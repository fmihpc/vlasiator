#include "../grid.h"
#include "../object_wrapper.h"
#include "../memoryallocation.h"
#include "vec.h"
#include "cpu_1d_ppm_nonuniform.hpp"
//#include "cpu_1d_ppm_nonuniform_conserving.hpp"

#include "gpu_trans_map_amr.hpp"
#include "cpu_trans_pencils.hpp"
#include "../arch/gpu_base.hpp"

// indices in padded source block, which is of type Vec with VECL
// elements in each vector.
//#define i_trans_ps_blockv_pencil(planeVectorIndex, planeIndex, blockIndex, lengthOfPencil) ( (blockIndex)  +  ( (planeVectorIndex) + (planeIndex) * VEC_PER_PLANE ) * ( lengthOfPencil) )
#define i_trans_ps_blockv_pencil(planeIndex, blockIndex, lengthOfPencil) ( (blockIndex)  +  ( (planeIndex) * VEC_PER_PLANE ) * ( lengthOfPencil) )

// Skip remapping if whole stencil for all vector elements consists of zeroes
__host__ __device__ inline bool check_skip_remapping(Vec* values, uint vectorindex) {
   for (int index=-VLASOV_STENCIL_WIDTH; index<VLASOV_STENCIL_WIDTH+1; ++index) {
      if (values[index][vectorindex] > 0) return false;
   }
   return true;
}

/* Propagate a given velocity block in all spatial cells of a pencil by a time step dt using a PPM reconstruction.
 *
 * @param pencilDZ Width of spatial cells in the direction of the pencil, vector datatype
 * @param values Density values of the block, vector datatype
 * @param dimension Satial dimension
 * @param blockGID Global ID of the velocity block.
 * @param dt Time step
 * @param vmesh Velocity mesh object
 * @param lengthOfPencil Number of cells in the pencil
 */
/* Copy the pencil source data to the temporary values array, so that the
 * dimensions are correctly swapped.
 *
 * This function must be thread-safe.
 *
 * @param blockDataPointer Vector of pre-prepared pointers to input (cell) block data
 * @param start Index from blockDataPointer to start at
 * @param int lengthOfPencil Number of spatial cells in pencil (not inclusive 2*VLASOV_STENCIL_WIDTH
 * @param values Vector into which the data should be loaded
 * @param vcell_transpose
 * @param popID ID of the particle species.
 */

//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor, maxBlocksPerCluster)
// assume never more than 8 threads/CPUs per GPU
__global__ void __launch_bounds__(WID3, 4) translation_kernel(
   const uint dimension,
   const unsigned int* const vcell_transpose,
   const Realv dt,
   uint* pencilLengths,
   uint* pencilStarts,
   //const vmesh::GlobalID blockGID, // which GID this thread is working on
   vmesh::GlobalID *allBlocks, // List of all blocks
   const uint nAllBlocks, // size of list of blocks which we won't exceed
   const uint startingBlockIndex, // First block index for this kernel invocation
   const uint blockIndexIncrement, // How much each kernel invocation should jump ahead
   const uint nPencils, // Number of total pencils (constant)
   const uint sumOfLengths, // sum of all pencil lengths (constant)
   const Realv threshold, // used by slope limiters
   split::SplitVector<vmesh::VelocityMesh*> *allPencilsMeshes, // Pointers to velocity meshes
   split::SplitVector<vmesh::VelocityBlockContainer*> *allPencilsContainers, // pointers to BlockContainers
   Realf** pencilBlockData, // pointers into cell block data, both written and read
   Vec* pencilOrderedSource, // Vec-ordered block data values for pencils
   Realf* pencilDZ,
   Realf* pencilRatios, // Vector holding target ratios
   uint* pencilBlocksCount // store how many non-empty blocks each pencil has for this GID
   ) {

   //const int gpuBlocks = gridDim.x;
   //const int blocki = blockIdx.x;
   //const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   // This is launched with block size (WID2,WID,1) assuming that VECL==WID2
   const vmesh::LocalID ti = threadIdx.y*blockDim.x + threadIdx.x;

   // offsets so this block of the kernel uses the correct part of temp arrays
   const uint pencilBlockDataOffset = blockIdx.x * sumOfLengths;
   const uint pencilOrderedSourceOffset = blockIdx.x * sumOfLengths * (WID3/VECL);
   const uint pencilBlocksCountOffset = blockIdx.x * nPencils;

   vmesh::VelocityMesh** pencilMeshes = allPencilsMeshes->data();
   vmesh::VelocityBlockContainer** pencilContainers = allPencilsContainers->data();
   vmesh::VelocityMesh* vmesh = pencilMeshes[0]; // just some vmesh
   const Realv dvz = vmesh->getCellSize()[dimension];
   const Realv vz_min = vmesh->getMeshMinLimits()[dimension];

   // Acting on velocity block blockGID, now found from array
   for (uint thisBlockIndex = startingBlockIndex + blockIdx.x; thisBlockIndex < nAllBlocks; thisBlockIndex += blockIndexIncrement) {
      const uint blockGID = allBlocks[thisBlockIndex];
      // First read data in
      for (uint pencili=0; pencili<nPencils; pencili++) {
         const uint lengthOfPencil = pencilLengths[pencili];
         const uint start = pencilStarts[pencili];
         // Get pointer to temprary buffer of VEC-ordered data for this kernel
         Vec* thisPencilOrderedSource = pencilOrderedSource + pencilOrderedSourceOffset + start * WID3/VECL;
         uint nonEmptyBlocks = 0;
         // Go over pencil length, gather cellblock data into aligned pencil source data
         for (uint celli = 0; celli < lengthOfPencil; celli++) {
            vmesh::VelocityMesh* vmesh = pencilMeshes[start + celli];
            // GPUTODO: Should we use the accelerated Hashinator interface to prefetch all LID-GID-pairs?
            // const vmesh::LocalID blockLID = vmesh->getLocalID(blockGID);
            // Now using warp accessor.
            const vmesh::LocalID blockLID = vmesh->warpGetLocalID(blockGID,ti);
            // Store block data pointer for both loading of data and writing back to the cell
            if (blockLID == vmesh->invalidLocalID()) {
               if (ti==0) {
                  pencilBlockData[pencilBlockDataOffset + start + celli] = NULL;
               }
               __syncthreads();
               // Non-existing block, push in zeroes
               thisPencilOrderedSource[i_trans_ps_blockv_pencil(threadIdx.y, celli, lengthOfPencil)][threadIdx.x] = 0.0;
            } else {
               #ifdef DEBUG_VLASIATOR
               const vmesh::LocalID meshSize = vmesh->size();
               const vmesh::LocalID VBCSize = pencilContainers[start + celli]->size();
               if ((blockLID>=meshSize) || (blockLID>=VBCSize)) {
                  if (ti==0) {
                     printf("Error in translation: trying to access LID %ul but sizes are vmesh %ul VBC %ul\n",blockLID,meshSize,VBCSize);
                  }
               }
               #endif
               if (ti==0) {
                  pencilBlockData[pencilBlockDataOffset + start + celli] = pencilContainers[start + celli]->getData(blockLID);
                  nonEmptyBlocks++;
               }
               __syncthreads();
               // Valid block, store values in Vec-order for efficient reading in propagation
               // Transpose block values so that mapping is along k direction.
               thisPencilOrderedSource[i_trans_ps_blockv_pencil(threadIdx.y, celli, lengthOfPencil)][threadIdx.x]
                  = (pencilBlockData[pencilBlockDataOffset + start + celli])[vcell_transpose[ti]];
            }
         } // End loop over this pencil
         if (ti==0) {
            pencilBlocksCount[pencilBlocksCountOffset + pencili] = nonEmptyBlocks;
         }
         __syncthreads();
      } // end loop over all pencils

      __syncthreads();
      // Now we reset target blocks
      for (uint celli=0; celli<sumOfLengths; celli++) {
         if (pencilRatios[celli] != 0) {
            // Is a target cell, needs to be reset
            vmesh::VelocityMesh* vmesh = pencilMeshes[celli];
            //const vmesh::LocalID blockLID = vmesh->getLocalID(blockGID);
            const vmesh::LocalID blockLID = vmesh->warpGetLocalID(blockGID,ti);
            if (blockLID != vmesh->invalidLocalID()) {
               // This block exists for this cell, reset
               (pencilBlockData[pencilBlockDataOffset + celli])[ti] = 0.0;
            }
         }
      } // end loop over all cells

      // Now we propagate the pencils and write data back to the block data containers
      // Get velocity data from vmesh that we need later to calculate the translation
      __syncthreads();
      vmesh::LocalID blockIndicesD = 0;
      if (dimension==0) {
         vmesh->getIndicesX(blockGID, blockIndicesD);
      } else if (dimension==1) {
         vmesh->getIndicesY(blockGID, blockIndicesD);
      } else if (dimension==2) {
         vmesh->getIndicesZ(blockGID, blockIndicesD);
      }

      // Assuming 1 neighbor in the target array because of the CFL condition
      // In fact propagating to > 1 neighbor will give an error
      // Also defined in the calling function for the allocation of targetValues
      // const uint nTargetNeighborsPerPencil = 1;
      for (uint pencili=0; pencili<nPencils; pencili++) {
         if (pencilBlocksCount[pencilBlocksCountOffset + pencili] == 0) {
            continue;
         }
         const uint lengthOfPencil = pencilLengths[pencili];
         const uint start = pencilStarts[pencili];
         Vec* thisPencilOrderedSource = pencilOrderedSource + pencilOrderedSourceOffset + start * WID3/VECL;

         // Go over length of propagated cells
         for (uint i = VLASOV_STENCIL_WIDTH; i < lengthOfPencil-VLASOV_STENCIL_WIDTH; i++){
            // Get pointers to block data used for output.
            // GPUTODO: use blockGID to get pointers here

            Realf* block_data_m1 = pencilBlockData[pencilBlockDataOffset + start + i - 1];
            Realf* block_data =    pencilBlockData[pencilBlockDataOffset + start + i];
            Realf* block_data_p1 = pencilBlockData[pencilBlockDataOffset + start + i + 1];

            // Cells which shouldn't be written to (e.g. sysboundary cells) have a targetRatio of 0
            // Also need to check if pointer is valid, because a cell can be missing an elsewhere propagated block
            Realf areaRatio_m1 = pencilRatios[start + i - 1];
            Realf areaRatio =    pencilRatios[start + i];
            Realf areaRatio_p1 = pencilRatios[start + i + 1];

            // (no longer loop over) planes (threadIdx.y) and vectors within planes (just 1 by construction)
            const Realf cell_vz = (blockIndicesD * WID + threadIdx.y + 0.5) * dvz + vz_min; //cell centered velocity
            const Realf z_translation = cell_vz * dt / pencilDZ[i]; // how much it moved in time dt (reduced units)

            // Determine direction of translation
            // part of density goes here (cell index change along spatial direcion)
            bool positiveTranslationDirection = (z_translation > 0.0);

            // Calculate normalized coordinates in current cell.
            // The coordinates (scaled units from 0 to 1) between which we will
            // integrate to put mass in the target  neighboring cell.
            // Normalize the coordinates to the origin cell. Then we scale with the difference
            // in volume between target and origin later when adding the integrated value.
            Realf z_1,z_2;
            z_1 = positiveTranslationDirection ? 1.0 - z_translation : 0.0;
            z_2 = positiveTranslationDirection ? 1.0 : - z_translation;

            #ifdef DEBUG_VLASIATOR
            if( horizontal_or(abs(z_1) > Vec(1.0)) || horizontal_or(abs(z_2) > Vec(1.0)) ) {
               assert( 0 && "Error in translation, CFL condition violated.");
            }
            #endif

            // Check if all values are 0:
            if (!check_skip_remapping(thisPencilOrderedSource
                                      + i_trans_ps_blockv_pencil(threadIdx.y, i, lengthOfPencil),threadIdx.x)) {
               // Note, you cannot request to sync threads within this code block or you risk deadlock.

               // Compute polynomial coefficients
               Realf a[3];
               // Silly indexing into coefficient calculation necessary due to
               // built-in assumptions of unsigned indexing.
               compute_ppm_coeff_nonuniform(pencilDZ + i - VLASOV_STENCIL_WIDTH,
                                            thisPencilOrderedSource
                                            + i_trans_ps_blockv_pencil(threadIdx.y, i, lengthOfPencil)
                                            - VLASOV_STENCIL_WIDTH,
                                            h4, VLASOV_STENCIL_WIDTH, a, threshold, threadIdx.x);
               // Compute integral
               const Realf ngbr_target_density =
                  z_2 * ( a[0] + z_2 * ( a[1] + z_2 * a[2] ) ) -
                  z_1 * ( a[0] + z_1 * ( a[1] + z_1 * a[2] ) );

               // Store mapped density in two target cells
               // in the current original cells we will put the rest of the original density
               // Now because each GPU block handles all pencils for an unique GID, we shouldn't need atomic additions here.

               // NOTE: not using atomic operations causes huge diffs (as if self contribution was neglected)! 11.01.2024 MB
               if (areaRatio && block_data) {
                  const Realf selfContribution = (thisPencilOrderedSource[i_trans_ps_blockv_pencil(threadIdx.y, i, lengthOfPencil)][threadIdx.x] - ngbr_target_density) * areaRatio;
                  atomicAdd(&block_data[vcell_transpose[ti]],selfContribution);
                  //block_data[vcell_transpose[ti]] += selfContribution;
               }
               if (areaRatio_p1 && block_data_p1) {
                  const Realf p1Contribution = (positiveTranslationDirection ? ngbr_target_density
                                                * pencilDZ[i] / pencilDZ[i + 1] : 0.0) * areaRatio_p1;
                  atomicAdd(&block_data_p1[vcell_transpose[ti]],p1Contribution);
                  //block_data_p1[vcell_transpose[ti]] += p1Contribution;
               }
               if (areaRatio_m1 && block_data_m1) {
                  const Realf m1Contribution = (!positiveTranslationDirection ? ngbr_target_density
                                                * pencilDZ[i] / pencilDZ[i - 1] : 0.0) * areaRatio_m1;
                  atomicAdd(&block_data_m1[vcell_transpose[ti]],m1Contribution);
                  //block_data_m1[vcell_transpose[ti]] += m1Contribution;
               }
            } // Did not skip remapping
            __syncthreads();
         } // end loop over this pencil
      } // end loop over all pencils
      __syncthreads();
   } // end loop over blocks
}

/* Mini-kernel for looping over all available velocity meshes and gathering
 * the union of all existing blocks.
 *
 * @param unionOfBlocksSet Hashmap, where keys are those blocks which are in the union of all blocks
 * @param allVmeshPointer Vector of pointers to velocitymeshes, used for gathering active blocks
 * @param nAllCells count of cells to read from allVmeshPointer
 */
__global__ void  __launch_bounds__(GPUTHREADS, 4) gather_union_of_blocks_kernel(
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet,
   split::SplitVector<vmesh::VelocityMesh*> *allVmeshPointer,
   const uint nAllCells)
{
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   //const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const vmesh::LocalID ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   for (vmesh::LocalID cellIndex=blocki; cellIndex<nAllCells; cellIndex += gpuBlocks) {
      vmesh::VelocityMesh* thisVmesh = allVmeshPointer->at(cellIndex);
      vmesh::LocalID nBlocks = thisVmesh->size();
      // for (vmesh::LocalID blockIndex=ti; blockIndex<nBlocks; blockIndex += warpSize) {
      //    const vmesh::GlobalID GID = thisVmesh->getGlobalID(blockIndex);
      //    unionOfBlocksSet->set_element(GID,GID);
      // }
      // Now with warp accessors
      for (vmesh::LocalID blockIndex=0; blockIndex<nBlocks; blockIndex++) {
         const vmesh::GlobalID GID = thisVmesh->getGlobalID(blockIndex);
         // warpInsert<true> only inserts if key does not yet exist
         unionOfBlocksSet->warpInsert<true>(GID, (vmesh::LocalID)GID, ti % GPUTHREADS);
      }
   }
}

/* Map velocity blocks in all local cells forward by one time step in one spatial dimension.
 * This function uses 1-cell wide pencils to update cells in-place to avoid allocating large
 * temporary buffers.
 *
 * @param [in] mpiGrid DCCRG grid object
 * @param [in] localPropagatedCells List of local cells that get propagated
 * ie. not boundary or DO_NOT_COMPUTE
 * @param [in] remoteTargetCells List of non-local target cells
 * @param [in] dimension Spatial dimension
 * @param [in] dt Time step
 * @param [in] popId Particle population ID
 */
bool gpu_trans_map_1d_amr(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                      const vector<CellID>& localPropagatedCells,
                      const vector<CellID>& remoteTargetCells,
                      std::vector<uint>& nPencils,
                      const uint dimension,
                      const Realv dt,
                      const uint popID) {

   phiprof::Timer setupTimer {"trans-amr-setup"};

   // return if there's no cells to propagate
   if(localPropagatedCells.size() == 0) {
      cout << "Returning because of no cells" << endl;
      return false;
   }

   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
   unsigned int vcell_transpose[WID3]; /*< defines the transpose for the solver internal (transposed) id: i + j*WID + k*WID2 to actual one*/
   // Fiddle indices x,y,z in VELOCITY SPACE
   switch (dimension) {
   case 0:
      // set values in array that is used to convert block indices
      // to global ID using a dot product.
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
   case 1:
      // set values in array that is used to convert block indices
      // to global ID using a dot product
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
   case 2:
      // set values in array that is used to convert block indices
      // to global id using a dot product.
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   default:
      cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
      abort();
      break;
   }

   // init vcell_transpose (moved here to take advantage of the omp parallel region)
#pragma omp parallel for collapse(3)
   for (uint k=0; k<WID; ++k) {
      for (uint j=0; j<WID; ++j) {
         for (uint i=0; i<WID; ++i) {
            const uint cell =
               i * cell_indices_to_id[0] +
               j * cell_indices_to_id[1] +
               k * cell_indices_to_id[2];
            vcell_transpose[ i + j * WID + k * WID2] = cell;
         }
      }
   }
   // Copy indexing information to device.
   gpuStream_t bgStream = gpu_getStream(); // uses stream assigned to thread 0, not the blocking default stream
   int device = gpu_getDevice();
   CHK_ERR( gpuMemcpyAsync(gpu_vcell_transpose, vcell_transpose, WID3*sizeof(uint), gpuMemcpyHostToDevice,bgStream) );

   // Vector with all cell ids
   vector<CellID> allCells(localPropagatedCells);
   allCells.insert(allCells.end(), remoteTargetCells.begin(), remoteTargetCells.end());
   const uint nAllCells = allCells.size();

   // Vectors of pointers to the cell structs
   std::vector<SpatialCell*> allCellsPointer(nAllCells);

   // Ensure GPU data has sufficient allocations/sizes, perform prefetches to CPU
   cuint sumOfLengths = DimensionPencils[dimension].sumOfLengths;
   gpu_vlasov_allocate(sumOfLengths);
   gpu_trans_allocate(nAllCells,sumOfLengths,0,0);
   // Initialize allCellsPointer. Find maximum mesh size.
   uint largestFoundMeshSize = 0;
   #pragma omp parallel
   {
      uint thread_largestFoundMeshSize = 0;
      #pragma omp for
      for(uint celli = 0; celli < nAllCells; celli++){
         allCellsPointer[celli] = mpiGrid[allCells[celli]];
         mpiGrid[allCells[celli]]->dev_upload_population(popID);
         allVmeshPointer->at(celli) = mpiGrid[allCells[celli]]->dev_get_velocity_mesh(popID);
         const uint thisMeshSize = mpiGrid[allCells[celli]]->get_velocity_mesh(popID)->size();
         thread_largestFoundMeshSize = thisMeshSize > thread_largestFoundMeshSize ? thisMeshSize : thread_largestFoundMeshSize;
         // Prefetches (in fact all data should already reside in device memory)
         // allCellsPointer[celli]->get_velocity_mesh(popID)->gpu_prefetchDevice();
         // allCellsPointer[celli]->get_velocity_blocks(popID)->gpu_prefetchDevice();
      }
      #pragma omp critical
      {
         largestFoundMeshSize = largestFoundMeshSize > thread_largestFoundMeshSize ? largestFoundMeshSize : thread_largestFoundMeshSize;
      }
   }
   // Prefetch vector of vmesh pointers to GPU
   allVmeshPointer->optimizeGPU(bgStream);

   // Reserve size for unionOfBlocksSet
   gpu_trans_allocate(0, 0, largestFoundMeshSize, 0);

   // Gather cell weights for load balancing
   if (Parameters::prepareForRebalance == true) {
      for (uint i=0; i<localPropagatedCells.size(); i++) {
         cuint myPencilCount = std::count(DimensionPencils[dimension].ids.begin(), DimensionPencils[dimension].ids.end(), localPropagatedCells[i]);
         nPencils[i] += myPencilCount;
         nPencils[nPencils.size()-1] += myPencilCount;
      }
   }

   phiprof::Timer buildTimer {"trans-amr-buildBlockList"};
   // Get a unique unsorted list of blockids that are in any of the
   // propagated cells. We launch this kernel, and do host-side pointer
   // gathering in parallel with it.
   const uint nGpuBlocks = nAllCells > GPUBLOCKS ? GPUBLOCKS : nAllCells;
   gather_union_of_blocks_kernel<<<nGpuBlocks, GPUTHREADS, 0, bgStream>>> (
      unionOfBlocksSet,
      allVmeshPointer,
      nAllCells
      );
   CHK_ERR( gpuPeekAtLastError() );
   buildTimer.stop();

   phiprof::Timer gatherPointerTimer {"trans-amr-gather-meshpointers"};
   // For each cellid listed in the pencils for this dimension, store the pointer to the vmesh.
   // At the same time, we could accumulate a list of unique cells included, but we already
   // get these from vlasovmover. This has to be on the host, as SpatialCells reside in host memory.
   // Host-side hashmaps are not threadsafe, so the allCellsPointer list cannot feasibly be gathered here.
   #pragma omp parallel for
   for (uint pencili = 0; pencili < DimensionPencils[dimension].N; ++pencili) {
      int L = DimensionPencils[dimension].lengthOfPencils[pencili];
      int start = DimensionPencils[dimension].idsStart[pencili];
      // Loop over cells in pencil
      for (int i = 0; i < L; i++) {
         const CellID thisCell = DimensionPencils[dimension].ids[start+i];
         allPencilsMeshes->at(start+i) = mpiGrid[thisCell]->dev_get_velocity_mesh(popID);
         allPencilsContainers->at(start+i) = mpiGrid[thisCell]->dev_get_velocity_blocks(popID);
      }
   }
   // Prefetch data back to GPU
   allPencilsMeshes->optimizeGPU(bgStream);
   allPencilsContainers->optimizeGPU(bgStream);

   // Extract pointers to data in managed memory
   uint* pencilLengths = DimensionPencils[dimension].gpu_lengthOfPencils->data();
   uint* pencilStarts = DimensionPencils[dimension].gpu_idsStart->data();
   Realf* pencilDZ = DimensionPencils[dimension].gpu_sourceDZ->data();
   Realf* pencilRatios = DimensionPencils[dimension].gpu_targetRatios->data();
   gatherPointerTimer.stop();

   phiprof::Timer buildTimer2 {"trans-amr-buildBlockList-2"};
   // Now we ensure the union of blocks gathering is complete and extract the union of blocks into a vector
   CHK_ERR( gpuStreamSynchronize(bgStream) );
   const vmesh::LocalID unionOfBlocksSetSize = unionOfBlocksSet->size();
   gpu_trans_allocate(0,0,0,unionOfBlocksSetSize);
   const uint nAllBlocks = unionOfBlocksSet->extractAllKeys(*unionOfBlocks,bgStream);
   CHK_ERR( gpuStreamSynchronize(bgStream) );
   vmesh::GlobalID *allBlocks = unionOfBlocks->data();
   // This threshold value is used by slope limiters.
   Realv threshold = mpiGrid[DimensionPencils[dimension].ids[VLASOV_STENCIL_WIDTH]]->getVelocityBlockMinValue(popID);
   buildTimer2.stop();

   /***********************/
   setupTimer.stop();
   /***********************/
   int bufferId {phiprof::initializeTimer("prepare buffers")};
   int mappingId {phiprof::initializeTimer("trans-amr-mapping")};
   #pragma omp parallel
   {
      // Thread id used for persistent device memory pointers
      #ifdef _OPENMP
      const uint cpuThreadID = omp_get_thread_num();
      const uint maxThreads = omp_get_max_threads();
      #else
      const uint cpuThreadID = 0;
      const uint maxThreads = 1;
      #endif
      gpuStream_t stream = gpu_getStream();

      phiprof::Timer bufferTimer {bufferId};
      // Vector of pointers to cell block data, used for both reading and writing
      // GPUTODO: pre-allocate one per thread, here verify sufficient size
      cuint nPencils = DimensionPencils[dimension].N;
      cuint currentAllocation = gpu_vlasov_getAllocation(); // in blocks, per thread buffer
      Vec* pencilOrderedSource = gpu_blockDataOrdered[cpuThreadID]; // pre-allocated temp buffer
      // How many blocks can each thread manage in parallel with this existing temp buffer?
      cuint nBlocksPerThread = currentAllocation / sumOfLengths;
      uint nGpuBlocks  = nBlocksPerThread > GPUBLOCKS ? GPUBLOCKS : nBlocksPerThread;

      Realf** pencilBlockData; // Array of pointers into actual block data
      uint* pencilBlocksCount; // Array of counters if pencil needs to be propagated for this block or not

      stringstream ss;
      ss<<" thread "<<cpuThreadID<<" malloc "<<sumOfLengths<<" * "<<nGpuBlocks<<" * sizeof(Realf*) + ";
      ss<<nPencils<<" * "<<nGpuBlocks<<" * sizeof(uint) = "<<sumOfLengths*nGpuBlocks*sizeof(Realf*)+nPencils*nGpuBlocks*sizeof(uint)<<std::endl;
      std::cerr<<ss.str();
#pragma omp barrier
      CHK_ERR( gpuMallocAsync((void**)&pencilBlockData, sumOfLengths*nGpuBlocks*sizeof(Realf*), stream) );
      CHK_ERR( gpuMallocAsync((void**)&pencilBlocksCount, nPencils*nGpuBlocks*sizeof(uint), stream) );
      CHK_ERR( gpuStreamSynchronize(stream) );
      bufferTimer.stop();

      // Loop over velocity space blocks (threaded, multi-stream, and multi-block parallel, but not using a for-loop)
      phiprof::Timer mappingTimer {mappingId}; // mapping (top-level)
      const uint startingBlockIndex = cpuThreadID*nGpuBlocks;
      const uint blockIndexIncrement = maxThreads*nGpuBlocks;
      // This thread, using its own stream, will launch nGpuBlocks instances of the below kernel, where each instance
      // propagates all pencils for the block in question.
      dim3 block(WID2,WID,1); // assumes VECL==WID2
#pragma omp barrier
      translation_kernel<<<nGpuBlocks, block, 0, stream>>> (
         dimension,
         gpu_vcell_transpose,
         dt,
         pencilLengths,
         pencilStarts,
         //blockGID, // which GID this thread is working on
         allBlocks, // List of all blocks
         nAllBlocks, // size of list of blocks which we won't exceed
         startingBlockIndex, // First block index for this kernel invocation
         blockIndexIncrement, // How much each kernel invocation should jump ahead
         nPencils, // Number of total pencils (constant)
         sumOfLengths, // sum of all pencil lengths (constant)
         threshold,
         allPencilsMeshes, // Pointers to velocity meshes
         allPencilsContainers, // pointers to BlockContainers
         pencilBlockData, // pointers into cell block data, both written and read
         pencilOrderedSource, // Vec-ordered block data values for pencils
         pencilDZ,
         pencilRatios, // Vector holding target ratios
         pencilBlocksCount // store how many non-empty blocks each pencil has for this GID
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuStreamSynchronize(stream) );
      mappingTimer.stop(); // mapping (top-level)
      CHK_ERR( gpuFree(pencilBlockData) );
      CHK_ERR( gpuFree(pencilBlocksCount) );

   } // closes pragma omp parallel
   return true;
}


/* Get an index that identifies which cell in the list of sibling cells this cell is.
 *
 * @param mpiGrid DCCRG grid object
 * @param cellid DCCRG id of this cell
 */
int get_sibling_index(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const CellID& cellid) {

   const int NO_SIBLINGS = 0;
   if(mpiGrid.get_refinement_level(cellid) == 0) {
      return NO_SIBLINGS;
   }

   //CellID parent = mpiGrid.mapping.get_parent(cellid);
   CellID parent = mpiGrid.get_parent(cellid);

   if (parent == INVALID_CELLID) {
      std::cerr<<"Invalid parent id"<<std::endl;
      abort();
   }

   // get_all_children returns an array instead of a vector now, need to map it to a vector for find and distance
   // std::array<uint64_t, 8> siblingarr = mpiGrid.mapping.get_all_children(parent);
   // vector<CellID> siblings(siblingarr.begin(), siblingarr.end());
   vector<CellID> siblings = mpiGrid.get_all_children(parent);
   auto location = std::find(siblings.begin(),siblings.end(),cellid);
   auto index = std::distance(siblings.begin(), location);
   if (index>7) {
      std::cerr<<"Invalid parent id"<<std::endl;
      abort();
   }
   return index;

}

/** GPU kernel for incrementing one set of blocks with another
    used for accruing remote neighbor contributions
*/
__global__ static void __launch_bounds__(WID3,4) remote_increment_kernel (
   Realf* blockData,
   Realf* neighborData,
   vmesh::LocalID nBlocks
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;
   // loop over whole velocity space
   for (uint blockLID=blocki; blockLID<nBlocks; blockLID += gpuBlocks) {
      // Increment value
      atomicAdd(&blockData[blockLID * WID3 + ti],neighborData[blockLID * WID3 + ti]);

      //blockData[blockLID * WID3 + ti] += neighborData[blockLID * WID3 + ti];
      // Note: this is not an atomic operation, so only one kernel per cell can be active at a time.
   }
}

/* This function communicates the mapping on process boundaries, and then updates the data to their correct values.
 * When sending data between neighbors of different refinement levels, special care has to be taken to ensure that
 * The sending and receiving ranks allocate the correct size arrays for neighbor_block_data.
 * This is partially due to DCCRG defining neighborhood size relative to the host cell. For details, see
 * https://github.com/fmihpc/dccrg/issues/12
 *
 * @param mpiGrid DCCRG grid object
 * @param dimension Spatial dimension
 * @param direction Direction of communication (+ or -)
 * @param popId Particle population ID
 */
void gpu_update_remote_mapping_contribution_amr(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const uint dimension,
   int direction,
   const uint popID) {

   // GPUTODO: First attempts at using managed memory for remote neighbours
   // Should move to re-using managed memory buffers and ensuring size is suitable?
   // If that path is taken, it should also check for any local cells *not* on process
   // boundary and free the buffers from those cells.
   int device = gpu_getDevice();

   int neighborhood = 0;
   int both_neighborhood = 0;

   //normalize and set neighborhoods
   if(direction > 0) {
      direction = 1;
      switch (dimension) {
      case 0:
         neighborhood = SHIFT_P_X_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = SHIFT_P_Y_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = SHIFT_P_Z_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID;
         break;
      }
   }
   if(direction < 0) {
      direction = -1;
      switch (dimension) {
      case 0:
         neighborhood = SHIFT_M_X_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = SHIFT_M_Y_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = SHIFT_M_Z_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID;
         break;
      }
   }

   //const vector<CellID>& local_cells = getLocalCells();
   const vector<CellID>& local_cells = mpiGrid.get_local_cells_on_process_boundary(both_neighborhood);
   const vector<CellID> remote_cells = mpiGrid.get_remote_cells_on_process_boundary(both_neighborhood);

   // Fast return if no cells to process
   if (local_cells.size() == remote_cells.size() == 0) {
      return;
   }

   vector<CellID> receive_cells;
   set<CellID> send_cells;
   vector<CellID> receive_origin_cells;
   vector<uint> receive_origin_index;

   phiprof::Timer updateRemoteTimerPre {"trans-amr-remotes-setup-getcells"};
   // Initialize remote cells
   #pragma omp parallel for
   for (auto rc : remote_cells) {
      SpatialCell *ccell = mpiGrid[rc];
      // Initialize number of blocks to 0 and block data to a default value.
      // We need the default for 1 to 1 communications
      if(ccell) {
         for (uint i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
            ccell->neighbor_block_data[i] = ccell->get_data(popID);
            ccell->neighbor_number_of_blocks[i] = 0;
         }
      }
   }

   // Initialize local cells
   #pragma omp parallel for
   for (auto lc : local_cells) {
      SpatialCell *ccell = mpiGrid[lc];
      if(ccell) {
         // Initialize number of blocks to 0 and neighbor block data pointer to the local block data pointer
         for (uint i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
            ccell->neighbor_block_data[i] = ccell->get_data(popID);
            ccell->neighbor_number_of_blocks[i] = 0;
         }
      }
   }
   updateRemoteTimerPre.stop();

   vector<Realf*> receiveBuffers;
   vector<Realf*> sendBuffers;

   phiprof::Timer updateRemoteTimer0 {"trans-amr-remotes-setup-localcells"};
   for (auto c : local_cells) {
      SpatialCell *ccell = mpiGrid[c];
      if (!ccell) continue;
      vector<CellID> p_nbrs;
      vector<CellID> n_nbrs;
      for (const auto& nbr : mpiGrid.get_face_neighbors_of(c)) {
         if(nbr.second == ((int)dimension + 1) * direction) {
            p_nbrs.push_back(nbr.first);
         }
         if(nbr.second == -1 * ((int)dimension + 1) * direction) {
            n_nbrs.push_back(nbr.first);
         }
      }

      uint sendIndex = 0;
      uint recvIndex = 0;
      int mySiblingIndex = get_sibling_index(mpiGrid,c);
      // Set up sends if any neighbor cells in p_nbrs are non-local.
      if (!all_of(p_nbrs.begin(), p_nbrs.end(), [&mpiGrid](CellID i){return mpiGrid.is_local(i);})) {
         phiprof::Timer updateRemoteTimer1 {"trans-amr-remotes-setup-sends"};
         // ccell adds a neighbor_block_data block for each neighbor in the positive direction to its local data
         for (const auto nbr : p_nbrs) {
            //Send data in nbr target array that we just mapped to, if
            // 1) it is a valid target,
            // 2) the source cell in center was translated,
            // 3) Cell is remote.
            if(nbr != INVALID_CELLID && do_translate_cell(ccell) && !mpiGrid.is_local(nbr)) {
               /*
                 Select the index to the neighbor_block_data and neighbor_number_of_blocks arrays
                 1) Ref_c == Ref_nbr == 0, index = 0
                 2) Ref_c == Ref_nbr != 0, index = c sibling index
                 3) Ref_c >  Ref_nbr     , index = c sibling index
                 4) Ref_c <  Ref_nbr     , index = nbr sibling index
                */
               if(mpiGrid.get_refinement_level(c) >= mpiGrid.get_refinement_level(nbr)) {
                  sendIndex = mySiblingIndex;
               } else {
                  sendIndex = get_sibling_index(mpiGrid,nbr);
               }
               SpatialCell *pcell = mpiGrid[nbr];
               // 4) it exists and is not a boundary cell,
               if(pcell && pcell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {

                  ccell->neighbor_number_of_blocks.at(sendIndex) = pcell->get_number_of_velocity_blocks(popID);

                  if(send_cells.find(nbr) == send_cells.end()) {
                     // 5 We have not already sent data from this rank to this cell.
                     ccell->neighbor_block_data.at(sendIndex) = pcell->get_data(popID);
                     send_cells.insert(nbr);
                  } else {
                     // The receiving cell can't know which cell is sending the data from this rank.
                     // Therefore, we have to send 0's from other cells in the case where multiple cells
                     // from one rank are sending to the same remote cell so that all sent cells can be
                     // summed for the correct result.

                     if (ccell->neighbor_number_of_blocks.at(sendIndex) == 0) {
                        continue;
                     }
                     // GPUTODO: This is now unified memory. With GPU-aware MPI it could be on-device.
                     CHK_ERR( gpuMallocManaged((void**)&ccell->neighbor_block_data.at(sendIndex), ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf)) );
                     CHK_ERR( gpuMemPrefetchAsync(ccell->neighbor_block_data.at(sendIndex),ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf),device,0) );
                     CHK_ERR( gpuMemset(ccell->neighbor_block_data.at(sendIndex), 0, ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf)) );
                     sendBuffers.push_back(ccell->neighbor_block_data.at(sendIndex));
                  } // closes if(send_cells.find(nbr) == send_cells.end())
               } // closes if(pcell && pcell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)
            } // closes if(nbr != INVALID_CELLID && do_translate_cell(ccell) && !mpiGrid.is_local(nbr))
         } // closes for(uint i_nbr = 0; i_nbr < nbrs_to.size(); ++i_nbr)
      } // closes if(!all_of(nbrs_to.begin(), nbrs_to.end(),[&mpiGrid](CellID i){return mpiGrid.is_local(i);}))

      // Set up receives if any neighbor cells in n_nbrs are non-local.
      if (!all_of(n_nbrs.begin(), n_nbrs.end(), [&mpiGrid](CellID i){return mpiGrid.is_local(i);})) {
         phiprof::Timer updateRemoteTimer2 {"trans-amr-remotes-setup-receives"};
         // ccell adds a neighbor_block_data block for each neighbor in the positive direction to its local data
         for (const auto nbr : n_nbrs) {
            if (nbr != INVALID_CELLID && !mpiGrid.is_local(nbr) &&
                ccell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
               //Receive data that ncell mapped to this local cell data array,
               //if 1) ncell is a valid source cell, 2) center cell is to be updated (normal cell) 3) ncell is remote
               SpatialCell *ncell = mpiGrid[nbr];
               // Check for null pointer
               if(!ncell) {
                  continue;
               }
               /*
                 Select the index to the neighbor_block_data and neighbor_number_of_blocks arrays
                 1) Ref_nbr == Ref_c == 0, index = 0
                 2) Ref_nbr == Ref_c != 0, index = nbr sibling index
                 3) Ref_nbr >  Ref_c     , index = nbr sibling index
                 4) Ref_nbr <  Ref_c     , index = c   sibling index
                */
               if(mpiGrid.get_refinement_level(nbr) >= mpiGrid.get_refinement_level(c)) {
                  // Allocate memory for one sibling at recvIndex.
                  recvIndex = get_sibling_index(mpiGrid,nbr);
                  ncell->neighbor_number_of_blocks.at(recvIndex) = ccell->get_number_of_velocity_blocks(popID);
                  if (ncell->neighbor_number_of_blocks.at(recvIndex) == 0) {
                     continue;
                  }

                  // GPUTODO: This is now unified memory. With GPU-aware MPI it could be on-device.
                  CHK_ERR( gpuMallocManaged((void**)&ncell->neighbor_block_data.at(recvIndex), ncell->neighbor_number_of_blocks.at(recvIndex) * WID3 * sizeof(Realf)) );
                  CHK_ERR( gpuMemPrefetchAsync(ncell->neighbor_block_data.at(recvIndex), ncell->neighbor_number_of_blocks.at(recvIndex) * WID3 * sizeof(Realf), device,0) );
                  receiveBuffers.push_back(ncell->neighbor_block_data.at(recvIndex));
               } else {
                  recvIndex = mySiblingIndex;
                  // std::array<uint64_t, 8> siblingarr = mpiGrid.mapping.get_all_children(mpiGrid.mapping.get_parent(c));
                  // vector<CellID> mySiblings(siblingarr.begin(), siblingarr.end());
                  auto mySiblings = mpiGrid.get_all_children(mpiGrid.get_parent(c));
                  auto myIndices = mpiGrid.mapping.get_indices(c);

                  // Allocate memory for each sibling to receive all the data sent by coarser ncell.
                  // only allocate blocks for face neighbors.
                  for (uint i_sib = 0; i_sib < MAX_NEIGHBORS_PER_DIM; ++i_sib) {
                     auto sibling = mySiblings.at(i_sib);
                     auto sibIndices = mpiGrid.mapping.get_indices(sibling);
                     auto* scell = mpiGrid[sibling];
                     // Only allocate siblings that are remote face neighbors to ncell
                     // Also take care to have these consistent with the sending process neighbor checks!
                     if(sibling != INVALID_CELLID
                        && scell
                        && mpiGrid.get_process(sibling) != mpiGrid.get_process(nbr)
                        && myIndices.at(dimension) == sibIndices.at(dimension)
                        && ncell->neighbor_number_of_blocks.at(i_sib) != scell->get_number_of_velocity_blocks(popID)
                        && scell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {

                        ncell->neighbor_number_of_blocks.at(i_sib) = scell->get_number_of_velocity_blocks(popID);
                        if (ncell->neighbor_number_of_blocks.at(i_sib) == 0) {
                           continue;
                        }

                        // GPUTODO: This is now unified memory. With GPU-aware MPI it could be on-device.
                        CHK_ERR( gpuMallocManaged((void**)&ncell->neighbor_block_data.at(i_sib), ncell->neighbor_number_of_blocks.at(i_sib) * WID3 * sizeof(Realf)) );
                        CHK_ERR( gpuMemPrefetchAsync(ncell->neighbor_block_data.at(i_sib), ncell->neighbor_number_of_blocks.at(i_sib) * WID3 * sizeof(Realf), device,0) );
                        receiveBuffers.push_back(ncell->neighbor_block_data.at(i_sib));
                     }
                  }
               }
               receive_cells.push_back(c);
               receive_origin_cells.push_back(nbr);
               receive_origin_index.push_back(recvIndex);
            } // closes (nbr != INVALID_CELLID && !mpiGrid.is_local(nbr) && ...)
         } // closes for(uint i_nbr = 0; i_nbr < nbrs_of.size(); ++i_nbr)
      } // closes if(!all_of(nbrs_of.begin(), nbrs_of.end(),[&mpiGrid](CellID i){return mpiGrid.is_local(i);}))
   } // closes for (auto c : local_cells) {
   updateRemoteTimer0.stop();

   MPI_Barrier(MPI_COMM_WORLD);
   phiprof::Timer updateRemoteTimer3 {"trans-amr-remotes-MPI"};

   // Do communication
   SpatialCell::setCommunicatedSpecies(popID);
   SpatialCell::set_mpi_transfer_type(Transfer::NEIGHBOR_VEL_BLOCK_DATA);
   mpiGrid.update_copies_of_remote_neighbors(neighborhood);
   updateRemoteTimer3.stop();

   MPI_Barrier(MPI_COMM_WORLD);

   // Reduce data: sum received data in the data array to
   // the target grid in the temporary block container
   // #pragma omp parallel
   if (receive_cells.size() != 0) {
      phiprof::Timer updateRemoteTimerIncrement {"trans-amr-remotes-increment"};
      for (size_t c = 0; c < receive_cells.size(); ++c) {
         SpatialCell* receive_cell = mpiGrid[receive_cells[c]];
         SpatialCell* origin_cell = mpiGrid[receive_origin_cells[c]];
         if(!receive_cell || !origin_cell) {
            continue;
         }

         Realf *blockData = receive_cell->get_data(popID);
         Realf *neighborData = origin_cell->neighbor_block_data[receive_origin_index[c]];
         vmesh::LocalID nBlocks = receive_cell->get_number_of_velocity_blocks(popID);
         const uint nGpuBlocks = nBlocks > GPUBLOCKS ? GPUBLOCKS : nBlocks;
         #ifdef _OPENMP
         const uint maxThreads = omp_get_max_threads();
         #else
         const uint maxThreads = 1;
         #endif
         // Increment needs to be parallel-safe, so use modulo of cellid as stream number
         gpuStream_t cellStream = gpuStreamList[receive_cells[c] % maxThreads];
         if (nGpuBlocks>0) {
            dim3 block(WID,WID,WID);
            remote_increment_kernel<<<nGpuBlocks, block, 0, cellStream>>> (
               blockData,
               neighborData,
               nBlocks
               );
            CHK_ERR( gpuPeekAtLastError() );
            //CHK_ERR( gpuStreamSynchronize(cellStream) );
         }
      }
      // Since all increment kernel streams were launched from outside an openmp region, use device sync here.
      CHK_ERR( gpuDeviceSynchronize() );

      // send cell data is set to zero. This is to avoid double copy if
      // one cell is the neighbor on bot + and - side to the same process
      vector<CellID> send_cells_vector(send_cells.begin(), send_cells.end());
      #pragma omp parallel for
      for (uint c = 0; c < send_cells_vector.size(); c++) {
         SpatialCell* send_cell = mpiGrid[send_cells_vector[c]];
         gpuStream_t stream = gpu_getStream();
         Realf* blockData = send_cell->get_data(popID);
         CHK_ERR( gpuMemsetAsync(blockData, 0, WID3*send_cell->get_number_of_velocity_blocks(popID)*sizeof(Realf),stream) );
      }
   }

   phiprof::Timer updateRemoteTimerFree {"trans-amr-remotes-free"};
   for (auto p : receiveBuffers) {
      CHK_ERR( gpuFree(p) );
   }
   for (auto p : sendBuffers) {
      CHK_ERR( gpuFree(p) );
   }
   updateRemoteTimerFree.stop();
}
