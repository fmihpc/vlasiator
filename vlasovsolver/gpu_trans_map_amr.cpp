/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

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
#define i_trans_ps_blockv_pencil(planeIndex, blockIndex, lengthOfPencil) ( (blockIndex)  +  ( (planeIndex) * VEC_PER_PLANE ) * ( lengthOfPencil) )

// Skip remapping if whole stencil for all vector elements consists of zeroes
__host__ __device__ inline bool check_skip_remapping(Vec* values, uint vectorindex) {
   for (int index=-VLASOV_STENCIL_WIDTH; index<VLASOV_STENCIL_WIDTH+1; ++index) {
      if (values[index][vectorindex] > 0) return false;
   }
   return true;
}

/* Propagate a given velocity block in all spatial cells of a pencil by a time step dt using a PPM reconstruction.
   Includes preparing intermediate buffers.
 *
 * @param  dimension  Cartesian direction of propagation
 * @param  dt  Time step length
 * @param  pencilLengths  pointer to buffer of lengths of all pencils to propagate
 * @param  pencilStarts  pointer to buffer of indexes of first cells of pencils
 * @param  allBlocks  pointer to list of all GIDs to propagate
 * @param  nAllBlocks  how many blocks exist in total
 * @param  startingBlockIndex  First block index for this kernel invocation
 * @param  blockIndexIncrement  How much each kernel invocation should jump ahead
 * @param  nPencils  Number of total pencils (constant)
 * @param  sumOfLengths  sum of all pencil lengths (constant)
 * @param  threshold  sparsity threshold, used by slope limiters
 * @param  allPencilsMeshes  Pointer to vector of pointers to velocity meshes
 * @param  allPencilsContainers  Pointer to vector of pointers to BlockContainers
 * @param  pencilBlockData  Pointer to buffer of pointers into cell block data, both written and read
 * @param  pencilOrderedSource  Pointer to Vec-ordered buffer used as interim values
 * @param  pencilDZ  Pointer into buffer of pencil cell sizes
 * @param  pencilRatios  Pointer into buffer with pencil target ratios (due to AMR)
 * @param  pencilBlocksCount  Pointer into buffer for storing how many non-empty blocks each pencil has for current GID
 */

//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor, maxBlocksPerCluster)
__global__ void __launch_bounds__(WID3, 4) translation_kernel(
   const uint dimension,
   const Realv dt,
   uint* pencilLengths,
   uint* pencilStarts,
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
   // This is launched with block size (WID,WID,WID) assuming that VECL==WID2 and VEC_PER_BLOCK=WID
   const vmesh::LocalID ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const uint vecIndex = threadIdx.y*blockDim.x + threadIdx.x;
   const uint planeIndex = threadIdx.z;

   // Translation direction
   uint vz_index;
   switch (dimension) {
   case 0:
      vz_index = threadIdx.x;
      break;
   case 1:
      vz_index = threadIdx.y;
      break;
   case 2:
      vz_index = threadIdx.z;
      break;
   default:
      printf(" Wrong dimension, abort\n");
      break;
   }

   // offsets so this block of the kernel uses the correct part of temp arrays
   const uint pencilBlockDataOffset = blockIdx.x * sumOfLengths;
   const uint pencilOrderedSourceOffset = blockIdx.x * sumOfLengths * (WID3/VECL);
   const uint pencilBlocksCountOffset = blockIdx.x * nPencils;

   vmesh::VelocityMesh** pencilMeshes = allPencilsMeshes->data();
   vmesh::VelocityBlockContainer** pencilContainers = allPencilsContainers->data();
   vmesh::VelocityMesh* randovmesh = pencilMeshes[0]; // just some vmesh
   const Realv dvz = randovmesh->getCellSize()[dimension];
   const Realv vz_min = randovmesh->getMeshMinLimits()[dimension];

   // Acting on velocity block blockGID, now found from array
   for (uint thisBlockIndex = startingBlockIndex + blockIdx.x; thisBlockIndex < nAllBlocks; thisBlockIndex += blockIndexIncrement) {
      if (thisBlockIndex >= nAllBlocks) {
         break;
      }
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
               thisPencilOrderedSource[i_trans_ps_blockv_pencil(planeIndex, celli, lengthOfPencil)][vecIndex] = 0.0;
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
               thisPencilOrderedSource[i_trans_ps_blockv_pencil(planeIndex, celli, lengthOfPencil)][vecIndex]
                  = (pencilBlockData[pencilBlockDataOffset + start + celli])[ti];
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
            if (pencilBlockData[pencilBlockDataOffset + celli]) {
               (pencilBlockData[pencilBlockDataOffset + celli])[ti] = 0.0;
            }
            // vmesh::VelocityMesh* vmesh = pencilMeshes[celli];
            // //const vmesh::LocalID blockLID = vmesh->getLocalID(blockGID);
            // const vmesh::LocalID blockLID = vmesh->warpGetLocalID(blockGID,ti);
            // if (blockLID != vmesh->invalidLocalID()) {
            //    // This block exists for this cell, reset
            //    (pencilBlockData[pencilBlockDataOffset + celli])[ti] = 0.0;
            // }
         }
      } // end loop over all cells

      // Now we propagate the pencils and write data back to the block data containers
      // Get velocity data from vmesh that we need later to calculate the translation
      __syncthreads();
      vmesh::LocalID blockIndicesD = 0;
      if (dimension==0) {
         randovmesh->getIndicesX(blockGID, blockIndicesD);
      } else if (dimension==1) {
         randovmesh->getIndicesY(blockGID, blockIndicesD);
      } else if (dimension==2) {
         randovmesh->getIndicesZ(blockGID, blockIndicesD);
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
            Realf* block_data_m1 = pencilBlockData[pencilBlockDataOffset + start + i - 1];
            Realf* block_data =    pencilBlockData[pencilBlockDataOffset + start + i];
            Realf* block_data_p1 = pencilBlockData[pencilBlockDataOffset + start + i + 1];

            // Cells which shouldn't be written to (e.g. sysboundary cells) have a targetRatio of 0
            // Also need to check if pointer is valid, because a cell can be missing an elsewhere propagated block
            Realf areaRatio_m1 = pencilRatios[start + i - 1];
            Realf areaRatio =    pencilRatios[start + i];
            Realf areaRatio_p1 = pencilRatios[start + i + 1];

            // (no longer loop over) planes (threadIdx.z) and vectors within planes (just 1 by construction)
            const Realf cell_vz = (blockIndicesD * WID + vz_index + 0.5) * dvz + vz_min; //cell centered velocity
            const Realf z_translation = cell_vz * dt / pencilDZ[start + i]; // how much it moved in time dt (reduced units)

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
                                      + i_trans_ps_blockv_pencil(planeIndex, i, lengthOfPencil),vecIndex)) {
               // Note, you cannot request to sync threads within this code block or you risk deadlock.

               // Compute polynomial coefficients
               Realf a[3];
               // Silly indexing into coefficient calculation necessary due to
               // built-in assumptions of unsigned indexing.
               compute_ppm_coeff_nonuniform(pencilDZ + start + i - VLASOV_STENCIL_WIDTH,
                                            thisPencilOrderedSource
                                            + i_trans_ps_blockv_pencil(planeIndex, i, lengthOfPencil)
                                            - VLASOV_STENCIL_WIDTH,
                                            h4, VLASOV_STENCIL_WIDTH, a, threshold, vecIndex);
               // Compute integral
               const Realf ngbr_target_density =
                  z_2 * ( a[0] + z_2 * ( a[1] + z_2 * a[2] ) ) -
                  z_1 * ( a[0] + z_1 * ( a[1] + z_1 * a[2] ) );

               // Store mapped density in two target cells
               // in the current original cells we will put the rest of the original density
               // Now because each GPU block handles all pencils for an unique GID, we shouldn't need atomic additions here.

               // NOTE: not using atomic operations causes huge diffs (as if self contribution was neglected)! 11.01.2024 MB
               if (areaRatio && block_data) {
                  const Realf selfContribution = (thisPencilOrderedSource[i_trans_ps_blockv_pencil(planeIndex, i, lengthOfPencil)][vecIndex] - ngbr_target_density) * areaRatio;
                  //atomicAdd(&block_data[ti],selfContribution);
                  block_data[ti] += selfContribution;
               }
               if (areaRatio_p1 && block_data_p1) {
                  const Realf p1Contribution = (positiveTranslationDirection ? ngbr_target_density
                                                * pencilDZ[start + i] / pencilDZ[start + i + 1] : 0.0) * areaRatio_p1;
                  //atomicAdd(&block_data_p1[ti],p1Contribution);
                  block_data_p1[ti] += p1Contribution;
               }
               if (areaRatio_m1 && block_data_m1) {
                  const Realf m1Contribution = (!positiveTranslationDirection ? ngbr_target_density
                                                * pencilDZ[start + i] / pencilDZ[start + i - 1] : 0.0) * areaRatio_m1;
                  //atomicAdd(&block_data_m1[ti],m1Contribution);
                  block_data_m1[ti] += m1Contribution;
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
__global__ void  gather_union_of_blocks_kernel(
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet,
   split::SplitVector<vmesh::VelocityMesh*> *allVmeshPointer,
   const uint nAllCells)
{
   const int ti = threadIdx.x; // [0,GPUTHREADS)
   const int indexInBlock = threadIdx.y; // [0,WARPSPERBLOCK)
   const uint cellIndex = blockIdx.x;
   vmesh::VelocityMesh* thisVmesh = allVmeshPointer->at(cellIndex);
   uint thisVmeshSize = thisVmesh->size();
   for (uint blockIndex = indexInBlock; blockIndex < thisVmeshSize; blockIndex += WARPSPERBLOCK) {
      if (blockIndex < thisVmeshSize) {
         // Now with warp accessors
         const vmesh::GlobalID GID = thisVmesh->getGlobalID(blockIndex);
         // warpInsert<true> only inserts if key does not yet exist
         unionOfBlocksSet->warpInsert<true>(GID, (vmesh::LocalID)GID, ti);
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
                      std::vector<uint>& nPencilsLB,
                      const uint dimension,
                      const Realv dt,
                      const uint popID) {

   phiprof::Timer setupTimer {"trans-amr-setup"};

   // return if there's no cells to propagate
   if(localPropagatedCells.size() == 0) {
      return false;
   }
   gpuStream_t bgStream = gpu_getStream(); // uses stream assigned to thread 0, not the blocking default stream
   int device = gpu_getDevice();

   // Vector with all cell ids
   vector<CellID> allCells(localPropagatedCells);
   allCells.insert(allCells.end(), remoteTargetCells.begin(), remoteTargetCells.end());
   const uint nAllCells = allCells.size();

   phiprof::Timer allocateTimer {"trans-amr-allocs"};
   // Ensure GPU data has sufficient allocations/sizes, perform prefetches to CPU
   const uint sumOfLengths = DimensionPencils[dimension].sumOfLengths;
   gpu_vlasov_allocate(sumOfLengths);
   // Resize allVmeshPointer, allPencilsMeshes, allPencilsContainers
   gpu_trans_allocate(nAllCells,sumOfLengths,0,0);
   allocateTimer.stop();

   // Find maximum mesh size.
   phiprof::Timer maxMeshSizeTimer {"trans-amr-find-maxmesh"};
   uint largestFoundMeshSize = 0;
   int checkMeshId {phiprof::initializeTimer("trans-amr-checkMesh")};
   #pragma omp parallel
   {
      uint thread_largestFoundMeshSize = 0;
      #pragma omp for
      for(uint celli = 0; celli < nAllCells; celli++){
         allVmeshPointer->at(celli) = mpiGrid[allCells[celli]]->dev_get_velocity_mesh(popID); // GPU-side vmesh
         const uint thisMeshSize = mpiGrid[allCells[celli]]->get_velocity_mesh(popID)->size(); // get cached size from CPU side
         thread_largestFoundMeshSize = thisMeshSize > thread_largestFoundMeshSize ? thisMeshSize : thread_largestFoundMeshSize;
         #ifdef DEBUG_VLASIATOR
         phiprof::Timer checkMeshTimer {checkMeshId};
         if (!mpiGrid[allCells[celli]]->checkMesh(popID)) {
            printf("GPU TRANS MAP AMR check of mesh for popID %d cell %lu failed!\n",popID,allCells[celli]);
         }
         #endif
      }
      #pragma omp critical
      {
         largestFoundMeshSize = largestFoundMeshSize > thread_largestFoundMeshSize ? largestFoundMeshSize : thread_largestFoundMeshSize;
      }
   }
   maxMeshSizeTimer.stop();

   // return if there's no blocks to propagate
   if(largestFoundMeshSize == 0) {
      return false;
   }

   allocateTimer.start();
   // Prefetch vector of vmesh pointers to GPU
   allVmeshPointer->optimizeGPU(bgStream);
   // Reserve size for unionOfBlocksSet
   gpu_trans_allocate(0, 0, largestFoundMeshSize, 0);
   allocateTimer.stop();

   // Gather cell weights for load balancing
   phiprof::Timer pencilCountTimer {"trans-amr-count-pencils"};
   if (Parameters::prepareForRebalance == true) {
      for (uint i=0; i<localPropagatedCells.size(); i++) {
         const uint myPencilCount = std::count(DimensionPencils[dimension].ids.begin(), DimensionPencils[dimension].ids.end(), localPropagatedCells[i]);
         nPencilsLB[i] += myPencilCount;
         nPencilsLB[nPencilsLB.size()-1] += myPencilCount;
      }
   }
   pencilCountTimer.stop();

   phiprof::Timer buildTimer {"trans-amr-buildBlockList"};
   // Get a unique unsorted list of blockids that are in any of the
   // propagated cells. We launch this kernel, and do host-side pointer
   // gathering in parallel with it.
   dim3 gatherdims(GPUTHREADS,WARPSPERBLOCK,1);
   gather_union_of_blocks_kernel<<<nAllCells, gatherdims, 0, bgStream>>> (
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
   const uint nPencils = DimensionPencils[dimension].N;
   #pragma omp parallel for
   for (uint pencili = 0; pencili < nPencils; ++pencili) {
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
   //phiprof::Timer gpuPrefetchTimer {"trans-amr-prefetch"};
   allPencilsMeshes->optimizeGPU(bgStream);
   allPencilsContainers->optimizeGPU(bgStream);
   //gpuPrefetchTimer.stop();

   // Extract pointers to data in managed memory
   uint* pencilLengths = DimensionPencils[dimension].gpu_lengthOfPencils->data();
   uint* pencilStarts = DimensionPencils[dimension].gpu_idsStart->data();
   Realf* pencilDZ = DimensionPencils[dimension].gpu_sourceDZ->data();
   Realf* pencilRatios = DimensionPencils[dimension].gpu_targetRatios->data();
   gatherPointerTimer.stop();

   //gpuPrefetchTimer.start();
   //CHK_ERR( gpuStreamSynchronize(bgStream) );
   //gpuPrefetchTimer.stop();

   allocateTimer.start();
   const vmesh::LocalID unionOfBlocksSetSize = unionOfBlocksSet->size();
   gpu_trans_allocate(0,0,0,unionOfBlocksSetSize);
   allocateTimer.stop();

   phiprof::Timer buildTimer2 {"trans-amr-buildBlockList-2"};
   // Now we ensure the union of blocks gathering is complete and extract the union of blocks into a vector
   const uint nAllBlocks = unionOfBlocksSet->extractAllKeys<false>(*unionOfBlocks,bgStream);
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

   // How many blocks worth of pre-allocated buffer do we have for each thread?
   const uint maxThreads = gpu_getMaxThreads();
   const uint currentAllocation = gpu_vlasov_getSmallestAllocation() * TRANSLATION_BUFFER_ALLOCATION_FACTOR;
   // GPUTODO: make tmep buffer allocation a config parameter? Current approach is not necessarily
   // good if average block count vs grid size are mismatcheds

   // How many block GIDs could each thread manage in parallel with this existing temp buffer? // floor int division
   const uint nBlocksPerThread = currentAllocation / sumOfLengths;

   // And how many block GIDs will we actually manage?
   const uint totalPerThread =  1 + ((nAllBlocks - 1) / maxThreads); // ceil int division, no more than this
   const uint nGpuBlocks  = nBlocksPerThread  < totalPerThread ? nBlocksPerThread : totalPerThread;

   #pragma omp parallel
   {
      // Thread id used for persistent device memory pointers.
      const uint cpuThreadID = gpu_getThread();
      gpuStream_t stream = gpu_getStream();

      phiprof::Timer bufferTimer {bufferId};
      Vec* pencilOrderedSource = gpu_blockDataOrdered[cpuThreadID];

      // Two temporary buffers, used in-kernel for both reading and writing
      // GPUTODO: pre-allocate one per thread, here verify sufficient size
      Realf** pencilBlockData; // Array of pointers into actual block data
      uint* pencilBlocksCount; // Array of counters if pencil needs to be propagated for this block or not
      CHK_ERR( gpuMallocAsync((void**)&pencilBlockData, sumOfLengths*nGpuBlocks*sizeof(Realf*), stream) );
      CHK_ERR( gpuMallocAsync((void**)&pencilBlocksCount, nPencils*nGpuBlocks*sizeof(uint), stream) );
      CHK_ERR( gpuStreamSynchronize(stream) );
      bufferTimer.stop();

      // Loop over velocity space blocks (threaded, each thread using one stream, and
      // each stream launching a kernel running several blocks in parallel, but not using a for-loop)
      phiprof::Timer mappingTimer {mappingId}; // mapping (top-level)
      const uint startingBlockIndex = cpuThreadID*nGpuBlocks;
      const uint blockIndexIncrement = maxThreads*nGpuBlocks;
      // Each thread, using its own stream, will launch nGpuBlocks instances of the below kernel, where each instance
      // propagates all pencils for the block in question.
      dim3 block(WID,WID,WID); // assumes VECL==WID2
      translation_kernel<<<nGpuBlocks, block, 0, stream>>> (
         dimension,
         dt,
         pencilLengths,
         pencilStarts,
         allBlocks, // List of all block GIDs
         nAllBlocks, // size of list of block GIDs which we won't exceed
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
__global__ static void remote_increment_kernel (
   Realf* blockData,
   Realf* neighborData,
   vmesh::LocalID nBlocks
   ) {
   //const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x; // ==LID
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;
   // Increment value
   // atomicAdd(&blockData[blocki * WID3 + ti],neighborData[blocki * WID3 + ti]);
   // As each target block has its own GPU stream, we ensure that we don't write concurrently from
   // several different threads or kernels, and thus don't need to use atomic operations.
   blockData[blocki * WID3 + ti] += neighborData[blocki * WID3 + ti];
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

   // Fast return if no cells to process
   int mpiProcs;
   MPI_Comm_size(MPI_COMM_WORLD,&mpiProcs);
   if (mpiProcs == 1) {
      return;
   }

   // GPUTODO: First attempts at using managed memory for remote neighbours
   // Should move to re-using managed memory buffers and ensuring size is suitable?
   // If that path is taken, it should also check for any local cells *not* on process
   // boundary and free the buffers from those cells.
   int device = gpu_getDevice();

   int neighborhood = 0;
   //normalize and set neighborhoods
   if(direction > 0) {
      direction = 1;
      switch (dimension) {
      case 0:
         neighborhood = SHIFT_P_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = SHIFT_P_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = SHIFT_P_Z_NEIGHBORHOOD_ID;
         break;
      }
   }
   if(direction < 0) {
      direction = -1;
      switch (dimension) {
      case 0:
         neighborhood = SHIFT_M_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = SHIFT_M_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = SHIFT_M_Z_NEIGHBORHOOD_ID;
         break;
      }
   }

   //const vector<CellID>& local_cells = getLocalCells();
   const vector<CellID> local_cells = mpiGrid.get_local_cells_on_process_boundary(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   const vector<CellID> remote_cells = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_NEIGHBORHOOD_ID);

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
      for (const auto& [neighbor, dir] : mpiGrid.get_face_neighbors_of(c)) {
         if(dir == ((int)dimension + 1) * direction) {
            p_nbrs.push_back(neighbor);
         }
         if(dir == -1 * ((int)dimension + 1) * direction) {
            n_nbrs.push_back(neighbor);
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
                        ccell->neighbor_block_data.at(sendIndex) = 0;
                        sendBuffers.push_back(0);
                     } else {
                        // GPUTODO: This is now unified memory. With GPU-aware MPI it could be on-device.
                        CHK_ERR( gpuMallocManaged((void**)&ccell->neighbor_block_data.at(sendIndex), ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf)) );
                        // CHK_ERR( gpuMemPrefetchAsync(ccell->neighbor_block_data.at(sendIndex),ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf),device,0) );
                        CHK_ERR( gpuMemset(ccell->neighbor_block_data.at(sendIndex), 0, ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf)) );
                        sendBuffers.push_back(ccell->neighbor_block_data.at(sendIndex));
                     }
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
                     receiveBuffers.push_back(0);
                  } else {
                     // GPUTODO: This is now unified memory. With GPU-aware MPI it could be on-device.
                     CHK_ERR( gpuMallocManaged((void**)&ncell->neighbor_block_data.at(recvIndex), ncell->neighbor_number_of_blocks.at(recvIndex) * WID3 * sizeof(Realf)) );
                     CHK_ERR( gpuMemPrefetchAsync(ncell->neighbor_block_data.at(recvIndex), ncell->neighbor_number_of_blocks.at(recvIndex) * WID3 * sizeof(Realf), device,0) );
                     receiveBuffers.push_back(ncell->neighbor_block_data.at(recvIndex));
                  }
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
                           receiveBuffers.push_back(0);
                        } else {
                           // GPUTODO: This is now unified memory. With GPU-aware MPI it could be on-device.
                           CHK_ERR( gpuMallocManaged((void**)&ncell->neighbor_block_data.at(i_sib), ncell->neighbor_number_of_blocks.at(i_sib) * WID3 * sizeof(Realf)) );
                           CHK_ERR( gpuMemPrefetchAsync(ncell->neighbor_block_data.at(i_sib), ncell->neighbor_number_of_blocks.at(i_sib) * WID3 * sizeof(Realf), device,0) );
                           receiveBuffers.push_back(ncell->neighbor_block_data.at(i_sib));
                        }
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
   if (receive_cells.size() != 0) {
      phiprof::Timer updateRemoteTimerIncrement {"trans-amr-remotes-increment"};
      for (size_t c = 0; c < receive_cells.size(); ++c) {
         SpatialCell* receive_cell = mpiGrid[receive_cells[c]];
         SpatialCell* origin_cell = mpiGrid[receive_origin_cells[c]];
         if (!receive_cell || !origin_cell) {
            continue;
         }

         Realf *blockData = receive_cell->get_data(popID);
         Realf *neighborData = origin_cell->neighbor_block_data[receive_origin_index[c]];
         vmesh::LocalID nBlocks = receive_cell->get_number_of_velocity_blocks(popID);
         const uint maxThreads = gpu_getMaxThreads();
         // Increment needs to be parallel-safe, so use modulo of cellid as stream number
         gpuStream_t cellStream = gpuStreamList[receive_cells[c] % maxThreads];
         if (nBlocks>0) {
            dim3 block(WID,WID,WID);
            remote_increment_kernel<<<nBlocks, block, 0, cellStream>>> (
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
