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

#include "gpu_1d_ppm_nonuniform.hpp"
//#include "gpu_1d_ppm_nonuniform_conserving.hpp"

#include "gpu_trans_map_amr.hpp"
#include "cpu_trans_pencils.hpp"
#include "../arch/gpu_base.hpp"

#ifdef USE_WARPACCESSORS
 #define USE_TRANS_WARPACCESSORS
#endif

// Skip remapping for this stencil, if no blocks exist
__device__ inline bool check_skip_blocks(const Realf* __restrict__ const *pencilBlockData, const uint centerOffset) {
   for (int index=-VLASOV_STENCIL_WIDTH; index<VLASOV_STENCIL_WIDTH+1; ++index) {
      if (pencilBlockData[centerOffset + index] != NULL) {
         return false;
      }
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
 * @param  nPencils  Number of total pencils (constant)
 * @param  sumOfLengths  sum of all pencil lengths (constant)
 * @param  threshold  sparsity threshold, used by slope limiters
 * @param  dev_allPencilsMeshes  Pointer to buffer of pointers to velocity meshes
 * @param  dev_allPencilsContainers  Pointer to buffer of pointers to BlockContainers
 * @param  pencilBlockData  Pointer to buffer of pointers into cell block data, both written and read
 * @param  dev_blockDataOrdered  Pointer to aligned buffer used as interim values
 * @param  pencilDZ  Pointer into buffer of pencil cell sizes
 * @param  pencilRatios  Pointer into buffer with pencil target ratios (due to AMR)
 * @param  pencilBlocksCount  Pointer into buffer for storing how many non-empty blocks each pencil has for current GID
 */

// GPUTODO: The translation kernel may need splitting up into one read/prep kernel and another translate/write kernel,
// so that pointers to each part can be declared const __restrict__ in turn. A quick attempt at this
// was actually slower, and e.g. Realf** pencilBlockData could not be const __restrict__ anyway. Using
// const_cast in some loops resulted in data corruption.

//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor, maxBlocksPerCluster)
__global__ void __launch_bounds__(WID3) translation_kernel(
   const uint dimension,
   const Realf dt,
   const uint* __restrict__ pencilLengths,
   const uint* __restrict__ pencilStarts,
   const vmesh::GlobalID* __restrict__ allBlocks, // List of all blocks
   const uint nAllBlocks, // size of list of blocks which we won't exceed
   const uint nPencils, // Number of total pencils (constant)
   const uint sumOfLengths, // sum of all pencil lengths (constant)
   const Realf threshold, // used by slope limiters
   const vmesh::VelocityMesh* __restrict__ const *dev_allPencilsMeshes, // Pointers to velocity meshes
   vmesh::VelocityBlockContainer* *dev_allPencilsContainers, // pointers to BlockContainers
   Realf** pencilBlockData, // pointers into cell block data, both written and read
   Realf** dev_blockDataOrdered, // buffer of pointers to mapping input data
   const Realf* __restrict__ pencilDZ,
   const Realf* __restrict__ pencilRatios, // buffer holding target ratios
   uint* pencilBlocksCount, // store how many non-empty blocks each pencil has for this GID
   uint *dev_pencilsInBin,
   uint *dev_binStart,
   uint *dev_binSize,
   const uint numberOfBins
   ) {
   // This is launched with grid size (nGpuBlocks,nAllocations,1)
   // where nGpuBlocks is the count of blocks which fit in the smallest temp buffer at once
   // and nAllocations is the number of temp GPU buffers to use.
   const uint startingBlockIndex = blockIdx.y*gridDim.x;
   const uint blockIndexIncrement = gridDim.y*gridDim.x;
   Realf* pencilOrderedSource = dev_blockDataOrdered[blockIdx.y];

   // This is launched with block size (WID,WID,WID)
   const vmesh::LocalID ti = (threadIdx.z)*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

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
   const uint pencilBlockDataOffset = (blockIdx.x * sumOfLengths) + (blockIdx.y * sumOfLengths * gridDim.x);
   const uint pencilOrderedSourceOffset = blockIdx.x * sumOfLengths * WID3;
   const uint pencilBlocksCountOffset = (blockIdx.x * nPencils) + (blockIdx.y * nPencils * gridDim.x);
   const vmesh::VelocityMesh* __restrict__ randovmesh = dev_allPencilsMeshes[0]; // just some vmesh
   const Realf dvz = randovmesh->getCellSize()[dimension];
   const Realf vz_min = randovmesh->getMeshMinLimits()[dimension];

   // Acting on velocity block blockGID, now found from array
   for (uint thisBlockIndex = startingBlockIndex + blockIdx.x; thisBlockIndex < nAllBlocks; thisBlockIndex += blockIndexIncrement) {
      
      const uint blockGID = allBlocks[thisBlockIndex];
      // First read data in
      uint nBin = blockIdx.z;

      for (uint pencilIndex = 0; pencilIndex < dev_binSize[nBin]; pencilIndex++) {
         const uint pencili = dev_pencilsInBin[dev_binStart[nBin]+pencilIndex];
         const uint lengthOfPencil = pencilLengths[pencili];
         const uint start = pencilStarts[pencili];
         // Get pointer to temprary buffer of VEC-ordered data for this kernel
         Realf* thisPencilOrderedSource = pencilOrderedSource + pencilOrderedSourceOffset + start * WID3;
         uint nonEmptyBlocks = 0;
         // Go over pencil length, gather cellblock data into pencil source data
         for (uint celli = 0; celli < lengthOfPencil; celli++) {
            const vmesh::VelocityMesh* __restrict__ vmesh = dev_allPencilsMeshes[start + celli];
            vmesh::VelocityBlockContainer* cellContainer = dev_allPencilsContainers[start + celli];
            // Now using warp accessor.
            #ifdef USE_TRANS_WARPACCESSORS
            const vmesh::LocalID blockLID = vmesh->warpGetLocalID(blockGID,ti);
            #else
            const vmesh::LocalID blockLID = vmesh->getLocalID(blockGID);
            #endif
            // Store block data pointer for both loading of data and writing back to the cell
            if (blockLID != vmesh->invalidLocalID()) {
               #ifdef DEBUG_VLASIATOR
               const vmesh::LocalID meshSize = vmesh->size();
               const vmesh::LocalID VBCSize = cellContainer->size();
               if ((blockLID>=meshSize) || (blockLID>=VBCSize)) {
                  if (ti==0) {
                     printf("Error in translation: trying to access LID %ul but sizes are vmesh %ul VBC %ul\n",blockLID,meshSize,VBCSize);
                  }
               }
               #endif
               if (ti==0) {
                  pencilBlockData[pencilBlockDataOffset + start + celli] = cellContainer->getData(blockLID);
                  nonEmptyBlocks++;
               }
               // Valid block, store values in contiguous data
               thisPencilOrderedSource[celli * WID3 + ti]
                  = (cellContainer->getData(blockLID))[ti];
            } else {
               if (ti==0) {
                  pencilBlockData[pencilBlockDataOffset + start + celli] = NULL;
               }
               // Non-existing block, push in zeroes
               thisPencilOrderedSource[celli * WID3 + ti] = (Realf)(0.0);
            }
            __syncthreads();
         } // End loop over this pencil
         if (ti==0) {
            pencilBlocksCount[pencilBlocksCountOffset + pencili] = nonEmptyBlocks;
         }
         __syncthreads();
      } // end loop over pencils in this bin

      __syncthreads();
      // Now we reset target blocks
      for (uint pencilIndex = 0; pencilIndex < dev_binSize[nBin]; pencilIndex++) {
         const uint pencili = dev_pencilsInBin[dev_binStart[nBin]+pencilIndex];
         const uint lengthOfPencil = pencilLengths[pencili];
         const uint start = pencilStarts[pencili];
         for (uint celli = 0; celli < lengthOfPencil; celli++) {
            if (pencilRatios[start + celli] != 0) {
               // Is a target cell, needs to be reset
               if (pencilBlockData[pencilBlockDataOffset + start + celli]) {
                  (pencilBlockData[pencilBlockDataOffset + start + celli])[ti] = (Realf)(0.0);
               }
            }
         } // end loop over this pencil
      } // end loop over pencils in this bin

      __syncthreads();

      // Now we propagate the pencils and write data back to the block data containers
      // Get velocity data from vmesh that we need later to calculate the translation
   
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
      for (uint pencilIndex = 0; pencilIndex < dev_binSize[nBin]; pencilIndex++) {
         const uint pencili = dev_pencilsInBin[dev_binStart[nBin]+pencilIndex];
         if (pencilBlocksCount[pencilBlocksCountOffset + pencili] == 0) {
            continue;
         }
         const uint lengthOfPencil = pencilLengths[pencili];
         const uint start = pencilStarts[pencili];
         const Realf* __restrict__ thisPencilOrderedSource = pencilOrderedSource + pencilOrderedSourceOffset + start * WID3;

         // Go over length of propagated cells
         for (uint i = VLASOV_STENCIL_WIDTH; i < lengthOfPencil-VLASOV_STENCIL_WIDTH; i++){
            // Get pointers to block data used for output.
            Realf* block_data_m1 = pencilBlockData[pencilBlockDataOffset + start + i - 1];
            Realf* block_data =    pencilBlockData[pencilBlockDataOffset + start + i];
            Realf* block_data_p1 = pencilBlockData[pencilBlockDataOffset + start + i + 1];

            // Cells which shouldn't be written to (e.g. sysboundary cells) have a targetRatio of 0
            // Also need to check if pointer is valid, because a cell can be missing an elsewhere propagated block
            const Realf areaRatio_m1 = pencilRatios[start + i - 1];
            const Realf areaRatio =    pencilRatios[start + i];
            const Realf areaRatio_p1 = pencilRatios[start + i + 1];

            // (no longer loop over) planes (threadIdx.z) and vectors within planes (just 1 by construction)
            const Realf cell_vz = (blockIndicesD * WID + vz_index + (Realf)(0.5)) * dvz + vz_min; //cell centered velocity
            const Realf z_translation = cell_vz * dt / pencilDZ[start + i]; // how much it moved in time dt (reduced units)

            // Determine direction of translation
            // part of density goes here (cell index change along spatial direcion)
            const bool positiveTranslationDirection = (z_translation > (Realf)(0.0));

            // Calculate normalized coordinates in current cell.
            // The coordinates (scaled units from 0 to 1) between which we will
            // integrate to put mass in the target  neighboring cell.
            // Normalize the coordinates to the origin cell. Then we scale with the difference
            // in volume between target and origin later when adding the integrated value.
            Realf z_1,z_2;
            z_1 = positiveTranslationDirection ? (Realf)(1.0) - z_translation : (Realf)(0.0);
            z_2 = positiveTranslationDirection ? (Realf)(1.0) : - z_translation;

            #ifdef DEBUG_VLASIATOR
            if ( abs(z_1) > (Realf)(1.0) || abs(z_2) > (Realf)(1.0) ) {
               assert( 0 && "Error in translation, CFL condition violated.");
            }
            #endif

            // If no blocks exist for this mapping, skip forward.
            if (!check_skip_blocks(pencilBlockData,pencilBlockDataOffset + start + i)) {
               // Compute polynomial coefficients
               Realf a[3];
               // Silly indexing into coefficient calculation necessary due to
               // built-in assumptions of unsigned indexing.
               compute_ppm_coeff_nonuniform(pencilDZ + start + i - VLASOV_STENCIL_WIDTH,
                                          thisPencilOrderedSource + (i - VLASOV_STENCIL_WIDTH) * WID3,
                                          h4, VLASOV_STENCIL_WIDTH, a, threshold, ti, WID3);
               // Compute integral
               const Realf ngbr_target_density =
                  z_2 * ( a[0] + z_2 * ( a[1] + z_2 * a[2] ) ) -
                  z_1 * ( a[0] + z_1 * ( a[1] + z_1 * a[2] ) );

               // Store mapped density in two target cells
               // in the current original cells we will put the rest of the original density
               // Now because each GPU block handles all pencils for an unique GID, we shouldn't need atomic additions here.

               // NOTE: not using atomic operations causes huge diffs (as if self contribution was neglected)! 11.01.2024 MB
               if (areaRatio && block_data) {
                  const Realf selfContribution = (thisPencilOrderedSource[i * WID3 + ti] - ngbr_target_density) * areaRatio;
                  //atomicAdd(&block_data[ti],selfContribution);
                  block_data[ti] += selfContribution;
               }
               if (areaRatio_p1 && block_data_p1) {
                  const Realf p1Contribution = (positiveTranslationDirection ? ngbr_target_density
                                                * pencilDZ[start + i] / pencilDZ[start + i + 1] : (Realf)(0.0)) * areaRatio_p1;
                  //atomicAdd(&block_data_p1[ti],p1Contribution);
                  block_data_p1[ti] += p1Contribution;
               }
               if (areaRatio_m1 && block_data_m1) {
                  const Realf m1Contribution = (!positiveTranslationDirection ? ngbr_target_density
                                                * pencilDZ[start + i] / pencilDZ[start + i - 1] : (Realf)(0.0)) * areaRatio_m1;
                  //atomicAdd(&block_data_m1[ti],m1Contribution);
                  block_data_m1[ti] += m1Contribution;
               }
            } // Did not skip remapping
            __syncthreads();
         } // end loop over this pencil
      } // end loop over pencils in this bin
      __syncthreads();
   } // end loop over blocks
}

/* Mini-kernel for looping over all available velocity meshes and gathering
 * the union of all existing blocks.
 *
 * @param unionOfBlocksSet Hashmap, where keys are those blocks which are in the union of all blocks
 * @param allVmeshPointer Buffer of pointers to velocitymeshes, used for gathering active blocks
 * @param nAllCells count of cells to read from allVmeshPointer
 */
#ifdef USE_WARPACCESSORS
__global__ void __launch_bounds__(GPUTHREADS*WARPSPERBLOCK) gather_union_of_blocks_kernel_WA(
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet,
   const vmesh::VelocityMesh* __restrict__ const *dev_vmeshes,
   const uint nAllCells)
{
   const int ti = threadIdx.x; // [0,GPUTHREADS)
   const int indexInBlock = threadIdx.y; // [0,WARPSPERBLOCK)
   const uint cellIndex = blockIdx.x;
   const uint blockIndexBase = blockIdx.y * WARPSPERBLOCK;
   const vmesh::VelocityMesh* __restrict__ thisVmesh = dev_vmeshes[cellIndex];
   const uint thisVmeshSize = thisVmesh->size();
   const uint blockIndex = blockIndexBase + indexInBlock;
   if (blockIndex < thisVmeshSize) {
      // Now with warp accessors
      const vmesh::GlobalID GID = thisVmesh->getGlobalID(blockIndex);
      // warpInsert<true> only inserts if key does not yet exist
      unionOfBlocksSet->warpInsert<true>(GID, (vmesh::LocalID)GID, ti);
   }
}
#else
__global__ void __launch_bounds__(GPUTHREADS*WARPSPERBLOCK) gather_union_of_blocks_kernel(
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet,
   const vmesh::VelocityMesh* __restrict__ const *dev_vmeshes,
   const uint nAllCells)
{
   const int indexInBlock = threadIdx.x; // [0,WARPSPERBLOCK*GPUTHREADS)
   const uint cellIndex = blockIdx.x;
   const uint blockIndexBase = blockIdx.y * WARPSPERBLOCK * GPUTHREADS;
   const vmesh::VelocityMesh* __restrict__ thisVmesh = dev_vmeshes[cellIndex];
   const uint thisVmeshSize = thisVmesh->size();
   const uint blockIndex = blockIndexBase + indexInBlock;
   if (blockIndex < thisVmeshSize) {
      const vmesh::GlobalID GID = thisVmesh->getGlobalID(blockIndex);
      // <true> only inserts if key does not yet exist
      unionOfBlocksSet->set_element<true>(GID, (vmesh::LocalID)GID);
   }
}
#endif

/* Map velocity blocks in all local cells forward by one time step in one spatial dimension.
 * This function uses 1-cell wide pencils to update cells in-place to avoid allocating large
 * temporary buffers.
 *
 * @param [in] mpiGrid DCCRG grid object
 * @param [in] localPropagatedCells List of local cells that get propagated
 * ie. not boundary or DO_NOT_COMPUTE
 * @param [in] remoteTargetCells List of non-local target cells
 * @param [in] nPencilsLB vector where the number of active pencils for each cell can be stored, to be used by load balance
 * @param [in] dimension Spatial dimension
 * @param [in] dt Time step
 * @param [in] popId Particle population ID
 */
bool trans_map_1d_amr(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                      const vector<CellID>& localPropagatedCells,
                      const vector<CellID>& remoteTargetCells,
                      std::vector<uint>& nPencilsLB,
                      const uint dimension,
                      const Realf dt,
                      const uint popID) {

   phiprof::Timer setupTimer {"trans-amr-setup"};

   // return if there's no cells to propagate
   if(localPropagatedCells.size() == 0) {
      return false;
   }
   gpuStream_t bgStream = gpu_getStream(); // uses stream assigned to thread 0, not the blocking default stream

   // Vector with all cell ids
   vector<CellID> allCells(localPropagatedCells);
   allCells.insert(allCells.end(), remoteTargetCells.begin(), remoteTargetCells.end());
   const uint nAllCells = allCells.size();

   phiprof::Timer allocateTimer {"trans-amr-allocs"};
   // Ensure GPU data has sufficient allocations/sizes
   const uint sumOfLengths = DimensionPencils[dimension].sumOfLengths;
   gpu_vlasov_allocate(sumOfLengths);
   // Ensure allocation for allPencilsMeshes, allPencilsContainers
   gpu_trans_allocate(nAllCells,sumOfLengths,0,0,0,0);
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
         host_vmeshes[celli] = mpiGrid[allCells[celli]]->dev_get_velocity_mesh(popID); // GPU-side vmesh
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
   // Copy vmesh pointers to GPU
   CHK_ERR( gpuMemcpy(dev_vmeshes, host_vmeshes, nAllCells*sizeof(vmesh::VelocityMesh*), gpuMemcpyHostToDevice) );
   // Reserve size for unionOfBlocksSet
   gpu_trans_allocate(0,0,largestFoundMeshSize,0,0,0);
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
   // propagated cells. We could do host-side pointer
   // gathering in parallel with it, but that leads to potentially incorrect phiprof output..
#ifdef USE_WARPACCESSORS
   const uint maxBlocksPerCell =  1 + ((largestFoundMeshSize - 1) / WARPSPERBLOCK); // ceil int division
   dim3 gatherdims_blocks(nAllCells,maxBlocksPerCell,1);
   dim3 gatherdims_threads(GPUTHREADS,WARPSPERBLOCK,1);
   gather_union_of_blocks_kernel_WA<<<gatherdims_blocks, gatherdims_threads, 0, bgStream>>> (
#else
   const uint maxBlocksPerCell =  1 + ((largestFoundMeshSize - 1) / (WARPSPERBLOCK*GPUTHREADS)); // ceil int division
   dim3 gatherdims_blocks(nAllCells,maxBlocksPerCell,1);
   dim3 gatherdims_threads(GPUTHREADS*WARPSPERBLOCK,1,1);
   gather_union_of_blocks_kernel<<<gatherdims_blocks, gatherdims_threads, 0, bgStream>>> (
#endif
      dev_unionOfBlocksSet,
      dev_vmeshes,
      nAllCells
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(bgStream) ); // So we get phiprof data of this gathering time
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
         host_allPencilsMeshes[start+i] = mpiGrid[thisCell]->dev_get_velocity_mesh(popID);
         host_allPencilsContainers[start+i] = mpiGrid[thisCell]->dev_get_velocity_blocks(popID);
      }
   }
   // Copy pencil meshes and VBCs to GPU
   CHK_ERR( gpuMemcpy(dev_allPencilsMeshes, host_allPencilsMeshes, sumOfLengths*sizeof(vmesh::VelocityMesh*), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuMemcpy(dev_allPencilsContainers, host_allPencilsContainers, sumOfLengths*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice) );

   // Extract pointers to data in unified memory
   uint* pencilLengths = DimensionPencils[dimension].gpu_lengthOfPencils;
   uint* pencilStarts = DimensionPencils[dimension].gpu_idsStart;
   Realf* pencilDZ = DimensionPencils[dimension].gpu_sourceDZ;
   Realf* pencilRatios = DimensionPencils[dimension].gpu_targetRatios;
   gatherPointerTimer.stop();

   // Now we ensure the union of blocks gathering is complete and find the size of it. Use it to ensure allocations.
   allocateTimer.start();
   // Use non-pagefaulting fetching of metadata
   Hashinator::Info mapInfo;
   unionOfBlocksSet->copyMetadata(&mapInfo, bgStream);
   CHK_ERR( gpuStreamSynchronize(bgStream) );
   const vmesh::LocalID unionOfBlocksSetSize = mapInfo.fill;
   gpu_trans_allocate(0,0,0,unionOfBlocksSetSize,0,0);
   allocateTimer.stop();

   phiprof::Timer buildTimer2 {"trans-amr-buildBlockList-2"};
   // Extract the union of blocks into a vector
   unionOfBlocksSet->extractAllKeysLoop(*dev_unionOfBlocks,bgStream);
   split::SplitInfo unionInfo;
   unionOfBlocks->copyMetadata(&unionInfo, bgStream);
   CHK_ERR( gpuStreamSynchronize(bgStream) );
   const uint nAllBlocks = unionInfo.size;
   const uint numberOfBins = DimensionPencils[dimension].activeBins.size();
   vmesh::GlobalID *allBlocks = unionOfBlocks->data();
   // This threshold value is used by slope limiters.
   Realf threshold = mpiGrid[DimensionPencils[dimension].ids[VLASOV_STENCIL_WIDTH]]->getVelocityBlockMinValue(popID);
   buildTimer2.stop();

   // GPUTODO: Improve config parameter use for temp buffer allocation, consolidate buffers, etc.
   // Current approach is not necessarily good if average block count vs grid size are mismatched.
   // Best would be to have one large buffer from which sub-buffers are provided to the various Vlasov
   // solvers as necessary.

   // How many blocks worth of pre-allocated buffer do we have for each thread?
   const uint currentAllocation = gpu_vlasov_getSmallestAllocation();
   // How many block GIDs could each thread manage in parallel with this existing temp buffer? // floor int division
   // Note: we no longer launch from several threads, but some buffers are still identified via threads.
   const uint nBlocksPerAllocation = currentAllocation / sumOfLengths;
   // And how many block GIDs will we actually manage at once?
   const uint numAllocations = gpu_getAllocationCount();
   const uint totalPerAllocation =  1 + ((nAllBlocks - 1) / numAllocations); // ceil int division
   // no more than this per allocation
   const uint nGpuBlocks = std::min(nBlocksPerAllocation,totalPerAllocation);
   // Limit is either how many blocks exist, or how many fit in buffer.

   phiprof::Timer bufferTimer {"trans-amr-buffers"};
   // Two temporary buffers, used in-kernel for both reading and writing
   // (dev_pencilBlockData and dev_pencilBlocksCount)
   allocateTimer.start();
   gpu_trans_allocate(0,sumOfLengths,0,0,nGpuBlocks,nPencils);
   allocateTimer.stop();
   bufferTimer.stop();

   /***********************/
   setupTimer.stop();
   /***********************/

   // Loop over velocity space blocks
   phiprof::Timer mappingTimer {"trans-amr-mapping"};
   // Launch 2D grid: First dimension is how many blocks fit in one temp buffer, second one
   // is which temp buffer allocation index to use. (GPUTODO: simplify together with buffer consolidation)
   dim3 grid(nGpuBlocks,numAllocations,numberOfBins);
   dim3 block(WID,WID,WID);
   translation_kernel<<<grid, block, 0, bgStream>>> (
      dimension,
      dt,
      pencilLengths,
      pencilStarts,
      allBlocks, // List of all block GIDs
      nAllBlocks, // size of list of block GIDs which we won't exceed
      nPencils, // Number of total pencils (constant)
      sumOfLengths, // sum of all pencil lengths (constant)
      threshold,
      dev_allPencilsMeshes, // Pointers to velocity meshes
      dev_allPencilsContainers, // pointers to BlockContainers
      dev_pencilBlockData, // pointers into cell block data, both written and read
      dev_blockDataOrdered, // buffer of pointers to ordered buffer data
      pencilDZ,
      pencilRatios, // buffer tor holding target ratios
      dev_pencilBlocksCount, // store how many non-empty blocks each pencil has for this GID
      gpuMemoryManager.getPointer<uint>(DimensionPencils[dimension].dev_pencilsInBin),
      gpuMemoryManager.getPointer<uint>(DimensionPencils[dimension].dev_binStart),
      gpuMemoryManager.getPointer<uint>(DimensionPencils[dimension].dev_binSize),
      numberOfBins
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(bgStream) );
   mappingTimer.stop();

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
void update_remote_mapping_contribution_amr(
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

   /**
      GPUTODO: First attempts at using unified memory for remote neighbours
      Should probably:
      a) Switch to using pure device buffers. Some MPI-CUDA implementations do not work well with UB buffers.
      b) transition to re-using buffers and ensuring size is suitable?
      If that path is taken, it should also check for any local cells *not* on process
      boundary (anymore, due to LB) and free the buffers from those cells... so gets tricky.
   */
   int device = gpu_getDevice();

   int neighborhood = 0;
   //normalize and set neighborhoods
   if(direction > 0) {
      direction = 1;
      switch (dimension) {
      case 0:
         neighborhood = Neighborhoods::SHIFT_P_X;
         break;
      case 1:
         neighborhood = Neighborhoods::SHIFT_P_Y;
         break;
      case 2:
         neighborhood = Neighborhoods::SHIFT_P_Z;
         break;
      }
   }
   if(direction < 0) {
      direction = -1;
      switch (dimension) {
      case 0:
         neighborhood = Neighborhoods::SHIFT_M_X;
         break;
      case 1:
         neighborhood = Neighborhoods::SHIFT_M_Y;
         break;
      case 2:
         neighborhood = Neighborhoods::SHIFT_M_Z;
         break;
      }
   }

   //const vector<CellID>& local_cells = getLocalCells();
   const vector<CellID> local_cells = mpiGrid.get_local_cells_on_process_boundary(Neighborhoods::VLASOV_SOLVER);
   const vector<CellID> remote_cells = mpiGrid.get_remote_cells_on_process_boundary(Neighborhoods::VLASOV_SOLVER);

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
      if (!ccell) {
         continue;
      }
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
      // one cell is the neighbor on both + and - side to the same process
      vector<CellID> send_cells_vector(send_cells.begin(), send_cells.end());
      for (uint c = 0; c < send_cells_vector.size(); c++) {
         SpatialCell* send_cell = mpiGrid[send_cells_vector[c]];
         gpuStream_t stream = gpu_getStream();
         Realf* blockData = send_cell->get_data(popID);
         CHK_ERR( gpuMemsetAsync(blockData, 0, WID3*send_cell->get_number_of_velocity_blocks(popID)*sizeof(Realf),stream) );
      }
      CHK_ERR( gpuDeviceSynchronize() );
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
