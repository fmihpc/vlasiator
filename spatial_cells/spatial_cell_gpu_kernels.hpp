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

#ifndef VLASIATOR_SPATIAL_CELL_KERNELS_HPP
#define VLASIATOR_SPATIAL_CELL_KERNELS_HPP

#ifdef USE_WARPACCESSORS
 #define USE_BATCH_WARPACCESSORS
#endif

// GPUTODO: Make error-checking functions to be called inside kernels instead of duplicating so much code.
// Same for gathering mass loss.

/*! Note: these kernels are called only when adjusting blocks for a single spatial cell. Usually, these
  are performed as batch operations, found in block_adjust_gpu.cpp and block_adjust_gpu_kernels.hpp
*/

/** GPU kernel for identifying which blocks have relevant content */
__global__ void __launch_bounds__(WID3,4) update_velocity_block_content_lists_kernel (
   const vmesh::VelocityMesh* __restrict__ vmesh,
   const vmesh::VelocityBlockContainer* __restrict__ blockContainer,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwcl_map,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwncl_map,
   const Real velocity_block_min_value
   ) {
   // const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const uint ti = threadIdx.x;

   // Each GPU block / workunit can manage several Vlasiator velocity blocks at once.
   const uint vlasiBlocksPerWorkUnit = 1;
   const uint workUnitIndex = 0; // [0,vlasiBlocksPerWorkUnit)
   // const uint vlasiBlocksPerWorkUnit = WARPSPERBLOCK * GPUTHREADS / WID3;
   // const uint workUnitIndex = ti / WID3; // [0,vlasiBlocksPerWorkUnit)
   const uint b_tid = ti % WID3; // [0,WID3)
   const uint blockLID = blocki * vlasiBlocksPerWorkUnit + workUnitIndex; // [0,nBlocksToChange)

   __shared__ int has_content[WARPSPERBLOCK * GPUTHREADS];
   // maps can only be cleared from host
   const uint nBlocks = vmesh->size();
   if (blockLID < nBlocks) {
      const vmesh::GlobalID blockGID = vmesh->getGlobalID(blockLID);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == vmesh->invalidGlobalID()) {
         if (b_tid==0) {
            printf("Invalid GID encountered in update_velocity_block_content_lists_kernel!\n");
         }
         return;
      }
      if (blockLID == vmesh->invalidLocalID()) {
         if (b_tid==0) {
            printf("Invalid LID encountered in update_velocity_block_content_lists_kernel!\n");
         }
         return;
      }
      #endif
      // Check each velocity cell if it is above the threshold
      const Realf* avgs = blockContainer->getData(blockLID);
      has_content[ti] = avgs[b_tid] >= velocity_block_min_value ? 1 : 0;
      __syncthreads(); // THIS SYNC IS CRUCIAL!
      // Implemented just a simple non-optimized thread OR
      // GPUTODO reductions via warp voting

      // Perform loop only until first value fulfills condition
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (has_content[0]) {
            break;
         }
         if (b_tid < s) {
            has_content[ti] = has_content[ti] || has_content[ti + s];
         }
         __syncthreads();
      }
      #ifdef USE_WARPACCESSORS
      // Insert into map only from threads 0...WARPSIZE
      if (b_tid < GPUTHREADS) {
         if (has_content[0]) {
            vbwcl_map->warpInsert(blockGID,blockLID,b_tid);
         } else {
            vbwncl_map->warpInsert(blockGID,blockLID,b_tid);
         }
      }
      #else
      // Insert into map only from thread 0
      if (b_tid == 0) {
         if (has_content[0]) {
            vbwcl_map->set_element(blockGID,blockLID);
         } else {
            vbwncl_map->set_element(blockGID,blockLID);
         }
      }
      #endif
      __syncthreads();
   }
}

/** Gpu Kernel to quickly gather the v-space halo of local content blocks
    Halo of 1 in each direction adds up to 26 neighbours.
    For NVIDIA/CUDA, we dan do 26 neighbours and 32 threads per warp in a single block.
    For AMD/HIP, we dan do 13 neighbours and 64 threads per warp in a single block, meaning two loops per cell.
    In either case, we launch blocks equal to velocity_block_with_content_list_size
*/
//__launch_bounds__(GPUTHREADS,4)
__global__ void update_velocity_halo_kernel (
   const vmesh::VelocityMesh* __restrict__ vmesh,
   const vmesh::LocalID velocity_block_with_content_list_size, // actually not used
   const vmesh::GlobalID* __restrict__ velocity_block_with_content_list_data,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* dev_velocity_block_with_content_map,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* dev_velocity_block_with_no_content_map
   ) {
   //const int gpuBlocks = gridDim.x; // Equal to VB with content list size (or at least 1)
   const int blocki = blockIdx.x;
   const int ti = threadIdx.x;
   //const int blockSize = blockDim.x; // should be 26*32 or 13*64
   const int offsetIndex1 = ti / GPUTHREADS; // [0,26) (NVIDIA) or [0,13) (AMD)
   const int w_tid = ti % GPUTHREADS; // [0,WARPSIZE)
   // Assumes addWidthV = 1
   #ifdef __CUDACC__
   const int max_i=1;
   #endif
   #ifdef __HIP_PLATFORM_HCC___
   const int max_i=2;
   #endif
   for (int i=0; i<max_i; i++) {
      int offsetIndex = offsetIndex1 + 13*i;
      // nudge latter half in order to exclude self
      if (offsetIndex > 12) {
         offsetIndex++;
      }
      const int offset_vx = (offsetIndex % 3) - 1;
      const int offset_vy = ((offsetIndex / 3) % 3) - 1;
      const int offset_vz = (offsetIndex / 9) - 1;
      // Offsets verified in python
      const vmesh::GlobalID GID = velocity_block_with_content_list_data[blocki];
      vmesh::LocalID ind0,ind1,ind2;
      vmesh->getIndices(GID,ind0,ind1,ind2);
      const int nind0 = ind0 + offset_vx;
      const int nind1 = ind1 + offset_vy;
      const int nind2 = ind2 + offset_vz;
      const vmesh::GlobalID nGID
         = vmesh->getGlobalID(nind0,nind1,nind2);
      #ifdef USE_WARPACCESSORS
      // Does block already exist in mesh?
      const vmesh::LocalID LID = vmesh->warpGetLocalID(nGID, w_tid);
      // Try adding this nGID to velocity_block_with_content_map. If it exists, do not overwrite.
      const bool newlyadded = dev_velocity_block_with_content_map->warpInsert_V<true>(nGID,LID, w_tid);
      if (newlyadded) {
         // Block did not previously exist in velocity_block_with_content_map
         if ( LID != vmesh->invalidLocalID()) {
            // Block exists in mesh, ensure it won't get deleted:
            dev_velocity_block_with_no_content_map->warpErase(nGID, w_tid);
         }
         // else:
         // Block does not yet exist in mesh at all. Needs adding!
         // Identified as invalidLID entries in velocity_block_with_content_map.
      }
      #else
      if (w_tid==0) {
         // Does block already exist in mesh?
         const vmesh::LocalID LID = vmesh->getLocalID(nGID);
         // Add this nGID to velocity_block_with_content_map.
         const bool newEntry = dev_velocity_block_with_content_map->set_element<true>(nGID,LID);
         if (newEntry) {
            // Block did not previously exist in velocity_block_with_content_map
            if ( LID != vmesh->invalidLocalID()) {
               // Block exists in mesh, ensure it won't get deleted:
               dev_velocity_block_with_no_content_map->device_erase(nGID);
            }
         }
      }
      #endif
      __syncthreads();
   }
}

/** Gpu Kernel to quickly gather the spatial halo of neighbour content blocks
*/
//__launch_bounds__(GPUTHREADS,4)
__global__ void update_neighbour_halo_kernel (
   const vmesh::VelocityMesh* __restrict__ vmesh,
   const uint neighbour_count,
   const vmesh::GlobalID* __restrict__ const *dev_neigh_vbwcls,
   const vmesh::LocalID* __restrict__ dev_neigh_Nvbwcls,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* dev_velocity_block_with_content_map,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* dev_velocity_block_with_no_content_map
   ) {
   //const int blockSize = blockDim.x; // should be 32*32 or 16*64
   //const int gpuBlocks = gridDim.x; // Equal to count of neighbour content blocks divided by (warps/block)
   const int ti = threadIdx.x; // [0,blockSize)
   const int w_tid = ti % GPUTHREADS; // [0,WARPSIZE)
   int myindex = blockIdx.x * WARPSPERBLOCK + ti / GPUTHREADS; // Starts off as total index
   // Find which neighbour we should access
   uint neigh_i = 0;
   for (uint i=0; i<neighbour_count; i++) {
      if (myindex < (int)dev_neigh_Nvbwcls[neigh_i]) {
         break;
      }
      myindex -= dev_neigh_Nvbwcls[neigh_i];
      neigh_i++;
   } // Now we know which neighbour buffer to access
   // Access neighbour GID from buffer and act on it
   if (neigh_i < neighbour_count) {
      const vmesh::GlobalID nGID = (dev_neigh_vbwcls[neigh_i])[myindex];
      #ifdef USE_WARPACCESSORS
      // Does block already exist in mesh?
      const vmesh::LocalID LID = vmesh->warpGetLocalID(nGID, w_tid);
      // Try adding this nGID to velocity_block_with_content_map. If it exists, do not overwrite.
      const bool newlyadded = dev_velocity_block_with_content_map->warpInsert_V<true>(nGID,LID, w_tid);
      if (newlyadded) {
         // Block did not previously exist in velocity_block_with_content_map
         if ( LID != vmesh->invalidLocalID()) {
            // Block exists in mesh, ensure it won't get deleted:
            dev_velocity_block_with_no_content_map->warpErase(nGID, w_tid);
         }
         // else:
         // Block does not yet exist in mesh at all. Needs adding!
         // Identified as invalidLID entries in velocity_block_with_content_map.
      }
      #else
      if (w_tid==0) {
         // Does block already exist in mesh?
         const vmesh::LocalID LID = vmesh->getLocalID(nGID);
         // Add this nGID to velocity_block_with_content_map.
         const bool newEntry = dev_velocity_block_with_content_map->set_element<true>(nGID,LID);
         if (newEntry) {
            // Block did not previously exist in velocity_block_with_content_map
            if ( LID != vmesh->invalidLocalID()) {
               // Block exists in mesh, ensure it won't get deleted:
               dev_velocity_block_with_no_content_map->device_erase(nGID);
            }
         }
      }
      #endif
      __syncthreads();
   }
}

/** GPU kernel for resizing vmesh GlobalToLocalMap and blockcontainer
    based on vmesh localToGlobalMap contents (and size).
    Also populates block parameters.
    Assumes capacity is already large enough.
 */
__global__ void update_vmesh_and_blockparameters_kernel (
   vmesh::VelocityMesh *dev_vmesh,
   vmesh::VelocityBlockContainer *dev_blockContainer,
   const vmesh::LocalID nLIDs // newSize
   ) {
   //const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int blockSize = blockDim.x*blockDim.y*blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   // if (dev_blockContainer->size() != nLIDs) {
   //    dev_blockContainer->setNewSize(nLIDs);
   // }

   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *map = dev_vmesh->gpu_expose_map();
   split::SplitVector<vmesh::GlobalID> *list = dev_vmesh->getGrid();

   #ifdef USE_WARPACCESSORS
   const uint w_tid = ti % GPUTHREADS;
   const vmesh::LocalID LID = ti/GPUTHREADS + blocki*blockSize;
   #else
   const uint w_tid = 0;
   const vmesh::LocalID LID = ti + blocki*blockSize;
   #endif

   if (LID < nLIDs) {
      // Set velocity block parameters:
      const vmesh::GlobalID GID = (*list)[LID];

      #ifdef USE_WARPACCESSORS
      // Add this GID to the velocity_block_with_content_map.
      map->warpInsert_V<true>(GID,LID, w_tid);
      #else
      map->set_element<true>(GID,LID);
      #endif
      // Write in block parameters
      if (w_tid==0) {
         Real* blockParameters = dev_blockContainer->getParameters(LID);
         dev_vmesh->getBlockInfo(GID, blockParameters + BlockParams::VXCRD);
      }
   }
}

/** Mini-kernel for checking list sizes and attempting to adjust vmesh and VBC size on-device */
__global__ void resize_vbc_kernel_pre(
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   const split::SplitVector<vmesh::GlobalID>* __restrict__ list_with_replace_new,
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ list_delete,
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ list_to_replace,
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ list_with_replace_old,
   vmesh::LocalID* returnLID, // return values: nbefore, nafter, nblockstochange, resize success
   Realf* gpu_rhoLossAdjust // mass loss, set to zero
   ) {
   const vmesh::LocalID nBlocksBeforeAdjust = vmesh->size();
   const vmesh::LocalID nToAdd = list_with_replace_new->size();
   const vmesh::LocalID nToRemove = list_delete->size() + list_to_replace->size();
   const vmesh::LocalID nBlocksAfterAdjust = nBlocksBeforeAdjust + nToAdd - nToRemove;
   const vmesh::LocalID nBlocksToChange = nToAdd > nToRemove ? nToAdd : nToRemove;

   gpu_rhoLossAdjust[0] = 0.0;
   returnLID[0] = nBlocksBeforeAdjust;
   returnLID[1] = nBlocksAfterAdjust;
   returnLID[2] = nBlocksToChange;
   // Should we grow the size?
   if (nBlocksAfterAdjust > nBlocksBeforeAdjust) {
      if ((nBlocksAfterAdjust <= vmesh->capacity()) && (nBlocksAfterAdjust <= blockContainer->capacity())) {
         returnLID[3] = 1; // Resize on-device will work.
         vmesh->device_setNewSize(nBlocksAfterAdjust);
         blockContainer->setNewSize(nBlocksAfterAdjust);
      } else {
         returnLID[3] = 0; // Need to recapacitate and resize from host
      }
   } else {
      // No error as no resize.
      returnLID[3] = 1;
   }
}

/** Mini-kernel for adjusting vmesh and VBC size on-device aftewards (shrink only) */
__global__ void resize_vbc_kernel_post(
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   const vmesh::LocalID nBlocksAfterAdjust
   ) {
   vmesh->device_setNewSize(nBlocksAfterAdjust);
   blockContainer->setNewSize(nBlocksAfterAdjust);
}

/** GPU kernel for updating blocks based on generated lists */
__global__ void __launch_bounds__(WID3,4) update_velocity_blocks_kernel(
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   const split::SplitVector<vmesh::GlobalID>* __restrict__ list_with_replace_new,
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ list_delete,
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ list_to_replace,
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ list_with_replace_old,
   const vmesh::LocalID nBlocksBeforeAdjust,
   const vmesh::LocalID nBlocksToChange,
   const vmesh::LocalID nBlocksAfterAdjust,
   Realf* gpu_rhoLossAdjust
   ) {
   //const int gpuBlocks = gridDim.x;
   //const int blocki = blockIdx.x;
   //const int blockSize = blockDim.x; // WID3
   const uint ti = threadIdx.x; // [0,blockSize)

   // Each GPU block / workunit could manage several Vlasiator velocity blocks at once.
   // However, thread syncs inside the kernel prevent this.
   //const uint vlasiBlocksPerWorkUnit = WARPSPERBLOCK * GPUTHREADS / WID3;
   //const uint workUnitIndex = ti / WID3; // [0,vlasiBlocksPerWorkUnit)
   //const uint index = blocki * vlasiBlocksPerWorkUnit + workUnitIndex; // [0,nBlocksToChange)
   //const uint vlasiBlocksPerWorkUnit = 1;
   //const uint workUnitIndex = 1;

   // This index into vectors can be adjusted along the way
   uint index = (uint)blockIdx.x;

   const int b_tid = ti % WID3; // [0,WID3)

   const vmesh::LocalID n_with_replace_new = list_with_replace_new->size();
   const vmesh::LocalID n_delete = list_delete->size();
   const vmesh::LocalID n_to_replace = list_to_replace->size();
   const vmesh::LocalID n_with_replace_old = list_with_replace_old->size();
   // For tracking mass-loss
   //__shared__ Realf massloss[blockSize];
   __shared__ Realf massloss[WID3];

   // Each block / workunit Processes one block from the lists.

   /*********
       Check if should delete item from end of vmesh.
       For this, we get both GID and LID from the vector.
   **/
   if (index < n_delete) {
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID rmGID = (list_delete->at(index)).first;
      const vmesh::GlobalID rmLID = (list_delete->at(index)).second;
      #else
      const vmesh::GlobalID rmGID = ((*list_delete)[index]).first;
      const vmesh::GlobalID rmLID = ((*list_delete)[index]).second;
      #endif

      #ifdef DEBUG_SPATIAL_CELL
      if (rmGID == vmesh->invalidGlobalID()) {
         if (rmLID != vmesh->invalidLocalID()) {
            // Valid LID but invalid GID: only remove from vmesh localToGlobal?
            if (b_tid==0) {
               printf("Removing blocks: Valid LID %u but invalid GID!\n",rmLID);
            }
         } else {
            if (b_tid==0) {
               printf("Removing blocks: Invalid LID and GID!\n");
            }
         }
         return;
      }
      if (rmLID == vmesh->invalidLocalID()) {
         if (rmGID != vmesh->invalidGlobalID()) {
            // Valid GID but invalid LID: only remove from vmesh globalToLocal?
            if (b_tid==0) {
               printf("Removing blocks: Valid GID %ul but invalid LID!\n",rmGID);
            }
         }
         return;
      }
      if ((unsigned long)rmLID >= (unsigned long)nBlocksBeforeAdjust) {
         if (b_tid==0) {
            printf("Trying to outright remove block which has LID %ul >= nBlocksBeforeAdjust %ul!\n",rmLID,nBlocksBeforeAdjust);
         }
         return;
      }
      if ((unsigned long)rmLID < (unsigned long)nBlocksAfterAdjust) {
         if (b_tid==0) {
            printf("Trying to outright remove block which has LID %u smaller than nBlocksAfterAdjust %u!\n",rmLID,nBlocksAfterAdjust);
         }
         return;
      }
      #endif

      // Track mass loss:
      Realf* rm_avgs = blockContainer->getData(rmLID);
      Real* rm_block_parameters = blockContainer->getParameters(rmLID);
      const Real rm_DV3 = rm_block_parameters[BlockParams::DVX]
         * rm_block_parameters[BlockParams::DVY]
         * rm_block_parameters[BlockParams::DVZ];
      // thread-sum for rho
      massloss[ti] = rm_avgs[b_tid]*rm_DV3;
      __syncthreads();
      // Implemented just a simple non-optimized thread sum
      for (int s=WID3/2; s>0; s>>=1) {
         if (b_tid < s) {
            massloss[ti] += massloss[ti + s];
         }
         __syncthreads();
      }
      // Bookkeeping only by one thread
      if (b_tid==0) {
         Realf old = atomicAdd(gpu_rhoLossAdjust, massloss[ti]);
      }
      __syncthreads();

      // Delete from vmesh
      #ifdef USE_WARPACCESSORS
      vmesh->warpDeleteBlock(rmGID,rmLID,b_tid);
      #else
      if (b_tid==0) {
         vmesh->deleteBlock(rmGID,rmLID);
      }
      #endif
      // GPUTODO debug checks
      return;
   }
   index -= n_delete;

   /*********
       Check if should replace existing block with either
       existing block from end of vmesh or new block
   **/
   if (index < n_to_replace) {
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID rmGID = (list_to_replace->at(index)).first;
      const vmesh::GlobalID rmLID = (list_to_replace->at(index)).second;
      #else
      const vmesh::GlobalID rmGID = ((*list_to_replace)[index]).first;
      const vmesh::GlobalID rmLID = ((*list_to_replace)[index]).second;
      #endif
      //const vmesh::LocalID rmLID = vmesh->warpGetLocalID(rmGID,b_tid);

      #ifdef DEBUG_SPATIAL_CELL
      if (rmGID == vmesh->invalidGlobalID()) {
         if (rmLID != vmesh->invalidLocalID()) {
            // Valid LID but invalid GID: only remove from vmesh localToGlobal?
            if (b_tid==0) {
               printf("Replacing blocks: Valid LID %u but invalid GID!\n",rmLID);
            }
         } else {
            if (b_tid==0) {
               printf("Replacing blocks: Invalid LID and GID!\n");
            }
         }
         return;
      }
      if (rmLID == vmesh->invalidLocalID()) {
         if (rmGID != vmesh->invalidGlobalID()) {
            // Valid GID but invalid LID: only remove from vmesh globalToLocal?
            if (b_tid==0) {
               printf("Replacing blocks: Valid GID %ul but invalid LID!\n",rmGID);
            }
         }
         return;
      }
      if (rmLID >= nBlocksBeforeAdjust) {
         if (b_tid==0) {
            printf("Trying to replace block which has LID %ul >= nBlocksBeforeAdjust %ul!\n",rmLID,nBlocksBeforeAdjust);
         }
         return;
      }
      #endif

      // Track mass loss:
      Realf* rm_avgs = blockContainer->getData(rmLID);
      Real* rm_block_parameters = blockContainer->getParameters(rmLID);
      const Real rm_DV3 = rm_block_parameters[BlockParams::DVX]
         * rm_block_parameters[BlockParams::DVY]
         * rm_block_parameters[BlockParams::DVZ];
      // thread-sum for rho
      massloss[ti] = rm_avgs[b_tid]*rm_DV3;
      __syncthreads();
      // Implemented just a simple non-optimized thread sum
      for (int s=WID3/2; s>0; s>>=1) {
         if (b_tid < s) {
            massloss[ti] += massloss[ti + s];
         }
         __syncthreads();
      }
      // Bookkeeping only by one thread
      if (b_tid==0) {
         Realf old = atomicAdd(gpu_rhoLossAdjust, massloss[ti]);
      }
      __syncthreads();

      // Figure out what to use as replacement
      vmesh::GlobalID replaceGID;
      vmesh::LocalID replaceLID;

      // First option: replace with existing block from end of vmesh
      if (index < n_with_replace_old) {
         #ifdef DEBUG_SPATIAL_CELL
         replaceGID = (list_with_replace_old->at(index)).first;
         replaceLID = (list_with_replace_old->at(index)).second;
         #else
         replaceGID = ((*list_with_replace_old)[index]).first;
         replaceLID = ((*list_with_replace_old)[index]).second;
         #endif

         Realf* repl_avgs = blockContainer->getData(replaceLID);
         Real*  repl_block_parameters = blockContainer->getParameters(replaceLID);
         rm_avgs[b_tid] = repl_avgs[b_tid];
         if (b_tid < BlockParams::N_VELOCITY_BLOCK_PARAMS) {
            rm_block_parameters[b_tid] = repl_block_parameters[b_tid];
         }
         __syncthreads();

      } else {
         // Second option: add new block instead
         #ifdef DEBUG_SPATIAL_CELL
         replaceGID = list_with_replace_new->at(index - n_with_replace_old);
         #else
         replaceGID = (*list_with_replace_new)[index - n_with_replace_old];
         #endif
         replaceLID = vmesh->invalidLocalID();

         rm_avgs[b_tid] = 0;
         if (b_tid==0) {
            // Write in block parameters
            vmesh->getBlockInfo(replaceGID, rm_block_parameters+BlockParams::VXCRD);
         }
         __syncthreads();
      }
      // Remove hashmap entry for removed block, add instead created block
      #ifdef USE_WARPACCESSORS
      vmesh->warpReplaceBlock(rmGID,rmLID,replaceGID,b_tid);
      #else
      if (b_tid==0) {
         vmesh->replaceBlock(rmGID,rmLID,replaceGID);
      }
      #endif

      #ifdef DEBUG_SPATIAL_CELL
      __syncthreads();
      if (vmesh->getGlobalID(rmLID) != replaceGID) {
         if (b_tid==0) {
            printf("Error! Replacing did not result in wanted GID at old LID in update_velocity_blocks_kernel! \n");
         }
         __syncthreads();
      }
      if (vmesh->getLocalID(replaceGID) != rmLID) {
         if (b_tid==0) {
            printf("Error! Replacing did not result in old LID at replaced GID in update_velocity_blocks_kernel! \n");
         }
         __syncthreads();
      }
      #endif

      return;
   }
   index -= n_to_replace;

   /*********
       Finally check if we should add new block after end of current vmesh
       We have reserved/used some entries from the beginning of the list_with_replace_new
       for the previous section, so now we access that with a different index.
   **/
   const uint add_index = index + (n_to_replace - n_with_replace_old);
   if (add_index < n_with_replace_new) {
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID addGID = list_with_replace_new->at(add_index);
      if (b_tid==0) {
         if (vmesh->getLocalID(addGID) != vmesh->invalidLocalID()) {
            printf("Trying to add new GID %u to mesh which already contains it! index=%u addindex=%u\n",addGID,index,add_index);
         }
      }
      #else
      const vmesh::GlobalID addGID = (*list_with_replace_new)[add_index];
      #endif

      // We need to add the data of addGID to a new LID. Here we still use the regular index.
      const vmesh::LocalID addLID = nBlocksBeforeAdjust + index;
      Realf* add_avgs = blockContainer->getData(addLID);
      #ifdef DEBUG_SPATIAL_CELL
      if (addGID == vmesh->invalidGlobalID()) {
         printf("Error! invalid addGID!\n");
         return;
      }
      if (addLID == vmesh->invalidLocalID()) {
         printf("Error! invalid addLID!\n");
         return;
      }
      #endif
      Real* add_block_parameters = blockContainer->getParameters(addLID);
      // Zero out blockdata
      add_avgs[b_tid] = 0;
      if (b_tid==0) {
         // Write in block parameters
         vmesh->getBlockInfo(addGID, add_block_parameters+BlockParams::VXCRD);
      }
      __syncthreads();

      // Insert new hashmap entry into vmesh
      #ifdef USE_WARPACCESSORS
      vmesh->warpPlaceBlock(addGID,addLID,b_tid);
      #else
      if (b_tid==0) {
         vmesh->placeBlock(addGID,addLID);
      }
      #endif

      #ifdef DEBUG_SPATIAL_CELL
      __syncthreads();
      if (vmesh->getGlobalID(addLID) == vmesh->invalidGlobalID()) {
         printf("Error! invalid GID after add from addLID!\n");
      }
      if (vmesh->getLocalID(addGID) == vmesh->invalidLocalID()) {
         printf("Error! invalid LID after add from addGID!\n");
      }
      #endif
      return;
   }

   // Fall-through error!
   if (b_tid==0) {
      printf("Error! Fall through in update_velocity_blocks_kernel! index %u nBlocksBeforeAdjust %u nBlocksAfterAdjust %u \n",
             index,nBlocksBeforeAdjust,nBlocksAfterAdjust);
   }
   __syncthreads();
}

#endif
