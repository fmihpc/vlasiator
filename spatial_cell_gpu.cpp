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

#include <unordered_set>

#include "spatial_cell_gpu.hpp"
#include "arch/gpu_base.hpp"
#include "object_wrapper.h"
#include "velocity_mesh_parameters.h"

using namespace std;

// GPUTODO: Make error-checking functions to be called inside kernels instead of duplicating so much code.
// Same for gathering mass loss.

/** GPU kernel for identifying which blocks have relevant content */
__global__ void __launch_bounds__(WID3,4) update_velocity_block_content_lists_kernel (
   vmesh::VelocityMesh *vmesh,
   //const uint nBlocks,
   vmesh::VelocityBlockContainer *blockContainer,
   //vmesh::GlobalID* vbwcl_gather,
   //vmesh::GlobalID* vbwncl_gather,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwcl_map,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwncl_map,
   Real velocity_block_min_value
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
   //for (uint blockLID=blocki; blockLID<nBlocks; blockLID += gpuBlocks) {
   if (blockLID < nBlocks) {
      const vmesh::GlobalID blockGID = vmesh->getGlobalID(blockLID);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == vmesh->invalidGlobalID()) {
         if (b_tid==0) printf("Invalid GID encountered in update_velocity_block_content_lists_kernel!\n");
         return;
      }
      if (blockLID == vmesh->invalidLocalID()) {
         if (b_tid==0) printf("Invalid LID encountered in update_velocity_block_content_lists_kernel!\n");
         return;
      }
      #endif
      // Check each velocity cell if it is above the threshold
      const Realf* avgs = blockContainer->getData(blockLID);
      has_content[ti] = avgs[b_tid] >= velocity_block_min_value ? 1 : 0;
      __syncthreads(); // THIS SYNC IS CRUCIAL!
      // Implemented just a simple non-optimized thread OR
      // GPUTODO reductions via warp voting
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (b_tid < s) {
            has_content[ti] = has_content[ti] || has_content[ti + s];
         }
         __syncthreads();
      }
      // Insert into map only from threads 0...WARPSIZE
      if (b_tid < GPUTHREADS) {
         if (has_content[0]) {
            vbwcl_map->warpInsert(blockGID,blockLID,b_tid);
         } else {
            vbwncl_map->warpInsert(blockGID,blockLID,b_tid);
         }
      }
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
   vmesh::VelocityMesh *vmesh,
   vmesh::LocalID velocity_block_with_content_list_size, // actually not used
   vmesh::GlobalID* velocity_block_with_content_list_data,
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
      // Does block already exist in mesh?
      const vmesh::LocalID LID = vmesh->warpGetLocalID(nGID, w_tid);
      // Try adding this nGID to velocity_block_with_content_map. If it exists, do not overwrite.
      const bool newlyadded = dev_velocity_block_with_content_map->warpInsert_V<true>(nGID,LID, w_tid);
      if (newlyadded) {
         // Block did not previously exist in velocity_block_with_content_map
         if ( LID != vmesh->invalidLocalID()) {
            // Block exists in mesh, ensure it won't get deleted:
            // try deleting from no_content map
            dev_velocity_block_with_no_content_map->warpErase(nGID, w_tid);
         }
         // else:
         // Block does not yet exist in mesh at all. Needs adding!
         // Identified as invalidLID entries in velocity_block_with_content_map.
      }
      __syncthreads();
   }
}

/** Gpu Kernel to quickly gather the spatial halo of neighbour content blocks
*/
//__launch_bounds__(GPUTHREADS,4)
__global__ void update_neighbour_halo_kernel (
   vmesh::VelocityMesh *vmesh,
   uint neighbour_count,
   vmesh::GlobalID **dev_neigh_vbwcls,
   vmesh::LocalID *dev_neigh_Nvbwcls,
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
      // Does block already exist in mesh?
      const vmesh::LocalID LID = vmesh->warpGetLocalID(nGID, w_tid);
      // Try adding this nGID to velocity_block_with_content_map. If it exists, do not overwrite.
      const bool newlyadded = dev_velocity_block_with_content_map->warpInsert_V<true>(nGID,LID, w_tid);
      if (newlyadded) {
         // Block did not previously exist in velocity_block_with_content_map
         if ( LID != vmesh->invalidLocalID()) {
            // Block exists in mesh, ensure it won't get deleted:
            // try deleting from no_content map
            dev_velocity_block_with_no_content_map->warpErase(nGID, w_tid);
         }
         // else:
         // Block does not yet exist in mesh at all. Needs adding!
         // Identified as invalidLID entries in velocity_block_with_content_map.
      }
      __syncthreads();
   }
}

/** GPU kernel for quickly filling block parameters.
 */
__global__ void update_blockparameters_kernel (
   vmesh::VelocityMesh *dev_vmesh,
   vmesh::VelocityBlockContainer *dev_blockContainer,
   vmesh::LocalID nLIDs
   ) {
   //const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const vmesh::LocalID LID = ti + blocki*warpSize;
   if (LID < nLIDs) {
      // Set velocity block parameters:
      const vmesh::GlobalID GID = dev_vmesh->getGlobalID(LID);
      Real* blockParameters = dev_blockContainer->getParameters(LID);
      // Write in block parameters
      dev_vmesh->getBlockInfo(GID, blockParameters + BlockParams::VXCRD);
   }
}

/** Mini-kernel for checking list sizes and attempting to adjust vmesh and VBC size on-device */
__global__ void resize_vbc_kernel_pre(
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   split::SplitVector<vmesh::GlobalID>* list_with_replace_new,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_delete,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_to_replace,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_with_replace_old,
   vmesh::LocalID* returnLID, // return values: nbefore, nafter, nblockstochange, resize success
   Realf* gpu_rhoLossAdjust // mass loss, set to zero
   ) {
   const vmesh::LocalID n_to_replace = list_to_replace->size(); // replace these blocks
   const vmesh::LocalID n_with_replace_new = list_with_replace_new->size(); // use to replace, or add at end
   const vmesh::LocalID n_with_replace_old = list_with_replace_old->size(); // use to replace
   const vmesh::LocalID n_to_delete = list_delete->size(); // delete from end

   const vmesh::LocalID nBlocksBeforeAdjust = vmesh->size();
   //const vmesh::LocalID nBlocksAfterAdjust = nBlocksBeforeAdjust + n_with_replace_new - n_to_delete;
   //const vmesh::LocalID nBlocksToChange = n_with_replace_new + n_with_replace_old + n_to_delete;

   // DEBUG
   const vmesh::LocalID nToAdd = list_with_replace_new->size();
   const vmesh::LocalID nToRemove = list_delete->size() + list_to_replace->size();
   const vmesh::LocalID nBlocksAfterAdjust = nBlocksBeforeAdjust + nToAdd - nToRemove;
   //const vmesh::LocalID nBlocksToChange = n_with_replace_new + n_with_replace_old + n_to_delete;
   //OR?
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
   vmesh::LocalID nBlocksAfterAdjust
   ) {
   vmesh->device_setNewSize(nBlocksAfterAdjust);
   blockContainer->setNewSize(nBlocksAfterAdjust);
}

/** GPU kernel for updating blocks based on generated lists */
__global__ void __launch_bounds__(WID3,4) update_velocity_blocks_kernel(
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   split::SplitVector<vmesh::GlobalID>* list_with_replace_new,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_delete,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_to_replace,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_with_replace_old,
   vmesh::LocalID nBlocksBeforeAdjust,
   vmesh::LocalID nBlocksToChange,
   vmesh::LocalID nBlocksAfterAdjust,
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
      vmesh->warpDeleteBlock(rmGID,rmLID,b_tid);
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
      vmesh->warpReplaceBlock(rmGID,rmLID,replaceGID,b_tid);
      #ifdef DEBUG_SPATIAL_CELL
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
      vmesh->warpPlaceBlock(addGID,addLID,b_tid);
      #ifdef DEBUG_SPATIAL_CELL
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

namespace spatial_cell {
   int SpatialCell::activePopID = 0;
   uint64_t SpatialCell::mpi_transfer_type = 0;
   bool SpatialCell::mpiTransferAtSysBoundaries = false;
   bool SpatialCell::mpiTransferInAMRTranslation = false;
   int SpatialCell::mpiTransferXYZTranslation = 0;

   SpatialCell::SpatialCell() {
      // Block list and cache always have room for all blocks
      this->sysBoundaryLayer=0; // Default value, layer not yet initialized
      for (unsigned int i=0; i<WID3; ++i) {
         null_block_data[i] = 0.0;
      }

      // reset spatial cell parameters
      for (unsigned int i = 0; i < CellParams::N_SPATIAL_CELL_PARAMS; i++) {
         this->parameters[i]=0.0;
      }

      // reset BVOL derivatives
      for (unsigned int i = 0; i < bvolderivatives::N_BVOL_DERIVATIVES; i++) {
         this->derivativesBVOL[i]=0;
      }

      for (unsigned int i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
         this->neighbor_number_of_blocks[i] = 0;
         this->neighbor_block_data[i] = NULL;
      }

      //is transferred by default
      this->mpiTransferEnabled=true;

      // Set correct number of populations
      populations.resize(getObjectWrapper().particleSpecies.size());

      // Set velocity meshes
      for (uint popID=0; popID<populations.size(); ++popID) {
         const species::Species& spec = getObjectWrapper().particleSpecies[popID];
         populations[popID].vmesh->initialize(spec.velocityMesh);
         populations[popID].Upload();
         populations[popID].velocityBlockMinValue = spec.sparseMinValue;
         populations[popID].N_blocks = 0;
      }

      // SplitVectors and hashmaps via pointers for unified memory

      // create in host instead of unified memory, upload device copy
      void *buf0 = malloc(sizeof(split::SplitVector<vmesh::GlobalID>));
      velocity_block_with_content_list = ::new (buf0) split::SplitVector<vmesh::GlobalID>(1);
      //velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(1);
      dev_velocity_block_with_content_list = velocity_block_with_content_list->upload();
      velocity_block_with_content_list->clear();
      velocity_block_with_content_list_size=0;
      velocity_block_with_content_list_capacity=1;

      // create in host instead of unified memory, upload device copy
      void *buf1 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
      void *buf2 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
      velocity_block_with_content_map = ::new (buf1) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      velocity_block_with_no_content_map = ::new (buf2)Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);

      // velocity_block_with_content_map = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      // velocity_block_with_no_content_map = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      dev_velocity_block_with_content_map = velocity_block_with_content_map->upload();
      dev_velocity_block_with_no_content_map = velocity_block_with_no_content_map->upload();
      vbwcl_sizePower = 7;
      vbwncl_sizePower = 7;
   }

   SpatialCell::~SpatialCell() {
      gpu_destructor();
   }

   void SpatialCell::gpu_destructor() {
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         populations[popID].gpu_destructor();
      }
      if (velocity_block_with_content_list) {
         ::delete velocity_block_with_content_list;
         velocity_block_with_content_list = 0;
      }
      if (velocity_block_with_content_map) {
         ::delete velocity_block_with_content_map;
         velocity_block_with_content_map = 0;
      }
      if (velocity_block_with_no_content_map) {
         ::delete velocity_block_with_no_content_map;
         velocity_block_with_no_content_map = 0;
      }
      velocity_block_with_content_list_size=0;
      velocity_block_with_content_list_capacity=0;
      vbwcl_sizePower = 0;
      vbwncl_sizePower = 0;
   }

   SpatialCell::SpatialCell(const SpatialCell& other) {
      const uint reserveSize = other.velocity_block_with_content_list_capacity;

      // create in host instead of unified memory, upload device copy
      void *buf0 = malloc(sizeof(split::SplitVector<vmesh::GlobalID>));
      velocity_block_with_content_list = ::new (buf0) split::SplitVector<vmesh::GlobalID>(reserveSize);
      //velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(reserveSize);
      dev_velocity_block_with_content_list = velocity_block_with_content_list->upload();
      velocity_block_with_content_list->clear();
      velocity_block_with_content_list_size = 0;
      velocity_block_with_content_list_capacity = reserveSize;

      // create in host instead of unified memory, upload device copy
      void *buf1 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
      void *buf2 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
      velocity_block_with_content_map = ::new (buf1) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(other.vbwcl_sizePower);
      velocity_block_with_no_content_map = ::new (buf2)Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(other.vbwncl_sizePower);
      dev_velocity_block_with_content_map = velocity_block_with_content_map->upload();
      dev_velocity_block_with_no_content_map = velocity_block_with_no_content_map->upload();
      vbwcl_sizePower = other.vbwcl_sizePower;
      vbwncl_sizePower = other.vbwncl_sizePower;

      // Member variables
      ioLocalCellId = other.ioLocalCellId;
      sysBoundaryFlag = other.sysBoundaryFlag;
      sysBoundaryLayer = other.sysBoundaryLayer;
      sysBoundaryLayerNew = other.sysBoundaryLayerNew;
      initialized = other.initialized;
      mpiTransferEnabled = other.mpiTransferEnabled;
      for (unsigned int i=0; i<bvolderivatives::N_BVOL_DERIVATIVES; ++i) {
         derivativesBVOL[i] = other.derivativesBVOL[i];
      }
      for (unsigned int i=0; i<CellParams::N_SPATIAL_CELL_PARAMS; ++i) {
         parameters[i] = other.parameters[i];
      }
      for (unsigned int i=0; i<WID3; ++i) {
         null_block_data[i] = 0.0;
      }
      for (unsigned int i=0; i<MAX_NEIGHBORS_PER_DIM; ++i) {
         neighbor_block_data[i] = 0;
         neighbor_number_of_blocks[i] = 0;
      }
      face_neighbor_ranks.clear();
      // for (unsigned int i=0; i<MAX_NEIGHBORS_PER_DIM; ++i) {
      //    neighbor_block_data[i] = other.neighbor_block_data[i];
      //    neighbor_number_of_blocks[i] = other.neighbor_number_of_blocks[i];
      // }
      // if (other.face_neighbor_ranks.size()>0) {
      //    face_neighbor_ranks = std::map<int,std::set<int>>(other.face_neighbor_ranks);
      // }
      if (other.populations.size()>0) {
         populations = std::vector<spatial_cell::Population>(other.populations);
      }
   }

   const SpatialCell& SpatialCell::operator=(const SpatialCell& other) {
      const uint reserveSize = other.velocity_block_with_content_list_capacity;
      velocity_block_with_content_list->reserve(reserveSize);
      velocity_block_with_content_list->clear();
      velocity_block_with_content_list_size = 0;
      velocity_block_with_content_list_capacity=reserveSize;
      dev_velocity_block_with_content_list = velocity_block_with_content_list->upload();

      velocity_block_with_content_map->resize(other.vbwcl_sizePower);
      velocity_block_with_no_content_map->resize(other.vbwncl_sizePower);
      dev_velocity_block_with_content_map = velocity_block_with_content_map->upload();
      dev_velocity_block_with_no_content_map = velocity_block_with_no_content_map->upload();
      vbwcl_sizePower = other.vbwcl_sizePower;
      vbwncl_sizePower = other.vbwncl_sizePower;

      // Member variables
      ioLocalCellId = other.ioLocalCellId;
      sysBoundaryFlag = other.sysBoundaryFlag;
      sysBoundaryLayer = other.sysBoundaryLayer;
      sysBoundaryLayerNew = other.sysBoundaryLayerNew;
      initialized = other.initialized;
      mpiTransferEnabled = other.mpiTransferEnabled;
      for (unsigned int i=0; i<bvolderivatives::N_BVOL_DERIVATIVES; ++i) {
         derivativesBVOL[i] = other.derivativesBVOL[i];
      }
      for (unsigned int i=0; i<CellParams::N_SPATIAL_CELL_PARAMS; ++i) {
         parameters[i] = other.parameters[i];
      }
      for (unsigned int i=0; i<WID3; ++i) {
         null_block_data[i] = 0.0;
      }
      for (unsigned int i=0; i<MAX_NEIGHBORS_PER_DIM; ++i) {
         neighbor_block_data[i] = 0;
         neighbor_number_of_blocks[i] = 0;
      }
      face_neighbor_ranks.clear();
      // for (unsigned int i=0; i<MAX_NEIGHBORS_PER_DIM; ++i) {
      //    neighbor_block_data[i] = other.neighbor_block_data[i];
      //    neighbor_number_of_blocks[i] = other.neighbor_number_of_blocks[i];
      // }
      //face_neighbor_ranks = std::map<int,std::set<int>>(other.face_neighbor_ranks);
      populations = std::vector<spatial_cell::Population>(other.populations);

      return *this;
   }

   /** Sets a guidance counter so that vmesh adjustment vectors have sufficient size
    */
   void SpatialCell::setReservation(const uint popID, const vmesh::LocalID reservationsize, bool force) {
      if (force || (reservationsize > populations[popID].reservation)) {
         populations[popID].reservation = reservationsize;
      }
   }
   vmesh::LocalID SpatialCell::getReservation(const uint popID) const {
      return populations[popID].reservation;
   }
   /** Recapacitates local temporary vectors based on guidance counter
    */
   void SpatialCell::applyReservation(const uint popID) {
      const size_t reserveSize = populations[popID].reservation * BLOCK_ALLOCATION_FACTOR;
      size_t newReserve = populations[popID].reservation * BLOCK_ALLOCATION_PADDING;
      const vmesh::LocalID HashmapReqSize = ceil(log2(reserveSize));
      gpuStream_t stream = gpu_getStream();
      // Now uses host-cached values

      // loop extraction used to require extra buffer
      //newReserve = ((newReserve /(WARPSPERBLOCK*GPUTHREADS))+2) * WARPSPERBLOCK * GPUTHREADS;

      if (velocity_block_with_content_list_capacity < newReserve) {
         //phiprof::Timer vectorTimer {"with_vector"};
         velocity_block_with_content_list->reserve(newReserve,true);
         velocity_block_with_content_list_capacity = newReserve;
         dev_velocity_block_with_content_list = velocity_block_with_content_list->upload(stream);
         //velocity_block_with_content_list->optimizeGPU(stream); // included in upload()
         //vectorTimer.stop();
      }
      if (vbwcl_sizePower < HashmapReqSize) {
         //phiprof::Timer map1cTimer {"with_map1 clear"};
         //velocity_block_with_content_map->clear<false>(Hashinator::targets::device,stream);
         //map1cTimer.stop();
         vbwcl_sizePower = HashmapReqSize+2;
         //velocity_block_with_content_map->resize(vbwcl_sizePower, targets::device, stream);
         //phiprof::Timer map1rTimer {"with_map1 resize"};
         //velocity_block_with_content_map->resize(vbwcl_sizePower,Hashinator::targets::device, stream);
         ::delete velocity_block_with_content_map;
         void *buf = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
         velocity_block_with_content_map = ::new (buf) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(vbwcl_sizePower);
         // velocity_block_with_content_map = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(vbwcl_sizePower);
         dev_velocity_block_with_content_map = velocity_block_with_content_map->upload(stream);
         //velocity_block_with_content_map->optimizeGPU(stream); // included in upload()
         //map1cTimer.stop();
      }
      if (vbwncl_sizePower < HashmapReqSize) {
         //phiprof::Timer map2cTimer {"with_map2 clear"};
         //velocity_block_with_no_content_map->clear<false>(Hashinator::targets::device,stream);
         //map2cTimer.stop();
         vbwncl_sizePower = HashmapReqSize+2;
         //phiprof::Timer map2rTimer {"with_map2 resize"};
         //velocity_block_with_no_content_map->resize(vbwncl_sizePower,Hashinator::targets::device, stream);
         ::delete velocity_block_with_no_content_map;
         void *buf = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
         velocity_block_with_no_content_map = ::new (buf) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(vbwncl_sizePower);
         //velocity_block_with_no_content_map = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(vbwncl_sizePower);
         dev_velocity_block_with_no_content_map = velocity_block_with_no_content_map->upload(stream);
         //velocity_block_with_no_content_map->optimizeGPU(stream); // included in upload()
         //map2cTimer.stop();
      }
   }

   /** Adds "important" and removes "unimportant" velocity blocks
    * to/from this cell.
    *
    * velocity_block_with_content_list needs to be up to date in local and remote cells.
    * velocity_block_with_no_content_list needs to be up to date in local cells.
    *
    * update_velocity_block_with_content_lists() should have
    * been called with the current distribution function values, and then the contetn list transferred.
    *
    * Removes all velocity blocks from this spatial cell which don't
    * have content and don't have spatial or velocity neighbors with
    * content.  Adds neighbors for all velocity blocks which do have
    * content (including spatial neighbors).  All cells in
    * spatial_neighbors are assumed to be neighbors of this cell.
    *
    * This function is thread-safe when called for different cells
    * per thread. We need the block_has_content vector from
    * neighbouring cells, but these are not written to here. We only
    * modify local cell.*/

   void SpatialCell::adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors,
                                            const uint popID, bool doDeleteEmptyBlocks) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      if (this->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         return;
      }

      const uint cpuThreadID = gpu_getThread();
      const gpuStream_t stream = gpu_getStream();

      // Ensure allocation
      const uint allocationSize = (populations[popID].reservation > velocity_block_with_content_list_size)
         ? populations[popID].reservation : velocity_block_with_content_list_size;
      gpu_blockadjust_allocate_perthread(cpuThreadID,allocationSize);

      phiprof::Timer adjustBlocksTimer {"Adjust velocity blocks"};

      vmesh::VelocityMesh* host_vmesh    = populations[popID].vmesh;
      vmesh::VelocityMesh* dev_vmesh    = populations[popID].dev_vmesh;
      vmesh::GlobalID* _withContentData = velocity_block_with_content_list->data();

      split::SplitVector<vmesh::GlobalID> *list_with_replace_new = gpu_list_with_replace_new[cpuThreadID];
      split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_delete = gpu_list_delete[cpuThreadID];
      split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_to_replace = gpu_list_to_replace[cpuThreadID];
      split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_with_replace_old = gpu_list_with_replace_old[cpuThreadID];

      // Evaluate velocity halo for local content blocks
      if (velocity_block_with_content_list_size>0) {
         //phiprof::Timer blockHaloTimer {"Block halo kernel"};
         const int addWidthV = getObjectWrapper().particleSpecies[popID].sparseBlockAddWidthV;
         if (addWidthV!=1) {
            std::cerr<<"Warning! "<<__FILE__<<":"<<__LINE__<<" Halo extent is not 1, unsupported size."<<std::endl;
         }
         // Halo of 1 in each direction adds up to 26 velocity neighbours.
         // For NVIDIA/CUDA, we dan do 26 neighbours and 32 threads per warp in a single block.
         // For AMD/HIP, we dan do 13 neighbours and 64 threads per warp in a single block, meaning two loops per cell.
         // In either case, we launch blocks equal to velocity_block_with_content_list_size
         update_velocity_halo_kernel<<<velocity_block_with_content_list_size, 26*32, 0, stream>>> (
            dev_vmesh,
            velocity_block_with_content_list_size,
            _withContentData,
            dev_velocity_block_with_content_map,
            dev_velocity_block_with_no_content_map
            );
         CHK_ERR( gpuPeekAtLastError() );
         //CHK_ERR( gpuStreamSynchronize(stream) );
         //blockHaloTimer.stop();
      }

      // Gather pointers and counts from neighbours
      uint neighbours_count = spatial_neighbors.size();
      uint neighbours_blocks_count = 0;
      std::vector<vmesh::GlobalID*> neigh_vbwcls;
      std::vector<vmesh::LocalID> neigh_Nvbwcls;

      if (neighbours_count > 0) {
         for (std::vector<SpatialCell*>::const_iterator neighbor=spatial_neighbors.begin();
              neighbor != spatial_neighbors.end(); ++neighbor) {
            if ((*neighbor)->velocity_block_with_content_list_size > 0) {
               neigh_Nvbwcls.push_back((*neighbor)->velocity_block_with_content_list_size);
               neigh_vbwcls.push_back((*neighbor)->velocity_block_with_content_list->data());
               neighbours_blocks_count += (*neighbor)->velocity_block_with_content_list_size;
            }
         }
         neighbours_count = neigh_Nvbwcls.size(); // Only manage neighbours with content.
      }

      if (neighbours_count > 0) {
         // Upload pointers and counters for neighbours
         vmesh::GlobalID** dev_neigh_vbwcls;
         vmesh::GlobalID* dev_neigh_Nvbwcls;
         CHK_ERR( gpuMallocAsync((void**)&dev_neigh_vbwcls, neighbours_count*sizeof(vmesh::GlobalID*), stream) );
         CHK_ERR( gpuMallocAsync((void**)&dev_neigh_Nvbwcls, neighbours_count*sizeof(vmesh::LocalID), stream) );
         CHK_ERR( gpuMemcpyAsync(dev_neigh_vbwcls, neigh_vbwcls.data(), neighbours_count*sizeof(vmesh::GlobalID*), gpuMemcpyHostToDevice, stream) );
         CHK_ERR( gpuMemcpyAsync(dev_neigh_Nvbwcls, neigh_Nvbwcls.data(), neighbours_count*sizeof(vmesh::LocalID), gpuMemcpyHostToDevice, stream) );
         //phiprof::Timer neighHaloTimer {"Neighbour halo kernel"};
         // For NVIDIA/CUDA, we dan do 26 neighbours and 32 threads per warp in a single block.
         // For AMD/HIP, we dan do 13 neighbours and 64 threads per warp in a single block, meaning two loops per cell.
         // This is managed in-kernel.
         // In either case, we launch blocks equal to velocity_block_with_content_list_size
         // We always launch at least one block in order to clear splitvectors.
         // ceil int division
         uint launchBlocks = 1 + ((neighbours_blocks_count - 1) / WARPSPERBLOCK);
         if (launchBlocks < std::pow(2,31)) {
            update_neighbour_halo_kernel<<<launchBlocks, WARPSPERBLOCK*GPUTHREADS, 0, stream>>> (
               dev_vmesh,
               neighbours_count,
               dev_neigh_vbwcls,
               dev_neigh_Nvbwcls,
               dev_velocity_block_with_content_map,
               dev_velocity_block_with_no_content_map
               );
            CHK_ERR( gpuPeekAtLastError() );
         } else {
            // Too many launch blocks, call one by one (unlikely)
            for (uint neigh_i = 0; neigh_i < neighbours_count; neigh_i++) {
               uint launchBlocks = 1 + ((neigh_Nvbwcls[neigh_i] - 1) / WARPSPERBLOCK);
               update_neighbour_halo_kernel<<<launchBlocks, WARPSPERBLOCK*GPUTHREADS, 0, stream>>> (
                  dev_vmesh,
                  1,
                  dev_neigh_vbwcls+neigh_i,
                  dev_neigh_Nvbwcls+neigh_i,
                  dev_velocity_block_with_content_map,
                  dev_velocity_block_with_no_content_map
                  );
               CHK_ERR( gpuPeekAtLastError() );
            }
         }
         //CHK_ERR( gpuStreamSynchronize(stream) );
         //neighHaloTimer.stop();
      }

      // Now extract vectors to be used in actual block adjustment
      // Previous kernels may have added dummy (to be added) entries to
      // velocity_block_with_content_map with LID=vmesh->invalidLocalID()
      // Or non-content blocks which should be retained (with correct LIDs).

      /** Rules used in extracting keys or elements from hashmaps
          Now these include passing pointers to GPU memory in order to evaluate
          nBlocksAfterAdjust without going via host. Pointers are copied by value.
       */
      const vmesh::GlobalID EMPTYBUCKET = std::numeric_limits<vmesh::GlobalID>::max();
      const vmesh::GlobalID TOMBSTONE   = EMPTYBUCKET - 1;
      const vmesh::GlobalID invalidGID  = host_vmesh->invalidGlobalID();
      const vmesh::LocalID  invalidLID  = host_vmesh->invalidLocalID();
      // Required GIDs which do not yet exist in vmesh were stored in velocity_block_with_content_map with invalidLID
      auto rule_add = [EMPTYBUCKET, TOMBSTONE, invalidGID, invalidLID]
         __device__(const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval) -> bool {
                         return kval.first != EMPTYBUCKET &&
                            kval.first != TOMBSTONE &&
                            kval.first != invalidGID &&
                            kval.second == invalidLID; };
      velocity_block_with_content_map->extractKeysByPatternLoop(*list_with_replace_new, rule_add, stream);

      if (doDeleteEmptyBlocks) {
         Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *vbwncm = dev_velocity_block_with_no_content_map;
         split::SplitVector<vmesh::GlobalID> *d_list_add = list_with_replace_new;

         auto rule_delete_move = [EMPTYBUCKET, TOMBSTONE, vbwncm, d_list_add, dev_vmesh, invalidGID, invalidLID]
            __device__(const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval) -> bool {
                                    const vmesh::LocalID nBlocksAfterAdjust1 = dev_vmesh->size()
                                       + d_list_add->size() - vbwncm->size();
                                    return kval.first != EMPTYBUCKET &&
                                       kval.first != TOMBSTONE &&
                                       kval.first != invalidGID &&
                                       kval.second >= nBlocksAfterAdjust1 &&
                                       kval.second != invalidLID; };
         auto rule_to_replace = [EMPTYBUCKET, TOMBSTONE, vbwncm, d_list_add, dev_vmesh, invalidGID, invalidLID]
            __device__(const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval) -> bool {
                                   const vmesh::LocalID nBlocksAfterAdjust2 = dev_vmesh->size()
                                      + d_list_add->size() - vbwncm->size();
                                   return kval.first != EMPTYBUCKET &&
                                      kval.first != TOMBSTONE &&
                                      kval.first != invalidGID &&
                                      kval.second < nBlocksAfterAdjust2 &&
                                      kval.second != invalidGID; };

         velocity_block_with_content_map->extractPatternLoop(*list_with_replace_old, rule_delete_move, stream);
         velocity_block_with_no_content_map->extractPatternLoop(*list_delete, rule_delete_move, stream);
         velocity_block_with_no_content_map->extractPatternLoop(*list_to_replace, rule_to_replace, stream);
      } else {
         list_with_replace_old->clear();
         list_delete->clear();
         list_to_replace->clear();
      }
      // Note:list_with_replace_new contains both new GIDs to use for replacements and new GIDs to place at end of vmesh

      // Actual adjustment calling happens in separate function as it is also called from within acceleration
      vmesh::LocalID nBlocksAfterAdjust = adjust_velocity_blocks_caller(popID);

      // Perform hashmap cleanup here (instead of at acceleration mid-steps)
      phiprof::Timer cleanupTimer {"Hashinator cleanup"};
      populations[popID].vmesh->gpu_cleanHashMap(stream);
      populations[popID].Upload();
      cleanupTimer.stop();

      #ifdef DEBUG_SPATIAL_CELL
      const size_t vmeshSize = (populations[popID].vmesh)->size();
      const size_t vbcSize = (populations[popID].blockContainer)->size();
      if (vmeshSize != vbcSize) {
         printf("ERROR: population vmesh %zu and blockcontainer %zu sizes do not match!\n",vmeshSize,vbcSize);
      }
      #endif
      #ifdef DEBUG_VLASIATOR
      // This is a bit extreme
      if (!populations[popID].vmesh->check()) {
         printf("ERROR in vmesh check: %s at %d\n",__FILE__,__LINE__);
      }
      #endif
   }

   /**
      Call GPU kernel with all necessary information for creation and deletion of blocks.
   **/
   vmesh::LocalID SpatialCell::adjust_velocity_blocks_caller(const uint popID) {
      phiprof::Timer addRemoveTimer {"GPU add and remove blocks"};
      Realf host_rhoLossAdjust = 0;
      vmesh::LocalID host_returnLID[4];

      const uint cpuThreadID = gpu_getThread();
      const gpuStream_t stream = gpu_getStream();
      // populations[popID].vmesh->print();
      // Grow the vmesh and block container, if necessary. Try performing this on-device, if possible.
      phiprof::Timer preparationTimer {"GPU resize mesh on-device"};
      resize_vbc_kernel_pre<<<1, 1, 0, stream>>> (
         populations[popID].dev_vmesh,
         populations[popID].dev_blockContainer,
         gpu_list_with_replace_new[cpuThreadID],
         gpu_list_delete[cpuThreadID],
         gpu_list_to_replace[cpuThreadID],
         gpu_list_with_replace_old[cpuThreadID],
         returnLID[cpuThreadID], // return values: nbefore, nafter, nblockstochange, resize success
         returnRealf[cpuThreadID] // mass loss, set to zero
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuMemcpyAsync(&host_returnLID[0], returnLID[cpuThreadID], 4*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, stream) );
      CHK_ERR( gpuStreamSynchronize(stream) );
      // Grow mesh if necessary and on-device resize did not work??
      const vmesh::LocalID nBlocksBeforeAdjust = host_returnLID[0];
      const vmesh::LocalID nBlocksAfterAdjust = host_returnLID[1];
      const vmesh::LocalID nBlocksToChange = host_returnLID[2];
      const vmesh::LocalID resizeDevSuccess = host_returnLID[3];
      preparationTimer.stop();
      if ( (nBlocksAfterAdjust > nBlocksBeforeAdjust) && (resizeDevSuccess == 0)) {
         //GPUTODO is _FACTOR enough instead of _PADDING?
         phiprof::Timer preparationHostTimer {"GPU resize mesh on host"};
         populations[popID].vmesh->setNewCapacity(nBlocksAfterAdjust*BLOCK_ALLOCATION_PADDING);
         populations[popID].vmesh->setNewSize(nBlocksAfterAdjust);
         populations[popID].blockContainer->setNewCapacity(nBlocksAfterAdjust*BLOCK_ALLOCATION_PADDING);
         populations[popID].blockContainer->setNewSize(nBlocksAfterAdjust);
         populations[popID].Upload();
      }

      if (nBlocksToChange==0) {
         // Update vmesh cached size
         populations[popID].vmesh->setNewCachedSize(nBlocksAfterAdjust);
         return nBlocksAfterAdjust;
      }

      phiprof::Timer addRemoveKernelTimer {"GPU add and remove blocks kernel"};
      // Each GPU block / workunit could manage several Vlasiator velocity blocks at once.
      // However, thread syncs inside the kernel prevent this.
      // const uint vlasiBlocksPerWorkUnit = WARPSPERBLOCK * GPUTHREADS / WID3;
      const uint vlasiBlocksPerWorkUnit = 1;
      // ceil int division
      const uint launchBlocks = 1 + ((nBlocksToChange - 1) / vlasiBlocksPerWorkUnit);

      // Third argument specifies the number of bytes in *shared memory* that is
      // dynamically allocated per block for this call in addition to the statically allocated memory.
      //CHK_ERR( gpuStreamSynchronize(stream) );
      update_velocity_blocks_kernel<<<launchBlocks, vlasiBlocksPerWorkUnit * WID3, 0, stream>>> (
         populations[popID].dev_vmesh,
         populations[popID].dev_blockContainer,
         gpu_list_with_replace_new[cpuThreadID],
         gpu_list_delete[cpuThreadID],
         gpu_list_to_replace[cpuThreadID],
         gpu_list_with_replace_old[cpuThreadID],
         nBlocksBeforeAdjust,
         nBlocksToChange,
         nBlocksAfterAdjust,
         returnRealf[cpuThreadID] // mass loss
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuMemcpyAsync(&host_rhoLossAdjust, returnRealf[cpuThreadID], sizeof(Realf), gpuMemcpyDeviceToHost, stream) );

      // Shrink the vmesh and block container, if necessary
      if (nBlocksAfterAdjust < nBlocksBeforeAdjust) {
         // Should not re-allocate on shrinking, so do on-device
         resize_vbc_kernel_post<<<1, 1, 0, stream>>> (
            populations[popID].dev_vmesh,
            populations[popID].dev_blockContainer,
            nBlocksAfterAdjust
            );
         CHK_ERR( gpuPeekAtLastError() );
      }
      // Update vmesh cached size
      populations[popID].vmesh->setNewCachedSize(nBlocksAfterAdjust);

      CHK_ERR( gpuStreamSynchronize(stream) );
      this->populations[popID].RHOLOSSADJUST += host_rhoLossAdjust;
      addRemoveKernelTimer.stop();

      // DEBUG output after kernel
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::LocalID nAll = populations[popID].vmesh->size();
      if (nAll!=nBlocksAfterAdjust) {
         //phiprof::Timer debugTimer {"Vmesh and VBC debug output"};
         populations[popID].vmesh->gpu_prefetchHost();
         CHK_ERR( gpuStreamSynchronize(stream) );
         printf("after kernel, size is %d should be %d\n",nAll,nBlocksAfterAdjust);
         for (vmesh::LocalID m=0; m<nAll; ++m) {
            const vmesh::GlobalID GIDs = populations[popID].vmesh->getGlobalID(m);
            const vmesh::LocalID LIDs = populations[popID].vmesh->getLocalID(GIDs);
            printf("LID %d GID-solved %d LID-solved %d\n",m,GIDs,LIDs);
         }
         populations[popID].vmesh->gpu_prefetchDevice();
      }
      #endif
      return nBlocksAfterAdjust;
   }

   void SpatialCell::adjustSingleCellVelocityBlocks(const uint popID, bool doDeleteEmpty) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      //neighbor_ptrs is empty as we do not have any consistent
      //data in neighbours yet, adjustments done only based on velocity
      //space.
      std::vector<SpatialCell*> neighbor_ptrs;
      update_velocity_block_content_lists(popID);
      adjust_velocity_blocks(neighbor_ptrs,popID,doDeleteEmpty);
   }

   /** Update the two lists containing blocks with content, and blocks without content.
    * @see adjustVelocityBlocks */
   void SpatialCell::update_velocity_block_content_lists(const uint popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
      const gpuStream_t stream = gpu_getStream();
      // phiprof::Timer updateListsTimer {"GPU update spatial cell block lists"};

      phiprof::Timer reservationTimer {"GPU apply reservation"};
      applyReservation(popID);
      reservationTimer.stop();

      phiprof::Timer clearTimer {"GPU clear maps"};
      velocity_block_with_content_list_size = 0;
      velocity_block_with_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwcl_sizePower));
      velocity_block_with_no_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwncl_sizePower));
      CHK_ERR( gpuStreamSynchronize(stream) );
      clearTimer.stop();
      phiprof::Timer sizeTimer {"GPU size"};
      const uint nBlocks = populations[popID].vmesh->size();
      if (nBlocks==0) {
         return;
      }
      CHK_ERR( gpuStreamSynchronize(stream) );
      sizeTimer.stop();
      phiprof::Timer updateListsTimer {"GPU update spatial cell block lists"};
      const Real velocity_block_min_value = getVelocityBlockMinValue(popID);
      // Each GPU block / workunit can manage several Vlasiator velocity blocks at once. (TODO FIX)
      //const uint vlasiBlocksPerWorkUnit = WARPSPERBLOCK * GPUTHREADS / WID3;
      const uint vlasiBlocksPerWorkUnit = 1;
      // ceil int division
      const uint launchBlocks = 1 + ((nBlocks - 1) / vlasiBlocksPerWorkUnit);
      //std::cerr<<launchBlocks<<" "<<vlasiBlocksPerWorkUnit<<" vmesh "<<populations[popID].dev_vmesh<<" vbc "<<populations[popID].dev_blockContainer<<" wcm "<<dev_velocity_block_with_content_map<<" wncm "<<dev_velocity_block_with_no_content_map<<" minval "<<velocity_block_min_value<<std::endl;

      // Third argument specifies the number of bytes in *shared memory* that is
      // dynamically allocated per block for this call in addition to the statically allocated memory.
      //update_velocity_block_content_lists_kernel<<<launchBlocks, WID3, WID3*sizeof(int), stream>>> (
      update_velocity_block_content_lists_kernel<<<launchBlocks, (vlasiBlocksPerWorkUnit * WID3), 0, stream>>> (
         populations[popID].dev_vmesh,
         populations[popID].dev_blockContainer,
         dev_velocity_block_with_content_map,
         dev_velocity_block_with_no_content_map,
         velocity_block_min_value
         );
      CHK_ERR( gpuPeekAtLastError() );

      // Now extract values from the map
      velocity_block_with_content_map->extractAllKeysLoop(*dev_velocity_block_with_content_list,stream);
      // split::SplitInfo info;
      // Hashinator::MapInfo info_m;
      split::SplitInfo info;
      velocity_block_with_content_list->copyMetadata(&info, stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      velocity_block_with_content_list_size = info.size;
      //velocity_block_with_content_list_size = velocity_block_with_content_list->size();
   }

   void SpatialCell::prefetchDevice() {

      for (size_t p=0; p<populations.size(); ++p) {
         populations[p].vmesh->gpu_prefetchDevice();
         populations[p].blockContainer->gpu_prefetchDevice();
      }
   }
   void SpatialCell::prefetchHost() {
      for (size_t p=0; p<populations.size(); ++p) {
         populations[p].vmesh->gpu_prefetchHost();
         populations[p].blockContainer->gpu_prefetchHost();
      }
   }

   /** Get maximum translation timestep for the given species.
    * @param popID ID of the particle species.
    * @return Maximum timestep calculated by the Vlasov translation.*/
   const Real& SpatialCell::get_max_r_dt(const uint popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      return populations[popID].max_dt[species::MAXRDT];
   }

   /** Get maximum acceleration timestep for the given species.
    * @param popID ID of the particle species.
    * @return Maximum timestep calculated by Vlasov acceleration.*/
   const Real& SpatialCell::get_max_v_dt(const uint popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      return populations[popID].max_dt[species::MAXVDT];
   }

   /** Get MPI datatype for sending the cell data.
    * @param cellID Spatial cell (dccrg) ID.
    * @param sender_rank Rank of the MPI process sending data from this cell.
    * @param receiver_rank Rank of the MPI process receiving data to this cell.
    * @param receiving If true, this process is receiving data.
    * @param neighborhood Neighborhood ID.
    * @return MPI datatype that transfers the requested data.*/
   std::tuple<void*, int, MPI_Datatype> SpatialCell::get_mpi_datatype(
                                                                      const CellID cellID,
                                                                      const int sender_rank,
                                                                      const int receiver_rank,
                                                                      const bool receiving,
                                                                      const int neighborhood
      ) {

      std::vector<MPI_Aint> displacements;
      std::vector<int> block_lengths;

      // create datatype for actual data if we are in the first two
      // layers around a boundary, or if we send for the whole system
      // in AMR translation, only send the necessary cells
      if (this->mpiTransferEnabled && ((SpatialCell::mpiTransferAtSysBoundaries==false && SpatialCell::mpiTransferInAMRTranslation==false) ||
                                       (SpatialCell::mpiTransferAtSysBoundaries==true && (this->sysBoundaryLayer ==1 || this->sysBoundaryLayer ==2)) ||
                                       (SpatialCell::mpiTransferInAMRTranslation==true &&
                                        this->parameters[CellParams::AMR_TRANSLATE_COMM_X+SpatialCell::mpiTransferXYZTranslation]==true ))) {

         //add data to send/recv to displacement and block length lists
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE1) != 0) {
            //first copy values in case this is the send operation
            populations[activePopID].N_blocks = populations[activePopID].vmesh->size();

            // send velocity block list size
            displacements.push_back((uint8_t*) &(populations[activePopID].N_blocks) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE2) != 0) {
            // STAGE1 should have been done, otherwise we have problems...
            if (receiving) {
               // Set population size based on mpi_number_of_blocks transferred earlier,
               // Cleared to be ready to receive
               populations[activePopID].vmesh->clear();
               populations[activePopID].vmesh->setNewSize(populations[activePopID].N_blocks);
               // VBC resized in prepare_to_receive_blocks
               populations[activePopID].Upload();
            } else {
               //Ensure N_blocks is still correct
               populations[activePopID].N_blocks = populations[activePopID].vmesh->size();
            }

            // send velocity block list
            displacements.push_back((uint8_t*) &(populations[activePopID].vmesh->getGrid()[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::GlobalID) * populations[activePopID].N_blocks);
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE1) !=0) {
            //Communicate size of list so that buffers can be allocated on receiving side
            // if (!receiving) { // Already done during block evaluation
            //    this->velocity_block_with_content_list_size = velocity_block_with_content_list->size();
            // }
            displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list_size) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE2) !=0) {
            const gpuStream_t stream = gpu_getStream();
            if (receiving) {
               this->velocity_block_with_content_list->resize(this->velocity_block_with_content_list_size,true);
               this->velocity_block_with_content_list_capacity = this->velocity_block_with_content_list->capacity();
               this->dev_velocity_block_with_content_list = this->velocity_block_with_content_list->upload();
               //this->velocity_block_with_content_list->optimizeGPU(stream); // included in upload()
             }
            //velocity_block_with_content_list_size should first be updated, before this can be done (STAGE1)
            displacements.push_back((uint8_t*) this->velocity_block_with_content_list->data() - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::GlobalID)*this->velocity_block_with_content_list_size);
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_DATA) !=0) {
            displacements.push_back((uint8_t*) get_data(activePopID) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Realf) * WID3 * populations[activePopID].blockContainer->size());
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::NEIGHBOR_VEL_BLOCK_DATA) != 0) {
            /*We are actually transferring the data of a
            * neighbor. The values of neighbor_block_data
            * and neighbor_number_of_blocks should be set in
            * solver.*/

            // Send this data only to ranks that contain face neighbors
            // this->neighbor_number_of_blocks has been initialized to 0, on other ranks it can stay that way.
            const set<int>& ranks = this->face_neighbor_ranks[neighborhood];
            if ( P::amrMaxSpatialRefLevel == 0 || receiving || ranks.find(receiver_rank) != ranks.end()) {

               for ( int i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
                  displacements.push_back((uint8_t*) this->neighbor_block_data[i] - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Realf) * WID3 * this->neighbor_number_of_blocks[i]);
               }

            }
         }

         // send  spatial cell parameters
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PARAMETERS)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * CellParams::N_SPATIAL_CELL_PARAMS);
         }

         // send spatial cell dimensions and coordinates
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_DIMENSIONS)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::XCRD]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 6);
         }

         // send  BGBXVOL BGBYVOL BGBZVOL PERBXVOL PERBYVOL PERBZVOL
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBXVOL]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 6);
         }

         // send RHOM, VX, VY, VZ
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHOM_V)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHOM]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 4);
         }

         // send RHOM_DT2, VX_DT2, VY_DT2, VZ_DT2
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHOMDT2_VDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHOM_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 4);
         }

         // send RHOQ
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHOQ)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHOQ]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real));
         }

         // send RHOQ_DT2
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHOQDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHOQ_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real));
         }

         // send  spatial cell BVOL derivatives
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL_DERIVATIVES)!=0){
            displacements.push_back((uint8_t*) &(this->derivativesBVOL[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * bvolderivatives::N_BVOL_DERIVATIVES);
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_IOLOCALCELLID)!=0){
            displacements.push_back((uint8_t*) &(this->ioLocalCellId) - (uint8_t*) this);
            block_lengths.push_back(sizeof(uint64_t));
         }

         // send electron pressure gradient term components
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_GRADPE_TERM)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::EXGRADPE]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }


         // send P tensor diagonal components
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_P)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::P_11]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::P_11_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }

         // send  sysBoundaryFlag
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_SYSBOUNDARYFLAG)!=0){
            displacements.push_back((uint8_t*) &(this->sysBoundaryFlag) - (uint8_t*) this);
            block_lengths.push_back(sizeof(uint));
            displacements.push_back((uint8_t*) &(this->sysBoundaryLayer) - (uint8_t*) this);
            block_lengths.push_back(sizeof(uint));
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_PARAMETERS) !=0) {
            displacements.push_back((uint8_t*) get_block_parameters(activePopID) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * size(activePopID) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
         }
         // Copy particle species metadata
         if ((SpatialCell::mpi_transfer_type & Transfer::POP_METADATA) != 0) {
            for (uint popID=0; popID<populations.size(); ++popID) {
               displacements.push_back((uint8_t*) &(populations[popID].RHO) - (uint8_t*)this);
               block_lengths.push_back(offsetof(spatial_cell::Population, N_blocks));
            }
         }

         // Refinement parameters
         if ((SpatialCell::mpi_transfer_type & Transfer::REFINEMENT_PARAMETERS)){
            displacements.push_back(reinterpret_cast<uint8_t*>(this->parameters.data() + CellParams::AMR_ALPHA1) - reinterpret_cast<uint8_t*>(this));
            block_lengths.push_back(sizeof(Real) * (CellParams::AMR_ALPHA2 - CellParams::AMR_ALPHA1 + 1)); // This is just 2, but let's be explicit
         }

      }

      void* address = this;
      int count;
      MPI_Datatype datatype;

      if (displacements.size() > 0) {
         count = 1;
         MPI_Type_create_hindexed(
            displacements.size(),
            &block_lengths[0],
            &displacements[0],
            MPI_BYTE,
            &datatype
         );
      } else {
         count = 0;
         datatype = MPI_BYTE;
      }

      const bool printMpiDatatype = false;
      if(printMpiDatatype) {
         int mpiSize;
         int myRank;
         MPI_Type_size(datatype,&mpiSize);
         MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
         cout << myRank << " get_mpi_datatype: " << cellID << " " << sender_rank << " " << receiver_rank << " " << mpiSize << ", Nblocks = " << populations[activePopID].N_blocks << ", nbr Nblocks =";
         for (uint i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
            const set<int>& ranks = this->face_neighbor_ranks[neighborhood];
            if ( receiving || ranks.find(receiver_rank) != ranks.end()) {
               cout << " " << this->neighbor_number_of_blocks[i];
            } else {
               cout << " " << 0;
            }
         }
         cout << " face_neighbor_ranks =";
         for (const auto& rank : this->face_neighbor_ranks[neighborhood]) {
            cout << " " << rank;
         }
         cout << endl;
      }

      return std::make_tuple(address,count,datatype);
   }

  /**< Minimum value of distribution function in any phase space cell
    * of a velocity block for the block to be considered to have content.
    * @param popID ID of the particle species.
    * @return Sparse min value for this species.*/
   Real SpatialCell::getVelocityBlockMinValue(const uint popID) const {
      return populations[popID].velocityBlockMinValue;
   }

   /** Prepares this spatial cell to receive the velocity grid over MPI.
    * At this stage we have received a new block list over MPI into
    * mpi_velocity_block_list, but the rest of the cell structures
    * have not been adapted to this new list. Here we re-initialize
    * the cell with empty blocks based on the new list.*/
   void SpatialCell::prepare_to_receive_blocks(const uint popID) {
      phiprof::Timer setGridTimer {"GPU receive blocks: set grid"};
      populations[popID].vmesh->setGrid(); // Based on localToGlobalMap
      const vmesh::LocalID meshSize = populations[popID].vmesh->size();
      populations[popID].blockContainer->setNewSize(meshSize);
      populations[popID].Upload();
      // Set velocity block parameters:
      const gpuStream_t stream = gpu_getStream();
      CHK_ERR( gpuStreamSynchronize(stream) );
      if (meshSize>0) {
         // ceil int division
         uint launchBlocks = 1 + ((meshSize - 1) / (WARPSPERBLOCK*GPUTHREADS));
         update_blockparameters_kernel<<<launchBlocks, (WARPSPERBLOCK*GPUTHREADS), 0, stream>>> (
            populations[popID].dev_vmesh,
            populations[popID].dev_blockContainer,
            meshSize
            );
         CHK_ERR( gpuPeekAtLastError() );
         CHK_ERR( gpuStreamSynchronize(stream) );
      }
   }

   /** Set the particle species SpatialCell should use in functions that
    * use the velocity mesh.
    * @param popID Population ID.
    * @return If true, the new species is in use.*/
   bool SpatialCell::setCommunicatedSpecies(const uint popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= getObjectWrapper().particleSpecies.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds species.size() " << getObjectWrapper().particleSpecies.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      activePopID = popID;
      return true;
   }

   /** Set maximum translation timestep for a particle species.
    * This function is called during Vlasov translation.
    * @param popID ID of the particle species.
    * @param value New maximum timestep.*/
   void SpatialCell::set_max_r_dt(const uint popID,const Real& value) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      populations[popID].max_dt[species::MAXRDT] = value;
   }

   /** Set maximum acceleration timestep for a particle species.
    * This function is called during Vlasov acceleration.
    * @param popID ID of the particle species.
    * @param value New maximum timestep.*/
   void SpatialCell::set_max_v_dt(const uint popID,const Real& value) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      populations[popID].max_dt[species::MAXVDT] = value;
   }

   /**  Purges extra capacity from block vectors. It sets size to
    * num_blocks * block_allocation_factor (if capacity greater than this),
    * and also forces capacity to this new smaller value.
    * @return True on success.*/
   bool SpatialCell::shrink_to_fit() {
      bool success = true;
      return success;

      for (size_t p=0; p<populations.size(); ++p) {
         const uint64_t amount
            = 2 + populations[p].blockContainer->size()
            * populations[p].blockContainer->getBlockAllocationFactor();

         // Allow capacity to be a bit large than needed by number of blocks, shrink otherwise
         if (populations[p].blockContainer->capacity() > amount ) {
            if (populations[p].blockContainer->setNewCapacity(amount) == false) {
               success = false;
            }
         }
      }
      return success;
   }
   void SpatialCell::printMeshSizes() {
      cerr << "SC::printMeshSizes:" << endl;
      for (size_t p=0; p<populations.size(); ++p) {
         cerr << "\t pop " << p << " " << populations[p].vmesh->size() << ' ' << populations[p].blockContainer->size() << endl;
      }
   }

   /** Updates minValue based on algorithm value from parameters (see parameters.cpp).
    * @param popID ID of the particle species.*/
   void SpatialCell::updateSparseMinValue(const uint popID) {

      species::Species& population = getObjectWrapper().particleSpecies[popID];

      if ( population.sparseDynamicAlgorithm == 1 || population.sparseDynamicAlgorithm == 2 ) {
         // Linear algorithm for the minValue: y=kx+b
         const Real k = (population.sparseDynamicMinValue2 - population.sparseDynamicMinValue1) / (population.sparseDynamicBulkValue2 - population.sparseDynamicBulkValue1);
         const Real b = population.sparseDynamicMinValue1 - k * population.sparseDynamicBulkValue1;
         Real x;
         if ( population.sparseDynamicAlgorithm == 1 ) {
            x = this->populations[popID].RHO;
         } else {
            x = this->get_number_of_velocity_blocks(popID);
         }
         const Real newMinValue = k*x+b;
         if( newMinValue < population.sparseDynamicMinValue1 ) { // Compare against the min minValue
            populations[popID].velocityBlockMinValue = population.sparseDynamicMinValue1;
         } else if( newMinValue > population.sparseDynamicMinValue2 ) { // Compare against the max minValue
            populations[popID].velocityBlockMinValue = population.sparseDynamicMinValue2;
         } else {
            populations[popID].velocityBlockMinValue = newMinValue;
         }
         return;
      } else {
         populations[popID].velocityBlockMinValue = getObjectWrapper().particleSpecies[popID].sparseMinValue;
         return;
      }
   }

} // namespace spatial_cell
