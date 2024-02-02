/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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

#ifdef DEBUG_VLASIATOR
   #define DEBUG_SPATIAL_CELL
#endif

using namespace std;

/** GPU kernel for identifying which blocks have relevant content */
__global__ void __launch_bounds__(WID3,4) update_velocity_block_content_lists_kernel (
   vmesh::VelocityMesh *vmesh,
   const uint nBlocks,
   vmesh::VelocityBlockContainer *blockContainer,
   vmesh::GlobalID* vbwcl_gather,
   vmesh::GlobalID* vbwncl_gather,
   Real velocity_block_min_value
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;
   __shared__ int has_content[WID3];
   //const uint nBlocks = vmesh->size();
   const vmesh::GlobalID invalidGID = vmesh->invalidGlobalID();
   for (uint blockLID=blocki; blockLID<nBlocks; blockLID += gpuBlocks) {
      const vmesh::GlobalID blockGID = vmesh->getGlobalID(blockLID);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == vmesh->invalidGlobalID()) {
         if (ti==0) printf("Invalid GID encountered in update_velocity_block_content_lists_kernel!\n");
         __syncthreads();
         continue;
      }
      if (blockLID == vmesh->invalidLocalID()) {
         if (ti==0) printf("Invalid LID encountered in update_velocity_block_content_lists_kernel!\n");
         __syncthreads();
         continue;
      }
      #endif
      // Check each velocity cell if it is above the threshold
      const Realf* avgs = blockContainer->getData(blockLID);
      has_content[ti] = avgs[ti] >= velocity_block_min_value ? 1 : 0;
      __syncthreads();
      // Implemented just a simple non-optimized thread OR
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (ti < s) {
            has_content[ti] = has_content[ti] || has_content[ti + s];
         }
         __syncthreads();
      }
      // Increment vector only from thread zero
      if (ti==0) {
         if (has_content[0]) {
            // velocity_block_with_content_list->device_push_back(blockGID);
            vbwcl_gather[blockLID] = blockGID;
            vbwncl_gather[blockLID] = invalidGID;
         } else {
            // velocity_block_with_no_content_list->device_push_back(blockGID);
            vbwncl_gather[blockLID] = blockGID;
            vbwcl_gather[blockLID] = invalidGID;
         }
      }
      __syncthreads();
   }
}

/** Gpu Kernel to quickly gather blocks and their v-space halo */
__global__ void __launch_bounds__(GPUTHREADS,4) update_blocks_required_halo_kernel (
   vmesh::VelocityMesh *vmesh,
   vmesh::GlobalID* velocity_block_with_content_list_data,
   vmesh::GlobalID* blocks_required_list,
   const vmesh::LocalID localContentBlocks,
   const int addWidthV
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int warpSize = blockDim.x*blockDim.y*blockDim.z; // should be 27
   const int offset_vx = (int)threadIdx.x - addWidthV; // expected block size 3x3x3
   const int offset_vy = (int)threadIdx.y - addWidthV; // and addWidthV 1
   const int offset_vz = (int)threadIdx.z - addWidthV;
   const vmesh::LocalID ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (vmesh::LocalID index=blocki; index<localContentBlocks; index += gpuBlocks) {
      const vmesh::GlobalID GID = velocity_block_with_content_list_data[index];
      vmesh::LocalID ind0,ind1,ind2;
      vmesh->getIndices(GID,ind0,ind1,ind2);
      const int nind0 = ind0 + offset_vx;
      const int nind1 = ind1 + offset_vy;
      const int nind2 = ind2 + offset_vz;
      const vmesh::GlobalID nGID
         = vmesh->getGlobalID(nind0,nind1,nind2);
      blocks_required_list[index*warpSize + ti] = nGID;
   } // for blocks
}

/** GPU kernel for identifying which blocks need to be moved from end of vspace to earlier positions.
    This kernel now searches through the end of the v-space, and checks which of
    those blocks are on the deleted blocks list.
 */
__global__ void __launch_bounds__(GPUTHREADS,4) update_blocks_to_move_kernel_2 (
   const split::SplitVector<vmesh::GlobalID>* localToGlobalMap,
   const vmesh::LocalID nBlocksCurrent,
   const vmesh::LocalID nBlocksRequired,
   const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* BlocksDeleteMap,
   vmesh::GlobalID* BlocksToMove_buffer
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   //const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   const uint nEvaluate = nBlocksCurrent - nBlocksRequired; // assumed to be positive
   for (uint index=blocki; index<nEvaluate; index+=gpuBlocks) {
      // Evaluate all GIDs at tail end of vmesh
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID GID = localToGlobalMap->at(nBlocksRequired + index);
      assert((GID != vmesh::INVALID_LOCALID) && "invalid GID in update_blocks_to_move_kernel_2");
      #else
      const vmesh::GlobalID GID = (*localToGlobalMap)[nBlocksRequired + index];
      #endif
      // Search for this GID in the BlocksDeleteMap
      vmesh::LocalID retval = vmesh::INVALID_LOCALID;
      BlocksDeleteMap->warpFind(GID, retval, ti % GPUTHREADS);
      if (ti==0) {
         if (retval == vmesh::INVALID_LOCALID) {
            // Block does not exist in deletion set/map, thus must be moved
            BlocksToMove_buffer[index] = GID;
         } else {
            // Block exists in map/set of GIDs to delete - do not move
            BlocksToMove_buffer[index] = vmesh::INVALID_GLOBALID;
         }
      }
      __syncthreads();
   }
   // After this kernel is run, we perform streamcompaction on BlocksToMove_buffer.
}

/** GPU kernel for identifying which blocks need to be moved from end of vspace to earlier positions.
    This kernel may be non-optimized in itself, but use of it gets rid
    of the need of vmesh prefetching back and forth.
 */
__global__ void __launch_bounds__(GPUTHREADS,4) update_blocks_to_move_kernel (
   vmesh::VelocityMesh *vmesh,
   split::SplitVector<vmesh::GlobalID>* BlocksRequired,
   split::SplitVector<vmesh::GlobalID>* BlocksToMove,
   const uint nBlocksRequired
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   //const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (uint index=blocki; index<nBlocksRequired; index += gpuBlocks) {
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID GID = BlocksRequired->at(index);
      assert((GID != vmesh->invalidGlobalID()) && "invalid GID in update_blocks_to_move_kernel");
      #else
      const vmesh::GlobalID GID = (*BlocksRequired)[index];
      #endif
      const vmesh::LocalID LID = vmesh->warpGetLocalID(GID, ti);
      if (ti==0) {
         if ( (LID!=vmesh->invalidLocalID()) && (LID>=nBlocksRequired)) {
            // Block exists but within region of vmesh which shall be deleted - queue for moving
            BlocksToMove->device_push_back(GID);
         }
      }
      __syncthreads();
   }
}

/** GPU kernel for making sure that in the list of blocks to remove,
    any blocks which need to be replaced with other ones (having LID < nAfterAdjust)
    are at the start of the list.
 */
__global__ void __launch_bounds__(GPUTHREADS,4) sort_blocks_remove_kernel (
   const vmesh::VelocityMesh *vmesh,
   split::SplitVector<vmesh::GlobalID>* BlocksToRemove,
   split::SplitVector<vmesh::GlobalID>* BlocksList,
   const vmesh::LocalID nBlocksToRemove,
   const vmesh::LocalID nBlocksAfterAdjust
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   //const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (uint index=blocki; index<nBlocksToRemove; index += gpuBlocks) {
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID GID = BlocksToRemove->at(index);
      assert((GID != vmesh::INVALID_GLOBALID) && "invalid GID in sort_blocks_remove_kernel");
      #else
      const vmesh::GlobalID GID = (*BlocksToRemove)[index];
      #endif
      // Find the LID of this GID
      const vmesh::LocalID LID = vmesh->warpGetLocalID(GID, ti);
      #ifdef DEBUG_SPATIAL_CELL
      assert((LID != vmesh::INVALID_LOCALID) && "invalid LID in sort_blocks_remove_kernel");
      #endif
      if (ti==0) {
         if (LID < nBlocksAfterAdjust) {
            // goes in beginning
            (*BlocksList)[index] = GID;
            (*BlocksList)[nBlocksToRemove+index] = vmesh::INVALID_GLOBALID;
         } else {
            // goes at end
            (*BlocksList)[nBlocksToRemove+index] = GID;
            (*BlocksList)[index] = vmesh::INVALID_GLOBALID;
         }
      }
      __syncthreads();
   }
   // After this kernel is run, we perform streamcompaction on BlocksToMove_buffer.
}

/** GPU kernel for quickly filling block parameters.
 */
__global__ void __launch_bounds__(GPUTHREADS,4) update_blockparameters_kernel (
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   Real *blockParameters,
   vmesh::LocalID nLIDs
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (uint index=blocki*warpSize; index<nLIDs; index += gpuBlocks*warpSize) {
      const vmesh::LocalID LID = index+ti;
      if (LID < nLIDs) {
         // Set velocity block parameters:
         const vmesh::GlobalID GID = vmesh->getGlobalID(LID);
         // Write in block parameters
         vmesh->getBlockInfo(GID, blockParameters + LID * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD);
      }
   }
}

/** GPU kernel for updating blocks based on generated lists */
__global__ void __launch_bounds__(WID3,4) update_velocity_blocks_kernel(
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   split::SplitVector<vmesh::GlobalID>* BlocksToAdd,
   split::SplitVector<vmesh::GlobalID>* BlocksToRemove,
   split::SplitVector<vmesh::GlobalID>* BlocksToMove,
   vmesh::LocalID nBlocksBeforeAdjust,
   vmesh::LocalID nBlocksAfterAdjust,
   Realf* gpu_rhoLossAdjust
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const uint ti = threadIdx.z*WID2 + threadIdx.y*WID + threadIdx.x;

   const vmesh::LocalID nToAdd = BlocksToAdd->size();
   const vmesh::LocalID nToRemove = BlocksToRemove->size();
   const vmesh::LocalID nToMove = BlocksToMove->size();
   // For tracking mass-loss
   Realf local_rhoLoss = 0;
   __shared__ Realf massloss[WID3];

   // We loop either up to nToRemove or or nToAdd, whichever is larger.
   const vmesh::LocalID loopExtent = nToAdd > nToRemove ? nToAdd : nToRemove;

   for (vmesh::LocalID m=blocki; m<loopExtent; m += gpuBlocks) {
      if (m < nToRemove) {
         #ifdef DEBUG_SPATIAL_CELL
         const vmesh::GlobalID rmGID = BlocksToRemove->at(m);
         #else
         const vmesh::GlobalID rmGID = (*BlocksToRemove)[m];
         #endif
         const vmesh::LocalID rmLID = vmesh->warpGetLocalID(rmGID,ti);
         #ifdef DEBUG_SPATIAL_CELL
         if (rmGID == vmesh->invalidGlobalID()) {
            if (rmLID != vmesh->invalidLocalID()) {
               // Valid LID but invalid GID: only remove from vmesh localToGlobal?
               printf("Removing blocks: Valid LID %u but invalid GID!\n",rmLID);
            } else {
               printf("Removing blocks: Invalid LID and GID!\n");
            }
            continue;
         }
         if (rmLID == vmesh->invalidLocalID()) {
            if (rmGID != vmesh->invalidGlobalID()) {
               // Valid GID but invalid LID: only remove from vmesh globalToLocal?
               printf("Removing blocks: Valid GID %ul but invalid LID!\n",rmGID);
            }
            continue;
         }
         assert(rmLID < nBlocksBeforeAdjust && "Trying to remove block which has LID >= nBlocksBeforeAdjust!");
         #endif
         // Track mass loss:
         Realf* rm_avgs = blockContainer->getData(rmLID);
         Real* rm_block_parameters = blockContainer->getParameters(rmLID);
         const Real rm_DV3 = rm_block_parameters[BlockParams::DVX]
            * rm_block_parameters[BlockParams::DVY]
            * rm_block_parameters[BlockParams::DVZ];

         // thread-sum for rho
         massloss[ti] = rm_avgs[ti]*rm_DV3;
         __syncthreads();
         // Implemented just a simple non-optimized thread sum
         for (unsigned int s=WID3/2; s>0; s>>=1) {
            if (ti < s) {
               massloss[ti] += massloss[ti + s];
            }
            __syncthreads();
         }

         if (ti==0) {
            // Bookkeeping only by one thread
            local_rhoLoss += massloss[0];
         }
         __syncthreads();

         if (m < nToAdd) {
            // Replace with new value
            #ifdef DEBUG_SPATIAL_CELL
            const vmesh::GlobalID addGID = BlocksToAdd->at(m);
            #else
            const vmesh::GlobalID addGID = (*BlocksToAdd)[m];
            #endif
            rm_avgs[ti] = 0;
            if (ti==0) {
               // Write in block parameters
               vmesh->getBlockInfo(addGID, rm_block_parameters+BlockParams::VXCRD);
            }
            __syncthreads();
            vmesh->warpReplaceBlock(rmGID,rmLID,addGID,ti);
            #ifdef DEBUG_SPATIAL_CELL
            if (vmesh->getGlobalID(rmLID) == vmesh->invalidGlobalID()) {
               if (ti==0) printf("Error! Adding resulted in invalid Global ID in update_velocity_blocks_kernel! \n");
               __syncthreads();
               continue;
            }
            if (vmesh->getLocalID(addGID) == vmesh->invalidLocalID()) {
               if (ti==0) printf("Error! Adding resulted in invalid Local ID in update_velocity_blocks_kernel! \n");
               __syncthreads();
               continue;
            }
            #endif
            continue;
         } else if (m >= nToAdd && (m < (nToAdd+nToMove))) {
            // GPUTODO: All these need to have rmLID<nAfterAdjust! So we need to get a semi-sorted list of BlocksToRemove.

            // Move replacement from end
            #ifdef DEBUG_SPATIAL_CELL
            const vmesh::GlobalID replaceGID = BlocksToMove->at(m - nToAdd);
            #else
            const vmesh::GlobalID replaceGID = (*BlocksToMove)[m - nToAdd];
            #endif
            const vmesh::LocalID replaceLID = vmesh->warpGetLocalID(replaceGID,ti);

            Realf* repl_avgs = blockContainer->getData(replaceLID);
            Real*  repl_block_parameters = blockContainer->getParameters(replaceLID);
            rm_avgs[ti] = repl_avgs[ti];
            if (ti < BlockParams::N_VELOCITY_BLOCK_PARAMS) {
               rm_block_parameters[ti] = repl_block_parameters[ti];
            }
            __syncthreads();

            // Remove hashmap entry for removed block, add instead created block
            vmesh->warpReplaceBlock(rmGID,rmLID,replaceGID,ti);

            #ifdef DEBUG_SPATIAL_CELL
            if (vmesh->getGlobalID(rmLID) == vmesh->invalidGlobalID()) {
               if (ti==0) printf("Invalid GID encountered in update_velocity_blocks_kernel!\n");
               __syncthreads();
               continue;
            }
            if (vmesh->getLocalID(replaceGID) == vmesh->invalidLocalID()) {
               if (ti==0) printf("Invalid LID encountered in update_velocity_blocks_kernel!\n");
               __syncthreads();
               continue;
            }
            #endif
            continue;
         } else {
            // This LID so large that we just delete without replacing
            #ifdef DEBUG_SPATIAL_CELL
            assert(rmLID >= nBlocksAfterAdjust && "Trying to remove block which has LID smaller than nBlocksAfterAdjust!");
            #endif
            vmesh->warpDeleteBlock(rmGID,rmLID,ti);
            continue;
         }
      } else if (m < nToAdd && m >= nToRemove) {
         // Add new block
         #ifdef DEBUG_SPATIAL_CELL
         const vmesh::GlobalID addGID = BlocksToAdd->at(m);
         if (ti==0) {
            if (vmesh->getLocalID(addGID) != vmesh->invalidLocalID()) {
               printf("Trying to add new GID %u to mesh which already contains it! m=%u nBlocksBeforeAdjust=%u nBlocksAfterAdjust=%u nToAdd=%u nToRemove=%u nToMove=%u\n",addGID,m,nBlocksBeforeAdjust,nBlocksAfterAdjust,nToAdd,nToRemove,nToMove);
            }
         }
         //assert(( == vmesh->invalidLocalID()) && "Error! Trying to add new GID to mesh which already contains it!");
         #else
         const vmesh::GlobalID addGID = (*BlocksToAdd)[m];
         #endif

         // We need to add the data of addGID to a new LID:
         const vmesh::LocalID addLID = nBlocksBeforeAdjust + m - nToRemove;
         Realf* add_avgs = blockContainer->getData(addLID);
         #ifdef DEBUG_SPATIAL_CELL
         // Debug check: if we are adding elements, then nToMove should be zero
         assert((nToMove==0) && "nToMove should be zero when adding blocks!");
         assert((addGID != vmesh->invalidGlobalID()) && "Error! Trying to add invalid GID!");
         assert((addLID != vmesh->invalidLocalID()) && "Error! Trying to add GID to invalid LID position!");
         #endif
         Real* add_block_parameters = blockContainer->getParameters(addLID);
         // Zero out blockdata
         add_avgs[ti] = 0;
         if (ti==0) {
            // Write in block parameters
            vmesh->getBlockInfo(addGID, add_block_parameters+BlockParams::VXCRD);
         }
         __syncthreads();
         vmesh->warpPlaceBlock(addGID,addLID,ti);
         #ifdef DEBUG_SPATIAL_CELL
         assert((vmesh->getGlobalID(addLID) != vmesh->invalidGlobalID()) && "Error! Trying to add invalid GID!");
         assert((vmesh->getLocalID(addGID) != vmesh->invalidLocalID()) && "Error! Trying to add GID to invalid LID position!");
         #endif
         continue;
      }
      #ifdef DEBUG_SPATIAL_CELL
      // Fall-through error!
      if (ti==0) {
         printf("Error! Fall through in update_velocity_blocks_kernel! nToAdd %u nToRemove %u nToMove %u nBlocksBeforeAdjust %u nBlocksAfterAdjust %u \n",
                nToAdd,nToRemove,nToMove,nBlocksBeforeAdjust,nBlocksAfterAdjust);
      }
      __syncthreads();
      assert(0 && "Error! Fall through in update_velocity_blocks_kernel!");
      #endif
   }
   // Atomically update accumulated mass loss
   if (ti==0) {
      Realf old = atomicAdd(gpu_rhoLossAdjust, local_rhoLoss);
   }
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

      // SplitVectors via pointers for unified memory
      velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(1);
      velocity_block_with_no_content_list = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksRequired = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToAdd = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToRemove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToMove = new split::SplitVector<vmesh::GlobalID>(1);
      velocity_block_with_content_list->clear();
      velocity_block_with_no_content_list->clear();
      BlocksRequired->clear();
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();
      attachedStream=0;
      velocity_block_with_content_list_size=0;
      gpu_velocity_block_with_content_list_buffer=0;
      BlocksRequiredMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      BlocksDeleteMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      velocity_block_with_content_list_capacity=1;
      velocity_block_with_no_content_list_capacity=1;
      BlocksRequired_capacity=1;
      BlocksToAdd_capacity=1;
      BlocksToRemove_capacity=1;
      BlocksToMove_capacity=1;
      BlocksRequiredMap_sizepower=7;
      BlocksDeleteMap_sizepower=7;
   }

   SpatialCell::~SpatialCell() {
      gpu_destructor();
   }

   void SpatialCell::gpu_destructor() {
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         populations[popID].gpu_destructor();
      }
      if (velocity_block_with_content_list) {
         delete velocity_block_with_content_list;
         velocity_block_with_content_list = 0;
      }
      if (velocity_block_with_no_content_list) {
         delete velocity_block_with_no_content_list;
         velocity_block_with_no_content_list = 0;
      }
      if (BlocksRequired) {
         delete BlocksRequired;
         BlocksRequired = 0;
      }
      if (BlocksToAdd) {
         delete BlocksToAdd;
         BlocksToAdd = 0;
      }
      if (BlocksToRemove) {
         delete BlocksToRemove;
         BlocksToRemove = 0;
      }
      if (BlocksToMove) {
         delete BlocksToMove;
         BlocksToMove = 0;
      }
      if (BlocksRequiredMap) {
         delete BlocksRequiredMap;
         BlocksRequiredMap = 0;
      }
      if (BlocksDeleteMap) {
         delete BlocksDeleteMap;
         BlocksDeleteMap = 0;
      }
      if (gpu_velocity_block_with_content_list_buffer) {
         CHK_ERR( gpuFree(gpu_velocity_block_with_content_list_buffer) );
         gpu_velocity_block_with_content_list_buffer = 0;
      }
      velocity_block_with_content_list_capacity=0;
      velocity_block_with_no_content_list_capacity=0;
      BlocksRequired_capacity=0;
      BlocksToAdd_capacity=0;
      BlocksToRemove_capacity=0;
      BlocksToMove_capacity=0;
      BlocksRequiredMap_sizepower=0;
      BlocksDeleteMap_sizepower=0;
   }

   SpatialCell::SpatialCell(const SpatialCell& other) {
      velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(1);
      velocity_block_with_no_content_list = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksRequired = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToAdd = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToRemove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToMove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksRequired->clear();
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();
      velocity_block_with_content_list->clear();
      velocity_block_with_no_content_list->clear();

      BlocksRequiredMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      BlocksDeleteMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);

      // Make space reservation guesses based on popID 0
      const uint reserveSize = other.populations[0].vmesh->size()*BLOCK_ALLOCATION_PADDING;
      BlocksRequired->reserve(reserveSize,true);
      BlocksToAdd->reserve(reserveSize,true);
      BlocksToRemove->reserve(reserveSize,true);
      BlocksToMove->reserve(reserveSize,true);
      velocity_block_with_content_list->reserve(reserveSize,true);
      velocity_block_with_no_content_list->reserve(reserveSize,true);

      velocity_block_with_content_list_capacity=reserveSize;
      velocity_block_with_no_content_list_capacity=reserveSize;
      BlocksRequired_capacity=reserveSize;
      BlocksToAdd_capacity=reserveSize;
      BlocksToRemove_capacity=reserveSize;
      BlocksToMove_capacity=reserveSize;
      BlocksRequiredMap_sizepower=7;
      BlocksDeleteMap_sizepower=7;

      // Member variables
      ioLocalCellId = other.ioLocalCellId;
      sysBoundaryFlag = other.sysBoundaryFlag;
      sysBoundaryLayer = other.sysBoundaryLayer;
      sysBoundaryLayerNew = other.sysBoundaryLayerNew;
      velocity_block_with_content_list_size = other.velocity_block_with_content_list_size;
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
      attachedStream=0;
      gpu_velocity_block_with_content_list_buffer=0;
   }
   const SpatialCell& SpatialCell::operator=(const SpatialCell& other) {
      const uint reserveSize = (other.BlocksRequired)->capacity();
      BlocksRequired->clear();
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();
      velocity_block_with_content_list->clear();
      velocity_block_with_no_content_list->clear();
      delete BlocksRequiredMap;
      delete BlocksDeleteMap;

      BlocksRequired->reserve(reserveSize,true);
      BlocksToAdd->reserve(reserveSize,true);
      BlocksToRemove->reserve(reserveSize,true);
      BlocksToMove->reserve(reserveSize,true);
      velocity_block_with_content_list->reserve(reserveSize,true);
      velocity_block_with_no_content_list->reserve(reserveSize,true);

      const int reserveSize2 = reserveSize > 0 ? reserveSize : 1;
      vmesh::LocalID HashmapReqSize = ceil(log2(reserveSize2)) +2;
      HashmapReqSize = HashmapReqSize > 7 ? HashmapReqSize : 7;
      BlocksRequiredMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(HashmapReqSize);
      BlocksDeleteMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(HashmapReqSize);

      velocity_block_with_content_list_capacity=reserveSize;
      velocity_block_with_no_content_list_capacity=reserveSize;
      BlocksRequired_capacity=reserveSize;
      BlocksToAdd_capacity=reserveSize;
      BlocksToRemove_capacity=reserveSize;
      BlocksToMove_capacity=reserveSize;
      BlocksRequiredMap_sizepower=HashmapReqSize;
      BlocksDeleteMap_sizepower=HashmapReqSize;

      // Member variables
      ioLocalCellId = other.ioLocalCellId;
      sysBoundaryFlag = other.sysBoundaryFlag;
      sysBoundaryLayer = other.sysBoundaryLayer;
      sysBoundaryLayerNew = other.sysBoundaryLayerNew;
      velocity_block_with_content_list_size = other.velocity_block_with_content_list_size;
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

      attachedStream=0;
      return *this;
   }

   /** Advises unified memory subsystem on preferred location of memory
       gpuMemAdviseSetPreferredLocation
       gpuMemAdviseUnsetPreferredLocation
       gpuMemAdviseSetReadMostly
       gpuMemAdviceUnsetReadMostly
       gpuMemAdviseSetAccessedBy
       gpuMemAdviseUnsetAccessedBy
    */
   void SpatialCell::gpu_advise() {
      return;
      // AMD advise is slow
      int device = gpu_getDevice();
      gpuStream_t stream = gpu_getStream();
      BlocksRequired->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksToAdd->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksToRemove->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksToMove->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      velocity_block_with_content_list->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      velocity_block_with_no_content_list->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksRequiredMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksDeleteMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);

      BlocksRequired->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksToAdd->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksToRemove->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksToMove->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      velocity_block_with_content_list->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      velocity_block_with_no_content_list->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksRequiredMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksDeleteMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);

      // Loop over populations
      for (size_t p=0; p<populations.size(); ++p) {
         populations[p].blockContainer->gpu_memAdvise(device,stream);
         populations[p].vmesh->gpu_memAdvise(device,stream);
      }
   }

   /** Attaches or deattaches unified memory to a GPU stream
       When attached, a stream can access this unified memory without
       issues.
    */
   void SpatialCell::gpu_attachToStream(gpuStream_t stream) {
      // Return if attaching is not needed
      if (!needAttachedStreams) {
         return;
      }
      // Attach unified memory regions to streams
      gpuStream_t newStream;
      if (stream==0) {
         newStream = gpu_getStream();
      } else {
         newStream = stream;
      }
      if (newStream == attachedStream) {
         return;
      } else {
         attachedStream = newStream;
      }
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,velocity_block_with_content_list, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,velocity_block_with_no_content_list, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToRemove, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToAdd, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToMove, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksRequired, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksRequiredMap, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksDeleteMap, 0,gpuMemAttachSingle) );
      // Loop over populations
      for (size_t p=0; p<populations.size(); ++p) {
         populations[p].blockContainer->gpu_attachToStream(attachedStream);
         populations[p].vmesh->gpu_attachToStream(attachedStream);
      }
      // Also call attach functions on all splitvectors and hashmaps
      velocity_block_with_content_list->streamAttach(attachedStream);
      velocity_block_with_no_content_list->streamAttach(attachedStream);
      BlocksToRemove->streamAttach(attachedStream);
      BlocksToAdd->streamAttach(attachedStream);
      BlocksToMove->streamAttach(attachedStream);
      BlocksRequired->streamAttach(attachedStream);
      BlocksRequiredMap->streamAttach(attachedStream);
      BlocksDeleteMap->streamAttach(attachedStream);
      return;
   }
   void SpatialCell::gpu_detachFromStream() {
      // Return if attaching is not needed
      if (!needAttachedStreams) {
         return;
      }
      if (attachedStream == 0) {
         // Already detached
         return;
      }
      attachedStream = 0;
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,velocity_block_with_content_list, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,velocity_block_with_no_content_list, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToRemove, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToAdd, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToMove, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksRequired, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksRequiredMap, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksDeleteMap, 0,gpuMemAttachGlobal) );
      // Loop over populations
      for (size_t p=0; p<populations.size(); ++p) {
         populations[p].blockContainer->gpu_detachFromStream();
         populations[p].vmesh->gpu_detachFromStream();
      }
      // Also call detach functions on all splitvectors and hashmaps
      velocity_block_with_content_list->streamAttach(0,gpuMemAttachGlobal);
      velocity_block_with_no_content_list->streamAttach(0,gpuMemAttachGlobal);
      BlocksToRemove->streamAttach(0,gpuMemAttachGlobal);
      BlocksToAdd->streamAttach(0,gpuMemAttachGlobal);
      BlocksToMove->streamAttach(0,gpuMemAttachGlobal);
      BlocksRequired->streamAttach(0,gpuMemAttachGlobal);
      BlocksRequiredMap->streamAttach(0,gpuMemAttachGlobal);
      BlocksDeleteMap->streamAttach(0,gpuMemAttachGlobal);
      return;
   }

   /** Sends the contents of velocity_block_with_content_list into a device buffer so that it can be accessed
       from several streams at once.
    */
   void SpatialCell::gpu_uploadContentLists() {
      //phiprof::Timer timer {"Upload local content lists"};
      gpuStream_t stream = gpu_getStream();
      if (velocity_block_with_content_list_size==0) {
         return;
      }
      CHK_ERR( gpuMallocAsync((void**)&gpu_velocity_block_with_content_list_buffer, velocity_block_with_content_list_size*sizeof(vmesh::LocalID), stream) );
      CHK_ERR( gpuMemcpyAsync(gpu_velocity_block_with_content_list_buffer, velocity_block_with_content_list->data(), velocity_block_with_content_list_size*sizeof(vmesh::LocalID), gpuMemcpyDeviceToDevice, stream) );
      SSYNC;
   }
   /** Clears the device buffer for velocity_block_with_content_list
    */
   void SpatialCell::gpu_clearContentLists() {
      gpuStream_t stream = gpu_getStream();
      if (velocity_block_with_content_list_size==0) {
         return;
      }
      CHK_ERR( gpuFreeAsync(gpu_velocity_block_with_content_list_buffer, stream) );
      gpu_velocity_block_with_content_list_buffer = 0;
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
      size_t reserveSize = populations[popID].reservation * BLOCK_ALLOCATION_FACTOR;
      size_t newReserve = populations[popID].reservation * BLOCK_ALLOCATION_PADDING;
      gpuStream_t stream = gpu_getStream();
      // Now uses host-cached values
      if (BlocksRequired_capacity < reserveSize) {
         BlocksRequired->reserve(newReserve,true);
         BlocksRequired_capacity = newReserve;
      }
      if (BlocksToAdd_capacity < reserveSize) {
         BlocksToAdd->reserve(newReserve,true);
         BlocksToAdd_capacity = newReserve;
      }
      if (BlocksToRemove_capacity < reserveSize) {
         BlocksToRemove->reserve(newReserve,true);
         BlocksToRemove_capacity = newReserve;
      }
      if (BlocksToMove_capacity < reserveSize) {
         BlocksToMove->reserve(newReserve,true);
         BlocksToMove_capacity = newReserve;
      }
      if (velocity_block_with_content_list_capacity < reserveSize) {
         velocity_block_with_content_list->reserve(newReserve,true);
         velocity_block_with_content_list_capacity = newReserve;
      }
      if (velocity_block_with_no_content_list_capacity < reserveSize) {
         velocity_block_with_no_content_list->reserve(newReserve,true);
         velocity_block_with_no_content_list_capacity = newReserve;
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

      // stream etc
#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif
      gpuStream_t stream = gpu_getStream();
      int nGpuBlocks;
      BlocksRequired->clear();
      BlocksToRemove->clear();
      BlocksToAdd->clear();
      BlocksToMove->clear();

      split::SplitVector<vmesh::GlobalID> *BlocksList = blockLists[thread_id];
      phiprof::Timer adjustBlocksTimer {"Adjust velocity blocks"};
      vmesh::LocalID currSize = populations[popID].vmesh->size();
      // Use cached values
      const vmesh::LocalID localContentBlocks = velocity_block_with_content_list_size;
      const vmesh::LocalID localNoContentBlocks = velocity_block_with_no_content_list_size;
      //int BlocksRequiredMapSizePower = BlocksRequiredMap->getSizePower();
      int BlocksRequiredMapSizePower = BlocksRequiredMap_sizepower; // use host-cached value
      vmesh::GlobalID* _withContentData = velocity_block_with_content_list->data(); // stored for self blocks
      vmesh::GlobalID* _withNoContentData = velocity_block_with_no_content_list->data();

      // // Neighbour and own prefetches
      // if (doPrefetches) {
      //    //phiprof::Timer prefetchTimer {"Prefetch"};
      //    populations[popID].vmesh->gpu_prefetchDevice(); // Queries active stream internally
      //    velocity_block_with_content_list->optimizeGPU(stream);
      //    velocity_block_with_no_content_list->optimizeGPU(stream);
      // }

      // Evaluate halo size
      const int addWidthV = getObjectWrapper().particleSpecies[popID].sparseBlockAddWidthV;
      const int halosize1 = (1+2*addWidthV);
      const int halosize3 = halosize1*halosize1*halosize1;
      const dim3 haloblock(halosize1,halosize1,halosize1);
      // Resize required blocks buffer to be large enough
      vmesh::LocalID nRequiredBlocksListSize = localContentBlocks*halosize3;
      if (!doDeleteEmptyBlocks) {
          nRequiredBlocksListSize += localNoContentBlocks;
      }
      // from neighbours
      const uint neighbors_count = spatial_neighbors.size();
      if (neighbors_count > 0) {
         for (std::vector<SpatialCell*>::const_iterator neighbor=spatial_neighbors.begin();
              neighbor != spatial_neighbors.end(); ++neighbor) {
            nRequiredBlocksListSize += (*neighbor)->velocity_block_with_content_list_size;
         }
      }
      BlocksList->resize(nRequiredBlocksListSize,true);
      if (BlocksRequired_capacity < localContentBlocks+4*localNoContentBlocks) {
         BlocksRequired_capacity=localContentBlocks+4*localNoContentBlocks;
         BlocksRequired->reserve(BlocksRequired_capacity,true);
      }
      nGpuBlocks = localContentBlocks > GPUBLOCKS ? GPUBLOCKS : localContentBlocks;

      //phiprof::Timer blockHaloTimer {"Block halo kernel"};
      if (localContentBlocks > 0) {
         // First add all local content blocks and their v-space halo to the list
         update_blocks_required_halo_kernel<<<nGpuBlocks, haloblock, 0, stream>>> (
            populations[popID].dev_vmesh,
            _withContentData,
            BlocksList->data(), // Start off at beginning of buffer
            localContentBlocks,
            addWidthV
            );
         CHK_ERR( gpuPeekAtLastError() );
      }
      //CHK_ERR( gpuStreamSynchronize(stream) );
      //blockHaloTimer.stop();

      //phiprof::Timer neighbourTimer {"Neighbour block memcopies"};
      // Then add neighbour blocks
      vmesh::LocalID incrementPoint = halosize3*localContentBlocks;
      if (neighbors_count > 0) {
         for (std::vector<SpatialCell*>::const_iterator neighbor=spatial_neighbors.begin();
              neighbor != spatial_neighbors.end(); ++neighbor) {
            const vmesh::LocalID nNeighBlocks = (*neighbor)->velocity_block_with_content_list_size;
            if (nNeighBlocks>0) {
               // just memcpy
               CHK_ERR( gpuMemcpyAsync(BlocksList->data()+incrementPoint,(*neighbor)->gpu_velocity_block_with_content_list_buffer,
                                       nNeighBlocks*sizeof(vmesh::LocalID), gpuMemcpyDeviceToDevice, stream) );
               incrementPoint += nNeighBlocks;
            }
         }
      }

      // Do we leave empty blocks unaltered?
      if (!doDeleteEmptyBlocks && localNoContentBlocks>0) {
         CHK_ERR( gpuMemcpyAsync(BlocksList->data()+incrementPoint,_withNoContentData,
                                 localNoContentBlocks*sizeof(vmesh::LocalID), gpuMemcpyDeviceToDevice, stream) );
         incrementPoint += localNoContentBlocks;
      }
      //CHK_ERR( gpuStreamSynchronize(stream) );
      //neighbourTimer.stop();

      //phiprof::Timer resizeTimer {"BlocksRequired hashmap resize / clear"};
      // Estimate required size based on existing blocks
      int HashmapReqSize = 2;
      if (localContentBlocks+localNoContentBlocks > 0) {
         const int HashmapReqSize2 = (localContentBlocks+localNoContentBlocks) > 0 ? (localContentBlocks+localNoContentBlocks) : 1;
         HashmapReqSize += ceil(log2(HashmapReqSize2));
      }

      if (BlocksRequiredMapSizePower >= HashmapReqSize) {
         // Map is already large enough
         BlocksRequiredMap->clear(Hashinator::targets::device,stream,false);
      } else {
         // Need larger empty map
         BlocksRequiredMap->clear(Hashinator::targets::device,stream,false);
         BlocksRequiredMap->resize(HashmapReqSize,Hashinator::targets::device, stream);
         BlocksRequiredMapSizePower = HashmapReqSize;
      }
      //resizeTimer.stop();

      // Dump all required content blocks into set/map with a fast hashinator interface.
      // GPU TODO: In fact, could do sort + unique with a list
      //phiprof::Timer blockInsertTimer {"All blocks with content"};
      //phiprof::Timer insertTimer {"BlocksRequired insert"};
      // 0.5 is target load factor
      BlocksRequiredMap->insert(BlocksList->data(),BlocksList->data(),incrementPoint,0.5,stream,false);
      CHK_ERR( gpuPeekAtLastError() );
      //CHK_ERR( gpuStreamSynchronize(stream) );
      //insertTimer.stop();

      // Ensure allocation for extraction calls is sufficient
      size_t bytesNeeded = split::tools::estimateMemoryForCompaction((size_t)std::pow(2,BlocksRequiredMapSizePower+4));
      gpu_compaction_allocate_buf_perthread(thread_id, bytesNeeded);

      //phiprof::Timer extractRequiredTimer {"BlocksRequired extract"};
      const vmesh::LocalID nBlocksRequired = BlocksRequiredMap->extractAllKeys(*BlocksRequired,compaction_buffer[thread_id],bytesNeeded,stream,false);
      //const vmesh::LocalID nBlocksRequired = BlocksRequiredMap->extractAllKeys(*BlocksRequired,stream,false);
      //extractRequiredTimer.stop();

      vmesh::LocalID nBlocksToRemove = 0;
      if (doDeleteEmptyBlocks && localNoContentBlocks>0) {
         //phiprof::Timer resizeDeleteTimer {"Blocksdelete hashmap resize / clear"};
         // Estimate required size based on existing blocks
         int HashmapDeleteReqSize = 2;
         if (localNoContentBlocks > 0) {
            const int HashmapDeleteReqSize2 = (localNoContentBlocks) > 0 ? (localNoContentBlocks) : 1;
            HashmapDeleteReqSize += ceil(log2(HashmapDeleteReqSize2));
         }
         int BlocksDeleteMapSizePower = BlocksDeleteMap->getSizePower();
         if (BlocksDeleteMapSizePower >= HashmapDeleteReqSize) {
            // Map is already large enough
            BlocksDeleteMap->clear(Hashinator::targets::device,stream,false);
         } else {
            // Need larger empty map
            BlocksDeleteMap->clear(Hashinator::targets::device,stream,false);
            BlocksDeleteMap->resize(HashmapDeleteReqSize,Hashinator::targets::device, stream);
            BlocksDeleteMapSizePower = HashmapDeleteReqSize;
         }
         //resizeDeleteTimer.stop();

         //phiprof::Timer findDeleteTimer {"BlocksToRemove all"};
         // Build set of blocks to potentially delete from localNoContentBlocks
         // 0.5 is target load factor
         BlocksDeleteMap->insert(_withNoContentData,_withNoContentData,localNoContentBlocks,0.5,stream,false);
         CHK_ERR( gpuPeekAtLastError() );
         CHK_ERR( gpuStreamSynchronize(stream) );
         // then delete all those blocks from the set which are still required
         BlocksDeleteMap->erase(BlocksRequired->data(),nBlocksRequired,stream);
         CHK_ERR( gpuStreamSynchronize(stream) );
         // And finally extract the actual list of blocks to remove
         // This extraction re-uses the same buffer as for required
         if (BlocksDeleteMap->size() > 0) {
            nBlocksToRemove = BlocksDeleteMap->extractAllKeys(*BlocksToRemove,compaction_buffer[thread_id],bytesNeeded,stream,false);
         }
         //findDeleteTimer.stop();
      }

      // Now clean the blocks required set/map of all blocks which already exist (these two calls should be merged)
      //phiprof::Timer eraseTimer {"BlocksRequired erase existing"};
      BlocksRequiredMap->erase(_withContentData,localContentBlocks,stream);
      BlocksRequiredMap->erase(_withNoContentData,localNoContentBlocks,stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      //eraseTimer.stop();
      //phiprof::Timer toAddTimer {"BlocksToAdd extract"};
      const vmesh::LocalID nBlocksToAdd = BlocksRequiredMap->extractAllKeys(*BlocksToAdd,compaction_buffer[thread_id],bytesNeeded,stream,false);
      //toAddTimer.stop();

      // If we remove more blocks than we create, we gather a list of blocks to rescue
      if (nBlocksToRemove > nBlocksToAdd) {
         update_blocks_to_move_caller_2(popID, nBlocksRequired, currSize);
      }

      // On-device adjustment calling happens in separate function as it is also called from within acceleration
      adjust_velocity_blocks_caller(popID);

      // Perform hashmap cleanup here (instead of at acceleration mid-steps)
      //phiprof::Timer cleanupTimer {"Hashinator cleanup"};
      populations[popID].vmesh->gpu_cleanHashMap(stream);
      SSYNC;
      //cleanupTimer.stop();

      #ifdef DEBUG_SPATIAL_CELL
      const size_t vmeshSize = (populations[popID].vmesh)->size();
      const size_t vbcSize = (populations[popID].blockContainer)->size();
      if (vmeshSize != vbcSize) {
         printf("ERROR: population vmesh %zu and blockcontainer %zu sizes do not match!\n",vmeshSize,vbcSize);
      }
      #endif
      #ifdef DEBUG_VLASIATOR
      // This is a bit extreme
      populations[popID].vmesh->check();
      #endif
   }

   void SpatialCell::update_blocks_to_move_caller(const uint popID) {
      phiprof::Timer updateToMoveTimer {"update_blocks_to_move"};
      // This helper calls a kernel which figures out which blocks need
      // to be rescued from the end-space of the block data.
      // To be used by acceleration in the special case that we hit v-space boundaries
      gpuStream_t stream = gpu_getStream();
      // GPUTODO: use host-caching for size
      const int nBlocksRequired = BlocksRequired->size();
      const uint nGpuBlocks = nBlocksRequired > GPUBLOCKS ? GPUBLOCKS : nBlocksRequired;
      if (BlocksToMove_capacity < nBlocksRequired) {
         BlocksToMove_capacity = nBlocksRequired*BLOCK_ALLOCATION_PADDING;
         BlocksToMove->reserve(BlocksToMove_capacity,true);
         BlocksToMove->optimizeGPU(stream);
      }
      if (nBlocksRequired>0) {
         //CHK_ERR( gpuStreamSynchronize(stream) );
         //phiprof::Timer blockMoveTimer {"blocks_to_move_kernel"};
         update_blocks_to_move_kernel<<<nGpuBlocks, GPUTHREADS, 0, stream>>> (
            populations[popID].dev_vmesh,
            BlocksRequired,
            BlocksToMove,
            nBlocksRequired
            );
         CHK_ERR( gpuPeekAtLastError() );
      }
   }
   void SpatialCell::update_blocks_to_move_caller_2(const uint popID, const vmesh::LocalID nBlocksRequired, const vmesh::LocalID nBlocksCurrent) {
      phiprof::Timer updateToMoveTimer {"update_blocks_to_move"};
      #ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
      #else
      const uint thread_id = 0;
      #endif
      // This helper calls a kernel which figures out which blocks need to be rescued 
      // from the end-space of the block data. Uses the BlocksDeleteMap (set).
      // To be used during regular block adjustment.
      const gpuStream_t stream = gpu_getStream();
      const int nEvaluate = nBlocksCurrent - nBlocksRequired;
      if (nEvaluate <= 0) {
         //blocksToMove->clear(); // Already cleared by block adjustment main function
         return;
      }
      const uint nGpuBlocks = nEvaluate > GPUBLOCKS ? GPUBLOCKS : nEvaluate;
      gpu_compaction_allocate_vec_perthread(thread_id, nEvaluate);
      if (BlocksToMove_capacity < nEvaluate) {
         BlocksToMove_capacity = nEvaluate*BLOCK_ALLOCATION_PADDING;
         BlocksToMove->reserve(BlocksToMove_capacity,true);
         BlocksToMove->optimizeGPU(stream);
      }
      BlocksToMove->resize(nEvaluate,true);
      update_blocks_to_move_kernel_2<<<nGpuBlocks, GPUTHREADS, 0, stream>>> (
         &(populations[popID].vmesh->getGrid()),
         nBlocksCurrent,
         nBlocksRequired,
         BlocksDeleteMap,
         vbwcl_gather[thread_id]->data() //pass bare pointer to buffer
         );
      CHK_ERR( gpuPeekAtLastError() );
      // Now do stream compaction on BlocksToMove, returning only valid GIDs
      size_t bytesNeeded=split::tools::estimateMemoryForCompaction((size_t)nEvaluate);
      gpu_compaction_allocate_buf_perthread(thread_id, bytesNeeded);
      auto Predicate = [] __host__ __device__ (vmesh::LocalID i ){return i != vmesh::INVALID_LOCALID; };
      CHK_ERR( gpuStreamSynchronize(stream) ); // Ensure kernel has finished
      const size_t nToMove = split::tools::copy_if(vbwcl_gather[thread_id]->data(), BlocksToMove->data(), nEvaluate,
                            Predicate, compaction_buffer[thread_id], bytesNeeded, stream);
      BlocksToMove->resize(nToMove,true,stream);
   }

   void SpatialCell::adjust_velocity_blocks_caller(const uint popID) {
      /**
          Call GPU kernel with all necessary information for creation and deletion of blocks.
          Potential optimization: take the vector lengths as input parameters
          instead of having to call the size and then prefetch back to device.
      **/
      phiprof::Timer addRemoveTimer {"GPU add and remove blocks"};

#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif
      gpuStream_t stream = gpu_getStream();
      int nGpuBlocks;

      //phiprof::Timer sizesTimer {"Block lists sizes"};
      CHK_ERR( gpuStreamSynchronize(stream) ); // To ensure all previous kernels have finished
      const vmesh::LocalID nBlocksBeforeAdjust = populations[popID].vmesh->size();
      // GPUTODO: If column extent analysis during acceleration is changed to use
      // buffers and stream compaction, we can use host-cached size values here.
      const vmesh::LocalID nToAdd = BlocksToAdd->size();
      const vmesh::LocalID nToRemove = BlocksToRemove->size();
      //const vmesh::LocalID nToMove = BlocksToMove->size(); // not used
      const vmesh::LocalID nBlocksAfterAdjust = nBlocksBeforeAdjust + nToAdd - nToRemove;
      const int nBlocksToChange = nToAdd > nToRemove ? nToAdd : nToRemove;
      //sizesTimer.stop();

      // We need to quasi-sort the BlocksToRemove list so that
      // entries at the start have LIDs < nBlocksAfterAdjust. BlocksList
      // is thread-owned and usually initialized to be very large, so
      // here we don't worry about allocation, we just resize it.
      split::SplitVector<vmesh::GlobalID> *BlocksList = blockLists[thread_id];
      BlocksList->resize(nToRemove*2,true);
      nGpuBlocks = nToRemove > GPUBLOCKS ? GPUBLOCKS : nToRemove;
      if (nGpuBlocks>0) {
         sort_blocks_remove_kernel<<<nGpuBlocks, GPUTHREADS, 0, stream>>> (
            populations[popID].dev_vmesh,
            BlocksToRemove,
            BlocksList,
            nToRemove,
            nBlocksAfterAdjust
            );
         CHK_ERR( gpuPeekAtLastError() );
         // Now do stream compaction from BlocksList back to BlocksToRemove, returning the quasi-sorted list
         size_t bytesNeeded=split::tools::estimateMemoryForCompaction((size_t)2*nToRemove);
         gpu_compaction_allocate_buf_perthread(thread_id, bytesNeeded);
         auto Predicate = [] __host__ __device__ (vmesh::GlobalID i ){return i != vmesh::INVALID_GLOBALID; };
         CHK_ERR( gpuStreamSynchronize(stream) ); // To ensure kernel has finished
         const size_t nCompacted = split::tools::copy_if(BlocksList->data(), BlocksToRemove->data(), 2*nToRemove,
                                                      Predicate, compaction_buffer[thread_id], bytesNeeded, stream);
         #ifdef DEBUG_SPATIAL_CELL
         assert((nCompacted == nToRemove) && "invalid outcome number from BlocksToRemove quasi-sorting");
         #endif
         BlocksToRemove->resize(nCompacted,true,stream);
      }
      // populations[popID].vmesh->print();

      // Grow the vmesh and block container, if necessary
      if (nBlocksAfterAdjust > nBlocksBeforeAdjust) {
         populations[popID].vmesh->setNewSize(nBlocksAfterAdjust);
         populations[popID].blockContainer->setSize(nBlocksAfterAdjust);
      }

      nGpuBlocks = nBlocksToChange > GPUBLOCKS ? GPUBLOCKS : nBlocksToChange;
      if (nGpuBlocks>0) {
         //phiprof::Timer addRemoveKernelTimer {"GPU add and remove blocks kernel"};
         CHK_ERR( gpuMemsetAsync(returnRealf[thread_id], 0, sizeof(Realf), stream) );
         const dim3 block(WID,WID,WID);
         // Third argument specifies the number of bytes in *shared memory* that is
         // dynamically allocated per block for this call in addition to the statically allocated memory.
         update_velocity_blocks_kernel<<<nGpuBlocks, block, 0, stream>>> (
            populations[popID].dev_vmesh,
            populations[popID].dev_blockContainer,
            BlocksToAdd,
            BlocksToRemove,
            BlocksToMove,
            nBlocksBeforeAdjust,
            nBlocksAfterAdjust,
            returnRealf[thread_id] // mass loss
            );
         CHK_ERR( gpuPeekAtLastError() );
         Realf host_rhoLossAdjust = 0;
         CHK_ERR( gpuMemcpyAsync(&host_rhoLossAdjust, returnRealf[thread_id], sizeof(Realf), gpuMemcpyDeviceToHost, stream) );
         CHK_ERR( gpuStreamSynchronize(stream) );
         this->populations[popID].RHOLOSSADJUST += host_rhoLossAdjust;
      }

      // Shrink the vmesh and block container, if necessary
      if (nBlocksAfterAdjust <= nBlocksBeforeAdjust) {
         populations[popID].vmesh->setNewSize(nBlocksAfterAdjust);
         populations[popID].blockContainer->setSize(nBlocksAfterAdjust);
         SSYNC;
         if (doPrefetches) {
            //phiprof::Timer prefetchTimer {"Vmesh and VBC lists prefetch dev"};
            populations[popID].vmesh->gpu_prefetchDevice();
            populations[popID].blockContainer->gpu_prefetchDevice();
            SSYNC;
         }
      }

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

      // Don't return until everything is done?
      CHK_ERR( gpuStreamSynchronize(stream) );
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
            populations[activePopID].N_blocks = populations[activePopID].blockContainer->size();

            // send velocity block list size
            displacements.push_back((uint8_t*) &(populations[activePopID].N_blocks) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE2) != 0) {
            // STAGE1 should have been done, otherwise we have problems...
            if (receiving) {
               //mpi_number_of_blocks transferred earlier
               populations[activePopID].vmesh->setNewSize(populations[activePopID].N_blocks);
            } else {
                //resize to correct size (it will avoid reallocation if it is big enough, I assume)
                populations[activePopID].N_blocks = populations[activePopID].blockContainer->size();
            }

            // send velocity block list
            displacements.push_back((uint8_t*) &(populations[activePopID].vmesh->getGrid()[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::GlobalID) * populations[activePopID].vmesh->size(true));
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
            gpuStream_t stream = gpu_getStream();
            if (receiving) {
               this->velocity_block_with_content_list->resize(this->velocity_block_with_content_list_size,true);
               this->velocity_block_with_content_list->optimizeGPU(stream);
               // Re receive velocity block content lists only for remote cells (?) so no need to
               // attach to a stream at this point.
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
      populations[popID].vmesh->setGrid();
      const vmesh::LocalID meshSize = populations[popID].vmesh->size(true);
      populations[popID].blockContainer->setSize(meshSize);
      // Set velocity block parameters:
      Real* parameters = get_block_parameters(popID);
      gpuStream_t stream = gpu_getStream();
      const uint nGpuBlocks = (meshSize/GPUTHREADS) > GPUBLOCKS ? GPUBLOCKS : std::ceil((Real)meshSize/(Real)GPUTHREADS);
      CHK_ERR( gpuStreamSynchronize(stream) );
      if (nGpuBlocks>0) {
         update_blockparameters_kernel<<<nGpuBlocks, GPUTHREADS, 0, stream>>> (
            populations[popID].dev_vmesh,
            populations[popID].dev_blockContainer,
            parameters,
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
         if (populations[p].blockContainer->capacity() > amount )
            if (populations[p].blockContainer->recapacitate(amount) == false) success = false;

      }
      return success;
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
#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif

      phiprof::Timer updateListsTimer {"GPU update spatial cell block lists"};
      gpuStream_t stream = gpu_getStream();
      velocity_block_with_content_list_size = 0;
      velocity_block_with_no_content_list_size = 0;
      const vmesh::LocalID currSize = populations[popID].vmesh->size();

      // Immediate return if no blocks to process
      if (currSize == 0) {
         velocity_block_with_content_list->clear();
         velocity_block_with_no_content_list->clear();
         return;
      }

      // Ensure allocation of gathering vectors and buffers
      size_t bytesNeeded=split::tools::estimateMemoryForCompaction((size_t)currSize);
      gpu_compaction_allocate_buf_perthread(thread_id, bytesNeeded);
      gpu_compaction_allocate_vec_perthread(thread_id, currSize);

      const vmesh::LocalID currCapacity = velocity_block_with_content_list_capacity; //host-cached
      if (currCapacity < currSize) {
         const uint reserveSize = currSize * BLOCK_ALLOCATION_FACTOR;
         velocity_block_with_content_list->reserve(reserveSize,true, stream);
         velocity_block_with_no_content_list->reserve(reserveSize,true, stream);
         velocity_block_with_content_list_capacity=reserveSize;
         velocity_block_with_no_content_list_capacity=reserveSize;
      }
      // Set gathering vectors to correct size
      vbwcl_gather[thread_id]->resize(currSize,true,stream);
      vbwncl_gather[thread_id]->resize(currSize,true,stream);
      velocity_block_with_content_list->resize(currSize,true,stream);
      velocity_block_with_no_content_list->resize(currSize,true,stream);

      phiprof::Timer kernelTimer {"GPU block lists kernel"};
      const Real velocity_block_min_value = getVelocityBlockMinValue(popID);
      dim3 block(WID,WID,WID);
      // Third argument specifies the number of bytes in *shared memory* that is
      // dynamically allocated per block for this call in addition to the statically allocated memory.
      update_velocity_block_content_lists_kernel<<<GPUBLOCKS, block, WID3*sizeof(bool), stream>>> (
         populations[popID].dev_vmesh,
         currSize,
         populations[popID].dev_blockContainer,
         vbwcl_gather[thread_id]->data(), // Now pass bare temporary gathering buffers
         vbwncl_gather[thread_id]->data(),
         velocity_block_min_value
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuStreamSynchronize(stream) ); // This sync is required!
      kernelTimer.stop();

      phiprof::Timer compactionTimer {"GPU block lists streamcompaction"};
      // Now do stream compaction on those two vectors, returning only valid GIDs
      // into the actual vectors.
      auto Predicate = [] __host__ __device__ (vmesh::GlobalID i ){return i != vmesh::INVALID_GLOBALID; };
      velocity_block_with_content_list_size = split::tools::copy_if(vbwcl_gather[thread_id]->data(), velocity_block_with_content_list->data(), currSize,
                            Predicate, compaction_buffer[thread_id], bytesNeeded, stream);
      velocity_block_with_content_list->resize(velocity_block_with_content_list_size,true,stream);
      velocity_block_with_no_content_list_size = split::tools::copy_if(vbwncl_gather[thread_id]->data(), velocity_block_with_no_content_list->data(), currSize,
                            Predicate, compaction_buffer[thread_id], bytesNeeded, stream);
      velocity_block_with_no_content_list->resize(velocity_block_with_no_content_list_size,true,stream);
      //CHK_ERR( gpuStreamSynchronize(stream) );
      compactionTimer.stop();

      // Note: Content list is not uploaded to device-only buffer here, but rather
      // in grid.cpp adjustVelocityBlocks()
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

   void SpatialCell::printMeshSizes() {
      cerr << "SC::printMeshSizes:" << endl;
      for (size_t p=0; p<populations.size(); ++p) {
         cerr << "\t pop " << p << " " << populations[p].vmesh->size(true) << ' ' << populations[p].blockContainer->size(true) << endl;
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
