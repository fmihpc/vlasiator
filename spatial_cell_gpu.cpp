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

#ifndef NDEBUG
   #define DEBUG_SPATIAL_CELL
#endif

using namespace std;

/** GPU kernel for identifying which blocks have relevant content */
__global__ void __launch_bounds__(WID3,4) update_velocity_block_content_lists_kernel (
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   split::SplitVector<vmesh::GlobalID>* velocity_block_with_content_list,
   split::SplitVector<vmesh::GlobalID>* velocity_block_with_no_content_list,
   Realf velocity_block_min_value
   ) {

   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;
   __shared__ bool has_content[WID3];
   const uint nBlocks = vmesh->size();
   for (uint blockLID=blocki; blockLID<nBlocks; blockLID += gpuBlocks) {
      vmesh::GlobalID blockGID = vmesh->getGlobalID(blockLID);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == vmesh->invalidGlobalID()) {
         continue;
      }
      if (blockLID == vmesh->invalidLocalID()) {
         continue;
      }
      #endif
      Realf* avgs = blockContainer->getData(blockLID);
      has_content[ti] = avgs[ti] >= velocity_block_min_value ? true : false;
      __syncthreads();
      // Implemented just a simple non-optimized thread OR
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (ti < s) {
            has_content[ti] = has_content[ti] || has_content[ti + s];
         }
         __syncthreads();
      }
      if (ti==0) {
         if (has_content[0]) {
            velocity_block_with_content_list->device_push_back(blockGID);
         } else {
            velocity_block_with_no_content_list->device_push_back(blockGID);
         }
      }
      __syncthreads();
   }
}

/** Gpu Kernel to quickly gather blocks and their v-space halo */
__global__ void __launch_bounds__(GPUTHREADS,4) update_blocks_required_halo_kernel (
   vmesh::VelocityMesh *vmesh,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* BlocksRequiredMap,
   split::SplitVector<vmesh::GlobalID> *velocity_block_with_content_list,
   split::SplitVector<vmesh::GlobalID> *BlocksHalo,
   const int addWidthV,
   // The following 4 vectors are passed just to be able to clear them on-device
   split::SplitVector<vmesh::GlobalID>* BlocksRequired,
   split::SplitVector<vmesh::GlobalID>* BlocksToRemove,
   split::SplitVector<vmesh::GlobalID>* BlocksToAdd,
   split::SplitVector<vmesh::GlobalID>* BlocksToMove
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const vmesh::LocalID ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   const unsigned long localContentBlocks = velocity_block_with_content_list->size();

   vmesh::GlobalID nGIDs[27] = {0};

   for (vmesh::LocalID index=blocki*warpSize; index<localContentBlocks; index += gpuBlocks*warpSize) {
      if (index+ti < localContentBlocks) {
         vmesh::GlobalID GID = velocity_block_with_content_list->at(index+ti);
         //BlocksRequiredMap->set_element(GID,GID); // Already added via Hashinator interface
         vmesh::LocalID ind0,ind1,ind2;
         vmesh->getIndices(GID,ind0,ind1,ind2);

         // New algorithm: attempts to minimize thread divergence and use of atomics.
         // Gather existence of all 26 neighbours into 32-bit bitmask
         vmesh::LocalID exist = 0;
         int localCounter=0;
         for (int offset_vx=-addWidthV; offset_vx<=addWidthV; offset_vx++) {
            for (int offset_vy=-addWidthV; offset_vy<=addWidthV; offset_vy++) {
               for (int offset_vz=-addWidthV; offset_vz<=addWidthV; offset_vz++) {
                  const int nind0 = ind0 + offset_vx;
                  const int nind1 = ind1 + offset_vy;
                  const int nind2 = ind2 + offset_vz;
                  const vmesh::GlobalID nGID
                     = vmesh->getGlobalID(nind0,nind1,nind2);
                  if (nGID != vmesh->invalidGlobalID()) {
                     // NEW:
                     nGIDs[localCounter] = nGID;
                     if (BlocksRequiredMap->device_count(nGID) != 0) {
                        exist = (exist | (1 << localCounter));
                     }
                     // OLD:
                     // if (BlocksRequiredMap->device_count(nGID) == 0) {
                     //    //BlocksHalo->device_push_back(nGID);
                     //    BlocksRequiredMap->set_element(nGID,nGID);
                     // }
                  }
                  localCounter++;
               } // for vz
            } // for vy
         } // for vx

         // Quick check: if all halo neighbours already exist, skip to next loop.
         if (exist == 134217727) {
            continue;
         }
         /**
            Define order of responsibility. Blocks are added by:
            first their face neighbours in x
            then their face neighbours in y
            then their face neighbours in z
            then their edge neighbours
            then their corner neighbours

            Bitmask: (cell 13 is self)

            Z
            |    8      17      26
            |  5      14      23
            |2   7  11  16  20  25
            |  4      13      22
            |1   6  10  15  19  24
            |  3      12      21
            |0      9       18
            L__________________ X

         // First the 8 corners: add only if corner cell has no other existing neighours (that we know of)

         // Then the 12 blocks which are between two corners
         // Add only if the four face neighbours we know of are not there

         // Then face neighbours
         // Add z-neighbours only if we don't know of any x or y face neighbours
         // Add y-neighbours only if we don't know of any x face neighbours
         // Add all x-face-neighbours

Start python code:
def val(x,y,z):
    # input is distances as -1,0,+1
    # returns the index
    return (z+1) + (y+1)*3 + (x+1)*9

def coord(inde):
    z = (inde % 3)
    y = ((inde//3) % 3)
    x = ((inde//9) % 3)
    return (x-1,y-1,z-1)

def dims2(inde):
    (i0,j0,k0) = coord(inde)
    listi = []
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                di = abs(i-i0)
                dj = abs(j-j0)
                dk = abs(k-k0)
                if (di<2) and (dj<2) and (dk<2) and ((di+dj+dk)<3):
                    listi.append(val(i,j,k))
    return listi

def neigh_x(inde):
    # looking at x-directional neighbour
    # only returns self
    listi = [inde]
    return listi

def neigh_y(inde):
    # looking at y-directional neighbour
    # returns the target and any x-directional neighbours of it
    (i0,j0,k0) = coord(inde)
    listi = []
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                di = abs(i-i0)
                dj = abs(j-j0)
                dk = abs(k-k0)
                if (di<2) and (dj<1) and (dk<1):
                    listi.append(val(i,j,k))
    return listi

def neigh_z(inde):
    # looking at z-directional neighbour
    # returns the target and any x-directional or y-directional neighbours of it
    (i0,j0,k0) = coord(inde)
    listi = []
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                di = abs(i-i0)
                dj = abs(j-j0)
                dk = abs(k-k0)
                if (di<2) and (dj<2) and (dk<1) and ((di+dj+dk)<2):
                    listi.append(val(i,j,k))
    return listi

def dims1(inde):
    (i0,j0,k0) = coord(inde)
    listi = []
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                di = abs(i-i0)
                dj = abs(j-j0)
                dk = abs(k-k0)
                if (di<2) and (dj<2) and (dk<2) and ((di+dj+dk)<2):
                    listi.append(val(i,j,k))
    return listi

def findmask(start,indices):
    val  = 0
    strmoi=""
    for x in indices:
        val |= (1<<x)
        strmoi = strmoi+str(x)+" "
    #print(start,indices,val)
    print("         // "+str(start)+": "+strmoi+" -> "+str(val))
    print("         if ( ( (exist & (vmesh::LocalID)"+str(val)+") == 0) && (nGIDs["+str(start)+"] != 0)) {")
    print("            BlocksHalo->device_push_back(nGIDs["+str(start)+"]);")
    print("         }")

def findall():
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                di = abs(i)
                dj = abs(j)
                dk = abs(k)
                if (di+dj+dk==3):
                    # corner
                    thiss = val(i,j,k)
                    findmask(thiss,dims2(thiss))
                elif (di+dj+dk==2):
                    # mid-edge
                    thiss = val(i,j,k)
                    findmask(thiss,dims1(thiss))
                elif (di+dj+dk==0):
                    # own, do nothing
                    pass
                elif ((di==1) and (dj+dk==0)):
                    # x-directional neighbour
                    thiss = val(i,j,k)
                    findmask(thiss,neigh_x(thiss))
                elif ((dj==1) and (di+dk==0)):
                    # y-directional neighbour
                    thiss = val(i,j,k)
                    findmask(thiss,neigh_y(thiss))
                elif ((dk==1) and (di+dj==0)):
                    # z-directional neighbour
                    thiss = val(i,j,k)
                    findmask(thiss,neigh_z(thiss))
                else:
                    print("error passed through!")
         **/

         // 0: 0 1 3 4 9 10 12  -> 5659
         if ( ( (exist & (vmesh::LocalID)5659) == 0) && (nGIDs[0] != 0)) {
            BlocksHalo->device_push_back(nGIDs[0]);
         }
         // 1: 0 1 2 4 10  -> 1047
         if ( ( (exist & (vmesh::LocalID)1047) == 0) && (nGIDs[1] != 0)) {
            BlocksHalo->device_push_back(nGIDs[1]);
         }
         // 2: 1 2 4 5 10 11 14  -> 19510
         if ( ( (exist & (vmesh::LocalID)19510) == 0) && (nGIDs[2] != 0)) {
            BlocksHalo->device_push_back(nGIDs[2]);
         }
         // 3: 0 3 4 6 12  -> 4185
         if ( ( (exist & (vmesh::LocalID)4185) == 0) && (nGIDs[3] != 0)) {
            BlocksHalo->device_push_back(nGIDs[3]);
         }
         // 4: 4  -> 16
         if ( ( (exist & (vmesh::LocalID)16) == 0) && (nGIDs[4] != 0)) {
            BlocksHalo->device_push_back(nGIDs[4]);
         }
         // 5: 2 4 5 8 14  -> 16692
         if ( ( (exist & (vmesh::LocalID)16692) == 0) && (nGIDs[5] != 0)) {
            BlocksHalo->device_push_back(nGIDs[5]);
         }
         // 6: 3 4 6 7 12 15 16  -> 102616
         if ( ( (exist & (vmesh::LocalID)102616) == 0) && (nGIDs[6] != 0)) {
            BlocksHalo->device_push_back(nGIDs[6]);
         }
         // 7: 4 6 7 8 16  -> 66000
         if ( ( (exist & (vmesh::LocalID)66000) == 0) && (nGIDs[7] != 0)) {
            BlocksHalo->device_push_back(nGIDs[7]);
         }
         // 8: 4 5 7 8 14 16 17  -> 213424
         if ( ( (exist & (vmesh::LocalID)213424) == 0) && (nGIDs[8] != 0)) {
            BlocksHalo->device_push_back(nGIDs[8]);
         }
         // 9: 0 9 10 12 18  -> 267777
         if ( ( (exist & (vmesh::LocalID)267777) == 0) && (nGIDs[9] != 0)) {
            BlocksHalo->device_push_back(nGIDs[9]);
         }
         // 10: 1 10 19  -> 525314
         if ( ( (exist & (vmesh::LocalID)525314) == 0) && (nGIDs[10] != 0)) {
            BlocksHalo->device_push_back(nGIDs[10]);
         }
         // 11: 2 10 11 14 20  -> 1068036
         if ( ( (exist & (vmesh::LocalID)1068036) == 0) && (nGIDs[11] != 0)) {
            BlocksHalo->device_push_back(nGIDs[11]);
         }
         // 12: 3 9 12 15 21  -> 2134536
         if ( ( (exist & (vmesh::LocalID)2134536) == 0) && (nGIDs[12] != 0)) {
            BlocksHalo->device_push_back(nGIDs[12]);
         }
         // 14: 5 11 14 17 23  -> 8538144
         if ( ( (exist & (vmesh::LocalID)8538144) == 0) && (nGIDs[14] != 0)) {
            BlocksHalo->device_push_back(nGIDs[14]);
         }
         // 15: 6 12 15 16 24  -> 16879680
         if ( ( (exist & (vmesh::LocalID)16879680) == 0) && (nGIDs[15] != 0)) {
            BlocksHalo->device_push_back(nGIDs[15]);
         }
         // 16: 7 16 25  -> 33620096
         if ( ( (exist & (vmesh::LocalID)33620096) == 0) && (nGIDs[16] != 0)) {
            BlocksHalo->device_push_back(nGIDs[16]);
         }
         // 17: 8 14 16 17 26  -> 67322112
         if ( ( (exist & (vmesh::LocalID)67322112) == 0) && (nGIDs[17] != 0)) {
            BlocksHalo->device_push_back(nGIDs[17]);
         }
         // 18: 9 10 12 18 19 21 22  -> 7083520
         if ( ( (exist & (vmesh::LocalID)7083520) == 0) && (nGIDs[18] != 0)) {
            BlocksHalo->device_push_back(nGIDs[18]);
         }
         // 19: 10 18 19 20 22  -> 6030336
         if ( ( (exist & (vmesh::LocalID)6030336) == 0) && (nGIDs[19] != 0)) {
            BlocksHalo->device_push_back(nGIDs[19]);
         }
         // 20: 10 11 14 19 20 22 23  -> 14175232
         if ( ( (exist & (vmesh::LocalID)14175232) == 0) && (nGIDs[20] != 0)) {
            BlocksHalo->device_push_back(nGIDs[20]);
         }
         // 21: 12 18 21 22 24  -> 23334912
         if ( ( (exist & (vmesh::LocalID)23334912) == 0) && (nGIDs[21] != 0)) {
            BlocksHalo->device_push_back(nGIDs[21]);
         }
         // 22: 22  -> 4194304
         if ( ( (exist & (vmesh::LocalID)4194304) == 0) && (nGIDs[22] != 0)) {
            BlocksHalo->device_push_back(nGIDs[22]);
         }
         // 23: 14 20 22 23 26  -> 80756736
         if ( ( (exist & (vmesh::LocalID)80756736) == 0) && (nGIDs[23] != 0)) {
            BlocksHalo->device_push_back(nGIDs[23]);
         }
         // 24: 12 15 16 21 22 24 25  -> 56725504
         if ( ( (exist & (vmesh::LocalID)56725504) == 0) && (nGIDs[24] != 0)) {
            BlocksHalo->device_push_back(nGIDs[24]);
         }
         // 25: 16 22 24 25 26  -> 121700352
         if ( ( (exist & (vmesh::LocalID)121700352) == 0) && (nGIDs[25] != 0)) {
            BlocksHalo->device_push_back(nGIDs[25]);
         }
         // 26: 14 16 17 22 23 25 26  -> 113459200
         if ( ( (exist & (vmesh::LocalID)113459200) == 0) && (nGIDs[26] != 0)) {
            BlocksHalo->device_push_back(nGIDs[26]);
         }
      } // if index
   } // for blocks

   // One thread also sets the sizes of these vectors to zero
   if (blockIdx.x == blockIdx.y == blockIdx.z == threadIdx.x == threadIdx.y == threadIdx.z == 0) {
      BlocksRequired->clear();
      BlocksToRemove->clear();
      BlocksToAdd->clear();
      BlocksToMove->clear();
   }
}

/** Gpu Kernel to quickly add blocks which have spatial neighbours */
__global__ void __launch_bounds__(GPUTHREADS,4) update_neighbours_have_content_kernel (
   vmesh::VelocityMesh *vmesh,
   split::SplitVector<vmesh::GlobalID> *BlocksHalo,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* BlocksRequiredMap,
   vmesh::GlobalID *neighbor_velocity_block_with_content_list,
   const unsigned long neighborContentBlocks
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const vmesh::LocalID ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   for (vmesh::LocalID index=blocki*warpSize; index<neighborContentBlocks; index += gpuBlocks*warpSize) {
      if (index+ti < neighborContentBlocks) {
         vmesh::GlobalID GID = neighbor_velocity_block_with_content_list[index+ti];
         if (BlocksRequiredMap->device_count(GID) == 0) {
            //BlocksRequiredMap->set_element(GID,GID);
            BlocksHalo->device_push_back(GID);
         }
      }
   }
}

/** GPU kernel for selecting only non-existing blocks for addition
    This kernel may be non-optimized in itself, but use of it gets rid
    of the need of vmesh prefetching back and forth.
 */
__global__ void __launch_bounds__(GPUTHREADS,4) update_blocks_to_add_kernel (
   vmesh::VelocityMesh *vmesh,
   split::SplitVector<vmesh::GlobalID>* BlocksRequired,
   split::SplitVector<vmesh::GlobalID>* BlocksToAdd,
   split::SplitVector<vmesh::GlobalID>* BlocksToMove,
   const uint nBlocksRequired
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (uint index=blocki*warpSize; index<nBlocksRequired; index += gpuBlocks*warpSize) {
      if (index+ti < nBlocksRequired) {
         const vmesh::GlobalID GIDreq = BlocksRequired->at(index+ti);
         const vmesh::LocalID LIDreq = vmesh->getLocalID(GIDreq);
         if ( LIDreq == vmesh->invalidLocalID() ) {
            BlocksToAdd->device_push_back(GIDreq);
         } else if (LIDreq >= nBlocksRequired) {
            BlocksToMove->device_push_back(GIDreq);
         }
      }
   }
}

/** GPU kernel for identifying blocks for deletion.
    This kernel may be non-optimized in itself, but use of it gets rid
    of the need of vmesh prefetching back and forth.
 */
__global__ void __launch_bounds__(GPUTHREADS,4) update_blocks_to_remove_kernel (
   split::SplitVector<vmesh::GlobalID>* velocity_block_with_no_content_list,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* BlocksRequiredMap,
   split::SplitVector<vmesh::GlobalID>* BlocksToRemove,
   const uint localNoContentBlocks
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int warpSize = blockDim.x*blockDim.y*blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (uint index=blocki*warpSize; index<localNoContentBlocks; index += gpuBlocks*warpSize) {
      if (index+ti < localNoContentBlocks) {
         const vmesh::GlobalID GIDcandidate = velocity_block_with_no_content_list->at(index+ti);
         if (BlocksRequiredMap->device_count(GIDcandidate) == 0) {
            BlocksToRemove->device_push_back(GIDcandidate);
         }
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
   vmesh::LocalID *VectorIndex,
   Realf* gpu_rhoLossAdjust
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;
   vmesh::LocalID *addVectorIndex = VectorIndex;
   vmesh::LocalID *moveVectorIndex = VectorIndex+1;

   const vmesh::LocalID nToAdd = BlocksToAdd->size();
   const vmesh::LocalID nToRemove = BlocksToRemove->size();
   const vmesh::LocalID nToMove = BlocksToMove->size();
   const vmesh::LocalID nToCreate = nToAdd > nToRemove ? (nToAdd-nToRemove) : 0;
   Realf local_rhoLoss = 0;
   // For tracking mass-loss
   __shared__ Realf massloss[WID3];
   __shared__ vmesh::LocalID moveIndex;
   __shared__ vmesh::LocalID addIndex;

   for (vmesh::LocalID m=blocki; m<nToRemove; m += gpuBlocks) {
      // Go through all blocks which are to be removed.
      // If there is a corresponding block to be added, place that in its stead.
      // Otherwise, take the corresponding block from the moved list instead.
      const vmesh::GlobalID rmGID = BlocksToRemove->at(m);
      const vmesh::LocalID rmLID = vmesh->getLocalID(rmGID);

      #ifdef DEBUG_SPATIAL_CELL
      if (rmGID == vmesh->invalidGlobalID()) {
         continue;
      }
      if (rmLID == vmesh->invalidLocalID()) {
         continue;
      }
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

      // Figure out which GID to put here instead
      if (rmLID >= nBlocksAfterAdjust) {
         // Delete without replacing
         if (ti==0) {
            vmesh->deleteBlock(rmGID,rmLID);
         }
         __syncthreads();
         continue;
      }
      if (ti==0) moveIndex = atomicAdd(moveVectorIndex,1);
      __syncthreads();
      if (moveIndex<nToMove) {
         // Move from latter part of vmesh
         const vmesh::GlobalID replaceGID = BlocksToMove->at(moveIndex);
         const vmesh::LocalID replaceLID = vmesh->getLocalID(replaceGID);
         Realf* repl_avgs = blockContainer->getData(replaceLID);
         Real*  repl_block_parameters = blockContainer->getParameters(replaceLID);
         rm_avgs[ti] = repl_avgs[ti];
         if (ti < BlockParams::N_VELOCITY_BLOCK_PARAMS) {
            rm_block_parameters[ti] = repl_block_parameters[ti];
         }
         __syncthreads();
         if (ti==0) {
            // Remove hashmap entry for removed block, add instead created block
            vmesh->replaceBlock(rmGID,rmLID,replaceGID);
         }
         __syncthreads();
         #ifdef DEBUG_SPATIAL_CELL
         if (vmesh->getGlobalID(rmLID) == vmesh->invalidGlobalID()) {
            continue;
         }
         if (vmesh->getLocalID(replaceGID) == vmesh->invalidLocalID()) {
            continue;
         }
         #endif
         continue;
      }
      if (ti==0) addIndex = atomicAdd(addVectorIndex,1);
      __syncthreads();
      if (addIndex<nToAdd) {
         // New GID
         const vmesh::GlobalID addGID = BlocksToAdd->at(addIndex);
         #ifdef DEBUG_SPATIAL_CELL
         if (addGID == vmesh->invalidGlobalID()) {
            __syncthreads();
            continue;
         }
         #endif
         rm_avgs[ti] = 0;
         if (ti==0) {
            // Write in block parameters
            vmesh->getCellSize(addGID,&(rm_block_parameters[BlockParams::DVX]));
            vmesh->getBlockCoordinates(addGID,&(rm_block_parameters[BlockParams::VXCRD]));
            vmesh->replaceBlock(rmGID,rmLID,addGID);
         }
         __syncthreads();
         #ifdef DEBUG_SPATIAL_CELL
         if (vmesh->getGlobalID(rmLID) == vmesh->invalidGlobalID()) {
            continue;
         }
         if (vmesh->getLocalID(addGID) == vmesh->invalidLocalID()) {
            continue;
         }
         #endif
         continue;
      }
      #ifdef DEBUG_SPATIAL_CELL
      if (ti==0) {
         printf("Error! Fall through in update_velocity_blocks_kernel! \n");
      }
      #endif
      __syncthreads();
   }
   // Now, if we need to expand the size of the vmesh, let's add blocks.
   // For thread-safety,this assumes that the localToGlobalMap is already of sufficient size, as should be
   // the block_data and block_parameters vectors.
   for (vmesh::LocalID m=blocki; m<nToCreate; m += gpuBlocks) {
      // Debug check: if we are adding elements, then nToMove should be zero
      // We have already used up nToRemove entries from the addition vector.
      const vmesh::GlobalID addGID = BlocksToAdd->at(nToRemove+m);
      // We need to add the data of addGID to a new LID:
      const vmesh::LocalID addLID = nBlocksBeforeAdjust + m;
      Realf* add_avgs = blockContainer->getData(addLID);
      Real* add_block_parameters = blockContainer->getParameters(addLID);
      // Zero out blockdata
      add_avgs[ti] = 0;
      if (ti==0) {
         // Write in block parameters
         vmesh->getCellSize(addGID,&(add_block_parameters[BlockParams::DVX]));
         vmesh->getBlockCoordinates(addGID,&(add_block_parameters[BlockParams::VXCRD]));
         vmesh->placeBlock(addGID,addLID);
      }
      __syncthreads();
      #ifdef DEBUG_SPATIAL_CELL
      if (vmesh->getGlobalID(addLID) == vmesh->invalidGlobalID()) {
         continue;
      }
      if (vmesh->getLocalID(addGID) == vmesh->invalidLocalID()) {
         continue;
      }
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
         populations[popID].velocityBlockMinValue = spec.sparseMinValue;
         populations[popID].N_blocks = 0;
      }

      // SplitVectors via pointers for unified memory
      velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(1);
      velocity_block_with_no_content_list = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksHalo = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksRequired = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToAdd = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToRemove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToMove = new split::SplitVector<vmesh::GlobalID>(1);
      velocity_block_with_content_list->clear();
      velocity_block_with_no_content_list->clear();
      BlocksHalo->clear();
      BlocksRequired->clear();
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();
      attachedStream=0;
      BlocksRequiredMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);

      gpuMallocHost((void **) &info_vbwcl, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_vbwncl, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_toRemove, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_toAdd, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_toMove, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_Required, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_Halo, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_brm, sizeof(Hashinator::MapInfo));
   }

   SpatialCell::~SpatialCell() {
      delete velocity_block_with_content_list;
      delete velocity_block_with_no_content_list;
      delete BlocksHalo;
      delete BlocksRequired;
      delete BlocksToAdd;
      delete BlocksToRemove;
      delete BlocksToMove;
      delete BlocksRequiredMap;
      gpuFreeHost(info_vbwcl);
      gpuFreeHost(info_vbwncl);
      gpuFreeHost(info_toRemove);
      gpuFreeHost(info_toAdd);
      gpuFreeHost(info_toMove);
      gpuFreeHost(info_Required);
      gpuFreeHost(info_Halo);
      gpuFreeHost(info_brm);
   }

   SpatialCell::SpatialCell(const SpatialCell& other) {
      velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(*(other.velocity_block_with_content_list));
      velocity_block_with_no_content_list = new split::SplitVector<vmesh::GlobalID>(*(other.velocity_block_with_no_content_list));
      BlocksHalo = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksRequired = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToAdd = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToRemove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToMove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksHalo->clear();
      BlocksRequired->clear();
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();
      velocity_block_with_content_list->clear();
      velocity_block_with_no_content_list->clear();

      BlocksRequiredMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);

      // Make space reservation guesses based on popID 0
      const uint reserveSize = other.populations[0].vmesh->size()*BLOCK_ALLOCATION_PADDING;
      BlocksHalo->reserve(reserveSize,true);
      BlocksRequired->reserve(reserveSize,true);
      BlocksToAdd->reserve(reserveSize,true);
      BlocksToRemove->reserve(reserveSize,true);
      BlocksToMove->reserve(reserveSize,true);
      velocity_block_with_content_list->reserve(reserveSize,true);
      velocity_block_with_no_content_list->reserve(reserveSize,true);

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
         neighbor_block_data[i] = other.neighbor_block_data[i];
         neighbor_number_of_blocks[i] = other.neighbor_number_of_blocks[i];
      }

      if (other.face_neighbor_ranks.size()>0) {
         face_neighbor_ranks = std::map<int,std::set<int>>(other.face_neighbor_ranks);
      }
      if (other.populations.size()>0) {
         populations = std::vector<spatial_cell::Population>(other.populations);
      }
      attachedStream=0;
      gpuMallocHost((void **) &info_vbwcl, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_vbwncl, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_toRemove, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_toAdd, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_toMove, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_Required, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_Halo, sizeof(split::SplitInfo));
      gpuMallocHost((void **) &info_brm, sizeof(Hashinator::MapInfo));
   }
   const SpatialCell& SpatialCell::operator=(const SpatialCell& other) {
      const uint reserveSize = (other.BlocksRequired)->capacity();
      BlocksHalo->clear();
      BlocksRequired->clear();
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();
      velocity_block_with_content_list->clear();
      velocity_block_with_no_content_list->clear();
      delete BlocksRequiredMap;

      BlocksHalo->reserve(reserveSize,true);
      BlocksRequired->reserve(reserveSize,true);
      BlocksToAdd->reserve(reserveSize,true);
      BlocksToRemove->reserve(reserveSize,true);
      BlocksToMove->reserve(reserveSize,true);
      velocity_block_with_content_list->reserve(reserveSize,true);
      velocity_block_with_no_content_list->reserve(reserveSize,true);

      const vmesh::LocalID HashmapReqSize = ceil(log2(reserveSize)) +2;
      BlocksRequiredMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(HashmapReqSize);

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
         neighbor_block_data[i] = other.neighbor_block_data[i];
         neighbor_number_of_blocks[i] = other.neighbor_number_of_blocks[i];
      }

      face_neighbor_ranks = std::map<int,std::set<int>>(other.face_neighbor_ranks);
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
      // CHK_ERR( gpuMemAdvise(ptr, count, advise, deviceID) );
      // CHK_ERR( gpuMemAdvise(velocity_block_with_content_list, sizeof(velocity_block_with_content_list),gpuMemAdviseSetPreferredLocation, gpu_getDevice()) );
      // gpu_getDevice()
      int device = gpu_getDevice();
      gpuStream_t stream = gpu_getStream();
      BlocksHalo->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksRequired->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksToAdd->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksToRemove->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksToMove->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      velocity_block_with_content_list->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      velocity_block_with_no_content_list->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      BlocksRequiredMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);

      BlocksHalo->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksRequired->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksToAdd->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksToRemove->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksToMove->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      velocity_block_with_content_list->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      velocity_block_with_no_content_list->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      BlocksRequiredMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);

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
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksHalo, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToRemove, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToAdd, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToMove, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksRequired, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksRequiredMap, 0,gpuMemAttachSingle) );
      // Loop over populations
      for (size_t p=0; p<populations.size(); ++p) {
         populations[p].blockContainer->gpu_attachToStream(attachedStream);
         populations[p].vmesh->gpu_attachToStream(attachedStream);
      }
      // Also call attach functions on all splitvectors and hashmaps
      velocity_block_with_content_list->streamAttach(attachedStream);
      velocity_block_with_no_content_list->streamAttach(attachedStream);
      BlocksHalo->streamAttach(attachedStream);
      BlocksToRemove->streamAttach(attachedStream);
      BlocksToAdd->streamAttach(attachedStream);
      BlocksToMove->streamAttach(attachedStream);
      BlocksRequired->streamAttach(attachedStream);
      BlocksRequiredMap->streamAttach(attachedStream);
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
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksHalo, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToRemove, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToAdd, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksToMove, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksRequired, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksRequiredMap, 0,gpuMemAttachGlobal) );
      // Loop over populations
      for (size_t p=0; p<populations.size(); ++p) {
         populations[p].blockContainer->gpu_detachFromStream();
         populations[p].vmesh->gpu_detachFromStream();
      }
      // Also call detach functions on all splitvectors and hashmaps
      velocity_block_with_content_list->streamAttach(0,gpuMemAttachGlobal);
      velocity_block_with_no_content_list->streamAttach(0,gpuMemAttachGlobal);
      BlocksHalo->streamAttach(0,gpuMemAttachGlobal);
      BlocksToRemove->streamAttach(0,gpuMemAttachGlobal);
      BlocksToAdd->streamAttach(0,gpuMemAttachGlobal);
      BlocksToMove->streamAttach(0,gpuMemAttachGlobal);
      BlocksRequired->streamAttach(0,gpuMemAttachGlobal);
      BlocksRequiredMap->streamAttach(0,gpuMemAttachGlobal);
      return;
   }

   /** Sends the contents of velocity_block_with_content_list into a device buffer so that it can be accessed
       from several streams at once.
    */
   void SpatialCell::gpu_uploadContentLists() {
      phiprof::start("Upload local content lists");
      gpuStream_t stream = gpu_getStream();
      velocity_block_with_content_list->copyMetadata(info_vbwcl,stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      velocity_block_with_content_list_size = info_vbwcl->size;
      if (velocity_block_with_content_list_size==0) {
         return;
      }
      CHK_ERR( gpuMallocAsync((void**)&gpu_velocity_block_with_content_list_buffer, velocity_block_with_content_list_size*sizeof(vmesh::LocalID), stream) );
      CHK_ERR( gpuMemcpyAsync(gpu_velocity_block_with_content_list_buffer, velocity_block_with_content_list->data(), velocity_block_with_content_list_size*sizeof(vmesh::LocalID), gpuMemcpyDeviceToDevice, stream) );
      SSYNC;
      phiprof::stop("Upload local content lists");
   }
   /** Clears the device buffer for velocity_block_with_content_list
    */
   void SpatialCell::gpu_clearContentLists() {
      gpuStream_t stream = gpu_getStream();
      if (velocity_block_with_content_list_size==0) {
         return;
      }
      CHK_ERR( gpuFreeAsync(gpu_velocity_block_with_content_list_buffer, stream) );
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
                                            const uint popID,bool doDeleteEmptyBlocks) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      // stream etc
#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif
      gpuStream_t stream = gpu_getStream();
      int nGpuBlocks;

      phiprof::start("Adjust velocity blocks");
      velocity_block_with_content_list->copyMetadata(info_vbwcl,stream);
      velocity_block_with_no_content_list->copyMetadata(info_vbwncl,stream);
      BlocksRequired->copyMetadata(info_Required,stream);
      BlocksHalo->copyMetadata(info_Halo,stream);
      BlocksRequiredMap->copyMetadata(info_brm, stream);
      vmesh::LocalID currSize = populations[popID].vmesh->size(); // Includes stream sync
      const vmesh::LocalID localContentBlocks = info_vbwcl->size;
      const vmesh::LocalID localNoContentBlocks = info_vbwncl->size;
      const vmesh::LocalID BlocksRequiredCapacity = info_Required->capacity;
      const vmesh::LocalID BlocksRequiredMapSizePower = info_brm->sizePower;
      vmesh::LocalID BlocksHaloCapacity = info_Halo->capacity;

      // Neighbour and own prefetches
      if (doPrefetches) {
         phiprof::start("Prefetch");
         populations[popID].vmesh->gpu_prefetchDevice(); // Queries active stream internally
         velocity_block_with_content_list->optimizeGPU(stream);
         velocity_block_with_no_content_list->optimizeGPU(stream);
         phiprof::stop("Prefetch");
      }

      phiprof::start("BlocksRequired hashmap resize / clear");
      // Estimate required size based on existing blocks
      vmesh::LocalID HashmapReqSize = 2;
      if (localContentBlocks+localNoContentBlocks > 0) {
         HashmapReqSize += ceil(log2(localContentBlocks+localNoContentBlocks));
      }
      if (BlocksRequiredMapSizePower >= HashmapReqSize) {
         // Map is already large enough
         BlocksRequiredMap->clear(Hashinator::targets::device,stream,false);
      } else {
         // Need larger empty map
         delete BlocksRequiredMap;
         BlocksRequiredMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(HashmapReqSize);
         int device = gpu_getDevice();
         BlocksRequiredMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         BlocksRequiredMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
         if ((attachedStream != 0)&&(needAttachedStreams)) {
            CHK_ERR( gpuStreamAttachMemAsync(attachedStream,BlocksRequiredMap, 0,gpuMemAttachSingle) );
            BlocksRequiredMap->streamAttach(attachedStream);
         }
         BlocksRequiredMap->optimizeGPU(stream);
      }
      BlocksHalo->clear();
      phiprof::stop("BlocksRequired hashmap resize / clear");

      if (localContentBlocks > 0) {
         // First add all local content blocks with a fast hashinator interface
         phiprof::start("Self Blocks with content");
         // 0.5 is target load factor
         BlocksRequiredMap->insert(velocity_block_with_content_list->data(),velocity_block_with_content_list->data(),localContentBlocks,0.5,stream,false);
         CHK_ERR( gpuPeekAtLastError() );
         CHK_ERR( gpuStreamSynchronize(stream) );
         phiprof::stop("Self Blocks with content");

         // add velocity space neighbors to map. We loop over blocks
         // with content, and insert all its v-space neighbors (halo)
         // Ensure at least one launch block
         nGpuBlocks = (localContentBlocks/GPUTHREADS) > GPUBLOCKS ? GPUBLOCKS : std::ceil((Real)localContentBlocks/(Real)GPUTHREADS);
         if (nGpuBlocks>0) {
            phiprof::start("Halo gather");
            BlocksHalo->clear();
            // Extreme estimate?
            vmesh::LocalID haloSizeEstimate = 8 * (localContentBlocks+localNoContentBlocks);
            if (BlocksHaloCapacity < haloSizeEstimate * BLOCK_ALLOCATION_FACTOR) {
               BlocksHalo->reserve(haloSizeEstimate * BLOCK_ALLOCATION_PADDING,true);
               BlocksHalo->memAdvise(gpuMemAdviseSetPreferredLocation,gpu_getDevice());
               BlocksHalo->memAdvise(gpuMemAdviseSetAccessedBy,gpu_getDevice());
               BlocksHalo->optimizeGPU(stream);
            }
            int addWidthV = getObjectWrapper().particleSpecies[popID].sparseBlockAddWidthV;
            update_blocks_required_halo_kernel<<<nGpuBlocks, GPUTHREADS, 0, stream>>> (
               populations[popID].vmesh,
               BlocksRequiredMap,
               velocity_block_with_content_list,
               BlocksHalo, // Now gathers blocks to be added into this vector
               addWidthV,
               // The following 4 vectors are passed just to be able to clear them on-device
               BlocksRequired,
               BlocksToRemove,
               BlocksToAdd,
               BlocksToMove
               );
            CHK_ERR( gpuPeekAtLastError() );
            CHK_ERR( gpuStreamSynchronize(stream) );
            phiprof::stop("Halo gather");
         }
         phiprof::start("Halo insert");
         BlocksHalo->copyMetadata(info_Halo,stream);
         CHK_ERR( gpuStreamSynchronize(stream) );
         const uint nHalo = info_Halo->size;
         if (nHalo > 0) {
            // 0.5 is target load factor
            BlocksRequiredMap->insert(BlocksHalo->data(),BlocksHalo->data(),nHalo,0.5,stream,false);
            CHK_ERR( gpuPeekAtLastError() );
            CHK_ERR( gpuStreamSynchronize(stream) );
         }
         phiprof::stop("Halo insert");
      }

      // add neighbor content info for spatial space neighbors to map. We loop over
      // neighbor cell lists with existing blocks, and raise the
      // flag for the local block with same block id

      // Here also gather to a vector, and then call the interface.
      // Re-use the halo vector.
      const uint neighbors_count = spatial_neighbors.size();
      if (neighbors_count > 0) {
         phiprof::start("Neighbor content lists");
         BlocksHalo->clear();
         CHK_ERR( gpuStreamSynchronize(stream) );
         for (std::vector<SpatialCell*>::const_iterator neighbor=spatial_neighbors.begin();
              neighbor != spatial_neighbors.end(); ++neighbor) {
            const int nNeighBlocks = (*neighbor)->velocity_block_with_content_list_size;
            // Ensure at least one launch block, try to do many neighbors at once
            nGpuBlocks = (nNeighBlocks/GPUTHREADS/neighbors_count) > GPUBLOCKS ? GPUBLOCKS : std::ceil((Real)nNeighBlocks/(Real)GPUTHREADS/(Real)neighbors_count);
            if (nGpuBlocks>0) {
               update_neighbours_have_content_kernel<<<nGpuBlocks, GPUTHREADS, 0, stream>>> (
                  populations[popID].vmesh,
                  BlocksHalo,
                  BlocksRequiredMap,
                  (*neighbor)->gpu_velocity_block_with_content_list_buffer,
                  nNeighBlocks
                  );
               CHK_ERR( gpuPeekAtLastError() );
            }
         }
         CHK_ERR( gpuStreamSynchronize(stream) );
         phiprof::stop("Neighbor content lists");
         phiprof::start("Neighbour Halo insert");
         const uint nHalo = BlocksHalo->size();
         if (nHalo > 0) {
            // 0.5 is target load factor
            BlocksRequiredMap->insert(BlocksHalo->data(),BlocksHalo->data(),nHalo,0.5,stream,false);
            CHK_ERR( gpuPeekAtLastError() );
            CHK_ERR( gpuStreamSynchronize(stream) );
         }
         phiprof::stop("Neighbour Halo insert");
      }
      // Same capacity for all
      phiprof::start("BlocksToXXX reserve");
      vmesh::LocalID growthOfBlocks = currSize * BLOCK_ALLOCATION_FACTOR;
      if (BlocksRequiredCapacity < growthOfBlocks) {
         BlocksRequired->reserve(currSize * BLOCK_ALLOCATION_PADDING,true);
         BlocksToAdd->reserve(currSize * BLOCK_ALLOCATION_PADDING,true);
         BlocksToRemove->reserve(currSize * BLOCK_ALLOCATION_PADDING,true);
         BlocksToMove->reserve(currSize * BLOCK_ALLOCATION_PADDING,true);
         int device = gpu_getDevice();
         BlocksRequired->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         BlocksToAdd->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         BlocksToRemove->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         BlocksToMove->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         BlocksRequired->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
         BlocksToAdd->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
         BlocksToRemove->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
         BlocksToMove->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      }
      SSYNC;
      phiprof::stop("BlocksToXXX reserve");
      phiprof::start("BlocksToXXX prefetch");
      if (doPrefetches || (BlocksRequiredCapacity < growthOfBlocks)) {
         BlocksRequired->optimizeGPU(stream);
         BlocksToRemove->optimizeGPU(stream);
         BlocksToAdd->optimizeGPU(stream);
         BlocksToMove->optimizeGPU(stream);
      }
      SSYNC;
      phiprof::stop("BlocksToXXX prefetch");

      // Extract list and count of all required blocks (content or with neighbors in spatial or velocity space)
      phiprof::start("Gather blocks required");
      const vmesh::LocalID nBlocksRequired = BlocksRequiredMap->extractAllKeys(*BlocksRequired,stream,false);
      phiprof::stop("Gather blocks required");

      // Flag all blocks in this cell without content + without neighbors with content to be removed
      if (doDeleteEmptyBlocks) {
         phiprof::start("Gather blocks to remove");
         // Ensure at least one launch block
         nGpuBlocks = (localNoContentBlocks/GPUTHREADS) > GPUBLOCKS ? GPUBLOCKS : std::ceil((Real)localNoContentBlocks/(Real)GPUTHREADS);
         if (nGpuBlocks>0) {
            update_blocks_to_remove_kernel<<<nGpuBlocks, GPUTHREADS, 0, stream>>> (
               velocity_block_with_no_content_list,
               BlocksRequiredMap,
               BlocksToRemove,
               localNoContentBlocks
               );
            CHK_ERR( gpuPeekAtLastError() );
            SSYNC;
         }
         phiprof::stop("Gather blocks to remove");
      }

      // Only add blocks which don't yet exist to optimize gpu parallel memory management.
      // Find these with a kernel.
      // This kernel also figures out which blocks need to be rescued from the end-space of the block data
      // Ensure at least one launch block
      nGpuBlocks = (nBlocksRequired/GPUTHREADS) > GPUBLOCKS ? GPUBLOCKS : std::ceil((Real)nBlocksRequired/(Real)GPUTHREADS);
      if (nBlocksRequired>0) {
         phiprof::start("blocks_to_add_kernel");
         update_blocks_to_add_kernel<<<nGpuBlocks, GPUTHREADS, 0, stream>>> (
            populations[popID].vmesh,
            BlocksRequired,
            BlocksToAdd,
            BlocksToMove,
            nBlocksRequired
            );
         CHK_ERR( gpuPeekAtLastError() );
         SSYNC;
         phiprof::stop("blocks_to_add_kernel");
      }

      // On-device adjustment calling happens in separate function as it is also called from within acceleration
      adjust_velocity_blocks_caller(popID);

      // Perform hashmap cleanup here (instead of at acceleration mid-steps)
      phiprof::start("Hashinator cleanup");
      if (needAttachedStreams) {
         populations[popID].vmesh->gpu_attachToStream(stream);
      }
      if (doPrefetches) {
         populations[popID].vmesh->gpu_prefetchDevice(stream);
      }
      populations[popID].vmesh->gpu_cleanHashMap(stream);
      SSYNC;
      phiprof::stop("Hashinator cleanup");

      phiprof::stop("Adjust velocity blocks");
   }

   void SpatialCell::adjust_velocity_blocks_caller(const uint popID) {
      /**
          Call GPU kernel with all necessary information for creation and deletion of blocks.
          Potential optimization: take the vector lengths as input parameters
          instead of having to call the size and then prefetch back to device.
      **/
      phiprof::start("GPU add and remove blocks");

#ifdef _OPENMP
      const uint thread_id = omp_get_thread_num();
#else
      const uint thread_id = 0;
#endif
      gpuStream_t stream = gpu_getStream();
      int nGpuBlocks;

      CHK_ERR( gpuStreamSynchronize(stream) ); // To ensure all previous kernels have finished

      phiprof::start("Block lists sizes");
      // Use copymetadata for these
      BlocksToAdd->copyMetadata(info_toAdd,stream);
      BlocksToRemove->copyMetadata(info_toRemove,stream);
      BlocksToMove->copyMetadata(info_toMove,stream);
      const vmesh::LocalID nBlocksBeforeAdjust = populations[popID].vmesh->size(); // includes a stream sync for the above
      const vmesh::LocalID nToAdd = info_toAdd->size;
      const vmesh::LocalID nToRemove = info_toRemove->size;
      //const vmesh::LocalID nToMove = info_toMove->size; // not used

      const vmesh::LocalID nBlocksAfterAdjust = nBlocksBeforeAdjust + nToAdd - nToRemove;
      const int nBlocksToChange = nToAdd > nToRemove ? nToAdd : nToRemove;
      nGpuBlocks = nBlocksToChange > GPUBLOCKS ? GPUBLOCKS : nBlocksToChange;
      phiprof::stop("Block lists sizes");

      // Grow the vectors, if necessary
      if (nBlocksAfterAdjust > nBlocksBeforeAdjust) {
         phiprof::start("GPU modify vmesh and VBC size (pre)");
         // These functions now prefetch back to device if necessary.
         populations[popID].vmesh->setNewSize(nBlocksAfterAdjust);
         populations[popID].blockContainer->setSize(nBlocksAfterAdjust);
         SSYNC;
         phiprof::stop("GPU modify vmesh and VBC size (pre)");
      }

      phiprof::start("GPU add and remove blocks kernel");
      if (nGpuBlocks>0) {
         CHK_ERR( gpuMemsetAsync(returnRealf[thread_id], 0, sizeof(Realf), stream) );
         CHK_ERR( gpuMemsetAsync(returnLID[thread_id], 0, 2*sizeof(vmesh::LocalID), stream) );
         dim3 block(WID,WID,WID);
         // Third argument specifies the number of bytes in *shared memory* that is
         // dynamically allocated per block for this call in addition to the statically allocated memory.
         update_velocity_blocks_kernel<<<nGpuBlocks, block, 0, stream>>> (
            populations[popID].vmesh,
            populations[popID].blockContainer,
            BlocksToAdd,
            BlocksToRemove,
            BlocksToMove,
            nBlocksBeforeAdjust,
            nBlocksAfterAdjust,
            returnLID[thread_id],//gpu_addVectorIndex and gpu_moveVectorIndex use these arrays
            returnRealf[thread_id]
            );
         CHK_ERR( gpuPeekAtLastError() );
         Realf host_rhoLossAdjust = 0;
         CHK_ERR( gpuMemcpyAsync(&host_rhoLossAdjust, returnRealf[thread_id], sizeof(Realf), gpuMemcpyDeviceToHost, stream) );
         CHK_ERR( gpuStreamSynchronize(stream) );
         this->populations[popID].RHOLOSSADJUST += host_rhoLossAdjust;
      }
      phiprof::stop("GPU add and remove blocks kernel");

      // Shrink the vectors, if necessary
      if (nBlocksAfterAdjust < nBlocksBeforeAdjust) {
         phiprof::start("GPU modify vmesh and VBC size (post)");
         // These functions now prefetch back to device if necessary.
         populations[popID].vmesh->setNewSize(nBlocksAfterAdjust);
         populations[popID].blockContainer->setSize(nBlocksAfterAdjust);
         SSYNC;
         phiprof::stop("GPU modify vmesh and VBC size (post)");
         if (doPrefetches) {
            phiprof::start("Vmesh and VBC lists prefetch dev");
            populations[popID].vmesh->gpu_prefetchDevice();
            populations[popID].blockContainer->gpu_prefetchDevice();
            SSYNC;
            phiprof::stop("Vmesh and VBC lists prefetch dev");
         }
      }

      // DEBUG output after kernel
      #ifdef DEBUG_SPATIAL_CELL
      phiprof::start("Vmesh and VBC debug output");
      populations[popID].vmesh->gpu_prefetchHost();
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID nAll = populations[popID].vmesh->size();
      printf("after kernel, size is %d should be %d\n",nAll,nBlocksAfterAdjust);
      for (vmesh::LocalID m=0; m<nAll; ++m) {
         const vmesh::GlobalID GIDs = populations[popID].vmesh->getGlobalID(m);
         const vmesh::LocalID LIDs = populations[popID].vmesh->getLocalID(GIDs);
         printf("LID %d GID-solved %d LID-solved %d\n",m,GIDs,LIDs);
      }
      populations[popID].vmesh->gpu_prefetchDevice();
      phiprof::start("Vmesh and VBC debug output");
      #endif

      // Don't return until everything is done?
      SSYNC;
      //CHK_ERR( gpuStreamSynchronize(stream) );
      phiprof::stop("GPU add and remove blocks");
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
      //space. TODO: should this delete blocks or not? Now not
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
            block_lengths.push_back(sizeof(vmesh::GlobalID) * populations[activePopID].vmesh->size());
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE1) !=0) {
            //Communicate size of list so that buffers can be allocated on receiving side
            if (!receiving) this->velocity_block_with_content_list_size = this->velocity_block_with_content_list->size();
            displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list_size) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE2) !=0) {
            if (receiving) {
               this->velocity_block_with_content_list->resize(this->velocity_block_with_content_list_size,true);
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
      populations[popID].vmesh->setGrid();
      populations[popID].blockContainer->setSize(populations[popID].vmesh->size());

      Real* parameters = get_block_parameters(popID);

      // Set velocity block parameters:
      for (vmesh::LocalID blockLID=0; blockLID<size(popID); ++blockLID) {
         const vmesh::GlobalID blockGID = get_velocity_block_global_id(blockLID,popID);
         parameters[BlockParams::VXCRD] = get_velocity_block_vx_min(popID,blockGID);
         parameters[BlockParams::VYCRD] = get_velocity_block_vy_min(popID,blockGID);
         parameters[BlockParams::VZCRD] = get_velocity_block_vz_min(popID,blockGID);
         populations[popID].vmesh->getCellSize(blockGID,&(parameters[BlockParams::DVX]));
         parameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
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

      phiprof::start("GPU update spatial cell block lists");
      gpuStream_t stream = gpu_getStream();
      phiprof::start("VB content list prefetches and allocations");
      // No obvious non-pagefaulting method for clearing?
      vmesh::LocalID currSize = populations[popID].vmesh->size();
      vmesh::LocalID currCapacity = velocity_block_with_content_list->capacity();
      velocity_block_with_content_list->clear();
      velocity_block_with_no_content_list->clear();
      if (currCapacity < currSize * BLOCK_ALLOCATION_FACTOR) {
         velocity_block_with_content_list->reserve(currSize * BLOCK_ALLOCATION_PADDING,true);
         velocity_block_with_no_content_list->reserve(currSize * BLOCK_ALLOCATION_PADDING,true);
         int device = gpu_getDevice();
         velocity_block_with_content_list->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         velocity_block_with_no_content_list->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         velocity_block_with_content_list->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
         velocity_block_with_no_content_list->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      }
      if (doPrefetches || (currCapacity < currSize)) {
         velocity_block_with_content_list->optimizeGPU(stream);
         velocity_block_with_no_content_list->optimizeGPU(stream);
      }
      phiprof::stop("VB content list prefetches and allocations");

      const Real velocity_block_min_value = getVelocityBlockMinValue(popID);

      phiprof::start("GPU update spatial cell block lists kernel");
      dim3 block(WID,WID,WID);
      // Third argument specifies the number of bytes in *shared memory* that is
      // dynamically allocated per block for this call in addition to the statically allocated memory.
      update_velocity_block_content_lists_kernel<<<GPUBLOCKS, block, WID3*sizeof(bool), stream>>> (
         populations[popID].vmesh,
         populations[popID].blockContainer,
         velocity_block_with_content_list,
         velocity_block_with_no_content_list,
         velocity_block_min_value
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuStreamSynchronize(stream) ); // This sync is required!
      phiprof::stop("GPU update spatial cell block lists kernel");

      // Note: Content list is not uploaded to device-only buffer here, but rather
      // in grid.cpp adjustVelocityBlocks()
      phiprof::stop("GPU update spatial cell block lists");
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
