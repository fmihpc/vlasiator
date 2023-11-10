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

// Non-AMR implementation of the velocity space, still in use despite the filename

#ifndef VELOCITY_MESH_GPU_H
#define VELOCITY_MESH_GPU_H

#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdint.h>
#include <vector>
#include <cmath>

#include "velocity_mesh_parameters.h"

#include "include/hashinator/hashinator.h"
#include "include/splitvector/splitvec.h"

#include "arch/gpu_base.hpp"
//#include "arch/arch_device_api.h" // included in above

#ifdef DEBUG_VLASIATOR
   #define DEBUG_VMESH
#endif

namespace vmesh {

   class VelocityMesh {
    public:
      VelocityMesh();
      ~VelocityMesh();
      VelocityMesh(const VelocityMesh& other);
      const VelocityMesh& operator=(const VelocityMesh& other);
      void gpu_destructor();

      size_t capacityInBytes() const;
      ARCH_HOSTDEV bool check() const;
      void clear();
      ARCH_HOSTDEV bool move(const vmesh::LocalID& sourceLocalID,const vmesh::LocalID& targetLocalID);
      ARCH_DEV bool warpMove(const vmesh::LocalID& sourceLocalID,const vmesh::LocalID& targetLocalID, const size_t b_tid);
      ARCH_HOSTDEV size_t count(const vmesh::GlobalID& globalID) const;
      ARCH_DEV size_t warpCount(const vmesh::GlobalID& globalID, const size_t b_tid) const;
      ARCH_HOSTDEV vmesh::GlobalID findBlock(vmesh::GlobalID cellIndices[3]) const;
      ARCH_DEV vmesh::GlobalID warpFindBlock(vmesh::GlobalID cellIndices[3], const size_t b_tid) const;
      ARCH_HOSTDEV bool getBlockCoordinates(const vmesh::GlobalID& globalID,Real coords[3]) const;
      ARCH_HOSTDEV void getBlockInfo(const vmesh::GlobalID& globalID,Real* array) const;
      ARCH_HOSTDEV const Real* getBlockSize() const;
      ARCH_HOSTDEV bool getBlockSize(const vmesh::GlobalID& globalID,Real size[3]) const;
      ARCH_HOSTDEV const Real* getCellSize(const uint8_t& refLevel=0) const;
      ARCH_HOSTDEV bool getCellSize(const vmesh::GlobalID& globalID,Real size[3]) const;
      ARCH_HOSTDEV vmesh::GlobalID getGlobalID(const vmesh::LocalID& localID) const;
      ARCH_HOSTDEV vmesh::GlobalID getGlobalID(const Real* coords) const;
      ARCH_HOSTDEV vmesh::GlobalID getGlobalID(vmesh::LocalID indices[3]) const;
      ARCH_HOSTDEV vmesh::GlobalID getGlobalID(const vmesh::LocalID& i,const vmesh::LocalID& j,const vmesh::LocalID& k) const;
      ARCH_HOSTDEV split::SplitVector<vmesh::GlobalID>& getGrid();
      ARCH_HOSTDEV const vmesh::LocalID* getGridLength(const uint8_t& refLevel=0) const;
      ARCH_HOSTDEV void getIndices(const vmesh::GlobalID& globalID,const uint8_t& refLevel,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const;
      ARCH_HOSTDEV void getIndices(const vmesh::GlobalID& globalID,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const;
      ARCH_HOSTDEV void getIndicesX(const vmesh::GlobalID& globalID,vmesh::LocalID& i) const;
      ARCH_HOSTDEV void getIndicesY(const vmesh::GlobalID& globalID,vmesh::LocalID& j) const;
      ARCH_HOSTDEV void getIndicesZ(const vmesh::GlobalID& globalID,vmesh::LocalID& k) const;
      ARCH_HOSTDEV size_t getMesh() const;
      ARCH_HOSTDEV vmesh::LocalID getLocalID(const vmesh::GlobalID& globalID) const;
      ARCH_DEV vmesh::LocalID warpGetLocalID(const vmesh::GlobalID& globalID, const size_t b_tid) const;
      ARCH_HOSTDEV const Real* getMeshMaxLimits() const;
      ARCH_HOSTDEV const Real* getMeshMinLimits() const;
      ARCH_HOSTDEV bool initialize(const size_t& meshID);
      ARCH_HOSTDEV static vmesh::LocalID invalidBlockIndex();
      ARCH_HOSTDEV static vmesh::GlobalID invalidGlobalID();
      ARCH_HOSTDEV static vmesh::LocalID invalidLocalID();
      ARCH_HOSTDEV bool isInitialized() const;
      ARCH_HOSTDEV void pop();
      ARCH_DEV void warpPop(const size_t b_tid);
      ARCH_HOSTDEV bool push_back(const vmesh::GlobalID& globalID);
      ARCH_DEV bool warpPush_back(const vmesh::GlobalID& globalID, const size_t b_tid);
      vmesh::LocalID push_back(const std::vector<vmesh::GlobalID>& blocks);
      ARCH_HOSTDEV vmesh::LocalID push_back(split::SplitVector<vmesh::GlobalID>* blocks);
      ARCH_DEV vmesh::LocalID warpPush_back(const split::SplitVector<vmesh::GlobalID>& blocks, const size_t b_tid);
      ARCH_DEV void replaceBlock(const vmesh::GlobalID& GIDold,const vmesh::LocalID& LID,const vmesh::GlobalID& GIDnew);
      ARCH_DEV void warpReplaceBlock(const vmesh::GlobalID& GIDold,const vmesh::LocalID& LID,const vmesh::GlobalID& GIDnew, const size_t b_tid);
      ARCH_DEV void placeBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID);
      ARCH_DEV void warpPlaceBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID, const size_t b_tid);
      ARCH_DEV void deleteBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID);
      ARCH_DEV void warpDeleteBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID, const size_t b_tid);
      void setGrid();
      bool setGrid(const std::vector<vmesh::GlobalID>& globalIDs);
      bool setGrid(const split::SplitVector<vmesh::GlobalID>& globalIDs);
      bool setMesh(const size_t& meshID);
      size_t getMeshID();
      void setNewSize(const vmesh::LocalID& newSize);
      void setNewCapacity(const vmesh::LocalID& newCapacity);
      ARCH_HOSTDEV size_t size() const;
      ARCH_HOSTDEV size_t sizeInBytes() const;
      void swap(VelocityMesh& vm);

      void gpu_prefetchHost(gpuStream_t stream);
      void gpu_prefetchDevice(gpuStream_t stream);
      void gpu_memAdvise(int device, gpuStream_t stream);
      void gpu_cleanHashMap(gpuStream_t stream);
      void gpu_attachToStream(gpuStream_t stream);
      void gpu_detachFromStream();
      void print_addresses();

   private:
      size_t meshID;
      gpuStream_t attachedStream;

      //std::vector<vmesh::GlobalID> localToGlobalMap;
      //OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID> globalToLocalMap;
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *globalToLocalMap;
      split::SplitVector<vmesh::GlobalID> *localToGlobalMap;
   };

   // ***** DEFINITIONS OF MEMBER FUNCTIONS ***** //


   inline VelocityMesh::VelocityMesh() {
      meshID = std::numeric_limits<size_t>::max();
      // Set sizepower to 10 (1024 blocks) straight away so there's enough room to grow
      globalToLocalMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      localToGlobalMap = new split::SplitVector<vmesh::GlobalID>(1);
      localToGlobalMap->clear();
      attachedStream = 0;
   }

   inline VelocityMesh::~VelocityMesh() {
      gpu_destructor();
   }
   inline void VelocityMesh::gpu_destructor() {
      if (globalToLocalMap) delete globalToLocalMap;
      if (localToGlobalMap) delete localToGlobalMap;
      globalToLocalMap = NULL;
      localToGlobalMap = NULL;
   }

   inline VelocityMesh::VelocityMesh(const VelocityMesh& other) {
      meshID = other.meshID;
      if (other.localToGlobalMap->size() > 0) {
         globalToLocalMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(*(other.globalToLocalMap));
         localToGlobalMap = new split::SplitVector<vmesh::GlobalID>(*(other.localToGlobalMap));
      } else {
         globalToLocalMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
         localToGlobalMap = new split::SplitVector<vmesh::GlobalID>(1);
         localToGlobalMap->clear();
      }
      attachedStream = 0;
      // gpuStream_t stream = gpu_getStream();
      // localToGlobalMap->optimizeJustDataGPU(stream);
      // globalToLocalMap->optimizeJustDataGPU(stream);
      // localToGlobalMap->optimizeMetadataGPU(stream);
      // globalToLocalMap->optimizeMetadataGPU(stream);
   }

   inline const VelocityMesh& VelocityMesh::operator=(const VelocityMesh& other) {
      //gpuStream_t stream = gpu_getStream();
      meshID = other.meshID;
      *globalToLocalMap = *(other.globalToLocalMap);
      *localToGlobalMap = *(other.localToGlobalMap);
      attachedStream = 0;
      // localToGlobalMap->optimizeJustDataGPU(stream);
      // globalToLocalMap->optimizeJustDataGPU(stream);
      // localToGlobalMap->optimizeMetadataGPU(stream);
      // globalToLocalMap->optimizeMetadataGPU(stream);
      return *this;
   }

   inline size_t VelocityMesh::capacityInBytes() const {
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      localToGlobalMap->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const size_t currentCapacity =  info_1[thread_id]->capacity;
      const size_t currentBucketCount = globalToLocalMap->bucket_count();

      const size_t capacityInBytes = currentCapacity*sizeof(vmesh::GlobalID)
           + currentBucketCount*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
      return capacityInBytes;
   }

   ARCH_HOSTDEV inline bool VelocityMesh::check() const {
      bool ok = true;
      #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      localToGlobalMap->optimizeJustDataCPU(stream);
      globalToLocalMap->optimizeJustDataCPU(stream);
      #endif

      if (localToGlobalMap->size() != globalToLocalMap->size()) {
         printf("VMO ERROR: sizes differ, %lu vs %lu\n",localToGlobalMap->size(),globalToLocalMap->size());
         ok = false;
         assert(0 && "VM check ERROR: sizes differ");
      }

      for (size_t b=0; b<size(); ++b) {
         const vmesh::LocalID globalID = localToGlobalMap->at(b);
         #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
         auto it = globalToLocalMap->device_find(globalID);
         #else
         auto it = globalToLocalMap->find(globalID);
         #endif
         const vmesh::GlobalID localID = it->second;
         if (localID != b) {
            ok = false;
            printf("VMO ERROR: localToGlobalMap[%lu] = %u but ",b,globalID);
            printf("globalToLocalMap[%u] = %u\n",globalID,localID);
            assert(0 && "VM check ERROR");
         }
      }
      #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      localToGlobalMap->optimizeJustDataGPU(stream);
      globalToLocalMap->optimizeJustDataGPU(stream);
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
      #endif
      return ok;
   }

   inline void VelocityMesh::clear() {
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      localToGlobalMap->clear();
      globalToLocalMap->clear(Hashinator::targets::device,stream,false);
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
   }

   ARCH_HOSTDEV inline bool VelocityMesh::move(const vmesh::LocalID& sourceLID,const vmesh::LocalID& targetLID) {
      const vmesh::GlobalID moveGID = localToGlobalMap->at(sourceLID); // block to move (at the end of list)
      const vmesh::GlobalID removeGID = localToGlobalMap->at(targetLID); // removed block
      #ifdef DEBUG_VMESH
      if (sourceLID != size()-1) {
         assert( 0 && "Error! Moving velocity mesh entry from position which is not last LID!");
      }
      #endif

      // at-function will throw out_of_range exception for non-existing global ID:
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      globalToLocalMap->set_element(moveGID,targetLID);
      globalToLocalMap->device_erase(removeGID);
      #else
      globalToLocalMap->at(moveGID) = targetLID;
      globalToLocalMap->erase(removeGID);
      #endif
      localToGlobalMap->at(targetLID) = moveGID;
      localToGlobalMap->pop_back();
      return true;
   }
   ARCH_HOSTDEV inline size_t VelocityMesh::count(const vmesh::GlobalID& globalID) const {
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      return globalToLocalMap->device_count(globalID);
      #else
      return globalToLocalMap->count(globalID);
      #endif
   }
   ARCH_HOSTDEV inline vmesh::GlobalID VelocityMesh::findBlock(vmesh::GlobalID cellIndices[3]) const {
      // Calculate i/j/k indices of the block that would own the cell:
      vmesh::GlobalID i_block = cellIndices[0] / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockLength[0];
      vmesh::GlobalID j_block = cellIndices[1] / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockLength[1];
      vmesh::GlobalID k_block = cellIndices[2] / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockLength[2];

      // Calculate block global ID:
      vmesh::GlobalID blockGID = getGlobalID(i_block,j_block,k_block);

      // If the block exists, return it:
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      if (globalToLocalMap->device_count(blockGID) != 0) {
      #else
      if (globalToLocalMap->count(blockGID) != 0) {
      #endif
         return blockGID;
      } else {
         return invalidGlobalID();
      }
   }

   ARCH_HOSTDEV inline bool VelocityMesh::getBlockCoordinates(const vmesh::GlobalID& globalID,Real coords[3]) const {
      if (globalID == invalidGlobalID()) {
         for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
         return false;
      }

      vmesh::LocalID indices[3];
      getIndices(globalID,indices[0],indices[1],indices[2]);
      if (indices[0] == invalidBlockIndex()) {
         for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
         return false;
      }

      coords[0] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[0] + indices[0]*(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[0];
      coords[1] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[1] + indices[1]*(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[1];
      coords[2] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[2] + indices[2]*(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[2];
      return true;
   }
   ARCH_HOSTDEV inline void VelocityMesh::getBlockInfo(const vmesh::GlobalID& globalID,Real* array) const {
      #ifdef DEBUG_VMESH
      if (globalID == invalidGlobalID()) {
         for (int i=0; i<6; ++i) array[i] = std::numeric_limits<Real>::infinity();
      }
      #endif

      vmesh::LocalID indices[3];
      indices[0] = globalID % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0];
      indices[1] = (globalID / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]) % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1];
      indices[2] = globalID / ((*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0] * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]);

      // Indices 0-2 contain coordinates of the lower left corner.
      // The values are the same as if getBlockCoordinates(globalID,&(array[0])) was called
      array[0] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[0] + indices[0]*(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[0];
      array[1] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[1] + indices[1]*(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[1];
      array[2] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[2] + indices[2]*(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[2];

      // Indices 3-5 contain the cell size.
      // The values are the same as if getCellSize(globalID,&(array[3])) was called
      array[3] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize[0];
      array[4] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize[1];
      array[5] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize[2];
   }

   ARCH_HOSTDEV inline const Real* VelocityMesh::getBlockSize() const {
      return (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize;
   }

   ARCH_HOSTDEV inline bool VelocityMesh::getBlockSize(const vmesh::GlobalID& globalID,Real size[3]) const {
      size[0] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[0];
      size[1] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[1];
      size[2] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[2];
      return true;
   }

   ARCH_HOSTDEV inline const Real* VelocityMesh::getCellSize(const uint8_t& refLevel) const {
      return (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize;
   }

   ARCH_HOSTDEV inline bool VelocityMesh::getCellSize(const vmesh::GlobalID& globalID,Real size[3]) const {
      size[0] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize[0];
      size[1] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize[1];
      size[2] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize[2];
      return true;
   }

   ARCH_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID& localID) const {
      #ifdef DEBUG_VMESH
      if (localID >= localToGlobalMap->size()) {
         assert (0 && "ERROR invalid local id");
      }
      #endif

      return localToGlobalMap->at(localID);
   }

   ARCH_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(const Real* coords) const {
      if (coords[0] < (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[0] || coords[0] >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMaxLimits[0] ||
         (coords[1] < (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[1] || coords[1] >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMaxLimits[1] ||
          coords[2] < (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[2] || coords[2] >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMaxLimits[2])) {
         return invalidGlobalID();
      }

      const vmesh::LocalID indices[3] = {
         static_cast<vmesh::LocalID>(floor((coords[0] - (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[0]) / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[0])),
         static_cast<vmesh::LocalID>(floor((coords[1] - (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[1]) / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[1])),
         static_cast<vmesh::LocalID>(floor((coords[2] - (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits[2]) / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[2]))
      };

      return indices[2]*(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]*(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]
              + indices[1]*(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0] + indices[0];
   }

   ARCH_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(vmesh::LocalID indices[3]) const {
      if (indices[0] >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]) return invalidGlobalID();
      if (indices[1] >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]) return invalidGlobalID();
      if (indices[2] >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[2]) return invalidGlobalID();
      return indices[0] + indices[1] * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]
         + indices[2] * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]
         * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0];
   }

   ARCH_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID& i,const vmesh::LocalID& j,const vmesh::LocalID& k) const {
      if (i >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]) return invalidGlobalID();
      if (j >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]) return invalidGlobalID();
      if (k >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[2]) return invalidGlobalID();
      return i + j * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]
         + k * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]
         * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0];
   }

   ARCH_HOSTDEV inline split::SplitVector<vmesh::GlobalID>& VelocityMesh::getGrid() {
      return *localToGlobalMap;
   }

   ARCH_HOSTDEV inline const vmesh::LocalID* VelocityMesh::getGridLength(const uint8_t& refLevel) const {
      return (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength;
   }

   ARCH_HOSTDEV inline void VelocityMesh::getIndices(const vmesh::GlobalID& globalID,const uint8_t& refLevel,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const {
      getIndices(globalID,i,j,k);
   }

   ARCH_HOSTDEV inline void VelocityMesh::getIndices(const vmesh::GlobalID& globalID,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const {
      if (globalID >= invalidGlobalID()) {
         i = j = k = invalidBlockIndex();
      } else {
         i = globalID % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0];
         j = (globalID / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]) % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1];
         k = globalID / ((*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0] * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]);
      }
   }

   ARCH_HOSTDEV inline void VelocityMesh::getIndicesX(const vmesh::GlobalID& globalID,vmesh::LocalID& i) const {
      if (globalID >= invalidGlobalID()) {
         i = invalidBlockIndex();
      } else {
         i = globalID % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0];
      }
   }
   ARCH_HOSTDEV inline void VelocityMesh::getIndicesY(const vmesh::GlobalID& globalID,vmesh::LocalID& j) const {
      if (globalID >= invalidGlobalID()) {
         j = invalidBlockIndex();
      } else {
         j = (globalID / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]) % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1];
      }
   }
   ARCH_HOSTDEV inline void VelocityMesh::getIndicesZ(const vmesh::GlobalID& globalID,vmesh::LocalID& k) const {
      if (globalID >= invalidGlobalID()) {
         k = invalidBlockIndex();
      } else {
         k = globalID / ((*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0] * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]);
      }
   }

   ARCH_HOSTDEV inline vmesh::LocalID VelocityMesh::getLocalID(const vmesh::GlobalID& globalID) const {
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      auto it = globalToLocalMap->device_find(globalID);
      if (it != globalToLocalMap->device_end()) return it->second;
      #else
      auto it = globalToLocalMap->find(globalID);
      if (it != globalToLocalMap->end()) return it->second;
      #endif
      return invalidLocalID();
   }
   ARCH_HOSTDEV inline size_t VelocityMesh::getMesh() const {
      return meshID;
   }

   ARCH_HOSTDEV inline const Real* VelocityMesh::getMeshMaxLimits() const {
      return (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMaxLimits;
   }

   ARCH_HOSTDEV inline const Real* VelocityMesh::getMeshMinLimits() const {
      return (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].meshMinLimits;
   }

   ARCH_HOSTDEV inline bool VelocityMesh::initialize(const size_t& meshID) {
      this->meshID = meshID;
      return true;
   }

   ARCH_HOSTDEV inline vmesh::LocalID VelocityMesh::invalidBlockIndex() {
      return INVALID_VEL_BLOCK_INDEX;
   }

   ARCH_HOSTDEV inline vmesh::GlobalID VelocityMesh::invalidGlobalID() {
      return INVALID_GLOBALID;
   }

   ARCH_HOSTDEV inline vmesh::LocalID VelocityMesh::invalidLocalID() {
      return INVALID_LOCALID;
   }

   ARCH_HOSTDEV inline bool VelocityMesh::isInitialized() const {
      return (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].initialized;
   }

   ARCH_HOSTDEV inline void VelocityMesh::pop() {
      const size_t mySize = size();
      if (mySize == 0) return;

      const vmesh::LocalID lastLID = mySize-1;
      const vmesh::GlobalID lastGID = localToGlobalMap->at(lastLID);
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      auto last = globalToLocalMap->device_find(lastGID);
      globalToLocalMap->device_erase(last);
      #else
      auto last = globalToLocalMap->find(lastGID);
      globalToLocalMap->erase(last);
      #endif
      localToGlobalMap->pop_back();
   }
   /** Returns true if added velocity block was new, false if it couldn't be added or already existed.
    */
   ARCH_HOSTDEV inline bool VelocityMesh::push_back(const vmesh::GlobalID& globalID) {
      const size_t mySize = size();
      if (mySize >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) return false;
      if (globalID == invalidGlobalID()) return false;

      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      // device_insert is slower, returns true or false for whether inserted key was new
      // globalToLocalMap->set_element(globalID,(vmesh::LocalID)mySize);
      // localToGlobalMap->device_push_back(globalID);
      auto position
         = globalToLocalMap->device_insert(Hashinator::make_pair(globalID,(vmesh::LocalID)mySize));
      if (position.second == true) {
         localToGlobalMap->device_push_back(globalID);
      }
      #else
      auto position
         = globalToLocalMap->insert(Hashinator::make_pair(globalID,(vmesh::LocalID)mySize));
      if (position.second == true) {
         localToGlobalMap->push_back(globalID);
      }
      #endif
      return position.second;
   }
   inline vmesh::LocalID VelocityMesh::push_back(const std::vector<vmesh::GlobalID>& blocks) {
      const size_t mySize = size();
      const size_t blocksSize = blocks.size();
      if (mySize+blocksSize > (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) {
         printf("vmesh: too many blocks, current size is %lu",mySize);
         printf(", adding %lu blocks", blocksSize);
         printf(", max is %u\n",(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks);
         return false;
      }

      for (size_t b=0; b<blocksSize; ++b) {
         globalToLocalMap->insert(Hashinator::make_pair(blocks[b],(vmesh::LocalID)(mySize+b)));
      }
      localToGlobalMap->insert(localToGlobalMap->end(),blocks.begin(),blocks.end());
      return blocksSize;
   }

   ARCH_HOSTDEV inline vmesh::LocalID VelocityMesh::push_back(split::SplitVector<vmesh::GlobalID>* blocks) {
     #if !(defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      //blocks->optimizeMetadataCPU(stream); // done in spatial cell
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      localToGlobalMap->copyMetadata(info_1[thread_id],stream);
      blocks->copyMetadata(info_2[thread_id],stream);
      localToGlobalMap->optimizeJustDataCPU(stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const size_t mySize = info_1[thread_id]->size;
      const size_t blocksSize = info_2[thread_id]->size;
      #else
      const size_t mySize = localToGlobalMap->size();
      const size_t blocksSize = blocks->size();
      #endif
      if (mySize+blocksSize > (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) {
         printf("vmesh: too many blocks, current size is %lu",mySize);
         printf(", adding %lu blocks", blocksSize);
         printf(", max is %u\n",(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks);
         return false;
      }

      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      for (size_t b=0; b<blocksSize; ++b) {
         // device_insert is slower than set_element, returns true or false for whether inserted key was new
         auto position
            = globalToLocalMap->device_insert(Hashinator::make_pair((*blocks)[b],(vmesh::LocalID)(mySize+b)));
         if (position.second == false) {
            printf("vmesh: failed to push_back new block! %d of %d\n",(uint)b,(uint)blocksSize);
            return b;
         }
         //globalToLocalMap->set_element(blocks[b],(vmesh::LocalID)(mySize+b));
      }
      localToGlobalMap->device_insert(localToGlobalMap->end(),blocks->begin(),blocks->end());
      #else
      if (mySize==0) {
         // Fast insertion into empty mesh
         localToGlobalMap->insert(localToGlobalMap->end(),blocks->begin(),blocks->end());
         CHK_ERR( gpuStreamSynchronize(stream) );
         localToGlobalMap->optimizeJustDataGPU(stream);
         globalToLocalMap->optimizeJustDataGPU(stream);
         CHK_ERR( gpuStreamSynchronize(stream) );
         globalToLocalMap->insertIndex(localToGlobalMap->data(),blocksSize,0.5,stream,false);
         CHK_ERR( gpuStreamSynchronize(stream) );
         localToGlobalMap->optimizeMetadataGPU(stream);
         globalToLocalMap->optimizeMetadataGPU(stream);
      } else {
         // GPUTODO: do inside kernel
         for (size_t b=0; b<blocksSize; ++b) {
            auto position
               = globalToLocalMap->insert(Hashinator::make_pair((*blocks)[b],(vmesh::LocalID)(mySize+b)));
            if (position.second == false) {
               printf("vmesh: failed to push_back new block! %d of %d\n",(uint)b,(uint)blocksSize);
               return b;
            }
         }
         localToGlobalMap->insert(localToGlobalMap->end(),blocks->begin(),blocks->end());
      }
      #endif
      return blocksSize;
   }

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
   /****
         Device-only accessors to be called from a single GPU thread
    **/
   ARCH_DEV inline void VelocityMesh::replaceBlock(const vmesh::GlobalID& GIDold,const vmesh::LocalID& LID,const vmesh::GlobalID& GIDnew) {
      #ifdef DEBUG_VMESH
      if (LID > size()-1) printf("vmesh replaceBlock error: LID is too large!\n");
      vmesh::LocalID LIDold = invalidLocalID();
      auto it = globalToLocalMap->device_find(GIDold);
      if (it != globalToLocalMap->device_end()) LIDold = it->second;
      if (localToGlobalMap->at(LIDold) != GIDold) printf("vmesh replaceBlock error: oldGID and oldLID don't match!\n");
      #endif
      globalToLocalMap->device_erase(GIDold);
      globalToLocalMap->set_element(GIDnew,LID);
      localToGlobalMap->at(LID) = GIDnew;
   }
   // Note: this function does not adjust the vmesh size, as it is used from within a parallel kernel.
   ARCH_DEV inline void VelocityMesh::placeBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID) {
      #ifdef DEBUG_VMESH
      if (LID > size()-1) {
         assert(0 && "vmesh placeBlock error: LID is too large!");
      }
      #endif
      globalToLocalMap->set_element(GID,LID);
      localToGlobalMap->at(LID) = GID;
   }
   ARCH_DEV inline void VelocityMesh::deleteBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID) {
      #ifdef DEBUG_VMESH
      // Verify that GID and LID match
      if (GID==invalidGlobalID()) printf("vmesh deleteBlock error: GID is invalidGlobalID!\n");
      if (LID==invalidLocalID()) printf("vmesh deleteBlock error: LID is invalidLocalID!\n");
      auto it = globalToLocalMap->device_find(GID);
      if (it == globalToLocalMap->device_end()) {
         printf("vmesh deleteBlock error: GID does not exist!\n");
      } else {
         if (it->second != LID) printf("vmesh deleteBlock error: LID %ul found with GID %ul does not match privided LID %ul!\n",it->second,GID,LID);
      }
      if (localToGlobalMap->at(LID) != GID) {
         printf("vmesh deleteBlock error: GID %ul found with LID %ul does not match privided GID %ul!\n",localToGlobalMap->at(LID),LID,GID);
      }
      #endif
      globalToLocalMap->device_erase(GID);
      localToGlobalMap->at(LID) = invalidGlobalID();
   }

   /****
        Warp accessor functions to be called from within GPU kernels over several threads
   **/
   ARCH_DEV inline void VelocityMesh::warpPop(const size_t b_tid) {
      const size_t mySize = size();
      if (mySize == 0) return;
      const vmesh::LocalID lastLID = mySize-1;
      const vmesh::GlobalID lastGID = localToGlobalMap->at(lastLID);
      #ifdef DEBUG_VMESH
      const vmesh::LocalID mapSize = globalToLocalMap->size();
      #endif
      if (b_tid < GPUTHREADS) {
         globalToLocalMap->warpErase(lastGID, b_tid);
      }
      if (b_tid==0) {
         localToGlobalMap->pop_back();
      }
      #ifdef DEBUG_VMESH
      __syncthreads();
      const vmesh::LocalID newMapSize = globalToLocalMap->size();
      const vmesh::LocalID newVecSize = size();
      if (newMapSize != mapSize-1) {
         printf("warpError in VelocityMesh::warpPop: map size %u is not expected %u! (thread %u)\n",newMapSize,(vmesh::LocalID)(mapSize-1),(vmesh::LocalID)b_tid);
         assert(0);
      }
      if (newVecSize != mySize-1) {
         printf("warpError in VelocityMesh::warpPop: vector size %u is not expected %u! (thread %u)\n",newVecSize,(vmesh::LocalID)(mySize-1),(vmesh::LocalID)b_tid);
         assert(0);
      }
      #endif
      __syncthreads();
   }
   ARCH_DEV inline vmesh::LocalID VelocityMesh::warpGetLocalID(const vmesh::GlobalID& globalID, const size_t b_tid) const {
      vmesh::LocalID retval = invalidLocalID();
      globalToLocalMap->warpFind(globalID, retval, b_tid % GPUTHREADS);
      #ifdef DEBUG_VMESH
      auto it = globalToLocalMap->device_find(globalID);
      if (it == globalToLocalMap->device_end()) {
         if (retval != invalidLocalID()) {
            printf("Warp error in VelocityMesh::warpGetLocalID: thread %u search did not find entry, warp search found LID %u for GID %u\n",(vmesh::LocalID)b_tid,retval,globalID);
         }
      } else if (retval != it->second) {
         printf("Warp error in VelocityMesh::warpGetLocalID: LID %u (warp) != %u (thread %u) for GID %u\n",
                retval,it->second,(vmesh::LocalID)b_tid,globalID);
      }
      #endif
      __syncthreads();
      return retval;
   }
   ARCH_DEV inline bool VelocityMesh::warpMove(const vmesh::LocalID& sourceLID,const vmesh::LocalID& targetLID, const size_t b_tid) {
      const vmesh::GlobalID moveGID = localToGlobalMap->at(sourceLID); // block to move (must at the end of list)
      const vmesh::GlobalID removeGID = localToGlobalMap->at(targetLID); // removed block
      #ifdef DEBUG_VMESH
      if (sourceLID != size()-1) {
         assert( 0 && "Error! Moving velocity mesh entry from position which is not last LID!");
      }
      const vmesh::LocalID preMapSize = globalToLocalMap->size();
      const vmesh::LocalID preVecSize = localToGlobalMap->size();
      #endif

      // at-function will throw out_of_range exception for non-existing global ID:
      if (b_tid < GPUTHREADS) {
         globalToLocalMap->warpInsert(moveGID,targetLID,b_tid);
         globalToLocalMap->warpErase(removeGID,b_tid);
      }
      if (b_tid==0) {
         localToGlobalMap->at(targetLID) = moveGID;
         localToGlobalMap->pop_back();
      }
      #ifdef DEBUG_VMESH
      __syncthreads();
      const vmesh::LocalID postMapSize = globalToLocalMap->size();
      const vmesh::LocalID postVecSize = localToGlobalMap->size();
      const vmesh::LocalID postLID = getLocalID(moveGID);
      if (preMapSize-1 != postMapSize) {
         printf("warpError in VelocityMesh::warpMove: map size %u is not expected %u! (thread %u)\n",postMapSize,preMapSize,(vmesh::LocalID)b_tid);
         assert(0);
      }
      if (preVecSize-1 != postVecSize) {
         printf("warpError in VelocityMesh::warpMove: vector size %u is not expected %u! (thread %u)\n",postVecSize,preVecSize,(vmesh::LocalID)b_tid);
         assert(0);
      }
      if (targetLID != postLID) {
         printf("warpError in VelocityMesh::warpMove: Entry at GID %u is %u instead of expected %u! (thread %u)\n",moveGID,postLID,targetLID,(vmesh::LocalID)b_tid);
         assert(0);
      }
      if (count(removeGID) != 0) {
         const vmesh::LocalID removedLIDfound = getLocalID(removeGID);
         printf("warpError in VelocityMesh::warpMove: Deleted GID %u still found in map with LID %u! (used to be LID %u) (thread %u)\n",removeGID,removedLIDfound,targetLID,(vmesh::LocalID)b_tid);
         assert(0);
      }
      #endif
      __syncthreads();
      return true;
   }
   ARCH_DEV inline size_t VelocityMesh::warpCount(const vmesh::GlobalID& globalID, const size_t b_tid) const {
      vmesh::LocalID retval = invalidLocalID();
      globalToLocalMap->warpFind(globalID, retval, b_tid % GPUTHREADS);
      #ifdef DEBUG_VMESH
      __syncthreads();
      vmesh::LocalID verify = getLocalID(globalID);
      if (verify != retval) {
         printf("warpError in VelocityMesh::warpCount: Searched GID %u found in map with LID %u for thread %u, but warpFind returned %u!\n",globalID,verify,(vmesh::LocalID)b_tid,retval);
         assert(0);
      }
      #endif
      __syncthreads();
      if (retval == invalidLocalID()) {
         return 0;
      } else {
         return 1;
      }
   }
   ARCH_DEV inline vmesh::GlobalID VelocityMesh::warpFindBlock(vmesh::GlobalID cellIndices[3], const size_t b_tid) const {
      // Calculate i/j/k indices of the block that would own the cell:
      vmesh::GlobalID i_block = cellIndices[0] / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockLength[0];
      vmesh::GlobalID j_block = cellIndices[1] / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockLength[1];
      vmesh::GlobalID k_block = cellIndices[2] / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockLength[2];

      // Calculate block global ID:
      vmesh::GlobalID blockGID = getGlobalID(i_block,j_block,k_block);

      // If the block exists, return it:
      vmesh::LocalID retval = invalidLocalID();
      globalToLocalMap->warpFind(blockGID, retval, b_tid % GPUTHREADS);

      #ifdef DEBUG_VMESH
      __syncthreads();
      auto it = globalToLocalMap->device_find(blockGID);
      if (it == globalToLocalMap->device_end()) {
         printf("Warp error in VelocityMesh::warpFindBlock: single-thread %u search did not find entry, warp search found LID %u for GID %u\n",(vmesh::LocalID)b_tid,retval,blockGID);
      } else if (retval != it->second) {
         printf("Warp error in VelocityMesh::warpFindBlock: LID %u (warp) != %u (thread %u) for GID %u\n",
                retval,it->second,(vmesh::LocalID)b_tid,blockGID);
      }
      #endif

      __syncthreads();
      if (retval == invalidLocalID()) {
         return invalidGlobalID();
      } else {
         return blockGID;
      }
   }
   ARCH_DEV inline bool VelocityMesh::warpPush_back(const vmesh::GlobalID& globalID, const size_t b_tid) {
      const vmesh::LocalID mySize = size();
      if (mySize >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) return false;
      if (globalID == invalidGlobalID()) return false;
      #ifdef DEBUG_VMESH
      const vmesh::LocalID mapSize = globalToLocalMap->size();
      #endif

      __shared__ bool inserted;
      if (b_tid==0) {
         inserted = false;
      }
      __syncthreads();
      if (b_tid < GPUTHREADS) {
         inserted = globalToLocalMap->warpInsert_V(globalID,(vmesh::LocalID)mySize, b_tid);
         if (inserted == true && b_tid==0) {
            localToGlobalMap->device_push_back(globalID);
         }
      }
      #ifdef DEBUG_VMESH
      __syncthreads();
      const vmesh::LocalID newMapSize = globalToLocalMap->size();
      const vmesh::LocalID newVecSize = size();
      const vmesh::LocalID postLID = getLocalID(globalID);
      const vmesh::GlobalID postGID = getGlobalID(mySize);
      if (newMapSize != mapSize+1) {
         printf("warpError in VelocityMesh::warpPush_back: map size %u is not expected %u! (thread %u)\n",newMapSize,mapSize+1,(vmesh::LocalID)b_tid);
         assert(0);
      }
      if (newVecSize != mySize+1) {
         printf("warpError in VelocityMesh::warpPush_back: vector size %u is not expected %u! (thread %u)\n",newVecSize,mySize+1,(vmesh::LocalID)b_tid);
         assert(0);
      }
      if (postLID != mySize) {
         printf("warpError in VelocityMesh::warpPush_back: hashmap returns LID %u but expected %u! (thread %u)\n",postLID,mySize,(vmesh::LocalID)b_tid);
         assert(0);
      }
      if (postLID != mySize) {
         printf("warpError in VelocityMesh::warpPush_back: vector entry is GID %u but expected %u! (thread %u)\n",postGID,globalID,(vmesh::LocalID)b_tid);
         assert(0);
      }
      #endif
      __syncthreads();
      return inserted;
   }
   ARCH_DEV inline vmesh::LocalID VelocityMesh::warpPush_back(const split::SplitVector<vmesh::GlobalID>& blocks, const size_t b_tid) {
      //GPUTODO: ADD debugs
      const vmesh::LocalID mySize = size();
      const vmesh::LocalID blocksSize = blocks.size();
      if (mySize+blocksSize > (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) {
         if (b_tid==0) {
            printf("vmesh: too many blocks, current size is %u",mySize);
            printf(", adding %u blocks", blocksSize);
            printf(", max is %u\n",(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks);
         }
         return false;
      }
      __shared__ vmesh::LocalID nInserted;
      __shared__ bool inserted;
      if (b_tid==0) {
         inserted = true;
      }
      __syncthreads();
      if (b_tid < GPUTHREADS) {
         for (vmesh::LocalID b=0; b<blocksSize; ++b) { // GPUTODO parallelize
            inserted = inserted && globalToLocalMap->warpInsert_V(blocks[b],(vmesh::LocalID)(mySize+b), b_tid);
            if (inserted == true && b_tid == 0) {
               nInserted = b;
            }
         }
      }
      __syncthreads();
      if (inserted == false) {
         if (b_tid == 0) {
            printf("vmesh: failed to push_back new block! %d of %d\n",(uint)(nInserted+1),(uint)blocksSize);
         }
         return nInserted+1;
      }
      if (b_tid==0) {
         localToGlobalMap->device_insert(localToGlobalMap->end(),blocks.begin(),blocks.end());
      }
      __syncthreads();
      return blocksSize;
   }
   ARCH_DEV inline void VelocityMesh::warpReplaceBlock(const vmesh::GlobalID& GIDold,const vmesh::LocalID& LID,const vmesh::GlobalID& GIDnew, const size_t b_tid) {
      // Inserts a (possibly new) block into the vmesh at a given position, and removes the existing block from there.
      #ifdef DEBUG_VMESH
      if (b_tid==0) {
         if (LID > size()-1) {
            printf("vmesh replaceBlock error: LID %u is too large for size %u!\n",LID,(vmesh::LocalID)size());
            assert(0);
         }
      }
      __syncthreads();
      vmesh::LocalID LIDold = invalidLocalID();
      globalToLocalMap->warpFind(GIDold, LIDold, b_tid % GPUTHREADS);
      auto it = globalToLocalMap->device_find(GIDold);
      if (it == globalToLocalMap->device_end()) {
         if (LIDold != invalidLocalID()) {
            printf("Warp error in VelocityMesh::warpReplaceBlock: thread %u search did not find entry, warp search found LID %u for GID %u\n",(vmesh::LocalID)b_tid,LIDold,GIDold);
            assert(0);
         }
      } else if (LIDold != it->second) {
         printf("Warp error in VelocityMesh::warpReplaceBlock: LID %u (warp) != %u (thread %u) for GID %u\n",
                LIDold,it->second,(vmesh::LocalID)b_tid,GIDold);
         assert(0);
      }
      if (b_tid==0) {
         if (localToGlobalMap->at(LIDold) != GIDold) {
            printf("vmesh replaceBlock error: oldGID and oldLID don't match!\n");
            assert(0);
         }
      }
      __syncthreads();
      const vmesh::LocalID preVecSize = size();
      const vmesh::LocalID preMapSize = globalToLocalMap->size();
      #endif
      if (b_tid < GPUTHREADS) { // GPUTODO these in parallel?
         globalToLocalMap->warpErase(GIDold, b_tid);

         #ifdef DEBUG_VMESH
         if (globalToLocalMap->size() != preMapSize-1) {
            printf("Warp error in VelocityMesh::warpReplaceBlock: map size %u does not match expected %u for thread %u!\n",(vmesh::LocalID)globalToLocalMap->size(),(vmesh::LocalID)(preMapSize-1),(vmesh::LocalID)b_tid);
            assert(0);
         }
         auto it2 = globalToLocalMap->device_find(GIDold);
         if (it2 != globalToLocalMap->device_end()) {
            printf("Warp error in VelocityMesh::warpReplaceBlock: warp-erased GID %u LID %u but thread %u still finds LID %u associated with it!\n",GIDold,LID,(vmesh::LocalID)b_tid,it2->second);
            assert(0);
         }
         bool newlyadded = false;
         newlyadded = globalToLocalMap->warpInsert_V(GIDnew,LID, b_tid);
         vmesh::LocalID postMapSize = globalToLocalMap->size();
         int change = 0;
         if (!newlyadded) change = -1;
         if (postMapSize != preMapSize + change) {
            printf("Warp error in VelocityMesh::warpReplaceBlock: map size %u does not match expected %u for thread %u!\n",postMapSize,(vmesh::LocalID)(preMapSize+change),(vmesh::LocalID)b_tid);
            if (b_tid==0) globalToLocalMap->stats();
            assert(0);
         }
         auto it3 = globalToLocalMap->device_find(GIDnew);
         if (it3 == globalToLocalMap->device_end()) {
            int sizePower = globalToLocalMap->getSizePower();
            printf("Warp error in VelocityMesh::warpReplaceBlock: warp-inserted GID %u LID %u but thread %u cannot find it!\n",GIDnew,LID,(vmesh::LocalID)b_tid);
            if (b_tid==0) {
               if (newlyadded) {
                  printf("warpAccessor reported true for insertion!\n");
               } else {
                  printf("warpAccessor reported false for insertion!\n");
               }
               globalToLocalMap->stats();
               //globalToLocalMap->dump_buckets();
            }
            assert(0);
         } else if (it3->second != LID) {
            printf("Warp error in VelocityMesh::warpReplaceBlock: warp-inserted GID %u LID %u but thread %u instead finds LID %u!\n",GIDnew,LID,(vmesh::LocalID)b_tid,it3->second);
            if (b_tid==0) globalToLocalMap->stats();
            assert(0);
         }
         #else
         bool newlyadded = globalToLocalMap->warpInsert_V(GIDnew,LID, b_tid);
         //globalToLocalMap->warpInsert(GIDnew,LID,b_tid);
         #endif
      }
      //__syncthreads(); // not needed
      if (b_tid==0) {
         localToGlobalMap->at(LID) = GIDnew;
      }
      __syncthreads();
   }
   ARCH_DEV inline void VelocityMesh::warpPlaceBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID, const size_t b_tid) {
      // Places block GID into the mesh with value LID. Assumes localToGlobalMap has already been grown sufficiently.
      #ifdef DEBUG_VMESH
      if (b_tid==0) {
         if (LID > size()-1) printf("vmesh placeBlock error: LID is too large!\n");
      }
      __syncthreads();
      #endif
      if (b_tid < GPUTHREADS) {
         bool newlyadded = false;
         newlyadded = globalToLocalMap->warpInsert_V(GID,LID, b_tid);
         if (b_tid==0) {
            if (!newlyadded) {
               globalToLocalMap->stats();
               //globalToLocalMap->dump_buckets();
            }
            assert(newlyadded && "newlyAdded warpPlaceBlock");
            localToGlobalMap->at(LID) = GID;
         }
      }
      #ifdef DEBUG_VMESH
      __syncthreads();
      // Note: no size check possible.
      auto it = globalToLocalMap->device_find(GID);
      if (it == globalToLocalMap->device_end()) {
         printf("Warp error in VelocityMesh::warpPlaceBlock: single-thread %u search did not find inserted GID %u LID %u\n",(vmesh::LocalID)b_tid,GID,LID);
         if (b_tid==0) globalToLocalMap->stats();
         assert(0);
      } else if (LID != it->second) {
         printf("Warp error in VelocityMesh::warpPlaceBlock: LID %u (warp) != %u (thread %u) for GID %u\n",LID,it->second,(vmesh::LocalID)b_tid,GID);
         if (b_tid==0) globalToLocalMap->stats();
         assert(0);
      }
      #endif
      __syncthreads();
   }
   ARCH_DEV inline void VelocityMesh::warpDeleteBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID, const size_t b_tid) {
      #ifdef DEBUG_VMESH
      // Verify that GID and LID match
      if (b_tid==0) {
         if (GID==invalidGlobalID()) printf("vmesh deleteBlock error: GID is invalidGlobalID!\n");
         if (LID==invalidLocalID()) printf("vmesh deleteBlock error: LID is invalidLocalID!\n");
      }
      __syncthreads();
      vmesh::LocalID retval = invalidLocalID();
      globalToLocalMap->warpFind(GID, retval, b_tid % GPUTHREADS);
      if (b_tid==0) {
         if (retval == invalidLocalID()) {
            printf("vmesh deleteBlock error: GID does not exist (warpFind(!\n");
         } else {
            if (retval != LID) {
               printf("vmesh deleteBlock error: LID %ul warpFound with GID %ul does not match provided LID %ul!\n",
                      retval,GID,LID);
            }
         }
         if (localToGlobalMap->at(LID) != GID) {
            printf("vmesh deleteBlock error: GID %ul warpFound with LID %ul does not match privided GID %ul!\n",
                   localToGlobalMap->at(LID),LID,GID);
         }
      }
      const vmesh::LocalID preMapSize = globalToLocalMap->size();
      __syncthreads();
      #endif
      if (b_tid < GPUTHREADS) {
         globalToLocalMap->warpErase(GID, b_tid);
         if (b_tid==0) {
            localToGlobalMap->at(LID) = invalidGlobalID();
         }
      }
      __syncthreads();
      #ifdef DEBUG_VMESH
      const vmesh::LocalID postMapSize = globalToLocalMap->size();
      if (postMapSize != preMapSize-1) {
         printf("Warp error in VelocityMesh::warpDeleteBlock: map size %u does not match expected %u for thread %u!\n",postMapSize,preMapSize-1,(vmesh::LocalID)b_tid);
         assert(0);
      }
      auto it = globalToLocalMap->device_find(GID);
      if (it != globalToLocalMap->device_end()) {
         printf("Warp error in VelocityMesh::warpDeleteBlock: GID %u still found with LID %u for thread %u!\n",GID,it->second,(vmesh::LocalID)b_tid);
         assert(0);
      }
      __syncthreads();
      #endif
   }
#endif

   inline void VelocityMesh::setGrid() {
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->clear(Hashinator::targets::device,stream,false);
      size_t nBlocks = localToGlobalMap->size();
      globalToLocalMap->insertIndex(localToGlobalMap->data(),nBlocks,0.5,stream,false);
      CHK_ERR( gpuStreamSynchronize(stream) );
   }

   // GPUTODO: These accessors are still slow, but we don't actually use them at all.
   inline bool VelocityMesh::setGrid(const std::vector<vmesh::GlobalID>& globalIDs) {
      printf("Warning! Slow version of VelocityMesh::setGrid.\n");
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->clear(Hashinator::targets::device,stream,false);
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap->insert(Hashinator::make_pair(globalIDs[i],(vmesh::LocalID)i));
      }
      localToGlobalMap->clear();
      localToGlobalMap->insert(localToGlobalMap->end(),globalIDs.begin(),globalIDs.end());
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
      return true;
   }
   inline bool VelocityMesh::setGrid(const split::SplitVector<vmesh::GlobalID>& globalIDs) {
      printf("Warning! Slow version of VelocityMesh::setGrid.\n");
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->clear(Hashinator::targets::device,stream,false);
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap->insert(Hashinator::make_pair(globalIDs[i],(vmesh::LocalID)i));
      }
      localToGlobalMap->clear();
      localToGlobalMap->insert(localToGlobalMap->end(),globalIDs.begin(),globalIDs.end());
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
      return true;
   }

   inline bool VelocityMesh::setMesh(const size_t& meshID) {
      if (meshID >= vmesh::getMeshWrapper()->velocityMeshes->size()) return false;
      this->meshID = meshID;
      return true;
   }
   inline size_t VelocityMesh::getMeshID() {
      return this->meshID;
   }

   inline void VelocityMesh::setNewSize(const vmesh::LocalID& newSize) {
      // Needed by GPU block adjustment
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      const uint device = gpu_getDevice();
      gpuStream_t stream = gpuStreamList[thread_id];
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      localToGlobalMap->copyMetadata(info_1[thread_id],stream);
      globalToLocalMap->copyMetadata(info_m[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      vmesh::LocalID currentCapacity =  info_1[thread_id]->capacity;
      const int currentSizePower = info_m[thread_id]->sizePower;
      // Passing eco flag = true to resize tells splitvector we manage padding manually.
      localToGlobalMap->resize(newSize,true);

      // CHK_ERR( gpuMemPrefetchAsync(localToGlobalMap, sizeof(split::SplitVector<vmesh::GlobalID>), gpuCpuDeviceId, stream) );
      // CHK_ERR( gpuMemPrefetchAsync(globalToLocalMap, sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>), gpuCpuDeviceId, stream) );
      if (newSize > currentCapacity) {
         // Was allocated new memory
         CHK_ERR( gpuStreamSynchronize(stream) );
         localToGlobalMap->optimizeJustDataGPU(stream);
         // localToGlobalMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         // localToGlobalMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      }
      // Ensure also that the map is large enough
      const int newSize2 = newSize > 0 ? newSize : 1;
      const int HashmapReqSize = ceil(log2(newSize2)) +2; // Make it really large enough
      if (currentSizePower < HashmapReqSize) {
         globalToLocalMap->device_rehash(HashmapReqSize, stream);
         CHK_ERR( gpuStreamSynchronize(stream) );
         globalToLocalMap->optimizeJustDataGPU(stream);
         // globalToLocalMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         // globalToLocalMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      }
      // Re-attach stream if required
      if ((attachedStream != 0)&&(needAttachedStreams)) {
         globalToLocalMap->streamAttach(attachedStream);
         localToGlobalMap->streamAttach(attachedStream);
      }
      // CHK_ERR( gpuMemPrefetchAsync(localToGlobalMap, sizeof(split::SplitVector<vmesh::GlobalID>), device, stream) );
      // CHK_ERR( gpuMemPrefetchAsync(globalToLocalMap, sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>), device, stream) );
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
   }

   // Used in initialization
   inline void VelocityMesh::setNewCapacity(const vmesh::LocalID& newCapacity) {
      // Passing eco flag = true to resize tells splitvector we manage padding manually.
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->copyMetadata(info_m[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      localToGlobalMap->reserve(newCapacity,true);
      // Ensure also that the map is large enough
      const int newCapacity2 = newCapacity > 0 ? newCapacity : 1;
      const int HashmapReqSize = ceil(log2(newCapacity2)) +1;
      if (info_m[thread_id]->sizePower < HashmapReqSize) {
         globalToLocalMap->resize(HashmapReqSize, Hashinator::targets::device, stream);
      }
      CHK_ERR( gpuStreamSynchronize(stream) );
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
   }

   ARCH_HOSTDEV inline size_t VelocityMesh::size() const {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      return localToGlobalMap->size();
#else
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      localToGlobalMap->optimizeMetadataCPU(stream);
      localToGlobalMap->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      return info_1[thread_id]->size;
      localToGlobalMap->optimizeMetadataGPU(stream);
#endif
   }

   ARCH_HOSTDEV inline size_t VelocityMesh::sizeInBytes() const {
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      const size_t currentSize = localToGlobalMap->size();
      const size_t currentFill = globalToLocalMap->size();
      #else
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      localToGlobalMap->copyMetadata(info_1[thread_id],stream);
      globalToLocalMap->copyMetadata(info_m[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const size_t currentSize = info_1[thread_id]->size;
      const size_t currentFill = info_m[thread_id]->fill;
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
      #endif
      const size_t sizeInBytes = currentSize*sizeof(vmesh::GlobalID)
           + currentFill*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
      return sizeInBytes;
   }

   inline void VelocityMesh::swap(VelocityMesh& vm) {
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      vm.localToGlobalMap->optimizeMetadataCPU(stream);
      vm.globalToLocalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->swap(*(vm.globalToLocalMap));
      localToGlobalMap->swap(*(vm.localToGlobalMap));
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
      vm.localToGlobalMap->optimizeMetadataGPU(stream);
      vm.globalToLocalMap->optimizeMetadataGPU(stream);
   }

   inline void VelocityMesh::gpu_prefetchHost(gpuStream_t stream=0) {
      if (stream==0) {
         stream = gpu_getStream();
      }
      phiprof::Timer prefetchTimer1 {"metadata prefetch2host"};
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      prefetchTimer1.stop();
      phiprof::Timer prefetchTimer2 {"justdata prefetch2host"};
      localToGlobalMap->optimizeJustDataCPU(stream);
      globalToLocalMap->optimizeJustDataCPU(stream);
      prefetchTimer2.stop();
      return;
   }

   inline void VelocityMesh::gpu_prefetchDevice(gpuStream_t stream=0) {
      if (stream==0) {
         stream = gpu_getStream();
      }
      phiprof::Timer prefetchTimer1 {"metadata prefetch2device"};
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      prefetchTimer1.stop();
      phiprof::Timer prefetchTimer2 {"justdata prefetch2device"};
      localToGlobalMap->optimizeJustDataGPU(stream);
      globalToLocalMap->optimizeJustDataGPU(stream);
      prefetchTimer2.stop();
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
      return;
   }

   inline void VelocityMesh::gpu_memAdvise(int device, gpuStream_t stream) {
      // int device = gpu_getDevice();
      // globalToLocalMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      // localToGlobalMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      // globalToLocalMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      // localToGlobalMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      return;
   }

   inline void VelocityMesh::gpu_cleanHashMap(gpuStream_t stream = 0) {
      if (stream==0) {
         stream = gpu_getStream();
      }
      const uint device = gpu_getDevice();
      localToGlobalMap->optimizeMetadataCPU(stream);
      globalToLocalMap->optimizeMetadataCPU(stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      globalToLocalMap->performCleanupTasks(stream);
      globalToLocalMap->resize_to_lf(0.5, Hashinator::targets::device, stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      localToGlobalMap->optimizeMetadataGPU(stream);
      globalToLocalMap->optimizeMetadataGPU(stream);
      return;
   }

   inline void VelocityMesh::gpu_attachToStream(gpuStream_t stream = 0) {
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
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,globalToLocalMap, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,localToGlobalMap, 0,gpuMemAttachSingle) );
      globalToLocalMap->streamAttach(attachedStream);
      localToGlobalMap->streamAttach(attachedStream);
      return;
   }
   inline void VelocityMesh::gpu_detachFromStream() {
      // Return if attaching is not needed
      if (!needAttachedStreams) {
         return;
      }
      // Detach unified memory regions from streams
      if (attachedStream == 0) {
         // Already detached
         return;
      }
      attachedStream = 0;
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,globalToLocalMap, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,localToGlobalMap, 0,gpuMemAttachGlobal) );
      globalToLocalMap->streamAttach(0,gpuMemAttachGlobal);
      localToGlobalMap->streamAttach(0,gpuMemAttachGlobal);
      return;
   }
   inline void VelocityMesh::print_addresses() {
      printf("GPU localToGlobalMap %p\n GPU globalToLocalMap %p\n",localToGlobalMap,globalToLocalMap);
   }

} // namespace vmesh

#endif
