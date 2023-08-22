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

namespace vmesh {

   class VelocityMesh : public Managed {
    public:
      VelocityMesh();
      ~VelocityMesh();
      VelocityMesh(const VelocityMesh& other);
      const VelocityMesh& operator=(const VelocityMesh& other);

      ARCH_HOSTDEV size_t capacityInBytes() const;
      ARCH_HOSTDEV bool check() const;
      void clear();
      ARCH_HOSTDEV bool copy(const vmesh::LocalID& sourceLocalID,const vmesh::LocalID& targetLocalID);
      ARCH_HOSTDEV size_t count(const vmesh::GlobalID& globalID) const;
      ARCH_HOSTDEV vmesh::GlobalID findBlock(vmesh::GlobalID cellIndices[3]) const;
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
      ARCH_HOSTDEV const Real* getMeshMaxLimits() const;
      ARCH_HOSTDEV const Real* getMeshMinLimits() const;
      ARCH_HOSTDEV bool initialize(const size_t& meshID);
      ARCH_HOSTDEV static vmesh::LocalID invalidBlockIndex();
      ARCH_HOSTDEV static vmesh::GlobalID invalidGlobalID();
      ARCH_HOSTDEV static vmesh::LocalID invalidLocalID();
      ARCH_HOSTDEV bool isInitialized() const;
      ARCH_HOSTDEV void pop();
      ARCH_HOSTDEV bool push_back(const vmesh::GlobalID& globalID);
      bool push_back(const std::vector<vmesh::GlobalID>& blocks);
      ARCH_HOSTDEV bool push_back(const split::SplitVector<vmesh::GlobalID>& blocks);
      ARCH_DEV void replaceBlock(const vmesh::GlobalID& GIDold,const vmesh::LocalID& LID,const vmesh::GlobalID& GIDnew);
      ARCH_DEV void placeBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID);
      ARCH_DEV void deleteBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID);
      void setGrid();
      bool setGrid(const std::vector<vmesh::GlobalID>& globalIDs);
      bool setGrid(const split::SplitVector<vmesh::GlobalID>& globalIDs);
      bool setMesh(const size_t& meshID);
      void setNewSize(const vmesh::LocalID& newSize);
      ARCH_HOSTDEV size_t size() const;
      ARCH_HOSTDEV size_t sizeInBytes() const;
      ARCH_HOSTDEV void swap(VelocityMesh& vm);

      void gpu_prefetchHost(gpuStream_t stream);
      void gpu_prefetchDevice(gpuStream_t stream);
      void gpu_memAdvise(int device, gpuStream_t stream);
      void gpu_cleanHashMap(gpuStream_t stream);
      void gpu_attachToStream(gpuStream_t stream);
      void gpu_detachFromStream();

   private:
      size_t meshID;
      gpuStream_t attachedStream;

      //std::vector<vmesh::GlobalID> localToGlobalMap;
      //OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID> globalToLocalMap;
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *globalToLocalMap;
      split::SplitVector<vmesh::GlobalID> *localToGlobalMap;
      split::SplitInfo *info_ltgm;
   };

   // ***** DEFINITIONS OF MEMBER FUNCTIONS ***** //


   inline VelocityMesh::VelocityMesh() {
      meshID = std::numeric_limits<size_t>::max();
      // Set sizepower to 10 (1024 blocks) straight away so there's enough room to grow
      globalToLocalMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      localToGlobalMap = new split::SplitVector<vmesh::GlobalID>(1);
      localToGlobalMap->clear();
      attachedStream = 0;
      gpuMallocHost((void **) &info_ltgm, sizeof(split::SplitInfo));
   }

   inline VelocityMesh::~VelocityMesh() {
      delete globalToLocalMap;
      delete localToGlobalMap;
      gpuFreeHost(info_ltgm);
   }

   inline VelocityMesh::VelocityMesh(const VelocityMesh& other) {
      meshID = other.meshID;
      globalToLocalMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(*(other.globalToLocalMap));
      if (other.localToGlobalMap->size() > 0) {
         localToGlobalMap = new split::SplitVector<vmesh::GlobalID>(*(other.localToGlobalMap));
      } else {
         localToGlobalMap = new split::SplitVector<vmesh::GlobalID>(1);
         localToGlobalMap->clear();
      }
      attachedStream = 0;
      gpuMallocHost((void **) &info_ltgm, sizeof(split::SplitInfo));
   }

   inline const VelocityMesh& VelocityMesh::operator=(const VelocityMesh& other) {
      delete globalToLocalMap;
      delete localToGlobalMap;
      meshID = other.meshID;
      globalToLocalMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(*(other.globalToLocalMap));
      if (other.localToGlobalMap->size() > 0) {
         localToGlobalMap = new split::SplitVector<vmesh::GlobalID>(*(other.localToGlobalMap));
      } else {
         localToGlobalMap = new split::SplitVector<vmesh::GlobalID>(1);
         localToGlobalMap->clear();
      }
      attachedStream = 0;
      //gpuMallocHost((void **) &info_ltgm, sizeof(split::SplitInfo)); already exists
      return *this;
   }

   ARCH_HOSTDEV inline size_t VelocityMesh::capacityInBytes() const {
      return localToGlobalMap->capacity()*sizeof(vmesh::GlobalID)
           + globalToLocalMap->bucket_count()*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
   }

   ARCH_HOSTDEV inline bool VelocityMesh::check() const {
      bool ok = true;

      if (localToGlobalMap->size() != globalToLocalMap->size()) {
         printf("VMO ERROR: sizes differ, %lu vs %lu\n",localToGlobalMap->size(),globalToLocalMap->size());
         ok = false;
         exit(1);
      }

      for (size_t b=0; b<size(); ++b) {
         const vmesh::LocalID globalID = localToGlobalMap->at(b);
         auto it = globalToLocalMap->find(globalID);
         const vmesh::GlobalID localID = it->second;
         if (localID != b) {
            ok = false;
            printf("VMO ERROR: localToGlobalMap[%lu] = %u but ",b,globalID);
            printf("globalToLocalMap[%u] = %u\n",globalID,localID);
            exit(1);
         }
      }

      return ok;
   }

   inline void VelocityMesh::clear() {
      split::SplitVector<vmesh::GlobalID>().swap(*localToGlobalMap);
      globalToLocalMap->clear();
   }

   ARCH_HOSTDEV inline bool VelocityMesh::copy(const vmesh::LocalID& sourceLID,const vmesh::LocalID& targetLID) {
      const vmesh::GlobalID sourceGID = localToGlobalMap->at(sourceLID); // block at the end of list
      const vmesh::GlobalID targetGID = localToGlobalMap->at(targetLID); // removed block

      // at-function will throw out_of_range exception for non-existing global ID:
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      globalToLocalMap->set_element(sourceGID,targetLID);
      globalToLocalMap->set_element(targetGID,sourceLID);
      #else
      globalToLocalMap->at(sourceGID) = targetLID;
      globalToLocalMap->at(targetGID) = sourceLID;
      #endif
      localToGlobalMap->at(targetLID) = sourceGID;
      localToGlobalMap->at(sourceLID) = targetGID;
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
      #ifndef NDEBUG
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
      #ifndef NDEBUG
      if (localID >= localToGlobalMap->size()) {
         printf("ERROR invalid local id %lu\n",localID);
         exit(1);
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
      if (size() == 0) return;

      const vmesh::LocalID lastLID = size()-1;
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

   ARCH_HOSTDEV inline bool VelocityMesh::push_back(const vmesh::GlobalID& globalID) {
      if (size() >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) return false;
      if (globalID == invalidGlobalID()) return false;

      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      auto position
         = globalToLocalMap->device_insert(Hashinator::make_pair(globalID,(vmesh::LocalID)localToGlobalMap->size()));
      if (position.second == true) {
         localToGlobalMap->device_push_back(globalID);
      }
      #else
      auto position
         = globalToLocalMap->insert(Hashinator::make_pair(globalID,(vmesh::LocalID)localToGlobalMap->size()));
      if (position.second == true) {
         localToGlobalMap->push_back(globalID);
      }
      #endif
      return position.second;
   }

   inline bool VelocityMesh::push_back(const std::vector<vmesh::GlobalID>& blocks) {
      if (size()+blocks.size() > (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) {
         printf("vmesh: too many blocks, current size is %lu",size());
         printf(", adding %lu blocks", blocks.size());
         printf(", max is %u\n",(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks);
         return false;
      }

      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap->insert(Hashinator::make_pair(blocks[b],(vmesh::LocalID)(localToGlobalMap->size()+b)));
      }
      localToGlobalMap->insert(localToGlobalMap->end(),blocks.begin(),blocks.end());
      return true;
   }

   ARCH_HOSTDEV inline bool VelocityMesh::push_back(const split::SplitVector<vmesh::GlobalID>& blocks) {
      if (size()+blocks.size() > (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) {
         printf("vmesh: too many blocks, current size is %lu",size());
         printf(", adding %lu blocks", blocks.size());
         printf(", max is %u\n",(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks);
         return false;
      }

      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap->device_insert(Hashinator::make_pair(blocks[b],(vmesh::LocalID)(localToGlobalMap->size()+b)));
      }
      localToGlobalMap->device_insert(localToGlobalMap->end(),blocks.begin(),blocks.end());
      #else
      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap->insert(Hashinator::make_pair(blocks[b],(vmesh::LocalID)(localToGlobalMap->size()+b)));
      }
      localToGlobalMap->insert(localToGlobalMap->end(),blocks.begin(),blocks.end());
      #endif
      return true;
   }

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
   ARCH_DEV inline void VelocityMesh::replaceBlock(const vmesh::GlobalID& GIDold,const vmesh::LocalID& LID,const vmesh::GlobalID& GIDnew) {
      globalToLocalMap->device_erase(GIDold);
      globalToLocalMap->set_element(GIDnew,LID);
      localToGlobalMap->at(LID) = GIDnew;
   }
   ARCH_DEV inline void VelocityMesh::placeBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID) {
      globalToLocalMap->set_element(GID,LID);
      localToGlobalMap->at(LID) = GID;
   }
   ARCH_DEV inline void VelocityMesh::deleteBlock(const vmesh::GlobalID& GID,const vmesh::LocalID& LID) {
      globalToLocalMap->device_erase(GID);
      localToGlobalMap->at(LID) = invalidGlobalID();
   }
#endif

   inline void VelocityMesh::setGrid() {
      globalToLocalMap->clear();
      for (size_t i=0; i<localToGlobalMap->size(); ++i) {
         globalToLocalMap->insert(Hashinator::make_pair(localToGlobalMap->at(i),(vmesh::LocalID)i));
      }
   }

   inline bool VelocityMesh::setGrid(const std::vector<vmesh::GlobalID>& globalIDs) {
      globalToLocalMap->clear();
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap->insert(Hashinator::make_pair(globalIDs[i],(vmesh::LocalID)i));
      }
      localToGlobalMap->clear();
      localToGlobalMap->insert(localToGlobalMap->end(),globalIDs.begin(),globalIDs.end());
      return true;
   }
   inline bool VelocityMesh::setGrid(const split::SplitVector<vmesh::GlobalID>& globalIDs) {
      globalToLocalMap->clear();
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap->insert(Hashinator::make_pair(globalIDs[i],(vmesh::LocalID)i));
      }
      localToGlobalMap->clear();
      localToGlobalMap->insert(localToGlobalMap->end(),globalIDs.begin(),globalIDs.end());
      return true;
   }

   inline bool VelocityMesh::setMesh(const size_t& meshID) {
      if (meshID >= vmesh::getMeshWrapper()->velocityMeshes->size()) return false;
      this->meshID = meshID;
      return true;
   }

   inline void VelocityMesh::setNewSize(const vmesh::LocalID& newSize) {
      // Needed by GPU block adjustment
      // Passing eco flag = true to resize tells splitvector we manage padding manually.
      vmesh::LocalID currentCapacity = localToGlobalMap->capacity();
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap->resize(newSize,true);
      int device = gpu_getDevice();
      if (newSize > currentCapacity) {
         // Was allocated new memory
         CHK_ERR( gpuStreamSynchronize(stream) );
         localToGlobalMap->optimizeGPU(stream);
         localToGlobalMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         localToGlobalMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      }
      // Ensure also that the map is large enough
      const int HashmapReqSize = ceil(log2(newSize)) +2; // Make it really large enough
      if (globalToLocalMap->getSizePower() < HashmapReqSize) {
         globalToLocalMap->device_rehash(HashmapReqSize, stream);
         CHK_ERR( gpuStreamSynchronize(stream) );
         globalToLocalMap->optimizeGPU(stream);
         globalToLocalMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         globalToLocalMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      }
      // Re-attach stream if required
      if ((attachedStream != 0)&&(needAttachedStreams)) {
         globalToLocalMap->streamAttach(attachedStream);
         localToGlobalMap->streamAttach(attachedStream);
      }
   }

   ARCH_HOSTDEV inline size_t VelocityMesh::size() const {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      return localToGlobalMap->size();
#else
      // Host-side non-pagefaulting approach
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap->copyMetadata(info_ltgm,stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      return info_ltgm->size;
#endif
   }

   ARCH_HOSTDEV inline size_t VelocityMesh::sizeInBytes() const {
      return globalToLocalMap->size()*sizeof(vmesh::GlobalID)
           + localToGlobalMap->size()*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
   }

   ARCH_HOSTDEV inline void VelocityMesh::swap(VelocityMesh& vm) {
      globalToLocalMap->swap(*(vm.globalToLocalMap));
      localToGlobalMap->swap(*(vm.localToGlobalMap));
   }

   inline void VelocityMesh::gpu_prefetchHost(gpuStream_t stream=0) {
      //if (localToGlobalMap->size() == 0) return; // This size check in itself causes a page fault
      // In fact we only need to prefetch the buckets inside the hashmap to GPU, but use this call.
      //Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *uploaded = globalToLocalMap->upload(gpu_getStream());
      if (stream==0) {
         globalToLocalMap->optimizeCPU(gpu_getStream());
         localToGlobalMap->optimizeCPU(gpu_getStream());
      } else {
         globalToLocalMap->optimizeCPU(stream);
         localToGlobalMap->optimizeCPU(stream);
      }
      return;
   }

   inline void VelocityMesh::gpu_prefetchDevice(gpuStream_t stream=0) {
      //if (localToGlobalMap->size() == 0) return; // This size check in itself causes a page fault
      // In fact we only need to prefetch the buckets inside the hashmap to GPU, but use this call.
      //Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *uploaded = globalToLocalMap->upload(gpu_getStream());
      if (stream==0) {
         globalToLocalMap->optimizeGPU(gpu_getStream());
         localToGlobalMap->optimizeGPU(gpu_getStream());
      } else {
         globalToLocalMap->optimizeGPU(stream);
         localToGlobalMap->optimizeGPU(stream);
      }
      return;
   }

   inline void VelocityMesh::gpu_memAdvise(int device, gpuStream_t stream) {
      // int device = gpu_getDevice();
      globalToLocalMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      localToGlobalMap->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      globalToLocalMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      localToGlobalMap->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      return;
   }

   inline void VelocityMesh::gpu_cleanHashMap(gpuStream_t stream = 0) {
      if (stream==0) {
         globalToLocalMap->performCleanupTasks(gpu_getStream());
      } else {
         globalToLocalMap->performCleanupTasks(stream);
      }
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
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,this, 0,gpuMemAttachSingle) );
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
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,this, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,globalToLocalMap, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,localToGlobalMap, 0,gpuMemAttachGlobal) );
      globalToLocalMap->streamAttach(0,gpuMemAttachGlobal);
      localToGlobalMap->streamAttach(0,gpuMemAttachGlobal);
      return;
   }

} // namespace vmesh

#endif
