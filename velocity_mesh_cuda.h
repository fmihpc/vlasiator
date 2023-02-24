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

#ifndef VELOCITY_MESH_CUDA_H
#define VELOCITY_MESH_CUDA_H

#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdint.h>
#include <vector>
//#include <unordered_map>
//#include <set>
#include <cmath>

#include "velocity_mesh_parameters.h"

//#include "object_wrapper.h"
//#include "open_bucket_hashtable.h"
#include "include/hashinator/hashinator.h"
#include "include/splitvector/splitvec.h"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cuda_context.cuh" // managed class, CUDA_HOSTDEV
#include <cuda/std/utility> // cuda::std::pair

namespace vmesh {

   class VelocityMesh : public Managed {
    public:
      VelocityMesh();
      ~VelocityMesh();
      VelocityMesh(const VelocityMesh& other);
      const VelocityMesh& operator=(const VelocityMesh& other);

      CUDA_HOSTDEV size_t capacityInBytes() const;
      CUDA_HOSTDEV bool check() const;
      void clear();
      CUDA_HOSTDEV bool copy(const vmesh::LocalID& sourceLocalID,const vmesh::LocalID& targetLocalID);
      CUDA_HOSTDEV size_t count(const vmesh::GlobalID& globalID) const;
      CUDA_HOSTDEV vmesh::GlobalID findBlock(vmesh::GlobalID cellIndices[3]) const;
      CUDA_HOSTDEV bool getBlockCoordinates(const vmesh::GlobalID& globalID,Real coords[3]) const;
      CUDA_HOSTDEV void getBlockInfo(const vmesh::GlobalID& globalID,Real* array) const;
      CUDA_HOSTDEV const Real* getBlockSize() const;
      CUDA_HOSTDEV bool getBlockSize(const vmesh::GlobalID& globalID,Real size[3]) const;
      CUDA_HOSTDEV const Real* getCellSize(const uint8_t& refLevel=0) const;
      CUDA_HOSTDEV bool getCellSize(const vmesh::GlobalID& globalID,Real size[3]) const;
      CUDA_HOSTDEV vmesh::GlobalID getGlobalID(const vmesh::LocalID& localID) const;
      CUDA_HOSTDEV vmesh::GlobalID getGlobalID(const Real* coords) const;
      CUDA_HOSTDEV vmesh::GlobalID getGlobalID(vmesh::LocalID indices[3]) const;
      CUDA_HOSTDEV vmesh::GlobalID getGlobalID(const vmesh::LocalID& i,const vmesh::LocalID& j,const vmesh::LocalID& k) const;
      CUDA_HOSTDEV split::SplitVector<vmesh::GlobalID>& getGrid();
      CUDA_HOSTDEV const vmesh::LocalID* getGridLength(const uint8_t& refLevel=0) const;
      CUDA_HOSTDEV void getIndices(const vmesh::GlobalID& globalID,const uint8_t& refLevel,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const;
      CUDA_HOSTDEV void getIndices(const vmesh::GlobalID& globalID,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const;
      CUDA_HOSTDEV size_t getMesh() const;
      CUDA_HOSTDEV vmesh::LocalID getLocalID(const vmesh::GlobalID& globalID) const;
      CUDA_HOSTDEV const Real* getMeshMaxLimits() const;
      CUDA_HOSTDEV const Real* getMeshMinLimits() const;
      CUDA_HOSTDEV bool initialize(const size_t& meshID);
      CUDA_HOSTDEV static vmesh::LocalID invalidBlockIndex();
      CUDA_HOSTDEV static vmesh::GlobalID invalidGlobalID();
      CUDA_HOSTDEV static vmesh::LocalID invalidLocalID();
      CUDA_HOSTDEV bool isInitialized() const;
      CUDA_HOSTDEV void pop();
      CUDA_HOSTDEV bool push_back(const vmesh::GlobalID& globalID);
      bool push_back(const std::vector<vmesh::GlobalID>& blocks);
      CUDA_HOSTDEV bool push_back(const split::SplitVector<vmesh::GlobalID>& blocks);
      CUDA_DEV bool replaceBlock(const vmesh::GlobalID& globalID,const vmesh::GlobalID& localID);
      void setGrid();
      bool setGrid(const std::vector<vmesh::GlobalID>& globalIDs);
      bool setGrid(const split::SplitVector<vmesh::GlobalID>& globalIDs);
      bool setMesh(const size_t& meshID);
      void setNewSize(const vmesh::LocalID& newSize);
      CUDA_HOSTDEV size_t size() const;
      CUDA_HOSTDEV size_t sizeInBytes() const;
      CUDA_HOSTDEV void swap(VelocityMesh& vm);

      void dev_prefetchHost();
      void dev_prefetchDevice();

   private:
      size_t meshID;

      //std::vector<vmesh::GlobalID> localToGlobalMap;
      //OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID> globalToLocalMap;
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *globalToLocalMap;
      split::SplitVector<vmesh::GlobalID> *localToGlobalMap;
   };

   // ***** DEFINITIONS OF MEMBER FUNCTIONS ***** //


   inline VelocityMesh::VelocityMesh() {
      meshID = std::numeric_limits<size_t>::max();
      // Set sizepower to 10 (1024 blocks) straight away so there's enough room to grow
      globalToLocalMap = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(10);
      localToGlobalMap = new split::SplitVector<vmesh::GlobalID>(1);
      localToGlobalMap->clear();
   }

   inline VelocityMesh::~VelocityMesh() {
      delete globalToLocalMap;
      delete localToGlobalMap;
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
      return *this;
   }

   CUDA_HOSTDEV inline size_t VelocityMesh::capacityInBytes() const {
      return localToGlobalMap->capacity()*sizeof(vmesh::GlobalID)
           + globalToLocalMap->bucket_count()*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
   }

   CUDA_HOSTDEV inline bool VelocityMesh::check() const {
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

   CUDA_HOSTDEV inline bool VelocityMesh::copy(const vmesh::LocalID& sourceLID,const vmesh::LocalID& targetLID) {
      const vmesh::GlobalID sourceGID = localToGlobalMap->at(sourceLID); // block at the end of list
      const vmesh::GlobalID targetGID = localToGlobalMap->at(targetLID); // removed block

      // at-function will throw out_of_range exception for non-existing global ID:
      #ifdef __CUDA_ARCH__
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

   CUDA_HOSTDEV inline size_t VelocityMesh::count(const vmesh::GlobalID& globalID) const {
      #ifdef __CUDA_ARCH__
      return globalToLocalMap->device_count(globalID);
      #else
      return globalToLocalMap->count(globalID);
      #endif
   }

   CUDA_HOSTDEV inline vmesh::GlobalID VelocityMesh::findBlock(vmesh::GlobalID cellIndices[3]) const {
      // Calculate i/j/k indices of the block that would own the cell:
      vmesh::GlobalID i_block = cellIndices[0] / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockLength[0];
      vmesh::GlobalID j_block = cellIndices[1] / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockLength[1];
      vmesh::GlobalID k_block = cellIndices[2] / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockLength[2];

      // Calculate block global ID:
      vmesh::GlobalID blockGID = getGlobalID(i_block,j_block,k_block);

      // If the block exists, return it:
      if (globalToLocalMap->find(blockGID) != globalToLocalMap->end()) {
         return blockGID;
      } else {
         return invalidGlobalID();
      }
   }

   CUDA_HOSTDEV inline bool VelocityMesh::getBlockCoordinates(const vmesh::GlobalID& globalID,Real coords[3]) const {
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

      coords[0] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[0] + indices[0]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[0];
      coords[1] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[1] + indices[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[1];
      coords[2] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[2] + indices[2]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[2];
      return true;
   }

   CUDA_HOSTDEV inline void VelocityMesh::getBlockInfo(const vmesh::GlobalID& globalID,Real* array) const {
      #ifndef NDEBUG
      if (globalID == invalidGlobalID()) {
         for (int i=0; i<6; ++i) array[i] = std::numeric_limits<Real>::infinity();
      }
      #endif

      vmesh::LocalID indices[3];
      indices[0] = globalID % (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0];
      indices[1] = (globalID / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]) % (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1];
      indices[2] = globalID / ((*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] * (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]);

      // Indices 0-2 contain coordinates of the lower left corner.
      // The values are the same as if getBlockCoordinates(globalID,&(array[0])) was called
      array[0] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[0] + indices[0]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[0];
      array[1] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[1] + indices[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[1];
      array[2] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[2] + indices[2]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[2];

      // Indices 3-5 contain the cell size.
      // The values are the same as if getCellSize(globalID,&(array[3])) was called
      array[3] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[0];
      array[4] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[1];
      array[5] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[2];
   }

   CUDA_HOSTDEV inline const Real* VelocityMesh::getBlockSize() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize;
   }

   CUDA_HOSTDEV inline bool VelocityMesh::getBlockSize(const vmesh::GlobalID& globalID,Real size[3]) const {
      size[0] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[0];
      size[1] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[1];
      size[2] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[2];
      return true;
   }

   CUDA_HOSTDEV inline const Real* VelocityMesh::getCellSize(const uint8_t& refLevel) const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize;
   }

   CUDA_HOSTDEV inline bool VelocityMesh::getCellSize(const vmesh::GlobalID& globalID,Real size[3]) const {
      size[0] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[0];
      size[1] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[1];
      size[2] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[2];
      return true;
   }

   CUDA_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID& localID) const {
      #ifndef NDEBUG
      if (localID >= localToGlobalMap->size()) {
         printf("ERROR invalid local id %lu\n",localID);
         exit(1);
      }
      #endif

      return localToGlobalMap->at(localID);
   }

   CUDA_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(const Real* coords) const {
      if (coords[0] < (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[0] || coords[0] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMaxLimits[0] ||
         (coords[1] < (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[1] || coords[1] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMaxLimits[1] ||
          coords[2] < (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[2] || coords[2] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMaxLimits[2])) {
         return invalidGlobalID();
      }

      const vmesh::LocalID indices[3] = {
         static_cast<vmesh::LocalID>(floor((coords[0] - (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[0]) / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[0])),
         static_cast<vmesh::LocalID>(floor((coords[1] - (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[1]) / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[1])),
         static_cast<vmesh::LocalID>(floor((coords[2] - (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits[2]) / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[2]))
      };

      return indices[2]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]
              + indices[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] + indices[0];
   }

   CUDA_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(vmesh::LocalID indices[3]) const {
      if (indices[0] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]) return invalidGlobalID();
      if (indices[1] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]) return invalidGlobalID();
      if (indices[2] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[2]) return invalidGlobalID();
      return indices[2]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]
              + indices[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] + indices[0];
   }

   CUDA_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID& i,const vmesh::LocalID& j,const vmesh::LocalID& k) const {
      if (i >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] || j >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1] || k >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[2]) {
         return invalidGlobalID();
      }
      return i + j*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]
              + k*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1];
   }

   CUDA_HOSTDEV inline split::SplitVector<vmesh::GlobalID>& VelocityMesh::getGrid() {
      return *localToGlobalMap;
   }

   CUDA_HOSTDEV inline const vmesh::LocalID* VelocityMesh::getGridLength(const uint8_t& refLevel) const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength;
   }

   CUDA_HOSTDEV inline void VelocityMesh::getIndices(const vmesh::GlobalID& globalID,const uint8_t& refLevel,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const {
      getIndices(globalID,i,j,k);
   }

   CUDA_HOSTDEV inline void VelocityMesh::getIndices(const vmesh::GlobalID& globalID,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const {
      if (globalID >= invalidGlobalID()) {
         i = j = k = invalidBlockIndex();
      } else {
         i = globalID % (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0];
         j = (globalID / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]) % (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1];
         k = globalID / ((*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] * (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]);
      }
   }

   CUDA_HOSTDEV inline vmesh::LocalID VelocityMesh::getLocalID(const vmesh::GlobalID& globalID) const {
      #ifdef __CUDA_ARCH__
      auto it = globalToLocalMap->device_find(globalID);
      if (it != globalToLocalMap->device_end()) return it->second;
      #else
      auto it = globalToLocalMap->find(globalID);
      if (it != globalToLocalMap->end()) return it->second;
      #endif
      return invalidLocalID();
   }

   CUDA_HOSTDEV inline size_t VelocityMesh::getMesh() const {
      return meshID;
   }

   CUDA_HOSTDEV inline const Real* VelocityMesh::getMeshMaxLimits() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMaxLimits;
   }

   CUDA_HOSTDEV inline const Real* VelocityMesh::getMeshMinLimits() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits;
   }

   CUDA_HOSTDEV inline bool VelocityMesh::initialize(const size_t& meshID) {
      this->meshID = meshID;
      return true;
   }

   CUDA_HOSTDEV inline vmesh::LocalID VelocityMesh::invalidBlockIndex() {
      return INVALID_VEL_BLOCK_INDEX;
   }

   CUDA_HOSTDEV inline vmesh::GlobalID VelocityMesh::invalidGlobalID() {
      return INVALID_GLOBALID;
   }

   CUDA_HOSTDEV inline vmesh::LocalID VelocityMesh::invalidLocalID() {
      return INVALID_LOCALID;
   }

   CUDA_HOSTDEV inline bool VelocityMesh::isInitialized() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].initialized;
   }

   CUDA_HOSTDEV inline void VelocityMesh::pop() {
      if (size() == 0) return;

      const vmesh::LocalID lastLID = size()-1;
      const vmesh::GlobalID lastGID = localToGlobalMap->at(lastLID);
      #ifdef __CUDA_ARCH__
      auto last = globalToLocalMap->device_find(lastGID);
      globalToLocalMap->device_erase(last);
      #else
      auto last = globalToLocalMap->find(lastGID);
      globalToLocalMap->erase(last);
      #endif
      localToGlobalMap->pop_back();
   }

   CUDA_HOSTDEV inline bool VelocityMesh::push_back(const vmesh::GlobalID& globalID) {
      if (size() >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks) return false;
      if (globalID == invalidGlobalID()) return false;

      #ifdef __CUDA_ARCH__
      auto position
         = globalToLocalMap->device_insert(cuda::std::make_pair(globalID,localToGlobalMap->size()));
      if (position.second == true) {
         localToGlobalMap->device_push_back(globalID);
      }
      #else
      auto position
         = globalToLocalMap->insert(cuda::std::make_pair(globalID,localToGlobalMap->size()));
      if (position.second == true) {
         localToGlobalMap->push_back(globalID);
      }
      #endif
      return position.second;
   }

   inline bool VelocityMesh::push_back(const std::vector<vmesh::GlobalID>& blocks) {
      if (size()+blocks.size() > (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks) {
         printf("vmesh: too many blocks, current size is %lu",size());
         printf(", adding %lu blocks", blocks.size());
         printf(", max is %u\n",(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks);
         return false;
      }

      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap->insert(cuda::std::make_pair(blocks[b],localToGlobalMap->size()+b));
      }
      localToGlobalMap->insert(localToGlobalMap->end(),blocks.begin(),blocks.end());
      return true;
   }

   CUDA_HOSTDEV inline bool VelocityMesh::push_back(const split::SplitVector<vmesh::GlobalID>& blocks) {
      if (size()+blocks.size() > (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks) {
         printf("vmesh: too many blocks, current size is %lu",size());
         printf(", adding %lu blocks", blocks.size());
         printf(", max is %u\n",(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks);
         return false;
      }

      #ifdef __CUDA_ARCH__
      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap->device_insert(cuda::std::make_pair(blocks[b],localToGlobalMap->size()+b));
      }
      localToGlobalMap->device_insert(localToGlobalMap->end(),blocks.begin(),blocks.end());
      #else
      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap->insert(cuda::std::make_pair(blocks[b],localToGlobalMap->size()+b));
      }
      localToGlobalMap->insert(localToGlobalMap->end(),blocks.begin(),blocks.end());
      #endif
      return true;
   }

   CUDA_DEV inline bool VelocityMesh::replaceBlock(const vmesh::GlobalID& removeGID,const vmesh::GlobalID& addGID) {
      vmesh::LocalID removeLID = globalToLocalMap->read_element(removeGID);
      globalToLocalMap->device_erase(removeGID);
      globalToLocalMap->set_element(addGID,removeLID);
      localToGlobalMap->at(removeLID) = addGID;
   }

   inline void VelocityMesh::setGrid() {
      globalToLocalMap->clear();
      for (size_t i=0; i<localToGlobalMap->size(); ++i) {
         globalToLocalMap->insert(cuda::std::make_pair(localToGlobalMap->at(i),i));
      }
   }

   inline bool VelocityMesh::setGrid(const std::vector<vmesh::GlobalID>& globalIDs) {
      globalToLocalMap->clear();
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap->insert(cuda::std::make_pair(globalIDs[i],i));
      }
      localToGlobalMap->clear();
      localToGlobalMap->insert(localToGlobalMap->end(),globalIDs.begin(),globalIDs.end());
      return true;
   }
   inline bool VelocityMesh::setGrid(const split::SplitVector<vmesh::GlobalID>& globalIDs) {
      globalToLocalMap->clear();
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap->insert(cuda::std::make_pair(globalIDs[i],i));
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
      // CUDATODO: What's the point of this function? Garbage LIDs added?
      localToGlobalMap->resize(newSize);
   }

   CUDA_HOSTDEV inline size_t VelocityMesh::size() const {
      return localToGlobalMap->size();
   }

   CUDA_HOSTDEV inline size_t VelocityMesh::sizeInBytes() const {
      return globalToLocalMap->size()*sizeof(vmesh::GlobalID)
           + localToGlobalMap->size()*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
   }

   CUDA_HOSTDEV inline void VelocityMesh::swap(VelocityMesh& vm) {
      globalToLocalMap->swap(*(vm.globalToLocalMap));
      localToGlobalMap->swap(*(vm.localToGlobalMap));
   }

   inline void VelocityMesh::dev_prefetchHost() {
      if (localToGlobalMap->size() == 0) return;
      // In fact we only need to prefetch the buckets inside the hashmap to GPU, but use this call.
      //Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *uploaded = globalToLocalMap->upload(cuda_getStream());
      globalToLocalMap->optimizeCPU(cuda_getStream());
      localToGlobalMap->optimizeCPU(cuda_getStream());
      return;
   }

   inline void VelocityMesh::dev_prefetchDevice() {
      if (localToGlobalMap->size() == 0) return;
      // In fact we only need to prefetch the buckets inside the hashmap to GPU, but use this call.
      //Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *uploaded = globalToLocalMap->upload(cuda_getStream());
      globalToLocalMap->optimizeGPU(cuda_getStream());
      localToGlobalMap->optimizeGPU(cuda_getStream());
      return;
   }

} // namespace vmesh

#endif