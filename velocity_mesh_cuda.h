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

#include "object_wrapper.h"
//#include "open_bucket_hashtable.h"
#include "include/hashinator/hashinator.h"
#include "include/splitvector/splitvec.h"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

// Forward declaration
ObjectWrapper& getObjectWrapper();

namespace vmesh {

   template<typename GID,typename LID>
   class VelocityMesh {
    public:
      VelocityMesh();
      ~VelocityMesh();

      size_t capacityInBytes() const;
      bool check() const;
      void clear();
      bool copy(const LID& sourceLocalID,const LID& targetLocalID);
      size_t count(const GID& globalID) const;
      GID findBlock(GID cellIndices[3]) const;
      bool getBlockCoordinates(const GID& globalID,Real coords[3]) const;
      void getBlockInfo(const GID& globalID,Real* array) const;
      const Real* getBlockSize() const;
      bool getBlockSize(const GID& globalID,Real size[3]) const;
      const Real* getCellSize(const uint8_t& refLevel=0) const;
      bool getCellSize(const GID& globalID,Real size[3]) const;
      GID getGlobalID(const LID& localID) const;
      GID getGlobalID(const Real* coords) const;
      GID getGlobalID(LID indices[3]) const;
      GID getGlobalID(const LID& i,const LID& j,const LID& k) const;
      split::SplitVector<GID>& getGrid();
      const LID* getGridLength(const uint8_t& refLevel=0) const;
      void getIndices(const GID& globalID,const uint8_t& refLevel,LID& i,LID& j,LID& k) const;
      void getIndices(const GID& globalID,LID& i,LID& j,LID& k) const;
      size_t getMesh() const;
      LID getLocalID(const GID& globalID) const;
      const Real* getMeshMaxLimits() const;
      const Real* getMeshMinLimits() const;
      bool initialize(const size_t& meshID);
      static LID invalidBlockIndex();
      static GID invalidGlobalID();
      static LID invalidLocalID();
      bool isInitialized() const;
      void pop();
      bool push_back(const GID& globalID);
      bool push_back(const std::vector<GID>& blocks);
      bool push_back(const split::SplitVector<GID>& blocks);
      void setGrid();
      bool setGrid(const std::vector<GID>& globalIDs);
      bool setGrid(const split::SplitVector<GID>& globalIDs);
      bool setMesh(const size_t& meshID);
      void setNewSize(const LID& newSize);
      size_t size() const;
      size_t sizeInBytes() const;
      void swap(VelocityMesh& vm);

   private:
      size_t meshID;

      //std::vector<GID> localToGlobalMap;
      //OpenBucketHashtable<GID,LID> globalToLocalMap;
      Hashinator::Hashmap<GID,LID> globalToLocalMap = Hashinator::Hashmap<GID,LID>(10); // Start of with a bit larger size for sizePower?
      split::SplitVector<GID> localToGlobalMap;
   };

   // ***** DEFINITIONS OF TEMPLATE MEMBER FUNCTIONS ***** //

   template<typename GID,typename LID> inline
   VelocityMesh<GID,LID>::VelocityMesh() {
      meshID = std::numeric_limits<size_t>::max();
   }

   template<typename GID,typename LID> inline
   VelocityMesh<GID,LID>::~VelocityMesh() { }

   template<typename GID,typename LID> inline
   size_t VelocityMesh<GID,LID>::capacityInBytes() const {
      return localToGlobalMap.capacity()*sizeof(GID)
           + globalToLocalMap.bucket_count()*(sizeof(GID)+sizeof(LID));
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::check() const {
      bool ok = true;

      if (localToGlobalMap.size() != globalToLocalMap.size()) {
         std::cerr << "VMO ERROR: sizes differ, " << localToGlobalMap.size() << " vs " << globalToLocalMap.size() << std::endl;
         ok = false;
         exit(1);
      }

      for (size_t b=0; b<size(); ++b) {
         const LID globalID = localToGlobalMap[b];
         auto it = globalToLocalMap.find(globalID);
         const GID localID = it->second;
         if (localID != b) {
            ok = false;
            std::cerr << "VMO ERROR: localToGlobalMap[" << b << "] = " << globalID << " but ";
            std::cerr << "globalToLocalMap[" << globalID << "] = " << localID << std::endl;
            exit(1);
         }
      }

      return ok;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::clear() {
      split::SplitVector<GID>().swap(localToGlobalMap);
      globalToLocalMap.clear();
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::copy(const LID& sourceLID,const LID& targetLID) {
      const GID sourceGID = localToGlobalMap[sourceLID]; // block at the end of list
      const GID targetGID = localToGlobalMap[targetLID]; // removed block

      // at-function will throw out_of_range exception for non-existing global ID:
      globalToLocalMap.at(sourceGID) = targetLID;
      localToGlobalMap[targetLID]    = sourceGID;
      globalToLocalMap.at(targetGID) = sourceLID; // These are needed to make pop() work
      localToGlobalMap[sourceLID]    = targetGID;
      return true;
   }

   template<typename GID,typename LID> inline
   size_t VelocityMesh<GID,LID>::count(const GID& globalID) const {
      return globalToLocalMap.count(globalID);
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::findBlock(GID cellIndices[3]) const {
      // Calculate i/j/k indices of the block that would own the cell:
      GID i_block = cellIndices[0] / getObjectWrapper().velocityMeshes[meshID].blockLength[0];
      GID j_block = cellIndices[1] / getObjectWrapper().velocityMeshes[meshID].blockLength[1];
      GID k_block = cellIndices[2] / getObjectWrapper().velocityMeshes[meshID].blockLength[2];

      // Calculate block global ID:
      GID blockGID = getGlobalID(0,i_block,j_block,k_block);

      // If the block exists, return it:
      if (globalToLocalMap.find(blockGID) != globalToLocalMap.end()) {
         return blockGID;
      } else {
         return invalidGlobalID();
      }
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getBlockCoordinates(const GID& globalID,Real coords[3]) const {
      if (globalID == invalidGlobalID()) {
         for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
         return false;
      }

      LID indices[3];
      getIndices(globalID,indices[0],indices[1],indices[2]);
      if (indices[0] == invalidBlockIndex()) {
         for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
         return false;
      }

      coords[0] = getObjectWrapper().velocityMeshes[meshID].meshMinLimits[0] + indices[0]*getObjectWrapper().velocityMeshes[meshID].blockSize[0];
      coords[1] = getObjectWrapper().velocityMeshes[meshID].meshMinLimits[1] + indices[1]*getObjectWrapper().velocityMeshes[meshID].blockSize[1];
      coords[2] = getObjectWrapper().velocityMeshes[meshID].meshMinLimits[2] + indices[2]*getObjectWrapper().velocityMeshes[meshID].blockSize[2];
      return true;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getBlockInfo(const GID& globalID,Real* array) const {
      #ifndef NDEBUG
      if (globalID == invalidGlobalID()) {
         for (int i=0; i<6; ++i) array[i] = std::numeric_limits<Real>::infinity();
      }
      #endif

      LID indices[3];
      indices[0] = globalID % getObjectWrapper().velocityMeshes[meshID].gridLength[0];
      indices[1] = (globalID / getObjectWrapper().velocityMeshes[meshID].gridLength[0]) % getObjectWrapper().velocityMeshes[meshID].gridLength[1];
      indices[2] = globalID / (getObjectWrapper().velocityMeshes[meshID].gridLength[0] * getObjectWrapper().velocityMeshes[meshID].gridLength[1]);

      // Indices 0-2 contain coordinates of the lower left corner.
      // The values are the same as if getBlockCoordinates(globalID,&(array[0])) was called
      array[0] = getObjectWrapper().velocityMeshes[meshID].meshMinLimits[0] + indices[0]*getObjectWrapper().velocityMeshes[meshID].blockSize[0];
      array[1] = getObjectWrapper().velocityMeshes[meshID].meshMinLimits[1] + indices[1]*getObjectWrapper().velocityMeshes[meshID].blockSize[1];
      array[2] = getObjectWrapper().velocityMeshes[meshID].meshMinLimits[2] + indices[2]*getObjectWrapper().velocityMeshes[meshID].blockSize[2];

      // Indices 3-5 contain the cell size.
      // The values are the same as if getCellSize(globalID,&(array[3])) was called
      array[3] = getObjectWrapper().velocityMeshes[meshID].cellSize[0];
      array[4] = getObjectWrapper().velocityMeshes[meshID].cellSize[1];
      array[5] = getObjectWrapper().velocityMeshes[meshID].cellSize[2];
   }

   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getBlockSize() const {
      return getObjectWrapper().velocityMeshes[meshID].blockSize;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getBlockSize(const GID& globalID,Real size[3]) const {
      size[0] = getObjectWrapper().velocityMeshes[meshID].blockSize[0];
      size[1] = getObjectWrapper().velocityMeshes[meshID].blockSize[1];
      size[2] = getObjectWrapper().velocityMeshes[meshID].blockSize[2];
      return true;
   }

   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getCellSize(const uint8_t& refLevel) const {
      return getObjectWrapper().velocityMeshes[meshID].cellSize;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getCellSize(const GID& globalID,Real size[3]) const {
      size[0] = getObjectWrapper().velocityMeshes[meshID].cellSize[0];
      size[1] = getObjectWrapper().velocityMeshes[meshID].cellSize[1];
      size[2] = getObjectWrapper().velocityMeshes[meshID].cellSize[2];
      return true;
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const LID& localID) const {
      #ifndef NDEBUG
      if (localID >= localToGlobalMap.size()) {
         std::cerr << "ERROR invalid local id" << std::endl; exit(1);
      }
      #endif

      return localToGlobalMap[localID];
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const Real* coords) const {
      if (coords[0] < getObjectWrapper().velocityMeshes[meshID].meshMinLimits[0] || coords[0] >= getObjectWrapper().velocityMeshes[meshID].meshMaxLimits[0] ||
         (coords[1] < getObjectWrapper().velocityMeshes[meshID].meshMinLimits[1] || coords[1] >= getObjectWrapper().velocityMeshes[meshID].meshMaxLimits[1] ||
          coords[2] < getObjectWrapper().velocityMeshes[meshID].meshMinLimits[2] || coords[2] >= getObjectWrapper().velocityMeshes[meshID].meshMaxLimits[2])) {
         return invalidGlobalID();
      }

      const LID indices[3] = {
         static_cast<LID>(floor((coords[0] - getObjectWrapper().velocityMeshes[meshID].meshMinLimits[0]) / getObjectWrapper().velocityMeshes[meshID].blockSize[0])),
         static_cast<LID>(floor((coords[1] - getObjectWrapper().velocityMeshes[meshID].meshMinLimits[1]) / getObjectWrapper().velocityMeshes[meshID].blockSize[1])),
         static_cast<LID>(floor((coords[2] - getObjectWrapper().velocityMeshes[meshID].meshMinLimits[2]) / getObjectWrapper().velocityMeshes[meshID].blockSize[2]))
      };

      return indices[2]*getObjectWrapper().velocityMeshes[meshID].gridLength[1]*getObjectWrapper().velocityMeshes[meshID].gridLength[0]
              + indices[1]*getObjectWrapper().velocityMeshes[meshID].gridLength[0] + indices[0];
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(LID indices[3]) const {
      if (indices[0] >= getObjectWrapper().velocityMeshes[meshID].gridLength[0]) return invalidGlobalID();
      if (indices[1] >= getObjectWrapper().velocityMeshes[meshID].gridLength[1]) return invalidGlobalID();
      if (indices[2] >= getObjectWrapper().velocityMeshes[meshID].gridLength[2]) return invalidGlobalID();
      return indices[2]*getObjectWrapper().velocityMeshes[meshID].gridLength[1]*getObjectWrapper().velocityMeshes[meshID].gridLength[0]
              + indices[1]*getObjectWrapper().velocityMeshes[meshID].gridLength[0] + indices[0];
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const LID& i,const LID& j,const LID& k) const {
      if (i >= getObjectWrapper().velocityMeshes[meshID].gridLength[0] || j >= getObjectWrapper().velocityMeshes[meshID].gridLength[1] || k >= getObjectWrapper().velocityMeshes[meshID].gridLength[2]) {
         return invalidGlobalID();
      }
      return i + j*getObjectWrapper().velocityMeshes[meshID].gridLength[0]
              + k*getObjectWrapper().velocityMeshes[meshID].gridLength[0]*getObjectWrapper().velocityMeshes[meshID].gridLength[1];
   }

   template<typename GID,typename LID> inline
   split::SplitVector<GID>& VelocityMesh<GID,LID>::getGrid() {
      return localToGlobalMap;
   }

   template<typename GID,typename LID> inline
   const LID* VelocityMesh<GID,LID>::getGridLength(const uint8_t& refLevel) const {
      return getObjectWrapper().velocityMeshes[meshID].gridLength;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getIndices(const GID& globalID,const uint8_t& refLevel,LID& i,LID& j,LID& k) const {
      getIndices(globalID,i,j,k);
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getIndices(const GID& globalID,LID& i,LID& j,LID& k) const {
      if (globalID >= invalidGlobalID()) {
         i = j = k = invalidBlockIndex();
      } else {
         i = globalID % getObjectWrapper().velocityMeshes[meshID].gridLength[0];
         j = (globalID / getObjectWrapper().velocityMeshes[meshID].gridLength[0]) % getObjectWrapper().velocityMeshes[meshID].gridLength[1];
         k = globalID / (getObjectWrapper().velocityMeshes[meshID].gridLength[0] * getObjectWrapper().velocityMeshes[meshID].gridLength[1]);
      }
   }

   template<typename GID,typename LID> inline
   LID VelocityMesh<GID,LID>::getLocalID(const GID& globalID) const {
      auto it = globalToLocalMap.find(globalID);
      if (it != globalToLocalMap.end()) return it->second;
      return invalidLocalID();
   }

   template<typename GID,typename LID> inline
   size_t VelocityMesh<GID,LID>::getMesh() const {
      return meshID;
   }

   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getMeshMaxLimits() const {
      return getObjectWrapper().velocityMeshes[meshID].meshMaxLimits;
   }

   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getMeshMinLimits() const {
      return getObjectWrapper().velocityMeshes[meshID].meshMinLimits;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::initialize(const size_t& meshID) {
      this->meshID = meshID;
      return true;
   }

   template<typename GID,typename LID> inline
   LID VelocityMesh<GID,LID>::invalidBlockIndex() {
      return INVALID_VEL_BLOCK_INDEX;
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::invalidGlobalID() {
      return INVALID_GLOBALID;
   }

   template<typename GID,typename LID> inline
   LID VelocityMesh<GID,LID>::invalidLocalID() {
      return INVALID_LOCALID;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::isInitialized() const {
      return getObjectWrapper().velocityMeshes[meshID].initialized;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::pop() {
      if (size() == 0) return;

      const LID lastLID = size()-1;
      const GID lastGID = localToGlobalMap[lastLID];
      auto last = globalToLocalMap.find(lastGID);

      globalToLocalMap.erase(last);
      localToGlobalMap.pop_back();
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::push_back(const GID& globalID) {
      if (size() >= getObjectWrapper().velocityMeshes[meshID].max_velocity_blocks) return false;
      if (globalID == invalidGlobalID()) return false;

      auto position
         = globalToLocalMap.insert(cuda::std::make_pair(globalID,localToGlobalMap.size()));

      if (position.second == true) {
         localToGlobalMap.push_back(globalID);
      }
      return position.second;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::push_back(const std::vector<GID>& blocks) {
      if (size()+blocks.size() > getObjectWrapper().velocityMeshes[meshID].max_velocity_blocks) {
         std::cerr << "vmesh: too many blocks, current size is " << size();
         std::cerr << ", adding " << blocks.size() << " blocks";
         std::cerr << ", max is " << getObjectWrapper().velocityMeshes[meshID].max_velocity_blocks << std::endl;
         return false;
      }

      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap.insert(cuda::std::make_pair(blocks[b],localToGlobalMap.size()+b));
      }
      //localToGlobalMap.insert(localToGlobalMap.end(),blocks.begin(),blocks.end());
      for (size_t b=0; b<blocks.size(); ++b) {
         localToGlobalMap.push_back(blocks.at(b));
      }
      return true;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::push_back(const split::SplitVector<GID>& blocks) {
      if (size()+blocks.size() > getObjectWrapper().velocityMeshes[meshID].max_velocity_blocks) {
         std::cerr << "vmesh: too many blocks, current size is " << size();
         std::cerr << ", adding " << blocks.size() << " blocks";
         std::cerr << ", max is " << getObjectWrapper().velocityMeshes[meshID].max_velocity_blocks << std::endl;
         return false;
      }

      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap.insert(cuda::std::make_pair(blocks[b],localToGlobalMap.size()+b));
      }
      localToGlobalMap.insert(localToGlobalMap.end(),blocks.begin(),blocks.end());

      return true;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::setGrid() {
      globalToLocalMap.clear();
      for (size_t i=0; i<localToGlobalMap.size(); ++i) {
         globalToLocalMap.insert(cuda::std::make_pair(localToGlobalMap[i],i));
      }
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::setGrid(const std::vector<GID>& globalIDs) {
      globalToLocalMap.clear();
      for (LID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap.insert(cuda::std::make_pair(globalIDs[i],i));
      }
      localToGlobalMap = globalIDs;
      return true;
   }
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::setGrid(const split::SplitVector<GID>& globalIDs) {
      globalToLocalMap.clear();
      for (LID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap.insert(cuda::std::make_pair(globalIDs[i],i));
      }
      localToGlobalMap = globalIDs;
      return true;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::setMesh(const size_t& meshID) {
      if (meshID >= getObjectWrapper().velocityMeshes.size()) return false;
      this->meshID = meshID;
      return true;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::setNewSize(const LID& newSize) {
      localToGlobalMap.resize(newSize);
   }

   template<typename GID,typename LID> inline
   size_t VelocityMesh<GID,LID>::size() const {
      return localToGlobalMap.size();
   }

   template<typename GID,typename LID> inline
   size_t VelocityMesh<GID,LID>::sizeInBytes() const {
      return globalToLocalMap.size()*sizeof(GID)
           + localToGlobalMap.size()*(sizeof(GID)+sizeof(LID));
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::swap(VelocityMesh& vm) {
      globalToLocalMap.swap(vm.globalToLocalMap);
      localToGlobalMap.swap(vm.localToGlobalMap);
   }

} // namespace vmesh

#endif
