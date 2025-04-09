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

#ifndef VELOCITY_MESH_CPU_H
#define VELOCITY_MESH_CPU_H

#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <cmath>

//#include "object_wrapper.h"
#include "../open_bucket_hashtable.h"
#include "../velocity_mesh_parameters.h"

#if defined(DEBUG_VLASIATOR) || defined(DEBUG_SPATIAL_CELL)
   #ifndef DEBUG_VMESH
   #define DEBUG_VMESH
   #endif
#endif

namespace vmesh {

   class VelocityMesh {
    public:
      VelocityMesh();
      ~VelocityMesh();
      VelocityMesh(const VelocityMesh& other);
      const VelocityMesh& operator=(const VelocityMesh& other);

      size_t capacityInBytes() const;
      bool check() const;
      void clear(bool shrink=false);
      void clearMap(const vmesh::LocalID& newSize);
      bool move(const vmesh::LocalID& sourceLocalID,const vmesh::LocalID& targetLocalID);
      size_t count(const vmesh::GlobalID& globalID) const;
      vmesh::GlobalID findBlock(vmesh::GlobalID cellIndices[3]) const;
      bool getBlockCoordinates(const vmesh::GlobalID& globalID,Real coords[3]) const;
      void getBlockInfo(const vmesh::GlobalID& globalID,Real* array) const;
      const Real* getBlockSize() const;
      bool getBlockSize(const vmesh::GlobalID& globalID,Real size[3]) const;
      const Real* getCellSize() const;
      bool getCellSize(const vmesh::GlobalID& globalID,Real size[3]) const;
      vmesh::GlobalID getGlobalID(const vmesh::LocalID& localID) const;
      vmesh::GlobalID getGlobalID(const Real* coords) const;
      vmesh::GlobalID getGlobalID(vmesh::LocalID indices[3]) const;
      vmesh::GlobalID getGlobalID(const vmesh::LocalID& i,const vmesh::LocalID& j,const vmesh::LocalID& k) const;
      vmesh::GlobalID getGlobalIndexOffset();
      std::vector<vmesh::GlobalID>* getGrid();
      const vmesh::LocalID* getGridLength() const;
//      void     getNeighbors(const GlobalID& globalID,std::vector<GlobalID>& neighborIDs);
      void getIndices(const vmesh::GlobalID& globalID,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const;
      size_t getMesh() const;
      vmesh::LocalID getLocalID(const vmesh::GlobalID& globalID) const;
      vmesh::GlobalID getMaxVelocityBlocks() const;
      const Real* getMeshMaxLimits() const;
      const Real* getMeshMinLimits() const;
      bool initialize(const size_t& meshID);
      static vmesh::LocalID invalidBlockIndex();
      static vmesh::GlobalID invalidGlobalID();
      static vmesh::LocalID invalidLocalID();
      bool isInitialized() const;
      void pop();
      bool push_back(const vmesh::GlobalID& globalID);
      bool push_back(const std::vector<vmesh::GlobalID>& blocks);
      void setGrid();
      bool setGrid(const std::vector<vmesh::GlobalID>& globalIDs);
      bool setMesh(const size_t& meshID);
      void setNewSize(const vmesh::LocalID& newSize);
      void setNewCapacity(const vmesh::LocalID& newCapacity);
      size_t size(bool dummy=0) const;
      size_t sizeInBytes() const;
      // void swap(VelocityMesh& vm);

    private:
      size_t meshID;

      std::vector<vmesh::GlobalID> localToGlobalMap;
      OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID> globalToLocalMap;
      //std::unordered_map<vmesh::GlobalID,vmesh::LocalID> globalToLocalMap;
   };

   // ***** DEFINITIONS OF TEMPLATE MEMBER FUNCTIONS ***** //

   inline VelocityMesh::VelocityMesh() {
      meshID = std::numeric_limits<size_t>::max();
      globalToLocalMap = OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID>();
      localToGlobalMap = std::vector<vmesh::GlobalID>(1);
      localToGlobalMap.clear();
   }

   inline VelocityMesh::~VelocityMesh() { }

   inline VelocityMesh::VelocityMesh(const VelocityMesh& other) {
      meshID = other.meshID;
      globalToLocalMap = OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID>(other.globalToLocalMap);
      if (other.localToGlobalMap.size() > 0) {
         localToGlobalMap = std::vector<vmesh::GlobalID>(other.localToGlobalMap);
      } else {
         localToGlobalMap = std::vector<vmesh::GlobalID>(1);
         localToGlobalMap.clear();
      }
   }

   inline const VelocityMesh& VelocityMesh::operator=(const VelocityMesh& other) {
      meshID = other.meshID;
      globalToLocalMap = OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID>(other.globalToLocalMap);
      if (other.localToGlobalMap.size() > 0) {
         localToGlobalMap = std::vector<vmesh::GlobalID>(other.localToGlobalMap);
      } else {
         localToGlobalMap = std::vector<vmesh::GlobalID>(1);
         localToGlobalMap.clear();
      }
      return *this;
   }

   inline size_t VelocityMesh::capacityInBytes() const {
      return localToGlobalMap.capacity()*sizeof(vmesh::GlobalID)
           + globalToLocalMap.bucket_count()*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
   }

   inline bool VelocityMesh::check() const {
      bool ok = true;

      if (localToGlobalMap.size() != globalToLocalMap.size()) {
         std::cerr << "VMO ERROR: sizes differ, " << localToGlobalMap.size() << " vs " << globalToLocalMap.size() << std::endl;
         ok = false;
         exit(1);
      }

      for (size_t b=0; b<size(); ++b) {
         const vmesh::LocalID globalID = localToGlobalMap.at(b);
         auto it = globalToLocalMap.find(globalID);
         const vmesh::GlobalID localID = it->second;
         if (localID != b) {
            ok = false;
            std::cerr << "VMO ERROR: localToGlobalMap[" << b << "] = " << globalID << " but ";
            std::cerr << "globalToLocalMap[" << globalID << "] = " << localID << std::endl;
            exit(1);
         }
      }

      return ok;
   }

   inline void VelocityMesh::clear(bool shrink) {
      if (shrink) {
         globalToLocalMap = OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID>();
         localToGlobalMap = std::vector<vmesh::GlobalID>(1);
         localToGlobalMap.clear();
      } else {
         globalToLocalMap.clear();
         localToGlobalMap.clear();
      }
   }
   inline void VelocityMesh::clearMap(const vmesh::LocalID& newSize) {
      globalToLocalMap.clear();
      //globalToLocalMap.reserve(newSize); //OpenBucketHashTable does not have a reserve function
   }

   inline bool VelocityMesh::move(const vmesh::LocalID& sourceLID,const vmesh::LocalID& targetLID) {
      const vmesh::GlobalID moveGID = localToGlobalMap.at(sourceLID); // block to move (at the end of list)
      const vmesh::GlobalID removeGID = localToGlobalMap.at(targetLID); // removed block
      #ifdef DEBUG_VMESH
      if (sourceLID != size()-1) {
         printf("Warning! Moving velocity mesh entry from position which is not last LID!\n");
      }
      #endif
      // at-function will throw out_of_range exception for non-existing global ID:
      globalToLocalMap.at(moveGID) = targetLID;
      globalToLocalMap.erase(removeGID);
      localToGlobalMap.at(targetLID) = moveGID;
      localToGlobalMap.pop_back();
      return true;
   }

   inline size_t VelocityMesh::count(const vmesh::GlobalID& globalID) const {
      return globalToLocalMap.count(globalID);
   }

   inline vmesh::GlobalID VelocityMesh::findBlock(vmesh::GlobalID cellIndices[3]) const {
      // Calculate i/j/k indices of the block that would own the cell:
      vmesh::GlobalID i_block = cellIndices[0] / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockLength[0];
      vmesh::GlobalID j_block = cellIndices[1] / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockLength[1];
      vmesh::GlobalID k_block = cellIndices[2] / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockLength[2];

      // Calculate block global ID:
      vmesh::GlobalID blockGID = getGlobalID(i_block,j_block,k_block);

      // If the block exists, return it:
      if (globalToLocalMap.find(blockGID) != globalToLocalMap.end()) {
         return blockGID;
      } else {
         return invalidGlobalID();
      }
   }

/*
   inline const vmesh::GlobalID* VelocityMesh::getBaseGridLength() {
      return gridLength;
   }

   inline const Real* VelocityMesh::getBaseGridBlockSize() {
      return blockSize;
   }

   inline const Real* VelocityMesh::getBaseGridCellSize() {
      return cellSize;
   }
*/
   inline bool VelocityMesh::getBlockCoordinates(const vmesh::GlobalID& globalID,Real coords[3]) const {
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

   inline void VelocityMesh::getBlockInfo(const vmesh::GlobalID& globalID,Real* array) const {
      #ifdef DEBUG_VMESH
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

   inline const Real* VelocityMesh::getBlockSize() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize;
   }

   inline bool VelocityMesh::getBlockSize(const vmesh::GlobalID& globalID,Real size[3]) const {
      size[0] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[0];
      size[1] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[1];
      size[2] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[2];
      return true;
   }

   inline const Real* VelocityMesh::getCellSize() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize;
   }

   inline bool VelocityMesh::getCellSize(const vmesh::GlobalID& globalID,Real size[3]) const {
      size[0] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[0];
      size[1] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[1];
      size[2] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[2];
      return true;
   }

   inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID& localID) const {
      #ifdef DEBUG_VMESH
      if (localID >= localToGlobalMap.size()) {
         std::cerr << "ERROR invalid local id" << std::endl; exit(1);
      }
      #endif

      return localToGlobalMap.at(localID);
   }

   inline vmesh::GlobalID VelocityMesh::getGlobalID(const Real* coords) const {
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

   inline vmesh::GlobalID VelocityMesh::getGlobalID(vmesh::LocalID indices[3]) const {
      if (indices[0] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]) {
         return invalidGlobalID();
      }
      if (indices[1] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]) {
         return invalidGlobalID();
      }
      if (indices[2] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[2]) {
         return invalidGlobalID();
      }
      return indices[2]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]
              + indices[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] + indices[0];
   }

   inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID& i,const vmesh::LocalID& j,const vmesh::LocalID& k) const {
      if (i >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] || j >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1] || k >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[2]) {
         return invalidGlobalID();
      }
      return i + j*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]
              + k*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1];
   }

   inline vmesh::GlobalID VelocityMesh::getGlobalIndexOffset() {
      return 0;
   }

   inline std::vector<vmesh::GlobalID>* VelocityMesh::getGrid() {
      return &localToGlobalMap;
   }

   inline const vmesh::LocalID* VelocityMesh::getGridLength() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength;
   }

   inline void VelocityMesh::getIndices(const vmesh::GlobalID& globalID,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const {
      if (globalID >= invalidGlobalID()) {
         i = j = k = invalidBlockIndex();
      } else {
         i = globalID % (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0];
         j = (globalID / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]) % (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1];
         k = globalID / ((*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] * (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]);
      }
   }

   inline vmesh::LocalID VelocityMesh::getLocalID(const vmesh::GlobalID& globalID) const {
      auto it = globalToLocalMap.find(globalID);
      if (it != globalToLocalMap.end()) {
         return it->second;
      }
      return invalidLocalID();
   }

   inline vmesh::GlobalID VelocityMesh::getMaxVelocityBlocks() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks;
   }

   inline size_t VelocityMesh::getMesh() const {
      return meshID;
   }

   inline const Real* VelocityMesh::getMeshMaxLimits() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMaxLimits;
   }

   inline const Real* VelocityMesh::getMeshMinLimits() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].meshMinLimits;
   }

   inline bool VelocityMesh::initialize(const size_t& meshID) {
      this->meshID = meshID;
      return true;
   }

   inline vmesh::LocalID VelocityMesh::invalidBlockIndex() {
      return INVALID_VEL_BLOCK_INDEX;
   }

   inline vmesh::GlobalID VelocityMesh::invalidGlobalID() {
      return INVALID_GLOBALID;
   }

   inline vmesh::LocalID VelocityMesh::invalidLocalID() {
      return INVALID_LOCALID;
   }

   inline bool VelocityMesh::isInitialized() const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].initialized;
   }

   inline void VelocityMesh::pop() {
      if (size() == 0) {
         return;
      }

      const vmesh::LocalID lastLID = size()-1;
      const vmesh::GlobalID lastGID = localToGlobalMap.at(lastLID);
      auto last = globalToLocalMap.find(lastGID);

      globalToLocalMap.erase(last);
      localToGlobalMap.pop_back();
   }

   inline bool VelocityMesh::push_back(const vmesh::GlobalID& globalID) {
      if (size() >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks) {
         return false;
      }
      if (globalID == invalidGlobalID()) {
         return false;
      }

      auto position
        = globalToLocalMap.insert(std::make_pair(globalID,localToGlobalMap.size()));

      if (position.second == true) {
         localToGlobalMap.push_back(globalID);
      }

      return position.second;
   }

   inline bool VelocityMesh::push_back(const std::vector<vmesh::GlobalID>& blocks) {
      if (size()+blocks.size() > (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks) {
         std::cerr << "vmesh: too many blocks, current size is " << size();
         std::cerr << ", adding " << blocks.size() << " blocks";
         std::cerr << ", max is " << (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks << std::endl;
         return false;
      }

      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap.insert(std::make_pair(blocks[b],localToGlobalMap.size()+b));
      }
      localToGlobalMap.insert(localToGlobalMap.end(),blocks.begin(),blocks.end());

      return true;
   }

   inline void VelocityMesh::setGrid() {
      globalToLocalMap.clear();
      for (size_t i=0; i<localToGlobalMap.size(); ++i) {
         globalToLocalMap.insert(std::make_pair(localToGlobalMap.at(i),i));
      }
   }

   inline bool VelocityMesh::setGrid(const std::vector<vmesh::GlobalID>& globalIDs) {
      globalToLocalMap.clear();
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap.insert(std::make_pair(globalIDs[i],i));
      }
      localToGlobalMap.clear();
      localToGlobalMap.insert(localToGlobalMap.end(),globalIDs.begin(),globalIDs.end());
      return true;
   }

   inline bool VelocityMesh::setMesh(const size_t& meshID) {
      if (meshID >= vmesh::getMeshWrapper()->velocityMeshes->size()) {
         return false;
      }
      this->meshID = meshID;
      return true;
   }

   inline void VelocityMesh::setNewSize(const vmesh::LocalID& newSize) {
      localToGlobalMap.resize(newSize);
   }

   // Used in initialization
   inline void VelocityMesh::setNewCapacity(const vmesh::LocalID& newCapacity) {
      localToGlobalMap.reserve(newCapacity);
   }

   inline size_t VelocityMesh::size(bool dummy) const {
      (void) dummy; // GPU mesh compatibility
      return localToGlobalMap.size();
   }

   inline size_t VelocityMesh::sizeInBytes() const {
      return globalToLocalMap.size()*sizeof(vmesh::GlobalID)
           + localToGlobalMap.size()*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
   }

   // inline void VelocityMesh::swap(VelocityMesh& vm) {
   //    globalToLocalMap.swap(vm.globalToLocalMap);
   //    localToGlobalMap.swap(vm.localToGlobalMap);
   // }

} // namespace vmesh

#endif
