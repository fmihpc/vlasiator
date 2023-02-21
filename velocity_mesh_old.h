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

#ifndef VELOCITY_MESH_OLD_H
#define VELOCITY_MESH_OLD_H

#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <cmath>

//#include "object_wrapper.h"
#include "open_bucket_hashtable.h"
#include "velocity_mesh_parameters.h"

namespace vmesh {

   class VelocityMesh {
    public:      
      VelocityMesh();
      ~VelocityMesh();
      VelocityMesh(const VelocityMesh& other);
      const VelocityMesh& operator=(const VelocityMesh& other);

      size_t capacityInBytes() const;
      bool check() const;
      void clear();
      bool coarsenAllowed(const vmesh::GlobalID& globalID) const;
      bool copy(const vmesh::LocalID& sourceLocalID,const vmesh::LocalID& targetLocalID);
      size_t count(const vmesh::GlobalID& globalID) const;
      vmesh::GlobalID findBlockDown(uint8_t& refLevel,vmesh::GlobalID cellIndices[3]) const;
      vmesh::GlobalID findBlock(uint8_t& refLevel,vmesh::GlobalID cellIndices[3]) const;
      bool getBlockCoordinates(const vmesh::GlobalID& globalID,Real coords[3]) const;
      void getBlockInfo(const vmesh::GlobalID& globalID,Real* array) const;
      const Real* getBlockSize(const uint8_t& refLevel) const;
      bool getBlockSize(const vmesh::GlobalID& globalID,Real size[3]) const;
      const Real* getCellSize(const uint8_t& refLevel) const;
      bool getCellSize(const vmesh::GlobalID& globalID,Real size[3]) const;
      void getChildren(const vmesh::GlobalID& globalID,std::vector<vmesh::GlobalID>& children) const;
//      void     getChildren(const GlobalID& globalID,std::vector<GlobalID>& children);
      vmesh::GlobalID getGlobalID(const vmesh::LocalID& localID) const;
      vmesh::GlobalID getGlobalID(const uint8_t& refLevel,const Real* coords) const;
      vmesh::GlobalID getGlobalID(const uint8_t& refLevel,vmesh::LocalID indices[3]) const;
      vmesh::GlobalID getGlobalID(const uint32_t& refLevel,const vmesh::LocalID& i,const vmesh::LocalID& j,const vmesh::LocalID& k) const;
      vmesh::GlobalID getGlobalIndexOffset(const uint8_t& refLevel=0);
      std::vector<vmesh::GlobalID>& getGrid();
      const vmesh::LocalID* getGridLength(const uint8_t& refLevel) const;
//      void     getNeighbors(const GlobalID& globalID,std::vector<GlobalID>& neighborIDs);
      void getIndices(const vmesh::GlobalID& globalID,uint8_t& refLevel,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const;
      size_t getMesh() const;
      vmesh::LocalID getLocalID(const vmesh::GlobalID& globalID) const;
      uint8_t getMaxAllowedRefinementLevel() const;
      vmesh::GlobalID getMaxVelocityBlocks() const;
      const Real* getMeshMaxLimits() const;
      const Real* getMeshMinLimits() const;
      void getNeighborsAtSameLevel(const vmesh::GlobalID& globalID,std::vector<vmesh::GlobalID>& neighborIDs) const;
      void getNeighborsExistingAtOffset(const vmesh::GlobalID& globalID,const int& i,const int& j,
              const int& k,std::vector<vmesh::LocalID>& neighborLIDs,int32_t& refLevelDifference) const;
      int getOctant(const vmesh::GlobalID& globalID) const;
      vmesh::GlobalID getParent(const vmesh::GlobalID& globalID) const;
      uint8_t getRefinementLevel(const vmesh::GlobalID& globalID) const;
//      void     getSiblingNeighbors(const GlobalID& globalID,std::vector<GlobalID>& nbrs);
//      void     getSiblings(const GlobalID& globalID,GlobalID siblings[8]);
      void getSiblings(const vmesh::GlobalID& globalID,std::vector<vmesh::GlobalID>& siblings) const;
      bool hasChildren(const vmesh::GlobalID& globalID) const;
      vmesh::GlobalID hasGrandParent(const vmesh::GlobalID& globalID) const;
      bool initialize(const size_t& meshID);
      static vmesh::LocalID invalidBlockIndex();
      static vmesh::GlobalID invalidGlobalID();
      static vmesh::LocalID invalidLocalID();
      bool isInitialized() const;
      void pop();
      bool push_back(const vmesh::GlobalID& globalID);
      bool push_back(const std::vector<vmesh::GlobalID>& blocks);
      bool refine(const vmesh::GlobalID& globalID,std::set<vmesh::GlobalID>& erasedBlocks,std::map<vmesh::GlobalID,vmesh::LocalID>& insertedBlocks);
      void setGrid();
      bool setGrid(const std::vector<vmesh::GlobalID>& globalIDs);
      bool setMesh(const size_t& meshID);
      void setNewSize(const vmesh::LocalID& newSize);
      size_t size() const;
      size_t sizeInBytes() const;
      void swap(VelocityMesh& vm);

    private:
      size_t meshID;

      std::vector<vmesh::GlobalID> *localToGlobalMap;
      OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID> *globalToLocalMap;
      //std::unordered_map<vmesh::GlobalID,vmesh::LocalID> globalToLocalMap;
   };

   // ***** DEFINITIONS OF TEMPLATE MEMBER FUNCTIONS ***** //

   inline VelocityMesh::VelocityMesh() { 
      meshID = std::numeric_limits<size_t>::max();
      globalToLocalMap = new OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID>();
      localToGlobalMap = new std::vector<vmesh::GlobalID>(1);
      localToGlobalMap->clear();
   }
   
   inline VelocityMesh::~VelocityMesh() {
      delete globalToLocalMap;
      delete localToGlobalMap;
   }

   inline VelocityMesh::VelocityMesh(const VelocityMesh& other) {
      meshID = other.meshID;
      globalToLocalMap = new OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID>(*(other.globalToLocalMap));
      if (other.localToGlobalMap->size() > 0) {
         localToGlobalMap = new std::vector<vmesh::GlobalID>(*(other.localToGlobalMap));
      } else {
         localToGlobalMap = new std::vector<vmesh::GlobalID>(1);
         localToGlobalMap->clear();
      }
   }

   inline const VelocityMesh& VelocityMesh::operator=(const VelocityMesh& other) {
      delete globalToLocalMap;
      delete localToGlobalMap;
      meshID = other.meshID;
      globalToLocalMap = new OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID>(*(other.globalToLocalMap));
      if (other.localToGlobalMap->size() > 0) {
         localToGlobalMap = new std::vector<vmesh::GlobalID>(*(other.localToGlobalMap));
      } else {
         localToGlobalMap = new std::vector<vmesh::GlobalID>(1);
         localToGlobalMap->clear();
      }
      return *this;
   }
   
   inline size_t VelocityMesh::capacityInBytes() const {
      return localToGlobalMap->capacity()*sizeof(vmesh::GlobalID)
           + globalToLocalMap->bucket_count()*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
   }

   inline bool VelocityMesh::check() const {
      bool ok = true;

      if (localToGlobalMap->size() != globalToLocalMap->size()) {
         std::cerr << "VMO ERROR: sizes differ, " << localToGlobalMap->size() << " vs " << globalToLocalMap->size() << std::endl;
         ok = false;
         exit(1);	 
      }

      for (size_t b=0; b<size(); ++b) {
         const vmesh::LocalID globalID = localToGlobalMap->at(b);
         auto it = globalToLocalMap->find(globalID);
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
   
   inline void VelocityMesh::clear() {
      std::vector<vmesh::GlobalID>().swap(*localToGlobalMap);
      globalToLocalMap->clear();
   }
   
   inline bool VelocityMesh::coarsenAllowed(const vmesh::GlobalID& globalID) const {
      return false;
   }
   
   inline bool VelocityMesh::copy(const vmesh::LocalID& sourceLID,const vmesh::LocalID& targetLID) {
      const vmesh::GlobalID sourceGID = localToGlobalMap->at(sourceLID); // block at the end of list
      const vmesh::GlobalID targetGID = localToGlobalMap->at(targetLID); // removed block

      // at-function will throw out_of_range exception for non-existing global ID:
      globalToLocalMap->at(sourceGID) = targetLID;
      localToGlobalMap->at(targetLID) = sourceGID;
      globalToLocalMap->at(targetGID) = sourceLID; // These are needed to make pop() work
      localToGlobalMap->at(sourceLID) = targetGID;
      return true;
   }
   
   inline size_t VelocityMesh::count(const vmesh::GlobalID& globalID) const {
      return globalToLocalMap->count(globalID);
   }
   
   inline vmesh::GlobalID VelocityMesh::findBlockDown(uint8_t& refLevel,vmesh::GlobalID cellIndices[3]) const {
      // Calculate i/j/k indices of the block that would own the cell:
      vmesh::GlobalID i_block = cellIndices[0] / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockLength[0];
      vmesh::GlobalID j_block = cellIndices[1] / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockLength[1];
      vmesh::GlobalID k_block = cellIndices[2] / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockLength[2];
      
      // Calculate block global ID:
      vmesh::GlobalID blockGID = getGlobalID(0,i_block,j_block,k_block);
      
      // If the block exists, return it:
      if (globalToLocalMap->find(blockGID) != globalToLocalMap->end()) {
         return blockGID;
      } else {
         return invalidGlobalID();
      }
   }
    
   inline vmesh::GlobalID VelocityMesh::findBlock(uint8_t& refLevel,vmesh::GlobalID cellIndices[3]) const {
      return findBlockDown(refLevel,cellIndices);
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
      
      uint8_t refLevel;
      vmesh::LocalID indices[3];
      getIndices(globalID,refLevel,indices[0],indices[1],indices[2]);
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

   inline const Real* VelocityMesh::getBlockSize(const uint8_t& refLevel) const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize;
   }
   
   inline bool VelocityMesh::getBlockSize(const vmesh::GlobalID& globalID,Real size[3]) const {
      size[0] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[0];
      size[1] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[1];
      size[2] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].blockSize[2];
      return true;
   }
   
   inline const Real* VelocityMesh::getCellSize(const uint8_t& refLevel) const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize;
   }

   inline bool VelocityMesh::getCellSize(const vmesh::GlobalID& globalID,Real size[3]) const {
      size[0] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[0];
      size[1] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[1];
      size[2] = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].cellSize[2];
      return true;
   }
   
   inline void VelocityMesh::getChildren(const vmesh::GlobalID& globalID,std::vector<vmesh::GlobalID>& children) const {
      children.clear();
      return;
   }

   inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID& localID) const {
      #ifndef NDEBUG
      if (localID >= localToGlobalMap->size()) {
         std::cerr << "ERROR invalid local id" << std::endl; exit(1);
      }
      #endif
      
      return localToGlobalMap->at(localID);
   }

   inline vmesh::GlobalID VelocityMesh::getGlobalID(const uint8_t& refLevel,const Real* coords) const {
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
   
   inline vmesh::GlobalID VelocityMesh::getGlobalID(const uint8_t& refLevel,vmesh::LocalID indices[3]) const {
      if (indices[0] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]) return invalidGlobalID();
      if (indices[1] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]) return invalidGlobalID();
      if (indices[2] >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[2]) return invalidGlobalID();
      return indices[2]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] 
              + indices[1]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] + indices[0];
   }
   
   inline vmesh::GlobalID VelocityMesh::getGlobalID(const uint32_t& refLevel,const vmesh::LocalID& i,const vmesh::LocalID& j,const vmesh::LocalID& k) const {
      if (i >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] || j >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1] || k >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[2]) {
         return invalidGlobalID();
      }
      return i + j*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] 
              + k*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]*(*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1];
   }
   
   inline vmesh::GlobalID VelocityMesh::getGlobalIndexOffset(const uint8_t& refLevel) {
      return 0;
   }
   
   inline std::vector<vmesh::GlobalID>& VelocityMesh::getGrid() {
      return *localToGlobalMap;
   }

   inline const vmesh::LocalID* VelocityMesh::getGridLength(const uint8_t& refLevel) const {
      return (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength;
   }
   
   inline void VelocityMesh::getIndices(const vmesh::GlobalID& globalID,uint8_t& refLevel,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const {
      refLevel = 0;
      if (globalID >= invalidGlobalID()) {
         i = j = k = invalidBlockIndex();
      } else {
         i = globalID % (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0];
         j = (globalID / (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0]) % (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1];
         k = globalID / ((*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0] * (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1]);
      }
   }

   inline vmesh::LocalID VelocityMesh::getLocalID(const vmesh::GlobalID& globalID) const {
      auto it = globalToLocalMap->find(globalID);
      if (it != globalToLocalMap->end()) return it->second;
      return invalidLocalID();
   }
   
   inline uint8_t VelocityMesh::getMaxAllowedRefinementLevel() const {
      return 0;
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
   
   inline void VelocityMesh::getNeighborsAtSameLevel(const vmesh::GlobalID& globalID,std::vector<vmesh::GlobalID>& neighborIDs) const {
      neighborIDs.resize(27);
      
      // Calculate block refinement level and indices
      uint8_t refLevel;
      vmesh::LocalID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      
      // Calculate global IDs of all 27 blocks:
      const vmesh::LocalID Nx_max = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[0];
      const vmesh::LocalID Ny_max = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[1];
      const vmesh::LocalID Nz_max = (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].gridLength[2];
      
      int nbr = 0;
      for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) {
         if (i+i_off < Nx_max && (j+j_off < Ny_max && k+k_off < Nz_max)) neighborIDs[nbr] = getGlobalID(0,i+i_off,j+j_off,k+k_off);
         else neighborIDs[nbr] = invalidGlobalID();
         ++nbr;
      }
   }

   inline void VelocityMesh::getNeighborsExistingAtOffset(const vmesh::GlobalID& globalID,const int& i_off,const int& j_off,const int& k_off,std::vector<vmesh::LocalID>& neighborLocalIDs,int32_t& refLevelDifference) const {
      #ifndef NDEBUG
      if (abs(i_off) > 1 || (abs(j_off) > 1 || abs(k_off) > 1)) {
         std::stringstream ss;
         ss << "VelocityMesh ERROR: invalid offsets in getNeighborsExistingAtOffset " << i_off << ' ' << j_off << ' ' << k_off << std::endl;
         std::cerr << ss.str();
         exit(1);
      }
      #endif
      
      refLevelDifference = 0;
      neighborLocalIDs.clear();
      
      // Calculate block refinement level and indices
      uint8_t refLevel;
      vmesh::LocalID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      
      // Return the requested neighbor if it exists:
      vmesh::GlobalID nbrGlobalID = getGlobalID(0,i+i_off,j+j_off,k+k_off);
      if (nbrGlobalID == invalidGlobalID()) return;

      auto nbr = globalToLocalMap->find(nbrGlobalID);
      if (nbr != globalToLocalMap->end()) {
         neighborLocalIDs.push_back(nbr->second);
         refLevelDifference = 0;
         return;
      }
   }

   inline int VelocityMesh::getOctant(const vmesh::GlobalID& globalID) const {
      // Calculate block indices and refinement level
      uint8_t refLevel;
      vmesh::LocalID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      
      const int i_oct = i % 2;
      const int j_oct = j % 2;
      const int k_oct = k % 2;
      return k_oct*2*2 + j_oct*2 + i_oct;
   }

   inline vmesh::GlobalID VelocityMesh::getParent(const vmesh::GlobalID& globalID) const {
      return globalID;
   }
   
   inline uint8_t VelocityMesh::getRefinementLevel(const vmesh::GlobalID& globalID) const {
      return 0;
   }

   inline void VelocityMesh::getSiblings(const vmesh::GlobalID& globalID,std::vector<vmesh::GlobalID>& siblings) const {
      uint8_t refLevel;
      vmesh::LocalID i,j,k;
      getIndices(globalID,refLevel,i,j,k);

      siblings.resize(8);
      
      i -= (i % 2);
      j -= (j % 2);
      k -= (k % 2);
      
      siblings[0] = getGlobalID(refLevel,i  ,j  ,k  );
      siblings[1] = getGlobalID(refLevel,i+1,j  ,k  );
      siblings[2] = getGlobalID(refLevel,i  ,j+1,k  );
      siblings[3] = getGlobalID(refLevel,i+1,j+1,k  );
      siblings[4] = getGlobalID(refLevel,i  ,j  ,k+1);
      siblings[5] = getGlobalID(refLevel,i+1,j  ,k+1);
      siblings[6] = getGlobalID(refLevel,i  ,j+1,k+1);
      siblings[7] = getGlobalID(refLevel,i+1,j+1,k+1);
   }
   
   inline bool VelocityMesh::hasChildren(const vmesh::GlobalID& globalID) const {
      return false;
   }

   inline vmesh::GlobalID VelocityMesh::hasGrandParent(const vmesh::GlobalID& globalID) const {
      return invalidGlobalID();
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
      if (size() == 0) return;

      const vmesh::LocalID lastLID = size()-1;
      const vmesh::GlobalID lastGID = localToGlobalMap->at(lastLID);
      auto last = globalToLocalMap->find(lastGID);

      globalToLocalMap->erase(last);
      localToGlobalMap->pop_back();
   }

   inline bool VelocityMesh::push_back(const vmesh::GlobalID& globalID) {
      if (size() >= (*vmesh::getMeshWrapper()->velocityMeshes)[meshID].max_velocity_blocks) return false;
      if (globalID == invalidGlobalID()) return false;

      auto position
        = globalToLocalMap->insert(std::make_pair(globalID,localToGlobalMap->size()));

      if (position.second == true) {
         localToGlobalMap->push_back(globalID);
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
         globalToLocalMap->insert(std::make_pair(blocks[b],localToGlobalMap->size()+b));
      }
      localToGlobalMap->insert(localToGlobalMap->end(),blocks.begin(),blocks.end());

      return true;
   }
   
   inline bool VelocityMesh::refine(const vmesh::GlobalID& globalID,std::set<vmesh::GlobalID>& erasedBlocks,std::map<vmesh::GlobalID,vmesh::LocalID>& insertedBlocks) {
      return false;
   }

   inline void VelocityMesh::setGrid() {
      globalToLocalMap->clear();
      for (size_t i=0; i<localToGlobalMap->size(); ++i) {
         globalToLocalMap->insert(std::make_pair(localToGlobalMap->at(i),i));
      }
   }

   inline bool VelocityMesh::setGrid(const std::vector<vmesh::GlobalID>& globalIDs) {
      globalToLocalMap->clear();
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap->insert(std::make_pair(globalIDs[i],i));
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
      localToGlobalMap->resize(newSize);
   }

   inline size_t VelocityMesh::size() const {
      return localToGlobalMap->size();
   }
   
   inline size_t VelocityMesh::sizeInBytes() const {
      return globalToLocalMap->size()*sizeof(vmesh::GlobalID)
           + localToGlobalMap->size()*(sizeof(vmesh::GlobalID)+sizeof(vmesh::LocalID));
   }

   inline void VelocityMesh::swap(VelocityMesh& vm) {
      globalToLocalMap->swap(*(vm.globalToLocalMap));
      localToGlobalMap->swap(*(vm.localToGlobalMap));
   }
   
} // namespace vmesh

#endif
