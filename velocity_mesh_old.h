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
#include <signal.h>

#include "velocity_mesh_parameters.h"

namespace vmesh {

   // Open bucket power-of-two sized hash table with multiplicative fibonacci hashing
   template<typename GID, typename LID, int maxBucketOverflow = 8> class OpenBucketHashtable {
      private:
         int sizePower; // Logarithm (base two) of the size of the table
         size_t fill; // Number of filled buckets
         std::vector<std::pair<GID,LID>> buckets;

         // Fibonacci hash function for 32bit values
         uint32_t hash(GID in) const {
            in ^= in >> (32 - sizePower);
            uint32_t retval = (uint32_t)(in * 2654435769ul) >> (32 - sizePower);
            return retval;
         }
      public:
         OpenBucketHashtable() : sizePower(4),fill(0),buckets(1<<sizePower, std::pair<GID,LID>(INVALID_GLOBALID,INVALID_LOCALID)) {};
         OpenBucketHashtable(const OpenBucketHashtable<GID,LID>& other) : sizePower(other.sizePower),fill(other.fill),buckets(other.buckets) {};

         // Resize the table to fit more things. This is automatically invoked once
         // maxBucketOverflow has triggered.
         void rehash(int newSizePower) {
            if(newSizePower > 32) {
               throw std::out_of_range("OpenBucketHashtable ran into rehashing catastrophe and exceeded 32bit buckets.");
            }
            std::vector<std::pair<GID,LID>> newBuckets(1<<newSizePower, std::pair<GID,LID>(INVALID_LOCALID,INVALID_GLOBALID));
            sizePower = newSizePower;
            int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size

            // Iterate through all old elements and rehash them into the new array.
            for(auto& e : buckets) {
               // Skip empty buckets
               if(e.first == INVALID_LOCALID) {
                  continue;
               }

               uint32_t newHash = hash(e.first);
               bool found = false;
               for(int i=0; i<maxBucketOverflow; i++) {
                  std::pair<GID, LID>& candidate = newBuckets[(newHash + i)&bitMask];
                  if(candidate.first == INVALID_GLOBALID) {
                     // Found an empty bucket, assign that one.
                     candidate = e;
                     found = true;
                     break;
                  }
               }

               if(!found) {
                  // Having arrived here means that we unsuccessfully rehashed and
                  // are *still* overflowing our buckets. So we need to try again with a bigger one.
                  return rehash(newSizePower+1);
               }
            }

            // Replace our buckets with the new ones
            buckets = newBuckets;
         }

         // Element access (by reference). Nonexistent elements get created.
         LID& at(const GID& key) {
            int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size
            uint32_t hashIndex = hash(key);

            // Try to find the matching bucket.
            for(int i=0; i<maxBucketOverflow; i++) {
               std::pair<GID,LID>& candidate = buckets[(hashIndex+i)&bitMask];
               if(candidate.first == key) {
                  // Found a match, return that
                  return candidate.second;
               }
               if(candidate.first == INVALID_GLOBALID) {
                  // Found an empty bucket, assign and return that.
                  candidate.first = key;
                  fill++;
                  return candidate.second;
               }
            }

            // Not found, and we have no free slots to create a new one. So we need to rehash to a larger size.
            rehash(sizePower+1);
            return at(key); // Recursive tail call to try again with larger table.
         }


         // Typical array-like access with [] operator
         LID& operator[](const GID& key) {
            return at(key);
         }

         // For STL compatibility: size(), bucket_count(), count(GID), clear()
         size_t size() const {
            return fill;
         }

         size_t bucket_count() const {
            return buckets.size();
         }

         size_t count(const GID& key) const {
            if(find(key) != end()) {
               return 1;
            } else {
               return 0;
            } 
         }

         void clear() {
            buckets = std::vector<std::pair<GID,LID>>(1<<sizePower, {INVALID_GLOBALID,INVALID_LOCALID});
            fill = 0;
         }

         // Iterator type. Iterates through all non-empty buckets.
         class iterator : public std::iterator<std::random_access_iterator_tag, std::pair<GID,LID>> {
               OpenBucketHashtable<GID,LID>* hashtable;
               size_t index;
            public:
               iterator(OpenBucketHashtable<GID,LID>& hashtable, size_t index): hashtable(&hashtable),index(index) { }

               iterator& operator++() {
                  do {
                     index++;
                  } while(hashtable->buckets[index].first == INVALID_GLOBALID && index < hashtable->buckets.size());
                  return *this;
               }

               bool operator==(iterator other) const {
                  return &hashtable->buckets[index] == &other.hashtable->buckets[other.index];
               }
               bool operator!=(iterator other) const {
                  return &hashtable->buckets[index] != &other.hashtable->buckets[other.index];
               }
               std::pair<GID,LID>& operator*() const {
                  return hashtable->buckets[index];
               }
               std::pair<GID,LID>* operator->() const {
                  return &hashtable->buckets[index];
               }
               size_t getIndex() {
                  return index;
               }

         };


         // Const iterator.
         class const_iterator : public std::iterator<std::random_access_iterator_tag, std::pair<GID,LID>> {
               const OpenBucketHashtable<GID,LID>* hashtable;
               size_t index;
            public:
               explicit const_iterator(const OpenBucketHashtable<GID,LID>& hashtable, size_t index): hashtable(&hashtable),index(index) { }

               const_iterator& operator++() {
                  do {
                     index++;
                  } while(hashtable->buckets[index].first == INVALID_GLOBALID && index < hashtable->buckets.size());
                  return *this;
               }

               bool operator==(const_iterator other) const {
                  return &hashtable->buckets[index] == &other.hashtable->buckets[other.index];
               }
               bool operator!=(const_iterator other) const {
                  return &hashtable->buckets[index] != &other.hashtable->buckets[other.index];
               }
               const std::pair<GID,LID>& operator*() const {
                  return hashtable->buckets[index];
               }
               const std::pair<GID,LID>* operator->() const {
                  return &hashtable->buckets[index];
               }
               size_t getIndex() {
                  return index;
               }
         };

         iterator begin() {
            for(size_t i=0; i<buckets.size(); i++) {
               if(buckets[i].first != INVALID_GLOBALID) {
                  return iterator(*this, i);
               }
            }
            return end();
         }
         const_iterator begin() const {
            for(size_t i=0; i<buckets.size(); i++) {
               if(buckets[i].first != INVALID_GLOBALID) {
                  return const_iterator(*this, i);
               }
            }
            return end();
         }

         iterator end() {
            return iterator(*this, buckets.size());
         }
         const_iterator end() const {
            return const_iterator(*this, buckets.size());
         }


         // Element access by iterator
         iterator find(GID key) {
            int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size
            uint32_t hashIndex = hash(key);

            // Try to find the matching bucket.
            for(int i=0; i<maxBucketOverflow; i++) {
               const std::pair<GID,LID>& candidate = buckets[(hashIndex+i)&bitMask];
               if(candidate.first == key) {
                  // Found a match, return that
                  return iterator(*this, (hashIndex+i)&bitMask);
               }

               if(candidate.first == INVALID_GLOBALID) {
                  // Found an empty bucket. Return empty.
                  return end();
               }
            }

            // Not found
            return end();
         }

         const const_iterator find(GID key) const {
            int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size
            uint32_t hashIndex = hash(key);

            // Try to find the matching bucket.
            for(int i=0; i<maxBucketOverflow; i++) {
               const std::pair<GID,LID>& candidate = buckets[(hashIndex+i)&bitMask];
               if(candidate.first == key) {
                  // Found a match, return that
                  return const_iterator(*this, (hashIndex+i)&bitMask);
               }

               if(candidate.first == INVALID_GLOBALID) {
                  // Found an empty bucket. Return empty.
                  return end();
               }
            }

            // Not found
            return end();

         }

         // More STL compatibility implementations
         std::pair<iterator,bool> insert(std::pair<GID,LID> newEntry) {
            bool found = find(newEntry.first) != end();
            if(!found) {
               at(newEntry.first) = newEntry.second;
            }
            return std::pair<iterator,bool>(find(newEntry.first),!found); 
         }

         // Remove one element from the hash table.
         iterator erase(iterator keyPos) {
            // Due to overflowing buckets, this might require moving quite a bit of stuff around.
            size_t index = keyPos.getIndex();

            if(buckets[index].first != INVALID_GLOBALID) {
               // Decrease fill count if this spot wasn't empty already
               fill--;
            }
            // Clear the element itself.
            buckets[index]=std::pair<GID,LID>(INVALID_GLOBALID,INVALID_LOCALID);

            int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size
            GID nextBucket = buckets[(index+1) & bitMask].first;
            if(nextBucket == INVALID_GLOBALID) {
               // Easy case: if the next bucket is empty, we are done.
               ++keyPos;
               return keyPos;
            } else {
               // Othrwise, we need to renumber.
               // TODO: This is potentially quite slow. Are there more
               // efficient or elegant ways to resolve this?
               rehash(sizePower);

               // Find the next bucket member at its potentially new location.
               return find(nextBucket);
            }
         }

         void swap(OpenBucketHashtable<GID,LID>& other) {
            buckets.swap(other.buckets);
            int tempSizePower = sizePower;
            sizePower = other.sizePower;
            other.sizePower = tempSizePower;

            size_t tempFill = fill;
            fill = other.fill;
            other.fill = tempFill;
         }

   };

   template<typename GID,typename LID>
   class VelocityMesh {
    public:      
      VelocityMesh();
      ~VelocityMesh();

      size_t capacityInBytes() const;
      bool check() const;
      void clear();
      bool coarsenAllowed(const GID& globalID) const;
      bool copy(const LID& sourceLocalID,const LID& targetLocalID);
      size_t count(const GID& globalID) const;
      GID findBlockDown(uint8_t& refLevel,GID cellIndices[3]) const;
      GID findBlock(uint8_t& refLevel,GID cellIndices[3]) const;
      bool getBlockCoordinates(const GID& globalID,Real coords[3]) const;
      void getBlockInfo(const GID& globalID,Real* array) const;
      const Real* getBlockSize(const uint8_t& refLevel) const;
      bool getBlockSize(const GID& globalID,Real size[3]) const;
      const Real* getCellSize(const uint8_t& refLevel) const;
      bool getCellSize(const GID& globalID,Real size[3]) const;
      void getChildren(const GID& globalID,std::vector<GID>& children) const;
//      void     getChildren(const GlobalID& globalID,std::vector<GlobalID>& children);
      GID getGlobalID(const LID& localID) const;
      GID getGlobalID(const uint8_t& refLevel,const Real* coords) const;
      GID getGlobalID(const uint8_t& refLevel,LID indices[3]) const;
      GID getGlobalID(const uint32_t& refLevel,const LID& i,const LID& j,const LID& k) const;
      GID getGlobalIndexOffset(const uint8_t& refLevel=0);
      std::vector<GID>& getGrid();
      const LID* getGridLength(const uint8_t& refLevel) const;
//      void     getNeighbors(const GlobalID& globalID,std::vector<GlobalID>& neighborIDs);
      void getIndices(const GID& globalID,uint8_t& refLevel,LID& i,LID& j,LID& k) const;
      size_t getMesh() const;
      LID getLocalID(const GID& globalID) const;
      uint8_t getMaxAllowedRefinementLevel() const;
      GID getMaxVelocityBlocks() const;
      const Real* getMeshMaxLimits() const;
      const Real* getMeshMinLimits() const;
      void getNeighborsAtSameLevel(const GID& globalID,std::vector<GID>& neighborIDs) const;
      void getNeighborsExistingAtOffset(const GID& globalID,const int& i,const int& j,
              const int& k,std::vector<LID>& neighborLIDs,int32_t& refLevelDifference) const;
      int getOctant(const GID& globalID) const;
      GID getParent(const GID& globalID) const;
      uint8_t getRefinementLevel(const GID& globalID) const;
//      void     getSiblingNeighbors(const GlobalID& globalID,std::vector<GlobalID>& nbrs);
//      void     getSiblings(const GlobalID& globalID,GlobalID siblings[8]);
      void getSiblings(const GID& globalID,std::vector<GID>& siblings) const;
      bool hasChildren(const GID& globalID) const;
      GID hasGrandParent(const GID& globalID) const;
      bool initialize(const size_t& meshID,std::vector<vmesh::MeshParameters>& meshParameters);
      bool initialize(const size_t& meshID);
      static LID invalidBlockIndex();
      static GID invalidGlobalID();
      static LID invalidLocalID();
      bool isInitialized() const;
      void pop();
      bool push_back(const GID& globalID);
      bool push_back(const std::vector<GID>& blocks);
      bool refine(const GID& globalID,std::set<GID>& erasedBlocks,std::map<GID,LID>& insertedBlocks);
      void setGrid();
      bool setGrid(const std::vector<GID>& globalIDs);
      bool setMesh(const size_t& meshID);
      void setNewSize(const LID& newSize);
      size_t size() const;
      size_t sizeInBytes() const;
      void swap(VelocityMesh& vm);

    private:
      static std::vector<vmesh::MeshParameters> meshParameters;
      size_t meshID;

      std::vector<GID> localToGlobalMap;
      OpenBucketHashtable<GID,LID> globalToLocalMap; //
      //std::unordered_map<GID,LID> globalToLocalMap;
   };

   // ***** INITIALIZERS FOR STATIC MEMBER VARIABLES ***** //
   template<typename GID,typename LID> std::vector<vmesh::MeshParameters> VelocityMesh<GID,LID>::meshParameters;
   
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
      std::vector<GID>().swap(localToGlobalMap);
      globalToLocalMap.clear();
   }
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::coarsenAllowed(const GID& globalID) const {
      return false;
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
   GID VelocityMesh<GID,LID>::findBlockDown(uint8_t& refLevel,GID cellIndices[3]) const {
      // Calculate i/j/k indices of the block that would own the cell:
      GID i_block = cellIndices[0] / meshParameters[meshID].blockLength[0];
      GID j_block = cellIndices[1] / meshParameters[meshID].blockLength[1];
      GID k_block = cellIndices[2] / meshParameters[meshID].blockLength[2];
      
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
   GID VelocityMesh<GID,LID>::findBlock(uint8_t& refLevel,GID cellIndices[3]) const {
      return findBlockDown(refLevel,cellIndices);
   }

/*
   template<typename GID,typename LID> inline
   const GID* VelocityMesh<GID,LID>::getBaseGridLength() {
      return gridLength;
   }

   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getBaseGridBlockSize() {
      return blockSize;
   }
   
   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getBaseGridCellSize() {
      return cellSize;
   }
*/
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getBlockCoordinates(const GID& globalID,Real coords[3]) const {
      if (globalID == invalidGlobalID()) {
         for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
         return false;
      }
      
      uint8_t refLevel;
      LID indices[3];
      getIndices(globalID,refLevel,indices[0],indices[1],indices[2]);
      if (indices[0] == invalidBlockIndex()) {
         for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
         return false;
      }

      coords[0] = meshParameters[meshID].meshMinLimits[0] + indices[0]*meshParameters[meshID].blockSize[0];
      coords[1] = meshParameters[meshID].meshMinLimits[1] + indices[1]*meshParameters[meshID].blockSize[1];
      coords[2] = meshParameters[meshID].meshMinLimits[2] + indices[2]*meshParameters[meshID].blockSize[2];
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
      indices[0] = globalID % meshParameters[meshID].gridLength[0];
      indices[1] = (globalID / meshParameters[meshID].gridLength[0]) % meshParameters[meshID].gridLength[1];
      indices[2] = globalID / (meshParameters[meshID].gridLength[0] * meshParameters[meshID].gridLength[1]);

      // Indices 0-2 contain coordinates of the lower left corner.
      // The values are the same as if getBlockCoordinates(globalID,&(array[0])) was called
      array[0] = meshParameters[meshID].meshMinLimits[0] + indices[0]*meshParameters[meshID].blockSize[0];
      array[1] = meshParameters[meshID].meshMinLimits[1] + indices[1]*meshParameters[meshID].blockSize[1];
      array[2] = meshParameters[meshID].meshMinLimits[2] + indices[2]*meshParameters[meshID].blockSize[2];

      // Indices 3-5 contain the cell size.
      // The values are the same as if getCellSize(globalID,&(array[3])) was called
      array[3] = meshParameters[meshID].cellSize[0];
      array[4] = meshParameters[meshID].cellSize[1];
      array[5] = meshParameters[meshID].cellSize[2];
   }

   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getBlockSize(const uint8_t& refLevel) const {
      return meshParameters[meshID].blockSize;
   }
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getBlockSize(const GID& globalID,Real size[3]) const {
      size[0] = meshParameters[meshID].blockSize[0];
      size[1] = meshParameters[meshID].blockSize[1];
      size[2] = meshParameters[meshID].blockSize[2];
      return true;
   }
   
   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getCellSize(const uint8_t& refLevel) const {
      return meshParameters[meshID].cellSize;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getCellSize(const GID& globalID,Real size[3]) const {
      size[0] = meshParameters[meshID].cellSize[0];
      size[1] = meshParameters[meshID].cellSize[1];
      size[2] = meshParameters[meshID].cellSize[2];
      return true;
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getChildren(const GID& globalID,std::vector<GID>& children) const {
      children.clear();
      return;
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
   GID VelocityMesh<GID,LID>::getGlobalID(const uint8_t& refLevel,const Real* coords) const {
      if (coords[0] < meshParameters[meshID].meshMinLimits[0] || coords[0] >= meshParameters[meshID].meshMaxLimits[0] ||
         (coords[1] < meshParameters[meshID].meshMinLimits[1] || coords[1] >= meshParameters[meshID].meshMaxLimits[1] ||
          coords[2] < meshParameters[meshID].meshMinLimits[2] || coords[2] >= meshParameters[meshID].meshMaxLimits[2])) {
         return invalidGlobalID();
      }

      const LID indices[3] = {
         static_cast<LID>(floor((coords[0] - meshParameters[meshID].meshMinLimits[0]) / meshParameters[meshID].blockSize[0])),
         static_cast<LID>(floor((coords[1] - meshParameters[meshID].meshMinLimits[1]) / meshParameters[meshID].blockSize[1])),
         static_cast<LID>(floor((coords[2] - meshParameters[meshID].meshMinLimits[2]) / meshParameters[meshID].blockSize[2]))
      };

      return indices[2]*meshParameters[meshID].gridLength[1]*meshParameters[meshID].gridLength[0] 
              + indices[1]*meshParameters[meshID].gridLength[0] + indices[0];
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const uint8_t& refLevel,LID indices[3]) const {
      if (indices[0] >= meshParameters[meshID].gridLength[0]) return invalidGlobalID();
      if (indices[1] >= meshParameters[meshID].gridLength[1]) return invalidGlobalID();
      if (indices[2] >= meshParameters[meshID].gridLength[2]) return invalidGlobalID();
      return indices[2]*meshParameters[meshID].gridLength[1]*meshParameters[meshID].gridLength[0] 
              + indices[1]*meshParameters[meshID].gridLength[0] + indices[0];
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const uint32_t& refLevel,const LID& i,const LID& j,const LID& k) const {
      if (i >= meshParameters[meshID].gridLength[0] || j >= meshParameters[meshID].gridLength[1] || k >= meshParameters[meshID].gridLength[2]) {
         return invalidGlobalID();
      }
      return i + j*meshParameters[meshID].gridLength[0] 
              + k*meshParameters[meshID].gridLength[0]*meshParameters[meshID].gridLength[1];
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalIndexOffset(const uint8_t& refLevel) {
      return 0;
   }
   
   template<typename GID,typename LID> inline
   std::vector<GID>& VelocityMesh<GID,LID>::getGrid() {
      return localToGlobalMap;
   }

   template<typename GID,typename LID> inline
   const LID* VelocityMesh<GID,LID>::getGridLength(const uint8_t& refLevel) const {
      return meshParameters[meshID].gridLength;
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getIndices(const GID& globalID,uint8_t& refLevel,LID& i,LID& j,LID& k) const {
      refLevel = 0;
      if (globalID >= invalidGlobalID()) {
         i = j = k = invalidBlockIndex();
      } else {
         i = globalID % meshParameters[meshID].gridLength[0];
         j = (globalID / meshParameters[meshID].gridLength[0]) % meshParameters[meshID].gridLength[1];
         k = globalID / (meshParameters[meshID].gridLength[0] * meshParameters[meshID].gridLength[1]);
      }
   }

   template<typename GID,typename LID> inline
   LID VelocityMesh<GID,LID>::getLocalID(const GID& globalID) const {
      auto it = globalToLocalMap.find(globalID);
      if (it != globalToLocalMap.end()) return it->second;
      return invalidLocalID();
   }
   
   template<typename GID,typename LID> inline
   uint8_t VelocityMesh<GID,LID>::getMaxAllowedRefinementLevel() const {
      return 0;
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getMaxVelocityBlocks() const {
      return meshParameters[meshID].max_velocity_blocks;
   }
   
   template<typename GID,typename LID> inline
   size_t VelocityMesh<GID,LID>::getMesh() const {
      return meshID;
   }
   
   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getMeshMaxLimits() const {
      return meshParameters[meshID].meshMaxLimits;
   }
   
   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getMeshMinLimits() const {
      return meshParameters[meshID].meshMinLimits;
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getNeighborsAtSameLevel(const GID& globalID,std::vector<GID>& neighborIDs) const {
      neighborIDs.resize(27);
      
      // Calculate block refinement level and indices
      uint8_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      
      // Calculate global IDs of all 27 blocks:
      const LID Nx_max = meshParameters[meshID].gridLength[0];
      const LID Ny_max = meshParameters[meshID].gridLength[1];
      const LID Nz_max = meshParameters[meshID].gridLength[2];
      
      int nbr = 0;
      for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) {
         if (i+i_off < Nx_max && (j+j_off < Ny_max && k+k_off < Nz_max)) neighborIDs[nbr] = getGlobalID(0,i+i_off,j+j_off,k+k_off);
         else neighborIDs[nbr] = invalidGlobalID();
         ++nbr;
      }
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getNeighborsExistingAtOffset(const GID& globalID,const int& i_off,const int& j_off,const int& k_off,std::vector<LID>& neighborLocalIDs,int32_t& refLevelDifference) const {
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
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      
      // Return the requested neighbor if it exists:
      GID nbrGlobalID = getGlobalID(0,i+i_off,j+j_off,k+k_off);
      if (nbrGlobalID == invalidGlobalID()) return;

      auto nbr = globalToLocalMap.find(nbrGlobalID);
      if (nbr != globalToLocalMap.end()) {
         neighborLocalIDs.push_back(nbr->second);
         refLevelDifference = 0;
         return;
      }
   }

   template<typename GID,typename LID> inline
   int VelocityMesh<GID,LID>::getOctant(const GID& globalID) const {
      // Calculate block indices and refinement level
      uint8_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      
      const int i_oct = i % 2;
      const int j_oct = j % 2;
      const int k_oct = k % 2;
      return k_oct*2*2 + j_oct*2 + i_oct;
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getParent(const GID& globalID) const {
      return globalID;
   }
   
   template<typename GID,typename LID> inline
   uint8_t VelocityMesh<GID,LID>::getRefinementLevel(const GID& globalID) const {
      return 0;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getSiblings(const GID& globalID,std::vector<GID>& siblings) const {
      uint8_t refLevel;
      LID i,j,k;
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
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::hasChildren(const GID& globalID) const {
      return false;
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::hasGrandParent(const GID& globalID) const {
      return invalidGlobalID();
   }
      
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::initialize(const size_t& meshID) {
      this->meshID = meshID;
      return true;
   }
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::initialize(const size_t& meshID,std::vector<vmesh::MeshParameters>& meshParameters) {
      meshParameters[meshID].initialized = false;

      meshParameters[meshID].meshMinLimits[0] = meshParameters[meshID].meshLimits[0];
      meshParameters[meshID].meshMinLimits[1] = meshParameters[meshID].meshLimits[2];
      meshParameters[meshID].meshMinLimits[2] = meshParameters[meshID].meshLimits[4];
      meshParameters[meshID].meshMaxLimits[0] = meshParameters[meshID].meshLimits[1];
      meshParameters[meshID].meshMaxLimits[1] = meshParameters[meshID].meshLimits[3];
      meshParameters[meshID].meshMaxLimits[2] = meshParameters[meshID].meshLimits[5];
      
      // Calculate derived mesh parameters:
      meshParameters[meshID].gridSize[0] = meshParameters[meshID].meshMaxLimits[0] - meshParameters[meshID].meshMinLimits[0];
      meshParameters[meshID].gridSize[1] = meshParameters[meshID].meshMaxLimits[1] - meshParameters[meshID].meshMinLimits[1];
      meshParameters[meshID].gridSize[2] = meshParameters[meshID].meshMaxLimits[2] - meshParameters[meshID].meshMinLimits[2];
      
      meshParameters[meshID].blockSize[0] = meshParameters[meshID].gridSize[0] / meshParameters[meshID].gridLength[0];
      meshParameters[meshID].blockSize[1] = meshParameters[meshID].gridSize[1] / meshParameters[meshID].gridLength[1];
      meshParameters[meshID].blockSize[2] = meshParameters[meshID].gridSize[2] / meshParameters[meshID].gridLength[2];
      
      meshParameters[meshID].cellSize[0] = meshParameters[meshID].blockSize[0] / meshParameters[meshID].blockLength[0];
      meshParameters[meshID].cellSize[1] = meshParameters[meshID].blockSize[1] / meshParameters[meshID].blockLength[1];
      meshParameters[meshID].cellSize[2] = meshParameters[meshID].blockSize[2] / meshParameters[meshID].blockLength[2];

      meshParameters[meshID].max_velocity_blocks 
              = meshParameters[meshID].gridLength[0]
              * meshParameters[meshID].gridLength[1]
              * meshParameters[meshID].gridLength[2];
      meshParameters[meshID].initialized = true;

      vmesh::VelocityMesh<GID,LID>::meshParameters = meshParameters;
      return meshParameters[meshID].initialized;
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
      return meshParameters[meshID].initialized;
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
      if (size() >= meshParameters[meshID].max_velocity_blocks) return false;
      if (globalID == invalidGlobalID()) return false;

      auto position
        = globalToLocalMap.insert(std::make_pair(globalID,localToGlobalMap.size()));

      if (position.second == true) {
         localToGlobalMap.push_back(globalID);
      }

      return position.second;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::push_back(const std::vector<GID>& blocks) {
      if (size()+blocks.size() > meshParameters[meshID].max_velocity_blocks) {
         std::cerr << "vmesh: too many blocks, current size is " << size();
         std::cerr << ", adding " << blocks.size() << " blocks";
         std::cerr << ", max is " << meshParameters[meshID].max_velocity_blocks << std::endl;
         return false;
      }
         
      for (size_t b=0; b<blocks.size(); ++b) {
         globalToLocalMap.insert(std::make_pair(blocks[b],localToGlobalMap.size()+b));
      }
      localToGlobalMap.insert(localToGlobalMap.end(),blocks.begin(),blocks.end());

      return true;
   }
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::refine(const GID& globalID,std::set<GID>& erasedBlocks,std::map<GID,LID>& insertedBlocks) {
      return false;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::setGrid() {
      globalToLocalMap.clear();
      for (size_t i=0; i<localToGlobalMap.size(); ++i) {
         globalToLocalMap.insert(std::make_pair(localToGlobalMap[i],i));
      }
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::setGrid(const std::vector<GID>& globalIDs) {
      globalToLocalMap.clear();
      for (LID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap.insert(std::make_pair(globalIDs[i],i));
      }
      localToGlobalMap = globalIDs;
      return true;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::setMesh(const size_t& meshID) {
      if (meshID >= meshParameters.size()) return false;
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
