/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2014 Finnish Meteorological Institute
 */

#ifndef VELOCITY_MESH_AMR_H
#define VELOCITY_MESH_AMR_H

#include <stdint.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <cmath>
#include <algorithm>

#ifndef NDEBUG
   #define DEBUG_AMR_MESH
#endif

namespace vmesh {

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
      GID findBlockDown(uint8_t& refLevel,GID cellIndices[3]);
      GID findBlock(uint8_t& refLevel,GID cellIndices[3]);
      static bool getBlockCoordinates(const GID& globalID,Real coords[3]);
      static void getBlockInfo(const GID& globalID,Real* array);
      static const Real* getBlockSize(const uint8_t& refLevel);
      static bool getBlockSize(const GID& globalID,Real size[3]);
      static const Real* getCellSize(const uint8_t& refLevel);
      static bool getCellSize(const GID& globalID,Real size[3]);
      static void getChildren(const GID& globalID,std::vector<GID>& children);
      GID getGlobalID(const LID& localID) const;
      static GID getGlobalID(const uint8_t& refLevel,const Real* coords);
      static GID getGlobalID(const uint8_t& refLevel,LID indices[3]);
      static GID getGlobalID(const uint32_t& refLevel,const LID& i,const LID& j,const LID& k);
      static GID getGlobalIndexOffset(const uint8_t& refLevel);
      static const LID* getGridLength(const uint8_t& refLevel);
      static void getIndices(const GID& globalID,uint8_t& refLevel,LID& i,LID& j,LID& k);
      LID getLocalID(const GID& globalID) const;
      static uint8_t getMaxAllowedRefinementLevel();
      static GID getMaxVelocityBlocks();
      static const Real* getMeshMaxLimits();
      static const Real* getMeshMinLimits();
      static void getNeighborsAtSameLevel(const GID& globalID,std::vector<GID>& neighborIDs);
      static void getNeighborsExistingAtSameLevel(const GID& globalID,std::vector<GID>& neighborIDs);
      void getNeighborsExistingAtOffset(const GID& globalID,const int& i,const int& j,const int& k,std::vector<LID>& neighborLIDs,int32_t& refLevelDifference) const;
      static int getOctant(const GID& globalID);
      static GID getParent(const GID& globalID);
      static uint8_t getRefinementLevel(const GID& globalID);
      static void getSiblingNeighbors(const GID& globalID,std::vector<GID>& nbrs);
      static void getSiblings(const GID& globalID,std::vector<GID>& siblings);
      bool hasChildren(const GID& globalID) const;
      GID hasGrandParent(const GID& globalID) const;
      static bool initialize(Real meshLimits[6],LID gridLength[3],LID blockLength[3],uint8_t refLevelMaxAllowed,
                             LID maxVelocityBlocks=std::numeric_limits<LID>::max());
      static LID invalidBlockIndex();
      static GID invalidGlobalID();
      static LID invalidLocalID();
      void pop();
      bool push_back(const GID& globalID);
      uint8_t push_back(const std::vector<vmesh::GlobalID>& blocks);
      bool refine(const GID& globalID,std::set<GID>& erasedBlocks,std::map<GID,LID>& insertedBlocks);
      bool setGrid(const std::vector<GID>& globalIDs);
      size_t size() const;
      size_t sizeInBytes() const;
      void swap(VelocityMesh& vm);

    private:
      static LID max_velocity_blocks;                                           /**< Maximum valid block local ID.*/
      static LID blockLength[3];                                                /**< Number of cells in a block per coordinate.*/
      static Real blockSize[3];                                                 /**< Size of a block at base grid level.*/
      static Real cellSize[3];                                                  /**< Size of a cell in a block at base grid level.*/
      static Real gridSize[3];                                                  /**< Size of the grid.*/
      static LID gridLength[3];                                                 /**< Number of blocks in the grid.*/
      static bool initialized;
      static Real meshMinLimits[3];                                             /**< Minimum coordinate values of the grid bounding box.*/
      static Real meshMaxLimits[3];                                             /**< Maximum coordinate values of the grid bounding box.*/
      static uint8_t refLevelMaxAllowed;      
      static std::vector<GID> offsets;                                          /**< Block global ID offsets for each refinement level.*/
      static Real* blockSizes;
      static Real* cellSizes;
      static LID* gridLengths;
      static size_t referenceCount;

      static Real* blockSizes;
      static Real* cellSizes;
      static LID* gridLengths;
      static size_t referenceCount;

      std::vector<GID> localToGlobalMap;
      std::unordered_map<GID,LID> globalToLocalMap;
      
      bool checkChildren(const GID& globalID) const;
      bool checkParent(const GID& globalID) const;
   };

   // ***** INITIALIZERS FOR STATIC MEMBER VARIABLES ***** //

   template<typename GID,typename LID> LID VelocityMesh<GID,LID>::max_velocity_blocks = 0;
   template<typename GID,typename LID> LID VelocityMesh<GID,LID>::blockLength[3] = {0,0,0};
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::blockSize[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> Real* VelocityMesh<GID,LID>::blockSizes = NULL;
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::cellSize[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> Real* VelocityMesh<GID,LID>::cellSizes = NULL;
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::gridSize[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> LID VelocityMesh<GID,LID>::gridLength[3] = {0,0,0};
   template<typename GID,typename LID> LID* VelocityMesh<GID,LID>::gridLengths = NULL;
   template<typename GID,typename LID> bool VelocityMesh<GID,LID>::initialized = false;
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::meshMinLimits[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::meshMaxLimits[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> std::vector<GID> VelocityMesh<GID,LID>::offsets;
   template<typename GID,typename LID> uint8_t VelocityMesh<GID,LID>::refLevelMaxAllowed = 0;
   template<typename GID,typename LID> size_t VelocityMesh<GID,LID>::referenceCount = 0;
   
   // ***** DEFINITIONS OF TEMPLATE MEMBER FUNCTIONS ***** //
   
   template<typename GID,typename LID> inline
   VelocityMesh<GID,LID>::VelocityMesh() { 
      ++VelocityMesh<GID,LID>::referenceCount;
   }

   template<typename GID,typename LID> inline
   VelocityMesh<GID,LID>::~VelocityMesh() { 
      // Reduce reference count:
      --VelocityMesh<GID,LID>::referenceCount;
      
      // If the last class was deleted, delete dynamically allocated arrays:
      if (VelocityMesh<GID,LID>::referenceCount == 0) {
	 delete [] VelocityMesh<GID,LID>::blockSizes;  VelocityMesh<GID,LID>::blockSizes = NULL;
	 delete [] VelocityMesh<GID,LID>::cellSizes;   VelocityMesh<GID,LID>::cellSizes = NULL;
	 delete [] VelocityMesh<GID,LID>::gridLengths; VelocityMesh<GID,LID>::gridLengths = NULL;
      }
   }

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
      }
      
      for (size_t b=0; b<size(); ++b) {
         const LID globalID = localToGlobalMap[b];
         typename std::unordered_map<GID,LID>::const_iterator it = globalToLocalMap.find(globalID);
         const GID localID = it->second;
         if (localID != b) {
            ok = false;
            std::cerr << "VMO ERROR: localToGlobalMap[" << b << "] = " << globalID << " but ";
            std::cerr << "globalToLocalMap[" << globalID << "] = " << localID << std::endl;
         }
         
         // Recursively check that none of the possible children exist:
         std::vector<GID> children;
         getChildren(globalID,children);
         bool hasChildren = false;
         for (size_t c=0; c<children.size(); ++c) {
            if (children[c] == invalidGlobalID()) continue;
            if (checkChildren(children[c]) == true) hasChildren = true;
         }
         if (hasChildren == true) {
            std::cerr << "VM AMR ERROR: block " << globalID << " at refinement level " << (int)getRefinementLevel(globalID);
            std::cerr << " exists and has children or grandchildren" << std::endl;
            std::cerr << "\t ERROR from " << __FILE__ << ':' << __LINE__ << std::endl;
            ok = false;
         }
         /*
          // Recursively check that parents do not exist:
          bool hasParents = false;
          if (checkParent(globalID) == true) hasParents=true;
          if (hasParents == true) {
          std::cerr << "VM AMR ERROR: block " << globalID << " at refinement level " << (int)getRefinementLevel(globalID);
          std::cerr << " exists and has a parent or a grandparent" << std::endl;
          std::cerr << "\t ERROR from " << __FILE__ << ':' << __LINE__ << std::endl;
          ok = false;
          }*/
      }
      return ok;
   }
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::checkChildren(const GID& globalID) const {
      // Check that none of the children exist:
      std::vector<GID> children;
      getChildren(globalID,children);
      for (size_t c=0; c<children.size(); ++c) {
         if (getLocalID(children[c]) != invalidLocalID()) {
            return true;
         }
      }
      
      // Recursively check that none of the grandchildren exist:
      for (size_t c=0; c<children.size(); ++c) {
         if (checkChildren(children[c]) == true) return true;
      }
      return false;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::checkParent(const GID& globalID) const {
      GID parent = getParent(globalID);
      if (parent == globalID) return false;
      if (getLocalID(parent) != invalidLocalID()) return true;
      return checkParent(parent);
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::clear() {
      std::vector<GID>().swap(localToGlobalMap);
      std::unordered_map<GID,LID>().swap(globalToLocalMap);
   }
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::coarsenAllowed(const GID& globalID) const {
      if (globalToLocalMap.count(globalID) == 0) return false;

      uint8_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      if (refLevel == 0) return false;

      // Coarsen not allowed if any of the siblings have children
      std::vector<GID> siblings;
      getSiblings(globalID,siblings);
      for (size_t s=0; s<siblings.size(); ++s) {
         std::vector<GID> children;
         getChildren(siblings[s],children);
         for (size_t c=0; c<children.size(); ++c) {
            if (globalToLocalMap.count(children[c]) > 0) return false;
         }
      }

      // Get all neighbors of the blocks in this octant
      std::vector<GID> octantNeighbors;
      getSiblingNeighbors(globalID,octantNeighbors);

      // Block cannot be coarsened if any of the octant neighbors have children
      for (size_t n=0; n<octantNeighbors.size(); ++n) {
         std::vector<GlobalID> children;
         getChildren(octantNeighbors[n],children);
         for (size_t c=0; c<children.size(); ++c) {
            if (globalToLocalMap.count(children[c]) != 0) return false;
         }
      }

      // Block can be coarsened
      return true;
   }
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::copy(const LID& sourceLID,const LID& targetLID) {
      #ifdef DEBUG_AMR_MESH
      bool ok=true;
      if (sourceLID >= localToGlobalMap.size()) ok=false;
      if (targetLID >= localToGlobalMap.size()) ok=false;
      if (ok == false) {
         std::cerr << "ERROR in copy, source/target LIDs are " << sourceLID << ' ' << targetLID << std::endl;
         std::cerr << "\t exiting from " << __FILE__ << ':' << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      const GID sourceGID = localToGlobalMap[sourceLID]; // block at the end of list
      const GID targetGID = localToGlobalMap[targetLID]; // removed block

      #ifdef DEBUG_AMR_MESH
      if (globalToLocalMap.find(sourceGID) == globalToLocalMap.end()) ok=false;
      if (globalToLocalMap.find(targetGID) == globalToLocalMap.end()) ok=false;
      if (ok == false) {
         std::cerr << "ERROR in copy, source/target LIDs are " << sourceLID << ' ' << targetLID << std::endl;
         std::cerr << "                             GIDs are " << sourceGID << ' ' << targetGID << std::endl;
         std::cerr << "\t sourceGID valid? ";
         if (globalToLocalMap.find(sourceGID) != globalToLocalMap.end()) std::cerr << "yes" << std::endl;
         else std::cerr << "no" << std::endl;
         std::cerr << "\t targetGID valid? ";
         if (globalToLocalMap.find(targetGID) != globalToLocalMap.end()) std::cerr << "yes" << std::endl;
         else std::cerr << "no" << std::endl;
         std::cerr << "\t exiting from " << __FILE__ << ':' << __LINE__ << std::endl;
         exit(1);
      }
      #endif
      
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
   GID VelocityMesh<GID,LID>::findBlockDown(uint8_t& refLevel,GID cellIndices[3]) {
      GID rvalue = invalidGlobalID();
      int refMul = std::pow(2,getMaxAllowedRefinementLevel()-refLevel);

      while (refLevel <= getMaxAllowedRefinementLevel()) {
         // Calculate i/j/k indices of the block that would own the cell at refinement level r:
         GID i_block = cellIndices[0] / (blockLength[0] * refMul);
         GID j_block = cellIndices[1] / (blockLength[1] * refMul);
         GID k_block = cellIndices[2] / (blockLength[2] * refMul);
         
         // Calculate block global ID:
         GID blockGID = getGlobalID(refLevel,i_block,j_block,k_block);
         
         // If the block exists, return it:
         if (globalToLocalMap.find(blockGID) != globalToLocalMap.end()) {
            return blockGID;
         }
         
         refMul *= 2;
         --refLevel;
      }

      refLevel = 0;
      return invalidGlobalID();
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::findBlock(uint8_t& refLevel,GID cellIndices[3]) {      
      GID rvalue = invalidGlobalID();
      if (refLevel > getMaxAllowedRefinementLevel()) return invalidGlobalID();

      int refMul = std::pow(2,getMaxAllowedRefinementLevel()-refLevel);
      for (uint8_t r=refLevel; r<=getMaxAllowedRefinementLevel(); ++r) {
         // Calculate i/j/k indices of the block that would own the cell at refinement level r:
         GID i_block = cellIndices[0] / (blockLength[0] * refMul);
         GID j_block = cellIndices[1] / (blockLength[1] * refMul);
         GID k_block = cellIndices[2] / (blockLength[2] * refMul);
         
         // Calculate block global ID:
         GID blockGID = getGlobalID(r,i_block,j_block,k_block);
         
         // If the block exists, return it:
         if (globalToLocalMap.find(blockGID) != globalToLocalMap.end()) {
            refLevel = r;
            return blockGID;
         }
         refMul /= 2;
      }
      return invalidGlobalID();
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getBlockCoordinates(const GID& globalID,Real coords[3]) {
      if (globalID == invalidGlobalID()) {
         for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
         return false;
      }
      
      // Calculate block's refinement level and (i,j,k) indices:
      uint8_t refLevel;
      LID indices[3];
      getIndices(globalID,refLevel,indices[0],indices[1],indices[2]);
      if (indices[0] == invalidBlockIndex()) {
         for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
         return false;
      }
      
      // NUmber of blocks per coordinate at this refinement level:
      const GlobalID multiplier = std::pow(2,refLevel);
      coords[0] = meshMinLimits[0] + indices[0]*blockSize[0]/multiplier;
      coords[1] = meshMinLimits[1] + indices[1]*blockSize[1]/multiplier;
      coords[2] = meshMinLimits[2] + indices[2]*blockSize[2]/multiplier;
      return true;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getBlockInfo(const GID& globalID,Real* array) {
      #ifdef DEBUG_AMR_MESH
      if (globalID == invalidGlobalID()) {
         for (int i=0; i<6; ++i) array[i] = std::numeric_limits<Real>::infinity();
      }
      #endif

      getBlockCoordinates(globalID,array+0);
      getCellSize(globalID,array+3);
   }

   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getBlockSize(const uint8_t& refLevel) {
      return blockSizes + refLevel*3;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getBlockSize(const GID& globalID,Real size[3]) {
      // Calculate block's refinement level and (i,j,k) indices:
      uint8_t refLevel;
      LID indices[3];
      getIndices(globalID,refLevel,indices[0],indices[1],indices[2]);

      // NUmber of blocks per coordinate at this refinement level:
      const GlobalID multiplier = std::pow(2,refLevel);
      
      // Calculate the number of blocks in each coordinate direction at this refinement level:
      size[0] = blockSize[0] / multiplier;
      size[1] = blockSize[1] / multiplier;
      size[2] = blockSize[2] / multiplier;
      return true;
   }

   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getCellSize(const uint8_t& refLevel) {
      return cellSizes + refLevel*3;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getCellSize(const GID& globalID,Real size[3]) {
      // Calculate block size and divide it by the number of cells in block
      getBlockSize(globalID,size);
      size[0] /= WID;
      size[1] /= WID;
      size[2] /= WID;
      return true;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getChildren(const GID& globalID,std::vector<GID>& children) {
      children.clear();

      uint8_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      if (refLevel == refLevelMaxAllowed) return;

      i *= 2;
      j *= 2;
      k *= 2;

      children.push_back(getGlobalID(refLevel+1,i  ,j  ,k  ));
      children.push_back(getGlobalID(refLevel+1,i+1,j  ,k  ));
      children.push_back(getGlobalID(refLevel+1,i  ,j+1,k  ));
      children.push_back(getGlobalID(refLevel+1,i+1,j+1,k  ));
      children.push_back(getGlobalID(refLevel+1,i  ,j  ,k+1));
      children.push_back(getGlobalID(refLevel+1,i+1,j  ,k+1));
      children.push_back(getGlobalID(refLevel+1,i  ,j+1,k+1));
      children.push_back(getGlobalID(refLevel+1,i+1,j+1,k+1));      
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const LID& localID) const {
      if (localID >= localToGlobalMap.size()) return invalidGlobalID();
      
      #ifdef DEBUG_AMR_MESH
      bool ok=true;
      const GID globalID = localToGlobalMap[localID];
      if (globalToLocalMap.find(globalID) == globalToLocalMap.end()) ok=false;
      if (ok == false) {
         std::cerr << "ERROR in getGlobalID, localID=" << localID << " and globalID=" << globalID;
         std::cerr << " but globalToLocalMap does not contain the block" << std::endl;
         std::cerr << "\t exiting from " << __FILE__ << ':' << __LINE__ << std::endl;
         exit(1);
      }
      #endif
      
      return localToGlobalMap[localID];
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const uint8_t& refLevel,LID indices[3]) {
      const uint8_t multiplier = std::pow(2,refLevel);
      if (indices[0] >= gridLength[0]*multiplier) return invalidGlobalID();
      if (indices[1] >= gridLength[1]*multiplier) return invalidGlobalID();
      if (indices[2] >= gridLength[2]*multiplier) return invalidGlobalID();
      return offsets[refLevel] + indices[2]*gridLength[1]*gridLength[0]*multiplier*multiplier 
	   + indices[1]*gridLength[0]*multiplier + indices[0];
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const uint32_t& refLevel,const LID& i,const LID& j,const LID& k) {
      const uint32_t multiplier = std::pow(2,refLevel);
      if (i >= gridLength[0]*multiplier) return invalidGlobalID();
      if (j >= gridLength[1]*multiplier) return invalidGlobalID();
      if (k >= gridLength[2]*multiplier) return invalidGlobalID();
      return offsets[refLevel] + k*gridLength[1]*gridLength[0]*multiplier*multiplier + j*gridLength[0]*multiplier + i;
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const uint8_t& refLevel,const Real* coords) {
      if (coords[0] < meshMinLimits[0] || coords[0] >= meshMaxLimits[0] || 
         (coords[1] < meshMinLimits[1] || coords[1] >= meshMaxLimits[1] || 
          coords[2] < meshMinLimits[2] || coords[2] >= meshMaxLimits[2])) {
         return invalidGlobalID();
      }
      
      const LID multiplier = std::pow(2,refLevel);
      const LID indices[3] = {
         static_cast<LID>(floor((coords[0] - meshMinLimits[0])*multiplier / blockSize[0])),
         static_cast<LID>(floor((coords[1] - meshMinLimits[1])*multiplier / blockSize[1])),
         static_cast<LID>(floor((coords[2] - meshMinLimits[2])*multiplier / blockSize[2]))
      };

      return offsets[refLevel] + indices[2]*gridLength[1]*gridLength[0]*multiplier*multiplier 
	     + indices[1]*gridLength[0]*multiplier + indices[0];
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalIndexOffset(const uint8_t& refLevel) {
      return offsets[refLevel];
   }
   
   template<typename GID,typename LID> inline
   const LID* VelocityMesh<GID,LID>::getGridLength(const uint8_t& refLevel) {
      return gridLengths + refLevel*3;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getIndices(const GID& globalID,uint8_t& refLevel,LID& i,LID& j,LID& k) {
      refLevel = std::upper_bound(offsets.begin(),offsets.end(),globalID)-offsets.begin()-1;
      const GID cellOffset = offsets[refLevel];

      const GID multiplier = std::pow(2,refLevel);
      const GID Nx = gridLength[0] * multiplier;
      const GID Ny = gridLength[1] * multiplier;

      GID index = globalID - cellOffset;
      k = index / (Ny*Nx);
      index -= k*Ny*Nx;
      j = index / Nx;
      i = index - j*Nx;
   }
      
   template<typename GID,typename LID> inline
   LID VelocityMesh<GID,LID>::getLocalID(const GID& globalID) const {
      typename std::unordered_map<GID,LID>::const_iterator it = globalToLocalMap.find(globalID);
      if (it != globalToLocalMap.end()) return it->second;
      return invalidLocalID();
   }

   template<typename GID,typename LID> inline
   uint8_t VelocityMesh<GID,LID>::getMaxAllowedRefinementLevel() {
      return refLevelMaxAllowed;
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getMaxVelocityBlocks() {
      return max_velocity_blocks;
   }
   
   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getMeshMaxLimits() {
      return meshMaxLimits;
   }
   
   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getMeshMinLimits() {
      return meshMinLimits;
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getNeighborsAtSameLevel(const GID& globalID,std::vector<GID>& neighborIDs) {
      neighborIDs.resize(27);
      
      // Calculate block refinement level and indices
      uint8_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      
      // Calculate global IDs of all 27 blocks:
      const LID Nx_max = gridLength[0] * std::pow(2,refLevel);
      const LID Ny_max = gridLength[1] * std::pow(2,refLevel);
      const LID Nz_max = gridLength[2] * std::pow(2,refLevel);
      
      int nbr = 0;
      for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) {
         if (i+i_off < Nx_max && (j+j_off < Ny_max && k+k_off < Nz_max)) neighborIDs[nbr] = getGlobalID(refLevel,i+i_off,j+j_off,k+k_off);
         else neighborIDs[nbr] = invalidGlobalID();
         ++nbr;
      }
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getNeighborsExistingAtSameLevel(const GID& globalID,std::vector<GID>& neighborIDs) {
      neighborIDs.clear();

      // Calculate block refinement level and indices
      uint8_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);

      // Calculate global IDs of all 26 neighbors. Note that 
      // some of the neighbors may not actually exist.
      const LID Nx_max = gridLength[0] * std::pow(2,refLevel);
      const LID Ny_max = gridLength[1] * std::pow(2,refLevel);
      const LID Nz_max = gridLength[2] * std::pow(2,refLevel);
      for (int k_off=-1; k_off<2; ++k_off) {
         if (k+k_off >= Nz_max) continue;
         for (int j_off=-1; j_off<2; ++j_off) {
            if (j+j_off >= Ny_max) continue;
            for (int i_off=-1; i_off<2; ++i_off) {
               if (i+i_off >= Nx_max) continue;
               if (i_off == 0 && (j_off == 0 && k_off == 0)) continue;
               neighborIDs.push_back(getGlobalID(refLevel,i+i_off,j+j_off,k+k_off));
            }
         }
      }
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getNeighborsExistingAtOffset(const GID& globalID,const int& i_off,const int& j_off,const int& k_off,std::vector<LID>& neighborLocalIDs,int32_t& refLevelDifference) const {
      #ifdef DEBUG_AMR_MESH
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

      // First check if the neighbor at same refinement level exists
      typename std::unordered_map<GID,LID>::const_iterator nbr;
      GID nbrGlobalID = getGlobalID(refLevel,i+i_off,j+j_off,k+k_off);
      if (nbrGlobalID == invalidGlobalID()) return;

      nbr = globalToLocalMap.find(nbrGlobalID);
      if (nbr != globalToLocalMap.end()) {
         neighborLocalIDs.push_back(nbr->second);
         refLevelDifference = 0;
         return;
      }

      // Check if parent of the neighbor exists
      nbr = globalToLocalMap.find(getParent(nbrGlobalID));
      if (nbr != globalToLocalMap.end()) {
         neighborLocalIDs.push_back(nbr->second);
         refLevelDifference = -1;
         return;
      }

      // Exit if neighbor cannot have children:
      if (refLevel == refLevelMaxAllowed) return;
      refLevelDifference = +1;

      // Check if neighbor's children exist:
      std::vector<GID> nbrChildren;
      getChildren(nbrGlobalID,nbrChildren);

      // Block can have zero to four face neighbors because of the sparse grid.
      // We need to return four neighbors so that the code calling this function 
      // indexes neighborLocalIDs correctly, non-existing neighbors have an invalid local ID.
      neighborLocalIDs.resize(4);
      for (size_t n=0; n<4; ++n) neighborLocalIDs[n] = invalidLocalID();

      const int nbrTypeID = (k_off+1)*9+(j_off+1)*3+(i_off+1);
      switch (nbrTypeID) {
       case 4:  // -z face neighbor
         nbr = globalToLocalMap.find(nbrChildren[4]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[0] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[5]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[1] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[6]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[7]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
         break;
       case 10: // -y face neighbor
         nbr = globalToLocalMap.find(nbrChildren[2]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[0] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[3]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[1] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[6]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[7]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
         break;
       case 12: // -x face neighbor
         nbr = globalToLocalMap.find(nbrChildren[1]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[0] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[3]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[1] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[5]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[7]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
         break;
       case 14: // +x face neighbor
         nbr = globalToLocalMap.find(nbrChildren[0]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[0] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[2]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[1] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[4]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[6]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
         break;
       case 16: // +y face neighbor
         nbr = globalToLocalMap.find(nbrChildren[0]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[0] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[1]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[1] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[4]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[5]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
         break;
       case 22: // +z face neighbor
         nbr = globalToLocalMap.find(nbrChildren[0]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[0] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[1]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[1] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[2]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
         nbr = globalToLocalMap.find(nbrChildren[3]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
         break;
       default:
         break;
      }
   }

   template<typename GID,typename LID> inline
   int VelocityMesh<GID,LID>::getOctant(const GID& globalID) {
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
   GID VelocityMesh<GID,LID>::getParent(const GID& globalID) {
      // Calculate block indices and refinement level
      uint8_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);

      // Blocks at base grid level don't have a parent
      if (refLevel == 0) return globalID;

      // Calculate parent's global ID
      i /= 2;
      j /= 2;
      k /= 2;
      return getGlobalID(refLevel-1,i,j,k);
   }
   
   template<typename GID,typename LID> inline
   uint8_t VelocityMesh<GID,LID>::getRefinementLevel(const GID& globalID) {
      for (size_t i=1; i<refLevelMaxAllowed+1; ++i) {
	 if (globalID < offsets[i]) return i-1;
      }
      return refLevelMaxAllowed;
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getSiblingNeighbors(const GID& globalID,std::vector<GID>& nbrs) {
      nbrs.clear();
      
      // Get block indices
      uint8_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      
      // Shift indices to bottom lower left corner in the octant
      i -= (i % 2);
      j -= (j % 2);
      k -= (k % 2);
      
      // Calculate the global IDs of all common neighbors of the blocks in the octant
      for (LID k_off=-1; k_off<3; ++k_off) {
         for (LID j_off=-1; j_off<3; ++j_off) {
            for (LID i_off=-1; i_off<3; ++i_off) {
               int cntr=0;
               if (i_off == 0 || i_off == 1) ++cntr;
               if (j_off == 0 || j_off == 1) ++cntr;
               if (k_off == 0 || k_off == 1) ++cntr;
               if (cntr == 3) continue;
               nbrs.push_back(getGlobalID(refLevel,i+i_off,j+j_off,k+k_off));
            }
         }
      }
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getSiblings(const GID& globalID,std::vector<GID>& siblings) {
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
      std::vector<GID> children;
      getChildren(globalID,children);
      
      // If children.size() is zero, block is at max refinement level:
      if (children.size() == 0) return false;

      // Check all children, if even one exists return 'true':
      for (size_t c=0; c<children.size(); ++c) {
         if (getLocalID(children[c]) != invalidLocalID()) return true;
      }
      return false;
   }
   
   /** Check if the block has any existing grandparents.
    * @param globalID Global ID of the block.
    * @return INVALID_GLOBALID if block has no existing grandparents, otherwise 
    * the global ID of the grandparent is returned.*/
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::hasGrandParent(const GID& globalID) const {
       // Exit if the block cannot have grandparents
       uint8_t refLevel = getRefinementLevel(globalID);
       if (refLevel <= 1) return invalidGlobalID();
       
       // Calculate the global ID of the first possible grandparent
       GID parentGID = getParent(globalID);
       parentGID = getParent(parentGID);
       
       do {
           // Grandparent exists, return its global ID
           if (getLocalID(parentGID) != invalidLocalID()) return parentGID;
           
           // Move down one refinement level
           GID tmp = parentGID;
           parentGID = getParent(parentGID);
           if (tmp == parentGID) break;
       } while (true);

       return invalidGlobalID();
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::initialize(Real meshLimits[6],LID gridLength[3],LID blockLength[3],uint8_t refLevelMaxAllowed,
                                          LID maxVelocityBlocks) {
      if (initialized == true) return initialized;
      VelocityMesh<GID,LID>::refLevelMaxAllowed = refLevelMaxAllowed;
      
      meshMinLimits[0] = meshLimits[0];
      meshMinLimits[1] = meshLimits[2];
      meshMinLimits[2] = meshLimits[4];
      meshMaxLimits[0] = meshLimits[1];
      meshMaxLimits[1] = meshLimits[3];
      meshMaxLimits[2] = meshLimits[5];
	
      VelocityMesh<GID,LID>::gridLength[0] = gridLength[0];
      VelocityMesh<GID,LID>::gridLength[1] = gridLength[1];
      VelocityMesh<GID,LID>::gridLength[2] = gridLength[2];
      
      VelocityMesh<GID,LID>::blockLength[0] = blockLength[0];
      VelocityMesh<GID,LID>::blockLength[1] = blockLength[1];
      VelocityMesh<GID,LID>::blockLength[2] = blockLength[2];
      
      // Calculate derived mesh parameters:
      gridSize[0] = meshMaxLimits[0] - meshMinLimits[0];
      gridSize[1] = meshMaxLimits[1] - meshMinLimits[1];
      gridSize[2] = meshMaxLimits[2] - meshMinLimits[2];

      VelocityMesh<GID,LID>::blockSize[0] = gridSize[0] / gridLength[0];
      VelocityMesh<GID,LID>::blockSize[1] = gridSize[1] / gridLength[1];
      VelocityMesh<GID,LID>::blockSize[2] = gridSize[2] / gridLength[2];

      cellSize[0] = blockSize[0] / blockLength[0];
      cellSize[1] = blockSize[1] / blockLength[1];
      cellSize[2] = blockSize[2] / blockLength[2];

      max_velocity_blocks = maxVelocityBlocks;

      // Calculate global ID offsets:
      const GID N_blocks0 = gridLength[0]*gridLength[1]*gridLength[2];
      offsets.resize(refLevelMaxAllowed+1);
      offsets[0] = 0;
      for (size_t i=1; i<refLevelMaxAllowed+1; ++i) {
         offsets[i] = offsets[i-1] + N_blocks0 * std::pow(8,i-1);
      }
      
      gridLengths = new LID[3*(refLevelMaxAllowed+1)];
      blockSizes  = new Real[3*(refLevelMaxAllowed+1)];
      cellSizes   = new Real[3*(refLevelMaxAllowed+1)];
      uint32_t mul = 1;
      for (uint8_t r=0; r<refLevelMaxAllowed+1; ++r) {
         for (int i=0; i<3; ++i) gridLengths[3*r+i] = VelocityMesh<GID,LID>::gridLength[i] * mul;
         for (int i=0; i<3; ++i) blockSizes[3*r+i]  = VelocityMesh<GID,LID>::blockSize[i] / mul;
         for (int i=0; i<3; ++i) cellSizes[3*r+i]   = VelocityMesh<GID,LID>::blockSize[i] / (mul * WID);
         mul *= 2;
      }
      initialized = true;
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
   void VelocityMesh<GID,LID>::pop() {
      if (size() == 0) return;

      const vmesh::LocalID lastLID = size()-1;
      const vmesh::GlobalID lastGID = localToGlobalMap[lastLID];
	  
      #ifdef DEBUG_AMR_MESH
         bool ok = true;
         if (globalToLocalMap.find(lastGID) == globalToLocalMap.end()) ok = false;
         if (ok == false) {
            std::cerr << "ERROR in pop(): last localID=" << lastLID << " globalID=" << lastGID;
            std::cerr << " but globalToLocalMap does not contain that block" << std::endl;
            std::cerr << "\t exiting from " << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
         }
      #endif
	  
      typename std::unordered_map<GID,LID>::iterator last = globalToLocalMap.find(lastGID);
      globalToLocalMap.erase(last);
      localToGlobalMap.pop_back();
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::push_back(const GID& globalID) {
      if (globalID == invalidGlobalID()) return false;
      if (size() >= max_velocity_blocks) {
         std::cerr << "vmesh: too many blocks, current size is " << size() << " max " << max_velocity_blocks << std::endl;
         return false;
      }

      std::pair<typename std::unordered_map<GID,LID>::iterator,bool> position
        = globalToLocalMap.insert(std::make_pair(globalID,localToGlobalMap.size()));

      if (position.second == true) {	 
         localToGlobalMap.push_back(globalID);
      }

      return position.second;
   }
   
   template<typename GID,typename LID> inline
   uint8_t VelocityMesh<GID,LID>::push_back(const std::vector<vmesh::GlobalID>& blocks) {
      if (size()+blocks.size() >= max_velocity_blocks) {
         std::cerr << "vmesh: too many blocks, current size is " << size() << " max " << max_velocity_blocks << std::endl;
         return 0;
      }

      // Attempt to add the given blocks
      uint8_t adds=0;
      for (size_t b=0; b<blocks.size(); ++b) {
         const GID globalID = blocks[b];
         std::pair<typename std::unordered_map<GID,LID>::iterator,bool> position
           = globalToLocalMap.insert(std::make_pair(globalID,localToGlobalMap.size()+b));
         if (position.second == true) {
            localToGlobalMap.push_back(globalID);
            ++adds;
         }
      }

      return adds;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::refine(const GID& globalID,std::set<GID>& erasedBlocks,std::map<GID,LID>& insertedBlocks) {
      // Check that the block exists
      typename std::unordered_map<GID,LID>::iterator it = globalToLocalMap.find(globalID);
      if (it == globalToLocalMap.end()) {
         return false;
      }
      const LID localID = it->second;

      // Calculate block's current refinement level and indices
      uint8_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);

      // Check that the block can still be refined
      if (refLevel >= getMaxAllowedRefinementLevel()) {
         return false;
      }

      // Calculate childrens' global IDs
      i *= 2;
      j *= 2;
      k *= 2;
      GID childrenGlobalIDs[8];
      childrenGlobalIDs[0] = getGlobalID(refLevel+1,i  ,j  ,k  );
      childrenGlobalIDs[1] = getGlobalID(refLevel+1,i+1,j  ,k  );
      childrenGlobalIDs[2] = getGlobalID(refLevel+1,i  ,j+1,k  );
      childrenGlobalIDs[3] = getGlobalID(refLevel+1,i+1,j+1,k  );
      childrenGlobalIDs[4] = getGlobalID(refLevel+1,i  ,j  ,k+1);
      childrenGlobalIDs[5] = getGlobalID(refLevel+1,i+1,j  ,k+1);
      childrenGlobalIDs[6] = getGlobalID(refLevel+1,i  ,j+1,k+1);
      childrenGlobalIDs[7] = getGlobalID(refLevel+1,i+1,j+1,k+1);

      // Erase the refined block
      globalToLocalMap.erase(globalID);
      std::pair<typename std::set<GID>::iterator,bool> success = erasedBlocks.insert(globalID);

      // Insert new blocks:
      const size_t currentSize = localToGlobalMap.size();

      // Old block is replaced with the first child
      globalToLocalMap.insert(std::make_pair(childrenGlobalIDs[0],localID));
      localToGlobalMap[localID] = childrenGlobalIDs[0];
      insertedBlocks.insert(std::make_pair(childrenGlobalIDs[0],localID));

      for (int i=1; i<8; ++i) {
         push_back(childrenGlobalIDs[i]);
         insertedBlocks.insert(std::make_pair(childrenGlobalIDs[i],localToGlobalMap.size()-1));
      }

      // Enforce that neighbors have at maximum one refinement level difference
      std::vector<GID> nbrs;
      getNeighborsExistingAtSameLevel(globalID,nbrs);
      for (size_t n=0; n<nbrs.size(); ++n) {
         // Get parent of the neighbor:
         const GID parentID = getParent(nbrs[n]);

         // If the parent exists, it is at two refinement levels higher than 
         // the block that was refined above. Refine the parent of the neighbor
         if (parentID == nbrs[n]) continue;
         it = globalToLocalMap.find(parentID);
         if (it != globalToLocalMap.end()) {
            refine(parentID,erasedBlocks,insertedBlocks);
         }
      }

      return true;
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
   size_t VelocityMesh<GID,LID>::size() const {
      return globalToLocalMap.size();
   }
   
   template<typename GID,typename LID> inline
   size_t VelocityMesh<GID,LID>::sizeInBytes() const {
      return globalToLocalMap.size()*sizeof(GID)
	+ localToGlobalMap.size()*(sizeof(GID)+sizeof(LID))
	+ offsets.size()*sizeof(GID);      
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::swap(VelocityMesh& vm) {
      globalToLocalMap.swap(vm.globalToLocalMap);
      localToGlobalMap.swap(vm.localToGlobalMap);      
   }

} // namespace vmesh
   
#endif
