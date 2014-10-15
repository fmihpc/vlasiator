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

namespace vmesh {

   template<typename GID,typename LID>
   class VelocityMesh {
    public:
      VelocityMesh();
      VelocityMesh(const VelocityMesh& vm);
      ~VelocityMesh();
      
      size_t capacityInBytes() const;
      bool check() const;
      void clear();
      bool copy(const LID& sourceLocalID,const LID& targetLocalID);
      size_t count(const GID& globalID) const;
      static const GID* getBaseGridLength();
      static const Real* getBaseGridBlockSize();
      static const Real* getBaseGridCellSize();
      static bool getBlockCoordinates(const GID& globalID,Real coords[3]);
      static void getBlockInfo(const GID& globalID,Real* array);
      static bool getBlockSize(const GID& globalID,Real size[3]);
      static bool getCellSize(const GID& globalID,Real size[3]);
      static void getChildren(const GID& globalID,std::vector<GID>& children);
      GID getGlobalID(const LID& localID) const;
      static GID getGlobalID(const Real& x,const Real& y,const Real& z);
      static GID getGlobalID(const uint32_t& refLevel,const LID& i,const LID& j,const LID& k);
      static void getIndices(const GID& globalID,uint32_t& refLevel,LID& i,LID& j,LID& k);
      LID getLocalID(const GID& globalID) const;
      static uint8_t getMaxAllowedRefinementLevel();
      static GID getMaxVelocityBlocks();
      static const Real* getMeshMaxLimits();
      static const Real* getMeshMinLimits();
      static void getNeighborsAtSameLevel(const GID& globalID,std::vector<GID>& neighborIDs);
      void getNeighborsExistingAtOffset(const GID& globalID,const int& i,const int& j,const int& k,std::vector<LID>& neighborIDs,int32_t& refLevelDifference) const;
      static int getOctant(const GID& globalID);
      static GID getParent(const GID& globalID);
      void getSiblingNeighbors(const GID& globalID,std::vector<GID>& nbrs);
      void getSiblings(const GID& globalID,GID siblings[8]);
      void getSiblings(const GID& globalID,std::vector<GID>& siblings);
      static bool initialize(Real meshLimits[6],LID gridLength[3],LID blockLength[3],uint8_t refLevelMaxAllowed);
      static LID invalidBlockIndex();
      static GID invalidGlobalID();
      static LID invalidLocalID();
      void pop();
      bool push_back(const GID& globalID);
      bool refine(const GID& globalID,std::set<GID>& erasedBlocks,std::map<GID,LID>& insertedBlocks);
      bool setGrid(const std::vector<GID>& globalIDs);
      size_t size() const;
      size_t sizeInBytes() const;

    private:
      VelocityMesh& operator=(const VelocityMesh& vm);

      static LID max_velocity_blocks;                                           /**< Maximum valid block local ID.*/
      static LID blockLength[3];                                                /**< Number of cells in a block per coordinate.*/
      static Real blockSize[3];                                                 /**< Size of a block at base grid level.*/
      static Real cellSize[3];                                                  /**< Size of a cell in a block at base grid level.*/
      static Real gridSize[3];                                                  /**< Size of the grid.*/
      static GID gridLength[3];                                                 /**< Number of blocks in the grid.*/
      static Real meshMinLimits[3];                                             /**< Minimum coordinate values of the grid bounding box.*/
      static Real meshMaxLimits[3];                                             /**< Maximum coordinate values of the grid bounding box.*/
      static uint8_t refLevelMaxAllowed;      
      static std::vector<GID> offsets;                                          /**< Block global ID offsets for each refinement level.*/

      std::vector<GID> localToGlobalMap;
      std::unordered_map<GID,LID> globalToLocalMap;
   };

   // ***** INITIALIZERS FOR STATIC MEMBER VARIABLES ***** //
   
   template<typename GID,typename LID> LID VelocityMesh<GID,LID>::max_velocity_blocks = 0;
   template<typename GID,typename LID> LID VelocityMesh<GID,LID>::blockLength[3] = {0,0,0};
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::blockSize[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::cellSize[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::gridSize[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> GID VelocityMesh<GID,LID>::gridLength[3] = {0,0,0};
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::meshMinLimits[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> Real VelocityMesh<GID,LID>::meshMaxLimits[3] = {NAN,NAN,NAN};
   template<typename GID,typename LID> std::vector<GID> VelocityMesh<GID,LID>::offsets;
   template<typename GID,typename LID> uint8_t VelocityMesh<GID,LID>::refLevelMaxAllowed = 0;
   
   // ***** DEFINITIONS OF TEMPLATE MEMBER FUNCTIONS ***** //
   
   template<typename GID,typename LID> inline
   VelocityMesh<GID,LID>::VelocityMesh() { }

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
	 typename std::unordered_map<GID,LID>::const_iterator it = globalToLocalMap.find(globalID);
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
      std::unordered_map<GID,LID>().swap(globalToLocalMap);
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
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getBlockCoordinates(const GID& globalID,Real coords[3]) {
      if (globalID == invalidGlobalID()) {
	 for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
	 return false;
      }
      
      // Calculate block's refinement level and (i,j,k) indices:
      uint32_t refLevel;
      LID indices[3];
      getIndices(globalID,refLevel,indices[0],indices[1],indices[2]);
      if (indices[0] == invalidBlockIndex()) {
	 for (int i=0; i<3; ++i) coords[i] = std::numeric_limits<Real>::quiet_NaN();
	 return false;
      }
      
      // NUmber of blocks per coordinate at this refinement level:
      const GlobalID multiplier = pow(2,refLevel);
      
      coords[0] = meshMinLimits[0] + indices[0]*blockSize[0]/multiplier;
      coords[1] = meshMinLimits[1] + indices[1]*blockSize[1]/multiplier;
      coords[2] = meshMinLimits[2] + indices[2]*blockSize[2]/multiplier;
      
      //std::cerr << "\t block " << globalID << " " << indices[0] << " " << indices[1] << " " << indices[2] << " lvl " << refLevel << " crds: ";
      //std::cerr << coords[0] << '\t' << coords[1] << '\t' << coords[2] << std::endl;
      //std::cerr << "\t offsets " << offsets[0] << '\t' << offsets[1] << std::endl;
      return true;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getBlockInfo(const GID& globalID,Real* array) {
      #ifndef NDEBUG
      if (globalID == invalidGlobalID()) {
	 for (int i=0; i<6; ++i) array[i] = std::numeric_limits<Real>::infinity();
      }
      #endif

      getBlockCoordinates(globalID,array+0);
      getCellSize(globalID,array+3);
   }
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getBlockSize(const GID& globalID,Real size[3]) {
      // Calculate block's refinement level and (i,j,k) indices:
      uint32_t refLevel;
      LID indices[3];
      getIndices(globalID,refLevel,indices[0],indices[1],indices[2]);

      // NUmber of blocks per coordinate at this refinement level:
      const GlobalID multiplier = pow(2,refLevel);
      
      // Calculate the number of blocks in each coordinate direction at this refinement level:
      size[0] = blockSize[0] / multiplier;
      size[1] = blockSize[1] / multiplier;
      size[2] = blockSize[2] / multiplier;
      return true;
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

      uint32_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      if (refLevel == refLevelMaxAllowed) return;

      i *= 2;
      j *= 2;
      k *= 2;

      children.push_back(getGlobalID(refLevel+1,i  ,j  ,k  ));
      children.push_back(getGlobalID(refLevel+1,i+1,j  ,k  ));
      children.push_back(getGlobalID(refLevel+1,i+1,j+1,k  ));
      children.push_back(getGlobalID(refLevel+1,i  ,j+1,k  ));
      children.push_back(getGlobalID(refLevel+1,i  ,j  ,k+1));
      children.push_back(getGlobalID(refLevel+1,i+1,j  ,k+1));
      children.push_back(getGlobalID(refLevel+1,i+1,j+1,k+1));
      children.push_back(getGlobalID(refLevel+1,i  ,j+1,k+1));
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
   GID VelocityMesh<GID,LID>::getGlobalID(const uint32_t& refLevel,const LID& i,const LID& j,const LID& k) {
      const uint32_t multiplier = pow(2,refLevel);
      if (i >= gridLength[0]*multiplier) return invalidGlobalID();
      if (j >= gridLength[1]*multiplier) return invalidGlobalID();
      if (k >= gridLength[2]*multiplier) return invalidGlobalID();
      return offsets[refLevel] + k*gridLength[1]*gridLength[0]*multiplier*multiplier + j*gridLength[0]*multiplier + i;
   }

   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const Real& x,const Real& y,const Real& z) {
      #warning This version of getGlobalID returns global IDs on base grid level
      if (x < meshMinLimits[0] || x >= meshMaxLimits[0] || y < meshMinLimits[1] || y >= meshMaxLimits[1] || z < meshMinLimits[2] || z >= meshMaxLimits[2]) {
	 return invalidGlobalID();
      }
      
      const LID indices[3] = {
	 static_cast<LID>(floor((x - meshMinLimits[0]) / blockSize[0])),
	 static_cast<LID>(floor((y - meshMinLimits[1]) / blockSize[1])),
	 static_cast<LID>(floor((z - meshMinLimits[2]) / blockSize[2]))
      };

      return getGlobalID(0,indices[0],indices[1],indices[2]);
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getIndices(const GID& globalID,uint32_t& refLevel,LID& i,LID& j,LID& k) {
      refLevel = std::upper_bound(offsets.begin(),offsets.end(),globalID)-offsets.begin()-1;
      const GID cellOffset = offsets[refLevel];

      const GlobalID multiplier = pow(2,refLevel);
      const GlobalID Nx = gridLength[0] * multiplier;
      const GlobalID Ny = gridLength[1] * multiplier;
      
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
      neighborIDs.clear();

      // Calculate block refinement level and indices
      uint32_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);

      // Calculate global IDs of all 26 neighbors. Note that 
      // some of the neighbors may not actually exist.
      const LID Nx_max = gridLength[0] * pow(2,refLevel);
      const LID Ny_max = gridLength[1] * pow(2,refLevel);
      const LID Nz_max = gridLength[2] * pow(2,refLevel);
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
      uint32_t refLevel;
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
	 nbr = globalToLocalMap.find(nbrChildren[7]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[6]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
	 break;
       case 10: // -y face neighbor
	 nbr = globalToLocalMap.find(nbrChildren[3]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[0] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[2]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[1] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[7]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[6]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
	 break;
       case 12: // -x face neighbor
	 nbr = globalToLocalMap.find(nbrChildren[1]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[0] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[2]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[1] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[5]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[6]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
	 break;
       case 14: // +x face neighbor
	 nbr = globalToLocalMap.find(nbrChildren[0]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[0] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[3]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[1] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[4]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[7]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
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
	 nbr = globalToLocalMap.find(nbrChildren[3]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[2] = nbr->second;}
	 nbr = globalToLocalMap.find(nbrChildren[2]); if (nbr != globalToLocalMap.end()) {neighborLocalIDs[3] = nbr->second;}
	 break;
       default:
	 break;
      }
   }
   
   template<typename GID,typename LID> inline
   int VelocityMesh<GID,LID>::getOctant(const GID& globalID) {
      // Calculate block indices and refinement level
      uint32_t refLevel;
      LID i,j,k;
      getIndices(globalID,refLevel,i,j,k);
      if (refLevel == 0) return 0;

      const int i_oct = i % 2;
      const int j_oct = j % 2;
      const int k_oct = k % 2;
      return k_oct*2*2 + j_oct*2 + i_oct;
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getParent(const GID& globalID) {
      // Calculate block indices and refinement level
      uint32_t refLevel;
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
   bool VelocityMesh<GID,LID>::initialize(Real meshLimits[6],LID gridLength[3],LID blockLength[3],uint8_t refLevelMaxAllowed) {
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

      #warning FIXME maxblocks
      //max_velocity_blocks = gridLength[0]*gridLength[1]*gridLength[2];
      max_velocity_blocks = 50000;
      
      // Calculate global ID offsets:
      const GID N_blocks0 = gridLength[0]*gridLength[1]*gridLength[2];
      offsets.resize(refLevelMaxAllowed+1);
      offsets[0] = 0;
      for (size_t i=1; i<refLevelMaxAllowed+1; ++i) {
	 offsets[i] = offsets[i-1] + N_blocks0 * pow(8,i-1);
      }
      
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
      typename std::unordered_map<GID,LID>::iterator last = globalToLocalMap.find(lastGID);
	
      globalToLocalMap.erase(last);
      localToGlobalMap.pop_back();
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::push_back(const GID& globalID) {
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
   bool VelocityMesh<GID,LID>::refine(const GID& globalID,std::set<GID>& erasedBlocks,std::map<GID,LID>& insertedBlocks) {
      // Check that the block exists
      typename std::unordered_map<GID,LID>::iterator it = globalToLocalMap.find(globalID);
      if (it == globalToLocalMap.end()) {
	 return false;
      }
      const LID localID = it->second;

      // Calculate block's current refinement level and indices
      uint32_t refLevel;
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
      getNeighborsAtSameLevel(globalID,nbrs);
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

} // namespace vmesh
   
#endif
