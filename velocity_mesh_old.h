/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2014 Finnish Meteorological Institute
 */

#ifndef VELOCITY_MESH_OLD_H
#define VELOCITY_MESH_OLD_H

#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <cmath>

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
      //static const GID* getBaseGridLength();
      //static const Real* getBaseGridBlockSize();
      //static const Real* getBaseGridCellSize();
      static bool getBlockCoordinates(const GID& globalID,Real coords[3]);
      static void getBlockInfo(const GID& globalID,Real* array);
      static const Real* getBlockSize(const uint8_t& refLevel);
      static bool getBlockSize(const GID& globalID,Real size[3]);
      static const Real* getCellSize(const uint8_t& refLevel);
      static bool getCellSize(const GID& globalID,Real size[3]);
      static void getChildren(const GID& globalID,std::vector<GID>& children);
//      void     getChildren(const GlobalID& globalID,std::vector<GlobalID>& children);
      GID getGlobalID(const LID& localID) const;
      static GID getGlobalID(const uint8_t& refLevel,const Real* coords);
      static GID getGlobalID(const uint8_t& refLevel,GID indices[3]);
      static GID getGlobalID(const uint32_t& refLevel,const LID& i,const LID& j,const LID& k);
      static GID getGlobalIndexOffset(const uint8_t& refLevel=0);
      static const GID* getGridLength(const uint8_t& refLevel);
//      void     getNeighbors(const GlobalID& globalID,std::vector<GlobalID>& neighborIDs);
      static void getIndices(const GID& globalID,uint8_t& refLevel,LID& i,LID& j,LID& k);
      LID getLocalID(const GID& globalID) const;
      static uint8_t getMaxAllowedRefinementLevel();
      static GID getMaxVelocityBlocks();
      static const Real* getMeshMaxLimits();
      static const Real* getMeshMinLimits();
      static void getNeighborsAtSameLevel(const GID& globalID,std::vector<GID>& neighborIDs);
      void getNeighborsExistingAtOffset(const GID& globalID,const int& i,const int& j,const int& k,std::vector<LID>& neighborIDs,int32_t& refLevelDifference) const;
      static int getOctant(const GID& globalID);
      static GID getParent(const GID& globalID);
      static uint8_t getRefinementLevel(const GID& globalID);
//      void     getSiblingNeighbors(const GlobalID& globalID,std::vector<GlobalID>& nbrs);
//      void     getSiblings(const GlobalID& globalID,GlobalID siblings[8]);
      static void getSiblings(const GID& globalID,std::vector<GID>& siblings);
      bool hasChildren(const GID& globalID) const;
      static bool initialize(Real meshLimits[6],LID gridLength[3],LID blockLength[3],uint8_t refLevelMaxAllowed=0);
      static LID invalidBlockIndex();
      static GID invalidGlobalID();
      static LID invalidLocalID();
      void pop();
      bool push_back(const GID& globalID);
      bool push_back(const std::vector<vmesh::GlobalID>& blocks);
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
      static GID gridLength[3];                                                 /**< Number of blocks in the grid.*/
      static Real meshMinLimits[3];                                             /**< Minimum coordinate values of the grid bounding box.*/
      static Real meshMaxLimits[3];                                             /**< Maximum coordinate values of the grid bounding box.*/
      
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
   GID VelocityMesh<GID,LID>::findBlockDown(uint8_t& refLevel,GID cellIndices[3]) {
      // Calculate i/j/k indices of the block that would own the cell:
      GID i_block = cellIndices[0] / blockLength[0];
      GID j_block = cellIndices[1] / blockLength[1];
      GID k_block = cellIndices[2] / blockLength[2];
      
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
   GID VelocityMesh<GID,LID>::findBlock(uint8_t& refLevel,GID cellIndices[3]) {
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
   bool VelocityMesh<GID,LID>::getBlockCoordinates(const GID& globalID,Real coords[3]) {
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

      coords[0] = meshMinLimits[0] + indices[0]*blockSize[0];
      coords[1] = meshMinLimits[1] + indices[1]*blockSize[1];
      coords[2] = meshMinLimits[2] + indices[2]*blockSize[2];
      return true;
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getBlockInfo(const GID& globalID,Real* array) {
      #ifndef NDEBUG
      if (globalID == invalidGlobalID()) {
         for (int i=0; i<6; ++i) array[i] = std::numeric_limits<Real>::infinity();
      }
      #endif

      LID indices[3];
      indices[0] = globalID % gridLength[0];
      indices[1] = (globalID / gridLength[0]) % gridLength[1];
      indices[2] = globalID / (gridLength[0] * gridLength[1]);

      // Indices 0-2 contain coordinates of the lower left corner.
      // The values are the same as if getBlockCoordinates(globalID,&(array[0])) was called
      array[0] = meshMinLimits[0] + indices[0]*blockSize[0];
      array[1] = meshMinLimits[1] + indices[1]*blockSize[1];
      array[2] = meshMinLimits[2] + indices[2]*blockSize[2];

      // Indices 3-5 contain the cell size.
      // The values are the same as if getCellSize(globalID,&(array[3])) was called
      array[3] = cellSize[0];
      array[4] = cellSize[1];
      array[5] = cellSize[2];
   }

   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getBlockSize(const uint8_t& refLevel) {
      return blockSize;
   }
   
   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getBlockSize(const GID& globalID,Real size[3]) {
      size[0] = blockSize[0];
      size[1] = blockSize[1];
      size[2] = blockSize[2];
      return true;
   }
   
   template<typename GID,typename LID> inline
   const Real* VelocityMesh<GID,LID>::getCellSize(const uint8_t& refLevel) {
      return cellSize;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::getCellSize(const GID& globalID,Real size[3]) {
      size[0] = cellSize[0];
      size[1] = cellSize[1];
      size[2] = cellSize[2];
      return true;
   }
   
   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getChildren(const GID& globalID,std::vector<GID>& children) {
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
   GID VelocityMesh<GID,LID>::getGlobalID(const uint8_t& refLevel,const Real* coords) {
      if (coords[0] < meshMinLimits[0] || coords[0] >= meshMaxLimits[0] ||
         (coords[1] < meshMinLimits[1] || coords[1] >= meshMaxLimits[1] ||
          coords[2] < meshMinLimits[2] || coords[2] >= meshMaxLimits[2])) {
         return invalidGlobalID();
      }

      const LID indices[3] = {
         static_cast<LID>(floor((coords[0] - meshMinLimits[0]) / blockSize[0])),
         static_cast<LID>(floor((coords[1] - meshMinLimits[1]) / blockSize[1])),
         static_cast<LID>(floor((coords[2] - meshMinLimits[2]) / blockSize[2]))
      };

      return indices[2]*gridLength[1]*gridLength[0] + indices[1]*gridLength[0] + indices[0];
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const uint8_t& refLevel,GID indices[3]) {
      if (indices[0] >= gridLength[0]) return invalidGlobalID();
      if (indices[1] >= gridLength[1]) return invalidGlobalID();
      if (indices[2] >= gridLength[2]) return invalidGlobalID();
      return indices[2]*gridLength[1]*gridLength[0] + indices[1]*gridLength[0] + indices[0];
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalID(const uint32_t& refLevel,const LID& i,const LID& j,const LID& k) {
      if (i >= gridLength[0] || j >= gridLength[1] || k >= gridLength[2]) {
         return invalidGlobalID();
      }
      return i + j*gridLength[0] + k*gridLength[0]*gridLength[1];
   }
   
   template<typename GID,typename LID> inline
   GID VelocityMesh<GID,LID>::getGlobalIndexOffset(const uint8_t& refLevel) {
      return 0;
   }

   template<typename GID,typename LID> inline
   const GID* VelocityMesh<GID,LID>::getGridLength(const uint8_t& refLevel) {
      return gridLength;
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::getIndices(const GID& globalID,uint8_t& refLevel,LID& i,LID& j,LID& k) {
      refLevel = 0;
      if (globalID >= invalidGlobalID()) {
         i = j = k = invalidBlockIndex();
      } else {
         i = globalID % gridLength[0];
         j = (globalID / gridLength[0]) % gridLength[1];
         k = globalID / (gridLength[0] * gridLength[1]);
      }
   }

   template<typename GID,typename LID> inline
   LID VelocityMesh<GID,LID>::getLocalID(const GID& globalID) const {
      typename std::unordered_map<GID,LID>::const_iterator it = globalToLocalMap.find(globalID);
      if (it != globalToLocalMap.end()) return it->second;
      return invalidLocalID();
   }
   
   template<typename GID,typename LID> inline
   uint8_t VelocityMesh<GID,LID>::getMaxAllowedRefinementLevel() {
      return 0;
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
      const LID Nx_max = gridLength[0];
      const LID Ny_max = gridLength[1];
      const LID Nz_max = gridLength[2];
      
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
      typename std::unordered_map<GID,LID>::const_iterator nbr;
      GID nbrGlobalID = getGlobalID(0,i+i_off,j+j_off,k+k_off);
      if (nbrGlobalID == invalidGlobalID()) return;

      nbr = globalToLocalMap.find(nbrGlobalID);
      if (nbr != globalToLocalMap.end()) {
         neighborLocalIDs.push_back(nbr->second);
         refLevelDifference = 0;
         return;
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
      return globalID;
   }
   
   template<typename GID,typename LID> inline
   uint8_t VelocityMesh<GID,LID>::getRefinementLevel(const GID& globalID) {
      return 0;
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
      return false;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::initialize(Real meshLimits[6],LID gridLength[3],LID blockLength[3],uint8_t refLevelMaxAllowed) {
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

      max_velocity_blocks = gridLength[0]*gridLength[1]*gridLength[2];
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
      if (size() >= max_velocity_blocks) return false;
      if (globalID == invalidGlobalID()) return false;

      std::pair<typename std::unordered_map<GID,LID>::iterator,bool> position
        = globalToLocalMap.insert(std::make_pair(globalID,localToGlobalMap.size()));

      if (position.second == true) {
         localToGlobalMap.push_back(globalID);
      }

      return position.second;
   }

   template<typename GID,typename LID> inline
   bool VelocityMesh<GID,LID>::push_back(const std::vector<vmesh::GlobalID>& blocks) {
      if (size()+blocks.size() >= max_velocity_blocks) {
         std::cerr << "vmesh: too many blocks, current size is " << size() << " max " << max_velocity_blocks << std::endl;
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
           + localToGlobalMap.size()*(sizeof(GID)+sizeof(LID));
   }

   template<typename GID,typename LID> inline
   void VelocityMesh<GID,LID>::swap(VelocityMesh& vm) {
      globalToLocalMap.swap(vm.globalToLocalMap);
      localToGlobalMap.swap(vm.localToGlobalMap);
   }
   
} // namespace vmesh

#endif
