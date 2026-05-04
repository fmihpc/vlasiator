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

// Non-AMR implementation of the velocity space, still in use despite the filename

#ifndef VELOCITY_MESH_GPU_H
#define VELOCITY_MESH_GPU_H

#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdint.h>
#include <vector>
#include <cmath>

#include "../velocity_mesh_parameters.h"

#include "include/hashinator/hashinator.h"
#include "include/splitvector/splitvec.h"

#include "../arch/gpu_base.hpp"
//#include "../arch/arch_device_api.h" // included in above

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

      ARCH_HOSTDEV size_t capacity() const;
      size_t capacityInBytes() const;
      ARCH_HOSTDEV bool check(); // These are no longer const as they prefetch back and forth
      ARCH_HOSTDEV void print();
      void clear(bool shrink);
      ARCH_HOSTDEV bool move(const vmesh::LocalID sourceLocalID,const vmesh::LocalID targetLocalID);
      ARCH_DEV bool warpMove(const vmesh::LocalID sourceLocalID,const vmesh::LocalID targetLocalID, const size_t b_tid);
      ARCH_HOSTDEV size_t count(const vmesh::GlobalID globalID) const;
      ARCH_DEV size_t warpCount(const vmesh::GlobalID globalID, const size_t b_tid) const;
      ARCH_HOSTDEV vmesh::GlobalID findBlock(vmesh::GlobalID cellIndices[3]) const;
      ARCH_DEV vmesh::GlobalID warpFindBlock(vmesh::GlobalID cellIndices[3], const size_t b_tid) const;
      ARCH_HOSTDEV bool getBlockCoordinates(const vmesh::GlobalID globalID,Real coords[3]) const;
      ARCH_HOSTDEV void getBlockInfo(const vmesh::GlobalID globalID,Real* array) const;
      ARCH_HOSTDEV const Real* getBlockSize() const;
      ARCH_HOSTDEV bool getBlockSize(const vmesh::GlobalID globalID,Real size[3]) const;
      ARCH_HOSTDEV const Real* getCellSize() const;
      ARCH_HOSTDEV bool getCellSize(const vmesh::GlobalID globalID,Real size[3]) const;
      ARCH_HOSTDEV vmesh::GlobalID getGlobalID(const vmesh::LocalID localID) const;
      ARCH_HOSTDEV vmesh::GlobalID getGlobalID(const Real* coords) const;
      ARCH_HOSTDEV vmesh::GlobalID getGlobalID(const vmesh::LocalID indices[3]) const;
      ARCH_HOSTDEV vmesh::GlobalID getGlobalID(const vmesh::LocalID i,const vmesh::LocalID j,const vmesh::LocalID k) const;
      ARCH_HOSTDEV split::SplitVector<vmesh::GlobalID>* getGrid();
      ARCH_HOSTDEV const vmesh::LocalID* getGridLength() const;
      ARCH_HOSTDEV void getIndices(const vmesh::GlobalID globalID,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const;
      ARCH_HOSTDEV void getIndicesX(const vmesh::GlobalID globalID,vmesh::LocalID& i) const;
      ARCH_HOSTDEV void getIndicesY(const vmesh::GlobalID globalID,vmesh::LocalID& j) const;
      ARCH_HOSTDEV void getIndicesZ(const vmesh::GlobalID globalID,vmesh::LocalID& k) const;
      ARCH_HOSTDEV size_t getMesh() const;
      ARCH_HOSTDEV vmesh::LocalID getLocalID(const vmesh::GlobalID globalID) const;
      ARCH_DEV vmesh::LocalID warpGetLocalID(const vmesh::GlobalID globalID, const size_t b_tid) const;
      ARCH_HOSTDEV const Real* getMeshMaxLimits() const;
      ARCH_HOSTDEV const Real* getMeshMinLimits() const;
      ARCH_HOSTDEV bool initialize(const size_t meshID);
      ARCH_HOSTDEV static vmesh::LocalID invalidBlockIndex();
      ARCH_HOSTDEV static vmesh::GlobalID invalidGlobalID();
      ARCH_HOSTDEV static vmesh::LocalID invalidLocalID();
      ARCH_HOSTDEV bool isInitialized() const;
      ARCH_HOSTDEV void pop();
      ARCH_DEV void warpPop(const size_t b_tid);
      ARCH_HOSTDEV bool push_back(const vmesh::GlobalID globalID);
      ARCH_DEV bool warpPush_back(const vmesh::GlobalID globalID, const size_t b_tid);
      vmesh::LocalID push_back(const std::vector<vmesh::GlobalID>& blocks);
      ARCH_HOSTDEV vmesh::LocalID push_back(split::SplitVector<vmesh::GlobalID>* blocks);
      ARCH_DEV vmesh::LocalID warpPush_back(const split::SplitVector<vmesh::GlobalID>& blocks, const size_t b_tid);
      ARCH_DEV void replaceBlock(const vmesh::GlobalID GIDold,const vmesh::LocalID LID,const vmesh::GlobalID GIDnew);
      ARCH_DEV void warpReplaceBlock(const vmesh::GlobalID GIDold,const vmesh::LocalID LID,const vmesh::GlobalID GIDnew, const size_t b_tid);
      ARCH_DEV void placeBlock(const vmesh::GlobalID GID,const vmesh::LocalID LID);
      ARCH_DEV void warpPlaceBlock(const vmesh::GlobalID GID,const vmesh::LocalID LID, const size_t b_tid);
      ARCH_DEV void deleteBlock(const vmesh::GlobalID GID,const vmesh::LocalID LID);
      ARCH_DEV void warpDeleteBlock(const vmesh::GlobalID GID,const vmesh::LocalID LID, const size_t b_tid);
      void setGrid();
      bool setGrid(const std::vector<vmesh::GlobalID>& globalIDs);
      bool setGrid(const split::SplitVector<vmesh::GlobalID>& globalIDs);
      bool setMesh(const size_t meshID);
      size_t getMeshID();
      void setNewSize(const vmesh::LocalID newSize);
      // bool setNewSizeClear(const vmesh::LocalID& newSize, gpuStream_t stream);
      ARCH_DEV void device_setNewSize(const vmesh::LocalID newSize);
      void setNewCachedSize(const vmesh::LocalID newSize);
      void updateCachedSize();
      void updateCachedCapacity();
      bool setNewCapacity(const vmesh::LocalID newCapacity);
      ARCH_HOSTDEV size_t size() const;
      ARCH_HOSTDEV size_t sizeInBytes() const;
      // void swap(VelocityMesh& vm);

      void gpu_prefetchHost(gpuStream_t stream);
      void gpu_prefetchDevice(gpuStream_t stream);
      void gpu_cleanHashMap(gpuStream_t stream);
      void print_addresses();
      void print_sizes();

      ARCH_HOSTDEV Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* gpu_expose_map();

   private:
      size_t meshID;
      size_t ltg_size=0, ltg_capacity=0, gtl_sizepower=0; // host-cached values

      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> globalToLocalMap;
      split::SplitVector<vmesh::GlobalID> localToGlobalMap;
   };

   // ***** DEFINITIONS OF MEMBER FUNCTIONS ***** //


   inline VelocityMesh::VelocityMesh() {
      meshID = std::numeric_limits<size_t>::max();
      // Set sizepower to 10 (1024 blocks) straight away so there's enough room to grow?
      globalToLocalMap = Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(INIT_MAP_SIZE);
      localToGlobalMap = split::SplitVector<vmesh::GlobalID>(INIT_VMESH_SIZE);
      localToGlobalMap.clear();
      ltg_size = 0;
      ltg_capacity = INIT_VMESH_SIZE;
      gtl_sizepower = INIT_MAP_SIZE;
   }

   inline VelocityMesh::~VelocityMesh() {}

   inline VelocityMesh::VelocityMesh(const VelocityMesh& other) {
      gpuStream_t stream = gpu_getStream();
      meshID = other.meshID;
      if (other.localToGlobalMap.size() > 0) {
         globalToLocalMap = Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(other.globalToLocalMap);
         localToGlobalMap = split::SplitVector<vmesh::GlobalID>(other.localToGlobalMap.capacity());
         // Overwrite is like a copy assign but takes a stream
         localToGlobalMap.overwrite(other.localToGlobalMap,stream);
         ltg_size = other.ltg_size;
         ltg_capacity = localToGlobalMap.capacity();
         gtl_sizepower = globalToLocalMap.getSizePower();
      } else {
         globalToLocalMap = Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(INIT_MAP_SIZE);
         localToGlobalMap = split::SplitVector<vmesh::GlobalID>(INIT_VMESH_SIZE);
         localToGlobalMap.clear();
         ltg_size = 0;
         ltg_capacity = INIT_VMESH_SIZE;
         gtl_sizepower = INIT_MAP_SIZE;
      }
   }

   inline const VelocityMesh& VelocityMesh::operator=(const VelocityMesh& other) {
      gpuStream_t stream = gpu_getStream();
      meshID = other.meshID;
      setNewCapacity(other.ltg_size);
      // Overwrite is like a copy assign but takes a stream
      globalToLocalMap.overwrite(other.globalToLocalMap,stream);
      // This constructor is used e.g. for boundary cells, where we in fact
      // don't need any extra capacity, so let's not call this reserve.
      // localToGlobalMap->reserve((other.localToGlobalMap)->capacity(),true);
      localToGlobalMap.overwrite(other.localToGlobalMap,stream);
      ltg_size = other.ltg_size;
      return *this;
   }

   inline size_t VelocityMesh::capacityInBytes() const {
      #ifdef DEBUG_VMESH
      const size_t cap1 = localToGlobalMap.capacity();
      if (ltg_capacity != cap1) {
         printf("VMESH CAPACITY ERROR: LTG capacity %lu vs cached value %lu in %s : %d\n",cap1,ltg_capacity,__FILE__,__LINE__);
      }
      const size_t cap2 = globalToLocalMap.getSizePower();
      if (gtl_sizepower != cap2) {
         printf("VMESH CAPACITY ERROR: GTL sizePower %lu vs cached value %lu in %s : %d\n",cap2,gtl_sizepower,__FILE__,__LINE__);
      }
      #endif
      // *** Using cached values
      const size_t currentCapacity =  ltg_capacity;
      const size_t currentBucketCount = std::pow(2,gtl_sizepower);

      const size_t capacityInBytes = currentCapacity * sizeof(vmesh::GlobalID)
           + currentBucketCount * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      return capacityInBytes;
   }

   ARCH_HOSTDEV inline size_t VelocityMesh::capacity() const {
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      return localToGlobalMap.capacity();
      #else
      #ifdef DEBUG_VMESH
      const size_t cap1 = localToGlobalMap.capacity();
      if (ltg_capacity != cap1) {
         printf("VMESH CAPACITY ERROR: capacity %lu vs cached value %lu in %s : %d\n",cap1,ltg_capacity,__FILE__,__LINE__);
      }
      #endif
      return ltg_capacity;
      #endif
   }

   ARCH_HOSTDEV inline bool VelocityMesh::check() {
      bool ok = true;
      #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Prefetch to host for the in-detail check
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap.optimizeCPU(stream);
      globalToLocalMap.optimizeCPU(stream);
      //CHK_ERR( gpuStreamSynchronize(stream) );
      CHK_ERR( gpuDeviceSynchronize() );
      #endif
      const size_t size1 = localToGlobalMap.size();
      const size_t size2 = globalToLocalMap.size();
      #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Only check cached values on host
      const size_t cap1 = localToGlobalMap.capacity();
      const size_t cap2 = globalToLocalMap.getSizePower();
      if (cap1 != ltg_capacity) {
         printf("VMESH CHECK ERROR: capacity %lu vs cached value %lu in %s : %d\n",cap1,ltg_capacity,__FILE__,__LINE__);
         return false;
      }
      if (cap2 != gtl_sizepower) {
         printf("VMESH CHECK ERROR: sizepower %lu vs cached value %lu in %s : %d\n",cap2,gtl_sizepower,__FILE__,__LINE__);
         return false;
      }
      if (size1 != ltg_size) {
         printf("VMESH CHECK ERROR: size %lu vs cached value %lu in %s : %d\n",size1,ltg_size,__FILE__,__LINE__);
         return false;
      }
      #endif
      if (size1 != size2) {
         printf("VMESH CHECK ERROR: LTG size %lu vs GTL size %lu in %s : %d\n",size1,size2,__FILE__,__LINE__);
         return false;
         //assert(0 && "VM check ERROR: sizes differ");
      }
      size_t fail = 0;
      for (size_t b=0; b<size1; ++b) {
         const vmesh::LocalID globalID = localToGlobalMap.at(b);
         #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
         auto it = globalToLocalMap.device_find(globalID);
         if (it != globalToLocalMap.device_end()) {
         #else
         auto it = globalToLocalMap.find(globalID);
         if (it != globalToLocalMap.end()) {
         #endif
            const vmesh::GlobalID localID = it->second;
            if (localID != b) {
               ok = false;
               printf("VMESH CHECK ERROR: localToGlobalMap[%lu] = %u but globalToLocalMap[%u] = %u\n",b,globalID,globalID,localID);
               //assert(0 && "VM check ERROR");
               fail++;
            }
         } else {
            ok = false;
            printf("VMESH CHECK ERROR: localToGlobalMap[%lu] = %u but could not find in globalToLocalMap ",b,globalID);
            fail++;
         }
      }
      if (fail>0) {
         printf("VMESH CHECK ERROR Encountered %lu failures.\n",fail);
         return false;
      }
      #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      localToGlobalMap.optimizeGPU(stream);
      globalToLocalMap.optimizeGPU(stream);
      #endif
      return ok;
   }

   ARCH_HOSTDEV inline void VelocityMesh::print() {
      #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap.optimizeCPU(stream);
      globalToLocalMap.optimizeCPU(stream);
      #endif

      if (localToGlobalMap.size() != globalToLocalMap.size()) {
         printf("VMO ERROR: sizes differ, %lu vs %lu\n",localToGlobalMap.size(),globalToLocalMap.size());
         assert(0 && "VM check ERROR: sizes differ");
      }
      vmesh::LocalID thisSize = size();
      printf("VM Size: %u \n",thisSize);
      for (vmesh::LocalID b=0; b<thisSize; ++b) {
         const vmesh::LocalID globalID = localToGlobalMap.at(b);
         #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
         auto it = globalToLocalMap.device_find(globalID);
         #else
         auto it = globalToLocalMap.find(globalID);
         #endif
         const vmesh::GlobalID localID = it->second;
         printf("vmesh LID [%6u] => GID [%6u] => [%6u]\n",b,globalID,localID);
      }
      #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      localToGlobalMap.optimizeGPU(stream);
      globalToLocalMap.optimizeGPU(stream);
      #endif
   }

   inline void VelocityMesh::clear(bool shrink=true) {
      ltg_size = 0;
      if (shrink) {
         ltg_capacity = 1;
         gtl_sizepower = 4;
         globalToLocalMap = Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(ltg_capacity);
         localToGlobalMap = split::SplitVector<vmesh::GlobalID>(gtl_sizepower);
         localToGlobalMap.clear();
      } else {
         gpuStream_t stream = gpu_getStream();
         localToGlobalMap.clear();
         globalToLocalMap.clear<false>(Hashinator::targets::device,stream, std::pow(2,gtl_sizepower));
         CHK_ERR( gpuStreamSynchronize(stream) );
      }

      #ifdef DEBUG_VMESH
      if ((localToGlobalMap.size() != 0) || (globalToLocalMap.size() != 0)) {
         std::cerr<<"VMESH CLEAR FAILED"<<std::endl;
      }
      #endif
   }

   ARCH_HOSTDEV inline bool VelocityMesh::move(const vmesh::LocalID sourceLID,const vmesh::LocalID targetLID) {
      const vmesh::GlobalID moveGID = localToGlobalMap.at(sourceLID); // block to move (at the end of list)
      const vmesh::GlobalID removeGID = localToGlobalMap.at(targetLID); // removed block
      #ifdef DEBUG_VMESH
      if (sourceLID != size()-1) {
         assert( 0 && "Error! Moving velocity mesh entry from position which is not last LID!");
      }
      #endif

      // at-function will throw out_of_range exception for non-existing global ID:
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      globalToLocalMap.set_element(moveGID,targetLID);
      globalToLocalMap.device_erase(removeGID);
      #else
      globalToLocalMap.at(moveGID) = targetLID;
      globalToLocalMap.erase(removeGID);
      #endif
      localToGlobalMap.at(targetLID) = moveGID;
      localToGlobalMap.pop_back();
      // Update cached value
      ltg_size--;
      return true;
   }
   ARCH_HOSTDEV inline size_t VelocityMesh::count(const vmesh::GlobalID globalID) const {
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      return globalToLocalMap.device_count(globalID);
      #else
      return globalToLocalMap.count(globalID);
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
      if (globalToLocalMap.device_count(blockGID) != 0) {
      #else
      if (globalToLocalMap.count(blockGID) != 0) {
      #endif
         return blockGID;
      } else {
         return invalidGlobalID();
      }
   }

   ARCH_HOSTDEV inline bool VelocityMesh::getBlockCoordinates(const vmesh::GlobalID globalID,Real coords[3]) const {
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
   ARCH_HOSTDEV inline void VelocityMesh::getBlockInfo(const vmesh::GlobalID globalID, Real* array) const {
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

   ARCH_HOSTDEV inline bool VelocityMesh::getBlockSize(const vmesh::GlobalID globalID, Real size[3]) const {
      size[0] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[0];
      size[1] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[1];
      size[2] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].blockSize[2];
      return true;
   }

   ARCH_HOSTDEV inline const Real* VelocityMesh::getCellSize() const {
      return (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize;
   }

   ARCH_HOSTDEV inline bool VelocityMesh::getCellSize(const vmesh::GlobalID globalID, Real size[3]) const {
      size[0] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize[0];
      size[1] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize[1];
      size[2] = (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].cellSize[2];
      return true;
   }

   ARCH_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID localID) const {
      #ifdef DEBUG_VMESH
      if (localID >= localToGlobalMap.size()) {
         assert (0 && "ERROR invalid local id");
      }
      #endif

      return localToGlobalMap.at(localID);
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

   ARCH_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID indices[3]) const {
      if (indices[0] >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]) {
         return invalidGlobalID();
      }
      if (indices[1] >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]) {
         return invalidGlobalID();
      }
      if (indices[2] >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[2]) {
         return invalidGlobalID();
      }
      return indices[0] + indices[1] * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]
         + indices[2] * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]
         * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0];
   }

   ARCH_HOSTDEV inline vmesh::GlobalID VelocityMesh::getGlobalID(const vmesh::LocalID i,const vmesh::LocalID j,const vmesh::LocalID k) const {
      if (i >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]) {
         return invalidGlobalID();
      }
      if (j >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]) {
         return invalidGlobalID();
      }
      if (k >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[2]) {
         return invalidGlobalID();
      }
      return i + j * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]
         + k * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]
         * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0];
   }

   ARCH_HOSTDEV inline split::SplitVector<vmesh::GlobalID>* VelocityMesh::getGrid() {
      return &localToGlobalMap;
   }

   ARCH_HOSTDEV inline const vmesh::LocalID* VelocityMesh::getGridLength() const {
      return (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength;
   }

   ARCH_HOSTDEV inline void VelocityMesh::getIndices(const vmesh::GlobalID globalID,vmesh::LocalID& i,vmesh::LocalID& j,vmesh::LocalID& k) const {
      if (globalID >= invalidGlobalID()) {
         i = j = k = invalidBlockIndex();
      } else {
         i = globalID % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0];
         j = (globalID / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]) % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1];
         k = globalID / ((*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0] * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]);
      }
   }

   ARCH_HOSTDEV inline void VelocityMesh::getIndicesX(const vmesh::GlobalID globalID,vmesh::LocalID& i) const {
      if (globalID >= invalidGlobalID()) {
         i = invalidBlockIndex();
      } else {
         i = globalID % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0];
      }
   }
   ARCH_HOSTDEV inline void VelocityMesh::getIndicesY(const vmesh::GlobalID globalID,vmesh::LocalID& j) const {
      if (globalID >= invalidGlobalID()) {
         j = invalidBlockIndex();
      } else {
         j = (globalID / (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0]) % (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1];
      }
   }
   ARCH_HOSTDEV inline void VelocityMesh::getIndicesZ(const vmesh::GlobalID globalID,vmesh::LocalID& k) const {
      if (globalID >= invalidGlobalID()) {
         k = invalidBlockIndex();
      } else {
         k = globalID / ((*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[0] * (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].gridLength[1]);
      }
   }

   ARCH_HOSTDEV inline vmesh::LocalID VelocityMesh::getLocalID(const vmesh::GlobalID globalID) const {
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      auto it = globalToLocalMap.device_find(globalID);
      if (it != globalToLocalMap.device_end()) {
         return it->second;
      }
      #else
      auto it = globalToLocalMap.find(globalID);
      if (it != globalToLocalMap.end()) {
         return it->second;
      }
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

   ARCH_HOSTDEV inline bool VelocityMesh::initialize(const size_t meshID) {
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
      if (mySize == 0) {
         return;
      }
      const vmesh::LocalID lastLID = mySize-1;

      #ifdef DEBUG_VMESH
      const vmesh::GlobalID lastGID = localToGlobalMap.at(lastLID);
      #else
      const vmesh::GlobalID lastGID = localToGlobalMap[lastLID];
      #endif
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      auto last = globalToLocalMap.device_find(lastGID);
      globalToLocalMap.device_erase(last);
      #else
      auto last = globalToLocalMap.find(lastGID);
      globalToLocalMap.erase(last);
      #endif
      localToGlobalMap.pop_back();
      // Update cached value
      ltg_size--;
   }
   /** Returns true if added velocity block was new, false if it couldn't be added or already existed.
    */
   ARCH_HOSTDEV inline bool VelocityMesh::push_back(const vmesh::GlobalID globalID) {
      const size_t mySize = size();
      if (mySize >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) {
         return false;
      }
      if (globalID == invalidGlobalID()) {
         return false;
      }

      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      // device_insert is slower, returns iterator and true or false for whether inserted key was new
      const bool newEntry = globalToLocalMap.set_element<true>(globalID,(vmesh::LocalID)mySize);
      if (newEntry) {
         localToGlobalMap.device_push_back(globalID);
         ltg_size++; // Note: called from inside kernel, cached size must be updated separately
         //ltg_capacity = localToGlobalMap.capacity(); // on-device, no recapacitate
      }
      return newEntry;
      #else
      auto position
         = globalToLocalMap.insert(Hashinator::make_pair(globalID,(vmesh::LocalID)mySize));
      if (position.second == true) {
         localToGlobalMap.push_back(globalID); // May increase capacity
         ltg_size++;
         if (ltg_size > ltg_capacity) {
            ltg_capacity = localToGlobalMap.capacity();
         }
         gtl_sizepower = globalToLocalMap.getSizePower();
      }
      return position.second;
      #endif
   }
   inline vmesh::LocalID VelocityMesh::push_back(const std::vector<vmesh::GlobalID>& blocks) {
      gpuStream_t stream = gpu_getStream();
      const size_t blocksSize = blocks.size();
      if (ltg_size+blocksSize > (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) {
         printf("vmesh: too many blocks, current size is %lu",ltg_size);
         printf(", adding %lu blocks", blocksSize);
         printf(", max is %u\n",(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks);
         return false;
      }

      if (ltg_size==0) {
         // Fast insertion into empty mesh
         if (blocksSize > ltg_capacity) {
            ltg_capacity = blocksSize*BLOCK_ALLOCATION_FACTOR;
            localToGlobalMap.reserve(ltg_capacity,true,stream);
         }
         localToGlobalMap.insert(localToGlobalMap.end(),blocks.begin(),blocks.end());
         vmesh::GlobalID* _localToGlobalMapData = localToGlobalMap.data();
         localToGlobalMap.optimizeGPU(stream);
         globalToLocalMap.insertIndex<false>(_localToGlobalMapData,blocksSize,0.5,stream);
         ltg_size = blocksSize;
         gtl_sizepower = globalToLocalMap.getSizePower();
         return blocksSize;
      } else {
         // GPUTODO: do inside kernel? This function is used by SpatialCell::add_velocity_blocks.
         if (ltg_size+blocksSize > ltg_capacity) {
            ltg_capacity = (ltg_size+blocksSize)*BLOCK_ALLOCATION_FACTOR;
            localToGlobalMap.reserve(ltg_capacity,true,stream);
         }
         localToGlobalMap.resize(ltg_size+blocksSize,true,stream);
         size_t newElements = 0;
         for (size_t b=0; b<blocksSize; ++b) {
            auto position
               = globalToLocalMap.insert(Hashinator::make_pair(blocks[b],(vmesh::LocalID)(ltg_size+b)));
            // Verify insertion into map and update vector
            if (position.second) { // this is true if the element did not previously exist in the map
               localToGlobalMap.at(ltg_size+newElements) = blocks[b];
               newElements++;
            }
         }
         localToGlobalMap.resize(ltg_size+newElements,true,stream);
         ltg_size += newElements;
         gtl_sizepower = globalToLocalMap.getSizePower();
         localToGlobalMap.optimizeGPU(stream);
         return newElements;
      }
   }

   ARCH_HOSTDEV inline vmesh::LocalID VelocityMesh::push_back(split::SplitVector<vmesh::GlobalID>* blocks) {
      #if !(defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap.optimizeCPU(stream); // insert one-by-one on CPU
      #endif
      const size_t blocksSize = blocks->size();
      if (ltg_size+blocksSize > (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) {
         printf("vmesh: too many blocks, current size is %lu",ltg_size);
         printf(", adding %lu blocks", blocksSize);
         printf(", max is %u\n",(*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks);
         return false;
      }

      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      localToGlobalMap.device_resize(ltg_size+blocksSize, false); //construct=false don't construct or set to zero
      size_t newElements = 0;
      for (size_t b=0; b<blocksSize; ++b) {
         // device_insert is slower than set_element, returns iterator and true or false for whether inserted key was new
         const bool newEntry = globalToLocalMap.set_element<true>((*blocks)[b],(vmesh::LocalID)(ltg_size+b));
         // Verify insertion into map and update vector
         if (newEntry) { // this is true if the element did not previously exist in the map
            localToGlobalMap[ltg_size+newElements] = (*blocks)[b];
            newElements++;
         }
      }
      localToGlobalMap.device_resize(ltg_size+newElements); //only make smaller so no construct
      ltg_size += newElements; // Note: called from inside kernel, cached size must be updated separately
      //ltg_capacity = localToGlobalMap.capacity(); // on-device, no recapacitate
      return newElements;
      #else
      if (ltg_size==0) {
         // Fast insertion into empty mesh
         if (blocksSize > ltg_capacity) {
            ltg_capacity = blocksSize*BLOCK_ALLOCATION_FACTOR;
            localToGlobalMap.reserve(ltg_capacity,true,stream);
         }
         localToGlobalMap.insert(localToGlobalMap.end(),blocks->begin(),blocks->end());
         vmesh::GlobalID* _localToGlobalMapData = localToGlobalMap.data();
         localToGlobalMap.optimizeGPU(stream);
         globalToLocalMap.insertIndex<false>(_localToGlobalMapData,blocksSize,0.5,stream);
         ltg_size = blocksSize;
         ltg_capacity = localToGlobalMap.capacity();
         gtl_sizepower = globalToLocalMap.getSizePower();
         return blocksSize;
      } else {
         // GPUTODO: do inside kernel?
         if (ltg_size+blocksSize > ltg_capacity) {
            ltg_capacity = (ltg_size+blocksSize)*BLOCK_ALLOCATION_FACTOR;
            localToGlobalMap.reserve(ltg_capacity,true,stream);
         }
         localToGlobalMap.resize(ltg_size+blocksSize,true,stream);
         size_t newElements = 0;
         for (size_t b=0; b<blocksSize; ++b) {
            auto position
               = globalToLocalMap.insert(Hashinator::make_pair((*blocks)[b],(vmesh::LocalID)(ltg_size+b)));
            // Verify insertion into map and update vector
            if (position.second) { // this is true if the element did not previously exist in the map
               localToGlobalMap.at(ltg_size+newElements) = (*blocks)[b];
               newElements++;
            }
         }
         localToGlobalMap.resize(ltg_size+newElements,true,stream);
         localToGlobalMap.optimizeGPU(stream);
         ltg_size += newElements;
         ltg_capacity = localToGlobalMap.capacity();
         gtl_sizepower = globalToLocalMap.getSizePower();
         return newElements;
      }
      #endif
   }

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
   /****
         Device-only accessors to be called from a single GPU thread
    **/
   ARCH_DEV inline void VelocityMesh::replaceBlock(const vmesh::GlobalID GIDold,const vmesh::LocalID LID,const vmesh::GlobalID GIDnew) {
      #ifdef DEBUG_VMESH
      if (LID > size()-1) {
         printf("vmesh replaceBlock error: LID is too large!\n");
      }
      vmesh::LocalID LIDold = invalidLocalID();
      auto it = globalToLocalMap.device_find(GIDold);
      if (it != globalToLocalMap.device_end()) {
         LIDold = it->second;
      }
      if (localToGlobalMap.at(LIDold) != GIDold) {
         printf("vmesh replaceBlock error: oldGID and oldLID don't match!\n");
      }
      #endif
      globalToLocalMap.device_erase(GIDold);
      globalToLocalMap.set_element(GIDnew,LID);
      localToGlobalMap.at(LID) = GIDnew;
   }
   // Note: this function does not adjust the vmesh size, as it is used from within a parallel kernel.
   ARCH_DEV inline void VelocityMesh::placeBlock(const vmesh::GlobalID GID,const vmesh::LocalID LID) {
      #ifdef DEBUG_VMESH
      if (LID > size()-1) {
         assert(0 && "vmesh placeBlock error: LID is too large!");
      }
      bool newEntry = globalToLocalMap.set_element(GID,LID);
      if (!newEntry) {
         assert(0 && "vmesh placeBlock error: set_element not a new entry!");
      }
      localToGlobalMap.at(LID) = GID;
      #else
      globalToLocalMap.set_element(GID,LID);
      localToGlobalMap[LID] = GID;
      #endif
   }
   ARCH_DEV inline void VelocityMesh::deleteBlock(const vmesh::GlobalID GID,const vmesh::LocalID LID) {
      #ifdef DEBUG_VMESH
      // Verify that GID and LID match
      if (GID==invalidGlobalID()) {
         printf("vmesh deleteBlock error: GID is invalidGlobalID! (LID %ul)\n",LID);
      }
      if (LID==invalidLocalID()) {
         printf("vmesh deleteBlock error: LID is invalidLocalID! (GID %ul)\n",GID);
      }
      auto it = globalToLocalMap.device_find(GID);
      if (it == globalToLocalMap.device_end()) {
         printf("vmesh deleteBlock error: GID %ul does not exist! (LID %ul)\n",GID,LID);
      } else {
         if (it->second != LID) {
            printf("vmesh deleteBlock error: LID %ul found with GID %ul does not match provided LID %ul!\n",it->second,GID,LID);
         }
      }
      if (localToGlobalMap.at(LID) != GID) {
         printf("vmesh deleteBlock error: GID %ul found with LID %ul does not match provided GID %ul!\n",localToGlobalMap.at(LID),LID,GID);
      }
      #endif
      globalToLocalMap.device_erase(GID);
      #ifdef DEBUG_VMESH
      localToGlobalMap.at(LID) = invalidGlobalID();
      #else
      localToGlobalMap[LID] = invalidGlobalID();
      #endif
      ltg_size--;
   }

   /****
        Warp accessor functions to be called from within GPU kernels over several threads
   **/
   ARCH_DEV inline void VelocityMesh::warpPop(const size_t b_tid) {
      const size_t mySize = size();
      if (mySize == 0) {
         return;
      }
      const vmesh::LocalID lastLID = mySize-1;
      #ifdef DEBUG_VMESH
      const vmesh::GlobalID lastGID = localToGlobalMap.at(lastLID);
      const vmesh::LocalID mapSize = globalToLocalMap.size();
      #else
      const vmesh::GlobalID lastGID = localToGlobalMap[lastLID];
      #endif
      if (b_tid < GPUTHREADS) {
         globalToLocalMap.warpErase(lastGID, b_tid);
      }
      if (b_tid==0) {
         localToGlobalMap.pop_back();
      }
      #ifdef DEBUG_VMESH
      __syncthreads();
      const vmesh::LocalID newMapSize = globalToLocalMap.size();
      const vmesh::LocalID newVecSize = localToGlobalMap.size();
      if (newMapSize != mapSize-1) {
         printf("warpError in VelocityMesh::warpPop: map size %u is not expected %u! (thread %u)\n",newMapSize,(vmesh::LocalID)(mapSize-1),(vmesh::LocalID)b_tid);
         assert(0);
      }
      if (newVecSize != mySize-1) {
         printf("warpError in VelocityMesh::warpPop: vector size %u is not expected %u! (thread %u)\n",newVecSize,(vmesh::LocalID)(mySize-1),(vmesh::LocalID)b_tid);
         assert(0);
      }
      #endif
      if (b_tid==0) {
         ltg_size--;
      }
      __syncthreads();
   }
   ARCH_DEV inline vmesh::LocalID VelocityMesh::warpGetLocalID(const vmesh::GlobalID globalID, const size_t b_tid) const {
      vmesh::LocalID retval = invalidLocalID();
      globalToLocalMap.warpFind(globalID, retval, b_tid % GPUTHREADS);
      #ifdef DEBUG_VMESH
      auto it = globalToLocalMap.device_find(globalID);
      if (it == globalToLocalMap.device_end()) {
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
   ARCH_DEV inline bool VelocityMesh::warpMove(const vmesh::LocalID sourceLID,const vmesh::LocalID targetLID, const size_t b_tid) {
      #ifdef DEBUG_VMESH
      const vmesh::GlobalID moveGID = localToGlobalMap.at(sourceLID); // block to move (must at the end of list)
      const vmesh::GlobalID removeGID = localToGlobalMap.at(targetLID); // removed block
      if (sourceLID != size()-1) {
         assert( 0 && "Error! Moving velocity mesh entry from position which is not last LID!");
      }
      const vmesh::LocalID preMapSize = globalToLocalMap.size();
      const vmesh::LocalID preVecSize = localToGlobalMap.size();
      #else
      const vmesh::GlobalID moveGID = localToGlobalMap[sourceLID]; // block to move (must at the end of list)
      const vmesh::GlobalID removeGID = localToGlobalMap[targetLID]; // removed block
      #endif

      // at-function will throw out_of_range exception for non-existing global ID:
      if (b_tid < GPUTHREADS) {
         globalToLocalMap.warpInsert(moveGID,targetLID,b_tid); // will overwrite
         globalToLocalMap.warpErase(removeGID,b_tid);
      }
      if (b_tid==0) {
         localToGlobalMap.at(targetLID) = moveGID;
         localToGlobalMap.pop_back();
      }
      #ifdef DEBUG_VMESH
      __syncthreads();
      const vmesh::LocalID postMapSize = globalToLocalMap.size();
      const vmesh::LocalID postVecSize = localToGlobalMap.size();
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
      if (b_tid==0) {
         ltg_size--;
      }
      __syncthreads();
      return true;
   }
   ARCH_DEV inline size_t VelocityMesh::warpCount(const vmesh::GlobalID globalID, const size_t b_tid) const {
      vmesh::LocalID retval = invalidLocalID();
      globalToLocalMap.warpFind(globalID, retval, b_tid % GPUTHREADS);
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
      globalToLocalMap.warpFind(blockGID, retval, b_tid % GPUTHREADS);

      #ifdef DEBUG_VMESH
      __syncthreads();
      auto it = globalToLocalMap.device_find(blockGID);
      if (it == globalToLocalMap.device_end()) {
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
   ARCH_DEV inline bool VelocityMesh::warpPush_back(const vmesh::GlobalID globalID, const size_t b_tid) {
      __shared__ bool inserted;
      const vmesh::LocalID mySize = size();
      #ifdef DEBUG_VMESH
      if (mySize >= (*(vmesh::getMeshWrapper()->velocityMeshes))[meshID].max_velocity_blocks) {
         return false;
      }
      if (globalID == invalidGlobalID()) {
         return false;
      }
      const vmesh::LocalID mapSize = globalToLocalMap.size();
      #endif

      __syncthreads();
      if (b_tid < GPUTHREADS) {
         // If exists, do not overwrite
         inserted = globalToLocalMap.warpInsert_V<true>(globalID,(vmesh::LocalID)mySize, b_tid);
         if (inserted == true && b_tid==0) {
            localToGlobalMap.device_push_back(globalID);
            ltg_size++; // Note: called from inside kernel, cached size must be updated separately
         }
      }
      #ifdef DEBUG_VMESH
      __syncthreads();
      const vmesh::LocalID newMapSize = globalToLocalMap.size();
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
            // Do not overwrite
            inserted = inserted && globalToLocalMap.warpInsert_V<true>(blocks[b],(vmesh::LocalID)(mySize+b), b_tid);
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
         localToGlobalMap.device_insert(localToGlobalMap.end(),blocks.begin(),blocks.end());
      }
      if (b_tid==0) {
         ltg_size += nInserted; // Note: called from inside kernel, cached size must be updated separately
      }
      __syncthreads();
      return blocksSize;
   }
   ARCH_DEV inline void VelocityMesh::warpReplaceBlock(const vmesh::GlobalID GIDold,
                                                       const vmesh::LocalID LID,
                                                       const vmesh::GlobalID GIDnew,
                                                       const size_t b_tid) {
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
      globalToLocalMap.warpFind(GIDold, LIDold, b_tid % GPUTHREADS);
      auto it = globalToLocalMap.device_find(GIDold);
      if (it == globalToLocalMap.device_end()) {
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
         if (localToGlobalMap.at(LIDold) != GIDold) {
            printf("vmesh replaceBlock error: oldGID and oldLID don't match!\n");
            assert(0);
         }
      }
      __syncthreads();
      const vmesh::LocalID preVecSize = size();
      const vmesh::LocalID preMapSize = globalToLocalMap.size();
      #endif
      if (b_tid < GPUTHREADS) { // GPUTODO these in parallel?
         globalToLocalMap.warpErase(GIDold, b_tid);

         #ifdef DEBUG_VMESH
         // if (globalToLocalMap.size() != preMapSize-1) {
         //    printf("Warp error in VelocityMesh::warpReplaceBlock: map size %u does not match expected %u for thread %u!\n",(vmesh::LocalID)globalToLocalMap.size(),(vmesh::LocalID)(preMapSize-1),(vmesh::LocalID)b_tid);
         //    assert(0);
         // }
         auto it2 = globalToLocalMap.device_find(GIDold);
         if (it2 != globalToLocalMap.device_end()) {
            printf("Warp error in VelocityMesh::warpReplaceBlock: warp-erased GID %u LID %u but thread %u still finds LID %u associated with it!\n",GIDold,LID,(vmesh::LocalID)b_tid,it2->second);
            assert(0);
         }
         bool newlyadded = false;
         // Do not overwrite
         newlyadded = globalToLocalMap.warpInsert_V(GIDnew,LID, b_tid);
         vmesh::LocalID postMapSize = globalToLocalMap.size();
         // int change = 0;
         // if (!newlyadded) {
         //    change = -1;
         // }
         // if (postMapSize != preMapSize + change) {
         //    printf("Warp error in VelocityMesh::warpReplaceBlock: map size %u does not match expected %u for thread %u!\n",postMapSize,(vmesh::LocalID)(preMapSize+change),(vmesh::LocalID)b_tid);
         //    if (b_tid==0) {
         //       globalToLocalMap.stats();
         //    }
         //    assert(0);
         // }
         auto it3 = globalToLocalMap.device_find(GIDnew);
         if (it3 == globalToLocalMap.device_end()) {
            int sizePower = globalToLocalMap.getSizePower();
            printf("Warp error in VelocityMesh::warpReplaceBlock: warp-inserted GID %u LID %u but thread %u cannot find it!\n",GIDnew,LID,(vmesh::LocalID)b_tid);
            if (b_tid==0) {
               if (newlyadded) {
                  printf("warpAccessor reported true for insertion!\n");
               } else {
                  printf("warpAccessor reported false for insertion!\n");
               }
               globalToLocalMap.stats();
               //globalToLocalMap.dump_buckets();
            }
            assert(0);
         } else if (it3->second != LID) {
            printf("Warp error in VelocityMesh::warpReplaceBlock: warp-inserted GID %u LID %u but thread %u instead finds LID %u!\n",GIDnew,LID,(vmesh::LocalID)b_tid,it3->second);
            if (b_tid==0) {
               globalToLocalMap.stats();
            }
            assert(0);
         }
         #else
         bool newlyadded = globalToLocalMap.warpInsert_V(GIDnew,LID, b_tid);
         //globalToLocalMap.warpInsert(GIDnew,LID,b_tid);
         #endif
      }
      //__syncthreads(); // not needed
      if (b_tid==0) {
         localToGlobalMap.at(LID) = GIDnew;
      }
      __syncthreads();
   }
   ARCH_DEV inline void VelocityMesh::warpPlaceBlock(const vmesh::GlobalID GID,const vmesh::LocalID LID, const size_t b_tid) {
      // Places block GID into the mesh with value LID. Assumes localToGlobalMap has already been grown sufficiently.

      #ifdef DEBUG_VMESH
      if (b_tid==0) {
         if (LID > size()-1) printf("vmesh placeBlock error: LID is too large!\n");
      }
      __syncthreads();
      {
         auto it = globalToLocalMap.device_find(GID);
         if (it != globalToLocalMap.device_end()) {
            printf("Warp error in VelocityMesh::warpPlaceBlock: single-thread %u search found GID %u=%u LID %u before it was inserted!\n",(vmesh::LocalID)b_tid,GID,it->first,it->second);
            if (b_tid==0) {
               globalToLocalMap.stats();
            }
            __syncthreads();
            assert(0);
         }
      }
      __syncthreads();
      #endif
      if (b_tid < GPUTHREADS) {
         bool newlyadded = false;
         newlyadded = globalToLocalMap.warpInsert_V(GID,LID, b_tid);
         if (!newlyadded) {
            if (b_tid==0) {
               globalToLocalMap.stats();
               printf("warpPlaceBlock error GID %u LID %u reported as not newly added! Size %zu.\n",GID,LID,localToGlobalMap.size());
               //globalToLocalMap.dump_buckets();
            }
            assert(newlyadded && "newlyAdded warpPlaceBlock");
         }
         localToGlobalMap.at(LID) = GID;
      }
      #ifdef DEBUG_VMESH
      __syncthreads();
      // Note: no size check possible.
      auto it = globalToLocalMap.device_find(GID);
      if (it == globalToLocalMap.device_end()) {
         printf("Warp error in VelocityMesh::warpPlaceBlock: single-thread %u search did not find inserted GID %u LID %u\n",(vmesh::LocalID)b_tid,GID,LID);
         if (b_tid==0) {
            globalToLocalMap.stats();
         }
         assert(0);
      } else if (LID != it->second) {
         printf("Warp error in VelocityMesh::warpPlaceBlock: LID %u (warp) != %u (thread %u) for GID %u\n",LID,it->second,(vmesh::LocalID)b_tid,GID);
         if (b_tid==0) {
            globalToLocalMap.stats();
         }
         assert(0);
      }
      #endif
      __syncthreads();
   }
   ARCH_DEV inline void VelocityMesh::warpDeleteBlock(const vmesh::GlobalID GID,const vmesh::LocalID LID, const size_t b_tid) {
      #ifdef DEBUG_VMESH
      // Verify that GID and LID match
      if (b_tid==0) {
         if (GID==invalidGlobalID()) {
            printf("vmesh warpDeleteBlock error: GID is invalidGlobalID! (LID %ul)\n",LID);
         }
         if (LID==invalidLocalID()) {
            printf("vmesh warpDeleteBlock error: LID is invalidLocalID! (GID %ul)\n",GID);
         }
      }
      __syncthreads();
      vmesh::LocalID retval = invalidLocalID();
      globalToLocalMap.warpFind(GID, retval, b_tid % GPUTHREADS);
      if (b_tid==0) {
         if (retval == invalidLocalID()) {
            printf("vmesh warpDeleteBlock error: GID %ul does not exist! (LID %ul)\n",GID,LID);
         } else {
            if (retval != LID) {
               printf("vmesh warpDeleteBlock error: LID %ul warpFound with GID %ul does not match provided LID %ul!\n",
                      retval,GID,LID);
            }
         }
         if (localToGlobalMap.at(LID) != GID) {
            printf("vmesh warpDeleteBlock error: GID %ul warpFound with LID %ul does not match provided GID %ul!\n",
                   localToGlobalMap.at(LID),LID,GID);
         }
      }
      const vmesh::LocalID preMapSize = globalToLocalMap.size();
      __syncthreads();
      #endif
      if (b_tid < GPUTHREADS) {
         globalToLocalMap.warpErase(GID, b_tid);
         if (b_tid==0) {
            localToGlobalMap.at(LID) = invalidGlobalID();
         }
      }
      __syncthreads();
      #ifdef DEBUG_VMESH
      // const vmesh::LocalID postMapSize = globalToLocalMap.size();
      // if (postMapSize != preMapSize-1) {
      //    printf("Warp error in VelocityMesh::warpDeleteBlock: map size %u does not match expected %u for thread %u!\n",postMapSize,preMapSize-1,(vmesh::LocalID)b_tid);
      //    assert(0);
      // }
      auto it = globalToLocalMap.device_find(GID);
      if (it != globalToLocalMap.device_end()) {
         printf("Warp error in VelocityMesh::warpDeleteBlock: GID %u still found with LID %u for thread %u!\n",GID,it->second,(vmesh::LocalID)b_tid);
         assert(0);
      }
      __syncthreads();
      #endif
   }
#endif

   inline void VelocityMesh::setGrid() {
      // Assumes we have a valid localToGlobalMap from e.g. MPI communication,
      // populates globalToLocalMap based on it.
      gpuStream_t stream = gpu_getStream();
      globalToLocalMap.clear<false>(Hashinator::targets::device,stream,std::pow(2,gtl_sizepower));
      CHK_ERR( gpuStreamSynchronize(stream) );
      size_t nBlocks = localToGlobalMap.size();
      globalToLocalMap.insertIndex<false>(localToGlobalMap.data(),nBlocks,0.5,stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      ltg_size = nBlocks;
      ltg_capacity = localToGlobalMap.capacity();
      gtl_sizepower = globalToLocalMap.getSizePower();
   }

   // GPUTODO: These accessors are still slow, but we don't actually use them at all.
   inline bool VelocityMesh::setGrid(const std::vector<vmesh::GlobalID>& globalIDs) {
      printf("Warning! Slow version of VelocityMesh::setGrid.\n");
      gpuStream_t stream = gpu_getStream();
      globalToLocalMap.clear<false>(Hashinator::targets::device,stream,std::pow(2,gtl_sizepower));
      CHK_ERR( gpuStreamSynchronize(stream) );
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap.insert(Hashinator::make_pair(globalIDs[i],(vmesh::LocalID)i));
      }
      localToGlobalMap.clear();
      localToGlobalMap.insert(localToGlobalMap.end(),globalIDs.begin(),globalIDs.end());
      ltg_size = globalIDs.size();
      ltg_capacity = localToGlobalMap.capacity();
      gtl_sizepower = globalToLocalMap.getSizePower();
      return true;
   }
   inline bool VelocityMesh::setGrid(const split::SplitVector<vmesh::GlobalID>& globalIDs) {
      printf("Warning! Slow version of VelocityMesh::setGrid.\n");
      gpuStream_t stream = gpu_getStream();
      globalToLocalMap.clear<false>(Hashinator::targets::device,stream,std::pow(2,gtl_sizepower));
      CHK_ERR( gpuStreamSynchronize(stream) );
      for (vmesh::LocalID i=0; i<globalIDs.size(); ++i) {
         globalToLocalMap.insert(Hashinator::make_pair(globalIDs[i],(vmesh::LocalID)i));
      }
      localToGlobalMap.clear();
      localToGlobalMap.insert(localToGlobalMap.end(),globalIDs.begin(),globalIDs.end());
      ltg_size = globalIDs.size();
      ltg_capacity = localToGlobalMap.capacity();
      gtl_sizepower = globalToLocalMap.getSizePower();
      return true;
   }

   inline bool VelocityMesh::setMesh(const size_t meshID) {
      if (meshID >= vmesh::getMeshWrapper()->velocityMeshes->size()) {
         return false;
      }
      this->meshID = meshID;
      return true;
   }
   inline size_t VelocityMesh::getMeshID() {
      return this->meshID;
   }

   inline void VelocityMesh::setNewSize(const vmesh::LocalID newSize) {
      // Needed by GPU block adjustment
      setNewCapacity(newSize);
      gpuStream_t stream = gpu_getStream();
      localToGlobalMap.resize(newSize,true,stream);
      ltg_size = newSize;
   }

   ARCH_DEV inline void VelocityMesh::device_setNewSize(const vmesh::LocalID newSize) {
      if (newSize>0) {
         const vmesh::LocalID currentCapacity = localToGlobalMap.capacity();
         assert(newSize <= currentCapacity && "insufficient vector capacity in vmesh::device_setNewSize");
         const int currentSizePower = globalToLocalMap.getSizePower();
         assert(ceil(log2(newSize)) <= currentSizePower && "insufficient map capacity in vmesh::device_setNewSize");
      }
      localToGlobalMap.device_resize(newSize,false); //construct=false don't construct or set to zero
      ltg_size = newSize; // Remember to update size on host as well
   }
   inline void VelocityMesh::setNewCachedSize(const vmesh::LocalID newSize) {
      // Should only be used to update host-side size if resizing on device
      ltg_size = newSize;
   }
   inline void VelocityMesh::updateCachedSize() {
      // More secure page-faulting way to update cached size
      ltg_size = localToGlobalMap.size();
   }
   inline void VelocityMesh::updateCachedCapacity() {
      // Should not be needed, added as an optional safeguard
      ltg_capacity = localToGlobalMap.capacity();
      gtl_sizepower = globalToLocalMap.getSizePower();
   }

   // Used in initialization
   inline bool VelocityMesh::setNewCapacity(const vmesh::LocalID newCapacity) {
      // Ensure also that the map is large enough (newCapacity always greater than zero here)
      uint HashmapReqSize = (uint)ceil(log2(newCapacity)) +1;

      if (ltg_capacity >= newCapacity && gtl_sizepower >= HashmapReqSize) {
         return false; // Still have enough buffer
      }
      gpuStream_t stream = gpu_getStream();
      ltg_capacity = newCapacity*BLOCK_ALLOCATION_FACTOR;
      // Passing eco flag = true to resize tells splitvector we manage padding manually.
      localToGlobalMap.reserve(ltg_capacity,true, stream);
      if (gtl_sizepower < HashmapReqSize) {
         gtl_sizepower = HashmapReqSize;
         globalToLocalMap.resize(gtl_sizepower, Hashinator::targets::device, stream);
      }
      CHK_ERR( gpuStreamSynchronize(stream) );
      return true;
   }

   ARCH_HOSTDEV inline size_t VelocityMesh::size() const {
      #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      return localToGlobalMap.size();
      #else
      #ifdef DEBUG_VMESH
      const size_t size1 = localToGlobalMap.size();
      if (ltg_size != size1) {
         printf("VMESH SIZE ERROR: size %lu vs cached value %lu in %s : %d\n",size1,ltg_size,__FILE__,__LINE__);
      }
      #endif
      return ltg_size;
      #endif
   }

   ARCH_HOSTDEV inline size_t VelocityMesh::sizeInBytes() const {
      #ifdef DEBUG_VMESH
      const size_t size1 = localToGlobalMap.size();
      const size_t size2 = globalToLocalMap.size();
      if (ltg_size != size1) {
         printf("VMESH SIZE ERROR: LTG size %lu vs cached value %lu in %s : %d\n",size1,ltg_size,__FILE__,__LINE__);
      }
      if (ltg_size != size2) {
         printf("VMESH SIZE ERROR: GTL size %lu vs cached value %lu in %s : %d\n",size1,ltg_size,__FILE__,__LINE__);
      }
      #endif
      const size_t sizeInBytes = ltg_size * sizeof(vmesh::GlobalID)
           + ltg_size * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      return sizeInBytes;
   }

   // inline void VelocityMesh::swap(VelocityMesh& vm) {
   //    gpuStream_t stream = gpu_getStream();
   //    globalToLocalMap.swap(vm.globalToLocalMap);
   //    localToGlobalMap.swap(vm.localToGlobalMap);
   // }
   inline void VelocityMesh::gpu_prefetchHost(gpuStream_t stream=0) {
      if (stream==0) {
         stream = gpu_getStream();
      }
      localToGlobalMap.optimizeCPU(stream);
      globalToLocalMap.optimizeCPU(stream);
      return;
   }

   inline void VelocityMesh::gpu_prefetchDevice(gpuStream_t stream=0) {
      if (stream==0) {
         stream = gpu_getStream();
      }
      //phiprof::Timer vmeshPrefetchTimer {"prefetch Vmesh"};
      localToGlobalMap.optimizeGPU(stream);
      globalToLocalMap.optimizeGPU(stream);
      //CHK_ERR( gpuStreamSynchronize(stream) );
      return;
   }

   inline void VelocityMesh::gpu_cleanHashMap(gpuStream_t stream = 0) {
      if (stream==0) {
         stream = gpu_getStream();
      }

      // phiprof::Timer resizeTimer {"Hashinator resize"};
      // //globalToLocalMap.performCleanupTasks<false>(stream);
      // globalToLocalMap.resize_to_lf(0.5, Hashinator::targets::device, stream);
      // CHK_ERR( gpuStreamSynchronize(stream) );
      // resizeTimer.stop();

      phiprof::Timer cleanupTimer {"Hashinator tombstones"};
      // globalToLocalMap.clean_tombstones<false>(stream);
      globalToLocalMap.performCleanupTasks<false>(stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      cleanupTimer.stop();
      return;
   }

   ARCH_HOSTDEV inline Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* VelocityMesh::gpu_expose_map() {
      return &globalToLocalMap;
   }

   inline void VelocityMesh::print_addresses() {
      printf("GPU localToGlobalMap %p\n GPU globalToLocalMap %p\n",&localToGlobalMap,&globalToLocalMap);
      printf("GPU localToGlobalMap capacity %zu size %zu \n GPU globalToLocalMap size %zu bucket count %zu\n",localToGlobalMap.capacity(),localToGlobalMap.size(),globalToLocalMap.size(),globalToLocalMap.bucket_count());
      printf("GPU localToGlobalMap data %p\n",localToGlobalMap.data());
      //printf("GPU localToGlobalMap iterators begin %p end  %p\n",localToGlobalMap.begin(),localToGlobalMap.end());
   }
   inline void VelocityMesh::print_sizes() {
      printf("GPU localToGlobalMap size %lu capacity %lu cached size %lu cached capacity %lu\nGPU globalToLocalMap fill %lu sizePower %lu cached sizePower %lu\n",
             localToGlobalMap.size(),localToGlobalMap.capacity(),ltg_size,ltg_capacity,globalToLocalMap.size(),(size_t)globalToLocalMap.getSizePower(),gtl_sizepower);
   }

} // namespace vmesh

#endif
