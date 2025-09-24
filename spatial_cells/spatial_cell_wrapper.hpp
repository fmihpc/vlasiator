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

#ifndef VELOCITY_BLOCK_CONTAINER_H
#define VELOCITY_BLOCK_CONTAINER_H

#include <stdio.h>
#include <unistd.h>
#include <vector>

#include "../common.h"

#if defined(DEBUG_VLASIATOR) || defined(DEBUG_SPATIAL_CELL)
#ifndef DEBUG_VBC
   #define DEBUG_VBC
#endif
#endif

#ifdef DEBUG_VBC
#include <sstream>
#endif

#ifdef USE_GPU
#include "../arch/gpu_base.hpp"
   //#include "../arch/arch_device_api.h" // included in above
   // Place block data and parameters inside splitvectors utilizing unified memory
#include "include/splitvector/splitvec.h"
#else
// GPU allocation factors are stored in arch/gpu_base.hpp
static const double BLOCK_ALLOCATION_FACTOR = 1.1;
static const double BLOCK_ALLOCATION_PADDING = 1.3;
#endif

// INIT_VMESH_SIZE defined in arch/gpu_base.hpph

using namespace std;

namespace vmesh {

   class VelocityBlockContainer {
   public:
      VelocityBlockContainer();
      ~VelocityBlockContainer();
      VelocityBlockContainer(const VelocityBlockContainer& other);
      const VelocityBlockContainer& operator=(const VelocityBlockContainer& other);

      ARCH_HOSTDEV vmesh::LocalID capacity() const;
      ARCH_HOSTDEV size_t capacityInBytes() const;
      void clear(bool shrink = true);
      ARCH_HOSTDEV void move(const vmesh::LocalID source, const vmesh::LocalID target);
      ARCH_HOSTDEV static double getBlockAllocationFactor();
      ARCH_HOSTDEV Realf* getData();
      ARCH_HOSTDEV const Realf* getData() const;
      ARCH_HOSTDEV Realf* getData(const vmesh::LocalID blockLID);
      ARCH_HOSTDEV const Realf* getData(const vmesh::LocalID blockLID) const;
      ARCH_HOSTDEV Real* getParameters();
      ARCH_HOSTDEV const Real* getParameters() const;
      ARCH_HOSTDEV Real* getParameters(const vmesh::LocalID blockLID);
      ARCH_HOSTDEV const Real* getParameters(const vmesh::LocalID blockLID) const;
      ARCH_HOSTDEV void pop();
      ARCH_HOSTDEV vmesh::LocalID push_back();
      ARCH_HOSTDEV vmesh::LocalID push_back_and_zero();
      ARCH_HOSTDEV vmesh::LocalID push_back(const uint32_t N_blocks);
      ARCH_HOSTDEV vmesh::LocalID push_back_and_zero(const uint32_t N_blocks);
      ARCH_HOSTDEV void resize(const vmesh::LocalID newSize);
      ARCH_HOSTDEV bool setNewSize(const vmesh::LocalID newSize);
      bool setNewCapacityShrink(const vmesh::LocalID reqCapacity);
      ARCH_HOSTDEV vmesh::LocalID size() const;
      ARCH_HOSTDEV size_t sizeInBytes() const;
      // void swap(VelocityBlockContainer& vbc);

#ifdef USE_GPU
      // Functions for use with GPU branch
      void setNewCachedSize(const vmesh::LocalID newSize);
      void updateCachedSize();
      void updateCachedCapacity();
      bool setNewCapacity(const vmesh::LocalID capacity, gpuStream_t stream);
      void gpu_prefetchHost(gpuStream_t stream);
      void gpu_prefetchDevice(gpuStream_t stream);
      void print_addresses();
#else
      bool setNewCapacity(const vmesh::LocalID capacity);
#endif

#ifdef DEBUG_VBC
      const Realf& getData(const vmesh::LocalID blockLID, const unsigned int cell) const;
      const Real& getParameters(const vmesh::LocalID blockLID, const unsigned int i) const;
      void setData(const vmesh::LocalID blockLID, const unsigned int cell, const Realf value);
#endif

   private:
      void exitInvalidLocalID(const vmesh::LocalID localID, const std::string& funcName) const;
      ARCH_DEV void exitInvalidLocalID(const vmesh::LocalID localID) const;

#ifdef USE_GPU
      split::SplitVector<Realf> block_data;
      split::SplitVector<Real> parameters;
      size_t cachedCapacity;
      size_t cachedSize;
#else
      std::vector<Realf, aligned_allocator<Realf, WID3>> block_data;
      std::vector<Real, aligned_allocator<Real, BlockParams::N_VELOCITY_BLOCK_PARAMS>> parameters;
#endif
   };

   inline VelocityBlockContainer::VelocityBlockContainer() {
#ifdef USE_GPU
      block_data = split::SplitVector<Realf>(INIT_VMESH_SIZE * WID3);
      parameters = split::SplitVector<Real>(INIT_VMESH_SIZE * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      cachedCapacity = INIT_VMESH_SIZE;
      cachedSize = 0;
#else
      block_data = std::vector<Realf, aligned_allocator<Realf, WID3>>(WID3);
      parameters = std::vector<Real, aligned_allocator<Real, BlockParams::N_VELOCITY_BLOCK_PARAMS>>(BlockParams::N_VELOCITY_BLOCK_PARAMS);
      // cachedCapacity = 1;
#endif
      block_data.clear();
      parameters.clear();
      // gpuStream_t stream = gpu_getStream();
   }

   inline VelocityBlockContainer::~VelocityBlockContainer() {}

   inline VelocityBlockContainer::VelocityBlockContainer(const VelocityBlockContainer& other) {
#ifdef USE_GPU
      block_data = split::SplitVector<Realf>(other.cachedCapacity * WID3);
      parameters = split::SplitVector<Real>(other.cachedCapacity * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      // Overwrite is like a copy assign but takes a stream
      gpuStream_t stream = gpu_getStream();
      block_data.overwrite(other.block_data, stream);
      parameters.overwrite(other.parameters, stream);
      cachedSize = other.cachedSize;
      cachedCapacity = other.cachedCapacity;
#else
      block_data = std::vector<Realf, aligned_allocator<Realf, WID3>>(other.block_data);
      parameters = std::vector<Real, aligned_allocator<Real, BlockParams::N_VELOCITY_BLOCK_PARAMS>>(other.parameters);
      // block_data.reserve(other.capacity()*WID3);
      // parameters.reserve(other.capacity()*BlockParams::N_VELOCITY_BLOCK_PARAMS);
#endif
   }

   inline const VelocityBlockContainer& VelocityBlockContainer::operator=(const VelocityBlockContainer& other) {
#ifdef USE_GPU
      gpuStream_t stream = gpu_getStream();
      block_data.reserve(other.cachedCapacity * WID3, true, stream);
      parameters.reserve(other.cachedCapacity * BlockParams::N_VELOCITY_BLOCK_PARAMS, true, stream);
      // Overwrite is like a copy assign but takes a stream
      block_data.overwrite(other.block_data, stream);
      parameters.overwrite(other.parameters, stream);
      cachedSize = other.cachedSize;
      cachedCapacity = other.cachedCapacity;
#else
      block_data = other.block_data;
      parameters = other.parameters;
      // block_data.reserve(other.capacity()*WID3);
      // parameters.reserve(other.capacity()*BlockParams::N_VELOCITY_BLOCK_PARAMS);
#endif
      return *this;
   }

   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::capacity() const {
#ifdef USE_GPU
#ifdef DEBUG_VBC
      const size_t currentCapacity = block_data.capacity() / WID3;
      if (currentCapacity != cachedCapacity) {
         printf("VBC CHECK ERROR: capacity %lu vs cached value %lu in %s : %d\n", currentCapacity, cachedCapacity, __FILE__, __LINE__);
      }
      #endif
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      return cachedCapacity;
      #else
      return block_data.capacity() / WID3;
      #endif
#else
      return block_data.capacity() / WID3;
#endif
   }

   inline ARCH_HOSTDEV size_t VelocityBlockContainer::capacityInBytes() const {
#ifdef USE_GPU
#ifdef DEBUG_VBC
      const size_t currentCapacity = block_data.capacity() / WID3;
      if (currentCapacity != cachedCapacity) {
         printf("VBC CHECK ERROR: capacity, %lu vs cached value %lu in %s : %d\n", currentCapacity, cachedCapacity, __FILE__, __LINE__);
      }
      #endif
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      return cachedCapacity * WID3 * sizeof(Realf) + cachedCapacity * BlockParams::N_VELOCITY_BLOCK_PARAMS * sizeof(Real);
      #else
      return block_data.capacity() * sizeof(Realf) + parameters.capacity() * BlockParams::N_VELOCITY_BLOCK_PARAMS * sizeof(Real);
      #endif
#else
      return block_data.capacity() * sizeof(Realf) + parameters.capacity() * BlockParams::N_VELOCITY_BLOCK_PARAMS * sizeof(Real);
#endif
   }

   /** Clears VelocityBlockContainer data and (if shrinking) deallocates all memory
    * reserved for velocity blocks.*/
   inline void VelocityBlockContainer::clear(bool shrink) {
#ifdef USE_GPU
      cachedSize = 0;
      if (shrink) {
         cachedCapacity = 1;
         block_data = split::SplitVector<Realf>(WID3);
         parameters = split::SplitVector<Real>(BlockParams::N_VELOCITY_BLOCK_PARAMS);
      }
#else
      if (shrink) {
         block_data = std::vector<Realf, aligned_allocator<Realf, WID3>>(WID3);
         parameters = std::vector<Real, aligned_allocator<Real, BlockParams::N_VELOCITY_BLOCK_PARAMS>>(BlockParams::N_VELOCITY_BLOCK_PARAMS);
      }
#endif
      block_data.clear();
      parameters.clear();
      #ifdef DEBUG_VBC
      if ((block_data.size() != 0) || (parameters.size() != 0)) {
         std::cerr << "VBC CLEAR FAILED" << std::endl;
      }
      #endif
   }

   inline ARCH_HOSTDEV void VelocityBlockContainer::move(const vmesh::LocalID source, const vmesh::LocalID target) {
#ifdef USE_GPU
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      const vmesh::LocalID numberOfBlocks = cachedSize;
      #else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      #endif
#else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
#endif

#ifdef DEBUG_VBC
      bool ok = true;
      const vmesh::LocalID currentCapacity = block_data.capacity() / WID3;
      const vmesh::LocalID currentCapacityP = parameters.capacity() / BlockParams::N_VELOCITY_BLOCK_PARAMS;
      const vmesh::LocalID numberOfBlocksP = parameters.size() / BlockParams::N_VELOCITY_BLOCK_PARAMS;
      if (source >= numberOfBlocks)
         ok = false;
      if (source >= currentCapacity)
         ok = false;
      if (source >= numberOfBlocksP)
         ok = false;
      if (source >= currentCapacityP)
         ok = false;
      if (target >= numberOfBlocks)
         ok = false;
      if (target >= currentCapacity)
         ok = false;
      if (numberOfBlocks > currentCapacity)
         ok = false;
      if (source != numberOfBlocks - 1)
         ok = false; // only allows moving from last entry
      if (source != numberOfBlocksP - 1)
         ok = false;
         #ifdef USE_GPU
      if (cachedCapacity != currentCapacity)
         ok = false;
         #endif
      if (currentCapacityP != currentCapacity)
         ok = false;
      if (numberOfBlocksP != numberOfBlocks)
         ok = false;
      if (ok == false) {
            #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         std::stringstream ss;
            ss << "VBC ERROR: invalid source LID=" << source << " in copy, target=" << target << " #blocks=" << numberOfBlocks << " capacity=" << currentCapacity << std::endl;
         ss << "or sizes are wrong, data->size()=" << block_data.size() << " parameters.size()=" << parameters.size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
            #else
            printf("VBC error: invalid source LID=%u in copy, target=%u #blocks=%u capacity=%u \n or sizes are wrong, data->size()=%u parameters.size()=%u \n",
                source,
                target,
                numberOfBlocks,
                currentCapacity,
                (vmesh::LocalID)block_data.size(),
                (vmesh::LocalID)parameters.size());
         assert(0);
            #endif
      }
      #endif

      for (unsigned int i = 0; i < WID3; ++i) {
         block_data[target * WID3 + i] = block_data[source * WID3 + i];
      }
      for (int i = 0; i < BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) {
         parameters[target * BlockParams::N_VELOCITY_BLOCK_PARAMS + i] = parameters[source * BlockParams::N_VELOCITY_BLOCK_PARAMS + i];
      }
      // and remove last entry
      block_data.erase(block_data.begin() + WID3 * (numberOfBlocks - 1), block_data.begin() + WID3 * (numberOfBlocks));
      parameters.erase(parameters.begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS * (numberOfBlocks - 1), parameters.begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS * (numberOfBlocks));
#ifdef USE_GPU
      cachedSize--; // Note: if called from inside GPU kernel, cached size must be updated separately
#endif
   }

   inline void VelocityBlockContainer::exitInvalidLocalID(const vmesh::LocalID localID, const std::string& funcName) const {
      #ifdef DEBUG_VBC
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::stringstream ss;
      ss << "Process " << rank << ' ';
      ss << "Invalid localID " << localID << " used in function '" << funcName << "' max allowed value is " << numberOfBlocks << std::endl;
      std::cerr << ss.str();
      sleep(1);
      exit(1);
      #endif
   }
   inline ARCH_DEV void VelocityBlockContainer::exitInvalidLocalID(const vmesh::LocalID localID) const {
      #ifdef DEBUG_VBC
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      printf("Invalid localID %u used in VBC; max allowed value is %u\n", localID, numberOfBlocks);
      assert(0);
      #endif
   }

   inline ARCH_HOSTDEV double VelocityBlockContainer::getBlockAllocationFactor() { return BLOCK_ALLOCATION_FACTOR; }

   inline ARCH_HOSTDEV Realf* VelocityBlockContainer::getData() { return block_data.data(); }

   inline ARCH_HOSTDEV const Realf* VelocityBlockContainer::getData() const { return block_data.data(); }

   inline ARCH_HOSTDEV Realf* VelocityBlockContainer::getData(const vmesh::LocalID blockLID) {
      #ifdef DEBUG_VBC
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
         #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      if (blockLID >= numberOfBlocks) {
         exitInvalidLocalID(blockLID);
      }
         #else
      if (blockLID >= numberOfBlocks) {
         exitInvalidLocalID(blockLID, "getData");
      }
         #endif
#endif
      return block_data.data() + blockLID * WID3;
   }

   inline ARCH_HOSTDEV const Realf* VelocityBlockContainer::getData(const vmesh::LocalID blockLID) const {
      #ifdef DEBUG_VBC
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
         #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      if (blockLID >= numberOfBlocks) {
         exitInvalidLocalID(blockLID);
      }
         #else
      if (blockLID >= numberOfBlocks) {
         exitInvalidLocalID(blockLID, "const getData const");
      }
         #endif
#endif
      return block_data.data() + blockLID * WID3;
   }

   inline ARCH_HOSTDEV Real* VelocityBlockContainer::getParameters() { return parameters.data(); }

   inline ARCH_HOSTDEV const Real* VelocityBlockContainer::getParameters() const { return parameters.data(); }

   inline ARCH_HOSTDEV Real* VelocityBlockContainer::getParameters(const vmesh::LocalID blockLID) {
      #ifdef DEBUG_VBC
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
         #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      if (blockLID >= numberOfBlocks) {
         exitInvalidLocalID(blockLID);
      }
         #else
      if (blockLID >= numberOfBlocks) {
         exitInvalidLocalID(blockLID, "getParameters");
      }
      if (blockLID >= parameters.size() / BlockParams::N_VELOCITY_BLOCK_PARAMS) {
         exitInvalidLocalID(blockLID, "getParameters 2");
      }
         #endif
#endif
      return parameters.data() + blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   inline ARCH_HOSTDEV const Real* VelocityBlockContainer::getParameters(const vmesh::LocalID blockLID) const {
      #ifdef DEBUG_VBC
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
         #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      if (blockLID >= numberOfBlocks) {
         exitInvalidLocalID(blockLID);
      }
         #else
      if (blockLID >= numberOfBlocks) {
         exitInvalidLocalID(blockLID, "const getParameters const");
      }
      if (blockLID >= parameters.size() / BlockParams::N_VELOCITY_BLOCK_PARAMS) {
         exitInvalidLocalID(blockLID, "const getParameters const 2");
      }
         #endif
#endif
      return parameters.data() + blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   inline ARCH_HOSTDEV void VelocityBlockContainer::pop() {
#ifdef USE_GPU
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      const vmesh::LocalID numberOfBlocks = cachedSize;
      #else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      #endif
#else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
#endif

      if (numberOfBlocks == 0) {
         return;
      }
      block_data.erase(block_data.begin() + WID3 * (numberOfBlocks - 1), block_data.begin() + WID3 * (numberOfBlocks));
      parameters.erase(parameters.begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS * (numberOfBlocks - 1), parameters.begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS * (numberOfBlocks));
#ifdef USE_GPU
      cachedSize--; // Note: if called from inside kernel, cached size must be updated separately
#endif
   }

   /** Grows the size of the VBC, does not touch data */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back() {
#ifdef USE_GPU
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      const vmesh::LocalID numberOfBlocks = cachedSize;
      #else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      #endif
#else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
#endif
      vmesh::LocalID newIndex = numberOfBlocks;

      #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      const vmesh::LocalID currentCapacityD = block_data.capacity() / WID3;
      if (newIndex >= currentCapacityD) {
         assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back).");
      }
      block_data.device_resize((numberOfBlocks + 1) * WID3);
      parameters.device_resize((numberOfBlocks + 1) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #elif defined(USE_GPU)
      setNewCapacity(numberOfBlocks + 1, stream);
      block_data.resize((numberOfBlocks + 1) * WID3, true, stream);
      parameters.resize((numberOfBlocks + 1) * BlockParams::N_VELOCITY_BLOCK_PARAMS, true, stream);
      #else
      setNewCapacity(numberOfBlocks + 1);
      block_data.resize((numberOfBlocks + 1) * WID3, true);
      parameters.resize((numberOfBlocks + 1) * BlockParams::N_VELOCITY_BLOCK_PARAMS, true);
      #endif

#ifdef DEBUG_VBC
      const vmesh::LocalID currentCapacity = block_data.capacity() / WID3;
      const vmesh::LocalID currentCapacityP = parameters.capacity() / BlockParams::N_VELOCITY_BLOCK_PARAMS;
      if (newIndex >= currentCapacity || newIndex >= currentCapacityP) {
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         std::stringstream ss;
         ss << "VBC ERROR in push_back, LID=" << newIndex << " for new block is out of bounds" << std::endl;
         ss << "\t data->size()=" << block_data.size() << " parameters.size()=" << parameters.size() << std::endl;
         ss << "\t data->capacity()=" << block_data.capacity() << " parameters.capacity()=" << parameters.capacity() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
         #else
         printf("VBC ERROR in device push_back, LID=%u for new block is out of bounds\n  data->size()=%u parameters.size()=%u\n",
                newIndex,
                (vmesh::LocalID)block_data.size(),
                (vmesh::LocalID)parameters.size());
         assert(0);
         #endif
      }
      #endif

#ifdef USE_GPU
      cachedSize++; // Note: if called from inside GPU kernel, cached size must be updated separately
#endif
      return newIndex;
   }

   /** Grows the size of the VBC, sets data to zero */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back_and_zero() {
#ifdef USE_GPU
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      const vmesh::LocalID numberOfBlocks = cachedSize;
      #else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      #endif
#else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
#endif
      const vmesh::LocalID newIndex = numberOfBlocks;

      #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      const vmesh::LocalID currentCapacityD = block_data.capacity() / WID3;
      if (newIndex >= currentCapacityD) {
         assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back_and_zero).");
      }
      block_data.device_resize((numberOfBlocks + 1) * WID3, false);                                 // construct=false don't construct or set to zero (performed below)
      parameters.device_resize((numberOfBlocks + 1) * BlockParams::N_VELOCITY_BLOCK_PARAMS, false); // construct=false don't construct or set to zero (performed below)
      #elif defined(USE_GPU)
      setNewCapacity(numberOfBlocks + 1, stream);
      block_data.resize((numberOfBlocks + 1) * WID3, true, stream);
      parameters.resize((numberOfBlocks + 1) * BlockParams::N_VELOCITY_BLOCK_PARAMS, true, stream);
      #else
      setNewCapacity(numberOfBlocks + 1);
      block_data.resize((numberOfBlocks + 1) * WID3);
      parameters.resize((numberOfBlocks + 1) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

#ifdef DEBUG_VBC
      const vmesh::LocalID currentCapacity = block_data.capacity() / WID3;
      const vmesh::LocalID currentCapacityP = parameters.capacity() / BlockParams::N_VELOCITY_BLOCK_PARAMS;
      if (newIndex >= currentCapacity || newIndex >= currentCapacityP) {
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         std::stringstream ss;
         ss << "VBC ERROR in push_back_and_zero, LID=" << newIndex << " for new block is out of bounds" << std::endl;
         ss << "\t data->size()=" << block_data.size() << " parameters.size()=" << parameters.size() << std::endl;
         ss << "\t data->capacity()=" << block_data.capacity() << " parameters.capacity()=" << parameters.capacity() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
         #else
         printf("VBC ERROR in device push_back_and_zero, LID=%u for new block is out of bounds \n data->size()=%u parameters.size()=%u \n",
                newIndex,
                (vmesh::LocalID)block_data.size(),
                (vmesh::LocalID)parameters.size());
         assert(0);
         #endif
      }
      #endif

      for (size_t i = 0; i < WID3; ++i) {
         block_data[newIndex * WID3 + i] = 0.0;
      }
      for (size_t i = 0; i < BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) {
         parameters[newIndex * BlockParams::N_VELOCITY_BLOCK_PARAMS + i] = 0.0;
      }
#ifdef USE_GPU
      cachedSize++; // Note: if called from inside GPU kernel, cached size must be updated separately
#endif
      return newIndex;
   }

   /** Grows the size of the VBC, does not touch data */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back(const uint32_t N_blocks) {
#ifdef USE_GPU
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      const vmesh::LocalID numberOfBlocks = cachedSize;
      const vmesh::LocalID currentCapacity = cachedCapacity;
      #else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      const vmesh::LocalID currentCapacity = block_data.capacity() / WID3;
      #endif
#else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      const vmesh::LocalID currentCapacity = block_data.capacity() / WID3;
#endif
      const vmesh::LocalID newIndex = numberOfBlocks;

      #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      if (newIndex + N_blocks >= currentCapacity - 1) {
         assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back N_blocks).");
      }
      block_data.device_resize((numberOfBlocks + N_blocks) * WID3);
      parameters.device_resize((numberOfBlocks + N_blocks) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #elif defined(USE_GPU)
      setNewCapacity(numberOfBlocks + N_blocks, stream);
      block_data.resize((numberOfBlocks + N_blocks) * WID3, true, stream);
      parameters.resize((numberOfBlocks + N_blocks) * BlockParams::N_VELOCITY_BLOCK_PARAMS, true, stream);
      #else
      setNewCapacity(numberOfBlocks + N_blocks);
      block_data.resize((numberOfBlocks + N_blocks) * WID3);
      parameters.resize((numberOfBlocks + N_blocks) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

#ifdef USE_GPU
      cachedSize += N_blocks; // Note: if called from inside GPU kernel, cached size must be updated separately
#endif
      return newIndex;
   }

   /** Grows the size of the VBC, sets data to zero */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back_and_zero(const uint32_t N_blocks) {
#ifdef USE_GPU
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      const vmesh::LocalID numberOfBlocks = cachedSize;
      const vmesh::LocalID currentCapacity = cachedCapacity;
      #else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      const vmesh::LocalID currentCapacity = block_data.capacity() / WID3;
      #endif
#else
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      const vmesh::LocalID currentCapacity = block_data.capacity() / WID3;
#endif
      const vmesh::LocalID newIndex = numberOfBlocks;

      #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      if (newIndex + N_blocks >= currentCapacity - 1) {
         assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back_and_zero N_blocks).");
      }
      block_data.device_resize((numberOfBlocks + N_blocks) * WID3, false);                                 // construct=false don't construct or set to zero (performed below)
      parameters.device_resize((numberOfBlocks + N_blocks) * BlockParams::N_VELOCITY_BLOCK_PARAMS, false); // construct=false don't construct or set to zero (performed below)
      #elif defined(USE_GPU)
      setNewCapacity(numberOfBlocks + N_blocks, stream);
      block_data.resize((numberOfBlocks + N_blocks) * WID3, true, stream);
      parameters.resize((numberOfBlocks + N_blocks) * BlockParams::N_VELOCITY_BLOCK_PARAMS, true, stream);
      #else
      setNewCapacity(numberOfBlocks + N_blocks);
      block_data.resize((numberOfBlocks + N_blocks) * WID3);
      parameters.resize((numberOfBlocks + N_blocks) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

#if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Clear velocity block data to zero values
      Realf* zero_blocks = block_data.data();
      Real* zero_parameters = parameters.data();
      // block_data.optimizeGPU(stream);
      // parameters.optimizeGPU(stream);
      CHK_ERR(gpuMemsetAsync(zero_blocks + newIndex * WID3, 0, WID3 * N_blocks * sizeof(Realf), stream));
      CHK_ERR(gpuMemsetAsync(zero_parameters + newIndex * BlockParams::N_VELOCITY_BLOCK_PARAMS, 0, BlockParams::N_VELOCITY_BLOCK_PARAMS * N_blocks * sizeof(Real), stream));
      CHK_ERR(gpuStreamSynchronize(stream));
      #else
      // Clear velocity block data to zero values
      for (size_t i = 0; i < WID3 * N_blocks; ++i) {
         block_data[newIndex * WID3 + i] = 0.0;
      }
      for (size_t i = 0; i < BlockParams::N_VELOCITY_BLOCK_PARAMS * N_blocks; ++i) {
         parameters[newIndex * BlockParams::N_VELOCITY_BLOCK_PARAMS + i] = 0.0;
      }
      #endif

#ifdef USE_GPU
      cachedSize += N_blocks; // Note: if called from inside GPU kernel, cached size must be updated separately
#endif
      return newIndex;
   }

#ifdef USE_GPU
   inline bool VelocityBlockContainer::setNewCapacity(const vmesh::LocalID reqCapacity, gpuStream_t stream = 0) {
      if (stream == 0) {
         gpuStream_t stream = gpu_getStream();
      }
      const vmesh::LocalID numberOfBlocks = cachedSize;
      const vmesh::LocalID currentCapacity = cachedCapacity;
#else
   inline bool VelocityBlockContainer::setNewCapacity(const vmesh::LocalID reqCapacity) {
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      const vmesh::LocalID currentCapacity = block_data.capacity() / WID3;
#endif
      // Note: No longer ever recapacitates down in size.

      // Reallocate so that free space is current * block_allocation_factor blocks,
      // and at least two in case of having zero blocks.
      vmesh::LocalID newCapacity = (reqCapacity > 2) ? reqCapacity : 2;
      if (currentCapacity >= newCapacity) {
         return false; // Still have enough buffer
      }
      if (newCapacity < numberOfBlocks) {
         std::cerr << " ERROR! Trying to recapacitate to " << newCapacity << " when VBC already contains " << numberOfBlocks << " blocks!" << std::endl;
         return false;
      }
      newCapacity = (vmesh::LocalID)((double)newCapacity * BLOCK_ALLOCATION_FACTOR);
      #ifdef USE_GPU
      // Passing eco flag = true to reserve tells splitvector we manage padding manually.
      block_data.reserve(newCapacity * WID3, true, stream);
      parameters.reserve(newCapacity * BlockParams::N_VELOCITY_BLOCK_PARAMS, true, stream);
      #else
      block_data.reserve(newCapacity * WID3);
      parameters.reserve(newCapacity * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

#ifdef USE_GPU
      cachedCapacity = newCapacity;
#endif
      return true;
   }

   inline bool VelocityBlockContainer::setNewCapacityShrink(const vmesh::LocalID reqCapacity) {
      #ifdef USE_HIP
      return setNewCapacity(reqCapacity);
      #endif
      // Reallocate so that capacity matches requested value.
      vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      if (reqCapacity < numberOfBlocks) {
         std::cerr << " ERROR! Trying to recapacitate to " << reqCapacity << " when VBC already contains " << numberOfBlocks << " blocks!" << std::endl;
         // Could enforce a minimum value here, but better to catch errors in call logic.
         return false;
      }
      vmesh::LocalID newCapacity = std::max(reqCapacity, (vmesh::LocalID)2); // At least 2 blocks
#ifdef USE_GPU
      gpuStream_t stream = gpu_getStream();
      // Overwrite/swap causing data corruption on LUMI-G, use reallocate method instead.
      block_data.reallocate(newCapacity * WID3, stream);
      parameters.reallocate(newCapacity * BlockParams::N_VELOCITY_BLOCK_PARAMS, stream);
      CHK_ERR(gpuStreamSynchronize(stream));
      cachedCapacity = newCapacity;
#else
      // Create with larger size (capacity), then resize down to actual size
      std::vector<Realf, aligned_allocator<Realf, WID3>> block_data_new(newCapacity * WID3);
      std::vector<Real, aligned_allocator<Real, BlockParams::N_VELOCITY_BLOCK_PARAMS>> parameters_new(newCapacity * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      block_data_new.resize(numberOfBlocks * WID3);
      parameters_new.resize(numberOfBlocks * BlockParams::N_VELOCITY_BLOCK_PARAMS);
      for (size_t i = 0; i < numberOfBlocks * WID3; ++i) {
         block_data_new[i] = block_data[i];
      }
      for (size_t i = 0; i < numberOfBlocks * BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) {
         parameters_new[i] = parameters[i];
      }
      block_data_new.swap(block_data);
      parameters_new.swap(parameters);
#endif
      return true;
   }

   inline ARCH_HOSTDEV void VelocityBlockContainer::resize(const vmesh::LocalID newSize) { setNewSize(newSize); }

   inline ARCH_HOSTDEV bool VelocityBlockContainer::setNewSize(const vmesh::LocalID newSize) {
      // Does not set new added blocks to zero
#ifdef USE_GPU
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      setNewCapacity(newSize, stream);
      parameters.resize((newSize)*BlockParams::N_VELOCITY_BLOCK_PARAMS, true, stream);
      block_data.resize((newSize)*WID3, true, stream);
      #else
      const vmesh::LocalID currentCapacity = block_data.capacity() / WID3;
      assert(newSize <= currentCapacity && "ERROR! Attempting to grow block container on-device beyond capacity (::setNewSize).");
      block_data.device_resize((newSize)*WID3, false);                                 // construct=false don't construct or set to zero
      parameters.device_resize((newSize)*BlockParams::N_VELOCITY_BLOCK_PARAMS, false); // construct=false don't construct or set to zero
      #endif
#else
      setNewCapacity(newSize);
      block_data.resize((newSize)*WID3);
      parameters.resize((newSize)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
#endif

#ifdef USE_GPU
      cachedSize = newSize; // Note: if called from inside GPU kernel, cached size must be updated separately
#endif
      return true;
   }

   /** Return the number of existing velocity blocks.
    * @return Number of existing velocity blocks.*/
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::size() const {
#ifdef USE_GPU
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      return block_data.size() / WID3;
      #else
#ifdef DEBUG_VBC
      const size_t currentSize = block_data.size() / WID3;
      if (currentSize != cachedSize) {
         printf("VBC CHECK ERROR: cached size mismatch, %lu vs %lu in %s : %d\n", currentSize, cachedSize, __FILE__, __LINE__);
      }
      #endif
      return cachedSize;
      #endif
#else
      return block_data.size() / WID3;
#endif
   }

   inline ARCH_HOSTDEV size_t VelocityBlockContainer::sizeInBytes() const {
#ifdef USE_GPU
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      return block_data.size() * sizeof(Realf) + parameters.size() * sizeof(Real);
      #else
#ifdef DEBUG_VBC
      const size_t currentSize = block_data.size() / WID3;
      if (currentSize != cachedSize) {
         printf("VBC CHECK ERROR: cached size mismatch, %lu vs %lu in %s : %d\n", currentSize, cachedSize, __FILE__, __LINE__);
      }
      #endif
      return cachedSize * WID3 * sizeof(Realf) + cachedSize * BlockParams::N_VELOCITY_BLOCK_PARAMS * sizeof(Real);
      #endif
#else
      return block_data.size() * sizeof(Realf) + parameters.size() * sizeof(Real);
#endif
   }

   // inline void VelocityBlockContainer::swap(VelocityBlockContainer& vbc) {
   //    block_data.swap(vbc.block_data);
   //    parameters.swap(vbc.parameters);
   // }

#ifdef USE_GPU
   /** GPU-ONLY FUNCTIONS */
   inline void VelocityBlockContainer::setNewCachedSize(const vmesh::LocalID newSize) {
      // Should only be used to update host-side size if resizing on device
      cachedSize = newSize;
   }
   inline void VelocityBlockContainer::updateCachedSize() {
      // More secure page-faulting way to update cached size
      cachedSize = block_data.size() / WID3;
   }
   inline void VelocityBlockContainer::updateCachedCapacity() {
      // Should not be needed, added as an optional safeguard
      cachedCapacity = block_data.capacity() / WID3;
   }
   inline void VelocityBlockContainer::print_addresses() { printf("GPU block_data %p\n GPU parameters %p\n", &block_data, &parameters); }
   inline void VelocityBlockContainer::gpu_prefetchHost(gpuStream_t stream = 0) {
      if (stream == 0) {
         gpuStream_t stream = gpu_getStream();
      }
      block_data.optimizeCPU(stream);
      parameters.optimizeCPU(stream);
      return;
   }
   inline void VelocityBlockContainer::gpu_prefetchDevice(gpuStream_t stream = 0) {
      if (stream == 0) {
         gpuStream_t stream = gpu_getStream();
      }
      block_data.optimizeGPU(stream);
      parameters.optimizeGPU(stream);
      return;
   }
#endif

#ifdef DEBUG_VBC
   inline const Realf& VelocityBlockContainer::getData(const vmesh::LocalID blockLID, const unsigned int cell) const {
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      bool ok = true;
      if (cell >= WID3) {
         ok = false;
      }
      if (blockLID >= numberOfBlocks) {
         ok = false;
      }
      if (blockLID * WID3 + cell >= block_data.size()) {
         ok = false;
      }
      if (ok == false) {
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in getData, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " data->size()=" << block_data.size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }
      return block_data[blockLID * WID3 + cell];
   }

   inline const Real& VelocityBlockContainer::getParameters(const vmesh::LocalID blockLID, const unsigned int cell) const {
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      bool ok = true;
      if (cell >= BlockParams::N_VELOCITY_BLOCK_PARAMS) {
         ok = false;
      }
      if (blockLID >= numberOfBlocks) {
         ok = false;
      }
      if (blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS + cell >= parameters.size()) {
         ok = false;
      }
      if (ok == false) {
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in getParameters, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " parameters.size()=" << parameters.size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }
      return parameters[blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS + cell];
   }

   inline void VelocityBlockContainer::setData(const vmesh::LocalID blockLID, const unsigned int cell, const Realf value) {
      const vmesh::LocalID numberOfBlocks = block_data.size() / WID3;
      bool ok = true;
      if (cell >= WID3) {
         ok = false;
      }
      if (blockLID >= numberOfBlocks) {
         ok = false;
      }
      if (blockLID * WID3 + cell >= block_data.size()) {
         ok = false;
      }
      if (ok == false) {
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in setData, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " data->size()=" << block_data.size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }

      block_data[blockLID * WID3 + cell] = value;
   }
#endif //debug VBC

} // namespace vmesh

#endif
