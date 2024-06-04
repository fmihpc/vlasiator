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

#include <vector>
#include <stdio.h>

#include "common.h"
#include "unistd.h"

#ifdef DEBUG_VLASIATOR
   #ifndef DEBUG_VBC
   #define DEBUG_VBC
   #endif
#endif

#ifdef DEBUG_VBC
#include <sstream>
#endif

#include "arch/arch_device_api.h"
#ifdef USE_GPU
   #include "arch/gpu_base.hpp"
   // Place block data and parameters inside splitvectors utilizing unified memory
   #include "include/splitvector/splitvec.h"
#else
   // GPU allocation factors are stored in arch/gpu_base.hpp
   static const double BLOCK_ALLOCATION_FACTOR = 1.1;
   static const double BLOCK_ALLOCATION_PADDING = 1.3;
#endif

using namespace std;

namespace vmesh {

   class VelocityBlockContainer {
   public:

      VelocityBlockContainer();
      ~VelocityBlockContainer();
      VelocityBlockContainer(const VelocityBlockContainer& other);
      const VelocityBlockContainer& operator=(const VelocityBlockContainer& other);
      void gpu_destructor();

      ARCH_HOSTDEV vmesh::LocalID capacity() const;
      ARCH_HOSTDEV size_t capacityInBytes() const;
      void clear(bool shrink=true);
      ARCH_HOSTDEV void move(const vmesh::LocalID& source,const vmesh::LocalID& target);
      ARCH_HOSTDEV static double getBlockAllocationFactor();
      ARCH_HOSTDEV Realf* getData();
      ARCH_HOSTDEV const Realf* getData() const;
      ARCH_HOSTDEV Realf* getData(const vmesh::LocalID& blockLID);
      ARCH_HOSTDEV const Realf* getData(const vmesh::LocalID& blockLID) const;
      ARCH_HOSTDEV Real* getParameters();
      ARCH_HOSTDEV const Real* getParameters() const;
      ARCH_HOSTDEV Real* getParameters(const vmesh::LocalID& blockLID);
      ARCH_HOSTDEV const Real* getParameters(const vmesh::LocalID& blockLID) const;
      ARCH_HOSTDEV void pop();
      ARCH_HOSTDEV vmesh::LocalID push_back();
      ARCH_HOSTDEV vmesh::LocalID push_back_and_zero();
      ARCH_HOSTDEV vmesh::LocalID push_back(const uint32_t& N_blocks);
      ARCH_HOSTDEV vmesh::LocalID push_back_and_zero(const uint32_t& N_blocks);
      ARCH_HOSTDEV void resize(const vmesh::LocalID& newSize);
      ARCH_HOSTDEV bool setNewSize(const vmesh::LocalID& newSize);
      ARCH_HOSTDEV vmesh::LocalID size() const;
      ARCH_HOSTDEV size_t sizeInBytes() const;
      // ARCH_HOSTDEV void swap(VelocityBlockContainer& vbc);

#ifdef USE_GPU // for GPU version
      bool setNewCapacity(const vmesh::LocalID& capacity, gpuStream_t stream);
      void gpu_prefetchHost(gpuStream_t stream);
      void gpu_prefetchDevice(gpuStream_t stream);
      void print_addresses();
#else
      bool setNewCapacity(const vmesh::LocalID& capacity);
#endif

      #ifdef DEBUG_VBC
      const Realf& getData(const vmesh::LocalID& blockLID,const unsigned int& cell) const;
      const Real& getParameters(const vmesh::LocalID& blockLID,const unsigned int& i) const;
      void setData(const vmesh::LocalID& blockLID,const unsigned int& cell,const Realf& value);
      #endif
    private:
      void exitInvalidLocalID(const vmesh::LocalID& localID,const std::string& funcName) const;
      ARCH_DEV void exitInvalidLocalID(const vmesh::LocalID& localID) const;

#ifdef USE_GPU
      split::SplitVector<Realf> *block_data;
      split::SplitVector<Real> *parameters;
#else
      std::vector<Realf,aligned_allocator<Realf,WID3> > *block_data;
      std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > *parameters;
#endif

   };

   inline VelocityBlockContainer::VelocityBlockContainer() {
#ifdef USE_GPU
      block_data= new split::SplitVector<Realf>(WID3);
      parameters= new split::SplitVector<Real>(BlockParams::N_VELOCITY_BLOCK_PARAMS);
#else
      block_data = new std::vector<Realf,aligned_allocator<Realf,WID3>>(WID3);
      parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS>>(BlockParams::N_VELOCITY_BLOCK_PARAMS);
#endif
      block_data->clear();
      parameters->clear();
      // gpuStream_t stream = gpu_getStream();
   }

   inline VelocityBlockContainer::~VelocityBlockContainer() {
      gpu_destructor();
   }
   inline void VelocityBlockContainer::gpu_destructor() {
      if (block_data) {
         delete block_data;
         block_data = 0;
      }
      if (parameters) {
         delete parameters;
         parameters = 0;
      }
   }

   inline VelocityBlockContainer::VelocityBlockContainer(const VelocityBlockContainer& other) {
#ifdef USE_GPU
      block_data = new split::SplitVector<Realf>(other.block_data->capacity());
      parameters = new split::SplitVector<Real>(other.parameters->capacity());
      // Overwrite is like a copy assign but takes a stream
      gpuStream_t stream = gpu_getStream();
      block_data->overwrite(*(other.block_data),stream);
      parameters->overwrite(*(other.parameters),stream);
#else
      block_data = new std::vector<Realf,aligned_allocator<Realf,WID3>>(*(other.block_data));
      parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS>>(*(other.parameters));
      block_data->reserve(other.block_data->capacity());
      parameters->reserve(other.parameters->capacity());
#endif
   }

   inline const VelocityBlockContainer& VelocityBlockContainer::operator=(const VelocityBlockContainer& other) {
      #ifdef USE_GPU
      gpuStream_t stream = gpu_getStream();
      block_data->reserve(other.block_data->capacity(), true, stream);
      parameters->reserve(other.parameters->capacity(), true, stream);
      // Overwrite is like a copy assign but takes a stream
      block_data->overwrite(*(other.block_data),stream);
      parameters->overwrite(*(other.parameters),stream);
      #else
      *block_data = *(other.block_data);
      *parameters = *(other.parameters);
      block_data->reserve(other.block_data->capacity());
      parameters->reserve(other.parameters->capacity());
      #endif
      return *this;
   }

   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::capacity() const {
      const vmesh::LocalID currentCapacity = block_data->capacity() / WID3;
      return currentCapacity;
   }

   inline ARCH_HOSTDEV size_t VelocityBlockContainer::capacityInBytes() const {
      const vmesh::LocalID currentCapacity = block_data->capacity();
      const vmesh::LocalID parametersCapacity = parameters->capacity();
      return currentCapacity*sizeof(Realf) + parametersCapacity*sizeof(Real);
   }

   /** Clears VelocityBlockContainer data and deallocates all memory
    * reserved for velocity blocks.*/
   inline void VelocityBlockContainer::clear(bool shrink) {
      // GPU DEBUG: For some reason, calling just clear seems broken?
      size_t capacity = block_data->capacity()/WID3;
      if (shrink) {
         capacity = 1;
      }
      delete block_data;
      delete parameters;
#ifdef USE_GPU
      block_data = new split::SplitVector<Realf>(capacity*WID3);
      parameters = new split::SplitVector<Real>(capacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
#else
      block_data = new std::vector<Realf,aligned_allocator<Realf,WID3>>(capacity*WID3);
      parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS>>(capacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
#endif
      block_data->clear();
      parameters->clear();
      if ((block_data->size() != 0) || (parameters->size() != 0)) {
         std::cerr<<"VBC CLEAR FAILED"<<std::endl;
      }
   }

   inline ARCH_HOSTDEV void VelocityBlockContainer::move(const vmesh::LocalID& source,const vmesh::LocalID& target) {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      #endif
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;

      #ifdef DEBUG_VBC
         bool ok = true;
         const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
         const vmesh::LocalID currentCapacityP = parameters->capacity()/BlockParams::N_VELOCITY_BLOCK_PARAMS;
         const vmesh::LocalID numberOfBlocksP = parameters->size()/BlockParams::N_VELOCITY_BLOCK_PARAMS;
         if (source >= numberOfBlocks) ok = false;
         if (source >= currentCapacity) ok = false;
         if (source >= numberOfBlocksP) ok = false;
         if (source >= currentCapacityP) ok = false;
         if (target >= numberOfBlocks) ok = false;
         if (target >= currentCapacity) ok = false;
         if (numberOfBlocks > currentCapacity) ok = false;
         if (source != numberOfBlocks-1) ok = false; // only allows moving from last entry
         if (source != numberOfBlocksP-1) ok = false;
         if (currentCapacityP != currentCapacity) ok = false;
         if (numberOfBlocksP != numberOfBlocks) ok = false;
         if (ok == false) {
            #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
            std::stringstream ss;
            ss << "VBC ERROR: invalid source LID=" << source << " in copy, target=" << target << " #blocks=" << numberOfBlocks << " capacity=" << currentCapacity << std::endl;
            ss << "or sizes are wrong, data->size()=" << block_data->size() << " parameters->size()=" << parameters->size() << std::endl;
            std::cerr << ss.str();
            sleep(1);
            exit(1);
            #else
            printf("VBC error: invalid source LID=%u in copy, target=%u #blocks=%u capacity=%u \n or sizes are wrong, data->size()=%u parameters->size()=%u \n",
                   source,target,numberOfBlocks,currentCapacity, (vmesh::LocalID)block_data->size(),(vmesh::LocalID)parameters->size());
            assert(0);
            #endif
         }
      #endif

      for (unsigned int i=0; i<WID3; ++i) {
         (*block_data)[target*WID3+i] = (*block_data)[source*WID3+i];
      }
      for (int i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) {
         (*parameters)[target*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = (*parameters)[source*BlockParams::N_VELOCITY_BLOCK_PARAMS+i];
      }
      // and remove last entry
      block_data->erase(block_data->begin() + WID3*(numberOfBlocks-1),
                        block_data->begin() + WID3*(numberOfBlocks));
      parameters->erase(parameters->begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS*(numberOfBlocks-1),
                        parameters->begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS*(numberOfBlocks));

   }

   inline void VelocityBlockContainer::exitInvalidLocalID(const vmesh::LocalID& localID,const std::string& funcName) const {
      #ifdef DEBUG_VBC
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      std::stringstream ss;
      ss << "Process " << rank << ' ';
      ss << "Invalid localID " << localID << " used in function '" << funcName << "' max allowed value is " << numberOfBlocks << std::endl;
      std::cerr << ss.str();
      sleep(1);
      exit(1);
      #endif
   }
   inline ARCH_DEV void VelocityBlockContainer::exitInvalidLocalID(const vmesh::LocalID& localID) const {
      #ifdef DEBUG_VBC
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      printf("Invalid localID %u used in VBC; max allowed value is %u\n",localID,numberOfBlocks);
      assert(0);
      #endif
   }

   inline ARCH_HOSTDEV double VelocityBlockContainer::getBlockAllocationFactor() {
      return BLOCK_ALLOCATION_FACTOR;
   }

   inline ARCH_HOSTDEV Realf* VelocityBlockContainer::getData() {
      return block_data->data();
   }

   inline ARCH_HOSTDEV const Realf* VelocityBlockContainer::getData() const {
      return block_data->data();
   }

   inline ARCH_HOSTDEV Realf* VelocityBlockContainer::getData(const vmesh::LocalID& blockLID) {
      #ifdef DEBUG_VBC
         const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
         #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID);
         #else
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"getData");
         #endif
      #endif
      return block_data->data() + blockLID*WID3;
   }

   inline ARCH_HOSTDEV const Realf* VelocityBlockContainer::getData(const vmesh::LocalID& blockLID) const {
      #ifdef DEBUG_VBC
         const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
         #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID);
         #else
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"const getData const");
         #endif
      #endif
      return block_data->data() + blockLID*WID3;
   }

#ifdef USE_GPU
   inline void VelocityBlockContainer::print_addresses() {
      printf("GPU block_data %p\n GPU parameters %p\n",block_data,parameters);
   }
   inline void VelocityBlockContainer::gpu_prefetchHost(gpuStream_t stream=0) {
      if (stream==0) {
         gpuStream_t stream = gpu_getStream();
      }
      block_data->optimizeCPU(stream);
      parameters->optimizeCPU(stream);
      return;
   }
   inline void VelocityBlockContainer::gpu_prefetchDevice(gpuStream_t stream=0) {
      if (stream==0) {
         gpuStream_t stream = gpu_getStream();
      }
      block_data->optimizeGPU(stream);
      parameters->optimizeGPU(stream);
      return;
   }
#endif

   inline ARCH_HOSTDEV Real* VelocityBlockContainer::getParameters() {
      return parameters->data();
   }

   inline ARCH_HOSTDEV const Real* VelocityBlockContainer::getParameters() const {
      return parameters->data();
   }

   inline ARCH_HOSTDEV Real* VelocityBlockContainer::getParameters(const vmesh::LocalID& blockLID) {
      #ifdef DEBUG_VBC
         const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
         #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID);
         #else
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"getParameters");
         if (blockLID >= parameters->size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) exitInvalidLocalID(blockLID,"getParameters 2");
         #endif
      #endif
      return parameters->data() + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   inline ARCH_HOSTDEV const Real* VelocityBlockContainer::getParameters(const vmesh::LocalID& blockLID) const {
      #ifdef DEBUG_VBC
         const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
         #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID);
         #else
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"const getParameters const");
         if (blockLID >= parameters->size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) exitInvalidLocalID(blockLID,"const getParameters const 2");
         #endif
      #endif
      return parameters->data() + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   inline ARCH_HOSTDEV void VelocityBlockContainer::pop() {
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;

      if (numberOfBlocks == 0) {
         return;
      }
      block_data->erase(block_data->begin() + WID3*(numberOfBlocks-1),
                        block_data->begin() + WID3*(numberOfBlocks));
      parameters->erase(parameters->begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS*(numberOfBlocks-1),
                        parameters->begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS*(numberOfBlocks));
   }

   /** Grows the size of the VBC, does not touch data */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back() {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      #endif
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;

      vmesh::LocalID newIndex = numberOfBlocks;
      #ifdef DEBUG_VBC
      const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
      const vmesh::LocalID currentCapacityP = parameters->capacity()/BlockParams::N_VELOCITY_BLOCK_PARAMS;
      if (newIndex >= currentCapacity || newIndex >= currentCapacityP) {
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         std::stringstream ss;
         ss << "VBC ERROR in push_back, LID=" << newIndex << " for new block is out of bounds" << std::endl;
         ss << "\t data->size()=" << block_data->size()  << " parameters->size()=" << parameters->size() << std::endl;
         ss << "\t data->capacity()=" << block_data->capacity()  << " parameters->capacity()=" << parameters->capacity() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
         #else
         printf("VBC ERROR in device push_back, LID=%u for new block is out of bounds\n  data->size()=%u parameters->size()=%u\n",
                newIndex,(vmesh::LocalID)block_data->size(),(vmesh::LocalID)parameters->size());
         assert(0);
         #endif
      }
      #endif

      #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      const vmesh::LocalID currentCapacityD = block_data->capacity()/WID3;
      if (newIndex >= currentCapacityD) {
         assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back).");
      }
      block_data->device_resize((numberOfBlocks+1)*WID3);
      parameters->device_resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #elif defined(USE_GPU)
      setNewCapacity(numberOfBlocks+1,stream);
      block_data->resize((numberOfBlocks+1)*WID3,true,stream);
      parameters->resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true,stream);
      #else
      setNewCapacity(numberOfBlocks+1);
      block_data->resize((numberOfBlocks+1)*WID3,true);
      parameters->resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true);
      #endif

      return newIndex;
   }

   /** Grows the size of the VBC, sets data to zero */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back_and_zero() {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      #endif
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;

      const vmesh::LocalID newIndex = numberOfBlocks;
      #ifdef DEBUG_VBC
      const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
      const vmesh::LocalID currentCapacityP = parameters->capacity()/BlockParams::N_VELOCITY_BLOCK_PARAMS;
      if (newIndex >= currentCapacity || newIndex >= currentCapacityP) {
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         std::stringstream ss;
         ss << "VBC ERROR in push_back_and_zero, LID=" << newIndex << " for new block is out of bounds" << std::endl;
         ss << "\t data->size()=" << block_data->size()  << " parameters->size()=" << parameters->size() << std::endl;
         ss << "\t data->capacity()=" << block_data->capacity()  << " parameters->capacity()=" << parameters->capacity() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
         #else
         printf("VBC ERROR in device push_back_and_zero, LID=%u for new block is out of bounds \n data->size()=%u parameters->size()=%u \n",
                newIndex,(vmesh::LocalID)block_data->size(),(vmesh::LocalID)parameters->size());
         assert(0);
         #endif
      }
      #endif

      #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      const vmesh::LocalID currentCapacityD = block_data->capacity()/WID3;
      if (newIndex >= currentCapacityD) {
         assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back_and_zero).");
      }
      block_data->device_resize((numberOfBlocks+1)*WID3, false); //construct=false don't construct or set to zero (performed below)
      parameters->device_resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS, false); //construct=false don't construct or set to zero (performed below)
      #elif defined(USE_GPU)
      setNewCapacity(numberOfBlocks+1,stream);
      block_data->resize((numberOfBlocks+1)*WID3,true,stream);
      parameters->resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true,stream);
      #else
      setNewCapacity(numberOfBlocks+1);
      block_data->resize((numberOfBlocks+1)*WID3);
      parameters->resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

      for (size_t i=0; i<WID3; ++i) {
         (*block_data)[newIndex*WID3+i] = 0.0;
      }
      for (size_t i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) {
         (*parameters)[newIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = 0.0;
      }

      return newIndex;
   }

   /** Grows the size of the VBC, does not touch data */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back(const uint32_t& N_blocks) {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      #endif
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;

      const vmesh::LocalID newIndex = numberOfBlocks;

      #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      if (newIndex + N_blocks >= currentCapacity-1) {
         assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back N_blocks).");
      }
      block_data->device_resize((numberOfBlocks+N_blocks)*WID3);
      parameters->device_resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #elif defined(USE_GPU)
      setNewCapacity(numberOfBlocks+N_blocks,stream);
      block_data->resize((numberOfBlocks+N_blocks)*WID3,true,stream);
      parameters->resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true,stream);
      #else
      setNewCapacity(numberOfBlocks+N_blocks);
      block_data->resize((numberOfBlocks+N_blocks)*WID3);
      parameters->resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

      return newIndex;
   }

   /** Grows the size of the VBC, sets data to zero */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back_and_zero(const uint32_t& N_blocks) {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      gpuStream_t stream = gpu_getStream();
      #endif
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;

      const vmesh::LocalID newIndex = numberOfBlocks;

      #if defined(USE_GPU) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      if (newIndex + N_blocks >= currentCapacity-1) {
         assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back_and_zero N_blocks).");
      }
      block_data->device_resize((numberOfBlocks+N_blocks)*WID3, false); //construct=false don't construct or set to zero (performed below)
      parameters->device_resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS, false); //construct=false don't construct or set to zero (performed below)
      #elif defined(USE_GPU)
      setNewCapacity(numberOfBlocks+N_blocks,stream);
      block_data->resize((numberOfBlocks+N_blocks)*WID3,true,stream);
      parameters->resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true,stream);
      #else
      setNewCapacity(numberOfBlocks+N_blocks);
      block_data->resize((numberOfBlocks+N_blocks)*WID3);
      parameters->resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Clear velocity block data to zero values
      Realf* zero_blocks = block_data->data();
      Real* zero_parameters = parameters->data();
      // block_data->optimizeGPU(stream);
      // parameters->optimizeGPU(stream);
      CHK_ERR( gpuMemsetAsync(zero_blocks + newIndex*WID3, 0, WID3*N_blocks*sizeof(Realf), stream) );
      CHK_ERR( gpuMemsetAsync(zero_parameters + newIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS, 0, BlockParams::N_VELOCITY_BLOCK_PARAMS*N_blocks*sizeof(Real), stream) );
      CHK_ERR( gpuStreamSynchronize(stream) );
      #else
      // Clear velocity block data to zero values
      for (size_t i=0; i<WID3*N_blocks; ++i) {
         (*block_data)[newIndex*WID3+i] = 0.0;
      }
      for (size_t i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS*N_blocks; ++i) {
         (*parameters)[newIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = 0.0;
      }
      #endif

      return newIndex;
   }

#ifdef USE_GPU
   inline bool VelocityBlockContainer::setNewCapacity(const vmesh::LocalID& reqCapacity, gpuStream_t stream=0) {
      if (stream==0) {
         gpuStream_t stream = gpu_getStream();
      }
#else
   inline bool VelocityBlockContainer::setNewCapacity(const vmesh::LocalID& reqCapacity) {
#endif
      // Note: No longer ever recapacitates down in size.
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;

      // Reallocate so that free space is current * block_allocation_padding blocks,
      // and at least two in case of having zero blocks.
      vmesh::LocalID newCapacity = (reqCapacity > 2) ? reqCapacity : 2;
      if (currentCapacity > BLOCK_ALLOCATION_FACTOR * newCapacity) {
         return false; // Still have enough buffer
      }
      if (newCapacity < numberOfBlocks) {
         std::cerr<<" ERROR! Trying to recapacitate to "<<newCapacity<<" when VBC already contains "<<numberOfBlocks<<" blocks!"<<std::endl;
         return false;
      }
      newCapacity *= BLOCK_ALLOCATION_PADDING;
      #ifdef USE_GPU
      // Passing eco flag = true to reserve tells splitvector we manage padding manually.
      block_data->reserve(newCapacity*WID3, true, stream);
      parameters->reserve(newCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS, true, stream);
      #else
      block_data->reserve(newCapacity*WID3);
      parameters->reserve(newCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

      // #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // block_data->optimizeGPU(stream);
      // parameters->optimizeGPU(stream);
      // #endif
      return true;
   }

   inline ARCH_HOSTDEV void VelocityBlockContainer::resize(const vmesh::LocalID& newSize) {
      setNewSize(newSize);
   }

   inline ARCH_HOSTDEV bool VelocityBlockContainer::setNewSize(const vmesh::LocalID& newSize) {
      // Does not set new added blocks to zero
      #ifdef USE_GPU
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         gpuStream_t stream = gpu_getStream();
         vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
         setNewCapacity(newSize,stream);
         parameters->resize((newSize)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true,stream);
         block_data->resize((newSize)*WID3,true,stream);
         #else
         const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
         assert(newSize <= currentCapacity && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back N_blocks).");
         block_data->device_resize((newSize)*WID3,false); //construct=false don't construct or set to zero
         parameters->device_resize((newSize)*BlockParams::N_VELOCITY_BLOCK_PARAMS,false); //construct=false don't construct or set to zero
         #endif
      #else
      setNewCapacity(newSize);
      block_data->resize((newSize)*WID3);
      parameters->resize((newSize)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif
      return true;
   }

   /** Return the number of existing velocity blocks.
    * @return Number of existing velocity blocks.*/
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::size() const {
      return block_data->size()/WID3;
   }

   inline ARCH_HOSTDEV size_t VelocityBlockContainer::sizeInBytes() const {
      const vmesh::LocalID currentSize = block_data->size();
      const vmesh::LocalID parametersSize = parameters->size();
      return currentSize*sizeof(Realf) + parametersSize*sizeof(Real);
   }

   // inline ARCH_HOSTDEV void VelocityBlockContainer::swap(VelocityBlockContainer& vbc) {
   //    block_data->swap(*(vbc.block_data));
   //    parameters->swap(*(vbc.parameters));
   // }

#ifdef DEBUG_VBC
   inline const Realf& VelocityBlockContainer::getData(const vmesh::LocalID& blockLID,const unsigned int& cell) const {
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      bool ok = true;
      if (cell >= WID3) ok = false;
      if (blockLID >= numberOfBlocks) ok = false;
      if (blockLID*WID3+cell >= block_data->size()) ok = false;
      if (ok == false) {
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in getData, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " data->size()=" << block_data->size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }
      return (*block_data)[blockLID*WID3+cell];
   }

   inline const Real& VelocityBlockContainer::getParameters(const vmesh::LocalID& blockLID,const unsigned int& cell) const {
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      bool ok = true;
      if (cell >= BlockParams::N_VELOCITY_BLOCK_PARAMS) ok = false;
      if (blockLID >= numberOfBlocks) ok = false;
      if (blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+cell >= parameters->size()) ok = false;
      if (ok == false) {
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in getParameters, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " parameters->size()=" << parameters->size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }
      return (*parameters)[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+cell];
   }

   inline void VelocityBlockContainer::setData(const vmesh::LocalID& blockLID,const unsigned int& cell,const Realf& value) {
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      bool ok = true;
      if (cell >= WID3) ok = false;
      if (blockLID >= numberOfBlocks) ok = false;
      if (blockLID*WID3+cell >= block_data->size()) ok = false;
      if (ok == false) {
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in setData, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " data->size()=" << block_data->size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }

      (*block_data)[blockLID*WID3+cell] = value;
   }
#endif //debug VBC

} // namespace block_cont

#endif
