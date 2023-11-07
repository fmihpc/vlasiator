/*
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
   #define DEBUG_VBC
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
      void clear();
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
      bool recapacitate(const vmesh::LocalID& capacity);
      ARCH_HOSTDEV bool setSize(const vmesh::LocalID& newSize);
      ARCH_HOSTDEV vmesh::LocalID size() const;
      ARCH_HOSTDEV size_t sizeInBytes() const;
      ARCH_HOSTDEV void swap(VelocityBlockContainer& vbc);

#ifdef USE_GPU // for GPU version
      void gpu_Allocate(vmesh::LocalID size);
      void gpu_Allocate();
      void gpu_prefetchHost(gpuStream_t stream);
      void gpu_prefetchDevice(gpuStream_t stream);
      void gpu_attachToStream(gpuStream_t stream);
      void gpu_detachFromStream();
      void print_addresses();
      void gpu_memAdvise(int device, gpuStream_t stream);
#endif

      #ifdef DEBUG_VBC
      const Realf& getData(const vmesh::LocalID& blockLID,const unsigned int& cell) const;
      const Real& getParameters(const vmesh::LocalID& blockLID,const unsigned int& i) const;
      void setData(const vmesh::LocalID& blockLID,const unsigned int& cell,const Realf& value);
      #endif
    private:
      void exitInvalidLocalID(const vmesh::LocalID& localID,const std::string& funcName) const;
      ARCH_DEV void exitInvalidLocalID(const vmesh::LocalID& localID) const;
      void resize(vmesh::LocalID add=1);

#ifdef USE_GPU
      gpuStream_t attachedStream;
      split::SplitVector<Realf> *block_data;
      split::SplitVector<Real> *parameters;
#else
      std::vector<Realf,aligned_allocator<Realf,WID3> > *block_data;
      std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > *parameters;
#endif

   };

   inline VelocityBlockContainer::VelocityBlockContainer() {
#ifdef USE_GPU
      block_data= new split::SplitVector<Realf>(1);
      parameters= new split::SplitVector<Real>(1);
      attachedStream = 0;
#else
      block_data = new std::vector<Realf,aligned_allocator<Realf,WID3>>(1);
      parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS>>(1);
#endif
      block_data->clear();
      parameters->clear();
   }

   inline VelocityBlockContainer::~VelocityBlockContainer() {
      gpu_destructor();
   }
   inline void VelocityBlockContainer::gpu_destructor() {
      if (block_data) delete block_data;
      if (parameters) delete parameters;
      block_data = NULL;
      parameters = NULL;
      #ifdef USE_GPU
      attachedStream = 0;
      #endif
   }

   inline VelocityBlockContainer::VelocityBlockContainer(const VelocityBlockContainer& other) {
#ifdef USE_GPU
      attachedStream = 0;
      block_data= new split::SplitVector<Realf>(*(other.block_data));
      parameters= new split::SplitVector<Real>(*(other.parameters));
      block_data->reserve(other.block_data->capacity(), true);
      parameters->reserve(other.parameters->capacity(), true);
#else
      block_data = new std::vector<Realf,aligned_allocator<Realf,WID3>>(*(other.block_data));
      parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS>>(*(other.parameters));
      block_data->reserve(other.block_data->capacity());
      parameters->reserve(other.parameters->capacity());
#endif
   }

   inline const VelocityBlockContainer& VelocityBlockContainer::operator=(const VelocityBlockContainer& other) {
      *block_data = *(other.block_data);
      *parameters = *(other.parameters);
      #ifdef USE_GPU
      // gpuStream_t stream = gpu_getStream();
      attachedStream = other.attachedStream;
      block_data->reserve(other.block_data->capacity(), true);
      parameters->reserve(other.parameters->capacity(), true);
      #else
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
      return (block_data->capacity())*sizeof(Realf) + parameters->capacity()*sizeof(Real);
   }

   /** Clears VelocityBlockContainer data and deallocates all memory
    * reserved for velocity blocks.*/
   inline void VelocityBlockContainer::clear() {
#ifdef USE_GPU
      block_data->resize(1,true);
      block_data->shrink_to_fit();
      parameters->resize(1,true);
      parameters->shrink_to_fit();
#else
      std::vector<Realf,aligned_allocator<Realf,WID3> > *dummy_data = new std::vector<Realf,aligned_allocator<Realf,WID3> >(1);
      std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > *dummy_parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> >(1);
      // initialization with zero capacity returns null pointers
      block_data->swap(*dummy_data);
      parameters->swap(*dummy_parameters);
      block_data->clear();
      parameters->clear();
      delete dummy_data;
      delete dummy_parameters;
#endif
   }

   inline ARCH_HOSTDEV void VelocityBlockContainer::move(const vmesh::LocalID& source,const vmesh::LocalID& target) {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      #endif

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
         if (target >= numberOfBlocksP) ok = false;
         if (target >= currentCapacityP) ok = false;
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
   inline void VelocityBlockContainer::gpu_Allocate(vmesh::LocalID size) {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      vmesh::LocalID currentCapacity = info_1[thread_id]->capacity / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
      #endif

      vmesh::LocalID requirement = (size > numberOfBlocks) ? size : numberOfBlocks;
      if (currentCapacity > BLOCK_ALLOCATION_FACTOR * requirement) {
         return; // Still have enough buffer
      }

      if (numberOfBlocks > size) {
         std::cerr<<"Error in "<<__FILE__<<" line "<<__LINE__<<": attempting to allocate less than numberOfBlocks"<<std::endl;
         abort();
      }
      if (numberOfBlocks > BLOCK_ALLOCATION_FACTOR * size) {
         std::cerr<<"Warning in "<<__FILE__<<" line "<<__LINE__<<": attempting to allocate less than safety margins"<<std::endl;
      }
      // Passing eco flag = true to reserve tells splitvector we manage padding manually.
      currentCapacity = BLOCK_ALLOCATION_PADDING * requirement;
      block_data->reserve(currentCapacity * WID3, true);
      parameters->reserve(currentCapacity * BlockParams::N_VELOCITY_BLOCK_PARAMS, true);
      return;
   }

   inline void VelocityBlockContainer::gpu_Allocate() {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      vmesh::LocalID currentCapacity = info_1[thread_id]->capacity / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
      #endif

      if (currentCapacity > BLOCK_ALLOCATION_FACTOR * numberOfBlocks) {
         return; // Still have enough buffer
      }
      // Passing eco flag = true to reserve tells splitvector we manage padding manually.
      currentCapacity = BLOCK_ALLOCATION_PADDING * numberOfBlocks;
      block_data->resize(currentCapacity * WID3, true);
      parameters->resize(currentCapacity * BlockParams::N_VELOCITY_BLOCK_PARAMS, true);
      return;
   }

   inline void VelocityBlockContainer::gpu_prefetchHost(gpuStream_t stream=0) {
      if (stream==0) {
         stream = gpu_getStream();
      }
      block_data->optimizeCPU(gpu_getStream());
      parameters->optimizeCPU(gpu_getStream());
      return;
   }

   inline void VelocityBlockContainer::gpu_prefetchDevice(gpuStream_t stream=0) {
      if (stream==0) {
         stream = gpu_getStream();
      }
      block_data->optimizeGPU(stream);
      parameters->optimizeGPU(stream);
      return;
   }

   inline void VelocityBlockContainer::gpu_memAdvise(int device, gpuStream_t stream) {
      // AMD advise is slow
      // // int device = gpu_getDevice();
      // block_data->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      // parameters->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      // block_data->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      // parameters->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      return;
   }

   inline void VelocityBlockContainer::gpu_attachToStream(gpuStream_t stream) {
      // Return if attaching is not needed
      if (!needAttachedStreams) {
         return;
      }
      // Attach unified memory regions to streams
      gpuStream_t newStream;
      if (stream==0) {
         newStream = gpu_getStream();
      } else {
         newStream = stream;
      }
      if (newStream == attachedStream) {
         return;
      } else {
         attachedStream = newStream;
      }
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,block_data, 0,gpuMemAttachSingle) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,parameters, 0,gpuMemAttachSingle) );
      block_data->streamAttach(attachedStream);
      parameters->streamAttach(attachedStream);
      return;
   }
   inline void VelocityBlockContainer::gpu_detachFromStream() {
      // Return if attaching is not needed
      if (!needAttachedStreams) {
         return;
      }
      if (attachedStream == 0) {
         // Already detached
         return;
      }
      attachedStream = 0;
      // Detach unified memory regions from streams
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,block_data, 0,gpuMemAttachGlobal) );
      CHK_ERR( gpuStreamAttachMemAsync(attachedStream,parameters, 0,gpuMemAttachGlobal) );
      block_data->streamAttach(0,gpuMemAttachGlobal);
      parameters->streamAttach(0,gpuMemAttachGlobal);
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
         #endif
      #endif
      return parameters->data() + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   inline ARCH_HOSTDEV void VelocityBlockContainer::pop() {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      #endif

      if (numberOfBlocks == 0) return;
      block_data->erase(block_data->begin() + WID3*(numberOfBlocks-1),
                        block_data->begin() + WID3*(numberOfBlocks));
      parameters->erase(parameters->begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS*(numberOfBlocks-1),
                        parameters->begin() + BlockParams::N_VELOCITY_BLOCK_PARAMS*(numberOfBlocks));
   }

   /** Grows the size of the VBC, does not touch data */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back() {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      #endif

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
      #ifdef USE_GPU
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         resize();
         block_data->resize((numberOfBlocks+1)*WID3,true);
         parameters->resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true);
         #else
         const vmesh::LocalID currentCapacityD = block_data->capacity()/WID3;
         if (newIndex >= currentCapacityD) {
            assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back).");
         }
         block_data->device_resize((numberOfBlocks+1)*WID3);
         parameters->device_resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
         #endif
      #else
      resize();
      block_data->resize((numberOfBlocks+1)*WID3);
      parameters->resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif
      return newIndex;
   }

   /** Grows the size of the VBC, sets data to zero */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back_and_zero() {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      #endif

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

      #ifdef USE_GPU
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         resize();
         block_data->resize((numberOfBlocks+1)*WID3,true);
         parameters->resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true);
         #else
         const vmesh::LocalID currentCapacityD = block_data->capacity()/WID3;
         if (newIndex >= currentCapacityD) {
            assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back_and_zero).");
         }
         block_data->device_resize((numberOfBlocks+1)*WID3);
         parameters->device_resize((numberOfBlocks+1)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
         #endif
      #else
      resize();
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
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      const vmesh::LocalID currentCapacity = info_1[thread_id]->capacity / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
      #endif

      const vmesh::LocalID newIndex = numberOfBlocks;

      #ifdef USE_GPU
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         resize(N_blocks);
         block_data->resize((numberOfBlocks+N_blocks)*WID3,true);
         parameters->resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true);
         #else
         if (newIndex + N_blocks >= currentCapacity-1) {
            assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back N_blocks).");
         }
         block_data->device_resize((numberOfBlocks+N_blocks)*WID3);
         parameters->device_resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
         #endif
      #else
      resize(N_blocks);
      block_data->resize((numberOfBlocks+N_blocks)*WID3);
      parameters->resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif
      return newIndex;
   }

   /** Grows the size of the VBC, sets data to zero */
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back_and_zero(const uint32_t& N_blocks) {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      const vmesh::LocalID currentCapacity = info_1[thread_id]->capacity / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
      #endif

      const vmesh::LocalID newIndex = numberOfBlocks;

      #ifdef USE_GPU
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         resize(N_blocks);
         block_data->resize((numberOfBlocks+N_blocks)*WID3,true);
         parameters->resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true);
         #else
         if (newIndex + N_blocks >= currentCapacity-1) {
            assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back_and_zero N_blocks).");
         }
         block_data->device_resize((numberOfBlocks+N_blocks)*WID3);
         parameters->device_resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
         #endif
      #else
      resize(N_blocks);
      block_data->resize((numberOfBlocks+N_blocks)*WID3);
      parameters->resize((numberOfBlocks+N_blocks)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Clear velocity block data to zero values
      CHK_ERR( gpuMemsetAsync(&(block_data[newIndex*WID3]), 0, WID3*N_blocks*sizeof(Realf), stream) );
      CHK_ERR( gpuMemsetAsync(&(parameters[newIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS]), 0, BlockParams::N_VELOCITY_BLOCK_PARAMS*N_blocks*sizeof(Real), stream) );
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

   inline bool VelocityBlockContainer::recapacitate(const vmesh::LocalID& newCapacity) {
      // Note: No longer ever recapacitates down in size.
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      const vmesh::LocalID currentCapacity = info_1[thread_id]->capacity / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
      #endif

      if (newCapacity < numberOfBlocks) {
         std::cerr<<" ERROR! Trying to recapacitate to "<<newCapacity<<" when VBC already contains "<<numberOfBlocks<<" blocks!"<<std::endl;
         return false;
      }
      #ifdef USE_GPU
      block_data->reserve(newCapacity*WID3, true);
      parameters->reserve(newCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS, true);
      #else
      block_data->reserve(newCapacity*WID3);
      parameters->reserve(newCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      #endif

      return true;
   }

   inline void VelocityBlockContainer::resize(const vmesh::LocalID add) {
      //This actually only alters capacity, not size.
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      vmesh::LocalID currentCapacity = info_1[thread_id]->capacity / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
      #endif

#ifdef USE_GPU
      if ((numberOfBlocks+add) * BLOCK_ALLOCATION_FACTOR >= currentCapacity) {
         // Resize so that free space is current * block_allocation_padding blocks,
         // and at least two in case of having zero blocks.
         // The order of velocity blocks is unaltered.
         currentCapacity = 2 + (numberOfBlocks+add) * BLOCK_ALLOCATION_PADDING;
         // reallocate() doesn't change vector size, would lead to data loss
         // when further reallocating due to lack of data copy
         // Passing eco flag = true to resize tells splitvector we manage padding manually.
         block_data->reserve(currentCapacity*WID3, true);
         parameters->reserve(currentCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS, true);
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         if ((attachedStream != 0)&&(needAttachedStreams)) {
            block_data->streamAttach(attachedStream);
            parameters->streamAttach(attachedStream);
         }
         gpuStream_t stream = gpu_getStream();
         CHK_ERR( gpuStreamSynchronize(stream) );
         block_data->optimizeGPU(stream);
         parameters->optimizeGPU(stream);
         // int device = gpu_getDevice();
         // block_data->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         // parameters->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
         // block_data->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
         // parameters->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
         #endif
      }
#else
      if ((numberOfBlocks+add) >= currentCapacity) {
         // Resize so that free space is current * block_allocation_factor blocks,
         // and at least two in case of having zero blocks.
         // The order of velocity blocks is unaltered.
         currentCapacity = 2 + (numberOfBlocks+add) * BLOCK_ALLOCATION_FACTOR;
         block_data->reserve(currentCapacity*WID3);
         parameters->reserve(currentCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      }
#endif
   }

   inline ARCH_HOSTDEV bool VelocityBlockContainer::setSize(const vmesh::LocalID& newSize) {
      #ifdef USE_GPU
         #if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
         block_data->resize((newSize)*WID3,true);
         parameters->resize((newSize)*BlockParams::N_VELOCITY_BLOCK_PARAMS,true);
         // Ensure buffer afterwards
         resize();
         #else
         const vmesh::LocalID currentCapacity = block_data->capacity()/WID3;
         if (newSize > currentCapacity) {
            assert(0 && "ERROR! Attempting to grow block container on-device beyond capacity (::push_back N_blocks).");
         }
         block_data->device_resize((newSize)*WID3);
         parameters->device_resize((newSize)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
         #endif
      #else
      block_data->resize((newSize)*WID3);
      parameters->resize((newSize)*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      // Ensure buffer afterwards
      resize();
      #endif

      return true;
   }

   /** Return the number of existing velocity blocks.
    * @return Number of existing velocity blocks.*/
   inline ARCH_HOSTDEV vmesh::LocalID VelocityBlockContainer::size() const {
      #if defined(USE_GPU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
      // Host-side non-pagefaulting approach
      const uint thread_id = gpu_getThread();
      gpuStream_t stream = gpuStreamList[thread_id];
      block_data->copyMetadata(info_1[thread_id],stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      const vmesh::LocalID numberOfBlocks = info_1[thread_id]->size / WID3;
      #else
      const vmesh::LocalID numberOfBlocks = block_data->size()/WID3;
      #endif
      return numberOfBlocks;
   }

   inline ARCH_HOSTDEV size_t VelocityBlockContainer::sizeInBytes() const {
      return block_data->capacity()*sizeof(Realf) + parameters->capacity()*sizeof(Real);
   }

   inline ARCH_HOSTDEV void VelocityBlockContainer::swap(VelocityBlockContainer& vbc) {
      block_data->swap(*(vbc.block_data));
      parameters->swap(*(vbc.parameters));
   }

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
