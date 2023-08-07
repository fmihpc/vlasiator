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

#ifdef DEBUG_VBC
#include <sstream>
#endif

#ifdef __CUDACC__
#endif

#ifdef USE_CUDA
   #include "cuda_context.cuh"
   // Place block data and parameters inside splitvectors utilizing unified memory
   #include "include/splitvector/splitvec.h"
#else
   #define CUDA_HOSTDEV
   // CUDA allocation factors are stored in cuda_context.cuh
   static const double BLOCK_ALLOCATION_FACTOR = 1.1;
#endif

using namespace std;

namespace vmesh {

#ifdef USE_CUDA
   class VelocityBlockContainer : public Managed {
#else
   class VelocityBlockContainer {
#endif
   public:

      VelocityBlockContainer();
      ~VelocityBlockContainer();
      VelocityBlockContainer(const VelocityBlockContainer& other);
      const VelocityBlockContainer& operator=(const VelocityBlockContainer& other);

      CUDA_HOSTDEV vmesh::LocalID capacity() const;
      CUDA_HOSTDEV size_t capacityInBytes() const;
      void clear();
      CUDA_HOSTDEV void copy(const vmesh::LocalID& source,const vmesh::LocalID& target);
      CUDA_HOSTDEV static double getBlockAllocationFactor();
      CUDA_HOSTDEV Realf* getData();
      CUDA_HOSTDEV const Realf* getData() const;
      CUDA_HOSTDEV Realf* getData(const vmesh::LocalID& blockLID);
      CUDA_HOSTDEV const Realf* getData(const vmesh::LocalID& blockLID) const;
      CUDA_HOSTDEV Real* getParameters();
      CUDA_HOSTDEV const Real* getParameters() const;
      CUDA_HOSTDEV Real* getParameters(const vmesh::LocalID& blockLID);
      CUDA_HOSTDEV const Real* getParameters(const vmesh::LocalID& blockLID) const;
      CUDA_HOSTDEV void pop();
      CUDA_HOSTDEV vmesh::LocalID push_back();
      CUDA_HOSTDEV vmesh::LocalID push_back(const uint32_t& N_blocks);
      bool recapacitate(const vmesh::LocalID& capacity);
      CUDA_HOSTDEV bool setSize(const vmesh::LocalID& newSize);
      CUDA_HOSTDEV vmesh::LocalID size() const;
      CUDA_HOSTDEV size_t sizeInBytes() const;
      CUDA_HOSTDEV void swap(VelocityBlockContainer& vbc);

#ifdef USE_CUDA // for CUDA version
      CUDA_HOSTDEV vmesh::LocalID dev_push_back();
      CUDA_HOSTDEV vmesh::LocalID dev_push_back(const uint32_t& N_blocks);
      void dev_Allocate(vmesh::LocalID size);
      void dev_Allocate();
      void dev_prefetchHost();
      void dev_prefetchDevice();
      void dev_attachToStream(cudaStream_t stream=0);
      void dev_detachFromStream();
      void dev_memAdvise(int device);
#endif

      #ifdef DEBUG_VBC
      const Realf& getData(const vmesh::LocalID& blockLID,const unsigned int& cell) const;
      const Real& getParameters(const vmesh::LocalID& blockLID,const unsigned int& i) const;
      void setData(const vmesh::LocalID& blockLID,const unsigned int& cell,const Realf& value);
      #endif
    private:
      void exitInvalidLocalID(const vmesh::LocalID& localID,const std::string& funcName) const;
      void resize();

      vmesh::LocalID currentCapacity;
      vmesh::LocalID numberOfBlocks;

#ifdef USE_CUDA
      cudaStream_t attachedStream;
      split::SplitVector<Realf> *block_data;
      split::SplitVector<Real> *parameters;
#else
      std::vector<Realf,aligned_allocator<Realf,WID3> > *block_data;
      std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > *parameters;
#endif

   };

   inline VelocityBlockContainer::VelocityBlockContainer() : currentCapacity {0}, numberOfBlocks {0} {
      // initialization with zero capacity would return null pointers
      const uint capacity = currentCapacity > 0 ? currentCapacity : 1;
#ifdef USE_CUDA
      block_data= new split::SplitVector<Realf>(capacity);
      parameters= new split::SplitVector<Real>(capacity);
      attachedStream = 0;
#else
      block_data = new std::vector<Realf,aligned_allocator<Realf,WID3>>(capacity);
      parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS>>(capacity);
#endif
      block_data->clear();
      parameters->clear();
   }

   inline VelocityBlockContainer::~VelocityBlockContainer() {
      delete block_data;
      delete parameters;
   }

   inline VelocityBlockContainer::VelocityBlockContainer(const VelocityBlockContainer& other) {
#ifdef USE_CUDA
      attachedStream = 0;
      block_data= new split::SplitVector<Realf>(*(other.block_data));
      parameters= new split::SplitVector<Real>(*(other.parameters));
#else
      block_data = new std::vector<Realf,aligned_allocator<Realf,WID3>>(*(other.block_data));
      parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS>>(*(other.parameters));
#endif
      currentCapacity = other.currentCapacity;
      numberOfBlocks = other.numberOfBlocks;
   }

   inline const VelocityBlockContainer& VelocityBlockContainer::operator=(const VelocityBlockContainer& other) {
      // Delete old vectors
      delete block_data;
      delete parameters;
#ifdef USE_CUDA
      attachedStream = 0;
      block_data= new split::SplitVector<Realf>(*(other.block_data));
      parameters= new split::SplitVector<Real>(*(other.parameters));
#else
      block_data = new std::vector<Realf,aligned_allocator<Realf,WID3>>(*(other.block_data));
      parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS>>(*(other.parameters));
#endif
      currentCapacity = other.currentCapacity;
      numberOfBlocks = other.numberOfBlocks;
      return *this;
   }

   inline CUDA_HOSTDEV vmesh::LocalID VelocityBlockContainer::capacity() const {
      return currentCapacity;
   }

   inline CUDA_HOSTDEV size_t VelocityBlockContainer::capacityInBytes() const {
      return (block_data->capacity())*sizeof(Realf) + parameters->capacity()*sizeof(Real);
   }

   /** Clears VelocityBlockContainer data and deallocates all memory
    * reserved for velocity blocks.*/
   inline void VelocityBlockContainer::clear() {
#ifdef USE_CUDA
      block_data->clear();
      block_data->shrink_to_fit();
      parameters->clear();
      parameters->shrink_to_fit();
 #else
      std::vector<Realf,aligned_allocator<Realf,WID3> > *dummy_data = new std::vector<Realf,aligned_allocator<Realf,WID3> >(1);
      std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > *dummy_parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> >(1);
      // initialization with zero capacity returns null pointers
      block_data->swap(*dummy_data);
      parameters->swap(*dummy_parameters);
      block_data->clear();
      parameters->clear();
      int device = cuda_getDevice();
      block_data->memAdvise(cudaMemAdviseSetPreferredLocation,device);
      parameters->memAdvise(cudaMemAdviseSetPreferredLocation,device);
      block_data->memAdvise(cudaMemAdviseSetAccessedBy,device);
      parameters->memAdvise(cudaMemAdviseSetAccessedBy,device);
      delete dummy_data;
      delete dummy_parameters;
#endif
      currentCapacity = 0;
      numberOfBlocks = 0;
   }

   inline CUDA_HOSTDEV void VelocityBlockContainer::copy(const vmesh::LocalID& source,const vmesh::LocalID& target) {
      #ifdef DEBUG_VBC
         bool ok = true;
         if (source >= numberOfBlocks) ok = false;
         if (source >= currentCapacity) ok = false;
         if (target >= numberOfBlocks) ok = false;
         if (target >= currentCapacity) ok = false;
         if (numberOfBlocks >= currentCapacity) ok = false;
         if (source != numberOfBlocks-1) ok = false;
         if (source*WID3+WID3-1 >= block_data->size()) ok = false;
         if (source*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::N_VELOCITY_BLOCK_PARAMS-1 >= parameters->size()) ok = false;
         if (target*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::N_VELOCITY_BLOCK_PARAMS-1 >= parameters->size()) ok = false;
         if (parameters->size()/BlockParams::N_VELOCITY_BLOCK_PARAMS != block_data->size()/WID3) ok = false;
         if (ok == false) {
            std::stringstream ss;
            ss << "VBC ERROR: invalid source LID=" << source << " in copy, target=" << target << " #blocks=" << numberOfBlocks << " capacity=" << currentCapacity << std::endl;
            ss << "or sizes are wrong, data.size()=" << block_data->size() << " parameters->size()=" << parameters->size() << std::endl;
            std::cerr << ss.str();
            sleep(1);
            exit(1);
         }
      #endif

      for (unsigned int i=0; i<WID3; ++i) {
         (*block_data)[target*WID3+i] = (*block_data)[source*WID3+i];
      }
      for (int i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) {
         (*parameters)[target*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = (*parameters)[source*BlockParams::N_VELOCITY_BLOCK_PARAMS+i];
      }
   }

   inline void VelocityBlockContainer::exitInvalidLocalID(const vmesh::LocalID& localID,const std::string& funcName) const {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);

      #ifdef DEBUG_VBC
      std::stringstream ss;
      ss << "Process " << rank << ' ';
      ss << "Invalid localID " << localID << " used in function '" << funcName << "' max allowed value is " << numberOfBlocks << std::endl;
      std::cerr << ss.str();
      sleep(1);
      #endif
      exit(1);
   }

   inline CUDA_HOSTDEV double VelocityBlockContainer::getBlockAllocationFactor() {
      return BLOCK_ALLOCATION_FACTOR;
   }

   inline CUDA_HOSTDEV Realf* VelocityBlockContainer::getData() {
      return block_data->data();
   }

   inline CUDA_HOSTDEV const Realf* VelocityBlockContainer::getData() const {
      return block_data->data();
   }

   inline CUDA_HOSTDEV Realf* VelocityBlockContainer::getData(const vmesh::LocalID& blockLID) {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"getData");
         if (blockLID >= block_data->size()/WID3) exitInvalidLocalID(blockLID,"const getData const");
      #endif
      return block_data->data() + blockLID*WID3;
   }

   inline CUDA_HOSTDEV const Realf* VelocityBlockContainer::getData(const vmesh::LocalID& blockLID) const {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"const getData const");
         if (blockLID >= block_data->size()/WID3) exitInvalidLocalID(blockLID,"const getData const");
      #endif
      return block_data->data() + blockLID*WID3;
   }

#ifdef USE_CUDA
   inline void VelocityBlockContainer::dev_Allocate(vmesh::LocalID size) {
      vmesh::LocalID requirement = (size > numberOfBlocks) ? size : numberOfBlocks;
      if (block_data->capacity() > BLOCK_ALLOCATION_FACTOR * requirement) {
         return; // Still have enough buffer
      }

      if (numberOfBlocks > size) {
         std::cerr<<"Error in "<<__FILE__<<" line "<<__LINE__<<": attempting to allocate less than numberOfBlocks"<<std::endl;
         abort();
      }
      if (numberOfBlocks > BLOCK_ALLOCATION_FACTOR * size) {
         std::cerr<<"Warning in "<<__FILE__<<" line "<<__LINE__<<": attempting to allocate less than safety margins"<<std::endl;
      }
      block_data->reallocate(BLOCK_ALLOCATION_PADDING * size);
      parameters->reallocate(BLOCK_ALLOCATION_PADDING * size);
      return;
   }

   inline void VelocityBlockContainer::dev_Allocate() {
      if (block_data->capacity() > BLOCK_ALLOCATION_FACTOR * numberOfBlocks) {
         return; // Still have enough buffer
      }
      block_data->reallocate(BLOCK_ALLOCATION_PADDING * numberOfBlocks);
      parameters->reallocate(BLOCK_ALLOCATION_PADDING * numberOfBlocks);
      return;
   }

   inline void VelocityBlockContainer::dev_prefetchHost() {
      //if (numberOfBlocks==0) return; // This size check in itself causes a page fault
      block_data->optimizeCPU(cuda_getStream());
      parameters->optimizeCPU(cuda_getStream());
      return;
   }

   inline void VelocityBlockContainer::dev_prefetchDevice() {
      //if (numberOfBlocks==0) return; // This size check in itself causes a page fault
      block_data->optimizeGPU(cuda_getStream());
      parameters->optimizeGPU(cuda_getStream());
      return;
   }

   inline void VelocityBlockContainer::dev_memAdvise(int device) {
      // int device = cuda_getDevice();
      block_data->memAdvise(cudaMemAdviseSetPreferredLocation,device);
      parameters->memAdvise(cudaMemAdviseSetPreferredLocation,device);
      block_data->memAdvise(cudaMemAdviseSetAccessedBy,device);
      parameters->memAdvise(cudaMemAdviseSetAccessedBy,device);
      return;
   }

   inline void VelocityBlockContainer::dev_attachToStream(cudaStream_t stream) {
      // Return if attaching is not needed
      if (!needAttachedStreams) {
         return;
      }
      // Attach unified memory regions to streams
      cudaStream_t newStream;
      if (stream==0) {
         newStream = cuda_getStream();
      } else {
         newStream = stream;
      }
      if (newStream == attachedStream) {
         return;
      } else {
         attachedStream = newStream;
      }
      HANDLE_ERROR( cudaStreamAttachMemAsync(attachedStream,this, 0,cudaMemAttachSingle) );
      HANDLE_ERROR( cudaStreamAttachMemAsync(attachedStream,block_data, 0,cudaMemAttachSingle) );
      HANDLE_ERROR( cudaStreamAttachMemAsync(attachedStream,parameters, 0,cudaMemAttachSingle) );
      block_data->streamAttach(attachedStream);
      parameters->streamAttach(attachedStream);
      return;
   }
   inline void VelocityBlockContainer::dev_detachFromStream() {
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
      HANDLE_ERROR( cudaStreamAttachMemAsync(attachedStream,this, 0,cudaMemAttachGlobal) );
      HANDLE_ERROR( cudaStreamAttachMemAsync(attachedStream,block_data, 0,cudaMemAttachGlobal) );
      HANDLE_ERROR( cudaStreamAttachMemAsync(attachedStream,parameters, 0,cudaMemAttachGlobal) );
      block_data->streamAttach(0,cudaMemAttachGlobal);
      parameters->streamAttach(0,cudaMemAttachGlobal);
      return;
   }
#endif

   inline CUDA_HOSTDEV Real* VelocityBlockContainer::getParameters() {
      return parameters->data();
   }

   inline CUDA_HOSTDEV const Real* VelocityBlockContainer::getParameters() const {
      return parameters->data();
   }

   inline CUDA_HOSTDEV Real* VelocityBlockContainer::getParameters(const vmesh::LocalID& blockLID) {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"getParameters");
         if (blockLID >= parameters->size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) exitInvalidLocalID(blockLID,"getParameters");
      #endif
      return parameters->data() + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   inline CUDA_HOSTDEV const Real* VelocityBlockContainer::getParameters(const vmesh::LocalID& blockLID) const {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"const getParameters const");
         if (blockLID >= parameters->size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) exitInvalidLocalID(blockLID,"getParameters");
      #endif
      return parameters->data() + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   inline CUDA_HOSTDEV void VelocityBlockContainer::pop() {
      if (numberOfBlocks == 0) return;
      --numberOfBlocks;
   }

   inline CUDA_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back() {
      vmesh::LocalID newIndex = numberOfBlocks;
      if (newIndex >= currentCapacity) {
#ifdef USE_CUDA
#pragma nv_diag_suppress 20011,20014
#endif
         resize();
      }

      #ifdef DEBUG_VBC
      if (newIndex >= block_data->size()/WID3 || newIndex >= parameters->size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) {
         std::stringstream ss;
         ss << "VBC ERROR in push_back, LID=" << newIndex << " for new block is out of bounds" << std::endl;
         ss << "\t data.size()=" << block_data->size()  << " parameters->size()=" << parameters->size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }
      #endif

      // Clear velocity block data to zero values
      for (size_t i=0; i<WID3; ++i) {
         (*block_data)[newIndex*WID3+i] = 0.0;
      }
      for (size_t i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) {
         (*parameters)[newIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = 0.0;
      }
      ++numberOfBlocks;
      return newIndex;
   }

   inline CUDA_HOSTDEV vmesh::LocalID VelocityBlockContainer::push_back(const uint32_t& N_blocks) {
      const vmesh::LocalID newIndex = numberOfBlocks;
      numberOfBlocks += N_blocks;
      if (numberOfBlocks > currentCapacity) {
#ifdef USE_CUDA
#pragma nv_diag_suppress 20014
#endif
         resize();
      }

      // Clear velocity block data to zero values
      for (size_t i=0; i<WID3*N_blocks; ++i) {
         (*block_data)[newIndex*WID3+i] = 0.0;
      }
      for (size_t i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS*N_blocks; ++i) {
         (*parameters)[newIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = 0.0;
      }

      return newIndex;
   }

/** PUSH BACK METHODS FOR GPU CODE - bookkeeping only
    Zeroing out block data and parameters must be done from caller kernel **/
#ifdef USE_CUDA
   inline CUDA_HOSTDEV vmesh::LocalID VelocityBlockContainer::dev_push_back() {
      vmesh::LocalID newIndex = numberOfBlocks;
      if (newIndex >= currentCapacity) {
         printf("Error: pushing back to VBC from device without enough capacity!\n");
         return 0;
      }
      ++numberOfBlocks;
      return newIndex;
   }
   inline CUDA_HOSTDEV vmesh::LocalID VelocityBlockContainer::dev_push_back(const uint32_t& N_blocks) {
      const vmesh::LocalID newIndex = numberOfBlocks;
      numberOfBlocks += N_blocks;
      if (numberOfBlocks > currentCapacity) {
         printf("Error: pushing back to VBC from device without enough capacity!\n");
         return 0;
      }
      return newIndex;
   }
#endif

   inline bool VelocityBlockContainer::recapacitate(const vmesh::LocalID& newCapacity) {
      if (newCapacity < numberOfBlocks) return false;
      {
#ifdef USE_CUDA
         split::SplitVector<Realf> *dummy_data = new split::SplitVector<Realf>(newCapacity*WID3);
#else
         std::vector<Realf,aligned_allocator<Realf,WID3> > *dummy_data = new std::vector<Realf,aligned_allocator<Realf,WID3> >(newCapacity*WID3);
#endif
         for (size_t i=0; i<numberOfBlocks*WID3; ++i) (*dummy_data)[i] = (*block_data)[i];
         dummy_data->swap(*block_data);
         delete dummy_data;
      }
      {
#ifdef USE_CUDA
         split::SplitVector<Real> *dummy_parameters = new split::SplitVector<Real>(newCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
#else
         std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > *dummy_parameters = new std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> >(newCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
#endif
         for (size_t i=0; i<numberOfBlocks*BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) (*dummy_parameters)[i] = (*parameters)[i];
         dummy_parameters->swap(*parameters);
         delete dummy_parameters;
      }
      #ifdef USE_CUDA
      int device = cuda_getDevice();
      block_data->memAdvise(cudaMemAdviseSetPreferredLocation,device);
      parameters->memAdvise(cudaMemAdviseSetPreferredLocation,device);
      block_data->memAdvise(cudaMemAdviseSetAccessedBy,device);
      parameters->memAdvise(cudaMemAdviseSetAccessedBy,device);
      #endif
      currentCapacity = newCapacity;
   return true;
   }

   inline void VelocityBlockContainer::resize() {
#ifdef USE_CUDA
      if ((numberOfBlocks+1)*BLOCK_ALLOCATION_FACTOR >= currentCapacity) {
         // Resize so that free space is block_allocation_chunk blocks,
         // and at least two in case of having zero blocks.
         // The order of velocity blocks is unaltered.
         currentCapacity = 2 + numberOfBlocks * BLOCK_ALLOCATION_PADDING;
         // reallocate() doesn't change vector size, leads to data loss when further reallocating due to lack of data copy
         // Passing eco flag = true to resize tells splitvector we manage padding manually.
         block_data->resize(currentCapacity*WID3, true);
         parameters->resize(currentCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS, true);
#ifndef __CUDA_ARCH__
         if ((attachedStream != 0)&&(needAttachedStreams)) {
            block_data->streamAttach(attachedStream);
            parameters->streamAttach(attachedStream);
         }
#endif
         cudaStream_t stream = cuda_getStream();
         HANDLE_ERROR( cudaStreamSynchronize(stream) );
         block_data->optimizeGPU(stream);
         parameters->optimizeGPU(stream);
         int device = cuda_getDevice();
         block_data->memAdvise(cudaMemAdviseSetPreferredLocation,device);
         parameters->memAdvise(cudaMemAdviseSetPreferredLocation,device);
         block_data->memAdvise(cudaMemAdviseSetAccessedBy,device);
         parameters->memAdvise(cudaMemAdviseSetAccessedBy,device);
      }
#else
      if ((numberOfBlocks+1) >= currentCapacity) {
         // Resize so that free space is block_allocation_chunk blocks,
         // and at least two in case of having zero blocks.
         // The order of velocity blocks is unaltered.
         currentCapacity = 2 + numberOfBlocks * BLOCK_ALLOCATION_FACTOR;
         block_data->resize(currentCapacity*WID3);
         parameters->resize(currentCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      }
#endif
   }

   inline CUDA_HOSTDEV bool VelocityBlockContainer::setSize(const vmesh::LocalID& newSize) {
      numberOfBlocks = newSize;
      if (newSize > currentCapacity) {
#ifdef USE_CUDA
#pragma nv_diag_suppress 20014
#endif
         resize();
      }
      return true;
   }

   /** Return the number of existing velocity blocks.
    * @return Number of existing velocity blocks.*/
   inline CUDA_HOSTDEV vmesh::LocalID VelocityBlockContainer::size() const {
      return numberOfBlocks;
   }

   inline CUDA_HOSTDEV size_t VelocityBlockContainer::sizeInBytes() const {
      return block_data->size()*sizeof(Realf) + parameters->size()*sizeof(Real);
   }

   inline CUDA_HOSTDEV void VelocityBlockContainer::swap(VelocityBlockContainer& vbc) {
      block_data->swap(*(vbc.block_data));
      parameters->swap(*(vbc.parameters));

      vmesh::LocalID dummy = currentCapacity;
      currentCapacity = vbc.currentCapacity;
      vbc.currentCapacity = dummy;

      dummy = numberOfBlocks;
      numberOfBlocks = vbc.numberOfBlocks;
      vbc.numberOfBlocks = dummy;
   }

   #ifdef DEBUG_VBC

   inline const Realf& VelocityBlockContainer::getData(const vmesh::LocalID& blockLID,const unsigned int& cell) const {
      bool ok = true;
      if (cell >= WID3) ok = false;
      if (blockLID >= numberOfBlocks) ok = false;
      if (blockLID*WID3+cell >= block_data->size()) ok = false;
      if (ok == false) {
         #ifdef DEBUG_VBC
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in getData, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " data.size()=" << block_data->size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         #endif
         exit(1);
      }

      return (*block_data)[blockLID*WID3+cell];
   }

   inline const Real& VelocityBlockContainer::getParameters(const vmesh::LocalID& blockLID,const unsigned int& cell) const {
      bool ok = true;
      if (cell >= BlockParams::N_VELOCITY_BLOCK_PARAMS) ok = false;
      if (blockLID >= numberOfBlocks) ok = false;
      if (blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+cell >= parameters->size()) ok = false;
      if (ok == false) {
         #ifdef DEBUG_VBC
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in getParameters, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " parameters->size()=" << parameters->size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         #endif
         exit(1);
      }

      return (*parameters)[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+cell];
   }

   inline void VelocityBlockContainer::setData(const vmesh::LocalID& blockLID,const unsigned int& cell,const Realf& value) {
      bool ok = true;
      if (cell >= WID3) ok = false;
      if (blockLID >= numberOfBlocks) ok = false;
      if (blockLID*WID3+cell >= block_data->size()) ok = false;
      if (ok == false) {
         #ifdef DEBUG_VBC
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in setData, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " data.size()=" << block_data->size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         #endif
         exit(1);
      }

      (*block_data)[blockLID*WID3+cell] = value;
   }

#endif //debug VBC

} // namespace block_cont

#endif
