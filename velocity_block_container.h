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

#ifndef VELOCITY_BLOCK_CONTAINER_H
#define VELOCITY_BLOCK_CONTAINER_H

#include <vector>
#include <stdio.h>

#include "common.h"
#include "unistd.h"

#ifdef DEBUG_VBC
   #include <sstream>
#endif

#ifdef USE_CUDA
   #include "cuda_context.cuh"
#endif

namespace vmesh {

   static const double BLOCK_ALLOCATION_FACTOR = 1.1;
   static const double CUDA_BLOCK_SAFECTY_FACTOR = 1.6;
   static const double CUDA_BLOCK_ALLOCATION_FACTOR = 2.0;

   template<typename LID>
   class VelocityBlockContainer {
    public:

      VelocityBlockContainer();
#ifdef USE_CUDA
      ~VelocityBlockContainer();
#endif
      LID capacity() const;
      size_t capacityInBytes() const;
      void clear();
      void copy(const LID& source,const LID& target);
      static double getBlockAllocationFactor();
      Realf* getData();
      const Realf* getData() const;
      Realf* getData(const LID& blockLID);
      const Realf* getData(const LID& blockLID) const;
      Realf* getNullData();
      Real* getParameters();
      const Real* getParameters() const;
      Real* getParameters(const LID& blockLID);
      const Real* getParameters(const LID& blockLID) const;
      void pop();
      LID push_back();
      LID push_back(const uint32_t& N_blocks);
      bool recapacitate(const LID& capacity);
      bool setSize(const LID& newSize);
      LID size() const;
      size_t sizeInBytes() const;
      void swap(VelocityBlockContainer& vbc);

#ifdef USE_CUDA // for CUDA version
      Realf* dev_getData();
      const Realf* dev_getData() const;
      Realf* dev_getData(const LID& blockLID);
      const Realf* dev_getData(const LID& blockLID) const;
      Real* dev_getParameters();
      const Real* dev_getParameters() const;
      Real* dev_getParameters(const LID& blockLID);
      const Real* dev_getParameters(const LID& blockLID) const;

      void dev_Deallocate();
      void dev_Allocate(LID size);
      void dev_Allocate();
      void dev_syncBlocksToHost();
      void dev_syncBlocksToDevice();
      void dev_syncParametersToHost();
      void dev_syncParametersToDevice();
      void dev_pinBlocks();
      void dev_unpinBlocks();
      void dev_pinParameters();
      void dev_unpinParameters();
      // Also add CUDA-capable version of vmesh when CUDA openhashmap is available

      bool dev_needsUpdatingBlocks;
      bool dev_needsUpdatingParameters;
#endif

      #ifdef DEBUG_VBC
      const Realf& getData(const LID& blockLID,const unsigned int& cell) const;
      const Real& getParameters(const LID& blockLID,const unsigned int& i) const;
      void setData(const LID& blockLID,const unsigned int& cell,const Realf& value);
      #endif

    private:
      void exitInvalidLocalID(const LID& localID,const std::string& funcName) const;
      void resize();

      std::vector<Realf,aligned_allocator<Realf,WID3> > block_data;
      Realf null_block_data[WID3];
      LID currentCapacity;
      LID numberOfBlocks;
      std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > parameters;

#ifdef USE_CUDA
      LID dev_allocatedSize;
      Realf *dev_block_data;
      Real *dev_parameters;
      bool pinnedBlocks;
      bool pinnedParameters;
#endif

   };

   template<typename LID> inline
   VelocityBlockContainer<LID>::VelocityBlockContainer() {
      currentCapacity = 0;
      numberOfBlocks = 0;
#ifdef USE_CUDA
      dev_allocatedSize = 0;
      dev_block_data = new Realf();
      dev_parameters = new Real();
      pinnedBlocks = false;
      pinnedParameters = false;
      dev_needsUpdatingBlocks = true;
      dev_needsUpdatingParameters = true;
#endif
   }

#ifdef USE_CUDA
   template<typename LID> inline
   VelocityBlockContainer<LID>::~VelocityBlockContainer() {
      dev_Deallocate();
      delete[] dev_block_data;
      delete[] dev_parameters;
      dev_unpinBlocks();
      dev_unpinParameters();

      block_data.~vector();
      parameters.~vector();
   }
#endif

   template<typename LID> inline
   LID VelocityBlockContainer<LID>::capacity() const {
      return currentCapacity;
   }

   template<typename LID> inline
   size_t VelocityBlockContainer<LID>::capacityInBytes() const {
      return (block_data.capacity())*sizeof(Realf) + parameters.capacity()*sizeof(Real);
   }

   /** Clears VelocityBlockContainer data and deallocates all memory
    * reserved for velocity blocks.*/
   template<typename LID> inline
   void VelocityBlockContainer<LID>::clear() {
      std::vector<Realf,aligned_allocator<Realf,WID3> > dummy_data;
      std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > dummy_parameters;

      block_data.swap(dummy_data);
      parameters.swap(dummy_parameters);
#ifdef USE_CUDA
      dev_Deallocate();
      dev_unpinBlocks();
      dev_unpinParameters();
#endif
      currentCapacity = 0;
      numberOfBlocks = 0;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::copy(const LID& source,const LID& target) {
      #ifdef DEBUG_VBC
         bool ok = true;
         if (source >= numberOfBlocks) ok = false;
         if (source >= currentCapacity) ok = false;
         if (target >= numberOfBlocks) ok = false;
         if (target >= currentCapacity) ok = false;
         if (numberOfBlocks >= currentCapacity) ok = false;
         if (source != numberOfBlocks-1) ok = false;
         if (source*WID3+WID3-1 >= block_data.size()) ok = false;
         if (source*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::N_VELOCITY_BLOCK_PARAMS-1 >= parameters.size()) ok = false;
         if (target*BlockParams::N_VELOCITY_BLOCK_PARAMS+BlockParams::N_VELOCITY_BLOCK_PARAMS-1 >= parameters.size()) ok = false;
         if (parameters.size()/BlockParams::N_VELOCITY_BLOCK_PARAMS != block_data.size()/WID3) ok = false;
         if (ok == false) {
            std::stringstream ss;
            ss << "VBC ERROR: invalid source LID=" << source << " in copy, target=" << target << " #blocks=" << numberOfBlocks << " capacity=" << currentCapacity << std::endl;
            ss << "or sizes are wrong, data.size()=" << block_data.size() << " parameters.size()=" << parameters.size() << std::endl;
            std::cerr << ss.str();
            sleep(1);
            exit(1);
         }
      #endif

      for (unsigned int i=0; i<WID3; ++i) block_data[target*WID3+i] = block_data[source*WID3+i];
      for (int i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) {
         parameters[target*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = parameters[source*BlockParams::N_VELOCITY_BLOCK_PARAMS+i];
      }
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::exitInvalidLocalID(const LID& localID,const std::string& funcName) const {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);

      std::stringstream ss;
      ss << "Process " << rank << ' ';
      ss << "Invalid localID " << localID << " used in function '" << funcName << "' max allowed value is " << numberOfBlocks << std::endl;
      std::cerr << ss.str();
      sleep(1);
      exit(1);
   }

   template<typename LID> inline
   double VelocityBlockContainer<LID>::getBlockAllocationFactor() {
      return BLOCK_ALLOCATION_FACTOR;
   }

   template<typename LID> inline
   Realf* VelocityBlockContainer<LID>::getData() {
      return block_data.data();
   }

   template<typename LID> inline
   const Realf* VelocityBlockContainer<LID>::getData() const {
      return block_data.data();
   }

   template<typename LID> inline
   Realf* VelocityBlockContainer<LID>::getData(const LID& blockLID) {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"getData");
         if (blockLID >= block_data.size()/WID3) exitInvalidLocalID(blockLID,"const getData const");
      #endif
      return block_data.data() + blockLID*WID3;
   }

   template<typename LID> inline
   const Realf* VelocityBlockContainer<LID>::getData(const LID& blockLID) const {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"const getData const");
         if (blockLID >= block_data.size()/WID3) exitInvalidLocalID(blockLID,"const getData const");
      #endif
      return block_data.data() + blockLID*WID3;
   }

#ifdef USE_CUDA
   template<typename LID> inline
   Realf* VelocityBlockContainer<LID>::dev_getData() {
      return dev_block_data;
   }

   template<typename LID> inline
   const Realf* VelocityBlockContainer<LID>::dev_getData() const {
      return dev_block_data;
   }

   template<typename LID> inline
   Realf* VelocityBlockContainer<LID>::dev_getData(const LID& blockLID) {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"getData");
         if (blockLID >= block_data.size()/WID3) exitInvalidLocalID(blockLID,"const getData const");
      #endif
      return dev_block_data + blockLID*WID3;
   }

   template<typename LID> inline
   const Realf* VelocityBlockContainer<LID>::dev_getData(const LID& blockLID) const {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"const getData const");
         if (blockLID >= block_data.size()/WID3) exitInvalidLocalID(blockLID,"const getData const");
      #endif
      return dev_block_data + blockLID*WID3;
   }

   template<typename LID> inline
   Real* VelocityBlockContainer<LID>::dev_getParameters() {
      return dev_parameters;
   }

   template<typename LID> inline
   const Real* VelocityBlockContainer<LID>::dev_getParameters() const {
      return dev_parameters;
   }

   template<typename LID> inline
   Real* VelocityBlockContainer<LID>::dev_getParameters(const LID& blockLID) {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"getParameters");
         if (blockLID >= parameters.size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) exitInvalidLocalID(blockLID,"getParameters");
      #endif
      return dev_parameters + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   template<typename LID> inline
   const Real* VelocityBlockContainer<LID>::dev_getParameters(const LID& blockLID) const {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"const getParameters const");
         if (blockLID >= parameters.size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) exitInvalidLocalID(blockLID,"getParameters");
      #endif
      return dev_parameters + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_Deallocate() {
      if (dev_allocatedSize > 0) cudaDeallocateBlockData(&dev_block_data, &dev_parameters);
      dev_allocatedSize = 0;
      return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_Allocate(LID size) {
      if (dev_allocatedSize > 0) {
         if ( (dev_allocatedSize > CUDA_BLOCK_SAFECTY_FACTOR * size) &&
              (dev_allocatedSize > numberOfBlocks) ) {
            return; // Still have enough buffer
         }
         cudaDeallocateBlockData(&dev_block_data, &dev_parameters);
      }

      if (numberOfBlocks > CUDA_BLOCK_ALLOCATION_FACTOR  *size) {
         std::cerr<<"Error in "<<__FILE__<<" line "<<__LINE__<<": attempting to allocate less than numberOfBlocks"<<std::endl;
         abort();
      }
      if (numberOfBlocks > CUDA_BLOCK_SAFECTY_FACTOR * CUDA_BLOCK_ALLOCATION_FACTOR * size) {
         std::cerr<<"Warning in "<<__FILE__<<" line "<<__LINE__<<": attempting to allocate less than safety margins"<<std::endl;
      }
      cudaAllocateBlockData(&dev_block_data, &dev_parameters, CUDA_BLOCK_ALLOCATION_FACTOR*size);
      dev_allocatedSize = CUDA_BLOCK_ALLOCATION_FACTOR*size;
      return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_Allocate() {
      if (dev_allocatedSize > 0) {
         if ( dev_allocatedSize > CUDA_BLOCK_SAFECTY_FACTOR * numberOfBlocks) {
            return; // Still have enough buffer
         }
         cudaDeallocateBlockData(&dev_block_data, &dev_parameters);
      }
      cudaAllocateBlockData(&dev_block_data, &dev_parameters, CUDA_BLOCK_ALLOCATION_FACTOR*numberOfBlocks);
      dev_allocatedSize = CUDA_BLOCK_ALLOCATION_FACTOR*numberOfBlocks;
      return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_syncBlocksToHost() {
      if (numberOfBlocks==0) return;
      if (dev_allocatedSize == 0) {
         std::cerr<<"Error in "<<__FILE__<<" line "<<__LINE__<<": attempting to syncToHost without allocated GPU memory"<<std::endl;
         abort();
      }
      cuda_DtoH_BlockData(dev_block_data, block_data.data(), numberOfBlocks);
      return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_syncBlocksToDevice() {
      if (numberOfBlocks==0) return;
      if (!dev_needsUpdatingBlocks) return;
      if (dev_allocatedSize == 0) {
         std::cerr<<"Error in "<<__FILE__<<" line "<<__LINE__<<": attempting to syncToDevice without allocated GPU memory"<<std::endl;
         abort();
      }
      cuda_HtoD_BlockData(dev_block_data, block_data.data(), numberOfBlocks);
      dev_needsUpdatingBlocks = false;
      return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_syncParametersToHost() {
      if (numberOfBlocks==0) return;
      if (dev_allocatedSize == 0) {
         std::cerr<<"Error in "<<__FILE__<<" line "<<__LINE__<<": attempting to syncToHost without allocated GPU memory"<<std::endl;
         abort();
      }
      cuda_DtoH_BlockParameters(dev_parameters, parameters.data(), numberOfBlocks);
      return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_syncParametersToDevice() {
      if (numberOfBlocks==0) return;
      if (!dev_needsUpdatingParameters) return;
      if (dev_allocatedSize == 0) {
         std::cerr<<"Error in "<<__FILE__<<" line "<<__LINE__<<": attempting to syncToDevice without allocated GPU memory"<<std::endl;
         abort();
      }
      cuda_HtoD_BlockParameters(dev_parameters, parameters.data(), numberOfBlocks);
      dev_needsUpdatingParameters = false;
      return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_pinBlocks() {
      //cuda_register_BlockData(block_data.data(), numberOfBlocks);
      if (pinnedBlocks) {
         cuda_unregister_BlockData(block_data.data());
      }
      cuda_register_BlockData(block_data.data(), currentCapacity);
      pinnedBlocks = true;
   return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_unpinBlocks() {
      if (pinnedBlocks) {
         cuda_unregister_BlockData(block_data.data());
      }
      pinnedBlocks = false;
      return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_pinParameters() {
      //cuda_register_BlockParameters(parameters.data(), numberOfBlocks);
      if (pinnedParameters) {
         cuda_unregister_BlockParameters(parameters.data());
      }
      cuda_register_BlockParameters(parameters.data(), currentCapacity);
      pinnedParameters = true;
      return;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::dev_unpinParameters() {
      if (pinnedParameters) {
         cuda_unregister_BlockParameters(parameters.data());
      }
      pinnedParameters = false;
      return;
   }

#endif

   template<typename LID> inline
   Realf* VelocityBlockContainer<LID>::getNullData() {
       return null_block_data;
   }

   template<typename LID> inline
   Real* VelocityBlockContainer<LID>::getParameters() {
      return parameters.data();
   }

   template<typename LID> inline
   const Real* VelocityBlockContainer<LID>::getParameters() const {
      return parameters.data();
   }

   template<typename LID> inline
   Real* VelocityBlockContainer<LID>::getParameters(const LID& blockLID) {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"getParameters");
         if (blockLID >= parameters.size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) exitInvalidLocalID(blockLID,"getParameters");
      #endif
      return parameters.data() + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   template<typename LID> inline
   const Real* VelocityBlockContainer<LID>::getParameters(const LID& blockLID) const {
      #ifdef DEBUG_VBC
         if (blockLID >= numberOfBlocks) exitInvalidLocalID(blockLID,"const getParameters const");
         if (blockLID >= parameters.size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) exitInvalidLocalID(blockLID,"getParameters");
      #endif
      return parameters.data() + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::pop() {
      if (numberOfBlocks == 0) return;
      --numberOfBlocks;
   }

   template<typename LID> inline
   LID VelocityBlockContainer<LID>::push_back() {
      LID newIndex = numberOfBlocks;
      if (newIndex >= currentCapacity) resize();

      #ifdef DEBUG_VBC
      if (newIndex >= block_data.size()/WID3 || newIndex >= parameters.size()/BlockParams::N_VELOCITY_BLOCK_PARAMS) {
         std::stringstream ss;
         ss << "VBC ERROR in push_back, LID=" << newIndex << " for new block is out of bounds" << std::endl;
         ss << "\t data.size()=" << block_data.size()  << " parameters.size()=" << parameters.size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }
      #endif

      // Clear velocity block data to zero values
      for (size_t i=0; i<WID3; ++i) block_data[newIndex*WID3+i] = 0.0;
      for (size_t i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i)
         parameters[newIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = 0.0;

      ++numberOfBlocks;
      return newIndex;
   }

   template<typename LID> inline
   LID VelocityBlockContainer<LID>::push_back(const uint32_t& N_blocks) {
      const LID newIndex = numberOfBlocks;
      numberOfBlocks += N_blocks;
      resize();

      // Clear velocity block data to zero values
      for (size_t i=0; i<WID3*N_blocks; ++i) block_data[newIndex*WID3+i] = 0.0;
      for (size_t i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS*N_blocks; ++i)
	parameters[newIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = 0.0;

      return newIndex;
   }

   template<typename LID> inline
   bool VelocityBlockContainer<LID>::recapacitate(const LID& newCapacity) {
      if (newCapacity < numberOfBlocks) return false;
#ifdef USE_CUDA
      dev_unpinBlocks();
      dev_unpinParameters();
#endif
      {
         std::vector<Realf,aligned_allocator<Realf,WID3> > dummy_data(newCapacity*WID3);
         for (size_t i=0; i<numberOfBlocks*WID3; ++i) dummy_data[i] = block_data[i];
         dummy_data.swap(block_data);
      }
      {
         std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > dummy_parameters(newCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
         for (size_t i=0; i<numberOfBlocks*BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) dummy_parameters[i] = parameters[i];
         dummy_parameters.swap(parameters);
      }
      currentCapacity = newCapacity;
#ifdef USE_CUDA
      dev_pinBlocks();
      dev_pinParameters();
#endif
   return true;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::resize() {
      if ((numberOfBlocks+1) >= currentCapacity) {
         // Resize so that free space is block_allocation_chunk blocks,
         // and at least two in case of having zero blocks.
         // The order of velocity blocks is unaltered.
#ifdef USE_CUDA
         dev_unpinBlocks();
         dev_unpinParameters();
#endif
         currentCapacity = 2 + numberOfBlocks * BLOCK_ALLOCATION_FACTOR;
         block_data.resize(currentCapacity*WID3);
         parameters.resize(currentCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
#ifdef USE_CUDA
         dev_pinBlocks();
         dev_pinParameters();
#endif
      }
   }

   template<typename LID> inline
   bool VelocityBlockContainer<LID>::setSize(const LID& newSize) {
      numberOfBlocks = newSize;
      if (newSize > currentCapacity) resize();
      return true;
   }

   /** Return the number of existing velocity blocks.
    * @return Number of existing velocity blocks.*/
   template<typename LID> inline
   LID VelocityBlockContainer<LID>::size() const {
      return numberOfBlocks;
   }

   template<typename LID> inline
   size_t VelocityBlockContainer<LID>::sizeInBytes() const {
      return block_data.size()*sizeof(Realf) + parameters.size()*sizeof(Real);
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::swap(VelocityBlockContainer& vbc) {
#ifdef USE_CUDA
      dev_unpinBlocks();
      dev_unpinParameters();
#endif
      block_data.swap(vbc.block_data);
      parameters.swap(vbc.parameters);

      LID dummy = currentCapacity;
      currentCapacity = vbc.currentCapacity;
      vbc.currentCapacity = dummy;

      dummy = numberOfBlocks;
      numberOfBlocks = vbc.numberOfBlocks;
      vbc.numberOfBlocks = dummy;

#ifdef USE_CUDA
      uint32_t dummy2 = *dev_block_data;
      *dev_block_data = *vbc.dev_block_data;
      *vbc.dev_block_data = dummy2;

      dummy2 = *dev_parameters;
      *dev_parameters = *vbc.dev_parameters;
      *vbc.dev_parameters = dummy2;

      dev_pinBlocks();
      dev_pinParameters();
#endif
   }

   #ifdef DEBUG_VBC

   template<typename LID> inline
   const Realf& VelocityBlockContainer<LID>::getData(const LID& blockLID,const unsigned int& cell) const {
      bool ok = true;
      if (cell >= WID3) ok = false;
      if (blockLID >= numberOfBlocks) ok = false;
      if (blockLID*WID3+cell >= block_data.size()) ok = false;
      if (ok == false) {
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in getData, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " data.size()=" << block_data.size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }

      return block_data[blockLID*WID3+cell];
   }

   template<typename LID> inline
   const Real& VelocityBlockContainer<LID>::getParameters(const LID& blockLID,const unsigned int& cell) const {
      bool ok = true;
      if (cell >= BlockParams::N_VELOCITY_BLOCK_PARAMS) ok = false;
      if (blockLID >= numberOfBlocks) ok = false;
      if (blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+cell >= parameters.size()) ok = false;
      if (ok == false) {
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in getParameters, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " parameters.size()=" << parameters.size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }

      return parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS+cell];
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::setData(const LID& blockLID,const unsigned int& cell,const Realf& value) {
      bool ok = true;
      if (cell >= WID3) ok = false;
      if (blockLID >= numberOfBlocks) ok = false;
      if (blockLID*WID3+cell >= block_data.size()) ok = false;
      if (ok == false) {
         std::stringstream ss;
         ss << "VBC ERROR: out of bounds in setData, LID=" << blockLID << " cell=" << cell << " #blocks=" << numberOfBlocks << " data.size()=" << block_data.size() << std::endl;
         std::cerr << ss.str();
         sleep(1);
         exit(1);
      }

      block_data[blockLID*WID3+cell] = value;
   }

   #endif

} // namespace block_cont

#endif
