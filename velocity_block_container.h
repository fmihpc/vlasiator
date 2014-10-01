/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2014 Finnish Meteorological Institute
 */

#ifndef VELOCITY_BLOCK_CONTAINER_H
#define VELOCITY_BLOCK_CONTAINER_H

#include <vector>

#include "common.h"

namespace vmesh {

   static const double BLOCK_ALLOCATION_FACTOR = 1.1;

   template<typename LID>
   class VelocityBlockContainer {
    public:

      VelocityBlockContainer();
      LID capacity() const;
      size_t capacityInBytes() const;
      void clear();
      void copy(const LID& source,const LID& target);
      Realf* getData();
      const Realf* getData() const;
      Realf* getData(const LID& blockLID);
      const Realf* getData(const LID& blockLID) const;
      Realf* getFx();
      const Realf* getFx() const;
      Realf* getFx(const LID& blockLID);
      const Realf* getFx(const LID& blockLID) const;
      Real* getParameters();
      const Real* getParameters() const;
      Real* getParameters(const LID& blockLID);      
      const Real* getParameters(const LID& blockLID) const;
      void pop();
      LID push_back();
      bool recapacitate(const LID& capacity);
      bool setSize(const LID& newSize);
      LID size() const;
      size_t sizeInBytes() const;

    private:
      void resize();
      
      std::vector<Realf,aligned_allocator<Realf,WID3> > block_data;
      std::vector<Realf,aligned_allocator<Realf,WID3> > block_fx;
      LID currentCapacity;
      LID numberOfBlocks;
      std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > parameters;
   };
   
   template<typename LID> inline
   VelocityBlockContainer<LID>::VelocityBlockContainer() {
      currentCapacity = 0;
      numberOfBlocks = 0;
   }
   
   template<typename LID> inline
   LID VelocityBlockContainer<LID>::capacity() const {
      return currentCapacity;
   }
   
   template<typename LID> inline
   size_t VelocityBlockContainer<LID>::capacityInBytes() const {
      return (block_data.capacity()+block_fx.capacity())*sizeof(Realf) + parameters.capacity()*sizeof(Real);
   }

   /** Clears VelocityBlockContainer data and deallocates all memory 
    * reserved for velocity blocks.*/
   template<typename LID> inline
   void VelocityBlockContainer<LID>::clear() {
      std::vector<Realf,aligned_allocator<Realf,WID3> > dummy_data;
      std::vector<Realf,aligned_allocator<Realf,WID3> > dummy_fx;
      std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > dummy_parameters;
      
      block_data.swap(dummy_data);
      block_fx.swap(dummy_fx);
      parameters.swap(dummy_parameters);
      
      currentCapacity = 0;
      numberOfBlocks = 0;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::copy(const LID& source,const LID& target) {      
      for (int i=0; i<WID3; ++i) block_data[target*WID3+i] = block_data[source*WID3+i];
      for (int i=0; i<WID3; ++i) block_fx[target*WID3+i]   = block_fx[source*WID3+i];
      for (int i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) {
	 parameters[target*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = 
	   parameters[source*BlockParams::N_VELOCITY_BLOCK_PARAMS+i];
      }
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
      return block_data.data() + blockLID*WID3;
   }
   
   template<typename LID> inline
   const Realf* VelocityBlockContainer<LID>::getData(const LID& blockLID) const {
      return block_data.data() + blockLID*WID3;
   }

   template<typename LID> inline
   Realf* VelocityBlockContainer<LID>::getFx() {
      return block_fx.data();
   }
   
   template<typename LID> inline
   const Realf* VelocityBlockContainer<LID>::getFx() const {
      return block_fx.data();
   }
   
   template<typename LID> inline
   Realf* VelocityBlockContainer<LID>::getFx(const LID& blockLID) {
      return block_fx.data() + blockLID*WID3;
   }
   
   template<typename LID> inline
   const Realf* VelocityBlockContainer<LID>::getFx(const LID& blockLID) const {
      return block_fx.data() + blockLID*WID3;
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
      return parameters.data() + blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
   }
   
   template<typename LID> inline
   const Real* VelocityBlockContainer<LID>::getParameters(const LID& blockLID) const {
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
      
      // Clear velocity block data to zero values
      for (size_t i=0; i<WID3; ++i) block_data[newIndex*WID3+i] = 0.0;
      for (size_t i=0; i<WID3; ++i) block_fx[newIndex*WID3+i] = 0.0;
      for (size_t i=0; i<BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) 
	parameters[newIndex*BlockParams::N_VELOCITY_BLOCK_PARAMS+i] = 0.0;

      ++numberOfBlocks;
      return newIndex;
   }

   template<typename LID> inline
   bool VelocityBlockContainer<LID>::recapacitate(const LID& newCapacity) {
      if (newCapacity < numberOfBlocks) return false;
	{
	   std::vector<Realf,aligned_allocator<Realf,WID3> > dummy_data(newCapacity*WID3);
	   for (size_t i=0; i<numberOfBlocks*WID3; ++i) dummy_data[i] = block_data[i];
	   dummy_data.swap(block_data);
	}
	{
	   std::vector<Realf,aligned_allocator<Realf,WID3> > dummy_fx(newCapacity*WID3);
	   for (size_t i=0; i<numberOfBlocks*WID3; ++i) dummy_fx[i] = block_fx[i];
	   dummy_fx.swap(block_fx);
	}
	{
	   std::vector<Real,aligned_allocator<Real,BlockParams::N_VELOCITY_BLOCK_PARAMS> > dummy_parameters(newCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
	   for (size_t i=0; i<numberOfBlocks*BlockParams::N_VELOCITY_BLOCK_PARAMS; ++i) dummy_parameters[i] = parameters[i];
	   dummy_parameters.swap(parameters);
	}
      currentCapacity = newCapacity;
      return true;
   }

   template<typename LID> inline
   void VelocityBlockContainer<LID>::resize() {
      if ((numberOfBlocks+1) >= currentCapacity) {
	 // Resize so that free space is block_allocation_chunk blocks, 
	 // and at least two in case of having zero blocks.
	 // The order of velocity blocks is unaltered.
	 currentCapacity = 2 + numberOfBlocks * BLOCK_ALLOCATION_FACTOR;
	 block_data.resize(currentCapacity*WID3);
	 block_fx.resize(currentCapacity*WID3);
	 parameters.resize(currentCapacity*BlockParams::N_VELOCITY_BLOCK_PARAMS);
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
      return block_data.size()*sizeof(Realf) + block_fx.size()*sizeof(Realf) + parameters.size()*sizeof(Real);
   }

} // namespace block_cont

#endif
