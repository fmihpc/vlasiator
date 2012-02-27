/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

#include "mpilogger.h"
#include "arrayallocator.h"

extern MPILogger mpilogger;

using namespace std;

/** Constructor for class ArrayAllocator. The contructor is empty - one must 
 * call ArrayAllocator::initialize member function to initialize the 
 * ArrayAllocator.
 */
ArrayAllocator::ArrayAllocator() {
   nextFree = 0;
   initialized = false;
   array = NULL;
}

/** Destructor for class ArrayAllocator. Deallocates memory that 
 * has been allocated with ArrayAllocator::initialize.
 */
ArrayAllocator::~ArrayAllocator() {
   freeDeviceArray();
}

bool ArrayAllocator::finalize() {
   freeDeviceArray();
   nextFree = 0;
   array = NULL;
   initialized = false;
   allocations.clear();
   freeSlots.clear();
   return true;
}

/** Initialize ArrayAllocator. 
 * @param maxSize Number of array elements the ArrayAllocator should contain.
 * @param byteSize Byte size of an array element.
 * @return If true, DeviceGrid initialized successfully.
 */
bool ArrayAllocator::initialize(cuint& maxSize,cuint& byteSize) {
   initialized = allocateDeviceArray(maxSize*byteSize);
   if (initialized == false) {
      array = NULL;
   } else {
      this->maxSize = maxSize;
      this->byteSize = byteSize;
   }
   
   return initialized;
}

/** Allocate an array on the device with the given size.
 * @byteSize Byte size of the array to allocate.
 * @return If true, an array was allocated successfully.
 */
bool ArrayAllocator::allocateDeviceArray(cuint& byteSize) {
   bool success = true;
   array = new char[byteSize];
   return success;
}

/** Deallocate the device array that has been previously allocated 
 * with a call to DeviceGrid::initialize.
 * @return If true, the array was deallocated successfully.
 */
bool ArrayAllocator::freeDeviceArray() {
   bool success = true;
   delete array;
   array = NULL;
   return success;
}

/** Get an offset into DeviceGrid for a spatial cell which has reserved memory. 
 * @param cellIndex Global ID of the spatial cell.
 * @return Offset into DeviceGrid corresponding to the cell's memory allocation, or 
 * numeric_limits<uint>::max() if the cell does not have a memory allocation.
 */
uint ArrayAllocator::getOffset(cuint& cellIndex) const {
   map<uint,ArrayAllocator::Allocation>::const_iterator it = allocations.find(cellIndex);
   if (it == allocations.end()) return numeric_limits<uint>::max();
   return it->second.offset;
}

/** Release a memory allocation of the given spatial cell.
 * @param cellIndex Global ID of a spatial cell that has reserved memory from DeviceGrid.
 */
void ArrayAllocator::releaseArray(cuint& cellIndex) {
   // Check that an allocation exists:
   map<uint,ArrayAllocator::Allocation>::iterator it = allocations.find(cellIndex);
   if (it == allocations.end()) return;
   
   // Reduce reference count to released array. If the last reference was removed, add 
   // an empty hole to freeAvgs:
   --(it->second.references);
   if (it->second.references > 0) return;   
   pair<map<uint,uint>::iterator,bool> position = freeSlots.insert(make_pair(it->second.offset,it->second.size));
   
   // Combine adjacent holes, if possible. Note that this includes 
   // testing for a hole before and after this allocation, and 
   // combining those holes with the new one is the memory holes are continuous.
   bool removeNext = false;
   bool removeThis = false;
   --position.first;
   map<uint,uint>::iterator prev = position.first; // Memory hole before the new hole
   ++position.first;
   ++position.first;
   map<uint,uint>::iterator next = position.first; // Memory hole after the new hole
   --position.first;
   
   // Attempt to combine the new hole with the next hole. This will not succeed if 
   // the iterator next points to the beginning of freeSlots.
   if (next != freeSlots.end() && next != position.first) {
     if (position.first->first + position.first->second == next->first) {
	position.first->second += next->second;
	removeNext = true;
     }
   } else {

   }
   
   // Attempt to combine the new hole with the previous hole. This will not 
   // succeed if the iterator prev points to the beginning of freeSlots.
   if (prev != position.first) {
     if (prev->first + prev->second == position.first->first) {
	prev->second += position.first->second;
	removeThis = true;
     }
   } else {

   }
   
   if (removeNext == true) freeSlots.erase(next);
   if (removeThis == true) freeSlots.erase(position.first);
   
   // Erase the allocation for the given array:
   allocations.erase(it);
}

/** Allocate a piece of continuous memory for the given spatial cell. 
 * A spatial cell can only have a single allocation at any given time.
 * @param cellIndex Global ID of the spatial cell.
 * @param size Number of array elements to allocate.
 * @return Offset into DeviceArray corresponding to the reserved memory. If 
 * numeric_limits<uint>::max() is returned, DeviceGrid was unable to allocate 
 * the requested memory.
 */
uint ArrayAllocator::reserveArray(cuint& cellIndex,cuint& size) {
   // Check that memory has not already been allocated for the given cell:
   map<uint,Allocation>::iterator tmp = allocations.find(cellIndex);
   if (tmp != allocations.end()) {
      ++(tmp->second.references);
      return tmp->second.offset;
   }
   
   uint position = numeric_limits<uint>::max();
   
   // Allocate requested array from one of the empty slots in array avgs, if possible:
   if (freeSlots.size() > 0) {
      for (map<uint,uint>::iterator it = freeSlots.begin(); it != freeSlots.end(); ++it) {
	 if (it->second >= size) {
	    cuint difference = it->second - size;
	    position = it->first;
	    
	    if (difference > 0) 
	      freeSlots.insert(make_pair(position+size,difference));

	    freeSlots.erase(it);
	    break;
	 }
      }
   } 
   // If array could not be allocated from one of the holes, attempt to allocate it 
   // from the end of array avgs:
   if (position == numeric_limits<uint>::max()) {
      if (nextFree + size > maxSize) position = numeric_limits<uint>::max();
      else {
	 position = nextFree;
	 nextFree += size;
      }
   }
   
   // Check for out-of-memory:
   if (position == numeric_limits<uint>::max()) {return position;}
   
   // Create a new allocation:
   ArrayAllocator::Allocation alloc;
   alloc.size = size;
   alloc.offset = position;
   alloc.references = 1;
   allocations.insert(make_pair(cellIndex,alloc));
   return position;
}




