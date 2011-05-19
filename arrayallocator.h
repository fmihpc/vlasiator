#ifndef ARRAYALLOCATOR_H
#define ARRAYALLOCATOR_H

#include <map>
#include "definitions.h"

/** Class DeviceGrid manages allocation and deallocation of arrays on a 
 * graphics card device (GPU). DeviceGrid creates a single array on the 
 * device. A user can request (and free) a continuous piece of memory from the 
 * DeviceGrid. 
 */
class ArrayAllocator {
 public:

   ArrayAllocator();
   ~ArrayAllocator();
   
   bool finalize();
   template<typename T> T* getArray() const;
   template<typename T> T* getArray(cuint& offset) const;
   uint getOffset(cuint& cellIndex) const;
   bool initialize(cuint& maxSize,cuint& byteSize);
   
   /** Query if the ArrayAllocator is initialized. 
    * @return If true, all memory was allocated successfully and ArrayAllocator is ready for use.
    * If false, an error has occurred and the simulation should be aborted.
    */
   bool isInitialized() const {return initialized;}
   
   void releaseArray(cuint& cellIndex);
   uint reserveArray(cuint& cellIndex,cuint& size);
   
 private:
   /** Struct Allocation is used to keep track of memory allocations 
    * on the device. It contains all necessary information of a single 
    * memory allocation.*/
   struct Allocation {
      uint size;                           /**< Number of allocated elements.*/
      uint offset;                         /**< Offset into DeviceGrid::array, indicates the starting position of this allocation.*/
      int references;                      /**< Number of references to the allocation.*/
   };

   std::map<uint,Allocation> allocations;  /**< Container for all current array allocations.*/
   char* array;                            /**< Pointer to array in GPU memory.*/
   uint byteSize;                          /**< Byte size of an array element.*/
   std::map<uint,uint> freeSlots;          /**< List of empty slots in array (offset + number of elements).*/
   bool initialized;                       /**< If true, DeviceGrid was initialized successfully.*/
   uint maxSize;                           /**< Maximum number of array elements that can be allocated on the device.*/
   uint nextFree;                          /**< Offset of next unallocated element.*/
   
   bool allocateDeviceArray(cuint& byteSize);
   bool freeDeviceArray();
};

/** Get a pointer to the device array that the DeviceGrid is using.
 * @return Pointer to an array on device.
 */
template<typename T> T* ArrayAllocator::getArray() const {
   T* ptr = reinterpret_cast<T*>(array);
   return ptr;
}

/** Get a pointer to a memory allocation.
 * @param offset Offset of the memory allocation, returned by member 
 * functions DeviceGrid::getOffset and DeviceGrid::reserveArray.
 * @return Pointer to an array on device.
 */
template<typename T> T* ArrayAllocator::getArray(cuint& offset) const {
   T* ptr = reinterpret_cast<T*>(array);
   return ptr + offset;
}

#endif
