/*
 * This file is part of Vlasiator.
 * Copyright 2010-2025 Finnish Meteorological Institute and University of Helsinki
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

#ifndef GPU_BASE_H
#define GPU_BASE_H

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "arch_device_api.h"

#include <stdio.h>
#include <mutex>
#include "include/splitvector/splitvec.h"
#include "include/hashinator/hashinator.h"
#include "../definitions.h"
#include "../vlasovsolver/vec.h"
#include "../velocity_mesh_parameters.h"
#include <phiprof.hpp>

#ifndef THREADS_PER_MP
#define THREADS_PER_MP 2048
#endif
#ifndef REGISTERS_PER_MP
#define REGISTERS_PER_MP 65536
#endif

// Device properties
extern int gpuMultiProcessorCount;
extern int blocksPerMP;
extern int threadsPerMP;

// Magic multipliers used to make educated guesses for initial allocations
// and for managing dynamic increases in allocation sizes. Some of these are
// scaled based on WID value for better guesses,
static const uint VLASOV_BUFFER_MINBLOCKS = 32768/WID3;
static const uint VLASOV_BUFFER_MINCOLUMNS = 2000/WID;
static const uint INIT_VMESH_SIZE (32768/WID3);
static const uint INIT_MAP_SIZE (16 - WID);
static const double BLOCK_ALLOCATION_PADDING = 1.2;
static const double BLOCK_ALLOCATION_FACTOR = 1.1;

// Used in acceleration column construction. The flattened version of the
// probe cube must store (5) counters / offsets, see vlasovsolver/gpu_acc_map.cpp for details.
static const int GPU_PROBEFLAT_N = 5;

// buffers need to be larger for translation to allow proper parallelism
// GPUTODO: Get rid of this multiplier and consolidate buffer allocations.
// WARNING: Simply removing this factor led to diffs in Flowthrough_trans_periodic, indicating that
// there is somethign wrong with the evaluation of buffers! To be investigated.
static const int TRANSLATION_BUFFER_ALLOCATION_FACTOR = 5;

#define MAXCPUTHREADS 512 // hypothetical max size for some allocation arrays

void gpu_init_device();
void gpu_clear_device();
gpuStream_t gpu_getStream();
gpuStream_t gpu_getPriorityStream();
uint gpu_getThread();
uint gpu_getMaxThreads();
int gpu_getDevice();
uint gpu_getAllocationCount();
int gpu_reportMemory(const size_t local_cap=0, const size_t ghost_cap=0, const size_t local_size=0, const size_t ghost_size=0);

unsigned int nextPowerOfTwo(unsigned int n);

void gpu_vlasov_allocate(uint maxBlockCount);
void gpu_calculateProbeAllocation(uint maxBlockCount);
void gpu_vlasov_deallocate();
void gpu_vlasov_allocate_perthread(uint cpuThreadID, uint maxBlockCount);
uint gpu_vlasov_getSmallestAllocation();

void gpu_batch_allocate(uint nCells=0, uint maxNeighbours=0);

void gpu_acc_allocate(uint maxBlockCount);
void gpu_acc_allocate_perthread(uint cpuThreadID, uint firstAllocationCount, uint columnSetAllocationCount=0);
void gpu_acc_deallocate();

void gpu_trans_allocate(cuint nAllCells=0,
                        cuint largestVmesh=0,
                        cuint unionSetSize=0);
void gpu_trans_deallocate();

extern gpuStream_t gpuStreamList[];
extern gpuStream_t gpuPriorityStreamList[];

// Struct used by Vlasov Acceleration semi-Lagrangian solver
struct ColumnOffsets {
   split::SplitVector<uint> setColumnOffsets; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
   split::SplitVector<uint> setNumColumns; // how many columns in set of columns (length nColumnSets)

   split::SplitVector<uint> columnBlockOffsets; // indexes where columns start (in blocks, length totalColumns)
   split::SplitVector<uint> columnNumBlocks; // length of column (in blocks, length totalColumns)
   split::SplitVector<int> minBlockK,maxBlockK;
   split::SplitVector<int> kBegin;
   split::SplitVector<int> i,j;
   uint colSize = 0;
   uint colSetSize = 0;
   uint colCapacity = 0;
   uint colSetCapacity = 0;

   ColumnOffsets(uint nColumns=1, uint nColumnSets=1) {
      gpuStream_t stream = gpu_getStream();
      setColumnOffsets.resize(nColumnSets);
      setNumColumns.resize(nColumnSets);
      columnBlockOffsets.resize(nColumns);
      columnNumBlocks.resize(nColumns);
      minBlockK.resize(nColumns);
      maxBlockK.resize(nColumns);
      kBegin.resize(nColumns);
      i.resize(nColumns);
      j.resize(nColumns);
      // These vectors themselves are not in unified memory, just their content data
      setColumnOffsets.optimizeGPU(stream);
      setNumColumns.optimizeGPU(stream);
      columnBlockOffsets.optimizeGPU(stream);
      columnNumBlocks.optimizeGPU(stream);
      minBlockK.optimizeGPU(stream);
      maxBlockK.optimizeGPU(stream);
      kBegin.optimizeGPU(stream);
      i.optimizeGPU(stream);
      j.optimizeGPU(stream);
      // Cached values
      colSize = nColumns;
      colSetSize = nColumnSets;
      colCapacity = columnBlockOffsets.capacity(); // Uses this as an example
      colSetCapacity = setNumColumns.capacity(); // Uses this as an example
   }
   void prefetchDevice(gpuStream_t stream) {
      setColumnOffsets.optimizeGPU(stream);
      setNumColumns.optimizeGPU(stream);
      columnBlockOffsets.optimizeGPU(stream);
      columnNumBlocks.optimizeGPU(stream);
      minBlockK.optimizeGPU(stream);
      maxBlockK.optimizeGPU(stream);
      kBegin.optimizeGPU(stream);
      i.optimizeGPU(stream);
      j.optimizeGPU(stream);
   }
   __host__ size_t sizeCols() const {
      return colSize;
   }
   __host__ size_t capacityCols() const {
      return colCapacity;
   }
   __host__ size_t capacityColSets() const {
      return colSetCapacity;
   }
   __device__ size_t dev_sizeCols() const {
      return columnBlockOffsets.size(); // Uses this as an example
   }
   __device__ size_t dev_sizeColSets() const {
      return setNumColumns.size(); // Uses this as an example
   }
   __device__ size_t dev_capacityCols() const {
      return columnBlockOffsets.capacity(); // Uses this as an example
   }
   __device__ size_t dev_capacityColSets() const {
      return setNumColumns.capacity(); // Uses this as an example
   }
   size_t capacityInBytes() const {
      return colCapacity * (2*sizeof(uint)+5*sizeof(int))
         + colSetCapacity * (2*sizeof(uint))
         + 4 * sizeof(split::SplitVector<uint>)
         + 5 * sizeof(split::SplitVector<int>);
   }
   void setSizes(size_t nCols=0, size_t nColSets=0) {
      // Ensure capacities are handled with cached values
      setCapacities(nCols,nColSets);
      // Only then resize
      setColumnOffsets.resize(nColSets,true);
      setNumColumns.resize(nColSets,true);
      columnBlockOffsets.resize(nCols,true);
      columnNumBlocks.resize(nCols,true);
      minBlockK.resize(nCols,true);
      maxBlockK.resize(nCols,true);
      kBegin.resize(nCols,true);
      i.resize(nCols,true);
      j.resize(nCols,true);
      colSize = nCols;
      colSetSize = nColSets;
   }
   __device__ void device_setSizes(size_t nCols=0, size_t nColSets=0) {
      // Cannot recapacitate
      setColumnOffsets.device_resize(nColSets);
      setNumColumns.device_resize(nColSets);
      columnBlockOffsets.device_resize(nCols);
      columnNumBlocks.device_resize(nCols);
      minBlockK.device_resize(nCols);
      maxBlockK.device_resize(nCols);
      kBegin.device_resize(nCols);
      i.device_resize(nCols);
      j.device_resize(nCols);
      colSize = nCols;
      colSetSize = nColSets;
   }
   void setCapacities(size_t nCols=0, size_t nColSets=0) {
      // check cached capacities to prevent page faults if not necessary
      if (nCols > colCapacity) {
         // Recapacitate column vectors
         colCapacity = nCols * BLOCK_ALLOCATION_PADDING;
         columnBlockOffsets.reallocate(colCapacity);
         columnNumBlocks.reallocate(colCapacity);
         minBlockK.reallocate(colCapacity);
         maxBlockK.reallocate(colCapacity);
         kBegin.reallocate(colCapacity);
         i.reallocate(colCapacity);
         j.reallocate(colCapacity);
      }
      if (nColSets > colSetCapacity) {
         // Recapacitate columnSet vectors
         colSetCapacity = nColSets * BLOCK_ALLOCATION_PADDING;
         setColumnOffsets.reallocate(colSetCapacity);
         setNumColumns.reallocate(colSetCapacity);
      }
   }
};

/*
Usage options:
(a):
1. Create your pointer using one of the macros for creating pointers. Pass the "name" of your
pointer to that macro. A corresponding variable in the memory manager will be automatically generated.
2. Allocate memory to that pointer by passing the "name" to one of the allocation macros
3. Get the pointer by passing the "name" to one of the get pointer methods macros
(b):
1. Define an index variable for your pointer and pass that to one of the methods for crate pointers
2. Allocate memory to that pointer by passing the index to one of the allocation methods
3. Get the pointer by passing the index to one of the get pointer methods
*/
struct GPUMemoryManager {
   // Store pointers and their allocation sizes
   std::vector<void*> gpuMemoryPointers;
   std::vector<size_t> allocationSizes;
   std::vector<uint> pointerDevice;
   std::vector<void*> sessionPointers;
   std::vector<size_t> sessionPointerOffset;
   std::vector<uint> sessionPointerDevice;
   std::vector<size_t> sessionAllocationSizes;
   std::mutex memoryMutex;
   bool sessionOn = false;
   size_t dev_sessionSize = 0;
   size_t host_sessionSize = 0;
   size_t dev_sessionAllocationSize = 0;
   size_t host_sessionAllocationSize = 0;
   size_t dev_previousSessionSize = 0;
   size_t host_previousSessionSize = 0;
   size_t maxPointerIndex = 0;
   size_t maxSessionPointerIndex = 0;

   // Indices for session pointers
   size_t dev_sessionPointer = 0;
   size_t host_sessionPointer = 0;

   #define NO_POINTER_DEVICE 0
   #define DEVICE_POINTER 1
   #define HOST_POINTER 2

   // Useful for passing typenames with commas to other macros
   #define SINGLE_ARG(...) __VA_ARGS__

   /* 
   Definitions for global pointer indices
   Run the updateGpuMemoryPointerList.sh script to automatically update the list
   All "names" for pointersinitialized with CREATE_UNIQUE_POINTER, CREATE_SUBPOINTERS, SESSION_HOST_ALLOCATE
   , or SESSION_HOST_ALLOCATE are listed here.
   */
   #define DEFINITIONS_HERE
   size_t dev_allMaps, dev_allPencilsContainers, dev_allPencilsMeshes, dev_blockDataOrdered, dev_bulkVX, dev_bulkVY, dev_bulkVZ, dev_bValues, dev_cellIdxArray, dev_cellIdxKeys, dev_cellIdxStartCutoff, dev_columnOffsetData, dev_Ddt, dev_densityPostAdjust, dev_densityPreAdjust, dev_dfdt_mu, dev_dxdydz, dev_fcount, dev_fmu, dev_intersections, dev_lists_delete, dev_lists_to_replace, dev_lists_with_replace_new, dev_lists_with_replace_old, dev_mass, dev_massLoss, dev_max_dt, dev_minValues, dev_moments1, dev_moments2, dev_nAfter, dev_nBefore, dev_nBlocksToChange, dev_nColumns, dev_nColumnSets, dev_nu0Values, dev_nWithContent, dev_overflownElements, dev_pencilBlockData, dev_pencilBlocksCount, dev_potentialDdtValues, dev_probeCubeData, dev_remappedCellIdxArray, dev_resizeSuccess, dev_smallCellIdxArray, dev_sparsity, dev_VBC, dev_VBCs, dev_vbwcl_neigh, dev_vbwcl_vec, dev_velocityIdxArray, dev_vmeshes, gpu_block_indices_to_id, gpu_block_indices_to_probe, gpu_cell_indices_to_id, host_allMaps, host_allPencilsContainers, host_allPencilsMeshes, host_blockDataOrdered, host_bulkVX, host_bulkVY, host_bulkVZ, host_bValues, host_cellIdxStartCutoff, host_Ddt, host_dxdydz, host_intersections, host_lists_delete, host_lists_to_replace, host_lists_with_replace_new, host_lists_with_replace_old, host_mass, host_massLoss, host_max_dt, host_minValues, host_moments1, host_moments2, host_nAfter, host_nBefore, host_nBlocksToChange, host_nColumns, host_nColumnSets, host_nu0Values, host_nWithContent, host_overflownElements, host_remappedCellIdxArray, host_resizeSuccess, host_returnLID, host_returnReal, host_returnRealf, host_smallCellIdxArray, host_sparsity, host_VBC, host_VBCs, host_vbwcl_neigh, host_vbwcl_vec, host_vmeshes, member, my_test_pointer, returnLID, returnReal, returnRealf;
   #define DEFINITIONS_END
   #undef DEFINITIONS_HERE
   #undef DEFINITIONS_END

   // Macro for using a global "name" for a pointer without manually adding the pointer to a list
   // Other methods have similar macros for using pointers with "names"
   #define CREATE_UNIQUE_POINTER(object, member) object.createPointer(object.member)
   // Create a new pointer and replace the given index with the index corresponding to the pointer
   bool createPointer(size_t& pointerIndex) {
      if (pointerIndex != 0){
         return false;
      }

      std::lock_guard<std::mutex> lock(memoryMutex);

      // Leave first pointer empty
      if (maxPointerIndex == 0) {
         gpuMemoryPointers.push_back(nullptr);
         allocationSizes.push_back((size_t)(0));
         pointerDevice.push_back(NO_POINTER_DEVICE);
      }

      gpuMemoryPointers.push_back(nullptr);
      allocationSizes.push_back((size_t)(0));
      pointerDevice.push_back(NO_POINTER_DEVICE);

      maxPointerIndex++;
      pointerIndex = maxPointerIndex;
      
      return true;
   }

   #define CREATE_SUBPOINTERS(object, member, amount) object.createSubPointers(object.member, amount)
   // Create a given number of subpointers
   bool createSubPointers(size_t& firstPointerIndex, const size_t numberOfPointers) {
      if (firstPointerIndex != 0 || numberOfPointers == 0){
         return false;
      }

      std::lock_guard<std::mutex> lock(memoryMutex);

      // Leave first pointer empty
      if (maxPointerIndex == 0) {
         gpuMemoryPointers.push_back(nullptr);
         allocationSizes.push_back((size_t)(0));
         pointerDevice.push_back(NO_POINTER_DEVICE);
      }

      firstPointerIndex = maxPointerIndex + 1;

      for (size_t i = 0; i < numberOfPointers; i++) {
         maxPointerIndex++;
         gpuMemoryPointers.push_back(nullptr);
         allocationSizes.push_back((size_t)(0));
         pointerDevice.push_back(NO_POINTER_DEVICE);
      }
      
      return true;
   }

   // Start a session, where we have one big session pointer which can be split into multiple smaller pointers
   bool startSession(size_t dev_bytes, size_t host_bytes){
      // Ensure that the session pointers are at least as big as the largest session so far
      size_t host_requiredSessionSize = max(host_previousSessionSize, host_bytes);
      size_t dev_requiredSessionSize = max(dev_previousSessionSize, dev_bytes);

      createPointer(dev_sessionPointer);
      createPointer(host_sessionPointer);

      if(sessionOn){
         std::cerr << "Concurrent sessions not supported. Please end previous session before starting a new one.\n";
         return false;
      }
      sessionOn = true;

      // Reallocate session pointers if the required size increases
      if(dev_requiredSessionSize > dev_sessionAllocationSize){
         allocate(dev_sessionPointer, dev_requiredSessionSize);
         dev_sessionAllocationSize = dev_requiredSessionSize;
      }
      dev_sessionSize = 0;
      pointerDevice[dev_sessionPointer] = DEVICE_POINTER;

      if(host_requiredSessionSize > host_sessionAllocationSize){
         hostAllocate(host_sessionPointer, host_requiredSessionSize);
         host_sessionAllocationSize = host_requiredSessionSize;
      }
      host_sessionSize = 0;
      pointerDevice[host_sessionPointer] = HOST_POINTER;

      return true;
   }

   // Ending a session wipes the divisions of the session pointer
   bool endSession(){
      if(!sessionOn){
         std::cerr << "No session is currently on. Please start a session before ending it.\n";
         return false;
      }
      sessionOn = false;
      dev_previousSessionSize = dev_sessionSize;
      dev_sessionSize = 0;
      host_previousSessionSize = host_sessionSize;
      host_sessionSize = 0;

      // Free the pointers that did not fit into the session pointer
      freeSessionPointers();

      maxSessionPointerIndex = 0;

      return true;
   }

   #define ALLOCATE_GPU(object, member, bytes) object.allocate(object.member, bytes)
   // Allocate memory to a pointer by index
   bool allocate(const size_t& pointerIndex, size_t bytes) {
      if (pointerIndex == 0) {
         std::cerr << "Error: Pointer not found in 'gpuMemoryManager.allocate'.\n";
         return false;
      }

      std::lock_guard<std::mutex> lock(memoryMutex);

      if (allocationSizes[pointerIndex] >= bytes) {
         //No need to reallocate
         return false;
      }
      
      if (gpuMemoryPointers[pointerIndex] != nullptr) {
         CHK_ERR( gpuFree(gpuMemoryPointers[pointerIndex]) );
      }

      CHK_ERR( gpuMalloc(&gpuMemoryPointers[pointerIndex], bytes) );
      allocationSizes[pointerIndex] = bytes;
      pointerDevice[pointerIndex] = DEVICE_POINTER;

      return true;
   }

   #define ALLOCATE_WITH_BUFFER(object, member, bytes, buffer) object.allocateWithBuffer(object.member, bytes, buffer)
   // Allocate memory to a pointer by index with an extra buffer added
   bool allocateWithBuffer(const size_t& pointerIndex, size_t bytes, size_t buffer) {
      if (pointerIndex == 0) {
         std::cerr << "Error: Pointer not found in 'gpuMemoryManager.allocateWithBuffer'.\n";
         return false;
      }

      std::lock_guard<std::mutex> lock(memoryMutex);

      if (allocationSizes[pointerIndex] >= bytes) {
         //No need to reallocate
         return false;
      }
      
      if (gpuMemoryPointers[pointerIndex] != nullptr) {
         CHK_ERR( gpuFree(gpuMemoryPointers[pointerIndex]) );
      }

      CHK_ERR( gpuMalloc(&gpuMemoryPointers[pointerIndex], bytes*buffer) );
      allocationSizes[pointerIndex] = bytes*buffer;
      pointerDevice[pointerIndex] = DEVICE_POINTER;

      return true;
   }

   #define HOST_ALLOCATE_GPU(object, member, bytes) object.hostAllocate(object.member, bytes)
   // Allocate pinned host memory to a pointer by index
   bool hostAllocate(const size_t& pointerIndex, size_t bytes) {
      if (pointerIndex == 0) {
         std::cerr << "Error: Pointer not found in 'gpuMemoryManager.hostAllocate'.\n";
         return false;
      }

      std::lock_guard<std::mutex> lock(memoryMutex);

      if (allocationSizes[pointerIndex] >= bytes) {
         //No need to reallocate
         return false;
      }

      if (gpuMemoryPointers[pointerIndex] != nullptr) {
         CHK_ERR( gpuFreeHost(gpuMemoryPointers[pointerIndex]) );
      }

      CHK_ERR( gpuMallocHost(&gpuMemoryPointers[pointerIndex], bytes) );
      allocationSizes[pointerIndex] = bytes;
      pointerDevice[pointerIndex] = HOST_POINTER;

      return true;
   }

   #define HOST_ALLOCATE_WITH_BUFFER(object, member, bytes, buffer) object.hostAllocateWithBuffer(object.member, bytes, buffer)
   // Allocate pinned host memory to a pointer by index with an extra buffer
   bool hostAllocateWithBuffer(const size_t& pointerIndex, size_t bytes, size_t buffer) {
      if (pointerIndex == 0) {
         std::cerr << "Error: Pointer not found in 'gpuMemoryManager.hostAllocateWithBuffer'.\n";
         return false;
      }

      std::lock_guard<std::mutex> lock(memoryMutex);

      if (allocationSizes[pointerIndex] >= bytes) {
         //No need to reallocate
         return false;
      }

      if (gpuMemoryPointers[pointerIndex] != nullptr) {
         CHK_ERR( gpuFreeHost(gpuMemoryPointers[pointerIndex]) );
      }

      CHK_ERR( gpuMallocHost(&gpuMemoryPointers[pointerIndex], bytes*buffer) );
      allocationSizes[pointerIndex] = bytes*buffer;
      pointerDevice[pointerIndex] = HOST_POINTER;

      return true;
   }

   #define SUBPOINTER_ALLOCATE(object, member, index, bytes) object.subPointerAllocate(object.member, index, bytes)
   // Allocate memory to a sub pointer by base index and index
   bool subPointerAllocate(const size_t& firstPointerIndex, const uint index, size_t bytes) {

      size_t pointerIndex = firstPointerIndex + index;

      if (firstPointerIndex == 0 || pointerIndex > maxPointerIndex) {
         std::cerr << "Error: Pointer not found in 'gpuMemoryManager.subPointerAllocate': firstPointerIndex = " << firstPointerIndex << ", pointerIndex = " << pointerIndex << ", maxPointerIndex = " << maxPointerIndex << "\n";
         return false;
      }

      return allocate(pointerIndex, bytes);
   }

   // Allocate memory to a sub pointer by base index and index with an extra buffer
   bool subPointerAllocateWithBuffer(const size_t& firstPointerIndex, const uint index, size_t bytes, size_t buffer) {

      size_t pointerIndex = firstPointerIndex + index;

      if (firstPointerIndex == 0 || pointerIndex > maxPointerIndex) {
         std::cerr << "Error: Pointer not found in 'gpuMemoryManager.subPointerAllocateWithBuffer'.\n";
         return false;
      }

      return allocateWithBuffer(pointerIndex, bytes, buffer);
   }

   #define SUBPOINTER_HOST_ALLOCATE(object, member, index, bytes) object.subPointerHostAllocate(object.member, index, bytes)
   // Allocate host memory to a sub pointer by base index and index
   bool subPointerHostAllocate(const size_t& firstPointerIndex, const uint index, size_t bytes) {

      size_t pointerIndex = firstPointerIndex + index;

      if (firstPointerIndex == 0 || pointerIndex > maxPointerIndex) {
         std::cerr << "Error: Pointer not found in 'gpuMemoryManager.subPointerHostAllocate'.\n";
         return false;
      }

      return hostAllocate(pointerIndex, bytes);
   }

   // Calculate the offset required to align a pointer
   template<typename T>
   size_t alignOffset(void* base, size_t offset) {
      uintptr_t fullAddress = reinterpret_cast<uintptr_t>(base) + offset;
      size_t alignment = std::max(alignof(T),size_t(256)); // Align to at least 256 bits
      size_t alignedAddress = (fullAddress + alignment - 1) & ~(alignment - 1);
      return alignedAddress - reinterpret_cast<uintptr_t>(base);
   }

   #define SESSION_ALLOCATE(object, type, member, bytes) object.sessionAllocate<type>(object.member, bytes)
   // Create a new pointer by partitioning the session pointer
   template<typename T>
   bool sessionAllocate(size_t& pointerIndex, size_t bytes){
      if(!sessionOn){
         std::cerr << "No session is currently on. Please start a session before allocating to it.\n";
         return false;
      }

      void *sessionPointer = getPointer<void>(dev_sessionPointer);
      size_t offset = alignOffset<T>(sessionPointer, dev_sessionSize);

      std::lock_guard<std::mutex> lock(memoryMutex);

      // Leave first pointer empty
      if (maxSessionPointerIndex == 0) {
         sessionPointerOffset.push_back((size_t)(0));
         sessionPointers.push_back(nullptr);
         sessionPointerDevice.push_back(NO_POINTER_DEVICE);
         sessionAllocationSizes.push_back((size_t)(0));
      }

      maxSessionPointerIndex++;
      pointerIndex = maxSessionPointerIndex;

      int padding = offset - dev_sessionSize;
      dev_sessionSize += bytes + padding;
      
      sessionPointerOffset.push_back(offset);
      sessionPointers.push_back(nullptr);
      sessionPointerDevice.push_back(DEVICE_POINTER);
      sessionAllocationSizes.push_back((size_t)(0));

      // If the pointer does not fit into the session pointer, create a new pointer
      if (dev_sessionSize > dev_sessionAllocationSize){
         sessionPointerOffset[pointerIndex] = dev_sessionSize;
         sessionAllocationSizes[pointerIndex] = bytes;
         CHK_ERR( gpuMalloc(&sessionPointers[pointerIndex], bytes) );
      }

      return true;
   }

   #define SESSION_HOST_ALLOCATE(object, type, member, bytes) object.sessionHostAllocate<type>(object.member, bytes)
   // Create a new pointer by partitioning the host session pointer
   template<typename T>
   bool sessionHostAllocate(size_t& pointerIndex, size_t bytes){
      if(!sessionOn){
         std::cerr << "No session is currently on. Please start a session before allocating to it.\n";
         return false;
      }

      void *sessionPointer = getPointer<void>(host_sessionPointer);
      size_t offset = alignOffset<T>(sessionPointer, host_sessionSize);

      std::lock_guard<std::mutex> lock(memoryMutex);

      // Leave first pointer empty
      if (maxSessionPointerIndex == 0) {
         sessionPointerOffset.push_back((size_t)(0));
         sessionPointers.push_back(nullptr);
         sessionPointerDevice.push_back(NO_POINTER_DEVICE);
         sessionAllocationSizes.push_back((size_t)(0));
      }

      maxSessionPointerIndex++;
      pointerIndex = maxSessionPointerIndex;

      int padding = offset - host_sessionSize;
      host_sessionSize += bytes + padding;
      
      sessionPointerOffset.push_back(offset);
      sessionPointers.push_back(nullptr);
      sessionPointerDevice.push_back(HOST_POINTER);
      sessionAllocationSizes.push_back((size_t)(0));

      // If the pointer does not fit into the session pointer, create a new pointer
      if (host_sessionSize > host_sessionAllocationSize){
         sessionPointerOffset[pointerIndex] = host_sessionSize;
         sessionAllocationSizes[pointerIndex] = bytes;
         CHK_ERR( gpuMallocHost(&sessionPointers[pointerIndex], bytes) );
      }

      return true;
   }

   // Get allocated size for a pointer
   size_t getSize(const size_t& pointerIndex) const {
      if (pointerIndex != 0){
         return allocationSizes[pointerIndex];
      }
      return 0;
   }

   // Get the total amount of GPU memory allocated with the memory manager
   size_t totalGpuAllocation(){
      size_t total = 0;

      for (size_t pointerIndex = 0; pointerIndex < gpuMemoryPointers.size(); pointerIndex++) {
         if (gpuMemoryPointers[pointerIndex] != nullptr) {
            if (pointerDevice[pointerIndex] == DEVICE_POINTER){
               total += allocationSizes[pointerIndex];
            }
         }
      }

      for (size_t pointerIndex = 0; pointerIndex < sessionPointers.size(); pointerIndex++) {
         if (sessionPointers[pointerIndex] != nullptr) {
            if (sessionPointerDevice[pointerIndex] == DEVICE_POINTER){
               total += sessionAllocationSizes[pointerIndex];
            }
         }
      }

      return total;
   }

   // get the total amount of host memory allocated with the memory manager
   size_t totalCpuAllocation(){
      size_t total = 0;

      for (size_t pointerIndex = 0; pointerIndex < gpuMemoryPointers.size(); pointerIndex++) {
         if (gpuMemoryPointers[pointerIndex] != nullptr) {
            if (pointerDevice[pointerIndex] == HOST_POINTER){
               total += allocationSizes[pointerIndex];
            }
         }
      }

      for (size_t pointerIndex = 0; pointerIndex < sessionPointers.size(); pointerIndex++) {
         if (sessionPointers[pointerIndex] != nullptr) {
            if (sessionPointerDevice[pointerIndex] == HOST_POINTER){
               total += sessionAllocationSizes[pointerIndex];
            }
         }
      }

      return total;
   }

   // Free all allocated pointers from session, besides the sessionpointer itself
   void freeSessionPointers() {
      for (size_t pointerIndex = 0; pointerIndex < sessionPointers.size(); pointerIndex++) {
         if (sessionPointers[pointerIndex] != nullptr) {
            if (sessionPointerDevice[pointerIndex] == DEVICE_POINTER){
               CHK_ERR( gpuFree(sessionPointers[pointerIndex]) );
            }else if (sessionPointerDevice[pointerIndex] == HOST_POINTER){
               CHK_ERR( gpuFreeHost(sessionPointers[pointerIndex]) );
            }
         }
      }

      sessionPointers.clear();
      sessionPointerOffset.clear();
      sessionPointerDevice.clear();
      sessionAllocationSizes.clear();
   }

   // Free all allocated GPU memory
   void freeAll() {
      for (size_t pointerIndex = 0; pointerIndex < gpuMemoryPointers.size(); pointerIndex++) {
         if (gpuMemoryPointers[pointerIndex] != nullptr) {
            if (allocationSizes[pointerIndex] > 0){
               if (pointerDevice[pointerIndex] == DEVICE_POINTER){
                  CHK_ERR( gpuFree(gpuMemoryPointers[pointerIndex]) );
               }else if (pointerDevice[pointerIndex] == HOST_POINTER){
                  CHK_ERR( gpuFreeHost(gpuMemoryPointers[pointerIndex]) );
               }
            }
         }
      }
      freeSessionPointers();

      gpuMemoryPointers.clear();
      allocationSizes.clear();
      pointerDevice.clear();
      sessionPointerOffset.clear();
      sessionOn = false;
      dev_sessionSize = 0;
      host_sessionSize = 0;
      dev_sessionAllocationSize = 0;
      host_sessionAllocationSize = 0;
      dev_previousSessionSize = 0;
      host_previousSessionSize = 0;
      maxPointerIndex = 0;
      maxSessionPointerIndex = 0;
      dev_sessionPointer = 0;
      host_sessionPointer = 0;
   }

   #define GET_POINTER(object, type, member) object.getPointer<type>(object.member)
   // Get typed pointer
   template <typename T>
   T* getPointer(const size_t& pointerIndex) const {
      if (pointerIndex == 0){
         throw std::runtime_error("Unknown pointer name at gpuMemoryManager.getPointer!\n");
      }
      return static_cast<T*>(gpuMemoryPointers[pointerIndex]);
   }

   #define GET_SUBPOINTER(object, type, member, index) object.getSubPointer<type>(object.member, index)
   // Get typed subpointer with an index
   template <typename T>
   T* getSubPointer(const size_t& firstPointerIndex, const uint index) const {

      size_t pointerIndex = firstPointerIndex + index;

      if (firstPointerIndex == 0 || pointerIndex > maxPointerIndex) {
         throw std::runtime_error("Unknown pointer name at gpuMemoryManager.getSubPointer!\n");
      }

      return getPointer<T>(pointerIndex);
   }

   #define GET_SESSION_POINTER(object, type, member) object.getSessionPointer<type>(object.member)
   // Get pointer in a session
   template <typename T>
   T* getSessionPointer(const size_t& pointerIndex) const {
      if (pointerIndex == 0 || pointerIndex > maxSessionPointerIndex){
         throw std::runtime_error("Unknown pointer name at gpuMemoryManager.getSessionPointer!\n");
      }

      // The pointer is usually contained in the session pointer and distinguished by an offset
      char *sessionPointer = static_cast<char*>(gpuMemoryPointers[dev_sessionPointer]);
      size_t offset = sessionPointerOffset[pointerIndex];

      // If the pointer did not fit inside the session pointer, retrieve it from the separate map
      if (offset > dev_sessionAllocationSize){
         return static_cast<T*>(sessionPointers[pointerIndex]);
      }
      
      return reinterpret_cast<T*>(sessionPointer + offset);
   }

   #define GET_SESSION_HOST_POINTER(object, type, member) object.getSessionHostPointer<type>(object.member)
   // Get host pointer in a session
   template <typename T>
   T* getSessionHostPointer(const size_t& pointerIndex) const {
      if (pointerIndex == 0 || pointerIndex > maxSessionPointerIndex){
         throw std::runtime_error("Unknown pointer name at gpuMemoryManager.getSessionHostPointer!\n");
      }

      // The pointer is usually contained in the session pointer and distinguished by an offset
      char *sessionPointer = static_cast<char*>(gpuMemoryPointers[host_sessionPointer]);
      size_t offset = sessionPointerOffset[pointerIndex];

      // If the pointer did not fit inside the session pointer, retrieve it from the separate map
      if (offset > host_sessionAllocationSize){
         return static_cast<T*>(sessionPointers[pointerIndex]);
      }
      
      return reinterpret_cast<T*>(sessionPointer + offset);
   }

   // Update a new pointer to the memory manager
   void updatePointer(const size_t& pointerIndex, void* newPtr) {
      std::lock_guard<std::mutex> lock(memoryMutex);
      gpuMemoryPointers[pointerIndex] = newPtr;
   }

   #define SET_SUBPOINTER(object, type, member, index, subPointerIndex) object.setSubPointer<type>(object.member, index, subPointerIndex)
   // Set the value of the base pointer at a given index to be the corresponding sub pointer
   template <typename T>
   void setSubPointer(const size_t basePointerIndex, const size_t index, const size_t subPointerIndex){
      if (basePointerIndex == 0 || subPointerIndex == 0) {
         throw std::runtime_error("Error: Pointer not found in 'gpuMemoryManager.setSubPointer'.");
      }

      T** basePointer = static_cast<T**>(gpuMemoryPointers[basePointerIndex]);
      T* subPointer  = static_cast<T*>(gpuMemoryPointers[subPointerIndex]);
      basePointer[index] = subPointer;
   }
};

extern GPUMemoryManager gpuMemoryManager;

extern ColumnOffsets *host_columnOffsetData;
extern uint gpu_largest_columnCount;
extern size_t gpu_probeFullSize, gpu_probeFlattenedSize, gpu_probeStride;

// Hash map and splitvectors buffers used in block adjustment are declared in block_adjust_gpu.hpp
// Vector and set for use in translation are declared in vlasovsolver/gpu_trans_map_amr.hpp

// Counters used in allocations
extern std::vector<uint> gpu_vlasov_allocatedSize;
extern uint gpu_acc_allocatedColumns;
extern uint gpu_acc_foundColumnsCount;

#endif
