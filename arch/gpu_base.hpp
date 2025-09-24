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
void gpu_vlasov_deallocate();
void gpu_vlasov_allocate_perthread(uint cpuThreadID, uint maxBlockCount);
void gpu_vlasov_deallocate_perthread(uint cpuThreadID);
uint gpu_vlasov_getSmallestAllocation();

void gpu_batch_allocate(uint nCells=0, uint maxNeighbours=0);
void gpu_batch_deallocate(bool first=true, bool second=true);

void gpu_acc_allocate(uint maxBlockCount);
void gpu_acc_allocate_perthread(uint cpuThreadID, uint firstAllocationCount, uint columnSetAllocationCount=0);
void gpu_acc_deallocate();

void gpu_trans_allocate(cuint nAllCells=0,
                        cuint sumOfLengths=0,
                        cuint largestVmesh=0,
                        cuint unionSetSize=0,
                        cuint transGpuBlocks=0,
                        cuint nPencils=0);
void gpu_trans_deallocate();

void gpu_pitch_angle_diffusion_allocate(size_t numberOfLocalCells, int nbins_v, int nbins_mu, int blocksPerSpatialCell, int totalNumberOfVelocityBlocks);
void gpu_pitch_angle_diffusion_deallocate();

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

struct GPUMemoryManager {
   // Store pointers and their allocation sizes
   std::unordered_map<std::string, void*> gpuMemoryPointers;
   std::unordered_map<std::string, size_t> allocationSizes;
   std::unordered_map<std::string, int> nameCounters;
   std::unordered_map<std::string, std::string> pointerDevice;
   std::mutex memoryMutex;

   // Create a new pointer with a base name, ensure unique name
   bool createPointer(const std::string& baseName, std::string &uniqueName) {
      if (uniqueName != "null"){
         return false;
      }

      std::lock_guard<std::mutex> lock(memoryMutex);

      if (gpuMemoryPointers.count(baseName)) {
         int& counter = nameCounters[baseName];
         uniqueName = baseName + "_" + std::to_string(++counter);
      } else {
         nameCounters[baseName] = 0;
         uniqueName = baseName;
      }

      gpuMemoryPointers[uniqueName] = nullptr;
      allocationSizes[uniqueName] = (size_t)(0);
      pointerDevice[uniqueName] = "None";
      return true;
   }

   // Allocate memory to a pointer by name
   bool allocate(const std::string& name, size_t bytes) {
      //TODO: only allocate if needs to increase in size
      std::lock_guard<std::mutex> lock(memoryMutex);
      if (gpuMemoryPointers.count(name) == 0) {
         std::cerr << "Error: Pointer name '" << name << "' not found.\n";
         return false;
      }
      
      if (gpuMemoryPointers[name] != nullptr) {
         CHK_ERR( gpuFree(gpuMemoryPointers[name]) );
      }

      CHK_ERR( gpuMalloc(&gpuMemoryPointers[name], bytes) );
      allocationSizes[name] = bytes;
      pointerDevice[name] = "dev";
      return true;
   }

   // Allocate pinned host memory to a pointer by name
   bool hostAllocate(const std::string& name, size_t bytes) {
      //TODO: only allocate if needs to increase in size
      std::lock_guard<std::mutex> lock(memoryMutex);
      if (gpuMemoryPointers.count(name) == 0) {
         std::cerr << "Error: Pointer name '" << name << "' not found.\n";
         return false;
      }

      if (gpuMemoryPointers[name] != nullptr) {
         CHK_ERR( gpuFreeHost(gpuMemoryPointers[name]) );
      }

      CHK_ERR( gpuMallocHost(&gpuMemoryPointers[name], bytes) );
      allocationSizes[name] = bytes;
      pointerDevice[name] = "host";
      return true;
   }

   // Get allocated size for a pointer
   size_t getSize(const std::string& name) const {
      if (allocationSizes.count(name)) return allocationSizes.at(name);
      return 0;
   }

   // Free all allocated GPU memory
   void freeAll() {
      for (auto& pair : gpuMemoryPointers) {
         if (pair.second != nullptr) {
            std::string name = pair.first;
            if (allocationSizes[name] > 0){
               if (pointerDevice[name] == "dev"){
                  CHK_ERR( gpuFree(pair.second) );
               }else if (pointerDevice[name] == "host"){
                  CHK_ERR( gpuFreeHost(pair.second) );
               }
            }
            pair.second = nullptr;
         }
      }

      gpuMemoryPointers.clear();
      allocationSizes.clear();
      nameCounters.clear();
      pointerDevice.clear();
   }

   // Get typed pointer
   template <typename T>
   T* getPointer(const std::string& name) const {
      if (!gpuMemoryPointers.count(name)) throw std::runtime_error("Unknown pointer name");
      return static_cast<T*>(gpuMemoryPointers.at(name));
   }
};

extern GPUMemoryManager gpuMemoryManager;

// Device data variables, to be allocated in good time. Made into an array so that each thread has their own pointer.
extern Realf **host_blockDataOrdered;
extern Realf **dev_blockDataOrdered;
extern uint *gpu_cell_indices_to_id;
extern uint *gpu_block_indices_to_id;
extern uint *gpu_block_indices_to_probe;

extern Realf** dev_pencilBlockData;
extern uint* dev_pencilBlocksCount;

extern Real *returnReal[];
extern Realf *returnRealf[];
extern vmesh::LocalID *returnLID[];
extern Real *host_returnReal[];
extern Realf *host_returnRealf[];
extern vmesh::LocalID *host_returnLID[];

extern ColumnOffsets *host_columnOffsetData;
extern ColumnOffsets *dev_columnOffsetData;
extern uint gpu_largest_columnCount;
extern size_t gpu_probeFullSize, gpu_probeFlattenedSize;

// Hash map and splitvectors buffers used in block adjustment are declared in block_adjust_gpu.hpp
// Vector and set for use in translation are declared in vlasovsolver/gpu_trans_map_amr.hpp

// Counters used in allocations
extern uint gpu_vlasov_allocatedSize[];
extern uint gpu_acc_allocatedColumns;
extern uint gpu_acc_foundColumnsCount;

// Pointers used in pitch angle diffusion
// Host pointers
extern Real *host_bValues, *host_nu0Values, *host_bulkVX, *host_bulkVY, *host_bulkVZ, *host_Ddt;
extern Realf *host_sparsity, *dev_densityPreAdjust, *dev_densityPostAdjust;
extern size_t *host_cellIdxStartCutoff, *host_smallCellIdxArray, *host_remappedCellIdxArray; // remappedCellIdxArray tells the position of the cell index in the sequence instead of the actual index
// Device pointers
extern Real *dev_bValues, *dev_nu0Values, *dev_bulkVX, *dev_bulkVY, *dev_bulkVZ, *dev_Ddt, *dev_potentialDdtValues;
extern Realf *dev_fmu, *dev_dfdt_mu, *dev_sparsity;
extern int *dev_fcount, *dev_cellIdxKeys;
extern size_t *dev_smallCellIdxArray, *dev_remappedCellIdxArray, *dev_cellIdxStartCutoff, *dev_cellIdxArray, *dev_velocityIdxArray;
// Counters
extern size_t latestNumberOfLocalCellsPitchAngle;
extern int latestNumberOfVelocityCellsPitchAngle;
extern bool memoryHasBeenAllocatedPitchAngle;

#endif
