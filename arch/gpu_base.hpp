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

#ifndef GPU_BASE_H
#define GPU_BASE_H

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "arch_device_api.h"

// Extra profiling stream synchronizations?
#define SSYNC CHK_ERR( gpuStreamSynchronize(stream) )
//#define SSYNC

#include <stdio.h>
#include "include/splitvector/splitvec.h"
#include "include/hashinator/hashinator.h"
#include "../definitions.h"
#include "../vlasovsolver/vec.h"
#include "../velocity_mesh_parameters.h"
#include <phiprof.hpp>

static const double BLOCK_ALLOCATION_PADDING = 2.5;
static const double BLOCK_ALLOCATION_FACTOR = 1.8;

// Extern flags
extern bool needAttachedStreams;
extern bool doPrefetches;

#define DIMS 1
#define MAXCPUTHREADS 64

void gpu_init_device();
void gpu_clear_device();
gpuStream_t gpu_getStream();
gpuStream_t gpu_getPriorityStream();
uint gpu_getThread();
int gpu_getDevice();
void gpu_vlasov_allocate(uint maxBlockCount);
void gpu_vlasov_deallocate();
uint gpu_vlasov_getAllocation();
void gpu_vlasov_allocate_perthread(uint cpuThreadID, uint blockAllocationCount);
void gpu_vlasov_deallocate_perthread(uint cpuThreadID);
void gpu_acc_allocate(uint maxBlockCount);
void gpu_acc_deallocate();
void gpu_acc_allocate_perthread(uint cpuThreadID, uint columnAllocationCount);
void gpu_acc_deallocate_perthread(uint cpuThreadID);

void gpu_compaction_deallocate();
void gpu_compaction_allocate(const uint vectorLength, const size_t bytesNeeded);
void gpu_compaction_allocate_vec_perthread(const uint cpuThreadID, const uint vectorLength);
void gpu_compaction_allocate_buf_perthread(const uint cpuThreadID, const size_t bytesNeeded);

void gpu_trans_allocate(cuint nAllCells=0, cuint sumOfLengths=0, cuint largestVmesh=0, cuint unionSetSize=0);
void gpu_trans_deallocate();

extern gpuStream_t gpuStreamList[];
extern gpuStream_t gpuPriorityStreamList[];

// Unified memory class for inheritance
class Managed {
public:
   void *operator new(size_t len) {
      void *ptr;
      gpuMallocManaged(&ptr, len);
      gpuDeviceSynchronize();
      return ptr;
   }

   void operator delete(void *ptr) {
      gpuDeviceSynchronize();
      gpuFree(ptr);
   }

   void* operator new[] (size_t len) {
      void *ptr;
      gpuMallocManaged(&ptr, len);
      gpuDeviceSynchronize();
      return ptr;
   }

   void operator delete[] (void* ptr) {
      gpuDeviceSynchronize();
      gpuFree(ptr);
   }

};

// Structs used by Vlasov Acceleration semi-Lagrangian solver
struct Column {
   int valuesOffset;                              // Source data values
   size_t targetBlockOffsets[MAX_BLOCKS_PER_DIM]; // Target data array offsets
   int nblocks;                                   // Number of blocks in this column
   int minBlockK,maxBlockK;                       // Column parallel coordinate limits
   int kBegin;                                    // Actual un-sheared starting block index
   int i,j;                                       // Blocks' perpendicular coordinates
};

struct ColumnOffsets {
   split::SplitVector<uint> columnBlockOffsets; // indexes where columns start (in blocks, length totalColumns)
   split::SplitVector<uint> columnNumBlocks; // length of column (in blocks, length totalColumns)
   split::SplitVector<uint> setColumnOffsets; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
   split::SplitVector<uint> setNumColumns; // how many columns in set of columns (length nColumnSets)

   ColumnOffsets(uint nColumns) {
      columnBlockOffsets.resize(nColumns);
      columnNumBlocks.resize(nColumns);
      setColumnOffsets.resize(nColumns);
      setNumColumns.resize(nColumns);
      columnBlockOffsets.clear();
      columnNumBlocks.clear();
      setColumnOffsets.clear();
      setNumColumns.clear();
      // These vectors themselves are not in unified memory, just their content data,
      // so we call a regular optimizeGPU without the unified flag.
      gpuStream_t stream = gpu_getStream();
      columnBlockOffsets.optimizeGPU(stream);
      columnNumBlocks.optimizeGPU(stream);
      setColumnOffsets.optimizeGPU(stream);
      setNumColumns.optimizeGPU(stream);
   }
   void prefetchDevice(gpuStream_t stream) {
      columnBlockOffsets.optimizeGPU(stream);
      columnNumBlocks.optimizeGPU(stream);
      setColumnOffsets.optimizeGPU(stream);
      setNumColumns.optimizeGPU(stream);
   }
};

// Device data variables, to be allocated in good time. Made into an array so that each thread has their own pointer.
extern vmesh::LocalID *gpu_GIDlist[];
extern vmesh::LocalID *gpu_LIDlist[];
extern vmesh::GlobalID *gpu_BlocksID_mapped[];
extern vmesh::GlobalID *gpu_BlocksID_mapped_sorted[];
extern vmesh::GlobalID *gpu_LIDlist_unsorted[];
extern vmesh::LocalID *gpu_columnNBlocks[];

extern Vec *gpu_blockDataOrdered[];
extern uint *gpu_cell_indices_to_id[];
extern uint *gpu_block_indices_to_id[];
extern uint *gpu_vcell_transpose;

extern void *gpu_RadixSortTemp[];
extern uint gpu_acc_RadixSortTempSize[];

extern Real *returnReal[];
extern Realf *returnRealf[];
extern vmesh::LocalID *returnLID[];

extern Column *gpu_columns[];
extern ColumnOffsets *cpu_columnOffsetData[];
extern ColumnOffsets *gpu_columnOffsetData[];

// SplitVectors and buffers for use in stream compaction
extern split::SplitVector<vmesh::GlobalID> *vbwcl_gather[];
extern split::SplitVector<vmesh::GlobalID> *vbwncl_gather[];
extern void *compaction_buffer[];

// SplitVector information structs for use in fetching sizes and capacities without page faulting
extern split::SplitInfo *info_1[];
extern split::SplitInfo *info_2[];
extern split::SplitInfo *info_3[];
extern split::SplitInfo *info_4[];
extern Hashinator::MapInfo *info_m[];

// Vectors and set for use in translation, actually declared in vlasovsolver/gpu_trans_map_amr.hpp
// to sidestep compilation errors
// extern split::SplitVector<vmesh::VelocityMesh*> *allVmeshPointer;
// extern split::SplitVector<vmesh::VelocityMesh*> *allPencilsMeshes;
// extern split::SplitVector<vmesh::VelocityBlockContainer*> *allPencilsContainers;
// extern split::SplitVector<vmesh::GlobalID> *unionOfBlocks;
// extern Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet;

// Counters used in allocations
extern uint gpu_vlasov_allocatedSize;
extern uint gpu_acc_allocatedColumns;
extern uint gpu_acc_columnContainerSize;
extern uint gpu_acc_foundColumnsCount;

#endif