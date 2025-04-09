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

//Scale sizes based on WID value
#define INIT_VMESH_SIZE (32768/WID3)
#define INIT_MAP_SIZE (16 - WID)

static const uint VLASOV_BUFFER_MINBLOCKS = 32768/WID3;
static const uint VLASOV_BUFFER_MINCOLUMNS = 2000/WID;
static const double BLOCK_ALLOCATION_PADDING = 1.2;
static const double BLOCK_ALLOCATION_FACTOR = 1.1;

// buffers need to be larger for translation to allow proper parallelism
static const int TRANSLATION_BUFFER_ALLOCATION_FACTOR = 5;

#define DIMS 1
#define MAXCPUTHREADS 64

void gpu_init_device();
void gpu_clear_device();
gpuStream_t gpu_getStream();
gpuStream_t gpu_getPriorityStream();
uint gpu_getThread();
uint gpu_getMaxThreads();
int gpu_getDevice();
int gpu_reportMemory(const size_t local_cap=0, const size_t ghost_cap=0, const size_t local_size=0, const size_t ghost_size=0);

void gpu_vlasov_allocate(uint maxBlockCount);
void gpu_vlasov_deallocate();
void gpu_vlasov_allocate_perthread(uint cpuThreadID, uint blockAllocationCount);
void gpu_vlasov_deallocate_perthread(uint cpuThreadID);
uint gpu_vlasov_getAllocation();
uint gpu_vlasov_getSmallestAllocation();

void gpu_batch_allocate(uint nCells=0, uint maxNeighbours=0);
void gpu_batch_deallocate(bool first=true, bool second=true);

void gpu_acc_allocate(uint maxBlockCount);
void gpu_acc_allocate_perthread(uint cpuThreadID, uint columnAllocationCount);
void gpu_acc_deallocate();
void gpu_acc_deallocate_perthread(uint cpuThreadID);

void gpu_trans_allocate(cuint nAllCells=0,
                        cuint sumOfLengths=0,
                        cuint largestVmesh=0,
                        cuint unionSetSize=0,
                        cuint transGpuBlocks=0,
                        cuint nPencils=0);
void gpu_trans_deallocate();

extern gpuStream_t gpuStreamList[];
extern gpuStream_t gpuPriorityStreamList[];

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
      // These vectors themselves are not in unified memory, just their content data
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
   int capacity() {
      return columnBlockOffsets.capacity()
         + columnNumBlocks.capacity()
         + setColumnOffsets.capacity()
         + setNumColumns.capacity();
   }
   int capacityInBytes() {
      return columnBlockOffsets.capacity() * sizeof(uint)
         + columnNumBlocks.capacity() * sizeof(uint)
         + setColumnOffsets.capacity() * sizeof(uint)
         + setNumColumns.capacity() * sizeof(uint)
         + 4 * sizeof(split::SplitVector<uint>);
   }
};

// Device data variables, to be allocated in good time. Made into an array so that each thread has their own pointer.
extern vmesh::GlobalID *gpu_GIDlist[];
extern vmesh::LocalID *gpu_LIDlist[];
extern vmesh::GlobalID *gpu_BlocksID_mapped[];
extern vmesh::GlobalID *gpu_BlocksID_mapped_sorted[];
extern vmesh::LocalID *gpu_LIDlist_unsorted[];
extern vmesh::LocalID *gpu_columnNBlocks[];

extern Vec *gpu_blockDataOrdered[];
extern uint *gpu_cell_indices_to_id[];
extern uint *gpu_block_indices_to_id[];
extern uint *gpu_vcell_transpose;

extern Vec** host_pencilOrderedPointers;
extern Vec** dev_pencilOrderedPointers;
extern Realf** dev_pencilBlockData;
extern uint* dev_pencilBlocksCount;

extern void *gpu_RadixSortTemp[];
extern size_t gpu_acc_RadixSortTempSize[];

extern Real *returnReal[];
extern Realf *returnRealf[];
extern vmesh::LocalID *returnLID[];
extern Real *host_returnReal[];
extern Realf *host_returnRealf[];
extern vmesh::LocalID *host_returnLID[];
extern vmesh::GlobalID *invalidGIDpointer;

extern Column *gpu_columns[];
extern ColumnOffsets *cpu_columnOffsetData[];
extern ColumnOffsets *gpu_columnOffsetData[];

// Hash map and splitvectors buffers used in block adjustment, actually declared in block_adjust_gpu.hpp
// to sidestep compilation errors
// extern vmesh::VelocityMesh** host_vmeshes, **dev_vmeshes;
// extern vmesh::VelocityBlockContainer** host_VBCs, **dev_VBCs;
// extern Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** host_allMaps, **dev_allMaps;
// extern split::SplitVector<vmesh::GlobalID> ** host_vbwcl_vec, **dev_vbwcl_vec;
// extern split::SplitVector<vmesh::GlobalID> ** host_lists_with_replace_new, **dev_lists_with_replace_new;
// extern split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_delete, **dev_lists_delete;
// extern split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_to_replace, **dev_lists_to_replace;
// extern split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_with_replace_old, **dev_lists_with_replace_old;
// extern split::SplitVector<vmesh::GlobalID> ** host_vbwcl_neigh, **dev_vbwcl_neigh;
// extern vmesh::LocalID* host_contentSizes, *dev_contentSizes;
// extern Real* host_minValues, *dev_minValues;
// extern Real* host_massLoss, *dev_massLoss;
// extern Real* host_mass, *dev_mass;

// SplitVector information structs for use in fetching sizes and capacities without page faulting
// extern split::SplitInfo *info_1[];
// extern split::SplitInfo *info_2[];
// extern split::SplitInfo *info_3[];
// extern split::SplitInfo *info_4[];
// extern Hashinator::MapInfo *info_m[];

// Vectors and set for use in translation, actually declared in vlasovsolver/gpu_trans_map_amr.hpp
// to sidestep compilation errors
// extern split::SplitVector<vmesh::VelocityMesh*> *allVmeshPointer;
// extern split::SplitVector<vmesh::VelocityMesh*> *allPencilsMeshes;
// extern split::SplitVector<vmesh::VelocityBlockContainer*> *allPencilsContainers;
// extern split::SplitVector<vmesh::GlobalID> *unionOfBlocks;
// extern Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet;

// Counters used in allocations
extern uint gpu_vlasov_allocatedSize[];
extern uint gpu_acc_allocatedColumns;
extern uint gpu_acc_columnContainerSize;
extern uint gpu_acc_foundColumnsCount;

#endif
