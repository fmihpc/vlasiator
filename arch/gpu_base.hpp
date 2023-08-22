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
int gpu_getDevice();
void gpu_vlasov_allocate (uint maxBlockCount);
uint gpu_vlasov_getAllocation();
void gpu_vlasov_allocate_perthread (uint cpuThreadID, uint blockAllocationCount);
void gpu_vlasov_deallocate_perthread (uint cpuThreadID);
void gpu_acc_allocate (uint maxBlockCount);
void gpu_acc_allocate_perthread (uint cpuThreadID, uint columnAllocationCount);
void gpu_acc_deallocate_perthread (uint cpuThreadID);


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

struct ColumnOffsets : public Managed {
   split::SplitVector<uint> columnBlockOffsets; // indexes where columns start (in blocks, length totalColumns)
   split::SplitVector<uint> columnNumBlocks; // length of column (in blocks, length totalColumns)
   split::SplitVector<uint> setColumnOffsets; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
   split::SplitVector<uint> setNumColumns; // how many columns in set of columns (length nColumnSets)
   gpuStream_t attachedStream;

   ColumnOffsets(uint nColumns) {
      columnBlockOffsets.resize(nColumns);
      columnNumBlocks.resize(nColumns);
      setColumnOffsets.resize(nColumns);
      setNumColumns.resize(nColumns);
      columnBlockOffsets.clear();
      columnNumBlocks.clear();
      setColumnOffsets.clear();
      setNumColumns.clear();
      attachedStream=0;
   }
   void gpu_attachToStream(gpuStream_t stream = 0) {
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
      CHK_ERR( gpuStreamAttachMemAsync(stream,this, 0,gpuMemAttachSingle) );
      columnBlockOffsets.streamAttach(stream);
      columnNumBlocks.streamAttach(stream);
      setColumnOffsets.streamAttach(stream);
      setNumColumns.streamAttach(stream);
   }
   void gpu_detachFromStream() {
      // Return if attaching is not needed
      if (!needAttachedStreams) {
         return;
      }
      if (attachedStream == 0) {
         // Already detached
         return;
      }
      attachedStream = 0;
      CHK_ERR( gpuStreamAttachMemAsync(0,this, 0,gpuMemAttachGlobal) );
      columnBlockOffsets.streamAttach(0,gpuMemAttachGlobal);
      columnNumBlocks.streamAttach(0,gpuMemAttachGlobal);
      setColumnOffsets.streamAttach(0,gpuMemAttachGlobal);
      setNumColumns.streamAttach(0,gpuMemAttachGlobal);
   }
   void gpu_advise() {
      int device = gpu_getDevice();
      gpuStream_t stream = gpu_getStream();
      columnBlockOffsets.memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      columnNumBlocks.memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      setColumnOffsets.memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      setNumColumns.memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      columnBlockOffsets.memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      columnNumBlocks.memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      setColumnOffsets.memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      setNumColumns.memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
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
extern uint *gpu_vcell_transpose[];

extern void *gpu_RadixSortTemp[];
extern uint gpu_acc_RadixSortTempSize[];

extern Real *returnReal[];
extern Realf *returnRealf[];
extern vmesh::LocalID *returnLID[];

extern Column *gpu_columns[];

// Unified (managed) memory variables
extern ColumnOffsets *unif_columnOffsetData[];

// Counters used in allocations
extern uint gpu_vlasov_allocatedSize;
extern uint gpu_acc_allocatedColumns;
extern uint gpu_acc_columnContainerSize;
extern uint gpu_acc_foundColumnsCount;

#endif
