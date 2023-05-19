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

#ifndef CUDA_CONTEXT_H
#define CUDA_CONTEXT_H

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#ifdef __CUDACC__
#define CUDA_DEV __device__
#else
#define CUDA_DEV
#endif

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "include/splitvector/splitvec.h"
#include "include/hashinator/hashinator.h"
#include "definitions.h"
#include "vlasovsolver/vec.h"
#include "velocity_mesh_parameters.h"

#include <stdio.h>

// Assumes CUDA Hardware with warpsize 32
#define FULL_MASK 0xffffffff

static const double BLOCK_ALLOCATION_PADDING = 2.5;
static const double BLOCK_ALLOCATION_FACTOR = 1.8;

// CUDA Error checking
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ));
static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

// Unified memory class for inheritance
class Managed {
public:
   void *operator new(size_t len) {
      void *ptr;
      cudaMallocManaged(&ptr, len);
      cudaDeviceSynchronize();
      return ptr;
   }

   void operator delete(void *ptr) {
      cudaDeviceSynchronize();
      cudaFree(ptr);
   }

   void* operator new[] (size_t len) {
      void *ptr;
      cudaMallocManaged(&ptr, len);
      cudaDeviceSynchronize();
      return ptr;
   }

   void operator delete[] (void* ptr) {
      cudaDeviceSynchronize();
      cudaFree(ptr);
   }

};

// Predicate rule for hashinator for extracting all valid elements
template <typename T, typename U>
struct Rule{
   Rule(){}
   __host__ __device__
   inline bool operator()( Hashinator::hash_pair<T,U>& element)const{
      if (element.first!=std::numeric_limits<vmesh::GlobalID>::max() && element.first!=std::numeric_limits<vmesh::GlobalID>::max()-1  ){return true;}
      //KEY_TYPE EMPTYBUCKET = std::numeric_limits<KEY_TYPE>::max(),
      //   KEY_TYPE TOMBSTONE = EMPTYBUCKET - 1,
      //if (element.first< 1000 ){return true;}
      return false;
   }
};

#define DIMS 1
#ifndef CUDABLOCKS
#  define CUDABLOCKS (108)
#endif
#ifndef CUDATHREADS
#  define CUDATHREADS (32) // NVIDIA: 32 AMD: 64
#endif

#define MAXCPUTHREADS 64

void cuda_init_device();
void cuda_set_device();
void cuda_clear_device();
cudaStream_t cuda_getStream();
cudaStream_t cuda_getPriorityStream();
int cuda_getDevice();
void cuda_vlasov_allocate (uint maxBlockCount);
void cuda_vlasov_allocate_perthread (uint cpuThreadID, uint blockAllocationCount);
void cuda_vlasov_deallocate_perthread (uint cpuThreadID);
void cuda_acc_allocate (uint maxBlockCount);
void cuda_acc_allocate_perthread (uint cpuThreadID, uint columnAllocationCount);
void cuda_acc_deallocate_perthread (uint cpuThreadID);

// Extern flag for stream attaching
extern cudaStream_t cudaStreamList[];
extern cudaStream_t cudaPriorityStreamList[];
//extern cudaStream_t cudaBaseStream;
extern bool needAttachedStreams;

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
   cudaStream_t attachedStream;

   ColumnOffsets(uint nColumns) {
      columnBlockOffsets.resize(nColumns);
      columnNumBlocks.resize(nColumns);
      setColumnOffsets.resize(nColumns);
      setNumColumns.resize(nColumns);
      columnBlockOffsets.clear();
      columnNumBlocks.clear();
      setColumnOffsets.clear();
      setNumColumns.clear();
      //HANDLE_ERROR( cudaMemAdvise(this, sizeof(ColumnOffsets), cudaMemAdviseSetPreferredLocation, cuda_getDevice()) );
      attachedStream=0;
   }
   void dev_attachToStream(cudaStream_t stream = 0) {
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
      HANDLE_ERROR( cudaStreamAttachMemAsync(stream,this, 0,cudaMemAttachSingle) );
      columnBlockOffsets.streamAttach(stream);
      columnNumBlocks.streamAttach(stream);
      setColumnOffsets.streamAttach(stream);
      setNumColumns.streamAttach(stream);
   }
   void dev_detachFromStream() {
      // Return if attaching is not needed
      if (!needAttachedStreams) {
         return;
      }
      if (attachedStream == 0) {
         // Already detached
         return;
      }
      attachedStream = 0;
      HANDLE_ERROR( cudaStreamAttachMemAsync(0,this, 0,cudaMemAttachGlobal) );
      columnBlockOffsets.streamAttach(0,cudaMemAttachGlobal);
      columnNumBlocks.streamAttach(0,cudaMemAttachGlobal);
      setColumnOffsets.streamAttach(0,cudaMemAttachGlobal);
      setNumColumns.streamAttach(0,cudaMemAttachGlobal);
   }
};

// Device data variables, to be allocated in good time. Made into an array so that each thread has their own pointer.
extern vmesh::LocalID *dev_GIDlist[];
extern vmesh::LocalID *dev_LIDlist[];
extern vmesh::GlobalID *dev_BlocksID_mapped[];
extern vmesh::GlobalID *dev_BlocksID_mapped_sorted[];
extern vmesh::GlobalID *dev_LIDlist_unsorted[];
extern vmesh::LocalID *dev_columnNBlocks[];

extern Vec *dev_blockDataOrdered[];

extern Realf *returnRealf[];
extern uint *dev_cell_indices_to_id[];
extern uint *dev_block_indices_to_id[];
extern uint *dev_vcell_transpose[];

// Unified (managed) memory variables
extern ColumnOffsets *unif_columnOffsetData[];
extern Column *unif_columns[];

// Counters used in allocations
extern uint cuda_vlasov_allocatedSize;
extern uint cuda_acc_allocatedColumns;
extern uint cuda_acc_columnContainerSize;
extern uint cuda_acc_foundColumnsCount;

#endif
