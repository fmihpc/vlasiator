/*
 * This file is part of Vlasiator.
 * Copyright 2010-2022 Finnish Meteorological Institute and University of Helsinki
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <iostream>

#include "common.h"
#include "mpi.h"

#include "cuda_context.cuh"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

// #define MAXCPUTHREADS 64 now in cuda_context.hpp

int myDevice;
int myRank;

// Allocate pointers for per-thread memory regions
cudaStream_t cudaStreamList[MAXCPUTHREADS];
cudaStream_t cudaPriorityStreamList[MAXCPUTHREADS];

Real *returnReal[MAXCPUTHREADS];
Realf *returnRealf[MAXCPUTHREADS];
vmesh::LocalID *returnLID[MAXCPUTHREADS];

bool needAttachedStreams = false;
bool doPrefetches=true; // only non-crucial prefetches are behind this check
// Note: disabling prefetches brings in strange memory errors and crashes (June 2023)

uint *dev_cell_indices_to_id[MAXCPUTHREADS];
uint *dev_block_indices_to_id[MAXCPUTHREADS];
uint *dev_vcell_transpose[MAXCPUTHREADS];

ColumnOffsets *unif_columnOffsetData[MAXCPUTHREADS];
Column *dev_columns[MAXCPUTHREADS];

Vec *dev_blockDataOrdered[MAXCPUTHREADS];

vmesh::LocalID *dev_GIDlist[MAXCPUTHREADS];
vmesh::LocalID *dev_LIDlist[MAXCPUTHREADS];
vmesh::GlobalID *dev_BlocksID_mapped[MAXCPUTHREADS];
vmesh::GlobalID *dev_BlocksID_mapped_sorted[MAXCPUTHREADS];
vmesh::GlobalID *dev_LIDlist_unsorted[MAXCPUTHREADS];
vmesh::LocalID *dev_columnNBlocks[MAXCPUTHREADS];
void *dev_RadixSortTemp[MAXCPUTHREADS];
uint cuda_acc_RadixSortTempSize[MAXCPUTHREADS] = {0};

// Memory allocation flags and values.
uint cuda_vlasov_allocatedSize = 0;
uint cuda_acc_allocatedColumns = 0;
uint cuda_acc_columnContainerSize = 0;
uint cuda_acc_foundColumnsCount = 0;

__host__ void cuda_init_device() {

#ifdef _OPENMP
   const uint maxThreads = omp_get_max_threads();
#else
   const uint maxThreads = 1
#endif

   int deviceCount;
   //HANDLE_ERROR( cudaFree(0));
   HANDLE_ERROR( cudaGetDeviceCount(&deviceCount) );
   printf("CUDA device count %d with %d threads/streams\n",deviceCount,maxThreads);

   /* Create communicator with one rank per compute node to identify which GPU to use */
   int amps_size;
   int amps_rank;
   int amps_node_rank;
   int amps_node_size;
   // int amps_write_rank;
   // int amps_write_size;
   MPI_Comm amps_CommWorld = MPI_COMM_NULL;
   MPI_Comm amps_CommNode = MPI_COMM_NULL;

   MPI_Comm_dup(MPI_COMM_WORLD, &amps_CommWorld);
   MPI_Comm_size(amps_CommWorld, &amps_size);
   MPI_Comm_rank(amps_CommWorld, &amps_rank);

   /* Create communicator with one rank per compute node */
#if MPI_VERSION >= 3
   MPI_Comm_split_type(amps_CommWorld, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &amps_CommNode);
#else
   /* Split the node level communicator based on Adler32 hash keys of processor name */
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int namelen;
   MPI_Get_processor_name(processor_name, &namelen);
   uint32_t checkSum = Adler32((unsigned char*)processor_name, namelen);
   /* Comm split only accepts non-negative numbers */
   /* Not super great for hashing purposes but hoping MPI-3 code will be used on most cases */
   checkSum &= INT_MAX;
   MPI_Comm_split(amps_CommWorld, checkSum, amps_rank, &amps_CommNode);
#endif

   MPI_Comm_rank(amps_CommNode, &amps_node_rank);
   MPI_Comm_size(amps_CommNode, &amps_node_size);
   std::cerr << "(Grid) rank " << amps_rank << " is noderank "<< amps_node_rank << " of "<< amps_node_size << std::endl;
   myRank = amps_node_rank;

   if (amps_node_rank >= deviceCount) {
      std::cerr<<"Error, attempting to use CUDA device beyond available count!"<<std::endl;
      abort();
   }
   if (amps_node_size > deviceCount) {
      std::cerr<<"Error, MPI tasks per node exceeds available CUDA device count!"<<std::endl;
      abort();
   }
   HANDLE_ERROR( cudaSetDevice(amps_node_rank) );
   HANDLE_ERROR( cudaDeviceSynchronize() );
   HANDLE_ERROR( cudaGetDevice(&myDevice) );

   // Query device capabilities
   int supportedMode;
   HANDLE_ERROR( cudaDeviceGetAttribute (&supportedMode, cudaDevAttrConcurrentManagedAccess, myDevice) );
   if (supportedMode==0) {
      printf("Warning! Current CUDA device does not support concurrent managed memory access from several streams.\n");
      needAttachedStreams = true;
   }

   // Pre-generate streams, allocate return pointers
   int *leastPriority = new int; // likely 0
   int *greatestPriority = new int; // likely -1
   HANDLE_ERROR( cudaDeviceGetStreamPriorityRange (leastPriority, greatestPriority) );
   if (*leastPriority==*greatestPriority) {
      printf("Warning when initializing CUDA streams: minimum and maximum stream priority are identical! %d == %d \n",*leastPriority, *greatestPriority);
   }
   for (uint i=0; i<maxThreads; ++i) {
      HANDLE_ERROR( cudaStreamCreateWithPriority(&(cudaStreamList[i]), cudaStreamDefault, *leastPriority) );
      HANDLE_ERROR( cudaStreamCreateWithPriority(&(cudaPriorityStreamList[i]), cudaStreamDefault, *greatestPriority) );
      HANDLE_ERROR( cudaMalloc((void**)&returnReal[i], 8*sizeof(Real)) );
      HANDLE_ERROR( cudaMalloc((void**)&returnRealf[i], 8*sizeof(Realf)) );
      HANDLE_ERROR( cudaMalloc((void**)&returnLID[i], 8*sizeof(vmesh::LocalID)) );
   }

   // Using just a single context for whole MPI task
}

__host__ void cuda_clear_device() {
   // Destroy streams
#ifdef _OPENMP
   const uint maxThreads = omp_get_max_threads();
#else
   const uint maxThreads = 1
#endif
   for (uint i=0; i<maxThreads; ++i) {
      HANDLE_ERROR( cudaStreamDestroy(cudaStreamList[i]) );
      HANDLE_ERROR( cudaStreamDestroy(cudaPriorityStreamList[i]) );
      HANDLE_ERROR( cudaFree(returnReal[i]) );
      HANDLE_ERROR( cudaFree(returnRealf[i]) );
      HANDLE_ERROR( cudaFree(returnLID[i]) );
   }
}

__host__ cudaStream_t cuda_getStream() {
#ifdef _OPENMP
   const uint thread_id = omp_get_thread_num();
#else
   const uint thread_id = 0;
#endif
   return cudaStreamList[thread_id];
}

__host__ cudaStream_t cuda_getPriorityStream() {
#ifdef _OPENMP
   const uint thread_id = omp_get_thread_num();
#else
   const uint thread_id = 0;
#endif
   return cudaPriorityStreamList[thread_id];
}

__host__ int cuda_getDevice() {
   int device;
   HANDLE_ERROR( cudaGetDevice(&device) );
   return device;
}

/*
   Top-level GPU memory allocation function.
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void cuda_vlasov_allocate (
   uint maxBlockCount
   ) {
   // Always prepare for at least 2500 blocks (affects also translation parallelism)
   const uint maxBlocksPerCell = maxBlockCount > 2500 ? maxBlockCount : 2500;
   // Check if we already have allocated enough memory?
   if (cuda_vlasov_allocatedSize > maxBlocksPerCell * BLOCK_ALLOCATION_FACTOR) {
      return;
   }
   // If not, add extra padding
   const uint newSize = maxBlocksPerCell * BLOCK_ALLOCATION_PADDING;

   // Deallocate before allocating new memory
#ifdef _OPENMP
   const uint maxNThreads = omp_get_max_threads();
#else
   const uint maxNThreads = 1;
#endif
   for (uint i=0; i<maxNThreads; ++i) {
      if (cuda_vlasov_allocatedSize > 0) {
         cuda_vlasov_deallocate_perthread(i);
      }
      cuda_vlasov_allocate_perthread(i, newSize);
   }
   cuda_vlasov_allocatedSize = newSize;
}

__host__ uint cuda_vlasov_getAllocation() {
   return cuda_vlasov_allocatedSize;
}

__host__ void cuda_vlasov_allocate_perthread (
   uint cpuThreadID,
   uint blockAllocationCount
   ) {
   // Mallocs should be in increments of 256 bytes. WID3 is at least 64, and len(Realf) is at least 4, so this is true.
   HANDLE_ERROR( cudaMalloc((void**)&dev_cell_indices_to_id[cpuThreadID], 3*sizeof(uint)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_block_indices_to_id[cpuThreadID], 3*sizeof(uint)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_vcell_transpose[cpuThreadID], WID3*sizeof(uint)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_blockDataOrdered[cpuThreadID], blockAllocationCount * (WID3 / VECL) * sizeof(Vec)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_BlocksID_mapped[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_BlocksID_mapped_sorted[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_LIDlist_unsorted[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_LIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::LocalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_GIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::LocalID)) );
}

__host__ void cuda_vlasov_deallocate_perthread (
   uint cpuThreadID
   ) {
   cuda_vlasov_allocatedSize = 0;
   HANDLE_ERROR( cudaFree(dev_cell_indices_to_id[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_block_indices_to_id[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_vcell_transpose[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_blockDataOrdered[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_BlocksID_mapped[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_BlocksID_mapped_sorted[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_LIDlist_unsorted[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_LIDlist[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_GIDlist[cpuThreadID]) );
}

/*
   Top-level GPU memory allocation function for acceleration-specific column data.
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void cuda_acc_allocate (
   uint maxBlockCount
   ) {
   uint requiredColumns;
   // Has the acceleration solver already figured out how many columns we have?
   if (cuda_acc_foundColumnsCount > 0) {
      // Always prepare for at least 500 columns
      requiredColumns = cuda_acc_foundColumnsCount > 500 ? cuda_acc_foundColumnsCount : 500;
   } else {
      // The worst case scenario for columns is with every block having content but no neighbours, creating up
      // to maxBlockCount columns with each needing three blocks (one value plus two for padding).

      // Use column count estimate from size of v-space? (quite large)
      // const uint estimatedColumns = 3 * std::pow(
      //    (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[0]
      //    * (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[1]
      //    * (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[2], 0.667);

      // Use column count estimate from blocks count
      const uint estimatedColumns = 3 * std::pow(maxBlockCount, 0.667);
      // Always prepare for at least 500 columns
      requiredColumns = estimatedColumns > 500 ? estimatedColumns : 500;
   }
   // Check if we already have allocated enough memory?
   if (cuda_acc_allocatedColumns > requiredColumns * BLOCK_ALLOCATION_FACTOR) {
      //std::cerr<<"CUDA_ACC_ALLOCATE: return "<<std::endl;
      return;
   }
   // If not, add extra padding
   const uint newSize = requiredColumns * BLOCK_ALLOCATION_PADDING;

   // Deallocate before allocating new memory
#ifdef _OPENMP
   const uint maxNThreads = omp_get_max_threads();
#else
   const uint maxNThreads = 1;
#endif
   if (cuda_acc_allocatedColumns > 0) {
      for (uint i=0; i<maxNThreads; ++i) {
         cuda_acc_deallocate_perthread(i);
      }
   }
   for (uint i=0; i<maxNThreads; ++i) {
      cuda_acc_allocate_perthread(i,newSize);
   }
   cuda_acc_allocatedColumns = newSize;
}

__host__ void cuda_acc_allocate_perthread (
   uint cpuThreadID,
   uint columnAllocationCount
   ) {
   // Unified memory; columndata contains several splitvectors.
   unif_columnOffsetData[cpuThreadID] = new ColumnOffsets(columnAllocationCount); // inherits managed
   HANDLE_ERROR( cudaMalloc((void**)&dev_columns[cpuThreadID], columnAllocationCount*sizeof(Column)) );

   // Potential ColumnSet block count container
   const uint c0 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[0];
   const uint c1 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[1];
   const uint c2 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[2];
   std::array<uint, 3> s = {c0,c1,c2};
   std::sort(s.begin(), s.end());
   cuda_acc_columnContainerSize = c2*c1;
   HANDLE_ERROR( cudaMalloc((void**)&dev_columnNBlocks[cpuThreadID], cuda_acc_columnContainerSize*sizeof(vmesh::LocalID)) );
}

__host__ void cuda_acc_deallocate_perthread (
   uint cpuThreadID
   ) {
   cuda_acc_allocatedColumns = 0;
   cuda_acc_columnContainerSize = 0;
   HANDLE_ERROR( cudaFree(dev_columnNBlocks[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_columns[cpuThreadID]) );
   delete unif_columnOffsetData[cpuThreadID];
}
