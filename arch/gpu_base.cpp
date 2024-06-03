/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
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

#include <stdio.h>
#include <iostream>
#include "common.h"
#include "mpi.h"

#include "gpu_base.hpp"
#include "../vlasovsolver/gpu_moments.h"
#include "../velocity_mesh_gpu.h"
#include "../velocity_block_container.h"
#include "../vlasovsolver/cpu_trans_pencils.hpp"

// #define MAXCPUTHREADS 64 now in gpu_base.hpp

int myDevice;
int myRank;

// Allocate pointers for per-thread memory regions
gpuStream_t gpuStreamList[MAXCPUTHREADS];
gpuStream_t gpuPriorityStreamList[MAXCPUTHREADS];

Real *returnReal[MAXCPUTHREADS];
Realf *returnRealf[MAXCPUTHREADS];
vmesh::LocalID *returnLID[MAXCPUTHREADS];

uint *gpu_cell_indices_to_id[MAXCPUTHREADS];
uint *gpu_block_indices_to_id[MAXCPUTHREADS];
uint *gpu_vcell_transpose; // only one needed, not one per thread

// Pointers to buffers used in acceleration
ColumnOffsets *cpu_columnOffsetData[MAXCPUTHREADS];
ColumnOffsets *gpu_columnOffsetData[MAXCPUTHREADS];
Column *gpu_columns[MAXCPUTHREADS];
Vec *gpu_blockDataOrdered[MAXCPUTHREADS];
vmesh::LocalID *gpu_GIDlist[MAXCPUTHREADS];
vmesh::LocalID *gpu_LIDlist[MAXCPUTHREADS];
vmesh::GlobalID *gpu_BlocksID_mapped[MAXCPUTHREADS];
vmesh::GlobalID *gpu_BlocksID_mapped_sorted[MAXCPUTHREADS];
vmesh::GlobalID *gpu_LIDlist_unsorted[MAXCPUTHREADS];
vmesh::LocalID *gpu_columnNBlocks[MAXCPUTHREADS];
vmesh::GlobalID *invalidGIDpointer = 0;
void *gpu_RadixSortTemp[MAXCPUTHREADS];
uint gpu_acc_RadixSortTempSize[MAXCPUTHREADS] = {0};

// Hash map and splitvectors used in block adjustment
split::SplitVector<vmesh::GlobalID> *gpu_list_with_replace_new[MAXCPUTHREADS];
split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> *gpu_list_delete[MAXCPUTHREADS];
split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> *gpu_list_to_replace[MAXCPUTHREADS];
split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> *gpu_list_with_replace_old[MAXCPUTHREADS];

// Vectors and set for use in translation
split::SplitVector<vmesh::VelocityMesh*> *allVmeshPointer;
split::SplitVector<vmesh::VelocityMesh*> *allPencilsMeshes;
split::SplitVector<vmesh::VelocityBlockContainer*> *allPencilsContainers;
split::SplitVector<vmesh::GlobalID> *unionOfBlocks;
Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet;
// counters for allocated sizes in translation
uint gpu_allocated_nAllCells = 0;
uint gpu_allocated_sumOfLengths = 0;
uint gpu_allocated_largestVmesh = 0;
uint gpu_allocated_unionSetSize = 0;

// Memory allocation flags and values (TODO make per-thread?).
uint gpu_blockadjust_allocatedSize[MAXCPUTHREADS] = {0};
uint gpu_vlasov_allocatedSize[MAXCPUTHREADS] = {0};
uint gpu_acc_allocatedColumns = 0;
uint gpu_acc_columnContainerSize = 0;
uint gpu_acc_foundColumnsCount = 0;

// SplitVector information structs for use in fetching sizes and capacities without page faulting
// split::SplitInfo *info_1[MAXCPUTHREADS];
// split::SplitInfo *info_2[MAXCPUTHREADS];
// split::SplitInfo *info_3[MAXCPUTHREADS];
// split::SplitInfo *info_4[MAXCPUTHREADS];
// Hashinator::MapInfo *info_m[MAXCPUTHREADS];

__host__ uint gpu_getThread() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}
__host__ uint gpu_getMaxThreads() {
#ifdef _OPENMP
   return omp_get_max_threads();
#else
   return 1;
#endif
}

__host__ void gpu_init_device() {
   const uint maxNThreads = gpu_getMaxThreads();
   // int deviceCount;
   // CHK_ERR( gpuFree(0));
   // CHK_ERR( gpuGetDeviceCount(&deviceCount) );
   // printf("GPU device count %d with %d threads/streams\n",deviceCount,maxThreads);

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

   // if (amps_node_rank >= deviceCount) {
   //    std::cerr<<"Error, attempting to use GPU device beyond available count!"<<std::endl;
   //    abort();
   // }
   // if (amps_node_size > deviceCount) {
   //    std::cerr<<"Error, MPI tasks per node exceeds available GPU device count!"<<std::endl;
   //    abort();
   // }
   // CHK_ERR( gpuSetDevice(amps_node_rank) );
   // CHK_ERR( gpuDeviceSynchronize() );
   CHK_ERR( gpuGetDevice(&myDevice) );

   // Query device capabilities (only for CUDA, not needed for HIP)
   #if defined(USE_GPU) && defined(__CUDACC__)
   int supportedMode;
   CHK_ERR( cudaDeviceGetAttribute (&supportedMode, cudaDevAttrConcurrentManagedAccess, myDevice) );
   if (supportedMode==0) {
      printf("Warning! Current GPU device does not support concurrent managed memory access from several streams.\n");
      abort();
   }
   #endif

   // Pre-generate streams, allocate return pointers
   int *leastPriority = new int; // likely 0
   int *greatestPriority = new int; // likely -1
   CHK_ERR( gpuDeviceGetStreamPriorityRange (leastPriority, greatestPriority) );
   if (*leastPriority==*greatestPriority) {
      printf("Warning when initializing GPU streams: minimum and maximum stream priority are identical! %d == %d \n",*leastPriority, *greatestPriority);
   }
   for (uint i=0; i<maxNThreads; ++i) {
      CHK_ERR( gpuStreamCreateWithPriority(&(gpuStreamList[i]), gpuStreamDefault, *leastPriority) );
      CHK_ERR( gpuStreamCreateWithPriority(&(gpuPriorityStreamList[i]), gpuStreamDefault, *greatestPriority) );
      CHK_ERR( gpuMalloc((void**)&returnReal[i], 8*sizeof(Real)) );
      CHK_ERR( gpuMalloc((void**)&returnRealf[i], 8*sizeof(Realf)) );
      CHK_ERR( gpuMalloc((void**)&returnLID[i], 8*sizeof(vmesh::LocalID)) );
      // CHK_ERR( gpuMallocHost((void **) &info_1[i], sizeof(split::SplitInfo)) );
      // CHK_ERR( gpuMallocHost((void **) &info_2[i], sizeof(split::SplitInfo)) );
      // CHK_ERR( gpuMallocHost((void **) &info_3[i], sizeof(split::SplitInfo)) );
      // CHK_ERR( gpuMallocHost((void **) &info_4[i], sizeof(split::SplitInfo)) );
      // CHK_ERR( gpuMallocHost((void **) &info_m[i], sizeof(Hashinator::MapInfo)) );
   }
   CHK_ERR( gpuMalloc((void**)&gpu_vcell_transpose, WID3*sizeof(uint)) );
   CHK_ERR( gpuMalloc((void**)&invalidGIDpointer, sizeof(vmesh::GlobalID)) );
   vmesh::GlobalID invalidGIDvalue = vmesh::INVALID_GLOBALID;
   CHK_ERR( gpuMemcpy(invalidGIDpointer, &invalidGIDvalue, sizeof(vmesh::LocalID), gpuMemcpyHostToDevice) );
   CHK_ERR( gpuDeviceSynchronize() );

   // Using just a single context for whole MPI task
}

__host__ void gpu_clear_device() {
   // Deallocate temporary buffers
   gpu_acc_deallocate();
   gpu_vlasov_deallocate();
   gpu_trans_deallocate();
   gpu_moments_deallocate();
   // Destroy streams
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      CHK_ERR( gpuStreamDestroy(gpuStreamList[i]) );
      CHK_ERR( gpuStreamDestroy(gpuPriorityStreamList[i]) );
      CHK_ERR( gpuFree(returnReal[i]) );
      CHK_ERR( gpuFree(returnRealf[i]) );
      CHK_ERR( gpuFree(returnLID[i]) );
   }
   CHK_ERR( gpuFree(gpu_vcell_transpose) );
   CHK_ERR( gpuDeviceSynchronize() );
}

__host__ gpuStream_t gpu_getStream() {
   return gpuStreamList[gpu_getThread()];
}

__host__ gpuStream_t gpu_getPriorityStream() {
   const uint thread_id = gpu_getThread();
   return gpuPriorityStreamList[thread_id];
}

__host__ int gpu_getDevice() {
   int device;
   CHK_ERR( gpuGetDevice(&device) );
   return device;
}

/*
   Top-level GPU memory allocation function.
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void gpu_vlasov_allocate(
   uint maxBlockCount // Largest found vmesh size
   ) {
   // Always prepare for at least 2500 blocks
   const uint maxBlocksPerCell = maxBlockCount > 2500 ? maxBlockCount : 2500;
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_vlasov_allocate_perthread(i, maxBlocksPerCell);
   }
}

/* Deallocation at end of simulation */
__host__ void gpu_vlasov_deallocate() {
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_vlasov_deallocate_perthread(i);
   }
}

__host__ uint gpu_vlasov_getAllocation() {
   const uint cpuThreadID = gpu_getThread();
   return gpu_vlasov_allocatedSize[cpuThreadID];
}
__host__ uint gpu_vlasov_getSmallestAllocation() {
   uint smallestAllocation = std::numeric_limits<uint>::max();
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      const uint threadAllocation = gpu_vlasov_allocatedSize[i];
      if (threadAllocation < smallestAllocation) {
         smallestAllocation = threadAllocation;
      }
   }
   return smallestAllocation;
}

__host__ void gpu_vlasov_allocate_perthread(
   uint cpuThreadID,
   uint blockAllocationCount
   ) {
   // Check if we already have allocated enough memory?
   if (gpu_vlasov_allocatedSize[cpuThreadID] > blockAllocationCount * BLOCK_ALLOCATION_FACTOR) {
      return;
   }
   // Potential new allocation with extra padding
   const uint newSize = blockAllocationCount * BLOCK_ALLOCATION_PADDING;
   // Deallocate before new allocation
   gpu_vlasov_deallocate_perthread(cpuThreadID);

   // Mallocs should be in increments of 256 bytes. WID3 is at least 64, and len(Realf) is at least 4, so this is true.
   CHK_ERR( gpuMalloc((void**)&gpu_cell_indices_to_id[cpuThreadID], 3*sizeof(uint)) );
   CHK_ERR( gpuMalloc((void**)&gpu_block_indices_to_id[cpuThreadID], 3*sizeof(uint)) );
   CHK_ERR( gpuMalloc((void**)&gpu_blockDataOrdered[cpuThreadID], newSize * TRANSLATION_BUFFER_ALLOCATION_FACTOR * (WID3 / VECL) * sizeof(Vec)) );
   CHK_ERR( gpuMalloc((void**)&gpu_BlocksID_mapped[cpuThreadID], newSize*sizeof(vmesh::GlobalID)) );
   CHK_ERR( gpuMalloc((void**)&gpu_BlocksID_mapped_sorted[cpuThreadID], newSize*sizeof(vmesh::GlobalID)) );
   CHK_ERR( gpuMalloc((void**)&gpu_LIDlist_unsorted[cpuThreadID], newSize*sizeof(vmesh::GlobalID)) );
   CHK_ERR( gpuMalloc((void**)&gpu_LIDlist[cpuThreadID], newSize*sizeof(vmesh::LocalID)) );
   CHK_ERR( gpuMalloc((void**)&gpu_GIDlist[cpuThreadID], newSize*sizeof(vmesh::LocalID)) );
   // Store size of new allocation
   gpu_vlasov_allocatedSize[cpuThreadID] = newSize;
}

__host__ void gpu_vlasov_deallocate_perthread (
   uint cpuThreadID
   ) {
   if (gpu_vlasov_allocatedSize[cpuThreadID] == 0) {
      return;
   }
   CHK_ERR( gpuFree(gpu_cell_indices_to_id[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_block_indices_to_id[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_blockDataOrdered[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_BlocksID_mapped[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_BlocksID_mapped_sorted[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_LIDlist_unsorted[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_LIDlist[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_GIDlist[cpuThreadID]) );
   gpu_vlasov_allocatedSize[cpuThreadID] = 0;
}

/*
   Top-level GPU memory allocation function for acceleration-specific column data.
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void gpu_acc_allocate(
   uint maxBlockCount
   ) {
   uint requiredColumns;
   // Has the acceleration solver already figured out how many columns we have?
   if (gpu_acc_foundColumnsCount > 0) {
      // Always prepare for at least 500 columns
      requiredColumns = gpu_acc_foundColumnsCount > 500 ? gpu_acc_foundColumnsCount : 500;
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
   if (gpu_acc_allocatedColumns > requiredColumns * BLOCK_ALLOCATION_FACTOR) {
      return;
   }
   // If not, add extra padding
   const uint newSize = requiredColumns * BLOCK_ALLOCATION_PADDING;

   // Deallocate before allocating new memory
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_acc_deallocate_perthread(i);
   }
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_acc_allocate_perthread(i,newSize);
   }
   gpu_acc_allocatedColumns = newSize;
}

/* Deallocation at end of simulation */
__host__ void gpu_acc_deallocate() {
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_acc_deallocate_perthread(i);
   }
}

__host__ void gpu_acc_allocate_perthread(
   uint cpuThreadID,
   uint columnAllocationCount
   ) {
   // columndata contains several splitvectors. columnData is host/device, but splitvector contents are unified.
   if (columnAllocationCount > 0) {
      cpu_columnOffsetData[cpuThreadID] = new ColumnOffsets(columnAllocationCount);
      CHK_ERR( gpuMalloc((void**)&gpu_columnOffsetData[cpuThreadID], sizeof(ColumnOffsets)) );
      CHK_ERR( gpuMemcpy(gpu_columnOffsetData[cpuThreadID], cpu_columnOffsetData[cpuThreadID], sizeof(ColumnOffsets), gpuMemcpyHostToDevice));
      CHK_ERR( gpuMalloc((void**)&gpu_columns[cpuThreadID], columnAllocationCount*sizeof(Column)) );
   }
   // Potential ColumnSet block count container
   const uint c0 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[0];
   const uint c1 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[1];
   const uint c2 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[2];
   std::array<uint, 3> s = {c0,c1,c2};
   std::sort(s.begin(), s.end());
   gpu_acc_columnContainerSize = c2*c1;
   CHK_ERR( gpuMalloc((void**)&gpu_columnNBlocks[cpuThreadID], gpu_acc_columnContainerSize*sizeof(vmesh::LocalID)) );
}

__host__ void gpu_acc_deallocate_perthread(
   uint cpuThreadID
   ) {
   gpu_acc_allocatedColumns = 0;
   gpu_acc_columnContainerSize = 0;
   if (gpu_acc_allocatedColumns > 0) {
      CHK_ERR( gpuFree(gpu_columns[cpuThreadID]) );
      delete cpu_columnOffsetData[cpuThreadID];
      CHK_ERR( gpuFree(gpu_columnOffsetData[cpuThreadID]) );
   }
   CHK_ERR( gpuFree(gpu_columnNBlocks[cpuThreadID]) );
}

/*
   Top-level GPU memory allocation function for block adjustment containers.
 */
__host__ void gpu_blockadjust_allocate(
   uint maxBlockCount
   ) {
   // Always prepare for at least 100 blocks
   const uint maxBlocksPerCell = maxBlockCount > 100 ? maxBlockCount : 100;
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_blockadjust_allocate_perthread(i, maxBlocksPerCell);
   }
}

/* Deallocation at end of simulation */
__host__ void gpu_blockadjust_deallocate() {
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_blockadjust_deallocate_perthread(i);
   }
}
__host__ void gpu_blockadjust_allocate_perthread(
   uint cpuThreadID,
   uint blockAllocationCount
   ) {
   gpuStream_t stream = gpu_getStream();
   // Check if we already have allocated enough memory?
   if (gpu_blockadjust_allocatedSize[cpuThreadID] > blockAllocationCount * BLOCK_ALLOCATION_FACTOR) {
      //printf("Early return from allocation, thread %lu with blockAllocationCount %lu existing allocation counter %lu and vector capacity %lu \n",(long unsigned)cpuThreadID,(long unsigned)blockAllocationCount,(long unsigned)gpu_blockadjust_allocatedSize[cpuThreadID],(long unsigned)gpu_list_with_replace_new[cpuThreadID]->capacity());
      return;
   }

   // Potential new allocation with extra padding
   const uint newSize = blockAllocationCount * BLOCK_ALLOCATION_PADDING;
   // Deallocate before new allocation
   gpu_blockadjust_deallocate_perthread(cpuThreadID);

   gpu_list_with_replace_new[cpuThreadID] = new split::SplitVector<vmesh::GlobalID>(newSize);
   gpu_list_delete[cpuThreadID] = new split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(newSize);
   gpu_list_to_replace[cpuThreadID] = new split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(newSize);
   gpu_list_with_replace_old[cpuThreadID] = new split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(newSize);

   // Store size of new allocation
   gpu_blockadjust_allocatedSize[cpuThreadID] = newSize;
   //printf("Allocated buffers for thread %lu with newSize %lu and vector length %lu = %lu \n",(long unsigned)cpuThreadID,(long unsigned)newSize,(long unsigned)loopReserve,(long unsigned)gpu_list_with_replace_new[cpuThreadID]->capacity());
}

__host__ void gpu_blockadjust_deallocate_perthread (
   uint cpuThreadID
   ) {
   if (gpu_blockadjust_allocatedSize[cpuThreadID] == 0) {
      return;
   }
   delete gpu_list_with_replace_new[cpuThreadID];
   delete gpu_list_delete[cpuThreadID];
   delete gpu_list_to_replace[cpuThreadID];
   delete gpu_list_with_replace_old[cpuThreadID];
   gpu_blockadjust_allocatedSize[cpuThreadID] = 0;
}

/*
   Top-level GPU memory allocation function for translation-specific vectors
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void gpu_trans_allocate(
   cuint nAllCells,
   cuint sumOfLengths,
   cuint largestVmesh,
   cuint unionSetSize
   ) {
   gpuStream_t stream = gpu_getStream();
   // Vectors with one entry per cell (prefetch to host)
   if (nAllCells != 0) {
      if (gpu_allocated_nAllCells == 0) {
         // New allocation
         allVmeshPointer = new split::SplitVector<vmesh::VelocityMesh*>(nAllCells);
      } else {
         // Resize
         allVmeshPointer->clear();
         allVmeshPointer->optimizeCPU(stream);
         allVmeshPointer->resize(nAllCells,true);
      }
      // Leave on CPU
      gpu_allocated_nAllCells = nAllCells;
   }
   // Vectors with one entry per pencil cell (prefetch to host)
   if (sumOfLengths != 0) {
      if (gpu_allocated_sumOfLengths == 0) {
         // New allocations
         allPencilsMeshes = new split::SplitVector<vmesh::VelocityMesh*>(sumOfLengths);
         allPencilsContainers = new split::SplitVector<vmesh::VelocityBlockContainer*>(sumOfLengths);
      } else {
         // Resize
         // allPencilsMeshes->optimizeCPU(stream);
         // allPencilsContainers->optimizeCPU(stream);
         allPencilsMeshes->resize(sumOfLengths,true);
         allPencilsContainers->resize(sumOfLengths,true);
      }
      // Leave on CPU
      gpu_allocated_sumOfLengths = sumOfLengths;
   }
   // Set for collecting union of blocks (prefetched to device)
   if (largestVmesh != 0) {
      const vmesh::LocalID HashmapReqSize = ceil(log2((int)largestVmesh)) +3;
      if (gpu_allocated_largestVmesh == 0) {
         // New allocation
         unionOfBlocksSet = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(HashmapReqSize);
         unionOfBlocksSet->optimizeGPU(stream);
      } else {
         // Ensure allocation
         const uint currSizePower = unionOfBlocksSet->getSizePower();
         if (currSizePower < HashmapReqSize) {
            unionOfBlocksSet->resize(HashmapReqSize);
            unionOfBlocksSet->optimizeGPU(stream);
         }
         // Ensure map is empty
         unionOfBlocksSet->clear(Hashinator::targets::device,stream,false);
      }
      gpu_allocated_largestVmesh = largestVmesh;
   }
   // Vector into which the set contents are read (prefetched to device)
   if (unionSetSize != 0) {
      if (gpu_allocated_unionSetSize == 0) {
         // New allocation
         unionOfBlocks = new split::SplitVector<vmesh::GlobalID>(unionSetSize);
         unionOfBlocks->clear();
         unionOfBlocks->optimizeGPU(stream);
      } else {
         if (unionOfBlocks->capacity() < unionSetSize) {
            // Recapacitate, clear, and prefetch
            unionOfBlocks->reserve(unionSetSize);
            unionOfBlocks->clear();
            unionOfBlocks->optimizeGPU(stream);
         } else {
            // Clear is enough
            unionOfBlocks->clear();
         }
      }
      gpu_allocated_unionSetSize = unionSetSize;
   }
   CHK_ERR( gpuStreamSynchronize(stream) );
}

/* Deallocation at end of simulation */
__host__ void gpu_trans_deallocate() {
   // Deallocate any translation vectors or sets which exist
   if (gpu_allocated_nAllCells != 0) {
      delete allVmeshPointer;
   }
   if (gpu_allocated_sumOfLengths != 0) {
      delete allPencilsMeshes;
      delete allPencilsContainers;
   }
   if (gpu_allocated_largestVmesh != 0) {
      delete unionOfBlocksSet;
   }
   if (gpu_allocated_unionSetSize != 0) {
      delete unionOfBlocks;
   }
   gpu_allocated_nAllCells = 0;
   gpu_allocated_sumOfLengths = 0;
   gpu_allocated_largestVmesh = 0;
   gpu_allocated_unionSetSize = 0;
   // Delete also the vectors for pencils for each dimension
   for (uint dimension=0; dimension<3; dimension++) {
      if (DimensionPencils[dimension].gpu_allocated) {
         delete DimensionPencils[dimension].gpu_lengthOfPencils;
         delete DimensionPencils[dimension].gpu_idsStart;
         delete DimensionPencils[dimension].gpu_sourceDZ;
         delete DimensionPencils[dimension].gpu_targetRatios;
      }
   }
}
