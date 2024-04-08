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

#include "gpu_base.hpp"
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

bool needAttachedStreams = false;
bool doPrefetches=false; // only non-crucial prefetches are behind this check
// Note: disabling prefetches brings in strange memory errors and crashes (June 2023)

uint *gpu_cell_indices_to_id[MAXCPUTHREADS];
uint *gpu_block_indices_to_id[MAXCPUTHREADS];
uint *gpu_vcell_transpose; // only one needed, not one per thread

ColumnOffsets *unif_columnOffsetData[MAXCPUTHREADS];
Column *gpu_columns[MAXCPUTHREADS];

Vec *gpu_blockDataOrdered[MAXCPUTHREADS];

vmesh::LocalID *gpu_GIDlist[MAXCPUTHREADS];
vmesh::LocalID *gpu_LIDlist[MAXCPUTHREADS];
vmesh::GlobalID *gpu_BlocksID_mapped[MAXCPUTHREADS];
vmesh::GlobalID *gpu_BlocksID_mapped_sorted[MAXCPUTHREADS];
vmesh::GlobalID *gpu_LIDlist_unsorted[MAXCPUTHREADS];
vmesh::LocalID *gpu_columnNBlocks[MAXCPUTHREADS];
void *gpu_RadixSortTemp[MAXCPUTHREADS];
uint gpu_acc_RadixSortTempSize[MAXCPUTHREADS] = {0};

// Vectors and set for use in translation
split::SplitVector<vmesh::VelocityMesh*> *allVmeshPointer;
split::SplitVector<vmesh::VelocityMesh*> *allPencilsMeshes;
split::SplitVector<vmesh::VelocityBlockContainer*> *allPencilsContainers;
split::SplitVector<vmesh::GlobalID> *unionOfBlocks;
Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet;
// counters for allocated sizes
uint gpu_allocated_nAllCells = 0;
uint gpu_allocated_sumOfLengths = 0;
uint gpu_allocated_largestVmesh = 0;
uint gpu_allocated_unionSetSize = 0;

// Memory allocation flags and values.
uint gpu_vlasov_allocatedSize = 0;
uint gpu_acc_allocatedColumns = 0;
uint gpu_acc_columnContainerSize = 0;
uint gpu_acc_foundColumnsCount = 0;

// SplitVectors and buffers for use in stream compaction
uint gpu_compaction_vectorsize[MAXCPUTHREADS] = {0};
uint gpu_compaction_buffersize[MAXCPUTHREADS] = {0};
split::SplitVector<vmesh::GlobalID> *vbwcl_gather[MAXCPUTHREADS];
split::SplitVector<vmesh::GlobalID> *vbwncl_gather[MAXCPUTHREADS];
void *compaction_buffer[MAXCPUTHREADS];

__host__ void gpu_init_device() {

#ifdef _OPENMP
   const uint maxThreads = omp_get_max_threads();
#else
   const uint maxThreads = 1;
#endif

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
      needAttachedStreams = true;
   }
   #endif

   // Init Umpire memory manager
   #if defined(USE_UMPIRE)
     auto& rm = umpire::ResourceManager::getInstance();
     umpire::Allocator allocator = rm.getAllocator("DEVICE");
     auto pooled_allocator = rm.makeAllocator<umpire::strategy::QuickPool>("DEV_POOL", allocator, 1024, 1024);
   #endif

   // Pre-generate streams, allocate return pointers
   int *leastPriority = new int; // likely 0
   int *greatestPriority = new int; // likely -1
   CHK_ERR( gpuDeviceGetStreamPriorityRange (leastPriority, greatestPriority) );
   if (*leastPriority==*greatestPriority) {
      printf("Warning when initializing GPU streams: minimum and maximum stream priority are identical! %d == %d \n",*leastPriority, *greatestPriority);
   }
   for (uint i=0; i<maxThreads; ++i) {
      CHK_ERR( gpuStreamCreateWithPriority(&(gpuStreamList[i]), gpuStreamDefault, *leastPriority) );
      CHK_ERR( gpuStreamCreateWithPriority(&(gpuPriorityStreamList[i]), gpuStreamDefault, *greatestPriority) );
      CHK_ERR( gpuMalloc((void**)&returnReal[i], 8*sizeof(Real)) );
      CHK_ERR( gpuMalloc((void**)&returnRealf[i], 8*sizeof(Realf)) );
      CHK_ERR( gpuMalloc((void**)&returnLID[i], 8*sizeof(vmesh::LocalID)) );
   }
   CHK_ERR( gpuDeviceSynchronize() );

   // Using just a single context for whole MPI task
}

__host__ void gpu_clear_device() {
   // Deallocate temporary buffers
   gpu_acc_deallocate();
   gpu_vlasov_deallocate();
   gpu_trans_deallocate();
   gpu_compaction_deallocate();
   // Destroy streams
#ifdef _OPENMP
   const uint maxThreads = omp_get_max_threads();
#else
   const uint maxThreads = 1;
#endif
   for (uint i=0; i<maxThreads; ++i) {
      CHK_ERR( gpuStreamDestroy(gpuStreamList[i]) );
      CHK_ERR( gpuStreamDestroy(gpuPriorityStreamList[i]) );
      CHK_ERR( gpuFree(returnReal[i]) );
      CHK_ERR( gpuFree(returnRealf[i]) );
      CHK_ERR( gpuFree(returnLID[i]) );
   }
   CHK_ERR( gpuDeviceSynchronize() );
}

__host__ gpuStream_t gpu_getStream() {
#ifdef _OPENMP
   const uint thread_id = omp_get_thread_num();
#else
   const uint thread_id = 0;
#endif
   return gpuStreamList[thread_id];
}

__host__ gpuStream_t gpu_getPriorityStream() {
#ifdef _OPENMP
   const uint thread_id = omp_get_thread_num();
#else
   const uint thread_id = 0;
#endif
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
   // Always prepare for at least 2500 blocks (affects also translation parallelism)
   const uint maxBlocksPerCell = maxBlockCount > 2500 ? maxBlockCount : 2500;
   // Check if we already have allocated enough memory?
   if (gpu_vlasov_allocatedSize > maxBlocksPerCell * BLOCK_ALLOCATION_FACTOR) {
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
      if (gpu_vlasov_allocatedSize > 0) {
         gpu_vlasov_deallocate_perthread(i);
         CHK_ERR( gpuFree(gpu_vcell_transpose) );
      }
      gpu_vlasov_allocate_perthread(i, newSize);
      CHK_ERR( gpuMalloc((void**)&gpu_vcell_transpose, WID3*sizeof(uint)) );
   }
   gpu_vlasov_allocatedSize = newSize;
}

/* Deallocation at end of simulation */
__host__ void gpu_vlasov_deallocate() {
   if (gpu_vlasov_allocatedSize == 0) return;
#ifdef _OPENMP
   const uint maxNThreads = omp_get_max_threads();
#else
   const uint maxNThreads = 1;
#endif
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_vlasov_deallocate_perthread(i);
   }
}

__host__ uint gpu_vlasov_getAllocation() {
   return gpu_vlasov_allocatedSize;
}

__host__ void gpu_vlasov_allocate_perthread(
   uint cpuThreadID,
   uint blockAllocationCount
   ) {
   // Mallocs should be in increments of 256 bytes. WID3 is at least 64, and len(Realf) is at least 4, so this is true.
   CHK_ERR( gpuMalloc((void**)&gpu_cell_indices_to_id[cpuThreadID], 3*sizeof(uint)) );
   CHK_ERR( gpuMalloc((void**)&gpu_block_indices_to_id[cpuThreadID], 3*sizeof(uint)) );
   CHK_ERR( gpuMalloc((void**)&gpu_blockDataOrdered[cpuThreadID], blockAllocationCount * (WID3 / VECL) * sizeof(Vec)) );
   CHK_ERR( gpuMalloc((void**)&gpu_BlocksID_mapped[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );
   CHK_ERR( gpuMalloc((void**)&gpu_BlocksID_mapped_sorted[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );
   CHK_ERR( gpuMalloc((void**)&gpu_LIDlist_unsorted[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );
   CHK_ERR( gpuMalloc((void**)&gpu_LIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::LocalID)) );
   CHK_ERR( gpuMalloc((void**)&gpu_GIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::LocalID)) );
}

__host__ void gpu_vlasov_deallocate_perthread (
   uint cpuThreadID
   ) {
   gpu_vlasov_allocatedSize = 0;
   CHK_ERR( gpuFree(gpu_cell_indices_to_id[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_block_indices_to_id[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_blockDataOrdered[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_BlocksID_mapped[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_BlocksID_mapped_sorted[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_LIDlist_unsorted[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_LIDlist[cpuThreadID]) );
   CHK_ERR( gpuFree(gpu_GIDlist[cpuThreadID]) );
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
#ifdef _OPENMP
   const uint maxNThreads = omp_get_max_threads();
#else
   const uint maxNThreads = 1;
#endif
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
#ifdef _OPENMP
   const uint maxNThreads = omp_get_max_threads();
#else
   const uint maxNThreads = 1;
#endif
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_acc_deallocate_perthread(i);
   }
}

__host__ void gpu_acc_allocate_perthread(
   uint cpuThreadID,
   uint columnAllocationCount
   ) {
   // Unified memory; columndata contains several splitvectors.
   if (columnAllocationCount > 0) {
      unif_columnOffsetData[cpuThreadID] = new ColumnOffsets(columnAllocationCount); // inherits managed
      unif_columnOffsetData[cpuThreadID]->gpu_advise();
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
      delete unif_columnOffsetData[cpuThreadID];
   }
   CHK_ERR( gpuFree(gpu_columnNBlocks[cpuThreadID]) );
}

/*
   Top-level GPU memory allocation function for stream compaction buffers.
 */
__host__ void gpu_compaction_allocate(
   const uint vectorLength,
   const size_t bytesNeeded
   ) {
#ifdef _OPENMP
   const uint maxNThreads = omp_get_max_threads();
#else
   const uint maxNThreads = 1;
#endif
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_compaction_allocate_vec_perthread(i,vectorLength);
      gpu_compaction_allocate_buf_perthread(i,bytesNeeded);
   }
}

/* Deallocation at end of simulation */
__host__ void gpu_compaction_deallocate() {
#ifdef _OPENMP
   const uint maxNThreads = omp_get_max_threads();
#else
   const uint maxNThreads = 1;
#endif
   for (uint i=0; i<maxNThreads; ++i) {
      gpu_compaction_allocate_vec_perthread(i,0);
      gpu_compaction_allocate_buf_perthread(i,0);
   }
}

__host__ void gpu_compaction_allocate_vec_perthread(
   const uint cpuThreadID,
   const uint vectorLength
   ) {
   gpuStream_t stream = gpu_getStream();
   // Deallocate if: buffer exists, and: requesting zero size or size is not large enough.
   if ( (gpu_compaction_vectorsize[cpuThreadID] > 0)
        && ( (vectorLength == 0) || vectorLength > gpu_compaction_vectorsize[cpuThreadID])
      ) {
      delete vbwcl_gather[cpuThreadID];
      delete vbwncl_gather[cpuThreadID];
      gpu_compaction_vectorsize[cpuThreadID] = 0;
   }
   // If buffer isn't large enough, allocate.
   if (vectorLength > gpu_compaction_vectorsize[cpuThreadID]) {
      const uint paddedSize = BLOCK_ALLOCATION_FACTOR * vectorLength;
      vbwcl_gather[cpuThreadID] = new split::SplitVector<vmesh::GlobalID>(paddedSize);
      vbwncl_gather[cpuThreadID] = new split::SplitVector<vmesh::GlobalID>(paddedSize);
      // int device = gpu_getDevice();
      // vbwcl_gather[cpuThreadID]->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      // vbwncl_gather[cpuThreadID]->memAdvise(gpuMemAdviseSetPreferredLocation,device,stream);
      // vbwcl_gather[cpuThreadID]->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      // vbwncl_gather[cpuThreadID]->memAdvise(gpuMemAdviseSetAccessedBy,device,stream);
      vbwcl_gather[cpuThreadID]->optimizeGPU(stream);
      vbwncl_gather[cpuThreadID]->optimizeGPU(stream);
      gpu_compaction_vectorsize[cpuThreadID] = paddedSize;
   }
}
__host__ void gpu_compaction_allocate_buf_perthread(
   const uint cpuThreadID,
   const size_t bytesNeeded
   ) {
   gpuStream_t stream = gpu_getStream();
   // Deallocate if: buffer exists, and: requesting zero size or size is not large enough.
   if ( (gpu_compaction_buffersize[cpuThreadID] > 0)
        && ( (bytesNeeded == 0) || bytesNeeded > gpu_compaction_buffersize[cpuThreadID])
      ) {
      CHK_ERR( gpuFreeAsync(compaction_buffer[cpuThreadID],stream) );
      gpu_compaction_buffersize[cpuThreadID] = 0;
   }
   // If buffer isn't large enough, allocate.
   if (bytesNeeded > gpu_compaction_buffersize[cpuThreadID]) {
      const uint paddedSize = BLOCK_ALLOCATION_FACTOR * bytesNeeded;
      CHK_ERR( gpuMallocAsync((void**)&compaction_buffer[cpuThreadID],paddedSize, stream) );
      gpu_compaction_buffersize[cpuThreadID] = paddedSize;
   }
   CHK_ERR( gpuStreamSynchronize(stream) );
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
         allVmeshPointer->resize(nAllCells,true);
         allVmeshPointer->optimizeCPU(stream);
      }
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
         allPencilsMeshes->resize(sumOfLengths,true);
         allPencilsContainers->resize(sumOfLengths,true);
         allPencilsMeshes->optimizeCPU(stream);
         allPencilsContainers->optimizeCPU(stream);
      }
      gpu_allocated_sumOfLengths = sumOfLengths;
   }
   // Set for collecting union of blocks (prefetched to device)
   if (largestVmesh != 0) {
      const vmesh::LocalID HashmapReqSize = ceil(log2(largestVmesh)) +3;
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
   CHK_ERR( gpuFree(gpu_vcell_transpose) );
}
