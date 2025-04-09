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

#include <stdio.h>
#include <iostream>
#include "common.h"
#include "mpi.h"

#include "gpu_base.hpp"
#include "velocity_mesh_parameters.h"
#include "../vlasovsolver/gpu_moments.h"
#include "../spatial_cells/velocity_mesh_gpu.h"
#include "../spatial_cells/velocity_block_container.h"
#include "../vlasovsolver/cpu_trans_pencils.hpp"
#include "../vlasovsolver/gpu_moments.h"
#include "../vlasovsolver/gpu_acc_sort_blocks.hpp"

#include "logger.h"

// #define MAXCPUTHREADS 64 now in gpu_base.hpp

extern Logger logFile;
int myDevice;
int myRank;

// Allocate pointers for per-thread memory regions
gpuStream_t gpuStreamList[MAXCPUTHREADS];
gpuStream_t gpuPriorityStreamList[MAXCPUTHREADS];

Real *returnReal[MAXCPUTHREADS];
Realf *returnRealf[MAXCPUTHREADS];
vmesh::LocalID *returnLID[MAXCPUTHREADS];
Real *host_returnReal[MAXCPUTHREADS];
Realf *host_returnRealf[MAXCPUTHREADS];
vmesh::LocalID *host_returnLID[MAXCPUTHREADS];

uint *gpu_cell_indices_to_id[MAXCPUTHREADS];
uint *gpu_block_indices_to_id[MAXCPUTHREADS];
uint *gpu_vcell_transpose; // only one needed, not one per thread

// Pointers to buffers used in acceleration
ColumnOffsets *cpu_columnOffsetData[MAXCPUTHREADS] = {0};
ColumnOffsets *gpu_columnOffsetData[MAXCPUTHREADS] = {0};
Column *gpu_columns[MAXCPUTHREADS];
Vec *gpu_blockDataOrdered[MAXCPUTHREADS] = {0};
vmesh::GlobalID *gpu_GIDlist[MAXCPUTHREADS];
vmesh::LocalID *gpu_LIDlist[MAXCPUTHREADS];
vmesh::GlobalID *gpu_BlocksID_mapped[MAXCPUTHREADS];
vmesh::GlobalID *gpu_BlocksID_mapped_sorted[MAXCPUTHREADS];
vmesh::LocalID *gpu_LIDlist_unsorted[MAXCPUTHREADS];
vmesh::LocalID *gpu_columnNBlocks[MAXCPUTHREADS] = {0};
vmesh::GlobalID *invalidGIDpointer = 0;
void *gpu_RadixSortTemp[MAXCPUTHREADS] = {0};
size_t gpu_acc_RadixSortTempSize[MAXCPUTHREADS] = {0};

// Hash map and splitvectors buffers used in block adjustment
vmesh::VelocityMesh** host_vmeshes, **dev_vmeshes;
vmesh::VelocityBlockContainer** host_VBCs, **dev_VBCs;
Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** host_allMaps, **dev_allMaps;
split::SplitVector<vmesh::GlobalID> ** host_vbwcl_vec, **dev_vbwcl_vec;
split::SplitVector<vmesh::GlobalID> ** host_lists_with_replace_new, **dev_lists_with_replace_new;
split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_delete, **dev_lists_delete;
split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_to_replace, **dev_lists_to_replace;
split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_with_replace_old, **dev_lists_with_replace_old;
split::SplitVector<vmesh::GlobalID> ** host_vbwcl_neigh, **dev_vbwcl_neigh;
vmesh::LocalID* host_contentSizes, *dev_contentSizes;
Real* host_minValues, *dev_minValues;
Real* host_massLoss, *dev_massLoss;
Real* host_mass, *dev_mass;

// Vectors and set for use in translation (and in vlasovsolver/gpu_dt.cpp)
split::SplitVector<vmesh::VelocityMesh*> *allVmeshPointer=0, *dev_allVmeshPointer=0;
split::SplitVector<vmesh::VelocityMesh*> *allPencilsMeshes=0, *dev_allPencilsMeshes=0;
split::SplitVector<vmesh::VelocityBlockContainer*> *allPencilsContainers=0, *dev_allPencilsContainers=0;
split::SplitVector<vmesh::GlobalID> *unionOfBlocks=0, *dev_unionOfBlocks=0;
Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet=0, *dev_unionOfBlocksSet=0;

// pointers for translation
Vec** host_pencilOrderedPointers;
Vec** dev_pencilOrderedPointers;
Realf** dev_pencilBlockData; // Array of pointers into actual block data
uint* dev_pencilBlocksCount; // Array of counters if pencil needs to be propagated for this block or not


// counters for allocated sizes in translation
uint gpu_allocated_nAllCells = 0;
uint gpu_allocated_sumOfLengths = 0;
uint gpu_allocated_largestVmeshSizePower = 0;
uint gpu_allocated_unionSetSize = 0;
uint gpu_allocated_trans_pencilBlockData = 0;
uint gpu_allocated_trans_pencilBlocksCount = 0;

// batch counters
uint gpu_allocated_batch_nCells = 0;
uint gpu_allocated_batch_maxNeighbours = 0;
// Memory allocation flags and values (TODO make per-thread?).
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
   //std::cerr << "(Grid) rank " << amps_rank << " is noderank "<< amps_node_rank << " of "<< amps_node_size << std::endl;
   myRank = amps_rank;

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
      printf("Error! Current GPU device does not support concurrent managed memory access from several streams.\n");
      printf("Please switch to a more recent CUDA compute architecture.\n");
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
      CHK_ERR( gpuMallocHost((void **) &host_returnReal[i], 8*sizeof(Real)) );
      CHK_ERR( gpuMallocHost((void **) &host_returnRealf[i], 8*sizeof(Realf)) );
      CHK_ERR( gpuMallocHost((void **) &host_returnLID[i], 8*sizeof(vmesh::LocalID)) );
      // CHK_ERR( gpuMallocHost((void **) &info_1[i], sizeof(split::SplitInfo)) );
      // CHK_ERR( gpuMallocHost((void **) &info_2[i], sizeof(split::SplitInfo)) );
      // CHK_ERR( gpuMallocHost((void **) &info_3[i], sizeof(split::SplitInfo)) );
      // CHK_ERR( gpuMallocHost((void **) &info_4[i], sizeof(split::SplitInfo)) );
      // CHK_ERR( gpuMallocHost((void **) &info_m[i], sizeof(Hashinator::MapInfo)) );
   }
   CHK_ERR( gpuMalloc((void**)&dev_pencilOrderedPointers, maxNThreads*sizeof(Vec*)) );
   CHK_ERR( gpuMallocHost((void **)&host_pencilOrderedPointers, maxNThreads*sizeof(Vec*)) );

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
   gpu_batch_deallocate();
   // Destroy streams
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      CHK_ERR( gpuStreamDestroy(gpuStreamList[i]) );
      CHK_ERR( gpuStreamDestroy(gpuPriorityStreamList[i]) );
      CHK_ERR( gpuFree(returnReal[i]) );
      CHK_ERR( gpuFree(returnRealf[i]) );
      CHK_ERR( gpuFree(returnLID[i]) );
      if (gpu_RadixSortTemp[i]) {
         CHK_ERR( gpuFree(gpu_RadixSortTemp[i]) );
         gpu_RadixSortTemp[i] = 0;
      }
      CHK_ERR( gpuFreeHost(host_returnReal[i]) );
      CHK_ERR( gpuFreeHost(host_returnRealf[i]) );
      CHK_ERR( gpuFreeHost(host_returnLID[i]) );
   }
   CHK_ERR( gpuFree(dev_pencilOrderedPointers) );
   CHK_ERR( gpuFreeHost(host_pencilOrderedPointers) );
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


/* Memory reporting function
 */
int gpu_reportMemory(const size_t local_cells_capacity, const size_t ghost_cells_capacity,
                     const size_t local_cells_size, const size_t ghost_cells_size) {
   /* Gather total CPU and GPU buffer sizes. Rank 0 reports details,
      all ranks return sum.
   */
   uint maxNThreads = gpu_getMaxThreads();

   size_t miniBuffers = maxNThreads * (
      8*sizeof(Real) // returnReal
      + 8*sizeof(Realf) // returnRealf
      + 8*sizeof(vmesh::LocalID) // returnLID
      + sizeof(Vec*) // dev_pencilOrderedPointers
      )
      + WID3*sizeof(uint) // gpu_vcell_transpose
      + sizeof(vmesh::GlobalID) // invalidGIDpointer
      + sizeof(std::array<vmesh::MeshParameters,MAX_VMESH_PARAMETERS_COUNT>) // velocityMeshes_upload
      + sizeof(vmesh::MeshWrapper) // MWdev
      + gpu_allocated_moments*sizeof(vmesh::VelocityBlockContainer*) // gpu_moments dev_VBC
      + gpu_allocated_moments*4*sizeof(Real)  // gpu_moments dev_moments1
      + gpu_allocated_moments*3*sizeof(Real); // gpu_moments dev_moments2
   // DT reduction buffers are deallocated every step (GPUTODO)

   size_t vlasovBuffers = 0;
   for (uint i=0; i<maxNThreads; ++i) {
      vlasovBuffers += 6*sizeof(uint) // gpu_cell_indices_to_id[cpuThreadID], gpu_block_indices_to_id[cpuThreadID]
         + gpu_vlasov_allocatedSize[i] * (
            TRANSLATION_BUFFER_ALLOCATION_FACTOR * (WID3 / VECL) * sizeof(Vec) // gpu_blockDataOrdered[cpuThreadID]
            + 3*sizeof(vmesh::GlobalID) // gpu_BlocksID_mapped, gpu_BlocksID_mapped_sorted, gpu_GIDlist
            + 2*sizeof(vmesh::LocalID) ); // gpu_LIDlist_unsorted, gpu_LIDlist
   }

   size_t batchBuffers = gpu_allocated_batch_nCells * (
      sizeof(vmesh::VelocityMesh*) // dev_vmeshes
      + sizeof(vmesh::VelocityBlockContainer*) // dev_VBCs
      + 2 * sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*) // dev_allMaps
      + 2 * sizeof(split::SplitVector<vmesh::GlobalID>*) // dev_vbwcl_vec, dev_lists_with_replace_new
      + 3 * sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*) // dev_lists_delete, dev_lists_to_replace, dev_lists_with_replace_old
      + 5 * sizeof(vmesh::LocalID) // dev_contentSizes
      + 3 * sizeof(Real) // dev_minValue, dev_massLoss, dev_mass
      );
   batchBuffers += gpu_allocated_batch_maxNeighbours * sizeof(split::SplitVector<vmesh::GlobalID>*); // dev_vbwcl_neigh

   size_t accBuffers = 0;
   for (uint i=0; i<maxNThreads; ++i) {
      accBuffers += gpu_acc_allocatedColumns * sizeof(Column); // gpu_columns[cpuThreadID]
      accBuffers += sizeof(ColumnOffsets); // gpu_columnOffsetData[cpuThreadID]
      if (cpu_columnOffsetData[i]) {
         accBuffers += cpu_columnOffsetData[i]->capacityInBytes(); // struct contents
      }
      accBuffers += gpu_acc_columnContainerSize*sizeof(vmesh::LocalID); // column id counts, maximal vmesh
      accBuffers += gpu_acc_RadixSortTempSize[i]; // gpu_RadixSortTemp[cpuThreadID]
   }

   size_t transBuffers = 0;
   if (allVmeshPointer) {
      transBuffers += sizeof(split::SplitVector<vmesh::VelocityMesh*>);
      transBuffers += allVmeshPointer->capacity() * sizeof(vmesh::VelocityMesh*);
   }
   if (allPencilsMeshes) {
      transBuffers += sizeof(split::SplitVector<vmesh::VelocityMesh*>);
      transBuffers += allPencilsMeshes->capacity() * sizeof(vmesh::VelocityMesh*);
   }
   if (allPencilsContainers) {
      transBuffers += sizeof(split::SplitVector<vmesh::VelocityBlockContainer*>);
      transBuffers += allPencilsContainers->capacity() * sizeof(vmesh::VelocityBlockContainer*);
   }
   if (unionOfBlocksSet) {
      transBuffers += sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>);
      transBuffers += unionOfBlocksSet->bucket_count() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
   }
   if (unionOfBlocks) {
      transBuffers += sizeof(split::SplitVector<vmesh::GlobalID>);
      transBuffers += unionOfBlocks->capacity() * sizeof(vmesh::GlobalID);
   }
   transBuffers += gpu_allocated_trans_pencilBlockData * sizeof(Realf*); //dev_pencilBlockData
   transBuffers += gpu_allocated_trans_pencilBlocksCount * sizeof(uint); // dev_pencilBlocksCount
   // Pencils:
   for (uint dimension=0; dimension<3; ++dimension) {
      transBuffers += DimensionPencils[dimension].gpu_allocated_sumOfLengths * 2 * sizeof(Realf); // gpu_lengthOfPencils, gpu_idsStart
      transBuffers += DimensionPencils[dimension].gpu_allocated_N * 2 * sizeof(uint); // gpu_sourceDZ, gpu_targetRatios
   }
   // Remote neighbor contribution buffers are in unified memory but deallocated after each use

   size_t free_byte ;
   size_t total_byte ;
   CHK_ERR( gpuMemGetInfo( &free_byte, &total_byte) );
   size_t used_mb = (total_byte-free_byte)/(1024*1024);
   size_t sum_mb = (miniBuffers+batchBuffers+vlasovBuffers+accBuffers+transBuffers+local_cells_capacity+ghost_cells_capacity)/(1024*1024);
   size_t local_req_mb = local_cells_size/(1024*1024);
   size_t ghost_req_mb = ghost_cells_size/(1024*1024);
   if (myRank==0) {
      logFile<<" =================================="<<std::endl;
      logFile<<" GPU Memory report"<<std::endl;
      logFile<<"     mini-buffers:          "<<miniBuffers/(1024*1024)<<" Mbytes"<<std::endl;
      logFile<<"     Batch buffers:         "<<batchBuffers/(1024*1024)<<" Mbytes"<<std::endl;
      logFile<<"     Vlasov buffers:        "<<vlasovBuffers/(1024*1024)<<" Mbytes"<<std::endl;
      logFile<<"     Acceleration buffers:  "<<accBuffers/(1024*1024)<<" Mbytes"<<std::endl;
      logFile<<"     Translation buffers:   "<<transBuffers/(1024*1024)<<" Mbytes"<<std::endl;
      logFile<<"     Local cells:           "<<local_cells_capacity/(1024*1024)<<" Mbytes"<<std::endl;
      logFile<<"     Ghost cells:           "<<ghost_cells_capacity/(1024*1024)<<" Mbytes"<<std::endl;
      if (local_req_mb || ghost_req_mb) {
         logFile<<"     Local cells required:  "<<local_req_mb<<" Mbytes"<<std::endl;
         logFile<<"     Ghost cells required:  "<<ghost_req_mb<<" Mbytes"<<std::endl;
      }
      logFile<<"   Total:                   "<<sum_mb<<" Mbytes"<<std::endl;
      logFile<<"   Reported Hardware use:   "<<used_mb<<" Mbytes"<<std::endl;
      logFile<<" =================================="<<std::endl;
   }
   return sum_mb;
}

/*
   Top-level GPU memory allocation function.
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void gpu_vlasov_allocate(
   uint maxBlockCount // Largest found vmesh size
   ) {
   // Always prepare for at least VLASOV_BUFFER_MINBLOCKS blocks
   const uint maxBlocksPerCell = max(VLASOV_BUFFER_MINBLOCKS, maxBlockCount);
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
   gpuStream_t stream = gpu_getStream();

   // Mallocs should be in increments of 256 bytes. WID3 is at least 64, and len(Realf) is at least 4, so this is true.
   CHK_ERR( gpuMallocAsync((void**)&gpu_cell_indices_to_id[cpuThreadID], 3*sizeof(uint), stream) );
   CHK_ERR( gpuMallocAsync((void**)&gpu_block_indices_to_id[cpuThreadID], 3*sizeof(uint), stream) );
   CHK_ERR( gpuMallocAsync((void**)&gpu_blockDataOrdered[cpuThreadID], newSize * TRANSLATION_BUFFER_ALLOCATION_FACTOR * (WID3 / VECL) * sizeof(Vec), stream) );
   CHK_ERR( gpuMallocAsync((void**)&gpu_BlocksID_mapped[cpuThreadID], newSize*sizeof(vmesh::GlobalID), stream) );
   CHK_ERR( gpuMallocAsync((void**)&gpu_BlocksID_mapped_sorted[cpuThreadID], newSize*sizeof(vmesh::GlobalID), stream) );
   CHK_ERR( gpuMallocAsync((void**)&gpu_LIDlist_unsorted[cpuThreadID], newSize*sizeof(vmesh::LocalID), stream) );
   CHK_ERR( gpuMallocAsync((void**)&gpu_LIDlist[cpuThreadID], newSize*sizeof(vmesh::LocalID), stream) );
   CHK_ERR( gpuMallocAsync((void**)&gpu_GIDlist[cpuThreadID], newSize*sizeof(vmesh::GlobalID), stream) );
   // Store size of new allocation
   gpu_vlasov_allocatedSize[cpuThreadID] = newSize;
}

__host__ void gpu_vlasov_deallocate_perthread (
   uint cpuThreadID
   ) {
   if (gpu_vlasov_allocatedSize[cpuThreadID] == 0) {
      return;
   }
   gpuStream_t stream = gpu_getStream();
   CHK_ERR( gpuFreeAsync(gpu_cell_indices_to_id[cpuThreadID],stream) );
   CHK_ERR( gpuFreeAsync(gpu_block_indices_to_id[cpuThreadID],stream) );
   CHK_ERR( gpuFreeAsync(gpu_blockDataOrdered[cpuThreadID],stream) );
   CHK_ERR( gpuFreeAsync(gpu_BlocksID_mapped[cpuThreadID],stream) );
   CHK_ERR( gpuFreeAsync(gpu_BlocksID_mapped_sorted[cpuThreadID],stream) );
   CHK_ERR( gpuFreeAsync(gpu_LIDlist_unsorted[cpuThreadID],stream) );
   CHK_ERR( gpuFreeAsync(gpu_LIDlist[cpuThreadID],stream) );
   CHK_ERR( gpuFreeAsync(gpu_GIDlist[cpuThreadID],stream) );
   gpu_vlasov_allocatedSize[cpuThreadID] = 0;
}

/** Allocation and deallocation for pointers used by batch operations in block adjustment */
__host__ void gpu_batch_allocate(uint nCells, uint maxNeighbours) {
   if (nCells > gpu_allocated_batch_nCells) {
      gpu_batch_deallocate(true, false);
      gpu_allocated_batch_nCells = nCells * BLOCK_ALLOCATION_FACTOR;

      CHK_ERR( gpuMallocHost((void**)&host_vmeshes,gpu_allocated_batch_nCells*sizeof(vmesh::VelocityMesh*)) );
      CHK_ERR( gpuMallocHost((void**)&host_VBCs,gpu_allocated_batch_nCells*sizeof(vmesh::VelocityBlockContainer*)) );
      CHK_ERR( gpuMallocHost((void**)&host_allMaps, 2*gpu_allocated_batch_nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*)) ); // note double size
      CHK_ERR( gpuMallocHost((void**)&host_vbwcl_vec, gpu_allocated_batch_nCells*sizeof(split::SplitVector<vmesh::GlobalID>*)) );
      CHK_ERR( gpuMallocHost((void**)&host_lists_with_replace_new, gpu_allocated_batch_nCells*sizeof(split::SplitVector<vmesh::GlobalID>*)) );
      CHK_ERR( gpuMallocHost((void**)&host_lists_delete, gpu_allocated_batch_nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*)) );
      CHK_ERR( gpuMallocHost((void**)&host_lists_to_replace, gpu_allocated_batch_nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*)) );
      CHK_ERR( gpuMallocHost((void**)&host_lists_with_replace_old, gpu_allocated_batch_nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*)) );
      CHK_ERR( gpuMallocHost((void**)&host_contentSizes,gpu_allocated_batch_nCells*5*sizeof(vmesh::LocalID)) ); // note quadruple size
      CHK_ERR( gpuMallocHost((void**)&host_minValues, gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMallocHost((void**)&host_massLoss, gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMallocHost((void**)&host_mass, gpu_allocated_batch_nCells*sizeof(Real)) );

      CHK_ERR( gpuMalloc((void**)&dev_vmeshes,gpu_allocated_batch_nCells*sizeof(vmesh::VelocityMesh*)) );
      CHK_ERR( gpuMalloc((void**)&dev_VBCs,gpu_allocated_batch_nCells*sizeof(vmesh::VelocityBlockContainer*)) );
      CHK_ERR( gpuMalloc((void**)&dev_allMaps, 2*gpu_allocated_batch_nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_vbwcl_vec, gpu_allocated_batch_nCells*sizeof(split::SplitVector<vmesh::GlobalID>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_lists_with_replace_new, gpu_allocated_batch_nCells*sizeof(split::SplitVector<vmesh::GlobalID>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_lists_delete, gpu_allocated_batch_nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_lists_to_replace, gpu_allocated_batch_nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_lists_with_replace_old, gpu_allocated_batch_nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_contentSizes,gpu_allocated_batch_nCells*5*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMalloc((void**)&dev_minValues,gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMalloc((void**)&dev_massLoss, gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMalloc((void**)&dev_mass, gpu_allocated_batch_nCells*sizeof(Real)) );
   }

   if (maxNeighbours*gpu_allocated_batch_nCells > gpu_allocated_batch_maxNeighbours) {
      gpu_batch_deallocate(false, true);
      gpu_allocated_batch_maxNeighbours = maxNeighbours * gpu_allocated_batch_nCells * BLOCK_ALLOCATION_FACTOR;
      CHK_ERR( gpuMallocHost((void**)&host_vbwcl_neigh, gpu_allocated_batch_maxNeighbours*sizeof(split::SplitVector<vmesh::GlobalID>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_vbwcl_neigh, gpu_allocated_batch_maxNeighbours*sizeof(split::SplitVector<vmesh::GlobalID>*)) );
   }
}
__host__ void gpu_batch_deallocate(bool first, bool second) {
   if (gpu_allocated_batch_nCells != 0 && first) {
      gpu_allocated_batch_nCells = 0;
      CHK_ERR( gpuFreeHost(host_vmeshes));
      CHK_ERR( gpuFreeHost(host_VBCs));
      CHK_ERR( gpuFreeHost(host_allMaps));
      CHK_ERR( gpuFreeHost(host_vbwcl_vec));
      CHK_ERR( gpuFreeHost(host_lists_with_replace_new));
      CHK_ERR( gpuFreeHost(host_lists_delete));
      CHK_ERR( gpuFreeHost(host_lists_to_replace));
      CHK_ERR( gpuFreeHost(host_lists_with_replace_old));
      CHK_ERR( gpuFreeHost(host_contentSizes));
      CHK_ERR( gpuFreeHost(host_minValues));
      CHK_ERR( gpuFreeHost(host_massLoss));
      CHK_ERR( gpuFreeHost(host_mass));
      CHK_ERR( gpuFree(dev_vmeshes));
      CHK_ERR( gpuFree(dev_VBCs));
      CHK_ERR( gpuFree(dev_allMaps));
      CHK_ERR( gpuFree(dev_vbwcl_vec));
      CHK_ERR( gpuFree(dev_lists_with_replace_new));
      CHK_ERR( gpuFree(dev_lists_delete));
      CHK_ERR( gpuFree(dev_lists_to_replace));
      CHK_ERR( gpuFree(dev_lists_with_replace_old));
      CHK_ERR( gpuFree(dev_contentSizes));
      CHK_ERR( gpuFree(dev_minValues));
      CHK_ERR( gpuFree(dev_massLoss));
      CHK_ERR( gpuFree(dev_mass));
   }
   if (gpu_allocated_batch_maxNeighbours != 0 && second) {
      gpu_allocated_batch_maxNeighbours = 0;
      CHK_ERR( gpuFreeHost(host_vbwcl_neigh));
      CHK_ERR( gpuFree(dev_vbwcl_neigh));
   }
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
      // Always prepare for at least VLASOV_BUFFER_MINCOLUMNS columns
      requiredColumns = max(VLASOV_BUFFER_MINCOLUMNS,gpu_acc_foundColumnsCount);
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
   gpuStream_t stream = gpu_getStream();
   if (columnAllocationCount > 0) {
      // Pointer to host memory struct, contains splitvectors with unified memory data
      cpu_columnOffsetData[cpuThreadID] = new ColumnOffsets(columnAllocationCount);
      CHK_ERR( gpuMallocAsync((void**)&gpu_columnOffsetData[cpuThreadID], sizeof(ColumnOffsets), stream) );
      CHK_ERR( gpuMemcpyAsync(gpu_columnOffsetData[cpuThreadID], cpu_columnOffsetData[cpuThreadID], sizeof(ColumnOffsets), gpuMemcpyHostToDevice, stream));
      CHK_ERR( gpuMallocAsync((void**)&gpu_columns[cpuThreadID], columnAllocationCount*sizeof(Column), stream) );
   }
   // Potential ColumnSet block count container
   const uint c0 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[0];
   const uint c1 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[1];
   const uint c2 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[2];
   std::array<uint, 3> s = {c0,c1,c2};
   std::sort(s.begin(), s.end());
   gpu_acc_columnContainerSize = s[1]*s[2];
   CHK_ERR( gpuMallocAsync((void**)&gpu_columnNBlocks[cpuThreadID], gpu_acc_columnContainerSize*sizeof(vmesh::LocalID), stream) );
}

__host__ void gpu_acc_deallocate_perthread(
   uint cpuThreadID
   ) {
   gpu_acc_allocatedColumns = 0;
   gpu_acc_columnContainerSize = 0;
   gpuStream_t stream = gpu_getStream();
   if (gpu_acc_allocatedColumns > 0) {
      CHK_ERR( gpuFreeAsync(gpu_columns[cpuThreadID],stream) );
      delete cpu_columnOffsetData[cpuThreadID];
      cpu_columnOffsetData[cpuThreadID] = 0;
      CHK_ERR( gpuFreeAsync(gpu_columnOffsetData[cpuThreadID],stream) );
      gpu_columnOffsetData[cpuThreadID] = 0;
   }
   if (gpu_columnNBlocks[cpuThreadID]) {
      CHK_ERR( gpuFreeAsync(gpu_columnNBlocks[cpuThreadID],stream) );
      gpu_columnNBlocks[cpuThreadID] = 0;
   }
}

/*
   Top-level GPU memory allocation function for translation-specific vectors
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void gpu_trans_allocate(
   cuint nAllCells,
   cuint sumOfLengths,
   cuint largestVmesh,
   cuint unionSetSize,
   cuint transGpuBlocks,
   cuint nPencils
   ) {
   gpuStream_t stream = gpu_getStream();
   // Vectors with one entry per cell (prefetch to host)
   if (nAllCells > 0) {
      // Note: this buffer used also in vlasovsolver/gpu_dt.cpp
      if (gpu_allocated_nAllCells == 0) {
         // New allocation
         void *buf0 = malloc(sizeof(split::SplitVector<vmesh::VelocityMesh*>));
         allVmeshPointer = ::new (buf0) split::SplitVector<vmesh::VelocityMesh*>(nAllCells);
         dev_allVmeshPointer = allVmeshPointer->upload<false>(stream);
      } else {
         // Resize
         allVmeshPointer->clear();
         allVmeshPointer->optimizeCPU(stream);
         allVmeshPointer->resize(nAllCells,true);
         dev_allVmeshPointer = allVmeshPointer->upload<false>(stream);
      }
      // Leave on CPU
      gpu_allocated_nAllCells = nAllCells;
   }
   // Vectors with one entry per pencil cell (prefetch to host)
   if (sumOfLengths > 0) {
      if (gpu_allocated_sumOfLengths == 0) {
         // New allocations
         void *buf0 = malloc(sizeof(split::SplitVector<vmesh::VelocityMesh*>));
         void *buf1 = malloc(sizeof(split::SplitVector<vmesh::VelocityBlockContainer*>));
         allPencilsMeshes = ::new (buf0) split::SplitVector<vmesh::VelocityMesh*>(sumOfLengths);
         allPencilsContainers = ::new (buf1) split::SplitVector<vmesh::VelocityBlockContainer*>(sumOfLengths);
         dev_allPencilsMeshes = allPencilsMeshes->upload<false>(stream);
         dev_allPencilsContainers = allPencilsContainers->upload<false>(stream);
      } else {
         // Resize
         allPencilsMeshes->optimizeCPU(stream);
         allPencilsContainers->optimizeCPU(stream);
         allPencilsMeshes->resize(sumOfLengths,true);
         allPencilsContainers->resize(sumOfLengths,true);
         dev_allPencilsMeshes = allPencilsMeshes->upload<false>(stream);
         dev_allPencilsContainers = allPencilsContainers->upload<false>(stream);
      }
      // Leave on CPU
      gpu_allocated_sumOfLengths = sumOfLengths;
   }
   // Set for collecting union of blocks (prefetched to device)
   if (largestVmesh > 0) {
      const vmesh::LocalID HashmapReqSize = ceil(log2((int)largestVmesh)) +2;
      if (gpu_allocated_largestVmeshSizePower == 0) {
         // New allocation
         void *buf0 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
         unionOfBlocksSet = ::new (buf0) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(HashmapReqSize);
         dev_unionOfBlocksSet = unionOfBlocksSet->upload<true>(stream); // <true> == optimize to GPU
         gpu_allocated_largestVmeshSizePower = HashmapReqSize;
      } else {
         // Ensure allocation
         if (HashmapReqSize > gpu_allocated_largestVmeshSizePower) {
            ::delete unionOfBlocksSet;
            void *buf0 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
            unionOfBlocksSet = ::new (buf0) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(HashmapReqSize);
            dev_unionOfBlocksSet = unionOfBlocksSet->upload<true>(stream); // <true> == optimize to GPU
            gpu_allocated_largestVmeshSizePower = HashmapReqSize;
         } else {
            // Ensure map is empty
            unionOfBlocksSet->clear<false>(Hashinator::targets::device,stream, std::pow(2,gpu_allocated_largestVmeshSizePower));
         }
      }
   }
   // Vector into which the set contents are read (prefetched to device)
   if (unionSetSize > 0) {
      if (gpu_allocated_unionSetSize == 0) {
         // New allocation
         void *buf0 = malloc(sizeof(split::SplitVector<vmesh::GlobalID>));
         unionOfBlocks = ::new (buf0) split::SplitVector<vmesh::GlobalID>(unionSetSize);
         unionOfBlocks->clear();
         //unionOfBlocks->optimizeGPU(stream);
         dev_unionOfBlocks = unionOfBlocks->upload<true>(stream); // <true> == optimize to GPU
      } else {
         // Clear is enough
         unionOfBlocks->clear();
         unionOfBlocks->reserve(unionSetSize);
         //unionOfBlocks->optimizeGPU(stream);
         dev_unionOfBlocks = unionOfBlocks->upload<true>(stream); // <true> == optimize to GPU
      }
      gpu_allocated_unionSetSize = unionSetSize;
   }
   // Two temporary buffers, used in-kernel for both reading and writing
   if (transGpuBlocks != 0) {
      if (nPencils == 0) {
         printf("Calling gpu_trans_allocate with transGpuBlocks but without nPencils is not supported!\n");
         abort();
      }
      const uint maxNThreads = gpu_getMaxThreads();
      if (gpu_allocated_trans_pencilBlockData < sumOfLengths*transGpuBlocks*maxNThreads) {
         // Need larger allocation
         if (gpu_allocated_trans_pencilBlockData != 0) {
            CHK_ERR( gpuFree(dev_pencilBlockData) ); // Free old
         }
         // New allocations
         gpu_allocated_trans_pencilBlockData = sumOfLengths*transGpuBlocks*maxNThreads * BLOCK_ALLOCATION_FACTOR;
         CHK_ERR( gpuMalloc((void**)&dev_pencilBlockData, gpu_allocated_trans_pencilBlockData*sizeof(Realf*)) );
      }
      if (gpu_allocated_trans_pencilBlocksCount < nPencils*transGpuBlocks*maxNThreads) {
         // Need larger allocation
         if (gpu_allocated_trans_pencilBlocksCount  != 0) {
            CHK_ERR( gpuFree(dev_pencilBlocksCount) ); // Free old
         }
         // New allocations
         gpu_allocated_trans_pencilBlocksCount = nPencils*transGpuBlocks*maxNThreads * BLOCK_ALLOCATION_FACTOR;
         CHK_ERR( gpuMalloc((void**)&dev_pencilBlocksCount, gpu_allocated_trans_pencilBlocksCount*sizeof(uint)) );
      }
   }
   CHK_ERR( gpuStreamSynchronize(stream) );
}

/* Deallocation at end of simulation */
__host__ void gpu_trans_deallocate() {
   // Deallocate any translation vectors or sets which exist
   if (gpu_allocated_nAllCells != 0) {
      ::delete allVmeshPointer;
   }
   if (gpu_allocated_sumOfLengths != 0) {
      ::delete allPencilsMeshes;
      ::delete allPencilsContainers;
   }
   if (gpu_allocated_largestVmeshSizePower != 0) {
      ::delete unionOfBlocksSet;
   }
   if (gpu_allocated_unionSetSize != 0) {
      ::delete unionOfBlocks;
   }
   if (gpu_allocated_trans_pencilBlockData != 0) {
      CHK_ERR( gpuFree(dev_pencilBlockData) );
   }
   if (gpu_allocated_trans_pencilBlocksCount != 0) {
      CHK_ERR( gpuFree(dev_pencilBlocksCount) );
   }
   gpu_allocated_nAllCells = 0;
   gpu_allocated_sumOfLengths = 0;
   gpu_allocated_largestVmeshSizePower = 0;
   gpu_allocated_unionSetSize = 0;
   gpu_allocated_trans_pencilBlockData = 0;
   gpu_allocated_trans_pencilBlocksCount = 0;
   // Delete also the vectors for pencils for each dimension
   for (uint dimension=0; dimension<3; dimension++) {
      if (DimensionPencils[dimension].gpu_allocated_N) {
         CHK_ERR( gpuFree(DimensionPencils[dimension].gpu_lengthOfPencils) );
         CHK_ERR( gpuFree(DimensionPencils[dimension].gpu_idsStart) );
      }
      if (DimensionPencils[dimension].gpu_allocated_sumOfLengths) {
         CHK_ERR( gpuFree(DimensionPencils[dimension].gpu_sourceDZ) );
         CHK_ERR( gpuFree(DimensionPencils[dimension].gpu_targetRatios) );
      }
   }
}
