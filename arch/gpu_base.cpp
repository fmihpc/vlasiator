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
#include "object_wrapper.h"
#include "../vlasovsolver/gpu_moments.h"
#include "../spatial_cells/velocity_mesh_gpu.h"
#include "../spatial_cells/velocity_block_container.h"
#include "../vlasovsolver/cpu_trans_pencils.hpp"
#include "../vlasovsolver/gpu_moments.h"

#include "logger.h"

// #define MAXCPUTHREADS 64 now in gpu_base.hpp

// Device properties
int gpuMultiProcessorCount = 0;
int blocksPerMP = 0;
int threadsPerMP = 0;

extern Logger logFile;
int myDevice;
int myRank;

// Allocate pointers for per-thread memory regions
gpuStream_t gpuStreamList[MAXCPUTHREADS];
gpuStream_t gpuPriorityStreamList[MAXCPUTHREADS];

// Pointers to buffers used in acceleration
ColumnOffsets *host_columnOffsetData = NULL;
// Counts used in acceleration
size_t gpu_probeFullSize = 0, gpu_probeFlattenedSize = 0, gpu_probeStride = 0;

// Buffers, Vector and set for use in translation
split::SplitVector<vmesh::GlobalID> *unionOfBlocks=NULL, *dev_unionOfBlocks=NULL;
Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet=NULL, *dev_unionOfBlocksSet=NULL;

// Memory manager
GPUMemoryManager gpuMemoryManager;

// Counter for how many parallel vlasov buffers are allocated
uint allocationCount = 0;
// Counter for how large each allocation is
std::vector<uint> gpu_vlasov_allocatedSize;
std::vector<size_t> gpu_vlasov_subPointers;

// counters for allocated sizes in translation
uint gpu_allocated_largestVmeshSizePower = 0;
uint gpu_allocated_unionSetSize = 0;
uint gpu_largest_columnCount = 0;

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

unsigned int nextPowerOfTwo(unsigned int n) {
   if (n == 0) return 1;
   n--; // Handle exact powers of two
   n |= n >> 1;
   n |= n >> 2;
   n |= n >> 4;
   n |= n >> 8;
   n |= n >> 16;
   return n + 1;
}

__host__ void gpu_init_device() {
   const uint maxNThreads = gpu_getMaxThreads();
   int deviceCount;
   // CHK_ERR( gpuFree(0));
   CHK_ERR( gpuGetDeviceCount(&deviceCount) );
   //printf("GPU device count %d with %d threads/streams\n",deviceCount,maxNThreads);

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
   myRank = amps_rank;

   // if only one visible device, assume MPI system handles device visibility and just use the only visible one.
   if (amps_rank == MASTER_RANK) {
      if (deviceCount > 1) {
         // If more than one device is visible, issue warning to user along with suggestion to use SLURM options.
         std::cout << "(Node 0) WARNING! MPI ranks see "<<deviceCount<<" GPU devices each." << std::endl;
         std::cout << "         Recommended usage is to utilize SLURM for showing only single GPU device per MPI rank:" << std::endl;
         std::cout << "         export CUDA_VISIBLE_DEVICES=\\$SLURM_LOCALID" << std::endl;
         std::cout << "         or" << std::endl;
         std::cout << "         export ROCR_VISIBLE_DEVICES=\\$SLURM_LOCALID" << std::endl;
      } else {
         std::cout << "(Node 0) MPI ranks see single GPU device each." << std::endl;
      }
   }
   CHK_ERR( gpuDeviceSynchronize() );
   CHK_ERR( gpuGetDevice(&myDevice) );

   // Decide on number of allocations to prepare
   const uint nBaseCells = P::xcells_ini * P::ycells_ini * P::zcells_ini;
   allocationCount = (nBaseCells == 1) ? 1 : P::GPUallocations;

   // Get device properties
   gpuDeviceProp prop;
   CHK_ERR( gpuGetDeviceProperties(&prop, myDevice) );
   gpuMultiProcessorCount = prop.multiProcessorCount;
   threadsPerMP = prop.maxThreadsPerMultiProcessor;
   #if defined(USE_GPU) && defined(__CUDACC__)
   CHK_ERR( gpuDeviceGetAttribute(&blocksPerMP, gpuDevAttrMaxBlocksPerMultiprocessor, myDevice) );
   #endif
   #if defined(USE_GPU) && defined(__HIP_PLATFORM_HCC___)
   blocksPerMP = threadsPerMP/GPUTHREADS; // This should be the maximum number of wavefronts per CU
   #endif


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
   CREATE_SUBPOINTERS(gpuMemoryManager, returnReal, maxNThreads);
   CREATE_SUBPOINTERS(gpuMemoryManager, returnRealf, maxNThreads);
   CREATE_SUBPOINTERS(gpuMemoryManager, returnLID, maxNThreads);
   CREATE_SUBPOINTERS(gpuMemoryManager, host_returnReal, maxNThreads);
   CREATE_SUBPOINTERS(gpuMemoryManager, host_returnRealf, maxNThreads);
   CREATE_SUBPOINTERS(gpuMemoryManager, host_returnLID, maxNThreads);


   int *leastPriority = new int; // likely 0
   int *greatestPriority = new int; // likely -1
   CHK_ERR( gpuDeviceGetStreamPriorityRange (leastPriority, greatestPriority) );
   if (*leastPriority==*greatestPriority) {
      printf("Warning when initializing GPU streams: minimum and maximum stream priority are identical! %d == %d \n",*leastPriority, *greatestPriority);
   }
   for (uint i=0; i<maxNThreads; ++i) {
      CHK_ERR( gpuStreamCreateWithPriority(&(gpuStreamList[i]), gpuStreamDefault, *leastPriority) );
      CHK_ERR( gpuStreamCreateWithPriority(&(gpuPriorityStreamList[i]), gpuStreamDefault, *greatestPriority) );

      SUBPOINTER_HOST_ALLOCATE(gpuMemoryManager, host_returnReal, i, 8*sizeof(Real));
      SUBPOINTER_HOST_ALLOCATE(gpuMemoryManager, host_returnRealf, i, 8*sizeof(Realf));
      SUBPOINTER_HOST_ALLOCATE(gpuMemoryManager, host_returnLID, i, 8*sizeof(vmesh::LocalID));

      SUBPOINTER_ALLOCATE(gpuMemoryManager, returnReal, i, 8*sizeof(Real));
      SUBPOINTER_ALLOCATE(gpuMemoryManager, returnRealf, i, 8*sizeof(Realf));
      SUBPOINTER_ALLOCATE(gpuMemoryManager, returnLID, i, 8*sizeof(vmesh::LocalID));
   }

   CREATE_UNIQUE_POINTER(gpuMemoryManager, gpu_cell_indices_to_id);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, gpu_block_indices_to_id);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, gpu_block_indices_to_probe);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, my_test_pointer);

   ALLOCATE_GPU(gpuMemoryManager, gpu_cell_indices_to_id, 3*sizeof(uint));
   ALLOCATE_GPU(gpuMemoryManager, gpu_block_indices_to_id, 3*sizeof(uint));
   ALLOCATE_GPU(gpuMemoryManager, gpu_block_indices_to_probe, 3*sizeof(uint));
   CHK_ERR( gpuDeviceSynchronize() );

   // Using just a single context for whole MPI task
}

__host__ void gpu_clear_device() {
   // Deallocate temporary buffers
   gpu_acc_deallocate();
   gpu_vlasov_deallocate();
   gpu_trans_deallocate();
   // Destroy streams
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      CHK_ERR( gpuStreamDestroy(gpuStreamList[i]) );
      CHK_ERR( gpuStreamDestroy(gpuPriorityStreamList[i]) );
   }
   CHK_ERR( gpuDeviceSynchronize() );
   gpuMemoryManager.freeAll();
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

__host__ uint gpu_getAllocationCount() {
   return allocationCount;
}

/*
   Memory reporting function
*/
int gpu_reportMemory(const size_t local_cells_capacity, const size_t ghost_cells_capacity,
                     const size_t local_cells_size, const size_t ghost_cells_size) {
   /* Gather total CPU and GPU buffer sizes. Rank 0 reports details,
      all ranks return sum.
   */
   uint maxNThreads = gpu_getMaxThreads();

   size_t miniBuffers = 
      sizeof(std::array<vmesh::MeshParameters,MAX_VMESH_PARAMETERS_COUNT>) // velocityMeshes_upload
      + sizeof(vmesh::MeshWrapper); // MWdev
   // DT reduction buffers are deallocated every step (GPUTODO, make persistent)

   size_t vlasovBuffers = 0;
   size_t batchBuffers = 0;

   size_t accBuffers = 0;
   for (uint i=0; i<allocationCount; ++i) {
      if (host_columnOffsetData) {
         accBuffers += host_columnOffsetData[i].capacityInBytes(); // struct contents
      }
   }

   size_t transBuffers = 0;
   if (unionOfBlocksSet) {
      transBuffers += sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>);
      transBuffers += unionOfBlocksSet->bucket_count() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
   }
   if (unionOfBlocks) {
      transBuffers += sizeof(split::SplitVector<vmesh::GlobalID>);
      transBuffers += unionOfBlocks->capacity() * sizeof(vmesh::GlobalID);
   }
   // Remote neighbor contribution buffers are in unified memory but deallocated after each use

   size_t memoryManagerCapacity = gpuMemoryManager.totalGpuAllocation();

   size_t free_byte ;
   size_t total_byte ;
   CHK_ERR( gpuMemGetInfo( &free_byte, &total_byte) );
   size_t used_mb = (total_byte-free_byte)/(1024*1024);
   size_t sum_mb = (miniBuffers+batchBuffers+vlasovBuffers+accBuffers+transBuffers+local_cells_capacity+ghost_cells_capacity+memoryManagerCapacity)/(1024*1024);
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
      logFile<<"     Memory manager:        "<<memoryManagerCapacity/(1024*1024)<<" Mbytes"<<std::endl;
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
   const uint maxBlockCount // Largest found vmesh size
   ) {
   // Always prepare for at least VLASOV_BUFFER_MINBLOCKS blocks
   const uint maxBlocksPerCell = max(VLASOV_BUFFER_MINBLOCKS, maxBlockCount);
   
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_blockDataOrdered);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_blockDataOrdered);
   ALLOCATE_GPU(gpuMemoryManager, dev_blockDataOrdered, allocationCount*sizeof(Realf*));
   HOST_ALLOCATE_GPU(gpuMemoryManager, host_blockDataOrdered, allocationCount*sizeof(Realf*));

   // per-buffer allocations
   for (uint i=0; i<allocationCount; ++i) {
      gpu_vlasov_allocate_perthread(i, maxBlocksPerCell);
   }

   // Above function stores buffer pointers in host_blockDataOrdered, copy pointers to dev_blockDataOrdered
   CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, Realf*, dev_blockDataOrdered), GET_POINTER(gpuMemoryManager, Realf*, host_blockDataOrdered), allocationCount*sizeof(Realf*), gpuMemcpyHostToDevice) );
}

/*
   Top-level GPU memory allocation function.
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void gpu_calculateProbeAllocation(
   const uint maxBlockCount // Largest found vmesh size
   ) {
   // Always prepare for at least VLASOV_BUFFER_MINBLOCKS blocks
   const uint maxBlocksPerCell = max(VLASOV_BUFFER_MINBLOCKS, maxBlockCount);
   
   // Evaluate required size for acceleration probe cube (based on largest population)
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      const uint c0 = (*vmesh::getMeshWrapper()->velocityMeshes)[popID].gridLength[0];
      const uint c1 = (*vmesh::getMeshWrapper()->velocityMeshes)[popID].gridLength[1];
      const uint c2 = (*vmesh::getMeshWrapper()->velocityMeshes)[popID].gridLength[2];
      std::array<uint, 3> s = {c0,c1,c2};
      std::sort(s.begin(), s.end());
      // Round values up to nearest 2*Hashinator::defaults::MAX_BLOCKSIZE
      size_t probeCubeExtentsFull = s[0]*s[1]*s[2];
      probeCubeExtentsFull = 2*Hashinator::defaults::MAX_BLOCKSIZE * (1 + ((probeCubeExtentsFull - 1) / (2*Hashinator::defaults::MAX_BLOCKSIZE)));
      gpu_probeFullSize = std::max(gpu_probeFullSize,probeCubeExtentsFull);
      size_t probeCubeExtentsFlat = s[1]*s[2];
      probeCubeExtentsFlat = 2*Hashinator::defaults::MAX_BLOCKSIZE * (1 + ((probeCubeExtentsFlat - 1) / (2*Hashinator::defaults::MAX_BLOCKSIZE)));
      gpu_probeFlattenedSize = std::max(gpu_probeFlattenedSize,probeCubeExtentsFlat);
   }

   size_t probeAllocation = sizeof(vmesh::LocalID)*gpu_probeFullSize + gpu_probeFlattenedSize*GPU_PROBEFLAT_N*sizeof(vmesh::LocalID);
   /*
   CUDA C Programming Guide
   6.3.2. Device Memory Accesses (June 2025)
   "Any address of a variable residing in global memory or returned by one of the memory allocation routines from the driver or
   runtime API is always aligned to at least 256 bytes."

   ROCm documentation
   HIP 6.4.43483 Documentation for hipMallocPitch
   "Currently the alignment is set to 128 bytes"

   Thus, our mallocs should be in increments of 256 bytes. WID3 is at least 64, and len(Realf) is at least 4, so this is true in all
   cases. Still, let us ensure (just to be sure) that probe cube addressing does not break alignment.
   And in fact let's use the block memory size as the stride.
   */
   probeAllocation = (1 + ((probeAllocation - 1) / (WID3 * sizeof(Realf)))) * (WID3 * sizeof(Realf));
   gpu_probeStride = max(gpu_probeStride,probeAllocation);
}

/* Deallocation at end of simulation */
__host__ void gpu_vlasov_deallocate() {
   while(gpu_vlasov_allocatedSize.size() < allocationCount){ //Make sure the gpu_vlasov_allocatedSize has enough elements
      gpu_vlasov_allocatedSize.push_back(0);
   }
   for (uint i=0; i<allocationCount; ++i) {
      gpu_vlasov_allocatedSize[i] = 0;
   }
}

__host__ uint gpu_vlasov_getSmallestAllocation() {
   uint smallestAllocation = std::numeric_limits<uint>::max();
   while(gpu_vlasov_allocatedSize.size() < allocationCount){ //Make sure the gpu_vlasov_allocatedSize has enough elements
      gpu_vlasov_allocatedSize.push_back(0);
   }
   for (uint i=0; i<allocationCount; ++i) {
      smallestAllocation = std::min(smallestAllocation,gpu_vlasov_allocatedSize[i]);
   }
   return smallestAllocation;
}

__host__ void gpu_vlasov_allocate_perthread(
   uint allocID,
   uint blockAllocationCount
   ) {
   while(gpu_vlasov_allocatedSize.size() < allocationCount){ //Make sure the gpu_vlasov_allocatedSize has enough elements
      gpu_vlasov_allocatedSize.push_back(0);
   }
   while(gpu_vlasov_subPointers.size() < allocationCount){ //Make sure the gpu_vlasov_subPointers has enough elements
      gpu_vlasov_subPointers.push_back(0);
   }

   // Dual use of blockDataOrdered: use also for acceleration probe cube and its flattened version.
   // Calculate required size
   size_t blockDataAllocation = blockAllocationCount * WID3 * sizeof(Realf);
   /*
     CUDA C Programming Guide
     6.3.2. Device Memory Accesses (June 2025)
     "Any address of a variable residing in global memory or returned by one of the memory allocation routines from the driver or
     runtime API is always aligned to at least 256 bytes."

     ROCm documentation
     HIP 6.4.43483 Documentation for hipMallocPitch
     "Currently the alignment is set to 128 bytes"

     Thus, our mallocs should be in increments of 256 bytes. WID3 is at least 64, and len(Realf) is at least 4, so this is true in all
     cases. Still, let us ensure (just to be sure) that probe cube addressing does not break alignment.
     And in fact let's use the block memory size as the stride.
   */
   blockDataAllocation = (1 + ((blockDataAllocation - 1) / (WID3 * sizeof(Realf)))) * (WID3 * sizeof(Realf));

   gpuMemoryManager.createPointer(gpu_vlasov_subPointers[allocID]);
   bool reAllocated = gpuMemoryManager.allocate(gpu_vlasov_subPointers[allocID], blockDataAllocation);
   SET_SUBPOINTER(gpuMemoryManager, Realf, host_blockDataOrdered, allocID,  gpu_vlasov_subPointers[allocID]);

   // Store size of new allocation (in units blocks)
   if (reAllocated) {
      gpu_vlasov_allocatedSize[allocID] = blockDataAllocation / (WID3 * sizeof(Realf));
   }
}

/** Allocation and deallocation for pointers used by batch operations in block adjustment */
__host__ void gpu_batch_allocate(uint nCells, uint maxNeighbours) {

   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_vmeshes);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_VBCs);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_allMaps);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_vbwcl_vec);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_lists_with_replace_new);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_lists_delete);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_lists_to_replace);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_lists_with_replace_old);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_nBefore);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_nAfter);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_nBlocksToChange);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_resizeSuccess);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_overflownElements);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_minValues);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_massLoss);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, host_intersections);

   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_vmeshes, nCells*sizeof(vmesh::VelocityMesh*), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_VBCs, nCells*sizeof(vmesh::VelocityBlockContainer*), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_allMaps, 2*nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), BLOCK_ALLOCATION_FACTOR); // note double size
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_vbwcl_vec, nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_lists_with_replace_new, nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_lists_delete, nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_lists_to_replace, nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_lists_with_replace_old, nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_nBefore, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_nAfter, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_nBlocksToChange, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_resizeSuccess, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_overflownElements, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_minValues, nCells*sizeof(Real), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_massLoss, nCells*sizeof(Real), BLOCK_ALLOCATION_FACTOR);
   HOST_ALLOCATE_WITH_BUFFER(gpuMemoryManager, host_intersections, nCells*4*sizeof(Realf), BLOCK_ALLOCATION_FACTOR);

   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_vmeshes);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_VBCs);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_allMaps);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_vbwcl_vec);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_lists_with_replace_new);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_lists_delete);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_lists_to_replace);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_lists_with_replace_old);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_nBefore);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_nAfter);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_nBlocksToChange);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_resizeSuccess);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_overflownElements);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_minValues);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_massLoss);
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_intersections);

   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_vmeshes, nCells*sizeof(vmesh::VelocityMesh*), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_VBCs, nCells*sizeof(vmesh::VelocityBlockContainer*), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_allMaps, 2*nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_vbwcl_vec, nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_lists_with_replace_new, nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_lists_delete, nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_lists_to_replace, nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_lists_with_replace_old, nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_nBefore, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_nAfter, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_nBlocksToChange, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_resizeSuccess, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_overflownElements, nCells*sizeof(vmesh::LocalID), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_minValues, nCells*sizeof(Real), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_massLoss, nCells*sizeof(Real), BLOCK_ALLOCATION_FACTOR);
   ALLOCATE_WITH_BUFFER(gpuMemoryManager, dev_intersections, nCells*4*sizeof(Realf), BLOCK_ALLOCATION_FACTOR);
}

/*
   Top-level GPU memory allocation function for acceleration-specific column data.
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void gpu_acc_allocate(
   uint maxBlockCount
   ) {
   if (host_columnOffsetData == NULL) {
      // This would be preferable as would use pinned memory but fails on exit
      void *buf;
      CHK_ERR( gpuMallocHost((void**)&buf,allocationCount*sizeof(ColumnOffsets)) );
      host_columnOffsetData = new (buf) ColumnOffsets[allocationCount];
   }
   CREATE_UNIQUE_POINTER(gpuMemoryManager, dev_columnOffsetData);
   ALLOCATE_GPU(gpuMemoryManager, dev_columnOffsetData, allocationCount*sizeof(ColumnOffsets));
   for (uint i=0; i<allocationCount; ++i) {
      gpu_acc_allocate_perthread(i,maxBlockCount);
   }
   // Above function stores buffer pointers in host_blockDataOrdered, copy pointers to dev_blockDataOrdered
   CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, ColumnOffsets, dev_columnOffsetData), host_columnOffsetData, allocationCount*sizeof(ColumnOffsets), gpuMemcpyHostToDevice) );
}

/* Deallocation at end of simulation */
__host__ void gpu_acc_deallocate() {
   if (host_columnOffsetData != NULL) {
      // delete[] host_columnOffsetData;
   }
   host_columnOffsetData = NULL;
}

/*
   Mid-level GPU memory allocation function for acceleration-specific column data.
   Supports calling within threaded regions and async operations.
 */
__host__ void gpu_acc_allocate_perthread(
   uint allocID,
   uint firstAllocationCount, // This is treated as maxBlockCount, unless the next
   //value is nonzero, in which case it is the column allocation count
   uint columnSetAllocationCount
   ) {
   uint columnAllocationCount;
   if (columnSetAllocationCount==0) {
      /*
        Estimate column count from maxblockcount, non-critical if ends up being too small.
        This makes a rough guess that we have a cubic velocity space domain, and thus one edge
        of it is the cubic root of the blocks count, and the area is the square of that. Thus,
        we take the two-thirds power of the block count, and multiply by the padding multiplier
        to be a bit on the safer side.
      */
      columnAllocationCount = BLOCK_ALLOCATION_PADDING * std::pow(firstAllocationCount,0.666);
      // Ensure a minimum value.
      columnAllocationCount = std::max(columnAllocationCount,(uint)VLASOV_BUFFER_MINCOLUMNS);
      columnSetAllocationCount = columnAllocationCount;
   } else {
      columnAllocationCount = firstAllocationCount;
      // Update tracker
      gpu_largest_columnCount = std::max(columnAllocationCount,gpu_largest_columnCount);
   }

   // columndata contains several splitvectors. columnData is host/device, but splitvector contents are unified.
   gpuStream_t stream = gpu_getStream();
   // Reallocate if necessary
   if ( (columnAllocationCount > host_columnOffsetData[allocID].capacityCols()) ||
        (columnSetAllocationCount > host_columnOffsetData[allocID].capacityColSets()) ) {
      // Also set size to match input
      host_columnOffsetData[allocID].setSizes(columnAllocationCount*BLOCK_ALLOCATION_PADDING, columnSetAllocationCount*BLOCK_ALLOCATION_PADDING);
      CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, ColumnOffsets, dev_columnOffsetData)+allocID, host_columnOffsetData+allocID, sizeof(ColumnOffsets), gpuMemcpyHostToDevice, stream));
   }
}

/*
   Top-level GPU memory allocation function for translation-specific vectors
   This is called from within non-threaded regions so does not perform async.
 */
__host__ void gpu_trans_allocate(
   cuint nAllCells,
   cuint largestVmesh,
   cuint unionSetSize
   ) {
   gpuStream_t stream = gpu_getStream();
   // Vectors with one entry per cell (prefetch to host)
   if (nAllCells > 0) {
      // Use batch allocation
      gpu_batch_allocate(nAllCells);
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
   CHK_ERR( gpuStreamSynchronize(stream) );
}

/* Deallocation at end of simulation */
__host__ void gpu_trans_deallocate() {
   // Deallocate any translation vectors or sets which exist
   if (gpu_allocated_largestVmeshSizePower != 0) {
      ::delete unionOfBlocksSet;
      gpu_allocated_largestVmeshSizePower = 0;
   }
   if (gpu_allocated_unionSetSize != 0) {
      ::delete unionOfBlocks;
      gpu_allocated_unionSetSize = 0;
   }
}
