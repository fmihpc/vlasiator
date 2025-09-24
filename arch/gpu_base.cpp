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
// Return parameters from kernels, only used for single-cell operations
Real *returnReal[MAXCPUTHREADS];
Realf *returnRealf[MAXCPUTHREADS];
vmesh::LocalID *returnLID[MAXCPUTHREADS];
Real *host_returnReal[MAXCPUTHREADS];
Realf *host_returnRealf[MAXCPUTHREADS];
vmesh::LocalID *host_returnLID[MAXCPUTHREADS];

// Transpose indices for solvers
uint *gpu_cell_indices_to_id;
uint *gpu_block_indices_to_id;
uint *gpu_block_indices_to_probe;

// Pointers to buffers used in acceleration
ColumnOffsets *host_columnOffsetData = NULL, *dev_columnOffsetData = NULL;
Realf **host_blockDataOrdered = NULL, **dev_blockDataOrdered = NULL;
// Counts used in acceleration
size_t gpu_probeFullSize = 0, gpu_probeFlattenedSize = 0;

// Hash map and splitvectors buffers used in batch operations (block adjustment, acceleration)
vmesh::VelocityMesh** host_vmeshes, **dev_vmeshes;
vmesh::VelocityBlockContainer** host_VBCs, **dev_VBCs;
Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** host_allMaps, **dev_allMaps;
split::SplitVector<vmesh::GlobalID> ** host_vbwcl_vec, **dev_vbwcl_vec;
split::SplitVector<vmesh::GlobalID> ** host_lists_with_replace_new, **dev_lists_with_replace_new;
split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_delete, **dev_lists_delete;
split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_to_replace, **dev_lists_to_replace;
split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_with_replace_old, **dev_lists_with_replace_old;
split::SplitVector<vmesh::GlobalID> ** host_vbwcl_neigh, **dev_vbwcl_neigh;

// Counters used in batch operations (block adjustment, acceleration)
vmesh::LocalID* host_nWithContent, *dev_nWithContent;
vmesh::LocalID* host_nBefore, *dev_nBefore;
vmesh::LocalID* host_nAfter, *dev_nAfter;
vmesh::LocalID* host_nBlocksToChange, *dev_nBlocksToChange;
vmesh::LocalID* host_resizeSuccess, *dev_resizeSuccess;
vmesh::LocalID* host_overflownElements, *dev_overflownElements;
vmesh::LocalID* host_nColumns, *dev_nColumns;
vmesh::LocalID* host_nColumnSets, *dev_nColumnSets;
Real* host_minValues, *dev_minValues;
Real* host_massLoss, *dev_massLoss;
Real* host_mass, *dev_mass;
Realf* host_intersections, *dev_intersections; // acceleration only

// Buffers, Vector and set for use in translation
vmesh::VelocityMesh **host_allPencilsMeshes=NULL, **dev_allPencilsMeshes=NULL;
vmesh::VelocityBlockContainer **host_allPencilsContainers=NULL, **dev_allPencilsContainers=NULL;
split::SplitVector<vmesh::GlobalID> *unionOfBlocks=NULL, *dev_unionOfBlocks=NULL;
Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet=NULL, *dev_unionOfBlocksSet=NULL;

// pointers for translation
Realf** dev_pencilBlockData; // Array of pointers into actual block data
uint* dev_pencilBlocksCount; // Array of counters if pencil needs to be propagated for this block or not
GPUMemoryManager gpuMemoryManager;

// Counter for how many parallel vlasov buffers are allocated
uint allocationCount = 0;
// Counter for how large each allocation is
uint gpu_vlasov_allocatedSize[MAXCPUTHREADS] = {0};

// counters for allocated sizes in translation
uint gpu_allocated_sumOfLengths = 0;
uint gpu_allocated_largestVmeshSizePower = 0;
uint gpu_allocated_unionSetSize = 0;
uint gpu_allocated_trans_pencilBlockData = 0;
uint gpu_allocated_trans_pencilBlocksCount = 0;
uint gpu_largest_columnCount = 0;

// batch buffer allocation counters
uint gpu_allocated_batch_nCells = 0;
uint gpu_allocated_batch_maxNeighbours = 0;

// Pointers used in pitch angle diffusion
// Host pointers
Real *host_bValues = nullptr, *host_nu0Values = nullptr, *host_bulkVX = nullptr, *host_bulkVY = nullptr, *host_bulkVZ = nullptr, *host_Ddt = nullptr;
Realf *host_sparsity = nullptr, *dev_densityPreAdjust = nullptr, *dev_densityPostAdjust = nullptr;
size_t *host_cellIdxStartCutoff = nullptr, *host_smallCellIdxArray = nullptr, *host_remappedCellIdxArray = nullptr; // remappedCellIdxArray tells the position of the cell index in the sequence instead of the actual index
// Device pointers
Real *dev_bValues = nullptr, *dev_nu0Values = nullptr, *dev_bulkVX = nullptr, *dev_bulkVY = nullptr, *dev_bulkVZ = nullptr,
   *dev_Ddt = nullptr, *dev_potentialDdtValues = nullptr;
Realf *dev_fmu = nullptr, *dev_dfdt_mu = nullptr, *dev_sparsity = nullptr;
int *dev_fcount = nullptr, *dev_cellIdxKeys = nullptr;
size_t *dev_smallCellIdxArray = nullptr, *dev_remappedCellIdxArray = nullptr, *dev_cellIdxStartCutoff = nullptr, *dev_cellIdxArray = nullptr, *dev_velocityIdxArray = nullptr;
// Counters
size_t latestNumberOfLocalCellsPitchAngle = 0;
int latestNumberOfVelocityCellsPitchAngle = 0;
bool memoryHasBeenAllocatedPitchAngle = false;

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
   }

   CHK_ERR( gpuMalloc((void**)&gpu_cell_indices_to_id, 3*sizeof(uint)) );
   CHK_ERR( gpuMalloc((void**)&gpu_block_indices_to_id, 3*sizeof(uint)) );
   CHK_ERR( gpuMalloc((void**)&gpu_block_indices_to_probe, 3*sizeof(uint)) );
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
   gpu_pitch_angle_diffusion_deallocate();
   // Destroy streams
   const uint maxNThreads = gpu_getMaxThreads();
   for (uint i=0; i<maxNThreads; ++i) {
      CHK_ERR( gpuStreamDestroy(gpuStreamList[i]) );
      CHK_ERR( gpuStreamDestroy(gpuPriorityStreamList[i]) );
      CHK_ERR( gpuFree(returnReal[i]) );
      CHK_ERR( gpuFree(returnRealf[i]) );
      CHK_ERR( gpuFree(returnLID[i]) );
      CHK_ERR( gpuFreeHost(host_returnReal[i]) );
      CHK_ERR( gpuFreeHost(host_returnRealf[i]) );
      CHK_ERR( gpuFreeHost(host_returnLID[i]) );
   }
   CHK_ERR( gpuFree(gpu_cell_indices_to_id) );
   CHK_ERR( gpuFree(gpu_block_indices_to_id) );
   CHK_ERR( gpuFree(gpu_block_indices_to_probe) );
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

   size_t miniBuffers = maxNThreads * (
      8*sizeof(Real) // returnReal
      + 8*sizeof(Realf) // returnRealf
      + 8*sizeof(vmesh::LocalID) // returnLID
      )
      + 9*sizeof(uint) // gpu_cell_indices_to_id, gpu_cell_indices_to_probe, gpu_block_indices_to_id
      + sizeof(std::array<vmesh::MeshParameters,MAX_VMESH_PARAMETERS_COUNT>) // velocityMeshes_upload
      + sizeof(vmesh::MeshWrapper) // MWdev
      + gpu_allocated_moments*sizeof(vmesh::VelocityBlockContainer*) // gpu_moments dev_VBC
      + gpu_allocated_moments*4*sizeof(Real)  // gpu_moments dev_moments1
      + gpu_allocated_moments*3*sizeof(Real); // gpu_moments dev_moments2
   // DT reduction buffers are deallocated every step (GPUTODO, make persistent)

   size_t vlasovBuffers = 0;
   for (uint i=0; i<allocationCount; ++i) {
      vlasovBuffers += gpu_vlasov_allocatedSize[i]
         * WID3 * sizeof(Realf); // gpu_blockDataOrdered[cpuThreadID] // sizes of actual buffers
      vlasovBuffers += sizeof(Realf*); // dev_blockDataOrdered // buffer of pointers to above
   }

   size_t batchBuffers = gpu_allocated_batch_nCells * (
      sizeof(vmesh::VelocityMesh*) // dev_vmeshes
      + sizeof(vmesh::VelocityBlockContainer*) // dev_VBCs
      + 2 * sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*) // dev_allMaps
      + 2 * sizeof(split::SplitVector<vmesh::GlobalID>*) // dev_vbwcl_vec, dev_lists_with_replace_new
      + 3 * sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*) // dev_lists_delete, dev_lists_to_replace, dev_lists_with_replace_old
      + 8 * sizeof(vmesh::LocalID) // dev_nWithContent, dev_nBefore, dev_nAfter, dev_nBlocksToChange, dev_resizeSuccess, dev_overflownElements, dev_nColumns, dev_nColumnSets
      + 3 * sizeof(Real) // dev_minValue, dev_massLoss, dev_mass
      );
   batchBuffers += gpu_allocated_batch_maxNeighbours * sizeof(split::SplitVector<vmesh::GlobalID>*); // dev_vbwcl_neigh

   size_t accBuffers = 0;
   for (uint i=0; i<allocationCount; ++i) {
      accBuffers += sizeof(ColumnOffsets); // dev_columnOffsetData[cpuThreadID]
      if (host_columnOffsetData) {
         accBuffers += host_columnOffsetData[i].capacityInBytes(); // struct contents
      }
   }

   size_t transBuffers = 0;
   if (dev_allPencilsMeshes) {
      transBuffers += gpu_allocated_sumOfLengths * sizeof(vmesh::VelocityMesh*);
   }
   if (dev_allPencilsContainers) {
      transBuffers += gpu_allocated_sumOfLengths * sizeof(vmesh::VelocityBlockContainer*);
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
   const uint maxBlockCount // Largest found vmesh size
   ) {
   // Always prepare for at least VLASOV_BUFFER_MINBLOCKS blocks
   const uint maxBlocksPerCell = max(VLASOV_BUFFER_MINBLOCKS, maxBlockCount);
   if (host_blockDataOrdered == NULL) {
      CHK_ERR( gpuMallocHost((void**)&host_blockDataOrdered,allocationCount*sizeof(Realf*)) );
   }
   if (dev_blockDataOrdered == NULL) {
      CHK_ERR( gpuMalloc((void**)&dev_blockDataOrdered,allocationCount*sizeof(Realf*)) );
   }
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
   // per-buffer allocations
   for (uint i=0; i<allocationCount; ++i) {
      gpu_vlasov_allocate_perthread(i, maxBlocksPerCell);
   }
   // Above function stores buffer pointers in host_blockDataOrdered, copy pointers to dev_blockDataOrdered
   CHK_ERR( gpuMemcpy(dev_blockDataOrdered, host_blockDataOrdered, allocationCount*sizeof(Realf*), gpuMemcpyHostToDevice) );
}

/* Deallocation at end of simulation */
__host__ void gpu_vlasov_deallocate() {
   for (uint i=0; i<allocationCount; ++i) {
      gpu_vlasov_deallocate_perthread(i);
   }
   if (host_blockDataOrdered != NULL) {
      CHK_ERR( gpuFreeHost(host_blockDataOrdered));
   }
   if (dev_blockDataOrdered != NULL) {
      CHK_ERR( gpuFree(dev_blockDataOrdered));
   }
   host_blockDataOrdered = dev_blockDataOrdered = NULL;
}

__host__ uint gpu_vlasov_getSmallestAllocation() {
   uint smallestAllocation = std::numeric_limits<uint>::max();
   for (uint i=0; i<allocationCount; ++i) {
      smallestAllocation = std::min(smallestAllocation,gpu_vlasov_allocatedSize[i]);
   }
   return smallestAllocation;
}

__host__ void gpu_vlasov_allocate_perthread(
   uint allocID,
   uint blockAllocationCount
   ) {
   // Check if we already have allocated enough memory?
   const size_t probeRequirement = sizeof(vmesh::LocalID)*gpu_probeFullSize + gpu_probeFlattenedSize*GPU_PROBEFLAT_N*sizeof(vmesh::LocalID);
   if ( (gpu_vlasov_allocatedSize[allocID] > blockAllocationCount * BLOCK_ALLOCATION_FACTOR) &&
        (gpu_vlasov_allocatedSize[allocID]*WID3*sizeof(Realf) >= probeRequirement) ) {
      return;
   }
   // Potential new allocation with extra padding (including translation multiplier - GPUTODO get rid of this)
   uint newSize = blockAllocationCount * BLOCK_ALLOCATION_PADDING * TRANSLATION_BUFFER_ALLOCATION_FACTOR;
   // Deallocate before new allocation
   gpu_vlasov_deallocate_perthread(allocID);
   gpuStream_t stream = gpu_getStream();

   // Dual use of blockDataOrdered: use also for acceleration probe cube and its flattened version.
   // Calculate required size
   size_t blockDataAllocation = newSize * WID3 * sizeof(Realf);
   // minimum allocation size:
   blockDataAllocation = std::max(blockDataAllocation,probeRequirement);
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
   CHK_ERR( gpuMallocAsync((void**)&host_blockDataOrdered[allocID], blockDataAllocation, stream) );
   // Store size of new allocation (in units blocks)
   gpu_vlasov_allocatedSize[allocID] = blockDataAllocation / (WID3 * sizeof(Realf));
}

__host__ void gpu_vlasov_deallocate_perthread (
   uint allocID
   ) {
   if (gpu_vlasov_allocatedSize[allocID] == 0) {
      return;
   }
   gpuStream_t stream = gpu_getStream();
   CHK_ERR( gpuFreeAsync(host_blockDataOrdered[allocID],stream) );
   gpu_vlasov_allocatedSize[allocID] = 0;
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
      CHK_ERR( gpuMallocHost((void**)&host_nWithContent,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMallocHost((void**)&host_nBefore,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMallocHost((void**)&host_nAfter,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMallocHost((void**)&host_nBlocksToChange,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMallocHost((void**)&host_resizeSuccess,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMallocHost((void**)&host_overflownElements,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMallocHost((void**)&host_nColumns,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMallocHost((void**)&host_nColumnSets,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMallocHost((void**)&host_minValues, gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMallocHost((void**)&host_massLoss, gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMallocHost((void**)&host_mass, gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMallocHost((void**)&host_intersections, gpu_allocated_batch_nCells*4*sizeof(Realf)) );

      CHK_ERR( gpuMalloc((void**)&dev_vmeshes,gpu_allocated_batch_nCells*sizeof(vmesh::VelocityMesh*)) );
      CHK_ERR( gpuMalloc((void**)&dev_VBCs,gpu_allocated_batch_nCells*sizeof(vmesh::VelocityBlockContainer*)) );
      CHK_ERR( gpuMalloc((void**)&dev_allMaps, 2*gpu_allocated_batch_nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_vbwcl_vec, gpu_allocated_batch_nCells*sizeof(split::SplitVector<vmesh::GlobalID>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_lists_with_replace_new, gpu_allocated_batch_nCells*sizeof(split::SplitVector<vmesh::GlobalID>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_lists_delete, gpu_allocated_batch_nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_lists_to_replace, gpu_allocated_batch_nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_lists_with_replace_old, gpu_allocated_batch_nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*)) );
      CHK_ERR( gpuMalloc((void**)&dev_nWithContent,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMalloc((void**)&dev_nBefore,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMalloc((void**)&dev_nAfter,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMalloc((void**)&dev_nBlocksToChange,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMalloc((void**)&dev_resizeSuccess,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMalloc((void**)&dev_overflownElements,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMalloc((void**)&dev_nColumns,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMalloc((void**)&dev_nColumnSets,gpu_allocated_batch_nCells*sizeof(vmesh::LocalID)) );
      CHK_ERR( gpuMalloc((void**)&dev_minValues,gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMalloc((void**)&dev_massLoss, gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMalloc((void**)&dev_mass, gpu_allocated_batch_nCells*sizeof(Real)) );
      CHK_ERR( gpuMalloc((void**)&dev_intersections, gpu_allocated_batch_nCells*4*sizeof(Realf)) );
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
      CHK_ERR( gpuFreeHost(host_nWithContent));
      CHK_ERR( gpuFreeHost(host_nBefore));
      CHK_ERR( gpuFreeHost(host_nAfter));
      CHK_ERR( gpuFreeHost(host_nBlocksToChange));
      CHK_ERR( gpuFreeHost(host_resizeSuccess));
      CHK_ERR( gpuFreeHost(host_overflownElements));
      CHK_ERR( gpuFreeHost(host_nColumns));
      CHK_ERR( gpuFreeHost(host_nColumnSets));
      CHK_ERR( gpuFreeHost(host_minValues));
      CHK_ERR( gpuFreeHost(host_massLoss));
      CHK_ERR( gpuFreeHost(host_mass));
      CHK_ERR( gpuFreeHost(host_intersections));
      CHK_ERR( gpuFree(dev_vmeshes));
      CHK_ERR( gpuFree(dev_VBCs));
      CHK_ERR( gpuFree(dev_allMaps));
      CHK_ERR( gpuFree(dev_vbwcl_vec));
      CHK_ERR( gpuFree(dev_lists_with_replace_new));
      CHK_ERR( gpuFree(dev_lists_delete));
      CHK_ERR( gpuFree(dev_lists_to_replace));
      CHK_ERR( gpuFree(dev_lists_with_replace_old));
      CHK_ERR( gpuFree(dev_nWithContent));
      CHK_ERR( gpuFree(dev_nBefore));
      CHK_ERR( gpuFree(dev_nAfter));
      CHK_ERR( gpuFree(dev_nBlocksToChange));
      CHK_ERR( gpuFree(dev_resizeSuccess));
      CHK_ERR( gpuFree(dev_overflownElements));
      CHK_ERR( gpuFree(dev_nColumns));
      CHK_ERR( gpuFree(dev_nColumnSets));
      CHK_ERR( gpuFree(dev_minValues));
      CHK_ERR( gpuFree(dev_massLoss));
      CHK_ERR( gpuFree(dev_mass));
      CHK_ERR( gpuFree(dev_intersections));
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
   if (host_columnOffsetData == NULL) {
      // This would be preferable as would use pinned memory but fails on exit
      void *buf;
      CHK_ERR( gpuMallocHost((void**)&buf,allocationCount*sizeof(ColumnOffsets)) );
      host_columnOffsetData = new (buf) ColumnOffsets[allocationCount];
   }
   if (dev_columnOffsetData == NULL) {
      CHK_ERR( gpuMalloc((void**)&dev_columnOffsetData,allocationCount*sizeof(ColumnOffsets)) );
   }
   for (uint i=0; i<allocationCount; ++i) {
      gpu_acc_allocate_perthread(i,maxBlockCount);
   }
   // Above function stores buffer pointers in host_blockDataOrdered, copy pointers to dev_blockDataOrdered
   CHK_ERR( gpuMemcpy(dev_columnOffsetData, host_columnOffsetData, allocationCount*sizeof(ColumnOffsets), gpuMemcpyHostToDevice) );
}

/* Deallocation at end of simulation */
__host__ void gpu_acc_deallocate() {
   if (host_columnOffsetData != NULL) {
      // delete[] host_columnOffsetData;
   }
   if (dev_columnOffsetData != NULL) {
      CHK_ERR( gpuFree(dev_columnOffsetData));
   }
   host_columnOffsetData = dev_columnOffsetData = NULL;
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
      CHK_ERR( gpuMemcpyAsync(dev_columnOffsetData+allocID, host_columnOffsetData+allocID, sizeof(ColumnOffsets), gpuMemcpyHostToDevice, stream));
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
      // Use batch allocation
      gpu_batch_allocate(nAllCells);
   }
   // Vectors with one entry per pencil cell (prefetch to host)
   if (sumOfLengths > 0) {
      if (gpu_allocated_sumOfLengths == 0) {
         // New allocations
         CHK_ERR( gpuMalloc((void**)&dev_allPencilsMeshes,sumOfLengths*sizeof(vmesh::VelocityMesh*)) );
         CHK_ERR( gpuMallocHost((void**)&host_allPencilsMeshes,sumOfLengths*sizeof(vmesh::VelocityMesh*)) );
         CHK_ERR( gpuMalloc((void**)&dev_allPencilsContainers,sumOfLengths*sizeof(vmesh::VelocityBlockContainer*)) );
         CHK_ERR( gpuMallocHost((void**)&host_allPencilsContainers,sumOfLengths*sizeof(vmesh::VelocityBlockContainer*)) );
         gpu_allocated_sumOfLengths = sumOfLengths;
      } else if (sumOfLengths > gpu_allocated_sumOfLengths) {
         // Free old
         CHK_ERR( gpuFree(dev_allPencilsMeshes));
         CHK_ERR( gpuFreeHost(host_allPencilsMeshes));
         CHK_ERR( gpuFree(dev_allPencilsContainers));
         CHK_ERR( gpuFreeHost(host_allPencilsContainers));
         // Allocate neow
         CHK_ERR( gpuMalloc((void**)&dev_allPencilsMeshes,sumOfLengths*sizeof(vmesh::VelocityMesh*)) );
         CHK_ERR( gpuMallocHost((void**)&host_allPencilsMeshes,sumOfLengths*sizeof(vmesh::VelocityMesh*)) );
         CHK_ERR( gpuMalloc((void**)&dev_allPencilsContainers,sumOfLengths*sizeof(vmesh::VelocityBlockContainer*)) );
         CHK_ERR( gpuMallocHost((void**)&host_allPencilsContainers,sumOfLengths*sizeof(vmesh::VelocityBlockContainer*)) );
         gpu_allocated_sumOfLengths = sumOfLengths;
      }
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
      if (gpu_allocated_trans_pencilBlockData < sumOfLengths*transGpuBlocks*allocationCount) {
         // Need larger allocation
         if (gpu_allocated_trans_pencilBlockData != 0) {
            CHK_ERR( gpuFree(dev_pencilBlockData) ); // Free old
         }
         // New allocations
         gpu_allocated_trans_pencilBlockData = sumOfLengths*transGpuBlocks*allocationCount * BLOCK_ALLOCATION_FACTOR;
         CHK_ERR( gpuMalloc((void**)&dev_pencilBlockData, gpu_allocated_trans_pencilBlockData*sizeof(Realf*)) );
      }
      if (gpu_allocated_trans_pencilBlocksCount < nPencils*transGpuBlocks*allocationCount) {
         // Need larger allocation
         if (gpu_allocated_trans_pencilBlocksCount  != 0) {
            CHK_ERR( gpuFree(dev_pencilBlocksCount) ); // Free old
         }
         // New allocations
         gpu_allocated_trans_pencilBlocksCount = nPencils*transGpuBlocks*allocationCount * BLOCK_ALLOCATION_FACTOR;
         CHK_ERR( gpuMalloc((void**)&dev_pencilBlocksCount, gpu_allocated_trans_pencilBlocksCount*sizeof(uint)) );
      }
   }
   CHK_ERR( gpuStreamSynchronize(stream) );
}

/* Deallocation at end of simulation */
__host__ void gpu_trans_deallocate() {
   // Deallocate any translation vectors or sets which exist
   if (gpu_allocated_sumOfLengths != 0) {
      CHK_ERR( gpuFree(dev_allPencilsMeshes));
      CHK_ERR( gpuFreeHost(host_allPencilsMeshes));
      CHK_ERR( gpuFree(dev_allPencilsContainers));
      CHK_ERR( gpuFreeHost(host_allPencilsContainers));
      dev_allPencilsMeshes = NULL;
      host_allPencilsMeshes = NULL;
      dev_allPencilsContainers = NULL;
      host_allPencilsContainers = NULL;
      gpu_allocated_sumOfLengths = 0;
   }
   if (gpu_allocated_largestVmeshSizePower != 0) {
      ::delete unionOfBlocksSet;
      gpu_allocated_largestVmeshSizePower = 0;
   }
   if (gpu_allocated_unionSetSize != 0) {
      ::delete unionOfBlocks;
      gpu_allocated_unionSetSize = 0;
   }
   if (gpu_allocated_trans_pencilBlockData != 0) {
      CHK_ERR( gpuFree(dev_pencilBlockData) );
      gpu_allocated_trans_pencilBlockData = 0;
   }
   if (gpu_allocated_trans_pencilBlocksCount != 0) {
      CHK_ERR( gpuFree(dev_pencilBlocksCount) );
      gpu_allocated_trans_pencilBlocksCount = 0;
   }
   // Delete also the vectors for pencils for each dimension
   for (uint dimension=0; dimension<3; dimension++) {
      if (DimensionPencils[dimension].gpu_allocated_N) {
         CHK_ERR( gpuFree(DimensionPencils[dimension].gpu_lengthOfPencils) );
         CHK_ERR( gpuFree(DimensionPencils[dimension].gpu_idsStart) );
         DimensionPencils[dimension].gpu_allocated_N = 0;
      }
      if (DimensionPencils[dimension].gpu_allocated_sumOfLengths) {
         CHK_ERR( gpuFree(DimensionPencils[dimension].gpu_sourceDZ) );
         CHK_ERR( gpuFree(DimensionPencils[dimension].gpu_targetRatios) );
         DimensionPencils[dimension].gpu_allocated_sumOfLengths = 0;
      }
   }
}

void gpu_pitch_angle_diffusion_allocate(size_t numberOfLocalCells, int nbins_v, int nbins_mu, int blocksPerSpatialCell, int totalNumberOfVelocityBlocks) {
   if (numberOfLocalCells <= latestNumberOfLocalCellsPitchAngle && totalNumberOfVelocityBlocks <= latestNumberOfVelocityCellsPitchAngle) {
      return;
   }

   latestNumberOfVelocityCellsPitchAngle = totalNumberOfVelocityBlocks;

   // Allocate device memory
   CHK_ERR( gpuMalloc((void**)&dev_cellIdxArray, totalNumberOfVelocityBlocks*sizeof(size_t)) );
   CHK_ERR( gpuMalloc((void**)&dev_velocityIdxArray, totalNumberOfVelocityBlocks*sizeof(size_t)) );

   if (numberOfLocalCells <= latestNumberOfLocalCellsPitchAngle) {
      return;
   }

   latestNumberOfLocalCellsPitchAngle = numberOfLocalCells;

   if(memoryHasBeenAllocatedPitchAngle){
      gpu_pitch_angle_diffusion_deallocate();
   }

   // Allocate host memory
   CHK_ERR( gpuHostAlloc(&host_bValues, 3*numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuHostAlloc(&host_nu0Values, numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuHostAlloc(&host_sparsity, numberOfLocalCells*sizeof(Realf)) );
   CHK_ERR( gpuHostAlloc(&host_bulkVX, numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuHostAlloc(&host_bulkVY, numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuHostAlloc(&host_bulkVZ, numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuHostAlloc(&host_cellIdxStartCutoff, numberOfLocalCells*sizeof(size_t)) );
   CHK_ERR( gpuHostAlloc(&host_smallCellIdxArray, numberOfLocalCells*sizeof(size_t)) );
   CHK_ERR( gpuHostAlloc(&host_remappedCellIdxArray, numberOfLocalCells*sizeof(size_t)) );
   CHK_ERR( gpuHostAlloc(&host_Ddt, numberOfLocalCells*sizeof(Real)) );

   // Allocate device memory
   CHK_ERR( gpuMalloc((void**)&dev_bValues, 3*numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_nu0Values, numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_sparsity, numberOfLocalCells*sizeof(Realf)) );
   CHK_ERR( gpuMalloc((void**)&dev_dfdt_mu, numberOfLocalCells*nbins_v*nbins_mu*sizeof(Realf)) );
   CHK_ERR( gpuMalloc((void**)&dev_fcount, numberOfLocalCells*nbins_v*nbins_mu*sizeof(int)) );
   CHK_ERR( gpuMalloc((void**)&dev_fmu, numberOfLocalCells*nbins_v*nbins_mu*sizeof(Realf)) );
   CHK_ERR( gpuMalloc((void**)&dev_bulkVX, numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_bulkVY, numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_bulkVZ, numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_densityPreAdjust, numberOfLocalCells*sizeof(Realf)) );
   CHK_ERR( gpuMalloc((void**)&dev_densityPostAdjust, numberOfLocalCells*sizeof(Realf)) );
   CHK_ERR( gpuMalloc((void**)&dev_cellIdxStartCutoff, numberOfLocalCells*sizeof(size_t)) );
   CHK_ERR( gpuMalloc((void**)&dev_smallCellIdxArray, numberOfLocalCells*sizeof(size_t)) );
   CHK_ERR( gpuMalloc((void**)&dev_remappedCellIdxArray, numberOfLocalCells*sizeof(size_t)) );
   CHK_ERR( gpuMalloc((void**)&dev_Ddt, numberOfLocalCells*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_potentialDdtValues, numberOfLocalCells*blocksPerSpatialCell*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_cellIdxKeys, numberOfLocalCells*blocksPerSpatialCell*sizeof(int)) );

   memoryHasBeenAllocatedPitchAngle = true;
}

void gpu_pitch_angle_diffusion_deallocate() {
   // Free memory
   if (dev_bValues) {
      CHK_ERR( gpuFree(dev_bValues) );
   }
   if (dev_nu0Values) {
      CHK_ERR( gpuFree(dev_nu0Values) );
   }
   if (dev_sparsity) {
      CHK_ERR( gpuFree(dev_sparsity) );
   }
   if (dev_dfdt_mu) {
      CHK_ERR( gpuFree(dev_dfdt_mu) );
   }
   if (dev_fcount) {
      CHK_ERR( gpuFree(dev_fcount) );
   }
   if (dev_fmu) {
      CHK_ERR( gpuFree(dev_fmu) );
   }
   if (dev_bulkVX) {
      CHK_ERR( gpuFree(dev_bulkVX) );
   }
   if (dev_bulkVY) {
      CHK_ERR( gpuFree(dev_bulkVY) );
   }
   if (dev_bulkVZ) {
      CHK_ERR( gpuFree(dev_bulkVZ) );
   }
   if (dev_densityPreAdjust) {
      CHK_ERR( gpuFree(dev_densityPreAdjust) );
   }
   if (dev_densityPostAdjust) {
      CHK_ERR( gpuFree(dev_densityPostAdjust) );
   }
   if (dev_cellIdxStartCutoff) {
      CHK_ERR( gpuFree(dev_cellIdxStartCutoff) );
   }
   if (dev_smallCellIdxArray) {
      CHK_ERR( gpuFree(dev_smallCellIdxArray) );
   }
   if (dev_remappedCellIdxArray) {
      CHK_ERR( gpuFree(dev_remappedCellIdxArray) );
   }
   if (dev_Ddt) {
      CHK_ERR( gpuFree(dev_Ddt) );
   }
   if (dev_potentialDdtValues) {
      CHK_ERR( gpuFree(dev_potentialDdtValues) );
   }
   if (dev_cellIdxKeys) {
      CHK_ERR( gpuFree(dev_cellIdxKeys) );
   }
   if (host_bValues) {
      CHK_ERR( gpuFreeHost(host_bValues) );
   }
   if (host_nu0Values) {
      CHK_ERR( gpuFreeHost(host_nu0Values) );
   }
   if (host_sparsity) {
      CHK_ERR( gpuFreeHost(host_sparsity) );
   }
   if (host_bulkVX) {
      CHK_ERR( gpuFreeHost(host_bulkVX) );
   }
   if (host_bulkVY) {
      CHK_ERR( gpuFreeHost(host_bulkVY) );
   }
   if (host_bulkVZ) {
      CHK_ERR( gpuFreeHost(host_bulkVZ) );
   }
   if (dev_velocityIdxArray) {
      CHK_ERR( gpuFree(dev_velocityIdxArray) );
   }
   if (dev_cellIdxArray) {
      CHK_ERR( gpuFree(dev_cellIdxArray) );
   }
   if (host_cellIdxStartCutoff) {
      CHK_ERR( gpuFreeHost(host_cellIdxStartCutoff) );
   }
   if (host_smallCellIdxArray) {
      CHK_ERR( gpuFreeHost(host_smallCellIdxArray) );
   }
   if (host_remappedCellIdxArray) {
      CHK_ERR( gpuFreeHost(host_remappedCellIdxArray) );
   }
   if (host_Ddt) {
      CHK_ERR( gpuFreeHost(host_Ddt) );
   }
}
