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

#ifndef VLASIATOR_BLOCK_ADJUST_KERNELS_HPP
#define VLASIATOR_BLOCK_ADJUST_KERNELS_HPP

#ifdef USE_WARPACCESSORS
 #define USE_BATCH_WARPACCESSORS
#endif

// __launch_bounds__(MAX_THREADS_PER_BLOCK, MIN_BLOCKS_PER_MP)
#ifdef __CUDACC__
#define WARPS_PER_MP 64
#define FULLBLOCKS_PER_MP THREADS_PER_MP/Hashinator::defaults::MAX_BLOCKSIZE
#define WID3S_PER_MP (2048/WID3)
#endif
#ifdef __HIP_PLATFORM_HCC___
#define WARPS_PER_MP 8
#define FULLBLOCKS_PER_MP 1
#define WID3S_PER_MP 7 // because of batch_update_velocity_blocks_kernel
#endif

/** GPU kernel for identifying which blocks have relevant content */
__global__ void __launch_bounds__(WID3,WID3S_PER_MP) batch_update_velocity_block_content_lists_kernel (
   const vmesh::VelocityMesh*  __restrict__  const *vmeshes,
   const vmesh::VelocityBlockContainer*  __restrict__  const *blockContainers,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* *allMaps,
   const Real* __restrict__ velocity_block_min_values,
   const bool gatherMass,
   Real* dev_mass
   ) {
   //const uint nCells = gridDim.y;
   const int cellIndex = blockIdx.y;
   const int blockiStart = blockIdx.x;

   const uint ti = threadIdx.x;

   const vmesh::VelocityMesh* __restrict__ vmesh = vmeshes[cellIndex];
   const vmesh::VelocityBlockContainer* __restrict__ blockContainer = blockContainers[cellIndex];
   const Real velocity_block_min_value = velocity_block_min_values[cellIndex];
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwcl_map = allMaps[2*cellIndex];
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwncl_map = allMaps[2*cellIndex+1];

   #define warpsPerBlockBatchContent WID3/GPUTHREADS

   __shared__ int has_content[warpsPerBlockBatchContent];
   __shared__ Real gathered_mass[WID3];
   const uint nBlocks = vmesh->size();
   #ifdef DEBUG_SPATIAL_CELL
   if (nBlocks != blockContainer->size()) {
      if (ti==0) {
         printf("VBC and vmesh size mismatch in batch_update_velocity_block_content_lists_kernel!\n");
      }
      assert(0);
   }
   #endif
   const uint blockLID = blockiStart;
   {
      if (blockLID >= nBlocks) {
         return;
      }
      // Check each velocity cell if it is above the threshold
      const Realf* __restrict__ avgs = blockContainer->getData(blockLID);

      const vmesh::GlobalID blockGID = vmesh->getGlobalID(blockLID);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == vmesh->invalidGlobalID()) {
         if (ti==0) {
            printf("Invalid GID encountered in batch_update_velocity_block_content_lists_kernel!\n");
         }
         assert(0);
      }
      if (blockLID == vmesh->invalidLocalID()) {
         if (ti==0) {
            printf("Invalid LID encountered in batch_update_velocity_block_content_lists_kernel!\n");
         }
         assert(0);
      }
      #endif

      bool hasContentThread = (avgs[ti] >= velocity_block_min_value);

      if (gatherMass) {
         gathered_mass[ti] = avgs[ti];
         __syncthreads();
         // Perform loop over all elements to gather total mass
         for (unsigned int s=WID3/2; s>0; s>>=1) {
            if (ti < s) {
               gathered_mass[ti] += gathered_mass[ti + s];
            }
            __syncthreads();
         }
      }
      
      const int indexInsideWarp = ti % GPUTHREADS;
      const int warpIndex = ti / GPUTHREADS;

      // Check for content with two consecutive warp votes
      hasContentThread = gpuKernelAny(0xFFFFFFFF, hasContentThread);

      if (WID3 > GPUTHREADS) {
         if (indexInsideWarp == 0) {
            has_content[warpIndex] = hasContentThread;
         }
         __syncthreads();

         if (warpIndex == 0) {
            hasContentThread = (indexInsideWarp < warpsPerBlockBatchContent) ? has_content[indexInsideWarp] : false;
            hasContentThread = gpuKernelAny(0xFFFFFFFF, hasContentThread);
         }
      }

      #ifdef USE_BATCH_WARPACCESSORS
      // Insert into map only from threads 0...WARPSIZE
      if (ti < GPUTHREADS) {
         if (ti == 0) {
            has_content[0] = hasContentThread;
         }
         __syncthreads();
         if (hasContentThread) {
            vbwcl_map->warpInsert(blockGID,blockLID,ti);
         } else {
            vbwncl_map->warpInsert(blockGID,blockLID,ti);
         }
      }
      #else
      // Insert into map only from thread 0
      if (ti == 0) {
         if (hasContentThread) {
            vbwcl_map->set_element(blockGID,blockLID);
         } else {
            vbwncl_map->set_element(blockGID,blockLID);
         }
      }
      #endif
      // Store gathered mass as atomic from one thread per block
      if (gatherMass && (ti == 0)) {
         Real old = atomicAdd(&dev_mass[cellIndex], gathered_mass[0]);
      }
   }
}

/*
 * Resets all elements in all provided hashmaps to EMPTY, VAL_TYPE()
 */
__global__ void __launch_bounds__(Hashinator::defaults::MAX_BLOCKSIZE, FULLBLOCKS_PER_MP) batch_reset_all_to_empty(
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>**maps
   ) {
   //launch parameters: dim3 grid(blocksNeeded,nCells,2);
   const size_t hashmapIndex = blockIdx.y * 2 + blockIdx.z;
   const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
   const size_t stride = gridDim.x * blockDim.x;
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* thisMap = maps[hashmapIndex];
   const size_t len = thisMap->bucket_count();
   const vmesh::GlobalID emptybucket = thisMap->get_emptybucket();
   Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>* dst = thisMap->expose_bucketdata<false>();

   for (size_t bucketIndex = tid; bucketIndex < len; bucketIndex += stride) {
      dst[bucketIndex].first = emptybucket;
   }

   //Thread 0 resets fill
   if (tid==0) {
      Hashinator::Info *info = thisMap->expose_mapinfo<false>();
      info->fill=0;
   }
}

/*
 * Reads sizes of hashmaps, compares with capacities of provided vectors, and sets the provided buffer
 * to indicate if the vector needs recapacitating. Assumes the required_capacities buffer has been
 * memset to zero before this kernel is called.
 */
__global__ void __launch_bounds__(Hashinator::defaults::MAX_BLOCKSIZE, FULLBLOCKS_PER_MP) check_vector_capacities(
   const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* __restrict__ const *maps,
   const split::SplitVector<vmesh::GlobalID>* __restrict__ const *vecs,
   vmesh::LocalID *required_capacities
   ) {
   const size_t index = threadIdx.x + blockIdx.x * blockDim.x;
   const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*  __restrict__ thisMap = maps[2*index];
   const split::SplitVector<vmesh::GlobalID>* __restrict__ thisVec = vecs[index];
   const size_t mapSize = thisMap->size();
   if (mapSize > thisVec->capacity()) {
      required_capacities[index] = mapSize;
   }
}

/*
 * Extracts keys (GIDs, if firstonly is true) or key-value pairs (GID-LID pairs)
 * from all provided hashmaps to provided splitvectors, and stores the vector size in an array.
 */
template <typename Rule, typename ELEMENT, bool FIRSTONLY=false>
__global__ void __launch_bounds__(Hashinator::defaults::MAX_BLOCKSIZE, FULLBLOCKS_PER_MP) extract_GIDs_kernel(
   const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* __restrict__ const *input_maps, // buffer of pointers to source maps
   split::SplitVector<ELEMENT> **output_vecs,
   vmesh::LocalID* output_sizes,
   Rule rule,
   const vmesh::VelocityMesh* __restrict__ const *rule_meshes, // buffer of pointers to vmeshes, sizes used by rules
   const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* __restrict__ const *rule_maps,
   const split::SplitVector<vmesh::GlobalID>* __restrict__ const *rule_vectors
   ) {
   //launch parameters: dim3 grid(nMaps,1,1); // As this is a looping reduction
   const size_t cellIndex = blockIdx.x;
   const size_t hashmapIndex = 2*blockIdx.x; // Assumes maps are with a stride of two due to allMaps buffer holding two for each cell
   if (input_maps[hashmapIndex]==0) {
      return; // Early return for invalid cells
   }
   const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* __restrict__ thisMap = input_maps[hashmapIndex];
   split::SplitVector<ELEMENT> *outputVec = output_vecs[cellIndex];

   // Threshold value used by some rules
   const vmesh::LocalID threshold = rule_meshes[cellIndex]->size()
      + rule_vectors[cellIndex]->size() - rule_maps[hashmapIndex]->size();

   const vmesh::LocalID  invalidLID = rule_meshes[cellIndex]->invalidLocalID();
   const vmesh::GlobalID invalidGID = rule_meshes[cellIndex]->invalidGlobalID();

   // This must be equal to at least both WARPLENGTH and MAX_BLOCKSIZE/WARPLENGTH
   __shared__ uint32_t warpSums[WARPLENGTH];
   __shared__ uint32_t outputCount;
   const int tid = threadIdx.x;
   const int wid = tid / WARPLENGTH;
   const int w_tid = tid % WARPLENGTH;
   //const int warpsPerBlock = BLOCKSIZE / WARPLENGTH;
   const size_t warpsPerBlock = blockDim.x / WARPLENGTH;
   // zero init shared buffer
   if (wid == 0) {
      warpSums[w_tid] = 0;
   }
   __syncthreads();
   // full warp votes for rule-> mask = [01010101010101010101010101010101]
   int64_t remaining = thisMap->bucket_count();
   const uint capacity = outputVec->capacity();
   uint32_t outputSize = 0;
   uint32_t inputOffset = 0;
   // Initial pointers into data
   //Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID> *input = thisMap->expose_bucketdata<false>();
   ELEMENT* output = outputVec->data();
   // Start loop
   while (remaining > 0) {
      const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>* __restrict__ input = thisMap->expose_bucketdata<false>();
      const int current = remaining > blockDim.x ? blockDim.x : remaining;
      __syncthreads();
      const int active = (tid < current) ? rule(thisMap, input[inputOffset + tid], threshold, invalidLID, invalidGID) : false;
      const auto mask = split::s_warpVote(active == 1, SPLIT_VOTING_MASK);
      const auto warpCount = split::s_pop_count(mask);
      if (w_tid == 0) {
         warpSums[wid] = warpCount;
      }
      __syncthreads();
      // Figure out the total here because we overwrite shared mem later
      if (wid == 0) {
         // ceil int division
         int activeWARPS = nextPow2(1 + ((current - 1) / WARPLENGTH));
         auto reduceCounts = [activeWARPS](int localCount) -> int {
                                for (int i = activeWARPS / 2; i > 0; i = i / 2) {
                                   localCount += split::s_shuffle_down(localCount, i, SPLIT_VOTING_MASK);
                                }
                                return localCount;
                             };
         auto localCount = warpSums[w_tid];
         const int totalCount = reduceCounts(localCount);
         if (w_tid == 0) {
            outputCount = totalCount;
            outputSize += totalCount;
            assert((outputSize <= capacity) && "extract_GIDs_kernel ran out of capacity!");
            outputVec->device_resize(outputSize);
         }
      }
      // Prefix scan WarpSums on the first warp
      if (wid == 0) {
         auto value = warpSums[w_tid];
         for (uint d = 1; d < warpsPerBlock; d = 2 * d) {
            int res = split::s_shuffle_up(value, (int)d, SPLIT_VOTING_MASK);
            if (tid % warpsPerBlock >= d) {
               value += res;
            }
         }
         warpSums[w_tid] = value;
      }
      __syncthreads();
      auto offset = (wid == 0) ? 0 : warpSums[wid - 1];
      auto pp = split::s_pop_count(mask & ((ONE << w_tid) - ONE));
      const auto warpTidWriteIndex = offset + pp;
      if (active) {
         if constexpr (FIRSTONLY) {
            output[warpTidWriteIndex] = input[inputOffset + tid].first;
         } else {
            output[warpTidWriteIndex] = input[inputOffset + tid];
         }
      }
      // Next loop iteration:
      //input += current;
      inputOffset += current;
      output += outputCount;
      remaining -= current;
   }
   __syncthreads();
   if (tid == 0) {
      // Resize to final correct output size.
      outputVec->device_resize(outputSize);
      if (output_sizes) {// Only store lengths if output buffer is not null
         output_sizes[cellIndex] = outputSize;
      }
   }
}

template <typename Rule, typename ELEMENT, bool FIRSTONLY=false>
void extract_GIDs_kernel_launcher(
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** input_maps,
   split::SplitVector<ELEMENT> **output_vecs,
   vmesh::LocalID* output_sizes,
   Rule rule,
   vmesh::VelocityMesh** rule_meshes,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** rule_maps,
   split::SplitVector<vmesh::GlobalID>** rule_vectors,
   const uint nCells,
   gpuStream_t stream
   ) {
   extract_GIDs_kernel<Rule,ELEMENT,FIRSTONLY><<<nCells, Hashinator::defaults::MAX_BLOCKSIZE, 0, stream>>>(
      input_maps,
      output_vecs,
      output_sizes,
      rule,
      rule_meshes,
      rule_maps,
      rule_vectors
      );
   CHK_ERR( gpuPeekAtLastError() );
}

/*
 * Extracts key-value (GID-LID) pairs matching the given rule
 * from the hashmaps of all provided velocity meshes,
 * stores them in provided splitvectors, and
 * clears all tombstones and matched elements.
 */
template <typename Rule>
__global__ void __launch_bounds__(Hashinator::defaults::MAX_BLOCKSIZE, FULLBLOCKS_PER_MP) extract_overflown_kernel(
   vmesh::VelocityMesh **vmeshes, // buffer of pointers to vmeshes, contain hashmaps
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **output_vecs,
   vmesh::LocalID* output_sizes,
   Rule rule
   ) {
   //launch parameters: dim3 grid(nMaps,1,1); // As this is a looping reduction
   const size_t vmeshIndex = blockIdx.x;
   if (vmeshes[vmeshIndex]==0) {
      return; // Early return for invalid cells
   }
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* thisMap = vmeshes[vmeshIndex]->gpu_expose_map();
   Hashinator::Info *info = thisMap->expose_mapinfo<false>();
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> *outputVec = output_vecs[vmeshIndex];

   if (info->tombstoneCounter == 0) {
      // If there are no tombstones, then also any overflown elements will be minimally overflown.
      outputVec->device_resize(0);
      return;
   }
   // This must be equal to at least both WARPLENGTH and MAX_BLOCKSIZE/WARPLENGTH
   __shared__ uint32_t warpSums[WARPLENGTH];
   __shared__ uint32_t outputCount;
   const int tid = threadIdx.x;
   const int wid = tid / WARPLENGTH;
   const int w_tid = tid % WARPLENGTH;
   //const int warpsPerBlock = BLOCKSIZE / WARPLENGTH;
   const uint warpsPerBlock = blockDim.x / WARPLENGTH;
   // zero init shared buffer
   if (wid == 0) {
      warpSums[w_tid] = 0;
   }
   __syncthreads();
   // full warp votes for rule-> mask = [01010101010101010101010101010101]
   int64_t remaining = thisMap->bucket_count();
   const uint capacity = outputVec->capacity();
   uint32_t outputSize = 0;
   // Initial pointers into data
   Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID> *input = thisMap->expose_bucketdata<false>();
   Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>* output = outputVec->data();
   const vmesh::GlobalID emptybucket = thisMap->get_emptybucket();
   // Start loop
   while (remaining > 0) {
      int current = remaining > blockDim.x ? blockDim.x : remaining;
      __syncthreads();
      const int active = (tid < current) ? rule(thisMap, input[tid]) : false;
      const auto mask = split::s_warpVote(active == 1, SPLIT_VOTING_MASK);
      const auto warpCount = split::s_pop_count(mask);
      if (w_tid == 0) {
         warpSums[wid] = warpCount;
      }
      __syncthreads();
      // Figure out the total here because we overwrite shared mem later
      if (wid == 0) {
         // ceil int division
         int activeWARPS = nextPow2(1 + ((current - 1) / WARPLENGTH));
         auto reduceCounts = [activeWARPS](int localCount) -> int {
                                for (int i = activeWARPS / 2; i > 0; i = i / 2) {
                                   localCount += split::s_shuffle_down(localCount, i, SPLIT_VOTING_MASK);
                                }
                                return localCount;
                             };
         auto localCount = warpSums[w_tid];
         int totalCount = reduceCounts(localCount);
         if (w_tid == 0) {
            outputCount = totalCount;
            outputSize += totalCount;
            assert((outputSize <= capacity) && "extract_overflown_kernel ran out of capacity!");
            outputVec->device_resize(outputSize);
         }
      }
      // Prefix scan WarpSums on the first warp
      if (wid == 0) {
         auto value = warpSums[w_tid];
         for (uint d = 1; d < warpsPerBlock; d = 2 * d) {
            int res = split::s_shuffle_up(value, (int)d, SPLIT_VOTING_MASK);
            if (tid % warpsPerBlock >= d) {
               value += res;
            }
         }
         warpSums[w_tid] = value;
      }
      __syncthreads();
      auto offset = (wid == 0) ? 0 : warpSums[wid - 1];
      auto pp = split::s_pop_count(mask & ((ONE << w_tid) - ONE));
      const auto warpTidWriteIndex = offset + pp;
      if (active) {
         output[warpTidWriteIndex] = input[tid];
         // Now also delete this entry. Must edit fill count at end of kernel.
         input[tid].first = emptybucket;
      }
      // Next loop iteration:
      input += current;
      output += outputCount;
      remaining -= current;
   }
   __syncthreads();
   if (tid == 0) {
      // Resize to final correct output size.
      outputVec->device_resize(outputSize);
      output_sizes[vmeshIndex] = outputSize;
      // Update mapInfo
      info->currentMaxBucketOverflow = Hashinator::defaults::BUCKET_OVERFLOW;
      info->fill -= outputSize; // subtract deleted (overflown) elements
      info->tombstoneCounter = 0;
   }
}

template <typename Rule>
void clean_tombstones_launcher(
   vmesh::VelocityMesh** vmeshes,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **overflown_elements,
   vmesh::LocalID* output_sizes,
   Rule rule,
   const uint nCells,
   gpuStream_t stream
   ) {
   // Extract overflown elements into temporary vector
   extract_overflown_kernel<Rule><<<nCells, Hashinator::defaults::MAX_BLOCKSIZE, 0, stream>>>(
      vmeshes,
      overflown_elements,
      output_sizes,
      rule
      );
   CHK_ERR( gpuPeekAtLastError() );
}

/*
 * Mini-kernel for inserting previously extracted overflown elements
 */
__global__ void __launch_bounds__(GPUTHREADS, WARPS_PER_MP) batch_insert_kernel(
   vmesh::VelocityMesh **vmeshes, // buffer of pointers to vmeshes, contain hashmaps
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ const *input_vecs
   ) {
   //launch parameters: dim3 grid(largestOverflow,nCells,1);
   const uint ti = threadIdx.x; // [0,blockSize)
   const int b_tid = ti % GPUTHREADS; // [0,GPUTHREADS)
   // GPUTODO: several entries in parallel per block
   const size_t vmeshIndex = blockIdx.y;
   const size_t blockIndex = blockIdx.x;
   if (vmeshes[vmeshIndex]==0) {
      return; // Early return for invalid cells
   }
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* thisMap = vmeshes[vmeshIndex]->gpu_expose_map();
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ inputVec = input_vecs[vmeshIndex];

   const size_t inputVecSize = inputVec->size();
   if (inputVecSize == 0 || blockIndex >= inputVecSize) {
      // No elements to insert
      return;
   }

   #ifdef USE_BATCH_WARPACCESSORS
   // Insert into map only from threads 0...WARPSIZE
   if (b_tid < GPUTHREADS) {
      #ifdef DEBUG_SPATIAL_CELL
      thisMap->warpInsert((inputVec->at(blockIndex)).first,(inputVec->at(blockIndex)).second,b_tid);
      #else
      thisMap->warpInsert(((*inputVec)[blockIndex]).first,((*inputVec)[blockIndex]).second,b_tid);
      #endif

   }
   #else
   // Insert into map only from thread 0
   if (b_tid == 0) {
      #ifdef DEBUG_SPATIAL_CELL
      thisMap->set_element((inputVec->at(blockIndex)).first,(inputVec->at(blockIndex)).second);
      #else
      thisMap->set_element(((*inputVec)[blockIndex]).first,((*inputVec)[blockIndex]).second);
      #endif
   }
   #endif
}

#ifdef USE_BATCH_WARPACCESSORS
/** Gpu Kernel to quickly gather the v-space halo of local content blocks
    Halo of 1 in each direction adds up to 26 neighbours.
    For NVIDIA/CUDA, we dan do 26 neighbours and 32 threads per warp in a single block.
    For AMD/HIP, we dan do 13 neighbours and 64 threads per warp in a single block, meaning two loops per cell.
    In either case, we launch blocks equal to nCells * max_velocity_block_with_content_list_size
*/
__global__ void __launch_bounds__(26*32, FULLBLOCKS_PER_MP) batch_update_velocity_halo_kernel (
   const vmesh::VelocityMesh* __restrict__ const *vmeshes,
   const split::SplitVector<vmesh::GlobalID>* __restrict__ const *velocity_block_with_content_lists,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** allMaps
   ) {
   // launch grid dim3 grid(launchBlocks,nCells,1);
   // Each block manages a single GID at a time, all velocity neighbours
   const uint nCells = gridDim.y;
   const uint cellIndex = blockIdx.y;
   const uint blockiStart = blockIdx.x;
   //const int blockSize = blockDim.x; // should be 26*32 or 13*64
   const uint ti = threadIdx.x;

   // Cells such as DO_NOT_COMPUTE are identified with a zero in the vmeshes pointer buffer
   if (vmeshes[cellIndex] == 0) {
      return;
   }
   const vmesh::VelocityMesh* __restrict__ vmesh = vmeshes[cellIndex];
   const split::SplitVector<vmesh::GlobalID>* __restrict__ velocity_block_with_content_list = velocity_block_with_content_lists[cellIndex];
   const vmesh::GlobalID* __restrict__ velocity_block_with_content_list_data = velocity_block_with_content_list->data();
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwcl_map = allMaps[2*cellIndex];
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwncl_map = allMaps[2*cellIndex+1];
   const vmesh::LocalID nBlocks = velocity_block_with_content_list->size();

   const vmesh::LocalID blocki = blockiStart;
   {
      // Return if we are beyond the size of the list for this cell
      if (blocki >= nBlocks) {
         return;
      }
      // Which spatial neighbour to consider out of the 26 face, edge, or corner neighbors
      const int offsetIndex1 = ti / GPUTHREADS; // [0,26) (NVIDIA) or [0,13) (AMD)
      const int w_tid = ti % GPUTHREADS; // [0,WARPSIZE)

      // Assumes addWidthV = 1
      #ifdef __CUDACC__
      const int max_i=1;
      #endif
      #ifdef __HIP_PLATFORM_HCC___
      const int max_i=2;
      #endif
      for (int i=0; i<max_i; i++) {
         int offsetIndex = offsetIndex1 + 13*i;
         // nudge latter half in order to exclude self
         if (offsetIndex > 12) {
            offsetIndex++;
         }
         const int offset_vx = (offsetIndex % 3) - 1;
         const int offset_vy = ((offsetIndex / 3) % 3) - 1;
         const int offset_vz = (offsetIndex / 9) - 1;
         // Offsets verified in python
         #ifdef DEBUG_SPATIAL_CELL
         const vmesh::GlobalID GID = velocity_block_with_content_list->at(blocki);
         #else
         const vmesh::GlobalID GID = velocity_block_with_content_list_data[blocki];
         #endif
         vmesh::LocalID ind0,ind1,ind2;
         vmesh->getIndices(GID,ind0,ind1,ind2);
         const int nind0 = ind0 + offset_vx;
         const int nind1 = ind1 + offset_vy;
         const int nind2 = ind2 + offset_vz;
         const vmesh::GlobalID nGID
            = vmesh->getGlobalID(nind0,nind1,nind2);
         if (nGID != vmesh->invalidGlobalID()) {
            // Does block already exist in mesh?
            const vmesh::LocalID LID = vmesh->warpGetLocalID(nGID, w_tid);
            // Try adding this nGID to velocity_block_with_content_map. If it exists, do not overwrite.
            const bool newlyadded = vbwcl_map->warpInsert_V<true>(nGID,LID, w_tid);
            if (newlyadded) {
               // Block did not previously exist in velocity_block_with_content_map
               if ( LID != vmesh->invalidLocalID()) {
                  // Block exists in mesh, ensure it won't get deleted:
                  vbwncl_map->warpErase(nGID, w_tid);
               }
               // else:
               // Block does not yet exist in mesh at all. Needs adding!
               // Identified as invalidLID entries in velocity_block_with_content_map.
            }
         }
      }
      __syncthreads();
   }
}
#else // if not using warp accessors
/** Gpu Kernel to quickly gather the v-space halo of local content blocks
    Halo of 1 in each direction adds up to 26 neighbours.
    This kernel does not use warp accessors so always does all 26 neighbors in a single block.
*/
__global__ void batch_update_velocity_halo_kernel (
   const vmesh::VelocityMesh* __restrict__ const *vmeshes,
   const split::SplitVector<vmesh::GlobalID>* __restrict__ const *velocity_block_with_content_lists,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** allMaps,
   const uint warpsPerBlockBatchHalo
   ) {
   // launch grid dim3 grid(launchBlocks,nCells,1);
   // Each block manages a single GID at a time, all velocity neighbours
   //const uint nCells = gridDim.y;
   const uint cellIndex = blockIdx.y;
   const uint blockiStart = blockIdx.x*warpsPerBlockBatchHalo+threadIdx.y; // launch grid block index inside number of velocity blocks
   const uint ti = threadIdx.x; // Thread index inside warp / wavefront acting on single LID

   // Cells such as DO_NOT_COMPUTE are identified with a zero in the vmeshes pointer buffer
   if (vmeshes[cellIndex] == 0) {
      return;
   }
   // Only act on first 26 threads of each warp / wavefront
   if (ti >= 26) {
      return; // Note: this prevents use of syncthreads!
   }

   const vmesh::VelocityMesh* __restrict__ vmesh = vmeshes[cellIndex];
   const split::SplitVector<vmesh::GlobalID>* __restrict__ velocity_block_with_content_list = velocity_block_with_content_lists[cellIndex];
   const vmesh::GlobalID* __restrict__ velocity_block_with_content_list_data = velocity_block_with_content_list->data();
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwcl_map = allMaps[2*cellIndex];
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwncl_map = allMaps[2*cellIndex+1];
   const vmesh::LocalID nBlocks = velocity_block_with_content_list->size();

   const vmesh::LocalID blocki = blockiStart;
   {
      // Return if we are beyond the size of the list for this cell
      if (blocki >= nBlocks) {
         return; // Disallows use of __syncthreads() in this kernel
      }
      int offsetIndex = ti;
      // nudge latter half in order to exclude self
      if (offsetIndex > 12) {
         offsetIndex++;
      }
      const int offset_vx = (offsetIndex % 3) - 1;
      const int offset_vy = ((offsetIndex / 3) % 3) - 1;
      const int offset_vz = (offsetIndex / 9) - 1;
      // Offsets verified in python
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID GID = velocity_block_with_content_list->at(blocki);
      #else
      const vmesh::GlobalID GID = velocity_block_with_content_list_data[blocki];
      #endif
      vmesh::LocalID ind0,ind1,ind2;
      vmesh->getIndices(GID,ind0,ind1,ind2);
      const int nind0 = ind0 + offset_vx;
      const int nind1 = ind1 + offset_vy;
      const int nind2 = ind2 + offset_vz;
      const vmesh::GlobalID nGID
         = vmesh->getGlobalID(nind0,nind1,nind2);
      if (nGID != vmesh->invalidGlobalID()) {
         // Does block already exist in mesh?
         const vmesh::LocalID LID = vmesh->getLocalID(nGID);
         // Add this nGID to velocity_block_with_content_map.
         const bool newlyadded = vbwcl_map->set_element<true>(nGID,LID);
         if (newlyadded) {
            // Block did not previously exist in velocity_block_with_content_map
            if ( LID != vmesh->invalidLocalID()) {
               // Block exists in mesh, ensure it won't get deleted:
               vbwncl_map->device_erase(nGID);
            }
         }
      }
      //__syncthreads(); // Not allowed due to early thread returns
   }
}
#endif // end if warp accessors


#ifdef USE_BATCH_WARPACCESSORS
/** Gpu Kernel to quickly gather the spatial halo of neighbour content blocks
*/
__global__ void __launch_bounds__(GPUTHREADS*WARPSPERBLOCK, FULLBLOCKS_PER_MP) batch_update_neighbour_halo_kernel (
   const vmesh::VelocityMesh* __restrict__ const *vmeshes,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** allMaps,
   const split::SplitVector<vmesh::GlobalID>* __restrict__ const *neigh_velocity_block_with_content_lists
   ) {
   const uint nCells = gridDim.y;
   const uint maxNeighbours = gridDim.z;
   const uint cellIndex = blockIdx.y;
   const uint neighIndex = blockIdx.y * maxNeighbours + blockIdx.z;

   // Cells such as DO_NOT_COMPUTE are identified with a zero in the vmeshes pointer buffer
   if (vmeshes[cellIndex] == 0) {
      return;
   }
   // Early return for non-existing neighbour indexes
   if (neigh_velocity_block_with_content_lists[neighIndex] == 0) {
      return;
   }

   const int ti = threadIdx.x; // [0,blockSize)
   const int w_tid = ti % GPUTHREADS; // [0,WARPSIZE)
   const int w_id = ti / GPUTHREADS; // [0,WARPSPERBLOCK)

   const int blockWidth = WARPSPERBLOCK; // how many GIDs each GPU block manages at once (in parallel)
   const int blockiStart = blockIdx.x * blockWidth;

   const split::SplitVector<vmesh::GlobalID>* __restrict__ velocity_block_with_content_list = neigh_velocity_block_with_content_lists[neighIndex];
   const int nBlocks = velocity_block_with_content_list->size();

   for (int blocki = blockiStart + w_id; blocki < blockiStart+blockWidth; blocki += blockWidth) {
      // Skip to sync if we are beyond the size of the list for this cell
      if (blocki < nBlocks) {
         const vmesh::VelocityMesh* __restrict__ vmesh = vmeshes[cellIndex];
         const vmesh::GlobalID* __restrict__ velocity_block_with_content_list_data = velocity_block_with_content_list->data();
         Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwcl_map = allMaps[2*cellIndex];
         Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwncl_map = allMaps[2*cellIndex+1];

         #ifdef DEBUG_SPATIAL_CELL
         const vmesh::GlobalID nGID = velocity_block_with_content_list->at(blocki);
         #else
         const vmesh::GlobalID nGID = velocity_block_with_content_list_data[blocki];
         #endif
         // Does block already exist in mesh?
         const vmesh::LocalID LID = vmesh->warpGetLocalID(nGID, w_tid);
         // Try adding this nGID to velocity_block_with_content_map. If it exists, do not overwrite.
         const bool newlyadded = vbwcl_map->warpInsert_V<true>(nGID,LID, w_tid);
         if (newlyadded) {
            // Block did not previously exist in velocity_block_with_content_map
            if ( LID != vmesh->invalidLocalID()) {
               // Block exists in mesh, ensure it won't get deleted:
               vbwncl_map->warpErase(nGID, w_tid);
            }
            // else:
            // Block does not yet exist in mesh at all. Needs adding!
            // Identified as invalidLID entries in velocity_block_with_content_map.
         }
      }
      __syncthreads();
   }
}
#else // if not using warp accessors
/** Gpu Kernel to quickly gather the spatial halo of neighbour content blocks
*/
__global__ void __launch_bounds__(GPUTHREADS*WARPSPERBLOCK, FULLBLOCKS_PER_MP) batch_update_neighbour_halo_kernel (
   const vmesh::VelocityMesh* __restrict__ const *vmeshes,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** allMaps,
   const split::SplitVector<vmesh::GlobalID>* __restrict__ const *neigh_velocity_block_with_content_lists
   ) {

   //const uint nCells = gridDim.y;
   const uint maxNeighbours = gridDim.z;
   const uint cellIndex = blockIdx.y;
   const uint neighIndex = blockIdx.y * maxNeighbours + blockIdx.z;

   const vmesh::VelocityMesh* __restrict__ vmeshCellIndex = vmeshes[cellIndex];
   const split::SplitVector<vmesh::GlobalID>* __restrict__ velocity_block_with_content_list = neigh_velocity_block_with_content_lists[neighIndex];

   // Cells such as DO_NOT_COMPUTE are identified with a zero in the vmeshes pointer buffer
   if (vmeshCellIndex == 0) {
      return;
   }
   // Early return for non-existing neighbour indexes
   if (velocity_block_with_content_list == 0) {
      return;
   }

   const int blockWidth = blockDim.x; // how many GIDs each GPU block manages at once (in parallel)
   const int ti = threadIdx.x; // [0,blockSize)

   const int nBlocks = velocity_block_with_content_list->size();

   {
      const int blocki = blockIdx.x * blockWidth + ti;
      // Return if we are beyond the size of the list for this cell
      if (blocki >= nBlocks) {
         return; // Disallows use of __syncthreads() in this kernel
      }
      const vmesh::GlobalID* __restrict__ velocity_block_with_content_list_data = velocity_block_with_content_list->data();
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwcl_map = allMaps[2*cellIndex];
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwncl_map = allMaps[2*cellIndex+1];

      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID nGID = velocity_block_with_content_list->at(blocki);
      #else
      const vmesh::GlobalID nGID = velocity_block_with_content_list_data[blocki];
      #endif
      // Does block already exist in mesh?
      const vmesh::LocalID LID = vmeshCellIndex->getLocalID(nGID);
      // Add this nGID to velocity_block_with_content_map.
      const bool newlyadded = vbwcl_map->set_element<true>(nGID,LID);
      if (newlyadded) {
         // Block did not previously exist in velocity_block_with_content_map
         if ( LID != vmeshCellIndex->invalidLocalID()) {
            // Block exists in mesh, ensure it won't get deleted:
            vbwncl_map->device_erase(nGID);
         }
      }
      //__syncthreads(); // Not allowed due to early thread returns
   }
}
#endif

/** Mini-kernel for checking list sizes and attempting to adjust vmesh and VBC size on-device */
__global__ void batch_resize_vbc_kernel_pre(
   vmesh::VelocityMesh **vmeshes,
   vmesh::VelocityBlockContainer **blockContainers,
   split::SplitVector<vmesh::GlobalID>** dev_list_with_replace_new,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>** dev_list_delete,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>** dev_list_to_replace,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>** dev_list_with_replace_old,
   // return values: nbefore, nafter, nblockstochange, resize success
   vmesh::LocalID* dev_nBefore,
   vmesh::LocalID* dev_nAfter,
   vmesh::LocalID* dev_nBlocksToChange,
   vmesh::LocalID* dev_resizeSuccess,
   Real* dev_rhoLossAdjust // mass loss, set to zero
   ) {
   const size_t cellIndex = blockIdx.x;
   if (vmeshes[cellIndex]==0) {
      return; // Early return for invalid cells
   }
   vmesh::VelocityMesh *vmesh = vmeshes[cellIndex];
   vmesh::VelocityBlockContainer *blockContainer = blockContainers[cellIndex];
   split::SplitVector<vmesh::GlobalID>* list_with_replace_new = dev_list_with_replace_new[cellIndex];
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_delete = dev_list_delete[cellIndex];
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_to_replace = dev_list_to_replace[cellIndex];
   //split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_with_replace_old = dev_list_with_replace_old[cellIndex];

   const vmesh::LocalID nBlocksBeforeAdjust = vmesh->size();
   const vmesh::LocalID nToAdd = list_with_replace_new->size();
   const vmesh::LocalID nToRemove = list_delete->size() + list_to_replace->size();
   const vmesh::LocalID nBlocksAfterAdjust = nBlocksBeforeAdjust + nToAdd - nToRemove;
   const vmesh::LocalID nBlocksToChange = nToAdd > nToRemove ? nToAdd : nToRemove;

   dev_rhoLossAdjust[cellIndex] = 0.0;
   dev_nBefore[cellIndex] = nBlocksBeforeAdjust;
   dev_nAfter[cellIndex] = nBlocksAfterAdjust;
   dev_nBlocksToChange[cellIndex] = nBlocksToChange;
   // Should we grow the size?
   if (nBlocksAfterAdjust > nBlocksBeforeAdjust) {
      if ((nBlocksAfterAdjust <= vmesh->capacity()) && (nBlocksAfterAdjust <= blockContainer->capacity())) {
         dev_resizeSuccess[cellIndex] = 1; // Resize on-device will work.
         vmesh->device_setNewSize(nBlocksAfterAdjust);
         blockContainer->setNewSize(nBlocksAfterAdjust);
      } else {
         dev_resizeSuccess[cellIndex] = 0; // Need to recapacitate and resize from host
      }
   } else {
      // No error as no resize.
      dev_resizeSuccess[cellIndex] = 2;
   }
}

/** Mini-kernel for adjusting vmesh and VBC size on-device aftewards (shrink only) */
__global__ void batch_resize_vbc_kernel_post(
   vmesh::VelocityMesh **vmeshes,
   vmesh::VelocityBlockContainer **blockContainers,
   vmesh::LocalID* dev_nAfter
   ) {
   const size_t cellIndex = blockIdx.x;
   if (vmeshes[cellIndex]==0) {
      return; // Early return for invalid cells
   }
   vmesh::VelocityMesh *vmesh = vmeshes[cellIndex];
   vmesh::VelocityBlockContainer *blockContainer = blockContainers[cellIndex];
   const vmesh::LocalID nBlocksAfterAdjust = dev_nAfter[cellIndex];
   vmesh->device_setNewSize(nBlocksAfterAdjust);
   blockContainer->setNewSize(nBlocksAfterAdjust);
}


/** GPU kernel for updating blocks based on generated lists */
__global__ void __launch_bounds__(WID3, WID3S_PER_MP) batch_update_velocity_blocks_kernel(
   vmesh::VelocityMesh **vmeshes,
   vmesh::VelocityBlockContainer **blockContainers,
   const split::SplitVector<vmesh::GlobalID>* __restrict__ const *dev_list_with_replace_new,
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ const *dev_list_delete,
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ const *dev_list_to_replace,
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ const *dev_list_with_replace_old,
   vmesh::LocalID* dev_nBefore,
   vmesh::LocalID* dev_nAfter,
   vmesh::LocalID* dev_nBlocksToChange,
   Real* dev_rhoLossAdjust // mass loss, gather from deleted blocks
   ) {
   // launch griddim3 grid(launchBlocks,nCells,1);
   const size_t cellIndex = blockIdx.y;
   if (vmeshes[cellIndex]==0) {
      return; // Early return for invalid cells
   }
   vmesh::VelocityMesh *vmesh = vmeshes[cellIndex];
   vmesh::VelocityBlockContainer *blockContainer = blockContainers[cellIndex];
   const split::SplitVector<vmesh::GlobalID>* __restrict__ list_with_replace_new = dev_list_with_replace_new[cellIndex];
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ list_delete = dev_list_delete[cellIndex];
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ list_to_replace = dev_list_to_replace[cellIndex];
   const split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* __restrict__ list_with_replace_old = dev_list_with_replace_old[cellIndex];

   const vmesh::LocalID nBlocksBeforeAdjust = dev_nBefore[cellIndex];
   const vmesh::LocalID nBlocksAfterAdjust = dev_nAfter[cellIndex];
   const vmesh::LocalID nBlocksToChange = dev_nBlocksToChange[cellIndex];

   if (blockIdx.x >= nBlocksToChange) {
      return; // Early return if outside list of blocks to change
   }
   const uint ti = threadIdx.x; // [0,blockSize)

   // This index into vectors can be adjusted along the way
   uint index = (uint)blockIdx.x;

   const int b_tid = ti % WID3; // [0,WID3)

   const vmesh::LocalID n_with_replace_new = list_with_replace_new->size();
   const vmesh::LocalID n_delete = list_delete->size();
   const vmesh::LocalID n_to_replace = list_to_replace->size();
   const vmesh::LocalID n_with_replace_old = list_with_replace_old->size();
   // For tracking mass-loss
   __shared__ Real massloss[WID3];

   // Each block Processes one block from the lists.

   /*********
       Check if should delete item from end of vmesh.
       For this, we get both GID and LID from the vector.
   **/
   if (index < n_delete) {
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID rmGID = (list_delete->at(index)).first;
      const vmesh::GlobalID rmLID = (list_delete->at(index)).second;
      #else
      const vmesh::GlobalID rmGID = ((*list_delete)[index]).first;
      const vmesh::GlobalID rmLID = ((*list_delete)[index]).second;
      #endif

      #ifdef DEBUG_SPATIAL_CELL
      if (rmGID == vmesh->invalidGlobalID()) {
         if (rmLID != vmesh->invalidLocalID()) {
            // Valid LID but invalid GID: only remove from vmesh localToGlobal?
            if (b_tid==0) {
               printf("Removing blocks: Valid LID %u but invalid GID!\n",rmLID);
            }
         } else {
            if (b_tid==0) {
               printf("Removing blocks: Invalid LID and GID!\n");
            }
         }
         assert(0);
      }
      if (rmLID == vmesh->invalidLocalID()) {
         if (rmGID != vmesh->invalidGlobalID()) {
            // Valid GID but invalid LID: only remove from vmesh globalToLocal?
            if (b_tid==0) {
               printf("Removing blocks: Valid GID %ul but invalid LID!\n",rmGID);
            }
         }
         assert(0);
      }
      if ((unsigned long)rmLID >= (unsigned long)nBlocksBeforeAdjust) {
         if (b_tid==0) {
            printf("Trying to outright remove block which has LID %ul >= nBlocksBeforeAdjust %ul!\n",rmLID,nBlocksBeforeAdjust);
         }
         assert(0);
      }
      if ((unsigned long)rmLID < (unsigned long)nBlocksAfterAdjust) {
         if (b_tid==0) {
            printf("Trying to outright remove block which has LID %u smaller than nBlocksAfterAdjust %u!\n",rmLID,nBlocksAfterAdjust);
         }
         assert(0);
      }
      #endif

      // Track mass loss:
      Realf* rm_avgs = blockContainer->getData(rmLID);
      Real* rm_block_parameters = blockContainer->getParameters(rmLID);
      const Real rm_DV3 = rm_block_parameters[BlockParams::DVX]
         * rm_block_parameters[BlockParams::DVY]
         * rm_block_parameters[BlockParams::DVZ];
      // thread-sum for rho
      massloss[ti] = rm_avgs[b_tid]*rm_DV3;
      __syncthreads();
      // Implemented just a simple non-optimized thread sum
      for (int s=WID3/2; s>0; s>>=1) {
         if (b_tid < s) {
            massloss[ti] += massloss[ti + s];
         }
         __syncthreads();
      }
      // Bookkeeping only by one thread per block
      if (b_tid==0) {
         Real old = atomicAdd(&dev_rhoLossAdjust[cellIndex], massloss[ti]);
      }
      __syncthreads();

      // Delete from vmesh
      #ifdef USE_BATCH_WARPACCESSORS
      vmesh->warpDeleteBlock(rmGID,rmLID,b_tid);
      #else
      if (b_tid==0) {
         vmesh->deleteBlock(rmGID,rmLID);
      }
      #endif
      // GPUTODO debug checks
      return;
   }
   index -= n_delete;

   /*********
       Check if should replace existing block with either
       existing block from end of vmesh or new block
   **/
   if (index < n_to_replace) {
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID rmGID = (list_to_replace->at(index)).first;
      const vmesh::GlobalID rmLID = (list_to_replace->at(index)).second;
      #else
      const vmesh::GlobalID rmGID = ((*list_to_replace)[index]).first;
      const vmesh::GlobalID rmLID = ((*list_to_replace)[index]).second;
      #endif
      //const vmesh::LocalID rmLID = vmesh->warpGetLocalID(rmGID,b_tid);

      #ifdef DEBUG_SPATIAL_CELL
      if (rmGID == vmesh->invalidGlobalID()) {
         if (rmLID != vmesh->invalidLocalID()) {
            // Valid LID but invalid GID: only remove from vmesh localToGlobal?
            if (b_tid==0) {
               printf("Replacing blocks: Valid LID %u but invalid GID!\n",rmLID);
            }
         } else {
            if (b_tid==0) {
               printf("Replacing blocks: Invalid LID and GID!\n");
            }
         }
         assert(0);
      }
      if (rmLID == vmesh->invalidLocalID()) {
         if (rmGID != vmesh->invalidGlobalID()) {
            // Valid GID but invalid LID: only remove from vmesh globalToLocal?
            if (b_tid==0) {
               printf("Replacing blocks: Valid GID %ul but invalid LID!\n",rmGID);
            }
         }
         assert(0);
      }
      if (rmLID >= nBlocksBeforeAdjust) {
         if (b_tid==0) {
            printf("Trying to replace block which has LID %ul >= nBlocksBeforeAdjust %ul!\n",rmLID,nBlocksBeforeAdjust);
         }
         assert(0);
      }
      #endif

      // Track mass loss:
      Realf* rm_avgs = blockContainer->getData(rmLID);
      Real* rm_block_parameters = blockContainer->getParameters(rmLID);
      const Real rm_DV3 = rm_block_parameters[BlockParams::DVX]
         * rm_block_parameters[BlockParams::DVY]
         * rm_block_parameters[BlockParams::DVZ];
      // thread-sum for rho
      massloss[ti] = rm_avgs[b_tid]*rm_DV3;
      __syncthreads();
      // Implemented just a simple non-optimized thread sum
      for (int s=WID3/2; s>0; s>>=1) {
         if (b_tid < s) {
            massloss[ti] += massloss[ti + s];
         }
         __syncthreads();
      }
      // Bookkeeping only by one thread per block
      if (b_tid==0) {
         Real old = atomicAdd(&dev_rhoLossAdjust[cellIndex], massloss[ti]);
      }
      __syncthreads();

      // Figure out what to use as replacement
      vmesh::GlobalID replaceGID;
      vmesh::LocalID replaceLID;

      // First option: replace with existing block from end of vmesh
      if (index < n_with_replace_old) {
         #ifdef DEBUG_SPATIAL_CELL
         replaceGID = (list_with_replace_old->at(index)).first;
         replaceLID = (list_with_replace_old->at(index)).second;
         #else
         replaceGID = ((*list_with_replace_old)[index]).first;
         replaceLID = ((*list_with_replace_old)[index]).second;
         #endif

         Realf* repl_avgs = blockContainer->getData(replaceLID);
         Real*  repl_block_parameters = blockContainer->getParameters(replaceLID);
         rm_avgs[b_tid] = repl_avgs[b_tid];
         if (b_tid < BlockParams::N_VELOCITY_BLOCK_PARAMS) {
            rm_block_parameters[b_tid] = repl_block_parameters[b_tid];
         }
         __syncthreads();

      } else {
         // Second option: add new block instead
         #ifdef DEBUG_SPATIAL_CELL
         replaceGID = list_with_replace_new->at(index - n_with_replace_old);
         #else
         replaceGID = (*list_with_replace_new)[index - n_with_replace_old];
         #endif
         replaceLID = vmesh->invalidLocalID();

         rm_avgs[b_tid] = 0;
         if (b_tid==0) {
            // Write in block parameters
            vmesh->getBlockInfo(replaceGID, rm_block_parameters+BlockParams::VXCRD);
         }
         __syncthreads();
      }
      // Remove hashmap entry for removed block, add instead created block
      #ifdef USE_BATCH_WARPACCESSORS
      vmesh->warpReplaceBlock(rmGID,rmLID,replaceGID,b_tid);
      #else
      if (b_tid==0) {
         vmesh->replaceBlock(rmGID,rmLID,replaceGID);
      }
      #endif
      #ifdef DEBUG_SPATIAL_CELL
      __syncthreads();
      if (vmesh->getGlobalID(rmLID) != replaceGID) {
         if (b_tid==0) {
            printf("Error! Replacing did not result in wanted GID at old LID in update_velocity_blocks_kernel! \n");
         }
         assert(0);
      }
      if (vmesh->getLocalID(replaceGID) != rmLID) {
         if (b_tid==0) {
            printf("Error! Replacing did not result in old LID at replaced GID in update_velocity_blocks_kernel! \n");
         }
         assert(0);
      }
      #endif

      return;
   }
   index -= n_to_replace;

   /*********
       Finally check if we should add new block after end of current vmesh
       We have reserved/used some entries from the beginning of the list_with_replace_new
       for the previous section, so now we access that with a different index.
   **/
   const uint add_index = index + (n_to_replace - n_with_replace_old);
   if (add_index < n_with_replace_new) {
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::GlobalID addGID = list_with_replace_new->at(add_index);
      if (vmesh->getLocalID(addGID) != vmesh->invalidLocalID()) {
         if (b_tid==0) {
            printf("Trying to add new GID %u to mesh which already contains it! index=%u addindex=%u\n",addGID,index,add_index);
         }
         assert(0);
      }
      __syncthreads();
      #else
      const vmesh::GlobalID addGID = (*list_with_replace_new)[add_index];
      #endif

      // We need to add the data of addGID to a new LID. Here we still use the regular index.
      const vmesh::LocalID addLID = nBlocksBeforeAdjust + index;
      Realf* add_avgs = blockContainer->getData(addLID);
      #ifdef DEBUG_SPATIAL_CELL
      if (addGID == vmesh->invalidGlobalID()) {
         printf("Error! invalid addGID!\n");
         assert(0);
      }
      if (addLID == vmesh->invalidLocalID()) {
         printf("Error! invalid addLID!\n");
         assert(0);
      }
      #endif
      Real* add_block_parameters = blockContainer->getParameters(addLID);
      // Zero out blockdata
      add_avgs[b_tid] = 0;
      if (b_tid==0) {
         // Write in block parameters
         vmesh->getBlockInfo(addGID, add_block_parameters+BlockParams::VXCRD);
      }
      __syncthreads();

      // Insert new hashmap entry into vmesh
      #ifdef USE_BATCH_WARPACCESSORS
      vmesh->warpPlaceBlock(addGID,addLID,b_tid);
      #else
      if (b_tid==0) {
         vmesh->placeBlock(addGID,addLID);
      }
      #endif
      #ifdef DEBUG_SPATIAL_CELL
      __syncthreads();
      if (vmesh->getGlobalID(addLID) == vmesh->invalidGlobalID()) {
         printf("Error! invalid GID after add from addLID!\n");
         assert(0);
      }
      if (vmesh->getLocalID(addGID) == vmesh->invalidLocalID()) {
         printf("Error! invalid LID after add from addGID!\n");
         assert(0);
      }
      #endif
      return;
   }

   // Fall-through error!
   if (b_tid==0) {
      printf("Error! Fall through in batch_update_velocity_blocks_kernel! index %u nBlocksBeforeAdjust %u nBlocksAfterAdjust %u \n",
             index,nBlocksBeforeAdjust,nBlocksAfterAdjust);
   }
   __syncthreads();
}

/** GPU kernel for batch-scaling particle populations
 */
__global__ void __launch_bounds__(WID3, WID3S_PER_MP) batch_population_scale_kernel (
   vmesh::VelocityBlockContainer **blockContainers,
   Real* dev_mass_scale
   ) {
   // launch griddim3 grid(launchBlocks,nCells,1);
   const int cellIndex = blockIdx.y;
   const int blocki = blockIdx.x;
   const uint ti = threadIdx.x;

   vmesh::VelocityBlockContainer* blockContainer = blockContainers[cellIndex];
   const Real cell_mass_scale = dev_mass_scale[cellIndex];

   const uint b_tid = ti % WID3; // [0,WID3)
   const uint blockLID = blocki; // [0,nBlocksToChange)

   const uint VBC_size = blockContainer->size();
   if (blockLID > VBC_size || cell_mass_scale <= 0) {
      return;
   }
   // Pointer to target block data
   Realf* data = blockContainer->getData(blockLID);
   // Scale value
   data[b_tid] = data[b_tid] * cell_mass_scale;
}

#endif
