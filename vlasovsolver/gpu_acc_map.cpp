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

#include "gpu_acc_map.hpp"
#include "../spatial_cells/block_adjust_gpu.hpp"


/* These macros are used for bank collision reduction on CUDA hardware in the
   scan_probe kernel.
   NUM_BANKS and LOG_NUM_BANKS are defined in splitvector headers (32 and 5)
   TODO: Which one provides best bank conflict avoidance? Depends on hardware?
   On some hardware this gives the warning #63-D: shift count is too large, yet works.
*/
#define LOG_BANKS 4
#ifdef USE_CUDA // Nvidia hardware
#define BANK_OFFSET(n)                        \
  ((n) >> (LOG_BANKS) + (n) >> (2 * LOG_BANKS))
//#define BANK_OFFSET(n) ((n) >> LOG_BANKS) // segfaults, do not use
#else // AMD hardware
#define BANK_OFFSET(n) 0 // Reduces to no bank conflict elimination
#endif

/*!
  \brief GPU kernel which fills the target probe cube with the invalid value for vmesh::LocalID

   Also clears provided vectors which are needed as empty for future kernels. Clearing on-device
   prevents page faulting of unified memory bookkeeping data.
   Resizes the per-cell velocity blocks with content list vector as it'll be re-used for a LID list.

   @param vmeshes pointer to buffer of pointers to vmeshes of all cells on this process
   @param dev_probeCubeData pointer to buffer of pointers, which point to allocated
    temporary arrays, which are recast and used as the probe cube.
   @param flatExtent the product of the non-accelerated dimensions of the probe cube, i.e. the
    extent of it when it's flattend, rounded up to provide memory alignment.
   @param nTot number of total entries in probe cube
   @param invalidLID value to be used as the value for invalid vmesh::LocalID
   @param lists_with_replace_new pointer to buffer of pointers to SplitVectors to be cleared
   @param lists_delete pointer to buffer of pointers to SplitVectors to be cleared
   @param lists_to_replace pointer to buffer of pointers to SplitVectors to be cleared
   @param lists_with_replace_old pointer to buffer of pointers to SplitVectors to be cleared
   @param dev_vbwcl_vec pointer to buffer of pointers to SplitVectors, re-sized and re-cast to use as a LID list
   @param cumulativeOffset the current cumulative offset at which to index the aforementioned buffers, using
    the grid index of this kernel on top of this provided offset.
*/
__global__ void prefill_probe_kernel(
   vmesh::VelocityMesh** __restrict__ vmeshes,
   vmesh::LocalID *dev_probeCubeData,
   const uint flatExtent,
   const size_t Dacc,
   const size_t Dother,
   const vmesh::LocalID invalidLID,
   // Pass these for emptying
   split::SplitVector<vmesh::GlobalID>* *lists_with_replace_new,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* *lists_delete,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* *lists_to_replace,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* *lists_with_replace_old,
   // This one is resized and re-used as a LIDlist
   split::SplitVector<vmesh::GlobalID> ** dev_vbwcl_vec,
   const uint cumulativeOffset,
   const size_t gpu_probeStride
   ) {
   const size_t ind = blockIdx.x * blockDim.x + threadIdx.x;
   const uint parallelOffsetIndex = blockIdx.y;
   const uint cellOffset = parallelOffsetIndex + cumulativeOffset;
   const size_t nTot = Dacc*Dother;
   vmesh::LocalID* probeFlattened = dev_probeCubeData + parallelOffsetIndex * gpu_probeStride;
   //vmesh::LocalID* probeCube = probeFlattened + flatExtent*GPU_PROBEFLAT_N;

   if (ind < flatExtent*GPU_PROBEFLAT_N) {
      // Flattened probe region
      probeFlattened[ind] = 0;
   } else if (ind < flatExtent*GPU_PROBEFLAT_N + nTot) {
      // Probe cube region
      probeFlattened[ind] = invalidLID;
   }
   // Device clears from single thread per cell
   if (ind==0) {
      lists_with_replace_new[cellOffset]->clear();
      lists_delete[cellOffset]->clear();
      lists_to_replace[cellOffset]->clear();
      lists_with_replace_old[cellOffset]->clear();
      const vmesh::VelocityMesh* __restrict__ vmesh = vmeshes[cellOffset];
      const vmesh::LocalID nBlocks = vmesh->size();
      dev_vbwcl_vec[cellOffset]->device_resize(nBlocks,false); // false: do not construct / reset new entries
   }
}

/*!
   GPU kernel which fills the target block container with zeroes

   @param blockContainers pointer to buffer of pointers to per-cell velocity block containers
   @param cumulativeOffset the current cumulative offset at which to index the aforementioned buffers, using
    the grid index of this kernel on top of this provided offset.
*/
__global__ void fill_VBC_zero_kernel(
   vmesh::VelocityBlockContainer** blockContainers,
   const uint cumulativeOffset
   ) {
   const uint parallelOffsetIndex = blockIdx.y;
   const uint cellOffset = parallelOffsetIndex + cumulativeOffset;
   vmesh::VelocityBlockContainer* blockContainer = blockContainers[cellOffset];
   const size_t VBC_size = blockContainer->size() * WID3;
   Realf *blockData = blockContainer->getData();

   const size_t ind = blockIdx.x * blockDim.x + threadIdx.x;
   for (size_t i = ind; i < VBC_size; i += gridDim.x * blockDim.x) {
      blockData[i] = 0;
   }
}

/*!
   \brief GPU kernel for taking the contents of a vmesh and placign the existing
   velocity blocks in a probe cube for further reduction by a flattening kernel.

   One option for a probe cube would be to reduce (with __ballot_sync) along the
   direction of propagation (k) to get the number of blocks in a column. However, this
   does not have an obvious way to support gathering columnsets (several columns at
   one set of perpendicular i,j indices).

   Thus, instead, for following analysis of the probe cube, we want the warp/wavefront
   to read a dimension *not* propagating along, because then the wavefront can loop
   over the dimension to propagate along (size Dacc), and each thread processes
   one potential columnSet. To even better parallelize and simplify, the two
   non-propagated dimensions (i,j) are merged int a dimension of size Dother
   (so it isn't an actual cube).

   Since we read the data (existing blocks) in LID order, writes will be jumbled
   anyhow, so we can take this opportunity to make future
   accesses to the probe cube efficient by writing in a smart order.

   TODO: ensure the first and second dimensions are powers of two for
   optimized reads? Then stepping will be based on array edge sizes instead of D0/1/2.

   @param vmeshes pointer to buffer of pointers to vmeshes of all cells on this process
   @param dev_probeCubeData pointer to buffer of pointers, which point to allocated
    temporary arrays, which are recast and used as the probe cube.
   @param flatExtent the product of the non-accelerated dimensions of the probe cube, i.e. the
    extent of it when it's flattend, rounded up to provide memory alignment.
   @param gpu_block_indices_to_probe 3-element array used for converting block GID to location
    within probe cube
   @param cumulativeOffset the current cumulative offset at which to index the aforementioned buffers, using
    the grid index of this kernel on top of this provided offset.
 */
__global__ void fill_probe_ordered(
   vmesh::VelocityMesh** __restrict__ vmeshes,
   vmesh::LocalID *dev_probeCubeData, // recast to vmesh::LocalID *probeCube
   const uint flatExtent,
   const uint* __restrict__ gpu_block_indices_to_probe,
   const uint cumulativeOffset,
   const size_t gpu_probeStride
   ) {
   const int ti = threadIdx.x; // [0,Hashinator::defaults::MAX_BLOCKSIZE)
   const vmesh::LocalID LID = blockDim.x * blockIdx.x + ti;
   const uint parallelOffsetIndex = blockIdx.y;
   const uint cellOffset = parallelOffsetIndex + cumulativeOffset;
   const vmesh::VelocityMesh* __restrict__ vmesh = vmeshes[cellOffset];
   const vmesh::LocalID nBlocks = vmesh->size();

   if (LID >= nBlocks) {
      return;
   }
   vmesh::LocalID* probeFlattened = dev_probeCubeData + parallelOffsetIndex * gpu_probeStride;
   vmesh::LocalID* probeCube = probeFlattened + flatExtent*GPU_PROBEFLAT_N;
   // Store in probe cube with ordering so that reading will be fast
   const vmesh::GlobalID GID = vmesh->getGlobalID(LID);
   vmesh::LocalID indices[3];
   vmesh->getIndices(GID,indices[0],indices[1],indices[2]);

   // Use pre-calculated probe indices
   const int target = indices[0] * gpu_block_indices_to_probe[0]
      + indices[1] * gpu_block_indices_to_probe[1]
      + indices[2] * gpu_block_indices_to_probe[2];

   probeCube[target] = LID;
}

/*!
   \brief GPU kernel which flattens the probe cube into two reduction results
   (counters): how many columns and how many blocks exist per columnset.

   In the probe cube, there's the dimension of acceleration (size Dacc)
   and the other two dimensions, merged (size Dother).
   For each index in the other two dimensions (i,j), we have a position in
   transcerse v-space associated with the possibility for constructing columns.
   Thus, that position is now termed a potential column position or potColumn.

   This kernel loops over Dacc to find:
   (1) How many acceleration columns were found for each potColumn
   (2) How many blocks were found for each potColumn

   @param dev_probeCubeData pointer to buffer of pointers, which point to allocated
    temporary arrays, which are recast and used as the probe cube.
   @param Dacc Extent of probe cube in accelerated dimension
   @param Dother Sum of extents of probe cube in transverse dimensions
   @param flatExtent the product of the non-accelerated dimensions of the probe cube, i.e. the
    extent of it when it's flattend, rounded up to provide memory alignment.
   @param invalidLID value to be used as the value for invalid vmesh::LocalID
*/
__global__ void flatten_probe_cube(
   vmesh::LocalID *dev_probeCubeData, // recast to vmesh::LocalID *probeCube, *probeFlattened
   const vmesh::LocalID Dacc,
   const vmesh::LocalID Dother,
   const size_t flatExtent,
   const vmesh::LocalID invalidLID,
   const size_t gpu_probeStride
   ) {
   // Probe cube contents have been ordered based on acceleration dimesion
   // so this kernel always reads in the same way.

   const int ti = threadIdx.x; // [0,Hashinator::defaults::MAX_BLOCKSIZE)
   const vmesh::LocalID ind = blockDim.x * blockIdx.x + ti;
   const uint parallelOffsetIndex = blockIdx.y;

   vmesh::LocalID* probeFlattened = dev_probeCubeData + parallelOffsetIndex * gpu_probeStride;
   vmesh::LocalID* probeCube = probeFlattened + flatExtent*GPU_PROBEFLAT_N;

   if (ind < Dother) {
      // Per-thread counters
      vmesh::LocalID foundBlocks = 0;
      vmesh::LocalID foundCols = 0;
      bool inCol = false;

      for (vmesh::LocalID j = 0; j < Dacc; j++) {
         if (probeCube[j*Dother + ind] == invalidLID) {
            // No block at this index.
            if (inCol) {
               // finish current column
               foundCols++;
               inCol = false;
            }
         } else {
            // Valid block found at this index
            foundBlocks++;
            if (!inCol) {
               // start new column
               inCol = true;
            }
         }
      }
      // Finished loop. If we are "still in a colum", count that.
      if (inCol) {
         foundCols++;
      }
      // Store values in global memory array
      probeFlattened[ind] = foundCols;
      probeFlattened[flatExtent + ind] = foundBlocks;
   }
}

/*!
   \brief GPU kernel which performs exclusive prefix scans of the flattened probe cube,
   providing cumulative sums up to each index for use in further kernels.

   This scan produces:

   (1) the cumulative sum of columns up to the beginning of each potColumn
   (2) the cumulative sum of column sets up to the beginning of each potColumn
   (3) the cumulative sum of blocks up to the beginning of each potColumn

   This is not a fully optimized scan as it uses only a single block
   in order to perform it all in one kernel. As the fixed-size shared memory
   buffer gets overwritten each cycle, we store the actual prefix scan results
   into the third and fourth entries in the probeFlattened buffer and keep
   track of the accumulated offset to the prefix.

   The count of columns to be evaluated is estimated to remain somewhat small-ish,
   so there should not be all that many loops of the cycle to deal with.

   @param vmeshes pointer to buffer of pointers to vmeshes of all cells on this process
   @param dev_probeCubeData pointer to buffer of pointers, which point to allocated
    temporary arrays, which are recast and used as the probe cube.
   @param Dacc Extent of probe cube in accelerated dimension
   @param Dother Sum of extents of probe cube in transverse dimensions
   @param flatExtent the product of the non-accelerated dimensions of the probe cube, i.e. the
    extent of it when it's flattend, rounded up to provide memory alignment.
   @param dev_numCols buffer for storing the total count of columns per cell (or probe cube)
   @param dev_numColSets buffer for storing the total count of columnSets per cell (or probe cube)
   @param dev_resizeSuccess buffer for storing a flag, if the columnData container re-sizing to
    match the counts of columns and columnSets was suffessful or not. If failed, re-sizing needs
    to occur on host.
   @param dev_columnOffsetData pointer to a buffer of ColumnOffsets structs (see arch/gpu_base.hpp) where
    per-cell information on all built columns and columnSets is stored
   @param cumulativeOffset the current cumulative offset at which to index the aforementioned buffers, using
    the grid index of this kernel on top of this provided offset.
*/

__global__ void scan_probe(
   vmesh::VelocityMesh** __restrict__ vmeshes,
   vmesh::LocalID *dev_probeCubeData, // recast to vmesh::LocalID *probeFlattened
   const vmesh::LocalID Dacc,
   const vmesh::LocalID Dother,
   const size_t flatExtent,
   vmesh::LocalID *dev_numCols,
   vmesh::LocalID *dev_numColSets,
   vmesh::LocalID *dev_resizeSuccess,
   ColumnOffsets* dev_columnOffsetData,
   const uint cumulativeOffset,
   const size_t gpu_probeStride
   ) {
   const uint parallelOffsetIndex = blockIdx.y;
   const uint cellOffset = parallelOffsetIndex + cumulativeOffset;
   vmesh::LocalID* probeFlattened = dev_probeCubeData + parallelOffsetIndex * gpu_probeStride;
   //vmesh::LocalID* probeCube = probeFlattened + flatExtent*GPU_PROBEFLAT_N;

   const vmesh::VelocityMesh* __restrict__ vmesh = vmeshes[cellOffset];
   const vmesh::LocalID nBlocks = vmesh->size(); // For early exit

   // Per-thread counters in shared memory for reduction. Double size buffer for better bank conflict avoidance.
   const int n = 2*Hashinator::defaults::MAX_BLOCKSIZE;
   __shared__ vmesh::LocalID reductionA[2*Hashinator::defaults::MAX_BLOCKSIZE]; // columns
   __shared__ vmesh::LocalID reductionB[2*Hashinator::defaults::MAX_BLOCKSIZE]; // columnsets
   __shared__ vmesh::LocalID reductionC[2*Hashinator::defaults::MAX_BLOCKSIZE]; // blocks
   __shared__ vmesh::LocalID offsetA;
   __shared__ vmesh::LocalID offsetB;
   __shared__ vmesh::LocalID offsetC;

   const int ti = threadIdx.x;
   if (ti==0) { // Cumulative result gathered per cycle
      offsetA = 0;
      offsetB = 0;
      offsetC = 0;
   }
   __syncthreads();
   ColumnOffsets* columnData = dev_columnOffsetData + parallelOffsetIndex;
   size_t majorOffset = 0;
   // Utilizes bank conflict avoidance scheme. To simplify handling, the input buffer
   // is enforced to be a multiple of 2*Hashinator::defaults::MAX_BLOCKSIZE in size.
   while ((majorOffset < flatExtent) && (offsetC<nBlocks)) {
      int offset = 1;
      // Load input into shared memory
      int ai = ti;
      int bi = ti + (n/2);
      int bankOffsetA = BANK_OFFSET(ai);
      int bankOffsetB = BANK_OFFSET(bi);
      reductionA[ai + bankOffsetA] = probeFlattened[majorOffset + ai];
      reductionA[bi + bankOffsetB] = probeFlattened[majorOffset + bi];
      reductionB[ai + bankOffsetA] = (probeFlattened[majorOffset + ai] != 0 ? 1 : 0);
      reductionB[bi + bankOffsetB] = (probeFlattened[majorOffset + bi] != 0 ? 1 : 0);
      reductionC[ai + bankOffsetA] = probeFlattened[flatExtent + majorOffset + ai];
      reductionC[bi + bankOffsetB] = probeFlattened[flatExtent + majorOffset + bi];

      // build sum in place up the tree
      for (int d = n>>1; d > 0; d >>= 1) {
         __syncthreads();
         if (ti < d) {
            int ai = offset*(2*ti+1)-1;
            int bi = offset*(2*ti+2)-1;
            ai += BANK_OFFSET(ai);
            bi += BANK_OFFSET(bi);
            reductionA[bi] += reductionA[ai];
            reductionB[bi] += reductionB[ai];
            reductionC[bi] += reductionC[ai];
         }
         offset *= 2;
      }
      // Clear the last element
      if (ti==0) {
         reductionA[n - 1 + BANK_OFFSET(n - 1)] = 0;
         reductionB[n - 1 + BANK_OFFSET(n - 1)] = 0;
         reductionC[n - 1 + BANK_OFFSET(n - 1)] = 0;
      }

      // traverse down tree & build scan
      for (int d = 1; d < n; d *= 2) {
         offset >>= 1;
         __syncthreads();
         if (ti < d) {
            int ai = offset*(2*ti+1)-1;
            int bi = offset*(2*ti+2)-1;
            ai += BANK_OFFSET(ai);
            bi += BANK_OFFSET(bi);

            vmesh::LocalID t = reductionA[ai];
            reductionA[ai] = reductionA[bi];
            reductionA[bi] += t;
            t = reductionB[ai];
            reductionB[ai] = reductionB[bi];
            reductionB[bi] += t;
            t = reductionC[ai];
            reductionC[ai] = reductionC[bi];
            reductionC[bi] += t;
         }
      }
      __syncthreads();

      // write results to device memory, increment majorOffset and offsetA/B/C.
      // Remember:
      // The flattened version must store:
      // 1) how many columns per potential column position (potColumn) (input for this kernel)
      // 2) how many blocks per potColumn (input for this kernel)
      // 3) cumulative offset into columns per potColumn (output for this kernel)
      // 4) cumulative offset into columnSets per potColumn (output for this kernel)
      // 5) cumulative offset into blocks per potColumn (output for this kernel)

      probeFlattened[2*flatExtent + majorOffset + ai] = reductionA[ai + bankOffsetA] + offsetA;
      probeFlattened[2*flatExtent + majorOffset + bi] = reductionA[bi + bankOffsetB] + offsetA;
      probeFlattened[3*flatExtent + majorOffset + ai] = reductionB[ai + bankOffsetA] + offsetB;
      probeFlattened[3*flatExtent + majorOffset + bi] = reductionB[bi + bankOffsetB] + offsetB;
      probeFlattened[4*flatExtent + majorOffset + ai] = reductionC[ai + bankOffsetA] + offsetC;
      probeFlattened[4*flatExtent + majorOffset + bi] = reductionC[bi + bankOffsetB] + offsetC;
      // Advance to reading next section of input buffer
      majorOffset += n;
      // Increment cumulative offset (exclusive sum result of last bin + contents of that one)
      __syncthreads();
      if (ti==0) {
         offsetA += reductionA[n-1] + probeFlattened[majorOffset-1];
         offsetB += reductionB[n-1] + (probeFlattened[majorOffset-1] != 0 ? 1 : 0);
         offsetC += reductionC[n-1] + probeFlattened[flatExtent + majorOffset-1];
      }
      __syncthreads();
   }
   __syncthreads();
   if (ti == 0) {
      // Store reduction results
      const vmesh::LocalID numCols = offsetA;
      const vmesh::LocalID numColSets = offsetB;
      //printf("found columns %u and columnsets %u, %u blocks vs %d\n",numCols,numColSets,offsetC,nBlocks);
      dev_numCols[parallelOffsetIndex] = numCols; // Total number of columns
      dev_numColSets[parallelOffsetIndex] = numColSets; // Total number of column sets
      // Resize device-side column offset container vectors. First verify capacity.
      // set dev_resizeSuccess to unity to indicate if re-capacitate on host is needed.
      if ( (columnData->dev_capacityCols() < numCols) ||
           (columnData->dev_capacityColSets() < numColSets) ) {
         dev_resizeSuccess[cellOffset] = 1;
         return;
      } else {
         dev_resizeSuccess[cellOffset] = 0;
      }
      columnData->device_setSizes(numCols,numColSets);
   }

   // Todo: unrolling  of reduction loops to get even more performance.
   // Memos:
   // Perform all-prefix-sum to gather offsets
   // Look at e.g.
   // https://developer.nvidia.com/gpugems/gpugems3/part-vi-gpu-computing/chapter-39-parallel-prefix-sum-scan-cuda
   // Example 39-2 onwards
   // See also  splitvector's stream compaction mechanism and
   // Credits to  https://www.eecs.umich.edu/courses/eecs570/hw/parprefix.pdf
   // Should also be made to work with arbitrary size buffers, not just powers-of-two
}

/*!
   \brief GPU kernel for building the columns and columnSets for each spatial cell, and storing
   offsets and lengths into the ColumnOffsets struct.

   Utilizing the previously computed offsets (with the probe cube), this parallel kernel
   builds the offsets required for columns. It read the contents of the probe cube and outputs the stored GIDs
   and LIDs into provided buffers.

   @param vmeshes pointer to buffer of pointers to vmeshes of all cells on this process
   @param dev_probeCubeData pointer to buffer of pointers, which point to allocated
    temporary arrays, which are recast and used as the probe cube.
   @param D0 Vmesh maximal extents in X-direction
   @param D1 Vmesh maximal extents in Y-direction
   @param D2 Vmesh maximal extents in Z-direction
   @param dimension direction of acceleration
   @param flatExtent the product of the non-accelerated dimensions of the probe cube, i.e. the
    extent of it when it's flattend, rounded up to provide memory alignment.
   @param invalidLID value to be used as the value for invalid vmesh::LocalID
   @param dev_columnOffsetData pointer to a buffer of ColumnOffsets structs (see arch/gpu_base.hpp) where
    per-cell information on all built columns and columnSets is stored
   @param dev_vbwcl_vec pointer to buffer of pointers to SplitVectors, re-cast to use as an (ordered) LID list
   @param cumulativeOffset the current cumulative offset at which to index the aforementioned buffers, using
    the grid index of this kernel on top of this provided offset.
*/
__global__ void build_column_offsets(
   vmesh::VelocityMesh** __restrict__ vmeshes,
   vmesh::LocalID* dev_probeCubeData, // recast to vmesh::LocalID *probeCube, *probeFlattened
   const vmesh::LocalID D0,
   const vmesh::LocalID D1,
   const vmesh::LocalID D2,
   const int dimension,
   const size_t flatExtent,
   const vmesh::LocalID invalidLID,
   ColumnOffsets* dev_columnOffsetData,
   split::SplitVector<vmesh::GlobalID> ** dev_vbwcl_vec, // use as LIDlist
   const uint cumulativeOffset,
   const size_t gpu_probeStride
   ) {
   // Probe cube contents have been ordered based on acceleration dimesion
   // so this kernel always reads in the same way.

   const int ti = threadIdx.x; // [0,Hashinator::defaults::MAX_BLOCKSIZE)
   const vmesh::LocalID ind = blockDim.x * blockIdx.x + ti;

   const uint parallelOffsetIndex = blockIdx.y;
   const uint cellOffset = parallelOffsetIndex + cumulativeOffset;

   // Caller function verified this cast is safe
   vmesh::LocalID* LIDlist = reinterpret_cast<vmesh::LocalID*>(dev_vbwcl_vec[cellOffset]->data());
   ColumnOffsets* columnData = dev_columnOffsetData + parallelOffsetIndex;

   vmesh::LocalID* probeFlattened = dev_probeCubeData + parallelOffsetIndex * gpu_probeStride;
   vmesh::LocalID* probeCube = probeFlattened + flatExtent*GPU_PROBEFLAT_N;

   // definition: potColumn is a potential column(set), i.e. a stack from the probe cube.
   // potColumn indexes/offsets into columnData and LIDlist
   const vmesh::LocalID N_cols = probeFlattened[ind];
   const vmesh::LocalID offset_cols = probeFlattened[2*flatExtent + ind];
   const vmesh::LocalID offset_colsets = probeFlattened[3*flatExtent + ind];
   const vmesh::LocalID offset_blocks = probeFlattened[4*flatExtent + ind];

   // Here we use ind to back-calculate the transverse "x" and "y" indices (i,j) of the column(set).
   // which is by agreement propagated in the "z"-direction. TODO: This could probably be done through
   // multiplication of indices and pre-computed multipliers as is done in many other kernels, but this
   // works well enough.
   vmesh::LocalID i,j;
   vmesh::LocalID Dacc, Dother;
   switch (dimension) {
      case 0:
         // propagate along x
         Dacc = D0;
         Dother = D1*D2;
         i = ind % D1; // Z (last dimension)
         j = ind / D1; // Y
         break;
      case 1:
         // propagate along y
         Dacc = D1;
         Dother = D0*D2;
         i = ind / D2; // X
         j = ind % D2; // Z (last dimension)
         break;
      case 2:
         // propagate along z
         Dacc = D2;
         Dother = D0*D1;
         i = ind / D1; // X
         j = ind % D1; // Y (last dimension)
         break;
      default:
         assert("ERROR! incorrect dimension!\n");
         return;
   }
   if (ind < Dother) {
      if (N_cols != 0) {
         // Update values in columnSets vector
         columnData->setColumnOffsets[offset_colsets] = offset_cols;
         columnData->setNumColumns[offset_colsets] = N_cols;
      }
      // Per-thread counters
      vmesh::LocalID foundBlocks = 0;
      vmesh::LocalID foundBlocksThisCol = 0;
      vmesh::LocalID foundCols = 0;
      bool inCol = false;

      // Loop through acceleration dimension of cube
      for (vmesh::LocalID k = 0; k < Dacc; k++) {
         // Early return when all columns have been completed
         if (foundCols >= N_cols) {
            return;
         }
         const vmesh::LocalID LID = probeCube[k*Dother + ind];
         if (LID == invalidLID) {
            // No block at this index.
            if (inCol) {
               // finish current column
               columnData->columnNumBlocks[offset_cols + foundCols] = foundBlocksThisCol;
               foundCols++;
               inCol = false;
            }
         } else {
            // Valid block found at this index!
            // Store LID into buffer
            LIDlist[offset_blocks + foundBlocks] = LID;
            if (!inCol) {
               // start new column
               inCol = true;
               foundBlocksThisCol = 0;
               columnData->columnBlockOffsets[offset_cols + foundCols] = offset_blocks + foundBlocks;
               columnData->i[offset_cols + foundCols] = i;
               columnData->j[offset_cols + foundCols] = j;
               columnData->kBegin[offset_cols + foundCols] = k;
            }
            foundBlocks++;
            foundBlocksThisCol++;
         }
      }
      // Finished loop. If we are "still in a colum", count that.
      if (inCol) {
         columnData->columnNumBlocks[offset_cols + foundCols] = foundBlocksThisCol;
      }
   }
}
/*!
   \brief GPU kernel for reading block data in from the spatial cell VelocityBlockContainers,
   transposing each block so the acceleration direction is the k-index, and storing each block
   into ordered buffers which will then be fed to the acceleration kernel.

   Reads block data in from spatial cells and places it in the ordered container.
   Each block is in column-based order. Works column-per-column and adds the necessary
   one empty block at each end

   Block contents are transposed so that values 0...WID3-1 (i.e. each block)
   are so that every WID3 consequtive elements form one "slice" perpendicular to the
   direction of acceleration. This is achieved by swapping the indexing order for two
   dimensions (if needed), in accordance with historic acceleration solver tradition.

   Example for WID=4:            Z
                                 ^   Y
         60   61   62   63       |  /
       56   57   58   59         | /
     52   53   54   55           |/
   48   49   50   51             *-----> X
         44   45   46   47
       40   41   42   43
     36   37   38   39
   32   33   34   35
         28   29   30   31
       24   25   26   27
     20   21   22   23
   16   17   18   19
         12   13   1    15
       8    9    10   11
     4    5    6    7
   0    1    2    3

   now the first 16 elements of the reodered buffer would correspond
   to source indices as follows:

   acceleration along X: (swap X and Z directions)
   0 16 32 48 4 20 36 52 8 24 40 12 28 44 60

   acceleration along Y: (swap Y and Z directions)
   0 1 2 3 16 17 18 19 32 33 34 35 48 49 50 51

   acceleration along Z: (no swaps)
   0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15

   This way the acceleration kernel can access cells from the input buffer
   such that there is a constant offset WID2 between cells along the same
   acceleration direction.

   @param blockContainers pointer to buffer of pointers to VelocityBlockContainers of all cells on this process
   @param dev_blockDataOrdered pointer to buffer of pointers, which point to allocated
    temporary arrays used for storing the phase-space data to be propagated.
   @param gpu_cell_indices_to_id buffer of 3 values used for converting block indices between directions
   @param dev_vbwcl_vec pointer to buffer of pointers to SplitVectors, re-cast to use as an (ordered) LID list
   @param dev_columnOffsetData pointer to a buffer of ColumnOffsets structs (see arch/gpu_base.hpp) where
    per-cell information on all built columns and columnSets is stored
   @param dev_nColumns pointer to buffer of data indicating how many columns each cell has
   @param cumulativeOffset the current cumulative offset at which to index the aforementioned buffers, using
    the grid index of this kernel on top of this provided offset.
 */
__global__ void __launch_bounds__(WID3) reorder_blocks_by_dimension_kernel(
   vmesh::VelocityBlockContainer** __restrict__ blockContainers,
   Realf** dev_blockDataOrdered,
   const uint* __restrict__ gpu_cell_indices_to_id,
   split::SplitVector<vmesh::GlobalID> ** dev_vbwcl_vec, // use as LIDlist
   ColumnOffsets* dev_columnOffsetData,
   vmesh::LocalID* dev_nColumns,
   const uint cumulativeOffset
   ) {
   // This is launched with block size (WID,WID,WID)
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   // Acceleration direction becomes "z"
   const uint sourcei = threadIdx.x*gpu_cell_indices_to_id[0]
      + threadIdx.y*gpu_cell_indices_to_id[1]
      + threadIdx.z*gpu_cell_indices_to_id[2];

   const uint iColumn = blockIdx.x;
   const uint parallelOffsetIndex = blockIdx.y;
   const uint cellOffset = parallelOffsetIndex + cumulativeOffset;
   // Early return if already dealt with all columns
   if (iColumn >= dev_nColumns[parallelOffsetIndex]) {
      return;
   }
   const vmesh::VelocityBlockContainer* __restrict__ blockContainer = blockContainers[cellOffset];
   ColumnOffsets* columnData = dev_columnOffsetData + parallelOffsetIndex;
   Realf *gpu_blockDataOrdered = dev_blockDataOrdered[parallelOffsetIndex];

   // Caller function verified this cast is safe
   vmesh::LocalID* LIDlist = reinterpret_cast<vmesh::LocalID*>(dev_vbwcl_vec[cellOffset]->data());

   // Each gpuBlock deals with one column.
   const uint inputOffset = columnData->columnBlockOffsets[iColumn];
   const uint outputOffset = (inputOffset + 2 * iColumn) * WID3;
   const uint columnLength = columnData->columnNumBlocks[iColumn];

   // Loop over column blocks
   for (uint b = 0; b < columnLength; b++) {
      #ifdef DEBUG_ACC
      assert((inputOffset + b) < blockContainer->size() && "reorder_blocks_by_dimension_kernel too large LID");
      #endif
      const vmesh::LocalID LID = LIDlist[inputOffset + b];
      const Realf* __restrict__ gpu_blockData = blockContainer->getData(LID);
      // Transpose block so that propagation direction becomes last dimension (z)
      gpu_blockDataOrdered[outputOffset + (1 + b) * WID3 + ti] = gpu_blockData[sourcei];
   }
   // Set first and last blocks to zero
   gpu_blockDataOrdered[outputOffset + ti] = 0.0;
   gpu_blockDataOrdered[outputOffset + (columnLength + 1) * WID3 + ti] = 0.0;
   // Note: this kernel does not memset gpu_blockData to zero, there is a separate kernel for that.
}

/*!
   \brief GPU Kernel for evaluating the constructed columns, which blocks need to be deleted from
   a velocity mesh, and which ones need to be added to it.

   For each column, evaluates which blocks are source blocks, which are target blocks,
   and which ones are both. Then goes through this list and pushes values into two sets
   and one vector.

   TODO: This kernel could probaly be streamlined, but it is fast to execute so hasn't been a priority.

   Note: max_v_length, v_min, dv could all be queried from any random vmesh instead of passing from host.

   @param dimension direction of acceleration
   @param vmeshes pointer to buffer of pointers to vmeshes of all cells on this process
   @param dev_columnOffsetData pointer to a buffer of ColumnOffsets structs (see arch/gpu_base.hpp) where
    per-cell information on all built columns and columnSets is stored
   @param lists_with_replace_new pointer to buffer of pointers to splitvectors, where newly added blocks are listed.
   @param allMaps pointer to buffer of pointers to maps, two per spatial cell, used for gathering deleted and required blocks.
   @param gpu_block_indices_to_probe 3-element array used for converting block GID to location
    within probe cube
   @param dev_intersections pointer to buffer of per-cell intersection values pre-calculated for this acceleration direction.
   @param bailout_velocity_space_wall_margin integer for how close to v-space walls we allow blocks to be created before bailout
   @param max_v_length maximum extent of v-space for this direction of acceleration
   @param v_min start of v-space for this direction of acceleration
   @param dv size of v for each velocity cell in this direction of acceleration
   @param dev_resizeSuccess pointer to buffer where the amount by which the list_with_replace_new capacity was exceeded
   @param dev_overflownElements pointer to buffer for bailout flag for touching v-space walls
   @param cumulativeOffset the current cumulative offset at which to index the aforementioned buffers, using
    the grid index of this kernel on top of this provided offset.
*/
   __global__ void __launch_bounds__(GPUTHREADS,4) evaluate_column_extents_kernel(
   const uint dimension,
   vmesh::VelocityMesh** __restrict__ vmeshes,
   ColumnOffsets* dev_columnOffsetData,
   split::SplitVector<vmesh::GlobalID>* *lists_with_replace_new,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* *allMaps,
   const uint* __restrict__ gpu_block_indices_to_id,
   const Realf* dev_intersections,
   const int bailout_velocity_space_wall_margin,
   const int max_v_length,
   const Realf v_min,
   const Realf dv,
   vmesh::LocalID *dev_resizeSuccess, // bailout flag: splitvector list_with_replace_new capacity error
   vmesh::LocalID *dev_overflownElements, // bailout flag: touching velspace wall
   const uint cumulativeOffset
   ) {
   const uint warpSize = blockDim.x;
   const uint setIndex = blockIdx.x;
   const uint ti = threadIdx.x;

   const uint parallelOffsetIndex = blockIdx.y;
   const uint cellOffset = parallelOffsetIndex + cumulativeOffset;

   const Realf intersection = dev_intersections[cellOffset*4+0];
   const Realf intersection_di = dev_intersections[cellOffset*4+1];
   const Realf intersection_dj = dev_intersections[cellOffset*4+2];
   const Realf intersection_dk = dev_intersections[cellOffset*4+3];
   ColumnOffsets* columnData = dev_columnOffsetData + parallelOffsetIndex;

   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *dev_map_require = allMaps[2*cellOffset];
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *dev_map_remove = allMaps[2*cellOffset+1];
   split::SplitVector<vmesh::GlobalID> *list_with_replace_new = lists_with_replace_new[cellOffset];
   const vmesh::VelocityMesh* __restrict__ vmesh = vmeshes[cellOffset];
   // Shared within all threads in one block (one columnSet)
   __shared__ int isTargetBlock[MAX_BLOCKS_PER_DIM];
   __shared__ int isSourceBlock[MAX_BLOCKS_PER_DIM];

   if (setIndex < columnData->setColumnOffsets.size()) {

      // Clear flags used for this columnSet
      for(uint tti = 0; tti < MAX_BLOCKS_PER_DIM; tti += warpSize ) {
         const uint index = tti + ti;
         if (index < MAX_BLOCKS_PER_DIM) {
            isTargetBlock[index] = 0;
            isSourceBlock[index] = 0;
         }
      }
      __syncthreads();

      /*need x,y coordinate of this column set */
      const vmesh::LocalID set_i = columnData->i[columnData->setColumnOffsets[setIndex]];
      const vmesh::LocalID set_j = columnData->j[columnData->setColumnOffsets[setIndex]];

      /* Compute the maximum starting point of the lagrangian (target) grid
         within the 4 corner cells in this block. Needed for computing
         maximum extent of target column.
      */

      Realf intersectionMins[4];
      intersectionMins[0] = intersection + (set_i * WID + 0) * intersection_di +
         (set_j * WID + 0) * intersection_dj;
      intersectionMins[1] = intersection + (set_i * WID + 0) * intersection_di +
         (set_j * WID + WID - 1) * intersection_dj;
      intersectionMins[2] = intersection + (set_i * WID + WID - 1) * intersection_di +
         (set_j * WID + 0) * intersection_dj;
      intersectionMins[3] = intersection + (set_i * WID + WID - 1) * intersection_di +
         (set_j * WID + WID - 1) * intersection_dj;

      Realf min_intersectionMin = std::min(std::min(intersectionMins[0],intersectionMins[1]),
                                           std::min(intersectionMins[2],intersectionMins[3]));
      Realf max_intersectionMin = std::max(std::max(intersectionMins[0],intersectionMins[1]),
                                           std::max(intersectionMins[2],intersectionMins[3]));

      // Now record which blocks are target blocks
      for (uint columnIndex = columnData->setColumnOffsets[setIndex];
           columnIndex < columnData->setColumnOffsets[setIndex] + columnData->setNumColumns[setIndex] ;
           ++columnIndex) {
         // Not parallelizing this at this level; not going to be many columns within a set
         // (and we want to manage each columnSet within one block)

         const vmesh::LocalID n_cblocks = columnData->columnNumBlocks[columnIndex];
         const vmesh::LocalID kBegin = columnData->kBegin[columnIndex];
         const vmesh::LocalID kEnd = kBegin + n_cblocks -1;

         /* firstBlockV is in z the minimum velocity value of the lower
          *  edge in source grid.
          * lastBlockV is in z the maximum velocity value of the upper
          *  edge in source grid. */
         const Realf firstBlockMinV = (WID * kBegin) * dv + v_min;
         const Realf lastBlockMaxV = (WID * (kEnd + 1)) * dv + v_min;

         /* gk is now the k value in terms of cells in target
            grid. This distance between max_intersectionMin (so lagrangian
            plan, well max value here) and V of source grid, divided by
            intersection_dk to find out how many grid cells that is*/
         const int firstBlock_gk = (int)((firstBlockMinV - max_intersectionMin)/intersection_dk);
         const int lastBlock_gk = (int)((lastBlockMaxV - min_intersectionMin)/intersection_dk);

         int firstBlockIndexK = firstBlock_gk/WID;
         int lastBlockIndexK = lastBlock_gk/WID;

         // now enforce mesh limits for target column blocks (and check if we are
         // too close to the velocity space boundaries)
         firstBlockIndexK = (firstBlockIndexK >= 0)            ? firstBlockIndexK : 0;
         firstBlockIndexK = (firstBlockIndexK < max_v_length ) ? firstBlockIndexK : max_v_length - 1;
         lastBlockIndexK  = (lastBlockIndexK  >= 0)            ? lastBlockIndexK  : 0;
         lastBlockIndexK  = (lastBlockIndexK  < max_v_length ) ? lastBlockIndexK  : max_v_length - 1;
         if(firstBlockIndexK < bailout_velocity_space_wall_margin
            || firstBlockIndexK >= max_v_length - bailout_velocity_space_wall_margin
            || lastBlockIndexK < bailout_velocity_space_wall_margin
            || lastBlockIndexK >= max_v_length - bailout_velocity_space_wall_margin
            ) {
            // Pass bailout (hitting the wall) flag back to host
            if (ti==0) {
               dev_overflownElements[cellOffset] = 1;
            }
         }

         //store source blocks
         for (uint blockK = kBegin; blockK <= kEnd; blockK +=warpSize){
            if ((blockK+ti) <= kEnd) {
               isSourceBlock[blockK+ti] = 1; // Does not need to be atomic, as long as it's no longer zero
            }
         }
         __syncthreads();

         //store target blocks
         for (uint blockK = (uint)firstBlockIndexK; blockK <= (uint)lastBlockIndexK; blockK+=warpSize){
            if ((blockK+ti) <= (uint)lastBlockIndexK) {
               isTargetBlock[blockK+ti] = 1; // Does not need to be atomic, as long as it's no longer zero
            }
         }
         __syncthreads();

         if (ti==0) {
            // Store for each column firstBlockIndexK, and lastBlockIndexK
            columnData->minBlockK[columnIndex] = firstBlockIndexK;
            columnData->maxBlockK[columnIndex] = lastBlockIndexK;
         }
      } // end loop over columns in set
      __syncthreads();

      for (uint blockT = 0; blockT < MAX_BLOCKS_PER_DIM; blockT +=warpSize) {
         const uint blockK = blockT + ti;
         // Not using warp accessors, as each thread has different block
         if (blockK < MAX_BLOCKS_PER_DIM) {
            if (isTargetBlock[blockK] != 0) {
               const int targetBlock =
                  set_i  * gpu_block_indices_to_id[0] +
                  set_j  * gpu_block_indices_to_id[1] +
                  blockK * gpu_block_indices_to_id[2];
               // Templated parameter: do not overwrite existing values
               dev_map_require->set_element<true>(targetBlock, vmesh->getLocalID(targetBlock));
            }
            if (isTargetBlock[blockK] !=0 && isSourceBlock[blockK] == 0 )  {
               const int targetBlock =
                  set_i  * gpu_block_indices_to_id[0] +
                  set_j  * gpu_block_indices_to_id[1] +
                  blockK * gpu_block_indices_to_id[2];
               if (!list_with_replace_new->device_push_back(targetBlock)) {
                  // out of capacity, bailout and gather how much capacity needs to grow
                  atomicAdd(&dev_resizeSuccess[cellOffset],1);
               }
            }
            if (isTargetBlock[blockK] == 0 && isSourceBlock[blockK] != 0 )  {
               const int targetBlock =
                  set_i  * gpu_block_indices_to_id[0] +
                  set_j  * gpu_block_indices_to_id[1] +
                  blockK * gpu_block_indices_to_id[2];
               // Templated parameter: do not overwrite existing values
               dev_map_remove->set_element<true>(targetBlock, vmesh->getLocalID(targetBlock));
            }
         } // block within MAX_BLOCKS_PER_DIM
      } // loop over all potential blocks
   } // if valid setIndex
}

// Use max 2048 per MP threads due to register usage limitations
#if THREADS_PER_MP < (REGISTERS_PER_MP/64 + 1)
  #define ACCELERATION_KERNEl_MIN_BLOCKS THREADS_PER_MP/(WID3)
#else
  #define ACCELERATION_KERNEl_MIN_BLOCKS (REGISTERS_PER_MP/64)/(WID3)
#endif
/*!
   \brief GPU kernel for main task of semi-Lagrangian acceleration. Reads data in from buffer,
   performs polynomial reconstruction and does piecewise integration of contribution,
   and stores results directly into VelocityBlockContainer.

   Note: v_min, dv could all be queried from any random vmesh instead of passing from host.

   @param vmeshes pointer to buffer of pointers to vmeshes of all cells on this process
   @param blockContainers pointer to buffer of pointers to VelocityBlockContainers of all cells on this process
   @param dev_blockDataOrdered pointer to buffer of pointers, which point to allocated
    temporary arrays used for storing the phase-space data to be propagated.
   @param gpu_cell_indices_to_id buffer of 3 values used for converting block indices between directions
   @param gpu_block_indices_to_probe 3-element array used for converting propageted block indexes bacl to GID
   @param dev_columnOffsetData pointer to a buffer of ColumnOffsets structs (see arch/gpu_base.hpp) where
    per-cell information on all built columns and columnSets is stored
   @param dev_intersections pointer to buffer of per-cell intersection values pre-calculated for this acceleration direction.
   @param v_min start of v-space for this direction of acceleration
   @param i_dv unity per size of v for each velocity cell in this direction of acceleration
   @param dv size of v for each velocity cell in this direction of acceleration
   @param dev_minvalues pointer to buffer of stored sparsity minValues for each cell, used by slope limiters
   @param invalidLID value to be used as the value for invalid vmesh::LocalID
   @param cumulativeOffset the current cumulative offset at which to index the aforementioned buffers, using
    the grid index of this kernel on top of this provided offset.
 */
__global__ void __launch_bounds__(WID3, ACCELERATION_KERNEl_MIN_BLOCKS) acceleration_kernel(
   vmesh::VelocityMesh** __restrict__ vmeshes, // indexing: cellOffset
   vmesh::VelocityBlockContainer **blockContainers, // indexing: cellOffset
   Realf** __restrict__ dev_blockDataOrdered, //indexing: blockIdx.y
   const uint* __restrict__ gpu_cell_indices_to_id,
   const uint* __restrict__ gpu_block_indices_to_id,
   ColumnOffsets* __restrict__ dev_columnOffsetData, //indexing: blockIdx.y
   const Realf *dev_intersections, // indexing: cellOffset
   const Realf v_min,
   const Realf i_dv,
   const Realf dv,
   const Real *dev_minValues, // indexing: cellOffset
   const size_t invalidLID,
   const uint cumulativeOffset
) {
   const uint parallelOffsetIndex = blockIdx.y; // which vlasov buffer allocation to access
   const uint cellOffset = parallelOffsetIndex + cumulativeOffset;

   // This is launched with block size (WID,WID,WID)
   // Indexes into transposed data blocks
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z; // Acceleration direction
   const int ij = threadIdx.x + threadIdx.y * blockDim.x; // transverse index
   const int ti = ij + k*blockDim.x*blockDim.y;

   const Realf* __restrict__ gpu_blockDataOrdered = dev_blockDataOrdered[parallelOffsetIndex];
   const ColumnOffsets* __restrict__ columnData = dev_columnOffsetData + parallelOffsetIndex;

   const Realf intersection = dev_intersections[cellOffset*4+0];
   const Realf intersection_di = dev_intersections[cellOffset*4+1];
   const Realf intersection_dj = dev_intersections[cellOffset*4+2];
   const Realf intersection_dk = dev_intersections[cellOffset*4+3];

   const vmesh::VelocityMesh* __restrict__ vmesh = vmeshes[cellOffset];
   vmesh::VelocityBlockContainer *blockContainer = blockContainers[cellOffset];
   Realf *gpu_blockData = blockContainer->getData();

   // Load minvalues to shared memory
   __shared__ Realf minValue;
   __shared__ uint setColumnOffset;
   __shared__ uint numColumns;

   // shared memory buffer for reducing looping count per block
   __shared__ int loopN[WID3/GPUTHREADS];

   {
      const uint setIndex = blockIdx.x;

      if (setIndex >= columnData->dev_sizeColSets()) {
         return;
      }

      if (ti == 0) {
         minValue = (Realf)dev_minValues[cellOffset];
         setColumnOffset = columnData->setColumnOffsets[setIndex];
         numColumns = columnData->setNumColumns[setIndex];
      }
   }

   __syncthreads();

   // Kernel must loop over all columns in set to ensure correct writes
   for (uint columnIndex = 0;
        columnIndex < numColumns;
        ++columnIndex) {
      
      const uint column = setColumnOffset + columnIndex;

      const Realf v_r0 = ( (Realf)(WID * columnData->kBegin[column]) * dv + v_min);
      const int nBlocks = columnData->columnNumBlocks[column];
      const int col_i = columnData->i[column];
      const int col_j = columnData->j[column];
      // Target block-k values for column
      const int col_mink = columnData->minBlockK[column];
      const int col_maxk = columnData->maxBlockK[column];
      const size_t stencilDataOffset = (columnData->columnBlockOffsets[column] + 2*column) * WID3;

      // Column index contribution for adjusted velocity block container at correct target GID/LID
      const vmesh::GlobalID columnGID = col_i  * gpu_block_indices_to_id[0] + col_j  * gpu_block_indices_to_id[1];

      // Intersection for this cell
      const Realf intersection_min =
         intersection
         + intersection_di * (Realf)(col_i * WID + i)
         + intersection_dj * (Realf)(col_j * WID + j);
      // Pre-computed constant target offset contribution
      const int target_cell_index_common = i * gpu_cell_indices_to_id[0]
         + j * gpu_cell_indices_to_id[1];

      // Loop over blocks in column
      for (int b = 0; b < nBlocks; b++) {
         const int blockOffset = WID * b; // in units k

         int minGk;

         {
            // Min/max Velocity coordinates in acceleration direction for this block
            const Realf min_lagrangian_v_l = v_r0 + blockOffset * dv;
            const Realf max_lagrangian_v_r = v_r0 + (blockOffset + WID) * dv;

            // Sub-column (single i and j) target k-index extent
            const int subcolumnMinGk = int(trunc((min_lagrangian_v_l - intersection_min)/intersection_dk));
            const int subcolumnMaxGk = int(trunc((max_lagrangian_v_r - intersection_min)/intersection_dk));

            // Truncate to possible output block values
            // min-value decreased by (WID-1) so even last slice in sub-column gets to calculate from first gk-index
            minGk = std::max(subcolumnMinGk, col_mink * WID) - (WID-1);
            const int maxGk = std::min(subcolumnMaxGk, (col_maxk + 1) * WID - 1);

            // Reduce Gk loop count
            int indexInsideWarp = ti % GPUTHREADS;
            int warpIndex = ti / GPUTHREADS;

            int val = maxGk - minGk + 1;

            for (int offset = GPUTHREADS/2; offset > 0; offset /= 2) {
               val = max(val, gpuKernelShflDown(val, offset));
            }

            if (indexInsideWarp == 0) {
               loopN[warpIndex] = val;
            }

            __syncthreads();

            if (warpIndex == 0) {
               val = (indexInsideWarp < WID3/GPUTHREADS) ? loopN[indexInsideWarp] : INT_MIN;
               for (int offset = GPUTHREADS/2; offset > 0; offset /= 2) {
                  val = max(val, gpuKernelShflDown(val, offset));
               }

               if (indexInsideWarp == 0) {
                  loopN[0] = val; // Store final result
               }
            }
         }

         __syncthreads();

         // Velocity coordinate in acceleration direction for this cell
         const Realf v_l = v_r0 + (blockOffset + k) * dv;
         const Realf v_r = v_r0 + (blockOffset + k + 1) * dv;
         // Target k-indexing for this cell
         const int lagrangian_gk_l = std::trunc((v_l-intersection_min)/intersection_dk);
         const int lagrangian_gk_r = std::trunc((v_r-intersection_min)/intersection_dk);

         // Compute reconstruction coefficients using WID2 as stride per slice
         // read from the offset for this column + the count of source blocks + 1 for an empty source block to begin with
         const size_t valuesOffset = stencilDataOffset + (b + 1) * WID3;
#ifdef ACC_SEMILAG_PLM
         Realf a[2];
         compute_plm_coeff(gpu_blockDataOrdered + valuesOffset, k, a, minValue, ij, WID2);
#endif
#ifdef ACC_SEMILAG_PPM
         Realf a[3];
         compute_ppm_coeff(gpu_blockDataOrdered + valuesOffset, h4, k, a, minValue, ij, WID2);
#endif
#ifdef ACC_SEMILAG_PQM
         Realf a[5];
         compute_pqm_coeff(gpu_blockDataOrdered + valuesOffset, h8, k, a, minValue, ij, WID2);
#endif

         // set the initial value for the integrand at the boundary at v = 0
         // (in reduced cell units), this will be shifted to target_density_1, see below.
         Realf target_density_r = (Realf)(0.0);

         // Perform the polynomial reconstruction for all cells the mapping streches into
         for(int loopgk = 0; loopgk < loopN[0]; loopgk++) {
            // Each cell within the subcolumn needs to consider a different gk value so writes don't overlap
            const int gk = minGk + loopgk + k;
            // // Does this cell need to consider this target gk?
            if (gk >= lagrangian_gk_l && gk <= lagrangian_gk_r) {
               const int blockK = gk/WID;
               const int gk_mod_WID = (gk - blockK * WID);
               // the velocities between which we will integrate, in order to put mass
               // into the target cell. If both v_r and v_l are in same cell
               // then v_1,v_2 should be between v_l and v_r.
               // v_1 and v_2 normalized to be between 0 and 1 in the cell.
               const Realf v_norm_r = (  std::min(  std::max( (gk + 1) * intersection_dk + intersection_min, v_l), v_r) - v_l) * i_dv;

               // shift, old right integrand is new left integrand
               const Realf target_density_l = target_density_r;

               // compute right integrand using FMA
               #ifdef ACC_SEMILAG_PLM
               target_density_r = a[1];
               target_density_r = a[0] + v_norm_r * target_density_r;
               target_density_r = v_norm_r * target_density_r;
               #endif
               #ifdef ACC_SEMILAG_PPM
               target_density_r = a[2];
               target_density_r = a[1] + v_norm_r * target_density_r;
               target_density_r = a[0] + v_norm_r * target_density_r;
               target_density_r = v_norm_r * target_density_r;
               #endif
               #ifdef ACC_SEMILAG_PQM
               target_density_r = a[4];
               target_density_r = a[3] + v_norm_r * target_density_r;
               target_density_r = a[2] + v_norm_r * target_density_r;
               target_density_r = a[1] + v_norm_r * target_density_r;
               target_density_r = a[0] + v_norm_r * target_density_r;
               target_density_r = v_norm_r * target_density_r;

               //target_density_r = v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * ( a[2] + v_norm_r * ( a[3] + v_norm_r * a[4] ) ) ) );
               #endif

               // integral area between the two integrands
               Realf tval = target_density_r - target_density_l;

               // Store directly into adjusted velocity block container at correct target GID/LID
               const vmesh::GlobalID targetGID = columnGID + blockK * gpu_block_indices_to_id[2];
               const vmesh::LocalID targetLID = vmesh->getLocalID(targetGID);
               // The target velocity cell within the target bloxk
               const int tcell = target_cell_index_common
                               + gk_mod_WID * gpu_cell_indices_to_id[2];
               // Write values into block data
               if (isfinite(tval) && (tval>(Realf)(0.0)) && (targetLID != invalidLID) ) {
                  // gpu_blockData[targetLID * WID3 + tcell] += tval;

                  // We use atomicAdd to avoid the need for sync
                  // It shouldn't be any slower if there is no competition
                  atomicAdd(&gpu_blockData[targetLID*WID3+tcell],tval);
               }
            } // end check if gk valid for this thread
         } // for loop over target k-indices
      } // for-loop over source blocks
   } // End this column
} // end semilag acc kernel


/*!
  \brief This function performs the semi-Lagrangian acceleration for a provided list of
  spatial cells, for one popID, for one dimension. See gpu_acc_semilag.cpp for
  Information on the calling structure.

  First, kernels are called to construct the existing columns and columnsets for the
  sparse velocity space data contained in this cell. A parallel launch and probe cubes
  are used for this, along with some unified memory column containers. Data from the
  current spatial cell is then read into an intermediate buffer in column-aligned order.

  Then, utilizing pre-calculated SLICE-3D intersections, a kernel is launched to evaluate
  the extents of velocity space after each column(set) has been accelerated as requested.
  Block adjustment functions are called to batch-update the velocity space of the cell
  in question to match the new requirements (adding and removing cells as needed). Then,
  the velocity block container is cleared.

  Finally the acceleration kernel itself is called. It reads the aligned velocity space data
  from the intermediate buffer, advects the columns according to the SLICE-3D intersections,
  and stores the resultant phase space density back into the velocity block container.

 * @param mpiGrid DCCRG container of spatial cells
 * @param launchCells vector of cells for which to perform acceleration in this "chunk"
 * @param popID ID of the accelerated particle species.
 * @param dimension Velocity dimension for acceleration (VX, VY, or VZ)
 * @param Dacc Maximal velocity block extent in accelerated dimension
 * @param Dother Product of maximal velocity block extents in non-accelerated dimensions
 * @param cumulativeOffset running counter for offset of indexing of launchCells into device buffers
*/
__host__ bool gpu_acc_map_1d(
   dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   vector<CellID> &launchCells,
   const uint popID,
   const uint dimension,
   const int Dacc, // velocity block max dimension, direction of acceleration
   const int Dother, // Product of other two dimensions (max blocks)
   const size_t cumulativeOffset
   ) {

   phiprof::Timer prepTimer {"preparation"};

   // Empty meshes are already excluded in gpu_acc_semilag.cpp
   // Sample vmesh from first cell
   vmesh::VelocityMesh* sampleVmesh    = mpiGrid[launchCells[0]]->get_velocity_mesh(popID);
   // These are constant for all cells included in this launch
   const vmesh::LocalID D0 = sampleVmesh->getGridLength()[0];
   const vmesh::LocalID D1 = sampleVmesh->getGridLength()[1];
   const vmesh::LocalID D2 = sampleVmesh->getGridLength()[2];
   const Realf dv    = sampleVmesh->getCellSize()[dimension];
   const Realf v_min = sampleVmesh->getMeshMinLimits()[dimension];
   const int max_v_length  = (int)sampleVmesh->getGridLength()[dimension];
   const Realf i_dv = 1.0/dv;
   const vmesh::LocalID invalidLocalID = sampleVmesh->invalidLocalID();

   const gpuStream_t baseStream = gpu_getStream();

   /** New merged kernel approach without sorting for columns

       First, we generate a "probe cube". It started off as an actual
       cube, but the transverse dimensions are considered as one.
       One dimension is that of the current acceleration (Dacc), and the other
       dimension is the product of the other two maximal velocity block
       domain extents (Dother). The "flattened" version is one where data is gathered
       over the acceleration direction into a single reduced value.

       The flattened version of the probe cube must store:
       1) how many columns per potential column position (potColumn)
       2) how many blocks per potColumn
       3) cumulative offset into columns per potColumn
       4) cumulative offset into columnSets per potColumn
       5) cumulative offset into blocks per potColumn

       For reductions, each slice of the flattened array should have a size a multiple of 2*MAX_BLOCKSIZE:
   */
   //const size_t flatExtent = 2*Hashinator::defaults::MAX_BLOCKSIZE * (1 + ((Dother - 1) / (2*Hashinator::defaults::MAX_BLOCKSIZE)));
   // This is now pre-computed in gpu_base.cpp
   const size_t flatExtent = gpu_probeFlattenedSize;

   /**
      For the gathered LIDlist, we re-use the allocation of spatial_cell->dev_velocity_block_with_content_list,
      It which contains variables of type vmesh::GlobalID (which should be the same as vmesh::LocalID, uint_32t).
      To ensure this static_cast is safe, we verify the sizes.
   */
   if constexpr (sizeof(vmesh::LocalID) != sizeof(vmesh::GlobalID)) {
      string message = " ERROR! vmesh::LocalID and vmesh::GlobalID are of different sizes, and thus";
      message += " the acceleration solver cannot safely use the spatial_cell->dev_list_delete";
      message += " Hashinator::splitVector object for storing a list of LIDs.";
      bailout(true, message, __FILE__, __LINE__);
   }

   const uint nLaunchCells = launchCells.size();
   size_t largestSizePower = 0;
   size_t largestNBefore = 0;
   vmesh::LocalID largest_totalColumns = 0;
   vmesh::LocalID largest_totalColumnSets = 0;
   vmesh::LocalID largest_nAfter = 0;

   for (size_t cellIndex = 0; cellIndex < nLaunchCells; cellIndex++) {
      const CellID cid = launchCells[cellIndex];
      const SpatialCell* SC = mpiGrid[cid];
      largestSizePower = std::max(largestSizePower, (size_t)SC->vbwcl_sizePower);
      largestSizePower = std::max(largestSizePower, (size_t)SC->vbwncl_sizePower);
      largestNBefore = std::max(largestNBefore, (size_t)SC->get_number_of_velocity_blocks(popID));
   }
   prepTimer.stop();

   phiprof::Timer clearTimer {"clear and prepare probe buffers"};
   // Clear hash maps used to evaluate block updates
   clear_maps_caller(nLaunchCells,largestSizePower,0,cumulativeOffset);

   gpu_calculateProbeAllocation(nLaunchCells);
   gpuMemoryManager.startSession(0,0);
   
   SESSION_ALLOCATE(gpuMemoryManager, vmesh::LocalID, dev_probeCubeData, gpu_getAllocationCount()*gpu_probeStride*sizeof(vmesh::LocalID));

   // probe cube and flattened version now re-use gpu_probeCubeData[cpuThreadID].
   // Due to alignment, Flattened version is at start of buffer, followed by the cube.
   // Required allocation sizes are calculated in gpu_base.cpp
   // Solve launch grid
   const size_t probeCombinedSize = gpu_probeFullSize + gpu_probeFlattenedSize * GPU_PROBEFLAT_N;
   const size_t n_prefill = 1 + ((probeCombinedSize - 1) / Hashinator::defaults::MAX_BLOCKSIZE);
   // This kernel fills the probe cube with invalid values and the flattened one with zeroes.
   const dim3 grid_prefill_probe(n_prefill,nLaunchCells,1);
   prefill_probe_kernel<<<grid_prefill_probe,Hashinator::defaults::MAX_BLOCKSIZE,0,baseStream>>>(
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
      GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_probeCubeData), // recast to vmesh::LocalID *probeCube
      flatExtent,
      Dacc,
      Dother,
      invalidLocalID,
      // Pass vectors for clearing
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new),
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_delete),
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_to_replace),
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_with_replace_old),
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec), // dev_velocity_block_with_content_list, // Resize to use as LIDlist
      cumulativeOffset,
      gpu_probeStride
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   clearTimer.stop();

   phiprof::Timer fillTimer {"fill probe cube"};
   // Read in GID list from vmesh, store LID values into probe cube in correct order
   // Launch params, fast ceil for positive ints
   const size_t n_fill_ord = 1 + ((largestNBefore - 1) / Hashinator::defaults::MAX_BLOCKSIZE);
   const dim3 grid_fill_ord(n_fill_ord,nLaunchCells,1);
   fill_probe_ordered<<<grid_fill_ord,Hashinator::defaults::MAX_BLOCKSIZE,0,baseStream>>>(
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
      GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_probeCubeData), // recast to vmesh::LocalID *probeCube
      flatExtent,
      GET_POINTER(gpuMemoryManager, uint, gpu_block_indices_to_probe),
      cumulativeOffset,
      gpu_probeStride
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   fillTimer.stop();

   // Now we perform reductions / flattenings / scans of the probe cube.
   // The kernel loops over the acceleration direction (Dacc).
   phiprof::Timer flattenTimer {"flatten probe cube"};
   const size_t n_grid_cube = 1 + ((Dother - 1) / Hashinator::defaults::MAX_BLOCKSIZE);
   const dim3 grid_cube(n_grid_cube,nLaunchCells,1);
   flatten_probe_cube<<<grid_cube,Hashinator::defaults::MAX_BLOCKSIZE,0,baseStream>>>(
      GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_probeCubeData), // recast to vmesh::LocalID *probeCube, *probeFlattened
      Dacc,
      Dother,
      flatExtent,
      invalidLocalID,
      gpu_probeStride
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   flattenTimer.stop();

   /*
     This kernel performs an exclusive prefix scan to get offsets for storing
     data from potential columns into the columnData container. Also gives us the total
     counts of columns, columnsets, and blocks, and uses the first two to resize
     our splitvector containers inside columnData.

     A proper prefix scan needs to be a two-phase process, thus two kernels,
     but here we do an iterative loop processing MAX_BLOCKSIZE elements at once.
     Not as efficient but simpler, and will be parallelized over spatial cells.
   */
   SESSION_HOST_ALLOCATE(gpuMemoryManager, vmesh::LocalID, host_nColumns, nLaunchCells*sizeof(vmesh::LocalID));
   SESSION_HOST_ALLOCATE(gpuMemoryManager, vmesh::LocalID, host_nColumnSets, nLaunchCells*sizeof(vmesh::LocalID));
   SESSION_ALLOCATE(gpuMemoryManager, vmesh::LocalID, dev_nColumns, nLaunchCells*sizeof(vmesh::LocalID));
   SESSION_ALLOCATE(gpuMemoryManager, vmesh::LocalID, dev_nColumnSets, nLaunchCells*sizeof(vmesh::LocalID));

   phiprof::Timer scanTimer {"scan probe cube"};
   const dim3 grid_scan(1,nLaunchCells,1);
   scan_probe<<<grid_scan,Hashinator::defaults::MAX_BLOCKSIZE,0,baseStream>>>(
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
      GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_probeCubeData), // recast to vmesh::LocalID *probeFlattened
      Dacc,
      Dother,
      flatExtent,
      GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nColumns),
      GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nColumnSets),
      GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess),
      GET_POINTER(gpuMemoryManager, ColumnOffsets, dev_columnOffsetData),
      cumulativeOffset,
      gpu_probeStride
      );
   CHK_ERR( gpuPeekAtLastError() );

   // Copy back to host sizes of found columns etc
   CHK_ERR( gpuMemcpyAsync(GET_SESSION_HOST_POINTER(gpuMemoryManager, vmesh::LocalID, host_nColumns), GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nColumns), nLaunchCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_SESSION_HOST_POINTER(gpuMemoryManager, vmesh::LocalID, host_nColumnSets), GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nColumnSets), nLaunchCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_resizeSuccess)+cumulativeOffset, GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess)+cumulativeOffset, nLaunchCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   scanTimer.stop();

   phiprof::Timer allocTimer {"ensure allocations"};
   // Ensure allocations (faster without threading)
   for (size_t cellIndex = 0; cellIndex < nLaunchCells; cellIndex++) {
      uint cellOffset = cellIndex + cumulativeOffset;
      // Read count of columns and columnsets, calculate required size of buffers
      vmesh::LocalID host_totalColumns = (GET_SESSION_HOST_POINTER(gpuMemoryManager, vmesh::LocalID, host_nColumns))[cellIndex];
      vmesh::LocalID host_totalColumnSets = (GET_SESSION_HOST_POINTER(gpuMemoryManager, vmesh::LocalID, host_nColumnSets))[cellIndex];
      vmesh::LocalID host_recapacitateVectors = (GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_resizeSuccess))[cellOffset]; // resize of columnData vectors
      largest_totalColumns = std::max(largest_totalColumns,host_totalColumns);
      largest_totalColumnSets = std::max(largest_totalColumnSets,host_totalColumnSets);
      if (host_recapacitateVectors) {
         // Can't call CPU reallocation directly as then copies go out of sync.
         // This function updates both CPU and GPU copies correctly.
         gpu_acc_allocate_perthread(cellIndex, host_totalColumns, host_totalColumnSets);
      }
   } // end parallel region
   allocTimer.stop();

   // Now we have gathered all the required offsets into probeFlattened, and can
   // now launch a kernel which constructs the columns offsets in parallel.
   phiprof::Timer columnsTimer {"build columns"};
   build_column_offsets<<<grid_cube,Hashinator::defaults::MAX_BLOCKSIZE,0,baseStream>>>(
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
      GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_probeCubeData), // recast to vmesh::LocalID *probeCube, *probeFlattened
      D0,D1,D2,
      dimension,
      flatExtent,
      invalidLocalID,
      GET_POINTER(gpuMemoryManager, ColumnOffsets, dev_columnOffsetData),
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec), //dev_velocity_block_with_content_list, // use as LIDlist
      cumulativeOffset,
      gpu_probeStride
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   columnsTimer.stop();

   phiprof::Timer allocTimer2 {"ensure vlasov allocations"};
   // Ensure allocations
   for (size_t cellIndex = 0; cellIndex < nLaunchCells; cellIndex++) {
      uint cellOffset = cellIndex + cumulativeOffset;
      // Read count of columns and columnsets, calculate required size of buffers
      vmesh::LocalID host_totalColumns = (GET_SESSION_HOST_POINTER(gpuMemoryManager, vmesh::LocalID, host_nColumns))[cellIndex];

      const CellID cid = launchCells[cellIndex];
      SpatialCell *SC = mpiGrid[cid];
      const vmesh::VelocityMesh *thisVmesh = SC->get_velocity_mesh(popID);
      const vmesh::LocalID nBlocks = thisVmesh->size();
      gpu_vlasov_allocate_perthread(cellIndex, 2*host_totalColumns+nBlocks);
   } // end parallel region

   CHK_ERR( gpuMemcpy(GET_POINTER(gpuMemoryManager, Realf*, dev_blockDataOrdered), GET_POINTER(gpuMemoryManager, Realf*, host_blockDataOrdered), gpu_getAllocationCount()*sizeof(Realf*), gpuMemcpyHostToDevice) );
   
   allocTimer2.stop();


   // Launch kernels for transposing and ordering velocity space data into columns
   phiprof::Timer reorderTimer {"reorder blocks"};
   const dim3 grid_reorder(largest_totalColumns,nLaunchCells,1);
   const dim3 block_reorder(WID,WID,WID);
   reorder_blocks_by_dimension_kernel<<<grid_reorder, block_reorder, 0, baseStream>>> (
      GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs),
      GET_POINTER(gpuMemoryManager, Realf*, dev_blockDataOrdered),
      GET_POINTER(gpuMemoryManager, uint, gpu_cell_indices_to_id),
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec), //dev_velocity_block_with_content_list, // use as LIDlist
      GET_POINTER(gpuMemoryManager, ColumnOffsets, dev_columnOffsetData),
      GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nColumns),
      cumulativeOffset
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   reorderTimer.stop();

   gpuMemoryManager.endSession();

   phiprof::Timer extentsTimer {"column extents"};
   // Reset counters used for verifying sufficient vector capacities and not overflowing v-space
   CHK_ERR( gpuMemset(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess)+cumulativeOffset, 0, nLaunchCells*sizeof(vmesh::LocalID)) );
   CHK_ERR( gpuMemset(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements)+cumulativeOffset, 0, nLaunchCells*sizeof(vmesh::LocalID)) );

   // Calculate target column extents
   const dim3 grid_column_extents(largest_totalColumnSets,nLaunchCells,1);
   evaluate_column_extents_kernel<<<grid_column_extents, GPUTHREADS, 0, baseStream>>> (
      dimension,
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
      GET_POINTER(gpuMemoryManager, ColumnOffsets, dev_columnOffsetData),
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new),
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps),
      GET_POINTER(gpuMemoryManager, uint, gpu_block_indices_to_id),
      GET_POINTER(gpuMemoryManager, Realf, dev_intersections),
      Parameters::bailout_velocity_space_wall_margin,
      max_v_length,
      v_min,
      dv,
      GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess), // bailout flag: splitvector list_with_replace_new capacity error
      GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements), // bailout flag: touching velspace wall
      cumulativeOffset
      );
   CHK_ERR( gpuPeekAtLastError() );
   // Check whether we exceeded the column data splitVectors on the way or if we need to bailout due to hitting v-space edge
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_resizeSuccess)+cumulativeOffset, GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess)+cumulativeOffset, nLaunchCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_overflownElements)+cumulativeOffset, GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements)+cumulativeOffset, nLaunchCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   extentsTimer.stop();

   phiprof::Timer extents2Timer {"column extents 2"};
   bool needSecondLaunchColumnExtents = false;
   // Faster without threading
   for (size_t cellIndex = 0; cellIndex < nLaunchCells; cellIndex++) {
      SpatialCell* SC = mpiGrid[launchCells[cellIndex]];
      const uint cellOffset = cellIndex + cumulativeOffset;
      if ((GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_resizeSuccess))[cellOffset] != 0) {
         needSecondLaunchColumnExtents = true;
         // counter indicates how many vector additions failed due to out-of-capacity.
         // Recapacitate with added safety factor and gather extents again.
         size_t newCapacity = (size_t)((SC->getReservation(popID)+(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_resizeSuccess))[cellOffset])*BLOCK_ALLOCATION_FACTOR);
         SC->setReservation(popID, newCapacity);
         SC->applyReservation(popID);
         // Clear the vector which receives push_backs. The maps do not need to be cleared.
         SC->list_with_replace_new->clear();
      }
   } // end parallel for

   if (needSecondLaunchColumnExtents) {
      // Reset counters, upload new pointers to splitvectors
      CHK_ERR( gpuMemsetAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess)+cumulativeOffset, 0, nLaunchCells*sizeof(vmesh::LocalID), baseStream) );
      CHK_ERR( gpuMemsetAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements)+cumulativeOffset, 0, nLaunchCells*sizeof(vmesh::LocalID), baseStream) );
      // Think this might not be actually needed, but let's play safe
      CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new)+cumulativeOffset, GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_lists_with_replace_new)+cumulativeOffset, nLaunchCells*sizeof(split::SplitVector<vmesh::GlobalID>*), gpuMemcpyHostToDevice, baseStream) );
      // Launch kernel a second time (now capacity should be sufficient)
      evaluate_column_extents_kernel<<<grid_column_extents, GPUTHREADS, 0, baseStream>>> (
         dimension,
         GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
         GET_POINTER(gpuMemoryManager, ColumnOffsets, dev_columnOffsetData),
         GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new),
         GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps),
         GET_POINTER(gpuMemoryManager, uint, gpu_block_indices_to_id),
         GET_POINTER(gpuMemoryManager, Realf, dev_intersections),
         Parameters::bailout_velocity_space_wall_margin,
         max_v_length,
         v_min,
         dv,
         GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess), // bailout flag: splitvector list_with_replace_new capacity error
         GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements), // bailout flag: touching velspace wall
         cumulativeOffset
         );
      CHK_ERR( gpuPeekAtLastError() );
      // Check whether we exceeded the column data splitVectors on the way and ensure capacity was now sufficient.
      CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_resizeSuccess)+cumulativeOffset, GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess)+cumulativeOffset, nLaunchCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
      CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_overflownElements)+cumulativeOffset, GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements)+cumulativeOffset, nLaunchCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
      CHK_ERR( gpuStreamSynchronize(baseStream) );
   }
   extents2Timer.stop();

   // Bailout checks (faster without threading)
   for (size_t cellIndex = 0; cellIndex < nLaunchCells; cellIndex++) {
      SpatialCell* SC = mpiGrid[launchCells[cellIndex]];
      const uint cellOffset = cellIndex + cumulativeOffset;
      // Check if we need to bailout due to hitting v-space edge
      if ((GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_overflownElements))[cellOffset] != 0) { //host_wallspace_margin_bailout_flag
         string message = "Some target blocks in acceleration are going to be less than ";
         message += std::to_string(Parameters::bailout_velocity_space_wall_margin);
         message += " blocks away from the current velocity space walls for population ";
         message += getObjectWrapper().particleSpecies[popID].name;
         message += " at CellID ";
         message += std::to_string((uint)SC->parameters[CellParams::CELLID]);
         message += ". Consider expanding velocity space for that population.";
         bailout(true, message, __FILE__, __LINE__);
      }
      // Also bail out if recapacitation was insufficient.
      if ((GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_resizeSuccess))[cellOffset] != 0) {
         string message = "Recapacitation of added velocity blocks vector for population ";
         message += " blocks away from the current velocity space walls for population ";
         message += getObjectWrapper().particleSpecies[popID].name;
         message += " at CellID ";
         message += std::to_string((uint)SC->parameters[CellParams::CELLID]);
         message += " failed. This should not happen.";
         bailout(true, message, __FILE__, __LINE__);
      }
   } // end parallel for

   /** Use block adjustment callers / lambda rules for extracting required map contents,
       building up vectors to use for parallel adjustment.
   */

   phiprof::Timer extractTimer {"extract block adjust vectors"};
   // TODO: Launch these three extracts in parallel from different streams?
   // Finds Blocks (GID,LID) to be rescued from end of v-space
   extract_to_delete_or_move_caller(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+2*cumulativeOffset, //dev_has_content_maps, // input maps
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_with_replace_old)+cumulativeOffset, // output vecs
      NULL, // pass null to not store vector lengths
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes)+cumulativeOffset, // rule_meshes
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+2*cumulativeOffset+1, //dev_has_no_content_maps// rule_maps
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new)+cumulativeOffset, // rule_vectors
      nLaunchCells,
      baseStream
      );
   // Find Blocks (GID,LID) to be outright deleted
   extract_to_delete_or_move_caller(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+2*cumulativeOffset+1,//dev_has_no_content_maps, // input maps
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_delete)+cumulativeOffset, // output vecs
      NULL, // pass null to not store vector lengths
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes)+cumulativeOffset, // rule_meshes
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+2*cumulativeOffset+1, //dev_has_no_content_maps, // rule_maps
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new)+cumulativeOffset, // rule_vectors
      nLaunchCells,
      baseStream
      );
   // Find Blocks (GID,LID) to be replaced with new ones
   extract_to_replace_caller(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+2*cumulativeOffset+1,//dev_has_no_content_maps, // input maps
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_to_replace)+cumulativeOffset, // output vecs
      NULL, // pass null to not store vector lengths
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes)+cumulativeOffset, // rule_meshes
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+2*cumulativeOffset+1,//dev_has_no_content_maps, // rule_maps
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new)+cumulativeOffset, // rule_vectors
      nLaunchCells,
      baseStream
      );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   extractTimer.stop();

   // Note: in this call, unless hitting v-space walls, we only grow the vspace size
   // and thus do not delete blocks or replace with old blocks.
   // The call now uses the batch block adjust interface.
   phiprof::Timer adjustTimer {"block adjust caller"};
   uint largestBlocksToChange; // Not needed
   uint largestBlocksBeforeOrAfter; // Not needed
   batch_adjust_blocks_caller(
      mpiGrid,
      launchCells,
      cumulativeOffset,
      largestBlocksToChange,
      largestBlocksBeforeOrAfter,
      popID);
   // This caller function updates values in host_nAfter
   // Velocity space has now all extra blocks added and/or removed for the transform target
   // and will not change shape anymore.
   adjustTimer.stop();

   // Track of largest vmesh size, evaluate launch parameters for zeroing kernel
   phiprof::Timer alloc2Timer {"ensure allocations 2"};
   for (size_t cellIndex = 0; cellIndex < nLaunchCells; cellIndex++) {
      SpatialCell* SC = mpiGrid[launchCells[cellIndex]];
      const uint cellOffset = cellIndex + cumulativeOffset;
      // The function batch_adjust_blocks_caller updates host_nAfter
      const vmesh::LocalID nBlocksAfterAdjust = (GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_nAfter))[cellOffset];
      SC->largestvmesh = SC->largestvmesh > nBlocksAfterAdjust ? SC->largestvmesh : nBlocksAfterAdjust;
      largest_nAfter = std::max(largest_nAfter,nBlocksAfterAdjust);
   } // end parallel region
   alloc2Timer.stop();

   /* Zero out target data on device (unified). We could call a separate gpuMemSet for each
      blockContainer, but gpuMemSets will just end up calling a kernel under the hood anyway, and
      with this kernel of our own we can use one call for all spatial cells at once, and the pointers
      to the velocity block containers already reside on-device.
    */
   phiprof::Timer zeroTimer {"zero target data"};
   const size_t n_fill_VBC_zero = 1 + ((largest_nAfter*WID3 - 1) / Hashinator::defaults::MAX_BLOCKSIZE);
   const dim3 grid_fill_VBC_zero(n_fill_VBC_zero,nLaunchCells,1);
   fill_VBC_zero_kernel<<<grid_fill_VBC_zero,Hashinator::defaults::MAX_BLOCKSIZE,0,baseStream>>>(
      GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs), // indexing: cellOffset
      cumulativeOffset
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   zeroTimer.stop();

   // Launch actual acceleration kernel performing Semi-Lagrangian re-mapping
   phiprof::Timer accTimer {"acceleration kernel"};
   const dim3 grid_acc(largest_totalColumnSets,nLaunchCells,1);
   const dim3 block_acc(WID,WID,WID); // Calculates a whole block at a time
   acceleration_kernel<<<grid_acc, block_acc, 0, baseStream>>> (
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), // indexing: cellOffset
      GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs), // indexing: cellOffset
      GET_POINTER(gpuMemoryManager, Realf*, dev_blockDataOrdered), //indexing: blockIdx.y
      GET_POINTER(gpuMemoryManager, uint, gpu_cell_indices_to_id),
      GET_POINTER(gpuMemoryManager, uint, gpu_block_indices_to_id),
      GET_POINTER(gpuMemoryManager, ColumnOffsets, dev_columnOffsetData), //indexing: blockIdx.y
      GET_POINTER(gpuMemoryManager, Realf, dev_intersections), // indexing: cellOffset
      v_min,
      i_dv,
      dv,
      GET_POINTER(gpuMemoryManager, Real, dev_minValues), // indexing: cellOffset, used by slope limiters
      invalidLocalID,
      cumulativeOffset
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   accTimer.stop();

   return true;
}
