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


#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "gpu_acc_sort_blocks.hpp"

using namespace std;
using namespace spatial_cell;

#ifdef DEBUG_VLASIATOR
   #ifndef DEBUG_ACC
   #define DEBUG_ACC
   #endif
#endif

__host__ void gpu_acc_allocate_radix_sort (
   const uint temp_storage_bytes,
   const uint cpuThreadID,
   const gpuStream_t stream
   ) {
   if (temp_storage_bytes * BLOCK_ALLOCATION_FACTOR > gpu_acc_RadixSortTempSize[cpuThreadID]) {
      if (gpu_acc_RadixSortTempSize[cpuThreadID] > 0) {
         CHK_ERR( gpuFreeAsync(gpu_RadixSortTemp[cpuThreadID], stream) );
      }
      gpu_acc_RadixSortTempSize[cpuThreadID] = temp_storage_bytes * BLOCK_ALLOCATION_PADDING;
      CHK_ERR( gpuMallocAsync((void**)&gpu_RadixSortTemp[cpuThreadID], gpu_acc_RadixSortTempSize[cpuThreadID], stream) );
   }
}
// This memory deallocated befor exit in gpu_clear_device() in arch/gpu_base.cpp

//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor, maxBlocksPerCluster)

// Kernels for converting GIDs to dimension-sorted indices
__global__ void __launch_bounds__(GPUTHREADS,4) blocksID_mapped_dim0_kernel(
   const vmesh::VelocityMesh* __restrict__ vmesh,
   vmesh::GlobalID *blocksID_mapped,
   vmesh::LocalID *blocksLID_unsorted,
   const uint nBlocks
   ) {
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   #ifdef DEBUG_ACC
   if (nBlocks != vmesh->size()) {
      printf("Error inside blocksID_mapped_dim0_kernel: nBlocks %u does not match vmesh size %lu!\n",nBlocks,vmesh->size());
   }
   #endif
   const vmesh::LocalID LID = blocki * warpSize + ti;
   if (LID < nBlocks) {
      blocksID_mapped[LID] = vmesh->getGlobalID(LID);
      blocksLID_unsorted[LID] = LID;
   }
}

__global__ void __launch_bounds__(GPUTHREADS,4) blocksID_mapped_dim1_kernel(
   const vmesh::VelocityMesh* __restrict__ vmesh,
   vmesh::GlobalID *blocksID_mapped,
   vmesh::LocalID *blocksLID_unsorted,
   const uint nBlocks
   ) {
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   #ifdef DEBUG_ACC
   if (nBlocks != vmesh->size()) {
      printf("Error inside blocksID_mapped_dim1_kernel: nBlocks %u does not match vmesh size %lu!\n",nBlocks,vmesh->size());
   }
   #endif
   const vmesh::LocalID D0 = vmesh->getGridLength()[0];
   const vmesh::LocalID D1 = vmesh->getGridLength()[1];
   // const vmesh::LocalID D2 = vmesh->getGridLength()[2];
   const vmesh::LocalID LID = blocki * warpSize + ti;
   if (LID < nBlocks) {
      const vmesh::GlobalID GID = vmesh->getGlobalID(LID);
      const vmesh::LocalID x_index = GID % D0;
      const vmesh::LocalID y_index = (GID / D0) % D1;
      blocksID_mapped[LID] = GID - (x_index + y_index*D0) + y_index + x_index * D1;
      blocksLID_unsorted[LID] = LID;
   }
}

__global__ void __launch_bounds__(GPUTHREADS,4) blocksID_mapped_dim2_kernel(
   const vmesh::VelocityMesh* __restrict__ vmesh,
   vmesh::GlobalID *blocksID_mapped,
   vmesh::LocalID *blocksLID_unsorted,
   const uint nBlocks
   ) {
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   #ifdef DEBUG_ACC
   if (nBlocks != vmesh->size()) {
      printf("Error inside blocksID_mapped_dim2_kernel: nBlocks %u does not match vmesh size %lu!\n",nBlocks,vmesh->size());
   }
   #endif
   const vmesh::LocalID D0 = vmesh->getGridLength()[0];
   const vmesh::LocalID D1 = vmesh->getGridLength()[1];
   const vmesh::LocalID D2 = vmesh->getGridLength()[2];
   const vmesh::LocalID LID = blocki * warpSize + ti;
   if (LID < nBlocks) {
      const vmesh::GlobalID GID = vmesh->getGlobalID(LID);
      const vmesh::LocalID x_index = GID % D0;
      const vmesh::LocalID y_index = (GID / D0) % D1;
      const vmesh::LocalID z_index = (GID / (D0*D1));
      blocksID_mapped[LID] = z_index + y_index*D2 + x_index*D1*D2;
      blocksLID_unsorted[LID] = LID;
   }
}

// LIDs are already in order.
// Now also order GIDS. (can be ridiculously parallel, minus memory access patterns)
// also used for scanning columnsets for block counts
__global__ void __launch_bounds__(GPUTHREADS,4) order_GIDs_kernel(
   const vmesh::VelocityMesh* __restrict__ vmesh,
   const vmesh::GlobalID* __restrict__ blocksLID,
   vmesh::GlobalID *blocksGID,
   const uint dimension,
   const vmesh::GlobalID* __restrict__ blocksID_mapped_sorted,
   vmesh::LocalID *gpu_columnNBlocks,
   const uint nBlocks,
   ColumnOffsets* columnData // passed just for resetting
   ) {
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   vmesh::LocalID DX;
   switch (dimension) {
      case 0:
         DX = vmesh->getGridLength()[0];
         break;
      case 1:
         DX = vmesh->getGridLength()[1];
         break;
      case 2:
         DX = vmesh->getGridLength()[2];
         break;
      default:
         printf("Incorrect dimension in __FILE__ __LINE__\n");
   }
   #ifdef DEBUG_ACC
   if (nBlocks != vmesh->size()) {
      printf("Error inside order_GIDs_kernel: nBlocks %u does not match vmesh size %lu!\n",nBlocks,vmesh->size());
   }
   #endif
   const vmesh::LocalID index = blocki * warpSize + ti;
   if (index < nBlocks) {
      #ifdef DEBUG_ACC
      if (blocksLID[index] >=  nBlocks ) {
         printf(" Too large LID %u!  (size %u)\n",blocksLID[index],nBlocks);
      }
      #endif
      blocksGID[index] = vmesh->getGlobalID(blocksLID[index]);

      const vmesh::LocalID column_id = blocksID_mapped_sorted[index] / DX;
      // Increment number of blocks in column
      const vmesh::LocalID old  = atomicAdd(&gpu_columnNBlocks[column_id],1);
      // // Evaluate smallest GID in column
      // old = atomicMin(&columnMinBlock[columnid],GID);
      // // Evaluate largest GID in colum
      // old = atomicMax(&columnMaxBlock[columnid],GID);
   }
   if (blockIdx.x == blockIdx.y == blockIdx.z == threadIdx.x == threadIdx.y == threadIdx.z == 0) {
      columnData->columnBlockOffsets.clear();
      columnData->columnNumBlocks.clear();
      columnData->setColumnOffsets.clear();
      columnData->setNumColumns.clear();
   }
}

/*** Kernel for constructing columns
 Checks if all blocks in a columnset belong to a single column, and
 can quickly jump through the whole column.
 For columnsets containing several columns, it trials blocks by scanning
 warpSize GIDs at a time.

 Still probably room for memory optimization.
**/

__global__ void __launch_bounds__(GPUTHREADS,4) construct_columns_kernel(
   const vmesh::VelocityMesh* __restrict__ vmesh,
   const uint dimension,
   const vmesh::GlobalID* __restrict__ blocksID_mapped_sorted,
   const vmesh::LocalID* __restrict__ gpu_columnNBlocks,
   ColumnOffsets* columnData,
   const uint nBlocks
   ) {
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   //const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   #ifdef DEBUG_ACC
   const int gpuBlocks = gridDim.x * gridDim.y * gridDim.z;
   if (gpuBlocks!=1) {
      printf("Error in construct_columns_kernel; unsafe gridDim\n");
      return;
   }
   #endif
   vmesh::LocalID DX;
   switch (dimension) {
      case 0:
         DX = vmesh->getGridLength()[0];
         break;
      case 1:
         DX = vmesh->getGridLength()[1];
         break;
      case 2:
         DX = vmesh->getGridLength()[2];
         break;
      default:
         printf("Incorrect dimension in __FILE__ __LINE__\n");
   }
   vmesh::LocalID prev_column_id, prev_dimension_id;

   __shared__ vmesh::LocalID i;
   __shared__ vmesh::LocalID blocks_in_columnset;
   if (ti==0) {
      i = 0;
      blocks_in_columnset = 0;
      // Put in the sorted blocks, and also compute column offsets and lengths:
      columnData->columnBlockOffsets.device_push_back(0); //first offset
      columnData->setColumnOffsets.device_push_back(0); //first offset
   }
   __syncthreads();
   while (i < nBlocks) {
      // identifies a particular column
      const vmesh::LocalID column_id = blocksID_mapped_sorted[i] / DX;
      // identifies a particular block in a column (along the dimension)
      const vmesh::LocalID dimension_id = blocksID_mapped_sorted[i] % DX;
      // How many blocks in this (new) column(set)?
      if ((ti==0) && (blocks_in_columnset==0)) {
         blocks_in_columnset = gpu_columnNBlocks[column_id];
      }
      // Trial: new column?
      if ( (ti==0) && (i > 0) &&  ( (column_id != prev_column_id) || (dimension_id != (prev_dimension_id + 1) ))) {
         //encountered new column! For i=0, we already entered the correct offset (0).
         //We also identify it as a new column if there is a break in the column (e.g., gap between two populations)
         //add offset where the next column will begin
         columnData->columnBlockOffsets.device_push_back(i);
         //add length of the current column that now ended
         columnData->columnNumBlocks.device_push_back(columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-1] - columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-2]);

         if (column_id != prev_column_id ){
            //encountered new set of columns, add offset to new set starting at present column
            columnData->setColumnOffsets.device_push_back(columnData->columnBlockOffsets.size() - 1);
            //add length of the previous column set that ended
            columnData->setNumColumns.device_push_back(columnData->setColumnOffsets[columnData->setColumnOffsets.size()-1] - columnData->setColumnOffsets[columnData->setColumnOffsets.size()-2]);
         }
      }
      __syncthreads();
      // Trial if only one column in columnset?
      if ( ( (blocksID_mapped_sorted[i+blocks_in_columnset-1] % DX) == (dimension_id + blocks_in_columnset - 1) ) &&
           ( (blocksID_mapped_sorted[i+blocks_in_columnset-1] / DX) == column_id ) ) {
         // skip to end of column
         if (ti==0) {
            i += blocks_in_columnset;
            blocks_in_columnset = 0;
         }
         // push_back to vectors happens at next loop
         __syncthreads();
      } else {
         // More than one column in columnset
         vmesh::LocalID this_col_length = 0;
         // Now trial by warpSize to see where column ends
         for (vmesh::LocalID ci=0; ci<blocks_in_columnset; ci += warpSize) {
            int notInColumn = 1;
            if (ci+ti < blocks_in_columnset) {
               // This evaluates if the block at the target point is no longer within the same column
               if ( (blocksID_mapped_sorted[i+ci+ti] % DX) == (dimension_id + ci+ti) &&
                    ( (blocksID_mapped_sorted[i+ci+ti] / DX) == column_id ) ) {
                  notInColumn = 0;
               }
            }
            // Warp vote to find first index (potentially) outside old column            
            //unsigned ballot_result = __ballot_sync(FULL_MASK, notInColumn);
            unsigned ballot_result = gpuKernelBallot(FULL_MASK, notInColumn);
            vmesh::LocalID minstep = __ffs(ballot_result); // Find first significant
             if (minstep==0) {
               minstep +=32; // no value found, jump whole warpSize
            } else {
               minstep--; // give actual index
            }
            this_col_length += minstep;
            // Exit this for loop if we reached the end of a column within a set
            if (minstep!=warpSize) {
               if (ti==0) {
                  // skip to end of column
                  i += this_col_length;
                  // Decrease number of free blocks in column(set)
                  blocks_in_columnset -= this_col_length;
               }
               __syncthreads();
               // push_back to vectors happens at next loop
               break;
            }
         }
      }
      prev_column_id = column_id;
      prev_dimension_id = dimension_id;
      __syncthreads();
   }
   // Add offsets for final column
   if (ti==0) {
      columnData->columnNumBlocks.device_push_back(nBlocks - columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-1]);
      columnData->setNumColumns.device_push_back(columnData->columnNumBlocks.size() - columnData->setColumnOffsets[columnData->setColumnOffsets.size()-1]);
   }
}

/*
   This function returns a sorted list of blocks in a cell.

   The sorted list is sorted according to the location, along the given dimension.
   This version uses triplets internally and also returns the LIDs of the sorted blocks.
*/
void sortBlocklistByDimension( vmesh::VelocityMesh* vmesh, //on-device vmesh
                               const vmesh::LocalID nBlocks,
                               const uint dimension,
                               vmesh::GlobalID *blocksID_mapped,
                               vmesh::GlobalID *blocksID_mapped_sorted,
                               vmesh::GlobalID *blocksGID,
                               vmesh::LocalID *blocksLID_unsorted,
                               vmesh::LocalID *blocksLID,
                               vmesh::LocalID *gpu_columnNBlocks,
                               ColumnOffsets* columnData,
                               const uint cpuThreadID,
                               gpuStream_t stream
   ) {

   // Ensure at least one launch block, ceil int division
   //const uint maxThreads = WARPSPERBLOCK*GPUTHREADS;
   // For some reason, a direct increase of launch threads breaks things. Stick with GPUTHREADS for now.
   const uint maxThreads = GPUTHREADS;
   const uint launchBlocks = 1 + ((nBlocks - 1) / (maxThreads));

   #ifdef DEBUG_ACC
   CHK_ERR( gpuMemsetAsync(blocksGID, 0, nBlocks*sizeof(vmesh::GlobalID), stream) );
   CHK_ERR( gpuMemsetAsync(blocksID_mapped, 0, nBlocks*sizeof(vmesh::GlobalID), stream) );
   CHK_ERR( gpuMemsetAsync(blocksID_mapped_sorted, 0, nBlocks*sizeof(vmesh::GlobalID), stream) );
   CHK_ERR( gpuMemsetAsync(blocksLID_unsorted, 0, nBlocks*sizeof(vmesh::LocalID), stream) );
   CHK_ERR( gpuMemsetAsync(blocksLID, 0, nBlocks*sizeof(vmesh::LocalID), stream) );
   #endif

   // Map blocks to new dimensionality
   switch( dimension ) {
      case 0: {
         blocksID_mapped_dim0_kernel<<<launchBlocks, maxThreads, 0, stream>>> (
            vmesh,
            blocksID_mapped,
            blocksLID_unsorted,
            nBlocks
            );
         break;
      }
      case 1: {
         blocksID_mapped_dim1_kernel<<<launchBlocks, maxThreads, 0, stream>>> (
            vmesh,
            blocksID_mapped,
            blocksLID_unsorted,
            nBlocks
            );
         break;
      }
      case 2: {
         blocksID_mapped_dim2_kernel<<<launchBlocks, maxThreads, 0, stream>>> (
            vmesh,
            blocksID_mapped,
            blocksLID_unsorted,
            nBlocks
            );
         break;
      }
      default:
         printf("Incorrect dimension in gpu_acc_sort_blocks.cpp\n");
   }
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(stream) );
   //SSYNC;

   // Determine temporary device storage requirements
   void     *temp_storage_null = NULL;
   size_t   temp_storage_bytes = 0;
   //GPUTODO: HIPIFY-option via arch
   #ifdef __CUDACC__
   cub::DeviceRadixSort::SortPairs(temp_storage_null, temp_storage_bytes,
                                   blocksID_mapped, blocksID_mapped_sorted,
                                   blocksLID_unsorted, blocksLID, nBlocks,
                                   0, sizeof(vmesh::GlobalID)*8, stream);
   #else
   // HIPCUB
   hipcub::DeviceRadixSort::SortPairs(temp_storage_null, temp_storage_bytes,
                                   blocksID_mapped, blocksID_mapped_sorted,
                                   blocksLID_unsorted, blocksLID, nBlocks,
                                   0, sizeof(vmesh::GlobalID)*8, stream);
   #endif
   CHK_ERR( gpuPeekAtLastError() );

   gpu_acc_allocate_radix_sort(temp_storage_bytes,cpuThreadID,stream);

   // Now sort
   #ifdef __CUDACC__
   cub::DeviceRadixSort::SortPairs(gpu_RadixSortTemp[cpuThreadID], temp_storage_bytes,
                                   blocksID_mapped, blocksID_mapped_sorted,
                                   blocksLID_unsorted, blocksLID, nBlocks,
                                   0, sizeof(vmesh::GlobalID)*8, stream);
   #else
   // HIPCUB
   hipcub::DeviceRadixSort::SortPairs(gpu_RadixSortTemp[cpuThreadID], temp_storage_bytes,
                                   blocksID_mapped, blocksID_mapped_sorted,
                                   blocksLID_unsorted, blocksLID, nBlocks,
                                   0, sizeof(vmesh::GlobalID)*8, stream);
   #endif
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(stream) );

   // Gather GIDs in order
   order_GIDs_kernel<<<launchBlocks, maxThreads, 0, stream>>> (
      vmesh,
      blocksLID,
      blocksGID,
      dimension,
      blocksID_mapped,
      gpu_columnNBlocks,
      nBlocks,
      columnData // Pass this just to clear it on device
      );
   CHK_ERR( gpuPeekAtLastError() );

   // Construct columns. To ensure order,
   // these are done serially, but still form within a kernel.
   // Optimizing launch threads for this kernel as well needs a bit more checking
   construct_columns_kernel<<<1, GPUTHREADS, 0, stream>>> (
      vmesh,
      dimension,
      blocksID_mapped_sorted,
      gpu_columnNBlocks,
      columnData,
      nBlocks
      );
   CHK_ERR( gpuPeekAtLastError() );
   //SSYNC;

}
