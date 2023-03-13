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


#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "cuda_acc_sort_blocks.hpp"

#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/device_ptr.h>

// Ensure printing of CUDA runtime errors to console
// #define CUB_STDERR
// #include <cub/util_allocator.cuh>
// #include <cub/device/device_radix_sort.cuh>

using namespace std;
using namespace spatial_cell;

// Comparator function for sorting vector of pairs
inline bool tripletcomparator( const std::pair<std::pair<uint, uint>, uint> & l, const std::pair<std::pair<uint, uint>, uint> & r ) {
   return l.first.first < r.first.first;
}

// Kernels for converting GIDs to dimension-sorted indices
__global__ void blocksID_mapped_dim0_kernel(
   const vmesh::VelocityMesh* vmesh,
   vmesh::GlobalID *blocksID_mapped,
   vmesh::LocalID *blocksLID,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z; 
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (vmesh::LocalID index=blocki*warpSize; index<nBlocks; index += cudaBlocks*warpSize) {
      const vmesh::LocalID LID = (index+ti);
      if (LID < nBlocks) {
         blocksID_mapped[LID] = vmesh->getGlobalID(LID);
         blocksLID[LID]=LID;
      }
   }
}

__global__ void blocksID_mapped_dim1_kernel(
   const vmesh::VelocityMesh* vmesh,
   vmesh::GlobalID *blocksID_mapped,
   vmesh::LocalID *blocksLID,
   vmesh::LocalID D0,
   vmesh::LocalID D1,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z; 
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (vmesh::LocalID index=blocki*warpSize; index<nBlocks; index += cudaBlocks*warpSize) {
      const vmesh::LocalID LID = (index+ti);
      if (LID < nBlocks) {
         const vmesh::GlobalID GID = vmesh->getGlobalID(LID);
         const vmesh::LocalID x_index = GID % D0;
         const vmesh::LocalID y_index = (GID / D0) % D1;
         blocksID_mapped[LID] = GID - (x_index + y_index*D0) + y_index + x_index * D1;
         blocksLID[LID]=LID;
      }
   }
}

__global__ void blocksID_mapped_dim2_kernel(
   const vmesh::VelocityMesh* vmesh,
   vmesh::GlobalID *blocksID_mapped,
   vmesh::LocalID *blocksLID,
   vmesh::LocalID D0,
   vmesh::LocalID D1,
   vmesh::LocalID D2,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z; 
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (vmesh::LocalID index=blocki*warpSize; index<nBlocks; index += cudaBlocks*warpSize) {
      const vmesh::LocalID LID = (index+ti);
      if (LID < nBlocks) {
         const vmesh::GlobalID GID = vmesh->getGlobalID(LID);
         const vmesh::LocalID x_index = GID % D0;
         const vmesh::LocalID y_index = (GID / D0) % D1;
         const vmesh::LocalID z_index = (GID / (D0*D1));
         blocksID_mapped[LID] = z_index + y_index*D2 + x_index*D1*D2;
         blocksLID[LID]=LID;
      }
   }
}

// LIDs are already in order.
// Now also order GIDS. (can be ridiculously parallel, minus memory access patterns)
__global__ void order_GIDs_kernel(
   const vmesh::VelocityMesh* vmesh,
   vmesh::GlobalID *blocksLID,
   vmesh::GlobalID *blocksGID,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z; 
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   for (vmesh::LocalID index=blocki*warpSize; index<nBlocks; index += cudaBlocks*warpSize) {
      const vmesh::LocalID i = (index+ti);
      if (i < nBlocks) {
         blocksGID[i]=vmesh->getGlobalID(blocksLID[i]);
      }
   }
}

__global__ void construct_columns_kernel(
   const vmesh::VelocityMesh* vmesh,
   vmesh::GlobalID *blocksID_mapped,
   ColumnOffsets* columnData,
   vmesh::LocalID DX,
   const uint nBlocks
   ) {
   const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z; 
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   // const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   // const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   if ((cudaBlocks!=1) || (warpSize!=1)) {
      printf("Error in construct_columns_kernel; unsafe warpSize or gridDim\n");
      return;
   }

   // Put in the sorted blocks, and also compute column offsets and lengths:
   columnData->columnBlockOffsets.device_push_back(0); //first offset
   columnData->setColumnOffsets.device_push_back(0); //first offset
   vmesh::LocalID prev_column_id, prev_dimension_id;
   
   for (vmesh::LocalID i=0; i<nBlocks; ++i) {
      // identifies a particular column
      vmesh::LocalID column_id = blocksID_mapped[i] / DX;
      // identifies a particular block in a column (along the dimension)
      vmesh::LocalID dimension_id = blocksID_mapped[i] % DX;

      if ( i > 0 &&  ( (column_id != prev_column_id) || (dimension_id != (prev_dimension_id + 1) ))) {
         //encountered new column! For i=0, we already entered the correct offset (0).
         //We also identify it as a new column if there is a break in the column (e.g., gap between two populations)
         /*add offset where the next column will begin*/
         columnData->columnBlockOffsets.device_push_back(i);
         /*add length of the current column that now ended*/
         columnData->columnNumBlocks.device_push_back(columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-1] - columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-2]);

         if (column_id != prev_column_id ){
            //encountered new set of columns, add offset to new set starting at present column
            columnData->setColumnOffsets.device_push_back(columnData->columnBlockOffsets.size() - 1);
            /*add length of the previous column set that ended*/
            columnData->setNumColumns.device_push_back(columnData->setColumnOffsets[columnData->setColumnOffsets.size()-1] - columnData->setColumnOffsets[columnData->setColumnOffsets.size()-2]);
         }
      }
      prev_column_id = column_id;
      prev_dimension_id = dimension_id;
   }

   columnData->columnNumBlocks.device_push_back(nBlocks - columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-1]);
   columnData->setNumColumns.device_push_back(columnData->columnNumBlocks.size() - columnData->setColumnOffsets[columnData->setColumnOffsets.size()-1]);
}

/*
   This function returns a sorted list of blocks in a cell.

   The sorted list is sorted according to the location, along the given dimension.
   This version uses triplets internally and also returns the LIDs of the sorted blocks.
*/
void sortBlocklistByDimension( //const spatial_cell::SpatialCell* spatial_cell,
                               const vmesh::VelocityMesh* vmesh,
                               const uint dimension,
                               vmesh::GlobalID *blocksGID,
                               vmesh::GlobalID *blocksID_mapped,
                               vmesh::LocalID *blocksLID,
                               ColumnOffsets* columnData,
   // split::SplitVector<uint> columnBlockOffsets; // indexes where columns start (in blocks, length totalColumns)
   // split::SplitVector<uint> columnNumBlocks; // length of column (in blocks, length totalColumns)
   // split::SplitVector<uint> setColumnOffsets; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
   // split::SplitVector<uint> setNumColumns; // how many columns in set of columns (length nColumnSets)
                               const uint cuda_async_queue_id,
                               cudaStream_t stream
   ) {

   columnData->columnBlockOffsets.clear();
   columnData->columnNumBlocks.clear();
   columnData->setColumnOffsets.clear();
   columnData->setNumColumns.clear();

#ifndef UNDEFINED
   //thrust::sort(thrust::device, data, data + size);
   // thrust::sort_by_key(keys, keys + N, values);
   // auto dptr = thrust::device_ptr<uint64_t>(data);
   // thrust::sort(dptr , dptr + size);
   // Could use CUB instead?

   const vmesh::LocalID nBlocks = vmesh->size();
   const uint refL=0; //vAMR
   const vmesh::LocalID D0 = vmesh->getGridLength(refL)[0];
   const vmesh::LocalID D1 = vmesh->getGridLength(refL)[1];
   const vmesh::LocalID D2 = vmesh->getGridLength(refL)[2];
   vmesh::LocalID DX[3] = {D0,D1,D2};
   
   uint nCudaBlocks  = (nBlocks/CUDATHREADS) > CUDABLOCKS ? CUDABLOCKS : (nBlocks/CUDATHREADS);
   dim3 grid(nCudaBlocks,1,1);
   dim3 block(CUDATHREADS,1,1);
   // thrust stream set already here
   // thrust::cuda::par.on(stream)
   
   // Map blocks to new dimensionality
   switch( dimension ) {
      case 0: {
         blocksID_mapped_dim0_kernel<<<grid, block, 0, stream>>> (
            vmesh,
            blocksID_mapped,
            blocksLID,
            nBlocks
            );
         cudaStreamSynchronize(stream);
         thrust::device_ptr<vmesh::GlobalID> dptrMapped(blocksID_mapped);
         thrust::device_ptr<vmesh::LocalID> dptrLID(blocksLID);
         thrust::sort_by_key(dptrMapped, dptrMapped + nBlocks, dptrLID);
         //cudaStreamSynchronize(stream);
         cudaDeviceSynchronize();
         break;
      }
      case 1: {
         blocksID_mapped_dim1_kernel<<<grid, block, 0, stream>>> (
            vmesh,
            blocksID_mapped,
            blocksLID,
            D0,D1,
            nBlocks
            );
         cudaStreamSynchronize(stream);
         // Sort (with thrust)
         // vmesh::GlobalID *dptrMapped = thrust::device_ptr<vmesh::GlobalID>(blocksID_mapped);
         // vmesh::LocalID *dptrLID = thrust::device_ptr<vmesh::LocalID>(blocksLID);
         //thrust::sort_by_key(thrust::cuda::par(stream), dptrMapped, dptrMapped + nBlocks, dptrLID);
         thrust::device_ptr<vmesh::GlobalID> dptrMapped(blocksID_mapped);
         thrust::device_ptr<vmesh::LocalID> dptrLID(blocksLID);
         thrust::sort_by_key(dptrMapped, dptrMapped + nBlocks, dptrLID);
         //cudaStreamSynchronize(stream);
         cudaDeviceSynchronize();
         // static CUB_RUNTIME_FUNCTION cudaError_t cub::DeviceRadixSort::SortPairs
         //    ( void *  d_temp_storage,
         //      size_t &  temp_storage_bytes,
         //      const KeyT *  d_keys_in,
         //      KeyT *  d_keys_out,
         //      const ValueT *  d_values_in,
         //      ValueT *  d_values_out,
         //      NumItemsT  num_items,
         //      int  begin_bit = 0,
         //      int  end_bit = sizeof(KeyT) * 8,
         //      cudaStream_t  stream = 0,
         //      bool  debug_synchronous = false
         //       )            
         break;
      }
      case 2: {
         blocksID_mapped_dim2_kernel<<<grid, block, 0, stream>>> (
            vmesh,
            blocksID_mapped,
            blocksLID,
            D0,D1,D2,
            nBlocks
            );
         cudaStreamSynchronize(stream);
         // Sort (with thrust)
         // vmesh::GlobalID *dptrMapped = thrust::device_ptr<vmesh::GlobalID>(blocksID_mapped);
         // vmesh::LocalID *dptrLID = thrust::device_ptr<vmesh::LocalID>(blocksLID);
         //thrust::sort_by_key(thrust::cuda::par(stream), dptrMapped, dptrMapped + nBlocks, dptrLID);
         thrust::device_ptr<vmesh::GlobalID> dptrMapped(blocksID_mapped);
         thrust::device_ptr<vmesh::LocalID> dptrLID(blocksLID);
         thrust::sort_by_key(dptrMapped, dptrMapped + nBlocks, dptrLID);
         //cudaStreamSynchronize(stream);
         cudaDeviceSynchronize();
         break;
      }
      default:
         printf("Incorrect dimension in cuda_acc_sort_blocks.cpp\n");
   }
   // Gather GIDs in order
   order_GIDs_kernel<<<grid, block, 0, stream>>> (
      vmesh,
      blocksLID,
      blocksGID,
      nBlocks
      );
   // Parallel with this, we can also construct columns. To ensure order,
   // these are done serially, but still form within a kernel.
   dim3 grid1(1,1,1);
   dim3 block1(1,1,1);
   printf("Start construct columns kernel\n");
   construct_columns_kernel<<<grid1, block1, 0, stream>>> (
      vmesh,
      blocksID_mapped,
      columnData,
      DX[dimension],
      nBlocks
      );

   // And then we synchronize
   cudaStreamSynchronize(stream);
   
#else

   //const uint nBlocks = spatial_cell->get_number_of_velocity_blocks(); // Number of blocks
   const vmesh::LocalID nBlocks = vmesh->size();

   // Velocity mesh refinement level, has no effect here
   // but is needed in some vmesh::VelocityMesh function calls.
   const uint8_t REFLEVEL = 0;
   // Copy block data to vector
   std::vector<std::pair<std::pair<vmesh::GlobalID,vmesh::GlobalID>,vmesh::LocalID> > block_triplets;
   block_triplets.resize( nBlocks );
   for (vmesh::LocalID i = 0; i < nBlocks; ++i ) {
      //const vmesh::GlobalID block = spatial_cell->get_velocity_block_global_id(i);
      const vmesh::GlobalID block = vmesh->getGlobalID(i);
      switch( dimension ) {
       case 0: {
          const vmesh::GlobalID blockId_mapped = block; // Mapping the block id to different coordinate system if dimension is not zero:
          block_triplets[i] = std::make_pair( std::make_pair( blockId_mapped, block ), i);
       }
         break;
       case 1: {
          // Do operation:
          //   block = x + y*x_max + z*y_max*x_max
          //=> block' = block - (x + y*x_max) + y + x*y_max = x + y*x_max + z*y_max*x_max - (x + y*x_max) + y + x*y_max
          //          = y + x*y_max + z*y_max*x_max
          //const uint x_indice = block%SpatialCell::get_velocity_grid_length()[0];
          //const uint y_indice = (block/SpatialCell::get_velocity_grid_length()[0])%SpatialCell::SpatialCell::get_velocity_grid_length()[1];
          const vmesh::LocalID x_index = block % vmesh->getGridLength(REFLEVEL)[0];
          const vmesh::LocalID y_index = (block / vmesh->getGridLength(REFLEVEL)[0]) % vmesh->getGridLength(REFLEVEL)[1];

          // Mapping the block id to different coordinate system if dimension is not zero:
          //const uint blockId_mapped
          //        = block - (x_indice + y_indice*SpatialCell::get_velocity_grid_length()[0])
          //        + y_indice
          //        + x_indice * SpatialCell::SpatialCell::get_velocity_grid_length()[1];
          const vmesh::GlobalID blockId_mapped
                  = block - (x_index + y_index*vmesh->getGridLength(REFLEVEL)[0])
                  + y_index
                  + x_index * vmesh->getGridLength(REFLEVEL)[1];
          block_triplets[i] = std::make_pair( std::make_pair( blockId_mapped, block ), i);
       }
         break;
       case 2: {
          // Do operation:
          //   block = x + y*x_max + z*y_max*x_max
          //=> block' = z + y*z_max + x*z_max*y_max
          //const uint x_indice = block%SpatialCell::get_velocity_grid_length()[0];
          //const uint y_indice = (block/SpatialCell::get_velocity_grid_length()[0])%SpatialCell::SpatialCell::get_velocity_grid_length()[1];
          //const uint z_indice =  (block/(SpatialCell::get_velocity_grid_length()[0]*SpatialCell::SpatialCell::get_velocity_grid_length()[1]));
          const vmesh::LocalID x_index = block % vmesh->getGridLength(REFLEVEL)[0];
          const vmesh::LocalID y_index = (block / vmesh->getGridLength(REFLEVEL)[0]) % vmesh->getGridLength(REFLEVEL)[1];
          const vmesh::LocalID z_index = (block / (vmesh->getGridLength(REFLEVEL)[0]*vmesh->getGridLength(REFLEVEL)[1]));

          // Mapping the block id to different coordinate system if dimension is not zero:
          //const uint blockId_mapped
          //  = z_indice
          //  + y_indice * SpatialCell::SpatialCell::get_velocity_grid_length()[2]
          //  + x_indice*SpatialCell::SpatialCell::get_velocity_grid_length()[1]*SpatialCell::SpatialCell::get_velocity_grid_length()[2];
          const vmesh::GlobalID blockId_mapped
            = z_index
            + y_index*vmesh->getGridLength(REFLEVEL)[2]
            + x_index*vmesh->getGridLength(REFLEVEL)[1]*vmesh->getGridLength(REFLEVEL)[2];
          block_triplets[i] = std::make_pair(std::make_pair( blockId_mapped, block ), i);
       }
         break;
      }
   }
   // Sort the list:
   std::sort( block_triplets.begin(), block_triplets.end(), tripletcomparator );

   // Put in the sorted blocks, and also compute column offsets and lengths:
   columnData->columnBlockOffsets.push_back(0); //first offset
   columnData->setColumnOffsets.push_back(0); //first offset
   uint prev_column_id, prev_dimension_id;

   for (vmesh::LocalID i=0; i<nBlocks; ++i) {
       // identifies a particular column
      vmesh::LocalID column_id = block_triplets[i].first.first / vmesh->getGridLength(REFLEVEL)[dimension];

       // identifies a particular block in a column (along the dimension)
       vmesh::LocalID dimension_id = block_triplets[i].first.first % vmesh->getGridLength(REFLEVEL)[dimension];

       //sorted lists
       blocksGID[i] = block_triplets[i].first.second;
       blocksLID[i] = vmesh->getLocalID(block_triplets[i].first.second);

       if ( i > 0 &&  ( column_id != prev_column_id || dimension_id != (prev_dimension_id + 1) )){
         //encountered new column! For i=0, we already entered the correct offset (0).
         //We also identify it as a new column if there is a break in the column (e.g., gap between two populations)
         /*add offset where the next column will begin*/
         columnData->columnBlockOffsets.push_back(i);
         /*add length of the current column that now ended*/
         columnData->columnNumBlocks.push_back(columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-1] - columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-2]);

         if (column_id != prev_column_id ){
            //encountered new set of columns, add offset to new set starting at present column
            columnData->setColumnOffsets.push_back(columnData->columnBlockOffsets.size() - 1);
            /*add length of the previous column set that ended*/
            columnData->setNumColumns.push_back(columnData->setColumnOffsets[columnData->setColumnOffsets.size()-1] - columnData->setColumnOffsets[columnData->setColumnOffsets.size()-2]);
         }
      }
      prev_column_id = column_id;
      prev_dimension_id = dimension_id;
   }

   columnData->columnNumBlocks.push_back(nBlocks - columnData->columnBlockOffsets[columnData->columnBlockOffsets.size()-1]);
   columnData->setNumColumns.push_back(columnData->columnNumBlocks.size() - columnData->setColumnOffsets[columnData->setColumnOffsets.size()-1]);

#endif
    printf("\n Output for dimension %d ",dimension);
    printf("\nColumnBlockOffsets %d\n", columnData->columnBlockOffsets.size());
   // for (auto i : columnData->columnBlockOffsets) printf("%d ",i);
    printf("\ncolumnNumBlocks %d\n", columnData->columnNumBlocks.size());
   // for (auto i : columnData->columnNumBlocks) printf("%d ",i);
    printf("\nsetColumnOffsets %d\n", columnData->setColumnOffsets.size());
   // for (auto i : columnData->setColumnOffsets) printf("%d ",i);
    printf("\nsetNumColumns %d\n", columnData->setNumColumns.size());
   // for (auto i : columnData->setNumColumns) printf("%d ",i);
    printf("\n \n",dimension);

}
