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

using namespace std;
using namespace spatial_cell;

// Comparator function for sorting vector of pairs
inline bool tripletcomparator( const std::pair<std::pair<uint, uint>, uint> & l, const std::pair<std::pair<uint, uint>, uint> & r ) {
   return l.first.first < r.first.first;
}

// // Kernel for vector reset
// __global__ void reset_column_search_vectors_kernel(
//    const vmesh::VelocityMesh* vmesh,
//    split::SplitVector<vmesh::LocalID> *columnNBlocks,
//    split::SplitVector<vmesh::LocalID> *columnMinBlock,
//    split::SplitVector<vmesh::LocalID> *columnMaxBlock
//    ) {
//    const int cudaBlocks = gridDim.x; //assumes 1D grid
//    const int blocki = blockIdx.x;
//    const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
//    const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

//    const vmesh::LocalID vl1 = columnNBlocks->size();
//    for (vmesh::LocalID vi=blocki*warpSize; vi<vl1; vi += cudaBlocks*warpSize) {
//       if (vi+ti<vl1) {
//          columnNBlocks->at(vi)=0;
//       }
//    }
//    const vmesh::LocalID vl2 = columnMinBlock->size();
//    for (vmesh::LocalID vi=blocki*warpSize; vi<vl2; vi += cudaBlocks*warpSize) {
//       if (vi+ti<vl2) {
//          columnMinBlock->at(vi)=vmesh->invalidLocalID();
//       }
//    }
//    const vmesh::LocalID vl3 = columnMaxBlock->size();
//    for (vmesh::LocalID vi=blocki*warpSize; vi<vl3; vi += cudaBlocks*warpSize) {
//       if (vi+ti<vl3) {
//          columnMaxBlock->at(vi)=0;
//       }
//    }
// }

// Kernel for scanning colums
__global__ void scan_blocks_for_columns_kernel(
   const vmesh::VelocityMesh* vmesh,
   const uint dimension,
   vmesh::LocalID *columnNBlocks,
   vmesh::LocalID *columnMinBlock,
   vmesh::LocalID *columnMaxBlock,
   const vmesh::LocalID D0,
   const vmesh::LocalID D1,
   const vmesh::LocalID D2
   ) {
   const int cudaBlocks = gridDim.x; //assumes 1D grid
   const int blocki = blockIdx.x;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   const uint nBlocks = vmesh->size();
   for (vmesh::LocalID LID=blocki*warpSize; LID<nBlocks; LID += cudaBlocks*warpSize) {
      if (LID+ti < nBlocks) {
         vmesh::GlobalID GID = vmesh->getGlobalID(LID+ti);
         // Solve column for GID
         vmesh::LocalID ind0,ind1,ind2;
         vmesh->getIndices(GID,ind0,ind1,ind2);
         uint columnid;
         switch (dimension) {
            case 0:
               columnid = ind1 + D1*ind2;
               break;
            case 1:
               columnid = ind0 + D1*ind2;
               break;
            case 2:
               columnid = ind0 + D1*ind1;
               break;
            default:
               printf(" Invalid dimension inside kernel!\n");
               return;
         }
         // Increment number of blocks in column
         unsigned long int old = atomicAdd(&columnNBlocks[columnid],1);
         // Evaluate smallest GID in column
         old = atomicMin(&columnMinBlock[columnid],GID);
         // Evaluate largest GID in colum
         old = atomicMax(&columnMaxBlock[columnid],GID);
      }
   }
}

// Kernel for sorting blocks
__global__ void sort_blocks_by_dimension_kernel(
   const vmesh::VelocityMesh* vmesh,
   const uint dimension,
   vmesh::LocalID *columnNBlocks,
   vmesh::LocalID *columnMinBlock,
   vmesh::LocalID *columnMaxBlock,
   const vmesh::LocalID D0,
   const vmesh::LocalID D1,
   const vmesh::LocalID D2,
   uint* blocksGID,
   uint* blocksLID,
   ColumnOffsets* columnData,
   vmesh::LocalID *dev_ColumnIndex,
   vmesh::LocalID *dev_ColumnSetIndex,
   vmesh::LocalID *dev_ColumnBlockIndex
   // columnData contains
   // split::SplitVector<uint> columnBlockOffsets; // indexes where columns start (in blocks, length totalColumns)
   // split::SplitVector<uint> columnNumBlocks; // length of column (in blocks, length totalColumns)
   // split::SplitVector<uint> setColumnOffsets; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
   // split::SplitVector<uint> setNumColumns; // how many columns in set of columns (length nColumnSets)
   ) {

   const int cudaBlocks = gridDim.x; //assumes 1D grid
   const int blocki = blockIdx.x;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   if ((warpSize!=1) || (ti!=0)) {
      printf("Error, kernel at __FILE__ __LINE__ called with unsupported launch size!\n");
   }
   // Only supports 1 thread for now
   uint stride;
   switch (dimension) {
      case 0: // X-direction, stride 1
         stride=1;
         break;
      case 1: // Y-direction, stride is gridsize[0]==D2
         stride=D2;
         break;
      case 2: // Z-direction, stride is gridsize[0]*gridsize[1] = D1*D2
         stride=D1*D2;
         break;
      default:
         printf(" Invalid dimension inside kernel!\n");
         return;
   }

   // Find own column from blocki, D1..2 and cudaBlocks
   const vmesh::LocalID Ctot = D1*D2;
   for (vmesh::LocalID Cind=blocki; Cind<Ctot; Cind += cudaBlocks) {
      if (columnNBlocks[Cind]==0) {
         // no active blocks along this column
         continue;
      }
      if (columnMinBlock[Cind] + stride*(columnNBlocks[Cind]-1) == columnMaxBlock[Cind]) {
         // Simple case: set contains only one column
         const uint ColumnSetIndex = atomicAdd(dev_ColumnSetIndex,1);
         const uint ColumnIndex = atomicAdd(dev_ColumnIndex,1);
         const uint BlockIndex = atomicAdd(dev_ColumnBlockIndex,columnNBlocks[Cind]);
         columnData->setNumColumns[ColumnSetIndex] = 1; // how many columns in set of columns (length nColumnSets)
         columnData->columnNumBlocks[ColumnIndex] = columnNBlocks[Cind]; // length of column (in blocks, length totalColumns)
         columnData->columnBlockOffsets[ColumnIndex] = BlockIndex; // indexes where columns start (in blocks, length totalColumns)
         columnData->setColumnOffsets[ColumnSetIndex]= BlockIndex; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
         // Now we proceed to place the GID,LID pairs into the allocated space
         for (uint bi = 0; bi < columnNBlocks[Cind]; bi++) {
            const vmesh::GlobalID GID = columnMinBlock[Cind] + bi * stride;
            blocksGID[BlockIndex+bi] = GID;
            blocksLID[BlockIndex+bi] = vmesh->getLocalID(GID);
         }
      } else {
         // Several columns in set; We need to gather all data before incrementing indices.
         // We do this by looping twice:
         vmesh::GlobalID nColumnsInSet = 1;
         bool inColumn = true;
         for (vmesh::GlobalID GID = columnMinBlock[Cind]; GID <= columnMaxBlock[Cind]; GID+=stride) {
            vmesh::LocalID LID = vmesh->getLocalID(GID);
            if ( inColumn && (LID == vmesh->invalidLocalID())) {
               // end of column
               inColumn = false;
            } else if (!inColumn && (LID != vmesh->invalidLocalID())) {
               // start of new column
               nColumnsInSet++;
               inColumn = true;
            }
         }
         const uint ColumnSetIndex = atomicAdd(dev_ColumnSetIndex,1);
         uint ColumnIndex = atomicAdd(dev_ColumnIndex,nColumnsInSet);
         uint BlockIndex = atomicAdd(dev_ColumnIndex,columnNBlocks[Cind]);
         columnData->setNumColumns[ColumnSetIndex] = nColumnsInSet; // how many columns in set of columns (length nColumnSets)
         columnData->setColumnOffsets[ColumnSetIndex]= BlockIndex; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
         // Now we know how many columns exist so we can store the data.
         vmesh::GlobalID blocksInThisColumn=0;
         vmesh::GlobalID blocksOffsetForThisColumn=BlockIndex;
         vmesh::GlobalID columnIndexInSet=0;
         inColumn = true;
         for (vmesh::GlobalID GID = columnMinBlock[Cind]; GID <= columnMaxBlock[Cind]; GID+=stride) {
            vmesh::LocalID LID = vmesh->getLocalID(GID);
            if (LID != vmesh->invalidLocalID()) {
               blocksGID[BlockIndex] = GID;
               blocksLID[BlockIndex] = LID;
               BlockIndex++;
               blocksInThisColumn++;
            }
            if ( inColumn && (LID == vmesh->invalidLocalID())) {
               // end of column
               columnData->columnNumBlocks[ColumnIndex+columnIndexInSet] = blocksInThisColumn; // length of column (in blocks, length totalColumns)
               columnData->columnBlockOffsets[ColumnIndex+columnIndexInSet] = blocksOffsetForThisColumn; // indexes where columns start (in blocks, length totalColumns)
               blocksOffsetForThisColumn += blocksInThisColumn;
               blocksInThisColumn=0;
               inColumn = false;
            } else if (!inColumn && (LID != vmesh->invalidLocalID())) {
               // start of new column
               inColumn = true;
            }
         }
         // Also write data for final column in this set
         columnData->columnNumBlocks[ColumnIndex+columnIndexInSet] = blocksInThisColumn; // length of column (in blocks, length totalColumns)
         columnData->columnBlockOffsets[ColumnIndex+columnIndexInSet] = blocksOffsetForThisColumn; // indexes where columns start (in blocks, length totalColumns)
      }
   }
}

/*
   This function returns a sorted list of blocks in a cell.

   The sorted list is sorted according to the location, along the given dimension.
   This version uses triplets internally and also returns the LIDs of the sorted blocks.
*/
void sortBlocklistByDimension( //const spatial_cell::SpatialCell* spatial_cell,
                               const vmesh::VelocityMesh* vmesh,
                               const uint dimension,
                               vmesh::GlobalID* blocksGID,
                               vmesh::LocalID* blocksLID,
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

#ifdef UNDEFINED
   // Call kernel instead of CPU code

   // Reset search vectors
      phiprof::start("Reset column search vectors");
      cudaMemsetAsync(columnNBlocks[cuda_async_queue_id], 0, sizeof(vmesh::LocalID)*cuda_acc_columnContainerSize, stream);
      cudaMemsetAsync(columnMinBlock[cuda_async_queue_id], 0xFF, sizeof(vmesh::LocalID)*cuda_acc_columnContainerSize, stream);
      cudaMemsetAsync(columnMaxBlock[cuda_async_queue_id], 0, sizeof(vmesh::LocalID)*cuda_acc_columnContainerSize, stream);
      cudaStreamSynchronize(stream);
      phiprof::stop("Reset column search vectors");

      phiprof::start("Scan blocks for required columns kernel");
      // Can't probably realistically launch 200*200 blocks...
      const uint refL=0; //vAMR
      const vmesh::LocalID D0 = vmesh->getGridLength(refL)[dimension];
      const vmesh::LocalID D1 = vmesh->getGridLength(refL)[(dimension+1)%3];
      const vmesh::LocalID D2 = vmesh->getGridLength(refL)[(dimension+2)%3];

      // First pass: go over all existing blocks. Accumulate how many blocks are in each column and store the
      // Smallest and largest GIDs in each column.
      uint nCudaBlocks  = (D1*D2/CUDATHREADS) > CUDABLOCKS ? CUDABLOCKS : (D1*D2/CUDATHREADS);
      dim3 grid1(nCudaBlocks,1,1);
      dim3 block1(CUDATHREADS,1,1);
      scan_blocks_for_columns_kernel<<<grid1, block1, 0, stream>>> (
         vmesh,
         dimension,
         columnNBlocks[cuda_async_queue_id],
         columnMinBlock[cuda_async_queue_id],
         columnMaxBlock[cuda_async_queue_id],
         D0,D1,D2
         );
      cudaStreamSynchronize(stream);
      phiprof::stop("Scan blocks for required columns kernel");

      phiprof::start("Sort blocks by dimension kernel");
      // Set up counters for gathering columns
      vmesh::LocalID *dev_ColumnIndex;
      vmesh::LocalID *dev_ColumnSetIndex;
      vmesh::LocalID *dev_ColumnBlockIndex;
      vmesh::LocalID host_ColumnIndex = 0;
      vmesh::LocalID host_ColumnSetIndex = 0;
      vmesh::LocalID host_ColumnBlockIndex = 0;
      HANDLE_ERROR( cudaMallocAsync((void**)&dev_ColumnIndex, sizeof(vmesh::LocalID), stream) );
      HANDLE_ERROR( cudaMallocAsync((void**)&dev_ColumnSetIndex, sizeof(vmesh::LocalID), stream) );
      HANDLE_ERROR( cudaMallocAsync((void**)&dev_ColumnBlockIndex, sizeof(vmesh::LocalID), stream) );
      HANDLE_ERROR( cudaMemcpyAsync(dev_ColumnIndex, &host_ColumnIndex, sizeof(vmesh::LocalID), cudaMemcpyHostToDevice, stream) );
      HANDLE_ERROR( cudaMemcpyAsync(dev_ColumnSetIndex, &host_ColumnSetIndex, sizeof(vmesh::LocalID), cudaMemcpyHostToDevice, stream) );
      HANDLE_ERROR( cudaMemcpyAsync(dev_ColumnBlockIndex, &host_ColumnBlockIndex, sizeof(vmesh::LocalID), cudaMemcpyHostToDevice, stream) );

      /** Now we launch a kernel which goes over every potential column. It skips columns which do not have any blocks.
          Within non-empty columns, it starts from the smallest GID and finishes at the largest GID with a
          dimension-specific stride. It then searches if those blocks exist, and places them at the suitable offset
          within the blocksGID and blocksLID lists.

          Do we actually need "sets" of columns or could we just have columns?
      **/
      nCudaBlocks = (vmesh->getGridLength(refL)[0]) > CUDABLOCKS ? CUDABLOCKS : (vmesh->getGridLength(refL)[0]);
      dim3 grid2(nCudaBlocks,1,1);
      //dim3 block1(CUDATHREADS,1,1);
      dim3 block2(1,1,1); // Currently only supports 1 thread
      sort_blocks_by_dimension_kernel<<<grid2, block2, 0, stream>>> (
         vmesh,
         dimension,
         columnNBlocks[cuda_async_queue_id],
         columnMinBlock[cuda_async_queue_id],
         columnMaxBlock[cuda_async_queue_id],
         D0,D1,D2,
         blocksGID,
         blocksLID,
         columnData,
         dev_ColumnIndex, // make into one 3-element pointer
         dev_ColumnSetIndex,
         dev_ColumnBlockIndex
         );
      cudaStreamSynchronize(stream);
      HANDLE_ERROR( cudaFree(dev_ColumnIndex) );
      HANDLE_ERROR( cudaFree(dev_ColumnSetIndex) );
      HANDLE_ERROR( cudaFree(dev_ColumnBlockIndex) );
      phiprof::stop("Sort blocks by dimension kernel");
      return;
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
}
