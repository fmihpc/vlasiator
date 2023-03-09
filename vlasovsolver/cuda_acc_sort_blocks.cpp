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

#include "cuda_acc_sort_blocks.hpp"

using namespace std;
using namespace spatial_cell;

// Comparator function for sorting vector of pairs
inline bool tripletcomparator( const std::pair<std::pair<uint, uint>, uint> & l, const std::pair<std::pair<uint, uint>, uint> & r ) {
   return l.first.first < r.first.first;
}

// Kernel
__global__ void sort_blocks_by_dimension_kernel(
   const vmesh::VelocityMesh* vmesh,
   const uint dimension,
   uint* blocksGID,
   uint* blocksLID,
   ColumnOffsets* columnData
   // contains
   // split::SplitVector<uint> columnBlockOffsets; // indexes where columns start (in blocks, length totalColumns)
   // split::SplitVector<uint> columnNumBlocks; // length of column (in blocks, length totalColumns)
   // split::SplitVector<uint> setColumnOffsets; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
   // split::SplitVector<uint> setNumColumns; // how many columns in set of columns (length nColumnSets)
   ) {
   const int index = threadIdx.x;
   const int trial1 = blockIdx.x;
   const int trial2 = blockIdx.y;


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
                               cudaStream_t stream
   ) {

   columnData->columnBlockOffsets.clear();
   columnData->columnNumBlocks.clear();
   columnData->setColumnOffsets.clear();
   columnData->setNumColumns.clear();

#ifdef UNDEFINED
      // Call kernel instead of CPU code
      phiprof::start("Sort blocks by dimension kernel");
      const uint refL=0; // vAMR
      switch (dimension) {
         case 0:
            dim3 grid(vmesh->getGridLength(refL)[1],vmesh->getGridLength(refL)[2],1);
            break;
         case 1:
            dim3 grid(vmesh->getGridLength(refL)[0],vmesh->getGridLength(refL)[2],1);
            break;
         case 2:
            dim3 grid(vmesh->getGridLength(refL)[0],vmesh->getGridLength(refL)[1],1);
            break:
         default:
            printf("Invalid dimension in cuda sortBlockListByDimension!\n");
            abort();
      }
      //dim3 block(WID,WID,WID);
      dim3 block1(1,1,1);
      //const int nCudaBlocks = nBlocks > CUDABLOCKS ? CUDABLOCKS : nBlocks;
      sort_blocks_by_dimension_kernel<<<grid, block1, 0, stream>>> (
         vmesh,
         dimension,
         blocksGID,
         blocksLID,
         columnData
         );
      cudaStreamSynchronize(stream);
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