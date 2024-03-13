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

#include "cpu_acc_sort_blocks.hpp"

using namespace std;
using namespace spatial_cell;

// Comparator function for sorting vector of pairs
inline bool paircomparator( const std::pair<uint, uint> & l, const std::pair<uint, uint> & r ) {
   return l.first < r.first;
}

/*
   This function returns a sorted list of blocks in a cell.

   The sorted list is sorted according to the location, along the given dimension.
   
*/
// TODO unfinished documentation
void sortBlocklistByDimension( //const spatial_cell::SpatialCell* spatial_cell,
                               const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
                               const uint dimension,
                               uint* blocks,
                               std::vector<uint> & columnBlockOffsets,
                               std::vector<uint> & columnNumBlocks,
                               std::vector<uint> & setColumnOffsets,
                               std::vector<uint> & setNumColumns) {
   //const uint nBlocks = spatial_cell->get_number_of_velocity_blocks(); // Number of blocks
   const vmesh::LocalID nBlocks = vmesh.size();

   // Velocity mesh refinement level, has no effect here
   // but is needed in some vmesh::VelocityMesh function calls.
   const uint8_t REFLEVEL = 0;
   
   // Copy block data to vector
   std::vector<std::pair<vmesh::GlobalID,vmesh::GlobalID> > block_pairs;
   block_pairs.resize( nBlocks );
   for (vmesh::LocalID i = 0; i < nBlocks; ++i ) {
      //const vmesh::GlobalID block = spatial_cell->get_velocity_block_global_id(i);
      const vmesh::GlobalID block = vmesh.getGlobalID(i);
      switch( dimension ) {
       case 0: {
          const vmesh::GlobalID blockId_mapped = block; // Mapping the block id to different coordinate system if dimension is not zero:
          block_pairs[i] = std::make_pair( blockId_mapped, block );
       }
         break;
       case 1: {
          // Do operation: 
          //   block = x + y*x_max + z*y_max*x_max 
          //=> block' = block - (x + y*x_max) + y + x*y_max = x + y*x_max + z*y_max*x_max - (x + y*x_max) + y + x*y_max
          //          = y + x*y_max + z*y_max*x_max
          //const uint x_indice = block%SpatialCell::get_velocity_grid_length()[0];
          //const uint y_indice = (block/SpatialCell::get_velocity_grid_length()[0])%SpatialCell::SpatialCell::get_velocity_grid_length()[1];
          const vmesh::LocalID x_index = block % vmesh.getGridLength(REFLEVEL)[0];
          const vmesh::LocalID y_index = (block / vmesh.getGridLength(REFLEVEL)[0]) % vmesh.getGridLength(REFLEVEL)[1];

          // Mapping the block id to different coordinate system if dimension is not zero:
          //const uint blockId_mapped 
          //        = block - (x_indice + y_indice*SpatialCell::get_velocity_grid_length()[0]) 
          //        + y_indice 
          //        + x_indice * SpatialCell::SpatialCell::get_velocity_grid_length()[1];
          const vmesh::GlobalID blockId_mapped 
                  = block - (x_index + y_index*vmesh.getGridLength(REFLEVEL)[0])
                  + y_index 
                  + x_index * vmesh.getGridLength(REFLEVEL)[1];
          block_pairs[i] = std::make_pair( blockId_mapped, block );
       }
         break;
       case 2: {
          // Do operation: 
          //   block = x + y*x_max + z*y_max*x_max 
          //=> block' = z + y*z_max + x*z_max*y_max
          //const uint x_indice = block%SpatialCell::get_velocity_grid_length()[0];
          //const uint y_indice = (block/SpatialCell::get_velocity_grid_length()[0])%SpatialCell::SpatialCell::get_velocity_grid_length()[1];
          //const uint z_indice =  (block/(SpatialCell::get_velocity_grid_length()[0]*SpatialCell::SpatialCell::get_velocity_grid_length()[1]));
          const vmesh::LocalID x_index = block % vmesh.getGridLength(REFLEVEL)[0];
          const vmesh::LocalID y_index = (block / vmesh.getGridLength(REFLEVEL)[0]) % vmesh.getGridLength(REFLEVEL)[1];
          const vmesh::LocalID z_index = (block / (vmesh.getGridLength(REFLEVEL)[0]*vmesh.getGridLength(REFLEVEL)[1]));

          // Mapping the block id to different coordinate system if dimension is not zero:
          //const uint blockId_mapped 
          //  = z_indice 
          //  + y_indice * SpatialCell::SpatialCell::get_velocity_grid_length()[2] 
          //  + x_indice*SpatialCell::SpatialCell::get_velocity_grid_length()[1]*SpatialCell::SpatialCell::get_velocity_grid_length()[2];
          const vmesh::GlobalID blockId_mapped 
            = z_index 
            + y_index*vmesh.getGridLength(REFLEVEL)[2]
            + x_index*vmesh.getGridLength(REFLEVEL)[1]*vmesh.getGridLength(REFLEVEL)[2];
          block_pairs[i] = std::make_pair( blockId_mapped, block );
       }
         break;
      }
   }
   // Sort the list:
   std::sort( block_pairs.begin(), block_pairs.end(), paircomparator );

   // Put in the sorted blocks, and also compute column offsets and lengths:
   columnBlockOffsets.push_back(0); //first offset
   setColumnOffsets.push_back(0); //first offset   
   uint prev_column_id, prev_dimension_id;

   for (vmesh::LocalID i=0; i<nBlocks; ++i) {
       // identifies a particular column
       vmesh::LocalID column_id = block_pairs[i].first / vmesh.getGridLength(REFLEVEL)[dimension];     
       
       // identifies a particular block in a column (along the dimension)
       vmesh::LocalID dimension_id = block_pairs[i].first % vmesh.getGridLength(REFLEVEL)[dimension];
      
       //sorted list
       blocks[i] = block_pairs[i].second;

      if ( i > 0 &&  ( column_id != prev_column_id || dimension_id != (prev_dimension_id + 1) )){
         //encountered new column! For i=0, we already entered the correct offset (0).
         //We also identify it as a new column if there is a break in the column (e.g., gap between two populations)
         /*add offset where the next column will begin*/
         columnBlockOffsets.push_back(i); 
         /*add length of the current column that now ended*/
         columnNumBlocks.push_back(columnBlockOffsets[columnBlockOffsets.size()-1] - columnBlockOffsets[columnBlockOffsets.size()-2]);

         if (column_id != prev_column_id ){
            //encountered new set of columns, add offset to new set starting at present column
            setColumnOffsets.push_back(columnBlockOffsets.size() - 1);
            /*add length of the previous column set that ended*/
            setNumColumns.push_back(setColumnOffsets[setColumnOffsets.size()-1] - setColumnOffsets[setColumnOffsets.size()-2]);
         }
      }      
      prev_column_id = column_id;
      prev_dimension_id = dimension_id;                        
   }
   
   columnNumBlocks.push_back(nBlocks - columnBlockOffsets[columnBlockOffsets.size()-1]);
   setNumColumns.push_back(columnNumBlocks.size() - setColumnOffsets[setColumnOffsets.size()-1]);
}
