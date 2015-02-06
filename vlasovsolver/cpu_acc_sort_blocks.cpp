
/*
This file is part of Vlasiator.
Copyright 2013-2015 Finnish Meteorological Institute

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
#warning "unfinished documentation"
void sortBlocklistByDimension( const spatial_cell::SpatialCell* spatial_cell, 
                               const uint dimension,
                               uint* blocks,
                               std::vector<uint> & columnBlockOffsets,
                               std::vector<uint> & columnNumBlocks,
                               std::vector<uint> & setColumnOffsets,
                               std::vector<uint> & setNumColumns) {
   const uint nBlocks = spatial_cell->get_number_of_velocity_blocks(); // Number of blocks
   // Copy block data to vector
   std::vector<std::pair<uint, uint> > block_pairs;
   block_pairs.resize( nBlocks );
   for (vmesh::LocalID i = 0; i < nBlocks; ++i ) {
      const vmesh::GlobalID block = spatial_cell->get_velocity_block_global_id(i);
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
          const uint x_indice = block%SpatialCell::get_velocity_grid_length()[0];
          const uint y_indice = (block/SpatialCell::get_velocity_grid_length()[0])%SpatialCell::SpatialCell::get_velocity_grid_length()[1];
          // Mapping the block id to different coordinate system if dimension is not zero:
          const uint blockId_mapped = block - (x_indice + y_indice*SpatialCell::get_velocity_grid_length()[0]) + y_indice + x_indice * SpatialCell::SpatialCell::get_velocity_grid_length()[1];
          block_pairs[i] = std::make_pair( blockId_mapped, block );
       }
         break;
       case 2: {
          // Do operation: 
          //   block = x + y*x_max + z*y_max*x_max 
          //=> block' = z + y*z_max + x*z_max*y_max
          const uint x_indice = block%SpatialCell::get_velocity_grid_length()[0];
          const uint y_indice = (block/SpatialCell::get_velocity_grid_length()[0])%SpatialCell::SpatialCell::get_velocity_grid_length()[1];
          const uint z_indice =  (block/(SpatialCell::get_velocity_grid_length()[0]*SpatialCell::SpatialCell::get_velocity_grid_length()[1]));
          // Mapping the block id to different coordinate system if dimension is not zero:
          const uint blockId_mapped = z_indice + y_indice * SpatialCell::SpatialCell::get_velocity_grid_length()[2] + x_indice*SpatialCell::SpatialCell::get_velocity_grid_length()[1]*SpatialCell::SpatialCell::get_velocity_grid_length()[2];
          block_pairs[i] = std::make_pair( blockId_mapped, block );
       }
         break;
      }
   }
   // Sort the list:
   std::sort( block_pairs.begin(), block_pairs.end(), paircomparator );
   
   // Put in the sorted blocks, and also compute columnoffsets, and column lengths:
   columnBlockOffsets.push_back(0); //first offset
   setColumnOffsets.push_back(0); //first offset   
   uint prev_column_id, prev_dimension_id;

   for (vmesh::LocalID i = 0; i < nBlocks; ++i ) {
      uint column_id; /* identifies a particlular column*/
      uint dimension_id; /*identifies a particular block in a column (along the dimension)*/
      blocks[i] = block_pairs[i].second; //sorted list
      switch( dimension ) {
         case 0:
            column_id = block_pairs[i].first / SpatialCell::get_velocity_grid_length()[0];
            dimension_id = block_pairs[i].first % SpatialCell::get_velocity_grid_length()[0];
            break;
         case 1:
            column_id = block_pairs[i].first / SpatialCell::SpatialCell::get_velocity_grid_length()[1];
            dimension_id = block_pairs[i].first % SpatialCell::SpatialCell::get_velocity_grid_length()[1];            
            break;
         case 2:
            column_id = block_pairs[i].first / SpatialCell::SpatialCell::get_velocity_grid_length()[2];
            dimension_id = block_pairs[i].first % SpatialCell::SpatialCell::get_velocity_grid_length()[2];            
            break;
      }
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
