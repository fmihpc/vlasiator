
/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_SORT_BLOCKS_FOR_ACC_H
#define CPU_SORT_BLOCKS_FOR_ACC_H

#include "algorithm"
#include "cmath"
#include "utility"
#include "common.h"
#include "spatial_cell.hpp"




// Comparator function for sorting vector of pairs
inline bool paircomparator( const pair<uint, uint> & l, const pair<uint, uint> & r ) {
   return l.first < r.first;
}

/*
   This function copies the block list from spatial cell to blocks and sorts it

   Note: blocks must be allocated
*/
static void sort_blocklist_by_dimension( const SpatialCell* spatial_cell, 
                                         const uint dimension,
                                         uint* blocks,
                                         std::vector<uint> & block_column_offsets ) {
  const uint nBlocks = spatial_cell->number_of_blocks; // Number of blocks
  // Copy block data to vector
  vector<pair<uint, uint> > block_pairs;
  block_pairs.resize( nBlocks );
  for( uint i = 0; i < nBlocks; ++i ) {
    const uint block = spatial_cell->velocity_block_list[i];
    switch( dimension ) {
    case 0:
      {
      const uint blockId_mapped = block; // Mapping the block id to different coordinate system if dimension is not zero:
      block_pairs[i] = make_pair( blockId_mapped, block );
      break;
      }
    case 1:
      {
      // Do operation: 
      //   block = x + y*x_max + z*y_max*x_max 
      //=> block' = block - (x + y*x_max) + y + x*y_max = x + y*x_max + z*y_max*x_max - (x + y*x_max) + y + x*y_max
      //          = y + x*y_max + z*y_max*x_max
      const uint x_indice = block%SpatialCell::vx_length;
      const uint y_indice = (block/SpatialCell::vx_length)%SpatialCell::vy_length;
      // Mapping the block id to different coordinate system if dimension is not zero:
      const uint blockId_mapped = block - (x_indice + y_indice*SpatialCell::vx_length) + y_indice + x_indice * SpatialCell::vy_length;
      block_pairs[i] = make_pair( blockId_mapped, block );
      break;
      }
    case 2:
      {
      // Do operation: 
      //   block = x + y*x_max + z*y_max*x_max 
      //=> block' = z + y*z_max + x*z_max*y_max
      const uint x_indice = block%SpatialCell::vx_length;
      const uint y_indice = (block/SpatialCell::vx_length)%SpatialCell::vy_length;
      const uint z_indice =  (block/(SpatialCell::vx_length*SpatialCell::vy_length));
      // Mapping the block id to different coordinate system if dimension is not zero:
      const uint blockId_mapped = z_indice + y_indice * SpatialCell::vz_length + x_indice*SpatialCell::vy_length*SpatialCell::vz_length;
      block_pairs[i] = make_pair( blockId_mapped, block );
      break;
      }
    }
  }
  // Sort the list:
  sort( block_pairs.begin(), block_pairs.end(), paircomparator );
  // Put in the sorted blocks, and also compute columnoffsets:
  block_column_offset.push_back(0); //first offset
  uint current_column_id
  for( uint i = 0; i < nBlocks; ++i ) {
     uint column_id;
     blocks[i] = block_pairs[i].second;
     switch( dimension ) {
         case 0:
            column_id = block_pairs[i].first / SpatialCell::vx_length;
            break;
         case 0:
            column_id = block_pairs[i].first / SpatialCell::vy_length;
            break;
         case 0:
            column_id = block_pairs[i].first / SpatialCell::vz_length;
            break;
     }
     if ( i > 0 && column_id != current_column_id ){
        //encountered new column! For i=0, we already entered the correct offset (0)
        block_column_offset.push_back(i); //first offset      
     }
     current_column_id = column_id;
  }
  return;
}


/*
   This function copies the block list from spatial cell to blocks and sorts it

   Note: blocks must be allocated
*/
static void sort_blocklist_by_dimension( const SpatialCell* spatial_cell, 
                                         const uint dimension,
                                         uint* blocks ) {
  const uint nBlocks = spatial_cell->number_of_blocks; // Number of blocks
  // Copy block data to vector
  vector<pair<uint, uint> > block_pairs;
  block_pairs.resize( nBlocks );
  for( uint i = 0; i < nBlocks; ++i ) {
    const uint block = spatial_cell->velocity_block_list[i];
    switch( dimension ) {
    case 0:
      {
      const uint blockId_mapped = block; // Mapping the block id to different coordinate system if dimension is not zero:
      block_pairs[i] = make_pair( blockId_mapped, block );
      break;
      }
    case 1:
      {
      // Do operation: 
      //   block = x + y*x_max + z*y_max*x_max 
      //=> block' = block - (x + y*x_max) + y + x*y_max = x + y*x_max + z*y_max*x_max - (x + y*x_max) + y + x*y_max
      //          = y + x*y_max + z*y_max*x_max
      const uint x_indice = block%SpatialCell::vx_length;
      const uint y_indice = (block/SpatialCell::vx_length)%SpatialCell::vy_length;
      // Mapping the block id to different coordinate system if dimension is not zero:
      const uint blockId_mapped = block - (x_indice + y_indice*SpatialCell::vx_length) + y_indice + x_indice * SpatialCell::vy_length;
      block_pairs[i] = make_pair( blockId_mapped, block );
      break;
      }
    case 2:
      {
      // Do operation: 
      //   block = x + y*x_max + z*y_max*x_max 
      //=> block' = z + y*z_max + x*z_max*y_max
      const uint x_indice = block%SpatialCell::vx_length;
      const uint y_indice = (block/SpatialCell::vx_length)%SpatialCell::vy_length;
      const uint z_indice =  (block/(SpatialCell::vx_length*SpatialCell::vy_length));
      // Mapping the block id to different coordinate system if dimension is not zero:
      const uint blockId_mapped = z_indice + y_indice * SpatialCell::vz_length + x_indice*SpatialCell::vy_length*SpatialCell::vz_length;
      block_pairs[i] = make_pair( blockId_mapped, block );
      break;
      }
    }
  }
  // Sort the list:
  sort( block_pairs.begin(), block_pairs.end(), paircomparator );
  // Put in the sorted blocks:
  for( uint i = 0; i < nBlocks; ++i ) {
     blocks[i] = block_pairs[i].second;
  }
  return;
}
