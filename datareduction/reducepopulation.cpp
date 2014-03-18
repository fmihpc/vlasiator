#include <iostream>
#include <unordered_set>
#include <parallel/algorithm>
#include <algorithm>
#include "reducepopulation.h"



using namespace std;
using namespace spatial_cell;

// Note: This is done to save memory (hopefully)
class Velocity_Cell {
   private:
      const SpatialCell * cell;
   public:
      //uintVELOCITY_BLOCK_LENGTH_t index; //Index could be uint32_t is enough
      //uint32_t blockId;
      const Velocity_Block * block
      uint16_t vCellId;

      inline void set_data( const SpatialCell * input_cell, const Velocity_Block input_block, const uint16_t input_vCellId ) {
         cell = input_cell;
         block = input_block;
         vCellId = input_vCellId;
      }

      // Compare values
      bool operator<( const Velocity_Cell & other ) const {
         return block->data[vCellId] < other.block->data[other.vCellId];
      }

//      // Compare equality
//      bool operator==(  const Velocity_Cell & other ) const {
//         return (blockId == other.blockId) && (vCellId == other.vCellId);
//      }

      // Function for getting the avgs value
      inline Realf get_avgs() const {
         return block->data[vCellId];
      }
};

// Create neighbors:
//static array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> local_vcell_neighbors;
//static array< vector< pair<uint16_t, vector<uint16_t> > >  , VELOCITY_BLOCK_LENGTH> remote_vcell_neighbors; // Note: This contains both velocity block index and velocity cell index


void set_local_and_remote_velocity_cell_neighbors( 
       array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
       array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                                                 ) {
   // Go through every velocity cell
   for( uint i = 0; i < WID; ++i ) for( uint j = 0; j < WID; ++j ) for( uint k = 0; k < WID; ++k ) {
      const uint16_t vCellId = i + j * block_vx_length + k * block_vx_length * block_vy_length;
      // Get the offsets:
      for( int i_offset = -1; i_offset <= 1; ++i_offset ) for( int j_offset = -1; j_offset <= 1; ++j_offset ) for( int k_offset = -1; k_offset <= 1; ++k_offset ) {
         // if i=j=k=0 then we're looking at the velocity cell itself, not neighbor
         if( i_offset == 0 && j_offset == 0 && k_offset == 0 ) { continue; }
         // Get the new indices:
         const int numberOfDirections = 3;
         const int neighbor_indices[numberOfDirections] = {
                                                          i + i_offset,
                                                          j + j_offset,
                                                          k + k_offset,
                                                          };
         bool isRemoteVCell = false;
         // Check if the indices are out of bounds:
         for( uint index = 0; index < numberOfDirections; ++index ) {
            if( neighbor_indices[index] < 0 || neighbor_indices[index] >= WID ) {
               // Out of bounds -> this is a remove velocity cell
               isRemoteVCell = true;
            }
         }
         if( isRemoteVCell ) {
            // Do stuff for remote_vcell_neighbors
            int neighbor_block_direction[numberOfDirections] = {1,1,1};
            // Go through every direction
            for( uint index = 0; index < numberOfDirections; ++index ) {
               if( neighbor_indices[index] < 0 ) {
                  // Move to neighbor block
                  neighbor_block_direction[index] -= 1;
                  // Move to neighbor indices in the neighbor block
                  neighbor_indices[index] = WID-1;
               } else if( neighbor_indices[index] >= WID ) {
                  // Move to neighbor block
                  neighbor_block_direction[index] += 1;
                  // Move to neighbor indices in the neighbor block
                  neighbor_indices[index] = 0;
               }
            }
            // Now get the neighbor block index: Note the 3 and 9 are so that we can get block indices with % and / operators
            const uint16_t neighbor_index = neighbor_block_direction[0]
                                    + neighbor_block_direction[1] * 3
                                    + neighbor_block_direction[2] * 9;
            const uint16_t neighbor_vCellId = neighbor_indices[0] 
                                            + neighbor_indices[1] * block_vx_length
                                            + neighbor_indices[2] * block_vx_length * block_vy_length;
            // Add the neighbor to remote velocity cell neighbors:
            // First check if the velocity block is already within the vector
            int index = -1;
            int iterator = 0;
            for( vector< pair<uint16_t, vector<uint16_t> > >::iterator it = remote_vcell_neighbors[vCellId].begin();
                 it != remote_vcell_neighbors[vCellId].end();
                 ++it ) {
               // Check for velocity block
               const uint16_t iterated_neighbor_index = get<0>(*it);
               if( iterated_neighbor_index == neighbor_index ) {
                  // Found the neighbor index:
                  index = iterator;
               }
               ++iterator;
            }
            // Check if the velocity block was found
            if( index == -1 ) {
               // Velocity block was not found so add it to the list
               vector<uint16_t> neighbor_vcells;
               neighbor_vcells.reserve(1);
               neighbor_vcells.push_back( neighbor_vCellId );
               pair<uint16_t, vector<uint16_t> > blockAndVCell = make_pair( neighbor_index, neighbor_vcells );
               // Add pair to the remote velocity cells
               remote_vcell_neighbors[vCellId].reserve( remote_vcell_neighbors[vCellId].size() + 1 )
               remote_vcell_neighbors[vCellId].push_back( blockAndVCell );
            } else {
               // Get the pair
               pair<uint16_t, vector<uint16_t> > & blockAndVCell = remote_vcell_neighbors[vCellId][index]
               // Add velocity cell:
               vector<uint16_t> & neighbor_vcells = get<1>( blockAndVCell );
               neighbor_vcells.reserve( neighbor_vcells.size() + 1 );
               neighbor_vcells.push_back( neighbor_vCellId );
            }
         } else {
            // This is not a remote velocity cell (meaning this velocity cell neighbor is within the same velocity block as the vCellId)
            // Add to local vcell neighbors:
            const uint16_t neighbor_vCellId = neighbor_indices[0] + neighbor_indices[1] * block_vx_length + neighbor_indices[2] * block_vx_length * block_vy_length;
            // Reserve space
            local_vcell_neighbors[vCellId].reserve( local_vcell_neighbors[vCellId].size() + 1 );
            local_vcell_neighbors[vCellId].push_back( neighbor_vCellId );
         }
      }
   }
}
//indices[0] + indices[1] * block_vx_length + indices[2] * block_vx_length * block_vy_length;
//for(int offset_vx=-1;offset_vx<=1;offset_vx++)
//   for(int offset_vy=-1;offset_vy<=1;offset_vy++)
//      for(int offset_vz=-1;offset_vz<=1;offset_vz++){
//         // location of neighbor pointer in this block's neighbor list
//         int neighbor_index=(1+offset_vx)+(1+offset_vy)*3+(1+offset_vz)*9;
//
//         if (neighbor_block == error_velocity_block) {
//            block_ptr->neighbors[neighbor_index] = NULL;
//         }


// Returns neighbors of the velociy cell
//Velocity_Cell * neighbors( const SpatialCell * cell, const Velocity_Cell & vCell ) {
   // Get the velocity cell neighbors:
   
//}

//// Returns neighbors of the velociy cell
//Velocity_Cell * neighbors( const SpatialCell * cell, const Velocity_Cell & vCell ) {
//   // Get the velocity cell block id:
//   const uint32_t blockId = vCell.blockId;
//   // Get the velocity cell id:
//   const uint16_t vCellId = vCell.vCellId;
//   // Transform to coordinates:
//   const uint32_t indices[3] = { 
//                        WID * (blockId % SpatialCell::vx_length) + vCellId % block_vx_length,
//                        WID * ((blockId / SpatialCell::vx_length) % SpatialCell::vy_length) + (vCellId / block_vx_length) % block_vy_length,
//                        WID * (blockId / (SpatialCell::vx_length * SpatialCell::vy_length)) + vCellId / (block_vx_length * block_vy_length)
//                               };
//   // Reserve space for the neighbors:
//   Velocity_Cell * neighbor_vCells = new Velocity_Cell[block_vx_length*block_vz_length*block_vy_length - 1];
//   // Get the neighbors:
//   // Go through neighbor indices e.g. i, j, k (-1,0,1)
//   int iterator = 0;
//   for( int i = -1; i <= 1; ++i ) for( unsigned int j = -1; j <= 1; ++j ) for( unsigned int k = -1; k <=1; ++k ) {
//      // Do not include self as neighbor
//      if( i == 0 && j == 0 && k == 0 ) { continue; }
//      // Get the neighbor's indices
//      const uint32_t neighbor_indices[3] = {
//                                           indices[0] + i,
//                                           indices[0] + j,
//                                           indices[0] + k,
//                                           };
//      // Transform into blockid and velocity cell id:
//      const uint32_t neighbor_blockId = 
//                             neighbor_indices[0] / WID
//                             + (neighbor_indices[1] / WID) * SpatialCell::vx_length
//                             + (neighbor_indices[2] / WID) * SpatialCell::vx_length * SpatialCell::vy_length;
//      const uint16_t neighbor_vCellId = 
//                             (neighbor_indices[0] % WID)
//                           + (neighbor_indices[1] % WID) * block_vx_length
//                           + (neighbor_indices[2] % WID) * block_vx_length * block_vy_length;
//      // Get index of the neighbor cell:
//      
//      neighbor_vCells[iterator].set_data( cell, neighbor_blockId, neighbor_vCellId );
//      iterator++;
//   }
//   return neighbor_vCells;
//}
//


//Fast implementation
Real evaluate_speed( 
                const SpatialCell * cell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                   ) {
   // Sort list of avgs values:
   vector<Velocity_Cell> velocityCells;
   // Get the block values
   const vector<Realf,aligned_allocator<Realf,VELOCITY_BLOCK_LENGTH> > * block_data = &(cell->block_data);
   // Initialize avgs values vector:
   velocityCells.resize( cell->number_of_blocks * VELOCITY_BLOCK_LENGTH );
   for( unsigned int i = 0; i < cell->number_of_blocks; ++i ) {
      // Create a new velocity cell
      Velocity_Cell input_cell;
      const uint32_t blockId = cell->velocity_block_list[i];
      for( unsigned int vCellId = 0; vCellId < VELOCITY_BLOCK_LENGTH; ++vCellId ) {
         // Input the block data
         input_cell.set_data( cell, cell->at(blockId), vCellId);
         // Input the velocity cell into the vector
         velocityCells[i * VELOCITY_BLOCK_LENGTH + vCellId] = input_cell;
      }
   }
   // Sort the list:
   sort(velocityCells.begin(), velocityCells.end());
   // Return value:
   Real value_to_return = 0;
   for( unsigned int i = 0; i < block_data->size(); ++i ) {
      if( i%2 == 0 ) {
         value_to_return += (Real)(velocityCells[i].get_avgs());
      } else {
         value_to_return -= (Real)(velocityCells[i].get_avgs());
      }
   }
   return value_to_return;
}

//Fast implementation
Real evaluate_speed_parallel(
                const SpatialCell * cell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                            ) {
   // Sort list of avgs values:
   vector<Velocity_Cell> velocityCells;
   // Get the block values
   const vector<Realf,aligned_allocator<Realf,VELOCITY_BLOCK_LENGTH> > * block_data = &(cell->block_data);
   // Initialize avgs values vector:
   velocityCells.resize( cell->number_of_blocks * VELOCITY_BLOCK_LENGTH );
   #pragma omp parallel for
   for( unsigned int i = 0; i < cell->number_of_blocks; ++i ) {
      // Create a new velocity cell
      Velocity_Cell input_cell;
      const uint32_t blockId = cell->velocity_block_list[i];
      for( unsigned int vCellId = 0; vCellId < VELOCITY_BLOCK_LENGTH; ++vCellId ) {
         // Input the block data
         input_cell.set_data( cell, cell->at(blockId), vCellId);
         // Input the velocity cell into the vector
         velocityCells[i * VELOCITY_BLOCK_LENGTH + vCellId] = input_cell;
      }
   }
   // Sort the list:
   __gnu_parallel::sort(velocityCells.begin(), velocityCells.end());
   // Return value:
   Real value_to_return = 0;
   #pragma omp parallel
   {
      Real tmp_value = 0;
      #pragma omp for
      for( unsigned int i = 0; i < block_data->size(); ++i ) {
         if( i%2 == 0 ) {
            tmp_value += (Real)(velocityCells[i].get_avgs());
         } else {
            tmp_value -= (Real)(velocityCells[i].get_avgs());
         }
      }
      #pragma omp atomic
      value_to_return += tmp_value;
   }
   return value_to_return;
}



