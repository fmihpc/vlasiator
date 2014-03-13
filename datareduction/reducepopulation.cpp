#include <iostream>
#include <vector>
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
      //uint64_t index; //Index could be uint32_t is enough
      uint32_t blockId;
      uint16_t vCellId;

      inline void set_data( const SpatialCell * input_cell, const uint32_t index ) {
         cell = input_cell;
         blockId = cell->velocity_block_list[index/64];
         vCellId = index%64;
      }

      inline void set_data( const SpatialCell * input_cell, const uint32_t input_blockId, const uint16_t input_vCellId ) {
         cell = input_cell;
         blockId = input_blockId;
         vCellId = input_vCellId;
      }

      // Compare values
      bool operator<( const Velocity_Cell & other ) const {
         return cell->at(blockId)->data[vCellId] < cell->at(blockId)->data[other.vCellId];
      }

      // Compare equality
      bool operator==(  const Velocity_Cell & other ) const {
         return (blockId == other.blockId) && (vCellId == other.vCellId);
      }

      // Function for getting the avgs value
      inline Realf get_avgs() const {
         return cell->at(blockId)->data[vCellId];
      }
};


// Returns neighbors of the velociy cell
Velocity_Cell * neighbors( const SpatialCell * cell, const Velocity_Cell & vCell ) {
   // Get the velocity cell block id:
   const uint32_t blockId = vCell.blockId;
   // Get the velocity cell id:
   const uint16_t vCellId = vCell.vCellId;
   // Transform to coordinates:
   const uint32_t indices[3] = { 
                        WID * (blockId % SpatialCell::vx_length) + vCellId % block_vx_length,
                        WID * ((blockId / SpatialCell::vx_length) % SpatialCell::vy_length) + (vCellId / block_vx_length) % block_vy_length,
                        WID * (blockId / (SpatialCell::vx_length * SpatialCell::vy_length)) + vCellId / (block_vx_length * block_vy_length)
                               };
   // Reserve space for the neighbors:
   Velocity_Cell * neighbor_vCells = new Velocity_Cell[block_vx_length*block_vz_length*block_vy_length - 1];
   // Get the neighbors:
   // Go through neighbor indices e.g. i, j, k (-1,0,1)
   unsigned int iterator = 0;
   for( unsigned int i = -1; i <= 1; ++i ) for( unsigned int j = -1; j <= 1; ++j ) for( unsigned int k = -1; k <=1; ++k ) {
      // Do not include self as neighbor
      if( i == 0 && j == 0 && k == 0 ) { continue; }
      // Get the neighbor's indices
      const uint32_t neighbor_indices[3] = {
                                           indices[0] + i,
                                           indices[0] + j,
                                           indices[0] + k,
                                           };
      // Transform into blockid and velocity cell id:
      const uint32_t neighbor_blockId = 
                             neighbor_indices[0] / WID
                             + (neighbor_indices[1] / WID) * SpatialCell::vx_length
                             + (neighbor_indices[2] / WID) * SpatialCell::vx_length * SpatialCell::vy_length;
      const uint16_t neighbor_vCellId = 
                             (neighbor_indices[0] % WID)
                           + (neighbor_indices[1] % WID) * block_vx_length
                           + (neighbor_indices[2] % WID) * block_vx_length * block_vy_length;
      // Get index of the neighbor cell:
      
      neighbor_vCells[iterator].set_data( cell, neighbor_blockId, neighbor_vCellId );
      iterator++;
   }
   return neighbor_vCells;
}



//Fast implementation
Real evaluate_speed( const SpatialCell * cell ) {
   // Sort list of avgs values:
   vector<Velocity_Cell> velocityCells;
   // Get the block values
   const vector<Realf,aligned_allocator<Realf,64> > * block_data = &(cell->block_data);
   // Initialize avgs values vector:
   velocityCells.resize( block_data->size() );
   for( unsigned int i = 0; i < block_data->size(); ++i ) {
      // Create a new velocity cell
      Velocity_Cell input_cell;
      // Input the block data
      input_cell.set_data( cell, i);
      // Input the velocity cell into the vector
      velocityCells[i] = input_cell;
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
Real evaluate_speed_parallel( const SpatialCell * cell ) {
   // Sort list of avgs values:
   vector<Velocity_Cell> velocityCells;
   // Get the block values
   const vector<Realf,aligned_allocator<Realf,64> > * block_data = &(cell->block_data);
   // Initialize avgs values vector:
   velocityCells.resize( block_data->size() );
   #pragma omp parallel for
   for( unsigned int i = 0; i < block_data->size(); ++i ) {
      // Create a new velocity cell
      Velocity_Cell input_cell;
      // Input the block data
      input_cell.set_data( cell, i);
      // Input the velocity cell into the vector
      velocityCells[i] = input_cell;
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



