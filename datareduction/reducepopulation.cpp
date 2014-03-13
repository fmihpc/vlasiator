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
      uint64_t index; //Index could be uint32_t is enough

      inline void set_data( const SpatialCell * input_cell, const uint32_t input_index ) {
         cell = input_cell;
         index = input_index;
      }

      // Returns block id
      inline uint32_t blockId() const {
         return cell->velocity_block_list[index/64];
      }

      //Returns velocity cell id
      inline uint16_t vCellId() const {
         return index%64;
      }

      // Compare values
      bool operator<( const Velocity_Cell & other ) const {
         return cell->block_data[index] < cell->block_data[other.index];
      }

      // Compare equality
      bool operator==(  const Velocity_Cell & other ) const {
         return index == other.index;
      }

      // Function for getting the avgs value
      inline Realf get_avgs() const {
         return cell->block_data[index];
      }
};

// Returns neighbors of the velociy cell
Velocity_Cell * neighbors( const Velocity_Cell & vCell ) {
   // Get the velocity cell block id:
   const uint32_t blockId = vCell.blockId();
   // Get the velocity cell id:
   const uint16_t vCellId = vCell.vCellId();
   // Transform to coordinates:
   
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



