#include <iostream>
#include <vector>
#include <unordered_set>
#include <parallel/algorithm>
#include "reducepopulation.h"

using namespace std;

// Note: This is done to save memory
class Velocity_Cell {
   private:
      const SpatialCell * cell;
   public:
      Velocity_Cell(); // Constructor

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
         return block_data[index] < block_data[other.index];
      }
      // Function for getting the avgs value
      inline Realf get_avgs() const {
         return block_data[index];
      }
};

////Fast implementation for clustering
//Real evaluate_speed_one( const SpatialCell * cell ) {
//   
//}

//Fast implementation
Real evaluate_speed( const SpatialCell * cell ) {
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
   for( unsigned int i = 0; i < block_data->size(); ++i ) {
      if( i%2 == 0 ) {
         value_to_return += (Real)(velocityCells[i].get_avgs());
      } else {
         value_to_return -= (Real)(velocityCells[i].get_avgs());
      }
   }
   return value_to_return;
}




















