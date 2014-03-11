#include <vector>
#include "reducepopulation.h"

using namespace std;

// Note: This is done to save memory
class Velocity_Cell {
   private:
      const vector<Realf,aligned_allocator<Realf,64> > * block_data;
   public:
      Velocity_Cell();
      inline void set_block_data( const vector<Realf,aligned_allocator<Realf,64> > * input_block_data ) {
         block_data = input_block_data;
      }
      unsigned int index;
      // Compare values
      bool operator<( const Velocity_Cell & other ) const {
         return block_data[index] < block_data[other.index];
      }
      // Function for getting the avgs value
      inline Realf get_avgs() const {
         return block_data[index];
      }
};

//Fast implementation
Real evaluate_speed( const SpatialCell * cell ) {
   // Sort list of avgs values:
   vector<Velocity_Cell> velocityCells;
   // Get the block values
   const vector<Realf,aligned_allocator<Realf,64> > * input_block_data = &(cell->block_data);
   // Initialize avgs values vector:
   velocityCells.resize( input_block_data->size() );
   #pragma omp parallel for
   for( unsigned int i = 0; i < input_block_data->size(); ++i ) {
      // Create a new velocity cell
      Velocity_Cell input_cell;
      // Input the block data
      input_cell.set_block_data( input_block_data );
      // Input the velocity cell into the vector
      velocityCells[i] = input_cell;
   }
   
}




















