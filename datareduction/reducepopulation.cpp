#include <iostream>
#include <unordered_set>
#include <parallel/algorithm>
#include <algorithm>
#include "phiprof.hpp"
#include "reducepopulation.h"



using namespace std;
using namespace spatial_cell;

// Note: This is done to save memory (hopefully)
class Velocity_Cell {
   private:
      //const SpatialCell * cell;
   public:
      //uintVELOCITY_BLOCK_LENGTH_t index; //Index could be uint32_t is enough
      //uint32_t blockId;
      const Velocity_Block * block;
      uint16_t vCellId;

      inline void set_data( const Velocity_Block * input_block, const uint16_t input_vCellId ) {
         block = input_block;
         vCellId = input_vCellId;
      }

      // Compare values
      bool operator<( const Velocity_Cell & other ) const {
         return block->data[vCellId] < other.block->data[other.vCellId];
      }

      // Compare equality
      bool operator==(  const Velocity_Cell & other ) const {
         // Compare memory addresses of blocks and velocity cell id
         return (block == other.block) && (vCellId == other.vCellId);
      }

      // Function for returning hash values for the velocity cell
      inline size_t hash( const size_t starting_point ) const {
         // Return the address of block's data in size_t form plus the vCellId
         return ((size_t)(static_cast<void*>(block->data)) - starting_point)/sizeof(Realf) + (size_t)(vCellId);
      }

      // Function for getting the avgs value
      inline Realf get_avgs() const {
         return block->data[vCellId];
      }
};

// Create neighbors:
//static array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> local_vcell_neighbors;
//static array< vector< pair<uint16_t, vector<uint16_t> > >  , VELOCITY_BLOCK_LENGTH> remote_vcell_neighbors; // Note: This contains both velocity block index and velocity cell index

//A function for creating neighbors for velocity cells.
// Note: local_vcell_neighbors[vCellId] gives the local neighbors of a velocity cell (neighbors that are within the same block)
// Note: remote_vcell_neighbors[vCellId] Gives a vector containing remote neighbors of the velocity cell (neighbors that are outside the block) in vector< pair<uint16_t, vector<uint16_t> > > format. The pair's first index gives the neighbor block index (check spatial_cell.hpp for more information) and the second index gives the local velocity cells withing that neighbor block.
void set_local_and_remote_velocity_cell_neighbors( 
       array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
       array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                                                 ) {
   // Go through every velocity cell in a block
   for( uint i = 0; i < WID; ++i ) for( uint j = 0; j < WID; ++j ) for( uint k = 0; k < WID; ++k ) {
      const uint16_t vCellId = i + j * block_vx_length + k * block_vx_length * block_vy_length;
      // Get the offsets:
      for( int i_offset = -1; i_offset <= 1; ++i_offset ) for( int j_offset = -1; j_offset <= 1; ++j_offset ) for( int k_offset = -1; k_offset <= 1; ++k_offset ) {
         // if i=j=k=0 then we're looking at the velocity cell itself, not neighbor
         if( i_offset == 0 && j_offset == 0 && k_offset == 0 ) { continue; }
         // Get the new indices:
         const int numberOfDirections = 3;
         int neighbor_indices[numberOfDirections] = {
                                                    (int)(i) + i_offset,
                                                    (int)(j) + j_offset,
                                                    (int)(k) + k_offset,
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
               remote_vcell_neighbors[vCellId].reserve( remote_vcell_neighbors[vCellId].size() + 1 );
               remote_vcell_neighbors[vCellId].push_back( blockAndVCell );
            } else {
               // Get the pair
               pair<uint16_t, vector<uint16_t> > & blockAndVCell = remote_vcell_neighbors[vCellId][index];
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



// A function for retrieving the velocity cell neighbors of a given velocity cell. Note: There are alwayus 3*3*3-1=26 neighbors
const static int numberOfVCellNeighbors = 26;
static inline void get_neighbors( 
                vector<Velocity_Cell> & neighbors,
                const Velocity_Cell & vCell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                                                                        ) {
   // Get the local neighbors:
   const uint16_t vCellId = vCell.vCellId;
   for( vector<uint16_t>::const_iterator it = local_vcell_neighbors[vCellId].begin(); it != local_vcell_neighbors[vCellId].end(); ++it ) {
      const uint16_t neighbor_vCellId = *it;
      // Add to the next neighbor:
      Velocity_Cell neighbor;
      // Still in the same block , new vCellId
      neighbor.set_data( vCell.block, neighbor_vCellId );
      neighbors.push_back( neighbor );
   }
   // Get the remote neighbors
   for( vector< pair<uint16_t, vector<uint16_t> > >::const_iterator it = remote_vcell_neighbors[vCellId].begin(); it != remote_vcell_neighbors[vCellId].end(); ++it ) {
      // Get the neighbor block's index
      const uint16_t neighbor_block_index = get<0>(*it);
      // Go through the local cells in the neighbor block:
      for( vector<uint16_t>::const_iterator jt = local_vcell_neighbors[vCellId].begin(); jt != local_vcell_neighbors[vCellId].end(); ++jt ) {
         const uint16_t neighbor_vCellId = *jt;
         // Check for null pointer:
         if( vCell.block->neighbors[neighbor_block_index] != NULL ) {
            Velocity_Cell neighbor;
            neighbor.set_data( vCell.block->neighbors[neighbor_block_index], neighbor_vCellId );
            neighbors.push_back( neighbor );
         }
      }
   }
   return;
}

//TODO: FINISH!
static inline void cluster( 
                  const vector<Velocity_Cell> & velocityCells,
                  const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                  const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
                  SpatialCell * cell
                          ) {
   // Start putting neighbors in a list:
   // Note: This space is already reserved but it is not being used currently:
   cerr << __LINE__ << endl;
   uint32_t * hashtable = reinterpret_cast<uint32_t*>(cell->block_fx.data());
   // Put everything to zero:
   for( unsigned int i = 0; i < velocityCells.size(); ++i ) {
      hashtable[i] = 0;
   }
   cerr << __LINE__ << endl;


   const size_t startingPoint = (size_t)(static_cast<const void*>((cell->block_data).data()));


   vector<Velocity_Cell> neighbors;
   neighbors.reserve( numberOfVCellNeighbors );
   // Do the algorithm:
   unsigned int clusterId = 1; // The largest cluster number. This gets updated as the algorithm progresses
   // Put the velocity cell with the highest value into the hash table:
   const int last_index = velocityCells.size()-2;
   hashtable[velocityCells[last_index].hash( startingPoint )] = clusterId;
   cerr << __LINE__ << endl;
   for( int i = velocityCells.size()-2; i >= 0; --i ) {
      const Velocity_Cell & vCell = velocityCells[i];
      // Get the hashed index of this velocity cell:
      const size_t index = vCell.hash( startingPoint );
      // Get the neighbors:
      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors );
      // Truth value for whether the velocity cell is a part of any cluster
      bool found_cluster = false;
      for( unsigned int j = 0; j < neighbors.size(); ++j ) {
         // Get the hashed neighbor_index:
         const size_t neighbor_index = neighbors[j].hash( startingPoint );
         // Check if the neighbor is a part of some cluster -- returns 0 if not
         const unsigned int neighbor_clusterId = hashtable[neighbor_index];
   cerr << __LINE__ << endl;
         if( neighbor_clusterId != 0 ) {
            if( found_cluster == true ) {
               // Already part of some cluster, so now we can link it to other clusters too:
               //TODO
               // For now, just print a to get the idea of how many clusters there might be
               //cerr << "HIT" << endl;
            } else {
               // Found out that the velocity cell is a part of a cluster -> don't update cluster number (because this is not a new cluster)
               // The velocity cell is a part of a cluster now because its neighbor is a part of a cluster, so set the cluster id equal to its neighbor's
               hashtable[index] = neighbor_clusterId;
               found_cluster = true;
            }
         }
      }
   cerr << __LINE__ << endl;
      // Check if the velocity cell was a part of a cluster, if not then make it a new cluster:
      if( found_cluster == false ) {
         // New cluster, so new cluster id
         ++clusterId;
         // Set the cluster id
         hashtable[index] = clusterId;
      }
   }
   cerr << __LINE__ << endl;
   // Print how many cluster ids:
   cerr << "CLUSTER IDS: " <<  clusterId << endl;
   return;
}


static void test_neighbor(
                SpatialCell * cell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                  ) {
   cerr << "TESTING NEIGHBOR" << endl;
   // Print out velocity cell neighbors:
   for( uint i = 0; i < local_vcell_neighbors.size(); ++i ) {
      // Velocity cell:
      cerr << "Velocity cell: " << i << " - ";
      cerr << "LOCALS: " << endl;
      // Its indices:
      velocity_cell_indices_t indices = cell->get_velocity_cell_indices( i );
      for( uint j = 0; j < 3; ++j ) {
         cerr << indices[j] << " ";
      }
      cerr << endl;
      for( vector<uint16_t>::const_iterator it = local_vcell_neighbors[i].begin(); it != local_vcell_neighbors[i].end(); ++it ) {
         cerr << "Neighbor: " << *it << " - ";
         velocity_cell_indices_t indices2 = cell->get_velocity_cell_indices( *it );
         for( uint j = 0; j < 3; ++j ) {
            cerr << indices2[j] << " ";
         }

         cerr << endl;
      }
      cerr << "REMOTES: " << endl;

      for( vector< pair<uint16_t, vector<uint16_t> > >::const_iterator it = remote_vcell_neighbors[i].begin();
      it != remote_vcell_neighbors[i].end(); ++it ) {
         for( vector<uint16_t>::const_iterator jt = get<1>(*it).begin(); jt != get<1>(*it).end(); ++jt ) {
            cerr << "Neighbor2: " << *jt << " - ";
            velocity_cell_indices_t indices2 = cell->get_velocity_cell_indices( *jt );
            for( uint j = 0; j < 3; ++j ) {
               cerr << indices2[j] << " ";
            }
            cerr << endl;
         }
      }

      
   }

}


//   cerr << "REMOTES: " << endl;
//   // Same for remote neighbors:
//   // Print out velocity cell neighbors:
//   for( uint i = 0; i < local_vcell_neighbors.size(); ++i ) {
//      // Velocity cell:
//      cerr << "Velocity cell: " << i << " - ";
//      // Its indices:
//      velocity_cell_indices_t indices = cell->get_velocity_cell_indices( i );
//      for( uint j = 0; j < 3; ++j ) {
//         cerr << indices[j] << " ";
//      }
//      cerr << endl;
//      for( vector< pair<uint16_t, vector<uint16_t> > >::const_iterator it = remote_vcell_neighbors[i].begin();
//      it != remote_vcell_neighbors[i].end(); ++it ) {
//         for( vector<uint16_t>::const_iterator jt = get<1>(*it).begin(); jt != get<1>(*it).end(); ++jt ) {
//            cerr << "Neighbor2: " << *jt << " - ";
//            velocity_cell_indices_t indices2 = cell->get_velocity_cell_indices( *jt );
//            for( uint j = 0; j < 3; ++j ) {
//               cerr << indices2[j] << " ";
//            }
//            cerr << endl;
//         }
//      }
//   }



//Fast implementation
Real evaluate_speed( 
                SpatialCell * cell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                   ) {
   // Test neighbor functioning:
   //test_neighbor( cell, local_vcell_neighbors,remote_vcell_neighbors );
   //cerr << "TESTED!" << endl;

   // Sort list of avgs values:
   vector<Velocity_Cell> velocityCells;
//   // Get the block values
//   const vector<Realf,aligned_allocator<Realf,VELOCITY_BLOCK_LENGTH> > * block_data = &(cell->block_data);
   // Get the pointer address of block data:
   const size_t startingPoint = (size_t)(static_cast<const void*>((cell->block_data).data()));
   // Initialize avgs values vector:
   velocityCells.resize( cell->number_of_blocks * VELOCITY_BLOCK_LENGTH );
   for( unsigned int i = 0; i < cell->number_of_blocks; ++i ) {
      // Create a new velocity cell
      Velocity_Cell input_cell;
      const uint32_t blockId = cell->velocity_block_list[i];
      for( uint16_t vCellId = 0; vCellId < VELOCITY_BLOCK_LENGTH; ++vCellId ) {
         // Input the block data
         input_cell.set_data( cell->at(blockId), vCellId);
         // Input the velocity cell into the vector
         velocityCells[i * VELOCITY_BLOCK_LENGTH + vCellId] = input_cell;
      }
   }
   // Sort the list:
   sort(velocityCells.begin(), velocityCells.end());
   // Return value:
   Real value_to_return = 0;
   for( unsigned int i = 0; i < velocityCells.size(); ++i ) {
      if( i%2 == 0 ) {
         value_to_return += (Real)(velocityCells[i].get_avgs());
      } else {
         value_to_return -= (Real)(velocityCells[i].get_avgs());
      }
   }


   // Start getting velocity neighbors
   vector<Velocity_Cell> neighbors;
   neighbors.reserve( numberOfVCellNeighbors );
   
   for( int i = velocityCells.size()-1; i >= 0; --i ) {
      const Velocity_Cell & vCell = velocityCells[i];
      // Get the neighbors:
      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors );
      for( vector<Velocity_Cell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
         value_to_return += it->vCellId;
      }
      neighbors.clear();
   }

   cluster( velocityCells, local_vcell_neighbors, remote_vcell_neighbors, cell );

   return value_to_return;
}

//Fast implementation
Real evaluate_speed_parallel(
                SpatialCell * cell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
                            ) {
//   // Sort list of avgs values:
//   vector<Velocity_Cell> velocityCells;
////   // Get the block values
////   const vector<Realf,aligned_allocator<Realf,VELOCITY_BLOCK_LENGTH> > * block_data = &(cell->block_data);
//   // Get the pointer address of block data:
//   const size_t startingPoint = (size_t)(static_cast<const void*>((cell->block_data).data()));
//   // Initialize avgs values vector:
//   velocityCells.resize( cell->number_of_blocks * VELOCITY_BLOCK_LENGTH );
//   #pragma omp parallel for
//   for( unsigned int i = 0; i < cell->number_of_blocks; ++i ) {
//      // Create a new velocity cell
//      Velocity_Cell input_cell;
//      const uint32_t blockId = cell->velocity_block_list[i];
//      for( unsigned int vCellId = 0; vCellId < VELOCITY_BLOCK_LENGTH; ++vCellId ) {
//         // Input the block data
//         input_cell.set_data( cell->at(blockId), vCellId);
//         // Input the velocity cell into the vector
//         velocityCells[i * VELOCITY_BLOCK_LENGTH + vCellId] = input_cell;
//      }
//   }
//   // Sort the list:
//   __gnu_parallel::sort(velocityCells.begin(), velocityCells.end());
//   // Return value:
//   Real value_to_return = 0;
//   #pragma omp parallel
//   {
//      Real tmp_value = 0;
//      #pragma omp for
//      for( unsigned int i = 0; i < velocityCells.size(); ++i ) {
//         if( i%2 == 0 ) {
//            tmp_value += (Real)(velocityCells[i].get_avgs());
//         } else {
//            tmp_value -= (Real)(velocityCells[i].get_avgs());
//         }
//      }
//      #pragma omp atomic
//      value_to_return += tmp_value;
//   }
//
//   // Start getting velocity neighbors
//   #pragma omp parallel for
//   for( int i = velocityCells.size()-1; i >= 0; --i ) {
//      array<Velocity_Cell, numberOfVCellNeighbors> neighbors;
//      const Velocity_Cell & vCell = velocityCells[i];
//      // Get the neighbors:
//      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors );
//      for( unsigned int j = 0; j < neighbors.size(); ++j ) {
//         if( neighbors[j].block ) {
//            value_to_return += neighbors[j].vCellId;
//         }
//      }
//   }
//

   //return value_to_return;
   return 0;
}



