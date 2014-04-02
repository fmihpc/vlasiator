#include <iostream>
#include <utility>
#include <unordered_set>
#include <parallel/algorithm>
#include <algorithm>
#include <deque>
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
         return ((size_t)(block->data) - starting_point)/sizeof(Realf) + (size_t)(vCellId);
      }

      // Function for getting the avgs value
      inline Realf get_avgs() const {
         return block->data[vCellId];
      }
};


namespace std {
   template <> struct hash<std::pair<uint32_t, uint32_t>> {
      inline size_t operator()(const std::pair<uint32_t, uint32_t> &v) const {
          std::hash<uint32_t> uint32_t_hasher;
          return uint32_t_hasher(v.first) ^ uint32_t_hasher(v.second);
      }
   };
}



namespace std {
   template <>
   struct hash<Velocity_Cell> {
      inline size_t operator() (const Velocity_Cell & vcell) const {
         std::hash<size_t> size_t_hasher;
         return size_t_hasher( (size_t)(vcell.block)/sizeof(Realf)+vcell.vCellId );
      }
   };
}

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
                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
                const SpatialCell * cell
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
         // Make sure the neighbor is valid
         if( vCell.block->neighbors[neighbor_block_index] && vCell.block->neighbors[neighbor_block_index]->data != cell->null_block_data.data() ) {
            Velocity_Cell neighbor;
            neighbor.set_data( vCell.block->neighbors[neighbor_block_index], neighbor_vCellId );
            neighbors.push_back( neighbor );
         }
      }
   }
   return;
}


//// A function for comparing clusters
//static inline bool same_cluster( const uint32_t id, const uint32_t neighbor_id, uint32_t *** clusterIds ) {
//   // Make sure that the pointers are NOT null:
//   phiprof_assert( clusterIds[id] && clusterIds[neighbor_id] );
//   // Make sure that the pointers are NOT null:
//   phiprof_assert( *clusterIds[id] && *clusterIds[neighbor_id] );
//   return (*clusterIds[id] == *clusterIds[neighbor_id]);
//}
//
//// A function for checking if cluster exists
//static inline bool cluster_exists( const uint32_t id, uint32_t *** clusterIds ) {
//   return clusterIds[id];
//}
//
//// A function for merging clusters
//static inline void merge_clusters( const uint32_t id, const uint32_t neighbor_id, uint32_t *** clusterIds ) {
//   // Make sure that the pointers are NOT null:
//   phiprof_assert( clusterIds[id] && clusterIds[neighbor_id] );
//   uint32_t ** cluster_one = clusterIds[id];
//   uint32_t ** cluster_two = clusterIds[neighbor_id];
//   const uint32_t numberOfCellsInCluster_one = **cluster_one;
//   const uint32_t numberOfCellsInCluser_two = **cluster_two;
//
//   // Increase the value of cluster two:
//   **cluster_two += numberOfCellsInCluster_one;
//   // Delete the value in cluster one:
//   delete *cluster_one;
//   // Delete the cluster two as well
//   // Set the cluster one to point at cluster two:
//   *clusterIds[id] = *cluster_two;
//   
//}
//
//// A function for creating a new cluster
//static inline void create_new_cluster( const uint32_t id, uint32_t *** clusterIds ) {
//   phiprof_assert( clusterIds );
//   // Reserve space for the cluster
//   clusterIds[id] = new uint32_t*;
//   uint32_t * new_cluster = new uint32_t;
//   // Input number of velocity cells in this cluster:
//   *new_cluster = 1;
//   // Set the velocity cell to point at the new cluster
//   *clusterIds[id] = new_cluster;
//}
//
//// A function for adding a velocity cell to an existing cluster
//static inline void add_to_neighbor_cluster( const uint32_t id, const uint32_t neighbor_id, uint32_t *** clusterIds ) {
//   // Make sure the pointer at id is null
//   phiprof_assert( !cluster_exists(id, clusterIds) );
//   // Note: We are assuming that clusterIds[id] is not pointing to anything yet
//   clusterIds[id] = clusterIds[neighbor_id];
//   // Now that there's an extra cluster, increase the cluster number:
//   **clusterIds[id] += 1;
//}
//
//
//// Delete cluster id
//static inline void delete_cluster( const uint32_t id, uint32_t *** clusterIds ) {
//   // Check for null pointers:
//   if( clusterIds ) {
//      if( clusterIds[id] ) {
//         if( *clusterIds[id] ) {
//            delete *clusterIds[id];
//            delete[] clusterIds[id];
//         } else {
//            delete[] clusterIds[id];
//         }
//      }
//   }
//}
//
//
//
//
////TODO: FINISH!
//static inline void cluster( 
//                                const vector<Velocity_Cell> & velocityCells,
//                                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
//                                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
//                                SpatialCell * cell
//                          ) {
//   // Reserve a table for clusters:
//   uint32_t *** clusterIds = new uint32_t**[velocityCells.size()];
//
//   // Set the first velocity cell to first cluster:
//   uint32_t last_vCell = velocityCells.size()-1;
//   create_new_cluster( last_vCell, clusterIds);
//
//
//   // Start getting velocity neighbors
//   vector<Velocity_Cell> neighbors;
//   neighbors.reserve( numberOfVCellNeighbors );
//
//   const size_t startingPoint = (size_t)(cell->block_data.data());
//
//   uint32_t merges = 0;
//
//   for( int i = velocityCells.size()-1; i >= 0; --i ) {
//      const Velocity_Cell & vCell = velocityCells[i];
//      // Get the neighbors:
//      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
//      // Get the id of the velocity cell
//      const uint32_t id = vCell.hash( startingPoint );
//      if( vCell.hash( startingPoint ) < 0 ) {
//         cerr << "SOMETHING WRONG " << vCell.hash( startingPoint ) << endl;
//      } else if( vCell.hash( startingPoint ) >= velocityCells.size() ) {
//         cerr << "SOMETHING WRONG2 " << vCell.hash( startingPoint ) << " "  << velocityCells.size() << endl;
//      }
//      for( vector<Velocity_Cell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
//         // Get the id of the neighbor:
//         const uint32_t neighbor_id = it->hash( startingPoint );
//
//         if( neighbor_id < 0 ) {
//            cerr << "SOMETHING WRONG " << neighbor_id << " " << it->vCellId << endl;
//         } else if( neighbor_id >= velocityCells.size() ) {
//            cerr << "SOMETHING WRONG2 " << neighbor_id << " "  << velocityCells.size() << " " << it->vCellId << endl;
//         }
//
//         // Set the id to equal the same as neighbors' if the neighbor is part of a cluster
//         if( cluster_exists( neighbor_id, clusterIds ) ) {
//            // If the cluster id has not been set yet, set it now:
//            if( !cluster_exists( id, clusterIds ) ) {
//               // Cluster id has not been set yet
//               add_to_neighbor_cluster( id, neighbor_id, clusterIds );
//            } else if( clusterIds[id] != clusterIds[neighbor_id] ) {
//               // id is a part of a cluster already, merge the clusterIds:
//               ++merges;
//               merge_clusters( id, neighbor_id, clusterIds );
//            }
//         }
//      }
//      // If the cell does not have any neighbors that are part of a cluster then this is a new cluster:
//      if( !cluster_exists( id, clusterIds ) ) {
//         create_new_cluster( id, clusterIds );
//      }
//      neighbors.clear();
//   }
//   cerr << cell->block_data.size() << " " << velocityCells.size() << " " << cell->block_fx.capacity() << " " << cell->block_fx.size() << endl;
//   if( cell->block_fx.size() < cell->block_data.size() ) { cerr << "WARNING" << endl; }
//   //Realf * blockdata = cell->block_data.data();
//   for( uint i = 0; i < velocityCells.size(); ++i ) {
//      phiprof_assert( cluster_exists( id, clusterIds ) );
//      const Realf value = **clusterIds[i]+0.1;
//      cell->block_data.at(i) = value;
//   }
//   // Print out the number of clusterIds:
//   cerr << "Merges: " << merges << endl;
//   cerr << "Sizeof: " << sizeof(Realf) << endl;
//
//   // Delete clusters:
//   for( uint i = 0; i < velocityCells.size(); ++i ) {
//      delete_cluster( i, clusterIds );
//   }
//   delete[] clusterIds;
//   return;
//}

//TODO: FINISH!
//static inline void cluster( 
//                                const vector<Velocity_Cell> & velocityCells,
//                                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
//                                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
//                                SpatialCell * cell
//                          ) {
//   // Reserve a table for clusters:
//   uint32_t * clusterIds = new uint32_t[velocityCells.size()];
//   
//   // Initialize to be part of no clusters:
//   const uint32_t noCluster = 1;
//   for( uint i = 0; i < velocityCells.size(); ++i ) {
//      clusterIds[i] = noCluster;
//   }
//
//   // Id for separating clusterIds
//   uint32_t clusterId=2;
//
//   // Set the first velocity cell to cluster one
//   uint32_t last_vCell = velocityCells.size()-1;
//   clusterIds[last_vCell] = clusterId;
//
//   // Start getting velocity neighbors
//   vector<Velocity_Cell> neighbors;
//   neighbors.reserve( numberOfVCellNeighbors );
//
////   
//
//   const size_t startingPoint = (size_t)(cell->block_data.data());
//
//   uint32_t merges = 0;
//
//   for( int i = velocityCells.size()-1; i >= 0; --i ) {
//      const Velocity_Cell & vCell = velocityCells[i];
//      // Get the neighbors:
//      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
//      // Get the id of the velocity cell
//      const uint32_t id = vCell.hash( startingPoint );
//      if( vCell.hash( startingPoint ) < 0 ) {
//         cerr << "SOMETHING WRONG " << vCell.hash( startingPoint ) << endl;
//      } else if( vCell.hash( startingPoint ) >= velocityCells.size() ) {
//         cerr << "SOMETHING WRONG2 " << vCell.hash( startingPoint ) << " "  << velocityCells.size() << endl;
//      }
//      for( vector<Velocity_Cell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
//         // Get the id of the neighbor:
//         const uint32_t neighbor_id = it->hash( startingPoint );
//
//         if( neighbor_id < 0 ) {
//            cerr << "SOMETHING WRONG " << neighbor_id << " " << it->vCellId << endl;
//         } else if( neighbor_id >= velocityCells.size() ) {
//            cerr << "SOMETHING WRONG2 " << neighbor_id << " "  << velocityCells.size() << " " << it->vCellId << endl;
//         }
//
//         // Set the id to equal the same as neighbors' if the neighbor is part of a cluster
//         if( clusterIds[neighbor_id] != noCluster ) {
//            // If the cluster id has not been set yet, set it now:
//            if( clusterIds[id] == noCluster ) {
//               // Cluster id has not been set yet
//               clusterIds[id] = clusterIds[neighbor_id];
//            } else if( clusterIds[id] != clusterIds[neighbor_id] ) {
//               // id is a part of a cluster already, merge the clusterIds:
//               ++merges;
//               //merge_clusterIds( id, neighbor_id,  );
//            }
//         }
//      }
//      // If the cell does not have any neighbors that are part of a cluster then this is a new cluster:
//      if( clusterIds[id] == noCluster ) {
//         ++clusterId;
//         clusterIds[id] = clusterId;
//      }
//      neighbors.clear();
//   }
//   cerr << cell->block_data.size() << " " << velocityCells.size() << " " << cell->block_fx.capacity() << " " << cell->block_fx.size() << endl;
//   if( cell->block_fx.size() < cell->block_data.size() ) { cerr << "WARNING" << endl; }
//   //Realf * blockdata = cell->block_data.data();
//   for( uint i = 0; i < velocityCells.size(); ++i ) {
//      const Realf value = clusterIds[i]+0.1;
//      cell->block_data.at(i) = value;
//   }
//   // Print out the number of clusterIds:
//   cerr << "Clusters: " << clusterId << endl;
//   cerr << "Merges: " << merges << endl;
//   cerr << "Sizeof: " << sizeof(Realf) << endl;
//
////   
//
//   delete[] clusterIds;
//   return;
//}


static inline void cluster_simple( 
                                vector<Velocity_Cell> & velocityCells,
                                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
                                SpatialCell * cell
                          ) {

   if( velocityCells.size() == 0 ) { return; }

   // Info on which are the current neighbors
   deque<Velocity_Cell> all_neighbors;
   //all_neighbors.reserve( 0.1 * velocityCells.size() );

   // Vector for holding already processed velocity cells:
   const uint16_t noCluster = 1;
   uint16_t * all_vCells = new uint16_t[velocityCells.size()];
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      // By default velocity cells are not a part of the cluster
      all_vCells[i] = noCluster;
   }

   // This is used for hashing velocity cells
   const size_t startingPoint = (size_t)(cell->block_data.data());

   
   uint16_t clusterId = noCluster + 1; // The first cluster id

   // Put the first velocity cell as processed
   {
      // Get the last velocity cell:
      const Velocity_Cell & vCell = velocityCells[velocityCells.size() - 1];
      const uint32_t id = vCell.hash( startingPoint );
      phiprof_assert( id < velocityCells.size() );
      all_vCells[ id ] = clusterId;
      // Add the first velocity cell neighbor
      all_neighbors.push_back( vCell );
   }

   const Realf threshold = -100;

   // Iterator that will be used for iterating velocity cells from largest value to smallest value
   vector<Velocity_Cell>::const_reverse_iterator rit = velocityCells.crbegin();
   ++rit; // Skip the first cell because we already put it in the velocity cell


   // Start getting velocity neighbors
   vector<Velocity_Cell> tmp_neighbors;
   tmp_neighbors.reserve( numberOfVCellNeighbors );


   // Start inputting velocity cells:
   for( uint i = 0; i < velocityCells.size(); ++i ) {

      // Get the velocity cell:
      const Velocity_Cell & vCell = all_neighbors.front();

      // Get neighbors
      get_neighbors( tmp_neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
      for( uint j = 0; j < tmp_neighbors.size(); ++j ) {
         const uint32_t id = tmp_neighbors[j].hash( startingPoint );

         phiprof_assert( id >= 0 && id < velocityCells.size() );

         // If the velocity cell is not a part of a cluster, add it
         // NOTE: This is probably not needed
         if( all_vCells[ id ] == noCluster ) {
            // Make sure the threshold holds:
            const Realf avgs_value = tmp_neighbors[j].get_avgs();
            if( avgs_value >= threshold ) {
               // The velocity cell has a bigger value than the threshold so add it
               all_vCells[ id ] = clusterId;
               all_neighbors.push_back( tmp_neighbors[j] );
            }
         }
      }

      tmp_neighbors.clear();

      // The front-most velocity cell is no longer needed
      phiprof_assert( all_neighbors.size() > 0 );
      all_neighbors.pop_front();

      // Check if any neighbors were found that passed the threshold qualification:
      if( all_neighbors.size() == 0 ) {
         // None were found, so this cluster has been processed:
         ++clusterId;
         // Find the next velocity cell:
         for( ; rit != velocityCells.rend(); ++rit ) {
            const uint32_t id = rit->hash( startingPoint );
            if( all_vCells[ id ] == noCluster ) {
               const Realf avgs_value = rit->get_avgs();
               if( avgs_value >= threshold ) {
                  all_vCells[ id ] = clusterId;
                  all_neighbors.push_back( *rit );
                  break;
               } else {
                  rit = velocityCells.rend();
               }
            }
         }
         // If no velocity cell was found then break the algorithm:
         if( all_neighbors.size() == 0 ) { break; }
      }
   }

   phiprof_assert( cell->block_fx.size() >= velocityCells.size() );

   for( uint i = 0; i < velocityCells.size(); ++i ) {
      cell->block_fx.at(i) = all_vCells[i];
   }

   // Print out the number of clusterIds:
//   cerr << "Merges: " << merges << endl;
//   cerr << "Sizeof: " << sizeof(Realf) << endl;

//   

   //delete[] clusterIds;

   return;
}





static inline void cluster( 
                                const vector<Velocity_Cell> & velocityCells,
                                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                                const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
                                SpatialCell * cell
                          ) {

   if( velocityCells.size() == 0 ) { return; }
   // Reserve a table for clusters:
   uint32_t * clusterIds = new uint32_t[velocityCells.size()];
   
   // Initialize to be part of no clusters:
   const uint32_t noCluster = 1;
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      clusterIds[i] = noCluster;
   }


   // Id for separating clusterIds
   uint32_t clusterId=2;

   // Set the first velocity cell to cluster one
   uint32_t last_vCell = velocityCells.size()-1;
   clusterIds[last_vCell] = clusterId;


   // Start getting velocity neighbors
   vector<Velocity_Cell> neighbors;
   neighbors.reserve( numberOfVCellNeighbors );

   const size_t startingPoint = (size_t)(cell->block_data.data());

   uint32_t merges = 0;


   for( int i = velocityCells.size()-1; i >= 0; --i ) {
      const Velocity_Cell & vCell = velocityCells[i];
      // Get the neighbors:
      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
      // Get the id of the velocity cell
      const uint32_t id = vCell.hash( startingPoint );

      phiprof_assert( id >= 0 && id < velocityCells.size() );

      for( vector<Velocity_Cell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
         // Get the id of the neighbor:
         const uint32_t neighbor_id = it->hash( startingPoint );

         phiprof_assert( neighbor_id >= 0 && neighbor_id < velocityCells.size() );

         // Set the id to equal the same as neighbors' if the neighbor is part of a cluster
         if( clusterIds[neighbor_id] != noCluster ) {
            // If the cluster id has not been set yet, set it now:
            if( clusterIds[id] == noCluster ) {
               // Cluster id has not been set yet
               clusterIds[id] = clusterIds[neighbor_id];
            } else if( clusterIds[id] != clusterIds[neighbor_id] ) {
               // Merge the clusters:
               ++merges;
            }
         }
      }
      // If the cell does not have any neighbors that are part of a cluster then this is a new cluster:
      if( clusterIds[id] == noCluster ) {
         ++clusterId;
         clusterIds[id] = clusterId;
      }
      neighbors.clear();
   }

   phiprof_assert( cell->block_fx.size() >= velocityCells.size() );
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      const Realf value = clusterIds[i];
      cell->block_fx.at(i) = value;
   }

   // Print out the number of clusterIds:
   cerr << "Clusters: " << clusterId << endl;
   cerr << "Merges: " << merges << endl;
//   cerr << "Sizeof: " << sizeof(Realf) << endl;

//   

   delete[] clusterIds;

   return;
}




////TODO: FINISH!
//static inline void cluster_backup( 
//                  const vector<Velocity_Cell> & velocityCells,
//                  const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
//                  const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
//                  SpatialCell * cell
//                          ) {
//   // Start putting neighbors in a list:
//   // Note: This space is already reserved but it is not being used currently:
//   
//   //uint32_t * hashtable = reinterpret_cast<uint32_t*>(cell->block_fx.data());
//   uint32_t * hashtable = new uint32_t(velocityCells.size());
//   // Put everything to zero:
//   for( unsigned int i = 0; i < velocityCells.size(); ++i ) {
//      hashtable[i] = 0;
//   }
//   
//
//
//   const size_t startingPoint = (size_t)(cell->block_data.data());
//
//
//   vector<Velocity_Cell> neighbors;
//   neighbors.reserve( numberOfVCellNeighbors );
//   // Do the algorithm:
//   unsigned int clusterId = 1; // The largest cluster number. This gets updated as the algorithm progresses
//   // Put the velocity cell with the highest value into the hash table:
//   const int last_index = velocityCells.size()-2;
//   hashtable[velocityCells[last_index].hash( startingPoint )] = clusterId;
//   
//   for( int i = velocityCells.size()-2; i >= 0; --i ) {
//      const Velocity_Cell & vCell = velocityCells[i];
//      // Get the hashed index of this velocity cell:
//      const size_t index = vCell.hash( startingPoint );
//      // Get the neighbors:
//      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
//      // Truth value for whether the velocity cell is a part of any cluster
//      bool found_cluster = false;
//      for( unsigned int j = 0; j < neighbors.size(); ++j ) {
//         // Get the hashed neighbor_index:
//         const size_t neighbor_index = neighbors[j].hash( startingPoint );
//         // Check if the neighbor is a part of some cluster -- returns 0 if not
//         const unsigned int neighbor_clusterId = hashtable[neighbor_index];
//         if( neighbor_clusterId != 0 ) {
//            if( found_cluster == true ) {
//               // Already part of some cluster, so now we can link it to other clusters too:
//               //TODO
//               // For now, just print a to get the idea of how many clusters there might be
//               //cerr << "HIT" << endl;
//            } else {
//               // Found out that the velocity cell is a part of a cluster -> don't update cluster number (because this is not a new cluster)
//               // The velocity cell is a part of a cluster now because its neighbor is a part of a cluster, so set the cluster id equal to its neighbor's
//               hashtable[index] = neighbor_clusterId;
//               found_cluster = true;
//            }
//         }
//      }
//      // Check if the velocity cell was a part of a cluster, if not then make it a new cluster:
//      if( found_cluster == false ) {
//         // New cluster, so new cluster id
//         ++clusterId;
//         // Set the cluster id
//         hashtable[index] = clusterId;
//      }
//   }
//   
//   // Print how many cluster ids:
//   cerr << "CLUSTER IDS: " <<  clusterId << endl;
//   return;
//}


//static inline void cluster_fast( 
//                  const vector<Velocity_Cell> & velocityCells,
//                  const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
//                  const array< vector< pair<uint16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
//                  SpatialCell * cell
//                               ) {
//   // Start putting neighbors in a list:
//   // Note: This space is already reserved but it is not being used currently:
//   
//   unordered_set<Velocity_Cell, Hash> processed_vcells;
//   processed_vcells.reserve( velocityCells.size() );
//   
//
//
//
//   vector<Velocity_Cell> neighbors;
//   neighbors.reserve( numberOfVCellNeighbors );
//   // Do the algorithm:
//   unsigned int clusterId = 1; // The largest cluster number. This gets updated as the algorithm progresses
//   // Put the velocity cell with the highest value into the hash table:
//   const int last_index = velocityCells.size()-2;
//   hashtable[velocityCells[last_index].hash( startingPoint )] = clusterId;
//   
//   for( int i = velocityCells.size()-2; i >= 0; --i ) {
//      const Velocity_Cell & vCell = velocityCells[i];
//      // Get the hashed index of this velocity cell:
//      const size_t index = vCell.hash( startingPoint );
//      // Get the neighbors:
//      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors );
//      // Truth value for whether the velocity cell is a part of any cluster
//      bool found_cluster = false;
//      for( unsigned int j = 0; j < neighbors.size(); ++j ) {
//         // Get the hashed neighbor_index:
//         const size_t neighbor_index = neighbors[j].hash( startingPoint );
//         // Check if the neighbor is a part of some cluster -- returns 0 if not
//         const unsigned int neighbor_clusterId = hashtable[neighbor_index];
//   
//         if( neighbor_clusterId != 0 ) {
//            if( found_cluster == true ) {
//               // Already part of some cluster, so now we can link it to other clusters too:
//               //TODO
//               // For now, just print a to get the idea of how many clusters there might be
//               //cerr << "HIT" << endl;
//            } else {
//               // Found out that the velocity cell is a part of a cluster -> don't update cluster number (because this is not a new cluster)
//               // The velocity cell is a part of a cluster now because its neighbor is a part of a cluster, so set the cluster id equal to its neighbor's
//               hashtable[index] = neighbor_clusterId;
//               found_cluster = true;
//            }
//         }
//      }
//   
//      // Check if the velocity cell was a part of a cluster, if not then make it a new cluster:
//      if( found_cluster == false ) {
//         // New cluster, so new cluster id
//         ++clusterId;
//         // Set the cluster id
//         hashtable[index] = clusterId;
//      }
//   }
//   
//   // Print how many cluster ids:
//   cerr << "CLUSTER IDS: " <<  clusterId << endl;
//   return;
//}



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
   //const size_t startingPoint = (size_t)(static_cast<const void*>((cell->block_data).data()));
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
//   // Return value:
//   Real value_to_return = 0;
//   for( unsigned int i = 0; i < velocityCells.size(); ++i ) {
//      if( i%2 == 0 ) {
//         value_to_return += (Real)(velocityCells[i].get_avgs());
//      } else {
//         value_to_return -= (Real)(velocityCells[i].get_avgs());
//      }
//   }


//   // Start getting velocity neighbors
//   vector<Velocity_Cell> neighbors;
//   neighbors.reserve( numberOfVCellNeighbors );
//   
//
//   const size_t startingPoint = (size_t)(cell->block_data.data());
//
//   
//   for( int i = velocityCells.size()-1; i >= 0; --i ) {
//      const Velocity_Cell & vCell = velocityCells[i];
//      // Get the neighbors:
//      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
//      if( vCell.hash( startingPoint ) == 0 ) {
//         cerr << "HIT" << endl;
//      } else if( vCell.hash( startingPoint ) == velocityCells.size()-1 ) {
//         cerr << "HIT2" << endl;
//      } else if( vCell.hash( startingPoint ) < 0 ) {
//         cerr << "SOMETHING WRONG " << vCell.hash( startingPoint ) << endl;
//      } else if( vCell.hash( startingPoint ) >= velocityCells.size() ) {
//         cerr << "SOMETHING WRONG2 " << vCell.hash( startingPoint ) << " "  << velocityCells.size() << endl;
//      }
//      for( vector<Velocity_Cell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
//         
//         value_to_return += it->hash(startingPoint)*0.01;
//         //cerr << it->hash(startingPoint);
//      }
//      neighbors.clear();
//   }
//   

   //cluster( velocityCells, local_vcell_neighbors, remote_vcell_neighbors, cell );
   cluster_simple( velocityCells, local_vcell_neighbors, remote_vcell_neighbors, cell );

   return 0;
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



