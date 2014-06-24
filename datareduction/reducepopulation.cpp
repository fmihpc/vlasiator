#include <iostream>
#include <utility>
#include <unordered_set>
#include <unordered_map>
#include <parallel/algorithm>
#include <algorithm>
#include <deque>
#include "phiprof.hpp"
#include "reducepopulation.h"

using namespace std;
using namespace spatial_cell;
using namespace vlsv;


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

class Cluster {
   private:
   public:
      unordered_set<uint32_t> neighbor_clusters; // Neighbor cluster indexes
      uint32_t members; // Number of members in this cluster
      uint32_t clusterId;
};

class Cluster_Fast {
   private:
   public:
      unordered_set<uint32_t> * neighbor_clusters; // Neighbor cluster indexes
      uint32_t * members; // Number of members in this cluster
      uint32_t * clusterId; // This cluster's cluster id

      // Constructor
      Cluster_Fast() {
         neighbor_clusters = NULL;
         members = NULL;
         clusterId = NULL;
      }
      // Destructor
      ~Cluster_Fast() {
         neighbor_clusters = NULL;
         members = NULL;
         clusterId = NULL;
      }

      // Compare values
      bool operator<( const Cluster_Fast & other ) const {
         return *clusterId < *other.clusterId;
      }

      // Pseudo-Constructor:
      inline void set_data( const uint32_t __clusterId ) {
         // Reserve space
         neighbor_clusters = new unordered_set<uint32_t>;
         members = new uint32_t;
         clusterId = new uint32_t;
         // Construct:
         // One member by default
         *members = 1;
         // Input cluster id:
         *clusterId = __clusterId;
         neighbor_clusters->insert( __clusterId );
      }

      // Pseudo-Destructor:
      inline void remove_data( vector<Cluster_Fast> & clusters ) {
         // Check if there's anything to delete:
         if( neighbor_clusters && clusterId && members ) {
            // Get pointers to be deleted:
            uint32_t * members_to_delete = members;
            uint32_t * clusterId_to_delete = clusterId;
            unordered_set<uint32_t> * neighbor_clusters_to_delete = neighbor_clusters;
   
            // Go through every member:
            for( unordered_set<uint32_t>::const_iterator it = neighbor_clusters->begin(); it != neighbor_clusters->end(); ++it ) {
               const uint index = *it;
               // Fix the members to null
               Cluster_Fast & cluster = clusters[index];
               cluster.members = NULL;
               cluster.clusterId = NULL;
               cluster.neighbor_clusters = NULL;
            }

            // Remove data:
            delete members_to_delete;
            delete clusterId_to_delete;
            delete neighbor_clusters_to_delete;
         }
      }

      inline bool find( const uint32_t __clusterId ) {
         phiprof_assert( neighbor_clusters );
         return neighbor_clusters->find( __clusterId ) == neighbor_clusters->end();
      }

      inline void merge( 
                  Cluster_Fast & cluster_neighbor,
                  vector<Cluster_Fast> & clusters
                            ) {

         // Make sure this class has no null pointers:
         phiprof_assert( neighbor_clusters && members && clusterId );

         phiprof_assert( cluster_neighbor.neighbor_clusters && cluster_neighbor.members && cluster_neighbor.clusterId );

         phiprof_assert( cluster_neighbor.clusterId != clusterId );

         // Optimize:
         if( cluster_neighbor.neighbor_clusters->size() > neighbor_clusters->size() ) {
            // It's faster to use this function the other way around if the neighbor is a bigger cluster:
            cluster_neighbor.merge( *this, clusters );
            return;
         }
  
         // CREATE NEW VALUES THAT WILL BE SHARED BY BOTH CLUSTERS
         //*****************
         // Create new cluster neighbors:
         unordered_set<uint32_t> * new_neighbor_clusters = neighbor_clusters;

         // Append values
         new_neighbor_clusters->insert( cluster_neighbor.neighbor_clusters->begin(), cluster_neighbor.neighbor_clusters->end() );

         // Create new cluster id:
         uint32_t * new_cluster_id = clusterId;
         *new_cluster_id = min( *clusterId, *cluster_neighbor.clusterId );

         // Get new members:
         uint32_t * new_members = members;
         *new_members = *members + *cluster_neighbor.members;
         //*****************



         // Get the pointers for the old values ( to be deleted ):
         //*****************
         unordered_set<uint32_t> * old_neighbor_clusters_neighbor = cluster_neighbor.neighbor_clusters;
         uint32_t * old_cluster_id_neighbor = cluster_neighbor.clusterId;
         uint32_t * old_members_neighbor = cluster_neighbor.members;
         //*****************





         // Update the cluster neighbors for every cluster:
         //----------------------------------------------------
         for( unordered_set<uint32_t>::const_iterator jt = old_neighbor_clusters_neighbor->begin(); jt != old_neighbor_clusters_neighbor->end(); ++jt ) {
            const uint iterator_index = *jt;
 
            phiprof_assert( iterator_index < clusters.size() );
 
            // Update the neighbors
            clusters[iterator_index].neighbor_clusters = new_neighbor_clusters;
            // Update members:
            clusters[iterator_index].members = new_members;
            // Update cluster ids:
            clusters[iterator_index].clusterId = new_cluster_id;
         }
         //----------------------------------------------------



         // Delete the old values:
         delete old_neighbor_clusters_neighbor;
         old_neighbor_clusters_neighbor = NULL;
         delete old_cluster_id_neighbor;
         old_cluster_id_neighbor = NULL;
         delete old_members_neighbor;
         old_members_neighbor = NULL;
      }

      inline void append( const int32_t numberOfMembersToAppend ) {
         *members = *members + numberOfMembersToAppend;
      }
};


class Cluster_Fast_Two {
   private:
   public:
      unordered_map<uint32_t, uint32_t> * neighbor_clusters; // Neighbor cluster indexes, connectivity
      uint32_t * members; // Number of members in this cluster
      uint32_t * clusterId; // This cluster's cluster id

      // Constructor
      Cluster_Fast_Two() {
         neighbor_clusters = NULL;
         members = NULL;
         clusterId = NULL;
      }
      // Destructor
      ~Cluster_Fast_Two() {
         neighbor_clusters = NULL;
         members = NULL;
         clusterId = NULL;
      }

      // Pseudo-Constructor:
      inline void set_data( const uint32_t __clusterId ) {
         // Reserve space
         neighbor_clusters = new unordered_map<uint32_t, uint32_t>;
         members = new uint32_t;
         clusterId = new uint32_t;
         // Construct:
         // One member by default
         *members = 1;
         // Input cluster id:
         *clusterId = __clusterId;
         neighbor_clusters->insert( make_pair(__clusterId, 0) );
      }

      // Pseudo-Destructor:
      inline void remove_data( vector<Cluster_Fast_Two> & clusters ) {
         // Check if there's anything to delete:
         if( neighbor_clusters && clusterId && members ) {
            // Get pointers to be deleted:
            uint32_t * members_to_delete = members;
            uint32_t * clusterId_to_delete = clusterId;
            unordered_map<uint32_t, uint32_t> * neighbor_clusters_to_delete = neighbor_clusters;
   
            // Go through every member:
            for( unordered_map<uint32_t, uint32_t>::const_iterator it = neighbor_clusters->begin(); it != neighbor_clusters->end(); ++it ) {
               const uint index = it->first;
               // Fix the members to null
               Cluster_Fast_Two & cluster = clusters[index];
               cluster.members = NULL;
               cluster.clusterId = NULL;
               cluster.neighbor_clusters = NULL;
            }

            // Remove data:
            delete members_to_delete;
            delete clusterId_to_delete;
            delete neighbor_clusters_to_delete;
         }
      }

      inline bool find( const uint32_t __clusterId ) {
         phiprof_assert( neighbor_clusters );
         return neighbor_clusters->find( __clusterId ) == neighbor_clusters->end();
      }

      inline void merge( 
                  Cluster_Fast_Two & cluster_neighbor,
                  vector<Cluster_Fast_Two> & clusters
                            ) {

         // Make sure this class has no null pointers:
         phiprof_assert( neighbor_clusters && members && clusterId );

         phiprof_assert( cluster_neighbor.neighbor_clusters && cluster_neighbor.members && cluster_neighbor.clusterId );

         phiprof_assert( cluster_neighbor.clusterId != clusterId );

         // Optimize:
         if( cluster_neighbor.neighbor_clusters->size() > neighbor_clusters->size() ) {
            // It's faster to use this function the other way around if the neighbor is a bigger cluster:
            cluster_neighbor.merge( *this, clusters );
            return;
         }
  
         // CREATE NEW VALUES THAT WILL BE SHARED BY BOTH CLUSTERS
         //*****************
         // Create new cluster neighbors:
         unordered_map<uint32_t, uint32_t> * new_neighbor_clusters = neighbor_clusters;

         // Append values
         new_neighbor_clusters->insert( cluster_neighbor.neighbor_clusters->begin(), cluster_neighbor.neighbor_clusters->end() );

         // Create new cluster id:
         uint32_t * new_cluster_id = clusterId;
         *new_cluster_id = min( *clusterId, *cluster_neighbor.clusterId );

         // Get new members:
         uint32_t * new_members = members;
         *new_members = *members + *cluster_neighbor.members;
         //*****************



         // Get the pointers for the old values ( to be deleted ):
         //*****************
         unordered_map<uint32_t, uint32_t> * old_neighbor_clusters_neighbor = cluster_neighbor.neighbor_clusters;
         uint32_t * old_cluster_id_neighbor = cluster_neighbor.clusterId;
         uint32_t * old_members_neighbor = cluster_neighbor.members;
         //*****************





         // Update the cluster neighbors for every cluster:
         //----------------------------------------------------
         for( unordered_map<uint32_t, uint32_t>::const_iterator jt = old_neighbor_clusters_neighbor->begin(); jt != old_neighbor_clusters_neighbor->end(); ++jt ) {
            const uint iterator_index = jt->first;
 
            phiprof_assert( iterator_index < clusters.size() );
 
            // Update the neighbors
            clusters[iterator_index].neighbor_clusters = new_neighbor_clusters;
            // Update members:
            clusters[iterator_index].members = new_members;
            // Update cluster ids:
            clusters[iterator_index].clusterId = new_cluster_id;
         }
         //----------------------------------------------------



         // Delete the old values:
         delete old_neighbor_clusters_neighbor;
         old_neighbor_clusters_neighbor = NULL;
         delete old_cluster_id_neighbor;
         old_cluster_id_neighbor = NULL;
         delete old_members_neighbor;
         old_members_neighbor = NULL;
      }

      inline void halfmerge( 
                  Cluster_Fast_Two & cluster_neighbor,
                  vector<Cluster_Fast_Two> & clusters
                            ) {

         // Make sure this class has no null pointers:
         phiprof_assert( neighbor_clusters && members && clusterId );

         phiprof_assert( cluster_neighbor.neighbor_clusters && cluster_neighbor.members && cluster_neighbor.clusterId );

         phiprof_assert( cluster_neighbor.clusterId != clusterId );

         // Check whether the neighbor cluster has already been added as a neighbor
         if( neighbor_clusters->find( *cluster_neighbor.clusterId ) == neighbor_clusters->end() ) {
            neighbor_clusters->insert( make_pair( *cluster_neighbor.clusterId, 0 ) );
            cluster_neighbor.neighbor_clusters->insert( make_pair( *clusterId, 0 ) );
         }
         // Append more neighbors: Basically what happens here is that we add neighbor to this cluster's neighbor list and vice versa for the neighbour cluster
         neighbor_clusters->operator[](*cluster_neighbor.clusterId) = neighbor_clusters->operator[](*cluster_neighbor.clusterId) + 1;
         cluster_neighbor.neighbor_clusters->operator[](*clusterId) = cluster_neighbor.neighbor_clusters->operator[](*clusterId) + 1;
      }

      inline void append( const uint32_t numberOfMembersToAppend ) {
         *members = *members + numberOfMembersToAppend;
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


/*! A function for calculating local and remote velocity cells for each velocity cell within a block. Remote velocity cells are the cells that would be outside a block and local velocity cells are the cells that are inside a block. Note: This is currently only used in the population_algorithm as a way to optimize the code. This way we don't need to calculate the neighbors for each velocity cell every time, and the memory consumption is negligible. After running this one can retrieve a given velocity cell's neighbors by the velocity cell id (0..63) with local_vcell_neighbors[vCellId] and remote neighbors likewise with remote_vcell_neighbors[vCellId]. The remote velocity cells are stored in pairs of blockid and list of neighbors within that block. So, for instance if remote_vcell_neighbors[vCellId] returned this: [<2,[0,4]> <10,[3,4]>] where <,> signals a pair, then the vCellId at block blockId would have remtoe neighbors at block blockId+2 at velocity cells 0 and 4 and at block blockId+10 at velocity cells 3 and 4

 \param local_vcell_neighbors                 Empty vector where local velocity cells get stored.  local_vcell_neighbors[vCellId] gives the local neighbors of a velocity cell (neighbors that are within the same block)
 \param remote_vcell_neighbors                Empty vector where remote velocity cells get stored. remote_vcell_neighbors[vCellId] Gives a vector containing remote neighbors of the velocity cell (neighbors that are outside the block) in vector< pair<int16_t, vector<uint16_t> > > format. The pair's first index gives the neighbor block index (check spatial_cell.hpp for more information) and the second index gives the local velocity cells withing that neighbor block.

 */
void set_local_and_remote_velocity_cell_neighbors( 
       array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
       array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
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
            int neighbor_block_direction[numberOfDirections] = {0,0,0};
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
            // Now get the neighbor block index:
            const int16_t neighbor_block_index = neighbor_block_direction[0]
                                    + neighbor_block_direction[1] * SpatialCell::vx_length
                                    + neighbor_block_direction[2] * SpatialCell::vx_length * SpatialCell::vy_length;
            const uint16_t neighbor_vCellId = neighbor_indices[0] 
                                            + neighbor_indices[1] * block_vx_length
                                            + neighbor_indices[2] * block_vx_length * block_vy_length;
            // Add the neighbor to remote velocity cell neighbors:
            // First check if the velocity block is already within the vector
            int index = -1;
            int iterator = 0;
            for( vector< pair<int16_t, vector<uint16_t> > >::iterator it = remote_vcell_neighbors[vCellId].begin();
                 it != remote_vcell_neighbors[vCellId].end();
                 ++it ) {
               // Check for velocity block
               const uint16_t iterated_neighbor_block_index = get<0>(*it);
               if( iterated_neighbor_block_index == neighbor_block_index ) {
                  // Found the neighbor index:
                  index = iterator;
               }
               ++iterator;
            }
            // Check if the velocity block was found, and if it was add one more velocity cell into the block
            if( index == -1 ) {
               // Velocity block was not found so add it to the list
               vector<uint16_t> neighbor_vcells;
               neighbor_vcells.reserve(1);
               neighbor_vcells.push_back( neighbor_vCellId );
               pair<int16_t, vector<uint16_t> > blockAndVCell = make_pair( neighbor_block_index, neighbor_vcells );
               // Add pair to the remote velocity cells
               remote_vcell_neighbors[vCellId].reserve( remote_vcell_neighbors[vCellId].size() + 1 );
               remote_vcell_neighbors[vCellId].push_back( blockAndVCell );
            } else {
               // Get the pair
               pair<int16_t, vector<uint16_t> > & blockAndVCell = remote_vcell_neighbors[vCellId][index];
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



//// A function for retrieving the velocity cell neighbors of a given velocity cell. Note: There are alwayus 3*3*3-1=26 neighbors
//const static int numberOfVCellNeighbors = 26;
//static inline void get_neighbors(
//                vector<Velocity_Cell> & neighbors,
//                const Velocity_Cell & vCell,
//                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
//                const array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
//                const SpatialCell * cell
//                                ) {
//   // Get the local neighbors:
//   const uint16_t vCellId = vCell.vCellId;
//   for( vector<uint16_t>::const_iterator it = local_vcell_neighbors[vCellId].begin(); it != local_vcell_neighbors[vCellId].end(); ++it ) {
//      const uint16_t neighbor_vCellId = *it;
//      // Add to the next neighbor:
//      Velocity_Cell neighbor;
//      // Still in the same block , new vCellId
//      neighbor.set_data( vCell.block, neighbor_vCellId );
//      neighbors.push_back( neighbor );
//   }
//   // Get the remote neighbors
//   for( vector< pair<int16_t, vector<uint16_t> > >::const_iterator it = remote_vcell_neighbors[vCellId].begin(); it != remote_vcell_neighbors[vCellId].end(); ++it ) {
//      // Get the neighbor block's index
//      const uint16_t neighbor_block_index = get<0>(*it);
//      // Go through the local cells in the neighbor block:
//      for( vector<uint16_t>::const_iterator jt = local_vcell_neighbors[vCellId].begin(); jt != local_vcell_neighbors[vCellId].end(); ++jt ) {
//         const uint16_t neighbor_vCellId = *jt;
//         // Make sure the neighbor is valid
//         if( vCell.block->neighbors[neighbor_block_index] && vCell.block->neighbors[neighbor_block_index]->data != cell->null_block_data.data() ) {
//            Velocity_Cell neighbor;
//            neighbor.set_data( vCell.block->neighbors[neighbor_block_index], neighbor_vCellId );
//            neighbors.push_back( neighbor );
//         }
//      }
//   }
//   return;
//}

// A function for retrieving the velocity cell neighbors of a given velocity cell. Note: There are alwayus 3*3*3-1=26 neighbors
const static int numberOfVCellNeighbors = 26;
static inline void get_neighbors(
                vector<Velocity_Cell> & neighbors,
                const Velocity_Cell & vCell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
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
   for( vector< pair<int16_t, vector<uint16_t> > >::const_iterator it = remote_vcell_neighbors[vCellId].begin(); it != remote_vcell_neighbors[vCellId].end(); ++it ) {
      // Get the neighbor block's index
      const uint32_t blockid = cell->get_velocity_block( vCell.block->parameters[BlockParams::VXCRD], 
                                                         vCell.block->parameters[BlockParams::VYCRD],
                                                         vCell.block->parameters[BlockParams::VZCRD] );
      const uint32_t neighbor_block_id = blockid + get<0>(*it);
      // Go through the local cells in the neighbor block:
      for( vector<uint16_t>::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt ) {
         const uint16_t neighbor_vCellId = *jt;
         const Velocity_Block * neighbor_block = cell->at(neighbor_block_id);
         // Make sure the neighbor is valid
         if( cell->at(neighbor_block_id) && neighbor_block->data != cell->null_block_data.data() ) {
            Velocity_Cell neighbor;
            neighbor.set_data( neighbor_block, neighbor_vCellId );
            neighbors.push_back( neighbor );
         }
      }
   }
   return;
}


// A function for clustering velocity cells into separate populations.
// This algorithm starts from the highest-valued velocity cell and starts clustering until it hits a certain threshold value, after which it starts to find new populations.
static inline void cluster_simple( 
                                vector<Velocity_Cell> & velocityCells,
                                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                                const array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
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

   uint debugging_iterator = 0;

   // Start inputting velocity cells:
   for( uint i = 0; i < velocityCells.size(); ++i ) {

      // Get the velocity cell:
      const Velocity_Cell & vCell = all_neighbors.front();

      // Get neighbors
      get_neighbors( tmp_neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
      for( uint j = 0; j < tmp_neighbors.size(); ++j ) {
         const uint32_t id = tmp_neighbors[j].hash( startingPoint );

         phiprof_assert( id < velocityCells.size() );

         // If the velocity cell is not a part of a cluster, add it
         // NOTE: This is probably not needed
         if( all_vCells[ id ] == noCluster ) {
            // Make sure the threshold holds:
            const Realf avgs_value = tmp_neighbors[j].get_avgs();
            if( avgs_value >= threshold ) {
               // The velocity cell has a bigger value than the threshold so add it
               all_vCells[ id ] = clusterId;
               all_neighbors.push_back( tmp_neighbors[j] );
               if( debugging_iterator >= velocityCells.size() ) {
                  cerr << "BUG FOUND2" << endl;
                  exit(1);
               }
               ++debugging_iterator;
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
         for( ; rit != velocityCells.crend(); ++rit ) {
            const uint32_t id = rit->hash( startingPoint );
            if( all_vCells[ id ] == noCluster ) {
               const Realf avgs_value = rit->get_avgs();
               if( avgs_value >= threshold ) {
                  all_vCells[ id ] = clusterId;
                  all_neighbors.push_back( *rit );
                  if( debugging_iterator >= velocityCells.size() ) {
                     cerr << "BUG FOUND" << endl;
                     exit(1);
                  }
                  ++debugging_iterator;
                  break;
               } else {
                  //rit = velocityCells.rend();
               }
            }
         }
         // If no velocity cell was found then break the algorithm:
         if( all_neighbors.size() <= 0 ) { break; }
      }
   }

   phiprof_assert( cell->block_fx.size() >= velocityCells.size() );

   for( uint i = 0; i < velocityCells.size(); ++i ) {
      cell->block_fx.at(i) = all_vCells[i];
   }


   delete[] all_vCells;

   return;
}




// A function for clustering velocity cells into separate populations.
// This algorithm uses Cluster class to save different clusters. It's identical to cluster_advanced but this one does not use pointers, which is why this one is slower.
static inline void cluster( 
                                const vector<Velocity_Cell> & velocityCells,
                                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                                const array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
                                SpatialCell * cell
                          ) {

   if( velocityCells.size() == 0 ) { return; }
   // Reserve a table for clusters:
   uint32_t * clusterIds = new uint32_t[velocityCells.size()];
   
   // Initialize to be part of no clusters:
   const uint32_t noCluster = 0;
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      clusterIds[i] = noCluster;
   }


   // Id for separating clusterIds
   uint32_t clusterId = noCluster + 1;

   // Set the first velocity cell to cluster one
   uint32_t last_vCell = velocityCells.size()-1;
   clusterIds[last_vCell] = clusterId;


   // Start getting velocity neighbors
   vector<Velocity_Cell> neighbors;
   neighbors.reserve( numberOfVCellNeighbors );

   const size_t startingPoint = (size_t)(cell->block_data.data());

   uint32_t merges = 0;

   //----------------------------------------------
   // Create clusters
   vector<Cluster> clusters;
   const uint estimated_clusters = 200;
   clusters.reserve( estimated_clusters );
   {
      // Create the no-cluster
      Cluster no_cluster;
      no_cluster.members = 0;
      no_cluster.clusterId = noCluster;
      no_cluster.neighbor_clusters.insert( noCluster );
      clusters.push_back( no_cluster );
      // Create the first cluster:
      Cluster first_cluster;
      first_cluster.members = 1;
      first_cluster.clusterId = clusterId;
      first_cluster.neighbor_clusters.insert( clusterId );
      clusters.push_back( first_cluster );
   }
   //----------------------------------------------


   for( int i = velocityCells.size()-1; i >= 0; --i ) {
      const Velocity_Cell & vCell = velocityCells[i];
      // Get the neighbors:
      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
      // Get the id of the velocity cell
      const uint32_t id = vCell.hash( startingPoint );

      phiprof_assert( id < velocityCells.size() );

      for( vector<Velocity_Cell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
         // Get the id of the neighbor:
         const uint32_t neighbor_id = it->hash( startingPoint );

         phiprof_assert( neighbor_id < velocityCells.size() );

         phiprof_assert( clusterIds[id] < clusters.size() && clusterIds[neighbor_id] < clusters.size() );

         // Set the id to equal the same as neighbors' if the neighbor is part of a cluster
         if( clusterIds[neighbor_id] != noCluster ) {
            const uint index = clusterIds[id];
            const uint index_neighbor = clusterIds[neighbor_id];

            phiprof_assert( index < clusters.size() && index_neighbor < clusters.size() );

            // If the cluster id has not been set yet, set it now:
            if( clusterIds[id] == noCluster ) {
               // Cluster id has not been set yet
               clusterIds[id] = clusterIds[neighbor_id];

               //--------------------------------------------------------------
               // Increase the amount of members in the cluster by one:
               unordered_set<uint32_t> & neighbor_clusters = clusters[index_neighbor].neighbor_clusters;
               // TODO: OPTIMIZE THIS WITH POINTERS WITH POINTERS O(1)
               for( unordered_set<uint32_t>::const_iterator jt = neighbor_clusters.begin(); jt != neighbor_clusters.end(); ++jt ) {
                  const uint iterator_index = *jt;

                  phiprof_assert( iterator_index < clusters.size() );

                  clusters[iterator_index].members += 1;
               }
               //--------------------------------------------------------------


            } else if( clusters[index].neighbor_clusters.find( index_neighbor ) == clusters[index].neighbor_clusters.end() ) {
               // Clusters are separate clusters

               // Merge the clusters:

               //Set clusters to be the same
               //----------------------------------------------------
               const uint index = clusterIds[id];
               Cluster & cluster = clusters[index];
               const uint index_neighbor = clusterIds[neighbor_id];
               Cluster & cluster_neighbor = clusters[index_neighbor];

               // Create new cluster neighbors:
               unordered_set<uint32_t> neighbor_clusters;
               neighbor_clusters.insert( cluster.neighbor_clusters.begin(), cluster.neighbor_clusters.end() );
               neighbor_clusters.insert( cluster_neighbor.neighbor_clusters.begin(), cluster_neighbor.neighbor_clusters.end() );

               // NOTE: IF THE CLUSTERS HAVE BEEN MERGED THEY WILL ALWAYS HAVE THE SAME MEMBERS AND NEIGHBOR_CLUSTERS SO POINTERS SHOULD BE USED


               // Update the cluster neighbors for every cluster:
               // TODO: OPTIMIZE, THIS ONLY NEEDS TO BE DONE ONCE WITH POINTERS (NO COPY CONSTRUCTION), USE POINTERS
               // Get the new amount of members in the merged cluster:
               const uint members = cluster.members + cluster_neighbor.members;
               // Get the minimum cluster id and set it as the cluster id for every cluster
               const uint32_t new_cluster_id = *min_element( neighbor_clusters.begin(), neighbor_clusters.end() );
               for( unordered_set<uint32_t>::const_iterator jt = neighbor_clusters.begin(); jt != neighbor_clusters.end(); ++jt ) {
                  const uint iterator_index = *jt;

                  phiprof_assert( iterator_index < clusters.size() );

                  // Update the neighbors
                  clusters[iterator_index].neighbor_clusters = neighbor_clusters;
                  // Update members:
                  clusters[iterator_index].members = members;
                  // Update cluster ids:
                  clusters[iterator_index].clusterId = new_cluster_id;
               }
               //----------------------------------------------------

               ++merges;
            }
         }
      }
      // If the cell does not have any neighbors that are part of a cluster then this is a new cluster:
      if( clusterIds[id] == noCluster ) {
         ++clusterId;

         // Create a new cluster
         //------------------------------------------------
         Cluster new_cluster;
         new_cluster.members = 1;
         new_cluster.clusterId = clusterId;
         new_cluster.neighbor_clusters.insert( clusterId );
         clusters.push_back( new_cluster );
         //------------------------------------------------

         clusterIds[id] = clusterId;
      }
      neighbors.clear();
   }

   phiprof_assert( cell->block_fx.size() >= velocityCells.size() );
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      const Realf value = clusters[clusterIds[i]].clusterId;
      cell->block_fx.at(i) = value;
   }

   // Print out the number of clusterIds:
//   cerr << "Clusters: " << clusterId << endl;
//   cerr << "Merges: " << merges << endl;
////   cerr << "Sizeof: " << sizeof(Realf) << endl;

//   

   delete[] clusterIds;
   return;
}


// A function for clustering velocity cells into separate populations.
// This algorithm uses Cluster_Fast class to save different clusters. The algorithm starts from the highest-valued velocity cell, moves to the second-highest and puts the second-highest into the same cluster as its neighbor, if its neighbor is a part of a cluster.
static inline void cluster_advanced( 
                                const vector<Velocity_Cell> & velocityCells,
                                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                                const array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
                                SpatialCell * cell,
                                const Real resolution_threshold
                                  ) {
   if( velocityCells.size() == 0 ) { return; }
   const uint numberOfVCells = velocityCells.size();
   // Reserve a table for clusters:
   //uint32_t * clusterIds = new uint32_t[numberOfVCells];
   vector<Realf,aligned_allocator<Realf,64> > & clusterIds = cell->block_fx;

   // Initialize to be part of no clusters:
   const uint32_t noCluster = 0;
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      clusterIds[i] = noCluster;
   }

   // Id for separating clusterIds
   uint32_t clusterId = noCluster + 1;

   // Set the first velocity cell to cluster one
   uint32_t last_vCell = velocityCells.size()-1;
   clusterIds[last_vCell] = clusterId;

   // Start getting velocity neighbors
   vector<Velocity_Cell> neighbors;
   neighbors.reserve( numberOfVCellNeighbors );

   const size_t startingPoint = (size_t)(cell->block_data.data());

   uint32_t merges = 0;

   //----------------------------------------------
   // Create clusters
   vector<Cluster_Fast> clusters;
   const uint estimated_clusters = 200;
   clusters.reserve( estimated_clusters );
   {
      // Create the no-cluster
      Cluster_Fast no_cluster;
      no_cluster.set_data( noCluster );
      clusters.push_back( no_cluster );
      // Create the first cluster:
      Cluster_Fast first_cluster;
      first_cluster.set_data( clusterId );
      clusters.push_back( first_cluster );
   }
   //----------------------------------------------

   const uint resolution = P::vxblocks_ini * P::vyblocks_ini * P::vzblocks_ini;

   for( int i = velocityCells.size()-1; i >= 0; --i ) {
      const Velocity_Cell & vCell = velocityCells[i];
      // Get the neighbors:
      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
      // Get the id of the velocity cell
      const uint32_t id = vCell.hash( startingPoint );

      phiprof_assert( id < velocityCells.size() );

      // Keep track of the highest avgs of the neighbor:
      //TODO: Remove highest avgs neighbor
      Realf highest_avgs_neighbor = 0;

      for( vector<Velocity_Cell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
         // Get the id of the neighbor:
         const uint32_t neighbor_id = it->hash( startingPoint );

         phiprof_assert( neighbor_id < velocityCells.size() );

         phiprof_assert( clusterIds[id] < clusters.size() && clusterIds[neighbor_id] < clusters.size() );


         // Set the id to equal the same as neighbors' if the neighbor is part of a cluster
         if( clusterIds[neighbor_id] != noCluster ) {

            const uint index = clusterIds[id];
            const uint index_neighbor = clusterIds[neighbor_id];

            phiprof_assert( index < clusters.size() && index_neighbor < clusters.size() );
            phiprof_assert( clusters[index].neighbor_clusters );


            // If the cluster id has not been set yet, set it now:
            if( clusterIds[id] == noCluster ) {

               // Cluster id has not been set yet
               clusterIds[id] = clusterIds[neighbor_id];

               //--------------------------------------------------------------
               // Increase the amount of members in the cluster by one:
               clusters[index_neighbor].append(1);
               //--------------------------------------------------------------

               // Update the highest avgs value counter:
               highest_avgs_neighbor = it->get_avgs();


            } else if( clusters[index].find( index_neighbor ) ) {

               // Clusters are separate clusters

               // Merge the clusters:

               //Set clusters to be the same (Merge them)
               //----------------------------------------------------
               Cluster_Fast & cluster_neighbor = clusters[index_neighbor];
               Cluster_Fast & cluster = clusters[index];

               // Merge clusters if the clusters are ok:
               const uint32_t clusterMembers = *cluster.members;
               const uint32_t neighborClusterMembers = *cluster_neighbor.members;
               // Check if the clusters should be merged: (If not, then check if the velocity cell should belong to the other cluster instead)
               if( clusterMembers < resolution*resolution_threshold || neighborClusterMembers < resolution*resolution_threshold ) {
                  // The other cluster is small enough, so merge the clusters:
                  cluster_neighbor.merge( cluster, clusters );
                  ++merges;
                  break;
               }
               else if( highest_avgs_neighbor < it->get_avgs() ) {
                  // The velocity cell should belong to the other cluster, because the other cluster has a highest avgs value:
                  // Remove it from the previous cluster
                  clusters[index].append(-1);
                  // Set to be the same clusters:
                  clusterIds[id] = clusterIds[neighbor_id];
                  // Increase the amount of members in the cluster by one
                  clusters[index_neighbor].append(1);

                  // Update the highest avgs value counter:
                  highest_avgs_neighbor = it->get_avgs();
               }
            }
         }
      }

      // If the cell does not have any neighbors that are part of a cluster then this is a new cluster:
      if( clusterIds[id] == noCluster ) {
         ++clusterId;

         // Create a new cluster
         //------------------------------------------------
         Cluster_Fast new_cluster;
         new_cluster.set_data( clusterId );
         clusters.push_back( new_cluster );
         //------------------------------------------------

         clusterIds[id] = clusterId;
      }

      neighbors.clear();
   }

   

   // Fix cluster ids so that they start from 0 and move up: now the cluster list (cluster ids) can look like this: [1 5 7 3 6 8 8 3]
   // Afterwards, it should look like this:                                                                         [1 3 5 2 4 6 6 2]
   {
      unordered_map<uint32_t, uint32_t> ids; // (id, rank)
      {
         // get sorted cluster ids:
         vector<uint32_t> tmp_clusterIds; tmp_clusterIds.reserve(clusters.size());
         {
            // Used to make sure we don't re-insert any cluster ids twice:
            unordered_set<uint32_t> tmp_clusterIds_set;
            
            for( uint i = 0; i < clusters.size(); ++i ) {
               const uint32_t id = *clusters[i].clusterId;
               if( id == 0) {continue;}
               if( tmp_clusterIds_set.find( id ) == tmp_clusterIds_set.end() ) { tmp_clusterIds_set.insert(id); tmp_clusterIds.push_back(id); }
            }
         }
         // Sort:
         sort( tmp_clusterIds.begin(), tmp_clusterIds.end() );
         // Map the ids:
         for( uint i = 0; i < tmp_clusterIds.size(); ++i ) {
            // Before inserting, check if id has already been inserted:
            const uint32_t id = tmp_clusterIds[i];
            // Map the id:
            ids.insert( make_pair( id, i ) );
         }
      }

      vector<uint32_t> clusterids_copy; clusterids_copy.resize(clusters.size());

      for( uint i = 0; i < clusters.size(); ++i ) {
         clusterids_copy[i] = *clusters[i].clusterId;
      }

      // Fix the clusterids according to the mapped values
      for( uint i = 0; i < clusters.size(); ++i ) {
         const uint32_t id = clusterids_copy[i];
         if( id == 0) {continue;}
         *clusters[i].clusterId = ids.at(id);
      }

      // Input number of populations:
      cell->number_of_populations = ids.size();
   }

   

   // Set values of clusters and calculate number of clusters
   phiprof_assert( cell->block_fx.size() >= velocityCells.size() );
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      const Realf value = *clusters[clusterIds[i]].clusterId;
      cell->block_fx[i] = value;
   }


   // Print out the number of clusterIds:
//   cerr << "Clusters: " << clusterId << endl;
//   cerr << "Merges: " << merges << endl;
//   cerr << "Sizeof: " << sizeof(Realf) << endl;


   //delete[] clusterIds;
   //clusterIds = NULL;
   // Free memory from clusters:
   for( uint i = 0; i < clusters.size(); ++i ) {
      clusters[i].remove_data( clusters );
   }
   return;
}

// A function for clustering velocity cells into separate populations.
// This algorithm uses Cluster_Fast class to save different clusters. The algorithm starts from the highest-valued velocity cell, moves to the second-highest and puts the second-highest into the same cluster as its neighbor, if its neighbor is a part of a cluster.
static inline void cluster_advanced_two( 
                                const vector<Velocity_Cell> & velocityCells,
                                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                                const array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
                                SpatialCell * cell
                                       ) {

   if( velocityCells.size() == 0 ) { return; }
   const uint numberOfVCells = velocityCells.size();
   // Reserve a table for clusters:
   uint32_t * clusterIds = new uint32_t[numberOfVCells];

   // Initialize to be part of no clusters:
   const uint32_t noCluster = 0;
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      clusterIds[i] = noCluster;
   }

   // Id for separating clusterIds
   uint32_t clusterId = noCluster + 1;

   // Set the first velocity cell to cluster one
   uint32_t last_vCell = velocityCells.size()-1;
   clusterIds[last_vCell] = clusterId;

   // Start getting velocity neighbors
   vector<Velocity_Cell> neighbors;
   neighbors.reserve( numberOfVCellNeighbors );

   const size_t startingPoint = (size_t)(cell->block_data.data());

   uint32_t merges = 0;

   //----------------------------------------------
   // Create clusters
   vector<Cluster_Fast_Two> clusters;
   const uint estimated_clusters = 200;
   clusters.reserve( estimated_clusters );
   {
      // Create the no-cluster
      Cluster_Fast_Two no_cluster;
      no_cluster.set_data( noCluster );
      clusters.push_back( no_cluster );
      // Create the first cluster:
      Cluster_Fast_Two first_cluster;
      first_cluster.set_data( clusterId );
      clusters.push_back( first_cluster );
   }
   //----------------------------------------------


   for( int i = velocityCells.size()-1; i >= 0; --i ) {
      const Velocity_Cell & vCell = velocityCells[i];
      // Get the neighbors:
      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
      // Get the id of the velocity cell
      const uint32_t id = vCell.hash( startingPoint );

      phiprof_assert( id < velocityCells.size() );

      for( vector<Velocity_Cell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {

         // Get the id of the neighbor:
         const uint32_t neighbor_id = it->hash( startingPoint );

         phiprof_assert( neighbor_id < velocityCells.size() );

         phiprof_assert( clusterIds[id] < clusters.size() && clusterIds[neighbor_id] < clusters.size() );


         // Set the id to equal the same as neighbors' if the neighbor is part of a cluster
         if( clusterIds[neighbor_id] != noCluster ) {

            const uint index = clusterIds[id];
            const uint index_neighbor = clusterIds[neighbor_id];

            phiprof_assert( index < clusters.size() && index_neighbor < clusters.size() );
            phiprof_assert( clusters[index].neighbor_clusters );

            // If the cluster id has not been set yet, set it now:
            if( clusterIds[id] == noCluster ) {

               // Cluster id has not been set yet
               clusterIds[id] = clusterIds[neighbor_id];

               //--------------------------------------------------------------
               // Increase the amount of members in the cluster by one:
               clusters[index_neighbor].append(1);
               //--------------------------------------------------------------

            } else {

               // Clusters are separate clusters

               // Merge the clusters:

               //Set clusters to be the same (Merge them)
               //----------------------------------------------------
               Cluster_Fast_Two & cluster_neighbor = clusters[index_neighbor];
               Cluster_Fast_Two & cluster = clusters[index];

               // Merge clusters if the clusters are ok:
               const uint32_t clusterMembers = *cluster.members;
               const uint32_t neighborClusterMembers = *cluster_neighbor.members;
               if( *(cluster.clusterId) != *(cluster_neighbor.clusterId) ) {
                  cluster_neighbor.halfmerge( cluster, clusters );
                  //----------------------------------------------------
                  ++merges;
               }
            }
         }
      }

      // If the cell does not have any neighbors that are part of a cluster then this is a new cluster:
      if( clusterIds[id] == noCluster ) {
         ++clusterId;

         // Create a new cluster
         //------------------------------------------------
         Cluster_Fast_Two new_cluster;
         new_cluster.set_data( clusterId );
         clusters.push_back( new_cluster );
         //------------------------------------------------

         clusterIds[id] = clusterId;
      }

      neighbors.clear();
   }


   // Do the merging:
   for( vector<Cluster_Fast_Two>::iterator it = clusters.begin(); it != clusters.end(); ++it ) {
      bool loopAgain = true;
      while( loopAgain ) {
         loopAgain = false;
         const double members = *(it->members);
         // Get exponent of members for evaluating if the connectivity is good enough to do merging
         const double M23 = pow(members, 2/3.);
         // Get measures of connectivity
         unordered_map<uint32_t, uint32_t> * connectivity = it->neighbor_clusters;
         // Loop through neighbors:
         for( unordered_map<uint32_t, uint32_t>::const_iterator jt = connectivity->begin(); jt != connectivity->end(); ++jt ) {
            // Get the connectivity:
            const double N = jt->second;
            // If connectivity is good enough, merge clusters: (also make sure that the clusters are not the same clusters)
            cerr << N/M23 << " " << 2.5 << endl;
            cerr << N << endl;
            const uint index = jt->first;
            Cluster_Fast_Two & cluster_neighbor = clusters[index];
            Cluster_Fast_Two & cluster = *it;
            if( N/M23 > 2.5 && *(cluster_neighbor.clusterId) != *(cluster.clusterId) ) {
               // Get neighbor cluster
               cluster_neighbor.merge( cluster, clusters );
               loopAgain = true;
               break;
            }
         }
      }
   }


   // Set values of clusters
   phiprof_assert( cell->block_fx.size() >= velocityCells.size() );
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      const Realf value = *clusters[clusterIds[i]].clusterId;
      cell->block_fx[i] = value;
   }


   // Print out the number of clusterIds:
//   cerr << "Clusters: " << clusterId << endl;
//   cerr << "Merges: " << merges << endl;
//   cerr << "Sizeof: " << sizeof(Realf) << endl;


   delete[] clusterIds;
   clusterIds = NULL;
   // Free memory from clusters:
   for( uint i = 0; i < clusters.size(); ++i ) {
      clusters[i].remove_data( clusters );
   }
   return;
}



//// Simplest case of clustering. Starts from highest valued velocity cell moves to the second-highest, third-highest, etc. This algorithm does not have cluster merging.
//static inline void cluster_simple_two( 
//                                const vector<Velocity_Cell> & velocityCells,
//                                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
//                                const array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
//                                SpatialCell * cell
//                          ) {
//
//   if( velocityCells.size() == 0 ) { return; }
//   // Reserve a table for clusters:
//   uint32_t * clusterIds = new uint32_t[velocityCells.size()];
//   
//   // Initialize to be part of no clusters:
//   const uint32_t noCluster = 1;
//   for( uint i = 0; i < velocityCells.size(); ++i ) {
//      clusterIds[i] = noCluster;
//   }
//
//
//   // Id for separating clusterIds
//   uint32_t clusterId=2;
//
//   // Set the first velocity cell to cluster one
//   uint32_t last_vCell = velocityCells.size()-1;
//   clusterIds[last_vCell] = clusterId;
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
//
//   for( int i = velocityCells.size()-1; i >= 0; --i ) {
//      const Velocity_Cell & vCell = velocityCells[i];
//      // Get the neighbors:
//      get_neighbors( neighbors, vCell, local_vcell_neighbors, remote_vcell_neighbors, cell );
//      // Get the id of the velocity cell
//      const uint32_t id = vCell.hash( startingPoint );
//
//      phiprof_assert( id < velocityCells.size() );
//
//      for( vector<Velocity_Cell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
//         // Get the id of the neighbor:
//         const uint32_t neighbor_id = it->hash( startingPoint );
//
//         phiprof_assert( neighbor_id < velocityCells.size() );
//
//         // Set the id to equal the same as neighbors' if the neighbor is part of a cluster
//         if( clusterIds[neighbor_id] != noCluster ) {
//            // If the cluster id has not been set yet, set it now:
//            if( clusterIds[id] == noCluster ) {
//               // Cluster id has not been set yet
//               clusterIds[id] = clusterIds[neighbor_id];
//            } else if( clusterIds[id] != clusterIds[neighbor_id] ) {
//               // Merge the clusters:
//               ++merges;
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
//
//   phiprof_assert( cell->block_fx.size() >= velocityCells.size() );
//   for( uint i = 0; i < velocityCells.size(); ++i ) {
//      const Realf value = clusterIds[i];
//      cell->block_fx.at(i) = value;
//   }
//
//   // Print out the number of clusterIds:
//   //cerr << "Clusters: " << clusterId << endl;
//   //cerr << "Merges: " << merges << endl;
////   cerr << "Sizeof: " << sizeof(Realf) << endl;
//
////   
//
//   delete[] clusterIds;
//
//   return;
//}



// Function for debugging the validity of cell neighbors function
static void test_neighbor(
                SpatialCell * cell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors
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

      for( vector< pair<int16_t, vector<uint16_t> > >::const_iterator it = remote_vcell_neighbors[i].begin();
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



/*! Function for calculating the different populations in distribution for a given spatial cell. Note that the results are currently saved in block_fx array within the SpatialCell.
 \param cell                                A velocity cell from the DCCRG grid
 \param local_vcell_neighbors               List of velocity cell neighbors within a block for each velocity cell within a block, note that this can be retrieved with the function set_local_and_remote_velocity_cell_neighbors
 \param remote_vcell_neighbors              List of velocity cell neighbors outside of a block for each velocity cell within a block, note that this can be retrieved with the function set_local_and_remote_velocity_cell_neighbors
 \param resolution_threshold                A value for determining how large at minimum we want our populations to be. 0.006 seems ok unless there's a reason to believe otherwise.

 */
void population_algorithm( 
                SpatialCell * cell,
                const array<vector<uint16_t>, VELOCITY_BLOCK_LENGTH> & local_vcell_neighbors,
                const array< vector< pair<int16_t, vector<uint16_t> > > , VELOCITY_BLOCK_LENGTH> & remote_vcell_neighbors,
                const Real resolution_threshold
                   ) {
   // Sort list of avgs values:

   // Vector for holding velocity cells:
   vector<Velocity_Cell> velocityCells;
   // Initialize avgs values vector:
   velocityCells.resize( cell->number_of_blocks * VELOCITY_BLOCK_LENGTH );


   // Input data
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

   // Test clustering:

   //cluster( velocityCells, local_vcell_neighbors, remote_vcell_neighbors, cell );
   //cluster_simple( velocityCells, local_vcell_neighbors, remote_vcell_neighbors, cell );
   //cluster_simple_two( velocityCells, local_vcell_neighbors, remote_vcell_neighbors, cell );
   cluster_advanced( velocityCells, local_vcell_neighbors, remote_vcell_neighbors, cell, resolution_threshold );
   //cluster_advanced_two( velocityCells, local_vcell_neighbors, remote_vcell_neighbors, cell );

   return;
}

// Function for getting the max number of populations:
static uint max_number_of_populations( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) {
   uint max_populations = 0;
   // Fetch cells:
   vector<uint64_t> local_cells = mpiGrid.get_cells();
   // Loop through cells:
   for( vector<uint64_t>::const_iterator it = local_cells.begin(); it != local_cells.end(); ++it ) {
      // Get cell:
      const uint64_t cellid = *it;
      // Get the cell's spatial cell class:
      SpatialCell * cell = mpiGrid[cellid];
      // Check number of populations:
      if( max_populations < cell->number_of_populations ) { max_populations = cell->number_of_populations; }
   }
   return max_populations;
}

// Function for writing out rho for different populations:
static inline bool write_rho( const uint max_populations, const vector<uint64_t> & local_cells, const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, Writer & vlsvWriter ) {
   //Declare attributes
   map<string, string> xmlAttributes;
   //We received mesh name as a parameter: MOST LIKELY THIS IS SpatialGrid!
   xmlAttributes["mesh"] = "SpatialGrid";
   xmlAttributes["name"] = "PopulationRho";
   const unsigned int arraySize = local_cells.size();
   const unsigned int vectorSize = max_populations;
   bool success = true;
//   if( arraySize == 0 ) {
//      const Real dummy_attribute = 0;
//      const uint64_t dataSize = sizeof(Real);
//      const string dataType = "float";
//      cerr << "WRITING ARRAY" << endl;
//      return vlsvWriter.writeArray("VARIABLE", xmlAttributes, dataType, arraySize, vectorSize, dataSize, (char*)(&dummy_attribute));
//   }
//   // Reserve space for writing out rho:
//   vector<Real> population_rho;
//   population_rho.resize(max_populations * local_cells.size());
//   // Initialize to zero;
//   for( uint i = 0; i < max_populations * local_cells.size(); ++i ) {
//      population_rho[i] = 0;
//   }
//   // Loop through cells:
//   for( uint i = 0; i < local_cells.size(); ++i ) {
//      // Get cell:
//      const uint64_t cellid = local_cells[i];
//      // Get the cell's spatial cell class:
//      const SpatialCell * cell = mpiGrid[cellid];
//      const Velocity_Block* block_ptr = cell->at(cell->velocity_block_list[0]);
//      const Real DV3 = block_ptr->parameters[BlockParams::DVX] * block_ptr->parameters[BlockParams::DVY] * block_ptr->parameters[BlockParams::DVZ];
//      // Loop through velocity cells and write out the population:
//      for( uint v = 0; v < cell->number_of_blocks * WID3; ++v ) {
//         population_rho[max_populations*i + (uint)(cell->block_fx[v])] += cell->block_data[v] * DV3;
//      }
//   }
//   // Print rho for debugging
//   for( uint i = 0; i < max_populations * local_cells.size(); ++i ) {
//      population_rho[i] += 4.17253e-11;
//   }
   // Write out the population:
   //success = vlsvWriter.writeArray("VARIABLE", xmlAttributes, arraySize, 1, population_rho.data());
   return success;
}

/*! Writes the population variables, so variables for different populations.

 \param mpiGrid                 The DCCRG grid with spatial cells
 \param vlsvWriter              The VLSV writer class for writing VLSV files, note that the file must have been opened already
 */
bool write_population_variables( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, Writer & vlsvWriter ) {
cerr << __LINE__ << endl;
   // Get max number of populations (to be used in determining how many arrays to write)
   const uint max_populations = max_number_of_populations( mpiGrid );
cerr << __LINE__ << endl;
   // Fetch cells:
   const vector<uint64_t> local_cells = mpiGrid.get_cells();
cerr << __LINE__ << endl;
   // Write rho:
   bool success = write_rho( max_populations, local_cells, mpiGrid, vlsvWriter );
cerr << __LINE__ << endl;
   return success;
}




/*! Writes the distribution function for different populations. Note that population_algorithm must have been called before this!

 \param mpiGrid                 The DCCRG grid with spatial cells
 \param vlsvWriter              The VLSV writer class for writing VLSV files, note that the file must have been opened already
 */
bool write_population_distribution( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const vector<uint64_t> & local_cells, vlsv::Writer & vlsvWriter ) {
   string name = "BLOCKVARIABLEPOPULATION";

   //Compute totalBlocks
   uint64_t totalBlocks = 0;  
   for(size_t cell=0;cell<local_cells.size();++cell){
      totalBlocks+=mpiGrid[local_cells[cell]]->number_of_blocks;
   }

   const string datatype_avgs = "float";
   const uint64_t arraySize_avgs = totalBlocks;
   const uint64_t vectorSize_avgs = WID3; // There are 64 elements in every velocity block
   const uint64_t dataSize_avgs = sizeof(Realf);
   map<string,string> attribs;
   attribs["mesh"] = "SpatialGrid";
   attribs["name"] = "avgs"; // Name of the velocity space distribution is written avgs

   // Start multi write
   vlsvWriter.startMultiwrite(datatype_avgs,arraySize_avgs,vectorSize_avgs,dataSize_avgs);
   // Loop over local_cells
   for( size_t cell = 0; cell < local_cells.size(); ++cell ) {
      // Get the spatial cell
      SpatialCell* SC = mpiGrid[local_cells[cell]];
      // Get the number of blocks in this cell
      const uint64_t arrayElements = SC->number_of_blocks;
      char * arrayToWrite = reinterpret_cast<char*>(SC->block_fx.data());
      // Add a subarray to write
      vlsvWriter.addMultiwriteUnit(arrayToWrite, arrayElements); // Note: We told beforehands that the vectorsize = WID3 = 64
   }
   if( local_cells.size() == 0 ) {
      // Write dummy array to avoid MPI blocking
      vlsvWriter.addMultiwriteUnit(NULL, 0);
   }
   // Write the subarrays
   vlsvWriter.endMultiwrite(name, attribs);
   return true;
}


