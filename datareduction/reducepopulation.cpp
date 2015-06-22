#include <iostream>
#include <utility>
#include <unordered_set>
#include <unordered_map>
#include <parallel/algorithm>
#include <algorithm>
#include <deque>
#include <sys/socket.h>
#include "phiprof.hpp"
#include "reducepopulation.h"

using namespace std;
using namespace spatial_cell;
using namespace vlsv;

const static uint32_t nullCluster = 0;


class VelocityCell {
   private:
      //const SpatialCell * cell;
   public:
      //uintVELOCITY_BLOCK_LENGTH_t index; //Index could be uint32_t is enough
      //uint32_t blockId;
      vmesh::GlobalID block;
      uint16_t vCellId;
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> * vmesh;
      vmesh::VelocityBlockContainer<vmesh::LocalID> * blockContainer;

      // Constructor
      inline void set_data( vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> & input_vmesh, 
                            vmesh::VelocityBlockContainer<vmesh::LocalID> & input_blockContainer,
                            const vmesh::GlobalID input_block, 
                            const uint16_t input_vCellId ) {
         block = input_block;
         vCellId = input_vCellId;
         vmesh = &input_vmesh;
         blockContainer = &input_blockContainer;
      }

      // Compare values
      bool operator<( const VelocityCell & other ) const {
         return get_avgs() < other.get_avgs();
      }

      // Compare equality
      bool operator==(  const VelocityCell & other ) const {
         // Compare memory addresses of blocks and velocity cell id
         return (block == other.block) && (vCellId == other.vCellId);
      }

      // Function for returning hash values for the velocity cell
      inline size_t hash() const {
         const size_t startingPoint = (size_t)(blockContainer->getData());
         // Return the address of block's data in size_t form plus the vCellId
         vmesh::LocalID blockLocalId = vmesh->getLocalID(block);
         blockContainer->getData(blockLocalId);
         return ((size_t)(blockContainer->getData(blockLocalId)) - startingPoint)/sizeof(Realf) + (size_t)(vCellId);
      }

      // Function for getting the avgs value
      inline Realf get_avgs() const {
         vmesh::LocalID blockLocalId = vmesh->getLocalID(block);
         Realf * array = blockContainer->getData(blockLocalId);
         return array[vCellId];
      }
};

class Cluster {
   private:
   public:
      unordered_set<uint32_t> * neighbor_clusters; // Neighbor cluster indexes
      uint32_t * members; // Number of members in this cluster
      uint32_t * clusterId; // This cluster's cluster id
      bool set_data_called = false;

      // Constructor
      Cluster() {
         neighbor_clusters = NULL;
         members = NULL;
         clusterId = NULL;
      }
      // Destructor
      ~Cluster() {
         neighbor_clusters = NULL;
         members = NULL;
         clusterId = NULL;
      }

      // Compare values
      bool operator<( const Cluster & other ) const {
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
	 
	 set_data_called = true;
      }

      // Pseudo-Destructor:
      inline void remove_data( vector<Cluster> & clusters ) {
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
               Cluster & cluster = clusters[index];
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

      // Check if a cluster id exists:
      inline bool find( const uint32_t __clusterId ) {
         phiprof_assert( neighbor_clusters );
         return neighbor_clusters->find( __clusterId ) == neighbor_clusters->end();
      }
      
      inline bool isNeighbor( const uint32_t __clusterId ) {
         phiprof_assert( neighbor_clusters );
	 return neighbor_clusters->find(__clusterId) != neighbor_clusters->end();
      }

      inline void merge( 
                  Cluster & cluster_neighbor,
                  vector<Cluster> & clusters
                            ) {
         if( neighbor_clusters == cluster_neighbor.neighbor_clusters ) return;
         // Make sure this class has no null pointers:
         phiprof_assert( neighbor_clusters != cluster_neighbor.neighbor_clusters );
         phiprof_assert( neighbor_clusters && members && clusterId );

         phiprof_assert( cluster_neighbor.neighbor_clusters && cluster_neighbor.members && cluster_neighbor.clusterId );

         phiprof_assert( cluster_neighbor.clusterId != clusterId );

         phiprof_assert( clusterId != NULL );

         phiprof_assert( cluster_neighbor.clusterId != NULL );

         phiprof_assert( *clusterId != nullCluster );

         phiprof_assert( *cluster_neighbor.clusterId != nullCluster );

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
         phiprof_assert( new_neighbor_clusters );

         // Append values
         new_neighbor_clusters->insert( cluster_neighbor.neighbor_clusters->begin(), cluster_neighbor.neighbor_clusters->end() );

         // Create new cluster id ( choose minimum of the two):
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
         
         // Just to make sure, put the neighbor_clusters to equal each other:
      }

      inline void append( const int32_t numberOfMembersToAppend ) {
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
   struct hash<VelocityCell> {
      inline size_t operator() (const VelocityCell & vcell) const {
         std::hash<size_t> size_t_hasher;
         return size_t_hasher( (size_t)(vcell.block)/sizeof(Realf)+vcell.vCellId );
      }
   };
}

inline uint16_t getVelocityCellid( const vmesh::LocalID globalVelocityCellIndices[3] ) {
  return (uint16_t)(globalVelocityCellIndices[0] % WID + (globalVelocityCellIndices[1] % WID) * WID + (globalVelocityCellIndices[2] % WID) * WID2);
}

inline void getLocalVelocityCellIndices( const uint16_t velocityCellid, vmesh::LocalID localVelocityCellIndices[3] ) {
   localVelocityCellIndices[0] = velocityCellid % WID;
   localVelocityCellIndices[1] = (velocityCellid / WID) % WID;
   localVelocityCellIndices[2] = (velocityCellid / (WID*WID) );
   phiprof_assert( velocityCellid == localVelocityCellIndices[0] + localVelocityCellIndices[1] * WID + localVelocityCellIndices[2] * WID * WID );
}

inline void getGlobalVelocityCellIndices( vmesh::VelocityMesh< vmesh::GlobalID, vmesh::LocalID > & vmesh, const VelocityCell & vCell, vmesh::LocalID globalVelocityCellIndices[3] ) {
   uint8_t refinement = 0;
   // Get blockid:
   const vmesh::GlobalID blockid = vCell.block;
   const uint16_t velocityCellid = vCell.vCellId;
   // Get block indices:
   vmesh::LocalID blockIndices[3], velocityCellIndices[3];
   vmesh.getIndices( blockid, refinement, blockIndices[0], blockIndices[1], blockIndices[2] );
   // Get velocity cell indices:
   getLocalVelocityCellIndices( velocityCellid, velocityCellIndices );
   // Input global indices:
   for( int direction = 0; direction < 3; ++direction ) {
      globalVelocityCellIndices[direction] = blockIndices[direction] * WID + velocityCellIndices[direction];
   }

   // Debugging:
   phiprof_assert( getVelocityCellid(globalVelocityCellIndices) == velocityCellid );
   phiprof_assert( vmesh.findBlockDown(refinement, globalVelocityCellIndices) );
}

/*! A function for retrieving the velocity cell neighbors of a given velocity cell.

 \param                      neighbors                           Empty vector where the neighbors of vCell gets stored
 \param                      vCell                               The velocity cell whose neighbors will be fetched
 \param                      local_vcell_neighbors               List of velocity cell neighbors within a block for each velocity cell within a block, note that this can be retrieved with the function set_local_and_remote_velocity_cell_neighbors
 \param                      remote_vcell_neighbors              List of velocity cell neighbors outside of a block for each velocity cell within a block, note that this can be retrieved with the function set_local_and_remote_velocity_cell_neighbors
 \param                      cell

 */
static inline void getNeighbors(
                vector<VelocityCell> & neighbors,
                const VelocityCell & vCell,
                SpatialCell * cell
                                ) {
   uint8_t refinement = 0;
   neighbors.clear();
   // Get the vmesh:
#warning no AMR
   const size_t notSet = 0;
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> & vmesh = cell->get_velocity_mesh(notSet);
   vmesh::VelocityBlockContainer<vmesh::LocalID> & blockContainer = cell->get_velocity_blocks(notSet);

   vmesh::LocalID globalVelocityCellIndices[3];
   getGlobalVelocityCellIndices( vmesh, vCell, globalVelocityCellIndices );
   int offsets[3];
   for( offsets[0] = -1; offsets[0] <= 1; ++offsets[0] ) for( offsets[1] = -1; offsets[1] <= 1; ++offsets[1] ) for( offsets[2] = -1; offsets[2] <= 1; ++offsets[2] ) {
      if( offsets[0] == 0 && offsets[1] == 0 && offsets[2] == 0 ) { continue; }
      // Out of bounds checking:
      bool outOfBounds = false;
      for( int direction = 0; direction < 3; ++direction ) {
         const vmesh::LocalID * gridLength = vmesh.getGridLength(0);
         if( (int)globalVelocityCellIndices[direction] + offsets[direction] < 0 ) { 
            outOfBounds = true; continue;
         }
         if( (int)globalVelocityCellIndices[direction] + offsets[direction] > gridLength[direction]*WID ) {
            outOfBounds = true; continue;
         }
      }
      if( outOfBounds ) { continue; }

      vmesh::LocalID neighborGlobalVelocityCellIndices[3];
      // Get the neighbor indices:
      for( int direction = 0; direction < 3; ++direction ) {
         neighborGlobalVelocityCellIndices[direction] = globalVelocityCellIndices[direction] + offsets[direction];
      }

      // Get neighbor data:
      const uint16_t neighborVelocityCellid = getVelocityCellid(neighborGlobalVelocityCellIndices);
      const vmesh::LocalID neighborBlock = vmesh.findBlockDown(refinement, neighborGlobalVelocityCellIndices);
      if( neighborBlock == vmesh.invalidGlobalID() ) { continue; }

      // Add the neighbor:
      VelocityCell newNeighbor;
      newNeighbor.set_data( vmesh, blockContainer, neighborBlock, neighborVelocityCellid );
      neighbors.push_back( newNeighbor );
   }
}




/*! A function that clusters velocity cells into separate populations. This algorithm uses the class Cluster to save different clusters and then writes them out into a spatial cell's block_fx values. It should be noted that this function changes spatial cell's block_fx value, so one should be careful when using this function. The algorithm starts from the highest-valued velocity cell, moves to the second-highest and puts the second highest into the same cluster as its neighbor, if its neighbor is part of a cluster.

Pseudocode:
\code{.cpp}
pre-processing:
// Sort velocity cells
velocity_cells = sort(velocity_cells)

cluster_advanced()
  // Start with no populations (or clusters, whichever we want to call them)
  populations=[]
  // Go through each velocity cell starting from highest avgs value, second highest, third highest, etc..
  for velocity_cell in velocity_cells
    neighbors = getNeighbors(velocity_cell) // Gets all of velocity cell's neighbors, in total 6 neighbors from each face
    for neighbor in neighbors // Go through every neighbor
      if neighbor in populations // Check if the neighbor has been put into some population
        if velocity_cell in populations // Check if the velocity cell is a part of a population
          // Now both velocity_cell and neighbor are in populations! Get the populations!
          velocity_cell_population = population( velocity_cell )
          neighbor_population = population( neighbor )
          // If one of the populations is really small, then merge the populations into one big population
          if neighbor_population < really_small or velocity_cell_population < really_small
             merge( velocity_cell_population and neighbor_population ) // Now all members of velocity_cell population are also members of neighbor_population, and vice-versa
        else
          // velocity_cell is not a part of a population, but the neighbor is! So add velocity_cell to the neighbor's population
          neighbor_population = population( neighbor )
          add( velocity_cell, neighbor_population ) // Now velocity cell is part of a population too
      else
        do_nothing()
   //if velocity_cell is still not part of a population, it means that none of its neighbors are part of a population, so create a new population for it!
   if velocity_cell not in populations
      new_population = create_new_population()
      populations.add( new_population )
      add( velocity_cell, new_population ) // Now velocity cell is part of a new population of velocity cells
      
\endcode

 \param velocityCells                         List of velocity cells within the spatial cell, note: this should be sorted based on the avgs value
 \param local_vcell_neighbors                 local_vcell_neighbors[vCellId] gives the local neighbors of a velocity cell (neighbors that are within the same block)
 \param remote_vcell_neighbors                remote_vcell_neighbors[vCellId] Gives a vector containing remote neighbors of the velocity cell (neighbors that are outside the block) in vector< pair<int16_t, vector<uint16_t> > > format. The pair's first index gives the neighbor block index (check spatial_cell.hpp for more information) and the second index gives the local velocity cells within that neighbor block.
 \param cell                                  The spatial cell whose populatios will be calculated
 \param resolution_threshold                A value for determining how large at minimum we want our populations to be. 0.006 seems ok unless there's a reason to believe otherwise.

 */
static inline void cluster_advanced( 
                                const vector<VelocityCell> & velocityCells,
                                SpatialCell * cell,
                                const Real resolution_threshold
                                  ) {
   
   cout << "Clustering.." << endl;
   if( velocityCells.size() == 0 ) { return; }
   const uint numberOfVCells = velocityCells.size();
   // Reserve a table for clusters:
   //uint32_t * clusterIds = new uint32_t[numberOfVCells];
   vector<Realf,aligned_allocator<Realf,64> > clusterIds; clusterIds.resize(numberOfVCells);

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
   vector<VelocityCell> neighbors;
   const size_t numberOfVCellNeighbors = 27;
   neighbors.reserve( numberOfVCellNeighbors );

   const size_t startingPoint = (size_t)(cell->get_data());

   uint32_t merges = 0;

   //----------------------------------------------
   // Create clusters
   vector<Cluster> clusters;
   const uint estimated_clusters = 200;
   clusters.reserve( estimated_clusters );
   {
      // Create the no-cluster
      Cluster no_cluster;
      no_cluster.set_data( noCluster );
      clusters.push_back( no_cluster );
      // Create the first cluster:
      Cluster first_cluster;
      first_cluster.set_data( clusterId );
      clusters.push_back( first_cluster );
      const VelocityCell & vCell = velocityCells[velocityCells.size()-1];
      const uint32_t id = vCell.hash(  );
      // Add the highest valued velocity cell to the first cluster:
      clusterIds[id] = clusterId;
   }
   //----------------------------------------------

   const uint resolution = P::vxblocks_ini * P::vyblocks_ini * P::vzblocks_ini;

   for( int i = velocityCells.size()-1; i >= 0; --i ) {
      const VelocityCell & vCell = velocityCells[i];
      // Get the neighbors:
      getNeighbors( neighbors, vCell, cell );
      // Get the id of the velocity cell
      const uint32_t id = vCell.hash(  );

      phiprof_assert( id < velocityCells.size() );

      // Keep track of the highest avgs of the neighbor:
      //TODO: Remove highest avgs neighbor
      Realf highest_avgs_neighbor = 0;

      for( vector<VelocityCell>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
         // Get the id of the neighbor:
         const uint32_t neighbor_id = it->hash(  );

         phiprof_assert( neighbor_id < velocityCells.size() );

         phiprof_assert( clusterIds[id] < clusters.size() && clusterIds[neighbor_id] < clusters.size() );


         // Set the id to equal the same as neighbors' if the neighbor is part of a cluster
         if( clusterIds[neighbor_id] != noCluster ) {

            const uint index = clusterIds[id];
            const uint index_neighbor = clusterIds[neighbor_id];

            phiprof_assert( index < clusters.size() && index_neighbor < clusters.size() );
            phiprof_assert( clusters[index].neighbor_clusters );


            // If the cluster id has not been set yet, set it now:
            if( index == noCluster ) {

               // Cluster id has not been set yet
               clusterIds[id] = clusterIds[neighbor_id];

               //--------------------------------------------------------------
               // Increase the amount of members in the cluster by one:
               clusters[index_neighbor].append(1);
               //--------------------------------------------------------------

               // Update the highest avgs value counter:
               if( it->get_avgs() > highest_avgs_neighbor) {highest_avgs_neighbor = it->get_avgs();}


            } else if( clusters[index].find( index_neighbor ) ) {

               phiprof_assert( index_neighbor != 0 );
               phiprof_assert( index != 0 );

               // Clusters are separate clusters

               // Merge the clusters:

               //Set clusters to be the same (Merge them)
               //----------------------------------------------------
               Cluster & cluster_neighbor = clusters[index_neighbor];
               Cluster & cluster = clusters[index];

               phiprof_assert( *cluster_neighbor.clusterId != 0 );
               phiprof_assert( *cluster.clusterId != 0 );

               // Merge clusters if the clusters are ok:
               const uint32_t clusterMembers = *cluster.members;
               const uint32_t neighborClusterMembers = *cluster_neighbor.members;
               // Check if the clusters should be merged: (If not, then check if the velocity cell should belong to the other cluster instead)
               if( clusterMembers < resolution*resolution_threshold || neighborClusterMembers < resolution*resolution_threshold ) {
                  cout << "Merge clusters!" << endl;
                  // The other cluster is small enough, so merge the clusters:
                  cluster_neighbor.merge( cluster, clusters );
                  ++merges;
                  break;
               } else if( highest_avgs_neighbor < it->get_avgs() ) {
                  cout << "Dont merge" << endl;
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

         phiprof_assert(clusterId != 0);

         // Create a new cluster
         //------------------------------------------------
         Cluster new_cluster;
         new_cluster.set_data( clusterId );
         clusters.push_back( new_cluster );
         //------------------------------------------------

         clusterIds[id] = clusterId;
      }

      neighbors.clear();
   }

   
#ifndef NDEBUG
int num_noclusters = 0;
for( vector<Cluster>::const_iterator it = clusters.begin(); it != clusters.end(); ++it ) {
   if( *(it->clusterId) == 0 ) { num_noclusters++; }
}
phiprof_assert( num_noclusters == 1 );
#endif

   // Fix cluster ids so that they start from 0 and move up: now the cluster list (cluster ids) can look like this: [0 1 5 7 3 6 8 8 3]
   // Afterwards, it should look like this:                                                                         [0 1 3 5 2 4 6 6 2]
   {
      unordered_map<uint32_t, uint32_t> ids; // (id, rank)
      {
         // get sorted cluster ids:
         vector<uint32_t> tmp_clusterIds; tmp_clusterIds.reserve(clusters.size());
         {
            // Used to make sure we don't re-insert any cluster ids twice: So after this tmp_clusterids should look like this: [0 1 5 7 3 6 8]
            unordered_set<uint32_t> tmp_clusterIds_set;
            
            for( uint i = 0; i < clusters.size(); ++i ) {
               const uint32_t id = *clusters[i].clusterId;
               //if( id == 0) {continue;}
               if( tmp_clusterIds_set.find( id ) == tmp_clusterIds_set.end() ) { tmp_clusterIds_set.insert(id); tmp_clusterIds.push_back(id); }
            }
         }
         // Sort:
         sort( tmp_clusterIds.begin(), tmp_clusterIds.end() ); //tmp_clusterids should look like this: [0 1 3 5 6 7 8]
         // Map the ids: after this, the ids should look like this: {0: 0, 1:1; 3:2, 5:3, 6:4, 7:5, 8:6}
         for( uint i = 0; i < tmp_clusterIds.size(); ++i ) {
            // Get id, and insert
            const uint32_t id = tmp_clusterIds[i];
            // Map the id:
            ids.insert( make_pair( id, i ) );
         }
      }

      // Make a temp copy of clusterids; clusterids_copy looks like [0 1 5 7 3 6 8 8 3]
      vector<uint32_t> clusterids_copy; clusterids_copy.resize(clusters.size());

      for( uint i = 0; i < clusters.size(); ++i ) {
         clusterids_copy[i] = *clusters[i].clusterId;
      }

      // Fix the clusterids according to the mapped values; afterwards clusters should look like this: [0 1 3 5 2 4 6 6 2]
      for( uint i = 0; i < clusters.size(); ++i ) {
         const uint32_t id = clusterids_copy[i];
         //if( id == 0) {continue;}
         *clusters[i].clusterId = ids.at(id);
      }

      // Input number of populations:
      cell->number_of_populations = ids.size();
   }

#ifndef NDEBUG
num_noclusters = 0;
for( vector<Cluster>::const_iterator it = clusters.begin(); it != clusters.end(); ++it ) {
   if( *(it->clusterId) == 0 ) { num_noclusters++; }
}
phiprof_assert( num_noclusters == 1 );
#endif


   // Set values of clusters and calculate number of clusters
   size_t notSet = 0;
   vmesh::VelocityBlockContainer<vmesh::LocalID> & blockContainer = cell->get_velocity_blocks(notSet);
   phiprof_assert( cell->get_number_of_velocity_blocks() * WID3 == velocityCells.size() );
   for( uint i = 0; i < velocityCells.size(); ++i ) {
      phiprof_assert( clusterIds[i] != 0 );
      const Realf value = *clusters[clusterIds[i]].clusterId;
      //cell->block_fx[i] = value-1;
   }


   // Print out the number of clusterIds:
   cout << "Clusters: " << clusterId << endl;
   cout << "Merges: " << merges << endl;
   cout << "Sizeof: " << sizeof(Realf) << endl;


   //delete[] clusterIds;
   //clusterIds = NULL;
   // Free memory from clusters:
   for( uint i = 0; i < clusters.size(); ++i ) {
      clusters[i].remove_data( clusters );
   }
   return;
}


namespace test {
   // A function for testing if VelocityCell class is functioning properly
   void testVelocityCells( SpatialCell * cell ) {
      size_t notSet = 0;
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> & vmesh = cell->get_velocity_mesh(notSet);
      vmesh::VelocityBlockContainer<vmesh::LocalID> & blockContainer = cell->get_velocity_blocks(notSet);
     
      vector<VelocityCell> velocityCells;
      velocityCells.resize(cell->get_number_of_velocity_blocks()*WID3);
      
      // Input data
      phiprof_assert( velocityCells.size() == cell->get_number_of_velocity_blocks() * VELOCITY_BLOCK_LENGTH );
      for( vmesh::LocalID i = 0; i < cell->get_number_of_velocity_blocks(); ++i ) {
         const vmesh::GlobalID blockId = cell->get_velocity_block_global_id(i);
         for( uint16_t vCellId = 0; vCellId < VELOCITY_BLOCK_LENGTH; ++vCellId ) {
            // Input the velocity cell into the vector
            velocityCells[i * VELOCITY_BLOCK_LENGTH + vCellId].set_data(vmesh, blockContainer, blockId, vCellId);
         }
      }
      
      
      Realf * data = blockContainer.getData();
      phiprof_assert( velocityCells.size() == blockContainer.size()*VELOCITY_BLOCK_LENGTH );
      for( int i = 0; i < cell->get_number_of_velocity_blocks(); ++i ) {
         for( uint16_t vCellId = 0; vCellId < VELOCITY_BLOCK_LENGTH; ++vCellId ) {
           VelocityCell & velocityCell = velocityCells[i*WID3 + vCellId];
           phiprof_assert(velocityCell.get_avgs() == data[i*WID3 + vCellId] );
         }
      }
   }
   
   void testIfSorted( vector<VelocityCell> velocityCells ) {
     if( velocityCells.size() == 0 ) {return;}
     phiprof_assert( velocityCells.size() != 0 );
     Realf last_value = -1e10;
     for( int i = 0; i < velocityCells.size(); ++i ) {
       VelocityCell & velocityCell = velocityCells[i];
       const Realf current_value = velocityCell.get_avgs();
       phiprof_assert( current_value >= last_value );
     }
   }
   
   // Function for testing that the hash works properly
   void testHash(  vector<VelocityCell> velocityCells,
                    SpatialCell * cell,
                    const Real resolution_threshold) {
     vmesh::VelocityBlockContainer< vmesh::LocalID > & blockContainer = cell->get_velocity_blocks(0);
     const int numberOfVCells = velocityCells.size();
     vector<Realf,aligned_allocator<Realf,64> > clusterIds; clusterIds.resize(numberOfVCells);
     for( int i = 0; i < numberOfVCells; ++i ) {
       clusterIds[i] = 0;
     }
     // Test hashing e.g. every id is unique:
     for( int i = 0; i < numberOfVCells; ++i ) {
       VelocityCell & vCell = velocityCells[i];
       const int id = vCell.hash();
       phiprof_assert( id >= 0 );
       phiprof_assert( id < numberOfVCells );
       phiprof_assert( clusterIds[id] == 0 );
       clusterIds[id] = 1;
     }
     return;
   }
   
   void getVelocityCellIndices(
                                 vmesh::LocalID vCell_indices[3],
                                 VelocityCell & vCell,
                                 SpatialCell * cell
                                 ) {
     // Get the velocity mesh:
#warning no AMR
     vmesh::VelocityMesh< vmesh::GlobalID, vmesh::LocalID > & vmesh = cell->get_velocity_mesh(0);
     // Get the block indices:
     vmesh::LocalID block_indices[3];
     uint8_t refinement_level = 0;
     vmesh.getIndices(vCell.block, refinement_level, block_indices[0], block_indices[1], block_indices[2]);
   
     // Get the local velocity cell indices:
     const uint16_t local_vCell_id = vCell.vCellId;
     const uint16_t local_vCell_indices[3] = { (uint16_t)(local_vCell_id%WID), (uint16_t)((local_vCell_id / WID)%WID), (uint16_t)(local_vCell_id / WID2) };
     
     // Input global velocity cell indices (unique indice for every velocity cell)
     for( int i = 0; i < 3; ++i ) {
       vCell_indices[i] = block_indices[i] * 4 + local_vCell_indices[i];
     }
   }
   
   uint64_t uniqueId( const vmesh::LocalID blockid, const uint16_t vCellid ) {
     return blockid * WID3 + vCellid;
   }
   
   // Function for testing that the fetching of neighbors works properly
   void testNeighbor( vector<VelocityCell> velocityCells,
                       SpatialCell * cell,
                       const Real resolution_threshold) {
     // the test should make sure that every velocity cell neighbor is
     // a) a distance of 1 away from the velocityCell
     // b) Unique (no duplicate neighbors)
     // c) Within the grid lengths
     const int numberOfVCells = velocityCells.size();
     for( int i = 0; i < numberOfVCells; ++i ) {
       VelocityCell & vCell = velocityCells[i];
       // Get the neighbors:
       vector<VelocityCell> neighbors; neighbors.reserve(30);
       getNeighbors(neighbors, vCell, cell); // Fetches neighbors
       phiprof_assert(neighbors.size() != 0 );
       
       // Loop through neighbors
       unordered_set<uint64_t> neighbor_duplicates_check;
       neighbor_duplicates_check.clear();
       for( int j = 0; j < neighbors.size(); ++j ) {
         VelocityCell & neighbor_vCell = neighbors[j];
         vmesh::LocalID vCell_indices[3], neighbor_vCell_indices[3];
         getVelocityCellIndices( vCell_indices, vCell, cell );
         getVelocityCellIndices( neighbor_vCell_indices, neighbor_vCell, cell );
 
         // Make sure the dist is one:
         int sum = 0;
         for( int dir = 0; dir < 3; ++dir ) {
           phiprof_assert( abs((int)vCell_indices[dir] - (int)neighbor_vCell_indices[dir]) <= 1 );
           sum += abs((int)vCell_indices[dir] - (int)neighbor_vCell_indices[dir]);
         }
         phiprof_assert( sum != 0 );
         
         // Make sure every neighbor is unique:
         //cout << "UNIQUE ID: " << uniqueId(neighbor_vCell.block, neighbor_vCell.vCellId) << endl;
         //for( unordered_set<uint64_t>::iterator it = neighbor_duplicates_check.begin(); it != neighbor_duplicates_check.end(); ++it ) {
         //   cout << "ALREADY HERE " <<  *it << endl;
         //}
         phiprof_assert( neighbor_duplicates_check.find( uniqueId(neighbor_vCell.block, neighbor_vCell.vCellId) ) == neighbor_duplicates_check.end() );
         neighbor_duplicates_check.insert( uniqueId(neighbor_vCell.block, neighbor_vCell.vCellId) );
       }
       
     
     }
     return;
   }
   
   void testCluster( vector<VelocityCell> velocityCells,
                      SpatialCell * cell,
                      const Real resolution_threshold ) {
     vector<Cluster> clusters;
     clusters.resize(4);
     // Create three clusters:

     uint32_t clusterId = nullCluster;
     clusters[clusterId].set_data( clusterId );
     clusterId++;
     clusters[clusterId].set_data( clusterId );
     clusterId++;
     clusters[clusterId].set_data( clusterId );
     clusterId++;
     clusters[clusterId].set_data( clusterId );
     // Make sure clusters 1 and 2 are not merged:
     phiprof_assert( !clusters[1].isNeighbor(*clusters[2].clusterId) );
     phiprof_assert( !clusters[2].isNeighbor(*clusters[1].clusterId) );
     // Merge cluster 1 and 2
     clusters[1].merge(clusters[2], clusters);
     //clusters[1].merge(clusters[0], clusters);
     // Check if cluster 1 and 2 are merged:
     phiprof_assert( clusters[1].isNeighbor(*clusters[2].clusterId) );
     phiprof_assert( clusters[2].isNeighbor(*clusters[1].clusterId) );
     // Make sure that cluster 1 and cluster 2 are not merged:
     phiprof_assert( !clusters[1].find(*clusters[2].clusterId) );
     // Merge cluster 2 and 3
     clusters[2].merge(clusters[3], clusters);
     // Now all 3 should be merged so check to make sure
     phiprof_assert(clusters[1].isNeighbor(*clusters[2].clusterId));
     phiprof_assert(clusters[1].isNeighbor(*clusters[3].clusterId));
     phiprof_assert(clusters[2].isNeighbor(*clusters[3].clusterId));
     // Deallocate data:
     clusters[0].remove_data(clusters);
     clusters[1].remove_data(clusters);
     clusters[2].remove_data(clusters);
     
     return;
   }
   
   void testNeighborSpeed(vector<VelocityCell> velocityCells,
                          SpatialCell * cell,
                          const Real resolution_threshold) {
     uint testValue = 0; vector <VelocityCell > neighbors;
     neighbors.reserve(30);
     for( int i = 0; i < velocityCells.size(); ++i ) {
       getNeighbors(neighbors, velocityCells[i], cell);
       int r = rand()%neighbors.size();
       VelocityCell neighborVelocityCell = neighbors[r];
       testValue += neighborVelocityCell.vCellId*i;
     }
     cell->number_of_populations = testValue;
   }
} // Namespace test

// Returns the total number of clusters:
int numberOfClusters( vector<Cluster> & clusters ) {
   if( clusters.size() == 0 ) { return 0; }
   unordered_set<uint32_t> uniqueClusterIds;
   for( int i = 1; i < clusters.size(); ++i ) {
      uniqueClusterIds.insert( *clusters[i].clusterId );
   }
   return uniqueClusterIds.size();
}

inline 
bool mergeClusters( Cluster & a, Cluster & b, VelocityCell & aVelocityCell, VelocityCell & bVelocityCell, const SpatialCell * cell, const Real minVolume, const Real minAvgs) {
  const Real * parameters = cell->parameters;
  return (min(*a.members, *b.members)*parameters[BlockParams::DVX]*parameters[BlockParams::DVX]*parameters[BlockParams::DVX] < minVolume) ||
         (min(aVelocityCell.get_avgs(), aVelocityCell.get_avgs()) > minAvgs );
}

void deallocateClusters( vector<Cluster> clusters ) {
   for( int i = 0; i < clusters.size(); ++i ) {
     clusters[i].remove_data( clusters );
   }
}

void saveClusterData( SpatialCell * cell, const vector<uint32_t> & velocityCellClusterId, const vector<Cluster> & clusters ) {
   vmesh::VelocityBlockContainer< vmesh::LocalID > & tmpBlockContainer = cell->get_velocity_blocks_temporary();
   tmpBlockContainer.clear();
   tmpBlockContainer.setSize(cell->get_number_of_velocity_blocks());
   Realf * tmpData = tmpBlockContainer.getData();

   for( int i = 0; i < velocityCellClusterId.size(); ++i ) {
      const Cluster & cluster = clusters[velocityCellClusterId[i]];

      tmpData[i] = (Realf)(*cluster.clusterId);
   }
}

void simpleCluster( vector<VelocityCell> & velocityCells, SpatialCell * cell, const Real totalVolumeThreshold, const Real minAvgs ) {
   if( velocityCells.size() == 0 ) { return; }
   phiprof_assert( velocityCells.size() > 0 );
   phiprof_assert( cell );
   
#warning no AMR
   const uint totalVolume = velocityCells.size();

   // For statistics:
   int merges = 0;
   int totalClusters = 0;

   // Create clusters:
   vector<Cluster> clusters;
   clusters.reserve(200);
   {
      // Create a null cluster
      Cluster newCluster;
      newCluster.set_data(nullCluster);
      clusters.push_back( newCluster );
   }
   uint32_t maxClusterIndex = nullCluster;

   // Create cluster ids for every velocity cell:
   vector<uint32_t> velocityCellClusterId;
   velocityCellClusterId.resize( velocityCells.size() );
   for( int i = 0; i < velocityCellClusterId.size(); ++i ) {
      velocityCellClusterId[i] = nullCluster;
   }

   // go through every velocity cell:
   for( int i = velocityCells.size()-1; i >= 0; --i ) {
      VelocityCell & velocityCell = velocityCells[i];
      const size_t id = velocityCell.hash();
      // Get the neighbors of the velocity cell:
      vector<VelocityCell> neighbors;
      neighbors.reserve(27);
      VelocityCell & vCell = velocityCells[i];
      getNeighbors(
                   neighbors,
                   vCell,
                   cell
                   );
      phiprof_assert( neighbors.size() != 0 );
      // go through the neighbors, check if any neighbor is not nullCluster:
      for( int k = 0; k < neighbors.size(); ++k ) {
         VelocityCell & velocityCellNeighbor = neighbors[k];
         const size_t neighborId = velocityCellNeighbor.hash();
         const uint32_t neighborClusterIndex = velocityCellClusterId[neighborId];
         // Check if the neighbor belongs to any cluster
         if( neighborClusterIndex != nullCluster ) {
            uint32_t & clusterIndex = velocityCellClusterId[id];
            // If the velocity cell is not already a part of a cluster, put it in the neighbor:
            if( clusterIndex == nullCluster ) {
               clusterIndex = neighborClusterIndex;
               phiprof_assert( clusters.size() > clusterIndex );
               clusters[clusterIndex].append(1); // Increase number of members by one
            }
            // If the velocity cell is already a part of a cluster (and they are not the same cluster), merge the cluster and the neighbor cluster:
            else if (clusterIndex != nullCluster && clusters[clusterIndex].clusterId != clusters[neighborClusterIndex].clusterId ) {
               phiprof_assert( clusters.size() > clusterIndex && clusters.size() > neighborClusterIndex );
               const Real avgsThreshold = 1.0e-15;
               const Real minVolume = 3000*cell->parameters[BlockParams::DVX]*cell->parameters[BlockParams::DVY]*cell->parameters[BlockParams::DVZ];
               if( mergeClusters(clusters[clusterIndex], clusters[neighborClusterIndex], velocityCell, velocityCellNeighbor, cell, minVolume, avgsThreshold) ) {
                 clusters[clusterIndex].merge( clusters[neighborClusterIndex], clusters );
                 ++merges;
               }
            }
         }
      }

      // If this velocity cell was not put into any cluster, then create a new cluster for it:
      uint32_t & clusterIndex = velocityCellClusterId[id];
      if( clusterIndex == nullCluster ) {
         maxClusterIndex++;
         clusterIndex = maxClusterIndex;
         Cluster newCluster;
         newCluster.set_data( maxClusterIndex );
         clusters.push_back( newCluster );
      }
   }
   // Print the number of clusters:
   cout << "Number of merges: " << merges << endl;
   cout << "Number of clusters: " << numberOfClusters( clusters ) << endl;
   cell->number_of_populations = numberOfClusters( clusters );
   
   saveClusterData(cell, velocityCellClusterId, clusters);
   
   deallocateClusters(clusters);
}



/*! Function for calculating the different populations in distribution for a given spatial cell. Note that the results are currently saved in block_fx array within the SpatialCell.

 \param cell                                A cell from the DCCRG grid
 \param resolution_threshold                A value for determining how large at minimum we want our populations to be. 0.006 seems ok unless there's a reason to believe otherwise.

 */
void populationAlgorithm( 
                SpatialCell * cell,
                const Real resolution_threshold
                   ) {
   // Get velocity mesh:
#warning vmesh: popId not set in get_velocity_mesh
   const size_t notSet = 0;
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> & vmesh = cell->get_velocity_mesh(notSet);
   vmesh::VelocityBlockContainer<vmesh::LocalID> & blockContainer = cell->get_velocity_blocks(notSet);

   
   // Vector for holding velocity cells:
   vector<VelocityCell> velocityCells;
   // Initialize avgs values vector:
   velocityCells.resize( cell->get_number_of_velocity_blocks() * VELOCITY_BLOCK_LENGTH );

   // Input data
   for( vmesh::LocalID i = 0; i < cell->get_number_of_velocity_blocks(); ++i ) {
      // Create a new velocity cell
      VelocityCell input_cell;
      const vmesh::GlobalID blockId = cell->get_velocity_block_global_id(i);
      for( uint16_t vCellId = 0; vCellId < VELOCITY_BLOCK_LENGTH; ++vCellId ) {
         // Input the block data
         input_cell.set_data( vmesh, blockContainer, blockId, vCellId);
         // Input the velocity cell into the vector
         velocityCells[i * VELOCITY_BLOCK_LENGTH + vCellId] = input_cell;
      }
   }

   phiprof::start("sort_velocity_space");
   // Sort the list:
   sort(velocityCells.begin(), velocityCells.end());
   phiprof::stop("sort_velocity_space");

   // Retrieve populations and save into block_fx array inside spatial cell:
   //cluster_advanced( velocityCells, cell, resolution_threshold );
   //testNeighbor_speed( velocityCells, local_vcell_neighbors, remote_vcell_neighbors, cell, resolution_threshold );

#ifndef NDEBUG
   // Check the functionality of the different algorithm components:
   test::testHash(velocityCells, cell, resolution_threshold);
   test::testNeighbor(velocityCells, cell, resolution_threshold);
   test::testCluster(velocityCells, cell, resolution_threshold);
   test::testVelocityCells(cell);
   test::testIfSorted(velocityCells);
#endif

   phiprof::start("getNeighbors");
   test::testNeighborSpeed( velocityCells, cell, resolution_threshold);
   phiprof::stop("getNeighbors");
   
   // Do clustering:
   const Real minAvgs = 1.0e-10;
   simpleCluster( velocityCells, cell, resolution_threshold, minAvgs );
   //connectivityCluster(velocityCells, cell, resolution_threshold);


   return;
}