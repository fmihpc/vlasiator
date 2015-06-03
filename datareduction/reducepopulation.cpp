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

         // Make sure this class has no null pointers:
         phiprof_assert( neighbor_clusters != cluster_neighbor.neighbor_clusters );
         phiprof_assert( neighbor_clusters && members && clusterId );

         phiprof_assert( cluster_neighbor.neighbor_clusters && cluster_neighbor.members && cluster_neighbor.clusterId );

         phiprof_assert( cluster_neighbor.clusterId != clusterId );

         phiprof_assert( clusterId != NULL );

         phiprof_assert( cluster_neighbor.clusterId != NULL );

         phiprof_assert( *clusterId != 0 );

         phiprof_assert( *cluster_neighbor.clusterId != 0 );

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
     uint32_t nullCluster = 0;
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
}

// Returns the total number of clusters:
int numberOfClusters( vector<Cluster> & clusters ) {
   if( clusters.size() == 0 ) { return 0; }
   unordered_set<uint32_t> uniqueClusterIds;
   for( int i = 1; i < clusters.size(); ++i ) {
      uniqueClusterIds.insert( *clusters[i].clusterId );
   }
   return uniqueClusterIds.size();
}

void simpleCluster( vector<VelocityCell> & velocityCells, SpatialCell * cell, const Real resolution_threshold ) {
   phiprof_assert( velocityCells.size() > 0 );
   phiprof_assert( cell );

   // For statistics:
   int merges = 0;

   // Create clusters:
   vector<Cluster> clusters;
   clusters.reserve(200);
   uint32_t nullCluster = 0;
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
   for( int i = 0; i < velocityCells.size(); ++i ) {
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
               if( min(*clusters[clusterIndex].members, *clusters[neighborClusterIndex].members) < resolution_threshold*velocityCells.size() ) {
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
   //cout << "Number of merges: " << merges << endl;
   //cout << "Number of clusters: " << numberOfClusters( clusters ) << endl;
   cell->number_of_populations = numberOfClusters( clusters );
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
#endif
   
   // Do clustering:
   simpleCluster( velocityCells, cell, resolution_threshold );


   return;
}

/*! Function for getting the max number of populations Note: the populations must have been saved

 \param mpiGrid                 The DCCRG grid with spatial cells
 */
//static uint max_number_of_populations( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) {
//   uint local_max_populations = 0;
//   // Fetch cells:
//   vector<uint64_t> local_cells = mpiGrid.get_cells();
//   // Loop through cells:
//   for( vector<uint64_t>::const_iterator it = local_cells.begin(); it != local_cells.end(); ++it ) {
//      // Get cell:
//      const uint64_t cellid = *it;
//      // Get the cell's spatial cell class:
//      SpatialCell * cell = mpiGrid[cellid];
//      // Check number of populations:
//      if( local_max_populations < cell->number_of_populations ) { local_max_populations = cell->number_of_populations; }
//   }
//   // Get the absolute maximum population out of all the processes:
//   uint max_populations = 0;
//   {
//      const int count = 1;
//      MPI_Allreduce( &local_max_populations, &max_populations, count, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD );
//   }
//   return max_populations;
//}

/*! Function for writing out rho for different populations:

 \param max_populations         Max number of populations in all cells
 \param local_cells             Cells for which to write the population
 \param mpiGrid                 The DCCRG grid with spatial cells
 \param vlsvWriter              The VLSV writer class for writing VLSV files, note that the file must have been opened already

 */
//static inline bool write_rho( const uint max_populations, const vector<uint64_t> & local_cells, const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, Writer & vlsvWriter ) {
//   bool success = true;
//   // Write out the population:
//   //success = vlsvWriter.writeArray("VARIABLE", mathematicsxmlAttributes, arraySize, 1, population_rho.data());
//   if( local_cells.size() == 0 ) {
//      const uint64_t arraySize = 0;
//      const uint64_t vectorSize = max_populations; // Population is Real, so a scalar
//      const string name = "PopulationRho";
//      map<string, string> xmlAttributes;
//      xmlAttributes["name"] = name;
//      xmlAttributes["mesh"] = "SpatialGrid";
//      Real dummy_data = 0;
//      return vlsvWriter.writeArray( "VARIABLE", xmlAttributes, arraySize, vectorSize, &dummy_data );
//   }
//
//
//
//   vector<Real> population_rho;
//   population_rho.resize( max_populations * local_cells.size() );
//
//   for(uint i = 0; i < population_rho.size(); ++i ) {
//      population_rho[i] = 0;
//   }
//
//   // Fetch different population_rho:
//   for( unsigned int i = 0; i < local_cells.size(); ++i ) {
//      const uint64_t cellId = local_cells[i];
//      SpatialCell * cell = mpiGrid[cellId];
//      if( cell->get_number_of_velocity_blocks() > 0 ) {
//         const Velocity_Block* block_ptr = cell->at(cell->velocity_block_list[0]);
//         const Real DV3 = block_ptr->parameters[BlockParams::DVX] * block_ptr->parameters[BlockParams::DVY] * block_ptr->parameters[BlockParams::DVZ];
//         vector<Real> thread_population_rho;
//         thread_population_rho.resize( max_populations );
//         // Set to zero:
//         for( uint j = 0; j < max_populations; ++j ) { thread_population_rho[j] = 0; }
//         #pragma omp parallel
//         {
//            vector<Real> thread_population_rho;
//            thread_population_rho.resize( max_populations );
//            // Set to zero:
//            for( uint j = 0; j < max_populations; ++j ) { thread_population_rho[j] = 0; }
//            // Loop through velocity cells:
//            #pragma omp for
//            for( uint v = 0; v < cell->number_of_blocks * WID3; ++v ) {
//               phiprof_assert( max_populations*i + (uint)(cell->block_fx[v]) < max_populations * local_cells.size() );
//               phiprof_assert( max_populations > (uint)(cell->block_fx[v]) );
//               thread_population_rho[(uint)(cell->block_fx[v])] += cell->block_data[v] * DV3;
//            }
//            const uint cell_index = max_populations*i;
//            for( uint j = 0; j < max_populations; ++j ) {
//               #pragma omp atomic
//               population_rho[cell_index + j] += thread_population_rho[j];
//            }
//         }
//      }
//   }
//
//   // Write the data out, we need arraySizze, vectorSize and name to do this
//   const uint64_t arraySize = local_cells.size();
//   const uint64_t vectorSize = max_populations; // Population is Real, so a scalar
//   const string name = "PopulationRho";
//
//   map<string, string> xmlAttributes;
//   xmlAttributes["name"] = name;
//   xmlAttributes["mesh"] = "SpatialGrid";
//
//   // Write the array and return false if the writing fails
//
//   if( vlsvWriter.writeArray( "VARIABLE", xmlAttributes, arraySize, vectorSize, population_rho.data() ) == false ) {
//      success = false;
//   }
//   return success;
//}



/*! Function for writing out rho_v for different populations:

 \param max_populations         Max number of populations in all cells
 \param local_cells             Cells for which to write the population
 \param mpiGrid                 The DCCRG grid with spatial cells
 \param vlsvWriter              The VLSV writer class for writing VLSV files, note that the file must have been opened already

 */
//static inline bool write_rho_v( const uint max_populations, const vector<uint64_t> & local_cells, const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, Writer & vlsvWriter ) {
//   bool success = true;
//
//   if( local_cells.size() == 0 ) {
//      const uint64_t arraySize = 0;
//      const uint64_t vectorSize = max_populations; // Population is Real, so a scalar
//      Real dummy_data = 0;
//
//
//      map<string, string> xmlAttributes;
//      xmlAttributes["name"] = "PopulationRhoVX";
//      xmlAttributes["mesh"] = "SpatialGrid";
//      // Write the array and return false if the writing fails
//      if( vlsvWriter.writeArray( "VARIABLE", xmlAttributes, arraySize, vectorSize, &dummy_data ) == false ) {
//         success = false;
//      }
//   
//      xmlAttributes.clear();
//      xmlAttributes["name"] = "PopulationRhoVY";
//      xmlAttributes["mesh"] = "SpatialGrid";
//      // Write the array and return false if the writing fails
//      if( vlsvWriter.writeArray( "VARIABLE", xmlAttributes, arraySize, vectorSize, &dummy_data ) == false ) {
//         success = false;
//      }
//   
//      xmlAttributes.clear();
//      xmlAttributes["name"] = "PopulationRhoVZ";
//      xmlAttributes["mesh"] = "SpatialGrid";
//        // Write the array and return false if the writing fails
//      if( vlsvWriter.writeArray( "VARIABLE", xmlAttributes, arraySize, vectorSize, &dummy_data ) == false ) {
//         success = false;
//      }
//      return success;
//   }
//
//   // Write out the population:
//   //success = vlsvWriter.writeArray("VARIABLE", xmlAttributes, arraySize, 1, population_rho.data());
//   vector<Real> population_rho_v_x;
//   population_rho_v_x.resize( max_populations * local_cells.size() );
//   vector<Real> population_rho_v_y;
//   population_rho_v_y.resize( max_populations * local_cells.size() );
//   vector<Real> population_rho_v_z;
//   population_rho_v_z.resize( max_populations * local_cells.size() );
//
//   const Real HALF = 0.5;
//
//   for(uint i = 0; i < population_rho_v_x.size(); ++i ) {
//      population_rho_v_x[i] = 0;
//      population_rho_v_y[i] = 0;
//      population_rho_v_z[i] = 0;
//   }
//
//   // Fetch different population_rho:
//   for( unsigned int i = 0; i < local_cells.size(); ++i ) {
//      const uint64_t cellId = local_cells[i];
//      SpatialCell * cell = mpiGrid[cellId];
//      const Velocity_Block* block_ptr = cell->at(cell->velocity_block_list[0]);
//      const Real DV3 = block_ptr->parameters[BlockParams::DVX] * block_ptr->parameters[BlockParams::DVY] * block_ptr->parameters[BlockParams::DVZ];
//      // Loop through blocks:
//      #pragma omp parallel
//      {
//         vector<Real> thread_population_rho_v_x;
//         thread_population_rho_v_x.resize( max_populations * local_cells.size() );
//         vector<Real> thread_population_rho_v_y;
//         thread_population_rho_v_y.resize( max_populations * local_cells.size() );
//         vector<Real> thread_population_rho_v_z;
//         thread_population_rho_v_z.resize( max_populations * local_cells.size() );
//
//         for( uint j = 0; j < max_populations; ++j ) {
//            thread_population_rho_v_x[j] = 0;
//            thread_population_rho_v_y[j] = 0;
//            thread_population_rho_v_z[j] = 0;
//         }
//
//         #pragma omp for
//         for( uint b = 0; b < cell->number_of_blocks; ++b ) {
//            const unsigned int blockId = cell->velocity_block_list[b];
//            const Velocity_Block* block = cell->at( blockId );
//            for( uint vi = 0; vi < WID; ++vi ) for( uint vj = 0; vj < WID; ++vj ) for( uint vk = 0; vk < WID; ++vk ) {
//               const uint vCellId = vi + vj*WID + vk*WID*WID;
//               const uint v = b*WID3 + vCellId;
//               phiprof_assert( max_populations*i + (uint)(cell->block_fx[v]) < max_populations * local_cells.size() );
//               phiprof_assert( max_populations > (uint)(cell->block_fx[v]) );
//               const Real VX = block->parameters[BlockParams::VXCRD] + (vi+HALF) * block->parameters[BlockParams::DVX];
//               const Real VY = block->parameters[BlockParams::VYCRD] + (vj+HALF) * block->parameters[BlockParams::DVY];
//               const Real VZ = block->parameters[BlockParams::VZCRD] + (vk+HALF) * block->parameters[BlockParams::DVZ];
//               // Input values into different populations
//               thread_population_rho_v_x[(uint)(cell->block_fx[v])] += cell->block_data[v] * DV3 * VX;
//               thread_population_rho_v_y[(uint)(cell->block_fx[v])] += cell->block_data[v] * DV3 * VY;
//               thread_population_rho_v_z[(uint)(cell->block_fx[v])] += cell->block_data[v] * DV3 * VZ;
//            }
//         }
//         #pragma omp critical
//         {
//            for( uint j = 0; j < max_populations; ++j ) {
//               population_rho_v_x[max_populations*i + j] += thread_population_rho_v_x[j];
//               population_rho_v_y[max_populations*i + j] += thread_population_rho_v_y[j];
//               population_rho_v_z[max_populations*i + j] += thread_population_rho_v_z[j];
//            }
//         }
//      }
//   }
//
//   // Write the data out, we need arraySizze, vectorSize and name to do this
//   const uint64_t arraySize = local_cells.size();
//   const uint64_t vectorSize = max_populations; // Population is Real, so a scalar
//
//   map<string, string> xmlAttributes;
//   xmlAttributes["name"] = "PopulationRhoVX";
//   xmlAttributes["mesh"] = "SpatialGrid";
//
//   // Write the array and return false if the writing fails
//   if( vlsvWriter.writeArray( "VARIABLE", xmlAttributes, arraySize, vectorSize, population_rho_v_x.data() ) == false ) {
//      success = false;
//   }
//
//   xmlAttributes.clear();
//   xmlAttributes["name"] = "PopulationRhoVY";
//   xmlAttributes["mesh"] = "SpatialGrid";
//
//   // Write the array and return false if the writing fails
//   if( vlsvWriter.writeArray( "VARIABLE", xmlAttributes, arraySize, vectorSize, population_rho_v_y.data() ) == false ) {
//      success = false;
//   }
//
//   xmlAttributes.clear();
//   xmlAttributes["name"] = "PopulationRhoVZ";
//   xmlAttributes["mesh"] = "SpatialGrid";
//
//   // Write the array and return false if the writing fails
//   if( vlsvWriter.writeArray( "VARIABLE", xmlAttributes, arraySize, vectorSize, population_rho_v_z.data() ) == false ) {
//      success = false;
//   }
//
//   return success;
//}


/*! Writes the population variables, so variables for different populations.

 \param mpiGrid                 The DCCRG grid with spatial cells
 \param vlsvWriter              The VLSV writer class for writing VLSV files, note that the file must have been opened already
 */
//bool write_population_variables( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, Writer & vlsvWriter ) {
//   // Get max number of populations (to be used in determining how many arrays to write)
//   const uint max_populations = max_number_of_populations( mpiGrid );
////   // Fetch cells:
////   const vector<uint64_t> local_cells = mpiGrid.get_cells();
//   // Write rho:
//   vector<uint64_t> local_cells = mpiGrid.get_cells();
//   bool success = write_rho( max_populations, local_cells, mpiGrid, vlsvWriter );
//   if( success == true ) { success == write_rho_v( max_populations, local_cells, mpiGrid, vlsvWriter ); }
//   return success;
//}




/*! Writes the distribution function for different populations. Note that population_algorithm must have been called before this!
    
 \param mpiGrid                 The DCCRG grid with spatial cells
 \param vlsvWriter              The VLSV writer class for writing VLSV files, note that the file must have been opened already
 \param local_cells             The cells for which we write out the velocity space
 */ 
//bool write_population_distribution( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const vector<uint64_t> & local_cells, vlsv::Writer & vlsvWriter ) {
//   string name = "BLOCKVARIABLEPOPULATION";
//
//   //Compute totalBlocks
//   uint64_t totalBlocks = 0;  
//   for(size_t cell=0;cell<local_cells.size();++cell){
//      totalBlocks+=mpiGrid[local_cells[cell]]->number_of_blocks;
//   }
//
//   const string datatype_avgs = "float";
//   const uint64_t arraySize_avgs = totalBlocks;
//   const uint64_t vectorSize_avgs = WID3; // There are 64 elements in every velocity block
//   const uint64_t dataSize_avgs = sizeof(Realf);
//   map<string,string> attribs;
//   attribs["mesh"] = "SpatialGrid";
//   attribs["name"] = "avgs"; // Name of the velocity space distribution is written avgs
//
//   // Start multi write
//   vlsvWriter.startMultiwrite(datatype_avgs,arraySize_avgs,vectorSize_avgs,dataSize_avgs);
//   // Loop over local_cells
//   for( size_t cell = 0; cell < local_cells.size(); ++cell ) {
//      // Get the spatial cell
//      SpatialCell* SC = mpiGrid[local_cells[cell]];
//      // Get the number of blocks in this cell
//      const uint64_t arrayElements = SC->get_number_of_velocity_blocks();
//      char * arrayToWrite = reinterpret_cast<char*>(SC->get_data());
//      // Add a subarray to write
//      vlsvWriter.addMultiwriteUnit(arrayToWrite, arrayElements); // Note: We told beforehands that the vectorsize = WID3 = 64
//   }
//   if( local_cells.size() == 0 ) {
//      // Write dummy array to avoid MPI blocking
//      vlsvWriter.addMultiwriteUnit(NULL, 0);
//   }
//   // Write the subarrays
//   vlsvWriter.endMultiwrite(name, attribs);
//   return true;
//}


