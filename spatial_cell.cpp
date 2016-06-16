/*!
Space for static variables of spatial cell class for Vlasiator.

Copyright 2011-2015 Finnish Meteorological Institute


*/

#include <unordered_set>
#include <vectorclass.h>

#include "spatial_cell.hpp"
#include "velocity_blocks.h"
#include "object_wrapper.h"

#ifndef NDEBUG
   #define DEBUG_SPATIAL_CELL
#endif

using namespace std;

namespace spatial_cell {
   int SpatialCell::activePopID = -1;
   uint64_t SpatialCell::mpi_transfer_type = 0;
   bool SpatialCell::mpiTransferAtSysBoundaries = false;

   SpatialCell::SpatialCell() {
      // Block list and cache always have room for all blocks
      this->sysBoundaryLayer=0; // Default value, layer not yet initialized
      for (unsigned int i=0; i<WID3; ++i) null_block_data[i] = 0.0;

      // reset spatial cell parameters
      for (unsigned int i = 0; i < CellParams::N_SPATIAL_CELL_PARAMS; i++) {
         this->parameters[i]=0.0;
      }
      
      // reset spatial cell derivatives
      for (unsigned int i = 0; i < fieldsolver::N_SPATIAL_CELL_DERIVATIVES; i++) {
         this->derivatives[i]=0;
      }
      
      // reset BVOL derivatives
      for (unsigned int i = 0; i < bvolderivatives::N_BVOL_DERIVATIVES; i++) {
         this->derivativesBVOL[i]=0;
      }
      //is transferred by default
      this->mpiTransferEnabled=true;
      
      // Set correct number of populations
      populations.resize(getObjectWrapper().particleSpecies.size());
      
      // Set velocity meshes
      for (int popID=0; popID<populations.size(); ++popID) {
         const species::Species& spec = getObjectWrapper().particleSpecies[popID];
         populations[popID].vmesh.initialize(spec.velocityMesh);
         populations[popID].velocityBlockMinValue = spec.sparseMinValue;
      }
   }

   SpatialCell::SpatialCell(const SpatialCell& other):
     initialized(other.initialized),
     mpiTransferEnabled(other.mpiTransferEnabled),
     velocity_block_with_content_list(other.velocity_block_with_content_list),
     velocity_block_with_no_content_list(other.velocity_block_with_no_content_list),
     sysBoundaryFlag(other.sysBoundaryFlag),
     sysBoundaryLayer(other.sysBoundaryLayer),
     sysBoundaryLayerNew(other.sysBoundaryLayerNew),
     populations(other.populations) {

        //copy parameters
        for(unsigned int i=0;i< CellParams::N_SPATIAL_CELL_PARAMS;i++){
           parameters[i]=other.parameters[i];
        }
        //copy derivatives
        for(unsigned int i=0;i< fieldsolver::N_SPATIAL_CELL_DERIVATIVES;i++){
           derivatives[i]=other.derivatives[i];
        }
        //copy BVOL derivatives
        for(unsigned int i=0;i< bvolderivatives::N_BVOL_DERIVATIVES;i++){
           derivativesBVOL[i]=other.derivativesBVOL[i];
        }
        
        //set null block data
        for (unsigned int i=0; i<WID3; ++i) null_block_data[i] = 0.0;
   }


   /** Adds "important" and removes "unimportant" velocity blocks
    * to/from this cell.
    * 
    * velocity_block_with_content_list needs to be up to date in local and remote cells.
    * velocity_block_with_no_content_list needs to be up to date in local cells.
    *         
    * update_velocity_block_with_content_lists() should have
    * been called with the current distribution function values, and then the contetn list transferred.
    * 
    * Removes all velocity blocks from this spatial cell which don't
    * have content and don't have spatial or velocity neighbors with
    * content.  Adds neighbors for all velocity blocks which do have
    * content (including spatial neighbors).  All cells in
    * spatial_neighbors are assumed to be neighbors of this cell.
    * 
    * This function is thread-safe when called for different cells
    * per thread. We need the block_has_content vector from
    * neighbouring cells, but these are not written to here. We only
    * modify local cell.
    * 
    * NOTE: The AMR mesh must be valid, otherwise this function will
    * remove some blocks that should not be removed.*/
   #ifndef AMR
   void SpatialCell::adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors,
                                            const int& popID,bool doDeleteEmptyBlocks) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      //  This set contains all those cellids which have neighbors in any
      //  of the 6-dimensions Actually, we would only need to add
      //  local blocks with no content here, as blocks with content
      //  do not need to be created and also will not be removed as
      //  we only check for removal for blocks with no content
      std::unordered_set<vmesh::GlobalID> neighbors_have_content;

      //add neighbor content info for velocity space neighbors to map. We loop over blocks
      //with content and raise the neighbors_have_content for
      //itself, and for all its neighbors
      for (vmesh::LocalID block_index=0; block_index<velocity_block_with_content_list.size(); ++block_index) {
         vmesh::GlobalID block = velocity_block_with_content_list[block_index];

         const uint8_t refLevel=0;
         const velocity_block_indices_t indices = SpatialCell::get_velocity_block_indices(popID,block);
         neighbors_have_content.insert(block); //also add the cell itself
         
         for (int offset_vx=-P::sparseBlockAddWidthV;offset_vx<=P::sparseBlockAddWidthV;offset_vx++) {
            for (int offset_vy=-P::sparseBlockAddWidthV;offset_vy<=P::sparseBlockAddWidthV;offset_vy++) {
               for (int offset_vz=-P::sparseBlockAddWidthV;offset_vz<=P::sparseBlockAddWidthV;offset_vz++) {
                  const vmesh::GlobalID neighbor_block 
                     = get_velocity_block(popID,{{indices[0]+offset_vx,indices[1]+offset_vy,indices[2]+offset_vz}},refLevel);
                  neighbors_have_content.insert(neighbor_block); //add all potential ngbrs of this block with content
               }
            }
         }
      }

      //add neighbor content info for spatial space neighbors to map. We loop over
      //neighbor cell lists with existing blocks, and raise the
      //flag for the local block with same block id
      for (std::vector<SpatialCell*>::const_iterator neighbor=spatial_neighbors.begin();
           neighbor != spatial_neighbors.end(); ++neighbor) {
         for (vmesh::LocalID block_index=0; block_index<(*neighbor)->velocity_block_with_content_list.size(); ++block_index) {
            vmesh::GlobalID block = (*neighbor)->velocity_block_with_content_list[block_index];
            neighbors_have_content.insert(block);
         }
      }

      // REMOVE all blocks in this cell without content + without neighbors with content
      // better to do it in the reverse order, as then blocks at the
      // end are removed first, and we may avoid copying extra data.
      if (doDeleteEmptyBlocks) {
         for (int block_index= this->velocity_block_with_no_content_list.size()-1; block_index>=0; --block_index) {
            const vmesh::GlobalID blockGID = velocity_block_with_no_content_list[block_index];
            #ifdef DEBUG_SPATIAL_CELL
            if (blockGID == invalid_global_id()) {
               cerr << "Got invalid block at " << __FILE__ << ' ' << __LINE__ << endl; 
               exit(1);
            }             
            #endif
            const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID,popID);            
            #ifdef DEBUG_SPATIAL_CELL
            if (blockLID == invalid_local_id()) {
               cerr << "Could not find block in " << __FILE__ << ' ' << __LINE__ << endl; 
               exit(1);
            }
            #endif
            
            bool removeBlock = false;
            std::unordered_set<vmesh::GlobalID>::iterator it = neighbors_have_content.find(blockGID);
            if (it == neighbors_have_content.end()) removeBlock = true;

            if (removeBlock == true) {
               //No content, and also no neighbor have content -> remove
               //and increment rho loss counters
               const Real* block_parameters = get_block_parameters(popID)+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
               const Real DV3 = block_parameters[BlockParams::DVX]
                 * block_parameters[BlockParams::DVY]
                 * block_parameters[BlockParams::DVZ];
               Real sum=0;
               for (unsigned int i=0; i<WID3; ++i) sum += get_data(popID)[blockLID*SIZE_VELBLOCK+i];
               this->parameters[CellParams::RHOLOSSADJUST] += DV3*sum;
	       
               // and finally remove block
               this->remove_velocity_block(blockGID,popID);
            }
         }
      }

      // ADD all blocks with neighbors in spatial or velocity space (if it exists then the block is unchanged)
      for (std::unordered_set<vmesh::GlobalID>::iterator it=neighbors_have_content.begin(); it != neighbors_have_content.end(); ++it) {
         this->add_velocity_block(*it,popID);
      }
   }

   #else       // AMR version

   void SpatialCell::adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors,
                                            const int& popID,bool doDeleteEmptyBlocks) {
      //  This set contains all those cell ids which have neighbors in any
      //  of the 6-dimensions Actually, we would only need to add
      //  local blocks with no content here, as blocks with content
      //  do not need to be created and also will not be removed as
      //  we only check for removal for blocks with no content
      std::unordered_set<vmesh::GlobalID> neighbors_have_content;

      for (vmesh::LocalID block_index=0; block_index<velocity_block_with_content_list.size(); ++block_index) {
         vmesh::GlobalID blockGID = velocity_block_with_content_list[block_index];
         vector<vmesh::GlobalID> neighborGIDs;
         populations[popID].vmesh.getNeighborsExistingAtSameLevel(blockGID,neighborGIDs);
         neighbors_have_content.insert(neighborGIDs.begin(),neighborGIDs.end());
         neighbors_have_content.insert(blockGID);
      }

      //add neighbor content info for spatial space neighbors to map. We loop over
      //neighbor cell lists with existing blocks, and raise the
      //flag for the local block with same block id
      unordered_set<vmesh::GlobalID> spat_nbr_has_content;
      for (std::vector<SpatialCell*>::const_iterator neighbor=spatial_neighbors.begin();
           neighbor != spatial_neighbors.end(); ++neighbor) {
         for (vmesh::LocalID block_index=0; block_index<(*neighbor)->velocity_block_with_content_list.size(); ++block_index) {
            vmesh::GlobalID block = (*neighbor)->velocity_block_with_content_list[block_index];
            spat_nbr_has_content.insert(block);
         }
      }      

      // REMOVE all blocks in this cell without content + without neighbors with content
      // better to do it in the reverse order, as then blocks at the
      // end are removed first, and we may avoid copying extra data.
      if (doDeleteEmptyBlocks) {
         for (int block_index= this->velocity_block_with_no_content_list.size()-1; block_index>=0; --block_index) {
            const vmesh::GlobalID blockGID = this->velocity_block_with_no_content_list[block_index];
            #ifdef DEBUG_SPATIAL_CELL
               if (blockGID == invalid_global_id()) {
                  cerr << "Got invalid block at " << __FILE__ << ' ' << __LINE__ << endl;
                  exit(1); 
               }
            #endif
            const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID,popID);
            #ifdef DEBUG_SPATIAL_CELL
               if (blockLID == invalid_local_id()) {
                  cerr << "Could not find block in " << __FILE__ << ' ' << __LINE__ << endl;
                  exit(1);
               }
            #endif
            
            bool removeBlock = false;

            // Check this block in the neighbor cells
            if (neighbors_have_content.find(blockGID) != neighbors_have_content.end()) continue;
            
            // Check the parent of this block in the neighbor cells
            if (neighbors_have_content.find(populations[popID].vmesh.getParent(blockGID)) != neighbors_have_content.end()) continue;
            
            // Check all the children of this block in the neighbor cells
            std::vector<vmesh::GlobalID> children;
            populations[popID].vmesh.getChildren(blockGID,children);
            int counter = 0;
            for (size_t c=0; c<children.size(); ++c) {
               if (neighbors_have_content.find(children[c]) != neighbors_have_content.end()) ++counter;
            }
            if (counter > 0) continue;
            
            // It is safe to remove this block
            removeBlock = true;

            if (removeBlock == true) {
               //No content, and also no neighbor have content -> remove
               //and increment rho loss counters
               const Real* block_parameters = get_block_parameters(popID)+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS;
               const Real DV3 = block_parameters[BlockParams::DVX]
                 * block_parameters[BlockParams::DVY]
                 * block_parameters[BlockParams::DVZ];
               Real sum=0;
               for (unsigned int i=0; i<WID3; ++i) sum += get_data(popID)[blockLID*SIZE_VELBLOCK+i];
               this->parameters[CellParams::RHOLOSSADJUST] += DV3*sum;
	       
               // and finally remove block
               this->remove_velocity_block(blockGID,popID);
            }
         }
      }

      // Filter the spat_nbr_has_content list so that it doesn't
      // contain overlapping blocks
      unordered_set<vmesh::GlobalID> ghostBlockList;
      vector<vector<vmesh::GlobalID> > sorted(populations[popID].vmesh.getMaxAllowedRefinementLevel()+1);

      // First sort the list according to refinement levels. This allows
      // us to skip checking the existence of children and grandchildren below.
      for (unordered_set<vmesh::GlobalID>::const_iterator it=spat_nbr_has_content.begin(); 
           it!=spat_nbr_has_content.end(); ++it) {
         sorted[populations[popID].vmesh.getRefinementLevel(*it)].push_back(*it);         
      }

      // Iterate through the sorted block list, and determine the no-content
      // blocks and their correct refinement levels that must exist
      for (int r=sorted.size()-1; r >= 0; --r) {
         for (size_t b=0; b<sorted[r].size(); ++b) {
            // If parent exists, all siblings must exist
            if (r > 0) {
               if (spat_nbr_has_content.find(populations[popID].vmesh.getParent(sorted[r][b])) != spat_nbr_has_content.end()) {
                  vector<vmesh::GlobalID> siblings;
                  populations[popID].vmesh.getSiblings(sorted[r][b],siblings);
                  ghostBlockList.insert(siblings.begin(),siblings.end());
                  continue;
               }
            }

            // If grandparent exists, parent octant must exist
            if (r > 1) {
               vmesh::GlobalID grandParentGID = populations[popID].vmesh.getParent(populations[popID].vmesh.getParent(sorted[r][b]));
               if (spat_nbr_has_content.find(grandParentGID) != spat_nbr_has_content.end()) {
                  vector<vmesh::GlobalID> siblings;
                  populations[popID].vmesh.getSiblings(populations[popID].vmesh.getParent(sorted[r][b]),siblings);
                  ghostBlockList.insert(siblings.begin(),siblings.end());
                  continue;
               }
            }

            // Parent or grandparent does not exist, add this block
            ghostBlockList.insert(sorted[r][b]);
         }         
      }

      // Add missing no-content blocks
      for (unordered_set<vmesh::GlobalID>::const_iterator it=ghostBlockList.begin();
           it != ghostBlockList.end(); ++it) {
         // Parent already exists
         if (populations[popID].vmesh.getLocalID(populations[popID].vmesh.getParent(*it)) != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) continue;

         // If any children exist, make sure they all exist
         std::vector<vmesh::GlobalID> children;
         populations[popID].vmesh.getChildren(*it,children);
         bool childrensExist = false;
         for (size_t c=0; c<children.size(); ++c) {
            if (populations[popID].vmesh.getLocalID(children[c]) != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
               childrensExist = true;
               break;
            }  
         }
         if (childrensExist == true) {
            // Attempt to add all children, only succeeds if the 
            // children does not exist
            for (size_t c=0; c<children.size(); ++c) add_velocity_block(children[c],popID);
            continue;
         }

         add_velocity_block(*it,popID);
      }
   }

   #endif

   void SpatialCell::adjustSingleCellVelocityBlocks(const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      //neighbor_ptrs is empty as we do not have any consistent
      //data in neighbours yet, adjustments done only based on velocity
      //space. TODO: should this delete blocks or not? Now not
      std::vector<SpatialCell*> neighbor_ptrs;
      update_velocity_block_content_lists(popID);
      adjust_velocity_blocks(neighbor_ptrs,popID,false);
   }

   void SpatialCell::coarsen_block(const vmesh::GlobalID& parent,const std::vector<vmesh::GlobalID>& children,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif

      // First create the parent (coarse) block and grab pointer to its data.
      // add_velocity_block initializes data to zero values.
      if (add_velocity_block(parent,popID) == false) return;
      vmesh::LocalID parentLID = populations[popID].vmesh.getLocalID(parent);
      Realf* parent_data = get_data(popID)+parentLID*SIZE_VELBLOCK;

      // Calculate children (fine) block local IDs, some of the children may not exist
      for (size_t c=0; c<children.size(); ++c) {
         vmesh::LocalID childrenLID = populations[popID].vmesh.getLocalID(children[c]);
         if (childrenLID == SpatialCell::invalid_local_id()) continue;
         Realf* data = get_data(popID)+childrenLID*SIZE_VELBLOCK;

         const int i_oct = c % 2;
         const int j_oct = (c/2) % 2;
         const int k_oct = c / 4;

         /*for (int k=0; k<WID; k+=2) for (int j=0; j<WID; j+=2) for (int i=0; i<WID; i+=2) {
            cerr << "\t" << i_oct*2+i/2 << ' ' << j_oct*2+j/2 << ' ' << k_oct*2+k/2 << " gets values from" << endl;
            
            // Sum the values in 8 cells that correspond to the same call in parent block
            Realf sum = 0;
            for (int kk=0; kk<2; ++kk) for (int jj=0; jj<2; ++jj) for (int ii=0; ii<2; ++ii) {
               cerr << "\t\t" << i+ii << ' ' << j+jj << ' ' << k+kk << endl;
               sum += data[vblock::index(i+ii,j+jj,k+kk)];
            }

            parent_data[vblock::index(i_oct*2+i/2,j_oct*2+j/2,k_oct*2+k/2)] = sum/8;
         }*/

         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
            parent_data[vblock::index(i_oct*2+i/2,j_oct*2+j/2,k_oct*2+k/2)] += data[vblock::index(i,j,k)]/8.0;
         }
      }

      // Remove the children
      for (size_t c=0; c<children.size(); ++c) {
         remove_velocity_block(children[c],popID);
      }
   }

   void SpatialCell::coarsen_blocks(amr_ref_criteria::Base* refCriterion,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      // Sort blocks according to their refinement levels
      vector<vector<vmesh::GlobalID> > blocks(populations[popID].vmesh.getMaxAllowedRefinementLevel()+1);

      for (vmesh::LocalID blockLID=0; blockLID<get_number_of_velocity_blocks(popID); ++blockLID) {
         vmesh::GlobalID blockGID = populations[popID].vmesh.getGlobalID(blockLID);
         uint8_t r = populations[popID].vmesh.getRefinementLevel(blockGID);
         blocks[r].push_back(blockGID);
      }

      // This is how much neighbor data we use when evaluating refinement criteria
      const int PAD=1;
      Realf array[(WID+2*PAD)*(WID+2*PAD)*(WID+2*PAD)];

      // Evaluate refinement criterion for velocity blocks, starting from 
      // the highest refinement level blocks
      for (size_t r=blocks.size()-1; r>=1; --r) {
         // List of blocks that can be coarsened
         //vector<vmesh::GlobalID> coarsenList;
         unordered_set<vmesh::GlobalID> coarsenList;

         // Evaluate refinement criterion for all blocks
         for (size_t b=0; b<blocks[r].size(); ++b) {
            const vmesh::GlobalID blockGID = blocks[r][b];
            fetch_data<1>(blockGID,populations[popID].vmesh,get_data(popID),array);
            if (refCriterion->evaluate(array,popID) < Parameters::amrCoarsenLimit) coarsenList.insert(blockGID);
         }

         // List of blocks created and removed during the coarsening. The first element (=key) 
         // is the global ID of the parent block that is created, the vector (=value) contains 
         // the global IDs of the children that will be removed
         unordered_map<vmesh::GlobalID,vector<vmesh::GlobalID> > allowCoarsen;

         // A block can refined only if all blocks in the octant allow it, 
         // and only if the octant (in which the block belongs to) neighbors 
         // do not have children
         for (unordered_set<vmesh::GlobalID>::const_iterator it=coarsenList.begin(); it!=coarsenList.end(); ++it) {
            vector<vmesh::GlobalID> siblings;
            populations[popID].vmesh.getSiblings(*it,siblings);
            bool allows=true;
            for (size_t s=0; s<siblings.size(); ++s) {
               // Skip non-existing blocks
               if (populations[popID].vmesh.getLocalID(siblings[s]) == invalid_local_id()) {
                  continue;
               }

               // Check that the sibling allows coarsening
               if (coarsenList.find(siblings[s]) == coarsenList.end()) {
                  allows = false;
                  break;
               }
            }

            // Check that the mesh structure allows coarsening
            if (populations[popID].vmesh.coarsenAllowed(*it) == false) continue;

            // If all siblings allow coarsening, add block to coarsen list
            if (allows == true) {
               allowCoarsen.insert(make_pair(populations[popID].vmesh.getParent(*it),siblings));
            }
         }

         cerr << "ref level " << r << " has " << allowCoarsen.size() << " blocks for coarsening" << endl;
         for (unordered_map<vmesh::GlobalID,vector<vmesh::GlobalID> >::const_iterator it=allowCoarsen.begin(); it!=allowCoarsen.end(); ++it) {
            coarsen_block(it->first,it->second,popID);
         }
      }
   }

   /*!
    Returns true if given velocity block has enough of a distribution function.
    Returns false if the value of the distribution function is too low in every
    sense in given block.
    Also returns false if given block doesn't exist or is an error block.
    */
   bool SpatialCell::compute_block_has_content(const vmesh::GlobalID& blockGID,const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
                                              
      if (blockGID == invalid_global_id()) return false;
      const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID,popID);
      if (blockLID == invalid_local_id()) return false;
            
      bool has_content = false;
      const Real velocity_block_min_value = getVelocityBlockMinValue(popID);
      const Realf* block_data = populations[popID].blockContainer.getData(blockLID);
      for (unsigned int i=0; i<VELOCITY_BLOCK_LENGTH; ++i) {
         if (block_data[i] >= velocity_block_min_value) {
            has_content = true;
            break;
         }
      }
      
      return has_content;
   }
   
   /** Get maximum translation timestep for the given species.
    * @param popID ID of the particle species.
    * @return Maximum timestep calculated by the Vlasov translation.*/
   const Real& SpatialCell::get_max_r_dt(const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      return populations[popID].max_dt[species::MAXRDT];
   }
   
   /** Get maximum acceleration timestep for the given species.
    * @param popID ID of the particle species.
    * @return Maximum timestep calculated by Vlasov acceleration.*/
   const Real& SpatialCell::get_max_v_dt(const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      return populations[popID].max_dt[species::MAXVDT];
   }

   /** Get MPI datatype for sending the cell data.
    * @param cellID Spatial cell (dccrg) ID.
    * @param sender_rank Rank of the MPI process sending data from this cell.
    * @param receiver_rank Rank of the MPI process receiving data to this cell.
    * @param receiving If true, this process is receiving data.
    * @param neighborhood Neighborhood ID.
    * @return MPI datatype that transfers the requested data.*/
   std::tuple<void*, int, MPI_Datatype> SpatialCell::get_mpi_datatype(
                                                                      const CellID cellID,
                                                                      const int sender_rank,
                                                                      const int receiver_rank,
                                                                      const bool receiving,
                                                                      const int neighborhood
      ) {

      std::vector<MPI_Aint> displacements;
      std::vector<int> block_lengths;
      vmesh::LocalID block_index = 0;

      // create datatype for actual data if we are in the first two 
      // layers around a boundary, or if we send for the whole system
      if (this->mpiTransferEnabled && (SpatialCell::mpiTransferAtSysBoundaries==false || this->sysBoundaryLayer ==1 || this->sysBoundaryLayer ==2 )) {
         //add data to send/recv to displacement and block length lists
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE1) != 0) {
            //first copy values in case this is the send operation
            populations[activePopID].N_blocks = populations[activePopID].blockContainer.size();

            // send velocity block list size
            displacements.push_back((uint8_t*) &(populations[activePopID].N_blocks) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE2) != 0) {
            // STAGE1 should have been done, otherwise we have problems...
            if (receiving) {
               //mpi_number_of_blocks transferred earlier
               populations[activePopID].vmesh.setNewSize(populations[activePopID].N_blocks);
            } else {
                //resize to correct size (it will avoid reallocation if it is big enough, I assume)
                populations[activePopID].N_blocks = populations[activePopID].blockContainer.size();
            }

            // send velocity block list
            displacements.push_back((uint8_t*) &(populations[activePopID].vmesh.getGrid()[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::GlobalID) * populations[activePopID].vmesh.size());
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE1) !=0) {
            //Communicate size of list so that buffers can be allocated on receiving side
            if (!receiving) this->velocity_block_with_content_list_size = this->velocity_block_with_content_list.size();
            displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list_size) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE2) !=0) {
            if (receiving) {
               this->velocity_block_with_content_list.resize(this->velocity_block_with_content_list_size);
            }

            //velocity_block_with_content_list_size should first be updated, before this can be done (STAGE1)
            displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::GlobalID)*this->velocity_block_with_content_list_size);
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_DATA) !=0) {
            displacements.push_back((uint8_t*) get_data(activePopID) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Realf) * VELOCITY_BLOCK_LENGTH * populations[activePopID].blockContainer.size());
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::NEIGHBOR_VEL_BLOCK_DATA) != 0) {
            /*We are actually transferring the data of a
            * neighbor. The values of neighbor_block_data
            * and neighbor_number_of_blocks should be set in
            * solver.*/               
            displacements.push_back((uint8_t*) this->neighbor_block_data - (uint8_t*) this);               
            block_lengths.push_back(sizeof(Realf) * VELOCITY_BLOCK_LENGTH* this->neighbor_number_of_blocks);
         }

         // send  spatial cell parameters
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PARAMETERS)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * CellParams::N_SPATIAL_CELL_PARAMS);
         }
         
         // send  spatial cell dimensions
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_DIMENSIONS)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::DX]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }
         
         // send  BGBX BGBY BGBZ and all edge-averaged BGBs
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BGB)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBX]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBX_000_010]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 24);
         }
         
         // send  BGBXVOL BGBYVOL BGBZVOL PERBXVOL PERBYVOL PERBZVOL
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBXVOL]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 6);
         }
         
         // send  EX, EY EZ
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_E)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::EX]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }
         
         // send  EX_DT2, EY_DT2, EZ_DT2
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_EDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::EX_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }
         
         // send  PERBX, PERBY, PERBZ
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PERB)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::PERBX]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }
         
         // send  PERBX_DT2, PERBY_DT2, PERBZ_DT2
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PERBDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::PERBX_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }
         
         // send RHO, RHOVX, RHOVY, RHOVZ
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHO_RHOV)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHO]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 4);
         }
         
         // send RHO_DT2, RHOVX_DT2, RHOVY_DT2, RHOVZ_DT2
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHODT2_RHOVDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHO_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 4);
         }
         
         // send  spatial cell derivatives
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_DERIVATIVES)!=0){
            displacements.push_back((uint8_t*) &(this->derivatives[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * fieldsolver::N_SPATIAL_CELL_DERIVATIVES);
         }
         
         // send  spatial cell BVOL derivatives
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL_DERIVATIVES)!=0){
            displacements.push_back((uint8_t*) &(this->derivativesBVOL[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * bvolderivatives::N_BVOL_DERIVATIVES);
         }
         
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_IOLOCALCELLID)!=0){
            displacements.push_back((uint8_t*) &(this->ioLocalCellId) - (uint8_t*) this);
            block_lengths.push_back(sizeof(uint64_t));
         }
         
         // send Hall term components
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_HALL_TERM)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::EXHALL_000_100]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 12);
         }
         // send electron pressure gradient term components
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_GRADPE_TERM)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::EXGRADPE]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }

         
         // send P tensor diagonal components
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_P)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::P_11]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }
         
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::P_11_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }
         
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHOQ_TOT) != 0) {
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHOQ_TOT]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real));
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PHI) != 0) {
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::PHI]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real)*2);
         }
         
         // send  sysBoundaryFlag
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_SYSBOUNDARYFLAG)!=0){
            displacements.push_back((uint8_t*) &(this->sysBoundaryFlag) - (uint8_t*) this);
            block_lengths.push_back(sizeof(uint));
            displacements.push_back((uint8_t*) &(this->sysBoundaryLayer) - (uint8_t*) this);
            block_lengths.push_back(sizeof(uint));
            displacements.push_back((uint8_t*) &(this->sysBoundaryLayerNew) - (uint8_t*) this);
            block_lengths.push_back(sizeof(int));
         }
         
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_PARAMETERS) !=0) {
            displacements.push_back((uint8_t*) get_block_parameters(activePopID) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * size(activePopID) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
         }
         // Copy particle species metadata
         if ((SpatialCell::mpi_transfer_type & Transfer::POP_METADATA) != 0) {
            for (unsigned int popID=0; popID<populations.size(); ++popID) {
               displacements.push_back((uint8_t*) &(populations[popID].max_dt) - (uint8_t*)this);
               block_lengths.push_back(species::SIZE_DT_ELEMENTS*sizeof(Real));
            }
         }
         // Copy random number generator state variables
         //if ((SpatialCell::mpi_transfer_type & Transfer::RANDOMGEN) != 0) {
         //   displacements.push_back((uint8_t*)get_rng_state_buffer() - (uint8_t*)this);
         //   block_lengths.push_back(256/8);
         //   displacements.push_back((uint8_t*)get_rng_data_buffer() - (uint8_t*)this);
         //   block_lengths.push_back(sizeof(random_data));
         //}
      }

      void* address = this;
      int count;
      MPI_Datatype datatype;
      
      if (displacements.size() > 0) {
         count = 1;
         MPI_Type_create_hindexed(
            displacements.size(),
            &block_lengths[0],
            &displacements[0],
            MPI_BYTE,
            &datatype
         );
      } else {
         count = 0;
         datatype = MPI_BYTE;
      }

      return std::make_tuple(address,count,datatype);
   }
   
   /** Get random number generator data buffer.
    * @return Random number generator data buffer.*/
   //random_data* SpatialCell::get_rng_data_buffer() {
   //   return &rngDataBuffer;
   //}

   /** Get random number generator state buffer.
    * @return Random number generator state buffer.*/
   //char* SpatialCell::get_rng_state_buffer() {
   //   return rngStateBuffer;
   //}

   /**< Minimum value of distribution function in any phase space cell 
    * of a velocity block for the block to be considered to have content.
    * @param popID ID of the particle species.
    * @return Sparse min value for this species.*/
   Real SpatialCell::getVelocityBlockMinValue(const int& popID) const {
      return populations[popID].velocityBlockMinValue;
   }

   void SpatialCell::merge_values_recursive(const int& popID,vmesh::GlobalID parentGID,vmesh::GlobalID blockGID,
                                            uint8_t refLevel,bool recursive,const Realf* data,
					    std::set<vmesh::GlobalID>& blockRemovalList) {
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == invalid_global_id()) {
         cerr << "merge_values_recursive called for GID=" << blockGID << " at r=" << (int)refLevel << " parent=" << parentGID << endl;
      }
      #endif

      // Get all possible children:
      vector<vmesh::GlobalID> childrenGIDs;
      populations[popID].vmesh.getChildren(blockGID,childrenGIDs);

      #warning FIXME activePopID is not correct here
      
      // Check if any of block's children exist:
      bool hasChildren = false;
      for (size_t c=0; c<childrenGIDs.size(); ++c) {
         if (populations[activePopID].vmesh.getLocalID(childrenGIDs[c]) != invalid_local_id()) {
            hasChildren = true;
            break;
         }
      }
      
      /*      // TEST
       for (size_t c=0; c<childrenGIDs.size(); ++c) {
       bool hasGrandChildren=false;
       vector<vmesh::GlobalID> grandChildren;
       if (vmesh.getLocalID(childrenGIDs[c]) != vmesh.invalidLocalID()) continue;
       vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getChildren(childrenGIDs[c],grandChildren);
       for (size_t cc=0; cc<grandChildren.size(); ++cc) {
	    if (vmesh.getLocalID(grandChildren[cc]) != vmesh.invalidLocalID()) {
       hasGrandChildren = true;
	       break;
	    }
       }
       if (hasGrandChildren==true) {
	    std::cerr << "block at r=" << (int)refLevel << " has lost grand children" << std::endl;
       }
       }
       // END TEST */

      // No children, try to merge values to this block:
      if (hasChildren == false) {
         vmesh::LocalID blockLID = populations[activePopID].vmesh.getLocalID(blockGID);
	 
         #ifdef DEBUG_SPATIAL_CELL
         if (blockLID == invalid_local_id()) {
            cerr << "ERROR: Failed to merge values, block does not exist!" << endl;
            cerr << "\t exiting from " << __FILE__ << ':' << __LINE__ << endl;
            exit(1);
         }
         #endif

         Realf* myData = populations[activePopID].blockContainer.getData(blockLID);
         if (parentGID == blockGID) {
            // If we enter here, the block is at the lowest refinement level.
            // If the block does not have enough content, flag it for removal
            //#warning REMOVED for debugging
            /*for (unsigned int i=0; i<WID3; ++i) {
             if (myData[i] >= SpatialCell::velocity_block_min_value) return;
             }
             blockRemovalList.insert(blockGID);*/
         } else {
            // Merge values to this block
            for (int i=0; i<WID3; ++i) myData[i] += data[i];
         }
         return;
      }
      
      // Iterate over all octants, each octant corresponds to a different child:
      bool removeBlock = false;
      for (int k_oct=0; k_oct<2; ++k_oct) for (int j_oct=0; j_oct<2; ++j_oct) for (int i_oct=0; i_oct<2; ++i_oct) {
         // Copy data belonging to the octant to a temporary array:
         Realf array[WID3];
         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
            array[vblock::index(i,j,k)] =   data[vblock::index(i_oct*2+i/2,j_oct*2+j/2,k_oct*2+k/2)];
         }
         
         // Send the data to the child:
         const int octant = k_oct*4 + j_oct*2 + i_oct;
         merge_values_recursive(popID,blockGID,childrenGIDs[octant],refLevel+1,true,array,blockRemovalList);
      }

      // Data merged to children, block can be removed
      blockRemovalList.insert(blockGID);
   }

   void SpatialCell::merge_values(const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      const uint8_t maxRefLevel = populations[popID].vmesh.getMaxAllowedRefinementLevel();

      for (int i=0; i<WID3; ++i) null_block_data[i] = 0;
      
      // Sort blocks according to their refinement levels:
      vector<vector<vmesh::GlobalID> > blocks(maxRefLevel+1);
      for (vmesh::LocalID blockLID=0; blockLID<get_number_of_velocity_blocks(popID); ++blockLID) {
         const vmesh::GlobalID blockGID = populations[popID].vmesh.getGlobalID(blockLID);
         
         if (blockGID == SpatialCell::invalid_global_id()) {
            cerr << "got invalid global id from mesh!" << endl;
            continue;
         }

         uint8_t refLevel = populations[popID].vmesh.getRefinementLevel(blockGID);
         blocks[refLevel].push_back(blockGID);
      }
      
      set<vmesh::GlobalID> blockRemovalList;
      for (uint8_t refLevel=0; refLevel<blocks.size()-1; ++refLevel) {
         for (size_t b=0; b<blocks[refLevel].size(); ++b) {
            const vmesh::GlobalID blockGID = blocks[refLevel][b];
            const vmesh::LocalID  blockLID = populations[popID].vmesh.getLocalID(blockGID);
            if (blockLID == invalid_local_id()) continue;
            
            const Realf* data = populations[popID].blockContainer.getData(blockLID);
            merge_values_recursive(popID,blockGID,blockGID,refLevel,true,data,blockRemovalList);
         }
      }
      
      cerr << "should remove " << blockRemovalList.size() << " blocks" << endl;
      for (set<vmesh::GlobalID>::const_iterator it=blockRemovalList.begin(); it!=blockRemovalList.end(); ++it) {
         //remove_velocity_block(*it);
      }
   }
   
   void SpatialCell::add_values(const vmesh::GlobalID& targetGID,
                                std::unordered_map<vmesh::GlobalID,Realf[(WID+2)*(WID+2)*(WID+2)]>& sourceData,
                                const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      vmesh::LocalID targetLID = get_velocity_block_local_id(targetGID,popID);
      if (targetLID == SpatialCell::invalid_local_id()) {
         std::cerr << "error has occurred" << std::endl;
         return;
      }

      Realf* targetData = get_data(popID)+targetLID*SIZE_VELBLOCK;

      // Add data from all same level blocks
      vector<vmesh::GlobalID> neighborIDs;
      populations[popID].vmesh.getNeighborsAtSameLevel(targetGID,neighborIDs);

      std::unordered_map<vmesh::GlobalID,Realf[(WID+2)*(WID+2)*(WID+2)]>::iterator it;

      it = sourceData.find(neighborIDs[vblock::nbrIndex( 0, 0, 0)]); // This block
      if (it != sourceData.end()) {
         Realf* source = it->second;
         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
            targetData[vblock::index(i,j,k)] += source[vblock::padIndex<1>(i+1,j+1,k+1)];
         }
      }
      
      // ***** face neighbors ***** //
      for (int i_off=-1; i_off<2; i_off+=2) {
         it = sourceData.find(neighborIDs[vblock::nbrIndex(i_off,0,0)]);
         if (it == sourceData.end()) continue;
         
         int i_trgt = 0;
         int i_src  = WID+1;
         if (i_off > 0) {i_trgt = WID-1; i_src=0;}
         Realf* source = it->second;
         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) {
            targetData[vblock::index(i_trgt,j,k)] += source[vblock::padIndex<1>(i_src,j+1,k+1)];
         }
      }
      
      for (int j_off=-1; j_off<2; j_off+=2) {
         it = sourceData.find(neighborIDs[vblock::nbrIndex(0,j_off,0)]);
         if (it == sourceData.end()) continue;

         int j_trgt = 0;
         int j_src  = WID+1;
         if (j_off > 0) {j_trgt = WID-1; j_src=0;}
         Realf* source = it->second;
         for (int k=0; k<WID; ++k) for (int i=0; i<WID; ++i) {
            targetData[vblock::index(i,j_trgt,k)] += source[vblock::padIndex<1>(i+1,j_src,k+1)];
         }
      }
      
      for (int k_off=-1; k_off<2; k_off+=2) {
         it = sourceData.find(neighborIDs[vblock::nbrIndex(0,0,k_off)]);
         if (it == sourceData.end()) continue;
         
         int k_trgt = 0;
         int k_src  = WID+1;
         if (k_off > 0) {k_trgt = WID-1; k_src=0;}
         Realf* source = it->second;
         for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
            targetData[vblock::index(i,j,k_trgt)] += source[vblock::padIndex<1>(i+1,j+1,k_src)];
         }
      }

      // ***** edge neighbors ***** //
      for (int j_off=-1; j_off<2; j_off+=2) for (int i_off=-1; i_off<2; i_off+=2) {
         it = sourceData.find(neighborIDs[vblock::nbrIndex(i_off,j_off,0)]);
         if (it == sourceData.end()) continue;
         Realf* source = it->second;
         
         const int i_trgt = max(i_off,0) * (WID-1);
         const int j_trgt = max(j_off,0) * (WID-1);
         const int i_src  = max(-i_off,0) * (WID+1);
         const int j_src  = max(-j_off,0) * (WID+1);
         for (int k=0; k<WID; ++k) {
            targetData[vblock::index(i_trgt,j_trgt,k)] += source[vblock::padIndex<1>(i_src,j_src,k+1)];
         }
      }
      
      for (int k_off=-1; k_off<2; k_off+=2) for (int i_off=-1; i_off<2; i_off+=2) {
         it = sourceData.find(neighborIDs[vblock::nbrIndex(i_off,0,k_off)]);
         if (it == sourceData.end()) continue;
         Realf* source = it->second;
         
         const int i_trgt = max(i_off,0) * (WID-1);
         const int k_trgt = max(k_off,0) * (WID-1);
         const int i_src  = max(-i_off,0) * (WID+1);
         const int k_src  = max(-k_off,0) * (WID+1);
         for (int j=0; j<WID; ++j) {
            targetData[vblock::index(i_trgt,j,k_trgt)] += source[vblock::padIndex<1>(i_src,j+1,k_src)];
         }
      }
      
      for (int k_off=-1; k_off<2; k_off+=2) for (int j_off=-1; j_off<2; j_off+=2) {
         it = sourceData.find(neighborIDs[vblock::nbrIndex(0,j_off,k_off)]);
         if (it == sourceData.end()) continue;
         Realf* source = it->second;
         
         const int j_trgt = max(j_off,0) * (WID-1);
         const int k_trgt = max(k_off,0) * (WID-1);
         const int j_src  = max(-j_off,0) * (WID+1);
         const int k_src  = max(-k_off,0) * (WID+1);
         for (int i=0; i<WID; ++i) {
            targetData[vblock::index(i,j_trgt,k_trgt)] += source[vblock::padIndex<1>(i+1,j_src,k_src)];
         }
      }

      // ***** corner neighbors ***** //
      for (int k_off=-1; k_off<2; k_off+=2) for (int j_off=-1; j_off<2; j_off+=2) for (int i_off=-1; i_off<2; i_off+=2) {
         it = sourceData.find(neighborIDs[vblock::nbrIndex(i_off,j_off,k_off)]);
         if (it == sourceData.end()) continue;

         const int i_trgt = max(i_off,0) * (WID-1);
         const int j_trgt = max(j_off,0) * (WID-1);
         const int k_trgt = max(k_off,0) * (WID-1);
         const int i_src  = max(-i_off,0) * (WID+1);
         const int j_src  = max(-j_off,0) * (WID+1);
         const int k_src  = max(-k_off,0) * (WID+1);
         Realf* source = it->second;
         targetData[vblock::index(i_trgt,j_trgt,k_trgt)] += source[vblock::padIndex<1>(i_src,j_src,k_src)];
         
         source[vblock::padIndex<1>(i_src,j_src,k_src)]=0;
      }
      
      // Exit if the block is at base grid level
      if (populations[popID].vmesh.getRefinementLevel(targetGID) == 0) {
         return;
      }
      const int octant = populations[popID].vmesh.getOctant(targetGID);
      
      int parentIndex[3];
      parentIndex[2] = 1 + 2*(octant / 4);
      parentIndex[1] = 1 + 2*((octant - 4*(octant/4))/2);
      parentIndex[0] = 1 + 2*(octant % 2);
      
      // Add data from all coarser blocks
      const Realf face_mul = 2.0;
      const Realf edge_mul = 4.0;
      const Realf corn_mul = 8.0;

      const vmesh::GlobalID targetParentGID = populations[popID].vmesh.getParent( neighborIDs[vblock::nbrIndex( 0, 0, 0)] );
      vmesh::GlobalID parentGID = targetParentGID;
      
      it = sourceData.find(parentGID);
      if (it != sourceData.end()) {
         Realf* source = it->second;
         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
            targetData[vblock::index(i,j,k)] += source[vblock::padIndex<1>(parentIndex[0]+i/2,parentIndex[1]+j/2,parentIndex[2]+k/2)];
         }
         
         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) source[vblock::padIndex<1>(parentIndex[0]+i/2,parentIndex[1]+j/2,parentIndex[2]+k/2)]=0;
      }
      
      // ***** face neighbors ***** //
      for (int i_off=-1; i_off<2; i_off+=2) {
         parentGID = populations[popID].vmesh.getParent( neighborIDs[vblock::nbrIndex(i_off,0,0)] );
         if (parentGID == targetParentGID) continue;
         it = sourceData.find(parentGID);
         if (it == sourceData.end()) continue;

         int i_trgt = max( i_off,0) * (WID-1);
         int i_src  = max(-i_off,0) * (WID+1);
         
         Realf* source = it->second;
         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) {
            targetData[vblock::index(i_trgt,j,k)] += face_mul*source[vblock::padIndex<1>(i_src,parentIndex[1]+j/2,parentIndex[2]+k/2)];
         }
      }

      for (int j_off=-1; j_off<2; j_off+=2) {
         parentGID = populations[popID].vmesh.getParent( neighborIDs[vblock::nbrIndex(0,j_off,0)] );
         if (parentGID == targetParentGID) continue;
         it = sourceData.find(parentGID);
         if (it == sourceData.end()) continue;
         
         int j_trgt = max( j_off,0) * (WID-1);
         int j_src  = max(-j_off,0) * (WID+1);
         
         Realf* source = it->second;
         for (int k=0; k<WID; ++k) for (int i=0; i<WID; ++i) {
            targetData[vblock::index(i,j_trgt,k)] += face_mul*source[vblock::padIndex<1>(parentIndex[0]+i/2,j_src,parentIndex[2]+k/2)];
         }
      }
      
      for (int k_off=-1; k_off<2; k_off+=2) {
         parentGID = populations[popID].vmesh.getParent( neighborIDs[vblock::nbrIndex(0,0,k_off)] );
         if (parentGID == targetParentGID) continue;
         it = sourceData.find(parentGID);
         if (it == sourceData.end()) continue;
         
         int k_trgt = max( k_off,0) * (WID-1);
         int k_src  = max(-k_off,0) * (WID+1);
         
         Realf* source = it->second;
         for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
            targetData[vblock::index(i,j,k_trgt)] += face_mul*source[vblock::padIndex<1>(parentIndex[0]+i/2,parentIndex[1]+j/2,k_src)];
         }
      }

      // ***** edge neighbors ***** //
      // A refined block in an octant is only allowed copy 
      // values from one coarser edge neighbor
      int i_off,j_off,k_off;
      if (octant == 0 || octant == 4) {i_off=-1; j_off=-1;}
      if (octant == 1 || octant == 5) {i_off=+1; j_off=-1;}
      if (octant == 2 || octant == 6) {i_off=-1; j_off=+1;}
      if (octant == 3 || octant == 7) {i_off=+1; j_off=+1;}
      Realf* source = NULL;
      parentGID = populations[popID].vmesh.getParent( neighborIDs[vblock::nbrIndex(i_off,j_off,0)] );
      if (parentGID != targetParentGID) {
         it = sourceData.find(parentGID);
         if (it != sourceData.end()) {
            source = it->second;
            const int i_trgt = max(i_off,0) * (WID-1);
            const int j_trgt = max(j_off,0) * (WID-1);
            const int i_src  = max(-i_off,0) * (WID+1);
            const int j_src  = max(-j_off,0) * (WID+1);
            for (int k=0; k<WID; ++k) {
               targetData[vblock::index(i_trgt,j_trgt,k)] += edge_mul*source[vblock::padIndex<1>(i_src,j_src,parentIndex[2]+k/2)];
            }
         }
      }
      
      if (octant == 0 || octant == 2) {i_off=-1; k_off=-1;}
      if (octant == 1 || octant == 3) {i_off=+1; k_off=-1;}
      if (octant == 4 || octant == 6) {i_off=-1; k_off=+1;}
      if (octant == 5 || octant == 7) {i_off=+1; k_off=+1;}
      source = NULL;
      parentGID = populations[popID].vmesh.getParent( neighborIDs[vblock::nbrIndex(i_off,0,k_off)] );
      if (parentGID != targetParentGID) {
         it = sourceData.find(parentGID);
         if (it != sourceData.end()) {
            source = it->second;
            const int i_trgt = max(i_off,0) * (WID-1);
            const int k_trgt = max(k_off,0) * (WID-1);
            const int i_src  = max(-i_off,0) * (WID+1);
            const int k_src  = max(-k_off,0) * (WID+1);
            for (int j=0; j<WID; ++j) {
               targetData[vblock::index(i_trgt,j,k_trgt)] += edge_mul*source[vblock::padIndex<1>(i_src,parentIndex[1]+j/2,k_src)];
            }
         }
      }

      if (octant == 0 || octant == 1) {j_off=-1; k_off=-1;}
      if (octant == 2 || octant == 3) {j_off=+1; k_off=-1;}
      if (octant == 4 || octant == 5) {j_off=-1; k_off=+1;}
      if (octant == 6 || octant == 7) {j_off=+1; k_off=+1;}
      source = NULL;
      parentGID = populations[popID].vmesh.getParent( neighborIDs[vblock::nbrIndex(0,j_off,k_off)] );
      if (parentGID != targetParentGID) {
         it = sourceData.find(parentGID);
         if (it != sourceData.end()) {
            source = it->second;
            const int j_trgt = max(j_off,0) * (WID-1);
            const int k_trgt = max(k_off,0) * (WID-1);
            const int j_src  = max(-j_off,0) * (WID+1);
            const int k_src  = max(-k_off,0) * (WID+1);
            for (int i=0; i<WID; ++i) {
               targetData[vblock::index(i,j_trgt,k_trgt)] += edge_mul*source[vblock::padIndex<1>(parentIndex[0]+i/2,j_src,k_src)];
            }
         }
      }
      
      // ***** corner neighbors ***** //
      // A refined block in an octant is allowed to copy
      // values from a single coarser corner neighbor
      switch (octant) {
       case 0: i_off=-1; j_off=-1; k_off=-1; break;
       case 1: i_off=+1; j_off=-1; k_off=-1; break;
       case 2: i_off=-1; j_off=+1; k_off=-1; break;
       case 3: i_off=+1; j_off=+1; k_off=-1; break;
       case 4: i_off=-1; j_off=-1; k_off=+1; break;
       case 5: i_off=+1; j_off=-1; k_off=+1; break;
       case 6: i_off=-1; j_off=+1; k_off=+1; break;
       case 7: i_off=+1; j_off=+1; k_off=+1; break;
      }
      source = NULL;
      parentGID = populations[popID].vmesh.getParent( neighborIDs[vblock::nbrIndex(i_off,j_off,k_off)] );
      if (parentGID != targetParentGID) {
         it = sourceData.find(parentGID);
         if (it != sourceData.end()) {
            source = it->second;
            const int i_trgt = max(i_off,0) * (WID-1);
            const int j_trgt = max(j_off,0) * (WID-1);
            const int k_trgt = max(k_off,0) * (WID-1);
            const int i_src  = max(-i_off,0) * (WID+1);
            const int j_src  = max(-j_off,0) * (WID+1);
            const int k_src  = max(-k_off,0) * (WID+1);
            targetData[vblock::index(i_trgt,j_trgt,k_trgt)] += corn_mul*source[vblock::padIndex<1>(i_src,j_src,k_src)];
         }
      }
   }
   
   /** Prepares this spatial cell to receive the velocity grid over MPI.
    * At this stage we have received a new block list over MPI into
    * mpi_velocity_block_list, but the rest of the cell structures
    * have not been adapted to this new list. Here we re-initialize
    * the cell with empty blocks based on the new list.*/
   void SpatialCell::prepare_to_receive_blocks(const int& popID) {
      populations[popID].vmesh.setGrid();
      populations[popID].blockContainer.setSize(populations[popID].vmesh.size());

      Real* parameters = get_block_parameters(popID);
      
      // Set velocity block parameters:
      for (vmesh::LocalID blockLID=0; blockLID<size(popID); ++blockLID) {
         const vmesh::GlobalID blockGID = get_velocity_block_global_id(blockLID,popID);
         parameters[BlockParams::VXCRD] = get_velocity_block_vx_min(popID,blockGID);
         parameters[BlockParams::VYCRD] = get_velocity_block_vy_min(popID,blockGID);
         parameters[BlockParams::VZCRD] = get_velocity_block_vz_min(popID,blockGID);
         populations[popID].vmesh.getCellSize(blockGID,&(parameters[BlockParams::DVX]));
         parameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      }
   }

   void SpatialCell::refine_block(const vmesh::GlobalID& blockGID,std::map<vmesh::GlobalID,vmesh::LocalID>& insertedBlocks,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == invalid_global_id()) {
         std::cerr << "invalid global ID, skip refinement" << std::endl;
         return;
      }
      #endif

      // Tell mesh to refine the given block. In return we get the erased 
      // and inserted blocks. Note that multiple blocks can be removed 
      // (and inserted) due to induced refinement. There are eight entries 
      // in newInserted (the children) for each entry in erasedBlocks.
      std::set<vmesh::GlobalID> erasedBlocks;
      std::map<vmesh::GlobalID,vmesh::LocalID> newInserted;
      if (populations[popID].vmesh.refine(blockGID,erasedBlocks,newInserted) == false) {
         return;
      }

      // Resize the block container, this preserves old data.
      const size_t newBlocks = newInserted.size()-erasedBlocks.size();
      populations[popID].blockContainer.setSize(populations[popID].blockContainer.size() + newBlocks);

      std::map<vmesh::GlobalID,vmesh::LocalID>::const_iterator ins=newInserted.begin();
      for (std::set<vmesh::GlobalID>::const_iterator er=erasedBlocks.begin(); er!=erasedBlocks.end(); ++er) {
         for (int child=0; child<8; ++child) {
            // Copy / interpolate data from old (coarse) block to new refined blocks.
            
            
            // Set refined block parameters
            Real* blockParams = populations[popID].blockContainer.getParameters(ins->second);
            populations[popID].vmesh.getBlockCoordinates(ins->first,blockParams);
            populations[popID].vmesh.getCellSize(ins->first,blockParams+3);
            
            ++ins;
         }
      }
      
      for (std::map<vmesh::GlobalID,vmesh::LocalID>::iterator it=newInserted.begin(); it!=newInserted.end(); ++it) {
         // Set refined block parameters
         Real* blockParams = populations[popID].blockContainer.getParameters(it->second);
         populations[popID].vmesh.getBlockCoordinates(it->first,blockParams);
         populations[popID].vmesh.getCellSize(it->first,blockParams+3);
         
      }

      insertedBlocks.insert(newInserted.begin(),newInserted.end());
   }

   /** Set the particle species SpatialCell should use in functions that 
    * use the velocity mesh.
    * @param popID Population ID.
    * @return If true, the new species is in use.*/
   bool SpatialCell::setCommunicatedSpecies(const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= getObjectWrapper().particleSpecies.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds species.size() " << getObjectWrapper().particleSpecies.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      activePopID = popID;
      return true;
   }
   
   /** Set maximum translation timestep for a particle species.
    * This function is called during Vlasov translation.
    * @param popID ID of the particle species.
    * @param value New maximum timestep.*/
   void SpatialCell::set_max_r_dt(const int& popID,const Real& value) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      populations[popID].max_dt[species::MAXRDT] = value;
   }

   /** Set maximum acceleration timestep for a particle species.
    * This function is called during Vlasov acceleration.
    * @param popID ID of the particle species.
    * @param value New maximum timestep.*/
   void SpatialCell::set_max_v_dt(const int& popID,const Real& value) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      populations[popID].max_dt[species::MAXVDT] = value;
   }

   /**  Purges extra capacity from block vectors. It sets size to
    * num_blocks * block_allocation_factor (if capacity greater than this), 
    * and also forces capacity to this new smaller value.
    * @return True on success.*/
   bool SpatialCell::shrink_to_fit() {
      bool success = true;
      for (size_t p=0; p<populations.size(); ++p) {
         const size_t amount 
            = 2 + populations[p].blockContainer.size() 
            * populations[p].blockContainer.getBlockAllocationFactor();

         // Allow capacity to be a bit large than needed by number of blocks, shrink otherwise
         if (populations[p].blockContainer.capacity() > amount * VELOCITY_BLOCK_LENGTH) 
            if (populations[p].blockContainer.recapacitate(amount) == false) success = false;
      }
      return success;
   }

   /** Update the two lists containing blocks with content, and blocks without content.
    * @see adjustVelocityBlocks */
   void SpatialCell::update_velocity_block_content_lists(const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      velocity_block_with_content_list.clear();
      velocity_block_with_no_content_list.clear();
      
      for (vmesh::LocalID block_index=0; block_index<populations[popID].vmesh.size(); ++block_index) {
         const vmesh::GlobalID globalID = populations[popID].vmesh.getGlobalID(block_index);
         if (compute_block_has_content(globalID,popID)){
            velocity_block_with_content_list.push_back(globalID);
         } else {
            velocity_block_with_no_content_list.push_back(globalID);
         }
      }
   }
   
   void SpatialCell::printMeshSizes() {
      cerr << "SC::printMeshSizes:" << endl;
      for (size_t p=0; p<populations.size(); ++p) {
         cerr << "\t pop " << p << " " << populations[p].vmesh.size() << ' ' << populations[p].blockContainer.size() << endl;  
      }
      cerr << "\t temp sizes are " << vmeshTemp.size() << ' ' << blockContainerTemp.size() << endl;      
   }

   /** Initialize the velocity mesh of the chosen particle population.
    * @param popID Population ID.
    * @param v_limits Velocity mesh min/max limits.
    * @param meshSize Number of blocks in each coordinate direction at base grid level.
    * @param blockSize Number of velocity cells in a block per coordinate direction.
    * @param f_min Sparse mesh threshold value.
    * @param maxRefLevel Maximum allowed mesh refinement level.*/
   bool SpatialCell::initialize_mesh() {
      bool success = true;
      for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         const species::Species& spec = getObjectWrapper().particleSpecies[popID];
         if (populations[popID].vmesh.initialize(spec.velocityMesh,getObjectWrapper().velocityMeshes) == false) {
            success = false;
         }
      }
      
      return success;
   }
   
   /** Updates minValue based on algorithm value from parameters (see parameters.cpp).
    * @param popID ID of the particle species.*/
   void SpatialCell::updateSparseMinValue(const int& popID) {
      if ( P::sparseDynamicAlgorithm == 1 || P::sparseDynamicAlgorithm == 2 ) {
         // Linear algorithm for the minValue: y=kx+b
         const Real k = (P::sparseDynamicMinValue2 - P::sparseDynamicMinValue1) / (P::sparseDynamicBulkValue2 - P::sparseDynamicBulkValue1);
         const Real b = P::sparseDynamicMinValue1 - k * P::sparseDynamicBulkValue1;
         Real x;
         if ( P::sparseDynamicAlgorithm == 1 ) {
            x = this->parameters[CellParams::RHO];
         } else {
            x = this->get_number_of_velocity_blocks(popID);
         }
         const Real newMinValue = k*x+b;
         if( newMinValue < P::sparseDynamicMinValue1 ) { // Compare against the min minValue
            populations[popID].velocityBlockMinValue = P::sparseDynamicMinValue1;
         } else if( newMinValue > P::sparseDynamicMinValue2 ) { // Compare against the max minValue
            populations[popID].velocityBlockMinValue = P::sparseDynamicMinValue2;
         } else {
            populations[popID].velocityBlockMinValue = newMinValue;
         }
         return;
      } else {
         populations[popID].velocityBlockMinValue = getObjectWrapper().particleSpecies[popID].sparseMinValue;
         return;
      }
      return;
   }

} // namespace spatial_cell
