/*!
Space for static variables of spatial cell class for Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute


*/

#include <unordered_set>

#include "spatial_cell.hpp"
#include "velocity_blocks.h"

using namespace std;

namespace spatial_cell {
   Real SpatialCell::velocity_block_min_value = 0;    
   uint64_t SpatialCell::mpi_transfer_type = 0;
   bool SpatialCell::mpiTransferAtSysBoundaries = false;

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
    * modify local cell.*/
   void SpatialCell::adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors, bool doDeleteEmptyBlocks) {
      //  This set contains all those cellids which have neighbors in any
      //  of the 6-dimensions Actually, we would only need to add
      //  local blocks with no content here, as blocks with content
      //  do not need to be created and also will not be removed as
      //  we only check for removal for blocks with no content
      //boost::unordered_set<vmesh::GlobalID> neighbors_have_content;
      std::unordered_set<vmesh::GlobalID> neighbors_have_content;
      
      #ifdef AMR
      for (vmesh::LocalID block_index=0; block_index<velocity_block_with_content_list.size(); ++block_index) {
	 vmesh::GlobalID blockGID = velocity_block_with_content_list[block_index];
	 vector<vmesh::GlobalID> neighborGIDs;
	 vmesh.getNeighborsAtSameLevel(blockGID,neighborGIDs);
	 neighbors_have_content.insert(neighborGIDs.begin(),neighborGIDs.end());
	 neighbors_have_content.insert(blockGID);
      }
      #else
      //add neighbor content info for velocity space neighbors to map. We loop over blocks
      //with content and raise the neighbors_have_content for
      //itself, and for all its neighbors
      for (vmesh::LocalID block_index=0; block_index<velocity_block_with_content_list.size(); ++block_index) {
	 vmesh::GlobalID block = velocity_block_with_content_list[block_index];
	 
	 const velocity_block_indices_t indices = get_velocity_block_indices(block);
	 neighbors_have_content.insert(block); //also add the cell itself
	          
	 for (int offset_vx=-P::sparseBlockAddWidthV;offset_vx<=P::sparseBlockAddWidthV;offset_vx++) {
	    for (int offset_vy=-P::sparseBlockAddWidthV;offset_vy<=P::sparseBlockAddWidthV;offset_vy++) {
	       for (int offset_vz=-P::sparseBlockAddWidthV;offset_vz<=P::sparseBlockAddWidthV;offset_vz++) {
		  const vmesh::GlobalID neighbor_block = get_velocity_block({{indices[0] + offset_vx, indices[1] + offset_vy, indices[2] + offset_vz}});
		  neighbors_have_content.insert(neighbor_block); //add all potential ngbrs of this block with content
	       }
	    }
	 }
      }
      #endif
      
      //add neighbor content info for spatial space neighbors to map. We loop over
      //neighbor cell lists with existing blocks, and raise the
      //flag for the local block with same block id
      for (std::vector<SpatialCell*>::const_iterator neighbor = spatial_neighbors.begin();
	   neighbor != spatial_neighbors.end(); neighbor++ ) {
	 for (vmesh::LocalID block_index=0;block_index< (*neighbor)->velocity_block_with_content_list.size();block_index++){
	    vmesh::GlobalID block = (*neighbor)->velocity_block_with_content_list[block_index];
	    neighbors_have_content.insert(block);
	 }
      }

      // REMOVE all blocks in this cell without content + without neighbors with content
      // better to do it in the reverse order, as then blocks at the
      // end are removed first, and we may avoid copying extra
      // data.
      if (doDeleteEmptyBlocks) {
	 for (int block_index= this->velocity_block_with_no_content_list.size()-1; block_index>=0; --block_index) {
	    const vmesh::GlobalID blockGID = this->velocity_block_with_no_content_list[block_index];
	    const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);

	    bool removeBlock = false;
	    #ifdef AMR
	    if (neighbors_have_content.find(blockGID) != neighbors_have_content.end()) continue;
	    if (neighbors_have_content.find(vmesh.getParent(blockGID)) != neighbors_have_content.end()) continue;
	    
	    std::vector<vmesh::GlobalID> children;
	    vmesh.getChildren(blockGID,children);
	    int counter = 0;
	    for (size_t c=0; c<children.size(); ++c) {
	       if (neighbors_have_content.find(children[c]) != neighbors_have_content.end()) ++counter;
	    }
	    if (counter > 0) continue;
	    removeBlock = true;
	    #else
	    std::unordered_set<vmesh::GlobalID>::iterator it = neighbors_have_content.find(blockGID);
	    if (it == neighbors_have_content.end()) removeBlock = true;
	    #endif
	    if (removeBlock == true) {
	       //No content, and also no neighbor have content -> remove
	       //increment rho loss counters
	       const Real* block_parameters = get_block_parameters(blockLID);
	       const Real DV3 = block_parameters[BlockParams::DVX]
		 * block_parameters[BlockParams::DVY]
		 * block_parameters[BlockParams::DVZ];
	       Real sum=0;
	       for (unsigned int i=0; i<WID3; ++i) sum += get_data(blockLID)[i];
	       this->parameters[CellParams::RHOLOSSADJUST] += DV3*sum;
	       
	       // and finally remove block
	       this->remove_velocity_block(blockGID);
	    }
	 }
      }
    
      #ifdef AMR

      #else
      // ADD all blocks with neighbors in spatial or velocity space (if it exists then the block is unchanged)
      for (boost::unordered_set<vmesh::GlobalID>::iterator it=neighbors_have_content.begin(); it != neighbors_have_content.end(); ++it) {
	 this->add_velocity_block(*it);
      }
      #endif
   }

   void SpatialCell::refine_block(const vmesh::GlobalID& blockGID,std::map<vmesh::GlobalID,vmesh::LocalID>& insertedBlocks) {
      if (blockGID == invalid_global_id()) {
	 std::cerr << "invalid global ID, skip refinement" << std::endl;
	 return;
      }
      
      std::set<vmesh::GlobalID> erasedBlocks;
      std::map<vmesh::GlobalID,vmesh::LocalID> newInserted;
      if (vmesh.refine(blockGID,erasedBlocks,newInserted) == false) {
	 return;
      }
      
      const size_t newBlocks = newInserted.size()-erasedBlocks.size();
      blockContainer.setSize(blockContainer.size() + newBlocks);
      
      for (std::map<vmesh::GlobalID,vmesh::LocalID>::iterator it=newInserted.begin(); it!=newInserted.end(); ++it) {
	 // Set refined block parameters
	 Real* blockParams = blockContainer.getParameters(it->second);
	 vmesh.getBlockCoordinates(it->first,blockParams);
	 vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getCellSize(it->first,blockParams+3);
      }
      
      insertedBlocks.insert(newInserted.begin(),newInserted.end());
   }

} // namespace spatial_cell
