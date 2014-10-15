/*!
Space for static variables of spatial cell class for Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute


*/

#include <unordered_set>
#include <vectorclass.h>

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
    * modify local cell.
    * 
    * NOTE: The AMR mesh must be valid, otherwise this function will
    * remove some blocks that should not be removed.*/
   void SpatialCell::adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors, bool doDeleteEmptyBlocks) {
      #ifdef AMR
//         return;
      #endif

      //  This set contains all those cellids which have neighbors in any
      //  of the 6-dimensions Actually, we would only need to add
      //  local blocks with no content here, as blocks with content
      //  do not need to be created and also will not be removed as
      //  we only check for removal for blocks with no content
      std::unordered_set<vmesh::GlobalID> neighbors_have_content;

      #ifdef AMR
      for (vmesh::LocalID block_index=0; block_index<velocity_block_with_content_list.size(); ++block_index) {
         vmesh::GlobalID blockGID = velocity_block_with_content_list[block_index];
         vector<vmesh::GlobalID> neighborGIDs;
         vmesh.getNeighborsExistingAtSameLevel(blockGID,neighborGIDs);
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
            const vmesh::GlobalID blockGID = this->velocity_block_with_no_content_list[block_index];
            #ifdef DEBUG_SPATIAL_CELL
               if (blockGID == invalid_global_id())
                  cerr << "Got invalid block at " << __FILE__ << ' ' << __LINE__ << endl; exit(1);               
            #endif
            const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);
            #ifdef DEBUG_SPATIAL_CELL
               if (blockLID == invalid_local_id())
                  cerr << "Could not find block in " << __FILE__ << ' ' << __LINE__ << endl; exit(1);               
            #endif
            
            bool removeBlock = false;
            #ifdef AMR
            // Check this block in the neighbor cells
            if (neighbors_have_content.find(blockGID) != neighbors_have_content.end()) continue;
            
            // Check the parent of this block in the neighbor cells
            if (neighbors_have_content.find(vmesh.getParent(blockGID)) != neighbors_have_content.end()) continue;
            
            // Check all the children of this block in the neighbor cells
            std::vector<vmesh::GlobalID> children;
            vmesh.getChildren(blockGID,children);
            int counter = 0;
            for (size_t c=0; c<children.size(); ++c) {
               if (neighbors_have_content.find(children[c]) != neighbors_have_content.end()) ++counter;
            }
            if (counter > 0) continue;
            
            // It is safe to remove this block
            removeBlock = true;
            #else
            std::unordered_set<vmesh::GlobalID>::iterator it = neighbors_have_content.find(blockGID);
            if (it == neighbors_have_content.end()) removeBlock = true;
            #endif

            if (removeBlock == true) {
               //No content, and also no neighbor have content -> remove
               //and increment rho loss counters
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
      for (std::unordered_set<vmesh::GlobalID>::iterator it=neighbors_have_content.begin(); it != neighbors_have_content.end(); ++it) {
         this->add_velocity_block(*it);
      }
      #endif
   }

   void SpatialCell::coarsen_block(const vmesh::GlobalID& parent,const std::vector<vmesh::GlobalID>& children) {
      // First create the parent (coarse) block and grab pointer to its data.
      // add_velocity_block initializes data to zero values.
      if (add_velocity_block(parent) == false) return;
      vmesh::LocalID parentLID = vmesh.getLocalID(parent);
      Realf* parent_data = get_data(parentLID);

      // Calculate children (fine) block local IDs, some of the children may not exist
      for (size_t c=0; c<children.size(); ++c) {
         vmesh::LocalID childrenLID = vmesh.getLocalID(children[c]);
         if (childrenLID == vmesh.invalidLocalID()) continue;
         Realf* data = get_data(childrenLID);

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
         remove_velocity_block(children[c]);
      }
   }

   void SpatialCell::coarsen_blocks(amr_ref_criteria::Base* refCriterion) {
      // Sort blocks according to their refinement levels
      vector<vector<vmesh::GlobalID> > blocks(vmesh.getMaxAllowedRefinementLevel()+1);

      for (vmesh::LocalID blockLID=0; blockLID<get_number_of_velocity_blocks(); ++blockLID) {
         vmesh::GlobalID blockGID = vmesh.getGlobalID(blockLID);
         uint8_t r = vmesh.getRefinementLevel(blockGID);
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
            fetch_data<1>(blockGID,vmesh,get_data(),array);
            if (refCriterion->evaluate(array) < Parameters::amrCoarsenLimit) coarsenList.insert(blockGID);
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
            vmesh.getSiblings(*it,siblings);
            bool allows=true;
            for (size_t s=0; s<siblings.size(); ++s) {
               // Skip non-existing blocks
               if (vmesh.getLocalID(siblings[s]) == vmesh.invalidLocalID()) {
                  
                  continue;
               }

               // Check that the sibling allows coarsening
               if (coarsenList.find(siblings[s]) == coarsenList.end()) {
                  allows = false;
                  break;
               }
            }

            // Check that the mesh structure allows coarsening
            if (vmesh.coarsenAllowed(*it) == false) continue;

            // If all siblings allow coarsening, add block to coarsen list
            if (allows == true) {
               allowCoarsen.insert(make_pair(vmesh.getParent(*it),siblings));
            }
         }

         cerr << "ref level " << r << " has " << allowCoarsen.size() << " blocks for coarsening" << endl;
         for (unordered_map<vmesh::GlobalID,vector<vmesh::GlobalID> >::const_iterator it=allowCoarsen.begin(); it!=allowCoarsen.end(); ++it) {
            coarsen_block(it->first,it->second);
         }
      }
   }

   void SpatialCell::merge_values_recursive(vmesh::GlobalID parentGID,vmesh::GlobalID blockGID,uint8_t refLevel,bool recursive,const Realf* data,
					    std::set<vmesh::GlobalID>& blockRemovalList) {
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == invalid_global_id()) {
         cerr << "merge_values_recursive called for GID=" << blockGID << " at r=" << (int)refLevel << " parent=" << parentGID << endl;
      }
      #endif

      // Get all possible children:
      vector<vmesh::GlobalID> childrenGIDs;
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getChildren(blockGID,childrenGIDs);

      // Check if any of block's children exist:
      bool hasChildren = false;
      for (size_t c=0; c<childrenGIDs.size(); ++c) {
         if (vmesh.getLocalID(childrenGIDs[c]) != vmesh.invalidLocalID()) {
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
         vmesh::LocalID blockLID = vmesh.getLocalID(blockGID);
	 
         #ifdef DEBUG_SPATIAL_CELL
         if (blockLID == invalid_local_id()) {
            cerr << "ERROR: Failed to merge values, block does not exist!" << endl;
            cerr << "\t exiting from " << __FILE__ << ':' << __LINE__ << endl;
            exit(1);
         }
         #endif

         Realf* myData = blockContainer.getData(blockLID);
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
         merge_values_recursive(blockGID,childrenGIDs[octant],refLevel+1,true,array,blockRemovalList);
      }

      // Data merged to children, block can be removed
      blockRemovalList.insert(blockGID);
   }

   void SpatialCell::merge_values() {
      const uint8_t maxRefLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getMaxAllowedRefinementLevel();

      for (int i=0; i<WID3; ++i) null_block_data[i] = 0;
      
      // Sort blocks according to their refinement levels:
      vector<vector<vmesh::GlobalID> > blocks(maxRefLevel+1);
      for (vmesh::LocalID blockLID=0; blockLID<get_number_of_velocity_blocks(); ++blockLID) {
         const vmesh::GlobalID blockGID = vmesh.getGlobalID(blockLID);
         
         if (blockGID == vmesh.invalidGlobalID()) {
            cerr << "got invalid global id from mesh!" << endl;
            continue;
         }

         uint8_t refLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getRefinementLevel(blockGID);
         blocks[refLevel].push_back(blockGID);
      }
      
      set<vmesh::GlobalID> blockRemovalList;
      for (uint8_t refLevel=0; refLevel<blocks.size()-1; ++refLevel) {
         for (size_t b=0; b<blocks[refLevel].size(); ++b) {
            const vmesh::GlobalID blockGID = blocks[refLevel][b];
            const vmesh::LocalID  blockLID = vmesh.getLocalID(blockGID);
            if (blockLID == vmesh.invalidLocalID()) continue;
            
            const Realf* data = blockContainer.getData(blockLID);
            merge_values_recursive(blockGID,blockGID,refLevel,true,data,blockRemovalList);
         }
      }
      
      cerr << "should remove " << blockRemovalList.size() << " blocks" << endl;
      for (set<vmesh::GlobalID>::const_iterator it=blockRemovalList.begin(); it!=blockRemovalList.end(); ++it) {
         //remove_velocity_block(*it);
      }
   }
   
   void SpatialCell::add_values(const vmesh::GlobalID& targetGID,
                                std::unordered_map<vmesh::GlobalID,Realf[(WID+2)*(WID+2)*(WID+2)]>& sourceData) {
      vmesh::LocalID targetLID = get_velocity_block_local_id( targetGID );
      if (targetLID == invalid_local_id()) {
         std::cerr << "error has occurred" << std::endl;
         return;
      }

      Realf* targetData = get_data(targetLID);

      // Add data from all same level blocks
      vector<vmesh::GlobalID> neighborIDs;
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getNeighborsAtSameLevel(targetGID,neighborIDs);

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
      if (vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getRefinementLevel(targetGID) == 0) {
         return;
      }
      const int octant = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getOctant(targetGID);
      
      int parentIndex[3];
      parentIndex[2] = 1 + 2*(octant / 4);
      parentIndex[1] = 1 + 2*((octant - 4*(octant/4))/2);
      parentIndex[0] = 1 + 2*(octant % 2);
      
      // Add data from all coarser blocks
      const Realf face_mul = 2.0;
      const Realf edge_mul = 4.0;
      const Realf corn_mul = 8.0;

      const vmesh::GlobalID targetParentGID = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getParent( neighborIDs[vblock::nbrIndex( 0, 0, 0)] );
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
         parentGID = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getParent( neighborIDs[vblock::nbrIndex(i_off,0,0)] );
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
         parentGID = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getParent( neighborIDs[vblock::nbrIndex(0,j_off,0)] );
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
         parentGID = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getParent( neighborIDs[vblock::nbrIndex(0,0,k_off)] );
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
      parentGID = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getParent( neighborIDs[vblock::nbrIndex(i_off,j_off,0)] );
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
      parentGID = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getParent( neighborIDs[vblock::nbrIndex(i_off,0,k_off)] );
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
      parentGID = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getParent( neighborIDs[vblock::nbrIndex(0,j_off,k_off)] );
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
      parentGID = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getParent( neighborIDs[vblock::nbrIndex(i_off,j_off,k_off)] );
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

   void SpatialCell::refine_block(const vmesh::GlobalID& blockGID,std::map<vmesh::GlobalID,vmesh::LocalID>& insertedBlocks) {
      if (blockGID == invalid_global_id()) {
         std::cerr << "invalid global ID, skip refinement" << std::endl;
         return;
      }

      // Tell mesh to refine the given block. In return we get the erased 
      // and inserted blocks. Note that multiple blocks can be removed 
      // (and inserted) due to induced refinement. There are eight entries 
      // in newInserted (the children) for each entry in erasedBlocks.
      std::set<vmesh::GlobalID> erasedBlocks;
      std::map<vmesh::GlobalID,vmesh::LocalID> newInserted;
      if (vmesh.refine(blockGID,erasedBlocks,newInserted) == false) {
         return;
      }

      // Resize the block container, this preserves old data.
      const size_t newBlocks = newInserted.size()-erasedBlocks.size();
      blockContainer.setSize(blockContainer.size() + newBlocks);

      std::map<vmesh::GlobalID,vmesh::LocalID>::const_iterator ins=newInserted.begin();
      for (std::set<vmesh::GlobalID>::const_iterator er=erasedBlocks.begin(); er!=erasedBlocks.end(); ++er) {
         for (int child=0; child<8; ++child) {
            // Copy / interpolate data from old (coarse) block to new refined blocks.
            
            
            // Set refined block parameters
            Real* blockParams = blockContainer.getParameters(ins->second);
            vmesh.getBlockCoordinates(ins->first,blockParams);
            vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getCellSize(ins->first,blockParams+3);
            
            ++ins;
         }
      }
      
      for (std::map<vmesh::GlobalID,vmesh::LocalID>::iterator it=newInserted.begin(); it!=newInserted.end(); ++it) {
         // Set refined block parameters
         Real* blockParams = blockContainer.getParameters(it->second);
         vmesh.getBlockCoordinates(it->first,blockParams);
         vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getCellSize(it->first,blockParams+3);
         
         Realf* fx = blockContainer.getFx(it->second);
         for (int i=0; i<WID3; ++i) fx[i] = 0;
      }

      insertedBlocks.insert(newInserted.begin(),newInserted.end());
   }

} // namespace spatial_cell
