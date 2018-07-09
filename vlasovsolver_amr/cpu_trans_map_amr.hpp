/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef CPU_TRANS_MAP_H
#define CPU_TRANS_MAP_H

#ifndef NDEBUG
   #define DEBUG_VLASOV_SOLVER
#endif

#include "vec.h"
#include "algorithm"
#include "cmath"
#include "utility"
#include "common.h"
#include "spatial_cell.hpp"
#include "cpu_1d_plm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_pqm.hpp"
#include "grid.h"

using namespace std;
using namespace spatial_cell;


// indices in padded block. b_k is the block index in z direction in
// ordinary space (- VLASOV_STENCIL_WIDTH to VLASOV_STENCIL_WIDTH)
//, i,j,k are the cell ids inside on block.
#define i_trans_pblockv(b_k, j, k)  ( (b_k + VLASOV_STENCIL_WIDTH ) + ( (j) + (k) * WID ) * ( 1 + 2 * VLASOV_STENCIL_WIDTH) )

// indices in padded target block, which has Vec4 elements. b_k is the
// block index in z direction in ordinary space, i,j,k are the cell
// ids inside on block (i in vector elements).
#define i_trans_ptblockv(b_k,j,k)  ( (j) + (k) * WID +((b_k) + 1 ) * WID2)

const int PAD=1;
static Realf tempSource[(WID+2*PAD)*(WID+2*PAD)*(WID+2*PAD)];

//Is cell translated? It is not translated if DO_NO_COMPUTE or if it is sysboundary cell and not in first sysboundarylayer
bool do_translate_cell(SpatialCell* SC) {
   if (SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
      (SC->sysBoundaryLayer != 1 && SC->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY))
     return false;
   else
     return true;
}

/*
 * return INVALID_CELLID if the spatial neighbor does not exist, or if
 * it is a cell that is not computed. If the
 * include_first_boundary_layer flag is set, then also first boundary
 * layer is inlcuded (does not return INVALID_CELLID).
 * This does not use dccrg's get_neighbor_of function as it does not support computing neighbors for remote cells

 TODO: not needed anymore as we do not need to compute ngbrs for remote cells
 */
CellID get_spatial_neighbor(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                            const CellID& cellID,
                            const bool include_first_boundary_layer,
                            const int spatial_di,
                            const int spatial_dj,
                            const int spatial_dk ) {
   dccrg::Types<3>::indices_t indices_unsigned = mpiGrid.mapping.get_indices(cellID);
   int64_t indices[3];
   dccrg::Grid_Length::type length = mpiGrid.mapping.length.get();
   
   //boost::array<double, 3> cell_length = mpiGrid.geometry.get_length(cells[i]);
   
   //compute raw new indices
   indices[0] = spatial_di + indices_unsigned[0];
   indices[1] = spatial_dj + indices_unsigned[1];
   indices[2] = spatial_dk + indices_unsigned[2];
   
   //take periodicity into account
   for(uint i = 0; i<3; i++) {
      if(mpiGrid.topology.is_periodic(i)) {
         while(indices[i] < 0 )
            indices[i] += length[i];
         while(indices[i] >= length[i] )
            indices[i] -= length[i];
      }
      
   }
   //return INVALID_CELLID for cells outside system (non-periodic)
   for(uint i = 0; i<3; i++) {
      if(indices[i]< 0)
         return INVALID_CELLID;
      if(indices[i]>=length[i])
         return INVALID_CELLID;
   }

   //store nbr indices into the correct datatype
   for(uint i = 0; i<3; i++) {
      indices_unsigned[i] = indices[i];
   }
   //get nbrID
   CellID nbrID =  mpiGrid.mapping.get_cell_from_indices(indices_unsigned,0);

   if (nbrID == dccrg::error_cell ) {
      std::cerr << __FILE__ << ":" << __LINE__
                << " No neighbor for cell?" << cellID
                << " at offsets " << spatial_di << ", " << spatial_dj << ", " << spatial_dk
                << std::endl;
      
      abort();
   }

   // not existing cell or do not compute       
   if( mpiGrid[nbrID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE)
      return INVALID_CELLID; 

   //cell on boundary, but not first layer and we want to include
   //first layer (e.g. when we compute source cells)
   if( include_first_boundary_layer &&
       mpiGrid[nbrID]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
       mpiGrid[nbrID]->sysBoundaryLayer != 1 ) {
      return INVALID_CELLID;
   }
   //cell on boundary, and we want none of the layers,
   //invalid.(e.g. when we compute targets)
   if( !include_first_boundary_layer &&
       mpiGrid[nbrID]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY){
      return INVALID_CELLID;
   }

   return nbrID; //no AMR
}

template<int DIR> inline
void addUpstreamBlocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,CellID nbrID,
                       int dim,vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh) {
   if (nbrID == INVALID_CELLID) return;
   
   SpatialCell* cellNbr = mpiGrid[nbrID];
   for (vmesh::LocalID blockLID=0; blockLID<cellNbr->get_number_of_velocity_blocks(); ++blockLID) {
      Real* blockParams = cellNbr->get_block_parameters(blockLID);
      
      #ifdef DEBUG_VLASOV_SOLVER
      if (blockParams == NULL) {
         std::cerr << "ERROR, cell " << cellID << " got NULL blockParams in " << __FILE__ << ' ' << __LINE__ << std::endl;
         std::cerr << "\t blockLID=" << blockLID << " nbr=" << cells[0] << std::endl;
         exit(1);
      }
      #endif

      switch (DIR) {
       case -1: {
         Real V = std::max(blockParams[dim],blockParams[dim] + WID*blockParams[BlockParams::DVX+dim]);
         if (V <= 0) continue; }
         break;
       case +1: {
         Real V = std::min(blockParams[dim],blockParams[dim] + WID*blockParams[BlockParams::DVX+dim]);
         if (V >= 0) continue; }
         break;
       default:
         std::cerr << "ERROR, invalid DIR in " << __FILE__ << ' ' << __LINE__ << std::endl; exit(1);
         break;
      }

      vmesh::GlobalID nbrGID = cellNbr->get_velocity_block_global_id(blockLID);
      if (vmesh.getLocalID(nbrGID) != vmesh.invalidLocalID()) {
         // The block exists in this cell
         continue;
      } else if (vmesh.getLocalID(vmesh.getParent(nbrGID)) != vmesh.invalidLocalID()) {
         // Parent block exists in this cell, need to refine
         std::set<vmesh::GlobalID> erased;
         std::map<vmesh::GlobalID,vmesh::LocalID> inserted;
         vmesh.refine(vmesh.getParent(nbrGID),erased,inserted);
      } else if (vmesh.hasChildren(nbrGID) == true) {
         // Children block(s) exist in this cell. The whole octant 
         // may not exist, however, so create the missing blocks.
         std::vector<vmesh::GlobalID> children;
         vmesh.getChildren(nbrGID,children);
         vmesh.push_back(children);
      } else {
         // Block, its parent or none of the children exist in this cell.
         // Need to create the block.
         vmesh.push_back(nbrGID);
      }
   }
}

void createTargetMesh(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,CellID cellID,int dim,
                      bool isRemoteCell) {
   const size_t popID = 0;

   // Get the immediate spatial face neighbors of this cell 
   // in the direction of propagation
   CellID cells[3];
   switch (dim) {
    case 0:
      cells[0] = get_spatial_neighbor(mpiGrid,cellID,true,-1,0,0);
      cells[1] = cellID;
      cells[2] = get_spatial_neighbor(mpiGrid,cellID,true,+1,0,0);
      break;
    case 1:
      cells[0] = get_spatial_neighbor(mpiGrid,cellID,true,0,-1,0);
      cells[1] = cellID;
      cells[2] = get_spatial_neighbor(mpiGrid,cellID,true,0,+1,0);
      break;
    case 2:
      cells[0] = get_spatial_neighbor(mpiGrid,cellID,true,0,0,-1);
      cells[1] = cellID;
      cells[2] = get_spatial_neighbor(mpiGrid,cellID,true,0,0,+1);
      break;
    default:
      std::cerr << "create error" << std::endl;
      exit(1);
      break;
   }

   // Remote (buffered) cells do not consider other remote cells as source cells,
   // i.e., only cells local to this process are translated
   if (isRemoteCell == true) {
      if (mpiGrid.is_local(cells[0]) == false) cells[0] = INVALID_CELLID;
      if (mpiGrid.is_local(cells[2]) == false) cells[2] = INVALID_CELLID;
   }

   SpatialCell* spatial_cell = mpiGrid[cellID];
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = spatial_cell->get_velocity_mesh_temporary();
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks_temporary();

   // At minimum the target mesh will be an identical copy of the existing mesh
   if (isRemoteCell == false) vmesh = spatial_cell->get_velocity_mesh(popID);
   else vmesh.clear();
   
   // Add or refine blocks arriving from the upstream
   addUpstreamBlocks<-1>(mpiGrid,cells[0],dim,vmesh);
   addUpstreamBlocks<+1>(mpiGrid,cells[2],dim,vmesh);

   // Target mesh generated, set block parameters
   blockContainer.setSize(vmesh.size());
   for (size_t b=0; b<vmesh.size(); ++b) {
      vmesh::GlobalID blockGID = vmesh.getGlobalID(b);
      Real* blockParams = blockContainer.getParameters(b);
      blockParams[BlockParams::VXCRD] = spatial_cell->get_velocity_block_vx_min(blockGID);
      blockParams[BlockParams::VYCRD] = spatial_cell->get_velocity_block_vy_min(blockGID);
      blockParams[BlockParams::VZCRD] = spatial_cell->get_velocity_block_vz_min(blockGID);
      vmesh.getCellSize(blockGID,&(blockParams[BlockParams::DVX]));
   }
}

/*compute spatial neighbors for source stencil with a size of 2*
 * VLASOV_STENCIL_WIDTH + 1, cellID at VLASOV_STENCIL_WIDTH. First
 * bondary layer included. Invalid cells are replaced by closest good
 * cells (i.e. boundary condition uses constant extrapolation for the
 * stencil values at boundaries*/ /*
void compute_spatial_source_neighbors(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,
                                      const uint dimension,
                                      CellID *neighbors){
   for (int i = -VLASOV_STENCIL_WIDTH; i <= VLASOV_STENCIL_WIDTH; i++) {
      switch (dimension){
       case 0:
         neighbors[i + VLASOV_STENCIL_WIDTH] = get_spatial_neighbor(mpiGrid, cellID, true, i, 0, 0);
         break;
       case 1:
         neighbors[i + VLASOV_STENCIL_WIDTH] = get_spatial_neighbor(mpiGrid, cellID, true, 0, i, 0);
         break;
       case 2:
         neighbors[i + VLASOV_STENCIL_WIDTH] = get_spatial_neighbor(mpiGrid, cellID, true, 0, 0, i);
         break;             
      }             
   }
   
   CellID last_good_cellID = cellID;
   //loop to neative side and replace all invalid cells with the closest good cell
   for(int i = -1; i>=-VLASOV_STENCIL_WIDTH; i--) {
      if (neighbors[i + VLASOV_STENCIL_WIDTH] == INVALID_CELLID) 
        neighbors[i + VLASOV_STENCIL_WIDTH] = last_good_cellID;
      else
        last_good_cellID = neighbors[i + VLASOV_STENCIL_WIDTH];
   }
   
   last_good_cellID = cellID;
   //loop to positive side and replace all invalid cells with the closest good cell
   for(int i=1; i<=VLASOV_STENCIL_WIDTH; i++) {
      if (neighbors[i + VLASOV_STENCIL_WIDTH] == INVALID_CELLID) 
        neighbors[i + VLASOV_STENCIL_WIDTH] = last_good_cellID;
      else
        last_good_cellID = neighbors[i + VLASOV_STENCIL_WIDTH];
   }
}*/

/*compute spatial target neighbors, stencil has a size of 3. No boundary cells are included*/
void compute_spatial_target_neighbors(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,
                                      const uint dimension,
                                      CellID *neighbors) {
   for (int i=-1; i<=1; ++i) {
      switch (dimension) {
       case 0:
         neighbors[i+1] = get_spatial_neighbor(mpiGrid,cellID,false,i,0,0);
         break;
       case 1:
         neighbors[i+1] = get_spatial_neighbor(mpiGrid,cellID,false,0,i,0);
         break;
       case 2:
         neighbors[i+1] = get_spatial_neighbor(mpiGrid,cellID,false,0,0,i);
         break;             
      }             
   }  
}

/* Copy the fx data to the temporary values array, so that the
 * dimensions are correctly swapped. Also, copy the same block for
 * then neighboring spatial cells (in the dimension). neighbors
 * generated with compute_spatial_neighbors_wboundcond)*/
/*
inline void copy_trans_block_data(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                  const CellID cellID,
                                  const CellID* source_neighbors,
                                  const vmesh::GlobalID blockGID,Vec4* values,int dimension) {
   uint cell_indices_to_id[3]={};
   switch (dimension) {
    case 0:
      // i and k coordinates have been swapped
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
    case 1:
      // j and k coordinates have been swapped
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
    case 2:
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   }
   // Copy volume averages of this block from all spatial cells:
   for (int b=-VLASOV_STENCIL_WIDTH; b<=VLASOV_STENCIL_WIDTH; ++b) {
      const CellID srcCell = source_neighbors[b + VLASOV_STENCIL_WIDTH];

      Realf* block_fx;
      const vmesh::LocalID blockLID = mpiGrid[srcCell]->get_velocity_block_local_id(blockGID);
      if (blockLID == mpiGrid[srcCell]->invalid_local_id()) {
         block_fx = mpiGrid[srcCell]->null_block_fx;
      } else {
         block_fx = mpiGrid[srcCell]->get_fx(blockLID);
      }
      
      // Copy fx table, spatial source_neighbors already taken care of when
      //   creating source_neighbors table. If a normal spatial cell does not
      //   simply have the block, its value will be its null_block which
      //   is fine. This null_block has a value of zero in fx, and that
      //   is thus the velocity space boundary
      for (uint k=0; k<WID; ++k) {
         for (uint j=0; j<WID; ++j) {
            for (uint i=0; i<WID; ++i) {
               const uint cell =
                 i * cell_indices_to_id[0] +
                 j * cell_indices_to_id[1] +
                 k * cell_indices_to_id[2];
               // copy data, when reading data from fx we swap dimensions using cell_indices_to_id
               values[i_trans_pblockv(b,j,k)].insert(i,(Real)block_fx[cell]);
            }
         }
      }
   }
}*/

/*!
  Store values to fx array from the new target data we have computed
  
  For dimension=0  we have rotated data
  i -> k
  j -> j
  k -> i
  For dimension=1   we have rotated data
  i -> i
  j -> k
  k -> j
  
*/ /*
inline void store_trans_block_data(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                   const CellID cellID, const CellID *target_neighbors, 
                                   const vmesh::GlobalID blockGID,
                                   Vec4 * __restrict__ target_values,int dimension) {
   uint cell_indices_to_id[3];
  
   switch (dimension) {
    case 0:
      // i and k coordinates have been swapped
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
    case 1:
      // j and k coordinates have been swapped
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
    case 2:
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
      
    default:
      //same as for dimension 2, mostly here to get rid of compiler warning
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=1;
      cell_indices_to_id[2]=1;
      cerr << "Dimension argument wrong: " << dimension << " at " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
      break;
   }
   
   //Store volume averages in target blocks:
   for (int b=-1; b<=1; ++b) {
      if (target_neighbors[b + 1] == INVALID_CELLID) {
         continue; //do not store to boundary cells or otherwise invalid cells
      }
      SpatialCell* spatial_cell = mpiGrid[target_neighbors[b + 1]];
      const vmesh::LocalID blockLID = spatial_cell->get_velocity_block_local_id(blockGID);
      if (blockLID == spatial_cell->invalid_local_id()) {
         // block does not exist. If so, we do not create it and add stuff to it here.
         // We have already created blocks around blocks with content in
         // spatial sense, so we have no need to create even more blocks here
         // TODO add loss counter
         continue;
      }
      
      Realf* block_data = spatial_cell->get_data(blockLID);
      for (uint k=0; k<WID; ++k) {
         for (uint j=0; j<WID; ++j) {
            for (uint i=0; i<WID; ++i) {
               const uint cell =
                 i * cell_indices_to_id[0] +
                 j * cell_indices_to_id[1] +
                 k * cell_indices_to_id[2];
               //store data, when reading data from  data we swap dimensions using cell_indices_to_id
               block_data[cell] += target_values[i_trans_ptblockv(b,j,k)][i];
            }
         }
      }
   }
}*/

/*
  For local cells that are not boundary cells  block data is copied from data to fx, and data is
  set to zero, if boundary cell then   we copy from data to fx, but do not
  touch data. FOr remote cells fx is already up to data as we receive there.  
*/
/*
bool trans_prepare_block_data(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const CellID cellID){
   bool return_value=false;
   SpatialCell* spatial_cell = mpiGrid[cellID];   
   // if we are on boundary then we do not set the data values to zero as these cells should not be updated
   const bool is_boundary = (spatial_cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY);
   // if the cell is remote, then we do no copy data to the fx table, it should already have been set there
   const bool is_local = mpiGrid.is_local(cellID);
   
   if (is_local && !is_boundary) {
      #pragma omp for nowait
      for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(); cell++) {
         //copy data to fx for solvers, and set data to zero as we will map new values there
         spatial_cell->get_fx()[cell] = spatial_cell->get_data()[cell];
         spatial_cell->get_data()[cell] = 0.0;
         return_value=true;
      }      
   } else if(!is_local && !is_boundary) {
      #pragma omp for nowait
      for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(); cell++) {
         //fx already up to date as we received to fx. data
         //needs to be reset as the updates we collect there will
         //be sent to other processes
         spatial_cell->get_data()[cell] = 0.0;
         return_value=true;
      }
   } else if(is_local && is_boundary) {
      #pragma omp for nowait
      for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(); cell++) {
         //data values are up to date, copy to fx for solvers. Do
         //not reset data as we will not propagate stuff there
         spatial_cell->get_fx()[cell] = spatial_cell->get_data()[cell];
         return_value=true;
      }
   } else if(!is_local && is_boundary) {
      #pragma omp for nowait
      for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(); cell++) {
         //fx already up to date as we received to fx. We copy to data, even if this is not needed...
         spatial_cell->get_data()[cell] = spatial_cell->get_fx()[cell];
         return_value=true;
      }
   }

   return return_value;
}*/

void getTargetArrays(const vmesh::GlobalID targetGID,const int& dim,SpatialCell* targetCell,std::vector<Realf*>& targetBlockData) {
   // Fetch pointers to neighbor block data arrays.
   // There will be either 1 or 8 pointers, depending on 
   // the difference in block refinement levels.
   targetBlockData.clear();
   if (targetCell->get_velocity_block_local_id(targetGID) != targetCell->invalid_global_id()) {
      // Neighbor at same refinement level, just a single pointer
      vmesh::LocalID nbrLID = targetCell->get_velocity_block_local_id(targetGID);
      targetBlockData.push_back(targetCell->get_data(nbrLID));
   } else {
      // Neighbor refined, eight pointers
      std::vector<vmesh::LocalID> nbrBlockLIDs;
      targetCell->get_velocity_block_children_local_ids(targetGID,nbrBlockLIDs);

      /*
      if (nbrBlockLIDs.size() != 8) {
          std::cerr << "ERROR occurred, children.size() " << nbrBlockLIDs.size() << std::endl;
          std::cerr << "targetGID " << targetGID << " ref level " << std::endl;
          exit(1);
      }*/
      if (nbrBlockLIDs.size() != 8) {
          //std::cerr << "error failed to find target block for GID " << targetGID;
          //std::cerr << " r=" << (int) targetCell->get_velocity_block_ref_level(targetGID) << std::endl;
          return;
      }
      
      targetBlockData.resize(8);
      switch (dim) {
       case 0:
         targetBlockData[0] = targetCell->get_data(nbrBlockLIDs[0]);
         targetBlockData[1] = targetCell->get_data(nbrBlockLIDs[1]);
         targetBlockData[2] = targetCell->get_data(nbrBlockLIDs[2]);
         targetBlockData[3] = targetCell->get_data(nbrBlockLIDs[3]);
         targetBlockData[4] = targetCell->get_data(nbrBlockLIDs[4]);
         targetBlockData[5] = targetCell->get_data(nbrBlockLIDs[5]);
         targetBlockData[6] = targetCell->get_data(nbrBlockLIDs[6]);
         targetBlockData[7] = targetCell->get_data(nbrBlockLIDs[7]);
         break;
       case 1:
         targetBlockData[0] = targetCell->get_data(nbrBlockLIDs[0]);
         targetBlockData[1] = targetCell->get_data(nbrBlockLIDs[1]);
         targetBlockData[2] = targetCell->get_data(nbrBlockLIDs[4]);
         targetBlockData[3] = targetCell->get_data(nbrBlockLIDs[5]);
         targetBlockData[4] = targetCell->get_data(nbrBlockLIDs[2]);
         targetBlockData[5] = targetCell->get_data(nbrBlockLIDs[3]);
         targetBlockData[6] = targetCell->get_data(nbrBlockLIDs[6]);
         targetBlockData[7] = targetCell->get_data(nbrBlockLIDs[7]);
         break;
       case 2:
         targetBlockData[0] = targetCell->get_data(nbrBlockLIDs[0]);
         targetBlockData[1] = targetCell->get_data(nbrBlockLIDs[1]);
         targetBlockData[2] = targetCell->get_data(nbrBlockLIDs[2]);
         targetBlockData[3] = targetCell->get_data(nbrBlockLIDs[3]);
         targetBlockData[4] = targetCell->get_data(nbrBlockLIDs[4]);
         targetBlockData[5] = targetCell->get_data(nbrBlockLIDs[5]);
         targetBlockData[6] = targetCell->get_data(nbrBlockLIDs[6]);
         targetBlockData[7] = targetCell->get_data(nbrBlockLIDs[7]);
         break;
       default:
         std::cerr << "ERROR in translation, incorrect dimension in " << __FILE__ << ' ' << __LINE__ << std::endl;
         exit(1);
         break;
      }
   }   
}
/*
template<typename REAL> inline
void depositToNeighbor(std::vector<Realf*>& targetBlockData,const REAL& amount,const int& i,const int& j,const int& k) {
   switch (targetBlockData.size()) {
       case 0:
           break;
       case 1:
//           targetBlockData[0][vblock::index(i,j,k)] += amount;
           targetBlockData[0][vblock::index(k,i,j)] += amount;
           break;
       case 8: {
            int i_trgt,j_trgt,k_trgt;
            int octant = vblock::refIndex(i,j,k,i_trgt,j_trgt,k_trgt);
            //targetBlockData[octant][vblock::index(i_trgt  ,j_trgt  ,k_trgt  )] += amount;
            //targetBlockData[octant][vblock::index(i_trgt+1,j_trgt  ,k_trgt  )] += amount;
            //targetBlockData[octant][vblock::index(i_trgt  ,j_trgt+1,k_trgt  )] += amount;
            //targetBlockData[octant][vblock::index(i_trgt+1,j_trgt+1,k_trgt  )] += amount;
            //targetBlockData[octant][vblock::index(i_trgt  ,j_trgt  ,k_trgt+1)] += amount;
            //targetBlockData[octant][vblock::index(i_trgt+1,j_trgt  ,k_trgt+1)] += amount;
            //targetBlockData[octant][vblock::index(i_trgt  ,j_trgt+1,k_trgt+1)] += amount;
            //targetBlockData[octant][vblock::index(i_trgt+1,j_trgt+1,k_trgt+1)] += amount;
            targetBlockData[octant][vblock::index(k_trgt  ,i_trgt  ,j_trgt  )] += amount;
            targetBlockData[octant][vblock::index(k_trgt  ,i_trgt+1,j_trgt  )] += amount;
            targetBlockData[octant][vblock::index(k_trgt  ,i_trgt  ,j_trgt+1)] += amount;
            targetBlockData[octant][vblock::index(k_trgt  ,i_trgt+1,j_trgt+1)] += amount;
            targetBlockData[octant][vblock::index(k_trgt+1,i_trgt  ,j_trgt  )] += amount;
            targetBlockData[octant][vblock::index(k_trgt+1,i_trgt+1,j_trgt  )] += amount;
            targetBlockData[octant][vblock::index(k_trgt+1,i_trgt  ,j_trgt+1)] += amount;
            targetBlockData[octant][vblock::index(k_trgt+1,i_trgt+1,j_trgt+1)] += amount;
        }
           break;
       default:
           break;
   }
}*/

template<typename REAL> inline
void depositToNeighbor(std::vector<Realf*>& targetBlockData,const REAL* amount,
                       const int& i,const int& j,const int& k) {
    switch (targetBlockData.size()) {
       case 0:
           break;
       case 1:
           //targetBlockData[0][vblock::index(i,j,k)] += amount;
           targetBlockData[0][vblock::index(k,i,j)] += amount[0];
           break;
       case 8: {
            int i_trgt,j_trgt,k_trgt;
            //int octant = vblock::refIndex(i,j,k,i_trgt,j_trgt,k_trgt);
            //targetBlockData[octant][vblock::index(i_trgt  ,j_trgt  ,k_trgt  )] += amount;
            //targetBlockData[octant][vblock::index(i_trgt+1,j_trgt  ,k_trgt  )] += amount;
            //targetBlockData[octant][vblock::index(i_trgt  ,j_trgt+1,k_trgt  )] += amount;
            //targetBlockData[octant][vblock::index(i_trgt+1,j_trgt+1,k_trgt  )] += amount;
            //targetBlockData[octant][vblock::index(i_trgt  ,j_trgt  ,k_trgt+1)] += amount;
            //targetBlockData[octant][vblock::index(i_trgt+1,j_trgt  ,k_trgt+1)] += amount;
            //targetBlockData[octant][vblock::index(i_trgt  ,j_trgt+1,k_trgt+1)] += amount;
            //targetBlockData[octant][vblock::index(i_trgt+1,j_trgt+1,k_trgt+1)] += amount;
            int octant = vblock::refIndex(k,i,j,k_trgt,i_trgt,j_trgt);
            targetBlockData[octant][vblock::index(k_trgt  ,i_trgt  ,j_trgt  )] += amount[0];
            targetBlockData[octant][vblock::index(k_trgt  ,i_trgt+1,j_trgt  )] += amount[1];
            targetBlockData[octant][vblock::index(k_trgt  ,i_trgt  ,j_trgt+1)] += amount[2];
            targetBlockData[octant][vblock::index(k_trgt  ,i_trgt+1,j_trgt+1)] += amount[3];
            targetBlockData[octant][vblock::index(k_trgt+1,i_trgt  ,j_trgt  )] += amount[4];
            targetBlockData[octant][vblock::index(k_trgt+1,i_trgt+1,j_trgt  )] += amount[5];
            targetBlockData[octant][vblock::index(k_trgt+1,i_trgt  ,j_trgt+1)] += amount[6];
            targetBlockData[octant][vblock::index(k_trgt+1,i_trgt+1,j_trgt+1)] += amount[7];
        }
           break;
       default:
           break;
   }
}

/* 
 Here we map from the current time step grid, to a target grid which
 is the lagrangian departure grid (so the grid at timestep +dt,
 tracked backwards by -dt). This is done in ordinary space in the translation step

 This function can, and should be, safely called in a parallel
 OpenMP region (as long as it does only one dimension per parallel
 refion). It is safe as each thread only computes certain blocks (blockID%tnum_threads = thread_num */
bool trans_map_1d(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const CellID cellID,const uint dimension,const Real dt) {
   // Compute target cells (this cell and its face neighbors)
   CellID targetCellIDs[3];
   compute_spatial_target_neighbors(mpiGrid,cellID,dimension,targetCellIDs);

   SpatialCell* targetCells[3];
   for (int i=0; i<3; ++i) targetCells[i] = mpiGrid[targetCellIDs[i]];

   // Get the source mesh (stored in the temporary mesh)
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = targetCells[1]->get_velocity_mesh_temporary();
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = targetCells[1]->get_velocity_blocks_temporary();

   vector<vector<Realf*> > targetBlocks(3);

   for (vmesh::LocalID blockLID=0; blockLID<vmesh.size(); ++blockLID) {
      const vmesh::GlobalID blockGID = vmesh.getGlobalID(blockLID);
      const Real* blockParams = blockContainer.getParameters(blockLID);
      const Real* cellParams  = targetCells[1]->get_cell_parameters();
      Realf* dataSource       = blockContainer.getData(blockLID);
      const Real DZ = cellParams[CellParams::DX+dimension];

      // Note: neighbor may not exist
      getTargetArrays(blockGID,dimension,targetCells[0],targetBlocks[0]);
      getTargetArrays(blockGID,dimension,targetCells[1],targetBlocks[1]);
      getTargetArrays(blockGID,dimension,targetCells[2],targetBlocks[2]);

      bool reconstruct=false;
      if (targetBlocks[0].size() > 1) reconstruct=true;
      if (targetBlocks[1].size() > 1) reconstruct=true;
      if (targetBlocks[2].size() > 1) reconstruct=true;

      if (reconstruct == true) {
          // Load padded data
          SpatialCell::fetch_data<PAD>(blockGID,vmesh,blockContainer.getData(),tempSource);

          for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
              // Calculate reconstruction coefficients
              Real a;
              Real f_lft = tempSource[vblock::padIndex<PAD>(i  ,j+1,k+1)];
              Real f_cen = tempSource[vblock::padIndex<PAD>(i+1,j+1,k+1)];
              Real f_rgt = tempSource[vblock::padIndex<PAD>(i+2,j+1,k+1)];
              //reconstruct_plm(f_lft,f_cen,f_rgt,a);
              a=0;

              Real b;
              f_lft = tempSource[vblock::padIndex<PAD>(i+1,j  ,k+1)];
              f_rgt = tempSource[vblock::padIndex<PAD>(i+1,j+2,k+1)];
              //reconstruct_plm(f_lft,f_cen,f_rgt,b);
              b=0;

              Real c;
              f_lft = tempSource[vblock::padIndex<PAD>(i+1,j+1,k  )];
              f_rgt = tempSource[vblock::padIndex<PAD>(i+1,j+1,k+2)];
              //reconstruct_plm(f_lft,f_cen,f_rgt,c);
              c=0;

              // Reconstructed values in each cell octant
              Real avgsRec[8];
              avgsRec[0] = dataSource[vblock::index(i,j,k)] - 0.25*a - 0.25*b - 0.25*c;
              avgsRec[1] = dataSource[vblock::index(i,j,k)] + 0.25*a - 0.25*b - 0.25*c;
              avgsRec[2] = dataSource[vblock::index(i,j,k)] - 0.25*a + 0.25*b - 0.25*c;
              avgsRec[3] = dataSource[vblock::index(i,j,k)] + 0.25*a + 0.25*b - 0.25*c;
              avgsRec[4] = dataSource[vblock::index(i,j,k)] - 0.25*a - 0.25*b + 0.25*c;
              avgsRec[5] = dataSource[vblock::index(i,j,k)] + 0.25*a - 0.25*b + 0.25*c;
              avgsRec[6] = dataSource[vblock::index(i,j,k)] - 0.25*a + 0.25*b + 0.25*c;
              avgsRec[7] = dataSource[vblock::index(i,j,k)] + 0.25*a + 0.25*b + 0.25*c;
/*
              Real V_bot = (blockParams[dimension] + 0.25*blockParams[BlockParams::DVX+dimension])*dt / DZ;
              if (V_bot < 0) depositToNeighbor(targetBlocks[0],avgsRec,-V_bot,i,j,k);
              else           depositToNeighbor(targetBlocks[2],avgsRec, V_bot,i,j,k);

              Real V_top = (blockParams[dimension] + 0.75*blockParams[BlockParams::DVX+dimension])*dt / DZ;
              if (V_top < 0) depositToNeighbor(targetBlocks[0],avgsRec,-V_top,i,j,k);
              else           depositToNeighbor(targetBlocks[2],avgsRec, V_top,i,j,k);

              #warning This needs to be transposed
              V_bot = 1 - fabs(V_bot);
              V_top = 1 - fabs(V_top);
              avgsRec[0] *= V_bot;
              avgsRec[1] *= V_top;
              avgsRec[2] *= V_bot;
              avgsRec[3] *= V_top;
              avgsRec[4] *= V_bot;
              avgsRec[5] *= V_top;
              avgsRec[6] *= V_bot;
              avgsRec[7] *= V_top;
*/
              /*
              switch (targetBlocks[1].size()) {
                  case 0:
                      break;
                  case 1:
                      targetBlocks[1][0][vblock::index(i,j,k)] += dataSource[vblock::index(i,j,k)];
                      break;
                  case 8: {
                      int i_trgt,j_trgt,k_trgt;
                      int octant = vblock::refIndex(i,j,k,i_trgt,j_trgt,k_trgt);
                      targetBlocks[1][octant][vblock::index(i_trgt  ,j_trgt  ,k_trgt  )] += avgsRec[0];
                      targetBlocks[1][octant][vblock::index(i_trgt+1,j_trgt  ,k_trgt  )] += avgsRec[1];
                      targetBlocks[1][octant][vblock::index(i_trgt  ,j_trgt+1,k_trgt  )] += avgsRec[2];
                      targetBlocks[1][octant][vblock::index(i_trgt+1,j_trgt+1,k_trgt  )] += avgsRec[3];
                      targetBlocks[1][octant][vblock::index(i_trgt  ,j_trgt  ,k_trgt+1)] += avgsRec[4];
                      targetBlocks[1][octant][vblock::index(i_trgt+1,j_trgt  ,k_trgt+1)] += avgsRec[5];
                      targetBlocks[1][octant][vblock::index(i_trgt  ,j_trgt+1,k_trgt+1)] += avgsRec[6];
                      targetBlocks[1][octant][vblock::index(i_trgt+1,j_trgt+1,k_trgt+1)] += avgsRec[7];
                      }
                      break;
                  default:
                      break;                      
              }*/
              
              depositToNeighbor(targetBlocks[1],avgsRec,i,j,k);
          }

          continue;
      }

      for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
         // Target block can be at the same refinement level or at +1 refinement level.
         // Note that the target block in this cell can also be at higher refinement level.
         //Real V_bot = (blockParams[dimension] + 0.25*blockParams[BlockParams::DVX+dimension])*dt / DZ;
         //Real V_top = (blockParams[dimension] + 0.75*blockParams[BlockParams::DVX+dimension])*dt / DZ;

         int k_trgt;
         Real removedMass = 0.0;
         //Real V_norm_min,V_norm_max;
         /*
         if (V_bot < 0) {
            V_norm_min = 0.0;
            V_norm_max = -V_bot;

            // This should be integrated using reconstructed values
            Real amount = 0.5*dataSource[vblock::index(i,j,k)]*(V_norm_max-V_norm_min);

            depositToNeighbor(targetBlocks[0],amount,i,j,k);
            removedMass += amount;
         } else {
            V_norm_min = 1.0 - V_top;
            V_norm_max = 1.0;

            // This should be integrated using reconstructed values
            Real amount = 0.5*dataSource[vblock::index(i,j,k)]*(V_norm_max-V_norm_min);

            depositToNeighbor(targetBlocks[2],amount,i,j,k);
            removedMass += amount;
         }*/

         const Real amount = dataSource[vblock::index(i,j,k)]-removedMass;
         //depositToNeighbor(targetBlocks[1],amount,i,j,k);
         targetBlocks[1][0][vblock::index(k,i,j)] += amount;
      }
   }

   return true;
}

/*!

  This function communicates the mapping on process boundaries, and then updates the data to their correct values.
  TODO, this could be inside an openmp region, in which case some m ore barriers and masters should be added

  \par dimension: 0,1,2 for x,y,z
  \par direction: 1 for + dir, -1 for - dir
*/
/*
void update_remote_mapping_contribution(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const uint dimension, int direction) {
   const vector<CellID> local_cells = mpiGrid.get_cells();
   const vector<CellID> remote_cells = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   vector<CellID> receive_cells;
   vector<CellID> send_cells;
   
   //normalize
   if(direction > 0)
      direction = 1;
   if(direction < 0)
      direction = -1;

   for (size_t c=0; c<remote_cells.size(); ++c) {
      SpatialCell *ccell = mpiGrid[remote_cells[c]];
      //default values, to avoid any extra sends and receives
      ccell->neighbor_block_data = &(ccell->get_data()[0]);
      ccell->neighbor_number_of_blocks = 0;
   }
   
   //prepare arrays
   for (size_t c=0; c<local_cells.size(); ++c) {
      SpatialCell *ccell = mpiGrid[local_cells[c]];
      //default values, to avoid any extra sends and receives
      ccell->neighbor_block_data = &(ccell->get_data()[0]);
      ccell->neighbor_number_of_blocks = 0;
      CellID p_ngbr,m_ngbr;
      
      switch(dimension) {
       case 0:
         p_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], false,  direction, 0, 0); //p_ngbr is target, if in boundaries then it is not updated
         m_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], true, -direction, 0, 0); //m_ngbr is source, first boundary layer is propagated so that it flows into system
         break;
       case 1:
         p_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], false, 0, direction, 0); //p_ngbr is target, if in boundaries then it is not update
         m_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], true, 0, -direction, 0); //m_ngbr is source, first boundary layer is propagated so that it flows into system
         break;
       case 2:
         p_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], false, 0, 0, direction); //p_ngbr is target, if in boundaries then it is not update
         m_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], true,  0, 0, -direction); //m_ngbr is source, first boundary layer is propagated so that it flows into system
         break;
       default:
         cerr << "Dimension wrong at (impossible!) "<< __FILE__ <<":" << __LINE__<<endl;
         exit(1);
      }

      if (mpiGrid.is_local(p_ngbr) && mpiGrid.is_local(m_ngbr)) continue; //internal cell, not much to do
            
      SpatialCell* pcell = NULL;
      if (p_ngbr != INVALID_CELLID) pcell = mpiGrid[p_ngbr];
      SpatialCell* mcell = NULL;
      if (m_ngbr != INVALID_CELLID) mcell = mpiGrid[m_ngbr];

      if (p_ngbr != INVALID_CELLID && !mpiGrid.is_local(p_ngbr) && do_translate_cell(ccell)) {
         //Send data in p_ngbr data array that we just
         //mapped to if 1) it is a valid target,
         //2) is remote cell, 3) if the source cell in center was
         //translated
         ccell->neighbor_block_data = &(pcell->get_data()[0]);
         ccell->neighbor_number_of_blocks = pcell->get_number_of_velocity_blocks();
         send_cells.push_back(p_ngbr);
      }
      
      if (m_ngbr != INVALID_CELLID && !mpiGrid.is_local(m_ngbr) && ccell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         //Receive data that mcell mapped to ccell to this local cell
         //fx array, if 1) m is a valid source cell, 2) center cell is to be updated (normal cell) 3)  m is remote
         mcell->neighbor_block_data = &(ccell->get_fx()[0]);
         mcell->neighbor_number_of_blocks = ccell->get_number_of_velocity_blocks();
         receive_cells.push_back(local_cells[c]);
      }

   }
   //Do communication
   SpatialCell::set_mpi_transfer_type(Transfer::NEIGHBOR_VEL_BLOCK_FLUXES);

   switch (dimension) {
    case 0:
      if (direction > 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_P_X_NEIGHBORHOOD_ID);  
      if (direction < 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_M_X_NEIGHBORHOOD_ID);  
      break;
    case 1:
      if (direction > 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_P_Y_NEIGHBORHOOD_ID);  
      if (direction < 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_M_Y_NEIGHBORHOOD_ID);  
      break;
    case 2:
      if (direction > 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_P_Z_NEIGHBORHOOD_ID);  
      if (direction < 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_M_Z_NEIGHBORHOOD_ID);  
      break;
   }

   #pragma omp parallel
   {
      //reduce data: sum received fx to data
      for (size_t c=0; c < receive_cells.size(); ++c) {
         SpatialCell *spatial_cell = mpiGrid[receive_cells[c]];      
         #pragma omp for nowait
         for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(); cell++) {
            //copy data to fx for solvers, and set data to zero as we will map new values there
            spatial_cell->get_data()[cell] += spatial_cell->get_fx()[cell];
         }
      }

      // send cell data is set to zero. This is to avoid double copy if
      // one cell is the neighbor on bot + and - side to the same
      // process
      for (size_t c=0; c < send_cells.size(); ++c) {
         SpatialCell *spatial_cell = mpiGrid[send_cells[c]];      
         #pragma omp for nowait
         for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(); cell++) {
            spatial_cell->get_data()[cell] = 0.0;
         }
      }
   }
}
*/ 

#endif   
