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

#include <algorithm>
#include <cmath>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../grid.h"
#include "../object_wrapper.h"
#include "vec.h"
#include "cpu_1d_plm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_ppm_nonuniform.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_trans_map.hpp"

using namespace std;
using namespace spatial_cell;

void compute_spatial_source_neighbors(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,const uint dimension,SpatialCell **neighbors);
void compute_spatial_target_neighbors(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,const uint dimension,SpatialCell **neighbors);
void copy_trans_block_data(SpatialCell** source_neighbors,const vmesh::GlobalID blockGID,
                           Vec* values,const unsigned char* const cellid_transpose,const uint popID);
CellID get_spatial_neighbor(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                            const CellID& cellID,const bool include_first_boundary_layer,
                            const int spatial_di,const int spatial_dj,const int spatial_dk);
SpatialCell* get_spatial_neighbor_pointer(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                          const CellID& cellID,const bool include_first_boundary_layer,
                                          const int spatial_di,const int spatial_dj,const int spatial_dk);
void store_trans_block_data(SpatialCell** target_neighbors,const vmesh::GlobalID blockGID,
                            Vec* __restrict__ target_values,
                            const unsigned char* const cellid_transpose,const uint popID);


// indices in padded source block, which is of type Vec with VECL
// element sin each vector. b_k is the block index in z direction in
// ordinary space [- VLASOV_STENCIL_WIDTH to VLASOV_STENCIL_WIDTH],
// i,j,k are the cell ids inside on block (i in vector elements).
// Vectors with same i,j,k coordinates, but in different spatial cells, are consequtive
//#define i_trans_ps_blockv(j, k, b_k)  ( (b_k + VLASOV_STENCIL_WIDTH ) + ( (((j) * WID + (k) * WID2)/VECL)  * ( 1 + 2 * VLASOV_STENCIL_WIDTH) ) )
#define i_trans_ps_blockv(planeVectorIndex, planeIndex, blockIndex) ( (blockIndex) + VLASOV_STENCIL_WIDTH  +  ( (planeVectorIndex) + (planeIndex) * VEC_PER_PLANE ) * ( 1 + 2 * VLASOV_STENCIL_WIDTH)  )

// indices in padded target block, which is of type Vec with VECL
// element sin each vector. b_k is the block index in z direction in
// ordinary space, i,j,k are the cell ids inside on block (i in vector
// elements).
//#define i_trans_pt_blockv(j, k, b_k) ( ( (j) * WID + (k) * WID2 + ((b_k) + 1 ) * WID3) / VECL )
#define i_trans_pt_blockv(planeVectorIndex, planeIndex, blockIndex)  ( planeVectorIndex + planeIndex * VEC_PER_PLANE + (blockIndex + 1) * VEC_PER_BLOCK)

#define i_trans_ps_blockv_pencil(planeVectorIndex, planeIndex, blockIndex, lengthOfPencil) ( (blockIndex) + VLASOV_STENCIL_WIDTH  +  ( (planeVectorIndex) + (planeIndex) * VEC_PER_PLANE ) * ( lengthOfPencil + 2 * VLASOV_STENCIL_WIDTH) )

//Is cell translated? It is not translated if DO_NO_COMPUTE or if it is sysboundary cell and not in first sysboundarylayer
bool do_translate_cell(SpatialCell* SC){
   if(SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
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
   CellID nbrID = mpiGrid.mapping.get_cell_from_indices(indices_unsigned,0);
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
                                      

/*
 * return NULL if the spatial neighbor does not exist, or if
 * it is a cell that is not computed. If the
 * include_first_boundary_layer flag is set, then also first boundary
 * layer is inlcuded (does not return INVALID_CELLID).
 * This does not use dccrg's get_neighbor_of function as it does not support computing neighbors for remote cells


 */

SpatialCell* get_spatial_neighbor_pointer(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                          const CellID& cellID,
                                          const bool include_first_boundary_layer,
                                          const int spatial_di,
                                          const int spatial_dj,
                                          const int spatial_dk ) {
   CellID nbrID=get_spatial_neighbor(mpiGrid, cellID, include_first_boundary_layer, spatial_di, spatial_dj, spatial_dk);

   if(nbrID!=INVALID_CELLID)
      return mpiGrid[nbrID];
   else
      return NULL;
}

/*compute spatial neighbors for source stencil with a size of 2*
 * VLASOV_STENCIL_WIDTH + 1, cellID at VLASOV_STENCIL_WIDTH. First
 * bondary layer included. Invalid cells are replaced by closest good
 * cells (i.e. boundary condition uses constant extrapolation for the
 * stencil values at boundaries*/
  
void compute_spatial_source_neighbors(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,
                                      const uint dimension,
                                      SpatialCell **neighbors){
   for(int i = -VLASOV_STENCIL_WIDTH; i <= VLASOV_STENCIL_WIDTH; i++){
      switch (dimension){
      case 0:
         neighbors[i + VLASOV_STENCIL_WIDTH] = get_spatial_neighbor_pointer(mpiGrid, cellID, true, i, 0, 0);
         break;
      case 1:
         neighbors[i + VLASOV_STENCIL_WIDTH] = get_spatial_neighbor_pointer(mpiGrid, cellID, true, 0, i, 0);
         break;
      case 2:
         neighbors[i + VLASOV_STENCIL_WIDTH] = get_spatial_neighbor_pointer(mpiGrid, cellID, true, 0, 0, i);
         break;             
      }             
   }

   SpatialCell* last_good_cell = mpiGrid[cellID];
   /*loop to neative side and replace all invalid cells with the closest good cell*/
   for(int i = -1;i>=-VLASOV_STENCIL_WIDTH;i--){
      if(neighbors[i + VLASOV_STENCIL_WIDTH] == NULL) 
         neighbors[i + VLASOV_STENCIL_WIDTH] = last_good_cell;
      else
         last_good_cell = neighbors[i + VLASOV_STENCIL_WIDTH];
   }

   last_good_cell = mpiGrid[cellID];
   /*loop to positive side and replace all invalid cells with the closest good cell*/
   for(int i = 1; i <= VLASOV_STENCIL_WIDTH; i++){
      if(neighbors[i + VLASOV_STENCIL_WIDTH] == NULL) 
         neighbors[i + VLASOV_STENCIL_WIDTH] = last_good_cell;
      else
         last_good_cell = neighbors[i + VLASOV_STENCIL_WIDTH];
   }
}

void compute_spatial_source_cells_for_pencil(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                             setOfPencils pencils,
                                             const uint iPencil,
                                             const uint dimension,
                                             SpatialCell **sourceCells){

   // L = length of the pencil iPencil
   int L = pencils.lengthOfPencils[iPencil];
   vector<CellID> ids = pencils.getIds(iPencil);
      
   for (int iCell = -VLASOV_STENCIL_WIDTH; iCell < L + VLASOV_STENCIL_WIDTH; iCell++) {
      CellID cellID = ids[min(max(iCell, 0), L - 1)];
      
      int i = 0;
      if(iCell <= 0) i = iCell;
      if(iCell >= L) i = iCell - (L - 1);
                                                       
      switch (dimension) {
      case 0:
         sourceCells[iCell + VLASOV_STENCIL_WIDTH] = get_spatial_neighbor_pointer(mpiGrid, cellID, false, i, 0, 0);
         break;
      case 1:
         sourceCells[iCell + VLASOV_STENCIL_WIDTH] = get_spatial_neighbor_pointer(mpiGrid, cellID, false, 0, i, 0);
         break;
      case 2:
         sourceCells[iCell + VLASOV_STENCIL_WIDTH] = get_spatial_neighbor_pointer(mpiGrid, cellID, false, 0, 0, i);
         break;
      }
   }

   SpatialCell* last_good_cell = mpiGrid[ids.front()];
   /*loop to neative side and replace all invalid cells with the closest good cell*/
   for(int i = -1;i>=-VLASOV_STENCIL_WIDTH;i--){
      if(sourceCells[i + VLASOV_STENCIL_WIDTH] == NULL) 
         sourceCells[i + VLASOV_STENCIL_WIDTH] = last_good_cell;
      else
         last_good_cell = sourceCells[i + VLASOV_STENCIL_WIDTH];
   }

   last_good_cell = mpiGrid[ids.back()];
   /*loop to positive side and replace all invalid cells with the closest good cell*/
   for(int i = L + 1; i <= L + VLASOV_STENCIL_WIDTH; i++){
      if(sourceCells[i + VLASOV_STENCIL_WIDTH] == NULL) 
         sourceCells[i + VLASOV_STENCIL_WIDTH] = last_good_cell;
      else
         last_good_cell = sourceCells[i + VLASOV_STENCIL_WIDTH];
   }
}


/*compute spatial target neighbors, stencil has a size of 3. No boundary cells are included*/
void compute_spatial_target_neighbors(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,
                                      const uint dimension,
                                      SpatialCell **neighbors){

   for(int i = -1; i <= 1; i++){
      switch (dimension){
      case 0:
         neighbors[i + 1] = get_spatial_neighbor_pointer(mpiGrid, cellID, false, i, 0, 0);
         break;
      case 1:
         neighbors[i + 1] = get_spatial_neighbor_pointer(mpiGrid, cellID, false, 0, i, 0);
         break;
      case 2:
         neighbors[i + 1] = get_spatial_neighbor_pointer(mpiGrid, cellID, false, 0, 0, i);
         break;             
      }             
   }

}

/*compute spatial target neighbors for pencils of size N. No boundary cells are included*/
void compute_spatial_target_cells_for_pencils(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                              setOfPencils& pencils,
                                              const uint dimension,
                                              SpatialCell **targetCells){

   uint GID = 0;
   for(uint iPencil = 0; iPencil < pencils.N; iPencil++){
      // L = length of the pencil iPencil
      int L = pencils.lengthOfPencils[iPencil];

      vector<CellID> ids = pencils.getIds(iPencil);
      for (int iCell = -1; iCell <= L; iCell++) {
         CellID cellID = ids[min(max(iCell,0),L - 1)];
         
         int i = 0;
         if(iCell == -1) i = -1;
         if(iCell == L)  i = 1;
                                                          
         switch (dimension) {
         case 0:
            targetCells[GID + iCell + 1] = get_spatial_neighbor_pointer(mpiGrid, cellID, false, i, 0, 0);
            break;
         case 1:
            targetCells[GID + iCell + 1] = get_spatial_neighbor_pointer(mpiGrid, cellID, false, 0, i, 0);
            break;
         case 2:
            targetCells[GID + iCell + 1] = get_spatial_neighbor_pointer(mpiGrid, cellID, false, 0, 0, i);
            break;
         }
      }
      GID += (L + 2);
   }
}


/* Copy the data to the temporary values array, so that the
 * dimensions are correctly swapped. Also, copy the same block for
 * then neighboring spatial cells (in the dimension). neighbors
 * generated with compute_spatial_neighbors_wboundcond).
 * 
 * This function must be thread-safe.
 *
 * @param source_neighbors Array containing the VLASOV_STENCIL_WIDTH closest 
 * spatial neighbors of this cell in the propagated dimension.
 * @param blockGID Global ID of the velocity block.
 * @param values Vector where loaded data is stored.
 * @param cellid_transpose
 * @param popID ID of the particle species.
 */
void copy_trans_block_data(
    SpatialCell** source_neighbors,
    const vmesh::GlobalID blockGID,
    Vec* values,
    const unsigned char* const cellid_transpose,
    const uint popID) { 

   /*load pointers to blocks and prefetch them to L1*/
   Realf* blockDatas[VLASOV_STENCIL_WIDTH * 2 + 1];
   for (int b = -VLASOV_STENCIL_WIDTH; b <= VLASOV_STENCIL_WIDTH; ++b) {
      SpatialCell* srcCell = source_neighbors[b + VLASOV_STENCIL_WIDTH];
      const vmesh::LocalID blockLID = srcCell->get_velocity_block_local_id(blockGID,popID);
      if (blockLID != srcCell->invalid_local_id()) {
         blockDatas[b + VLASOV_STENCIL_WIDTH] = srcCell->get_data(blockLID,popID);
         //prefetch storage pointers to L1
         _mm_prefetch((char *)(blockDatas[b + VLASOV_STENCIL_WIDTH]), _MM_HINT_T0);
         _mm_prefetch((char *)(blockDatas[b + VLASOV_STENCIL_WIDTH]) + 64, _MM_HINT_T0);
         _mm_prefetch((char *)(blockDatas[b + VLASOV_STENCIL_WIDTH]) + 128, _MM_HINT_T0);
         _mm_prefetch((char *)(blockDatas[b + VLASOV_STENCIL_WIDTH]) + 192, _MM_HINT_T0);
         if(VPREC  == 8) {
            //prefetch storage pointers to L1
            _mm_prefetch((char *)(blockDatas[b + VLASOV_STENCIL_WIDTH]) + 256, _MM_HINT_T0);
            _mm_prefetch((char *)(blockDatas[b + VLASOV_STENCIL_WIDTH]) + 320, _MM_HINT_T0);
            _mm_prefetch((char *)(blockDatas[b + VLASOV_STENCIL_WIDTH]) + 384, _MM_HINT_T0);
            _mm_prefetch((char *)(blockDatas[b + VLASOV_STENCIL_WIDTH]) + 448, _MM_HINT_T0);
         }
      }
      else{
         blockDatas[b + VLASOV_STENCIL_WIDTH] = NULL;
      }
   }
 
   //  Copy volume averages of this block from all spatial cells:
   for (int b = -VLASOV_STENCIL_WIDTH; b <= VLASOV_STENCIL_WIDTH; ++b) {
      if(blockDatas[b + VLASOV_STENCIL_WIDTH] != NULL) {
         Realv blockValues[WID3];
         const Realf* block_data = blockDatas[b + VLASOV_STENCIL_WIDTH];
         // Copy data to a temporary array and transpose values so that mapping is along k direction.
         // spatial source_neighbors already taken care of when
         // creating source_neighbors table. If a normal spatial cell does not
         // simply have the block, its value will be its null_block which
         // is fine. This null_block has a value of zero in data, and that
         // is thus the velocity space boundary
         for (uint i=0; i<WID3; ++i) {
            blockValues[i] = block_data[cellid_transpose[i]];
         }
      
         // now load values into the actual values table..
         uint offset =0;
         for (uint k=0; k<WID; ++k) {
            for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
               // store data, when reading data from data we swap dimensions 
               // using precomputed plane_index_to_id and cell_indices_to_id
               values[i_trans_ps_blockv(planeVector, k, b)].load(blockValues + offset);
               offset += VECL;
            }
         }
      } else {
         uint cellid=0;
         for (uint k=0; k<WID; ++k) {
            for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++) {
               values[i_trans_ps_blockv(planeVector, k, b)] = Vec(0);
            }
         }
      }
   }
}

/* 
   Here we map from the current time step grid, to a target grid which
   is the lagrangian departure grid (so th grid at timestep +dt,
   tracked backwards by -dt). This is done in ordinary space in the translation step

   This function can, and should be, safely called in a parallel
   OpenMP region (as long as it does only one dimension per parallel
   refion). It is safe as each thread only computes certain blocks (blockID%tnum_threads = thread_num */

bool trans_map_1d(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const vector<CellID>& localPropagatedCells,
                  const vector<CellID>& remoteTargetCells,
                  const uint dimension,
                  const Realv dt,
                  const uint popID) {
   // values used with an stencil in 1 dimension, initialized to 0. 
   // Contains a block, and its spatial neighbours in one dimension.
   Realv dz,z_min, dvz,vz_min;
   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
   unsigned char  cellid_transpose[WID3]; /*< defines the transpose for the solver internal (transposed) id: i + j*WID + k*WID2 to actual one*/

   if(localPropagatedCells.size() == 0) 
      return true; 
//vector with all cells
   vector<CellID> allCells(localPropagatedCells);
   allCells.insert(allCells.end(), remoteTargetCells.begin(), remoteTargetCells.end());
   
   const uint nSourceNeighborsPerCell = 1 + 2 * VLASOV_STENCIL_WIDTH;
   std::vector<SpatialCell*> allCellsPointer(allCells.size());
   std::vector<SpatialCell*> sourceNeighbors(localPropagatedCells.size() * nSourceNeighborsPerCell);
   std::vector<SpatialCell*> targetNeighbors(3 * localPropagatedCells.size() );

   
#pragma omp parallel for
   for(uint celli = 0; celli < allCells.size(); celli++){         
      allCellsPointer[celli] = mpiGrid[allCells[celli]];
   }
   
   
#pragma omp parallel for
   for(uint celli = 0; celli < localPropagatedCells.size(); celli++){         
         // compute spatial neighbors, separately for targets and source. In
         // source cells we have a wider stencil and take into account
         // boundaries. For targets we only have actual cells as we do not
         // want to propagate boundary cells (array may contain
         // INVALID_CELLIDs at boundaries).
      compute_spatial_source_neighbors(mpiGrid, localPropagatedCells[celli], dimension, sourceNeighbors.data() + celli * nSourceNeighborsPerCell);
      compute_spatial_target_neighbors(mpiGrid, localPropagatedCells[celli], dimension, targetNeighbors.data() + celli * 3);
   }
   
    
   //Get a unique sorted list of blockids that are in any of the
   // propagated cells. First use set for this, then add to vector (may not
   // be the most nice way to do this and in any case we could do it along
   // dimension for data locality reasons => copy acc map column code, TODO: FIXME
   std::unordered_set<vmesh::GlobalID> unionOfBlocksSet;
   


   for(uint celli = 0; celli < allCellsPointer.size(); celli++) {
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = allCellsPointer[celli]->get_velocity_mesh(popID);
      for (vmesh::LocalID block_i=0; block_i< vmesh.size(); ++block_i) {
         unionOfBlocksSet.insert(vmesh.getGlobalID(block_i));
      }
   }
   
   std::vector<vmesh::GlobalID> unionOfBlocks;
   unionOfBlocks.reserve(unionOfBlocksSet.size());
   for(const auto blockGID:  unionOfBlocksSet) {
      unionOfBlocks.push_back(blockGID);
   }
   
    

   
   const uint8_t REFLEVEL=0;
   const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = allCellsPointer[0]->get_velocity_mesh(popID);
   // set cell size in dimension direction
   dvz = vmesh.getCellSize(REFLEVEL)[dimension];
   vz_min = vmesh.getMeshMinLimits()[dimension];
   switch (dimension) {
   case 0:
      dz = P::dx_ini;
      z_min = P::xmin;      
      // set values in array that is used to convert block indices 
      // to global ID using a dot product.
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
   case 1:
      dz = P::dy_ini;
      z_min = P::ymin;
      // set values in array that is used to convert block indices 
      // to global ID using a dot product
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
   case 2:
      dz = P::dz_ini;
      z_min = P::zmin;
      // set values in array that is used to convert block indices
      // to global id using a dot product.
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   default:
      cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
      abort();
      break;
   }
         
   // init plane_index_to_id
   for (uint k=0; k<WID; ++k) {
      for (uint j=0; j<WID; ++j) {
         for (uint i=0; i<WID; ++i) {
            const uint cell =
               i * cell_indices_to_id[0] +
               j * cell_indices_to_id[1] +
               k * cell_indices_to_id[2];
            cellid_transpose[ i + j * WID + k * WID2] = cell;
         }
      }
   }

   const Realv i_dz=1.0/dz;
   
   int t1 = phiprof::initializeTimer("mapping");
   int t2 = phiprof::initializeTimer("store");
   
   
#pragma omp parallel 
   {      
      std::vector<Realf> targetBlockData(3 * localPropagatedCells.size() * WID3);
      std::vector<bool> targetsValid(localPropagatedCells.size());
      std::vector<vmesh::LocalID> allCellsBlockLocalID(allCells.size());

      
      
#pragma omp for schedule(guided)
      for(uint blocki = 0; blocki < unionOfBlocks.size(); blocki++){
         vmesh::GlobalID blockGID = unionOfBlocks[blocki];
         phiprof::start(t1);
         
         for(uint celli = 0; celli < allCellsPointer.size(); celli++){
            allCellsBlockLocalID[celli] = allCellsPointer[celli]->get_velocity_block_local_id(blockGID, popID);
         }

      
         for(uint celli = 0; celli < localPropagatedCells.size(); celli++){
            SpatialCell *spatial_cell = allCellsPointer[celli];
            const CellID cellID =  localPropagatedCells[celli];
            const vmesh::LocalID blockLID = allCellsBlockLocalID[celli];
            
            //Reset list of valid targets, will be set to true later for those
            //that are valid
            targetsValid[celli] = false;
            
            if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID() ||
                get_spatial_neighbor(mpiGrid, cellID, true, 0, 0, 0) == INVALID_CELLID) {
               //do nothing if it is not a normal cell, or a cell that is in the
               //first boundary layer, or the block does not exist in this
               //spatial cell
               continue;
            }

          
            // Vector buffer where we write data, initialized to 0*/
            Vec targetVecValues[3 * WID3 / VECL];
            // init target_values
            for (uint i = 0; i< 3 * WID3 / VECL; ++i) {
               targetVecValues[i] = Vec(0.0);
            }
          
            // buffer where we read in source data. i index vectorized
            Vec values[(1 + 2 * VLASOV_STENCIL_WIDTH) * WID3 / VECL];
            copy_trans_block_data(sourceNeighbors.data() + celli * nSourceNeighborsPerCell, blockGID, values, cellid_transpose, popID);
            velocity_block_indices_t block_indices;
            uint8_t refLevel;
            vmesh.getIndices(blockGID,refLevel, block_indices[0], block_indices[1], block_indices[2]);
          
            //i,j,k are now relative to the order in which we copied data to the values array. 
            //After this point in the k,j,i loops there should be no branches based on dimensions
            //
            //Note that the i dimension is vectorized, and thus there are no loops over i
            for (uint k=0; k<WID; ++k) {
               const Realv cell_vz = (block_indices[dimension] * WID + k + 0.5) * dvz + vz_min; //cell centered velocity
               const Realv z_translation = cell_vz * dt * i_dz; // how much it moved in time dt (reduced units)
               const int target_scell_index = (z_translation > 0) ? 1: -1; //part of density goes here (cell index change along spatial direcion)
             
               //the coordinates (scaled units from 0 to 1) between which we will
               //integrate to put mass in the target  neighboring cell. 
               //As we are below CFL<1, we know
               //that mass will go to two cells: current and the new one.
               Realv z_1,z_2;
               if ( z_translation < 0 ) {
                  z_1 = 0;
                  z_2 = -z_translation; 
               } else {
                  z_1 = 1.0 - z_translation;
                  z_2 = 1.0;
               }
               
               for (uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++) {         
                  //compute reconstruction
#ifdef TRANS_SEMILAG_PLM
                  Vec a[3];
                  compute_plm_coeff(values + i_trans_ps_blockv(planeVector, k, -VLASOV_STENCIL_WIDTH), VLASOV_STENCIL_WIDTH, a);
#endif
#ifdef TRANS_SEMILAG_PPM
                  Vec a[3];
                  //Check that stencil width VLASOV_STENCIL_WIDTH in grid.h corresponds to order of face estimates  (h4 & h5 =2, H6=3, h8=4)
                  compute_ppm_coeff(values + i_trans_ps_blockv(planeVector, k, -VLASOV_STENCIL_WIDTH), h4, VLASOV_STENCIL_WIDTH, a);
#endif
#ifdef TRANS_SEMILAG_PQM
                  Vec a[5];
                  //Check that stencil width VLASOV_STENCIL_WIDTH in grid.h corresponds to order of face estimates (h4 & h5 =2, H6=3, h8=4)
                  compute_pqm_coeff(values + i_trans_ps_blockv(planeVector, k, -VLASOV_STENCIL_WIDTH), h6, VLASOV_STENCIL_WIDTH, a);
#endif
          
#ifdef TRANS_SEMILAG_PLM
                  const Vec ngbr_target_density =
                     z_2 * ( a[0] + z_2 * a[1] ) -
                     z_1 * ( a[0] + z_1 * a[1] );
#endif
#ifdef TRANS_SEMILAG_PPM
                  const Vec ngbr_target_density =
                     z_2 * ( a[0] + z_2 * ( a[1] + z_2 * a[2] ) ) -
                     z_1 * ( a[0] + z_1 * ( a[1] + z_1 * a[2] ) );
#endif
#ifdef TRANS_SEMILAG_PQM
                  const Vec ngbr_target_density =
                     z_2 * ( a[0] + z_2 * ( a[1] + z_2 * ( a[2] + z_2 * ( a[3] + z_2 * a[4] ) ) ) ) -
                     z_1 * ( a[0] + z_1 * ( a[1] + z_1 * ( a[2] + z_1 * ( a[3] + z_1 * a[4] ) ) ) );
#endif
                  targetVecValues[i_trans_pt_blockv(planeVector, k, target_scell_index)] +=  ngbr_target_density;                     //in the current original cells we will put this density        
                  targetVecValues[i_trans_pt_blockv(planeVector, k, 0)] +=  values[i_trans_ps_blockv(planeVector, k, 0)] - ngbr_target_density; //in the current original cells we will put the rest of the original density
               }
            }
         
            //Store final vector data in temporary data for all target blocks,
            //and mark that this celli produced valid targets
         
            targetsValid[celli] = true;
            for (int b = -1; b< 2 ; ++b) {
               Realv vector[VECL];
               for (uint k=0; k<WID; ++k) {
                  for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
                     targetVecValues[i_trans_pt_blockv(planeVector, k, b)].store(vector);
#pragma ivdep
#pragma GCC ivdep
                     for(uint i = 0; i< VECL; i++){
                        // store data, when reading data from data we swap
                        // dimensions 
                        // using precomputed plane_index_to_id and
                        // cell_indices_to_id
                        targetBlockData[(celli * 3 + b + 1) * WID3 +  cellid_transpose[i + planeVector * VECL + k * WID2]] = 
                           vector[i];
                     }
                  }
               }
            }
         }
      
         phiprof::stop(t1);
         phiprof::start(t2);
               
         //reset blocks in all non-sysboundary spatial cells for this block id
         for(uint celli = 0; celli < allCellsPointer.size(); celli++){
            SpatialCell* spatial_cell = allCellsPointer[celli];
            if(spatial_cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
               const vmesh::LocalID blockLID = allCellsBlockLocalID[celli];
               if (blockLID != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
                  Realf* blockData = spatial_cell->get_data(blockLID, popID);
                  for(int i = 0; i < WID3; i++) {
                     blockData[i] = 0.0;
                  }
               }
            }
         }
      
         //store values from target_values array to the actual blocks
         for(uint celli = 0; celli < localPropagatedCells.size(); celli++){
            if(targetsValid[celli]) {
               for(uint ti = 0; ti < 3; ti++) {
                  SpatialCell* spatial_cell = targetNeighbors[celli * 3 + ti];
                  if(spatial_cell ==NULL) {
                     //invalid target spatial cell
                     continue;
                  }
               
                  const vmesh::LocalID blockLID = spatial_cell->get_velocity_block_local_id(blockGID, popID);
                  if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
                     // block does not exist. If so, we do not create it and add stuff to it here.
                     // We have already created blocks around blocks with content in
                     // spatial sense, so we have no need to create even more blocks here
                     // TODO add loss counter
                     continue;
                  }
                  Realf* blockData = spatial_cell->get_data(blockLID, popID);
                  for(int i = 0; i < WID3 ; i++) {
                     blockData[i] += targetBlockData[(celli * 3 + ti) * WID3 + i];
                  }
               }
            }
         
         }
         phiprof::stop(t2);

      
      } //loop over set of blocks on process
   }
   

   return true;
}

/*!

  This function communicates the mapping on process boundaries, and then updates the data to their correct values.
  TODO, this could be inside an openmp region, in which case some m ore barriers and masters should be added

  \par dimension: 0,1,2 for x,y,z
  \par direction: 1 for + dir, -1 for - dir
*/
void update_remote_mapping_contribution(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const uint dimension,
   int direction,
   const uint popID) {
   
   const vector<CellID> local_cells = mpiGrid.get_cells();
   const vector<CellID> remote_cells = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   vector<CellID> receive_cells;
   vector<CellID> send_cells;
   vector<Realf*> receiveBuffers;
   
   //normalize
   if(direction > 0) direction = 1;
   if(direction < 0) direction = -1;
   for (size_t c=0; c<remote_cells.size(); ++c) {
      SpatialCell *ccell = mpiGrid[remote_cells[c]];
      //default values, to avoid any extra sends and receives
      ccell->neighbor_block_data = ccell->get_data(popID);
      ccell->neighbor_number_of_blocks = 0;
   }

   //TODO: prepare arrays, make parallel by avoidin push_back and by checking also for other stuff
   for (size_t c=0; c<local_cells.size(); ++c) {
      SpatialCell *ccell = mpiGrid[local_cells[c]];
      //default values, to avoid any extra sends and receives
      ccell->neighbor_block_data = ccell->get_data(popID);
      ccell->neighbor_number_of_blocks = 0;
      CellID p_ngbr,m_ngbr;
      switch (dimension) {
      case 0:
         p_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], false, direction, 0, 0); //p_ngbr is target, if in boundaries then it is not updated
         m_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], true, -direction, 0, 0); //m_ngbr is source, first boundary layer is propagated so that it flows into system
         break;
      case 1:
         p_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], false, 0, direction, 0); //p_ngbr is target, if in boundaries then it is not update
         m_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], true, 0, -direction, 0); //m_ngbr is source, first boundary layer is propagated so that it flows into system
         break;
      case 2:
         p_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], false, 0, 0, direction); //p_ngbr is target, if in boundaries then it is not update
         m_ngbr=get_spatial_neighbor(mpiGrid, local_cells[c], true, 0, 0, -direction); //m_ngbr is source, first boundary layer is propagated so that it flows into system
         break;
      default:
         cerr << "Dimension wrong at (impossible!) "<< __FILE__ <<":" << __LINE__<<endl;
         exit(1);
         break;
      }
      //internal cell, not much to do
      if (mpiGrid.is_local(p_ngbr) && mpiGrid.is_local(m_ngbr)) continue;

      SpatialCell *pcell = NULL;
      if (p_ngbr != INVALID_CELLID) pcell = mpiGrid[p_ngbr];
      SpatialCell *mcell = NULL;
      if (m_ngbr != INVALID_CELLID) mcell = mpiGrid[m_ngbr];
      if (p_ngbr != INVALID_CELLID && pcell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) 
         if (!mpiGrid.is_local(p_ngbr) && do_translate_cell(ccell)) {
            //if (p_ngbr != INVALID_CELLID && !mpiGrid.is_local(p_ngbr) && do_translate_cell(ccell)) {
            //Send data in p_ngbr target array that we just
            //mapped to if 1) it is a valid target,
            //2) is remote cell, 3) if the source cell in center was
            //translated
            ccell->neighbor_block_data = pcell->get_data(popID);
            ccell->neighbor_number_of_blocks = pcell->get_number_of_velocity_blocks(popID);
            send_cells.push_back(p_ngbr);
         }
      if (m_ngbr != INVALID_CELLID &&
          !mpiGrid.is_local(m_ngbr) &&
          ccell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         //Receive data that mcell mapped to ccell to this local cell
         //data array, if 1) m is a valid source cell, 2) center cell is to be updated (normal cell) 3) m is remote
         //we will here allocate a receive buffer, since we need to aggregate values
         mcell->neighbor_number_of_blocks = ccell->get_number_of_velocity_blocks(popID);
         mcell->neighbor_block_data = (Realf*) aligned_malloc(mcell->neighbor_number_of_blocks * WID3 * sizeof(Realf), 64);
         
         receive_cells.push_back(local_cells[c]);
         receiveBuffers.push_back(mcell->neighbor_block_data);
      }
   }
    
   // Do communication
   SpatialCell::setCommunicatedSpecies(popID);
   SpatialCell::set_mpi_transfer_type(Transfer::NEIGHBOR_VEL_BLOCK_DATA);
   switch(dimension) {
   case 0:
      if(direction > 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_P_X_NEIGHBORHOOD_ID);
      if(direction < 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_M_X_NEIGHBORHOOD_ID);
      break;
   case 1:
      if(direction > 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_P_Y_NEIGHBORHOOD_ID);
      if(direction < 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_M_Y_NEIGHBORHOOD_ID);
      break;
   case 2:
      if(direction > 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_P_Z_NEIGHBORHOOD_ID);
      if(direction < 0) mpiGrid.update_copies_of_remote_neighbors(SHIFT_M_Z_NEIGHBORHOOD_ID);
      break;
   }
   
#pragma omp parallel
   {
      //reduce data: sum received data in the data array to 
      // the target grid in the temporary block container
      for (size_t c=0; c < receive_cells.size(); ++c) {
         SpatialCell* spatial_cell = mpiGrid[receive_cells[c]];
         Realf *blockData = spatial_cell->get_data(popID);
          
#pragma omp for 
         for(unsigned int cell = 0; cell<VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(popID); ++cell) {
            blockData[cell] += receiveBuffers[c][cell];
         }
      }
       
      // send cell data is set to zero. This is to avoid double copy if
      // one cell is the neighbor on bot + and - side to the same
      // process
      for (size_t c=0; c<send_cells.size(); ++c) {
         SpatialCell* spatial_cell = mpiGrid[send_cells[c]];
         Realf * blockData = spatial_cell->get_data(popID);
           
#pragma omp for nowait
         for(unsigned int cell = 0; cell< VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(popID); ++cell) {
            // copy received target data to temporary array where target data is stored.
            blockData[cell] = 0;
         }
      }
   }
    
   //and finally free temporary receive buffer
   for (size_t c=0; c < receiveBuffers.size(); ++c) {
      aligned_free(receiveBuffers[c]);
   }
}

CellID selectNeighbor(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> &grid,
                      CellID id, int dimension = 0, uint path = 0) {

   const auto neighbors = grid.get_face_neighbors_of(id);
   const int myProcess = grid.get_process(id);
   
   vector < CellID > myNeighbors;
   // Collect neighbor ids in the positive direction of the chosen dimension,
   // that are on the same process as the origin.
   // Note that dimension indexing starts from 1 (of course it does)
   for (const auto cell : neighbors) {
      if (cell.second == dimension + 1 && grid.get_process(cell.first) == myProcess)
         myNeighbors.push_back(cell.first);
   }

   CellID neighbor;
  
   switch( myNeighbors.size() ) {
      // Since refinement can only increase by 1 level the only possibilities
      // Should be 0 neighbors, 1 neighbor or 4 neighbors.
   case 0 : {
      // did not find neighbors
      neighbor = INVALID_CELLID;
      break;
   }
   case 1 : {
      neighbor = myNeighbors[0];
      break;
   }
   case 4 : {
      neighbor = myNeighbors[path];
      break;
   }
   default: {
      // something is wrong
      neighbor = INVALID_CELLID;
      throw "Invalid neighbor count!";
      break;
   }
   }

   return neighbor;
  
}

setOfPencils buildPencilsWithNeighbors( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> &grid, 
					setOfPencils &pencils, CellID startingId,
					vector<CellID> ids, uint dimension, 
					vector<uint> path) {

   const bool debug = false;
   CellID nextNeighbor;
   uint id = startingId;
   uint startingRefLvl = grid.get_refinement_level(id);

   if( ids.size() == 0 )
      ids.push_back(startingId);

   // If the cell where we start is refined, we need to figure out which path
   // to follow in future refined cells. This is a bit hacky but we have to
   // use the order or the children of the parent cell to figure out which
   // corner we are in.
   // Maybe you could use physical coordinates here?
   if( startingRefLvl > path.size() ) {
      for ( uint i = path.size(); i < startingRefLvl; i++) {
         auto parent = grid.get_parent(id);
         auto children = grid.get_all_children(parent);
         auto it = std::find(children.begin(),children.end(),id);
         auto index = std::distance(children.begin(),it);
         auto index2 = index;
      
         switch( dimension ) {
         case 0: {
            index2 = index / 2;
            break;
         }
         case 1: {
            index2 = index - index / 2;
            break;
         }
         case 2: {
            index2 = index % 4;
            break;
         }
         }
         path.insert(path.begin(),index2);
         id = parent;
      }
   }

   id = startingId;

   bool periodic;
   
   while (id > 0) {

      periodic = false;

      // Find the refinement level in the neighboring cell. Any neighbor will do
      // since refinement level can only increase by 1 between neighbors.
      nextNeighbor = selectNeighbor(grid,id,dimension);

      // If there are no neighbors, we can stop.
      if (nextNeighbor == 0)
         break;
    
      uint refLvl = grid.get_refinement_level(nextNeighbor);

      if (refLvl > 0) {
    
         // If we have encountered this refinement level before and stored
         // the path this builder follows, we will just take the same path
         // again.
         if ( path.size() >= refLvl ) {
      
            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << " ";
               std::cout << std::endl;
            }
	
            nextNeighbor = selectNeighbor(grid,id,dimension,path[refLvl-1]);      
	
         } else {
	
            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have NOT seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << ' ';
               std::cout << std::endl;
            }
	
            // New refinement level, create a path through each neighbor cell
            for ( uint i : {0,1,2,3} ) {
	  
               vector < uint > myPath = path;
               myPath.push_back(i);
	  
               nextNeighbor = selectNeighbor(grid,id,dimension,myPath.back());
	  
               if ( i == 3 ) {
	    
                  // This builder continues with neighbor 3
                  ids.push_back(nextNeighbor);
                  path = myPath;
	    
               } else {
	    
                  // Spawn new builders for neighbors 0,1,2
                  buildPencilsWithNeighbors(grid,pencils,id,ids,dimension,myPath);
	    
               }
	  
            }
	
         }

      } else {
         if(debug) {
            std::cout << "I am cell " << id << ". ";
            std::cout << " I am on refinement level 0." << std::endl;
         }
      }// Closes if (refLvl == 0)

      // If we found a neighbor, add it to the list of ids for this pencil.
      if(nextNeighbor != INVALID_CELLID) {
         if (debug) {
            std::cout << " Next neighbor is " << nextNeighbor << "." << std::endl;
         }         

         for (auto id : ids) {
            if (nextNeighbor == id) {
               periodic = true;
            }
         }
         if (periodic) {
            // Exit the while loop
            id = -1;
         } else {
            ids.push_back(nextNeighbor);
            // Move to the next cell.
            id = nextNeighbor;
         }
      }
    
   } // Closes while loop

   // Get the x,y - coordinates of the pencil (in the direction perpendicular to the pencil)
   const auto coordinates = grid.get_center(ids[0]);
   double x,y;
   uint ix,iy,iz;
      
   switch(dimension) {
   case 0: 
      ix = 1;
      iy = 2;
      iz = 0;
      break;
   
   case 1: 
      ix = 2;
      iy = 0;
      iz = 1;
      break;
   
   case 2: 
      ix = 0;
      iy = 1;
      iz = 2;
      break;
   
   default: 
      ix = 0;
      iy = 1;
      iz = 2;
      break;   
   }
   
   x = coordinates[ix];
   y = coordinates[iy];

   pencils.addPencil(ids,x,y,periodic);
   
   return pencils;
  
}

//void propagatePencil(Vec dr[], Vec values, Vec z_translation, uint blocks_per_dim ) {
void propagatePencil(Vec* dz, Vec* values, const uint dimension, const uint blockGID, const Realv dt,
                     const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> &vmesh, const uint lengthOfPencil, bool debugflag) {

   // Get velocity data from vmesh that we need later to calculate the translation
   velocity_block_indices_t block_indices;
   uint8_t refLevel;
   vmesh.getIndices(blockGID,refLevel, block_indices[0], block_indices[1], block_indices[2]);
   Realv dvz = vmesh.getCellSize(refLevel)[dimension];
   Realv vz_min = vmesh.getMeshMinLimits()[dimension];
   
   // Assuming 1 neighbor in the target array because of the CFL condition
   // In fact propagating to > 1 neighbor will give an error
   const uint nTargetNeighborsPerPencil = 1;

   // Vector buffer where we write data, initialized to 0*/
   Vec targetValues[(lengthOfPencil + 2 * nTargetNeighborsPerPencil) * WID3 / VECL];
   
   for (uint i = 0; i < (lengthOfPencil + 2 * nTargetNeighborsPerPencil) * WID3 / VECL; i++) {
      
      // init target_values
      targetValues[i] = Vec(0.0);
      
   }
   // Go from 0 to length here to propagate all the cells in the pencil
   for (uint i = 0; i < lengthOfPencil; i++){

      // The source array is padded by VLASOV_STENCIL_WIDTH on both sides.
      uint i_source   = i + VLASOV_STENCIL_WIDTH;
      
      for (uint k = 0; k < WID; ++k) {

         const Realv cell_vz = (block_indices[dimension] * WID + k + 0.5) * dvz + vz_min; //cell centered velocity
         const Vec z_translation = cell_vz * dt / dz[i_source]; // how much it moved in time dt (reduced units)

         // Determine direction of translation
         // part of density goes here (cell index change along spatial direcion)
         Vecb positiveTranslationDirection = (z_translation > Vec(0.0));
         
         // Calculate normalized coordinates in current cell.
         // The coordinates (scaled units from 0 to 1) between which we will
         // integrate to put mass in the target  neighboring cell.
         // Normalize the coordinates to the origin cell. Then we scale with the difference
         // in volume between target and origin later when adding the integrated value.
         Vec z_1,z_2;
         z_1 = select(positiveTranslationDirection, 1.0 - z_translation, 0.0);
         z_2 = select(positiveTranslationDirection, 1.0, - z_translation);

         if( horizontal_or(abs(z_1) > Vec(1.0)) || horizontal_or(abs(z_2) > Vec(1.0)) ) {
            std::cout << "Error, CFL condition violated\n";
            std::cout << "Exiting\n";
            std::exit(1);
         }
         
         for (uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++) {   
      
            // Compute polynomial coefficients
            Vec a[3];
            // Dz: is a padded array, pointer can point to the beginning, i + VLASOV_STENCIL_WIDTH will get the right cell.
            // values: transpose function adds VLASOV_STENCIL_WIDTH to the block index, therefore we substract it here, then
            // i + VLASOV_STENCIL_WIDTH will point to the right cell. Complicated! Why! Sad! MVGA!
            compute_ppm_coeff_nonuniform(dz,
                                         values + i_trans_ps_blockv_pencil(planeVector, k, i-VLASOV_STENCIL_WIDTH, lengthOfPencil),
                                         h4, VLASOV_STENCIL_WIDTH, a);
            
            // Compute integral
            const Vec ngbr_target_density =
               z_2 * ( a[0] + z_2 * ( a[1] + z_2 * a[2] ) ) -
               z_1 * ( a[0] + z_1 * ( a[1] + z_1 * a[2] ) );
            // Store mapped density in two target cells
            // in the neighbor cell we will put this density
            targetValues[i_trans_pt_blockv(planeVector, k, i + 1)] += select( positiveTranslationDirection,
                                                  ngbr_target_density * dz[i_source] / dz[i_source + 1],Vec(0.0));
            targetValues[i_trans_pt_blockv(planeVector, k, i - 1 )] += select(!positiveTranslationDirection,
                                                 ngbr_target_density * dz[i_source] / dz[i_source - 1],Vec(0.0));
            
            // in the current original cells we will put the rest of the original density
            targetValues[i_trans_pt_blockv(planeVector, k, i)] +=
               values[i_trans_ps_blockv_pencil(planeVector, k, i, lengthOfPencil)] - ngbr_target_density;
         }
      }
   }

   // Write target data into source data
   // VLASOV_STENCIL_WIDTH >= nTargetNeighborsPerPencil is required (default 2 >= 1)

   for (int i = 0; i < lengthOfPencil + 2 * nTargetNeighborsPerPencil; i++) {

      for (uint k = 0; k < WID; ++k) {
         
         for (uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++) {

            values[i_trans_ps_blockv_pencil(planeVector, k, i - 1, lengthOfPencil)] =
               targetValues[i_trans_pt_blockv(planeVector, k, i - 1)];
            
         }
      }
   }  
}

void get_seed_ids(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const vector<CellID> &localPropagatedCells,
                  const uint dimension,
                  vector<CellID> &seedIds) {

   const bool debug = false;
   
   //#pragma omp parallel for      
   for(auto celli: localPropagatedCells) {
      int myProcess = mpiGrid.get_process(celli);
      // Collect a list of cell ids that do not have a neighbor in the negative direction
      // These are the seed ids for the pencils.
      vector<CellID> negativeNeighbors;
      // Returns all neighbors as (id, direction-dimension) pairs.
      for ( const auto neighbor : mpiGrid.get_face_neighbors_of(celli ) ) {

         if ( mpiGrid.get_process(neighbor.first) == myProcess ) {
            // select the neighbor in the negative dimension of the propagation
            if (neighbor.second == - (static_cast<int>(dimension) + 1)) {
               
               // add the id of the neighbor to a list if it's on the same process
               negativeNeighbors.push_back(neighbor.first);
               
            }

         }
      }
      //cout << endl;
      // if no neighbors were found in the negative direction, add this cell id to the seed cells
      if (negativeNeighbors.size() == 0)
         seedIds.push_back(celli);    
   }

   // If no seed ids were found, let's assume we have a periodic boundary and
   // a single process in the dimension of propagation. In this case we start from
   // the first cells of the plane perpendicular to the propagation dimension
   if (seedIds.size() == 0) {
      for (uint ix = 0; ix < P::xcells_ini; ix++) {
         for (uint iy = 0; iy < P::ycells_ini; iy++) {
            for (uint iz = 0; iz < P::zcells_ini; iz++) {
               CellID seedId;
               switch (dimension) {
               case 0:
                  // yz - plane
                  if(ix == 0) {
                     seedId = P::xcells_ini * P::ycells_ini * iz + P::xcells_ini * iy + 1;
                     if(mpiGrid.get_process(seedId) == mpiGrid.get_process(localPropagatedCells[0]))
                        seedIds.push_back(seedId);

                  }
                  break;
               case 1:                  
                  // xz - plane
                  if(iy == 0) {
                     seedId = P::xcells_ini * P::ycells_ini * iz + ix + 1;
                     if(mpiGrid.get_process(seedId) == mpiGrid.get_process(localPropagatedCells[0]))
                        seedIds.push_back(seedId);

                  }
                  break;
               case 2:
                  // xy - plane
                  if(iz == 0) {
                     seedId = P::xcells_ini * iy + ix + 1;
                     if(mpiGrid.get_process(seedId) == mpiGrid.get_process(localPropagatedCells[0]))
                        seedIds.push_back(seedId);
                  }
                  break;
               }
            }
         }
      }
   }

   if(debug) {
      cout << "Number of seed ids is " << seedIds.size() << endl;
      cout << "Seed ids are: ";
      for (const auto seedId : seedIds) {
         cout << seedId << " ";
      }
      cout << endl;
   }   
}




/* Copy the data to the temporary values array, so that the
 * dimensions are correctly swapped. Also, copy the same block for
 * then neighboring spatial cells (in the dimension). neighbors
 * generated with compute_spatial_neighbors_wboundcond).
 * 
 * This function must be thread-safe.
 *
 * @param source_neighbors Array containing the VLASOV_STENCIL_WIDTH closest 
 * spatial neighbors of this cell in the propagated dimension.
 * @param blockGID Global ID of the velocity block.
 * @param int lengthOfPencil Number of spatial cells in pencil
 * @param values Vector where loaded data is stored.
 * @param cellid_transpose
 * @param popID ID of the particle species.
 */
void copy_trans_block_data_amr(
    SpatialCell** source_neighbors,
    const vmesh::GlobalID blockGID,
    int lengthOfPencil,
    Vec* values,
    const unsigned char* const cellid_transpose,
    const uint popID) { 

   // Allocate data pointer for all blocks in pencil. Pad on both ends by VLASOV_STENCIL_WIDTH
   Realf* blockDataPointer[lengthOfPencil + 2 * VLASOV_STENCIL_WIDTH];   

   for (int b = -VLASOV_STENCIL_WIDTH; b < lengthOfPencil + VLASOV_STENCIL_WIDTH; b++) {
      // Get cell pointer and local block id
      SpatialCell* srcCell = source_neighbors[b + VLASOV_STENCIL_WIDTH];
         
      const vmesh::LocalID blockLID = srcCell->get_velocity_block_local_id(blockGID,popID);
      if (blockLID != srcCell->invalid_local_id()) {
         // Get data pointer
         blockDataPointer[b + VLASOV_STENCIL_WIDTH] = srcCell->get_data(blockLID,popID);
         // //prefetch storage pointers to L1
         // _mm_prefetch((char *)(blockDataPointer[b + VLASOV_STENCIL_WIDTH]), _MM_HINT_T0);
         // _mm_prefetch((char *)(blockDataPointer[b + VLASOV_STENCIL_WIDTH]) + 64, _MM_HINT_T0);
         // _mm_prefetch((char *)(blockDataPointer[b + VLASOV_STENCIL_WIDTH]) + 128, _MM_HINT_T0);
         // _mm_prefetch((char *)(blockDataPointer[b + VLASOV_STENCIL_WIDTH]) + 192, _MM_HINT_T0);
         // if(VPREC  == 8) {
         //   //prefetch storage pointers to L1
         //   _mm_prefetch((char *)(blockDataPointer[b + VLASOV_STENCIL_WIDTH]) + 256, _MM_HINT_T0);
         //   _mm_prefetch((char *)(blockDataPointer[b + VLASOV_STENCIL_WIDTH]) + 320, _MM_HINT_T0);
         //   _mm_prefetch((char *)(blockDataPointer[b + VLASOV_STENCIL_WIDTH]) + 384, _MM_HINT_T0);
         //   _mm_prefetch((char *)(blockDataPointer[b + VLASOV_STENCIL_WIDTH]) + 448, _MM_HINT_T0);
         // }
         
      } else {
         blockDataPointer[b + VLASOV_STENCIL_WIDTH] = NULL;
      }
   }
   
   //  Copy volume averages of this block from all spatial cells:
   for (int b = -VLASOV_STENCIL_WIDTH; b < lengthOfPencil + VLASOV_STENCIL_WIDTH; b++) {
      if(blockDataPointer[b + VLASOV_STENCIL_WIDTH] != NULL) {
         Realv blockValues[WID3];
         const Realf* block_data = blockDataPointer[b + VLASOV_STENCIL_WIDTH];
         // Copy data to a temporary array and transpose values so that mapping is along k direction.
         // spatial source_neighbors already taken care of when
         // creating source_neighbors table. If a normal spatial cell does not
         // simply have the block, its value will be its null_block which
         // is fine. This null_block has a value of zero in data, and that
         // is thus the velocity space boundary
         for (uint i=0; i<WID3; ++i) {
            blockValues[i] = block_data[cellid_transpose[i]];
         }
         
         // now load values into the actual values table..
         uint offset =0;
         for (uint k=0; k<WID; k++) {
            for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
               // store data, when reading data from data we swap dimensions 
               // using precomputed plane_index_to_id and cell_indices_to_id
               values[i_trans_ps_blockv_pencil(planeVector, k, b, lengthOfPencil)].load(blockValues + offset);
               offset += VECL;
            }
         }
      } else {
         for (uint k=0; k<WID; ++k) {
            for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++) {
               values[i_trans_ps_blockv_pencil(planeVector, k, b, lengthOfPencil)] = Vec(0);
            }
         }
      }
   }
}


bool trans_map_1d_amr(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                      const vector<CellID>& localPropagatedCells,
                      const vector<CellID>& remoteTargetCells,
                      const uint dimension,
                      const Realv dt,
                      const uint popID) {
   
   Realv dvz,vz_min;  
   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
   unsigned char  cellid_transpose[WID3]; /*< defines the transpose for the solver internal (transposed) id: i + j*WID + k*WID2 to actual one*/
   const uint blocks_per_dim = 1;
   // return if there's no cells to propagate
   if(localPropagatedCells.size() == 0) {
      cout << "Returning because of no cells" << endl;
      return false;
   }
  
   // Vector with all cell ids
   vector<CellID> allCells(localPropagatedCells);
   allCells.insert(allCells.end(), remoteTargetCells.begin(), remoteTargetCells.end());  

   // Vectors of pointers to the cell structs
   std::vector<SpatialCell*> allCellsPointer(allCells.size());  
   
   // Initialize allCellsPointer
   //#pragma omp parallel for
   for(uint celli = 0; celli < allCells.size(); celli++){
      allCellsPointer[celli] = mpiGrid[allCells[celli]];
   }

   // Fiddle indices x,y,z in VELOCITY SPACE
   switch (dimension) {
   case 0:
      // set values in array that is used to convert block indices 
      // to global ID using a dot product.
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
   case 1:
      // set values in array that is used to convert block indices 
      // to global ID using a dot product
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
   case 2:
      // set values in array that is used to convert block indices
      // to global id using a dot product.
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   default:
      cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
      abort();
      break;
   }
           
   // init cellid_transpose
   for (uint k=0; k<WID; ++k) {
      for (uint j=0; j<WID; ++j) {
         for (uint i=0; i<WID; ++i) {
            const uint cell =
               i * cell_indices_to_id[0] +
               j * cell_indices_to_id[1] +
               k * cell_indices_to_id[2];
            cellid_transpose[ i + j * WID + k * WID2] = cell;
         }
      }
   }
   // ****************************************************************************

   // compute pencils => set of pencils (shared datastructure)

   vector<CellID> seedIds;
   get_seed_ids(mpiGrid, localPropagatedCells, dimension, seedIds);
   
   // Empty vectors for internal use of buildPencilsWithNeighbors. Could be default values but
   // default vectors are complicated. Should overload buildPencilsWithNeighbors like suggested here
   // https://stackoverflow.com/questions/3147274/c-default-argument-for-vectorint
   vector<CellID> ids;
   vector<uint> path;

   // Output vectors for ready pencils
   setOfPencils pencils;
   vector<setOfPencils> pencilSets;

   for (const auto seedId : seedIds) {
      // Construct pencils from the seedIds into a set of pencils.
      pencils = buildPencilsWithNeighbors(mpiGrid, pencils, seedId, ids, dimension, path);
   }

   // Print out ids of pencils (if needed for debugging)
   if (false) {
      uint ibeg = 0;
      uint iend = 0;
      std::cout << "I have created " << pencils.N << " pencils along dimension " << dimension << ":\n";
      std::cout << "(x, y): indices " << std::endl;
      std::cout << "-----------------------------------------------------------------" << std::endl;
      for (uint i = 0; i < pencils.N; i++) {
         iend += pencils.lengthOfPencils[i];
         std::cout << "(" << pencils.x[i] << ", " << pencils.y[i] << "): ";
         for (auto j = pencils.ids.begin() + ibeg; j != pencils.ids.begin() + iend; ++j) {
            std::cout << *j << " ";
         }
         ibeg  = iend;
         std::cout << std::endl;
      }
   }

   
   // Add the final set of pencils to the pencilSets - vector.
   // Only one set is created for now but we retain support for multiple sets
   pencilSets.push_back(pencils);
   // ****************************************************************************
     
   const uint8_t VMESH_REFLEVEL = 0;
   
   // Get a pointer to the velocity mesh of the first spatial cell
   const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = allCellsPointer[0]->get_velocity_mesh(popID);   
   
   // set cell size in dimension direction
   dvz = vmesh.getCellSize(VMESH_REFLEVEL)[dimension];
   vz_min = vmesh.getMeshMinLimits()[dimension];
   
   // Get a unique sorted list of blockids that are in any of the
   // propagated cells. First use set for this, then add to vector (may not
   // be the most nice way to do this and in any case we could do it along
   // dimension for data locality reasons => copy acc map column code, TODO: FIXME
   // TODO: Do this separately for each pencil?
   std::unordered_set<vmesh::GlobalID> unionOfBlocksSet;    
   
   for(auto cell : allCellsPointer) {
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = cell->get_velocity_mesh(popID);
      for (vmesh::LocalID block_i=0; block_i< vmesh.size(); ++block_i) {
         unionOfBlocksSet.insert(vmesh.getGlobalID(block_i));
      }
   }
   
   std::vector<vmesh::GlobalID> unionOfBlocks;
   unionOfBlocks.reserve(unionOfBlocksSet.size());
   for(const auto blockGID:  unionOfBlocksSet) {
      unionOfBlocks.push_back(blockGID);
   }
   // ****************************************************************************
   
   int t1 = phiprof::initializeTimer("mappingAndStore");
   
   //#pragma omp parallel
   {
      //#pragma omp for schedule(guided)
      // Loop over velocity space blocks. Thread this loop (over vspace blocks) with OpenMP.    
      for(uint blocki = 0; blocki < unionOfBlocks.size(); blocki++){
         
         phiprof::start(t1);

         // Get global id of the velocity block
         vmesh::GlobalID blockGID = unionOfBlocks[blocki];

         bool debugflag = false;
         CellID debugcell;
         uint allCellsPointerIndex = 16;

         const vmesh::LocalID debugLID = allCellsPointer[allCellsPointerIndex]->get_velocity_block_local_id(blockGID, popID);
         Realf* data = allCellsPointer[allCellsPointerIndex]->get_data(debugLID,popID);
         //for (uint i = 0; i < WID3; i++) if (data[i] != 0) debugflag = true;
         velocity_block_indices_t block_indices;
         uint8_t vRefLevel;
         vmesh.getIndices(blockGID,vRefLevel, block_indices[0],
                          block_indices[1], block_indices[2]);      
      
         // Loop over sets of pencils
         // This loop only has one iteration for now
         for ( auto pencils: pencilSets ) {

            std::vector<Realf> targetBlockData((pencils.sumOfLengths + 2 * pencils.N) * WID3);
            // Allocate vectorized targetvecdata sum(lengths of pencils)*WID3 / VECL)
            // Add padding by 2 for each pencil
            Vec targetVecData[(pencils.sumOfLengths + 2 * pencils.N) * WID3 / VECL];

            // Initialize targetvecdata to 0
            for( uint i = 0; i < (pencils.sumOfLengths + 2 * pencils.N) * WID3 / VECL; i++ ) {
               targetVecData[i] = Vec(0.0);
            }

            // TODO: There's probably a smarter way to keep track of where we are writing
            //       in the target data array.
            uint targetDataIndex = 0;
            
            // Compute spatial neighbors for target cells.
            // For targets we need the local cells, plus a padding of 1 cell at both ends
            std::vector<SpatialCell*> targetCells(pencils.sumOfLengths + pencils.N * 2 );
            
            compute_spatial_target_cells_for_pencils(mpiGrid, pencils, dimension, targetCells.data());           
            
            // Loop over pencils
            uint totalTargetLength = 0;
            for(uint pencili = 0; pencili < pencils.N; pencili++){

               vector<CellID> pencilIds = pencils.getIds(pencili);
               int L = pencils.lengthOfPencils[pencili];
               uint targetLength = L + 2;
               uint sourceLength = L + 2 * VLASOV_STENCIL_WIDTH;
               
               // Compute spatial neighbors for the source cells of the pencil. In
               // source cells we have a wider stencil and take into account boundaries.
               std::vector<SpatialCell*> sourceCells(sourceLength);
               compute_spatial_source_cells_for_pencil(mpiGrid, pencils, pencili, dimension, sourceCells.data());


               Vec dz[sourceCells.size()];
               uint i = 0;
               for(auto neighbor: sourceCells) {
                  switch (dimension) {
                  case(0):
                     dz[i] = neighbor->SpatialCell::parameters[CellParams::DX];
                     break;                  
                  case(1):
                     dz[i] = neighbor->SpatialCell::parameters[CellParams::DY];
                     break;                  
                  case(2):
                     dz[i] = neighbor->SpatialCell::parameters[CellParams::DZ];
                     break;
                  }
                  i++;
               }
               
               // Allocate source data: sourcedata<length of pencil * WID3)
               // Add padding by 2 * VLASOV_STENCIL_WIDTH
               Vec sourceVecData[(sourceLength) * WID3 / VECL];                              
               
               // load data(=> sourcedata) / (proper xy reconstruction in future)
               copy_trans_block_data_amr(sourceCells.data(), blockGID, L, sourceVecData,
                                         cellid_transpose, popID);                                 
                             
               // Dz and sourceVecData are both padded by VLASOV_STENCIL_WIDTH
               // Dz has 1 value/cell, sourceVecData has WID3 values/cell
               propagatePencil(dz, sourceVecData, dimension, blockGID, dt, vmesh, L, debugflag);

               // sourcedata => targetdata[this pencil])
               for (int i = 0; i < targetLength; i++) {
                  for (uint k=0; k<WID; k++) {
                     for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
                        
                        targetVecData[i_trans_pt_blockv(planeVector, k, totalTargetLength + i - 1)] =
                           sourceVecData[i_trans_ps_blockv_pencil(planeVector, k, i - 1, L)];
                     
                     }
                  }
               }
               totalTargetLength += targetLength;
               // dealloc source data -- Should be automatic since it's declared in this iteration?
            }

            // reset blocks in all non-sysboundary neighbor spatial cells for this block id
            // At this point the data is saved in targetVecData so we can reset the spatial cells
            for (auto *spatial_cell: targetCells) {
               // Check for system boundary
               if(spatial_cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  // Get local velocity block id
                  const vmesh::LocalID blockLID = spatial_cell->get_velocity_block_local_id(blockGID, popID);
                  // Check for invalid id
                  if (blockLID != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
                     // Get a pointer to the block data
                     Realf* blockData = spatial_cell->get_data(blockLID, popID);
                     // Loop over velocity block cells
                     for(int i = 0; i < WID3; i++) {
                        blockData[i] = 0.0;
                     }
                  }
               }
            }
            
            // store_data(target_data => targetCells)  :Aggregate data for blockid to original location 
            // Loop over pencils again
            totalTargetLength = 0;
            for(uint pencili = 0; pencili < pencils.N; pencili++){
               
               int L = pencils.lengthOfPencils[pencili];
               int targetLength = L + 2;
               vector<CellID> pencilIds = pencils.getIds(pencili);

               bool debugPencilFlag = false;
               
               // Unpack the vector data

               // Loop over cells in pencil +- 1 padded cell
               for ( uint celli = 0; celli < targetLength; ++celli ) {

                  // // If the pencil is periodic, we do not write the ghost cells because
                  // // They are copies of cells that are already in the pencil
                  // // - It seems that doing this was wrong. Investigate!
                  // if(pencils.periodic[pencili] && (celli == 0 || celli == targetLength - 1))
                  //   continue;
                  
                  if(celli > 0 && pencilIds[celli - 1] == debugcell) debugPencilFlag = true;
                  
                  Realv vector[VECL];
                  // Loop over 1st vspace dimension
                  for (uint k = 0; k < WID; ++k) {
                     // Loop over 2nd vspace dimension
                     for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++) {
                        targetVecData[i_trans_pt_blockv(planeVector, k, totalTargetLength + celli - 1)].store(vector);
                        // Loop over 3rd (vectorized) vspace dimension
                        for (uint i = 0; i < VECL; i++) {
                           targetBlockData[(totalTargetLength + celli) * WID3 +
                                           cellid_transpose[i + planeVector * VECL + k * WID2]]
                              = vector[i];
                        }
                     }
                  }
               }

               // store values from targetBlockData array to the actual blocks
               // Loop over cells in the pencil, including the padded cells of the target array
               for ( uint celli = 0; celli < targetLength; celli++ ) {
                  
                  uint GID = celli + totalTargetLength; 
                  SpatialCell* spatial_cell = targetCells[GID];
                  
                  if(spatial_cell == NULL) {
                     // Invalid target spatial cell
                     continue;
                  }
                  
                  const vmesh::LocalID blockLID = spatial_cell->get_velocity_block_local_id(blockGID, popID);
                  if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
                     // Invalid local id.
                     continue;
                  }

                  Realf* blockData = spatial_cell->get_data(blockLID, popID);
                  for(int i = 0; i < WID3 ; i++) {
                     blockData[i] += targetBlockData[GID * WID3 + i];
                  }
               }

               if(debugflag && debugPencilFlag) cout << "TotalTargetLength = " << totalTargetLength << endl;
               totalTargetLength += targetLength;
               
               // dealloc target data -- Should be automatic again?
            }
         }

         data = allCellsPointer[allCellsPointerIndex]->get_data(debugLID,popID);         
      }
   }   

   //cout << "I am at line " << __LINE__ << " of " << __FILE__ << endl;
  
   return true;
}
