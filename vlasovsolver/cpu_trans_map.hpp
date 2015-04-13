/*
  This file is part of Vlasiator.
  Copyright 2014 Finnish Meteorological Institute
*/
#ifndef CPU_TRANS_MAP_H
#define CPU_TRANS_MAP_H
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


//Is cell translated? It is not translated if DO_NO_COMPUTE or if it is sysboundary cell and not in first sysboundarylayer
bool do_translate_cell(SpatialCell* SC){
   if(SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
      (SC->sysBoundaryLayer != 1 && SC->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY))
      return false;
   else
      return true;
}



//Create temporary target grid where we write the mapped values for
//all cells in cells vector. In this non-AMR version it contains the
//same blocks as the normal grid.
void createTargetGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,   const vector<CellID>& cells){
   phiprof::start("create-target-grid");
#pragma omp  parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell *spatial_cell = mpiGrid[cells[c]];
      /*get target mesh & blocks (in temporary arrays)*/

      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = spatial_cell->get_velocity_mesh_temporary();
      //copy mesh (population 0)
      vmesh = spatial_cell->get_velocity_mesh(0);

      vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks_temporary();
      // allocate space in block container and set block parameters
      blockContainer.clear();
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
   phiprof::stop("create-target-grid");
}

//Clear temporary target grid  for all cells in cells vector.
void clearTargetGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,   const vector<CellID>& cells){
   phiprof::start("clear-target-grid");
#pragma omp  parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell *spatial_cell = mpiGrid[cells[c]];
      spatial_cell->get_velocity_mesh_temporary().clear();
      spatial_cell->get_velocity_blocks_temporary().clear();
   }
   phiprof::stop("clear-target-grid");

         
}

//Set all values in the temporary target grid to zero (0.0), for all cells in cells vector.   
void zeroTargetGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,   const vector<CellID>& cells){
   phiprof::start("zero-target-grid");      
#pragma omp  parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell *spatial_cell = mpiGrid[cells[c]];
      vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks_temporary();
      for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * blockContainer.size(); cell++) {
         blockContainer.getData()[cell] = 0;
      }
   }
   phiprof::stop("zero-target-grid");
}


//Swap temporary target grid and normal grid. This is cheap as values are not copied.
void swapTargetSourceGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,   const vector<CellID>& cells){
   phiprof::start("swap-target-grid");
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell *spatial_cell = mpiGrid[cells[c]];
      spatial_cell->swap(spatial_cell->get_velocity_mesh_temporary(), spatial_cell->get_velocity_blocks_temporary());
   }
   phiprof::stop("swap-target-grid");

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

/* Copy the data to the temporary values array, so that the
 * dimensions are correctly swapped. Also, copy the same block for
 * then neighboring spatial cells (in the dimension). neighbors
 * generated with compute_spatial_neighbors_wboundcond)*/

inline void copy_trans_block_data(
   SpatialCell** source_neighbors,
   const vmesh::GlobalID blockGID,
   Vec* values,
   const unsigned char * const cellid_transpose) {


   //  Copy volume averages of this block from all spatial cells:
   for (int b = -VLASOV_STENCIL_WIDTH; b <= VLASOV_STENCIL_WIDTH; ++b) {
      SpatialCell* srcCell = source_neighbors[b + VLASOV_STENCIL_WIDTH];
      const vmesh::LocalID blockLID = srcCell->get_velocity_block_local_id(blockGID);
      if (blockLID != srcCell->invalid_local_id()) {
         Realv blockValues[WID3];
         Realf* block_data = srcCell->get_data(blockLID);
         /* Copy data to a temporary array and transpose values so that mapping is along k direction.
          * spatial source_neighbors already taken care of when
          * creating source_neighbors table. If a normal spatial cell does not
          * simply have the block, its value will be its null_block which
          * is fine. This null_block has a value of zero in data, and that
          * is thus the velocity space boundary*/
         for (uint i=0; i<WID3; ++i) {
            blockValues[i] = block_data[cellid_transpose[i]];
         }

         /*now load values into the actual values table..*/
         uint offset =0;
         for (uint k=0; k<WID; ++k) {
            for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
               /*store data, when reading data from  data we swap dimensions using precomputed plane_index_to_id and cell_indices_to_id*/
               values[i_trans_ps_blockv(planeVector, k, b)].load(blockValues + offset);
               offset += VECL;
            }
         }
      }
      else {
         uint cellid=0;
         for (uint k=0; k<WID; ++k) {
            for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
               values[i_trans_ps_blockv(planeVector, k, b)] = Vec(0);
            }
         }
      }
   }
}


/*!
  Store values to target grid from the new target data we have computed
  
  For dimension=0  we have rotated data
  i -> k
  j -> j
  k -> i
  For dimension=1   we have rotated data
  i -> i
  j -> k
  k -> j
  
*/
inline void store_trans_block_data(
   SpatialCell **target_neighbors,
   const vmesh::GlobalID blockGID,
   Vec * __restrict__ target_values,
   const unsigned char * const cellid_transpose
) {
   //Store volume averages in target blocks:
   for (int b=-1; b<=1; ++b) {
      if (target_neighbors[b + 1] == INVALID_CELLID) {
         continue; //do not store to boundary cells or otherwise invalid cells
      }
      SpatialCell* spatial_cell = target_neighbors[b + 1];
      const vmesh::LocalID blockLID = spatial_cell->get_velocity_block_local_id(blockGID);
      if (blockLID == spatial_cell->invalid_local_id()) {
         // block does not exist. If so, we do not create it and add stuff to it here.
         // We have already created blocks around blocks with content in
         // spatial sense, so we have no need to create even more blocks here
         // TODO add loss counter
         continue;
      }

      /*get block container for target cells*/
      vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks_temporary();
      Realf* block_data = blockContainer.getData(blockLID);

      Realv blockValues[VECL];
      uint cellid=0;
      for (uint k=0; k<WID; ++k) {
         for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
            target_values[i_trans_pt_blockv(planeVector, k, b)].store(blockValues);
            for(uint i = 0; i< VECL; i++){
               /*store data, when reading data from  data we swap dimensions using precomputed plane_index_to_id and cell_indices_to_id*/
               block_data[cellid_transpose[cellid++]] += blockValues[i];
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


bool trans_map_1d(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const CellID cellID,const uint dimension, const Realv dt) {
   /*values used with an stencil in 1 dimension, initialized to 0. Contains a block, and its spatial neighbours in one dimension */  
   Realv dz,z_min, dvz,vz_min;
   SpatialCell* spatial_cell = mpiGrid[cellID];
   if (spatial_cell->get_number_of_velocity_blocks() == 0) return true;
   
   uint block_indices_to_id[3]; /*< used when computing id of target block */
   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
   unsigned char  cellid_transpose[WID3]; /*< defines the transpose for the solver internal (transposed) id: i + j*WID + k*WID2 to actual one*/
   uint thread_id = 0;  //thread id. Default value for serial case
   uint num_threads = 1; //Number of threads. Default value for serial case
   #ifdef _OPENMP
   //get actual values if OpenMP is enabled
   thread_id = omp_get_thread_num();
   num_threads = omp_get_num_threads(); 
   #endif
   //do nothing if it is not a normal cell, or a cell that is in the first boundary layer
   if (get_spatial_neighbor(mpiGrid, cellID, true, 0, 0, 0) == INVALID_CELLID)
      return true; 
   
   /*compute spatial neighbors, separately for targets and source. In
    * source cells we have a wider stencil and take into account
    * boundaries. For targets we only have actual cells as we do not
    * want to propagate boundary cells (array may contain
    * INVALID_CELLIDs at boundaries)*/
   SpatialCell* source_neighbors[1 + 2 * VLASOV_STENCIL_WIDTH];
   SpatialCell* target_neighbors[3];
   compute_spatial_source_neighbors(mpiGrid,cellID,dimension,source_neighbors);
   compute_spatial_target_neighbors(mpiGrid,cellID,dimension,target_neighbors);
   
   /*set cell size in dimension direction*/
   switch (dimension){
       case 0:
          dz = P::dx_ini;
          z_min = P::xmin;
          dvz = SpatialCell::get_velocity_grid_cell_size()[0];
          vz_min = SpatialCell::get_velocity_grid_min_limits()[0];
          block_indices_to_id[0]=SpatialCell::get_velocity_grid_length()[0]*SpatialCell::get_velocity_grid_length()[1];
          block_indices_to_id[1]=SpatialCell::get_velocity_grid_length()[0];
          block_indices_to_id[2]=1;
      
          /*set values in array that is used to transfer blockindices to id using a dot product*/
          cell_indices_to_id[0]=WID2;
          cell_indices_to_id[1]=WID;
          cell_indices_to_id[2]=1;

          break;
       case 1:
          dz = P::dy_ini;
          z_min = P::ymin;
          dvz = SpatialCell::get_velocity_grid_cell_size()[1];
          vz_min = SpatialCell::get_velocity_grid_min_limits()[1];
      
          /*set values in array that is used to transfer blockindices to id using a dot product*/
          block_indices_to_id[0]=1;
          block_indices_to_id[1]=SpatialCell::get_velocity_grid_length()[0]*SpatialCell::get_velocity_grid_length()[1];
          block_indices_to_id[2]=SpatialCell::get_velocity_grid_length()[0];

          /*set values in array that is used to transfer blockindices to id using a dot product*/
          cell_indices_to_id[0]=1;
          cell_indices_to_id[1]=WID2;
          cell_indices_to_id[2]=WID;

          break;
       case 2:
          dz = P::dz_ini;
          z_min = P::zmin;
          dvz = SpatialCell::get_velocity_grid_cell_size()[2];
          vz_min = SpatialCell::get_velocity_grid_min_limits()[2];    
      
          /*set values in array that is used to transfer blockindices to id using a dot product*/
          block_indices_to_id[0]=1;
          block_indices_to_id[1]=SpatialCell::get_velocity_grid_length()[0];
          block_indices_to_id[2]=SpatialCell::get_velocity_grid_length()[0]*SpatialCell::get_velocity_grid_length()[1];

          /*set values in array that is used to transfer blockindices to id using a dot product*/
          cell_indices_to_id[0]=1;
          cell_indices_to_id[1]=WID;
          cell_indices_to_id[2]=WID2;
          break;

       default:
          cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
          abort();
          break;
   }

   /*init plane_index_to_id*/
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

   /*Loop over blocks in spatial cell. In ordinary space the number of
    * blocks in this spatial cell does not change*/
#pragma omp for
   for (vmesh::LocalID block_i=0; block_i<spatial_cell->get_number_of_velocity_blocks(); ++block_i) {
      const vmesh::GlobalID blockGID = spatial_cell->get_velocity_block_global_id(block_i);

      //Each thread only computes a certain non-overlapping subset of blocks
      //      if (blockGID % num_threads != thread_id) continue;

      /*buffer where we write data, initialized to 0*/
      Vec target_values[3 * WID3 / VECL];

      //init target_values
      for (uint i = 0; i< 3 * WID3 / VECL; ++i) {
         target_values[i] = Vec(0.0);
      }

      /*buffer where we read in source data. i index vectorized*/
      Vec values[(1 + 2 * VLASOV_STENCIL_WIDTH) * WID3 / VECL];
      copy_trans_block_data(source_neighbors, blockGID, values, cellid_transpose);
      velocity_block_indices_t block_indices = SpatialCell::get_velocity_block_indices(blockGID);
      
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
         }
         else {
            z_1 = 1.0 - z_translation;
            z_2 = 1.0;
         }
         for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){         
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
            target_values[i_trans_pt_blockv(planeVector, k, target_scell_index)] +=  ngbr_target_density;                     //in the current original cells we will put this density        
            target_values[i_trans_pt_blockv(planeVector, k, 0)] +=  values[i_trans_ps_blockv(planeVector, k, 0)] - ngbr_target_density; //in the current original cells we will put the rest of the original density
         }
      }
      
      //store values from target_values array to the actual blocks
      store_trans_block_data(target_neighbors,blockGID,target_values,cellid_transpose);
   }
   
   return true;
}

/*!

  This function communicates the mapping on process boundaries, and then updates the data to their correct values.
  TODO, this could be inside an openmp region, in which case some more barriers and masters should be added

  \par dimension: 0,1,2 for x,y,z
  \par direction: 1 for + dir, -1 for - dir
*/






void update_remote_mapping_contribution(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const uint dimension, int direction) {
   const vector<CellID>& local_cells = getLocalCells();
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
      }
      if(mpiGrid.is_local(p_ngbr) && mpiGrid.is_local(m_ngbr))
         continue; //internal cell, not much to do
      SpatialCell *pcell = NULL;
      if(p_ngbr != INVALID_CELLID)
         pcell = mpiGrid[p_ngbr];
      SpatialCell *mcell = NULL;
      if(m_ngbr != INVALID_CELLID)
         mcell = mpiGrid[m_ngbr];
      if(p_ngbr != INVALID_CELLID &&
         !mpiGrid.is_local(p_ngbr) &&
         do_translate_cell(ccell)
         ) {
         //Send data in p_ngbr temporary target array that we just
         //mapped to if 1) it is a valid target,
         //2) is remote cell, 3) if the source cell in center was
         //translated
         vmesh::VelocityBlockContainer<vmesh::LocalID>& pcellBlockContainer = pcell->get_velocity_blocks_temporary();
         ccell->neighbor_block_data = &(pcellBlockContainer.getData()[0]);
         ccell->neighbor_number_of_blocks = pcellBlockContainer.size();
         send_cells.push_back(p_ngbr);
      }
      if(m_ngbr != INVALID_CELLID &&
         !mpiGrid.is_local(m_ngbr) &&
         ccell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY
         ){
         //Receive data that mcell mapped to ccell to this local cell
         //data array, if 1) m is a valid source cell, 2) center cell is to be updated (normal cell) 3) m is remote
         //we can reuse the normal data array as we do not anymore need the original distribution function values
         mcell->neighbor_block_data = &(ccell->get_data()[0]);
         mcell->neighbor_number_of_blocks = ccell->get_number_of_velocity_blocks();
         receive_cells.push_back(local_cells[c]);
      }
   }
   
   //Do communication
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
      //reduce data: sum received data in the data array to the target grid in the temporary block contaniner
      for (size_t c=0; c < receive_cells.size(); ++c) {
         SpatialCell *spatial_cell = mpiGrid[receive_cells[c]];
         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks_temporary();
         #pragma omp for nowait
         for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(); cell++) {
            //copy received target data to temporary array where target data is stored.
            blockContainer.getData()[cell] += spatial_cell->get_data()[cell];
         }
      }
      /*send cell data is set to zero. This is to avoid double copy if
       * one cell is the neighbor on bot + and - side to the same
       * process*/
      for (size_t c=0; c < send_cells.size(); ++c) {
         SpatialCell *spatial_cell = mpiGrid[send_cells[c]];
         vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks_temporary();
         #pragma omp for nowait
         for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * blockContainer.size(); cell++) {
            blockContainer.getData()[cell] = 0.0;
         }
      }
   }
}



#endif
