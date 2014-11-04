/*
  This file is part of Vlasiator.
  Copyright 2014 Finnish Meteorological Institute
*/
#ifndef CPU_TRANS_MAP_H
#define CPU_TRANS_MAP_H
#include "vec4.h"
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

/*compute spatial neighbors for source stencil with a size of 2*
 * VLASOV_STENCIL_WIDTH + 1, cellID at VLASOV_STENCIL_WIDTH. First
 * bondary layer included. Invalid cells are replaced by closest good
 * cells (i.e. boundary condition uses constant extrapolation for the
 * stencil values at boundaries*/
  
void compute_spatial_source_neighbors(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,
                                      const uint dimension,
                                      CellID *neighbors){
   for(int i = -VLASOV_STENCIL_WIDTH; i <= VLASOV_STENCIL_WIDTH; i++){
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
   /*loop to neative side and replace all invalid cells with the closest good cell*/
   for(int i = -1;i>=-VLASOV_STENCIL_WIDTH;i--){
      if(neighbors[i + VLASOV_STENCIL_WIDTH] == INVALID_CELLID) 
         neighbors[i + VLASOV_STENCIL_WIDTH] = last_good_cellID;
      else
         last_good_cellID = neighbors[i + VLASOV_STENCIL_WIDTH];
   }

   last_good_cellID = cellID;
   /*loop to positive side and replace all invalid cells with the closest good cell*/
   for(int i = 1; i <= VLASOV_STENCIL_WIDTH; i++){
      if(neighbors[i + VLASOV_STENCIL_WIDTH] == INVALID_CELLID) 
         neighbors[i + VLASOV_STENCIL_WIDTH] = last_good_cellID;
      else
         last_good_cellID = neighbors[i + VLASOV_STENCIL_WIDTH];
   }
}

/*compute spatial target neighbors, stencil has a size of 3. No boundary cells are included*/
void compute_spatial_target_neighbors(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,
                                      const uint dimension,
                                      CellID *neighbors){

   for(int i = -1; i <= 1; i++){
      switch (dimension){
          case 0:
             neighbors[i + 1] = get_spatial_neighbor(mpiGrid, cellID, false, i, 0, 0);
             break;
          case 1:
             neighbors[i + 1] = get_spatial_neighbor(mpiGrid, cellID, false, 0, i, 0);
             break;
          case 2:
             neighbors[i + 1] = get_spatial_neighbor(mpiGrid, cellID, false, 0, 0, i);
             break;             
      }             
   }

}

/* Copy the fx data to the temporary values array, so that the
 * dimensions are correctly swapped. Also, copy the same block for
 * then neighboring spatial cells (in the dimension). neighbors
 * generated with compute_spatial_neighbors_wboundcond)*/

inline void copy_trans_block_data(
   const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const CellID cellID,
   const CellID* source_neighbors,
   const vmesh::GlobalID blockGID,Vec4* values,int dimension
) {
   uint cell_indices_to_id[3]={};
   switch (dimension){
       case 0:
          /* i and k coordinates have been swapped*/
          cell_indices_to_id[0]=WID2;
          cell_indices_to_id[1]=WID;
          cell_indices_to_id[2]=1;
          break;
       case 1:
          /* j and k coordinates have been swapped*/
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
   for (int b = -VLASOV_STENCIL_WIDTH; b <= VLASOV_STENCIL_WIDTH; ++b) {
      const CellID srcCell = source_neighbors[b + VLASOV_STENCIL_WIDTH];

      Realf* block_fx;
      const vmesh::LocalID blockLID = mpiGrid[srcCell]->get_velocity_block_local_id(blockGID);
      if (blockLID == mpiGrid[srcCell]->invalid_local_id()) {
         block_fx = mpiGrid[srcCell]->null_block_fx;
      } else {
         block_fx = mpiGrid[srcCell]->get_fx(blockLID);
      }

      /* Copy fx table, spatial source_neighbors already taken care of when
       *   creating source_neighbors table. If a normal spatial cell does not
       *   simply have the block, its value will be its null_block which
       *   is fine. This null_block has a value of zero in fx, and that
       *   is thus the velocity space boundary*/
      for (uint k=0; k<WID; ++k) {
         for (uint j=0; j<WID; ++j) {
            for (uint i=0; i<WID; ++i) {
               const uint cell =
                  i * cell_indices_to_id[0] +
                  j * cell_indices_to_id[1] +
                  k * cell_indices_to_id[2];
               /*copy data, when reading data from fx we swap dimensions using cell_indices_to_id*/
               values[i_trans_pblockv(b,j,k)].insert(i,(Real)block_fx[cell]);
            }
         }
      }
   }
}


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
  
*/
inline void store_trans_block_data(
   const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const CellID cellID, const CellID *target_neighbors,
   const vmesh::GlobalID blockGID,
   Vec4 * __restrict__ target_values,int dimension
) {
   uint cell_indices_to_id[3];
   
   switch (dimension){
       case 0:
          /* i and k coordinates have been swapped*/
          cell_indices_to_id[0]=WID2;
          cell_indices_to_id[1]=WID;
          cell_indices_to_id[2]=1;
          break;
       case 1:
          /* j and k coordinates have been swapped*/
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
               /*store data, when reading data from  data we swap dimensions using cell_indices_to_id*/
               block_data[cell] += target_values[i_trans_ptblockv(b,j,k)][i];
            }
         }
      }
   }
}

/*
  For local cells that are not boundary cells  block data is copied from data to fx, and data is
  set to zero, if boundary cell then   we copy from data to fx, but do not
  touch data. FOr remote cells fx is already up to data as we receive there.
*/
bool trans_prepare_block_data(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const CellID cellID){
   bool return_value=false;
   SpatialCell* spatial_cell = mpiGrid[cellID];   
   /*if we are on boundary then we do not set the data values to zero as these cells should not be updated*/
   const bool is_boundary = (spatial_cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY);
   /*if the cell is remote, then we do no copy data to the fx table, it should already have been set there*/
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
}





/* 
   Here we map from the current time step grid, to a target grid which
   is the lagrangian departure grid (so th grid at timestep +dt,
   tracked backwards by -dt). This is done in ordinary space in the translation step

   This function can, and should be, safely called in a parallel
OpenMP region (as long as it does only one dimension per parallel
refion). It is safe as each thread only computes certain blocks (blockID%tnum_threads = thread_num */

bool trans_map_1d(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const CellID cellID,const uint dimension, const Real dt) {
   /*values used with an stencil in 1 dimension, initialized to 0. Contains a block, and its spatial neighbours in one dimension */  
   Real dz,z_min, dvz,vz_min;
   SpatialCell* spatial_cell = mpiGrid[cellID];
   uint block_indices_to_id[3]; /*< used when computing id of target block */
   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/   
   uint thread_id = 0;  //thread id. Default value for serial case
   uint num_threads = 1; //Number of threads. Default value for serial case
   #ifdef _OPENMP
   //get actual values if OpenMP is enabled
   thread_id = omp_get_thread_num();
   num_threads = omp_get_num_threads(); 
   #endif
   //do nothing if it is not a normal cell, or a cell that is in the first boundary layer
   if (get_spatial_neighbor(mpiGrid, cellID, true, 0, 0, 0) == INVALID_CELLID) return true; 
   
   /*compute spatial neighbors, separately for targets and source. In
    * source cells we have a wider stencil and take into account
    * boundaries. For targets we only have actual cells as we do not
    * want to propagate boundary cells (array may contain
    * INVALID_CELLIDs at boundaries)*/
   CellID source_neighbors[1 + 2 * VLASOV_STENCIL_WIDTH];
   CellID target_neighbors[3];
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
   const Real i_dz=1.0/dz;

   /*Loop over blocks in spatial cell. In ordinary space the number of
    * blocks in this spatial cell does not change*/
   for (vmesh::LocalID block_i=0; block_i<spatial_cell->get_number_of_velocity_blocks(); ++block_i) {
      const vmesh::GlobalID blockGID = spatial_cell->get_velocity_block_global_id(block_i);

      //Each thread only computes a certain non-overlapping subset of blocks
      if (blockGID % num_threads != thread_id) continue;

      /*buffer where we write data, initialized to 0*/
      Vec4 target_values[3 * WID2];

      //init target_values
      for (uint i = 0; i<3*WID2; ++i) target_values[i] = Vec4(0.0, 0.0, 0.0, 0.0);

      /*buffer where we read in source data. i index vectorized*/
      Vec4 values[(1 + 2 * VLASOV_STENCIL_WIDTH) * WID3];
      copy_trans_block_data(mpiGrid, cellID, source_neighbors, blockGID, values, dimension);
      velocity_block_indices_t block_indices = SpatialCell::get_velocity_block_indices(blockGID);

      //i,j,k are now relative to the order in which we copied data to the values array. 
      //After this point in the k,j,i loops there should be no branches based on dimensions
      //
      //Note that the i dimension is vectorized, and thus there are no loops over i
      for (uint k=0; k<WID; ++k) {
         const Real cell_vz = (block_indices[dimension] * WID + k + 0.5) * dvz + vz_min; //cell centered velocity
         const Real z_translation = cell_vz * dt * i_dz; // how much it moved in time dt (reduced units)
         const int target_scell_index = (z_translation > 0) ? 1: -1; //part of density goes here (cell index change along spatial direcion)
         
         //the coordinates (scaled units from 0 to 1) between which we will
         //integrate to put mass in the target  neighboring cell. 
         //As we are below CFL<1, we know
         //that mass will go to two cells: current and the new one.
         Real z_1,z_2;
         if ( z_translation < 0 ) {
            z_1 = 0;
            z_2 = -z_translation; 
         }
         else {
            z_1 = 1.0 - z_translation;
            z_2 = 1.0;
         }
         
         for (uint j = 0; j < WID; ++j){ 
            //compute reconstruction
#ifdef TRANS_SEMILAG_PLM
            Vec4 a[3];
            compute_plm_coeff(values + i_trans_pblockv(-VLASOV_STENCIL_WIDTH , j, k), VLASOV_STENCIL_WIDTH, a);
#endif
#ifdef TRANS_SEMILAG_PPM
            Vec4 a[3];
            //Check that stencil width VLASOV_STENCIL_WIDTH in grid.h corresponds to order of face estimates  (h4 & h5 =2, H6=3, h8=4)
            compute_ppm_coeff(values + i_trans_pblockv(-VLASOV_STENCIL_WIDTH , j, k), h4, VLASOV_STENCIL_WIDTH, a);
#endif
#ifdef TRANS_SEMILAG_PQM
            Vec4 a[5];
            //Check that stencil width VLASOV_STENCIL_WIDTH in grid.h corresponds to order of face estimates (h4 & h5 =2, H6=3, h8=4)
            compute_pqm_coeff(values + i_trans_pblockv(-VLASOV_STENCIL_WIDTH , j, k), h6, VLASOV_STENCIL_WIDTH, a);
#endif
          
#ifdef TRANS_SEMILAG_PLM
            const Vec4 ngbr_target_density =
               (z_2 - z_1) * a[0] +
               (z_2 * z_2 - z_1 * z_1) * a[1];
#endif
#ifdef TRANS_SEMILAG_PPM
            const Vec4 ngbr_target_density =
               (z_2 - z_1) * a[0] +
               (z_2 * z_2 - z_1 * z_1) * a[1] +
               (z_2 * z_2 * z_2 - z_1 * z_1 * z_1) * a[2];
#endif
#ifdef TRANS_SEMILAG_PQM
            const Vec4 ngbr_target_density =
               (z_2                         - z_1                        ) * a[0] +
               (z_2 * z_2                   - z_1 * z_1                  ) * a[1] +
               (z_2 * z_2 * z_2             - z_1 * z_1 * z_1            ) * a[2] +
               (z_2 * z_2 * z_2 * z_2       - z_1 * z_1 * z_1 * z_1      ) * a[3] +
               (z_2 * z_2 * z_2 * z_2 * z_2 - z_1 * z_1 * z_1 * z_1 * z_1) * a[4];
#endif
            target_values[i_trans_ptblockv(target_scell_index,j,k)] +=  ngbr_target_density; //in the current original cells we will put this density        
            target_values[i_trans_ptblockv(0,j,k)] +=  values[i_trans_pblockv(0,j,k)] - ngbr_target_density; //in the current original cells we will put the rest of the original density
         }
      }

      //store values from target_values array to the actual blocks
      store_trans_block_data(mpiGrid,cellID,target_neighbors,blockGID,target_values,dimension);
   }

   return true;
}

/*!

  This function communicates the mapping on process boundaries, and then updates the data to their correct values.
  TODO, this could be inside an openmp region, in which case some m ore barriers and masters should be added

  \par dimension: 0,1,2 for x,y,z
  \par direction: 1 for + dir, -1 for - dir
*/
  
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
         //Send data in p_ngbr data array that we just
         //mapped to if 1) it is a valid target,
         //2) is remote cell, 3) if the source cell in center was
         //translated
         ccell->neighbor_block_data = &(pcell->get_data()[0]);
         ccell->neighbor_number_of_blocks = pcell->get_number_of_velocity_blocks();
         send_cells.push_back(p_ngbr);
      }
      
      if(m_ngbr != INVALID_CELLID &&
         !mpiGrid.is_local(m_ngbr) &&
         ccell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY 
         ){
         //Receive data that mcell mapped to ccell to this local cell
         //fx array, if 1) m is a valid source cell, 2) center cell is to be updated (normal cell) 3)  m is remote
         mcell->neighbor_block_data = &(ccell->get_fx()[0]);
         mcell->neighbor_number_of_blocks = ccell->get_number_of_velocity_blocks();
         receive_cells.push_back(local_cells[c]);
      }

   }
   //Do communication
   SpatialCell::set_mpi_transfer_type(Transfer::NEIGHBOR_VEL_BLOCK_FLUXES);

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
      //reduce data: sum received fx to data
      for (size_t c=0; c < receive_cells.size(); ++c) {
         SpatialCell *spatial_cell = mpiGrid[receive_cells[c]];      
#pragma omp for nowait
         for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(); cell++) {
            //copy data to fx for solvers, and set data to zero as we will map new values there
            spatial_cell->get_data()[cell] += spatial_cell->get_fx()[cell];
         }
      }

   
      /*send cell data is set to zero. This is to avoid double copy if
       * one cell is the neighbor on bot + and - side to the same
       * process*/
      for (size_t c=0; c < send_cells.size(); ++c) {
         SpatialCell *spatial_cell = mpiGrid[send_cells[c]];      
#pragma omp for nowait
         for(unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH * spatial_cell->get_number_of_velocity_blocks(); cell++) {
            spatial_cell->get_data()[cell] = 0.0;
         }
      }
   }
}

#endif
