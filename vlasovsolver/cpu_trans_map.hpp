
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
#include "cpu_1d_interpolations.hpp"

#ifdef TRANS_SEMILAG_PCONSTM
#define  TRANS_STENCIL_WIDTH 0
const int
#endif
#ifdef TRANS_SEMILAG_PLM
#define  TRANS_STENCIL_WIDTH 0
#endif
#if TRANS_SEMILAG_PPM
#define  TRANS_STENCIL_WIDTH 0
#endif



using namespace std;
using namespace spatial_cell;

//size of padded block where we store the block, and its spatial neighbors in one dimension
#define trans_pblock_xw  WID
#define trans_pblock_yw  WID
#define trans_pblock_zblockw  (WID * (1 + 2 * TRANS_STENCIL_WIDTH))    //how many blocks in z directions
#define trans_pblock_zw  (WID * trans_pblock_zblockw)  //how many cells in z direction

// indices in padded block. b_k is the block index in z direction in ordinary space, i,j,k are the cell ids inside on block. Data order i,b_k,j,k.
#define i_trans_pblock(b_k,i,j,k) ( (i) + ((b_k) + TRANS_STENCIL_WIDTH) * trans_pblock_xw + (j) * trans_pblock_xw * trans_pblock_zblockw  + (k) * trans_pblock_xw * trans_pblock_yw * trans_pblock_zblockw)
#define i_trans_pblockv(b_k,j,k)  ( i_trans_pblock(b_k,0,j,k) )

//size of xadded block where we write the target data
#define trans_ptblock_yw  WID
#define trans_ptblock_zblockw  3  //how many blocks in z directions
#define trans_ptblock_zw  (WID * trans_ptblock_zblockw)  //how many cells in z direction
// indices in padded  target block, which has Vec4 elements. b_k is the block index in z direction in ordinary space, i,j,k are the cell ids inside on block (i in vector elements). Data order i,b_k,j,k.
#define i_trans_ptblockv(b_k,j,k)  ( ((b_k) + 1) + (j) *  trans_ptblock_zblockw  + (k) * trans_ptblock_yw * trans_ptblock_zblockw)




/*return INVALID_CELLID if the spatial neighbor does not exist, or if
 * it is a cell that is not computed. If the
 * include_first_boundary_layer flag is set, then also first boundary
 * layer is inlcuded (does not return INVALID_CELLID).


 */

uint get_spatial_neighbor(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                          const CellID& cellID,
                          const bool include_first_boundary_layer,
                          const int spatial_di,
                          const int spatial_dj,
                          const int spatial_dk ) {
   if(spatial_di == 0 && spatial_dj == 0 && spatial_dk == 0 )
      return cellID;
   
   
   const std::vector<CellID> neighbors = mpiGrid.get_neighbors_of(cellID,spatial_di,spatial_dj,spatial_dk);
   if (neighbors.size() == 0) {
      std::cerr << __FILE__ << ":" << __LINE__
                << " No neighbor for cell " << cellID
                << " at offsets " << spatial_di << ", " << spatial_dj << ", " << spatial_dk
                << std::endl;
      abort();
  }

   // not existing cell or do not compute       
   if( neighbors[0] == INVALID_CELLID  || mpiGrid[neighbors[0]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE)
      return INVALID_CELLID; 

   //cell on boundary, but not first layer and we want to include first layer (e.g. when we compute source cells)
   if( include_first_boundary_layer &&
       mpiGrid[neighbors[0]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
       mpiGrid[neighbors[0]]->sysBoundaryLayer != 1 ) {
      return INVALID_CELLID;
   }
   //cell on boundary, and we want none of the layers, invalid.(e.g. when we compute targets)
   if( !include_first_boundary_layer &&
       mpiGrid[neighbors[0]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY){
      return INVALID_CELLID;
   }
   
   return neighbors[0]; //no AMR

}

/*compute spatial neighbors for source stencil with a size of 2*
 * TRANS_STENCIL_WIDTH + 1, cellID at TRANS_STENCIL_WIDTH. First
 * bondary layer included. Invalid cells are replaced by closest good
 * cells (i.e. boundary condition uses constant extrapolation for the
 * stencil values at boundaries*/
  
void compute_spatial_source_neighbors(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                               const CellID& cellID,
                               const uint dimension,
                               CellID *neighbors){
   for(int i = -TRANS_STENCIL_WIDTH; i <= TRANS_STENCIL_WIDTH; i++){
      switch (dimension){
          case 0:
             neighbors[i + TRANS_STENCIL_WIDTH] = get_spatial_neighbor(mpiGrid, cellID, true, i, 0, 0);
             break;
          case 1:
             neighbors[i + TRANS_STENCIL_WIDTH] = get_spatial_neighbor(mpiGrid, cellID, true, 0, i, 0);
             break;
          case 2:
             neighbors[i + TRANS_STENCIL_WIDTH] = get_spatial_neighbor(mpiGrid, cellID, true, 0, 0, i);
             break;             
      }             
   }

   CellID last_good_cellID = cellID;
   /*loop to neative side and replace all invalid cells with the closest good cell*/
   for(int i = -1;i>=-TRANS_STENCIL_WIDTH;i--){
      if(neighbors[i + TRANS_STENCIL_WIDTH] == INVALID_CELLID) 
         neighbors[i + TRANS_STENCIL_WIDTH] = last_good_cellID;
      else
         last_good_cellID = neighbors[i + TRANS_STENCIL_WIDTH];
   }

   last_good_cellID = cellID;
   /*loop to positive side and replace all invalid cells with the closest good cell*/
   for(int i = 1; i <= TRANS_STENCIL_WIDTH; i++){
      if(neighbors[i + TRANS_STENCIL_WIDTH] == INVALID_CELLID) 
         neighbors[i + TRANS_STENCIL_WIDTH] = last_good_cellID;
      else
         last_good_cellID = neighbors[i + TRANS_STENCIL_WIDTH];
   }

}

/*compute spatial target neighbors, stencil has a size of 3. No boundary cells are included*/
void compute_spatial_target_neighbors(const dccrg::Dccrg<SpatialCell>& mpiGrid,
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

inline void copy_trans_block_data(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CellID cellID, const CellID *source_neighbors, const uint blockID, Real * __restrict__ values, int dimension){
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
  }
  // Copy volume averages of this block from all spatial cells:
   for(int b = -TRANS_STENCIL_WIDTH; b <= TRANS_STENCIL_WIDTH; b++){
      Velocity_Block * block = mpiGrid[source_neighbors[b]]->at(blockID);
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
             values[i_trans_pblock(b,i,j,k)] = block->fx[cell];
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
inline void store_trans_block_data(const dccrg::Dccrg<SpatialCell>& mpiGrid, const CellID cellID, const CellID *target_neighbors, const uint blockID, Vec4 * __restrict__ target_values, int dimension){
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
  for (uint b=0; b<3; ++b) {
     if (target_neighbors[b] == INVALID_CELLID)
        continue; //do not store to boundary cells or otherwise invalid cells
     SpatialCell* spatial_cell = mpiGrid[target_neighbors[b]];
     Velocity_Block *block = spatial_cell->at(blockID);
     if (spatial_cell->is_null_block(block)) {
        /*block does not exist. If so, we do not create it and add stuff to it here. We
         * have already created blocks around blocks with content in
         * spatial sense, so we have no need to create even more blocks here*/
        /*TODO add loss counter*/
        continue;
     }
     for (uint k=0; k<WID; ++k) {
        for (uint j=0; j<WID; ++j) {
           for (uint i=0; i<WID; ++i) {
              const uint cell =
                 i * cell_indices_to_id[0] +
                 j * cell_indices_to_id[1] +
                 k * cell_indices_to_id[2];
              /*store data, when reading data from  data we swap dimensions using cell_indices_to_id*/
              block->data[cell] += target_values[i_trans_ptblockv(b,j,k)][i];
           }
        }
     }
  }
}


/*
  For cells that are not boundary cells  block data is copied from data to fx, and data is
  set to zero, if boundary cell then   we copy from data to fx, but do not
  touch data.

  TODO: MPI communication and boundary conditions could be made smarter to avoid these extra copies.
  TODO: parallelize with OpenMP
*/
bool trans_prepare_block_data(const dccrg::Dccrg<SpatialCell>& mpiGrid, const CellID cellID){
  /* 
     Move densities from data to fx and clear data, to prepare for mapping
  */
   SpatialCell* spatial_cell = mpiGrid[cellID];   
   const bool clear_data = (spatial_cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY);

   
   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
      const unsigned int block = spatial_cell->velocity_block_list[block_i];
      Velocity_Block * __restrict__ block_ptr = spatial_cell->at(block);
      Real * __restrict__ fx = block_ptr->fx;
      Real * __restrict__ data = block_ptr->data;
      
      for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
         fx[cell] = data[cell];
         if(clear_data)
            data[cell] = 0.0;
      }
   }

}





/* 
   Here we map from the current time step grid, to a target grid which
   is the lagrangian departure grid (so th grid at timestep +dt,
   tracked backwards by -dt). This is done in ordinary space in the translation step

   TODO: parallelize with OpenMP over blocks
*/

bool trans_map_1d(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CellID cellID,const uint dimension, const Real dt) {
   uint count = 0;
//   cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
   /*values used with an stencil in 1 dimension, initialized to 0. Contains a block, and its spatial neighbours in one dimension */  
   Real dz,z_min, dvz,vz_min;
   SpatialCell* spatial_cell = mpiGrid[cellID];
   uint block_indices_to_id[3]; /*< used when computing id of target block */
   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/   


   if(dimension>2)
      return false; //not possible


   /*compute spatial neighbors, separately for targets and source. In
    * source cells we have a wider stencil and take into account
    * boundaries. For targets we only have actual cells as we do not
    * want to propagate boundary cells (array may contain
    * INVALID_CELLIDs at boundaries)*/
   CellID source_neighbors[1 + 2 * TRANS_STENCIL_WIDTH];
   CellID target_neighbors[3];
//   cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
   compute_spatial_source_neighbors(mpiGrid,cellID,dimension,source_neighbors);
//   cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
   compute_spatial_target_neighbors(mpiGrid,cellID,dimension,target_neighbors); 
//   cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;   
  /*set cell size in dimension direction*/  
   switch (dimension){
      case 0:
    dz = P::dx_ini;
    z_min = P::xmin;
    dvz = SpatialCell::cell_dvx;
    vz_min = SpatialCell::vx_min;
    block_indices_to_id[0]=SpatialCell::vx_length * SpatialCell::vy_length;
    block_indices_to_id[1]=SpatialCell::vx_length;
    block_indices_to_id[2]=1;
    /*set values in array that is used to transfer blockindices to id using a dot product*/
    cell_indices_to_id[0]=WID2;
    cell_indices_to_id[1]=WID;
    cell_indices_to_id[2]=1;

    break;
  case 1:
    dz = P::dy_ini;
    z_min = P::ymin;
    dvz = SpatialCell::cell_dvy;
    vz_min = SpatialCell::vy_min;
    /*set values in array that is used to transfer blockindices to id using a dot product*/
    block_indices_to_id[0]=1;
    block_indices_to_id[1]=SpatialCell::vx_length * SpatialCell::vy_length;
    block_indices_to_id[2]=SpatialCell::vx_length;
    /*set values in array that is used to transfer blockindices to id using a dot product*/
    cell_indices_to_id[0]=1;
    cell_indices_to_id[1]=WID2;
    cell_indices_to_id[2]=WID;

    break;
  case 2:
    dz = P::dz_ini;
    z_min = P::zmin;
    dvz = SpatialCell::cell_dvz;
    vz_min = SpatialCell::vz_min;
    /*set values in array that is used to transfer blockindices to id using a dot product*/
    block_indices_to_id[0]=1;
    block_indices_to_id[1]=SpatialCell::vx_length;
    block_indices_to_id[2]=SpatialCell::vx_length * SpatialCell::vy_length;
    /*set values in array that is used to transfer blockindices to id using a dot product*/
    cell_indices_to_id[0]=1;
    cell_indices_to_id[1]=WID;
    cell_indices_to_id[2]=WID2;

    
    break;
  }
  const Real i_dz=1.0/dz;

  
  /*Loop over blocks in spatial cell. In ordinary space the number of
   * blocks in this spatial cell does not change, blocks in
   * neighboring cells might*/
//  cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
  
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
//     cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;  
     const unsigned int blockID = spatial_cell->velocity_block_list[block_i];
    Velocity_Block * __restrict__ block = spatial_cell->at(blockID);

    Real values[trans_pblock_xw * trans_pblock_yw * trans_pblock_zw];
    Vec4 target_values[trans_ptblock_yw * trans_ptblock_zw]={}; /*buffer where we write data, initialized to 0*/
//    cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;  
    copy_trans_block_data(mpiGrid,cellID,source_neighbors,blockID,values,dimension);
//   cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
    velocity_block_indices_t block_indices=SpatialCell::get_velocity_block_indices(blockID);
    
    /*i,j,k are now relative to the order in which we copied data to the values array. 
      After this point in the k,j,i loops there should be no branches based on dimensions
      
      Note that the i dimension is vectorized, and thus there are no loops over i
    */

    for (uint k=0; k<WID; ++k){
       const Real cell_vz = (block_indices[dimension] * WID + k + 0.5) * dvz + vz_min; //cell centered velocity
       const Real z_translation = - cell_vz * dt * i_dz; // how much it moved in time - dt (departure grid) (in units of dz)
       const int target_scell_index = (z_translation > 0) ? 1: -1; //part of density goes here (cell index change along spatial direcion)

       for (uint j = 0; j < WID; ++j){ 
          /*compute reconstruction*/
//          cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
#ifdef TRANS_SEMILAG_PLM
          Vec4 a[2];
          Vec4 mv,cv,pv;
          mv.load(values + i_trans_pblockv(-1,j,k));
          cv.load(values + i_trans_pblockv(0,j,k));
          pv.load(values + i_trans_pblockv(1,j,k));
          compute_plm_coeff(mv,cv,pv,a);
#endif
#ifdef TRANS_SEMILAG_PPM
          Vec4 a[3];
          Vec4 mmv,mv,cv,pv,ppv;
          mmv.load(values + i_trans_pblockv(-2,j,k));
          mv.load(values + i_trans_pblockv(-1,j,k));
          cv.load(values + i_trans_pblockv(0,j,k));
          pv.load(values + i_trans_pblockv(1,j,k));
          ppv.load(values + i_trans_pblockv(2,j,k));
          compute_ppm_coeff(mmv,mv,cv,pv,ppv,a);
#endif
//             cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
	  
	  //the velocity between which we will integrate to put mass
	  //in the targe cell. z_translation defines the departure grid. As we are below CFL<1, we know that mass will go to two cells: current and the new one. 
          Real z_1,z_2;
          if ( z_translation < 0 ) {
             z_1 = 0.5 + z_translation;
             z_2 = 0.5;
          }
          else {
             z_1 = -0.5;
             z_2 = -0.5 + z_translation;
          }
          
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
//          cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
          target_values[i_trans_ptblockv(target_scell_index,j,k)] +=  ngbr_target_density; //in the current original cells we will put this density        
          target_values[i_trans_ptblockv(0,j,k)] +=  cv - ngbr_target_density; //in the current original cells we will put the rest of the original density
//          cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
       }
    }
//    cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
    //store values from target_values array to the actual blocks
    store_trans_block_data(mpiGrid, cellID, target_neighbors, blockID,target_values,dimension);
//    cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
  }

//  cout<< "Count "<< count++<<" at " << __FILE__ << ":"<<__LINE__<<endl;
  return true;
}


#endif   
