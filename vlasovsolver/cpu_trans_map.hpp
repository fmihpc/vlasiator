
/*
This file is part of Vlasiator.
Copyright 2012 Finnish Meteorological Institute

*/

#ifndef CPU_TRANS_MAP_H
#define CPU_TRANS_MAP_H

#include  "vec4.h"
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
#define pblock_xw  WID
#define pblock_yw  WID
#define pblock_zw  (1 + 2 * TRANS_STENCIL_WIDTH)
// indices in padded block. b_k is the block index in z direction in ordinary space, i,j,k are the cell ids inside on block. Data order i,b_k,j,k.
#define i_pblock(b_k,i,j,k) ( (i) + ((b_k) + TRANS_STENCIL_WIDTH) * pblock_xw + (j) * pblock_xw * pblock_zw  + (k) * pblock_xw * pblock_yw * pblock_zw)
#define i_pblockv(b_k,j,k)  ( ((b_k) + TRANS_STENCIL_WIDTH) * pblock_xw + (j) * pblock_xw * pblock_zw  + (k) * pblock_xw * pblock_yw * pblock_zw)

//size of padded block where we write the target data
#define ptblock_xw  WID 
#define ptblock_yw  WID
#define ptblock_zw  3
// indices in padded  target block, which has Vec4 elements. b_k is the block index in z direction in ordinary space, i,j,k are the cell ids inside on block (i in vector elements). Data order i,b_k,j,k.
#define i_ptblockv(b_k,j,k)  ( ((b_k) + 1) + (j) *  ptblock_zw  + (k) * ptblock_yw * ptblock_zw)



Velocity_Block * getSpatialNeighbourBlock(
			 const dccrg::Dccrg<SpatialCell>& mpiGrid,
			 const CellID& cellID,
			 const uint& blockID,
			 const int spatial_di,
			 const int spatial_dj,
			 const int spatial_dk ) {
  
   const std::vector<CellID> neighbors = mpiGrid.get_neighbors_of(cellID,spatial_di,spatial_dj,spatial_dk);
   if (neighbors.size() == 0) {
      std::cerr << __FILE__ << ":" << __LINE__
                << " No neighbor for cell " << cellID
                << " at offsets " << int(i) - 2 << ", " << int(j) - 2 << ", " << int(k) - 2
	      << std::endl;
      abort();
  }
  
  if( neighbors[0] == INVALID_CELLID  ||
      mpiGrid[neighbors[0]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
      (mpiGrid[neighbors[0]]->sysBoundaryLayer != 1  &&
       mpiGrid[neighbors[0]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY )
      ) {
     return NULL; //use whatever ngbr condition we need to do in copy function
  } 
  else {
    const CellID ngbrCellID=neighbors[0]; //no AMR
    SpatialCell* spatial_cell = mpiGrid[ngbrCellID];
    return spatial_cell->at(blockID); /*if it does not exits, its null_block will be returned*/
  }
}

/*!
  For dimension=0 data copy  we have rotated data
  i -> k
  j -> j
  k -> i
For dimension=1 data copy  we have rotated data
  i -> i
  j -> k
  k -> j

*/

/*FIXME, use stencil and not hardcoded eighborhood*/
inline void copy_trans_block_data(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CellID cellID,const uint blockID,Real * __restrict__ values, int dimension){
  SpatialCell* spatial_cell = mpiGrid[cellID];
  Velocity_Block *blocks[pblock_zw];
  uint cell_indices_to_id[3];
  switch (dimension){
      case 0:
    /* i and k coordinates have been swapped*/
         cell_indices_to_id[0]=WID2;
         cell_indices_to_id[1]=WID;
         cell_indices_to_id[2]=1;
         for (int bk=-TRANS_STENCIL_WIDTH; bk <= TRANS_STENCIL_WIDTH;bk++) {
            if(bk!=0)
               blocks[bk + TRANS_STENCIL_WIDTH] = getSpatialNeighbourBlock(mpiGrid, cellID, blockID, bk, 0, 0);
            else
               blocks[bk + TRANS_STENCIL_WIDTH] = spatial_cell->at(blockID);
         }
         break;
      case 1:
         /* j and k coordinates have been swapped*/
         cell_indices_to_id[0]=1;
         cell_indices_to_id[1]=WID2;
         cell_indices_to_id[2]=WID;
         for (int bk=-TRANS_STENCIL_WIDTH; bk <= TRANS_STENCIL_WIDTH;bk++) {
            if(bk!=0)
          blocks[bk + TRANS_STENCIL_WIDTH] = getSpatialNeighbourBlock(mpiGrid, cellID, blockID, 0, bk, 0);
            else
               blocks[bk + TRANS_STENCIL_WIDTH] = spatial_cell->at(blockID);
         }
         break;
      case 2:
         cell_indices_to_id[0]=1;
         cell_indices_to_id[1]=WID;
         cell_indices_to_id[2]=WID2;
         for (int bk=-TRANS_STENCIL_WIDTH; bk <= TRANS_STENCIL_WIDTH;bk++) {
            if(bk!=0)
               blocks[bk + TRANS_STENCIL_WIDTH] = getSpatialNeighbourBlock(mpiGrid, cellID, blockID, 0, 0, bk);
            else
               blocks[bk + TRANS_STENCIL_WIDTH] = spatial_cell->at(blockID);
         }
         break;
  }
  
  // Copy volume averages of this block:
  for (uint b=0; b<pblock_zw; ++b) {
    Real * __restrict__ fx;

    /*get closest possible block (on boundaries some cells have no
     * blocks). If a normal spatial cell does not simply have the block, its value
     * will be its null_block which is fine. This null_block has a
     * value of zero in fx, and that is thus the velocity space
     * boundary*/
    if (b == STENCIL_WIDTH)
       fx = blocks[b]->fx;
    else if (b<STENCIL_WIDTH)
       /*copy closest existing block for spatial cell on the negative side*/
       for (db = 0; b + db <= STENCIL_WIDTH;db ++) {
          if ( blocks[b + db] !=NULL) {
             fx = blocks[b+ db]->fx;
             break;
          }
       }

    else if (b<STENCIL_WIDTH)
       /*copy closest existing block for spatial cell on the positive side*/
       for (db = 0; b - db <= STENCIL_WIDTH;db ++) {
          if ( blocks[b - db] !=NULL) {
             fx = blocks[b - db]->fx;
             break;
          }
       }
    //we have our fx table from the closest possible block, copy data now
    for (uint k=0; k<WID; ++k) {
      for (uint j=0; j<WID; ++j) {
	for (uint i=0; i<WID; ++i) {
	  const uint cell =
	    i * cell_indices_to_id[0] +
	    j * cell_indices_to_id[1] +
	    k * cell_indices_to_id[2];
	  /*copy data, when reading data from fx we swap dimensions using cell_indices_to_id*/
	  values[i_pblock(b,i,j,k)] = fx[cell];
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
inline void store_trans_block_data(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CellID cellID,const uint blockID,Vec4 * __restrict__ target_values, int dimension){
  SpatialCell* spatial_cell = mpiGrid[cellID];
  Velocity_Block *blocks[3];
  uint cell_indices_to_id[3];
  Velocity_Block * block_p,block_m,block_pp,block_mm;
  
  blocks[1] = spatial_cell->at(blockID);
  switch (dimension){
  case 0:
    /* i and k coordinates have been swapped*/
    cell_indices_to_id[0]=WID2;
    cell_indices_to_id[1]=WID;
    cell_indices_to_id[2]=1;
    blocks[0] = getSpatialNeighbourBlock(mpiGrid, cellID, blockID, -1, 0, 0);
    blocks[2] = getSpatialNeighbourBlock(mpiGrid, cellID, blockID, 1, 0, 0);
    break;
  case 1:
    /* j and k coordinates have been swapped*/
    cell_indices_to_id[0]=1;
    cell_indices_to_id[1]=WID2;
    cell_indices_to_id[2]=WID;
    blocks[0] = getSpatialNeighbourBlock(mpiGrid, cellID, blockID, 0, -1, 0);
    blocks[2] = getSpatialNeighbourBlock(mpiGrid, cellID, blockID, 0, 1, 0);
    break;
  case 2:
    cell_indices_to_id[0]=1;
    cell_indices_to_id[1]=WID;
    cell_indices_to_id[2]=WID2;
    blocks[0] = getSpatialNeighbourBlock(mpiGrid, cellID, blockID, 0, 0, -1);
    blocks[2] = getSpatialNeighbourBlock(mpiGrid, cellID, blockID, 0, 0, 1);
    break;
  }
  
  
  //Store volume averages in target blocks:
  for (uint b=0; b<3; ++b) {
    Real * __restrict__ data;
    if (blocks[b] != NULL) {
       
       if (spatial_cell->count(tblock) == 0) {
          //count faster since the add_velocity_block call is more expensive
          spatial_cell->add_velocity_block(tblock);
       }
       Velocity_Block* block_ptr = spatial_cell->at_fast(tblock);
       target_block_data=block_ptr->data;
       
          
       for (uint k=0; k<WID; ++k) {
          for (uint j=0; j<WID; ++j) {
             for (uint i=0; i<WID; ++i) {
                const uint cell =
                   i * cell_indices_to_id[0] +
                   j * cell_indices_to_id[1] +
                   k * cell_indices_to_id[2];
                /*copy data, when reading data from    data we swap dimensions using cell_indices_to_id*/
                
                blocks[b]->data[cell] = target_values[i_pblock(b,j,k)][i];
             }
          }
       }
    }
  }
}


/* 
   Here we map from the current time step grid, to a target grid which
   is the lagrangian departure grid (so th grid at timestep +dt,
   tracked backwards by -dt). This is done in ordinary space in the translation step
*/

bool map_trans_1d(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CellID cellID,const uint dimension, const Real dt) {
   /*values used with an stencil in 1 dimension, initialized to 0. Contains a block, and its spatial neighbours in one dimension */  
   Real dz,z_min, dvz,vz_min;
   SpatialCell* spatial_cell = mpiGrid[cellID];
   uint block_indices_to_id[3]; /*< used when computing id of target block */
   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/   
   if(dimension>2)
      return false; //not possible
   
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

  
  /* 
     Move densities from data to fx and clear data, to prepare for mapping
  */
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
    const unsigned int block = spatial_cell->velocity_block_list[block_i];
    Velocity_Block * __restrict__ block_ptr = spatial_cell->at(block);
    Real * __restrict__ fx = block_ptr->fx;
    Real * __restrict__ data = block_ptr->data;
    
    for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
      fx[cell] = data[cell];
      data[cell] = 0.0;
    }
  }
  
  
  

  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
    const unsigned int blockID = spatial_cell->velocity_block_list[block_i];
    Velocity_Block * __restrict__ block = spatial_cell->at(block);

    Real values[pblock_xw * pblock_yw * pblock_zw];
    Vec4 target_values[ptblock_yw * ptblock_zw]={}; /*buffer where we write data, initialized to 0*/

    copy_block_data(mpiGrid,cellID,blockID,values,dimension);

    velocity_block_indices_t block_indices=SpatialCell::get_velocity_block_indices(blockID);
    
    /*i,j,k are now relative to the order in which we copied data to the values array. 
      After this point in the k,j,i loops there should be no branches based on dimensions
      
      Note that the i dimension is vectorized, and thus there are no loops over i
    */

    for (uint k=0; k<WID; ++k){
       const Real cell_vz = (block_indices[dimension] * WID + k + 0.5) * dvz + vz_min; //cell centered velocity
       const Real z_translation = - cell_vz * dt * i_dz; // how much it moved in time - dt (departure grid) (in units of dz)
       const Real target_scell_index = (z_translation > 0) ? 1: -1; //part of density goes here (cell index change along spatial direcion)

       for (uint j = 0; j < WID; ++j){ 
          /*compute reconstruction*/
#ifdef TRANS_SEMILAG_PLM
          Vec4 a[2];
          Vec4 mv,cv,pv;
          mv.load(values + i_pblockv(-1,j,k));
          cv.load(values + i_pblockv(0,j,k));
          pv.load(values + i_pblockv(1,j,k));
          compute_plm_coeff(mv,cv,pv,a);
#endif
#ifdef TRANS_SEMILAG_PPM
          Vec4 a[3];
          Vec4 mmv,mv,cv,pv,ppv;
          mmv.load(values + i_pblockv(-2,j,k));
          mv.load(values + i_pblockv(-1,j,k));
          cv.load(values + i_pblockv(0,j,k));
          pv.load(values + i_pblockv(1,j,k));
          ppv.load(values + i_pblockv(2,j,k));
          compute_ppm_coeff(mmv,mv,cv,pv,ppv,a);
#endif
          
	  
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
	    (v_2 - v_1) * a[0] +
	    (v_2 * v_2 - v_1 * v_1) * a[1];
#endif
#ifdef TRANS_SEMILAG_PPM
	  const Vec4 ngbr_target_density=
             (v_2 - v_1) * a[0] +
             (v_2 * v_2 - v_1 * v_1) * a[1] +
             (v_2 * v_2 * v_2 - v_1 * v_1 * v_1) * a[2];
#endif
          target_values[i_ptblockv(target_scell_index,j,k)] + =  ngbr_target_density; //in the current original cells we will put this density        
          target_values[i_ptblockv(0,j,k)] + =  cv - ngbr_target_density; //in the current original cells we will put the rest of the original density
       }
    }

    
    //store values from target_values array to the actual , one element at a time
    for(uint target_i = 0; target_i < 4;target_i ++ ){
       
       if(previous_target_block != tblock) {
          previous_target_block = tblock;
		//not the same block as last time, lets create it if we
		//need to and fetch its data array pointer and store it in target_block_data.
		if (spatial_cell->count(tblock) == 0) {
                   //count faster since the add_velocity_block call is more expensive
		  spatial_cell->add_velocity_block(tblock);
		}
		Velocity_Block* block_ptr = spatial_cell->at_fast(tblock);
		target_block_data=block_ptr->data;
	      }
	      target_block_data[tcell] += tval;
	    }
	  }
	  gk++; //next iteration in while loop
	}
      }
    }
  }
  delete [] blocks;
  return true;
}


#endif   
