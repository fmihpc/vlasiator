
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
const int STENCIL_WIDTH=0;
#endif
#ifdef TRANS_SEMILAG_PLM
const int STENCIL_WIDTH=1;
#endif
#if TRANS_SEMILAG_PPM
const int STENCIL_WIDTH=2;
#endif



using namespace std;
using namespace spatial_cell;

// indices in padded block. b_k is the block index in z direction in ordinary space, i,j,k are the cell ids inside on block. Data order i,b_k,j,k.
// b_k range is  -STENCIL_WIDTH -  + STENCil_WIDTH, so in total 1 + 2 * STENCIL_WIDTH values 
#define i_pblock(b_k,i,j,k) ( (i) + ((b_k) + STENCIL_WIDTH) * WID + (j) * WID * (1 + 2 * STENCIL_WIDTH) + (k) * WID2 * (1 + 2 * STENCIL_WIDTH))
#define i_pblockv(b_k,j,k) ( ((b_k) + STENCIL_WIDTH) * WID + (j) * WID * (1 + 2 * STENCIL_WIDTH) + (k) * WID2 * (1 + 2 * STENCIL_WIDTH))



CellID getNeighbourID(
   const dccrg::Dccrg<SpatialCell>& mpiGrid,
   const CellID& cellID,
   const uchar& i,
   const uchar& j,
   const uchar& k ) {
  const std::vector<CellID> neighbors = mpiGrid.get_neighbors_of(cellID, int(i) - 2, int(j) - 2, int(k) - 2);
  if (neighbors.size() == 0) {
    std::cerr << __FILE__ << ":" << __LINE__
	      << " No neighbor for cell " << cellID
	      << " at offsets " << int(i) - 2 << ", " << int(j) - 2 << ", " << int(k) - 2
                 << std::endl;
    abort();
  }
  // TODO support spatial refinement
  if( neighbors[0] == INVALID_CELLID  ||
      mpiGrid[neighbors[0]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
      (mpiGrid[neighbors[0]]->sysBoundaryLayer != 1  &&
       mpiGrid[neighbors[0]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY )
      ) {
     return INVALID_CELLID;
  } 
  else {
    return neighbors[0];
  }
}

/*!
  value array should be initialized to zero
  
  For dimension=0 data copy  we have rotated data
  i -> k
  j -> j
  k -> i
For dimension=1 data copy  we have rotated data
  i -> i
  j -> k
  k -> j
For dimension=0 data copy 
*/
inline void copy_trans_block_data(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CellID cellID,uint dimension ) {

  Velocity_Block *block=spatial_cell->at(blockID);
  Velocity_Block *nbrBlock;
    uint cell_indices_to_id[3];
    uint block_P1,block_M1;
    switch (dimension){
        case 0:
           /* i and k coordinates have been swapped*/
           cell_indices_to_id[0]=WID2;
           cell_indices_to_id[1]=WID;
           cell_indices_to_id[2]=1;
           block_P1=velocity_neighbor::XP1_YCC_ZCC;
           block_M1=velocity_neighbor::XM1_YCC_ZCC;
           break;
        case 1:
           /* i and k coordinates have been swapped*/
           cell_indices_to_id[0]=1;
           cell_indices_to_id[1]=WID2;
           cell_indices_to_id[2]=WID;
           block_P1=velocity_neighbor::XCC_YP1_ZCC;
           block_M1=velocity_neighbor::XCC_YM1_ZCC;
           break;
        case 2:
           cell_indices_to_id[0]=1;
           cell_indices_to_id[1]=WID;
           cell_indices_to_id[2]=WID2;
           block_P1=velocity_neighbor::XCC_YCC_ZP1;
           block_M1=velocity_neighbor::XCC_YCC_ZM1;
           break;
    }
    
    // Construct values
    // Copy averages from -1 neighbour if it exists (if not, array is initialized to zero)
    nbrBlock = block->neighbors[block_M1];
    if ( nbrBlock != NULL) {
       Real * __restrict__ ngbr_fx = nbrBlock->fx;
       for (int k=-STENCIL_WIDTH; k<0; ++k) {
          for (uint j=0; j<WID; ++j) {
             for (uint i=0; i<WID; ++i) {
                const uint cell =
                   i * cell_indices_to_id[0] +
                   j * cell_indices_to_id[1] +
                   (k + WID) * cell_indices_to_id[2];
                values[i_pblock(i,j,k)] = ngbr_fx[cell];
             }
          }
       }
    }
    else {
       for (int k=-STENCIL_WIDTH; k<0; ++k) {
          for (uint j=0; j<WID; ++j) {
             for (uint i=0; i<WID; ++i) {
	       values[i_pblock(i,j,k)] = 0.0;
             }
          }
       }
    }


    Real * __restrict__ fx = block->fx;
    // Copy volume averages of this block:
    for (uint k=0; k<WID; ++k) {
       for (uint j=0; j<WID; ++j) {
          for (uint i=0; i<WID; ++i) {
             const uint cell =
                i * cell_indices_to_id[0] +
                j * cell_indices_to_id[1] +
                k * cell_indices_to_id[2];
             values[i_pblock(i,j,k)] = fx[cell];
          }
       }
    }
    
    // Copy averages from +1 neighbour if it exists (if not, array is initialized to zero)
    nbrBlock = block->neighbors[block_P1];
    if ( nbrBlock != NULL) {
       Real * __restrict__ ngbr_fx = nbrBlock->fx;
       for (uint k=WID; k<WID+STENCIL_WIDTH; ++k) {             
          for (uint j=0; j<WID; ++j) {
             for (uint i=0; i<WID; ++i) {
                const uint cell =
                   i * cell_indices_to_id[0] +
                   j * cell_indices_to_id[1] +
                   (k-WID) * cell_indices_to_id[2];
                values[i_pblock(i,j,k)] = ngbr_fx[cell];
             }
          }
       }
    }
    else {
      //no neighbor, set neighbor values to zero
      for (uint k=WID; k<WID+STENCIL_WIDTH; ++k) {             
	for (uint j=0; j<WID; ++j) {
	  for (uint i=0; i<WID; ++i) {
	    values[i_pblock(i,j,k)] = 0.0;
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

bool map_trans_1d(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CellID cellID,uint dimension ) {
  /*values used with an stencil in 1 dimension, initialized to 0. Contains a block, and its spatial neighbours in one dimension */  
  Real values[WID3*(1+2*STENCIL_WIDTH)];
  Real dz,z_min;
  SpatialCell* spatial_cell = mpiGrid[cellID];

  if(dimension>2)
    return false; //not possible

  /*set cell size in dimension direction*/  
  switch (dimension){
  case 0:
    dz = P::dx_ini;
    z_min = P::xmin;
    break;
  case 1:
    dz = P::dy_ini;
    z_min = P::ymin;
    break;
  case 2:
    dz = P::dz_ini;
    z_min = P::zmin;
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
    const unsigned int block = spatial_cell->velocity_block_list[block_i];
    Velocity_Block * __restrict__ block_ptr = spatial_cell->at(block);
    
    copy_block_data(block_ptr,values,dimension); 
    velocity_block_indices_t block_indices=SpatialCell::get_velocity_block_indices(block);
    
    
    uint temp;
    //Switch block indices according to dimensions, the alogrithm has
    //been written for integrating along z.
    switch (dimension){
        case 0:
	  /*i and k coordinates have been swapped*/
	  temp=block_indices[2];
	  block_indices[2]=block_indices[0];
	  block_indices[0]=temp;
	  break;
    case 1:
           /*in values j and k coordinates have been swapped*/
           temp=block_indices[2];
           block_indices[2]=block_indices[1];
           block_indices[1]=temp;
           break;
        case 2:
           break;
    }


    /*i,j,k are now relative to the order in which we copied data to the values array. 
      After this point in the k,j,i loops there should be no branches based on dimensions
      
      Note that the i dimension is vectorized, and thus there are no loops over i
    */
    for (uint j = 0; j < WID; ++j){ 
      /*target cell/block index contribution not dependent on k index*/
      const Vec4i target_cell_index_common = j*cell_indices_to_id[1] + Vec4i(0, cell_indices_to_id[0], 2 * cell_indices_to_id[0], 3 * cell_indices_to_id[0]);
      const int target_block_index_common(block_indices[0]*block_indices_to_id[0]+block_indices[1]*block_indices_to_id[1]);
      /* 
	 intersection_min is the intersection z coordinate (z after
	 swaps that is) of the lowest possible z plane for each i,j
	 index (i in vector)
      */	 
      const Real intersection_min_base = intersection +
	(block_indices[0]*WID)*intersection_di + 
	(block_indices[1]*WID+j)*intersection_dj;
      const Vec4 intersection_min(intersection_min_base,
				  intersection_min_base + intersection_di,
				  intersection_min_base + 2.0 * intersection_di,
				  intersection_min_base + 3.0 * intersection_di);
      
      
      for (uint k=0; k<WID; ++k){ 
	/*v_l, v_c, v_r are the left, mid and right velocity coordinates of source cell*/
	Vec4 v_l((WID * block_indices[2] + k) * dv + v_min);
	Vec4 v_c(v_l+0.5*dv);
	Vec4 v_r(v_l+dv);
	//left(l) and right(r) k values (global index) in the target
	//lagrangian grid, the intersecting cells
	const Vec4i lagrangian_gk_l=truncate_to_int((v_l-intersection_min)/intersection_dk);
	const Vec4i lagrangian_gk_r=truncate_to_int((v_r-intersection_min)/intersection_dk);
	
#ifdef ACC_SEMILAG_PLM
	Vec4 a[2];
	Vec4 mv,cv,pv;
	mv.load(values + i_pblockv(j,k-1));
	cv.load(values + i_pblockv(j,k));
	pv.load(values + i_pblockv(j,k+1));
	compute_plm_coeff(mv,cv,pv,a);

#endif
#ifdef ACC_SEMILAG_PPM
	Vec4 a[3];
	Vec4 mmv,mv,cv,pv,ppv;
	mmv.load(values + i_pblockv(j,k-2));
	mv.load(values + i_pblockv(j,k-1));
	cv.load(values + i_pblockv(j,k));
	pv.load(values + i_pblockv(j,k+1));
	ppv.load(values + i_pblockv(j,k+2));
	compute_ppm_coeff(mmv,mv,cv,pv,ppv,a);
#endif
	Vec4i gk(lagrangian_gk_l);	
	while (horizontal_or(gk <= lagrangian_gk_r)){
	  const Vec4i gk_div_WID = gk/WID;
	  const Vec4i gk_mod_WID = (gk - gk_div_WID * WID);
	  //the block of the lagrangian cell to which we map
	  const Vec4i target_block(target_block_index_common + gk_div_WID * block_indices_to_id[2]);
	  //cell index in the target block 
	  const Vec4i target_cell(target_cell_index_common + gk_mod_WID * cell_indices_to_id[2]);
	  
	  //the velocity between which we will integrate to put mass
	  //in the targe cell. If both v_r and v_l are in same cell
	  //then v_1,v_2 should be between v_l and v_r.
	  //v_1 and v_2 normalized to be between 0 and 1 in the cell.
	  //For vector elements where gk is already larger than needed (lagrangian_gk_r), v_2=v_1=v_r and thus the value is zero.
	  const Vec4 v_1 = (min(max(to_double(gk) * intersection_dk + intersection_min, v_l), v_r) - v_l) * i_dv;
	  const Vec4 v_2 = (min(to_double(gk + 1) * intersection_dk + intersection_min,       v_r) - v_l) * i_dv;
#ifdef ACC_SEMILAG_PCONSTM
	  Vec4 cv;	    
	  cv.load(values + i_pblockv(j,k));
	  const Vec4 target_density=(v_2 - v_1) *  cv;
#endif
#ifdef ACC_SEMILAG_PLM	    
	  const Vec4 target_density=
	    (v_2 - v_1) * a[0] +
	    (v_2 * v_2 - v_1 * v_1) * a[1];
#endif
#ifdef ACC_SEMILAG_PPM
	  const Vec4 target_density=
	    (v_2 - v_1) * a[0] +
	    (v_2 * v_2 - v_1 * v_1) * a[1] +
	    (v_2 * v_2 * v_2 - v_1 * v_1 * v_1) * a[2];
#endif
	  
	  //store values, one element at a time
	  for(uint target_i = 0; target_i < 4;target_i ++ ){
	    const uint tblock=target_block[target_i];
	    const uint tcell=target_cell[target_i];
	    const Real tval=target_density[target_i];
	    if (tblock < SpatialCell::max_velocity_blocks && tval != 0.0) {	      
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
