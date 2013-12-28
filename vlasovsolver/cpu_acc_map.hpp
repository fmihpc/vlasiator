
/*
This file is part of Vlasiator.
Copyright 2012 Finnish Meteorological Institute

*/

#ifndef CPU_ACC_MAP_H
#define CPU_ACC_MAP_H

#include  "vec4.h"
#include "algorithm"
#include "cmath"
#include "utility"
#include "common.h"
#include "spatial_cell.hpp"

#ifdef ACC_SEMILAG_PCONSTM
const int STENCIL_WIDTH=0;
#endif
#ifdef ACC_SEMILAG_PLM
const int STENCIL_WIDTH=1;
#endif
#if ACC_SEMILAG_PPM
const int STENCIL_WIDTH=2;
#endif



using namespace std;
using namespace spatial_cell;

const Vec4 one_sixth(1.0/6.0);
const Vec4 one_twelfth(1.0/12.0);
const Vec4 seven_twelfth(7.0/12.0);
const Vec4 one_third(1.0/3.0);

// indices in padded z block
//#define i_pblock(i,j,k) ( ((k) + STENCIL_WIDTH ) * WID2 + (j) * WID + (i) ) // x first, then y and z
#define i_pblock(i,j,k) ( ((k) + STENCIL_WIDTH ) * WID + (j) * WID * (WID + 2* STENCIL_WIDTH) + (i) )
#define i_pblockv(j,k) ( ((k) + STENCIL_WIDTH ) * WID + (j) * WID * (WID + 2* STENCIL_WIDTH) )


/*!
  MC slope limiter
*/

inline Vec4 slope_limiter(const Vec4& l,const Vec4& m, const Vec4& r) {
  const Vec4 two(2.0);
  const Vec4 half(0.5);
  const Vec4 zero(0.0);
  Vec4 sign;
  Vec4 a=r-m;
  Vec4 b=m-l; 
  Vec4 minval=min(two*abs(a),two*abs(b));
  minval=min(minval,half*abs(a+b));
  
  //check for extrema
  Vec4 output = select(a*b < 0,zero,minval);
  
  //set sign
  return select(a + b < 0,-output,output);
}


/*!
 Compute PLM coefficients
f(v) = a[0] + a[1]/2.0*t 
t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
The factor 2.0 is in the polynom to ease integration, then integral is a[0]*t + a[1]*t**2
*/

inline void compute_plm_coeff(Vec4 mv,Vec4 cv, Vec4 pv,
			      Vec4 * __restrict__ a){
  const Vec4 d_cv=slope_limiter(mv,cv,pv);
  a[0] = cv - d_cv * 0.5;
  a[1] = d_cv * 0.5;
}




inline void compute_ppm_coeff(Vec4 mmv,Vec4 mv, Vec4 cv, Vec4 pv,Vec4 ppv,
			      Vec4 * __restrict__ a){
   Vec4 p_face;
   Vec4 m_face;

   //Compute estimations of the face values. Limited estimates, like in Coella 1984. 
   const Vec4 d_cv=slope_limiter(mv,cv,pv);
   const Vec4 d_pv=slope_limiter(cv,pv,ppv);
   const Vec4 d_mv=slope_limiter(mmv,mv,cv);

   p_face=0.5*(pv+cv) + one_sixth * (d_cv-d_pv);
   m_face=0.5*(cv+mv) + one_sixth * (d_mv-d_cv);        

   //Coella1984 eq. 1.10, detect extrema
   m_face = select((p_face - cv) * (cv - m_face) < 0, m_face, cv);
   p_face = select((p_face - cv) * (cv - m_face) < 0, p_face, cv);
   
   //Coella et al, check for monotonicity   
   m_face = select((p_face-m_face)*(cv-0.5*(m_face+p_face))>(p_face-m_face)*(p_face-m_face)*one_sixth,
		  3*cv-2*p_face,
		  m_face);
   p_face = select(-(p_face-m_face)*(p_face-m_face)*one_sixth > (p_face-m_face)*(cv-0.5*(m_face+p_face)),
		   3*cv-2*m_face,
		   p_face);

   //Fit a second order polynomial for reconstruction see, e.g., White
   //2008 (PQM article) (note additional integration factors built in,
   //contrary to White (2008) eq. 4
   a[0]=m_face;
   a[1]=3.0*cv-2.0*m_face-p_face;
   a[2]=(m_face+p_face-2.0*cv);
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
inline void copy_block_data(Velocity_Block *block,Real * __restrict__ values, int dimension){
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
   tracked backwards by -dt)
*/

bool map_1d(SpatialCell* spatial_cell,   
	    Real intersection, Real intersection_di, Real intersection_dj,Real intersection_dk,
	    uint dimension ) {

  /*values used with an stencil in 1 dimension, initialized to 0 */  
  Real values[WID*WID*(WID+2*STENCIL_WIDTH)]={};
  // Make a copy of the blocklist, the blocklist will change during this algorithm
  uint*  blocks=new uint[spatial_cell->number_of_blocks];
  const uint nblocks=spatial_cell->number_of_blocks;
  /* 
     Move densities from data to fx and clear data, to prepare for mapping
     Also copy blocklist since the velocity block list in spatial cell changes when we add values
  */
  for (unsigned int block_i = 0; block_i < nblocks; block_i++) {
    const unsigned int block = spatial_cell->velocity_block_list[block_i];
    blocks[block_i] = block; 
    Velocity_Block * __restrict__ block_ptr = spatial_cell->at(block);
    Real * __restrict__ fx = block_ptr->fx;
    Real * __restrict__ data = block_ptr->data;

    for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
      fx[cell] = data[cell];
      data[cell] = 0.0;
    }
  }
  if(dimension>2)
    return false; //not possible

  
  Real dv,v_min;
  Real is_temp;
  uint block_indices_to_id[3]; /*< used when computing id of target block */
  uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
  switch (dimension){
  case 0:
    /* i and k coordinates have been swapped*/
    /*set cell size in dimension direction*/
    dv=SpatialCell::cell_dvx; 
    v_min=SpatialCell::vx_min; 
    /*swap intersection i and k coordinates*/
    is_temp=intersection_di;
    intersection_di=intersection_dk;
    intersection_dk=is_temp;
    /*set values in array that is used to transfer blockindices to id using a dot product*/
    block_indices_to_id[0]=SpatialCell::vx_length * SpatialCell::vy_length;
    block_indices_to_id[1]=SpatialCell::vx_length;
    block_indices_to_id[2]=1;
    /*set values in array that is used to transfer blockindices to id using a dot product*/
    cell_indices_to_id[0]=WID2;
    cell_indices_to_id[1]=WID;
    cell_indices_to_id[2]=1;
    break;
  case 1:
    /* j and k coordinates have been swapped*/
    /*set cell size in dimension direction*/
    dv=SpatialCell::cell_dvy;
    v_min=SpatialCell::vy_min; 
    /*swap intersection j and k coordinates*/
    is_temp=intersection_dj;
    intersection_dj=intersection_dk;
    intersection_dk=is_temp;
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
    /*set cell size in dimension direction*/
    dv=SpatialCell::cell_dvz;
    v_min=SpatialCell::vz_min; 
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
  const Real i_dv=1.0/dv;

  
  for (unsigned int block_i = 0; block_i < nblocks; block_i++) {
    Velocity_Block *block=spatial_cell->at(blocks[block_i]);
    velocity_block_indices_t block_indices=SpatialCell::get_velocity_block_indices(blocks[block_i]);
    uint temp;
    //Switch block indices according to dimensions, the alogirthm has
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
    copy_block_data(block,values,dimension); 

    /*i,j,k are now relative to the order in which we copied data to the values array. 
      After this point in the k,j,i loops there should be no branches based on dimensions
      
      Note that the i dimension is vectorized, and thus there are no loops over i
    */
    for (uint j = 0; j < WID; ++j){ 
      const Real intersection_min_base = intersection +
	(block_indices[0]*WID)*intersection_di + 
	(block_indices[1]*WID+j)*intersection_dj;
      /* 
	 intersection_min is the intersection z coordinate (z after
	 swaps that is) of the lowest possible z plane for each i,j
	 index (i in vector)
      */	 
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
	mmv.load(values + i_pblockv(j,k-1));
	mv.load(values + i_pblockv(j,k-1));
	cv.load(values + i_pblockv(j,k));
	pv.load(values + i_pblockv(j,k+1));
	ppv.load(values + i_pblockv(j,k+2));
	compute_ppm_coeff(mmv,mv,cv,pv,ppv,a);
#endif
	Vec4i gk = lagrangian_gk_l;	
	while (horizontal_or(gk <= lagrangian_gk_r)){
	  const Vec4i gk_div_WID = gk/WID;
	  const Vec4i gk_mod_WID = (gk - gk_div_WID * WID);
	  //the block of the lagrangian cell to which we map
	  const Vec4i target_block = 	
	    block_indices[0]*block_indices_to_id[0]+
	    block_indices[1]*block_indices_to_id[1]+
	    gk_div_WID*block_indices_to_id[2];
	  //cell index in the target block 
	  const Vec4i target_cell_base = j*cell_indices_to_id[1] + gk_mod_WID*cell_indices_to_id[2];
	  const Vec4i target_cell(target_cell_base + Vec4i(0, cell_indices_to_id[0], 2 * cell_indices_to_id[0], 3 * cell_indices_to_id[0]));
	  
	  //the velocity between which we will integrate to put mass
	  //in the targe cell. If both v_r and v_l are in same cell
	  //then v_1,v_2 should be between v_l and v_r.
	  //v_1 and v_2 normalized to be between 0 and 1 in the cell.
	  //For vector elements where gk is already larger than needed (lagrangian_gk_r), v_2=v_1=v_r and thus the value is zero.
	  const Vec4 v_1 = (min(max(to_double(gk) * intersection_dk + intersection_min, v_l), v_r) - v_l) * i_dv;
	  const Vec4 v_2 = (min(to_double(gk + 1) * intersection_dk + intersection_min, v_r) - v_l) * i_dv;
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
	  for(uint k = 0; k < WID;k ++ ){
	    if (target_block[k] < SpatialCell::max_velocity_blocks) 
	      spatial_cell->increment_value(target_block[k],target_cell[k],target_density[k]);
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
