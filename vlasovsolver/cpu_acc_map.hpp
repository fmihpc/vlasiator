/*
This file is part of Vlasiator.

Copyright 2012 Finnish Meteorological Institute

*/

#ifndef CPU_ACC_MAP_H
#define CPU_ACC_MAP_H

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

const Real one_sixth=1.0/6.0;
const Real one_twelfth=1.0/12.0;
const Real seven_twelfth=7.0/12.0;
const Real one_third=1.0/3.0;

// indices in padded z block
//#define i_pblock(i,j,k) ( ((k) + STENCIL_WIDTH ) * WID2 + (j) * WID + (i) ) // x first, then y and z
#define i_pblock(i,j,k) ( ((k) + STENCIL_WIDTH ) * WID + (j) * WID * (WID + 2* STENCIL_WIDTH) + (i) )


/*!
  MC slope limiter
*/

template<typename T> inline T slope_limiter(const T& l,const T& m, const T& r) {
  T sign;
  T a=r-m;
  T b=m-l; 
  if (a*b < 0.0)
    return 0.0;
  T minval=min(2.0*fabs(a),2.0*fabs(b));
  minval=min(minval,0.5*fabs(a+b));
  
  if(a+b<0)
    return -minval;
  else
    return minval;
}


/*!
 Compute PLM coefficients
f(v) = a[0] + a[1]/2.0*t 
t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
The factor 2.0 is in the polynom to ease integration, then integral is a[0]*t + a[1]*t**2
*/

inline void compute_plm_coeff(Real mv,Real cv, Real pv,
			      Real * __restrict__ a){
  const Real d_cv=slope_limiter(mv,cv,pv);
  a[0] = cv - d_cv * 0.5;
  a[1] = d_cv * 0.5;
}




inline void compute_ppm_coeff(Real mmv,Real mv, Real cv, Real pv,Real ppv,
			      Real * __restrict__ a){
   Real p_face;
   Real m_face;

   //Compute estimations of the face values. 

#ifdef PPM_FACEEST_H4_POS
   //h4 estimates (fourth order), white 2008. No negative values, no other filtering here
   p_face=seven_twelfth*(pv+cv)-one_twelfth*(ppv+mv);
   m_face=seven_twelfth*(cv+mv)-one_twelfth*(pv+mmv);
   if(p_face<0) p_face=0;
   if(m_face<0) m_face=0;

#endif

#ifdef PPM_FACEEST_LIMITED
   //limited estimates, like in Coella 1984. This is (in my opinion), the best one. 
   const Real d_cv=slope_limiter(mv,cv,pv);
   const Real d_pv=slope_limiter(cv,pv,ppv);
   const Real d_mv=slope_limiter(mmv,mv,cv);

   p_face=0.5*(pv+cv) + one_sixth * (d_cv-d_pv);
   m_face=0.5*(cv+mv) + one_sixth * (d_mv-d_cv);        
#endif



  //Coella1984 eq. 1.10, detect extream
   if( (p_face-cv)*(cv-m_face) <0) {
     //set to constant approximation
     p_face=cv;
     m_face=cv;
   }

   //Coella et al, check for monotonicity
   else if( (p_face-m_face)*(cv-0.5*(m_face+p_face))>(p_face-m_face)*(p_face-m_face)*one_sixth){
     m_face=3*cv-2*p_face;
   }
   else if( -(p_face-m_face)*(p_face-m_face)*one_sixth > (p_face-m_face)*(cv-0.5*(m_face+p_face))) {
     p_face=3*cv-2*m_face;
   }

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
  uint block_indices_to_id[3];
  uint cell_indices_to_id[3];
  Real integration_v[4]; //v boundaries for integration when mapping
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
     */
    //TODO, the data layout is not optimal anymore, in memory we shoul dhave i first, k then and then j
    
    for (uint j = 0; j < WID; ++j){ 
      for (uint i = 0; i < WID; ++i){ //TODO, use i to vectorize computation
	const Real intersection_min=intersection +
	  (block_indices[0]*WID+i)*intersection_di + 
	  (block_indices[1]*WID+j)*intersection_dj;
	for (uint k=0; k<WID; ++k){ 
	  const Real v_l=(WID*block_indices[2]+k)*dv+v_min;
	  const Real v_c=v_l+0.5*dv;
	  const Real v_r=v_l+dv;
	  //left(l) and right(r) k values (global index) in the target
	  //lagrangian grid, the intersecting cells
	  const uint lagrangian_gk_l=(v_l-intersection_min)/intersection_dk; 
	  const uint lagrangian_gk_r=(v_r-intersection_min)/intersection_dk;

#ifdef ACC_SEMILAG_PLM
          Real a[2];
          compute_plm_coeff(values[i_pblock(i,j,k-1)],
			    values[i_pblock(i,j,k  )],
			    values[i_pblock(i,j,k+1)],
			    a);
#endif
#ifdef ACC_SEMILAG_PPM
          Real a[3];
	  compute_ppm_coeff(values[i_pblock(i,j,k-2)],
			    values[i_pblock(i,j,k-1)],
			    values[i_pblock(i,j,k  )],
			    values[i_pblock(i,j,k+1)],
			    values[i_pblock(i,j,k+2)],
			    a);
#endif
	  //add values to target cell
          for(uint gk=lagrangian_gk_l;gk<=lagrangian_gk_r;gk++){
             //the blocks of the lagrangian cell to which we map
             const uint target_block = 
	      block_indices[0]*block_indices_to_id[0]+
	      block_indices[1]*block_indices_to_id[1]+
	      (gk/WID)*block_indices_to_id[2];
	    //cell index in the target block 
	    const uint target_cell = 
               i*cell_indices_to_id[0] + 
	       j*cell_indices_to_id[1] +
	       (gk%WID)*cell_indices_to_id[2];
	    //the velocity between which we will integrate to put mass
	    //in the targe cell. If both v_r and v_l are in same cell
	    //then v_1,v_2 should be between v_l and v_r.
            //v_1 and v_2 normalized to be between 0 and 1 in the cell
	    const Real v_1 = (max(gk * intersection_dk + intersection_min, v_l)-v_l)*i_dv;
	    const Real v_2 = (min((gk+1) * intersection_dk + intersection_min, v_r)-v_l)*i_dv;
#ifdef ACC_SEMILAG_PCONSTM	    
            const Real target_density=
	      (v_2 - v_1) *  values[i_pblock(i,j,k)];
#endif
#ifdef ACC_SEMILAG_PLM	    
            const Real target_density=
	      (v_2 - v_1) * a[0] +
	      (v_2 * v_2 - v_1 * v_1) * a[1];
#endif
#ifdef ACC_SEMILAG_PPM
            const Real target_density=
               (v_2 - v_1) * a[0] +
               (v_2 * v_2 - v_1 * v_1) * a[1] +
               (v_2 * v_2 * v_2 - v_1 * v_1 * v_1) * a[2];
#endif
	    if (target_block < SpatialCell::max_velocity_blocks) 
	      spatial_cell->increment_value(target_block,target_cell,target_density);
	  }
	}
      }
    }
  }
  delete [] blocks;
  return true;
}


#endif   
