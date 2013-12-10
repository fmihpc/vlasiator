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

#ifdef ACC_SEMILAG_PLM
const int STENCIL_WIDTH=1;
#else
const int STENCIL_WIDTH=2;
#endif

using namespace std;
using namespace spatial_cell;

const Real one_sixth=1.0/6.0;
const Real one_twelfth=1.0/12.0;
const Real seven_twelfth=7.0/12.0;
const Real one_third=1.0/3.0;

// indices in padded z block
#define i_pblock(i,j,k) ( ((k) + STENCIL_WIDTH ) * WID2 + (j) * WID + (i) )



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
  
  if(a<0)
    return -minval;
  else
    return minval;
}

/*!
Compute PPM coefficients as in:
 Coella 1984
 Carpenter et al. 1990
For v=-dv/2 .. dv/2 in a cell we compute
f(v)= cv + A * v + B * ( BB - v**2 )
*/

inline void compute_ppm_coeff(Real dv, Real mmv, Real mv, Real cv, Real pv, Real ppv,
                          Real &A, Real &B, Real &BB){

  //compute p_face,m_face. 
  const Real d_pv=slope_limiter(cv,pv,ppv);
  const Real d_cv=slope_limiter(mv,cv,pv);
  const Real d_mv=slope_limiter(mmv,mv,cv);
  Real p_face=0.5*(pv+cv) + one_sixth * (d_cv-d_pv);
  Real m_face=0.5*(cv+mv) + one_sixth * (d_mv-d_cv);
  
  /*
    Coella1984 eq. 1.09, claimed to give same resoult as above but that is not the case.
    Does not include the slope limiter
    
    Real p_face=seven_twelfth*(pv+cv)-one_twelfth*(ppv+mv);
    Real m_face=seven_twelfth*(cv+mv)-one_twelfth*(pv+mmv);
  */

  
  //Coella1984 eq. 1.10
  if( (p_face-cv)*(cv-m_face) <0) {
     //Extrema, cv higher/lowe than both face values. This is the
     //same as setting B=0 and A=0, so constant approximation
     p_face=cv;
     m_face=cv;
     
  }
  else if( (p_face-m_face)*(cv-0.5*(m_face+p_face))>(p_face-m_face)*(p_face-m_face)*one_sixth){
    m_face=3*cv-2*p_face;
  }
  else if( -(p_face-m_face)*(p_face-m_face)*one_sixth > (p_face-m_face)*(cv-0.5*(m_face+p_face))) {
    p_face=3*cv-2*m_face;
  }

  //Fit a second order polynomial for reconstruction
  // f(v)= cv + A * v + B * ( BB - v**2 )
  // f(-dv/2) = m_face
  // f(+dv/2) = p_face
  //see carpenter et al. 1990
  
  A=(p_face-m_face)/dv;
  B=(6*cv - 3*(p_face+m_face))/(dv*dv);    
  BB=dv*dv/12.0;
}

/*!
Piecewise cubic method as in

Zerroukat et al., SLICE: A Semi-Lagrangian Inherently Conserving and Ef cient scheme for
transport problems, Q. J. R. Meteorol. Soc. (2002), 128, pp. 2801–2820

f(v) = a[0] + a[1]*t + a[2]*t**2 + a[3]*t**3
t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
*/

inline void compute_pcm_coeff(Real mmv, Real mv, Real cv, Real pv, Real ppv,
			      Real * __restrict__ a){
  //compute rho_R,rho_L exactly as in PPM. 
  /*
  const Real d_pv=slope_limiter(cv,pv,ppv);
  const Real d_cv=slope_limiter(mv,cv,pv);
  const Real d_mv=slope_limiter(mmv,mv,cv);
  Real rho_R=0.5*(pv+cv) + one_sixth * (d_cv-d_pv);
  Real rho_L=0.5*(cv+mv) + one_sixth * (d_mv-d_cv);
  */
  /*
    Coella1984 eq. 1.09, claimed to give same resoult as above but that is not the case.
    Does not include the slope limiter
  */  
  Real rho_R=seven_twelfth*(pv+cv)-one_twelfth*(ppv+mv);
  Real rho_L=seven_twelfth*(cv+mv)-one_twelfth*(pv+mmv);
  //derivative, eq (17) in Zerroukat 2002. See notebook sheet for computation details
  // Real d_rho = (-12.0*cv - 22.0*mv + 5.0*mmv + 22.0*pv + 7.0*ppv)/48.0;
  Real d_rho=seven_twelfth*(pv-mv) + 1.0/24.0 * ( mmv - ppv );
  //Coella1984 eq. 1.10
  if( (rho_R-cv)*(cv-rho_L) <0) {
    //Maxima, constant 
    rho_R=cv;
    rho_L=cv;
    d_rho=0;
  }

  /*
    else if( (rho_R-rho_L)*(cv-0.5*(rho_L+rho_R))>(rho_R-rho_L)*(rho_R-rho_L)*one_sixth){
    rho_L=3*cv-2*rho_R;
    }
    else if( -(rho_R-rho_L)*(rho_R-rho_L)*one_sixth > (rho_R-rho_L)*(cv-0.5*(rho_L+rho_R))) {
    rho_R=3*cv-2*rho_L;
    }
  */

  //Coefficients as in eq. (13) 
  a[0] = rho_L;
  a[1] = 6.0*(cv - rho_L)                     - 2.0 * d_rho;
  a[2] = 3.0*(3.0 * rho_L - 2.0 * cv - rho_R) + 6.0 * d_rho;
  a[3] = 4.0*(rho_R - rho_L)                  - 4.0 * d_rho;
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
    for (uint k=0; k<WID; ++k){ 
      //todo, could also compute z index and compute the start velocity of block
      const Real v_l=(WID*block_indices[2]+k)*dv+v_min;
      const Real v_c=v_l+0.5*dv;
      const Real v_r=v_l+dv;
      for (uint j = 0; j < WID; ++j){
         for (uint i = 0; i < WID; ++i){
	  //f is mean value
	  const Real cv = values[i_pblock(i,j,k)];
          
#ifdef ACC_SEMILAG_PLM
	  //A is slope of linear approximation
	  const Real A = slope_limiter(values[i_pblock(i,j,k-1)],
				       cv,values[i_pblock(i,j,k+1)])*i_dv;
#endif
#ifdef ACC_SEMILAG_PPM
          Real A,B,BB;
          compute_ppm_coeff(dv,values[i_pblock(i,j,k-2)],values[i_pblock(i,j,k-1)],
			    cv,values[i_pblock(i,j,k+1)],values[i_pblock(i,j,k+2)],
			    A,B,BB);
#endif
#ifdef ACC_SEMILAG_PCM
	  Real a[4];
          compute_pcm_coeff(values[i_pblock(i,j,k-2)],values[i_pblock(i,j,k-1)],
			    cv,values[i_pblock(i,j,k+1)],values[i_pblock(i,j,k+2)],a);
#endif	  

	  //left(l) and right(r) k values (global index) in the target
	  //lagrangian grid, the intersecting cells
	  const Real intersection_min=intersection +
	    (block_indices[0]*WID+i)*intersection_di + 
             (block_indices[1]*WID+j)*intersection_dj;
	  const uint lagrangian_gk_l=(v_l-intersection_min)/intersection_dk; 
	  const uint lagrangian_gk_r=(v_r-intersection_min)/intersection_dk;
	  //todo, annoying loop. probably pretty bad for vectorization
	  //since each i,j pair can have a different amount
	  //here... based on distances one knows beforhand if the loop
	  //can have 1-2, or 2-3 iterations, we could split the whole
	  //thing..
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



#ifdef ACC_SEMILAG_PLM	    
	    //Note that we also have shifted velocity to have origo at v_c
	    //(center velocity of euclidian cell)
	    const Real v_1 = max(gk * intersection_dk + intersection_min, v_l)-v_c;
	    const Real v_2 = min((gk+1) * intersection_dk + intersection_min, v_r)-v_c;
	    //target mass is value in center of intersecting length,
	    //times length (missing x,y, but they would be cancelled
	    //anyway when we divide to get density
	    const Real target_density = (cv + A * 0.5*(v_1+v_2))*(v_2 - v_1)*i_dv;
#endif
#ifdef ACC_SEMILAG_PPM
	    //Note that we also have shifted velocity to have origo at v_c
	    //(center velocity of euclidian cell)
	    const Real v_1 = max(gk * intersection_dk + intersection_min, v_l)-v_c;
	    const Real v_2 = min((gk+1) * intersection_dk + intersection_min, v_r)-v_c;
	    //todo, we recompute integrals, could do some reuse at least over gk loop (same with v also)
            const Real integral_1=v_1*(cv+0.5*A*v_1+B*(BB-v_1*v_1*one_third));
            const Real integral_2=v_2*(cv+0.5*A*v_2+B*(BB-v_2*v_2*one_third));
            const Real target_density=(integral_2-integral_1)*i_dv;
#endif
#ifdef ACC_SEMILAG_PCM
	    //in pcm we normalize by dv and in the center cell the v coordinate goes from 0 to 1
	    const Real v_1 = (max(gk * intersection_dk + intersection_min, v_l)-v_l)/dv;
	    const Real v_2 = (min((gk+1) * intersection_dk + intersection_min, v_r)-v_l)/dv;
	    //todo, we recompute integrals, could do some reuse at least over gk loop (same with v also)
            const Real integral_1=v_1*(a[0]+0.5*a[1]*v_1+one_third*a[2]*v_1*v_1+0.25*a[3]*v_1*v_1*v_1);
            const Real integral_2=v_2*(a[0]+0.5*a[1]*v_2+one_third*a[2]*v_2*v_2+0.25*a[3]*v_2*v_2*v_2);
            const Real target_density=integral_2-integral_1;

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
