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

#ifdef ACC_SEMILAG_ORDER_1
const int STENCIL_WIDTH=1;
#else
const int STENCIL_WIDTH=2;
#endif

using namespace std;
using namespace spatial_cell;





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


  


inline void fit_poly2_rec(Real dv, Real mmv, Real mv, Real cv, Real pv, Real ppv,
                          Real &A, Real &B, Real &BB){
   // First fit a second order         polynomial for faces(it would acually be better with higher order?)
    // f(x)=a+b*x+c*x*x
    // f(0) = cv
    // f(-1) = mv 
    // f(+1) = pv
    // m_face=f(-0.5)
    // p_face=f(0.5)
    Real a=cv;
    Real b=(pv - mv)*0.5;
    Real c=(pv + mv)*0.5 - cv;
    Real m_face=a-b*0.5+c*0.25; 
    Real p_face=a+b*0.5+c*0.25;

    //face values should not cause new extremas, use linear average if it is above, below either cell value
    Real m_face_max=max(mv,cv);
    Real m_face_min=min(mv,cv);
    if(m_face < m_face_min || m_face > m_face_max)
      m_face=(mv+cv)*0.5; 
    Real p_face_max=max(pv,cv);
    Real p_face_min=min(pv,cv);
    if(p_face < p_face_min || p_face > p_face_max)
       p_face=(pv+cv)*0.5; 
    
    //Fit a second order polynomial for reconstruction, set rec_par[A], rec_par[B]
    // f(v)= cv + A * v + B * ( BB - v**2 )
    // f(-dv/2) = m_face
    // f(+dv/2) = p_face
    //see carpenter et al. 1990
    A=(p_face-m_face)/dv;
    B=(6*cv - 3*(p_face+m_face))/(dv*dv);    
    BB=dv*dv/12.0;
    
    //Hack to see if we are monotonic in cell, if not compute a linear approximation
    Real extrema_v=A/(2*B);
    if(extrema_v>-0.5*dv && extrema_v < 0.5*dv) {
       A=slope_limiter(mv,cv,pv)/dv;
       B=0.0;
    }

  }



// indices in padded z block
#define i_pblock(i,j,k) ( ((k) + STENCIL_WIDTH ) * WID2 + (j) * WID + (i) )
// indices in normal block
#define i_block(i,j,k) ( (k) * WID2 + (j) * WID + (i) )

/*
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


bool map_1d(SpatialCell* spatial_cell,   Real intersection, Real intersection_di, Real intersection_dj,Real intersection_dk,
	   uint dimension ) {
   const Real onethird=1.0/3.0;
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
          
#ifdef ACC_SEMILAG_ORDER_1
	  //A is slope of linear approximation
	  const Real A = slope_limiter(values[i_pblock(i,j,k-1)],
				       cv,values[i_pblock(i,j,k+1)])*i_dv;
#endif
#ifdef ACC_SEMILAG_ORDER_2
          Real A,B,BB;
          fit_poly2_rec(dv,values[i_pblock(i,j,k-2)],values[i_pblock(i,j,k-1)],
                        cv,values[i_pblock(i,j,k+1)],values[i_pblock(i,j,k+2)],
                        A,B,BB);
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
	    
	    //the velocity between which we will integrate to put mass int the targe cell. If both v_r and
	    //v_l are in same cell then v will be small than v_l, set it then to v_l
	    //note that we also have shifted velocity to have origo at v_c (center velocity of euclidian cell)
	    const Real v_1 = max(gk * intersection_dk + intersection_min, v_l)-v_c;
	    const Real v_2 = min((gk+1) * intersection_dk + intersection_min, v_r)-v_c;

#ifdef ACC_SEMILAG_ORDER_1	    
	    //target mass is value in center of intersecting length,
	    //times length (missing x,y, but they would be cancelled
	    //anyway when we divide to get density
	    const Real target_mass = (cv + A * 0.5*(v_1+v_2)*(v_2 - v_1);
#endif
#ifdef ACC_SEMILAG_ORDER_2
            //todo, we recompute integrals, could do some reuse at least over gk loop (same with v also)
            const Real integral_1=v_1*(cv+0.5*A*v_1+B*(BB-v_1*v_1*onethird));
            const Real integral_2=v_2*(cv+0.5*A*v_2+B*(BB-v_2*v_2*onethird));
            const Real target_mass=integral_2-integral_1;
#endif
	    if (target_block < SpatialCell::max_velocity_blocks) 
	      spatial_cell->increment_value(target_block,target_cell,target_mass*i_dv);
	    
	  }
	}
      }
    }
  }
  delete [] blocks;
  return true;
}


#endif   
