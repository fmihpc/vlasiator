/*
This file is part of Vlasiator.

Copyright 2012 Finnish Meteorological Institute

*/

#ifndef CPU_ACC_MAP_H
#define CPU_ACC_MAP_H

#include "algorithm"
#include "cmath"
#include "utility"

/*TODO - replace with standard library c++11 functions as soon as possible*/
#include "boost/array.hpp"
#include "boost/unordered_map.hpp"

#include "common.h"
#include "spatial_cell.hpp"

const int STENCIL_WIDTH=2; //at most equalt to WID

using namespace std;
using namespace spatial_cell;

enum InterpolationType { CONSTANT,LINEAR,PARABOLIC };


/*
  f(v)=f_cell_average +   // CONSTANT 
  A*(v-v0) +         // LINEAR
  B*(BB-(v-v0)**2)   // PARABOLIC  
*/

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

// indices in padded z block
template<typename T> inline T i_pblock(const T& i,const T& j,const T& k) {return (k+2)*WID2+j*WID+i;}
// indices in normal block
template<typename T> inline T i_block(const T& i,const T& j,const T& k) {return k*WID2+j*WID+i;}


inline void copy_block_data_z(Velocity_Block *block,Real *values){
    Velocity_Block *nbrBlock;  
    // Construct values
    // Copy averages from -z neighbour if it exists (if not, array is initialized to zero)
    nbrBlock = block->neighbors[velocity_neighbor::XCC_YCC_ZM1];
    if ( nbrBlock != NULL) {
      for (int k=-STENCIL_WIDTH; k<0; ++k) 
	for (int j=0; j<WID; ++j) 
	  for (int i=0; i<WID; ++i) {
	    values[i_pblock(i,j,k)] = nbrBlock->fx[i_blockx(i,j,k+WID)];
	  }
    }   
    // Copy volume averages of this block:
    for (int k=0; k<WID; ++k) 
      for (int j=0; j<WID; ++j) 
	for (int i=0; i<WID; ++i) {
	  values[i_pblock(i,j,k)] = block->fx[i_blockx(i,j,k)];
	}    
    // Copy averages from +z neighbour if it exists (if not, array is initialized to zero)
    nbrBlock = block->neighbors[velocity_neighbor::XCC_YCC_ZP1];
    if ( nbrBlock != NULL) {
      for (int k=WID; k<WID+STENCIL_WIDTH; ++k) 
	for (int j=0; j<WID; ++j) 
	  for (int i=0; i<WID; ++i) {
	    values[i_pblock(i,j,k)] = nbrBlock->fx[i_blockx(i,j,k-WID)];
	  }
    }
}

void map_z(SpatialCell* spatial_cell,   Real intersection, Real intersection_di, Real intersection_dj,Real intersection_dk) {
  //TODO, align arrays

  /*values used with an stencil in 1 dimension, initialized to 0 */  
  Real values[WID*WID*(WID+2*STENCIL_WIDTH)]={};
  // Make a copy of the blocklist, the blocklist will change during this algorithm
  uint*  blocks=new uint[spatial_cell->number_of_blocks];
  const int nblocks=spatial_cell->number_of_blocks;
  /* 
     Move densities from data to fx and clear data, to prepare for mapping
     Also copy blocklist since the velocity block list in spatial cell changes when we add values
  */
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
    const unsigned int block = spatial_cell->velocity_block_list[block_i];
    blocks[block_i] = block; 
    Velocity_Block* block_ptr = spatial_cell->at(block);
    for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
      block_ptr->fx[cell] = block_ptr->data[cell];
      block_ptr->data[cell] = 0.0;
    }
  }



  const Real i_dz=1.0/SpatialCell::cell_dvz;
  
  for (unsigned int block_i = 0; block_i < nblocks; block_i++) {
    Velocity_Block *block=spatial_cell->at(blocks[block_i]);
    velocity_block_indices_t block_indices=get_velocity_block_indices(blocks[block_i]);
    copy_block_data_z(block,values);
    for (int k=0; k<WID; ++k){ 
      //todo, could also compute z index and compute the start velocity of block
      const Real vl=(WID*block_indices[2]+k)*SpatialCell::cell_dvz;
      const Real vc=vl+0.5*SpatialCell::cell_dvz;
      const Real vr=vl+SpatialCell::cell_dvz;
      for (int j = 0; j < WID; ++j){
	for (int i = 0; i < WID; ++i){
	  const intersection_min=intersection +
	    (block_indices[0]*WID+i)*intersection_di + 
	    (block_indices[1]*WID+j)*intersection_dj;
	  //f is mean value
	  const real f = value[i_pblock(i,j,k)];
	  //A is slope of linear approximation
	  const Real A = slope_limiter(value[i_pblock(i,j,k-1)],
				       f,value[i_pblock(i,j,k+1)])*i_dz;
	  //left(l) and right(r) k values (global index) in the target
	  //lagrangian grid, the intersecting cells
	  int lagrangian_k_l=(vl-intersection_min)/intersection_dk; 
	  int lagrangian_k_r=(vr-intersection_min)/intersection_dk;
	  
	  //the blocks of the two lagrangian cells
	  int target_block_l= block_indices[0] + 
	    block_indices[1] * SpatialCell::vx_length + 
	    lagrangian_k_l/WID * SpatialCell::vx_length * SpatialCell::vy_length;
	  int target_block_r= block_indices[0] + 
	    block_indices[1] * SpatialCell::vx_length + 
	    lagrangian_k_r/WID * SpatialCell::vx_length * SpatialCell::vy_length;
	  //cell index in the block of the l,r target cells
	  int target_cell_l=i+j*WID+(lagrangian_k_l%WID)*WID2;
	  int target_cell_r=i+j*WID+(lagrangian_k_r%WID)*WID2;
	  
	  //the velocity between the two target cells
	  Real lagrangian_v = lagrangian_k_r * intersection_dk + intersection_min;
	  
	  if(lagrangian_k_l==lagrangian_k_r){
	    //the original cell is completely inside one lagrangian
	    //cell, let's add all to that one
	    spatial_cell->increment_value(target_block_l,target_cell_l,f);
	  }
	  else {
	    //target mass is value in center of intersecting length,
	    //times length (missing x,y, but they would be cancelled
	    //anyway when we divide to get density
	    Real target_mass_l = (f + A * (0.5*(lagrangian_v+v_l)-v_c))*(lagrangian_v - v_l);
	    Real target_mass_r = SpatialCell::cell_dvz*f-target_mass_l; //the rest
	    spatial_cell->increment_value(target_block_l,target_cell_l,target_mass_l/SpatialCell::cell_dvz);
	    spatial_cell->increment_value(target_block_r,target_cell_r,target_mass_r/SpatialCell::cell_dvz);
	  }
	}
      }
    }
  }
  
  delete [] blocks;
}
void map_x(SpatialCell* spatial_cell,   Real intersection, Real intersection_di, Real intersection_dj,Real intersection_dk) {}
void map_y(SpatialCell* spatial_cell,   Real intersection, Real intersection_di, Real intersection_dj,Real intersection_dk) {}

#endif   
