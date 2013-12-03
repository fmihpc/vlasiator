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

using namespace std;
using namespace spatial_cell;

enum InterpolationType { CONSTANT,LINEAR,PARABOLIC };
enum ReconstructionParams {V0,A,B,NUM_REC_PARAMS };
enum CommonReconstructionParams {BB,NUM_COMMON_REC_PARAMS};

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

/*!
  1D mapping from source to target. 
  TODO - needs to bew essentially redone
  - add stuff to normal blocks, indexing just means a differnet thing, ew can keep on swapping between fx and data...
  - do stuff one block at a time
  - vectorize as now over the dims perpendicular to 1d mapping
  
  \param source_indices Indices of cells in source along dimension
  \param source_value   Value in source cells
  \param source_min_coordinate  Min coordinate for source cells, source cell i goes from l=source_indices[i]*source_delta + source_min_coordinate  to r=l+source_delta
  \param source_delta Size of source cells.
  \param target_coordinates 
  \param mapped_values
*/

void mapping_1d(const InterpolationType interpolationType,
		const std::vector< uint >& source_indices,
		const std::vector< Real >& source_value,
		const Real source_min_coordinate,
		const Real source_delta,
		const Real target_min_coordinate,
		const Real target_delta,
		std::vector< uint >& target_indices,
		std::vector< Real >& target_values) {  
  switch(interpolationType) {
  case CONSTANT:
    for(uint i=0;i<source_indices.size();i++) {
      const Real c=source_value[i];
    }
    break;
    
  case LINEAR:
    for(uint i = 0; i < source_indices.size(); i++ ) {
      const Real c = source_value[i];
      const Real l = (i>0 && source_indices[i-1]==source_indices[i]-1)?source_value[i-1]:0.0;
      const Real r = (i<source_indices.size()-1 && source_indices[i+1]==source_indices[i]+1)?source_value[i+1]:0.0;
      const Real A = slope_limiter(l,c,r)/source_delta;
      const Real vl = (source_indices[i])*source_delta+source_min_coordinate;
      const Real vc = vl+0.5*source_delta;
      const Real vr = vl+source_delta;
            
      /*compute the two overlapping target cells */
      const int target_index_1 = (vl-target_min_coordinate)/target_delta; //first target cell
      const int target_index_2 = target_index_1+1;                        //second target cell
      const Real target_v_mid = target_index_2*target_delta+target_min_coordinate;  // mid betweem target cells
      
      const Real mass_1 = c + (0.5*(target_v_mid - vl) - vc ) * A;
      const Real mass_2 = c + (0.5*(vr-target_v_mid) - vc ) * A;

      if(target_indices.back() == target_index_1 )
	target_values[target_indices.size()-1]+=mass_1;
      else{
	target_values.push_back(mass_1);
	target_indices.push_back(target_index_1);
      }
      target_values.push_back(mass_2);
      target_indices.push_back(target_index_2);
    }
  }
}

void map_z(SpatialCell* spatial_cell,   Real intersection, Real intersection_di, Real intersection_dj,Real intersection_dk) {}
void map_x(SpatialCell* spatial_cell,   Real intersection, Real intersection_di, Real intersection_dj,Real intersection_dk) {}
void map_y(SpatialCell* spatial_cell,   Real intersection, Real intersection_di, Real intersection_dj,Real intersection_dk) {}

#endif   
