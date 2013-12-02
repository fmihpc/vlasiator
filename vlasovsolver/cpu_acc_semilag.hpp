/*
This file is part of Vlasiator.

Copyright 2012 Finnish Meteorological Institute

*/

#ifndef CPU_ACC_SEMILAG_H
#define CPU_ACC_SEMILAG_H

#include "algorithm"
#include "cmath"
#include "utility"

/*TODO - replace with standard library c++11 functions as soon as possible*/
#include "boost/array.hpp"
#include "boost/unordered_map.hpp"

#include "common.h"
#include "spatial_cell.hpp"

#include <Eigen/Geometry>
#include <Eigen/Core>

#include "vlasovsolver/cpu_acc_compute_downstream_blocks.hpp"
#include "vlasovsolver/cpu_acc_transform.hpp"
#include "vlasovsolver/cpu_acc_intersections.hpp"

using namespace std;
using namespace spatial_cell;
using namespace Eigen;


/*!  

Propagates the distribution function in velocity space of given real
space cell.

Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
three‐dimensional monotone and conservative semi‐Lagrangian scheme
(SLICE‐3D) for transport problems." Quarterly Journal of the Royal
Meteorological Society 138.667 (2012): 1640-1651.

*/

void cpu_accelerate_cell(SpatialCell* spatial_cell,const Real dt) {

   /*compute transform, forward in time and backward in time*/
   phiprof::start("compute-transform");
   //compute the transform performed in this acceleration
   Transform<Real,3,Affine> fwd_transform= compute_acceleration_transformation(spatial_cell,dt);
   Transform<Real,3,Affine> bwd_transform= fwd_transform.inverse();
   phiprof::stop("compute-transform");

 
  // Make a copy of the blocklist, these are the current Eulerian cells
  std::vector<unsigned int> upstream_blocks;
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
    upstream_blocks.push_back(spatial_cell->velocity_block_list[block_i]);
  }
  
  /* Compute masses based on densities, and store it in fx. data is
     cleared to make way for the actual accelerated values.
  */
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
    const unsigned int block = spatial_cell->velocity_block_list[block_i];
    Velocity_Block* block_ptr = spatial_cell->at(block);
    const Real volume=block_ptr->parameters[BlockParams::DVX]*
      block_ptr->parameters[BlockParams::DVY]*
         block_ptr->parameters[BlockParams::DVZ];
    for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
      block_ptr->fx[cell] = block_ptr->data[cell]*volume;
      block_ptr->data[cell] = 0.0;
    }
  }
  
  
 
  //TODO, we might want to compute a index instead of using an array, let's keep it for now and benchmark later on
  //intersections_x  key: y (euclidian),z (lagrangian) index of cell, value: vector of pair(x index of cell (Lagrangian), x where cell starts)
  Real intersection_x;
  Real intersection_x_di;
  Real intersection_x_dj;
  Real intersection_x_dk;
    
  
  compute_intersections_x(spatial_cell, bwd_transform, fwd_transform,
			  intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk);
  
  //   compute_mapping
  exit(1);
  
   
   
   
 }


   

#endif

