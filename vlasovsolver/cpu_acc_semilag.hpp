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

#include "vlasovsolver/cpu_acc_transform.hpp"
#include "vlasovsolver/cpu_acc_intersections.hpp"
#include "vlasovsolver/cpu_acc_map.hpp"

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

 
  

  
  Real intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk;  
  Real intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk;  
  Real intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk;  
  compute_intersections_z(spatial_cell, bwd_transform, fwd_transform,
			  intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk);
  compute_intersections_x(spatial_cell, bwd_transform, fwd_transform,
			  intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk);
  compute_intersections_y(spatial_cell, bwd_transform, fwd_transform,
			  intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk);
  
  
  map_z(spatial_cell, intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk); /*< map along z, mass from fx->data, clear fx*/
  map_x(spatial_cell, intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk); /*< map along z, mass from fx->data, clear fx*/
  map_y(spatial_cell, intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk); /*< map along z, mass from fx->data, clear fx*/
  

  /* 
     Compute densities from the mass we have created
  */
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
    const unsigned int block = spatial_cell->velocity_block_list[block_i];
    Velocity_Block* block_ptr = spatial_cell->at(block);
    const Real volume=block_ptr->parameters[BlockParams::DVX]*
      block_ptr->parameters[BlockParams::DVY]*
      block_ptr->parameters[BlockParams::DVZ];
    for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
      block_ptr->data[cell]/=volume;
    }
  }
  
  //   compute_mapping
  exit(1);
  
   
   
   
 }


   

#endif

