/*
  This file is part of Vlasiator.
  Copyright 2013,2014 Finnish Meteorological Institute
*/

#ifndef CPU_ACC_SEMILAG_H
#define CPU_ACC_SEMILAG_H

#include "algorithm"
#include "cmath"
#include "utility"

/*TODO - replace with standard library c++11 functions*/
#include "boost/array.hpp"
#include "boost/unordered_map.hpp"

#include "common.h"
#include "spatial_cell.hpp"

#include <Eigen/Geometry>
#include <Eigen/Core>

#include "vlasovsolver_amr/cpu_acc_transform.hpp"
#include "vlasovsolver_amr/cpu_acc_intersections.hpp"
#include "vlasovsolver_amr/cpu_acc_map.hpp"

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

void cpu_accelerate_cell(SpatialCell* spatial_cell, uint map_order, const Real dt) {
   double t1=MPI_Wtime();
   /*compute transform, forward in time and backward in time*/
   phiprof::start("compute-transform");

   //compute the transform performed in this acceleration (ok for AMR)
   Transform<Real,3,Affine> fwd_transform= compute_acceleration_transformation(spatial_cell,dt);
   Transform<Real,3,Affine> bwd_transform= fwd_transform.inverse();
   phiprof::stop("compute-transform");

   map_order=0;
   static int dim = map_order;

   static int counter=0;

   if (counter % 2 == 0) {
   switch (map_order) {
    case 0:
      if (dim == 0) map_1d(spatial_cell, fwd_transform, bwd_transform,0,0);
      if (dim == 1) map_1d(spatial_cell, fwd_transform, bwd_transform,1,1);
      //map_1d(spatial_cell, fwd_transform, bwd_transform,2,2);
      break;
    case 1:
      if (dim == 1) map_1d(spatial_cell, fwd_transform, bwd_transform,1,0);
      //map_1d(spatial_cell, fwd_transform, bwd_transform,2,1);
      //map_1d(spatial_cell, fwd_transform, bwd_transform,0,2);
      if (dim == 0) map_1d(spatial_cell, fwd_transform, bwd_transform,0,1);
      break;
    case 2:
      map_1d(spatial_cell, fwd_transform, bwd_transform,2,0);
      map_1d(spatial_cell, fwd_transform, bwd_transform,0,1);
      map_1d(spatial_cell, fwd_transform, bwd_transform,1,2);
      break;
    default:
      map_1d(spatial_cell, fwd_transform, bwd_transform,2,0);
      map_1d(spatial_cell, fwd_transform, bwd_transform,0,1);
      map_1d(spatial_cell, fwd_transform, bwd_transform,1,2);
      break;
   }
   }

   if (counter % 2 != 0) {
      phiprof::start("mesh coarsening");
      amr_ref_criteria::Base* refCriterion = getObjectWrapper().amrVelRefCriteria.create(Parameters::amrVelRefCriterion);
      if (refCriterion != NULL) {
         refCriterion->initialize("");
         spatial_cell->coarsen_blocks(refCriterion);
         delete refCriterion;
      }
      phiprof::stop("mesh coarsening");
   } else {
      ++dim;
      if (dim == 2) dim=0;
   }
   ++counter;

   double t2=MPI_Wtime();
   spatial_cell->parameters[CellParams::LBWEIGHTCOUNTER] += t2 - t1;
}

#endif

