/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
   phiprof::Timer computeTransformTimer {"compute-transform"};

   //compute the transform performed in this acceleration
   Transform<Real,3,Affine> fwd_transform= compute_acceleration_transformation(spatial_cell,dt);
   Transform<Real,3,Affine> bwd_transform= fwd_transform.inverse();
   computeTransformTimer.stop();

   // NOTE: This is now in a debugging / testing state. The propagator 
   // only does one thing each time step. Currently this does
   // step=0 accel vx
   //      1 coarsen mesh
   //      2 accel vy
   //      3 coarsen mesh
   //      4 accel vz
   //      5 coarsen mesh
   // (repeat)
   
   // It is then easy to see the effect of each step in the output vlsv files.
   
   // BEGIN TEST
   map_order=0;
   static int dim = map_order;
   static int counter=0;
   if (counter % 2 == 0) {
   // END TEST

   switch (map_order) {
    case 0: // x -> y -> z
      // BEGIN TEST
      if (dim == 0) map_1d(spatial_cell, fwd_transform, bwd_transform,0,0);
      if (dim == 1) map_1d(spatial_cell, fwd_transform, bwd_transform,1,1);
      if (dim == 2) map_1d(spatial_cell, fwd_transform, bwd_transform,2,2);      
      // END TEST
      //map_1d(spatial_cell, fwd_transform, bwd_transform,0,0);
      //map_1d(spatial_cell, fwd_transform, bwd_transform,1,1);
      //map_1d(spatial_cell, fwd_transform, bwd_transform,2,2);
      break;
    case 1: // y -> z -> x
      map_1d(spatial_cell, fwd_transform, bwd_transform,1,0);
      map_1d(spatial_cell, fwd_transform, bwd_transform,2,1);
      map_1d(spatial_cell, fwd_transform, bwd_transform,0,2);
      break;
    case 2: // z -> x -> y
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
   // BEGIN TEST
   }
   // END TEST

   // NOTE: Mesh coarsening might be needed after each acceleration substep 
   
   // BEGIN TEST
   if (counter % 2 != 0) {
      phiprof::Timer timer {"mesh coarsening"};
      vamr_ref_criteria::Base* refCriterion = getObjectWrapper().amrVelRefCriteria.create(Parameters::vamrVelRefCriterion);
      if (refCriterion != NULL) {
         refCriterion->initialize("");
         spatial_cell->coarsen_blocks(refCriterion);
         delete refCriterion;
      }
   } else {
      ++dim;
      if (dim == 2) dim=0;
   }
   ++counter;
   // END TEST
   
   double t2=MPI_Wtime();
   spatial_cell->parameters[CellParams::LBWEIGHTCOUNTER] += t2 - t1;
}

#endif

