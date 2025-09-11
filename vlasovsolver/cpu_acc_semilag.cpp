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

#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <phiprof.hpp>
#include "../definitions.h"

#include "cpu_acc_semilag.hpp"
#include "cpu_acc_intersections.hpp"
#include "cpu_acc_map.hpp"

#include "cpu_acc_transform.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
  Calls semi-lagrangian acceleration routines for the provided list of cells

 * @param mpiGrid DCCRG container of spatial cells
 * @param acceleratedCells vector of cells for which to perform acceleration
 * @param popID ID of the accelerated particle species.
 * @param map_order Order in which vx,vy,vz mappings are performed.
*/

void cpu_accelerate_cells(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const std::vector<CellID>& acceleratedCells,
                          const uint popID,
                          const uint map_order
   ) {
   int timerId {phiprof::initializeTimer("cell-semilag-acc")};
   int intersections_id {phiprof::initializeTimer("cell-compute-intersections")};

   #pragma omp parallel // Launch workshare region
   {
      // Calculate intersections (should be constant cost per cell)
      #pragma omp for schedule(static,1)
      for (size_t c=0; c<acceleratedCells.size(); ++c) {
         const CellID cellID = acceleratedCells[c];
         SpatialCell* SC = mpiGrid[cellID];
         Population& pop = SC->get_population(popID);

         // compute transform, forward in time and backward in time, performed in this acceleration
         pop.fwd_transform = compute_acceleration_transformation(SC, popID, pop.subcycleDt);
         pop.bwd_transform = pop.fwd_transform.inverse();
         // compute_cell_intersections(SC, popID, map_order, pop.subcycleDt, intersections_id);
      }
      #pragma omp barrier
      // Semi-Lagrangian acceleration for all cells active in this subcycle,
      // dimension-by-dimension. Dynamic cost due to varying block counts.
      #pragma omp for schedule(dynamic,1)
      for (size_t c=0; c<acceleratedCells.size(); ++c) {
         const CellID cellID = acceleratedCells[c];
         SpatialCell* SC = mpiGrid[cellID];

         phiprof::Timer semilagAccTimer {timerId};
         cpu_accelerate_cell(SC,popID,map_order);
         semilagAccTimer.stop();
      }
   }
}

/*!
  Propagates the distribution function in velocity space of given real
  space cell.

  Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
  three‐dimensional monotone and conservative semi‐Lagrangian scheme
  (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
  Meteorological Society 138.667 (2012): 1640-1651.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
 * @param map_order Order in which vx,vy,vz mappings are performed.
*/

void cpu_accelerate_cell(SpatialCell* spatial_cell,
                         const uint popID,
                         const uint map_order
   ) {

   Population& pop = spatial_cell->get_population(popID);
   switch(map_order){
      case 0: {
         //Map order XYZ
         map_1d(spatial_cell, popID, pop.intersection_x, pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk, map_order, 0); // map along x
         map_1d(spatial_cell, popID, pop.intersection_y, pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk, map_order, 1); // map along y
         map_1d(spatial_cell, popID, pop.intersection_z, pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk, map_order, 2); // map along z
         break;
      }
      case 1: {
         //Map order YZX
         map_1d(spatial_cell, popID, pop.intersection_y, pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk, map_order, 1); // map along y
         map_1d(spatial_cell, popID, pop.intersection_z, pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk, map_order, 2); // map along z
         map_1d(spatial_cell, popID, pop.intersection_x, pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk, map_order, 0); // map along x
         break;
      }
      case 2: {
         //Map order Z X Y
         map_1d(spatial_cell, popID, pop.intersection_z, pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk, map_order, 2); // map along z
         map_1d(spatial_cell, popID, pop.intersection_x, pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk, map_order, 0); // map along x
         map_1d(spatial_cell, popID, pop.intersection_y, pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk, map_order, 1); // map along y
         break;
      }
   }
}
