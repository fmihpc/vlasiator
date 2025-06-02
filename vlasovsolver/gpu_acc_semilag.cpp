/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
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

#include "gpu_acc_semilag.hpp"
#include "cpu_acc_intersections.hpp"
#include "gpu_acc_map.hpp"

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

void gpu_accelerate_cells(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const std::vector<CellID>& acceleratedCells,
                          const uint popID,
                          const uint map_order
   ) {
   int intersections_id {phiprof::initializeTimer("cell-compute-intersections")};
   uint gpuMaxBlockCount = 0;
   // Calculate intersections (should be constant cost per cell)
   #pragma omp parallel
   {
      uint threadGpuMaxBlockCount = 0;
      #pragma omp for schedule(static,1)
      for (size_t c=0; c<acceleratedCells.size(); ++c) {
         const CellID cellID = acceleratedCells[c];
         SpatialCell* SC = mpiGrid[cellID];
         Population& pop = SC->get_population(popID);
         compute_cell_intersections(SC, popID, map_order, pop.subcycleDt, intersections_id);

         const vmesh::VelocityMesh* vmesh = SC->get_velocity_mesh(popID);
         const uint blockCount = vmesh->size();
         threadGpuMaxBlockCount = std::max(threadGpuMaxBlockCount,blockCount);
      }
      #pragma omp critical
      {
         gpuMaxBlockCount = std::max(gpuMaxBlockCount,threadGpuMaxBlockCount);
      }
   }

   // Ensure accelerator has enough temporary memory allocated
   phiprof::Timer verificationTimer {"gpu ACC allocation verifications"};
   gpu_vlasov_allocate(gpuMaxBlockCount);
   gpu_acc_allocate(gpuMaxBlockCount);
   verificationTimer.stop();

   // Semi-Lagrangian acceleration for all cells active in this subcycle,
   // dimension-by-dimension. Dynamic cost due to varying block counts.
   int timerId {phiprof::initializeTimer("cell-semilag-acc")};
   #pragma omp parallel for schedule(dynamic,1)
   for (size_t c=0; c<acceleratedCells.size(); ++c) {
      const CellID cellID = acceleratedCells[c];
      SpatialCell* SC = mpiGrid[cellID];

      phiprof::Timer semilagAccTimer {timerId};
      gpu_accelerate_cell(SC,popID,map_order);
      semilagAccTimer.stop();
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

void gpu_accelerate_cell(SpatialCell* spatial_cell,
                         const uint popID,
                         const uint map_order
   ) {

   Population& pop = spatial_cell->get_population(popID);
   switch(map_order){
      case 0: {
         //Map order XYZ
         gpu_acc_map_1d(spatial_cell, popID, pop.intersection_x,
                pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk,0); // map along x
         gpu_acc_map_1d(spatial_cell, popID, pop.intersection_y,
                pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk,1); // map along y
         gpu_acc_map_1d(spatial_cell, popID, pop.intersection_z,
                pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk,2); // map along z
         break;
      }
      case 1: {
         //Map order YZX
         gpu_acc_map_1d(spatial_cell, popID, pop.intersection_y,
                pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk,1); // map along y
         gpu_acc_map_1d(spatial_cell, popID, pop.intersection_z,
                pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk,2); // map along z
         gpu_acc_map_1d(spatial_cell, popID, pop.intersection_x,
                pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk,0); // map along x
         break;
      }
      case 2: {
         //Map order Z X Y
         gpu_acc_map_1d(spatial_cell, popID, pop.intersection_z,
                pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk,2); // map along z
         gpu_acc_map_1d(spatial_cell, popID, pop.intersection_x,
                pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk,0); // map along x
         gpu_acc_map_1d(spatial_cell, popID, pop.intersection_y,
                pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk,1); // map along y
         break;
      }
   }
}
