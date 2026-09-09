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

// #include <dccrg.hpp>
// #include <dccrg_cartesian_geometry.hpp>
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
                          const std::vector<AccelerationPayload>& acceleratedCells,
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
         const CellID cellID = acceleratedCells[c].cellptr->get_cellid();
         Population& pop = acceleratedCells[c].cellptr->get_population(popID, acceleratedCells[c].timeclass);
         compute_cell_intersections(acceleratedCells[c].cellptr, popID, map_order, pop.subcycleDt, intersections_id, acceleratedCells[c].timeclass);
      }
      #pragma omp barrier
      // Semi-Lagrangian acceleration for all cells active in this subcycle,
      // dimension-by-dimension. Dynamic cost due to varying block counts.
      #pragma omp for schedule(dynamic,1)
      for (size_t c=0; c<acceleratedCells.size(); ++c) {
         const CellID cellID = acceleratedCells[c].cellptr->get_cellid();
         phiprof::Timer semilagAccTimer {timerId};
         cpu_accelerate_cell(acceleratedCells[c].cellptr,popID,map_order,acceleratedCells[c].dt,acceleratedCells[c].timeclass);
         semilagAccTimer.stop();
      }
   }
}


/*!
  Compute the number of subcycles needed for the acceleration of the particle
  species in the spatial cell. Note that one should first prepare to
  accelerate the cell with prepareAccelerateCell.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
*/

// uint getAccelerationSubcycles(SpatialCell* spatial_cell, Real dt, const uint popID)
// {
//    //return max( convert<uint>(ceil(dt*spatial_cell->CellParams[CELLPARAMS::TIMECLASSDT] / spatial_cell->get_max_v_dt(popID))), 1u);
//    return max( convert<uint>(ceil(dt / spatial_cell->get_max_v_dt(popID))), 1u);
// }
// /*!
//   Compute the number of subcycles needed from maxVdt and target dt.

//  * @param spatial_cell Spatial cell containing the accelerated population.
//  * @param popID ID of the accelerated particle species.
// */

// uint getAccelerationSubcycles(Real maxVdt, Real dt)
// {
//    //return max( convert<uint>(ceil(dt*spatial_cell->CellParams[CELLPARAMS::TIMECLASSDT] / spatial_cell->get_max_v_dt(popID))), 1u);
//    return max( convert<uint>(ceil(dt / maxVdt)), 1u);
// }

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
                         const uint map_order,
                         const Real& dt,
                         int timeclass) {
   //double t1 = MPI_Wtime();
   vmesh::VelocityMesh* vmeshPtr = NULL; 
   vmesh::VelocityBlockContainer* blockContainerPtr = NULL;
   // Main branch: check if time-ghost data required, handle that and recurse as needed
   // .... of course this is now substepping unnecessarily...
   // std::cout << spatial_cell->parameters[CellParams::CELLID] << "c Accelerate, initial refs" << " vmesh " << &vmesh << " blockContainer " << &blockContainer << "\n";

   vmeshPtr = spatial_cell->get_velocity_mesh(popID,timeclass);
//    std::cout << "Accelerating a vmesh with meshID " << vmeshPtr->getMesh() << " at "<< spatial_cell->get_cellid() << " size: " << vmeshPtr->size() << "\n";
   blockContainerPtr = spatial_cell->get_velocity_blocks(popID,timeclass);
    
   vmesh::VelocityMesh& vmesh = *vmeshPtr;
   vmesh::VelocityBlockContainer& blockContainer = *blockContainerPtr;


   Population& pop = spatial_cell->get_population(popID, timeclass);
   switch(map_order){
      case 0: {
         //Map order XYZ
         map_1d(spatial_cell, popID, vmeshPtr, blockContainerPtr, pop.intersection_x,
                pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk,0, timeclass); // map along x
         map_1d(spatial_cell, popID, vmeshPtr, blockContainerPtr, pop.intersection_y,
                pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk,1, timeclass); // map along y
         map_1d(spatial_cell, popID, vmeshPtr, blockContainerPtr, pop.intersection_z,
                pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk,2, timeclass); // map along z
         break;
      }
      case 1: {
         //Map order YZX
         map_1d(spatial_cell, popID, vmeshPtr, blockContainerPtr, pop.intersection_y,
                pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk,1, timeclass); // map along y
         map_1d(spatial_cell, popID, vmeshPtr, blockContainerPtr, pop.intersection_z,
                pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk,2, timeclass); // map along z
         map_1d(spatial_cell, popID, vmeshPtr, blockContainerPtr, pop.intersection_x,
                pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk,0, timeclass); // map along x
         break;
      }
      case 2: {
         //Map order Z X Y
         map_1d(spatial_cell, popID, vmeshPtr, blockContainerPtr, pop.intersection_z,
                pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk,2, timeclass); // map along z
         map_1d(spatial_cell, popID, vmeshPtr, blockContainerPtr, pop.intersection_x,
                pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk,0, timeclass); // map along x
         map_1d(spatial_cell, popID, vmeshPtr, blockContainerPtr, pop.intersection_y,
                pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk,1, timeclass); // map along y
         break;
      }
   }
}
