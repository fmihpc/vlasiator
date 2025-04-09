/*
 * This file is part of Vlasiator.
 * Copyright 2024-2025 University of Helsinki 
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

#include <vector>
#include "../definitions.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include "../object_wrapper.h"
//#include <stdint.h>
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

//using namespace std;
using namespace spatial_cell;

void reduce_vlasov_dt(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                      const vector<CellID>& cells,
                      Real (&dtMaxLocal)[3]) {

   phiprof::Timer computeTimestepTimer {"compute-vlasov-timestep"};
   const Real HALF = 0.5;

   for (vector<CellID>::const_iterator cell_id = cells.begin(); cell_id != cells.end(); ++cell_id) {
      SpatialCell* cell = mpiGrid[*cell_id];
      const Real dx = cell->parameters[CellParams::DX];
      const Real dy = cell->parameters[CellParams::DY];
      const Real dz = cell->parameters[CellParams::DZ];
      cell->parameters[CellParams::MAXRDT] = numeric_limits<Real>::max();

      for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
         cell->set_max_r_dt(popID, numeric_limits<Real>::max());
         const Real EPS = numeric_limits<Real>::min() * 1000;

         const uint nBlocks = cell->get_number_of_velocity_blocks(popID);
         if (nBlocks==0) {
            continue;
         }
         #ifdef USE_GPU
         const vmesh::VelocityBlockContainer *blockContainer = cell->dev_get_velocity_blocks(popID);
         #else
         const vmesh::VelocityBlockContainer *blockContainer = cell->get_velocity_blocks(popID);
         #endif

         Real threadMin = std::numeric_limits<Real>::max();
         arch::parallel_reduce<arch::min>({2, nBlocks},
            ARCH_LOOP_LAMBDA (uint i, const uint blockLID, Real *lthreadMin) -> void{
               i = i * (WID - 1); // ie, i == 0, i == WID - 1
               const Real* blockParams = blockContainer->getParameters(popID);
               const Real Vx =
                   blockParams[blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] +
                   (i + HALF) * blockParams[blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX] + EPS;
               const Real Vy =
                   blockParams[blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] +
                   (i + HALF) * blockParams[blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY] + EPS;
               const Real Vz =
                   blockParams[blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] +
                   (i + HALF) * blockParams[blockLID * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ] + EPS;

               const Real dt_max_cell = min({dx / fabs(Vx), dy / fabs(Vy), dz / fabs(Vz)});
               lthreadMin[0] = min(dt_max_cell,lthreadMin[0]);
         }, threadMin);
         cell->set_max_r_dt(popID, threadMin);
         cell->parameters[CellParams::MAXRDT] = min(cell->get_max_r_dt(popID), cell->parameters[CellParams::MAXRDT]);
      } // end loop over popID

      if (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
          (cell->sysBoundaryLayer == 1 && cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY)) {
         // spatial fluxes computed also for boundary cells
         dtMaxLocal[0] = min(dtMaxLocal[0], cell->parameters[CellParams::MAXRDT]);
      }

      if (cell->parameters[CellParams::MAXVDT] != 0 &&
          (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
           (P::vlasovAccelerateMaxwellianBoundaries && cell->sysBoundaryFlag == sysboundarytype::MAXWELLIAN))) {
         // acceleration only done on non-boundary cells
         dtMaxLocal[1] = min(dtMaxLocal[1], cell->parameters[CellParams::MAXVDT]);
      }
   }
}
