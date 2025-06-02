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

#include "spatial_cell_cpu.hpp"
#include "block_adjust_cpu.hpp"
#include "../object_wrapper.h"
#include "../velocity_mesh_parameters.h"

namespace spatial_cell {
/** Bulk call over listed cells of spatial grid
    Prepares the content / no-content velocity block lists
    for all requested cells, for the requested popID
**/
   void update_velocity_block_content_lists(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const vector<CellID>& cells,
      const uint popID) {

      if (cells.size()==0) {
         return;
      }

      // int computeId {phiprof::initializeTimer("Compute with_content_list")};
#pragma omp parallel
      {
         // phiprof::Timer timer {computeId};
#pragma omp for schedule(dynamic)
         for (uint i=0; i<cells.size(); ++i) {
            mpiGrid[cells[i]]->updateSparseMinValue(popID);
            mpiGrid[cells[i]]->update_velocity_block_content_lists(popID);
         }
         //  timer.stop();
      } // end parallel region
   }

   void adjust_velocity_blocks_in_cells(
      dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const vector<CellID>& cellsToAdjust,
      const uint popID
      ) {

      const size_t n_cells = cellsToAdjust.size();
      if (n_cells==0) {
         return;
      }

//      int adjustId {phiprof::initializeTimer("Adjusting blocks")};
#pragma omp parallel
      {
//         phiprof::Timer timer {adjustId};
#pragma omp for schedule(dynamic)
         for (size_t i=0; i < n_cells; ++i) {
            Real density_pre_adjust=0.0;
            Real density_post_adjust=0.0;
            CellID cell_id=cellsToAdjust[i];
            SpatialCell* cell = mpiGrid[cell_id];
            vector<SpatialCell*> neighbor_ptrs;
            // gather spatial neighbor list and gather vector with pointers to cells
            const auto* neighbors = mpiGrid.get_neighbors_of(cell_id, Neighborhoods::NEAREST);
            // Note: at AMR refinement boundaries this can cause blocks to propagate further
            // than absolutely required. Face neighbours, however, are not enough as we must
            // account for diagonal propagation.

            std::unordered_set<CellID> uniqueNeighbors;
            uniqueNeighbors.reserve(neighbors->size());
            // find only unique neighbor cells
            for ( const auto& [neighbor_id, dir] : *neighbors) {
               if ((neighbor_id != 0) && (neighbor_id != cell_id)) {
                  uniqueNeighbors.insert(neighbor_id);
               }
            }
            neighbor_ptrs.reserve(uniqueNeighbors.size());
            for ( const CellID neighbor_id : uniqueNeighbors) {
               neighbor_ptrs.push_back(mpiGrid[neighbor_id]);
            }

            if (getObjectWrapper().particleSpecies[popID].sparse_conserve_mass) {
               for (size_t i=0; i<cell->get_number_of_velocity_blocks(popID)*WID3; ++i) {
                  density_pre_adjust += cell->get_data(popID)[i];
               }
            }

            cell->adjust_velocity_blocks(neighbor_ptrs,popID);

            if (getObjectWrapper().particleSpecies[popID].sparse_conserve_mass) {
               for (size_t i=0; i<cell->get_number_of_velocity_blocks(popID)*WID3; ++i) {
                  density_post_adjust += cell->get_data(popID)[i];
               }
               if (density_post_adjust != 0.0) {
                  for (size_t i=0; i<cell->get_number_of_velocity_blocks(popID)*WID3; ++i) {
                     cell->get_data(popID)[i] *= density_pre_adjust/density_post_adjust;
                  }
               }
            }
         }
//         timer.stop();
      } // end parallel region
   }

} //namespace
