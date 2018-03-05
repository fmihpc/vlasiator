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

#include <cstdlib>
#include <unordered_map>
#include <map>

#include "fs_cache.h"
#include "fs_common.h"

using namespace std;

extern map<CellID,uint> existingCellsFlags;

// ***** Definitions for static fs_cache::CacheContainer members ***** //

fs_cache::CacheContainer cacheContainer;
long int fs_cache::CacheContainer::cacheCalculatedStep = -1;
vector<fs_cache::CellCache> fs_cache::CacheContainer::localCellsCache;
vector<uint16_t> fs_cache::CacheContainer::boundaryCellsWithLocalNeighbours;
vector<uint16_t> fs_cache::CacheContainer::boundaryCellsWithRemoteNeighbours;
vector<uint16_t> fs_cache::CacheContainer::cellsWithLocalNeighbours;
vector<uint16_t> fs_cache::CacheContainer::cellsWithRemoteNeighbours;
vector<uint16_t> fs_cache::CacheContainer::local_NOT_DO_NOT_COMPUTE;
vector<uint16_t> fs_cache::CacheContainer::local_NOT_SYSBOUND_DO_NOT_COMPUTE;

namespace fs_cache {
           
   CellCache::CellCache() {
      for (int i=0; i<27; ++i) cells[i] = NULL;
   }
   
   void calculateCache(
      dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID>& cells
   ) {
      if (Parameters::tstep == cacheContainer.cacheCalculatedStep) return;

      cacheContainer.clear();
      unordered_map<CellID,uint16_t> globalToLocalMap;

      for (size_t c=0; c<cells.size(); ++c) {
         const CellID cellID = cells[c];
         globalToLocalMap[cellID] = c;
         spatial_cell::SpatialCell* cell = mpiGrid[cellID];

         if (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            cacheContainer.local_NOT_SYSBOUND_DO_NOT_COMPUTE.push_back(c);
            cacheContainer.local_NOT_DO_NOT_COMPUTE.push_back(c);
         } else if (cell->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE) {
            cacheContainer.local_NOT_DO_NOT_COMPUTE.push_back(c);
         }

         CellCache cellCache;
         cellCache.cellID             = cellID;
         cellCache.existingCellsFlags = existingCellsFlags[cellID];
         cellCache.sysBoundaryFlag    = cell->sysBoundaryFlag;

         // Get pointers to all cells withing 3x3x3 centered around this cell
         // NOTE: the version below does the same thing as the commented out code block
         // but it is MUCH faster!
         //for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
         //   if (i == 0 && (j == 0 && k == 0)) {
         //      cell = mpiGrid[cellID];
         //   } else {
         //      const CellID nbrID = mpiGrid.get_neighbors_of_at_offset(cellID,i,j,k)[0];
         //      cell = mpiGrid[nbrID];
         //   }
         //   cellCache.cells[fs_cache::calculateNbrID(1+i,1+j,1+k)] = cell;
         //}
         const vector<CellID>* nbrs = mpiGrid.get_neighbors_of(cellID,FIELD_SOLVER_NEIGHBORHOOD_ID);         
         int i_nbrs = 0;
         for (int n=0; n<27; ++n) {
            if (n == fs_cache::calculateNbrID(1,1,1)) {
               cellCache.cells[n] = mpiGrid[cellID];
               continue;
            }
            cellCache.cells[n] = mpiGrid[(*nbrs)[i_nbrs]];
            ++i_nbrs;
         }
         cacheContainer.localCellsCache.push_back(cellCache);
      }

      // Calculate caches needed in propagateMagneticFieldSimple
      vector<CellID> temp;
      getBoundaryCellList(mpiGrid,
                          mpiGrid.get_local_cells_not_on_process_boundary(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID),
                          temp);
      for (size_t c=0; c<temp.size(); ++c) {
         cacheContainer.boundaryCellsWithLocalNeighbours.push_back(globalToLocalMap[temp[c]]);
      }

      temp.clear();
      getBoundaryCellList(mpiGrid,
                          mpiGrid.get_local_cells_on_process_boundary(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID),
                          temp);
      for (size_t c=0; c<temp.size(); ++c) {
         cacheContainer.boundaryCellsWithRemoteNeighbours.push_back(globalToLocalMap[temp[c]]);
      }

      // Calculate caches needed in calculateBVOLDerivativesSimple (derivatives.cpp)
      // NOTE: DO_NOT_COMPUTE cells have been removed
      const vector<uint64_t> cellsWithLocalNeighbours
        = mpiGrid.get_local_cells_not_on_process_boundary(FIELD_SOLVER_NEIGHBORHOOD_ID);
      for (size_t c=0; c<cellsWithLocalNeighbours.size(); ++c) {
         if (mpiGrid[cellsWithLocalNeighbours[c]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         cacheContainer.cellsWithLocalNeighbours.push_back(globalToLocalMap[cellsWithLocalNeighbours[c]]);
      }

      const vector<uint64_t> cellsWithRemoteNeighbours
        = mpiGrid.get_local_cells_on_process_boundary(FIELD_SOLVER_NEIGHBORHOOD_ID);
      for (size_t c=0; c<cellsWithRemoteNeighbours.size(); ++c) {
         if (mpiGrid[cellsWithRemoteNeighbours[c]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         cacheContainer.cellsWithRemoteNeighbours.push_back(globalToLocalMap[cellsWithRemoteNeighbours[c]]);
      }

      cacheContainer.cacheCalculatedStep = Parameters::tstep;
   }
   
   void CacheContainer::clear() {
      vector<fs_cache::CellCache>().swap(localCellsCache);
      vector<uint16_t>().swap(local_NOT_DO_NOT_COMPUTE);
      vector<uint16_t>().swap(local_NOT_SYSBOUND_DO_NOT_COMPUTE);
      vector<uint16_t>().swap(boundaryCellsWithLocalNeighbours);
      vector<uint16_t>().swap(boundaryCellsWithRemoteNeighbours);
      vector<uint16_t>().swap(cellsWithRemoteNeighbours);
      vector<uint16_t>().swap(cellsWithLocalNeighbours);
   }
   
   CacheContainer& getCache() {return cacheContainer;}
      
} // namespace fs_cache
