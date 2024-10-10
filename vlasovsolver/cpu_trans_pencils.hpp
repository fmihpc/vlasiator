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
#ifndef CPU_TRANS_PENCILS_H
#define CPU_TRANS_PENCILS_H

#include <vector>
#include "vec.h"
#include "../common.h"
#include "../spatial_cell_wrapper.hpp"

struct setOfPencils {

   uint N; // Number of pencils in the set
   uint sumOfLengths;
   std::vector< uint > lengthOfPencils; // Lengths of pencils (including stencil cells)
   std::vector< CellID > ids; // List of pencil cells (including stencil cells)
   std::vector< uint > idsStart; // List of where a pencil's CellIDs start in the ids array
   std::vector< Realf > sourceDZ; // Widths of source cells
   std::vector< Realf > targetRatios; // Pencil to target cell area ratios of target cells
   std::vector< Real > x,y; // x,y - position
   std::vector< bool > periodic;
   std::vector< std::vector<uint> > path; // Path taken through refinement levels

   setOfPencils() {
      N = 0;
      sumOfLengths = 0;
   }

   void removeAllPencils() {
      N = 0;
      sumOfLengths = 0;
      sumOfLengths = 0;
      lengthOfPencils.clear();
      idsStart.clear();
      ids.clear();
      sourceDZ.clear();
      targetRatios.clear();
      x.clear();
      y.clear();
      periodic.clear();
      path.clear();
   }

   void addPencil(std::vector<CellID> idsIn, Real xIn, Real yIn, bool periodicIn, std::vector<uint> pathIn) {
      N++;
      // If necessary, add the zero cells to the beginning and end
      if (idsIn.front() != 0) {
            idsIn.insert(idsIn.begin(),VLASOV_STENCIL_WIDTH,0);
      }
      if (idsIn.back() != 0) {
         for (int i = 0; i < VLASOV_STENCIL_WIDTH; i++)
         {
            idsIn.push_back(0);
         }
      }
      sumOfLengths += idsIn.size();
      lengthOfPencils.push_back(idsIn.size());
      idsStart.push_back(ids.size());
      ids.insert(ids.end(),idsIn.begin(),idsIn.end());
      sourceDZ.resize(sumOfLengths);
      targetRatios.resize(sumOfLengths);
      x.push_back(xIn);
      y.push_back(yIn);
      periodic.push_back(periodicIn);
      path.push_back(pathIn);
   }

   void removePencil(const uint pencilId) {
      x.erase(x.begin() + pencilId);
      y.erase(y.begin() + pencilId);
      periodic.erase(periodic.begin() + pencilId);
      path.erase(path.begin() + pencilId);

      uint ibeg = idsStart[pencilId];
      ids.erase(ids.begin() + ibeg, ids.begin() + ibeg + lengthOfPencils[pencilId] + 2*VLASOV_STENCIL_WIDTH);
      targetRatios.erase(targetRatios.begin() + ibeg, targetRatios.begin() + ibeg + lengthOfPencils[pencilId] + 2*VLASOV_STENCIL_WIDTH);
      sourceDZ.erase(sourceDZ.begin() + ibeg, sourceDZ.begin() + ibeg + lengthOfPencils[pencilId] + 2*VLASOV_STENCIL_WIDTH);
      idsStart.erase(idsStart.begin() + pencilId);

      N--;
      sumOfLengths -= lengthOfPencils[pencilId];
      lengthOfPencils.erase(lengthOfPencils.begin() + pencilId);
   }

   std::vector<CellID> getIds(const uint pencilId) const {
      if (pencilId >= N) {
         std::vector<CellID> idsEmpty;
         return idsEmpty;
      }
      // Use vector range constructor. Only return actual pencil ids, not the stencils at the ends
      std::vector<CellID>::const_iterator ibeg = ids.begin() + idsStart[pencilId] + VLASOV_STENCIL_WIDTH;
      std::vector<CellID>::const_iterator iend = ibeg + lengthOfPencils[pencilId] - 2*VLASOV_STENCIL_WIDTH;
      std::vector<CellID> idsOut(ibeg, iend);
      return idsOut;
   }

   // Split one pencil into up to four pencils covering the same space.
   // dx and dy are the dimensions of the original pencil.
   void split(const uint myPencilId, const Real dx, const Real dy) {
      auto myIds = this->getIds(myPencilId);

      // Find paths that members of this pencil may have in other pencils (can happen)
      // so that we don't add duplicates.
      std::vector<int> existingSteps;

#pragma omp parallel for
      for (uint theirPencilId = 0; theirPencilId < this->N; ++theirPencilId) {
         if(theirPencilId == myPencilId) continue;
         auto theirIds = this->getIds(theirPencilId);
         for (auto theirId : theirIds) {
            for (auto myId : myIds) {
               if (myId == theirId) {
                  std::vector<uint> theirPath = this->path.at(theirPencilId);
                  std::vector<uint> myPath = this->path.at(myPencilId);
                  if(theirPath.size() > myPath.size()) {
                     bool samePath = true;
                     for (uint i = 0; i < myPath.size(); ++i) {
                        if(myPath.at(i) != theirPath.at(i)) {
                           samePath = false;
                        }
                     }

                     if(samePath) {
                        uint theirStep = theirPath.at(myPath.size());
#pragma omp critical
                        {
                           existingSteps.push_back(theirStep);
                        }
                     }
                  }
               }
            }
         }
      } // end parallel region

      bool firstPencil = true;
      const auto copy_of_path = path.at(myPencilId);
      const auto copy_of_x = x.at(myPencilId);
      const auto copy_of_y = y.at(myPencilId);

      // Add those pencils whose steps dont already exist in the pencils struct
      for (int step = 0; step < 4; ++step) {
         if (std::any_of(existingSteps.begin(), existingSteps.end(), [step](int i){return step == i;})) {
            continue;
         }

         Real signX = 1.0;
         Real signY = 1.0;

         if(step < 2) {
            signY = -1.0;
         }

         if(step % 2 == 0) {
            signX = -1.0;
         }

         auto myX = copy_of_x + signX * 0.25 * dx;
         auto myY = copy_of_y + signY * 0.25 * dy;

         if (firstPencil) {
            //TODO: set x and y correctly. Right now they are not used anywhere.
            path.at(myPencilId).push_back(step);
            x.at(myPencilId) = myX;
            y.at(myPencilId) = myY;
            firstPencil = false;
         } else {
            auto myPath = copy_of_path;
            myPath.push_back(step);
            addPencil(myIds, myX, myY, periodic.at(myPencilId), myPath);
         }
      }
   }
};
// Note: Splitting does not handle target or source cells, as those are computed after all pencil splitting has concluded.

bool do_translate_cell(spatial_cell::SpatialCell* SC);

// grid.cpp calls this function to both find seed cells and build pencils for all dimensions
void prepareSeedIdsAndPencils(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
// find seed cells and build pencils for one dimension
void prepareSeedIdsAndPencils(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              const uint dimension);

// pencils used for AMR translation
extern std::array<setOfPencils,3> DimensionPencils;
extern std::array<std::unordered_set<CellID>,3> DimensionTargetCells;

// Ghost translation cell lists (no interim comms)
void prepareGhostTranslationCellLists(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      const std::vector<CellID>& localPropagatedCells);

// defined in cpu_trans_map_amr.cpp
extern std::unordered_set<CellID> ghostTranslate_sources_x;
extern std::unordered_set<CellID> ghostTranslate_sources_y;
extern std::unordered_set<CellID> ghostTranslate_sources_z;
extern std::unordered_set<CellID> ghostTranslate_active_x;
extern std::unordered_set<CellID> ghostTranslate_active_y;
extern std::unordered_set<CellID> ghostTranslate_active_z;

#endif
