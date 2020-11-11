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
#ifndef CPU_TRANS_MAP_AMR_H
#define CPU_TRANS_MAP_AMR_H

#include <vector>

#include "vec.h"
#include "../common.h"
#include "../spatial_cell.hpp"


struct setOfPencils {

   uint N; // Number of pencils in the set
   uint sumOfLengths;
   std::vector< uint > lengthOfPencils; // Lengths of pencils
   std::vector< CellID > ids; // List of cells
   std::vector< uint > idsStart; // List of where a pencil's CellIDs start in the ids array
   std::vector< Realv > x,y; // x,y - position
   std::vector< bool > periodic;
   std::vector< std::vector<uint> > path; // Path taken through refinement levels

   setOfPencils() {
      
      N = 0;
      sumOfLengths = 0;
   }

   void addPencil(std::vector<CellID> idsIn, Real xIn, Real yIn, bool periodicIn, std::vector<uint> pathIn) {

      N++;
      sumOfLengths += idsIn.size();
      lengthOfPencils.push_back(idsIn.size());
      idsStart.push_back(ids.size());
      ids.insert(ids.end(),idsIn.begin(),idsIn.end());
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
      ids.erase(ids.begin() + ibeg, ids.begin() + ibeg + lengthOfPencils[pencilId]);
      idsStart.erase(idsStart.begin() + pencilId);

      N--;
      sumOfLengths -= lengthOfPencils[pencilId];
      lengthOfPencils.erase(lengthOfPencils.begin() + pencilId);
         
   }

   std::vector<CellID> getIds(const uint pencilId) const {
      
      if (pencilId > N) {
         std::vector<CellID> idsEmpty;
         return idsEmpty;
      }
      
      // Use vector range constructor
      std::vector<CellID>::const_iterator ibeg = ids.begin() + idsStart[pencilId];
      std::vector<CellID>::const_iterator iend = ibeg + lengthOfPencils[pencilId];
      std::vector<CellID> idsOut(ibeg, iend);
      
      return idsOut;
   }

   // Split one pencil into up to four pencils covering the same space.
   // dx and dy are the dimensions of the original pencil.
   void split(const uint myPencilId, const Realv dx, const Realv dy) {

      auto myIds = this->getIds(myPencilId);

      // Find paths that members of this pencil may have in other pencils (can happen)
      // so that we don't add duplicates.
      std::vector<int> existingSteps;
      
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
                        existingSteps.push_back(theirStep);
                     }
                  }
               }
            }
         }
      }

      bool firstPencil = true;
      const auto copy_of_path = path.at(myPencilId);
      const auto copy_of_x = x.at(myPencilId);
      const auto copy_of_y = y.at(myPencilId);

      // Add those pencils whose steps dont already exist in the pencils struct
      for (int step = 0; step < 4; ++step) {
         if (std::any_of(existingSteps.begin(), existingSteps.end(), [step](int i){return step == i;})) {
            continue;
         }

         Realv signX = 1.0;
         Realv signY = 1.0;
         
         if(step < 2) {
            signY = -1.0;
         }

         if(step % 2 == 0) {
            signX = -1.0;
         }
         
         auto myX = copy_of_x + signX * 0.25 * dx;
         auto myY = copy_of_y + signY * 0.25 * dy;

         if(firstPencil) {
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


bool trans_map_1d_amr(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const std::vector<CellID>& localPropagatedCells,
                  const std::vector<CellID>& remoteTargetCells,
                  std::vector<uint>& nPencils,
                  const uint dimension,
                  const Realv dt,
                  const uint popID);


void update_remote_mapping_contribution_amr(dccrg::Dccrg<spatial_cell::SpatialCell,
                                            dccrg::Cartesian_Geometry>& mpiGrid,
                                            const uint dimension,
                                            int direction,
                                            const uint popID);


#endif
