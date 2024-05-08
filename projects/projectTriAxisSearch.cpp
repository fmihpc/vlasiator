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
#include "projectTriAxisSearch.h"
#include "../object_wrapper.h"

using namespace std;
using namespace spatial_cell;

using namespace std;

namespace projects {
   /*!
    * WARNING This assumes that the velocity space is isotropic (same resolution in vx, vy, vz).
    */
   std::vector<vmesh::GlobalID> TriAxisSearch::findBlocksToInitialize(SpatialCell* cell,const uint popID) const {
      set<vmesh::GlobalID> blocksToInitialize;
      bool search;
      unsigned int counterX, counterY, counterZ;

      creal minValue = cell->getVelocityBlockMinValue(popID);
      // And how big a buffer do we add to the edges?
      const uint buffer = 2;
      // How much below the sparsity can a cell be to still be included?
      creal tolerance = 0.1;

      creal x = cell->parameters[CellParams::XCRD];
      creal y = cell->parameters[CellParams::YCRD];
      creal z = cell->parameters[CellParams::ZCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      creal dvxBlock = cell->get_velocity_grid_block_size(popID)[0];
      creal dvyBlock = cell->get_velocity_grid_block_size(popID)[1];
      creal dvzBlock = cell->get_velocity_grid_block_size(popID)[2];
      creal dvxCell = cell->get_velocity_grid_cell_size(popID)[0];
      creal dvyCell = cell->get_velocity_grid_cell_size(popID)[1];
      creal dvzCell = cell->get_velocity_grid_cell_size(popID)[2];

      const size_t vxblocks_ini = cell->get_velocity_grid_length(popID)[0];
      const size_t vyblocks_ini = cell->get_velocity_grid_length(popID)[1];
      const size_t vzblocks_ini = cell->get_velocity_grid_length(popID)[2];

      const vector<std::array<Real, 3>> V0 = this->getV0(x+0.5*dx, y+0.5*dy, z+0.5*dz, popID);
      for (vector<std::array<Real, 3>>::const_iterator it = V0.begin(); it != V0.end(); it++) {
         // VX search
         search = true;
         counterX = 0;
         while (search) {
            if ( (tolerance * minValue >
                  calcPhaseSpaceDensity(x, y, z, dx, dy, dz,
                                        it->at(0) + counterX*dvxBlock, it->at(1), it->at(2),
                                        dvxCell, dvyCell, dvzCell, popID)
                  || counterX > vxblocks_ini ) ) {
               search = false;
            }
            counterX++;
         }
         counterX+=buffer;
         Real vRadiusSquared = (Real)counterX*(Real)counterX*dvxBlock*dvxBlock;

         // VY search
         search = true;
         counterY = 0;
         while(search) {
            if ( (tolerance * minValue >
                  calcPhaseSpaceDensity(x, y, z, dx, dy, dz,
                                        it->at(0), it->at(1) + counterY*dvyBlock, it->at(2),
                                        dvxCell, dvyCell, dvzCell, popID)
                  || counterY > vyblocks_ini ) ) {
               search = false;
            }
            counterY++;
         }
         counterY+=buffer;
         vRadiusSquared = max(vRadiusSquared, (Real)counterY*(Real)counterY*dvyBlock*dvyBlock);

         // VZ search
         search = true;
         counterZ = 0;
         while(search) {
            if ( (tolerance * minValue >
                  calcPhaseSpaceDensity(x, y, z, dx, dy, dz,
                                        it->at(0), it->at(1), it->at(2) + counterZ*dvzBlock,
                                        dvxCell, dvyCell, dvzCell, popID)
                  || counterZ > vzblocks_ini ) ) {
               search = false;
            }
            counterZ++;
         }
         counterZ+=buffer;
         vRadiusSquared = max(vRadiusSquared, (Real)counterZ*(Real)counterZ*dvzBlock*dvzBlock);

         // Block listing
         Real V_crds[3];
         for (uint kv=0; kv<vzblocks_ini; ++kv) {
            for (uint jv=0; jv<vyblocks_ini; ++jv) {
               for (uint iv=0; iv<vxblocks_ini; ++iv) {
                  vmesh::GlobalID blockIndices[3];
                  blockIndices[0] = iv;
                  blockIndices[1] = jv;
                  blockIndices[2] = kv;
                  const vmesh::GlobalID blockGID = cell->get_velocity_block(popID,blockIndices);

                  cell->get_velocity_block_coordinates(popID,blockGID,V_crds);
                  V_crds[0] += (0.5*dvxBlock - it->at(0) );
                  V_crds[1] += (0.5*dvyBlock - it->at(1) );
                  V_crds[2] += (0.5*dvzBlock - it->at(2) );
                  Real R2 = ((V_crds[0])*(V_crds[0])
                             + (V_crds[1])*(V_crds[1])
                             + (V_crds[2])*(V_crds[2]));

                  if (R2 < vRadiusSquared) {
                     blocksToInitialize.insert(blockGID);
                  }
               } // vxblocks_ini
            } // vyblocks_ini
         } // vzblocks_ini
      } // iteration over V0's

      vector<vmesh::GlobalID> returnVector;
      for (set<vmesh::GlobalID>::const_iterator it=blocksToInitialize.begin(); it!=blocksToInitialize.end(); ++it) {
         returnVector.push_back(*it);
      }
      return returnVector;
   }

} // namespace projects
