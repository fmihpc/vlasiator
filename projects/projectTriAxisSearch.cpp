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
    * This assumes that the velocity space is isotropic (same resolution in vx, vy, vz).
    */
   uint TriAxisSearch::findBlocksToInitialize(SpatialCell* cell,const uint popID) const {
      vmesh::VelocityMesh *vmesh = cell->get_velocity_mesh(popID);

      vmesh::GlobalID *GIDbuffer;
      const auto& vblocks_ini = cell->get_velocity_grid_length(popID);
      const uint blocksCount = vblocks_ini[0]*vblocks_ini[1]*vblocks_ini[2];
      #ifdef USE_GPU
      // Host-pinned memory buffer, max possible size
      CHK_ERR( gpuMallocHost((void**)&GIDbuffer,blocksCount*sizeof(vmesh::GlobalID)) );
      #endif
      // Non-GPU: insert directly into vmesh

      std::set<vmesh::GlobalID> singleSet;
      bool search;
      unsigned int counterX, counterY, counterZ;

      creal minValue = cell->getVelocityBlockMinValue(popID);
      // And how big a buffer do we add to the edges?
      uint buffer = 2;
      if (WID > 4 && blocksCount > 8) {
         // With WID8 a two-block buffer increases memory requirements too much.
         // However, we allow extra buffer for very minimal v-spaces.
         buffer = 1;
      }
      // How much below the sparsity can a cell be to still be included?
      creal tolerance = 0.1;

      creal x = cell->parameters[CellParams::XCRD];
      creal y = cell->parameters[CellParams::YCRD];
      creal z = cell->parameters[CellParams::ZCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      // creal dvxCell = cell->get_velocity_grid_cell_size(popID)[0];
      // creal dvyCell = cell->get_velocity_grid_cell_size(popID)[1];
      // creal dvzCell = cell->get_velocity_grid_cell_size(popID)[2];

      const size_t vxblocks_ini = cell->get_velocity_grid_length(popID)[0];
      const size_t vyblocks_ini = cell->get_velocity_grid_length(popID)[1];
      const size_t vzblocks_ini = cell->get_velocity_grid_length(popID)[2];

      vmesh::LocalID LID = 0;
      const vector<std::array<Real, 3>> V0 = this->getV0(x+0.5*dx, y+0.5*dy, z+0.5*dz, popID);
      const bool singlePeak = ( V0.size() == 1 );
      // Loop over possible V peaks
      for (int peak = 0; peak < V0.size(); ++peak) {
         std::array<Real, 3> currentV0 {V0[peak]};
         std::array<Real, 3> vRadiusSquared {probePhaseSpaceInv(cell, popID, tolerance * minValue, peak)};
         std::array<Real, 3> vRadius {sqrt(vRadiusSquared[0]), sqrt(vRadiusSquared[1]), sqrt(vRadiusSquared[2])};

         // Assuming here blocks are the smallest around V0
         Real dV[3];
         vmesh::GlobalID block0 {vmesh->getGlobalID(currentV0.data())};
         vmesh->getBlockSize(block0, dV);
         Real counterX = (vRadius[0] / dV[0]);
         Real counterY = (vRadius[1] / dV[1]);
         Real counterZ = (vRadius[2] / dV[2]);
         
         #ifndef USE_GPU // non-GPU mesh resizing
         // sphere volume is 4/3 pi r^3, approximate that 5*counterX*counterY*counterZ is enough.
         vmesh::LocalID currentMaxSize = LID + 5*counterX*counterY*counterZ;
         vmesh->setNewSize(currentMaxSize);
         GIDbuffer = vmesh->getGrid()->data();
         #endif

         // Block listing
         Real V_crds[3];
         for (uint kv=0; kv<vzblocks_ini; ++kv) {
            for (uint jv=0; jv<vyblocks_ini; ++jv) {
               for (uint iv=0; iv<vxblocks_ini; ++iv) {
                  const vmesh::GlobalID GID = vmesh->getGlobalID(iv,jv,kv);
                  vmesh->getBlockCoordinates(GID,V_crds);
                  vmesh->getBlockSize(GID, dV);

                  // Check block center point
                  for (int i = 0; i < 3; ++i) {
                     // TODO why 2 * dV??
                     V_crds[i] += 0.5 * dV[i] - currentV0[i];
                  }

                  #ifndef USE_GPU // non-GPU mesh resizing
                  if (LID >= currentMaxSize) {
                     currentMaxSize = LID + counterX*counterY*counterZ;
                     vmesh->setNewSize(currentMaxSize);
                     GIDbuffer = vmesh->getGrid()->data();
                  }
                  #endif
                  // Add this block only if it doesn't exist yet
                  if (V_crds[0] * V_crds[0] / vRadiusSquared[0] + V_crds[1] * V_crds[1] / vRadiusSquared[1] + V_crds[2] * V_crds[2] / vRadiusSquared[2] < 1 && singleSet.count(GID)==0) {
                     singleSet.insert(GID);
                     GIDbuffer[LID] = GID;
                     LID++;
                  }
               } // vxblocks_ini
            } // vyblocks_ini
         } // vzblocks_ini
      } // iteration over V0's
      // Set final size of vmesh
      cell->get_population(popID).N_blocks = LID;

      #ifdef USE_GPU
      // Copy data from CPU to GPU
      cell->dev_resize_vmesh(popID,LID);
      vmesh::GlobalID *GIDtarget = vmesh->getGrid()->data();
      gpuStream_t stream = gpu_getStream();
      CHK_ERR( gpuMemcpyAsync(GIDtarget, GIDbuffer, LID*sizeof(vmesh::GlobalID), gpuMemcpyHostToDevice, stream));
      CHK_ERR( gpuStreamSynchronize(stream) );
      CHK_ERR( gpuFreeHost(GIDbuffer));
      #else
      // Resize vmesh down to final size
      vmesh->setNewSize(LID);
      #endif

      return LID;
   }

} // namespace projects
