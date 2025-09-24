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

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"
#include "Diffusion.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   Diffusion::Diffusion() : Project() {}
   Diffusion::~Diffusion() {}

   bool Diffusion::initialize(void) { return Project::initialize(); }

   void Diffusion::addParameters() {
      typedef Readparameters RP;
      RP::add("Diffusion.B0", "Background field value (T)", 1.0e-9);

      // Per-population parameters
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_Diffusion.rho", "Number density (m^-3)", 1.0e7);
         RP::add(pop + "_Diffusion.Temperature", "Temperature (K)", 2.0e6);
         RP::add(pop + "_Diffusion.Scale_x", "Scale length in x (m)", 100000.0);
         RP::add(pop + "_Diffusion.Scale_y", "Scale length in y (m)", 100000.0);
      }
   }

   void Diffusion::getParameters() {
      Project::getParameters();

      typedef Readparameters RP;
      RP::get("Diffusion.B0", this->B0);

      // Per-population parameters
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         DiffusionSpeciesParameters sP;

         RP::get(pop + "_Diffusion.rho", sP.DENSITY);
         RP::get(pop + "_Diffusion.Temperature", sP.TEMPERATURE);
         RP::get(pop + "_Diffusion.Scale_x", sP.SCA_X);
         RP::get(pop + "_Diffusion.Scale_y", sP.SCA_Y);

         speciesParams.push_back(sP);
      }
   }

   Realf Diffusion::fillPhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, const uint nRequested) const {
      const DiffusionSpeciesParameters& sP = speciesParams[popID];
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      // Fetch spatial cell center coordinates
      const Real x = cell->parameters[CellParams::XCRD] + 0.5 * cell->parameters[CellParams::DX];
      const Real y = cell->parameters[CellParams::YCRD] + 0.5 * cell->parameters[CellParams::DY];
      const Real z = cell->parameters[CellParams::ZCRD] + 0.5 * cell->parameters[CellParams::DZ];
      const Real initV0X = 0;
      const Real initV0Y = 0;
      const Real initV0Z = 0;
      const Real kb = physicalconstants::K_B;
      const Real initRho = sP.DENSITY;
      const Real initT = sP.TEMPERATURE;
      const Real initScaY = sP.SCA_Y;
      const Real initScaX = sP.SCA_X;

      #ifdef USE_GPU
      vmesh::VelocityMesh* vmesh = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh* vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      #endif
      // Loop over blocks
      Realf rhosum = 0;
      arch::parallel_reduce<arch::null>({WID, WID, WID, nRequested},
                                        ARCH_LOOP_LAMBDA(const uint i, const uint j, const uint k, const uint initIndex, Realf* lsum) {
                                           vmesh::GlobalID* GIDlist = vmesh->getGrid()->data();
                                           Realf* bufferData = VBC->getData();
                                           const vmesh::GlobalID blockGID = GIDlist[initIndex];
                                           // Calculate parameters for new block
                                           Real blockCoords[6];
                                           vmesh->getBlockInfo(blockGID, &blockCoords[0]);
                                           creal vxBlock = blockCoords[0];
                                           creal vyBlock = blockCoords[1];
                                           creal vzBlock = blockCoords[2];
                                           creal dvxCell = blockCoords[3];
                                           creal dvyCell = blockCoords[4];
                                           creal dvzCell = blockCoords[5];
                                           ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
                                              creal vx = vxBlock + (i + 0.5) * dvxCell - initV0X;
                                              creal vy = vyBlock + (j + 0.5) * dvyCell - initV0Y;
                                              creal vz = vzBlock + (k + 0.5) * dvzCell - initV0Z;
                                              const Realf value = initRho * pow(mass / (2.0 * M_PI * kb * initT), 1.5) *
                                                                  (5.0 * exp(-(pow(x, 2.0) / pow(initScaX, 2.0) + pow(y, 2.0) / pow(initScaY, 2.0))) *
                                                                       exp(-mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * kb * initT)) +
                                                                   exp(-mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * kb * initT)));
                                              bufferData[initIndex * WID3 + k * WID2 + j * WID + i] = value;
                                              // lsum[0] += value;
                                           };
                                        },
                                        rhosum);
      return rhosum;
   }

   void Diffusion::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

   void Diffusion::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid, FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
      ConstantField bgField;
      bgField.initialize(0, 0, this->B0); // bg bx, by,bz
      setBackgroundField(bgField, BgBGrid);
   }
} // namespace projects
