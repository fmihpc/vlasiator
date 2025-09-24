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
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"
#include "../../velocity_mesh_parameters.h"

#include "Shock.h"

namespace projects {
   Shock::Shock() : Project() {}
   Shock::~Shock() {}

   bool Shock::initialize(void) { return Project::initialize(); }

   void Shock::addParameters() {
      typedef Readparameters RP;
      RP::add("Shock.BX0", "Background field value (T)", 1.0e-9);
      RP::add("Shock.BY0", "Background field value (T)", 2.0e-9);
      RP::add("Shock.BZ0", "Background field value (T)", 3.0e-9);
      RP::add("Shock.EX0", "Background electric field", 0.0);
      RP::add("Shock.VX0", "Bulk velocity in x", 0.0);
      RP::add("Shock.VY0", "Bulk velocity in y", 0.0);
      RP::add("Shock.VZ0", "Bulk velocuty in z", 0.0);
      RP::add("Shock.rho", "Number density (m^-3)", 1.0e7);
      RP::add("Shock.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Shock.magPertAmp", "Amplitude of the magnetic perturbation", 1.0e-9);
      RP::add("Shock.densityPertAmp", "Amplitude factor of the density perturbation", 0.1);
      RP::add("Shock.velocityPertAmp", "Amplitude of the velocity perturbation", 1.0e6);
      RP::add("Shock.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
      RP::add("Shock.Scale_x", "Scale length in x (m)", 2.0e6);
      RP::add("Shock.Scale_y", "Scale length in y (m)", 2.0e6);
      RP::add("Shock.Sharp_Y", "Sharpness of tannh", 0.1);
   }

   void Shock::getParameters() {
      Project::getParameters();
      typedef Readparameters RP;

      if (getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("Shock.BX0", this->BX0);
      RP::get("Shock.BY0", this->BY0);
      RP::get("Shock.BZ0", this->BZ0);
      RP::get("Shock.EX0", this->EX0);
      RP::get("Shock.VX0", this->VX0);
      RP::get("Shock.VY0", this->VY0);
      RP::get("Shock.VZ0", this->VZ0);
      RP::get("Shock.rho", this->DENSITY);
      RP::get("Shock.Temperature", this->TEMPERATURE);
      RP::get("Shock.magPertAmp", this->magPertAmp);
      RP::get("Shock.densityPertAmp", this->densityPertAmp);
      RP::get("Shock.velocityPertAmp", this->velocityPertAmp);
      RP::get("Shock.maxwCutoff", this->maxwCutoff);
      RP::get("Shock.Scale_x", this->SCA_X);
      RP::get("Shock.Scale_y", this->SCA_Y);
      RP::get("Shock.Sharp_Y", this->Sharp_Y);
   }

   Realf Shock::fillPhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, const uint nRequested) const {
      // const speciesParameters& sP = this->speciesParams[popID];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = this->DENSITY;
      Real initT = this->TEMPERATURE;
      const Real initV0X = this->VX0;
      const Real initV0Y = this->VY0;
      const Real initV0Z = this->VZ0;

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
                                              const Realf value = MaxwellianPhaseSpaceDensity(vx, vy, vz, initT, initRho, mass);
                                              bufferData[initIndex * WID3 + k * WID2 + j * WID + i] = value;
                                              // lsum[0] += value;
                                           };
                                        },
                                        rhosum);
      return rhosum;
   }

   void Shock::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

   void Shock::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid, FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
      setBackgroundFieldToZero(BgBGrid);

      if (!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();

#pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);

                  cell->at(fsgrids::bfield::PERBX) = 0.0;
                  cell->at(fsgrids::bfield::PERBY) = 0.0;
                  cell->at(fsgrids::bfield::PERBZ) = this->BZ0 * (3.0 + 2.0 * tanh((xyz[1] - Parameters::ymax / 2.0) / (this->Sharp_Y * Parameters::ymax)));
               }
            }
         }
      }
   }
} // namespace projects
