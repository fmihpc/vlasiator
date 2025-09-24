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
#include "../../velocity_mesh_parameters.h"

#include "Fluctuations.h"

using namespace spatial_cell;

Real projects::Fluctuations::rndRho, projects::Fluctuations::rndVel[3];

namespace projects {
   Fluctuations::Fluctuations() : TriAxisSearch() {}
   Fluctuations::~Fluctuations() {}
   bool Fluctuations::initialize(void) { return Project::initialize(); }

   void Fluctuations::addParameters() {
      typedef Readparameters RP;
      RP::add("Fluctuations.BX0", "Background field value (T)", 1.0e-9);
      RP::add("Fluctuations.BY0", "Background field value (T)", 2.0e-9);
      RP::add("Fluctuations.BZ0", "Background field value (T)", 3.0e-9);
      RP::add("Fluctuations.magXPertAbsAmp", "Amplitude of the magnetic perturbation along x", 1.0e-9);
      RP::add("Fluctuations.magYPertAbsAmp", "Amplitude of the magnetic perturbation along y", 1.0e-9);
      RP::add("Fluctuations.magZPertAbsAmp", "Amplitude of the magnetic perturbation along z", 1.0e-9);

      // Per-population parameters
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_Fluctuations.rho", "Number density (m^-3)", 1.0e7);
         RP::add(pop + "_Fluctuations.TemperatureX", "Temperature (K)", 2.0e6);
         RP::add(pop + "_Fluctuations.TemperatureY", "Temperature (K)", 2.0e6);
         RP::add(pop + "_Fluctuations.TemperatureZ", "Temperature (K)", 2.0e6);
         RP::add(pop + "_Fluctuations.densityPertRelAmp", "Amplitude factor of the density perturbation", 0.1);
         RP::add(pop + "_Fluctuations.velocityPertAbsAmp", "Amplitude of the velocity perturbation", 1.0e6);
         RP::add(pop + "_Fluctuations.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
      }
   }

   void Fluctuations::getParameters() {
      Project::getParameters();
      typedef Readparameters RP;
      Project::getParameters();
      RP::get("Fluctuations.BX0", this->BX0);
      RP::get("Fluctuations.BY0", this->BY0);
      RP::get("Fluctuations.BZ0", this->BZ0);
      RP::get("Fluctuations.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("Fluctuations.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("Fluctuations.magZPertAbsAmp", this->magZPertAbsAmp);

      // Per-population parameters
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         FluctuationsSpeciesParameters sP;
         RP::get(pop + "_Fluctuations.rho", sP.DENSITY);
         RP::get(pop + "_Fluctuations.TemperatureX", sP.TEMPERATUREX);
         RP::get(pop + "_Fluctuations.TemperatureY", sP.TEMPERATUREY);
         RP::get(pop + "_Fluctuations.TemperatureZ", sP.TEMPERATUREZ);
         RP::get(pop + "_Fluctuations.densityPertRelAmp", sP.densityPertRelAmp);
         RP::get(pop + "_Fluctuations.velocityPertAbsAmp", sP.velocityPertAbsAmp);
         RP::get(pop + "_Fluctuations.maxwCutoff", sP.maxwCutoff);

         speciesParams.push_back(sP);
      }
   }

   Realf Fluctuations::fillPhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, const uint nRequested) const {
      const FluctuationsSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      // const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = sP.DENSITY * (1.0 + sP.densityPertRelAmp * (0.5 - rndRho));
      Real initTx = sP.TEMPERATUREX;
      Real initTy = sP.TEMPERATUREY;
      Real initTz = sP.TEMPERATUREZ;
      const Real initV0X = sP.velocityPertAbsAmp * (0.5 - rndVel[0]);
      const Real initV0Y = sP.velocityPertAbsAmp * (0.5 - rndVel[1]);
      const Real initV0Z = sP.velocityPertAbsAmp * (0.5 - rndVel[2]);

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
                                              const Realf value = TriMaxwellianPhaseSpaceDensity(vx, vy, vz, initTx, initTy, initTz, initRho, mass);
                                              bufferData[initIndex * WID3 + k * WID2 + j * WID + i] = value;
                                              // lsum[0] += value;
                                           };
                                        },
                                        rhosum);
      return rhosum;
   }

   /* Evaluates local SpatialCell properties for the project and population,
      then evaluates the phase-space density at the given coordinates.
      Used as a probe for projectTriAxisSearch.
   */
   Realf Fluctuations::probePhaseSpace(spatial_cell::SpatialCell* cell, const uint popID, Real vx_in, Real vy_in, Real vz_in) const {
      const FluctuationsSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      // const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initRho = sP.DENSITY * (1.0 + sP.densityPertRelAmp * (0.5 - rndRho));
      Real initTx = sP.TEMPERATUREX;
      Real initTy = sP.TEMPERATUREY;
      Real initTz = sP.TEMPERATUREZ;
      const Real initV0X = sP.velocityPertAbsAmp * (0.5 - rndVel[0]);
      const Real initV0Y = sP.velocityPertAbsAmp * (0.5 - rndVel[1]);
      const Real initV0Z = sP.velocityPertAbsAmp * (0.5 - rndVel[2]);
      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;
      const Realf value = TriMaxwellianPhaseSpaceDensity(vx, vy, vz, initTx, initTy, initTz, initRho, mass);
      return value;
   }

   void Fluctuations::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {

      std::default_random_engine rndState;
      setRandomCellSeed(cell, rndState);

      this->rndRho = getRandomNumber(rndState);
      this->rndVel[0] = getRandomNumber(rndState);
      this->rndVel[1] = getRandomNumber(rndState);
      this->rndVel[2] = getRandomNumber(rndState);
   }

   void Fluctuations::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid, FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                       FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
      ConstantField bgField;
      bgField.initialize(this->BX0, this->BY0, this->BZ0);

      setBackgroundField(bgField, BgBGrid);

      if (!P::isRestart) {
         const auto localSize = BgBGrid.getLocalSize().data();

         #pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  const int64_t cellid = perBGrid.GlobalIDForCoords(x, y, z);

                  std::default_random_engine rndState;
                  setRandomSeed(cellid, rndState);

                  cell->at(fsgrids::bfield::PERBX) = this->magXPertAbsAmp * (0.5 - getRandomNumber(rndState));
                  cell->at(fsgrids::bfield::PERBY) = this->magYPertAbsAmp * (0.5 - getRandomNumber(rndState));
                  cell->at(fsgrids::bfield::PERBZ) = this->magZPertAbsAmp * (0.5 - getRandomNumber(rndState));
               }
            }
         }
      }
   }

   std::vector<std::array<Real, 3>> Fluctuations::getV0(creal x, creal y, creal z, const uint popID) const {
      std::array<Real, 3> V0{{0.0, 0.0, 0.0}};
      std::vector<std::array<Real, 3>> centerPoints;
      centerPoints.push_back(V0);
      return centerPoints;
   }

} // namespace projects
