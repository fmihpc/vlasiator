/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * 2017-2025 University of Helsinki
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
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../object_wrapper.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"

#include "LossCone.h"

using namespace spatial_cell;

Real projects::LossCone::rndRho, projects::LossCone::rndVel[3];


namespace projects {
   LossCone::LossCone(): TriAxisSearch() { }
   LossCone::~LossCone() { }
   bool LossCone::initialize(void) {return Project::initialize();}

   void LossCone::addParameters() {
      typedef Readparameters RP;
      RP::add("LossCone.BX0", "Background field value (T)", 1.0e-9);
      RP::add("LossCone.BY0", "Background field value (T)", 2.0e-9);
      RP::add("LossCone.BZ0", "Background field value (T)", 3.0e-9);
      RP::add("LossCone.magXPertAbsAmp", "Amplitude of the magnetic perturbation along x", 1.0e-9);
      RP::add("LossCone.magYPertAbsAmp", "Amplitude of the magnetic perturbation along y", 1.0e-9);
      RP::add("LossCone.magZPertAbsAmp", "Amplitude of the magnetic perturbation along z", 1.0e-9);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_LossCone.rho", "Number density (m^-3)", 1.0e7);
         RP::add(pop + "_LossCone.TemperatureX", "Temperature (K)", 2.0e6);
         RP::add(pop + "_LossCone.TemperatureY", "Temperature (K)", 2.0e6);
         RP::add(pop + "_LossCone.TemperatureZ", "Temperature (K)", 2.0e6);
         RP::add(pop + "_LossCone.densityPertRelAmp", "Amplitude factor of the density perturbation", 0.1);
         RP::add(pop + "_LossCone.VX0", "Initial bulk velocity in x-direction", 0.0);
         RP::add(pop + "_LossCone.VY0", "Initial bulk velocity in y-direction", 0.0);
         RP::add(pop + "_LossCone.VZ0", "Initial bulk velocity in z-direction", 0.0);
         RP::add(pop + "_LossCone.velocityPertAbsAmp", "Amplitude of the velocity perturbation", 1.0e6);
         RP::add(pop + "_LossCone.muLimit", "Cutoff value for pitch-cosine mu positive and negative)", 0.5);
      }
   }

   void LossCone::getParameters() {
      Project::getParameters();
      typedef Readparameters RP;
      Project::getParameters();
      RP::get("LossCone.BX0", this->BX0);
      RP::get("LossCone.BY0", this->BY0);
      RP::get("LossCone.BZ0", this->BZ0);
      RP::get("LossCone.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("LossCone.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("LossCone.magZPertAbsAmp", this->magZPertAbsAmp);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         LossConeSpeciesParameters sP;
         RP::get(pop + "_LossCone.rho", sP.DENSITY);
         RP::get(pop + "_LossCone.VX0", sP.V0[0]);
         RP::get(pop + "_LossCone.VY0", sP.V0[1]);
         RP::get(pop + "_LossCone.VZ0", sP.V0[2]);
         RP::get(pop + "_LossCone.TemperatureX", sP.TEMPERATUREX);
         RP::get(pop + "_LossCone.TemperatureY", sP.TEMPERATUREY);
         RP::get(pop + "_LossCone.TemperatureZ", sP.TEMPERATUREZ);
         RP::get(pop + "_LossCone.densityPertRelAmp", sP.densityPertRelAmp);
         RP::get(pop + "_LossCone.velocityPertAbsAmp", sP.velocityPertAbsAmp);
         RP::get(pop + "_LossCone.muLimit", sP.muLimit);
         speciesParams.push_back(sP);
      }
   }

   Realf LossCone::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      const LossConeSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      // const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real initRho = sP.DENSITY * (1.0 + sP.densityPertRelAmp * (0.5 - rndRho));
      const Real initTx = sP.TEMPERATUREX;
      const Real initTy = sP.TEMPERATUREY;
      const Real initTz = sP.TEMPERATUREZ;
      const Real initV0X = sP.V0[0] + sP.velocityPertAbsAmp * (0.5 - rndVel[0] );
      const Real initV0Y = sP.V0[1] + sP.velocityPertAbsAmp * (0.5 - rndVel[1] );
      const Real initV0Z = sP.V0[2] + sP.velocityPertAbsAmp * (0.5 - rndVel[2] );
      const Real muLimit = abs(sP.muLimit);

      #ifdef USE_GPU
      vmesh::VelocityMesh *vmesh = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh *vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      #endif
      // Loop over blocks
      Realf rhosum = 0;
      arch::parallel_reduce<arch::null>(
         {WID, WID, WID, nRequested},
         ARCH_LOOP_LAMBDA (const uint i, const uint j, const uint k, const uint initIndex, Realf *lsum ) {
            vmesh::GlobalID *GIDlist = vmesh->getGrid()->data();
            Realf* bufferData = VBC->getData();
            const vmesh::GlobalID blockGID = GIDlist[initIndex];
            // Calculate parameters for new block
            Real blockCoords[6];
            vmesh->getBlockInfo(blockGID,&blockCoords[0]);
            creal vxBlock = blockCoords[0];
            creal vyBlock = blockCoords[1];
            creal vzBlock = blockCoords[2];
            creal dvxCell = blockCoords[3];
            creal dvyCell = blockCoords[4];
            creal dvzCell = blockCoords[5];
            ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
               creal vx = vxBlock + (i+0.5)*dvxCell - initV0X;
               creal vy = vyBlock + (j+0.5)*dvyCell - initV0Y;
               creal vz = vzBlock + (k+0.5)*dvzCell - initV0Z;

               // TODO: Use Eigen vectors, get magnetic field as well and calculate components from that
               Real vpara = vx;
               // Real vperp = sqrt(vy*vy + vz*vz);
               Real modv = sqrt(vx*vx + vy*vy + vz*vz);
               Real mu    = vpara / modv;

               Real value = 0;
               // Only fill outside losscone
               if (mu > -muLimit && mu < muLimit) {
                  value += TriMaxwellianPhaseSpaceDensity(vx,vy,vz,initTx,initTy,initTz,initRho,mass);
               }
               bufferData[initIndex*WID3 + k*WID2 + j*WID + i] = value;
               //lsum[0] += value;
            };
         }, rhosum);
      return rhosum;
   }

   /* Evaluates local SpatialCell properties for the project and population,
      then evaluates the phase-space density at the given coordinates.
      Used as a probe for projectTriAxisSearch.
   */
   Realf LossCone::probePhaseSpace(spatial_cell::SpatialCell *cell,
                                        const uint popID,
                                        Real vx_in, Real vy_in, Real vz_in
      ) const {
      const LossConeSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      // const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      // const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      // const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const Real initRho = sP.DENSITY * (1.0 + sP.densityPertRelAmp * (0.5 - rndRho));
      const Real initTx = sP.TEMPERATUREX;
      const Real initTy = sP.TEMPERATUREY;
      const Real initTz = sP.TEMPERATUREZ;
      const Real initV0X = sP.V0[0] + sP.velocityPertAbsAmp * (0.5 - rndVel[0] );
      const Real initV0Y = sP.V0[1] + sP.velocityPertAbsAmp * (0.5 - rndVel[1] );
      const Real initV0Z = sP.V0[2] + sP.velocityPertAbsAmp * (0.5 - rndVel[2] );
      const Real muLimit = abs(sP.muLimit);
      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;

      // Probe function should not account for mu Limit
      Real value = TriMaxwellianPhaseSpaceDensity(vx,vy,vz,initTx,initTy,initTz,initRho,mass);
      return value;
   }

   void LossCone::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      creal z = cellParams[CellParams::ZCRD];
      creal dz = cellParams[CellParams::DZ];

      CellID cellID = (int) ((x - Parameters::xmin) / dx) +
         (int) ((y - Parameters::ymin) / dy) * Parameters::xcells_ini +
         (int) ((z - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;

      std::default_random_engine rndState;
      setRandomCellSeed(cell,rndState);

      this->rndRho=getRandomNumber(rndState);
      this->rndVel[0]=getRandomNumber(rndState);
      this->rndVel[1]=getRandomNumber(rndState);
      this->rndVel[2]=getRandomNumber(rndState);
   }

   void LossCone::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(this->BX0,
                         this->BY0,
                         this->BZ0);

      setBackgroundField(bgField, BgBGrid);

      if(!P::isRestart) {
         const auto localSize = BgBGrid.getLocalSize().data();

         #pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  const int64_t cellid = perBGrid.GlobalIDForCoords(x, y, z);

                  std::default_random_engine rndState;
                  setRandomSeed(cellid,rndState);

                  cell->at(fsgrids::bfield::PERBX) = this->magXPertAbsAmp * (0.5 - getRandomNumber(rndState));
                  cell->at(fsgrids::bfield::PERBY) = this->magYPertAbsAmp * (0.5 - getRandomNumber(rndState));
                  cell->at(fsgrids::bfield::PERBZ) = this->magZPertAbsAmp * (0.5 - getRandomNumber(rndState));
               }
            }
         }
      }
   }

   std::vector<std::array<Real, 3> > LossCone::getV0(
      creal x,
      creal y,
      creal z,
      const uint popID
   ) const {
      std::array<Real, 3> V0 {{0.0, 0.0, 0.0}};
      std::vector<std::array<Real, 3> > centerPoints;
      centerPoints.push_back(V0);
      return centerPoints;
   }

} // namespace projects
