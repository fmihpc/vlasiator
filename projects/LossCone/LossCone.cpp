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
         RP::add(pop + "_LossCone.velocityPertAbsAmp", "Amplitude of the velocity perturbation", 1.0e6);
         RP::add(pop + "_LossCone.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
         RP::add(pop + "_LossCone.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
         RP::add(pop + "_LossCone.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
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
         RP::get(pop + "_LossCone.TemperatureX", sP.TEMPERATUREX);
         RP::get(pop + "_LossCone.TemperatureY", sP.TEMPERATUREY);
         RP::get(pop + "_LossCone.TemperatureZ", sP.TEMPERATUREZ);
         RP::get(pop + "_LossCone.densityPertRelAmp", sP.densityPertRelAmp);
         RP::get(pop + "_LossCone.velocityPertAbsAmp", sP.velocityPertAbsAmp);
         RP::get(pop + "_LossCone.nSpaceSamples", sP.nSpaceSamples);
         RP::get(pop + "_LossCone.nVelocitySamples", sP.nVelocitySamples);
         RP::get(pop + "_LossCone.maxwCutoff", sP.maxwCutoff);

         speciesParams.push_back(sP);
      }
   }
   
   Real LossCone::getDistribValue(creal& vx,creal& vy, creal& vz, const uint popID) const {
      const LossConeSpeciesParameters& sP = speciesParams[popID];

      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;
      Realf vpara = vx;
      Realf vperp = sqrt(vy*vy + vz*vz);
      Realf mu    = vpara / sqrt(vpara*vpara + vperp*vperp);
      if (mu <= -0.5 or mu >= 0.5) {return 0.0;}
      else {return exp((- mass / (2.0 * kb)) * ((vx*vx) / sP.TEMPERATUREX + (vy*vy) / sP.TEMPERATUREY + (vz*vz) / sP.TEMPERATUREZ));}
   }

   Real LossCone::calcPhaseSpaceDensity(
      creal& x, creal& y, creal& z,
      creal& dx, creal& dy, creal& dz,
      creal& vx, creal& vy, creal& vz,
      creal& dvx, creal& dvy, creal& dvz,const uint popID
   ) const {
      const LossConeSpeciesParameters& sP = speciesParams[popID];
      const size_t meshID = getObjectWrapper().particleSpecies[popID].velocityMesh;
      vmesh::MeshParameters& meshParams = getObjectWrapper().velocityMeshes[meshID];
      if (vx < meshParams.meshMinLimits[0] + 0.5*dvx ||
          vy < meshParams.meshMinLimits[1] + 0.5*dvy ||
          vz < meshParams.meshMinLimits[2] + 0.5*dvz ||
          vx > meshParams.meshMaxLimits[0] - 1.5*dvx ||
          vy > meshParams.meshMaxLimits[1] - 1.5*dvy ||
          vz > meshParams.meshMaxLimits[2] - 1.5*dvz) {
         return 0.0;
      }
      
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;
      
      creal d_vx = dvx / (sP.nVelocitySamples-1);
      creal d_vy = dvy / (sP.nVelocitySamples-1);
      creal d_vz = dvz / (sP.nVelocitySamples-1);
      Real avg = 0.0;
      
      for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
         for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
            for (uint vk=0; vk<sP.nVelocitySamples; ++vk)
            {
               avg += getDistribValue(
                  vx+vi*d_vx - sP.velocityPertAbsAmp * (0.5 - rndVel[0] ),
                  vy+vj*d_vy - sP.velocityPertAbsAmp * (0.5 - rndVel[1] ),
                  vz+vk*d_vz - sP.velocityPertAbsAmp * (0.5 - rndVel[2] ), popID);
            }
      
      creal result = avg *
         sP.DENSITY * (1.0 + sP.densityPertRelAmp * (0.5 - rndRho)) *
         pow(mass / (2.0 * M_PI * kb ), 1.5) / ( sqrt(sP.TEMPERATUREX) * sqrt(sP.TEMPERATUREY) * sqrt(sP.TEMPERATUREZ) ) /
         (sP.nVelocitySamples*sP.nVelocitySamples*sP.nVelocitySamples );
      
      if(result < sP.maxwCutoff) {
         return 0.0;
      } else {
         return result;
      }
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
      
      setRandomSeed(cellID);
      
      this->rndRho=getRandomNumber();
      this->rndVel[0]=getRandomNumber();
      this->rndVel[1]=getRandomNumber();
      this->rndVel[2]=getRandomNumber();
   }

   void LossCone::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
      FsGrid< fsgrids::technical, 2>& technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(this->BX0,
                         this->BY0,
                         this->BZ0);

      setBackgroundField(bgField, BgBGrid);
      
      if(!P::isRestart) {
         const auto localSize = BgBGrid.getLocalSize().data();
         
         #pragma omp parallel for collapse(3)
         for (int x = 0; x < localSize[0]; ++x) {
            for (int y = 0; y < localSize[1]; ++y) {
               for (int z = 0; z < localSize[2]; ++z) {
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  const int64_t cellid = perBGrid.GlobalIDForCoords(x, y, z);
                  
                  setRandomSeed(cellid);
                  
                  cell->at(fsgrids::bfield::PERBX) = this->magXPertAbsAmp * (0.5 - getRandomNumber());
                  cell->at(fsgrids::bfield::PERBY) = this->magYPertAbsAmp * (0.5 - getRandomNumber());
                  cell->at(fsgrids::bfield::PERBZ) = this->magZPertAbsAmp * (0.5 - getRandomNumber());
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
