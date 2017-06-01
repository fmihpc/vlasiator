/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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

#include "Fluctuations.h"

using namespace spatial_cell;

Real projects::Fluctuations::rndRho, projects::Fluctuations::rndVel[3];


namespace projects {
   Fluctuations::Fluctuations(): TriAxisSearch() { }
   Fluctuations::~Fluctuations() { }
   bool Fluctuations::initialize(void) {return Project::initialize();}
   
   void Fluctuations::addParameters() {
      typedef Readparameters RP;
      RP::add("Fluctuations.BX0", "Background field value (T)", 1.0e-9);
      RP::add("Fluctuations.BY0", "Background field value (T)", 2.0e-9);
      RP::add("Fluctuations.BZ0", "Background field value (T)", 3.0e-9);
      RP::add("Fluctuations.magXPertAbsAmp", "Amplitude of the magnetic perturbation along x", 1.0e-9);
      RP::add("Fluctuations.magYPertAbsAmp", "Amplitude of the magnetic perturbation along y", 1.0e-9);
      RP::add("Fluctuations.magZPertAbsAmp", "Amplitude of the magnetic perturbation along z", 1.0e-9);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_Fluctuations.rho", "Number density (m^-3)", 1.0e7);
         RP::add(pop + "_Fluctuations.Temperature", "Temperature (K)", 2.0e6);
         RP::add(pop + "_Fluctuations.densityPertRelAmp", "Amplitude factor of the density perturbation", 0.1);
         RP::add(pop + "_Fluctuations.velocityPertAbsAmp", "Amplitude of the velocity perturbation", 1.0e6);
         RP::add(pop + "_Fluctuations.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
         RP::add(pop + "_Fluctuations.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
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
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         FluctuationsSpeciesParameters sP;
         RP::get(pop + "_Fluctuations.rho", sP.DENSITY);
         RP::get(pop + "_Fluctuations.Temperature", sP.TEMPERATURE);
         RP::get(pop + "_Fluctuations.densityPertRelAmp", sP.densityPertRelAmp);
         RP::get(pop + "_Fluctuations.velocityPertAbsAmp", sP.velocityPertAbsAmp);
         RP::get(pop + "_Fluctuations.nSpaceSamples", sP.nSpaceSamples);
         RP::get(pop + "_Fluctuations.nVelocitySamples", sP.nVelocitySamples);
         RP::get(pop + "_Fluctuations.maxwCutoff", sP.maxwCutoff);

         speciesParams.push_back(sP);
      }
   }
   
   Real Fluctuations::getDistribValue(creal& vx,creal& vy, creal& vz, const uint popID) const {
      const FluctuationsSpeciesParameters& sP = speciesParams[popID];

      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;
      return exp(- mass * (vx*vx + vy*vy + vz*vz) / (2.0 * kb * sP.TEMPERATURE));
   }

   Real Fluctuations::calcPhaseSpaceDensity(
      creal& x, creal& y, creal& z,
      creal& dx, creal& dy, creal& dz,
      creal& vx, creal& vy, creal& vz,
      creal& dvx, creal& dvy, creal& dvz,const uint popID
   ) const {
      const FluctuationsSpeciesParameters& sP = speciesParams[popID];
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
         pow(mass / (2.0 * M_PI * kb * sP.TEMPERATURE), 1.5) /
         (sP.nVelocitySamples*sP.nVelocitySamples*sP.nVelocitySamples);
      
      if(result < sP.maxwCutoff) {
         return 0.0;
      } else {
         return result;
      }
   }
   
   void Fluctuations::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
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
      
      setRandomSeed(cell,cellID);

      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      
      this->rndRho=getRandomNumber(cell);
      this->rndVel[0]=getRandomNumber(cell);
      this->rndVel[1]=getRandomNumber(cell);
      this->rndVel[2]=getRandomNumber(cell);
      
      cellParams[CellParams::PERBX] = this->magXPertAbsAmp * (0.5 - getRandomNumber(cell));
      cellParams[CellParams::PERBY] = this->magYPertAbsAmp * (0.5 - getRandomNumber(cell));
      cellParams[CellParams::PERBZ] = this->magZPertAbsAmp * (0.5 - getRandomNumber(cell));
   }

   void Fluctuations::setCellBackgroundField(SpatialCell* cell) const {
      ConstantField bgField;
      bgField.initialize(this->BX0,
                         this->BY0,
                         this->BZ0);
      
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
   }
   
   std::vector<std::array<Real, 3> > Fluctuations::getV0(
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
