/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2011, 2012 Finnish Meteorological Institute
 * 
 * Vlasiator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 * 
 * Vlasiator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"

#include "Dispersion.h"

#ifndef _AIX
int32_t projects::Dispersion::rndRho, projects::Dispersion::rndVel[3];
#else
int64_t projects::Dispersion::rndRho, projects::Dispersion::rndVel[3];
#endif

namespace projects {
   Dispersion::Dispersion(): Project() { }
   Dispersion::~Dispersion() { }
   
   bool Dispersion::initialize(void) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      
      memset(&(this->rngDataBuffer), 0, sizeof(this->rngDataBuffer));
      #ifndef _AIX
      initstate_r(this->seed*myRank, &(this->rngStateBuffer[0]), 256, &(this->rngDataBuffer));
      #else
      initstate_r(this->seed*myRank, &(this->rngStateBuffer[0]), 256, NULL, &(this->rngDataBuffer));
      #endif
      return true;
   }
   
   void Dispersion::addParameters() {
      typedef Readparameters RP;
      RP::add("Dispersion.B0", "Guide magnetic field strength (T)", 1.0e-9);
      RP::add("Dispersion.angleXY", "Orientation of the guide magnetic field with respect to the x-axis in x-y plane (rad)", 0.001);
      RP::add("Dispersion.angleXZ", "Orientation of the guide magnetic field with respect to the x-axis in x-z plane (rad)", 0.001);
      RP::add("Dispersion.rho", "Number density (m^-3)", 1.0e7);
      RP::add("Dispersion.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Dispersion.magXPertAbsAmp", "Absolute amplitude of the magnetic perturbation along x (T)", 1.0e-9);
      RP::add("Dispersion.magYPertAbsAmp", "Absolute amplitude of the magnetic perturbation along y (T)", 1.0e-9);
      RP::add("Dispersion.magZPertAbsAmp", "Absolute amplitude of the magnetic perturbation along z (T)", 1.0e-9);
      RP::add("Dispersion.densityPertRelAmp", "Relative amplitude of the density perturbation", 0.1);
      RP::add("Dispersion.velocityPertAbsAmp", "Absolute amplitude of the velocity perturbation", 1.0e6);
      RP::add("Dispersion.seed", "Seed for the RNG", 42);
      RP::add("Dispersion.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Dispersion.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
      RP::add("Dispersion.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
   }
   
   void Dispersion::getParameters() {
      typedef Readparameters RP;
      RP::get("Dispersion.B0", this->B0);
      RP::get("Dispersion.angleXY", this->angleXY);
      RP::get("Dispersion.angleXZ", this->angleXZ);
      RP::get("Dispersion.rho", this->DENSITY);
      RP::get("Dispersion.Temperature", this->TEMPERATURE);
      RP::get("Dispersion.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("Dispersion.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("Dispersion.magZPertAbsAmp", this->magZPertAbsAmp);
      RP::get("Dispersion.densityPertRelAmp", this->densityPertRelAmp);
      RP::get("Dispersion.velocityPertAbsAmp", this->velocityPertAbsAmp);
      RP::get("Dispersion.seed", this->seed);
      RP::get("Dispersion.nSpaceSamples", this->nSpaceSamples);
      RP::get("Dispersion.nVelocitySamples", this->nVelocitySamples);
      RP::get("Dispersion.maxwCutoff", this->maxwCutoff);
   }
   
   Real Dispersion::getDistribValue(creal& vx,creal& vy, creal& vz) {
      creal k = 1.3806505e-23; // Boltzmann
      creal mass = 1.67262171e-27; // m_p in kg
      return exp(- mass * (vx*vx + vy*vy + vz*vz) / (2.0 * k * this->TEMPERATURE));
   }
   
   Real Dispersion::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      if(vx < Parameters::vxmin + 0.5 * dvx ||
         vy < Parameters::vymin + 0.5 * dvy ||
         vz < Parameters::vzmin + 0.5 * dvz ||
         vx > Parameters::vxmax - 1.5 * dvx ||
         vy > Parameters::vymax - 1.5 * dvy ||
         vz > Parameters::vzmax - 1.5 * dvz
      ) return 0.0;
      
      creal mass = Parameters::m;
      creal q = Parameters::q;
      creal k = 1.3806505e-23; // Boltzmann
      creal mu0 = 1.25663706144e-6; // mu_0
      
      creal d_vx = dvx / (this->nVelocitySamples-1);
      creal d_vy = dvy / (this->nVelocitySamples-1);
      creal d_vz = dvz / (this->nVelocitySamples-1);
      Real avg = 0.0;
      
      for (uint vi=0; vi<this->nVelocitySamples; ++vi)
         for (uint vj=0; vj<this->nVelocitySamples; ++vj)
            for (uint vk=0; vk<this->nVelocitySamples; ++vk)
            {
               avg += getDistribValue(
                  vx+vi*d_vx - this->velocityPertAbsAmp * (0.5 - (double)(this->rndVel[0]) / (double)RAND_MAX),
                  vy+vj*d_vy - this->velocityPertAbsAmp * (0.5 - (double)(this->rndVel[1]) / (double)RAND_MAX),
                  vz+vk*d_vz - this->velocityPertAbsAmp * (0.5 - (double)(this->rndVel[2]) / (double)RAND_MAX)
               );
            }
            
            creal result = avg *
            this->DENSITY * (1.0 + this->densityPertRelAmp * (0.5 - (double)(Dispersion::rndRho) / (double)RAND_MAX)) *
            pow(mass / (2.0 * M_PI * k * this->TEMPERATURE), 1.5) /
            //            (Parameters::vzmax - Parameters::vzmin) / 
            (this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
            if(result < this->maxwCutoff) {
               return 0.0;
            } else {
               return result;
            }
   }
   
   void Dispersion::calcCellParameters(Real* cellParams,creal& t) {
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      creal z = cellParams[CellParams::ZCRD];
      creal dz = cellParams[CellParams::DZ];
      
      uint cellID = (int) ((x - Parameters::xmin) / dx) +
      (int) ((y - Parameters::ymin) / dy) * Parameters::xcells_ini +
      (int) ((z - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
      
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      
      #ifndef _AIX
      int32_t rndBuffer[3];
      random_r(&rngDataBuffer, &rndBuffer[0]);
      random_r(&rngDataBuffer, &rndBuffer[1]);
      random_r(&rngDataBuffer, &rndBuffer[2]);
      random_r(&rngDataBuffer, &(this->rndRho));
      random_r(&rngDataBuffer, &(this->rndVel[0]));
      random_r(&rngDataBuffer, &(this->rndVel[1]));
      random_r(&rngDataBuffer, &(this->rndVel[2]));
      #else
      int64_t rndBuffer[3];
      random_r(&rndBuffer[0], &rngDataBuffer);
      random_r(&rndBuffer[1], &rngDataBuffer);
      random_r(&rndBuffer[2], &rngDataBuffer);
      random_r(&(this->rndRho), &rngDataBuffer);
      random_r(&(this->rndVel[0]), &rngDataBuffer);
      random_r(&(this->rndVel[1]), &rngDataBuffer);
      random_r(&(this->rndVel[2]), &rngDataBuffer);
      #endif

      cellParams[CellParams::PERBX] = this->magXPertAbsAmp * (0.5 - (double)rndBuffer[0] / (double)RAND_MAX);
      cellParams[CellParams::PERBY] = this->magYPertAbsAmp * (0.5 - (double)rndBuffer[1] / (double)RAND_MAX);
      cellParams[CellParams::PERBZ] = this->magZPertAbsAmp * (0.5 - (double)rndBuffer[2] / (double)RAND_MAX);
      
      cellParams[CellParams::BGBX] = this->B0 * cos(this->angleXY) * cos(this->angleXZ);
      cellParams[CellParams::BGBY] = this->B0 * sin(this->angleXY) * cos(this->angleXZ);
      cellParams[CellParams::BGBZ] = this->B0 * sin(this->angleXZ);
   }
} // namespace projects
