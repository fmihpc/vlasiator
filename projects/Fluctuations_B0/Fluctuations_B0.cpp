/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"

typedef fluctuationsParameters FlucP;
Real FlucP::BX0 = NAN;
Real FlucP::BY0 = NAN;
Real FlucP::BZ0 = NAN;
Real FlucP::DENSITY = NAN;
Real FlucP::TEMPERATURE = NAN;
Real FlucP::magPertAmp = NAN;
Real FlucP::densityPertAmp = NAN;
Real FlucP::velocityPertAmp = NAN;
Real FlucP::maxwCutoff = NAN;
uint FlucP::nSpaceSamples = 0;
uint FlucP::nVelocitySamples = 0;

bool initializeProject(void) { return true; }

bool addProjectParameters() {
   typedef Readparameters RP;
   RP::add("Fluctuations.BX0", "Background field value (T)", 1.0e-9);
   RP::add("Fluctuations.BY0", "Background field value (T)", 2.0e-9);
   RP::add("Fluctuations.BZ0", "Background field value (T)", 3.0e-9);
   RP::add("Fluctuations.rho", "Number density (m^-3)", 1.0e7);
   RP::add("Fluctuations.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Fluctuations.magPertAmp", "Amplitude of the magnetic perturbation", 1.0e-9);
   RP::add("Fluctuations.densityPertAmp", "Amplitude factor of the density perturbation", 0.1);
   RP::add("Fluctuations.velocityPertAmp", "Amplitude of the velocity perturbation", 1.0e6);
   RP::add("Fluctuations.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Fluctuations.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   RP::add("Fluctuations.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
   return true;
}

bool getProjectParameters() {
   typedef Readparameters RP;
   RP::get("Fluctuations.BX0", FlucP::BX0);
   RP::get("Fluctuations.BY0", FlucP::BY0);
   RP::get("Fluctuations.BZ0", FlucP::BZ0);
   RP::get("Fluctuations.rho", FlucP::DENSITY);
   RP::get("Fluctuations.Temperature", FlucP::TEMPERATURE);
   RP::get("Fluctuations.magPertAmp", FlucP::magPertAmp);
   RP::get("Fluctuations.densityPertAmp", FlucP::densityPertAmp);
   RP::get("Fluctuations.velocityPertAmp", FlucP::velocityPertAmp);
   RP::get("Fluctuations.nSpaceSamples", FlucP::nSpaceSamples);
   RP::get("Fluctuations.nVelocitySamples", FlucP::nVelocitySamples);
   RP::get("Fluctuations.maxwCutoff", FlucP::maxwCutoff);
   return true;
}


bool cellParametersChanged(creal& t) {return false;}


Real getDistribValue(creal& vx,creal& vy, creal& vz) {
   creal k = 1.3806505e-23; // Boltzmann
   creal mass = 1.67262171e-27; // m_p in kg
   return exp(- mass * (vx*vx + vy*vy + vz*vz) / (2.0 * k * FlucP::TEMPERATURE));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
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
   
   static int rndRho = 0;
   static int rndVel[3] = {0};
   int cellID = (int) (x / dx) +
                (int) (y / dy) * Parameters::xcells_ini +
                (int) (z / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
   srand(cellID);
   rndRho = rand();
   rndVel[0] =rand();
   rndVel[1] =rand();
   rndVel[2] =rand();
   
   creal d_vx = dvx / (FlucP::nVelocitySamples-1);
   creal d_vy = dvy / (FlucP::nVelocitySamples-1);
   creal d_vz = dvz / (FlucP::nVelocitySamples-1);
   Real avg = 0.0;
   
   for (uint vi=0; vi<FlucP::nVelocitySamples; ++vi)
      for (uint vj=0; vj<FlucP::nVelocitySamples; ++vj)
	 for (uint vk=0; vk<FlucP::nVelocitySamples; ++vk)
         {
	    avg += getDistribValue(
	       vx+vi*d_vx - FlucP::velocityPertAmp * (0.5 - (double)rndVel[0] / (double)RAND_MAX),
	       vy+vj*d_vy - FlucP::velocityPertAmp * (0.5 - (double)rndVel[1] / (double)RAND_MAX),
	       vz+vk*d_vz - FlucP::velocityPertAmp * (0.5 - (double)rndVel[2] / (double)RAND_MAX));
         }
   
   creal result = avg *
	          FlucP::DENSITY * (1.0 + FlucP::densityPertAmp * (0.5 - (double)rndRho / (double)RAND_MAX)) *
	          pow(mass / (2.0 * M_PI * k * FlucP::TEMPERATURE), 1.5) /
//   	            (Parameters::vzmax - Parameters::vzmin) / 
                  (FlucP::nVelocitySamples*FlucP::nVelocitySamples*FlucP::nVelocitySamples);
   if(result < FlucP::maxwCutoff) {
      return 0.0;
   } else {
      return result;
   }
}
      
void calcBlockParameters(Real* blockParams) {}

void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD];
   creal dy = cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD];
   creal dz = cellParams[CellParams::DZ];
   
   int cellID = (int) (x / dx) +
   (int) (y / dy) * Parameters::xcells_ini +
   (int) (z / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
   srand(cellID);

   // default initial values
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = 0.0;
   cellParams[CellParams::BY   ] = 0.0;
   cellParams[CellParams::BZ   ] = 0.0;
   cellParams[CellParams::BXFACEX0] = 0;
   cellParams[CellParams::BYFACEX0] = 0;
   cellParams[CellParams::BZFACEX0] = 0;
   cellParams[CellParams::BXFACEY0] = 0;
   cellParams[CellParams::BYFACEY0] = 0;
   cellParams[CellParams::BZFACEY0] = 0;
   cellParams[CellParams::BXFACEZ0] = 0;
   cellParams[CellParams::BYFACEZ0] = 0;
   cellParams[CellParams::BZFACEZ0] = 0;

   // user supplied initial values
   cellParams[CellParams::BXFACEX0] = FlucP::BX0;
   cellParams[CellParams::BYFACEY0] = FlucP::BY0;
   cellParams[CellParams::BZFACEZ0] = FlucP::BZ0;
   cellParams[CellParams::BY   ] = FlucP::magPertAmp * (0.5 - (double)rand() / (double)RAND_MAX);
   cellParams[CellParams::BZ   ] = FlucP::magPertAmp * (0.5 - (double)rand() / (double)RAND_MAX);
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   const std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

