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
along with this program. If not, see <http://www.gnu.org/licenses/>.
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

typedef dispersionParameters DispP;
Real DispP::BX0 = NAN;
Real DispP::BY0 = NAN;
Real DispP::BZ0 = NAN;
Real DispP::DENSITY = NAN;
Real DispP::TEMPERATURE = NAN;
Real DispP::magPertAmp = NAN;
Real DispP::densityPertAmp = NAN;
Real DispP::velocityPertAmp = NAN;
uint DispP::seed = 0;
uint DispP::sectorSize = 0;
uint DispP::nSpaceSamples = 0;
uint DispP::nVelocitySamples = 0;

bool initializeProject(void) {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   srand(DispP::seed*rank);
   return true;
}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("Dispersion.BX0", "Background field value (T)", 1.0e-9);
   RP::add("Dispersion.BY0", "Background field value (T)", 2.0e-9);
   RP::add("Dispersion.BZ0", "Background field value (T)", 3.0e-9);
   RP::add("Dispersion.rho", "Number density (m^-3)", 1.0e7);
   RP::add("Dispersion.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Dispersion.magPertAmp", "Amplitude of the magnetic perturbation", 1.0e-9);
   RP::add("Dispersion.densityPertAmp", "Amplitude factor of the density perturbation", 0.1);
   RP::add("Dispersion.velocityPertAmp", "Amplitude of the velocity perturbation", 1.0e6);
   RP::add("Dispersion.seed","Seed integer for the srand() function multiplied by rank", 42);
   RP::add("Dispersion.sectorSize", "Maximal size of the sectors for randomising", 10);
   RP::add("Dispersion.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Dispersion.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("Dispersion.BX0", DispP::BX0);
   RP::get("Dispersion.BY0", DispP::BY0);
   RP::get("Dispersion.BZ0", DispP::BZ0);
   RP::get("Dispersion.rho", DispP::DENSITY);
   RP::get("Dispersion.Temperature", DispP::TEMPERATURE);
   RP::get("Dispersion.magPertAmp", DispP::magPertAmp);
   RP::get("Dispersion.densityPertAmp", DispP::densityPertAmp);
   RP::get("Dispersion.velocityPertAmp", DispP::velocityPertAmp);
   RP::get("Dispersion.seed", DispP::seed);
   RP::get("Dispersion.sectorSize", DispP::sectorSize);
   RP::get("Dispersion.nSpaceSamples", DispP::nSpaceSamples);
   RP::get("Dispersion.nVelocitySamples", DispP::nVelocitySamples);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

Real getDistribValue(creal& vx,creal& vy, creal& vz) {
   creal k = 1.3806505e-23; // Boltzmann
   creal mass = 1.67262171e-27; // m_p in kg
   return exp(- mass * (vx*vx + vy*vy + vz*vz) / (2.0 * k * DispP::TEMPERATURE));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal mass = 1.67262171e-27; // m_p in kg
   creal k = 1.3806505e-23; // Boltzmann
   creal mu0 = 1.25663706144e-6; // mu_0
   creal q = 1.60217653e-19; // q_i
   
   static int spaceIndexOld[3] = {std::numeric_limits<int>::min(),
                                  std::numeric_limits<int>::min(),
                                  std::numeric_limits<int>::min()};
   static int spaceIndex[3] = {0};
   static int rndRho = 0;
   static int rndVel[3] = {0};
#pragma omp critical
   {
      //critical region since srand not thread-safe. A nicer fix would be to use a thread-safe rng
      //e.g., http://linux.die.net/man/3/random_r (FIXME)
      static int rndRhoSector = rand()%DispP::sectorSize+1;
      static int rndVelSector = rand()%DispP::sectorSize+1;
      static int cptRhoSector = 0;
      static int cptVelSector = 0;
      
      //static variables should be threadprivate
#pragma omp threadprivate(spaceIndexOld,spaceIndex,rndRho,rndVel,rndRhoSector,rndVelSector,cptRhoSector,cptVelSector)
   
      spaceIndex[0] = (int) (x / dx);
      spaceIndex[1] = (int) (y / dy);
      spaceIndex[2] = (int) (z / dz);
      if(spaceIndex[0] != spaceIndexOld[0] ||
	 spaceIndex[1] != spaceIndexOld[1] ||
	 spaceIndex[2] != spaceIndexOld[2]) {
	 if(cptRhoSector++%rndRhoSector == 0)
	 {
	    rndRho = rand();
	    rndRhoSector = rand()%DispP::sectorSize+1;
	 }
	 if(cptVelSector++%rndVelSector == 0)
	 {
	    rndVel = {rand(), rand(), rand()};
	    rndVelSector = rand()%DispP::sectorSize+1;
	 }
      }
   } // end of omp critical region
   spaceIndexOld[0] = spaceIndex[0];
   spaceIndexOld[1] = spaceIndex[1];
   spaceIndexOld[2] = spaceIndex[2];
   
   creal d_vx = dvx / (DispP::nVelocitySamples-1);
   creal d_vy = dvy / (DispP::nVelocitySamples-1);
   creal d_vz = dvz / (DispP::nVelocitySamples-1);
   Real avg = 0.0;
   for (uint vi=0; vi<DispP::nVelocitySamples; ++vi)
      for (uint vj=0; vj<DispP::nVelocitySamples; ++vj)
	 for (uint vk=0; vk<DispP::nVelocitySamples; ++vk)
         {
	    avg += getDistribValue(
	       vx+vi*d_vx - DispP::velocityPertAmp * (0.5 - (double)rndVel[0] / (double)RAND_MAX),
	       vy+vj*d_vy - DispP::velocityPertAmp * (0.5 - (double)rndVel[1] / (double)RAND_MAX),
	       vz+vk*d_vz - DispP::velocityPertAmp * (0.5 - (double)rndVel[2] / (double)RAND_MAX));
         }
   
   return avg *
   DispP::DENSITY * (1.0 + DispP::densityPertAmp * (0.5 - (double)rndRho / (double)RAND_MAX)) *
   pow(mass / (2.0 * M_PI * k * DispP::TEMPERATURE), 1.5) /
//   (Parameters::vzmax - Parameters::vzmin) / 
(DispP::nVelocitySamples*DispP::nVelocitySamples*DispP::nVelocitySamples);
}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   //creal x = cellParams[CellParams::XCRD];
   //creal dx = cellParams[CellParams::DX];
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = DispP::BX0;
#pragma omp critical
   {
      //critical region since srand not thread-safe. A nicer fix would be to use a thread-safe rng
      //e.g., http://linux.die.net/man/3/random_r (FIXME)
      cellParams[CellParams::BY   ] = DispP::magPertAmp * (0.5 - (double)rand() / (double)RAND_MAX);
      cellParams[CellParams::BZ   ] = DispP::magPertAmp * (0.5 - (double)rand() / (double)RAND_MAX);
   } // end of omp critical region
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   const std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

