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

using namespace std;

typedef harrisParameters HP;
Real HP::SCA_LAMBDA = NAN;
Real HP::B0 = NAN;
Real HP::TEMPERATURE = NAN;
Real HP::DENSITY = NAN;

bool initializeProject(void) {
   return true;
}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("Harris.Scale_size", "Harris sheet scale size (m)", 150000.0);
   RP::add("Harris.B0", "Magnetic field at infinity (T)", 8.33061003094e-8);
   RP::add("Harris.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Harris.rho", "Number density at infinity (m^-3)", 1.0e7);
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("Harris.Scale_size", HP::SCA_LAMBDA);
   RP::get("Harris.B0", HP::B0);
   RP::get("Harris.Temperature", HP::TEMPERATURE);
   RP::get("Harris.rho", HP::DENSITY);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
  
  creal mass = 1.67262171e-27; // m_p in kg
  creal k = 1.3806505e-23; // Boltzmann
  creal mu0 = 1.25663706144e-6; // mu_0
  creal q = 1.60217653e-19; // q_i
  
  creal Vy0 = 0.0;
  
  return HP::DENSITY * pow(mass / (2.0 * M_PI * k * HP::TEMPERATURE), 1.5) * (
     5.0 / pow(cosh((x + 0.5 * dx) / (HP::SCA_LAMBDA)), 2.0) * exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy - Vy0, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * k * HP::TEMPERATURE))
    +
    exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * k * HP::TEMPERATURE)));
}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
  creal x = cellParams[CellParams::XCRD];
  creal dx = cellParams[CellParams::DX];
  
  cellParams[CellParams::EX   ] = 0.0;
  cellParams[CellParams::EY   ] = 0.0;
  cellParams[CellParams::EZ   ] = 0.0;
  cellParams[CellParams::BX   ] = 0.0;
  cellParams[CellParams::BY   ] = 0.0;
  cellParams[CellParams::BZ   ] = HP::B0 * tanh((x + 0.5 * dx) / HP::SCA_LAMBDA);
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

