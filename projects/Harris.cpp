/*
This file is part of Vlasiator.

Copyright 2011 Finnish Meteorological Institute

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

using namespace std;

bool initializeProject(void) {
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

/*Real calcPhaseSpaceDensity(creal& z,creal& x,creal& y,creal& dz,creal& dx,creal& dy,
			   creal& vz,creal& vx,creal& vy,creal& dvz,creal& dvx,creal& dvy) {*/
Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
  
  creal mass = 1.67262171e-27; // m_p in kg
  creal k = 1.3806505e-23; // Boltzmann
  creal mu0 = 1.25663706144e-6; // mu_0
  creal q = 1.60217653e-19; // q_i
  
  creal Vy0 = 0.0;
  
  return DENSITY * pow(mass / (2.0 * M_PI * k * TEMPERATURE), 1.5) * (
    5.0 / pow(cosh((x + 0.5 * dx) / (SCA_LAMBDA)), 2.0) * exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy - Vy0, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * k * TEMPERATURE))
    +
    exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * k * TEMPERATURE)));
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
  cellParams[CellParams::BZ   ] = B0 * tanh((x + 0.5 * dx) / SCA_LAMBDA);
}

/*void calcCellParameters(Real* cellParams,creal& t) {
   creal y = cellParams[CellParams::YCRD];
   creal dy = cellParams[CellParams::DY];
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BY   ] = 0.0;
   cellParams[CellParams::BZ   ] = 0.0;
   cellParams[CellParams::BX   ] = B0 * tanh((y + 0.5 * dy) / SCA_LAMBDA);
}*/

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
#ifndef PARGRID
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
#else
void calcSimParameters(ParGrid<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<ID::type> cells;
   mpiGrid.getCells(cells);
#endif
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->cpu_cellParams, t);
   }
}

