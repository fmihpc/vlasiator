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
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "cell_spatial.h"
#include "common.h"
#include "project.h"
#include "parameters.h"

using namespace std;

#ifndef PARGRID
bool initializeProject(dccrg<SpatialCell>& mpiGrid) {
#else
bool initializeProject(ParGrid<SpatialCell>& mpiGrid) {
#endif
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

/*Real calcPhaseSpaceDensity(creal& z,creal& x,creal& y,creal& dz,creal& dx,creal& dy,
			   creal& vz,creal& vx,creal& vy,creal& dvz,creal& dvx,creal& dvy) {*/
Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
  
  creal mass = 1.67262171e-27; // m_p in kg
  creal k = 1.3806505e-23; // Boltzmann
//  creal mu0 = 1.25663706144e-6; // mu_0
//  creal q = 1.60217653e-19; // q_i

  creal ksi = ((x + 0.5 * dx)  * cos(ALPHA) + (y + 0.5 * dy) * sin(ALPHA)) / WAVELENGTH;
  creal Vx = A_VEL * ALFVEN_VEL * sin(ALPHA) * sin(2.0 * M_PI * ksi);
  creal Vy = - A_VEL * ALFVEN_VEL * cos(ALPHA) * sin(2.0 * M_PI * ksi);
  creal Vz = - A_VEL * ALFVEN_VEL * cos(2.0 * M_PI * ksi);
  
  creal den = DENSITY * pow(mass / (2.0 * M_PI * k * TEMPERATURE), 1.5) *
    exp(- mass * (pow(vx + 0.5 * dvx - Vx, 2.0) + pow(vy + 0.5 * dvy - Vy, 2.0) + pow(vz + 0.5 * dvz - Vz, 2.0)) / (2.0 * k * TEMPERATURE));
  return den;
}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
  creal x = cellParams[CellParams::XCRD];
  creal dx = cellParams[CellParams::DX];
  creal y = cellParams[CellParams::YCRD];
  creal dy = cellParams[CellParams::DY];
  creal ksi = ((x + 0.5 * dx)  * cos(ALPHA) + (y + 0.5 * dy) * sin(ALPHA)) / WAVELENGTH;
  
  cellParams[CellParams::EX   ] = 0.0;
  cellParams[CellParams::EY   ] = 0.0;
  cellParams[CellParams::EZ   ] = 0.0;
  cellParams[CellParams::BX   ] = B0 * cos(ALPHA) - A_MAG * B0 * sin(ALPHA) * sin(2.0 * M_PI * ksi);
  cellParams[CellParams::BY   ] = B0 * sin(ALPHA) + A_MAG * B0 * cos(ALPHA) * sin(2.0 * M_PI * ksi);
  cellParams[CellParams::BZ   ] = B0 * A_MAG * cos(2.0 * M_PI * ksi);
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
#ifndef PARGRID
void calcSimParameters(dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
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

