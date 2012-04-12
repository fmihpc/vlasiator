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

typedef larmorParameters LarmP;
Real LarmP::BX0 = NAN;
Real LarmP::BY0 = NAN;
Real LarmP::BZ0 = NAN;
Real LarmP::VX0 = NAN;
Real LarmP::VY0 = NAN;
Real LarmP::VZ0 = NAN;
Real LarmP::X0 = NAN;
Real LarmP::Y0 = NAN;
Real LarmP::Z0 = NAN;
Real LarmP::DENSITY = NAN;


bool initializeProject(void) {return true;}

bool addProjectParameters() {
   typedef Readparameters RP;
   RP::add("Larmor.BX0", "Background field value (T)", 0.0);
   RP::add("Larmor.BY0", "Background field value (T)", 0.0);
   RP::add("Larmor.BZ0", "Background field value (T)", 0.0);
   RP::add("Larmor.VX0", "Bulk velocity in x", 0.0);
   RP::add("Larmor.VY0", "Bulk velocity in y", 0.0);
   RP::add("Larmor.VZ0", "Bulk velocity in z", 0.0);
   RP::add("Larmor.X0", "Initial Position", 0.0);
   RP::add("Larmor.Y0", "Initial Position", 0.0);
   RP::add("Larmor.Z0", "Initial Position", 0.0);
   RP::add("Larmor.rho", "Number density (m^-3)", 1.0e7);
   return true;
}

bool getProjectParameters() {
   typedef Readparameters RP;
   RP::get("Larmor.BX0", LarmP::BX0);
   RP::get("Larmor.BY0", LarmP::BY0);
   RP::get("Larmor.BZ0", LarmP::BZ0);
   RP::get("Larmor.VX0", LarmP::VX0);
   RP::get("Larmor.VY0", LarmP::VY0);
   RP::get("Larmor.VZ0", LarmP::VZ0);
   RP::get("Larmor.X0", LarmP::X0);
   RP::get("Larmor.Y0", LarmP::Y0);
   RP::get("Larmor.Z0", LarmP::Z0);
   RP::get("Larmor.rho", LarmP::DENSITY);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}




Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {

   static bool isSet=false;

   if(vx < Parameters::vxmin + 0.5 * dvx ||
      vy < Parameters::vymin + 0.5 * dvy ||
      vz < Parameters::vzmin + 0.5 * dvz ||
      vx > Parameters::vxmax - 1.5 * dvx ||
      vy > Parameters::vymax - 1.5 * dvy ||
      vz > Parameters::vzmax - 1.5 * dvz
   ) return 0.0;

   if(isSet)
      return 0.0; //exactly one value to be set

   
   creal mass = Parameters::m;
   creal q = Parameters::q;
   creal k = 1.3806505e-23; // Boltzmann
   creal mu0 = 1.25663706144e-6; // mu_0

   if( fabs(vx-LarmP::VX0)<dvx &&
       fabs(vy-LarmP::VY0)<dvy &&
       fabs(vz-LarmP::VZ0)<dvz &&
       fabs(x-LarmP::X0)<dx &&
       fabs(y-LarmP::Y0)<dy &&
       fabs(z-LarmP::Z0)<dz){
      isSet=true;
      return LarmP::DENSITY/(dvx*dvy*dvz);
   }

   return 0.0;
   
}
      
void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

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
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   
   cellParams[CellParams::BX   ] = LarmP::BX0;
   cellParams[CellParams::BY   ] = LarmP::BY0;
   cellParams[CellParams::BZ   ] = LarmP::BZ0;
   
   cellParams[CellParams::BXVOL   ] = LarmP::BX0;
   cellParams[CellParams::BYVOL   ] = LarmP::BY0;
   cellParams[CellParams::BZVOL   ] = LarmP::BZ0;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   const std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

