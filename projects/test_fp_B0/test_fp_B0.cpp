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
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"

enum cases {BXCASE,BYCASE,BZCASE};

// SET THIS TO BXCASE,BYCASE,OR BZCASE TO SELECT ONE OF THE THREE CASES

static int CASE = BZCASE;

using namespace std;

bool initializeProject(void) {return true;}
bool addProjectParameters(){return true;}
bool getProjectParameters(){return true;}

bool cellParametersChanged(creal& t) {return false;}

Real sign(creal value)
{
   if(abs(value) < 1e-5) return 0.0;
   else return value / abs(value);
}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
   const Real T  = 1e-6;
   const Real n  = 1.0e7;
   creal ALPHA = 7.0 * M_PI / 4.0;

   Real VX,VY,VZ,ksi,eta;
   switch (CASE) {
      case BXCASE:
      ksi = ((y + 0.5 * dy)  * cos(ALPHA) + (z + 0.5 * dz) * sin(ALPHA)) / (2.0 * sqrt(2.0));
      eta = (-(y + 0.5 * dy)  * sin(ALPHA) + (z + 0.5 * dz) * cos(ALPHA)) / (2.0 * sqrt(2.0));
      VX = 0.0;
      VY = sign(cos(ALPHA)) * 0.5 + 0.1*cos(ALPHA) * sin(2.0 * M_PI * eta);
      VZ = sign(sin(ALPHA)) * 0.5 + 0.1*sin(ALPHA) * sin(2.0 * M_PI * eta); 
      break;
    case BYCASE:
      ksi = ((z + 0.5 * dz)  * cos(ALPHA) + (x + 0.5 * dx) * sin(ALPHA)) / (2.0 * sqrt(2.0));
      eta = (-(z + 0.5 * dz)  * sin(ALPHA) + (x + 0.5 * dx) * cos(ALPHA)) / (2.0 * sqrt(2.0));
      VX = sign(sin(ALPHA)) * 0.5 + 0.1*sin(ALPHA) * sin(2.0 * M_PI * eta);
      VY = 0.0;
      VZ = sign(cos(ALPHA)) * 0.5 + 0.1*cos(ALPHA) * sin(2.0 * M_PI * eta);
      break;
    case BZCASE:
      ksi = ((x + 0.5 * dx)  * cos(ALPHA) + (y + 0.5 * dy) * sin(ALPHA)) / (2.0 * sqrt(2.0));
      eta = (-(x + 0.5 * dx)  * sin(ALPHA) + (y + 0.5 * dy) * cos(ALPHA)) / (2.0 * sqrt(2.0));
      VX = sign(cos(ALPHA)) * 0.5 + 0.1*cos(ALPHA) * sin(2.0 * M_PI * eta);
      VY = sign(sin(ALPHA)) * 0.5 + 0.1*sin(ALPHA) * sin(2.0 * M_PI * eta);
      VZ = 0.0;
      break;
   }
   
   creal VX2 = (vx+0.5*dvx-VX)*(vx+0.5*dvx-VX);
   creal VY2 = (vy+0.5*dvy-VY)*(vy+0.5*dvy-VY);
   creal VZ2 = (vz+0.5*dvz-VZ)*(vz+0.5*dvz-VZ);
   
   creal CONST = physicalconstants::MASS_PROTON / 2.0 / physicalconstants::K_B / T;
   Real NORM = (physicalconstants::MASS_PROTON / 2.0 / M_PI / physicalconstants::K_B / T);
   NORM = n * pow(NORM,1.5);
   
   return NORM*exp(-CONST*(VX2+VY2+VZ2));
}

void calcBlockParameters(Real* blockParams) { }

void calcCellParameters(Real* cellParams,creal& t) {
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
   
   typedef Parameters P;
   creal x = cellParams[CellParams::XCRD] + 0.5 * cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD] + 0.5 * cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD] + 0.5 * cellParams[CellParams::DZ];

   creal B0 = 1e-10;
   creal B1 = 1e-9;
   
   switch (CASE) {
    case BXCASE:
      if (y >= -0.2 && y <= 0.2)
	if (z >= -0.2 && z <= 0.2)
	  cellParams[CellParams::BX] = B1;
	  cellParams[CellParams::BXFACEX0] = B0;
      break;
    case BYCASE:
      if (x >= -0.2 && x <= 0.2)
	if (z >= -0.2 && z <= 0.2)
	  cellParams[CellParams::BY] = B1;
	  cellParams[CellParams::BYFACEY0] = B0;
      break;
    case BZCASE:
      if (x >= -0.2 && x <= 0.2)
	if (y >= -0.2 && y <= 0.2)
	  cellParams[CellParams::BZ] = B1;
	  cellParams[CellParams::BZFACEZ0] = B0;
      break;
   }
}



void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}



