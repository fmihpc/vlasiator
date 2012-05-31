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

using namespace std;

typedef test_fpParameters tfP;
Real tfP::B0 = NAN;
Real tfP::DENSITY = NAN;
Real tfP::TEMPERATURE = NAN;
Real tfP::ALPHA = NAN;
int  tfP::CASE = 5;
bool tfP::shear = false;

bool initializeProject(void) {
   tfP::ALPHA *= M_PI / 4.0;
   return true;
}

bool addProjectParameters(void){
   typedef Readparameters RP;
   RP::add("test_fp.B0", "Background field value (T)", 1.0e-9);
   RP::add("test_fp.rho", "Number density (m^-3)", 1.0e7);
   RP::add("test_fp.Temperature", "Temperature (K)", 1.0e-6);
   RP::add("test_fp.angle", "Orientation of the propagation expressed in pi/4", 0.0);
   RP::add("test_fp.Bdirection", "Direction of the magnetic field (0:x, 1:y, 2:z)", 0);
   RP::add("test_fp.shear", "Add a shear (if false, V=0.5 everywhere).", true);
   return true;
}

bool getProjectParameters(void){
   typedef Readparameters RP;
   RP::get("test_fp.B0", tfP::B0);
   RP::get("test_fp.rho", tfP::DENSITY);
   RP::get("test_fp.Temperature", tfP::TEMPERATURE);
   RP::get("test_fp.angle", tfP::ALPHA);
   RP::get("test_fp.Bdirection", tfP::CASE);
   RP::get("test_fp.shear", tfP::shear);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

Real sign(creal value)
{
   if(abs(value) < 1e-5) return 0.0;
   else return value / abs(value);
}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
   Real VX,VY,VZ;
   if (tfP::shear == true)
   {
      Real ksi,eta;
      switch (tfP::CASE) {
	 case BXCASE:
	 ksi = ((y + 0.5 * dy)  * cos(tfP::ALPHA) + (z + 0.5 * dz) * sin(tfP::ALPHA)) / (2.0 * sqrt(2.0));
	 eta = (-(y + 0.5 * dy)  * sin(tfP::ALPHA) + (z + 0.5 * dz) * cos(tfP::ALPHA)) / (2.0 * sqrt(2.0));
	 VX = 0.0;
	 VY = sign(cos(tfP::ALPHA)) * 0.5 + 0.1*cos(tfP::ALPHA) * sin(2.0 * M_PI * eta);
	 VZ = sign(sin(tfP::ALPHA)) * 0.5 + 0.1*sin(tfP::ALPHA) * sin(2.0 * M_PI * eta); 
	 break;
      case BYCASE:
	 ksi = ((z + 0.5 * dz)  * cos(tfP::ALPHA) + (x + 0.5 * dx) * sin(tfP::ALPHA)) / (2.0 * sqrt(2.0));
	 eta = (-(z + 0.5 * dz)  * sin(tfP::ALPHA) + (x + 0.5 * dx) * cos(tfP::ALPHA)) / (2.0 * sqrt(2.0));
	 VX = sign(sin(tfP::ALPHA)) * 0.5 + 0.1*sin(tfP::ALPHA) * sin(2.0 * M_PI * eta);
	 VY = 0.0;
	 VZ = sign(cos(tfP::ALPHA)) * 0.5 + 0.1*cos(tfP::ALPHA) * sin(2.0 * M_PI * eta);
	 break;
      case BZCASE:
	 ksi = ((x + 0.5 * dx)  * cos(tfP::ALPHA) + (y + 0.5 * dy) * sin(tfP::ALPHA)) / (2.0 * sqrt(2.0));
	 eta = (-(x + 0.5 * dx)  * sin(tfP::ALPHA) + (y + 0.5 * dy) * cos(tfP::ALPHA)) / (2.0 * sqrt(2.0));
	 VX = sign(cos(tfP::ALPHA)) * 0.5 + 0.1*cos(tfP::ALPHA) * sin(2.0 * M_PI * eta);
	 VY = sign(sin(tfP::ALPHA)) * 0.5 + 0.1*sin(tfP::ALPHA) * sin(2.0 * M_PI * eta);
	 VZ = 0.0;
	 break;
      }
   } else {
      switch (tfP::CASE) {
	 case BXCASE:
	    VX = 0.0;
	    VY = cos(tfP::ALPHA) * 0.5;
	    VZ = sin(tfP::ALPHA) * 0.5; 
	    break;
	 case BYCASE:
	    VX = sin(tfP::ALPHA) * 0.5;
	    VY = 0.0;
	    VZ = cos(tfP::ALPHA) * 0.5;
	    break;
	 case BZCASE:
	    VX = cos(tfP::ALPHA) * 0.5;
	    VY = sin(tfP::ALPHA) * 0.5;
	    VZ = 0.0;
	    break;
      }
   }
   
   creal VX2 = (vx+0.5*dvx-VX)*(vx+0.5*dvx-VX);
   creal VY2 = (vy+0.5*dvy-VY)*(vy+0.5*dvy-VY);
   creal VZ2 = (vz+0.5*dvz-VZ)*(vz+0.5*dvz-VZ);
   
   creal CONST = physicalconstants::MASS_PROTON / 2.0 / physicalconstants::K_B / tfP::TEMPERATURE;
   Real NORM = (physicalconstants::MASS_PROTON / 2.0 / M_PI / physicalconstants::K_B / tfP::TEMPERATURE);
   NORM = tfP::DENSITY * pow(NORM,1.5);
   
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
   
   typedef Parameters P;
   creal x = cellParams[CellParams::XCRD] + 0.5 * cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD] + 0.5 * cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD] + 0.5 * cellParams[CellParams::DZ];
   
   switch (tfP::CASE) {
    case BXCASE:
      if (y >= -0.2 && y <= 0.2)
	if (z >= -0.2 && z <= 0.2)
	   cellParams[CellParams::BX] = tfP::B0;
      break;
    case BYCASE:
      if (x >= -0.2 && x <= 0.2)
	if (z >= -0.2 && z <= 0.2)
	   cellParams[CellParams::BY] = tfP::B0;
      break;
    case BZCASE:
      if (x >= -0.2 && x <= 0.2)
	if (y >= -0.2 && y <= 0.2)
	   cellParams[CellParams::BZ] = tfP::B0;
      break;
   }
}



void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}



