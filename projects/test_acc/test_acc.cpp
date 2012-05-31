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

using namespace std;

typedef test_accParameters taP;
Real taP::SPEED = NAN;
Real taP::v_min = NAN;
Real taP::v_max = NAN;

bool initializeProject(void) {return true;}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("test_acc.SPEED", "Acceleration value", 0.5);
   RP::add("test_acc.v_min", "Inner corner of initial blobs", 0.2);
   RP::add("test_acc.v_max", "Outer corner of initial blobs", 0.3);
   return true;
}
bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("test_acc.SPEED", taP::SPEED);
   RP::get("test_acc.v_min", taP::v_min);
   RP::get("test_acc.v_max", taP::v_max);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
   Real Vx = vx + 0.5 * dvx;
   Real Vy = vy + 0.5 * dvy;
   Real Vz = vz + 0.5 * dvz;
   if (Vz >= taP::v_min && Vz < taP::v_max) {
      if (Vx >= taP::v_min && Vx < taP::v_max) {
	 if (Vy >= -taP::v_max && Vy < -taP::v_min) return 1.0;
	 if (Vy >= taP::v_min && Vy < taP::v_max) return 1.0;
      }
      if (Vx >= -taP::v_max && Vx < -taP::v_min) {
	 if (Vy >= -taP::v_max && Vy < -taP::v_min) return 1.0;
	 if (Vy >= taP::v_min && Vy < taP::v_max) return 1.0;
      }
   }
   if (Vz >= -taP::v_max && Vz < -taP::v_min) {
      if (Vx >= taP::v_min && Vx < taP::v_max) {
	 if (Vy >= -taP::v_max && Vy < -taP::v_min) return 1.0;
	 if (Vy >= taP::v_min && Vy < taP::v_max) return 1.0;
      }
      if (Vx >= -taP::v_max && Vx < -taP::v_min) {
	 if (Vy >= -taP::v_max && Vy < -taP::v_min) return 1.0;
	 if (Vy >= taP::v_min && Vy < taP::v_max) return 1.0;
      }
   }
   return 0.0;
}

void calcBlockParameters(Real* blockParams) { }

void calcCellParameters(Real* cellParams,creal& t) {
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = 0.0;
   cellParams[CellParams::BY   ] = 0.0;
   cellParams[CellParams::BZ   ] = 0.0;
   
   cellParams[CellParams::EXVOL] = 0.0;
   cellParams[CellParams::EYVOL] = 0.0;
   cellParams[CellParams::EZVOL] = 0.0;
   cellParams[CellParams::BXVOL] = 0.0;
   cellParams[CellParams::BYVOL] = 0.0;
   cellParams[CellParams::BZVOL] = 0.0;
}



void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}



