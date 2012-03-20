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

using namespace std;

bool initializeProject(void) {
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
   if (vz >= 0.2 && vz < 0.2+dvz) {
      if (vx >= 0.2 && vx < 0.2+dvx) {
	 if (vy >= -0.2-dvy*1.001 && vy < -0.2) return 1.0;
	 if (vy >= 0.2 && vy < 0.2+dvy) return 1.0;
      }
      if (vx >= -0.2-dvx*1.001 && vx < -0.2) {
	 if (vy >= -0.2-dvy*1.001 && vy < -0.2) return 1.0;
	 if (vy >= 0.2 && vy < 0.2+dvy) return 1.0;
      }
   }
   if (vz >= -0.2-dvz*1.001 && vz < -0.2) {
      if (vx >= 0.2 && vx < 0.2+dvx) {
	 if (vy >= -0.2-dvy*1.001 && vy < -0.2) return 1.0;
	 if (vy >= 0.2 && vy < 0.2+dvy) return 1.0;
      }
      if (vx >= -0.2-dvx*1.001 && vx < -0.2) {
	 if (vy >= -0.2-dvy*1.001 && vy < -0.2) return 1.0;
	 if (vy >= 0.2 && vy < 0.2+dvy) return 1.0;
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



