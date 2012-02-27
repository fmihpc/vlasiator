/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

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
#include <cmath>
#include <omp.h>

#include "cpu_acc_kt.h"

using namespace std;

bool cpu_acceleration(SpatialCell& cell) {
   // Clear spatial cell velocity moments:
   cell.cpu_cellParams[CellParams::RHO]   = 0.0;
   cell.cpu_cellParams[CellParams::RHOVX] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVY] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVZ] = 0.0;
   
   bool success = true;
   #pragma omp parallel 
     {
	creal DT = Parameters::dt;
	// Calculate derivatives in vx,vy,vz:
        #pragma omp for
	for (uint block=0; block<cell.N_blocks; ++block) {
	   cpu_calcVelDerivs<Real>(cell,block);
	}
	
	// Calculate velocity fluxes:
        #pragma omp for
	for (uint block=0; block<cell.N_blocks; ++block) {
	   cpu_calcVelFluxesX<Real>(cell,block);
	   cpu_calcVelFluxesY<Real>(cell,block);
	   cpu_calcVelFluxesZ<Real>(cell,block);
	}
	
	// Propagate volume averages in velocity space:
        #pragma omp for
	for (uint block=0; block<cell.N_blocks; ++block) {
	   cpu_propagateVel<Real>(cell,block,DT);
	}
     }
   return success;
}

