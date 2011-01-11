#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "cpu_acc.h"

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

