#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "cpu_acc.h"

using namespace std;

bool cpu_acceleration(SpatialCell& cell) {
   bool success = true;
   #pragma omp parallel 
     {
	creal DT = Parameters::dt;
	// Calculate derivatives in vx,vy,vz:
        #pragma omp for nowait
	for (uint block=0; block<cell.N_blocks; ++block) {
	   cpu_calcVelDerivs<real>(cell,block);
	}
	
	// Calculate velocity fluxes:
        #pragma omp for nowait
	for (uint block=0; block<cell.N_blocks; ++block) {
	   cpu_calcVelFluxesX<real>(cell,block);
	   cpu_calcVelFluxesY<real>(cell,block);
	   cpu_calcVelFluxesZ<real>(cell,block);
	}

	// Propagate volume averages in velocity space:
        #pragma omp for nowait
	for (uint block=0; block<cell.N_blocks; ++block) {
	   cpu_propagateVel<real>(cell,block,DT);
	}
     }
   return success;
}

