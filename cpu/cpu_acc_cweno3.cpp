#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "cpu_acc_cweno3.h"

using namespace std;

bool cpu_acceleration(SpatialCell& cell) {
   // Clear spatial cell velocity moments:
   cell.cpu_cellParams[CellParams::RHO]   = 0.0;
   cell.cpu_cellParams[CellParams::RHOVX] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVY] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVZ] = 0.0;
   
   bool success = true;
   Real DT = Parameters::dt;

   copyAverages(cell);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelDerivsX<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelFluxesX<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_propagateVelX<Real>(cell,block,cell.cpu_avgs,DT);
   //for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelDerivsX<Real>(cell,block);
   //for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelFluxesX<Real>(cell,block);
   //for (uint block=0; block<cell.N_blocks; ++block) cpu_propagateVelX<Real>(cell,block,cell.cpu_avgs,DT);
   sumAverages(cell);
   
   /*
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelDerivsY<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelFluxesY<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_propagateVelY<Real>(cell,block,cell.cpu_avgs,1.0*DT);
   */
   /*
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelDerivsX<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelFluxesX<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_propagateVelX<Real>(cell,block,cell.cpu_avgs,0.5*DT);
    */
}


