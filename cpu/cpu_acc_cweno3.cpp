/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

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


