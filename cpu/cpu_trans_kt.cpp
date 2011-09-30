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

#include "cpu_trans_kt.h"

using namespace std;

bool cpu_translation1(SpatialCell& cell,const std::vector<const SpatialCell*>& spatNbrs) {
   #pragma omp parallel for
   for (uint block=0; block<cell.N_blocks; ++block) {
      cpu_calcSpatDerivs<Real>(cell,block,spatNbrs);
   }
   return true;
}

bool cpu_translation2(SpatialCell& cell,const std::vector<const SpatialCell*>& spatNbrs) {
   #pragma omp parallel for
   for (uint block=0; block<cell.N_blocks; ++block) {
      cpu_calcSpatFluxesX<Real>(cell,block,spatNbrs);
      cpu_calcSpatFluxesY<Real>(cell,block,spatNbrs);
      cpu_calcSpatFluxesZ<Real>(cell,block,spatNbrs);
   }
   return true;
}

bool cpu_translation3(SpatialCell& cell,const std::vector<const SpatialCell*>& spatNbrs) {
   creal DT = Parameters::dt;
   #pragma omp parallel for
   for (uint block=0; block<cell.N_blocks; ++block) {
      cpu_propagateSpat(cell,block,spatNbrs,DT);
   }
   return true;
}

/** Calculate velocity moments for given spatial cell. This function is only called 
 * when plotting the initial state of the simulation, solver calculates the moments 
 * afterwards.
 * @param cell The spatial cell whose velocity moments are to be calculated.
 * @return If true, the velocity moments were calculated successfully.
 */
bool cpu_calcVelocityMoments(SpatialCell& cell) {
   cell.cpu_cellParams[CellParams::RHO  ] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVX] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVY] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVZ] = 0.0;

   #pragma omp parallel for
   for (uint block=0; block<cell.N_blocks; ++block) {
      creal* const blockParams = cell.cpu_blockParams + block*SIZE_BLOCKPARAMS;
      const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
      Real sum_n = 0.0;
      Real sum_nvx = 0.0;
      Real sum_nvy = 0.0;
      Real sum_nvz = 0.0;
      for (uint k=0; k<WID; ++k) {
	 const Real VZ = blockParams[BlockParams::VZCRD] + (k+convert<Real>(0.5))*blockParams[BlockParams::DVZ];
	 for (uint j=0; j<WID; ++j) {
	    const Real VY = blockParams[BlockParams::VYCRD] + (j+convert<Real>(0.5))*blockParams[BlockParams::DVY];
	    for (uint i=0; i<WID; ++i) {
	       const Real VX = blockParams[BlockParams::VXCRD] + (i+convert<Real>(0.5))*blockParams[BlockParams::DVX];
	       sum_n   += cell.cpu_avgs[block*SIZE_VELBLOCK + k*WID2+j*WID+i];
	       sum_nvx += cell.cpu_avgs[block*SIZE_VELBLOCK + k*WID2+j*WID+i] * VX;
	       sum_nvy += cell.cpu_avgs[block*SIZE_VELBLOCK + k*WID2+j*WID+i] * VY;
	       sum_nvz += cell.cpu_avgs[block*SIZE_VELBLOCK + k*WID2+j*WID+i] * VZ;
	    }
	 }
      }
      #pragma omp atomic
      cell.cpu_cellParams[CellParams::RHO  ] += sum_n  *DV3;
      #pragma omp atomic
      cell.cpu_cellParams[CellParams::RHOVX] += sum_nvx*DV3;
      #pragma omp atomic
      cell.cpu_cellParams[CellParams::RHOVY] += sum_nvy*DV3;
      #pragma omp atomic
      cell.cpu_cellParams[CellParams::RHOVZ] += sum_nvz*DV3;
   }
   return true;
}


