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
   creal volume = cell.cpu_cellParams[CellParams::DX] * cell.cpu_cellParams[CellParams::DY] * cell.cpu_cellParams[CellParams::DZ];
   cell.cpu_cellParams[CellParams::RHO  ] /= volume;
   cell.cpu_cellParams[CellParams::RHOVX] /= volume;
   cell.cpu_cellParams[CellParams::RHOVY] /= volume;
   cell.cpu_cellParams[CellParams::RHOVZ] /= volume;

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

   creal volume = cell.cpu_cellParams[CellParams::DX] * cell.cpu_cellParams[CellParams::DY] * cell.cpu_cellParams[CellParams::DZ];
   cell.cpu_cellParams[CellParams::RHO  ] /= volume;
   cell.cpu_cellParams[CellParams::RHOVX] /= volume;
   cell.cpu_cellParams[CellParams::RHOVY] /= volume;
   cell.cpu_cellParams[CellParams::RHOVZ] /= volume;
   return true;
}
