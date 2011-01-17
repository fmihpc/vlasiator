#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "cpu_trans.h"

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

