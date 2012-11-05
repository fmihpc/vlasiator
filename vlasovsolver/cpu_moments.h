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

#ifndef CPU_MOMENTS_H
#define CPU_MOMENTS_H

#include <vector>
#include <limits>
#include <dccrg.hpp>

#include "../definitions.h"
#include "../common.h"
#include "leveque_common.h"
#include "spatial_cell.hpp"
using namespace spatial_cell;

template<typename REAL> void cpu_blockVelocityMoments(const Real* const avgs,const REAL* const blockParams,REAL* const cellParams,
                                                      const int cp_rho, const int cp_rhovx, const int cp_rhovy, const int cp_rhovz) {
   const REAL HALF = 0.5;
   
   REAL n_sum = 0.0;
   REAL nvx_sum = 0.0;
   REAL nvy_sum = 0.0;
   REAL nvz_sum = 0.0;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      const REAL VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
      const REAL VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
      const REAL VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
      
      n_sum   += avgs[cellIndex(i,j,k)];
      nvx_sum += avgs[cellIndex(i,j,k)]*VX;
      nvy_sum += avgs[cellIndex(i,j,k)]*VY;
      nvz_sum += avgs[cellIndex(i,j,k)]*VZ;
   }
   
   // Accumulate contributions coming from this velocity block to the 
   // spatial cell velocity moments. If multithreading / OpenMP is used, 
   // these updates need to be atomic:
   const REAL DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
   cellParams[cp_rho] += n_sum * DV3;
   cellParams[cp_rhovx] += nvx_sum * DV3;
   cellParams[cp_rhovy] += nvy_sum * DV3;
   cellParams[cp_rhovz] += nvz_sum * DV3;
}





template<typename UINT> void cpu_calcVelocityMoments(SpatialCell *cell,const UINT blockId,
                                                     const int cp_rho, const int cp_rhovx, const int cp_rhovy, const int cp_rhovz) {
   Velocity_Block* block=cell->at(blockId); //returns a reference to block            
   // Calculate velocity moments:
   cpu_blockVelocityMoments(block->data,block->parameters,cell->parameters,cp_rho,cp_rhovx,cp_rhovy,cp_rhovz);

}




#endif
      


 
