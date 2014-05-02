/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute

*/

#ifndef CPU_MOMENTS_H
#define CPU_MOMENTS_H

#include <vector>
#include <limits>
#include <dccrg.hpp>

#include "../definitions.h"
#include "../common.h"
#include "spatial_cell.hpp"
using namespace spatial_cell;

template<typename REAL> void cpu_blockVelocityFirstMoments(
   const Realf* const avgs,
   const REAL* const blockParams,
   REAL* const cellParams,
   const int cp_rho,
   const int cp_rhovx,
   const int cp_rhovy,
   const int cp_rhovz
) {
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

template<typename REAL> void cpu_blockVelocitySecondMoments(
   const Realf* const avgs,
   const REAL* const blockParams,
   REAL* const cellParams,
   const int cp_rho,
   const int cp_rhovx,
   const int cp_rhovy,
   const int cp_rhovz,
   const int cp_p11,
   const int cp_p22,
   const int cp_p33
) {
   const REAL HALF = 0.5;
   
   REAL averageVX = cellParams[cp_rhovx] / cellParams[cp_rho];
   REAL averageVY = cellParams[cp_rhovy] / cellParams[cp_rho];
   REAL averageVZ = cellParams[cp_rhovz] / cellParams[cp_rho];
   REAL nvx2_sum = 0.0;
   REAL nvy2_sum = 0.0;
   REAL nvz2_sum = 0.0;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      const REAL VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
      const REAL VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
      const REAL VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
      
      nvx2_sum += avgs[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX);
      nvy2_sum += avgs[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY);
      nvz2_sum += avgs[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ);
   }
   
   // Accumulate contributions coming from this velocity block to the 
   // spatial cell velocity moments. If multithreading / OpenMP is used, 
   // these updates need to be atomic:
   const REAL mpDV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ]*physicalconstants::MASS_PROTON;
   cellParams[cp_p11] += nvx2_sum * mpDV3;
   cellParams[cp_p22] += nvy2_sum * mpDV3;
   cellParams[cp_p33] += nvz2_sum * mpDV3;
}



template<typename UINT> void cpu_calcVelocityFirstMoments(
   SpatialCell *cell,
   const UINT blockId,
   const int cp_rho,
   const int cp_rhovx,
   const int cp_rhovy,
   const int cp_rhovz
) {
   Velocity_Block* block=cell->at(blockId); //returns a reference to block            
   // Calculate velocity moments:
   cpu_blockVelocityFirstMoments(block->data,block->parameters,cell->parameters,cp_rho,cp_rhovx,cp_rhovy,cp_rhovz);
}

template<typename UINT> void cpu_calcVelocitySecondMoments(
   SpatialCell *cell,
   const UINT blockId,
   const int cp_rho,
   const int cp_rhovx,
   const int cp_rhovy,
   const int cp_rhovz,
   const int cp_p11,
   const int cp_p22,
   const int cp_p33
) {
   Velocity_Block* block=cell->at(blockId); //returns a reference to block
   // Calculate velocity moments:
   cpu_blockVelocitySecondMoments(block->data,block->parameters,cell->parameters,cp_rho,cp_rhovx,cp_rhovy,cp_rhovz,cp_p11,cp_p22,cp_p33);
}


#endif
      


 
