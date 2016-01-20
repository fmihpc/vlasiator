/*
This file is part of Vlasiator.

Copyright 2010-2013,2015 Finnish Meteorological Institute

*/

#ifndef CPU_MOMENTS_H
#define CPU_MOMENTS_H

#include <vector>
#include <limits>
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include "../definitions.h"
#include "../common.h"
#include "../spatial_cell.hpp"

using namespace spatial_cell;

// ***** FUNCTION DECLARATIONS ***** //

template<typename REAL> 
void blockVelocityFirstMoments(const Realf* avgs,const Real* blockParams,
                               const Real& massRatio,REAL* array);

template<typename REAL> 
void blockVelocitySecondMoments(const Realf* avgs,const Real* blockParams,const Real* cellParams,
                                const int cp_rho,const int cp_rhovx,const int cp_rhovy,const int cp_rhovz,
                                REAL* array);

void calculateMoments_R_maxdt(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              const std::vector<CellID>& cells,
                              const bool& computeSecond);

void calculateMoments_V(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                        const std::vector<CellID>& cells,
                        const bool& computeSecond);



// ***** TEMPLATE FUNCTION DEFINITIONS ***** //

/** Calculate the zeroth and first velocity moments for the given 
 * velocity block and add results to 'array', which must have at 
 * least size four. After this function returns, the contents of 
 * 'array' are as follows: array[0]=n; array[1]=n*Vx; array[2]=nVy;
 * array[3]=nVz; Here n is the scaled number density, i.e., number density 
 * times population mass / proton mass. This function is AMR safe.
 * @param avgs Distribution function.
 * @param blockParams Parameters for the given velocity block.
 * @param massRatio Population mass / proton mass.
 * @param array Array of at least size four where the calculated moments are added.*/
template<typename REAL> inline
void blockVelocityFirstMoments(
        const Realf* avgs,
        const Real* blockParams,
        const Real& massRatio,REAL* array) {

    const Real HALF = 0.5;

   Real n_sum = 0.0;
   Real nvx_sum = 0.0;
   Real nvy_sum = 0.0;
   Real nvz_sum = 0.0;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      const REAL VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
      const REAL VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
      const REAL VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
      
      n_sum   += avgs[cellIndex(i,j,k)];
      nvx_sum += avgs[cellIndex(i,j,k)]*VX;
      nvy_sum += avgs[cellIndex(i,j,k)]*VY;
      nvz_sum += avgs[cellIndex(i,j,k)]*VZ;        
   }
   
   const Real mrDV3 = massRatio * blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
   array[0] += n_sum   * mrDV3;
   array[1] += nvx_sum * mrDV3;
   array[2] += nvy_sum * mrDV3;
   array[3] += nvz_sum * mrDV3;
}

/** Calculate the second velocity moments for the given velocity block, and add 
 * results to 'array', which must have at least size three. After this function 
 * returns, the contents of 'array' are as follows: array[0]=n(Vx-Vx0); 
 * array[1]=n(Vy-Vy0); array[2]=n(Vz-Vz0); Here Vx0,Vy0,Vz0 are the components 
 * of the bulk velocity (calculated over all species). This function is AMR safe.
 * @param avgs Distribution function.
 * @param blockParams Parameters for the given velocity block.
 * @param cellParams Parameters for the spatial cell containing the given velocity block.
 * @param rho Index into cellParams, used to read bulk number density.
 * @param rhovx Index into cellParams, used to read bulk Vx times number density.
 * @param rhovy Index into cellParams, used to read bulk Vy times number density.
 * @param rhovz Index into cellParams, used to read bulk Vz times number density.
 * @param array Array where the calculated moments are added.*/
template<typename REAL> inline
void blockVelocitySecondMoments(
        const Realf* avgs,
        const Real* blockParams,
        const Real* cellParams,
        const int cp_rho,
        const int cp_rhovx,
        const int cp_rhovy,
        const int cp_rhovz,
        REAL* array) {

   const Real HALF = 0.5;

   const Real RHO = std::max(cellParams[cp_rho], std::numeric_limits<REAL>::min());
   const Real averageVX = cellParams[cp_rhovx] / RHO;
   const Real averageVY = cellParams[cp_rhovy] / RHO;
   const Real averageVZ = cellParams[cp_rhovz] / RHO;
   Real nvx2_sum = 0.0;
   Real nvy2_sum = 0.0;
   Real nvz2_sum = 0.0;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      const Real VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
      const Real VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
      const Real VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
      
      nvx2_sum += avgs[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX);
      nvy2_sum += avgs[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY);
      nvz2_sum += avgs[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ);
   }
   
   const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
   array[0] += nvx2_sum * DV3;
   array[1] += nvy2_sum * DV3;
   array[2] += nvz2_sum * DV3;
}

#endif
