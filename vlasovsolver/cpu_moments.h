/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef CPU_MOMENTS_H
#define CPU_MOMENTS_H

#include <vector>
#include <limits>
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include "../definitions.h"
#include "../common.h"
#include "../spatial_cell_wrapper.hpp"

using namespace spatial_cell;

// ***** FUNCTION DECLARATIONS ***** //

template<typename REAL> 
void blockVelocityFirstMoments(const Realf* avgs,const Real* blockParams,
                               REAL* array);

template<typename REAL> 
void blockVelocitySecondMoments(const Realf* avgs,const Real* blockParams,
                                const REAL v[3],
                                std::vector<Real>& array);

void calculateMoments_R(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
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
 * times population mass / proton mass. This function is VAMR safe.
 * @param avgs Distribution function.
 * @param blockParams Parameters for the given velocity block.
 * @param array Array of at least size four where the calculated moments are added.*/
template<typename REAL> inline
void blockVelocityFirstMoments(
        const Realf* avgs,
        const Real* blockParams,
        REAL* array) {

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
   
   const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
   array[0] += n_sum   * DV3;
   array[1] += nvx_sum * DV3;
   array[2] += nvy_sum * DV3;
   array[3] += nvz_sum * DV3;
}

/** Calculate the second velocity moments for the given velocity block, and add 
 * results to 'array', which must have at least size three. After this function 
 * returns, the contents of 'array' are as follows: array[0]=n(Vx-Vx0); 
 * array[1]=n(Vy-Vy0); array[2]=n(Vz-Vz0); Here Vx0,Vy0,Vz0 are the components 
 * of the bulk velocity (calculated over all species). This function is VAMR safe.
 * @param avgs Distribution function.
 * @param blockParams Parameters for the given velocity block.
 * @param averageVX Bulk velocity x
 * @param averageVY Bulk velocity y
 * @param averageVZ Bulk velocity z
 * @param array Array where the calculated moments are added.*/
template<typename REAL> inline
void blockVelocitySecondMoments(
        const Realf* avgs,
        const Real* blockParams,
        const REAL averageVX,
        const REAL averageVY,
        const REAL averageVZ,
        std::vector<Real>& array) {

   const Real HALF = 0.5;

   Real nvx2_sum = 0.0;
   Real nvy2_sum = 0.0;
   Real nvz2_sum = 0.0;
   Real nvyvz_sum = 0.0;
   Real nvxvz_sum = 0.0;
   Real nvxvy_sum = 0.0;
   for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
      const Real VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
      const Real VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
      const Real VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
      
      nvx2_sum += avgs[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX);
      nvy2_sum += avgs[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY);
      nvz2_sum += avgs[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ);

      nvyvz_sum += avgs[cellIndex(i, j, k)] * (VY - averageVY) * (VZ - averageVZ);
      nvxvz_sum += avgs[cellIndex(i, j, k)] * (VX - averageVX) * (VZ - averageVZ);
      nvxvy_sum += avgs[cellIndex(i, j, k)] * (VX - averageVX) * (VY - averageVY);
   }
   
   const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
   array[0] += nvx2_sum * DV3;
   array[1] += nvy2_sum * DV3;
   array[2] += nvz2_sum * DV3;
   array[3] += nvyvz_sum * DV3;
   array[4] += nvxvz_sum * DV3;
   array[5] += nvxvy_sum * DV3;
}

#endif
