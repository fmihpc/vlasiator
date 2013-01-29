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

#ifndef CPU_LORENTZ_H
#define CPU_LORENTZ_H

#include "definitions.h"
#include "common.h"
#include "spatial_cell.hpp"
#include "projects/project.h"
#include "leveque_common.h"
#include <algorithm>

using namespace spatial_cell;



//  Calculate bulk components of acceleration due to three-dimensional Lorentz force
//  cellParams Array containing spatial cell parameters.
//  cellBVOLDerivatives Array containing the BVOL derivatives.

template<typename REAL> void lorentzForceFaceBulk(
   REAL a_bulk[3],
   const REAL* const cellParams,
   const REAL* const cellBVOLDerivatives
) {

   
   const REAL UX = ((cellParams[CellParams::RHO_V] <= 0) ? 0.0 : cellParams[CellParams::RHOVX_V]/cellParams[CellParams::RHO_V]);
   const REAL UY = ((cellParams[CellParams::RHO_V] <= 0) ? 0.0 : cellParams[CellParams::RHOVY_V]/cellParams[CellParams::RHO_V]);
   const REAL UZ = ((cellParams[CellParams::RHO_V] <= 0) ? 0.0 : cellParams[CellParams::RHOVZ_V]/cellParams[CellParams::RHO_V]);
   
   const REAL BX = cellParams[CellParams::PERBXVOL] + cellParams[CellParams::BGBXVOL];
   const REAL BY = cellParams[CellParams::PERBYVOL] + cellParams[CellParams::BGBYVOL];
   const REAL BZ = cellParams[CellParams::PERBZVOL] + cellParams[CellParams::BGBZVOL];
   const REAL B = sqrt(BX*BX+BY*BY+BZ*BZ);

   const REAL EX = cellParams[CellParams::EXVOL];
   const REAL EY = cellParams[CellParams::EYVOL];
   const REAL EZ = cellParams[CellParams::EZVOL];
   
   const REAL dBXdy = cellBVOLDerivatives[bvolderivatives::dPERBXVOLdy]/cellParams[CellParams::DY];
   const REAL dBXdz = cellBVOLDerivatives[bvolderivatives::dPERBXVOLdz]/cellParams[CellParams::DZ];
   const REAL dBYdx = cellBVOLDerivatives[bvolderivatives::dPERBYVOLdx]/cellParams[CellParams::DX];
   const REAL dBYdz = cellBVOLDerivatives[bvolderivatives::dPERBYVOLdz]/cellParams[CellParams::DZ];
   const REAL dBZdx = cellBVOLDerivatives[bvolderivatives::dPERBZVOLdx]/cellParams[CellParams::DX];
   const REAL dBZdy = cellBVOLDerivatives[bvolderivatives::dPERBZVOLdy]/cellParams[CellParams::DY];
   
   const REAL prefactor =  ((cellParams[CellParams::RHO_V] <= 0) ? 0.0 : 1.0 / (physicalconstants::MU_0 * cellParams[CellParams::RHO_V] * Parameters::q )); 

   
   if(Parameters::lorentzUseFieldSolverE) {
      a_bulk[0] = Parameters::q_per_m*EX;
      a_bulk[1] = Parameters::q_per_m*EY;
      a_bulk[2] = Parameters::q_per_m*EZ;
   }
   else {
      a_bulk[0] = Parameters::q_per_m*( BY*UZ - BZ*UY);
      a_bulk[1] = Parameters::q_per_m*( BZ*UX - BX*UZ);
      a_bulk[2] = Parameters::q_per_m*( BX*UY - BY*UX);
   }


   //resisitivity term E=n J
   a_bulk[0]+=  Parameters::q_per_m*Parameters::resistivity * B *prefactor *(dBZdy-dBYdz);
   a_bulk[1]+=  Parameters::q_per_m*Parameters::resistivity * B *prefactor *(dBXdz-dBZdx);
   a_bulk[2]+=  Parameters::q_per_m*Parameters::resistivity * B *prefactor *(dBYdx-dBXdy);
   
   /*
     Hall term  (see QUESPACE-281)  (JxB = curl(B) x B)
   */
   
   if(Parameters::lorentzHallTerm) {
      a_bulk[0] += Parameters::q_per_m * prefactor * (BZ*dBXdz - BZ*dBZdx - BY*dBYdx + BY*dBXdy);
      a_bulk[1] += Parameters::q_per_m * prefactor * (BX*dBYdx - BX*dBXdy - BZ*dBZdy + BZ*dBYdz);
      a_bulk[2] += Parameters::q_per_m * prefactor * (BY*dBZdy - BY*dBYdz - BX*dBXdz + BX*dBZdx);
   }
}




//  Calculate conponents of acceleration due to three-dimensional Lorentz force caused by phase space in velocity space
//  cellParams Array containing spatial cell parameters.
//  cellBVOLDerivatives Array containing the BVOL derivatives.


template<typename UINT,typename REAL> void lorentzForceFaceVelspace(
   REAL a_phasespacecell[9],
   const UINT& I, const UINT& J, const UINT& K, 
   const REAL* const cellParams, 
   const REAL* const blockParams) {
   
   const REAL BX = cellParams[CellParams::PERBXVOL] + cellParams[CellParams::BGBXVOL];
   const REAL BY = cellParams[CellParams::PERBYVOL] + cellParams[CellParams::BGBYVOL];
   const REAL BZ = cellParams[CellParams::PERBZVOL] + cellParams[CellParams::BGBZVOL];
   
   
   // I The i-index of the cell, within a velocity block, in which the computed face is stored.
   // J The j-index of the cell, within a velocity block, in which the computed face is stored.
   // K The k-index of the cell, within a velocity block, in which the computed face is stored.
   
   // Calculate block components of acceleration due to three-dimensional Lorentz force at face with normal to vx-direction.  
   REAL VX = blockParams[BlockParams::VXCRD] + I*blockParams[BlockParams::DVX];
   REAL VY = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   REAL VZ = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   
   // ax=a_phasespacecell[0], ay=a_phasespacecell[1], az=a_phasespacecell[2]
   a_phasespacecell[0] = Parameters::q_per_m*(VY*BZ - VZ*BY);
   a_phasespacecell[1] = Parameters::q_per_m*(VZ*BX - VX*BZ);
   a_phasespacecell[2] = Parameters::q_per_m*(VX*BY - VY*BX);
   
   // Calculate block components of acceleration due to three-dimensional Lorentz force at face with normal to vy-direction.
   
   VX = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   VY = blockParams[BlockParams::VYCRD] + J*blockParams[BlockParams::DVY];
   VZ = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
      
   // ax=a_phasespacecell[3], ay=a_phasespacecell[4], az=a_phasespacecell[5]
   a_phasespacecell[3] = Parameters::q_per_m*(VY*BZ - VZ*BY);
   a_phasespacecell[4] = Parameters::q_per_m*(VZ*BX - VX*BZ);
   a_phasespacecell[5] = Parameters::q_per_m*(VX*BY - VY*BX);
   
   //Calculate block components of acceleration due to three-dimensional Lorentz force at face with normal to vz-direction.
   VX = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   VY = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   VZ = blockParams[BlockParams::VZCRD] + K*blockParams[BlockParams::DVZ];
   
   // ax=a_phasespacecell[6], ay=a_phasespacecell[7], az=a_phasespacecell[8]
   a_phasespacecell[6] = Parameters::q_per_m*(VY*BZ - VZ*BY);
   a_phasespacecell[7] = Parameters::q_per_m*(VZ*BX - VX*BZ);
   a_phasespacecell[8] = Parameters::q_per_m*(VX*BY - VY*BX);
}





#endif
