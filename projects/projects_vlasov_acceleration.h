/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PROJECTS_ACCELERATION_COMMON_H
#define PROJECTS_ACCELERATION_COMMON_H

// *********************************
// ***** TEMPLATE DECLARATIONS *****
// *********************************

template<typename UINT,typename REAL> void lorentzForceFaceX(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
);

template<typename UINT,typename REAL> void lorentzForceFaceY(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
);

template<typename UINT,typename REAL> void lorentzForceFaceZ(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
);

// ********************************
// ***** TEMPLATE DEFINITIONS *****
// ********************************

/** Calculate acceleration due to three-dimensional Lorentz force
 * at face with normal to vx-direction.
 * @param ax Variable in which the calculated x-component of acceleration is to be written.
 * @param ay Variable in which the calculated y-component of acceleration is to be written.
 * @param az Variable in which the calculated z-component of acceleration is to be written.
 * @param I The i-index of the cell, within a velocity block, in which the computed face is stored.
 * @param J The j-index of the cell, within a velocity block, in which the computed face is stored.
 * @param K The k-index of the cell, within a velocity block, in which the computed face is stored.
 * @param cellParams Array containing spatial cell parameters.
 * @param blockParams Array containing velocity block parameters.
 * @param cellBVOLDerivatives Array containing the BVOL derivatives.
 */
template<typename UINT,typename REAL> void lorentzForceFaceX(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
) {
   const REAL VX = blockParams[BlockParams::VXCRD] + I*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   const REAL UX = cellParams[CellParams::RHO] == 0 ? 0.0 : cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   const REAL UY = cellParams[CellParams::RHO] == 0 ? 0.0 : cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   const REAL UZ = cellParams[CellParams::RHO] == 0 ? 0.0 : cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   const REAL BX = cellParams[CellParams::BXVOL];
   const REAL BY = cellParams[CellParams::BYVOL];
   const REAL BZ = cellParams[CellParams::BZVOL];
   const REAL dBXdy = cellBVOLDerivatives[bvolderivatives::dBXVOLdy];
   const REAL dBXdz = cellBVOLDerivatives[bvolderivatives::dBXVOLdz];
   const REAL dBYdx = cellBVOLDerivatives[bvolderivatives::dBYVOLdx];
   const REAL dBYdz = cellBVOLDerivatives[bvolderivatives::dBYVOLdz];
   const REAL dBZdx = cellBVOLDerivatives[bvolderivatives::dBZVOLdx];
   const REAL dBZdy = cellBVOLDerivatives[bvolderivatives::dBZVOLdy];
   const REAL prefactor = cellParams[CellParams::RHO] == 0 ? 0.0 : 1.0 / (physicalconstants::MU_0 * cellParams[CellParams::RHO]);
   ax = Parameters::q_per_m*(prefactor * (BZ*dBXdz - BZ*dBZdx - BY*dBYdx + BY*dBXdy) + BY*(VZ-UZ) - BZ*(VY-UY));
   ay = Parameters::q_per_m*(prefactor * (BX*dBYdx - BX*dBXdy - BZ*dBZdy + BZ*dBYdz) + BZ*(VX-UX) - BX*(VZ-UZ));
   az = Parameters::q_per_m*(prefactor * (BY*dBZdy - BY*dBYdz - BX*dBXdz + BX*dBZdx) + BX*(VY-UY) - BY*(VX-UX));
}

/** Calculate acceleration due to three-dimensional Lorentz force
 * at face with normal to vy-direction.
 * @param ax Variable in which the calculated x-component of acceleration is to be written.
 * @param ay Variable in which the calculated y-component of acceleration is to be written.
 * @param az Variable in which the calculated z-component of acceleration is to be written.
 * @param I The i-index of the cell, within a velocity block, in which the computed face is stored.
 * @param J The j-index of the cell, within a velocity block, in which the computed face is stored.
 * @param K The k-index of the cell, within a velocity block, in which the computed face is stored.
 * @param cellParams Array containing spatial cell parameters.
 * @param blockParams Array containing velocity block parameters.
 * @param cellBVOLDerivatives Array containing the BVOL derivatives.
 */
template<typename UINT,typename REAL> void lorentzForceFaceY(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
) {
   const REAL VX = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + J*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   const REAL UX = cellParams[CellParams::RHO] == 0 ? 0.0 : cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   const REAL UY = cellParams[CellParams::RHO] == 0 ? 0.0 : cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   const REAL UZ = cellParams[CellParams::RHO] == 0 ? 0.0 : cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   const REAL BX = cellParams[CellParams::BXVOL];
   const REAL BY = cellParams[CellParams::BYVOL];
   const REAL BZ = cellParams[CellParams::BZVOL];
   const REAL dBXdy = cellBVOLDerivatives[bvolderivatives::dBXVOLdy];
   const REAL dBXdz = cellBVOLDerivatives[bvolderivatives::dBXVOLdz];
   const REAL dBYdx = cellBVOLDerivatives[bvolderivatives::dBYVOLdx];
   const REAL dBYdz = cellBVOLDerivatives[bvolderivatives::dBYVOLdz];
   const REAL dBZdx = cellBVOLDerivatives[bvolderivatives::dBZVOLdx];
   const REAL dBZdy = cellBVOLDerivatives[bvolderivatives::dBZVOLdy];
   const REAL prefactor = cellParams[CellParams::RHO] == 0 ? 0.0 : 1.0 / (physicalconstants::MU_0 * cellParams[CellParams::RHO]);
   ax = Parameters::q_per_m*(prefactor * (BZ*dBXdz - BZ*dBZdx - BY*dBYdx + BY*dBXdy) + BY*(VZ-UZ) - BZ*(VY-UY));
   ay = Parameters::q_per_m*(prefactor * (BX*dBYdx - BX*dBXdy - BZ*dBZdy + BZ*dBYdz) + BZ*(VX-UX) - BX*(VZ-UZ));
   az = Parameters::q_per_m*(prefactor * (BY*dBZdy - BY*dBYdz - BX*dBXdz + BX*dBZdx) + BX*(VY-UY) - BY*(VX-UX));
}

/** Calculate acceleration due to three-dimensional Lorentz force
 * at face with normal to vz-direction.
 * @param ax Variable in which the calculated x-component of acceleration is to be written.
 * @param ay Variable in which the calculated y-component of acceleration is to be written.
 * @param az Variable in which the calculated z-component of acceleration is to be written.
 * @param I The i-index of the cell, within a velocity block, in which the computed face is stored.
 * @param J The j-index of the cell, within a velocity block, in which the computed face is stored.
 * @param K The k-index of the cell, within a velocity block, in which the computed face is stored.
 * @param cellParams Array containing spatial cell parameters.
 * @param blockParams Array containing velocity block parameters.
 * @param cellBVOLDerivatives Array containing the BVOL derivatives.
 */
template<typename UINT,typename REAL> void lorentzForceFaceZ(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
) {
   const REAL VX = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + K*blockParams[BlockParams::DVZ];
   const REAL UX = cellParams[CellParams::RHO] == 0 ? 0.0 : cellParams[CellParams::RHOVX]/cellParams[CellParams::RHO];
   const REAL UY = cellParams[CellParams::RHO] == 0 ? 0.0 : cellParams[CellParams::RHOVY]/cellParams[CellParams::RHO];
   const REAL UZ = cellParams[CellParams::RHO] == 0 ? 0.0 : cellParams[CellParams::RHOVZ]/cellParams[CellParams::RHO];
   const REAL BX = cellParams[CellParams::BXVOL];
   const REAL BY = cellParams[CellParams::BYVOL];
   const REAL BZ = cellParams[CellParams::BZVOL];
   const REAL dBXdy = cellBVOLDerivatives[bvolderivatives::dBXVOLdy];
   const REAL dBXdz = cellBVOLDerivatives[bvolderivatives::dBXVOLdz];
   const REAL dBYdx = cellBVOLDerivatives[bvolderivatives::dBYVOLdx];
   const REAL dBYdz = cellBVOLDerivatives[bvolderivatives::dBYVOLdz];
   const REAL dBZdx = cellBVOLDerivatives[bvolderivatives::dBZVOLdx];
   const REAL dBZdy = cellBVOLDerivatives[bvolderivatives::dBZVOLdy];
   const REAL prefactor = cellParams[CellParams::RHO] == 0 ? 0.0 : 1.0 / (physicalconstants::MU_0 * cellParams[CellParams::RHO]);
   ax = Parameters::q_per_m*(prefactor * (BZ*dBXdz - BZ*dBZdx - BY*dBYdx + BY*dBXdy) + BY*(VZ-UZ) - BZ*(VY-UY));
   ay = Parameters::q_per_m*(prefactor * (BX*dBYdx - BX*dBXdy - BZ*dBZdy + BZ*dBYdz) + BZ*(VX-UX) - BX*(VZ-UZ));
   az = Parameters::q_per_m*(prefactor * (BY*dBZdy - BY*dBYdz - BX*dBXdz + BX*dBZdx) + BX*(VY-UY) - BY*(VX-UX));
}
#endif
