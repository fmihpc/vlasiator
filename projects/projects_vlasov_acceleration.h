#ifndef PROJECTS_ACCELERATION_COMMON_H
#define PROJECTS_ACCELERATION_COMMON_H

// *********************************
// ***** TEMPLATE DECLARATIONS *****
// *********************************

template<typename UINT,typename REAL> void lorentzForceFaceX(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,
							     const REAL* const cellParams,const REAL* const blockParams);

template<typename UINT,typename REAL> void lorentzForceFaceY(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,
							     const REAL* const cellParams,const REAL* const blockParams);

template<typename UINT,typename REAL> void lorentzForceFaceZ(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,
							     const REAL* const cellParams,const REAL* const blockParams);

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
 */
template<typename UINT,typename REAL> void lorentzForceFaceX(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,
							const REAL* const cellParams,const REAL* const blockParams) {
   const REAL VX = blockParams[BlockParams::VXCRD] + I*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   const REAL EX = cellParams[CellParams::EX];
   const REAL EY = cellParams[CellParams::EY];
   const REAL EZ = cellParams[CellParams::EZ];
   const REAL BX = cellParams[CellParams::BX];
   const REAL BY = cellParams[CellParams::BY];
   const REAL BZ = cellParams[CellParams::BZ];
   ax = Parameters::q_per_m*(EX + VY*BZ - VZ*BY);
   ay = Parameters::q_per_m*(EY + VZ*BX - VX*BZ);
   az = Parameters::q_per_m*(EZ + VX*BY - VY*BX);
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
 */
template<typename UINT,typename REAL> void lorentzForceFaceY(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,
							const REAL* const cellParams,const REAL* const blockParams) {
   const REAL VX = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + J*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   const REAL EX = cellParams[CellParams::EX];
   const REAL EY = cellParams[CellParams::EY];
   const REAL EZ = cellParams[CellParams::EZ];
   const REAL BX = cellParams[CellParams::BX];
   const REAL BY = cellParams[CellParams::BY];
   const REAL BZ = cellParams[CellParams::BZ];
   ax = Parameters::q_per_m*(EX + VY*BZ - VZ*BY);
   ay = Parameters::q_per_m*(EY + VZ*BX - VX*BZ);
   az = Parameters::q_per_m*(EZ + VX*BY - VY*BX);
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
 */
template<typename UINT,typename REAL> void lorentzForceFaceZ(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,
							const REAL* const cellParams,const REAL* const blockParams) {
   const REAL VX = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + K*blockParams[BlockParams::DVZ];
   const REAL EX = cellParams[CellParams::EX];
   const REAL EY = cellParams[CellParams::EY];
   const REAL EZ = cellParams[CellParams::EZ];
   const REAL BX = cellParams[CellParams::BX];
   const REAL BY = cellParams[CellParams::BY];
   const REAL BZ = cellParams[CellParams::BZ];
   ax = Parameters::q_per_m*(EX + VY*BZ - VZ*BY);
   ay = Parameters::q_per_m*(EY + VZ*BX - VX*BZ);
   az = Parameters::q_per_m*(EZ + VX*BY - VY*BX);
}

#endif
