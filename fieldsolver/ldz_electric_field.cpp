/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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

#include <cstdlib>

#include "fs_common.h"
#include "fs_cache.h"
#include "ldz_electric_field.hpp"

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

namespace fs = fieldsolver;
namespace pc = physicalconstants;
using namespace std;

extern map<CellID,uint> existingCellsFlags;

/*! \brief Low-level helper function.
 *
 * Computes the correct combination of speeds to determine the CFL limits.
 *
 * It should be in-plane, but we use the complete wave speeds from calculateWaveSpeed??().
 * 
 * At the moment it computes the geometric mean of both bulk velocity components
 * and takes the maximum of that plus either the magnetosonic or the whistler speed.
 * 
 * \sa calculateWaveSpeedYZ calculateWaveSpeedXY calculateWaveSpeedZX
 *
 * \param v0 Flow in first direction
 * \param v1 Flow in second direction
 * \param vA Alfven speed
 * \param vS Sound speed
 * \param vW Whistler speed
 */
Real calculateCflSpeed(
   const Real& v0,
   const Real& v1,
   const Real& vA,
   const Real& vS,
   const Real& vW
) {
   const Real v = sqrt(v0*v0 + v1*v1);
   const Real vMS = sqrt(vA*vA + vS*vS);
   return max(v + vMS, v + vW);
}

/*! \brief Low-level helper function.
 * 
 * Computes the magnetosonic speed in the YZ plane. Used in upwinding the electric field X component.
 * 
 * Selects the RHO/RHO_DT2 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * 
 * If fields are not propagated, returns 0.0 as there is no information propagating.
 * 
 * \param cp Curent cell's parameters
 * \param derivs Curent cell's derivatives
 * \param nbr_cp Neighbor cell's parameters
 * \param nbr_derivs Neighbor cell's derivatives
 * \param By Current cell's By
 * \param Bz Current cell's Bz
 * \param dBydx dBydx derivative
 * \param dBydz dBydz derivative
 * \param dBzdx dBzdx derivative
 * \param dBzdy dBzdy derivative
 * \param ydir +1 or -1 depending on the interpolation direction in y
 * \param zdir +1 or -1 depending on the interpolation direction in z
 * \param minRhom Minimum mass density allowed from the neighborhood
 * \param maxRhom Maximum mass density allowed from the neighborhood
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param ret_vA Alfven speed returned
 * \param ret_vS Sound speed returned
 * \param ret_vW Whistler speed returned
 */
void calculateWaveSpeedYZ(
   const Real* cp,
   const Real* derivs,
   const Real* nbr_cp,
   const Real* nbr_derivs,
   const Real& By,
   const Real& Bz,
   const Real& dBydx,
   const Real& dBydz,
   const Real& dBzdx,
   const Real& dBzdy,
   const Real& ydir,
   const Real& zdir,
   const Real& minRhom,
   const Real& maxRhom,
   cint& RKCase,
   Real& ret_vA,
   Real& ret_vS,
   Real& ret_vW
) {
   if (Parameters::propagateField == false) {
      return;
   }

   Real A_0, A_X, rhom, p11, p22, p33;
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      A_0  = HALF*(nbr_cp[CellParams::PERBX] + nbr_cp[CellParams::BGBX] + cp[CellParams::PERBX] + cp[CellParams::BGBX]);
      A_X  = (nbr_cp[CellParams::PERBX] + nbr_cp[CellParams::BGBX]) - (cp[CellParams::PERBX] + cp[CellParams::BGBX]);
      rhom = cp[CellParams::RHOM] + ydir*HALF*derivs[fs::drhomdy] + zdir*HALF*derivs[fs::drhomdz];
      p11 = cp[CellParams::P_11] + ydir*HALF*derivs[fs::dp11dy] + zdir*HALF*derivs[fs::dp11dz];
      p22 = cp[CellParams::P_22] + ydir*HALF*derivs[fs::dp22dy] + zdir*HALF*derivs[fs::dp22dz];
      p33 = cp[CellParams::P_33] + ydir*HALF*derivs[fs::dp33dy] + zdir*HALF*derivs[fs::dp33dz];
   } else { // RKCase == RK_ORDER2_STEP1
      A_0  = HALF*(nbr_cp[CellParams::PERBX_DT2] + nbr_cp[CellParams::BGBX] + cp[CellParams::PERBX_DT2] + cp[CellParams::BGBX]);
      A_X  = (nbr_cp[CellParams::PERBX_DT2] + nbr_cp[CellParams::BGBX]) - (cp[CellParams::PERBX_DT2] + cp[CellParams::BGBX]);
      rhom = cp[CellParams::RHOM_DT2] + ydir*HALF*derivs[fs::drhomdy] + zdir*HALF*derivs[fs::drhomdz];
      p11 = cp[CellParams::P_11_DT2] + ydir*HALF*derivs[fs::dp11dy] + zdir*HALF*derivs[fs::dp11dz];
      p22 = cp[CellParams::P_22_DT2] + ydir*HALF*derivs[fs::dp22dy] + zdir*HALF*derivs[fs::dp22dz];
      p33 = cp[CellParams::P_33_DT2] + ydir*HALF*derivs[fs::dp33dy] + zdir*HALF*derivs[fs::dp33dz];
   }
   if (rhom < minRhom) {
      rhom = minRhom;
   } else if (rhom > maxRhom) {
      rhom = maxRhom;
   }
   
   const Real A_Y  = nbr_derivs[fs::dPERBxdy] + nbr_derivs[fs::dBGBxdy] + derivs[fs::dPERBxdy] + derivs[fs::dBGBxdy];
   const Real A_XY = nbr_derivs[fs::dPERBxdy] + nbr_derivs[fs::dBGBxdy] - (derivs[fs::dPERBxdy] + derivs[fs::dBGBxdy]);
   const Real A_Z  = nbr_derivs[fs::dPERBxdz] + nbr_derivs[fs::dBGBxdz] + derivs[fs::dPERBxdz] + derivs[fs::dBGBxdz];
   const Real A_XZ = nbr_derivs[fs::dPERBxdz] + nbr_derivs[fs::dBGBxdz] - (derivs[fs::dPERBxdz] + derivs[fs::dBGBxdz]);

   const Real Bx2  = (A_0 + ydir*HALF*A_Y + zdir*HALF*A_Z)*(A_0 + ydir*HALF*A_Y + zdir*HALF*A_Z)
     + TWELWTH*(A_X + ydir*HALF*A_XY + zdir*HALF*A_XZ)*(A_X + ydir*HALF*A_XY + zdir*HALF*A_XZ); // OK
   const Real By2  = (By + zdir*HALF*dBydz)*(By + zdir*HALF*dBydz) + TWELWTH*dBydx*dBydx; // OK
   const Real Bz2  = (Bz + ydir*HALF*dBzdy)*(Bz + ydir*HALF*dBzdy) + TWELWTH*dBzdx*dBzdx; // OK
   
   const Real Bmag2 = Bx2 + By2 + Bz2;
   
   p11 = p11 < 0.0 ? 0.0 : p11;
   p22 = p22 < 0.0 ? 0.0 : p22;
   p33 = p33 < 0.0 ? 0.0 : p33;

   const Real vA2 = divideIfNonZero(Bmag2, pc::MU_0*rhom); // Alfven speed
   const Real vS2 = divideIfNonZero(p11+p22+p33, 2.0*rhom); // sound speed, adiabatic coefficient 3/2, P=1/3*trace in sound speed
#warning Which ion species to take into whistler speed?
   const Real vW = Parameters::ohmHallTerm > 0 ? divideIfNonZero(2.0*M_PI*vA2*pc::MASS_PROTON, cp[CellParams::DX]*pc::CHARGE*sqrt(Bmag2)) : 0.0; // whistler speed
   
   ret_vA = sqrt(vA2);
   ret_vS = sqrt(vS2);
   ret_vW = vW;
}

/*! \brief Low-level helper function.
 * 
 * Computes the magnetosonic speed in the XZ plane. Used in upwinding the electric field Y component.
 * 
 * Selects the RHO/RHO_DT2 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * 
 * If fields are not propagated, returns 0.0 as there is no information propagating.
 * 
 * \param cp Curent cell's parameters
 * \param derivs Curent cell's derivatives
 * \param nbr_cp Neighbor cell's parameters
 * \param nbr_derivs Neighbor cell's derivatives
 * \param Bx Current cell's Bx
 * \param Bz Current cell's Bz
 * \param dBxdy dBxdy derivative
 * \param dBxdz dBxdz derivative
 * \param dBzdx dBzdx derivative
 * \param dBzdy dBzdy derivative
 * \param xdir +1 or -1 depending on the interpolation direction in x
 * \param zdir +1 or -1 depending on the interpolation direction in z
 * \param minRhom Minimum mass density allowed from the neighborhood
 * \param maxRhom Maximum mass density allowed from the neighborhood
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param ret_vA Alfven speed returned
 * \param ret_vS Sound speed returned
 * \param ret_vW Whistler speed returned
 */
void calculateWaveSpeedXZ(
   const Real* cp,
   const Real* derivs,
   const Real* nbr_cp,
   const Real* nbr_derivs,
   const Real& Bx,
   const Real& Bz,
   const Real& dBxdy,
   const Real& dBxdz,
   const Real& dBzdx,
   const Real& dBzdy,
   const Real& xdir,
   const Real& zdir,
   const Real& minRhom,
   const Real& maxRhom,
   cint& RKCase,
   Real& ret_vA,
   Real& ret_vS,
   Real& ret_vW
) {
   if (Parameters::propagateField == false) {
      return;
   }

   Real B_0, B_Y, rhom, p11, p22, p33;
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      B_0  = HALF*(nbr_cp[CellParams::PERBY] + nbr_cp[CellParams::BGBY] + cp[CellParams::PERBY] + cp[CellParams::BGBY]);
      B_Y  = (nbr_cp[CellParams::PERBY] + nbr_cp[CellParams::BGBY]) - (cp[CellParams::PERBY] + cp[CellParams::BGBY]);
      rhom = cp[CellParams::RHOM] + xdir*HALF*derivs[fs::drhomdx] + zdir*HALF*derivs[fs::drhomdz];
      p11 = cp[CellParams::P_11] + xdir*HALF*derivs[fs::dp11dx] + zdir*HALF*derivs[fs::dp11dz];
      p22 = cp[CellParams::P_22] + xdir*HALF*derivs[fs::dp22dx] + zdir*HALF*derivs[fs::dp22dz];
      p33 = cp[CellParams::P_33] + xdir*HALF*derivs[fs::dp33dx] + zdir*HALF*derivs[fs::dp33dz];
   } else { // RKCase == RK_ORDER2_STEP1
      B_0  = HALF*(nbr_cp[CellParams::PERBY_DT2] + nbr_cp[CellParams::BGBY] + cp[CellParams::PERBY_DT2] + cp[CellParams::BGBY]);
      B_Y  = (nbr_cp[CellParams::PERBY_DT2] + nbr_cp[CellParams::BGBY]) - (cp[CellParams::PERBY_DT2] + cp[CellParams::BGBY]);
      rhom = cp[CellParams::RHOM_DT2] + xdir*HALF*derivs[fs::drhomdx] + zdir*HALF*derivs[fs::drhomdz];
      p11 = cp[CellParams::P_11_DT2] + xdir*HALF*derivs[fs::dp11dx] + zdir*HALF*derivs[fs::dp11dz];
      p22 = cp[CellParams::P_22_DT2] + xdir*HALF*derivs[fs::dp22dx] + zdir*HALF*derivs[fs::dp22dz];
      p33 = cp[CellParams::P_33_DT2] + xdir*HALF*derivs[fs::dp33dx] + zdir*HALF*derivs[fs::dp33dz];
   }
   if (rhom < minRhom) {
      rhom = minRhom;
   } else if (rhom > maxRhom) {
      rhom = maxRhom;
   }
   
   const Real B_X  = nbr_derivs[fs::dPERBydx] + nbr_derivs[fs::dBGBydx] + derivs[fs::dPERBydx] + derivs[fs::dBGBydx];
   const Real B_XY = nbr_derivs[fs::dPERBydx] + nbr_derivs[fs::dBGBydx] - (derivs[fs::dPERBydx] + derivs[fs::dBGBydx]);
   const Real B_Z  = nbr_derivs[fs::dPERBydz] + nbr_derivs[fs::dBGBydz] + derivs[fs::dPERBydz] + derivs[fs::dBGBydz];
   const Real B_YZ = nbr_derivs[fs::dPERBydz] + nbr_derivs[fs::dBGBydz] - (derivs[fs::dPERBydz] + derivs[fs::dBGBydz]);
      
   const Real By2  = (B_0 + xdir*HALF*B_X + zdir*HALF*B_Z)*(B_0 + xdir*HALF*B_X + zdir*HALF*B_Z)
     + TWELWTH*(B_Y + xdir*HALF*B_XY + zdir*HALF*B_YZ)*(B_Y + xdir*HALF*B_XY + zdir*HALF*B_YZ); // OK
   const Real Bx2  = (Bx + zdir*HALF*dBxdz)*(Bx + zdir*HALF*dBxdz) + TWELWTH*dBxdy*dBxdy; // OK
   const Real Bz2  = (Bz + xdir*HALF*dBzdx)*(Bz + xdir*HALF*dBzdx) + TWELWTH*dBzdy*dBzdy; // OK
   
   const Real Bmag2 = Bx2 + By2 + Bz2;
   
   p11 = p11 < 0.0 ? 0.0 : p11;
   p22 = p22 < 0.0 ? 0.0 : p22;
   p33 = p33 < 0.0 ? 0.0 : p33;
   
   const Real vA2 = divideIfNonZero(Bmag2, pc::MU_0*rhom); // Alfven speed
   const Real vS2 = divideIfNonZero(p11+p22+p33, 2.0*rhom); // sound speed, adiabatic coefficient 3/2, P=1/3*trace in sound speed
#warning Which ion species to take into whistler speed?
   const Real vW = Parameters::ohmHallTerm > 0 ? divideIfNonZero(2.0*M_PI*vA2*pc::MASS_PROTON, cp[CellParams::DX]*pc::CHARGE*sqrt(Bmag2)) : 0.0; // whistler speed
   
   ret_vA = sqrt(vA2);
   ret_vS = sqrt(vS2);
   ret_vW = vW;
}

/*! \brief Low-level helper function.
 * 
 * Computes the magnetosonic speed in the XY plane. Used in upwinding the electric field Z component.
 * 
 * Selects the RHO/RHO_DT2 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * 
 * If fields are not propagated, returns 0.0 as there is no information propagating.
 * 
 * \param cp Curent cell's parameters
 * \param derivs Curent cell's derivatives
 * \param nbr_cp Neighbor cell's parameters
 * \param nbr_derivs Neighbor cell's derivatives
 * \param Bx Current cell's Bx
 * \param By Current cell's By
 * \param dBxdy dBxdy derivative
 * \param dBxdz dBxdz derivative
 * \param dBydx dBydx derivative
 * \param dBydz dBydz derivative
 * \param xdir +1 or -1 depending on the interpolation direction in x
 * \param ydir +1 or -1 depending on the interpolation direction in y
 * \param minRhom Minimum mass density allowed from the neighborhood
 * \param maxRhom Maximum mass density allowed from the neighborhood
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param ret_vA Alfven speed returned
 * \param ret_vS Sound speed returned
 * \param ret_vW Whistler speed returned
 */
void calculateWaveSpeedXY(
   const Real* cp,
   const Real* derivs,
   const Real* nbr_cp,
   const Real* nbr_derivs,
   const Real& Bx,
   const Real& By,
   const Real& dBxdy,
   const Real& dBxdz,
   const Real& dBydx,
   const Real& dBydz,
   const Real& xdir,
   const Real& ydir,
   const Real& minRhom,
   const Real& maxRhom,
   cint& RKCase,
   Real& ret_vA,
   Real& ret_vS,
   Real& ret_vW
) {
   if (Parameters::propagateField == false) {
      return;
   }

   Real C_0, C_Z, rhom, p11, p22, p33;
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      C_0  = HALF*(nbr_cp[CellParams::PERBZ] + nbr_cp[CellParams::BGBZ] + cp[CellParams::PERBZ] + cp[CellParams::BGBZ]);
      C_Z  = (nbr_cp[CellParams::PERBZ] + nbr_cp[CellParams::BGBZ]) - (cp[CellParams::PERBZ] + cp[CellParams::BGBZ]);
      rhom = cp[CellParams::RHOM] + xdir*HALF*derivs[fs::drhomdx] + ydir*HALF*derivs[fs::drhomdy];
      p11 = cp[CellParams::P_11] + xdir*HALF*derivs[fs::dp11dx] + ydir*HALF*derivs[fs::dp11dy];
      p22 = cp[CellParams::P_22] + xdir*HALF*derivs[fs::dp22dx] + ydir*HALF*derivs[fs::dp22dy];
      p33 = cp[CellParams::P_33] + xdir*HALF*derivs[fs::dp33dx] + ydir*HALF*derivs[fs::dp33dy];
   } else { // RKCase == RK_ORDER2_STEP1
      C_0  = HALF*(nbr_cp[CellParams::PERBZ_DT2] + nbr_cp[CellParams::BGBZ] + cp[CellParams::PERBZ_DT2] + cp[CellParams::BGBZ]);
      C_Z  = (nbr_cp[CellParams::PERBZ_DT2] + nbr_cp[CellParams::BGBZ]) - (cp[CellParams::PERBZ_DT2] + cp[CellParams::BGBZ]);
      rhom = cp[CellParams::RHOM_DT2] + xdir*HALF*derivs[fs::drhomdx] + ydir*HALF*derivs[fs::drhomdy];
      p11 = cp[CellParams::P_11_DT2] + xdir*HALF*derivs[fs::dp11dx] + ydir*HALF*derivs[fs::dp11dy];
      p22 = cp[CellParams::P_22_DT2] + xdir*HALF*derivs[fs::dp22dx] + ydir*HALF*derivs[fs::dp22dy];
      p33 = cp[CellParams::P_33_DT2] + xdir*HALF*derivs[fs::dp33dx] + ydir*HALF*derivs[fs::dp33dy];
   }
   if (rhom < minRhom) {
      rhom = minRhom;
   } else if (rhom > maxRhom) {
      rhom = maxRhom;
   }
   
   const Real C_X  = nbr_derivs[fs::dPERBzdx] + nbr_derivs[fs::dBGBzdx] + derivs[fs::dPERBzdx] + derivs[fs::dBGBzdx];
   const Real C_XZ = nbr_derivs[fs::dPERBzdx] + nbr_derivs[fs::dBGBzdx] - (derivs[fs::dPERBzdx] + derivs[fs::dBGBzdx]);
   const Real C_Y  = nbr_derivs[fs::dPERBzdy] + nbr_derivs[fs::dBGBzdy] + derivs[fs::dPERBzdy] + derivs[fs::dBGBzdy];
   const Real C_YZ = nbr_derivs[fs::dPERBzdy] + nbr_derivs[fs::dBGBzdy] - (derivs[fs::dPERBzdy] + derivs[fs::dBGBzdy]);
   
   const Real Bz2  = (C_0 + xdir*HALF*C_X + ydir*HALF*C_Y)*(C_0 + xdir*HALF*C_X + ydir*HALF*C_Y)
     + TWELWTH*(C_Z + xdir*HALF*C_XZ + ydir*HALF*C_YZ)*(C_Z + xdir*HALF*C_XZ + ydir*HALF*C_YZ);
   const Real Bx2  = (Bx + ydir*HALF*dBxdy)*(Bx + ydir*HALF*dBxdy) + TWELWTH*dBxdz*dBxdz;
   const Real By2  = (By + xdir*HALF*dBydx)*(By + xdir*HALF*dBydx) + TWELWTH*dBydz*dBydz;
   
   const Real Bmag2 = Bx2 + By2 + Bz2;
   
   p11 = p11 < 0.0 ? 0.0 : p11;
   p22 = p22 < 0.0 ? 0.0 : p22;
   p33 = p33 < 0.0 ? 0.0 : p33;
      
   const Real vA2 = divideIfNonZero(Bmag2, pc::MU_0*rhom); // Alfven speed
   const Real vS2 = divideIfNonZero(p11+p22+p33, 2.0*rhom); // sound speed, adiabatic coefficient 3/2, P=1/3*trace in sound speed
#warning Which ion species to take into whistler speed?
   const Real vW = Parameters::ohmHallTerm > 0 ? divideIfNonZero(2.0*M_PI*vA2*pc::MASS_PROTON, cp[CellParams::DX]*pc::CHARGE*sqrt(Bmag2)) : 0.0; // whistler speed
   
   ret_vA = sqrt(vA2);
   ret_vS = sqrt(vS2);
   ret_vW = vW;
}

/*! \brief Low-level electric field propagation function.
 * 
 * Computes the upwinded electric field X component along the cell's corresponding edge as the cross product of B and V in the YZ plane. Also includes the calculation of the maximally allowed time step.
 * 
 * Selects the RHO/RHO_DT2 RHOV[XYZ]/RHOV[XYZ]1 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * 
 * Note that the background B field is excluded from the diffusive term calculations because they are equivalent to a current term and the background field is curl-free.
 * 
 * \param cache Field solver cell cache
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
void calculateEdgeElectricFieldX(
   fs_cache::CellCache& cache,
   cint& RKCase
) {
   #ifdef DEBUG_FSOLVER
   bool ok = true;
   if (cache.cells[fs_cache::calculateNbrID(1  ,1  ,1  )] == NULL) ok = false;
   if (cache.cells[fs_cache::calculateNbrID(1  ,1-1,1  )] == NULL) ok = false;
   if (cache.cells[fs_cache::calculateNbrID(1  ,1-1,1-1)] == NULL) ok = false;
   if (cache.cells[fs_cache::calculateNbrID(1  ,1  ,1-1)] == NULL) ok = false;
   if (ok == false) {
      cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << std::endl;
      exit(1);
   }
   #endif

   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.
   Real ay_pos,ay_neg;              // Max. characteristic velocities to y-direction
   Real az_pos,az_neg;              // Max. characteristic velocities to z-direction
   Real Vy0,Vz0;                    // Reconstructed V
   Real vA, vS, vW;                 // Alfven, sound, whistler speed
   Real maxV = 0.0;                 // Max velocity for CFL purposes
   Real c_y, c_z;                   // Wave speeds to yz-directions

   // Get read-only pointers to NE,NW,SE,SW states (SW is rw, result is written there):
   Real*  const cp_SW = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
   creal* const cp_SE = cache.cells[fs_cache::calculateNbrID(1  ,1-1,1  )]->parameters;
   creal* const cp_NE = cache.cells[fs_cache::calculateNbrID(1  ,1-1,1-1)]->parameters;
   creal* const cp_NW = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1-1)]->parameters;
   creal* const derivs_SW = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->derivatives;
   creal* const derivs_SE = cache.cells[fs_cache::calculateNbrID(1  ,1-1,1  )]->derivatives;
   creal* const derivs_NE = cache.cells[fs_cache::calculateNbrID(1  ,1-1,1-1)]->derivatives;
   creal* const derivs_NW = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1-1)]->derivatives;

   Real By_S, Bz_W, Bz_E, By_N, perBy_S, perBz_W, perBz_E, perBy_N, rhoq_S;
   Real minRhom = std::numeric_limits<Real>::max();
   Real maxRhom = std::numeric_limits<Real>::min();
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      By_S = cp_SW[CellParams::PERBY]+cp_SW[CellParams::BGBY];
      Bz_W = cp_SW[CellParams::PERBZ]+cp_SW[CellParams::BGBZ];
      Bz_E = cp_SE[CellParams::PERBZ]+cp_SE[CellParams::BGBZ];
      By_N = cp_NW[CellParams::PERBY]+cp_NW[CellParams::BGBY];
      perBy_S = cp_SW[CellParams::PERBY];
      perBz_W = cp_SW[CellParams::PERBZ];
      perBz_E = cp_SE[CellParams::PERBZ];
      perBy_N = cp_NW[CellParams::PERBY];
      Vy0  = divideIfNonZero(cp_SW[CellParams::RHOMVY], cp_SW[CellParams::RHOM]);
      Vz0  = divideIfNonZero(cp_SW[CellParams::RHOMVZ], cp_SW[CellParams::RHOM]);
      rhoq_S = FOURTH*(cp_SW[CellParams::RHOQ] + cp_SE[CellParams::RHOQ] + cp_NW[CellParams::RHOQ] + cp_NE[CellParams::RHOQ]);
      minRhom = min(minRhom,
                  min(cp_SW[CellParams::RHOM],
                     min(cp_SE[CellParams::RHOM],
                        min(cp_NW[CellParams::RHOM],
                            cp_NE[CellParams::RHOM])
                        )
                     )
                  );
      maxRhom = max(maxRhom,
                  max(cp_SW[CellParams::RHOM],
                     max(cp_SE[CellParams::RHOM],
                        max(cp_NW[CellParams::RHOM],
                            cp_NE[CellParams::RHOM])
                        )
                     )
                  );
   } else { // RKCase == RK_ORDER2_STEP1
      By_S = cp_SW[CellParams::PERBY_DT2]+cp_SW[CellParams::BGBY];
      Bz_W = cp_SW[CellParams::PERBZ_DT2]+cp_SW[CellParams::BGBZ];
      Bz_E = cp_SE[CellParams::PERBZ_DT2]+cp_SE[CellParams::BGBZ];
      By_N = cp_NW[CellParams::PERBY_DT2]+cp_NW[CellParams::BGBY];
      perBy_S = cp_SW[CellParams::PERBY_DT2];
      perBz_W = cp_SW[CellParams::PERBZ_DT2];
      perBz_E = cp_SE[CellParams::PERBZ_DT2];
      perBy_N = cp_NW[CellParams::PERBY_DT2];
      Vy0  = divideIfNonZero(cp_SW[CellParams::RHOMVY_DT2], cp_SW[CellParams::RHOM_DT2]);
      Vz0  = divideIfNonZero(cp_SW[CellParams::RHOMVZ_DT2], cp_SW[CellParams::RHOM_DT2]);
      rhoq_S = FOURTH*(cp_SW[CellParams::RHOQ_DT2] + cp_SE[CellParams::RHOQ_DT2] + cp_NW[CellParams::RHOQ_DT2] + cp_NE[CellParams::RHOQ_DT2]);
      minRhom = min(minRhom,
                  min(cp_SW[CellParams::RHOM_DT2],
                     min(cp_SE[CellParams::RHOM_DT2],
                        min(cp_NW[CellParams::RHOM_DT2],
                            cp_NE[CellParams::RHOM_DT2])
                        )
                     )
                  );
      maxRhom = max(maxRhom,
                  max(cp_SW[CellParams::RHOM_DT2],
                     max(cp_SE[CellParams::RHOM_DT2],
                        max(cp_NW[CellParams::RHOM_DT2],
                            cp_NE[CellParams::RHOM_DT2])
                        )
                     )
                  );
   }
   rhoq_S = (rhoq_S <= Parameters::hallMinimumRhoq) ? Parameters::hallMinimumRhoq : rhoq_S ;
   
   creal dBydx_S = derivs_SW[fs::dPERBydx] + derivs_SW[fs::dBGBydx];
   creal dBydz_S = derivs_SW[fs::dPERBydz] + derivs_SW[fs::dBGBydz];
   creal dBzdx_W = derivs_SW[fs::dPERBzdx] + derivs_SW[fs::dBGBzdx];
   creal dBzdy_W = derivs_SW[fs::dPERBzdy] + derivs_SW[fs::dBGBzdy];
   creal dBzdx_E = derivs_SE[fs::dPERBzdx] + derivs_SE[fs::dBGBzdx];
   creal dBzdy_E = derivs_SE[fs::dPERBzdy] + derivs_SE[fs::dBGBzdy];
   creal dBydx_N = derivs_NW[fs::dPERBydx] + derivs_NW[fs::dBGBydx];
   creal dBydz_N = derivs_NW[fs::dPERBydz] + derivs_NW[fs::dBGBydz];
   creal dperBydz_S = derivs_SW[fs::dPERBydz];
   creal dperBydz_N = derivs_NW[fs::dPERBydz];
   creal dperBzdy_W = derivs_SW[fs::dPERBzdy];
   creal dperBzdy_E = derivs_SE[fs::dPERBzdy];

   // Ex and characteristic speeds on this cell:
   // 1st order terms:
   Real Ex_SW = By_S*Vz0 - Bz_W*Vy0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ex_SW += Parameters::resistivity *
     sqrt((cp_SW[CellParams::BGBX]+cp_SW[CellParams::PERBX])*
          (cp_SW[CellParams::BGBX]+cp_SW[CellParams::PERBX]) +
          (cp_SW[CellParams::BGBY]+cp_SW[CellParams::PERBY])*
          (cp_SW[CellParams::BGBY]+cp_SW[CellParams::PERBY]) +
          (cp_SW[CellParams::BGBZ]+cp_SW[CellParams::PERBZ])*
          (cp_SW[CellParams::BGBZ]+cp_SW[CellParams::PERBZ])
         ) /
     cp_SW[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_SW[fs::dPERBzdy]/cp_SW[CellParams::DY] - derivs_SW[fs::dPERBydz]/cp_SW[CellParams::DZ]);
   
   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ex_SW += cp_SW[CellParams::EXHALL_000_100] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ex_SW += cp_SW[CellParams::EXGRADPE];
   }

   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_SW += +HALF*((By_S - HALF*dBydz_S)*(-derivs_SW[fs::dVzdy] - derivs_SW[fs::dVzdz]) - dBydz_S*Vz0 + SIXTH*dBydx_S*derivs_SW[fs::dVzdx]);
      Ex_SW += -HALF*((Bz_W - HALF*dBzdy_W)*(-derivs_SW[fs::dVydy] - derivs_SW[fs::dVydz]) - dBzdy_W*Vy0 + SIXTH*dBzdx_W*derivs_SW[fs::dVydx]);
   #endif

   creal* const nbr_cp_SW     = cache.cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->parameters;
   creal* const nbr_derivs_SW = cache.cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->derivatives;
   
   calculateWaveSpeedYZ(cp_SW, derivs_SW, nbr_cp_SW, nbr_derivs_SW, By_S, Bz_W, dBydx_S, dBydz_S, dBzdx_W, dBzdy_W, MINUS, MINUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_y = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_z = c_y;
   ay_neg   = max(ZERO,-Vy0 + c_y);
   ay_pos   = max(ZERO,+Vy0 + c_y);
   az_neg   = max(ZERO,-Vz0 + c_z);
   az_pos   = max(ZERO,+Vz0 + c_z);
   maxV = max(maxV, calculateCflSpeed(Vy0, Vz0, vA, vS, vW));

   // Ex and characteristic speeds on j-1 neighbour:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vy0  = divideIfNonZero(cp_SE[CellParams::RHOMVY], cp_SE[CellParams::RHOM]);
      Vz0  = divideIfNonZero(cp_SE[CellParams::RHOMVZ], cp_SE[CellParams::RHOM]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vy0  = divideIfNonZero(cp_SE[CellParams::RHOMVY_DT2], cp_SE[CellParams::RHOM_DT2]);
      Vz0  = divideIfNonZero(cp_SE[CellParams::RHOMVZ_DT2], cp_SE[CellParams::RHOM_DT2]);
   }

   // 1st order terms:
   Real Ex_SE = By_S*Vz0 - Bz_E*Vy0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ex_SE += Parameters::resistivity *
     sqrt((cp_SE[CellParams::BGBX]+cp_SE[CellParams::PERBX])*
          (cp_SE[CellParams::BGBX]+cp_SE[CellParams::PERBX]) +
          (cp_SE[CellParams::BGBY]+cp_SE[CellParams::PERBY])*
          (cp_SE[CellParams::BGBY]+cp_SE[CellParams::PERBY]) +
          (cp_SE[CellParams::BGBZ]+cp_SE[CellParams::PERBZ])*
          (cp_SE[CellParams::BGBZ]+cp_SE[CellParams::PERBZ])
         ) /
     cp_SE[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_SE[fs::dPERBzdy]/cp_SE[CellParams::DY] - derivs_SE[fs::dPERBydz]/cp_SE[CellParams::DZ]);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ex_SE += cp_SE[CellParams::EXHALL_010_110] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ex_SE += cp_SE[CellParams::EXGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_SE += +HALF*((By_S - HALF*dBydz_S)*(+derivs_SE[fs::dVzdy] - derivs_SE[fs::dVzdz]) - dBydz_S*Vz0 + SIXTH*dBydx_S*derivs_SE[fs::dVzdx]);
      Ex_SE += -HALF*((Bz_E + HALF*dBzdy_E)*(+derivs_SE[fs::dVydy] - derivs_SE[fs::dVydz]) + dBzdy_E*Vy0 + SIXTH*dBzdx_E*derivs_SE[fs::dVydx]);
   #endif

   creal* const nbr_cp_SE     = cache.cells[fs_cache::calculateNbrID(1+1,1-1,1  )]->parameters;
   creal* const nbr_derivs_SE = cache.cells[fs_cache::calculateNbrID(1+1,1-1,1  )]->derivatives;
   
   calculateWaveSpeedYZ(cp_SE, derivs_SE, nbr_cp_SE, nbr_derivs_SE, By_S, Bz_E, dBydx_S, dBydz_S, dBzdx_E, dBzdy_E, PLUS, MINUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_y = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_z = c_y;
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   maxV = max(maxV, calculateCflSpeed(Vy0, Vz0, vA, vS, vW));

   // Ex and characteristic speeds on k-1 neighbour:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vy0  = divideIfNonZero(cp_NW[CellParams::RHOMVY], cp_NW[CellParams::RHOM]);
      Vz0  = divideIfNonZero(cp_NW[CellParams::RHOMVZ], cp_NW[CellParams::RHOM]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vy0  = divideIfNonZero(cp_NW[CellParams::RHOMVY_DT2], cp_NW[CellParams::RHOM_DT2]);
      Vz0  = divideIfNonZero(cp_NW[CellParams::RHOMVZ_DT2], cp_NW[CellParams::RHOM_DT2]);
   }

   // 1st order terms:
   Real Ex_NW    = By_N*Vz0 - Bz_W*Vy0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ex_NW += Parameters::resistivity *
     sqrt((cp_NW[CellParams::BGBX]+cp_NW[CellParams::PERBX])*
          (cp_NW[CellParams::BGBX]+cp_NW[CellParams::PERBX]) +
          (cp_NW[CellParams::BGBY]+cp_NW[CellParams::PERBY])*
          (cp_NW[CellParams::BGBY]+cp_NW[CellParams::PERBY]) +
          (cp_NW[CellParams::BGBZ]+cp_NW[CellParams::PERBZ])*
          (cp_NW[CellParams::BGBZ]+cp_NW[CellParams::PERBZ])
         ) /
     cp_NW[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_NW[fs::dPERBzdy]/cp_NW[CellParams::DY] - derivs_NW[fs::dPERBydz]/cp_NW[CellParams::DZ]);
   
   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ex_NW += cp_NW[CellParams::EXHALL_001_101]  / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ex_NW += cp_NW[CellParams::EXGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_NW += +HALF*((By_N + HALF*dBydz_N)*(-derivs_NW[fs::dVzdy] + derivs_NW[fs::dVzdz]) + dBydz_N*Vz0 + SIXTH*dBydx_N*derivs_NW[fs::dVzdx]);
      Ex_NW += -HALF*((Bz_W - HALF*dBzdy_W)*(-derivs_NW[fs::dVydy] + derivs_NW[fs::dVydz]) - dBzdy_W*Vy0 + SIXTH*dBzdx_W*derivs_NW[fs::dVydx]);
   #endif
   
   creal* const nbr_cp_NW     = cache.cells[fs_cache::calculateNbrID(1+1,1  ,1-1)]->parameters;
   creal* const nbr_derivs_NW = cache.cells[fs_cache::calculateNbrID(1+1,1  ,1-1)]->derivatives;
   
   calculateWaveSpeedYZ(cp_NW, derivs_NW, nbr_cp_NW, nbr_derivs_NW, By_N, Bz_W, dBydx_N, dBydz_N, dBzdx_W, dBzdy_W, MINUS, PLUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_y = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_z = c_y;
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   maxV = max(maxV, calculateCflSpeed(Vy0, Vz0, vA, vS, vW));

   // Ex and characteristic speeds on j-1,k-1 neighbour:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vy0 = divideIfNonZero(cp_NE[CellParams::RHOMVY], cp_NE[CellParams::RHOM]);
      Vz0 = divideIfNonZero(cp_NE[CellParams::RHOMVZ], cp_NE[CellParams::RHOM]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vy0 = divideIfNonZero(cp_NE[CellParams::RHOMVY_DT2], cp_NE[CellParams::RHOM_DT2]);
      Vz0 = divideIfNonZero(cp_NE[CellParams::RHOMVZ_DT2], cp_NE[CellParams::RHOM_DT2]);
   }
   
   // 1st order terms:
   Real Ex_NE    = By_N*Vz0 - Bz_E*Vy0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ex_NE += Parameters::resistivity *
            sqrt((cp_NE[CellParams::BGBX]+cp_NE[CellParams::PERBX])*
                 (cp_NE[CellParams::BGBX]+cp_NE[CellParams::PERBX]) +
                 (cp_NE[CellParams::BGBY]+cp_NE[CellParams::PERBY])*
                 (cp_NE[CellParams::BGBY]+cp_NE[CellParams::PERBY]) +
                 (cp_NE[CellParams::BGBZ]+cp_NE[CellParams::PERBZ])*
                 (cp_NE[CellParams::BGBZ]+cp_NE[CellParams::PERBZ])
                ) /
               cp_NE[CellParams::RHOQ] /
            physicalconstants::MU_0 *
            (derivs_NE[fs::dPERBzdy]/cp_NE[CellParams::DY] - derivs_NE[fs::dPERBydz]/cp_NE[CellParams::DZ]);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ex_NE += cp_NE[CellParams::EXHALL_011_111] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ex_NE += cp_NE[CellParams::EXGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_NE += +HALF*((By_N + HALF*dBydz_N)*(+derivs_NE[fs::dVzdy] + derivs_NE[fs::dVzdz]) + dBydz_N*Vz0 + SIXTH*dBydx_N*derivs_NE[fs::dVzdx]);
      Ex_NE += -HALF*((Bz_E + HALF*dBzdy_E)*(+derivs_NE[fs::dVydy] + derivs_NE[fs::dVydz]) + dBzdy_E*Vy0 + SIXTH*dBzdx_E*derivs_NE[fs::dVydx]);
   #endif
   
   creal* const nbr_cp_NE     = cache.cells[fs_cache::calculateNbrID(1+1,1-1,1-1)]->parameters;
   creal* const nbr_derivs_NE = cache.cells[fs_cache::calculateNbrID(1+1,1-1,1-1)]->derivatives;
   
   calculateWaveSpeedYZ(cp_NE, derivs_NE, nbr_cp_NE, nbr_derivs_NE, By_N, Bz_E, dBydx_N, dBydz_N, dBzdx_E, dBzdy_E, PLUS, PLUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_y = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_z = c_y;
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   maxV = max(maxV, calculateCflSpeed(Vy0, Vz0, vA, vS, vW));

   // Calculate properly upwinded edge-averaged Ex:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      cp_SW[CellParams::EX]  = ay_pos*az_pos*Ex_NE + ay_pos*az_neg*Ex_SE + ay_neg*az_pos*Ex_NW + ay_neg*az_neg*Ex_SW;
      cp_SW[CellParams::EX] /= ((ay_pos+ay_neg)*(az_pos+az_neg)+EPS);
      if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
         // 1st order diffusive terms:
         cp_SW[CellParams::EX] -= az_pos*az_neg/(az_pos+az_neg+EPS)*(perBy_S-perBy_N);
         cp_SW[CellParams::EX] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(perBz_W-perBz_E);
#else
         // 2nd     order diffusive terms
         cp_SW[CellParams::EX] -= az_pos*az_neg/(az_pos+az_neg+EPS)*((perBy_S-HALF*dperBydz_S) - (perBy_N+HALF*dperBydz_N));
         cp_SW[CellParams::EX] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((perBz_W-HALF*dperBzdy_W) - (perBz_E+HALF*dperBzdy_E));
#endif
      }
   } else { // RKCase == RK_ORDER2_STEP1
      cp_SW[CellParams::EX_DT2]  = ay_pos*az_pos*Ex_NE + ay_pos*az_neg*Ex_SE + ay_neg*az_pos*Ex_NW + ay_neg*az_neg*Ex_SW;
      cp_SW[CellParams::EX_DT2] /= ((ay_pos+ay_neg)*(az_pos+az_neg)+EPS);
      if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
         // 1st order diffusive terms:
         cp_SW[CellParams::EX_DT2] -= az_pos*az_neg/(az_pos+az_neg+EPS)*(perBy_S-perBy_N);
         cp_SW[CellParams::EX_DT2] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(perBz_W-perBz_E);
#else
         // 2nd order diffusive terms
         cp_SW[CellParams::EX_DT2] -= az_pos*az_neg/(az_pos+az_neg+EPS)*((perBy_S-HALF*dperBydz_S) - (perBy_N+HALF*dperBydz_N));
         cp_SW[CellParams::EX_DT2] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((perBz_W-HALF*dperBzdy_W) - (perBz_E+HALF*dperBzdy_E));
#endif
      }
   }
   
   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      //compute maximum timestep for fieldsolver in this cell (CFL=1)
      Real min_dx=std::numeric_limits<Real>::max();
      min_dx=min(min_dx,cp_SW[CellParams::DY]);
      min_dx=min(min_dx,cp_SW[CellParams::DZ]);
      //update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (maxV != ZERO) cp_SW[CellParams::MAXFDT] = min(cp_SW[CellParams::MAXFDT],min_dx/maxV);
   }
}

/*! \brief Low-level electric field propagation function.
 * 
 * Computes the upwinded electric field Y component along the cell's corresponding edge as the cross product of B and V in the XZ plane. Also includes the calculation of the maximally allowed time step.
 * 
 * Selects the RHO/RHO_DT2 RHOV[XYZ]/RHOV[XYZ]1 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * 
 * Note that the background B field is excluded from the diffusive term calculations because they are equivalent to a current term and the background field is curl-free.
 * 
 * \param cache Field solver cell cache
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
void calculateEdgeElectricFieldY(
   fs_cache::CellCache& cache,
   cint& RKCase
) {
   #ifdef DEBUG_FSOLVER
   bool ok = true;
   if (cache.cells[fs_cache::calculateNbrID(1  ,1  ,1  )] == NULL) ok = false;
   if (cache.cells[fs_cache::calculateNbrID(1  ,1  ,1-1)] == NULL) ok = false;
   if (cache.cells[fs_cache::calculateNbrID(1-1,1  ,1  )] == NULL) ok = false;
   if (cache.cells[fs_cache::calculateNbrID(1-1,1  ,1-1)] == NULL) ok = false;
   if (ok == false) {
      cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << std::endl;
      exit(1);
   }
   #endif
   
   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.
   Real ax_pos,ax_neg;              // Max. characteristic velocities to x-direction
   Real az_pos,az_neg;              // Max. characteristic velocities to z-direction
   Real Vx0,Vz0;                    // Reconstructed V
   Real vA, vS, vW;                 // Alfven, sound, whistler speed
   Real maxV = 0.0;                 // Max velocity for CFL purposes
   Real c_x,c_z;                    // Wave speeds to xz-directions

   // Get read-only pointers to NE,NW,SE,SW states (SW is rw, result is written there):
   Real* const  cp_SW = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
   creal* const cp_SE = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1-1)]->parameters;
   creal* const cp_NW = cache.cells[fs_cache::calculateNbrID(1-1,1  ,1  )]->parameters;
   creal* const cp_NE = cache.cells[fs_cache::calculateNbrID(1-1,1  ,1-1)]->parameters;
   creal* const derivs_SW = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->derivatives;
   creal* const derivs_SE = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1-1)]->derivatives;
   creal* const derivs_NW = cache.cells[fs_cache::calculateNbrID(1-1,1  ,1  )]->derivatives;
   creal* const derivs_NE = cache.cells[fs_cache::calculateNbrID(1-1,1  ,1-1)]->derivatives;

   // Fetch required plasma parameters:
   Real Bz_S, Bx_W, Bx_E, Bz_N, perBz_S, perBx_W, perBx_E, perBz_N, rhoq_S;
   Real minRhom = std::numeric_limits<Real>::max();
   Real maxRhom = std::numeric_limits<Real>::min();
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Bz_S = cp_SW[CellParams::PERBZ]+cp_SW[CellParams::BGBZ];
      Bx_W = cp_SW[CellParams::PERBX]+cp_SW[CellParams::BGBX];
      Bx_E = cp_SE[CellParams::PERBX]+cp_SE[CellParams::BGBX];
      Bz_N = cp_NW[CellParams::PERBZ]+cp_NW[CellParams::BGBZ];
      perBz_S = cp_SW[CellParams::PERBZ];
      perBx_W = cp_SW[CellParams::PERBX];
      perBx_E = cp_SE[CellParams::PERBX];
      perBz_N = cp_NW[CellParams::PERBZ];
      Vx0  = divideIfNonZero(cp_SW[CellParams::RHOMVX], cp_SW[CellParams::RHOM]);
      Vz0  = divideIfNonZero(cp_SW[CellParams::RHOMVZ], cp_SW[CellParams::RHOM]);
      rhoq_S = FOURTH*(cp_SW[CellParams::RHOQ] + cp_SE[CellParams::RHOQ] + cp_NW[CellParams::RHOQ] + cp_NE[CellParams::RHOQ]);
      minRhom = min(minRhom,
                  min(cp_SW[CellParams::RHOM],
                     min(cp_SE[CellParams::RHOM],
                        min(cp_NW[CellParams::RHOM],
                            cp_NE[CellParams::RHOM])
                        )
                     )
                  );
      maxRhom = max(maxRhom,
                  max(cp_SW[CellParams::RHOM],
                     max(cp_SE[CellParams::RHOM],
                        max(cp_NW[CellParams::RHOM],
                            cp_NE[CellParams::RHOM])
                        )
                     )
                  );
   } else { // RKCase == RK_ORDER2_STEP1
      Bz_S = cp_SW[CellParams::PERBZ_DT2]+cp_SW[CellParams::BGBZ];
      Bx_W = cp_SW[CellParams::PERBX_DT2]+cp_SW[CellParams::BGBX];
      Bx_E = cp_SE[CellParams::PERBX_DT2]+cp_SE[CellParams::BGBX];
      Bz_N = cp_NW[CellParams::PERBZ_DT2]+cp_NW[CellParams::BGBZ];
      perBz_S = cp_SW[CellParams::PERBZ_DT2];
      perBx_W = cp_SW[CellParams::PERBX_DT2];
      perBx_E = cp_SE[CellParams::PERBX_DT2];
      perBz_N = cp_NW[CellParams::PERBZ_DT2];
      Vx0  = divideIfNonZero(cp_SW[CellParams::RHOMVX_DT2], cp_SW[CellParams::RHOM_DT2]);
      Vz0  = divideIfNonZero(cp_SW[CellParams::RHOMVZ_DT2], cp_SW[CellParams::RHOM_DT2]);
      rhoq_S = FOURTH*(cp_SW[CellParams::RHOQ_DT2] + cp_SE[CellParams::RHOQ_DT2] + cp_NW[CellParams::RHOQ_DT2] + cp_NE[CellParams::RHOQ_DT2]);
      minRhom = min(minRhom,
                  min(cp_SW[CellParams::RHOM_DT2],
                     min(cp_SE[CellParams::RHOM_DT2],
                        min(cp_NW[CellParams::RHOM_DT2],
                            cp_NE[CellParams::RHOM_DT2])
                        )
                     )
                  );
      maxRhom = max(maxRhom,
                  max(cp_SW[CellParams::RHOM_DT2],
                     max(cp_SE[CellParams::RHOM_DT2],
                        max(cp_NW[CellParams::RHOM_DT2],
                            cp_NE[CellParams::RHOM_DT2])
                        )
                     )
                  );
   }
   rhoq_S = (rhoq_S <= Parameters::hallMinimumRhoq) ? Parameters::hallMinimumRhoq : rhoq_S ;
   
   creal dBxdy_W = derivs_SW[fs::dPERBxdy] + derivs_SW[fs::dBGBxdy];
   creal dBxdz_W = derivs_SW[fs::dPERBxdz] + derivs_SW[fs::dBGBxdz];
   creal dBzdx_S = derivs_SW[fs::dPERBzdx] + derivs_SW[fs::dBGBzdx];
   creal dBzdy_S = derivs_SW[fs::dPERBzdy] + derivs_SW[fs::dBGBzdy];
   creal dBxdy_E = derivs_SE[fs::dPERBxdy] + derivs_SE[fs::dBGBxdy];
   creal dBxdz_E = derivs_SE[fs::dPERBxdz] + derivs_SE[fs::dBGBxdz];
   creal dBzdx_N = derivs_NW[fs::dPERBzdx] + derivs_NW[fs::dBGBzdx];
   creal dBzdy_N = derivs_NW[fs::dPERBzdy] + derivs_NW[fs::dBGBzdy];
   creal dperBzdx_S = derivs_SW[fs::dPERBzdx];
   creal dperBzdx_N = derivs_NW[fs::dPERBzdx];
   creal dperBxdz_W = derivs_SW[fs::dPERBxdz];
   creal dperBxdz_E = derivs_SE[fs::dPERBxdz];
   
   // Ey and characteristic speeds on this cell:
   // 1st order terms:
   Real Ey_SW  = Bz_S*Vx0 - Bx_W*Vz0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ey_SW += Parameters::resistivity *
     sqrt((cp_SW[CellParams::BGBX]+cp_SW[CellParams::PERBX])*
          (cp_SW[CellParams::BGBX]+cp_SW[CellParams::PERBX]) +
          (cp_SW[CellParams::BGBY]+cp_SW[CellParams::PERBY])*
          (cp_SW[CellParams::BGBY]+cp_SW[CellParams::PERBY]) +
          (cp_SW[CellParams::BGBZ]+cp_SW[CellParams::PERBZ])*
          (cp_SW[CellParams::BGBZ]+cp_SW[CellParams::PERBZ])
         ) /
     cp_SW[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_SW[fs::dPERBxdz]/cp_SW[CellParams::DZ] - derivs_SW[fs::dPERBzdx]/cp_SW[CellParams::DX]);
   
   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ey_SW += cp_SW[CellParams::EYHALL_000_010] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ey_SW += cp_SW[CellParams::EYGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms
      Ey_SW += +HALF*((Bz_S - HALF*dBzdx_S)*(-derivs_SW[fs::dVxdx] - derivs_SW[fs::dVxdz]) - dBzdx_S*Vx0 + SIXTH*dBzdy_S*derivs_SW[fs::dVxdy]);
      Ey_SW += -HALF*((Bx_W - HALF*dBxdz_W)*(-derivs_SW[fs::dVzdx] - derivs_SW[fs::dVzdz]) - dBxdz_W*Vz0 + SIXTH*dBxdy_W*derivs_SW[fs::dVzdy]);
   #endif
   
   creal* const nbr_cp_SW     = cache.cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->parameters;
   creal* const nbr_derivs_SW = cache.cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->derivatives;
   
   calculateWaveSpeedXZ(cp_SW, derivs_SW, nbr_cp_SW, nbr_derivs_SW, Bx_W, Bz_S, dBxdy_W, dBxdz_W, dBzdx_S, dBzdy_S, MINUS, MINUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_z = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_x = c_z;
   az_neg   = max(ZERO,-Vz0 + c_z);
   az_pos   = max(ZERO,+Vz0 + c_z);
   ax_neg   = max(ZERO,-Vx0 + c_x);
   ax_pos   = max(ZERO,+Vx0 + c_x);
   maxV = max(maxV, calculateCflSpeed(Vz0, Vx0, vA, vS, vW));

   // Ey and characteristic speeds on k-1 neighbour:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vx0  = divideIfNonZero(cp_SE[CellParams::RHOMVX], cp_SE[CellParams::RHOM]);
      Vz0  = divideIfNonZero(cp_SE[CellParams::RHOMVZ], cp_SE[CellParams::RHOM]);
   } else { //RKCase == RK_ORDER2_STEP1
      Vx0  = divideIfNonZero(cp_SE[CellParams::RHOMVX_DT2], cp_SE[CellParams::RHOM_DT2]);
      Vz0  = divideIfNonZero(cp_SE[CellParams::RHOMVZ_DT2], cp_SE[CellParams::RHOM_DT2]);
   }

   // 1st order terms:
   Real Ey_SE    = Bz_S*Vx0 - Bx_E*Vz0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ey_SE += Parameters::resistivity *
     sqrt((cp_SE[CellParams::BGBX]+cp_SE[CellParams::PERBX])*
          (cp_SE[CellParams::BGBX]+cp_SE[CellParams::PERBX]) +
          (cp_SE[CellParams::BGBY]+cp_SE[CellParams::PERBY])*
          (cp_SE[CellParams::BGBY]+cp_SE[CellParams::PERBY]) +
          (cp_SE[CellParams::BGBZ]+cp_SE[CellParams::PERBZ])*
          (cp_SE[CellParams::BGBZ]+cp_SE[CellParams::PERBZ])
         ) /
     cp_SE[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_SE[fs::dPERBxdz]/cp_SE[CellParams::DZ] - derivs_SE[fs::dPERBzdx]/cp_SE[CellParams::DX]);

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ey_SE += cp_SE[CellParams::EYHALL_001_011] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ey_SE += cp_SE[CellParams::EYGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ey_SE += +HALF*((Bz_S - HALF*dBzdx_S)*(-derivs_SE[fs::dVxdx] + derivs_SE[fs::dVxdz]) - dBzdx_S*Vx0 + SIXTH*dBzdy_S*derivs_SE[fs::dVxdy]);
      Ey_SE += -HALF*((Bx_E + HALF*dBxdz_E)*(-derivs_SE[fs::dVzdx] + derivs_SE[fs::dVzdz]) + dBxdz_E*Vz0 + SIXTH*dBxdy_E*derivs_SE[fs::dVzdy]);
   #endif
   
   creal* const nbr_cp_SE     = cache.cells[fs_cache::calculateNbrID(1  ,1+1,1-1)]->parameters;
   creal* const nbr_derivs_SE = cache.cells[fs_cache::calculateNbrID(1  ,1+1,1-1)]->derivatives;
   
   calculateWaveSpeedXZ(cp_SE, derivs_SE, nbr_cp_SE, nbr_derivs_SE, Bx_E, Bz_S, dBxdy_E, dBxdz_E, dBzdx_S, dBzdy_S, MINUS, PLUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_z = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_x = c_z;
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   maxV = max(maxV, calculateCflSpeed(Vz0, Vx0, vA, vS, vW));
   
   // Ey and characteristic speeds on i-1 neighbour:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vz0  = divideIfNonZero(cp_NW[CellParams::RHOMVZ], cp_NW[CellParams::RHOM]);
      Vx0  = divideIfNonZero(cp_NW[CellParams::RHOMVX], cp_NW[CellParams::RHOM]);
   } else { //RKCase == RK_ORDER2_STEP1
      Vz0  = divideIfNonZero(cp_NW[CellParams::RHOMVZ_DT2], cp_NW[CellParams::RHOM_DT2]);
      Vx0  = divideIfNonZero(cp_NW[CellParams::RHOMVX_DT2], cp_NW[CellParams::RHOM_DT2]);
   }
   
   // 1st order terms:
   Real Ey_NW    = Bz_N*Vx0 - Bx_W*Vz0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ey_NW += Parameters::resistivity *
     sqrt((cp_NW[CellParams::BGBX]+cp_NW[CellParams::PERBX])*
          (cp_NW[CellParams::BGBX]+cp_NW[CellParams::PERBX]) +
          (cp_NW[CellParams::BGBY]+cp_NW[CellParams::PERBY])*
          (cp_NW[CellParams::BGBY]+cp_NW[CellParams::PERBY]) +
          (cp_NW[CellParams::BGBZ]+cp_NW[CellParams::PERBZ])*
          (cp_NW[CellParams::BGBZ]+cp_NW[CellParams::PERBZ])
         ) /
     cp_NW[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_NW[fs::dPERBxdz]/cp_NW[CellParams::DZ] - derivs_NW[fs::dPERBzdx]/cp_NW[CellParams::DX]);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ey_NW += cp_NW[CellParams::EYHALL_100_110] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ey_NW += cp_NW[CellParams::EYGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ey_NW += +HALF*((Bz_N + HALF*dBzdx_N)*(+derivs_NW[fs::dVxdx] - derivs_NW[fs::dVxdz]) + dBzdx_N*Vx0 + SIXTH*dBzdy_N*derivs_NW[fs::dVxdy]);
      Ey_NW += -HALF*((Bx_W - HALF*dBxdz_W)*(+derivs_NW[fs::dVzdx] - derivs_NW[fs::dVzdz]) - dBxdz_W*Vz0 + SIXTH*dBxdy_W*derivs_NW[fs::dVzdy]);
   #endif
   
   creal* const nbr_cp_NW     = cache.cells[fs_cache::calculateNbrID(1-1,1+1,1  )]->parameters;
   creal* const nbr_derivs_NW = cache.cells[fs_cache::calculateNbrID(1-1,1+1,1  )]->derivatives;
   
   calculateWaveSpeedXZ(cp_NW, derivs_NW, nbr_cp_NW, nbr_derivs_NW, Bx_W, Bz_N, dBxdy_W, dBxdz_W, dBzdx_N, dBzdy_N, PLUS, MINUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_z = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_x = c_z;
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   maxV = max(maxV, calculateCflSpeed(Vz0, Vx0, vA, vS, vW));

   // Ey and characteristic speeds on i-1,k-1 neighbour:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vz0 = divideIfNonZero(cp_NE[CellParams::RHOMVZ], cp_NE[CellParams::RHOM]);
      Vx0 = divideIfNonZero(cp_NE[CellParams::RHOMVX], cp_NE[CellParams::RHOM]);
   } else { //RKCase == RK_ORDER2_STEP1
      Vz0 = divideIfNonZero(cp_NE[CellParams::RHOMVZ_DT2], cp_NE[CellParams::RHOM_DT2]);
      Vx0 = divideIfNonZero(cp_NE[CellParams::RHOMVX_DT2], cp_NE[CellParams::RHOM_DT2]);
   }
   
   // 1st order terms:
   Real Ey_NE    = Bz_N*Vx0 - Bx_E*Vz0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ey_NE += Parameters::resistivity *
     sqrt((cp_NE[CellParams::BGBX]+cp_NE[CellParams::PERBX])*
          (cp_NE[CellParams::BGBX]+cp_NE[CellParams::PERBX]) +
          (cp_NE[CellParams::BGBY]+cp_NE[CellParams::PERBY])*
          (cp_NE[CellParams::BGBY]+cp_NE[CellParams::PERBY]) +
          (cp_NE[CellParams::BGBZ]+cp_NE[CellParams::PERBZ])*
          (cp_NE[CellParams::BGBZ]+cp_NE[CellParams::PERBZ])
         ) /
     cp_NE[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_NE[fs::dPERBxdz]/cp_NE[CellParams::DZ] - derivs_NE[fs::dPERBzdx]/cp_NE[CellParams::DX]);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ey_NE += cp_NE[CellParams::EYHALL_101_111] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ey_NE += cp_NE[CellParams::EYGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ey_NE += +HALF*((Bz_N + HALF*dBzdx_N)*(+derivs_NE[fs::dVxdx] + derivs_NE[fs::dVxdz]) + dBzdx_N*Vx0 + SIXTH*dBzdy_N*derivs_NE[fs::dVxdy]);
      Ey_NE += -HALF*((Bx_E + HALF*dBxdz_E)*(+derivs_NE[fs::dVzdx] + derivs_NE[fs::dVzdz]) + dBxdz_E*Vz0 + SIXTH*dBxdy_E*derivs_NE[fs::dVzdy]);
   #endif

   creal* const nbr_cp_NE     = cache.cells[fs_cache::calculateNbrID(1-1,1+1,1-1)]->parameters;
   creal* const nbr_derivs_NE = cache.cells[fs_cache::calculateNbrID(1-1,1+1,1-1)]->derivatives;
   
   calculateWaveSpeedXZ(cp_NE, derivs_NE, nbr_cp_NE, nbr_derivs_NE, Bx_E, Bz_N, dBxdy_E, dBxdz_E, dBzdx_N, dBzdy_N, PLUS, PLUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_z = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_x = c_z;
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   maxV = max(maxV, calculateCflSpeed(Vz0, Vx0, vA, vS, vW));

   // Calculate properly upwinded edge-averaged Ey:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      cp_SW[CellParams::EY]  = az_pos*ax_pos*Ey_NE + az_pos*ax_neg*Ey_SE + az_neg*ax_pos*Ey_NW + az_neg*ax_neg*Ey_SW;
      cp_SW[CellParams::EY] /= ((az_pos+az_neg)*(ax_pos+ax_neg)+EPS);

      if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
         cp_SW[CellParams::EY] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(perBz_S-perBz_N);
         cp_SW[CellParams::EY] += az_pos*az_neg/(az_pos+az_neg+EPS)*(perBx_W-perBx_E);
#else
         cp_SW[CellParams::EY] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((perBz_S-HALF*dperBzdx_S) - (perBz_N+HALF*dperBzdx_N));
         cp_SW[CellParams::EY] += az_pos*az_neg/(az_pos+az_neg+EPS)*((perBx_W-HALF*dperBxdz_W) - (perBx_E+HALF*dperBxdz_E));
#endif
      }
   } else { // RKCase == RK_ORDER2_STEP1
      cp_SW[CellParams::EY_DT2]  = az_pos*ax_pos*Ey_NE + az_pos*ax_neg*Ey_SE + az_neg*ax_pos*Ey_NW + az_neg*ax_neg*Ey_SW;
      cp_SW[CellParams::EY_DT2] /= ((az_pos+az_neg)*(ax_pos+ax_neg)+EPS);
      if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
         cp_SW[CellParams::EY_DT2] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(perBz_S-perBz_N);
         cp_SW[CellParams::EY_DT2] += az_pos*az_neg/(az_pos+az_neg+EPS)*(perBx_W-perBx_E);
#else
         cp_SW[CellParams::EY_DT2] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((perBz_S-HALF*dperBzdx_S) - (perBz_N+HALF*dperBzdx_N));
         cp_SW[CellParams::EY_DT2] += az_pos*az_neg/(az_pos+az_neg+EPS)*((perBx_W-HALF*dperBxdz_W) - (perBx_E+HALF*dperBxdz_E));
#endif
      }
   }
   
   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      //compute maximum timestep for fieldsolver in this cell (CFL=1)      
      Real min_dx=std::numeric_limits<Real>::max();;
      min_dx=min(min_dx,cp_SW[CellParams::DX]);
      min_dx=min(min_dx,cp_SW[CellParams::DZ]);
      //update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (maxV!=ZERO) cp_SW[CellParams::MAXFDT]=min(cp_SW[CellParams::MAXFDT],min_dx/maxV);
   }
}

/*! \brief Low-level electric field propagation function.
 *
 * Computes the upwinded electric field Z component along the cell's corresponding edge as the cross product of B and V in the XY plane. Also includes the calculation of the maximally allowed time step.
 * 
 * Selects the RHO/RHO_DT2 RHOV[XYZ]/RHOV[XYZ]1 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * 
 * Note that the background B field is excluded from the diffusive term calculations because they are equivalent to a current term and the background field is curl-free.
 * 
 * \param cache Field solver cell cache
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
void calculateEdgeElectricFieldZ(
   fs_cache::CellCache& cache,
   cint& RKCase
) {
   #ifdef DEBUG_FSOLVER
   bool ok = true;
   if (cache.cells[fs_cache::calculateNbrID(1  ,1  ,1  )] == NULL) ok = false;
   if (cache.cells[fs_cache::calculateNbrID(1-1,1  ,1  )] == NULL) ok = false;
   if (cache.cells[fs_cache::calculateNbrID(1-1,1-1,1  )] == NULL) ok = false;
   if (cache.cells[fs_cache::calculateNbrID(1  ,1-1,1  )] == NULL) ok = false;
   if (ok == false) {
      cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << std::endl;
      exit(1);
   }
   #endif

   // An edge has four neighbouring spatial cells. Calculate 
   // electric field in each of the four cells per edge.
   Real ax_pos,ax_neg;              // Max. characteristic velocities to x-direction
   Real ay_pos,ay_neg;              // Max. characteristic velocities to y-direction
   Real Vx0,Vy0;                    // Reconstructed V
   Real vA, vS, vW;                         // Alfven, sound, whistler speed
   Real maxV = 0.0;                 // Max velocity for CFL purposes
   Real c_x,c_y;                    // Characteristic speeds to xy-directions
   
   // Get read-only pointers to NE,NW,SE,SW states (SW is rw, result is written there):
   Real* const cp_SW  = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
   creal* const cp_SE = cache.cells[fs_cache::calculateNbrID(1-1,1  ,1  )]->parameters;
   creal* const cp_NE = cache.cells[fs_cache::calculateNbrID(1-1,1-1,1  )]->parameters;
   creal* const cp_NW = cache.cells[fs_cache::calculateNbrID(1  ,1-1,1  )]->parameters;
   
   creal* const derivs_SW = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->derivatives;
   creal* const derivs_SE = cache.cells[fs_cache::calculateNbrID(1-1,1  ,1  )]->derivatives;
   creal* const derivs_NE = cache.cells[fs_cache::calculateNbrID(1-1,1-1,1  )]->derivatives;
   creal* const derivs_NW = cache.cells[fs_cache::calculateNbrID(1  ,1-1,1  )]->derivatives;

   // Fetch needed plasma parameters/derivatives from the four cells:
   Real Bx_S, By_W, By_E, Bx_N, perBx_S, perBy_W, perBy_E, perBx_N, rhoq_S;
   Real minRhom = std::numeric_limits<Real>::max();
   Real maxRhom = std::numeric_limits<Real>::min();
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Bx_S    = cp_SW[CellParams::PERBX] + cp_SW[CellParams::BGBX];
      By_W    = cp_SW[CellParams::PERBY] + cp_SW[CellParams::BGBY];
      By_E    = cp_SE[CellParams::PERBY] + cp_SE[CellParams::BGBY];
      Bx_N    = cp_NW[CellParams::PERBX] + cp_NW[CellParams::BGBX];
      perBx_S    = cp_SW[CellParams::PERBX];
      perBy_W    = cp_SW[CellParams::PERBY];
      perBy_E    = cp_SE[CellParams::PERBY];
      perBx_N    = cp_NW[CellParams::PERBX];
      Vx0  = divideIfNonZero(cp_SW[CellParams::RHOMVX], cp_SW[CellParams::RHOM]);
      Vy0  = divideIfNonZero(cp_SW[CellParams::RHOMVY], cp_SW[CellParams::RHOM]);
      rhoq_S = FOURTH*(cp_SW[CellParams::RHOQ] + cp_SE[CellParams::RHOQ] + cp_NW[CellParams::RHOQ] + cp_NE[CellParams::RHOQ]);
      minRhom = min(minRhom,
                  min(cp_SW[CellParams::RHOM],
                     min(cp_SE[CellParams::RHOM],
                        min(cp_NW[CellParams::RHOM],
                            cp_NE[CellParams::RHOM])
                        )
                     )
                  );
      maxRhom = max(maxRhom,
                  max(cp_SW[CellParams::RHOM],
                     max(cp_SE[CellParams::RHOM],
                        max(cp_NW[CellParams::RHOM],
                            cp_NE[CellParams::RHOM])
                        )
                     )
                  );
   } else { // RKCase == RK_ORDER2_STEP1
      Bx_S    = cp_SW[CellParams::PERBX_DT2] + cp_SW[CellParams::BGBX];
      By_W    = cp_SW[CellParams::PERBY_DT2] + cp_SW[CellParams::BGBY];
      By_E    = cp_SE[CellParams::PERBY_DT2] + cp_SE[CellParams::BGBY];
      Bx_N    = cp_NW[CellParams::PERBX_DT2] + cp_NW[CellParams::BGBX];
      perBx_S    = cp_SW[CellParams::PERBX_DT2];
      perBy_W    = cp_SW[CellParams::PERBY_DT2];
      perBy_E    = cp_SE[CellParams::PERBY_DT2];
      perBx_N    = cp_NW[CellParams::PERBX_DT2];
      Vx0  = divideIfNonZero(cp_SW[CellParams::RHOMVX_DT2], cp_SW[CellParams::RHOM_DT2]);
      Vy0  = divideIfNonZero(cp_SW[CellParams::RHOMVY_DT2], cp_SW[CellParams::RHOM_DT2]);
      rhoq_S = FOURTH*(cp_SW[CellParams::RHOQ_DT2] + cp_SE[CellParams::RHOQ_DT2] + cp_NW[CellParams::RHOQ_DT2] + cp_NE[CellParams::RHOQ_DT2]);
      minRhom = min(minRhom,
                  min(cp_SW[CellParams::RHOM_DT2],
                     min(cp_SE[CellParams::RHOM_DT2],
                        min(cp_NW[CellParams::RHOM_DT2],
                            cp_NE[CellParams::RHOM_DT2])
                        )
                     )
                  );
      maxRhom = max(maxRhom,
                  max(cp_SW[CellParams::RHOM_DT2],
                     max(cp_SE[CellParams::RHOM_DT2],
                        max(cp_NW[CellParams::RHOM_DT2],
                            cp_NE[CellParams::RHOM_DT2])
                        )
                     )
                  );
   }
   rhoq_S = (rhoq_S <= Parameters::hallMinimumRhoq) ? Parameters::hallMinimumRhoq : rhoq_S ;
   
   creal dBxdy_S = derivs_SW[fs::dPERBxdy] + derivs_SW[fs::dBGBxdy];
   creal dBxdz_S = derivs_SW[fs::dPERBxdz] + derivs_SW[fs::dBGBxdz];
   creal dBydx_W = derivs_SW[fs::dPERBydx] + derivs_SW[fs::dBGBydx];
   creal dBydz_W = derivs_SW[fs::dPERBydz] + derivs_SW[fs::dBGBydz];
   creal dBydx_E = derivs_SE[fs::dPERBydx] + derivs_SE[fs::dBGBydx];
   creal dBydz_E = derivs_SE[fs::dPERBydz] + derivs_SE[fs::dBGBydz];
   creal dBxdy_N = derivs_NW[fs::dPERBxdy] + derivs_NW[fs::dBGBxdy];
   creal dBxdz_N = derivs_NW[fs::dPERBxdz] + derivs_NW[fs::dBGBxdz];
   creal dperBxdy_S = derivs_SW[fs::dPERBxdy];
   creal dperBxdy_N = derivs_NW[fs::dPERBxdy];
   creal dperBydx_W = derivs_SW[fs::dPERBydx];
   creal dperBydx_E = derivs_SE[fs::dPERBydx];
   
   // Ez and characteristic speeds on SW cell:
   // 1st order terms:
   Real Ez_SW = Bx_S*Vy0 - By_W*Vx0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ez_SW += Parameters::resistivity *
     sqrt((cp_SW[CellParams::BGBX]+cp_SW[CellParams::PERBX])*
          (cp_SW[CellParams::BGBX]+cp_SW[CellParams::PERBX]) +
          (cp_SW[CellParams::BGBY]+cp_SW[CellParams::PERBY])*
          (cp_SW[CellParams::BGBY]+cp_SW[CellParams::PERBY]) +
          (cp_SW[CellParams::BGBZ]+cp_SW[CellParams::PERBZ])*
          (cp_SW[CellParams::BGBZ]+cp_SW[CellParams::PERBZ])
         ) /
     cp_SW[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_SW[fs::dPERBydx]/cp_SW[CellParams::DX] - derivs_SW[fs::dPERBxdy]/cp_SW[CellParams::DY]);
   
   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ez_SW += cp_SW[CellParams::EZHALL_000_001] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ez_SW += cp_SW[CellParams::EZGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_SW  += +HALF*((Bx_S - HALF*dBxdy_S)*(-derivs_SW[fs::dVydx] - derivs_SW[fs::dVydy]) - dBxdy_S*Vy0 + SIXTH*dBxdz_S*derivs_SW[fs::dVydz]);
      Ez_SW  += -HALF*((By_W - HALF*dBydx_W)*(-derivs_SW[fs::dVxdx] - derivs_SW[fs::dVxdy]) - dBydx_W*Vx0 + SIXTH*dBydz_W*derivs_SW[fs::dVxdz]);
   #endif
   
   // Calculate maximum wave speed (fast magnetosonic speed) on SW cell. In order 
   // to get Alfven speed we need to calculate some reconstruction coeff. for Bz:
   creal* const nbr_cp_SW     = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->parameters;
   creal* const nbr_derivs_SW = cache.cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->derivatives;
   
   calculateWaveSpeedXY(cp_SW, derivs_SW, nbr_cp_SW, nbr_derivs_SW, Bx_S, By_W, dBxdy_S, dBxdz_S, dBydx_W, dBydz_W, MINUS, MINUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_x = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_y = c_x;
   ax_neg   = max(ZERO,-Vx0 + c_x);
   ax_pos   = max(ZERO,+Vx0 + c_x);
   ay_neg   = max(ZERO,-Vy0 + c_y);
   ay_pos   = max(ZERO,+Vy0 + c_y);
   maxV = max(maxV, calculateCflSpeed(Vx0, Vy0, vA, vS, vW));

   // Ez and characteristic speeds on SE (i-1) cell:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vx0  = divideIfNonZero(cp_SE[CellParams::RHOMVX], cp_SE[CellParams::RHOM]);
      Vy0  = divideIfNonZero(cp_SE[CellParams::RHOMVY], cp_SE[CellParams::RHOM]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vx0  = divideIfNonZero(cp_SE[CellParams::RHOMVX_DT2], cp_SE[CellParams::RHOM_DT2]);
      Vy0  = divideIfNonZero(cp_SE[CellParams::RHOMVY_DT2], cp_SE[CellParams::RHOM_DT2]);
   }
   
   // 1st order terms:
   Real Ez_SE = Bx_S*Vy0 - By_E*Vx0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ez_SE += Parameters::resistivity *
     sqrt((cp_SE[CellParams::BGBX]+cp_SE[CellParams::PERBX])*
          (cp_SE[CellParams::BGBX]+cp_SE[CellParams::PERBX]) +
          (cp_SE[CellParams::BGBY]+cp_SE[CellParams::PERBY])*
          (cp_SE[CellParams::BGBY]+cp_SE[CellParams::PERBY]) +
          (cp_SE[CellParams::BGBZ]+cp_SE[CellParams::PERBZ])*
          (cp_SE[CellParams::BGBZ]+cp_SE[CellParams::PERBZ])
         ) /
     cp_SE[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_SE[fs::dPERBydx]/cp_SE[CellParams::DX] - derivs_SE[fs::dPERBxdy]/cp_SE[CellParams::DY]);
   
   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ez_SE += cp_SE[CellParams::EZHALL_100_101] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ez_SE += cp_SE[CellParams::EZGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_SE  += +HALF*((Bx_S - HALF*dBxdy_S)*(+derivs_SE[fs::dVydx] - derivs_SE[fs::dVydy]) - dBxdy_S*Vy0 + SIXTH*dBxdz_S*derivs_SE[fs::dVydz]);
      Ez_SE  += -HALF*((By_E + HALF*dBydx_E)*(+derivs_SE[fs::dVxdx] - derivs_SE[fs::dVxdy]) + dBydx_E*Vx0 + SIXTH*dBydz_E*derivs_SE[fs::dVxdz]);
   #endif
   
   creal* const nbr_cp_SE     = cache.cells[fs_cache::calculateNbrID(1-1,1  ,1+1)]->parameters;
   creal* const nbr_derivs_SE = cache.cells[fs_cache::calculateNbrID(1-1,1  ,1+1)]->derivatives;
   
   calculateWaveSpeedXY(cp_SE, derivs_SE, nbr_cp_SE, nbr_derivs_SE, Bx_S, By_E, dBxdy_S, dBxdz_S, dBydx_E, dBydz_E, PLUS, MINUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_x = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_y = c_x;
   ax_neg = max(ax_neg,-Vx0 + c_x);
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);
   maxV = max(maxV, calculateCflSpeed(Vx0, Vy0, vA, vS, vW));

   // Ez and characteristic speeds on NW (j-1) cell:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vx0  = divideIfNonZero(cp_NW[CellParams::RHOMVX], cp_NW[CellParams::RHOM]);
      Vy0  = divideIfNonZero(cp_NW[CellParams::RHOMVY], cp_NW[CellParams::RHOM]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vx0  = divideIfNonZero(cp_NW[CellParams::RHOMVX_DT2], cp_NW[CellParams::RHOM_DT2]);
      Vy0  = divideIfNonZero(cp_NW[CellParams::RHOMVY_DT2], cp_NW[CellParams::RHOM_DT2]);
   }

   // 1st order terms:
   Real Ez_NW = Bx_N*Vy0 - By_W*Vx0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ez_NW += Parameters::resistivity *
     sqrt((cp_NW[CellParams::BGBX]+cp_NW[CellParams::PERBX])*
          (cp_NW[CellParams::BGBX]+cp_NW[CellParams::PERBX]) +
          (cp_NW[CellParams::BGBY]+cp_NW[CellParams::PERBY])*
          (cp_NW[CellParams::BGBY]+cp_NW[CellParams::PERBY]) +
          (cp_NW[CellParams::BGBZ]+cp_NW[CellParams::PERBZ])*
          (cp_NW[CellParams::BGBZ]+cp_NW[CellParams::PERBZ])
         ) /
     cp_NW[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_NW[fs::dPERBydx]/cp_NW[CellParams::DX] - derivs_NW[fs::dPERBxdy]/cp_NW[CellParams::DY]);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ez_NW += cp_NW[CellParams::EZHALL_010_011] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ez_NW += cp_NW[CellParams::EZGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_NW  += +HALF*((Bx_N + HALF*dBxdy_N)*(-derivs_NW[fs::dVydx] + derivs_NW[fs::dVydy]) + dBxdy_N*Vy0 + SIXTH*dBxdz_N*derivs_NW[fs::dVydz]);
      Ez_NW  += -HALF*((By_W - HALF*dBydx_W)*(-derivs_NW[fs::dVxdx] + derivs_NW[fs::dVxdy]) - dBydx_W*Vx0 + SIXTH*dBydz_W*derivs_NW[fs::dVxdz]);
   #endif
   
   creal* const nbr_cp_NW     = cache.cells[fs_cache::calculateNbrID(1  ,1-1,1+1)]->parameters;
   creal* const nbr_derivs_NW = cache.cells[fs_cache::calculateNbrID(1  ,1-1,1+1)]->derivatives;
   
   calculateWaveSpeedXY(cp_NW, derivs_NW, nbr_cp_NW, nbr_derivs_NW, Bx_N, By_W, dBxdy_N, dBxdz_N, dBydx_W, dBydz_W, MINUS, PLUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_x = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_y = c_x;
   ax_neg = max(ax_neg,-Vx0 + c_x); 
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);
   maxV = max(maxV, calculateCflSpeed(Vx0, Vy0, vA, vS, vW));
   
   // Ez and characteristic speeds on NE (i-1,j-1) cell:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vx0  = divideIfNonZero(cp_NE[CellParams::RHOMVX], cp_NE[CellParams::RHOM]);
      Vy0  = divideIfNonZero(cp_NE[CellParams::RHOMVY], cp_NE[CellParams::RHOM]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vx0  = divideIfNonZero(cp_NE[CellParams::RHOMVX_DT2], cp_NE[CellParams::RHOM_DT2]);
      Vy0  = divideIfNonZero(cp_NE[CellParams::RHOMVY_DT2], cp_NE[CellParams::RHOM_DT2]);
   }
   
   // 1st order terms:
   Real Ez_NE = Bx_N*Vy0 - By_E*Vx0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ez_NE += Parameters::resistivity *
     sqrt((cp_NE[CellParams::BGBX]+cp_NE[CellParams::PERBX])*
          (cp_NE[CellParams::BGBX]+cp_NE[CellParams::PERBX]) +
          (cp_NE[CellParams::BGBY]+cp_NE[CellParams::PERBY])*
          (cp_NE[CellParams::BGBY]+cp_NE[CellParams::PERBY]) +
          (cp_NE[CellParams::BGBZ]+cp_NE[CellParams::PERBZ])*
          (cp_NE[CellParams::BGBZ]+cp_NE[CellParams::PERBZ])
         ) /
     cp_NE[CellParams::RHOQ] /
     physicalconstants::MU_0 *
     (derivs_NE[fs::dPERBydx]/cp_NE[CellParams::DX] - derivs_NE[fs::dPERBxdy]/cp_NE[CellParams::DY]);
   
   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ez_NE += cp_NE[CellParams::EZHALL_110_111] / (rhoq_S*physicalconstants::MU_0);
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ez_NE += cp_NE[CellParams::EZGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_NE  += +HALF*((Bx_N + HALF*dBxdy_N)*(+derivs_NE[fs::dVydx] + derivs_NE[fs::dVydy]) + dBxdy_N*Vy0 + SIXTH*dBxdz_N*derivs_NE[fs::dVydz]);
      Ez_NE  += -HALF*((By_E + HALF*dBydx_E)*(+derivs_NE[fs::dVxdx] + derivs_NE[fs::dVxdy]) + dBydx_E*Vx0 + SIXTH*dBydz_E*derivs_NE[fs::dVxdz]);
   #endif
   
   creal* const nbr_cp_NE     = cache.cells[fs_cache::calculateNbrID(1-1,1-1,1+1)]->parameters;
   creal* const nbr_derivs_NE = cache.cells[fs_cache::calculateNbrID(1-1,1-1,1+1)]->derivatives;
   
   calculateWaveSpeedXY(cp_NE, derivs_NE, nbr_cp_NE, nbr_derivs_NE, Bx_N, By_E, dBxdy_N, dBxdz_N, dBydx_E, dBydz_E, PLUS, PLUS, minRhom, maxRhom, RKCase, vA, vS, vW);
   c_x = min(Parameters::maxWaveVelocity,sqrt(vA*vA + vS*vS));
   c_y = c_x;
   ax_neg = max(ax_neg,-Vx0 + c_x);
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);
   maxV = max(maxV, calculateCflSpeed(Vx0, Vy0, vA, vS, vW));

   // Calculate properly upwinded edge-averaged Ez:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      cp_SW[CellParams::EZ] = ax_pos*ay_pos*Ez_NE + ax_pos*ay_neg*Ez_SE + ax_neg*ay_pos*Ez_NW + ax_neg*ay_neg*Ez_SW;
      cp_SW[CellParams::EZ] /= ((ax_pos+ax_neg)*(ay_pos+ay_neg)+EPS);

      if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
         cp_SW[CellParams::EZ] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(perBx_S-perBx_N);
         cp_SW[CellParams::EZ] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(perBy_W-perBy_E);
#else
         cp_SW[CellParams::EZ] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((perBx_S-HALF*dperBxdy_S) - (perBx_N+HALF*dperBxdy_N));
         cp_SW[CellParams::EZ] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((perBy_W-HALF*dperBydx_W) - (perBy_E+HALF*dperBydx_E));
#endif
      }
   } else { // RKCase == RK_ORDER2_STEP1
      cp_SW[CellParams::EZ_DT2] = ax_pos*ay_pos*Ez_NE + ax_pos*ay_neg*Ez_SE + ax_neg*ay_pos*Ez_NW + ax_neg*ay_neg*Ez_SW;
      cp_SW[CellParams::EZ_DT2] /= ((ax_pos+ax_neg)*(ay_pos+ay_neg)+EPS);

      if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
         cp_SW[CellParams::EZ_DT2] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(perBx_S-perBx_N);
         cp_SW[CellParams::EZ_DT2] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(perBy_W-perBy_E);
#else
         cp_SW[CellParams::EZ_DT2] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((perBx_S-HALF*dperBxdy_S) - (perBx_N+HALF*dperBxdy_N));
         cp_SW[CellParams::EZ_DT2] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((perBy_W-HALF*dperBydx_W) - (perBy_E+HALF*dperBydx_E));
#endif
      }
   }
   
   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      //compute maximum timestep for fieldsolver in this cell (CFL=1)      
      Real min_dx=std::numeric_limits<Real>::max();;
      min_dx=min(min_dx,cp_SW[CellParams::DX]);
      min_dx=min(min_dx,cp_SW[CellParams::DY]);
      //update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if(maxV!=ZERO) cp_SW[CellParams::MAXFDT]=min(cp_SW[CellParams::MAXFDT],min_dx/maxV);
   }
}

/*! \brief Electric field propagation function.
 * 
 * Calls the general or the system boundary electric field propagation functions.
 * 
 * \param mpiGrid Grid
 * \param cellCache Field solver cell cache
 * \param cells Vector of cells to process
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa calculateUpwindedElectricFieldSimple calculateEdgeElectricFieldX calculateEdgeElectricFieldY calculateEdgeElectricFieldZ
 * 
 */
void calculateElectricField(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   std::vector<fs_cache::CellCache>& cellCache,
   const std::vector<uint16_t>& cells,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const uint16_t localID = cells[c];
      fs_cache::CellCache& cache = cellCache[localID];
      const CellID cellID = cache.cellID;

      if (cache.sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;

      cuint fieldSolverSysBoundaryFlag = cache.existingCellsFlags;
      cuint cellSysBoundaryFlag        = cache.sysBoundaryFlag;
      cuint cellSysBoundaryLayer       = cache.cells[fs_cache::calculateNbrID(1,1,1)]->sysBoundaryLayer;

      if ((fieldSolverSysBoundaryFlag & CALCULATE_EX) == CALCULATE_EX) {
         if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) &&
             (cellSysBoundaryLayer != 1)
            ) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->
              fieldSolverBoundaryCondElectricField(mpiGrid, cellID, RKCase, 0);
         } else {
            calculateEdgeElectricFieldX(cache,RKCase);
         }
      }

      if ((fieldSolverSysBoundaryFlag & CALCULATE_EY) == CALCULATE_EY) {
         if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) &&
            (cellSysBoundaryLayer != 1)
           ) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->
              fieldSolverBoundaryCondElectricField(mpiGrid, cellID, RKCase, 1);
         } else {
            calculateEdgeElectricFieldY(cache,RKCase);
         }
      }

      if ((fieldSolverSysBoundaryFlag & CALCULATE_EZ) == CALCULATE_EZ) {
         if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) &&
            (cellSysBoundaryLayer != 1)
           ) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->
              fieldSolverBoundaryCondElectricField(mpiGrid, cellID, RKCase, 2);
         } else {
            calculateEdgeElectricFieldZ(cache,RKCase);
         }
      }
   } // for-loop over spatial cells
}

/*! \brief High-level electric field computation function.
 * 
 * Transfers the derivatives, calculates the edge electric fields and transfers the new electric fields.
 * 
 * \param mpiGrid Grid
 * \param sysBoundaries System boundary conditions existing
 * \param localCells Vector of local cells to process
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa calculateElectricField calculateEdgeElectricFieldX calculateEdgeElectricFieldY calculateEdgeElectricFieldZ
 */
void calculateUpwindedElectricFieldSimple(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase
) {
   namespace fs = fieldsolver;
   int timer;
   phiprof::start("Calculate upwinded electric field");
   uint64_t transferMask = 0;
   if(P::ohmHallTerm > 0) {
      transferMask = transferMask | Transfer::CELL_HALL_TERM;
   }
   if(P::ohmGradPeTerm > 0) {
      transferMask = transferMask | Transfer::CELL_GRADPE_TERM;
   }
   if(P::ohmHallTerm == 0 && P::ohmGradPeTerm == 0) {
      transferMask = Transfer::CELL_DERIVATIVES;
   }
   SpatialCell::set_mpi_transfer_type(transferMask);
   
   timer=phiprof::initializeTimer("Start communication in calculateUpwindedElectricFieldSimple","MPI");
   phiprof::start(timer);
   mpiGrid.start_remote_neighbor_copy_updates(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);
   
   // Calculate upwinded electric field on inner cells
   timer=phiprof::initializeTimer("Compute inner cells");
   phiprof::start(timer);
   calculateElectricField(mpiGrid,fs_cache::getCache().localCellsCache,
                          fs_cache::getCache().cellsWithLocalNeighbours,
                          sysBoundaries,RKCase);
   phiprof::stop(timer,fs_cache::getCache().cellsWithLocalNeighbours.size(),"Spatial Cells");
   
   timer=phiprof::initializeTimer("Wait for receives","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_receives(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);

   // Calculate upwinded electric field on boundary cells:
   timer=phiprof::initializeTimer("Compute boundary cells");
   phiprof::start(timer);
   calculateElectricField(mpiGrid,fs_cache::getCache().localCellsCache,   
                          fs_cache::getCache().cellsWithRemoteNeighbours,
                          sysBoundaries,RKCase);
   phiprof::stop(timer,fs_cache::getCache().cellsWithRemoteNeighbours.size(),"Spatial Cells");


   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_sends();
   phiprof::stop(timer);
   
   // Exchange electric field with neighbouring processes
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_E);
   } else { // RKCase == RK_ORDER2_STEP1
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_EDT2);
   }
   timer=phiprof::initializeTimer("Communicate electric fields","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.update_copies_of_remote_neighbors(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);

   const size_t N_cells = fs_cache::getCache().cellsWithLocalNeighbours.size() 
     + fs_cache::getCache().cellsWithRemoteNeighbours.size();
   
   phiprof::stop("Calculate upwinded electric field",N_cells,"Spatial Cells");
}
