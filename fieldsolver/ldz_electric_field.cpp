/*
This file is part of Vlasiator.

Copyright 2015 Finnish Meteorological Institute
*/

#include <cstdlib>

#include "fs_common.h"
#include "ldz_electric_field.hpp"

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

namespace fs = fieldsolver;
namespace pc = physicalconstants;
using namespace std;

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
 * \param minRho Minimum density allowed from the neighborhood
 * \param maxRho Maximum density allowed from the neighborhood
 */
Real calculateWaveSpeedYZ(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   cint i,
   cint j,
   cint k,
   cint nbi,
   cint nbj,
   cint nbk,
   const Real& By,
   const Real& Bz,
   const Real& dBydx,
   const Real& dBydz,
   const Real& dBzdx,
   const Real& dBzdy,
   const Real& ydir,
   const Real& zdir,
   const Real& minRho,
   const Real& maxRho
) {
   if (Parameters::propagateField == false) return 0.0;
   
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb = perBGrid.get(i,j,k);
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const nbr_perb = perBGrid.get(nbi,nbj,nbk);
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments = momentsGrid.get(i,j,k);
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments = dMomentsGrid.get(i,j,k);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb = dPerBGrid.get(i,j,k);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const nbr_dperb = dPerBGrid.get(nbi,nbj,nbk);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb = BgBGrid.get(i,j,k);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const  nbr_bgb = BgBGrid.get(nbi,nbj,nbk);
   
   Real A_0, A_X, rhom, p11, p22, p33;
   A_0  = HALF*(nbr_perb[fsgrids::bfield::PERBX] + nbr_bgb[fsgrids::bgbfield::BGBX] + perb[fsgrids::bfield::PERBX] + bgb[fsgrids::bgbfield::BGBX]);
   A_X  = (nbr_perb[fsgrids::bfield::PERBX] + nbr_bgb[fsgrids::bgbfield::BGBX]) - (perb[fsgrids::bfield::PERBX] + bgb[fsgrids::bgbfield::BGBX]);
   rhom = moments[fsgrids::moments::RHO] + ydir*HALF*dmoments[fsgrids::dmoments::drhody] + zdir*HALF*dmoments[fsgrids::dmoments::drhodz];
   p11 = moments[fsgrids::moments::P_11] + ydir*HALF*dmoments[fsgrids::dmoments::dp11dy] + zdir*HALF*dmoments[fsgrids::dmoments::dp11dz];
   p22 = moments[fsgrids::moments::P_22] + ydir*HALF*dmoments[fsgrids::dmoments::dp22dy] + zdir*HALF*dmoments[fsgrids::dmoments::dp22dz];
   p33 = moments[fsgrids::moments::P_33] + ydir*HALF*dmoments[fsgrids::dmoments::dp33dy] + zdir*HALF*dmoments[fsgrids::dmoments::dp33dz];
   
   if (rhom < minRho) {
      rhom = minRho;
   } else if (rhom > maxRho) {
      rhom = maxRho;
   }
   rhom *= pc::MASS_PROTON;
   
   const Real A_Y  = nbr_dperb[fsgrids::dperb::dPERBxdy] + nbr_dbgb[fsgrids::bgbfield::dBGBxdy] + dperb[fsgrids::dperb::dPERBxdy] + dbgb[fsgrids::dbgb::dBGBxdy];
   const Real A_XY = nbr_dperb[fsgrids::dperb::dPERBxdy] + nbr_dbgb[fsgrids::bgbfield::dBGBxdy] - (dperb[fsgrids::dperb::dPERBxdy] + dbgb[fsgrids::dbgb::dBGBxdy]);
   const Real A_Z  = nbr_dperb[fsgrids::dperb::dPERBxdz] + nbr_dbgb[fsgrids::bgbfield::dBGBxdz] + dperb[fsgrids::dperb::dPERBxdz] + dbgb[fsgrids::dbgb::dBGBxdz];
   const Real A_XZ = nbr_dperb[fsgrids::dperb::dPERBxdz] + nbr_dbgb[fsgrids::bgbfield::dBGBxdz] - (dperb[fsgrids::dperb::dPERBxdz] + dbgb[fsgrids::dbgb::dBGBxdz]);

   const Real Bx2  = (A_0 + ydir*HALF*A_Y + zdir*HALF*A_Z)*(A_0 + ydir*HALF*A_Y + zdir*HALF*A_Z)
     + TWELWTH*(A_X + ydir*HALF*A_XY + zdir*HALF*A_XZ)*(A_X + ydir*HALF*A_XY + zdir*HALF*A_XZ); // OK
   const Real By2  = (By + zdir*HALF*dBydz)*(By + zdir*HALF*dBydz) + TWELWTH*dBydx*dBydx; // OK
   const Real Bz2  = (Bz + ydir*HALF*dBzdy)*(Bz + ydir*HALF*dBzdy) + TWELWTH*dBzdx*dBzdx; // OK

   p11 = p11 < 0.0 ? 0.0 : p11;
   p22 = p22 < 0.0 ? 0.0 : p22;
   p33 = p33 < 0.0 ? 0.0 : p33;

   const Real vA2 = divideIfNonZero(Bx2+By2+Bz2, pc::MU_0*rhom); // Alfven speed
   const Real vS2 = divideIfNonZero(p11+p22+p33, 2.0*rhom); // sound speed, adiabatic coefficient 3/2, P=1/3*trace in sound speed
   const Real vW = Parameters::ohmHallTerm > 0 ? divideIfNonZero(2.0*M_PI*vA2*pc::MASS_PROTON, 
                                                                 perBGrid.DX*pc::CHARGE*sqrt(Bx2+By2+Bz2)) : 0.0; // whistler speed

   return min(Parameters::maxWaveVelocity,sqrt(vA2 + vS2) + vW);
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
 * \param minRho Minimum density allowed from the neighborhood
 * \param maxRho Maximum density allowed from the neighborhood
 */
Real calculateWaveSpeedXZ(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   cint i,
   cint j,
   cint k,
   cint nbi,
   cint nbj,
   cint nbk,
   const Real& Bx,
   const Real& Bz,
   const Real& dBxdy,
   const Real& dBxdz,
   const Real& dBzdx,
   const Real& dBzdy,
   const Real& xdir,
   const Real& zdir,
   const Real& minRho,
   const Real& maxRho
) {
   if (Parameters::propagateField == false) return 0.0;
   
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb = perBGrid.get(i,j,k);
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const nbr_perb = perBGrid.get(nbi,nbj,nbk);
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments = momentsGrid.get(i,j,k);
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments = dMomentsGrid.get(i,j,k);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb = dPerBGrid.get(i,j,k);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const nbr_dperb = dPerBGrid.get(nbi,nbj,nbk);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb = BgBGrid.get(i,j,k);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const  nbr_bgb = BgBGrid.get(nbi,nbj,nbk);
   
   Real B_0, B_Y, rhom, p11, p22, p33;
   B_0  = HALF*(nbr_perb[fsgrids::bfield::PERBY] + nbr_bgb[fsgrids::bgbfield::BGBY] + perb[fsgrids::bfield::PERBY] + bgb[fsgrids::bgbfield::BGBY]);
   B_Y  = (nbr_perb[fsgrids::bfield::PERBY] + nbr_bgb[fsgrids::bgbfield::BGBY]) - (perb[fsgrids::bfield::PERBY] + bgb[fsgrids::bgbfield::BGBY]);
   rhom = moments[fsgrids::moments::RHO] + xdir*HALF*dmoments[fsgrids::dmoments::drhodx] + zdir*HALF*dmoments[fsgrids::dmoments::drhodz];
   p11 = moments[fsgrids::moments::P_11] + xdir*HALF*dmoments[fsgrids::dmoments::dp11dx] + zdir*HALF*dmoments[fsgrids::dmoments::dp11dz];
   p22 = moments[fsgrids::moments::P_22] + xdir*HALF*dmoments[fsgrids::dmoments::dp22dx] + zdir*HALF*dmoments[fsgrids::dmoments::dp22dz];
   p33 = moments[fsgrids::moments::P_33] + xdir*HALF*dmoments[fsgrids::dmoments::dp33dx] + zdir*HALF*dmoments[fsgrids::dmoments::dp33dz];
   
   if (rhom < minRho) {
      rhom = minRho;
   } else if (rhom > maxRho) {
      rhom = maxRho;
   }
   rhom *= pc::MASS_PROTON;
   
   const Real B_X  = nbr_dperb[fsgrids::dperb::dPERBydx] + nbr_dbgb[fsgrids::bgbfield::dBGBydx] + dperb[fsgrids::dperb::dPERBydx] + dbgb[fsgrids::dbgb::dBGBydx];
   const Real B_XY = nbr_dperb[fsgrids::dperb::dPERBydx] + nbr_dbgb[fsgrids::bgbfield::dBGBydx] - (dperb[fsgrids::dperb::dPERBydx] + dbgb[fsgrids::dbgb::dBGBydx]);
   const Real B_Z  = nbr_dperb[fsgrids::dperb::dPERBydz] + nbr_dbgb[fsgrids::bgbfield::dBGBydz] + dperb[fsgrids::dperb::dPERBydz] + dbgb[fsgrids::dbgb::dBGBydz];
   const Real B_YZ = nbr_dperb[fsgrids::dperb::dPERBydz] + nbr_dbgb[fsgrids::bgbfield::dBGBydz] - (dperb[fsgrids::dperb::dPERBydz] + dbgb[fsgrids::dbgb::dBGBydz]);
      
   const Real By2  = (B_0 + xdir*HALF*B_X + zdir*HALF*B_Z)*(B_0 + xdir*HALF*B_X + zdir*HALF*B_Z)
     + TWELWTH*(B_Y + xdir*HALF*B_XY + zdir*HALF*B_YZ)*(B_Y + xdir*HALF*B_XY + zdir*HALF*B_YZ); // OK
   const Real Bx2  = (Bx + zdir*HALF*dBxdz)*(Bx + zdir*HALF*dBxdz) + TWELWTH*dBxdy*dBxdy; // OK
   const Real Bz2  = (Bz + xdir*HALF*dBzdx)*(Bz + xdir*HALF*dBzdx) + TWELWTH*dBzdy*dBzdy; // OK
   
   p11 = p11 < 0.0 ? 0.0 : p11;
   p22 = p22 < 0.0 ? 0.0 : p22;
   p33 = p33 < 0.0 ? 0.0 : p33;
   
   const Real vA2 = divideIfNonZero(Bx2+By2+Bz2, pc::MU_0*rhom); // Alfven speed
   const Real vS2 = divideIfNonZero(p11+p22+p33, 2.0*rhom); // sound speed, adiabatic coefficient 3/2, P=1/3*trace in sound speed
   const Real vW = Parameters::ohmHallTerm > 0 ? divideIfNonZero(2.0*M_PI*vA2*pc::MASS_PROTON, perBGrid.DX*pc::CHARGE*sqrt(Bx2+By2+Bz2)) : 0.0; // whistler speed
      
   return min(Parameters::maxWaveVelocity,sqrt(vA2 + vS2) + vW);
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
 * \param minRho Minimum density allowed from the neighborhood
 * \param maxRho Maximum density allowed from the neighborhood
 */
Real calculateWaveSpeedXY(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   cint i,
   cint j,
   cint k,
   cint nbi,
   cint nbj,
   cint nbk,
   const Real& Bx,
   const Real& By,
   const Real& dBxdy,
   const Real& dBxdz,
   const Real& dBydx,
   const Real& dBydz,
   const Real& xdir,
   const Real& ydir,
   const Real& minRho,
   const Real& maxRho
) {
   if (Parameters::propagateField == false) return 0.0;
   
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb = perBGrid.get(i,j,k);
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const nbr_perb = perBGrid.get(nbi,nbj,nbk);
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments = momentsGrid.get(i,j,k);
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments = dMomentsGrid.get(i,j,k);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb = dPerBGrid.get(i,j,k);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const nbr_dperb = dPerBGrid.get(nbi,nbj,nbk);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb = BgBGrid.get(i,j,k);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const  nbr_bgb = BgBGrid.get(nbi,nbj,nbk);
   
   Real C_0, C_Z, rhom, p11, p22, p33;
   C_0  = HALF*(nbr_perb[fsgrids::bfield::PERBZ] + nbr_bgb[fsgrids::bgbfield::BGBZ] + perb[fsgrids::bfield::PERBZ] + bgb[fsgrids::bgbfield::BGBZ]);
   C_Z  = (nbr_perb[fsgrids::bfield::PERBZ] + nbr_bgb[fsgrids::bgbfield::BGBZ]) - (perb[fsgrids::bfield::PERBZ] + bgb[fsgrids::bgbfield::BGBZ]);
   rhom = moments[fsgrids::moments::RHO] + xdir*HALF*dmoments[fsgrids::dmoments::drhodx] + ydir*HALF*dmoments[fsgrids::dmoments::drhody];
   p11 = moments[fsgrids::moments::P_11] + xdir*HALF*dmoments[fsgrids::dmoments::dp11dx] + ydir*HALF*dmoments[fsgrids::dmoments::dp11dy];
   p22 = moments[fsgrids::moments::P_22] + xdir*HALF*dmoments[fsgrids::dmoments::dp22dx] + ydir*HALF*dmoments[fsgrids::dmoments::dp22dy];
   p33 = moments[fsgrids::moments::P_33] + xdir*HALF*dmoments[fsgrids::dmoments::dp33dx] + ydir*HALF*dmoments[fsgrids::dmoments::dp33dy];
   
   if (rhom < minRho) {
      rhom = minRho;
   } else if (rhom > maxRho) {
      rhom = maxRho;
   }
   rhom *= pc::MASS_PROTON;
   
   const Real C_X  = nbr_dperb[fsgrids::dperb::dPERBzdx] + nbr_dbgb[fsgrids::bgbfield::dBGBzdx] + dperb[fsgrids::dperb::dPERBzdx] + dbgb[fsgrids::dbgb::dBGBzdx];
   const Real C_XZ = nbr_dperb[fsgrids::dperb::dPERBzdx] + nbr_dbgb[fsgrids::bgbfield::dBGBzdx] - (dperb[fsgrids::dperb::dPERBzdx] + dbgb[fsgrids::dbgb::dBGBzdx]);
   const Real C_Y  = nbr_dperb[fsgrids::dperb::dPERBzdy] + nbr_dbgb[fsgrids::bgbfield::dBGBzdy] + dperb[fsgrids::dperb::dPERBzdy] + dbgb[fsgrids::dbgb::dBGBzdy];
   const Real C_YZ = nbr_dperb[fsgrids::dperb::dPERBzdy] + nbr_dbgb[fsgrids::bgbfield::dBGBzdy] - (dperb[fsgrids::dperb::dPERBzdy] + dbgb[fsgrids::dbgb::dBGBzdy]);
   
   const Real Bz2  = (C_0 + xdir*HALF*C_X + ydir*HALF*C_Y)*(C_0 + xdir*HALF*C_X + ydir*HALF*C_Y)
     + TWELWTH*(C_Z + xdir*HALF*C_XZ + ydir*HALF*C_YZ)*(C_Z + xdir*HALF*C_XZ + ydir*HALF*C_YZ);
   const Real Bx2  = (Bx + ydir*HALF*dBxdy)*(Bx + ydir*HALF*dBxdy) + TWELWTH*dBxdz*dBxdz;
   const Real By2  = (By + xdir*HALF*dBydx)*(By + xdir*HALF*dBydx) + TWELWTH*dBydz*dBydz;
   
   p11 = p11 < 0.0 ? 0.0 : p11;
   p22 = p22 < 0.0 ? 0.0 : p22;
   p33 = p33 < 0.0 ? 0.0 : p33;
      
   const Real vA2 = divideIfNonZero(Bx2+By2+Bz2, pc::MU_0*rhom); // Alfven speed
   const Real vS2 = divideIfNonZero(p11+p22+p33, 2.0*rhom); // sound speed, adiabatic coefficient 3/2, P=1/3*trace in sound speed
   const Real vW = Parameters::ohmHallTerm > 0 ? divideIfNonZero(2.0*M_PI*vA2*pc::MASS_PROTON, perBGrid.DX*pc::CHARGE*sqrt(Bx2+By2+Bz2)) : 0.0; // whistler speed

   return min(Parameters::maxWaveVelocity,sqrt(vA2 + vS2) + vW);
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
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   i,
   j,
   k,
   cint& RKCase
) {
   #ifdef DEBUG_FSOLVER
   bool ok = true;
   if (technicalGrid.get(i,j,k) == NULL) ok = false;
   if (technicalGrid.get(i,j-1,k) == NULL) ok = false;
   if (technicalGrid.get(i,j-1,k-1) == NULL) ok = false;
   if (technicalGrid.get(i,j,k-1) == NULL) ok = false;
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
   Real c_y, c_z;                   // Wave speeds to yz-directions
   
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_SW = perBGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_SE = perBGrid.get(i  ,j-1,k  );
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_NE = perBGrid.get(i  ,j-1,k-1);
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_NW = perBGrid.get(i  ,j  ,k-1);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_SW = BgBGrid.get(i,j  ,k  );
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_SE = BgBGrid.get(i,j-1,k  );
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_NE = BgBGrid.get(i,j-1,k-1);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_NW = BgBGrid.get(i,j  ,k-1);
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_SW = momentsGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_SE = momentsGrid.get(i  ,j-1,k  );
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_NE = momentsGrid.get(i  ,j-1,k-1);
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_NW = momentsGrid.get(i  ,j  ,k-1);
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_SW = dMomentsGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_SE = dMomentsGrid.get(i  ,j-1,k  );
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_NE = dMomentsGrid.get(i  ,j-1,k-1);
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_NW = dMomentsGrid.get(i  ,j  ,k-1);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_SW = dPerBGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_SE = dPerBGrid.get(i  ,j-1,k  );
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_NE = dPerBGrid.get(i  ,j-1,k-1);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_NW = dPerBGrid.get(i  ,j  ,k-1);
   
   const std::array<Real, fsgrids::efield::N_EFIELD> * efield_SW = EGrid.get(i,j,k);
   
   Real By_S, Bz_W, Bz_E, By_N, perBy_S, perBz_W, perBz_E, perBy_N;
   Real minRho = std::numeric_limits<Real>::max();
   Real maxRho = std::numeric_limits<Real>::min();
   By_S = perb_SW[fsgrids::bfield::PERBY]+bgb_SW[fsgrids::bgbfield::BGBY];
   Bz_W = perb_SW[fsgrids::bfield::PERBZ]+bgb_SW[fsgrids::bgbfield::BGBZ];
   Bz_E = perb_SE[fsgrids::bfield::PERBZ]+bgb_SE[fsgrids::bgbfield::BGBZ];
   By_N = perb_NW[fsgrids::bfield::PERBY]+bgb_NW[fsgrids::bgbfield::BGBY];
   perBy_S = perb_SW[fsgrids::bfield::PERBY];
   perBz_W = perb_SW[fsgrids::bfield::PERBZ];
   perBz_E = perb_SE[fsgrids::bfield::PERBZ];
   perBy_N = perb_NW[fsgrids::bfield::PERBY];
   Vy0  = divideIfNonZero(moments_SW[fsgrids::moments::RHOVY], moments_SW[fsgrids::moments::RHO]);
   Vz0  = divideIfNonZero(moments_SW[fsgrids::moments::RHOVZ], moments_SW[fsgrids::moments::RHO]);
   minRho = min(minRho,
               min(moments_SW[fsgrids::moments::RHO],
                  min(moments_SE[fsgrids::moments::RHO],
                     min(moments_NW[fsgrids::moments::RHO],
                           moments_NE[fsgrids::moments::RHO])
                     )
                  )
               );
   maxRho = max(maxRho,
               max(moments_SW[fsgrids::moments::RHO],
                  max(moments_SE[fsgrids::moments::RHO],
                     max(moments_NW[fsgrids::moments::RHO],
                           moments_NE[fsgrids::moments::RHO])
                     )
                  )
               );

   creal dBydx_S = dperb_SW[fsgrids::dperb::dPERBydx] + dbgb_SW[fsgrids::bgbfield::dBGBydx];
   creal dBydz_S = dperb_SW[fsgrids::dperb::dPERBydz] + dbgb_SW[fsgrids::bgbfield::dBGBydz];
   creal dBzdx_W = dperb_SW[fsgrids::dperb::dPERBzdx] + dbgb_SW[fsgrids::bgbfield::dBGBzdx];
   creal dBzdy_W = dperb_SW[fsgrids::dperb::dPERBzdy] + dbgb_SW[fsgrids::bgbfield::dBGBzdy];
   creal dBzdx_E = dperb_SE[fsgrids::dperb::dPERBzdx] + dbgb_SE[fsgrids::bgbfield::dBGBzdx];
   creal dBzdy_E = dperb_SE[fsgrids::dperb::dPERBzdy] + dbgb_SE[fsgrids::bgbfield::dBGBzdy];
   creal dBydx_N = dperb_NW[fsgrids::dperb::dPERBydx] + dbgb_NW[fsgrids::bgbfield::dBGBydx];
   creal dBydz_N = dperb_NW[fsgrids::dperb::dPERBydz] + dbgb_NW[fsgrids::bgbfield::dBGBydz];
   creal dperBydz_S = dperb_SW[fsgrids::dperb::dPERBydz];
   creal dperBydz_N = dperb_NW[fsgrids::dperb::dPERBydz];
   creal dperBzdy_W = dperb_SW[fsgrids::dperb::dPERBzdy];
   creal dperBzdy_E = dperb_SE[fsgrids::dperb::dPERBzdy];

   // Ex and characteristic speeds on this cell:
   // 1st order terms:
   Real Ex_SW = By_S*Vz0 - Bz_W*Vy0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ex_SW += Parameters::resistivity *
     sqrt((bgb_SW[fsgrids::bgbfield::BGBX]+perb_SW[fsgrids::bfield::PERBX])*
          (bgb_SW[fsgrids::bgbfield::BGBX]+perb_SW[fsgrids::bfield::PERBX]) +
          (bgb_SW[fsgrids::bgbfield::BGBY]+perb_SW[fsgrids::bfield::PERBY])*
          (bgb_SW[fsgrids::bgbfield::BGBY]+perb_SW[fsgrids::bfield::PERBY]) +
          (bgb_SW[fsgrids::bgbfield::BGBZ]+perb_SW[fsgrids::bfield::PERBZ])*
          (bgb_SW[fsgrids::bgbfield::BGBZ]+perb_SW[fsgrids::bfield::PERBZ])
         ) /
     (moments_SW[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_SW[fsgrids::dperb::dPERBzdy]/technicalGrid.DY - dperb_SW[fsgrids::dperb::dPERBydz]/technicalGrid.DZ);
   
   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ex_SW += EHallGrid.get(i,j,k)[fsgrids::ehall::EXHALL_000_100];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ex_SW += EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EXGRADPE];
   }

   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_SW += +HALF*((By_S - HALF*dBydz_S)*(-dmoments_SW[fsgrids::dmoments::dVzdy] - dmoments_SW[fsgrids::dmoments::dVzdz]) - dBydz_S*Vz0 + SIXTH*dBydx_S*dmoments_SW[fsgrids::dmoments::dVzdx]);
      Ex_SW += -HALF*((Bz_W - HALF*dBzdy_W)*(-dmoments_SW[fsgrids::dmoments::dVydy] - dmoments_SW[fsgrids::dmoments::dVydz]) - dBzdy_W*Vy0 + SIXTH*dBzdx_W*dmoments_SW[fsgrids::dmoments::dVydx]);
   #endif
   
   c_y = calculateWaveSpeedYZ(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i  , j, k,
      i+1, j, k,
      By_S, Bz_W, dBydx_S, dBydz_S, dBzdx_W, dBzdy_W, MINUS, MINUS, minRho, maxRho);
   c_z = c_y;
   ay_neg   = max(ZERO,-Vy0 + c_y);
   ay_pos   = max(ZERO,+Vy0 + c_y);
   az_neg   = max(ZERO,-Vz0 + c_z);
   az_pos   = max(ZERO,+Vz0 + c_z);

   // Ex and characteristic speeds on j-1 neighbour:
   Vy0  = divideIfNonZero(moments_SE[fsgrids::moments::RHOVY], moments_SE[fsgrids::moments::RHO]);
   Vz0  = divideIfNonZero(moments_SE[fsgrids::moments::RHOVZ], moments_SE[fsgrids::moments::RHO]);

   // 1st order terms:
   Real Ex_SE = By_S*Vz0 - Bz_E*Vy0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ex_SE += Parameters::resistivity *
     sqrt((bgb_SE[fsgrids::bgbfield::BGBX]+perb_SE[fsgrids::bfield::PERBX])*
          (bgb_SE[fsgrids::bgbfield::BGBX]+perb_SE[fsgrids::bfield::PERBX]) +
          (bgb_SE[fsgrids::bgbfield::BGBY]+perb_SE[fsgrids::bfield::PERBY])*
          (bgb_SE[fsgrids::bgbfield::BGBY]+perb_SE[fsgrids::bfield::PERBY]) +
          (bgb_SE[fsgrids::bgbfield::BGBZ]+perb_SE[fsgrids::bfield::PERBZ])*
          (bgb_SE[fsgrids::bgbfield::BGBZ]+perb_SE[fsgrids::bfield::PERBZ])
         ) /
     (moments_SE[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_SE[fsgrids::dperb::dPERBzdy]/technicalGrid.DY - dperb_SE[fsgrids::dperb::dPERBydz]/technicalGrid.DZ);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ex_SE += EHallGrid.get(i,j-1,k)[fsgrids::ehall::EXHALL_010_110];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ex_SE += EGradPeGrid.get(i,j-1,k)[fsgrids::egradpe::EXGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_SE += +HALF*((By_S - HALF*dBydz_S)*(+dmoments_SE[fsgrids::dmoments::dVzdy] - dmoments_SE[fsgrids::dmoments::dVzdz]) - dBydz_S*Vz0 + SIXTH*dBydx_S*dmoments_SE[fsgrids::dmoments::dVzdx]);
      Ex_SE += -HALF*((Bz_E + HALF*dBzdy_E)*(+dmoments_SE[fsgrids::dmoments::dVydy] - dmoments_SE[fsgrids::dmoments::dVydz]) + dBzdy_E*Vy0 + SIXTH*dBzdx_E*dmoments_SE[fsgrids::dmoments::dVydx]);
   #endif

   c_y = calculateWaveSpeedYZ(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i  ,   j, k,
      i+1, j-1, k,
      By_S, Bz_E, dBydx_S, dBydz_S, dBzdx_E, dBzdy_E, PLUS, MINUS, minRho, maxRho);
   c_z = c_y;
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);

   // Ex and characteristic speeds on k-1 neighbour:
   Vy0  = divideIfNonZero(moments_NW[fsgrids::moments::RHOVY], moments_NW[fsgrids::moments::RHO]);
   Vz0  = divideIfNonZero(moments_NW[fsgrids::moments::RHOVZ], moments_NW[fsgrids::moments::RHO]);

   // 1st order terms:
   Real Ex_NW    = By_N*Vz0 - Bz_W*Vy0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ex_NW += Parameters::resistivity *
     sqrt((bgb_NW[fsgrids::bgbfield::BGBX]+perb_NW[fsgrids::bfield::PERBX])*
          (bgb_NW[fsgrids::bgbfield::BGBX]+perb_NW[fsgrids::bfield::PERBX]) +
          (bgb_NW[fsgrids::bgbfield::BGBY]+perb_NW[fsgrids::bfield::PERBY])*
          (bgb_NW[fsgrids::bgbfield::BGBY]+perb_NW[fsgrids::bfield::PERBY]) +
          (bgb_NW[fsgrids::bgbfield::BGBZ]+perb_NW[fsgrids::bfield::PERBZ])*
          (bgb_NW[fsgrids::bgbfield::BGBZ]+perb_NW[fsgrids::bfield::PERBZ])
         ) /
     (moments_NW[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_NW[fsgrids::dperb::dPERBzdy]/technicalGrid.DY - dperb_NW[fsgrids::dperb::dPERBydz]/technicalGrid.DZ);
   
   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ex_NW += EHallGrid.get(i,j,k-1)[fsgrids::ehall::EXHALL_001_101];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ex_NW += EGradPeGrid.get(i,j,k-1)[fsgrids::egradpe::EXGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_NW += +HALF*((By_N + HALF*dBydz_N)*(-dmoments_NW[fsgrids::dmoments::dVzdy] + dmoments_NW[fsgrids::dmoments::dVzdz]) + dBydz_N*Vz0 + SIXTH*dBydx_N*dmoments_NW[fsgrids::dmoments::dVzdx]);
      Ex_NW += -HALF*((Bz_W - HALF*dBzdy_W)*(-dmoments_NW[fsgrids::dmoments::dVydy] + dmoments_NW[fsgrids::dmoments::dVydz]) - dBzdy_W*Vy0 + SIXTH*dBzdx_W*dmoments_NW[fsgrids::dmoments::dVydx]);
   #endif
   
   c_y = calculateWaveSpeedYZ(
      perBGrid,
      perBDt2Grid,
      momentsGrid,
      momentsDt2Grid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i+1, j, k-1,
      By_N, Bz_W, dBydx_N, dBydz_N, dBzdx_W, dBzdy_W, MINUS, PLUS, minRho, maxRho);
   c_z = c_y;
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);

   // Ex and characteristic speeds on j-1,k-1 neighbour:
   Vy0 = divideIfNonZero(moments_NE[fsgrids::moments::RHOVY], moments_NE[fsgrids::moments::RHO]);
   Vz0 = divideIfNonZero(moments_NE[fsgrids::moments::RHOVZ], moments_NE[fsgrids::moments::RHO]);
   
   // 1st order terms:
   Real Ex_NE    = By_N*Vz0 - Bz_E*Vy0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ex_NE += Parameters::resistivity *
            sqrt((bgb_NE[fsgrids::bgbfield::BGBX]+perb_NE[fsgrids::bfield::PERBX])*
                 (bgb_NE[fsgrids::bgbfield::BGBX]+perb_NE[fsgrids::bfield::PERBX]) +
                 (bgb_NE[fsgrids::bgbfield::BGBY]+perb_NE[fsgrids::bfield::PERBY])*
                 (bgb_NE[fsgrids::bgbfield::BGBY]+perb_NE[fsgrids::bfield::PERBY]) +
                 (bgb_NE[fsgrids::bgbfield::BGBZ]+perb_NE[fsgrids::bfield::PERBZ])*
                 (bgb_NE[fsgrids::bgbfield::BGBZ]+perb_NE[fsgrids::bfield::PERBZ])
                ) /
               (moments_NE[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
            physicalconstants::MU_0 *
            (dperb_NE[fsgrids::dperb::dPERBzdy]/technicalGrid.DY - dperb_NE[fsgrids::dperb::dPERBydz]/technicalGrid.DZ);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ex_NE += EHallGrid.get(i,j-1,k-1)[fsgrids::ehall::EXHALL_011_111];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ex_NE += EGradPeGrid.get(i,j-1,k-1)[fsgrids::egradpe::EXGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_NE += +HALF*((By_N + HALF*dBydz_N)*(+dmoments_NE[fsgrids::dmoments::dVzdy] + dmoments_NE[fsgrids::dmoments::dVzdz]) + dBydz_N*Vz0 + SIXTH*dBydx_N*dmoments_NE[fsgrids::dmoments::dVzdx]);
      Ex_NE += -HALF*((Bz_E + HALF*dBzdy_E)*(+dmoments_NE[fsgrids::dmoments::dVydy] + dmoments_NE[fsgrids::dmoments::dVydz]) + dBzdy_E*Vy0 + SIXTH*dBzdx_E*dmoments_NE[fsgrids::dmoments::dVydx]);
   #endif
   
   c_y = calculateWaveSpeedYZ(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i+1,j-1,k-1,
      By_N, Bz_E, dBydx_N, dBydz_N, dBzdx_E, dBzdy_E, PLUS, PLUS, minRho, maxRho);
   c_z = c_y;
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   
   // Calculate properly upwinded edge-averaged Ex:
   efield_SW[fsgrids::efield::EX]  = ay_pos*az_pos*Ex_NE + ay_pos*az_neg*Ex_SE + ay_neg*az_pos*Ex_NW + ay_neg*az_neg*Ex_SW;
   efield_SW[fsgrids::efield::EX] /= ((ay_pos+ay_neg)*(az_pos+az_neg)+EPS);
   if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
      // 1st order diffusive terms:
      efield_SW[fsgrids::efield::EX] -= az_pos*az_neg/(az_pos+az_neg+EPS)*(perBy_S-perBy_N);
      efield_SW[fsgrids::efield::EX] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(perBz_W-perBz_E);
#else
      // 2nd     order diffusive terms
      efield_SW[fsgrids::efield::EX] -= az_pos*az_neg/(az_pos+az_neg+EPS)*((perBy_S-HALF*dperBydz_S) - (perBy_N+HALF*dperBydz_N));
      efield_SW[fsgrids::efield::EX] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((perBz_W-HALF*dperBzdy_W) - (perBz_E+HALF*dperBzdy_E));
#endif
   }
   
   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      //compute maximum timestep for fieldsolver in this cell (CFL=1)
      Real max_a=ZERO;
      max_a=max(fabs(az_neg),max_a); 
      max_a=max(fabs(az_pos),max_a);
      max_a=max(fabs(ay_neg),max_a);
      max_a=max(fabs(ay_pos),max_a);
      Real min_dx=std::numeric_limits<Real>::max();
      min_dx=min(min_dx,technicalGrid.DY);
      min_dx=min(min_dx,technicalGrid.DZ);
      //update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (max_a != ZERO) technicalGrid.get(i,j,k)->maxFsDt = min(technicalGrid.get(i,j,k)->maxFsDt,min_dx/max_a);
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
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   i,
   j,
   k,
   cint& RKCase
) {
   #ifdef DEBUG_FSOLVER
   bool ok = true;
   if (technicalGrid.get(i,j,k) == NULL) ok = false;
   if (technicalGrid.get(i,j,k-1) == NULL) ok = false;
   if (technicalGrid.get(i-1,j,k) == NULL) ok = false;
   if (technicalGrid.get(i-1,j,k-1) == NULL) ok = false;
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
   Real c_x,c_z;                    // Wave speeds to xz-directions
   
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_SW = perBGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_SE = perBGrid.get(i  ,j  ,k-1);
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_NW = perBGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_NE = perBGrid.get(i-1,j  ,k-1);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_SW = BgBGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_SE = BgBGrid.get(i  ,j  ,k-1);
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_NW = BgBGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_NE = BgBGrid.get(i-1,j  ,k-1);
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_SW = momentsGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_SE = momentsGrid.get(i  ,j  ,k-1);
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_NW = momentsGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_NE = momentsGrid.get(i-1,j  ,k-1);
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_SW = dMomentsGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_SE = dMomentsGrid.get(i  ,j  ,k-1);
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_NW = dMomentsGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_NE = dMomentsGrid.get(i-1,j  ,k-1);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_SW = dPerBGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_SE = dPerBGrid.get(i  ,j  ,k-1);
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_NW = dPerBGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_NE = dPerBGrid.get(i-1,j  ,k-1);
   
   const std::array<Real, fsgrids::efield::N_EFIELD> * efield_SW = EGrid.get(i,j,k);
   
   // Fetch required plasma parameters:
   Real Bz_S, Bx_W, Bx_E, Bz_N, perBz_S, perBx_W, perBx_E, perBz_N;
   Real minRho = std::numeric_limits<Real>::max();
   Real maxRho = std::numeric_limits<Real>::min();
   Bz_S = perb_SW[fsgrids::bfield::PERBZ]+bgb_SW[fsgrids::bgbfield::BGBZ];
   Bx_W = perb_SW[fsgrids::bfield::PERBX]+bgb_SW[fsgrids::bgbfield::BGBX];
   Bx_E = perb_SE[fsgrids::bfield::PERBX]+bgb_SE[fsgrids::bgbfield::BGBX];
   Bz_N = perb_NW[fsgrids::bfield::PERBZ]+bgb_NW[fsgrids::bgbfield::BGBZ];
   perBz_S = perb_SW[fsgrids::bfield::PERBZ];
   perBx_W = perb_SW[fsgrids::bfield::PERBX];
   perBx_E = perb_SE[fsgrids::bfield::PERBX];
   perBz_N = perb_NW[fsgrids::bfield::PERBZ];
   Vx0  = divideIfNonZero(moments_SW[fsgrids::moments::RHOVX], moments_SW[fsgrids::moments::RHO]);
   Vz0  = divideIfNonZero(moments_SW[fsgrids::moments::RHOVZ], moments_SW[fsgrids::moments::RHO]);
   minRho = min(minRho,
               min(moments_SW[fsgrids::moments::RHO],
                  min(moments_SE[fsgrids::moments::RHO],
                     min(moments_NW[fsgrids::moments::RHO],
                           moments_NE[fsgrids::moments::RHO])
                     )
                  )
               );
   maxRho = max(maxRho,
               max(moments_SW[fsgrids::moments::RHO],
                  max(moments_SE[fsgrids::moments::RHO],
                     max(moments_NW[fsgrids::moments::RHO],
                           moments_NE[fsgrids::moments::RHO])
                     )
                  )
               );
   
   creal dBxdy_W = dperb_SW[fsgrids::dperb::dPERBxdy] + dbgb_SW[fsgrids::bgbfield::dBGBxdy];
   creal dBxdz_W = dperb_SW[fsgrids::dperb::dPERBxdz] + dbgb_SW[fsgrids::bgbfield::dBGBxdz];
   creal dBzdx_S = dperb_SW[fsgrids::dperb::dPERBzdx] + dbgb_SW[fsgrids::bgbfield::dBGBzdx];
   creal dBzdy_S = dperb_SW[fsgrids::dperb::dPERBzdy] + dbgb_SW[fsgrids::bgbfield::dBGBzdy];
   creal dBxdy_E = dperb_SE[fsgrids::dperb::dPERBxdy] + dbgb_SE[fsgrids::bgbfield::dBGBxdy];
   creal dBxdz_E = dperb_SE[fsgrids::dperb::dPERBxdz] + dbgb_SE[fsgrids::bgbfield::dBGBxdz];
   creal dBzdx_N = dperb_NW[fsgrids::dperb::dPERBzdx] + dbgb_NW[fsgrids::bgbfield::dBGBzdx];
   creal dBzdy_N = dperb_NW[fsgrids::dperb::dPERBzdy] + dbgb_NW[fsgrids::bgbfield::dBGBzdy];
   creal dperBzdx_S = dperb_SW[fsgrids::dperb::dPERBzdx];
   creal dperBzdx_N = dperb_NW[fsgrids::dperb::dPERBzdx];
   creal dperBxdz_W = dperb_SW[fsgrids::dperb::dPERBxdz];
   creal dperBxdz_E = dperb_SE[fsgrids::dperb::dPERBxdz];
   
   // Ey and characteristic speeds on this cell:
   // 1st order terms:
   Real Ey_SW  = Bz_S*Vx0 - Bx_W*Vz0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ey_SW += Parameters::resistivity *
     sqrt((bgb_SW[fsgrids::bgbfield::BGBX]+perb_SW[fsgrids::bfield::PERBX])*
          (bgb_SW[fsgrids::bgbfield::BGBX]+perb_SW[fsgrids::bfield::PERBX]) +
          (bgb_SW[fsgrids::bgbfield::BGBY]+perb_SW[fsgrids::bfield::PERBY])*
          (bgb_SW[fsgrids::bgbfield::BGBY]+perb_SW[fsgrids::bfield::PERBY]) +
          (bgb_SW[fsgrids::bgbfield::BGBZ]+perb_SW[fsgrids::bfield::PERBZ])*
          (bgb_SW[fsgrids::bgbfield::BGBZ]+perb_SW[fsgrids::bfield::PERBZ])
         ) /
     (moments_SW[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_SW[fsgrids::dperb::dPERBxdz]/technicalGrid.DZ - dperb_SW[fsgrids::dperb::dPERBzdx]/technicalGrid.DX);
   
   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ey_SW += EHallGrid.get(i,j,k)[fsgrids::ehall::EYHALL_000_010];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ey_SW += EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EYGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms
      Ey_SW += +HALF*((Bz_S - HALF*dBzdx_S)*(-dmoments_SW[fsgrids::dmoments::dVxdx] - dmoments_SW[fsgrids::dmoments::dVxdz]) - dBzdx_S*Vx0 + SIXTH*dBzdy_S*dmoments_SW[fsgrids::dmoments::dVxdy]);
      Ey_SW += -HALF*((Bx_W - HALF*dBxdz_W)*(-dmoments_SW[fsgrids::dmoments::dVzdx] - dmoments_SW[fsgrids::dmoments::dVzdz]) - dBxdz_W*Vz0 + SIXTH*dBxdy_W*dmoments_SW[fsgrids::dmoments::dVzdy]);
   #endif
   
   c_z = calculateWaveSpeedXZ(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i, j+1, k,
      Bx_W, Bz_S, dBxdy_W, dBxdz_W, dBzdx_S, dBzdy_S, MINUS, MINUS, minRho, maxRho);
   c_x = c_z;
   az_neg   = max(ZERO,-Vz0 + c_z);
   az_pos   = max(ZERO,+Vz0 + c_z);
   ax_neg   = max(ZERO,-Vx0 + c_x);
   ax_pos   = max(ZERO,+Vx0 + c_x);

   // Ey and characteristic speeds on k-1 neighbour:
   Vx0  = divideIfNonZero(moments_SE[fsgrids::moments::RHOVX], moments_SE[fsgrids::moments::RHO]);
   Vz0  = divideIfNonZero(moments_SE[fsgrids::moments::RHOVZ], moments_SE[fsgrids::moments::RHO]);

   // 1st order terms:
   Real Ey_SE    = Bz_S*Vx0 - Bx_E*Vz0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ey_SE += Parameters::resistivity *
     sqrt((bgb_SE[fsgrids::bgbfield::BGBX]+perb_SE[fsgrids::bfield::PERBX])*
          (bgb_SE[fsgrids::bgbfield::BGBX]+perb_SE[fsgrids::bfield::PERBX]) +
          (bgb_SE[fsgrids::bgbfield::BGBY]+perb_SE[fsgrids::bfield::PERBY])*
          (bgb_SE[fsgrids::bgbfield::BGBY]+perb_SE[fsgrids::bfield::PERBY]) +
          (bgb_SE[fsgrids::bgbfield::BGBZ]+perb_SE[fsgrids::bfield::PERBZ])*
          (bgb_SE[fsgrids::bgbfield::BGBZ]+perb_SE[fsgrids::bfield::PERBZ])
         ) /
     (moments_SE[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_SE[fsgrids::dperb::dPERBxdz]/technicalGrid.DZ - dperb_SE[fsgrids::dperb::dPERBzdx]/technicalGrid.DX);

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ey_SE += EHallGrid.get(i,j,k-1)[fsgrids::ehall::EYHALL_001_011];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ey_SE += EGradPeGrid.get(i,j,k-1)[fsgrids::egradpe::EYGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ey_SE += +HALF*((Bz_S - HALF*dBzdx_S)*(-dmoments_SE[fsgrids::dmoments::dVxdx] + dmoments_SE[fsgrids::dmoments::dVxdz]) - dBzdx_S*Vx0 + SIXTH*dBzdy_S*dmoments_SE[fsgrids::dmoments::dVxdy]);
      Ey_SE += -HALF*((Bx_E + HALF*dBxdz_E)*(-dmoments_SE[fsgrids::dmoments::dVzdx] + dmoments_SE[fsgrids::dmoments::dVzdz]) + dBxdz_E*Vz0 + SIXTH*dBxdy_E*dmoments_SE[fsgrids::dmoments::dVzdy]);
   #endif
   
   c_z = calculateWaveSpeedXZ(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i, j+1, k-1,
      Bx_E, Bz_S, dBxdy_E, dBxdz_E, dBzdx_S, dBzdy_S, MINUS, PLUS, minRho, maxRho);
   c_x = c_z;
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   
   // Ey and characteristic speeds on i-1 neighbour:
   Vz0  = divideIfNonZero(moments_NW[fsgrids::moments::RHOVZ], moments_NW[fsgrids::moments::RHO]);
   Vx0  = divideIfNonZero(moments_NW[fsgrids::moments::RHOVX], moments_NW[fsgrids::moments::RHO]);
   
   // 1st order terms:
   Real Ey_NW    = Bz_N*Vx0 - Bx_W*Vz0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ey_NW += Parameters::resistivity *
     sqrt((bgb_NW[fsgrids::bgbfield::BGBX]+perb_NW[fsgrids::bfield::PERBX])*
          (bgb_NW[fsgrids::bgbfield::BGBX]+perb_NW[fsgrids::bfield::PERBX]) +
          (bgb_NW[fsgrids::bgbfield::BGBY]+perb_NW[fsgrids::bfield::PERBY])*
          (bgb_NW[fsgrids::bgbfield::BGBY]+perb_NW[fsgrids::bfield::PERBY]) +
          (bgb_NW[fsgrids::bgbfield::BGBZ]+perb_NW[fsgrids::bfield::PERBZ])*
          (bgb_NW[fsgrids::bgbfield::BGBZ]+perb_NW[fsgrids::bfield::PERBZ])
         ) /
     (moments_NW[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_NW[fsgrids::dperb::dPERBxdz]/technicalGrid.DZ - dperb_NW[fsgrids::dperb::dPERBzdx]/technicalGrid.DX);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ey_NW += EHallGrid.get(i-1,j,k)[fsgrids::ehall::EYHALL_100_110];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ey_NW += EGradPeGrid.get(i-1,j,k)[fsgrids::egradpe::EYGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ey_NW += +HALF*((Bz_N + HALF*dBzdx_N)*(+dmoments_NW[fsgrids::dmoments::dVxdx] - dmoments_NW[fsgrids::dmoments::dVxdz]) + dBzdx_N*Vx0 + SIXTH*dBzdy_N*dmoments_NW[fsgrids::dmoments::dVxdy]);
      Ey_NW += -HALF*((Bx_W - HALF*dBxdz_W)*(+dmoments_NW[fsgrids::dmoments::dVzdx] - dmoments_NW[fsgrids::dmoments::dVzdz]) - dBxdz_W*Vz0 + SIXTH*dBxdy_W*dmoments_NW[fsgrids::dmoments::dVzdy]);
   #endif
   
   c_z = calculateWaveSpeedXZ(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i-1,j+1,k,
      Bx_W, Bz_N, dBxdy_W, dBxdz_W, dBzdx_N, dBzdy_N, PLUS, MINUS, minRho, maxRho);
   c_x = c_z;
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);

   // Ey and characteristic speeds on i-1,k-1 neighbour:
   Vz0 = divideIfNonZero(moments_NE[fsgrids::moments::RHOVZ], moments_NE[fsgrids::moments::RHO]);
   Vx0 = divideIfNonZero(moments_NE[fsgrids::moments::RHOVX], moments_NE[fsgrids::moments::RHO]);
   
   // 1st order terms:
   Real Ey_NE    = Bz_N*Vx0 - Bx_E*Vz0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ey_NE += Parameters::resistivity *
     sqrt((bgb_NE[fsgrids::bgbfield::BGBX]+perb_NE[fsgrids::bfield::PERBX])*
          (bgb_NE[fsgrids::bgbfield::BGBX]+perb_NE[fsgrids::bfield::PERBX]) +
          (bgb_NE[fsgrids::bgbfield::BGBY]+perb_NE[fsgrids::bfield::PERBY])*
          (bgb_NE[fsgrids::bgbfield::BGBY]+perb_NE[fsgrids::bfield::PERBY]) +
          (bgb_NE[fsgrids::bgbfield::BGBZ]+perb_NE[fsgrids::bfield::PERBZ])*
          (bgb_NE[fsgrids::bgbfield::BGBZ]+perb_NE[fsgrids::bfield::PERBZ])
         ) /
     (moments_NE[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_NE[fsgrids::dperb::dPERBxdz]/technicalGrid.DZ - dperb_NE[fsgrids::dperb::dPERBzdx]/technicalGrid.DX);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ey_NE += EHallGrid.get(i-1,j,k-1)[fsgrids::ehall::EYHALL_101_111];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ey_NE += EGradPeGrid.get(i-1,j,k-1)[fsgrids::egradpe::EYGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ey_NE += +HALF*((Bz_N + HALF*dBzdx_N)*(+dmoments_NE[fsgrids::dmoments::dVxdx] + dmoments_NE[fsgrids::dmoments::dVxdz]) + dBzdx_N*Vx0 + SIXTH*dBzdy_N*dmoments_NE[fsgrids::dmoments::dVxdy]);
      Ey_NE += -HALF*((Bx_E + HALF*dBxdz_E)*(+dmoments_NE[fsgrids::dmoments::dVzdx] + dmoments_NE[fsgrids::dmoments::dVzdz]) + dBxdz_E*Vz0 + SIXTH*dBxdy_E*dmoments_NE[fsgrids::dmoments::dVzdy]);
   #endif
   
   c_z = calculateWaveSpeedXZ(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i-1,j+1,k-1,
      Bx_E, Bz_N, dBxdy_E, dBxdz_E, dBzdx_N, dBzdy_N, PLUS, PLUS, minRho, maxRho);
   c_x = c_z;
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   
   // Calculate properly upwinded edge-averaged Ey:
   efield_SW[fsgrids::efield::EY]  = az_pos*ax_pos*Ey_NE + az_pos*ax_neg*Ey_SE + az_neg*ax_pos*Ey_NW + az_neg*ax_neg*Ey_SW;
   efield_SW[fsgrids::efield::EY] /= ((az_pos+az_neg)*(ax_pos+ax_neg)+EPS);

   if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
      efield_SW[fsgrids::efield::EY] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(perBz_S-perBz_N);
      efield_SW[fsgrids::efield::EY] += az_pos*az_neg/(az_pos+az_neg+EPS)*(perBx_W-perBx_E);
#else
      efield_SW[fsgrids::efield::EY] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((perBz_S-HALF*dperBzdx_S) - (perBz_N+HALF*dperBzdx_N));
      efield_SW[fsgrids::efield::EY] += az_pos*az_neg/(az_pos+az_neg+EPS)*((perBx_W-HALF*dperBxdz_W) - (perBx_E+HALF*dperBxdz_E));
#endif
   }
   
   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      //compute maximum timestep for fieldsolver in this cell (CFL=1)      
      Real max_a=ZERO;
      max_a=max(fabs(az_neg),max_a);
      max_a=max(fabs(az_pos),max_a);
      max_a=max(fabs(ax_neg),max_a);
      max_a=max(fabs(ax_pos),max_a);
      Real min_dx=std::numeric_limits<Real>::max();;
      min_dx=min(min_dx,technicalGrid.DX);
      min_dx=min(min_dx,technicalGrid.DZ);
      //update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (max_a!=ZERO) technicalGrid.get(i,j,k)->maxFsDt=min(technicalGrid.get(i,j,k)->maxFsDt,min_dx/max_a);
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
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   i,
   j,
   k,
   cint& RKCase
) {
   #ifdef DEBUG_FSOLVER
   bool ok = true;
   if (technicalGrid.get(i,j,k) == NULL) ok = false;
   if (technicalGrid.get(i-1,j,k) == NULL) ok = false;
   if (technicalGrid.get(i-1,j-1,k) == NULL) ok = false;
   if (technicalGrid.get(i,j-1,k) == NULL) ok = false;
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
   Real c_x,c_y;                    // Characteristic speeds to xy-directions
   
   // Get read-only pointers to NE,NW,SE,SW states (SW is rw, result is written there):
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_SW = perBGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_SE = perBGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_NE = perBGrid.get(i-1,j-1,k  );
   const std::array<Real, fsgrids::bfield::N_BFIELD> * const perb_NW = perBGrid.get(i  ,j-1,k  );
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_SW = BgBGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_SE = BgBGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_NE = BgBGrid.get(i-1,j-1,k  );
   const std::array<Real, fsgrids::bgbfield::N_BGB> * const bgb_NW = BgBGrid.get(i  ,j-1,k  );
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_SW = momentsGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_SE = momentsGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_NE = momentsGrid.get(i-1,j-1,k  );
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const moments_NW = momentsGrid.get(i  ,j-1,k  );
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_SW = dMomentsGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_SE = dMomentsGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_NE = dMomentsGrid.get(i-1,j-1,k  );
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dmoments_NW = dMomentsGrid.get(i  ,j-1,k  );
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_SW = dPerBGrid.get(i  ,j  ,k  );
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_SE = dPerBGrid.get(i-1,j  ,k  );
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_NE = dPerBGrid.get(i-1,j-1,k  );
   const std::array<Real, fsgrids::dperb::N_DPERB> * const dperb_NW = dPerBGrid.get(i  ,j-1,k  );
   
   const std::array<Real, fsgrids::efield::N_EFIELD> * efield_SW = EGrid.get(i,j,k);
   
   // Fetch needed plasma parameters/derivatives from the four cells:
   Real Bx_S, By_W, By_E, Bx_N, perBx_S, perBy_W, perBy_E, perBx_N;
   Real minRho = std::numeric_limits<Real>::max();
   Real maxRho = std::numeric_limits<Real>::min();
   Bx_S    = perb_SW[fsgrids::bfield::PERBX] + bgb_SW[fsgrids::bgbfield::BGBX];
   By_W    = perb_SW[fsgrids::bfield::PERBY] + bgb_SW[fsgrids::bgbfield::BGBY];
   By_E    = perb_SE[fsgrids::bfield::PERBY] + bgb_SE[fsgrids::bgbfield::BGBY];
   Bx_N    = perb_NW[fsgrids::bfield::PERBX] + bgb_NW[fsgrids::bgbfield::BGBX];
   perBx_S    = perb_SW[fsgrids::bfield::PERBX];
   perBy_W    = perb_SW[fsgrids::bfield::PERBY];
   perBy_E    = perb_SE[fsgrids::bfield::PERBY];
   perBx_N    = perb_NW[fsgrids::bfield::PERBX];
   Vx0  = divideIfNonZero(moments_SW[fsgrids::moments::RHOVX], moments_SW[fsgrids::moments::RHO]);
   Vy0  = divideIfNonZero(moments_SW[fsgrids::moments::RHOVY], moments_SW[fsgrids::moments::RHO]);
   minRho = min(minRho,
               min(moments_SW[fsgrids::moments::RHO],
                  min(moments_SE[fsgrids::moments::RHO],
                     min(moments_NW[fsgrids::moments::RHO],
                           moments_NE[fsgrids::moments::RHO])
                     )
                  )
               );
   maxRho = max(maxRho,
               max(moments_SW[fsgrids::moments::RHO],
                  max(moments_SE[fsgrids::moments::RHO],
                     max(moments_NW[fsgrids::moments::RHO],
                           moments_NE[fsgrids::moments::RHO])
                     )
                  )
               );
   
   creal dBxdy_S = dperb_SW[fsgrids::dperb::dPERBxdy] + dbgb_SW[fsgrids::bgbfield::dBGBxdy];
   creal dBxdz_S = dperb_SW[fsgrids::dperb::dPERBxdz] + dbgb_SW[fsgrids::bgbfield::dBGBxdz];
   creal dBydx_W = dperb_SW[fsgrids::dperb::dPERBydx] + dbgb_SW[fsgrids::bgbfield::dBGBydx];
   creal dBydz_W = dperb_SW[fsgrids::dperb::dPERBydz] + dbgb_SW[fsgrids::bgbfield::dBGBydz];
   creal dBydx_E = dperb_SE[fsgrids::dperb::dPERBydx] + dbgb_SE[fsgrids::bgbfield::dBGBydx];
   creal dBydz_E = dperb_SE[fsgrids::dperb::dPERBydz] + dbgb_SE[fsgrids::bgbfield::dBGBydz];
   creal dBxdy_N = dperb_NW[fsgrids::dperb::dPERBxdy] + dbgb_NW[fsgrids::bgbfield::dBGBxdy];
   creal dBxdz_N = dperb_NW[fsgrids::dperb::dPERBxdz] + dbgb_NW[fsgrids::bgbfield::dBGBxdz];
   creal dperBxdy_S = dperb_SW[fsgrids::dperb::dPERBxdy];
   creal dperBxdy_N = dperb_NW[fsgrids::dperb::dPERBxdy];
   creal dperBydx_W = dperb_SW[fsgrids::dperb::dPERBydx];
   creal dperBydx_E = dperb_SE[fsgrids::dperb::dPERBydx];
   
   // Ez and characteristic speeds on SW cell:
   // 1st order terms:
   Real Ez_SW = Bx_S*Vy0 - By_W*Vx0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ez_SW += Parameters::resistivity *
     sqrt((bgb_SW[fsgrids::bgbfield::BGBX]+perb_SW[fsgrids::bfield::PERBX])*
          (bgb_SW[fsgrids::bgbfield::BGBX]+perb_SW[fsgrids::bfield::PERBX]) +
          (bgb_SW[fsgrids::bgbfield::BGBY]+perb_SW[fsgrids::bfield::PERBY])*
          (bgb_SW[fsgrids::bgbfield::BGBY]+perb_SW[fsgrids::bfield::PERBY]) +
          (bgb_SW[fsgrids::bgbfield::BGBZ]+perb_SW[fsgrids::bfield::PERBZ])*
          (bgb_SW[fsgrids::bgbfield::BGBZ]+perb_SW[fsgrids::bfield::PERBZ])
         ) /
     (moments_SW[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_SW[fsgrids::dperb::dPERBydx]/technicalGrid.DX - dperb_SW[fsgrids::dperb::dPERBxdy]/technicalGrid.DY);
   
   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ez_SW += EHallGrid.get(i,j,k)[fsgrids::ehall::EZHALL_000_001];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ez_SW += EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EZGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_SW  += +HALF*((Bx_S - HALF*dBxdy_S)*(-dmoments_SW[fsgrids::dmoments::dVydx] - dmoments_SW[fsgrids::dmoments::dVydy]) - dBxdy_S*Vy0 + SIXTH*dBxdz_S*dmoments_SW[fsgrids::dmoments::dVydz]);
      Ez_SW  += -HALF*((By_W - HALF*dBydx_W)*(-dmoments_SW[fsgrids::dmoments::dVxdx] - dmoments_SW[fsgrids::dmoments::dVxdy]) - dBydx_W*Vx0 + SIXTH*dBydz_W*dmoments_SW[fsgrids::dmoments::dVxdz]);
   #endif
   
   // Calculate maximum wave speed (fast magnetosonic speed) on SW cell. In order 
   // to get Alfven speed we need to calculate some reconstruction coeff. for Bz:
   c_x = calculateWaveSpeedXY(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i, j, k+1,
      Bx_S, By_W, dBxdy_S, dBxdz_S, dBydx_W, dBydz_W, MINUS, MINUS, minRho, maxRho);
   c_y = c_x;
   ax_neg   = max(ZERO,-Vx0 + c_x);
   ax_pos   = max(ZERO,+Vx0 + c_x);
   ay_neg   = max(ZERO,-Vy0 + c_y);
   ay_pos   = max(ZERO,+Vy0 + c_y);

   // Ez and characteristic speeds on SE (i-1) cell:
   Vx0  = divideIfNonZero(moments_SE[fsgrids::moments::RHOVX], moments_SE[fsgrids::moments::RHO]);
   Vy0  = divideIfNonZero(moments_SE[fsgrids::moments::RHOVY], moments_SE[fsgrids::moments::RHO]);
   
   // 1st order terms:
   Real Ez_SE = Bx_S*Vy0 - By_E*Vx0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ez_SE += Parameters::resistivity *
     sqrt((bgb_SE[fsgrids::bgbfield::BGBX]+perb_SE[fsgrids::bfield::PERBX])*
          (bgb_SE[fsgrids::bgbfield::BGBX]+perb_SE[fsgrids::bfield::PERBX]) +
          (bgb_SE[fsgrids::bgbfield::BGBY]+perb_SE[fsgrids::bfield::PERBY])*
          (bgb_SE[fsgrids::bgbfield::BGBY]+perb_SE[fsgrids::bfield::PERBY]) +
          (bgb_SE[fsgrids::bgbfield::BGBZ]+perb_SE[fsgrids::bfield::PERBZ])*
          (bgb_SE[fsgrids::bgbfield::BGBZ]+perb_SE[fsgrids::bfield::PERBZ])
         ) /
     (moments_SE[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_SE[fsgrids::dperb::dPERBydx]/technicalGrid.DX - dperb_SE[fsgrids::dperb::dPERBxdy]/technicalGrid.DY);
   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ez_SE += EHallGrid.get(i-1,j,k)[fsgrids::ehall::EZHALL_100_101];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ez_SE += EGradPeGrid.get(i-1,j,k)[fsgrids::egradpe::EZGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_SE  += +HALF*((Bx_S - HALF*dBxdy_S)*(+dmoments_SE[fsgrids::dmoments::dVydx] - dmoments_SE[fsgrids::dmoments::dVydy]) - dBxdy_S*Vy0 + SIXTH*dBxdz_S*dmoments_SE[fsgrids::dmoments::dVydz]);
      Ez_SE  += -HALF*((By_E + HALF*dBydx_E)*(+dmoments_SE[fsgrids::dmoments::dVxdx] - dmoments_SE[fsgrids::dmoments::dVxdy]) + dBydx_E*Vx0 + SIXTH*dBydz_E*dmoments_SE[fsgrids::dmoments::dVxdz]);
   #endif
   
   c_x = calculateWaveSpeedXY(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i-1,j,k+1,
      Bx_S, By_E, dBxdy_S, dBxdz_S, dBydx_E, dBydz_E, PLUS, MINUS, minRho, maxRho);
   c_y = c_x;
   ax_neg = max(ax_neg,-Vx0 + c_x);
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);

   // Ez and characteristic speeds on NW (j-1) cell:
   Vx0  = divideIfNonZero(moments_NW[fsgrids::moments::RHOVX], moments_NW[fsgrids::moments::RHO]);
   Vy0  = divideIfNonZero(moments_NW[fsgrids::moments::RHOVY], moments_NW[fsgrids::moments::RHO]);

   // 1st order terms:
   Real Ez_NW = Bx_N*Vy0 - By_W*Vx0;

   // Resistive term
   // FIXME this does not include RK stepping
   Ez_NW += Parameters::resistivity *
     sqrt((bgb_NW[fsgrids::bgbfield::BGBX]+perb_NW[fsgrids::bfield::PERBX])*
          (bgb_NW[fsgrids::bgbfield::BGBX]+perb_NW[fsgrids::bfield::PERBX]) +
          (bgb_NW[fsgrids::bgbfield::BGBY]+perb_NW[fsgrids::bfield::PERBY])*
          (bgb_NW[fsgrids::bgbfield::BGBY]+perb_NW[fsgrids::bfield::PERBY]) +
          (bgb_NW[fsgrids::bgbfield::BGBZ]+perb_NW[fsgrids::bfield::PERBZ])*
          (bgb_NW[fsgrids::bgbfield::BGBZ]+perb_NW[fsgrids::bfield::PERBZ])
         ) /
     (moments_NW[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_NW[fsgrids::dperb::dPERBydx]/technicalGrid.DX - dperb_NW[fsgrids::dperb::dPERBxdy]/technicalGrid.DY);

   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ez_NW += EHallGrid.get(i,j-1,k)[fsgrids::ehall::EZHALL_010_011];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ez_NW += EGradPeGrid.get(i,j-1,k)[fsgrids::egradpe::EZGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_NW  += +HALF*((Bx_N + HALF*dBxdy_N)*(-dmoments_NW[fsgrids::dmoments::dVydx] + dmoments_NW[fsgrids::dmoments::dVydy]) + dBxdy_N*Vy0 + SIXTH*dBxdz_N*dmoments_NW[fsgrids::dmoments::dVydz]);
      Ez_NW  += -HALF*((By_W - HALF*dBydx_W)*(-dmoments_NW[fsgrids::dmoments::dVxdx] + dmoments_NW[fsgrids::dmoments::dVxdy]) - dBydx_W*Vx0 + SIXTH*dBydz_W*dmoments_NW[fsgrids::dmoments::dVxdz]);
   #endif
   
   c_x = calculateWaveSpeedXY(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i, j-1, k+1,
      Bx_N, By_W, dBxdy_N, dBxdz_N, dBydx_W, dBydz_W, MINUS, PLUS, minRho, maxRho);
   c_y = c_x;
   ax_neg = max(ax_neg,-Vx0 + c_x); 
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);
   
   // Ez and characteristic speeds on NE (i-1,j-1) cell:
   Vx0  = divideIfNonZero(moments_NE[fsgrids::moments::RHOVX], moments_NE[fsgrids::moments::RHO]);
   Vy0  = divideIfNonZero(moments_NE[fsgrids::moments::RHOVY], moments_NE[fsgrids::moments::RHO]);
   
   // 1st order terms:
   Real Ez_NE = Bx_N*Vy0 - By_E*Vx0;
   
   // Resistive term
   // FIXME this does not include RK stepping
   Ez_NE += Parameters::resistivity *
     sqrt((bgb_NE[fsgrids::bgbfield::BGBX]+perb_NE[fsgrids::bfield::PERBX])*
          (bgb_NE[fsgrids::bgbfield::BGBX]+perb_NE[fsgrids::bfield::PERBX]) +
          (bgb_NE[fsgrids::bgbfield::BGBY]+perb_NE[fsgrids::bfield::PERBY])*
          (bgb_NE[fsgrids::bgbfield::BGBY]+perb_NE[fsgrids::bfield::PERBY]) +
          (bgb_NE[fsgrids::bgbfield::BGBZ]+perb_NE[fsgrids::bfield::PERBZ])*
          (bgb_NE[fsgrids::bgbfield::BGBZ]+perb_NE[fsgrids::bfield::PERBZ])
         ) /
     (moments_NE[fsgrids::moments::RHO]*physicalconstants::CHARGE) /
     physicalconstants::MU_0 *
     (dperb_NE[fsgrids::dperb::dPERBydx]/technicalGrid.DX - dperb_NE[fsgrids::dperb::dPERBxdy]/technicalGrid.DY);
   
   // Hall term
   if(Parameters::ohmHallTerm > 0) {
      Ez_NE += EHallGrid.get(i-1,j-1,k)[fsgrids::ehall::EZHALL_110_111];
   }
   
   // Electron pressure gradient term
   if(Parameters::ohmGradPeTerm > 0) {
      Ez_NE += EGradPeGrid.get(i-1,j-1,k)[fsgrids::egradpe::EZGRADPE];
   }
   
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_NE  += +HALF*((Bx_N + HALF*dBxdy_N)*(+dmoments_NE[fsgrids::dmoments::dVydx] + dmoments_NE[fsgrids::dmoments::dVydy]) + dBxdy_N*Vy0 + SIXTH*dBxdz_N*dmoments_NE[fsgrids::dmoments::dVydz]);
      Ez_NE  += -HALF*((By_E + HALF*dBydx_E)*(+dmoments_NE[fsgrids::dmoments::dVxdx] + dmoments_NE[fsgrids::dmoments::dVxdy]) + dBydx_E*Vx0 + SIXTH*dBydz_E*dmoments_NE[fsgrids::dmoments::dVxdz]);
   #endif
   
   c_x = calculateWaveSpeedXY(
      perBGrid,
      momentsGrid,
      dPerBGrid,
      dMomentsGrid,
      BgBGrid,
      i, j, k,
      i-1, j-1, k+1,
      Bx_N, By_E, dBxdy_N, dBxdz_N, dBydx_E, dBydz_E, PLUS, PLUS, minRho, maxRho);
   c_y = c_x;
   ax_neg = max(ax_neg,-Vx0 + c_x);
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);

   // Calculate properly upwinded edge-averaged Ez:
   efield_SW[fsgrids::efield::EZ] = ax_pos*ay_pos*Ez_NE + ax_pos*ay_neg*Ez_SE + ax_neg*ay_pos*Ez_NW + ax_neg*ay_neg*Ez_SW;
   efield_SW[fsgrids::efield::EZ] /= ((ax_pos+ax_neg)*(ay_pos+ay_neg)+EPS);

   if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
      efield_SW[fsgrids::efield::EZ] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(perBx_S-perBx_N);
      efield_SW[fsgrids::efield::EZ] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(perBy_W-perBy_E);
#else
      efield_SW[fsgrids::efield::EZ] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((perBx_S-HALF*dperBxdy_S) - (perBx_N+HALF*dperBxdy_N));
      efield_SW[fsgrids::efield::EZ] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((perBy_W-HALF*dperBydx_W) - (perBy_E+HALF*dperBydx_E));
#endif
   }
   
   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      //compute maximum timestep for fieldsolver in this cell (CFL=1)      
      Real max_a=ZERO;
      max_a=max(fabs(ay_neg),max_a);
      max_a=max(fabs(ay_pos),max_a);
      max_a=max(fabs(ax_neg),max_a);
      max_a=max(fabs(ax_pos),max_a);
      Real min_dx=std::numeric_limits<Real>::max();;
      min_dx=min(min_dx,technicalGrid.DX);
      min_dx=min(min_dx,technicalGrid.DY);
      //update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if(max_a!=ZERO) technicalGrid.get(i,j,k)->maxFsDt=min(technicalGrid.get(i,j,k)->maxFsDt,min_dx/max_a);
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
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   cint i,
   cint j,
   cint k,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {
   cuint cellSysBoundaryFlag = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   
   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
   
   cuint cellSysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;
   
   if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && (cellSysBoundaryLayer != 1)) {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, 0);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, 1);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, 2);
   } else {
      calculateEdgeElectricFieldX(
         perBGrid,
         EGrid,
         EHallGrid,
         EGradPeGrid,
         momentsGrid,
         dPerBGrid,
         dMomentsGrid,
         BgBGrid,
         technicalGrid,
         i,
         j,
         k,
         RKCase
      );
      calculateEdgeElectricFieldY(
         perBGrid,
         EGrid,
         EHallGrid,
         EGradPeGrid,
         momentsGrid,
         dPerBGrid,
         dMomentsGrid,
         BgBGrid,
         technicalGrid,
         i,
         j,
         k,
         RKCase
      );
      calculateEdgeElectricFieldZ(
         perBGrid,
         EGrid,
         EHallGrid,
         EGradPeGrid,
         momentsGrid,
         dPerBGrid,
         dMomentsGrid,
         BgBGrid,
         technicalGrid,
         i,
         j,
         k,
         RKCase
      );
   }
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
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {
   int timer;
   const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::start("Calculate upwinded electric field");
   if(P::ohmHallTerm > 0) {
      EHallGrid.updateGhostCells();
   }
   if(P::ohmGradPeTerm > 0) {
      EGradPeGrid.updateGhostCells();
   }
   if(P::ohmHallTerm == 0 && P::ohmGradPeTerm == 0) {
      dPerBGrid.updateGhostCells();
      dMomentsGrid.updateGhostCells();
   }
   
   // Calculate upwinded electric field on inner cells
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);
   #pragma omp parallel for collapse(3)
   for (uint k=0; k<gridDims[2]; k++) {
      for (uint j=0; j<gridDims[1]; j++) {
         for (uint i=0; i<gridDims[0]; i++) {
            if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
               calculateElectricField(
                  perBGrid,
                  EGrid,
                  EHallGrid,
                  EGradPeGrid,
                  momentsGrid,
                  dPerBGrid,
                  dMomentsGrid,
                  BgBGrid,
                  technicalGrid,
                  i,
                  j,
                  k,
                  sysBoundaries,
                  RKCase
               );
            } else { // RKCase == RK_ORDER2_STEP1
               calculateElectricField(
                  perBDt2Grid,
                  EDt2Grid,
                  EHallGrid,
                  EGradPeGrid,
                  momentsDt2Grid,
                  dPerBGrid,
                  dMomentsGrid,
                  BgBGrid,
                  technicalGrid,
                  i,
                  j,
                  k,
                  sysBoundaries,
                  RKCase
               );
            }
         }
      }
   }
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   // Exchange electric field with neighbouring processes
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      EGrid.updateGhostCells();
   } else { 
      EDt2Grid.updateGhostCells();
   }
   
   phiprof::stop("Calculate upwinded electric field",N_cells,"Spatial Cells");
}
