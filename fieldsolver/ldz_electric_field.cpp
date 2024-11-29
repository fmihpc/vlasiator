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

#include <algorithm>
#include <cstdlib>

#include "fs_common.h"
#include "ldz_electric_field.hpp"

#ifdef DEBUG_VLASIATOR
#define DEBUG_FSOLVER
#endif

namespace pc = physicalconstants;
using namespace std;

struct Wavespeeds {
   Real alfven = 0.0;
   Real sound = 0.0;
   Real whistler = 0.0;

   // Effective wave speeds for advection and CFL calculation
   // Note that these are calculated as if the plasma is purely made up of hydrogen, which
   // is a reasonable approximation if it is proton-dominant.
   // Simulations which predominantly contain heavier ion species will have to change this!
   //
   // See
   // T E Stringer, Low-frequency waves in an unbounded plasma
   // Journal of Nuclear Energy. Part C, Plasma Physics, Accelerators, Thermonuclear Research
   // Volume 5, Number 2, page 89
   // Section BC in Fig. 1 and Table 1, referenced later e.g. by Gary, for the current one,
   // and the below for the heavy ions
   // https://www.ann-geophys.net/26/1605/2008/  (Whistler waves)
   // and
   // http://iopscience.iop.org/article/10.1088/0253-6102/43/2/026/meta (Alfven waves)
   // for details.
   Wavespeeds(Real bmag2, Real rhom, Real p11, Real p22, Real p33, const std::array<Real, 3>& gridSpacing)
       : alfven(sqrt(divideIfNonZero(bmag2, pc::MU_0 * rhom))),
         sound(sqrt(divideIfNonZero(p11 + p22 + p33, 2.0 * rhom))),
         whistler(Parameters::ohmHallTerm > 0
                      ? alfven * (1 + divideIfNonZero(2 * M_PI * M_PI * pc::MASS_PROTON * pc::MASS_PROTON,
                                                      gridSpacing[0] * gridSpacing[0] * rhom * pc::CHARGE * pc::CHARGE *
                                                          pc::MU_0) /
                                          sqrt(1 + divideIfNonZero(M_PI * M_PI * pc::MASS_PROTON * pc::MASS_PROTON,
                                                                   gridSpacing[0] * gridSpacing[0] * rhom * pc::CHARGE *
                                                                       pc::CHARGE * pc::MU_0)))
                      : 0.0) {}

   Real minVelocity() const {
      return min(Parameters::maxWaveVelocity, sqrt(alfven * alfven + sound * sound) + whistler);
   }

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
    */
   Real cflSpeed(Real v0, Real v1) const {
      const Real v = sqrt(v0 * v0 + v1 * v1);
      const Real vMS = sqrt(alfven * alfven + sound * sound);
      return max(v + vMS, v + whistler);
   }
};

struct Limits {
   Real min = 0.0;
   Real max = 0.0;
};

Real clampNegativeToZero(Real v) { return std::clamp(v, 0.0, v); }

/*! \brief Low-level helper function.
 *
 * Computes the magnetosonic speed in the YZ plane. Used in upwinding the electric field X component,
 * at the interface between cell (i,j,k) and (nbi,nbj,nbk).
 *
 * Expects that the correct RHO and B fields are being passed, depending on the
 * stage of the Runge-Kutta time stepping method.
 *
 * If fields are not propagated, returns 0.0 as there is no information propagating.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param nbi,nbj,nbk fsGrid cell coordinates for the adjacent cell
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
 * \param ret_vA Alfven speed returned
 * \param ret_vS Sound speed returned
 * \param ret_vW Whistler speed returned
 */
Wavespeeds
calculateWaveSpeedYZ(fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                     const std::array<Real, 3>& gridSpacing, cint i, cint j, cint k, cint nbi, cint nbj, cint nbk,
                     const Real& By, const Real& Bz, const Real& dBydx, const Real& dBydz, const Real& dBzdx,
                     const Real& dBzdy, const Real& ydir, const Real& zdir, Limits rhomLimits) {
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb = perBGrid.get(i, j, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* nbr_perb = perBGrid.get(nbi, nbj, nbk);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments = momentsGrid.get(i, j, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments = dMomentsGrid.get(i, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb = dPerBGrid.get(i, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* nbr_dperb = dPerBGrid.get(nbi, nbj, nbk);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb = BgBGrid.get(i, j, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* nbr_bgb = BgBGrid.get(nbi, nbj, nbk);

   const Real A_0 = HALF * (nbr_perb->at(fsgrids::bfield::PERBX) + nbr_bgb->at(fsgrids::bgbfield::BGBX) +
                            perb->at(fsgrids::bfield::PERBX) + bgb->at(fsgrids::bgbfield::BGBX));
   const Real A_X = (nbr_perb->at(fsgrids::bfield::PERBX) + nbr_bgb->at(fsgrids::bgbfield::BGBX)) -
                    (perb->at(fsgrids::bfield::PERBX) + bgb->at(fsgrids::bgbfield::BGBX));
   const Real rhom =
       std::clamp(moments->at(fsgrids::moments::RHOM) + ydir * HALF * dmoments->at(fsgrids::dmoments::drhomdy) +
                      zdir * HALF * dmoments->at(fsgrids::dmoments::drhomdz),
                  rhomLimits.min, rhomLimits.max);
   const Real p11 =
       clampNegativeToZero(moments->at(fsgrids::moments::P_11) + ydir * HALF * dmoments->at(fsgrids::dmoments::dp11dy) +
                           zdir * HALF * dmoments->at(fsgrids::dmoments::dp11dz));
   const Real p22 =
       clampNegativeToZero(moments->at(fsgrids::moments::P_22) + ydir * HALF * dmoments->at(fsgrids::dmoments::dp22dy) +
                           zdir * HALF * dmoments->at(fsgrids::dmoments::dp22dz));
   const Real p33 =
       clampNegativeToZero(moments->at(fsgrids::moments::P_33) + ydir * HALF * dmoments->at(fsgrids::dmoments::dp33dy) +
                           zdir * HALF * dmoments->at(fsgrids::dmoments::dp33dz));

   const Real A_Y = nbr_dperb->at(fsgrids::dperb::dPERBxdy) + nbr_bgb->at(fsgrids::bgbfield::dBGBxdy) +
                    dperb->at(fsgrids::dperb::dPERBxdy) + bgb->at(fsgrids::bgbfield::dBGBxdy);
   const Real A_XY = nbr_dperb->at(fsgrids::dperb::dPERBxdy) + nbr_bgb->at(fsgrids::bgbfield::dBGBxdy) -
                     (dperb->at(fsgrids::dperb::dPERBxdy) + bgb->at(fsgrids::bgbfield::dBGBxdy));
   const Real A_Z = nbr_dperb->at(fsgrids::dperb::dPERBxdz) + nbr_bgb->at(fsgrids::bgbfield::dBGBxdz) +
                    dperb->at(fsgrids::dperb::dPERBxdz) + bgb->at(fsgrids::bgbfield::dBGBxdz);
   const Real A_XZ = nbr_dperb->at(fsgrids::dperb::dPERBxdz) + nbr_bgb->at(fsgrids::bgbfield::dBGBxdz) -
                     (dperb->at(fsgrids::dperb::dPERBxdz) + bgb->at(fsgrids::bgbfield::dBGBxdz));

   const Real bx2 = (A_0 + ydir * HALF * A_Y + zdir * HALF * A_Z) * (A_0 + ydir * HALF * A_Y + zdir * HALF * A_Z) +
                    TWELWTH * (A_X + ydir * HALF * A_XY + zdir * HALF * A_XZ) *
                        (A_X + ydir * HALF * A_XY + zdir * HALF * A_XZ);                               // OK
   const Real by2 = (By + zdir * HALF * dBydz) * (By + zdir * HALF * dBydz) + TWELWTH * dBydx * dBydx; // OK
   const Real bz2 = (Bz + ydir * HALF * dBzdy) * (Bz + ydir * HALF * dBzdy) + TWELWTH * dBzdx * dBzdx; // OK

   return Wavespeeds(bx2 + by2 + bz2, rhom, p11, p22, p33, gridSpacing);
}

/*! \brief Low-level helper function.
 *
 * Computes the magnetosonic speed in the XZ plane. Used in upwinding the electric field Y component,
 * at the interface between cell (i,j,k) and (nbi,nbj,nbk).
 *
 * Expects that the correct RHO and B fields are being passed, depending on the stage of the Runge-Kutta time stepping
 * method.
 *
 * If fields are not propagated, returns 0.0 as there is no information propagating.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param nbi,nbj,nbk fsGrid cell coordinates for the adjacent cell
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
 * \param ret_vA Alfven speed returned
 * \param ret_vS Sound speed returned
 * \param ret_vW Whistler speed returned
 */
Wavespeeds
calculateWaveSpeedXZ(fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                     const std::array<Real, 3>& gridSpacing, cint i, cint j, cint k, cint nbi, cint nbj, cint nbk,
                     const Real& Bx, const Real& Bz, const Real& dBxdy, const Real& dBxdz, const Real& dBzdx,
                     const Real& dBzdy, const Real& xdir, const Real& zdir, Limits rhomLimits) {
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb = perBGrid.get(i, j, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* nbr_perb = perBGrid.get(nbi, nbj, nbk);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments = momentsGrid.get(i, j, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments = dMomentsGrid.get(i, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb = dPerBGrid.get(i, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* nbr_dperb = dPerBGrid.get(nbi, nbj, nbk);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb = BgBGrid.get(i, j, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* nbr_bgb = BgBGrid.get(nbi, nbj, nbk);

   const Real B_0 = HALF * (nbr_perb->at(fsgrids::bfield::PERBY) + nbr_bgb->at(fsgrids::bgbfield::BGBY) +
                            perb->at(fsgrids::bfield::PERBY) + bgb->at(fsgrids::bgbfield::BGBY));
   const Real B_Y = (nbr_perb->at(fsgrids::bfield::PERBY) + nbr_bgb->at(fsgrids::bgbfield::BGBY)) -
                    (perb->at(fsgrids::bfield::PERBY) + bgb->at(fsgrids::bgbfield::BGBY));
   const Real rhom =
       std::clamp(moments->at(fsgrids::moments::RHOM) + xdir * HALF * dmoments->at(fsgrids::dmoments::drhomdx) +
                      zdir * HALF * dmoments->at(fsgrids::dmoments::drhomdz),
                  rhomLimits.min, rhomLimits.max);
   const Real p11 =
       clampNegativeToZero(moments->at(fsgrids::moments::P_11) + xdir * HALF * dmoments->at(fsgrids::dmoments::dp11dx) +
                           zdir * HALF * dmoments->at(fsgrids::dmoments::dp11dz));
   const Real p22 =
       clampNegativeToZero(moments->at(fsgrids::moments::P_22) + xdir * HALF * dmoments->at(fsgrids::dmoments::dp22dx) +
                           zdir * HALF * dmoments->at(fsgrids::dmoments::dp22dz));
   const Real p33 =
       clampNegativeToZero(moments->at(fsgrids::moments::P_33) + xdir * HALF * dmoments->at(fsgrids::dmoments::dp33dx) +
                           zdir * HALF * dmoments->at(fsgrids::dmoments::dp33dz));

   const Real B_X = nbr_dperb->at(fsgrids::dperb::dPERBydx) + nbr_bgb->at(fsgrids::bgbfield::dBGBydx) +
                    dperb->at(fsgrids::dperb::dPERBydx) + bgb->at(fsgrids::bgbfield::dBGBydx);
   const Real B_XY = nbr_dperb->at(fsgrids::dperb::dPERBydx) + nbr_bgb->at(fsgrids::bgbfield::dBGBydx) -
                     (dperb->at(fsgrids::dperb::dPERBydx) + bgb->at(fsgrids::bgbfield::dBGBydx));
   const Real B_Z = nbr_dperb->at(fsgrids::dperb::dPERBydz) + nbr_bgb->at(fsgrids::bgbfield::dBGBydz) +
                    dperb->at(fsgrids::dperb::dPERBydz) + bgb->at(fsgrids::bgbfield::dBGBydz);
   const Real B_YZ = nbr_dperb->at(fsgrids::dperb::dPERBydz) + nbr_bgb->at(fsgrids::bgbfield::dBGBydz) -
                     (dperb->at(fsgrids::dperb::dPERBydz) + bgb->at(fsgrids::bgbfield::dBGBydz));

   const Real by2 = (B_0 + xdir * HALF * B_X + zdir * HALF * B_Z) * (B_0 + xdir * HALF * B_X + zdir * HALF * B_Z) +
                    TWELWTH * (B_Y + xdir * HALF * B_XY + zdir * HALF * B_YZ) *
                        (B_Y + xdir * HALF * B_XY + zdir * HALF * B_YZ);                               // OK
   const Real bx2 = (Bx + zdir * HALF * dBxdz) * (Bx + zdir * HALF * dBxdz) + TWELWTH * dBxdy * dBxdy; // OK
   const Real bz2 = (Bz + xdir * HALF * dBzdx) * (Bz + xdir * HALF * dBzdx) + TWELWTH * dBzdy * dBzdy; // OK

   return Wavespeeds(bx2 + by2 + bz2, rhom, p11, p22, p33, gridSpacing);
}

/*! \brief Low-level helper function.
 *
 * Computes the magnetosonic speed in the XY plane. Used in upwinding the electric field Z component,
 * at the interface between cell (i,j,k) and (nbi,nbj,nbk).
 *
 * Expects that the correct RHO and B fields are being passed, depending on the stage of the Runge-Kutta time stepping
 * method.
 *
 * If fields are not propagated, returns 0.0 as there is no information propagating.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param nbi,nbj,nbk fsGrid cell coordinates for the adjacent cell
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
 * \param ret_vA Alfven speed returned
 * \param ret_vS Sound speed returned
 * \param ret_vW Whistler speed returned
 */
Wavespeeds
calculateWaveSpeedXY(fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
                     fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                     const std::array<Real, 3>& gridSpacing, cint i, cint j, cint k, cint nbi, cint nbj, cint nbk,
                     const Real& Bx, const Real& By, const Real& dBxdy, const Real& dBxdz, const Real& dBydx,
                     const Real& dBydz, const Real& xdir, const Real& ydir, Limits rhomLimits) {
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb = perBGrid.get(i, j, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* nbr_perb = perBGrid.get(nbi, nbj, nbk);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments = momentsGrid.get(i, j, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments = dMomentsGrid.get(i, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb = dPerBGrid.get(i, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* nbr_dperb = dPerBGrid.get(nbi, nbj, nbk);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb = BgBGrid.get(i, j, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* nbr_bgb = BgBGrid.get(nbi, nbj, nbk);

   const Real C_0 = HALF * (nbr_perb->at(fsgrids::bfield::PERBZ) + nbr_bgb->at(fsgrids::bgbfield::BGBZ) +
                            perb->at(fsgrids::bfield::PERBZ) + bgb->at(fsgrids::bgbfield::BGBZ));
   const Real C_Z = (nbr_perb->at(fsgrids::bfield::PERBZ) + nbr_bgb->at(fsgrids::bgbfield::BGBZ)) -
                    (perb->at(fsgrids::bfield::PERBZ) + bgb->at(fsgrids::bgbfield::BGBZ));
   const Real rhom =
       std::clamp(moments->at(fsgrids::moments::RHOM) + xdir * HALF * dmoments->at(fsgrids::dmoments::drhomdx) +
                      ydir * HALF * dmoments->at(fsgrids::dmoments::drhomdy),
                  rhomLimits.min, rhomLimits.max);
   const Real p11 =
       clampNegativeToZero(moments->at(fsgrids::moments::P_11) + xdir * HALF * dmoments->at(fsgrids::dmoments::dp11dx) +
                           ydir * HALF * dmoments->at(fsgrids::dmoments::dp11dy));
   const Real p22 =
       clampNegativeToZero(moments->at(fsgrids::moments::P_22) + xdir * HALF * dmoments->at(fsgrids::dmoments::dp22dx) +
                           ydir * HALF * dmoments->at(fsgrids::dmoments::dp22dy));
   const Real p33 =
       clampNegativeToZero(moments->at(fsgrids::moments::P_33) + xdir * HALF * dmoments->at(fsgrids::dmoments::dp33dx) +
                           ydir * HALF * dmoments->at(fsgrids::dmoments::dp33dy));

   const Real C_X = nbr_dperb->at(fsgrids::dperb::dPERBzdx) + nbr_bgb->at(fsgrids::bgbfield::dBGBzdx) +
                    dperb->at(fsgrids::dperb::dPERBzdx) + bgb->at(fsgrids::bgbfield::dBGBzdx);
   const Real C_XZ = nbr_dperb->at(fsgrids::dperb::dPERBzdx) + nbr_bgb->at(fsgrids::bgbfield::dBGBzdx) -
                     (dperb->at(fsgrids::dperb::dPERBzdx) + bgb->at(fsgrids::bgbfield::dBGBzdx));
   const Real C_Y = nbr_dperb->at(fsgrids::dperb::dPERBzdy) + nbr_bgb->at(fsgrids::bgbfield::dBGBzdy) +
                    dperb->at(fsgrids::dperb::dPERBzdy) + bgb->at(fsgrids::bgbfield::dBGBzdy);
   const Real C_YZ = nbr_dperb->at(fsgrids::dperb::dPERBzdy) + nbr_bgb->at(fsgrids::bgbfield::dBGBzdy) -
                     (dperb->at(fsgrids::dperb::dPERBzdy) + bgb->at(fsgrids::bgbfield::dBGBzdy));

   const Real bz2 =
       (C_0 + xdir * HALF * C_X + ydir * HALF * C_Y) * (C_0 + xdir * HALF * C_X + ydir * HALF * C_Y) +
       TWELWTH * (C_Z + xdir * HALF * C_XZ + ydir * HALF * C_YZ) * (C_Z + xdir * HALF * C_XZ + ydir * HALF * C_YZ);
   const Real bx2 = (Bx + ydir * HALF * dBxdy) * (Bx + ydir * HALF * dBxdy) + TWELWTH * dBxdz * dBxdz;
   const Real by2 = (By + xdir * HALF * dBydx) * (By + xdir * HALF * dBydx) + TWELWTH * dBydz * dBydz;

   return Wavespeeds(bx2 + by2 + bz2, rhom, p11, p22, p33, gridSpacing);
}

void fsdebugCheck([[maybe_unused]] fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
                  [[maybe_unused]] int32_t i, [[maybe_unused]] int32_t j, [[maybe_unused]] int32_t k) {
#ifdef DEBUG_FSOLVER
   const bool ok = technicalGrid.get(i, j, k) != NULL && technicalGrid.get(i - 1, j, k) != NULL &&
                   technicalGrid.get(i - 1, j - 1, k) != NULL && technicalGrid.get(i, j - 1, k) != NULL;
   if (!ok) {
      cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << std::endl;
      exit(1);
   }
#endif
}

Limits getRhomLimits(const std::array<std::array<Real, fsgrids::moments::N_MOMENTS>, 4>& moments) {
   Real minRhom = std::numeric_limits<Real>::max();
   Real maxRhom = std::numeric_limits<Real>::min();
   for (const auto& m : moments) {
      const auto rhom = m[fsgrids::moments::RHOM];
      minRhom = min(minRhom, rhom);
      maxRhom = max(maxRhom, rhom);
   }
   return {
       minRhom,
       maxRhom,
   };
}

Real resistiveTerm(const auto& bgb, const auto& perb, const auto& moments, const auto& dperb,
                   std::array<Real, 2> indices, std::array<Real, 2> spacing) {
   const auto ax = bgb[fsgrids::bgbfield::BGBX];
   const auto bx = perb[fsgrids::bfield::PERBX];
   const auto ay = bgb[fsgrids::bgbfield::BGBY];
   const auto by = perb[fsgrids::bfield::PERBY];
   const auto az = bgb[fsgrids::bgbfield::BGBZ];
   const auto bz = perb[fsgrids::bfield::PERBZ];

   return Parameters::resistivity * sqrt((ax + bx) * (ax + bx) + (ay + by) * (ay + by) + (az + bz) * (az + bz)) /
          moments[fsgrids::moments::RHOQ] / physicalconstants::MU_0 *
          (dperb[indices[0]] / spacing[0] - dperb[indices[1]] / spacing[1]);
}

struct FieldValues {
   Real E_NE = 0.0;
   Real E_SE = 0.0;
   Real E_NW = 0.0;
   Real E_SW = 0.0;
   Real apos = 0.0;
   Real aneg = 0.0;
   Real bpos = 0.0;
   Real bneg = 0.0;
   Real perB_S = 0.0;
   Real perB_N = 0.0;
   Real perB_W = 0.0;
   Real perB_E = 0.0;
   Real dperB_S = 0.0;
   Real dperB_N = 0.0;
   Real dperB_W = 0.0;
   Real dperB_E = 0.0;

   Real field() const {
      Real efield = apos * bpos * E_NE + apos * bneg * E_SE + aneg * bpos * E_NW + aneg * bneg * E_SW;
      efield /= ((apos + aneg) * (bpos + bneg) + EPS);
      if (Parameters::fieldSolverDiffusiveEterms) {
#ifdef FS_1ST_ORDER_SPACE
         // 1st order diffusive terms:
         efield -= bpos * bneg / (bpos + bneg + EPS) * (perB_S - perB_N);
         efield += apos * aneg / (apos + aneg + EPS) * (perB_W - perB_E);
#else
         // 2nd     order diffusive terms
         efield -= bpos * bneg / (bpos + bneg + EPS) * ((perB_S - HALF * dperB_S) - (perB_N + HALF * dperB_N));
         efield += apos * aneg / (apos + aneg + EPS) * ((perB_W - HALF * dperB_W) - (perB_E + HALF * dperB_E));
#endif
      }

      return efield;
   }
};

/*! \brief Low-level electric field propagation function.
 *
 * Computes the upwinded electric field X component along the cell's corresponding edge as the cross product of B and V
 * in the YZ plane. Also includes the calculation of the maximally allowed time step.
 *
 * Expects that the correct RHO and B fields are being passed, depending on the stage of the Runge-Kutta time stepping
 * method.
 *
 * Note that the background B field is excluded from the diffusive term calculations because they are equivalent to a
 * current term and the background field is curl-free.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param EGrid fsGrid holding the electric field
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param EGradPeGrid fsGrid holding the electron pressure gradient E field
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
void calculateEdgeElectricFieldX(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH>& EGradPeGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, cint i, cint j, cint k, cint& RKCase) {
   fsdebugCheck(technicalGrid, i, j, k);

   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.
   Real ay_pos, ay_neg; // Max. characteristic velocities to y-direction
   Real az_pos, az_neg; // Max. characteristic velocities to z-direction
   Real maxV = 0.0;     // Max velocity for CFL purposes
   Real c_y, c_z;       // Wave speeds to yz-directions

   // Get values at all four neighbours, result is written to SW.
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_SW = perBGrid.get(i, j, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_SE = perBGrid.get(i, j - 1, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_NE = perBGrid.get(i, j - 1, k - 1);
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_NW = perBGrid.get(i, j, k - 1);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_SW = BgBGrid.get(i, j, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_SE = BgBGrid.get(i, j - 1, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_NE = BgBGrid.get(i, j - 1, k - 1);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_NW = BgBGrid.get(i, j, k - 1);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_SW = momentsGrid.get(i, j, k);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_SE = momentsGrid.get(i, j - 1, k);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_NE = momentsGrid.get(i, j - 1, k - 1);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_NW = momentsGrid.get(i, j, k - 1);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_SW = dMomentsGrid.get(i, j, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_SE = dMomentsGrid.get(i, j - 1, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_NE = dMomentsGrid.get(i, j - 1, k - 1);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_NW = dMomentsGrid.get(i, j, k - 1);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_SW = dPerBGrid.get(i, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_SE = dPerBGrid.get(i, j - 1, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_NE = dPerBGrid.get(i, j - 1, k - 1);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_NW = dPerBGrid.get(i, j, k - 1);

   const auto rhomLimits = getRhomLimits({
       *moments_SW,
       *moments_SE,
       *moments_NW,
       *moments_NE,
   });

   const Real By_S = perb_SW->at(fsgrids::bfield::PERBY) + bgb_SW->at(fsgrids::bgbfield::BGBY);
   const Real Bz_W = perb_SW->at(fsgrids::bfield::PERBZ) + bgb_SW->at(fsgrids::bgbfield::BGBZ);
   const Real Bz_E = perb_SE->at(fsgrids::bfield::PERBZ) + bgb_SE->at(fsgrids::bgbfield::BGBZ);
   const Real By_N = perb_NW->at(fsgrids::bfield::PERBY) + bgb_NW->at(fsgrids::bgbfield::BGBY);
   const Real perBy_S = perb_SW->at(fsgrids::bfield::PERBY);
   const Real perBz_W = perb_SW->at(fsgrids::bfield::PERBZ);
   const Real perBz_E = perb_SE->at(fsgrids::bfield::PERBZ);
   const Real perBy_N = perb_NW->at(fsgrids::bfield::PERBY);
   Real Vy0 = moments_SW->at(fsgrids::moments::VY);
   Real Vz0 = moments_SW->at(fsgrids::moments::VZ);

   creal dBydx_S = dperb_SW->at(fsgrids::dperb::dPERBydx) + bgb_SW->at(fsgrids::bgbfield::dBGBydx);
   creal dBydz_S = dperb_SW->at(fsgrids::dperb::dPERBydz) + bgb_SW->at(fsgrids::bgbfield::dBGBydz);
   creal dBzdx_W = dperb_SW->at(fsgrids::dperb::dPERBzdx) + bgb_SW->at(fsgrids::bgbfield::dBGBzdx);
   creal dBzdy_W = dperb_SW->at(fsgrids::dperb::dPERBzdy) + bgb_SW->at(fsgrids::bgbfield::dBGBzdy);
   creal dBzdx_E = dperb_SE->at(fsgrids::dperb::dPERBzdx) + bgb_SE->at(fsgrids::bgbfield::dBGBzdx);
   creal dBzdy_E = dperb_SE->at(fsgrids::dperb::dPERBzdy) + bgb_SE->at(fsgrids::bgbfield::dBGBzdy);
   creal dBydx_N = dperb_NW->at(fsgrids::dperb::dPERBydx) + bgb_NW->at(fsgrids::bgbfield::dBGBydx);
   creal dBydz_N = dperb_NW->at(fsgrids::dperb::dPERBydz) + bgb_NW->at(fsgrids::bgbfield::dBGBydz);
   creal dperBydz_S = dperb_SW->at(fsgrids::dperb::dPERBydz);
   creal dperBydz_N = dperb_NW->at(fsgrids::dperb::dPERBydz);
   creal dperBzdy_W = dperb_SW->at(fsgrids::dperb::dPERBzdy);
   creal dperBzdy_E = dperb_SE->at(fsgrids::dperb::dPERBzdy);

   // Ex and characteristic speeds on this cell:
   // 1st order terms:
   Real Ex_SW = By_S * Vz0 - Bz_W * Vy0;

   const auto& gridSpacing = technicalGrid.getGridSpacing();
   // Resistive term
   if (Parameters::resistivity > 0) {
      Ex_SW += resistiveTerm(*bgb_SW, *perb_SW, *moments_SW, *dperb_SW,
                             {fsgrids::dperb::dPERBzdy, fsgrids::dperb::dPERBydz}, {gridSpacing[1], gridSpacing[2]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ex_SW += EHallGrid.get(i, j, k)->at(fsgrids::ehall::EXHALL_000_100);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ex_SW += EGradPeGrid.get(i, j, k)->at(fsgrids::egradpe::EXGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ex_SW += +HALF * ((By_S - HALF * dBydz_S) *
                         (-dmoments_SW->at(fsgrids::dmoments::dVzdy) - dmoments_SW->at(fsgrids::dmoments::dVzdz)) -
                     dBydz_S * Vz0 + SIXTH * dBydx_S * dmoments_SW->at(fsgrids::dmoments::dVzdx));
   Ex_SW += -HALF * ((Bz_W - HALF * dBzdy_W) *
                         (-dmoments_SW->at(fsgrids::dmoments::dVydy) - dmoments_SW->at(fsgrids::dmoments::dVydz)) -
                     dBzdy_W * Vy0 + SIXTH * dBzdx_W * dmoments_SW->at(fsgrids::dmoments::dVydx));
#endif
   auto wavespeeds = calculateWaveSpeedYZ(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                          technicalGrid.getGridSpacing(), i, j, k, i + 1, j, k, By_S, Bz_W, dBydx_S,
                                          dBydz_S, dBzdx_W, dBzdy_W, MINUS, MINUS, rhomLimits);
   c_y = wavespeeds.minVelocity();
   c_z = c_y;
   ay_neg = max(ZERO, -Vy0 + c_y);
   ay_pos = max(ZERO, +Vy0 + c_y);
   az_neg = max(ZERO, -Vz0 + c_z);
   az_pos = max(ZERO, +Vz0 + c_z);
   maxV = max(maxV, wavespeeds.cflSpeed(Vy0, Vz0));

   // Ex and characteristic speeds on j-1 neighbour:
   Vy0 = moments_SE->at(fsgrids::moments::VY);
   Vz0 = moments_SE->at(fsgrids::moments::VZ);

   // 1st order terms:
   Real Ex_SE = By_S * Vz0 - Bz_E * Vy0;

   // Resistive term
   if (Parameters::resistivity > 0) {
      Ex_SE += resistiveTerm(*bgb_SE, *perb_SE, *moments_SE, *dperb_SE,
                             {fsgrids::dperb::dPERBzdy, fsgrids::dperb::dPERBydz}, {gridSpacing[1], gridSpacing[2]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ex_SE += EHallGrid.get(i, j - 1, k)->at(fsgrids::ehall::EXHALL_010_110);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ex_SE += EGradPeGrid.get(i, j - 1, k)->at(fsgrids::egradpe::EXGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ex_SE += +HALF * ((By_S - HALF * dBydz_S) *
                         (+dmoments_SE->at(fsgrids::dmoments::dVzdy) - dmoments_SE->at(fsgrids::dmoments::dVzdz)) -
                     dBydz_S * Vz0 + SIXTH * dBydx_S * dmoments_SE->at(fsgrids::dmoments::dVzdx));
   Ex_SE += -HALF * ((Bz_E + HALF * dBzdy_E) *
                         (+dmoments_SE->at(fsgrids::dmoments::dVydy) - dmoments_SE->at(fsgrids::dmoments::dVydz)) +
                     dBzdy_E * Vy0 + SIXTH * dBzdx_E * dmoments_SE->at(fsgrids::dmoments::dVydx));
#endif

   wavespeeds = calculateWaveSpeedYZ(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                     technicalGrid.getGridSpacing(), i, j - 1, k, i + 1, j - 1, k, By_S, Bz_E, dBydx_S,
                                     dBydz_S, dBzdx_E, dBzdy_E, PLUS, MINUS, rhomLimits);
   c_y = wavespeeds.minVelocity();
   c_z = c_y;
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   maxV = max(maxV, wavespeeds.cflSpeed(Vy0, Vz0));

   // Ex and characteristic speeds on k-1 neighbour:
   Vy0 = moments_NW->at(fsgrids::moments::VY);
   Vz0 = moments_NW->at(fsgrids::moments::VZ);

   // 1st order terms:
   Real Ex_NW = By_N * Vz0 - Bz_W * Vy0;

   // Resistive term
   if (Parameters::resistivity > 0) {
      Ex_NW += resistiveTerm(*bgb_NW, *perb_NW, *moments_NW, *dperb_NW,
                             {fsgrids::dperb::dPERBzdy, fsgrids::dperb::dPERBydz}, {gridSpacing[1], gridSpacing[2]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ex_NW += EHallGrid.get(i, j, k - 1)->at(fsgrids::ehall::EXHALL_001_101);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ex_NW += EGradPeGrid.get(i, j, k - 1)->at(fsgrids::egradpe::EXGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ex_NW += +HALF * ((By_N + HALF * dBydz_N) *
                         (-dmoments_NW->at(fsgrids::dmoments::dVzdy) + dmoments_NW->at(fsgrids::dmoments::dVzdz)) +
                     dBydz_N * Vz0 + SIXTH * dBydx_N * dmoments_NW->at(fsgrids::dmoments::dVzdx));
   Ex_NW += -HALF * ((Bz_W - HALF * dBzdy_W) *
                         (-dmoments_NW->at(fsgrids::dmoments::dVydy) + dmoments_NW->at(fsgrids::dmoments::dVydz)) -
                     dBzdy_W * Vy0 + SIXTH * dBzdx_W * dmoments_NW->at(fsgrids::dmoments::dVydx));
#endif

   wavespeeds = calculateWaveSpeedYZ(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                     technicalGrid.getGridSpacing(), i, j, k - 1, i + 1, j, k - 1, By_N, Bz_W, dBydx_N,
                                     dBydz_N, dBzdx_W, dBzdy_W, MINUS, PLUS, rhomLimits);
   c_y = wavespeeds.minVelocity();
   c_z = c_y;
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   maxV = max(maxV, wavespeeds.cflSpeed(Vy0, Vz0));

   // Ex and characteristic speeds on j-1,k-1 neighbour:
   Vy0 = moments_NE->at(fsgrids::moments::VY);
   Vz0 = moments_NE->at(fsgrids::moments::VZ);

   // 1st order terms:
   Real Ex_NE = By_N * Vz0 - Bz_E * Vy0;

   // Resistive term
   if (Parameters::resistivity > 0) {
      Ex_NE += resistiveTerm(*bgb_NE, *perb_NE, *moments_NE, *dperb_NE,
                             {fsgrids::dperb::dPERBzdy, fsgrids::dperb::dPERBydz}, {gridSpacing[1], gridSpacing[2]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ex_NE += EHallGrid.get(i, j - 1, k - 1)->at(fsgrids::ehall::EXHALL_011_111);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ex_NE += EGradPeGrid.get(i, j - 1, k - 1)->at(fsgrids::egradpe::EXGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ex_NE += +HALF * ((By_N + HALF * dBydz_N) *
                         (+dmoments_NE->at(fsgrids::dmoments::dVzdy) + dmoments_NE->at(fsgrids::dmoments::dVzdz)) +
                     dBydz_N * Vz0 + SIXTH * dBydx_N * dmoments_NE->at(fsgrids::dmoments::dVzdx));
   Ex_NE += -HALF * ((Bz_E + HALF * dBzdy_E) *
                         (+dmoments_NE->at(fsgrids::dmoments::dVydy) + dmoments_NE->at(fsgrids::dmoments::dVydz)) +
                     dBzdy_E * Vy0 + SIXTH * dBzdx_E * dmoments_NE->at(fsgrids::dmoments::dVydx));
#endif

   wavespeeds = calculateWaveSpeedYZ(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                     technicalGrid.getGridSpacing(), i, j - 1, k - 1, i + 1, j - 1, k - 1, By_N, Bz_E,
                                     dBydx_N, dBydz_N, dBzdx_E, dBzdy_E, PLUS, PLUS, rhomLimits);
   c_y = wavespeeds.minVelocity();
   c_z = c_y;
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   maxV = max(maxV, wavespeeds.cflSpeed(Vy0, Vz0));

   // Calculate properly upwinded edge-averaged Ex:
   const FieldValues f(Ex_NE, Ex_SE, Ex_NW, Ex_SW, ay_pos, ay_neg, az_pos, az_neg, perBy_S, perBy_N, perBz_W, perBz_E,
                       dperBydz_S, dperBydz_N, dperBzdy_W, dperBzdy_E);
   EGrid.get(i, j, k)->at(fsgrids::efield::EX) = f.field();

   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      // compute maximum timestep for fieldsolver in this cell (CFL=1)
      Real min_dx = std::numeric_limits<Real>::max();
      min_dx = min(min_dx, gridSpacing[1]);
      min_dx = min(min_dx, gridSpacing[2]);
      // update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (maxV != ZERO)
         technicalGrid.get(i, j, k)->maxFsDt = min(technicalGrid.get(i, j, k)->maxFsDt, min_dx / maxV);
   }
}

/*! \brief Low-level electric field propagation function.
 *
 * Computes the upwinded electric field Y component along the cell's corresponding edge as the cross product of B and V
 * in the XZ plane. Also includes the calculation of the maximally allowed time step.
 *
 * Expects that the correct RHO and B fields are being passed, depending on the stage of the Runge-Kutta time stepping
 * method.
 *
 * Note that the background B field is excluded from the diffusive term calculations because they are equivalent to a
 * current term and the background field is curl-free.
 *
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
void calculateEdgeElectricFieldY(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH>& EGradPeGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, cint i, cint j, cint k, cint& RKCase) {
   fsdebugCheck(technicalGrid, i, j, k);

   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.
   Real ax_pos, ax_neg; // Max. characteristic velocities to x-direction
   Real az_pos, az_neg; // Max. characteristic velocities to z-direction
   Real maxV = 0.0;     // Max velocity for CFL purposes
   Real c_x, c_z;       // Wave speeds to xz-directions

   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_SW = perBGrid.get(i, j, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_SE = perBGrid.get(i, j, k - 1);
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_NW = perBGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_NE = perBGrid.get(i - 1, j, k - 1);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_SW = BgBGrid.get(i, j, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_SE = BgBGrid.get(i, j, k - 1);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_NW = BgBGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_NE = BgBGrid.get(i - 1, j, k - 1);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_SW = momentsGrid.get(i, j, k);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_SE = momentsGrid.get(i, j, k - 1);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_NW = momentsGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_NE = momentsGrid.get(i - 1, j, k - 1);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_SW = dMomentsGrid.get(i, j, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_SE = dMomentsGrid.get(i, j, k - 1);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_NW = dMomentsGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_NE = dMomentsGrid.get(i - 1, j, k - 1);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_SW = dPerBGrid.get(i, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_SE = dPerBGrid.get(i, j, k - 1);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_NW = dPerBGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_NE = dPerBGrid.get(i - 1, j, k - 1);

   const auto rhomLimits = getRhomLimits({
       *moments_SW,
       *moments_SE,
       *moments_NW,
       *moments_NE,
   });
   // Fetch required plasma parameters:
   const Real Bz_S = perb_SW->at(fsgrids::bfield::PERBZ) + bgb_SW->at(fsgrids::bgbfield::BGBZ);
   const Real Bx_W = perb_SW->at(fsgrids::bfield::PERBX) + bgb_SW->at(fsgrids::bgbfield::BGBX);
   const Real Bx_E = perb_SE->at(fsgrids::bfield::PERBX) + bgb_SE->at(fsgrids::bgbfield::BGBX);
   const Real Bz_N = perb_NW->at(fsgrids::bfield::PERBZ) + bgb_NW->at(fsgrids::bgbfield::BGBZ);
   const Real perBz_S = perb_SW->at(fsgrids::bfield::PERBZ);
   const Real perBx_W = perb_SW->at(fsgrids::bfield::PERBX);
   const Real perBx_E = perb_SE->at(fsgrids::bfield::PERBX);
   const Real perBz_N = perb_NW->at(fsgrids::bfield::PERBZ);
   Real Vx0 = moments_SW->at(fsgrids::moments::VX);
   Real Vz0 = moments_SW->at(fsgrids::moments::VZ);

   creal dBxdy_W = dperb_SW->at(fsgrids::dperb::dPERBxdy) + bgb_SW->at(fsgrids::bgbfield::dBGBxdy);
   creal dBxdz_W = dperb_SW->at(fsgrids::dperb::dPERBxdz) + bgb_SW->at(fsgrids::bgbfield::dBGBxdz);
   creal dBzdx_S = dperb_SW->at(fsgrids::dperb::dPERBzdx) + bgb_SW->at(fsgrids::bgbfield::dBGBzdx);
   creal dBzdy_S = dperb_SW->at(fsgrids::dperb::dPERBzdy) + bgb_SW->at(fsgrids::bgbfield::dBGBzdy);
   creal dBxdy_E = dperb_SE->at(fsgrids::dperb::dPERBxdy) + bgb_SE->at(fsgrids::bgbfield::dBGBxdy);
   creal dBxdz_E = dperb_SE->at(fsgrids::dperb::dPERBxdz) + bgb_SE->at(fsgrids::bgbfield::dBGBxdz);
   creal dBzdx_N = dperb_NW->at(fsgrids::dperb::dPERBzdx) + bgb_NW->at(fsgrids::bgbfield::dBGBzdx);
   creal dBzdy_N = dperb_NW->at(fsgrids::dperb::dPERBzdy) + bgb_NW->at(fsgrids::bgbfield::dBGBzdy);
   creal dperBzdx_S = dperb_SW->at(fsgrids::dperb::dPERBzdx);
   creal dperBzdx_N = dperb_NW->at(fsgrids::dperb::dPERBzdx);
   creal dperBxdz_W = dperb_SW->at(fsgrids::dperb::dPERBxdz);
   creal dperBxdz_E = dperb_SE->at(fsgrids::dperb::dPERBxdz);

   // Ey and characteristic speeds on this cell:
   // 1st order terms:
   Real Ey_SW = Bz_S * Vx0 - Bx_W * Vz0;

   const auto& gridSpacing = technicalGrid.getGridSpacing();
   // Resistive term
   if (Parameters::resistivity > 0) {
      Ey_SW += resistiveTerm(*bgb_SW, *perb_SW, *moments_SW, *dperb_SW,
                             {fsgrids::dperb::dPERBxdz, fsgrids::dperb::dPERBzdx}, {gridSpacing[2], gridSpacing[0]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ey_SW += EHallGrid.get(i, j, k)->at(fsgrids::ehall::EYHALL_000_010);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ey_SW += EGradPeGrid.get(i, j, k)->at(fsgrids::egradpe::EYGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms
   Ey_SW += +HALF * ((Bz_S - HALF * dBzdx_S) *
                         (-dmoments_SW->at(fsgrids::dmoments::dVxdx) - dmoments_SW->at(fsgrids::dmoments::dVxdz)) -
                     dBzdx_S * Vx0 + SIXTH * dBzdy_S * dmoments_SW->at(fsgrids::dmoments::dVxdy));
   Ey_SW += -HALF * ((Bx_W - HALF * dBxdz_W) *
                         (-dmoments_SW->at(fsgrids::dmoments::dVzdx) - dmoments_SW->at(fsgrids::dmoments::dVzdz)) -
                     dBxdz_W * Vz0 + SIXTH * dBxdy_W * dmoments_SW->at(fsgrids::dmoments::dVzdy));
#endif

   auto wavespeeds = calculateWaveSpeedXZ(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                          technicalGrid.getGridSpacing(), i, j, k, i, j + 1, k, Bx_W, Bz_S, dBxdy_W,
                                          dBxdz_W, dBzdx_S, dBzdy_S, MINUS, MINUS, rhomLimits);
   c_z = wavespeeds.minVelocity();
   c_x = c_z;
   az_neg = max(ZERO, -Vz0 + c_z);
   az_pos = max(ZERO, +Vz0 + c_z);
   ax_neg = max(ZERO, -Vx0 + c_x);
   ax_pos = max(ZERO, +Vx0 + c_x);
   maxV = max(maxV, wavespeeds.cflSpeed(Vz0, Vx0));

   // Ey and characteristic speeds on k-1 neighbour:
   Vx0 = moments_SE->at(fsgrids::moments::VX);
   Vz0 = moments_SE->at(fsgrids::moments::VZ);

   // 1st order terms:
   Real Ey_SE = Bz_S * Vx0 - Bx_E * Vz0;

   // Resistive term
   if (Parameters::resistivity > 0) {
      Ey_SE += resistiveTerm(*bgb_SE, *perb_SE, *moments_SE, *dperb_SE,
                             {fsgrids::dperb::dPERBxdz, fsgrids::dperb::dPERBzdx}, {gridSpacing[2], gridSpacing[0]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ey_SE += EHallGrid.get(i, j, k - 1)->at(fsgrids::ehall::EYHALL_001_011);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ey_SE += EGradPeGrid.get(i, j, k - 1)->at(fsgrids::egradpe::EYGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ey_SE += +HALF * ((Bz_S - HALF * dBzdx_S) *
                         (-dmoments_SE->at(fsgrids::dmoments::dVxdx) + dmoments_SE->at(fsgrids::dmoments::dVxdz)) -
                     dBzdx_S * Vx0 + SIXTH * dBzdy_S * dmoments_SE->at(fsgrids::dmoments::dVxdy));
   Ey_SE += -HALF * ((Bx_E + HALF * dBxdz_E) *
                         (-dmoments_SE->at(fsgrids::dmoments::dVzdx) + dmoments_SE->at(fsgrids::dmoments::dVzdz)) +
                     dBxdz_E * Vz0 + SIXTH * dBxdy_E * dmoments_SE->at(fsgrids::dmoments::dVzdy));
#endif

   wavespeeds = calculateWaveSpeedXZ(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                     technicalGrid.getGridSpacing(), i, j, k - 1, i, j + 1, k - 1, Bx_E, Bz_S, dBxdy_E,
                                     dBxdz_E, dBzdx_S, dBzdy_S, MINUS, PLUS, rhomLimits);
   c_z = wavespeeds.minVelocity();
   c_x = c_z;
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   maxV = max(maxV, wavespeeds.cflSpeed(Vz0, Vx0));

   // Ey and characteristic speeds on i-1 neighbour:
   Vz0 = moments_NW->at(fsgrids::moments::VZ);
   Vx0 = moments_NW->at(fsgrids::moments::VX);

   // 1st order terms:
   Real Ey_NW = Bz_N * Vx0 - Bx_W * Vz0;

   // Resistive term
   if (Parameters::resistivity > 0) {
      Ey_NW += resistiveTerm(*bgb_NW, *perb_NW, *moments_NW, *dperb_NW,
                             {fsgrids::dperb::dPERBxdz, fsgrids::dperb::dPERBzdx}, {gridSpacing[2], gridSpacing[0]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ey_NW += EHallGrid.get(i - 1, j, k)->at(fsgrids::ehall::EYHALL_100_110);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ey_NW += EGradPeGrid.get(i - 1, j, k)->at(fsgrids::egradpe::EYGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ey_NW += +HALF * ((Bz_N + HALF * dBzdx_N) *
                         (+dmoments_NW->at(fsgrids::dmoments::dVxdx) - dmoments_NW->at(fsgrids::dmoments::dVxdz)) +
                     dBzdx_N * Vx0 + SIXTH * dBzdy_N * dmoments_NW->at(fsgrids::dmoments::dVxdy));
   Ey_NW += -HALF * ((Bx_W - HALF * dBxdz_W) *
                         (+dmoments_NW->at(fsgrids::dmoments::dVzdx) - dmoments_NW->at(fsgrids::dmoments::dVzdz)) -
                     dBxdz_W * Vz0 + SIXTH * dBxdy_W * dmoments_NW->at(fsgrids::dmoments::dVzdy));
#endif

   wavespeeds = calculateWaveSpeedXZ(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                     technicalGrid.getGridSpacing(), i - 1, j, k, i - 1, j + 1, k, Bx_W, Bz_N, dBxdy_W,
                                     dBxdz_W, dBzdx_N, dBzdy_N, PLUS, MINUS, rhomLimits);
   c_z = wavespeeds.minVelocity();
   c_x = c_z;
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   maxV = max(maxV, wavespeeds.cflSpeed(Vz0, Vx0));

   // Ey and characteristic speeds on i-1,k-1 neighbour:
   Vz0 = moments_NE->at(fsgrids::moments::VZ);
   Vx0 = moments_NE->at(fsgrids::moments::VX);

   // 1st order terms:
   Real Ey_NE = Bz_N * Vx0 - Bx_E * Vz0;

   // Resistive term
   if (Parameters::resistivity > 0) {
      Ey_NE += resistiveTerm(*bgb_NE, *perb_NE, *moments_NE, *dperb_NE,
                             {fsgrids::dperb::dPERBxdz, fsgrids::dperb::dPERBzdx}, {gridSpacing[2], gridSpacing[0]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ey_NE += EHallGrid.get(i - 1, j, k - 1)->at(fsgrids::ehall::EYHALL_101_111);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ey_NE += EGradPeGrid.get(i - 1, j, k - 1)->at(fsgrids::egradpe::EYGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ey_NE += +HALF * ((Bz_N + HALF * dBzdx_N) *
                         (+dmoments_NE->at(fsgrids::dmoments::dVxdx) + dmoments_NE->at(fsgrids::dmoments::dVxdz)) +
                     dBzdx_N * Vx0 + SIXTH * dBzdy_N * dmoments_NE->at(fsgrids::dmoments::dVxdy));
   Ey_NE += -HALF * ((Bx_E + HALF * dBxdz_E) *
                         (+dmoments_NE->at(fsgrids::dmoments::dVzdx) + dmoments_NE->at(fsgrids::dmoments::dVzdz)) +
                     dBxdz_E * Vz0 + SIXTH * dBxdy_E * dmoments_NE->at(fsgrids::dmoments::dVzdy));
#endif

   wavespeeds = calculateWaveSpeedXZ(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                     technicalGrid.getGridSpacing(), i - 1, j, k - 1, i - 1, j + 1, k - 1, Bx_E, Bz_N,
                                     dBxdy_E, dBxdz_E, dBzdx_N, dBzdy_N, PLUS, PLUS, rhomLimits);
   c_z = wavespeeds.minVelocity();
   c_x = c_z;
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   maxV = max(maxV, wavespeeds.cflSpeed(Vz0, Vx0));

   // Calculate properly upwinded edge-averaged Ey:
   const FieldValues f(Ey_NE, Ey_SE, Ey_NW, Ey_SW, az_pos, az_neg, ax_pos, ax_neg, perBz_S, perBz_N, perBx_W, perBx_E,
                       dperBzdx_S, dperBzdx_N, dperBxdz_W, dperBxdz_E);
   EGrid.get(i, j, k)->at(fsgrids::efield::EY) = f.field();

   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      // compute maximum timestep for fieldsolver in this cell (CFL=1)
      Real min_dx = std::numeric_limits<Real>::max();
      ;
      min_dx = min(min_dx, gridSpacing[0]);
      min_dx = min(min_dx, gridSpacing[2]);
      // update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (maxV != ZERO)
         technicalGrid.get(i, j, k)->maxFsDt = min(technicalGrid.get(i, j, k)->maxFsDt, min_dx / maxV);
   }
}

/*! \brief Low-level electric field propagation function.
 *
 * Computes the upwinded electric field Z component along the cell's corresponding edge as the cross product of B and V
 * in the XY plane. Also includes the calculation of the maximally allowed time step.
 *
 * Expects that the correct RHO and B fields are being passed, depending on the stage of the Runge-Kutta time stepping
 * method.
 *
 * Note that the background B field is excluded from the diffusive term calculations because they are equivalent to a
 * current term and the background field is curl-free.
 *
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
void calculateEdgeElectricFieldZ(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH>& EGradPeGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, cint i, cint j, cint k, cint& RKCase) {
   fsdebugCheck(technicalGrid, i, j, k);

   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.
   Real ax_pos, ax_neg; // Max. characteristic velocities to x-direction
   Real ay_pos, ay_neg; // Max. characteristic velocities to y-direction
   Real maxV = 0.0;     // Max velocity for CFL purposes
   Real c_x, c_y;       // Characteristic speeds to xy-directions

   // Get read-only pointers to NE,NW,SE,SW states (SW is rw, result is written there):
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_SW = perBGrid.get(i, j, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_SE = perBGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_NE = perBGrid.get(i - 1, j - 1, k);
   std::array<Real, fsgrids::bfield::N_BFIELD>* perb_NW = perBGrid.get(i, j - 1, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_SW = BgBGrid.get(i, j, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_SE = BgBGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_NE = BgBGrid.get(i - 1, j - 1, k);
   std::array<Real, fsgrids::bgbfield::N_BGB>* bgb_NW = BgBGrid.get(i, j - 1, k);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_SW = momentsGrid.get(i, j, k);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_SE = momentsGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_NE = momentsGrid.get(i - 1, j - 1, k);
   std::array<Real, fsgrids::moments::N_MOMENTS>* moments_NW = momentsGrid.get(i, j - 1, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_SW = dMomentsGrid.get(i, j, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_SE = dMomentsGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_NE = dMomentsGrid.get(i - 1, j - 1, k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>* dmoments_NW = dMomentsGrid.get(i, j - 1, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_SW = dPerBGrid.get(i, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_SE = dPerBGrid.get(i - 1, j, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_NE = dPerBGrid.get(i - 1, j - 1, k);
   std::array<Real, fsgrids::dperb::N_DPERB>* dperb_NW = dPerBGrid.get(i, j - 1, k);

   const auto rhomLimits = getRhomLimits({
       *moments_SW,
       *moments_SE,
       *moments_NW,
       *moments_NE,
   });

   // Fetch needed plasma parameters/derivatives from the four cells:
   const Real Bx_S = perb_SW->at(fsgrids::bfield::PERBX) + bgb_SW->at(fsgrids::bgbfield::BGBX);
   const Real By_W = perb_SW->at(fsgrids::bfield::PERBY) + bgb_SW->at(fsgrids::bgbfield::BGBY);
   const Real By_E = perb_SE->at(fsgrids::bfield::PERBY) + bgb_SE->at(fsgrids::bgbfield::BGBY);
   const Real Bx_N = perb_NW->at(fsgrids::bfield::PERBX) + bgb_NW->at(fsgrids::bgbfield::BGBX);
   const Real perBx_S = perb_SW->at(fsgrids::bfield::PERBX);
   const Real perBy_W = perb_SW->at(fsgrids::bfield::PERBY);
   const Real perBy_E = perb_SE->at(fsgrids::bfield::PERBY);
   const Real perBx_N = perb_NW->at(fsgrids::bfield::PERBX);
   Real Vx0 = moments_SW->at(fsgrids::moments::VX);
   Real Vy0 = moments_SW->at(fsgrids::moments::VY);

   creal dBxdy_S = dperb_SW->at(fsgrids::dperb::dPERBxdy) + bgb_SW->at(fsgrids::bgbfield::dBGBxdy);
   creal dBxdz_S = dperb_SW->at(fsgrids::dperb::dPERBxdz) + bgb_SW->at(fsgrids::bgbfield::dBGBxdz);
   creal dBydx_W = dperb_SW->at(fsgrids::dperb::dPERBydx) + bgb_SW->at(fsgrids::bgbfield::dBGBydx);
   creal dBydz_W = dperb_SW->at(fsgrids::dperb::dPERBydz) + bgb_SW->at(fsgrids::bgbfield::dBGBydz);
   creal dBydx_E = dperb_SE->at(fsgrids::dperb::dPERBydx) + bgb_SE->at(fsgrids::bgbfield::dBGBydx);
   creal dBydz_E = dperb_SE->at(fsgrids::dperb::dPERBydz) + bgb_SE->at(fsgrids::bgbfield::dBGBydz);
   creal dBxdy_N = dperb_NW->at(fsgrids::dperb::dPERBxdy) + bgb_NW->at(fsgrids::bgbfield::dBGBxdy);
   creal dBxdz_N = dperb_NW->at(fsgrids::dperb::dPERBxdz) + bgb_NW->at(fsgrids::bgbfield::dBGBxdz);
   creal dperBxdy_S = dperb_SW->at(fsgrids::dperb::dPERBxdy);
   creal dperBxdy_N = dperb_NW->at(fsgrids::dperb::dPERBxdy);
   creal dperBydx_W = dperb_SW->at(fsgrids::dperb::dPERBydx);
   creal dperBydx_E = dperb_SE->at(fsgrids::dperb::dPERBydx);

   // Ez and characteristic speeds on SW cell:
   // 1st order terms:
   Real Ez_SW = Bx_S * Vy0 - By_W * Vx0;

   const auto& gridSpacing = technicalGrid.getGridSpacing();
   // Resistive term
   if (Parameters::resistivity > 0) {
      Ez_SW += resistiveTerm(*bgb_SW, *perb_SW, *moments_SW, *dperb_SW,
                             {fsgrids::dperb::dPERBydx, fsgrids::dperb::dPERBxdy}, {gridSpacing[0], gridSpacing[1]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ez_SW += EHallGrid.get(i, j, k)->at(fsgrids::ehall::EZHALL_000_001);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ez_SW += EGradPeGrid.get(i, j, k)->at(fsgrids::egradpe::EZGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ez_SW += +HALF * ((Bx_S - HALF * dBxdy_S) *
                         (-dmoments_SW->at(fsgrids::dmoments::dVydx) - dmoments_SW->at(fsgrids::dmoments::dVydy)) -
                     dBxdy_S * Vy0 + SIXTH * dBxdz_S * dmoments_SW->at(fsgrids::dmoments::dVydz));
   Ez_SW += -HALF * ((By_W - HALF * dBydx_W) *
                         (-dmoments_SW->at(fsgrids::dmoments::dVxdx) - dmoments_SW->at(fsgrids::dmoments::dVxdy)) -
                     dBydx_W * Vx0 + SIXTH * dBydz_W * dmoments_SW->at(fsgrids::dmoments::dVxdz));
#endif

   // Calculate maximum wave speed (fast magnetosonic speed) on SW cell. In order
   // to get Alfven speed we need to calculate some reconstruction coeff. for Bz:
   auto wavespeeds = calculateWaveSpeedXY(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                          technicalGrid.getGridSpacing(), i, j, k, i, j, k + 1, Bx_S, By_W, dBxdy_S,
                                          dBxdz_S, dBydx_W, dBydz_W, MINUS, MINUS, rhomLimits);
   c_x = wavespeeds.minVelocity();
   c_y = c_x;
   ax_neg = max(ZERO, -Vx0 + c_x);
   ax_pos = max(ZERO, +Vx0 + c_x);
   ay_neg = max(ZERO, -Vy0 + c_y);
   ay_pos = max(ZERO, +Vy0 + c_y);
   maxV = max(maxV, wavespeeds.cflSpeed(Vx0, Vy0));

   // Ez and characteristic speeds on SE (i-1) cell:
   Vx0 = moments_SE->at(fsgrids::moments::VX);
   Vy0 = moments_SE->at(fsgrids::moments::VY);

   // 1st order terms:
   Real Ez_SE = Bx_S * Vy0 - By_E * Vx0;

   // Resistive term
   if (Parameters::resistivity > 0) {
      Ez_SE += resistiveTerm(*bgb_SE, *perb_SE, *moments_SE, *dperb_SE,
                             {fsgrids::dperb::dPERBydx, fsgrids::dperb::dPERBxdy}, {gridSpacing[0], gridSpacing[1]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ez_SE += EHallGrid.get(i - 1, j, k)->at(fsgrids::ehall::EZHALL_100_101);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ez_SE += EGradPeGrid.get(i - 1, j, k)->at(fsgrids::egradpe::EZGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ez_SE += +HALF * ((Bx_S - HALF * dBxdy_S) *
                         (+dmoments_SE->at(fsgrids::dmoments::dVydx) - dmoments_SE->at(fsgrids::dmoments::dVydy)) -
                     dBxdy_S * Vy0 + SIXTH * dBxdz_S * dmoments_SE->at(fsgrids::dmoments::dVydz));
   Ez_SE += -HALF * ((By_E + HALF * dBydx_E) *
                         (+dmoments_SE->at(fsgrids::dmoments::dVxdx) - dmoments_SE->at(fsgrids::dmoments::dVxdy)) +
                     dBydx_E * Vx0 + SIXTH * dBydz_E * dmoments_SE->at(fsgrids::dmoments::dVxdz));
#endif

   wavespeeds = calculateWaveSpeedXY(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                     technicalGrid.getGridSpacing(), i - 1, j, k, i - 1, j, k + 1, Bx_S, By_E, dBxdy_S,
                                     dBxdz_S, dBydx_E, dBydz_E, PLUS, MINUS, rhomLimits);
   c_x = wavespeeds.minVelocity();
   c_y = c_x;
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   maxV = max(maxV, wavespeeds.cflSpeed(Vx0, Vy0));

   // Ez and characteristic speeds on NW (j-1) cell:
   Vx0 = moments_NW->at(fsgrids::moments::VX);
   Vy0 = moments_NW->at(fsgrids::moments::VY);

   // 1st order terms:
   Real Ez_NW = Bx_N * Vy0 - By_W * Vx0;

   // Resistive term
   if (Parameters::resistivity > 0) {
      Ez_NW += resistiveTerm(*bgb_NW, *perb_NW, *moments_NW, *dperb_NW,
                             {fsgrids::dperb::dPERBydx, fsgrids::dperb::dPERBxdy}, {gridSpacing[0], gridSpacing[1]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ez_NW += EHallGrid.get(i, j - 1, k)->at(fsgrids::ehall::EZHALL_010_011);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ez_NW += EGradPeGrid.get(i, j - 1, k)->at(fsgrids::egradpe::EZGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ez_NW += +HALF * ((Bx_N + HALF * dBxdy_N) *
                         (-dmoments_NW->at(fsgrids::dmoments::dVydx) + dmoments_NW->at(fsgrids::dmoments::dVydy)) +
                     dBxdy_N * Vy0 + SIXTH * dBxdz_N * dmoments_NW->at(fsgrids::dmoments::dVydz));
   Ez_NW += -HALF * ((By_W - HALF * dBydx_W) *
                         (-dmoments_NW->at(fsgrids::dmoments::dVxdx) + dmoments_NW->at(fsgrids::dmoments::dVxdy)) -
                     dBydx_W * Vx0 + SIXTH * dBydz_W * dmoments_NW->at(fsgrids::dmoments::dVxdz));
#endif

   wavespeeds = calculateWaveSpeedXY(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                     technicalGrid.getGridSpacing(), i, j - 1, k, i, j - 1, k + 1, Bx_N, By_W, dBxdy_N,
                                     dBxdz_N, dBydx_W, dBydz_W, MINUS, PLUS, rhomLimits);
   c_x = wavespeeds.minVelocity();
   c_y = c_x;
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   maxV = max(maxV, wavespeeds.cflSpeed(Vx0, Vy0));

   // Ez and characteristic speeds on NE (i-1,j-1) cell:
   Vx0 = moments_NE->at(fsgrids::moments::VX);
   Vy0 = moments_NE->at(fsgrids::moments::VY);

   // 1st order terms:
   Real Ez_NE = Bx_N * Vy0 - By_E * Vx0;

   // Resistive term
   if (Parameters::resistivity > 0) {
      Ez_NE += resistiveTerm(*bgb_NE, *perb_NE, *moments_NE, *dperb_NE,
                             {fsgrids::dperb::dPERBydx, fsgrids::dperb::dPERBxdy}, {gridSpacing[0], gridSpacing[1]});
   }

   // Hall term
   if (Parameters::ohmHallTerm > 0) {
      Ez_NE += EHallGrid.get(i - 1, j - 1, k)->at(fsgrids::ehall::EZHALL_110_111);
   }

   // Electron pressure gradient term
   if (Parameters::ohmGradPeTerm > 0) {
      Ez_NE += EGradPeGrid.get(i - 1, j - 1, k)->at(fsgrids::egradpe::EZGRADPE);
   }

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ez_NE += +HALF * ((Bx_N + HALF * dBxdy_N) *
                         (+dmoments_NE->at(fsgrids::dmoments::dVydx) + dmoments_NE->at(fsgrids::dmoments::dVydy)) +
                     dBxdy_N * Vy0 + SIXTH * dBxdz_N * dmoments_NE->at(fsgrids::dmoments::dVydz));
   Ez_NE += -HALF * ((By_E + HALF * dBydx_E) *
                         (+dmoments_NE->at(fsgrids::dmoments::dVxdx) + dmoments_NE->at(fsgrids::dmoments::dVxdy)) +
                     dBydx_E * Vx0 + SIXTH * dBydz_E * dmoments_NE->at(fsgrids::dmoments::dVxdz));
#endif

   wavespeeds = calculateWaveSpeedXY(perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid,
                                     technicalGrid.getGridSpacing(), i - 1, j - 1, k, i - 1, j - 1, k + 1, Bx_N, By_E,
                                     dBxdy_N, dBxdz_N, dBydx_E, dBydz_E, PLUS, PLUS, rhomLimits);
   c_x = wavespeeds.minVelocity();
   c_y = c_x;
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   maxV = max(maxV, wavespeeds.cflSpeed(Vx0, Vy0));

   // Calculate properly upwinded edge-averaged Ez:
   const FieldValues f(Ez_NE, Ez_SE, Ez_NW, Ez_SW, ax_pos, ax_neg, ay_pos, ay_neg, perBx_S, perBx_N, perBy_W, perBy_E,
                       dperBxdy_S, dperBxdy_N, dperBydx_W, dperBydx_E);
   EGrid.get(i, j, k)->at(fsgrids::efield::EZ) = f.field();

   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      // compute maximum timestep for fieldsolver in this cell (CFL=1)
      Real min_dx = std::numeric_limits<Real>::max();
      min_dx = min(min_dx, gridSpacing[0]);
      min_dx = min(min_dx, gridSpacing[1]);
      // update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (maxV != ZERO)
         technicalGrid.get(i, j, k)->maxFsDt = min(technicalGrid.get(i, j, k)->maxFsDt, min_dx / maxV);
   }
}

/*! \brief Electric field propagation function.
 *
 * Calls the general or the system boundary electric field propagation functions.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param EGrid fsGrid holding the electric field
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param EGradPeGrid fsGrid holding the electron pressure gradient E field
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 *
 * \sa calculateUpwindedElectricFieldSimple calculateEdgeElectricFieldX calculateEdgeElectricFieldY
 * calculateEdgeElectricFieldZ
 *
 */
void calculateElectricField(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH>& EGradPeGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, cint i, cint j, cint k,
    SysBoundary& sysBoundaries, cint& RKCase) {
   cuint cellSysBoundaryFlag = technicalGrid.get(i, j, k)->sysBoundaryFlag;

   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
       cellSysBoundaryFlag == sysboundarytype::OUTER_BOUNDARY_PADDING)
      return;

   cuint bitfield = technicalGrid.get(i, j, k)->SOLVE;

   if ((bitfield & compute::EX) == compute::EX) {
      calculateEdgeElectricFieldX(perBGrid, EGrid, EHallGrid, EGradPeGrid, momentsGrid, dPerBGrid, dMomentsGrid,
                                  BgBGrid, technicalGrid, i, j, k, RKCase);
   } else {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, 0);
   }

   if ((bitfield & compute::EY) == compute::EY) {
      calculateEdgeElectricFieldY(perBGrid, EGrid, EHallGrid, EGradPeGrid, momentsGrid, dPerBGrid, dMomentsGrid,
                                  BgBGrid, technicalGrid, i, j, k, RKCase);
   } else {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, 1);
   }

   if ((bitfield & compute::EZ) == compute::EZ) {
      calculateEdgeElectricFieldZ(perBGrid, EGrid, EHallGrid, EGradPeGrid, momentsGrid, dPerBGrid, dMomentsGrid,
                                  BgBGrid, technicalGrid, i, j, k, RKCase);
   } else {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, 2);
   }
}

/*! \brief High-level electric field computation function.
 *
 * Transfers the derivatives, calculates the edge electric fields and transfers the new electric fields.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities at runge-kutta t=0
 * \param perBDt2Grid fsGrid holding the perturbed B quantities at runge-kutta t=0.5
 * \param EGrid fsGrid holding the Electric field quantities at runge-kutta t=0
 * \param EDt2Grid fsGrid holding the Electric field quantities at runge-kutta t=0.5
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param EGradPeGrid fsGrid holding the electron pressure gradient E field
 * \param momentsGrid fsGrid holding the moment quantities at runge-kutta t=0
 * \param momentsDt2Grid fsGrid holding the moment quantities at runge-kutta t=0.5
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derivatives of moments
 * \param dMomentsDt2Grid fsGrid holding the derivatives of moments at runge-kutta t=0.5
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 *
 * \sa calculateElectricField calculateEdgeElectricFieldX calculateEdgeElectricFieldY calculateEdgeElectricFieldZ
 */
void calculateUpwindedElectricFieldSimple(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH>& EGradPeGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH>& EGradPeDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, SysBoundary& sysBoundaries, cint& RKCase,
    const bool communicateEGradPeOrMomentsDerivatives) {
   // const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const fsgrid::FsIndex_t* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0] * gridDims[1] * gridDims[2];
   phiprof::Timer upwindedETimer{"Calculate upwinded electric field"};
   int computeTimerID{phiprof::initializeTimer("Electric field compute cells")};

   phiprof::Timer mpiTimer{"Electric field ghost updates MPI", {"MPI"}};
   // Update ghosts if necessary, unless previous terms have already updated them
   if (P::ohmHallTerm > 0) {
      EHallGrid.updateGhostCells();
   }
   if (P::ohmGradPeTerm > 0 && communicateEGradPeOrMomentsDerivatives) {
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         EGradPeGrid.updateGhostCells();
      } else {
         EGradPeDt2Grid.updateGhostCells();
      }
   }
   if (P::ohmHallTerm == 0) {
      dPerBGrid.updateGhostCells();
   }
   if (P::ohmHallTerm == 0 && P::ohmGradPeTerm == 0 && communicateEGradPeOrMomentsDerivatives) {
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         dMomentsGrid.updateGhostCells();
      } else {
         dMomentsDt2Grid.updateGhostCells();
      }
   }

   mpiTimer.stop();

// Calculate upwinded electric field on inner cells
#pragma omp parallel
   {
      phiprof::Timer computeTimer{computeTimerID};
#pragma omp for collapse(2)
      for (fsgrid::FsIndex_t k = 0; k < gridDims[2]; k++) {
         for (fsgrid::FsIndex_t j = 0; j < gridDims[1]; j++) {
            for (fsgrid::FsIndex_t i = 0; i < gridDims[0]; i++) {
               if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
                  calculateElectricField(perBGrid, EGrid, EHallGrid, EGradPeGrid, momentsGrid, dPerBGrid, dMomentsGrid,
                                         BgBGrid, technicalGrid, i, j, k, sysBoundaries, RKCase);
               } else { // RKCase == RK_ORDER2_STEP1
                  calculateElectricField(perBDt2Grid, EDt2Grid, EHallGrid, EGradPeDt2Grid, momentsDt2Grid, dPerBGrid,
                                         dMomentsDt2Grid, BgBGrid, technicalGrid, i, j, k, sysBoundaries, RKCase);
               }
            }
         }
      }
      computeTimer.stop(N_cells, "Spatial Cells");
   }

   mpiTimer.start();
   // Exchange electric field with neighbouring processes
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      EGrid.updateGhostCells();
   } else {
      EDt2Grid.updateGhostCells();
   }
   mpiTimer.stop();

   upwindedETimer.stop(N_cells, "Spatial Cells");
}
