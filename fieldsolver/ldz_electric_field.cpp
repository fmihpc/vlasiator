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
private:
   Real alfvenSq = 0.0;
   Real soundSq = 0.0;
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
public:
   Wavespeeds(Real bmag2, Real rhom, Real p11, Real p22, Real p33, const std::array<Real, 3>& gridSpacing)
       : alfvenSq(divideIfNonZero(bmag2, pc::MU_0 * rhom)), soundSq(divideIfNonZero(p11 + p22 + p33, 2.0 * rhom)),
         whistler(Parameters::ohmHallTerm > 0
                      ? sqrt(alfvenSq) *
                            (1 + divideIfNonZero(2 * M_PI * M_PI * pc::MASS_PROTON * pc::MASS_PROTON,
                                                 gridSpacing[0] * gridSpacing[0] * rhom * pc::CHARGE * pc::CHARGE *
                                                     pc::MU_0) /
                                     sqrt(1 + divideIfNonZero(M_PI * M_PI * pc::MASS_PROTON * pc::MASS_PROTON,
                                                              gridSpacing[0] * gridSpacing[0] * rhom * pc::CHARGE *
                                                                  pc::CHARGE * pc::MU_0)))
                      : 0.0) {}

   Real minVelocity() const { return min(Parameters::maxWaveVelocity, sqrt(alfvenSq + soundSq) + whistler); }

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
      const Real vMS = sqrt(alfvenSq + soundSq);
      return max(v + vMS, v + whistler);
   }
};

struct Limits {
   Real min = 0.0;
   Real max = 0.0;
};

/*!< Provide linearly-interpolated moments (limited to min/max of neighborign values for rhom)
 *   as well as Balsara-reconstructed magnetic field components used in the determination of
 *   characteristic wave speeds, used in turn in the upwind constrained transport Londrillo &
 *   del Zanna algorithm.
 *   
 *   Balsara DS (2009) Divergence-free reconstruction of magnetic fields and WENO schemes for
 *   magnetohydrodynamics. J Comput Phys 228:5040–5056. https://doi.org/10.1016/j.jcp.2009.03.038
 *   
 *   Londrillo P, Del Zanna L (2004) On the divergence-free condition in Godunov-type schemes for
 *   ideal magnetohydrodynamics: the upwind constrained transport method. J Comput Phys 195:17–48.
 *   https://doi.org/10.1016/j.jcp.2003.09.016
 */
struct Reconstructions {
private:
   const std::array<Real, fsgrids::bfield::N_BFIELD>& perb;
   const std::array<Real, fsgrids::bfield::N_BFIELD>& nbr_perb;
   const std::array<Real, fsgrids::dperb::N_DPERB>& dperb;
   const std::array<Real, fsgrids::dperb::N_DPERB>& nbr_dperb;
   const std::array<Real, fsgrids::bgbfield::N_BGB>& bgb;
   const std::array<Real, fsgrids::bgbfield::N_BGB>& nbr_bgb;
   const std::array<Real, fsgrids::moments::N_MOMENTS>& moment;
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS>& dmoment;
   const Limits& rhomLimits;

   template <size_t N>
   std::tuple<Real, Real> compute(const std::array<Real, N>& nbr_arr, const std::array<Real, N>& arr, size_t i,
                                  size_t j, Real mul0, Real mul1) const {
      const Real a = nbr_arr[i] + nbr_bgb[j];
      const Real b = arr[i] + bgb[j];

      return {mul0 * (a + b), mul1 * (a - b)};
   }

   /*!<
    * Interpolate the volume-averaged, i.e. cell-centred moment value by HALF a cell to the given cell edge.
    * \param dir0 +/-1 in the first direction
    * \param dir1 +/-1 in the second direction
    * \param i index into moment array
    * \param j index into dmoment array for derivative of moment i in dir0
    * \param k index into dmoment array for derivative of moment i in dir1
    * \param min limit inteprolated result to no smaller than this
    * \param max limit interpolated result to no larger than this
    */
   Real interpolateAndLimitMoments(Real dir0, Real dir1, size_t i, size_t j, size_t k, Real min, Real max) const {
      return std::clamp(moment[i] + HALF * (dir0 * dmoment[j] + dir1 * dmoment[k]), min, max);
   }

public:
   Reconstructions(fsgrids::perbspan perB,
                     fsgrids::constdperbspan dPerB,
                     fsgrids::constbgbspan BgB,
                     fsgrids::constmomentsspan moments,
                     fsgrids::constdmomentsspan dMoments, size_t self, size_t nbr,
                     const Limits& rhomLimits)
       : perb(perB[self]), nbr_perb(perB[nbr]), dperb(dPerB[self]), nbr_dperb(dPerB[nbr]), bgb(BgB[self]),
         nbr_bgb(BgB[nbr]), moment(moments[self]), dmoment(dMoments[self]), rhomLimits(rhomLimits) {}

   std::tuple<Real, Real> perBCoeffs(size_t i, size_t j) const { return compute(nbr_perb, perb, i, j, HALF, 1.0); }

   std::tuple<Real, Real> dPerBCoeffs(size_t i, size_t j) const { return compute(nbr_dperb, dperb, i, j, 1.0, 1.0); }

   Real rhom(Real dir0, Real dir1, size_t i, size_t j, size_t k) const {
      return interpolateAndLimitMoments(dir0, dir1, i, j, k, rhomLimits.min, rhomLimits.max);
   }

   Real p(Real dir0, Real dir1, size_t i, size_t j, size_t k) const {
      return interpolateAndLimitMoments(dir0, dir1, i, j, k, 0.0, std::numeric_limits<Real>::max());
   }

   static Real squared(Real a, Real b) { return a * a + TWELWTH * b * b; }
};

/*! \brief Low-level helper function.
 *
 * Computes the magnetosonic speed in the YZ plane. Used in upwinding the electric field X component,
 * at the interface between cells self and nbr.
 *
 * Expects that the correct RHO and B fields are being passed, depending on the
 * stage of the Runge-Kutta time stepping method.
 *
 * If fields are not propagated, returns 0.0 as there is no information propagating.
 *
 * \param perB fsGrid holding the perturbed B quantities
 * \param moments fsGrid holding the moment quantities
 * \param dPerB fsGrid holding the derivatives of perturbed B
 * \param dMoments fsGrid holding the derviatives of moments
 * \param BgB fsGrid holding the background B quantities
 * \param gridSpacing grid cell size in x,y,z
 * \param rhomLimits allowed min and max density when interpolating
 * \param self Current cell index
 * \param nbr Neighbor cell index
 * \param By Current cell's By
 * \param Bz Current cell's Bz
 * \param dBydx dBydx derivative
 * \param dBydz dBydz derivative
 * \param dBzdx dBzdx derivative
 * \param dBzdy dBzdy derivative
 * \param ydir +1 or -1 depending on the interpolation direction in y
 * \param zdir +1 or -1 depending on the interpolation direction in z
 */
Wavespeeds calculateWaveSpeedYZ(fsgrids::perbspan perB,
                                fsgrids::constmomentsspan moments,
                                fsgrids::constdperbspan dPerB,
                                fsgrids::constdmomentsspan dMoments,
                                fsgrids::constbgbspan BgB,
                                const std::array<Real, 3>& gridSpacing, const Limits& rhomLimits, size_t self,
                                size_t nbr, Real By, Real Bz, Real dBydx, Real dBydz, Real dBzdx, Real dBzdy, Real ydir,
                                Real zdir) {
   const Reconstructions rec(perB, dPerB, BgB, moments, dMoments, self, nbr, rhomLimits);

   const Real rhom =
       rec.rhom(ydir, zdir, fsgrids::moments::RHOM, fsgrids::dmoments::drhomdy, fsgrids::dmoments::drhomdz);
   const Real p11 = rec.p(ydir, zdir, fsgrids::moments::P_11, fsgrids::dmoments::dp11dy, fsgrids::dmoments::dp11dz);
   const Real p22 = rec.p(ydir, zdir, fsgrids::moments::P_22, fsgrids::dmoments::dp22dy, fsgrids::dmoments::dp22dz);
   const Real p33 = rec.p(ydir, zdir, fsgrids::moments::P_33, fsgrids::dmoments::dp33dy, fsgrids::dmoments::dp33dz);

   const auto [A_0, A_X] = rec.perBCoeffs(fsgrids::bfield::PERBX, fsgrids::bgbfield::BGBX);
   const auto [A_Y, A_XY] = rec.dPerBCoeffs(fsgrids::dperb::dPERBxdy, fsgrids::bgbfield::dBGBxdy);
   const auto [A_Z, A_XZ] = rec.dPerBCoeffs(fsgrids::dperb::dPERBxdz, fsgrids::bgbfield::dBGBxdz);

   const Real bx2 =
       Reconstructions::squared(A_0 + HALF * (ydir * A_Y + zdir * A_Z), A_X + HALF * (ydir * A_XY + zdir * A_XZ));
   const Real by2 = Reconstructions::squared(By + zdir * HALF * dBydz, dBydx);
   const Real bz2 = Reconstructions::squared(Bz + ydir * HALF * dBzdy, dBzdx);

   return Wavespeeds(bx2 + by2 + bz2, rhom, p11, p22, p33, gridSpacing);
}

/*! \brief Low-level helper function.
 *
 * Computes the magnetosonic speed in the XZ plane. Used in upwinding the electric field Y component,
 * at the interface between cells self and nbr.
 *
 * Expects that the correct RHO and B fields are being passed, depending on the stage of the Runge-Kutta time stepping
 * method.
 *
 * If fields are not propagated, returns 0.0 as there is no information propagating.
 *
 * \param perB fsGrid holding the perturbed B quantities
 * \param moments fsGrid holding the moment quantities
 * \param dPerB fsGrid holding the derivatives of perturbed B
 * \param dMoments fsGrid holding the derviatives of moments
 * \param BgB fsGrid holding the background B quantities
 * \param gridSpacing grid cell size in x,y,z
 * \param rhomLimits allowed min and max density when interpolating
 * \param self Current cell index
 * \param nbr Neighbor cell index
 * \param Bx Current cell's Bx
 * \param Bz Current cell's Bz
 * \param dBxdy dBxdy derivative
 * \param dBxdz dBxdz derivative
 * \param dBzdx dBzdx derivative
 * \param dBzdy dBzdy derivative
 * \param xdir +1 or -1 depending on the interpolation direction in x
 * \param zdir +1 or -1 depending on the interpolation direction in z
 */
Wavespeeds calculateWaveSpeedXZ(fsgrids::perbspan perB,
                                fsgrids::constmomentsspan moments,
                                fsgrids::constdperbspan dPerB,
                                fsgrids::constdmomentsspan dMoments,
                                fsgrids::constbgbspan BgB,
                                const std::array<Real, 3>& gridSpacing, const Limits& rhomLimits, size_t self,
                                size_t nbr, Real Bx, Real Bz, Real dBxdy, Real dBxdz, Real dBzdx, Real dBzdy, Real xdir,
                                Real zdir) {
   const Reconstructions rec(perB, dPerB, BgB, moments, dMoments, self, nbr, rhomLimits);

   const Real rhom =
       rec.rhom(xdir, zdir, fsgrids::moments::RHOM, fsgrids::dmoments::drhomdx, fsgrids::dmoments::drhomdz);
   const Real p11 = rec.p(xdir, zdir, fsgrids::moments::P_11, fsgrids::dmoments::dp11dx, fsgrids::dmoments::dp11dz);
   const Real p22 = rec.p(xdir, zdir, fsgrids::moments::P_22, fsgrids::dmoments::dp22dx, fsgrids::dmoments::dp22dz);
   const Real p33 = rec.p(xdir, zdir, fsgrids::moments::P_33, fsgrids::dmoments::dp33dx, fsgrids::dmoments::dp33dz);

   const auto [B_0, B_Y] = rec.perBCoeffs(fsgrids::bfield::PERBY, fsgrids::bgbfield::BGBY);
   const auto [B_X, B_XY] = rec.dPerBCoeffs(fsgrids::dperb::dPERBydx, fsgrids::bgbfield::dBGBydx);
   const auto [B_Z, B_YZ] = rec.dPerBCoeffs(fsgrids::dperb::dPERBydz, fsgrids::bgbfield::dBGBydz);

   const Real by2 =
       Reconstructions::squared(B_0 + HALF * (xdir * B_X + zdir * B_Z), B_Y + HALF * (xdir * B_XY + zdir * B_YZ));
   const Real bx2 = Reconstructions::squared(Bx + zdir * HALF * dBxdz, dBxdy);
   const Real bz2 = Reconstructions::squared(Bz + xdir * HALF * dBzdx, dBzdy);

   return Wavespeeds(bx2 + by2 + bz2, rhom, p11, p22, p33, gridSpacing);
}

/*! \brief Low-level helper function.
 *
 * Computes the magnetosonic speed in the XY plane. Used in upwinding the electric field Z component,
 * at the interface between cells self and nbr.
 *
 * Expects that the correct RHO and B fields are being passed, depending on the stage of the Runge-Kutta time stepping
 * method.
 *
 * If fields are not propagated, returns 0.0 as there is no information propagating.
 *
 * \param perB fsGrid holding the perturbed B quantities
 * \param moments fsGrid holding the moment quantities
 * \param dPerB fsGrid holding the derivatives of perturbed B
 * \param dMoments fsGrid holding the derviatives of moments
 * \param BgB fsGrid holding the background B quantities
 * \param gridSpacing grid cell size in x,y,z
 * \param rhomLimits allowed min and max density when interpolating
 * \param self Current cell index
 * \param nbr Neighbor cell index
 * \param Bx Current cell's Bx
 * \param By Current cell's By
 * \param dBxdy dBxdy derivative
 * \param dBxdz dBxdz derivative
 * \param dBydx dBydx derivative
 * \param dBydz dBydz derivative
 * \param xdir +1 or -1 depending on the interpolation direction in x
 * \param ydir +1 or -1 depending on the interpolation direction in y
 */
Wavespeeds calculateWaveSpeedXY(fsgrids::perbspan perB,
                                fsgrids::constmomentsspan moments,
                                fsgrids::constdperbspan dPerB,
                                fsgrids::constdmomentsspan dMoments,
                                fsgrids::constbgbspan BgB,
                                const std::array<Real, 3>& gridSpacing, const Limits& rhomLimits, size_t self,
                                size_t nbr, Real Bx, Real By, Real dBxdy, Real dBxdz, Real dBydx, Real dBydz, Real xdir,
                                Real ydir) {
   const Reconstructions rec(perB, dPerB, BgB, moments, dMoments, self, nbr, rhomLimits);

   const Real rhom =
       rec.rhom(xdir, ydir, fsgrids::moments::RHOM, fsgrids::dmoments::drhomdx, fsgrids::dmoments::drhomdy);
   const Real p11 = rec.p(xdir, ydir, fsgrids::moments::P_11, fsgrids::dmoments::dp11dx, fsgrids::dmoments::dp11dy);
   const Real p22 = rec.p(xdir, ydir, fsgrids::moments::P_22, fsgrids::dmoments::dp22dx, fsgrids::dmoments::dp22dy);
   const Real p33 = rec.p(xdir, ydir, fsgrids::moments::P_33, fsgrids::dmoments::dp33dx, fsgrids::dmoments::dp33dy);

   const auto [C_0, C_Z] = rec.perBCoeffs(fsgrids::bfield::PERBZ, fsgrids::bgbfield::BGBZ);
   const auto [C_X, C_XZ] = rec.dPerBCoeffs(fsgrids::dperb::dPERBzdx, fsgrids::bgbfield::dBGBzdx);
   const auto [C_Y, C_YZ] = rec.dPerBCoeffs(fsgrids::dperb::dPERBzdy, fsgrids::bgbfield::dBGBzdy);

   const Real bz2 =
       Reconstructions::squared(C_0 + HALF * (xdir * C_X + ydir * C_Y), C_Z + HALF * (xdir * C_XZ + ydir * C_YZ));
   const Real bx2 = Reconstructions::squared(Bx + ydir * HALF * dBxdy, dBxdz);
   const Real by2 = Reconstructions::squared(By + xdir * HALF * dBydx, dBydz);

   return Wavespeeds(bx2 + by2 + bz2, rhom, p11, p22, p33, gridSpacing);
}

void fsdebugCheck([[maybe_unused]] const fsgrid::FsStencil& stencil, [[maybe_unused]] size_t len,
                  [[maybe_unused]] const char* file, [[maybe_unused]] uint32_t line) {
#ifdef DEBUG_FSOLVER
   const bool ok = stencil.ooo() < len && stencil.oom() < len && stencil.omo() < len && stencil.omm() < len &&
                   stencil.moo() < len && stencil.mom() < len && stencil.mmo() < len;

   if (!ok) {
      cerr << "Out-of-bounds access in " << file << ":" << line << std::endl;
      exit(1);
   }
#endif
}

/*!< Get min and max density of the four neighbors, used for limiting/clamping the results of interpolation. */
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

Real resistiveTerm(const auto& bgb, const auto& perb, const auto& dperb, Real rhoq, std::array<size_t, 2> indices,
                   std::array<Real, 2> spacing) {
   const auto x = bgb[fsgrids::bgbfield::BGBX] + perb[fsgrids::bfield::PERBX];
   const auto y = bgb[fsgrids::bgbfield::BGBY] + perb[fsgrids::bfield::PERBY];
   const auto z = bgb[fsgrids::bgbfield::BGBZ] + perb[fsgrids::bfield::PERBZ];

   return Parameters::resistivity * sqrt(x * x + y * y + z * z) / rhoq / physicalconstants::MU_0 *
          (dperb[indices[0]] / spacing[0] - dperb[indices[1]] / spacing[1]);
}

struct UpwindField {
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

   Real operator()() const {
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

struct CardinalIndices {
   size_t sw = 0;
   size_t se = 0;
   size_t nw = 0;
   size_t ne = 0;
};

struct DataArrays {
   const std::array<Real, fsgrids::bfield::N_BFIELD>& perb;
   const std::array<Real, fsgrids::dperb::N_DPERB>& dperb;
   const std::array<Real, fsgrids::moments::N_MOMENTS>& moments;
   const std::array<Real, fsgrids::dmoments::N_DMOMENTS>& dmoments;
   const std::array<Real, fsgrids::bgbfield::N_BGB>& bgb;

   DataArrays(fsgrids::perbspan perb,
              fsgrids::constdperbspan dperb,
              fsgrids::constmomentsspan moments,
              fsgrids::constdmomentsspan dmoments,
              fsgrids::constbgbspan bgb, size_t index)
       : perb(perb[index]), dperb(dperb[index]), moments(moments[index]), dmoments(dmoments[index]), bgb(bgb[index]) {}
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
 * \param perb fsGrid holding the perturbed B quantities
 * \param dperb fsGrid holding the derivatives of perturbed B
 * \param e fsGrid holding the electric field
 * \param ehall fsGrid holding the Hall contributions to the electric field
 * \param egradpe fsGrid holding the electron pressure gradient E field
 * \param moments fsGrid holding the moment quantities
 * \param dmoments fsGrid holding the derviatives of moments
 * \param bgb fsGrid holding the background B quantities
 * \param technical fsGrid holding technical information (such as boundary types)
 * \param stencil current cell's fsgrid stencil
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param gridSpacing fsgrid cell size in x,y,z
 */
void calculateEdgeElectricFieldX(fsgrids::perbspan perb,
                                 fsgrids::constdperbspan dperb,
                                 fsgrids::efieldspan e,
                                 fsgrids::constehallspan ehall,
                                 fsgrids::constegradpespan egradpe,
                                 fsgrids::constmomentsspan moments,
                                 fsgrids::constdmomentsspan dmoments,
                                 fsgrids::constbgbspan bgb,
                                 fsgrids::technicalspan technical, const fsgrid::FsStencil& stencil,
                                 int32_t RKCase, const std::array<Real, 3>& gridSpacing) {
   fsdebugCheck(stencil, perb.size(), __FILE__, __LINE__);

   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.
   Real ay_pos, ay_neg; // Max. characteristic velocities to y-direction
   Real az_pos, az_neg; // Max. characteristic velocities to z-direction
   Real maxV = 0.0;     // Max velocity for CFL purposes
   Real c_y, c_z;       // Wave speeds to yz-directions

   const CardinalIndices ci{stencil.ooo(), stencil.omo(), stencil.oom(), stencil.omm()};
   const DataArrays sw{perb, dperb, moments, dmoments, bgb, ci.sw};
   const DataArrays se{perb, dperb, moments, dmoments, bgb, ci.se};
   const DataArrays nw{perb, dperb, moments, dmoments, bgb, ci.nw};
   const DataArrays ne{perb, dperb, moments, dmoments, bgb, ci.ne};

   const auto rhomLimits = getRhomLimits({
       sw.moments,
       se.moments,
       nw.moments,
       ne.moments,
   });

   const Real By_S = sw.perb[fsgrids::bfield::PERBY] + sw.bgb[fsgrids::bgbfield::BGBY];
   const Real Bz_W = sw.perb[fsgrids::bfield::PERBZ] + sw.bgb[fsgrids::bgbfield::BGBZ];
   const Real Bz_E = se.perb[fsgrids::bfield::PERBZ] + se.bgb[fsgrids::bgbfield::BGBZ];
   const Real By_N = nw.perb[fsgrids::bfield::PERBY] + nw.bgb[fsgrids::bgbfield::BGBY];
   const Real perBy_S = sw.perb[fsgrids::bfield::PERBY];
   const Real perBz_W = sw.perb[fsgrids::bfield::PERBZ];
   const Real perBz_E = se.perb[fsgrids::bfield::PERBZ];
   const Real perBy_N = nw.perb[fsgrids::bfield::PERBY];
   Real Vy0 = sw.moments[fsgrids::moments::VY];
   Real Vz0 = sw.moments[fsgrids::moments::VZ];

   creal dBydx_S = sw.dperb[fsgrids::dperb::dPERBydx] + sw.bgb[fsgrids::bgbfield::dBGBydx];
   creal dBydz_S = sw.dperb[fsgrids::dperb::dPERBydz] + sw.bgb[fsgrids::bgbfield::dBGBydz];
   creal dBzdx_W = sw.dperb[fsgrids::dperb::dPERBzdx] + sw.bgb[fsgrids::bgbfield::dBGBzdx];
   creal dBzdy_W = sw.dperb[fsgrids::dperb::dPERBzdy] + sw.bgb[fsgrids::bgbfield::dBGBzdy];
   creal dBzdx_E = se.dperb[fsgrids::dperb::dPERBzdx] + se.bgb[fsgrids::bgbfield::dBGBzdx];
   creal dBzdy_E = se.dperb[fsgrids::dperb::dPERBzdy] + se.bgb[fsgrids::bgbfield::dBGBzdy];
   creal dBydx_N = nw.dperb[fsgrids::dperb::dPERBydx] + nw.bgb[fsgrids::bgbfield::dBGBydx];
   creal dBydz_N = nw.dperb[fsgrids::dperb::dPERBydz] + nw.bgb[fsgrids::bgbfield::dBGBydz];
   creal dperBydz_S = sw.dperb[fsgrids::dperb::dPERBydz];
   creal dperBydz_N = nw.dperb[fsgrids::dperb::dPERBydz];
   creal dperBzdy_W = sw.dperb[fsgrids::dperb::dPERBzdy];
   creal dperBzdy_E = se.dperb[fsgrids::dperb::dPERBzdy];

   // Ex and characteristic speeds on this cell:
   // 1st order terms:
   Real Ex_SW = By_S * Vz0 - Bz_W * Vy0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ex_SW += +HALF * ((By_S - HALF * dBydz_S) *
                         (-sw.dmoments[fsgrids::dmoments::dVzdy] - sw.dmoments[fsgrids::dmoments::dVzdz]) -
                     dBydz_S * Vz0 + SIXTH * dBydx_S * sw.dmoments[fsgrids::dmoments::dVzdx]);
   Ex_SW += -HALF * ((Bz_W - HALF * dBzdy_W) *
                         (-sw.dmoments[fsgrids::dmoments::dVydy] - sw.dmoments[fsgrids::dmoments::dVydz]) -
                     dBzdy_W * Vy0 + SIXTH * dBzdx_W * sw.dmoments[fsgrids::dmoments::dVydx]);
#endif
   size_t self = stencil.ooo();
   size_t nbr = stencil.poo();
   auto wavespeeds = calculateWaveSpeedYZ(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, By_S,
                                          Bz_W, dBydx_S, dBydz_S, dBzdx_W, dBzdy_W, MINUS, MINUS);
   c_y = wavespeeds.minVelocity();
   c_z = c_y;
   ay_neg = max(ZERO, -Vy0 + c_y);
   ay_pos = max(ZERO, +Vy0 + c_y);
   az_neg = max(ZERO, -Vz0 + c_z);
   az_pos = max(ZERO, +Vz0 + c_z);
   maxV = max(maxV, wavespeeds.cflSpeed(Vy0, Vz0));

   // Ex and characteristic speeds on SE neighbour:
   Vy0 = se.moments[fsgrids::moments::VY];
   Vz0 = se.moments[fsgrids::moments::VZ];

   // 1st order terms:
   Real Ex_SE = By_S * Vz0 - Bz_E * Vy0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ex_SE += +HALF * ((By_S - HALF * dBydz_S) *
                         (+se.dmoments[fsgrids::dmoments::dVzdy] - se.dmoments[fsgrids::dmoments::dVzdz]) -
                     dBydz_S * Vz0 + SIXTH * dBydx_S * se.dmoments[fsgrids::dmoments::dVzdx]);
   Ex_SE += -HALF * ((Bz_E + HALF * dBzdy_E) *
                         (+se.dmoments[fsgrids::dmoments::dVydy] - se.dmoments[fsgrids::dmoments::dVydz]) +
                     dBzdy_E * Vy0 + SIXTH * dBzdx_E * se.dmoments[fsgrids::dmoments::dVydx]);
#endif

   self = stencil.omo();
   nbr = stencil.pmo();
   wavespeeds = calculateWaveSpeedYZ(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, By_S,
                                     Bz_E, dBydx_S, dBydz_S, dBzdx_E, dBzdy_E, PLUS, MINUS);
   c_y = wavespeeds.minVelocity();
   c_z = c_y;
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   maxV = max(maxV, wavespeeds.cflSpeed(Vy0, Vz0));

   // Ex and characteristic speeds on NW neighbour:
   Vy0 = nw.moments[fsgrids::moments::VY];
   Vz0 = nw.moments[fsgrids::moments::VZ];

   // 1st order terms:
   Real Ex_NW = By_N * Vz0 - Bz_W * Vy0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ex_NW += +HALF * ((By_N + HALF * dBydz_N) *
                         (-nw.dmoments[fsgrids::dmoments::dVzdy] + nw.dmoments[fsgrids::dmoments::dVzdz]) +
                     dBydz_N * Vz0 + SIXTH * dBydx_N * nw.dmoments[fsgrids::dmoments::dVzdx]);
   Ex_NW += -HALF * ((Bz_W - HALF * dBzdy_W) *
                         (-nw.dmoments[fsgrids::dmoments::dVydy] + nw.dmoments[fsgrids::dmoments::dVydz]) -
                     dBzdy_W * Vy0 + SIXTH * dBzdx_W * nw.dmoments[fsgrids::dmoments::dVydx]);
#endif

   self = stencil.oom();
   nbr = stencil.pom();
   wavespeeds = calculateWaveSpeedYZ(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, By_N,
                                     Bz_W, dBydx_N, dBydz_N, dBzdx_W, dBzdy_W, MINUS, PLUS);
   c_y = wavespeeds.minVelocity();
   c_z = c_y;
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   maxV = max(maxV, wavespeeds.cflSpeed(Vy0, Vz0));

   // Ex and characteristic speeds on NE neighbour:
   Vy0 = ne.moments[fsgrids::moments::VY];
   Vz0 = ne.moments[fsgrids::moments::VZ];

   // 1st order terms:
   Real Ex_NE = By_N * Vz0 - Bz_E * Vy0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ex_NE += +HALF * ((By_N + HALF * dBydz_N) *
                         (+ne.dmoments[fsgrids::dmoments::dVzdy] + ne.dmoments[fsgrids::dmoments::dVzdz]) +
                     dBydz_N * Vz0 + SIXTH * dBydx_N * ne.dmoments[fsgrids::dmoments::dVzdx]);
   Ex_NE += -HALF * ((Bz_E + HALF * dBzdy_E) *
                         (+ne.dmoments[fsgrids::dmoments::dVydy] + ne.dmoments[fsgrids::dmoments::dVydz]) +
                     dBzdy_E * Vy0 + SIXTH * dBzdx_E * ne.dmoments[fsgrids::dmoments::dVydx]);
#endif

   self = stencil.omm();
   nbr = stencil.pmm();
   wavespeeds = calculateWaveSpeedYZ(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, By_N,
                                     Bz_E, dBydx_N, dBydz_N, dBzdx_E, dBzdy_E, PLUS, PLUS);
   c_y = wavespeeds.minVelocity();
   c_z = c_y;
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   maxV = max(maxV, wavespeeds.cflSpeed(Vy0, Vz0));

   // Resistive terms
   if (Parameters::resistivity > 0) {
      const std::array indices = {static_cast<size_t>(fsgrids::dperb::dPERBzdy),
                                  static_cast<size_t>(fsgrids::dperb::dPERBydz)};
      const std::array spacing = {gridSpacing[1], gridSpacing[2]};

      Ex_SW += resistiveTerm(sw.bgb, sw.perb, sw.dperb, sw.moments[fsgrids::moments::RHOQ], indices, spacing);
      Ex_SE += resistiveTerm(se.bgb, se.perb, se.dperb, se.moments[fsgrids::moments::RHOQ], indices, spacing);
      Ex_NW += resistiveTerm(nw.bgb, nw.perb, nw.dperb, nw.moments[fsgrids::moments::RHOQ], indices, spacing);
      Ex_NE += resistiveTerm(ne.bgb, ne.perb, ne.dperb, ne.moments[fsgrids::moments::RHOQ], indices, spacing);
   }

   // Hall terms
   if (Parameters::ohmHallTerm > 0) {
      Ex_SW += ehall[ci.sw][fsgrids::ehall::EXHALL_000_100];
      Ex_SE += ehall[ci.se][fsgrids::ehall::EXHALL_010_110];
      Ex_NW += ehall[ci.nw][fsgrids::ehall::EXHALL_001_101];
      Ex_NE += ehall[ci.ne][fsgrids::ehall::EXHALL_011_111];
   }

   // Electron pressure gradient terms
   if (Parameters::ohmGradPeTerm > 0) {
      Ex_SW += egradpe[ci.sw][fsgrids::egradpe::EXGRADPE];
      Ex_SE += egradpe[ci.se][fsgrids::egradpe::EXGRADPE];
      Ex_NW += egradpe[ci.nw][fsgrids::egradpe::EXGRADPE];
      Ex_NE += egradpe[ci.ne][fsgrids::egradpe::EXGRADPE];
   }

   // Calculate properly upwinded edge-averaged Ex:
   const UpwindField f(Ex_NE, Ex_SE, Ex_NW, Ex_SW, ay_pos, ay_neg, az_pos, az_neg, perBy_S, perBy_N, perBz_W, perBz_E,
                       dperBydz_S, dperBydz_N, dperBzdy_W, dperBzdy_E);
   e[ci.sw][fsgrids::efield::EX] = f();

   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      // compute maximum timestep for fieldsolver in this cell (CFL=1)
      Real min_dx = std::numeric_limits<Real>::max();
      min_dx = min(min_dx, gridSpacing[1]);
      min_dx = min(min_dx, gridSpacing[2]);
      // update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (maxV != ZERO) {
         auto& maxFsDt = technical[stencil.ooo()].maxFsDt;
         maxFsDt = min(maxFsDt, min_dx / maxV);
      }
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
 * \param perb fsGrid holding the perturbed B quantities
 * \param dperb fsGrid holding the derivatives of perturbed B
 * \param e fsGrid holding the electric field
 * \param ehall fsGrid holding the Hall contributions to the electric field
 * \param egradpe fsGrid holding the electron pressure gradient E field
 * \param moments fsGrid holding the moment quantities
 * \param dmoments fsGrid holding the derviatives of moments
 * \param bgb fsGrid holding the background B quantities
 * \param technical fsGrid holding technical information (such as boundary types)
 * \param stencil current cell's fsgrid stencil
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param gridSpacing fsgrid cell size in x,y,z
 */
void calculateEdgeElectricFieldY(fsgrids::perbspan perb,
                                 fsgrids::constdperbspan dperb,
                                 fsgrids::efieldspan e,
                                 fsgrids::constehallspan ehall,
                                 fsgrids::constegradpespan egradpe,
                                 fsgrids::constmomentsspan moments,
                                 fsgrids::constdmomentsspan dmoments,
                                 fsgrids::constbgbspan bgb,
                                 fsgrids::technicalspan technical, const fsgrid::FsStencil& stencil,
                                 int32_t RKCase, const std::array<Real, 3>& gridSpacing) {
   fsdebugCheck(stencil, perb.size(), __FILE__, __LINE__);

   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.
   Real ax_pos, ax_neg; // Max. characteristic velocities to x-direction
   Real az_pos, az_neg; // Max. characteristic velocities to z-direction
   Real maxV = 0.0;     // Max velocity for CFL purposes
   Real c_x, c_z;       // Wave speeds to xz-directions

   const CardinalIndices ci{stencil.ooo(), stencil.oom(), stencil.moo(), stencil.mom()};
   const DataArrays sw{perb, dperb, moments, dmoments, bgb, ci.sw};
   const DataArrays se{perb, dperb, moments, dmoments, bgb, ci.se};
   const DataArrays nw{perb, dperb, moments, dmoments, bgb, ci.nw};
   const DataArrays ne{perb, dperb, moments, dmoments, bgb, ci.ne};

   const auto rhomLimits = getRhomLimits({
       sw.moments,
       se.moments,
       nw.moments,
       ne.moments,
   });
   // Fetch required plasma parameters:
   const Real Bz_S = sw.perb[fsgrids::bfield::PERBZ] + sw.bgb[fsgrids::bgbfield::BGBZ];
   const Real Bx_W = sw.perb[fsgrids::bfield::PERBX] + sw.bgb[fsgrids::bgbfield::BGBX];
   const Real Bx_E = se.perb[fsgrids::bfield::PERBX] + se.bgb[fsgrids::bgbfield::BGBX];
   const Real Bz_N = nw.perb[fsgrids::bfield::PERBZ] + nw.bgb[fsgrids::bgbfield::BGBZ];
   const Real perBz_S = sw.perb[fsgrids::bfield::PERBZ];
   const Real perBx_W = sw.perb[fsgrids::bfield::PERBX];
   const Real perBx_E = se.perb[fsgrids::bfield::PERBX];
   const Real perBz_N = nw.perb[fsgrids::bfield::PERBZ];
   Real Vx0 = sw.moments[fsgrids::moments::VX];
   Real Vz0 = sw.moments[fsgrids::moments::VZ];

   creal dBxdy_W = sw.dperb[fsgrids::dperb::dPERBxdy] + sw.bgb[fsgrids::bgbfield::dBGBxdy];
   creal dBxdz_W = sw.dperb[fsgrids::dperb::dPERBxdz] + sw.bgb[fsgrids::bgbfield::dBGBxdz];
   creal dBzdx_S = sw.dperb[fsgrids::dperb::dPERBzdx] + sw.bgb[fsgrids::bgbfield::dBGBzdx];
   creal dBzdy_S = sw.dperb[fsgrids::dperb::dPERBzdy] + sw.bgb[fsgrids::bgbfield::dBGBzdy];
   creal dBxdy_E = se.dperb[fsgrids::dperb::dPERBxdy] + se.bgb[fsgrids::bgbfield::dBGBxdy];
   creal dBxdz_E = se.dperb[fsgrids::dperb::dPERBxdz] + se.bgb[fsgrids::bgbfield::dBGBxdz];
   creal dBzdx_N = nw.dperb[fsgrids::dperb::dPERBzdx] + nw.bgb[fsgrids::bgbfield::dBGBzdx];
   creal dBzdy_N = nw.dperb[fsgrids::dperb::dPERBzdy] + nw.bgb[fsgrids::bgbfield::dBGBzdy];
   creal dperBzdx_S = sw.dperb[fsgrids::dperb::dPERBzdx];
   creal dperBzdx_N = nw.dperb[fsgrids::dperb::dPERBzdx];
   creal dperBxdz_W = sw.dperb[fsgrids::dperb::dPERBxdz];
   creal dperBxdz_E = se.dperb[fsgrids::dperb::dPERBxdz];

   // Ey and characteristic speeds on this cell:
   // 1st order terms:
   Real Ey_SW = Bz_S * Vx0 - Bx_W * Vz0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms
   Ey_SW += +HALF * ((Bz_S - HALF * dBzdx_S) *
                         (-sw.dmoments[fsgrids::dmoments::dVxdx] - sw.dmoments[fsgrids::dmoments::dVxdz]) -
                     dBzdx_S * Vx0 + SIXTH * dBzdy_S * sw.dmoments[fsgrids::dmoments::dVxdy]);
   Ey_SW += -HALF * ((Bx_W - HALF * dBxdz_W) *
                         (-sw.dmoments[fsgrids::dmoments::dVzdx] - sw.dmoments[fsgrids::dmoments::dVzdz]) -
                     dBxdz_W * Vz0 + SIXTH * dBxdy_W * sw.dmoments[fsgrids::dmoments::dVzdy]);
#endif

   size_t self = stencil.ooo();
   size_t nbr = stencil.opo();
   auto wavespeeds = calculateWaveSpeedXZ(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, Bx_W,
                                          Bz_S, dBxdy_W, dBxdz_W, dBzdx_S, dBzdy_S, MINUS, MINUS);
   c_z = wavespeeds.minVelocity();
   c_x = c_z;
   az_neg = max(ZERO, -Vz0 + c_z);
   az_pos = max(ZERO, +Vz0 + c_z);
   ax_neg = max(ZERO, -Vx0 + c_x);
   ax_pos = max(ZERO, +Vx0 + c_x);
   maxV = max(maxV, wavespeeds.cflSpeed(Vz0, Vx0));

   // Ey and characteristic speeds on SE neighbour:
   Vx0 = se.moments[fsgrids::moments::VX];
   Vz0 = se.moments[fsgrids::moments::VZ];

   // 1st order terms:
   Real Ey_SE = Bz_S * Vx0 - Bx_E * Vz0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ey_SE += +HALF * ((Bz_S - HALF * dBzdx_S) *
                         (-se.dmoments[fsgrids::dmoments::dVxdx] + se.dmoments[fsgrids::dmoments::dVxdz]) -
                     dBzdx_S * Vx0 + SIXTH * dBzdy_S * se.dmoments[fsgrids::dmoments::dVxdy]);
   Ey_SE += -HALF * ((Bx_E + HALF * dBxdz_E) *
                         (-se.dmoments[fsgrids::dmoments::dVzdx] + se.dmoments[fsgrids::dmoments::dVzdz]) +
                     dBxdz_E * Vz0 + SIXTH * dBxdy_E * se.dmoments[fsgrids::dmoments::dVzdy]);
#endif

   self = stencil.oom();
   nbr = stencil.opm();
   wavespeeds = calculateWaveSpeedXZ(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, Bx_E,
                                     Bz_S, dBxdy_E, dBxdz_E, dBzdx_S, dBzdy_S, MINUS, PLUS);
   c_z = wavespeeds.minVelocity();
   c_x = c_z;
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   maxV = max(maxV, wavespeeds.cflSpeed(Vz0, Vx0));

   // Ey and characteristic speeds on NW neighbour:
   Vz0 = nw.moments[fsgrids::moments::VZ];
   Vx0 = nw.moments[fsgrids::moments::VX];

   // 1st order terms:
   Real Ey_NW = Bz_N * Vx0 - Bx_W * Vz0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ey_NW += +HALF * ((Bz_N + HALF * dBzdx_N) *
                         (+nw.dmoments[fsgrids::dmoments::dVxdx] - nw.dmoments[fsgrids::dmoments::dVxdz]) +
                     dBzdx_N * Vx0 + SIXTH * dBzdy_N * nw.dmoments[fsgrids::dmoments::dVxdy]);
   Ey_NW += -HALF * ((Bx_W - HALF * dBxdz_W) *
                         (+nw.dmoments[fsgrids::dmoments::dVzdx] - nw.dmoments[fsgrids::dmoments::dVzdz]) -
                     dBxdz_W * Vz0 + SIXTH * dBxdy_W * nw.dmoments[fsgrids::dmoments::dVzdy]);
#endif

   self = stencil.moo();
   nbr = stencil.mpo();
   wavespeeds = calculateWaveSpeedXZ(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, Bx_W,
                                     Bz_N, dBxdy_W, dBxdz_W, dBzdx_N, dBzdy_N, PLUS, MINUS);
   c_z = wavespeeds.minVelocity();
   c_x = c_z;
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   maxV = max(maxV, wavespeeds.cflSpeed(Vz0, Vx0));

   // Ey and characteristic speeds on NE neighbour:
   Vz0 = ne.moments[fsgrids::moments::VZ];
   Vx0 = ne.moments[fsgrids::moments::VX];

   // 1st order terms:
   Real Ey_NE = Bz_N * Vx0 - Bx_E * Vz0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ey_NE += +HALF * ((Bz_N + HALF * dBzdx_N) *
                         (+ne.dmoments[fsgrids::dmoments::dVxdx] + ne.dmoments[fsgrids::dmoments::dVxdz]) +
                     dBzdx_N * Vx0 + SIXTH * dBzdy_N * ne.dmoments[fsgrids::dmoments::dVxdy]);
   Ey_NE += -HALF * ((Bx_E + HALF * dBxdz_E) *
                         (+ne.dmoments[fsgrids::dmoments::dVzdx] + ne.dmoments[fsgrids::dmoments::dVzdz]) +
                     dBxdz_E * Vz0 + SIXTH * dBxdy_E * ne.dmoments[fsgrids::dmoments::dVzdy]);
#endif

   self = stencil.mom();
   nbr = stencil.mpm();
   wavespeeds = calculateWaveSpeedXZ(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, Bx_E,
                                     Bz_N, dBxdy_E, dBxdz_E, dBzdx_N, dBzdy_N, PLUS, PLUS);
   c_z = wavespeeds.minVelocity();
   c_x = c_z;
   az_neg = max(az_neg, -Vz0 + c_z);
   az_pos = max(az_pos, +Vz0 + c_z);
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   maxV = max(maxV, wavespeeds.cflSpeed(Vz0, Vx0));

   // Resistive terms
   if (Parameters::resistivity > 0) {
      const std::array indices = {static_cast<size_t>(fsgrids::dperb::dPERBxdz),
                                  static_cast<size_t>(fsgrids::dperb::dPERBzdx)};
      const std::array spacing = {gridSpacing[2], gridSpacing[0]};

      Ey_SW += resistiveTerm(sw.bgb, sw.perb, sw.dperb, sw.moments[fsgrids::moments::RHOQ], indices, spacing);
      Ey_SE += resistiveTerm(se.bgb, se.perb, se.dperb, se.moments[fsgrids::moments::RHOQ], indices, spacing);
      Ey_NW += resistiveTerm(nw.bgb, nw.perb, nw.dperb, nw.moments[fsgrids::moments::RHOQ], indices, spacing);
      Ey_NE += resistiveTerm(ne.bgb, ne.perb, ne.dperb, ne.moments[fsgrids::moments::RHOQ], indices, spacing);
   }

   // Hall terms
   if (Parameters::ohmHallTerm > 0) {
      Ey_SW += ehall[ci.sw][fsgrids::ehall::EYHALL_000_010];
      Ey_SE += ehall[ci.se][fsgrids::ehall::EYHALL_001_011];
      Ey_NW += ehall[ci.nw][fsgrids::ehall::EYHALL_100_110];
      Ey_NE += ehall[ci.ne][fsgrids::ehall::EYHALL_101_111];
   }

   // Electron pressure gradient terms
   if (Parameters::ohmGradPeTerm > 0) {
      Ey_SW += egradpe[ci.sw][fsgrids::egradpe::EYGRADPE];
      Ey_SE += egradpe[ci.se][fsgrids::egradpe::EYGRADPE];
      Ey_NW += egradpe[ci.nw][fsgrids::egradpe::EYGRADPE];
      Ey_NE += egradpe[ci.ne][fsgrids::egradpe::EYGRADPE];
   }

   // Calculate properly upwinded edge-averaged Ey:
   const UpwindField f(Ey_NE, Ey_SE, Ey_NW, Ey_SW, az_pos, az_neg, ax_pos, ax_neg, perBz_S, perBz_N, perBx_W, perBx_E,
                       dperBzdx_S, dperBzdx_N, dperBxdz_W, dperBxdz_E);
   e[ci.sw][fsgrids::efield::EY] = f();

   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      // compute maximum timestep for fieldsolver in this cell (CFL=1)
      Real min_dx = std::numeric_limits<Real>::max();
      min_dx = min(min_dx, gridSpacing[0]);
      min_dx = min(min_dx, gridSpacing[2]);
      // update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (maxV != ZERO) {
         auto& maxFsDt = technical[stencil.ooo()].maxFsDt;
         maxFsDt = min(maxFsDt, min_dx / maxV);
      }
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
 * \param perb fsGrid holding the perturbed B quantities
 * \param dperb fsGrid holding the derivatives of perturbed B
 * \param e fsGrid holding the electric field
 * \param ehall fsGrid holding the Hall contributions to the electric field
 * \param egradpe fsGrid holding the electron pressure gradient E field
 * \param moments fsGrid holding the moment quantities
 * \param dmoments fsGrid holding the derviatives of moments
 * \param bgb fsGrid holding the background B quantities
 * \param technical fsGrid holding technical information (such as boundary types)
 * \param stencil current cell's fsgrid stencil
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param gridSpacing fsgrid cell size in x,y,z
 */
void calculateEdgeElectricFieldZ(fsgrids::perbspan perb,
                                 fsgrids::constdperbspan dperb,
                                 fsgrids::efieldspan e,
                                 fsgrids::constehallspan ehall,
                                 fsgrids::constegradpespan egradpe,
                                 fsgrids::constmomentsspan moments,
                                 fsgrids::constdmomentsspan dmoments,
                                 fsgrids::constbgbspan bgb,
                                 fsgrids::technicalspan technical, const fsgrid::FsStencil& stencil,
                                 int32_t RKCase, const std::array<Real, 3>& gridSpacing) {
   fsdebugCheck(stencil, perb.size(), __FILE__, __LINE__);

   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.
   Real ax_pos, ax_neg; // Max. characteristic velocities to x-direction
   Real ay_pos, ay_neg; // Max. characteristic velocities to y-direction
   Real maxV = 0.0;     // Max velocity for CFL purposes
   Real c_x, c_y;       // Characteristic speeds to xy-directions

   const CardinalIndices ci{stencil.ooo(), stencil.moo(), stencil.omo(), stencil.mmo()};
   const DataArrays sw{perb, dperb, moments, dmoments, bgb, ci.sw};
   const DataArrays se{perb, dperb, moments, dmoments, bgb, ci.se};
   const DataArrays nw{perb, dperb, moments, dmoments, bgb, ci.nw};
   const DataArrays ne{perb, dperb, moments, dmoments, bgb, ci.ne};

   const auto rhomLimits = getRhomLimits({
       sw.moments,
       se.moments,
       nw.moments,
       ne.moments,
   });

   // Fetch needed plasma parameters/derivatives from the four cells:
   const Real Bx_S = sw.perb[fsgrids::bfield::PERBX] + sw.bgb[fsgrids::bgbfield::BGBX];
   const Real By_W = sw.perb[fsgrids::bfield::PERBY] + sw.bgb[fsgrids::bgbfield::BGBY];
   const Real By_E = se.perb[fsgrids::bfield::PERBY] + se.bgb[fsgrids::bgbfield::BGBY];
   const Real Bx_N = nw.perb[fsgrids::bfield::PERBX] + nw.bgb[fsgrids::bgbfield::BGBX];
   const Real perBx_S = sw.perb[fsgrids::bfield::PERBX];
   const Real perBy_W = sw.perb[fsgrids::bfield::PERBY];
   const Real perBy_E = se.perb[fsgrids::bfield::PERBY];
   const Real perBx_N = nw.perb[fsgrids::bfield::PERBX];
   Real Vx0 = sw.moments[fsgrids::moments::VX];
   Real Vy0 = sw.moments[fsgrids::moments::VY];

   creal dBxdy_S = sw.dperb[fsgrids::dperb::dPERBxdy] + sw.bgb[fsgrids::bgbfield::dBGBxdy];
   creal dBxdz_S = sw.dperb[fsgrids::dperb::dPERBxdz] + sw.bgb[fsgrids::bgbfield::dBGBxdz];
   creal dBydx_W = sw.dperb[fsgrids::dperb::dPERBydx] + sw.bgb[fsgrids::bgbfield::dBGBydx];
   creal dBydz_W = sw.dperb[fsgrids::dperb::dPERBydz] + sw.bgb[fsgrids::bgbfield::dBGBydz];
   creal dBydx_E = se.dperb[fsgrids::dperb::dPERBydx] + se.bgb[fsgrids::bgbfield::dBGBydx];
   creal dBydz_E = se.dperb[fsgrids::dperb::dPERBydz] + se.bgb[fsgrids::bgbfield::dBGBydz];
   creal dBxdy_N = nw.dperb[fsgrids::dperb::dPERBxdy] + nw.bgb[fsgrids::bgbfield::dBGBxdy];
   creal dBxdz_N = nw.dperb[fsgrids::dperb::dPERBxdz] + nw.bgb[fsgrids::bgbfield::dBGBxdz];
   creal dperBxdy_S = sw.dperb[fsgrids::dperb::dPERBxdy];
   creal dperBxdy_N = nw.dperb[fsgrids::dperb::dPERBxdy];
   creal dperBydx_W = sw.dperb[fsgrids::dperb::dPERBydx];
   creal dperBydx_E = se.dperb[fsgrids::dperb::dPERBydx];

   // Ez and characteristic speeds on SW cell:
   // 1st order terms:
   Real Ez_SW = Bx_S * Vy0 - By_W * Vx0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ez_SW += +HALF * ((Bx_S - HALF * dBxdy_S) *
                         (-sw.dmoments[fsgrids::dmoments::dVydx] - sw.dmoments[fsgrids::dmoments::dVydy]) -
                     dBxdy_S * Vy0 + SIXTH * dBxdz_S * sw.dmoments[fsgrids::dmoments::dVydz]);
   Ez_SW += -HALF * ((By_W - HALF * dBydx_W) *
                         (-sw.dmoments[fsgrids::dmoments::dVxdx] - sw.dmoments[fsgrids::dmoments::dVxdy]) -
                     dBydx_W * Vx0 + SIXTH * dBydz_W * sw.dmoments[fsgrids::dmoments::dVxdz]);
#endif

   size_t self = stencil.ooo();
   size_t nbr = stencil.oop();
   auto wavespeeds = calculateWaveSpeedXY(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, Bx_S,
                                          By_W, dBxdy_S, dBxdz_S, dBydx_W, dBydz_W, MINUS, MINUS);
   c_x = wavespeeds.minVelocity();
   c_y = c_x;
   ax_neg = max(ZERO, -Vx0 + c_x);
   ax_pos = max(ZERO, +Vx0 + c_x);
   ay_neg = max(ZERO, -Vy0 + c_y);
   ay_pos = max(ZERO, +Vy0 + c_y);
   maxV = max(maxV, wavespeeds.cflSpeed(Vx0, Vy0));

   // Ez and characteristic speeds on SE cell:
   Vx0 = se.moments[fsgrids::moments::VX];
   Vy0 = se.moments[fsgrids::moments::VY];

   // 1st order terms:
   Real Ez_SE = Bx_S * Vy0 - By_E * Vx0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ez_SE += +HALF * ((Bx_S - HALF * dBxdy_S) *
                         (+se.dmoments[fsgrids::dmoments::dVydx] - se.dmoments[fsgrids::dmoments::dVydy]) -
                     dBxdy_S * Vy0 + SIXTH * dBxdz_S * se.dmoments[fsgrids::dmoments::dVydz]);
   Ez_SE += -HALF * ((By_E + HALF * dBydx_E) *
                         (+se.dmoments[fsgrids::dmoments::dVxdx] - se.dmoments[fsgrids::dmoments::dVxdy]) +
                     dBydx_E * Vx0 + SIXTH * dBydz_E * se.dmoments[fsgrids::dmoments::dVxdz]);
#endif

   self = stencil.moo();
   nbr = stencil.mop();
   wavespeeds = calculateWaveSpeedXY(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, Bx_S,
                                     By_E, dBxdy_S, dBxdz_S, dBydx_E, dBydz_E, PLUS, MINUS);
   c_x = wavespeeds.minVelocity();
   c_y = c_x;
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   maxV = max(maxV, wavespeeds.cflSpeed(Vx0, Vy0));

   // Ez and characteristic speeds on NW cell:
   Vx0 = nw.moments[fsgrids::moments::VX];
   Vy0 = nw.moments[fsgrids::moments::VY];

   // 1st order terms:
   Real Ez_NW = Bx_N * Vy0 - By_W * Vx0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ez_NW += +HALF * ((Bx_N + HALF * dBxdy_N) *
                         (-nw.dmoments[fsgrids::dmoments::dVydx] + nw.dmoments[fsgrids::dmoments::dVydy]) +
                     dBxdy_N * Vy0 + SIXTH * dBxdz_N * nw.dmoments[fsgrids::dmoments::dVydz]);
   Ez_NW += -HALF * ((By_W - HALF * dBydx_W) *
                         (-nw.dmoments[fsgrids::dmoments::dVxdx] + nw.dmoments[fsgrids::dmoments::dVxdy]) -
                     dBydx_W * Vx0 + SIXTH * dBydz_W * nw.dmoments[fsgrids::dmoments::dVxdz]);
#endif

   self = stencil.omo();
   nbr = stencil.omp();
   wavespeeds = calculateWaveSpeedXY(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, Bx_N,
                                     By_W, dBxdy_N, dBxdz_N, dBydx_W, dBydz_W, MINUS, PLUS);
   c_x = wavespeeds.minVelocity();
   c_y = c_x;
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   maxV = max(maxV, wavespeeds.cflSpeed(Vx0, Vy0));

   // Ez and characteristic speeds on NE cell:
   Vx0 = ne.moments[fsgrids::moments::VX];
   Vy0 = ne.moments[fsgrids::moments::VY];

   // 1st order terms:
   Real Ez_NE = Bx_N * Vy0 - By_E * Vx0;

#ifndef FS_1ST_ORDER_SPACE
   // 2nd order terms:
   Ez_NE += +HALF * ((Bx_N + HALF * dBxdy_N) *
                         (+ne.dmoments[fsgrids::dmoments::dVydx] + ne.dmoments[fsgrids::dmoments::dVydy]) +
                     dBxdy_N * Vy0 + SIXTH * dBxdz_N * ne.dmoments[fsgrids::dmoments::dVydz]);
   Ez_NE += -HALF * ((By_E + HALF * dBydx_E) *
                         (+ne.dmoments[fsgrids::dmoments::dVxdx] + ne.dmoments[fsgrids::dmoments::dVxdy]) +
                     dBydx_E * Vx0 + SIXTH * dBydz_E * ne.dmoments[fsgrids::dmoments::dVxdz]);
#endif

   self = stencil.mmo();
   nbr = stencil.mmp();
   wavespeeds = calculateWaveSpeedXY(perb, moments, dperb, dmoments, bgb, gridSpacing, rhomLimits, self, nbr, Bx_N,
                                     By_E, dBxdy_N, dBxdz_N, dBydx_E, dBydz_E, PLUS, PLUS);
   c_x = wavespeeds.minVelocity();
   c_y = c_x;
   ax_neg = max(ax_neg, -Vx0 + c_x);
   ax_pos = max(ax_pos, +Vx0 + c_x);
   ay_neg = max(ay_neg, -Vy0 + c_y);
   ay_pos = max(ay_pos, +Vy0 + c_y);
   maxV = max(maxV, wavespeeds.cflSpeed(Vx0, Vy0));

   // Resistive terms
   if (Parameters::resistivity > 0) {
      const std::array indices = {static_cast<size_t>(fsgrids::dperb::dPERBydx),
                                  static_cast<size_t>(fsgrids::dperb::dPERBxdy)};
      const std::array spacing = {gridSpacing[0], gridSpacing[1]};

      Ez_SW += resistiveTerm(sw.bgb, sw.perb, sw.dperb, sw.moments[fsgrids::moments::RHOQ], indices, spacing);
      Ez_SE += resistiveTerm(se.bgb, se.perb, se.dperb, se.moments[fsgrids::moments::RHOQ], indices, spacing);
      Ez_NW += resistiveTerm(nw.bgb, nw.perb, nw.dperb, nw.moments[fsgrids::moments::RHOQ], indices, spacing);
      Ez_NE += resistiveTerm(ne.bgb, ne.perb, ne.dperb, ne.moments[fsgrids::moments::RHOQ], indices, spacing);
   }

   // Hall terms
   if (Parameters::ohmHallTerm > 0) {
      Ez_SW += ehall[ci.sw][fsgrids::ehall::EZHALL_000_001];
      Ez_SE += ehall[ci.se][fsgrids::ehall::EZHALL_100_101];
      Ez_NW += ehall[ci.nw][fsgrids::ehall::EZHALL_010_011];
      Ez_NE += ehall[ci.ne][fsgrids::ehall::EZHALL_110_111];
   }

   // Electron pressure gradient terms
   if (Parameters::ohmGradPeTerm > 0) {
      Ez_SW += egradpe[ci.sw][fsgrids::egradpe::EZGRADPE];
      Ez_SE += egradpe[ci.se][fsgrids::egradpe::EZGRADPE];
      Ez_NW += egradpe[ci.nw][fsgrids::egradpe::EZGRADPE];
      Ez_NE += egradpe[ci.ne][fsgrids::egradpe::EZGRADPE];
   }

   // Calculate properly upwinded edge-averaged Ez:
   const UpwindField f(Ez_NE, Ez_SE, Ez_NW, Ez_SW, ax_pos, ax_neg, ay_pos, ay_neg, perBx_S, perBx_N, perBy_W, perBy_E,
                       dperBxdy_S, dperBxdy_N, dperBydx_W, dperBydx_E);
   e[ci.sw][fsgrids::efield::EZ] = f();

   if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
      // compute maximum timestep for fieldsolver in this cell (CFL=1)
      Real min_dx = std::numeric_limits<Real>::max();
      min_dx = min(min_dx, gridSpacing[0]);
      min_dx = min(min_dx, gridSpacing[1]);
      // update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
      if (maxV != ZERO) {
         auto& maxFsDt = technical[stencil.ooo()].maxFsDt;
         maxFsDt = min(maxFsDt, min_dx / maxV);
      }
   }
}

/*! \brief Electric field propagation function.
 *
 * Calls the general or the system boundary electric field propagation functions.
 *
 * \param perb fsGrid holding the perturbed B quantities
 * \param dperb fsGrid holding the derivatives of perturbed B
 * \param e fsGrid holding the electric field
 * \param ehall fsGrid holding the Hall contributions to the electric field
 * \param egradpe fsGrid holding the electron pressure gradient E field
 * \param moments fsGrid holding the moment quantities
 * \param dmoments fsGrid holding the derviatives of moments
 * \param bgb fsGrid holding the background B quantities
 * \param technical fsGrid holding technical information (such as boundary types)
 * \param stencil current cell's fsgrid stencil
 * \param gridSpacing fsgrid cell size in x,y,z
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 *
 * \sa calculateUpwindedElectricFieldSimple calculateEdgeElectricFieldX calculateEdgeElectricFieldY
 * calculateEdgeElectricFieldZ
 *
 */
void calculateElectricField(fsgrids::perbspan perb,
                            fsgrids::constdperbspan dperb,
                            fsgrids::efieldspan e,
                            fsgrids::constehallspan ehall,
                            fsgrids::constegradpespan egradpe,
                            fsgrids::constmomentsspan moments,
                            fsgrids::constdmomentsspan dmoments,
                            fsgrids::constbgbspan bgb,
                            fsgrids::technicalspan technical, const fsgrid::FsStencil& stencil,
                            const std::array<Real, 3>& gridSpacing, SysBoundary& sysBoundaries, int32_t RKCase) {
   cuint cellSysBoundaryFlag = technical[stencil.ooo()].sysBoundaryFlag;
   cuint bitfield = technical[stencil.ooo()].SOLVE;

   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
       cellSysBoundaryFlag == sysboundarytype::OUTER_BOUNDARY_PADDING) {
      return;
   }

   if ((bitfield & compute::EX) == compute::EX) {
      calculateEdgeElectricFieldX(perb, dperb, e, ehall, egradpe, moments, dmoments, bgb, technical, stencil, RKCase,
                                  gridSpacing);
   } else {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondElectricField(e, stencil, 0);
   }

   if ((bitfield & compute::EY) == compute::EY) {
      calculateEdgeElectricFieldY(perb, dperb, e, ehall, egradpe, moments, dmoments, bgb, technical, stencil, RKCase,
                                  gridSpacing);
   } else {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondElectricField(e, stencil, 1);
   }

   if ((bitfield & compute::EZ) == compute::EZ) {
      calculateEdgeElectricFieldZ(perb, dperb, e, ehall, egradpe, moments, dmoments, bgb, technical, stencil, RKCase,
                                  gridSpacing);
   } else {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondElectricField(e, stencil, 2);
   }
}

/*! \brief High-level electric field computation function.
 *
 * Transfers the derivatives, calculates the edge electric fields and transfers the new electric fields.
 *
 * \param perb fsGrid holding the perturbed B quantities at runge-kutta t=0
 * \param perbdt2 fsGrid holding the perturbed B quantities at runge-kutta t=0.5
 * \param e fsGrid holding the Electric field quantities at runge-kutta t=0
 * \param edt2 fsGrid holding the Electric field quantities at runge-kutta t=0.5
 * \param ehall fsGrid holding the Hall contributions to the electric field
 * \param egradpe fsGrid holding the electron pressure gradient E field
 * \param egradpedt2 fsGrid holding the electron pressure gradient E field
 * \param moments fsGrid holding the moment quantities at runge-kutta t=0
 * \param momentsdt2 fsGrid holding the moment quantities at runge-kutta t=0.5
 * \param dperb fsGrid holding the derivatives of perturbed B
 * \param dmoments fsGrid holding the derivatives of moments
 * \param dmomentsdt2 fsGrid holding the derivatives of moments at runge-kutta t=0.5
 * \param bgb fsGrid holding the background B quantities
 * \param technical fsGrid holding technical information (such as boundary types)
 * \param fsgrid fsgrids container
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param communicateEGradPeOrMomentsDerivatives Boolean flag whether grad(Pe) electric field or moments need a ghost update
 *
 * \sa calculateElectricField calculateEdgeElectricFieldX calculateEdgeElectricFieldY calculateEdgeElectricFieldZ
 */
void calculateUpwindedElectricFieldSimple(fsgrids::perbspan perb,
                                          fsgrids::perbspan perbdt2,
                                          fsgrids::efieldspan e,
                                          fsgrids::efieldspan edt2,
                                          fsgrids::ehallspan ehall,
                                          fsgrids::egradpespan egradpe,
                                          fsgrids::egradpespan egradpedt2,
                                          fsgrids::momentsspan moments,
                                          fsgrids::momentsspan momentsdt2,
                                          fsgrids::dperbspan dperb,
                                          fsgrids::dmomentsspan dmoments,
                                          fsgrids::dmomentsspan dmomentsdt2,
                                          fsgrids::bgbspan bgb,
                                          fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
                                          SysBoundary& sysBoundaries, int32_t RKCase,
                                          const bool communicateEGradPeOrMomentsDerivatives) {
   const size_t numCells = fsgrid.getNumCells();

   if (not(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2)) {
      perb = perbdt2;
      e = edt2;
      egradpe = egradpedt2;
      moments = momentsdt2;
      dmoments = dmomentsdt2;
   }

   phiprof::Timer upwindedETimer{"Calculate upwinded electric field"};

   phiprof::Timer mpiTimer{"Electric field ghost updates MPI", {"MPI"}};

   // Update ghosts if necessary, unless previous terms have already updated them
   if (P::ohmHallTerm > 0) {
      fsgrid.updateGhostCells(ehall);
   }

   if (P::ohmGradPeTerm > 0 && communicateEGradPeOrMomentsDerivatives) {
      fsgrid.updateGhostCells(egradpe);
   }

   if (P::ohmHallTerm == 0) {
      fsgrid.updateGhostCells(dperb);
   }

   if (P::ohmHallTerm == 0 && P::ohmGradPeTerm == 0 && communicateEGradPeOrMomentsDerivatives) {
      fsgrid.updateGhostCells(dmoments);
   }

   mpiTimer.stop();
   // Calculate upwinded electric field
   fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                       phiprof::initializeTimer("Electric field compute cells"), technical,
                       [=, &sysBoundaries](const fsgrid::Coordinates &coordinates, const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
                          calculateElectricField(perb, dperb, e, ehall, egradpe, moments, dmoments, bgb, technical, stencil,
                                                 coordinates.physicalGridSpacing, sysBoundaries, RKCase);
                       });

   mpiTimer.start();
   // Exchange electric field with neighbouring processes
   fsgrid.updateGhostCells(e);
   mpiTimer.stop();

   upwindedETimer.stop(numCells, "Spatial Cells");
}
