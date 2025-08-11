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

#include "fs_common.h"
#include "../fieldtracing/fieldtracing.h"

/*! \brief Helper function
 *
 * Divides the first value by the second or returns zero if the denominator is zero.
 *
 * \param numerator Numerator
 * \param denominator Denominator
 */
Real divideIfNonZero(creal numerator, creal denominator) {
   if (denominator <= 0.0) {
      return 0.0;
   } else {
      return numerator / denominator;
   }
}

/*! \brief Low-level helper function.
 *
 * Computes the reconstruction coefficients used for field component reconstruction.
 * Only implemented for 2nd and 3rd order.
 *
 * \param perb fsGrid holding the perturbed B quantities
 * \param dperb fsGrid holding the derivatives of perturbed B
 * \param stencil fsgrid stencil for the current cell
 * \param reconstructionOrder Reconstruction order of the fields after Balsara 2009, 2 used for BVOL, 3 used for
 * 2nd-order Hall term calculations.
 */
std::array<Real, Rec::N_REC_COEFFICIENTS>
reconstructionCoefficients(fsgrids::perbspan perb,
                           fsgrids::constdperbspan dperb,
                           const fsgrid::FsStencil& stencil, Real reconstructionOrder) {
   std::array<Real, Rec::N_REC_COEFFICIENTS> perturbedResult;
   const auto ooo = stencil.ooo();
   const auto poo = stencil.poo();
   const auto opo = stencil.opo();
   const auto oop = stencil.oop();
   const std::array<Real, fsgrids::dperb::N_DPERB>& der_i1j1k1 = dperb[ooo];
   const std::array<Real, fsgrids::bfield::N_BFIELD>& cep_i1j1k1 = perb[ooo];

   const bool pooOk = stencil.cellExists(1, 0, 0);
   const bool opoOk = stencil.cellExists(0, 1, 0);
   const bool oopOk = stencil.cellExists(0, 0, 1);
   const std::array<Real, fsgrids::bfield::N_BFIELD>& cep_i2j1k1 = pooOk ? perb[poo] : cep_i1j1k1;
   const std::array<Real, fsgrids::bfield::N_BFIELD>& cep_i1j2k1 = opoOk ? perb[opo] : cep_i1j1k1;
   const std::array<Real, fsgrids::bfield::N_BFIELD>& cep_i1j1k2 = oopOk ? perb[oop] : cep_i1j1k1;

#ifndef FS_1ST_ORDER_SPACE
   // Create a dummy array for containing zero values for derivatives on non-existing cells:
   std::array<Real, fsgrids::dperb::N_DPERB> dummyDerivatives;
   dummyDerivatives.fill(0.0);

   // Fetch neighbour cell derivatives, or in case the neighbour does not
   // exist, use dummyDerivatives array:
   const std::array<Real, fsgrids::dperb::N_DPERB>& der_i2j1k1 = pooOk ? dperb[poo] : dummyDerivatives;
   const std::array<Real, fsgrids::dperb::N_DPERB>& der_i1j2k1 = opoOk ? dperb[opo] : dummyDerivatives;
   const std::array<Real, fsgrids::dperb::N_DPERB>& der_i1j1k2 = oopOk ? dperb[oop] : dummyDerivatives;

   // Calculate 3rd order reconstruction coefficients:
   if (reconstructionOrder == 2) {
      perturbedResult[Rec::a_yy] = 0.0;
      perturbedResult[Rec::a_zz] = 0.0;
      perturbedResult[Rec::a_yz] = 0.0;
      perturbedResult[Rec::a_xxx] = 0.0;
      perturbedResult[Rec::a_xxy] = 0.0;
      perturbedResult[Rec::a_xxz] = 0.0;
      perturbedResult[Rec::a_xyy] = 0.0;
      perturbedResult[Rec::a_xyz] = 0.0;
      perturbedResult[Rec::a_xzz] = 0.0;
      perturbedResult[Rec::b_xx] = 0.0;
      perturbedResult[Rec::b_xz] = 0.0;
      perturbedResult[Rec::b_zz] = 0.0;
      perturbedResult[Rec::b_xxy] = 0.0;
      perturbedResult[Rec::b_xyy] = 0.0;
      perturbedResult[Rec::b_xyz] = 0.0;
      perturbedResult[Rec::b_yyy] = 0.0;
      perturbedResult[Rec::b_yyz] = 0.0;
      perturbedResult[Rec::b_yzz] = 0.0;
      perturbedResult[Rec::c_xx] = 0.0;
      perturbedResult[Rec::c_xy] = 0.0;
      perturbedResult[Rec::c_yy] = 0.0;
      perturbedResult[Rec::c_xxz] = 0.0;
      perturbedResult[Rec::c_xyz] = 0.0;
      perturbedResult[Rec::c_xzz] = 0.0;
      perturbedResult[Rec::c_yyz] = 0.0;
      perturbedResult[Rec::c_yzz] = 0.0;
      perturbedResult[Rec::c_zzz] = 0.0;
   } else if (reconstructionOrder == 3) {
      perturbedResult[Rec::a_yy] =
          HALF * (der_i2j1k1[fsgrids::dperb::dPERBxdyy] + der_i1j1k1[fsgrids::dperb::dPERBxdyy]);
      perturbedResult[Rec::a_zz] =
          HALF * (der_i2j1k1[fsgrids::dperb::dPERBxdzz] + der_i1j1k1[fsgrids::dperb::dPERBxdzz]);
      perturbedResult[Rec::a_yz] =
          HALF * (der_i2j1k1[fsgrids::dperb::dPERBxdyz] + der_i1j1k1[fsgrids::dperb::dPERBxdyz]);
      perturbedResult[Rec::a_xyy] = (der_i2j1k1[fsgrids::dperb::dPERBxdyy] - der_i1j1k1[fsgrids::dperb::dPERBxdyy]);
      perturbedResult[Rec::a_xyz] = (der_i2j1k1[fsgrids::dperb::dPERBxdyz] - der_i1j1k1[fsgrids::dperb::dPERBxdyz]);
      perturbedResult[Rec::a_xzz] = (der_i2j1k1[fsgrids::dperb::dPERBxdzz] - der_i1j1k1[fsgrids::dperb::dPERBxdzz]);

      perturbedResult[Rec::b_xx] =
          HALF * (der_i1j2k1[fsgrids::dperb::dPERBydxx] + der_i1j1k1[fsgrids::dperb::dPERBydxx]);
      perturbedResult[Rec::b_xz] =
          HALF * (der_i1j2k1[fsgrids::dperb::dPERBydxz] + der_i1j1k1[fsgrids::dperb::dPERBydxz]);
      perturbedResult[Rec::b_zz] =
          HALF * (der_i1j2k1[fsgrids::dperb::dPERBydzz] + der_i1j1k1[fsgrids::dperb::dPERBydzz]);
      perturbedResult[Rec::b_xxy] = (der_i1j2k1[fsgrids::dperb::dPERBydxx] - der_i1j1k1[fsgrids::dperb::dPERBydxx]);
      perturbedResult[Rec::b_xyz] = (der_i1j2k1[fsgrids::dperb::dPERBydxz] - der_i1j1k1[fsgrids::dperb::dPERBydxz]);
      perturbedResult[Rec::b_yzz] = (der_i1j2k1[fsgrids::dperb::dPERBydzz] - der_i1j1k1[fsgrids::dperb::dPERBydzz]);

      perturbedResult[Rec::c_xx] =
          HALF * (der_i1j1k2[fsgrids::dperb::dPERBzdxx] + der_i1j1k1[fsgrids::dperb::dPERBzdxx]);
      perturbedResult[Rec::c_xy] =
          HALF * (der_i1j1k2[fsgrids::dperb::dPERBzdxy] + der_i1j1k1[fsgrids::dperb::dPERBzdxy]);
      perturbedResult[Rec::c_yy] =
          HALF * (der_i1j1k2[fsgrids::dperb::dPERBzdyy] + der_i1j1k1[fsgrids::dperb::dPERBzdyy]);
      perturbedResult[Rec::c_xxz] = (der_i1j1k2[fsgrids::dperb::dPERBzdxx] - der_i1j1k1[fsgrids::dperb::dPERBzdxx]);
      perturbedResult[Rec::c_xyz] = (der_i1j1k2[fsgrids::dperb::dPERBzdxy] - der_i1j1k1[fsgrids::dperb::dPERBzdxy]);
      perturbedResult[Rec::c_yyz] = (der_i1j1k2[fsgrids::dperb::dPERBzdyy] - der_i1j1k1[fsgrids::dperb::dPERBzdyy]);

      perturbedResult[Rec::a_xxx] = -THIRD * (perturbedResult[Rec::b_xxy] + perturbedResult[Rec::c_xxz]);
      perturbedResult[Rec::a_xxy] = -FOURTH * perturbedResult[Rec::c_xyz];
      perturbedResult[Rec::a_xxz] = -FOURTH * perturbedResult[Rec::b_xyz];

      perturbedResult[Rec::b_xyy] = -FOURTH * perturbedResult[Rec::c_xyz];
      perturbedResult[Rec::b_yyy] = -THIRD * (perturbedResult[Rec::c_yyz] + perturbedResult[Rec::a_xyy]);
      perturbedResult[Rec::b_yyz] = -FOURTH * perturbedResult[Rec::a_xyz];

      perturbedResult[Rec::c_xzz] = -FOURTH * perturbedResult[Rec::b_xyz];
      perturbedResult[Rec::c_yzz] = -FOURTH * perturbedResult[Rec::a_xyz];
      perturbedResult[Rec::c_zzz] = -THIRD * (perturbedResult[Rec::a_xzz] + perturbedResult[Rec::b_yzz]);
   } else {
      cerr << __FILE__ << ":" << __LINE__ << ":"
           << " Not coded yet!" << endl;
      abort();
   }

   // Calculate 2nd order reconstruction coefficients:
   perturbedResult[Rec::a_xy] = der_i2j1k1[fsgrids::dperb::dPERBxdy] - der_i1j1k1[fsgrids::dperb::dPERBxdy];
   perturbedResult[Rec::a_xz] = der_i2j1k1[fsgrids::dperb::dPERBxdz] - der_i1j1k1[fsgrids::dperb::dPERBxdz];
   perturbedResult[Rec::a_y] = HALF * (der_i2j1k1[fsgrids::dperb::dPERBxdy] + der_i1j1k1[fsgrids::dperb::dPERBxdy]) -
                               SIXTH * perturbedResult[Rec::a_xxy];
   perturbedResult[Rec::a_z] = HALF * (der_i2j1k1[fsgrids::dperb::dPERBxdz] + der_i1j1k1[fsgrids::dperb::dPERBxdz]) -
                               SIXTH * perturbedResult[Rec::a_xxz];

   perturbedResult[Rec::b_xy] = der_i1j2k1[fsgrids::dperb::dPERBydx] - der_i1j1k1[fsgrids::dperb::dPERBydx];
   perturbedResult[Rec::b_yz] = der_i1j2k1[fsgrids::dperb::dPERBydz] - der_i1j1k1[fsgrids::dperb::dPERBydz];
   perturbedResult[Rec::b_x] = HALF * (der_i1j2k1[fsgrids::dperb::dPERBydx] + der_i1j1k1[fsgrids::dperb::dPERBydx]) -
                               SIXTH * perturbedResult[Rec::b_xyy];
   perturbedResult[Rec::b_z] = HALF * (der_i1j2k1[fsgrids::dperb::dPERBydz] + der_i1j1k1[fsgrids::dperb::dPERBydz]) -
                               SIXTH * perturbedResult[Rec::b_yyz];

   perturbedResult[Rec::c_xz] = der_i1j1k2[fsgrids::dperb::dPERBzdx] - der_i1j1k1[fsgrids::dperb::dPERBzdx];
   perturbedResult[Rec::c_yz] = der_i1j1k2[fsgrids::dperb::dPERBzdy] - der_i1j1k1[fsgrids::dperb::dPERBzdy];
   perturbedResult[Rec::c_x] = HALF * (der_i1j1k2[fsgrids::dperb::dPERBzdx] + der_i1j1k1[fsgrids::dperb::dPERBzdx]) -
                               SIXTH * perturbedResult[Rec::c_xzz];
   perturbedResult[Rec::c_y] = HALF * (der_i1j1k2[fsgrids::dperb::dPERBzdy] + der_i1j1k1[fsgrids::dperb::dPERBzdy]) -
                               SIXTH * perturbedResult[Rec::c_yzz];

   perturbedResult[Rec::a_xx] = -HALF * (perturbedResult[Rec::b_xy] + perturbedResult[Rec::c_xz]);
   perturbedResult[Rec::b_yy] = -HALF * (perturbedResult[Rec::a_xy] + perturbedResult[Rec::c_yz]);
   perturbedResult[Rec::c_zz] = -HALF * (perturbedResult[Rec::a_xz] + perturbedResult[Rec::b_yz]);

   perturbedResult[Rec::a_x] =
       cep_i2j1k1[fsgrids::bfield::PERBX] - cep_i1j1k1[fsgrids::bfield::PERBX] - TENTH * perturbedResult[Rec::a_xxx];
   perturbedResult[Rec::b_y] =
       cep_i1j2k1[fsgrids::bfield::PERBY] - cep_i1j1k1[fsgrids::bfield::PERBY] - TENTH * perturbedResult[Rec::b_yyy];
   perturbedResult[Rec::c_z] =
       cep_i1j1k2[fsgrids::bfield::PERBZ] - cep_i1j1k1[fsgrids::bfield::PERBZ] - TENTH * perturbedResult[Rec::c_zzz];

#else
   perturbedResult.fill(0.0);
#endif

   // Calculate 1st order reconstruction coefficients:
   perturbedResult[Rec::a_0] = HALF * (cep_i2j1k1[fsgrids::bfield::PERBX] + cep_i1j1k1[fsgrids::bfield::PERBX]) -
                               SIXTH * perturbedResult[Rec::a_xx];
   perturbedResult[Rec::b_0] = HALF * (cep_i1j2k1[fsgrids::bfield::PERBY] + cep_i1j1k1[fsgrids::bfield::PERBY]) -
                               SIXTH * perturbedResult[Rec::b_yy];
   perturbedResult[Rec::c_0] = HALF * (cep_i1j1k2[fsgrids::bfield::PERBZ] + cep_i1j1k1[fsgrids::bfield::PERBZ]) -
                               SIXTH * perturbedResult[Rec::c_zz];

   return perturbedResult;
}

/*! Interpolate perturbed B to arbitrary x,y,z in cell
 *  Uses the reconstruction coefficients and equations from
 *  Divergence-free reconstruction of magnetic fields and WENO schemes for magnetohydrodynamics
 *  D.S. Balsara, J. Comp. Phys., 228, 2009
 *  doi:10.1016/j.jcp.2009.03.038
 *
 * \param perb perturbed B fsGrid
 * \param dperb perturbed B derivatives fsGrid
 * \param technical technical fsGrid
 * \param reconstructionCoefficientsCache map for cached values
 * \param i local fsGrid x-index
 * \param j local fsGrid y-index
 * \param k local fsGrid z-index
 * \param x 3D global simulation x,y,z coordinates of point to interpolate to
 */
std::array<Real, 3> interpolatePerturbedB(
    fsgrids::perbspan perb,
    fsgrids::constdperbspan dperb,
    fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
    std::map<std::array<int, 3>, std::array<Real, Rec::N_REC_COEFFICIENTS>>& reconstructionCoefficientsCache, cint i,
    cint j, cint k, const std::array<Real, 3> x) {
   const auto stencil = fsgrid.makeStencil(i, j, k);

   cuint cellSysBoundaryFlag = technical[stencil.ooo()].sysBoundaryFlag;
   if (cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
      return {0, 0, 0};
   }

   // Balsara reconstruction formulas: x,y,z are in [-1/2, 1/2] local coordinates
   std::array<Real, 3> xLocal = fsgrid.physicalToCellFractional(x[0], x[1], x[2]);
   xLocal[0] -= 0.5;
   xLocal[1] -= 0.5;
   xLocal[2] -= 0.5;

   if (fabs(xLocal[0]) > 0.5 || fabs(xLocal[1]) > 0.5 || fabs(xLocal[2]) > 0.5) {
      cerr << __FILE__ << ":" << __LINE__ << ": Coordinate (" << xLocal[0] << "," << xLocal[1] << "," << xLocal[2]
           << ")  outside of this cell!" << endl;
      abort();
   }

   const std::array<Real, Rec::N_REC_COEFFICIENTS> rc = [&stencil, &reconstructionCoefficientsCache, &perb, &dperb]() {
      const std::array<int, 3> cellIds = {stencil.i, stencil.j, stencil.k};
      if (FieldTracing::fieldTracingParameters.useCache) {
#pragma omp critical
         {
            // Reconstruction order of the fields after Balsara 2009, 2 used for general B, but 3 used here to
            // allow for cache reuse, see interpolatePerturbedJ below
            if (reconstructionCoefficientsCache.find(cellIds) == reconstructionCoefficientsCache.end()) {
               const auto rc = reconstructionCoefficients(perb, dperb, stencil, 3);
               reconstructionCoefficientsCache.insert({cellIds, rc});
            }
         }

         return reconstructionCoefficientsCache.at(cellIds);
      } else {
         return reconstructionCoefficients(perb, dperb, stencil, 3);
      }
   }();

   return {
       // Eq. (7) Balsara 2009
       rc[Rec::a_0] + rc[Rec::a_x] * xLocal[0] + rc[Rec::a_y] * xLocal[1] + rc[Rec::a_z] * xLocal[2] +
           rc[Rec::a_xx] * (xLocal[0] * xLocal[0] - TWELWTH) + rc[Rec::a_xy] * xLocal[0] * xLocal[1] +
           rc[Rec::a_xz] * xLocal[0] * xLocal[2],
       // Eq. (8) Balsara 2009
       rc[Rec::b_0] + rc[Rec::b_x] * xLocal[0] + rc[Rec::b_y] * xLocal[1] + rc[Rec::b_z] * xLocal[2] +
           rc[Rec::b_yy] * (xLocal[1] * xLocal[1] - TWELWTH) + rc[Rec::b_xy] * xLocal[0] * xLocal[1] +
           rc[Rec::b_yz] * xLocal[1] * xLocal[2],
       // Eq. (9) Balsara 2009
       rc[Rec::c_0] + rc[Rec::c_x] * xLocal[0] + rc[Rec::c_y] * xLocal[1] + rc[Rec::c_z] * xLocal[2] +
           rc[Rec::c_zz] * (xLocal[2] * xLocal[2] - TWELWTH) + rc[Rec::c_xz] * xLocal[0] * xLocal[2] +
           rc[Rec::c_yz] * xLocal[1] * xLocal[2],
   };
}

/*! Interpolate curl(perturbed B) to arbitrary x,y,z in cell
 *  Uses the reconstruction coefficients and equations from
 *  Divergence-free reconstruction of magnetic fields and WENO schemes for magnetohydrodynamics
 *  D.S. Balsara, J. Comp. Phys., 228, 2009
 *  doi:10.1016/j.jcp.2009.03.038
 *  and the wxMaxima file at
 *  doc/fieldsolver/Balsara_curlB_at_arbitrary_xyz.wxmx
 *
 * \param perb perturbed B fsGrid
 * \param dperb perturbed B derivatives fsGrid
 * \param technical technical fsGrid
 * \param fsgrid fsgfrid container
 * \param i local fsGrid x-index
 * \param j local fsGrid y-index
 * \param k local fsGrid z-index
 * \param x 3D global simulation x,y,z coordinates of point to interpolate to
 */
std::array<Real, 3> interpolateCurlB(
    fsgrids::perbspan perb,
    fsgrids::constdperbspan dperb,
    fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
    std::map<std::array<int, 3>, std::array<Real, Rec::N_REC_COEFFICIENTS>>& reconstructionCoefficientsCache, cint i,
    cint j, cint k, const std::array<Real, 3> x) {
   const auto stencil = fsgrid.makeStencil(i, j, k);

   cuint cellSysBoundaryFlag = technical[stencil.ooo()].sysBoundaryFlag;
   if (cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
      return {0, 0, 0};
   }

   // Balsara reconstruction formulas: x,y,z are in [-1/2, 1/2] local coordinates
   std::array<Real, 3> xLocal = fsgrid.physicalToCellFractional(x[0], x[1], x[2]);
   xLocal[0] -= 0.5;
   xLocal[1] -= 0.5;
   xLocal[2] -= 0.5;

   if (fabs(xLocal[0]) > 0.5 || fabs(xLocal[1]) > 0.5 || fabs(xLocal[2]) > 0.5) {
      cerr << __FILE__ << ":" << __LINE__ << ": Coordinate (" << xLocal[0] << "," << xLocal[1] << "," << xLocal[2]
           << ")  outside of this cell!" << endl;
      abort();
   }

   const std::array<Real, Rec::N_REC_COEFFICIENTS> rc = [&stencil, &reconstructionCoefficientsCache, &perb, &dperb]() {
      std::array<int, 3> cellIds = {stencil.i, stencil.j, stencil.k};
      // Actual use of the coefficient cache has proven not to be thread safe. But it appears to be reasonably fast even
      // without it.
      if (FieldTracing::fieldTracingParameters.useCache) {
#pragma omp critical
         {
            // Reconstruction order of the fields after Balsara 2009, 3 used to obtain 2nd
            // order curl(B) and allows for cache reuse, see interpolatePerturbedB above
            if (reconstructionCoefficientsCache.find(cellIds) == reconstructionCoefficientsCache.end()) {
               const auto rc = reconstructionCoefficients(perb, dperb, stencil, 3);
               reconstructionCoefficientsCache.insert({cellIds, rc});
            }
         }
         return reconstructionCoefficientsCache.at(cellIds);
      } else {
         return reconstructionCoefficients(perb, dperb, stencil, 3);
      }
   }();

   return {
       (12 * rc[Rec::c_yzz] * xLocal[2] * xLocal[2] + 24 * rc[Rec::c_yyz] * xLocal[1] * xLocal[2] -
        24 * rc[Rec::b_yzz] * xLocal[1] * xLocal[2] + 12 * rc[Rec::c_xyz] * xLocal[0] * xLocal[2] +
        12 * rc[Rec::c_yz] * xLocal[2] - 24 * rc[Rec::b_zz] * xLocal[2] - 12 * rc[Rec::b_yyz] * xLocal[1] * xLocal[1] -
        12 * rc[Rec::b_xyz] * xLocal[0] * xLocal[1] + 24 * rc[Rec::c_yy] * xLocal[1] - 12 * rc[Rec::b_yz] * xLocal[1] +
        12 * rc[Rec::c_xy] * xLocal[0] - 12 * rc[Rec::b_xz] * xLocal[0] - rc[Rec::c_yzz] + 12 * rc[Rec::c_y] -
        12 * rc[Rec::b_z] + rc[Rec::b_yyz]) /
           12,
       // See that minus if you ever copy again from wxMaxima!
       -(12 * rc[Rec::c_xzz] * xLocal[2] * xLocal[2] + 12 * rc[Rec::c_xyz] * xLocal[1] * xLocal[2] +
         24 * rc[Rec::c_xxz] * xLocal[0] * xLocal[2] - 24 * rc[Rec::a_xzz] * xLocal[0] * xLocal[2] +
         12 * rc[Rec::c_xz] * xLocal[2] - 24 * rc[Rec::a_zz] * xLocal[2] - 12 * rc[Rec::a_xyz] * xLocal[0] * xLocal[1] +
         12 * rc[Rec::c_xy] * xLocal[1] - 12 * rc[Rec::a_yz] * xLocal[1] - 12 * rc[Rec::a_xxz] * xLocal[0] * xLocal[0] +
         24 * rc[Rec::c_xx] * xLocal[0] - 12 * rc[Rec::a_xz] * xLocal[0] - rc[Rec::c_xzz] + 12 * rc[Rec::c_x] -
         12 * rc[Rec::a_z] + rc[Rec::a_xxz]) /
           12,
       (12 * rc[Rec::b_xyz] * xLocal[1] * xLocal[2] - 12 * rc[Rec::a_xyz] * xLocal[0] * xLocal[2] +
        12 * rc[Rec::b_xz] * xLocal[2] - 12 * rc[Rec::a_yz] * xLocal[2] + 12 * rc[Rec::b_xyy] * xLocal[1] * xLocal[1] +
        24 * rc[Rec::b_xxy] * xLocal[0] * xLocal[1] - 24 * rc[Rec::a_xyy] * xLocal[0] * xLocal[1] +
        12 * rc[Rec::b_xy] * xLocal[1] - 24 * rc[Rec::a_yy] * xLocal[1] - 12 * rc[Rec::a_xxy] * xLocal[0] * xLocal[0] +
        24 * rc[Rec::b_xx] * xLocal[0] - 12 * rc[Rec::a_xy] * xLocal[0] - rc[Rec::b_xyy] + 12 * rc[Rec::b_x] -
        12 * rc[Rec::a_y] + rc[Rec::a_xxy]) /
           12,
   };
}
