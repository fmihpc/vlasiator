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

#include <cstdlib>

#include "derivatives.hpp"
#include "fs_common.h"
#include "fs_limiters.h"
#include <Eigen/Geometry>
#include <span>

template <typename T, size_t N> struct DerivativesData {
   const std::array<T, N>& ooo = {};
   const std::array<T, N>& poo = {};
   const std::array<T, N>& moo = {};
   const std::array<T, N>& opo = {};
   const std::array<T, N>& omo = {};
   const std::array<T, N>& oop = {};
   const std::array<T, N>& oom = {};
};

void computeMomentsDerivatives(fsgrids::momentsspan moments,
                    std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                    const fsgrid::FsStencil& stencil, const bool atSysBoundary) {
   using dmo = fsgrids::dmoments;
   using mom = fsgrids::moments;

   std::array<Real, dmo::N_DMOMENTS>& dMoments = dmoments[stencil.ooo()];

   auto computeDiff = [](const auto& i, const auto& right, const auto& left) { return 0.5 * (right[i] - left[i]); };

   auto computeLimiter = [](const auto& i, const auto& right, const auto& left, const auto& center) {
      return limiter(left[i], center[i], right[i]);
   };

   // Constants for electron pressure derivatives, see also ldz_gradpe.cpp
   // Calculate anchor point constants: First the pressure, then a derived constant.
   const Real Pe_anchor = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   const Real Pe_const = Pe_anchor * pow(Parameters::electronDensity, -Parameters::electronPTindex);
   auto computeGradPeLimiter = [=](const auto& right, const auto& left, const auto& center) {
      // pres_e = const * np.power(rho_e, index)
      return Pe_const * limiter(pow(left[mom::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex),
                               pow(center[mom::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex),
                               pow(right[mom::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex));
   };

   auto computeGradPeDiff = [=](const auto& right, const auto& left) {
      return Pe_const * 0.5 * (pow(right[mom::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex) - pow(left[mom::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex));
   };

   const DerivativesData momData{
       moments[stencil.ooo()], moments[stencil.poo()], moments[stencil.moo()], moments[stencil.opo()],
       moments[stencil.omo()],   moments[stencil.oop()],  moments[stencil.oom()],
   };

   {
#ifdef DEBUG_SOLVERS
      const auto& cv = momData.ooo[mom::RHOM];
      if (cv <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << (cv < 0 ? " Negative" : " Zero") << " density in fsgrid cell id " << stencil.indexFromOffset(0,0,0) << std::endl;
         abort();
      }

      const auto& lv = momData.moo[mom::RHOM];
      if (lv <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << (lv < 0 ? " Negative" : " Zero") << " density in fsgrid cell " << stencil.indexFromOffset(-1,0,0)
                   << std::endl;
         abort();
      }

      const auto& rv = momData.poo[mom::RHOM];
      if (rv <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << (rv < 0 ? " Negative" : " Zero") << " density in fsgrid cell " << stencil.indexFromOffset(1,0,0)
                   << std::endl;
         abort();
      }
#endif
   }

   static constexpr std::array moms{mom::RHOM, mom::RHOQ, mom::P_11, mom::P_22, mom::P_33, mom::VX, mom::VY, mom::VZ};
   static constexpr std::array dmix{
       dmo::drhomdx, dmo::drhoqdx, dmo::dp11dx, dmo::dp22dx, dmo::dp33dx, dmo::dVxdx, dmo::dVydx, dmo::dVzdx,
   };
   static constexpr std::array dmiy{
       dmo::drhomdy, dmo::drhoqdy, dmo::dp11dy, dmo::dp22dy, dmo::dp33dy, dmo::dVxdy, dmo::dVydy, dmo::dVzdy,
   };
   static constexpr std::array dmiz{
       dmo::drhomdz, dmo::drhoqdz, dmo::dp11dz, dmo::dp22dz, dmo::dp33dz, dmo::dVxdz, dmo::dVydz, dmo::dVzdz,
   };

   if (atSysBoundary) {
      for (size_t i = 0; i < moms.size(); i++) {
         dMoments[dmix[i]] = computeDiff(moms[i], momData.poo, momData.moo);
         dMoments[dmiy[i]] = computeDiff(moms[i], momData.opo, momData.omo);
         dMoments[dmiz[i]] = computeDiff(moms[i], momData.oop, momData.oom);
      }
      dMoments[dmo::dPedx] = computeGradPeDiff(momData.poo, momData.moo);
      dMoments[dmo::dPedy] = computeGradPeDiff(momData.opo, momData.omo);
      dMoments[dmo::dPedz] = computeGradPeDiff(momData.oop, momData.oom);
   } else {
      for (size_t i = 0; i < moms.size(); i++) {
         dMoments[dmix[i]] = computeLimiter(moms[i], momData.poo, momData.moo, momData.ooo);
         dMoments[dmiy[i]] = computeLimiter(moms[i], momData.opo, momData.omo, momData.ooo);
         dMoments[dmiz[i]] = computeLimiter(moms[i], momData.oop, momData.oom, momData.ooo);
      }
      dMoments[dmo::dPedx] = computeGradPeLimiter(momData.poo, momData.moo, momData.ooo);
      dMoments[dmo::dPedy] = computeGradPeLimiter(momData.opo, momData.omo, momData.ooo);
      dMoments[dmo::dPedz] = computeGradPeLimiter(momData.oop, momData.oom, momData.ooo);
   }


}

void computePerbDerivatives(std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                 std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb, const fsgrid::FsStencil& stencil,
                 bool dontCompute2ndDerivatives, bool atSysBoundary, cuint sysBoundaryFlag) {
   using dpb = fsgrids::dperb;
   using bfi = fsgrids::bfield;
   std::array<Real, dpb::N_DPERB>& dPerB = dperb[stencil.ooo()];

   auto computeDiff = [](const auto& i, const auto& right, const auto& left) { return 0.5 * (right[i] - left[i]); };

   auto computeLimiter = [](const auto& i, const auto& right, const auto& left, const auto& center) {
      return limiter(left[i], center[i], right[i]);
   };
   const DerivativesData perbData{
       perb[stencil.ooo()], perb[stencil.poo()], perb[stencil.moo()], perb[stencil.opo()],
       perb[stencil.omo()],   perb[stencil.oop()],  perb[stencil.oom()],
   };

   if (atSysBoundary) {
      dPerB[dpb::dPERBydx] = computeDiff(bfi::PERBY, perbData.poo, perbData.moo);
      dPerB[dpb::dPERBzdx] = computeDiff(bfi::PERBZ, perbData.poo, perbData.moo);
      dPerB[dpb::dPERBxdy] = computeDiff(bfi::PERBX, perbData.opo, perbData.omo);
      dPerB[dpb::dPERBzdy] = computeDiff(bfi::PERBZ, perbData.opo, perbData.omo);
      dPerB[dpb::dPERBxdz] = computeDiff(bfi::PERBX, perbData.oop, perbData.oom);
      dPerB[dpb::dPERBydz] = computeDiff(bfi::PERBY, perbData.oop, perbData.oom);
   } else {
      dPerB[dpb::dPERBydx] = computeLimiter(bfi::PERBY, perbData.poo, perbData.moo, perbData.ooo);
      dPerB[dpb::dPERBzdx] = computeLimiter(bfi::PERBZ, perbData.poo, perbData.moo, perbData.ooo);
      dPerB[dpb::dPERBxdy] = computeLimiter(bfi::PERBX, perbData.opo, perbData.omo, perbData.ooo);
      dPerB[dpb::dPERBzdy] = computeLimiter(bfi::PERBZ, perbData.opo, perbData.omo, perbData.ooo);
      dPerB[dpb::dPERBxdz] = computeLimiter(bfi::PERBX, perbData.oop, perbData.oom, perbData.ooo);
      dPerB[dpb::dPERBydz] = computeLimiter(bfi::PERBY, perbData.oop, perbData.oom, perbData.ooo);
   }

   if (dontCompute2ndDerivatives) {
      dPerB[dpb::dPERBydxx] = 0.0;
      dPerB[dpb::dPERBzdxx] = 0.0;
      dPerB[dpb::dPERBxdyy] = 0.0;
      dPerB[dpb::dPERBzdyy] = 0.0;
      dPerB[dpb::dPERBxdzz] = 0.0;
      dPerB[dpb::dPERBydzz] = 0.0;
      dPerB[dpb::dPERBxdyz] = 0.0;
      dPerB[dpb::dPERBydxz] = 0.0;
      dPerB[dpb::dPERBzdxy] = 0.0;
   } else {
      auto compute2ndDerivative = [](auto i, const auto& right, const auto& left, const auto& center) {
         return left[i] + right[i] - 2.0 * center[i];
      };
      dPerB[dpb::dPERBydxx] = compute2ndDerivative(bfi::PERBY, perbData.poo, perbData.moo, perbData.ooo);
      dPerB[dpb::dPERBzdxx] = compute2ndDerivative(bfi::PERBZ, perbData.poo, perbData.moo, perbData.ooo);
      dPerB[dpb::dPERBxdyy] = compute2ndDerivative(bfi::PERBX, perbData.opo, perbData.omo, perbData.ooo);
      dPerB[dpb::dPERBzdyy] = compute2ndDerivative(bfi::PERBZ, perbData.opo, perbData.omo, perbData.ooo);
      dPerB[dpb::dPERBxdzz] = compute2ndDerivative(bfi::PERBX, perbData.oop, perbData.oom, perbData.ooo);
      dPerB[dpb::dPERBydzz] = compute2ndDerivative(bfi::PERBY, perbData.oop, perbData.oom, perbData.ooo);

      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         auto compute = [&perb](auto bl, auto br, auto tl, auto tr, auto i) {
            const auto& botLeft = perb[bl];
            const auto& botRght = perb[br];
            const auto& topLeft = perb[tl];
            const auto& topRght = perb[tr];
            return FOURTH * (botLeft[i] + topRght[i] - botRght[i] - topLeft[i]);
         };

         dPerB[dpb::dPERBxdyz] =
             compute(stencil.omm(), stencil.opm(), stencil.omp(), stencil.opp(), bfi::PERBX);
         dPerB[dpb::dPERBydxz] =
             compute(stencil.mom(), stencil.pom(), stencil.mop(), stencil.pop(), bfi::PERBY);
         dPerB[dpb::dPERBzdxy] =
             compute(stencil.mmo(), stencil.pmo(), stencil.mpo(), stencil.ppo(), bfi::PERBZ);
      }
   }
}

/*! \brief Low-level spatial derivatives calculation.
 *
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in
 * project.h. Uses RHO, V[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the
 * second-order method, and RHO_DT2, V[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param fsgrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * \param doMoments Bool telling whether the derivatives for moments need updating too.
 *
 * \sa calculateDerivativesSimple calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivatives(std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                          fsgrids::momentsspan moments,
                          std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb,
                          std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                          const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer,
                          const bool doMoments) {
   /*
    * For sysBoundaryLayer 1 or 2, we are near a boundary, and we wish to use regular centered differences instead of
    * slope limiter-adjusted values. This is to minimize oscillations as a smooth behaviour is required near artificial
    * boundaries, unlike at boundaries and shocks inside the simulation domain.
    */
   const bool atSysBoundary =
       sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY);
   const bool dontCompute2ndDerivatives = Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1;

   if (doMoments) {
      computeMomentsDerivatives(moments, dmoments, stencil, atSysBoundary);
   }

   computePerbDerivatives(perb, dperb, stencil, dontCompute2ndDerivatives, atSysBoundary, sysBoundaryFlag);

   if (sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
      SBC::SysBoundaryCondition::setCellDerivativesToZero(dperb, dmoments, stencil, 3);
      SBC::SysBoundaryCondition::setCellDerivativesToZero(dperb, dmoments, stencil, 4);
      SBC::SysBoundaryCondition::setCellDerivativesToZero(dperb, dmoments, stencil, 5);
   }
}

/*! \brief High-level derivative calculation wrapper function.
 *

 * B has to be updated because after the system boundary update in propagateMagneticFieldSimple there is no consistent
 state of B yet everywhere.
 *
 * Then the derivatives are calculated.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param fsgrid fsGrid holding technical information (such as boundary types)
 * \param communicateMoments If true, the derivatives of moments (rho, V, P) are communicated to neighbours.

 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivativesSimple(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                std::span<std::array<Real, fsgrids::moments::N_MOMENTS>> moments,
                                std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb,
                                std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                                std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid,
                                const bool doMoments) {
   phiprof::Timer derivativesTimer{"Calculate face derivatives"};
   const size_t numCells = fsgrid.getNumCells();

   phiprof::Timer mpiTimer{"FS derivatives ghost updates MPI", {"MPI"}};
   fsgrid.updateGhostCells(perb);
   if (doMoments) {
      fsgrid.updateGhostCells(moments);
   }
   mpiTimer.stop();

   // Calculate derivatives
   fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                       phiprof::initializeTimer("FS derivatives compute cells"), technical,
                       [=](const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
                          calculateDerivatives(perb, moments, dperb, dmoments, stencil, sysBoundaryFlag,
                                               sysBoundaryLayer, doMoments);
                       });

   derivativesTimer.stop(numCells, "Spatial Cells");
}

/*! \brief Low-level spatial derivatives calculation.
 *
 * Calculate the spatial derivatives of BVOL or set them to zero.
 *
 * For sysBoundaryLayer 1 or 2, we are near a boundary, and we wish to use regular centered differences instead of slope
 * limiter-adjusted values. This is to minimize oscillations as a smooth behaviour is required near artificial
 * boundaries, unlike at boundaries and shocks inside the simulation domain.
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param fsgrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateBVOLDerivatives(std::span<std::array<Real, fsgrids::volfields::N_VOL>> vol,
                              std::span<const fsgrids::technical> technical, const fsgrid::FsStencil& stencil) {
   const auto& tech = technical[stencil.ooo()];

   cuint sysBoundaryFlag = tech.sysBoundaryFlag;
   cuint sysBoundaryLayer = tech.sysBoundaryLayer;
   const bool notSysBoundary =
       sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY);

   using vf = fsgrids::volfields;

   auto computeDiff = [](const auto& i, const auto& right, const auto& left) { return 0.5 * (right[i] - left[i]); };

   auto computeLimiter = [](const auto& i, const auto& right, const auto& left, const auto& center) {
      return limiter(left[i], center[i], right[i]);
   };

   auto& volCenter = vol[stencil.ooo()];
   const auto poo = stencil.poo();
   const auto moo = stencil.moo();
   if (notSysBoundary) {
      volCenter[vf::dPERBXVOLdx] = computeDiff(vf::PERBXVOL, vol[poo], vol[moo]);
      volCenter[vf::dPERBYVOLdx] = computeDiff(vf::PERBYVOL, vol[poo], vol[moo]);
      volCenter[vf::dPERBZVOLdx] = computeDiff(vf::PERBZVOL, vol[poo], vol[moo]);
   } else {
      volCenter[vf::dPERBXVOLdx] = computeLimiter(vf::PERBXVOL, vol[poo], vol[moo], volCenter);
      volCenter[vf::dPERBYVOLdx] = computeLimiter(vf::PERBYVOL, vol[poo], vol[moo], volCenter);
      volCenter[vf::dPERBZVOLdx] = computeLimiter(vf::PERBZVOL, vol[poo], vol[moo], volCenter);
   }

   const auto opo = stencil.opo();
   const auto omo = stencil.omo();
   if (notSysBoundary) {
      volCenter[vf::dPERBXVOLdy] = computeDiff(vf::PERBXVOL, vol[opo], vol[omo]);
      volCenter[vf::dPERBYVOLdy] = computeDiff(vf::PERBYVOL, vol[opo], vol[omo]);
      volCenter[vf::dPERBZVOLdy] = computeDiff(vf::PERBZVOL, vol[opo], vol[omo]);
   } else {
      volCenter[vf::dPERBXVOLdy] = computeLimiter(vf::PERBXVOL, vol[opo], vol[omo], volCenter);
      volCenter[vf::dPERBYVOLdy] = computeLimiter(vf::PERBYVOL, vol[opo], vol[omo], volCenter);
      volCenter[vf::dPERBZVOLdy] = computeLimiter(vf::PERBZVOL, vol[opo], vol[omo], volCenter);
   }

   const auto oop = stencil.oop();
   const auto oom = stencil.oom();
   if (notSysBoundary) {
      volCenter[vf::dPERBXVOLdz] = computeDiff(vf::PERBXVOL, vol[oop], vol[oom]);
      volCenter[vf::dPERBYVOLdz] = computeDiff(vf::PERBYVOL, vol[oop], vol[oom]);
      volCenter[vf::dPERBZVOLdz] = computeDiff(vf::PERBZVOL, vol[oop], vol[oom]);
   } else {
      volCenter[vf::dPERBXVOLdz] = computeLimiter(vf::PERBXVOL, vol[oop], vol[oom], volCenter);
      volCenter[vf::dPERBYVOLdz] = computeLimiter(vf::PERBYVOL, vol[oop], vol[oom], volCenter);
      volCenter[vf::dPERBZVOLdz] = computeLimiter(vf::PERBZVOL, vol[oop], vol[oom], volCenter);
   }
}

/*! \brief High-level derivative calculation wrapper function.
 *
 * BVOL has been calculated locally by calculateVolumeAveragedFields but not communicated.
 * For the acceleration step one needs the cross-derivatives of BVOL
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param fsgrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateBVOLDerivativesSimple(std::span<std::array<Real, fsgrids::volfields::N_VOL>> vol,
                                    std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid) {
   phiprof::Timer derivsTimer{"Calculate volume derivatives"};
   const size_t numCells = fsgrid.getNumCells();

   phiprof::Timer commTimer{"BVOL derivatives ghost updates MPI", {"MPI"}};
   fsgrid.updateGhostCells(vol);
   commTimer.stop(numCells, "Spatial Cells");

   // Calculate derivatives
   fsgrid.parallel_for([](int timerId) ->  phiprof::Timer { return phiprof::Timer{timerId}; },
                       phiprof::initializeTimer("FS derivatives BVOL compute cells"), technical,
                       [=](const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
                          calculateBVOLDerivatives(vol, technical, stencil);
                      });

   derivsTimer.stop(numCells, "Spatial Cells");
}

/*! \brief Low-level curvature calculation.
 *
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param bgbGrid fsGrid holding the background fields
 * \param fsgrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 *
 * http://fusionwiki.ciemat.es/wiki/Magnetic_curvature
 *
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateCurvature(std::span<std::array<Real, fsgrids::volfields::N_VOL>> vol,
                        std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                        const fsgrid::FsStencil& stencil, const std::array<Real, 3>& gridSpacing) {
   auto compute = [&vol, &bgb](auto i) -> std::array<Real, 3> {
      const auto& b = bgb[i];
      const auto& v = vol[i];
      const Real bx = b[fsgrids::bgbfield::BGBXVOL] + v[fsgrids::volfields::PERBXVOL];
      const Real by = b[fsgrids::bgbfield::BGBYVOL] + v[fsgrids::volfields::PERBYVOL];
      const Real bz = b[fsgrids::bgbfield::BGBZVOL] + v[fsgrids::volfields::PERBZVOL];
      const Real bnorm = sqrt(bx * bx + by * by + bz * bz);

      return {
          bx / bnorm,
          by / bnorm,
          bz / bnorm,
      };
   };

   const auto [bx, by, bz] = compute(stencil.ooo());
   const auto [left_x_bx, left_x_by, left_x_bz] = compute(stencil.moo());
   const auto [rght_x_bx, rght_x_by, rght_x_bz] = compute(stencil.poo());
   const auto [left_y_bx, left_y_by, left_y_bz] = compute(stencil.omo());
   const auto [rght_y_bx, rght_y_by, rght_y_bz] = compute(stencil.opo());
   const auto [left_z_bx, left_z_by, left_z_bz] = compute(stencil.oom());
   const auto [rght_z_bx, rght_z_by, rght_z_bz] = compute(stencil.oop());

   auto& volCenter = vol[stencil.ooo()];
   volCenter[fsgrids::volfields::CURVATUREX] = bx * 0.5 * (rght_x_bx - left_x_bx) / gridSpacing[0] +
                                               by * 0.5 * (rght_y_bx - left_y_bx) / gridSpacing[1] +
                                               bz * 0.5 * (rght_z_bx - left_z_bx) / gridSpacing[2];
   volCenter[fsgrids::volfields::CURVATUREY] = bx * 0.5 * (rght_x_by - left_x_by) / gridSpacing[0] +
                                               by * 0.5 * (rght_y_by - left_y_by) / gridSpacing[1] +
                                               bz * 0.5 * (rght_z_by - left_z_by) / gridSpacing[2];
   volCenter[fsgrids::volfields::CURVATUREZ] = bx * 0.5 * (rght_x_bz - left_x_bz) / gridSpacing[0] +
                                               by * 0.5 * (rght_y_bz - left_y_bz) / gridSpacing[1] +
                                               bz * 0.5 * (rght_z_bz - left_z_bz) / gridSpacing[2];
}

/*! \brief High-level curvature calculation wrapper function.
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param bgbGrid fsGrid holding the background fields
 * \param fsgrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateCurvatureSimple(std::span<std::array<Real, fsgrids::volfields::N_VOL>> vol,
                              std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                              std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid) {
   phiprof::Timer curvatureTimer{"Calculate curvature"};
   const auto& gridSpacing = fsgrid.getGridSpacing();
   const size_t numCells = fsgrid.getNumCells();

   phiprof::Timer commTimer{"Calculate curvature ghost updates MPI", {"MPI"}};
   fsgrid.updateGhostCells(vol);
   commTimer.stop(numCells, "Spatial Cells");

   fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                       phiprof::initializeTimer("Calculate curvature compute cells"), technical,
                       [=](const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
                          const bool compute = (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY &&
                                     sysBoundaryLayer != 1 && sysBoundaryLayer != 2);
                          if (compute) {
                             calculateCurvature(vol, bgb, stencil, gridSpacing);
                          }
                       });
   curvatureTimer.stop(numCells, "Spatial Cells");
}

/*! \brief Returns perturbed volumetric B of cell
 *
 */
[[maybe_unused]] static std::array<Real, 3> getPerBVol(SpatialCell* cell) {
   return std::array<Real, 3>{{cell->parameters[CellParams::PERBXVOL], cell->parameters[CellParams::PERBYVOL],
                               cell->parameters[CellParams::PERBZVOL]}};
}

/*! \brief Returns volumetric B of cell
 *
 */
static std::array<Real, 3> getBVol(SpatialCell* cell) {
   return std::array<Real, 3>{{cell->parameters[CellParams::BGBXVOL] + cell->parameters[CellParams::PERBXVOL],
                               cell->parameters[CellParams::BGBYVOL] + cell->parameters[CellParams::PERBYVOL],
                               cell->parameters[CellParams::BGBZVOL] + cell->parameters[CellParams::PERBZVOL]}};
}

/*! \brief Calculates momentum density of cell
 *
 */
static std::array<Real, 3> getMomentumDensity(SpatialCell* cell) {
   Real rho = cell->parameters[CellParams::RHOM];
   return std::array<Real, 3>{{rho * cell->parameters[CellParams::VX], rho * cell->parameters[CellParams::VY],
                               rho * cell->parameters[CellParams::VZ]}};
}

/*! \brief Calculates energy density for spatial cell
 *
 */
static Real calculateU(SpatialCell* cell) {
   Real rho = cell->parameters[CellParams::RHOM];
   std::array<Real, 3> p = getMomentumDensity(cell);
   std::array<Real, 3> B = getBVol(cell);
   return (pow(B[0], 2) + pow(B[1], 2) + pow(B[2], 2)) / (2.0 * physicalconstants::MU_0) + // Magnetic field energy
          (rho > EPS ? (pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2)) / (2.0 * cell->parameters[CellParams::RHOM])
                     : 0.0); // Kinetic energy
}

/*! \brief Calculates pressure anistotropy from B and P
 *  \param P elements of pressure order in order: P_11, P_22, P_33, P_23, P_13, P_12
 */
static Real calculateAnisotropy(const Eigen::Matrix3d& rot, const std::array<Real, 6>& P) {
   // Now, rotation matrix to get parallel and perpendicular pressure
   // Eigen::Quaterniond q {Quaterniond::FromTwoVectors(Eigen::vector3d{0, 0, 1}, Eigen::vector3d{myB[0], myB[1],
   // myB[2]})}; Eigen::Matrix3d rot = q.toRotationMatrix();
   Eigen::Matrix3d Ptensor{
       {P[0], P[5], P[4]},
       {P[5], P[1], P[3]},
       {P[4], P[3], P[2]},
   };

   Eigen::Matrix3d transposerot = rot.transpose();
   Eigen::Matrix3d Pprime = rot * Ptensor * transposerot;

   Real Panisotropy{0.0};
   if (Pprime(2, 2) > EPS) {
      Panisotropy = (Pprime(0, 0) + Pprime(1, 1)) / (2 * Pprime(2, 2));
   }

   return Panisotropy;
}

/*! \brief Low-level scaled gradients calculation
 *
 * For the SpatialCell* cell and its neighbors, calculate scaled gradients and their maximum alpha
 * The gradients are the same as in the GUMICS simulation, see
 * Janhunen, P., Palmroth, M., Laitinen, T., Honkonen, I., Juusola, L., Facsko, G., & Pulkkinen, T. I. (2012). The
 * GUMICS-4 global MHD magnetosphere-ionosphere coupling simulation. Journal of Atmospheric and Solar - Terrestrial
 * Physics, 80, 48-59. https://doi.org/10.1016/j.jastp.2012.03.006
 *
 */
void calculateScaledDeltas(SpatialCell* cell, std::vector<SpatialCell*>& neighbors) {
   Real dRho{0};
   Real dU{0};
   Real dPsq{0};
   Real dBsq{0};
   Real dB{0};

   Real myRho{cell->parameters[CellParams::RHOM]};
   Real myU{calculateU(cell)};
   Real myV{std::sqrt(std::pow(cell->parameters[CellParams::VX], 2) + std::pow(cell->parameters[CellParams::VY], 2) +
                      std::pow(cell->parameters[CellParams::VZ], 2))};
   Real maxV{myV};
   std::array<Real, 3> myP = getMomentumDensity(cell);
   std::array<Real, 3> myB = getBVol(cell);
   for (SpatialCell* neighbor : neighbors) {
      Real otherRho = neighbor->parameters[CellParams::RHOM];
      Real otherU = calculateU(neighbor);
      Real otherV{std::sqrt(std::pow(neighbor->parameters[CellParams::VX], 2) +
                            std::pow(neighbor->parameters[CellParams::VY], 2) +
                            std::pow(neighbor->parameters[CellParams::VZ], 2))};
      std::array<Real, 3> otherP = getMomentumDensity(neighbor);
      std::array<Real, 3> otherB = getBVol(neighbor);
      Real deltaBsq = pow(myB[0] - otherB[0], 2) + pow(myB[1] - otherB[1], 2) + pow(myB[2] - otherB[2], 2);

      if (myV < EPS) {
         maxV = std::max(maxV, otherV);
      }

      Real maxRho = std::max(myRho, otherRho);
      if (maxRho > EPS) {
         dRho = std::max(fabs(myRho - otherRho) / maxRho, dRho);
      }
      Real maxU = std::max(myU, otherU);
      if (maxU > EPS) {
         dU = std::max(fabs(myU - otherU) / maxU, dU);
         dBsq = std::max(deltaBsq / (2 * physicalconstants::MU_0 * maxU), dBsq);
         if (myRho > EPS) {
            dPsq = std::max((pow(myP[0] - otherP[0], 2) + pow(myP[1] - otherP[1], 2) + pow(myP[2] - otherP[2], 2)) /
                                (2 * myRho * maxU),
                            dPsq);
         }
      }
      Real maxB = sqrt(std::max(pow(myB[0], 2) + pow(myB[1], 2) + pow(myB[2], 2),
                                pow(otherB[0], 2) + pow(otherB[1], 2) + pow(otherB[2], 2)));
      if (maxB > EPS) {
         dB = std::max(sqrt(deltaBsq) / maxB, dB);
      }
   }

   Real alpha{0.0};
   alpha = std::max(alpha, dRho * P::alphaDRhoWeight);
   alpha = std::max(alpha, dU * P::alphaDUWeight);
   alpha = std::max(alpha, dPsq * P::alphaDPSqWeight);
   alpha = std::max(alpha, dBsq * P::alphaDBSqWeight);
   alpha = std::max(alpha, dB * P::alphaDBWeight);

   Real dBXdy{cell->derivativesBVOL[bvolderivatives::dPERBXVOLdy]};
   Real dBXdz{cell->derivativesBVOL[bvolderivatives::dPERBXVOLdz]};
   Real dBYdx{cell->derivativesBVOL[bvolderivatives::dPERBYVOLdx]};
   Real dBYdz{cell->derivativesBVOL[bvolderivatives::dPERBYVOLdz]};
   Real dBZdx{cell->derivativesBVOL[bvolderivatives::dPERBZVOLdx]};
   Real dBZdy{cell->derivativesBVOL[bvolderivatives::dPERBZVOLdy]};

   // Note missing factor of mu_0, since we want B and J in same units later
   myB = getBVol(cell); // Redundant, but this makes sure we use total B here
   std::array<Real, 3> myJ = {dBZdy - dBYdz, dBXdz - dBZdx, dBYdx - dBXdy};
   Real BdotJ{0.0};
   Real Bsq{0.0};
   Real J{0.0};
   for (int i = 0; i < 3; ++i) {
      BdotJ += myB[i] * myJ[i];
      Bsq += myB[i] * myB[i];
      J += myJ[i] * myJ[i];
   }
   J = std::sqrt(J);

   Real Bperp{0.0};
   if (Bsq > EPS) {
      for (int i = 0; i < 3; ++i) {
         Bperp += std::pow(myB[i] * (1 - BdotJ / Bsq), 2);
      }
      Bperp = std::sqrt(Bperp);
   }

   // Vorticity
   Real dVxdy{cell->derivativesV[vderivatives::dVxdy]};
   Real dVxdz{cell->derivativesV[vderivatives::dVxdz]};
   Real dVydx{cell->derivativesV[vderivatives::dVydx]};
   Real dVydz{cell->derivativesV[vderivatives::dVydz]};
   Real dVzdx{cell->derivativesV[vderivatives::dVzdx]};
   Real dVzdy{cell->derivativesV[vderivatives::dVzdy]};
   Real vorticity{std::sqrt(std::pow(dVxdy - dVydz, 2) + std::pow(dVxdz - dVzdx, 2) + std::pow(dVydx - dVxdy, 2))};
   // Real vA {std::sqrt(Bsq / (physicalconstants::MU_0 * myRho))};
   Real amr_vorticity{-1.0}; // Error value
   if (maxV > EPS) {
      amr_vorticity = vorticity * cell->parameters[CellParams::DX] / maxV;
   }

   std::array<Real, 6> myPressure{cell->parameters[CellParams::P_11], cell->parameters[CellParams::P_22],
                                  cell->parameters[CellParams::P_33], cell->parameters[CellParams::P_23],
                                  cell->parameters[CellParams::P_13], cell->parameters[CellParams::P_12]};
   Eigen::Matrix3d rot =
       Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d{myB[0], myB[1], myB[2]}, Eigen::Vector3d{0, 0, 1})
           .normalized()
           .toRotationMatrix();
   Real Panisotropy{calculateAnisotropy(rot, myPressure)};
   for (const auto& pop : cell->get_populations()) {
      // TODO I hate this. Change all this crap to std::vectors?
      std::array<Real, 6> popP{pop.P[0], pop.P[1], pop.P[2], pop.P[3], pop.P[4], pop.P[5]};
      Real popPanisotropy{calculateAnisotropy(rot, popP)};
      // low value refines
      Panisotropy = std::min(Panisotropy, popPanisotropy);
   }

   cell->parameters[CellParams::AMR_DRHO] = dRho;
   cell->parameters[CellParams::AMR_DU] = dU;
   cell->parameters[CellParams::AMR_DPSQ] = dPsq;
   cell->parameters[CellParams::AMR_DBSQ] = dBsq;
   cell->parameters[CellParams::AMR_DB] = dB;
   cell->parameters[CellParams::AMR_ALPHA1] = alpha;
   cell->parameters[CellParams::AMR_ALPHA2] =
       cell->parameters[CellParams::DX] * J / (Bperp + EPS); // Epsilon in denominator so we don't get infinities
   cell->parameters[CellParams::P_ANISOTROPY] = Panisotropy;
   // Experimental, current scaling is bulk velocity
   cell->parameters[CellParams::AMR_VORTICITY] = amr_vorticity;
}

/*! \brief High-level scaled gradient calculation wrapper function.
 *
 * Calculates gradients needed for alpha everywhere in the grid
 *
 */

void calculateScaledDeltasSimple(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
   const vector<CellID>& cells = getLocalCells();
   int N_cells = cells.size();
   phiprof::Timer gradientsTimer{"Calculate volume gradients"};
   int computeTimerId{phiprof::initializeTimer("Calculate volume gradients compute cells")};

   phiprof::Timer commTimer{"Calculate volume gradients ghost updates MPI", {"MPI"}};
   // We only need nearest neighbourhood and spatial data here
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);
   mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::NEAREST);
   commTimer.stop(N_cells, "Spatial Cells");

// Calculate derivatives
#pragma omp parallel
   {
      phiprof::Timer computeTimer{computeTimerId};
#pragma omp for
      for (uint i = 0; i < cells.size(); ++i) {
         CellID id = cells[i];
         SpatialCell* cell = mpiGrid[id];
         std::vector<SpatialCell*> neighbors;
         for (const auto& [neighbor, dir] : mpiGrid.get_face_neighbors_of(id)) {
            neighbors.push_back(mpiGrid[neighbor]);
         }
         calculateScaledDeltas(cell, neighbors);
      }
      computeTimer.stop(N_cells, "Spatial Cells");
   }

   gradientsTimer.stop(N_cells, "Spatial Cells");
}
