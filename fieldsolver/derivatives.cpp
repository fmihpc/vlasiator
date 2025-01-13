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
   const std::array<T, N>& center = {};
   const std::array<T, N>& right = {};
   const std::array<T, N>& left = {};
   const std::array<T, N>& up = {};
   const std::array<T, N>& down = {};
   const std::array<T, N>& near = {};
   const std::array<T, N>& far = {};
};

void computeMoments(std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments,
                    std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                    const fsgrid::FsStencil& stencil, const bool notSysBoundary) {
   using dmo = fsgrids::dmoments;
   using mom = fsgrids::moments;

   std::array<Real, dmo::N_DMOMENTS>& dMoments = dmoments[stencil.center()];

   auto computeDiff = [](const auto& i, const auto& right, const auto& left) { return 0.5 * (right[i] - left[i]); };

   auto computeLimiter = [](const auto& i, const auto& right, const auto& left, const auto& center) {
      return limiter(left[i], center[i], right[i]);
   };

   const DerivativesData momData{
       moments[stencil.center()], moments[stencil.right()], moments[stencil.left()], moments[stencil.up()],
       moments[stencil.down()],   moments[stencil.near()],  moments[stencil.far()],
   };

   {
#ifdef DEBUG_SOLVERS
      const auto& cv = momData.center[mom::RHOM];
      if (cv <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << (cv < 0 ? " Negative" : " Zero") << " density in spatial cell at ("
                   << i << " " << j << " " << k << ")" << std::endl;
         abort();
      }

      const auto& lv = momData.left[mom::RHOM];
      if (lv <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << (lv < 0 ? " Negative" : " Zero") << " density in spatial cell"
                   << std::endl;
         abort();
      }

      const auto& rv = momData.right[mom::RHOM];
      if (rv <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << (rv < 0 ? " Negative" : " Zero") << " density in spatial cell"
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

   if (notSysBoundary) {
      for (size_t i = 0; i < moms.size(); i++) {
         dMoments[dmix[i]] = computeDiff(moms[i], momData.right, momData.left);
      }
      for (size_t i = 0; i < moms.size(); i++) {
         dMoments[dmiy[i]] = computeDiff(moms[i], momData.up, momData.down);
      }
      for (size_t i = 0; i < moms.size(); i++) {
         dMoments[dmiz[i]] = computeDiff(moms[i], momData.near, momData.far);
      }
   } else {
      for (size_t i = 0; i < moms.size(); i++) {
         dMoments[dmix[i]] = computeLimiter(moms[i], momData.right, momData.left, momData.center);
      }
      for (size_t i = 0; i < moms.size(); i++) {
         dMoments[dmiy[i]] = computeLimiter(moms[i], momData.up, momData.down, momData.center);
      }
      for (size_t i = 0; i < moms.size(); i++) {
         dMoments[dmiz[i]] = computeLimiter(moms[i], momData.near, momData.far, momData.center);
      }
   }

   // electron pressure
   // Constants for electron pressure derivatives
   // Upstream pressure
   const Real Peupstream = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   const Real Peconst = Peupstream * pow(Parameters::electronDensity, -Parameters::electronPTindex);
   auto computePresE = [&dMoments, &Peconst](const auto& right, const auto& left, const auto& center) {
      // pres_e = const * np.power(rho_e, index)
      return Peconst * limiter(pow(left[mom::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex),
                               pow(center[mom::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex),
                               pow(right[mom::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex));
   };

   dMoments[dmo::dPedx] = computePresE(momData.right, momData.left, momData.center);
   dMoments[dmo::dPedy] = computePresE(momData.up, momData.down, momData.center);
   dMoments[dmo::dPedz] = computePresE(momData.near, momData.far, momData.center);
}

void computePerb(std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                 std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb, const fsgrid::FsStencil& stencil,
                 bool dontCompute2ndDerivatives, bool notSysBoundary, cuint sysBoundaryFlag) {
   using dpb = fsgrids::dperb;
   using bfi = fsgrids::bfield;
   std::array<Real, dpb::N_DPERB>& dPerB = dperb[stencil.center()];

   auto computeDiff = [](const auto& i, const auto& right, const auto& left) { return 0.5 * (right[i] - left[i]); };

   auto computeLimiter = [](const auto& i, const auto& right, const auto& left, const auto& center) {
      return limiter(left[i], center[i], right[i]);
   };
   const DerivativesData perbData{
       perb[stencil.center()], perb[stencil.right()], perb[stencil.left()], perb[stencil.up()],
       perb[stencil.down()],   perb[stencil.near()],  perb[stencil.far()],
   };

   if (notSysBoundary) {
      dPerB[dpb::dPERBydx] = computeDiff(bfi::PERBY, perbData.right, perbData.left);
      dPerB[dpb::dPERBzdx] = computeDiff(bfi::PERBZ, perbData.right, perbData.left);
      dPerB[dpb::dPERBxdy] = computeDiff(bfi::PERBX, perbData.up, perbData.down);
      dPerB[dpb::dPERBzdy] = computeDiff(bfi::PERBZ, perbData.up, perbData.down);
      dPerB[dpb::dPERBxdz] = computeDiff(bfi::PERBX, perbData.near, perbData.far);
      dPerB[dpb::dPERBydz] = computeDiff(bfi::PERBY, perbData.near, perbData.far);
   } else {
      dPerB[dpb::dPERBydx] = computeLimiter(bfi::PERBY, perbData.right, perbData.left, perbData.center);
      dPerB[dpb::dPERBzdx] = computeLimiter(bfi::PERBZ, perbData.right, perbData.left, perbData.center);
      dPerB[dpb::dPERBxdy] = computeLimiter(bfi::PERBX, perbData.up, perbData.down, perbData.center);
      dPerB[dpb::dPERBzdy] = computeLimiter(bfi::PERBZ, perbData.up, perbData.down, perbData.center);
      dPerB[dpb::dPERBxdz] = computeLimiter(bfi::PERBX, perbData.near, perbData.far, perbData.center);
      dPerB[dpb::dPERBydz] = computeLimiter(bfi::PERBY, perbData.near, perbData.far, perbData.center);
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
      dPerB[dpb::dPERBydxx] = compute2ndDerivative(bfi::PERBY, perbData.right, perbData.left, perbData.center);
      dPerB[dpb::dPERBzdxx] = compute2ndDerivative(bfi::PERBZ, perbData.right, perbData.left, perbData.center);
      dPerB[dpb::dPERBxdyy] = compute2ndDerivative(bfi::PERBX, perbData.up, perbData.down, perbData.center);
      dPerB[dpb::dPERBzdyy] = compute2ndDerivative(bfi::PERBZ, perbData.up, perbData.down, perbData.center);
      dPerB[dpb::dPERBxdzz] = compute2ndDerivative(bfi::PERBX, perbData.near, perbData.far, perbData.center);
      dPerB[dpb::dPERBydzz] = compute2ndDerivative(bfi::PERBY, perbData.near, perbData.far, perbData.center);

      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         auto compute = [&perb](auto bl, auto br, auto tl, auto tr, auto i) {
            const auto& botLeft = perb[bl];
            const auto& botRght = perb[br];
            const auto& topLeft = perb[tl];
            const auto& topRght = perb[tr];
            return FOURTH * (botLeft[i] + topRght[i] - botRght[i] - topLeft[i]);
         };

         dPerB[dpb::dPERBxdyz] =
             compute(stencil.downfar(), stencil.upfar(), stencil.downnear(), stencil.upnear(), bfi::PERBX);
         dPerB[dpb::dPERBydxz] =
             compute(stencil.leftfar(), stencil.rightfar(), stencil.leftnear(), stencil.rightnear(), bfi::PERBY);
         dPerB[dpb::dPERBzdxy] =
             compute(stencil.leftdown(), stencil.rightdown(), stencil.leftup(), stencil.rightup(), bfi::PERBZ);
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
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * \param doMoments Bool telling whether the derivatives for moments need updating too.
 *
 * \sa calculateDerivativesSimple calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivatives(std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                          std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments,
                          std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb,
                          std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                          const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer,
                          const bool doMoments) {
   /*
    * For sysBoundaryLayer 1 or 2, we are near a boundary, and we wish to use regular centered differences instead of
    * slope limiter-adjusted values. This is to minimize oscillations as a smooth behaviour is required near artificial
    * boundaries, unlike at boundaries and shocks inside the simulation domain.
    */
   const bool notSysBoundary =
       sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY);
   const bool dontCompute2ndDerivatives = Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1;

   if (doMoments) {
      computeMoments(moments, dmoments, stencil, notSysBoundary);
   }

   computePerb(perb, dperb, stencil, dontCompute2ndDerivatives, notSysBoundary, sysBoundaryFlag);

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
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param communicateMoments If true, the derivatives of moments (rho, V, P) are communicated to neighbours.

 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivativesSimple(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                std::span<std::array<Real, fsgrids::moments::N_MOMENTS>> moments,
                                std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb,
                                std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                                fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
                                const bool doMoments) {
   phiprof::Timer derivativesTimer{"Calculate face derivatives"};
   const size_t numCells = technicalGrid.getNumCells();

   phiprof::Timer mpiTimer{"FS derivatives ghost updates MPI", {"MPI"}};
   technicalGrid.updateGhostCells(perb);
   if (doMoments) {
      technicalGrid.updateGhostCells(moments);
   }
   mpiTimer.stop();

   // Calculate derivatives
   technicalGrid.parallel_for(
       [=](const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
          calculateDerivatives(perb, moments, dperb, dmoments, stencil, sysBoundaryFlag, sysBoundaryLayer, doMoments);
       },
       [](int timerId) { return phiprof::Timer{timerId}; }, phiprof::initializeTimer("FS derivatives compute cells"));

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
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateBVOLDerivatives(std::span<std::array<Real, fsgrids::volfields::N_VOL>> vol,
                              std::span<const fsgrids::technical> technical, const fsgrid::FsStencil& stencil) {
   const auto& tech = technical[stencil.center()];

   cuint sysBoundaryFlag = tech.sysBoundaryFlag;
   cuint sysBoundaryLayer = tech.sysBoundaryLayer;
   const bool notSysBoundary =
       sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY);

   using vf = fsgrids::volfields;

   auto computeDiff = [](const auto& i, const auto& right, const auto& left) { return 0.5 * (right[i] - left[i]); };

   auto computeLimiter = [](const auto& i, const auto& right, const auto& left, const auto& center) {
      return limiter(left[i], center[i], right[i]);
   };

   auto& volCenter = vol[stencil.center()];
   const auto right = stencil.right();
   const auto left = stencil.left();
   if (notSysBoundary) {
      volCenter[vf::dPERBXVOLdx] = computeDiff(vf::PERBXVOL, vol[right], vol[left]);
      volCenter[vf::dPERBYVOLdx] = computeDiff(vf::PERBYVOL, vol[right], vol[left]);
      volCenter[vf::dPERBZVOLdx] = computeDiff(vf::PERBZVOL, vol[right], vol[left]);
   } else {
      volCenter[vf::dPERBXVOLdx] = computeLimiter(vf::PERBXVOL, vol[right], vol[left], volCenter);
      volCenter[vf::dPERBYVOLdx] = computeLimiter(vf::PERBYVOL, vol[right], vol[left], volCenter);
      volCenter[vf::dPERBZVOLdx] = computeLimiter(vf::PERBZVOL, vol[right], vol[left], volCenter);
   }

   const auto up = stencil.up();
   const auto down = stencil.down();
   if (notSysBoundary) {
      volCenter[vf::dPERBXVOLdy] = computeDiff(vf::PERBXVOL, vol[up], vol[down]);
      volCenter[vf::dPERBYVOLdy] = computeDiff(vf::PERBYVOL, vol[up], vol[down]);
      volCenter[vf::dPERBZVOLdy] = computeDiff(vf::PERBZVOL, vol[up], vol[down]);
   } else {
      volCenter[vf::dPERBXVOLdy] = computeLimiter(vf::PERBXVOL, vol[up], vol[down], volCenter);
      volCenter[vf::dPERBYVOLdy] = computeLimiter(vf::PERBYVOL, vol[up], vol[down], volCenter);
      volCenter[vf::dPERBZVOLdy] = computeLimiter(vf::PERBZVOL, vol[up], vol[down], volCenter);
   }

   const auto near = stencil.near();
   const auto far = stencil.far();
   if (notSysBoundary) {
      volCenter[vf::dPERBXVOLdz] = computeDiff(vf::PERBXVOL, vol[near], vol[far]);
      volCenter[vf::dPERBYVOLdz] = computeDiff(vf::PERBYVOL, vol[near], vol[far]);
      volCenter[vf::dPERBZVOLdz] = computeDiff(vf::PERBZVOL, vol[near], vol[far]);
   } else {
      volCenter[vf::dPERBXVOLdz] = computeLimiter(vf::PERBXVOL, vol[near], vol[far], volCenter);
      volCenter[vf::dPERBYVOLdz] = computeLimiter(vf::PERBYVOL, vol[near], vol[far], volCenter);
      volCenter[vf::dPERBZVOLdz] = computeLimiter(vf::PERBZVOL, vol[near], vol[far], volCenter);
   }
}

/*! \brief High-level derivative calculation wrapper function.
 *
 * BVOL has been calculated locally by calculateVolumeAveragedFields but not communicated.
 * For the acceleration step one needs the cross-derivatives of BVOL
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateBVOLDerivativesSimple(std::span<std::array<Real, fsgrids::volfields::N_VOL>> vol,
                                    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   std::span<const fsgrids::technical> technical = technicalGrid.getData();
   const auto* localSize = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = localSize[0] * localSize[1] * localSize[2];

   phiprof::Timer derivsTimer{"Calculate volume derivatives"};
   int computeTimerId{phiprof::initializeTimer("FS derivatives BVOL compute cells")};

   phiprof::Timer commTimer{"BVOL derivatives ghost updates MPI", {"MPI"}};
   technicalGrid.updateGhostCells(vol);
   commTimer.stop(N_cells, "Spatial Cells");

// Calculate derivatives
#pragma omp parallel
   {
      phiprof::Timer computeTimer{computeTimerId};
#pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const auto stencil = technicalGrid.makeStencil(i, j, k);
               calculateBVOLDerivatives(vol, technical, stencil);
            }
         }
      }
      computeTimer.stop(N_cells, "Spatial Cells");
   }

   derivsTimer.stop(N_cells, "Spatial Cells");
}

/*! \brief Low-level curvature calculation.
 *
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param bgbGrid fsGrid holding the background fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
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

   const auto [bx, by, bz] = compute(stencil.center());
   const auto [left_x_bx, left_x_by, left_x_bz] = compute(stencil.left());
   const auto [rght_x_bx, rght_x_by, rght_x_bz] = compute(stencil.right());
   const auto [left_y_bx, left_y_by, left_y_bz] = compute(stencil.down());
   const auto [rght_y_bx, rght_y_by, rght_y_bz] = compute(stencil.up());
   const auto [left_z_bx, left_z_by, left_z_bz] = compute(stencil.far());
   const auto [rght_z_bx, rght_z_by, rght_z_bz] = compute(stencil.near());

   auto& volCenter = vol[stencil.center()];
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
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateCurvatureSimple(std::span<std::array<Real, fsgrids::volfields::N_VOL>> vol,
                              std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                              fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   const auto& gridSpacing = technicalGrid.getGridSpacing();
   const auto* localSize = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = localSize[0] * localSize[1] * localSize[2];

   std::span<const fsgrids::technical> technical = technicalGrid.getData();

   phiprof::Timer curvatureTimer{"Calculate curvature"};
   int computeTimerId{phiprof::initializeTimer("Calculate curvature compute cells")};

   phiprof::Timer commTimer{"Calculate curvature ghost updates MPI", {"MPI"}};
   technicalGrid.updateGhostCells(vol);
   commTimer.stop(N_cells, "Spatial Cells");

#pragma omp parallel
   {
      phiprof::Timer computeTimer{computeTimerId};
#pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const auto stencil = technicalGrid.makeStencil(i, j, k);
               const auto& tech = technical[stencil.center()];

               const bool compute = (tech.sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY &&
                                     tech.sysBoundaryLayer != 1 && tech.sysBoundaryLayer != 2);
               if (compute) {
                  calculateCurvature(vol, bgb, stencil, gridSpacing);
               }
            }
         }
      }
      computeTimer.stop(N_cells, "Spatial Cells");
   }

   curvatureTimer.stop(N_cells, "Spatial Cells");
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
