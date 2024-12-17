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

auto computeDerivative(bool notSysBoundary, const auto& i, const auto& right, const auto& left, const auto& center) {
   return notSysBoundary ? 0.5 * (right[i] - left[i]) : limiter(left[i], center[i], right[i]);
};

/*! \brief Low-level spatial derivatives calculation.
 *
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in
 * project.h. Uses RHO, V[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the
 * second-order method, and RHO_DT2, V[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 *
 * For sysBoundaryLayer 1 or 2, we are near a boundary, and we wish to use regular centered differences instead of slope
 * limiter-adjusted values. This is to minimize oscillations as a smooth behaviour is required near artificial
 * boundaries, unlike at boundaries and shocks inside the simulation domain.
 *
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * \param shouldCalculateMoments Bool telling whether the derivatives for moments need updating too.
 *
 * \sa calculateDerivativesSimple calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivatives(std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                          std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments,
                          std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb,
                          std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                          std::span<const fsgrids::technical> technical, const fsgrid::FsStencil& stencil,
                          const bool shouldCalculateMoments) {
   std::array<Real, fsgrids::dperb::N_DPERB>& dPerB = dperb[stencil.center()];
   std::array<Real, fsgrids::dmoments::N_DMOMENTS>& dMoments = dmoments[stencil.center()];
   const auto& tech = technical[stencil.center()];

   // Get boundary flag for the cell:
   cuint sysBoundaryFlag = tech.sysBoundaryFlag;
   cuint sysBoundaryLayer = tech.sysBoundaryLayer;
   const bool notSysBoundary =
       sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY);
   const bool dontCompute2ndDerivatives = Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1;

   // Constants for electron pressure derivatives
   // Upstream pressure
   const Real Peupstream = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   const Real Peconst = Peupstream * pow(Parameters::electronDensity, -Parameters::electronPTindex);

   // clang-format off
   static constexpr std::array<std::array<fsgrids::dmoments, 8>, 3> dmomentsIndices = {
       std::array {
           fsgrids::dmoments::drhomdx,
           fsgrids::dmoments::drhoqdx,
           fsgrids::dmoments::dp11dx,
           fsgrids::dmoments::dp22dx,
           fsgrids::dmoments::dp33dx,
           fsgrids::dmoments::dVxdx,
           fsgrids::dmoments::dVydx,
           fsgrids::dmoments::dVzdx
       },
       std::array {
           fsgrids::dmoments::drhomdy,
           fsgrids::dmoments::drhoqdy,
           fsgrids::dmoments::dp11dy,
           fsgrids::dmoments::dp22dy,
           fsgrids::dmoments::dp33dy,
           fsgrids::dmoments::dVxdy,
           fsgrids::dmoments::dVydy,
           fsgrids::dmoments::dVzdy,
       },
       std::array {
         fsgrids::dmoments::drhomdz,
         fsgrids::dmoments::drhoqdz,
         fsgrids::dmoments::dp11dz,
         fsgrids::dmoments::dp22dz,
         fsgrids::dmoments::dp33dz,
         fsgrids::dmoments::dVxdz,
         fsgrids::dmoments::dVydz,
         fsgrids::dmoments::dVzdz,
       },
   };

   static constexpr std::array momentsIndices = {
       fsgrids::moments::RHOM,
       fsgrids::moments::RHOQ,
       fsgrids::moments::P_11,
       fsgrids::moments::P_22,
       fsgrids::moments::P_33,
       fsgrids::moments::VX,
       fsgrids::moments::VY,
       fsgrids::moments::VZ,
   };

   static constexpr std::array perBIndices = {
       std::array {
           fsgrids::bfield::PERBY,
           fsgrids::bfield::PERBZ,
       },
       std::array {
           fsgrids::bfield::PERBX,
           fsgrids::bfield::PERBZ,
       },
       std::array {
           fsgrids::bfield::PERBX,
           fsgrids::bfield::PERBY,
       },
   };

   static constexpr std::array dPerBIndices = {
       std::array {
          fsgrids::dperb::dPERBydx,
          fsgrids::dperb::dPERBzdx,
          fsgrids::dperb::dPERBydxx,
          fsgrids::dperb::dPERBzdxx,
       },
       std::array {
          fsgrids::dperb::dPERBxdy,
          fsgrids::dperb::dPERBzdy,
          fsgrids::dperb::dPERBxdyy,
          fsgrids::dperb::dPERBzdyy,
       },
       std::array {
          fsgrids::dperb::dPERBxdz,
          fsgrids::dperb::dPERBydz,
          fsgrids::dperb::dPERBxdzz,
          fsgrids::dperb::dPERBydzz,
       },
   };

   static constexpr std::array presEIndices = {
    fsgrids::dmoments::dPedx,
    fsgrids::dmoments::dPedy,
    fsgrids::dmoments::dPedz,
   };

   const std::array cellIndices = {
       std::array {
           stencil.left(),
           stencil.right(),
       },
       std::array {
           stencil.down(),
           stencil.up(),
       },
       std::array {
           stencil.far(),
           stencil.near(),
       },
   };
   // clang-format on

   auto computeMoments = [&shouldCalculateMoments, &dMoments, &notSysBoundary](auto component, const auto& right,
                                                                               const auto& left, const auto& center) {
      if (shouldCalculateMoments) {
         for (size_t i = 0; i < momentsIndices.size(); i++) {
            dMoments[dmomentsIndices[component][i]] =
                computeDerivative(notSysBoundary, momentsIndices[i], right, left, center);
         }
      }
   };

   auto computePresE = [&shouldCalculateMoments, &dMoments, &Peconst](auto component, const auto& right,
                                                                      const auto& left, const auto& center) {
      if (shouldCalculateMoments) {
         // pres_e = const * np.power(rho_e, index)
         dMoments[presEIndices[component]] =
             Peconst *
             limiter(pow(left[fsgrids::moments::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex),
                     pow(center[fsgrids::moments::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex),
                     pow(right[fsgrids::moments::RHOQ] / physicalconstants::CHARGE, Parameters::electronPTindex));
      }
   };

   auto computePerB = [&dPerB, &notSysBoundary, &dontCompute2ndDerivatives](auto component, const auto& right,
                                                                            const auto& left, const auto& center) {
      for (auto i = 0; i < 2; i++) {
         const auto j = dPerBIndices[component][i];
         const auto k = perBIndices[component][i];
         dPerB[j] = computeDerivative(notSysBoundary, k, right, left, center);
      }

      for (auto i = 0; i < 2; i++) {
         const auto j = dPerBIndices[component][i + 2];
         const auto k = perBIndices[component][i];

         dPerB[j] = dontCompute2ndDerivatives ? 0.0 : left[k] + right[k] - 2.0 * center[k];
      }
   };

   // Compute moments
   const std::array<Real, fsgrids::moments::N_MOMENTS>& centerMoments = moments[stencil.center()];
#ifdef DEBUG_SOLVERS
   {
      const auto& cv = centerMoments[fsgrids::moments::RHOM];
      if (cv <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << (cv < 0 ? " Negative" : " Zero") << " density in spatial cell at ("
                   << i << " " << j << " " << k << ")" << std::endl;
         abort();
      }
   }
#endif

   for (auto component = 0; component < 3; component++) {
      const auto& ci = cellIndices[component];
      const auto& left = moments[ci[0]];
      const auto& right = moments[ci[1]];
#ifdef DEBUG_SOLVERS
      {
         const auto& lv = left[fsgrids::moments::RHOM];
         if (lv <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__ << (lv < 0 ? " Negative" : " Zero") << " density in spatial cell"
                      << std::endl;
            abort();
         }

         const auto& rv = right[fsgrids::moments::RHOM];
         if (rv <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__ << (rv < 0 ? " Negative" : " Zero") << " density in spatial cell"
                      << std::endl;
            abort();
         }
      }
#endif

      computeMoments(component, right, left, centerMoments);
      computePresE(component, right, left, centerMoments);
   }

   // Compute perb
   const std::array<Real, fsgrids::bfield::N_BFIELD>& centerPerB = perb[stencil.center()];
   for (auto component = 0; component < 3; component++) {
      const auto& ci = cellIndices[component];
      const auto& left = perb[ci[0]];
      const auto& right = perb[ci[1]];
      computePerB(component, right, left, centerPerB);
   }

   if (dontCompute2ndDerivatives) {
      dPerB[fsgrids::dperb::dPERBxdyz] = 0.0;
      dPerB[fsgrids::dperb::dPERBydxz] = 0.0;
      dPerB[fsgrids::dperb::dPERBzdxy] = 0.0;
   } else if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      // clang-format off
      const std::array cellIndices = {
          std::array {
              stencil.leftdown(),
              stencil.rightdown(),
              stencil.leftup(),
              stencil.rightup(),
          },
          std::array {
              stencil.leftfar(),
              stencil.rightfar(),
              stencil.leftnear(),
              stencil.rightnear(),
          },
          std::array {
              stencil.downfar(),
              stencil.upfar(),
              stencil.downnear(),
              stencil.upnear(),
          },
      };
      // clang-format on

      static constexpr std::array dPerBIndices = {
          fsgrids::dperb::dPERBzdxy,
          fsgrids::dperb::dPERBydxz,
          fsgrids::dperb::dPERBxdyz,
      };

      static constexpr std::array perBIndices = {
          fsgrids::bfield::PERBZ,
          fsgrids::bfield::PERBY,
          fsgrids::bfield::PERBX,
      };

      for (size_t component = 0; component < dPerBIndices.size(); component++) {
         const auto& ci = cellIndices[component];
         const auto& botLeft = perb[ci[0]];
         const auto& botRght = perb[ci[1]];
         const auto& topLeft = perb[ci[2]];
         const auto& topRght = perb[ci[3]];

         const auto& i = dPerBIndices[component];
         const auto& j = perBIndices[component];
         dPerB[i] = FOURTH * (botLeft[j] + topRght[j] - botRght[j] - topLeft[j]);
      }
   } else {
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
void calculateDerivativesSimple(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, const bool doMoments) {
   const auto& localSize = technicalGrid.getLocalSize();
   const size_t N_cells = localSize[0] * localSize[1] * localSize[2];

   std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perb = perBGrid.getData();
   std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments = momentsGrid.getData();
   std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb = dPerBGrid.getData();
   std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments = dMomentsGrid.getData();
   std::span<const fsgrids::technical> technical = technicalGrid.getData();

   phiprof::Timer derivativesTimer{"Calculate face derivatives"};
   int computeTimerId{phiprof::initializeTimer("FS derivatives compute cells")};

   phiprof::Timer mpiTimer{"FS derivatives ghost updates MPI", {"MPI"}};

   perBGrid.updateGhostCells();
   if (doMoments) {
      momentsGrid.updateGhostCells();
   }

   mpiTimer.stop();

// Calculate derivatives
#pragma omp parallel
   {
      phiprof::Timer computeTimer{computeTimerId};
#pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const auto stencil = technicalGrid.makeStencil(i, j, k);
               calculateDerivatives(perb, moments, dperb, dmoments, technical, stencil, doMoments);
            }
         }
      }
      computeTimer.stop(N_cells, "Spatial Cells");
   }

   derivativesTimer.stop(N_cells, "Spatial Cells");
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
   const auto right = stencil.right();
   const auto left = stencil.left();
   const auto up = stencil.up();
   const auto down = stencil.down();
   const auto near = stencil.near();
   const auto far = stencil.far();

   auto& volCenter = vol[stencil.center()];
   volCenter[vf::dPERBXVOLdx] = computeDerivative(notSysBoundary, vf::PERBXVOL, vol[right], vol[left], volCenter);
   volCenter[vf::dPERBYVOLdx] = computeDerivative(notSysBoundary, vf::PERBYVOL, vol[right], vol[left], volCenter);
   volCenter[vf::dPERBZVOLdx] = computeDerivative(notSysBoundary, vf::PERBZVOL, vol[right], vol[left], volCenter);
   volCenter[vf::dPERBXVOLdy] = computeDerivative(notSysBoundary, vf::PERBXVOL, vol[up], vol[down], volCenter);
   volCenter[vf::dPERBYVOLdy] = computeDerivative(notSysBoundary, vf::PERBYVOL, vol[up], vol[down], volCenter);
   volCenter[vf::dPERBZVOLdy] = computeDerivative(notSysBoundary, vf::PERBZVOL, vol[up], vol[down], volCenter);
   volCenter[vf::dPERBXVOLdz] = computeDerivative(notSysBoundary, vf::PERBXVOL, vol[near], vol[far], volCenter);
   volCenter[vf::dPERBYVOLdz] = computeDerivative(notSysBoundary, vf::PERBYVOL, vol[near], vol[far], volCenter);
   volCenter[vf::dPERBZVOLdz] = computeDerivative(notSysBoundary, vf::PERBZVOL, vol[near], vol[far], volCenter);
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
void calculateBVOLDerivativesSimple(
    fsgrid::FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH>& volGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   std::span<std::array<Real, fsgrids::volfields::N_VOL>> vol = volGrid.getData();
   std::span<const fsgrids::technical> technical = technicalGrid.getData();
   const auto& localSize = technicalGrid.getLocalSize();
   const size_t N_cells = localSize[0] * localSize[1] * localSize[2];

   phiprof::Timer derivsTimer{"Calculate volume derivatives"};
   int computeTimerId{phiprof::initializeTimer("FS derivatives BVOL compute cells")};

   phiprof::Timer commTimer{"BVOL derivatives ghost updates MPI", {"MPI"}};
   volGrid.updateGhostCells();
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
void calculateCurvatureSimple(fsgrid::FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH>& volGrid,
                              fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& bgbGrid,
                              fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   const auto& gridSpacing = technicalGrid.getGridSpacing();
   const auto& localSize = technicalGrid.getLocalSize();
   const size_t N_cells = localSize[0] * localSize[1] * localSize[2];

   std::span<std::array<Real, fsgrids::volfields::N_VOL>> vol = volGrid.getData();
   std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb = bgbGrid.getData();
   std::span<const fsgrids::technical> technical = technicalGrid.getData();

   phiprof::Timer curvatureTimer{"Calculate curvature"};
   int computeTimerId{phiprof::initializeTimer("Calculate curvature compute cells")};

   phiprof::Timer commTimer{"Calculate curvature ghost updates MPI", {"MPI"}};
   volGrid.updateGhostCells();
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
