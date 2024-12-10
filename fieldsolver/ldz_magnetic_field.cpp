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

#include <span>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "ldz_magnetic_field.hpp"

void propagateMagneticFieldComponent(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                     std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perbdt2,
                                     std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> e,
                                     std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> edt2,
                                     const fsgrid::FsStencil& stencil, Real dt0, Real dt1, int32_t RKCase,
                                     size_t perbIdx, size_t i, size_t j, const std::array<size_t, 3>& idx) {
   const auto& e0 = RKCase == RK_ORDER2_STEP2 ? edt2[idx[0]] : e[idx[0]];
   const auto& e1 = RKCase == RK_ORDER2_STEP2 ? edt2[idx[1]] : e[idx[1]];
   const auto& e2 = RKCase == RK_ORDER2_STEP2 ? edt2[idx[2]] : e[idx[2]];
   auto& pb = perb[stencil.center()][perbIdx];
   const Real v = dt0 * (e2[i] - e0[i]) + dt1 * (e0[j] - e1[j]);

   switch (RKCase) {
   case RK_ORDER1 | RK_ORDER2_STEP2:
      pb += v;
      break;
   case RK_ORDER2_STEP1:
      perbdt2[stencil.center()][perbIdx] = pb + 0.5 * v;
      break;
   default:
      std::cerr << __FILE__ << ":" << __LINE__ << ":"
                << "Invalid RK case." << std::endl;
      abort();
   }
}

/*! \brief Low-level magnetic field propagation function.
 *
 * Propagates the cell's face-averaged magnetic field components by
 * using Faraday's law on the face edges. Depending on the time order
 * of accuracy it is done in one stage or in two stages using the
 * intermediate E1 components for the first stage of the second-order
 * Runge-Kutta method and E for the other cases.
 *
 *
 * \param perBGrid fsGrid holding the perturbed B quantities at runge-kutta t=0
 * \param perBDt2Grid fsGrid holding the perturbed B quantities at runge-kutta t=0.5
 * \param EGrid fsGrid holding the Electric field quantities at runge-kutta t=0
 * \param EDt2Grid fsGrid holding the Electric field quantities at runge-kutta t=0.5
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param dt Length of the time step
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param doX If true, compute the x component (default true).
 * \param doY If true, compute the y component (default true).
 * \param doZ If true, compute the z component (default true).
 */
void propagateMagneticField(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                            std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perbdt2,
                            std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> e,
                            std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> edt2,
                            const fsgrid::FsStencil& stencil, Real dt, int32_t RKCase, bool doX, bool doY, bool doZ,
                            const std::array<Real, 3>& gridSpacing) {
   creal dtdx = dt / gridSpacing[0];
   creal dtdy = dt / gridSpacing[1];
   creal dtdz = dt / gridSpacing[2];

   if (doX == true) {
      propagateMagneticFieldComponent(perb, perbdt2, e, edt2, stencil, dtdz, dtdy, RKCase, fsgrids::bfield::PERBX,
                                      fsgrids::efield::EY, fsgrids::efield::EZ,
                                      {stencil.center(), stencil.up(), stencil.near()});
   }

   if (doY == true) {
      propagateMagneticFieldComponent(perb, perbdt2, e, edt2, stencil, dtdx, dtdz, RKCase, fsgrids::bfield::PERBY,
                                      fsgrids::efield::EZ, fsgrids::efield::EX,
                                      {stencil.center(), stencil.near(), stencil.right()});
   }

   if (doZ == true) {
      propagateMagneticFieldComponent(perb, perbdt2, e, edt2, stencil, dtdy, dtdx, RKCase, fsgrids::bfield::PERBZ,
                                      fsgrids::efield::EX, fsgrids::efield::EY,
                                      {stencil.center(), stencil.right(), stencil.up()});
   }
}

/*! \brief Low-level magnetic field propagation function.
 *
 * Propagates the magnetic field according to the system boundary conditions.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities at runge-kutta t=0
 * \param perBDt2Grid fsGrid holding the perturbed B quantities at runge-kutta t=0.5
 * \param bgbGrid fsGrid holding the background field B quantities
 * \param EGrid fsGrid holding the Electric field quantities at runge-kutta t=0
 * \param EDt2Grid fsGrid holding the Electric field quantities at runge-kutta t=0.5
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 * \param dt Length of the time step
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 *
 * \sa propagateMagneticFieldSimple propagateMagneticField
 */
void propagateSysBoundaryMagneticField(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& bgbGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, int32_t i, int32_t j, int32_t k,
    SysBoundary& sysBoundaries, Real dt, int32_t RKCase, uint32_t component) {
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      perBGrid.get(i, j, k)->at(fsgrids::bfield::PERBX + component) =
          sysBoundaries.getSysBoundary(technicalGrid.get(i, j, k)->sysBoundaryFlag)
              ->fieldSolverBoundaryCondMagneticField(perBGrid, bgbGrid, technicalGrid, i, j, k, dt, component);
   } else {
      perBDt2Grid.get(i, j, k)->at(fsgrids::bfield::PERBX + component) =
          sysBoundaries.getSysBoundary(technicalGrid.get(i, j, k)->sysBoundaryFlag)
              ->fieldSolverBoundaryCondMagneticField(perBDt2Grid, bgbGrid, technicalGrid, i, j, k, dt, component);
   }
}

/*! \brief High-level magnetic field propagation function.
 *
 * Propagates the magnetic field and applies the field boundary conditions defined in project.h where needed.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities at runge-kutta t=0
 * \param perBDt2Grid fsGrid holding the perturbed B quantities at runge-kutta t=0.5
 * \param EGrid fsGrid holding the Electric field quantities at runge-kutta t=0
 * \param EDt2Grid fsGrid holding the Electric field quantities at runge-kutta t=0.5
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * \param dt Length of the time step
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 *
 * \sa propagateMagneticField propagateSysBoundaryMagneticField
 */
void propagateMagneticFieldSimple(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& bgbGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EDt2Grid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, SysBoundary& sysBoundaries, creal& dt,
    cint& RKCase) {
   const auto& localSize = technicalGrid.getLocalSize();
   const auto& gridSpacing = technicalGrid.getGridSpacing();
   const size_t N_cells = localSize[0] * localSize[1] * localSize[2];

   std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb = perBGrid.getData();
   std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perbdt2 = perBDt2Grid.getData();
   std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb = bgbGrid.getData();
   std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> e = EGrid.getData();
   std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> edt2 = EDt2Grid.getData();
   std::span<const fsgrids::technical> technical = technicalGrid.getData();

   phiprof::Timer propagateBTimer{"Propagate magnetic field"};
   int computeTimerId{phiprof::initializeTimer("Magnetic Field compute cells")};
   int sysBoundaryTimerId{phiprof::initializeTimer("Magnetic Field compute sysboundary cells")};
#pragma omp parallel
   {
      phiprof::Timer computeTimer{computeTimerId};
#pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const fsgrid::FsStencil stencil = technicalGrid.makeStencil(i, j, k);
               cuint bitfield = technical[stencil.center()].SOLVE;
               propagateMagneticField(
                   perb, perbdt2, e, edt2, stencil, dt, RKCase, ((bitfield & compute::BX) == compute::BX),
                   ((bitfield & compute::BY) == compute::BY), ((bitfield & compute::BZ) == compute::BZ), gridSpacing);
            }
         }
      }
   }

   // This communication is needed for boundary conditions, in practice almost all
   // of the communication is going to be redone in calculateDerivativesSimple
   // TODO: do not transfer if there are no field boundaryconditions
   phiprof::Timer mpiTimer{"MPI", {"MPI"}};
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ with neighbours
      perBGrid.updateGhostCells();
   } else { // RKCase == RK_ORDER2_STEP1
      // Exchange PERBX_DT2,PERBY_DT2,PERBZ_DT2 with neighbours
      perBDt2Grid.updateGhostCells();
   }
   mpiTimer.stop();

// Propagate B on system boundary/process inner cells
#pragma omp parallel
   {
      phiprof::Timer sysBoundaryTimer{sysBoundaryTimerId};
// L1 pass
#pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const fsgrid::FsStencil stencil = technicalGrid.makeStencil(i, j, k);
               const auto tech = technical[stencil.center()];
               cuint bitfield = tech.SOLVE;
               // L1 pass
               if (tech.sysBoundaryLayer == 1) {
                  if ((bitfield & compute::BX) != compute::BX) {
                     propagateSysBoundaryMagneticField(perBGrid, perBDt2Grid, bgbGrid, technicalGrid, i, j, k,
                                                       sysBoundaries, dt, RKCase, 0);
                  }
                  if ((bitfield & compute::BY) != compute::BY) {
                     propagateSysBoundaryMagneticField(perBGrid, perBDt2Grid, bgbGrid, technicalGrid, i, j, k,
                                                       sysBoundaries, dt, RKCase, 1);
                  }
                  if ((bitfield & compute::BZ) != compute::BZ) {
                     propagateSysBoundaryMagneticField(perBGrid, perBDt2Grid, bgbGrid, technicalGrid, i, j, k,
                                                       sysBoundaries, dt, RKCase, 2);
                  }
               }
            }
         }
      }
   }

   mpiTimer.start();
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ with neighbours
      perBGrid.updateGhostCells();
   } else { // RKCase == RK_ORDER2_STEP1
      // Exchange PERBX_DT2,PERBY_DT2,PERBZ_DT2 with neighbours
      perBDt2Grid.updateGhostCells();
   }
   mpiTimer.stop();

#pragma omp parallel
   {
      phiprof::Timer sysBoundaryTimer{sysBoundaryTimerId};
// L2 pass
#pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const fsgrid::FsStencil stencil = technicalGrid.makeStencil(i, j, k);
               const auto tech = technical[stencil.center()];
               if (tech.sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY && tech.sysBoundaryLayer == 2) {
                  for (uint component = 0; component < 3; component++) {
                     propagateSysBoundaryMagneticField(perBGrid, perBDt2Grid, bgbGrid, technicalGrid, i, j, k,
                                                       sysBoundaries, dt, RKCase, component);
                  }
               }
            }
         }
      }
   }
   propagateBTimer.stop(N_cells, "Spatial Cells");
}
