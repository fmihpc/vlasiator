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

void propagateMagneticField(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                            std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perbdt2,
                            std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> e,
                            std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> edt2,
                            const fsgrid::FsStencil& stencil, Real dt, int32_t RKCase, bool doX, bool doY, bool doZ,
                            const std::array<Real, 3>& gridSpacing) {
   creal dtdx = dt / gridSpacing[0];
   creal dtdy = dt / gridSpacing[1];
   creal dtdz = dt / gridSpacing[2];

   std::array<Real, fsgrids::bfield::N_BFIELD>& perBGrid0 = perb[stencil.center()];

   if (doX == true) {
      switch (RKCase) {
      case RK_ORDER1: {
         const auto& EGrid0 = e[stencil.center()]; // i,j,k;
         const auto& EGrid1 = e[stencil.up()];     // i,j+1,k;
         const auto& EGrid2 = e[stencil.near()];   // i,j,k+1;
         perBGrid0[fsgrids::bfield::PERBX] += dtdz * (EGrid2[fsgrids::efield::EY] - EGrid0[fsgrids::efield::EY]) +
                                              dtdy * (EGrid0[fsgrids::efield::EZ] - EGrid1[fsgrids::efield::EZ]);
         break;
      }

      case RK_ORDER2_STEP1: {
         auto& perBDt2Grid0 = perbdt2[stencil.center()];       // i,j,k;
         const auto& EGrid0 = e[stencil.center()];             // i,j,k;
         const auto& EGrid1 = e[stencil.up()];                 // i,j+1,k;
         const auto& EGrid2 = e[stencil.near()];               // i,j,k+1;
         perBDt2Grid0[fsgrids::bfield::PERBX] =
             perBGrid0[fsgrids::bfield::PERBX] +
             0.5 * (dtdz * (EGrid2[fsgrids::efield::EY] - EGrid0[fsgrids::efield::EY]) +
                    dtdy * (EGrid0[fsgrids::efield::EZ] - EGrid1[fsgrids::efield::EZ]));
         break;
      }

      case RK_ORDER2_STEP2: {
         const auto& EGrid0 = edt2[stencil.center()]; // i,j,k;
         const auto& EGrid1 = edt2[stencil.up()];     // i,j+1,k;
         const auto& EGrid2 = edt2[stencil.near()];   // i,j,k+1;
         perBGrid0[fsgrids::bfield::PERBX] += (dtdz * (EGrid2[fsgrids::efield::EY] - EGrid0[fsgrids::efield::EY]) +
                                               dtdy * (EGrid0[fsgrids::efield::EZ] - EGrid1[fsgrids::efield::EZ]));
         break;
      }

      default:
         std::cerr << __FILE__ << ":" << __LINE__ << ":"
                   << "Invalid RK case." << std::endl;
         abort();
      }
   }

   if (doY == true) {
      switch (RKCase) {
      case RK_ORDER1: {
         const auto& EGrid0 = e[stencil.center()]; // i,j,k;
         const auto& EGrid1 = e[stencil.near()];   // i,j,k+1;
         const auto& EGrid2 = e[stencil.right()];  // i+1,j,k;
         perBGrid0[fsgrids::bfield::PERBY] += dtdx * (EGrid2[fsgrids::efield::EZ] - EGrid0[fsgrids::efield::EZ]) +
                                              dtdz * (EGrid0[fsgrids::efield::EX] - EGrid1[fsgrids::efield::EX]);
         break;
      }
      case RK_ORDER2_STEP1: {
         auto& perBDt2Grid0 = perbdt2[stencil.center()];       // i,j,k;
         const auto& EGrid0 = e[stencil.center()];             // i,j,k;
         const auto& EGrid1 = e[stencil.near()];               // i,j,k+1;
         const auto& EGrid2 = e[stencil.right()];              // i+1,j,k;
         perBDt2Grid0[fsgrids::bfield::PERBY] =
             perBGrid0[fsgrids::bfield::PERBY] +
             0.5 * (dtdx * (EGrid2[fsgrids::efield::EZ] - EGrid0[fsgrids::efield::EZ]) +
                    dtdz * (EGrid0[fsgrids::efield::EX] - EGrid1[fsgrids::efield::EX]));
         break;
      }
      case RK_ORDER2_STEP2: {
         const auto& EGrid0 = edt2[stencil.center()]; // i,j,k;
         const auto& EGrid1 = edt2[stencil.near()];   // i,j,k+1;
         const auto& EGrid2 = edt2[stencil.right()];  // i+1,j,k;
         perBGrid0[fsgrids::bfield::PERBY] += (dtdx * (EGrid2[fsgrids::efield::EZ] - EGrid0[fsgrids::efield::EZ]) +
                                               dtdz * (EGrid0[fsgrids::efield::EX] - EGrid1[fsgrids::efield::EX]));
         break;
      }
      default:
         std::cerr << __FILE__ << ":" << __LINE__ << ":"
                   << "Invalid RK case." << std::endl;
         abort();
      }
   }

   if (doZ == true) {
      switch (RKCase) {
      case RK_ORDER1: {
         const auto& EGrid0 = e[stencil.center()]; // i,j,k;
         const auto& EGrid1 = e[stencil.right()];  // i+1,j,k;
         const auto& EGrid2 = e[stencil.up()];     // i,j+1,k;
         perBGrid0[fsgrids::bfield::PERBZ] += dtdy * (EGrid2[fsgrids::efield::EX] - EGrid0[fsgrids::efield::EX]) +
                                              dtdx * (EGrid0[fsgrids::efield::EY] - EGrid1[fsgrids::efield::EY]);
         break;
      }
      case RK_ORDER2_STEP1: {
         auto& perBDt2Grid0 = perbdt2[stencil.center()];       // i,j,k;
         const auto& EGrid0 = e[stencil.center()];             // i,j,k;
         const auto& EGrid1 = e[stencil.right()];              // i+1,j,k;
         const auto& EGrid2 = e[stencil.up()];                 // i,j+1,k;
         perBDt2Grid0[fsgrids::bfield::PERBZ] =
             perBGrid0[fsgrids::bfield::PERBZ] +
             0.5 * (dtdy * (EGrid2[fsgrids::efield::EX] - EGrid0[fsgrids::efield::EX]) +
                    dtdx * (EGrid0[fsgrids::efield::EY] - EGrid1[fsgrids::efield::EY]));
         break;
      }
      case RK_ORDER2_STEP2: {
         const auto& EGrid0 = edt2[stencil.center()]; // i,j,k;
         const auto& EGrid1 = edt2[stencil.right()];  // i+1,j,k;
         const auto& EGrid2 = edt2[stencil.up()];     // i,j+1,k;
         perBGrid0[fsgrids::bfield::PERBZ] += (dtdy * (EGrid2[fsgrids::efield::EX] - EGrid0[fsgrids::efield::EX]) +
                                               dtdx * (EGrid0[fsgrids::efield::EY] - EGrid1[fsgrids::efield::EY]));
         break;
      }
      default:
         std::cerr << __FILE__ << ":" << __LINE__ << ":"
                   << "Invalid RK case." << std::endl;
         abort();
      }
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
void propagateSysBoundaryMagneticField(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                       std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perbdt2,
                                       std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                                       std::span<const fsgrids::technical> technical,
                                       const std::array<Real, 3>& gridSpacing,
                                       const std::array<fsgrid::FsSize_t, 3>& globalCoordinates,
                                       const fsgrid::FsStencil& stencil, SysBoundary& sysBoundaries, int32_t RKCase,
                                       uint32_t component) {
   const bool case0 = RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2;
   auto& out = case0 ? perb[stencil.center()] : perbdt2[stencil.center()];
   const auto& pb = case0 ? perb : perbdt2;

   out[fsgrids::bfield::PERBX + component] =
       sysBoundaries.getSysBoundary(technical[stencil.center()].sysBoundaryFlag)
           ->fieldSolverBoundaryCondMagneticField(pb, bgb, technical, gridSpacing, globalCoordinates, stencil,
                                                  component);
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
void propagateMagneticFieldSimple(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                  std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perbdt2,
                                  std::span<std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                                  std::span<std::array<Real, fsgrids::efield::N_EFIELD>> e,
                                  std::span<std::array<Real, fsgrids::efield::N_EFIELD>> edt2,
                                  fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
                                  SysBoundary& sysBoundaries, creal& dt, cint& RKCase) {
   const auto* localSize = &technicalGrid.getLocalSize()[0];
   const auto& gridSpacing = technicalGrid.getGridSpacing();
   const size_t numCells = technicalGrid.getNumCells();

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
      computeTimer.stop(numCells, "Spatial Cells");
   }

   // This communication is needed for boundary conditions, in practice almost all
   // of the communication is going to be redone in calculateDerivativesSimple
   // TODO: do not transfer if there are no field boundaryconditions
   phiprof::Timer mpiTimer{"MPI", {"MPI"}};
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ with neighbours
      technicalGrid.updateGhostCells(perb);
   } else { // RKCase == RK_ORDER2_STEP1
      // Exchange PERBX_DT2,PERBY_DT2,PERBZ_DT2 with neighbours
      technicalGrid.updateGhostCells(perbdt2);
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
               const auto globalCoordinates = technicalGrid.localToGlobal(stencil.i, stencil.j, stencil.k);
               const auto tech = technical[stencil.center()];
               cuint bitfield = tech.SOLVE;
               // L1 pass
               if (tech.sysBoundaryLayer == 1) {
                  if ((bitfield & compute::BX) != compute::BX) {
                     propagateSysBoundaryMagneticField(perb, perbdt2, bgb, technical, gridSpacing, globalCoordinates,
                                                       stencil, sysBoundaries, RKCase, 0);
                  }
                  if ((bitfield & compute::BY) != compute::BY) {
                     propagateSysBoundaryMagneticField(perb, perbdt2, bgb, technical, gridSpacing, globalCoordinates,
                                                       stencil, sysBoundaries, RKCase, 1);
                  }
                  if ((bitfield & compute::BZ) != compute::BZ) {
                     propagateSysBoundaryMagneticField(perb, perbdt2, bgb, technical, gridSpacing, globalCoordinates,
                                                       stencil, sysBoundaries, RKCase, 2);
                  }
               }
            }
         }
      }
      sysBoundaryTimer.stop(numCells, "Spatial Cells");
   }

   mpiTimer.start();
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ with neighbours
      technicalGrid.updateGhostCells(perb);
   } else { // RKCase == RK_ORDER2_STEP1
      // Exchange PERBX_DT2,PERBY_DT2,PERBZ_DT2 with neighbours
      technicalGrid.updateGhostCells(perbdt2);
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
               const auto globalCoordinates = technicalGrid.localToGlobal(stencil.i, stencil.j, stencil.k);
               const auto tech = technical[stencil.center()];
               if (tech.sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY && tech.sysBoundaryLayer == 2) {
                  for (uint component = 0; component < 3; component++) {
                     propagateSysBoundaryMagneticField(perb, perbdt2, bgb, technical, gridSpacing, globalCoordinates,
                                                       stencil, sysBoundaries, RKCase, component);
                  }
               }
            }
         }
      }
      sysBoundaryTimer.stop(numCells, "Spatial Cells");
   }
   propagateBTimer.stop(numCells, "Spatial Cells");
}
