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

// clang-format off
#include "fs_common.h"
#include <limits>
#include "ldz_gradpe.hpp"
// clang-format on

#ifdef DEBUG_VLASIATOR
#define DEBUG_FSOLVER
#endif

using namespace std;

/** Calculate the electron pressure gradient term on all given cells.
 * @param sysBoundaries System boundary condition functions.
 */
void calculateGradPeTerm(std::span<std::array<Real, fsgrids::egradpe::N_EGRADPE>> egradpes,
                         std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments,
                         std::span<const std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                         std::span<const fsgrids::technical> technical, const fsgrid::FsStencil& stencil,
                         const auto& gridSpacing, SysBoundary& sysBoundaries) {
#ifdef DEBUG_FSOLVER
   if (stencil.center() >= moments.size()) {
      cerr << "Out-of-bounds access in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
#endif
   auto& egradpe = egradpes[stencil.center()];
   const auto& moment = moments[stencil.center()];
   const auto& dmoment = dmoments[stencil.center()];
   const auto& tech = technical[stencil.center()];

   cuint cellSysBoundaryFlag = tech.sysBoundaryFlag;
   cuint cellSysBoundaryLayer = tech.sysBoundaryLayer;

   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
       cellSysBoundaryFlag == sysboundarytype::OUTER_BOUNDARY_PADDING) {
      return;
   }

   if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && (cellSysBoundaryLayer != 1)) {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)
          ->fieldSolverBoundaryCondGradPeElectricField(egradpes, stencil, 0);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)
          ->fieldSolverBoundaryCondGradPeElectricField(egradpes, stencil, 1);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)
          ->fieldSolverBoundaryCondGradPeElectricField(egradpes, stencil, 2);
   } else {
      if (Parameters::ohmGradPeTerm == 0) {
         cerr << __FILE__ << __LINE__
              << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0."
              << endl;
      } else if (Parameters::ohmGradPeTerm == 1) {
         auto calculateEdgeGradPeTermComponent = [&egradpe, &moment, &dmoment](size_t i, size_t j, Real spacing) {
            const Real rhoq = moment[fsgrids::moments::RHOQ];
            const Real min = Parameters::hallMinimumRhoq;
            const Real max = std::numeric_limits<Real>::max();
            const Real hallRhoq = std::clamp(rhoq, min, max);
            egradpe[i] = -dmoment[j] / (hallRhoq * spacing);
         };
         calculateEdgeGradPeTermComponent(fsgrids::egradpe::EXGRADPE, fsgrids::dmoments::dPedx, gridSpacing[0]);
         calculateEdgeGradPeTermComponent(fsgrids::egradpe::EYGRADPE, fsgrids::dmoments::dPedy, gridSpacing[1]);
         calculateEdgeGradPeTermComponent(fsgrids::egradpe::EZGRADPE, fsgrids::dmoments::dPedz, gridSpacing[2]);
      } else {
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms."
              << endl;
      }
   }
}

void calculateGradPeTermSimple(
    fsgrid::FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH>& EGradPeGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH>& EGradPeDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsDt2Grid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, SysBoundary& sysBoundaries, cint& RKCase) {
   const auto& gridSpacing = technicalGrid.getGridSpacing();
   const auto& localSize = technicalGrid.getLocalSize();
   const size_t N_cells = localSize[0] * localSize[1] * localSize[2];
   const bool case0 = RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2;

   std::span<std::array<Real, fsgrids::egradpe::N_EGRADPE>> egradpes = EGradPeGrid.getData();
   std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments = momentsGrid.getData();
   std::span<const std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments = dMomentsGrid.getData();
   std::span<const fsgrids::technical> technical = technicalGrid.getData();
   if (case0) {
      egradpes = EGradPeDt2Grid.getData();
      moments = momentsDt2Grid.getData();
      dmoments = dMomentsDt2Grid.getData();
   }

   phiprof::Timer gradPeTimer{"Calculate GradPe term"};
   int computeTimerId{phiprof::initializeTimer("EgradPe compute cells")};

   phiprof::Timer mpiTimer{"EgradPe field update ghosts MPI", {"MPI"}};
   if (case0) {
      dMomentsGrid.updateGhostCells();
   } else {
      dMomentsDt2Grid.updateGhostCells();
   }
   mpiTimer.stop();

// Calculate GradPe term
#pragma omp parallel
   {
      phiprof::Timer computeTimer{computeTimerId};
#pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const auto& stencil = technicalGrid.makeStencil(i, j, k);
               calculateGradPeTerm(egradpes, moments, dmoments, technical, stencil, gridSpacing, sysBoundaries);
            }
         }
      }
      computeTimer.stop(N_cells, "Spatial Cells");
   }

   gradPeTimer.stop(N_cells, "Spatial Cells");
}
