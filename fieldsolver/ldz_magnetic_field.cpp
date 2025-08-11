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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ldz_magnetic_field.hpp"

void propagateMagneticField(fsgrids::perbspan perb,
                            fsgrids::perbspan perbdt2,
                            fsgrids::constefieldspan e,
                            fsgrids::constefieldspan edt2,
                            const fsgrid::FsStencil& stencil, Real dt, int32_t RKCase, bool doX, bool doY, bool doZ,
                            const std::array<Real, 3>& gridSpacing) {
   creal dtdx = dt / gridSpacing[0];
   creal dtdy = dt / gridSpacing[1];
   creal dtdz = dt / gridSpacing[2];

   std::array<Real, fsgrids::bfield::N_BFIELD>& perBGrid0 = perb[stencil.ooo()];

   if (doX == true) {
      switch (RKCase) {
      case RK_ORDER1: {
         const auto& EGrid0 = e[stencil.ooo()]; // i,j,k;
         const auto& EGrid1 = e[stencil.opo()];     // i,j+1,k;
         const auto& EGrid2 = e[stencil.oop()];   // i,j,k+1;
         perBGrid0[fsgrids::bfield::PERBX] += dtdz * (EGrid2[fsgrids::efield::EY] - EGrid0[fsgrids::efield::EY]) +
                                              dtdy * (EGrid0[fsgrids::efield::EZ] - EGrid1[fsgrids::efield::EZ]);
         break;
      }

      case RK_ORDER2_STEP1: {
         auto& perBDt2Grid0 = perbdt2[stencil.ooo()];       // i,j,k;
         const auto& EGrid0 = e[stencil.ooo()];             // i,j,k;
         const auto& EGrid1 = e[stencil.opo()];                 // i,j+1,k;
         const auto& EGrid2 = e[stencil.oop()];               // i,j,k+1;
         perBDt2Grid0[fsgrids::bfield::PERBX] =
             perBGrid0[fsgrids::bfield::PERBX] +
             0.5 * (dtdz * (EGrid2[fsgrids::efield::EY] - EGrid0[fsgrids::efield::EY]) +
                    dtdy * (EGrid0[fsgrids::efield::EZ] - EGrid1[fsgrids::efield::EZ]));
         break;
      }

      case RK_ORDER2_STEP2: {
         const auto& EGrid0 = edt2[stencil.ooo()]; // i,j,k;
         const auto& EGrid1 = edt2[stencil.opo()];     // i,j+1,k;
         const auto& EGrid2 = edt2[stencil.oop()];   // i,j,k+1;
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
         const auto& EGrid0 = e[stencil.ooo()]; // i,j,k;
         const auto& EGrid1 = e[stencil.oop()];   // i,j,k+1;
         const auto& EGrid2 = e[stencil.poo()];  // i+1,j,k;
         perBGrid0[fsgrids::bfield::PERBY] += dtdx * (EGrid2[fsgrids::efield::EZ] - EGrid0[fsgrids::efield::EZ]) +
                                              dtdz * (EGrid0[fsgrids::efield::EX] - EGrid1[fsgrids::efield::EX]);
         break;
      }
      case RK_ORDER2_STEP1: {
         auto& perBDt2Grid0 = perbdt2[stencil.ooo()];       // i,j,k;
         const auto& EGrid0 = e[stencil.ooo()];             // i,j,k;
         const auto& EGrid1 = e[stencil.oop()];               // i,j,k+1;
         const auto& EGrid2 = e[stencil.poo()];              // i+1,j,k;
         perBDt2Grid0[fsgrids::bfield::PERBY] =
             perBGrid0[fsgrids::bfield::PERBY] +
             0.5 * (dtdx * (EGrid2[fsgrids::efield::EZ] - EGrid0[fsgrids::efield::EZ]) +
                    dtdz * (EGrid0[fsgrids::efield::EX] - EGrid1[fsgrids::efield::EX]));
         break;
      }
      case RK_ORDER2_STEP2: {
         const auto& EGrid0 = edt2[stencil.ooo()]; // i,j,k;
         const auto& EGrid1 = edt2[stencil.oop()];   // i,j,k+1;
         const auto& EGrid2 = edt2[stencil.poo()];  // i+1,j,k;
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
         const auto& EGrid0 = e[stencil.ooo()]; // i,j,k;
         const auto& EGrid1 = e[stencil.poo()];  // i+1,j,k;
         const auto& EGrid2 = e[stencil.opo()];     // i,j+1,k;
         perBGrid0[fsgrids::bfield::PERBZ] += dtdy * (EGrid2[fsgrids::efield::EX] - EGrid0[fsgrids::efield::EX]) +
                                              dtdx * (EGrid0[fsgrids::efield::EY] - EGrid1[fsgrids::efield::EY]);
         break;
      }
      case RK_ORDER2_STEP1: {
         auto& perBDt2Grid0 = perbdt2[stencil.ooo()];       // i,j,k;
         const auto& EGrid0 = e[stencil.ooo()];             // i,j,k;
         const auto& EGrid1 = e[stencil.poo()];              // i+1,j,k;
         const auto& EGrid2 = e[stencil.opo()];                 // i,j+1,k;
         perBDt2Grid0[fsgrids::bfield::PERBZ] =
             perBGrid0[fsgrids::bfield::PERBZ] +
             0.5 * (dtdy * (EGrid2[fsgrids::efield::EX] - EGrid0[fsgrids::efield::EX]) +
                    dtdx * (EGrid0[fsgrids::efield::EY] - EGrid1[fsgrids::efield::EY]));
         break;
      }
      case RK_ORDER2_STEP2: {
         const auto& EGrid0 = edt2[stencil.ooo()]; // i,j,k;
         const auto& EGrid1 = edt2[stencil.poo()];  // i+1,j,k;
         const auto& EGrid2 = edt2[stencil.opo()];     // i,j+1,k;
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
 * \param perb fsGrid holding the perturbed B quantities at runge-kutta t=0
 * \param perbdt2 fsGrid holding the perturbed B quantities at runge-kutta t=0.5
 * \param bgb fsGrid holding the background field B quantities
 * \param technical fsgrid holding the technical parameters
 * \param gridSpacing cell size in x,y,z
 * \param globalCoordinates cell global grid coordinate indices
 * \param stencil fsgrid stencil for cell
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param component component to compute
 *
 * \sa propagateMagneticFieldSimple propagateMagneticField
 */
void propagateSysBoundaryMagneticField(fsgrids::perbspan perb,
                                       fsgrids::perbspan perbdt2,
                                       fsgrids::constbgbspan bgb,
                                       fsgrids::consttechnicalspan technical,
                                       const std::array<Real, 3>& gridSpacing,
                                       const std::array<fsgrid::FsSize_t, 3>& globalCoordinates,
                                       const fsgrid::FsStencil& stencil, SysBoundary& sysBoundaries, int32_t RKCase,
                                       uint32_t component) {
   const bool case0 = RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2;
   auto& out = case0 ? perb[stencil.ooo()] : perbdt2[stencil.ooo()];
   const auto& pb = case0 ? perb : perbdt2;

   out[fsgrids::bfield::PERBX + component] =
       sysBoundaries.getSysBoundary(technical[stencil.ooo()].sysBoundaryFlag)
           ->fieldSolverBoundaryCondMagneticField(pb, bgb, technical, gridSpacing, globalCoordinates, stencil,
                                                  component);
}

/*! \brief High-level magnetic field propagation function.
 *
 * Propagates the magnetic field and applies the field boundary conditions defined in project.h where needed.
 *
 * \param perb fsGrid holding the perturbed B quantities at runge-kutta t=0
 * \param perbdt2 fsGrid holding the perturbed B quantities at runge-kutta t=0.5
 * \param bgb fsGrid holding the background B quantities
 * \param e fsGrid holding the Electric field quantities at runge-kutta t=0
 * \param edt2 fsGrid holding the Electric field quantities at runge-kutta t=0.5
 * \param technical fsGrid holding technical information (such as boundary types)
 * \param fsgrid container of all fsgrids
 * \param sysBoundaries System boundary conditions existing
 * \param dt Length of the time step
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 *
 * \sa propagateMagneticField propagateSysBoundaryMagneticField
 */
void propagateMagneticFieldSimple(fsgrids::perbspan perb,
                                  fsgrids::perbspan perbdt2,
                                  fsgrids::bgbspan bgb,
                                  fsgrids::efieldspan e,
                                  fsgrids::efieldspan edt2,
                                  fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
                                  SysBoundary& sysBoundaries, creal& dt, cint& RKCase) {
   phiprof::Timer propagateBTimer{"Propagate magnetic field"};
   const auto* localSize = &fsgrid.getLocalSize()[0];
   const auto& gridSpacing = fsgrid.getGridSpacing();
   const size_t numCells = fsgrid.getNumCells();

   int sysBoundaryTimerId{phiprof::initializeTimer("Magnetic Field compute sysboundary cells")};
   fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                       phiprof::initializeTimer("Magnetic Field compute cells"), technical,
                       [=](const fsgrid::Coordinates &coordinates, const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
                          cuint bitfield = technical[stencil.ooo()].SOLVE;
                          propagateMagneticField(
                             perb, perbdt2, e, edt2, stencil, dt, RKCase, ((bitfield & compute::BX) == compute::BX),
                             ((bitfield & compute::BY) == compute::BY), ((bitfield & compute::BZ) == compute::BZ), coordinates.physicalGridSpacing);
                       });

   // This communication is needed for boundary conditions, in practice almost all
   // of the communication is going to be redone in calculateDerivativesSimple
   // TODO: do not transfer if there are no field boundaryconditions
   phiprof::Timer mpiTimer{"MPI", {"MPI"}};
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ with neighbours
      fsgrid.updateGhostCells(perb);
   } else { // RKCase == RK_ORDER2_STEP1
      // Exchange PERBX_DT2,PERBY_DT2,PERBZ_DT2 with neighbours
      fsgrid.updateGhostCells(perbdt2);
   }
   mpiTimer.stop();

   // Propagate B on system boundary/process inner cells
   phiprof::Timer sysBoundaryTimer {sysBoundaryTimerId};
   // L1 pass, gather which faces to solve
   std::vector<std::array<int,4>> L1Solve;
#pragma omp parallel
   {
      std::vector<std::array<int,4>> threadL1Solve;
// L1 pass
#pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const fsgrid::FsStencil stencil = fsgrid.makeStencil(i, j, k);
               const auto tech = technical[stencil.ooo()];
               cuint bitfield = tech.SOLVE;
               // L1 pass
               if (tech.sysBoundaryLayer == 1) {
                  if ((bitfield & compute::BX) != compute::BX) {
                     threadL1Solve.push_back({i,j,k,0});
                  }
                  if ((bitfield & compute::BY) != compute::BY) {
                     threadL1Solve.push_back({i,j,k,1});
                  }
                  if ((bitfield & compute::BZ) != compute::BZ) {
                     threadL1Solve.push_back({i,j,k,2});
                  }
               }
            }
         }
      }
      #pragma omp critical
      {
         L1Solve.insert(L1Solve.end(), threadL1Solve.begin(), threadL1Solve.end());
      }
      #pragma omp barrier
      //for (auto [i,j,k,dir] : L1Solve) { // not supported with OpenMP on old CLANG
      #pragma omp for // default i.e. schedule(static,1)
      for (uint entry=0; entry<L1Solve.size(); ++entry) {
         const fsgrid::FsStencil stencil = fsgrid.makeStencil(L1Solve.at(entry)[0], L1Solve.at(entry)[1], L1Solve.at(entry)[2]);
         cint component = L1Solve.at(entry)[3];
         const auto globalCoordinates = fsgrid.localToGlobal(stencil.i, stencil.j, stencil.k);
         const auto tech = technical[stencil.ooo()];
         propagateSysBoundaryMagneticField(perb, perbdt2, bgb, technical, gridSpacing, globalCoordinates, stencil, sysBoundaries, RKCase, component);
      }
   }
   sysBoundaryTimer.stop();

   mpiTimer.start();
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ with neighbours
      fsgrid.updateGhostCells(perb);
   } else { // RKCase == RK_ORDER2_STEP1
      // Exchange PERBX_DT2,PERBY_DT2,PERBZ_DT2 with neighbours
      fsgrid.updateGhostCells(perbdt2);
   }
   mpiTimer.stop();

   sysBoundaryTimer.start();
   // L2 pass, gather which faces to solve
   std::vector<std::array<int,4>> L2Solve;
   #pragma omp parallel
   {
      std::vector<std::array<int,4>> threadL2Solve;
      #pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const fsgrid::FsStencil stencil = fsgrid.makeStencil(i, j, k);
               const auto tech = technical[stencil.ooo()];
               if(tech.sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
                  tech.sysBoundaryLayer == 2
                  ) {
                  for (int component = 0; component < 3; component++) {
                     threadL2Solve.push_back({i,j,k,component});
                  }
               }
            }
         }
      }
      #pragma omp critical
      {
         L2Solve.insert(L2Solve.end(), threadL2Solve.begin(), threadL2Solve.end());
      }
      #pragma omp barrier
      //for (auto [i,j,k,dir] : L2Solve) { // not supported with OpenMP on old CLANG
      #pragma omp for // default i.e. schedule(static,1)
      for (uint entry=0; entry<L2Solve.size(); ++entry) {
         const fsgrid::FsStencil stencil = fsgrid.makeStencil(L2Solve.at(entry)[0], L2Solve.at(entry)[1], L2Solve.at(entry)[2]);
         cint component = L2Solve.at(entry)[3];
         const auto globalCoordinates = fsgrid.localToGlobal(stencil.i, stencil.j, stencil.k);
         const auto tech = technical[stencil.ooo()];
         propagateSysBoundaryMagneticField(perb, perbdt2, bgb, technical, gridSpacing, globalCoordinates, stencil, sysBoundaries, RKCase, component);
      }
   }
   propagateBTimer.stop(numCells, "Spatial Cells");
}
