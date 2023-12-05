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
void propagateMagneticField(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EDt2Grid,
   cint i,
   cint j,
   cint k,
   creal& dt,
   cint& RKCase,
   const bool doX, //=true (default)
   const bool doY, //=true (default)
   const bool doZ  //=true (default)
) {
   creal dx = perBGrid.DX;
   creal dy = perBGrid.DY;
   creal dz = perBGrid.DZ;

   std::array<Real, fsgrids::bfield::N_BFIELD> * perBGrid0 = perBGrid.get(i,j,k);
   std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0;
   std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1;
   std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2;
   std::array<Real, fsgrids::bfield::N_BFIELD> * perBDt2Grid0;

   if (doX == true) {
      switch (RKCase) {
         case RK_ORDER1:
            EGrid0 = EGrid.get(i,j,k);
            EGrid1 = EGrid.get(i,j+1,k);
            EGrid2 = EGrid.get(i,j,k+1);
            perBGrid0->at(fsgrids::bfield::PERBX) += dt/dz*(EGrid2->at(fsgrids::efield::EY) - EGrid0->at(fsgrids::efield::EY)) + dt/dy*(EGrid0->at(fsgrids::efield::EZ) - EGrid1->at(fsgrids::efield::EZ));
            break;

         case RK_ORDER2_STEP1:
            perBDt2Grid0 = perBDt2Grid.get(i,j,k);
            EGrid0 = EGrid.get(i,j,k);
            EGrid1 = EGrid.get(i,j+1,k);
            EGrid2 = EGrid.get(i,j,k+1);
            perBDt2Grid0->at(fsgrids::bfield::PERBX) = perBGrid0->at(fsgrids::bfield::PERBX) + 0.5*dt*(1.0/dz*(EGrid2->at(fsgrids::efield::EY) - EGrid0->at(fsgrids::efield::EY)) + 1.0/dy*(EGrid0->at(fsgrids::efield::EZ) - EGrid1->at(fsgrids::efield::EZ)));
            break;

         case RK_ORDER2_STEP2:
            EGrid0 = EDt2Grid.get(i,j,k);
            EGrid1 = EDt2Grid.get(i,j+1,k);
            EGrid2 = EDt2Grid.get(i,j,k+1);
            perBGrid0->at(fsgrids::bfield::PERBX) += dt * (1.0/dz*(EGrid2->at(fsgrids::efield::EY) - EGrid0->at(fsgrids::efield::EY)) + 1.0/dy*(EGrid0->at(fsgrids::efield::EZ) - EGrid1->at(fsgrids::efield::EZ)));
            break;

         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "Invalid RK case." << std::endl;
            abort();
      }
   }

   if (doY == true) {
      switch (RKCase) {
         case RK_ORDER1:
            EGrid0 = EGrid.get(i,j,k);
            EGrid1 = EGrid.get(i,j,k+1);
            EGrid2 = EGrid.get(i+1,j,k);
            perBGrid0->at(fsgrids::bfield::PERBY) += dt/dx*(EGrid2->at(fsgrids::efield::EZ) - EGrid0->at(fsgrids::efield::EZ)) + dt/dz*(EGrid0->at(fsgrids::efield::EX) - EGrid1->at(fsgrids::efield::EX));
            break;
         case RK_ORDER2_STEP1:
            perBDt2Grid0 = perBDt2Grid.get(i,j,k);
            EGrid0 = EGrid.get(i,j,k);
            EGrid1 = EGrid.get(i,j,k+1);
            EGrid2 = EGrid.get(i+1,j,k);
            perBDt2Grid0->at(fsgrids::bfield::PERBY) = perBGrid0->at(fsgrids::bfield::PERBY) + 0.5*dt*(1.0/dx*(EGrid2->at(fsgrids::efield::EZ) - EGrid0->at(fsgrids::efield::EZ)) + 1.0/dz*(EGrid0->at(fsgrids::efield::EX) - EGrid1->at(fsgrids::efield::EX)));
            break;
         case RK_ORDER2_STEP2:
            EGrid0 = EDt2Grid.get(i,j,k);
            EGrid1 = EDt2Grid.get(i,j,k+1);
            EGrid2 = EDt2Grid.get(i+1,j,k);
            perBGrid0->at(fsgrids::bfield::PERBY) += dt * (1.0/dx*(EGrid2->at(fsgrids::efield::EZ) - EGrid0->at(fsgrids::efield::EZ)) + 1.0/dz*(EGrid0->at(fsgrids::efield::EX) - EGrid1->at(fsgrids::efield::EX)));
            break;
         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "Invalid RK case." << std::endl;
            abort();
      }
   }

   if (doZ == true) {
      switch (RKCase) {
         case RK_ORDER1:
            EGrid0 = EGrid.get(i,j,k);
            EGrid1 = EGrid.get(i+1,j,k);
            EGrid2 = EGrid.get(i,j+1,k);
            perBGrid0->at(fsgrids::bfield::PERBZ) += dt/dy*(EGrid2->at(fsgrids::efield::EX) - EGrid0->at(fsgrids::efield::EX)) + dt/dx*(EGrid0->at(fsgrids::efield::EY) - EGrid1->at(fsgrids::efield::EY));
            break;
         case RK_ORDER2_STEP1:
            perBDt2Grid0 = perBDt2Grid.get(i,j,k);
            EGrid0 = EGrid.get(i,j,k);
            EGrid1 = EGrid.get(i+1,j,k);
            EGrid2 = EGrid.get(i,j+1,k);
            perBDt2Grid0->at(fsgrids::bfield::PERBZ) = perBGrid0->at(fsgrids::bfield::PERBZ) + 0.5*dt*(1.0/dy*(EGrid2->at(fsgrids::efield::EX) - EGrid0->at(fsgrids::efield::EX)) + 1.0/dx*(EGrid0->at(fsgrids::efield::EY) - EGrid1->at(fsgrids::efield::EY)));
            break;
         case RK_ORDER2_STEP2:
            EGrid0 = EDt2Grid.get(i,j,k);
            EGrid1 = EDt2Grid.get(i+1,j,k);
            EGrid2 = EDt2Grid.get(i,j+1,k);
            perBGrid0->at(fsgrids::bfield::PERBZ) += dt  * (1.0/dy*(EGrid2->at(fsgrids::efield::EX) - EGrid0->at(fsgrids::efield::EX)) + 1.0/dx*(EGrid0->at(fsgrids::efield::EY) - EGrid1->at(fsgrids::efield::EY)));
            break;
         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "Invalid RK case." << std::endl;
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
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EDt2Grid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   cint i,
   cint j,
   cint k,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& RKCase,
   cuint component
) {
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBX + component) = sysBoundaries.getSysBoundary(technicalGrid.get(i,j,k)->sysBoundaryFlag)->fieldSolverBoundaryCondMagneticField(perBGrid, technicalGrid, i, j, k, dt, component);
   } else {
      perBDt2Grid.get(i,j,k)->at(fsgrids::bfield::PERBX + component) = sysBoundaries.getSysBoundary(technicalGrid.get(i,j,k)->sysBoundaryFlag)->fieldSolverBoundaryCondMagneticField(perBDt2Grid, technicalGrid, i, j, k, dt, component);
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
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EDt2Grid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& RKCase
) {
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const int* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::Timer propagateBTimer {"Propagate magnetic field"};

   int computeTimerId {phiprof::initializeTimer("Magnetic Field compute cells")};
   int sysBoundaryTimerId {phiprof::initializeTimer("Magnetic Field compute sysboundary cells")};
   #pragma omp parallel
   {
      phiprof::Timer computeTimer {computeTimerId};
      #pragma omp for collapse(2)
      for (int k=0; k<gridDims[2]; k++) {
         for (int j=0; j<gridDims[1]; j++) {
            for (int i=0; i<gridDims[0]; i++) {
               cuint bitfield = technicalGrid.get(i,j,k)->SOLVE;
               propagateMagneticField(perBGrid, perBDt2Grid, EGrid, EDt2Grid, i, j, k, dt, RKCase, ((bitfield & compute::BX) == compute::BX), ((bitfield & compute::BY) == compute::BY), ((bitfield & compute::BZ) == compute::BZ));
            }
         }
      }
   }

   //This communication is needed for boundary conditions, in practice almost all
   //of the communication is going to be redone in calculateDerivativesSimple
   //TODO: do not transfer if there are no field boundaryconditions
   phiprof::Timer mpiTimer {"MPI", {"MPI"}};
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
      phiprof::Timer sysBoundaryTimer {sysBoundaryTimerId};
      // L1 pass
      #pragma omp for collapse(2)
      for (int k=0; k<gridDims[2]; k++) {
         for (int j=0; j<gridDims[1]; j++) {
            for (int i=0; i<gridDims[0]; i++) {
               cuint bitfield = technicalGrid.get(i,j,k)->SOLVE;
               // L1 pass
               if (technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) {
                  if ((bitfield & compute::BX) != compute::BX) {
                     propagateSysBoundaryMagneticField(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, i, j, k, sysBoundaries, dt, RKCase, 0);
                  }
                  if ((bitfield & compute::BY) != compute::BY) {
                     propagateSysBoundaryMagneticField(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, i, j, k, sysBoundaries, dt, RKCase, 1);
                  }
                  if ((bitfield & compute::BZ) != compute::BZ) {
                     propagateSysBoundaryMagneticField(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, i, j, k, sysBoundaries, dt, RKCase, 2);
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
      phiprof::Timer sysBoundaryTimer {sysBoundaryTimerId};
      // L2 pass
      #pragma omp for collapse(2)
      for (int k=0; k<gridDims[2]; k++) {
         for (int j=0; j<gridDims[1]; j++) {
            for (int i=0; i<gridDims[0]; i++) {
               if(technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
                  technicalGrid.get(i,j,k)->sysBoundaryLayer == 2
                  ) {
                  for (uint component = 0; component < 3; component++) {
                     propagateSysBoundaryMagneticField(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, i, j, k, sysBoundaries, dt, RKCase, component);
                  }
               }
            }
         }
      }
   }
   propagateBTimer.stop(N_cells,"Spatial Cells");
}
