/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute
 * 
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
 * \param cellCache Field solver cell cache
 * \param cells Vector of cells to process
 * \param dt Length of the time step
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param doX If true, compute the x component (default true).
 * \param doY If true, compute the y component (default true).
 * \param doZ If true, compute the z component (default true).
 */
void propagateMagneticField(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
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
   
   const std::array<Real, fsgrids::bfield::N_BFIELD> * perBGrid0 = perBGrid.get(i,j,k);
   
   if (doX == true) {
      switch (RKCase) {
         case RK_ORDER1:
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EGrid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1 = EGrid.get(i,j+1,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2 = EGrid.get(i,j,k+1);
            perBGrid0[fsgrids::bfield::PERBX] += dt/dz*(EGrid2[fsgrids::efield::EY] - EGrid0[fsgrids::efield::EY]) + dt/dy*(EGrid0[fsgrids::efield::EZ] - EGrid1[fsgrids::efield::EZ]);
            break;
         case RK_ORDER2_STEP1:
            const std::array<Real, fsgrids::bfield::N_BFIELD> * perBDt2Grid0 = perBDt2Grid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EGrid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1 = EGrid.get(i,j+1,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2 = EGrid.get(i,j,k+1);
            perBDt2Grid0[fsgrids::bfield::PERBX] = perBGrid0[fsgrids::bfield::PERBX] + 0.5*dt*(1.0/dz*(EGrid2[fsgrids::efield::EY] - EGrid0[fsgrids::efield::EY]) + 1.0/dy*(EGrid0[fsgrids::efield::EZ] - EGrid1[fsgrids::efield::EZ]));
            break;
         case RK_ORDER2_STEP2:
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EDt2Grid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1 = EDt2Grid.get(i,j+1,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2 = EDt2Grid.get(i,j,k+1);
            perBGrid0[fsgrids::bfield::PERBX] += dt * (1.0/dz*(EGrid2[fsgrids::efield::EY] - EGrid0[fsgrids::efield::EY]) + 1.0/dy*(EGrid0[fsgrids::efield::EZ] - EGrid1[fsgrids::efield::EZ]));
            break;
         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "Invalid RK case." << std::endl;
            abort();
      }
   }
   
   if (doY == true) {
      switch (RKCase) {
         case RK_ORDER1:
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EGrid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1 = EGrid.get(i,j,k+1);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2 = EGrid.get(i+1,j,k);
            perBGrid0[fsgrids::bfield::PERBY] += dt/dx*(EGrid2[fsgrids::efield::EZ] - EGrid0[fsgrids::efield::EZ]) + dt/dz*(EGrid0[fsgrids::efield::EX] - EGrid1[fsgrids::efield::EX]);
            break;
         case RK_ORDER2_STEP1:
            const std::array<Real, fsgrids::bfield::N_BFIELD> * perBDt2Grid0 = perBDt2Grid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EGrid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1 = EGrid.get(i,j,k+1);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2 = EGrid.get(i+1,j,k);
            perBDt2Grid0[fsgrids::bfield::PERBY] = perBGrid0[fsgrids::bfield::PERBY] + 0.5*dt*(1.0/dx*(EGrid2[fsgrids::efield::EZ] - EGrid0[fsgrids::efield::EZ]) + 1.0/dz*(EGrid0[fsgrids::efield::EX] - EGrid1[fsgrids::efield::EX]));
            break;
         case RK_ORDER2_STEP2:
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EDt2Grid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1 = EDt2Grid.get(i,j,k+1);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2 = EDt2Grid.get(i+1,j,k);
            perBGrid0[fsgrids::bfield::PERBY] += dt * (1.0/dx*(EGrid2[fsgrids::efield::EZ] - EGrid0[fsgrids::efield::EZ]) + 1.0/dz*(EGrid0[fsgrids::efield::EX] - EGrid1[fsgrids::efield::EX]));
            break;
         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "Invalid RK case." << std::endl;
            abort();
      }
   }
   
   if (doZ == true) {
      switch (RKCase) {
         case RK_ORDER1:
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EGrid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1 = EGrid.get(i+1,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2 = EGrid.get(i,j+1,k);
            perBGrid0[fsgrids::bfield::PERBZ] += dt/dy*(EGrid2[fsgrids::efield::EX] - EGrid0[fsgrids::efield::EX]) + dt/dx*(EGrid0[fsgrids::efield::EY] - EGrid1[fsgrids::efield::EY]);
            break;
         case RK_ORDER2_STEP1:
            const std::array<Real, fsgrids::bfield::N_BFIELD> * perBDt2Grid0 = perBDt2Grid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EGrid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1 = EGrid.get(i+1,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2 = EGrid.get(i,j+1,k);
            perBDt2Grid0[fsgrids::bfield::PERBZ] = perBGrid0[fsgrids::bfield::PERBZ] + 0.5*dt*(1.0/dy*(EGrid2[fsgrids::efield::EX] - EGrid0[fsgrids::efield::EX]) + 1.0/dx*(EGrid0[fsgrids::efield::EY] - EGrid1[fsgrids::efield::EY]));
            break;
         case RK_ORDER2_STEP2:
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EDt2Grid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid0 = EDt2Grid.get(i,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid1 = EDt2Grid.get(i+1,j,k);
            const std::array<Real, fsgrids::efield::N_EFIELD> * EGrid2 = EDt2Grid.get(i,j+1,k);
            perBGrid0[fsgrids::bfield::PERBZ] += dt  * (1.0/dy*(EGrid2[fsgrids::efield::EX] - EGrid0[fsgrids::efield::EX]) + 1.0/dx*(EGrid0[fsgrids::efield::EY] - EGrid1[fsgrids::efield::EY]));
            break;
         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "Invalid RK case." << std::endl;
            abort();
      }
   }
}

/*! \brief High-level magnetic field propagation function.
 * 
 * Propagates the magnetic field and applies the field boundary conditions defined in project.h where needed.
 * 
 * \param mpiGrid Grid
 * \param sysBoundaries System boundary conditions existing
 * \param dt Length of the time step
 * \param localCells Vector of local cells to process
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa propagateMagneticField propagateSysBoundaryMagneticField
 */
void propagateMagneticFieldSimple(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& RKCase
) {
   int timer;
   const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   
   phiprof::start("Propagate magnetic field");
   
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);
   
   #pragma omp parallel for collapse(3)
   for (uint k=0; k<gridDims[2]; k++) {
      for (uint j=0; j<gridDims[1]; j++) {
         for (uint i=0; i<gridDims[0]; i++) {
            if(technicalGrid.get(i,j,k)->sysBoundaryFlag != NOT_SYSBOUNDARY) continue;
            // Propagate B on all local cells:
            propagateMagneticField(perBGrid, perBDt2Grid, EGrid, EDt2Grid, i, j, k, dt, RKCase);
         }
      }
   }
   
   //phiprof::stop("propagate not sysbound",localCells.size(),"Spatial Cells");
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   //This communication is needed for boundary conditions, in practice almost all
   //of the communication is going to be redone in calculateDerivativesSimple
   //TODO: do not transfer if there are no field boundaryconditions
   phiprof::start("MPI");
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ with neighbours
      perBGrid.updateGhostCells();
   } else { // RKCase == RK_ORDER2_STEP1
      // Exchange PERBX_DT2,PERBY_DT2,PERBZ_DT2 with neighbours
      perBDt2Grid.updateGhostCells();
   }
   
   phiprof::stop("MPI");
   
   // Propagate B on system boundary/process inner cells
   timer=phiprof::initializeTimer("Compute system boundary cells");
   phiprof::start(timer);
   #pragma omp parallel for collapse(3)
   for (uint k=0; k<gridDims[2]; k++) {
      for (uint j=0; j<gridDims[1]; j++) {
         for (uint i=0; i<gridDims[0]; i++) {
            if(technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundaries::NOT_SYSBOUNDARY) {
               propagateSysBoundaryMagneticField(perBGrid, perBDt2Grid, technicalGrid, i, j, k, sysBoundaries, dt, RKCase);
            }
         }
      }
   }
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   phiprof::stop("Propagate magnetic field",N_cells,"Spatial Cells");
}

/*! \brief Low-level magnetic field propagation function.
 * 
 * Propagates the magnetic field according to the system boundary conditions.
 * 
 * \param mpiGrid Grid
 * \param cellCache Field solver cell cache
 * \param localID Field solver cache local cell ID
 * \param sysBoundaries System boundary conditions existing
 * \param dt Length of the time step
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa propagateMagneticFieldSimple propagateMagneticField
 */
void propagateSysBoundaryMagneticField(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   cint i,
   cint j,
   cint k,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& RKCase
) {
   const FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> * bGrid;
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      bGrid = & perBGrid.get(i,j,k);
   } else {
      bGrid = & perBDt2Grid.get(i,j,k);
   }
   cuint sysBoundaryFlag = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   for (uint component = 0; component < 3; component++) {
      bGrid[fsgrids::bfield::PERBX + component] = sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondMagneticField(perBGrid, perBDt2Grid, technicalGrid, i, j, k, dt, RKCase, component);
   }
}
