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

/*! \file ldz_main.cpp
 * \brief Londrillo -- Del Zanna upwind constrained transport field solver.
 * 
 * On the divergence-free condition in Godunov-type schemes for
 * ideal magnetohydrodynamics: the upwind constrained transport method,
 * P. Londrillo and L. Del Zanna, J. Comp. Phys., 195, 2004.
 * http://dx.doi.org/10.1016/j.jcp.2003.09.016
 *
 * Reconstructions taken from:
 * Efficient, high accuracy ADER-WENO schemes for hydrodynamics and
 * divergence-free magnetohydrodynamics, D. S. Balsara, T. Rumpf,
 * M. Dumbser, C.-D. Munz, J. Comp. Phys, 228, 2480-2516, 2009.
 * http://dx.doi.org/10.1016/j.jcp.2008.12.003
 * and
 * Divergence-free reconstruction of magnetic fields and WENO
 * schemes for magnetohydrodynamics, D. S. Balsara, J. Comp. Phys.,
 * 228, 5040-5056, 2009.
 * http://dx.doi.org/10.1016/j.jcp.2009.03.038
 * 
 * *****  NOTATION USED FOR VARIABLES FOLLOWS THE ONES USED  *****\n
 * *****      IN THE ABOVEMENTIONED PUBLICATION(S)           *****
 */

#include "ldz_electric_field.hpp"
#include "ldz_magnetic_field.hpp"
#include "ldz_hall.hpp"
#include "ldz_gradpe.hpp"
#include "ldz_volume.hpp"
#include "fs_common.h"
#include "derivatives.hpp"
#include "fs_limiters.h"
#include "mpiconversion.h"

/*! Re-initialize field propagator after rebalance. E, BGB, RHO, RHO_V,
 cell_dimensions, sysboundaryflag need to be up to date for the
 extended neighborhood
 */
bool initializeFieldPropagatorAfterRebalance() {
   // Assume static background field, they are not communicated here
   // but are assumed to be ok after each load balance as that
   // communicates all spatial data
   
   return true;
}

bool finalizeFieldPropagator() {
   return true;
}

/*! \brief Top-level field propagation function.
 * 
 * Propagates the magnetic field, computes the derivatives and the upwinded
 * electric field, then computes the volume-averaged field values. Takes care
 * of the Runge-Kutta iteration at the top level, the functions called get as
 * an argument the element from the enum defining the current stage and handle
 * their job correspondingly.
 * 
 * \param dt Length of the time step
 * \param subcycles Number of subcycles to compute.
 * 
 * \sa propagateMagneticFieldSimple calculateDerivativesSimple calculateUpwindedElectricFieldSimple calculateVolumeAveragedFields calculateBVOLDerivativesSimple
 * 
 */
bool propagateFields(
    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, 2> &perBGrid,
    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, 2> &perBDt2Grid,
    FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, 2> &EGrid,
    FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, 2> &EDt2Grid,
    FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, 2> &EHallGrid,
    FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> &EGradPeGrid,
    FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, 2> &momentsGrid,
    FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, 2> &momentsDt2Grid,
    FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, 2> &dPerBGrid,
    FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> &dMomentsGrid,
    FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, 2> &BgBGrid,
    FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, 2> &volGrid,
    FsGrid<fsgrids::technical, 2> &technicalGrid,
    FsGrid<std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,
    SysBoundary &sysBoundaries,
    creal &dt,
    cuint subcycles)
{

   if(subcycles == 0) {
      cerr << "Field solver subcycles cannot be 0." << endl;
      exit(1);
   }
   
   bool hallTermCommunicateDerivatives = true;
   
   const int* gridDims = &technicalGrid.getLocalSize()[0];
   
   #pragma omp parallel for collapse(3)
   for (int k=0; k<gridDims[2]; k++) {
      for (int j=0; j<gridDims[1]; j++) {
         for (int i=0; i<gridDims[0]; i++) {
            technicalGrid.get(i,j,k)->maxFsDt=std::numeric_limits<Real>::max();
         }
      }
   }
   
   
   if (subcycles == 1) {
      #ifdef FS_1ST_ORDER_TIME
      propagateMagneticFieldSimple(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, sysBoundaries, pmlGrid, dt, RK_ORDER1);
      calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER1, true);
      if(P::ohmGradPeTerm > 0){
         calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER1);
         hallTermCommunicateDerivatives = false;
      }
      if(P::ohmHallTerm > 0) {
         calculateHallTermSimple(
            perBGrid,
            perBDt2Grid,
            EHallGrid,
            momentsGrid,
            momentsDt2Grid,
            dPerBGrid,
            dMomentsGrid,
            BgBGrid,
            technicalGrid,
            sysBoundaries,
            RK_ORDER1,
            hallTermCommunicateDerivatives
         );
      }
      calculateUpwindedElectricFieldSimple(
         perBGrid,
         perBDt2Grid,
         EGrid,
         EDt2Grid,
         EHallGrid,
         EGradPeGrid,
         momentsGrid,
         momentsDt2Grid,
         dPerBGrid,
         dMomentsGrid,
         BgBGrid,
         technicalGrid,
         pmlGrid,
         sysBoundaries,
         RK_ORDER1
      );
      #else
      propagateMagneticFieldSimple(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, sysBoundaries,pmlGrid, dt, RK_ORDER2_STEP1);
      calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER2_STEP1, true);
      if(P::ohmGradPeTerm > 0) {
         calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER2_STEP1);
         hallTermCommunicateDerivatives = false;
      }
      if(P::ohmHallTerm > 0) {
         calculateHallTermSimple(
            perBGrid,
            perBDt2Grid,
            EHallGrid,
            momentsGrid,
            momentsDt2Grid,
            dPerBGrid,
            dMomentsGrid,
            BgBGrid,
            technicalGrid,
            sysBoundaries,
            RK_ORDER2_STEP1,
            hallTermCommunicateDerivatives
         );
      }
      calculateUpwindedElectricFieldSimple(
         perBGrid,
         perBDt2Grid,
         EGrid,
         EDt2Grid,
         EHallGrid,
         EGradPeGrid,
         momentsGrid,
         momentsDt2Grid,
         dPerBGrid,
         dMomentsGrid,
         BgBGrid,
         technicalGrid,
         pmlGrid,
         sysBoundaries,
         RK_ORDER2_STEP1
      );
      
      propagateMagneticFieldSimple(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, sysBoundaries,pmlGrid, dt, RK_ORDER2_STEP2);
      calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER2_STEP2, true);
      if(P::ohmGradPeTerm > 0) {
         calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER2_STEP2);
         hallTermCommunicateDerivatives = false;
      }
      if(P::ohmHallTerm > 0) {
         calculateHallTermSimple(
            perBGrid,
            perBDt2Grid,
            EHallGrid,
            momentsGrid,
            momentsDt2Grid,
            dPerBGrid,
            dMomentsGrid,
            BgBGrid,
            technicalGrid,
            sysBoundaries,
            RK_ORDER2_STEP2,
            hallTermCommunicateDerivatives
         );
      }
      calculateUpwindedElectricFieldSimple(
         perBGrid,
         perBDt2Grid,
         EGrid,
         EDt2Grid,
         EHallGrid,
         EGradPeGrid,
         momentsGrid,
         momentsDt2Grid,
         dPerBGrid,
         dMomentsGrid,
         BgBGrid,
         technicalGrid,
         pmlGrid,
         sysBoundaries,
         RK_ORDER2_STEP2
      );
      #endif
   } else {
      Real subcycleDt = dt/convert<Real>(subcycles);
      Real subcycleT = P::t;
      creal targetT = P::t + dt;
      uint subcycleCount = 0;
      uint maxSubcycleCount = std::numeric_limits<uint>::max();
      int myRank = perBGrid.getRank();
      
      while (subcycleCount < maxSubcycleCount ) {         
         // In case of subcycling, we decided to go for a blunt Runge-Kutta subcycling even though e.g. moments are not going along.
         // Result of the Summer of Debugging 2016, the behaviour in wave dispersion was much improved with this.
         propagateMagneticFieldSimple(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, sysBoundaries,pmlGrid, subcycleDt, RK_ORDER2_STEP1);
         
         // We need to calculate derivatives of the moments at every substep, but they only
         // need to be communicated in the first one.
         calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER2_STEP1, (subcycleCount==0));
         if(P::ohmGradPeTerm > 0 && subcycleCount==0) {
            calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER2_STEP1);
            hallTermCommunicateDerivatives = false;
         }
         if(P::ohmHallTerm > 0) {
            calculateHallTermSimple(
               perBGrid,
               perBDt2Grid,
               EHallGrid,
               momentsGrid,
               momentsDt2Grid,
               dPerBGrid,
               dMomentsGrid,
               BgBGrid,
               technicalGrid,
               sysBoundaries,
               RK_ORDER2_STEP1,
               hallTermCommunicateDerivatives
            );
         }
         calculateUpwindedElectricFieldSimple(
            perBGrid,
            perBDt2Grid,
            EGrid,
            EDt2Grid,
            EHallGrid,
            EGradPeGrid,
            momentsGrid,
            momentsDt2Grid,
            dPerBGrid,
            dMomentsGrid,
            BgBGrid,
            technicalGrid,
            pmlGrid,
            sysBoundaries,
            RK_ORDER2_STEP1
         );
         
         propagateMagneticFieldSimple(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, sysBoundaries,pmlGrid, subcycleDt, RK_ORDER2_STEP2);
         
         // We need to calculate derivatives of the moments at every substep, but they only
         // need to be communicated in the first one.
         calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER2_STEP2, (subcycleCount==0));
         if(P::ohmGradPeTerm > 0 && subcycleCount==0) {
            calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, technicalGrid, sysBoundaries, RK_ORDER2_STEP2);
            hallTermCommunicateDerivatives = false;
         }
         if(P::ohmHallTerm > 0) {
            calculateHallTermSimple(
               perBGrid,
               perBDt2Grid,
               EHallGrid,
               momentsGrid,
               momentsDt2Grid,
               dPerBGrid,
               dMomentsGrid,
               BgBGrid,
               technicalGrid,
               sysBoundaries,
               RK_ORDER2_STEP2,
               hallTermCommunicateDerivatives
            );
         }
         calculateUpwindedElectricFieldSimple(
            perBGrid,
            perBDt2Grid,
            EGrid,
            EDt2Grid,
            EHallGrid,
            EGradPeGrid,
            momentsGrid,
            momentsDt2Grid,
            dPerBGrid,
            dMomentsGrid,
            BgBGrid,
            technicalGrid,
            pmlGrid,
            sysBoundaries,
            RK_ORDER2_STEP2
         );
         
         phiprof::start("FS subcycle stuff");
         subcycleT += subcycleDt; 
         subcycleCount++;

         if( subcycleT >= targetT || subcycleCount >= maxSubcycleCount  ) {
            //we are done
            if( subcycleT > targetT ) {
               //due to roundoff we might hit this, should add delta
               std::cerr << "subcycleT > targetT, should not happen! (values: subcycleT " << subcycleT << ", subcycleDt " << subcycleDt << ", targetT " << targetT << ")" << std::endl;
            }

            // Make sure the phiprof group is closed when leaving the loop
            phiprof::stop("FS subcycle stuff");

            break;
         }


         
         // Reassess subcycle dt
         Real dtMaxLocal;
         Real dtMaxGlobal;
         dtMaxLocal=std::numeric_limits<Real>::max();

         std::array<int32_t, 3>& localSize = technicalGrid.getLocalSize();
         for(int z=0; z<localSize[2]; z++) {
            for(int y=0; y<localSize[1]; y++) {
               for(int x=0; x<localSize[0]; x++) {
                  fsgrids::technical* cell = technicalGrid.get(x,y,z);
                  if ( cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
                        (cell->sysBoundaryLayer == 1 && cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY )) {
                     dtMaxLocal=min(dtMaxLocal, cell->maxFsDt);
                  }
               }
            }
         }

         phiprof::start("MPI_Allreduce");
         technicalGrid.Allreduce(&(dtMaxLocal), &(dtMaxGlobal), 1, MPI_Type<Real>(), MPI_MIN);
         phiprof::stop("MPI_Allreduce");
         
         //reduce dt if it is too high
         if( subcycleDt > dtMaxGlobal * P::fieldSolverMaxCFL ) {
            creal meanFieldsCFL = 0.5*(P::fieldSolverMaxCFL+ P::fieldSolverMinCFL);
            subcycleDt = meanFieldsCFL * dtMaxGlobal;
            if ( myRank == MASTER_RANK ) {
               logFile << "(TIMESTEP) New field solver subcycle dt = " << subcycleDt << " computed on step " <<  P::tstep << " and substep " << subcycleCount << " at " << P::t << " s" << std::endl;
            }
         }

         // Readjust the dt to hit targetT. Try to avoid having a very
         // short delta step at the end, instead 2 more normal ones
         if( subcycleT + 1.5 * subcycleDt  > targetT ) {
            subcycleDt = targetT - subcycleT;
            maxSubcycleCount = subcycleCount + 1; // 1 more steps
            //check that subcyclDt has correct CFL, take 2 if not
            if(subcycleDt > dtMaxGlobal * P::fieldSolverMaxCFL ) {
               subcycleDt = (targetT - subcycleT)/2;
               maxSubcycleCount = subcycleCount + 2; 
            }
         }
         
         phiprof::stop("FS subcycle stuff");
      }
      
      
      if( subcycles != subcycleCount && myRank == MASTER_RANK) {
         logFile << "Effective field solver subcycles were " << subcycleCount << " instead of " << P::fieldSolverSubcycles << " on step " <<  P::tstep << std::endl;
      }
   }
   
   calculateVolumeAveragedFields(perBGrid,EGrid,dPerBGrid,volGrid,technicalGrid);
   calculateBVOLDerivativesSimple(volGrid, technicalGrid, sysBoundaries);
   return true;
}
