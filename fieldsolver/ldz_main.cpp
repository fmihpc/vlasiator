/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

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

/*! Re-initialize field propagator after rebalance. E, BGB, RHO, RHO_V,
 cell_dimensions, sysboundaryflag need to be up to date for the
 extended neighborhood
 */
bool initializeFieldPropagatorAfterRebalance(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
) {
   // Assume static background field, they are not communicated here
   // but are assumed to be ok after each load balance as that
   // communicates all spatial data
   
   return true;
}

/*! Calculates bit masks used in the field solver and computes the initial edge electric fields from the initial magnetic fields. Then computes the initial volume averages.
 */
bool initializeFieldPropagator(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 3, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 3, 2> & BgBGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 3, 2> & volGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries
) {
   const vector<uint64_t>& localCells = getLocalCells();
   
   // Force recalculate of cell caches
   phiprof::start("Calculate Caches");
   fs_cache::calculateCache(mpiGrid,localCells);
   phiprof::stop("Calculate Caches",localCells.size(),"Spatial Cells");

   // Checking that spatial cells are cubic, otherwise field solver is incorrect (cf. derivatives in E, Hall term)
   if((abs((P::dx_ini-P::dy_ini)/P::dx_ini) > 0.001) ||
      (abs((P::dx_ini-P::dz_ini)/P::dx_ini) > 0.001) ||
      (abs((P::dy_ini-P::dz_ini)/P::dy_ini) > 0.001)) {
      std::cerr << "WARNING: Your spatial cells seem not to be cubic. However the field solver is assuming them to be. Use at your own risk and responsibility!" << std::endl;
   }
   
   // Assume static background field, they are not communicated here
   // but are assumed to be ok after each load balance as that
   // communicates all spatial data
   
   // Calculate derivatives and upwinded edge-E. Exchange derivatives 
   // and edge-E:s between neighbouring processes and calculate 
   // face-averaged E,B fields.
   bool hallTermCommunicateDerivatives = true;
   calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER1, true);
   if(P::ohmGradPeTerm > 0) {
      calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER1);
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
         localCells,
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
      sysBoundaries,
      localCells,
      RK_ORDER1
   );
   calculateVolumeAveragedFields(perBGrid,EGrid,dPerBGrid,volGrid,fs_cache::getCache().localCellsCache,fs_cache::getCache().local_NOT_DO_NOT_COMPUTE);
   calculateBVOLDerivativesSimple(volGrid, technicalGrid, sysBoundaries, localCells);
   
   return true;
}

bool finalizeFieldPropagator() {
   return true;
}

/*! \brief Top-level field propagation function.
 * 
 * Propagates the magnetic field, computes the derivatives and the upwinded electric field, then computes the volume-averaged field values. Takes care of the Runge-Kutta iteration at the top level, the functions called get as an argument the element from the enum defining the current stage and handle their job correspondingly.
 * 
 * \param mpiGrid Grid
 * \param dt Length of the time step
 * \param subcycles Number of subcycles to compute.
 * 
 * \sa propagateMagneticFieldSimple calculateDerivativesSimple calculateUpwindedElectricFieldSimple calculateVolumeAveragedFields calculateBVOLDerivativesSimple
 * 
 */
bool propagateFields(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 3, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 3, 2> & BgBGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 3, 2> & volGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& subcycles
) {
   
   if(subcycles == 0) {
      cerr << "Field solver subcycles cannot be 0." << endl;
      exit(1);
   }
   
   // Reserve memory for derivatives for all cells on this process:
   const vector<CellID>& localCells = getLocalCells();
   bool hallTermCommunicateDerivatives = true;

   if (Parameters::meshRepartitioned == true) {
      phiprof::start("Calculate Caches");
      fs_cache::calculateCache(mpiGrid,localCells);
      phiprof::stop("Calculate Caches");
   }

   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      mpiGrid[cellID]->parameters[CellParams::MAXFDT]=std::numeric_limits<Real>::max();
   }


   if (subcycles == 1) {
      #ifdef FS_1ST_ORDER_TIME
      propagateMagneticFieldSimple(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, sysBoundaries, dt, localCells, RK_ORDER1);
      calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER1, true);
      if(P::ohmGradPeTerm > 0){
         calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER1);
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
            localCells,
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
         sysBoundaries,
         localCells,
         RK_ORDER1
      );
      #else
      propagateMagneticFieldSimple(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, sysBoundaries, dt, localCells, RK_ORDER2_STEP1);
      calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER2_STEP1, true);
      if(P::ohmGradPeTerm > 0) {
         calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER2_STEP1);
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
            localCells,
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
         sysBoundaries,
         localCells,
         RK_ORDER2_STEP1
      );
      
      propagateMagneticFieldSimple(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, sysBoundaries, dt, localCells, RK_ORDER2_STEP2);
      calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER2_STEP2, true);
      if(P::ohmGradPeTerm > 0) {
         calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER2_STEP2);
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
            localCells,
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
         sysBoundaries,
         localCells,
         RK_ORDER2_STEP2
      );
      #endif
   } else {
      for (uint i=0; i<subcycles; i++) {
         propagateMagneticFieldSimple(perBGrid, perBDt2Grid, EGrid, EDt2Grid, technicalGrid, sysBoundaries, dt/convert<Real>(subcycles), localCells, RK_ORDER1);
         // If we are at the first subcycle we need to update the derivatives of the moments, 
         // otherwise only B changed and those derivatives need to be updated.
         calculateDerivativesSimple(perBGrid, perBDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER1, (i==0));
         if(P::ohmGradPeTerm > 0 && i==0) {
            calculateGradPeTermSimple(EGradPeGrid, momentsGrid, momentsDt2Grid, dMomentsGrid, sysBoundaries, localCells, RK_ORDER1);
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
               localCells,
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
            sysBoundaries,
            localCells,
            RK_ORDER1
         );
      }
   }
   
   calculateVolumeAveragedFields(perBGrid,EGrid,dPerBGrid,volGrid,fs_cache::getCache().localCellsCache,fs_cache::getCache().local_NOT_DO_NOT_COMPUTE);
   calculateBVOLDerivativesSimple(volGrid, technicalGrid, sysBoundaries, localCells);
   return true;
}
