/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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

#include "fs_cache.h"
#include "ldz_electric_field.hpp"
#include "ldz_magnetic_field.hpp"
#include "ldz_hall.hpp"
#include "ldz_gradpe.hpp"
#include "ldz_volume.hpp"
#include "fs_common.h"
#include "derivatives.hpp"
#include "fs_limiters.h"
#include "mpiconversion.h"

extern map<CellID,uint> existingCellsFlags; /**< Defined in fs_common.cpp */

extern Logger logFile; //, diagnostic; can be used also later

void calculateExistingCellsFlags(
                                 dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                 const vector<CellID>& localCells
) {
   
   existingCellsFlags.clear();
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      
      // Raise the bit for each existing cell within a 3x3 cube of 
      // spatial cells. This cell sits at the center of the cube.
      uint existingCellsFlag = (1 << calcNbrNumber(1,1,1)); // The cell itself exists (bit 13 set to 1)
      
      for (int k=-1; k<2; ++k)
         for (int j=-1; j<2; ++j)
            for (int i=-1; i<2; ++i) {
               if (i == 0 && (j == 0 && k == 0)) continue;
               const CellID nbr = getNeighbourID(mpiGrid, cellID, 2 + i, 2 + j, 2 + k);
               if (nbr == INVALID_CELLID)
                  continue;
               if (mpiGrid[nbr]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE)
                  continue;
               existingCellsFlag = existingCellsFlag | (1 << calcNbrNumber(1+i,1+j,1+k));
            }
      existingCellsFlags[cellID] = existingCellsFlag;
   }
}

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
   
   const vector<uint64_t>& localCells = getLocalCells();
   //calculateSysBoundaryFlags(mpiGrid,localCells);
   calculateExistingCellsFlags(mpiGrid,localCells);
   return true;
}

/*! Calculates bit masks used in the field solver and computes the initial edge electric fields from the initial magnetic fields. Then computes the initial volume averages.
 */
bool initializeFieldPropagator(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
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

   //calculateSysBoundaryFlags(mpiGrid,localCells);
   calculateExistingCellsFlags(mpiGrid,localCells);

   // Calculate bit masks used for if-statements by field propagator. 
   // These are used to test whether or not certain combination of 
   // neighbours exists for a cell. These can be replaced by honest 
   // if-statements, but you will just end up needing very many of them 
   // as each bit mask tests the existence of several neighbours at once.
   // Existence of neighbours would also need to be queried from the 
   // parallel grid, i.e. using if-statements is likely to be much 
   // slower.

   // x-derivatives are calculated if x-face neighbours exist:
   CALCULATE_DX = 0;
   CALCULATE_DX = CALCULATE_DX | (1 << calcNbrNumber(0,1,1));
   CALCULATE_DX = CALCULATE_DX | (1 << calcNbrNumber(2,1,1));
   
   // y-derivatives are calculated if y-face neighbours exist:
   CALCULATE_DY = 0;
   CALCULATE_DY = CALCULATE_DY | (1 << calcNbrNumber(1,0,1));
   CALCULATE_DY = CALCULATE_DY | (1 << calcNbrNumber(1,2,1));
   
   // z-derivatives are calculated if z-face neighbours exist:
   CALCULATE_DZ = 0;
   CALCULATE_DZ = CALCULATE_DZ | (1 << calcNbrNumber(1,1,0));
   CALCULATE_DZ = CALCULATE_DZ | (1 << calcNbrNumber(1,1,2));
   
   // xy mixed derivatives are calculated if +/-x,+/-y diagonal neighbours exist
   CALCULATE_DXY = 0;
   if(P::ohmHallTerm > 1) {
      CALCULATE_DXY = CALCULATE_DXY | (1 << calcNbrNumber(0,0,1));
      CALCULATE_DXY = CALCULATE_DXY | (1 << calcNbrNumber(2,2,1));
      CALCULATE_DXY = CALCULATE_DXY | (1 << calcNbrNumber(0,2,1));
      CALCULATE_DXY = CALCULATE_DXY | (1 << calcNbrNumber(2,0,1));
   }
   
   // xz mixed derivatives are calculated if +/-x,+/-z diagonal neighbours exist
   CALCULATE_DXZ = 0;
   if(P::ohmHallTerm > 1) {
      CALCULATE_DXZ = CALCULATE_DXZ | (1 << calcNbrNumber(0,1,0));
      CALCULATE_DXZ = CALCULATE_DXZ | (1 << calcNbrNumber(2,1,2));
      CALCULATE_DXZ = CALCULATE_DXZ | (1 << calcNbrNumber(0,1,2));
      CALCULATE_DXZ = CALCULATE_DXZ | (1 << calcNbrNumber(2,1,0));
   }
   // yz mixed derivatives are calculated if +/-y,+/-z diagonal neighbours exist
   CALCULATE_DYZ = 0;
   if(P::ohmHallTerm > 1) {
      CALCULATE_DYZ = CALCULATE_DYZ | (1 << calcNbrNumber(1,0,0));
      CALCULATE_DYZ = CALCULATE_DYZ | (1 << calcNbrNumber(1,2,2));
      CALCULATE_DYZ = CALCULATE_DYZ | (1 << calcNbrNumber(1,0,2));
      CALCULATE_DYZ = CALCULATE_DYZ | (1 << calcNbrNumber(1,2,0));
   }
   
   // Edge Ex is calculated if -y,-z,+/-x neighbours exist:
   CALCULATE_EX = 0;
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(1,0,1)); // -y nbr
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(1,1,0)); // -z nbr
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(0,1,1)); // -x nbr
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(2,1,1)); // +x nbr
      
   // Edge Ey is calculated if -x,-z,+/-y neighbours exist:
   CALCULATE_EY = 0;
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(0,1,1)); // -x nbr
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(1,1,0)); // -z nbr
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(1,0,1)); // -y nbr
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(1,2,1)); // +y nbr
      
   // Edge Ez is calculated      if -x,-y,+/-z neighbours exist:
   CALCULATE_EZ = 0;
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(0,1,1)); // -x nbr
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(1,0,1)); // -y nbr
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(1,1,0)); // -z nbr
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(1,1,2)); // +z nbr

   // Bx is propagated if -x,+/-y,+/-z neighbours exist:
   PROPAGATE_BX = 0;
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(0,1,1)); // -x nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,0,1)); // -y nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,2,1)); // +y nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,1,0)); // -z nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,1,2)); // +z nbr
   
   // By is propagated if -y,+/-x,+/-z neighbours exist:
   PROPAGATE_BY = 0;
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(1,0,1)); // -y nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(0,1,1)); // -x nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(2,1,1)); // +x nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(1,1,0)); // -z nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(1,1,2)); // +z nbr
   
   // Bz is propagated if -z,+/-x,+/-y neighbours exist:
   PROPAGATE_BZ = 0;
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(1,1,0)); // -z nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(0,1,1)); // -x nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(2,1,1)); // +x nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(1,0,1)); // -y nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(1,2,1)); // +y nbr
   
   // Assume static background field, they are not communicated here
   // but are assumed to be ok after each load balance as that
   // communicates all spatial data
   
   // Calculate derivatives and upwinded edge-E. Exchange derivatives 
   // and edge-E:s between neighbouring processes and calculate 
   // face-averaged E,B fields.
   bool hallTermCommunicateDerivatives = true;
   calculateDerivativesSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER1, true);
   if(P::ohmGradPeTerm > 0) {
      calculateGradPeTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER1);
      hallTermCommunicateDerivatives = false;
   }
   if(P::ohmHallTerm > 0) {
      calculateHallTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER1, hallTermCommunicateDerivatives);
   }
   calculateUpwindedElectricFieldSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER1);
   calculateVolumeAveragedFields(mpiGrid,fs_cache::getCache().localCellsCache,fs_cache::getCache().local_NOT_DO_NOT_COMPUTE);
   calculateBVOLDerivativesSimple(mpiGrid, sysBoundaries, localCells);
   
   return true;
}

bool finalizeFieldPropagator(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
) {
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
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   cuint subcycles
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
      propagateMagneticFieldSimple(mpiGrid, sysBoundaries, dt, localCells, RK_ORDER1);
      calculateDerivativesSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER1, true);
      if(P::ohmGradPeTerm > 0){
         calculateGradPeTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER1);
         hallTermCommunicateDerivatives = false;
      }
      if(P::ohmHallTerm > 0) {
         calculateHallTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER1, hallTermCommunicateDerivatives);
      }
      calculateUpwindedElectricFieldSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER1);
      #else
      propagateMagneticFieldSimple(mpiGrid, sysBoundaries, dt, localCells, RK_ORDER2_STEP1);
      calculateDerivativesSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP1, true);
      if(P::ohmGradPeTerm > 0) {
         calculateGradPeTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP1);
         hallTermCommunicateDerivatives = false;
      }
      if(P::ohmHallTerm > 0) {
         calculateHallTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP1, hallTermCommunicateDerivatives);
      }
      calculateUpwindedElectricFieldSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP1);
      
      propagateMagneticFieldSimple(mpiGrid, sysBoundaries, dt, localCells, RK_ORDER2_STEP2);
      calculateDerivativesSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP2, true);
      if(P::ohmGradPeTerm > 0) {
         calculateGradPeTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP2);
         hallTermCommunicateDerivatives = false;
      }
      if(P::ohmHallTerm > 0) {
         calculateHallTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP2, hallTermCommunicateDerivatives);
      }
      calculateUpwindedElectricFieldSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP2);
      #endif
   } else {
      const vector<CellID> cells = mpiGrid.get_cells();
      Real subcycleDt = dt/convert<Real>(subcycles);
      Real subcycleT = P::t;
      creal targetT = P::t + dt;
      uint subcycleCount = 0;
      uint maxSubcycleCount = std::numeric_limits<uint>::max();
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      
      while (subcycleCount < maxSubcycleCount ) {         
         // In case of subcycling, we decided to go for a blunt Runge-Kutta subcycling even though e.g. moments are not going along.
         // Result of the Summer of Debugging 2016, the behaviour in wave dispersion was much improved with this.
         propagateMagneticFieldSimple(mpiGrid, sysBoundaries, subcycleDt, localCells, RK_ORDER2_STEP1);
         // We need to calculate derivatives of the moments at every substep, but they only
         // need to be communicated in the first one.
         calculateDerivativesSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP1, (subcycleCount==0));
         if(P::ohmGradPeTerm > 0 && subcycleCount==0) {
            calculateGradPeTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP1);
            hallTermCommunicateDerivatives = false;
         }
         if(P::ohmHallTerm > 0) {
            calculateHallTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP1, hallTermCommunicateDerivatives);
         }
         calculateUpwindedElectricFieldSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP1);
         
         propagateMagneticFieldSimple(mpiGrid, sysBoundaries, subcycleDt, localCells, RK_ORDER2_STEP2);
         // We need to calculate derivatives of the moments at every substep, but they only
         // need to be communicated in the first one.
         calculateDerivativesSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP2, (subcycleCount==0));
         if(P::ohmGradPeTerm > 0 && subcycleCount==0) {
            calculateGradPeTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP2);
            hallTermCommunicateDerivatives = false;
         }
         if(P::ohmHallTerm > 0) {
            calculateHallTermSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP2, hallTermCommunicateDerivatives);
         }
         calculateUpwindedElectricFieldSimple(mpiGrid, sysBoundaries, localCells, RK_ORDER2_STEP2);
         
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

         for (std::vector<uint64_t>::const_iterator cell_id = cells.begin(); cell_id != cells.end(); ++cell_id) {
            SpatialCell* cell = mpiGrid[*cell_id];
            if ( cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
               (cell->sysBoundaryLayer == 1 && cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY )) {
               dtMaxLocal=min(dtMaxLocal, cell->parameters[CellParams::MAXFDT]);
            }
         }
         phiprof::start("MPI_Allreduce");
         MPI_Allreduce(&(dtMaxLocal), &(dtMaxGlobal), 1, MPI_Type<Real>(), MPI_MIN, MPI_COMM_WORLD);
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
   
   calculateVolumeAveragedFields(mpiGrid,fs_cache::getCache().localCellsCache,fs_cache::getCache().local_NOT_DO_NOT_COMPUTE);
   calculateBVOLDerivativesSimple(mpiGrid, sysBoundaries, localCells);
   return true;
}
