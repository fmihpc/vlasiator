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

#include "fs_common.h"
#include "fs_cache.h"
#include "ldz_gradpe.hpp"

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

using namespace std;

// // X
// template<typename REAL> inline
// REAL gradPeX_000_100(
//    const REAL* const pC,
//    creal BGBY,
//    creal BGBZ,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// template<typename REAL> inline
// REAL gradPeX_010_110(
//    const REAL* const pC,
//    creal BGBY,
//    creal BGBZ,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// template<typename REAL> inline
// REAL gradPeX_001_101(
//    const REAL* const pC,
//    creal BGBY,
//    creal BGBZ,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// template<typename REAL> inline
// REAL gradPeX_011_111(
//    const REAL* const pC,
//    creal BGBY,
//    creal BGBZ,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// // Y
// template<typename REAL> inline
// REAL gradPeY_000_010(
//    const REAL* const pC,
//    creal BGBX,
//    creal BGBZ,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// template<typename REAL> inline
// REAL gradPeY_100_110(
//    const REAL* const pC,
//    creal BGBX,
//    creal BGBZ,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// template<typename REAL> inline
// REAL gradPeY_001_011(
//    const REAL* const pC,
//    creal BGBX,
//    creal BGBZ,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// template<typename REAL> inline
// REAL gradPeY_101_111(
//    const REAL* const pC,
//    creal BGBX,
//    creal BGBZ,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// // Z
// template<typename REAL> inline
// REAL gradPeZ_000_001(
//    const REAL* const pC,
//    creal BGBX,
//    creal BGBY,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// template<typename REAL> inline
// REAL gradPeZ_100_101(
//    const REAL* const pC,
//    creal BGBX,
//    creal BGBY,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// template<typename REAL> inline
// REAL gradPeZ_010_011(
//    const REAL* const pC,
//    creal BGBX,
//    creal BGBY,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }
// 
// template<typename REAL> inline
// REAL gradPeZ_110_111(
//    const REAL* const pC,
//    creal BGBX,
//    creal BGBY,
//    creal dx,
//    creal dy,
//    creal dz
// ) {
//    return;
// }

void calculateEdgeGradPeTermXComponents(
   Real* cp,
   Real* derivs,
   cint& RKCase
) {
   Real hallRhoq = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            hallRhoq =  (cp[CellParams::RHOQ] <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : cp[CellParams::RHOQ] ;
         }
         if (RKCase == RK_ORDER2_STEP1) {
            hallRhoq =  (cp[CellParams::RHOQ_DT2] <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : cp[CellParams::RHOQ_DT2] ;
         }
         cp[CellParams::EXGRADPE] = -physicalconstants::K_B*Parameters::electronTemperature*derivs[fieldsolver::drhoqdx] / (hallRhoq*cp[CellParams::DX]);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermYComponents(
   Real* cp,
   Real* derivs,
   cint& RKCase
) {
   Real hallRhoq = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            hallRhoq =  (cp[CellParams::RHOQ] <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : cp[CellParams::RHOQ] ;
         }
         if (RKCase == RK_ORDER2_STEP1) {
            hallRhoq =  (cp[CellParams::RHOQ_DT2] <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : cp[CellParams::RHOQ_DT2] ;
         }
         cp[CellParams::EYGRADPE] = -physicalconstants::K_B*Parameters::electronTemperature*derivs[fieldsolver::drhoqdy] / (hallRhoq*physicalconstants::CHARGE*cp[CellParams::DY]);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermZComponents(
   Real* cp,
   Real* derivs,
   cint& RKCase
) {
   Real hallRhoq = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            hallRhoq =  (cp[CellParams::RHOQ] <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : cp[CellParams::RHOQ] ;
         }
         if (RKCase == RK_ORDER2_STEP1) {
            hallRhoq =  (cp[CellParams::RHOQ_DT2] <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : cp[CellParams::RHOQ_DT2] ;
         }
         cp[CellParams::EZGRADPE] = -physicalconstants::K_B*Parameters::electronTemperature*derivs[fieldsolver::drhoqdz] / (hallRhoq*physicalconstants::CHARGE*cp[CellParams::DZ]);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

/** Calculate the electron pressure gradient term on all given cells.
 * @param sysBoundaries System boundary condition functions.
 * @param cache Cache for local cells.
 * @param cells Local IDs of calculated cells, one of the vectors in fs_cache::CacheContainer.
 * @param RKCase
 */
void calculateGradPeTerm(
   SysBoundary& sysBoundaries,
   std::vector<fs_cache::CellCache>& cache,
   const std::vector<uint16_t>& cells,
   cint& RKCase
) {
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) { // DO_NOT_COMPUTE cells already removed
      const uint16_t localID = cells[c];

      #ifdef DEBUG_FSOLVER
      if (localID >= cache.size()) {
         cerr << "local index out of bounds in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }
      if (cache[localID].cells[fs_cache::calculateNbrID(1,1,1)] == NULL) {
         cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }
      if (cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->parameters == NULL) {
         cerr << "NULL cell parameters in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }
      if (cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->derivatives == NULL) {
         cerr << "NULL derivatives in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }
      #endif

      cuint fieldSolverSysBoundaryFlag = cache[localID].existingCellsFlags;
      cuint cellSysBoundaryFlag        = cache[localID].sysBoundaryFlag;
      cuint cellSysBoundaryLayer       = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->sysBoundaryLayer;

      Real perturbedCoefficients[Rec::N_REC_COEFFICIENTS];

      if ((fieldSolverSysBoundaryFlag & CALCULATE_EX) == CALCULATE_EX) {
         if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) &&
             (cellSysBoundaryLayer != 1)) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(cache[localID],RKCase,0);
         } else {
            Real* cp     = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
            Real* derivs = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->derivatives;
            calculateEdgeGradPeTermXComponents(cp,derivs,RKCase);
         }
      }
      if ((fieldSolverSysBoundaryFlag & CALCULATE_EY) == CALCULATE_EY) {
         if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) &&
             (cellSysBoundaryLayer != 1)) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(cache[localID],RKCase,1);
         } else {
            Real* cp     = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
            Real* derivs = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->derivatives;
            calculateEdgeGradPeTermYComponents(cp,derivs,RKCase);
         }
      }
      if ((fieldSolverSysBoundaryFlag & CALCULATE_EZ) == CALCULATE_EZ) {
         if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) &&
             (cellSysBoundaryLayer != 1)) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(cache[localID],RKCase,2);
         } else {
            Real* cp     = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
            Real* derivs = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->derivatives;
            calculateEdgeGradPeTermZComponents(cp,derivs,RKCase);
         }
      }
   }
}

void calculateGradPeTermSimple(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase
) {
   namespace fs = fieldsolver;
   int timer;
   phiprof::start("Calculate GradPe term");
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_DERIVATIVES);

   fs_cache::CacheContainer& cacheContainer = fs_cache::getCache();

   timer=phiprof::initializeTimer("Start communication of derivatives","MPI");
   phiprof::start(timer);
   mpiGrid.start_remote_neighbor_copy_updates(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);

   // Calculate GradPe term on inner cells
   timer=phiprof::initializeTimer("Compute inner cells");
   phiprof::start(timer);
   calculateGradPeTerm(sysBoundaries,cacheContainer.localCellsCache,cacheContainer.cellsWithLocalNeighbours,RKCase);
   phiprof::stop(timer,cacheContainer.cellsWithLocalNeighbours.size(),"Spatial Cells");

   timer=phiprof::initializeTimer("Wait for receives","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_receives(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);
   
   // Calculate GradPe term on boundary cells:
   timer=phiprof::initializeTimer("Compute boundary cells");
   phiprof::start(timer);
   calculateGradPeTerm(sysBoundaries,cacheContainer.localCellsCache,cacheContainer.cellsWithRemoteNeighbours,RKCase);
   phiprof::stop(timer,cacheContainer.cellsWithRemoteNeighbours.size(),"Spatial Cells");

   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_sends();
   phiprof::stop(timer);

   const size_t N_cells = cacheContainer.cellsWithRemoteNeighbours.size()
     + cacheContainer.cellsWithLocalNeighbours.size();

   phiprof::stop("Calculate GradPe term",N_cells,"Spatial Cells");
}
