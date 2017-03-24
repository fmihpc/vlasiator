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
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   #warning Particles (charge) assumed to be protons here
   Real hallRho = 0.0;
   Real rho = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         rho = momentsGrid.get(i,j,k)->at(fsgrids::moments::RHO);
         hallRho = (rho <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : rho ;
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE) = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::drhodx) / (hallRho*physicalconstants::CHARGE*EGradPeGrid.DX);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermYComponents(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   #warning Particles (charge) assumed to be protons here
   Real hallRho = 0.0;
   Real rho = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         rho = momentsGrid.get(i,j,k)->at(fsgrids::moments::RHO);
         hallRho = (rho <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : rho ;
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EYGRADPE) = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::drhody) / (hallRho*physicalconstants::CHARGE*EGradPeGrid.DY);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermZComponents(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
  #warning Particles (charge) assumed to be protons here
   Real hallRho = 0.0;
   Real rho = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         rho = momentsGrid.get(i,j,k)->at(fsgrids::moments::RHO);
         hallRho = (rho <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : rho ;
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EZGRADPE) = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::drhodz) / (hallRho*physicalconstants::CHARGE*EGradPeGrid.DZ);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

/** Calculate the electron pressure gradient term on all given cells.
 * @param sysBoundaries System boundary condition functions.
 */
void calculateGradPeTerm(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   cint i,
   cint j,
   cint k,
   SysBoundary& sysBoundaries
) {
   #ifdef DEBUG_FSOLVER
   if (technicalGrid.get(i,j,k) == NULL) {
      cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
   #endif
   
   cuint cellSysBoundaryFlag = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   
   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) return;
   
   cuint cellSysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;
   
   if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && (cellSysBoundaryLayer != 1)) {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,i,j,k,0);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,i,j,k,1);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,i,j,k,2);
   } else {
      calculateEdgeGradPeTermXComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
      calculateEdgeGradPeTermXComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
      calculateEdgeGradPeTermYComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
   }
}

void calculateGradPeTermSimple(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {
   int timer;
   const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::start("Calculate GradPe term");

   timer=phiprof::initializeTimer("MPI","MPI");
   phiprof::start(timer);
   dMomentsGrid.updateGhostCells();
   phiprof::stop(timer);

   // Calculate GradPe term
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);
   #pragma omp parallel for collapse(3)
   for (int k=0; k<gridDims[2]; k++) {
      for (int j=0; j<gridDims[1]; j++) {
         for (int i=0; i<gridDims[0]; i++) {
            if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
               calculateGradPeTerm(EGradPeGrid, momentsGrid, dMomentsGrid, technicalGrid, i, j, k, sysBoundaries);
            } else {
               calculateGradPeTerm(EGradPeGrid, momentsDt2Grid, dMomentsGrid, technicalGrid, i, j, k, sysBoundaries);
            }
         }
      }
   }
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   phiprof::stop("Calculate GradPe term",N_cells,"Spatial Cells");
}
