/*
This file is part of Vlasiator.

Copyright 2015 Finnish Meteorological Institute

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
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   #warning Particles (charge) assumed to be protons here
   Real hallRho;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         creal rho = momentsGrid.get(i,j,k)[fsgrids::moments::RHO];
         hallRho = (rho <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : rho ;
         EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EXGRADPE] = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid,get(i,j,k)[fsgrids::dmoments::drhodx] / (hallRho*physicalconstants::CHARGE*EGradPeGrid.DX);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermYComponents(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   #warning Particles (charge) assumed to be protons here
   Real hallRho;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         creal rho = momentsGrid.get(i,j,k)[fsgrids::moments::RHO];
         hallRho = (rho <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : rho ;
         EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EYGRADPE] = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid,get(i,j,k)[fsgrids::dmoments::drhody] / (hallRho*physicalconstants::CHARGE*EGradPeGrid.DY);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermZComponents(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
  #warning Particles (charge) assumed to be protons here
   Real hallRho;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         creal rho = momentsGrid.get(i,j,k)[fsgrids::moments::RHO];
         hallRho = (rho <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : rho ;
         EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EZGRADPE] = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid,get(i,j,k)[fsgrids::dmoments::drhodz] / (hallRho*physicalconstants::CHARGE*EGradPeGrid.DZ);
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
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   const FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   const FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   const FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   const int i,
   const int j,
   const int k,
   SysBoundary& sysBoundaries
) {
   #ifdef DEBUG_FSOLVER
   if (technicalGrid.get(i,j,k) == NULL) {
      cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
   #endif
   
   cuint cellSysBoundaryFlag = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   
   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
   
   cuint cellSysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;
   
   if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && (cellSysBoundaryLayer != 1)) {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,cache[localID],RKCase,0);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,cache[localID],RKCase,1);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,cache[localID],RKCase,2);
   } else {
      calculateEdgeGradPeTermXComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
      calculateEdgeGradPeTermXComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
      calculateEdgeGradPeTermYComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
   }
}

void calculateGradPeTermSimple(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase
) {
   namespace fs = fieldsolver;
   int timer;
   const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::start("Calculate GradPe term");

   timer=phiprof::initializeTimer("Start communication of derivatives","MPI");
   phiprof::start(timer);
   dMomentsGrid.updateGhostCells();
   phiprof::stop(timer);

   // Calculate GradPe term
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);
   #pragma omp parallel for collapse(3)
   for (uint k=0; k<gridDims[2]; k++) {
      for (uint j=0; j<gridDims[1]; j++) {
         for (uint i=0; i<gridDims[0]; i++) {
            if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
               calculateGradPeTerm(EGradPeGrid, momentsGrid, dMomentsGrid, technicalGrid, i, j, k, sysBoundaries);
            }
            if (RKCase == RK_ORDER2_STEP1) {
               calculateGradPeTerm(EGradPeGrid, momentsDt2Grid, dMomentsGrid, technicalGrid, i, j, k, sysBoundaries);
            }
         }
      }
   }
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   phiprof::stop("Calculate GradPe term",N_cells,"Spatial Cells");
}
