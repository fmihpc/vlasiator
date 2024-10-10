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

#include "fs_common.h"
#include "ldz_gradpe.hpp"

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

using namespace std;

void calculateEdgeGradPeTermXComponents(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   Real hallRhoq = 0.0;
   Real rhoq = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;

      case 1:
         rhoq = momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ);
         hallRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
         //EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE) = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::drhoqdx) / (hallRhoq*EGradPeGrid.DX);
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE) = - dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::dPedx) / (hallRhoq*EGradPeGrid.DX);
	 break;

      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermYComponents(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   Real hallRhoq = 0.0;
   Real rhoq = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;

      case 1:
         rhoq = momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ);
         hallRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
         //EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EYGRADPE) = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::drhoqdy) / (hallRhoq*EGradPeGrid.DY);
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EYGRADPE) = - dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::dPedy) / (hallRhoq*EGradPeGrid.DY);
         break;

      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermZComponents(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   Real hallRhoq = 0.0;
   Real rhoq = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;

      case 1:
         rhoq = momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ);
         hallRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
         //EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EZGRADPE) = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::drhoqdz) / (hallRhoq*EGradPeGrid.DZ);
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EZGRADPE) = - dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::dPedz) / (hallRhoq*EGradPeGrid.DZ);
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
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
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

   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE || cellSysBoundaryFlag == sysboundarytype::OUTER_BOUNDARY_PADDING) return;

   cuint cellSysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;

   if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && (cellSysBoundaryLayer != 1)) {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,i,j,k,0);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,i,j,k,1);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,i,j,k,2);
   } else {
      calculateEdgeGradPeTermXComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
      calculateEdgeGradPeTermYComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
      calculateEdgeGradPeTermZComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
   }
}

void calculateGradPeTermSimple(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeDt2Grid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsDt2Grid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const FsGridTools::FsIndex_t* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::Timer gradPeTimer {"Calculate GradPe term"};
   int computeTimerId {phiprof::initializeTimer("EgradPe compute cells")};

   phiprof::Timer mpiTimer {"EgradPe field update ghosts MPI", {"MPI"}};
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      dMomentsGrid.updateGhostCells();
   } else {
      dMomentsDt2Grid.updateGhostCells();
   }
   mpiTimer.stop();

   // Calculate GradPe term
   #pragma omp parallel
   {
      phiprof::Timer computeTimer {computeTimerId};
      #pragma omp for collapse(2)
      for (FsGridTools::FsIndex_t k=0; k<gridDims[2]; k++) {
         for (FsGridTools::FsIndex_t j=0; j<gridDims[1]; j++) {
            for (FsGridTools::FsIndex_t i=0; i<gridDims[0]; i++) {
               if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
                  calculateGradPeTerm(EGradPeGrid, momentsGrid, dMomentsGrid, technicalGrid, i, j, k, sysBoundaries);
               } else {
                  calculateGradPeTerm(EGradPeDt2Grid, momentsDt2Grid, dMomentsDt2Grid, technicalGrid, i, j, k, sysBoundaries);
               }
            }
         }
      }
      computeTimer.stop(N_cells,"Spatial Cells");
   }

   gradPeTimer.stop(N_cells,"Spatial Cells");
}
