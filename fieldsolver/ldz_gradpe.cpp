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

#ifdef DEBUG_VLASIATOR
   #define DEBUG_FSOLVER
#endif

using namespace std;

/**
   Vlasiator considers electrons an inertia-free charge-neutralizing fluid. To model the effects of electron pressure on
   plasma dynamics (such as the cross-shock potential, an electric field induced by electron pressure changes), we consider the
   electron pressure gradient term in the generalized Ohm's law:

   - nabla dot P_e /(n_e e)

   To model this, we can assume some equation of state for electrons. We now assume that the electron fluid has an isotropic
   pressure tensor (so pressure is a scalar), and that it is governed by a polytropic process, where P * V^k is constant. 
   Here P is the pressure, V is the volume of an unit of fluid, and k is a polytropic index describing the process.

   k = 0     : Isobaric process (pressure is constant)
   k = 1     : Isothermal process (temperature is constant)
   k = gamma : Adiabatic process (no energy transfer)

   here gamma is the ratio of specific heats, that is the heat capacity at constant pressure (C_P) divided by the 
   heat capacity at constant volume (C_V). For an ideal monoatomic gas, gamma = 5/3.

   Now because P * V^k is constant, we can evaluate that constant for the whole simulation based on some anchor point,
   usually considered the upstream inflow conditions. These anchor point values are provided by the user as
   config parameters. The anchor point electron density should match the density of charge from ions, but It is left up
   to the user to decide if anchor point electron temperature should be equal to the proton temperature or if there should be
   a correction constant (e.g. 1/4).

   In derivatives.cpp we first calculate the electron pressure at this anchor point:
   Real Pe_anchor = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   And then calculate the constant used further on:
   Real Pe_const = Pe_anchor * pow(Parameters::electronDensity, -Parameters::electronPTindex);

   Remembering that P * V^k = P * (n_e)^-k is constant, we can calculate that constant at the anchor point and use that for solving the
   electron pressure at any given position as 
   P(r) = Pe_const * (n_e(r))^k
   where n_e is the electron number density (inverse of the volume of a given parcel of electron fluid), which is solved by
   dividing the ion charge density fsgrids::moments::RHOQ with the elementary charge.

   Thus, the gradient of electron pressure is calculated in derivatives.cpp using a standard slope limiter, for each Cartesian direction. 

   Note:
   The value stored in e.g. dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::drhoqdx) is a difference, not a derivative. Thus, in the next step,
   evaluating the electron pressure gradient term in the general Ohm's law as
   E_gradPe(x,y,z) = - dPe/d(x,y,z) / (n_e * e * Delta(x,y,z))
   requires adding EGradPeGrid.DX/DY/DZ in the denominator.

   Similar to the Hall term, we use the Parameters::hallMinimumRhoq value as a lower bound for electron number density in order to
   prevent runaway electric field terms at depletion zones.

 */


void calculateEdgeGradPeTermXComponents(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   Real limitedRhoq = 0.0;
   Real rhoq = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;

      case 1:
         rhoq = momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ);
         limitedRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE) = - dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::dPedx) / (limitedRhoq*EGradPeGrid.DX);
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
   Real limitedRhoq = 0.0;
   Real rhoq = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;

      case 1:
         rhoq = momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ);
         limitedRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EYGRADPE) = - dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::dPedy) / (limitedRhoq*EGradPeGrid.DY);
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
   Real limitedRhoq = 0.0;
   Real rhoq = 0.0;
   switch (Parameters::ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0." << endl;
         break;

      case 1:
         rhoq = momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ);
         limitedRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EZGRADPE) = - dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::dPedz) / (limitedRhoq*EGradPeGrid.DZ);
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
