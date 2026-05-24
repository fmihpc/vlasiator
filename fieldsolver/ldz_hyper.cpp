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
#include "ldz_hyper.hpp"

#ifdef DEBUG_VLASIATOR
   #define DEBUG_FSOLVER
#endif

using namespace std;

/*
   Functions for calculating the contribution of hyperresistivity to the electric field:
   E_hyper = -eta_H nabla^2 (curl(B)) / mu_0

   Currently, the calculations are done fully on the Yee lattice, without any Balsara interpolation or volume averaging.
 */
Real calculateSecondDerivativeOfCurl(
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   cint i,
   cint j,
   cint k,
   cint curlComp,
   cint derComp
) {
   std::array<Real, fsgrids::dperb::N_DPERB> * centdPerB = dPerBGrid.get(i,j,k);
   std::array<Real, fsgrids::dperb::N_DPERB> * middPerB = NULL;
   std::array<Real, fsgrids::dperb::N_DPERB> * leftdPerB = NULL;
   std::array<Real, fsgrids::dperb::N_DPERB> * rightdPerB = NULL;
   Real SecondDerOfCurl = 0.0;

   switch (derComp) {

      case 0:
         leftdPerB = dPerBGrid.get(i-1,j,k);
         rightdPerB = dPerBGrid.get(i+1,j,k);
         break;

      case 1:
         leftdPerB = dPerBGrid.get(i,j-1,k);
         rightdPerB = dPerBGrid.get(i,j+1,k);
         break;

      case 2:
         leftdPerB = dPerBGrid.get(i,j,k-1);
         rightdPerB = dPerBGrid.get(i,j,k+1);
         break;

      default:
         cerr << __FILE__ << ":" << __LINE__ << "Derivative component should be 0 (xx), 1 (yy), or 2 (zz)." << endl;
         break;
   }

   middPerB = centdPerB;
   if (leftdPerB == NULL || rightdPerB == NULL) {
      return 0.0;
   }

   switch (curlComp) {
      case 0:
         SecondDerOfCurl = (leftdPerB->at(fsgrids::dperb::dPERBzdy) + rightdPerB->at(fsgrids::dperb::dPERBzdy) - 2 * middPerB->at(fsgrids::dperb::dPERBzdy)) / (dPerBGrid.DX * dPerBGrid.DX * dPerBGrid.DX) - 
         (leftdPerB->at(fsgrids::dperb::dPERBydz) + rightdPerB->at(fsgrids::dperb::dPERBydz) - 2 * middPerB->at(fsgrids::dperb::dPERBydz)) / (dPerBGrid.DX * dPerBGrid.DX * dPerBGrid.DX);
         break;

      case 1:
         SecondDerOfCurl = (leftdPerB->at(fsgrids::dperb::dPERBxdz) + rightdPerB->at(fsgrids::dperb::dPERBxdz) - 2 * middPerB->at(fsgrids::dperb::dPERBxdz)) / (dPerBGrid.DX * dPerBGrid.DX * dPerBGrid.DX) - 
         (leftdPerB->at(fsgrids::dperb::dPERBzdx) + rightdPerB->at(fsgrids::dperb::dPERBzdx) - 2 * middPerB->at(fsgrids::dperb::dPERBzdx)) / (dPerBGrid.DX * dPerBGrid.DX * dPerBGrid.DX);
         break;

      case 2:
         SecondDerOfCurl = (leftdPerB->at(fsgrids::dperb::dPERBydx) + rightdPerB->at(fsgrids::dperb::dPERBydx) - 2 * middPerB->at(fsgrids::dperb::dPERBydx)) / (dPerBGrid.DX * dPerBGrid.DX * dPerBGrid.DX) - 
         (leftdPerB->at(fsgrids::dperb::dPERBxdy) + rightdPerB->at(fsgrids::dperb::dPERBxdy) - 2 * middPerB->at(fsgrids::dperb::dPERBxdy)) / (dPerBGrid.DX * dPerBGrid.DX * dPerBGrid.DX);
         break;
      
      default:
         cerr << __FILE__ << ":" << __LINE__ << "Curl component should be 0 (x), 1 (y), or 2 (z)." << endl;
         break;
   }

   return SecondDerOfCurl;

}

void calculateEdgeHyperTermXComponents(
   FsGrid< std::array<Real, fsgrids::ehyper::N_EHYPER>, FS_STENCIL_WIDTH> & EHyperGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   cint i,
   cint j,
   cint k
) {
   Real limitedRhoq = 0.0;
   Real rhoq = 0.0;
   Real EHyperX = 0.0;
   if (Parameters::ohmHyperTerm == 0) {
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a hyperresistivity term function if Parameters::ohmHyperTerm == 0." << endl;
   }
   std::array<Real, fsgrids::bfield::N_BFIELD> * centPerB = perBGrid.get(i,j,k);
   std::array<Real, fsgrids::bgbfield::N_BGB> * centBgB = BgBGrid.get(i,j,k);
   std::array<Real, fsgrids::moments::N_MOMENTS> * centMom = momentsGrid.get(i,j,k);

   Real Bmag = sqrt((centBgB->at(fsgrids::bgbfield::BGBX)+centPerB->at(fsgrids::bfield::PERBX))*
            (centBgB->at(fsgrids::bgbfield::BGBX)+centPerB->at(fsgrids::bfield::PERBX)) +
            (centBgB->at(fsgrids::bgbfield::BGBY)+centPerB->at(fsgrids::bfield::PERBY))*
            (centBgB->at(fsgrids::bgbfield::BGBY)+centPerB->at(fsgrids::bfield::PERBY)) +
            (centBgB->at(fsgrids::bgbfield::BGBZ)+centPerB->at(fsgrids::bfield::PERBZ))*
            (centBgB->at(fsgrids::bgbfield::BGBZ)+centPerB->at(fsgrids::bfield::PERBZ))
           );
   rhoq = centMom->at(fsgrids::moments::RHOQ);
   limitedRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;

   // -eta_H/mu_0 calculated with center values
   Real hyperres_coeff = -1.0 / 4.0 / M_PI / M_PI *
           dPerBGrid.DX * dPerBGrid.DX *
           Bmag / limitedRhoq / physicalconstants::MU_0;

   // xx-derivative
   EHyperX += P::ohmHyperFactor * hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,0,0);

   // yy-derivative
   EHyperX += P::ohmHyperFactor * hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,0,1);

   // zz-derivative
   EHyperX += P::ohmHyperFactor * hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,0,2);

   EHyperGrid.get(i,j,k)->at(fsgrids::ehyper::EXHYPER) = EHyperX;
}


void calculateEdgeHyperTermYComponents(
   FsGrid< std::array<Real, fsgrids::ehyper::N_EHYPER>, FS_STENCIL_WIDTH> & EHyperGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   cint i,
   cint j,
   cint k
) {
   Real limitedRhoq = 0.0;
   Real rhoq = 0.0;
   Real EHyperY = 0.0;
   if (Parameters::ohmHyperTerm == 0) {
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a hyperresistivity term function if Parameters::ohmHyperTerm == 0." << endl;
   }
   std::array<Real, fsgrids::bfield::N_BFIELD> * centPerB = perBGrid.get(i,j,k);
   std::array<Real, fsgrids::bgbfield::N_BGB> * centBgB = BgBGrid.get(i,j,k);
   std::array<Real, fsgrids::moments::N_MOMENTS> * centMom = momentsGrid.get(i,j,k);

   Real Bmag = sqrt((centBgB->at(fsgrids::bgbfield::BGBX)+centPerB->at(fsgrids::bfield::PERBX))*
            (centBgB->at(fsgrids::bgbfield::BGBX)+centPerB->at(fsgrids::bfield::PERBX)) +
            (centBgB->at(fsgrids::bgbfield::BGBY)+centPerB->at(fsgrids::bfield::PERBY))*
            (centBgB->at(fsgrids::bgbfield::BGBY)+centPerB->at(fsgrids::bfield::PERBY)) +
            (centBgB->at(fsgrids::bgbfield::BGBZ)+centPerB->at(fsgrids::bfield::PERBZ))*
            (centBgB->at(fsgrids::bgbfield::BGBZ)+centPerB->at(fsgrids::bfield::PERBZ))
           );
   rhoq = centMom->at(fsgrids::moments::RHOQ);
   limitedRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;

   // -eta_H/mu_0 calculated with center values
   Real hyperres_coeff = -1.0 / 4.0 / M_PI / M_PI *
           dPerBGrid.DY * dPerBGrid.DY *
           Bmag / limitedRhoq / physicalconstants::MU_0;

   // xx-derivative
   EHyperY += P::ohmHyperFactor * hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,1,0);

   // yy-derivative
   EHyperY += P::ohmHyperFactor * hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,1,1);

   // zz-derivative
   EHyperY += P::ohmHyperFactor * hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,1,2);

   EHyperGrid.get(i,j,k)->at(fsgrids::ehyper::EYHYPER) = EHyperY;
}

void calculateEdgeHyperTermZComponents(
   FsGrid< std::array<Real, fsgrids::ehyper::N_EHYPER>, FS_STENCIL_WIDTH> & EHyperGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   cint i,
   cint j,
   cint k
) {
   Real limitedRhoq = 0.0;
   Real rhoq = 0.0;
   Real EHyperZ = 0.0;
   if (Parameters::ohmHyperTerm == 0) {
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a hyperresistivity term function if Parameters::ohmHyperTerm == 0." << endl;
   }
   std::array<Real, fsgrids::bfield::N_BFIELD> * centPerB = perBGrid.get(i,j,k);
   std::array<Real, fsgrids::bgbfield::N_BGB> * centBgB = BgBGrid.get(i,j,k);
   std::array<Real, fsgrids::moments::N_MOMENTS> * centMom = momentsGrid.get(i,j,k);

   Real Bmag = sqrt((centBgB->at(fsgrids::bgbfield::BGBX)+centPerB->at(fsgrids::bfield::PERBX))*
            (centBgB->at(fsgrids::bgbfield::BGBX)+centPerB->at(fsgrids::bfield::PERBX)) +
            (centBgB->at(fsgrids::bgbfield::BGBY)+centPerB->at(fsgrids::bfield::PERBY))*
            (centBgB->at(fsgrids::bgbfield::BGBY)+centPerB->at(fsgrids::bfield::PERBY)) +
            (centBgB->at(fsgrids::bgbfield::BGBZ)+centPerB->at(fsgrids::bfield::PERBZ))*
            (centBgB->at(fsgrids::bgbfield::BGBZ)+centPerB->at(fsgrids::bfield::PERBZ))
           );
   rhoq = centMom->at(fsgrids::moments::RHOQ);
   limitedRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;

   // -eta_H/mu_0 calculated with center values
   Real hyperres_coeff = -1.0 / 4.0 / M_PI / M_PI *
           dPerBGrid.DZ * dPerBGrid.DZ *
           Bmag / limitedRhoq / physicalconstants::MU_0;

   // xx-derivative
   EHyperZ += P::ohmHyperFactor * hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,2,0);

   // yy-derivative
   EHyperZ += P::ohmHyperFactor * hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,2,1);

   // zz-derivative
   EHyperZ += P::ohmHyperFactor * hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,2,2);

   EHyperGrid.get(i,j,k)->at(fsgrids::ehyper::EZHYPER) = EHyperZ;
}

/** Calculate the hyperresistivity term on all given cells.
 * @param sysBoundaries System boundary condition functions.
 */
void calculateHyperTerm(
   FsGrid< std::array<Real, fsgrids::ehyper::N_EHYPER>, FS_STENCIL_WIDTH> & EHyperGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
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
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondHyperElectricField(EHyperGrid,i,j,k,0);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondHyperElectricField(EHyperGrid,i,j,k,1);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondHyperElectricField(EHyperGrid,i,j,k,2);
   } else {
      calculateEdgeHyperTermXComponents(EHyperGrid,momentsGrid,dPerBGrid,perBGrid,BgBGrid,technicalGrid,i,j,k);
      calculateEdgeHyperTermYComponents(EHyperGrid,momentsGrid,dPerBGrid,perBGrid,BgBGrid,technicalGrid,i,j,k);
      calculateEdgeHyperTermZComponents(EHyperGrid,momentsGrid,dPerBGrid,perBGrid,BgBGrid,technicalGrid,i,j,k);
   }
}

void calculateHyperTermSimple(
   FsGrid< std::array<Real, fsgrids::ehyper::N_EHYPER>, FS_STENCIL_WIDTH> & EHyperGrid,
   FsGrid< std::array<Real, fsgrids::ehyper::N_EHYPER>, FS_STENCIL_WIDTH> & EHyperDt2Grid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const FsGridTools::FsIndex_t* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::Timer hyperTimer {"Calculate Hyper term"};
   int computeTimerId {phiprof::initializeTimer("Ehyper compute cells")};

   phiprof::Timer mpiTimer {"Hyper field update ghosts MPI", {"MPI"}};
   if (Parameters::ohmHallTerm == 0) {
      dPerBGrid.updateGhostCells();
   }
   mpiTimer.stop();

   // Calculate Hyper term
   #pragma omp parallel
   {
      phiprof::Timer computeTimer {computeTimerId};
      #pragma omp for collapse(2)
      for (FsGridTools::FsIndex_t k=0; k<gridDims[2]; k++) {
         for (FsGridTools::FsIndex_t j=0; j<gridDims[1]; j++) {
            for (FsGridTools::FsIndex_t i=0; i<gridDims[0]; i++) {
               if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
                  calculateHyperTerm(EHyperGrid, momentsGrid, dPerBGrid, perBGrid, BgBGrid, technicalGrid, i, j, k, sysBoundaries);
               } else {
                  calculateHyperTerm(EHyperDt2Grid, momentsDt2Grid, dPerBGrid, perBDt2Grid, BgBGrid, technicalGrid, i, j, k, sysBoundaries);
               }
            }
         }
      }
      computeTimer.stop(N_cells,"Spatial Cells");
   }

   hyperTimer.stop(N_cells,"Spatial Cells");
}
