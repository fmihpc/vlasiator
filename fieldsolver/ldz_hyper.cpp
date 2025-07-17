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

# define M_PI		3.14159265358979323846	/* pi */

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

   Using this equation, for any position in the simulation, the electron temperature can be evaluated based on the local
   electron number density (assumed to be equal to the ion charge density fsgrids::moments::RHOQ divided by the elementary
   charge) and some normalization or anchoring values.

   The polytropic index can be, e.g.:
   k = 0     : Isobaric process (pressure is constant)
   k = 1     : Isothermal process (temperature is constant)
   k = gamma : Adiabatic process (no energy transfer)

   Here gamma is the ratio of specific heats, that is the heat capacity at constant pressure (C_P) divided by the
   heat capacity at constant volume (C_V). For an ideal monoatomic gas, gamma = 5/3.

   Now because P * V^k is constant, we can evaluate that constant for the whole simulation based on the values at
   some hypothetical anchor point, with values usually declared based on the upstream inflow conditions. These values are
   only used for normalizing the relationship between electron density and electron temperature for the given polytropic
   relation. These anchor point values are provided by the user as config parameters - the anchor point is not an actual
   physical point in the simulation, but just imagined for normalization purposes..

   Example in a cfg-script:
   [fieldsolver]
   ohmGradPeTerm = 1 # active
   electronPTindex = 1.666667 # adiabatic, 3/2
   electronDensity = 1.0e6 # anchor point (here: inflow solar wind) electron density value
   electronTemperature = 0.125e6 # electron temperature associated with the anchor point density value. Here
   # the value is set to 1/4 of the solar wind proton temperature.

   The user may, for example, use inflow solar wind electron values as the anchor point values. A good rule of thumb is
   that electron number density times elementary charge should match the ion charge density in the solar wind. A zeroth-order
   approximation is to set electron temperature to equal that of solar wind ions, but alternatively an observation-based factor
   (e.g. 0.25) may be applied to decrease electron temperatures at the anchor point.

   In derivatives.cpp we first calculate the electron pressure at this anchor point:
   Real Pe_anchor = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   And then calculate the constant used further on:
   Real Pe_const = Pe_anchor * pow(Parameters::electronDensity, -Parameters::electronPTindex);

   Remembering that P * V^k = P * (n_e)^-k is constant, we can calculate that constant from the anchor point values,
   and use that for solving the electron pressure at any given position as
   P(r) = Pe_const * (n_e(r))^k

   Thus, the gradient of electron pressure is calculated in derivatives.cpp using a standard slope limiter, for each Cartesian direction.

   Note:
   The value stored in e.g. dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::drhoqdx) is a difference, not a derivative. Thus, in the next step,
   evaluating the electron pressure gradient term in the general Ohm's law as
   E_gradPe(x,y,z) = - dPe/d(x,y,z) / (n_e * e * Delta(x,y,z))
   requires adding EGradPeGrid.DX/DY/DZ in the denominator.

   Similar to the Hall term, we use the Parameters::hallMinimumRhoq value as a lower bound for electron number density in order to
   prevent runaway electric field terms at depletion zones.

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
   // std::array<Real, fsgrids::dperb::N_DPERB> * leftleftdPerB = NULL;
   std::array<Real, fsgrids::dperb::N_DPERB> * rightdPerB = NULL;
   // std::array<Real, fsgrids::dperb::N_DPERB> * rightrightdPerB = NULL;
   Real SecondDerOfCurl = 0.0;

   switch (derComp) {

      case 0:
         leftdPerB = dPerBGrid.get(i-1,j,k);
         // leftleftdPerB = dPerBGrid.get(i-2,j,k);
         rightdPerB = dPerBGrid.get(i+1,j,k);
         // rightrightdPerB = dPerBGrid.get(i+2,j,k);
         break;

      case 1:
         leftdPerB = dPerBGrid.get(i,j-1,k);
         // leftleftdPerB = dPerBGrid.get(i,j-2,k);
         rightdPerB = dPerBGrid.get(i,j+1,k);
         // rightrightdPerB = dPerBGrid.get(i,j+2,k);
         break;

      case 2:
         leftdPerB = dPerBGrid.get(i,j,k-1);
         // leftleftdPerB = dPerBGrid.get(i,j,k-2);
         rightdPerB = dPerBGrid.get(i,j,k+1);
         // rightrightdPerB = dPerBGrid.get(i,j,k+2);
         break;

      default:
         cerr << __FILE__ << ":" << __LINE__ << "Derivative component should be 0 (xx), 1 (yy), or 2 (zz)." << endl;
         break;
   }

   middPerB = centdPerB;
   if (leftdPerB == NULL || rightdPerB == NULL) {
      return 0.0;
   }

   // if (!(leftdPerB == NULL || rightdPerB == NULL)) {
   //    middPerB = centdPerB;
   // } else if (leftdPerB == NULL && !(rightdPerB == NULL || rightrightdPerB == NULL) ) {
   //    leftdPerB = centdPerB;
   //    middPerB = rightdPerB;
   //    rightdPerB = rightrightdPerB;
   // } else if (rightdPerB == NULL && !(leftdPerB == NULL || leftleftdPerB == NULL) ) {
   //    leftdPerB = leftleftdPerB;
   //    middPerB = leftdPerB;
   //    rightdPerB = centdPerB;
   // } else {
   //    leftdPerB =
   //    middPerB =
   //    rightdPerB = centdPerB;
   // }

   switch (curlComp) {
      case 0:
         SecondDerOfCurl = (leftdPerB->at(fsgrids::dperb::dPERBzdy) + rightdPerB->at(fsgrids::dperb::dPERBzdy) - 2 * middPerB->at(fsgrids::dperb::dPERBzdy)) / (dPerBGrid.DX * dPerBGrid.DX) - 
         (leftdPerB->at(fsgrids::dperb::dPERBydz) + rightdPerB->at(fsgrids::dperb::dPERBydz) - 2 * middPerB->at(fsgrids::dperb::dPERBydz)) / (dPerBGrid.DX * dPerBGrid.DX);
         break;

      case 1:
         SecondDerOfCurl = (leftdPerB->at(fsgrids::dperb::dPERBxdz) + rightdPerB->at(fsgrids::dperb::dPERBxdz) - 2 * middPerB->at(fsgrids::dperb::dPERBxdz)) / (dPerBGrid.DX * dPerBGrid.DX) - 
         (leftdPerB->at(fsgrids::dperb::dPERBzdx) + rightdPerB->at(fsgrids::dperb::dPERBzdx) - 2 * middPerB->at(fsgrids::dperb::dPERBzdx)) / (dPerBGrid.DX * dPerBGrid.DX);
         break;

      case 2:
         SecondDerOfCurl = (leftdPerB->at(fsgrids::dperb::dPERBydx) + rightdPerB->at(fsgrids::dperb::dPERBydx) - 2 * middPerB->at(fsgrids::dperb::dPERBydx)) / (dPerBGrid.DX * dPerBGrid.DX) - 
         (leftdPerB->at(fsgrids::dperb::dPERBxdy) + rightdPerB->at(fsgrids::dperb::dPERBxdy) - 2 * middPerB->at(fsgrids::dperb::dPERBxdy)) / (dPerBGrid.DX * dPerBGrid.DX);
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
   Real hyperres_coeff = -1.0 * 4.0 * M_PI * M_PI *
           dPerBGrid.DX * dPerBGrid.DX *
           Bmag / limitedRhoq / physicalconstants::MU_0;

   // xx-derivative
   EHyperX += hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,0,0);

   // yy-derivative
   EHyperX += hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,0,1);

   // zz-derivative
   EHyperX += hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,0,2);

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
   Real hyperres_coeff = -1.0 * 4.0 * M_PI * M_PI *
           dPerBGrid.DY * dPerBGrid.DY *
           Bmag / limitedRhoq / physicalconstants::MU_0;

   // xx-derivative
   EHyperY += hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,1,0);

   // yy-derivative
   EHyperY += hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,1,1);

   // zz-derivative
   EHyperY += hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,1,2);

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
   Real hyperres_coeff = -1.0 * 4.0 * M_PI * M_PI *
           dPerBGrid.DZ * dPerBGrid.DZ *
           Bmag / limitedRhoq / physicalconstants::MU_0;

   // xx-derivative
   EHyperZ += hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,2,0);

   // yy-derivative
   EHyperZ += hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,2,1);

   // zz-derivative
   EHyperZ += hyperres_coeff * calculateSecondDerivativeOfCurl(dPerBGrid,i,j,k,2,2);

   EHyperGrid.get(i,j,k)->at(fsgrids::ehyper::EZHYPER) = EHyperZ;
}

/** Calculate the electron pressure gradient term on all given cells.
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

   if (P::ohmHallTerm == 0) {
      phiprof::Timer mpiTimer {"Hyper field update ghosts MPI", {"MPI"}};
      dPerBGrid.updateGhostCells();
      mpiTimer.stop();
   }

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
