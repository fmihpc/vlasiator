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

#include <cstdlib>

#include "fs_common.h"
#include "derivatives.hpp"
#include "fs_limiters.h"

/*! \brief Low-level spatial derivatives calculation.
 *
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in project.h.
 * Uses RHO, V[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the second-order method,
 * and RHO_DT2, V[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 *
 * For sysBoundaryLayer 1 or 2, we are near a boundary, and we wish to use regular centered differences instead of slope limiter-adjusted values.
 * This is to minimize oscillations as a smooth behaviour is required near artificial boundaries,
 * unlike at boundaries and shocks inside the simulation domain.
 *
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa calculateDerivativesSimple calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivatives(
   cint i,
   cint j,
   cint k,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {
   std::array<Real, fsgrids::dperb::N_DPERB> * dPerB = dPerBGrid.get(i,j,k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS> * dMoments = dMomentsGrid.get(i,j,k);

   // Get boundary flag for the cell:
   cuint sysBoundaryFlag  = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   cuint sysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;

   // Constants for electron pressure derivatives
   // Upstream pressure
   Real Peupstream = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   Real Peconst = Peupstream * pow(Parameters::electronDensity, -Parameters::electronPTindex);

   std::array<Real, fsgrids::moments::N_MOMENTS> * leftMoments = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD> * leftPerB = NULL;
   std::array<Real, fsgrids::moments::N_MOMENTS> * centMoments = momentsGrid.get(i,j,k);
   std::array<Real, fsgrids::bfield::N_BFIELD> * centPerB = perBGrid.get(i,j,k);
   #ifdef DEBUG_SOLVERS
   if (centMoments->at(fsgrids::moments::RHOM) <= 0) {
      std::cerr << __FILE__ << ":" << __LINE__
         << (centMoments->at(fsgrids::moments::RHOM) < 0 ? " Negative" : " Zero") << " density in spatial cell at (" << i << " " << j << " " << k << ")"
         << std::endl;
      abort();
   }
   #endif
   std::array<Real, fsgrids::moments::N_MOMENTS> * rghtMoments = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * rghtPerB = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * botLeft = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * botRght = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * topLeft = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * topRght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      leftPerB = perBGrid.get(i-1,j,k);
      rghtPerB = perBGrid.get(i+1,j,k);
      leftMoments = momentsGrid.get(i-1,j,k);
      rghtMoments = momentsGrid.get(i+1,j,k);
      #ifdef DEBUG_SOLVERS
      if (leftMoments->at(fsgrids::moments::RHOM) <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (leftMoments->at(fsgrids::moments::RHOM) < 0 ? " Negative" : " Zero") << " density in spatial cell " //<< leftNbrID
            << std::endl;
         abort();
      }
      if (rghtMoments->at(fsgrids::moments::RHOM) <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (rghtMoments->at(fsgrids::moments::RHOM) < 0 ? " Negative" : " Zero") << " density in spatial cell " //<< rightNbrID
            << std::endl;
         abort();
      }
      #endif
      
      if(sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         dMoments->at(fsgrids::dmoments::drhomdx) = (rghtMoments->at(fsgrids::moments::RHOM)-leftMoments->at(fsgrids::moments::RHOM))/2;
         dMoments->at(fsgrids::dmoments::drhoqdx) = (rghtMoments->at(fsgrids::moments::RHOQ)-leftMoments->at(fsgrids::moments::RHOQ))/2;
         dMoments->at(fsgrids::dmoments::dp11dx) = (rghtMoments->at(fsgrids::moments::P_11)-leftMoments->at(fsgrids::moments::P_11))/2;
         dMoments->at(fsgrids::dmoments::dp22dx) = (rghtMoments->at(fsgrids::moments::P_22)-leftMoments->at(fsgrids::moments::P_22))/2;
         dMoments->at(fsgrids::dmoments::dp33dx) = (rghtMoments->at(fsgrids::moments::P_33)-leftMoments->at(fsgrids::moments::P_33))/2;

         dMoments->at(fsgrids::dmoments::dVxdx)  = (rghtMoments->at(fsgrids::moments::VX)-leftMoments->at(fsgrids::moments::VX))/2;
         dMoments->at(fsgrids::dmoments::dVydx)  = (rghtMoments->at(fsgrids::moments::VY)-leftMoments->at(fsgrids::moments::VY))/2;
         dMoments->at(fsgrids::dmoments::dVzdx)  = (rghtMoments->at(fsgrids::moments::VZ)-leftMoments->at(fsgrids::moments::VZ))/2;
         dPerB->at(fsgrids::dperb::dPERBydx)  = (rghtPerB->at(fsgrids::bfield::PERBY)-leftPerB->at(fsgrids::bfield::PERBY))/2;
         dPerB->at(fsgrids::dperb::dPERBzdx)  = (rghtPerB->at(fsgrids::bfield::PERBZ)-leftPerB->at(fsgrids::bfield::PERBZ))/2;
      } else {
         dMoments->at(fsgrids::dmoments::drhomdx) = limiter(leftMoments->at(fsgrids::moments::RHOM),centMoments->at(fsgrids::moments::RHOM),rghtMoments->at(fsgrids::moments::RHOM));
         dMoments->at(fsgrids::dmoments::drhoqdx) = limiter(leftMoments->at(fsgrids::moments::RHOQ),centMoments->at(fsgrids::moments::RHOQ),rghtMoments->at(fsgrids::moments::RHOQ));
         dMoments->at(fsgrids::dmoments::dp11dx) = limiter(leftMoments->at(fsgrids::moments::P_11),centMoments->at(fsgrids::moments::P_11),rghtMoments->at(fsgrids::moments::P_11));
         dMoments->at(fsgrids::dmoments::dp22dx) = limiter(leftMoments->at(fsgrids::moments::P_22),centMoments->at(fsgrids::moments::P_22),rghtMoments->at(fsgrids::moments::P_22));
         dMoments->at(fsgrids::dmoments::dp33dx) = limiter(leftMoments->at(fsgrids::moments::P_33),centMoments->at(fsgrids::moments::P_33),rghtMoments->at(fsgrids::moments::P_33));

         dMoments->at(fsgrids::dmoments::dVxdx)  = limiter(leftMoments->at(fsgrids::moments::VX), centMoments->at(fsgrids::moments::VX), rghtMoments->at(fsgrids::moments::VX));
         dMoments->at(fsgrids::dmoments::dVydx)  = limiter(leftMoments->at(fsgrids::moments::VY), centMoments->at(fsgrids::moments::VY), rghtMoments->at(fsgrids::moments::VY));
         dMoments->at(fsgrids::dmoments::dVzdx)  = limiter(leftMoments->at(fsgrids::moments::VZ), centMoments->at(fsgrids::moments::VZ), rghtMoments->at(fsgrids::moments::VZ));
         dPerB->at(fsgrids::dperb::dPERBydx)  = limiter(leftPerB->at(fsgrids::bfield::PERBY),centPerB->at(fsgrids::bfield::PERBY),rghtPerB->at(fsgrids::bfield::PERBY));
         dPerB->at(fsgrids::dperb::dPERBzdx)  = limiter(leftPerB->at(fsgrids::bfield::PERBZ),centPerB->at(fsgrids::bfield::PERBZ),rghtPerB->at(fsgrids::bfield::PERBZ));
      }

      // pres_e = const * np.power(rho_e, index)
      dMoments->at(fsgrids::dmoments::dPedx) = Peconst * limiter(pow(leftMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex));      
      
      if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
        dPerB->at(fsgrids::dperb::dPERBydxx) = 0.0;
        dPerB->at(fsgrids::dperb::dPERBzdxx) = 0.0;
      } else {
        dPerB->at(fsgrids::dperb::dPERBydxx) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
        dPerB->at(fsgrids::dperb::dPERBzdxx) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
      }
   } else {
      SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 0);
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      leftPerB = perBGrid.get(i,j-1,k);
      rghtPerB = perBGrid.get(i,j+1,k);
      leftMoments = momentsGrid.get(i,j-1,k);
      rghtMoments = momentsGrid.get(i,j+1,k);

      if(sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         dMoments->at(fsgrids::dmoments::drhomdy) = (rghtMoments->at(fsgrids::moments::RHOM)-leftMoments->at(fsgrids::moments::RHOM))/2;
         dMoments->at(fsgrids::dmoments::drhoqdy) = (rghtMoments->at(fsgrids::moments::RHOQ)-leftMoments->at(fsgrids::moments::RHOQ))/2;
         dMoments->at(fsgrids::dmoments::dp11dy) = (rghtMoments->at(fsgrids::moments::P_11)-leftMoments->at(fsgrids::moments::P_11))/2;
         dMoments->at(fsgrids::dmoments::dp22dy) = (rghtMoments->at(fsgrids::moments::P_22)-leftMoments->at(fsgrids::moments::P_22))/2;
         dMoments->at(fsgrids::dmoments::dp33dy) = (rghtMoments->at(fsgrids::moments::P_33)-leftMoments->at(fsgrids::moments::P_33))/2;
         dMoments->at(fsgrids::dmoments::dVxdy)  = (rghtMoments->at(fsgrids::moments::VX)-leftMoments->at(fsgrids::moments::VX))/2;
         dMoments->at(fsgrids::dmoments::dVydy)  = (rghtMoments->at(fsgrids::moments::VY)-leftMoments->at(fsgrids::moments::VY))/2;
         dMoments->at(fsgrids::dmoments::dVzdy)  = (rghtMoments->at(fsgrids::moments::VZ)-leftMoments->at(fsgrids::moments::VZ))/2;

         dPerB->at(fsgrids::dperb::dPERBxdy)  = (rghtPerB->at(fsgrids::bfield::PERBX)-leftPerB->at(fsgrids::bfield::PERBX))/2;
         dPerB->at(fsgrids::dperb::dPERBzdy)  = (rghtPerB->at(fsgrids::bfield::PERBZ)-leftPerB->at(fsgrids::bfield::PERBZ))/2;
      } else {
         dMoments->at(fsgrids::dmoments::drhomdy) = limiter(leftMoments->at(fsgrids::moments::RHOM),centMoments->at(fsgrids::moments::RHOM),rghtMoments->at(fsgrids::moments::RHOM));
         dMoments->at(fsgrids::dmoments::drhoqdy) = limiter(leftMoments->at(fsgrids::moments::RHOQ),centMoments->at(fsgrids::moments::RHOQ),rghtMoments->at(fsgrids::moments::RHOQ));
         dMoments->at(fsgrids::dmoments::dp11dy) = limiter(leftMoments->at(fsgrids::moments::P_11),centMoments->at(fsgrids::moments::P_11),rghtMoments->at(fsgrids::moments::P_11));
         dMoments->at(fsgrids::dmoments::dp22dy) = limiter(leftMoments->at(fsgrids::moments::P_22),centMoments->at(fsgrids::moments::P_22),rghtMoments->at(fsgrids::moments::P_22));
         dMoments->at(fsgrids::dmoments::dp33dy) = limiter(leftMoments->at(fsgrids::moments::P_33),centMoments->at(fsgrids::moments::P_33),rghtMoments->at(fsgrids::moments::P_33));
         dMoments->at(fsgrids::dmoments::dVxdy)  = limiter(leftMoments->at(fsgrids::moments::VX), centMoments->at(fsgrids::moments::VX), rghtMoments->at(fsgrids::moments::VX));
         dMoments->at(fsgrids::dmoments::dVydy)  = limiter(leftMoments->at(fsgrids::moments::VY), centMoments->at(fsgrids::moments::VY), rghtMoments->at(fsgrids::moments::VY));
         dMoments->at(fsgrids::dmoments::dVzdy)  = limiter(leftMoments->at(fsgrids::moments::VZ), centMoments->at(fsgrids::moments::VZ), rghtMoments->at(fsgrids::moments::VZ));

         dPerB->at(fsgrids::dperb::dPERBxdy)  = limiter(leftPerB->at(fsgrids::bfield::PERBX),centPerB->at(fsgrids::bfield::PERBX),rghtPerB->at(fsgrids::bfield::PERBX));
         dPerB->at(fsgrids::dperb::dPERBzdy)  = limiter(leftPerB->at(fsgrids::bfield::PERBZ),centPerB->at(fsgrids::bfield::PERBZ),rghtPerB->at(fsgrids::bfield::PERBZ));
      }

      // pres_e = const * np.power(rho_e, index)
      dMoments->at(fsgrids::dmoments::dPedy) = Peconst * limiter(pow(leftMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex));      

      if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
         dPerB->at(fsgrids::dperb::dPERBxdyy) = 0.0;
         dPerB->at(fsgrids::dperb::dPERBzdyy) = 0.0;
      } else {
         dPerB->at(fsgrids::dperb::dPERBxdyy) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
         dPerB->at(fsgrids::dperb::dPERBzdyy) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
      }
      
   } else {
      SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 1);
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      leftPerB = perBGrid.get(i,j,k-1);
      rghtPerB = perBGrid.get(i,j,k+1);
      leftMoments = momentsGrid.get(i,j,k-1);
      rghtMoments = momentsGrid.get(i,j,k+1);
      if(sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         dMoments->at(fsgrids::dmoments::drhomdz) = (rghtMoments->at(fsgrids::moments::RHOM)-leftMoments->at(fsgrids::moments::RHOM))/2;
         dMoments->at(fsgrids::dmoments::drhoqdz) = (rghtMoments->at(fsgrids::moments::RHOQ)-leftMoments->at(fsgrids::moments::RHOQ))/2;
         dMoments->at(fsgrids::dmoments::dp11dz) = (rghtMoments->at(fsgrids::moments::P_11)-leftMoments->at(fsgrids::moments::P_11))/2;
         dMoments->at(fsgrids::dmoments::dp22dz) = (rghtMoments->at(fsgrids::moments::P_22)-leftMoments->at(fsgrids::moments::P_22))/2;
         dMoments->at(fsgrids::dmoments::dp33dz) = (rghtMoments->at(fsgrids::moments::P_33)-leftMoments->at(fsgrids::moments::P_33))/2;
         dMoments->at(fsgrids::dmoments::dVxdz)  = (rghtMoments->at(fsgrids::moments::VX)-leftMoments->at(fsgrids::moments::VX))/2;
         dMoments->at(fsgrids::dmoments::dVydz)  = (rghtMoments->at(fsgrids::moments::VY)-leftMoments->at(fsgrids::moments::VY))/2;
         dMoments->at(fsgrids::dmoments::dVzdz)  = (rghtMoments->at(fsgrids::moments::VZ)-leftMoments->at(fsgrids::moments::VZ))/2;
         
         dPerB->at(fsgrids::dperb::dPERBxdz)  = (rghtPerB->at(fsgrids::bfield::PERBX)-leftPerB->at(fsgrids::bfield::PERBX))/2;
         dPerB->at(fsgrids::dperb::dPERBydz)  = (rghtPerB->at(fsgrids::bfield::PERBY)-leftPerB->at(fsgrids::bfield::PERBY))/2;
      } else {
         dMoments->at(fsgrids::dmoments::drhomdz) = limiter(leftMoments->at(fsgrids::moments::RHOM),centMoments->at(fsgrids::moments::RHOM),rghtMoments->at(fsgrids::moments::RHOM));
         dMoments->at(fsgrids::dmoments::drhoqdz) = limiter(leftMoments->at(fsgrids::moments::RHOQ),centMoments->at(fsgrids::moments::RHOQ),rghtMoments->at(fsgrids::moments::RHOQ));
         dMoments->at(fsgrids::dmoments::dp11dz) = limiter(leftMoments->at(fsgrids::moments::P_11),centMoments->at(fsgrids::moments::P_11),rghtMoments->at(fsgrids::moments::P_11));
         dMoments->at(fsgrids::dmoments::dp22dz) = limiter(leftMoments->at(fsgrids::moments::P_22),centMoments->at(fsgrids::moments::P_22),rghtMoments->at(fsgrids::moments::P_22));
         dMoments->at(fsgrids::dmoments::dp33dz) = limiter(leftMoments->at(fsgrids::moments::P_33),centMoments->at(fsgrids::moments::P_33),rghtMoments->at(fsgrids::moments::P_33));
         dMoments->at(fsgrids::dmoments::dVxdz)  = limiter(leftMoments->at(fsgrids::moments::VX), centMoments->at(fsgrids::moments::VX), rghtMoments->at(fsgrids::moments::VX));
         dMoments->at(fsgrids::dmoments::dVydz)  = limiter(leftMoments->at(fsgrids::moments::VY), centMoments->at(fsgrids::moments::VY), rghtMoments->at(fsgrids::moments::VY));
         dMoments->at(fsgrids::dmoments::dVzdz)  = limiter(leftMoments->at(fsgrids::moments::VZ), centMoments->at(fsgrids::moments::VZ), rghtMoments->at(fsgrids::moments::VZ));
         
         dPerB->at(fsgrids::dperb::dPERBxdz)  = limiter(leftPerB->at(fsgrids::bfield::PERBX),centPerB->at(fsgrids::bfield::PERBX),rghtPerB->at(fsgrids::bfield::PERBX));
         dPerB->at(fsgrids::dperb::dPERBydz)  = limiter(leftPerB->at(fsgrids::bfield::PERBY),centPerB->at(fsgrids::bfield::PERBY),rghtPerB->at(fsgrids::bfield::PERBY));
      }

      // pres_e = const * np.power(rho_e, index)
      dMoments->at(fsgrids::dmoments::dPedz) = Peconst * limiter(pow(leftMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex));      

      if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
        dPerB->at(fsgrids::dperb::dPERBxdzz) = 0.0;
        dPerB->at(fsgrids::dperb::dPERBydzz) = 0.0;
      } else {
        dPerB->at(fsgrids::dperb::dPERBxdzz) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
        dPerB->at(fsgrids::dperb::dPERBydzz) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
      }
      
   } else {
      SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 2);
   }
   
   if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
      dPerB->at(fsgrids::dperb::dPERBxdyz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBydxz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBzdxy) = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         botLeft = perBGrid.get(i-1,j-1,k);
         botRght = perBGrid.get(i+1,j-1,k);
         topLeft = perBGrid.get(i-1,j+1,k);
         topRght = perBGrid.get(i+1,j+1,k);
         dPerB->at(fsgrids::dperb::dPERBzdxy) = FOURTH * (botLeft->at(fsgrids::bfield::PERBZ) + topRght->at(fsgrids::bfield::PERBZ) - botRght->at(fsgrids::bfield::PERBZ) - topLeft->at(fsgrids::bfield::PERBZ));
      } else {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 3);
      }
      
      // Calculate xz mixed derivatives:
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         botLeft = perBGrid.get(i-1,j,k-1);
         botRght = perBGrid.get(i+1,j,k-1);
         topLeft = perBGrid.get(i-1,j,k+1);
         topRght = perBGrid.get(i+1,j,k+1);
         dPerB->at(fsgrids::dperb::dPERBydxz) = FOURTH * (botLeft->at(fsgrids::bfield::PERBY) + topRght->at(fsgrids::bfield::PERBY) - botRght->at(fsgrids::bfield::PERBY) - topLeft->at(fsgrids::bfield::PERBY));
      } else {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 4);
      }
      
      // Calculate yz mixed derivatives:
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         botLeft = perBGrid.get(i,j-1,k-1);
         botRght = perBGrid.get(i,j+1,k-1);
         topLeft = perBGrid.get(i,j-1,k+1);
         topRght = perBGrid.get(i,j+1,k+1);
         dPerB->at(fsgrids::dperb::dPERBxdyz) = FOURTH * (botLeft->at(fsgrids::bfield::PERBX) + topRght->at(fsgrids::bfield::PERBX) - botRght->at(fsgrids::bfield::PERBX) - topLeft->at(fsgrids::bfield::PERBX));
      } else {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 5);
      }
   }
}


/*! \brief High-level derivative calculation wrapper function.
 * 

 * B has to be updated because after the system boundary update in propagateMagneticFieldSimple there is no consistent state of B yet everywhere.
 * 
 * Then the derivatives are calculated.
 * 
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param perBDt2Grid fsGrid holding the perturbed B quantities at runge-kutta t=0.5
 * \param momentsGrid fsGrid holding the moment quantities
 * \param momentsDt2Grid fsGrid holding the moment quantities at runge-kutta t=0.5
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param communicateMoments If true, the derivatives of moments (rho, V, P) are communicated to neighbours.
 
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivativesSimple(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase,
   const bool communicateMoments) {
   int timer;
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const int* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   
   phiprof::start("Calculate face derivatives");
   
   timer=phiprof::initializeTimer("MPI","MPI");
   phiprof::start(timer);
   
   switch (RKCase) {
    case RK_ORDER1:
      // Means initialising the solver as well as RK_ORDER1
      // standard case Exchange PERB* with neighbours
      // The update of PERB[XYZ] is needed after the system
      // boundary update of propagateMagneticFieldSimple.
       perBGrid.updateGhostCells();
       if(communicateMoments) {
         momentsGrid.updateGhostCells();
       }
       break;
    case RK_ORDER2_STEP1:
      // Exchange PERB*_DT2,RHO_DT2,V*_DT2 with neighbours The
      // update of PERB[XYZ]_DT2 is needed after the system
      // boundary update of propagateMagneticFieldSimple.
       perBDt2Grid.updateGhostCells();
       if(communicateMoments) {
         momentsDt2Grid.updateGhostCells();
       }
       break;
    case RK_ORDER2_STEP2:
      // Exchange PERB*,RHO,V* with neighbours The update of B
      // is needed after the system boundary update of
      // propagateMagneticFieldSimple.
       perBGrid.updateGhostCells();
       if(communicateMoments) {
         momentsGrid.updateGhostCells();
       }
      break;
    default:
      cerr << __FILE__ << ":" << __LINE__ << " Went through switch, this should not happen." << endl;
      abort();
   }
   
   phiprof::stop(timer);

   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);

   // Calculate derivatives
   #pragma omp parallel for collapse(3)
   for (int k=0; k<gridDims[2]; k++) {
      for (int j=0; j<gridDims[1]; j++) {
         for (int i=0; i<gridDims[0]; i++) {
            if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
            if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
               calculateDerivatives(i,j,k, perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RKCase);
            } else {
               calculateDerivatives(i,j,k, perBDt2Grid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RKCase);
            }
         }
      }
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   phiprof::stop("Calculate face derivatives",N_cells,"Spatial Cells");   
}

/*! \brief Low-level spatial derivatives calculation.
 *
 * Calculate the spatial derivatives of BVOL or set them to zero.
 *
 * For sysBoundaryLayer 1 or 2, we are near a boundary, and we wish to use regular centered differences instead of slope limiter-adjusted values.
 * This is to minimize oscillations as a smooth behaviour is required near artificial boundaries,
 * unlike at boundaries and shocks inside the simulation domain.
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateBVOLDerivatives(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   cint i,
   cint j,
   cint k,
   SysBoundary& sysBoundaries
) {
   std::array<Real, fsgrids::volfields::N_VOL> * array = volGrid.get(i,j,k);
   
   std::array<Real, fsgrids::volfields::N_VOL> * left = NULL;
   std::array<Real, fsgrids::volfields::N_VOL> * rght = NULL;
   
   cuint sysBoundaryFlag = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   cuint sysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;

   // Calculate x-derivatives (is not TVD for AMR mesh):
   if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY || sysBoundaryLayer == 1) {

      left = volGrid.get(i-1,j,k);
      rght = volGrid.get(i+1,j,k);
      
      if (sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         array->at(fsgrids::volfields::dPERBYVOLdx) = (rght->at(fsgrids::volfields::PERBYVOL)-left->at(fsgrids::volfields::PERBYVOL))/2;
         array->at(fsgrids::volfields::dPERBZVOLdx) = (rght->at(fsgrids::volfields::PERBZVOL)-left->at(fsgrids::volfields::PERBZVOL))/2;
      } else {
         array->at(fsgrids::volfields::dPERBYVOLdx) = limiter(left->at(fsgrids::volfields::PERBYVOL),array->at(fsgrids::volfields::PERBYVOL),rght->at(fsgrids::volfields::PERBYVOL));
         array->at(fsgrids::volfields::dPERBZVOLdx) = limiter(left->at(fsgrids::volfields::PERBZVOL),array->at(fsgrids::volfields::PERBZVOL),rght->at(fsgrids::volfields::PERBZVOL));
      }
   } else {
      SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 0);
   }
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
   if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY || sysBoundaryLayer == 1) {
      left = volGrid.get(i,j-1,k);
      rght = volGrid.get(i,j+1,k);
      
      if (sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         array->at(fsgrids::volfields::dPERBXVOLdy) = (rght->at(fsgrids::volfields::PERBXVOL)-left->at(fsgrids::volfields::PERBXVOL))/2;
         array->at(fsgrids::volfields::dPERBZVOLdy) = (rght->at(fsgrids::volfields::PERBZVOL)-left->at(fsgrids::volfields::PERBZVOL))/2;
      } else {
         array->at(fsgrids::volfields::dPERBXVOLdy) = limiter(left->at(fsgrids::volfields::PERBXVOL),array->at(fsgrids::volfields::PERBXVOL),rght->at(fsgrids::volfields::PERBXVOL));
         array->at(fsgrids::volfields::dPERBZVOLdy) = limiter(left->at(fsgrids::volfields::PERBZVOL),array->at(fsgrids::volfields::PERBZVOL),rght->at(fsgrids::volfields::PERBZVOL));
      }
   } else {
      SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 1);
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY || sysBoundaryLayer == 1) {
      left = volGrid.get(i,j,k-1);
      rght = volGrid.get(i,j,k+1);
      
      if (sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         array->at(fsgrids::volfields::dPERBXVOLdz) = (rght->at(fsgrids::volfields::PERBXVOL)-left->at(fsgrids::volfields::PERBXVOL))/2;
         array->at(fsgrids::volfields::dPERBYVOLdz) = (rght->at(fsgrids::volfields::PERBYVOL)-left->at(fsgrids::volfields::PERBYVOL))/2;
      } else {
         array->at(fsgrids::volfields::dPERBXVOLdz) = limiter(left->at(fsgrids::volfields::PERBXVOL),array->at(fsgrids::volfields::PERBXVOL),rght->at(fsgrids::volfields::PERBXVOL));
         array->at(fsgrids::volfields::dPERBYVOLdz) = limiter(left->at(fsgrids::volfields::PERBYVOL),array->at(fsgrids::volfields::PERBYVOL),rght->at(fsgrids::volfields::PERBYVOL));
      }
   } else {
      SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 2);
   }
}

/*! \brief High-level derivative calculation wrapper function.
 * 
 * BVOL has been calculated locally by calculateVolumeAveragedFields but not communicated.
 * For the acceleration step one needs the cross-derivatives of BVOL
 * 
 * \param volGrid fsGrid holding the volume averaged fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * 
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateBVOLDerivativesSimple(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries
) {
   int timer;
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const int* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   
   phiprof::start("Calculate volume derivatives");
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);
   volGrid.updateGhostCells();
   
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   
   // Calculate derivatives
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);
   
   #pragma omp parallel for collapse(3)
   for (int k=0; k<gridDims[2]; k++) {
      for (int j=0; j<gridDims[1]; j++) {
         for (int i=0; i<gridDims[0]; i++) {
            if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               continue;
            }
            calculateBVOLDerivatives(volGrid,technicalGrid,i,j,k,sysBoundaries);
         }
      }
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");

   phiprof::stop("Calculate volume derivatives",N_cells,"Spatial Cells");
}

/*! \brief Returns volumetric E of cell
 *
 */
static std::array<Real, 3> getE(SpatialCell* cell)
{
   return std::array<Real, 3> { {cell->parameters[CellParams::EXVOL], cell->parameters[CellParams::EYVOL], cell->parameters[CellParams::EZVOL]} };
}

/*! \brief Returns perturbed volumetric B of cell
 *
 */
static std::array<Real, 3> getPerB(SpatialCell* cell)
{
   return std::array<Real, 3> { {cell->parameters[CellParams::PERBXVOL], cell->parameters[CellParams::PERBYVOL], cell->parameters[CellParams::PERBZVOL]} };
}

/*! \brief Returns volumetric B of cell
 *
 */
static std::array<Real, 3> getB(SpatialCell* cell)
{
   return std::array<Real, 3> { 
      {
         cell->parameters[CellParams::BGBXVOL] + cell->parameters[CellParams::PERBXVOL], 
         cell->parameters[CellParams::BGBYVOL] + cell->parameters[CellParams::PERBYVOL], 
         cell->parameters[CellParams::BGBZVOL] + cell->parameters[CellParams::PERBZVOL]
      } 
   };
}

/*! \brief Calculates momentum density of cell
 *
 */
static std::array<Real, 3> getMomentumDensity(SpatialCell* cell)
{
   Real rho = cell->parameters[CellParams::RHOM];
   return std::array<Real, 3> { {rho * cell->parameters[CellParams::VX], rho * cell->parameters[CellParams::VY], rho * cell->parameters[CellParams::VZ]} };
}

/*! \brief Calculates energy density for spatial cell with only perturbated magnetic field
 *
 */
static Real calculateU1(SpatialCell* cell)
{
   std::array<Real, 3> p = getMomentumDensity(cell);
   std::array<Real, 3> B = getPerB(cell);
   return (pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2)) / (2.0 * cell->parameters[CellParams::RHOM]) + (pow(B[0], 2) + pow(B[1], 2) + pow(B[2], 2)) / (2.0 * physicalconstants::MU_0);
}

/*! \brief Low-level scaled gradients calculation
 * 
 * For the SpatialCell* cell and its neighbors, calculate scaled gradients and their maximum alpha
 * The gradients are the same as in the GUMICS simulation, see
 * Janhunen, P., Palmroth, M., Laitinen, T., Honkonen, I., Juusola, L., Facsko, G., & Pulkkinen, T. I. (2012). The GUMICS-4 global MHD magnetosphere-ionosphere coupling simulation. Journal of Atmospheric and Solar - Terrestrial Physics, 80, 48-59. https://doi.org/10.1016/j.jastp.2012.03.006
 *
 */
void calculateScaledDeltas(
   SpatialCell* cell,
   std::vector<SpatialCell*>& neighbors)
{
   Real dRho {0};
   Real dU {0};
   Real dPsq {0};
   Real dBsq {0};
   Real dB {0};

   Real myRho {cell->parameters[CellParams::RHOM]};
   Real myU {calculateU1(cell)};
   std::array<Real, 3> myP = getMomentumDensity(cell);
   std::array<Real, 3> myB = getPerB(cell);
   for (SpatialCell* neighbor : neighbors) {
      Real otherRho = neighbor->parameters[CellParams::RHOM];
      Real otherU = calculateU1(neighbor);
      std::array<Real, 3> otherP = getMomentumDensity(neighbor);
      std::array<Real, 3> otherB = getPerB(neighbor);
      Real deltaBsq = pow(myB[0] - otherB[0], 2) + pow(myB[1] - otherB[1], 2) + pow(myB[2] - otherB[2], 2);

      // Assignment intentional
      if (Real maxRho = std::max(myRho, otherRho)) {
         dRho = std::max(fabs(myRho - otherRho) / maxRho, dRho);
      }
      if (Real maxU = std::max(myU, otherU)) {
         dU = std::max(fabs(myU - otherU) / maxU, dU);
         dPsq = std::max((pow(myP[0] - otherP[0], 2) + pow(myP[1] - otherP[1], 2) + pow(myP[2] - otherP[2], 2)) / (2 * myRho * maxU), dPsq) / 4.0;
         dBsq = std::max(deltaBsq / (2 * physicalconstants::MU_0 * maxU), dBsq) / 4.0;
      }
      if(Real maxB = sqrt(std::max(pow(myB[0], 2) + pow(myB[1], 2) + pow(myB[2], 2), pow(otherB[0], 2) + pow(otherB[1], 2) + pow(otherB[2], 2)))) {
         dB = std::max(sqrt(deltaBsq) / maxB, dB) / 2.0;
      }
   }
   
   Real alpha = dRho;
   if (dU > alpha) {
      alpha = dU;
   }
   if (dPsq > alpha) {
      alpha = dPsq;
   }
   if (dBsq > alpha) {
      alpha = dBsq;
   }
   if (dB > alpha) {
      alpha = dB;
   }

   Real dBXdy {cell->derivativesBVOL[bvolderivatives::dPERBXVOLdy]};
   Real dBXdz {cell->derivativesBVOL[bvolderivatives::dPERBXVOLdz]};
   Real dBYdx {cell->derivativesBVOL[bvolderivatives::dPERBYVOLdx]};
   Real dBYdz {cell->derivativesBVOL[bvolderivatives::dPERBYVOLdz]};
   Real dBZdx {cell->derivativesBVOL[bvolderivatives::dPERBZVOLdx]};
   Real dBZdy {cell->derivativesBVOL[bvolderivatives::dPERBZVOLdy]};

   // Note missing factor of mu_0, since we want B and J in same units later
   std::array<Real, 3> myJ = {dBZdy - dBYdz, dBXdz - dBZdx, dBYdx - dBXdy};
   Real BdotJ {0.0};
   Real Bsq {0.0};
   for (int i = 0; i < 3; ++i) {
      BdotJ += myB[i] * myJ[i];
      Bsq += myB[i] * myB[i];
   }

   Real Bperp {0.0};
   Real J {0.0};
   for (int i = 0; i < 3; ++i) {
      Bperp += std::pow(myB[i] * (1 - BdotJ / Bsq), 2);
      J += myJ[i] * myJ[i];
   }
   Bperp = std::sqrt(Bperp);
   J = std::sqrt(J);

   cell->parameters[CellParams::AMR_DRHO] = dRho;
   cell->parameters[CellParams::AMR_DU] = dU;
   cell->parameters[CellParams::AMR_DPSQ] = dPsq;
   cell->parameters[CellParams::AMR_DBSQ] = dBsq;
   cell->parameters[CellParams::AMR_DB] = dB;
   cell->parameters[CellParams::AMR_ALPHA] = alpha;
   cell->parameters[CellParams::AMR_JPERB] = J / Bperp;
}

/*! \brief High-level scaled gradient calculation wrapper function.
 * 
 * Calculates gradients needed for alpha everywhere in the grid
 * 
 */

void calculateScaledDeltasSimple(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid)
{
   const vector<CellID>& cells = getLocalCells();
   int N_cells = cells.size();
   int timer;
   phiprof::start("Calculate volume gradients");
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);

   // We only need nearest neighbourhood and spatial data here
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);
   mpiGrid.update_copies_of_remote_neighbors(NEAREST_NEIGHBORHOOD_ID);
   
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   // Calculate derivatives
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);
   
   #pragma omp parallel for
   for (uint i = 0; i < cells.size(); ++i) {
   //for (CellID id : cells) {
      CellID id = cells[i];
      SpatialCell* cell = mpiGrid[id];
      std::vector<SpatialCell*> neighbors;
      for (auto neighPair : mpiGrid.get_face_neighbors_of(id)) {
         neighbors.push_back(mpiGrid[neighPair.first]);
      }
      calculateScaledDeltas(cell, neighbors);
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");

   phiprof::stop("Calculate volume gradients",N_cells,"Spatial Cells");
}
