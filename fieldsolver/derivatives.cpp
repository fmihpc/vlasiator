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
#include <Eigen/Geometry>

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
 * \param calculateMoments Bool telling whether the derivatives for moments need updating too.
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
   const bool calculateMoments
) {
   std::array<Real, fsgrids::dperb::N_DPERB> * dPerB = dPerBGrid.get(i,j,k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS> * dMoments = dMomentsGrid.get(i,j,k);

   // Get boundary flag for the cell:
   cuint sysBoundaryFlag  = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   cuint sysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;

   // Constants for electron pressure derivatives, see also ldz_gradpe.cpp
   // Calculate anchor point constants: First the pressure, then a derived constant.
   Real Pe_anchor = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   Real Pe_const = Pe_anchor * pow(Parameters::electronDensity, -Parameters::electronPTindex);

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
   leftPerB = perBGrid.get(i-1,j,k);
   rghtPerB = perBGrid.get(i+1,j,k);
   leftMoments = momentsGrid.get(i-1,j,k);
   rghtMoments = momentsGrid.get(i+1,j,k);
   if(leftPerB == NULL) {
      leftPerB = centPerB;
      leftMoments = centMoments;
   }
   if(rghtPerB == NULL) {
      rghtPerB = centPerB;
      rghtMoments = centMoments;
   }
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

   if(P::fieldSolverFiniteDifferencingAtBoundaries && (sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY))) {
      if(calculateMoments) {
         dMoments->at(fsgrids::dmoments::drhomdx) = (rghtMoments->at(fsgrids::moments::RHOM)-leftMoments->at(fsgrids::moments::RHOM))/2;
         dMoments->at(fsgrids::dmoments::drhoqdx) = (rghtMoments->at(fsgrids::moments::RHOQ)-leftMoments->at(fsgrids::moments::RHOQ))/2;
         dMoments->at(fsgrids::dmoments::dp11dx) = (rghtMoments->at(fsgrids::moments::P_11)-leftMoments->at(fsgrids::moments::P_11))/2;
         dMoments->at(fsgrids::dmoments::dp22dx) = (rghtMoments->at(fsgrids::moments::P_22)-leftMoments->at(fsgrids::moments::P_22))/2;
         dMoments->at(fsgrids::dmoments::dp33dx) = (rghtMoments->at(fsgrids::moments::P_33)-leftMoments->at(fsgrids::moments::P_33))/2;
         dMoments->at(fsgrids::dmoments::dVxdx)  = (rghtMoments->at(fsgrids::moments::VX)-leftMoments->at(fsgrids::moments::VX))/2;
         dMoments->at(fsgrids::dmoments::dVydx)  = (rghtMoments->at(fsgrids::moments::VY)-leftMoments->at(fsgrids::moments::VY))/2;
         dMoments->at(fsgrids::dmoments::dVzdx)  = (rghtMoments->at(fsgrids::moments::VZ)-leftMoments->at(fsgrids::moments::VZ))/2;
      }
      dPerB->at(fsgrids::dperb::dPERBydx)  = (rghtPerB->at(fsgrids::bfield::PERBY)-leftPerB->at(fsgrids::bfield::PERBY))/2;
      dPerB->at(fsgrids::dperb::dPERBzdx)  = (rghtPerB->at(fsgrids::bfield::PERBZ)-leftPerB->at(fsgrids::bfield::PERBZ))/2;
   } else {
      if(calculateMoments) {
         dMoments->at(fsgrids::dmoments::drhomdx) = limiter(leftMoments->at(fsgrids::moments::RHOM),centMoments->at(fsgrids::moments::RHOM),rghtMoments->at(fsgrids::moments::RHOM));
         dMoments->at(fsgrids::dmoments::drhoqdx) = limiter(leftMoments->at(fsgrids::moments::RHOQ),centMoments->at(fsgrids::moments::RHOQ),rghtMoments->at(fsgrids::moments::RHOQ));
         dMoments->at(fsgrids::dmoments::dp11dx) = limiter(leftMoments->at(fsgrids::moments::P_11),centMoments->at(fsgrids::moments::P_11),rghtMoments->at(fsgrids::moments::P_11));
         dMoments->at(fsgrids::dmoments::dp22dx) = limiter(leftMoments->at(fsgrids::moments::P_22),centMoments->at(fsgrids::moments::P_22),rghtMoments->at(fsgrids::moments::P_22));
         dMoments->at(fsgrids::dmoments::dp33dx) = limiter(leftMoments->at(fsgrids::moments::P_33),centMoments->at(fsgrids::moments::P_33),rghtMoments->at(fsgrids::moments::P_33));
         dMoments->at(fsgrids::dmoments::dVxdx)  = limiter(leftMoments->at(fsgrids::moments::VX), centMoments->at(fsgrids::moments::VX), rghtMoments->at(fsgrids::moments::VX));
         dMoments->at(fsgrids::dmoments::dVydx)  = limiter(leftMoments->at(fsgrids::moments::VY), centMoments->at(fsgrids::moments::VY), rghtMoments->at(fsgrids::moments::VY));
         dMoments->at(fsgrids::dmoments::dVzdx)  = limiter(leftMoments->at(fsgrids::moments::VZ), centMoments->at(fsgrids::moments::VZ), rghtMoments->at(fsgrids::moments::VZ));
      }
      dPerB->at(fsgrids::dperb::dPERBydx)  = limiter(leftPerB->at(fsgrids::bfield::PERBY),centPerB->at(fsgrids::bfield::PERBY),rghtPerB->at(fsgrids::bfield::PERBY));
      dPerB->at(fsgrids::dperb::dPERBzdx)  = limiter(leftPerB->at(fsgrids::bfield::PERBZ),centPerB->at(fsgrids::bfield::PERBZ),rghtPerB->at(fsgrids::bfield::PERBZ));
   }

   if(calculateMoments) {
      // pres_e = const * np.power(rho_e, index)
      dMoments->at(fsgrids::dmoments::dPedx) = Pe_const * limiter(pow(leftMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex));
   }

   if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
      dPerB->at(fsgrids::dperb::dPERBydxx) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBzdxx) = 0.0;
   } else {
      dPerB->at(fsgrids::dperb::dPERBydxx) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
      dPerB->at(fsgrids::dperb::dPERBzdxx) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
   leftPerB = perBGrid.get(i,j-1,k);
   rghtPerB = perBGrid.get(i,j+1,k);
   leftMoments = momentsGrid.get(i,j-1,k);
   rghtMoments = momentsGrid.get(i,j+1,k);
   if(leftPerB == NULL) {
      leftPerB = centPerB;
      leftMoments = centMoments;
   }
   if(rghtPerB == NULL) {
      rghtPerB = centPerB;
      rghtMoments = centMoments;
   }

   if(P::fieldSolverFiniteDifferencingAtBoundaries && (sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY))) {
      if(calculateMoments) {
         dMoments->at(fsgrids::dmoments::drhomdy) = (rghtMoments->at(fsgrids::moments::RHOM)-leftMoments->at(fsgrids::moments::RHOM))/2;
         dMoments->at(fsgrids::dmoments::drhoqdy) = (rghtMoments->at(fsgrids::moments::RHOQ)-leftMoments->at(fsgrids::moments::RHOQ))/2;
         dMoments->at(fsgrids::dmoments::dp11dy) = (rghtMoments->at(fsgrids::moments::P_11)-leftMoments->at(fsgrids::moments::P_11))/2;
         dMoments->at(fsgrids::dmoments::dp22dy) = (rghtMoments->at(fsgrids::moments::P_22)-leftMoments->at(fsgrids::moments::P_22))/2;
         dMoments->at(fsgrids::dmoments::dp33dy) = (rghtMoments->at(fsgrids::moments::P_33)-leftMoments->at(fsgrids::moments::P_33))/2;
         dMoments->at(fsgrids::dmoments::dVxdy)  = (rghtMoments->at(fsgrids::moments::VX)-leftMoments->at(fsgrids::moments::VX))/2;
         dMoments->at(fsgrids::dmoments::dVydy)  = (rghtMoments->at(fsgrids::moments::VY)-leftMoments->at(fsgrids::moments::VY))/2;
         dMoments->at(fsgrids::dmoments::dVzdy)  = (rghtMoments->at(fsgrids::moments::VZ)-leftMoments->at(fsgrids::moments::VZ))/2;
      }
      dPerB->at(fsgrids::dperb::dPERBxdy)  = (rghtPerB->at(fsgrids::bfield::PERBX)-leftPerB->at(fsgrids::bfield::PERBX))/2;
      dPerB->at(fsgrids::dperb::dPERBzdy)  = (rghtPerB->at(fsgrids::bfield::PERBZ)-leftPerB->at(fsgrids::bfield::PERBZ))/2;
   } else {
      if(calculateMoments) {
         dMoments->at(fsgrids::dmoments::drhomdy) = limiter(leftMoments->at(fsgrids::moments::RHOM),centMoments->at(fsgrids::moments::RHOM),rghtMoments->at(fsgrids::moments::RHOM));
         dMoments->at(fsgrids::dmoments::drhoqdy) = limiter(leftMoments->at(fsgrids::moments::RHOQ),centMoments->at(fsgrids::moments::RHOQ),rghtMoments->at(fsgrids::moments::RHOQ));
         dMoments->at(fsgrids::dmoments::dp11dy) = limiter(leftMoments->at(fsgrids::moments::P_11),centMoments->at(fsgrids::moments::P_11),rghtMoments->at(fsgrids::moments::P_11));
         dMoments->at(fsgrids::dmoments::dp22dy) = limiter(leftMoments->at(fsgrids::moments::P_22),centMoments->at(fsgrids::moments::P_22),rghtMoments->at(fsgrids::moments::P_22));
         dMoments->at(fsgrids::dmoments::dp33dy) = limiter(leftMoments->at(fsgrids::moments::P_33),centMoments->at(fsgrids::moments::P_33),rghtMoments->at(fsgrids::moments::P_33));
         dMoments->at(fsgrids::dmoments::dVxdy)  = limiter(leftMoments->at(fsgrids::moments::VX), centMoments->at(fsgrids::moments::VX), rghtMoments->at(fsgrids::moments::VX));
         dMoments->at(fsgrids::dmoments::dVydy)  = limiter(leftMoments->at(fsgrids::moments::VY), centMoments->at(fsgrids::moments::VY), rghtMoments->at(fsgrids::moments::VY));
         dMoments->at(fsgrids::dmoments::dVzdy)  = limiter(leftMoments->at(fsgrids::moments::VZ), centMoments->at(fsgrids::moments::VZ), rghtMoments->at(fsgrids::moments::VZ));
      }
      dPerB->at(fsgrids::dperb::dPERBxdy)  = limiter(leftPerB->at(fsgrids::bfield::PERBX),centPerB->at(fsgrids::bfield::PERBX),rghtPerB->at(fsgrids::bfield::PERBX));
      dPerB->at(fsgrids::dperb::dPERBzdy)  = limiter(leftPerB->at(fsgrids::bfield::PERBZ),centPerB->at(fsgrids::bfield::PERBZ),rghtPerB->at(fsgrids::bfield::PERBZ));
   }

   // pres_e = const * np.power(rho_e, index)
   dMoments->at(fsgrids::dmoments::dPedy) = Pe_const * limiter(pow(leftMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex));
   if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
      dPerB->at(fsgrids::dperb::dPERBxdyy) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBzdyy) = 0.0;
   } else {
      dPerB->at(fsgrids::dperb::dPERBxdyy) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
      dPerB->at(fsgrids::dperb::dPERBzdyy) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
   }

   // Calculate z-derivatives (is not TVD for AMR mesh):
   leftPerB = perBGrid.get(i,j,k-1);
   rghtPerB = perBGrid.get(i,j,k+1);
   leftMoments = momentsGrid.get(i,j,k-1);
   rghtMoments = momentsGrid.get(i,j,k+1);
   if(leftPerB == NULL) {
      leftPerB = centPerB;
      leftMoments = centMoments;
   }
   if(rghtPerB == NULL) {
      rghtPerB = centPerB;
      rghtMoments = centMoments;
   }

   if(P::fieldSolverFiniteDifferencingAtBoundaries && (sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY))) {
      if(calculateMoments) {
         dMoments->at(fsgrids::dmoments::drhomdz) = (rghtMoments->at(fsgrids::moments::RHOM)-leftMoments->at(fsgrids::moments::RHOM))/2;
         dMoments->at(fsgrids::dmoments::drhoqdz) = (rghtMoments->at(fsgrids::moments::RHOQ)-leftMoments->at(fsgrids::moments::RHOQ))/2;
         dMoments->at(fsgrids::dmoments::dp11dz) = (rghtMoments->at(fsgrids::moments::P_11)-leftMoments->at(fsgrids::moments::P_11))/2;
         dMoments->at(fsgrids::dmoments::dp22dz) = (rghtMoments->at(fsgrids::moments::P_22)-leftMoments->at(fsgrids::moments::P_22))/2;
         dMoments->at(fsgrids::dmoments::dp33dz) = (rghtMoments->at(fsgrids::moments::P_33)-leftMoments->at(fsgrids::moments::P_33))/2;
         dMoments->at(fsgrids::dmoments::dVxdz)  = (rghtMoments->at(fsgrids::moments::VX)-leftMoments->at(fsgrids::moments::VX))/2;
         dMoments->at(fsgrids::dmoments::dVydz)  = (rghtMoments->at(fsgrids::moments::VY)-leftMoments->at(fsgrids::moments::VY))/2;
         dMoments->at(fsgrids::dmoments::dVzdz)  = (rghtMoments->at(fsgrids::moments::VZ)-leftMoments->at(fsgrids::moments::VZ))/2;
      }
      dPerB->at(fsgrids::dperb::dPERBxdz)  = (rghtPerB->at(fsgrids::bfield::PERBX)-leftPerB->at(fsgrids::bfield::PERBX))/2;
      dPerB->at(fsgrids::dperb::dPERBydz)  = (rghtPerB->at(fsgrids::bfield::PERBY)-leftPerB->at(fsgrids::bfield::PERBY))/2;
   } else {
      if(calculateMoments) {
         dMoments->at(fsgrids::dmoments::drhomdz) = limiter(leftMoments->at(fsgrids::moments::RHOM),centMoments->at(fsgrids::moments::RHOM),rghtMoments->at(fsgrids::moments::RHOM));
         dMoments->at(fsgrids::dmoments::drhoqdz) = limiter(leftMoments->at(fsgrids::moments::RHOQ),centMoments->at(fsgrids::moments::RHOQ),rghtMoments->at(fsgrids::moments::RHOQ));
         dMoments->at(fsgrids::dmoments::dp11dz) = limiter(leftMoments->at(fsgrids::moments::P_11),centMoments->at(fsgrids::moments::P_11),rghtMoments->at(fsgrids::moments::P_11));
         dMoments->at(fsgrids::dmoments::dp22dz) = limiter(leftMoments->at(fsgrids::moments::P_22),centMoments->at(fsgrids::moments::P_22),rghtMoments->at(fsgrids::moments::P_22));
         dMoments->at(fsgrids::dmoments::dp33dz) = limiter(leftMoments->at(fsgrids::moments::P_33),centMoments->at(fsgrids::moments::P_33),rghtMoments->at(fsgrids::moments::P_33));
         dMoments->at(fsgrids::dmoments::dVxdz)  = limiter(leftMoments->at(fsgrids::moments::VX), centMoments->at(fsgrids::moments::VX), rghtMoments->at(fsgrids::moments::VX));
         dMoments->at(fsgrids::dmoments::dVydz)  = limiter(leftMoments->at(fsgrids::moments::VY), centMoments->at(fsgrids::moments::VY), rghtMoments->at(fsgrids::moments::VY));
         dMoments->at(fsgrids::dmoments::dVzdz)  = limiter(leftMoments->at(fsgrids::moments::VZ), centMoments->at(fsgrids::moments::VZ), rghtMoments->at(fsgrids::moments::VZ));
      }
      dPerB->at(fsgrids::dperb::dPERBxdz)  = limiter(leftPerB->at(fsgrids::bfield::PERBX),centPerB->at(fsgrids::bfield::PERBX),rghtPerB->at(fsgrids::bfield::PERBX));
      dPerB->at(fsgrids::dperb::dPERBydz)  = limiter(leftPerB->at(fsgrids::bfield::PERBY),centPerB->at(fsgrids::bfield::PERBY),rghtPerB->at(fsgrids::bfield::PERBY));
   }

   if(calculateMoments) {
      // pres_e = const * np.power(rho_e, index)
      dMoments->at(fsgrids::dmoments::dPedz) = Pe_const * limiter(pow(leftMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex));
   }

   if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
      dPerB->at(fsgrids::dperb::dPERBxdzz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBydzz) = 0.0;
   } else {
      dPerB->at(fsgrids::dperb::dPERBxdzz) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
      dPerB->at(fsgrids::dperb::dPERBydzz) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
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
 * \param dMomentsGridDt2 fsGrid holding the derviatives of moments at runge-kutta t=0.5
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
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsDt2Grid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase,
   const bool doMoments) {
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const FsGridTools::FsIndex_t* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::Timer derivativesTimer {"Calculate face derivatives"};
   int computeTimerId {phiprof::initializeTimer("FS derivatives compute cells")};

   phiprof::Timer mpiTimer {"FS derivatives ghost updates MPI", {"MPI"}};
   switch (RKCase) {
    case RK_ORDER1:
      // Means initialising the solver as well as RK_ORDER1
      // standard case Exchange PERB* with neighbours
      // The update of PERB[XYZ] is needed after the system
      // boundary update of propagateMagneticFieldSimple.
       perBGrid.updateGhostCells();
       if(doMoments) {
         momentsGrid.updateGhostCells();
       }
       break;
    case RK_ORDER2_STEP1:
      // Exchange PERB*_DT2,RHO_DT2,V*_DT2 with neighbours The
      // update of PERB[XYZ]_DT2 is needed after the system
      // boundary update of propagateMagneticFieldSimple.
       perBDt2Grid.updateGhostCells();
       if(doMoments) {
         momentsDt2Grid.updateGhostCells();
       }
       break;
    case RK_ORDER2_STEP2:
      // Exchange PERB*,RHO,V* with neighbours The update of B
      // is needed after the system boundary update of
      // propagateMagneticFieldSimple.
       perBGrid.updateGhostCells();
       if(doMoments) {
         momentsGrid.updateGhostCells();
       }
      break;
    default:
      cerr << __FILE__ << ":" << __LINE__ << " Went through switch, this should not happen." << endl;
      abort();
   }
   mpiTimer.stop();

   // Calculate derivatives
   #pragma omp parallel
   {
      phiprof::Timer computeTimer {computeTimerId};
      #pragma omp for collapse(2)
      for (FsGridTools::FsIndex_t k=0; k<gridDims[2]; k++) {
         for (FsGridTools::FsIndex_t j=0; j<gridDims[1]; j++) {
            for (FsGridTools::FsIndex_t i=0; i<gridDims[0]; i++) {
               if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
                  calculateDerivatives(i,j,k, perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, doMoments);
               } else {
                  calculateDerivatives(i,j,k, perBDt2Grid, momentsDt2Grid, dPerBGrid, dMomentsDt2Grid, technicalGrid, sysBoundaries, doMoments);
               }
            }
         }
      }
      computeTimer.stop(N_cells, "Spatial Cells");
   }

   derivativesTimer.stop(N_cells, "Spatial Cells");
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
   left = volGrid.get(i-1,j,k);
   rght = volGrid.get(i+1,j,k);

   if(left == NULL) {
      left = array;
   }
   if(rght == NULL) {
      rght = array;
   }

   if(P::fieldSolverFiniteDifferencingAtBoundaries && (sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY))) {
      array->at(fsgrids::volfields::dPERBXVOLdx) = (rght->at(fsgrids::volfields::PERBXVOL)-left->at(fsgrids::volfields::PERBXVOL))/2;
      array->at(fsgrids::volfields::dPERBYVOLdx) = (rght->at(fsgrids::volfields::PERBYVOL)-left->at(fsgrids::volfields::PERBYVOL))/2;
      array->at(fsgrids::volfields::dPERBZVOLdx) = (rght->at(fsgrids::volfields::PERBZVOL)-left->at(fsgrids::volfields::PERBZVOL))/2;
   } else {
      array->at(fsgrids::volfields::dPERBXVOLdx) = limiter(left->at(fsgrids::volfields::PERBXVOL),array->at(fsgrids::volfields::PERBXVOL),rght->at(fsgrids::volfields::PERBXVOL));
      array->at(fsgrids::volfields::dPERBYVOLdx) = limiter(left->at(fsgrids::volfields::PERBYVOL),array->at(fsgrids::volfields::PERBYVOL),rght->at(fsgrids::volfields::PERBYVOL));
      array->at(fsgrids::volfields::dPERBZVOLdx) = limiter(left->at(fsgrids::volfields::PERBZVOL),array->at(fsgrids::volfields::PERBZVOL),rght->at(fsgrids::volfields::PERBZVOL));
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
   left = volGrid.get(i,j-1,k);
   rght = volGrid.get(i,j+1,k);

   if(left == NULL) {
      left = array;
   }
   if(rght == NULL) {
      rght = array;
   }

   if(P::fieldSolverFiniteDifferencingAtBoundaries && (sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY))) {
      array->at(fsgrids::volfields::dPERBXVOLdy) = (rght->at(fsgrids::volfields::PERBXVOL)-left->at(fsgrids::volfields::PERBXVOL))/2;
      array->at(fsgrids::volfields::dPERBYVOLdy) = (rght->at(fsgrids::volfields::PERBYVOL)-left->at(fsgrids::volfields::PERBYVOL))/2;
      array->at(fsgrids::volfields::dPERBZVOLdy) = (rght->at(fsgrids::volfields::PERBZVOL)-left->at(fsgrids::volfields::PERBZVOL))/2;
   } else {
      array->at(fsgrids::volfields::dPERBXVOLdy) = limiter(left->at(fsgrids::volfields::PERBXVOL),array->at(fsgrids::volfields::PERBXVOL),rght->at(fsgrids::volfields::PERBXVOL));
      array->at(fsgrids::volfields::dPERBYVOLdy) = limiter(left->at(fsgrids::volfields::PERBYVOL),array->at(fsgrids::volfields::PERBYVOL),rght->at(fsgrids::volfields::PERBYVOL));
      array->at(fsgrids::volfields::dPERBZVOLdy) = limiter(left->at(fsgrids::volfields::PERBZVOL),array->at(fsgrids::volfields::PERBZVOL),rght->at(fsgrids::volfields::PERBZVOL));
   }

   // Calculate z-derivatives (is not TVD for AMR mesh):
   left = volGrid.get(i,j,k-1);
   rght = volGrid.get(i,j,k+1);

   if(left == NULL) {
      left = array;
   }
   if(rght == NULL) {
      rght = array;
   }

   if(P::fieldSolverFiniteDifferencingAtBoundaries && (sysBoundaryLayer == 1 || (sysBoundaryLayer == 2 && sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY))) {
      array->at(fsgrids::volfields::dPERBXVOLdz) = (rght->at(fsgrids::volfields::PERBXVOL)-left->at(fsgrids::volfields::PERBXVOL))/2;
      array->at(fsgrids::volfields::dPERBYVOLdz) = (rght->at(fsgrids::volfields::PERBYVOL)-left->at(fsgrids::volfields::PERBYVOL))/2;
      array->at(fsgrids::volfields::dPERBZVOLdz) = (rght->at(fsgrids::volfields::PERBZVOL)-left->at(fsgrids::volfields::PERBZVOL))/2;
   } else {
      array->at(fsgrids::volfields::dPERBXVOLdz) = limiter(left->at(fsgrids::volfields::PERBXVOL),array->at(fsgrids::volfields::PERBXVOL),rght->at(fsgrids::volfields::PERBXVOL));
      array->at(fsgrids::volfields::dPERBYVOLdz) = limiter(left->at(fsgrids::volfields::PERBYVOL),array->at(fsgrids::volfields::PERBYVOL),rght->at(fsgrids::volfields::PERBYVOL));
      array->at(fsgrids::volfields::dPERBZVOLdz) = limiter(left->at(fsgrids::volfields::PERBZVOL),array->at(fsgrids::volfields::PERBZVOL),rght->at(fsgrids::volfields::PERBZVOL));
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
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const FsGridTools::FsIndex_t* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::Timer derivsTimer {"Calculate volume derivatives"};
   int computeTimerId {phiprof::initializeTimer("FS derivatives BVOL compute cells")};

   phiprof::Timer commTimer {"BVOL derivatives ghost updates MPI", {"MPI"}};
   volGrid.updateGhostCells();
   commTimer.stop(N_cells,"Spatial Cells");

   // Calculate derivatives
   #pragma omp parallel
   {
      phiprof::Timer computeTimer {computeTimerId};
      #pragma omp for collapse(2)
      for (FsGridTools::FsIndex_t k=0; k<gridDims[2]; k++) {
         for (FsGridTools::FsIndex_t j=0; j<gridDims[1]; j++) {
            for (FsGridTools::FsIndex_t i=0; i<gridDims[0]; i++) {
               calculateBVOLDerivatives(volGrid,technicalGrid,i,j,k,sysBoundaries);
            }
         }
      }
      computeTimer.stop(N_cells,"Spatial Cells");
   }

   derivsTimer.stop(N_cells,"Spatial Cells");
}

/*! \brief Low-level curvature calculation.
 *
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param bgbGrid fsGrid holding the background fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 *
 * http://fusionwiki.ciemat.es/wiki/Magnetic_curvature
 *
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateCurvature(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & bgbGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   cint i,
   cint j,
   cint k,
   SysBoundary& sysBoundaries
) {
   if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY && technicalGrid.get(i,j,k)->sysBoundaryLayer != 1 && technicalGrid.get(i,j,k)->sysBoundaryLayer != 2) {
      std::array<Real, fsgrids::volfields::N_VOL> * vol = volGrid.get(i,j,k);
      std::array<Real, fsgrids::bgbfield::N_BGB> * bg = bgbGrid.get(i,j,k);

      std::array<Real, fsgrids::volfields::N_VOL> * vol_left_x = volGrid.get(i-1,j,k);
      std::array<Real, fsgrids::volfields::N_VOL> * vol_rght_x = volGrid.get(i+1,j,k);
      std::array<Real, fsgrids::volfields::N_VOL> * vol_left_y = volGrid.get(i,j-1,k);
      std::array<Real, fsgrids::volfields::N_VOL> * vol_rght_y = volGrid.get(i,j+1,k);
      std::array<Real, fsgrids::volfields::N_VOL> * vol_left_z = volGrid.get(i,j,k-1);
      std::array<Real, fsgrids::volfields::N_VOL> * vol_rght_z = volGrid.get(i,j,k+1);
      std::array<Real, fsgrids::bgbfield::N_BGB> * bg_left_x = bgbGrid.get(i-1,j,k);
      std::array<Real, fsgrids::bgbfield::N_BGB> * bg_rght_x = bgbGrid.get(i+1,j,k);
      std::array<Real, fsgrids::bgbfield::N_BGB> * bg_left_y = bgbGrid.get(i,j-1,k);
      std::array<Real, fsgrids::bgbfield::N_BGB> * bg_rght_y = bgbGrid.get(i,j+1,k);
      std::array<Real, fsgrids::bgbfield::N_BGB> * bg_left_z = bgbGrid.get(i,j,k-1);
      std::array<Real, fsgrids::bgbfield::N_BGB> * bg_rght_z = bgbGrid.get(i,j,k+1);

      Real bx = bg->at(fsgrids::bgbfield::BGBXVOL) + vol->at(fsgrids::volfields::PERBXVOL);
      Real by = bg->at(fsgrids::bgbfield::BGBYVOL) + vol->at(fsgrids::volfields::PERBYVOL);
      Real bz = bg->at(fsgrids::bgbfield::BGBZVOL) + vol->at(fsgrids::volfields::PERBZVOL);
      creal bnorm = sqrt(bx*bx + by*by + bz*bz);
      bx /= bnorm;
      by /= bnorm;
      bz /= bnorm;
      Real left_x_bx = bg_left_x->at(fsgrids::bgbfield::BGBXVOL) + vol_left_x->at(fsgrids::volfields::PERBXVOL);
      Real left_x_by = bg_left_x->at(fsgrids::bgbfield::BGBYVOL) + vol_left_x->at(fsgrids::volfields::PERBYVOL);
      Real left_x_bz = bg_left_x->at(fsgrids::bgbfield::BGBZVOL) + vol_left_x->at(fsgrids::volfields::PERBZVOL);
      creal left_x_bnorm = sqrt(left_x_bx*left_x_bx + left_x_by*left_x_by + left_x_bz*left_x_bz);
      left_x_bx /= left_x_bnorm;
      left_x_by /= left_x_bnorm;
      left_x_bz /= left_x_bnorm;

      Real rght_x_bx = bg_rght_x->at(fsgrids::bgbfield::BGBXVOL) + vol_rght_x->at(fsgrids::volfields::PERBXVOL);
      Real rght_x_by = bg_rght_x->at(fsgrids::bgbfield::BGBYVOL) + vol_rght_x->at(fsgrids::volfields::PERBYVOL);
      Real rght_x_bz = bg_rght_x->at(fsgrids::bgbfield::BGBZVOL) + vol_rght_x->at(fsgrids::volfields::PERBZVOL);
      creal rght_x_bnorm = sqrt(rght_x_bx*rght_x_bx + rght_x_by*rght_x_by + rght_x_bz*rght_x_bz);
      rght_x_bx /= rght_x_bnorm;
      rght_x_by /= rght_x_bnorm;
      rght_x_bz /= rght_x_bnorm;

      Real left_y_bx = bg_left_y->at(fsgrids::bgbfield::BGBXVOL) + vol_left_y->at(fsgrids::volfields::PERBXVOL);
      Real left_y_by = bg_left_y->at(fsgrids::bgbfield::BGBYVOL) + vol_left_y->at(fsgrids::volfields::PERBYVOL);
      Real left_y_bz = bg_left_y->at(fsgrids::bgbfield::BGBZVOL) + vol_left_y->at(fsgrids::volfields::PERBZVOL);
      creal left_y_bnorm = sqrt(left_y_bx*left_y_bx + left_y_by*left_y_by + left_y_bz*left_y_bz);
      left_y_bx /= left_y_bnorm;
      left_y_by /= left_y_bnorm;
      left_y_bz /= left_y_bnorm;

      Real rght_y_bx = bg_rght_y->at(fsgrids::bgbfield::BGBXVOL) + vol_rght_y->at(fsgrids::volfields::PERBXVOL);
      Real rght_y_by = bg_rght_y->at(fsgrids::bgbfield::BGBYVOL) + vol_rght_y->at(fsgrids::volfields::PERBYVOL);
      Real rght_y_bz = bg_rght_y->at(fsgrids::bgbfield::BGBZVOL) + vol_rght_y->at(fsgrids::volfields::PERBZVOL);
      creal rght_y_bnorm = sqrt(rght_y_bx*rght_y_bx + rght_y_by*rght_y_by + rght_y_bz*rght_y_bz);
      rght_y_bx /= rght_y_bnorm;
      rght_y_by /= rght_y_bnorm;
      rght_y_bz /= rght_y_bnorm;

      Real left_z_bx = bg_left_z->at(fsgrids::bgbfield::BGBXVOL) + vol_left_z->at(fsgrids::volfields::PERBXVOL);
      Real left_z_by = bg_left_z->at(fsgrids::bgbfield::BGBYVOL) + vol_left_z->at(fsgrids::volfields::PERBYVOL);
      Real left_z_bz = bg_left_z->at(fsgrids::bgbfield::BGBZVOL) + vol_left_z->at(fsgrids::volfields::PERBZVOL);
      creal left_z_bnorm = sqrt(left_z_bx*left_z_bx + left_z_by*left_z_by + left_z_bz*left_z_bz);
      left_z_bx /= left_z_bnorm;
      left_z_by /= left_z_bnorm;
      left_z_bz /= left_z_bnorm;

      Real rght_z_bx = bg_rght_z->at(fsgrids::bgbfield::BGBXVOL) + vol_rght_z->at(fsgrids::volfields::PERBXVOL);
      Real rght_z_by = bg_rght_z->at(fsgrids::bgbfield::BGBYVOL) + vol_rght_z->at(fsgrids::volfields::PERBYVOL);
      Real rght_z_bz = bg_rght_z->at(fsgrids::bgbfield::BGBZVOL) + vol_rght_z->at(fsgrids::volfields::PERBZVOL);
      creal rght_z_bnorm = sqrt(rght_z_bx*rght_z_bx + rght_z_by*rght_z_by + rght_z_bz*rght_z_bz);
      rght_z_bx /= rght_z_bnorm;
      rght_z_by /= rght_z_bnorm;
      rght_z_bz /= rght_z_bnorm;

      vol->at(fsgrids::volfields::CURVATUREX) = bx * 0.5*(rght_x_bx-left_x_bx) / technicalGrid.DX + by * 0.5*(rght_y_bx-left_y_bx) / technicalGrid.DY + bz * 0.5*(rght_z_bx-left_z_bx) / technicalGrid.DZ;
      vol->at(fsgrids::volfields::CURVATUREY) = bx * 0.5*(rght_x_by-left_x_by) / technicalGrid.DX + by * 0.5*(rght_y_by-left_y_by) / technicalGrid.DY + bz * 0.5*(rght_z_by-left_z_by) / technicalGrid.DZ;
      vol->at(fsgrids::volfields::CURVATUREZ) = bx * 0.5*(rght_x_bz-left_x_bz) / technicalGrid.DX + by * 0.5*(rght_y_bz-left_y_bz) / technicalGrid.DY + bz * 0.5*(rght_z_bz-left_z_bz) / technicalGrid.DZ;
   }
}

/*! \brief High-level curvature calculation wrapper function.
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param bgbGrid fsGrid holding the background fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateCurvatureSimple(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & bgbGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries
) {
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const FsGridTools::FsIndex_t* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::Timer curvatureTimer {"Calculate curvature"};
   int computeTimerId {phiprof::initializeTimer("Calculate curvature compute cells")};

   phiprof::Timer commTimer {"Calculate curvature ghost updates MPI", {"MPI"}};
   volGrid.updateGhostCells();
   commTimer.stop(N_cells,"Spatial Cells");

   #pragma omp parallel
   {
      phiprof::Timer computeTimer {computeTimerId};
      #pragma omp for collapse(2)
      for (FsGridTools::FsIndex_t k=0; k<gridDims[2]; k++) {
         for (FsGridTools::FsIndex_t j=0; j<gridDims[1]; j++) {
            for (FsGridTools::FsIndex_t i=0; i<gridDims[0]; i++) {
               if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
                   technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::OUTER_BOUNDARY_PADDING) {
                  continue;
               }
               calculateCurvature(volGrid,bgbGrid,technicalGrid,i,j,k,sysBoundaries);
            }
         }
      }
      computeTimer.stop(N_cells, "Spatial Cells");
   }

   curvatureTimer.stop(N_cells, "Spatial Cells");
}

/*! \brief Returns perturbed volumetric B of cell
 *
 */
[[maybe_unused]] static std::array<Real, 3> getPerBVol(SpatialCell* cell)
{
   return std::array<Real, 3> { {cell->parameters[CellParams::PERBXVOL], cell->parameters[CellParams::PERBYVOL], cell->parameters[CellParams::PERBZVOL]} };
}

/*! \brief Returns volumetric B of cell
 *
 */
static std::array<Real, 3> getBVol(SpatialCell* cell)
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
   return std::array<Real, 3> { 
      {
         rho * cell->parameters[CellParams::VX], 
         rho * cell->parameters[CellParams::VY], 
         rho * cell->parameters[CellParams::VZ]} 
   };
}

/*! \brief Calculates energy density for spatial cell
 *
 */
static Real calculateU(SpatialCell* cell)
{
   Real rho = cell->parameters[CellParams::RHOM];
   std::array<Real, 3> p = getMomentumDensity(cell);
   std::array<Real, 3> B = getBVol(cell);
   return (pow(B[0], 2) + pow(B[1], 2) + pow(B[2], 2)) / (2.0 * physicalconstants::MU_0) + // Magnetic field energy
      (rho > EPS ? (pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2)) / (2.0 * cell->parameters[CellParams::RHOM]) : 0.0); // Kinetic energy
}

/*! \brief Calculates pressure anistotropy from B and P
 *  \param P elements of pressure order in order: P_11, P_22, P_33, P_23, P_13, P_12
 */
static Real calculateAnisotropy(const Eigen::Matrix3d& rot, const std::array<Real, 6>& P)
{
   // Now, rotation matrix to get parallel and perpendicular pressure
   // Eigen::Quaterniond q {Quaterniond::FromTwoVectors(Eigen::vector3d{0, 0, 1}, Eigen::vector3d{myB[0], myB[1], myB[2]})};
   // Eigen::Matrix3d rot = q.toRotationMatrix();
   Eigen::Matrix3d Ptensor {
      {P[0], P[5], P[4]},
      {P[5], P[1], P[3]},
      {P[4], P[3], P[2]},
   };
   
   Eigen::Matrix3d transposerot = rot.transpose();
   Eigen::Matrix3d Pprime = rot * Ptensor * transposerot;

   Real Panisotropy {0.0};
   if (Pprime(2, 2) > EPS) {
      Panisotropy = (Pprime(0, 0) + Pprime(1, 1)) / (2 * Pprime(2, 2));
   }

   return Panisotropy;
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
   Real myU {calculateU(cell)};
   Real myV {std::sqrt(std::pow(cell->parameters[CellParams::VX], 2) + std::pow(cell->parameters[CellParams::VY], 2) + std::pow(cell->parameters[CellParams::VZ], 2))};
   Real maxV {myV};
   std::array<Real, 3> myP = getMomentumDensity(cell);
   std::array<Real, 3> myB = getBVol(cell);
   for (SpatialCell* neighbor : neighbors) {
      Real otherRho = neighbor->parameters[CellParams::RHOM];
      Real otherU = calculateU(neighbor);
      Real otherV {std::sqrt(std::pow(neighbor->parameters[CellParams::VX], 2) + std::pow(neighbor->parameters[CellParams::VY], 2) + std::pow(neighbor->parameters[CellParams::VZ], 2))};
      std::array<Real, 3> otherP = getMomentumDensity(neighbor);
      std::array<Real, 3> otherB = getBVol(neighbor);
      Real deltaBsq = pow(myB[0] - otherB[0], 2) + pow(myB[1] - otherB[1], 2) + pow(myB[2] - otherB[2], 2);

      if (myV < EPS) {
         maxV = std::max(maxV,otherV);
      }
      
      Real maxRho = std::max(myRho, otherRho);
      if (maxRho > EPS) {
         dRho = std::max(fabs(myRho - otherRho) / maxRho, dRho);
      }
      Real maxU = std::max(myU, otherU);
      if (maxU > EPS) {
         dU = std::max(fabs(myU - otherU) / maxU, dU);
         dBsq = std::max(deltaBsq / (2 * physicalconstants::MU_0 * maxU), dBsq);
         if (myRho > EPS) {
            dPsq = std::max((pow(myP[0] - otherP[0], 2) + pow(myP[1] - otherP[1], 2) + pow(myP[2] - otherP[2], 2)) / (2 * myRho * maxU), dPsq);
         }
      }
      Real maxB = sqrt(std::max(pow(myB[0], 2) + pow(myB[1], 2) + pow(myB[2], 2), pow(otherB[0], 2) + pow(otherB[1], 2) + pow(otherB[2], 2)));
      if (maxB > EPS) {
         dB = std::max(sqrt(deltaBsq) / maxB, dB);
      }
   }
   
   Real alpha {0.0};
   alpha = std::max(alpha, dRho * P::alphaDRhoWeight);
   alpha = std::max(alpha, dU * P::alphaDUWeight);
   alpha = std::max(alpha, dPsq * P::alphaDPSqWeight);
   alpha = std::max(alpha, dBsq * P::alphaDBSqWeight);
   alpha = std::max(alpha, dB * P::alphaDBWeight);

   Real dBXdy {cell->derivativesBVOL[bvolderivatives::dPERBXVOLdy]};
   Real dBXdz {cell->derivativesBVOL[bvolderivatives::dPERBXVOLdz]};
   Real dBYdx {cell->derivativesBVOL[bvolderivatives::dPERBYVOLdx]};
   Real dBYdz {cell->derivativesBVOL[bvolderivatives::dPERBYVOLdz]};
   Real dBZdx {cell->derivativesBVOL[bvolderivatives::dPERBZVOLdx]};
   Real dBZdy {cell->derivativesBVOL[bvolderivatives::dPERBZVOLdy]};

   // Note missing factor of mu_0, since we want B and J in same units later
   myB = getBVol(cell); // Redundant, but this makes sure we use total B here
   std::array<Real, 3> myJ = {dBZdy - dBYdz, dBXdz - dBZdx, dBYdx - dBXdy};
   Real BdotJ {0.0};
   Real Bsq {0.0};
   Real J {0.0};
   for (int i = 0; i < 3; ++i) {
      BdotJ += myB[i] * myJ[i];
      Bsq += myB[i] * myB[i];
      J += myJ[i] * myJ[i];
   }
   J = std::sqrt(J);

   Real Bperp {0.0};
   if (Bsq > EPS) {
      for (int i = 0; i < 3; ++i) {
         Bperp += std::pow(myB[i] * (1 - BdotJ / Bsq), 2);
      }
      Bperp = std::sqrt(Bperp);
   }

   // Vorticity
   Real dVxdy {cell->derivativesV[vderivatives::dVxdy]};
   Real dVxdz {cell->derivativesV[vderivatives::dVxdz]};
   Real dVydx {cell->derivativesV[vderivatives::dVydx]};
   Real dVydz {cell->derivativesV[vderivatives::dVydz]};
   Real dVzdx {cell->derivativesV[vderivatives::dVzdx]};
   Real dVzdy {cell->derivativesV[vderivatives::dVzdy]};
   Real vorticity {std::sqrt(std::pow(dVzdy - dVydz, 2) + std::pow(dVxdz - dVzdx, 2 ) + std::pow(dVydx - dVxdy, 2))};
   //Real vA {std::sqrt(Bsq / (physicalconstants::MU_0 * myRho))};
   Real amr_vorticity {-1.0}; // Error value
   if (maxV > EPS) {
      amr_vorticity = vorticity * cell->parameters[CellParams::DX] / maxV;
   }

   std::array<Real, 6> myPressure {cell->parameters[CellParams::P_11], cell->parameters[CellParams::P_22], cell->parameters[CellParams::P_33], cell->parameters[CellParams::P_23], cell->parameters[CellParams::P_13], cell->parameters[CellParams::P_12]};
   Eigen::Matrix3d rot = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d{myB[0], myB[1], myB[2]}, Eigen::Vector3d{0, 0, 1}).normalized().toRotationMatrix();
   Real Panisotropy {calculateAnisotropy(rot, myPressure)};
   for (const auto& pop : cell->get_populations()) {
      // TODO I hate this. Change all this crap to std::vectors?
      std::array<Real, 6> popP {pop.P[0], pop.P[1], pop.P[2], pop.P[3], pop.P[4], pop.P[5]};
      Real popPanisotropy {calculateAnisotropy(rot, popP)};
      // low value refines
      Panisotropy = std::min(Panisotropy, popPanisotropy);
   }

   cell->parameters[CellParams::AMR_DRHO] = dRho;
   cell->parameters[CellParams::AMR_DU] = dU;
   cell->parameters[CellParams::AMR_DPSQ] = dPsq;
   cell->parameters[CellParams::AMR_DBSQ] = dBsq;
   cell->parameters[CellParams::AMR_DB] = dB;
   cell->parameters[CellParams::AMR_ALPHA1] = alpha;
   cell->parameters[CellParams::AMR_ALPHA2] = cell->parameters[CellParams::DX] * J / (Bperp + EPS);   // Epsilon in denominator so we don't get infinities
   cell->parameters[CellParams::P_ANISOTROPY] = Panisotropy;
   // Experimental, current scaling is bulk velocity
   cell->parameters[CellParams::AMR_VORTICITY] = amr_vorticity;
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
   phiprof::Timer gradientsTimer {"Calculate volume gradients"};
   int computeTimerId {phiprof::initializeTimer("Calculate volume gradients compute cells")};

   phiprof::Timer commTimer {"Calculate volume gradients ghost updates MPI", {"MPI"}};
   // We only need nearest neighbourhood and spatial data here
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);
   mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::NEAREST);
   commTimer.stop(N_cells,"Spatial Cells");

   // Calculate derivatives
   #pragma omp parallel
   {
      phiprof::Timer computeTimer {computeTimerId};
      #pragma omp for
      for (uint i = 0; i < cells.size(); ++i) {
         CellID id = cells[i];
         SpatialCell* cell = mpiGrid[id];
         std::vector<SpatialCell*> neighbors;
         for (const auto& [neighbor, dir] : mpiGrid.get_face_neighbors_of(id)) {
            neighbors.push_back(mpiGrid[neighbor]);
         }
         calculateScaledDeltas(cell, neighbors);
      }
      computeTimer.stop(N_cells,"Spatial Cells");
   }

   

   gradientsTimer.stop(N_cells,"Spatial Cells");
}
