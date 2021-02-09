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
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in project.h. Uses RHO, V[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the second-order method, and RHO_DT2, V[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param boundaries Boundary conditions existing
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
   Boundary& boundaries,
   cint& RKCase
) {
   std::array<Real, fsgrids::dperb::N_DPERB> * dPerB = dPerBGrid.get(i,j,k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS> * dMoments = dMomentsGrid.get(i,j,k);

   // Get boundary flag for the cell:
   cuint boundaryFlag  = technicalGrid.get(i,j,k)->boundaryFlag;
   cuint boundaryLayer = technicalGrid.get(i,j,k)->boundaryLayer;

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
   if ((boundaryFlag == boundarytype::NOT_BOUNDARY) || (boundaryLayer == 1)) {
      
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

      // pres_e = const * np.power(rho_e, index)
      dMoments->at(fsgrids::dmoments::dPedx) = Peconst * limiter(pow(leftMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex));      
      
      if (Parameters::ohmHallTerm < 2 || boundaryLayer == 1) {
        dPerB->at(fsgrids::dperb::dPERBydxx) = 0.0;
        dPerB->at(fsgrids::dperb::dPERBzdxx) = 0.0;
      } else {
        dPerB->at(fsgrids::dperb::dPERBydxx) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
        dPerB->at(fsgrids::dperb::dPERBzdxx) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
      }
   } else {
      // Boundary conditions handle derivatives.
      if (boundaryFlag == boundarytype::NOT_BOUNDARY) {
         BC::BoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 0);
      } else {
         boundaries.getBoundary(boundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 0);
      }
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
   if ((boundaryFlag == boundarytype::NOT_BOUNDARY) || (boundaryLayer == 1)) {
      
      leftPerB = perBGrid.get(i,j-1,k);
      rghtPerB = perBGrid.get(i,j+1,k);
      leftMoments = momentsGrid.get(i,j-1,k);
      rghtMoments = momentsGrid.get(i,j+1,k);
      
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

      // pres_e = const * np.power(rho_e, index)
      dMoments->at(fsgrids::dmoments::dPedy) = Peconst * limiter(pow(leftMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex));      

      if (Parameters::ohmHallTerm < 2 || boundaryLayer == 1) {
         dPerB->at(fsgrids::dperb::dPERBxdyy) = 0.0;
         dPerB->at(fsgrids::dperb::dPERBzdyy) = 0.0;
      } else {
         dPerB->at(fsgrids::dperb::dPERBxdyy) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
         dPerB->at(fsgrids::dperb::dPERBzdyy) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
      }
      
   } else {
      // Boundary conditions handle derivatives.
      if (boundaryFlag == boundarytype::NOT_BOUNDARY) {
         BC::BoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 1);
      } else {
         boundaries.getBoundary(boundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if ((boundaryFlag == boundarytype::NOT_BOUNDARY) || (boundaryLayer == 1)) {
      
      leftPerB = perBGrid.get(i,j,k-1);
      rghtPerB = perBGrid.get(i,j,k+1);
      leftMoments = momentsGrid.get(i,j,k-1);
      rghtMoments = momentsGrid.get(i,j,k+1);
      
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

      // pres_e = const * np.power(rho_e, index)
      dMoments->at(fsgrids::dmoments::dPedz) = Peconst * limiter(pow(leftMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments->at(fsgrids::moments::RHOQ)/physicalconstants::CHARGE,Parameters::electronPTindex));      

      if (Parameters::ohmHallTerm < 2 || boundaryLayer == 1) {
        dPerB->at(fsgrids::dperb::dPERBxdzz) = 0.0;
        dPerB->at(fsgrids::dperb::dPERBydzz) = 0.0;
      } else {
        dPerB->at(fsgrids::dperb::dPERBxdzz) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
        dPerB->at(fsgrids::dperb::dPERBydzz) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
      }
      
   } else {
      // Boundary conditions handle derivatives.
      if (boundaryFlag == boundarytype::NOT_BOUNDARY) {
         BC::BoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 2);
      } else {
         boundaries.getBoundary(boundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 2);
      }
   }
   
   if (Parameters::ohmHallTerm < 2 || boundaryLayer == 1) {
      dPerB->at(fsgrids::dperb::dPERBxdyz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBydxz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBzdxy) = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      if ((boundaryFlag == boundarytype::NOT_BOUNDARY) || (boundaryLayer == 1)) {
         
         botLeft = perBGrid.get(i-1,j-1,k);
         botRght = perBGrid.get(i+1,j-1,k);
         topLeft = perBGrid.get(i-1,j+1,k);
         topRght = perBGrid.get(i+1,j+1,k);
         
         dPerB->at(fsgrids::dperb::dPERBzdxy) = FOURTH * (botLeft->at(fsgrids::bfield::PERBZ) + topRght->at(fsgrids::bfield::PERBZ) - botRght->at(fsgrids::bfield::PERBZ) - topLeft->at(fsgrids::bfield::PERBZ));
         
      } else {
         // Boundary conditions handle derivatives.
         if (boundaryFlag == boundarytype::NOT_BOUNDARY) {
            BC::BoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 3);
         } else {
            boundaries.getBoundary(boundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 3);
         }
      }
      
      // Calculate xz mixed derivatives:
      if ((boundaryFlag == boundarytype::NOT_BOUNDARY) || (boundaryLayer == 1)) {
         
         botLeft = perBGrid.get(i-1,j,k-1);
         botRght = perBGrid.get(i+1,j,k-1);
         topLeft = perBGrid.get(i-1,j,k+1);
         topRght = perBGrid.get(i+1,j,k+1);
         
         dPerB->at(fsgrids::dperb::dPERBydxz) = FOURTH * (botLeft->at(fsgrids::bfield::PERBY) + topRght->at(fsgrids::bfield::PERBY) - botRght->at(fsgrids::bfield::PERBY) - topLeft->at(fsgrids::bfield::PERBY));
         
      } else {
         // Boundary conditions handle derivatives.
         if (boundaryFlag == boundarytype::NOT_BOUNDARY) {
            BC::BoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 4);
         } else {
            boundaries.getBoundary(boundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 4);
         }
      }
      
      // Calculate yz mixed derivatives:
      if ((boundaryFlag == boundarytype::NOT_BOUNDARY) || (boundaryLayer == 1)) {
         
         botLeft = perBGrid.get(i,j-1,k-1);
         botRght = perBGrid.get(i,j+1,k-1);
         topLeft = perBGrid.get(i,j-1,k+1);
         topRght = perBGrid.get(i,j+1,k+1);
         
         dPerB->at(fsgrids::dperb::dPERBxdyz) = FOURTH * (botLeft->at(fsgrids::bfield::PERBX) + topRght->at(fsgrids::bfield::PERBX) - botRght->at(fsgrids::bfield::PERBX) - topLeft->at(fsgrids::bfield::PERBX));
         
      } else {
         // Boundary conditions handle derivatives.
         if (boundaryFlag == boundarytype::NOT_BOUNDARY) {
            BC::BoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 5);
         } else {
            boundaries.getBoundary(boundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 5);
         }
      }
   }
}


/*! \brief High-level derivative calculation wrapper function.
 * 

 * B has to be updated because after the boundary update in propagateMagneticFieldSimple there is no consistent state of B yet everywhere.
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
 * \param boundaries Boundary conditions existing
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
   Boundary& boundaries,
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
      // The update of PERB[XYZ] is needed after the
      // boundary update of propagateMagneticFieldSimple.
       perBGrid.updateGhostCells();
       if(communicateMoments) {
         momentsGrid.updateGhostCells();
       }
       break;
    case RK_ORDER2_STEP1:
      // Exchange PERB*_DT2,RHO_DT2,V*_DT2 with neighbours The
      // update of PERB[XYZ]_DT2 is needed after the
      // boundary update of propagateMagneticFieldSimple.
       perBDt2Grid.updateGhostCells();
       if(communicateMoments) {
         momentsDt2Grid.updateGhostCells();
       }
       break;
    case RK_ORDER2_STEP2:
      // Exchange PERB*,RHO,V* with neighbours The update of B
      // is needed after the boundary update of
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
            if (technicalGrid.get(i,j,k)->boundaryFlag == boundarytype::NOTHING) continue;
            if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
               calculateDerivatives(i,j,k, perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, technicalGrid, boundaries, RKCase);
            } else {
               calculateDerivatives(i,j,k, perBDt2Grid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, boundaries, RKCase);
            }
         }
      }
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   phiprof::stop("Calculate face derivatives",N_cells,"Spatial Cells");   
}

/*! \brief Low-level spatial derivatives calculation.
 * 
 * For the cell with ID cellID calculate the spatial derivatives of BVOL or apply the derivative boundary conditions defined in project.h.
 * 
 * \param volGrid fsGrid holding the volume averaged fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param boundaries Boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateBVOLDerivatives(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   cint i,
   cint j,
   cint k,
   Boundary& boundaries
) {
   std::array<Real, fsgrids::volfields::N_VOL> * array = volGrid.get(i,j,k);
   
   std::array<Real, fsgrids::volfields::N_VOL> * left = NULL;
   std::array<Real, fsgrids::volfields::N_VOL> * rght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if (technicalGrid.get(i,j,k)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
      left = volGrid.get(i-1,j,k);
      rght = volGrid.get(i+1,j,k);
      
      array->at(fsgrids::volfields::dPERBYVOLdx) = limiter(left->at(fsgrids::volfields::PERBYVOL),array->at(fsgrids::volfields::PERBYVOL),rght->at(fsgrids::volfields::PERBYVOL));
      array->at(fsgrids::volfields::dPERBZVOLdx) = limiter(left->at(fsgrids::volfields::PERBZVOL),array->at(fsgrids::volfields::PERBZVOL),rght->at(fsgrids::volfields::PERBZVOL));
   } else {
      if (technicalGrid.get(i,j,k)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
         BC::BoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 0);
      } else {
         boundaries.getBoundary(technicalGrid.get(i,j,k)->boundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 0);
      }
   }
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
   if (technicalGrid.get(i,j,k)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
      left = volGrid.get(i,j-1,k);
      rght = volGrid.get(i,j+1,k);
      
      array->at(fsgrids::volfields::dPERBXVOLdy) = limiter(left->at(fsgrids::volfields::PERBXVOL),array->at(fsgrids::volfields::PERBXVOL),rght->at(fsgrids::volfields::PERBXVOL));
      array->at(fsgrids::volfields::dPERBZVOLdy) = limiter(left->at(fsgrids::volfields::PERBZVOL),array->at(fsgrids::volfields::PERBZVOL),rght->at(fsgrids::volfields::PERBZVOL));
   } else {
      if (technicalGrid.get(i,j,k)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
         BC::BoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 1);
      } else {
         boundaries.getBoundary(technicalGrid.get(i,j,k)->boundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if (technicalGrid.get(i,j,k)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
      left = volGrid.get(i,j,k-1);
      rght = volGrid.get(i,j,k+1);
      
      array->at(fsgrids::volfields::dPERBXVOLdz) = limiter(left->at(fsgrids::volfields::PERBXVOL),array->at(fsgrids::volfields::PERBXVOL),rght->at(fsgrids::volfields::PERBXVOL));
      array->at(fsgrids::volfields::dPERBYVOLdz) = limiter(left->at(fsgrids::volfields::PERBYVOL),array->at(fsgrids::volfields::PERBYVOL),rght->at(fsgrids::volfields::PERBYVOL));
   } else {
      if (technicalGrid.get(i,j,k)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
         BC::BoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 2);
      } else {
         boundaries.getBoundary(technicalGrid.get(i,j,k)->boundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 2);
      }
   }
}

/*! \brief High-level derivative calculation wrapper function.
 * 
 * BVOL has been calculated locally by calculateVolumeAveragedFields but not communicated.
 * For the acceleration step one needs the cross-derivatives of BVOL
 * 
 * \param volGrid fsGrid holding the volume averaged fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param boundaries Boundary conditions existing
 * 
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateBVOLDerivativesSimple(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   Boundary& boundaries
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
            if (technicalGrid.get(i,j,k)->boundaryFlag == boundarytype::NOTHING) continue;
            
            calculateBVOLDerivatives(volGrid,technicalGrid,i,j,k,boundaries);
         }
      }
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");

   phiprof::stop("Calculate volume derivatives",N_cells,"Spatial Cells");
}
