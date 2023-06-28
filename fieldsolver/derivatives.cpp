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
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa calculateDerivativesSimple calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivatives(
   cint i,
   cint j,
   cint k,
   FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid<Real, fsgrids::moments::N_MOMENTS, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid<Real, fsgrids::dperb::N_DPERB, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid<Real, fsgrids::dmoments::N_DMOMENTS, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {
   Real* dPerB = dPerBGrid.get(i,j,k);
   Real* dMoments = dMomentsGrid.get(i,j,k);

   // Get boundary flag for the cell:
   cuint sysBoundaryFlag  = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   cuint sysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;

   // Constants for electron pressure derivatives
   // Upstream pressure
   Real Peupstream = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   Real Peconst = Peupstream * pow(Parameters::electronDensity, -Parameters::electronPTindex);

   Real* leftMoments = NULL;
   Real* leftPerB = NULL;
   Real* centMoments = momentsGrid.get(i,j,k);
   Real* centPerB = perBGrid.get(i,j,k);
   #ifdef DEBUG_SOLVERS
   if (centMoments[fsgrids::moments::RHOM] <= 0) {
      std::cerr << __FILE__ << ":" << __LINE__
         << (centMoments[fsgrids::moments::RHOM] < 0 ? " Negative" : " Zero") << " density in spatial cell at (" << i << " " << j << " " << k << ")"
         << std::endl;
      abort();
   }
   #endif
   Real* rghtMoments = NULL;
   Real* rghtPerB = NULL;
   Real* botLeft = NULL;
   Real* botRght = NULL;
   Real* topLeft = NULL;
   Real* topRght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      
      leftPerB = perBGrid.get(i-1,j,k);
      rghtPerB = perBGrid.get(i+1,j,k);
      leftMoments = momentsGrid.get(i-1,j,k);
      rghtMoments = momentsGrid.get(i+1,j,k);
      #ifdef DEBUG_SOLVERS
      if (leftMoments[fsgrids::moments::RHOM] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (leftMoments[fsgrids::moments::RHOM] < 0 ? " Negative" : " Zero") << " density in spatial cell " //<< leftNbrID
            << std::endl;
         abort();
      }
      if (rghtMoments[fsgrids::moments::RHOM] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (rghtMoments[fsgrids::moments::RHOM] < 0 ? " Negative" : " Zero") << " density in spatial cell " //<< rightNbrID
            << std::endl;
         abort();
      }
      #endif
      
      dMoments[fsgrids::dmoments::drhomdx] = limiter(leftMoments[fsgrids::moments::RHOM],centMoments[fsgrids::moments::RHOM],rghtMoments[fsgrids::moments::RHOM]);
      dMoments[fsgrids::dmoments::drhoqdx] = limiter(leftMoments[fsgrids::moments::RHOQ],centMoments[fsgrids::moments::RHOQ],rghtMoments[fsgrids::moments::RHOQ]);
      dMoments[fsgrids::dmoments::dp11dx] = limiter(leftMoments[fsgrids::moments::P_11],centMoments[fsgrids::moments::P_11],rghtMoments[fsgrids::moments::P_11]);
      dMoments[fsgrids::dmoments::dp22dx] = limiter(leftMoments[fsgrids::moments::P_22],centMoments[fsgrids::moments::P_22],rghtMoments[fsgrids::moments::P_22]);
      dMoments[fsgrids::dmoments::dp33dx] = limiter(leftMoments[fsgrids::moments::P_33],centMoments[fsgrids::moments::P_33],rghtMoments[fsgrids::moments::P_33]);

      dMoments[fsgrids::dmoments::dVxdx]  = limiter(leftMoments[fsgrids::moments::VX], centMoments[fsgrids::moments::VX], rghtMoments[fsgrids::moments::VX]);
      dMoments[fsgrids::dmoments::dVydx]  = limiter(leftMoments[fsgrids::moments::VY], centMoments[fsgrids::moments::VY], rghtMoments[fsgrids::moments::VY]);
      dMoments[fsgrids::dmoments::dVzdx]  = limiter(leftMoments[fsgrids::moments::VZ], centMoments[fsgrids::moments::VZ], rghtMoments[fsgrids::moments::VZ]);
      dPerB[fsgrids::dperb::dPERBydx]  = limiter(leftPerB[fsgrids::bfield::PERBY],centPerB[fsgrids::bfield::PERBY],rghtPerB[fsgrids::bfield::PERBY]);
      dPerB[fsgrids::dperb::dPERBzdx]  = limiter(leftPerB[fsgrids::bfield::PERBZ],centPerB[fsgrids::bfield::PERBZ],rghtPerB[fsgrids::bfield::PERBZ]);

      // pres_e = const * np.power(rho_e, index)
      dMoments[fsgrids::dmoments::dPedx] = Peconst * limiter(pow(leftMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex));      
      
      if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
        dPerB[fsgrids::dperb::dPERBydxx] = 0.0;
        dPerB[fsgrids::dperb::dPERBzdxx] = 0.0;
      } else {
        dPerB[fsgrids::dperb::dPERBydxx] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] - 2.0*centPerB[fsgrids::bfield::PERBY];
        dPerB[fsgrids::dperb::dPERBzdxx] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] - 2.0*centPerB[fsgrids::bfield::PERBZ];
      }
   } else {
      // Boundary conditions handle derivatives.
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 0);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 0);
      }
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      
      leftPerB = perBGrid.get(i,j-1,k);
      rghtPerB = perBGrid.get(i,j+1,k);
      leftMoments = momentsGrid.get(i,j-1,k);
      rghtMoments = momentsGrid.get(i,j+1,k);
      
       dMoments[fsgrids::dmoments::drhomdy] = limiter(leftMoments[fsgrids::moments::RHOM],centMoments[fsgrids::moments::RHOM],rghtMoments[fsgrids::moments::RHOM]);
       dMoments[fsgrids::dmoments::drhoqdy] = limiter(leftMoments[fsgrids::moments::RHOQ],centMoments[fsgrids::moments::RHOQ],rghtMoments[fsgrids::moments::RHOQ]);
       dMoments[fsgrids::dmoments::dp11dy] = limiter(leftMoments[fsgrids::moments::P_11],centMoments[fsgrids::moments::P_11],rghtMoments[fsgrids::moments::P_11]);
       dMoments[fsgrids::dmoments::dp22dy] = limiter(leftMoments[fsgrids::moments::P_22],centMoments[fsgrids::moments::P_22],rghtMoments[fsgrids::moments::P_22]);
       dMoments[fsgrids::dmoments::dp33dy] = limiter(leftMoments[fsgrids::moments::P_33],centMoments[fsgrids::moments::P_33],rghtMoments[fsgrids::moments::P_33]);
       dMoments[fsgrids::dmoments::dVxdy]  = limiter(leftMoments[fsgrids::moments::VX], centMoments[fsgrids::moments::VX], rghtMoments[fsgrids::moments::VX]);
       dMoments[fsgrids::dmoments::dVydy]  = limiter(leftMoments[fsgrids::moments::VY], centMoments[fsgrids::moments::VY], rghtMoments[fsgrids::moments::VY]);
       dMoments[fsgrids::dmoments::dVzdy]  = limiter(leftMoments[fsgrids::moments::VZ], centMoments[fsgrids::moments::VZ], rghtMoments[fsgrids::moments::VZ]);

      dPerB[fsgrids::dperb::dPERBxdy]  = limiter(leftPerB[fsgrids::bfield::PERBX],centPerB[fsgrids::bfield::PERBX],rghtPerB[fsgrids::bfield::PERBX]);
      dPerB[fsgrids::dperb::dPERBzdy]  = limiter(leftPerB[fsgrids::bfield::PERBZ],centPerB[fsgrids::bfield::PERBZ],rghtPerB[fsgrids::bfield::PERBZ]);

      // pres_e = const * np.power(rho_e, index)
      dMoments[fsgrids::dmoments::dPedy] = Peconst * limiter(pow(leftMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex));      

      if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
         dPerB[fsgrids::dperb::dPERBxdyy] = 0.0;
         dPerB[fsgrids::dperb::dPERBzdyy] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBxdyy] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] - 2.0*centPerB[fsgrids::bfield::PERBX];
         dPerB[fsgrids::dperb::dPERBzdyy] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] - 2.0*centPerB[fsgrids::bfield::PERBZ];
      }
      
   } else {
      // Boundary conditions handle derivatives.
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 1);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      
      leftPerB = perBGrid.get(i,j,k-1);
      rghtPerB = perBGrid.get(i,j,k+1);
      leftMoments = momentsGrid.get(i,j,k-1);
      rghtMoments = momentsGrid.get(i,j,k+1);
      
      dMoments[fsgrids::dmoments::drhomdz] = limiter(leftMoments[fsgrids::moments::RHOM],centMoments[fsgrids::moments::RHOM],rghtMoments[fsgrids::moments::RHOM]);
      dMoments[fsgrids::dmoments::drhoqdz] = limiter(leftMoments[fsgrids::moments::RHOQ],centMoments[fsgrids::moments::RHOQ],rghtMoments[fsgrids::moments::RHOQ]);
      dMoments[fsgrids::dmoments::dp11dz] = limiter(leftMoments[fsgrids::moments::P_11],centMoments[fsgrids::moments::P_11],rghtMoments[fsgrids::moments::P_11]);
      dMoments[fsgrids::dmoments::dp22dz] = limiter(leftMoments[fsgrids::moments::P_22],centMoments[fsgrids::moments::P_22],rghtMoments[fsgrids::moments::P_22]);
      dMoments[fsgrids::dmoments::dp33dz] = limiter(leftMoments[fsgrids::moments::P_33],centMoments[fsgrids::moments::P_33],rghtMoments[fsgrids::moments::P_33]);
      dMoments[fsgrids::dmoments::dVxdz]  = limiter(leftMoments[fsgrids::moments::VX], centMoments[fsgrids::moments::VX], rghtMoments[fsgrids::moments::VX]);
      dMoments[fsgrids::dmoments::dVydz]  = limiter(leftMoments[fsgrids::moments::VY], centMoments[fsgrids::moments::VY], rghtMoments[fsgrids::moments::VY]);
      dMoments[fsgrids::dmoments::dVzdz]  = limiter(leftMoments[fsgrids::moments::VZ], centMoments[fsgrids::moments::VZ], rghtMoments[fsgrids::moments::VZ]);
      
      dPerB[fsgrids::dperb::dPERBxdz]  = limiter(leftPerB[fsgrids::bfield::PERBX],centPerB[fsgrids::bfield::PERBX],rghtPerB[fsgrids::bfield::PERBX]);
      dPerB[fsgrids::dperb::dPERBydz]  = limiter(leftPerB[fsgrids::bfield::PERBY],centPerB[fsgrids::bfield::PERBY],rghtPerB[fsgrids::bfield::PERBY]);

      // pres_e = const * np.power(rho_e, index)
      dMoments[fsgrids::dmoments::dPedz] = Peconst * limiter(pow(leftMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex));      

      if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
        dPerB[fsgrids::dperb::dPERBxdzz] = 0.0;
        dPerB[fsgrids::dperb::dPERBydzz] = 0.0;
      } else {
        dPerB[fsgrids::dperb::dPERBxdzz] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] - 2.0*centPerB[fsgrids::bfield::PERBX];
        dPerB[fsgrids::dperb::dPERBydzz] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] - 2.0*centPerB[fsgrids::bfield::PERBY];
      }
      
   } else {
      // Boundary conditions handle derivatives.
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 2);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 2);
      }
   }
   
   if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
      dPerB[fsgrids::dperb::dPERBxdyz] = 0.0;
      dPerB[fsgrids::dperb::dPERBydxz] = 0.0;
      dPerB[fsgrids::dperb::dPERBzdxy] = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
         
         botLeft = perBGrid.get(i-1,j-1,k);
         botRght = perBGrid.get(i+1,j-1,k);
         topLeft = perBGrid.get(i-1,j+1,k);
         topRght = perBGrid.get(i+1,j+1,k);
         
         dPerB[fsgrids::dperb::dPERBzdxy] = FOURTH * (botLeft[fsgrids::bfield::PERBZ] + topRght[fsgrids::bfield::PERBZ] - botRght[fsgrids::bfield::PERBZ] - topLeft[fsgrids::bfield::PERBZ]);
         
      } else {
         // Boundary conditions handle derivatives.
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 3);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 3);
         }
      }
      
      // Calculate xz mixed derivatives:
      if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
         
         botLeft = perBGrid.get(i-1,j,k-1);
         botRght = perBGrid.get(i+1,j,k-1);
         topLeft = perBGrid.get(i-1,j,k+1);
         topRght = perBGrid.get(i+1,j,k+1);
         
         dPerB[fsgrids::dperb::dPERBydxz] = FOURTH * (botLeft[fsgrids::bfield::PERBY] + topRght[fsgrids::bfield::PERBY] - botRght[fsgrids::bfield::PERBY] - topLeft[fsgrids::bfield::PERBY]);
         
      } else {
         // Boundary conditions handle derivatives.
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 4);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 4);
         }
      }
      
      // Calculate yz mixed derivatives:
      if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
         
         botLeft = perBGrid.get(i,j-1,k-1);
         botRght = perBGrid.get(i,j+1,k-1);
         topLeft = perBGrid.get(i,j-1,k+1);
         topRght = perBGrid.get(i,j+1,k+1);
         
         dPerB[fsgrids::dperb::dPERBxdyz] = FOURTH * (botLeft[fsgrids::bfield::PERBX] + topRght[fsgrids::bfield::PERBX] - botRght[fsgrids::bfield::PERBX] - topLeft[fsgrids::bfield::PERBX]);
         
      } else {
         // Boundary conditions handle derivatives.
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 5);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 5);
         }
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
   FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> & perBDt2Grid,
   FsGrid<Real, fsgrids::moments::N_MOMENTS, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid<Real, fsgrids::moments::N_MOMENTS, FS_STENCIL_WIDTH> & momentsDt2Grid,
   FsGrid<Real, fsgrids::dperb::N_DPERB, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid<Real, fsgrids::dmoments::N_DMOMENTS, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
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
 * For the cell with ID cellID calculate the spatial derivatives of BVOL or apply the derivative boundary conditions defined in project.h.
 * 
 * \param volGrid fsGrid holding the volume averaged fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateBVOLDerivatives(
   FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH> & volGrid,
   FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
   cint i,
   cint j,
   cint k,
   SysBoundary& sysBoundaries
) {
   Real* array = volGrid.get(i,j,k);
   Real* left = NULL;
   Real* rght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      left = volGrid.get(i-1,j,k);
      rght = volGrid.get(i+1,j,k);
      
      array[fsgrids::volfields::dPERBYVOLdx] = limiter(left[fsgrids::volfields::PERBYVOL],array[fsgrids::volfields::PERBYVOL],rght[fsgrids::volfields::PERBYVOL]);
      array[fsgrids::volfields::dPERBZVOLdx] = limiter(left[fsgrids::volfields::PERBZVOL],array[fsgrids::volfields::PERBZVOL],rght[fsgrids::volfields::PERBZVOL]);
   } else {
      if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 0);
      } else {
         sysBoundaries.getSysBoundary(technicalGrid.get(i,j,k)->sysBoundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 0);
      }
   }
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
   if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      left = volGrid.get(i,j-1,k);
      rght = volGrid.get(i,j+1,k);
      
      array[fsgrids::volfields::dPERBXVOLdy] = limiter(left[fsgrids::volfields::PERBXVOL],array[fsgrids::volfields::PERBXVOL],rght[fsgrids::volfields::PERBXVOL]);
      array[fsgrids::volfields::dPERBZVOLdy] = limiter(left[fsgrids::volfields::PERBZVOL],array[fsgrids::volfields::PERBZVOL],rght[fsgrids::volfields::PERBZVOL]);
   } else {
      if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 1);
      } else {
         sysBoundaries.getSysBoundary(technicalGrid.get(i,j,k)->sysBoundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      left = volGrid.get(i,j,k-1);
      rght = volGrid.get(i,j,k+1);
      
      array[fsgrids::volfields::dPERBXVOLdz] = limiter(left[fsgrids::volfields::PERBXVOL],array[fsgrids::volfields::PERBXVOL],rght[fsgrids::volfields::PERBXVOL]);
      array[fsgrids::volfields::dPERBYVOLdz] = limiter(left[fsgrids::volfields::PERBYVOL],array[fsgrids::volfields::PERBYVOL],rght[fsgrids::volfields::PERBYVOL]);
   } else {
      if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 2);
      } else {
         sysBoundaries.getSysBoundary(technicalGrid.get(i,j,k)->sysBoundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 2);
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
 * \param sysBoundaries System boundary conditions existing
 * 
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateBVOLDerivativesSimple(
   FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH> & volGrid,
   FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
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
            if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
            
            calculateBVOLDerivatives(volGrid,technicalGrid,i,j,k,sysBoundaries);
         }
      }
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");

   phiprof::stop("Calculate volume derivatives",N_cells,"Spatial Cells");
}
