/*
This file is part of Vlasiator.

Copyright 2015 Finnish Meteorological Institute
*/

#include <cstdlib>

#include "fs_common.h"
#include "derivatives.hpp"
#include "fs_limiters.h"

/*! \brief Low-level spatial derivatives calculation.
 * 
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in project.h. Uses RHO, RHOV[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the second-order method, and RHO_DT2, RHOV[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 * \param mpiGrid Grid
 * \param cellCache Field solver cell cache
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param doMoments If true, the derivatives of moments (rho, V, P) are computed.
 * 
 * \sa calculateDerivativesSimple calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivatives(
   const int i,
   const int j,
   const int k,
   const FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   const FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   const FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   const FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase,
   const bool& doMoments
) {
   std::array<Real, fsgrids::dperb::N_DPERB> * const dPerB = & dPerBGrid->get(i,j,k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS> * const dMoments = & dMomentsGrid->get(i,j,k);

   // Get boundary flag for the cell:
   cuint sysBoundaryFlag  = technicalGrid->get(i,j,k).sysBoundaryFlag;
   cuint sysBoundaryLayer = technicalGrid->get(i,j,k).sysBoundaryLayer;
   
   std::array<Real, fsgrids::moments::N_MOMENTS> * const leftMoments = NULL;
   std::array<Real, fsgrids::moments::N_BFIELD> * const leftPerB = NULL;
   const std::array<Real, fsgrids::moments::N_MOMENTS> * const centMoments = & momentsGrid->get(i,j,k);
   const std::array<Real, fsgrids::moments::N_BFIELD> * const centPerB = & perBGrid->get(i,j,k);
   #ifdef DEBUG_SOLVERS
   if (centMoments[fsgrids::moments::RHO] <= 0) {
      std::cerr << __FILE__ << ":" << __LINE__
         << (centMoments[fsgrids::moments::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << cellID
         << std::endl;
      abort();
   }
   #endif
   std::array<Real, fsgrids::moments::N_MOMENTS> * const rghtMoments = NULL;
   std::array<Real, fsgrids::moments::N_BFIELD>  * const rghtPerB = NULL;
   std::array<Real, fsgrids::moments::N_BFIELD>  * const botLeft = NULL;
   std::array<Real, fsgrids::moments::N_BFIELD>  * const botRght = NULL;
   std::array<Real, fsgrids::moments::N_BFIELD>  * const topLeft = NULL;
   std::array<Real, fsgrids::moments::N_BFIELD>  * const topRght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         leftPerB = & perBGrid->get(i-1,j,k);
         rghtPerB = & perBGrid->get(i+1,j,k);
         if (doMoments) {
            leftMoments = & momentsGrid->get(i-1,j,k);
            rghtMoments = & momentsGrid->get(i+1,j,k);
         }
      }
      if (RKCase == RK_ORDER2_STEP1) {
         leftPerB = & perBGridDt2->get(i-1,j,k);
         rghtPerB = & perBGridDt2->get(i+1,j,k);
         if (doMoments) {
            leftMoments = & momentsGridDt2->get(i-1,j,k);
            rghtMoments = & momentsGridDt2->get(i+1,j,k);
         }
      }
      
      if (doMoments) {
         dMoments[fsgrids::dmoments::drhodx] = limiter(leftMoments[fsgrids::moments::RHO],centMoments[fsgrids::moments::RHO],rghtMoments[fsgrids::moments::RHO]);
         dMoments[fsgrids::dmoments::dp11dx] = limiter(leftMoments[fsgrids::moments::P_11],centMoments[fsgrids::moments::P_11],rghtMoments[fsgrids::moments::P_11]);
         dMoments[fsgrids::dmoments::dp22dx] = limiter(leftMoments[fsgrids::moments::P_22],centMoments[fsgrids::moments::P_22],rghtMoments[fsgrids::moments::P_22]);
         dMoments[fsgrids::dmoments::dp33dx] = limiter(leftMoments[fsgrids::moments::P_33],centMoments[fsgrids::moments::P_33],rghtMoments[fsgrids::moments::P_33]);

         dMoments[fsgrids::dmoments::dVxdx]  = limiter(leftMoments[fsgrids::moments::RHOVX], leftMoments[fsgrids::moments::RHO],
                                                       centMoments[fsgrids::moments::RHOVX], centMoments[fsgrids::moments::RHO],
                                                       rghtMoments[fsgrids::moments::RHOVX], rghtMoments[fsgrids::moments::RHO]);
         dMoments[fsgrids::dmoments::dVydx]  = limiter(leftMoments[fsgrids::moments::RHOVY], leftMoments[fsgrids::moments::RHO],
                                                       centMoments[fsgrids::moments::RHOVY], centMoments[fsgrids::moments::RHO],
                                                       rghtMoments[fsgrids::moments::RHOVY], rghtMoments[fsgrids::moments::RHO]);
         dMoments[fsgrids::dmoments::dVzdx]  = limiter(leftMoments[fsgrids::moments::RHOVZ], leftMoments[fsgrids::moments::RHO],
                                                       centMoments[fsgrids::moments::RHOVZ], centMoments[fsgrids::moments::RHO],
                                                       rghtMoments[fsgrids::moments::RHOVZ], rghtMoments[fsgrids::moments::RHO]);
         }
      dPerB[fsgrids::dperb::dPERBydx]  = limiter(leftPerB[fsgrids::bfield::PERBY],centPerB[fsgrids::bfield::PERBY],rghtPerB[fsgrids::bfield::PERBY]);
      dPerB[fsgrids::dperb::dPERBzdx]  = limiter(leftPerB[fsgrids::bfield::PERBZ],centPerB[fsgrids::bfield::PERBZ],rghtPerB[fsgrids::bfield::PERBZ]);
      if(Parameters::ohmHallTerm < 2) {
         dPerB[fsgrids::dperb::dPERBydxx] = 0.0;
         dPerB[fsgrids::dperb::dPERBzdxx] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBydxx] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] - 2.0*centPerB[fsgrids::bfield::PERBY];
         dPerB[fsgrids::dperb::dPERBzdxx] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] - 2.0*centPerB[fsgrids::bfield::PERBZ];
      }
   } else {
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 0);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 0);
      }
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         leftPerB = & perBGrid->get(i,j-1,k);
         rghtPerB = & perBGrid->get(i,j+1,k);
         if (doMoments) {
            leftMoments = & momentsGrid->get(i,j-1,k);
            rghtMoments = & momentsGrid->get(i,j+1,k);
         }
      }
      if (RKCase == RK_ORDER2_STEP1) {
         leftPerB = & perBGridDt2->get(i,j-1,k);
         rghtPerB = & perBGridDt2->get(i,j+1,k);
         if (doMoments) {
            leftMoments = & momentsGridDt2->get(i,j-1,k);
            rghtMoments = & momentsGridDt2->get(i,j+1,k);
         }
      }

      if (doMoments) {
         dMoments[fsgrids::dmoments::drhody] = limiter(leftMoments[fsgrids::moments::RHO],centMoments[fsgrids::moments::RHO],rghtMoments[fsgrids::moments::RHO]);
         dMoments[fsgrids::dmoments::dp11dy] = limiter(leftMoments[fsgrids::moments::P_11],centMoments[fsgrids::moments::P_11],rghtMoments[fsgrids::moments::P_11]);
         dMoments[fsgrids::dmoments::dp22dy] = limiter(leftMoments[fsgrids::moments::P_22],centMoments[fsgrids::moments::P_22],rghtMoments[fsgrids::moments::P_22]);
         dMoments[fsgrids::dmoments::dp33dy] = limiter(leftMoments[fsgrids::moments::P_33],centMoments[fsgrids::moments::P_33],rghtMoments[fsgrids::moments::P_33]);
         dMoments[fsgrids::dmoments::dVxdy]  = limiter(leftMoments[fsgrids::moments::RHOVX], leftMoments[fsgrids::moments::RHO],
                                                       centMoments[fsgrids::moments::RHOVX], centMoments[fsgrids::moments::RHO],
                                                       rghtMoments[fsgrids::moments::RHOVX], rghtMoments[fsgrids::moments::RHO]);
         dMoments[fsgrids::dmoments::dVydy]  = limiter(leftMoments[fsgrids::moments::RHOVY], leftMoments[fsgrids::moments::RHO],
                                                       centMoments[fsgrids::moments::RHOVY], centMoments[fsgrids::moments::RHO],
                                                       rghtMoments[fsgrids::moments::RHOVY], rghtMoments[fsgrids::moments::RHO]);
         dMoments[fsgrids::dmoments::dVzdy]  = limiter(leftMoments[fsgrids::moments::RHOVZ], leftMoments[fsgrids::moments::RHO],
                                                       centMoments[fsgrids::moments::RHOVZ], centMoments[fsgrids::moments::RHO],
                                                       rghtMoments[fsgrids::moments::RHOVZ], rghtMoments[fsgrids::moments::RHO]);
      }
      dPerB[fsgrids::dperb::dPERBxdy]  = limiter(leftPerB[fsgrids::bfield::PERBX],centPerB[fsgrids::bfield::PERBX],rghtPerB[fsgrids::bfield::PERBX]);
      dPerB[fsgrids::dperb::dPERBzdy]  = limiter(leftPerB[fsgrids::bfield::PERBZ],centPerB[fsgrids::bfield::PERBZ],rghtPerB[fsgrids::bfield::PERBZ]);

      if(Parameters::ohmHallTerm < 2) {
         dPerB[fsgrids::dperb::dPERBxdyy] = 0.0;
         dPerB[fsgrids::dperb::dPERBzdyy] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBxdyy] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] - 2.0*centPerB[fsgrids::bfield::PERBX];
         dPerB[fsgrids::dperb::dPERBzdyy] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] - 2.0*centPerB[fsgrids::bfield::PERBZ];
      }
   } else {
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 1);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         leftPerB = & perBGrid->get(i,j,k-1);
         rghtPerB = & perBGrid->get(i,j,k+1);
         if (doMoments) {
            leftMoments = & momentsGrid->get(i,j,k-1);
            rghtMoments = & momentsGrid->get(i,j,k+1);
         }
      }
      if (RKCase == RK_ORDER2_STEP1) {
         leftPerB = & perBGridDt2->get(i,j,k-1);
         rghtPerB = & perBGridDt2->get(i,j,k+1);
         if (doMoments) {
            leftMoments = & momentsGridDt2->get(i,j,k-1);
            rghtMoments = & momentsGridDt2->get(i,j,k+1);
         }
      }

      if (doMoments) {
         dMoments[fsgrids::dmoments::drhodz] = limiter(leftMoments[fsgrids::moments::RHO],centMoments[fsgrids::moments::RHO],rghtMoments[fsgrids::moments::RHO]);
         dMoments[fsgrids::dmoments::dp11dz] = limiter(leftMoments[fsgrids::moments::P_11],centMoments[fsgrids::moments::P_11],rghtMoments[fsgrids::moments::P_11]);
         dMoments[fsgrids::dmoments::dp22dz] = limiter(leftMoments[fsgrids::moments::P_22],centMoments[fsgrids::moments::P_22],rghtMoments[fsgrids::moments::P_22]);
         dMoments[fsgrids::dmoments::dp33dz] = limiter(leftMoments[fsgrids::moments::P_33],centMoments[fsgrids::moments::P_33],rghtMoments[fsgrids::moments::P_33]);
         dMoments[fsgrids::dmoments::dVxdz]  = limiter(leftMoments[fsgrids::moments::RHOVX], leftMoments[fsgrids::moments::RHO],
                                                       centMoments[fsgrids::moments::RHOVX], centMoments[fsgrids::moments::RHO],
                                                       rghtMoments[fsgrids::moments::RHOVX], rghtMoments[fsgrids::moments::RHO]);
         dMoments[fsgrids::dmoments::dVydz]  = limiter(leftMoments[fsgrids::moments::RHOVY], leftMoments[fsgrids::moments::RHO],
                                                       centMoments[fsgrids::moments::RHOVY], centMoments[fsgrids::moments::RHO],
                                                       rghtMoments[fsgrids::moments::RHOVY], rghtMoments[fsgrids::moments::RHO]);
         dMoments[fsgrids::dmoments::dVzdz]  = limiter(leftMoments[fsgrids::moments::RHOVZ], leftMoments[fsgrids::moments::RHO],
                                                       centMoments[fsgrids::moments::RHOVZ], centMoments[fsgrids::moments::RHO],
                                                       rghtMoments[fsgrids::moments::RHOVZ], rghtMoments[fsgrids::moments::RHO]);
      }
      dPerB[fsgrids::dperb::dPERBxdz]  = limiter(leftPerB[fsgrids::bfield::PERBX],centPerB[fsgrids::bfield::PERBX],rghtPerB[fsgrids::bfield::PERBX]);
      dPerB[fsgrids::dperb::dPERBydz]  = limiter(leftPerB[fsgrids::bfield::PERBY],centPerB[fsgrids::bfield::PERBY],rghtPerB[fsgrids::bfield::PERBY]);
      if(Parameters::ohmHallTerm < 2) {
         dPerB[fsgrids::dperb::dPERBxdzz] = 0.0;
         dPerB[fsgrids::dperb::dPERBydzz] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBxdzz] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] - 2.0*centPerB[fsgrids::bfield::PERBX];
         dPerB[fsgrids::dperb::dPERBydzz] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] - 2.0*centPerB[fsgrids::bfield::PERBY];
      }
   } else {
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 2);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 2);
      }
   }
   
   if (Parameters::ohmHallTerm < 2) {
      dPerB[fsgrids::dperb::dPERBxdyz] = 0.0;
      dPerB[fsgrids::dperb::dPERBydxz] = 0.0;
      dPerB[fsgrids::dperb::dPERBzdxy] = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
         
         if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            botLeft = & perBGrid->get(i-1,j-1,k);
            botRght = & perBGrid->get(i+1,j-1,k);
            topLeft = & perBGrid->get(i-1,j+1,k);
            topRght = & perBGrid->get(i+1,j+1,k);
         }
         if (RKCase == RK_ORDER2_STEP1) {
            botLeft = & perBGridDt2->get(i-1,j-1,k);
            botRght = & perBGridDt2->get(i+1,j-1,k);
            topLeft = & perBGridDt2->get(i-1,j+1,k);
            topRght = & perBGridDt2->get(i+1,j+1,k);
         }
         
         dPerB[fsgrids::dperb::dPERBzdxy] = FOURTH * (botLeft[fsgrids::bfield::PERBZ] + topRght[fsgrids::bfield::PERBZ] - botRght[fsgrids::bfield::PERBZ] - topLeft[fsgrids::bfield::PERBZ]);
         
      } else {
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 3);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 3);
         }
      }
      
      // Calculate xz mixed derivatives:
      if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
         
         if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            botLeft = & perBGrid->get(i-1,j,k-1);
            botRght = & perBGrid->get(i+1,j,k-1);
            topLeft = & perBGrid->get(i-1,j,k+1);
            topRght = & perBGrid->get(i+1,j,k+1);
         }
         if (RKCase == RK_ORDER2_STEP1) {
            botLeft = & perBGridDt2->get(i-1,j,k-1);
            botRght = & perBGridDt2->get(i+1,j,k-1);
            topLeft = & perBGridDt2->get(i-1,j,k+1);
            topRght = & perBGridDt2->get(i+1,j,k+1);
         }
         
         dPerB[fsgrids::dperb::dPERBydxz] = FOURTH * (botLeft[fsgrids::bfield::PERBY] + topRght[fsgrids::bfield::PERBY] - botRght[fsgrids::bfield::PERBY] - topLeft[fsgrids::bfield::PERBY]);
         
      } else {
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 4);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 4);
         }
      }
      
      // Calculate yz mixed derivatives:
      if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
         
         if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            botLeft = & perBGrid->get(i,j-1,k-1);
            botRght = & perBGrid->get(i,j+1,k-1);
            topLeft = & perBGrid->get(i,j-1,k+1);
            topRght = & perBGrid->get(i,j+1,k+1);
         }
         if (RKCase == RK_ORDER2_STEP1) {
            botLeft = & perBGridDt2->get(i,j-1,k-1);
            botRght = & perBGridDt2->get(i,j+1,k-1);
            topLeft = & perBGridDt2->get(i,j-1,k+1);
            topRght = & perBGridDt2->get(i,j+1,k+1);
         }
         
         dPerB[fsgrids::dperb::dPERBxdyz] = FOURTH * (botLeft[fsgrids::bfield::PERBX] + topRght[fsgrids::bfield::PERBX] - botRght[fsgrids::bfield::PERBX] - topLeft[fsgrids::bfield::PERBX]);
         
      } else {
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
 * \param mpiGrid Grid
 * \param sysBoundaries System boundary conditions existing
 * \param localCells Vector of local cells to process
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param doMoments If true, the derivatives of moments (rho, V, P) are computed.
 
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivativesSimple(
   const FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   const FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   const FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   const FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase,
   const bool& doMoments
) {
   int timer;
   const std::array<int, 3> gridDims = technicalGrid->getLocalSize();
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   
   phiprof::start("Calculate face derivatives");
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);
   
   switch (RKCase) {
      case RK_ORDER1:
         // Means initialising the solver as well as RK_ORDER1
         // standard case exchange PERB*,RHO,RHOV,P with neighbours
         // The update of PERB[XYZ] is needed after the system
         // boundary update of propagateMagneticFieldSimple.
         perBGrid.updateGhostCells();
         momentsGrid.updateGhostCells();
         break;
      case RK_ORDER2_STEP1:
         // Exchange PERB*_DT2,RHO_DT2,RHOV*_DT2,P*DT2 with neighbours The
         // update of PERB[XYZ]_DT2 is needed after the system
         // boundary update of propagateMagneticFieldSimple.
         perBDt2Grid.updateGhostCells();
         momentsDt2Grid.updateGhostCells();
         break;
      case RK_ORDER2_STEP2:
         // Exchange PERB*,RHO,RHOV*,P* with neighbours The update of B
         // is needed after the system boundary update of
         // propagateMagneticFieldSimple.
         perBGrid.updateGhostCells();
         momentsGrid.updateGhostCells();
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
   for (uint k=0; k<gridDims[2]; k++) {
      for (uint j=0; j<gridDims[1]; j++) {
         for (uint i=0; i<gridDims[0]; i++) {
            if (technicalGrid->get(i,j,k).sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
            
            calculateDerivatives(i,j,k, perBGrid, perBDt2Grid, EGrid, EDt2Grid, momentsGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RKCase, doMoments);
         }
      }
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   phiprof::stop("Calculate face derivatives",N_cells,"Spatial Cells");   
}

/*! \brief Low-level spatial derivatives calculation.
 * 
 * For the cell with ID cellID calculate the spatial derivatives of BVOL or apply the derivative boundary conditions defined in project.h.
 * \param mpiGrid Grid
 * \param cache Field solver cell cache
 * \param cells Vector of local cells to process
 * 
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateBVOLDerivatives(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 3, 2> & volGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   const int i,
   const int j,
   const int k,
   SysBoundary& sysBoundaries
) {
   Real* const array = volGrid->get(i,j,k);
   
   std::array<Real, fsgrids::volfields::N_VOL> * left = NULL;
   std::array<Real, fsgrids::volfields::N_VOL> * rght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if (technicalGrid->get(i,j,k).sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      left = & volGrid->get(i-1,j,k);
      rght = & volGrid->get(i+1,j,k);
      
      array[fsgrids::volfields::dPERBYVOLdx] = limiter(left[fsgrids::volfields::PERBYVOL],cent[fsgrids::volfields::PERBYVOL],rght[fsgrids::volfields::PERBYVOL]);
      array[fsgrids::volfields::dPERBZVOLdx] = limiter(left[fsgrids::volfields::PERBZVOL],cent[fsgrids::volfields::PERBZVOL],rght[fsgrids::volfields::PERBZVOL]);
   } else {
      if (technicalGrid->get(i,j,k).sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 0);
      } else {
         sysBoundaries.getSysBoundary(technicalGrid->get(i,j,k).sysBoundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 0);
      }
   }
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
   if (technicalGrid->get(i,j,k).sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      left = & volGrid->get(i,j-1,k);
      rght = & volGrid->get(i,j+1,k);
      
      array[fsgrids::volfields::dPERBXVOLdy] = limiter(left[fsgrids::volfields::PERBXVOL],cent[fsgrids::volfields::PERBXVOL],rght[fsgrids::volfields::PERBXVOL]);
      array[fsgrids::volfields::dPERBZVOLdy] = limiter(left[fsgrids::volfields::PERBZVOL],cent[fsgrids::volfields::PERBZVOL],rght[fsgrids::volfields::PERBZVOL]);
   } else {
      if (technicalGrid->get(i,j,k).sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 1);
      } else {
         sysBoundaries.getSysBoundary(technicalGrid->get(i,j,k).sysBoundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if (technicalGrid->get(i,j,k).sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      left = cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1-1)]->parameters;
      rght = cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->parameters;
      
      array[fsgrids::volfields::dPERBXVOLdz] = limiter(left[fsgrids::volfields::PERBXVOL],cent[fsgrids::volfields::PERBXVOL],rght[fsgrids::volfields::PERBXVOL]);
      array[fsgrids::volfields::dPERBYVOLdz] = limiter(left[fsgrids::volfields::PERBYVOL],cent[fsgrids::volfields::PERBYVOL],rght[fsgrids::volfields::PERBYVOL]);
   } else {
      if (technicalGrid->get(i,j,k).sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 2);
      } else {
         sysBoundaries.getSysBoundary(technicalGrid->get(i,j,k).sysBoundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 2);
      }
   }
}

/*! \brief High-level derivative calculation wrapper function.
 * 
 * BVOL has been calculated locally by calculateVolumeAveragedFields but not communicated.
 * For the acceleration step one needs the cross-derivatives of BVOL
 * 
 * \param mpiGrid Grid
 * \param sysBoundaries System boundary conditions existing
 * \param localCells Vector of local cells to process
 * 
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateBVOLDerivativesSimple(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 3, 2> & volGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells
) {
   int timer;
   const std::array<int, 3> gridDims = technicalGrid->getLocalSize();
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
   for (uint k=0; k<gridDims[2]; k++) {
      for (uint j=0; j<gridDims[1]; j++) {
         for (uint i=0; i<gridDims[0]; i++) {
            if (technicalGrid->get(i,j,k).sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
            
            calculateBVOLDerivatives(volGrid,i,j,k,sysBoundaries);
         }
      }
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");

   phiprof::stop("Calculate volume derivatives",N_cells,"Spatial Cells");
}
