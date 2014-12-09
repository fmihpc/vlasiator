/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef DERIVATIVES_HPP
#define DERIVATIVES_HPP

#include "fs_limiters.h"

/*! \brief Low-level spatial derivatives calculation.
 * 
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in project.h. Uses RHO, RHOV[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the second-order method, and RHO_DT2, RHOV[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 * \param cellID Index of the cell to process
 * \param mpiGrid Grid
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param doMoments If true, the derivatives of moments (rho, V, P) are computed.
 */
static void calculateDerivatives(
   const CellID& cellID,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase,
   const bool& doMoments
) {
   namespace cp = CellParams;
   namespace fs = fieldsolver;
   Real* const array       = mpiGrid[cellID]->derivatives;
   // Get boundary flag for the cell:
   #ifndef NDEBUG
      map<CellID,uint>::const_iterator it = sysBoundaryFlags.find(cellID);
      if (it == sysBoundaryFlags.end()) {cerr << "ERROR Could not find boundary flag for cell #" << cellID << endl; exit(1);}
      cuint existingCells = it->second;
   #else
      cuint existingCells = sysBoundaryFlags[cellID];
   #endif
   cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
   
   CellID leftNbrID,rghtNbrID;
   creal* left = NULL;
   creal* cent = mpiGrid[cellID   ]->parameters;
   #ifdef DEBUG_SOLVERS
   if (cent[cp::RHO] <= 0) {
      std::cerr << __FILE__ << ":" << __LINE__
         << (cent[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << cellID
         << std::endl;
      abort();
   }
   #endif
   creal* rght = NULL;
   CellID botLeftNbrID, botRghtNbrID, topLeftNbrID, topRghtNbrID;
   creal* botLeft = NULL;
   creal* botRght = NULL;
   creal* topLeft = NULL;
   creal* topRght = NULL;
   
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if (((existingCells & CALCULATE_DX) == CALCULATE_DX) &&
       ((mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ||
        (mpiGrid[cellID]->sysBoundaryLayer == 1))
      ) {
      leftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2  ,2  );
      rghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2  ,2  );
      left = mpiGrid[leftNbrID]->parameters;
      #ifdef DEBUG_SOLVERS
      if (left[cp::RHO] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (left[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << leftNbrID
            << std::endl;
         abort();
      }
      #endif
      rght = mpiGrid[rghtNbrID]->parameters;
      #ifdef DEBUG_SOLVERS
      if (rght[cp::RHO] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (rght[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << rghtNbrID
            << std::endl;
         abort();
      }
      #endif
      
      if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         if (doMoments) {
            array[fs::drhodx] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
            array[fs::dp11dx] = limiter(left[cp::P_11],cent[cp::P_11],rght[cp::P_11]);
            array[fs::dp22dx] = limiter(left[cp::P_22],cent[cp::P_22],rght[cp::P_22]);
            array[fs::dp33dx] = limiter(left[cp::P_33],cent[cp::P_33],rght[cp::P_33]);
            array[fs::dVxdx]  = limiter(divideIfNonZero(left[cp::RHOVX], left[cp::RHO]),
                                       divideIfNonZero(cent[cp::RHOVX], cent[cp::RHO]),
                                       divideIfNonZero(rght[cp::RHOVX], rght[cp::RHO]));
            array[fs::dVydx]  = limiter(divideIfNonZero(left[cp::RHOVY], left[cp::RHO]),
                                       divideIfNonZero(cent[cp::RHOVY], cent[cp::RHO]),
                                       divideIfNonZero(rght[cp::RHOVY], rght[cp::RHO]));
            array[fs::dVzdx]  = limiter(divideIfNonZero(left[cp::RHOVZ], left[cp::RHO]),
                                       divideIfNonZero(cent[cp::RHOVZ], cent[cp::RHO]),
                                       divideIfNonZero(rght[cp::RHOVZ], rght[cp::RHO]));
         }
         array[fs::dPERBydx]  = limiter(left[cp::PERBY],cent[cp::PERBY],rght[cp::PERBY]);
         array[fs::dPERBzdx]  = limiter(left[cp::PERBZ],cent[cp::PERBZ],rght[cp::PERBZ]);
         if(Parameters::ohmHallTerm < 2) {
            array[fs::dPERBydxx] = 0.0;
            array[fs::dPERBzdxx] = 0.0;
         } else {
            array[fs::dPERBydxx] = left[cp::PERBY] + rght[cp::PERBY] - 2.0*cent[cp::PERBY];
            array[fs::dPERBzdxx] = left[cp::PERBZ] + rght[cp::PERBZ] - 2.0*cent[cp::PERBZ];
         }
      }
      if (RKCase == RK_ORDER2_STEP1) {
         if (doMoments) {
            array[fs::drhodx] = limiter(left[cp::RHO_DT2],cent[cp::RHO_DT2],rght[cp::RHO_DT2]);
            array[fs::dp11dx] = limiter(left[cp::P_11_DT2],cent[cp::P_11_DT2],rght[cp::P_11_DT2]);
            array[fs::dp22dx] = limiter(left[cp::P_22_DT2],cent[cp::P_22_DT2],rght[cp::P_22_DT2]);
            array[fs::dp33dx] = limiter(left[cp::P_33_DT2],cent[cp::P_33_DT2],rght[cp::P_33_DT2]);
            array[fs::dVxdx]  = limiter(divideIfNonZero(left[cp::RHOVX_DT2], left[cp::RHO_DT2]),
                                       divideIfNonZero(cent[cp::RHOVX_DT2], cent[cp::RHO_DT2]),
                                       divideIfNonZero(rght[cp::RHOVX_DT2], rght[cp::RHO_DT2]));
            array[fs::dVydx]  = limiter(divideIfNonZero(left[cp::RHOVY_DT2], left[cp::RHO_DT2]),
                                       divideIfNonZero(cent[cp::RHOVY_DT2], cent[cp::RHO_DT2]),
                                       divideIfNonZero(rght[cp::RHOVY_DT2], rght[cp::RHO_DT2]));
            array[fs::dVzdx]  = limiter(divideIfNonZero(left[cp::RHOVZ_DT2], left[cp::RHO_DT2]),
                                       divideIfNonZero(cent[cp::RHOVZ_DT2], cent[cp::RHO_DT2]),
                                       divideIfNonZero(rght[cp::RHOVZ_DT2], rght[cp::RHO_DT2]));
         }
         array[fs::dPERBydx]  = limiter(left[cp::PERBY_DT2],cent[cp::PERBY_DT2],rght[cp::PERBY_DT2]);
         array[fs::dPERBzdx]  = limiter(left[cp::PERBZ_DT2],cent[cp::PERBZ_DT2],rght[cp::PERBZ_DT2]);
         if(Parameters::ohmHallTerm < 2) {
            array[fs::dPERBydxx] = 0.0;
            array[fs::dPERBzdxx] = 0.0;
         } else {
            array[fs::dPERBydxx] = left[cp::PERBY_DT2] + rght[cp::PERBY_DT2] - 2.0*cent[cp::PERBY_DT2];
            array[fs::dPERBzdxx] = left[cp::PERBZ_DT2] + rght[cp::PERBZ_DT2] - 2.0*cent[cp::PERBZ_DT2];
         }
      }
   } else {
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 0);
      } else {
         sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)
            ->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 0);
      }
   }
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
   if (((existingCells & CALCULATE_DY) == CALCULATE_DY) &&
       ((mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ||
        (mpiGrid[cellID]->sysBoundaryLayer == 1))) {
      leftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2-1,2  );
      rghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2+1,2  );
      
      left = mpiGrid[leftNbrID]->parameters;
      #ifdef DEBUG_SOLVERS
      if (left[cp::RHO] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (left[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << leftNbrID
            << " Zero density in spatial cell " << leftNbrID
            << std::endl;
         abort();
      }
      #endif
      
      rght = mpiGrid[rghtNbrID]->parameters;
      #ifdef DEBUG_SOLVERS
      if (rght[cp::RHO] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (rght[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << rghtNbrID
            << std::endl;
         abort();
      }
      #endif
      
      if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         if (doMoments) {
            array[fs::drhody] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
            array[fs::dp11dy] = limiter(left[cp::P_11],cent[cp::P_11],rght[cp::P_11]);
            array[fs::dp22dy] = limiter(left[cp::P_22],cent[cp::P_22],rght[cp::P_22]);
            array[fs::dp33dy] = limiter(left[cp::P_33],cent[cp::P_33],rght[cp::P_33]);
            array[fs::dVxdy]  = limiter(divideIfNonZero(left[cp::RHOVX], left[cp::RHO]),
                                       divideIfNonZero(cent[cp::RHOVX], cent[cp::RHO]),
                                       divideIfNonZero(rght[cp::RHOVX], rght[cp::RHO]));
            array[fs::dVydy]  = limiter(divideIfNonZero(left[cp::RHOVY], left[cp::RHO]),
                                       divideIfNonZero(cent[cp::RHOVY], cent[cp::RHO]),
                                       divideIfNonZero(rght[cp::RHOVY], rght[cp::RHO]));
            array[fs::dVzdy]  = limiter(divideIfNonZero(left[cp::RHOVZ], left[cp::RHO]),
                                       divideIfNonZero(cent[cp::RHOVZ], cent[cp::RHO]),
                                       divideIfNonZero(rght[cp::RHOVZ], rght[cp::RHO]));
         }
         array[fs::dPERBxdy]  = limiter(left[cp::PERBX],cent[cp::PERBX],rght[cp::PERBX]);
         array[fs::dPERBzdy]  = limiter(left[cp::PERBZ],cent[cp::PERBZ],rght[cp::PERBZ]);
         if(Parameters::ohmHallTerm < 2) {
            array[fs::dPERBxdyy] = 0.0;
            array[fs::dPERBzdyy] = 0.0;
         } else {
            array[fs::dPERBxdyy] = left[cp::PERBX] + rght[cp::PERBX] - 2.0*cent[cp::PERBX];
            array[fs::dPERBzdyy] = left[cp::PERBZ] + rght[cp::PERBZ] - 2.0*cent[cp::PERBZ];
         }
      }
      if (RKCase == RK_ORDER2_STEP1) {
         if (doMoments) {
            array[fs::drhody] = limiter(left[cp::RHO_DT2],cent[cp::RHO_DT2],rght[cp::RHO_DT2]);
            array[fs::dp11dy] = limiter(left[cp::P_11_DT2],cent[cp::P_11_DT2],rght[cp::P_11_DT2]);
            array[fs::dp22dy] = limiter(left[cp::P_22_DT2],cent[cp::P_22_DT2],rght[cp::P_22_DT2]);
            array[fs::dp33dy] = limiter(left[cp::P_33_DT2],cent[cp::P_33_DT2],rght[cp::P_33_DT2]);
            array[fs::dVxdy]  = limiter(divideIfNonZero(left[cp::RHOVX_DT2], left[cp::RHO_DT2]),
                                       divideIfNonZero(cent[cp::RHOVX_DT2], cent[cp::RHO_DT2]),
                                       divideIfNonZero(rght[cp::RHOVX_DT2], rght[cp::RHO_DT2]));
            array[fs::dVydy]  = limiter(divideIfNonZero(left[cp::RHOVY_DT2], left[cp::RHO_DT2]),
                                       divideIfNonZero(cent[cp::RHOVY_DT2], cent[cp::RHO_DT2]),
                                       divideIfNonZero(rght[cp::RHOVY_DT2], rght[cp::RHO_DT2]));
            array[fs::dVzdy]  = limiter(divideIfNonZero(left[cp::RHOVZ_DT2], left[cp::RHO_DT2]),
                                       divideIfNonZero(cent[cp::RHOVZ_DT2], cent[cp::RHO_DT2]),
                                       divideIfNonZero(rght[cp::RHOVZ_DT2], rght[cp::RHO_DT2]));
         }
         array[fs::dPERBxdy]  = limiter(left[cp::PERBX_DT2],cent[cp::PERBX_DT2],rght[cp::PERBX_DT2]);
         array[fs::dPERBzdy]  = limiter(left[cp::PERBZ_DT2],cent[cp::PERBZ_DT2],rght[cp::PERBZ_DT2]);
         if(Parameters::ohmHallTerm < 2) {
            array[fs::dPERBxdyy] = 0.0;
            array[fs::dPERBzdyy] = 0.0;
         } else {
            array[fs::dPERBxdyy] = left[cp::PERBX_DT2] + rght[cp::PERBX_DT2] - 2.0*cent[cp::PERBX_DT2];
            array[fs::dPERBzdyy] = left[cp::PERBZ_DT2] + rght[cp::PERBZ_DT2] - 2.0*cent[cp::PERBZ_DT2];
         }
      }
   } else {
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 1);
      } else {
         sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)
         ->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if (((existingCells & CALCULATE_DZ) == CALCULATE_DZ) &&
       ((mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ||
        (mpiGrid[cellID]->sysBoundaryLayer == 1))) {
      leftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2  ,2-1);
      rghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2  ,2+1);
      left = mpiGrid[leftNbrID]->parameters;
      
      #ifdef DEBUG_SOLVERS
      if (left[cp::RHO] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (left[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << leftNbrID
            << std::endl;
         abort();
      }
      #endif
      rght = mpiGrid[rghtNbrID]->parameters;
      #ifdef DEBUG_SOLVERS
      if (rght[cp::RHO] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (rght[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << rghtNbrID
            << std::endl;
         abort();
      }
      #endif
      
      if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         if (doMoments) {
            array[fs::drhodz] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
            array[fs::dp11dz] = limiter(left[cp::P_11],cent[cp::P_11],rght[cp::P_11]);
            array[fs::dp22dz] = limiter(left[cp::P_22],cent[cp::P_22],rght[cp::P_22]);
            array[fs::dp33dz] = limiter(left[cp::P_33],cent[cp::P_33],rght[cp::P_33]);
            array[fs::dVxdz]  = limiter(divideIfNonZero(left[cp::RHOVX], left[cp::RHO]),
                                       divideIfNonZero(cent[cp::RHOVX], cent[cp::RHO]),
                                       divideIfNonZero(rght[cp::RHOVX], rght[cp::RHO]));
            array[fs::dVydz]  = limiter(divideIfNonZero(left[cp::RHOVY], left[cp::RHO]),
                                       divideIfNonZero(cent[cp::RHOVY], cent[cp::RHO]),
                                       divideIfNonZero(rght[cp::RHOVY], rght[cp::RHO]));
            array[fs::dVzdz]  = limiter(divideIfNonZero(left[cp::RHOVZ], left[cp::RHO]),
                                       divideIfNonZero(cent[cp::RHOVZ], cent[cp::RHO]),
                                       divideIfNonZero(rght[cp::RHOVZ], rght[cp::RHO]));
         }
         array[fs::dPERBxdz]  = limiter(left[cp::PERBX],cent[cp::PERBX],rght[cp::PERBX]);
         array[fs::dPERBydz]  = limiter(left[cp::PERBY],cent[cp::PERBY],rght[cp::PERBY]);
         if(Parameters::ohmHallTerm < 2) {
            array[fs::dPERBxdzz] = 0.0;
            array[fs::dPERBydzz] = 0.0;
         } else {
            array[fs::dPERBxdzz] = left[cp::PERBX] + rght[cp::PERBX] - 2.0*cent[cp::PERBX];
            array[fs::dPERBydzz] = left[cp::PERBY] + rght[cp::PERBY] - 2.0*cent[cp::PERBY];
         }
      }
      if (RKCase == RK_ORDER2_STEP1) {
         if (doMoments) {
            array[fs::drhodz] = limiter(left[cp::RHO_DT2],cent[cp::RHO_DT2],rght[cp::RHO_DT2]);
            array[fs::dp11dz] = limiter(left[cp::P_11_DT2],cent[cp::P_11_DT2],rght[cp::P_11_DT2]);
            array[fs::dp22dz] = limiter(left[cp::P_22_DT2],cent[cp::P_22_DT2],rght[cp::P_22_DT2]);
            array[fs::dp33dz] = limiter(left[cp::P_33_DT2],cent[cp::P_33_DT2],rght[cp::P_33_DT2]);
            array[fs::dVxdz]  = limiter(divideIfNonZero(left[cp::RHOVX_DT2], left[cp::RHO_DT2]),
                                       divideIfNonZero(cent[cp::RHOVX_DT2], cent[cp::RHO_DT2]),
                                       divideIfNonZero(rght[cp::RHOVX_DT2], rght[cp::RHO_DT2]));
            array[fs::dVydz]  = limiter(divideIfNonZero(left[cp::RHOVY_DT2], left[cp::RHO_DT2]),
                                       divideIfNonZero(cent[cp::RHOVY_DT2], cent[cp::RHO_DT2]),
                                       divideIfNonZero(rght[cp::RHOVY_DT2], rght[cp::RHO_DT2]));
            array[fs::dVzdz]  = limiter(divideIfNonZero(left[cp::RHOVZ_DT2], left[cp::RHO_DT2]),
                                       divideIfNonZero(cent[cp::RHOVZ_DT2], cent[cp::RHO_DT2]),
                                       divideIfNonZero(rght[cp::RHOVZ_DT2], rght[cp::RHO_DT2]));
         }
         array[fs::dPERBxdz]  = limiter(left[cp::PERBX_DT2],cent[cp::PERBX_DT2],rght[cp::PERBX_DT2]);
         array[fs::dPERBydz]  = limiter(left[cp::PERBY_DT2],cent[cp::PERBY_DT2],rght[cp::PERBY_DT2]);
         if(Parameters::ohmHallTerm < 2) {
            array[fs::dPERBxdzz] = 0.0;
            array[fs::dPERBydzz] = 0.0;
         } else {
            array[fs::dPERBxdzz] = left[cp::PERBX_DT2] + rght[cp::PERBX_DT2] - 2.0*cent[cp::PERBX_DT2];
            array[fs::dPERBydzz] = left[cp::PERBY_DT2] + rght[cp::PERBY_DT2] - 2.0*cent[cp::PERBY_DT2];
         }
      }
   } else {
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 2);
      } else {
         sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)
         ->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 2);
      }
   }
   
   if(Parameters::ohmHallTerm < 2) {
      array[fs::dPERBxdyz] = 0.0;
      array[fs::dPERBydxz] = 0.0;
      array[fs::dPERBzdxy] = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      if (((existingCells & CALCULATE_DXY) == CALCULATE_DXY) &&
         ((mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ||
         (mpiGrid[cellID]->sysBoundaryLayer == 1))
      ) {
         botLeftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2-1,2  );
         botRghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2-1,2  );
         topLeftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2+1,2  );
         topRghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2+1,2  );
         botLeft = mpiGrid[botLeftNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (botLeft[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (botLeft[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << botLeftNbrID
            << std::endl;
            abort();
         }
         #endif
         botRght = mpiGrid[botRghtNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (botRght[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (botRght[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << botRghtNbrID
            << std::endl;
            abort();
         }
         #endif
         topLeft = mpiGrid[topLeftNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (topLeft[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (topLeft[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << topLeftNbrID
            << std::endl;
            abort();
         }
         #endif
         topRght = mpiGrid[topRghtNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (topRght[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (topRght[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << topRghtNbrID
            << std::endl;
            abort();
         }
         #endif
         
         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            array[fs::dPERBzdxy] = FOURTH * (botLeft[cp::PERBZ] + topRght[cp::PERBZ] - botRght[cp::PERBZ] - topLeft[cp::PERBZ]);
         }
         if (RKCase == RK_ORDER2_STEP1) {
            array[fs::dPERBzdxy] = FOURTH * (botLeft[cp::PERBZ_DT2] + topRght[cp::PERBZ_DT2] - botRght[cp::PERBZ_DT2] - topLeft[cp::PERBZ_DT2]);
         }
      } else {
         if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 3);
         } else {
            sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)
               ->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 3);
         }
      }
      
      // Calculate xz mixed derivatives:
      if (((existingCells & CALCULATE_DXZ) == CALCULATE_DXZ) &&
         ((mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ||
         (mpiGrid[cellID]->sysBoundaryLayer == 1))
      ) {
         botLeftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2  ,2-1);
         botRghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2  ,2-1);
         topLeftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2  ,2+1);
         topRghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2  ,2+1);
         botLeft = mpiGrid[botLeftNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (botLeft[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (botLeft[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << botLeftNbrID
            << std::endl;
            abort();
         }
         #endif
         botRght = mpiGrid[botRghtNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (botRght[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (botRght[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << botRghtNbrID
            << std::endl;
            abort();
         }
         #endif
         topLeft = mpiGrid[topLeftNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (topLeft[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (topLeft[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << topLeftNbrID
            << std::endl;
            abort();
         }
         #endif
         topRght = mpiGrid[topRghtNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (topRght[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (topRght[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << topRghtNbrID
            << std::endl;
            abort();
         }
         #endif
         
         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            array[fs::dPERBydxz] = FOURTH * (botLeft[cp::PERBY] + topRght[cp::PERBY] - botRght[cp::PERBY] - topLeft[cp::PERBY]);
         }
         if (RKCase == RK_ORDER2_STEP1) {
            array[fs::dPERBydxz] = FOURTH * (botLeft[cp::PERBY_DT2] + topRght[cp::PERBY_DT2] - botRght[cp::PERBY_DT2] - topLeft[cp::PERBY_DT2]);
         }
      } else {
         if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 4);
         } else {
            sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)
               ->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 4);
         }
      }
      
      // Calculate yz mixed derivatives:
      if (((existingCells & CALCULATE_DYZ) == CALCULATE_DYZ) &&
         ((mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ||
         (mpiGrid[cellID]->sysBoundaryLayer == 1))
      ) {
         botLeftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2-1,2-1);
         botRghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2+1,2-1);
         topLeftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2-1,2+1);
         topRghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2+1,2+1);
         botLeft = mpiGrid[botLeftNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (botLeft[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (botLeft[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << botLeftNbrID
            << std::endl;
            abort();
         }
         #endif
         botRght = mpiGrid[botRghtNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (botRght[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (botRght[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << botRghtNbrID
            << std::endl;
            abort();
         }
         #endif
         topLeft = mpiGrid[topLeftNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (topLeft[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (topLeft[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << topLeftNbrID
            << std::endl;
            abort();
         }
         #endif
         topRght = mpiGrid[topRghtNbrID]->parameters;
         #ifdef DEBUG_SOLVERS
         if (topRght[cp::RHO] <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__
            << (topRght[cp::RHO] < 0 ? " Negative" : " Zero") << " density in spatial cell " << topRghtNbrID
            << std::endl;
            abort();
         }
         #endif
         
         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            array[fs::dPERBxdyz] = FOURTH * (botLeft[cp::PERBX] + topRght[cp::PERBX] - botRght[cp::PERBX] - topLeft[cp::PERBX]);
         }
         if (RKCase == RK_ORDER2_STEP1) {
            array[fs::dPERBxdyz] = FOURTH * (botLeft[cp::PERBX_DT2] + topRght[cp::PERBX_DT2] - botRght[cp::PERBX_DT2] - topLeft[cp::PERBX_DT2]);
         }
      } else {
         if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 5);
         } else {
            sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)
               ->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 5);
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
 * \param localCells Vector of local cells to process
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param doMoments If true, the derivatives of moments (rho, V, P) are computed.
 
 * \sa calculateDerivatives
 */
void calculateDerivativesSimple(
       dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
       SysBoundary& sysBoundaries,
       const vector<CellID>& localCells,
       cint& RKCase,
       const bool& doMoments
) {
   int timer;
   namespace fs = fieldsolver;
   
   phiprof::start("Calculate face derivatives");
   
   switch(RKCase) {
      case RK_ORDER1:
         // Means initialising the solver as well as RK_ORDER1
         // standard case Exchange PERB* with neighbours
         // The update of PERB[XYZ] is needed after the system
         // boundary update of propagateMagneticFieldSimple.
         SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB | Transfer::CELL_RHO_RHOV | Transfer::CELL_P);
         break;
      case RK_ORDER2_STEP1:
         // Exchange PERB*_DT2,RHO_DT2,RHOV*_DT2 with neighbours The
         // update of PERB[XYZ]_DT2 is needed after the system
         // boundary update of propagateMagneticFieldSimple.
         SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERBDT2 | Transfer::CELL_RHODT2_RHOVDT2 | Transfer::CELL_PDT2);
         break;
      case RK_ORDER2_STEP2:
         // Exchange PERB*,RHO,RHOV* with neighbours The update of B
         // is needed after the system boundary update of
         // propagateMagneticFieldSimple.
         SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB | Transfer::CELL_RHO_RHOV | Transfer::CELL_P);
         break;
      default:
         cerr << __FILE__ << ":" << __LINE__ << " Went through switch, this should not happen." << endl;
         abort();
   }
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);
   mpiGrid.start_remote_neighbor_copy_updates(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Compute process inner cells");
   phiprof::start(timer);
   // Calculate derivatives on process inner cells
   const vector<uint64_t> cellsWithLocalNeighbours
      = mpiGrid.get_local_cells_not_on_process_boundary(FIELD_SOLVER_NEIGHBORHOOD_ID);
#pragma omp parallel for
   for(uint i=0;i<cellsWithLocalNeighbours.size();i++){
      const CellID cellID = cellsWithLocalNeighbours[i];
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      calculateDerivatives(cellID, mpiGrid, sysBoundaries, RKCase, true);
   }
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_receives(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);
   
   // Calculate derivatives on process boundary cells
   timer=phiprof::initializeTimer("Compute process boundary cells");
   phiprof::start(timer);
   const vector<uint64_t> cellsWithRemoteNeighbours
      = mpiGrid.get_local_cells_on_process_boundary(FIELD_SOLVER_NEIGHBORHOOD_ID);
#pragma omp parallel for
   for(uint i=0;i<cellsWithRemoteNeighbours.size();i++){
      const CellID cellID = cellsWithRemoteNeighbours[i];
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      calculateDerivatives(cellID, mpiGrid, sysBoundaries, RKCase, true);
   }
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_sends();
   phiprof::stop(timer);
   
   phiprof::stop("Calculate face derivatives");
}

/*! \brief Low-level spatial derivatives calculation.
 * 
 * For the cell with ID cellID calculate the spatial derivatives of BVOL or apply the derivative boundary conditions defined in project.h.
 * \param cellID Index of the cell to process
 * \param mpiGrid Grid
 */
static void calculateBVOLDerivatives(
   const CellID& cellID,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries
) {
   namespace cp = CellParams;
   namespace der = bvolderivatives;
   Real* const array = mpiGrid[cellID]->derivativesBVOL;
   // Get boundary flag for the cell:
   #ifndef NDEBUG
   map<CellID,uint>::const_iterator it = sysBoundaryFlags.find(cellID);
   if (it == sysBoundaryFlags.end()) {cerr << "ERROR Could not find boundary flag for cell #" << cellID << endl; exit(1);}
   cuint existingCells = it->second;
   #else
   cuint existingCells = sysBoundaryFlags[cellID];
   #endif
   cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
   
   CellID leftNbrID,rghtNbrID;
   creal* left = NULL;
   creal* cent = mpiGrid[cellID]->parameters;
   creal* rght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if (((existingCells & CALCULATE_DX) == CALCULATE_DX) &&
      (mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)) {
      leftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2  ,2  );
   rghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2  ,2  );
   left = mpiGrid[leftNbrID]->parameters;
   rght = mpiGrid[rghtNbrID]->parameters;
   
   array[der::dPERBYVOLdx] = limiter(left[cp::PERBYVOL],cent[cp::PERBYVOL],rght[cp::PERBYVOL]);
   array[der::dPERBZVOLdx] = limiter(left[cp::PERBZVOL],cent[cp::PERBZVOL],rght[cp::PERBZVOL]);
   
      } else {
         if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(mpiGrid, cellID, 0);
         } else {
            sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)
            ->fieldSolverBoundaryCondBVOLDerivatives(mpiGrid, cellID, 0);
         }
      }
      
      // Calculate y-derivatives (is not TVD for AMR mesh):
      if (((existingCells & CALCULATE_DY) == CALCULATE_DY) &&
         (mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)) {
         leftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2-1,2  );
      rghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2+1,2  );
      left = mpiGrid[leftNbrID]->parameters;
      rght = mpiGrid[rghtNbrID]->parameters;
      
      array[der::dPERBXVOLdy] = limiter(left[cp::PERBXVOL],cent[cp::PERBXVOL],rght[cp::PERBXVOL]);
      array[der::dPERBZVOLdy] = limiter(left[cp::PERBZVOL],cent[cp::PERBZVOL],rght[cp::PERBZVOL]);
      
         } else {
            if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
               SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(mpiGrid, cellID, 1);
            } else {
               sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)
               ->fieldSolverBoundaryCondBVOLDerivatives(mpiGrid, cellID, 1);
            }
         }
         
         // Calculate z-derivatives (is not TVD for AMR mesh):
         if (((existingCells & CALCULATE_DZ) == CALCULATE_DZ) &&
            (mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)) {
            leftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2  ,2-1);
         rghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2  ,2+1);
         left = mpiGrid[leftNbrID]->parameters;
         rght = mpiGrid[rghtNbrID]->parameters;
         
         array[der::dPERBXVOLdz] = limiter(left[cp::PERBXVOL],cent[cp::PERBXVOL],rght[cp::PERBXVOL]);
         array[der::dPERBYVOLdz] = limiter(left[cp::PERBYVOL],cent[cp::PERBYVOL],rght[cp::PERBYVOL]);
         
            } else {
               if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(mpiGrid, cellID, 2);
               } else {
                  sysBoundaries.getSysBoundary(mpiGrid[cellID]->sysBoundaryFlag)
                  ->fieldSolverBoundaryCondBVOLDerivatives(mpiGrid, cellID, 2);
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
 * \sa calculateDerivatives
 */
void calculateBVOLDerivativesSimple(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells
) {
   int timer;
   namespace fs = fieldsolver;
   
   phiprof::start("Calculate volume derivatives");
   
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_BVOL);
   mpiGrid.update_copies_of_remote_neighbors(FIELD_SOLVER_NEIGHBORHOOD_ID);
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);
   mpiGrid.start_remote_neighbor_copy_updates(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Compute process inner cells");
   phiprof::start(timer);
   // Calculate derivatives on process inner cells
   const vector<uint64_t> cellsWithLocalNeighbours
   = mpiGrid.get_local_cells_not_on_process_boundary(FIELD_SOLVER_NEIGHBORHOOD_ID);
   for (vector<uint64_t>::const_iterator cell = cellsWithLocalNeighbours.begin(); cell != cellsWithLocalNeighbours.end(); cell++) {
      if(mpiGrid[*cell]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      calculateBVOLDerivatives(*cell, mpiGrid, sysBoundaries);
   }
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_receives(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);
   
   // Calculate derivatives on process boundary cells
   timer=phiprof::initializeTimer("Compute process boundary cells");
   phiprof::start(timer);
   const vector<uint64_t> cellsWithRemoteNeighbours
   = mpiGrid.get_local_cells_on_process_boundary(FIELD_SOLVER_NEIGHBORHOOD_ID);
   for (vector<uint64_t>::const_iterator cell = cellsWithRemoteNeighbours.begin(); cell != cellsWithRemoteNeighbours.end(); cell++) {
      if(mpiGrid[*cell]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      calculateBVOLDerivatives(*cell, mpiGrid, sysBoundaries);
   }
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_sends();
   phiprof::stop(timer);
   
   phiprof::stop("Calculate volume derivatives");
}

#endif
