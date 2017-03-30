/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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
#include "fs_cache.h"

extern map<CellID,uint> existingCellsFlags; /**< Defined in fs_common.cpp */

/*! \brief Low-level spatial derivatives calculation.
 * 
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in project.h. Uses RHO, RHOV[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the second-order method, and RHO_DT2, RHOV[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 * \param mpiGrid Grid
 * \param cellCache Field solver cell cache
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa calculateDerivativesSimple calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivatives(
   dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   fs_cache::CellCache& cellCache,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {

   namespace cp = CellParams;
   namespace fs = fieldsolver;
   const CellID cellID = cellCache.cellID;
   
   #ifdef DEBUG_FSOLVER
   bool ok = true;
   if (cellCache.cells[fs_cache::calculateNbrID(1,1,1)] == NULL) ok = false;
   if (cellCache.cells[fs_cache::calculateNbrID(1,1,1)]->derivatives == NULL) ok = false;
   if (ok == false) {
      cerr << "ERROR, NULL pointer detected in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
   #endif
   
   Real* const array = cellCache.cells[fs_cache::calculateNbrID(1,1,1)]->derivatives;

   // Get boundary flag for the cell:
   cuint existingCells    = cellCache.existingCellsFlags;
   cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
   cuint sysBoundaryFlag  = cellCache.cells[fs_cache::calculateNbrID(1,1,1)]->sysBoundaryFlag;
   cuint sysBoundaryLayer = cellCache.cells[fs_cache::calculateNbrID(1,1,1)]->sysBoundaryLayer;
   
   CellID leftNbrID,rghtNbrID;
   creal* left = NULL;
   creal* cent = cellCache.cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
   #ifdef DEBUG_SOLVERS
   if (cent[cp::RHOM] <= 0) {
      std::cerr << __FILE__ << ":" << __LINE__
         << (cent[cp::RHOM] < 0 ? " Negative" : " Zero") << " density in spatial cell " << cellID
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
       ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1))) {
      left = cellCache.cells[fs_cache::calculateNbrID(1-1,1  ,1  )]->parameters;
      #ifdef DEBUG_SOLVERS
      if (left[cp::RHOM] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (left[cp::RHOM] < 0 ? " Negative" : " Zero") << " density in spatial cell " << leftNbrID
            << std::endl;
         abort();
      }
      #endif
      rght = cellCache.cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->parameters;
      #ifdef DEBUG_SOLVERS
      if (rght[cp::RHOM] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (rght[cp::RHOM] < 0 ? " Negative" : " Zero") << " density in spatial cell " << rghtNbrID
            << std::endl;
         abort();
      }
      #endif
      
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         array[fs::drhomdx] = limiter(left[cp::RHOM],cent[cp::RHOM],rght[cp::RHOM]);
         array[fs::drhoqdx] = limiter(left[cp::RHOQ],cent[cp::RHOQ],rght[cp::RHOQ]);
         array[fs::dp11dx] = limiter(left[cp::P_11],cent[cp::P_11],rght[cp::P_11]);
         array[fs::dp22dx] = limiter(left[cp::P_22],cent[cp::P_22],rght[cp::P_22]);
         array[fs::dp33dx] = limiter(left[cp::P_33],cent[cp::P_33],rght[cp::P_33]);

         array[fs::dVxdx]  = limiter(left[cp::RHOMVX], left[cp::RHOM],
                                     cent[cp::RHOMVX], cent[cp::RHOM],
                                     rght[cp::RHOMVX], rght[cp::RHOM]);
         array[fs::dVydx]  = limiter(left[cp::RHOMVY], left[cp::RHOM],
                                     cent[cp::RHOMVY], cent[cp::RHOM],
                                     rght[cp::RHOMVY], rght[cp::RHOM]);
         array[fs::dVzdx]  = limiter(left[cp::RHOMVZ], left[cp::RHOM],
                                     cent[cp::RHOMVZ], cent[cp::RHOM],
                                     rght[cp::RHOMVZ], rght[cp::RHOM]);
         
         array[fs::dPERBydx]  = limiter(left[cp::PERBY],cent[cp::PERBY],rght[cp::PERBY]);
         array[fs::dPERBzdx]  = limiter(left[cp::PERBZ],cent[cp::PERBZ],rght[cp::PERBZ]);
         if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
            array[fs::dPERBydxx] = 0.0;
            array[fs::dPERBzdxx] = 0.0;
         } else {
            array[fs::dPERBydxx] = left[cp::PERBY] + rght[cp::PERBY] - 2.0*cent[cp::PERBY];
            array[fs::dPERBzdxx] = left[cp::PERBZ] + rght[cp::PERBZ] - 2.0*cent[cp::PERBZ];
         }
      }
      if (RKCase == RK_ORDER2_STEP1) {
         array[fs::drhomdx] = limiter(left[cp::RHOM_DT2],cent[cp::RHOM_DT2],rght[cp::RHOM_DT2]);
         array[fs::drhoqdx] = limiter(left[cp::RHOQ_DT2],cent[cp::RHOQ_DT2],rght[cp::RHOQ_DT2]);
         array[fs::dp11dx] = limiter(left[cp::P_11_DT2],cent[cp::P_11_DT2],rght[cp::P_11_DT2]);
         array[fs::dp22dx] = limiter(left[cp::P_22_DT2],cent[cp::P_22_DT2],rght[cp::P_22_DT2]);
         array[fs::dp33dx] = limiter(left[cp::P_33_DT2],cent[cp::P_33_DT2],rght[cp::P_33_DT2]);
         array[fs::dVxdx]  = limiter(left[cp::RHOMVX_DT2], left[cp::RHOM_DT2],
                                     cent[cp::RHOMVX_DT2], cent[cp::RHOM_DT2],
                                     rght[cp::RHOMVX_DT2], rght[cp::RHOM_DT2]);
         array[fs::dVydx]  = limiter(left[cp::RHOMVY_DT2], left[cp::RHOM_DT2],
                                     cent[cp::RHOMVY_DT2], cent[cp::RHOM_DT2],
                                     rght[cp::RHOMVY_DT2], rght[cp::RHOM_DT2]);
         array[fs::dVzdx]  = limiter(left[cp::RHOMVZ_DT2], left[cp::RHOM_DT2],
                                     cent[cp::RHOMVZ_DT2], cent[cp::RHOM_DT2],
                                     rght[cp::RHOMVZ_DT2], rght[cp::RHOM_DT2]);
         
         array[fs::dPERBydx]  = limiter(left[cp::PERBY_DT2],cent[cp::PERBY_DT2],rght[cp::PERBY_DT2]);
         array[fs::dPERBzdx]  = limiter(left[cp::PERBZ_DT2],cent[cp::PERBZ_DT2],rght[cp::PERBZ_DT2]);
         if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
            array[fs::dPERBydxx] = 0.0;
            array[fs::dPERBzdxx] = 0.0;
         } else {
            array[fs::dPERBydxx] = left[cp::PERBY_DT2] + rght[cp::PERBY_DT2] - 2.0*cent[cp::PERBY_DT2];
            array[fs::dPERBzdxx] = left[cp::PERBZ_DT2] + rght[cp::PERBZ_DT2] - 2.0*cent[cp::PERBZ_DT2];
         }
      }
   } else {
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 0);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)
            ->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 0);
      }
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
   if (((existingCells & CALCULATE_DY) == CALCULATE_DY) &&
       ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1))) {
      left = cellCache.cells[fs_cache::calculateNbrID(1  ,1-1,1  )]->parameters;
      rght = cellCache.cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->parameters;

      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         array[fs::drhomdy] = limiter(left[cp::RHOM],cent[cp::RHOM],rght[cp::RHOM]);
         array[fs::drhoqdy] = limiter(left[cp::RHOQ],cent[cp::RHOQ],rght[cp::RHOQ]);
         array[fs::dp11dy] = limiter(left[cp::P_11],cent[cp::P_11],rght[cp::P_11]);
         array[fs::dp22dy] = limiter(left[cp::P_22],cent[cp::P_22],rght[cp::P_22]);
         array[fs::dp33dy] = limiter(left[cp::P_33],cent[cp::P_33],rght[cp::P_33]);
         array[fs::dVxdy]  = limiter(left[cp::RHOMVX], left[cp::RHOM],
                                     cent[cp::RHOMVX], cent[cp::RHOM],
                                     rght[cp::RHOMVX], rght[cp::RHOM]);
         array[fs::dVydy]  = limiter(left[cp::RHOMVY], left[cp::RHOM],
                                     cent[cp::RHOMVY], cent[cp::RHOM],
                                     rght[cp::RHOMVY], rght[cp::RHOM]);
         array[fs::dVzdy]  = limiter(left[cp::RHOMVZ], left[cp::RHOM],
                                     cent[cp::RHOMVZ], cent[cp::RHOM],
                                     rght[cp::RHOMVZ], rght[cp::RHOM]);
         
         array[fs::dPERBxdy]  = limiter(left[cp::PERBX],cent[cp::PERBX],rght[cp::PERBX]);
         array[fs::dPERBzdy]  = limiter(left[cp::PERBZ],cent[cp::PERBZ],rght[cp::PERBZ]);

         if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
            array[fs::dPERBxdyy] = 0.0;
            array[fs::dPERBzdyy] = 0.0;
         } else {
            array[fs::dPERBxdyy] = left[cp::PERBX] + rght[cp::PERBX] - 2.0*cent[cp::PERBX];
            array[fs::dPERBzdyy] = left[cp::PERBZ] + rght[cp::PERBZ] - 2.0*cent[cp::PERBZ];
         }
      }
      if (RKCase == RK_ORDER2_STEP1) {
         array[fs::drhomdy] = limiter(left[cp::RHOM_DT2],cent[cp::RHOM_DT2],rght[cp::RHOM_DT2]);
         array[fs::drhoqdy] = limiter(left[cp::RHOQ_DT2],cent[cp::RHOQ_DT2],rght[cp::RHOQ_DT2]);
         array[fs::dp11dy] = limiter(left[cp::P_11_DT2],cent[cp::P_11_DT2],rght[cp::P_11_DT2]);
         array[fs::dp22dy] = limiter(left[cp::P_22_DT2],cent[cp::P_22_DT2],rght[cp::P_22_DT2]);
         array[fs::dp33dy] = limiter(left[cp::P_33_DT2],cent[cp::P_33_DT2],rght[cp::P_33_DT2]);
         array[fs::dVxdy]  = limiter(left[cp::RHOMVX_DT2], left[cp::RHOM_DT2],
                                     cent[cp::RHOMVX_DT2], cent[cp::RHOM_DT2],
                                     rght[cp::RHOMVX_DT2], rght[cp::RHOM_DT2]);
         array[fs::dVydy]  = limiter(left[cp::RHOMVY_DT2], left[cp::RHOM_DT2],
                                     cent[cp::RHOMVY_DT2], cent[cp::RHOM_DT2],
                                     rght[cp::RHOMVY_DT2], rght[cp::RHOM_DT2]);
         array[fs::dVzdy]  = limiter(left[cp::RHOMVZ_DT2], left[cp::RHOM_DT2],
                                     cent[cp::RHOMVZ_DT2], cent[cp::RHOM_DT2],
                                     rght[cp::RHOMVZ_DT2], rght[cp::RHOM_DT2]);
         
         array[fs::dPERBxdy]  = limiter(left[cp::PERBX_DT2],cent[cp::PERBX_DT2],rght[cp::PERBX_DT2]);
         array[fs::dPERBzdy]  = limiter(left[cp::PERBZ_DT2],cent[cp::PERBZ_DT2],rght[cp::PERBZ_DT2]);
         if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
            array[fs::dPERBxdyy] = 0.0;
            array[fs::dPERBzdyy] = 0.0;
         } else {
            array[fs::dPERBxdyy] = left[cp::PERBX_DT2] + rght[cp::PERBX_DT2] - 2.0*cent[cp::PERBX_DT2];
            array[fs::dPERBzdyy] = left[cp::PERBZ_DT2] + rght[cp::PERBZ_DT2] - 2.0*cent[cp::PERBZ_DT2];
         }
      }
   } else {
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 1);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if (((existingCells & CALCULATE_DZ) == CALCULATE_DZ) &&
       ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1))) {

      left = cellCache.cells[fs_cache::calculateNbrID(1  ,1  ,1-1)]->parameters;
      rght = cellCache.cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->parameters;

      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         array[fs::drhomdz] = limiter(left[cp::RHOM],cent[cp::RHOM],rght[cp::RHOM]);
         array[fs::drhoqdz] = limiter(left[cp::RHOQ],cent[cp::RHOQ],rght[cp::RHOQ]);
         array[fs::dp11dz] = limiter(left[cp::P_11],cent[cp::P_11],rght[cp::P_11]);
         array[fs::dp22dz] = limiter(left[cp::P_22],cent[cp::P_22],rght[cp::P_22]);
         array[fs::dp33dz] = limiter(left[cp::P_33],cent[cp::P_33],rght[cp::P_33]);
         array[fs::dVxdz]  = limiter(left[cp::RHOMVX], left[cp::RHOM],
                                     cent[cp::RHOMVX], cent[cp::RHOM],
                                     rght[cp::RHOMVX], rght[cp::RHOM]);
         array[fs::dVydz]  = limiter(left[cp::RHOMVY], left[cp::RHOM],
                                     cent[cp::RHOMVY], cent[cp::RHOM],
                                     rght[cp::RHOMVY], rght[cp::RHOM]);
         array[fs::dVzdz]  = limiter(left[cp::RHOMVZ], left[cp::RHOM],
                                     cent[cp::RHOMVZ], cent[cp::RHOM],
                                     rght[cp::RHOMVZ], rght[cp::RHOM]);
         
         array[fs::dPERBxdz]  = limiter(left[cp::PERBX],cent[cp::PERBX],rght[cp::PERBX]);
         array[fs::dPERBydz]  = limiter(left[cp::PERBY],cent[cp::PERBY],rght[cp::PERBY]);
         if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
            array[fs::dPERBxdzz] = 0.0;
            array[fs::dPERBydzz] = 0.0;
         } else {
            array[fs::dPERBxdzz] = left[cp::PERBX] + rght[cp::PERBX] - 2.0*cent[cp::PERBX];
            array[fs::dPERBydzz] = left[cp::PERBY] + rght[cp::PERBY] - 2.0*cent[cp::PERBY];
         }
      }
      if (RKCase == RK_ORDER2_STEP1) {
         array[fs::drhomdz] = limiter(left[cp::RHOM_DT2],cent[cp::RHOM_DT2],rght[cp::RHOM_DT2]);
         array[fs::drhoqdz] = limiter(left[cp::RHOQ_DT2],cent[cp::RHOQ_DT2],rght[cp::RHOQ_DT2]);
         array[fs::dp11dz] = limiter(left[cp::P_11_DT2],cent[cp::P_11_DT2],rght[cp::P_11_DT2]);
         array[fs::dp22dz] = limiter(left[cp::P_22_DT2],cent[cp::P_22_DT2],rght[cp::P_22_DT2]);
         array[fs::dp33dz] = limiter(left[cp::P_33_DT2],cent[cp::P_33_DT2],rght[cp::P_33_DT2]);
         array[fs::dVxdz]  = limiter(left[cp::RHOMVX_DT2], left[cp::RHOM_DT2],
                                     cent[cp::RHOMVX_DT2], cent[cp::RHOM_DT2],
                                     rght[cp::RHOMVX_DT2], rght[cp::RHOM_DT2]);
         array[fs::dVydz]  = limiter(left[cp::RHOMVY_DT2], left[cp::RHOM_DT2],
                                     cent[cp::RHOMVY_DT2], cent[cp::RHOM_DT2],
                                     rght[cp::RHOMVY_DT2], rght[cp::RHOM_DT2]);
         array[fs::dVzdz]  = limiter(left[cp::RHOMVZ_DT2], left[cp::RHOM_DT2],
                                     cent[cp::RHOMVZ_DT2], cent[cp::RHOM_DT2],
                                     rght[cp::RHOMVZ_DT2], rght[cp::RHOM_DT2]);
         
         array[fs::dPERBxdz]  = limiter(left[cp::PERBX_DT2],cent[cp::PERBX_DT2],rght[cp::PERBX_DT2]);
         array[fs::dPERBydz]  = limiter(left[cp::PERBY_DT2],cent[cp::PERBY_DT2],rght[cp::PERBY_DT2]);
         if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
            array[fs::dPERBxdzz] = 0.0;
            array[fs::dPERBydzz] = 0.0;
         } else {
            array[fs::dPERBxdzz] = left[cp::PERBX_DT2] + rght[cp::PERBX_DT2] - 2.0*cent[cp::PERBX_DT2];
            array[fs::dPERBydzz] = left[cp::PERBY_DT2] + rght[cp::PERBY_DT2] - 2.0*cent[cp::PERBY_DT2];
         }
      }
   } else {
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 2);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 2);
      }
   }
   
   if (Parameters::ohmHallTerm < 2 || sysBoundaryLayer == 1) {
      array[fs::dPERBxdyz] = 0.0;
      array[fs::dPERBydxz] = 0.0;
      array[fs::dPERBzdxy] = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      if (((existingCells & CALCULATE_DXY) == CALCULATE_DXY) &&
          ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1))) {
         botLeft = cellCache.cells[fs_cache::calculateNbrID(1-1,1-1,1  )]->parameters;
         botRght = cellCache.cells[fs_cache::calculateNbrID(1+1,1-1,1  )]->parameters;
         topLeft = cellCache.cells[fs_cache::calculateNbrID(1-1,1+1,1  )]->parameters;
         topRght = cellCache.cells[fs_cache::calculateNbrID(1+1,1+1,1  )]->parameters;

         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            array[fs::dPERBzdxy] = FOURTH * (botLeft[cp::PERBZ] + topRght[cp::PERBZ] - botRght[cp::PERBZ] - topLeft[cp::PERBZ]);
         }
         if (RKCase == RK_ORDER2_STEP1) {
            array[fs::dPERBzdxy] = FOURTH * (botLeft[cp::PERBZ_DT2] + topRght[cp::PERBZ_DT2] - botRght[cp::PERBZ_DT2] - topLeft[cp::PERBZ_DT2]);
         }
      } else {
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 3);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 3);
         }
      }
      
      // Calculate xz mixed derivatives:
      if (((existingCells & CALCULATE_DXZ) == CALCULATE_DXZ) &&
          ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1))) {

         botLeft = cellCache.cells[fs_cache::calculateNbrID(1-1,1  ,1-1)]->parameters;
         botRght = cellCache.cells[fs_cache::calculateNbrID(1+1,1  ,1-1)]->parameters;
         topLeft = cellCache.cells[fs_cache::calculateNbrID(1-1,1  ,1+1)]->parameters;
         topRght = cellCache.cells[fs_cache::calculateNbrID(1+1,1  ,1+1)]->parameters;

         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            array[fs::dPERBydxz] = FOURTH * (botLeft[cp::PERBY] + topRght[cp::PERBY] - botRght[cp::PERBY] - topLeft[cp::PERBY]);
         }
         if (RKCase == RK_ORDER2_STEP1) {
            array[fs::dPERBydxz] = FOURTH * (botLeft[cp::PERBY_DT2] + topRght[cp::PERBY_DT2] - botRght[cp::PERBY_DT2] - topLeft[cp::PERBY_DT2]);
         }
      } else {
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 4);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 4);
         }
      }
      
      // Calculate yz mixed derivatives:
      if (((existingCells & CALCULATE_DYZ) == CALCULATE_DYZ) &&
          ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1))) {

         botLeft = cellCache.cells[fs_cache::calculateNbrID(1  ,1-1,1-1)]->parameters;
         botRght = cellCache.cells[fs_cache::calculateNbrID(1  ,1+1,1-1)]->parameters;
         topLeft = cellCache.cells[fs_cache::calculateNbrID(1  ,1-1,1+1)]->parameters;
         topRght = cellCache.cells[fs_cache::calculateNbrID(1  ,1+1,1+1)]->parameters;

         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            array[fs::dPERBxdyz] = FOURTH * (botLeft[cp::PERBX] + topRght[cp::PERBX] - botRght[cp::PERBX] - topLeft[cp::PERBX]);
         }
         if (RKCase == RK_ORDER2_STEP1) {
            array[fs::dPERBxdyz] = FOURTH * (botLeft[cp::PERBX_DT2] + topRght[cp::PERBX_DT2] - botRght[cp::PERBX_DT2] - topLeft[cp::PERBX_DT2]);
         }
      } else {
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(mpiGrid, cellID, 5);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(mpiGrid, cellID, RKCase, 5);
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
   dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase,
   const bool communicateMoments) {
   int timer;
   namespace fs = fieldsolver;

   #warning IMPROVE ME this function is not multithreaded
   
   phiprof::start("Calculate face derivatives");

   switch (RKCase) {
    case RK_ORDER1:
      // Means initialising the solver as well as RK_ORDER1
      // standard case Exchange PERB* with neighbours
      // The update of PERB[XYZ] is needed after the system
      // boundary update of propagateMagneticFieldSimple.
      if(communicateMoments) {
        spatial_cell::SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB | Transfer::CELL_RHOM_V | Transfer::CELL_RHOQ | Transfer::CELL_P);
      } else {
        spatial_cell::SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB);
      }
      break;
    case RK_ORDER2_STEP1:
      // Exchange PERB*_DT2,RHO_DT2,RHOV*_DT2 with neighbours The
      // update of PERB[XYZ]_DT2 is needed after the system
      // boundary update of propagateMagneticFieldSimple.
      if(communicateMoments) {
        spatial_cell::SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERBDT2 | Transfer::CELL_RHOMDT2_VDT2 | Transfer::CELL_RHOQDT2 | Transfer::CELL_PDT2);
      } else {
        spatial_cell::SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERBDT2);
      }
      break;
    case RK_ORDER2_STEP2:
      // Exchange PERB*,RHO,RHOV* with neighbours The update of B
      // is needed after the system boundary update of
      // propagateMagneticFieldSimple.
      if(communicateMoments) {
        spatial_cell::SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB | Transfer::CELL_RHOM_V | Transfer::CELL_RHOQ | Transfer::CELL_P);
      } else {
        spatial_cell::SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB);
      }
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
   for (size_t c=0; c<fs_cache::getCache().cellsWithLocalNeighbours.size(); ++c) {
      const uint16_t localID = fs_cache::getCache().cellsWithLocalNeighbours[c];
      fs_cache::CellCache cache = fs_cache::getCache().localCellsCache[localID];

      if (cache.sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      calculateDerivatives(mpiGrid,cache,sysBoundaries, RKCase);
   }
   phiprof::stop(timer,fs_cache::getCache().cellsWithLocalNeighbours.size(),"Spatial Cells");

   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_receives(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);
   
   // Calculate derivatives on process boundary cells
   timer=phiprof::initializeTimer("Compute process boundary cells");
   phiprof::start(timer);
   for (size_t c=0; c<fs_cache::getCache().cellsWithRemoteNeighbours.size(); ++c) {
      const uint16_t localID = fs_cache::getCache().cellsWithRemoteNeighbours[c];
      fs_cache::CellCache cache = fs_cache::getCache().localCellsCache[localID];

      if (cache.sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      calculateDerivatives(mpiGrid,cache,sysBoundaries, RKCase);
   }
   phiprof::stop(timer,fs_cache::getCache().cellsWithRemoteNeighbours.size(),"Spatial Cells");

   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_sends();
   phiprof::stop(timer);

   const size_t N_cells = fs_cache::getCache().cellsWithLocalNeighbours.size()
     + fs_cache::getCache().cellsWithRemoteNeighbours.size();
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

void calculateBVOLDerivatives(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              std::vector<fs_cache::CellCache>& cache,
                              const std::vector<uint16_t>& cells,
                              SysBoundary& sysBoundaries) {
   #warning This function still contains mpiGrid!

   namespace cp = CellParams;
   namespace der = bvolderivatives;

   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const uint16_t localID = cells[c];
      if (cache[localID].sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;

      Real* const array = cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->derivativesBVOL;

      // Get boundary flag for the cell:
      cuint existingCells = cache[localID].existingCellsFlags;
      cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());

      creal* left = NULL;
      creal* cent = cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
      creal* rght = NULL;
      
      // Calculate x-derivatives (is not TVD for AMR mesh):
      if (((existingCells & CALCULATE_DX) == CALCULATE_DX) &&
          (cache[localID].sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)) {

         left = cache[localID].cells[fs_cache::calculateNbrID(1-1,1  ,1  )]->parameters;
         rght = cache[localID].cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->parameters;

         array[der::dPERBYVOLdx] = limiter(left[cp::PERBYVOL],cent[cp::PERBYVOL],rght[cp::PERBYVOL]);
         array[der::dPERBZVOLdx] = limiter(left[cp::PERBZVOL],cent[cp::PERBZVOL],rght[cp::PERBZVOL]);
      } else {
         if (cache[localID].sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(mpiGrid, cache[localID].cellID, 0);
         } else {
            sysBoundaries.getSysBoundary(cache[localID].sysBoundaryFlag)
              ->fieldSolverBoundaryCondBVOLDerivatives(mpiGrid, cache[localID].cellID, 0);
         }
      }
      
      // Calculate y-derivatives (is not TVD for AMR mesh):
      if (((existingCells & CALCULATE_DY) == CALCULATE_DY) &&
          (cache[localID].sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)) {
         
         left = cache[localID].cells[fs_cache::calculateNbrID(1  ,1-1,1  )]->parameters;
         rght = cache[localID].cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->parameters;

         array[der::dPERBXVOLdy] = limiter(left[cp::PERBXVOL],cent[cp::PERBXVOL],rght[cp::PERBXVOL]);
         array[der::dPERBZVOLdy] = limiter(left[cp::PERBZVOL],cent[cp::PERBZVOL],rght[cp::PERBZVOL]);
      } else {
         if (cache[localID].sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(mpiGrid, cache[localID].cellID, 1);
         } else {
            sysBoundaries.getSysBoundary(cache[localID].sysBoundaryFlag)
              ->fieldSolverBoundaryCondBVOLDerivatives(mpiGrid, cache[localID].cellID, 1);
         }
      }
      
      // Calculate z-derivatives (is not TVD for AMR mesh):
      if (((existingCells & CALCULATE_DZ) == CALCULATE_DZ) &&
          (cache[localID].sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)) {

         left = cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1-1)]->parameters;
         rght = cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->parameters;

         array[der::dPERBXVOLdz] = limiter(left[cp::PERBXVOL],cent[cp::PERBXVOL],rght[cp::PERBXVOL]);
         array[der::dPERBYVOLdz] = limiter(left[cp::PERBYVOL],cent[cp::PERBYVOL],rght[cp::PERBYVOL]);
      } else {
         if (cache[localID].sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(mpiGrid, cache[localID].cellID, 2);
         } else {
            sysBoundaries.getSysBoundary(cache[localID].sysBoundaryFlag)
              ->fieldSolverBoundaryCondBVOLDerivatives(mpiGrid, cache[localID].cellID, 2);
         }
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
   dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells
) {
   int timer;
   namespace fs = fieldsolver;
   
   phiprof::start("Calculate volume derivatives");
   
   spatial_cell::SpatialCell::set_mpi_transfer_type(Transfer::CELL_BVOL);
   mpiGrid.update_copies_of_remote_neighbors(FIELD_SOLVER_NEIGHBORHOOD_ID);
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);
   mpiGrid.start_remote_neighbor_copy_updates(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);

   // Calculate derivatives on process inner cells
   timer=phiprof::initializeTimer("Compute process inner cells");
   phiprof::start(timer);
   calculateBVOLDerivatives(mpiGrid,fs_cache::getCache().localCellsCache,fs_cache::getCache().cellsWithLocalNeighbours,sysBoundaries);
   phiprof::stop(timer,fs_cache::getCache().cellsWithLocalNeighbours.size(),"Spatial Cells");

   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_receives(FIELD_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop(timer);

   // Calculate derivatives on process boundary cells
   timer=phiprof::initializeTimer("Compute process boundary cells");
   phiprof::start(timer);
   calculateBVOLDerivatives(mpiGrid,fs_cache::getCache().localCellsCache,fs_cache::getCache().cellsWithRemoteNeighbours,sysBoundaries);
   phiprof::stop(timer,fs_cache::getCache().cellsWithRemoteNeighbours.size(),"Spatial Cells");

   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_remote_neighbor_copy_update_sends();
   phiprof::stop(timer);

   CellID N_cells = fs_cache::getCache().cellsWithLocalNeighbours.size() 
                  + fs_cache::getCache().cellsWithRemoteNeighbours.size();
   phiprof::stop("Calculate volume derivatives",N_cells,"Spatial Cells");
}
