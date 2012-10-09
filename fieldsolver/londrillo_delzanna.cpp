/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/*! \file londrillo_delzanna.cpp
 * \brief Londrillo -- Del Zanna upwind constrained transport field solver.
 * 
 * On the divergence-free condition in Godunov-type schemes for
 * ideal magnetohydrodynamics: the upwind constrained transport method,
 * P. Londrillo and L. Del Zanna, J. Comp. Phys., 195, 2004.
 * http://dx.doi.org/10.1016/j.jcp.2003.09.016
 *
 * Reconstructions taken from:
 * Efficient, high accuracy ADER-WENO schemes for hydrodynamics and
 * divergence-free magnetohydrodynamics, D. S. Balsara, T. Rumpf,
 * M. Dumbser, C.-D. Munz, J. Comp. Phys, 228, 2480-2516, 2009.
 * http://dx.doi.org/10.1016/j.jcp.2008.12.003
 * and
 * Divergence-free reconstruction of magnetic fields and WENO
 * schemes for magnetohydrodynamics, D. S. Balsara, J. Comp. Phys.,
 * 228, 5040-5056, 2009.
 * http://dx.doi.org/10.1016/j.jcp.2009.03.038
 * 
 * *****  NOTATION USED FOR VARIABLES FOLLOWS THE ONES USED  *****\n
 * *****      IN THE ABOVEMENTIONED PUBLICATION(S)           *****
 */

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <list>
#include <set>

#include "../common.h"
#include "../parameters.h"
#include "../fieldsolver.h"
#include "limiters.h"
#include "../project.h"
#include "phiprof.hpp"

using namespace std;
using namespace fieldsolver;

#include <stdint.h>
typedef uint64_t CellID;

static creal EPS = 1.0e-30;

// TODO why this?
// static set<CellID> ghostCells;

static map<CellID,uint> sysBoundaryFlags;
/*!< Boundary status flags for all cells on this process. Here "boundary cell" 
 * means that the cell is at the physical boundary of the simulation volume, 
 * in some cases this condition means that this cell is a "ghost cell". However, 
 * this is algorithm-dependent, so one must be careful with ghost cell definition.
 * 
 * Consider a cell and its immediate neighbours (26 in total), i.e. a 3x3 cube 
 * of cells at base grid level. Considered cell is at the center of the cube. 
 * Number the cells with the usual C array numbering, k*9+j*3+i, where i,j,k 
 * are the cell indices in x,y,z directions. Each existing cell within the 
 * 3x3 cube has its bit (calculated with C array indexing) flipped to value 1.
 * The bit 13 is always set to unit value (considered cell always exists).
 * 
 * These boundary flags can be used to determine whether a numerical algorithm 
 * should be applied to a cell, for example, to calculate an edge electric field.
 * The boundary status can be checked with a single bitwise operation instead of 
 * N if-statements.
 * 
 * Note that this definition works with mesh refinement. The boundary flag 
 * should only change for a cell if some of its neighbours are deleted or 
 * created during the simulation.
 */


static uint CALCULATE_DX = 0; /**< Bit mask determining if x-derivatives can be calculated on a cell.*/
static uint CALCULATE_DY = 0; /**< Bit mask determining if y-derivatives can be calculated on a cell.*/
static uint CALCULATE_DZ = 0; /**< Bit mask determining if z-derivatives can be calculated on a cell.*/
static uint CALCULATE_EX = 0; /**< Bit mask determining if edge Ex can be calculated on a cell.*/
static uint CALCULATE_EY = 0; /**< Bit mask determining if edge Ey can be calculated on a cell.*/
static uint CALCULATE_EZ = 0; /**< Bit mask determining if edge Ez can be calculated on a cell.*/
static uint PROPAGATE_BX = 0; /**< Bit mask determining if face Bx is propagated on a cell.*/
static uint PROPAGATE_BY = 0; /**< Bit mask determining if face By is propagated on a cell.*/
static uint PROPAGATE_BZ = 0; /**< Bit mask determining if face Bz is propagated on a cell.*/

// Constants: not needed as such, but if field solver is implemented on GPUs 
// these force CPU to use float accuracy, which in turn helps to compare 
// CPU and GPU results.

const Real HALF    = 0.5;
const Real MINUS   = -1.0;
const Real PLUS    = +1.0;
const Real EIGTH   = 1.0/8.0;
const Real FOURTH  = 1.0/4.0;
const Real SIXTH   = 1.0/6.0;
const Real TWELWTH = 1.0/12.0;
const Real TWO     = 2.0;
const Real ZERO    = 0.0;

/*! Steps in the Runge-Kutta methods */
enum {RK_ORDER1,     	/*!< First order method, one step (and initialisation) */
      RK_ORDER2_STEP1,	/*!< Two-step second order method, first step */
      RK_ORDER2_STEP2	/*!< Two-step second order method, second step */
};

void calculateDerivativesSimple(dccrg::Dccrg<SpatialCell>& mpiGrid,const vector<CellID>& localCells, cint& RKCase);
void calculateUpwindedElectricFieldSimple(dccrg::Dccrg<SpatialCell>& mpiGrid,const vector<CellID>& localCells, cint& RKCase);


/*! \brief Calculate the neighbour number.
 * 
 * Calculate the neighbour number. For the inspected cell the (i,j,k) are (1,1,1). Add or 
 * reduce one from an index to get the "neighbour number" for the neighbour in that direction. 
 * For example, neighbour number for i-1,j-1,k neighbour is calculated with calcNbrNumber(1-1,1-1,1+0).
 * Thus, the cell in question has a "neighbour number" 13.
 * The purpose of this function (and neighbour numbers) is to indicate whether a cell has 
 * existing neighbours on a given direction. The neighbour existence status can be stored in 
 * a single 32bit word and tested with bitwise operations.
 */
inline uchar calcNbrNumber(const uchar& i,const uchar& j,const uchar& k) {return k*9+j*3+i;}

inline uchar calcNbrTypeID(const uchar& i,const uchar& j,const uchar& k) {return k*25+j*5+i;}

/*! Select the limiter to be used in the field solver. Limiters defined in limiters.h */
Real limiter(creal& left,creal& cent,creal& rght) {
   //const Real limited = minmod(left,cent,rght);
   const Real limited = MClimiter(left,cent,rght);
   //const Real limited = vanLeer(left,cent,rght);

   #ifdef DEBUG_SOLVERS
   if (limited != limited
   || limited * 0 != 0) {
      std::cerr << __FILE__ << ":" << __LINE__
         << " Limiter returned an invalid value " << limited
         << " with left, center, right: " << left << ", " << cent << ", " << rght
         << std::endl;
      abort();
   }
   #endif

   return limited;
}

CellID getNeighbourID(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   const CellID& cellID,
   const uchar& i,
   const uchar& j,
   const uchar& k
) {
   const std::vector<CellID> neighbors = mpiGrid.get_neighbors_of(cellID, int(i) - 2, int(j) - 2, int(k) - 2);
   if (neighbors.size() == 0) {
      cerr << __FILE__ << ":" << __LINE__
         << " No neighbor for cell " << cellID
         << " at offsets " << int(i) - 2 << ", " << int(j) - 2 << ", " << int(k) - 2
         << endl;
      abort();
   }
   // TODO support spatial refinement
   return neighbors[0];
}

static void calculateSysBoundaryFlags(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   const vector<CellID>& localCells
) {
   sysBoundaryFlags.clear();
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      
      // Raise the bit for each existing cell within a 3x3 cube of 
      // spatial cells. This cell sits at the center of the cube.
      uint sysBoundaryFlag = (1 << calcNbrNumber(1,1,1)); // The cell itself exists (bit 13 set to 1)
      
      for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
         if (i == 0 && (j == 0 && k == 0)) continue;
         const CellID nbr = getNeighbourID(mpiGrid, cellID, 2 + i, 2 + j, 2 + k);
         if (nbr == INVALID_CELLID) continue;
         sysBoundaryFlag = sysBoundaryFlag | (1 << calcNbrNumber(1+i,1+j,1+k));
      }
      sysBoundaryFlags[cellID] = sysBoundaryFlag;
   }
}

/*! \brief Low-level helper function.
 * 
 * Avoid crashes on zero density by returning V = rhoV = 0.0 if rho = 0.0.
 * 
 */
Real divideIfNonZero(creal rhoV, creal rho) {
   if(rho <= 0.0) {
      return 0.0;
   } else {
      return rhoV / rho;
   }
}

/*! \brief Low-level spatial derivatives calculation.
 * 
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in project.h. Uses RHO, RHOV[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the second-order method, and RHO1, RHOV[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 * \param cellID Index of the cell to process
 * \param mpiGrid Grid
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
static void calculateDerivatives(
   const CellID& cellID,
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   cint& RKCase
) {
   namespace cp = CellParams;
   namespace fs = fieldsolver;
   Real* const array       = mpiGrid[cellID]->derivatives;
   Real* const derivatives = array;
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
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if ((existingCells & CALCULATE_DX) == CALCULATE_DX) {
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
         array[fs::drhodx] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
         array[fs::dVxdx]  = limiter(divideIfNonZero(left[cp::RHOVX], left[cp::RHO]),
                                     divideIfNonZero(cent[cp::RHOVX], cent[cp::RHO]),
                                     divideIfNonZero(rght[cp::RHOVX], rght[cp::RHO]));
         array[fs::dVydx]  = limiter(divideIfNonZero(left[cp::RHOVY], left[cp::RHO]),
                                     divideIfNonZero(cent[cp::RHOVY], cent[cp::RHO]),
                                     divideIfNonZero(rght[cp::RHOVY], rght[cp::RHO]));
         array[fs::dVzdx]  = limiter(divideIfNonZero(left[cp::RHOVZ], left[cp::RHO]),
                                     divideIfNonZero(cent[cp::RHOVZ], cent[cp::RHO]),
                                     divideIfNonZero(rght[cp::RHOVZ], rght[cp::RHO]));
         array[fs::dBydx]  = limiter(left[cp::BGBY]+left[cp::PERBY],cent[cp::BGBY]+cent[cp::PERBY],rght[cp::BGBY]+rght[cp::PERBY]);
         array[fs::dBzdx]  = limiter(left[cp::BGBZ]+left[cp::PERBZ],cent[cp::BGBZ]+cent[cp::PERBZ],rght[cp::BGBZ]+rght[cp::PERBZ]);
      }
      if (RKCase == RK_ORDER2_STEP1) {
         array[fs::drhodx] = limiter(left[cp::RHO1],cent[cp::RHO1],rght[cp::RHO1]);
         array[fs::dVxdx]  = limiter(divideIfNonZero(left[cp::RHOVX1], left[cp::RHO1]),
                                     divideIfNonZero(cent[cp::RHOVX1], cent[cp::RHO1]),
                                     divideIfNonZero(rght[cp::RHOVX1], rght[cp::RHO1]));
         array[fs::dVydx]  = limiter(divideIfNonZero(left[cp::RHOVY1], left[cp::RHO1]),
                                     divideIfNonZero(cent[cp::RHOVY1], cent[cp::RHO1]),
                                     divideIfNonZero(rght[cp::RHOVY1], rght[cp::RHO1]));
         array[fs::dVzdx]  = limiter(divideIfNonZero(left[cp::RHOVZ1], left[cp::RHO1]),
                                     divideIfNonZero(cent[cp::RHOVZ1], cent[cp::RHO1]),
                                     divideIfNonZero(rght[cp::RHOVZ1], rght[cp::RHO1]));
         array[fs::dBydx]  = limiter(left[cp::BGBY]+left[cp::PERBY1],cent[cp::BGBY]+cent[cp::PERBY1],rght[cp::BGBY]+rght[cp::PERBY1]);
         array[fs::dBzdx]  = limiter(left[cp::BGBZ]+left[cp::PERBZ1],cent[cp::BGBZ]+cent[cp::PERBZ1],rght[cp::BGBZ]+rght[cp::PERBZ1]);
      }
   } else {
      fieldSolverBoundaryCondDerivX(cellID,array,existingCells,nonExistingCells,derivatives,mpiGrid);
   }
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
   if ((existingCells & CALCULATE_DY) == CALCULATE_DY) {
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
         array[fs::drhody] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
         array[fs::dVxdy]  = limiter(divideIfNonZero(left[cp::RHOVX], left[cp::RHO]),
                                     divideIfNonZero(cent[cp::RHOVX], cent[cp::RHO]),
                                     divideIfNonZero(rght[cp::RHOVX], rght[cp::RHO]));
         array[fs::dVydy]  = limiter(divideIfNonZero(left[cp::RHOVY], left[cp::RHO]),
                                     divideIfNonZero(cent[cp::RHOVY], cent[cp::RHO]),
                                     divideIfNonZero(rght[cp::RHOVY], rght[cp::RHO]));
         array[fs::dVzdy]  = limiter(divideIfNonZero(left[cp::RHOVZ], left[cp::RHO]),
                                     divideIfNonZero(cent[cp::RHOVZ], cent[cp::RHO]),
                                     divideIfNonZero(rght[cp::RHOVZ], rght[cp::RHO]));
         array[fs::dBxdy]  = limiter(left[cp::BGBX]+left[cp::PERBX],cent[cp::BGBX]+cent[cp::PERBX],rght[cp::BGBX]+rght[cp::PERBX]);
         array[fs::dBzdy]  = limiter(left[cp::BGBZ]+left[cp::PERBZ],cent[cp::BGBZ]+cent[cp::PERBZ],rght[cp::BGBZ]+rght[cp::PERBZ]);
      }
      if (RKCase == RK_ORDER2_STEP1) {
         array[fs::drhody] = limiter(left[cp::RHO1],cent[cp::RHO1],rght[cp::RHO1]);
         array[fs::dVxdy]  = limiter(divideIfNonZero(left[cp::RHOVX1], left[cp::RHO1]),
                                     divideIfNonZero(cent[cp::RHOVX1], cent[cp::RHO1]),
                                     divideIfNonZero(rght[cp::RHOVX1], rght[cp::RHO1]));
         array[fs::dVydy]  = limiter(divideIfNonZero(left[cp::RHOVY1], left[cp::RHO1]),
                                     divideIfNonZero(cent[cp::RHOVY1], cent[cp::RHO1]),
                                     divideIfNonZero(rght[cp::RHOVY1], rght[cp::RHO1]));
         array[fs::dVzdy]  = limiter(divideIfNonZero(left[cp::RHOVZ1], left[cp::RHO1]),
                                     divideIfNonZero(cent[cp::RHOVZ1], cent[cp::RHO1]),
                                     divideIfNonZero(rght[cp::RHOVZ1], rght[cp::RHO1]));
         array[fs::dBxdy]  = limiter(left[cp::BGBX]+left[cp::PERBX1],cent[cp::BGBX]+cent[cp::PERBX1],rght[cp::BGBX]+rght[cp::PERBX1]);
         array[fs::dBzdy]  = limiter(left[cp::BGBZ]+left[cp::PERBZ1],cent[cp::BGBZ]+cent[cp::PERBZ1],rght[cp::BGBZ]+rght[cp::PERBZ1]);
      }
   } else {
      fieldSolverBoundaryCondDerivY(cellID,array,existingCells,nonExistingCells,derivatives,mpiGrid);
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if ((existingCells & CALCULATE_DZ) == CALCULATE_DZ) {
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
         array[fs::drhodz] = limiter(left[cp::RHO],cent[cp::RHO],rght[cp::RHO]);
         array[fs::dVxdz]  = limiter(divideIfNonZero(left[cp::RHOVX], left[cp::RHO]),
                                     divideIfNonZero(cent[cp::RHOVX], cent[cp::RHO]),
                                     divideIfNonZero(rght[cp::RHOVX], rght[cp::RHO]));
         array[fs::dVydz]  = limiter(divideIfNonZero(left[cp::RHOVY], left[cp::RHO]),
                                     divideIfNonZero(cent[cp::RHOVY], cent[cp::RHO]),
                                     divideIfNonZero(rght[cp::RHOVY], rght[cp::RHO]));
         array[fs::dVzdz]  = limiter(divideIfNonZero(left[cp::RHOVZ], left[cp::RHO]),
                                     divideIfNonZero(cent[cp::RHOVZ], cent[cp::RHO]),
                                     divideIfNonZero(rght[cp::RHOVZ], rght[cp::RHO]));
         array[fs::dBxdz]  = limiter(left[cp::BGBX]+left[cp::PERBX],cent[cp::BGBX]+cent[cp::PERBX],rght[cp::BGBX]+rght[cp::PERBX]);
         array[fs::dBydz]  = limiter(left[cp::BGBY]+left[cp::PERBY],cent[cp::BGBY]+cent[cp::PERBY],rght[cp::BGBY]+rght[cp::PERBY]);
      }
      if (RKCase == RK_ORDER2_STEP1) {
         array[fs::drhodz] = limiter(left[cp::RHO1],cent[cp::RHO1],rght[cp::RHO1]);
         array[fs::dVxdz]  = limiter(divideIfNonZero(left[cp::RHOVX1], left[cp::RHO1]),
                                     divideIfNonZero(cent[cp::RHOVX1], cent[cp::RHO1]),
                                     divideIfNonZero(rght[cp::RHOVX1], rght[cp::RHO1]));
         array[fs::dVydz]  = limiter(divideIfNonZero(left[cp::RHOVY1], left[cp::RHO1]),
                                     divideIfNonZero(cent[cp::RHOVY1], cent[cp::RHO1]),
                                     divideIfNonZero(rght[cp::RHOVY1], rght[cp::RHO1]));
         array[fs::dVzdz]  = limiter(divideIfNonZero(left[cp::RHOVZ1], left[cp::RHO1]),
                                     divideIfNonZero(cent[cp::RHOVZ1], cent[cp::RHO1]),
                                     divideIfNonZero(rght[cp::RHOVZ1], rght[cp::RHO1]));
         array[fs::dBxdz]  = limiter(left[cp::BGBX]+left[cp::PERBX1],cent[cp::BGBX]+cent[cp::PERBX1],rght[cp::BGBX]+rght[cp::PERBX1]);
         array[fs::dBydz]  = limiter(left[cp::BGBY]+left[cp::PERBY1],cent[cp::BGBY]+cent[cp::PERBY1],rght[cp::BGBY]+rght[cp::PERBY1]);
      }
   } else {
      fieldSolverBoundaryCondDerivZ(cellID,array,existingCells,nonExistingCells,derivatives,mpiGrid);
   }
}

/*! \brief Low-level helper function.
 * 
 * Computes the magnetosonic speed in the YZ plane. Used in upwinding the electric field X component.
 * 
 * Selects the RHO/RHO1 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * 
 * If fields are not propagated, rewturns 0.0 as there is no information propagating.
 * 
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
template<typename REAL> REAL calculateFastMSspeedYZ(const REAL* cp, const REAL* derivs, const REAL* nbr_cp, const REAL* nbr_derivs, const REAL& By, const REAL& Bz, const REAL& dBydx, const REAL& dBydz, const REAL& dBzdx, const REAL& dBzdy, const REAL& ydir, const REAL& zdir, cint& RKCase
) {
   namespace fs = fieldsolver;
   namespace pc = physicalconstants;
   REAL A_0, A_X, rho;
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      A_0  = HALF*(nbr_cp[CellParams::PERBX] + nbr_cp[CellParams::BGBX] + cp[CellParams::PERBX] + cp[CellParams::BGBX]);
      A_X  = (nbr_cp[CellParams::PERBX] + nbr_cp[CellParams::BGBX]) - (cp[CellParams::PERBX] + cp[CellParams::BGBX]);
      rho = Parameters::m*(cp[CellParams::RHO] + ydir*HALF*derivs[fs::drhody] + zdir*HALF*derivs[fs::drhodz]);
   } else { // RKCase == RK_ORDER2_STEP1
      A_0  = HALF*(nbr_cp[CellParams::PERBX1] + nbr_cp[CellParams::BGBX] + cp[CellParams::PERBX1] + cp[CellParams::BGBX]);
      A_X  = (nbr_cp[CellParams::PERBX1] + nbr_cp[CellParams::BGBX]) - (cp[CellParams::PERBX1] + cp[CellParams::BGBX]);
      rho = Parameters::m*(cp[CellParams::RHO1] + ydir*HALF*derivs[fs::drhody] + zdir*HALF*derivs[fs::drhodz]);
   }
   const REAL A_Y  = nbr_derivs[fs::dBxdy]  + derivs[fs::dBxdy];
   const REAL A_XY = nbr_derivs[fs::dBxdy]  - derivs[fs::dBxdy];
   const REAL A_Z  = nbr_derivs[fs::dBxdz]  + derivs[fs::dBxdz];
   const REAL A_XZ = nbr_derivs[fs::dBxdz]  - derivs[fs::dBxdz];
   
   const REAL Bx2  = (A_0 + ydir*HALF*A_Y + zdir*HALF*A_Z)*(A_0 + ydir*HALF*A_Y + zdir*HALF*A_Z)
     + TWELWTH*(A_X + ydir*HALF*A_XY + zdir*HALF*A_XZ)*(A_X + ydir*HALF*A_XY + zdir*HALF*A_XZ); // OK
   const REAL By2  = (By + zdir*HALF*dBydz)*(By + zdir*HALF*dBydz) + TWELWTH*dBydx*dBydx; // OK
   const REAL Bz2  = (Bz + ydir*HALF*dBzdy)*(Bz + ydir*HALF*dBzdy) + TWELWTH*dBzdx*dBzdx; // OK
   
   if(!Parameters::propagateField) {
      return 0.0;
   } else {
      return sqrt((Bx2+By2+Bz2) / (pc::MU_0 * rho));
   }
}

/*! \brief Low-level helper function.
 * 
 * Computes the magnetosonic speed in the XZ plane. Used in upwinding the electric field Y component.
 * 
 * Selects the RHO/RHO1 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * 
 * If fields are not propagated, rewturns 0.0 as there is no information propagating.
 * 
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
template<typename REAL> REAL calculateFastMSspeedXZ(const REAL* cp, const REAL* derivs, const REAL* nbr_cp, const REAL* nbr_derivs, const REAL& Bx, const REAL& Bz, const REAL& dBxdy, const REAL& dBxdz, const REAL& dBzdx, const REAL& dBzdy, const REAL& xdir,const REAL& zdir, cint& RKCase
) {
   namespace fs = fieldsolver;
   namespace pc = physicalconstants;
   REAL B_0, B_Y, rho;
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      B_0  = HALF*(nbr_cp[CellParams::PERBY] + nbr_cp[CellParams::BGBY] + cp[CellParams::PERBY] + cp[CellParams::BGBY]);
      B_Y  = (nbr_cp[CellParams::PERBY] + nbr_cp[CellParams::BGBY]) - (cp[CellParams::PERBY] + cp[CellParams::BGBY]);
      rho = Parameters::m*(cp[CellParams::RHO] + xdir*HALF*derivs[fs::drhodx] + zdir*HALF*derivs[fs::drhodz]);
   } else { // RKCase == RK_ORDER2_STEP1
      B_0  = HALF*(nbr_cp[CellParams::PERBY1] + nbr_cp[CellParams::BGBY] + cp[CellParams::PERBY1] + cp[CellParams::BGBY]);
      B_Y  = (nbr_cp[CellParams::PERBY1] + nbr_cp[CellParams::BGBY]) - (cp[CellParams::PERBY1] + cp[CellParams::BGBY]);
      rho = Parameters::m*(cp[CellParams::RHO1] + xdir*HALF*derivs[fs::drhodx] + zdir*HALF*derivs[fs::drhodz]);
   }
   const REAL B_X  = nbr_derivs[fs::dBydx]  + derivs[fs::dBydx];
   const REAL B_XY = nbr_derivs[fs::dBydx]  - derivs[fs::dBydx];
   const REAL B_Z  = nbr_derivs[fs::dBydz]  + derivs[fs::dBydz];
   const REAL B_YZ = nbr_derivs[fs::dBydz]  - derivs[fs::dBydz];
   
   const REAL By2  = (B_0 + xdir*HALF*B_X + zdir*HALF*B_Z)*(B_0 + xdir*HALF*B_X + zdir*HALF*B_Z)
     + TWELWTH*(B_Y + xdir*HALF*B_XY + zdir*HALF*B_YZ)*(B_Y + xdir*HALF*B_XY + zdir*HALF*B_YZ); // OK
   const REAL Bx2  = (Bx + zdir*HALF*dBxdz)*(Bx + zdir*HALF*dBxdz) + TWELWTH*dBxdy*dBxdy; // OK
   const REAL Bz2  = (Bz + xdir*HALF*dBzdx)*(Bz + xdir*HALF*dBzdx) + TWELWTH*dBzdy*dBzdy; // OK
   
   if(!Parameters::propagateField) {
      return 0.0;
   } else {
      return sqrt((Bx2+By2+Bz2) / (pc::MU_0 * rho));
   }
}

/*! \brief Low-level helper function.
 * 
 * Computes the magnetosonic speed in the XY plane. Used in upwinding the electric field Z component.
 * 
 * Selects the RHO/RHO1 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * 
 * If fields are not propagated, rewturns 0.0 as there is no information propagating.
 * 
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
template<typename REAL> REAL calculateFastMSspeedXY(const REAL* cp, const REAL* derivs, const REAL* nbr_cp, const REAL* nbr_derivs, const REAL& Bx, const REAL& By, const REAL& dBxdy, const REAL& dBxdz, const REAL& dBydx, const REAL& dBydz, const REAL& xdir,const REAL& ydir, cint& RKCase
) {
   namespace fs = fieldsolver;
   namespace pc = physicalconstants;
   REAL C_0, C_Z, rho;
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      C_0  = HALF*(nbr_cp[CellParams::PERBZ] + nbr_cp[CellParams::BGBZ] + cp[CellParams::PERBZ] + cp[CellParams::BGBZ]);
      C_Z  = (nbr_cp[CellParams::PERBZ] + nbr_cp[CellParams::BGBZ]) - (cp[CellParams::PERBZ] + cp[CellParams::BGBZ]);
      rho = Parameters::m*(cp[CellParams::RHO] + xdir*HALF*derivs[fs::drhodx] + ydir*HALF*derivs[fs::drhody]);
   } else { // RKCase == RK_ORDER2_STEP1
      C_0  = HALF*(nbr_cp[CellParams::PERBZ1] + nbr_cp[CellParams::BGBZ] + cp[CellParams::PERBZ1] + cp[CellParams::BGBZ]);
      C_Z  = (nbr_cp[CellParams::PERBZ1] + nbr_cp[CellParams::BGBZ]) - (cp[CellParams::PERBZ1] + cp[CellParams::BGBZ]);
      rho = Parameters::m*(cp[CellParams::RHO1] + xdir*HALF*derivs[fs::drhodx] + ydir*HALF*derivs[fs::drhody]);
   }
   const REAL C_X  = nbr_derivs[fs::dBzdx]  + derivs[fs::dBzdx];
   const REAL C_XZ = nbr_derivs[fs::dBzdx]  - derivs[fs::dBzdx];
   const REAL C_Y  = nbr_derivs[fs::dBzdy]  + derivs[fs::dBzdy];
   const REAL C_YZ = nbr_derivs[fs::dBzdy]  - derivs[fs::dBzdy];
   
   const REAL Bz2  = (C_0 + xdir*HALF*C_X + ydir*HALF*C_Y)*(C_0 + xdir*HALF*C_X + ydir*HALF*C_Y)
     + TWELWTH*(C_Z + xdir*HALF*C_XZ + ydir*HALF*C_YZ)*(C_Z + xdir*HALF*C_XZ + ydir*HALF*C_YZ);
   const REAL Bx2  = (Bx + ydir*HALF*dBxdy)*(Bx + ydir*HALF*dBxdy) + TWELWTH*dBxdz*dBxdz;
   const REAL By2  = (By + xdir*HALF*dBydx)*(By + xdir*HALF*dBydx) + TWELWTH*dBydz*dBydz;
   
   if(!Parameters::propagateField) {
      return 0.0;
   } else {
      return sqrt((Bx2+By2+Bz2) / (pc::MU_0 * rho));
   }
}

/*! \brief Low-level electric field propagation function.
 * 
 * Computes the upwinded electric field X component along the cell's corresponding edge as the cross product of B and V in the YZ plane. Also includes the calculation of the maximally allowed time step.
 * Selects the RHO/RHO1 RHOV[XYZ]/RHOV[XYZ]1 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * \param cellID Index of the cell to process
 * \param mpiGrid Grid
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
static void calculateEdgeElectricFieldX(
      const CellID& cellID,
      dccrg::Dccrg<SpatialCell>& mpiGrid,
      cint& RKCase
) {
   namespace fs = fieldsolver;
   
   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge.   
   Real ay_pos,ay_neg;              // Max. characteristic velocities to y-direction
   Real az_pos,az_neg;              // Max. characteristic velocities to z-direction
   Real Vy0,Vz0;                    // Reconstructed V
   Real c_y, c_z;                   // Wave speeds to yz-directions
   
   // Get read-only pointers to NE,NW,SE,SW states (SW is rw, result is written there):
   const CellID nbr_SE = getNeighbourID(mpiGrid, cellID, 2  , 2-1, 2  );
   const CellID nbr_NE = getNeighbourID(mpiGrid, cellID, 2  , 2-1, 2-1);
   const CellID nbr_NW = getNeighbourID(mpiGrid, cellID, 2  , 2  , 2-1);
   #ifndef NDEBUG
      if (nbr_SE == INVALID_CELLID) {cerr << "ERROR: Could not find SE neighbour!" << endl; exit(1);}
      if (nbr_NE == INVALID_CELLID) {cerr << "ERROR: Could not find NE neighbour!" << endl; exit(1);}
      if (nbr_NW == INVALID_CELLID) {cerr << "ERROR: Could not find NW neighbour!" << endl; exit(1);}
   #endif
   
   Real*  const cp_SW = mpiGrid[cellID]->parameters;
   creal* const cp_SE = mpiGrid[nbr_SE]->parameters;
   creal* const cp_NE = mpiGrid[nbr_NE]->parameters;
   creal* const cp_NW = mpiGrid[nbr_NW]->parameters;
   
   creal* const derivs_SW = mpiGrid[cellID]->derivatives;
   creal* const derivs_SE = mpiGrid[nbr_SE]->derivatives;
   creal* const derivs_NE = mpiGrid[nbr_NE]->derivatives;
   creal* const derivs_NW = mpiGrid[nbr_NW]->derivatives;
   
   Real By_S, Bz_W, Bz_E, By_N;
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      By_S = cp_SW[CellParams::PERBY]+cp_SW[CellParams::BGBY];
      Bz_W = cp_SW[CellParams::PERBZ]+cp_SW[CellParams::BGBZ];
      Bz_E = cp_SE[CellParams::PERBZ]+cp_SE[CellParams::BGBZ];
      By_N = cp_NW[CellParams::PERBY]+cp_NW[CellParams::BGBY];
      Vy0  = divideIfNonZero(cp_SW[CellParams::RHOVY], cp_SW[CellParams::RHO]);
      Vz0  = divideIfNonZero(cp_SW[CellParams::RHOVZ], cp_SW[CellParams::RHO]);
   } else { // RKCase == RK_ORDER2_STEP1
      By_S = cp_SW[CellParams::PERBY1]+cp_SW[CellParams::BGBY];
      Bz_W = cp_SW[CellParams::PERBZ1]+cp_SW[CellParams::BGBZ];
      Bz_E = cp_SE[CellParams::PERBZ1]+cp_SE[CellParams::BGBZ];
      By_N = cp_NW[CellParams::PERBY1]+cp_NW[CellParams::BGBY];
      Vy0  = divideIfNonZero(cp_SW[CellParams::RHOVY1], cp_SW[CellParams::RHO1]);
      Vz0  = divideIfNonZero(cp_SW[CellParams::RHOVZ1], cp_SW[CellParams::RHO1]);
   }
   
   creal dBydx_S = derivs_SW[fs::dBydx];
   creal dBydz_S = derivs_SW[fs::dBydz];
   creal dBzdx_W = derivs_SW[fs::dBzdx];
   creal dBzdy_W = derivs_SW[fs::dBzdy];
   creal dBzdx_E = derivs_SE[fs::dBzdx];
   creal dBzdy_E = derivs_SE[fs::dBzdy];
   creal dBydx_N = derivs_NW[fs::dBydx];
   creal dBydz_N = derivs_NW[fs::dBydz];
   
   // Ex and characteristic speeds on this cell:
   // 1st order terms:
   Real Ex_SW = By_S*Vz0 - Bz_W*Vy0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_SW += +HALF*((By_S - HALF*dBydz_S)*(-derivs_SW[fs::dVzdy] - derivs_SW[fs::dVzdz]) - dBydz_S*Vz0 + SIXTH*dBydx_S*derivs_SW[fs::dVzdx]);
      Ex_SW += -HALF*((Bz_W - HALF*dBzdy_W)*(-derivs_SW[fs::dVydy] - derivs_SW[fs::dVydz]) - dBzdy_W*Vy0 + SIXTH*dBzdx_W*derivs_SW[fs::dVydx]);
   #endif
   
   const CellID nbrID_SW      = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2  );
   #ifndef NDEBUG
      if (nbrID_SW == INVALID_CELLID) {cerr << "ERROR: Could not find SW cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_SW     = mpiGrid[nbrID_SW]->parameters;
   creal* const nbr_derivs_SW = mpiGrid[nbrID_SW]->derivatives;
   c_y = calculateFastMSspeedYZ(cp_SW,derivs_SW,nbr_cp_SW,nbr_derivs_SW,By_S,Bz_W,dBydx_S,dBydz_S,dBzdx_W,dBzdy_W,MINUS,MINUS, RKCase);
   c_z = c_y;
   ay_neg   = max(ZERO,-Vy0 + c_y);
   ay_pos   = max(ZERO,+Vy0 + c_y);
   az_neg   = max(ZERO,-Vz0 + c_z);
   az_pos   = max(ZERO,+Vz0 + c_z);
   
   // Ex and characteristic speeds on j-1 neighbour:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vy0  = divideIfNonZero(cp_SE[CellParams::RHOVY], cp_SE[CellParams::RHO]);
      Vz0  = divideIfNonZero(cp_SE[CellParams::RHOVZ], cp_SE[CellParams::RHO]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vy0  = divideIfNonZero(cp_SE[CellParams::RHOVY1], cp_SE[CellParams::RHO1]);
      Vz0  = divideIfNonZero(cp_SE[CellParams::RHOVZ1], cp_SE[CellParams::RHO1]);
   }
   
   // 1st order terms:
   Real Ex_SE = By_S*Vz0 - Bz_E*Vy0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_SE += +HALF*((By_S - HALF*dBydz_S)*(+derivs_SE[fs::dVzdy] - derivs_SE[fs::dVzdz]) - dBydz_S*Vz0 + SIXTH*dBydx_S*derivs_SE[fs::dVzdx]);
      Ex_SE += -HALF*((Bz_E + HALF*dBzdy_E)*(+derivs_SE[fs::dVydy] - derivs_SE[fs::dVydz]) + dBzdy_E*Vy0 + SIXTH*dBzdx_E*derivs_SE[fs::dVydx]);
   #endif
   
   const CellID nbrID_SE      = getNeighbourID(mpiGrid, cellID, 2+1, 2-1, 2  );
   #ifndef NDEBUG
      if (nbrID_SE == INVALID_CELLID) {cerr << "ERROR: Could not find SE cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_SE     = mpiGrid[nbrID_SE]->parameters;
   creal* const nbr_derivs_SE = mpiGrid[nbrID_SE]->derivatives;
   c_y = calculateFastMSspeedYZ(cp_SE,derivs_SE,nbr_cp_SE,nbr_derivs_SE,By_S,Bz_E,dBydx_S,dBydz_S,dBzdx_E,dBzdy_E,PLUS,MINUS, RKCase);
   c_z = c_y;
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   
   // Ex and characteristic speeds on k-1 neighbour:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vy0  = divideIfNonZero(cp_NW[CellParams::RHOVY], cp_NW[CellParams::RHO]);
      Vz0  = divideIfNonZero(cp_NW[CellParams::RHOVZ], cp_NW[CellParams::RHO]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vy0  = divideIfNonZero(cp_NW[CellParams::RHOVY1], cp_NW[CellParams::RHO1]);
      Vz0  = divideIfNonZero(cp_NW[CellParams::RHOVZ1], cp_NW[CellParams::RHO1]);
   }
   
   // 1st order terms:
   Real Ex_NW    = By_N*Vz0 - Bz_W*Vy0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_NW += +HALF*((By_N + HALF*dBydz_N)*(-derivs_NW[fs::dVzdy] + derivs_NW[fs::dVzdz]) + dBydz_N*Vz0 + SIXTH*dBydx_N*derivs_NW[fs::dVzdx]);
      Ex_NW += -HALF*((Bz_W - HALF*dBzdy_W)*(-derivs_NW[fs::dVydy] + derivs_NW[fs::dVydz]) - dBzdy_W*Vy0 + SIXTH*dBzdx_W*derivs_NW[fs::dVydx]);
   #endif
   
   const CellID nbrID_NW      = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2-1);
   #ifndef NDEBUG
      if (nbrID_NW == INVALID_CELLID) {cerr << "ERROR: Could not find NW cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_NW     = mpiGrid[nbrID_NW]->parameters;
   creal* const nbr_derivs_NW = mpiGrid[nbrID_NW]->derivatives;
   c_y = calculateFastMSspeedYZ(cp_NW,derivs_NW,nbr_cp_NW,nbr_derivs_NW,By_N,Bz_W,dBydx_N,dBydz_N,dBzdx_W,dBzdy_W,MINUS,PLUS, RKCase);
   c_z = c_y;
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   
   // Ex and characteristic speeds on j-1,k-1 neighbour:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vy0 = divideIfNonZero(cp_NE[CellParams::RHOVY], cp_NE[CellParams::RHO]);
      Vz0 = divideIfNonZero(cp_NE[CellParams::RHOVZ], cp_NE[CellParams::RHO]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vy0 = divideIfNonZero(cp_NE[CellParams::RHOVY1], cp_NE[CellParams::RHO1]);
      Vz0 = divideIfNonZero(cp_NE[CellParams::RHOVZ1], cp_NE[CellParams::RHO1]);
   }
   
   // 1st order terms:
   Real Ex_NE    = By_N*Vz0 - Bz_E*Vy0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ex_NE += +HALF*((By_N + HALF*dBydz_N)*(+derivs_NE[fs::dVzdy] + derivs_NE[fs::dVzdz]) + dBydz_N*Vz0 + SIXTH*dBydx_N*derivs_NE[fs::dVzdx]);
      Ex_NE += -HALF*((Bz_E + HALF*dBzdy_E)*(+derivs_NE[fs::dVydy] + derivs_NE[fs::dVydz]) + dBzdy_E*Vy0 + SIXTH*dBzdx_E*derivs_NE[fs::dVydx]);
   #endif
   
   const CellID nbrID_NE      = getNeighbourID(mpiGrid, cellID, 2+1, 2-1, 2-1);
   #ifndef NDEBUG
      if (nbrID_NE == INVALID_CELLID) {cerr << "ERROR: Could not find NE cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_NE     = mpiGrid[nbrID_NE]->parameters;
   creal* const nbr_derivs_NE = mpiGrid[nbrID_NE]->derivatives;
   c_y = calculateFastMSspeedYZ(cp_NE,derivs_NE,nbr_cp_NE,nbr_derivs_NE,By_N,Bz_E,dBydx_N,dBydz_N,dBzdx_E,dBzdy_E,PLUS,PLUS, RKCase);
   c_z = c_y;
   ay_neg   = max(ay_neg,-Vy0 + c_y);
   ay_pos   = max(ay_pos,+Vy0 + c_y);
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   
   // Calculate properly upwinded edge-averaged Ex:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      cp_SW[CellParams::EX]  = ay_pos*az_pos*Ex_NE + ay_pos*az_neg*Ex_SE + ay_neg*az_pos*Ex_NW + ay_neg*az_neg*Ex_SW;
      cp_SW[CellParams::EX] /= ((ay_pos+ay_neg)*(az_pos+az_neg)+EPS);
      #ifdef FS_1ST_ORDER_SPACE
      // 1st order diffusive terms:
      cp_SW[CellParams::EX] -= az_pos*az_neg/(az_pos+az_neg+EPS)*(By_S-By_N);
      cp_SW[CellParams::EX] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(Bz_W-Bz_E);
      #else
      // 2nd order diffusive terms
      cp_SW[CellParams::EX] -= az_pos*az_neg/(az_pos+az_neg+EPS)*((By_S-HALF*dBydz_S) - (By_N+HALF*dBydz_N));
      cp_SW[CellParams::EX] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((Bz_W-HALF*dBzdy_W) - (Bz_E+HALF*dBzdy_E));
      #endif
   } else { // RKCase == RK_ORDER2_STEP1
      cp_SW[CellParams::EX1]  = ay_pos*az_pos*Ex_NE + ay_pos*az_neg*Ex_SE + ay_neg*az_pos*Ex_NW + ay_neg*az_neg*Ex_SW;
      cp_SW[CellParams::EX1] /= ((ay_pos+ay_neg)*(az_pos+az_neg)+EPS);
      #ifdef FS_1ST_ORDER_SPACE
      // 1st order diffusive terms:
      cp_SW[CellParams::EX1] -= az_pos*az_neg/(az_pos+az_neg+EPS)*(By_S-By_N);
      cp_SW[CellParams::EX1] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(Bz_W-Bz_E);
      #else
      // 2nd order diffusive terms
      cp_SW[CellParams::EX1] -= az_pos*az_neg/(az_pos+az_neg+EPS)*((By_S-HALF*dBydz_S) - (By_N+HALF*dBydz_N));
      cp_SW[CellParams::EX1] += ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((Bz_W-HALF*dBzdy_W) - (Bz_E+HALF*dBzdy_E));
      #endif
   }
   
   //compute maximum timestep for fieldsolver in this cell (CFL=1)      
   Real max_a=ZERO;
   max_a=max(fabs(az_neg),max_a); 
   max_a=max(fabs(az_pos),max_a);
   max_a=max(fabs(ay_neg),max_a);
   max_a=max(fabs(ay_pos),max_a);
   Real min_dx=std::numeric_limits<Real>::max();
   min_dx=min(min_dx,cp_SW[CellParams::DY]);
   min_dx=min(min_dx,cp_SW[CellParams::DZ]);
   //update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
   if(max_a!=ZERO) cp_SW[CellParams::MAXFDT]=min(cp_SW[CellParams::MAXFDT],min_dx/max_a);
}

/*! \brief Low-level electric field propagation function.
 * 
 * Computes the upwinded electric field Y component along the cell's corresponding edge as the cross product of B and V in the XZ plane. Also includes the calculation of the maximally allowed time step.
 * Selects the RHO/RHO1 RHOV[XYZ]/RHOV[XYZ]1 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * \param cellID Index of the cell to process
 * \param mpiGrid Grid
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
static void calculateEdgeElectricFieldY(
   const CellID& cellID,
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   cint& RKCase
) {
   // An edge has four neighbouring spatial cells. Calculate
   // electric field in each of the four cells per edge. 
   namespace fs = fieldsolver;
   
   Real ax_pos,ax_neg;              // Max. characteristic velocities to x-direction
   Real az_pos,az_neg;              // Max. characteristic velocities to z-direction
   Real Vx0,Vz0;                    // Reconstructed V
   Real c_x,c_z;                    // Wave speeds to xz-directions
   
   // Get read-only pointers to NE,NW,SE,SW states (SW is rw, result is written there):
   const CellID nbr_SE = getNeighbourID(mpiGrid, cellID, 2  , 2  , 2-1);
   const CellID nbr_NW = getNeighbourID(mpiGrid, cellID, 2-1, 2  , 2  );
   const CellID nbr_NE = getNeighbourID(mpiGrid, cellID, 2-1, 2  , 2-1);
   #ifndef NDEBUG
      if (nbr_SE == INVALID_CELLID) {cerr << "ERROR: Could not find SE neighbour!" << endl; exit(1);}
      if (nbr_NE == INVALID_CELLID) {cerr << "ERROR: Could not find NE neighbour!" << endl; exit(1);}
      if (nbr_NW == INVALID_CELLID) {cerr << "ERROR: Could not find NW neighbour!" << endl; exit(1);}
   #endif
   
   Real* const  cp_SW = mpiGrid[cellID]->parameters;
   creal* const cp_SE = mpiGrid[nbr_SE]->parameters;
   creal* const cp_NE = mpiGrid[nbr_NE]->parameters;
   creal* const cp_NW = mpiGrid[nbr_NW]->parameters;
   
   creal* const derivs_SW = mpiGrid[cellID]->derivatives;
   creal* const derivs_SE = mpiGrid[nbr_SE]->derivatives;
   creal* const derivs_NE = mpiGrid[nbr_NE]->derivatives;
   creal* const derivs_NW = mpiGrid[nbr_NW]->derivatives;
   
   // Fetch required plasma parameters:
   Real Bz_S, Bx_W, Bx_E, Bz_N;
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Bz_S = cp_SW[CellParams::PERBZ]+cp_SW[CellParams::BGBZ];
      Bx_W = cp_SW[CellParams::PERBX]+cp_SW[CellParams::BGBX];
      Bx_E = cp_SE[CellParams::PERBX]+cp_SE[CellParams::BGBX];
      Bz_N = cp_NW[CellParams::PERBZ]+cp_NW[CellParams::BGBZ];
      Vx0  = divideIfNonZero(cp_SW[CellParams::RHOVX], cp_SW[CellParams::RHO]);
      Vz0  = divideIfNonZero(cp_SW[CellParams::RHOVZ], cp_SW[CellParams::RHO]);
   } else { // RKCase == RK_ORDER2_STEP1
      Bz_S = cp_SW[CellParams::PERBZ1]+cp_SW[CellParams::BGBZ];
      Bx_W = cp_SW[CellParams::PERBX1]+cp_SW[CellParams::BGBX];
      Bx_E = cp_SE[CellParams::PERBX1]+cp_SE[CellParams::BGBX];
      Bz_N = cp_NW[CellParams::PERBZ1]+cp_NW[CellParams::BGBZ];
      Vx0  = divideIfNonZero(cp_SW[CellParams::RHOVX1], cp_SW[CellParams::RHO1]);
      Vz0  = divideIfNonZero(cp_SW[CellParams::RHOVZ1], cp_SW[CellParams::RHO1]);
   }
   
   creal dBxdy_W = derivs_SW[fs::dBxdy];
   creal dBxdz_W = derivs_SW[fs::dBxdz];
   creal dBzdx_S = derivs_SW[fs::dBzdx];
   creal dBzdy_S = derivs_SW[fs::dBzdy];
   creal dBxdy_E = derivs_SE[fs::dBxdy];
   creal dBxdz_E = derivs_SE[fs::dBxdz];
   creal dBzdx_N = derivs_NW[fs::dBzdx];
   creal dBzdy_N = derivs_NW[fs::dBzdy];
   
   // Ey and characteristic speeds on this cell:
   // 1st order terms:
   Real Ey_SW  = Bz_S*Vx0 - Bx_W*Vz0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms
      Ey_SW += +HALF*((Bz_S - HALF*dBzdx_S)*(-derivs_SW[fs::dVxdx] - derivs_SW[fs::dVxdz]) - dBzdx_S*Vx0 + SIXTH*dBzdy_S*derivs_SW[fs::dVxdy]);
      Ey_SW += -HALF*((Bx_W - HALF*dBxdz_W)*(-derivs_SW[fs::dVzdx] - derivs_SW[fs::dVzdz]) - dBxdz_W*Vz0 + SIXTH*dBxdy_W*derivs_SW[fs::dVzdy]);
   #endif
   
   const CellID nbrID_SW      = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2  );
   #ifndef NDEBUG
      if (nbrID_SW == INVALID_CELLID) {cerr << "ERROR: Could not find SW cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_SW     = mpiGrid[nbrID_SW]->parameters;
   creal* const nbr_derivs_SW = mpiGrid[nbrID_SW]->derivatives;
   c_z = calculateFastMSspeedXZ(cp_SW,derivs_SW,nbr_cp_SW,nbr_derivs_SW,Bx_W,Bz_S,dBxdy_W,dBxdz_W,dBzdx_S,dBzdy_S,MINUS,MINUS, RKCase);
   c_x = c_z;
   az_neg   = max(ZERO,-Vz0 + c_z);
   az_pos   = max(ZERO,+Vz0 + c_z);
   ax_neg   = max(ZERO,-Vx0 + c_x);
   ax_pos   = max(ZERO,+Vx0 + c_x);
   
   // Ey and characteristic speeds on k-1 neighbour:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vx0  = divideIfNonZero(cp_SE[CellParams::RHOVX], cp_SE[CellParams::RHO]);
      Vz0  = divideIfNonZero(cp_SE[CellParams::RHOVZ], cp_SE[CellParams::RHO]);
   } else { //RKCase == RK_ORDER2_STEP1
      Vx0  = divideIfNonZero(cp_SE[CellParams::RHOVX1], cp_SE[CellParams::RHO1]);
      Vz0  = divideIfNonZero(cp_SE[CellParams::RHOVZ1], cp_SE[CellParams::RHO1]);
   }
   
   // 1st order terms:
   Real Ey_SE    = Bz_S*Vx0 - Bx_E*Vz0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ey_SE += +HALF*((Bz_S - HALF*dBzdx_S)*(-derivs_SE[fs::dVxdx] + derivs_SE[fs::dVxdz]) - dBzdx_S*Vx0 + SIXTH*dBzdy_S*derivs_SE[fs::dVxdy]);
      Ey_SE += -HALF*((Bx_E + HALF*dBxdz_E)*(-derivs_SE[fs::dVzdx] + derivs_SE[fs::dVzdz]) + dBxdz_E*Vz0 + SIXTH*dBxdy_E*derivs_SE[fs::dVzdy]);
   #endif
   
   const CellID nbrID_SE      = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2-1);
   #ifndef NDEBUG
      if (nbrID_SE == INVALID_CELLID) {cerr << "ERROR: Could not find SE cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_SE     = mpiGrid[nbrID_SE]->parameters;
   creal* const nbr_derivs_SE = mpiGrid[nbrID_SE]->derivatives;
   c_z = calculateFastMSspeedXZ(cp_SE,derivs_SE,nbr_cp_SE,nbr_derivs_SE,Bx_E,Bz_S,dBxdy_E,dBxdz_E,dBzdx_S,dBzdy_S,MINUS,PLUS, RKCase);
   c_x = c_z;
   az_neg   = max(az_neg,-Vz0 - c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 - c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   
   // Ey and characteristic speeds on i-1 neighbour:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vz0  = divideIfNonZero(cp_NW[CellParams::RHOVZ], cp_NW[CellParams::RHO]);
      Vx0  = divideIfNonZero(cp_NW[CellParams::RHOVX], cp_NW[CellParams::RHO]);
   } else { //RKCase == RK_ORDER2_STEP1
      Vz0  = divideIfNonZero(cp_NW[CellParams::RHOVZ1], cp_NW[CellParams::RHO1]);
      Vx0  = divideIfNonZero(cp_NW[CellParams::RHOVX1], cp_NW[CellParams::RHO1]);
   }
   
   // 1st order terms:
   Real Ey_NW    = Bz_N*Vx0 - Bx_W*Vz0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ey_NW += +HALF*((Bz_N + HALF*dBzdx_N)*(+derivs_NW[fs::dVxdx] - derivs_NW[fs::dVxdz]) + dBzdx_N*Vx0 + SIXTH*dBzdy_N*derivs_NW[fs::dVxdy]);
      Ey_NW += -HALF*((Bx_W - HALF*dBxdz_W)*(+derivs_NW[fs::dVzdx] - derivs_NW[fs::dVzdz]) - dBxdz_W*Vz0 + SIXTH*dBxdy_W*derivs_NW[fs::dVzdy]);
   #endif
   
   const CellID nbrID_NW      = getNeighbourID(mpiGrid, cellID, 2-1, 2+1, 2  );
   #ifndef NDEBUG
      if (nbrID_NW == INVALID_CELLID) {cerr << "ERROR: Could not find NW cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_NW     = mpiGrid[nbrID_NW]->parameters;
   creal* const nbr_derivs_NW = mpiGrid[nbrID_NW]->derivatives;
   c_z = calculateFastMSspeedXZ(cp_NW,derivs_NW,nbr_cp_NW,nbr_derivs_NW,Bx_W,Bz_N,dBxdy_W,dBxdz_W,dBzdx_N,dBzdy_N,PLUS,MINUS, RKCase);
   c_x = c_z;
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   
   // Ey and characteristic speeds on i-1,k-1 neighbour:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vz0 = divideIfNonZero(cp_NE[CellParams::RHOVZ], cp_NE[CellParams::RHO]);
      Vx0 = divideIfNonZero(cp_NE[CellParams::RHOVX], cp_NE[CellParams::RHO]);
   } else { //RKCase == RK_ORDER2_STEP1
      Vz0 = divideIfNonZero(cp_NE[CellParams::RHOVZ1], cp_NE[CellParams::RHO1]);
      Vx0 = divideIfNonZero(cp_NE[CellParams::RHOVX1], cp_NE[CellParams::RHO1]);
   }
   
   // 1st order terms:
   Real Ey_NE    = Bz_N*Vx0 - Bx_E*Vz0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ey_NE += +HALF*((Bz_N + HALF*dBzdx_N)*(+derivs_NE[fs::dVxdx] + derivs_NE[fs::dVxdz]) + dBzdx_N*Vx0 + SIXTH*dBzdy_N*derivs_NE[fs::dVxdy]);
      Ey_NE += -HALF*((Bx_E + HALF*dBxdz_E)*(+derivs_NE[fs::dVzdx] + derivs_NE[fs::dVzdz]) + dBxdz_E*Vz0 + SIXTH*dBxdy_E*derivs_NE[fs::dVzdy]);
   #endif
   
   const CellID nbrID_NE      = getNeighbourID(mpiGrid, cellID, 2-1, 2+1, 2-1);
   #ifndef NDEBUG
      if (nbrID_NE == INVALID_CELLID) {cerr << "ERROR: Could not find NE cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_NE     = mpiGrid[nbrID_NE]->parameters;
   creal* const nbr_derivs_NE = mpiGrid[nbrID_NE]->derivatives;
   c_z = calculateFastMSspeedXZ(cp_NE,derivs_NE,nbr_cp_NE,nbr_derivs_NE,Bx_E,Bz_N,dBxdy_E,dBxdz_E,dBzdx_N,dBzdy_N,PLUS,PLUS, RKCase);
   c_x = c_z;
   az_neg   = max(az_neg,-Vz0 + c_z);
   az_pos   = max(az_pos,+Vz0 + c_z);
   ax_neg   = max(ax_neg,-Vx0 + c_x);
   ax_pos   = max(ax_pos,+Vx0 + c_x);
   
   // Calculate properly upwinded edge-averaged Ey:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      cp_SW[CellParams::EY]  = az_pos*ax_pos*Ey_NE + az_pos*ax_neg*Ey_SE + az_neg*ax_pos*Ey_NW + az_neg*ax_neg*Ey_SW;
      cp_SW[CellParams::EY] /= ((az_pos+az_neg)*(ax_pos+ax_neg)+EPS);
      #ifdef FS_1ST_ORDER_SPACE
      cp_SW[CellParams::EY] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(Bz_S-Bz_N);
      cp_SW[CellParams::EY] += az_pos*az_neg/(az_pos+az_neg+EPS)*(Bx_W-Bx_E);
      #else
      cp_SW[CellParams::EY] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((Bz_S-HALF*dBzdx_S) - (Bz_N+HALF*dBzdx_N));
      cp_SW[CellParams::EY] += az_pos*az_neg/(az_pos+az_neg+EPS)*((Bx_W-HALF*dBxdz_W) - (Bx_E+HALF*dBxdz_E));
      #endif
   } else { // RKCase == RK_ORDER2_STEP1
      cp_SW[CellParams::EY1]  = az_pos*ax_pos*Ey_NE + az_pos*ax_neg*Ey_SE + az_neg*ax_pos*Ey_NW + az_neg*ax_neg*Ey_SW;
      cp_SW[CellParams::EY1] /= ((az_pos+az_neg)*(ax_pos+ax_neg)+EPS);
      #ifdef FS_1ST_ORDER_SPACE
      cp_SW[CellParams::EY1] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(Bz_S-Bz_N);
      cp_SW[CellParams::EY1] += az_pos*az_neg/(az_pos+az_neg+EPS)*(Bx_W-Bx_E);
      #else
      cp_SW[CellParams::EY1] -= ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((Bz_S-HALF*dBzdx_S) - (Bz_N+HALF*dBzdx_N));
      cp_SW[CellParams::EY1] += az_pos*az_neg/(az_pos+az_neg+EPS)*((Bx_W-HALF*dBxdz_W) - (Bx_E+HALF*dBxdz_E));
      #endif
   }
   
   //compute maximum timestep for fieldsolver in this cell (CFL=1)      
   Real max_a=ZERO;
   max_a=max(fabs(az_neg),max_a);
   max_a=max(fabs(az_pos),max_a);
   max_a=max(fabs(ax_neg),max_a);
   max_a=max(fabs(ax_pos),max_a);
   Real min_dx=std::numeric_limits<Real>::max();;
   min_dx=min(min_dx,cp_SW[CellParams::DX]);
   min_dx=min(min_dx,cp_SW[CellParams::DZ]);
   //update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps 
   if(max_a!=ZERO) cp_SW[CellParams::MAXFDT]=min(cp_SW[CellParams::MAXFDT],min_dx/max_a);
}

/*! \brief Low-level electric field propagation function.
 *
 * Computes the upwinded electric field Z component along the cell's corresponding edge as the cross product of B and V in the XY plane. Also includes the calculation of the maximally allowed time step.
 * Selects the RHO/RHO1 RHOV[XYZ]/RHOV[XYZ]1 and B[XYZ]/B[XYZ]1 values depending on the stage of the Runge-Kutta time stepping method.
 * \param cellID Index of the cell to process
 * \param mpiGrid Grid
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
static void calculateEdgeElectricFieldZ(
   const CellID& cellID,
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   cint& RKCase
) {
   namespace fs = fieldsolver;
   namespace pc = physicalconstants;
   
   // An edge has four neighbouring spatial cells. Calculate 
   // electric field in each of the four cells per edge.
   Real ax_pos,ax_neg;              // Max. characteristic velocities to x-direction
   Real ay_pos,ay_neg;              // Max. characteristic velocities to y-direction
   Real Vx0,Vy0;                    // Reconstructed V
   Real c_x,c_y;                    // Characteristic speeds to xy-directions
   
   // Get read-only pointers to NE,NW,SE,SW states (SW is rw, result is written there):
   const CellID nbr_SE = getNeighbourID(mpiGrid, cellID, 2-1, 2  , 2  );
   const CellID nbr_NE = getNeighbourID(mpiGrid, cellID, 2-1, 2-1, 2  );
   const CellID nbr_NW = getNeighbourID(mpiGrid, cellID, 2  , 2-1, 2  );
   #ifndef NDEBUG
      if (nbr_SE == INVALID_CELLID) {cerr << "ERROR: Could not find SE neighbour!" << endl; exit(1);}
      if (nbr_NE == INVALID_CELLID) {cerr << "ERROR: Could not find NE neighbour!" << endl; exit(1);}
      if (nbr_NW == INVALID_CELLID) {cerr << "ERROR: Could not find NW neighbour!" << endl; exit(1);}
   #endif
   
   Real* const cp_SW  = mpiGrid[cellID]->parameters;
   creal* const cp_SE = mpiGrid[nbr_SE]->parameters;
   creal* const cp_NE = mpiGrid[nbr_NE]->parameters;
   creal* const cp_NW = mpiGrid[nbr_NW]->parameters;
   
   creal* const derivs_SW = mpiGrid[cellID]->derivatives;
   creal* const derivs_SE = mpiGrid[nbr_SE]->derivatives;
   creal* const derivs_NE = mpiGrid[nbr_NE]->derivatives;
   creal* const derivs_NW = mpiGrid[nbr_NW]->derivatives;
   
   // Fetch needed plasma parameters/derivatives from the four cells:
   Real Bx_S, By_W, By_E, Bx_N;
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Bx_S    = cp_SW[CellParams::PERBX] + cp_SW[CellParams::BGBX];
      By_W    = cp_SW[CellParams::PERBY] + cp_SW[CellParams::BGBY];
      By_E    = cp_SE[CellParams::PERBY] + cp_SE[CellParams::BGBY];
      Bx_N    = cp_NW[CellParams::PERBX] + cp_NW[CellParams::BGBX];
      Vx0  = divideIfNonZero(cp_SW[CellParams::RHOVX], cp_SW[CellParams::RHO]);
      Vy0  = divideIfNonZero(cp_SW[CellParams::RHOVY], cp_SW[CellParams::RHO]);
   } else { // RKCase == RK_ORDER2_STEP1
      Bx_S    = cp_SW[CellParams::PERBX1] + cp_SW[CellParams::BGBX];
      By_W    = cp_SW[CellParams::PERBY1] + cp_SW[CellParams::BGBY];
      By_E    = cp_SE[CellParams::PERBY1] + cp_SE[CellParams::BGBY];
      Bx_N    = cp_NW[CellParams::PERBX1] + cp_NW[CellParams::BGBX];
      Vx0  = divideIfNonZero(cp_SW[CellParams::RHOVX1], cp_SW[CellParams::RHO1]);
      Vy0  = divideIfNonZero(cp_SW[CellParams::RHOVY1], cp_SW[CellParams::RHO1]);
   }
   
   creal dBxdy_S = derivs_SW[fs::dBxdy];
   creal dBxdz_S = derivs_SW[fs::dBxdz];
   creal dBydx_W = derivs_SW[fs::dBydx];
   creal dBydz_W = derivs_SW[fs::dBydz];
   creal dBydx_E = derivs_SE[fs::dBydx];
   creal dBydz_E = derivs_SE[fs::dBydz];
   creal dBxdy_N = derivs_NW[fs::dBxdy];
   creal dBxdz_N = derivs_NW[fs::dBxdz];
   
   // Ez and characteristic speeds on SW cell:
   // 1st order terms:
   Real Ez_SW = Bx_S*Vy0 - By_W*Vx0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_SW  += +HALF*((Bx_S - HALF*dBxdy_S)*(-derivs_SW[fs::dVydx] - derivs_SW[fs::dVydy]) - dBxdy_S*Vy0 + SIXTH*dBxdz_S*derivs_SW[fs::dVydz]);
      Ez_SW  += -HALF*((By_W - HALF*dBydx_W)*(-derivs_SW[fs::dVxdx] - derivs_SW[fs::dVxdy]) - dBydx_W*Vx0 + SIXTH*dBydz_W*derivs_SW[fs::dVxdz]);
   #endif
   
   // Calculate maximum wave speed (fast magnetosonic speed) on SW cell. In order 
   // to get Alfven speed we need to calculate some reconstruction coeff. for Bz:
   const CellID nbrID_SW      = getNeighbourID(mpiGrid, cellID, 2  , 2  , 2+1);
   #ifndef NDEBUG
      if (nbrID_SW == INVALID_CELLID) {cerr << "ERROR: Could not find SW cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_SW     = mpiGrid[nbrID_SW]->parameters;
   creal* const nbr_derivs_SW = mpiGrid[nbrID_SW]->derivatives;
   c_x = calculateFastMSspeedXY(cp_SW,derivs_SW,nbr_cp_SW,nbr_derivs_SW,Bx_S,By_W,dBxdy_S,dBxdz_S,dBydx_W,dBydz_W,MINUS,MINUS, RKCase);
   c_y = c_x;
   ax_neg   = max(ZERO,-Vx0 + c_x);
   ax_pos   = max(ZERO,+Vx0 + c_x);
   ay_neg   = max(ZERO,-Vy0 + c_y);
   ay_pos   = max(ZERO,+Vy0 + c_y);
   
   // Ez and characteristic speeds on SE (i-1) cell:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vx0  = divideIfNonZero(cp_SE[CellParams::RHOVX], cp_SE[CellParams::RHO]);
      Vy0  = divideIfNonZero(cp_SE[CellParams::RHOVY], cp_SE[CellParams::RHO]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vx0  = divideIfNonZero(cp_SE[CellParams::RHOVX1], cp_SE[CellParams::RHO1]);
      Vy0  = divideIfNonZero(cp_SE[CellParams::RHOVY1], cp_SE[CellParams::RHO1]);
   }
   
   // 1st order terms:
   Real Ez_SE = Bx_S*Vy0 - By_E*Vx0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_SE  += +HALF*((Bx_S - HALF*dBxdy_S)*(+derivs_SE[fs::dVydx] - derivs_SE[fs::dVydy]) - dBxdy_S*Vy0 + SIXTH*dBxdz_S*derivs_SE[fs::dVydz]);
      Ez_SE  += -HALF*((By_E + HALF*dBydx_E)*(+derivs_SE[fs::dVxdx] - derivs_SE[fs::dVxdy]) + dBydx_E*Vx0 + SIXTH*dBydz_E*derivs_SE[fs::dVxdz]);
   #endif
   
   const CellID nbrID_SE      = getNeighbourID(mpiGrid, cellID, 2-1, 2  , 2+1);
   #ifndef NDEBUG
      if (nbrID_SE == INVALID_CELLID) {cerr << "ERROR: Could not find SE cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_SE     = mpiGrid[nbrID_SE]->parameters;
   creal* const nbr_derivs_SE = mpiGrid[nbrID_SE]->derivatives;
   c_x = calculateFastMSspeedXY(cp_SE,derivs_SE,nbr_cp_SE,nbr_derivs_SE,Bx_S,By_E,dBxdy_S,dBxdz_S,dBydx_E,dBydz_E,PLUS,MINUS, RKCase);
   c_y = c_x;   
   ax_neg = max(ax_neg,-Vx0 + c_x);
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);
   
   // Ez and characteristic speeds on NW (j-1) cell:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vx0  = divideIfNonZero(cp_NW[CellParams::RHOVX], cp_NW[CellParams::RHO]);
      Vy0  = divideIfNonZero(cp_NW[CellParams::RHOVY], cp_NW[CellParams::RHO]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vx0  = divideIfNonZero(cp_NW[CellParams::RHOVX1], cp_NW[CellParams::RHO1]);
      Vy0  = divideIfNonZero(cp_NW[CellParams::RHOVY1], cp_NW[CellParams::RHO1]);
   }
   
   // 1st order terms:
   Real Ez_NW = Bx_N*Vy0 - By_W*Vx0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_NW  += +HALF*((Bx_N + HALF*dBxdy_N)*(-derivs_NW[fs::dVydx] + derivs_NW[fs::dVydy]) + dBxdy_N*Vy0 + SIXTH*dBxdz_N*derivs_NW[fs::dVydz]);
      Ez_NW  += -HALF*((By_W - HALF*dBydx_W)*(-derivs_NW[fs::dVxdx] + derivs_NW[fs::dVxdy]) - dBydx_W*Vx0 + SIXTH*dBydz_W*derivs_NW[fs::dVxdz]);
   #endif
   
   const CellID nbrID_NW      = getNeighbourID(mpiGrid, cellID, 2  , 2-1, 2+1);
   #ifndef NDEBUG
      if (nbrID_NW == INVALID_CELLID) {cerr << "ERROR: Could not find NW cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_NW     = mpiGrid[nbrID_NW]->parameters;
   creal* const nbr_derivs_NW = mpiGrid[nbrID_NW]->derivatives;
   c_x = calculateFastMSspeedXY(cp_NW,derivs_NW,nbr_cp_NW,nbr_derivs_NW,Bx_N,By_W,dBxdy_N,dBxdz_N,dBydx_W,dBydz_W,MINUS,PLUS, RKCase);
   c_y = c_x;
   ax_neg = max(ax_neg,-Vx0 + c_x); 
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);
   
   // Ez and characteristic speeds on NE (i-1,j-1) cell:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      Vx0  = divideIfNonZero(cp_NE[CellParams::RHOVX], cp_NE[CellParams::RHO]);
      Vy0  = divideIfNonZero(cp_NE[CellParams::RHOVY], cp_NE[CellParams::RHO]);
   } else { // RKCase == RK_ORDER2_STEP1
      Vx0  = divideIfNonZero(cp_NE[CellParams::RHOVX1], cp_NE[CellParams::RHO1]);
      Vy0  = divideIfNonZero(cp_NE[CellParams::RHOVY1], cp_NE[CellParams::RHO1]);
   }
   
   // 1st order terms:
   Real Ez_NE = Bx_N*Vy0 - By_E*Vx0;
   #ifndef FS_1ST_ORDER_SPACE
      // 2nd order terms:
      Ez_NE  += +HALF*((Bx_N + HALF*dBxdy_N)*(+derivs_NE[fs::dVydx] + derivs_NE[fs::dVydy]) + dBxdy_N*Vy0 + SIXTH*dBxdz_N*derivs_NE[fs::dVydz]);
      Ez_NE  += -HALF*((By_E + HALF*dBydx_E)*(+derivs_NE[fs::dVxdx] + derivs_NE[fs::dVxdy]) + dBydx_E*Vx0 + SIXTH*dBydz_E*derivs_NE[fs::dVxdz]);
   #endif
   
   const CellID nbrID_NE      = getNeighbourID(mpiGrid, cellID, 2-1, 2-1, 2+1);
   #ifndef NDEBUG
      if (nbrID_NE == INVALID_CELLID) {cerr << "ERROR: Could not find NE cell!" << endl; exit(1);}
   #endif
   creal* const nbr_cp_NE     = mpiGrid[nbrID_NE]->parameters;
   creal* const nbr_derivs_NE = mpiGrid[nbrID_NE]->derivatives;
   c_x = calculateFastMSspeedXY(cp_NE,derivs_NE,nbr_cp_NE,nbr_derivs_NE,Bx_N,By_E,dBxdy_N,dBxdz_N,dBydx_E,dBydz_E,PLUS,PLUS, RKCase);
   c_y = c_x;
   ax_neg = max(ax_neg,-Vx0 + c_x);
   ax_pos = max(ax_pos,+Vx0 + c_x);
   ay_neg = max(ay_neg,-Vy0 + c_y);
   ay_pos = max(ay_pos,+Vy0 + c_y);
   
   // Calculate properly upwinded edge-averaged Ez:
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      cp_SW[CellParams::EZ] = ax_pos*ay_pos*Ez_NE + ax_pos*ay_neg*Ez_SE + ax_neg*ay_pos*Ez_NW + ax_neg*ay_neg*Ez_SW;
      CHECK_FLOAT(cp_SW[CellParams::EZ])
      cp_SW[CellParams::EZ] /= ((ax_pos+ax_neg)*(ay_pos+ay_neg)+EPS);
      CHECK_FLOAT(cp_SW[CellParams::EZ])
      #ifdef FS_1ST_ORDER_SPACE
      cp_SW[CellParams::EZ] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(Bx_S-Bx_N);
      CHECK_FLOAT(cp_SW[CellParams::EZ])
      cp_SW[CellParams::EZ] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(By_W-By_E);
      CHECK_FLOAT(cp_SW[CellParams::EZ])
      #else
      cp_SW[CellParams::EZ] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((Bx_S-HALF*dBxdy_S) - (Bx_N+HALF*dBxdy_N));
      CHECK_FLOAT(cp_SW[CellParams::EZ])
      cp_SW[CellParams::EZ] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((By_W-HALF*dBydx_W) - (By_E+HALF*dBydx_E));
      CHECK_FLOAT(cp_SW[CellParams::EZ])
      #endif
   } else { // RKCase == RK_ORDER2_STEP1
      cp_SW[CellParams::EZ1] = ax_pos*ay_pos*Ez_NE + ax_pos*ay_neg*Ez_SE + ax_neg*ay_pos*Ez_NW + ax_neg*ay_neg*Ez_SW;
      cp_SW[CellParams::EZ1] /= ((ax_pos+ax_neg)*(ay_pos+ay_neg)+EPS);
      #ifdef FS_1ST_ORDER_SPACE
      cp_SW[CellParams::EZ1] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*(Bx_S-Bx_N);
      cp_SW[CellParams::EZ1] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*(By_W-By_E);
      #else
      cp_SW[CellParams::EZ1] -= ay_pos*ay_neg/(ay_pos+ay_neg+EPS)*((Bx_S-HALF*dBxdy_S) - (Bx_N+HALF*dBxdy_N));
      cp_SW[CellParams::EZ1] += ax_pos*ax_neg/(ax_pos+ax_neg+EPS)*((By_W-HALF*dBydx_W) - (By_E+HALF*dBydx_E));
      #endif
   }
   
   
   //compute maximum timestep for fieldsolver in this cell (CFL=1)      
   Real max_a=ZERO;
   max_a=max(fabs(ay_neg),max_a);
   max_a=max(fabs(ay_pos),max_a);
   max_a=max(fabs(ax_neg),max_a);
   max_a=max(fabs(ax_pos),max_a);
   Real min_dx=std::numeric_limits<Real>::max();;
   min_dx=min(min_dx,cp_SW[CellParams::DX]);
   min_dx=min(min_dx,cp_SW[CellParams::DY]);
   //update max allowed timestep for field propagation in this cell, which is the minimum of CFL=1 timesteps
   if(max_a!=ZERO) cp_SW[CellParams::MAXFDT]=min(cp_SW[CellParams::MAXFDT],min_dx/max_a);
}

/*! \brief Low-level magnetic field propagation function.
 * 
 * Propagates the cell's face-averaged magnetic field components by
 * using Faraday's law on the face edges. Depending on the time order
 * of accuracy it is done in one stage or in two stages using the
 * intermediate E1 components for the first stage of the second-order
 * Runge-Kutta method and E for the other cases.
 * \param cellID Index of the cell to process
 * \param mpiGrid Grid
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
static void propagateMagneticField(
   const CellID& cellID,
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   creal& dt,
   cint& RKCase
) {
   CellID nbrID;
   Real* const cp0 = mpiGrid[cellID]->parameters;
   creal* cp1;
   creal* cp2;
   creal dx = cp0[CellParams::DX];
   CHECK_FLOAT(dx)
   creal dy = cp0[CellParams::DY];
   CHECK_FLOAT(dy)
   creal dz = cp0[CellParams::DZ];
   CHECK_FLOAT(dz)
   
   #ifndef NDEBUG
      map<CellID,uint>::const_iterator it = sysBoundaryFlags.find(cellID);
      if (it == sysBoundaryFlags.end()) {cerr << "Could not find boundary flags for cell #" << cellID << endl; exit(1);}
      cuint sysBoundaryFlag = it->second;
   #else
      cuint sysBoundaryFlag = sysBoundaryFlags[cellID];
   #endif
   
   // Propagate face-averaged Bx:
   if ((sysBoundaryFlag & PROPAGATE_BX) == PROPAGATE_BX) {
      nbrID = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2  );
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive y direction"
            << std::endl;
         abort();
      }
      cp1 = mpiGrid[nbrID]->parameters;

      nbrID = getNeighbourID(mpiGrid, cellID, 2  , 2  , 2+1);
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive z direction"
            << std::endl;
         abort();
      }
      cp2 = mpiGrid[nbrID]->parameters;

      # ifdef FS_1ST_ORDER_TIME
      cp0[CellParams::PERBX] += dt/dz*(cp2[CellParams::EY] - cp0[CellParams::EY]) + dt/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ]);
      # else
      if(RKCase == RK_ORDER2_STEP1) {
         cp0[CellParams::PERBX1] = cp0[CellParams::PERBX] +
         Parameters::RK_alpha*dt*(1.0/dz*(cp2[CellParams::EY] - cp0[CellParams::EY]) + 1.0/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ]));
      } else {
         cp0[CellParams::PERBX] += dt * ((1.0 - 0.5/Parameters::RK_alpha) * (1.0/dz*(cp2[CellParams::EY] -cp0[CellParams::EY]) +
                                                                            1.0/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ])) +
                                        0.5/Parameters::RK_alpha * (1.0/dz*(cp2[CellParams::EY1] - cp0[CellParams::EY1]) +
                                                                    1.0/dy*(cp0[CellParams::EZ1] - cp1[CellParams::EZ1]))
                                        );
      }
      # endif
   }
   
   // Propagate face-averaged By:
   if ((sysBoundaryFlag & PROPAGATE_BY) == PROPAGATE_BY) {
      nbrID = getNeighbourID(mpiGrid, cellID, 2  , 2  , 2+1);
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive z direction"
            << std::endl;
         abort();
      }
      cp1 = mpiGrid[nbrID]->parameters;

      nbrID = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2  );
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive x direction"
            << std::endl;
         abort();
      }
      cp2 = mpiGrid[nbrID]->parameters;

      # ifdef FS_1ST_ORDER_TIME
      cp0[CellParams::PERBY] += dt/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) + dt/dz*(cp0[CellParams::EX] - cp1[CellParams::EX]);
      # else
      if(RKCase == RK_ORDER2_STEP1) {
         cp0[CellParams::PERBY1] = cp0[CellParams::PERBY] +
         Parameters::RK_alpha*dt*(1.0/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) + 1.0/dz*(cp0[CellParams::EX] - cp1[CellParams::EX]));
      } else {
         cp0[CellParams::PERBY] += dt * ((1.0 - 0.5/Parameters::RK_alpha) * (1.0/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) +
                                                                            1.0/dz*(cp0[CellParams::EX] - cp1[CellParams::EX])) + 
                                        0.5/Parameters::RK_alpha * (1.0/dx*(cp2[CellParams::EZ1] - cp0[CellParams::EZ1]) +
                                                                    1.0/dz*(cp0[CellParams::EX1] - cp1[CellParams::EX1]))
                                        );
      }
      # endif
   }
      
   // Propagate face-averaged Bz:
   if ((sysBoundaryFlag & PROPAGATE_BZ) == PROPAGATE_BZ) {
      nbrID = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2  );
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive x direction"
            << std::endl;
         abort();
      }
      cp1 = mpiGrid[nbrID]->parameters;

      nbrID = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2  );
      #ifndef NDEBUG
         if (nbrID == INVALID_CELLID) {cerr << "Failed to get nbr pointer!" << endl; exit(1);}
      #endif

      if (mpiGrid[nbrID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << nbrID
            << " while solving cell " << cellID
            << " in positive y direction"
            << std::endl;
         abort();
      }
      cp2 = mpiGrid[nbrID]->parameters;
      
      # ifdef FS_1ST_ORDER_TIME
      cp0[CellParams::PERBZ] += dt/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) + dt/dx*(cp0[CellParams::EY] - cp1[CellParams::EY]);
      # else
      if(RKCase == RK_ORDER2_STEP1) {
         cp0[CellParams::PERBZ1] = cp0[CellParams::PERBZ] +
         Parameters::RK_alpha*dt*(1.0/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) + 1.0/dx*(cp0[CellParams::EY] - cp1[CellParams::EY]));
      } else {
         cp0[CellParams::PERBZ] += dt * ((1.0 - 0.5/Parameters::RK_alpha) * (1.0/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) +
                                                                          1.0/dx*(cp0[CellParams::EY] - cp1[CellParams::EY])) +
         0.5/Parameters::RK_alpha * (1.0/dy*(cp2[CellParams::EX1] - cp0[CellParams::EX1]) +
                                     1.0/dx*(cp0[CellParams::EY1] - cp1[CellParams::EY1])));
      }
# endif   
   }
}

bool initializeFieldPropagatorAfterRebalance(
        dccrg::Dccrg<SpatialCell>& mpiGrid
) {
   
   vector<uint64_t> localCells = mpiGrid.get_cells();

   calculateSysBoundaryFlags(mpiGrid,localCells);

   //need E when computing magnetic field later on
   // ASSUME STATIC background field, we do not later on explicitly transfer it!
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_E | Transfer::CELL_BGB);
   int timer=phiprof::initializeTimer("Communicate E and BGB","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.update_remote_neighbor_data();
   phiprof::stop(timer);
   
   return true;
}

/*! Calculates bit masks used in the field solver and computes the initial edge electric fields from the initial magnetic fields. Then computes the initial volume averages.
 */
bool initializeFieldPropagator(
        dccrg::Dccrg<SpatialCell>& mpiGrid
) {
   vector<uint64_t> localCells = mpiGrid.get_cells();
   
   calculateSysBoundaryFlags(mpiGrid,localCells);
   
   // Calculate bit masks used for if-statements by field propagator. 
   // These are used to test whether or not certain combination of 
   // neighbours exists for a cell. These can be replaced by honest 
   // if-statements, but you will just end up needing very many of them 
   // as each bit mask tests the existence of several neighbours at once.
   // Existence of neighbours would also need to be queried from the 
   // parallel grid, i.e. using if-statements is likely to be much 
   // slower.

   // x-derivatives are calculated if x-face neighbours exist:
   CALCULATE_DX = 0;
   CALCULATE_DX = CALCULATE_DX | (1 << calcNbrNumber(0,1,1));
   CALCULATE_DX = CALCULATE_DX | (1 << calcNbrNumber(2,1,1));
   
   // y-derivatives are calculated if y-face neighbours exist:
   CALCULATE_DY = 0;
   CALCULATE_DY = CALCULATE_DY | (1 << calcNbrNumber(1,0,1));
   CALCULATE_DY = CALCULATE_DY | (1 << calcNbrNumber(1,2,1));
   
   // z-derivatives are calculated if z-face neighbours exist:
   CALCULATE_DZ = 0;
   CALCULATE_DZ = CALCULATE_DZ | (1 << calcNbrNumber(1,1,0));
   CALCULATE_DZ = CALCULATE_DZ | (1 << calcNbrNumber(1,1,2));
   
   // Edge Ex is calculated i   f -y,-z,+/-x neighbours exist:
   CALCULATE_EX = 0;
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(1,0,1)); // -y nbr
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(1,1,0)); // -z nbr
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(0,1,1)); // -x nbr
   CALCULATE_EX = CALCULATE_EX | (1 << calcNbrNumber(2,1,1)); // +x nbr
      
   // Edge Ey is calculated if -x,-z,+/-y neighbours exist:
   CALCULATE_EY = 0;
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(0,1,1)); // -x nbr
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(1,1,0)); // -z nbr
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(1,0,1)); // -y nbr
   CALCULATE_EY = CALCULATE_EY | (1 << calcNbrNumber(1,2,1)); // +y nbr
      
   // Edge Ez is calculated      if -x,-y,+/-z neighbours exist:
   CALCULATE_EZ = 0;
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(0,1,1)); // -x nbr
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(1,0,1)); // -y nbr
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(1,1,0)); // -z nbr
   CALCULATE_EZ = CALCULATE_EZ | (1 << calcNbrNumber(1,1,2)); // +z nbr

   // Bx is propagated if -x,+/-y,+/-z neighbours exist:
   PROPAGATE_BX = 0;
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(0,1,1)); // -x nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,0,1)); // -y nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,2,1)); // +y nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,1,0)); // -z nbr
   PROPAGATE_BX = PROPAGATE_BX | (1 << calcNbrNumber(1,1,2)); // +z nbr
   
   // By is propagated if -y,+/-x,+/-z neighbours exist:
   PROPAGATE_BY = 0;
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(1,0,1)); // -y nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(0,1,1)); // -x nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(2,1,1)); // +x nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(1,1,0)); // -z nbr
   PROPAGATE_BY = PROPAGATE_BY | (1 << calcNbrNumber(1,1,2)); // +z nbr
   
   // Bz is propagated if -z,+/-x,+/-y neighbours exist:
   PROPAGATE_BZ = 0;
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(1,1,0)); // -z nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(0,1,1)); // -x nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(2,1,1)); // +x nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(1,0,1)); // -y nbr
   PROPAGATE_BZ = PROPAGATE_BZ | (1 << calcNbrNumber(1,2,1)); // +y nbr
   

   // ASSUME STATIC background field, we do not later on explicitly transfer it
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_BGB);
   int timer=phiprof::initializeTimer("Communicate BGB","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.update_remote_neighbor_data();
   phiprof::stop(timer);
   // Calculate derivatives and upwinded edge-E. Exchange derivatives 
   // and edge-E:s between neighbouring processes and calculate 
   // face-averaged E,B fields.
   calculateDerivativesSimple(mpiGrid,localCells, RK_ORDER1);
   calculateUpwindedElectricFieldSimple(mpiGrid,localCells, RK_ORDER1);
   calculateVolumeAveragedFields(mpiGrid);
   
   return true;
}

bool finalizeFieldPropagator(
   dccrg::Dccrg<SpatialCell>& mpiGrid
) {
   return true;
}

/*! \brief High-level derivative calculation wrapper function.
 * 
 * In the first stage of the second-order Runge-Kutta time stepping scheme, this linearly interpolates the velocity moments and puts them into the RHO1, RHOV[XYZ]1 variables.
 * 
 * Then the derivatives are calculated.
 * 
 * Finally the current values of the moments RHO, RHOV[XYZ] are put into RHO1, RHOV[XYZ]1 for future interpolation.
 * 
 * \param mpiGrid Grid
 * \param localCells Vector of local cells to process
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa calculateDerivatives
 */
void calculateDerivativesSimple(
       dccrg::Dccrg<SpatialCell>& mpiGrid,
       const vector<CellID>& localCells,
       cint& RKCase
) {
   int timer;
   namespace fs = fieldsolver;
   
   phiprof::start("Calculate derivatives");
   
   # ifndef FS_1ST_ORDER_TIME
   if(RKCase == RK_ORDER1) { // Means initialising the solver
      for (vector<uint64_t>::const_iterator cell = localCells.begin(); cell != localCells.end(); cell++) {
         if(mpiGrid[*cell]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         mpiGrid[*cell]->parameters[CellParams::RHO1] = mpiGrid[*cell]->parameters[CellParams::RHO];
         mpiGrid[*cell]->parameters[CellParams::RHOVX1] = mpiGrid[*cell]->parameters[CellParams::RHOVX];
         mpiGrid[*cell]->parameters[CellParams::RHOVY1] = mpiGrid[*cell]->parameters[CellParams::RHOVY];
         mpiGrid[*cell]->parameters[CellParams::RHOVZ1] = mpiGrid[*cell]->parameters[CellParams::RHOVZ];
      }
   }
   # endif
   if(RKCase == RK_ORDER2_STEP1) {
      // Interpolate linearly the moments.
      // The RHO(V?)1 fields contain previous value, the RHO(V?) the current one.
      // After this the RHO(V?)1 contain the interpolated value.
      for (vector<uint64_t>::const_iterator cell = localCells.begin(); cell != localCells.end(); cell++) {
         if(mpiGrid[*cell]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         mpiGrid[*cell]->parameters[CellParams::RHO1] = mpiGrid[*cell]->parameters[CellParams::RHO1] +
         Parameters::RK_alpha * (mpiGrid[*cell]->parameters[CellParams::RHO] - mpiGrid[*cell]->parameters[CellParams::RHO1]);
         mpiGrid[*cell]->parameters[CellParams::RHOVX1] = mpiGrid[*cell]->parameters[CellParams::RHOVX1] +
         Parameters::RK_alpha * (mpiGrid[*cell]->parameters[CellParams::RHOVX] - mpiGrid[*cell]->parameters[CellParams::RHOVX1]);
         mpiGrid[*cell]->parameters[CellParams::RHOVY1] = mpiGrid[*cell]->parameters[CellParams::RHOVY1] +
         Parameters::RK_alpha * (mpiGrid[*cell]->parameters[CellParams::RHOVY] - mpiGrid[*cell]->parameters[CellParams::RHOVY1]);
         mpiGrid[*cell]->parameters[CellParams::RHOVZ1] = mpiGrid[*cell]->parameters[CellParams::RHOVZ1] +
         Parameters::RK_alpha * (mpiGrid[*cell]->parameters[CellParams::RHOVZ] - mpiGrid[*cell]->parameters[CellParams::RHOVZ1]);
      }
   }
   
   if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ,BGBX,BGBY,BGBZ,RHO,RHOVX,RHOVY,RHOVZ with neighbours
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB_RHO_RHOV);
   } else { // RKCase == RK_ORDER2_STEP1        
      // Exchange BX1,BY1,BZ1,RHO1,RHOVX1,RHOVY1,RHOVZ1 with neighbours
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB1_RHO1_RHOV1);
   }
   
   timer=phiprof::initializeTimer("Start comm of PERB, BGB  and RHOV","MPI");
   phiprof::start(timer);
   mpiGrid.start_remote_neighbor_data_update();
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Compute inner cells");
   phiprof::start(timer);
   // Calculate derivatives on inner cells
   const vector<uint64_t> local_cells = mpiGrid.get_cells_with_local_neighbors();
   for (vector<uint64_t>::const_iterator cell = local_cells.begin(); cell != local_cells.end(); cell++) {
      if(mpiGrid[*cell]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      calculateDerivatives(*cell,mpiGrid, RKCase);
   }
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_neighbor_data_update_receives();
   phiprof::stop(timer);
   
   // Calculate derivatives on boundary cells
   timer=phiprof::initializeTimer("Compute boundary cells");
   phiprof::start(timer);
   const vector<uint64_t> boundary_cells = mpiGrid.get_cells_with_remote_neighbor();
   for (vector<uint64_t>::const_iterator cell = boundary_cells.begin(); cell != boundary_cells.end(); cell++) {
      if(mpiGrid[*cell]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      calculateDerivatives(*cell,mpiGrid, RKCase);
   }
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_neighbor_data_update_sends();
   phiprof::stop(timer);
   
   if(RKCase == RK_ORDER2_STEP2) { // Shift down the moments, now the RHO(V?)1 fields contain the current ie future previous values.
      for (vector<uint64_t>::const_iterator cell = localCells.begin(); cell != localCells.end(); cell++) {
         if(mpiGrid[*cell]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         mpiGrid[*cell]->parameters[CellParams::RHO1] = mpiGrid[*cell]->parameters[CellParams::RHO];
         mpiGrid[*cell]->parameters[CellParams::RHOVX1] = mpiGrid[*cell]->parameters[CellParams::RHOVX];
         mpiGrid[*cell]->parameters[CellParams::RHOVY1] = mpiGrid[*cell]->parameters[CellParams::RHOVY];
         mpiGrid[*cell]->parameters[CellParams::RHOVZ1] = mpiGrid[*cell]->parameters[CellParams::RHOVZ];
      }
   }
   
   phiprof::stop("Calculate derivatives");
}

/*! \brief High-level electric field computation function.
 * 
 * Transfers the derivatives, calculates the edge electric fields and transfers the new electric fields.
 * 
 * \param mpiGrid Grid
 * \param localCells Vector of local cells to process
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa calculateEdgeElectricFieldX calculateEdgeElectricFieldY calculateEdgeElectricFieldZ
 */
void calculateUpwindedElectricFieldSimple(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   const vector<CellID>& localCells,
   cint& RKCase
) {
   namespace fs = fieldsolver;
   int timer;
   phiprof::start("Calculate upwinded electric field");
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_DERIVATIVES);
   
   timer=phiprof::initializeTimer("Start communication of derivatives","MPI");
   phiprof::start(timer);
   mpiGrid.start_remote_neighbor_data_update();
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Compute inner cells");
   phiprof::start(timer);
   // Calculate upwinded electric field on inner cells
   const vector<uint64_t> local_cells = mpiGrid.get_cells_with_local_neighbors();
   for (vector<uint64_t>::const_iterator cell = local_cells.begin(); cell != local_cells.end(); cell++) {
      if(mpiGrid[*cell]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      cuint sysBoundaryFlag = sysBoundaryFlags[*cell];
      if ((sysBoundaryFlag & CALCULATE_EX) == CALCULATE_EX) calculateEdgeElectricFieldX(*cell,mpiGrid, RKCase);
      if ((sysBoundaryFlag & CALCULATE_EY) == CALCULATE_EY) calculateEdgeElectricFieldY(*cell,mpiGrid, RKCase);
      if ((sysBoundaryFlag & CALCULATE_EZ) == CALCULATE_EZ) calculateEdgeElectricFieldZ(*cell,mpiGrid, RKCase);
   }
   phiprof::stop(timer);
   timer=phiprof::initializeTimer("Wait for receives","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_neighbor_data_update_receives();
   phiprof::stop(timer);
   timer=phiprof::initializeTimer("Compute boundary cells");
   phiprof::start(timer);
   // Calculate upwinded electric field on boundary cells:
   const vector<uint64_t> boundary_cells = mpiGrid.get_cells_with_remote_neighbor();
   for (vector<uint64_t>::const_iterator cell = boundary_cells.begin(); cell != boundary_cells.end(); cell++) {
      if(mpiGrid[*cell]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      cuint sysBoundaryFlag = sysBoundaryFlags[*cell];
      if ((sysBoundaryFlag & CALCULATE_EX) == CALCULATE_EX) calculateEdgeElectricFieldX(*cell,mpiGrid, RKCase);
      if ((sysBoundaryFlag & CALCULATE_EY) == CALCULATE_EY) calculateEdgeElectricFieldY(*cell,mpiGrid, RKCase);
      if ((sysBoundaryFlag & CALCULATE_EZ) == CALCULATE_EZ) calculateEdgeElectricFieldZ(*cell,mpiGrid, RKCase);
   }
   phiprof::stop(timer);
   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_neighbor_data_update_sends();
   phiprof::stop(timer);
   
   // Exchange electric field with neighbouring processes
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_E);
   } else { // RKCase == RK_ORDER2_STEP1
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_E1);
   }
   timer=phiprof::initializeTimer("Communicate electric fields","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.update_remote_neighbor_data();
   phiprof::stop(timer);

   phiprof::stop("Calculate upwinded electric field");
}

/*! \brief High-level magnetic field propagation function.
 * 
 * Propagates the magnetic field and applies the field boundary conditions defined in project.h where needed.
 * 
 * \param mpiGrid Grid
 * \param dt Length of the time step
 * \param localCells Vector of local cells to process
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa propagateMagneticField
 */
static void propagateMagneticFieldSimple(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   creal& dt,
   const vector<CellID>& localCells,
   cint& RKCase
) {
   phiprof::start("Propagate magnetic field");
   // Propagate B on all local cells:
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      propagateMagneticField(cellID,mpiGrid,dt,RKCase);
   }
   
   // Calculate new B on faces outside the simulation domain using boundary conditions.
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      #ifndef NDEBUG
         const map<CellID,uint>::const_iterator it = sysBoundaryFlags.find(cellID);
         if (it == sysBoundaryFlags.end()) {cerr << "ERROR Could not find boundary flag for cell #" << cellID << endl; exit(1);}
         cuint existingCells = it->second;
      #else
         cuint existingCells = sysBoundaryFlags[cellID];
      #endif
      cuint nonExistingCells = (existingCells ^ numeric_limits<uint>::max());
      if (mpiGrid[cellID] == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " No data for cell " << cellID
            << std::endl;
         abort();
      }

      if ((existingCells & PROPAGATE_BX) != PROPAGATE_BX) {
         mpiGrid[cellID]->parameters[CellParams::PERBX] = fieldSolverBoundaryCondBx<CellID,uint,Real>(cellID,existingCells,nonExistingCells,mpiGrid);
         if(RKCase == RK_ORDER2_STEP1) { // WARNING not correct for time-varying boundary conditions
            mpiGrid[cellID]->parameters[CellParams::PERBX1] = mpiGrid[cellID]->parameters[CellParams::PERBX];
         }
      }
      if ((existingCells & PROPAGATE_BY) != PROPAGATE_BY) {
         mpiGrid[cellID]->parameters[CellParams::PERBY] = fieldSolverBoundaryCondBy<CellID,uint,Real>(cellID,existingCells,nonExistingCells,mpiGrid);
         if(RKCase == RK_ORDER2_STEP1) { // WARNING not correct for time-varying boundary conditions
            mpiGrid[cellID]->parameters[CellParams::PERBY1] = mpiGrid[cellID]->parameters[CellParams::PERBY];
         }
      }
      if ((existingCells & PROPAGATE_BZ) != PROPAGATE_BZ) {
         mpiGrid[cellID]->parameters[CellParams::PERBZ] = fieldSolverBoundaryCondBz<CellID,uint,Real>(cellID,existingCells,nonExistingCells,mpiGrid);
         if(RKCase == RK_ORDER2_STEP1) { // WARNING not correct for time-varying boundary conditions
            mpiGrid[cellID]->parameters[CellParams::PERBZ1] = mpiGrid[cellID]->parameters[CellParams::PERBZ];
         }
      }
   }
   phiprof::stop("Propagate magnetic field");
}

/*! \brief Top-level field propagation function.
 * 
 * Propagates the magnetic field, computes the derivatives and the upwinded electric field, then computes the volume-averaged field values. Takes care of the Runge-Kutta iteration at the top level, the functions called get as an argument the element from the enum defining the current stage and handle their job correspondingly.
 * 
 * \param mpiGrid Grid
 * \param dt Length of the time step
 * 
 * \sa propagateMagneticFieldSimple calculateDerivativesSimple calculateUpwindedElectricFieldSimple calculateVolumeAveragedFields
 * 
 */
bool propagateFields(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   creal& dt
) {
   // Reserve memory for derivatives for all cells on this process:
   vector<CellID> localCells = mpiGrid.get_cells();
   
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];   
      mpiGrid[cellID]->parameters[CellParams::MAXFDT]=std::numeric_limits<Real>::max();;
   }
   # ifdef FS_1ST_ORDER_TIME
   propagateMagneticFieldSimple(mpiGrid, dt, localCells, RK_ORDER1);
   calculateDerivativesSimple(mpiGrid, localCells, RK_ORDER1);
   calculateUpwindedElectricFieldSimple(mpiGrid, localCells, RK_ORDER1);
   # else
   propagateMagneticFieldSimple(mpiGrid, dt, localCells, RK_ORDER2_STEP1);
   calculateDerivativesSimple(mpiGrid, localCells, RK_ORDER2_STEP1);
   calculateUpwindedElectricFieldSimple(mpiGrid, localCells, RK_ORDER2_STEP1);
   
   propagateMagneticFieldSimple(mpiGrid, dt, localCells, RK_ORDER2_STEP2);
   calculateDerivativesSimple(mpiGrid, localCells, RK_ORDER2_STEP2);
   calculateUpwindedElectricFieldSimple(mpiGrid, localCells, RK_ORDER2_STEP2);
   # endif
   calculateVolumeAveragedFields(mpiGrid);
   return true;
}

/*! Namespace encompassing the enum defining the list of reconstruction coefficients used in field component reconstructions.*/
namespace Rec {
   /*! Enum defining the list of reconstruction coefficients used in field component reconstructions.*/
   enum Rec {a_0,a_x,a_y,a_z,a_xx,a_xy,a_xz,
   b_0,b_x,b_y,b_z,b_yx,b_yy,b_yz,
   c_0,c_x,c_y,c_z,c_zx,c_zy,c_zz
   };
}

/*! \brief Low-level helper function.
 *
 * Computes the reconstruction coefficients used for field component reconstruction.
 */
void reconstructionCoefficients(
   const CellID& cellID,
   const CellID& nbr_i2j1k1,
   const CellID& nbr_i1j2k1,
   const CellID& nbr_i1j1k2,
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   Real* result
) {
   // Do not calculate values for non-existing cells:
   if (cellID == INVALID_CELLID) {
      for (int i=0; i<Rec::c_zz+1; ++i) result[i] = 0.0;
      return;
   }
   
   namespace fs = fieldsolver;
   namespace cp = CellParams;
   
   Real* const cep_i1j1k1 = mpiGrid[cellID]->parameters;
   
   // Create a dummy array for containing zero values for cellParams on non-existing cells:
   Real dummyCellParams[CellParams::N_SPATIAL_CELL_PARAMS];
   for (uint i=0; i<CellParams::N_SPATIAL_CELL_PARAMS; ++i) dummyCellParams[i] = 0.0;
   
   Real* cep_i2j1k1 = NULL;
   Real* cep_i1j2k1 = NULL;
   Real* cep_i1j1k2 = NULL;
   if (nbr_i2j1k1 == INVALID_CELLID) cep_i2j1k1 = dummyCellParams;
   else cep_i2j1k1 = mpiGrid[nbr_i2j1k1]->parameters;
   if (nbr_i1j2k1 == INVALID_CELLID) cep_i1j2k1 = dummyCellParams;
   else cep_i1j2k1 = mpiGrid[nbr_i1j2k1]->parameters;
   if (nbr_i1j1k2 == INVALID_CELLID) cep_i1j1k2 = dummyCellParams;
   else cep_i1j1k2 = mpiGrid[nbr_i1j1k2]->parameters;

   #ifndef FS_1ST_ORDER_SPACE
      creal* const der_i1j1k1 = mpiGrid[cellID]->derivatives;

      // Create a dummy array for containing zero values for derivatives on non-existing cells:
      Real dummyDerivatives[N_SPATIAL_CELL_DERIVATIVES];
      for (uint i=0; i<N_SPATIAL_CELL_DERIVATIVES; ++i) dummyDerivatives[i] = 0.0;
   
      // Fetch neighbour cell derivatives, or in case the neighbour does not 
      // exist, use dummyDerivatives array:
      Real* der_i2j1k1 = NULL;
      Real* der_i1j2k1 = NULL;
      Real* der_i1j1k2 = NULL;
      if (nbr_i2j1k1 == INVALID_CELLID) der_i2j1k1 = dummyDerivatives;
      else der_i2j1k1 = mpiGrid[nbr_i2j1k1]->derivatives;
      if (nbr_i1j2k1 == INVALID_CELLID) der_i1j2k1 = dummyDerivatives;
      else der_i1j2k1 = mpiGrid[nbr_i1j2k1]->derivatives;
      if (nbr_i1j1k2 == INVALID_CELLID) der_i1j1k2 = dummyDerivatives;
      else der_i1j1k2 = mpiGrid[nbr_i1j1k2]->derivatives;

      // Calculate 2nd order reconstruction coefficients:
      result[Rec::a_xy] = der_i2j1k1[fs::dBxdy] - der_i1j1k1[fs::dBxdy];
      CHECK_FLOAT(result[Rec::a_xy])
      result[Rec::a_xz] = der_i2j1k1[fs::dBxdz] - der_i1j1k1[fs::dBxdz];
      CHECK_FLOAT(result[Rec::a_xz])
      result[Rec::a_x ] = (cep_i2j1k1[cp::PERBX] + cep_i2j1k1[cp::BGBX]) - (cep_i1j1k1[cp::PERBX] + cep_i1j1k1[cp::BGBX]);
      CHECK_FLOAT(result[Rec::a_x ])
      result[Rec::a_y ] = HALF*(der_i2j1k1[fs::dBxdy] + der_i1j1k1[fs::dBxdy]);
      CHECK_FLOAT(result[Rec::a_y ])
      result[Rec::a_z ] = HALF*(der_i2j1k1[fs::dBxdz] + der_i1j1k1[fs::dBxdz]);
      CHECK_FLOAT(result[Rec::a_z ])
   
      result[Rec::b_yx] = der_i1j2k1[fs::dBydx] - der_i1j1k1[fs::dBydx];
      CHECK_FLOAT(result[Rec::b_yx])
      result[Rec::b_yz] = der_i1j2k1[fs::dBydz] - der_i1j1k1[fs::dBydz];
      CHECK_FLOAT(result[Rec::b_yz])
      result[Rec::b_x ] = HALF*(der_i1j2k1[fs::dBydx] + der_i1j1k1[fs::dBydx]);
      CHECK_FLOAT(result[Rec::b_x ])
      result[Rec::b_y ] = (cep_i1j2k1[cp::PERBY] + cep_i1j2k1[cp::BGBY]) - (cep_i1j1k1[cp::PERBY] + cep_i1j1k1[cp::BGBY]);
      CHECK_FLOAT(result[Rec::b_y ])
      result[Rec::b_z ] = HALF*(der_i1j2k1[fs::dBydz] + der_i1j1k1[fs::dBydz]);
      CHECK_FLOAT(result[Rec::b_z ])
   
      result[Rec::c_zx] = der_i1j1k2[fs::dBzdx] - der_i1j1k1[fs::dBzdx];
      CHECK_FLOAT(result[Rec::c_zx])
      result[Rec::c_zy] = der_i1j1k2[fs::dBzdy] - der_i1j1k1[fs::dBzdy];
      CHECK_FLOAT(result[Rec::c_zy])
      result[Rec::c_x ] = HALF*(der_i1j1k2[fs::dBzdx] + der_i1j1k1[fs::dBzdx]);
      CHECK_FLOAT(result[Rec::c_x ])
      result[Rec::c_y ] = HALF*(der_i1j1k2[fs::dBzdy] + der_i1j1k1[fs::dBzdy]);
      CHECK_FLOAT(result[Rec::c_y ])
      result[Rec::c_z ] = (cep_i1j1k2[cp::PERBZ] + cep_i1j1k2[cp::BGBZ]) - (cep_i1j1k1[cp::PERBZ] + cep_i1j1k1[cp::BGBZ]);
      CHECK_FLOAT(result[Rec::c_z ])
   
      result[Rec::a_xx] = -HALF*(result[Rec::b_yx] + result[Rec::c_zx]);
      CHECK_FLOAT(result[Rec::a_xx])
      result[Rec::b_yy] = -HALF*(result[Rec::a_xy] + result[Rec::c_zy]);
      CHECK_FLOAT(result[Rec::b_yy])
      result[Rec::c_zz] = -HALF*(result[Rec::a_xz] + result[Rec::b_yz]);
      CHECK_FLOAT(result[Rec::c_zz])
   #else
      for (int i=0; i<Rec::c_zz+1; ++i) result[i] = 0.0;
   #endif
   
   // Calculate 1st order reconstruction coefficients:
   result[Rec::a_0 ] = HALF*(cep_i2j1k1[cp::PERBX] + cep_i2j1k1[cp::BGBX] + cep_i1j1k1[cp::PERBX] + cep_i1j1k1[cp::BGBX]) - SIXTH*result[Rec::a_xx];
   CHECK_FLOAT(result[Rec::a_0 ])
   result[Rec::b_0 ] = HALF*(cep_i1j2k1[cp::PERBY] + cep_i1j2k1[cp::BGBY] + cep_i1j1k1[cp::PERBY] + cep_i1j1k1[cp::BGBY]) - SIXTH*result[Rec::b_yy];
   CHECK_FLOAT(result[Rec::b_0 ])
   result[Rec::c_0 ] = HALF*(cep_i1j1k2[cp::PERBZ] + cep_i1j1k2[cp::BGBZ] + cep_i1j1k1[cp::PERBZ] + cep_i1j1k1[cp::BGBZ]) - SIXTH*result[Rec::c_zz];
   CHECK_FLOAT(result[Rec::c_0 ])
}

/*! \brief Top-level field averaging function.
 * 
 * Averages the electric and magnetic fields over the cell volumes.
 * 
 * \sa reconstructionCoefficients
 */
void calculateVolumeAveragedFields(
   dccrg::Dccrg<SpatialCell>& mpiGrid
) {
   phiprof::start("Calculate volume averaged fields");
   
   namespace fs = fieldsolver;
   namespace cp = CellParams;
                                   
   vector<uint64_t> localCells = mpiGrid.get_cells();
   
   Real coefficients[Rec::c_zz+1];
      
   cuint EX_CELLS = (1 << calcNbrNumber(1,1,1))
      | (1 << calcNbrNumber(1,2,1))
      | (1 << calcNbrNumber(1,1,2))
      | (1 << calcNbrNumber(1,2,2));
   cuint EY_CELLS = (1 << calcNbrNumber(1,1,1))
      | (1 << calcNbrNumber(2,1,1))
      | (1 << calcNbrNumber(1,1,2))
      | (1 << calcNbrNumber(2,1,2));
   cuint EZ_CELLS = (1 << calcNbrNumber(1,1,1))
      | (1 << calcNbrNumber(2,1,1))
      | (1 << calcNbrNumber(1,2,1))
      | (1 << calcNbrNumber(2,2,1));
   
   uint existingCells = 0;
   for (size_t cell=0; cell<localCells.size(); ++cell) {
      const CellID cellID = localCells[cell];
      
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
      
      // Get neighbour flags for the cell:
      map<CellID,uint>::const_iterator it = sysBoundaryFlags.find(cellID);
      if (it == sysBoundaryFlags.end()) existingCells = 0;
      else existingCells = it->second;
      
      // Calculate reconstruction coefficients for this cell:
      const CellID nbr_i2j1k1 = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2  );
      const CellID nbr_i1j2k1 = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2  );
      const CellID nbr_i1j1k2 = getNeighbourID(mpiGrid, cellID, 2  , 2  , 2+1);
      reconstructionCoefficients(cellID,nbr_i2j1k1,nbr_i1j2k1,nbr_i1j1k2,mpiGrid,coefficients);
      
      // Calculate volume average of B:
      Real* const cellParams = mpiGrid[cellID]->parameters;
      cellParams[cp::BXVOL] = coefficients[Rec::a_0]; //these include both background and 
      cellParams[cp::BYVOL] = coefficients[Rec::b_0];
      cellParams[cp::BZVOL] = coefficients[Rec::c_0];
      
      // Calculate volume average of E (NEEDS IMPROVEMENT):
      const CellID nbr_i1j2k2 = getNeighbourID(mpiGrid, cellID, 2  , 2+1, 2+1);
      const CellID nbr_i2j1k2 = getNeighbourID(mpiGrid, cellID, 2+1, 2  , 2+1);
      const CellID nbr_i2j2k1 = getNeighbourID(mpiGrid, cellID, 2+1, 2+1, 2  );
      creal* const cep_i1j1k1 = cellParams;
      
      if ((existingCells & EX_CELLS) == EX_CELLS) {
         creal* const cep_i1j2k1 = mpiGrid[nbr_i1j2k1]->parameters;
         creal* const cep_i1j1k2 = mpiGrid[nbr_i1j1k2]->parameters;
         creal* const cep_i1j2k2 = mpiGrid[nbr_i1j2k2]->parameters;

         CHECK_FLOAT(cep_i1j1k1[cp::EX])
         CHECK_FLOAT(cep_i1j2k1[cp::EX])
         CHECK_FLOAT(cep_i1j1k2[cp::EX])
         CHECK_FLOAT(cep_i1j2k2[cp::EX])
         cellParams[cp::EXVOL] = FOURTH*(cep_i1j1k1[cp::EX] + cep_i1j2k1[cp::EX] + cep_i1j1k2[cp::EX] + cep_i1j2k2[cp::EX]);
         CHECK_FLOAT(cellParams[cp::EXVOL])
      } else {
         cellParams[cp::EXVOL] = 0.0;
      }
      
      if ((existingCells & EY_CELLS) == EY_CELLS) {
         creal* const cep_i2j1k1 = mpiGrid[nbr_i2j1k1]->parameters;
         creal* const cep_i1j1k2 = mpiGrid[nbr_i1j1k2]->parameters;
         creal* const cep_i2j1k2 = mpiGrid[nbr_i2j1k2]->parameters;

         CHECK_FLOAT(cep_i1j1k1[cp::EY])
         CHECK_FLOAT(cep_i2j1k1[cp::EY])
         CHECK_FLOAT(cep_i1j1k2[cp::EY])
         CHECK_FLOAT(cep_i2j1k2[cp::EY])
         cellParams[cp::EYVOL] = FOURTH*(cep_i1j1k1[cp::EY] + cep_i2j1k1[cp::EY] + cep_i1j1k2[cp::EY] + cep_i2j1k2[cp::EY]);
         CHECK_FLOAT(cellParams[cp::EYVOL])
      } else {
         cellParams[cp::EYVOL] = 0.0;
      }
      
      if ((existingCells & EZ_CELLS) == EZ_CELLS) {
         creal* const cep_i2j1k1 = mpiGrid[nbr_i2j1k1]->parameters;
         creal* const cep_i1j2k1 = mpiGrid[nbr_i1j2k1]->parameters;
         creal* const cep_i2j2k1 = mpiGrid[nbr_i2j2k1]->parameters;

         CHECK_FLOAT(cep_i1j1k1[cp::EZ])
         CHECK_FLOAT(cep_i2j1k1[cp::EZ])
         CHECK_FLOAT(cep_i1j2k1[cp::EZ])
         CHECK_FLOAT(cep_i2j2k1[cp::EZ])
         cellParams[cp::EZVOL] = FOURTH*(cep_i1j1k1[cp::EZ] + cep_i2j1k1[cp::EZ] + cep_i1j2k1[cp::EZ] + cep_i2j2k1[cp::EZ]);
         CHECK_FLOAT(cellParams[cp::EZVOL])
      } else {
         cellParams[cp::EZVOL] = 0.0;
      }
   }
   phiprof::stop("Calculate volume averaged fields");
}
