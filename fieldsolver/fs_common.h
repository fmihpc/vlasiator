/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef FS_COMMON_H
#define FS_COMMON_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <list>
#include <set>
#include <stdint.h>

#include "../common.h"
#include "../grid.h"
#include "../parameters.h"
#include "../projects/project.h"
#include "phiprof.hpp"
#include "../sysboundary/sysboundary.h"
#include "../sysboundary/sysboundarycondition.h"

// Constants: not needed as such, but if field solver is implemented on GPUs 
// these force CPU to use float accuracy, which in turn helps to compare 
// CPU and GPU results.

const Real HALF    = 0.5;
const Real MINUS   = -1.0;
const Real PLUS    = +1.0;
const Real THIRD   = 1.0/3.0;
const Real FOURTH  = 1.0/4.0;
const Real SIXTH   = 1.0/6.0;
const Real EIGTH   = 1.0/8.0;
const Real TENTH   = 1.0/10.0;
const Real TWELWTH = 1.0/12.0;
const Real TWO     = 2.0;
const Real ZERO    = 0.0;

static creal EPS = 1.0e-30;

using namespace std;
using namespace fieldsolver;

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
static uint CALCULATE_DXY = 0; /**< Bit mask determining if xy mixed derivatives can be calculated on a cell.*/
static uint CALCULATE_DXZ = 0; /**< Bit mask determining if xz mixed derivatives can be calculated on a cell.*/
static uint CALCULATE_DYZ = 0; /**< Bit mask determining if yz mixed derivatives can be calculated on a cell.*/
static uint CALCULATE_EX = 0; /**< Bit mask determining if edge Ex can be calculated on a cell.*/
static uint CALCULATE_EY = 0; /**< Bit mask determining if edge Ey can be calculated on a cell.*/
static uint CALCULATE_EZ = 0; /**< Bit mask determining if edge Ez can be calculated on a cell.*/
static uint PROPAGATE_BX = 0; /**< Bit mask determining if face Bx is propagated on a cell.*/
static uint PROPAGATE_BY = 0; /**< Bit mask determining if face By is propagated on a cell.*/
static uint PROPAGATE_BZ = 0; /**< Bit mask determining if face Bz is propagated on a cell.*/

bool initializeFieldPropagator(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries
);
bool initializeFieldPropagatorAfterRebalance(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
bool finalizeFieldPropagator(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
bool propagateFields(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& subcycles
);

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

CellID getNeighbourID(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const CellID& cellID,
   const uchar& i,
   const uchar& j,
   const uchar& k
);

/*! \brief Low-level helper function.
 * 
 * Avoid crashes on zero density by returning V = rhoV = 0.0 if rho = 0.0.
 * 
 */
Real divideIfNonZero(creal rhoV, creal rho);

/*! Namespace encompassing the enum defining the list of reconstruction coefficients used in field component reconstructions.*/
namespace Rec {
   /*! Enum defining the list of reconstruction coefficients used in field component reconstructions.*/
   enum Rec {
      a_0, a_x, a_y, a_z, a_xx, a_yy, a_zz, a_xy, a_xz, a_yz, a_xxx, a_xxy, a_xyy, a_xxz, a_xzz, a_xyz,
      b_0, b_x, b_y, b_z, b_xx, b_yy, b_zz, b_xy, b_xz, b_yz, b_xxy, b_xyy, b_yyy, b_yyz, b_yzz, b_xyz,
      c_0, c_x, c_y, c_z, c_xx, c_yy, c_zz, c_xy, c_xz, c_yz, c_xxz, c_xzz, c_yyz, c_yzz, c_xyz, c_zzz,
      N_REC_COEFFICIENTS
   };
}

/*! \brief Low-level helper function.
 * 
 * Computes the reconstruction coefficients used for field component reconstruction.
 * Only implemented for 2nd and 3rd order.
 * 
 * \param reconstructionOrder Reconstruction order of the fields after Balsara 2009, 2 used for BVOL, 3 used for 2nd-order Hall term calculations.
 */
void reconstructionCoefficients(
   const CellID& cellID,
   const CellID& nbr_i2j1k1,
   const CellID& nbr_i1j2k1,
   const CellID& nbr_i1j1k2,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   Real* perturbedResult,
   creal& reconstructionOrder,
   cint& RKCase
);



#endif
