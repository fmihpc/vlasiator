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

#include "fs_common.h"
#include "fs_cache.h"

map<CellID,uint> existingCellsFlags;


/*! \brief Helper function
 * 
 * Uses mpiGrid/dccrg to retrieve the cellID of a neighbour, taking its offset. Does not use the field solver cache, so it is not recommended to use it in intensive code sections.
 * 
 * \param mpiGrid Grid
 * \param cellID ID of the current cell
 * \param i Offset in x, 2 is current cell.
 * \param j Offset in y, 2 is current cell.
 * \param k Offset in z, 2 is current cell.
 */
CellID getNeighbourID(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const CellID& cellID,
   const uchar& i,
   const uchar& j,
   const uchar& k
) {
   
   //#ifdef DEBUG_SOLVERS
   // check that requested neighbor is within one index
   if (i < 1) {
      cerr << __FILE__ << ":" << __LINE__ << endl;
      abort();
   }
   if (i > 3) {
      cerr << __FILE__ << ":" << __LINE__ << endl;
      abort();
   }
   
   if (j < 1) {
      cerr << __FILE__ << ":" << __LINE__ << endl;
      abort();
   }
   if (j > 3) {
      cerr << __FILE__ << ":" << __LINE__ << endl;
      abort();
   }
   
   if (k < 1) {
      cerr << __FILE__ << ":" << __LINE__ << endl;
      abort();
   }
   if (k > 3) {
      cerr << __FILE__ << ":" << __LINE__ << endl;
      abort();
   }
   
   // NOTE: The following dccrg call is really slow..
   const std::vector<CellID> neighbors = mpiGrid.get_neighbors_of_at_offset(cellID, int(i) - 2, int(j) - 2, int(k) - 2);

   if (neighbors.size() == 0) {
      stringstream ss;
      ss << __FILE__ << ":" << __LINE__
      << " No neighbor for cell " << cellID
      << " at offsets " << int(i) - 2 << ", " << int(j) - 2 << ", " << int(k) - 2
      << endl;
      cerr << ss.str();
      abort();
   }
   // TODO support spatial refinement
   return neighbors[0];
}

/*! \brief Helper function
 * 
 * Divides the first value by the second or returns zero if the denominator is zero.
 * 
 * \param numerator Numerator
 * \param denominator Denominator
 */
Real divideIfNonZero(
   creal numerator,
   creal denominator
) {
   if(denominator <= 0.0) {
      return 0.0;
   } else {
      return numerator / denominator;
   }
}

/*! \brief Low-level helper function.
 * 
 * Computes the reconstruction coefficients used for field component reconstruction.
 * Only implemented for 2nd and 3rd order.
 * 
 * \param cell Field solver cell cache
 * \param perturbedResult Array in which to store the coefficients.
 * \param reconstructionOrder Reconstruction order of the fields after Balsara 2009, 2 used for BVOL, 3 used for 2nd-order Hall term calculations.
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 */
void reconstructionCoefficients(
   fs_cache::CellCache& cell,
   Real* perturbedResult,
   creal& reconstructionOrder,
   cint& RKCase
) {
   
   namespace fs = fieldsolver;
   namespace cp = CellParams;

   Real*  const cep_i1j1k1 = cell.cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
   creal* const der_i1j1k1 = cell.cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->derivatives;

   // Create a dummy array for containing zero values for cellParams on non-existing cells:
   Real dummyCellParams[CellParams::N_SPATIAL_CELL_PARAMS];
   for (uint i=0; i<CellParams::N_SPATIAL_CELL_PARAMS; ++i) dummyCellParams[i] = 0.0;

   creal* cep_i2j1k1 = dummyCellParams;
   creal* cep_i1j2k1 = dummyCellParams;
   creal* cep_i1j1k2 = dummyCellParams;
   if (cell.cells[fs_cache::calculateNbrID(1+1,1  ,1  )] != NULL) cep_i2j1k1 = cell.cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->parameters;
   if (cell.cells[fs_cache::calculateNbrID(1  ,1+1,1  )] != NULL) cep_i1j2k1 = cell.cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->parameters;
   if (cell.cells[fs_cache::calculateNbrID(1  ,1  ,1+1)] != NULL) cep_i1j1k2 = cell.cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->parameters;

   #ifndef FS_1ST_ORDER_SPACE

   // Create a dummy array for containing zero values for derivatives on non-existing cells:
   Real dummyDerivatives[N_SPATIAL_CELL_DERIVATIVES];
   for (uint i=0; i<N_SPATIAL_CELL_DERIVATIVES; ++i) dummyDerivatives[i] = 0.0;

   // Fetch neighbour cell derivatives, or in case the neighbour does not 
   // exist, use dummyDerivatives array:
   creal* der_i2j1k1 = dummyDerivatives;
   creal* der_i1j2k1 = dummyDerivatives;
   creal* der_i1j1k2 = dummyDerivatives;
   if (cell.cells[fs_cache::calculateNbrID(1+1,1  ,1  )] != NULL) der_i2j1k1 = cell.cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->derivatives;
   if (cell.cells[fs_cache::calculateNbrID(1  ,1+1,1  )] != NULL) der_i1j2k1 = cell.cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->derivatives;
   if (cell.cells[fs_cache::calculateNbrID(1  ,1  ,1+1)] != NULL) der_i1j1k2 = cell.cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->derivatives;

   // Calculate 3rd order reconstruction coefficients:
   if (reconstructionOrder == 2) {
      perturbedResult[Rec::a_yy] = 0.0;
      perturbedResult[Rec::a_zz] = 0.0;
      perturbedResult[Rec::a_yz] = 0.0;
      perturbedResult[Rec::a_xxx] = 0.0;
      perturbedResult[Rec::a_xxy] = 0.0;
      perturbedResult[Rec::a_xxz] = 0.0;
      perturbedResult[Rec::a_xyy] = 0.0;
      perturbedResult[Rec::a_xyz] = 0.0;
      perturbedResult[Rec::a_xzz] = 0.0;
      perturbedResult[Rec::b_xx] = 0.0;
      perturbedResult[Rec::b_xz] = 0.0;
      perturbedResult[Rec::b_zz] = 0.0;
      perturbedResult[Rec::b_xxy] = 0.0;
      perturbedResult[Rec::b_xyy] = 0.0;
      perturbedResult[Rec::b_xyz] = 0.0;
      perturbedResult[Rec::b_yyy] = 0.0;
      perturbedResult[Rec::b_yyz] = 0.0;
      perturbedResult[Rec::b_yzz] = 0.0;
      perturbedResult[Rec::c_xx] = 0.0;
      perturbedResult[Rec::c_xy] = 0.0;
      perturbedResult[Rec::c_yy] = 0.0;
      perturbedResult[Rec::c_xxz] = 0.0;
      perturbedResult[Rec::c_xyz] = 0.0;
      perturbedResult[Rec::c_xzz] = 0.0;
      perturbedResult[Rec::c_yyz] = 0.0;
      perturbedResult[Rec::c_yzz] = 0.0;
      perturbedResult[Rec::c_zzz] = 0.0;
   } else if (reconstructionOrder == 3) {
      perturbedResult[Rec::a_yy] = HALF * (der_i2j1k1[fs::dPERBxdyy] + der_i1j1k1[fs::dPERBxdyy]);
      perturbedResult[Rec::a_zz] = HALF * (der_i2j1k1[fs::dPERBxdzz] + der_i1j1k1[fs::dPERBxdzz]);
      perturbedResult[Rec::a_yz] = HALF * (der_i2j1k1[fs::dPERBxdyz] + der_i1j1k1[fs::dPERBxdyz]);
      perturbedResult[Rec::a_xyy] = (der_i2j1k1[fs::dPERBxdyy] - der_i1j1k1[fs::dPERBxdyy]);
      perturbedResult[Rec::a_xyz] = (der_i2j1k1[fs::dPERBxdyz] - der_i1j1k1[fs::dPERBxdyz]);
      perturbedResult[Rec::a_xzz] = (der_i2j1k1[fs::dPERBxdzz] - der_i1j1k1[fs::dPERBxdzz]);
      
      perturbedResult[Rec::b_xx] = HALF * (der_i1j2k1[fs::dPERBydxx] + der_i1j1k1[fs::dPERBydxx]);
      perturbedResult[Rec::b_xz] = HALF * (der_i1j2k1[fs::dPERBydxz] + der_i1j1k1[fs::dPERBydxz]);
      perturbedResult[Rec::b_zz] = HALF * (der_i1j2k1[fs::dPERBydzz] + der_i1j1k1[fs::dPERBydzz]);
      perturbedResult[Rec::b_xxy] = (der_i1j2k1[fs::dPERBydxx] - der_i1j1k1[fs::dPERBydxx]);
      perturbedResult[Rec::b_xyz] = (der_i1j2k1[fs::dPERBydxz] - der_i1j1k1[fs::dPERBydxz]);
      perturbedResult[Rec::b_yzz] = (der_i1j2k1[fs::dPERBydzz] - der_i1j1k1[fs::dPERBydzz]);
      
      perturbedResult[Rec::c_xx] = HALF * (der_i1j1k2[fs::dPERBzdxx] + der_i1j1k1[fs::dPERBzdxx]);
      perturbedResult[Rec::c_xy] = HALF * (der_i1j1k2[fs::dPERBzdxy] + der_i1j1k1[fs::dPERBzdxy]);
      perturbedResult[Rec::c_yy] = HALF * (der_i1j1k2[fs::dPERBzdyy] + der_i1j1k1[fs::dPERBzdyy]);
      perturbedResult[Rec::c_xxz] = (der_i1j1k2[fs::dPERBzdxx] - der_i1j1k1[fs::dPERBzdxx]);
      perturbedResult[Rec::c_xyz] = (der_i1j1k2[fs::dPERBzdxy] - der_i1j1k1[fs::dPERBzdxy]);
      perturbedResult[Rec::c_yyz] = (der_i1j1k2[fs::dPERBzdyy] - der_i1j1k1[fs::dPERBzdyy]);
      
      perturbedResult[Rec::a_xxx] = -THIRD*(perturbedResult[Rec::b_xxy] + perturbedResult[Rec::c_xxz]);
      perturbedResult[Rec::a_xxy] = -FOURTH*perturbedResult[Rec::c_xyz];
      perturbedResult[Rec::a_xxz] = -FOURTH*perturbedResult[Rec::b_xyz];
      
      perturbedResult[Rec::b_xyy] = -FOURTH*perturbedResult[Rec::c_xyz];
      perturbedResult[Rec::b_yyy] = -THIRD*(perturbedResult[Rec::c_yyz] + perturbedResult[Rec::a_xyy]);
      perturbedResult[Rec::b_yyz] = -FOURTH*perturbedResult[Rec::a_xyz];
      
      perturbedResult[Rec::c_xzz] = -FOURTH*perturbedResult[Rec::b_xyz];
      perturbedResult[Rec::c_yzz] = -FOURTH*perturbedResult[Rec::a_xyz];
      perturbedResult[Rec::c_zzz] = -THIRD*(perturbedResult[Rec::a_xzz] + perturbedResult[Rec::b_yzz]);
   } else {
      cerr << __FILE__ << ":" << __LINE__ << ":" << " Not coded yet!" << endl;
      abort();
   }
   
   // Calculate 2nd order reconstruction coefficients:
   perturbedResult[Rec::a_xy] = der_i2j1k1[fs::dPERBxdy] - der_i1j1k1[fs::dPERBxdy];
   perturbedResult[Rec::a_xz] = der_i2j1k1[fs::dPERBxdz] - der_i1j1k1[fs::dPERBxdz];
   perturbedResult[Rec::a_y ] = HALF*(der_i2j1k1[fs::dPERBxdy] + der_i1j1k1[fs::dPERBxdy]) - SIXTH*perturbedResult[Rec::a_xxy];
   perturbedResult[Rec::a_z ] = HALF*(der_i2j1k1[fs::dPERBxdz] + der_i1j1k1[fs::dPERBxdz]) - SIXTH*perturbedResult[Rec::a_xxz];
   
   perturbedResult[Rec::b_xy] = der_i1j2k1[fs::dPERBydx] - der_i1j1k1[fs::dPERBydx];
   perturbedResult[Rec::b_yz] = der_i1j2k1[fs::dPERBydz] - der_i1j1k1[fs::dPERBydz];
   perturbedResult[Rec::b_x ] = HALF*(der_i1j2k1[fs::dPERBydx] + der_i1j1k1[fs::dPERBydx]) - SIXTH*perturbedResult[Rec::b_xyy];
   perturbedResult[Rec::b_z ] = HALF*(der_i1j2k1[fs::dPERBydz] + der_i1j1k1[fs::dPERBydz]) - SIXTH*perturbedResult[Rec::b_yyz];
   
   perturbedResult[Rec::c_xz] = der_i1j1k2[fs::dPERBzdx] - der_i1j1k1[fs::dPERBzdx];
   perturbedResult[Rec::c_yz] = der_i1j1k2[fs::dPERBzdy] - der_i1j1k1[fs::dPERBzdy];
   perturbedResult[Rec::c_x ] = HALF*(der_i1j1k2[fs::dPERBzdx] + der_i1j1k1[fs::dPERBzdx]) - SIXTH*perturbedResult[Rec::c_xzz];
   perturbedResult[Rec::c_y ] = HALF*(der_i1j1k2[fs::dPERBzdy] + der_i1j1k1[fs::dPERBzdy]) - SIXTH*perturbedResult[Rec::c_yzz];
   
   perturbedResult[Rec::a_xx] = -HALF*(perturbedResult[Rec::b_xy] + perturbedResult[Rec::c_xz]);
   perturbedResult[Rec::b_yy] = -HALF*(perturbedResult[Rec::a_xy] + perturbedResult[Rec::c_yz]);
   perturbedResult[Rec::c_zz] = -HALF*(perturbedResult[Rec::a_xz] + perturbedResult[Rec::b_yz]);
         
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      perturbedResult[Rec::a_x ] = cep_i2j1k1[cp::PERBX] - cep_i1j1k1[cp::PERBX] - TENTH*perturbedResult[Rec::a_xxx];
      perturbedResult[Rec::b_y ] = cep_i1j2k1[cp::PERBY] - cep_i1j1k1[cp::PERBY] - TENTH*perturbedResult[Rec::b_yyy];
      perturbedResult[Rec::c_z ] = cep_i1j1k2[cp::PERBZ] - cep_i1j1k1[cp::PERBZ] - TENTH*perturbedResult[Rec::c_zzz];
   }
   if (RKCase == RK_ORDER2_STEP1) {
      perturbedResult[Rec::a_x ] = cep_i2j1k1[cp::PERBX_DT2] - cep_i1j1k1[cp::PERBX_DT2] - TENTH*perturbedResult[Rec::a_xxx];
      perturbedResult[Rec::b_y ] = cep_i1j2k1[cp::PERBY_DT2] - cep_i1j1k1[cp::PERBY_DT2] - TENTH*perturbedResult[Rec::b_yyy];
      perturbedResult[Rec::c_z ] = cep_i1j1k2[cp::PERBZ_DT2] - cep_i1j1k1[cp::PERBZ_DT2] - TENTH*perturbedResult[Rec::c_zzz];
   }

   #else
   for (int i=0; i<Rec::N_REC_COEFFICIENTS; ++i) {
      perturbedResult[i] = 0.0;
   }
   #endif

   // Calculate 1st order reconstruction coefficients:
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      perturbedResult[Rec::a_0 ] = HALF*(cep_i2j1k1[cp::PERBX] + cep_i1j1k1[cp::PERBX]) - SIXTH*perturbedResult[Rec::a_xx];
      perturbedResult[Rec::b_0 ] = HALF*(cep_i1j2k1[cp::PERBY] + cep_i1j1k1[cp::PERBY]) - SIXTH*perturbedResult[Rec::b_yy];
      perturbedResult[Rec::c_0 ] = HALF*(cep_i1j1k2[cp::PERBZ] + cep_i1j1k1[cp::PERBZ]) - SIXTH*perturbedResult[Rec::c_zz];
   }
   if (RKCase == RK_ORDER2_STEP1) {
      perturbedResult[Rec::a_0 ] = HALF*(cep_i2j1k1[cp::PERBX_DT2] + cep_i1j1k1[cp::PERBX_DT2]) - SIXTH*perturbedResult[Rec::a_xx];
      perturbedResult[Rec::b_0 ] = HALF*(cep_i1j2k1[cp::PERBY_DT2] + cep_i1j1k1[cp::PERBY_DT2]) - SIXTH*perturbedResult[Rec::b_yy];
      perturbedResult[Rec::c_0 ] = HALF*(cep_i1j1k2[cp::PERBZ_DT2] + cep_i1j1k1[cp::PERBZ_DT2]) - SIXTH*perturbedResult[Rec::c_zz];
   }
}

