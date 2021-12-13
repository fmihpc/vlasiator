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

#include "fs_common.h"

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
 * \param perBGrid fsGrid holding the perturbed B quantities 
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param perturbedResult Array in which to store the coefficients.
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param reconstructionOrder Reconstruction order of the fields after Balsara 2009, 2 used for BVOL, 3 used for 2nd-order Hall term calculations.
 */
void reconstructionCoefficients(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   Real* perturbedResult,
   cint i,
   cint j,
   cint k,
   creal& reconstructionOrder
) {
   std::array<Real, fsgrids::bfield::N_BFIELD> * cep_i1j1k1 = NULL;
   std::array<Real, fsgrids::dperb::N_DPERB> * der_i1j1k1 = dPerBGrid.get(i,j,k);
   std::array<Real, fsgrids::bfield::N_BFIELD> * dummyCellParams = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD> * cep_i2j1k1 = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD> * cep_i1j2k1 = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD> * cep_i1j1k2 = NULL;
   
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> * params = & perBGrid;
   
   cep_i1j1k1 = params->get(i,j,k);
   dummyCellParams = cep_i1j1k1;
   cep_i2j1k1 = dummyCellParams;
   cep_i1j2k1 = dummyCellParams;
   cep_i1j1k2 = dummyCellParams;
   if (params->get(i+1,j,k) != NULL) cep_i2j1k1 = params->get(i+1,j,k);
   if (params->get(i,j+1,k) != NULL) cep_i1j2k1 = params->get(i,j+1,k);
   if (params->get(i,j,k+1) != NULL) cep_i1j1k2 = params->get(i,j,k+1);
   
   #ifndef FS_1ST_ORDER_SPACE

   // Create a dummy array for containing zero values for derivatives on non-existing cells:
   std::array<Real, fsgrids::dperb::N_DPERB> dummyDerivatives;
   for (int ii=0; ii<fsgrids::dperb::N_DPERB; ii++) {
      dummyDerivatives.at(ii) = 0.0;
   }
   
   // Fetch neighbour cell derivatives, or in case the neighbour does not 
   // exist, use dummyDerivatives array:
   std::array<Real, fsgrids::dperb::N_DPERB> * der_i2j1k1 = &dummyDerivatives;
   std::array<Real, fsgrids::dperb::N_DPERB> * der_i1j2k1 = &dummyDerivatives;
   std::array<Real, fsgrids::dperb::N_DPERB> * der_i1j1k2 = &dummyDerivatives;
   if (dPerBGrid.get(i+1,j,k) != NULL) der_i2j1k1 = dPerBGrid.get(i+1,j,k);
   if (dPerBGrid.get(i,j+1,k) != NULL) der_i1j2k1 = dPerBGrid.get(i,j+1,k);
   if (dPerBGrid.get(i,j,k+1) != NULL) der_i1j1k2 = dPerBGrid.get(i,j,k+1);
   
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
      perturbedResult[Rec::a_yy] = HALF * (der_i2j1k1->at(fsgrids::dperb::dPERBxdyy) + der_i1j1k1->at(fsgrids::dperb::dPERBxdyy));
      perturbedResult[Rec::a_zz] = HALF * (der_i2j1k1->at(fsgrids::dperb::dPERBxdzz) + der_i1j1k1->at(fsgrids::dperb::dPERBxdzz));
      perturbedResult[Rec::a_yz] = HALF * (der_i2j1k1->at(fsgrids::dperb::dPERBxdyz) + der_i1j1k1->at(fsgrids::dperb::dPERBxdyz));
      perturbedResult[Rec::a_xyy] = (der_i2j1k1->at(fsgrids::dperb::dPERBxdyy) - der_i1j1k1->at(fsgrids::dperb::dPERBxdyy));
      perturbedResult[Rec::a_xyz] = (der_i2j1k1->at(fsgrids::dperb::dPERBxdyz) - der_i1j1k1->at(fsgrids::dperb::dPERBxdyz));
      perturbedResult[Rec::a_xzz] = (der_i2j1k1->at(fsgrids::dperb::dPERBxdzz) - der_i1j1k1->at(fsgrids::dperb::dPERBxdzz));
      
      perturbedResult[Rec::b_xx] = HALF * (der_i1j2k1->at(fsgrids::dperb::dPERBydxx) + der_i1j1k1->at(fsgrids::dperb::dPERBydxx));
      perturbedResult[Rec::b_xz] = HALF * (der_i1j2k1->at(fsgrids::dperb::dPERBydxz) + der_i1j1k1->at(fsgrids::dperb::dPERBydxz));
      perturbedResult[Rec::b_zz] = HALF * (der_i1j2k1->at(fsgrids::dperb::dPERBydzz) + der_i1j1k1->at(fsgrids::dperb::dPERBydzz));
      perturbedResult[Rec::b_xxy] = (der_i1j2k1->at(fsgrids::dperb::dPERBydxx) - der_i1j1k1->at(fsgrids::dperb::dPERBydxx));
      perturbedResult[Rec::b_xyz] = (der_i1j2k1->at(fsgrids::dperb::dPERBydxz) - der_i1j1k1->at(fsgrids::dperb::dPERBydxz));
      perturbedResult[Rec::b_yzz] = (der_i1j2k1->at(fsgrids::dperb::dPERBydzz) - der_i1j1k1->at(fsgrids::dperb::dPERBydzz));
      
      perturbedResult[Rec::c_xx] = HALF * (der_i1j1k2->at(fsgrids::dperb::dPERBzdxx) + der_i1j1k1->at(fsgrids::dperb::dPERBzdxx));
      perturbedResult[Rec::c_xy] = HALF * (der_i1j1k2->at(fsgrids::dperb::dPERBzdxy) + der_i1j1k1->at(fsgrids::dperb::dPERBzdxy));
      perturbedResult[Rec::c_yy] = HALF * (der_i1j1k2->at(fsgrids::dperb::dPERBzdyy) + der_i1j1k1->at(fsgrids::dperb::dPERBzdyy));
      perturbedResult[Rec::c_xxz] = (der_i1j1k2->at(fsgrids::dperb::dPERBzdxx) - der_i1j1k1->at(fsgrids::dperb::dPERBzdxx));
      perturbedResult[Rec::c_xyz] = (der_i1j1k2->at(fsgrids::dperb::dPERBzdxy) - der_i1j1k1->at(fsgrids::dperb::dPERBzdxy));
      perturbedResult[Rec::c_yyz] = (der_i1j1k2->at(fsgrids::dperb::dPERBzdyy) - der_i1j1k1->at(fsgrids::dperb::dPERBzdyy));
      
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
   perturbedResult[Rec::a_xy] = der_i2j1k1->at(fsgrids::dperb::dPERBxdy) - der_i1j1k1->at(fsgrids::dperb::dPERBxdy);
   perturbedResult[Rec::a_xz] = der_i2j1k1->at(fsgrids::dperb::dPERBxdz) - der_i1j1k1->at(fsgrids::dperb::dPERBxdz);
   perturbedResult[Rec::a_y ] = HALF*(der_i2j1k1->at(fsgrids::dperb::dPERBxdy) + der_i1j1k1->at(fsgrids::dperb::dPERBxdy)) - SIXTH*perturbedResult[Rec::a_xxy];
   perturbedResult[Rec::a_z ] = HALF*(der_i2j1k1->at(fsgrids::dperb::dPERBxdz) + der_i1j1k1->at(fsgrids::dperb::dPERBxdz)) - SIXTH*perturbedResult[Rec::a_xxz];
   
   perturbedResult[Rec::b_xy] = der_i1j2k1->at(fsgrids::dperb::dPERBydx) - der_i1j1k1->at(fsgrids::dperb::dPERBydx);
   perturbedResult[Rec::b_yz] = der_i1j2k1->at(fsgrids::dperb::dPERBydz) - der_i1j1k1->at(fsgrids::dperb::dPERBydz);
   perturbedResult[Rec::b_x ] = HALF*(der_i1j2k1->at(fsgrids::dperb::dPERBydx) + der_i1j1k1->at(fsgrids::dperb::dPERBydx)) - SIXTH*perturbedResult[Rec::b_xyy];
   perturbedResult[Rec::b_z ] = HALF*(der_i1j2k1->at(fsgrids::dperb::dPERBydz) + der_i1j1k1->at(fsgrids::dperb::dPERBydz)) - SIXTH*perturbedResult[Rec::b_yyz];
   
   perturbedResult[Rec::c_xz] = der_i1j1k2->at(fsgrids::dperb::dPERBzdx) - der_i1j1k1->at(fsgrids::dperb::dPERBzdx);
   perturbedResult[Rec::c_yz] = der_i1j1k2->at(fsgrids::dperb::dPERBzdy) - der_i1j1k1->at(fsgrids::dperb::dPERBzdy);
   perturbedResult[Rec::c_x ] = HALF*(der_i1j1k2->at(fsgrids::dperb::dPERBzdx) + der_i1j1k1->at(fsgrids::dperb::dPERBzdx)) - SIXTH*perturbedResult[Rec::c_xzz];
   perturbedResult[Rec::c_y ] = HALF*(der_i1j1k2->at(fsgrids::dperb::dPERBzdy) + der_i1j1k1->at(fsgrids::dperb::dPERBzdy)) - SIXTH*perturbedResult[Rec::c_yzz];
   
   perturbedResult[Rec::a_xx] = -HALF*(perturbedResult[Rec::b_xy] + perturbedResult[Rec::c_xz]);
   perturbedResult[Rec::b_yy] = -HALF*(perturbedResult[Rec::a_xy] + perturbedResult[Rec::c_yz]);
   perturbedResult[Rec::c_zz] = -HALF*(perturbedResult[Rec::a_xz] + perturbedResult[Rec::b_yz]);
         
   perturbedResult[Rec::a_x ] = cep_i2j1k1->at(fsgrids::bfield::PERBX) - cep_i1j1k1->at(fsgrids::bfield::PERBX) - TENTH*perturbedResult[Rec::a_xxx];
   perturbedResult[Rec::b_y ] = cep_i1j2k1->at(fsgrids::bfield::PERBY) - cep_i1j1k1->at(fsgrids::bfield::PERBY) - TENTH*perturbedResult[Rec::b_yyy];
   perturbedResult[Rec::c_z ] = cep_i1j1k2->at(fsgrids::bfield::PERBZ) - cep_i1j1k1->at(fsgrids::bfield::PERBZ) - TENTH*perturbedResult[Rec::c_zzz];

   #else
   for (int i=0; i<Rec::N_REC_COEFFICIENTS; ++i) {
      perturbedResult[i] = 0.0;
   }
   #endif

   // Calculate 1st order reconstruction coefficients:
   perturbedResult[Rec::a_0 ] = HALF*(cep_i2j1k1->at(fsgrids::bfield::PERBX) + cep_i1j1k1->at(fsgrids::bfield::PERBX)) - SIXTH*perturbedResult[Rec::a_xx];
   perturbedResult[Rec::b_0 ] = HALF*(cep_i1j2k1->at(fsgrids::bfield::PERBY) + cep_i1j1k1->at(fsgrids::bfield::PERBY)) - SIXTH*perturbedResult[Rec::b_yy];
   perturbedResult[Rec::c_0 ] = HALF*(cep_i1j1k2->at(fsgrids::bfield::PERBZ) + cep_i1j1k1->at(fsgrids::bfield::PERBZ)) - SIXTH*perturbedResult[Rec::c_zz];
}

