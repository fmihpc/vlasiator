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
#include "ldz_hall.hpp"

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

using namespace std;

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBY Background By
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermXComponents
 *
 */
template<typename REAL> inline
REAL JXBX_000_100(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> pC,
   creal BGBY,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz-(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)+(pC[c_xz]*BGBZ)/(2*dx)-(pC[c_xyz]*BGBZ)/(4*dx)+(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx-(pC[a_yz]*BGBY)/(2*dy)-(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy+(pC[b_xz]*BGBY)/(2*dx)-(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)+(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx-(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)-
     (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)-(pC[a_z]*pC[c_z])/(2*dz)+(pC[a_yz]*pC[c_z])/(4*dz)-(pC[a_zz]*pC[c_yz])/(4*dz)+(pC[a_z]*pC[c_yz])/(4*dz)-(pC[a_yz]*pC[c_yz])/(8*dz)+(pC[a_zz]*pC[c_y])/(2*dz)-(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)-(pC[a_xz]*pC[c_xz])/(24*dz)+(pC[a_xyz]*pC[c_xz])/(48*dz)-(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)-(pC[a_xyz]*pC[c_x])/(24*dz)-(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz-(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)+
     (pC[a_yy]*pC[b_z])/(2*dy)-(pC[a_y]*pC[b_z])/(2*dy)-(pC[a_yz]*pC[b_yz])/(8*dy)-(pC[a_yy]*pC[b_yz])/(4*dy)+(pC[a_y]*pC[b_yz])/(4*dy)-(pC[a_yz]*pC[b_yy])/(12*dy)-(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)+(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)-(pC[a_y]*pC[b_y])/(2*dy)+(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)-(pC[a_xy]*pC[b_xy])/(24*dy)-(pC[a_xyz]*pC[b_x])/(24*dy)-(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)-(pC[a_yz]*pC[b_0])/(2*dy)-(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
     (pC[c_xzz]*pC[c_zz])/(36*dx)+(pC[c_xz]*pC[c_zz])/(12*dx)-(pC[c_xyz]*pC[c_zz])/(24*dx)+(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)+(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)+(pC[c_xyz]*pC[c_z])/(8*dx)-(pC[c_xy]*pC[c_z])/(4*dx)+(pC[c_x]*pC[c_z])/(2*dx)-(pC[c_xzz]*pC[c_yz])/(24*dx)+(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)+(pC[c_xy]*pC[c_yz])/(8*dx)-(pC[c_x]*pC[c_yz])/(4*dx)+(pC[c_xzz]*pC[c_y])/(12*dx)-(pC[c_xz]*pC[c_y])/(4*dx)+(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)+(pC[c_x]*pC[c_y])/(2*dx)-
     (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)+(pC[c_xx]*pC[c_xz])/(12*dx)+(pC[c_0]*pC[c_xz])/(2*dx)-(pC[c_0]*pC[c_xyz])/(4*dx)+(pC[c_0]*pC[c_xy])/(2*dx)+(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)+(pC[b_xyz]*pC[b_z])/(8*dx)+(pC[b_xyy]*pC[b_z])/(12*dx)-(pC[b_xy]*pC[b_z])/(4*dx)+(pC[b_x]*pC[b_z])/(2*dx)+(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)-(pC[b_xyy]*pC[b_yz])/(24*dx)+(pC[b_xy]*pC[b_yz])/(8*dx)-(pC[b_x]*pC[b_yz])/(4*dx)+(pC[b_xz]*pC[b_yy])/(12*dx)-
     (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)+(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)-(pC[b_xz]*pC[b_y])/(4*dx)+(pC[b_xyz]*pC[b_y])/(8*dx)+(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)+(pC[b_x]*pC[b_y])/(2*dx)+(pC[b_0]*pC[b_xz])/(2*dx)-(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)+(pC[b_xx]*pC[b_xy])/(12*dx)+(pC[b_0]*pC[b_xy])/(2*dx)+(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBY Background By
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermXComponents
 *
 */
template<typename REAL> inline
REAL JXBX_010_110(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBY,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz+(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)+(pC[c_xz]*BGBZ)/(2*dx)+(pC[c_xyz]*BGBZ)/(4*dx)-(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx-(pC[a_yz]*BGBY)/(2*dy)+(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy+(pC[b_xz]*BGBY)/(2*dx)+(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)-(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx-(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)+
   (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)-(pC[a_z]*pC[c_z])/(2*dz)-(pC[a_yz]*pC[c_z])/(4*dz)+(pC[a_zz]*pC[c_yz])/(4*dz)-(pC[a_z]*pC[c_yz])/(4*dz)-(pC[a_yz]*pC[c_yz])/(8*dz)-(pC[a_zz]*pC[c_y])/(2*dz)+(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)-(pC[a_xz]*pC[c_xz])/(24*dz)-(pC[a_xyz]*pC[c_xz])/(48*dz)-(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)+(pC[a_xyz]*pC[c_x])/(24*dz)-(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz+(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)-
   (pC[a_yy]*pC[b_z])/(2*dy)-(pC[a_y]*pC[b_z])/(2*dy)+(pC[a_yz]*pC[b_yz])/(8*dy)-(pC[a_yy]*pC[b_yz])/(4*dy)-(pC[a_y]*pC[b_yz])/(4*dy)-(pC[a_yz]*pC[b_yy])/(12*dy)+(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)-(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)+(pC[a_y]*pC[b_y])/(2*dy)-(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)+(pC[a_xy]*pC[b_xy])/(24*dy)-(pC[a_xyz]*pC[b_x])/(24*dy)+(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)-(pC[a_yz]*pC[b_0])/(2*dy)+(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
   (pC[c_xzz]*pC[c_zz])/(36*dx)+(pC[c_xz]*pC[c_zz])/(12*dx)+(pC[c_xyz]*pC[c_zz])/(24*dx)-(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)+(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)-(pC[c_xyz]*pC[c_z])/(8*dx)+(pC[c_xy]*pC[c_z])/(4*dx)+(pC[c_x]*pC[c_z])/(2*dx)+(pC[c_xzz]*pC[c_yz])/(24*dx)-(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)+(pC[c_xy]*pC[c_yz])/(8*dx)+(pC[c_x]*pC[c_yz])/(4*dx)-(pC[c_xzz]*pC[c_y])/(12*dx)+(pC[c_xz]*pC[c_y])/(4*dx)+(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)-(pC[c_x]*pC[c_y])/(2*dx)-
   (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)+(pC[c_xx]*pC[c_xz])/(12*dx)+(pC[c_0]*pC[c_xz])/(2*dx)+(pC[c_0]*pC[c_xyz])/(4*dx)-(pC[c_0]*pC[c_xy])/(2*dx)+(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)-(pC[b_xyz]*pC[b_z])/(8*dx)+(pC[b_xyy]*pC[b_z])/(12*dx)+(pC[b_xy]*pC[b_z])/(4*dx)+(pC[b_x]*pC[b_z])/(2*dx)-(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)+(pC[b_xyy]*pC[b_yz])/(24*dx)+(pC[b_xy]*pC[b_yz])/(8*dx)+(pC[b_x]*pC[b_yz])/(4*dx)+(pC[b_xz]*pC[b_yy])/(12*dx)+
   (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)-(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)+(pC[b_xz]*pC[b_y])/(4*dx)+(pC[b_xyz]*pC[b_y])/(8*dx)-(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)-(pC[b_x]*pC[b_y])/(2*dx)+(pC[b_0]*pC[b_xz])/(2*dx)+(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)-(pC[b_xx]*pC[b_xy])/(12*dx)-(pC[b_0]*pC[b_xy])/(2*dx)-(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBY Background By
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermXComponents
 *
 */
template<typename REAL> inline
REAL JXBX_001_101(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBY,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return(pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz-(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)-(pC[c_xz]*BGBZ)/(2*dx)+(pC[c_xyz]*BGBZ)/(4*dx)+(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx+(pC[a_yz]*BGBY)/(2*dy)-(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy-(pC[b_xz]*BGBY)/(2*dx)+(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)+(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx+(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)-
   (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)+(pC[a_z]*pC[c_z])/(2*dz)-(pC[a_yz]*pC[c_z])/(4*dz)-(pC[a_zz]*pC[c_yz])/(4*dz)-(pC[a_z]*pC[c_yz])/(4*dz)+(pC[a_yz]*pC[c_yz])/(8*dz)-(pC[a_zz]*pC[c_y])/(2*dz)-(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)+(pC[a_xz]*pC[c_xz])/(24*dz)-(pC[a_xyz]*pC[c_xz])/(48*dz)+(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)-(pC[a_xyz]*pC[c_x])/(24*dz)+(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz-(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)-
   (pC[a_yy]*pC[b_z])/(2*dy)+(pC[a_y]*pC[b_z])/(2*dy)-(pC[a_yz]*pC[b_yz])/(8*dy)+(pC[a_yy]*pC[b_yz])/(4*dy)-(pC[a_y]*pC[b_yz])/(4*dy)+(pC[a_yz]*pC[b_yy])/(12*dy)-(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)-(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)-(pC[a_y]*pC[b_y])/(2*dy)-(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)-(pC[a_xy]*pC[b_xy])/(24*dy)+(pC[a_xyz]*pC[b_x])/(24*dy)-(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)+(pC[a_yz]*pC[b_0])/(2*dy)-(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
   (pC[c_xzz]*pC[c_zz])/(36*dx)-(pC[c_xz]*pC[c_zz])/(12*dx)+(pC[c_xyz]*pC[c_zz])/(24*dx)+(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)-(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)+(pC[c_xyz]*pC[c_z])/(8*dx)+(pC[c_xy]*pC[c_z])/(4*dx)-(pC[c_x]*pC[c_z])/(2*dx)+(pC[c_xzz]*pC[c_yz])/(24*dx)+(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)-(pC[c_xy]*pC[c_yz])/(8*dx)+(pC[c_x]*pC[c_yz])/(4*dx)+(pC[c_xzz]*pC[c_y])/(12*dx)+(pC[c_xz]*pC[c_y])/(4*dx)-(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)+(pC[c_x]*pC[c_y])/(2*dx)-
   (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)-(pC[c_xx]*pC[c_xz])/(12*dx)-(pC[c_0]*pC[c_xz])/(2*dx)+(pC[c_0]*pC[c_xyz])/(4*dx)+(pC[c_0]*pC[c_xy])/(2*dx)-(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)+(pC[b_xyz]*pC[b_z])/(8*dx)-(pC[b_xyy]*pC[b_z])/(12*dx)+(pC[b_xy]*pC[b_z])/(4*dx)-(pC[b_x]*pC[b_z])/(2*dx)+(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)+(pC[b_xyy]*pC[b_yz])/(24*dx)-(pC[b_xy]*pC[b_yz])/(8*dx)+(pC[b_x]*pC[b_yz])/(4*dx)-(pC[b_xz]*pC[b_yy])/(12*dx)+
   (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)+(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)+(pC[b_xz]*pC[b_y])/(4*dx)-(pC[b_xyz]*pC[b_y])/(8*dx)+(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)+(pC[b_x]*pC[b_y])/(2*dx)-(pC[b_0]*pC[b_xz])/(2*dx)+(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)+(pC[b_xx]*pC[b_xy])/(12*dx)+(pC[b_0]*pC[b_xy])/(2*dx)+(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx ;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBY Background By
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermXComponents
 *
 */
template<typename REAL> inline
REAL JXBX_011_111(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBY,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return (pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz+(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)-(pC[c_xz]*BGBZ)/(2*dx)-(pC[c_xyz]*BGBZ)/(4*dx)-(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx+(pC[a_yz]*BGBY)/(2*dy)+(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy-(pC[b_xz]*BGBY)/(2*dx)-(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)-(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx+(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)+
   (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)+(pC[a_z]*pC[c_z])/(2*dz)+(pC[a_yz]*pC[c_z])/(4*dz)+(pC[a_zz]*pC[c_yz])/(4*dz)+(pC[a_z]*pC[c_yz])/(4*dz)+(pC[a_yz]*pC[c_yz])/(8*dz)+(pC[a_zz]*pC[c_y])/(2*dz)+(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)+(pC[a_xz]*pC[c_xz])/(24*dz)+(pC[a_xyz]*pC[c_xz])/(48*dz)+(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)+(pC[a_xyz]*pC[c_x])/(24*dz)+(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz+(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)+
   (pC[a_yy]*pC[b_z])/(2*dy)+(pC[a_y]*pC[b_z])/(2*dy)+(pC[a_yz]*pC[b_yz])/(8*dy)+(pC[a_yy]*pC[b_yz])/(4*dy)+(pC[a_y]*pC[b_yz])/(4*dy)+(pC[a_yz]*pC[b_yy])/(12*dy)+(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)+(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)+(pC[a_y]*pC[b_y])/(2*dy)+(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)+(pC[a_xy]*pC[b_xy])/(24*dy)+(pC[a_xyz]*pC[b_x])/(24*dy)+(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)+(pC[a_yz]*pC[b_0])/(2*dy)+(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
   (pC[c_xzz]*pC[c_zz])/(36*dx)-(pC[c_xz]*pC[c_zz])/(12*dx)-(pC[c_xyz]*pC[c_zz])/(24*dx)-(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)-(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)-(pC[c_xyz]*pC[c_z])/(8*dx)-(pC[c_xy]*pC[c_z])/(4*dx)-(pC[c_x]*pC[c_z])/(2*dx)-(pC[c_xzz]*pC[c_yz])/(24*dx)-(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)-(pC[c_xy]*pC[c_yz])/(8*dx)-(pC[c_x]*pC[c_yz])/(4*dx)-(pC[c_xzz]*pC[c_y])/(12*dx)-(pC[c_xz]*pC[c_y])/(4*dx)-(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)-(pC[c_x]*pC[c_y])/(2*dx)-
   (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)-(pC[c_xx]*pC[c_xz])/(12*dx)-(pC[c_0]*pC[c_xz])/(2*dx)-(pC[c_0]*pC[c_xyz])/(4*dx)-(pC[c_0]*pC[c_xy])/(2*dx)-(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)-(pC[b_xyz]*pC[b_z])/(8*dx)-(pC[b_xyy]*pC[b_z])/(12*dx)-(pC[b_xy]*pC[b_z])/(4*dx)-(pC[b_x]*pC[b_z])/(2*dx)-(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)-(pC[b_xyy]*pC[b_yz])/(24*dx)-(pC[b_xy]*pC[b_yz])/(8*dx)-(pC[b_x]*pC[b_yz])/(4*dx)-(pC[b_xz]*pC[b_yy])/(12*dx)-
   (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)-(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)-(pC[b_xz]*pC[b_y])/(4*dx)-(pC[b_xyz]*pC[b_y])/(8*dx)-(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)-(pC[b_x]*pC[b_y])/(2*dx)-(pC[b_0]*pC[b_xz])/(2*dx)-(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)-(pC[b_xx]*pC[b_xy])/(12*dx)-(pC[b_0]*pC[b_xy])/(2*dx)-(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx;
}

// Y
/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermYComponents
 *
 */
template<typename REAL> inline
REAL JXBY_000_010(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz-(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)+(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy-(pC[c_xyz]*BGBZ)/(4*dy)+(pC[c_xy]*BGBZ)/(2*dy)+(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy-(pC[a_xyz]*BGBX)/(4*dy)+(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)-(pC[b_xz]*BGBX)/(2*dx)-(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx-(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)-
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)-(pC[b_z]*pC[c_z])/(2*dz)+(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)-(pC[b_yz]*pC[c_yz])/(24*dz)+(pC[b_xyz]*pC[c_yz])/(48*dz)-(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)-(pC[b_xyz]*pC[c_y])/(24*dz)-(pC[b_zz]*pC[c_xz])/(4*dz)+(pC[b_z]*pC[c_xz])/(4*dz)-(pC[b_xz]*pC[c_xz])/(8*dz)+(pC[b_zz]*pC[c_x])/(2*dz)-(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)-(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz-(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)+
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)-(pC[c_xyz]*pC[c_zz])/(24*dy)+(pC[c_xy]*pC[c_zz])/(12*dy)+(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)+(pC[c_y]*pC[c_z])/(2*dy)+(pC[c_xyz]*pC[c_z])/(8*dy)-(pC[c_xy]*pC[c_z])/(4*dy)-(pC[c_xz]*pC[c_yzz])/(24*dy)+(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)+(pC[c_yy]*pC[c_yz])/(12*dy)+(pC[c_xz]*pC[c_yz])/(8*dy)-(pC[c_x]*pC[c_yz])/(4*dy)+(pC[c_0]*pC[c_yz])/(2*dy)+(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)-(pC[c_xz]*pC[c_y])/(4*dy)+
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)+(pC[c_xy]*pC[c_xz])/(8*dy)+(pC[c_x]*pC[c_xyz])/(8*dy)-(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)+(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)+(pC[a_y]*pC[a_z])/(2*dy)+(pC[a_xyz]*pC[a_z])/(8*dy)-(pC[a_xy]*pC[a_z])/(4*dy)+(pC[a_xxy]*pC[a_z])/(12*dy)+(pC[a_xz]*pC[a_yz])/(8*dy)+(pC[a_xx]*pC[a_yz])/(12*dy)-(pC[a_x]*pC[a_yz])/(4*dy)+(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)+(pC[a_xy]*pC[a_yy])/(12*dy)-(pC[a_xz]*pC[a_y])/(4*dy)+
     (pC[a_xyy]*pC[a_y])/(12*dy)-(pC[a_xx]*pC[a_y])/(6*dy)+(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)+(pC[a_xy]*pC[a_xz])/(8*dy)-(pC[a_xxy]*pC[a_xz])/(24*dy)-(pC[a_xx]*pC[a_xyz])/(24*dy)+(pC[a_x]*pC[a_xyz])/(8*dy)-(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)+(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)+(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)+(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)-(pC[a_xz]*pC[b_xz])/(8*dx)-
     (pC[a_xx]*pC[b_xz])/(12*dx)+(pC[a_x]*pC[b_xz])/(4*dx)-(pC[a_0]*pC[b_xz])/(2*dx)-(pC[a_y]*pC[b_xyz])/(24*dx)+(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)-(pC[a_xy]*pC[b_xy])/(24*dx)-(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)+(pC[a_z]*pC[b_xx])/(2*dx)-(pC[a_xz]*pC[b_xx])/(4*dx)-(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)-(pC[a_0]*pC[b_xx])/dx-(pC[a_z]*pC[b_x])/(2*dx)+(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)-(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermYComponents
 *
 */
template<typename REAL> inline
REAL JXBY_100_110(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz+(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)+(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy+(pC[c_xyz]*BGBZ)/(4*dy)-(pC[c_xy]*BGBZ)/(2*dy)+(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy+(pC[a_xyz]*BGBX)/(4*dy)-(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)-(pC[b_xz]*BGBX)/(2*dx)+(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx-(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)+
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)-(pC[b_z]*pC[c_z])/(2*dz)-(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)-(pC[b_yz]*pC[c_yz])/(24*dz)-(pC[b_xyz]*pC[c_yz])/(48*dz)-(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)+(pC[b_xyz]*pC[c_y])/(24*dz)+(pC[b_zz]*pC[c_xz])/(4*dz)-(pC[b_z]*pC[c_xz])/(4*dz)-(pC[b_xz]*pC[c_xz])/(8*dz)-(pC[b_zz]*pC[c_x])/(2*dz)+(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)-(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz+(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)+
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)+(pC[c_xyz]*pC[c_zz])/(24*dy)-(pC[c_xy]*pC[c_zz])/(12*dy)+(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)+(pC[c_y]*pC[c_z])/(2*dy)-(pC[c_xyz]*pC[c_z])/(8*dy)+(pC[c_xy]*pC[c_z])/(4*dy)+(pC[c_xz]*pC[c_yzz])/(24*dy)-(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)+(pC[c_yy]*pC[c_yz])/(12*dy)-(pC[c_xz]*pC[c_yz])/(8*dy)+(pC[c_x]*pC[c_yz])/(4*dy)+(pC[c_0]*pC[c_yz])/(2*dy)+(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)+(pC[c_xz]*pC[c_y])/(4*dy)-
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)+(pC[c_xy]*pC[c_xz])/(8*dy)+(pC[c_x]*pC[c_xyz])/(8*dy)+(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)-(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)+(pC[a_y]*pC[a_z])/(2*dy)-(pC[a_xyz]*pC[a_z])/(8*dy)+(pC[a_xy]*pC[a_z])/(4*dy)+(pC[a_xxy]*pC[a_z])/(12*dy)-(pC[a_xz]*pC[a_yz])/(8*dy)+(pC[a_xx]*pC[a_yz])/(12*dy)+(pC[a_x]*pC[a_yz])/(4*dy)+(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)-(pC[a_xy]*pC[a_yy])/(12*dy)+(pC[a_xz]*pC[a_y])/(4*dy)-
     (pC[a_xyy]*pC[a_y])/(12*dy)-(pC[a_xx]*pC[a_y])/(6*dy)-(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)+(pC[a_xy]*pC[a_xz])/(8*dy)+(pC[a_xxy]*pC[a_xz])/(24*dy)+(pC[a_xx]*pC[a_xyz])/(24*dy)+(pC[a_x]*pC[a_xyz])/(8*dy)+(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)-(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)-(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)-(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)+(pC[a_xz]*pC[b_xz])/(8*dx)-
     (pC[a_xx]*pC[b_xz])/(12*dx)-(pC[a_x]*pC[b_xz])/(4*dx)-(pC[a_0]*pC[b_xz])/(2*dx)-(pC[a_y]*pC[b_xyz])/(24*dx)-(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)+(pC[a_xy]*pC[b_xy])/(24*dx)+(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)-(pC[a_z]*pC[b_xx])/(2*dx)-(pC[a_xz]*pC[b_xx])/(4*dx)+(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)+(pC[a_0]*pC[b_xx])/dx-(pC[a_z]*pC[b_x])/(2*dx)-(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)+(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermYComponents
 *
 */
template<typename REAL> inline
REAL JXBY_001_011(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return (pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz-(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)-(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy+(pC[c_xyz]*BGBZ)/(4*dy)+(pC[c_xy]*BGBZ)/(2*dy)-(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy+(pC[a_xyz]*BGBX)/(4*dy)+(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)+(pC[b_xz]*BGBX)/(2*dx)-(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx+(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)-
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)+(pC[b_z]*pC[c_z])/(2*dz)-(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)+(pC[b_yz]*pC[c_yz])/(24*dz)-(pC[b_xyz]*pC[c_yz])/(48*dz)+(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)-(pC[b_xyz]*pC[c_y])/(24*dz)-(pC[b_zz]*pC[c_xz])/(4*dz)-(pC[b_z]*pC[c_xz])/(4*dz)+(pC[b_xz]*pC[c_xz])/(8*dz)-(pC[b_zz]*pC[c_x])/(2*dz)-(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)+(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz-(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)-
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)+(pC[c_xyz]*pC[c_zz])/(24*dy)+(pC[c_xy]*pC[c_zz])/(12*dy)-(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)-(pC[c_y]*pC[c_z])/(2*dy)+(pC[c_xyz]*pC[c_z])/(8*dy)+(pC[c_xy]*pC[c_z])/(4*dy)+(pC[c_xz]*pC[c_yzz])/(24*dy)+(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)-(pC[c_yy]*pC[c_yz])/(12*dy)+(pC[c_xz]*pC[c_yz])/(8*dy)+(pC[c_x]*pC[c_yz])/(4*dy)-(pC[c_0]*pC[c_yz])/(2*dy)-(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)+(pC[c_xz]*pC[c_y])/(4*dy)+
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)-(pC[c_xy]*pC[c_xz])/(8*dy)-(pC[c_x]*pC[c_xyz])/(8*dy)+(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)+(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)-(pC[a_y]*pC[a_z])/(2*dy)+(pC[a_xyz]*pC[a_z])/(8*dy)+(pC[a_xy]*pC[a_z])/(4*dy)-(pC[a_xxy]*pC[a_z])/(12*dy)+(pC[a_xz]*pC[a_yz])/(8*dy)-(pC[a_xx]*pC[a_yz])/(12*dy)+(pC[a_x]*pC[a_yz])/(4*dy)-(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)+(pC[a_xy]*pC[a_yy])/(12*dy)+(pC[a_xz]*pC[a_y])/(4*dy)+
     (pC[a_xyy]*pC[a_y])/(12*dy)-(pC[a_xx]*pC[a_y])/(6*dy)+(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)-(pC[a_xy]*pC[a_xz])/(8*dy)+(pC[a_xxy]*pC[a_xz])/(24*dy)+(pC[a_xx]*pC[a_xyz])/(24*dy)-(pC[a_x]*pC[a_xyz])/(8*dy)+(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)+(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)+(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)+(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)-(pC[a_xz]*pC[b_xz])/(8*dx)+
     (pC[a_xx]*pC[b_xz])/(12*dx)-(pC[a_x]*pC[b_xz])/(4*dx)+(pC[a_0]*pC[b_xz])/(2*dx)+(pC[a_y]*pC[b_xyz])/(24*dx)-(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)-(pC[a_xy]*pC[b_xy])/(24*dx)-(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)-(pC[a_z]*pC[b_xx])/(2*dx)+(pC[a_xz]*pC[b_xx])/(4*dx)-(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)-(pC[a_0]*pC[b_xx])/dx+(pC[a_z]*pC[b_x])/(2*dx)-(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)-(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBZ Background Bz
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermYComponents
 *
 */
template<typename REAL> inline
REAL JXBY_101_111(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBZ,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return (pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz+(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)-(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy-(pC[c_xyz]*BGBZ)/(4*dy)-(pC[c_xy]*BGBZ)/(2*dy)-(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy-(pC[a_xyz]*BGBX)/(4*dy)-(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)+(pC[b_xz]*BGBX)/(2*dx)+(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx+(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)+
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)+(pC[b_z]*pC[c_z])/(2*dz)+(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)+(pC[b_yz]*pC[c_yz])/(24*dz)+(pC[b_xyz]*pC[c_yz])/(48*dz)+(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)+(pC[b_xyz]*pC[c_y])/(24*dz)+(pC[b_zz]*pC[c_xz])/(4*dz)+(pC[b_z]*pC[c_xz])/(4*dz)+(pC[b_xz]*pC[c_xz])/(8*dz)+(pC[b_zz]*pC[c_x])/(2*dz)+(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)+(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz+(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)-
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)-(pC[c_xyz]*pC[c_zz])/(24*dy)-(pC[c_xy]*pC[c_zz])/(12*dy)-(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)-(pC[c_y]*pC[c_z])/(2*dy)-(pC[c_xyz]*pC[c_z])/(8*dy)-(pC[c_xy]*pC[c_z])/(4*dy)-(pC[c_xz]*pC[c_yzz])/(24*dy)-(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)-(pC[c_yy]*pC[c_yz])/(12*dy)-(pC[c_xz]*pC[c_yz])/(8*dy)-(pC[c_x]*pC[c_yz])/(4*dy)-(pC[c_0]*pC[c_yz])/(2*dy)-(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)-(pC[c_xz]*pC[c_y])/(4*dy)-
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)-(pC[c_xy]*pC[c_xz])/(8*dy)-(pC[c_x]*pC[c_xyz])/(8*dy)-(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)-(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)-(pC[a_y]*pC[a_z])/(2*dy)-(pC[a_xyz]*pC[a_z])/(8*dy)-(pC[a_xy]*pC[a_z])/(4*dy)-(pC[a_xxy]*pC[a_z])/(12*dy)-(pC[a_xz]*pC[a_yz])/(8*dy)-(pC[a_xx]*pC[a_yz])/(12*dy)-(pC[a_x]*pC[a_yz])/(4*dy)-(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)-(pC[a_xy]*pC[a_yy])/(12*dy)-(pC[a_xz]*pC[a_y])/(4*dy)-(pC[a_xyy]*pC[a_y])/(12*dy)-
     (pC[a_xx]*pC[a_y])/(6*dy)-(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)-(pC[a_xy]*pC[a_xz])/(8*dy)-(pC[a_xxy]*pC[a_xz])/(24*dy)-(pC[a_xx]*pC[a_xyz])/(24*dy)-(pC[a_x]*pC[a_xyz])/(8*dy)-(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)-(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)-(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)-(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)+(pC[a_xz]*pC[b_xz])/(8*dx)+(pC[a_xx]*pC[b_xz])/(12*dx)+(pC[a_x]*pC[b_xz])/(4*dx)+
     (pC[a_0]*pC[b_xz])/(2*dx)+(pC[a_y]*pC[b_xyz])/(24*dx)+(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)+(pC[a_xy]*pC[b_xy])/(24*dx)+(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)+(pC[a_z]*pC[b_xx])/(2*dx)+(pC[a_xz]*pC[b_xx])/(4*dx)+(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)+(pC[a_0]*pC[b_xx])/dx+(pC[a_z]*pC[b_x])/(2*dx)+(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)+(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

// Z
/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBY Background By
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermZComponents
 *
 */
template<typename REAL> inline
REAL JXBZ_000_001(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBY,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz+(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)+(pC[b_xz]*BGBY)/(2*dz)-(pC[b_xyz]*BGBY)/(4*dz)-(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy-(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz+(pC[a_yz]*BGBX)/(2*dz)+(pC[a_xz]*BGBX)/(2*dz)-(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)-(pC[c_xy]*BGBX)/(2*dx)-(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)+(pC[b_yz]*pC[b_zz])/(12*dz)+
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)+(pC[b_y]*pC[b_z])/(2*dz)-(pC[b_xy]*pC[b_z])/(4*dz)+(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)+(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)+(pC[b_xy]*pC[b_yz])/(8*dz)-(pC[b_x]*pC[b_yz])/(4*dz)+(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)+(pC[b_y]*pC[b_yyz])/(12*dz)-(pC[b_xy]*pC[b_yyz])/(24*dz)+(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)+(pC[b_xz]*pC[b_yy])/(12*dz)-(pC[b_xyz]*pC[b_yy])/(24*dz)-(pC[b_xz]*pC[b_y])/(4*dz)+
     (pC[b_xyz]*pC[b_y])/(8*dz)+(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)+(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)+(pC[b_x]*pC[b_xyz])/(8*dz)-(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)+(pC[a_xz]*pC[a_zz])/(12*dz)+(pC[a_y]*pC[a_z])/(2*dz)+(pC[a_xzz]*pC[a_z])/(12*dz)-(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)+(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)+(pC[a_xy]*pC[a_yz])/(8*dz)+(pC[a_xx]*pC[a_yz])/(12*dz)-(pC[a_x]*pC[a_yz])/(4*dz)+(pC[a_0]*pC[a_yz])/(2*dz)-
     (pC[a_xz]*pC[a_y])/(4*dz)+(pC[a_xyz]*pC[a_y])/(8*dz)+(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)+(pC[a_xy]*pC[a_xz])/(8*dz)+(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)+(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)-(pC[a_xx]*pC[a_xyz])/(24*dz)+(pC[a_x]*pC[a_xyz])/(8*dz)-(pC[a_0]*pC[a_xyz])/(4*dz)-(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)+(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)-(pC[b_yz]*pC[c_yz])/(24*dy)-(pC[b_z]*pC[c_yyz])/(12*dy)+
     (pC[b_yz]*pC[c_yyz])/(24*dy)-(pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)-(pC[b_xy]*pC[c_yy])/(4*dy)+(pC[b_x]*pC[c_yy])/(2*dy)-(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)-(pC[b_y]*pC[c_y])/(2*dy)+(pC[b_xy]*pC[c_y])/(4*dy)-(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy-(pC[b_z]*pC[c_xyz])/(24*dy)+(pC[b_yz]*pC[c_xyz])/(48*dy)-(pC[b_yy]*pC[c_xy])/(12*dy)+(pC[b_y]*pC[c_xy])/(4*dy)-(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)-(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)-(pC[a_xz]*pC[c_xz])/(24*dx)-
     (pC[a_z]*pC[c_xyz])/(24*dx)+(pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)-(pC[a_xy]*pC[c_xy])/(8*dx)-(pC[a_xx]*pC[c_xy])/(12*dx)+(pC[a_x]*pC[c_xy])/(4*dx)-(pC[a_0]*pC[c_xy])/(2*dx)-(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)+(pC[a_y]*pC[c_xx])/(2*dx)-(pC[a_xy]*pC[c_xx])/(4*dx)-(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)-(pC[a_0]*pC[c_xx])/dx-(pC[a_y]*pC[c_x])/(2*dx)+(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)-(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBY Background By
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermZComponents
 *
 */
template<typename REAL> inline
REAL JXBZ_100_101(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBY,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz+(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)-(pC[b_xz]*BGBY)/(2*dz)+(pC[b_xyz]*BGBY)/(4*dz)-(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy+(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz+(pC[a_yz]*BGBX)/(2*dz)-(pC[a_xz]*BGBX)/(2*dz)+(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)-(pC[c_xy]*BGBX)/(2*dx)+(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)+(pC[b_yz]*pC[b_zz])/(12*dz)+
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)+(pC[b_y]*pC[b_z])/(2*dz)+(pC[b_xy]*pC[b_z])/(4*dz)-(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)+(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)-(pC[b_xy]*pC[b_yz])/(8*dz)+(pC[b_x]*pC[b_yz])/(4*dz)+(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)+(pC[b_y]*pC[b_yyz])/(12*dz)+(pC[b_xy]*pC[b_yyz])/(24*dz)-(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)-(pC[b_xz]*pC[b_yy])/(12*dz)+(pC[b_xyz]*pC[b_yy])/(24*dz)+(pC[b_xz]*pC[b_y])/(4*dz)-
     (pC[b_xyz]*pC[b_y])/(8*dz)+(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)-(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)+(pC[b_x]*pC[b_xyz])/(8*dz)+(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)-(pC[a_xz]*pC[a_zz])/(12*dz)+(pC[a_y]*pC[a_z])/(2*dz)-(pC[a_xzz]*pC[a_z])/(12*dz)+(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)-(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)-(pC[a_xy]*pC[a_yz])/(8*dz)+(pC[a_xx]*pC[a_yz])/(12*dz)+(pC[a_x]*pC[a_yz])/(4*dz)+(pC[a_0]*pC[a_yz])/(2*dz)+(pC[a_xz]*pC[a_y])/(4*dz)
       -(pC[a_xyz]*pC[a_y])/(8*dz)+(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)+(pC[a_xy]*pC[a_xz])/(8*dz)-(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)-(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)+(pC[a_xx]*pC[a_xyz])/(24*dz)+(pC[a_x]*pC[a_xyz])/(8*dz)+(pC[a_0]*pC[a_xyz])/(4*dz)+(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)-(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)-(pC[b_yz]*pC[c_yz])/(24*dy)-(pC[b_z]*pC[c_yyz])/(12*dy)+(pC[b_yz]*pC[c_yyz])/(24*dy)-
     (pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)+(pC[b_xy]*pC[c_yy])/(4*dy)-(pC[b_x]*pC[c_yy])/(2*dy)-(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)-(pC[b_y]*pC[c_y])/(2*dy)-(pC[b_xy]*pC[c_y])/(4*dy)+(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy+(pC[b_z]*pC[c_xyz])/(24*dy)-(pC[b_yz]*pC[c_xyz])/(48*dy)+(pC[b_yy]*pC[c_xy])/(12*dy)-(pC[b_y]*pC[c_xy])/(4*dy)-(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)+(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)+(pC[a_xz]*pC[c_xz])/(24*dx)-(pC[a_z]*pC[c_xyz])/(24*dx)-
     (pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)+(pC[a_xy]*pC[c_xy])/(8*dx)-(pC[a_xx]*pC[c_xy])/(12*dx)-(pC[a_x]*pC[c_xy])/(4*dx)-(pC[a_0]*pC[c_xy])/(2*dx)+(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)-(pC[a_y]*pC[c_xx])/(2*dx)-(pC[a_xy]*pC[c_xx])/(4*dx)+(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)+(pC[a_0]*pC[c_xx])/dx-(pC[a_y]*pC[c_x])/(2*dx)-(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)+(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBY Background By
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermZComponents
 *
 */
template<typename REAL> inline
REAL JXBZ_010_011(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBY,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz-(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)+(pC[b_xz]*BGBY)/(2*dz)+(pC[b_xyz]*BGBY)/(4*dz)+(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy-(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz-(pC[a_yz]*BGBX)/(2*dz)+(pC[a_xz]*BGBX)/(2*dz)+(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)+(pC[c_xy]*BGBX)/(2*dx)-(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)-(pC[b_yz]*pC[b_zz])/(12*dz)-
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)-(pC[b_y]*pC[b_z])/(2*dz)+(pC[b_xy]*pC[b_z])/(4*dz)+(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)-(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)+(pC[b_xy]*pC[b_yz])/(8*dz)+(pC[b_x]*pC[b_yz])/(4*dz)-(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)-(pC[b_y]*pC[b_yyz])/(12*dz)+(pC[b_xy]*pC[b_yyz])/(24*dz)+(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)+(pC[b_xz]*pC[b_yy])/(12*dz)+(pC[b_xyz]*pC[b_yy])/(24*dz)+(pC[b_xz]*pC[b_y])/(4*dz)+
     (pC[b_xyz]*pC[b_y])/(8*dz)-(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)+(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)-(pC[b_x]*pC[b_xyz])/(8*dz)+(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)+(pC[a_xz]*pC[a_zz])/(12*dz)-(pC[a_y]*pC[a_z])/(2*dz)+(pC[a_xzz]*pC[a_z])/(12*dz)+(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)+(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)+(pC[a_xy]*pC[a_yz])/(8*dz)-(pC[a_xx]*pC[a_yz])/(12*dz)+(pC[a_x]*pC[a_yz])/(4*dz)-(pC[a_0]*pC[a_yz])/(2*dz)+(pC[a_xz]*pC[a_y])/(4*dz)
       +(pC[a_xyz]*pC[a_y])/(8*dz)-(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)-(pC[a_xy]*pC[a_xz])/(8*dz)+(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)+(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)+(pC[a_xx]*pC[a_xyz])/(24*dz)-(pC[a_x]*pC[a_xyz])/(8*dz)+(pC[a_0]*pC[a_xyz])/(4*dz)+(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)+(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)+(pC[b_yz]*pC[c_yz])/(24*dy)+(pC[b_z]*pC[c_yyz])/(12*dy)+(pC[b_yz]*pC[c_yyz])/(24*dy)
         +(pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)-(pC[b_xy]*pC[c_yy])/(4*dy)-(pC[b_x]*pC[c_yy])/(2*dy)+(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)+(pC[b_y]*pC[c_y])/(2*dy)-(pC[b_xy]*pC[c_y])/(4*dy)-(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy-(pC[b_z]*pC[c_xyz])/(24*dy)-(pC[b_yz]*pC[c_xyz])/(48*dy)-(pC[b_yy]*pC[c_xy])/(12*dy)-(pC[b_y]*pC[c_xy])/(4*dy)+(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)-(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)-(pC[a_xz]*pC[c_xz])/(24*dx)+(pC[a_z]*pC[c_xyz])/(24*dx)-
     (pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)-(pC[a_xy]*pC[c_xy])/(8*dx)+(pC[a_xx]*pC[c_xy])/(12*dx)-(pC[a_x]*pC[c_xy])/(4*dx)+(pC[a_0]*pC[c_xy])/(2*dx)-(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)-(pC[a_y]*pC[c_xx])/(2*dx)+(pC[a_xy]*pC[c_xx])/(4*dx)-(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)-(pC[a_0]*pC[c_xx])/dx+(pC[a_y]*pC[c_x])/(2*dx)-(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)-(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

/*! \brief Low-level Hall component computation
 *
 * Hall term computation following Balsara reconstruction, edge-averaged.
 *
 * \param pC Reconstruction coefficients
 * \param BGBX Background Bx
 * \param BGBY Background By
 * \param dx Cell dx
 * \param dy Cell dy
 * \param dz Cell dz
 *
 * \sa calculateEdgeHallTermZComponents
 *
 */
template<typename REAL> inline
REAL JXBZ_110_111(
   const std::array<REAL, Rec::N_REC_COEFFICIENTS> & pC,
   creal BGBX,
   creal BGBY,
   creal dx,
   creal dy,
   creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz-(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)-(pC[b_xz]*BGBY)/(2*dz)-(pC[b_xyz]*BGBY)/(4*dz)+(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy+(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz-(pC[a_yz]*BGBX)/(2*dz)-(pC[a_xz]*BGBX)/(2*dz)-(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)+(pC[c_xy]*BGBX)/(2*dx)+(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)-(pC[b_yz]*pC[b_zz])/(12*dz)-
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)-(pC[b_y]*pC[b_z])/(2*dz)-(pC[b_xy]*pC[b_z])/(4*dz)-(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)-(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)-(pC[b_xy]*pC[b_yz])/(8*dz)-(pC[b_x]*pC[b_yz])/(4*dz)-(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)-(pC[b_y]*pC[b_yyz])/(12*dz)-(pC[b_xy]*pC[b_yyz])/(24*dz)-(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)-(pC[b_xz]*pC[b_yy])/(12*dz)-(pC[b_xyz]*pC[b_yy])/(24*dz)-(pC[b_xz]*pC[b_y])/(4*dz)-
     (pC[b_xyz]*pC[b_y])/(8*dz)-(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)-(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)-(pC[b_x]*pC[b_xyz])/(8*dz)-(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)-(pC[a_xz]*pC[a_zz])/(12*dz)-(pC[a_y]*pC[a_z])/(2*dz)-(pC[a_xzz]*pC[a_z])/(12*dz)-(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)-(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)-(pC[a_xy]*pC[a_yz])/(8*dz)-(pC[a_xx]*pC[a_yz])/(12*dz)-(pC[a_x]*pC[a_yz])/(4*dz)-(pC[a_0]*pC[a_yz])/(2*dz)-(pC[a_xz]*pC[a_y])/(4*dz)-
     (pC[a_xyz]*pC[a_y])/(8*dz)-(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)-(pC[a_xy]*pC[a_xz])/(8*dz)-(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)-(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)-(pC[a_xx]*pC[a_xyz])/(24*dz)-(pC[a_x]*pC[a_xyz])/(8*dz)-(pC[a_0]*pC[a_xyz])/(4*dz)-(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)-(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)+(pC[b_yz]*pC[c_yz])/(24*dy)+(pC[b_z]*pC[c_yyz])/(12*dy)+(pC[b_yz]*pC[c_yyz])/(24*dy)+
     (pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)+(pC[b_xy]*pC[c_yy])/(4*dy)+(pC[b_x]*pC[c_yy])/(2*dy)+(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)+(pC[b_y]*pC[c_y])/(2*dy)+(pC[b_xy]*pC[c_y])/(4*dy)+(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy+(pC[b_z]*pC[c_xyz])/(24*dy)+(pC[b_yz]*pC[c_xyz])/(48*dy)+(pC[b_yy]*pC[c_xy])/(12*dy)+(pC[b_y]*pC[c_xy])/(4*dy)+(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)+(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)+(pC[a_xz]*pC[c_xz])/(24*dx)+(pC[a_z]*pC[c_xyz])/(24*dx)+
     (pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)+(pC[a_xy]*pC[c_xy])/(8*dx)+(pC[a_xx]*pC[c_xy])/(12*dx)+(pC[a_x]*pC[c_xy])/(4*dx)+(pC[a_0]*pC[c_xy])/(2*dx)+(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)+(pC[a_y]*pC[c_xx])/(2*dx)+(pC[a_xy]*pC[c_xx])/(4*dx)+(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)+(pC[a_0]*pC[c_xx])/dx+(pC[a_y]*pC[c_x])/(2*dx)+(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)+(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

/*! \brief Low-level function computing the Hall term numerator x components.
 *
 * Calls the lower-level inline templates and scales the components properly.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param perturbedCoefficients Reconstruction coefficients
 * \param i,j,k fsGrid cell coordinates for the current cell
 *
 * \sa calculateHallTerm JXBX_000_100 JXBX_001_101 JXBX_010_110 JXBX_011_111
 *
 */
void calculateEdgeHallTermXComponents(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   const std::array<Real, Rec::N_REC_COEFFICIENTS> & perturbedCoefficients,
   cint i,
   cint j,
   cint k
) {
   Real By = 0.0;
   Real Bz = 0.0;
   Real hallRhoq = 0.0;
   Real EXHall = 0.0;

   switch (Parameters::ohmHallTerm) {
    case 0:
      cerr << __FILE__ << __LINE__ << "You shouldn't be in a Hall term function if Parameters::ohmHallTerm == 0." << endl;
      break;

    case 1:
      By = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBY)+BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY);
      Bz = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBZ)+BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ);

      hallRhoq =  (momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ) <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ) ;
      EXHall = Bz*((BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBxdz)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBxdz)) / technicalGrid.DZ -
                  (BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBzdx)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBzdx)) / technicalGrid.DX) -
               By*((BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBydx)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBydx)) / technicalGrid.DX -
                  (BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBxdy)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBxdy)) / technicalGrid.DY);
      EXHall /= physicalconstants::MU_0 * hallRhoq;

      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EXHALL_000_100) =
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EXHALL_010_110) =
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EXHALL_001_101) =
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EXHALL_011_111) = EXHall;

      break;
    case 2:
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j-1,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j  ,k-1)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j-1,k-1)->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EXHALL_000_100) = JXBX_000_100(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j+1,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j  ,k-1)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j+1,k-1)->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EXHALL_010_110) = JXBX_010_110(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j-1,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j  ,k+1)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j-1,k+1)->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EXHALL_001_101) = JXBX_001_101(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j+1,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j  ,k+1)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j+1,k+1)->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EXHALL_011_111) = JXBX_011_111(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      break;

    default:
      cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
      break;
   }
}

/*! \brief Low-level function computing the Hall term numerator y components.
 *
 * Calls the lower-level inline templates and scales the components properly.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param perturbedCoefficients Reconstruction coefficients
 * \param i,j,k fsGrid cell coordinates for the current cell
 *
 * \sa calculateHallTerm JXBY_000_010 JXBY_001_011 JXBY_100_110 JXBY_101_111
 *
 */
void calculateEdgeHallTermYComponents(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   const std::array<Real, Rec::N_REC_COEFFICIENTS> & perturbedCoefficients,
   cint i,
   cint j,
   cint k
) {
   Real Bx = 0.0;
   Real Bz = 0.0;
   Real hallRhoq = 0.0;
   Real EYHall = 0.0;

   switch (Parameters::ohmHallTerm) {
    case 0:
      cerr << __FILE__ << __LINE__ << "You shouldn't be in a Hall term function if Parameters::ohmHallTerm == 0." << endl;
      break;

    case 1:
      Bx = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBX)+BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX);
      Bz = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBZ)+BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ);

      hallRhoq =  (momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ) <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ) ;
      EYHall = Bx*((BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBydx)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBydx)) / technicalGrid.DX -
                  (BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBxdy)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBxdy)) / technicalGrid.DY) -
               Bz*((BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBzdy)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBzdy)) / technicalGrid.DY -
                  (BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBydz)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBydz)) / technicalGrid.DZ);
      EYHall /= physicalconstants::MU_0 * hallRhoq;

      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EYHALL_000_010) =
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EYHALL_100_110) =
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EYHALL_101_111) =
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EYHALL_001_011) = EYHall;
      break;

    case 2:
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i-1,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j  ,k-1)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i-1,j  ,k-1)->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EYHALL_000_010) = JXBY_000_010(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i+1,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j  ,k-1)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i+1,j  ,k-1)->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EYHALL_100_110) = JXBY_100_110(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i-1,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j  ,k+1)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i-1,j  ,k+1)->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EYHALL_001_011) = JXBY_001_011(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i+1,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j  ,k+1)->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i+1,j  ,k+1)->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EYHALL_101_111) = JXBY_101_111(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBZ), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      break;

    default:
      cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
      break;
   }
}

/*! \brief Low-level function computing the Hall term numerator z components.
 *
 * Calls the lower-level inline templates and scales the components properly.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param perturbedCoefficients Reconstruction coefficients
 * \param i,j,k fsGrid cell coordinates for the current cell
 *
 * \sa calculateHallTerm JXBZ_000_001 JXBZ_010_011 JXBZ_100_101 JXBZ_110_111
 *
 */
void calculateEdgeHallTermZComponents(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   const std::array<Real, Rec::N_REC_COEFFICIENTS> & perturbedCoefficients,
   cint i,
   cint j,
   cint k
) {
   Real Bx = 0.0;
   Real By = 0.0;
   Real hallRhoq = 0.0;
   Real EZHall = 0.0;

   switch (Parameters::ohmHallTerm) {
   case 0:
     cerr << __FILE__ << __LINE__ << "You shouldn't be in a Hall term function if Parameters::ohmHallTerm == 0." << endl;
     break;

   case 1:
     Bx = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBX)+BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX);
     By = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBY)+BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY);

     hallRhoq =  (momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ) <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : momentsGrid.get(i,j,k)->at(fsgrids::moments::RHOQ) ;
     EZHall = By*((BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBzdy)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBzdy)) / technicalGrid.DY -
                 (BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBydz)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBydz)) / technicalGrid.DZ) -
              Bx*((BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBxdz)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBxdz)) / technicalGrid.DZ -
                 (BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::dBGBzdx)+dPerBGrid.get(i,j,k)->at(fsgrids::dperb::dPERBzdx)) / technicalGrid.DX);
     EZHall /= physicalconstants::MU_0 * hallRhoq;

     EHallGrid.get(i,j,k)->at(fsgrids::ehall::EZHALL_000_001) =
     EHallGrid.get(i,j,k)->at(fsgrids::ehall::EZHALL_100_101) =
     EHallGrid.get(i,j,k)->at(fsgrids::ehall::EZHALL_110_111) =
     EHallGrid.get(i,j,k)->at(fsgrids::ehall::EZHALL_010_011) = EZHall;
     break;

   case 2:
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i-1,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j-1,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i-1,j-1,k  )->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EZHALL_000_001) = JXBZ_000_001(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i+1,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j-1,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i+1,j-1,k  )->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EZHALL_100_101) = JXBZ_100_101(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i-1,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j+1,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i-1,j+1,k  )->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EZHALL_010_011) = JXBZ_010_011(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      hallRhoq = FOURTH * (
         momentsGrid.get(i  ,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i+1,j  ,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i  ,j+1,k  )->at(fsgrids::moments::RHOQ) +
         momentsGrid.get(i+1,j+1,k  )->at(fsgrids::moments::RHOQ)
      );
      hallRhoq =  (hallRhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : hallRhoq ;
      EHallGrid.get(i,j,k)->at(fsgrids::ehall::EZHALL_110_111) = JXBZ_110_111(perturbedCoefficients, BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBX), BgBGrid.get(i,j,k)->at(fsgrids::bgbfield::BGBY), technicalGrid.DX, technicalGrid.DY, technicalGrid.DZ) / (physicalconstants::MU_0 * hallRhoq);
      break;

    default:
      cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
      break;
   }
}

/** \brief Calculate the numerator of the Hall term on all given cells.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary condition functions.
 * \param i,j,k fsGrid cell coordinates for the current cell
 *
 * \sa calculateHallTermSimple calculateEdgeHallTermXComponents calculateEdgeHallTermYComponents calculateEdgeHallTermZComponents
 */
void calculateHallTerm(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint i,
   cint j,
   cint k
) {

   #ifdef DEBUG_FSOLVER
   if (technicalGrid.get(i,j,k) == NULL) {
      cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
   #endif

   cuint cellSysBoundaryFlag = technicalGrid.get(i,j,k)->sysBoundaryFlag;

   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) return;

   cuint cellSysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;

   std::array<Real, Rec::N_REC_COEFFICIENTS> perturbedCoefficients;

   reconstructionCoefficients(
      perBGrid,
      dPerBGrid,
      perturbedCoefficients,
      i,
      j,
      k,
      3 // Reconstruction order of the fields after Balsara 2009, 2 used for general B, 3 used here for 2nd-order Hall term
   );

   if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && (cellSysBoundaryLayer != 1)) {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, 0);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, 1);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, 2);
   } else {
      calculateEdgeHallTermXComponents(perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid, technicalGrid, perturbedCoefficients, i, j, k);
      calculateEdgeHallTermYComponents(perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid, technicalGrid, perturbedCoefficients, i, j, k);
      calculateEdgeHallTermZComponents(perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid, technicalGrid, perturbedCoefficients, i, j, k);
   }

}

/*! \brief High-level function computing the Hall term.
 *
 * Performs the communication before and after the computation as well as the computation of all Hall term numerator components.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param perBDt2Grid fsGrid holding the perturbed B quantities at runge-kutta half step
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param momentsGrid fsGrid holding the moment quantities
 * \param momentsDt2Grid fsGrid holding the moment quantities at runge-kutta half step
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary condition functions.
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param communicateMomentsDerivatives whether to communicate derivatves with the neighbour CPUs
 *
 * \sa calculateHallTerm
 */
void calculateHallTermSimple(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase
) {
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const int* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];

   phiprof::Timer hallTimer {"Calculate Hall term"};
   phiprof::Timer mpiTimer {"EHall ghost updates MPI", {"MPI"}};
   int computeTimerId {phiprof::initializeTimer("EHall compute cells")};
   dPerBGrid.updateGhostCells();
   if(P::ohmGradPeTerm == 0) {
      dMomentsGrid.updateGhostCells();
   }
   mpiTimer.stop();

   #pragma omp parallel
   {
      phiprof::Timer computeTimer {computeTimerId};
      #pragma omp for collapse(2)
      for (int k=0; k<gridDims[2]; k++) {
	 for (int j=0; j<gridDims[1]; j++) {
	    for (int i=0; i<gridDims[0]; i++) {
	       if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
		  calculateHallTerm(perBGrid, EHallGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid, technicalGrid,sysBoundaries, i, j, k);
	       } else {
		  calculateHallTerm(perBDt2Grid, EHallGrid, momentsDt2Grid, dPerBGrid, dMomentsGrid, BgBGrid, technicalGrid,sysBoundaries, i, j, k);
	       }
	    }
	 }
      }
      computeTimer.stop(N_cells,"Spatial Cells");
   }

   hallTimer.stop(N_cells, "Spatial Cells");
}
