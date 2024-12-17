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

// clang-format off
#include "fs_common.h"
#include "ldz_hall.hpp"
// clang-format on
#include <limits>
#include <span>

#ifdef DEBUG_VLASIATOR
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
 * \sa calculateEdgeHallTermComponents
 *
 */
template <typename REAL>
inline REAL JXBX_000_100(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBY, creal BGBZ,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return -(pC[a_zz] * BGBZ) / dz + (pC[a_z] * BGBZ) / dz - (pC[a_yz] * BGBZ) / (2 * dz) -
          (pC[c_xzz] * BGBZ) / (6 * dx) + (pC[c_xz] * BGBZ) / (2 * dx) - (pC[c_xyz] * BGBZ) / (4 * dx) +
          (pC[c_xy] * BGBZ) / (2 * dx) - (pC[c_x] * BGBZ) / dx - (pC[a_yz] * BGBY) / (2 * dy) - (pC[a_yy] * BGBY) / dy +
          (pC[a_y] * BGBY) / dy + (pC[b_xz] * BGBY) / (2 * dx) - (pC[b_xyz] * BGBY) / (4 * dx) -
          (pC[b_xyy] * BGBY) / (6 * dx) + (pC[b_xy] * BGBY) / (2 * dx) - (pC[b_x] * BGBY) / dx -
          (pC[a_zz] * pC[c_zz]) / (6 * dz) + (pC[a_z] * pC[c_zz]) / (6 * dz) - (pC[a_yz] * pC[c_zz]) / (12 * dz) +
          (pC[a_zz] * pC[c_z]) / (2 * dz) - (pC[a_z] * pC[c_z]) / (2 * dz) + (pC[a_yz] * pC[c_z]) / (4 * dz) -
          (pC[a_zz] * pC[c_yz]) / (4 * dz) + (pC[a_z] * pC[c_yz]) / (4 * dz) - (pC[a_yz] * pC[c_yz]) / (8 * dz) +
          (pC[a_zz] * pC[c_y]) / (2 * dz) - (pC[a_z] * pC[c_y]) / (2 * dz) + (pC[a_yz] * pC[c_y]) / (4 * dz) +
          (pC[a_xzz] * pC[c_xz]) / (24 * dz) - (pC[a_xz] * pC[c_xz]) / (24 * dz) + (pC[a_xyz] * pC[c_xz]) / (48 * dz) -
          (pC[a_xzz] * pC[c_x]) / (12 * dz) + (pC[a_xz] * pC[c_x]) / (12 * dz) - (pC[a_xyz] * pC[c_x]) / (24 * dz) -
          (pC[a_zz] * pC[c_0]) / dz + (pC[a_z] * pC[c_0]) / dz - (pC[a_yz] * pC[c_0]) / (2 * dz) +
          (pC[a_yz] * pC[b_z]) / (4 * dy) + (pC[a_yy] * pC[b_z]) / (2 * dy) - (pC[a_y] * pC[b_z]) / (2 * dy) -
          (pC[a_yz] * pC[b_yz]) / (8 * dy) - (pC[a_yy] * pC[b_yz]) / (4 * dy) + (pC[a_y] * pC[b_yz]) / (4 * dy) -
          (pC[a_yz] * pC[b_yy]) / (12 * dy) - (pC[a_yy] * pC[b_yy]) / (6 * dy) + (pC[a_y] * pC[b_yy]) / (6 * dy) +
          (pC[a_yz] * pC[b_y]) / (4 * dy) + (pC[a_yy] * pC[b_y]) / (2 * dy) - (pC[a_y] * pC[b_y]) / (2 * dy) +
          (pC[a_xyz] * pC[b_xy]) / (48 * dy) + (pC[a_xyy] * pC[b_xy]) / (24 * dy) - (pC[a_xy] * pC[b_xy]) / (24 * dy) -
          (pC[a_xyz] * pC[b_x]) / (24 * dy) - (pC[a_xyy] * pC[b_x]) / (12 * dy) + (pC[a_xy] * pC[b_x]) / (12 * dy) -
          (pC[a_yz] * pC[b_0]) / (2 * dy) - (pC[a_yy] * pC[b_0]) / dy + (pC[a_y] * pC[b_0]) / dy -
          (pC[c_xzz] * pC[c_zz]) / (36 * dx) + (pC[c_xz] * pC[c_zz]) / (12 * dx) - (pC[c_xyz] * pC[c_zz]) / (24 * dx) +
          (pC[c_xy] * pC[c_zz]) / (12 * dx) - (pC[c_x] * pC[c_zz]) / (6 * dx) + (pC[c_xzz] * pC[c_z]) / (12 * dx) -
          (pC[c_xz] * pC[c_z]) / (4 * dx) + (pC[c_xyz] * pC[c_z]) / (8 * dx) - (pC[c_xy] * pC[c_z]) / (4 * dx) +
          (pC[c_x] * pC[c_z]) / (2 * dx) - (pC[c_xzz] * pC[c_yz]) / (24 * dx) + (pC[c_xz] * pC[c_yz]) / (8 * dx) -
          (pC[c_xyz] * pC[c_yz]) / (16 * dx) + (pC[c_xy] * pC[c_yz]) / (8 * dx) - (pC[c_x] * pC[c_yz]) / (4 * dx) +
          (pC[c_xzz] * pC[c_y]) / (12 * dx) - (pC[c_xz] * pC[c_y]) / (4 * dx) + (pC[c_xyz] * pC[c_y]) / (8 * dx) -
          (pC[c_xy] * pC[c_y]) / (4 * dx) + (pC[c_x] * pC[c_y]) / (2 * dx) - (pC[c_0] * pC[c_xzz]) / (6 * dx) -
          (pC[c_xxz] * pC[c_xz]) / (24 * dx) + (pC[c_xx] * pC[c_xz]) / (12 * dx) + (pC[c_0] * pC[c_xz]) / (2 * dx) -
          (pC[c_0] * pC[c_xyz]) / (4 * dx) + (pC[c_0] * pC[c_xy]) / (2 * dx) + (pC[c_x] * pC[c_xxz]) / (12 * dx) -
          (pC[c_x] * pC[c_xx]) / (6 * dx) - (pC[c_0] * pC[c_x]) / dx - (pC[b_xz] * pC[b_z]) / (4 * dx) +
          (pC[b_xyz] * pC[b_z]) / (8 * dx) + (pC[b_xyy] * pC[b_z]) / (12 * dx) - (pC[b_xy] * pC[b_z]) / (4 * dx) +
          (pC[b_x] * pC[b_z]) / (2 * dx) + (pC[b_xz] * pC[b_yz]) / (8 * dx) - (pC[b_xyz] * pC[b_yz]) / (16 * dx) -
          (pC[b_xyy] * pC[b_yz]) / (24 * dx) + (pC[b_xy] * pC[b_yz]) / (8 * dx) - (pC[b_x] * pC[b_yz]) / (4 * dx) +
          (pC[b_xz] * pC[b_yy]) / (12 * dx) - (pC[b_xyz] * pC[b_yy]) / (24 * dx) - (pC[b_xyy] * pC[b_yy]) / (36 * dx) +
          (pC[b_xy] * pC[b_yy]) / (12 * dx) - (pC[b_x] * pC[b_yy]) / (6 * dx) - (pC[b_xz] * pC[b_y]) / (4 * dx) +
          (pC[b_xyz] * pC[b_y]) / (8 * dx) + (pC[b_xyy] * pC[b_y]) / (12 * dx) - (pC[b_xy] * pC[b_y]) / (4 * dx) +
          (pC[b_x] * pC[b_y]) / (2 * dx) + (pC[b_0] * pC[b_xz]) / (2 * dx) - (pC[b_0] * pC[b_xyz]) / (4 * dx) -
          (pC[b_0] * pC[b_xyy]) / (6 * dx) - (pC[b_xxy] * pC[b_xy]) / (24 * dx) + (pC[b_xx] * pC[b_xy]) / (12 * dx) +
          (pC[b_0] * pC[b_xy]) / (2 * dx) + (pC[b_x] * pC[b_xxy]) / (12 * dx) - (pC[b_x] * pC[b_xx]) / (6 * dx) -
          (pC[b_0] * pC[b_x]) / dx;
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
 * \sa calculateEdgeHallTermComponents
 *
 */
template <typename REAL>
inline REAL JXBX_010_110(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBY, creal BGBZ,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return -(pC[a_zz] * BGBZ) / dz + (pC[a_z] * BGBZ) / dz + (pC[a_yz] * BGBZ) / (2 * dz) -
          (pC[c_xzz] * BGBZ) / (6 * dx) + (pC[c_xz] * BGBZ) / (2 * dx) + (pC[c_xyz] * BGBZ) / (4 * dx) -
          (pC[c_xy] * BGBZ) / (2 * dx) - (pC[c_x] * BGBZ) / dx - (pC[a_yz] * BGBY) / (2 * dy) + (pC[a_yy] * BGBY) / dy +
          (pC[a_y] * BGBY) / dy + (pC[b_xz] * BGBY) / (2 * dx) + (pC[b_xyz] * BGBY) / (4 * dx) -
          (pC[b_xyy] * BGBY) / (6 * dx) - (pC[b_xy] * BGBY) / (2 * dx) - (pC[b_x] * BGBY) / dx -
          (pC[a_zz] * pC[c_zz]) / (6 * dz) + (pC[a_z] * pC[c_zz]) / (6 * dz) + (pC[a_yz] * pC[c_zz]) / (12 * dz) +
          (pC[a_zz] * pC[c_z]) / (2 * dz) - (pC[a_z] * pC[c_z]) / (2 * dz) - (pC[a_yz] * pC[c_z]) / (4 * dz) +
          (pC[a_zz] * pC[c_yz]) / (4 * dz) - (pC[a_z] * pC[c_yz]) / (4 * dz) - (pC[a_yz] * pC[c_yz]) / (8 * dz) -
          (pC[a_zz] * pC[c_y]) / (2 * dz) + (pC[a_z] * pC[c_y]) / (2 * dz) + (pC[a_yz] * pC[c_y]) / (4 * dz) +
          (pC[a_xzz] * pC[c_xz]) / (24 * dz) - (pC[a_xz] * pC[c_xz]) / (24 * dz) - (pC[a_xyz] * pC[c_xz]) / (48 * dz) -
          (pC[a_xzz] * pC[c_x]) / (12 * dz) + (pC[a_xz] * pC[c_x]) / (12 * dz) + (pC[a_xyz] * pC[c_x]) / (24 * dz) -
          (pC[a_zz] * pC[c_0]) / dz + (pC[a_z] * pC[c_0]) / dz + (pC[a_yz] * pC[c_0]) / (2 * dz) +
          (pC[a_yz] * pC[b_z]) / (4 * dy) - (pC[a_yy] * pC[b_z]) / (2 * dy) - (pC[a_y] * pC[b_z]) / (2 * dy) +
          (pC[a_yz] * pC[b_yz]) / (8 * dy) - (pC[a_yy] * pC[b_yz]) / (4 * dy) - (pC[a_y] * pC[b_yz]) / (4 * dy) -
          (pC[a_yz] * pC[b_yy]) / (12 * dy) + (pC[a_yy] * pC[b_yy]) / (6 * dy) + (pC[a_y] * pC[b_yy]) / (6 * dy) -
          (pC[a_yz] * pC[b_y]) / (4 * dy) + (pC[a_yy] * pC[b_y]) / (2 * dy) + (pC[a_y] * pC[b_y]) / (2 * dy) -
          (pC[a_xyz] * pC[b_xy]) / (48 * dy) + (pC[a_xyy] * pC[b_xy]) / (24 * dy) + (pC[a_xy] * pC[b_xy]) / (24 * dy) -
          (pC[a_xyz] * pC[b_x]) / (24 * dy) + (pC[a_xyy] * pC[b_x]) / (12 * dy) + (pC[a_xy] * pC[b_x]) / (12 * dy) -
          (pC[a_yz] * pC[b_0]) / (2 * dy) + (pC[a_yy] * pC[b_0]) / dy + (pC[a_y] * pC[b_0]) / dy -
          (pC[c_xzz] * pC[c_zz]) / (36 * dx) + (pC[c_xz] * pC[c_zz]) / (12 * dx) + (pC[c_xyz] * pC[c_zz]) / (24 * dx) -
          (pC[c_xy] * pC[c_zz]) / (12 * dx) - (pC[c_x] * pC[c_zz]) / (6 * dx) + (pC[c_xzz] * pC[c_z]) / (12 * dx) -
          (pC[c_xz] * pC[c_z]) / (4 * dx) - (pC[c_xyz] * pC[c_z]) / (8 * dx) + (pC[c_xy] * pC[c_z]) / (4 * dx) +
          (pC[c_x] * pC[c_z]) / (2 * dx) + (pC[c_xzz] * pC[c_yz]) / (24 * dx) - (pC[c_xz] * pC[c_yz]) / (8 * dx) -
          (pC[c_xyz] * pC[c_yz]) / (16 * dx) + (pC[c_xy] * pC[c_yz]) / (8 * dx) + (pC[c_x] * pC[c_yz]) / (4 * dx) -
          (pC[c_xzz] * pC[c_y]) / (12 * dx) + (pC[c_xz] * pC[c_y]) / (4 * dx) + (pC[c_xyz] * pC[c_y]) / (8 * dx) -
          (pC[c_xy] * pC[c_y]) / (4 * dx) - (pC[c_x] * pC[c_y]) / (2 * dx) - (pC[c_0] * pC[c_xzz]) / (6 * dx) -
          (pC[c_xxz] * pC[c_xz]) / (24 * dx) + (pC[c_xx] * pC[c_xz]) / (12 * dx) + (pC[c_0] * pC[c_xz]) / (2 * dx) +
          (pC[c_0] * pC[c_xyz]) / (4 * dx) - (pC[c_0] * pC[c_xy]) / (2 * dx) + (pC[c_x] * pC[c_xxz]) / (12 * dx) -
          (pC[c_x] * pC[c_xx]) / (6 * dx) - (pC[c_0] * pC[c_x]) / dx - (pC[b_xz] * pC[b_z]) / (4 * dx) -
          (pC[b_xyz] * pC[b_z]) / (8 * dx) + (pC[b_xyy] * pC[b_z]) / (12 * dx) + (pC[b_xy] * pC[b_z]) / (4 * dx) +
          (pC[b_x] * pC[b_z]) / (2 * dx) - (pC[b_xz] * pC[b_yz]) / (8 * dx) - (pC[b_xyz] * pC[b_yz]) / (16 * dx) +
          (pC[b_xyy] * pC[b_yz]) / (24 * dx) + (pC[b_xy] * pC[b_yz]) / (8 * dx) + (pC[b_x] * pC[b_yz]) / (4 * dx) +
          (pC[b_xz] * pC[b_yy]) / (12 * dx) + (pC[b_xyz] * pC[b_yy]) / (24 * dx) - (pC[b_xyy] * pC[b_yy]) / (36 * dx) -
          (pC[b_xy] * pC[b_yy]) / (12 * dx) - (pC[b_x] * pC[b_yy]) / (6 * dx) + (pC[b_xz] * pC[b_y]) / (4 * dx) +
          (pC[b_xyz] * pC[b_y]) / (8 * dx) - (pC[b_xyy] * pC[b_y]) / (12 * dx) - (pC[b_xy] * pC[b_y]) / (4 * dx) -
          (pC[b_x] * pC[b_y]) / (2 * dx) + (pC[b_0] * pC[b_xz]) / (2 * dx) + (pC[b_0] * pC[b_xyz]) / (4 * dx) -
          (pC[b_0] * pC[b_xyy]) / (6 * dx) - (pC[b_xxy] * pC[b_xy]) / (24 * dx) - (pC[b_xx] * pC[b_xy]) / (12 * dx) -
          (pC[b_0] * pC[b_xy]) / (2 * dx) - (pC[b_x] * pC[b_xxy]) / (12 * dx) - (pC[b_x] * pC[b_xx]) / (6 * dx) -
          (pC[b_0] * pC[b_x]) / dx;
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
 * \sa calculateEdgeHallTermComponents
 *
 */
template <typename REAL>
inline REAL JXBX_001_101(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBY, creal BGBZ,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return (pC[a_zz] * BGBZ) / dz + (pC[a_z] * BGBZ) / dz - (pC[a_yz] * BGBZ) / (2 * dz) -
          (pC[c_xzz] * BGBZ) / (6 * dx) - (pC[c_xz] * BGBZ) / (2 * dx) + (pC[c_xyz] * BGBZ) / (4 * dx) +
          (pC[c_xy] * BGBZ) / (2 * dx) - (pC[c_x] * BGBZ) / dx + (pC[a_yz] * BGBY) / (2 * dy) - (pC[a_yy] * BGBY) / dy +
          (pC[a_y] * BGBY) / dy - (pC[b_xz] * BGBY) / (2 * dx) + (pC[b_xyz] * BGBY) / (4 * dx) -
          (pC[b_xyy] * BGBY) / (6 * dx) + (pC[b_xy] * BGBY) / (2 * dx) - (pC[b_x] * BGBY) / dx +
          (pC[a_zz] * pC[c_zz]) / (6 * dz) + (pC[a_z] * pC[c_zz]) / (6 * dz) - (pC[a_yz] * pC[c_zz]) / (12 * dz) +
          (pC[a_zz] * pC[c_z]) / (2 * dz) + (pC[a_z] * pC[c_z]) / (2 * dz) - (pC[a_yz] * pC[c_z]) / (4 * dz) -
          (pC[a_zz] * pC[c_yz]) / (4 * dz) - (pC[a_z] * pC[c_yz]) / (4 * dz) + (pC[a_yz] * pC[c_yz]) / (8 * dz) -
          (pC[a_zz] * pC[c_y]) / (2 * dz) - (pC[a_z] * pC[c_y]) / (2 * dz) + (pC[a_yz] * pC[c_y]) / (4 * dz) +
          (pC[a_xzz] * pC[c_xz]) / (24 * dz) + (pC[a_xz] * pC[c_xz]) / (24 * dz) - (pC[a_xyz] * pC[c_xz]) / (48 * dz) +
          (pC[a_xzz] * pC[c_x]) / (12 * dz) + (pC[a_xz] * pC[c_x]) / (12 * dz) - (pC[a_xyz] * pC[c_x]) / (24 * dz) +
          (pC[a_zz] * pC[c_0]) / dz + (pC[a_z] * pC[c_0]) / dz - (pC[a_yz] * pC[c_0]) / (2 * dz) +
          (pC[a_yz] * pC[b_z]) / (4 * dy) - (pC[a_yy] * pC[b_z]) / (2 * dy) + (pC[a_y] * pC[b_z]) / (2 * dy) -
          (pC[a_yz] * pC[b_yz]) / (8 * dy) + (pC[a_yy] * pC[b_yz]) / (4 * dy) - (pC[a_y] * pC[b_yz]) / (4 * dy) +
          (pC[a_yz] * pC[b_yy]) / (12 * dy) - (pC[a_yy] * pC[b_yy]) / (6 * dy) + (pC[a_y] * pC[b_yy]) / (6 * dy) -
          (pC[a_yz] * pC[b_y]) / (4 * dy) + (pC[a_yy] * pC[b_y]) / (2 * dy) - (pC[a_y] * pC[b_y]) / (2 * dy) -
          (pC[a_xyz] * pC[b_xy]) / (48 * dy) + (pC[a_xyy] * pC[b_xy]) / (24 * dy) - (pC[a_xy] * pC[b_xy]) / (24 * dy) +
          (pC[a_xyz] * pC[b_x]) / (24 * dy) - (pC[a_xyy] * pC[b_x]) / (12 * dy) + (pC[a_xy] * pC[b_x]) / (12 * dy) +
          (pC[a_yz] * pC[b_0]) / (2 * dy) - (pC[a_yy] * pC[b_0]) / dy + (pC[a_y] * pC[b_0]) / dy -
          (pC[c_xzz] * pC[c_zz]) / (36 * dx) - (pC[c_xz] * pC[c_zz]) / (12 * dx) + (pC[c_xyz] * pC[c_zz]) / (24 * dx) +
          (pC[c_xy] * pC[c_zz]) / (12 * dx) - (pC[c_x] * pC[c_zz]) / (6 * dx) - (pC[c_xzz] * pC[c_z]) / (12 * dx) -
          (pC[c_xz] * pC[c_z]) / (4 * dx) + (pC[c_xyz] * pC[c_z]) / (8 * dx) + (pC[c_xy] * pC[c_z]) / (4 * dx) -
          (pC[c_x] * pC[c_z]) / (2 * dx) + (pC[c_xzz] * pC[c_yz]) / (24 * dx) + (pC[c_xz] * pC[c_yz]) / (8 * dx) -
          (pC[c_xyz] * pC[c_yz]) / (16 * dx) - (pC[c_xy] * pC[c_yz]) / (8 * dx) + (pC[c_x] * pC[c_yz]) / (4 * dx) +
          (pC[c_xzz] * pC[c_y]) / (12 * dx) + (pC[c_xz] * pC[c_y]) / (4 * dx) - (pC[c_xyz] * pC[c_y]) / (8 * dx) -
          (pC[c_xy] * pC[c_y]) / (4 * dx) + (pC[c_x] * pC[c_y]) / (2 * dx) - (pC[c_0] * pC[c_xzz]) / (6 * dx) -
          (pC[c_xxz] * pC[c_xz]) / (24 * dx) - (pC[c_xx] * pC[c_xz]) / (12 * dx) - (pC[c_0] * pC[c_xz]) / (2 * dx) +
          (pC[c_0] * pC[c_xyz]) / (4 * dx) + (pC[c_0] * pC[c_xy]) / (2 * dx) - (pC[c_x] * pC[c_xxz]) / (12 * dx) -
          (pC[c_x] * pC[c_xx]) / (6 * dx) - (pC[c_0] * pC[c_x]) / dx - (pC[b_xz] * pC[b_z]) / (4 * dx) +
          (pC[b_xyz] * pC[b_z]) / (8 * dx) - (pC[b_xyy] * pC[b_z]) / (12 * dx) + (pC[b_xy] * pC[b_z]) / (4 * dx) -
          (pC[b_x] * pC[b_z]) / (2 * dx) + (pC[b_xz] * pC[b_yz]) / (8 * dx) - (pC[b_xyz] * pC[b_yz]) / (16 * dx) +
          (pC[b_xyy] * pC[b_yz]) / (24 * dx) - (pC[b_xy] * pC[b_yz]) / (8 * dx) + (pC[b_x] * pC[b_yz]) / (4 * dx) -
          (pC[b_xz] * pC[b_yy]) / (12 * dx) + (pC[b_xyz] * pC[b_yy]) / (24 * dx) - (pC[b_xyy] * pC[b_yy]) / (36 * dx) +
          (pC[b_xy] * pC[b_yy]) / (12 * dx) - (pC[b_x] * pC[b_yy]) / (6 * dx) + (pC[b_xz] * pC[b_y]) / (4 * dx) -
          (pC[b_xyz] * pC[b_y]) / (8 * dx) + (pC[b_xyy] * pC[b_y]) / (12 * dx) - (pC[b_xy] * pC[b_y]) / (4 * dx) +
          (pC[b_x] * pC[b_y]) / (2 * dx) - (pC[b_0] * pC[b_xz]) / (2 * dx) + (pC[b_0] * pC[b_xyz]) / (4 * dx) -
          (pC[b_0] * pC[b_xyy]) / (6 * dx) - (pC[b_xxy] * pC[b_xy]) / (24 * dx) + (pC[b_xx] * pC[b_xy]) / (12 * dx) +
          (pC[b_0] * pC[b_xy]) / (2 * dx) + (pC[b_x] * pC[b_xxy]) / (12 * dx) - (pC[b_x] * pC[b_xx]) / (6 * dx) -
          (pC[b_0] * pC[b_x]) / dx;
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
 * \sa calculateEdgeHallTermComponents
 *
 */
template <typename REAL>
inline REAL JXBX_011_111(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBY, creal BGBZ,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return (pC[a_zz] * BGBZ) / dz + (pC[a_z] * BGBZ) / dz + (pC[a_yz] * BGBZ) / (2 * dz) -
          (pC[c_xzz] * BGBZ) / (6 * dx) - (pC[c_xz] * BGBZ) / (2 * dx) - (pC[c_xyz] * BGBZ) / (4 * dx) -
          (pC[c_xy] * BGBZ) / (2 * dx) - (pC[c_x] * BGBZ) / dx + (pC[a_yz] * BGBY) / (2 * dy) + (pC[a_yy] * BGBY) / dy +
          (pC[a_y] * BGBY) / dy - (pC[b_xz] * BGBY) / (2 * dx) - (pC[b_xyz] * BGBY) / (4 * dx) -
          (pC[b_xyy] * BGBY) / (6 * dx) - (pC[b_xy] * BGBY) / (2 * dx) - (pC[b_x] * BGBY) / dx +
          (pC[a_zz] * pC[c_zz]) / (6 * dz) + (pC[a_z] * pC[c_zz]) / (6 * dz) + (pC[a_yz] * pC[c_zz]) / (12 * dz) +
          (pC[a_zz] * pC[c_z]) / (2 * dz) + (pC[a_z] * pC[c_z]) / (2 * dz) + (pC[a_yz] * pC[c_z]) / (4 * dz) +
          (pC[a_zz] * pC[c_yz]) / (4 * dz) + (pC[a_z] * pC[c_yz]) / (4 * dz) + (pC[a_yz] * pC[c_yz]) / (8 * dz) +
          (pC[a_zz] * pC[c_y]) / (2 * dz) + (pC[a_z] * pC[c_y]) / (2 * dz) + (pC[a_yz] * pC[c_y]) / (4 * dz) +
          (pC[a_xzz] * pC[c_xz]) / (24 * dz) + (pC[a_xz] * pC[c_xz]) / (24 * dz) + (pC[a_xyz] * pC[c_xz]) / (48 * dz) +
          (pC[a_xzz] * pC[c_x]) / (12 * dz) + (pC[a_xz] * pC[c_x]) / (12 * dz) + (pC[a_xyz] * pC[c_x]) / (24 * dz) +
          (pC[a_zz] * pC[c_0]) / dz + (pC[a_z] * pC[c_0]) / dz + (pC[a_yz] * pC[c_0]) / (2 * dz) +
          (pC[a_yz] * pC[b_z]) / (4 * dy) + (pC[a_yy] * pC[b_z]) / (2 * dy) + (pC[a_y] * pC[b_z]) / (2 * dy) +
          (pC[a_yz] * pC[b_yz]) / (8 * dy) + (pC[a_yy] * pC[b_yz]) / (4 * dy) + (pC[a_y] * pC[b_yz]) / (4 * dy) +
          (pC[a_yz] * pC[b_yy]) / (12 * dy) + (pC[a_yy] * pC[b_yy]) / (6 * dy) + (pC[a_y] * pC[b_yy]) / (6 * dy) +
          (pC[a_yz] * pC[b_y]) / (4 * dy) + (pC[a_yy] * pC[b_y]) / (2 * dy) + (pC[a_y] * pC[b_y]) / (2 * dy) +
          (pC[a_xyz] * pC[b_xy]) / (48 * dy) + (pC[a_xyy] * pC[b_xy]) / (24 * dy) + (pC[a_xy] * pC[b_xy]) / (24 * dy) +
          (pC[a_xyz] * pC[b_x]) / (24 * dy) + (pC[a_xyy] * pC[b_x]) / (12 * dy) + (pC[a_xy] * pC[b_x]) / (12 * dy) +
          (pC[a_yz] * pC[b_0]) / (2 * dy) + (pC[a_yy] * pC[b_0]) / dy + (pC[a_y] * pC[b_0]) / dy -
          (pC[c_xzz] * pC[c_zz]) / (36 * dx) - (pC[c_xz] * pC[c_zz]) / (12 * dx) - (pC[c_xyz] * pC[c_zz]) / (24 * dx) -
          (pC[c_xy] * pC[c_zz]) / (12 * dx) - (pC[c_x] * pC[c_zz]) / (6 * dx) - (pC[c_xzz] * pC[c_z]) / (12 * dx) -
          (pC[c_xz] * pC[c_z]) / (4 * dx) - (pC[c_xyz] * pC[c_z]) / (8 * dx) - (pC[c_xy] * pC[c_z]) / (4 * dx) -
          (pC[c_x] * pC[c_z]) / (2 * dx) - (pC[c_xzz] * pC[c_yz]) / (24 * dx) - (pC[c_xz] * pC[c_yz]) / (8 * dx) -
          (pC[c_xyz] * pC[c_yz]) / (16 * dx) - (pC[c_xy] * pC[c_yz]) / (8 * dx) - (pC[c_x] * pC[c_yz]) / (4 * dx) -
          (pC[c_xzz] * pC[c_y]) / (12 * dx) - (pC[c_xz] * pC[c_y]) / (4 * dx) - (pC[c_xyz] * pC[c_y]) / (8 * dx) -
          (pC[c_xy] * pC[c_y]) / (4 * dx) - (pC[c_x] * pC[c_y]) / (2 * dx) - (pC[c_0] * pC[c_xzz]) / (6 * dx) -
          (pC[c_xxz] * pC[c_xz]) / (24 * dx) - (pC[c_xx] * pC[c_xz]) / (12 * dx) - (pC[c_0] * pC[c_xz]) / (2 * dx) -
          (pC[c_0] * pC[c_xyz]) / (4 * dx) - (pC[c_0] * pC[c_xy]) / (2 * dx) - (pC[c_x] * pC[c_xxz]) / (12 * dx) -
          (pC[c_x] * pC[c_xx]) / (6 * dx) - (pC[c_0] * pC[c_x]) / dx - (pC[b_xz] * pC[b_z]) / (4 * dx) -
          (pC[b_xyz] * pC[b_z]) / (8 * dx) - (pC[b_xyy] * pC[b_z]) / (12 * dx) - (pC[b_xy] * pC[b_z]) / (4 * dx) -
          (pC[b_x] * pC[b_z]) / (2 * dx) - (pC[b_xz] * pC[b_yz]) / (8 * dx) - (pC[b_xyz] * pC[b_yz]) / (16 * dx) -
          (pC[b_xyy] * pC[b_yz]) / (24 * dx) - (pC[b_xy] * pC[b_yz]) / (8 * dx) - (pC[b_x] * pC[b_yz]) / (4 * dx) -
          (pC[b_xz] * pC[b_yy]) / (12 * dx) - (pC[b_xyz] * pC[b_yy]) / (24 * dx) - (pC[b_xyy] * pC[b_yy]) / (36 * dx) -
          (pC[b_xy] * pC[b_yy]) / (12 * dx) - (pC[b_x] * pC[b_yy]) / (6 * dx) - (pC[b_xz] * pC[b_y]) / (4 * dx) -
          (pC[b_xyz] * pC[b_y]) / (8 * dx) - (pC[b_xyy] * pC[b_y]) / (12 * dx) - (pC[b_xy] * pC[b_y]) / (4 * dx) -
          (pC[b_x] * pC[b_y]) / (2 * dx) - (pC[b_0] * pC[b_xz]) / (2 * dx) - (pC[b_0] * pC[b_xyz]) / (4 * dx) -
          (pC[b_0] * pC[b_xyy]) / (6 * dx) - (pC[b_xxy] * pC[b_xy]) / (24 * dx) - (pC[b_xx] * pC[b_xy]) / (12 * dx) -
          (pC[b_0] * pC[b_xy]) / (2 * dx) - (pC[b_x] * pC[b_xxy]) / (12 * dx) - (pC[b_x] * pC[b_xx]) / (6 * dx) -
          (pC[b_0] * pC[b_x]) / dx;
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
template <typename REAL>
inline REAL JXBY_000_010(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBX, creal BGBZ,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return -(pC[b_zz] * BGBZ) / dz + (pC[b_z] * BGBZ) / dz - (pC[b_xz] * BGBZ) / (2 * dz) -
          (pC[c_yzz] * BGBZ) / (6 * dy) + (pC[c_yz] * BGBZ) / (2 * dy) - (pC[c_y] * BGBZ) / dy -
          (pC[c_xyz] * BGBZ) / (4 * dy) + (pC[c_xy] * BGBZ) / (2 * dy) + (pC[a_yz] * BGBX) / (2 * dy) -
          (pC[a_y] * BGBX) / dy - (pC[a_xyz] * BGBX) / (4 * dy) + (pC[a_xy] * BGBX) / (2 * dy) -
          (pC[a_xxy] * BGBX) / (6 * dy) - (pC[b_xz] * BGBX) / (2 * dx) - (pC[b_xx] * BGBX) / dx +
          (pC[b_x] * BGBX) / dx - (pC[b_zz] * pC[c_zz]) / (6 * dz) + (pC[b_z] * pC[c_zz]) / (6 * dz) -
          (pC[b_xz] * pC[c_zz]) / (12 * dz) + (pC[b_zz] * pC[c_z]) / (2 * dz) - (pC[b_z] * pC[c_z]) / (2 * dz) +
          (pC[b_xz] * pC[c_z]) / (4 * dz) + (pC[b_yzz] * pC[c_yz]) / (24 * dz) - (pC[b_yz] * pC[c_yz]) / (24 * dz) +
          (pC[b_xyz] * pC[c_yz]) / (48 * dz) - (pC[b_yzz] * pC[c_y]) / (12 * dz) + (pC[b_yz] * pC[c_y]) / (12 * dz) -
          (pC[b_xyz] * pC[c_y]) / (24 * dz) - (pC[b_zz] * pC[c_xz]) / (4 * dz) + (pC[b_z] * pC[c_xz]) / (4 * dz) -
          (pC[b_xz] * pC[c_xz]) / (8 * dz) + (pC[b_zz] * pC[c_x]) / (2 * dz) - (pC[b_z] * pC[c_x]) / (2 * dz) +
          (pC[b_xz] * pC[c_x]) / (4 * dz) - (pC[b_zz] * pC[c_0]) / dz + (pC[b_z] * pC[c_0]) / dz -
          (pC[b_xz] * pC[c_0]) / (2 * dz) - (pC[c_yzz] * pC[c_zz]) / (36 * dy) + (pC[c_yz] * pC[c_zz]) / (12 * dy) -
          (pC[c_y] * pC[c_zz]) / (6 * dy) - (pC[c_xyz] * pC[c_zz]) / (24 * dy) + (pC[c_xy] * pC[c_zz]) / (12 * dy) +
          (pC[c_yzz] * pC[c_z]) / (12 * dy) - (pC[c_yz] * pC[c_z]) / (4 * dy) + (pC[c_y] * pC[c_z]) / (2 * dy) +
          (pC[c_xyz] * pC[c_z]) / (8 * dy) - (pC[c_xy] * pC[c_z]) / (4 * dy) - (pC[c_xz] * pC[c_yzz]) / (24 * dy) +
          (pC[c_x] * pC[c_yzz]) / (12 * dy) - (pC[c_0] * pC[c_yzz]) / (6 * dy) - (pC[c_yyz] * pC[c_yz]) / (24 * dy) +
          (pC[c_yy] * pC[c_yz]) / (12 * dy) + (pC[c_xz] * pC[c_yz]) / (8 * dy) - (pC[c_x] * pC[c_yz]) / (4 * dy) +
          (pC[c_0] * pC[c_yz]) / (2 * dy) + (pC[c_y] * pC[c_yyz]) / (12 * dy) - (pC[c_y] * pC[c_yy]) / (6 * dy) -
          (pC[c_xz] * pC[c_y]) / (4 * dy) + (pC[c_x] * pC[c_y]) / (2 * dy) - (pC[c_0] * pC[c_y]) / dy -
          (pC[c_xyz] * pC[c_xz]) / (16 * dy) + (pC[c_xy] * pC[c_xz]) / (8 * dy) + (pC[c_x] * pC[c_xyz]) / (8 * dy) -
          (pC[c_0] * pC[c_xyz]) / (4 * dy) - (pC[c_x] * pC[c_xy]) / (4 * dy) + (pC[c_0] * pC[c_xy]) / (2 * dy) -
          (pC[a_yz] * pC[a_z]) / (4 * dy) + (pC[a_y] * pC[a_z]) / (2 * dy) + (pC[a_xyz] * pC[a_z]) / (8 * dy) -
          (pC[a_xy] * pC[a_z]) / (4 * dy) + (pC[a_xxy] * pC[a_z]) / (12 * dy) + (pC[a_xz] * pC[a_yz]) / (8 * dy) +
          (pC[a_xx] * pC[a_yz]) / (12 * dy) - (pC[a_x] * pC[a_yz]) / (4 * dy) + (pC[a_0] * pC[a_yz]) / (2 * dy) -
          (pC[a_y] * pC[a_yy]) / (6 * dy) + (pC[a_xy] * pC[a_yy]) / (12 * dy) - (pC[a_xz] * pC[a_y]) / (4 * dy) +
          (pC[a_xyy] * pC[a_y]) / (12 * dy) - (pC[a_xx] * pC[a_y]) / (6 * dy) + (pC[a_x] * pC[a_y]) / (2 * dy) -
          (pC[a_0] * pC[a_y]) / dy - (pC[a_xyz] * pC[a_xz]) / (16 * dy) + (pC[a_xy] * pC[a_xz]) / (8 * dy) -
          (pC[a_xxy] * pC[a_xz]) / (24 * dy) - (pC[a_xx] * pC[a_xyz]) / (24 * dy) + (pC[a_x] * pC[a_xyz]) / (8 * dy) -
          (pC[a_0] * pC[a_xyz]) / (4 * dy) - (pC[a_xy] * pC[a_xyy]) / (24 * dy) + (pC[a_xx] * pC[a_xy]) / (12 * dy) -
          (pC[a_x] * pC[a_xy]) / (4 * dy) + (pC[a_0] * pC[a_xy]) / (2 * dy) - (pC[a_xx] * pC[a_xxy]) / (36 * dy) +
          (pC[a_x] * pC[a_xxy]) / (12 * dy) - (pC[a_0] * pC[a_xxy]) / (6 * dy) + (pC[a_z] * pC[b_xz]) / (4 * dx) -
          (pC[a_xz] * pC[b_xz]) / (8 * dx) - (pC[a_xx] * pC[b_xz]) / (12 * dx) + (pC[a_x] * pC[b_xz]) / (4 * dx) -
          (pC[a_0] * pC[b_xz]) / (2 * dx) - (pC[a_y] * pC[b_xyz]) / (24 * dx) + (pC[a_xy] * pC[b_xyz]) / (48 * dx) +
          (pC[a_y] * pC[b_xy]) / (12 * dx) - (pC[a_xy] * pC[b_xy]) / (24 * dx) - (pC[a_y] * pC[b_xxy]) / (12 * dx) +
          (pC[a_xy] * pC[b_xxy]) / (24 * dx) + (pC[a_z] * pC[b_xx]) / (2 * dx) - (pC[a_xz] * pC[b_xx]) / (4 * dx) -
          (pC[a_xx] * pC[b_xx]) / (6 * dx) + (pC[a_x] * pC[b_xx]) / (2 * dx) - (pC[a_0] * pC[b_xx]) / dx -
          (pC[a_z] * pC[b_x]) / (2 * dx) + (pC[a_xz] * pC[b_x]) / (4 * dx) + (pC[a_xx] * pC[b_x]) / (6 * dx) -
          (pC[a_x] * pC[b_x]) / (2 * dx) + (pC[a_0] * pC[b_x]) / dx;
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
template <typename REAL>
inline REAL JXBY_100_110(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBX, creal BGBZ,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return -(pC[b_zz] * BGBZ) / dz + (pC[b_z] * BGBZ) / dz + (pC[b_xz] * BGBZ) / (2 * dz) -
          (pC[c_yzz] * BGBZ) / (6 * dy) + (pC[c_yz] * BGBZ) / (2 * dy) - (pC[c_y] * BGBZ) / dy +
          (pC[c_xyz] * BGBZ) / (4 * dy) - (pC[c_xy] * BGBZ) / (2 * dy) + (pC[a_yz] * BGBX) / (2 * dy) -
          (pC[a_y] * BGBX) / dy + (pC[a_xyz] * BGBX) / (4 * dy) - (pC[a_xy] * BGBX) / (2 * dy) -
          (pC[a_xxy] * BGBX) / (6 * dy) - (pC[b_xz] * BGBX) / (2 * dx) + (pC[b_xx] * BGBX) / dx +
          (pC[b_x] * BGBX) / dx - (pC[b_zz] * pC[c_zz]) / (6 * dz) + (pC[b_z] * pC[c_zz]) / (6 * dz) +
          (pC[b_xz] * pC[c_zz]) / (12 * dz) + (pC[b_zz] * pC[c_z]) / (2 * dz) - (pC[b_z] * pC[c_z]) / (2 * dz) -
          (pC[b_xz] * pC[c_z]) / (4 * dz) + (pC[b_yzz] * pC[c_yz]) / (24 * dz) - (pC[b_yz] * pC[c_yz]) / (24 * dz) -
          (pC[b_xyz] * pC[c_yz]) / (48 * dz) - (pC[b_yzz] * pC[c_y]) / (12 * dz) + (pC[b_yz] * pC[c_y]) / (12 * dz) +
          (pC[b_xyz] * pC[c_y]) / (24 * dz) + (pC[b_zz] * pC[c_xz]) / (4 * dz) - (pC[b_z] * pC[c_xz]) / (4 * dz) -
          (pC[b_xz] * pC[c_xz]) / (8 * dz) - (pC[b_zz] * pC[c_x]) / (2 * dz) + (pC[b_z] * pC[c_x]) / (2 * dz) +
          (pC[b_xz] * pC[c_x]) / (4 * dz) - (pC[b_zz] * pC[c_0]) / dz + (pC[b_z] * pC[c_0]) / dz +
          (pC[b_xz] * pC[c_0]) / (2 * dz) - (pC[c_yzz] * pC[c_zz]) / (36 * dy) + (pC[c_yz] * pC[c_zz]) / (12 * dy) -
          (pC[c_y] * pC[c_zz]) / (6 * dy) + (pC[c_xyz] * pC[c_zz]) / (24 * dy) - (pC[c_xy] * pC[c_zz]) / (12 * dy) +
          (pC[c_yzz] * pC[c_z]) / (12 * dy) - (pC[c_yz] * pC[c_z]) / (4 * dy) + (pC[c_y] * pC[c_z]) / (2 * dy) -
          (pC[c_xyz] * pC[c_z]) / (8 * dy) + (pC[c_xy] * pC[c_z]) / (4 * dy) + (pC[c_xz] * pC[c_yzz]) / (24 * dy) -
          (pC[c_x] * pC[c_yzz]) / (12 * dy) - (pC[c_0] * pC[c_yzz]) / (6 * dy) - (pC[c_yyz] * pC[c_yz]) / (24 * dy) +
          (pC[c_yy] * pC[c_yz]) / (12 * dy) - (pC[c_xz] * pC[c_yz]) / (8 * dy) + (pC[c_x] * pC[c_yz]) / (4 * dy) +
          (pC[c_0] * pC[c_yz]) / (2 * dy) + (pC[c_y] * pC[c_yyz]) / (12 * dy) - (pC[c_y] * pC[c_yy]) / (6 * dy) +
          (pC[c_xz] * pC[c_y]) / (4 * dy) - (pC[c_x] * pC[c_y]) / (2 * dy) - (pC[c_0] * pC[c_y]) / dy -
          (pC[c_xyz] * pC[c_xz]) / (16 * dy) + (pC[c_xy] * pC[c_xz]) / (8 * dy) + (pC[c_x] * pC[c_xyz]) / (8 * dy) +
          (pC[c_0] * pC[c_xyz]) / (4 * dy) - (pC[c_x] * pC[c_xy]) / (4 * dy) - (pC[c_0] * pC[c_xy]) / (2 * dy) -
          (pC[a_yz] * pC[a_z]) / (4 * dy) + (pC[a_y] * pC[a_z]) / (2 * dy) - (pC[a_xyz] * pC[a_z]) / (8 * dy) +
          (pC[a_xy] * pC[a_z]) / (4 * dy) + (pC[a_xxy] * pC[a_z]) / (12 * dy) - (pC[a_xz] * pC[a_yz]) / (8 * dy) +
          (pC[a_xx] * pC[a_yz]) / (12 * dy) + (pC[a_x] * pC[a_yz]) / (4 * dy) + (pC[a_0] * pC[a_yz]) / (2 * dy) -
          (pC[a_y] * pC[a_yy]) / (6 * dy) - (pC[a_xy] * pC[a_yy]) / (12 * dy) + (pC[a_xz] * pC[a_y]) / (4 * dy) -
          (pC[a_xyy] * pC[a_y]) / (12 * dy) - (pC[a_xx] * pC[a_y]) / (6 * dy) - (pC[a_x] * pC[a_y]) / (2 * dy) -
          (pC[a_0] * pC[a_y]) / dy - (pC[a_xyz] * pC[a_xz]) / (16 * dy) + (pC[a_xy] * pC[a_xz]) / (8 * dy) +
          (pC[a_xxy] * pC[a_xz]) / (24 * dy) + (pC[a_xx] * pC[a_xyz]) / (24 * dy) + (pC[a_x] * pC[a_xyz]) / (8 * dy) +
          (pC[a_0] * pC[a_xyz]) / (4 * dy) - (pC[a_xy] * pC[a_xyy]) / (24 * dy) - (pC[a_xx] * pC[a_xy]) / (12 * dy) -
          (pC[a_x] * pC[a_xy]) / (4 * dy) - (pC[a_0] * pC[a_xy]) / (2 * dy) - (pC[a_xx] * pC[a_xxy]) / (36 * dy) -
          (pC[a_x] * pC[a_xxy]) / (12 * dy) - (pC[a_0] * pC[a_xxy]) / (6 * dy) + (pC[a_z] * pC[b_xz]) / (4 * dx) +
          (pC[a_xz] * pC[b_xz]) / (8 * dx) - (pC[a_xx] * pC[b_xz]) / (12 * dx) - (pC[a_x] * pC[b_xz]) / (4 * dx) -
          (pC[a_0] * pC[b_xz]) / (2 * dx) - (pC[a_y] * pC[b_xyz]) / (24 * dx) - (pC[a_xy] * pC[b_xyz]) / (48 * dx) +
          (pC[a_y] * pC[b_xy]) / (12 * dx) + (pC[a_xy] * pC[b_xy]) / (24 * dx) + (pC[a_y] * pC[b_xxy]) / (12 * dx) +
          (pC[a_xy] * pC[b_xxy]) / (24 * dx) - (pC[a_z] * pC[b_xx]) / (2 * dx) - (pC[a_xz] * pC[b_xx]) / (4 * dx) +
          (pC[a_xx] * pC[b_xx]) / (6 * dx) + (pC[a_x] * pC[b_xx]) / (2 * dx) + (pC[a_0] * pC[b_xx]) / dx -
          (pC[a_z] * pC[b_x]) / (2 * dx) - (pC[a_xz] * pC[b_x]) / (4 * dx) + (pC[a_xx] * pC[b_x]) / (6 * dx) +
          (pC[a_x] * pC[b_x]) / (2 * dx) + (pC[a_0] * pC[b_x]) / dx;
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
template <typename REAL>
inline REAL JXBY_001_011(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBX, creal BGBZ,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return (pC[b_zz] * BGBZ) / dz + (pC[b_z] * BGBZ) / dz - (pC[b_xz] * BGBZ) / (2 * dz) -
          (pC[c_yzz] * BGBZ) / (6 * dy) - (pC[c_yz] * BGBZ) / (2 * dy) - (pC[c_y] * BGBZ) / dy +
          (pC[c_xyz] * BGBZ) / (4 * dy) + (pC[c_xy] * BGBZ) / (2 * dy) - (pC[a_yz] * BGBX) / (2 * dy) -
          (pC[a_y] * BGBX) / dy + (pC[a_xyz] * BGBX) / (4 * dy) + (pC[a_xy] * BGBX) / (2 * dy) -
          (pC[a_xxy] * BGBX) / (6 * dy) + (pC[b_xz] * BGBX) / (2 * dx) - (pC[b_xx] * BGBX) / dx +
          (pC[b_x] * BGBX) / dx + (pC[b_zz] * pC[c_zz]) / (6 * dz) + (pC[b_z] * pC[c_zz]) / (6 * dz) -
          (pC[b_xz] * pC[c_zz]) / (12 * dz) + (pC[b_zz] * pC[c_z]) / (2 * dz) + (pC[b_z] * pC[c_z]) / (2 * dz) -
          (pC[b_xz] * pC[c_z]) / (4 * dz) + (pC[b_yzz] * pC[c_yz]) / (24 * dz) + (pC[b_yz] * pC[c_yz]) / (24 * dz) -
          (pC[b_xyz] * pC[c_yz]) / (48 * dz) + (pC[b_yzz] * pC[c_y]) / (12 * dz) + (pC[b_yz] * pC[c_y]) / (12 * dz) -
          (pC[b_xyz] * pC[c_y]) / (24 * dz) - (pC[b_zz] * pC[c_xz]) / (4 * dz) - (pC[b_z] * pC[c_xz]) / (4 * dz) +
          (pC[b_xz] * pC[c_xz]) / (8 * dz) - (pC[b_zz] * pC[c_x]) / (2 * dz) - (pC[b_z] * pC[c_x]) / (2 * dz) +
          (pC[b_xz] * pC[c_x]) / (4 * dz) + (pC[b_zz] * pC[c_0]) / dz + (pC[b_z] * pC[c_0]) / dz -
          (pC[b_xz] * pC[c_0]) / (2 * dz) - (pC[c_yzz] * pC[c_zz]) / (36 * dy) - (pC[c_yz] * pC[c_zz]) / (12 * dy) -
          (pC[c_y] * pC[c_zz]) / (6 * dy) + (pC[c_xyz] * pC[c_zz]) / (24 * dy) + (pC[c_xy] * pC[c_zz]) / (12 * dy) -
          (pC[c_yzz] * pC[c_z]) / (12 * dy) - (pC[c_yz] * pC[c_z]) / (4 * dy) - (pC[c_y] * pC[c_z]) / (2 * dy) +
          (pC[c_xyz] * pC[c_z]) / (8 * dy) + (pC[c_xy] * pC[c_z]) / (4 * dy) + (pC[c_xz] * pC[c_yzz]) / (24 * dy) +
          (pC[c_x] * pC[c_yzz]) / (12 * dy) - (pC[c_0] * pC[c_yzz]) / (6 * dy) - (pC[c_yyz] * pC[c_yz]) / (24 * dy) -
          (pC[c_yy] * pC[c_yz]) / (12 * dy) + (pC[c_xz] * pC[c_yz]) / (8 * dy) + (pC[c_x] * pC[c_yz]) / (4 * dy) -
          (pC[c_0] * pC[c_yz]) / (2 * dy) - (pC[c_y] * pC[c_yyz]) / (12 * dy) - (pC[c_y] * pC[c_yy]) / (6 * dy) +
          (pC[c_xz] * pC[c_y]) / (4 * dy) + (pC[c_x] * pC[c_y]) / (2 * dy) - (pC[c_0] * pC[c_y]) / dy -
          (pC[c_xyz] * pC[c_xz]) / (16 * dy) - (pC[c_xy] * pC[c_xz]) / (8 * dy) - (pC[c_x] * pC[c_xyz]) / (8 * dy) +
          (pC[c_0] * pC[c_xyz]) / (4 * dy) - (pC[c_x] * pC[c_xy]) / (4 * dy) + (pC[c_0] * pC[c_xy]) / (2 * dy) -
          (pC[a_yz] * pC[a_z]) / (4 * dy) - (pC[a_y] * pC[a_z]) / (2 * dy) + (pC[a_xyz] * pC[a_z]) / (8 * dy) +
          (pC[a_xy] * pC[a_z]) / (4 * dy) - (pC[a_xxy] * pC[a_z]) / (12 * dy) + (pC[a_xz] * pC[a_yz]) / (8 * dy) -
          (pC[a_xx] * pC[a_yz]) / (12 * dy) + (pC[a_x] * pC[a_yz]) / (4 * dy) - (pC[a_0] * pC[a_yz]) / (2 * dy) -
          (pC[a_y] * pC[a_yy]) / (6 * dy) + (pC[a_xy] * pC[a_yy]) / (12 * dy) + (pC[a_xz] * pC[a_y]) / (4 * dy) +
          (pC[a_xyy] * pC[a_y]) / (12 * dy) - (pC[a_xx] * pC[a_y]) / (6 * dy) + (pC[a_x] * pC[a_y]) / (2 * dy) -
          (pC[a_0] * pC[a_y]) / dy - (pC[a_xyz] * pC[a_xz]) / (16 * dy) - (pC[a_xy] * pC[a_xz]) / (8 * dy) +
          (pC[a_xxy] * pC[a_xz]) / (24 * dy) + (pC[a_xx] * pC[a_xyz]) / (24 * dy) - (pC[a_x] * pC[a_xyz]) / (8 * dy) +
          (pC[a_0] * pC[a_xyz]) / (4 * dy) - (pC[a_xy] * pC[a_xyy]) / (24 * dy) + (pC[a_xx] * pC[a_xy]) / (12 * dy) -
          (pC[a_x] * pC[a_xy]) / (4 * dy) + (pC[a_0] * pC[a_xy]) / (2 * dy) - (pC[a_xx] * pC[a_xxy]) / (36 * dy) +
          (pC[a_x] * pC[a_xxy]) / (12 * dy) - (pC[a_0] * pC[a_xxy]) / (6 * dy) + (pC[a_z] * pC[b_xz]) / (4 * dx) -
          (pC[a_xz] * pC[b_xz]) / (8 * dx) + (pC[a_xx] * pC[b_xz]) / (12 * dx) - (pC[a_x] * pC[b_xz]) / (4 * dx) +
          (pC[a_0] * pC[b_xz]) / (2 * dx) + (pC[a_y] * pC[b_xyz]) / (24 * dx) - (pC[a_xy] * pC[b_xyz]) / (48 * dx) +
          (pC[a_y] * pC[b_xy]) / (12 * dx) - (pC[a_xy] * pC[b_xy]) / (24 * dx) - (pC[a_y] * pC[b_xxy]) / (12 * dx) +
          (pC[a_xy] * pC[b_xxy]) / (24 * dx) - (pC[a_z] * pC[b_xx]) / (2 * dx) + (pC[a_xz] * pC[b_xx]) / (4 * dx) -
          (pC[a_xx] * pC[b_xx]) / (6 * dx) + (pC[a_x] * pC[b_xx]) / (2 * dx) - (pC[a_0] * pC[b_xx]) / dx +
          (pC[a_z] * pC[b_x]) / (2 * dx) - (pC[a_xz] * pC[b_x]) / (4 * dx) + (pC[a_xx] * pC[b_x]) / (6 * dx) -
          (pC[a_x] * pC[b_x]) / (2 * dx) + (pC[a_0] * pC[b_x]) / dx;
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
template <typename REAL>
inline REAL JXBY_101_111(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBX, creal BGBZ,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return (pC[b_zz] * BGBZ) / dz + (pC[b_z] * BGBZ) / dz + (pC[b_xz] * BGBZ) / (2 * dz) -
          (pC[c_yzz] * BGBZ) / (6 * dy) - (pC[c_yz] * BGBZ) / (2 * dy) - (pC[c_y] * BGBZ) / dy -
          (pC[c_xyz] * BGBZ) / (4 * dy) - (pC[c_xy] * BGBZ) / (2 * dy) - (pC[a_yz] * BGBX) / (2 * dy) -
          (pC[a_y] * BGBX) / dy - (pC[a_xyz] * BGBX) / (4 * dy) - (pC[a_xy] * BGBX) / (2 * dy) -
          (pC[a_xxy] * BGBX) / (6 * dy) + (pC[b_xz] * BGBX) / (2 * dx) + (pC[b_xx] * BGBX) / dx +
          (pC[b_x] * BGBX) / dx + (pC[b_zz] * pC[c_zz]) / (6 * dz) + (pC[b_z] * pC[c_zz]) / (6 * dz) +
          (pC[b_xz] * pC[c_zz]) / (12 * dz) + (pC[b_zz] * pC[c_z]) / (2 * dz) + (pC[b_z] * pC[c_z]) / (2 * dz) +
          (pC[b_xz] * pC[c_z]) / (4 * dz) + (pC[b_yzz] * pC[c_yz]) / (24 * dz) + (pC[b_yz] * pC[c_yz]) / (24 * dz) +
          (pC[b_xyz] * pC[c_yz]) / (48 * dz) + (pC[b_yzz] * pC[c_y]) / (12 * dz) + (pC[b_yz] * pC[c_y]) / (12 * dz) +
          (pC[b_xyz] * pC[c_y]) / (24 * dz) + (pC[b_zz] * pC[c_xz]) / (4 * dz) + (pC[b_z] * pC[c_xz]) / (4 * dz) +
          (pC[b_xz] * pC[c_xz]) / (8 * dz) + (pC[b_zz] * pC[c_x]) / (2 * dz) + (pC[b_z] * pC[c_x]) / (2 * dz) +
          (pC[b_xz] * pC[c_x]) / (4 * dz) + (pC[b_zz] * pC[c_0]) / dz + (pC[b_z] * pC[c_0]) / dz +
          (pC[b_xz] * pC[c_0]) / (2 * dz) - (pC[c_yzz] * pC[c_zz]) / (36 * dy) - (pC[c_yz] * pC[c_zz]) / (12 * dy) -
          (pC[c_y] * pC[c_zz]) / (6 * dy) - (pC[c_xyz] * pC[c_zz]) / (24 * dy) - (pC[c_xy] * pC[c_zz]) / (12 * dy) -
          (pC[c_yzz] * pC[c_z]) / (12 * dy) - (pC[c_yz] * pC[c_z]) / (4 * dy) - (pC[c_y] * pC[c_z]) / (2 * dy) -
          (pC[c_xyz] * pC[c_z]) / (8 * dy) - (pC[c_xy] * pC[c_z]) / (4 * dy) - (pC[c_xz] * pC[c_yzz]) / (24 * dy) -
          (pC[c_x] * pC[c_yzz]) / (12 * dy) - (pC[c_0] * pC[c_yzz]) / (6 * dy) - (pC[c_yyz] * pC[c_yz]) / (24 * dy) -
          (pC[c_yy] * pC[c_yz]) / (12 * dy) - (pC[c_xz] * pC[c_yz]) / (8 * dy) - (pC[c_x] * pC[c_yz]) / (4 * dy) -
          (pC[c_0] * pC[c_yz]) / (2 * dy) - (pC[c_y] * pC[c_yyz]) / (12 * dy) - (pC[c_y] * pC[c_yy]) / (6 * dy) -
          (pC[c_xz] * pC[c_y]) / (4 * dy) - (pC[c_x] * pC[c_y]) / (2 * dy) - (pC[c_0] * pC[c_y]) / dy -
          (pC[c_xyz] * pC[c_xz]) / (16 * dy) - (pC[c_xy] * pC[c_xz]) / (8 * dy) - (pC[c_x] * pC[c_xyz]) / (8 * dy) -
          (pC[c_0] * pC[c_xyz]) / (4 * dy) - (pC[c_x] * pC[c_xy]) / (4 * dy) - (pC[c_0] * pC[c_xy]) / (2 * dy) -
          (pC[a_yz] * pC[a_z]) / (4 * dy) - (pC[a_y] * pC[a_z]) / (2 * dy) - (pC[a_xyz] * pC[a_z]) / (8 * dy) -
          (pC[a_xy] * pC[a_z]) / (4 * dy) - (pC[a_xxy] * pC[a_z]) / (12 * dy) - (pC[a_xz] * pC[a_yz]) / (8 * dy) -
          (pC[a_xx] * pC[a_yz]) / (12 * dy) - (pC[a_x] * pC[a_yz]) / (4 * dy) - (pC[a_0] * pC[a_yz]) / (2 * dy) -
          (pC[a_y] * pC[a_yy]) / (6 * dy) - (pC[a_xy] * pC[a_yy]) / (12 * dy) - (pC[a_xz] * pC[a_y]) / (4 * dy) -
          (pC[a_xyy] * pC[a_y]) / (12 * dy) - (pC[a_xx] * pC[a_y]) / (6 * dy) - (pC[a_x] * pC[a_y]) / (2 * dy) -
          (pC[a_0] * pC[a_y]) / dy - (pC[a_xyz] * pC[a_xz]) / (16 * dy) - (pC[a_xy] * pC[a_xz]) / (8 * dy) -
          (pC[a_xxy] * pC[a_xz]) / (24 * dy) - (pC[a_xx] * pC[a_xyz]) / (24 * dy) - (pC[a_x] * pC[a_xyz]) / (8 * dy) -
          (pC[a_0] * pC[a_xyz]) / (4 * dy) - (pC[a_xy] * pC[a_xyy]) / (24 * dy) - (pC[a_xx] * pC[a_xy]) / (12 * dy) -
          (pC[a_x] * pC[a_xy]) / (4 * dy) - (pC[a_0] * pC[a_xy]) / (2 * dy) - (pC[a_xx] * pC[a_xxy]) / (36 * dy) -
          (pC[a_x] * pC[a_xxy]) / (12 * dy) - (pC[a_0] * pC[a_xxy]) / (6 * dy) + (pC[a_z] * pC[b_xz]) / (4 * dx) +
          (pC[a_xz] * pC[b_xz]) / (8 * dx) + (pC[a_xx] * pC[b_xz]) / (12 * dx) + (pC[a_x] * pC[b_xz]) / (4 * dx) +
          (pC[a_0] * pC[b_xz]) / (2 * dx) + (pC[a_y] * pC[b_xyz]) / (24 * dx) + (pC[a_xy] * pC[b_xyz]) / (48 * dx) +
          (pC[a_y] * pC[b_xy]) / (12 * dx) + (pC[a_xy] * pC[b_xy]) / (24 * dx) + (pC[a_y] * pC[b_xxy]) / (12 * dx) +
          (pC[a_xy] * pC[b_xxy]) / (24 * dx) + (pC[a_z] * pC[b_xx]) / (2 * dx) + (pC[a_xz] * pC[b_xx]) / (4 * dx) +
          (pC[a_xx] * pC[b_xx]) / (6 * dx) + (pC[a_x] * pC[b_xx]) / (2 * dx) + (pC[a_0] * pC[b_xx]) / dx +
          (pC[a_z] * pC[b_x]) / (2 * dx) + (pC[a_xz] * pC[b_x]) / (4 * dx) + (pC[a_xx] * pC[b_x]) / (6 * dx) +
          (pC[a_x] * pC[b_x]) / (2 * dx) + (pC[a_0] * pC[b_x]) / dx;
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
template <typename REAL>
inline REAL JXBZ_000_001(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBX, creal BGBY,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return -(pC[b_z] * BGBY) / dz + (pC[b_yz] * BGBY) / (2 * dz) - (pC[b_yyz] * BGBY) / (6 * dz) +
          (pC[b_xz] * BGBY) / (2 * dz) - (pC[b_xyz] * BGBY) / (4 * dz) - (pC[c_yy] * BGBY) / dy +
          (pC[c_y] * BGBY) / dy - (pC[c_xy] * BGBY) / (2 * dy) - (pC[a_z] * BGBX) / dz + (pC[a_yz] * BGBX) / (2 * dz) +
          (pC[a_xz] * BGBX) / (2 * dz) - (pC[a_xyz] * BGBX) / (4 * dz) - (pC[a_xxz] * BGBX) / (6 * dz) -
          (pC[c_xy] * BGBX) / (2 * dx) - (pC[c_xx] * BGBX) / dx + (pC[c_x] * BGBX) / dx -
          (pC[b_z] * pC[b_zz]) / (6 * dz) + (pC[b_yz] * pC[b_zz]) / (12 * dz) + (pC[b_yzz] * pC[b_z]) / (12 * dz) -
          (pC[b_yy] * pC[b_z]) / (6 * dz) + (pC[b_y] * pC[b_z]) / (2 * dz) - (pC[b_xy] * pC[b_z]) / (4 * dz) +
          (pC[b_x] * pC[b_z]) / (2 * dz) - (pC[b_0] * pC[b_z]) / dz - (pC[b_yz] * pC[b_yzz]) / (24 * dz) +
          (pC[b_yy] * pC[b_yz]) / (12 * dz) - (pC[b_y] * pC[b_yz]) / (4 * dz) + (pC[b_xy] * pC[b_yz]) / (8 * dz) -
          (pC[b_x] * pC[b_yz]) / (4 * dz) + (pC[b_0] * pC[b_yz]) / (2 * dz) - (pC[b_yy] * pC[b_yyz]) / (36 * dz) +
          (pC[b_y] * pC[b_yyz]) / (12 * dz) - (pC[b_xy] * pC[b_yyz]) / (24 * dz) + (pC[b_x] * pC[b_yyz]) / (12 * dz) -
          (pC[b_0] * pC[b_yyz]) / (6 * dz) + (pC[b_xz] * pC[b_yy]) / (12 * dz) - (pC[b_xyz] * pC[b_yy]) / (24 * dz) -
          (pC[b_xz] * pC[b_y]) / (4 * dz) + (pC[b_xyz] * pC[b_y]) / (8 * dz) + (pC[b_xy] * pC[b_xz]) / (8 * dz) -
          (pC[b_x] * pC[b_xz]) / (4 * dz) + (pC[b_0] * pC[b_xz]) / (2 * dz) - (pC[b_xy] * pC[b_xyz]) / (16 * dz) +
          (pC[b_x] * pC[b_xyz]) / (8 * dz) - (pC[b_0] * pC[b_xyz]) / (4 * dz) - (pC[a_z] * pC[a_zz]) / (6 * dz) +
          (pC[a_xz] * pC[a_zz]) / (12 * dz) + (pC[a_y] * pC[a_z]) / (2 * dz) + (pC[a_xzz] * pC[a_z]) / (12 * dz) -
          (pC[a_xy] * pC[a_z]) / (4 * dz) - (pC[a_xx] * pC[a_z]) / (6 * dz) + (pC[a_x] * pC[a_z]) / (2 * dz) -
          (pC[a_0] * pC[a_z]) / dz - (pC[a_y] * pC[a_yz]) / (4 * dz) + (pC[a_xy] * pC[a_yz]) / (8 * dz) +
          (pC[a_xx] * pC[a_yz]) / (12 * dz) - (pC[a_x] * pC[a_yz]) / (4 * dz) + (pC[a_0] * pC[a_yz]) / (2 * dz) -
          (pC[a_xz] * pC[a_y]) / (4 * dz) + (pC[a_xyz] * pC[a_y]) / (8 * dz) + (pC[a_xxz] * pC[a_y]) / (12 * dz) -
          (pC[a_xz] * pC[a_xzz]) / (24 * dz) + (pC[a_xy] * pC[a_xz]) / (8 * dz) + (pC[a_xx] * pC[a_xz]) / (12 * dz) -
          (pC[a_x] * pC[a_xz]) / (4 * dz) + (pC[a_0] * pC[a_xz]) / (2 * dz) - (pC[a_xy] * pC[a_xyz]) / (16 * dz) -
          (pC[a_xx] * pC[a_xyz]) / (24 * dz) + (pC[a_x] * pC[a_xyz]) / (8 * dz) - (pC[a_0] * pC[a_xyz]) / (4 * dz) -
          (pC[a_xxz] * pC[a_xy]) / (24 * dz) - (pC[a_xx] * pC[a_xxz]) / (36 * dz) + (pC[a_x] * pC[a_xxz]) / (12 * dz) -
          (pC[a_0] * pC[a_xxz]) / (6 * dz) + (pC[b_z] * pC[c_yz]) / (12 * dy) - (pC[b_yz] * pC[c_yz]) / (24 * dy) -
          (pC[b_z] * pC[c_yyz]) / (12 * dy) + (pC[b_yz] * pC[c_yyz]) / (24 * dy) - (pC[b_yy] * pC[c_yy]) / (6 * dy) +
          (pC[b_y] * pC[c_yy]) / (2 * dy) - (pC[b_xy] * pC[c_yy]) / (4 * dy) + (pC[b_x] * pC[c_yy]) / (2 * dy) -
          (pC[b_0] * pC[c_yy]) / dy + (pC[b_yy] * pC[c_y]) / (6 * dy) - (pC[b_y] * pC[c_y]) / (2 * dy) +
          (pC[b_xy] * pC[c_y]) / (4 * dy) - (pC[b_x] * pC[c_y]) / (2 * dy) + (pC[b_0] * pC[c_y]) / dy -
          (pC[b_z] * pC[c_xyz]) / (24 * dy) + (pC[b_yz] * pC[c_xyz]) / (48 * dy) - (pC[b_yy] * pC[c_xy]) / (12 * dy) +
          (pC[b_y] * pC[c_xy]) / (4 * dy) - (pC[b_xy] * pC[c_xy]) / (8 * dy) + (pC[b_x] * pC[c_xy]) / (4 * dy) -
          (pC[b_0] * pC[c_xy]) / (2 * dy) + (pC[a_z] * pC[c_xz]) / (12 * dx) - (pC[a_xz] * pC[c_xz]) / (24 * dx) -
          (pC[a_z] * pC[c_xyz]) / (24 * dx) + (pC[a_xz] * pC[c_xyz]) / (48 * dx) + (pC[a_y] * pC[c_xy]) / (4 * dx) -
          (pC[a_xy] * pC[c_xy]) / (8 * dx) - (pC[a_xx] * pC[c_xy]) / (12 * dx) + (pC[a_x] * pC[c_xy]) / (4 * dx) -
          (pC[a_0] * pC[c_xy]) / (2 * dx) - (pC[a_z] * pC[c_xxz]) / (12 * dx) + (pC[a_xz] * pC[c_xxz]) / (24 * dx) +
          (pC[a_y] * pC[c_xx]) / (2 * dx) - (pC[a_xy] * pC[c_xx]) / (4 * dx) - (pC[a_xx] * pC[c_xx]) / (6 * dx) +
          (pC[a_x] * pC[c_xx]) / (2 * dx) - (pC[a_0] * pC[c_xx]) / dx - (pC[a_y] * pC[c_x]) / (2 * dx) +
          (pC[a_xy] * pC[c_x]) / (4 * dx) + (pC[a_xx] * pC[c_x]) / (6 * dx) - (pC[a_x] * pC[c_x]) / (2 * dx) +
          (pC[a_0] * pC[c_x]) / dx;
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
template <typename REAL>
inline REAL JXBZ_100_101(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBX, creal BGBY,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return -(pC[b_z] * BGBY) / dz + (pC[b_yz] * BGBY) / (2 * dz) - (pC[b_yyz] * BGBY) / (6 * dz) -
          (pC[b_xz] * BGBY) / (2 * dz) + (pC[b_xyz] * BGBY) / (4 * dz) - (pC[c_yy] * BGBY) / dy +
          (pC[c_y] * BGBY) / dy + (pC[c_xy] * BGBY) / (2 * dy) - (pC[a_z] * BGBX) / dz + (pC[a_yz] * BGBX) / (2 * dz) -
          (pC[a_xz] * BGBX) / (2 * dz) + (pC[a_xyz] * BGBX) / (4 * dz) - (pC[a_xxz] * BGBX) / (6 * dz) -
          (pC[c_xy] * BGBX) / (2 * dx) + (pC[c_xx] * BGBX) / dx + (pC[c_x] * BGBX) / dx -
          (pC[b_z] * pC[b_zz]) / (6 * dz) + (pC[b_yz] * pC[b_zz]) / (12 * dz) + (pC[b_yzz] * pC[b_z]) / (12 * dz) -
          (pC[b_yy] * pC[b_z]) / (6 * dz) + (pC[b_y] * pC[b_z]) / (2 * dz) + (pC[b_xy] * pC[b_z]) / (4 * dz) -
          (pC[b_x] * pC[b_z]) / (2 * dz) - (pC[b_0] * pC[b_z]) / dz - (pC[b_yz] * pC[b_yzz]) / (24 * dz) +
          (pC[b_yy] * pC[b_yz]) / (12 * dz) - (pC[b_y] * pC[b_yz]) / (4 * dz) - (pC[b_xy] * pC[b_yz]) / (8 * dz) +
          (pC[b_x] * pC[b_yz]) / (4 * dz) + (pC[b_0] * pC[b_yz]) / (2 * dz) - (pC[b_yy] * pC[b_yyz]) / (36 * dz) +
          (pC[b_y] * pC[b_yyz]) / (12 * dz) + (pC[b_xy] * pC[b_yyz]) / (24 * dz) - (pC[b_x] * pC[b_yyz]) / (12 * dz) -
          (pC[b_0] * pC[b_yyz]) / (6 * dz) - (pC[b_xz] * pC[b_yy]) / (12 * dz) + (pC[b_xyz] * pC[b_yy]) / (24 * dz) +
          (pC[b_xz] * pC[b_y]) / (4 * dz) - (pC[b_xyz] * pC[b_y]) / (8 * dz) + (pC[b_xy] * pC[b_xz]) / (8 * dz) -
          (pC[b_x] * pC[b_xz]) / (4 * dz) - (pC[b_0] * pC[b_xz]) / (2 * dz) - (pC[b_xy] * pC[b_xyz]) / (16 * dz) +
          (pC[b_x] * pC[b_xyz]) / (8 * dz) + (pC[b_0] * pC[b_xyz]) / (4 * dz) - (pC[a_z] * pC[a_zz]) / (6 * dz) -
          (pC[a_xz] * pC[a_zz]) / (12 * dz) + (pC[a_y] * pC[a_z]) / (2 * dz) - (pC[a_xzz] * pC[a_z]) / (12 * dz) +
          (pC[a_xy] * pC[a_z]) / (4 * dz) - (pC[a_xx] * pC[a_z]) / (6 * dz) - (pC[a_x] * pC[a_z]) / (2 * dz) -
          (pC[a_0] * pC[a_z]) / dz - (pC[a_y] * pC[a_yz]) / (4 * dz) - (pC[a_xy] * pC[a_yz]) / (8 * dz) +
          (pC[a_xx] * pC[a_yz]) / (12 * dz) + (pC[a_x] * pC[a_yz]) / (4 * dz) + (pC[a_0] * pC[a_yz]) / (2 * dz) +
          (pC[a_xz] * pC[a_y]) / (4 * dz) - (pC[a_xyz] * pC[a_y]) / (8 * dz) + (pC[a_xxz] * pC[a_y]) / (12 * dz) -
          (pC[a_xz] * pC[a_xzz]) / (24 * dz) + (pC[a_xy] * pC[a_xz]) / (8 * dz) - (pC[a_xx] * pC[a_xz]) / (12 * dz) -
          (pC[a_x] * pC[a_xz]) / (4 * dz) - (pC[a_0] * pC[a_xz]) / (2 * dz) - (pC[a_xy] * pC[a_xyz]) / (16 * dz) +
          (pC[a_xx] * pC[a_xyz]) / (24 * dz) + (pC[a_x] * pC[a_xyz]) / (8 * dz) + (pC[a_0] * pC[a_xyz]) / (4 * dz) +
          (pC[a_xxz] * pC[a_xy]) / (24 * dz) - (pC[a_xx] * pC[a_xxz]) / (36 * dz) - (pC[a_x] * pC[a_xxz]) / (12 * dz) -
          (pC[a_0] * pC[a_xxz]) / (6 * dz) + (pC[b_z] * pC[c_yz]) / (12 * dy) - (pC[b_yz] * pC[c_yz]) / (24 * dy) -
          (pC[b_z] * pC[c_yyz]) / (12 * dy) + (pC[b_yz] * pC[c_yyz]) / (24 * dy) - (pC[b_yy] * pC[c_yy]) / (6 * dy) +
          (pC[b_y] * pC[c_yy]) / (2 * dy) + (pC[b_xy] * pC[c_yy]) / (4 * dy) - (pC[b_x] * pC[c_yy]) / (2 * dy) -
          (pC[b_0] * pC[c_yy]) / dy + (pC[b_yy] * pC[c_y]) / (6 * dy) - (pC[b_y] * pC[c_y]) / (2 * dy) -
          (pC[b_xy] * pC[c_y]) / (4 * dy) + (pC[b_x] * pC[c_y]) / (2 * dy) + (pC[b_0] * pC[c_y]) / dy +
          (pC[b_z] * pC[c_xyz]) / (24 * dy) - (pC[b_yz] * pC[c_xyz]) / (48 * dy) + (pC[b_yy] * pC[c_xy]) / (12 * dy) -
          (pC[b_y] * pC[c_xy]) / (4 * dy) - (pC[b_xy] * pC[c_xy]) / (8 * dy) + (pC[b_x] * pC[c_xy]) / (4 * dy) +
          (pC[b_0] * pC[c_xy]) / (2 * dy) + (pC[a_z] * pC[c_xz]) / (12 * dx) + (pC[a_xz] * pC[c_xz]) / (24 * dx) -
          (pC[a_z] * pC[c_xyz]) / (24 * dx) - (pC[a_xz] * pC[c_xyz]) / (48 * dx) + (pC[a_y] * pC[c_xy]) / (4 * dx) +
          (pC[a_xy] * pC[c_xy]) / (8 * dx) - (pC[a_xx] * pC[c_xy]) / (12 * dx) - (pC[a_x] * pC[c_xy]) / (4 * dx) -
          (pC[a_0] * pC[c_xy]) / (2 * dx) + (pC[a_z] * pC[c_xxz]) / (12 * dx) + (pC[a_xz] * pC[c_xxz]) / (24 * dx) -
          (pC[a_y] * pC[c_xx]) / (2 * dx) - (pC[a_xy] * pC[c_xx]) / (4 * dx) + (pC[a_xx] * pC[c_xx]) / (6 * dx) +
          (pC[a_x] * pC[c_xx]) / (2 * dx) + (pC[a_0] * pC[c_xx]) / dx - (pC[a_y] * pC[c_x]) / (2 * dx) -
          (pC[a_xy] * pC[c_x]) / (4 * dx) + (pC[a_xx] * pC[c_x]) / (6 * dx) + (pC[a_x] * pC[c_x]) / (2 * dx) +
          (pC[a_0] * pC[c_x]) / dx;
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
template <typename REAL>
inline REAL JXBZ_010_011(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBX, creal BGBY,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return -(pC[b_z] * BGBY) / dz - (pC[b_yz] * BGBY) / (2 * dz) - (pC[b_yyz] * BGBY) / (6 * dz) +
          (pC[b_xz] * BGBY) / (2 * dz) + (pC[b_xyz] * BGBY) / (4 * dz) + (pC[c_yy] * BGBY) / dy +
          (pC[c_y] * BGBY) / dy - (pC[c_xy] * BGBY) / (2 * dy) - (pC[a_z] * BGBX) / dz - (pC[a_yz] * BGBX) / (2 * dz) +
          (pC[a_xz] * BGBX) / (2 * dz) + (pC[a_xyz] * BGBX) / (4 * dz) - (pC[a_xxz] * BGBX) / (6 * dz) +
          (pC[c_xy] * BGBX) / (2 * dx) - (pC[c_xx] * BGBX) / dx + (pC[c_x] * BGBX) / dx -
          (pC[b_z] * pC[b_zz]) / (6 * dz) - (pC[b_yz] * pC[b_zz]) / (12 * dz) - (pC[b_yzz] * pC[b_z]) / (12 * dz) -
          (pC[b_yy] * pC[b_z]) / (6 * dz) - (pC[b_y] * pC[b_z]) / (2 * dz) + (pC[b_xy] * pC[b_z]) / (4 * dz) +
          (pC[b_x] * pC[b_z]) / (2 * dz) - (pC[b_0] * pC[b_z]) / dz - (pC[b_yz] * pC[b_yzz]) / (24 * dz) -
          (pC[b_yy] * pC[b_yz]) / (12 * dz) - (pC[b_y] * pC[b_yz]) / (4 * dz) + (pC[b_xy] * pC[b_yz]) / (8 * dz) +
          (pC[b_x] * pC[b_yz]) / (4 * dz) - (pC[b_0] * pC[b_yz]) / (2 * dz) - (pC[b_yy] * pC[b_yyz]) / (36 * dz) -
          (pC[b_y] * pC[b_yyz]) / (12 * dz) + (pC[b_xy] * pC[b_yyz]) / (24 * dz) + (pC[b_x] * pC[b_yyz]) / (12 * dz) -
          (pC[b_0] * pC[b_yyz]) / (6 * dz) + (pC[b_xz] * pC[b_yy]) / (12 * dz) + (pC[b_xyz] * pC[b_yy]) / (24 * dz) +
          (pC[b_xz] * pC[b_y]) / (4 * dz) + (pC[b_xyz] * pC[b_y]) / (8 * dz) - (pC[b_xy] * pC[b_xz]) / (8 * dz) -
          (pC[b_x] * pC[b_xz]) / (4 * dz) + (pC[b_0] * pC[b_xz]) / (2 * dz) - (pC[b_xy] * pC[b_xyz]) / (16 * dz) -
          (pC[b_x] * pC[b_xyz]) / (8 * dz) + (pC[b_0] * pC[b_xyz]) / (4 * dz) - (pC[a_z] * pC[a_zz]) / (6 * dz) +
          (pC[a_xz] * pC[a_zz]) / (12 * dz) - (pC[a_y] * pC[a_z]) / (2 * dz) + (pC[a_xzz] * pC[a_z]) / (12 * dz) +
          (pC[a_xy] * pC[a_z]) / (4 * dz) - (pC[a_xx] * pC[a_z]) / (6 * dz) + (pC[a_x] * pC[a_z]) / (2 * dz) -
          (pC[a_0] * pC[a_z]) / dz - (pC[a_y] * pC[a_yz]) / (4 * dz) + (pC[a_xy] * pC[a_yz]) / (8 * dz) -
          (pC[a_xx] * pC[a_yz]) / (12 * dz) + (pC[a_x] * pC[a_yz]) / (4 * dz) - (pC[a_0] * pC[a_yz]) / (2 * dz) +
          (pC[a_xz] * pC[a_y]) / (4 * dz) + (pC[a_xyz] * pC[a_y]) / (8 * dz) - (pC[a_xxz] * pC[a_y]) / (12 * dz) -
          (pC[a_xz] * pC[a_xzz]) / (24 * dz) - (pC[a_xy] * pC[a_xz]) / (8 * dz) + (pC[a_xx] * pC[a_xz]) / (12 * dz) -
          (pC[a_x] * pC[a_xz]) / (4 * dz) + (pC[a_0] * pC[a_xz]) / (2 * dz) - (pC[a_xy] * pC[a_xyz]) / (16 * dz) +
          (pC[a_xx] * pC[a_xyz]) / (24 * dz) - (pC[a_x] * pC[a_xyz]) / (8 * dz) + (pC[a_0] * pC[a_xyz]) / (4 * dz) +
          (pC[a_xxz] * pC[a_xy]) / (24 * dz) - (pC[a_xx] * pC[a_xxz]) / (36 * dz) + (pC[a_x] * pC[a_xxz]) / (12 * dz) -
          (pC[a_0] * pC[a_xxz]) / (6 * dz) + (pC[b_z] * pC[c_yz]) / (12 * dy) + (pC[b_yz] * pC[c_yz]) / (24 * dy) +
          (pC[b_z] * pC[c_yyz]) / (12 * dy) + (pC[b_yz] * pC[c_yyz]) / (24 * dy) + (pC[b_yy] * pC[c_yy]) / (6 * dy) +
          (pC[b_y] * pC[c_yy]) / (2 * dy) - (pC[b_xy] * pC[c_yy]) / (4 * dy) - (pC[b_x] * pC[c_yy]) / (2 * dy) +
          (pC[b_0] * pC[c_yy]) / dy + (pC[b_yy] * pC[c_y]) / (6 * dy) + (pC[b_y] * pC[c_y]) / (2 * dy) -
          (pC[b_xy] * pC[c_y]) / (4 * dy) - (pC[b_x] * pC[c_y]) / (2 * dy) + (pC[b_0] * pC[c_y]) / dy -
          (pC[b_z] * pC[c_xyz]) / (24 * dy) - (pC[b_yz] * pC[c_xyz]) / (48 * dy) - (pC[b_yy] * pC[c_xy]) / (12 * dy) -
          (pC[b_y] * pC[c_xy]) / (4 * dy) + (pC[b_xy] * pC[c_xy]) / (8 * dy) + (pC[b_x] * pC[c_xy]) / (4 * dy) -
          (pC[b_0] * pC[c_xy]) / (2 * dy) + (pC[a_z] * pC[c_xz]) / (12 * dx) - (pC[a_xz] * pC[c_xz]) / (24 * dx) +
          (pC[a_z] * pC[c_xyz]) / (24 * dx) - (pC[a_xz] * pC[c_xyz]) / (48 * dx) + (pC[a_y] * pC[c_xy]) / (4 * dx) -
          (pC[a_xy] * pC[c_xy]) / (8 * dx) + (pC[a_xx] * pC[c_xy]) / (12 * dx) - (pC[a_x] * pC[c_xy]) / (4 * dx) +
          (pC[a_0] * pC[c_xy]) / (2 * dx) - (pC[a_z] * pC[c_xxz]) / (12 * dx) + (pC[a_xz] * pC[c_xxz]) / (24 * dx) -
          (pC[a_y] * pC[c_xx]) / (2 * dx) + (pC[a_xy] * pC[c_xx]) / (4 * dx) - (pC[a_xx] * pC[c_xx]) / (6 * dx) +
          (pC[a_x] * pC[c_xx]) / (2 * dx) - (pC[a_0] * pC[c_xx]) / dx + (pC[a_y] * pC[c_x]) / (2 * dx) -
          (pC[a_xy] * pC[c_x]) / (4 * dx) + (pC[a_xx] * pC[c_x]) / (6 * dx) - (pC[a_x] * pC[c_x]) / (2 * dx) +
          (pC[a_0] * pC[c_x]) / dx;
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
template <typename REAL>
inline REAL JXBZ_110_111(const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, creal BGBX, creal BGBY,
                         const std::array<Real, 3>& gridSpacing) {
   using namespace Rec;
   const auto dx = gridSpacing[0];
   const auto dy = gridSpacing[1];
   const auto dz = gridSpacing[2];
   return -(pC[b_z] * BGBY) / dz - (pC[b_yz] * BGBY) / (2 * dz) - (pC[b_yyz] * BGBY) / (6 * dz) -
          (pC[b_xz] * BGBY) / (2 * dz) - (pC[b_xyz] * BGBY) / (4 * dz) + (pC[c_yy] * BGBY) / dy +
          (pC[c_y] * BGBY) / dy + (pC[c_xy] * BGBY) / (2 * dy) - (pC[a_z] * BGBX) / dz - (pC[a_yz] * BGBX) / (2 * dz) -
          (pC[a_xz] * BGBX) / (2 * dz) - (pC[a_xyz] * BGBX) / (4 * dz) - (pC[a_xxz] * BGBX) / (6 * dz) +
          (pC[c_xy] * BGBX) / (2 * dx) + (pC[c_xx] * BGBX) / dx + (pC[c_x] * BGBX) / dx -
          (pC[b_z] * pC[b_zz]) / (6 * dz) - (pC[b_yz] * pC[b_zz]) / (12 * dz) - (pC[b_yzz] * pC[b_z]) / (12 * dz) -
          (pC[b_yy] * pC[b_z]) / (6 * dz) - (pC[b_y] * pC[b_z]) / (2 * dz) - (pC[b_xy] * pC[b_z]) / (4 * dz) -
          (pC[b_x] * pC[b_z]) / (2 * dz) - (pC[b_0] * pC[b_z]) / dz - (pC[b_yz] * pC[b_yzz]) / (24 * dz) -
          (pC[b_yy] * pC[b_yz]) / (12 * dz) - (pC[b_y] * pC[b_yz]) / (4 * dz) - (pC[b_xy] * pC[b_yz]) / (8 * dz) -
          (pC[b_x] * pC[b_yz]) / (4 * dz) - (pC[b_0] * pC[b_yz]) / (2 * dz) - (pC[b_yy] * pC[b_yyz]) / (36 * dz) -
          (pC[b_y] * pC[b_yyz]) / (12 * dz) - (pC[b_xy] * pC[b_yyz]) / (24 * dz) - (pC[b_x] * pC[b_yyz]) / (12 * dz) -
          (pC[b_0] * pC[b_yyz]) / (6 * dz) - (pC[b_xz] * pC[b_yy]) / (12 * dz) - (pC[b_xyz] * pC[b_yy]) / (24 * dz) -
          (pC[b_xz] * pC[b_y]) / (4 * dz) - (pC[b_xyz] * pC[b_y]) / (8 * dz) - (pC[b_xy] * pC[b_xz]) / (8 * dz) -
          (pC[b_x] * pC[b_xz]) / (4 * dz) - (pC[b_0] * pC[b_xz]) / (2 * dz) - (pC[b_xy] * pC[b_xyz]) / (16 * dz) -
          (pC[b_x] * pC[b_xyz]) / (8 * dz) - (pC[b_0] * pC[b_xyz]) / (4 * dz) - (pC[a_z] * pC[a_zz]) / (6 * dz) -
          (pC[a_xz] * pC[a_zz]) / (12 * dz) - (pC[a_y] * pC[a_z]) / (2 * dz) - (pC[a_xzz] * pC[a_z]) / (12 * dz) -
          (pC[a_xy] * pC[a_z]) / (4 * dz) - (pC[a_xx] * pC[a_z]) / (6 * dz) - (pC[a_x] * pC[a_z]) / (2 * dz) -
          (pC[a_0] * pC[a_z]) / dz - (pC[a_y] * pC[a_yz]) / (4 * dz) - (pC[a_xy] * pC[a_yz]) / (8 * dz) -
          (pC[a_xx] * pC[a_yz]) / (12 * dz) - (pC[a_x] * pC[a_yz]) / (4 * dz) - (pC[a_0] * pC[a_yz]) / (2 * dz) -
          (pC[a_xz] * pC[a_y]) / (4 * dz) - (pC[a_xyz] * pC[a_y]) / (8 * dz) - (pC[a_xxz] * pC[a_y]) / (12 * dz) -
          (pC[a_xz] * pC[a_xzz]) / (24 * dz) - (pC[a_xy] * pC[a_xz]) / (8 * dz) - (pC[a_xx] * pC[a_xz]) / (12 * dz) -
          (pC[a_x] * pC[a_xz]) / (4 * dz) - (pC[a_0] * pC[a_xz]) / (2 * dz) - (pC[a_xy] * pC[a_xyz]) / (16 * dz) -
          (pC[a_xx] * pC[a_xyz]) / (24 * dz) - (pC[a_x] * pC[a_xyz]) / (8 * dz) - (pC[a_0] * pC[a_xyz]) / (4 * dz) -
          (pC[a_xxz] * pC[a_xy]) / (24 * dz) - (pC[a_xx] * pC[a_xxz]) / (36 * dz) - (pC[a_x] * pC[a_xxz]) / (12 * dz) -
          (pC[a_0] * pC[a_xxz]) / (6 * dz) + (pC[b_z] * pC[c_yz]) / (12 * dy) + (pC[b_yz] * pC[c_yz]) / (24 * dy) +
          (pC[b_z] * pC[c_yyz]) / (12 * dy) + (pC[b_yz] * pC[c_yyz]) / (24 * dy) + (pC[b_yy] * pC[c_yy]) / (6 * dy) +
          (pC[b_y] * pC[c_yy]) / (2 * dy) + (pC[b_xy] * pC[c_yy]) / (4 * dy) + (pC[b_x] * pC[c_yy]) / (2 * dy) +
          (pC[b_0] * pC[c_yy]) / dy + (pC[b_yy] * pC[c_y]) / (6 * dy) + (pC[b_y] * pC[c_y]) / (2 * dy) +
          (pC[b_xy] * pC[c_y]) / (4 * dy) + (pC[b_x] * pC[c_y]) / (2 * dy) + (pC[b_0] * pC[c_y]) / dy +
          (pC[b_z] * pC[c_xyz]) / (24 * dy) + (pC[b_yz] * pC[c_xyz]) / (48 * dy) + (pC[b_yy] * pC[c_xy]) / (12 * dy) +
          (pC[b_y] * pC[c_xy]) / (4 * dy) + (pC[b_xy] * pC[c_xy]) / (8 * dy) + (pC[b_x] * pC[c_xy]) / (4 * dy) +
          (pC[b_0] * pC[c_xy]) / (2 * dy) + (pC[a_z] * pC[c_xz]) / (12 * dx) + (pC[a_xz] * pC[c_xz]) / (24 * dx) +
          (pC[a_z] * pC[c_xyz]) / (24 * dx) + (pC[a_xz] * pC[c_xyz]) / (48 * dx) + (pC[a_y] * pC[c_xy]) / (4 * dx) +
          (pC[a_xy] * pC[c_xy]) / (8 * dx) + (pC[a_xx] * pC[c_xy]) / (12 * dx) + (pC[a_x] * pC[c_xy]) / (4 * dx) +
          (pC[a_0] * pC[c_xy]) / (2 * dx) + (pC[a_z] * pC[c_xxz]) / (12 * dx) + (pC[a_xz] * pC[c_xxz]) / (24 * dx) +
          (pC[a_y] * pC[c_xx]) / (2 * dx) + (pC[a_xy] * pC[c_xx]) / (4 * dx) + (pC[a_xx] * pC[c_xx]) / (6 * dx) +
          (pC[a_x] * pC[c_xx]) / (2 * dx) + (pC[a_0] * pC[c_xx]) / dx + (pC[a_y] * pC[c_x]) / (2 * dx) +
          (pC[a_xy] * pC[c_x]) / (4 * dx) + (pC[a_xx] * pC[c_x]) / (6 * dx) + (pC[a_x] * pC[c_x]) / (2 * dx) +
          (pC[a_0] * pC[c_x]) / dx;
}

template <typename REAL>
inline REAL JXB(fsgrids::ehall term, const std::array<REAL, Rec::N_REC_COEFFICIENTS>& pC, Real BGBX, Real BGBY,
                Real BGBZ, const std::array<Real, 3>& gridSpacing) {
   switch (term) {
   case fsgrids::EXHALL_000_100: {
      return JXBX_000_100(pC, BGBY, BGBZ, gridSpacing);
   }
   case fsgrids::EXHALL_010_110: {
      return JXBX_010_110(pC, BGBY, BGBZ, gridSpacing);
   }
   case fsgrids::EXHALL_001_101: {
      return JXBX_001_101(pC, BGBY, BGBZ, gridSpacing);
   }
   case fsgrids::EXHALL_011_111: {
      return JXBX_011_111(pC, BGBY, BGBZ, gridSpacing);
   }
   case fsgrids::EYHALL_000_010: {
      return JXBY_000_010(pC, BGBX, BGBZ, gridSpacing);
   }
   case fsgrids::EYHALL_100_110: {
      return JXBY_100_110(pC, BGBX, BGBZ, gridSpacing);
   }
   case fsgrids::EYHALL_001_011: {
      return JXBY_001_011(pC, BGBX, BGBZ, gridSpacing);
   }
   case fsgrids::EYHALL_101_111: {
      return JXBY_101_111(pC, BGBX, BGBZ, gridSpacing);
   }
   case fsgrids::EZHALL_000_001: {
      return JXBZ_000_001(pC, BGBX, BGBY, gridSpacing);
   }
   case fsgrids::EZHALL_100_101: {
      return JXBZ_100_101(pC, BGBX, BGBY, gridSpacing);
   }
   case fsgrids::EZHALL_010_011: {
      return JXBZ_010_011(pC, BGBX, BGBY, gridSpacing);
   }
   case fsgrids::EZHALL_110_111: {
      return JXBZ_110_111(pC, BGBX, BGBY, gridSpacing);
   }
   case fsgrids::N_EHALL: {
      break;
   }
   }
   return 0.0;
}

/*! \brief Low-level function computing the Hall term numerator x components.
 *
 * Calls the lower-level inline templates and scales the components properly.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param perturbedCoefficients Reconstruction coefficients
 * \param i,j,k fsGrid cell coordinates for the current cell
 *
 * \sa calculateHallTerm JXBX_000_100 JXBX_001_101 JXBX_010_110 JXBX_011_111
 *
 */
void calculateEdgeHallTermComponents(std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perbs,
                                     std::span<std::array<Real, fsgrids::ehall::N_EHALL>> ehalls,
                                     std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments,
                                     std::span<const std::array<Real, fsgrids::dperb::N_DPERB>> dperbs,
                                     std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgbs,
                                     const std::array<Real, 3>& gridSpacing,
                                     const std::array<Real, Rec::N_REC_COEFFICIENTS>& perturbedCoefficients,
                                     const fsgrid::FsStencil& stencil) {
   const auto center = stencil.center();
   const auto& bgb = bgbs[center];
   const auto& perb = perbs[center];
   const auto& dperb = dperbs[center];
   const auto& moment = moments[center];
   auto& ehall = ehalls[center];

   const Real bgbx = bgb[fsgrids::bgbfield::BGBX];
   const Real bgby = bgb[fsgrids::bgbfield::BGBY];
   const Real bgbz = bgb[fsgrids::bgbfield::BGBZ];

   auto computeHallRhoq = [&moments, &moment](const std::array<size_t, 4>& indices) {
      const auto min = Parameters::hallMinimumRhoq;
      const auto max = std::numeric_limits<Real>::max();

      return std::clamp(
          Parameters::ohmHallTerm == 1
              ? moment[fsgrids::moments::RHOQ]
              : FOURTH * (moments[indices[0]][fsgrids::moments::RHOQ] + moments[indices[1]][fsgrids::moments::RHOQ] +
                          moments[indices[2]][fsgrids::moments::RHOQ] + moments[indices[3]][fsgrids::moments::RHOQ]),
          min, max);
   };

   switch (Parameters::ohmHallTerm) {
   case 0:
      cerr << __FILE__ << __LINE__ << "You shouldn't be in a Hall term function if Parameters::ohmHallTerm == 0."
           << endl;
      break;

   case 1: {
      const Real Bx = perb[fsgrids::bfield::PERBX] + bgbx;
      const Real By = perb[fsgrids::bfield::PERBY] + bgby;
      const Real Bz = perb[fsgrids::bfield::PERBZ] + bgbz;

      const Real invHallRhoqMU0 = 1.0 / (physicalconstants::MU_0 * computeHallRhoq({}));

      const Real ydx = (bgb[fsgrids::bgbfield::dBGBydx] + dperb[fsgrids::dperb::dPERBydx]) / gridSpacing[0];
      const Real zdx = (bgb[fsgrids::bgbfield::dBGBzdx] + dperb[fsgrids::dperb::dPERBzdx]) / gridSpacing[0];
      const Real xdy = (bgb[fsgrids::bgbfield::dBGBxdy] + dperb[fsgrids::dperb::dPERBxdy]) / gridSpacing[1];
      const Real zdy = (bgb[fsgrids::bgbfield::dBGBzdy] + dperb[fsgrids::dperb::dPERBzdy]) / gridSpacing[1];
      const Real xdz = (bgb[fsgrids::bgbfield::dBGBxdz] + dperb[fsgrids::dperb::dPERBxdz]) / gridSpacing[2];
      const Real ydz = (bgb[fsgrids::bgbfield::dBGBydz] + dperb[fsgrids::dperb::dPERBydz]) / gridSpacing[2];

      const Real EXHall = (Bz * (xdz - zdx) - By * (ydx - xdy)) * invHallRhoqMU0;
      ehall[fsgrids::ehall::EXHALL_000_100] = EXHall;
      ehall[fsgrids::ehall::EXHALL_010_110] = EXHall;
      ehall[fsgrids::ehall::EXHALL_001_101] = EXHall;
      ehall[fsgrids::ehall::EXHALL_011_111] = EXHall;

      const Real EYHall = (Bx * (ydx - xdy) - Bz * (zdy - ydz)) * invHallRhoqMU0;
      ehall[fsgrids::ehall::EYHALL_000_010] = EYHall;
      ehall[fsgrids::ehall::EYHALL_100_110] = EYHall;
      ehall[fsgrids::ehall::EYHALL_101_111] = EYHall;
      ehall[fsgrids::ehall::EYHALL_001_011] = EYHall;

      const Real EZHall = (By * (zdy - ydz) - Bx * (xdz - zdx)) * invHallRhoqMU0;
      ehall[fsgrids::ehall::EZHALL_000_001] = EZHall;
      ehall[fsgrids::ehall::EZHALL_100_101] = EZHall;
      ehall[fsgrids::ehall::EZHALL_110_111] = EZHall;
      ehall[fsgrids::ehall::EZHALL_010_011] = EZHall;

      break;
   }
   case 2: {
      auto computeEHall = [&perturbedCoefficients, &bgbx, &bgby, &bgbz, &gridSpacing](fsgrids::ehall term,
                                                                                      Real hallRhoq) {
         return JXB(term, perturbedCoefficients, bgbx, bgby, bgbz, gridSpacing) / (physicalconstants::MU_0 * hallRhoq);
      };

      const auto down = stencil.down();
      const auto up = stencil.up();
      const auto left = stencil.left();
      const auto right = stencil.right();
      const auto far = stencil.far();
      const auto near = stencil.near();
      // clang-format off
      const std::array<std::array<size_t, 4>, 12> indices = {
          std::array{
              center,
              down,
              far,
              stencil.downfar(),
          },
          std::array{
              center,
              left,
              far,
              stencil.leftfar(),
          },
          std::array{
              center,
              left,
              down,
              stencil.leftdown(),
          },
          std::array{
              center,
              right,
              far,
              stencil.rightfar(),
          },
          std::array{
              center,
              right,
              down,
              stencil.rightdown(),
          },
          std::array{
              center,
              up,
              far,
              stencil.upfar(),
          },
          std::array{
              center,
              left,
              up,
              stencil.leftup(),
          },
          std::array{
              center,
              right,
              up,
              stencil.rightup(),
          },
          std::array{
              center,
              down,
              near,
              stencil.downnear(),
          },
          std::array{
              center,
              left,
              near,
              stencil.leftnear(),
          },
          std::array{
              center,
              right,
              near,
              stencil.rightnear(),
          },
          std::array{
              center,
              up,
              near,
              stencil.upnear(),
          },
      };
      // clang-format on

      const std::array<fsgrids::ehall, 12> terms = {
          fsgrids::ehall::EXHALL_000_100, fsgrids::ehall::EYHALL_000_010, fsgrids::ehall::EZHALL_000_001,
          fsgrids::ehall::EYHALL_100_110, fsgrids::ehall::EZHALL_100_101, fsgrids::ehall::EXHALL_010_110,
          fsgrids::ehall::EZHALL_010_011, fsgrids::ehall::EZHALL_110_111, fsgrids::ehall::EXHALL_001_101,
          fsgrids::ehall::EYHALL_001_011, fsgrids::ehall::EYHALL_101_111, fsgrids::ehall::EXHALL_011_111,
      };

      for (size_t index = 0; index < terms.size(); index++) {
         ehall[terms[index]] = computeEHall(terms[index], computeHallRhoq(indices[index]));
      }

      break;
   }

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
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary condition functions.
 * \param i,j,k fsGrid cell coordinates for the current cell
 *
 * \sa calculateHallTermSimple calculateEdgeHallTermComponents
 */
void calculateHallTerm(std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                       std::span<std::array<Real, fsgrids::ehall::N_EHALL>> ehall,
                       std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments,
                       std::span<const std::array<Real, fsgrids::dperb::N_DPERB>> dperb,
                       std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                       std::span<const fsgrids::technical> technical, const fsgrid::FsStencil& stencil,
                       SysBoundary& sysBoundaries, const std::array<Real, 3>& gridSpacing) {
#ifdef DEBUG_FSOLVER
   if (!stencil.cellExists(0, 0, 0)) {
      cerr << "Out-of-bounds access in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
#endif

   const auto& tech = technical[stencil.center()];
   cuint cellSysBoundaryFlag = tech.sysBoundaryFlag;
   cuint cellSysBoundaryLayer = tech.sysBoundaryLayer;

   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
       cellSysBoundaryFlag == sysboundarytype::OUTER_BOUNDARY_PADDING) {
      return;
   }

   // Reconstruction order of the fields after Balsara 2009, 2 used for general B, 3 used
   // here for 2nd-order Hall term
   const std::array<Real, Rec::N_REC_COEFFICIENTS> perturbedCoefficients =
       reconstructionCoefficients(perb, dperb, stencil, 3);

   if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && (cellSysBoundaryLayer != 1)) {
      auto* const sb = sysBoundaries.getSysBoundary(cellSysBoundaryFlag);
      sb->fieldSolverBoundaryCondHallElectricField(ehall, stencil, 0);
      sb->fieldSolverBoundaryCondHallElectricField(ehall, stencil, 1);
      sb->fieldSolverBoundaryCondHallElectricField(ehall, stencil, 2);
   } else {
      calculateEdgeHallTermComponents(perb, ehall, moments, dperb, bgb, gridSpacing, perturbedCoefficients, stencil);
   }
}

/*! \brief High-level function computing the Hall term.
 *
 * Performs the communication before and after the computation as well as the computation of all Hall term numerator
 * components.
 *
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param perBDt2Grid fsGrid holding the perturbed B quantities at runge-kutta half step
 * \param EHallGrid fsGrid holding the Hall contributions to the electric field
 * \param momentsGrid fsGrid holding the moment quantities
 * \param momentsDt2Grid fsGrid holding the moment quantities at runge-kutta half step
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param BgBGrid fsGrid holding the background B quantities
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary condition functions.
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 *
 * \sa calculateHallTerm
 */
void calculateHallTermSimple(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsDt2Grid,
    fsgrid::FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, SysBoundary& sysBoundaries, int32_t RKCase,
    const bool communicateMomentsDerivatives) {

   std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> perb = perBGrid.getData();
   std::span<std::array<Real, fsgrids::ehall::N_EHALL>> ehall = EHallGrid.getData();
   std::span<const std::array<Real, fsgrids::moments::N_MOMENTS>> moments = momentsGrid.getData();
   std::span<const std::array<Real, fsgrids::dperb::N_DPERB>> dperb = dPerBGrid.getData();
   std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb = BgBGrid.getData();
   std::span<const fsgrids::technical> technical = technicalGrid.getData();

   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      perb = perBDt2Grid.getData();
      moments = momentsDt2Grid.getData();
   }

   const auto& gridSpacing = technicalGrid.getGridSpacing();
   const auto& localSize = technicalGrid.getLocalSize();
   const size_t N_cells = localSize[0] * localSize[1] * localSize[2];

   phiprof::Timer hallTimer{"Calculate Hall term"};

   phiprof::Timer mpiTimer{"EHall ghost updates MPI", {"MPI"}};
   int computeTimerId{phiprof::initializeTimer("EHall compute cells")};
   dPerBGrid.updateGhostCells();
   if (P::ohmGradPeTerm == 0 && communicateMomentsDerivatives) {
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         dMomentsGrid.updateGhostCells();
      } else {
         dMomentsDt2Grid.updateGhostCells();
      }
   }
   mpiTimer.stop();

#pragma omp parallel
   {
      phiprof::Timer computeTimer{computeTimerId};
#pragma omp for collapse(2)
      for (auto k = 0; k < localSize[2]; k++) {
         for (auto j = 0; j < localSize[1]; j++) {
            for (auto i = 0; i < localSize[0]; i++) {
               const auto stencil = technicalGrid.makeStencil(i, j, k);
               calculateHallTerm(perb, ehall, moments, dperb, bgb, technical, stencil, sysBoundaries, gridSpacing);
            }
         }
      }
      computeTimer.stop(N_cells, "Spatial Cells");
   }

   hallTimer.stop(N_cells, "Spatial Cells");
}
