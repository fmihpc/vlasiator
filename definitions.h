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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <limits>
#include <stdint.h>

// set floating point precision for storing the distribution function here. Default is single precision, use -DDPF to set double precision
#ifdef DPF
typedef double Realf;
#else
typedef float Realf;
#endif

// set general floating point precision here. Default is single precision, use -DDP to set double precision
#ifdef DP
typedef double Real;
typedef const double creal;
#else
typedef float Real;
typedef const float creal;
#endif

typedef const int cint;
typedef unsigned char uchar;
typedef const unsigned char cuchar;

typedef uint32_t uint;
typedef const uint32_t cuint;

typedef cuint csize;

typedef uint64_t CellID;

template <typename T> T convert(const T& number) { return number; }

namespace vmesh {
   typedef uint32_t GlobalID; /**< Datatype used for velocity block global IDs.*/
   typedef uint32_t LocalID;  /**< Datatype used for velocity block local IDs.*/

   /** Global ID of a non-existing or otherwise erroneous velocity block.*/
   static const GlobalID INVALID_GLOBALID = std::numeric_limits<GlobalID>::max();

   /** Local ID of a non-existing or otherwise erroneous velocity block.*/
   static const LocalID INVALID_LOCALID = std::numeric_limits<LocalID>::max();

   /** Block index of a non-existing or erroneous velocity block.*/
   static const LocalID INVALID_VEL_BLOCK_INDEX = INVALID_LOCALID;
} // namespace vmesh

// fieldsolver stencil.
#define FS_STENCIL_WIDTH 2

// Vlasov propagator stencils in ordinary space, velocity space may be
// higher. Assume H4 (or H5) for PPM, H6 for PQM
#ifdef TRANS_SEMILAG_PLM
#define VLASOV_STENCIL_WIDTH 1
#endif
#ifdef TRANS_SEMILAG_PPM
   #define  VLASOV_STENCIL_WIDTH 2
#endif
#ifdef TRANS_SEMILAG_PQM
#define VLASOV_STENCIL_WIDTH 3
#endif

// Max number of face neighbors per dimension with AMR
#define MAX_NEIGHBORS_PER_DIM 8
#define MAX_FACE_NEIGHBORS_PER_DIM 4

#endif
