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

#include <stdint.h>
#include <limits>

//set floating point precision for storing the distribution function here. Default is single precision, use -DDPF to set double precision
#ifdef DPF
typedef double Realf;
#else
typedef float Realf;
#endif

//set general floating point precision here. Default is single precision, use -DDP to set double precision
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

template<typename T> T convert(const T& number) {return number;}

/** Definition of simulation geometry, used to speed up 
 * computations in cases where all velocity coordinates are not used.*/
namespace geometry {
   enum Setup {
      XY4D,            /**< Simulation is 2D, only x,y,vx,vy are used.*/
      XZ4D,            /**< Simulation is 2D, only x,z,vx,vz are used.*/
      XY5D,            /**< Simulation is 5D, only x,y,vx,vy,vz are used.*/
      XZ5D,            /**< Simulation is 5D, only x,z,vx,vy,vz are used.*/
      XYZ6D            /**< Simulation is 6D (default).*/
   };
}

namespace vmesh {
   #ifndef VAMR
   typedef uint32_t GlobalID;              /**< Datatype used for velocity block global IDs.*/
   typedef uint32_t LocalID;               /**< Datatype used for velocity block local IDs.*/
   #else
   typedef uint32_t GlobalID;
   typedef uint32_t LocalID;
   #endif

   /** Global ID of a non-existing or otherwise erroneous velocity block.*/
   static const GlobalID INVALID_GLOBALID = std::numeric_limits<GlobalID>::max();

   /** Local ID of a non-existing or otherwise erroneous velocity block.*/
   static const LocalID INVALID_LOCALID  = std::numeric_limits<LocalID>::max();

   /** Block index of a non-existing or erroneous velocity block.*/
   static const LocalID INVALID_VEL_BLOCK_INDEX = INVALID_LOCALID;
}

/** Definition of a function that takes in a velocity block with neighbor data, 
 * and returns a number that is used to decide whether or not the block should 
 * be refined or coarsened.*/
typedef Realf (*VAmrVelRefinement)(const Realf* velBlock);

// neighborhoods, these are initialized in grid.cpp:initializeGrid

//#define FIELD_SOLVER_NEIGHBORHOOD_ID 1 // Now fields on FSgrid
#define VLASOV_SOLVER_NEIGHBORHOOD_ID 2   //up to second(PPM) neighbor in each face direction
#define VLASOV_SOLVER_X_NEIGHBORHOOD_ID 3 //up to second(PPM) neighbor in x face directions
#define VLASOV_SOLVER_Y_NEIGHBORHOOD_ID 4 //up to second(PPM) neighbor in y face directions
#define VLASOV_SOLVER_Z_NEIGHBORHOOD_ID 5 //up to second(PPM) neighbor in z face directions
#define VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID 6 //nearest neighbor in X face direction, f() can propagate to local cells in X dir, and are target for local cells
#define VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID 7 //nearest neighbor in Y face direction, f() can propagate to local cells in Y dir, and are target for local cells
#define VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID 8 //nearest neighbor in Z face direction, f() can propagate to local cells in Z dir, and are target for local cells
#define SYSBOUNDARIES_NEIGHBORHOOD_ID 9 // When classifying sysboundaries, all 26 nearest neighbors are included,
#define SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID 10 //Up to second nearest neighbors in all directions (also diagonals)
#define NEAREST_NEIGHBORHOOD_ID 11  //nearest neighbors
#define FULL_NEIGHBORHOOD_ID 12      //Up to second nearest neighbors in all directions (also diagonals) + vlasov solver neighborhood
#define DIST_FUNC_NEIGHBORHOOD_ID 13 //nearest neighbors in all directions (also diagonals) + vlasov solver neighborhood
#define SHIFT_P_X_NEIGHBORHOOD_ID 14 //Shift in +x direction
#define SHIFT_P_Y_NEIGHBORHOOD_ID 15 //Shift in +y direction
#define SHIFT_P_Z_NEIGHBORHOOD_ID 16 //Shift in +z direction
#define SHIFT_M_X_NEIGHBORHOOD_ID 17 //Shift in -x direction
#define SHIFT_M_Y_NEIGHBORHOOD_ID 18 //Shift in -y direction
#define SHIFT_M_Z_NEIGHBORHOOD_ID 19 //Shift in -z direction
#define VLASOV_SOLVER_X_GHOST_NEIGHBORHOOD_ID 20 //up to third(PPM+ghost) neighbor in x face directions
#define VLASOV_SOLVER_Y_GHOST_NEIGHBORHOOD_ID 21 //up to third(PPM+ghost) neighbor in y face directions
#define VLASOV_SOLVER_Z_GHOST_NEIGHBORHOOD_ID 22 //up to third(PPM+ghost) neighbor in z face directions
#define VLASOV_SOLVER_GHOST_NEIGHBORHOOD_ID 23 //all required neighbors for ghost translation
#define VLASOV_SOLVER_GHOST_REQNEIGH_NEIGHBORHOOD_ID 24 //all ghost translation neighbors which require own neighbor information

//fieldsolver stencil.
#define FS_STENCIL_WIDTH 2

//Vlasov propagator stencils in ordinary space, velocity space may be
//higher. Assume H4 (or H5) for PPM, H6 for PQM
#ifdef TRANS_SEMILAG_PLM
   #define  VLASOV_STENCIL_WIDTH 1
#endif
#ifdef TRANS_SEMILAG_PPM
   #define  VLASOV_STENCIL_WIDTH 2
#endif
#ifdef TRANS_SEMILAG_PQM
   #define  VLASOV_STENCIL_WIDTH 3
#endif

// Max number of face neighbors per dimension with AMR
#define MAX_NEIGHBORS_PER_DIM 8
#define MAX_FACE_NEIGHBORS_PER_DIM 4

#endif
