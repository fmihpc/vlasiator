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
 *
 * File:   vlsvextract.h
 * Author: sandroos
 *
 * Created on April 14, 2015, 2:08 PM
 * 
 */

#ifndef VLSVEXTRACT_H
#define	VLSVEXTRACT_H

#include <cstdlib>
#include <array>
#include <vector>

#include "definitions.h"

//A struct for holding info on cell structure (the grid)
struct CellStructure {
   uint64_t cell_bounds[3];     /**< The number of cells in x, y, z direction.*/
   Real cell_length[3];         /**< Length of a cell in x, y, z direction. */
   Real min_coordinates[3];     /**< x_min, y_min, z_min are stored here.*/

   uint64_t vcell_bounds[3];    /**< The number of cells in vx, vy, vz directions.*/
   Real vblock_length[3];       /**< Size (dvx,dvy,dvz) of a velocity block in vx, vy, vz directions.*/
   Real min_vcoordinates[3];    /**< vx_min, vy_min, vz_min are stored here.*/
   uint32_t maxVelRefLevel;     /**< Maximum refinement level of velocity meshes.*/

   int slicedCoords[3];
   Real slicedCoordValues[3];
};

template<typename REAL>
struct NodeCrd {
   static REAL EPS;
   REAL x;
   REAL y;
   REAL z;

   NodeCrd(const REAL& x, const REAL& y, const REAL& z);
   
   bool comp(const NodeCrd<REAL>& n) const;
};

struct NodeComp {
   bool operator()(const NodeCrd<double>& a, const NodeCrd<double>& b) const;
   bool operator()(const NodeCrd<float>& a, const NodeCrd<float>& b) const;
};

//A class for holding user options
class UserOptions {
public:
   bool getCellIdFromLine;
   bool getCellIdFromInput;
   bool getCellIdFromCoordinates;
   bool rotateVectors;
   bool plasmaFrame;
   uint64_t cellId;
   std::vector<uint64_t> cellIdList;
   uint32_t numberOfCoordinatesInALine;
   std::vector<std::string> outputDirectoryPath;
   std::array<Real, 3> coordinates;
   std::array<Real, 3> point1;
   std::array<Real, 3> point2;

   UserOptions() {
      getCellIdFromLine = false;
      getCellIdFromInput = false;
      getCellIdFromCoordinates = false;
      rotateVectors = false;
      plasmaFrame =false;
      cellId = std::numeric_limits<uint64_t>::max();
      numberOfCoordinatesInALine = 0;
   }

   ~UserOptions() {}
};

bool setVelocityMeshVariables(vlsv::Reader& vlsvReader,CellStructure& cellStruct);
bool setVelocityMeshVariables(vlsv::Reader& vlsvReader,CellStructure& cellStruct,
                              const std::string& popName);

template<typename REAL> inline
NodeCrd<REAL>::NodeCrd(const REAL& x,const REAL& y,const REAL& z): x(x),y(y),z(z) { }

template<typename REAL> inline
bool NodeCrd<REAL>::comp(const NodeCrd<REAL>& n) const {
   REAL EPS1, EPS2, EPS;
   EPS1 = 1.0e-6 * fabs(x);
   EPS2 = 1.0e-6 * fabs(n.x);
   if (x == 0.0) EPS1 = 1.0e-7;
   if (n.x == 0.0) EPS2 = 1.0e-7;
   EPS = max(EPS1, EPS2);
   if (fabs(x - n.x) > EPS) return false;

   EPS1 = 1.0e-6 * fabs(y);
   EPS2 = 1.0e-6 * fabs(n.y);
   if (y == 0.0) EPS1 = 1.0e-7;
   if (n.y == 0.0) EPS2 = 1.0e-7;
   EPS = max(EPS1, EPS2);
   if (fabs(y - n.y) > EPS) return false;

   EPS1 = 1.0e-6 * fabs(z);
   EPS2 = 1.0e-6 * fabs(n.z);
   if (z == 0.0) EPS1 = 1.0e-7;
   if (n.z == 0.0) EPS2 = 1.0e-7;
   EPS = max(EPS1, EPS2);
   if (fabs(z - n.z) > EPS) return false;
   return true;
}

#endif	// VLSVEXTRACT_H

