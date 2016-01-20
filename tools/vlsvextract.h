/* This file is part of Vlasiator.
 * 
 * File:   vlsvextract.h
 * Author: sandroos
 *
 * Created on April 14, 2015, 2:08 PM
 * 
 * Copyright 2015 Finnish Meteorological Institute
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

