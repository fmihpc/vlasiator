/*
This file is part of Vlasiator.

Copyright 2012 Finnish Meteorological Institute

*/

#ifndef CPU_TRANS_SEMILAG_H
#define CPU_TRANS_SEMILAG_H

#include "algorithm"
#include "cmath"
#include "utility"

/*TODO - replace with standard library c++11 functions as soon as possible*/
#include "boost/array.hpp"
#include "boost/unordered_map.hpp"

#include "common.h"
#include "spatial_cell.hpp"

#include "vlasovsolver/cpu_trans_intersections.hpp"

using namespace std;
using namespace spatial_cell;
using namespace Eigen;





void cpu_translate_cell(const dccrg::Dccrg<SpatialCell>& mpiGrid,const CellID cellID,const Real dt) {


  /*compute transform, forward in time and backward in time*/

 }


   





CellID getNeighbourID(
   const dccrg::Dccrg<SpatialCell>& mpiGrid,
   const CellID& cellID,
   const uchar& i,
   const uchar& j,
   const uchar& k ) {
  const std::vector<CellID> neighbors = mpiGrid.get_neighbors_of(cellID, int(i) - 2, int(j) - 2, int(k) - 2);
  if (neighbors.size() == 0) {
    std::cerr << __FILE__ << ":" << __LINE__
	      << " No neighbor for cell " << cellID
	      << " at offsets " << int(i) - 2 << ", " << int(j) - 2 << ", " << int(k) - 2
                 << std::endl;
    abort();
  }
  // TODO support spatial refinement
  if( neighbors[0] == INVALID_CELLID  ||
      mpiGrid[neighbors[0]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
      (mpiGrid[neighbors[0]]->sysBoundaryLayer != 1  &&
       mpiGrid[neighbors[0]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY )
      ) {
     return INVALID_CELLID;
  } 
  else {
    return neighbors[0];
  }
}


   

#endif

