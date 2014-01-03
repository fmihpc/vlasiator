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


   






   

#endif

