#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>

#include "cell_spatial.h"

using namespace std;

SpatialCell::SpatialCell() {
   gpuIndex = numeric_limits<uint>::max();
   
   cellType = Cell::UNINITIALIZED;
   cellIndex = numeric_limits<uint>::max();
   velBlockIndex = numeric_limits<uint>::max();
   i_ind = numeric_limits<uint>::max();
   j_ind = numeric_limits<uint>::max();
   k_ind = numeric_limits<uint>::max();

   cpu_nbrsSpa = NULL;
   cpu_cellParams = NULL;
   
   cpu_nbrsVel = NULL;
   cpu_avgs = NULL;
   cpu_blockParams = NULL;

   gpu_cellParams = NULL;
   
   gpu_nbrsVel = NULL;
   gpu_avgs = NULL;
   gpu_blockParams = NULL;
}

SpatialCell::SpatialCell(const SpatialCell& s) {
   gpuIndex = s.gpuIndex;
   
   cellType = s.cellType;
   cellIndex = s.cellIndex;
   velBlockIndex = s.velBlockIndex;
   i_ind = s.i_ind;
   j_ind = s.j_ind;
   k_ind = s.k_ind;
   
   cpu_nbrsSpa = s.cpu_nbrsSpa;
   cpu_cellParams = s.cpu_cellParams;
   
   cpu_nbrsVel = s.cpu_nbrsVel;
   cpu_avgs = s.cpu_avgs;
   cpu_blockParams = s.cpu_blockParams;
   
   gpu_cellParams = s.gpu_cellParams;
   
   gpu_nbrsVel = s.gpu_nbrsVel;
   gpu_avgs = s.gpu_avgs;
   gpu_blockParams = s.gpu_blockParams;
}

SpatialCell::~SpatialCell() {
   cpu_nbrsSpa = NULL;
   cpu_cellParams = NULL;
   
   cpu_nbrsVel = NULL;
   cpu_avgs = NULL;
   cpu_blockParams = NULL;
   
   gpu_cellParams = NULL;
   
   gpu_nbrsVel = NULL;
   gpu_avgs = NULL;
   gpu_blockParams = NULL;
}


