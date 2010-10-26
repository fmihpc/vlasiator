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
   cpu_fx = NULL;
   cpu_fy = NULL;
   cpu_fz = NULL;
   cpu_d1x = NULL;
   cpu_d1y = NULL;
   cpu_d1z = NULL;
   cpu_d2x = NULL;
   cpu_d2y = NULL;
   cpu_d2z = NULL;
   
   gpu_cellParams = NULL;
   
   gpu_nbrsVel = NULL;
   gpu_avgs = NULL;
   gpu_blockParams = NULL;
   gpu_fx = NULL;
   gpu_fy = NULL;
   gpu_fz = NULL;
   gpu_d1x = NULL;
   gpu_d1y = NULL;
   gpu_d1z = NULL;
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
   cpu_fx = s.cpu_fx;
   cpu_fy = s.cpu_fy;
   cpu_fz = s.cpu_fz;
   cpu_d1x = s.cpu_d1x;
   cpu_d1y = s.cpu_d1y;
   cpu_d1z = s.cpu_d1z;
   cpu_d2x = s.cpu_d2x;
   cpu_d2y = s.cpu_d2y;
   cpu_d2z = s.cpu_d2z;
   
   gpu_cellParams = s.gpu_cellParams;
   
   gpu_nbrsVel = s.gpu_nbrsVel;
   gpu_avgs = s.gpu_avgs;
   gpu_blockParams = s.gpu_blockParams;
   gpu_fx = s.gpu_fx;
   gpu_fy = s.gpu_fy;
   gpu_fz = s.gpu_fz;
   gpu_d1x = s.gpu_d1x;
   gpu_d1y = s.gpu_d1y;
   gpu_d1z = s.gpu_d1z;
}

SpatialCell::~SpatialCell() {
   cpu_nbrsSpa = NULL;
   cpu_cellParams = NULL;
   
   cpu_nbrsVel = NULL;
   cpu_avgs = NULL;
   cpu_blockParams = NULL;
   cpu_fx = NULL;
   cpu_fy = NULL;
   cpu_fz = NULL;
   cpu_d1x = NULL;
   cpu_d1y = NULL;
   cpu_d1z = NULL;
   cpu_d2x = NULL;
   cpu_d2y = NULL;
   cpu_d2z = NULL;
   
   gpu_cellParams = NULL;
   
   gpu_nbrsVel = NULL;
   gpu_avgs = NULL;
   gpu_blockParams = NULL;
   gpu_fx = NULL;
   gpu_fy = NULL;
   gpu_fz = NULL;
   gpu_d1x = NULL;
   gpu_d1y = NULL;
   gpu_d1z = NULL;
}


