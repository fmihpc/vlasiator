/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <cuda_runtime_api.h>

#include "common.h"
#include "spatial_cell.hpp"
#include "logger.h"

using namespace std;

extern Logger logger;

#ifndef NOCUDA
bool SpatialCell::devSync(const Cell::Array& array,const Cell::Dir direction) {
   // Set pointers to cpu and gpu arrays which are to be synced, 
   // and calculate the amount of bytes to transfer:
   void* cpuArray;
   void* gpuArray;
   size_t bytes;
   switch (array) {
    case Cell::Blocks:
      cpuArray = cpu_avgs;
      gpuArray = gpu_avgs;
      bytes = N_blocks*SIZE_VELBLOCK*sizeof(real);
      break;
    case Cell::BlockParams:
      cpuArray = cpu_blockParams;
      gpuArray = gpu_blockParams;
      bytes = N_blocks*SIZE_BLOCKPARAMS*sizeof(real);
      break;
    case Cell::NbrsVel:
      cpuArray = cpu_nbrsVel;
      gpuArray = gpu_nbrsVel;
      bytes = N_blocks*SIZE_NBRS_VEL*sizeof(uint);
      break;
    case Cell::CellParams:
      cpuArray = cpu_cellParams;
      gpuArray = gpu_cellParams;
      bytes = SIZE_CELLPARAMS*sizeof(real);
      break;
    case Cell::Fx:
      cpuArray = cpu_fx;
      gpuArray = gpu_fx;
      bytes = N_blocks*SIZE_FLUXS*sizeof(real);
      break;
    case Cell::Fy:
      cpuArray = cpu_fy;
      gpuArray = gpu_fy;
      bytes = N_blocks*SIZE_FLUXS*sizeof(real);
      break;
    case Cell::Fz:
      cpuArray = cpu_fz;
      gpuArray = gpu_fz;
      bytes = N_blocks*SIZE_FLUXS*sizeof(real);
      break;
    case Cell::D1x:
      cpuArray = cpu_d1x;
      gpuArray = gpu_d1x;
      bytes = N_blocks*SIZE_DERIV*sizeof(real);
      break;
    case Cell::D1y:
      cpuArray = cpu_d1y;
      gpuArray = gpu_d1y;
      bytes = N_blocks*SIZE_DERIV*sizeof(real);
      break;
    case Cell::D1z:
      cpuArray = cpu_d1z;
      gpuArray = gpu_d1z;
      bytes = N_blocks*SIZE_DERIV*sizeof(real);
      break;
   }
   // Set destination and source arrays:
   void* src;
   void* dst;
   cudaMemcpyKind copyDirection;
   switch (direction) {
    case Cell::CpuToDev:
      src = cpuArray;
      dst = gpuArray;
      copyDirection = cudaMemcpyHostToDevice;
      break;
    case Cell::DevToCpu:
      src = gpuArray;
      dst = cpuArray;
      copyDirection = cudaMemcpyDeviceToHost;
      break;
   }
   // Perform copy and return result:
   //cerr << "\t Copying from " << src << " to " << dst << endl;
   cudaError_t result = cudaMemcpy(dst,src,bytes,copyDirection);   
   if (result == cudaSuccess) return true;

   logger << "(CELLSYNC): Error occurred while copying data, error = '" << cudaGetErrorString(result) << "'" << endl;
   return false;
}
#endif



