#include <cstdlib>
#include <iostream>
#include <cuda_runtime_api.h>

#include "common.h"
#include "cell_spatial.h"
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



