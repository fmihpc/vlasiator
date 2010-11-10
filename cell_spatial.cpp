#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>

#include "common.h"
#include "cell_spatial.h"
#include "parameters.h"
#include "cudafuncs.h"

using namespace std;

Real* allocateReal(cuint& ELEMENTS,const bool& PAGELOCKED) {
   #ifdef NOCUDA
   return new Real[ELEMENTS];
   #else
   if (PAGELOCKED == true) {
      Real* ptr;
      if (deviceCreateArray(ptr,ELEMENTS*sizeof(Real)) == false) return NULL;
      return ptr;
   } else {
      return new Real[ELEMENTS];
   }
   #endif
}

uint* allocateUint(cuint& ELEMENTS,const bool& PAGELOCKED) {
   #ifdef NOCUDA
   return new uint[ELEMENTS];
   #else
   if (PAGELOCKED == true) {
      uint* ptr;
      if (deviceCreateArray(ptr,ELEMENTS*sizeof(uint)) == false) return NULL;
      return ptr;
   } else {
      return new uint[ELEMENTS];
   }
   #endif
}

void freeReal(Real*& ptr,const bool& PAGELOCKED) {
   #ifdef NOCUDA
   delete ptr;
   #else
   if (PAGELOCKED == true) {
      deviceDeleteArray(ptr);
   } else {
      delete ptr;
   }
   #endif
   ptr = NULL;
}

void freeUint(uint*& ptr,const bool& PAGELOCKED) {
   #ifdef NOCUDA
   delete ptr;
   #else
   if (PAGELOCKED == true) {
      deviceDeleteArray(ptr);
   } else {
      delete ptr;
   }
   #endif
   ptr = NULL;
} 

SpatialCell::SpatialCell() {
   typedef Parameters P;
   pageLocked = false; // Should determine this from somewhere
   N_blocks = P::vxblocks_ini*P::vyblocks_ini*P::vzblocks_ini;
   
   cpuIndex = numeric_limits<uint>::max();
   gpuIndex = numeric_limits<uint>::max();
   
   cellType = Cell::UNINITIALIZED;
   cellIndex = numeric_limits<uint>::max();
   velBlockIndex = numeric_limits<uint>::max();
   i_ind = numeric_limits<uint>::max();
   j_ind = numeric_limits<uint>::max();
   k_ind = numeric_limits<uint>::max();
   
   // Allocate memory for data in CPU memory:
   allocateMemory();

   // Pointers to corresponding data in GPU memory. Leave as NULL for now.
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
   //std::cerr << "SpatialCell copy constructor called" << std::endl;
   
   // Determine the amount of data and allocate arrays:
   N_blocks = s.N_blocks;
   pageLocked = s.pageLocked;
   allocateMemory();
   // Copy CPU data:
   copyMemory(s);
   /*
   for (uint i=0; i<SIZE_NBRS_SPA; ++i) cpu_nbrsSpa[i] = s.cpu_nbrsSpa[i];
   for (uint i=0; i<SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   for (uint i=0; i<N_blocks*SIZE_NBRS_VEL; ++i) cpu_nbrsVel[i] = s.cpu_nbrsVel[i];
   for (uint i=0; i<N_blocks*SIZE_BLOCKPARAMS; ++i) cpu_blockParams[i] = s.cpu_blockParams[i];
   for (uint i=0; i<N_blocks*SIZE_VELBLOCK; ++i) cpu_avgs[i] = s.cpu_avgs[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) cpu_fx[i] = s.cpu_fx[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) cpu_fy[i] = s.cpu_fy[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) cpu_fz[i] = s.cpu_fz[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d1x[i] = s.cpu_d1x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d1y[i] = s.cpu_d1y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d1z[i] = s.cpu_d1z[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d2x[i] = s.cpu_d2x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d2y[i] = s.cpu_d2y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d2z[i] = s.cpu_d2z[i];
   */
   copyVariables(s);
   /*
   // Copy the rest of member variables:
   cpuIndex = s.cpuIndex;
   gpuIndex = s.gpuIndex;
   
   cellType = s.cellType;
   cellIndex = s.cellIndex;
   velBlockIndex = s.velBlockIndex;
   i_ind = s.i_ind;
   j_ind = s.j_ind;
   k_ind = s.k_ind;

   // For now, just copy the pointers to GPU memory:
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
    */
}

SpatialCell::~SpatialCell() {
   //cerr << "Cell #" << cpuIndex << " freeing up memory" << endl;
   
   // Free CPU memory:
   freeMemory();
   
   // GPU data, CHECK THIS !!!
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

void SpatialCell::allocateMemory() {
   // Allocate memory for cell data:
   cpu_nbrsSpa    = allocateUint(SIZE_NBRS_SPA,pageLocked);
   cpu_cellParams = allocateReal(SIZE_CELLPARAMS,pageLocked);
   // Allocate memory for velocity grid block data:
   cpu_nbrsVel = allocateUint(N_blocks*SIZE_NBRS_VEL,pageLocked);
   cpu_blockParams = allocateReal(N_blocks*SIZE_BLOCKPARAMS,pageLocked);
   cpu_avgs = allocateReal(N_blocks*SIZE_VELBLOCK,pageLocked);
   cpu_fx   = allocateReal(N_blocks*SIZE_FLUXS,pageLocked);
   cpu_fy   = allocateReal(N_blocks*SIZE_FLUXS,pageLocked);
   cpu_fz   = allocateReal(N_blocks*SIZE_FLUXS,pageLocked);
   cpu_d1x  = allocateReal(N_blocks*SIZE_DERIV,pageLocked);
   cpu_d1y  = allocateReal(N_blocks*SIZE_DERIV,pageLocked);
   cpu_d1z  = allocateReal(N_blocks*SIZE_DERIV,pageLocked);
   cpu_d2x  = allocateReal(N_blocks*SIZE_DERIV,pageLocked);
   cpu_d2y  = allocateReal(N_blocks*SIZE_DERIV,pageLocked);
   cpu_d2z  = allocateReal(N_blocks*SIZE_DERIV,pageLocked);
}

void SpatialCell::freeMemory() {
   freeUint(cpu_nbrsSpa,pageLocked);
   freeReal(cpu_cellParams,pageLocked);

   freeUint(cpu_nbrsVel,pageLocked);
   freeReal(cpu_blockParams,pageLocked);
   freeReal(cpu_avgs,pageLocked);
   freeReal(cpu_fx  ,pageLocked);
   freeReal(cpu_fy  ,pageLocked);
   freeReal(cpu_fz  ,pageLocked);
   freeReal(cpu_d1x ,pageLocked);
   freeReal(cpu_d1y ,pageLocked);
   freeReal(cpu_d1z ,pageLocked);
   freeReal(cpu_d2x ,pageLocked);
   freeReal(cpu_d2y ,pageLocked);
   freeReal(cpu_d2z ,pageLocked);
}

void SpatialCell::copyMemory(const SpatialCell& s) {
   for (uint i=0; i<SIZE_NBRS_SPA; ++i) cpu_nbrsSpa[i] = s.cpu_nbrsSpa[i];
   for (uint i=0; i<SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   for (uint i=0; i<N_blocks*SIZE_NBRS_VEL; ++i) cpu_nbrsVel[i] = s.cpu_nbrsVel[i];
   for (uint i=0; i<N_blocks*SIZE_BLOCKPARAMS; ++i) cpu_blockParams[i] = s.cpu_blockParams[i];
   for (uint i=0; i<N_blocks*SIZE_VELBLOCK; ++i) cpu_avgs[i] = s.cpu_avgs[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) cpu_fx[i] = s.cpu_fx[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) cpu_fy[i] = s.cpu_fy[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) cpu_fz[i] = s.cpu_fz[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d1x[i] = s.cpu_d1x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d1y[i] = s.cpu_d1y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d1z[i] = s.cpu_d1z[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d2x[i] = s.cpu_d2x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d2y[i] = s.cpu_d2y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) cpu_d2z[i] = s.cpu_d2z[i];
}

void SpatialCell::copyVariables(const SpatialCell& s) {
   N_blocks = s.N_blocks;
   pageLocked = s.pageLocked;
   
   cpuIndex = s.cpuIndex;
   gpuIndex = s.gpuIndex;
   
   cellType = s.cellType;
   cellIndex = s.cellIndex;
   velBlockIndex = s.velBlockIndex;
   i_ind = s.i_ind;
   j_ind = s.j_ind;
   k_ind = s.k_ind;
   
   // For now, just copy GPU memory pointers:
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

void SpatialCell::swapMemoryType() {
   // Not implemented yet
}

SpatialCell& SpatialCell::operator=(const SpatialCell& s) {
   //std::cerr << "SpatialCell operator= called" << std::endl;

   // Deallocate old arrays, allocate new ones with the correct size, and copy:
   freeMemory();
   N_blocks = s.N_blocks;
   pageLocked = s.pageLocked;
   allocateMemory();
   copyMemory(s);
   
   // Copy the rest of member variables:
   copyVariables(s);
   
   return *this;
}


