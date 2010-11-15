#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>

#include "logger.h"
#include "common.h"
#include "cell_spatial.h"
#include "parameters.h"

using namespace std;

extern Grid grid;
extern Logger logger;

SpatialCell::SpatialCell() {
   typedef Parameters P;
   N_blocks = 0;
   
   cpuIndex = numeric_limits<uint>::max();
   cpu_cellParams = new Real[SIZE_CELLPARAMS];
   cpu_nbrsSpa    = new uint[SIZE_NBRS_SPA];
}

SpatialCell::SpatialCell(const SpatialCell& s) {
   // Copy variables related to the spatial cell:
   N_blocks       = s.N_blocks;
   cpu_cellParams = new Real[SIZE_CELLPARAMS];
   cpu_nbrsSpa    = new uint[SIZE_NBRS_SPA];
   for (uint i=0; i<SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   for (uint i=0; i<SIZE_NBRS_SPA; ++i  ) cpu_nbrsSpa[i]    = s.cpu_nbrsSpa[i];
   
   cpuIndex = s.cpuIndex;
   // If the SpatialCell to copy has allocated memory, increase the reference count 
   // and set pointers:
   if (cpuIndex == numeric_limits<uint>::max()) return;
   if (grid.addReference(cpuIndex) == false) {
      logger << "SpatialCell: reference increase failed, aborting." << endl << flush;
      exit(1);
   }
   cpu_nbrsVel     = s.cpu_nbrsVel;
   cpu_blockParams = s.cpu_blockParams;
   cpu_avgs        = s.cpu_avgs;
   cpu_fx          = s.cpu_fx;
   cpu_fy          = s.cpu_fy;
   cpu_fz          = s.cpu_fz;
   cpu_d1x         = s.cpu_d1x;
   cpu_d1y         = s.cpu_d1y;
   cpu_d1z         = s.cpu_d1z;
   cpu_d2x         = s.cpu_d2x;
   cpu_d2y         = s.cpu_d2y;
   cpu_d2z         = s.cpu_d2z;
}

SpatialCell& SpatialCell::operator=(const SpatialCell& s) {
   // Clear previous memory:
   if (cpuIndex != numeric_limits<uint>::max()) {
      if (grid.removeReference(cpuIndex) == false) {
	 cerr << "SpatialCell operator=: Failed to remove reference." << endl;
      }
   }
   // Copy variables related to the spatial cell:
   N_blocks = s.N_blocks;
   for (uint i=0; i<SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   for (uint i=0; i<SIZE_NBRS_SPA; ++i  ) cpu_nbrsSpa[i]    = s.cpu_nbrsSpa[i];
   
   // Copy variables related to the velocity grid:
   cpuIndex = s.cpuIndex;
   if (cpuIndex == numeric_limits<uint>::max()) return *this;
   
   if (grid.addReference(cpuIndex) == false) {
      logger << "SpatialCell: reference increase failed, aborting." << endl << flush;
      logger << "\t operator= got cpuIndex " << s.cpuIndex << " from copied SpatialCell." << endl << flush;
      exit(1);
   }
   cpu_nbrsVel     = s.cpu_nbrsVel;
   cpu_blockParams = s.cpu_blockParams;
   cpu_avgs        = s.cpu_avgs;
   cpu_fx          = s.cpu_fx;
   cpu_fy          = s.cpu_fy;
   cpu_fz          = s.cpu_fz;
   cpu_d1x         = s.cpu_d1x;
   cpu_d1y         = s.cpu_d1y;
   cpu_d1z         = s.cpu_d1z;
   cpu_d2x         = s.cpu_d2x;
   cpu_d2y         = s.cpu_d2y;
   cpu_d2z         = s.cpu_d2z;
   return *this;
}

SpatialCell::~SpatialCell() {
   // Free CPU memory:
   freeMemory();
}

bool SpatialCell::allocateMemory(const uint& BLOCKS) {
   N_blocks = BLOCKS;
   // Attempt to get pointers to available memory from Grid:
   cpuIndex = grid.getFreeMemory(N_blocks);
   if (cpuIndex == numeric_limits<uint>::max()) return false;
   cpu_nbrsVel     = grid.getNbrsVel()     + cpuIndex*SIZE_NBRS_VEL;
   cpu_blockParams = grid.getBlockParams() + cpuIndex*SIZE_BLOCKPARAMS;
   cpu_avgs        = grid.getAvgs()        + cpuIndex*SIZE_VELBLOCK;
   cpu_fx          = grid.getFx()          + cpuIndex*SIZE_FLUXS;
   cpu_fy          = grid.getFy()          + cpuIndex*SIZE_FLUXS;
   cpu_fz          = grid.getFz()          + cpuIndex*SIZE_FLUXS;
   cpu_d1x         = grid.getD1x()         + cpuIndex*SIZE_DERIV;
   cpu_d1y         = grid.getD1y()         + cpuIndex*SIZE_DERIV;
   cpu_d1z         = grid.getD1z()         + cpuIndex*SIZE_DERIV;
   cpu_d2x         = grid.getD2x()         + cpuIndex*SIZE_DERIV;
   cpu_d2y         = grid.getD2y()         + cpuIndex*SIZE_DERIV;
   cpu_d2z         = grid.getD2z()         + cpuIndex*SIZE_DERIV;
   return true;
}

bool SpatialCell::freeMemory() {
   if (cpuIndex != numeric_limits<uint>::max()) {
      if (grid.removeReference(cpuIndex) == false) {
	 logger << "SpatialCell ERROR: Reference removal failed" << endl << flush;
      }
      cpuIndex = numeric_limits<uint>::max();
   }
   delete cpu_cellParams;
   delete cpu_nbrsSpa;
   return true;
}

void SpatialCell::getMemInfo() {
   logger << "cpuIndex = " << cpuIndex << endl << flush;
}

bool SpatialCell::clone(const SpatialCell& s) {
   N_blocks = s.N_blocks;
   if (allocateMemory(N_blocks) == false) return false;
   // Copy cell contents to new memory locations:
   for (uint i=0; i<SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   for (uint i=0; i<SIZE_NBRS_SPA; ++i  ) cpu_nbrsSpa[i]    = s.cpu_nbrsSpa[i];
   
   for (uint i=0; i<N_blocks*SIZE_NBRS_VEL; ++i)    cpu_nbrsVel[i]     = s.cpu_nbrsVel[i];
   for (uint i=0; i<N_blocks*SIZE_BLOCKPARAMS; ++i) cpu_blockParams[i] = s.cpu_blockParams[i];
   for (uint i=0; i<N_blocks*SIZE_VELBLOCK; ++i)    cpu_avgs[i]        = s.cpu_avgs[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i)       cpu_fx[i]          = s.cpu_fx[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i)       cpu_fy[i]          = s.cpu_fy[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i)       cpu_fz[i]          = s.cpu_fz[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d1x[i]         = s.cpu_d1x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d1y[i]         = s.cpu_d1y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d1z[i]         = s.cpu_d1z[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d2x[i]         = s.cpu_d2x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d2y[i]         = s.cpu_d2y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d2z[i]         = s.cpu_d2z[i];
   return  true;
}

