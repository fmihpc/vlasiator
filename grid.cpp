#include <cstdlib>
#include <iostream>
#include <limits>

#include "common.h"
#include "grid.h"
#include "logger.h"

#ifndef NOCUDA
  #include "cudafuncs.h"
#endif

using namespace std;

extern Logger logger;

#ifndef NOCUDA
bool Grid::allocateArray(const string& name,const size_t& BYTES,Real*& arrptr) {
   if (deviceCreateArray(arrptr,BYTES) == false) {
      logger << "\t Failed to allocate page-locked memory for " << name << endl;
      arrptr = NULL;
      return false;
   } else {
      logger << "\t Allocated " << BYTES << " bytes of page-locked memory for array " << name << endl;
      return true;
   }
}

bool Grid::allocateArray(const string& name,const size_t& BYTES,uint*& arrptr) {
   if (deviceCreateArray(arrptr,BYTES) == false) {
      logger << "\t Failed to allocate page-locked memory for " << name << endl;
      arrptr = NULL;
      return false;
   } else {
      logger << "\t Allocated " << BYTES << " bytes of page-locked memory for array " << name << endl;
      return true;
   }
}
#endif

/** Creates a new Grid and allocates arrays for storing velocity grid blocks and other 
 * data structures in CPU memory.
 */
Grid::Grid() {
   logger << "(GRID): Initializing" << endl;
   initialized = true;

   // Allocate memory for spatial cells:
   nbrsSpa = new uint[MAX_SPA_CELLS*SIZE_NBRS_SPA];
   cells = new SpatialCell[MAX_SPA_CELLS];

   logger << "\t Allocated " << MAX_SPA_CELLS*SIZE_NBRS_SPA*sizeof(uint) << " bytes for array nbrsSpa" << endl;
   logger << "\t Allocated " << MAX_SPA_CELLS*sizeof(SpatialCell) << " bytes for array cells" << endl;
   
   // Allocate memory for velocity blocks:
   nbrsVel = NULL;
   blockParams = NULL;
   blockArray = NULL;
   #ifndef NOCUDA
     if (allocateArray("nbrsVel",MAX_VEL_BLOCKS*SIZE_NBRS_VEL*sizeof(uint),nbrsVel) == false) initialized = false;
     if (allocateArray("blockParams",MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS*sizeof(Real),blockParams) == false) initialized = false;
     if (allocateArray("blockArray",MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(Real),blockArray) == false) initialized = false;
     if (allocateArray("cellParams",MAX_SPA_CELLS*SIZE_CELLPARAMS*sizeof(Real),cellParams) == false) initialized = false;
     if (allocateArray("fx",MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(Real),fx) == false) initialized = false;
     if (allocateArray("fy",MAX_VEL_BLOCKS*SIZE_FLUXS*sizeof(Real),fy) == false) initialized = false;
     if (allocateArray("fz",MAX_VEL_BLOCKS*SIZE_FLUXS*sizeof(Real),fz) == false) initialized = false;
     if (allocateArray("d1x",MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(Real),d1x) == false) initialized = false;
     if (allocateArray("d1y",MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(Real),d1y) == false) initialized = false;
     if (allocateArray("d1z",MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(Real),d1z) == false) initialized = false;
     if (allocateArray("d2x",MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(Real),d2x) == false) initialized = false;
     if (allocateArray("d2y",MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(Real),d2y) == false) initialized = false;
     if (allocateArray("d2z",MAX_VEL_BLOCKS*SIZE_DERIV*sizeof(Real),d2z) == false) initialized = false;
   #else
     nbrsVel = new uint[MAX_VEL_BLOCKS*SIZE_NBRS_VEL];
     blockParams = new Real[MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS];
     blockArray = new Real[MAX_VEL_BLOCKS*SIZE_VELBLOCK];
     cellParams = new Real[MAX_SPA_CELLS*SIZE_CELLPARAMS];
     fx = new Real[MAX_VEL_BLOCKS*SIZE_FLUXS];
     fy = new Real[MAX_VEL_BLOCKS*SIZE_FLUXS];
     fz = new Real[MAX_VEL_BLOCKS*SIZE_FLUXS];
     d1x = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
     d1y = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
     d1z = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
     d2x = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
     d2y = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
     d2z = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
   #endif
   
   if (initialized == false) {
      logger << "\t Failed to allocate all requested memory!" << endl;
      #ifndef NOCUDA
        deviceDeleteArray(nbrsVel);
        deviceDeleteArray(blockParams);
        deviceDeleteArray(blockArray);
        deviceDeleteArray(cellParams);
      #endif
      delete nbrsSpa;
      delete cells;
      nbrsSpa = NULL;
      cellParams = NULL;
      cells = NULL;
      nextVelBlock = 0;
      nextSpaCell = 0;
   } else {
      nextVelBlock = 0;
      nextSpaCell = 0;
   }
}

/** Free up the memory used by the Grid. The destructor does not deallocate 
 * anything in the device memory.
 */
Grid::~Grid() {
   logger << "(GRID): Freeing up CPU grid." << endl;
   #ifndef NOCUDA
     deviceDeleteArray(nbrsVel);
     deviceDeleteArray(blockParams);
     deviceDeleteArray(blockArray);
     deviceDeleteArray(cellParams);
     deviceDeleteArray(fx);
     deviceDeleteArray(fy);
     deviceDeleteArray(fz);
     deviceDeleteArray(d1x);
     deviceDeleteArray(d1y);
     deviceDeleteArray(d1z);
     deviceDeleteArray(d2x);
     deviceDeleteArray(d2y);
     deviceDeleteArray(d2z);
   #else
     delete nbrsVel;
     delete blockParams;
     delete blockArray;
     delete cellParams;
     delete fx;
     delete fy;
     delete fz;
     delete d1x;
     delete d2x;
     delete d1y;
     delete d2y;
     delete d1z;
     delete d2z;
   #endif
   delete nbrsSpa;
   delete [] cells;
   
   nbrsSpa = NULL;
   cellParams = NULL;
   cells = NULL;
}

/** Initialize a spatial cell and reserve memory for the velocity blocks and 
 * other data structures.
 * @param velBlocks The number of velocity grid blocks to reserve for the cell.
 * @return The index of the initialized cell, or numeric_limits<uint>::max() if the 
 * requested cell could not be initialized.
 */
uint Grid::getSpatialCell(const uint& velBlocks) {
   if (initialized == false) return numeric_limits<uint>::max();
   if (nextVelBlock >= MAX_VEL_BLOCKS) return numeric_limits<uint>::max();
   if (nextSpaCell >= MAX_SPA_CELLS) return numeric_limits<uint>::max();
   
   cells[nextSpaCell].N_blocks = velBlocks;
   cells[nextSpaCell].cellIndex = nextSpaCell;
   cells[nextSpaCell].velBlockIndex = nextVelBlock;
   
   // Set pointers to arrays containing cell parameters (common for each velocity block):
   cells[nextSpaCell].cpu_cellParams = cellParams + nextSpaCell*SIZE_CELLPARAMS;
   cells[nextSpaCell].cpu_nbrsSpa    = nbrsSpa    + nextSpaCell*SIZE_NBRS_SPA;
   // Set pointers to arrays containing parameters of each velocity block:
   cells[nextSpaCell].cpu_nbrsVel     = nbrsVel     + nextVelBlock*SIZE_NBRS_VEL;
   cells[nextSpaCell].cpu_avgs        = blockArray  + nextVelBlock*SIZE_VELBLOCK;
   cells[nextSpaCell].cpu_blockParams = blockParams + nextVelBlock*SIZE_BLOCKPARAMS;
   cells[nextSpaCell].cpu_fx          = fx + nextVelBlock*SIZE_FLUXS;
   cells[nextSpaCell].cpu_fy          = fy + nextVelBlock*SIZE_FLUXS;
   cells[nextSpaCell].cpu_fz          = fz + nextVelBlock*SIZE_FLUXS;
   cells[nextSpaCell].cpu_d1x         = d1x + nextVelBlock*SIZE_DERIV;
   cells[nextSpaCell].cpu_d1y         = d1y + nextVelBlock*SIZE_DERIV;
   cells[nextSpaCell].cpu_d1z         = d1z + nextVelBlock*SIZE_DERIV;
   cells[nextSpaCell].cpu_d2x         = d2x + nextVelBlock*SIZE_DERIV;
   cells[nextSpaCell].cpu_d2y         = d2y + nextVelBlock*SIZE_DERIV;
   cells[nextSpaCell].cpu_d2z         = d2z + nextVelBlock*SIZE_DERIV;

   ++nextSpaCell;
   nextVelBlock += velBlocks;
   return nextSpaCell-1;
}

