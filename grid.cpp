#include <cstdlib>
#include <iostream>
#include <limits>

#include "common.h"
#include "grid.h"
#include "logger.h"
#include "cudafuncs.h"

using namespace std;

extern Logger logger;

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
   if (deviceCreateArray(nbrsVel,MAX_VEL_BLOCKS*SIZE_NBRS_VEL*sizeof(uint)) == false) {
      initialized = false;
      logger << "\t Failed to allocate page-locked memory for nbrsVel!" << endl;
   } else {
      logger << "\t Allocated " << MAX_VEL_BLOCKS*SIZE_NBRS_VEL*sizeof(uint) << " bytes of page-locked CPU memory for array nbrsVel" << endl;
   }

   if (deviceCreateArray(blockParams,MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS*sizeof(real)) == false) {
      initialized = false;
      logger << "\t Failed to allocate page-locked memory for blockParams!" << endl;
   } else {
     logger << "\t Allocated " << MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS*sizeof(real) << " bytes of page-locked CPU memory for array blockParams" << endl;
   }
   
   if (deviceCreateArray(blockArray,MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(real)) == false) {
      initialized = false;
      logger << "\t Failed to allocate page-locked memory for blockArray!" << endl;
   } else {
      logger << "\t Allocated " << MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(real) << " bytes of page-locked CPU memory for array blockArray" << endl;
   }

   if (deviceCreateArray(cellParams,MAX_SPA_CELLS*SIZE_CELLPARAMS*sizeof(real)) == false) {
      initialized = false;
      logger << "\t Failed to allocate page-locked memory for cellParams!" << endl;
   } else {
      logger << "\t Allocated " << MAX_SPA_CELLS*SIZE_CELLPARAMS*sizeof(real) << " bytes of page-locked CPU memory for array cellParams" << endl;
   }
   
   if (initialized == false) {
      logger << "\t Failed to allocate all requested memory!" << endl;
      deviceDeleteArray(nbrsVel);
      deviceDeleteArray(blockParams);
      deviceDeleteArray(blockArray);
      deviceDeleteArray(cellParams);
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
   
   deviceDeleteArray(nbrsVel);
   deviceDeleteArray(blockParams);
   deviceDeleteArray(blockArray);
   deviceDeleteArray(cellParams);
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
   
   ++nextSpaCell;
   nextVelBlock += velBlocks;
   return nextSpaCell-1;
}

