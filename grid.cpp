#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <utility>

#include "common.h"
#include "grid.h"
#include "logger.h"
#include "parameters.h"

using namespace std;

extern Logger logger;

static uint nextBlock;
static map<uint,int> referenceCount;
static list<uint> freePositions;

Grid::Grid() {
   typedef Parameters P;

   try {
      nbrsVel     = new uint[MAX_VEL_BLOCKS*SIZE_NBRS_VEL];
   }
   catch (exception& e) {
      cerr << "Couldn't allocate memory for nbrsVel: " << e.what() << endl;
      exit(EXIT_FAILURE);
   }

   try {
      blockParams = new Real[MAX_VEL_BLOCKS*SIZE_BLOCKPARAMS];
   }
   catch (exception& e) {
      cerr << "Couldn't allocate memory for blockParams: " << e.what() << endl;
      exit(EXIT_FAILURE);
   }
   /*
   avgs        = new Real[MAX_VEL_BLOCKS*SIZE_VELBLOCK];
   fx          = new Real[MAX_VEL_BLOCKS*SIZE_FLUXS];
   fy          = new Real[MAX_VEL_BLOCKS*SIZE_FLUXS];
   fz          = new Real[MAX_VEL_BLOCKS*SIZE_FLUXS];
   d1x         = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
   d1y         = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
   d1z         = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
   d2x         = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
   d2y         = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
   d2z         = new Real[MAX_VEL_BLOCKS*SIZE_DERIV];
   */
   cluint SIZE = 
     MAX_VEL_BLOCKS*SIZE_VELBLOCK
     + MAX_VEL_BLOCKS*SIZE_FLUXS
     + MAX_VEL_BLOCKS*SIZE_FLUXS
     + MAX_VEL_BLOCKS*SIZE_FLUXS
     + MAX_VEL_BLOCKS*SIZE_DERIV
     + MAX_VEL_BLOCKS*SIZE_DERIV
     + MAX_VEL_BLOCKS*SIZE_DERIV
     + MAX_VEL_BLOCKS*SIZE_DERIV
     + MAX_VEL_BLOCKS*SIZE_DERIV
     + MAX_VEL_BLOCKS*SIZE_DERIV;
   
   try {
      avgs = new Real[SIZE];
   }
   catch (exception& e) {
      cerr << "Couldn't allocate memory for avgs: " << e.what() << endl;
      exit(EXIT_FAILURE);
   }

   fx = avgs + MAX_VEL_BLOCKS*SIZE_VELBLOCK;
   fy = fx + MAX_VEL_BLOCKS*SIZE_FLUXS;
   fz = fy + MAX_VEL_BLOCKS*SIZE_FLUXS;
   d1x = fz + MAX_VEL_BLOCKS*SIZE_FLUXS;
   d1y = d1x + MAX_VEL_BLOCKS*SIZE_DERIV;
   d1z = d1y + MAX_VEL_BLOCKS*SIZE_DERIV;
   d2x = d1z + MAX_VEL_BLOCKS*SIZE_DERIV;
   d2y = d2x + MAX_VEL_BLOCKS*SIZE_DERIV;
   d2z = d2y + MAX_VEL_BLOCKS*SIZE_DERIV;
}

Grid::~Grid() {
   delete nbrsVel;
   delete blockParams;
   delete avgs;
   /*
   delete fx;
   delete fy;
   delete fz;
   delete d1x;
   delete d1y;
   delete d1z;
   delete d2x;
   delete d2y;
   delete d2z;
    */
   // Check that SpatialCells have no remaining references to the memory:
   map<uint,int>::iterator it = referenceCount.begin();
   while (it != referenceCount.end()) {
      if (it->second > 0) 
	logger << "Grid: Cell #" << it->first << " has " << it->second << " references remaining!" << endl << flush;
      ++it;
   }
}

bool Grid::addReference(cuint& INDEX) {
   #ifdef DEBUG
      logger << "Grid::addReference called for index " << INDEX;
      if (referenceCount.find(INDEX) != referenceCount.end()) {
	 logger << ", now has " << referenceCount[INDEX]+1 << " refs" << endl << flush;
      } else {
	 logger << ", now has 1 refs" << endl << flush;
      }
   #endif
   if (referenceCount.find(INDEX) == referenceCount.end()) return false;
   ++referenceCount[INDEX];
   return true;
}

bool Grid::removeReference(cuint& INDEX) {
   #ifdef DEBUG
      logger << "Grid: removeReference called for index " << INDEX;
      if (referenceCount.find(INDEX) != referenceCount.end()) {
	 logger << ", " << referenceCount[INDEX]-1 << " refs remaining" << endl << flush;
      } else {
	 logger << " CELL NOT FOUND!" << endl << flush;
      }
   #endif
   if (referenceCount.find(INDEX) == referenceCount.end()) return false;
   --referenceCount[INDEX];
   if (referenceCount[INDEX] <= 0) { // Mark the positions as free
      referenceCount[INDEX] = 0;
      freePositions.push_back(INDEX);
   }
   return true;
}

uint Grid::getFreeMemory(cuint& BLOCKS) {
   if (BLOCKS == 0) {
      logger << "Grid ERROR: getFreeMemory called with BLOCKS=0, aborting." << endl << flush;
      exit(1);
   }
   
   if (freePositions.size() > 0) {
      // Try to get available from holes in the arrays:
      uint rvalue = freePositions.front();
      freePositions.pop_front();
      ++referenceCount[rvalue];
      #ifdef DEBUG
         logger << "Grid getFreeMemory: Cell " << rvalue << " returned, now has " << referenceCount[rvalue] << " refs." << endl;
      #endif
      return rvalue;
   } else if (nextBlock+BLOCKS <= MAX_VEL_BLOCKS) {
      // Return the next free position:
      uint rvalue = nextBlock;
      nextBlock += BLOCKS;
      ++referenceCount[rvalue];
      #ifdef DEBUG
         logger << "Grid getFreeMemory: Cell " << rvalue << " returned, now has " << referenceCount[rvalue] << " refs." << endl;
      #endif
      return rvalue;
   } else {
      // Out of memory:
      //logger << "Out of memory" << endl << flush;
      return numeric_limits<uint>::max();
   }
}

void Grid::printReferences() {
   logger << "Grid: printReferences:" << endl;
   map<uint,int>::iterator it = referenceCount.begin();
   while (it != referenceCount.end()) {
      logger << "\t Cell " << it->first << " refs = " << it->second << endl;
      ++it;
   }
   logger << flush;
   
   logger << "     freePositions:" << endl;
   list<uint>::iterator itt = freePositions.begin();
   while (itt != freePositions.end()) {
      logger << "\t Position #" << *itt << " is free" << endl;
      ++itt;
   }
   logger << flush;
}

uint Grid::getTotalNumberOfBlocks() const {
   return nextBlock;
}



