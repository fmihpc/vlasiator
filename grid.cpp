#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <utility>

#include "common.h"
#include "grid.h"
#include "mpilogger.h"
#include "parameters.h"

using namespace std;

extern MPILogger mpilogger;

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
      for (unsigned long long int i=0; i<SIZE; ++i) avgs[i] = 0.0;
   }
   catch (exception& e) {
      cerr << "\tERROR Couldn't allocate memory for avgs: " << e.what() << endl;
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
	mpilogger << "Grid: Cell #" << it->first << " has " << it->second << " references remaining!" << endl << write;
      ++it;
   }
}

bool Grid::addReference(cuint& INDEX) {
   #ifdef DEBUG
      mpilogger << "Grid::addReference called for index " << INDEX;
      if (referenceCount.find(INDEX) != referenceCount.end()) {
	 mpilogger << ", now has " << referenceCount[INDEX]+1 << " refs" << endl << write;
      } else {
	 mpilogger << ", now has 1 refs" << endl << write;
      }
   #endif
   if (referenceCount.find(INDEX) == referenceCount.end()) return false;
   ++referenceCount[INDEX];
   return true;
}

uint* Grid::getNbrsVel(cuint& cpuIndex) const {return nbrsVel + cpuIndex*SIZE_NBRS_VEL;}
Real* Grid::getBlockParams(cuint& cpuIndex) const {return blockParams + cpuIndex*SIZE_BLOCKPARAMS;}
Real* Grid::getAvgs(cuint& cpuIndex) const {return avgs + cpuIndex*SIZE_VELBLOCK;}
Real* Grid::getFx(cuint& cpuIndex) const {return fx + cpuIndex*SIZE_FLUXS;}
Real* Grid::getFy(cuint& cpuIndex) const {return fy + cpuIndex*SIZE_FLUXS;}
Real* Grid::getFz(cuint& cpuIndex) const {return fz + cpuIndex*SIZE_FLUXS;}
Real* Grid::getD1x(cuint& cpuIndex) const {return d1x + cpuIndex*SIZE_DERIV;}
Real* Grid::getD1y(cuint& cpuIndex) const {return d1y + cpuIndex*SIZE_DERIV;}
Real* Grid::getD1z(cuint& cpuIndex) const {return d1z + cpuIndex*SIZE_DERIV;}
Real* Grid::getD2x(cuint& cpuIndex) const {return d2x + cpuIndex*SIZE_DERIV;}
Real* Grid::getD2y(cuint& cpuIndex) const {return d2y + cpuIndex*SIZE_DERIV;}
Real* Grid::getD2z(cuint& cpuIndex) const {return d2z + cpuIndex*SIZE_DERIV;}

bool Grid::removeReference(cuint& INDEX) {
   #ifdef DEBUG
      mpilogger << "Grid: removeReference called for index " << INDEX;
      if (referenceCount.find(INDEX) != referenceCount.end()) {
	 mpilogger << ", " << referenceCount[INDEX]-1 << " refs remaining" << endl << write;
      } else {
	 mpilogger << " CELL NOT FOUND!" << endl << write;
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
      mpilogger << "Grid ERROR: getFreeMemory called with BLOCKS=0, aborting." << endl << write;
      exit(1);
   }
   
   if (freePositions.size() > 0) {
      // Try to get available from holes in the arrays:
      uint rvalue = freePositions.front();
      freePositions.pop_front();
      ++referenceCount[rvalue];
      #ifdef DEBUG
         mpilogger << "Grid getFreeMemory: Cell " << rvalue << " returned, now has " << referenceCount[rvalue] << " refs." << endl << write;
      #endif
      return rvalue;
   } else if (nextBlock+BLOCKS <= MAX_VEL_BLOCKS) {
      // Return the next free position:
      uint rvalue = nextBlock;
      nextBlock += BLOCKS;
      ++referenceCount[rvalue];
      #ifdef DEBUG
         mpilogger << "Grid getFreeMemory: Cell " << rvalue << " returned, now has " << referenceCount[rvalue] << " refs." << endl << write;
      #endif
      return rvalue;
   } else {
      // Out of memory:
      //logger << "Out of memory" << endl << flush;
      return numeric_limits<uint>::max();
   }
}

void Grid::printReferences() {
   mpilogger << "Grid: printReferences:" << endl;
   map<uint,int>::iterator it = referenceCount.begin();
   while (it != referenceCount.end()) {
      mpilogger << "\t Cell " << it->first << " refs = " << it->second << endl;
      ++it;
   }
   
   mpilogger << "     freePositions:" << endl;
   list<uint>::iterator itt = freePositions.begin();
   while (itt != freePositions.end()) {
      mpilogger << "\t Position #" << *itt << " is free" << endl;
      ++itt;
   }
   mpilogger << write;
}

uint Grid::getTotalNumberOfBlocks() const {
   return nextBlock;
}



