#include <cstdlib>
#include <iostream>
#include <exception>
#include <limits>
#include <list>
#include <map>
#include <utility>

#include "common.h"
#include "parameters.h"
#include "grid.h"
#include "mpilogger.h"
#include "memalloc.h"

using namespace std;

extern MPILogger mpilogger;

static uint nextBlock;
static map<uint,int> referenceCount;
static list<uint> freePositions;

Grid::Grid() {
   typedef Parameters P;

   try {
      nbrsSpa     = new uint[MAX_VEL_BLOCKS*SIZE_NBRS_SPA];
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
   
   #ifdef CUDA
      cluint SIZE = MAX_VEL_BLOCKS*WID3;
   #else
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
   #endif
   /*
   try {
      avgs = new Real[SIZE];
      for (unsigned long long int i=0; i<SIZE; ++i) avgs[i] = 0.0;
   }
   catch (exception& e) {
      cerr << "\tERROR Couldn't allocate memory for avgs: " << e.what() << endl;
      exit(EXIT_FAILURE);
   }
   */
   allocateArray(&avgs,SIZE);

   #ifndef CUDA
      fx = avgs + MAX_VEL_BLOCKS*SIZE_VELBLOCK;
      fy = fx + MAX_VEL_BLOCKS*SIZE_FLUXS;
      fz = fy + MAX_VEL_BLOCKS*SIZE_FLUXS;
      d1x = fz + MAX_VEL_BLOCKS*SIZE_FLUXS;
      d1y = d1x + MAX_VEL_BLOCKS*SIZE_DERIV;
      d1z = d1y + MAX_VEL_BLOCKS*SIZE_DERIV;
      d2x = d1z + MAX_VEL_BLOCKS*SIZE_DERIV;
      d2y = d2x + MAX_VEL_BLOCKS*SIZE_DERIV;
      d2z = d2y + MAX_VEL_BLOCKS*SIZE_DERIV;
   #endif
}

Grid::~Grid() {
   delete nbrsSpa;
   delete nbrsVel;
   delete blockParams;
   //delete avgs;
   freeArray(avgs);
   
   // Check that SpatialCells have no remaining references to the memory:
   // This loop is commented out as it hangs sometimes, there is a bug somewhere.
   //It is only called after the program has finished, so disabling these cleanups is not that bad
   /*
     
   map<uint,int>::iterator it = referenceCount.begin();
   while (it != referenceCount.end()) {
      if (it->second > 0) 
	mpilogger << "Grid: Cell #" << it->first << " has " << it->second << " references remaining!" << endl << write;
      ++it;
   }
   */
   
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

uint* Grid::getNbrsSpa(cuint& cpuIndex) const {return nbrsSpa + cpuIndex*SIZE_NBRS_SPA;}
uint* Grid::getNbrsVel(cuint& cpuIndex) const {return nbrsVel + cpuIndex*SIZE_NBRS_VEL;}
Real* Grid::getBlockParams(cuint& cpuIndex) const {return blockParams + cpuIndex*SIZE_BLOCKPARAMS;}
Real* Grid::getAvgs(cuint& cpuIndex) const {return avgs + cpuIndex*SIZE_VELBLOCK;}

#ifndef CUDA
   Real* Grid::getFx(cuint& cpuIndex) const {return fx + cpuIndex*SIZE_FLUXS;}
   Real* Grid::getFy(cuint& cpuIndex) const {return fy + cpuIndex*SIZE_FLUXS;}
   Real* Grid::getFz(cuint& cpuIndex) const {return fz + cpuIndex*SIZE_FLUXS;}
   Real* Grid::getD1x(cuint& cpuIndex) const {return d1x + cpuIndex*SIZE_DERIV;}
   Real* Grid::getD1y(cuint& cpuIndex) const {return d1y + cpuIndex*SIZE_DERIV;}
   Real* Grid::getD1z(cuint& cpuIndex) const {return d1z + cpuIndex*SIZE_DERIV;}
   Real* Grid::getD2x(cuint& cpuIndex) const {return d2x + cpuIndex*SIZE_DERIV;}
   Real* Grid::getD2y(cuint& cpuIndex) const {return d2y + cpuIndex*SIZE_DERIV;}
   Real* Grid::getD2z(cuint& cpuIndex) const {return d2z + cpuIndex*SIZE_DERIV;}
#endif

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



