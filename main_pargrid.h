#ifdef PARGRID
#ifndef MAIN_H
#define MAIN_H

#include <cstdlib>
#include <iostream>
#include <vector>

#include "mpilogger.h"
#include "pargrid.h"
#include "definitions.h"
#include "cell_spatial.h"
#include "project.h"
#include "timer.h"

// NOTE: If preprocessor flag PROFILE is undefined, the compiler should optimize out Timer:: calls.

MPILogger mpilogger;

extern bool cpu_acceleration(SpatialCell& cell);
extern bool cpu_translation1(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation2(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation3(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);

namespace Main {
   std::vector<ID::type> cells;
   std::vector<const SpatialCell*> nbrPtrs(6,NULL);
   SpatialCell* cellPtr;
   
   cuint calcAcc           = Timer::create("Computing: vel. propagation  (total) : ");
   cuint calcSpatDerivs    = Timer::create("Computing: spat. derivatives (total) : ");
   cuint spatDerivsMPIRecv = Timer::create("MPI Recv : spat. derivs              : ");
   cuint spatDerivsMPISend = Timer::create("MPI Send : spat. derivs              : ");
   cuint calcSpatFluxes    = Timer::create("Computing: spat. fluxes      (total) : ");
   cuint spatFluxesMPIRecv = Timer::create("MPI Recv : spat. fluxes              : ");
   cuint spatFluxesMPISend = Timer::create("MPI Send : spat. fluxes              : ");
   cuint calcSpatProp      = Timer::create("Computing: spat. propag      (total) : ");
   cuint spatPropMPIRecv   = Timer::create("MPI Recv : spat. propag              : ");
   cuint spatPropMPISend   = Timer::create("MPI Send : spat. propag              : ");
}

void initialLoadBalance(ParGrid<SpatialCell>& mpiGrid) { }

bool findNeighbours(std::vector<const SpatialCell*>& nbrPtr,const ParGrid<SpatialCell>& mpiGrid,const ID::type& CELLID) {
   std::vector<ID::type> nbrs;
   mpiGrid.getNeighbours(nbrs,CELLID);
   if (nbrs.size() != 24) {
      mpilogger << "findNeighbours ERROR: Failed to get neighbour indices!" << std::endl << write;
      return false;
   }
   // Get a pointer to each neighbour. If the neighbour does not exists, 
   // a NULL pointer is inserted:
   if (nbrs[ 0] != std::numeric_limits<ID::type>::max()) nbrPtr[0] = mpiGrid[nbrs[ 0]];
   else nbrPtr[0] = NULL;
   if (nbrs[ 4] != std::numeric_limits<ID::type>::max()) nbrPtr[1] = mpiGrid[nbrs[ 4]];
   else nbrPtr[1] = NULL;
   if (nbrs[ 8] != std::numeric_limits<ID::type>::max()) nbrPtr[2] = mpiGrid[nbrs[ 8]];
   else nbrPtr[2] = NULL;
   if (nbrs[12] != std::numeric_limits<ID::type>::max()) nbrPtr[3] = mpiGrid[nbrs[12]];
   else nbrPtr[3] = NULL;
   if (nbrs[16] != std::numeric_limits<ID::type>::max()) nbrPtr[4] = mpiGrid[nbrs[16]];
   else nbrPtr[4] = NULL;
   if (nbrs[20] != std::numeric_limits<ID::type>::max()) nbrPtr[5] = mpiGrid[nbrs[20]];
   else nbrPtr[5] = NULL;
   return true;
}

void calculateSimParameters(ParGrid<SpatialCell>& mpiGrid, creal& t, Real& dt) {
   // TODO let the project function decide if something should really be calculated
   if (!cellParametersChanged(t)) {
   	return;
   }
   calcSimParameters(mpiGrid, t, dt);
}

void calculateCellParameters(ParGrid<SpatialCell>& mpiGrid, creal& t, ID::type cell) {
   // TODO let the project function decide if something should really be calculated
   if (!cellParametersChanged(t)) {
   	return;
   }
   calcCellParameters(mpiGrid[cell]->cpu_cellParams,t);
}

void calculateAcceleration(ParGrid<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcAcc);

   // Calculate acceleration for all cells (inner + boundary):
   mpiGrid.getCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (Main::cellPtr != NULL) cpu_acceleration(*Main::cellPtr);
   }

   Timer::stop(Main::calcAcc);
}

void calculateSpatialDerivatives(ParGrid<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcSpatDerivs);

   // Start neighbour data exchange:
   mpiGrid.startNeighbourExchange(0);
   // Calculate derivatives for inner cells:
   mpiGrid.getInnerCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write;
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }

   #ifdef PARGRID_WAITANY
      // Loop until all remote cell data has been received:
      ID::type readyCellID;
      while (mpiGrid.waitAnyReceive() == true) {
	 // Calculate all ready local cells:
	 while (mpiGrid.getReadyCell(readyCellID) == true) {
	    Main::cellPtr = mpiGrid[readyCellID];
	    if (findNeighbours(Main::nbrPtrs,mpiGrid,readyCellID) == false) {
	       mpilogger << "Failed to find neighbours." << std::endl << write;
	       continue;
	    }
	    if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
	 }
      }
   #elif PARGRID_WAITSOME
      // Loop until all remote cell data has been received:
      ID::type readyCellID;
      while (mpiGrid.waitSomeReceives() == true) {
	 // Calculate all ready local cells:
	 while (mpiGrid.getReadyCell(readyCellID) == true) {
	    Main::cellPtr = mpiGrid[readyCellID];
	    if (findNeighbours(Main::nbrPtrs,mpiGrid,readyCellID) == false) {
	       mpilogger << "Failed to find neighbours." << std::endl << write;
	       continue;
	    }
	    if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
	 }
      }
   #else
      // Calculate derivatives for boundary cells when transfers have completed:
      Timer::start(Main::spatDerivsMPIRecv);
      mpiGrid.waitAllReceives();
      Timer::stop(Main::spatDerivsMPIRecv);
   
      mpiGrid.getBoundaryCells(Main::cells);
      for (size_t c=0; c<Main::cells.size(); ++c) {
	 Main::cellPtr = mpiGrid[Main::cells[c]];
	 if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	    mpilogger << "Failed to find neighbours." << std::endl << write;
	    continue;
	 }
	 if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
      }
   #endif
   // Wait for all sends to complete:
   Timer::start(Main::spatDerivsMPISend);
   mpiGrid.waitAllSends();
   Timer::stop(Main::spatDerivsMPISend);
   Timer::stop(Main::calcSpatDerivs);
}

void calculateSpatialFluxes(ParGrid<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcSpatFluxes);

   // Start neighbour data exchange:
   mpiGrid.startNeighbourExchange(1);
   // Calculate fluxes for inner cells:
   mpiGrid.getInnerCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
   
   #ifdef PARGRID_WAITSOME
      // Loop until all remote cell data has been received:
      ID::type readyCellID;
      while (mpiGrid.waitSomeReceives() == true) {
	 // Calculate all ready local cells:
	 while (mpiGrid.getReadyCell(readyCellID) == true) {
	    Main::cellPtr = mpiGrid[readyCellID];
	    if (findNeighbours(Main::nbrPtrs,mpiGrid,readyCellID) == false) {
	       mpilogger << "Failed to find neighbours." << std::endl << write;
	       continue;
	    }
	    if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
	 }
      }
   #else
      // Calculate fluxes for boundary cells when transfers have completed:
      Timer::start(Main::spatFluxesMPIRecv);
      mpiGrid.waitAllReceives();
      Timer::stop(Main::spatFluxesMPIRecv);
   
      mpiGrid.getBoundaryCells(Main::cells);
      for (size_t c=0; c<Main::cells.size(); ++c) {
	 Main::cellPtr = mpiGrid[Main::cells[c]];
	 if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	    mpilogger << "Failed to find neighbours." << std::endl << write; 
	    continue;
	 }
	 if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
      }
   #endif
   Timer::start(Main::spatFluxesMPISend);
   mpiGrid.waitAllSends();
   Timer::stop(Main::spatFluxesMPISend);
   Timer::stop(Main::calcSpatFluxes);
}

void calculateSpatialPropagation(ParGrid<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcSpatProp);

   // Start neighbour data exchange:
   mpiGrid.startNeighbourExchange(2);
   // Start neighbour data exchange:
   mpiGrid.getInnerCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 mpilogger << "Failed to find neighbours." << std::endl << write;
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
   }
   #ifdef PARGRID_WAITSOME
      // Loop until all remote cell data has been received:
      ID::type readyCellID;
      while (mpiGrid.waitSomeReceives() == true) {
	 // Calculate all ready local cells:
	 while (mpiGrid.getReadyCell(readyCellID) == true) {
	    Main::cellPtr = mpiGrid[readyCellID];
	    if (findNeighbours(Main::nbrPtrs,mpiGrid,readyCellID) == false) {
	       mpilogger << "Failed to find neighbours." << std::endl << write;
	       continue;
	    }
	    if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
	 }
      }
   #else
      // Propagate boundary cells when transfers have completed:
      Timer::start(Main::spatPropMPIRecv);
      mpiGrid.waitAllReceives();
      Timer::stop(Main::spatPropMPIRecv);
   
      mpiGrid.getBoundaryCells(Main::cells);
      for (size_t c=0; c<Main::cells.size(); ++c) {
	 Main::cellPtr = mpiGrid[Main::cells[c]];
	 if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	    mpilogger << "Failed to find neighbours." << std::endl << write;
	    continue;
	 }
	 if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
      }
   #endif
   Timer::start(Main::spatPropMPISend);
   mpiGrid.waitAllSends();
   Timer::stop(Main::spatPropMPISend);
   Timer::stop(Main::calcSpatProp);
}

void writeTimers() {Timer::print();}

#endif // #ifndef MAIN_H
#endif // #ifdef PARGRID
