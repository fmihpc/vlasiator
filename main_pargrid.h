#ifdef PARGRID
#ifndef MAIN_H
#define MAIN_H

#include <cstdlib>
#include <iostream>
#include <vector>

#include "logger.h"
#include "pargrid.h"
#include "definitions.h"
#include "cell_spatial.h"
#include "project.h"
#include "timer.h"

// NOTE: If preprocessor flag PROFILE is undefined, the compiler should optimize out Timer:: calls.

Logger logger;

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
   if (nbrs.size() != 6) {
      logger << "findNeighbours ERROR: Failed to get neighbour indices!" << std::endl;
      return false;
   }
   // Get a pointer to each neighbour. If the neighbour does not exists, 
   // a NULL pointer is inserted:
   if (nbrs[0] != std::numeric_limits<ID::type>::max()) nbrPtr[0] = mpiGrid[nbrs[0]];
   else nbrPtr[0] = NULL;
   if (nbrs[1] != std::numeric_limits<ID::type>::max()) nbrPtr[1] = mpiGrid[nbrs[1]];
   else nbrPtr[1] = NULL;
   if (nbrs[2] != std::numeric_limits<ID::type>::max()) nbrPtr[2] = mpiGrid[nbrs[2]];
   else nbrPtr[2] = NULL;
   if (nbrs[3] != std::numeric_limits<ID::type>::max()) nbrPtr[3] = mpiGrid[nbrs[3]];
   else nbrPtr[3] = NULL;
   if (nbrs[4] != std::numeric_limits<ID::type>::max()) nbrPtr[4] = mpiGrid[nbrs[4]];
   else nbrPtr[4] = NULL;
   if (nbrs[5] != std::numeric_limits<ID::type>::max()) nbrPtr[5] = mpiGrid[nbrs[5]];
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
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }
   // Calculate derivatives for boundary cells when transfers have completed:
   Timer::start(Main::spatDerivsMPIRecv);
   mpiGrid.waitAllReceives();
   Timer::stop(Main::spatDerivsMPIRecv);
   
   mpiGrid.getBoundaryCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }
   
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
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
   // Calculate fluxes for boundary cells when transfers have completed:
   Timer::start(Main::spatFluxesMPIRecv);
   mpiGrid.waitAllReceives();
   Timer::stop(Main::spatFluxesMPIRecv);
   
   mpiGrid.getBoundaryCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
   
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
	 std::cerr << "Failed to find neighbours." << std::endl;
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
   }
   // Propagate boundary cells when transfers have completed:
   Timer::start(Main::spatPropMPIRecv);
   mpiGrid.waitAllReceives();
   Timer::stop(Main::spatPropMPIRecv);
   
   mpiGrid.getBoundaryCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl;
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
   }
   
   Timer::start(Main::spatPropMPISend);
   mpiGrid.waitAllSends();
   Timer::stop(Main::spatPropMPISend);
   Timer::stop(Main::calcSpatProp);
}

void writeTimers() {Timer::print();}

#endif // #ifndef MAIN_H
#endif // #ifdef PARGRID
