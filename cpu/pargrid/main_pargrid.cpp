#ifdef PARGRID

#include <cstdlib>
#include <iostream>
#include <vector>

#include <mpilogger.h>
#include <pargrid.h>
#include <definitions.h>
#include <cell_spatial.h>
#include <project.h>
#include <timer.h>
#include <vlsreader.h>
#include "vlasovmover.h"

using namespace std;

// NOTE: If preprocessor flag PROFILE is undefined, the compiler should optimize out Timer:: calls.

extern MPILogger mpilogger;

extern bool cpu_acceleration(SpatialCell& cell);
extern bool cpu_translation1(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation2(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation3(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_calcVelocityMoments(SpatialCell& cell);

inline uchar calcNbrTypeID(cuchar& i,cuchar& j,cuchar& k) {return k*25+j*5+i;}

namespace Main {
   std::vector<ID::type> cells;
   std::vector<const SpatialCell*> nbrPtrs(6,NULL);
   SpatialCell* cellPtr;
   
   uint calcAcc;
   uint calcSpatDerivs;
   uint spatDerivsMPIRecv;
   uint spatDerivsMPISend;
   uint calcSpatFluxes;
   uint spatFluxesMPIRecv;
   uint spatFluxesMPISend;
   uint calcSpatProp;
   uint spatPropMPIRecv;
   uint spatPropMPISend;
}

bool finalizeMover() {
   Timer::print();
   return true;
}

bool initializeMover(ParGrid<SpatialCell>& mpiGrid) {
   Main::calcAcc           = Timer::create("Computing: vel. propagation  (total) : ");
   Main::calcSpatDerivs    = Timer::create("Computing: spat. derivatives (total) : ");
   Main::spatDerivsMPIRecv = Timer::create("MPI Recv : spat. derivs              : ");
   Main::spatDerivsMPISend = Timer::create("MPI Send : spat. derivs              : ");
   Main::calcSpatFluxes    = Timer::create("Computing: spat. fluxes      (total) : ");
   Main::spatFluxesMPIRecv = Timer::create("MPI Recv : spat. fluxes              : ");
   Main::spatFluxesMPISend = Timer::create("MPI Send : spat. fluxes              : ");
   Main::calcSpatProp      = Timer::create("Computing: spat. propag      (total) : ");
   Main::spatPropMPIRecv   = Timer::create("MPI Recv : spat. propag              : ");
   Main::spatPropMPISend   = Timer::create("MPI Send : spat. propag              : ");

   // Populate spatial neighbour list:
   mpiGrid.getAllCells(Main::cells);
   for (size_t cell=0; cell<Main::cells.size(); ++cell) {
      cuint cellID = Main::cells[cell];
      uint* const nbrsSpa = mpiGrid[cellID]->cpu_nbrsSpa;

      // Get spatial neighbour IDs and store them into a vector:
      uint counter = 0;
      vector<uint> nbrIDs;
      for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
	 if (i == 0 & (j == 0 & k == 0)) nbrIDs.push_back(cellID); // in ParGrid cells do not consider themselves as their own neighbours
	 else nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(2+i,2+j,2+k)));
	 ++counter;
      }
      nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(0,2,2))); // i-1,j,k nbr, goes to index 27
      nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,0,2))); // i,j-2,k nbr, goes to index 28
      nbrIDs.push_back(mpiGrid.getNeighbour(cellID,calcNbrTypeID(2,2,0))); // i,j,k-2 nbr, goes to index 29

      // Store neighbour offsets into a vector:
      vector<uint> cellOffsets(nbrIDs.size());
      for (size_t i=0; i<nbrIDs.size(); ++i) {
	 if (nbrIDs[i] == numeric_limits<uint>::max()) cellOffsets[i] = numeric_limits<uint>::max();
	 else cellOffsets[i] = mpiGrid[nbrIDs[i]]->cpuIndex * SIZE_VELBLOCK;
      }
      
      // Create spatial neighbour list entry for each block:
      for (uint block=0; block<mpiGrid[cellID]->N_blocks; ++block) {
	 uint boundaryFlag = 0;
	 // Store offsets to each spatial neighbour of this block. Note that 
	 // the offset to this block is stored to index 13:
	 for (size_t i=0; i<nbrIDs.size(); ++i) {
	    if (cellOffsets[i] == numeric_limits<uint>::max()) {
	       nbrsSpa[block*SIZE_NBRS_SPA + i] = numeric_limits<uint>::max();
	    } else {
	       boundaryFlag = boundaryFlag | (1 << i);
	       nbrsSpa[block*SIZE_NBRS_SPA + i] = cellOffsets[i] + block*SIZE_VELBLOCK;
	    }
	 }
	 // Boundary flags are stored to the last position in nbrsSpa array:
	 nbrsSpa[block*SIZE_NBRS_SPA + 30] = boundaryFlag;
      }
   }
   
   return true;
}

void initialLoadBalance(ParGrid<SpatialCell>& mpiGrid) { 
   
}

bool findNeighbours(std::vector<const SpatialCell*>& nbrPtr,const ParGrid<SpatialCell>& mpiGrid,const ID::type& CELLID) {
   ID::type nbrID;
   for (int i=0; i<6; ++i) nbrPtr[i] = NULL;
   nbrID = mpiGrid.getNeighbour(CELLID,calcNbrTypeID(2-1,2  ,2  ));
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[0] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID,calcNbrTypeID(2+1,2  ,2  ));
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[1] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID,calcNbrTypeID(2  ,2-1,2  ));
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[2] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID,calcNbrTypeID(2  ,2+1,2  ));
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[3] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID,calcNbrTypeID(2  ,2  ,2-1));
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[4] = mpiGrid[nbrID];
   nbrID = mpiGrid.getNeighbour(CELLID,calcNbrTypeID(2  ,2  ,2+1));
   if (nbrID != std::numeric_limits<ID::type>::max()) nbrPtr[5] = mpiGrid[nbrID];
   
   return true;
}

void calculateVelocityMoments(ParGrid<SpatialCell>& mpiGrid) {
   mpiGrid.getCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (Main::cellPtr != NULL) cpu_calcVelocityMoments(*Main::cellPtr);
   }
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

void calculateSpatialPropagation(ParGrid<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs) {
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

#endif // #ifdef PARGRID
