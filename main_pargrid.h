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

Logger logger;

extern bool cpu_acceleration(SpatialCell& cell);
extern bool cpu_translation1(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation2(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation3(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);

namespace Main {
   std::vector<ID::type> cells;
   std::vector<const SpatialCell*> nbrPtrs(6,NULL);
   SpatialCell* cellPtr;
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

void calculateCellParameters(ParGrid<SpatialCell>& mpiGrid,creal& t) {
   // First check if the cell parameters really need to be recalculated:
   if (cellParametersChanged(t) == false) return;
   mpiGrid.getCells(Main::cells);
   for (size_t i=0; i<Main::cells.size(); ++i) calcCellParameters(mpiGrid[Main::cells[i]]->cpu_cellParams,t);
}

void calculateAcceleration(ParGrid<SpatialCell>& mpiGrid) {
   // Calculate acceleration for all cells (inner + boundary):
   mpiGrid.getCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (Main::cellPtr != NULL) cpu_acceleration(*Main::cellPtr);
   }
}

void calculateSpatialDerivatives(ParGrid<SpatialCell>& mpiGrid) {
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
   mpiGrid.waitAll();
   mpiGrid.getBoundaryCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }
}

void calculateSpatialFluxes(ParGrid<SpatialCell>& mpiGrid) {
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
   mpiGrid.waitAll();
   mpiGrid.getBoundaryCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
}

void calculateSpatialPropagation(ParGrid<SpatialCell>& mpiGrid) {
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
   mpiGrid.waitAll();
   mpiGrid.getBoundaryCells(Main::cells);
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl;
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
   }
}
   
#endif // #ifndef MAIN_H
#endif // #ifdef PARGRID
