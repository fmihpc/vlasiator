#ifndef PARGRID // Do not use these functions if ParGrid is used:

#ifndef MAIN_H
#define MAIN_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <dccrg.hpp>
#include <boost/mpi.hpp>
#include <zoltan.h>

#include "logger.h"
#include "definitions.h"
#include "parameters.h"
#include "cell_spatial.h"

namespace Main {
   std::vector<uint64_t> cells;
   std::vector<const SpatialCell*> nbrPtrs(6,NULL);
   SpatialCell* cellPtr;
}

Logger logger;

extern bool cpu_acceleration(SpatialCell& cell);
extern bool cpu_translation1(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation2(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation3(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);

void initialLoadBalance(dccrg<SpatialCell>& mpiGrid) {
   mpiGrid.balance_load();
   comm.barrier();
}

bool findNeighbours(std::vector<const SpatialCell*>& nbrPtr,dccrg<SpatialCell>& mpiGrid,const uint64_t& cell) {
   std::vector<uint64_t> nbrs;
   // Get a pointer to each neighbour. If the neighbour does not exists, 
   // a NULL pointer is inserted.
   nbrs = mpiGrid.get_neighbours_x(cell,-1.0);
   if (nbrs.size() > 0) nbrPtr[0] = mpiGrid[nbrs[0]]; // Get ptr to -x neighbour
   else nbrPtr[0] = NULL;
   nbrs = mpiGrid.get_neighbours_x(cell,+1.0);
   if (nbrs.size() > 0) nbrPtr[1] = mpiGrid[nbrs[0]]; //            +x
   else nbrPtr[1] = NULL;
   nbrs = mpiGrid.get_neighbours_y(cell,-1.0);
   if (nbrs.size() > 0) nbrPtr[2] = mpiGrid[nbrs[0]]; //            -y
   else nbrPtr[2] = NULL;
   nbrs = mpiGrid.get_neighbours_y(cell,+1.0);
   if (nbrs.size() > 0) nbrPtr[3] = mpiGrid[nbrs[0]]; //            +y
   else nbrPtr[3] = NULL;
   nbrs = mpiGrid.get_neighbours_z(cell,-1.0);
   if (nbrs.size() > 0) nbrPtr[4] = mpiGrid[nbrs[0]]; //            -z
   else nbrPtr[4] = NULL;
   nbrs = mpiGrid.get_neighbours_z(cell,+1.0);
   if (nbrs.size() > 0) nbrPtr[5] = mpiGrid[nbrs[0]]; //            +z
   else nbrPtr[5] = NULL;
   return true;
}

void calculateAcceleration(dccrg<SpatialCell>& mpiGrid) {
   // Calculate acceleration for all cells (inner + boundary):
   Main::cells = mpiGrid.get_cells();
   for (size_t i=0; i<Main::cells.size(); ++i) {
      Main::cellPtr = mpiGrid[Main::cells[i]];
      if (Main::cellPtr != NULL) cpu_acceleration(*Main::cellPtr);
   }   
}

void calculateSpatialDerivatives(dccrg<SpatialCell>& mpiGrid) {
   typedef Parameters P;
   // Start neighbour data exchange:
   P::transmit = Transmit::AVGS;
   mpiGrid.start_remote_neighbour_data_update();
   // Calculate derivatives for inner cells:
   Main::cells = mpiGrid.get_cells_with_local_neighbours();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }
   // Calculate derivatives for boundary cells when transfers have completed:
   mpiGrid.wait_neighbour_data_update();
   Main::cells = mpiGrid.get_cells_with_remote_neighbour();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	   std::cerr << "Failed to find neighbours." << std::endl; 
	   continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }
}

void calculateSpatialFluxes(dccrg<SpatialCell>& mpiGrid) {
   typedef Parameters P;
   // Start neighbour data exchange:
   P::transmit = Transmit::AVGS | Transmit::DERIV1 | Transmit::DERIV2;
   mpiGrid.start_remote_neighbour_data_update();
   // Calculate fluxes for inner cells:
   Main::cells = mpiGrid.get_cells_with_local_neighbours();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
   // Calculate fluxes for boundary cells when transfers have completed:
   mpiGrid.wait_neighbour_data_update();
   Main::cells = mpiGrid.get_cells_with_remote_neighbour();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
}

void calculateSpatialPropagation(dccrg<SpatialCell>& mpiGrid) {
   typedef Parameters P;
   // Start neighbour data exchange:
   P::transmit = Transmit::FLUXES;
   mpiGrid.start_remote_neighbour_data_update();
   // Propagate inner cells:
   Main::cells = mpiGrid.get_cells_with_local_neighbours();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
   }
   // Propagate boundary cells when transfers have completed:
   mpiGrid.wait_neighbour_data_update();
   Main::cells = mpiGrid.get_cells_with_remote_neighbour();
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
#endif // #ifndef PARGRID
