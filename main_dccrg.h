#ifndef PARGRID // Do not use these functions if ParGrid is used:

#ifndef MAIN_H
#define MAIN_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <boost/mpi.hpp>
#include <zoltan.h>

#define DCCRG_SEND_SINGLE_CELLS
#define DCCRG_CELL_DATA_SIZE_FROM_USER
#define DCCRG_USER_MPI_DATA_TYPE
#include <dccrg.hpp>

#include "logger.h"
#include "definitions.h"
#include "parameters.h"
#include "cell_spatial.h"
#include "project.h"
#include "timer.h"

namespace Main {
   std::vector<uint64_t> cells;
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

Logger logger;

extern bool cpu_acceleration(SpatialCell& cell);
extern bool cpu_translation1(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation2(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation3(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);

void initialLoadBalance(dccrg<SpatialCell>& mpiGrid) {
   mpiGrid.balance_load();
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

void calculateSimParameters(dccrg<SpatialCell>& mpiGrid, creal& t, Real& dt) {
   // TODO let the project function decide if something should really be calculated
   if (!cellParametersChanged(t)) {
   	return;
   }
   calcSimParameters(mpiGrid, t, dt);
}

void calculateCellParameters(dccrg<SpatialCell>& mpiGrid,creal& t, uint64_t cell) {
   // TODO let the project function decide if something should really be calculated
   if (!cellParametersChanged(t)) {
   	return;
   }
   calcCellParameters(mpiGrid[cell]->cpu_cellParams,t);
}

void calculateAcceleration(dccrg<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcAcc);
   
   // Calculate acceleration for all cells (inner + boundary):
   Main::cells = mpiGrid.get_cells();
   for (size_t i=0; i<Main::cells.size(); ++i) {
      Main::cellPtr = mpiGrid[Main::cells[i]];
      if (Main::cellPtr != NULL) cpu_acceleration(*Main::cellPtr);
   }
   
   Timer::stop(Main::calcAcc);
}

void calculateSpatialDerivatives(dccrg<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcSpatDerivs);
   
   typedef Parameters P;
   // Start neighbour data exchange:
   P::transmit = Transmit::AVGS;
   SpatialCell::base_address_identifier = 0;
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
   Timer::start(Main::spatDerivsMPIRecv);
   mpiGrid.wait_neighbour_data_update();
   Timer::stop(Main::spatDerivsMPIRecv);
   
   Main::cells = mpiGrid.get_cells_with_remote_neighbour();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	   std::cerr << "Failed to find neighbours." << std::endl; 
	   continue;
      }
      if (Main::cellPtr != NULL) cpu_translation1(*Main::cellPtr,Main::nbrPtrs);
   }
   
   Timer::start(Main::spatDerivsMPISend);
   #warning MPI sends cannot be timed with DCCRG, recv contains sends and receives
   Timer::stop(Main::spatDerivsMPISend);
   Timer::stop(Main::calcSpatDerivs);
}

void calculateSpatialFluxes(dccrg<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcSpatFluxes);
   
   typedef Parameters P;
   // Start neighbour data exchange:
   P::transmit = Transmit::AVGS | Transmit::DERIV1 | Transmit::DERIV2;
   SpatialCell::base_address_identifier = 1;
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
   Timer::start(Main::spatFluxesMPIRecv);
   mpiGrid.wait_neighbour_data_update();
   Timer::stop(Main::spatFluxesMPIRecv);
   
   Main::cells = mpiGrid.get_cells_with_remote_neighbour();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation2(*Main::cellPtr,Main::nbrPtrs);
   }
   
   Timer::start(Main::spatFluxesMPISend);
   #warning MPI sends cannot be timed with DCCRG, recv contains sends and receives
   Timer::stop(Main::spatFluxesMPISend);
   Timer::stop(Main::calcSpatFluxes);
}

void calculateSpatialPropagation(dccrg<SpatialCell>& mpiGrid) {
   Timer::start(Main::calcSpatProp);
   
   typedef Parameters P;
   // Start neighbour data exchange:
   P::transmit = Transmit::FLUXES;
   SpatialCell::base_address_identifier = 2;
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
   Timer::start(Main::spatPropMPIRecv);
   mpiGrid.wait_neighbour_data_update();
   Timer::stop(Main::spatPropMPIRecv);
   
   Main::cells = mpiGrid.get_cells_with_remote_neighbour();
   for (size_t c=0; c<Main::cells.size(); ++c) {
      Main::cellPtr = mpiGrid[Main::cells[c]];
      if (findNeighbours(Main::nbrPtrs,mpiGrid,Main::cells[c]) == false) {
	 std::cerr << "Failed to find neighbours." << std::endl; 
	 continue;
      }
      if (Main::cellPtr != NULL) cpu_translation3(*Main::cellPtr,Main::nbrPtrs);
   }
   
   Timer::start(Main::spatPropMPISend);
   #warning MPI sends cannot be timed with DCCRG, recv contains sends and receives
   Timer::stop(Main::spatPropMPISend);
   Timer::stop(Main::calcSpatProp);
}

void writeTimers() {Timer::print();}

#endif // #ifndef MAIN_H
#endif // #ifndef PARGRID
