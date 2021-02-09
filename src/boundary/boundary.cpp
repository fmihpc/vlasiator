/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/*! \file boundary.cpp
 * \brief Implementation of the class Boundary.
 */

#include <cstdlib>
#include <iostream>

#include "../grid.h"
#include "../object_wrapper.h"
#include "../vlasovsolver/cpu_moments.h"

#include "boundary.h"
#include "ionosphere.h"
#include "maxwellian.h"
#include "nothing.h"
#include "outflow.h"
#include "user.h"

using namespace std;
using namespace spatial_cell;

bool precedenceSort(const BC::BoundaryCondition *first, const BC::BoundaryCondition *second) {
   if (first->getPrecedence() < second->getPrecedence())
      return true;
   else
      return false;
}

// ************************************************************
// ***** DEFINITIONS FOR BOUNDARY CLASS *****
// ************************************************************

/*! Constructor for class Boundary.*/
Boundary::Boundary() : anyDynamic(false) {}

/*!\brief Destructor for class Boundary.
 *
 * Reduces the value of Boundary::nBoundaries by one, and if after the
 * destruction Boundary::nBoundaries equals zero all stored Boundaries are
 * deleted.
 */
Boundary::~Boundary() {
   // Call delete for each BoundaryCondition:
   for (list<BC::BoundaryCondition *>::iterator it = boundaries.begin(); it != boundaries.end(); it++) {
      delete *it;
      *it = NULL;
   }
}

/*!\brief Add its own and all existing BoundaryConditions' parameters.
 *
 * Adds the parameters specific to the Boundary class handling the list of
 * BoundaryConditions and then calls the static addParameters functions of all
 * BoundaryConditions implemented in the code in order to have them appear also
 * in the help.
 */
void Boundary::addParameters() {
   Readparameters::addComposing(
       "boundaries.boundary",
       "List of boundary condition (BC) types to be used. Each boundary condition to be used has to be on a new line "
       "boundary = YYY. Available options are: Outflow, Ionosphere, Maxwellian, User.");
   Readparameters::add("boundaries.periodic_x", "Set the grid periodicity in x-direction. true(default)/false.", true);
   Readparameters::add("boundaries.periodic_y", "Set the grid periodicity in y-direction. true(default)/false.", true);
   Readparameters::add("boundaries.periodic_z", "Set the grid periodicity in z-direction. true(default)/false.", true);

   // Call static addParameter functions in all BC's
   BC::Nothing::addParameters();
   BC::Ionosphere::addParameters();
   BC::Outflow::addParameters();
   BC::Maxwellian::addParameters();
   BC::User::addParameters();
}

/*!\brief Get this class' parameters.
 *
 * Get actually used boundary conditions, parameters of which are read from
 * each BoundaryCondition's initialization function.
 */
void Boundary::getParameters() {
   Readparameters::get("boundaries.boundary", boundaryCondList);
   Readparameters::get("boundaries.periodic_x", periodic[0]);
   Readparameters::get("boundaries.periodic_y", periodic[1]);
   Readparameters::get("boundaries.periodic_z", periodic[2]);
}

/*! Add a new BC::BoundaryCondition.
 * Boundary will take care of deleting it.
 * \param bc BoundaryCondition object
 * \param project Project object
 * \param t Current time
 */
void Boundary::addBoundary(BC::BoundaryCondition *bc, Project &project, creal t) {
   // Initialize the boundary condition
   bc->initBoundary(t, project);

   boundaries.push_back(bc);
   if (boundaries.size() > 1)
      boundaries.sort(precedenceSort);

   // This assumes that only one instance of each type is created.
   indexToBoundary[bc->getIndex()] = bc;
}

/*!\brief Initialise all boundary conditions actually used.
 *
 * This function loops through the list of boundary conditions listed as to be
 * used in the configuration file/command line arguments. For each of these it
 * adds the corresponding instance and updates the member anyDynamic to
 * determine whether any BoundaryCondition is dynamic in time.
 * \param project Project object
 * \param t Current time
 * \sa addBoundary
 */
void Boundary::initBoundaries(Project &project, creal t) {
   vector<string>::const_iterator it;

   if (boundaryCondList.size() == 0) {
      if (!periodic[0] && !Readparameters::helpRequested)
         BC::abort_mpi("Non-periodic in x but no boundary condtion loaded!");
      if (!periodic[1] && !Readparameters::helpRequested)
         BC::abort_mpi("Non-periodic in y but no boundary condtion loaded!");
      if (!periodic[2] && !Readparameters::helpRequested)
         BC::abort_mpi("Non-periodic in z but no boundary condtion loaded!");
   }

   for (it = boundaryCondList.begin(); it != boundaryCondList.end(); it++) {
      if (*it == "Outflow") {
         this->addBoundary(new BC::Outflow, project, t);

         bool faces[6];
         this->getBoundary(boundarytype::OUTFLOW)->getFaces(&faces[0]);

         if ((faces[0] || faces[1]) && periodic[0])
            BC::abort_mpi("Conflict: x boundaries set to periodic but found Outflow conditions!");

         if ((faces[2] || faces[3]) && periodic[1])
            BC::abort_mpi("Conflict: y boundaries set to periodic but found Outflow conditions!");

         if ((faces[4] || faces[5]) && periodic[2])
            BC::abort_mpi("Conflict: z boundaries set to periodic but found Outflow conditions!");

         if ((faces[0] || faces[1]) && P::xcells_ini < 5)
            BC::abort_mpi("Outflow condition loaded on x- or x+ face but not enough cells in x!");

         if ((faces[2] || faces[3]) && P::ycells_ini < 5)
            BC::abort_mpi("Outflow condition loaded on y- or y+ face but not enough cells in y!");

         if ((faces[4] || faces[5]) && P::zcells_ini < 5)
            BC::abort_mpi("Outflow condition loaded on z- or z+ face but not enough cells in z!");
      }
      if (*it == "Ionosphere") {
         this->addBoundary(new BC::Ionosphere, project, t);
         this->addBoundary(new BC::Nothing, project, t);
         anyDynamic = anyDynamic || this->getBoundary(boundarytype::IONOSPHERE)->isDynamic();
      }
      if (*it == "Maxwellian") {
         this->addBoundary(new BC::Maxwellian, project, t);
         anyDynamic = anyDynamic || this->getBoundary(boundarytype::MAXWELLIAN)->isDynamic();
         bool faces[6];
         this->getBoundary(boundarytype::MAXWELLIAN)->getFaces(&faces[0]);

         if ((faces[0] || faces[1]) && periodic[0])
            BC::abort_mpi("Conflict: x boundaries set to periodic but found Maxwellian also!");

         if ((faces[2] || faces[3]) && periodic[1])
            BC::abort_mpi("Conflict: y boundaries set to periodic but found Maxwellian also!");

         if ((faces[4] || faces[5]) && periodic[2])
            BC::abort_mpi("Conflict: z boundaries set to periodic but found Maxwellian also!");

         if ((faces[0] || faces[1]) && P::xcells_ini < 5)
            BC::abort_mpi("Maxwellian condition loaded on x- or x+ face but not enough cells in x!");

         if ((faces[2] || faces[3]) && P::ycells_ini < 5)
            BC::abort_mpi("Maxwellian condition loaded on y- or y+ face but not enough cells in y!");

         if ((faces[4] || faces[5]) && P::zcells_ini < 5)
            BC::abort_mpi("Maxwellian condition loaded on z- or z+ face but not enough cells in z!");
      }
   }

   list<BC::BoundaryCondition *>::iterator it2;
   for (it2 = boundaries.begin(); it2 != boundaries.end(); it2++)
      (*it2)->setPeriodicity(periodic);
}

/* Verifies that all cells within FULL_NEIGHBORHOOD_ID of L1 boundary cells are on the same refinement
 * level (one group for inner boundary, another for outer boundary). */
void Boundary::checkRefinement(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid) {
   // Set is used such that each cell will only be checked once.
   set<CellID> innerBoundaryCells, outerBoundaryCells;

   int innerBoundaryRefLvl = -1, outerBoundaryRefLvl = -1;

   // Collect cells by boundarytype
   for (auto cellId : mpiGrid.get_cells()) {
      SpatialCell *cell = mpiGrid[cellId];
      if (cell) {
         if (cell->boundaryFlag == boundarytype::IONOSPHERE) {
            innerBoundaryCells.insert(cellId);
            innerBoundaryRefLvl = mpiGrid.get_refinement_level(cellId);
            if (cell->boundaryLayer == 1) {
               // Add all stencil neighbors of layer 1 cells
               auto *nbrPairVector = mpiGrid.get_neighbors_of(cellId, BOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
               for (auto nbrPair : *nbrPairVector) {
                  if (nbrPair.first != INVALID_CELLID)
                     innerBoundaryCells.insert(nbrPair.first);
               }
            }
         } else if (cell->boundaryFlag != boundarytype::NOT_BOUNDARY &&
                    cell->boundaryFlag != boundarytype::NOTHING) {
            outerBoundaryCells.insert(cellId);
            outerBoundaryRefLvl = mpiGrid.get_refinement_level(cellId);
            // Add all stencil neighbors of outer boundary cells
            auto *nbrPairVector = mpiGrid.get_neighbors_of(cellId, BOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
            for (auto nbrPair : *nbrPairVector) {
               if (nbrPair.first != INVALID_CELLID)
                  outerBoundaryCells.insert(nbrPair.first);
            }
         }
      }
   }

   for (auto cellId : innerBoundaryCells) {
      if (cellId != INVALID_CELLID && mpiGrid.get_refinement_level(cellId) != innerBoundaryRefLvl)
         BC::abort_mpi("ERROR: inner boundary cells must have identical refinement level!");
   }

   for (auto cellId : outerBoundaryCells) {
      if (cellId != INVALID_CELLID && mpiGrid.get_refinement_level(cellId) != outerBoundaryRefLvl)
         BC::abort_mpi("ERROR: outer boundary cells must have identical refinement level!");
   }
}

bool belongsToLayer(const int layer, const int x, const int y, const int z,
                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid) {
   bool belongs = false;

   // loop through all neighbors (including diagonals)
   for (int ix = -1; ix <= 1; ++ix) {
      for (int iy = -1; iy <= 1; ++iy) {
         for (int iz = -1; iz <= 1; ++iz) {
            // not strictly necessary but logically we should not consider the cell itself
            // among its neighbors.
            if ((ix == 0 && iy == 0 && iz == 0) || !technicalGrid.get(x + ix, y + iy, z + iz))
               continue;

            if (layer == 1 && technicalGrid.get(x + ix, y + iy, z + iz)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
               // in the first layer, boundary cell belongs if it has a non-boundary neighbor
               belongs = true;
               return belongs;
            } else if (layer > 1 && technicalGrid.get(x + ix, y + iy, z + iz)->boundaryLayer == layer - 1) {
               // in all other layers, boundary cell belongs if it has a neighbor in the previous layer
               belongs = true;
               return belongs;
            }
         }
      }
   }

   return belongs;
}

/*!\brief Classify all simulation cells with respect to the boundary conditions.
 *
 * Loops through all cells and and for each assigns the correct boundaryFlag
 * depending on the return value of each BoundaryCondition's assignBoundary.
 * \param mpiGrid Grid
 */
void Boundary::classifyCells(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                             FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid) {
   vector<CellID> cells = mpiGrid.get_cells();
   auto localSize = technicalGrid.getLocalSize().data();

   // Set all cells to default value
   for (uint i = 0; i < cells.size(); i++)
      mpiGrid[cells[i]]->boundaryFlag = boundarytype::NOT_BOUNDARY;

#pragma omp parallel for collapse(3)
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            technicalGrid.get(x, y, z)->boundaryFlag = boundarytype::NOT_BOUNDARY;
            technicalGrid.get(x, y, z)->boundaryLayer = 0;
            technicalGrid.get(x, y, z)->maxFsDt = std::numeric_limits<Real>::max();
            // Set the fsgrid rank in the technical grid
            technicalGrid.get(x, y, z)->fsGridRank = technicalGrid.getRank();
         }
      }
   }

   /*
     loop through boundaries and let all boundaries set in local
   cells if they are part of which boundary (cell location needs to
   be updated by now. No remote data needed/available, so assignement
   has to be based individually on each cells location
   */
   list<BC::BoundaryCondition *>::iterator it;
   for (it = boundaries.begin(); it != boundaries.end(); it++)
      (*it)->assignBoundary(mpiGrid, technicalGrid);

   // communicate boundary assignments (boundaryFlag and
   // boundaryLayer communicated)
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_BOUNDARYFLAG);
   mpiGrid.update_copies_of_remote_neighbors(BOUNDARIES_NEIGHBORHOOD_ID);

   // set distance 1 cells to boundary cells, that have neighbors which are normal cells
   for (uint i = 0; i < cells.size(); i++) {
      mpiGrid[cells[i]]->boundaryLayer = 0; /*Initial value*/
      if (mpiGrid[cells[i]]->boundaryFlag != boundarytype::NOT_BOUNDARY) {
         const auto *nbrs = mpiGrid.get_neighbors_of(cells[i], BOUNDARIES_NEIGHBORHOOD_ID);
         for (uint j = 0; j < (*nbrs).size(); j++) {
            if ((*nbrs)[j].first != 0) {
               if (mpiGrid[(*nbrs)[j].first]->boundaryFlag == boundarytype::NOT_BOUNDARY) {
                  mpiGrid[cells[i]]->boundaryLayer = 1;
               }
            }
         }
      }
   }

   /*communicate which cells have Layer 1 set above for local cells (boundaryFlag
    * and boundaryLayer communicated)*/
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_BOUNDARYFLAG);

   mpiGrid.update_copies_of_remote_neighbors(BOUNDARIES_NEIGHBORHOOD_ID);

   /*Compute distances*/
   uint maxLayers = 3; // max(max(P::xcells_ini, P::ycells_ini), P::zcells_ini);
   for (uint layer = 1; layer < maxLayers; layer++) {
      for (uint i = 0; i < cells.size(); i++) {
         if (mpiGrid[cells[i]]->boundaryLayer == 0) {
            const auto *nbrs = mpiGrid.get_neighbors_of(cells[i], BOUNDARIES_NEIGHBORHOOD_ID);
            // Note: this distance calculation will be non-plateau monotonic only assuming that
            // Boundary::checkRefinement has been applied correctly and there are no refinement
            // level changes within BOUNDARIES_NEIGHBORHOOD_ID.
            for (uint j = 0; j < (*nbrs).size(); j++) {
               if ((*nbrs)[j].first != 0 && (*nbrs)[j].first != cells[i]) {
                  if (mpiGrid[(*nbrs)[j].first]->boundaryLayer == layer) {
                     mpiGrid[cells[i]]->boundaryLayer = layer + 1;
                     break;
                  }
               }
            }
         }
      }

      SpatialCell::set_mpi_transfer_type(Transfer::CELL_BOUNDARYFLAG);
      mpiGrid.update_copies_of_remote_neighbors(BOUNDARIES_NEIGHBORHOOD_ID);
   }

   /*set cells to NOTHING if they are on boundary, and are not
    * in the first two layers of the boundary*/
   for (uint i = 0; i < cells.size(); i++) {
      if (mpiGrid[cells[i]]->boundaryFlag != boundarytype::NOT_BOUNDARY && mpiGrid[cells[i]]->boundaryLayer != 1 &&
          mpiGrid[cells[i]]->boundaryLayer != 2) {
         mpiGrid[cells[i]]->boundaryFlag = boundarytype::NOTHING;
      }
   }

   // The following is done so that everyone knows their neighbour's layer
   // flags.  This is needed for the correct use of the boundary local
   // communication patterns.
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_BOUNDARYFLAG);
   mpiGrid.update_copies_of_remote_neighbors(FULL_NEIGHBORHOOD_ID);

   // Now the layers need to be set on fsgrid too
   // In dccrg initialization the max number of boundary layers is set to 3.
   const uint MAX_NUMBER_OF_BOUNDARY_LAYERS = 3 * pow(2, mpiGrid.get_maximum_refinement_level());

   technicalGrid.updateGhostCells();

   // loop through max number of layers
   for (uint layer = 1; layer <= MAX_NUMBER_OF_BOUNDARY_LAYERS; ++layer) {

// loop through all cells in grid
#pragma omp parallel for collapse(3)
      for (int x = 0; x < localSize[0]; ++x) {
         for (int y = 0; y < localSize[1]; ++y) {
            for (int z = 0; z < localSize[2]; ++z) {

               // for the first layer, consider all cells that belong to a boundary, for other layers
               // consider all cells that have not yet been labeled.
               if ((layer == 1 && technicalGrid.get(x, y, z)->boundaryFlag != boundarytype::NOT_BOUNDARY) ||
                   (layer > 1 && technicalGrid.get(x, y, z)->boundaryLayer == 0)) {

                  if (belongsToLayer(layer, x, y, z, technicalGrid)) {

                     technicalGrid.get(x, y, z)->boundaryLayer = layer;

                     if (layer > 2 && technicalGrid.get(x, y, z)->boundaryFlag != boundarytype::NOT_BOUNDARY) {
                        technicalGrid.get(x, y, z)->boundaryFlag = boundarytype::NOTHING;
                     }
                  }
               }
            }
         }
      }
      technicalGrid.updateGhostCells();
   }

// One more pass to make sure, in particular if the ionosphere is wide enough
// there is remaining cells of IONOSPHERE type inside the max layers gone through previously.
// This last pass now gets rid of them.
#pragma omp parallel for collapse(3)
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            if (technicalGrid.get(x, y, z)->boundaryLayer == 0 &&
                technicalGrid.get(x, y, z)->boundaryFlag == boundarytype::IONOSPHERE) {
               technicalGrid.get(x, y, z)->boundaryFlag = boundarytype::NOTHING;
            }
         }
      }
   }

   technicalGrid.updateGhostCells();

   const std::array<int, 3> fsGridDimensions = technicalGrid.getGlobalSize();

// One pass to setup the bit field to know which components the field solver should propagate.
#pragma omp parallel for collapse(3)
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            technicalGrid.get(x, y, z)->SOLVE = 0;

            std::array<int32_t, 3> globalIndices = technicalGrid.getGlobalIndices(x, y, z);

            if (((globalIndices[0] == 0 || globalIndices[0] == fsGridDimensions[0] - 1) && !this->isPeriodic(0)) ||
                ((globalIndices[1] == 0 || globalIndices[1] == fsGridDimensions[1] - 1) && !this->isPeriodic(1)) ||
                ((globalIndices[2] == 0 || globalIndices[2] == fsGridDimensions[2] - 1) && !this->isPeriodic(2))) {
               continue;
            }
            if (technicalGrid.get(x, y, z)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BX;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BY;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BZ;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EX;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EY;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EZ;
            } else {
               if (technicalGrid.get(x - 1, y, z)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BX;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EY;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EZ;
               }
               if (technicalGrid.get(x, y - 1, z)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BY;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EX;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EZ;
               }
               if (technicalGrid.get(x, y, z - 1)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BZ;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EX;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EY;
               }
               if (technicalGrid.get(x - 1, y - 1, z)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EZ;
               }
               if (technicalGrid.get(x - 1, y, z - 1)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EY;
               }
               if (technicalGrid.get(x, y - 1, z - 1)->boundaryFlag == boundarytype::NOT_BOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EX;
               }
            }
         }
      }
   }

   technicalGrid.updateGhostCells();
}

/*!\brief Apply the initial state to all boundary cells.
 * Loops through all BoundaryConditions and apply the initial state for all
 * existing particle species.
 * \param mpiGrid Grid
 */
void Boundary::applyInitialState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                 FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                                 Project &project) {
   list<BC::BoundaryCondition *>::iterator it;
   for (it = boundaries.begin(); it != boundaries.end(); it++) {
      // Skip when not restarting or not requested.
      if (Parameters::isRestart && !(*it)->doApplyUponRestart())
         continue;
      (*it)->applyInitialState(mpiGrid, perBGrid, project);
   }
}

void Boundary::updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                           FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid, creal t) {
   list<BC::BoundaryCondition *>::iterator it;
   if (anyDynamic) {
      for (it = boundaries.begin(); it != boundaries.end(); it++) {
         // Skip when not restarting or not requested.
         if (Parameters::isRestart && !(*it)->doApplyUponRestart())
            continue;
         if ((*it)->isDynamic())
            (*it)->updateState(mpiGrid, perBGrid, t);
      }
   }
}

/*!\brief Apply the Vlasov boundary conditions to all boundary cells at time t.
 *
 * Loops through all BoundaryConditions and calls the corresponding
 * vlasovBoundaryCondition() function. The boundary condition functions are
 * called for all particle species, one at a time.
 *
 * WARNING (see end of the function) Blocks are changed but lists not updated
 * now, if you need to use/communicate them before the next update is done, add
 * an update at the end of this function.
 *
 * \param mpiGrid Grid
 * \param t Current time
 * \param doCalcMomentsV Compute into _V moments if true or _R moments if false
 * so that the interpolated ones can be done.
 */
void Boundary::applyBoundaryVlasovConditions(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid, creal t,
                                             const bool doCalcMomentsV) {
   if (boundaries.size() == 0)
      return; // no boundaries

      /*Transfer along boundaries*/
      // First the small stuff without overlapping in an extended neighbourhood:
#warning TODO This now communicates in the wider neighbourhood for both layers, could be reduced to smaller neighbourhood for layer 1, larger neighbourhood for layer 2.
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_PARAMETERS | Transfer::POP_METADATA | Transfer::CELL_BOUNDARYFLAG,
                                      true);
   mpiGrid.update_copies_of_remote_neighbors(BOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);

   // Loop over existing particle species
   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
      SpatialCell::setCommunicatedSpecies(popID);
      // update lists in larger neighborhood
      updateRemoteVelocityBlockLists(mpiGrid, popID, BOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);

      // Then the block data in the reduced neighbourhood:
      int timer = phiprof::initializeTimer("Start comm of cell and block data", "MPI");
      phiprof::start(timer);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA, true);
      mpiGrid.start_remote_neighbor_copy_updates(BOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
      phiprof::stop(timer);

      timer = phiprof::initializeTimer("Compute process inner cells");
      phiprof::start(timer);

      // Compute Vlasov boundary condition on boundary/process inner cells
      vector<CellID> localCells;
      getBoundaryCellList(mpiGrid, mpiGrid.get_local_cells_not_on_process_boundary(BOUNDARIES_EXTENDED_NEIGHBORHOOD_ID),
                          localCells);

#pragma omp parallel for
      for (uint i = 0; i < localCells.size(); i++) {
         cuint boundaryType = mpiGrid[localCells[i]]->boundaryFlag;
         this->getBoundary(boundaryType)->vlasovBoundaryCondition(mpiGrid, localCells[i], popID, doCalcMomentsV);
      }
      if (doCalcMomentsV)
         calculateMoments_V(mpiGrid, localCells, true);
      else
         calculateMoments_R(mpiGrid, localCells, true);

      phiprof::stop(timer);

      timer = phiprof::initializeTimer("Wait for receives", "MPI", "Wait");
      phiprof::start(timer);
      mpiGrid.wait_remote_neighbor_copy_update_receives(BOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
      phiprof::stop(timer);

      // Compute vlasov boundary on boundary/process boundary cells
      timer = phiprof::initializeTimer("Compute process boundary cells");
      phiprof::start(timer);
      vector<CellID> boundaryCells;
      getBoundaryCellList(mpiGrid, mpiGrid.get_local_cells_on_process_boundary(BOUNDARIES_EXTENDED_NEIGHBORHOOD_ID),
                          boundaryCells);
#pragma omp parallel for
      for (uint i = 0; i < boundaryCells.size(); i++) {
         cuint boundaryType = mpiGrid[boundaryCells[i]]->boundaryFlag;
         this->getBoundary(boundaryType)->vlasovBoundaryCondition(mpiGrid, boundaryCells[i], popID, doCalcMomentsV);
      }
      if (doCalcMomentsV)
         calculateMoments_V(mpiGrid, boundaryCells, true);
      else
         calculateMoments_R(mpiGrid, boundaryCells, true);

      phiprof::stop(timer);

      timer = phiprof::initializeTimer("Wait for sends", "MPI", "Wait");
      phiprof::start(timer);
      mpiGrid.wait_remote_neighbor_copy_update_sends();
      phiprof::stop(timer);

      // WARNING Blocks are changed but lists not updated now, if you need to
      // use/communicate them before the next update is done, add an update
      // here. Reset lists in smaller default neighborhood.
      updateRemoteVelocityBlockLists(mpiGrid, popID);

   } // for-loop over populations
}

/*! Get a pointer to the BoundaryCondition of given index.
 * \param boundaryType Type of the boundary condition to return
 * \return Pointer to the instance of the BoundaryCondition
 */
BC::BoundaryCondition *Boundary::getBoundary(cuint boundaryType) const {
   auto it = indexToBoundary.find(boundaryType);
   if (it != indexToBoundary.end())
      return it->second;
   else
      BC::abort_mpi("ERROR: Boundary " + to_string(boundaryType) + " is invalid", 1);
}

/*! Get the number of BoundaryConditions stored in Boundary. */
unsigned int Boundary::size() const { return boundaries.size(); }

/*! Check whether any boundary condition is dynamic in time. */
bool Boundary::isDynamic() const { return anyDynamic; }

/*! Get a bool telling whether the system is periodic in the queried direction.
 * \param direction 0: x, 1: y, 2: z.
 * \retval periodic Is the system periodic in the queried direction.
 */
bool Boundary::isPeriodic(uint direction) const { return periodic[direction]; }

/*! Get a vector containing the cellID of all cells which are not NOTHING or
 * NOT_BOUNDARY in the vector of cellIDs passed to the function.
 * \param mpiGrid Grid
 * \param cellList Vector of cellIDs in which to look for boundary cells
 * \param boundaryCellList Vector of boundary the cells' cellIDs
 */
void getBoundaryCellList(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                         const vector<uint64_t> &cellList, vector<uint64_t> &boundaryCellList) {
   boundaryCellList.clear();
   for (size_t cell = 0; cell < cellList.size(); ++cell) {
      const CellID cellID = cellList[cell];
      if (mpiGrid[cellID]->boundaryFlag == boundarytype::NOTHING ||
          mpiGrid[cellID]->boundaryFlag == boundarytype::NOT_BOUNDARY)
         continue;
      boundaryCellList.push_back(cellID);
   }
}

/*! Updates all NonboundaryCells into an internal map. This should be called in
 * loadBalance.
 * \param mpiGrid The DCCRG grid
 */
void Boundary::updateBoundariesAfterLoadBalance(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid) {
   vector<uint64_t> local_cells_on_boundary;
   getBoundaryCellList(mpiGrid, mpiGrid.get_cells(), local_cells_on_boundary);
   // Loop over boundaries:
   for (std::list<BC::BoundaryCondition *>::iterator it = boundaries.begin(); it != boundaries.end(); ++it) {
      (*it)->updateBoundaryConditionsAfterLoadBalance(mpiGrid, local_cells_on_boundary);
   }
}
