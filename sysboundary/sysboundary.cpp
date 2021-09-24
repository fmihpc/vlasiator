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

/*! \file sysboundary.cpp
 * \brief Implementation of the class SysBoundary.
 */

#include <cstdlib>
#include <iostream>

#include "../grid.h"
#include "../object_wrapper.h"
#include "../vlasovsolver/cpu_moments.h"

#include "donotcompute.h"
#include "ionosphere.h"
#include "maxwellian.h"
#include "outflow.h"
#include "sysboundary.h"
#include "user.h"

using namespace std;
using namespace spatial_cell;

bool precedenceSort(const SBC::SysBoundaryCondition* first, const SBC::SysBoundaryCondition* second) {
   if (first->getPrecedence() < second->getPrecedence())
      return true;
   else
      return false;
}

// ************************************************************
// ***** DEFINITIONS FOR SYSBOUNDARY CLASS *****
// ************************************************************

SysBoundary::SysBoundary() : anyDynamic(false) {}

SysBoundary::~SysBoundary() {
   // Call delete for each SysBoundaryCondition
   for (list<SBC::SysBoundaryCondition*>::iterator it = sysBoundaries.begin(); it != sysBoundaries.end(); it++) {
      delete *it;
      *it = NULL;
   }
}

/*!\brief Add its own and all existing SysBoundaryConditions' parameters.
 *
 * Adds the parameters specific to the SysBoundary class handling the list of
 * SysBoundaryConditions and then calls the static addParameters functions of all
 * SysBoundaryConditions implemented in the code in order to have them appear also
 * in the help.
 */
void SysBoundary::addParameters() {
   Readparameters::addComposing(
       "boundaries.boundary",
       "List of boundary condition (BC) types to be used. Each boundary condition to be used has to be on a new line "
       "boundary = YYY. Available options are: Outflow, Ionosphere, Maxwellian, User.");
   Readparameters::add("boundaries.periodic_x", "Set the grid periodicity in x-direction. 'yes'(default)/'no'.", "yes");
   Readparameters::add("boundaries.periodic_y", "Set the grid periodicity in y-direction. 'yes'(default)/'no'.", "yes");
   Readparameters::add("boundaries.periodic_z", "Set the grid periodicity in z-direction. 'yes'(default)/'no'.", "yes");

   // Call static addParameter functions in all BC's
   SBC::DoNotCompute::addParameters();
   SBC::Ionosphere::addParameters();
   SBC::Outflow::addParameters();
   SBC::Maxwellian::addParameters();
   SBC::User::addParameters();
}

/*!\brief Get this class' parameters.
 *
 * Get actually used boundary conditions, parameters of which are read from
 * each SysBoundaryCondition's initialization function.
 */
void SysBoundary::getParameters() {
   string periodic_x, periodic_y, periodic_z;
   Readparameters::get("boundaries.boundary", sysBoundaryCondList);
   Readparameters::get("boundaries.periodic_x", periodic_x);
   Readparameters::get("boundaries.periodic_y", periodic_y);
   Readparameters::get("boundaries.periodic_z", periodic_z);

   periodic[0] = (periodic_x == "yes") ? true : false;
   periodic[1] = (periodic_y == "yes") ? true : false;
   periodic[2] = (periodic_z == "yes") ? true : false;
}

/*! Add a new SBC::SysBoundaryCondition.
 * SysBoundary will take care of deleting it.
 * \param bc SysBoundaryCondition object
 * \param project Project object
 * \param t Current time
 */
void SysBoundary::addBoundary(SBC::SysBoundaryCondition* bc, Project& project, creal t) {
   // Initialize the boundary condition
   bc->initBoundary(t, project);

   sysBoundaries.push_back(bc);
   if (sysBoundaries.size() > 1)
      sysBoundaries.sort(precedenceSort);

   // This assumes that only one instance of each type is created.
   indexToBoundary[bc->getIndex()] = bc;
}

/*!\brief Initialise all boundary conditions actually used.
 *
 * This function loops through the list of boundary conditions listed as to be
 * used in the configuration file/command line arguments. For each of these it
 * adds the corresponding instance and updates the member anyDynamic to
 * determine whether any SysBoundaryCondition is dynamic in time.
 * \param project Project object
 * \param t Current time
 * \sa addBoundary
 */
void SysBoundary::initBoundaries(Project& project, creal t) {
   vector<string>::const_iterator it;

   if (sysBoundaryCondList.size() == 0) {
      if (!periodic[0] && !Readparameters::helpRequested)
         SBC::abort_mpi("Non-periodic in x but no boundary condtion loaded!");
      if (!periodic[1] && !Readparameters::helpRequested)
         SBC::abort_mpi("Non-periodic in y but no boundary condtion loaded!");
      if (!periodic[2] && !Readparameters::helpRequested)
         SBC::abort_mpi("Non-periodic in z but no boundary condtion loaded!");
   }

   for (it = sysBoundaryCondList.begin(); it != sysBoundaryCondList.end(); it++) {
      if (*it == "Outflow") {
         this->addBoundary(new SBC::Outflow, project, t);

         bool faces[6];
         this->getBoundary(sysboundarytype::OUTFLOW)->getFaces(&faces[0]);

         if ((faces[0] || faces[1]) && periodic[0])
            SBC::abort_mpi("Conflict: x boundaries set to periodic but found Outflow conditions!");

         if ((faces[2] || faces[3]) && periodic[1])
            SBC::abort_mpi("Conflict: y boundaries set to periodic but found Outflow conditions!");

         if ((faces[4] || faces[5]) && periodic[2])
            SBC::abort_mpi("Conflict: z boundaries set to periodic but found Outflow conditions!");

         if ((faces[0] || faces[1]) && P::xcells_ini < 5)
            SBC::abort_mpi("Outflow condition loaded on x- or x+ face but not enough cells in x!");

         if ((faces[2] || faces[3]) && P::ycells_ini < 5)
            SBC::abort_mpi("Outflow condition loaded on y- or y+ face but not enough cells in y!");

         if ((faces[4] || faces[5]) && P::zcells_ini < 5)
            SBC::abort_mpi("Outflow condition loaded on z- or z+ face but not enough cells in z!");
      } else if (*it == "Ionosphere") {
         this->addBoundary(new SBC::Ionosphere, project, t);
         this->addBoundary(new SBC::DoNotCompute, project, t);
         anyDynamic = anyDynamic || this->getBoundary(sysboundarytype::IONOSPHERE)->isDynamic();
      } else if (*it == "Maxwellian") {
         this->addBoundary(new SBC::Maxwellian, project, t);
         anyDynamic = anyDynamic || this->getBoundary(sysboundarytype::MAXWELLIAN)->isDynamic();
         bool faces[6];
         this->getBoundary(sysboundarytype::MAXWELLIAN)->getFaces(&faces[0]);

         if ((faces[0] || faces[1]) && periodic[0])
            SBC::abort_mpi("Conflict: x boundaries set to periodic but found Maxwellian also!");

         if ((faces[2] || faces[3]) && periodic[1])
            SBC::abort_mpi("Conflict: y boundaries set to periodic but found Maxwellian also!");

         if ((faces[4] || faces[5]) && periodic[2])
            SBC::abort_mpi("Conflict: z boundaries set to periodic but found Maxwellian also!");

         if ((faces[0] || faces[1]) && P::xcells_ini < 5)
            SBC::abort_mpi("Maxwellian condition loaded on x- or x+ face but not enough cells in x!");

         if ((faces[2] || faces[3]) && P::ycells_ini < 5)
            SBC::abort_mpi("Maxwellian condition loaded on y- or y+ face but not enough cells in y!");

         if ((faces[4] || faces[5]) && P::zcells_ini < 5)
            SBC::abort_mpi("Maxwellian condition loaded on z- or z+ face but not enough cells in z!");
      } else {
         std::ostringstream msg;
         msg << "Unknown type of boundary read: " << *it;
         SBC::abort_mpi(msg.str());
      }
   }

   list<SBC::SysBoundaryCondition*>::iterator it2;
   for (it2 = sysBoundaries.begin(); it2 != sysBoundaries.end(); it2++)
      (*it2)->setPeriodicity(periodic);
}

/* Verifies that all cells within FULL_NEIGHBORHOOD_ID of L1 boundary cells are on the same refinement
 * level (one group for inner boundary, another for outer boundary). */
void SysBoundary::checkRefinement(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
   // Set is used such that each cell will only be checked once.
   set<CellID> innerBoundaryCells, outerBoundaryCells;

   int innerBoundaryRefLvl = -1, outerBoundaryRefLvl = -1;

   // Collect cells by sysboundarytype
   for (auto cellId : mpiGrid.get_cells()) {
      SpatialCell* cell = mpiGrid[cellId];
      if (cell) {
         if (cell->sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
            innerBoundaryCells.insert(cellId);
            innerBoundaryRefLvl = mpiGrid.get_refinement_level(cellId);
            if (cell->sysBoundaryLayer == 1) {
               // Add all stencil neighbors of layer 1 cells
               auto* nbrPairVector = mpiGrid.get_neighbors_of(cellId, SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
               for (auto nbrPair : *nbrPairVector) {
                  if (nbrPair.first != INVALID_CELLID)
                     innerBoundaryCells.insert(nbrPair.first);
               }
            }
         } else if (cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
                    cell->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE) {
            outerBoundaryCells.insert(cellId);
            outerBoundaryRefLvl = mpiGrid.get_refinement_level(cellId);
            // Add all stencil neighbors of outer boundary cells
            auto* nbrPairVector = mpiGrid.get_neighbors_of(cellId, SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
            for (auto nbrPair : *nbrPairVector) {
               if (nbrPair.first != INVALID_CELLID)
                  outerBoundaryCells.insert(nbrPair.first);
            }
         }
      }
   }

   for (auto cellId : innerBoundaryCells) {
      if (cellId != INVALID_CELLID && mpiGrid.get_refinement_level(cellId) != innerBoundaryRefLvl)
         SBC::abort_mpi("ERROR: inner boundary cells must have identical refinement level!");
   }

   for (auto cellId : outerBoundaryCells) {
      if (cellId != INVALID_CELLID && mpiGrid.get_refinement_level(cellId) != outerBoundaryRefLvl)
         SBC::abort_mpi("ERROR: outer boundary cells must have identical refinement level!");
   }
}

bool belongsToLayer(const int layer, const int x, const int y, const int z,
                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   bool belongs = false;

   // loop through all neighbors (including diagonals)
   for (int ix = -1; ix <= 1; ++ix) {
      for (int iy = -1; iy <= 1; ++iy) {
         for (int iz = -1; iz <= 1; ++iz) {
            // not strictly necessary but logically we should not consider the cell itself
            // among its neighbors.
            if ((ix == 0 && iy == 0 && iz == 0) || !technicalGrid.get(x + ix, y + iy, z + iz))
               continue;

            if (layer == 1 && technicalGrid.get(x + ix, y + iy, z + iz)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
               // in the first layer, boundary cell belongs if it has a non-boundary neighbor
               belongs = true;
               return belongs;
            } else if (layer > 1 && technicalGrid.get(x + ix, y + iy, z + iz)->sysBoundaryLayer == layer - 1) {
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
 * Loops through all cells and and for each assigns the correct sysBoundaryFlag
 * depending on the return value of each SysBoundaryCondition's assignBoundary.
 * \param mpiGrid Grid
 */
void SysBoundary::classifyCells(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   vector<CellID> cells = mpiGrid.get_cells();
   auto localSize = technicalGrid.getLocalSize().data();

   // Set all cells to default value
   for (uint i = 0; i < cells.size(); i++)
      mpiGrid[cells[i]]->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;

#pragma omp parallel for collapse(3)
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            technicalGrid.get(x, y, z)->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
            technicalGrid.get(x, y, z)->sysBoundaryLayer = 0;
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
   list<SBC::SysBoundaryCondition*>::iterator it;
   for (it = sysBoundaries.begin(); it != sysBoundaries.end(); it++)
      (*it)->assignBoundary(mpiGrid, technicalGrid);

   // communicate boundary assignments (sysBoundaryFlag and
   // sysBoundaryLayer communicated)
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
   mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_NEIGHBORHOOD_ID);

   // set distance 1 cells to boundary cells, that have neighbors which are normal cells
   for (uint i = 0; i < cells.size(); i++) {
      mpiGrid[cells[i]]->sysBoundaryLayer = 0; /*Initial value*/
      if (mpiGrid[cells[i]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
         const auto* nbrs = mpiGrid.get_neighbors_of(cells[i], SYSBOUNDARIES_NEIGHBORHOOD_ID);
         for (uint j = 0; j < (*nbrs).size(); j++) {
            if ((*nbrs)[j].first != 0) {
               if (mpiGrid[(*nbrs)[j].first]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  mpiGrid[cells[i]]->sysBoundaryLayer = 1;
               }
            }
         }
      }
   }

   /*communicate which cells have Layer 1 set above for local cells (sysBoundaryFlag
    * and sysBoundaryLayer communicated)*/
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);

   mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_NEIGHBORHOOD_ID);

   /*Compute distances*/
   uint maxLayers = 3; // max(max(P::xcells_ini, P::ycells_ini), P::zcells_ini);
   for (uint layer = 1; layer < maxLayers; layer++) {
      for (uint i = 0; i < cells.size(); i++) {
         if (mpiGrid[cells[i]]->sysBoundaryLayer == 0) {
            const auto* nbrs = mpiGrid.get_neighbors_of(cells[i], SYSBOUNDARIES_NEIGHBORHOOD_ID);
            // Note: this distance calculation will be non-plateau monotonic only assuming that
            // SysBoundary::checkRefinement has been applied correctly and there are no refinement
            // level changes within SYSBOUNDARIES_NEIGHBORHOOD_ID.
            for (uint j = 0; j < (*nbrs).size(); j++) {
               if ((*nbrs)[j].first != 0 && (*nbrs)[j].first != cells[i]) {
                  if (mpiGrid[(*nbrs)[j].first]->sysBoundaryLayer == layer) {
                     mpiGrid[cells[i]]->sysBoundaryLayer = layer + 1;
                     break;
                  }
               }
            }
         }
      }

      SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
      mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_NEIGHBORHOOD_ID);
   }

   /*set cells to DO_NOT_COMPUTE if they are on boundary, and are not
    * in the first two layers of the boundary*/
   for (uint i = 0; i < cells.size(); i++) {
      if (mpiGrid[cells[i]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
          mpiGrid[cells[i]]->sysBoundaryLayer != 1 && mpiGrid[cells[i]]->sysBoundaryLayer != 2) {
         mpiGrid[cells[i]]->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
      }
   }

   // The following is done so that everyone knows their neighbour's layer
   // flags.  This is needed for the correct use of the boundary local
   // communication patterns.
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
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
               if ((layer == 1 && technicalGrid.get(x, y, z)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) ||
                   (layer > 1 && technicalGrid.get(x, y, z)->sysBoundaryLayer == 0)) {

                  if (belongsToLayer(layer, x, y, z, technicalGrid)) {

                     technicalGrid.get(x, y, z)->sysBoundaryLayer = layer;

                     if (layer > 2 && technicalGrid.get(x, y, z)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
                        technicalGrid.get(x, y, z)->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
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
            if (technicalGrid.get(x, y, z)->sysBoundaryLayer == 0 &&
                technicalGrid.get(x, y, z)->sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
               technicalGrid.get(x, y, z)->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
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
            if (technicalGrid.get(x, y, z)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BX;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BY;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BZ;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EX;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EY;
               technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EZ;
            } else {
               if (technicalGrid.get(x - 1, y, z)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BX;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EY;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EZ;
               }
               if (technicalGrid.get(x, y - 1, z)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BY;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EX;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EZ;
               }
               if (technicalGrid.get(x, y, z - 1)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::BZ;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EX;
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EY;
               }
               if (technicalGrid.get(x - 1, y - 1, z)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EZ;
               }
               if (technicalGrid.get(x - 1, y, z - 1)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  technicalGrid.get(x, y, z)->SOLVE = technicalGrid.get(x, y, z)->SOLVE | compute::EY;
               }
               if (technicalGrid.get(x, y - 1, z - 1)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
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
void SysBoundary::applyInitialState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                    Project& project) {
   list<SBC::SysBoundaryCondition*>::iterator it;
   for (it = sysBoundaries.begin(); it != sysBoundaries.end(); it++) {
      // Skip when not restarting or not requested.
      if (Parameters::isRestart && !(*it)->doApplyUponRestart()) {
         continue;
      }
      (*it)->applyInitialState(mpiGrid, perBGrid, project);
   }
}

void SysBoundary::updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                              FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                              creal t) {
   list<SBC::SysBoundaryCondition*>::iterator it;
   if (anyDynamic) {
      for (it = sysBoundaries.begin(); it != sysBoundaries.end(); it++) {
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
void SysBoundary::applyBoundaryVlasovConditions(
    dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, creal t,
    const bool doCalcMomentsV) {

   if (sysBoundaries.size() == 0) {
      return; // no system boundaries
   }

      /*Transfer along boundaries*/
      // First the small stuff without overlapping in an extended neighbourhood:
#warning TODO This now communicates in the wider neighbourhood for both layers, could be reduced to smaller neighbourhood for layer 1, larger neighbourhood for layer 2.
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_PARAMETERS | Transfer::POP_METADATA | Transfer::CELL_SYSBOUNDARYFLAG, true);
   mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);

   // Loop over existing particle species
   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
      SpatialCell::setCommunicatedSpecies(popID);
      // update lists in larger neighborhood
      updateRemoteVelocityBlockLists(mpiGrid, popID, SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);

      // Then the block data in the reduced neighbourhood:
      int timer = phiprof::initializeTimer("Start comm of cell and block data", "MPI");
      phiprof::start(timer);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA, true);
      mpiGrid.start_remote_neighbor_copy_updates(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
      phiprof::stop(timer);

      timer = phiprof::initializeTimer("Compute process inner cells");
      phiprof::start(timer);

      // Compute Vlasov boundary condition on boundary/process inner cells
      vector<CellID> localCells;
      getBoundaryCellList(mpiGrid, mpiGrid.get_local_cells_not_on_process_boundary(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID), localCells);

#pragma omp parallel for
      for (uint i = 0; i < localCells.size(); i++) {
         cuint sysBoundaryType = mpiGrid[localCells[i]]->sysBoundaryFlag;
         this->getBoundary(sysBoundaryType)->vlasovBoundaryCondition(mpiGrid, localCells[i], popID, doCalcMomentsV);
      }
      if (doCalcMomentsV)
         calculateMoments_V(mpiGrid, localCells, true);
      else
         calculateMoments_R(mpiGrid, localCells, true);

      phiprof::stop(timer);

      timer = phiprof::initializeTimer("Wait for receives", "MPI", "Wait");
      phiprof::start(timer);
      mpiGrid.wait_remote_neighbor_copy_updates(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
      phiprof::stop(timer);

      // Compute vlasov boundary on boundary/process boundary cells
      timer = phiprof::initializeTimer("Compute process boundary cells");
      phiprof::start(timer);
      vector<CellID> boundaryCells;
      getBoundaryCellList(mpiGrid, mpiGrid.get_local_cells_on_process_boundary(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID),
                          boundaryCells);
#pragma omp parallel for
      for (uint i = 0; i < boundaryCells.size(); i++) {
         cuint sysBoundaryType = mpiGrid[boundaryCells[i]]->sysBoundaryFlag;
         this->getBoundary(sysBoundaryType)->vlasovBoundaryCondition(mpiGrid, boundaryCells[i], popID, doCalcMomentsV);
      }
      if (doCalcMomentsV)
         calculateMoments_V(mpiGrid, boundaryCells, true);
      else
         calculateMoments_R(mpiGrid, boundaryCells, true);

      phiprof::stop(timer);

      // WARNING Blocks are changed but lists not updated now, if you need to
      // use/communicate them before the next update is done, add an update
      // here. Reset lists in smaller default neighborhood.
      updateRemoteVelocityBlockLists(mpiGrid, popID);

   } // for-loop over populations
}

/*! Get a pointer to the SysBoundaryCondition of given index.
 * \param sysBoundaryType Type of the boundary condition to return
 * \return Pointer to the instance of the SysBoundaryCondition
 */
SBC::SysBoundaryCondition* SysBoundary::getBoundary(cuint sysBoundaryType) const {
   auto it = indexToBoundary.find(sysBoundaryType);
   if (it != indexToBoundary.end())
      return it->second;
   else
      SBC::abort_mpi("ERROR: Boundary " + to_string(sysBoundaryType) + " is invalid", 1);
}

/*! Get the number of SysBoundaryConditions stored in SysBoundary. */
unsigned int SysBoundary::size() const { return sysBoundaries.size(); }

/*! Check whether any boundary condition is dynamic in time. */
bool SysBoundary::isDynamic() const { return anyDynamic; }

/*! Get a bool telling whether the system is periodic in the queried direction.
 * \param direction 0: x, 1: y, 2: z.
 * \retval periodic Is the system periodic in the queried direction.
 */
bool SysBoundary::isPeriodic(uint direction) const { return periodic[direction]; }

/*! Get a vector containing the cellID of all cells which are not DO_NOT_COMPUTE or
 * NOT_SYSBOUNDARY in the vector of cellIDs passed to the function.
 * \param mpiGrid Grid
 * \param cellList Vector of cellIDs in which to look for boundary cells
 * \param boundaryCellList Vector of boundary the cells' cellIDs
 */
void getBoundaryCellList(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                         const vector<uint64_t>& cellList, vector<uint64_t>& boundaryCellList) {
   boundaryCellList.clear();
   for (size_t cell = 0; cell < cellList.size(); ++cell) {
      const CellID cellID = cellList[cell];
      if (mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
          mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         continue;
      }
      boundaryCellList.push_back(cellID);
   }
}

/*! Updates all NonboundaryCells into an internal map. This should be called in
 * loadBalance.
 * \param mpiGrid The DCCRG grid
 */
void SysBoundary::updateBoundariesAfterLoadBalance(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
   vector<uint64_t> local_cells_on_boundary;
   getBoundaryCellList(mpiGrid, mpiGrid.get_cells(), local_cells_on_boundary);
   // Loop over boundaries:
   for (std::list<SBC::SysBoundaryCondition*>::iterator it = sysBoundaries.begin(); it != sysBoundaries.end(); ++it) {
      (*it)->updateBoundaryConditionsAfterLoadBalance(mpiGrid, local_cells_on_boundary);
   }
}
