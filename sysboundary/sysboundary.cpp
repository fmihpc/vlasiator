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
#include "../vlasovsolver/arch_moments.h"

#include "donotcompute.h"
#include "ionosphere.h"
#include "copysphere.h"
#include "outflow.h"
#include "setmaxwellian.h"
#include "sysboundary.h"
#include "../fieldsolver/gridGlue.hpp"

using namespace std;
using namespace spatial_cell;

bool precedenceSort(const SBC::SysBoundaryCondition* first, const SBC::SysBoundaryCondition* second) {
   if (first->getPrecedence() < second->getPrecedence()) {
      return true;
   } else {
      return false;
   }
}

// ************************************************************
// ***** DEFINITIONS FOR BOUNDARY CLASS *****
// ************************************************************

SysBoundary::SysBoundary() : anyDynamic(false) {}

/*!\brief Destructor for class SysBoundary.
 *
 * Reduces the value of SysBoundary::nSysBoundaries by one,
 * and if after the destruction SysBoundary::nSysBoundaries equals zero all stored SysBoundaries are deleted.
 */
SysBoundary::~SysBoundary() {
   // Call delete for each SysBoundaryCondition:
   for (list<SBC::SysBoundaryCondition*>::iterator it = sysBoundaries.begin(); it != sysBoundaries.end(); it++) {
      delete *it;
      *it = NULL;
   }
}

/*!\brief Add its own and all existing SysBoundaryConditions' parameters.
 *
 * Adds the parameters specific to the SysBondary class handling the list of
 * SysBoundaryConditions and then calls the static addParameters functions of all
 * SysBoundaryConditions implemented in the code in order to have them appear also in the
 * help.
 */
void SysBoundary::addParameters() {
   Readparameters::addComposing(
       "boundaries.boundary",
       "List of boundary condition (BC) types to be used. Each boundary condition to be used has to be on a new line "
       "boundary = YYY. Available options are: Outflow, Ionosphere, Copysphere, Maxwellian.");
   Readparameters::add("boundaries.periodic_x", "Set the grid periodicity in x-direction. 'yes'(default)/'no'.", "yes");
   Readparameters::add("boundaries.periodic_y", "Set the grid periodicity in y-direction. 'yes'(default)/'no'.", "yes");
   Readparameters::add("boundaries.periodic_z", "Set the grid periodicity in z-direction. 'yes'(default)/'no'.", "yes");

   // call static addParameter functions in all bc's
   SBC::DoNotCompute::addParameters();
   SBC::Ionosphere::addParameters();
   SBC::Copysphere::addParameters();
   SBC::Outflow::addParameters();
   SBC::Maxwellian::addParameters();
}

/*!\brief Get this class' parameters.
 *
 * Get the parameters pertaining to this class.
 *
 * getParameters for each actually used system boundary condition is called by each
 * SysBoundaryCondition's initialization function.
 */
void SysBoundary::getParameters() {
   string periodic_x, periodic_y, periodic_z;

   Readparameters::get("boundaries.boundary", sysBoundaryCondList);
   Readparameters::get("boundaries.periodic_x", periodic_x);
   Readparameters::get("boundaries.periodic_y", periodic_y);
   Readparameters::get("boundaries.periodic_z", periodic_z);

   periodic[0] = (periodic_x == "yes");
   periodic[1] = (periodic_y == "yes");
   periodic[2] = (periodic_z == "yes");
}

/*! Add a new SBC::SysBoundaryCondition which has been created with new sysBoundary.
 * SysBoundary will take care of deleting it.
 *
 * \param bc SysBoundaryCondition object
 * \param project Project object
 * \param t Current time
 * \retval success If true, the given SBC::SysBoundaryCondition was added successfully.
 */
void SysBoundary::addSysBoundary(SBC::SysBoundaryCondition* bc, Project& project, creal& t) {
   // Initialize the boundary condition
   stringstream timername;
   timername<<"Initialize system boundary condition "<<bc->getName();
   phiprof::Timer timer {timername.str()};
   bc->initSysBoundary(t, project);
   timer.stop();

   sysBoundaries.push_back(bc);
   if (sysBoundaries.size() > 1) {
      sysBoundaries.sort(precedenceSort);
   }

   // This assumes that only one instance of each type is created.
   indexToSysBoundary[bc->getIndex()] = bc;
}

/*!\brief Initialise all system boundary conditions actually used.
 *
 * This function loops through the list of system boundary conditions listed as to be used
 * in the configuration file/command line arguments. For each of these it adds the
 * corresponding instance and updates the member isDynamic to determine whether any
 * SysBoundaryCondition is dynamic in time.
 *
 * \param project Project object
 * \param t Current time
 * \retval success If true, the initialisation of all system boundaries succeeded.
 * \sa addSysBoundary
 */
void SysBoundary::initSysBoundaries(Project& project, creal& t) {
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   vector<string>::const_iterator it;

   if (sysBoundaryCondList.size() == 0) {
      if (!periodic[0] && !Readparameters::helpRequested) {
         abort_mpi("Non-periodic in x but no boundary condtion loaded!");
      }
      if (!periodic[1] && !Readparameters::helpRequested) {
         abort_mpi("Non-periodic in y but no boundary condtion loaded!");
      }
      if (!periodic[2] && !Readparameters::helpRequested) {
         abort_mpi("Non-periodic in z but no boundary condtion loaded!");
      }
   }

   for (it = sysBoundaryCondList.begin(); it != sysBoundaryCondList.end(); it++) {
      if (*it == "Outflow" || *it == "outflow") {
         this->addSysBoundary(::new SBC::Outflow, project, t);

         anyDynamic = anyDynamic | this->getSysBoundary(sysboundarytype::OUTFLOW)->isDynamic();
         bool faces[6];
         this->getSysBoundary(sysboundarytype::OUTFLOW)->getFaces(&faces[0]);

         if ((faces[0] || faces[1]) && periodic[0]) {
            abort_mpi("Conflict: x boundaries set to periodic but found Outflow conditions!");
         }

         if ((faces[2] || faces[3]) && periodic[1]) {
            abort_mpi("Conflict: y boundaries set to periodic but found Outflow conditions!");
         }

         if ((faces[4] || faces[5]) && periodic[2]) {
            abort_mpi("Conflict: z boundaries set to periodic but found Outflow conditions!");
         }
         if ((faces[0] || faces[1]) && P::xcells_ini < 5) {
            abort_mpi("Outflow condition loaded on x- or x+ face but not enough cells in x!");
         }
         if ((faces[2] || faces[3]) && P::ycells_ini < 5) {
            abort_mpi("Outflow condition loaded on y- or y+ face but not enough cells in y!");
         }

         if ((faces[4] || faces[5]) && P::zcells_ini < 5) {
            abort_mpi("Outflow condition loaded on z- or z+ face but not enough cells in z!");
         }

      } else if (*it == "Ionosphere" || *it == "ionosphere") {
         this->addSysBoundary(::new SBC::Ionosphere, project, t);
         this->addSysBoundary(::new SBC::DoNotCompute, project, t);
         anyDynamic = anyDynamic | this->getSysBoundary(sysboundarytype::IONOSPHERE)->isDynamic();
      } else if(*it == "Copysphere" || *it == "copysphere") {
         this->addSysBoundary(::new SBC::Copysphere, project, t);
         this->addSysBoundary(::new SBC::DoNotCompute, project, t);
         anyDynamic = anyDynamic | this->getSysBoundary(sysboundarytype::COPYSPHERE)->isDynamic();
      } else if (*it == "Maxwellian" || *it == "maxwellian") {
         this->addSysBoundary(::new SBC::Maxwellian, project, t);
         anyDynamic = anyDynamic | this->getSysBoundary(sysboundarytype::MAXWELLIAN)->isDynamic();
         bool faces[6];
         this->getSysBoundary(sysboundarytype::MAXWELLIAN)->getFaces(&faces[0]);
         if ((faces[0] || faces[1]) && periodic[0]) {
            abort_mpi("Conflict: x boundaries set to periodic but found Maxwellian also!");
         }
         if ((faces[2] || faces[3]) && periodic[1]) {
            abort_mpi("Conflict: y boundaries set to periodic but found Maxwellian also!");
         }
         if ((faces[4] || faces[5]) && periodic[2]) {
            abort_mpi("Conflict: z boundaries set to periodic but found Maxwellian also!");
         }
         if ((faces[0] || faces[1]) && P::xcells_ini < 5) {
            abort_mpi("Maxwellian condition loaded on x- or x+ face but not enough cells in x!");
         }
         if ((faces[2] || faces[3]) && P::ycells_ini < 5) {
            abort_mpi("Maxwellian condition loaded on y- or y+ face but not enough cells in y!");
         }
         if ((faces[4] || faces[5]) && P::zcells_ini < 5) {
            abort_mpi("Maxwellian condition loaded on z- or z+ face but not enough cells in z!");
         }
      } else {
         std::ostringstream msg;
         msg << "Unknown type of boundary read: " << *it;
         abort_mpi(msg.str());
      }
   }

   for (auto& b : sysBoundaries)  {
      b->setPeriodicity(periodic);
   }
}

/*!\brief Boolean check if queried sysboundarycondition exists
 * Note: this queries against the parsed list of names, not against the actual existing boundaries.
 * This is so that the project configuration (which is handled before grid and boundary initialization)
 * can use this call.
 *
 * \param string name
 * \retval success if queried boundary exists
 * \sa addSysBoundary
 */
bool SysBoundary::existSysBoundary(std::string name) {
   vector<string>::const_iterator it;
   for (it = sysBoundaryCondList.begin(); it != sysBoundaryCondList.end(); it++) {
      if ((*it) == name) {
         return true;
      }
   }
   return false;
}

void SysBoundary::checkRefinement(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
   // Verifies that all cells within Neighborhoods::FULL of L1 boundary cells are on the same refinement
   // level (one group for inner boundary, another for outer boundary)

   // Set is used to avoid storing duplicates - each cell only needs to be checked once
   set<CellID> innerBoundaryCells;
   set<CellID> outerBoundaryCells;

   int innerBoundaryRefLvl = -1;
   int outerBoundaryRefLvl = -1;

   const vector<CellID>& local_cells = getLocalCells();
   // Collect cells by sysboundarytype
   for (auto cellId : local_cells) {
      SpatialCell* cell = mpiGrid[cellId];
      if (cell) {
         if (cell->sysBoundaryFlag == sysboundarytype::IONOSPHERE || cell->sysBoundaryFlag == sysboundarytype::COPYSPHERE) {
            innerBoundaryCells.insert(cellId);
            innerBoundaryRefLvl = mpiGrid.get_refinement_level(cellId);
            if (cell->sysBoundaryLayer == 1) {
               // Add all stencil neighbors of layer 1 cells
               auto* nbrPairVector = mpiGrid.get_neighbors_of(cellId, Neighborhoods::SYSBOUNDARIES);
               for (auto nbrPair : *nbrPairVector) {
                  if (nbrPair.first != INVALID_CELLID) {
                     innerBoundaryCells.insert(nbrPair.first);
                  }
               }
            }
         } else if (cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
                    cell->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE) {
            outerBoundaryCells.insert(cellId);
            outerBoundaryRefLvl = mpiGrid.get_refinement_level(cellId);
            // Add all stencil neighbors of outer boundary cells
            auto* nbrPairVector = mpiGrid.get_neighbors_of(cellId, Neighborhoods::SYSBOUNDARIES);
            for (auto nbrPair : *nbrPairVector) {
               if (nbrPair.first != INVALID_CELLID) {
                  outerBoundaryCells.insert(nbrPair.first);
               }
            }
         }
      }
   }

   for (auto cellId : innerBoundaryCells) {
      if (cellId != INVALID_CELLID && mpiGrid.get_refinement_level(cellId) != innerBoundaryRefLvl) {
         abort_mpi("ERROR: inner boundary cells must have identical refinement level!");
      }
   }

   for (auto cellId : outerBoundaryCells) {
      if (cellId != INVALID_CELLID && mpiGrid.get_refinement_level(cellId) != outerBoundaryRefLvl) {
         abort_mpi("ERROR: outer boundary cells must have identical refinement level!");
      }
   }
}

bool belongsToLayer(const int layer, const int x, const int y, const int z,
                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {

   bool belongs = false;

   // loop through all neighbors (including diagonals)
   for (int iz = -1; iz <= 1; ++iz) {
      for (int iy = -1; iy <= 1; ++iy) {
         for (int ix = -1; ix <= 1; ++ix) {

            // not strictly necessary but logically we should not consider the cell itself
            // among its neighbors.
            if ((ix == 0 && iy == 0 && iz == 0) || !technicalGrid.get(x + ix, y + iy, z + iz)) {
               continue;
            }

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

/*!\brief Classify all simulation cells with respect to the system boundary conditions.
 *
 * Loops through all cells and and for each assigns the correct sysBoundaryFlag depending on
 * the return value of each SysBoundaryCondition's assignSysBoundary.
 *
 * \param mpiGrid Grid
 */
void SysBoundary::classifyCells(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                 FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
                                 const vector<CellID>& cells = getLocalCells();
                                 auto localSize = technicalGrid.getLocalSize().data();

   /*set all cells to default value, not_sysboundary */
#pragma omp parallel for
   for (uint i = 0; i < cells.size(); i++) {
      mpiGrid[cells[i]]->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
   }
#pragma omp parallel for collapse(2)
   for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
      for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            //technicalGrid.get(x, y, z)->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
            // Here for debugging since boundarytype should be fed from MPIGrid
            technicalGrid.get(x, y, z)->sysBoundaryFlag = sysboundarytype::N_SYSBOUNDARY_CONDITIONS;
            technicalGrid.get(x, y, z)->sysBoundaryLayer = 0;
            // Function called on every refinement, we only want to reset dt on simulation start
            if (P::tstep == P::tstep_min) {
               technicalGrid.get(x, y, z)->maxFsDt = numeric_limits<Real>::max();
            }
            // Set the fsgrid rank in the technical grid
            technicalGrid.get(x, y, z)->fsGridRank = technicalGrid.getRank();
         }
      }
   }

   /*
     loop through sysboundaries and let all sysboundaries set in local
   cells if they are part of which sysboundary (cell location needs to
   be updated by now. No remote data needed/available, so assignement
   has to be based individually on each cells location
   */
   for(auto& b : sysBoundaries) {
      b->assignSysBoundary(mpiGrid, technicalGrid);
   }

   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
   mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::SYSBOUNDARIES);

   feedBoundaryIntoFsGrid(mpiGrid, cells, technicalGrid);

   // set distance 1 cells to boundary cells, that have neighbors which are normal cells
   for (CellID cell : cells) {
      mpiGrid[cell]->sysBoundaryLayer = 0; /*Initial value*/

      std::array<double, 3> dx = mpiGrid.geometry.get_length(cell);
      std::array<double, 3> x = mpiGrid.get_center(cell);
      if (!isPeriodic(0) && (x[0] > Parameters::xmax - dx[0] || x[0] < Parameters::xmin + dx[0])) {
         continue;
      } else if (!isPeriodic(1) && (x[1] > Parameters::ymax - dx[1] || x[1] < Parameters::ymin + dx[1])) {
         continue;
      } else if (!isPeriodic(2) && (x[2] > Parameters::zmax - dx[2] || x[2] < Parameters::zmin + dx[2])) {
         continue;
      }

      if (mpiGrid[cell]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
         // Cornerwise neighbor, i.e. cell must be in both neighbors_of and neighbors_to
         for (auto i : *mpiGrid.get_neighbors_of(cell, Neighborhoods::SYSBOUNDARIES)) {
            CellID neighbor = i.first;
            if (neighbor && mpiGrid[neighbor]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
               for (auto j : *mpiGrid.get_neighbors_to(cell, Neighborhoods::SYSBOUNDARIES)) {
                  if (j.first == neighbor) {
                     mpiGrid[cell]->sysBoundaryLayer = 1;
                  }
               }
            }
         }
      }
   }

   /*communicate which cells have Layer 1 set above for local cells (sysBoundaryFlag
    * and sysBoundaryLayer communicated)*/
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);

   mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::SYSBOUNDARIES);

   /*Compute distances*/
   uint maxLayers = 3; // max(max(P::xcells_ini, P::ycells_ini), P::zcells_ini);
   for (uint layer = 1; layer < maxLayers; layer++) {
      for (CellID cell : cells) {
         if (mpiGrid[cell]->sysBoundaryLayer == 0) {
            // Note: this distance calculation will be non-plateau monotonic only assuming that
            // SysBoundary::checkRefinement has been applied correctly and there are no refinement
            // level changes within Neighborhoods::SYSBOUNDARIES.
            for (auto i : *mpiGrid.get_neighbors_of(cell, Neighborhoods::SYSBOUNDARIES)) {
               CellID neighbor = i.first;
               if (neighbor && mpiGrid[neighbor]->sysBoundaryLayer == layer) {
                  for (auto j : *mpiGrid.get_neighbors_to(cell, Neighborhoods::SYSBOUNDARIES)) {
                     if (j.first == neighbor) {
                        mpiGrid[cell]->sysBoundaryLayer = layer + 1;
                     }
                  }
               }
            }
         }
      }

      SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
      mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::SYSBOUNDARIES);
   }

   /*set cells to DO_NOT_COMPUTE if they are on boundary, and are not
    * in the first two layers of the boundary*/
   for (uint i = 0; i < cells.size(); i++) {
      if (mpiGrid[cells[i]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
          mpiGrid[cells[i]]->sysBoundaryLayer != 1 && mpiGrid[cells[i]]->sysBoundaryLayer != 2) {
         mpiGrid[cells[i]]->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
      }
   }

   // The following is done so that everyone knows their neighbour's
   // layer flags.  This is needed for the correct use of the system
   // boundary local communication patterns.
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
   mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::FULL);

   SysBoundary& sysBoundaryContainer = getObjectWrapper().sysBoundaryContainer;
   Real ionosphereDownmapRadius = 0;
   if (sysBoundaryContainer.existSysBoundary("Ionosphere")) {
      Readparameters::get("ionosphere.downmapRadius", ionosphereDownmapRadius);
   }
   if(ionosphereDownmapRadius < 1000) {
      ionosphereDownmapRadius *= physicalconstants::R_E;
   }

   // Now the layers need to be set on fsgrid too
   // In dccrg initialization the max number of boundary layers is set to 3.
   const uint MAX_NUMBER_OF_BOUNDARY_LAYERS = 3 * pow(2, mpiGrid.get_maximum_refinement_level());

   technicalGrid.updateGhostCells();

   // loop through max number of layers
   for (uint layer = 1; layer <= MAX_NUMBER_OF_BOUNDARY_LAYERS; ++layer) {

// loop through all cells in grid
#pragma omp parallel for collapse(2)
      for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
         for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
            for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {

               // for the first layer, consider all cells that belong to a boundary, for other layers
               // consider all cells that have not yet been labeled.
               if ((layer == 1 && technicalGrid.get(x, y, z)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) ||
                   (layer > 1 && technicalGrid.get(x, y, z)->sysBoundaryLayer == 0)) {

                  if (belongsToLayer(layer, x, y, z, technicalGrid)) {

                     technicalGrid.get(x, y, z)->sysBoundaryLayer = layer;

                     if (layer > 2 && (technicalGrid.get(x,y,z)->sysBoundaryFlag == sysboundarytype::IONOSPHERE || 
                                       technicalGrid.get(x,y,z)->sysBoundaryFlag == sysboundarytype::COPYSPHERE)) {
                        technicalGrid.get(x, y, z)->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
                     } else if (layer > 2 && technicalGrid.get(x, y, z)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
                        technicalGrid.get(x, y, z)->sysBoundaryFlag = sysboundarytype::OUTER_BOUNDARY_PADDING;
                     }
                  }
               }
            }
         }
      }
      // This needs an update every iteration as belongsToLayer() needs up to date data.
      technicalGrid.updateGhostCells();
   }

// One more pass to make sure, in particular if the ionosphere is wide enough
// there is remaining cells of IONOSPHERE type inside the max layers gone through previously.
// This last pass now gets rid of them.
#pragma omp parallel for collapse(2)
   for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
      for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            if (technicalGrid.get(x,y,z)->sysBoundaryLayer == 0 && (
                technicalGrid.get(x,y,z)->sysBoundaryFlag == sysboundarytype::IONOSPHERE ||
                technicalGrid.get(x,y,z)->sysBoundaryFlag == sysboundarytype::COPYSPHERE)) {
               technicalGrid.get(x, y, z)->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
            }
         }
      }
   }

   technicalGrid.updateGhostCells();

   const array<FsGridTools::FsSize_t,3> fsGridDimensions = technicalGrid.getGlobalSize();

   // One pass to setup the bit field to know which components the field solver should propagate.
#pragma omp parallel for collapse(2)
   for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
      for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            technicalGrid.get(x, y, z)->SOLVE = 0;

            array<FsGridTools::FsIndex_t, 3> globalIndices = technicalGrid.getGlobalIndices(x, y, z);

            if (((globalIndices[0] == 0 || globalIndices[0] == (FsGridTools::FsIndex_t)fsGridDimensions[0] - 1) &&
                 !this->isPeriodic(0)) ||
                ((globalIndices[1] == 0 || globalIndices[1] == (FsGridTools::FsIndex_t)fsGridDimensions[1] - 1) &&
                 !this->isPeriodic(1)) ||
                ((globalIndices[2] == 0 || globalIndices[2] == (FsGridTools::FsIndex_t)fsGridDimensions[2] - 1) &&
                 !this->isPeriodic(2))) {
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

/*!\brief Apply the initial state to all system boundary cells.
 * Loops through all SysBoundaryConditions and calls the corresponding applyInitialState
 * function. This function must apply the initial state for all existing particle species.
 *
 * \param mpiGrid Grid
 * \retval success If true, the application of all system boundary states succeeded.
 */
void SysBoundary::applyInitialState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>&technicalGrid,
                                    FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                    FsGrid<array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                    Project& project) {

   list<SBC::SysBoundaryCondition*>::iterator it;
   for (it = sysBoundaries.begin(); it != sysBoundaries.end(); it++) {
      if (                                // This is to skip the reapplication
          Parameters::isRestart           // When not restarting
          && !(*it)->doApplyUponRestart() // When reapplicaiton is not requested
      ) {
         continue;
      }
      stringstream timername;
      timername<<"Apply system boundary condition "<<(*it)->getName()<<" initial state";
      phiprof::Timer timer {timername.str()};
      (*it)->applyInitialState(mpiGrid, technicalGrid, perBGrid, BgBGrid, project);
   }
}

void SysBoundary::updateState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                              FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
                              FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                              FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                              creal t) {
   if (isAnyDynamic()) {
      for(auto& b : sysBoundaries) {
         if (b->isDynamic()) {
            b->updateState(mpiGrid, technicalGrid, perBGrid, BgBGrid, t);
         }
      }
   }
}


void SysBoundary::setupL2OutflowAtRestart(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
   list<SBC::SysBoundaryCondition*>::iterator it;
   for (it = sysBoundaries.begin(); it != sysBoundaries.end(); it++) {
      if (Parameters::isRestart           // When restarting
         && !(*it)->doApplyUponRestart() // When reapplication is not requested
         && (*it)->getIndex() == sysboundarytype::OUTFLOW
      ) {
         (*it)->setupL2OutflowAtRestart(mpiGrid);
      }
   }
   // Needed after copying VDFs in L1 Outflow cells or LUMI (and our communications) breaks.
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      updateRemoteVelocityBlockLists(mpiGrid, popID);
   }
}


/*!\brief Apply the Vlasov system boundary conditions to all system boundary cells at time t.
 *
 * Loops through all SysBoundaryConditions and calls the corresponding vlasovBoundaryCondition() function.
 * The boundary condition functions are called for all particle species, one at a time.
 *
 * WARNING (see end of the function) Blocks are changed but lists not updated now,
 * if you need to use/communicate them before the next update is done, add an
 * update at the end of this function.
 *
 * \param mpiGrid Grid
 * \param t Current time
 * \param calculate_V_moments if true, compute into _V, false into _R moments so that the interpolated ones can be done
 */
void SysBoundary::applySysBoundaryVlasovConditions(
    dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, creal& t,
    const bool calculate_V_moments) {

   if (sysBoundaries.size() == 0) {
      return; // no system boundaries
   }

/*Transfer along boundaries*/
// First the small stuff without overlapping in an extended neighbourhood:
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_PARAMETERS | Transfer::POP_METADATA | Transfer::CELL_SYSBOUNDARYFLAG, true);
   mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::SYSBOUNDARIES_EXTENDED);

   // Loop over existing particle species
   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
      SpatialCell::setCommunicatedSpecies(popID);
      // update lists in neighborhood
      updateRemoteVelocityBlockLists(mpiGrid, popID, Neighborhoods::SYSBOUNDARIES);

      // Then the block data in the reduced neighbourhood:
      phiprof::Timer commTimer {"Start comm of cell and block data", {"MPI"}};
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA, true);
      mpiGrid.start_remote_neighbor_copy_updates(Neighborhoods::SYSBOUNDARIES);
      commTimer.stop();

      phiprof::Timer computeInnerTimer {"Compute process inner cells"};
      // Compute Vlasov boundary condition on system boundary/process inner cells
      vector<CellID> localCells;
      getBoundaryCellList(mpiGrid, mpiGrid.get_local_cells_not_on_process_boundary(Neighborhoods::SYSBOUNDARIES), localCells);

#pragma omp parallel for
      for (uint i = 0; i < localCells.size(); i++) {
         cuint sysBoundaryType = mpiGrid[localCells[i]]->sysBoundaryFlag;
         this->getSysBoundary(sysBoundaryType)->vlasovBoundaryCondition(mpiGrid, localCells[i], popID, calculate_V_moments);
      }
      if (popID==getObjectWrapper().particleSpecies.size()-1) {
         // Only calculate moments when handling last population
         if (calculate_V_moments) {
            calculateMoments_V(mpiGrid, localCells, true);
         } else {
            calculateMoments_R(mpiGrid, localCells, true);
         }
      }
      computeInnerTimer.stop();

      phiprof::Timer waitimer {"Wait for receives", {"MPI", "Wait"}};
      mpiGrid.wait_remote_neighbor_copy_updates(Neighborhoods::SYSBOUNDARIES);
      waitimer.stop();

      // Compute vlasov boundary on system boundary/process boundary cells
      phiprof::Timer computeBoundaryTimer {"Compute process boundary cells"};
      vector<CellID> boundaryCells;
      getBoundaryCellList(mpiGrid, mpiGrid.get_local_cells_on_process_boundary(Neighborhoods::SYSBOUNDARIES),
                          boundaryCells);
#pragma omp parallel for
      for (uint i = 0; i < boundaryCells.size(); i++) {
         cuint sysBoundaryType = mpiGrid[boundaryCells[i]]->sysBoundaryFlag;
         this->getSysBoundary(sysBoundaryType)->vlasovBoundaryCondition(mpiGrid, boundaryCells[i], popID, calculate_V_moments);
      }
      if (popID==getObjectWrapper().particleSpecies.size()-1) {
         // Only calculate moments when handling last population
         if (calculate_V_moments) {
            calculateMoments_V(mpiGrid, boundaryCells, true);
         } else {
            calculateMoments_R(mpiGrid, boundaryCells, true);
         }
      }
      computeBoundaryTimer.stop();

      // WARNING Blocks are changed but lists not updated now, if you need to use/communicate them before the next
      // update is done, add an update here. reset lists in smaller default neighborhood
      updateRemoteVelocityBlockLists(mpiGrid, popID);

   } // for-loop over populations
}

/*! Get a pointer to the SysBoundaryCondition of given index.
 * \param sysBoundaryType Type of the system boundary condition to return
 * \return Pointer to the instance of the SysBoundaryCondition. NULL if sysBoundaryType is invalid.
 */
SBC::SysBoundaryCondition* SysBoundary::getSysBoundary(cuint sysBoundaryType) const {
   auto it = indexToSysBoundary.find(sysBoundaryType);
   if (it != indexToSysBoundary.end()) {
      return it->second;
   } else {
      abort_mpi("ERROR: Boundary " + to_string(sysBoundaryType) + " is invalid", 1);
   }
}

/*! Get the number of SysBoundaryConditions stored in SysBoundary.
 * \retval size Number of SysBoundaryConditions stored in SysBoundary.
 */
unsigned int SysBoundary::size() const { return sysBoundaries.size(); }

/*! Get a bool telling whether any system boundary condition is dynamic in time (and thus needs updating).
 * \retval isDynamic Is any system boundary condition dynamic in time.
 */
bool SysBoundary::isAnyDynamic() const { return anyDynamic; }

/*! Get a bool telling whether the system is periodic in the queried direction.
 * \param direction 0: x, 1: y, 2: z.
 * \retval periodic Is the system periodic in the queried direction.
 */
bool SysBoundary::isPeriodic(uint direction) const { return periodic[direction]; }

/*! Get a vector containing the cellID of all cells which are not DO_NOT_COMPUTE or NOT_SYSBOUNDARY in the vector of
 * cellIDs passed to the function.
 *
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

/*! Updates all NonsysboundaryCells into an internal map. This should be called in loadBalance.
 * \param mpiGrid The DCCRG grid
 */
void SysBoundary::updateSysBoundariesAfterLoadBalance(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
   phiprof::Timer timer {"updateSysBoundariesAfterLoadBalance"};
   vector<uint64_t> local_cells_on_boundary;
   getBoundaryCellList(mpiGrid, mpiGrid.get_cells(), local_cells_on_boundary);
   // Loop over sysboundaries:
   for (list<SBC::SysBoundaryCondition*>::iterator it = sysBoundaries.begin(); it != sysBoundaries.end(); ++it) {
      (*it)->updateSysBoundaryConditionsAfterLoadBalance(mpiGrid, local_cells_on_boundary);
   }
}
