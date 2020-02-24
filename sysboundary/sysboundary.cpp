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

#include "sysboundary.h"
#include "donotcompute.h"
#include "ionosphere.h"
#include "outflow.h"
#include "setmaxwellian.h"

using namespace std;
using namespace spatial_cell;

bool precedenceSort(const SBC::SysBoundaryCondition* first,
                    const SBC::SysBoundaryCondition* second) {
   if(first->getPrecedence() < second->getPrecedence()) return true;
   else return false;
}

// ************************************************************
// ***** DEFINITIONS FOR BOUNDARY CLASS *****
// ************************************************************

/*! Constructor for class SysBoundary.*/
SysBoundary::SysBoundary() { }

/*!\brief Destructor for class SysBoundary.
 * 
 * Reduces the value of SysBoundary::nSysBoundaries by one,
 * and if after the destruction SysBoundary::nSysBoundaries equals zero all stored SysBoundaries are deleted.
 */
SysBoundary::~SysBoundary() {
   // Call delete for each SysBoundaryCondition:
   for (list<SBC::SysBoundaryCondition*>::iterator it=sysBoundaries.begin();
        it!=sysBoundaries.end();
        it++) {
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
   Readparameters::addComposing("boundaries.boundary", "List of boundary condition (BC) types to be used. Each boundary condition to be used has to be on a new line boundary = YYY. Available (20140113) are Outflow Ionosphere Maxwellian.");
   Readparameters::add("boundaries.periodic_x","If 'yes' the grid is periodic in x-direction. Defaults to 'no'.","no");
   Readparameters::add("boundaries.periodic_y","If 'yes' the grid is periodic in y-direction. Defaults to 'no'.","no");
   Readparameters::add("boundaries.periodic_z","If 'yes' the grid is periodic in z-direction. Defaults to 'no'.","no");
   
   //call static addParameter functions in all bc's
   SBC::DoNotCompute::addParameters();
   SBC::Ionosphere::addParameters();
   SBC::Outflow::addParameters();
   SBC::SetMaxwellian::addParameters();
}

/*!\brief Get this class' parameters.
 * 
 * Get the parameters pertaining to this class.
 * 
 * getParameters for each actually used system boundary condition is called by each
 * SysBoundaryCondition's initialization function.
 */
void SysBoundary::getParameters() {
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   if(!Readparameters::get("boundaries.boundary", sysBoundaryCondList)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   std::string periodic_x,periodic_y,periodic_z;
   if(!Readparameters::get("boundaries.periodic_x",periodic_x)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!Readparameters::get("boundaries.periodic_y",periodic_y)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   };
   if(!Readparameters::get("boundaries.periodic_z",periodic_z)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   };
   isPeriodic[0] = false;
   isPeriodic[1] = false;
   isPeriodic[2] = false;
   if (periodic_x == "yes") isPeriodic[0] = true;
   if (periodic_y == "yes") isPeriodic[1] = true;
   if (periodic_z == "yes") isPeriodic[2] = true;
}

/*! Add a new SBC::SysBoundaryCondition which has been created with new sysBoundary. 
 * SysBoundary will take care of deleting it.
 * 
 * \param bc SysBoundaryCondition object
 * \param project Project object
 * \param t Current time
 * \retval success If true, the given SBC::SysBoundaryCondition was added successfully.
 */
bool SysBoundary::addSysBoundary(
   SBC::SysBoundaryCondition* bc,
   Project &project,
   creal& t
) {
   // Initialize the boundary condition
   bool success = true;
   if (bc->initSysBoundary(t, project) == false) {
      cerr << "Failed to initialize system boundary condition '" << bc->getName() << "'" << endl;
      return false;
   }

   sysBoundaries.push_back(bc);
   if (sysBoundaries.size() > 1) {
      sysBoundaries.sort(precedenceSort);
   }

   // This assumes that only one instance of each type is created.
   indexToSysBoundary[bc->getIndex()] = bc;

   return success;
}

/*!\brief Initialise all system boundary conditions actually used.
 * 
 * This function loops through the list of system boundary conditions listed as to be used
 * in the configuration file/command line arguments. For each of these it adds the
 * corresponding instance and updates the member isThisDynamic to determine whether any
 * SysBoundaryCondition is dynamic in time.
 * 
 * \param project Project object
 * \param t Current time
 * \retval success If true, the initialisation of all system boundaries succeeded.
 * \sa addSysBoundary
 */
bool SysBoundary::initSysBoundaries(
                                    Project& project,
                                    creal& t
                                   ) {
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   bool success = true;
   vector<string>::const_iterator it;
   
   if (sysBoundaryCondList.size() == 0) {
      if(!isPeriodic[0] && !Readparameters::helpRequested) {
         if(myRank == MASTER_RANK) cerr << "You set boundaries.periodic_x = no but you didn't load any system boundary condition using the option boundaries.boundary, are you sure this is correct?" << endl;
      }
      if(!isPeriodic[1] && !Readparameters::helpRequested) {
         if(myRank == MASTER_RANK) cerr << "You set boundaries.periodic_y = no but you didn't load any system boundary condition using the option boundaries.boundary, are you sure this is correct?" << endl;
      }
      if(!isPeriodic[2] && !Readparameters::helpRequested) {
         if(myRank == MASTER_RANK) cerr << "You set boundaries.periodic_z = no but you didn't load any system boundary condition using the option boundaries.boundary, are you sure this is correct?" << endl;
      }
   }
   
   for (it = sysBoundaryCondList.begin();
        it != sysBoundaryCondList.end();
        it++) {
      if(*it == "Outflow") {
         if(this->addSysBoundary(new SBC::Outflow, project, t) == false) {
            if(myRank == MASTER_RANK) cerr << "Error in adding Outflow boundary." << endl;
            success = false;
         }
         isThisDynamic = isThisDynamic|this->getSysBoundary(sysboundarytype::OUTFLOW)->isDynamic();
         bool faces[6];
         this->getSysBoundary(sysboundarytype::OUTFLOW)->getFaces(&faces[0]);
         if((faces[0] || faces[1]) && isPeriodic[0]) {
            if(myRank == MASTER_RANK) cerr << "You set boundaries.periodic_x = yes and load Outflow system boundary conditions on the x+ or x- face, are you sure this is correct?" << endl;
         }
         if((faces[2] || faces[3]) && isPeriodic[1]) {
            if(myRank == MASTER_RANK) cerr << "You set boundaries.periodic_y = yes and load Outflow system boundary conditions on the y+ or y- face, are you sure this is correct?" << endl;
         }
         if((faces[4] || faces[5]) && isPeriodic[2]) {
            if(myRank == MASTER_RANK) cerr << "You set boundaries.periodic_z = yes and load Outflow system boundary conditions on the z+ or z- face, are you sure this is correct?" << endl;
         }
         if((faces[0] || faces[1]) && P::xcells_ini < 5) {
            if(myRank == MASTER_RANK) cerr << "You load Outflow system boundary conditions on the x+ or x- face but there is not enough cells in that direction to make sense." << endl;
            exit(1);
         }
         if((faces[2] || faces[3]) && P::ycells_ini < 5) {
            if(myRank == MASTER_RANK) cerr << "You load Outflow system boundary conditions on the y+ or y- face but there is not enough cells in that direction to make sense." << endl;
            exit(1);
         }
         if((faces[4] || faces[5]) && P::zcells_ini < 5) {
            if(myRank == MASTER_RANK) cerr << "You load Outflow system boundary conditions on the z+ or z- face but there is not enough cells in that direction to make sense." << endl;
            exit(1);
         }
      }
      if(*it == "Ionosphere") {
         if(this->addSysBoundary(new SBC::Ionosphere, project, t) == false) {
            if(myRank == MASTER_RANK) cerr << "Error in adding Ionosphere boundary." << endl;
            success = false;
         }
         if(this->addSysBoundary(new SBC::DoNotCompute, project, t) == false) {
            if(myRank == MASTER_RANK) cerr << "Error in adding DoNotCompute boundary (for Ionosphere)." << endl;
            success = false;
         }
         isThisDynamic = isThisDynamic|
         this->getSysBoundary(sysboundarytype::IONOSPHERE)->isDynamic();
      }
      if(*it == "Maxwellian") {
         if(this->addSysBoundary(new SBC::SetMaxwellian, project, t) == false) {
            if(myRank == MASTER_RANK) cerr << "Error in adding Maxwellian boundary." << endl;
            success = false;
         }
         isThisDynamic = isThisDynamic|
         this->getSysBoundary(sysboundarytype::SET_MAXWELLIAN)->isDynamic();
         bool faces[6];
         this->getSysBoundary(sysboundarytype::SET_MAXWELLIAN)->getFaces(&faces[0]);
         if((faces[0] || faces[1]) && isPeriodic[0]) {
            if(myRank == MASTER_RANK) cerr << "You set boundaries.periodic_x = yes and load Maxwellian system boundary conditions on the x+ or x- face, are you sure this is correct?" << endl;
         }
         if((faces[2] || faces[3]) && isPeriodic[1]) {
            if(myRank == MASTER_RANK) cerr << "You set boundaries.periodic_y = yes and load Maxwellian system boundary conditions on the y+ or y- face, are you sure this is correct?" << endl;
         }
         if((faces[4] || faces[5]) && isPeriodic[2]) {
            if(myRank == MASTER_RANK) cerr << "You set boundaries.periodic_z = yes and load Maxwellian system boundary conditions on the z+ or z- face, are you sure this is correct?" << endl;
         }
         if((faces[0] || faces[1]) && P::xcells_ini < 5) {
            if(myRank == MASTER_RANK) cerr << "You load Maxwellian system boundary conditions on the x+ or x- face but there is not enough cells in that direction to make sense." << endl;
            exit(1);
         }
         if((faces[2] || faces[3]) && P::ycells_ini < 5) {
            if(myRank == MASTER_RANK) cerr << "You load Maxwellian system boundary conditions on the y+ or y- face but there is not enough cells in that direction to make sense." << endl;
            exit(1);
         }
         if((faces[4] || faces[5]) && P::zcells_ini < 5) {
            if(myRank == MASTER_RANK) cerr << "You load Maxwellian system boundary conditions on the z+ or z- face but there is not enough cells in that direction to make sense." << endl;
            exit(1);
         }
      }
   }
   
   list<SBC::SysBoundaryCondition*>::iterator it2;
   for (it2 = sysBoundaries.begin();
        it2 != sysBoundaries.end();
        it2++
       ) {
      (*it2)->setPeriodicity(isPeriodic);
      }
   
   return success;
}

bool SysBoundary::checkRefinement(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

   // Set is used to avoid storing duplicates - each cell only needs to be checked once
   std::set<CellID> innerBoundaryCells;
   std::set<CellID> outerBoundaryCells;

   int innerBoundaryRefLvl = -1;
   int outerBoundaryRefLvl = -1;
   
   // Collect cells by sysboundarytype
   for (auto cellId : mpiGrid.get_cells()) {
      SpatialCell* cell = mpiGrid[cellId];
      if(cell) {
         if (cell->sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
            innerBoundaryCells.insert(cellId);
            innerBoundaryRefLvl = mpiGrid.get_refinement_level(cellId);
            if (cell->sysBoundaryLayer == 1) {
               // Add non-boundary neighbors of layer 1 cells
               auto* nbrPairVector = mpiGrid.get_neighbors_of(cellId,FULL_NEIGHBORHOOD_ID);
               for (auto nbrPair : *nbrPairVector) {
                  if(nbrPair.first != INVALID_CELLID) {
                     innerBoundaryCells.insert(nbrPair.first);
                  }
               }
            }
         } else if (cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
                    cell->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE) {
            outerBoundaryCells.insert(cellId);
            outerBoundaryRefLvl = mpiGrid.get_refinement_level(cellId);
            // Add non-boundary neighbors of outer boundary cells
            auto* nbrPairVector = mpiGrid.get_neighbors_of(cellId,FULL_NEIGHBORHOOD_ID);
            for (auto nbrPair : *nbrPairVector) {
               if(nbrPair.first != INVALID_CELLID) {
                  outerBoundaryCells.insert(nbrPair.first);
               }
            }
         }
      }
   }

   for (auto cellId : innerBoundaryCells) {
      if (cellId != INVALID_CELLID && mpiGrid.get_refinement_level(cellId) != innerBoundaryRefLvl) {
         return false;
      }
   }

   for (auto cellId : outerBoundaryCells) {
      if (cellId != INVALID_CELLID && mpiGrid.get_refinement_level(cellId) != outerBoundaryRefLvl) {
         // cout << "Failed refinement check " << cellId << " " << mpiGrid.get_refinement_level(cellId) << " "<< outerBoundaryRefLvl << endl;
         return false;
      }
   }
   
   return true;
}


bool belongsToLayer(const int layer, const int x, const int y, const int z,
                    FsGrid< fsgrids::technical, 2>& technicalGrid) {
   
   bool belongs = false;
   
   // loop through all neighbors (including diagonals)
   for (int ix = -1; ix <= 1; ++ix) {
      for (int iy = -1; iy <= 1; ++iy) {
         for (int iz = -1; iz <= 1; ++iz) {
            
            // not strictly necessary but logically we should not consider the cell itself
            // among its neighbors.
            if( ( ix == 0 && iy == 0 && iz == 0 ) || !technicalGrid.get(x+ix,y+iy,z+iz)) {
               continue;
            }
            
            if(layer == 1 && technicalGrid.get(x+ix,y+iy,z+iz)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
               // in the first layer, boundary cell belongs if it has a non-boundary neighbor
               belongs = true;
               return belongs;
               
            } else if (layer > 1 && technicalGrid.get(x+ix,y+iy,z+iz)->sysBoundaryLayer == layer - 1) {
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
bool SysBoundary::classifyCells(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                FsGrid< fsgrids::technical, 2> & technicalGrid) {
   bool success = true;
   vector<CellID> cells = mpiGrid.get_cells();
   auto localSize = technicalGrid.getLocalSize().data();
   
   /*set all cells to default value, not_sysboundary*/
   for(uint i=0; i<cells.size(); i++) {
      mpiGrid[cells[i]]->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
   }
   #pragma omp parallel for collapse(3)
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            technicalGrid.get(x,y,z)->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
            technicalGrid.get(x,y,z)->sysBoundaryLayer = 0;
            technicalGrid.get(x,y,z)->maxFsDt = std::numeric_limits<Real>::max();
         }
      }
   }
   
   /*
     loop through sysboundaries and let all sysboundaries set in local
   cells if they are part of which sysboundary (cell location needs to
   be updated by now. No remote data needed/available, so assignement
   has to be based individually on each cells location
   */
   list<SBC::SysBoundaryCondition*>::iterator it;
   for (it = sysBoundaries.begin(); it != sysBoundaries.end(); it++) {
      success = success && (*it)->assignSysBoundary(mpiGrid,technicalGrid);
   }
   

   // communicate boundary assignments (sysBoundaryFlag and
   // sysBoundaryLayer communicated)
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
   mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_NEIGHBORHOOD_ID);

   // set distance 1 cells to boundary cells, that have neighbors which are normal cells
   for(uint i=0; i<cells.size(); i++) {
      mpiGrid[cells[i]]->sysBoundaryLayer=0; /*Initial value*/
      if(mpiGrid[cells[i]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY ) {
         const auto* nbrs = mpiGrid.get_neighbors_of(cells[i],SYSBOUNDARIES_NEIGHBORHOOD_ID);
         for(uint j=0; j<(*nbrs).size(); j++) {
            if((*nbrs)[j].first!=0 ) {
               if(mpiGrid[(*nbrs)[j].first]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ) {
                  mpiGrid[cells[i]]->sysBoundaryLayer=1;
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
   uint maxLayers=22;//max(max(P::xcells_ini, P::ycells_ini), P::zcells_ini);
   for(uint layer=1;layer<maxLayers;layer++){
      for(uint i=0; i<cells.size(); i++) {
         if(mpiGrid[cells[i]]->sysBoundaryLayer==0){
            const auto* nbrs = mpiGrid.get_neighbors_of(cells[i],SYSBOUNDARIES_NEIGHBORHOOD_ID);
            for(uint j=0; j<(*nbrs).size(); j++) {
               if((*nbrs)[j].first!=0 && (*nbrs)[j].first!=cells[i] ) {
                  if(mpiGrid[(*nbrs)[j].first]->sysBoundaryLayer==layer) { 
                     mpiGrid[cells[i]]->sysBoundaryLayer=layer+1;
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
   for(uint i=0; i<cells.size(); i++) {
      if(mpiGrid[cells[i]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
         mpiGrid[cells[i]]->sysBoundaryLayer != 1 &&
         mpiGrid[cells[i]]->sysBoundaryLayer != 2
      ) {
         mpiGrid[cells[i]]->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
      }
   }

   // The following is done so that everyone knows their neighbour's
   // layer flags.  This is needed for the correct use of the system
   // boundary local communication patterns.
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
   mpiGrid.update_copies_of_remote_neighbors(FULL_NEIGHBORHOOD_ID);
   
   
   // Now the layers need to be set on fsgrid too
   // In dccrg initialization the max number of boundary layers is set to 3.
   const uint MAX_NUMBER_OF_BOUNDARY_LAYERS = 3 * pow(2,mpiGrid.get_maximum_refinement_level());
   
   technicalGrid.updateGhostCells();
   
   // loop through max number of layers
   for(uint layer = 1; layer <= MAX_NUMBER_OF_BOUNDARY_LAYERS; ++layer) {
      
      // loop through all cells in grid
      #pragma omp parallel for collapse(3)
      for (int x = 0; x < localSize[0]; ++x) {
         for (int y = 0; y < localSize[1]; ++y) {
            for (int z = 0; z < localSize[2]; ++z) {
               
               // for the first layer, consider all cells that belong to a boundary, for other layers
               // consider all cells that have not yet been labeled.
               if((layer == 1 && technicalGrid.get(x,y,z)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) ||
                  (layer > 1 && technicalGrid.get(x,y,z)->sysBoundaryLayer == 0)) {
                  
                  if (belongsToLayer(layer, x, y, z, technicalGrid)) {
                     
                     technicalGrid.get(x,y,z)->sysBoundaryLayer = layer;
                     
                     if (layer > 2 && technicalGrid.get(x,y,z)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
                        technicalGrid.get(x,y,z)->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
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
            if (technicalGrid.get(x,y,z)->sysBoundaryLayer == 0 && technicalGrid.get(x,y,z)->sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
               technicalGrid.get(x,y,z)->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
            }
         }
      }
   }
   
   technicalGrid.updateGhostCells();
   
   return success;
}

/*!\brief Apply the initial state to all system boundary cells.
 * Loops through all SysBoundaryConditions and calls the corresponding applyInitialState
 * function. This function must apply the initial state for all existing particle species.
 * 
 * \param mpiGrid Grid
 * \retval success If true, the application of all system boundary states succeeded.
 */
bool SysBoundary::applyInitialState(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   Project& project
) {
   bool success = true;
   
   list<SBC::SysBoundaryCondition*>::iterator it;
   for (it = sysBoundaries.begin();
        it != sysBoundaries.end();
        it++) {
      if(                                                        // This is to skip the reapplication
         Parameters::isRestart == true                           // When not restarting
         && (*it)->doApplyUponRestart() == false                 // When reapplicaiton is not requested
         && (*it)->getIndex() != sysboundarytype::IONOSPHERE     // But this is to force it when we have either IONOSPHERE
         && (*it)->getIndex() != sysboundarytype::SET_MAXWELLIAN // or SET_MAXWELLIAN as otherwise the POP_METADA are not properly set
      ) {
         continue;
      }
      if((*it)->applyInitialState(mpiGrid, perBGrid, project) == false) {
         cerr << "ERROR: " << (*it)->getName() << " system boundary condition initial state not applied correctly." << endl;
         success = false;
      }
   }

   return success;
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
 */
void SysBoundary::applySysBoundaryVlasovConditions(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   creal& t
) {
   if(sysBoundaries.size()==0) {
      return; //no system boundaries
   }

   /*Transfer along boundaries*/
   // First the small stuff without overlapping in an extended neighbourhood:
   SpatialCell::set_mpi_transfer_type(
      Transfer::CELL_PARAMETERS|
      Transfer::POP_METADATA|
      Transfer::CELL_SYSBOUNDARYFLAG,true);
   mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
   
   // Loop over existing particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      SpatialCell::setCommunicatedSpecies(popID);

      // Then the block data in the reduced neighbourhood:
      int timer=phiprof::initializeTimer("Start comm of cell and block data","MPI");
      phiprof::start(timer);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,true);
      mpiGrid.start_remote_neighbor_copy_updates(SYSBOUNDARIES_NEIGHBORHOOD_ID);
      phiprof::stop(timer);

      timer=phiprof::initializeTimer("Compute process inner cells");
      phiprof::start(timer);

      // Compute Vlasov boundary condition on system boundary/process inner cells
      vector<CellID> localCells;
      getBoundaryCellList(mpiGrid,mpiGrid.get_local_cells_not_on_process_boundary(SYSBOUNDARIES_NEIGHBORHOOD_ID),localCells);
   
      #pragma omp parallel for
      for (uint i=0; i<localCells.size(); i++) {
         cuint sysBoundaryType = mpiGrid[localCells[i]]->sysBoundaryFlag;
         this->getSysBoundary(sysBoundaryType)->vlasovBoundaryCondition(mpiGrid,localCells[i],popID);
      }
      phiprof::stop(timer);
   
      timer=phiprof::initializeTimer("Wait for receives","MPI","Wait");
      phiprof::start(timer);
      mpiGrid.wait_remote_neighbor_copy_update_receives(SYSBOUNDARIES_NEIGHBORHOOD_ID);
      phiprof::stop(timer);

      // Compute vlasov boundary on system boundary/process boundary cells
      timer=phiprof::initializeTimer("Compute process boundary cells");
      phiprof::start(timer);
      vector<CellID> boundaryCells;
      getBoundaryCellList(mpiGrid,mpiGrid.get_local_cells_on_process_boundary(SYSBOUNDARIES_NEIGHBORHOOD_ID),boundaryCells);
      #pragma omp parallel for
      for (uint i=0; i<boundaryCells.size(); i++) {
         cuint sysBoundaryType = mpiGrid[boundaryCells[i]]->sysBoundaryFlag;
         this->getSysBoundary(sysBoundaryType)->vlasovBoundaryCondition(mpiGrid, boundaryCells[i],popID);
      }
      phiprof::stop(timer);

      timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
      phiprof::start(timer);
      mpiGrid.wait_remote_neighbor_copy_update_sends();
      phiprof::stop(timer);

      // WARNING Blocks are changed but lists not updated now, if you need to use/communicate them before the next update is done, add an update here.
      updateRemoteVelocityBlockLists(mpiGrid, popID);

   } // for-loop over populations
}

/*! Get a pointer to the SysBoundaryCondition of given index.
 * \param sysBoundaryType Type of the system boundary condition to return
 * \return Pointer to the instance of the SysBoundaryCondition. NULL if sysBoundaryType is invalid.
 */
SBC::SysBoundaryCondition* SysBoundary::getSysBoundary(cuint sysBoundaryType) const {
   auto it = indexToSysBoundary.find(sysBoundaryType);
   if(it != indexToSysBoundary.end()){
      return it->second;
   }
   else{
      std::cerr << "System boundary "<< sysBoundaryType << " is invalid  " << __FILE__ << ":" << __LINE__ << std::endl;
      return NULL;
   }
}

/*! Get the number of SysBoundaryConditions stored in SysBoundary.
 * \retval size Number of SysBoundaryConditions stored in SysBoundary.
 */
unsigned int SysBoundary::size() const {return sysBoundaries.size();}

/*! Get a bool telling whether any system boundary condition is dynamic in time (and thus needs updating).
 * \retval isThisDynamic Is any system boundary condition dynamic in time.
 */
bool SysBoundary::isDynamic() const {return isThisDynamic;}

/*! Get a bool telling whether the system is periodic in the queried direction.
 * \param direction 0: x, 1: y, 2: z.
 * \retval isPeriodic Is the system periodic in the queried direction.
 */
bool SysBoundary::isBoundaryPeriodic(uint direction) const {return isPeriodic[direction];}

/*! Get a vector containing the cellID of all cells which are not DO_NOT_COMPUTE or NOT_SYSBOUNDARY in the vector of cellIDs passed to the function.
 * 
 * \param mpiGrid Grid
 * \param cellList Vector of cellIDs in which to look for boundary cells
 * \param boundaryCellList Vector of boundary the cells' cellIDs
 * \retval Returns true if the operation is successful
 */
bool getBoundaryCellList(
   const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const vector<uint64_t>& cellList,
   vector<uint64_t>& boundaryCellList
){
   boundaryCellList.clear();
   for (size_t cell=0; cell<cellList.size(); ++cell) {
      const CellID cellID = cellList[cell];
      if(mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
         mpiGrid[cellID]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) continue;
      boundaryCellList.push_back(cellID);
   }
   return true;
}

/*! Updates all NonsysboundaryCells into an internal map. This should be called in loadBalance.
 * \param mpiGrid The DCCRG grid
 * \retval Returns true if the operation is successful
 */
bool SysBoundary::updateSysBoundariesAfterLoadBalance(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   phiprof::start("updateSysBoundariesAfterLoadBalance");
   vector<uint64_t> local_cells_on_boundary;
   getBoundaryCellList(mpiGrid, mpiGrid.get_cells(), local_cells_on_boundary);
   // Loop over sysboundaries:
   for( std::list<SBC::SysBoundaryCondition*>::iterator it = sysBoundaries.begin(); it != sysBoundaries.end(); ++it ) {
      (*it)->updateSysBoundaryConditionsAfterLoadBalance(mpiGrid, local_cells_on_boundary);
   }

   phiprof::stop("updateSysBoundariesAfterLoadBalance");
   return true;
}

