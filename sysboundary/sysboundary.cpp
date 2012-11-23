/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/*! \file sysboundary.cpp
 * \brief Implementation of the class SysBoundary.
 */

#include <cstdlib>
#include <iostream>

#include "sysboundary.h"
# include "../grid.h"

using namespace std;

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
   Readparameters::addComposing("boundaries.boundary", "List of boundary condition (BC) types to be used. Each boundary condition to be used has to be on a new line boundary = YYY. Available (20120820) are outflow ionosphere maxwellian.");
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
 * \retval success If true, the given SBC::SysBoundaryCondition was added successfully.
 */
bool SysBoundary::addSysBoundary(
   SBC::SysBoundaryCondition* bc,
   Project &project,
   creal& t
) {
   sysBoundaries.push_back(bc);
   if(sysBoundaries.size() > 1) {
      sysBoundaries.sort(precedenceSort);
   }
   
   bool success = true;
   if(bc->initSysBoundary(t, project) == false) {
      success = false;
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
      }
   }
   
   return success;
}

/*!\brief Classify all simulation cells with respect to the system boundary conditions.
 * 
 * Loops through all cells and and for each assigns the correct sysBoundaryFlag depending on
 * the return value of each SysBoundaryCondition's assignSysBoundary.
 */
bool SysBoundary::classifyCells(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   bool success = true;
   vector<CellID> cells = mpiGrid.get_cells();
   
   for(uint i=0; i<cells.size(); i++) {
      mpiGrid[cells[i]]->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
   }
   
   list<SBC::SysBoundaryCondition*>::iterator it;
   for (it = sysBoundaries.begin();
        it != sysBoundaries.end();
        it++) {
      success = success && (*it)->assignSysBoundary(mpiGrid);
   }
   
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
   mpiGrid.update_remote_neighbor_data(SYSBOUNDARIES_NEIGHBORHOOD_ID);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   /*set distance 1 cells to boundary cells, that have neighbors which are normal cells */
   for(uint i=0; i<cells.size(); i++) {
      mpiGrid[cells[i]]->sysBoundaryLayer=0; /*Initial value*/
      if(mpiGrid[cells[i]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY ) {
         const std::vector<CellID>* nbrs = mpiGrid.get_neighbors(cells[i],SYSBOUNDARIES_NEIGHBORHOOD_ID);
         for(uint j=0; j<(*nbrs).size(); j++) {
            if((*nbrs)[j]!=0 ) {
               if(mpiGrid[(*nbrs)[j]]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ) {
                  mpiGrid[cells[i]]->sysBoundaryLayer=1;
               }
            }
         }
      }
   }
   
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
   mpiGrid.update_remote_neighbor_data(SYSBOUNDARIES_NEIGHBORHOOD_ID);
   
   /*Compute distances*/
   /*Construct vector with all cells, both local and remote*/
   uint maxLayers=10;
   for(uint layer=1;layer<maxLayers;layer++){
      for(uint i=0; i<cells.size(); i++) {
         if(mpiGrid[cells[i]]->sysBoundaryLayer==0){
            const std::vector<CellID>* nbrs = mpiGrid.get_neighbors(cells[i],SYSBOUNDARIES_NEIGHBORHOOD_ID);
            for(uint j=0; j<(*nbrs).size(); j++) {
               if((*nbrs)[j]!=0 && (*nbrs)[j]!=cells[i] ) {
                  if(mpiGrid[(*nbrs)[j]]->sysBoundaryLayer==layer) { 
                     mpiGrid[cells[i]]->sysBoundaryLayer=layer+1;
                     break;
                  }
               }
            }
         }
      }
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
      mpiGrid.update_remote_neighbor_data(SYSBOUNDARIES_NEIGHBORHOOD_ID);
   }
   
   for(uint i=0; i<cells.size(); i++) {
      if(mpiGrid[cells[i]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
         mpiGrid[cells[i]]->sysBoundaryLayer != 1 &&
         mpiGrid[cells[i]]->sysBoundaryLayer != 2
      ) {
         mpiGrid[cells[i]]->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
      }
   }
   
   return success;
}

/*!\brief Apply the initial state to all system boundary cells.
 * 
 * Loops through all SysBoundaryConditions and calls the corresponding applyInitialState
 * function.
 * \retval success If true, the application of all system boundary states succeeded.
 */
bool SysBoundary::applyInitialState(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   Project& project
) {
   bool success = true;
   using namespace sysboundarytype;
   
   list<SBC::SysBoundaryCondition*>::iterator it;
   for (it = sysBoundaries.begin();
        it != sysBoundaries.end();
        it++) {
      if((*it)->applyInitialState(mpiGrid, project) == false) {
         cerr << "ERROR: " << (*it)->getName() << " system boundary condition initial state not applied correctly." << endl;
         success = false;
      }
   }
   
   return success;
}

/*!\brief Apply the Vlasov system boundary conditions to all system boundary cells at time t.
 *
 * Loops through all SysBoundaryConditions and calls the corresponding vlasovBoundaryCondition()
 * function.
 */
void SysBoundary::applySysBoundaryVlasovConditions(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t) {
   if(sysBoundaries.size()==0) {
      return; //no system boundaries
   }

   /*Transfer along boundaries*/
   // First the small stuff without overlapping in an extended neighbourhood:
   SpatialCell::set_mpi_transfer_type(
      Transfer::CELL_PARAMETERS|
      Transfer::CELL_SYSBOUNDARYFLAG,true);
   mpiGrid.update_remote_neighbor_data(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
   // Then the block data in the reduced neighbourhood:
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,true);
   
   int timer=phiprof::initializeTimer("Start comm of cell and block data","MPI");
   phiprof::start(timer);
   mpiGrid.start_remote_neighbor_data_update(SYSBOUNDARIES_NEIGHBORHOOD_ID);
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Compute process inner cells");
   phiprof::start(timer);
   // Compute Vlasov boundary condition on system boundary/on process inner cells
   const vector<CellID> localCells = mpiGrid.get_cells_with_local_neighbors(SYSBOUNDARIES_NEIGHBORHOOD_ID);
# pragma omp parallel for
   for (uint i=0; i<localCells.size(); i++) {
      cuint sysBoundaryType = mpiGrid[localCells[i]]->sysBoundaryFlag;
      if(sysBoundaryType == sysboundarytype::DO_NOT_COMPUTE ||
         sysBoundaryType == sysboundarytype::NOT_SYSBOUNDARY) continue;
      this->getSysBoundary(sysBoundaryType)->vlasovBoundaryCondition(mpiGrid, localCells[i]);
   }
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Wait for receives","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_neighbor_data_update_receives(SYSBOUNDARIES_NEIGHBORHOOD_ID);
   phiprof::stop(timer);
   
   // Compute vlasov boundary on system boundary/process boundary cells
   timer=phiprof::initializeTimer("Compute process boundary cells");
   phiprof::start(timer);
   const vector<CellID> boundaryCells = mpiGrid.get_cells_with_remote_neighbor(SYSBOUNDARIES_NEIGHBORHOOD_ID);
# pragma omp parallel for
   for (uint i=0; i<boundaryCells.size(); i++) {
      cuint sysBoundaryType = mpiGrid[boundaryCells[i]]->sysBoundaryFlag;
      if(sysBoundaryType == sysboundarytype::DO_NOT_COMPUTE ||
         sysBoundaryType == sysboundarytype::NOT_SYSBOUNDARY) continue;
      this->getSysBoundary(sysBoundaryType)->vlasovBoundaryCondition(mpiGrid, boundaryCells[i]);
   }
   phiprof::stop(timer);
   
   timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
   phiprof::start(timer);
   mpiGrid.wait_neighbor_data_update_sends();
   phiprof::stop(timer);
   
   //FIXME, we should have a boolean that tells us if this is needed or not... 
   updateRemoteVelocityBlockLists(mpiGrid);
   adjustVelocityBlocks(mpiGrid);

   /*
  // NOTE, I do not think these are needed.  VEL_BLOCK_DATA is
  // transferred in the leveque solver, and the moments are needed when
  // computing fields -> should be handled there
   
   SpatialCell::set_mpi_transfer_type(
   Transfer::CELL_RHO_RHOV| Transfer::VEL_BLOCK_DATA );
   mpiGrid.update_remote_neighbor_data(xxx);
   */
}   

/*! Get a pointer to the SysBoundaryCondition of given index.
 * \retval ptr Pointer to the instance of the SysBoundaryCondition.
 */
SBC::SysBoundaryCondition* SysBoundary::getSysBoundary(cuint sysBoundaryType) const {
   return indexToSysBoundary.find(sysBoundaryType)->second;
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
