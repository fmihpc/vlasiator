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

#include <boost/assign/list_of.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip> // for setprecision()
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>

#include "grid.h"

#include "vlasovmover.h"
#include "definitions.h"
#include "mpiconversion.h"
#include "logger.h"
#include "parameters.h"

#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"
#include "transferstencil.h"

#include "vlsvwriter2.h" 
#include "fieldsolver.h"
#include "projects/project.h"
#include "fileio.h"


using namespace std;
using namespace phiprof;

extern Logger logFile, diagnostic;

void initVelocityGridGeometry();
void initSpatialCellCoordinates(dccrg::Dccrg<SpatialCell>& mpiGrid);
bool applyInitialState(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   Project& project
);
bool adjust_local_velocity_blocks(dccrg::Dccrg<SpatialCell>& mpiGrid);

void initializeGrid(
   int argn,
   char **argc,
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   SysBoundary& sysBoundaries,
   Project& project
) {
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   
   // initialize velocity grid of spatial cells before creating cells in dccrg.initialize
   initVelocityGridGeometry();

   // Init Zoltan:
   float zoltanVersion;
   if (Zoltan_Initialize(argn,argc,&zoltanVersion) != ZOLTAN_OK) {
      if(myRank == MASTER_RANK) cerr << "\t ERROR: Zoltan initialization failed." << endl;
      exit(1);
   } else {
      logFile << "\t Zoltan " << zoltanVersion << " initialized successfully" << std::endl << writeVerbose;
   }
   
   mpiGrid.set_geometry(
      P::xcells_ini, P::ycells_ini, P::zcells_ini,
      P::xmin, P::ymin, P::zmin,
      P::dx_ini, P::dy_ini, P::dz_ini
      );

   
   MPI_Comm comm = MPI_COMM_WORLD;
   mpiGrid.initialize(
      comm,
      &P::loadBalanceAlgorithm[0],
      // neighborhood size
      #ifdef SOLVER_KT
      1, // kt needs 0 but field volume average calculation needs 1
      #elif defined SOLVER_LEVEQUE
      2,
      #endif
      0, // maximum refinement level
      sysBoundaries.isBoundaryPeriodic(0),
      sysBoundaries.isBoundaryPeriodic(1),
      sysBoundaries.isBoundaryPeriodic(2)
   );

   // set reduced neighborhoods
   typedef dccrg::Types<3>::neighborhood_item_t neigh_t;
   
   // set a reduced neighborhood for field solver
   std::vector<neigh_t> nearestneighbor_neighborhood;
   for (int z = -1; z <= 1; z++) {
      for (int y = -1; y <= 1; y++) {
         for (int x = -1; x <= 1; x++) {
            if (x == 0 && y == 0 && z == 0) {
               continue;
            }
            
            neigh_t offsets = {{x, y, z}};
            nearestneighbor_neighborhood.push_back(offsets);
         }
      }
   }
   
   if (!mpiGrid.add_remote_update_neighborhood(FIELD_SOLVER_NEIGHBORHOOD_ID, nearestneighbor_neighborhood)) {
      std::cerr << __FILE__ << ":" << __LINE__
         << " Couldn't set field solver neighborhood"
         << std::endl;
      abort();
   }
   
   // set a reduced neighborhood for all possible communication in  vlasov solver
   //FIXME, the +2 neighbors can be removed as we do not receive from +2, do check though...
   const std::vector<neigh_t> vlasov_neighborhood
      = boost::assign::list_of<neigh_t>
      (boost::assign::list_of( 0)( 0)(-2))
      (boost::assign::list_of(-1)(-1)(-1))
      (boost::assign::list_of( 0)(-1)(-1))
      (boost::assign::list_of( 1)(-1)(-1))
      (boost::assign::list_of(-1)( 0)(-1))
      (boost::assign::list_of( 0)( 0)(-1))
      (boost::assign::list_of( 1)( 0)(-1))
      (boost::assign::list_of(-1)( 1)(-1))
      (boost::assign::list_of( 0)( 1)(-1))
      (boost::assign::list_of( 1)( 1)(-1))
      (boost::assign::list_of( 0)(-2)( 0))
      (boost::assign::list_of(-1)(-1)( 0))
      (boost::assign::list_of( 0)(-1)( 0))
      (boost::assign::list_of( 1)(-1)( 0))
      (boost::assign::list_of(-2)( 0)( 0))
      (boost::assign::list_of(-1)( 0)( 0))
      (boost::assign::list_of( 1)( 0)( 0))
      (boost::assign::list_of( 2)( 0)( 0))
      (boost::assign::list_of(-1)( 1)( 0))
      (boost::assign::list_of( 0)( 1)( 0))
      (boost::assign::list_of( 1)( 1)( 0))
      (boost::assign::list_of( 0)( 2)( 0))
      (boost::assign::list_of(-1)(-1)( 1))
      (boost::assign::list_of( 0)(-1)( 1))
      (boost::assign::list_of( 1)(-1)( 1))
      (boost::assign::list_of(-1)( 0)( 1))
      (boost::assign::list_of( 0)( 0)( 1))
      (boost::assign::list_of( 1)( 0)( 1))
      (boost::assign::list_of(-1)( 1)( 1))
      (boost::assign::list_of( 0)( 1)( 1))
      (boost::assign::list_of( 1)( 1)( 1))
      (boost::assign::list_of( 0)( 0)( 2));
   
   if (!mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_NEIGHBORHOOD_ID, vlasov_neighborhood)) {
      std::cerr << __FILE__ << ":" << __LINE__
                << " Couldn't set vlasov solver neighborhood"
                << std::endl;
      abort();
   }

   
   // A reduced neighborhood for vlasov distribution function receives
   const std::vector<neigh_t> vlasov_density_neighborhood
      = boost::assign::list_of<neigh_t>
      (boost::assign::list_of( 0)( 0)(-1))
      (boost::assign::list_of( 0)( 0)( 1))
      (boost::assign::list_of( 0)(-1)( 0))
      (boost::assign::list_of( 0)( 1)( 0))
      (boost::assign::list_of(-1)( 0)( 0))
      (boost::assign::list_of( 1)( 0)( 0))
      (boost::assign::list_of(-2)( 0)( 0))
      (boost::assign::list_of( 0)(-2)( 0))
      (boost::assign::list_of( 0)( 0)(-2));
   
   if (!mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_DENSITY_NEIGHBORHOOD_ID, vlasov_density_neighborhood)) {
      std::cerr << __FILE__ << ":" << __LINE__
                << " Couldn't set field solver neighborhood"
                << std::endl;
      abort();
   }


   if (!mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_FLUXES_NEIGHBORHOOD_ID, nearestneighbor_neighborhood)) {
      std::cerr << __FILE__ << ":" << __LINE__
                << " Couldn't set field solver neighborhood"
                << std::endl;
      abort();
   }

   if (!mpiGrid.add_remote_update_neighborhood(SYSBOUNDARIES_NEIGHBORHOOD_ID, nearestneighbor_neighborhood)) {
      std::cerr << __FILE__ << ":" << __LINE__
                << " Couldn't set system boundaries neighborhood"
                << std::endl;
      abort();
   }
   
   
   mpiGrid.set_partitioning_option("IMBALANCE_TOL", P::loadBalanceTolerance);
   phiprof::start("Initial load-balancing");
   if (myRank == MASTER_RANK) logFile << "(INIT): Starting initial load balance." << endl << writeVerbose;
   mpiGrid.balance_load();
   phiprof::stop("Initial load-balancing");
   
   if (myRank == MASTER_RANK) logFile << "(INIT): Set initial state." << endl << writeVerbose;
   phiprof::start("Set initial state");
   
   phiprof::start("Set spatial cell coordinates");
   initSpatialCellCoordinates(mpiGrid);
   phiprof::stop("Set spatial cell coordinates");
   
   phiprof::start("Initialize system boundary conditions");
   if(sysBoundaries.initSysBoundaries(project, P::t_min) == false) {
      if (myRank == MASTER_RANK) cerr << "Error in initialising the system boundaries." << endl;
      exit(1);
   }
   phiprof::stop("Initialize system boundary conditions");
   
   // Initialise system boundary conditions (they need the initialised positions!!)
   phiprof::start("Classify cells (sys boundary conditions)");
   if(sysBoundaries.classifyCells(mpiGrid) == false) {
      cerr << "(MAIN) ERROR: System boundary conditions were not set correctly." << endl;
      exit(1);
   }
   phiprof::stop("Classify cells (sys boundary conditions)");
   
   if (P::isRestart){
      logFile << "Restart from "<< P::restartFileName << std::endl << writeVerbose;
      phiprof::start("Read restart");
      if ( readGrid(mpiGrid,P::restartFileName)== false ){
         logFile << "(MAIN) ERROR: restarting failed"<<endl;
         exit(1);
      }
      phiprof::stop("Read restart");
   } else {
      // Go through every spatial cell on this CPU, and create the initial state:
      phiprof::start("Apply initial state");
      if(applyInitialState(mpiGrid, project) == false) {
         cerr << "(MAIN) ERROR: Initial state was not applied correctly." << endl;
         exit(1);
      }
      
      phiprof::stop("Apply initial state");
      phiprof::start("Apply system boundary conditions state");
      if(sysBoundaries.applyInitialState(mpiGrid, project) == false) {
         cerr << " (MAIN) ERROR: System boundary conditions initial state was not applied correctly." << endl;
         exit(1);
      }
      phiprof::stop("Apply system boundary conditions state");
   }

   updateRemoteVelocityBlockLists(mpiGrid);
   adjustVelocityBlocks(mpiGrid,false); // do not initialize mover, mover has not yet been initialized here
   
   phiprof::initializeTimer("Fetch Neighbour data","MPI");
   phiprof::start("Fetch Neighbour data");
   // update complete spatial cell data 
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_DATA);
   mpiGrid.update_remote_neighbor_data(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop("Fetch Neighbour data");
   
   phiprof::stop("Set initial state");
   
   balanceLoad(mpiGrid);
}

// initialize velocity grid of spatial cells before creating cells in dccrg.initialize
void initVelocityGridGeometry(){
   spatial_cell::SpatialCell::vx_length = P::vxblocks_ini;
   spatial_cell::SpatialCell::vy_length = P::vyblocks_ini;
   spatial_cell::SpatialCell::vz_length = P::vzblocks_ini;
   spatial_cell::SpatialCell::max_velocity_blocks = 
      spatial_cell::SpatialCell::vx_length * spatial_cell::SpatialCell::vy_length * spatial_cell::SpatialCell::vz_length;
   spatial_cell::SpatialCell::vx_min = P::vxmin;
   spatial_cell::SpatialCell::vx_max = P::vxmax;
   spatial_cell::SpatialCell::vy_min = P::vymin;
   spatial_cell::SpatialCell::vy_max = P::vymax;
   spatial_cell::SpatialCell::vz_min = P::vzmin;
   spatial_cell::SpatialCell::vz_max = P::vzmax;
   spatial_cell::SpatialCell::grid_dvx = spatial_cell::SpatialCell::vx_max - spatial_cell::SpatialCell::vx_min;
   spatial_cell::SpatialCell::grid_dvy = spatial_cell::SpatialCell::vy_max - spatial_cell::SpatialCell::vy_min;
   spatial_cell::SpatialCell::grid_dvz = spatial_cell::SpatialCell::vz_max - spatial_cell::SpatialCell::vz_min;
   spatial_cell::SpatialCell::block_dvx = spatial_cell::SpatialCell::grid_dvx / spatial_cell::SpatialCell::vx_length;
   spatial_cell::SpatialCell::block_dvy = spatial_cell::SpatialCell::grid_dvy / spatial_cell::SpatialCell::vy_length;
   spatial_cell::SpatialCell::block_dvz = spatial_cell::SpatialCell::grid_dvz / spatial_cell::SpatialCell::vz_length;
   spatial_cell::SpatialCell::cell_dvx = spatial_cell::SpatialCell::block_dvx / block_vx_length;
   spatial_cell::SpatialCell::cell_dvy = spatial_cell::SpatialCell::block_dvy / block_vy_length;
   spatial_cell::SpatialCell::cell_dvz = spatial_cell::SpatialCell::block_dvz / block_vz_length;
   spatial_cell::SpatialCell::velocity_block_min_value = P::sparseMinValue;
}



void initSpatialCellCoordinates(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   vector<uint64_t> cells = mpiGrid.get_cells();
#pragma omp parallel for
   for (uint i=0; i<cells.size(); ++i) {
      mpiGrid[cells[i]]->parameters[CellParams::XCRD] = mpiGrid.get_cell_x_min(cells[i]);
      mpiGrid[cells[i]]->parameters[CellParams::YCRD] = mpiGrid.get_cell_y_min(cells[i]);
      mpiGrid[cells[i]]->parameters[CellParams::ZCRD] = mpiGrid.get_cell_z_min(cells[i]);
      mpiGrid[cells[i]]->parameters[CellParams::DX  ] = mpiGrid.get_cell_x_size(cells[i]);
      mpiGrid[cells[i]]->parameters[CellParams::DY  ] = mpiGrid.get_cell_y_size(cells[i]);
      mpiGrid[cells[i]]->parameters[CellParams::DZ  ] = mpiGrid.get_cell_z_size(cells[i]);
   }
}

bool applyInitialState(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   Project& project
) {
   typedef Parameters P;
   using namespace sysboundarytype;
   
   vector<uint64_t> cells = mpiGrid.get_cells();
   
   //  Go through every cell on this node and initialize the pointers to 
   // cpu memory, physical parameters and volume averages for each phase space 
   // point in the velocity grid. Velocity block neighbour list is also 
   // constructed here:
   // Each initialization has to be independent to avoid threading problems 
#pragma omp parallel for schedule(dynamic)
   //WARNING no threading here if setProjectCell has threading
   for (uint i=0; i<cells.size(); ++i) {
      SpatialCell* cell = mpiGrid[cells[i]];
      if(cell->sysBoundaryFlag != NOT_SYSBOUNDARY) continue;
      project.setCell(cell);
   }
   return true;
}

void balanceLoad(dccrg::Dccrg<SpatialCell>& mpiGrid){

   //set weights based on each cells LB weight counter
   vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i=0; i<cells.size(); ++i){
      if(mpiGrid[cells[i]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         // No load in that case
         mpiGrid.set_cell_weight(cells[i], 0.0);
      } else {
         if(P::maxAccelerationSubsteps!=1) {
            //use time-metric from solvers
            mpiGrid.set_cell_weight(cells[i], mpiGrid[cells[i]]->parameters[CellParams::LBWEIGHTCOUNTER]);
         }
         else{
            //No substepping in acceleration step, use number of blocks instead as metric as that provides slightly better results
            mpiGrid.set_cell_weight(cells[i], mpiGrid[cells[i]]->number_of_blocks);
         }
      }
      
      mpiGrid[cells[i]]->parameters[CellParams::LBWEIGHTCOUNTER]=0.0; //zero counter
   }
   
// tell other processes which velocity blocks exist in remote spatial cells
   phiprof::initializeTimer("Balancing load", "Load balance");
   phiprof::start("Balancing load");
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_SIZE_AND_LIST);
   mpiGrid.prepare_to_balance_load();
   
   // reserve space for velocity block data in arriving remote cells
   phiprof::start("Preparing receives");
   const boost::unordered_set<uint64_t>* incoming_cells = mpiGrid.get_balance_added_cells();
   std::vector<uint64_t> incoming_cells_list (incoming_cells->begin(),incoming_cells->end()); 
   
#pragma omp parallel for
   for(unsigned int i=0;i<incoming_cells_list.size();i++){
      uint64_t cell_id=incoming_cells_list[i];
      SpatialCell* cell = mpiGrid[cell_id];
      if (cell == NULL) {
         cerr << "No data for spatial cell " << cell_id << endl;
         abort();
      }
      cell->prepare_to_receive_blocks();
   }
   phiprof::stop("Preparing receives", incoming_cells_list.size(), "Spatial cells");


   phiprof::start("balance load");
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_DATA);
   mpiGrid.balance_load(true);
   phiprof::stop("balance load");

   phiprof::start("update block lists");
   //new partition, re/initialize blocklists of remote cells.
   updateRemoteVelocityBlockLists(mpiGrid);
   phiprof::stop("update block lists");

   phiprof::start("Init solvers");
   //need to re-initialize stencils and neighbors in leveque solver
   if (initializeMover(mpiGrid) == false) {
      logFile << "(MAIN): Vlasov propagator did not initialize correctly!" << endl << writeVerbose;
      exit(1);
   }

      // Initialize field propagator:
   if (initializeFieldPropagatorAfterRebalance(mpiGrid) == false) {
       logFile << "(MAIN): Field propagator did not initialize correctly!" << endl << writeVerbose;
       exit(1);
   }
   phiprof::stop("Init solvers");
   
   phiprof::stop("Balancing load");
}

//Compute which blocks have content, adjust local velocity blocks, and
//make sure remote cells are up-to-date and ready to receive
//data. Solvers are also updated so that their internal structures are
//ready for the new number of blocks.

bool adjustVelocityBlocks(dccrg::Dccrg<SpatialCell>& mpiGrid, bool reInitMover) {
   phiprof::initializeTimer("re-adjust blocks","Block adjustment");
   phiprof::start("re-adjust blocks");

   vector<uint64_t> cells = mpiGrid.get_cells();
   phiprof::start("Check for content");
#pragma omp parallel for  
   for (uint i=0; i<cells.size(); ++i) 
      mpiGrid[cells[i]]->update_all_block_has_content();     
   phiprof::stop("Check for content");
   phiprof::initializeTimer("Transfer block data","MPI");
   phiprof::start("Transfer block data");
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_HAS_CONTENT );
   mpiGrid.update_remote_neighbor_data(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop("Transfer block data");
   
   adjust_local_velocity_blocks(mpiGrid);
   updateRemoteVelocityBlockLists(mpiGrid);

   //re-init vlasovmover
   if(reInitMover) {
      phiprof::start("InitMoverAfterBlockChange");
      initMoverAfterBlockChange(mpiGrid);
      phiprof::stop("InitMoverAfterBlockChange");
   }
   phiprof::stop("re-adjust blocks");
   return true;
}

/*!
Adjusts velocity blocks in local spatial cells.

Doesn't adjust velocity blocks in remote cells, this has to be fixed
by calling updateRemoteVelocityBlockLists()

*/
bool adjust_local_velocity_blocks(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   phiprof::start("Adjusting blocks");

   const vector<uint64_t> cells = mpiGrid.get_cells();

#pragma omp parallel for
   for(unsigned int i=0;i<cells.size();i++){
      uint64_t cell_id=cells[i];
      SpatialCell* cell = mpiGrid[cell_id];
      if (cell == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__
                   << " No data for spatial cell " << cell_id
                   << endl;
         abort();
      }
      // gather spatial neighbor list
      const vector<uint64_t>* neighbors
      	= mpiGrid.get_neighbors(cell_id, VLASOV_SOLVER_NEIGHBORHOOD_ID);
      vector<SpatialCell*> neighbor_ptrs;

      neighbor_ptrs.reserve(neighbors->size());
      
      for (vector<uint64_t>::const_iterator
           neighbor_id = neighbors->begin();
           neighbor_id != neighbors->end();
           ++neighbor_id
      ) {
         if (*neighbor_id == 0 || *neighbor_id == cell_id) {
            continue;
         }
         
         SpatialCell* neighbor = mpiGrid[*neighbor_id];
         if (neighbor == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
                      << " No data for neighbor " << *neighbor_id
                      << " of cell " << cell_id
                      << endl;
            abort();
         }

         neighbor_ptrs.push_back(neighbor);         
      }
      //is threadsafe
      cell->adjust_velocity_blocks(neighbor_ptrs);
   }
   phiprof::stop("Adjusting blocks");
   phiprof::start("Set cell weight");
   // set cells' weights based on adjusted number of velocity blocks
   for (std::vector<uint64_t>::const_iterator
        cell_id = cells.begin();
        cell_id != cells.end();
        ++cell_id
   ) {
      SpatialCell* cell = mpiGrid[*cell_id];
      if (cell == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__ << " No data for spatial cell " << *cell_id << endl;
         abort();
      }

      

   }
   phiprof::stop("Set cell weight");

   return true;
}

/*
Updates velocity block lists between remote neighbors and prepares local
copies of remote neighbors for receiving velocity block data.
*/
void updateRemoteVelocityBlockLists(dccrg::Dccrg<SpatialCell>& mpiGrid)
{
   // update velocity block lists
   // Faster to do it in one operation, and not by first sending size, then list.
   phiprof::initializeTimer("Velocity block list update","MPI");
   phiprof::start("Velocity block list update");
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_SIZE_AND_LIST);
   mpiGrid.update_remote_neighbor_data(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   phiprof::stop("Velocity block list update");

   /*      
   Prepare spatial cells for receiving velocity block data
   */
   
   phiprof::start("Preparing receives");
   std::vector<uint64_t> incoming_cells
      = mpiGrid.get_list_of_remote_cells_with_local_neighbors(VLASOV_SOLVER_NEIGHBORHOOD_ID);
#pragma omp parallel for
   for(unsigned int i=0;i<incoming_cells.size();i++){
      uint64_t cell_id=incoming_cells[i];
      SpatialCell* cell = mpiGrid[cell_id];
      if (cell == NULL) {
         cerr << __FILE__ << ":" << __LINE__
              << " No data for spatial cell " << cell_id
              << endl;
         abort();
      }
      
      cell->prepare_to_receive_blocks();
   }
   phiprof::stop("Preparing receives", incoming_cells.size(), "SpatialCells");
}
