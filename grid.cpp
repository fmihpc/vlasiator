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

#include "boost/mpi.hpp"
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

#include "transferstencil.h"

#include "fieldsolver.h"
#include "project.h"



using namespace std;
using namespace phiprof;

extern Logger logFile, diagnostic;

//Subroutine for initializing all local spatial cells
void initSpatialCells(dccrg::Dccrg<SpatialCell>& mpiGrid);

void initSpatialCellCoordinates(dccrg::Dccrg<SpatialCell>& mpiGrid);

//subroutine to adjust blocks of local cells; remove/add based on user-defined limits
bool adjust_local_velocity_blocks(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);



bool initializeGrid(int argn, char **argc,dccrg::Dccrg<SpatialCell>& mpiGrid){
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

   // initialize velocity grid of spatial cells before creating cells in dccrg.initialize
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
   spatial_cell::SpatialCell::velocity_block_min_avg_value = P::sparseMinAvgValue;
 
   
   
// Init Zoltan:
   float zoltanVersion;
   if (Zoltan_Initialize(argn,argc,&zoltanVersion) != ZOLTAN_OK) {
      logFile << "\t ERROR: Zoltan initialization failed, aborting." << std::endl << writeVerbose;
   } else {
      logFile << "\t Zoltan " << zoltanVersion << " initialized successfully" << std::endl << writeVerbose;
   }
   
   mpiGrid.set_geometry(
      P::xcells_ini, P::ycells_ini, P::zcells_ini,
      P::xmin, P::ymin, P::zmin,
      P::dx_ini, P::dy_ini, P::dz_ini
   );
   boost::mpi::communicator comm;//corresponds to MPI_COMM_WORLD
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
      P::periodic_x, P::periodic_y, P::periodic_z
      );
   
   
   
   
   mpiGrid.set_partitioning_option("IMBALANCE_TOL", P::loadBalanceTolerance);
   phiprof::start("Initial load-balancing");
   if (myrank == 0) logFile << "(INIT): Starting initial load balance." << endl << writeVerbose;
   mpiGrid.balance_load();
   phiprof::stop("Initial load-balancing");


   if (myrank == 0) logFile << "(INIT): Set initial state." << endl << writeVerbose;
   
   // Go through every spatial cell on this CPU, and create the initial state:


   phiprof::start("Set initial state");
   if (P::restartFileName==string("")){
      initSpatialCells(mpiGrid);
   }
   else{
      if (myrank == 0) logFile << "(INIT): Reading in state from restartfile" << endl << writeVerbose;
      readGrid(mpiGrid,P::restartFileName);
   }

//now we make sure all remote cells are set up to receive block data and that all sparse velocity blocks are up-to-date
   vector<uint64_t> cells = mpiGrid.get_cells();
   updateRemoteVelocityBlockLists(mpiGrid);
   
   //in principle not needed as that was done in initSpatialCell, but lets be safe and do it anyway as it does not  cost much
   for (uint i=0; i<cells.size(); ++i) 
      mpiGrid[cells[i]]->update_all_block_has_content();     
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_HAS_CONTENT );
    mpiGrid.update_remote_neighbor_data();
    
    adjust_local_velocity_blocks(mpiGrid);
    //velocity blocks adjusted, lets prepare again for new lists
    updateRemoteVelocityBlockLists(mpiGrid);

    phiprof::initializeTimer("Fetch Neighbour data","MPI");
    phiprof::start("Fetch Neighbour data");
    // update complete spatial cell data 
    SpatialCell::set_mpi_transfer_type(Transfer::ALL_DATA);
    mpiGrid.update_remote_neighbor_data();       
    phiprof::stop("Fetch Neighbour data");   


   
   phiprof::stop("Set initial state");


   
   balanceLoad(mpiGrid);
   return true;
}




void initSpatialCells(dccrg::Dccrg<SpatialCell>& mpiGrid){
    typedef Parameters P;
    phiprof::start("init cell values");
    vector<uint64_t> cells = mpiGrid.get_cells();


    //  Go through every cell on this node and initialize the pointers to 
    // cpu memory, physical parameters and volume averages for each phase space 
    // point in the velocity grid. Velocity block neighbour list is also 
    // constructed here:
    // Each initialization has to be independent to avoid threadingproblems 
#pragma omp parallel for schedule(dynamic)
    for (uint i=0; i<cells.size(); ++i) {
       Real xmin,ymin,zmin,dx,dy,dz;
       dx = mpiGrid.get_cell_x_size(cells[i]);
       dy = mpiGrid.get_cell_y_size(cells[i]);
       dz = mpiGrid.get_cell_z_size(cells[i]);
       xmin = mpiGrid.get_cell_x_min(cells[i]);
       ymin = mpiGrid.get_cell_y_min(cells[i]);
       zmin = mpiGrid.get_cell_z_min(cells[i]);
       initSpatialCell(*(mpiGrid[cells[i]]),xmin,ymin,zmin,dx,dy,dz,false);
    }
    phiprof::stop("init cell values");
}


/** Set up a spatial cell.
 */
bool initSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,
		      creal& zmin,creal& dx,creal& dy,creal& dz,
		     const bool& isRemote) {
   typedef Parameters P;
   // Set up cell parameters:
   cell.parameters[CellParams::XCRD] = xmin;
   cell.parameters[CellParams::YCRD] = ymin;
   cell.parameters[CellParams::ZCRD] = zmin;
   cell.parameters[CellParams::DX  ] = dx;
   cell.parameters[CellParams::DY  ] = dy;
   cell.parameters[CellParams::DZ  ] = dz;
   calcCellParameters(&(cell.parameters[0]),0.0);
   cell.parameters[CellParams::RHO  ] = 0.0;
   cell.parameters[CellParams::RHOVX] = 0.0;
   cell.parameters[CellParams::RHOVY] = 0.0;
   cell.parameters[CellParams::RHOVZ] = 0.0;
   cell.parameters[CellParams::RHOLOSSADJUST] = 0.0;
   cell.parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;

   // Go through each velocity block in the velocity phase space grid.
   // Set the initial state and block parameters:
   
   creal dvx_block = SpatialCell::block_dvx; // Size of a block in vx-direction
   creal dvy_block = SpatialCell::block_dvy; //                    vy
   creal dvz_block = SpatialCell::block_dvz; //                    vz
   creal dvx_blockCell = SpatialCell::cell_dvx;                 // Size of one cell in a block in vx-direction
   creal dvy_blockCell = SpatialCell::cell_dvy;                 //                                vy
   creal dvz_blockCell = SpatialCell::cell_dvz;                 //                                vz

   
   for (uint kv=0; kv<P::vzblocks_ini; ++kv) for (uint jv=0; jv<P::vyblocks_ini; ++jv) for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
      creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
      creal vy_block = P::vymin + jv*dvy_block; // vy-
      creal vz_block = P::vzmin + kv*dvz_block; // vz-

      if (isRemote == true) continue;
      // Calculate volume average of distrib. function for each cell in the block.
      for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
	 creal vx_cell = vx_block + ic*dvx_blockCell;
	 creal vy_cell = vy_block + jc*dvy_blockCell;
	 creal vz_cell = vz_block + kc*dvz_blockCell;
         Real average=calcPhaseSpaceDensity(xmin,ymin,zmin,dx,dy,dz,vx_cell,vy_cell,vz_cell,dvx_blockCell,dvy_blockCell,dvz_blockCell);

         if(average!=0.0){
            creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
            creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
            creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
            cell.set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
            // Add contributions to spatial cell velocity moments:
            creal dV = dvx_blockCell*dvy_blockCell*dvz_blockCell;  // Volume of one cell in a block      
            cell.parameters[CellParams::RHO  ] += average*dV;
            cell.parameters[CellParams::RHOVX] += average*vx_cell_center*dV;
            cell.parameters[CellParams::RHOVY] += average*vy_cell_center*dV;
            cell.parameters[CellParams::RHOVZ] += average*vz_cell_center*dV;
         }
      }
   }
   creal spatialVolume = cell.parameters[CellParams::DX]*cell.parameters[CellParams::DY]*cell.parameters[CellParams::DZ];
   cell.parameters[CellParams::RHO  ] /= spatialVolume;
   cell.parameters[CellParams::RHOVX] /= spatialVolume;
   cell.parameters[CellParams::RHOVY] /= spatialVolume;
   cell.parameters[CellParams::RHOVZ] /= spatialVolume;

   //lets get rid of blocks not fulfilling the criteria here to save
   //memory.  neighbor_ptrs is empty as we do not have any consistent
   //data in neighbours yet, adjustments done only based on velocity
   //space.
   vector<SpatialCell*> neighbor_ptrs;
   cell.update_all_block_has_content();
   cell.adjust_velocity_blocks(neighbor_ptrs);
   return true;
}


void balanceLoad(dccrg::Dccrg<SpatialCell>& mpiGrid){

   //set weights based on each cells LB weight counter
   vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i=0; i<cells.size(); ++i){
      if(P::maxAccelerationSubsteps!=1) {
         //use time-metric from solvers
         mpiGrid.set_cell_weight(cells[i], mpiGrid[cells[i]]->parameters[CellParams::LBWEIGHTCOUNTER]);
      }
      else{
         //No substepping in acceleration step, use number of blocks instead as metric as that provides slightly better results
         mpiGrid.set_cell_weight(cells[i], mpiGrid[cells[i]]->number_of_blocks);
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

bool adjustVelocityBlocks(dccrg::Dccrg<SpatialCell>& mpiGrid) {
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
   mpiGrid.update_remote_neighbor_data();
   phiprof::stop("Transfer block data");
   
   adjust_local_velocity_blocks(mpiGrid);

   updateRemoteVelocityBlockLists(mpiGrid);
   //re-init vlasovmover
   phiprof::start("InitMoverAfterBlockChange");
   initMoverAfterBlockChange(mpiGrid);
   phiprof::stop("InitMoverAfterBlockChange");
   
   phiprof::stop("re-adjust blocks");
   return true;
}


/*!
Adjusts velocity blocks in local spatial cells.

Doesn't adjust velocity blocks of copies of remote neighbors.
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
      const vector<uint64_t>* neighbors = mpiGrid.get_neighbors(cell_id);
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
   mpiGrid.update_remote_neighbor_data();
   phiprof::stop("Velocity block list update");

   /*      
   Prepare spatial cells for receiving velocity block data
   */
   
   phiprof::start("Preparing receives");
   std::vector<uint64_t> incoming_cells = mpiGrid.get_list_of_remote_cells_with_local_neighbors();
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
