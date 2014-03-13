/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
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
#include "fieldsolver.h"
#include "projects/project.h"
#include "iowrite.h"
#include "ioread.h"


using namespace std;
using namespace phiprof;

extern Logger logFile, diagnostic;

void initVelocityGridGeometry();
void initSpatialCellCoordinates(dccrg::Dccrg<SpatialCell>& mpiGrid);
void initializeStencils(dccrg::Dccrg<SpatialCell>& mpiGrid);

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
   int neighborhood_size = 2; // At least this needed by fieldsolver. It is also fine for vlasovsolver

   mpiGrid.initialize(
      comm,
      &P::loadBalanceAlgorithm[0],
      neighborhood_size, // neighborhood size
      0, // maximum refinement level
      sysBoundaries.isBoundaryPeriodic(0),
      sysBoundaries.isBoundaryPeriodic(1),
      sysBoundaries.isBoundaryPeriodic(2)
   );
   
   initializeStencils(mpiGrid);
   
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
      vector<uint64_t> cells = mpiGrid.get_cells();
      //set background field, FIXME should be read in from restart
#pragma omp parallel for schedule(dynamic)
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         project.setCellBackgroundField(cell);
      }
      
   }
   else {
      //Initial state based on project, background field in all cells
      //and other initial values in non-sysboundary cells
      phiprof::start("Apply initial state");
      // Go through every cell on this node and initialize the 
      //    -Background field on all cells
      //  -Perturbed fields and ion distribution function in non-sysboundary cells
      //    Each initialization has to be independent to avoid threading problems 
      vector<uint64_t> cells = mpiGrid.get_cells();
#pragma omp parallel for schedule(dynamic)
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         project.setCellBackgroundField(cell);
         if(cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)
            project.setCell(cell);
      }
      //initial state for sys-boundary cells
      phiprof::stop("Apply initial state");
      phiprof::start("Apply system boundary conditions state");
      if(sysBoundaries.applyInitialState(mpiGrid, project) == false) {
         cerr << " (MAIN) ERROR: System boundary conditions initial state was not applied correctly." << endl;
         exit(1);
      }      
      phiprof::stop("Apply system boundary conditions state");
      adjustVelocityBlocks(mpiGrid); // do not initialize mover, mover has not yet been initialized here
      shrink_to_fit_grid_data(mpiGrid); //get rid of excess data already here
   }
   
   //Balance load before we transfer all data below
   balanceLoad(mpiGrid);
   
   phiprof::initializeTimer("Fetch Neighbour data","MPI");
   phiprof::start("Fetch Neighbour data");
   // update complet cell spatial data for full stencil (
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);
   mpiGrid.update_remote_neighbor_data(FULL_NEIGHBORHOOD_ID);

   phiprof::stop("Fetch Neighbour data");
   phiprof::stop("Set initial state");

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
      mpiGrid[cells[i]]->parameters[CellParams::DX  ] = mpiGrid.get_cell_length_x(cells[i]);
      mpiGrid[cells[i]]->parameters[CellParams::DY  ] = mpiGrid.get_cell_length_y(cells[i]);
      mpiGrid[cells[i]]->parameters[CellParams::DZ  ] = mpiGrid.get_cell_length_z(cells[i]);
   }
}


   


void balanceLoad(dccrg::Dccrg<SpatialCell>& mpiGrid){
// tell other processes which velocity blocks exist in remote spatial cells
   phiprof::initializeTimer("Balancing load", "Load balance");
   phiprof::start("Balancing load");

   phiprof::start("deallocate boundary data");
   //deallocate blocks in remote cells to decrease memory load
   deallocateRemoteCellBlocks(mpiGrid);
   phiprof::stop("deallocate boundary data");
   //set weights based on each cells LB weight counter
   vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i=0; i<cells.size(); ++i){
      //weight set 
      mpiGrid[cells[i]]->parameters[CellParams::LBWEIGHTCOUNTER]=
         Parameters::loadBalanceGamma +
         mpiGrid[cells[i]]->number_of_blocks * Parameters::loadBalanceAlpha + 
         mpiGrid[cells[i]]->number_of_blocks * mpiGrid[cells[i]]->number_of_blocks * Parameters::loadBalanceBeta;
      mpiGrid.set_cell_weight(cells[i], mpiGrid[cells[i]]->parameters[CellParams::LBWEIGHTCOUNTER]);
   }
   phiprof::start("dccrg.initialize_balance_load");
   mpiGrid.initialize_balance_load(true);
   phiprof::stop("dccrg.initialize_balance_load");
   
   const boost::unordered_set<uint64_t>& incoming_cells = mpiGrid.get_balance_added_cells();
   std::vector<uint64_t> incoming_cells_list (incoming_cells.begin(),incoming_cells.end()); 

   const boost::unordered_set<uint64_t>& outgoing_cells = mpiGrid.get_balance_removed_cells();
   std::vector<uint64_t> outgoing_cells_list (outgoing_cells.begin(),outgoing_cells.end()); 
   
   /*transfer cells in parts to preserve memory*/
   phiprof::start("Data transfers");
   const uint64_t num_part_transfers=5;
   for(uint64_t transfer_part=0;transfer_part<num_part_transfers;transfer_part++){
     
     //Set transfers on/off for the incming cells in this transfer set and prepare for receive
     for(unsigned int i=0;i<incoming_cells_list.size();i++){
        uint64_t cell_id=incoming_cells_list[i];
        SpatialCell* cell = mpiGrid[cell_id];
        if(cell_id%num_part_transfers!=transfer_part) {
           cell->set_mpi_transfer_enabled(false);
        }
        else{
           cell->set_mpi_transfer_enabled(true);
        }
     }
     
     //Set transfers on/off for the outgoing cells in this transfer set
     for(unsigned int i=0;i<outgoing_cells_list.size();i++){
       uint64_t cell_id=outgoing_cells_list[i];
       SpatialCell* cell = mpiGrid[cell_id];
       if(cell_id%num_part_transfers!=transfer_part) {
          cell->set_mpi_transfer_enabled(false);
       }
       else {
          cell->set_mpi_transfer_enabled(true);
       }
     }
     
     //Transfer velocity block list
     SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE1);
     mpiGrid.continue_balance_load();

     SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE2);
     mpiGrid.continue_balance_load();

     for(unsigned int i=0;i<incoming_cells_list.size();i++){
       uint64_t cell_id=incoming_cells_list[i];
       SpatialCell* cell = mpiGrid[cell_id];
       if(cell_id%num_part_transfers==transfer_part) {
          phiprof::start("Preparing receives");
          // reserve space for velocity block data in arriving remote cells
          cell->prepare_to_receive_blocks();
          phiprof::stop("Preparing receives", incoming_cells_list.size(), "Spatial cells");

       }
     }
     
     //do the actual transfer of data for the set of cells to be transferred
     phiprof::start("transfer_all_data");     
     SpatialCell::set_mpi_transfer_type(Transfer::ALL_DATA);
     mpiGrid.continue_balance_load();
     phiprof::stop("transfer_all_data");     

     //Free memory for cells that have been sent (the blockdata)
     for(unsigned int i=0;i<outgoing_cells_list.size();i++){
       uint64_t cell_id=outgoing_cells_list[i];
       SpatialCell* cell = mpiGrid[cell_id];
       if(cell_id%num_part_transfers==transfer_part) 
	 cell->clear(); //free memory of this cell as it has already been transferred. It will not be used anymore
     }
   }
   phiprof::stop("Data transfers");
   //finish up load balancing
   phiprof::start("dccrg.finish_balance_load");
   mpiGrid.finish_balance_load();
   phiprof::stop("dccrg.finish_balance_load");
 
   //Make sure transfers are enabled for all cells
   cells = mpiGrid.get_cells();
   for (uint i=0; i<cells.size(); ++i) 
      mpiGrid[cells[i]]->set_mpi_transfer_enabled(true);

   // Communicate all spatial data for FULL neighborhood, which
   // includes all data with the exception of dist function data
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);
   mpiGrid.update_remote_neighbor_data(FULL_NEIGHBORHOOD_ID);


   phiprof::start("update block lists");
   //new partition, re/initialize blocklists of remote cells.
   updateRemoteVelocityBlockLists(mpiGrid);
   phiprof::stop("update block lists");
   
   
   phiprof::start("Init solvers");
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
//ready for the new number of blocks.  Blocks exist if they have
//contents, or if their nearest neighbor in spatial or velocity space
//have content. Note that block existence does not use vlasov stencil
//as it is important to also include diagonals to avoid massloss
bool adjustVelocityBlocks(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   phiprof::initializeTimer("re-adjust blocks","Block adjustment");
   phiprof::start("re-adjust blocks");
   const vector<uint64_t> cells = mpiGrid.get_cells();
   
   phiprof::start("Compute with_content_list");
#pragma omp parallel for  
   for (uint i=0; i<cells.size(); ++i) 
     mpiGrid[cells[i]]->update_velocity_block_content_lists();
   phiprof::stop("Compute with_content_list");
   
   phiprof::initializeTimer("Transfer with_content_list","MPI");
   phiprof::start("Transfer with_content_list");
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_WITH_CONTENT_STAGE1 );
   mpiGrid.update_remote_neighbor_data(NEAREST_NEIGHBORHOOD_ID);
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_WITH_CONTENT_STAGE2 );
   mpiGrid.update_remote_neighbor_data(NEAREST_NEIGHBORHOOD_ID);
   phiprof::stop("Transfer with_content_list");

   
   //Adjusts velocity blocks in local spatial cells, doesn't adjust velocity blocks in remote cells.
   phiprof::start("Adjusting blocks");
#pragma omp parallel for
   for(unsigned int i=0;i<cells.size();i++){
     Real density_pre_adjust=0.0;
     Real density_post_adjust=0.0;
     uint64_t cell_id=cells[i];
     SpatialCell* cell = mpiGrid[cell_id];
     
     // gather spatial neighbor list and create vector with pointers to neighbor spatial cells
     const vector<uint64_t>* neighbors = mpiGrid.get_neighbors(cell_id, NEAREST_NEIGHBORHOOD_ID);
     vector<SpatialCell*> neighbor_ptrs;
     neighbor_ptrs.reserve(neighbors->size());
     for (vector<uint64_t>::const_iterator neighbor_id = neighbors->begin();
	  neighbor_id != neighbors->end(); ++neighbor_id) {
       if (*neighbor_id == 0 || *neighbor_id == cell_id) {
	 continue;
         }
       neighbor_ptrs.push_back(mpiGrid[*neighbor_id]);
     }
     if(P::sparse_conserve_mass) {
       for (unsigned int block_i = 0; block_i < cell->number_of_blocks; block_i++) {
	 Velocity_Block * block_ptr = cell->at(cell->velocity_block_list[block_i]);
	 for(unsigned int cell_i = 0;cell_i<WID3;cell_i++){
	   density_pre_adjust+=block_ptr->data[cell_i]; 
	 }
       }
     }
     cell->adjust_velocity_blocks(neighbor_ptrs);
     
     if(P::sparse_conserve_mass) {
       for (unsigned int block_i = 0; block_i < cell->number_of_blocks; block_i++) {
	 Velocity_Block * block_ptr = cell->at(cell->velocity_block_list[block_i]);
	 for(unsigned int cell_i = 0;cell_i<WID3;cell_i++){
	   density_post_adjust+=block_ptr->data[cell_i]; 
	  }
       }
       if(density_post_adjust!=0.0){
	 for (unsigned int block_i = 0; block_i < cell->number_of_blocks; block_i++) {
	   Velocity_Block * block_ptr = cell->at(cell->velocity_block_list[block_i]);
	   for(unsigned int cell_i = 0;cell_i<WID3;cell_i++){
	     block_ptr->data[cell_i]*=density_pre_adjust/density_post_adjust;                           
	   }
	 }
       }    
     }    
   }
   phiprof::stop("Adjusting blocks");

   //Updated newly adjusted velocity block lists on remote cells, and
   //prepare to receive block data
   updateRemoteVelocityBlockLists(mpiGrid);
   phiprof::stop("re-adjust blocks");
   return true;
}

/*! Shirnk to fit velocity space data to save memory.
 * \param mpiGrid Spatial grid
 */
void shrink_to_fit_grid_data(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   const std::vector<uint64_t> remote_cells = mpiGrid.get_remote_cells_on_process_boundary();

   /*append remote cells to cells*/
   cells.insert(cells.end(),remote_cells.begin(),remote_cells.end());
#pragma omp parallel for
   for(unsigned int i=0;i<cells.size();i++){
      mpiGrid[cells[i]]->shrink_to_fit();
   }
}
/*! Estimates memory consumption and writes it into logfile. Collective operation on MPI_COMM_WORLD
 * \param mpiGrid Spatial grid
 */
   


void report_memory_consumption(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   /*now report memory consumption into logfile*/
   const vector<uint64_t> cells = mpiGrid.get_cells();
   const std::vector<uint64_t> remote_cells = mpiGrid.get_remote_cells_on_process_boundary(DIST_FUNC_NEIGHBORHOOD_ID);   
   int rank,n_procs;
   MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   /*compute memory consumption of the block data, double as MPI does
    * not define proper uint64_t datatypes for MAXLOCNot Real, as we
    * want double here not to loose accuracy.
    * Computed as number of blocks * 2 arrays with block data (fx, data) *  WID3 amount of cells per block *  each cell has a size of Real
    */

   /*report data for memory needed by blocks*/
   double mem[9] = {};
   double sum_mem[9];   


   for(unsigned int i=0;i<cells.size();i++){
      // Multiplied by two because there are two block lists
      mem[0] += mpiGrid[cells[i]]->number_of_blocks * WID3 * 2 * sizeof(Realf); 
      mem[3] += (mpiGrid[cells[i]]->block_data.size() +  mpiGrid[cells[i]]->block_fx.size()) * sizeof(Realf);
      mem[6] += (mpiGrid[cells[i]]->block_data.capacity() +  mpiGrid[cells[i]]->block_fx.capacity()) * sizeof(Realf);
   }

   for(unsigned int i=0;i<remote_cells.size();i++){
      // Multiplied by two because there are two block lists
      mem[1] += mpiGrid[remote_cells[i]]->number_of_blocks * WID3 * 2 * sizeof(Realf); 
      mem[4] += (mpiGrid[remote_cells[i]]->block_data.size() +  mpiGrid[remote_cells[i]]->block_fx.size()) * sizeof(Realf);
      mem[7] += (mpiGrid[remote_cells[i]]->block_data.capacity() +  mpiGrid[remote_cells[i]]->block_fx.capacity()) * sizeof(Realf);
   }
   
   mem[2] = mem[0] + mem[1];
   mem[5] = mem[3] + mem[4];
   mem[8] = mem[6] + mem[7];

   MPI_Reduce(mem, sum_mem, 9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   logFile << "(MEM) needed for blocks: " << sum_mem[2] << endl;
   logFile << "(MEM) needed for block vectors (size): " << sum_mem[5] << endl;   
   logFile << "(MEM) needed for block vectors (capacity) " << sum_mem[8] << endl;   
   
   struct {
      double val;
      int   rank;
   } max_mem[3],mem_usage_loc[3],min_mem[3];
   for(uint i = 0; i<3; i++){
      mem_usage_loc[i].val = mem[i + 6]; //report on capacity number
      mem_usage_loc[i].rank = rank;
   }
   
   MPI_Reduce(mem_usage_loc, max_mem, 3, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
   MPI_Reduce(mem_usage_loc, min_mem, 3, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
   
   logFile << "(MEM)   Average capacity: " << sum_mem[2]/n_procs << " local cells " << sum_mem[0]/n_procs << " remote cells " << sum_mem[1]/n_procs << endl;
   logFile << "(MEM)   Max capacity:     " << max_mem[2].val   << " on  process " << max_mem[2].rank << endl;
   logFile << "(MEM)   Min capacity:     " << min_mem[2].val   << " on  process " << min_mem[2].rank << endl;
   logFile << writeVerbose;
}        

/*!      Deallocates all block data in remote cells in order to save
 *  memory
 * \param mpiGrid Spatial grid
 */
 
void deallocateRemoteCellBlocks(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   const std::vector<uint64_t> incoming_cells
      = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   for(unsigned int i=0;i<incoming_cells.size();i++){
      uint64_t cell_id=incoming_cells[i];
      SpatialCell* cell = mpiGrid[cell_id];
      if (cell != NULL) {
         cell->clear();
      }
   }

}



/*
Updates velocity block lists between remote neighbors and prepares local
copies of remote neighbors for receiving velocity block data.
*/
void updateRemoteVelocityBlockLists(dccrg::Dccrg<SpatialCell>& mpiGrid)
{
   // update velocity block lists For small velocity spaces it is
   // faster to do it in one operation, and not by first sending size,
   // then list. For large we do it in two steps
   
   phiprof::initializeTimer("Velocity block list update","MPI");
   phiprof::start("Velocity block list update");
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE1);
   mpiGrid.update_remote_neighbor_data(DIST_FUNC_NEIGHBORHOOD_ID);
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE2);
   mpiGrid.update_remote_neighbor_data(DIST_FUNC_NEIGHBORHOOD_ID);

   phiprof::stop("Velocity block list update");

   /*      
   Prepare spatial cells for receiving velocity block data
   */
   
   phiprof::start("Preparing receives");
   const std::vector<uint64_t> incoming_cells
      = mpiGrid.get_remote_cells_on_process_boundary(DIST_FUNC_NEIGHBORHOOD_ID);
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
/*
  Set stencils. These are the stencils (in 2D, real ones in 3D of
  course). x are stencil neighbor to cell local cell o:

NEAREST FIELD_SOLVER  SYSBOUNDARIES  (nearest neighbor)
-----------
  xxx
  xox
  xxx
-----------

EXTENDED_SYSBOUNDARIES (second nearest neighbor, also in diagonal)
-----------
  xxxxx
  xxxxx
  xxoxx
  xxxxx
  xxxxx
-----------  

VLASOV
-----------  
    x
    x
  xxoxx
    x
    x
-----------    

VLASOV_{XYZ}
-----------
 xxoxxx
-----------

VLASOV_SOURCE
-----------
   x
  xox
   x

-----------

DIST_FUNC  (Includes all cells which should know about each others blocks and have space for them. VLASOV + SYSBOUNDARIES.
-----------  
    x
   xxx
  xxoxx
   xxx
    x
    
-----------    

   
FULL (Includes all possible communication)
-----------
  xxxxx
  xxxxx
  xxoxx
  xxxxx
  xxxxx

-----------

SHIFT_M_X    ox
SHIFT_P_X   xo
 Y, Z in the same way
*/

void initializeStencils(dccrg::Dccrg<SpatialCell>& mpiGrid){
   
#ifdef TRANS_SEMILAG_PLM
   const int vlasov_stencil_width=1;
#endif
#if TRANS_SEMILAG_PPM
   const int vlasov_stencil_width=2;
#endif
   
   // set reduced neighborhoods
   typedef dccrg::Types<3>::neighborhood_item_t neigh_t;
   
   // set a reduced neighborhood for field solver
   std::vector<neigh_t> neighborhood;
   for (int z = -1; z <= 1; z++) {
      for (int y = -1; y <= 1; y++) {
         for (int x = -1; x <= 1; x++) {
            if (x == 0 && y == 0 && z == 0) {
               continue;
            }            
            neigh_t offsets = {{x, y, z}};
            neighborhood.push_back(offsets);
         }
      }
   }
   mpiGrid.add_remote_update_neighborhood(FIELD_SOLVER_NEIGHBORHOOD_ID, neighborhood);
   mpiGrid.add_remote_update_neighborhood(NEAREST_NEIGHBORHOOD_ID, neighborhood);
   mpiGrid.add_remote_update_neighborhood(SYSBOUNDARIES_NEIGHBORHOOD_ID, neighborhood);

   neighborhood.clear();
   for (int z = -2; z <= 2; z++) {
      for (int y = -2; y <= 2; y++) {
         for (int x = -2; x <= 2; x++) {
            if (x == 0 && y == 0 && z == 0) {
               continue;
            }
            neigh_t offsets = {{x, y, z}};
            neighborhood.push_back(offsets);
         }
      }
   }
   mpiGrid.add_remote_update_neighborhood(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID, neighborhood);

   if(vlasov_stencil_width>=3) {
      /*add face neighbors if stencil width larger than 2*/
      neighborhood.push_back({{ vlasov_stencil_width, 0, 0}});
      neighborhood.push_back({{-vlasov_stencil_width, 0, 0}});
      neighborhood.push_back({{0, vlasov_stencil_width, 0}});
      neighborhood.push_back({{0,-vlasov_stencil_width, 0}});
      neighborhood.push_back({{0, 0, vlasov_stencil_width}});
      neighborhood.push_back({{0, 0,-vlasov_stencil_width}});     
   }
   /*all possible communication pairs*/
   mpiGrid.add_remote_update_neighborhood(FULL_NEIGHBORHOOD_ID, neighborhood);

   
   /*stencils for semilagrangian propagators*/ 
   neighborhood.clear();
   for (int d = -vlasov_stencil_width; d <= vlasov_stencil_width; d++) {
     if (d != 0) {
        neighborhood.push_back({{d, 0, 0}});
        neighborhood.push_back({{0, d, 0}});
        neighborhood.push_back({{0, 0, d}});
     }
   }
   mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_NEIGHBORHOOD_ID, neighborhood);

   // add remaining nearest neighbors for DIST_FUNC neighborhood
   for (int z = -1; z <= 1; z++) {
      for (int y = -1; y <= 1; y++) {
         for (int x = -1; x <= 1; x++) {
            //do not add cells already in neighborhood (vlasov solver)
            if (x == 0 && y == 0 ) continue;
            if (x == 0 && z == 0 ) continue;
            if (y == 0 && z == 0 ) continue;
            
            neigh_t offsets = {{x, y, z}};
            neighborhood.push_back(offsets);
         }
      }
   }
   mpiGrid.add_remote_update_neighborhood(DIST_FUNC_NEIGHBORHOOD_ID, neighborhood);
   
   neighborhood.clear();
   for (int d = -vlasov_stencil_width; d <= vlasov_stencil_width; d++) {
     if (d != 0) {
        neighborhood.push_back({{d, 0, 0}});
     }
   }
   mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_X_NEIGHBORHOOD_ID, neighborhood);

   
   neighborhood.clear();
   for (int d = -vlasov_stencil_width; d <= vlasov_stencil_width; d++) {
     if (d != 0) {
        neighborhood.push_back({{0, d, 0}});
     }
   }
   mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID, neighborhood);

   
   neighborhood.clear();
   for (int d = -vlasov_stencil_width; d <= vlasov_stencil_width; d++) {
     if (d != 0) {
        neighborhood.push_back({{0, 0, d}});
     }
   }
   mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID, neighborhood);

   neighborhood.clear();
   for (int d = -1; d <= 1; d++) {
     if (d != 0) {
        neighborhood.push_back({{d, 0, 0}});
     }
   }
   mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_SOURCE_X_NEIGHBORHOOD_ID, neighborhood);

   neighborhood.clear();
   for (int d = -1; d <= 1; d++) {
     if (d != 0) {
        neighborhood.push_back({{0, d, 0}});
     }
   }
   mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_SOURCE_Y_NEIGHBORHOOD_ID, neighborhood);

   neighborhood.clear();
   for (int d = -1; d <= 1; d++) {
     if (d != 0) {
        neighborhood.push_back({{0, 0, d}});
     }
   }
   mpiGrid.add_remote_update_neighborhood(VLASOV_SOLVER_SOURCE_Z_NEIGHBORHOOD_ID, neighborhood);


   neighborhood.clear();
   neighborhood.push_back({{1, 0, 0}});
   mpiGrid.add_remote_update_neighborhood(SHIFT_M_X_NEIGHBORHOOD_ID, neighborhood);
   neighborhood.clear();
   neighborhood.push_back({{0, 1, 0}});
   mpiGrid.add_remote_update_neighborhood(SHIFT_M_Y_NEIGHBORHOOD_ID, neighborhood);
   neighborhood.clear();
   neighborhood.push_back({{0, 0, 1}});
   mpiGrid.add_remote_update_neighborhood(SHIFT_M_Z_NEIGHBORHOOD_ID, neighborhood);
   neighborhood.clear();
   neighborhood.push_back({{-1, 0, 0}});
   mpiGrid.add_remote_update_neighborhood(SHIFT_P_X_NEIGHBORHOOD_ID, neighborhood);
   neighborhood.clear();
   neighborhood.push_back({{0, -1, 0}});
   mpiGrid.add_remote_update_neighborhood(SHIFT_P_Y_NEIGHBORHOOD_ID, neighborhood);
   neighborhood.clear();
   neighborhood.push_back({{0, 0, -1}});
   mpiGrid.add_remote_update_neighborhood(SHIFT_P_Z_NEIGHBORHOOD_ID, neighborhood);

}
