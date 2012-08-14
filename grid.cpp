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

#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"
#include "transferstencil.h"

#include "vlsvwriter2.h" 
#include "fieldsolver.h"
#include "project.h"



using namespace std;
using namespace phiprof;

extern Logger logfile, diagnostic;

//Subroutine for initializing all local spatial cells
void initSpatialCells(dccrg::Dccrg<SpatialCell>& mpiGrid,SysBoundary& sysBoundaries);



//subroutine to adjust blocks of local cells; remove/add based on user-defined limits
bool adjust_local_velocity_blocks(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);



bool initializeGrid(int argn,
                    char **argc,
                    dccrg::Dccrg<SpatialCell>& mpiGrid,
                    SysBoundary& sysBoundaries
                   ){
   bool success = true;
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
      logfile << "\t ERROR: Zoltan initialization failed." << std::endl << writeVerbose;
      success = false;
   } else {
      logfile << "\t Zoltan " << zoltanVersion << " initialized successfully" << std::endl << writeVerbose;
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
      sysBoundaries.isBoundaryPeriodic(0),
      sysBoundaries.isBoundaryPeriodic(1),
      sysBoundaries.isBoundaryPeriodic(2)
      );
   
   mpiGrid.set_partitioning_option("IMBALANCE_TOL", P::loadBalanceTolerance);
   phiprof::start("Initial load-balancing");
   if (myrank == MASTER_RANK) logfile << "(INIT): Starting initial load balance." << endl << writeVerbose;
   mpiGrid.balance_load();
   phiprof::stop("Initial load-balancing");
   
   if (myrank == MASTER_RANK) logfile << "(INIT): Set initial state." << endl << writeVerbose;
   // Go through every spatial cell on this CPU, and create the initial state:
   phiprof::start("Set initial state");
   if(sysBoundaries.initializeSysBoundaries(P::t_min) == false) {
      if (myrank == MASTER_RANK) cerr << "Error in initialising the system boundaries." << endl;
      success = false;
   }
   
   initSpatialCells(mpiGrid, sysBoundaries);

   /*
   //The new structure for the initializing cells instead of calling initSpatialCells above

   initSpatialCellCoordinates();  //just initialize coordinates 
   sysBoundaries.classifyCells(mpigrid); // Classify all local cells in mpigrid (potentially a collective MPI operation)
   applyInitialState(mpigrid); // Apply initial state, will only set it for non-sysboundary cells
   sysboundaries.setState(mpigrid); // set sysboundary states for all local sysboundary cells 
   updateSparseVelocityStuff(mpiGrid) // The blockupdates done now in end of initSpatialCells

   */    

   phiprof::stop("Set initial state");

   balanceLoad(mpiGrid);
   return success;
}


void initSpatialCells(dccrg::Dccrg<SpatialCell>& mpiGrid,
                      SysBoundary& sysBoundaries
){
   typedef Parameters P;
   phiprof::start("init cell values");
   vector<uint64_t> cells = mpiGrid.get_cells();
   
   
   //  Go through every cell on this node and initialize the pointers to 
   // cpu memory, physical parameters and volume averages for each phase space 
   // point in the velocity grid. Velocity block neighbour list is also 
   // constructed here:
   // Each initialization has to be independent to avoid threading problems 
//#pragma omp parallel for schedule(dynamic)
   for (uint i=0; i<cells.size(); ++i) {
      Real xmin,ymin,zmin,dx,dy,dz;
      dx = mpiGrid.get_cell_x_size(cells[i]);
      dy = mpiGrid.get_cell_y_size(cells[i]);
      dz = mpiGrid.get_cell_z_size(cells[i]);
      xmin = mpiGrid.get_cell_x_min(cells[i]);
      ymin = mpiGrid.get_cell_y_min(cells[i]);
      zmin = mpiGrid.get_cell_z_min(cells[i]);
      initSpatialCell(*(mpiGrid[cells[i]]),sysBoundaries,xmin,ymin,zmin,dx,dy,dz,false);
   }
   phiprof::stop("init cell values");
   updateRemoteVelocityBlockLists(mpiGrid);
   
   //in principle not needed as that was done in initSpatialCell, but lets be safe and do it anyway as it does not  cost much
   for (uint i=0; i<cells.size(); ++i)
   mpiGrid[cells[i]]->update_all_block_has_content();
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_HAS_CONTENT);
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
}


/*! Set up a spatial cell.*/
bool initSpatialCell(SpatialCell& cell,
                     SysBoundary& sysBoundaries,
                     creal& xmin,creal& ymin,creal& zmin,
                     creal& dx,creal& dy,creal& dz,
                     const bool& isRemote
) {
   typedef Parameters P;
   using namespace sysboundarytype;
   // Set up cell parameters:
   cell.parameters[CellParams::XCRD] = xmin;
   cell.parameters[CellParams::YCRD] = ymin;
   cell.parameters[CellParams::ZCRD] = zmin;
   cell.parameters[CellParams::DX  ] = dx;
   cell.parameters[CellParams::DY  ] = dy;
   cell.parameters[CellParams::DZ  ] = dz;
   
   // Initialise system boundary conditions (they need the initialised positions!!)
   if(sysBoundaries.assignSysBoundaryType(cell) == false) {
      cerr << "(MAIN) ERROR: System boundary conditions were not set correctly." << endl;
      exit(1);
   }
   
   switch(cell.sysBoundaryFlag) {
      case DO_NOT_COMPUTE:
         cell.parameters[CellParams::BX] = 0.0;
         cell.parameters[CellParams::BY] = 0.0;
         cell.parameters[CellParams::BZ] = 0.0;
         cell.parameters[CellParams::EX] = 0.0;
         cell.parameters[CellParams::EY] = 0.0;
         cell.parameters[CellParams::EZ] = 0.0;
         break;
      case SW:
         calcSysBoundaryCellParameters(sysBoundaries, cell.sysBoundaryFlag, &(cell.parameters[0]), 0.0);
         break;
      default:
         calcCellParameters(&(cell.parameters[0]), 0.0);
   }
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
   creal dvx_blockCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
   creal dvy_blockCell = SpatialCell::cell_dvy; //                                vy
   creal dvz_blockCell = SpatialCell::cell_dvz; //                                vz
   
   
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
         Real average;
         switch(cell.sysBoundaryFlag) {
            case DO_NOT_COMPUTE:
               average = 0.0;
               break;
            case SW:
               average = 
               calcSysBoundaryPhaseSpaceDensity(sysBoundaries,
                                                cell.sysBoundaryFlag,
                                                xmin,ymin,zmin,
                                                dx,dy,dz,
                                                vx_cell,vy_cell,vz_cell,
                                                dvx_blockCell,dvy_blockCell,dvz_blockCell);
               break;
            // NOT_SYSBOUNDARY, OUTFLOW, IONOSPHERE should be initialised in the same way:
            default:
               average = calcPhaseSpaceDensity(xmin,ymin,zmin,
                                               dx,dy,dz,
                                               vx_cell,vy_cell,vz_cell,
                                               dvx_blockCell,dvy_blockCell,dvz_blockCell);
         }

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
      logfile << "(MAIN): Vlasov propagator did not initialize correctly!" << endl << writeVerbose;
      exit(1);
   }

      // Initialize field propagator:
   if (initializeFieldPropagatorAfterRebalance(mpiGrid) == false) {
       logfile << "(MAIN): Field propagator did not initialize correctly!" << endl << writeVerbose;
       exit(1);
   }
   phiprof::stop("Init solvers");
   
   phiprof::stop("Balancing load");
}


bool writeGrid(const dccrg::Dccrg<SpatialCell>& mpiGrid,
               DataReducer& dataReducer,
               const string& name,
               const uint& index,
               const bool& writeRestart) {
    double allStart = MPI_Wtime();
    bool success = true;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(writeRestart)
        phiprof::start("writeGrid-restart");
    else
        phiprof::start("writeGrid-reduced");
    
   // Create a name for the output file and open it with VLSVWriter:
   stringstream fname;
   fname << name <<".";
   fname.width(7);
   fname.fill('0');
   fname << index << ".vlsv";
   
   VLSVWriter vlsvWriter;
   vlsvWriter.open(fname.str(),MPI_COMM_WORLD,0);
   
   // Get all local cell IDs and write to file:
   map<string,string> attribs;
   vector<uint64_t> cells = mpiGrid.get_cells();

   if (vlsvWriter.writeArray("MESH","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) {
      cerr << "Proc #" << myrank << " failed to write cell IDs!" << endl;
   }

   // Create a buffer for spatial cell coordinates. Copy all coordinates to 
   // buffer and write:
   Real* buffer = new Real[6*cells.size()];
   for (size_t i=0; i<cells.size(); ++i) {
      SpatialCell* SC = mpiGrid[cells[i]];
      for (int j=0; j<6; ++j) {
	 buffer[6*i+j] = SC->parameters[j];
      }
   }
   if (vlsvWriter.writeArray("COORDS","SpatialGrid",attribs,cells.size(),6,buffer) == false) {
      cerr << "Proc #" << myrank << " failed to write cell coords!" << endl;
   }
   delete[] buffer;

   // Write variables calculated by DataReductionOperators (DRO). We do not know how many 
   // numbers each DRO calculates, so a buffer has to be re-allocated for each DRO:
   char* varBuffer;
   attribs.clear();
   attribs["mesh"] = "SpatialGrid";
   for (uint i=0; i<dataReducer.size(); ++i) {
      string variableName,dataType;
      uint dataSize,vectorSize;
      variableName = dataReducer.getName(i);
      if (dataReducer.getDataVectorInfo(i,dataType,dataSize,vectorSize) == false) {
	 cerr << "ERROR when requesting info from DRO " << i << endl;
      }
      uint64_t arraySize = cells.size()*vectorSize*dataSize;
      
      // Request DataReductionOperator to calculate the reduced data for all local cells:
      varBuffer = new char[arraySize];
      for (uint64_t cell=0; cell<cells.size(); ++cell) {
	 if (dataReducer.reduceData(mpiGrid[cells[cell]],i,varBuffer + cell*vectorSize*dataSize) == false) success = false;
      }
      if (success == false) logfile << "(MAIN) writeGrid: ERROR datareductionoperator '" << dataReducer.getName(i) << "' returned false!" << endl << writeVerbose;
      
      // Write reduced data to file:
      if (vlsvWriter.writeArray("VARIABLE",variableName,attribs,cells.size(),vectorSize,dataType,dataSize,varBuffer) == false) success = false;
      if (success == false) logfile << "(MAIN) writeGrid: ERROR failed to write datareductionoperator data to file!" << endl << writeVerbose;
      delete[] varBuffer;
      varBuffer = NULL;
   }

   // If restart data is not written, exit here:
   if (writeRestart == false) {
      phiprof::initializeTimer("Barrier","MPI","Barrier");
      phiprof::start("Barrier");
      MPI_Barrier(MPI_COMM_WORLD);
      phiprof::stop("Barrier");
      vlsvWriter.close();
      phiprof::stop("writeGrid-reduced");
      return success;
   }

   attribs.clear();
//START TO WRITE RESTART
   
   // Write spatial cell parameters:
   Real* paramsBuffer = new Real[cells.size()*CellParams::N_SPATIAL_CELL_PARAMS];

   for (size_t i = 0; i < cells.size(); ++i)
      for (uint j = 0; j < CellParams::N_SPATIAL_CELL_PARAMS; ++j) {
      paramsBuffer[i*CellParams::N_SPATIAL_CELL_PARAMS+j] = mpiGrid[cells[i]]->parameters[j];
   }
   
   if (vlsvWriter.writeArray("CELLPARAMS","SpatialGrid",attribs,cells.size(),CellParams::N_SPATIAL_CELL_PARAMS,paramsBuffer) == false) {
      logfile << "(MAIN) writeGrid: ERROR failed to write spatial cell parameters!" << endl << writeVerbose;
      success = false;
   }
   delete[] paramsBuffer;
   
   //REMOVED-did not write out real data Write the number of spatial neighbours each cell has 

   
   // Write velocity blocks and related data. Which cells write velocity grids 
   // should be requested from a function, but for now we just write velocity grids for all cells.
   // First write global IDs of those cells which write velocity blocks (here: all cells):
   if (vlsvWriter.writeArray("CELLSWITHBLOCKS","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) success = false;
   if (success == false) logfile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;
   
   // Write the number of velocity blocks in each spatial cell. Again a temporary buffer is used:
   uint* N_blocks = new uint[cells.size()];
   uint64_t totalBlocks = 0;
   for (size_t cell=0; cell<cells.size(); ++cell) {
      N_blocks[cell] = mpiGrid[cells[cell]]->size();
      totalBlocks += mpiGrid[cells[cell]]->size();
   }
   if (vlsvWriter.writeArray("NBLOCKS","SpatialGrid",attribs,cells.size(),1,N_blocks) == false) success = false;
   if (success == false) logfile << "(MAIN) writeGrid: ERROR failed to write NBLOCKS to file!" << endl << writeVerbose;
   delete[] N_blocks;
   
   double start = MPI_Wtime();

   // Write velocity block coordinates.// TODO: add support for MPI_Datatype in startMultiwrite... or use normal writeArray as all data is collected already in one place

   std::vector<Real> velocityBlockParameters;
   velocityBlockParameters.reserve(totalBlocks*BlockParams::N_VELOCITY_BLOCK_PARAMS);

   // gather data for writing
   for (size_t cell=0; cell<cells.size(); ++cell) {
      int index=0;
      SpatialCell* SC = mpiGrid[cells[cell]];
      for (unsigned int block_i=0;block_i < SC->number_of_blocks;block_i++){
         unsigned int block = SC->velocity_block_list[block_i];
         Velocity_Block* block_data = SC->at(block);
         for(unsigned int p=0;p<BlockParams::N_VELOCITY_BLOCK_PARAMS;++p){
            velocityBlockParameters.push_back(block_data->parameters[p]);
         }
      }
      
   }

   if (vlsvWriter.writeArray("BLOCKCOORDINATES","SpatialGrid",attribs,totalBlocks,BlockParams::N_VELOCITY_BLOCK_PARAMS,&(velocityBlockParameters[0])) == false) success = false;
   if (success == false) logfile << "(MAIN) writeGrid: ERROR failed to write BLOCKCOORDINATES to file!" << endl << writeVerbose;
   velocityBlockParameters.clear();

   
   // Write values of distribution function:
   std::vector<Real> velocityBlockData;
   velocityBlockData.reserve(totalBlocks*SIZE_VELBLOCK);
   
   for (size_t cell=0; cell<cells.size(); ++cell) {
      int index=0;
      SpatialCell* SC = mpiGrid[cells[cell]];
      for (unsigned int block_i=0;block_i < SC->number_of_blocks;block_i++){
         unsigned int block = SC->velocity_block_list[block_i];
         Velocity_Block* block_data = SC->at(block);
         for(unsigned int vc=0;vc<SIZE_VELBLOCK;++vc){
            velocityBlockData.push_back(block_data->data[vc]);
         }
      }
   }

   attribs["mesh"] = "SpatialGrid";
   if (vlsvWriter.writeArray("BLOCKVARIABLE","avgs",attribs,totalBlocks,SIZE_VELBLOCK,&(velocityBlockData[0])) == false) success=false;
   if (success ==false)      logfile << "(MAIN) writeGrid: ERROR occurred when writing BLOCKVARIABLE avgs" << endl << writeVerbose;
   velocityBlockData.clear();
   
   double end = MPI_Wtime();


   vlsvWriter.close();

   double allEnd = MPI_Wtime();
   
   //double bytesWritten = totalBlocks*((SIZE_VELBLOCK+SIZE_BLOCKPARAMS)*sizeof(Real)+(SIZE_NBRS_VEL)*sizeof(uint));
   double secs = end-start;
   //FIXME, should be global info
   //logfile << "Wrote " << bytesWritten/1.0e6 << " MB of data in " << secs << " seconds, datarate is " << bytesWritten/secs/1.0e9 << " GB/s" << endl << writeVerbose;
   
   double allSecs = allEnd-allStart;

   phiprof::stop("writeGrid-restart");//,1.0e-6*bytesWritten,"MB");
   return success;
}


bool computeDiagnostic(const dccrg::Dccrg<SpatialCell>& mpiGrid,
		       DataReducer& dataReducer,
		       luint tstep)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   
   string dataType;
   uint dataSize, vectorSize;
   vector<uint64_t> cells = mpiGrid.get_cells();
   cuint nCells = cells.size();
   cuint nOps = dataReducer.size();
   vector<Real> localMin(nOps), localMax(nOps), localSum(nOps+1), localAvg(nOps),
               globalMin(nOps),globalMax(nOps),globalSum(nOps+1),globalAvg(nOps);
   localSum[0] = 1.0 * nCells;
   Real buffer;
   bool success = true;
   static bool printDiagnosticHeader = true;
   
   if (printDiagnosticHeader == true && myrank == 0) {
      diagnostic << "# Column 1 Step" << endl;
      for (uint i=0; i<nOps; ++i) {
	 diagnostic << "# Columns " << 2 + i*4 << " to " << 5 + i*4 << ": " << dataReducer.getName(i) << " min max sum average" << endl;
      }
      printDiagnosticHeader = false;
   }
   
   for (uint i=0; i<nOps; ++i) {
      
      if (dataReducer.getDataVectorInfo(i,dataType,dataSize,vectorSize) == false) {
	 cerr << "ERROR when requesting info from diagnostic DRO " << dataReducer.getName(i) << endl;
      }
      localMin[i] = std::numeric_limits<Real>::max();
      localMax[i] = std::numeric_limits<Real>::min();
      localSum[i+1] = 0.0;
      buffer = 0.0;
      
      // Request DataReductionOperator to calculate the reduced data for all local cells:
      for (uint64_t cell=0; cell<nCells; ++cell) {
	 success = true;
	 if (dataReducer.reduceData(mpiGrid[cells[cell]], i, &buffer) == false) success = false;
	 localMin[i] = min(buffer, localMin[i]);
	 localMax[i] = max(buffer, localMax[i]);
	 localSum[i+1] += buffer;
      }
      localAvg[i] = localSum[i+1];
      
      if (success == false) logfile << "(MAIN) computeDiagnostic: ERROR datareductionoperator '" << dataReducer.getName(i) << "' returned false!" << endl << writeVerbose;
   }
   
   MPI_Reduce(&localMin[0], &globalMin[0], nOps, MPI_Type<Real>(), MPI_MIN, 0, MPI_COMM_WORLD);
   MPI_Reduce(&localMax[0], &globalMax[0], nOps, MPI_Type<Real>(), MPI_MAX, 0, MPI_COMM_WORLD);
   MPI_Reduce(&localSum[0], &globalSum[0], nOps + 1, MPI_Type<Real>(), MPI_SUM, 0, MPI_COMM_WORLD);

   diagnostic << setprecision(12); 
   diagnostic << tstep << "\t";
   
   for (uint i=0; i<nOps; ++i) {
      if (globalSum[0] != 0.0) globalAvg[i] = globalSum[i+1] / globalSum[0];
      else globalAvg[i] = globalSum[i+1];
      if (myrank == 0) {
	 diagnostic << globalMin[i] << "\t" <<
	 globalMax[i] << "\t" <<
	 globalSum[i+1] << "\t" <<
	 globalAvg[i] << "\t";
      }
   }
   if (myrank == 0) diagnostic << endl << write;
   return true;
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
