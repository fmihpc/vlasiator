#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>

#include <dccrg.hpp>
#include <boost/mpi.hpp>
#include <zoltan.h>

#include "definitions.h"
#include "logger.h"
#include "parameters.h"
#include "grid.h"
#include "writevars.h"
#include "silowriter.h"

#ifndef NOCUDA
  #include "devicegrid.h"
  #include "cudafuncs.h"
#endif

using namespace std;

#ifndef NOCUDA
  extern void buildBaseGrid(Grid& grid,DeviceGrid& deviceGrid);
  extern bool acceleration(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
  extern bool translation_step1(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
  extern bool translation_step2(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
  extern bool translation_step3(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
#else
  extern void buildBaseGrid(Grid& grid);
#endif

extern void buildSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,creal& zmin,creal& dx,creal& dy,creal& dz);
extern bool cpu_acceleration(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT);
extern bool cpu_translation1(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT);
extern bool cpu_translation2(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT);
extern bool cpu_translation3(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT);

extern bool cpu_acceleration(SpatialCell& cell);
extern bool cpu_translation1(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation2(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);
extern bool cpu_translation3(SpatialCell& cell,const std::vector<const SpatialCell*>& nbrPtrs);

Logger logger;
Grid grid;

void initSpatialCells(dccrg<SpatialCell>& mpiGrid,boost::mpi::communicator& comm) {
   typedef Parameters P;
   
   // This can be replaced by an iterator.
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   // Go through every cell on this node and initialize the pointers to 
   // cpu memory, physical parameters and volume averages for each phase space 
   // point in the velocity grid. Velocity block neighbour list is also 
   // constructed here:
   Real xmin,ymin,zmin,dx,dy,dz;
   for (uint i=0; i<cells.size(); ++i) {
      dx = mpiGrid.get_cell_x_size(cells[i]);
      dy = mpiGrid.get_cell_y_size(cells[i]);
      dz = mpiGrid.get_cell_z_size(cells[i]);
      xmin = mpiGrid.get_cell_x(cells[i]);
      ymin = mpiGrid.get_cell_y(cells[i]);
      zmin = mpiGrid.get_cell_z(cells[i]);
      xmin -= 0.5*dx;
      ymin -= 0.5*dy;
      zmin -= 0.5*dz;

      //mpiGrid[cells[i]]->cpuIndex = cells[i];
      mpiGrid[cells[i]]->allocateMemory(P::vxblocks_ini*P::vyblocks_ini*P::vzblocks_ini);
      buildSpatialCell(*(mpiGrid[cells[i]]),xmin,ymin,zmin,dx,dy,dz);
   }
}

void writeSpatialCells(const boost::mpi::communicator& comm,dccrg<SpatialCell>& mpiGrid) {
   // This can be replaced by an iterator.
   std::vector<uint64_t> cells = mpiGrid.get_cells();

   std::stringstream fname;
   fname << "cells." << comm.rank() << '.';
   fname.width(5);
   fname.fill('0');
   fname << Parameters::tstep << ".silo";
   
   openOutputFile(fname.str(),"spatial_cells");
   reserveSpatialCells(cells.size());
   for (uint i=0; i<cells.size(); ++i) {
      Real* const avgs = mpiGrid[cells[i]]->cpu_avgs;
      if (avgs == NULL) {
	 std::cerr << "(MAIN) ERROR expected a pointer, got NULL" << std::endl;
	 continue;
      }
      Real n = 0.0;
      for (uint b=0; b<SIZE_VELBLOCK*mpiGrid[cells[i]]->N_blocks; ++b) n += avgs[b];
      addSpatialCell(mpiGrid[cells[i]]->cpu_cellParams,n);
   }
   writeSpatialCells("spatcells","n");
   closeOutputFile();
   freeCells();
}

void writeVelocityBlocks(const boost::mpi::communicator& comm,dccrg<SpatialCell>& mpiGrid) {
   // This can be replaced by an iterator.
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   
   std::stringstream fname;
   fname << "blocks." << comm.rank() << '.';
   fname.width(5);
   fname.fill('0');
   fname << Parameters::tstep << ".silo";
   
   // Write velocity grid
   openOutputFile(fname.str(),"vel_blocks");
   SpatialCell cell;
   //cell = *(mpiGrid[cells[0]]);
   cell.clone(*(mpiGrid[cells[0]]));
   
   reserveVelocityBlocks(cell.N_blocks);
   for (uint b=0; b<cell.N_blocks; ++b) {
      addVelocityGridBlock3D(cell.cpu_blockParams+b*SIZE_BLOCKPARAMS);
   }
   writeReservedBlocks("velgrid");
   freeBlocks();
   
   // Integrate phase space densities 
   for (uint i=0; i<cell.N_blocks*SIZE_VELBLOCK; ++i) cell.cpu_avgs[i] = 0.0;
   for (uint i=0; i<cell.N_blocks*SIZE_FLUXS; ++i) cell.cpu_fx[i] = 0.0;
   for (uint i=0; i<cell.N_blocks*SIZE_FLUXS; ++i) cell.cpu_fy[i] = 0.0;
   for (uint i=0; i<cell.N_blocks*SIZE_FLUXS; ++i) cell.cpu_fz[i] = 0.0;
   for (uint i=0; i<cell.N_blocks*SIZE_DERIV; ++i) cell.cpu_d1x[i] = 0.0;
   for (uint i=0; i<cell.N_blocks*SIZE_DERIV; ++i) cell.cpu_d1y[i] = 0.0;
   for (uint i=0; i<cell.N_blocks*SIZE_DERIV; ++i) cell.cpu_d1z[i] = 0.0;
   for (uint i=0; i<cell.N_blocks*SIZE_DERIV; ++i) cell.cpu_d2x[i] = 0.0;
   for (uint i=0; i<cell.N_blocks*SIZE_DERIV; ++i) cell.cpu_d2y[i] = 0.0;
   for (uint i=0; i<cell.N_blocks*SIZE_DERIV; ++i) cell.cpu_d2z[i] = 0.0;
   
   for (uint i=0; i<cells.size(); ++i) {
      creal* const avgs = mpiGrid[cells[i]]->cpu_avgs;
      if (avgs == NULL) continue;
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_avgs[j] += (mpiGrid[cells[i]])->cpu_avgs[j];
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_fx[j]   += (mpiGrid[cells[i]])->cpu_fx[j];
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_fy[j]   += (mpiGrid[cells[i]])->cpu_fy[j];
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_fz[j]   += (mpiGrid[cells[i]])->cpu_fz[j];
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_d1x[j]  += (mpiGrid[cells[i]])->cpu_d1x[j];
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_d1y[j]  += (mpiGrid[cells[i]])->cpu_d1y[j];
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_d1z[j]  += (mpiGrid[cells[i]])->cpu_d1z[j];
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_d2x[j]  += (mpiGrid[cells[i]])->cpu_d2x[j];
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_d2y[j]  += (mpiGrid[cells[i]])->cpu_d2y[j];
      for (uint j=0; j<cell.N_blocks*SIZE_VELBLOCK; ++j) cell.cpu_d2z[j]  += (mpiGrid[cells[i]])->cpu_d2z[j];
   }

   writeVelocityBlockGridScalar3D("f","velgrid",cell.N_blocks,cell.cpu_avgs);
   writeVelocityBlockGridScalar3D("Fx","velgrid",cell.N_blocks,cell.cpu_fx);
   writeVelocityBlockGridScalar3D("Fy","velgrid",cell.N_blocks,cell.cpu_fy);
   writeVelocityBlockGridScalar3D("Fz","velgrid",cell.N_blocks,cell.cpu_fz);
   writeVelocityBlockGridScalar3D("D1x","velgrid",cell.N_blocks,cell.cpu_d1x);
   writeVelocityBlockGridScalar3D("D1y","velgrid",cell.N_blocks,cell.cpu_d1y);
   writeVelocityBlockGridScalar3D("D1z","velgrid",cell.N_blocks,cell.cpu_d1z);
   writeVelocityBlockGridScalar3D("D2x","velgrid",cell.N_blocks,cell.cpu_d2x);
   writeVelocityBlockGridScalar3D("D2y","velgrid",cell.N_blocks,cell.cpu_d2y);
   writeVelocityBlockGridScalar3D("D2z","velgrid",cell.N_blocks,cell.cpu_d2z);
   closeOutputFile();
}

bool findNeighbours(std::vector<const SpatialCell*>& nbrPtr,dccrg<SpatialCell>& mpiGrid,const uint64_t& cell) {
   std::vector<uint64_t> nbrs;
   // Get a pointer to each neighbour. If the neighbour does not exists, insert a NULL ptr.
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
   
   // Periodic boundary conditions in z: FIX THIS!
   nbrPtr[4] = mpiGrid[cell];
   nbrPtr[5] = mpiGrid[cell];
   return true;
}

int main(int argn,char* args[]) {
   typedef Parameters P;
   boost::mpi::environment env(argn,args);
   boost::mpi::communicator comm;
     {
	std::stringstream ss;
	ss << "logfile." << comm.rank() << ".txt";
	logger.setOutputFile(ss.str());
     }
   
   logger << "(MAIN): Starting up." << std::endl;
   Parameters parameters;

   // Allocate memory for grids 
   #ifndef NOCUDA
   DeviceGrid deviceGrid;
   #endif

   // Check that requested memory was allocated:
   bool success = true;
   #ifndef NOCUDA
   if (deviceGrid.isInitialized() == false) {
      logger << "\t ERROR: DeviceGrid initialization failed." << std::endl;
      success = false;
   }
   #endif
      
   // Create parallel MPI grid and init Zoltan:
   float zoltanVersion;
   if (Zoltan_Initialize(argn,args,&zoltanVersion) != ZOLTAN_OK) {
      logger << "\t ERROR: Zoltan initialization failed, aborting." << std::endl;
      success = false;
   } else {
      logger << "\t Zoltan initialized successfully" << std::endl;
   }
   dccrg<SpatialCell> mpiGrid(comm,"RCB",P::xmin,P::ymin,P::zmin,P::dx_ini,P::xcells_ini,P::ycells_ini,P::zcells_ini,1,0);
   
   // If initialization was not successful, abort.
   if (success == false) {
      std::cerr << "An error has occurred, aborting. See logfile for details." << std::endl;
      logger << "Aborting" << std::endl;
      return 1;
   }

   // Do initial load balancing:
   mpiGrid.balance_load();
   comm.barrier();
   
   // Go through every spatial cell on this CPU, and create the initial state:
   initSpatialCells(mpiGrid,comm);
   comm.barrier();
   
   // Fetch neighbour data:
   mpiGrid.start_remote_neighbour_data_update(); // TEST
   mpiGrid.wait_neighbour_data_update();
   comm.barrier();
   
   // Write initial state:
   writeSpatialCells(comm,mpiGrid);
   writeVelocityBlocks(comm,mpiGrid);
   comm.barrier();

   /*
   // Create initial state to CPU memory:
   #ifndef NOCUDA
   buildBaseGrid(grid,deviceGrid);
   #else
   buildBaseGrid(grid);
   #endif
   
   #ifndef NOCUDA
     // Copy volume averages, spatial cell parameters, velocity grid block parameters, 
     // and neighbour lists to device:
     for (uint i=0; i<grid.size(); ++i) {
	if (grid[i]->devSync(Cell::Blocks,Cell::CpuToDev) == false) cerr << "Error while copying blockArray" << endl;
	if (grid[i]->devSync(Cell::BlockParams,Cell::CpuToDev) == false) cerr << "Error while copying blockParams" << endl;
	if (grid[i]->devSync(Cell::CellParams,Cell::CpuToDev) == false) cerr << "Error while copying cellParams" << endl;
	if (grid[i]->devSync(Cell::NbrsVel,Cell::CpuToDev) == false) cerr << "Error while copying nbrsVel" << endl;
     }
   #endif
    */

   // Main simulation loop:
   logger << "(MAIN): Starting main simulation loop." << std::endl;
   logger << "\t (TIME) reports value " << std::time(NULL) << std::endl;
   for (uint tstep=0; tstep < P::tsteps; ++tstep) {
      /*
      #ifndef NOCUDA
        acceleration(grid,deviceGrid,0.025);
        translation_step1(grid,deviceGrid,0.025);
        translation_step2(grid,deviceGrid,0.025);
        translation_step3(grid,deviceGrid,0.025);
      #else
        for (uint i=0; i<grid.size(); ++i) cpu_acceleration(i,*(grid[i]),grid,0.025);
        for (uint i=0; i<grid.size(); ++i) cpu_translation1(i,*(grid[i]),grid,0.025);
        for (uint i=0; i<grid.size(); ++i) cpu_translation2(i,*(grid[i]),grid,0.025);
        for (uint i=0; i<grid.size(); ++i) cpu_translation3(i,*(grid[i]),grid,0.025);
      #endif
      */
      std::vector<uint64_t> cells;
      std::vector<const SpatialCell*> nbrPtrs(6,NULL);
      SpatialCell* cellPtr;
      
      // Calculate acceleration for each cell on this process:
      cells = mpiGrid.get_cells();
      for (size_t c=0; c<cells.size(); ++c) {
	 cellPtr = mpiGrid[cells[c]];
	 if (cellPtr != NULL) cpu_acceleration(*cellPtr);
      }

      // Calculate derivatives in spatial coordinates:
      P::transmit = Transmit::AVGS;
      mpiGrid.start_remote_neighbour_data_update();
      cells = mpiGrid.get_cells_with_local_neighbours();
      for (size_t c=0; c<cells.size(); ++c) {
	 cellPtr = mpiGrid[cells[c]];
	 if (findNeighbours(nbrPtrs,mpiGrid,cells[c]) == false) 
	   {std::cerr << "Failed to find neighbours." << std::endl; continue;}
	 if (cellPtr != NULL) cpu_translation1(*cellPtr,nbrPtrs);
      }
      mpiGrid.wait_neighbour_data_update();
      cells = mpiGrid.get_cells_with_remote_neighbour();
      for (size_t c=0; c<cells.size(); ++c) {
	 cellPtr = mpiGrid[cells[c]];
	 if (findNeighbours(nbrPtrs,mpiGrid,cells[c]) == false) 
	   {std::cerr << "Failed to find neighbours." << std::endl; continue;}
	 if (cellPtr != NULL) cpu_translation1(*cellPtr,nbrPtrs);
      }
      
      // Calculate fluxes in spatial coordinates:
      P::transmit = Transmit::AVGS | Transmit::DERIV1 | Transmit::DERIV2;
      mpiGrid.start_remote_neighbour_data_update();
      cells = mpiGrid.get_cells_with_local_neighbours();
      for (size_t c=0; c<cells.size(); ++c) {
	 cellPtr = mpiGrid[cells[c]];
	 if (findNeighbours(nbrPtrs,mpiGrid,cells[c]) == false)
	   {std::cerr << "Failed to find neighbours." << std::endl; continue;}
	 if (cellPtr != NULL) cpu_translation2(*cellPtr,nbrPtrs);
      }
      mpiGrid.wait_neighbour_data_update();
      cells = mpiGrid.get_cells_with_remote_neighbour();
      for (size_t c=0; c<cells.size(); ++c) {
	 cellPtr = mpiGrid[cells[c]];
	 if (findNeighbours(nbrPtrs,mpiGrid,cells[c]) == false)
	   {std::cerr << "Failed to find neighbours." << std::endl; continue;}
	 if (cellPtr != NULL) cpu_translation2(*cellPtr,nbrPtrs);
      }

      // Calculate propagation in spatial coordinates:
      P::transmit = Transmit::FLUXES;
      mpiGrid.start_remote_neighbour_data_update();
      cells = mpiGrid.get_cells_with_local_neighbours();
      for (size_t c=0; c<cells.size(); ++c) {
	 cellPtr = mpiGrid[cells[c]];
	 if (findNeighbours(nbrPtrs,mpiGrid,cells[c]) == false)
	   {std::cerr << "Failed to find neighbours." << std::endl; continue;}
	 if (cellPtr != NULL) cpu_translation3(*cellPtr,nbrPtrs);
      }
      mpiGrid.wait_neighbour_data_update();
      cells = mpiGrid.get_cells_with_remote_neighbour();
      for (size_t c=0; c<cells.size(); ++c) {
	 cellPtr = mpiGrid[cells[c]];
	 if (findNeighbours(nbrPtrs,mpiGrid,cells[c]) == false)
	   {std::cerr << "Failed to find neighbours." << std::endl; continue;}
	 if (cellPtr != NULL) cpu_translation3(*cellPtr,nbrPtrs);
      }

      ++P::tstep;
      
      // Check if the full simulation state should be written to disk
      if (P::tstep % P::saveInterval == 0) {
	 logger << "(MAIN): Saving full state to disk at tstep = " << tstep << std::endl;
      }

      // Check if variables and derived quantities should be written to disk
      if (P::tstep % P::diagnInterval == 0) {
	 logger << "(MAIN): Saving variables to disk at tstep = " << tstep << std::endl;

	 #ifndef NOCUDA
	   // Load data from device before plotting:
	   for (uint cell=0; cell<grid.size(); ++cell) {
	      grid[cell]->devSync(Cell::Blocks,Cell::DevToCpu);
	   }
	 #endif
	 writeSpatialCells(comm,mpiGrid);
	 writeVelocityBlocks(comm,mpiGrid);
      }
      comm.barrier();
   }
   logger << "(MAIN): All timesteps calculated." << std::endl;
   logger << "\t (TIME) reports value " << std::time(NULL) << std::endl;

   // Write final state:
   writeSpatialCells(comm,mpiGrid);
   writeVelocityBlocks(comm,mpiGrid);

   logger << "(MAIN): Exiting." << std::endl;
   return 0;
}




