#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>

#ifndef PARGRID
   #include "main_dccrg.h"
#else
   #include "main_pargrid.h"
#endif

#include "definitions.h"
#include "logger.h"
#include "parameters.h"
#include "grid.h"
#include "silowriter.h"
#include "cell_spatial.h"
#include "writevars.h"

extern void buildSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,creal& zmin,creal& dx,creal& dy,creal& dz,const bool& isRemote);

Grid grid;

using namespace std;

#ifndef PARGRID
void initSpatialCells(dccrg<SpatialCell>& mpiGrid,boost::mpi::communicator& comm) {
#else
void initSpatialCells(const ParGrid<SpatialCell>& mpiGrid) {
#endif
   typedef Parameters P;
   
   // This can be replaced by an iterator.
   #ifndef PARGRID
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   #else
   std::vector<ID::type> cells;
   mpiGrid.getCells(cells);
   #endif
   
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
      
      buildSpatialCell(*(mpiGrid[cells[i]]),xmin,ymin,zmin,dx,dy,dz,false);
   }
   
   #ifdef PARGRID
     mpiGrid.getRemoteCells(cells);
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
	
	buildSpatialCell(*(mpiGrid[cells[i]]),xmin,ymin,zmin,dx,dy,dz,true);
     }
   #endif
}

#ifndef PARGRID
void writeSpatialCells(const boost::mpi::communicator& comm,dccrg<SpatialCell>& mpiGrid) {
#else
void writeSpatialCells(const ParGrid<SpatialCell>& mpiGrid) {
#endif
   // This can be replaced by an iterator.
   #ifndef PARGRID
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   #else
   std::vector<ID::type> cells;
   mpiGrid.getCells(cells);
   #endif

   std::stringstream fname;
   #ifndef PARGRID
   fname << "cells." << comm.rank() << '.';
   #else
   fname << "cells." << mpiGrid.rank() << '.';
   #endif
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

#ifndef PARGRID
   
#else
void writeRemoteCells(const ParGrid<SpatialCell>& mpiGrid) {
   #ifndef PARGRID
   
   #else
      std::vector<ID::type> cells;
      mpiGrid.getRemoteCells(cells);
   #endif
   if (cells.size() == 0) return;
   std::stringstream fname;
   #ifndef PARGRID
   
   #else
      fname << "remcells." << mpiGrid.rank() << '.';
   #endif
   fname.width(5);
   fname.fill('0');
   fname << Parameters::tstep << ".silo";
   
   openOutputFile(fname.str(),"remote_cells");
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
   writeSpatialCells("remotecells","n");
   closeOutputFile();
   freeCells();
}
#endif

#ifndef PARGRID
void writeVelocityBlocks(const boost::mpi::communicator& comm,dccrg<SpatialCell>& mpiGrid) {
#else
void writeVelocityBlocks(const ParGrid<SpatialCell>& mpiGrid) {
#endif
   // This can be replaced by an iterator.
   #ifndef PARGRID
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   #else
   std::vector<ID::type> cells;
   mpiGrid.getCells(cells);
   #endif
   
   std::stringstream fname;
   #ifndef PARGRID
   fname << "blocks." << comm.rank() << '.';
   #else
   fname << "blocks." << mpiGrid.rank() << '.';
   #endif
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

#ifndef PARGRID
void writeVelocityBlocks(const boost::mpi::communicator& comm, dccrg<SpatialCell>& mpiGrid, const uint64_t cell) {
#else
void writeVelocityBlocks(const ParGrid<SpatialCell>& mpiGrid, const ID::type cell) {
#endif
   std::stringstream fname;
   #ifndef PARGRID
   fname << "block_" << cell << "." << comm.rank() << '.';
   #else
   fname << "block_" << cell << "." << mpiGrid.rank() << '.';
   #endif
   fname.width(5);
   fname.fill('0');
   fname << Parameters::tstep << ".silo";

   // Write velocity grid
   openOutputFile(fname.str(),"vel_blocks");

   writeVelocityBlockGridScalar3D("f", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_avgs);
   writeVelocityBlockGridScalar3D("Fx", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_fx);
   writeVelocityBlockGridScalar3D("Fy", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_fy);
   writeVelocityBlockGridScalar3D("Fz", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_fz);
   writeVelocityBlockGridScalar3D("D1x", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_d1x);
   writeVelocityBlockGridScalar3D("D1y", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_d1y);
   writeVelocityBlockGridScalar3D("D1z", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_d1z);
   writeVelocityBlockGridScalar3D("D2x", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_d2x);
   writeVelocityBlockGridScalar3D("D2y", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_d2y);
   writeVelocityBlockGridScalar3D("D2z", "velgrid", mpiGrid[cell]->N_blocks, mpiGrid[cell]->cpu_d2z);
   closeOutputFile();
}

#ifndef PARGRID
void writeAllVelocityBlocks(const boost::mpi::communicator& comm, dccrg<SpatialCell>& mpiGrid) {
#else
void writeAllVelocityBlocks(const ParGrid<SpatialCell>& mpiGrid) {
#endif

   #ifndef PARGRID
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   #else
   std::vector<ID::type> cells;
   mpiGrid.getCells(cells);
   #endif

   for (uint i = 0; i < cells.size(); ++i) {
      #ifndef PARGRID
      writeVelocityBlocks(comm, mpiGrid, cells[i]);
      #else
      writeVelocityBlocks(mpiGrid, cells[i]);
      #endif
   }
}

int main(int argn,char* args[]) {
   bool success = true;
   
   typedef Parameters P;
   Parameters parameters(argn,args);

   #ifndef PARGRID // INITIALIZE USING DCCRG
      boost::mpi::environment env(argn,args);
      boost::mpi::communicator comm;
        {
	   std::stringstream ss;
	   ss << "logfile." << comm.rank() << ".txt";
	   logger.setOutputFile(ss.str());
	}
      logger << "(MAIN): Starting up." << std::endl;
   
      // Create parallel MPI grid and init Zoltan:
      float zoltanVersion;
      if (Zoltan_Initialize(argn,args,&zoltanVersion) != ZOLTAN_OK) {
	 logger << "\t ERROR: Zoltan initialization failed, aborting." << std::endl;
	 success = false;
      } else {
	 logger << "\t Zoltan initialized successfully" << std::endl;
      }
      dccrg<SpatialCell> mpiGrid(comm,"HYPERGRAPH",P::xmin,P::ymin,P::zmin,P::dx_ini,P::dy_ini,P::dz_ini,P::xcells_ini,P::ycells_ini,P::zcells_ini,1,0);
   
   #else           // INITIALIZE USING PARGRID
      ParGrid<SpatialCell> mpiGrid(P::xcells_ini,P::ycells_ini,P::zcells_ini,P::xmin,P::ymin,P::zmin,
				   P::xmax,P::ymax,P::zmax,Hypergraph,argn,args);
        {
	   std::stringstream ss;
	   ss << "logfile." << mpiGrid.rank() << ".txt";
	   logger.setOutputFile(ss.str());
	}
      logger << "(MAIN): Starting up." << std::endl;
   #endif
   
   // If initialization was not successful, abort.
   if (success == false) {
      std::cerr << "An error has occurred, aborting. See logfile for details." << std::endl;
      logger << "Aborting" << std::endl;
      return 1;
   }
   // Do initial load balancing:
   initialLoadBalance(mpiGrid);
   #ifndef PARGRID
      comm.barrier();
   #endif
   
   // Go through every spatial cell on this CPU, and create the initial state:
   #ifndef PARGRID
      initSpatialCells(mpiGrid,comm);
      comm.barrier();
   #else
      initSpatialCells(mpiGrid);
      writeCellDistribution(mpiGrid);
      mpiGrid.barrier();
   #endif

   // Fetch neighbour data:
   #ifndef PARGRID
      P::transmit = Transmit::AVGS;
      mpiGrid.start_remote_neighbour_data_update(); // TEST
      mpiGrid.wait_neighbour_data_update();
      comm.barrier();
   #else
      P::transmit = Transmit::AVGS;
      mpiGrid.startNeighbourExchange(1);
      mpiGrid.waitAll();
      mpiGrid.barrier();
   #endif
   logger << "(MAIN): Total no. reserved velocity blocks in Grid = ";
   #ifdef PARGRID
      logger << grid.getTotalNumberOfBlocks() << std::endl;
   #endif
   
   // Write initial state:
   #ifndef PARGRID
      writeSpatialCells(comm,mpiGrid);
      //writeVelocityBlocks(comm,mpiGrid);
      writeAllVelocityBlocks(comm, mpiGrid);
      comm.barrier();
   #else
      writeSpatialCells(mpiGrid);
      //writeRemoteCells(mpiGrid);
      //writeVelocityBlocks(mpiGrid);
      writeAllVelocityBlocks(mpiGrid);
      mpiGrid.barrier();
   #endif

   // Main simulation loop:
   logger << "(MAIN): Starting main simulation loop." << std::endl;
   time_t before = std::time(NULL);
   logger << "\t (TIME) reports value " << before << std::endl;
   for (uint tstep=0; tstep < P::tsteps; ++tstep) {
      // Check if spatial cell parameters need to be recalculated, i.e. if they 
      // are time-dependent:
      calculateCellParameters(mpiGrid,P::tstep*P::dt);
      
      // Propagate the state of simulation forward in time by dt:
      calculateAcceleration(mpiGrid);
      calculateSpatialDerivatives(mpiGrid);
      calculateSpatialFluxes(mpiGrid);
      calculateSpatialPropagation(mpiGrid);
      ++P::tstep;
      
      // Check if the full simulation state should be written to disk
      if (P::tstep % P::saveInterval == 0) {
	 logger << "(MAIN): Saving full state to disk at tstep = " << tstep << std::endl;
	 #ifdef PARGRID
	    logger << "\t # sends to other MPI processes      = " << mpiGrid.getNumberOfSends() << std::endl;
	    logger << "\t # receives from other MPI processes = " << mpiGrid.getNumberOfReceives() << std::endl;
	 #endif
      }

      // Check if variables and derived quantities should be written to disk
      if (P::tstep % P::diagnInterval == 0) {
	 logger << "(MAIN): Saving variables to disk at tstep = " << tstep << std::endl;
	 #ifndef PARGRID
	    writeSpatialCells(comm,mpiGrid);
	    //writeVelocityBlocks(comm,mpiGrid);
	 #else
	    writeSpatialCells(mpiGrid);
	    //writeRemoteCells(mpiGrid);
	    //writeVelocityBlocks(mpiGrid);
	 #endif
      }
      #ifndef PARGRID
         comm.barrier();
      #else
         mpiGrid.barrier();
      #endif
   }
   logger << "(MAIN): All timesteps calculated." << std::endl;
   time_t after = std::time(NULL);
   logger << "\t (TIME) reports value " << after << std::endl;
   logger << "\t (TIME) total run time " << after - before << " s" << std::endl;

   // Write final state:
   #ifndef PARGRID
      writeSpatialCells(comm,mpiGrid);
      //writeVelocityBlocks(comm,mpiGrid);
   #else
      writeSpatialCells(mpiGrid);
      //writeRemoteCells(mpiGrid);
      //writeVelocityBlocks(mpiGrid);
   #endif

   logger << "(MAIN): Exiting." << std::endl;
   return 0;
}




