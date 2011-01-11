#include "boost/mpi.hpp"
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

#include "datareducer.h"
#include "datareductionoperator.h"
#include "vlswriter.h"

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
      xmin = mpiGrid.get_cell_x_min(cells[i]);
      ymin = mpiGrid.get_cell_y_min(cells[i]);
      zmin = mpiGrid.get_cell_z_min(cells[i]);
      
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
   fname.width(7);
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
   double x = mpiGrid.get_cell_x(cell);
   double y = mpiGrid.get_cell_y(cell);
   double z = mpiGrid.get_cell_z(cell);
   fname << "block_" << x / 6.3712e6 << "_" << y / 6.3712e6 << "_" << z / 6.3712e6 << "_";
   fname.width(7);
   fname.fill('0');
   fname << Parameters::tstep << ".silo";

   // Write velocity grid
   openOutputFile(fname.str(),"vel_blocks");
   reserveVelocityBlocks(mpiGrid[cell]->N_blocks);
   for (uint i = 0; i < mpiGrid[cell]->N_blocks; ++i) {
      addVelocityGridBlock3D(mpiGrid[cell]->cpu_blockParams + i * SIZE_BLOCKPARAMS);
   }
   writeReservedBlocks("velgrid");
   freeBlocks();

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


#ifndef PARGRID
void writeSomeVelocityGrids(const boost::mpi::communicator& comm, dccrg<SpatialCell>& mpiGrid, const std::vector<Real> x, const std::vector<Real> y, const std::vector<Real> z) {
#else
void writeSomeVelocityGrids(const ParGrid<SpatialCell>& mpiGrid, const std::vector<Real> x, const std::vector<Real> y, const std::vector<Real> z) {
#endif

	#ifndef PARGRID
	std::vector<uint64_t> cells = mpiGrid.get_cells();
	#else
	std::vector<ID::type> cells;
	mpiGrid.getCells(cells);
	#endif

	if (x.size() != y.size() || x.size() != z.size()) {
		cerr << "writeSomeVelocityGrids: x, y and z sizes must be equal" << endl;
		exit(EXIT_FAILURE);
	}

   for (int i = 0; i < x.size(); i++) {
      for (unsigned int j = 0; j < cells.size(); j++) {
         #ifndef PARGRID
	 double cell_x = mpiGrid.get_cell_x(cells[j]);
	 double cell_y = mpiGrid.get_cell_y(cells[j]);
	 double cell_z = mpiGrid.get_cell_z(cells[j]);
	 double cell_dx = mpiGrid.get_cell_x_size(cells[j]);
	 double cell_dy = mpiGrid.get_cell_y_size(cells[j]);
	 double cell_dz = mpiGrid.get_cell_z_size(cells[j]);
	 
	 if (fabs(x[i] - cell_x) <= cell_dx / 2
	     && fabs(y[i] - cell_y) <= cell_dy / 2
	     && fabs(z[i] - cell_z) <= cell_dz / 2) {
	    writeVelocityBlocks(comm, mpiGrid, cells[j]);
	 }
         #else
         #warning writeSomeVelocityGrids not supported with PARGRID
	 //writeVelocityBlocks(mpiGrid, cells[i]);
         #endif
      }
   }
}

#ifdef PARGRID
bool writeSpatialCellData(const ParGrid<SpatialCell>& mpiGrid,VlsWriter& vlsWriter,DataReducer& dataReducer) {
#else
bool writeSpatialCellData(const dccrg<SpatialCell>& mpiGrid,VlsWriter& vlsWriter,DataReducer& dataReducer) {
#endif
   bool success = true;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

   stringstream fname;
   fname << "celldata.";
   fname.width(7);
   fname.fill('0');
   fname << Parameters::tstep << ".vlsv";
   
   // Open a file for writing using MPI I/O:
   if (vlsWriter.openWrite(MPI_COMM_WORLD,fname.str()) == false) {
      logger << "Error opening output file on process " << mpiGrid.rank() << endl;
      success = false;
   }
   
   // Master process (rank=0) within the given communicator writes the file header:
   if (vlsWriter.writeHeader(MPI_COMM_WORLD,0) == false) {
      logger << "Error writing header on process " << mpiGrid.rank() << endl;
      success = false;
   }
   
   // Master process writes description of static-size variables. Here "static-size" 
   // means that the size of variable data is the same for all spatial cells.
   if (vlsWriter.writeStaticVariableDesc(MPI_COMM_WORLD,0,&dataReducer) == false) {
      logger << "Error writing variable description on process " << mpiGrid.rank() << endl;
      success = false;
   }
   
   // Get the global IDs of all local cells:
   #ifdef PARGRID
      mpiGrid.getCells(Main::cells);
      ID::type cellGID;
   #else
      Main::cells = mpiGrid.get_cells();
      uint64_t cellGID;
   #endif
   
   // Write local cell data to file:
   SpatialCell* cellptr;
   for (size_t i=0; i<Main::cells.size(); ++i) {
      cellGID = Main::cells[i];
      cellptr = mpiGrid[cellGID];
      if (cellptr != NULL) {
	 if (vlsWriter.writeSpatCellCoordEntry(cellGID,*cellptr,&dataReducer) == false) {
	    logger << "Error writing spatial cell with global ID " << cellGID << " on process " << myrank << endl;
	    success = false;
	 }
      } else {
	 cerr << "Proc #" << myrank << " received NULL pointer, i = " << i << " global ID = " << Main::cells[i] << endl;
	 success = false;
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   
   // Close output file and exit:
   if (vlsWriter.close() == false) {
      logger << "Error closing file on process " << mpiGrid.rank() << endl;
      success = false;
   }
   return success;
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

   #ifdef PARGRID
   logger << "(MAIN): Total no. reserved velocity blocks in Grid = ";
   logger << grid.getTotalNumberOfBlocks() << std::endl;
   #endif
   
   // Initialize data reduction operators. This should be done elsewhere in order to initialize 
   // user-defined operators:
   DataReducer reducer;
   reducer.addOperator(new DRO::VariableB);
   reducer.addOperator(new DRO::VariableE);
   reducer.addOperator(new DRO::VariableRho);
   reducer.addOperator(new DRO::VariableRhoV);
   reducer.addOperator(new DRO::MPIrank);
   VlsWriter vlsWriter;

   // Write initial state:
   if (P::save_spatial_grid) {
      if (writeSpatialCellData(mpiGrid,vlsWriter,reducer) == false) {
	 logger << "(MAIN): ERROR occurred while writing data to file!" << endl;
      }
   }
   #ifndef PARGRID
      if (P::save_velocity_grid) {
	 writeAllVelocityBlocks(comm, mpiGrid);
      }
      writeSomeVelocityGrids(comm, mpiGrid, P::save_spatial_cells_x, P::save_spatial_cells_y, P::save_spatial_cells_z);
      comm.barrier();
   #else
      if (P::save_velocity_grid) {
	 writeAllVelocityBlocks(mpiGrid);
      }
      mpiGrid.barrier();
   #endif
   
   // Main simulation loop:
   logger << "(MAIN): Starting main simulation loop." << std::endl;
   time_t before = std::time(NULL);
   for (uint tstep=0; tstep < P::tsteps; ++tstep) {
      // Recalculate (maybe) spatial cell parameters
      calculateSimParameters(mpiGrid, P::t, P::dt);

      // use globally minimum timestep
      #ifndef PARGRID
      P::dt = all_reduce(comm, P::dt, boost::mpi::minimum<Real>());
      #else
      #warning No communicator for all_reduce when using PARGRID
      #endif

      // Propagate the state of simulation forward in time by dt:
      calculateAcceleration(mpiGrid);
      calculateSpatialDerivatives(mpiGrid);
      calculateSpatialFluxes(mpiGrid);
      calculateSpatialPropagation(mpiGrid);
      ++P::tstep;
      P::t += P::dt;
      
      // Check if the full simulation state should be written to disk
      if (P::tstep % P::saveRestartInterval == 0) {
         // TODO: implement full state saving
         logger << "(MAIN): Saving full state to disk at tstep = " << tstep << ", time = " << P::t << std::endl;
         #ifndef PARGRID
            logger << "\t # sends to other MPI processes      = " << mpiGrid.get_number_of_update_send_cells() << std::endl;
            logger << "\t # receives from other MPI processes = " << mpiGrid.get_number_of_update_receive_cells() << std::endl;
            logger << "\t # total = " << mpiGrid.get_number_of_update_send_cells() + mpiGrid.get_number_of_update_receive_cells() << std::endl;
         #else
            logger << "\t # sends to other MPI processes      = " << mpiGrid.getNumberOfSends() << std::endl;
            logger << "\t # receives from other MPI processes = " << mpiGrid.getNumberOfReceives() << std::endl;
            logger << "\t # total = " << mpiGrid.getNumberOfSends() + mpiGrid.getNumberOfReceives() << std::endl;
         #endif
      }

      // Check if variables and derived quantities should be written to disk
      if (P::tstep % P::diagnInterval == 0) {
         logger << "(MAIN): Saving variables to disk at tstep = " << tstep << ", time = " << P::t << std::endl;

	 if (P::save_spatial_grid) {
	    if (writeSpatialCellData(mpiGrid,vlsWriter,reducer) == false) {
	       logger << "(MAIN): ERROR occurred while writing data to file!" << endl;
	    }
	 }
	 
         #ifndef PARGRID
         if (P::save_velocity_grid) {
            writeAllVelocityBlocks(comm, mpiGrid);
         }
         writeSomeVelocityGrids(comm, mpiGrid, P::save_spatial_cells_x, P::save_spatial_cells_y, P::save_spatial_cells_z);

         #else
         if (P::save_velocity_grid) {
            writeAllVelocityBlocks(mpiGrid);
         }

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
   logger << "\t (TIME) total run time " << after - before << " s, total simulated time " << P::t << " s" << std::endl;
   logger << "\t (TIME) seconds per timestep " << double(after - before) / P::tsteps << ", seconds per simulated second " << double(after - before) / P::t << std::endl;

   // Write final state:
   if (P::save_spatial_grid) {
      if (writeSpatialCellData(mpiGrid,vlsWriter,reducer) == false) {
	 logger << "(MAIN): ERROR occurred while writing data to file!" << endl;
      }
   }


   #ifndef PARGRID
   if (P::save_velocity_grid) {
      writeAllVelocityBlocks(comm, mpiGrid);
   }
   writeSomeVelocityGrids(comm, mpiGrid, P::save_spatial_cells_x, P::save_spatial_cells_y, P::save_spatial_cells_z);

   #else

   if (P::save_velocity_grid) {
      writeAllVelocityBlocks(mpiGrid);
   }

   #endif

   logger << "(MAIN): Exiting." << std::endl;
   return 0;
}

