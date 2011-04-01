

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
#include "mpiconversion.h"
#include "mpilogger.h"
#include "parameters.h"
#include "grid.h"
#include "silowriter.h"
#include "cell_spatial.h"
#include "writevars.h"
#include "gridbuilder.h"
#include "datareducer.h"
#include "datareductionoperator.h"
#include "vlswriter.h"

#include "vlsvwriter2.h" // TEST

#ifdef CRAYPAT
//include craypat api headers if compiled with craypat on Cray XT/XE
#include "pat_api.h"
#endif 
Grid grid;

using namespace std;




#ifndef PARGRID
void initSpatialCells(const dccrg<SpatialCell>& mpiGrid,boost::mpi::communicator& comm) {
#else
void initSpatialCells(const ParGrid<SpatialCell>& mpiGrid) {
#endif
    
    typedef Parameters P;


   // This can be replaced by an iterator.
   #ifndef PARGRID
     Main::cells = mpiGrid.get_cells();
   #else
     mpiGrid.getCells(Main::cells);
   #endif
   
   // Go through every cell on this node and initialize the pointers to 
   // cpu memory, physical parameters and volume averages for each phase space 
   // point in the velocity grid. Velocity block neighbour list is also 
   // constructed here:
   Real xmin,ymin,zmin,dx,dy,dz;
   for (uint i=0; i<Main::cells.size(); ++i) {
      dx = mpiGrid.get_cell_x_size(Main::cells[i]);
      dy = mpiGrid.get_cell_y_size(Main::cells[i]);
      dz = mpiGrid.get_cell_z_size(Main::cells[i]);
      xmin = mpiGrid.get_cell_x_min(Main::cells[i]);
      ymin = mpiGrid.get_cell_y_min(Main::cells[i]);
      zmin = mpiGrid.get_cell_z_min(Main::cells[i]);
      buildSpatialCell(*(mpiGrid[Main::cells[i]]),xmin,ymin,zmin,dx,dy,dz,false);
   }
   #ifdef PARGRID
     // For ParGrid memory for remote cells needs to be allocated here:
     mpiGrid.getRemoteCells(Main::cells);
     for (uint i=0; i<Main::cells.size(); ++i) {
	dx = mpiGrid.get_cell_x_size(Main::cells[i]);
	dy = mpiGrid.get_cell_y_size(Main::cells[i]);
	dz = mpiGrid.get_cell_z_size(Main::cells[i]);
	xmin = mpiGrid.get_cell_x_min(Main::cells[i]);
	ymin = mpiGrid.get_cell_y_min(Main::cells[i]);
	zmin = mpiGrid.get_cell_z_min(Main::cells[i]);
	buildSpatialCell(*(mpiGrid[Main::cells[i]]),xmin,ymin,zmin,dx,dy,dz,true);
     }
   #endif
}

#ifndef PARGRID
void writeVelocityBlocks(const boost::mpi::communicator& comm, const dccrg<SpatialCell>& mpiGrid) {
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
void writeVelocityBlocks(const dccrg<SpatialCell>& mpiGrid, const uint64_t cell) {
#else
void writeVelocityBlocks(const ParGrid<SpatialCell>& mpiGrid, const ID::type cell) {
#endif

   if (mpiGrid[cell] == NULL) {
      return;
   }

   std::stringstream fname;
   fname.precision(1);
   double x = mpiGrid.get_cell_x(cell);
   double y = mpiGrid.get_cell_y(cell);
   double z = mpiGrid.get_cell_z(cell);
   fname << "block_" << std::fixed << x / 6.3712e6 << "_" << y / 6.3712e6 << "_" << z / 6.3712e6 << "_";
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
void writeAllVelocityBlocks(const dccrg<SpatialCell>& mpiGrid) {
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
      writeVelocityBlocks(mpiGrid, cells[i]);
   }
}


#ifndef PARGRID
void writeSomeVelocityGrids(const dccrg<SpatialCell>& mpiGrid, const std::vector<Real> x, const std::vector<Real> y, const std::vector<Real> z) {
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
			double cell_x = mpiGrid.get_cell_x(cells[j]);
			double cell_y = mpiGrid.get_cell_y(cells[j]);
			double cell_z = mpiGrid.get_cell_z(cells[j]);
			double cell_dx = mpiGrid.get_cell_x_size(cells[j]);
			double cell_dy = mpiGrid.get_cell_y_size(cells[j]);
			double cell_dz = mpiGrid.get_cell_z_size(cells[j]);

			if (fabs(x[i] - cell_x) <= cell_dx / 2
			&& fabs(y[i] - cell_y) <= cell_dy / 2
			&& fabs(z[i] - cell_z) <= cell_dz / 2) {
				writeVelocityBlocks(mpiGrid, cells[j]);
			}
		}
	}
}

#ifdef PARGRID
bool writeGrid(const ParGrid<SpatialCell>& mpiGrid,DataReducer& dataReducer) {
#else
bool writeGrid(const dccrg<SpatialCell>& mpiGrid,DataReducer& dataReducer) {
#endif
   double allStart = MPI_Wtime();
   bool success = true;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

   // Create a name for the output file and open it with VLSVWriter:
   stringstream fname;
   fname << "grid.";
   fname.width(7);
   fname.fill('0');
   fname << Parameters::tstep << ".vlsv";
   
   VLSVWriter vlsvWriter;
   vlsvWriter.open(fname.str(),MPI_COMM_WORLD,0);
   
   // Get all local cell IDs and write to file:
   map<string,string> attribs;
   #ifdef PARGRID
      mpiGrid.getCells(Main::cells);
   #else 
      Main::cells = mpiGrid.get_cells();
   #endif

   if (vlsvWriter.writeArray("MESH","SpatialGrid",attribs,Main::cells.size(),1,&(Main::cells[0])) == false) {
      cerr << "Proc #" << myrank << " failed to write cell IDs!" << endl;
   }

   // Create a buffer for spatial cell coordinates. Copy all coordinates to 
   // buffer and write:
   Real* buffer = new Real[6*Main::cells.size()];
   for (size_t i=0; i<Main::cells.size(); ++i) {
      SpatialCell* SC = mpiGrid[Main::cells[i]];
      for (int j=0; j<6; ++j) {
	 buffer[6*i+j] = SC->cpu_cellParams[j];
      }
   }
   if (vlsvWriter.writeArray("COORDS","SpatialGrid",attribs,Main::cells.size(),6,buffer) == false) {
      cerr << "Proc #" << myrank << " failed to write cell coords!" << endl;
   }
   delete buffer;


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
      uint64_t arraySize = Main::cells.size()*vectorSize*dataSize;
      
      // Request DataReductionOperator to calculate the reduced data for all local cells:
      varBuffer = new char[arraySize];
      for (uint64_t cell=0; cell<Main::cells.size(); ++cell) {
	 if (dataReducer.reduceData(mpiGrid[Main::cells[cell]],i,varBuffer + cell*vectorSize*dataSize) == false) success = false;
      }
      
      // Write reduced data to file:
      if (vlsvWriter.writeArray("VARIABLE",variableName,attribs,Main::cells.size(),vectorSize,dataType,dataSize,varBuffer) == false) success = false;      
      delete varBuffer;
      varBuffer = NULL;
   }

   // Write velocity blocks and related data. Which cells write velocity grids 
   // should be requested from a function, but for now we just write velocity grids for all cells:
   attribs.clear();
   
   // First write global IDs of those cells which write velocity blocks (here: all cells):
   if (vlsvWriter.writeArray("CELLSWITHBLOCKS","SpatialGrid",attribs,Main::cells.size(),1,&(Main::cells[0])) == false) success = false;
   
   // Write the number of velocity blocks in each spatial cell. Again a temporary buffer is used:
   uint* N_blocks = new uint[Main::cells.size()];
   uint64_t totalBlocks = 0;
   for (size_t cell=0; cell<Main::cells.size(); ++cell) {
      N_blocks[cell] = mpiGrid[Main::cells[cell]]->N_blocks;
      totalBlocks += mpiGrid[Main::cells[cell]]->N_blocks;
   }
   if (vlsvWriter.writeArray("NBLOCKS","SpatialGrid",attribs,Main::cells.size(),1,N_blocks) == false) success = false;

   double start = MPI_Wtime();
   
   // Write velocity block coordinates. At this point the data may get too large to be buffered, 
   // so a "multi-write" mode is used - coordinates are written one velocity grid at a time. Besides, 
   // the velocity block coordinates are already stored in a suitable format for writing in SpatialCells:
   if (vlsvWriter.startMultiwrite("float",totalBlocks,6,sizeof(Real)) == false) success = false;
   if (success == true) {
      uint64_t counter = 0;
      SpatialCell* SC;
      for (size_t cell=0; cell<Main::cells.size(); ++cell) {
	 SC = mpiGrid[Main::cells[cell]];
	 if (vlsvWriter.multiwriteArray(counter,N_blocks[cell],SC->cpu_blockParams) == false) success = false;
	 counter += N_blocks[cell];
      }
   }
   if (success == true) if (vlsvWriter.endMultiwrite("BLOCKCOORDINATES","SpatialGrid",attribs) == false) success = false;
   
   // Write values of distribution function:
   if (vlsvWriter.startMultiwrite("float",totalBlocks,64,sizeof(Real)) == false) success = false;
   if (success == true) {
      uint64_t counter = 0;
      SpatialCell* SC;
      for (size_t cell=0; cell<Main::cells.size(); ++cell) {
	 SC = mpiGrid[Main::cells[cell]];
	 if (vlsvWriter.multiwriteArray(counter,N_blocks[cell],SC->cpu_avgs) == false) success = false;
	 counter += N_blocks[cell];
      }
   }
   attribs["mesh"] = "SpatialGrid";
   if (success == true) if (vlsvWriter.endMultiwrite("BLOCKVARIABLE","avgs",attribs) == false) success = false;
   
   double end = MPI_Wtime();
    
   delete N_blocks;

   vlsvWriter.close();

   double allEnd = MPI_Wtime();
   
   double bytesWritten = Main::cells.size()*1000*4*(64+6);
   double secs = end-start;
   mpilogger << "Wrote " << bytesWritten/1.0e6 << " MB of data in " << secs << " seconds, datarate is " << bytesWritten/secs/1.0e9 << " GB/s" << endl << write;
   
   double allSecs = allEnd-allStart;
   if (myrank == 0) mpilogger << "All data written in " << allSecs << " seconds" << endl << write;
   
   return success;
}
   
#ifdef PARGRID
bool writeSpatialCellData(const ParGrid<SpatialCell>& mpiGrid,VlsWriter& vlsWriter,DataReducer& dataReducer) {
#else
bool writeSpatialCellData(const dccrg<SpatialCell>& mpiGrid,VlsWriter& vlsWriter,DataReducer& dataReducer) {
#endif
   bool success = true;
   bool writeSpatNbrLists = false;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   
   stringstream fname;
   fname << "celldata.";
   fname.width(7);
   fname.fill('0');
   fname << Parameters::tstep << ".vlsv";
   
   // Open a file for writing using MPI I/O:
   if (vlsWriter.open(MPI_COMM_WORLD,fname.str()) == false) {
      mpilogger << "Error opening output file on process " << myrank << endl << write;
      success = false;
   }

   // Set some parameters:
   #ifdef PARGRID
      vlsWriter.setBytesPerCellGID(sizeof(ID::type));
   #else
      vlsWriter.setBytesPerCellGID(sizeof(uint64_t));
   #endif
   vlsWriter.setWriteSpatNbrsLists(writeSpatNbrLists);
   
   // Master process (rank=0) within the given communicator writes the file header:
   if (vlsWriter.writeHeader(MPI_COMM_WORLD,0) == false) {
      mpilogger << "Error writing header on process " << myrank << endl << write;
      success = false;
   }

   // Master process writes description of static-size variables. Here "static-size" 
   // means that the size of variable data is the same for all spatial cells.
   if (vlsWriter.writeStaticVariableDesc(MPI_COMM_WORLD,0,&dataReducer) == false) {
      mpilogger << "Error writing variable description on process " << myrank << endl << write;
      success = false;
   }

   // Get the global IDs of all local cells:
   #ifdef PARGRID
      mpiGrid.getCells(Main::cells);
      ID::type cellGID;
      vector<ID::type> nbrs;
   #else
      Main::cells = mpiGrid.get_cells();
      uint64_t cellGID;
      vector<uint64_t> nbrs;
   #endif
   
   // Write local cell data to file using a buffer:
   if (vlsWriter.reserveSpatCellCoordBuffer(Main::cells.size(),&dataReducer) == false) {
      mpilogger << "ERROR: VlsWriter failed to reserve cell buffer!" << endl;
      success = false;
   }
   
   // If errors have come up, exit:
   if (success == false) {
      vlsWriter.close();
      return success;
   }
   
   // Write spatial cell data:
   SpatialCell* cellptr;
   unsigned char refLevel;
   for (size_t i=0; i<Main::cells.size(); ++i) {
      cellGID = Main::cells[i];
      cellptr = mpiGrid[cellGID];
      
      if (writeSpatNbrLists == true) {
         #ifdef PARGRID
            mpiGrid.getExistingNeighbours(nbrs,cellGID);
            refLevel = mpiGrid.getRefinementLevel(cellGID);
            if (refLevel == numeric_limits<unsigned char>::max())
	       mpilogger << "Received erroneous refinement level for cell with global ID " << cellGID << endl << write;
         #else
	    nbrs.clear();
            refLevel = 0;
         #endif
      }

      if (cellptr != NULL) {
	 if (vlsWriter.writeSpatCellCoordEntryBuffered(cellGID,*cellptr,&dataReducer,nbrs,refLevel) == false) {
	    mpilogger << "Error writing spatial cell with global ID " << cellGID << " on process " << myrank << endl << write;
	    success = false;
	 }
      } else {
	 mpilogger << "ERROR: Received NULL pointer, i = " << i << " global ID = " << Main::cells[i] << endl << write;
	 success = false;
      }
   }
   if (vlsWriter.flushBuffer() == false) {
      mpilogger << "ERROR: VlsWriter failed to flush buffer!" << endl << write;
      success = false;
   }
   
   // Wait until all processes have completed writing before inserting the end marker:
   MPI_Barrier(MPI_COMM_WORLD);
   if (vlsWriter.writeSpatCellCoordEntryEndMarker(MPI_COMM_WORLD,0) == false) {
      mpilogger << "ERROR: VlsWriter failed to write end marker!" << endl << write;
      success = false;
   }
   
   // Close output file and exit:
   if (vlsWriter.close() == false) {
      mpilogger << "Error closing file on process " << myrank << endl << write;
      success = false;
   }
   return success;
}

#ifdef PARGRID
void log_send_receive_info(const ParGrid<SpatialCell>& mpiGrid) {
#else
void log_send_receive_info(const dccrg<SpatialCell>& mpiGrid) {
#endif
   mpilogger << "Number of sends / receives:" << endl;
   #ifndef PARGRID
   mpilogger << "\tto other MPI processes   = " << mpiGrid.get_number_of_update_send_cells() << endl;
   mpilogger << "\tfrom other MPI processes = " << mpiGrid.get_number_of_update_receive_cells() << endl;
   mpilogger << "\ttotal = " << mpiGrid.get_number_of_update_send_cells() + mpiGrid.get_number_of_update_receive_cells() << endl;
   #else
   mpilogger << "\tto other MPI processes   = " << mpiGrid.getNumberOfSends() << endl;
   mpilogger << "\tfrom other MPI processes = " << mpiGrid.getNumberOfReceives() << endl;
   mpilogger << "\ttotal = " << mpiGrid.getNumberOfSends() + mpiGrid.getNumberOfReceives() << endl;
   #endif
   mpilogger << write;
}

int main(int argn,char* args[]) {
   bool success = true;
   double totTime;
   double initTime;

   
   // Init MPI:
   #ifndef PARGRID
      boost::mpi::environment env(argn,args);
      boost::mpi::communicator comm;
   #else
      if (MPI_Init(&argn,&args) != MPI_SUCCESS) {
	 cerr << "(MAIN): MPI init failed!" << endl;
	 exit(1);
      }
   #endif
   const int MASTER_RANK = 0;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

      
   totTime=MPI_Wtime();
   initTime=MPI_Wtime();    

   // Init parameter file reader:
   typedef Parameters P;
   typedef Readparameters RP;
   Readparameters readparameters(argn,args,MPI_COMM_WORLD);
   readparameters.parse();
   if (readparameters.isInitialized() == false) {
       success = false;
       cerr << "(MAIN) Readparameters failed to init, aborting!" << endl;
       return 1;
   }
   
   // Read in some parameters
   // Parameters related to solar wind simulations:
   // DEPRECATED: These will be moved to somewhere else in the future

   RP::add("solar_wind_file","Read solar wind data from the file arg","");        
   
   
   // Parameters related to saving data:
   // WARNING: Some of these parameters may become deprecated in the future.

   RP::add("save_interval", "Save the simulation every arg time steps",1);
   RP::add("restart_interval","Save the complete simulation every arg time steps",numeric_limits<uint>::max());
   RP::add("save_spatial_grid", "Save spatial cell averages for the whole simulation",true);
   RP::add("save_velocity_grid","Save velocity grid from every spatial cell in the simulation",false);
   RP::addComposing("save_spatial_cells_at_x,X","Save the velocity grid in spatial cells at these coordinates (x components, also give as many y and z components, values from command line, configuration files and environment variables are added together [short version only works on command line])");
   RP::addComposing("save_spatial_cells_at_y,Y","Save the velocity grid in spatial cells at these (y components, also give as many x and z components, values from command line, configuration files and environment variables are added together [short version only works on command line])");
   
   RP::addComposing("save_spatial_cells_at_z,Z","Save the velocity grid in spatial cells at these coordinates (z components, also give as many x and y components, values from command line, configuration files and environment variables are added together [short version only works on command line])");
   
   RP::parse();
   RP::get("solar_wind_file",P::solar_wind_file);
   RP::get("save_interval", P::diagnInterval);
   RP::get("restart_interval", P::saveRestartInterval);
   RP::get("save_spatial_grid", P::save_spatial_grid);
   RP::get("save_velocity_grid", P::save_velocity_grid);
   RP::get("save_spatial_cells_at_x,X", P::save_spatial_cells_x);
   RP::get("save_spatial_cells_at_y,Y", P::save_spatial_cells_y);
   RP::get("save_spatial_cells_at_z,Z", P::save_spatial_cells_z);

   
   // Sanity checks (DEPRECATED):
   if (P::save_spatial_cells_x.size() != P::save_spatial_cells_y.size()
       || P::save_spatial_cells_x.size() != P::save_spatial_cells_z.size()) {
       cerr << "Must have equal number of values for x, y and z components, but given: x " << P::save_spatial_cells_x.size();
       cerr << ", y " << P::save_spatial_cells_y.size() << ", z " << P::save_spatial_cells_z.size() << endl;
       MPI_Finalize();
       exit(1);
   }
   
   



#ifdef CRAYPAT
/*initialize variables for reducing sampling & tracing to a slice in iteration & process
  space.
*/
   char *envvar;
   int firstiter=-1;
   int lastiter=numeric_limits<int>::max();
   int firstrank=-1;
   int lastrank=numeric_limits<int>::max();       
#warning Including CrayPAT API 
   if(myrank==0){
       envvar=getenv("CRAYPAT_FIRSTITER");
       if(envvar!=NULL){
           //FIXME, happily we do no error checking here...
           firstiter=atoi(envvar);
       }
       envvar=getenv("CRAYPAT_LASTITER");
       if(envvar!=NULL){
       //FIXME, happily we do no error checking here...
           lastiter=atoi(envvar);
       }
       envvar=getenv("CRAYPAT_FIRSTRANK");
       if(envvar!=NULL){
           //FIXME, happily we do no error checking here...
           firstrank=atoi(envvar);
       }
       envvar=getenv("CRAYPAT_LASTRANK");
       if(envvar!=NULL){
           //FIXME, happily we do no error checking here...
           lastrank=atoi(envvar);
       }
   }
   
   MPI_Bcast(&firstiter,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&lastiter,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&firstrank,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&lastrank,1,MPI_INT,0,MPI_COMM_WORLD);

   //turn off craypat sampling & tracing based on the defined slices.
   if(firstiter>0){
       PAT_state(PAT_STATE_OFF);
   }
   else if(myrank<firstrank || myrank> lastrank){
       PAT_state(PAT_STATE_OFF);
   }
#endif

   
   // Init parallel logger:
   if (mpilogger.open(MPI_COMM_WORLD,MASTER_RANK,"logfile.txt") == false) {
      cerr << "(MAIN) ERROR: MPILogger failed to open output file!" << endl;
      exit(1);
   }
   
   #ifndef PARGRID // INITIALIZE USING DCCRG
      // Create parallel MPI grid and init Zoltan:
      float zoltanVersion;
      if (Zoltan_Initialize(argn,args,&zoltanVersion) != ZOLTAN_OK) {
	 mpilogger << "\t ERROR: Zoltan initialization failed, aborting." << std::endl << write;
	 success = false;
      } else {
	 mpilogger << "\t Zoltan initialized successfully" << std::endl << write;
      }
      if (buildGrid(MPI_COMM_WORLD,MASTER_RANK) == false) {
	 mpilogger << "(MAIN) Grid builder failed!" << endl << write;
	 success = false;
      } else {
	 mpilogger << "(MAIN) Grid built successfully" << endl << write;
      }

      
      dccrg<SpatialCell> mpiGrid(comm,"HIER",P::xmin,P::ymin,P::zmin,P::dx_ini,P::dy_ini,P::dz_ini,P::xcells_ini,P::ycells_ini,P::zcells_ini,0,0);
      //read in partitioning levels from input
      RP::addComposing("dccrg.partition_procs","Procs per load balance group");
      RP::addComposing("dccrg.partition_lb_method","Load balance method");
      RP::addComposing("dccrg.partition_imbalance_tol","Imbalance tolerance");

      RP::parse();
      vector<int> partitionProcs;
      vector<string> partitionLbMethod;
      vector<string> partitionImbalanceTol;
      RP::get("dccrg.partition_procs",partitionProcs);
      RP::get("dccrg.partition_lb_method",partitionLbMethod);
      RP::get("dccrg.partition_imbalance_tol",partitionImbalanceTol);

      //check that all options set for all levels
      if ( partitionProcs.size()!=partitionLbMethod.size() ||
           partitionProcs.size()!=partitionImbalanceTol.size()){
          if(myrank==0){
              cerr << "DCCRG partition levels not defined completely! Needed options per level are:" <<endl;
              cerr << " partition_procs "<<endl << " partition_lb_method " <<endl << " partition_imbalance_tol " <<endl;
          }
          MPI_Finalize();
          exit(1);
      }
       
      // set default values if nothing has been defined
      if ( partitionProcs.size()==0){
          partitionProcs.push_back(1);
          partitionLbMethod.push_back("HYPERGRAPH");
          partitionImbalanceTol.push_back("1.1");
      }

      //set options 
      for(int i=0;i<partitionProcs.size();i++){
          // set hierarchial partitioning parameters for first level
          if(myrank==0){
              mpilogger << "(MAIN) Partition parameters level "<<i <<" procs: " << partitionProcs[i] << " LB_METHOD: " <<
                  partitionLbMethod[i]<< " IMBALANCE_TOL: "<< partitionImbalanceTol[i] << endl << write;
          }
          
          mpiGrid.add_partitioning_level(partitionProcs[i]);
          mpiGrid.add_partitioning_option(i, "LB_METHOD", partitionLbMethod[i]);
          mpiGrid.add_partitioning_option(i, "IMBALANCE_TOL", partitionImbalanceTol[i]);
      }
      

  
      
   
   #else           // INITIALIZE USING PARGRID
      ParGrid<SpatialCell> mpiGrid(Hypergraph,argn,args);

      // Add all cells to mpiGrid:
      if (buildGrid(mpiGrid,MPI_COMM_WORLD,MASTER_RANK) == false) {
	 mpilogger << "(MAIN) Grid builder failed!" << endl << write;
	 success = false;
      } else {
	 mpilogger << "(MAIN) Grid built successfully" << endl << write;
      }
      // Load balance is most likely far from optimal. Do an 
      // initial load balance before reading cell data:
      //mpiGrid.initialize();
   #endif

   // If initialization was not successful, abort.
   if (success == false) {
      if (myrank == MASTER_RANK) {
	 std::cerr << "An error has occurred, aborting. See logfile for details." << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      mpilogger.close();
      return 1;
   }
   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Starting up." << endl << write;
   
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
      //initSpatialCells(mpiGrid);
      //mpiGrid.barrier();
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

   log_send_receive_info(mpiGrid);

   #ifdef PARGRID
   mpilogger << "(MAIN): Total no. reserved velocity blocks in Grid = ";
   mpilogger << grid.getTotalNumberOfBlocks() << std::endl << write;
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

   // ***********************************
   // ***** INITIALIZATION COMPLETE *****
   // ***********************************
   
   // Free up memory:
   readparameters.finalize();
   initTime=MPI_Wtime()-initTime;
   
   double initIoTime=MPI_Wtime();
   double loopIoTime=0;
   // Write initial state:
   if (P::save_spatial_grid) {
      #ifdef PARGRID
         mpiGrid.getCells(Main::cells);
      #else
         Main::cells = mpiGrid.get_cells();
      #endif
      for (size_t i=0; i<Main::cells.size(); ++i) {
	 cpu_calcVelocityMoments(*(mpiGrid[Main::cells[i]]));
      }
      if (myrank == MASTER_RANK) {
	 mpilogger << "(MAIN): Saving initial state of variables to disk." << endl << write;
      }
      //writeGrid(mpiGrid,reducer);
      if (writeSpatialCellData(mpiGrid,vlsWriter,reducer) == false) {
	 mpilogger << "(MAIN): ERROR occurred while writing data to file!" << endl << write;
      }
   }
   initIoTime=MPI_Wtime()-initIoTime;
   
   if (P::save_velocity_grid) {
      writeAllVelocityBlocks(mpiGrid);
   }
   writeSomeVelocityGrids(mpiGrid, P::save_spatial_cells_x, P::save_spatial_cells_y, P::save_spatial_cells_z);
   #ifndef PARGRID
      comm.barrier();
   #else
      mpiGrid.barrier();
   #endif
   
   // Main simulation loop:
   if (myrank == MASTER_RANK) 
     mpilogger << "(MAIN): Starting main simulation loop." << endl << write;

   double before = MPI_Wtime();
   for (luint tstep=P::tstep_min; tstep < P::tsteps; ++tstep) {
#ifdef CRAYPAT
//turn on & off sampling & tracing
       if(myrank>=firstrank && myrank<=lastrank){
           if(tstep==firstiter) PAT_state(PAT_STATE_ON);
           if(tstep>lastiter) PAT_state(PAT_STATE_OFF);
       }
#endif 
       // Recalculate (maybe) spatial cell parameters
      calculateSimParameters(mpiGrid, P::t, P::dt);

      // use globally minimum timestep
      #ifndef PARGRID
      P::dt = all_reduce(comm, P::dt, boost::mpi::minimum<Real>());
      #else
      #warning Cannot calculate minimum timestep when using PARGRID: no communicator for all_reduce
      #endif

      // Propagate the state of simulation forward in time by dt:
#ifdef CRAYPAT
      PAT_region_begin(1,"calculateAcceleration");
#endif 
      calculateAcceleration(mpiGrid);
#ifdef CRAYPAT
      PAT_region_end(1);
      PAT_region_begin(2,"calculateSpatialDerivatives");
#endif 
      calculateSpatialDerivatives(mpiGrid);
#ifdef CRAYPAT
      PAT_region_end(2);
      PAT_region_begin(3,"calculateSpatialFluxes");
#endif 
      calculateSpatialFluxes(mpiGrid);
#ifdef CRAYPAT
      PAT_region_end(3);
      PAT_region_begin(4,"calculateSpatialPropagation");
#endif 
      calculateSpatialPropagation(mpiGrid);
#ifdef CRAYPAT
      PAT_region_end(4);
#endif 
      ++P::tstep;
      P::t += P::dt;
      
      // Check if the full simulation state should be written to disk
      double t1=MPI_Wtime();
      if (P::tstep % P::saveRestartInterval == 0) {
         // TODO: implement full state saving
	 if (myrank == MASTER_RANK) {
	    mpilogger << "(MAIN): Writing restart files to disk at tstep = " << P::tstep << ", time = " << P::t << endl;
	    mpilogger << "\t NOT IMPLEMENTED YET" << endl << write;
	 }
      }

      // Check if variables and derived quantities should be written to disk
      if (P::tstep % P::diagnInterval == 0) {
	 if (myrank == MASTER_RANK) {
	    mpilogger << "(MAIN): Saving variables to disk at tstep = " << P::tstep << ", time = " << P::t << endl << write;
	 }
	 if (P::save_spatial_grid) {
	    if (writeSpatialCellData(mpiGrid,vlsWriter,reducer) == false) {
	       mpilogger << "(MAIN): ERROR occurred while writing data to file!" << endl << write;
	    }
	 }
	 
         if (P::save_velocity_grid) {
            writeAllVelocityBlocks(mpiGrid);
         }
         writeSomeVelocityGrids(mpiGrid, P::save_spatial_cells_x, P::save_spatial_cells_y, P::save_spatial_cells_z);
      }
      loopIoTime+=MPI_Wtime()-t1;
          
      #ifndef PARGRID
         comm.barrier();
      #else
         mpiGrid.barrier();
      #endif
   }
   double after = MPI_Wtime();

   double finalIoTime=MPI_Wtime();
   // Write final state:
   if (P::save_spatial_grid) {
      mpilogger << "(MAIN): Saving variables to disk at tstep = " << P::tstep << ", time = " << P::t << endl << write;
      if (writeSpatialCellData(mpiGrid,vlsWriter,reducer) == false) {
	 mpilogger << "(MAIN): ERROR occurred while writing data to file!" << endl << write;
      }
   }
   if (P::save_velocity_grid) {
      writeAllVelocityBlocks(mpiGrid);
   }
   writeSomeVelocityGrids(mpiGrid, P::save_spatial_cells_x, P::save_spatial_cells_y, P::save_spatial_cells_z);

   finalIoTime=MPI_Wtime()-finalIoTime;
   totTime=MPI_Wtime()-totTime;
   
   if (myrank == MASTER_RANK) {
      mpilogger << "(MAIN): All timesteps calculated." << endl;
      mpilogger << "(TIME) total simulated time " << P::t << " s" << endl;
      mpilogger << "(TIME) Total time " << totTime << " s" << endl;
      mpilogger << "(TIME)   Initialization time " << initTime << " s" << endl;
      mpilogger << "(TIME)   Initial IO time " << initIoTime << " s" << endl;
      mpilogger << "(TIME)   Main loop: total run time  " << after - before << " s" << endl;
      mpilogger << "(TIME)     Main loop: Compute+MPI " << after - before - loopIoTime << " s" << endl;
      mpilogger << "(TIME)     Main loop: IO " << loopIoTime << " s" << endl;
      mpilogger << "(TIME)   Final IO time " << finalIoTime << " s" << endl;
      mpilogger << write;
   }




   // Write the timer values (if timers have been defined):
   writeTimers();
   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Exiting." << endl << write;
   mpilogger.close();
   return 0;
}

