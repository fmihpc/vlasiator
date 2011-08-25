#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>

#include "vlasovmover.h"
#include "definitions.h"
#include "mpiconversion.h"
#include "mpilogger.h"
#include "parameters.h"
#include "grid.h"
#include "silowriter.h"
#include "cell_spatial.h"
#include "gridbuilder.h"
#include "datareducer.h"
#include "datareductionoperator.h"
#include "transferstencil.h"

#include "vlsvwriter2.h" // TEST
#include "fieldsolver.h"

#ifdef CRAYPAT
//include craypat api headers if compiled with craypat on Cray XT/XE
#include "pat_api.h"
#endif

#include "profile.h"

Grid grid;
MPILogger mpilogger;

bool inistate = true;

using namespace std;

#ifndef PARGRID
void initSpatialCells(const dccrg<SpatialCell>& mpiGrid,boost::mpi::communicator& comm) {
    typedef Parameters P;

   vector<uint64_t> cells = mpiGrid.get_cells();
   
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
      if (!mpiGrid[cells[i]]->initialize(mpiGrid[cells[i]]->N_blocks)) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " Failed to initialize spatial cell " << cells[i]
            << ", too small MAX_VEL_BLOCKS in parameters.h?"
            << std::endl;
         abort();
      }
      buildSpatialCell(*(mpiGrid[cells[i]]),xmin,ymin,zmin,dx,dy,dz,false);
   }
}
#endif

#ifdef PARGRID
bool writeGrid(const ParGrid<SpatialCell>& mpiGrid,DataReducer& dataReducer,const bool& writeRestart) {
#else
bool writeGrid(const dccrg<SpatialCell>& mpiGrid,DataReducer& dataReducer,const bool& writeRestart) {
#endif
    
    double allStart = MPI_Wtime();
    bool success = true;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(writeRestart)
        profile::start("writeGrid-restart");
    else
        profile::start("writeGrid-reduced");
    
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
      vector<ID::type> cells;
      mpiGrid.getCells(cells);
   #else 
      vector<uint64_t> cells = mpiGrid.get_cells();
   #endif

   if (vlsvWriter.writeArray("MESH","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) {
      cerr << "Proc #" << myrank << " failed to write cell IDs!" << endl;
   }

   // Create a buffer for spatial cell coordinates. Copy all coordinates to 
   // buffer and write:
   Real* buffer = new Real[6*cells.size()];
   for (size_t i=0; i<cells.size(); ++i) {
      SpatialCell* SC = mpiGrid[cells[i]];
      for (int j=0; j<6; ++j) {
	 buffer[6*i+j] = SC->cpu_cellParams[j];
      }
   }
   if (vlsvWriter.writeArray("COORDS","SpatialGrid",attribs,cells.size(),6,buffer) == false) {
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
      uint64_t arraySize = cells.size()*vectorSize*dataSize;
      
      // Request DataReductionOperator to calculate the reduced data for all local cells:
      varBuffer = new char[arraySize];
      for (uint64_t cell=0; cell<cells.size(); ++cell) {
	 if (dataReducer.reduceData(mpiGrid[cells[cell]],i,varBuffer + cell*vectorSize*dataSize) == false) success = false;
      }
      if (success == false) mpilogger << "(MAIN) writeGrid: ERROR datareductionoperator '" << dataReducer.getName(i) << "' returned false!" << endl << write;
      
      // Write reduced data to file:
      if (vlsvWriter.writeArray("VARIABLE",variableName,attribs,cells.size(),vectorSize,dataType,dataSize,varBuffer) == false) success = false;
      if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to write datareductionoperator data to file!" << endl << write;
      delete varBuffer;
      varBuffer = NULL;
   }

   // If restart data is not written, exit here:
   if (writeRestart == false) {
      MPI_Barrier(MPI_COMM_WORLD);
      vlsvWriter.close();
      profile::stop("writeGrid-reduced");
      return success;
   }
   
   attribs.clear();

   // Write spatial cell parameters:
   Real* paramsBuffer = new Real[cells.size()*SIZE_CELLPARAMS];
   for (size_t i=0; i<cells.size(); ++i) for (uint j=0; j<SIZE_CELLPARAMS; ++j) 
     paramsBuffer[i*SIZE_CELLPARAMS+j] = mpiGrid[cells[i]]->cpu_cellParams[j];   
   if (vlsvWriter.writeArray("CELLPARAMS","SpatialGrid",attribs,cells.size(),SIZE_CELLPARAMS,paramsBuffer) == false) {
      mpilogger << "(MAIN) writeGrid: ERROR failed to write spatial cell parameters!" << endl << write;
      success = false;
   }
   delete paramsBuffer;
   
   // Write the number of spatial neighbours each cell has:
   uchar* N_neighbours = new uchar[cells.size()];
   uint64_t neighbourSum = 0;
   for (size_t i=0; i<cells.size(); ++i) {
      #ifdef PARGRID
         N_neighbours[i] = mpiGrid.getNumberOfNeighbours(cells[i]);
      #else
      
      #endif
      neighbourSum += N_neighbours[i];
   }
   if (vlsvWriter.writeArray("NBRSUM","SpatialGrid",attribs,cells.size(),1,N_neighbours) == false) success = false;
   delete N_neighbours;
   
   // Write global IDs of spatial neighbours:
   #ifdef PARGRID
      ID::type* nbrIDs = new ID::type[neighbourSum];
      uchar* nbrTypeIDs = new uchar[neighbourSum];
      vector<uchar> tmpTypes;
      vector<uint> tmpIDs;
      size_t counter = 0;
      for (size_t i=0; i<cells.size(); ++i) {
	 mpiGrid.getExistingNeighbours(cells[i],tmpIDs,tmpTypes);
	 for (size_t j=0; j<tmpTypes.size(); ++j) {
	    nbrTypeIDs[counter] = tmpTypes[j];
	    nbrIDs[counter] = tmpIDs[j];
	    ++counter;
	 }
      }
      if (vlsvWriter.writeArray("NBRTYPES","SpatialGrid",attribs,neighbourSum,1,nbrTypeIDs) == false) success = false;
      if (vlsvWriter.writeArray("NBRIDS","SpatialGrid",attribs,neighbourSum,1,nbrIDs) == false) success = false;
      delete nbrIDs;
      delete nbrTypeIDs;
   #else

   #endif
   
   // Write velocity blocks and related data. Which cells write velocity grids 
   // should be requested from a function, but for now we just write velocity grids for all cells.
   // First write global IDs of those cells which write velocity blocks (here: all cells):
   if (vlsvWriter.writeArray("CELLSWITHBLOCKS","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) success = false;
   if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << write;
   
   // Write the number of velocity blocks in each spatial cell. Again a temporary buffer is used:
   uint* N_blocks = new uint[cells.size()];
   uint64_t totalBlocks = 0;
   for (size_t cell=0; cell<cells.size(); ++cell) {
      N_blocks[cell] = mpiGrid[cells[cell]]->N_blocks;
      totalBlocks += mpiGrid[cells[cell]]->N_blocks;
   }
   if (vlsvWriter.writeArray("NBLOCKS","SpatialGrid",attribs,cells.size(),1,N_blocks) == false) success = false;
   if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to write NBLOCKS to file!" << endl << write;

   double start = MPI_Wtime();
   
   // Write velocity block coordinates. At this point the data may get too large to be buffered, 
   // so a "multi-write" mode is used - coordinates are written one velocity grid at a time. Besides, 
   // the velocity block coordinates are already stored in a suitable format for writing in SpatialCells:
   if (vlsvWriter.startMultiwrite("float",totalBlocks,SIZE_BLOCKPARAMS,sizeof(Real)) == false) success = false;
   if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to start BLOCKCOORDINATES multiwrite!" << endl << write;
   if (success == true) {
      SpatialCell* SC;
      for (size_t cell=0; cell<cells.size(); ++cell) {
	 SC = mpiGrid[cells[cell]];
	 if (vlsvWriter.multiwriteArray(N_blocks[cell],SC->cpu_blockParams) == false) success = false;
      }
   }
   if (success == true) {
      if (vlsvWriter.endMultiwrite("BLOCKCOORDINATES","SpatialGrid",attribs) == false) success = false;
      if (success == false) mpilogger << "(MAIN) writeGrid: ERROR occurred when ending BLOCKCOORDINATES multiwrite!" << endl << write;
   }
   
   // Write values of distribution function:
   if (vlsvWriter.startMultiwrite("float",totalBlocks,SIZE_VELBLOCK,sizeof(Real)) == false) success = false;
   if (success == false) mpilogger << "(MAIN) writeGrid: ERROR failed to start BLOCKVARIABLE avgs multiwrite!" << endl << write;
   if (success == true) {
      uint64_t counter = 0;
      SpatialCell* SC;
      for (size_t cell=0; cell<cells.size(); ++cell) {
	 SC = mpiGrid[cells[cell]];
	 if (vlsvWriter.multiwriteArray(N_blocks[cell],SC->cpu_avgs) == false) success = false;
	 //if (vlsvWriter.multiwriteArray(counter,N_blocks[cell],SC->cpu_fx) == false) success = false;
	 counter += N_blocks[cell];
      }
   }
   attribs["mesh"] = "SpatialGrid";
   if (success == true) {
      if (vlsvWriter.endMultiwrite("BLOCKVARIABLE","avgs",attribs) == false) success = false;
      if (success == false) mpilogger << "(MAIN) writeGrid: ERROR occurred when ending BLOCKVARIABLE avgs multiwrite!" << endl << write;
   }
   
   // Write velocity block neighbour lists:
   if (vlsvWriter.startMultiwrite("uint",totalBlocks,SIZE_NBRS_VEL,sizeof(uint)) == false) {
      success = false;
      mpilogger << "(MAIN) writeGrid: ERROR failed to start BLOCKNBRS multiwrite!" << endl << write;
   } else {
      SpatialCell* SC;
      for (size_t cell=0; cell<cells.size(); ++cell) {
	 SC = mpiGrid[cells[cell]];
	 if (vlsvWriter.multiwriteArray(N_blocks[cell],SC->cpu_nbrsVel) == false) success = false;
      }
   }
   if (success == true) {
      if (vlsvWriter.endMultiwrite("BLOCKNBRS","SpatialGrid",attribs) == false) {
	 success = false;
	 mpilogger << "(MAIN) writeGrid: ERROR occurred when ending BLOCKNBRS multiwrite!" << endl << write;
      }
   }
   double end = MPI_Wtime();
   delete N_blocks;

   vlsvWriter.close();

   double allEnd = MPI_Wtime();
   
   double bytesWritten = totalBlocks*((SIZE_VELBLOCK+SIZE_BLOCKPARAMS)*sizeof(Real)+(SIZE_NBRS_VEL)*sizeof(uint));
   double secs = end-start;
   mpilogger << "Wrote " << bytesWritten/1.0e6 << " MB of data in " << secs << " seconds, datarate is " << bytesWritten/secs/1.0e9 << " GB/s" << endl << write;
   
   double allSecs = allEnd-allStart;
   if (myrank == 0) mpilogger << "All data written in " << allSecs << " seconds" << endl << write;

   profile::stop("writeGrid-restart",1.0e-6*bytesWritten,"MB");
   return success;
}

void exchangeVelocityGridMetadata(
	#ifdef PARGRID
	ParGrid<SpatialCell>& mpiGrid
	#else
	dccrg<SpatialCell>& mpiGrid
	#endif
) {
   #ifdef PARGRID
   vector<uchar> nbrTypes;
   for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
      if ((i == 0) & ((j== 0) & (k == 0))) continue;
      const int ii = i+2;
      const int jj = j+2;
      const int kk = k+2;
      nbrTypes.push_back(kk*25+jj*5+ii);
   }
   nbrTypes.push_back((2+0)*25+(2+0)*5+(2-2)); // --x 
   nbrTypes.push_back((2+0)*25+(2+0)*5+(2+2)); // ++x
   nbrTypes.push_back((2+0)*25+(2-2)*5+(2-0)); // --y
   nbrTypes.push_back((2+0)*25+(2+2)*5+(2+0)); // ++y
   nbrTypes.push_back((2-2)*25+(2+0)*5+(2-0)); // --z
   nbrTypes.push_back((2+2)*25+(2+0)*5+(2+0)); // ++z
   
   TransferStencil<ID::type> stencil(numeric_limits<ID::type>::max());
   stencil.addReceives(mpiGrid,nbrTypes);
   stencil.addSends(mpiGrid,nbrTypes);
   
   // Post receives for N_blocks
   mpiGrid.startSingleMode();
   for (map<pair<int,int>,ID::type>::const_iterator i=stencil.recvs.begin(); i!=stencil.recvs.end(); ++i) {
      const int host        = i->first.first;
      const int tag         = i->first.second;
      const ID::type cellID = i->second;
      char* buffer          = reinterpret_cast<char*>(&(mpiGrid[cellID]->N_blocks));
      mpiGrid.singleReceive(host,tag,sizeof(uint),buffer,cellID);
   }
   // Post sends for N_blocks
   for (multimap<ID::type,pair<int,int> >::const_iterator i=stencil.sends.begin(); i!=stencil.sends.end(); ++i) {
      const ID::type cellID = i->first;
      const int dest        = i->second.first;
      const int tag         = i->second.second;
      char* buffer          = reinterpret_cast<char*>(&(mpiGrid[cellID]->N_blocks));
      mpiGrid.singleSend(dest,tag,sizeof(uint),buffer,cellID);
   }
   mpiGrid.waitAllReceives();
   mpiGrid.waitAllSends();

   #else	// ifdef PARGRID

   SpatialCell::base_address_identifier = 5;
   mpiGrid.update_remote_neighbour_data();

   #endif	// ifdef PARGRID
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

   profile::start("main");
   profile::start("Initialization");

   profile::start("Read parameters");
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
   RP::add("propagate_field","Propagate magnetic field during the simulation",true);
   RP::add("propagate_vlasov","Propagate distribution functions during the simulation",true);
   
   RP::parse();
   RP::get("solar_wind_file",P::solar_wind_file);
   RP::get("save_interval", P::diagnInterval);
   RP::get("restart_interval", P::saveRestartInterval);
   RP::get("save_spatial_grid", P::save_spatial_grid);
   RP::get("save_velocity_grid", P::save_velocity_grid);
   RP::get("save_spatial_cells_at_x,X", P::save_spatial_cells_x);
   RP::get("save_spatial_cells_at_y,Y", P::save_spatial_cells_y);
   RP::get("save_spatial_cells_at_z,Z", P::save_spatial_cells_z);
   RP::get("propagate_field",P::propagateField);
   RP::get("propagate_vlasov",P::propagateVlasov);
   
   // Sanity checks (DEPRECATED):
   if (P::save_spatial_cells_x.size() != P::save_spatial_cells_y.size()
       || P::save_spatial_cells_x.size() != P::save_spatial_cells_z.size()) {
       cerr << "Must have equal number of values for x, y and z components, but given: x " << P::save_spatial_cells_x.size();
       cerr << ", y " << P::save_spatial_cells_y.size() << ", z " << P::save_spatial_cells_z.size() << endl;
       MPI_Finalize();
       exit(1);
   }

   profile::stop("Read parameters");

#ifdef CRAYPAT
/*initialize variables for reducing sampling & tracing to a slice of all iterations.
*/
   char *envvar;
   int firstiter=-1;
   int lastiter=numeric_limits<int>::max();
#warning Including CrayPAT API 
   if(myrank==0){
       envvar=getenv("CRAYPAT_FIRSTITER");
       if(envvar!=NULL){
           firstiter=atoi(envvar);
       }
       envvar=getenv("CRAYPAT_LASTITER");
       if(envvar!=NULL){
           lastiter=atoi(envvar);
       }
       
   }
   
   MPI_Bcast(&firstiter,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&lastiter,1,MPI_INT,0,MPI_COMM_WORLD);

   //turn off craypat sampling & tracing based on the defined slices.
   if(firstiter>0){
       PAT_state(PAT_STATE_OFF);
   }
   
#endif
   
   profile::start("open mpilogger");
   // Init parallel logger:
   if (mpilogger.open(MPI_COMM_WORLD,MASTER_RANK,"logfile.txt") == false) {
      cerr << "(MAIN) ERROR: MPILogger failed to open output file!" << endl;
      exit(1);
   }
   profile::stop("open mpilogger");
   profile::start("Initialize Grid");

#ifndef PARGRID // INITIALIZE USING DCCRG
      // Create parallel MPI grid and init Zoltan:
      float zoltanVersion;
      if (Zoltan_Initialize(argn,args,&zoltanVersion) != ZOLTAN_OK) {
	 mpilogger << "\t ERROR: Zoltan initialization failed, aborting." << std::endl << write;
	 success = false;
      } else {
	 mpilogger << "\t Zoltan " << zoltanVersion << " initialized successfully" << std::endl << write;
      }
      if (buildGrid(MPI_COMM_WORLD,MASTER_RANK) == false) {
	 mpilogger << "(MAIN) Grid builder failed!" << endl << write;
	 success = false;
      } else {
	 mpilogger << "(MAIN) Grid built successfully" << endl << write;
      }

      
      dccrg<SpatialCell> mpiGrid(
         comm,
         "HYPERGRAPH",
         P::xmin, P::ymin, P::zmin,
         P::dx_ini, P::dy_ini, P::dz_ini,
         P::xcells_ini, P::ycells_ini, P::zcells_ini,
#ifdef SOLVER_KT
         2, // neighborhood size
#elif defined SOLVER_LEVEQUE
         4, //neighborhood size
#endif
         0, // maximum refinement level
         P::periodic_x, P::periodic_y, P::periodic_z
      );
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
      for(unsigned int i=0;i<partitionProcs.size();i++){
          // set hierarchial partitioning parameters for first level
          if(myrank==0){
              mpilogger << "(MAIN) Partition parameters level "<<i <<" procs: " << partitionProcs[i] << " LB_METHOD: " <<
                  partitionLbMethod[i]<< " IMBALANCE_TOL: "<< partitionImbalanceTol[i] << endl << write;
          }
          
          mpiGrid.add_partitioning_level(partitionProcs[i]);
          mpiGrid.add_partitioning_option(i, "LB_METHOD", partitionLbMethod[i]);
          mpiGrid.add_partitioning_option(i, "IMBALANCE_TOL", partitionImbalanceTol[i]);
      }

      profile::start("Initial load-balancing");
      if (myrank == MASTER_RANK) mpilogger << "(MAIN): Starting initial load balance." << endl << write;
      initialLoadBalance(mpiGrid);
      profile::stop("Initial load-balancing");
   
   
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
      profile::stop("Initialize Grid");

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

   // Receive N_blocks from remote neighbours
   exchangeVelocityGridMetadata(mpiGrid);
   // reserve space for velocity blocks in local copies of remote neighbors
   #ifdef PARGRID
      vector<ID::type> remoteCells;
      mpiGrid.getRemoteCells(remoteCells);
      for (size_t i=0; i<remoteCells.size(); ++i) mpiGrid[remoteCells[i]]->initialize(mpiGrid[remoteCells[i]]->N_blocks);
   #else
   const boost::unordered_set<uint64_t>* incoming_cells = mpiGrid.get_remote_cells_with_local_neighbours();
   for (boost::unordered_set<uint64_t>::const_iterator
      cell_id = incoming_cells->begin();
      cell_id != incoming_cells->end();
      cell_id++
   ) {
      SpatialCell* cell = mpiGrid[*cell_id];
      if (cell == NULL) {
         cerr << "No data for spatial cell " << *cell_id << endl;
         abort();
      }
      cell->initialize(cell->N_blocks);
   }
   #endif
   
   // Go through every spatial cell on this CPU, and create the initial state:

#ifndef PARGRID
   profile::start("Set initial state");
   initSpatialCells(mpiGrid,comm);
   profile::stop("Set initial state");
#endif

  

   profile::start("Fetch Neighbour data");
   // Fetch neighbour data:
   #ifndef PARGRID
      // FIXME: add a mpi_datatype to send all cell data
      P::transmit = 5;
      mpiGrid.update_remote_neighbour_data(); // TEST
   #endif
   profile::stop("Fetch Neighbour data");
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
   reducer.addOperator(new DRO::VariableVolE);
   reducer.addOperator(new DRO::VariableVolB);
   
   //VlsWriter vlsWriter;
   profile::start("Init vlasov propagator");
   // Initialize Vlasov propagator:
   if (initializeMover(mpiGrid) == false) {
      mpilogger << "(MAIN): Vlasov propagator did not initialize correctly!" << endl << write;
      exit(1);
   }
   calculateVelocityMoments(mpiGrid);
   profile::stop("Init vlasov propagator");
   
   profile::start("Init field propagator");
   // Initialize field propagator:
   if (initializeFieldPropagator(mpiGrid) == false) {
      mpilogger << "(MAIN): Field propagator did not initialize correctly!" << endl << write;
      exit(1);
   }
   calculateVolumeAveragedFields(mpiGrid);
   profile::stop("Init field propagator");

   // ***********************************
   // ***** INITIALIZATION COMPLETE *****
   // ***********************************
   
   // Free up memory:
   readparameters.finalize();
   profile::start("Save initial state");
   // Write initial state:
   if (P::save_spatial_grid) {
      if (myrank == MASTER_RANK) {
	 mpilogger << "(MAIN): Saving initial state of variables to disk." << endl << write;
      }

      //writegrid has new vlsvwriter routines
      if (writeGrid(mpiGrid,reducer,true) == false) {
	 mpilogger << "(MAIN): ERROR occurred while writing spatial cell and restart data!" << endl << write;
      }
   }
   profile::stop("Save initial state");
   profile::stop("Initialization");
#ifndef PARGRID
   comm.barrier();
#else
   mpiGrid.barrier();
#endif
   inistate = false;
   // Main simulation loop:
   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Starting main simulation loop." << endl << write;

   double before = MPI_Wtime();
   unsigned int totalComputedSpatialCells=0;
   unsigned int computedSpatialCells=0;
   profile::start("Simulation");
   for (luint tstep=P::tstep_min; tstep < P::tsteps; ++tstep) {
#ifdef CRAYPAT //turn on & off sampling & tracing
       if(myrank>=firstrank && myrank<=lastrank){
           if(tstep==firstiter) PAT_state(PAT_STATE_ON);
           if(tstep>lastiter) PAT_state(PAT_STATE_OFF);
       }
#endif

       
       //compute how many spatial cells we solve for this step
#ifndef PARGRID
       computedSpatialCells=mpiGrid.get_cells().size();
#else
       vector<uint64_t> cells;
       mpiGrid.getCells(cells);
       computedSpatialCells=cells.size();
#endif
       
       totalComputedSpatialCells+=computedSpatialCells;
       profile::start("Propagate");
       // Recalculate (maybe) spatial cell parameters
       calculateSimParameters(mpiGrid, P::t, P::dt);

      // use globally minimum timestep
      #ifndef PARGRID
      P::dt = all_reduce(comm, P::dt, boost::mpi::minimum<Real>());
      #else
      #warning Cannot calculate minimum timestep when using PARGRID: no communicator for all_reduce
      #endif

      // Propagate the state of simulation forward in time by dt:      
      if (P::propagateVlasov == true) {
          profile::start("Propagate Vlasov");
          profile::start("First propagation");
          calculateSpatialDerivatives(mpiGrid);
          calculateSpatialFluxes(mpiGrid);
          calculateSpatialPropagation(mpiGrid,false,false);
          profile::stop("First propagation",computedSpatialCells,"SpatialCells");
          bool transferAvgs = false;
          if (P::tstep % P::saveRestartInterval == 0 ||
              P::tstep == P::tsteps-1) transferAvgs = true;
          
	  profile::start("Acceleration");
          calculateAcceleration(mpiGrid);
	  profile::stop("Acceleration",computedSpatialCells,"SpatialCells");
	 
          profile::start("Second propagation");
          calculateSpatialDerivatives(mpiGrid);
          calculateSpatialFluxes(mpiGrid);
          calculateSpatialPropagation(mpiGrid,true,transferAvgs);
          profile::stop("Second propagation",computedSpatialCells,"SpatialCells");
          profile::stop("Propagate Vlasov",computedSpatialCells,"SpatialCells");
      }

      // Propagate fields forward in time by dt. If field is not 
      // propagated self-consistently (test-Vlasov simulation), then 
      // re-calculate face-averaged E,B fields. This requires that 
      // edge-E and face-B have been shared with remote neighbours 
      // (not done by calculateFaceAveragedFields).

      if (P::propagateField == true) {
          profile::start("Propagate Fields");
          propagateFields(mpiGrid,P::dt);
          profile::stop("Propagate Fields",computedSpatialCells,"SpatialCells");
      } else {
	 calculateFaceAveragedFields(mpiGrid);
      }
      profile::stop("Propagate",computedSpatialCells,"SpatialCells");
      ++P::tstep;
      P::t += P::dt;
      
      // Check if data needs to be written to disk:
      if (P::tstep % P::saveRestartInterval == 0 || P::tstep % P::diagnInterval == 0) {
          profile::start("IO");
	 bool writeRestartData = false;
	 if (P::tstep % P::saveRestartInterval == 0) {
	   writeRestartData = true;
	   if (myrank == MASTER_RANK)
	   mpilogger << "(MAIN): Writing spatial cell and restart data to disk, tstep = " << P::tstep << " t = " << P::t << endl << write;
	 } else
	   if (myrank == MASTER_RANK)
	     mpilogger << "(MAIN): Writing spatial cell data to disk, tstep = " << P::tstep << " t = " << P::t << endl << write;
	 
	 calculateVolumeAveragedFields(mpiGrid);
	 if (writeGrid(mpiGrid,reducer,writeRestartData) == false) {
	    if (myrank == MASTER_RANK)
	      mpilogger << "(MAIN): ERROR occurred while writing spatial cell and restart data!" << endl << write;
	 }
         profile::stop("IO");
      }
      
      MPI_Barrier(MPI_COMM_WORLD);
   }
   double after = MPI_Wtime();
   profile::stop("Simulation",totalComputedSpatialCells,"SpatialCells");
   profile::start("Finalization");   
   finalizeMover();
#ifdef PARGRID
   finalizeFieldPropagator(mpiGrid);
#endif
   
   if (myrank == MASTER_RANK) {
       mpilogger << "(MAIN): All timesteps calculated." << endl;
       mpilogger << "\t (TIME) total run time " << after - before << " s, total simulated time " << P::t << " s" << endl;
       mpilogger << "\t (TIME) seconds per timestep " << double(after - before) / P::tsteps <<
           ", seconds per simulated second " << double(after - before) / P::t << endl;
       mpilogger << write;
   }
   
   profile::stop("Finalization");   
   profile::stop("main");
   profile::print(MPI_COMM_WORLD);
   

   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Exiting." << endl << write;
   mpilogger.close();
   return 0;
}

