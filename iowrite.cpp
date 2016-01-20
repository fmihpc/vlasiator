/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
*/

/*! \file vlsvdiff.cpp
 \brief File containing write IO for vlasiator. More info at: https://agora.fmi.fi/display/CORSAIR/VLSV+File+Format
*/

#include <cstdlib>
#include <iostream>
#include <iomanip> // for setprecision()
#include <cmath>
#include <sstream>
#include <ctime>
#include <array>
#include "iowrite.h"
#include "grid.h"
#include "phiprof.hpp"
#include "parameters.h"
#include "logger.h"
#include "vlasovmover.h"

using namespace std;
using namespace phiprof;
using namespace vlsv;

extern Logger logFile, diagnostic;

typedef Parameters P;


/*! Updates local ids across MPI to let other processes know in which order this process saves the local cell ids
 \param mpiGrid Vlasiator's MPI grid
 \param local_cells local cells on in the current process (no ghost cells included)
 \param comm The MPi comm
 \return Returns true if operation was successful
 */
bool updateLocalIds(  dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> & mpiGrid,
                      const vector<uint64_t> & local_cells,
                      MPI_Comm comm ) {
   int myRank;
   MPI_Comm_rank(comm,&myRank);

   //Declare an iterator for iterating though the cell ids
   vector<uint64_t>::const_iterator it;
   //Local ids for the process start from 0 (this is used in the iteration)
   uint64_t thisProcessLocalId = 0;
   //Iterate through local cells
   for( it = local_cells.begin(); it != local_cells.end(); ++it ) {
      //NOTE: (*it) = cellId
      //Set the local id
      mpiGrid[(*it)]->ioLocalCellId = thisProcessLocalId;
      //Increment the local id
      thisProcessLocalId++;
   }
   //Update the local ids (let the other processes know they've been updated)
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_IOLOCALCELLID);
   mpiGrid.update_copies_of_remote_neighbors(NEAREST_NEIGHBORHOOD_ID);

   return true;
}

/*! Checks the success of write IO
 \param success Parameter for determining whether the IO has been successful
 \param errorMessage The message to be shown if global success is false
 \param comm The MPI comm
 \return Returns true if all of the processes running this function had success value of true, false otherwise
 */
bool globalSuccess(bool success,string errorMessage,MPI_Comm comm){
   int successInt;
   int globalSuccessInt;
   if(success)
      successInt=1;
   else
      successInt=0;
   
   MPI_Allreduce(&successInt,&globalSuccessInt,1,MPI_INT,MPI_MIN,comm);
   
   if(globalSuccessInt==1) {
      return true;
   }
   else{
      logFile << errorMessage << endl<<write ;
      return false;
   }
}

/** Writes the velocity distribution into the file.
 @param vlsvWriter Some vlsv writer with a file open.
 @param mpiGrid Vlasiator's grid.
 @param cells Vector of local cells within this process (no ghost cells).
 @param comm The MPI communicator.
 @return Returns true if operation was successful.*/
bool writeVelocityDistributionData(Writer& vlsvWriter,
                                   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                   const vector<uint64_t> & cells,MPI_Comm comm) {
   // Write velocity blocks and related data. 
   // In restart we just write velocity grids for all cells.
   // First write global Ids of those cells which write velocity blocks (here: all cells):
   map<string,string> attribs;
   const string popName      = "avgs";
   const string spatMeshName = "SpatialGrid";
   attribs["name"] = popName;
   bool success=true;

   // Compute totalBlocks
   uint64_t totalBlocks = 0;  
   vector<uint> blocksPerCell;   
   for (size_t cell=0; cell<cells.size(); ++cell){
      totalBlocks+=mpiGrid[cells[cell]]->get_number_of_velocity_blocks();
      blocksPerCell.push_back(mpiGrid[cells[cell]]->get_number_of_velocity_blocks());
   }

   // The name of the mesh is "SpatialGrid"
   attribs["mesh"] = spatMeshName;
   const unsigned int vectorSize = 1;
   // Write the array:
   if (vlsvWriter.writeArray("CELLSWITHBLOCKS",attribs,cells.size(),vectorSize,cells.data()) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;
   // Write blocks per cell, this has to be in the same order as cellswitblocks so that extracting works
   if(vlsvWriter.writeArray("BLOCKSPERCELL",attribs,blocksPerCell.size(),vectorSize,blocksPerCell.data()) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;

   // Write (partial) velocity mesh data
   const uint8_t refLevel = 0; // Mesh refinement level, 0 here since velocity mesh is unrefined
   uint64_t bbox[6];
   bbox[0] = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getGridLength(refLevel)[0];
   bbox[1] = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getGridLength(refLevel)[1];
   bbox[2] = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getGridLength(refLevel)[2];
   bbox[3] = WID;
   bbox[4] = WID;
   bbox[5] = WID;

   attribs.clear();
   attribs["mesh"] = popName;
   attribs["type"] = vlsv::mesh::STRING_UCD_AMR;
   if (mpiGrid.get_rank() == MASTER_RANK) {
      if (vlsvWriter.writeArray("MESH_BBOX",attribs,6,1,bbox) == false) success = false;
      
      for (int crd=0; crd<3; ++crd) {
         const size_t N_nodes = bbox[crd]*bbox[crd+3]+1;
         Real* crds = new Real[N_nodes];
         const Real dV = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getCellSize(refLevel)[crd];

         for (size_t i=0; i<N_nodes; ++i) {
            crds[i] = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getMeshMinLimits()[crd] + i*dV;
         }
         
         if (crd == 0) {
            if (vlsvWriter.writeArray("MESH_NODE_CRDS_X",attribs,N_nodes,1,crds) == false) success = false;
         }
         if (crd == 1) {
            if (vlsvWriter.writeArray("MESH_NODE_CRDS_Y",attribs,N_nodes,1,crds) == false) success = false;
         }
         if (crd == 2) {
            if (vlsvWriter.writeArray("MESH_NODE_CRDS_Z",attribs,N_nodes,1,crds) == false) success = false;
         }
         delete [] crds; crds = NULL;
      }
   } else {
      if (vlsvWriter.writeArray("MESH_BBOX",attribs,0,1,bbox) == false) success = false;
      Real* crds = NULL;
      if (vlsvWriter.writeArray("MESH_NODE_CRDS_X",attribs,0,1,crds) == false) success = false;
      if (vlsvWriter.writeArray("MESH_NODE_CRDS_Y",attribs,0,1,crds) == false) success = false;
      if (vlsvWriter.writeArray("MESH_NODE_CRDS_Z",attribs,0,1,crds) == false) success = false;
   }
   
   // Write velocity block ids
   vector<vmesh::GlobalID> velocityBlockIds;
   try {
      velocityBlockIds.reserve( totalBlocks );
      // gather data for writing
      for (size_t cell=0; cell<cells.size(); ++cell) {
         SpatialCell* SC = mpiGrid[cells[cell]];
         for (vmesh::LocalID block_i=0; block_i<SC->get_number_of_velocity_blocks(); ++block_i) {
            vmesh::GlobalID block = SC->get_velocity_block_global_id(block_i);
            velocityBlockIds.push_back( block );
         }
      }
   } catch (...) {
      cerr << "FAILED TO WRITE VELOCITY BLOCK IDS AT: " << __FILE__ << " " << __LINE__ << endl;
      success=false;
   }

   if (globalSuccess(success,"(MAIN) writeGrid: ERROR: Failed to fill temporary array velocityBlockIds",MPI_COMM_WORLD) == false) {
      vlsvWriter.close();
      return false;
   }

   attribs.clear();
   attribs["mesh"] = spatMeshName;
   attribs["name"] = popName;
   if (vlsvWriter.writeArray("BLOCKIDS", attribs, totalBlocks, vectorSize, velocityBlockIds.data()) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write BLOCKIDS to file!" << endl << writeVerbose;
   {
      vector<vmesh::GlobalID>().swap(velocityBlockIds);
   }

   // Write the velocity space data
   // set everything that is needed for writing in data such as the array's name, size, data type, etc..
   attribs.clear();
   attribs["mesh"] = spatMeshName; // Name of the spatial mesh
   attribs["name"] = popName;      // Name of the velocity space distribution is written avgs
   const string datatype_avgs = "float";
   const uint64_t arraySize_avgs = totalBlocks;
   const uint64_t vectorSize_avgs = WID3; // There are 64 elements in every velocity block

   // Get the data size needed for writing in data
   uint64_t dataSize_avgs = sizeof(Realf);

   // Start multi write
   vlsvWriter.startMultiwrite(datatype_avgs,arraySize_avgs,vectorSize_avgs,dataSize_avgs);

   // Loop over cells
   for (size_t cell = 0; cell<cells.size(); ++cell) {
      // Get the spatial cell
      SpatialCell* SC = mpiGrid[cells[cell]];
      
      // Get the number of blocks in this cell
      const uint64_t arrayElements = SC->get_number_of_velocity_blocks();
      char* arrayToWrite = reinterpret_cast<char*>(SC->get_data());

      // Add a subarray to write
      vlsvWriter.addMultiwriteUnit(arrayToWrite, arrayElements); // Note: We told beforehands that the vectorsize = WID3 = 64
   }
   if (cells.size() == 0) {
      vlsvWriter.addMultiwriteUnit(NULL, 0); //Dummy write to avoid hang in end multiwrite
   }

   // Write the subarrays
   vlsvWriter.endMultiwrite("BLOCKVARIABLE", attribs);

   if (globalSuccess(success,"(MAIN) writeGrid: ERROR: Failed to fill temporary velocityBlockData array",MPI_COMM_WORLD) == false) {
      vlsvWriter.close();
      return false;
   }

   if (success ==false) {
      logFile << "(MAIN) writeGrid: ERROR occurred when writing BLOCKVARIABLE f" << endl << writeVerbose;
   }
   return success;
}

/*! Writes info received from data reducer. This function writes out the variable arrays into the file
 \param mpiGrid The Vlasiator's grid
 \param cells List of local cells (no ghost cells included)
 \param writeAsFloat If true, the data reducer writes variable arrays as float instead of double
 \param dataReducer The data reducer which contains the necessary functions for calculating variables
 \param dataReducerIndex Index in the data reducer (determines which variable to read) Note: size of the data reducer can be retrieved with dataReducer.size()
 \param vlsvWriter Some vlsv writer with a file open
 \return Returns true if operation was successful
 */
bool writeDataReducer(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                      const vector<uint64_t>& cells,
                      const bool writeAsFloat,
                      DataReducer& dataReducer,
                      int dataReducerIndex,
                      Writer& vlsvWriter){
   map<string,string> attribs;                      
   string variableName,dataType;
   bool success=true;

   //Get basic data on a variable:
   uint dataSize,vectorSize;
   attribs["mesh"] = "SpatialGrid";
   variableName = dataReducer.getName(dataReducerIndex);
   attribs["name"] = variableName;
   if (dataReducer.getDataVectorInfo(dataReducerIndex,dataType,dataSize,vectorSize) == false) {
      cerr << "ERROR when requesting info from DRO " << dataReducerIndex << endl;
      return false;
   }
  
   
   const uint64_t varBufferArraySize = cells.size()*vectorSize*dataSize;
   
   //Request DataReductionOperator to calculate the reduced data for all local cells:
   char* varBuffer = NULL;
   try {
      varBuffer = new char[varBufferArraySize];
   } catch( bad_alloc& ) {
      cerr << "ERROR, FAILED TO ALLOCATE MEMORY AT: " << __FILE__ << " " << __LINE__ << endl;
      logFile << "(MAIN) writeGrid: ERROR FAILED TO ALLOCATE MEMORY AT: " << __FILE__ << " " << __LINE__ << endl << writeVerbose;
      return false;
   }


   for (uint64_t cell=0; cell<cells.size(); ++cell) {
      //Reduce data ( return false if the operation fails )
      if (dataReducer.reduceData(mpiGrid[cells[cell]],dataReducerIndex,varBuffer + cell*vectorSize*dataSize) == false){
         success = false;
         logFile << "(MAIN) writeGrid: ERROR datareductionoperator '" << dataReducer.getName(dataReducerIndex) <<
            "' returned false!" << endl << writeVerbose;
      }
   }
   if( success ) {
      if( (writeAsFloat == true && dataType.compare("float") == 0) && dataSize == sizeof(double) ) {
         double * varBuffer_double = reinterpret_cast<double*>(varBuffer);
         //Declare smaller varbuffer:
         const uint64_t arraySize_smaller = cells.size();
         const uint32_t vectorSize_smaller = vectorSize;
         const uint32_t dataSize_smaller = sizeof(float);
         const string dataType_smaller = dataType;
         float * varBuffer_smaller = NULL;
         try {
            varBuffer_smaller = new float[arraySize_smaller * vectorSize_smaller];
         } catch( bad_alloc& ) {
            cerr << "ERROR, FAILED TO ALLOCATE MEMORY AT: " << __FILE__ << " " << __LINE__ << endl;
            logFile << "(MAIN) writeGrid: ERROR FAILED TO ALLOCATE MEMORY AT: " << __FILE__ << " " << __LINE__ << endl << writeVerbose;
            delete[] varBuffer;
            varBuffer = NULL;
            return false;
         }
         //Input varBuffer_double into varBuffer_smaller:
         for( uint64_t i = 0; i < arraySize_smaller * vectorSize_smaller; ++i ) {
            const double value = varBuffer_double[i];
            varBuffer_smaller[i] = (float)(value);
         }
         //Cast the varBuffer to char:
         char * varBuffer_smaller_char = reinterpret_cast<char*>(varBuffer_smaller);
         //Write the array:
         if (vlsvWriter.writeArray("VARIABLE", attribs, dataType_smaller, arraySize_smaller, vectorSize_smaller, dataSize_smaller, varBuffer_smaller_char) == false) {
            success = false;
            logFile << "(MAIN) writeGrid: ERROR failed to write datareductionoperator data to file!" << endl << writeVerbose;
         }
         delete[] varBuffer_smaller;
         varBuffer_smaller = NULL;
      } else {
         // Write  reduced data to file if DROP was successful:
         if (vlsvWriter.writeArray("VARIABLE",attribs, dataType, cells.size(), vectorSize, dataSize, varBuffer) == false) {
            success = false;
            logFile << "(MAIN) writeGrid: ERROR failed to write datareductionoperator data to file!" << endl << writeVerbose;
         }
      }
   }

   delete[] varBuffer;
   varBuffer = NULL;
   return success;
}




/*! Writes common grid data such as parameters (time steps, x_min, ..) as well as local cell ids as variables
 \param vlsvWriter Some vlsv writer with a file open
 \param mpiGrid Vlasiator's grid
 \param local_cells The local cell ids in this process
 \param fileIndex File index, file will be called "name.index.vlsv"
 \param comm The MPI comm
 \return Returns true if operation was successful
 */
bool writeCommonGridData(
   Writer& vlsvWriter,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const vector<uint64_t>& local_cells,
   const uint& fileIndex,
   MPI_Comm comm
) {
   // Writes parameters and cell ids into the VLSV file
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   const int masterProcessId = 0;
   //Write local cells into array as a variable:
   //Note: This needs to be done separately from the array MESH
   const short unsigned int vectorSize = 1;
   const uint32_t arraySize = local_cells.size();
   map<string, string> xmlAttributes;
   xmlAttributes["name"] = "CellID";
   xmlAttributes["mesh"] = "SpatialGrid";
   if( vlsvWriter.writeArray( "VARIABLE", xmlAttributes, arraySize, vectorSize, local_cells.data() ) == false ) {
      return false;
   }

   //Write parameters:
   if( vlsvWriter.writeParameter("time", &P::t) == false ) { return false; }
   if( vlsvWriter.writeParameter("dt", &P::dt) == false ) { return false; }
   if( vlsvWriter.writeParameter("timestep", &P::tstep) == false ) { return false; }
   if( vlsvWriter.writeParameter("fieldSolverSubcycles", &P::fieldSolverSubcycles) == false ) { return false; }
   if( vlsvWriter.writeParameter("fileIndex", &fileIndex) == false ) { return false; }
   if( vlsvWriter.writeParameter("xmin", &P::xmin) == false ) { return false; }
   if( vlsvWriter.writeParameter("xmax", &P::xmax) == false ) { return false; }
   if( vlsvWriter.writeParameter("ymin", &P::ymin) == false ) { return false; }
   if( vlsvWriter.writeParameter("ymax", &P::ymax) == false ) { return false; }
   if( vlsvWriter.writeParameter("zmin", &P::zmin) == false ) { return false; }
   if( vlsvWriter.writeParameter("zmax", &P::zmax) == false ) { return false; }
   if( vlsvWriter.writeParameter("xcells_ini", &P::xcells_ini) == false ) { return false; }
   if( vlsvWriter.writeParameter("ycells_ini", &P::ycells_ini) == false ) { return false; }
   if( vlsvWriter.writeParameter("zcells_ini", &P::zcells_ini) == false ) { return false; }
   if( vlsvWriter.writeParameter("vxmin", &P::vxmin) == false ) { return false; }
   if( vlsvWriter.writeParameter("vxmax", &P::vxmax) == false ) { return false; }
   if( vlsvWriter.writeParameter("vymin", &P::vymin) == false ) { return false; }
   if( vlsvWriter.writeParameter("vymax", &P::vymax) == false ) { return false; }
   if( vlsvWriter.writeParameter("vzmin", &P::vzmin) == false ) { return false; }
   if( vlsvWriter.writeParameter("vzmax", &P::vzmax) == false ) { return false; }
   if( vlsvWriter.writeParameter("vxblocks_ini", &P::vxblocks_ini) == false ) { return false; }
   if( vlsvWriter.writeParameter("vyblocks_ini", &P::vyblocks_ini) == false ) { return false; }
   if( vlsvWriter.writeParameter("vzblocks_ini", &P::vzblocks_ini) == false ) { return false; }
   if ( vlsvWriter.writeParameter("max_velocity_ref_level", &P::amrMaxVelocityRefLevel) == false) {return false;}

   //Mark the new version:
   float version = 1.00;
   if( vlsvWriter.writeParameter( "version", &version ) == false ) { return false; }
   return true; 
}


/*! Writes ghost cell ids into the file
 \param mpiGrid Vlasiator's grid
 \param vlsvWriter Some vlsv writer with a file open
 \param meshName Name of the mesh (If unsure, put SpatialGrid)
 \param ghost_cells List of cell ids on the process boundary (Ghost cells)
 \return Returns true if operation was successful
 \sa updateLocalIds
 */
bool writeGhostZoneDomainAndLocalIdNumbers( dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                              Writer & vlsvWriter,
                                              const string & meshName,
                                              const vector<uint64_t> & ghost_cells ) {
   //Declare vectors for storing data
   vector<uint64_t> ghostDomainIds;
   ghostDomainIds.reserve( ghost_cells.size() );
   vector<uint64_t> ghostLocalIds;
   ghostLocalIds.reserve( ghost_cells.size() );

   //Iterate through all ghost zones:
   vector<uint64_t>::const_iterator it;
   for( it = ghost_cells.begin(); it != ghost_cells.end(); ++it ) {
      //Domain id is the MPI process rank owning the ghost zone

      //get the local id of the zone in the process where THIS ghost zone is a local zone:
      //In order to do this we need MPI (Done in createZone)
      //Example:
      //Marking zones with a letter A, B, C, D etc and local ids are numbers above the zones
      //local id:                     0  1  2                  3 4 5
      //Process 1. has: local zones ( A, B, C ), ghost zones ( D E F )

      //local id:                     0  1                  2
      //Process 2. has: local zones ( D, G ), ghost zones ( H )
      //Now if we're in process 1. and our ghost zone is D, its domainId would be 2. because process 2. has D as local zone
      //In process 2, the local id of D is 0, so that's the local id we want now
      //The local id is being saved in createZone function     

      //Append to the vectors Note: Check updateLocalIds function
      ghostDomainIds.push_back( mpiGrid.get_process( *it ) );
      ghostLocalIds.push_back( mpiGrid[(*it)]->ioLocalCellId );
   }

   //We need the number of ghost zones for vlsvWriter:
   uint64_t numberOfGhosts = ghost_cells.size();

   //Write:
   map<string, string> xmlAttributes; //Used for writing in info
   //Note: should be "SpatialGrid"
   xmlAttributes["mesh"] = meshName;
   const unsigned int vectorSize = 1;
   //Write the in the number of ghost domains: (Returns false if writing fails)
   if( vlsvWriter.writeArray( "MESH_GHOST_DOMAINS", xmlAttributes, numberOfGhosts, vectorSize, ghostDomainIds.data() ) == false ) {
      cerr << "Error, failed to write MEST_GHOST_DOMAINS at: " << __FILE__ << " " << __LINE__ << endl;
      logFile << "(MAIN) writeGrid: ERROR failed to write MEST_GHOST_DOMAINS at: " << __FILE__ << " " << __LINE__ << endl << writeVerbose;
      return false;
   }
   //Write the in the number of ghost local ids: (Returns false if writing fails)
   if( vlsvWriter.writeArray( "MESH_GHOST_LOCALIDS", xmlAttributes, numberOfGhosts, vectorSize, ghostLocalIds.data()) == false ) {
      cerr << "Error, failed to write MEST_GHOST_LOCALIDS at: " << __FILE__ << " " << __LINE__ << endl;
      logFile << "(MAIN) writeGrid: ERROR failed to write MEST_GHOST_LOCALIDS at: " << __FILE__ << " " << __LINE__ << endl << writeVerbose;
      return false;
   }
   //Everything good
   return true;
}


/*! Writes domain sizes into the vlsv file, so the number of ghost and local cell ids in this process
 \param vlsvWriter Some vlsv writer with a file open
 \param meshName Name of the mesh (SpatialGrid used in the writeGrid function)
 \param numberOfLocalZones Number of local cells in this process
 \param numberOfGhostZones Number of ghost cells in this process ( Cells on the process boundary )
 \return Returns true if operation was successful
 */
bool writeDomainSizes( Writer & vlsvWriter,
                         const string & meshName,
                         const unsigned int & numberOfLocalZones,
                         const unsigned int & numberOfGhostZones ) {
   //Declare domainSize. There are two types of domain sizes -- ghost and local
   const unsigned int numberOfDomainTypes = 2;
   uint32_t domainSize[numberOfDomainTypes];
   domainSize[0] = numberOfLocalZones + numberOfGhostZones;
   domainSize[1] = numberOfGhostZones;

   //Write the array:
   map<string, string> xmlAttributes;
   //Put the meshName
   xmlAttributes["mesh"] = meshName;
   const unsigned int arraySize = 1;
   const unsigned int vectorSize = 2;
   //Write (writeArray does the writing) Note: Here the important part is "MESH_DOMAIN_SIZES" -- visit plugin needs this
   if( vlsvWriter.writeArray( "MESH_DOMAIN_SIZES", xmlAttributes, arraySize, vectorSize, domainSize ) == false ) {
      cerr << "Error at: " << __FILE__ << " " << __LINE__ << ", FAILED TO WRITE MESH_DOMAIN_SIZES" << endl;
      logFile << "(MAIN) writeGrid: ERROR FAILED TO WRITE MESH_DOMAIN_SIZES AT: " << __FILE__ << " " << __LINE__ << endl << writeVerbose;
      return false;
   }
   return true;
}



/*! Writes the zone global id numbers into the file. The vlsv file needs to know in which order the local cells + ghost cells are written. Local cells are first appended to a vector called global ids, after which the ghost cells are appended. The global ids vector will then be saved into a vlsv file
 \param mpiGrid Vlasiator's MPI grid
 \param vlsvWriter Some vlsv writer with a file open
 \param meshName Name of the mesh ("SpatialGrid" used in the writeGrid function and it should be the default)
 \param local_cells Vector containing the local cells of this process
 \param ghost_cells Vector containing the ghost cells of this process ( The cells on process boundary )
 \return Returns true if the operation was successful
 */
bool writeZoneGlobalIdNumbers( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                 Writer & vlsvWriter,
                                 const string & meshName,
                                 const vector<uint64_t> & local_cells,
                                 const vector<uint64_t> & ghost_cells ) {
   if( local_cells.empty() ) {
      if( !ghost_cells.empty() ) {
         //Something very wrong -- local zones should always have members when ghost zones has members
         cerr << "ERROR, LOCAL ZONES EMPTY BUT GHOST ZONES NOT AT " << __FILE__ << __LINE__ << endl;
         return false;
      }
   }

   //Get the cells in x, y, z direction right off the bat (for the sake of clarity):
   const unsigned int xCells = P::xcells_ini;
   const unsigned int yCells = P::ycells_ini;
   const unsigned int zCells = P::zcells_ini;

   vector<uint64_t> globalIds;
   globalIds.reserve( local_cells.size() + ghost_cells.size() );

   //Iterate through local_cells and store the values into globalIDs
   //Note: globalID is defined as follows: global ID = z*yCells*xCells + y*xCells + x
   vector<uint64_t>::const_iterator it;
   for( it = local_cells.begin(); it != local_cells.end(); ++it ) {
      if( (*it) == 0 ) {
         cerr << "ERROR, Invalid cell id at " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      //Add the global id:
      //Note: Unlike cell ids, global ids start from 0
      globalIds.push_back( (*it) - 1 );
   }
   //Do the same for ghost zones: (Append to the end of the list of global ids)
   for( it = ghost_cells.begin(); it != ghost_cells.end(); ++it ) {
      if( (*it) == 0 ) {
         cerr << "ERROR, Invalid cell id at " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      //Add the global id:
      globalIds.push_back( (*it) - 1 );
   }

   //Get the total number of zones:
   const uint64_t numberOfZones = globalIds.size();

   //Write the array:
   map<string, string> xmlAttributes;
   //The name of the mesh (user input -- should be "SpatialGrid")
   xmlAttributes["name"] = meshName;
   //A mandatory 'type' -- just something visit hopefully understands, because I dont (some of us do!) :)
   xmlAttributes["type"] = "multi_ucd";

   //Set periodicity:
   if( mpiGrid.topology.is_periodic( 0 ) ) { xmlAttributes["xperiodic"] = "yes"; } else { xmlAttributes["xperiodic"] = "no"; }
   if( mpiGrid.topology.is_periodic( 1 ) ) { xmlAttributes["yperiodic"] = "yes"; } else { xmlAttributes["yperiodic"] = "no"; }
   if( mpiGrid.topology.is_periodic( 2 ) ) { xmlAttributes["zperiodic"] = "yes"; } else { xmlAttributes["zperiodic"] = "no"; }
   //Write:
   if( numberOfZones == 0 ) {
      const uint64_t dummy_data = 0;
      const unsigned int dummy_array = 0;
      if( vlsvWriter.writeArray( "MESH", xmlAttributes, dummy_array, 1, &dummy_data ) == false ) {
         cerr << "Unsuccessful writing of MESH at: " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
   } else {
      if( vlsvWriter.writeArray( "MESH", xmlAttributes, numberOfZones, 1, globalIds.data() ) == false ) {
         cerr << "Unsuccessful writing of MESH at: " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
   }
   //Successfully wrote the array
   return true;
}

/*! Writes the node coordinates. This means basically every cell's node (Corner of each cell). Note: The grid is a structured one, so writing the nodes means starting from the corner of the grid and writing coordinates per every cell length until reaching the other corner of the grid
 \param vlsvWriter Some vlsv writer with a file open
 \param meshName Name of the mesh (SpatialGrid used in writeGrid and should be the default)
 \param masterRank The master rank (Vlasiator uses 0)
 \param comm The MPI comm
 \return Returns true if the operation was successful
 */
bool writeBoundingBoxNodeCoordinates ( Writer & vlsvWriter,
                                       const string & meshName,
                                       const int masterRank,
                                       MPI_Comm comm ) {

   //Create variables xCells, yCells, zCells which tell the number of zones in the given direction
   //Note: This is for the sake of clarity.
   const uint64_t & xCells = P::xcells_ini;
   const uint64_t & yCells = P::ycells_ini;
   const uint64_t & zCells = P::zcells_ini;

   //Create variables xmin, ymin, zmin for calculations
   const Real & xmin = (Real)P::xmin;
   const Real & ymin = (Real)P::ymin;
   const Real & zmin = (Real)P::zmin;

   //Create variables for cell lengths in x, y, z directions for calculations
   const Real & xCellLength = (Real)P::dx_ini;
   const Real & yCellLength = (Real)P::dy_ini;
   const Real & zCellLength = (Real)P::dz_ini;
   

   //Create node coordinates:
   //These are the coordinates for any given node in x y or z direction
   //Note: Nodes are basically the box coordinates
   vector<Real> xNodeCoordinates;
   xNodeCoordinates.reserve(xCells + 1);
   vector<Real> yNodeCoordinates;
   yNodeCoordinates.reserve(yCells + 1);
   vector<Real> zNodeCoordinates;
   zNodeCoordinates.reserve(zCells + 1);

   //Input the coordinates for the nodes:
   for( unsigned int i = 0; i < xCells + 1; ++i ) {
      //The x coordinate of the first node should be xmin, the second xmin + xCellLength and so on
      xNodeCoordinates.push_back(xmin + xCellLength * i);
   }
   for( unsigned int i = 0; i < yCells + 1; ++i ) {
      yNodeCoordinates.push_back(ymin + yCellLength * i);
   }
   for( unsigned int i = 0; i < zCells + 1; ++i ) {
      zNodeCoordinates.push_back(zmin + zCellLength * i);
   }

   //Write the arrays:
   map<string, string> xmlAttributes;
   //Note: meshName should be "SpatialGrid", probably
   xmlAttributes["mesh"] = meshName;
   //"success"'s value will be returned. By default it's true but if some of the vlsvWriter operations fail it will be false:
   bool success = true;
   //Depending on whether our rank is master rank or not the operation is slightly different, so let's get out rank from the MPI_Comm comm:
   int myRank;
   //Input myRank:
   MPI_Comm_rank(comm, &myRank);
   //Check the rank and write the arrays:
   const unsigned int vectorSize = 1;
   uint64_t arraySize;
   if( myRank == masterRank ) {
      //Save with the correct name "MESH_NODE_CRDS_X" -- writeArray returns false if something goes wrong
      arraySize = xCells + 1;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_X", xmlAttributes, arraySize, vectorSize, xNodeCoordinates.data()) == false ) success = false;
      arraySize = yCells + 1;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Y", xmlAttributes, arraySize, vectorSize, yNodeCoordinates.data()) == false ) success = false;
      arraySize = zCells + 1;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Z", xmlAttributes, arraySize, vectorSize, zNodeCoordinates.data()) == false ) success = false;
   } else {
      //Not a master process, so write empty:
      arraySize = 0;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_X", xmlAttributes, arraySize, vectorSize, xNodeCoordinates.data()) == false ) success = false;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Y", xmlAttributes, arraySize, vectorSize, yNodeCoordinates.data()) == false ) success = false;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Z", xmlAttributes, arraySize, vectorSize, zNodeCoordinates.data()) == false ) success = false;
   }
   //Free the memory
   xNodeCoordinates.clear();
   yNodeCoordinates.clear();
   zNodeCoordinates.clear();
   return success;
}


/*! Function for writing the bounding box. This writes only if the process running it is the master rank. This array contains info on the boundaries of the grid so for example the number of cells in x, y, z direction.
 \param vlsvWriter Some vlsv writer with a file open
 \param meshName Name of the mesh to write ("SpatialGrid" is used in writeGrid and it should be the default)
 \param masterRank The master process' id. Vlasiator uses 0 as the master process id, so by default should be 0
 \param comm MPI comm
 \return Returns true if operation was successful
 */
bool writeMeshBoundingBox( Writer & vlsvWriter, 
                           const string & meshName, 
                           const int masterRank,
                           MPI_Comm comm ) {
   //Get my rank from the MPI_Comm
   int myRank;
   MPI_Comm_rank(comm, &myRank);

   //Declare boundaryBox (writeArray expects it to tell the size of
   const unsigned int box_size = 6;
   const unsigned int notBlockBasedMesh = 1; // 1 because we are not interested in block based mesh
                                             //Note: If we were, the 3 last values in boundaryBox(below) would tell the
                                             //number of cells in blocks in x, y, z direction
   //Set the boundary box
   const uint64_t & numberOfXCells = P::xcells_ini;
   const uint64_t & numberOfYCells = P::ycells_ini;
   const uint64_t & numberOfZCells = P::zcells_ini;
   uint64_t boundaryBox[box_size] = { numberOfXCells, numberOfYCells, numberOfZCells, 
                                      notBlockBasedMesh, notBlockBasedMesh, notBlockBasedMesh };

   //Write:
   //Declare attributes
   map<string, string> xmlAttributes;
   //We received mesh name as a parameter: MOST LIKELY THIS IS SpatialGrid!
   xmlAttributes["mesh"] = meshName;

   //Write an array (NOTE: success will be returned and writeArray will return true or false depending on whether or not the write is successful)
   bool success;
   if( myRank == masterRank ) {
      //The visit plugin expects MESH_BBOX as a keyword
      //NOTE: writeArray writes boundaryBox
      const unsigned int arraySize = 6;
      const unsigned int vectorSize = 1;
      success = vlsvWriter.writeArray("MESH_BBOX", xmlAttributes, arraySize, vectorSize, boundaryBox);
   } else {
      const unsigned int arraySize = 0;
      const unsigned int vectorSize = 1;
      success = vlsvWriter.writeArray("MESH_BBOX", xmlAttributes, arraySize, vectorSize, boundaryBox);
   }
   return success;
}

/** This function writes the velocity space.
 * @param mpiGrid Vlasiator's grid.
 * @param vlsvWriter some vlsv writer with a file open.
 * @param index Index to call the correct member of the various parameter vectors.
 * @param cells Vector containing local cells of this process.
 * @return Returns true if the operation was successful.
 * @sa writeVelocityDistributionData. */
bool writeVelocitySpace(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                        Writer& vlsvWriter,int index,const vector<uint64_t>& cells) {
      //Compute which cells will write out their velocity space
      vector<uint64_t> velSpaceCells;
      int lineX, lineY, lineZ;
      for (uint i = 0; i < cells.size(); i++) {
         mpiGrid[cells[i]]->parameters[CellParams::ISCELLSAVINGF] = 0.0;
         // CellID stride selection
         if (P::systemWriteDistributionWriteStride[index] > 0 &&
             cells[i] % P::systemWriteDistributionWriteStride[index] == 0) {
            velSpaceCells.push_back(cells[i]);
            mpiGrid[cells[i]]->parameters[CellParams::ISCELLSAVINGF] = 1.0;
            continue; // Avoid double entries in case the cell also matches following conditions.
         }
         // Cell lines selection
         // Determine cellID's 3D indices
         lineX =  (cells[i]-1) % P::xcells_ini;
         lineY = ((cells[i]-1) / P::xcells_ini) % P::ycells_ini;
         lineZ = ((cells[i]-1) /(P::xcells_ini *  P::ycells_ini)) % P::zcells_ini;
         // Check that indices are in correct intersection at least in one plane
         if ((P::systemWriteDistributionWriteXlineStride[index] > 0 &&
              P::systemWriteDistributionWriteYlineStride[index] > 0 &&
              lineX % P::systemWriteDistributionWriteXlineStride[index] == 0 &&
              lineY % P::systemWriteDistributionWriteYlineStride[index] == 0)
             &&
             (P::systemWriteDistributionWriteYlineStride[index] > 0 &&
              P::systemWriteDistributionWriteZlineStride[index] > 0 &&
              lineY % P::systemWriteDistributionWriteYlineStride[index] == 0 &&
              lineZ % P::systemWriteDistributionWriteZlineStride[index] == 0)
             &&
             (P::systemWriteDistributionWriteZlineStride[index] > 0 &&
              P::systemWriteDistributionWriteXlineStride[index] > 0 &&
              lineZ % P::systemWriteDistributionWriteZlineStride[index] == 0 &&
              lineX % P::systemWriteDistributionWriteXlineStride[index] == 0)
         ) {
            velSpaceCells.push_back(cells[i]);
            mpiGrid[cells[i]]->parameters[CellParams::ISCELLSAVINGF] = 1.0;
         }
      }

      uint64_t numVelSpaceCells;
      uint64_t localNumVelSpaceCells;
      localNumVelSpaceCells=velSpaceCells.size();
      MPI_Allreduce(&localNumVelSpaceCells,&numVelSpaceCells,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
      //write out velocity space data NOTE: There is mpi communication in writeVelocityDistributionData
      if (writeVelocityDistributionData(vlsvWriter, mpiGrid, velSpaceCells, MPI_COMM_WORLD ) == false ) {
         cerr << "ERROR, FAILED TO WRITE VELOCITY DISTRIBUTION DATA AT " << __FILE__ << " " << __LINE__ << endl;
         logFile << "(MAIN) writeGrid: ERROR FAILED TO WRITE VELOCITY DISTRIBUTION DATA AT: " << __FILE__ << " " << __LINE__ << endl << writeVerbose;
      }
      return true;
}

/*! This function makes sure that local cells and ghost cells do not have any identical members (used for error checking)
 \param local_cells List of local cells within this process
 \param ghost_cells List of ghost cells within this process (cells on the process boundary)
 */
bool checkForSameMembers( const vector<uint64_t> local_cells, const vector<uint64_t> ghost_cells ) {
   //NOTE: VECTORS MUST BE SORTED
   //Make sure ghost cells and local cells don't have same members in them:
   vector<uint64_t>::const_iterator i = local_cells.begin();
   vector<uint64_t>::const_iterator j = ghost_cells.begin();
   while( i != local_cells.end() && j != ghost_cells.end() ) {
      if( (*i) < (*j) ) {
         ++i;
      } else if( (*i) > (*j) ) {
         ++j;
      } else {
         //Has a same member
         cerr << "ERROR SAME CELL ID " << *i << " -" << endl;
         logFile << "(MAIN) writeGrid: ERROR SAME CELL ID AT: " << __FILE__ << " " << __LINE__ << endl << writeVerbose;
         return true;
      }
   }
   return false;
}


/*!

\brief Write out system into a vlsv file

\param mpiGrid     The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute data that is added into file
\param index       Index to call the correct member of the various parameter vectors
\param writeGhosts If true, writes out ghost cells (cells that exist on the process boundary so other process' cells)
*/
bool writeGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
               DataReducer& dataReducer,
               const uint& index,
               const bool writeGhosts ) {
   double allStart = MPI_Wtime();
   bool success = true;
   int myRank;
   phiprof::initializeTimer("Barrier-entering-writegrid","MPI","Barrier");
   phiprof::start("Barrier-entering-writegrid");
   MPI_Barrier(MPI_COMM_WORLD);
   phiprof::stop("Barrier-entering-writegrid");


   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   phiprof::start("writeGrid-reduced");
   // Create a name for the output file and open it with VLSVWriter:
   stringstream fname;
   fname << P::systemWritePath.at(index) << "/" << P::systemWriteName.at(index) << ".";
   fname.width(7);
   fname.fill('0');
   fname << P::systemWrites.at(index) << ".vlsv";


   //Open the file with vlsvWriter:
   Writer vlsvWriter;
   const int masterProcessId = 0;
   MPI_Info MPIinfo = MPI_INFO_NULL;

   phiprof::start("open");
   vlsvWriter.open( fname.str(), MPI_COMM_WORLD, masterProcessId, MPIinfo );
   phiprof::stop("open");

   phiprof::start("metadataIO");

   // Get all local cell Ids 
   const vector<CellID>& local_cells = getLocalCells();
   
   //Declare ghost cells:
   vector<CellID> ghost_cells;
   if( writeGhosts ) {
      // Writing ghost cells:
      // Get all ghost cell Ids (NOTE: this works slightly differently depending on whether the grid is periodic or not)
      ghost_cells = mpiGrid.get_remote_cells_on_process_boundary( NEAREST_NEIGHBORHOOD_ID );
   }


   //Make sure the local cells and ghost cells are fetched properly
   if( local_cells.empty() ) {
      if( !ghost_cells.empty() ) {
         //Local cells empty but ghost cells not empty -- something very wrong
         cerr << "ERROR! LOCAL CELLS EMPTY BUT GHOST CELLS NOT AT: " << __FILE__ << " " << __LINE__ << endl;
      }
   }

   //The mesh name is "SpatialGrid" (This is used for writing in data)
   const string meshName = "SpatialGrid";

   //Write mesh boundaries: NOTE: master process only
   //Visit plugin needs to know the boundaries of the mesh so the number of cells in x, y, z direction
   if( writeMeshBoundingBox( vlsvWriter, meshName, masterProcessId, MPI_COMM_WORLD ) == false ) return false;

   //Write the node coordinates: NOTE: master process only
   if( writeBoundingBoxNodeCoordinates( vlsvWriter, meshName, masterProcessId, MPI_COMM_WORLD ) == false ) return false;

   //Write basic grid variables: NOTE: master process only
   if( writeCommonGridData(vlsvWriter, mpiGrid, local_cells, P::systemWrites[index], MPI_COMM_WORLD) == false ) return false;

   //Write zone global id numbers:
   if( writeZoneGlobalIdNumbers( mpiGrid, vlsvWriter, meshName, local_cells, ghost_cells ) == false ) return false;

   //Write domain sizes:
   if( writeDomainSizes( vlsvWriter, meshName, local_cells.size(), ghost_cells.size() ) == false ) return false;

   //Update local ids for cells:
   if( updateLocalIds( mpiGrid, local_cells, MPI_COMM_WORLD ) == false ) return false;

   //Write ghost zone domain and local id numbers ( VisIt plugin needs this for MPI )
   if( writeGhostZoneDomainAndLocalIdNumbers( mpiGrid, vlsvWriter, meshName, ghost_cells ) == false ) return false;
   phiprof::stop("metadataIO");
   phiprof::start("velocityspaceIO");
   if( writeVelocitySpace( mpiGrid, vlsvWriter, index, local_cells ) == false ) return false;
   phiprof::stop("velocityspaceIO");

   phiprof::start("reduceddataIO");
   //Write necessary variables:
   //Determines whether we write in floats or doubles
   for( uint i = 0; i < dataReducer.size(); ++i ) {
      if( writeDataReducer( mpiGrid, local_cells, (P::writeAsFloat==1), dataReducer, i, vlsvWriter ) == false ) return false;
   }
   phiprof::stop("reduceddataIO");

   phiprof::start("close");
   vlsvWriter.close();
   phiprof::stop("close");
   const uint64_t bytesWritten = vlsvWriter.getBytesWritten();
   phiprof::stop("writeGrid-reduced",bytesWritten*1e-9,"GB");
   return success;
}

/*!

\brief Write out a restart of the simulation into a vlsv file. All block data in remote cells will be reset.

\param mpiGrid   The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute data that is added into file
\param name       File name prefix, file will be called "name.index.vlsv"
\param fileIndex  File index, file will be called "name.index.vlsv"
*/
bool writeRestart(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  DataReducer& dataReducer,
                  const string& name,
                  const uint& fileIndex,
                  const int& stripe) {
   // Writes a restart
   double allStart = MPI_Wtime();
   bool success = true;
   int myRank;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   phiprof::initializeTimer("BarrierEnteringWriteRestart","MPI","Barrier");
   phiprof::start("BarrierEnteringWriteRestart");
   MPI_Barrier(MPI_COMM_WORLD);
   phiprof::stop("BarrierEnteringWriteRestart");

   phiprof::start("writeRestart");
   phiprof::start("DeallocateRemoteBlocks");
   //deallocate blocks in remote cells to decrease memory load
   deallocateRemoteCellBlocks(mpiGrid);
   phiprof::stop("DeallocateRemoteBlocks");
   
   // Get the current time.
   // Avoid different times on different processes!
   char currentDate[80];
   if(myRank == MASTER_RANK) {
      const time_t rawTime = time(NULL);
      const struct tm * timeInfo = localtime(&rawTime);
      strftime(currentDate, 80, "%F_%H-%M-%S", timeInfo);
   }
   MPI_Bcast(&currentDate,80,MPI_CHAR,MASTER_RANK,MPI_COMM_WORLD);
   
   // Create a name for the output file and open it with VLSVWriter:
   stringstream fname;
   fname << P::restartWritePath << "/" << name << ".";
   fname.width(7);
   fname.fill('0');
   fname << fileIndex << "." << currentDate << ".vlsv";

   phiprof::start("open");
   //Open the file with vlsvWriter:
   Writer vlsvWriter;
   const int masterProcessId = 0;
   MPI_Info MPIinfo; 
   if (stripe == 0 || stripe < -1){
      MPIinfo = MPI_INFO_NULL;
   } else {
      MPI_Info_create(&MPIinfo);
      char stripeChar[6];
      sprintf(stripeChar,"%d",stripe);
      /* no. of I/O devices to be used for file striping */
      char factor[] = "striping_factor";
      MPI_Info_set(MPIinfo, factor, stripeChar);
   }
   
   if( vlsvWriter.open( fname.str(), MPI_COMM_WORLD, masterProcessId, MPIinfo ) == false) return false;

   phiprof::stop("open");

   phiprof::start("metadataIO");
   
   // Get all local cell Ids 
   vector<CellID> local_cells = getLocalCells();
   //no order assumed so let's order cells here
   std::sort(local_cells.begin(), local_cells.end());
   
   //Note: No need to write ghost zones for write restart
   const vector<CellID> ghost_cells;
   
   //The mesh name is "SpatialGrid"
   const string meshName = "SpatialGrid";
   
   //Write mesh boundaries: NOTE: master process only
   //Visit plugin needs to know the boundaries of the mesh so the number of cells in x, y, z direction
   if( writeMeshBoundingBox( vlsvWriter, meshName, masterProcessId, MPI_COMM_WORLD ) == false ) return false;
   
   //Write the node coordinates: NOTE: master process only
   if( writeBoundingBoxNodeCoordinates( vlsvWriter, meshName, masterProcessId, MPI_COMM_WORLD ) == false ) return false;
   
   //Write basic grid parameters: NOTE: master process only ( I think )
   if( writeCommonGridData(vlsvWriter, mpiGrid, local_cells, fileIndex, MPI_COMM_WORLD) == false ) return false;
   
   //Write zone global id numbers:
   if( writeZoneGlobalIdNumbers( mpiGrid, vlsvWriter, meshName, local_cells, ghost_cells ) == false ) return false;
   phiprof::stop("metadataIO");
   phiprof::start("reduceddataIO");   
   //write out DROs we need for restarts
   DataReducer restartReducer;
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("background_B",CellParams::BGBX,3));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("perturbed_B",CellParams::PERBX,3));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("moments",CellParams::RHO,4));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("moments_dt2",CellParams::RHO_DT2,4));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("moments_r",CellParams::RHO_R,4));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("moments_v",CellParams::RHO_V,4));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("pressure",CellParams::P_11,3));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("pressure_dt2",CellParams::P_11_DT2,3));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("pressure_r",CellParams::P_11_R,3));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("pressure_v",CellParams::P_11_V,3));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("LB_weight",CellParams::LBWEIGHTCOUNTER,1));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("max_v_dt",CellParams::MAXVDT,1));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("max_r_dt",CellParams::MAXRDT,1));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("max_fields_dt",CellParams::MAXFDT,1));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("rho_loss_adjust",CellParams::RHOLOSSADJUST,1));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("rho_loss_velocity_boundary",CellParams::RHOLOSSVELBOUNDARY,1));
   
   restartReducer.addOperator(new DRO::MPIrank);
   restartReducer.addOperator(new DRO::BoundaryType);
   restartReducer.addOperator(new DRO::BoundaryLayer);
   
   //Write necessary variables:
   const bool writeAsFloat = false;
   for (uint i=0; i<restartReducer.size(); ++i) {
      writeDataReducer(mpiGrid, local_cells, writeAsFloat, restartReducer, i, vlsvWriter);
   }
   phiprof::stop("reduceddataIO");   
   //write the velocity distribution data -- note: it's expecting a vector of pointers:
   // Note: restart should always write double values to ensure the accuracy of the restart runs. 
   // In case of distribution data it is not as important as they are mainly used for visualization purpose
   phiprof::start("velocityspaceIO");
   writeVelocityDistributionData(vlsvWriter, mpiGrid, local_cells, MPI_COMM_WORLD);
   phiprof::stop("velocityspaceIO");

   phiprof::start("close");
   vlsvWriter.close();
   phiprof::stop("close");

   phiprof::start("updateRemoteBlocks");
   //Updated newly adjusted velocity block lists on remote cells, and
   //prepare to receive block data
   updateRemoteVelocityBlockLists(mpiGrid);
   phiprof::stop("updateRemoteBlocks");

   const uint64_t bytesWritten = vlsvWriter.getBytesWritten();
   phiprof::stop("writeRestart",bytesWritten*1e-9,"GB");
   
   return success;
}


/*!

\brief Write out simulation diagnostics into diagnostic.txt

\param mpiGrid   The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute diagnostic data
*/
bool writeDiagnostic(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                     DataReducer& dataReducer)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   
   string dataType;
   uint dataSize, vectorSize;
   const vector<CellID>& cells = getLocalCells();
   cuint nCells = cells.size();
   cuint nOps = dataReducer.size();
   
   // Exit if the user does not want any diagnostics output
   if (nOps == 0) return true;

   vector<Real> localMin(nOps), localMax(nOps), localSum(nOps+1), localAvg(nOps),
               globalMin(nOps),globalMax(nOps),globalSum(nOps+1),globalAvg(nOps);
   localSum[0] = 1.0 * nCells;
   Real buffer;
   bool success = true;
   static bool printDiagnosticHeader = true;
   
   if (printDiagnosticHeader == true && myRank == MASTER_RANK) {
      if (P::isRestart){
         diagnostic << "# ==== Restart from file "<< P::restartFileName << " ===="<<endl;
      }
      diagnostic << "# Column 1 Step" << endl;
      diagnostic << "# Column 2 Simulation time" << endl;
      diagnostic << "# Column 3 Time step dt" << endl;
      for (uint i=0; i<nOps; ++i) {
         diagnostic << "# Columns " << 4 + i*4 << " to " << 7 + i*4 << ": " << dataReducer.getName(i) << " min max sum average" << endl;
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
      
      if (success == false) logFile << "(MAIN) writeDiagnostic: ERROR datareductionoperator '" << dataReducer.getName(i) <<
                               "' returned false!" << endl << writeVerbose;
   }
   
   MPI_Reduce(&localMin[0], &globalMin[0], nOps, MPI_Type<Real>(), MPI_MIN, 0, MPI_COMM_WORLD);
   MPI_Reduce(&localMax[0], &globalMax[0], nOps, MPI_Type<Real>(), MPI_MAX, 0, MPI_COMM_WORLD);
   MPI_Reduce(&localSum[0], &globalSum[0], nOps + 1, MPI_Type<Real>(), MPI_SUM, 0, MPI_COMM_WORLD);
   
   diagnostic << setprecision(12); 
   diagnostic << Parameters::tstep << "\t";
   diagnostic << Parameters::t << "\t";
   diagnostic << Parameters::dt << "\t";
   
   for (uint i=0; i<nOps; ++i) {
      if (globalSum[0] != 0.0) globalAvg[i] = globalSum[i+1] / globalSum[0];
      else globalAvg[i] = globalSum[i+1];
      if (myRank == MASTER_RANK) {
         diagnostic << globalMin[i] << "\t" <<
         globalMax[i] << "\t" <<
         globalSum[i+1] << "\t" <<
         globalAvg[i] << "\t";
      }
   }
   if (myRank == MASTER_RANK) diagnostic << endl << write;
   return true;
}



