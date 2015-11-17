#include <cstdlib>
#include <iostream>
#include <iomanip> // for setprecision()
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>
#include "ioread.h"
#include "phiprof.hpp"
#include "parameters.h"
#include "logger.h"
#include "vlsv_reader_parallel.h"
#include "vlasovmover.h"

using namespace std;
using namespace phiprof;
using namespace vlsv;

extern Logger logFile, diagnostic;

typedef Parameters P;

/*!
 * \brief Checks for command files written to the local directory.
 * If a file STOP was written and is readable, then a bailout with restart writing is initiated.
 * If a file KILL was written and is readable, then a bailout without a restart is initiated.
 * If a file SAVE was written and is readable, then restart writing without a bailout is initiated.
 * To avoid bailing out upfront on a new run the files are renamed with the date to keep a trace.
 * The function should only be called by MASTER_RANK. This ensures that resetting P::bailout_write_restart works.
 */
void checkExternalCommands() {
   struct stat tempStat;
   if (stat("STOP", &tempStat) == 0) {
      bailout(true, "Received an external STOP command. Setting bailout.write_restart to true.");
      P::bailout_write_restart = true;
      char newName[80];
      // Get the current time.
      const time_t rawTime = time(NULL);
      const struct tm * timeInfo = localtime(&rawTime);
      strftime(newName, 80, "STOP_%F_%H-%M-%S", timeInfo);
      rename("STOP", newName);
      return;
   }
   if(stat("KILL", &tempStat) == 0) {
      bailout(true, "Received an external KILL command. Setting bailout.write_restart to false.");
      P::bailout_write_restart = false;
      char newName[80];
      // Get the current time.
      const time_t rawTime = time(NULL);
      const struct tm * timeInfo = localtime(&rawTime);
      strftime(newName, 80, "KILL_%F_%H-%M-%S", timeInfo);
      rename("KILL", newName);
      return;
   }
   if(stat("SAVE", &tempStat) == 0) {
      cerr << "Received an external SAVE command. Writing a restart file." << endl;
      globalflags::writeRestart = true;
      char newName[80];
      // Get the current time.
      const time_t rawTime = time(NULL);
      const struct tm * timeInfo = localtime(&rawTime);
      strftime(newName, 80, "SAVE_%F_%H-%M-%S", timeInfo);
      rename("SAVE", newName);
      return;
   }
}

/*!
  \brief Collective exit on error functions

  If any process(es) have a false success values then the program will
  abort and write out the message to the logfile
  \param success This process' value for successful execution
  \param message Error message to be shown upon failure
  \param comm MPI comm
  \return Returns true if the operation was successful
*/
bool exitOnError(bool success,string message,MPI_Comm comm) {
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
      logFile << message << endl<<write ;
      exit(1);
   }
}


/*!
 \brief Read cell ID's
 Read in cell ID's from file.
 \param file Vlsv reader with a file open
 \param fileCells Vector in whic to store the cell ids
 \param masterRank The simulation's master rank id (Vlasiator uses 0, which should be the default)
 \param comm MPI comm (MPI_COMM_WORLD should be the default)
*/
bool readCellIds(ParallelReader & file,
                 vector<uint64_t>& fileCells, const int masterRank,MPI_Comm comm){
   // Get info on array containing cell Ids:
   uint64_t arraySize = 0;
   uint64_t vectorSize;
   datatype::type dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   bool success=true;
   int rank;
   MPI_Comm_rank(comm,&rank);
   if(rank==masterRank){

      const short int readFromFirstIndex = 0;
      //let's let master read cellId's, we anyway have at max ~1e6 cells
      attribs.push_back(make_pair("name","CellID"));
      attribs.push_back(make_pair("mesh","SpatialGrid"));
      if (file.getArrayInfoMaster("VARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
         logFile << "(RESTART) ERROR: Failed to read cell ID array info!" << endl << write;
         return false;
      }

      //Make a routine error check:
      if( vectorSize != 1 ) {
         logFile << "(RESTART) ERROR: Bad vectorsize at " << __FILE__ << " " << __LINE__ << endl << write;
         return false;
      }
      
      //   Read cell Ids:
      char* IDbuffer = new char[arraySize*vectorSize*byteSize];
      if (file.readArrayMaster("VARIABLE",attribs,readFromFirstIndex,arraySize,IDbuffer) == false) {
         logFile << "(RESTART) ERROR: Failed to read cell Ids!" << endl << write;
         success = false;
      }
   
   // Convert global Ids into our local DCCRG 64 bit uints
      const uint64_t & numberOfCells = arraySize;
      fileCells.resize(numberOfCells);
      if (dataType == datatype::type::UINT && byteSize == 4) {
         uint32_t* ptr = reinterpret_cast<uint32_t*>(IDbuffer);
         //Input cell ids
         for (uint64_t i=0; i<numberOfCells; ++i) {
            const uint64_t cellID = ptr[i];
            fileCells[i] = cellID;
         }
      } else if (dataType == datatype::type::UINT && byteSize == 8) {
         uint64_t* ptr = reinterpret_cast<uint64_t*>(IDbuffer);
         for (uint64_t i=0; i<numberOfCells; ++i) {
            const uint64_t cellID = ptr[i];
            fileCells[i] = cellID;
         }
      } else {
         logFile << "(RESTART) ERROR: ParallelReader returned an unsupported datatype for cell Ids!" << endl << write;
         success = false;
      }
      delete[] IDbuffer;
   }
   //broadcast cellId's to everybody
   MPI_Bcast(&arraySize,1,MPI_UINT64_T,masterRank,comm);   
   fileCells.resize(arraySize);
   MPI_Bcast(&(fileCells[0]),arraySize,MPI_UINT64_T,masterRank,comm);
   
   return success;
}



/*!
 \brief Read number of blocks per cell
 \param file Vlsv reader with a file open 
 \param nBlocks Vector for holding information on cells and the number of blocks in them -- this function saves data here
 \param masterRank The master rank of this process (Vlasiator uses masterRank = 0 and so it should be the default)
 \param comm MPI comm
 \return Returns true if the operation was successful
 \sa readGrid
*/
bool readNBlocks( ParallelReader & file,
                 vector<unsigned int>& nBlocks, int masterRank,MPI_Comm comm){
   // Get info on array containing cell Ids:
   uint64_t arraySize;
   uint64_t vectorSize;
   datatype::type dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   bool success=true;
   int rank;
   MPI_Comm_rank(comm,&rank);
   if(rank==masterRank){
      //master reads data
      attribs.push_back(make_pair("mesh","SpatialGrid"));
      if (file.getArrayInfoMaster("BLOCKSPERCELL",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
         logFile << "(RESTART) ERROR: Failed to read number of blocks" << endl << write;
         success= false;
      }

      const short int readFileFromBeginning = 0;
      nBlocks.resize(vectorSize*arraySize);
      if (file.readArrayMaster("BLOCKSPERCELL",attribs,readFileFromBeginning,arraySize,(char*)&(nBlocks[0])) == false) {
         logFile << "(RESTART) ERROR: Failed to read number of blocks!" << endl << write;
         success = false;
      }
   }

   //now broadcast the data to everybody
   MPI_Bcast(&arraySize,1,MPI_UINT64_T,masterRank,comm);
   MPI_Bcast(&vectorSize,1,MPI_UINT64_T,masterRank,comm);
   nBlocks.resize(vectorSize*arraySize);
   MPI_Bcast(&(nBlocks[0]),vectorSize*arraySize,MPI_UNSIGNED,masterRank,comm);
   return success;
}

//Outputs the velocity block indices of some given block into indices
//Input:
//[0] cellStruct -- some cell structure that has been constructed properly
//[1] block -- some velocity block id
//Output:
//[0] indices -- the array where to store the indices

/*! Outputs given block's velocity min coordinates (the corner of the block) itn blockCoordinates
 \param block The block's id
 \param blockCoordinates An empty array where to store the block coordinates
 \sa readBlockData
 */
void getVelocityBlockCoordinates( const uint64_t & block, boost::array<Real, 3> & blockCoordinates ) {
   //Get indices:
   boost::array<uint64_t, 3> blockIndices;
   blockIndices[0] = block % P::vxblocks_ini;
   blockIndices[1] = (block / P::vxblocks_ini) % P::vyblocks_ini;
   blockIndices[2] = block / (P::vxblocks_ini * P::vyblocks_ini);
   //Store the coordinates:
   blockCoordinates[0] = P::vxmin + ((P::vxmax - P::vxmin) / P::vxblocks_ini) * blockIndices[0];
   blockCoordinates[1] = P::vymin + ((P::vymax - P::vymin) / P::vyblocks_ini) * blockIndices[1];
   blockCoordinates[2] = P::vzmin + ((P::vzmax - P::vzmin) / P::vzblocks_ini) * blockIndices[2];
   return;
}


/*!
  This reads in velocity space data
  \param file                  Vlsv reader with a file open
  \param fileCells             List of all cell ids
  \param localCellStartOffset  The offset from which to start reading cells ( This should be balanced so that every process has roughly the same amount of blocks to read )
  \param localCells            How many cells after the offset to read ( This should be balanced so that every process has roughly the same amount of blocks to read )
  \param localBlockStartOffset The offset from which to start reading block data for thisp process ( Calculated from nBlocks and localCellStartOffset )
  \param localBlocks           Total number of blocks forthis process ( Calculated from nBlocks and localCellStartOffset and localCellStartOffset )
  \param mpiGrid               Vlasiator's grid
 \sa readGrid
*/

   
bool readBlockData(
   ParallelReader & file,
   const vector<uint64_t>& fileCells,
   const uint64_t localCellStartOffset,
   const uint64_t localCells,
   const vector<uint>& nBlocks,
   const uint64_t localBlockStartOffset,
   const uint64_t localBlocks,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
) {
  uint64_t arraySize;
  uint64_t avgVectorSize;
  uint64_t cellParamsVectorSize;
  datatype::type dataType;
  uint64_t byteSize;
  list<pair<string,string> > avgAttribs;
  bool success=true;
   
  avgAttribs.push_back(make_pair("name","avgs"));
  avgAttribs.push_back(make_pair("mesh","SpatialGrid"));
  
  //Get block id array info and store them into blockIdAttribs, lockIdByteSize, blockIdDataType, blockIdVectorSize
  list<pair<string,string> > blockIdAttribs;
  uint64_t blockIdVectorSize, blockIdByteSize;
  datatype::type blockIdDataType;
  blockIdAttribs.push_back( make_pair("mesh", "SpatialGrid") );
  if (file.getArrayInfo("BLOCKIDS",blockIdAttribs,arraySize,blockIdVectorSize,blockIdDataType,blockIdByteSize) == false ){
     logFile << "(RESTART) ERROR: Failed to read BLOCKIDS array info " << endl << write;
     return false;
  }
  if(file.getArrayInfo("BLOCKVARIABLE", avgAttribs, arraySize, avgVectorSize, dataType, byteSize) == false ){
    logFile << "(RESTART) ERROR: Failed to read BLOCKVARIABLE array info " << endl << write;
    return false;
  }

   //Some routine error checks:
   if( avgVectorSize!=WID3 ){
      logFile << "(RESTART) ERROR: Blocksize does not match in restart file " << endl << write;
      return false;
   }
   if( sizeof(Realf) != byteSize ) {
      logFile << "(RESTART) ERROR: Bad avgs bytesize at " << __FILE__ << " " << __LINE__ << endl << write;
      return false;
   }
   
   if( blockIdByteSize  != sizeof(vmesh::GlobalID)) {
      logFile << "(RESTART) ERROR: BlockID data size does not match " << __FILE__ << " " << __LINE__ << endl << write;
      return false;
   }
   
   vmesh::GlobalID * blockIdBuffer = new vmesh::GlobalID[blockIdVectorSize * localBlocks]; //blockids of all cells
   
   //Read block ids and data
   file.readArray("BLOCKIDS", blockIdAttribs, localBlockStartOffset, localBlocks, (char*)blockIdBuffer );
   file.startMultiread("BLOCKVARIABLE", avgAttribs);
   uint64_t blockBufferOffset=0;
   //Go through all spatial cells     
   std::vector<vmesh::GlobalID> blockIdsInCell; //blockIds in a particular cell, temporary usage
   for(uint i=0;i<localCells;i++){
      CellID cell = fileCells[localCellStartOffset + i]; //spatial cell id 
      uint nBlocksInCell = nBlocks[localCellStartOffset + i];
      //copy blocks in this cell to vector blockIdsInCell, size of read in data has been checked earlier
      blockIdsInCell.assign(blockIdBuffer + blockBufferOffset,
                            blockIdBuffer + blockBufferOffset + nBlocksInCell);
      mpiGrid[cell]->add_velocity_blocks(blockIdsInCell); //allocate space for all blocks and create them
      //register read
      file.addMultireadUnit((char*)(mpiGrid[cell]->get_data()), nBlocksInCell);
      blockBufferOffset += nBlocksInCell; //jump to location of next local cell
   }
   
   file.endMultiread(localBlockStartOffset);

   delete[] blockIdBuffer;
   return success;
}




/*! Reads cell parameters from the file and saves them in the right place in mpiGrid
 \param file Some parallel vlsv reader with a file open
 \param fileCells List of all cell ids
 \param localCellStartOffset Offset in the fileCells list for this process ( calculated so that the amount of blocks is distributed somewhat evenly between processes)
 \param localCells The amount of cells to read in this process after localCellStartOffset
 \param cellParamsIndex The parameter of the cell index e.g. CellParams::RHO
 \param expectedVectorSize The amount of elements in the parameter (parameter can be a scalar or a vector of size N)
 \param mpiGrid Vlasiator's grid (the parameters are saved here)
 \return Returns true if the operation is successful
 */
template <typename fileReal>
static bool _readCellParamsVariable(
   ParallelReader & file,
   const vector<uint64_t>& fileCells,
   const uint64_t localCellStartOffset,
   const uint64_t localCells,
   const string& variableName,
   const size_t cellParamsIndex,
   const size_t expectedVectorSize,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
) {
   uint64_t arraySize;
   uint64_t vectorSize;
   datatype::type dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   fileReal *buffer;
   bool success=true;
   
   attribs.push_back(make_pair("name",variableName));
   attribs.push_back(make_pair("mesh","SpatialGrid"));
   
   
   if (file.getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      logFile << "(RESTART)  ERROR: Failed to read " << endl << write;
      return false;
   }

   if(vectorSize!=expectedVectorSize){
      logFile << "(RESTART)  vectorsize wrong " << endl << write;
      return false;
   }
   
   buffer=new fileReal[vectorSize*localCells];
   if(file.readArray("VARIABLE",attribs,localCellStartOffset,localCells,(char *)buffer) == false ) {
      logFile << "(RESTART)  ERROR: Failed to read " << variableName << endl << write;
      return false;
   }
   
   for(uint i=0;i<localCells;i++){
     uint cell=fileCells[localCellStartOffset+i];
     for(uint j=0;j<vectorSize;j++){
        mpiGrid[cell]->parameters[cellParamsIndex+j]=buffer[i*vectorSize+j];
     }
   }
   
   delete[] buffer;
   return success;
}

/*! Reads cell parameters from the file and saves them in the right place in mpiGrid
 \param file Some parallel vlsv reader with a file open
 \param fileCells List of all cell ids
 \param localCellStartOffset Offset in the fileCells list for this process ( calculated so that the amount of blocks is distributed somewhat evenly between processes)
 \param localCells The amount of cells to read in this process after localCellStartOffset
 \param cellParamsIndex The parameter of the cell index e.g. CellParams::RHO
 \param expectedVectorSize The amount of elements in the parameter (parameter can be a scalar or a vector of size N)
 \param mpiGrid Vlasiator's grid (the parameters are saved here)
 \return Returns true if the operation is successful
 */
bool readCellParamsVariable(
   ParallelReader & file,
   const vector<uint64_t>& fileCells,
   const uint64_t localCellStartOffset,
   const uint64_t localCells,
   const string& variableName,
   const size_t cellParamsIndex,
   const size_t expectedVectorSize,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
) {
   uint64_t arraySize;
   uint64_t vectorSize;
   datatype::type dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   
   attribs.push_back(make_pair("name",variableName));
   attribs.push_back(make_pair("mesh","SpatialGrid"));
   
   
   if (file.getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      logFile << "(RESTART)  ERROR: Failed to read " << endl << write;
      return false;
   }

   // Call _readCellParamsVariable
   if( dataType == datatype::type::FLOAT ) {
      switch (byteSize) {
         case sizeof(double):
            return _readCellParamsVariable<double>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
         case sizeof(float):
            return _readCellParamsVariable<float>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
      }
   } else if( dataType == datatype::type::UINT ) {
      switch (byteSize) {

         case sizeof(uint32_t):
            return _readCellParamsVariable<uint32_t>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
         case sizeof(uint64_t):
            return _readCellParamsVariable<uint64_t>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
      }
   } else if( dataType == datatype::type::INT ) {
      switch (byteSize) {
         case sizeof(int32_t):
            return _readCellParamsVariable<int32_t>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
         case sizeof(int64_t):
            return _readCellParamsVariable<int64_t>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
      }
   } else {
      logFile << "(RESTART)  ERROR: Failed to read data type at readCellParamsVariable" << endl << write;
      return false;
   }
}


/*! A function for reading parameters e.g. 'timestep'
 \param file Vlsv parallel reader with a file open
 \param name Name of the parameter
 \param value Variable in which to store the scalar variable (double, float, int .. )
 \param masterRank The master process' id (Vlasiator uses 0 so this should equal 0 by default)
 \param comm MPI comm (MPI_COMM_WORLD should be the default)
 \return Returns true if the operation is successful
 */
template <typename T>
bool readScalarParameter(ParallelReader & file, const string & name,T& value, int masterRank,MPI_Comm comm){
   return file.readParameter( name, value );
}

/*! A function for checking the scalar parameter
 \param file Some parallel vlsv reader with a file open
 \param name Name of the parameter
 \param correctValue The correct value of the parameter to compare to
 \param masterRank The master process' id (Vlasiator uses 0 so this should be 0 by default)
 \param comm MPI comm (Default should be MPI_COMM_WORLD)
 \return Returns true if the operation is successful
 */
template <typename T>
bool checkScalarParameter(ParallelReader & file, const string & name, T correctValue, int masterRank,MPI_Comm comm){
   T value;
   readScalarParameter(file,name,value,masterRank,comm);
   if(value!=correctValue){
      std::ostringstream s;
      s << "(RESTART) Parameter " << name << " has mismatching value.";
      s << " CFG value = " << correctValue;
      s << " Restart file value = " << value;
      exitOnError(false,s.str(),MPI_COMM_WORLD);
      return false;
   }
   else{
      exitOnError(true,"",MPI_COMM_WORLD);
      return true;
   }
}



/*!
\brief Read in state from a vlsv file in order to restart simulations
\param mpiGrid Vlasiator's grid
\param name Name of the restart file e.g. "restart.00052.vlsv"
 \return Returns true if the operation was successful
*/
bool readGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
              const std::string& name){

   vector<uint64_t> fileCells; /*< CellIds for all cells in file*/
   vector<uint> nBlocks;/*< Number of blocks for all cells in file*/
   bool success=true;
   int myRank,processes;
   
   
   // Attempt to open VLSV file for reading:
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   MPI_Comm_size(MPI_COMM_WORLD,&processes);
   
   
   phiprof::start("readGrid");

   ParallelReader file;
   MPI_Info mpiInfo = MPI_INFO_NULL;

   if (file.open(name,MPI_COMM_WORLD,MASTER_RANK,mpiInfo) == false) {
      success=false;
   }
   exitOnError(success,"(RESTART) Could not open file",MPI_COMM_WORLD);

   //Around May 2015 time was renamed from "t" to "time", we try to read both, new way is read first
   if(readScalarParameter(file,"time",P::t,MASTER_RANK,MPI_COMM_WORLD) == false)
      if(readScalarParameter(file,"t",P::t,MASTER_RANK,MPI_COMM_WORLD) == false)
         success=false;
   P::t_min=P::t;

   //Around May 2015 timestep was renamed from "tstep" to "timestep", we try to read both, new way is read first   
   if(readScalarParameter(file,"timestep",P::tstep,MASTER_RANK,MPI_COMM_WORLD) == false)
      if(readScalarParameter(file,"tstep",P::tstep,MASTER_RANK,MPI_COMM_WORLD) == false)
         success=false;
   P::tstep_min=P::tstep;

   if(readScalarParameter(file,"dt",P::dt,MASTER_RANK,MPI_COMM_WORLD) ==false) success=false;

   if(readScalarParameter(file,"fieldSolverSubcycles",P::fieldSolverSubcycles,MASTER_RANK,MPI_COMM_WORLD) ==false) {
      // Legacy restarts do not have this field, it "should" be safe for one or two steps...
      P::fieldSolverSubcycles = 1.0;
      std::cout << " No P::fieldSolverSubcycles found in restart, setting 1." << std::endl;
   }
   MPI_Bcast(&(P::fieldSolverSubcycles),1,MPI_Type<Real>(),MASTER_RANK,MPI_COMM_WORLD);
   
   checkScalarParameter(file,"xmin",P::xmin,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"ymin",P::ymin,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"zmin",P::zmin,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"xmax",P::xmax,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"ymax",P::ymax,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"zmax",P::zmax,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"xcells_ini",P::xcells_ini,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"ycells_ini",P::ycells_ini,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"zcells_ini",P::zcells_ini,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"vxmin",P::vxmin,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"vymin",P::vymin,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"vzmin",P::vzmin,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"vxmax",P::vxmax,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"vymax",P::vymax,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"vzmax",P::vzmax,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"vxblocks_ini",P::vxblocks_ini,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"vyblocks_ini",P::vyblocks_ini,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"vzblocks_ini",P::vzblocks_ini,MASTER_RANK,MPI_COMM_WORLD);

   phiprof::start("readDatalayout");
   if(success) { success=readCellIds(file,fileCells,MASTER_RANK,MPI_COMM_WORLD); }

   //check that the cellID lists are identical in file and grid
   if(myRank==0){
      vector<CellID> allGridCells=mpiGrid.get_all_cells();
      if(fileCells.size() != allGridCells.size()){
         success=false;
      }
   }
   
   exitOnError(success,"(RESTART) Wrong number of cells in restartfile",MPI_COMM_WORLD);
   if(success) {
      success = readNBlocks(file,nBlocks,MASTER_RANK,MPI_COMM_WORLD);
   }
   //make sure all cells are empty, we will anyway overwrite everything and in that case moving cells is easier...
   {
      const vector<CellID>& gridCells = getLocalCells();
      for(uint i=0;i<gridCells.size();i++){
         mpiGrid[gridCells[i]]->clear();
      }
   }
   
   uint64_t totalNumberOfBlocks=0;
   unsigned int numberOfBlocksPerProcess;

   for(uint i=0;i<nBlocks.size();i++){
      totalNumberOfBlocks+=nBlocks[i];
   }
   numberOfBlocksPerProcess=1+totalNumberOfBlocks/processes;



   uint64_t localCellStartOffset=0; /*<!this is where local cells start in file-list after migration*/
   uint64_t localCells=0;
   
   uint64_t numberOfBlocksCount=0;
   //pin local cells to remote processes, we try to balance number of blocks so that each process has the same amount of blocks, more or less

   for(uint i=0;i<fileCells.size();i++){
      numberOfBlocksCount+=nBlocks[i];
      int newCellProcess=numberOfBlocksCount/numberOfBlocksPerProcess;
      if(newCellProcess==myRank) {
         if(localCells==0)
            localCellStartOffset=i; //here local cells start
         localCells++;
      }
      if(mpiGrid.is_local(fileCells[i])){
          mpiGrid.pin(fileCells[i],newCellProcess);
      }
   }
   
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_DATA);
   mpiGrid.balance_load(false);
   //update list of local gridcells
   recalculateLocalCellsCache();
   //get new list of local gridcells
   const vector<CellID>& gridCells = getLocalCells();
   //unpin cells, otherwise we will never change this initial bad balance
   for(uint i=0;i<gridCells.size();i++){
      mpiGrid.unpin(gridCells[i]);
   }

   //check for errors, has migration succeeded
   if(localCells != gridCells.size() ){
      success=false;
   }
   if(success) {
      for(uint i=localCellStartOffset;i< localCellStartOffset+localCells;i++){
         if(!mpiGrid.is_local(fileCells[i])) {
            success=false;
         }
      }
   }

   exitOnError(success,"(RESTART) Cell migration failed",MPI_COMM_WORLD);

   //set cell coordinates based on cfg (mpigrid) information
   for(uint i=0;i<gridCells.size();i++){
      std::array<double, 3> cell_min = mpiGrid.geometry.get_min(gridCells[i]);
      std::array<double, 3> cell_length = mpiGrid.geometry.get_length(gridCells[i]);
      
      mpiGrid[gridCells[i]]->parameters[CellParams::XCRD] = cell_min[0];
      mpiGrid[gridCells[i]]->parameters[CellParams::YCRD] = cell_min[1];
      mpiGrid[gridCells[i]]->parameters[CellParams::ZCRD] = cell_min[2];
      mpiGrid[gridCells[i]]->parameters[CellParams::DX  ] = cell_length[0];
      mpiGrid[gridCells[i]]->parameters[CellParams::DY  ] = cell_length[1];
      mpiGrid[gridCells[i]]->parameters[CellParams::DZ  ] = cell_length[2];
   }
   //where local data start in the blocklists
   uint64_t localBlockStartOffset=0;
   for(uint i=0;i<localCellStartOffset;i++){
      localBlockStartOffset+=nBlocks[i];
   }
   uint64_t localBlocks=0;
   for(uint i=localCellStartOffset;i<localCellStartOffset+localCells;i++){
     localBlocks+=nBlocks[i];
   }

   phiprof::stop("readDatalayout");
   //todo, check file datatype, and do not just use double
   phiprof::start("readCellParameters");


   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"perturbed_B",CellParams::PERBX,3,mpiGrid); }
// This has to be set anyway, there are also the derivatives that should be written/read if we want to only read in background field
//   if(success)
//     success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"background_B",CellParams::BGBX,3,mpiGrid);
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"moments",CellParams::RHO,4,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"moments_dt2",CellParams::RHO_DT2,4,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"moments_r",CellParams::RHO_R,4,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"moments_v",CellParams::RHO_V,4,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"pressure",CellParams::P_11,3,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"pressure_dt2",CellParams::P_11_DT2,3,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"pressure_r",CellParams::P_11_R,3,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"pressure_v",CellParams::P_11_V,3,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"LB_weight",CellParams::LBWEIGHTCOUNTER,1,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"max_v_dt",CellParams::MAXVDT,1,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"max_r_dt",CellParams::MAXRDT,1,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"max_fields_dt",CellParams::MAXFDT,1,mpiGrid); }
   // Read rho losses Note: vector size = 1 (In the older versions the rho loss wasn't recorded)
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"rho_loss_adjust",CellParams::RHOLOSSADJUST,1,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"rho_loss_velocity_boundary",CellParams::RHOLOSSVELBOUNDARY,1,mpiGrid); }


   
   phiprof::stop("readCellParameters");
   phiprof::start("readBlockData");
   if(success) { success=readBlockData(file,fileCells,localCellStartOffset,localCells,nBlocks,localBlockStartOffset,localBlocks,mpiGrid); }
   phiprof::stop("readBlockData");
   if(success) { success=file.close(); }
   phiprof::stop("readGrid");

   exitOnError(success,"(RESTART) Other failure",MPI_COMM_WORLD);
   return success;

}


