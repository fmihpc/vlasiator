#include <cstdlib>
#include <iostream>
#include <iomanip> // for setprecision()
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>
#include "boost/array.hpp"
#include "ioread.h"
#include "phiprof.hpp"
#include "parameters.h"
#include "logger.h"
#include "vlsvreader2.h"
#include "vlsv_reader_parallel.h"
#include "vlasovmover.h"

using namespace std;
using namespace phiprof;
using namespace vlsv;

extern Logger logFile, diagnostic;

typedef Parameters P;

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
 Read in cell ID's from file. Note: Uses the older version of vlsv parallel reader
 \param file Some vlsv reader with a file open
 \param fileCells Vector in whic to store the cell ids
 \param masterRank The simulation's master rank id (Vlasiator uses 0, which should be the default)
 \param comm MPI comm (MPI_COMM_WORLD should be the default)
*/

bool readCellIds(VLSVParReader & file,
                 vector<uint64_t>& fileCells, int masterRank,MPI_Comm comm){
   // Get info on array containing cell Ids:
   uint64_t arraySize;
   uint64_t vectorSize;
   VLSV::datatype dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   bool success=true;
   int rank;
   MPI_Comm_rank(comm,&rank);
   if(rank==masterRank){
      //let's let master read cellId's, we anyway have at max ~1e6 cells
      attribs.push_back(make_pair("name","SpatialGrid"));
      if (file.getArrayInfoMaster("MESH",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
         logFile << "(RESTARTBUILDER) ERROR: Failed to read cell ID array info!" << endl << write;
         return false;
      }
      
      //   Read cell Ids:
      char* IDbuffer = new char[arraySize*vectorSize*byteSize];
      if (file.readArrayMaster("MESH",attribs,0,arraySize,IDbuffer) == false) {
         logFile << "(RESTARTBUILDER) ERROR: Failed to read cell Ids!" << endl << write;
         success = false;
      }
   
   // Convert global Ids into our local DCCRG 64 bit uints
      fileCells.resize(arraySize);
      int N_cells = arraySize;
      if (dataType == VLSV::UINT && byteSize == 4) {
         uint32_t* ptr = reinterpret_cast<uint32_t*>(IDbuffer);
         for (uint64_t i=0; i<arraySize; ++i) fileCells[i] = ptr[i];
      } else if (dataType == VLSV::UINT && byteSize == 8) {
         uint64_t* ptr = reinterpret_cast<uint64_t*>(IDbuffer);
         for (uint64_t i=0; i<arraySize; ++i) fileCells[i] = ptr[i];
      } else {
         logFile << "(RESTARTBUILDER) ERROR: VLSVParReader returned an unsupported datatype for cell Ids!" << endl << write;
         success = false;
      }
   }
   //broadcast cellId's to everybody
   MPI_Bcast(&arraySize,1,MPI_UINT64_T,masterRank,comm);   
   fileCells.resize(arraySize);
   MPI_Bcast(&(fileCells[0]),arraySize,MPI_UINT64_T,masterRank,comm);
   
   return success;
}

/*!
 \brief Read cell ID's
 Read in cell ID's from file. Note: Uses the newer version of vlsv parallel reader
 \param file Some vlsv reader with a file open
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
         logFile << "(RESTARTBUILDER) ERROR: Failed to read cell ID array info!" << endl << write;
         return false;
      }

      //Make a routine error check:
      if( vectorSize != 1 ) {
         logFile << "(RESTARTBUILDER) ERROR: Bad vectorsize at " << __FILE__ << " " << __LINE__ << endl << write;
         return false;
      }
      
      //   Read cell Ids:
      char* IDbuffer = new char[arraySize*vectorSize*byteSize];
      if (file.readArrayMaster("VARIABLE",attribs,readFromFirstIndex,arraySize,IDbuffer) == false) {
         logFile << "(RESTARTBUILDER) ERROR: Failed to read cell Ids!" << endl << write;
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
         logFile << "(RESTARTBUILDER) ERROR: VLSVParReader returned an unsupported datatype for cell Ids!" << endl << write;
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
 \param file Some vlsv reader with a file open (can be old or new vlsv reader)
 \param nBlocks Vector for holding information on cells and the number of blocks in them -- this function saves data here
 \param masterRank The master rank of this process (Vlasiator uses masterRank = 0 and so it should be the default)
 \param comm MPI comm
 \return Returns true if the operation was successful
 \sa exec_readGrid
*/
template <class T>
bool readNBlocks( T & file,
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
      if( typeid(T) == typeid(ParallelReader) ) {
         attribs.push_back(make_pair("mesh","SpatialGrid"));
      } else if( typeid(T) == typeid(VLSVParReader) ) {
         attribs.push_back(make_pair("name","SpatialGrid"));
      } else {
         cout << "(RESTARTBUILDER) ERROR: BAD TYPEID IN READNBLOCKS AT " << __FILE__ << " " << __LINE__ << endl;
         logFile << "(RESTARTBUILDER) ERROR: BAD TYPEID IN READNBLOCKS AT " << __FILE__ << " " << __LINE__ << endl << write;
         return false;
      }
      if (file.getArrayInfoMaster("BLOCKSPERCELL",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
         logFile << "(RESTARTBUILDER) ERROR: Failed to read number of blocks" << endl << write;
         success= false;
      }

      const short int readFileFromBeginning = 0;
      nBlocks.resize(vectorSize*arraySize);
      if (file.readArrayMaster("BLOCKSPERCELL",attribs,readFileFromBeginning,arraySize,(char*)&(nBlocks[0])) == false) {
         logFile << "(RESTARTBUILDER) ERROR: Failed to read number of blocks!" << endl << write;
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

/*! This is a function for calling templated functions whose template is determined by data type and size read by vlsv reader. Note: 
 \param
 \param

*/
bool executeWithTemplate( 7 ) {

}

/*! This reads in data one cell at a time. It is not the most efficient way but has the following benefits
 - For large datasets (distribution function), we avoid any problem with having to store all distribution functions twice in memory
 - Machinery in readvlsv does not at the moment support setting fileviews, this should be improved.
 The template stands for the file type so if one is reading doubles, fileReal should be double
 TODO: Get rid of fileReal (Was done once in a branch but it caused problems)
 \param file Some vlsv reader with a file open
 \param fileCells List of all cell ids
 \param localCellStartOffset The offset from which to start reading cells ( This should be balanced so that every process has roughly the same amount of blocks to read )
 \param localCells How many cells after the offset to read ( This should be balanced so that every process has roughly the same amount of blocks to read )
 \param localBlockStartOffset localCellStartOffset's corresponding block offset in this process ( Calculated from nBlocks and localCellStartOffset )
 \param localCells localCellStartOffset's corresponding block block amount in this process ( Calculated from nBlocks and localCellStartOffset and localCellStartOffset )
 \param mpiGrid Vlasiator's grid
 \sa exec_readGrid
*/

template <typename fileReal>
bool readBlockData(
   VLSVParReader & file,
   const vector<uint64_t>& fileCells,
   const uint64_t localCellStartOffset,
   const uint64_t localCells,
   const vector<uint>& nBlocks,
   const uint64_t localBlockStartOffset,
   const uint64_t localBlocks,
   dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid
) {
  uint64_t arraySize;
  uint64_t avgVectorSize;
  uint64_t coordVectorSize;
  uint64_t cellParamsVectorSize;
  VLSV::datatype dataType;
  uint64_t byteSize;
  list<pair<string,string> > avgAttribs;
  list<pair<string,string> > coordAttribs;
  fileReal *coordBuffer;
  fileReal *avgBuffer;
  bool success=true;
   
  coordAttribs.push_back(make_pair("name","SpatialGrid"));
  avgAttribs.push_back(make_pair("name","avgs"));
  avgAttribs.push_back(make_pair("mesh","SpatialGrid"));
  

  //Get array info for cell parameters 
  if (file.getArrayInfo("BLOCKCOORDINATES",coordAttribs,arraySize,coordVectorSize,dataType,byteSize) == false ){
    logFile << "(RESTARTBUILDER) ERROR: Failed to read BLOCKCOORDINATES array info " << endl << write;
    return false;
  }

  if(file.getArrayInfo("BLOCKVARIABLE",avgAttribs,arraySize,avgVectorSize,dataType,byteSize) == false ){
    logFile << "(RESTARTBUILDER) ERROR: Failed to read BLOCKVARIABLE array info " << endl << write;
    return false;
  }

  //todo: more errorchecks!!  
   if(avgVectorSize!=WID3){
      logFile << "(RESTARTBUILDER) ERROR: Blocksize does not match in restart file " << endl << write;
      return false;
   }
      
   coordBuffer=new fileReal[coordVectorSize*localBlocks];
   avgBuffer=new fileReal[avgVectorSize*localBlocks];
   
   file.readArray("BLOCKCOORDINATES",coordAttribs,localBlockStartOffset,localBlocks,(char*)coordBuffer);
   file.readArray("BLOCKVARIABLE",avgAttribs,localBlockStartOffset,localBlocks,(char*)avgBuffer);
   
   uint64_t bufferBlock=0;
   for(uint i=0;i<localCells;i++){
     uint cell=fileCells[localCellStartOffset+i];
     for (uint blockIndex=0;blockIndex<nBlocks[localCellStartOffset+i];blockIndex++){
        creal vx_block = coordBuffer[bufferBlock*coordVectorSize+BlockParams::VXCRD];
        creal vy_block = coordBuffer[bufferBlock*coordVectorSize+BlockParams::VYCRD];
        creal vz_block = coordBuffer[bufferBlock*coordVectorSize+BlockParams::VZCRD];
        creal dvx_blockCell = coordBuffer[bufferBlock*coordVectorSize+BlockParams::DVX];
        creal dvy_blockCell = coordBuffer[bufferBlock*coordVectorSize+BlockParams::DVY];
        creal dvz_blockCell = coordBuffer[bufferBlock*coordVectorSize+BlockParams::DVZ];
        // set    volume average of distrib. function for each cell in the block.
        for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
           creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
           creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
           creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
           //todo, use faster set_value interface
           mpiGrid[cell]->set_value(vx_cell_center,vy_cell_center,vz_cell_center,avgBuffer[bufferBlock*avgVectorSize+cellIndex(ic,jc,kc)]);
        }
        bufferBlock++; 
     }
   }


   delete(avgBuffer);
   delete(coordBuffer);
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


/*! This reads in data one cell at a time. It is not the most efficient way but has the following benefits
 - For large datasets (distribution function), we avoid any problem with having to store all distribution functions twice in memory
 - Machinery in readvlsv does not at the moment support setting fileviews, this should be improved.
 The template stands for the file type so if one is reading doubles, fileReal should be double
 TODO: Get rid of fileReal (Was done once in a branch but it caused problems)
 \param file Some vlsv reader with a file open
 \param fileCells List of all cell ids
 \param localCellStartOffset The offset from which to start reading cells ( This should be balanced so that every process has roughly the same amount of blocks to read )
 \param localCells How many cells after the offset to read ( This should be balanced so that every process has roughly the same amount of blocks to read )
 \param localBlockStartOffset localCellStartOffset's corresponding block offset in this process ( Calculated from nBlocks and localCellStartOffset )
 \param localCells localCellStartOffset's corresponding block block amount in this process ( Calculated from nBlocks and localCellStartOffset and localCellStartOffset )
 \param mpiGrid Vlasiator's grid
 \sa exec_readGrid
*/
template <typename fileReal>
bool readBlockData(
   ParallelReader & file,
   const vector<uint64_t>& fileCells,
   const uint64_t localCellStartOffset,
   const uint64_t localCells,
   const vector<uint>& nBlocks,
   const uint64_t localBlockStartOffset,
   const uint64_t localBlocks,
   dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid
) {
  uint64_t arraySize;
  uint64_t avgVectorSize;
  uint64_t cellParamsVectorSize;
  datatype::type dataType;
  uint64_t byteSize;
  list<pair<string,string> > avgAttribs;
  fileReal *avgBuffer;
  bool success=true;
   
  avgAttribs.push_back(make_pair("name","avgs"));
  avgAttribs.push_back(make_pair("mesh","SpatialGrid"));
  


  //Get block id array info and store them into blockIdAttribs, lockIdByteSize, blockIdDataType, blockIdVectorSize
  list<pair<string,string> > blockIdAttribs;
  uint64_t blockIdVectorSize, blockIdByteSize;
  datatype::type blockIdDataType;
  blockIdAttribs.push_back( make_pair("mesh", "SpatialGrid") );
  if (file.getArrayInfo("BLOCKIDS",blockIdAttribs,arraySize,blockIdVectorSize,blockIdDataType,blockIdByteSize) == false ){
    logFile << "(RESTARTBUILDER) ERROR: Failed to read BLOCKCOORDINATES array info " << endl << write;
    return false;
  }

  if(file.getArrayInfo("BLOCKVARIABLE",avgAttribs,arraySize,avgVectorSize,dataType,byteSize) == false ){
    logFile << "(RESTARTBUILDER) ERROR: Failed to read BLOCKVARIABLE array info " << endl << write;
    return false;
  }

   //Some routine error checks:
   if( avgVectorSize!=WID3 ){
      logFile << "(RESTARTBUILDER) ERROR: Blocksize does not match in restart file " << endl << write;
      return false;
   }
   if( byteSize != sizeof(fileReal) ) {
      logFile << "(RESTARTBUILDER) ERROR: Bad avgs bytesize at " << __FILE__ << " " << __LINE__ << endl << write;
      return false;
   }
      
   avgBuffer=new fileReal[avgVectorSize*localBlocks];

   //Create a buffer for writing in block ids
   char * blockIdBuffer_char = new char[blockIdVectorSize * localBlocks * blockIdByteSize / sizeof(char)];
  
   //Read block ids into blockIdBuffer_char
   file.readArray( "BLOCKIDS", blockIdAttribs, localBlockStartOffset, localBlocks, blockIdBuffer_char );

   file.readArray("BLOCKVARIABLE",avgAttribs,localBlockStartOffset,localBlocks,(char*)avgBuffer);

   //Get velocity cell lengths for iteration:
   const unsigned short int numberOfCellsInBlocksPerDirection = 4;
   const creal dvx_blockCell = ((P::vxmax - P::vxmin) / P::vxblocks_ini) / (creal)(numberOfCellsInBlocksPerDirection);
   const creal dvy_blockCell = ((P::vymax - P::vymin) / P::vyblocks_ini) / (creal)(numberOfCellsInBlocksPerDirection);
   const creal dvz_blockCell = ((P::vzmax - P::vzmin) / P::vzblocks_ini) / (creal)(numberOfCellsInBlocksPerDirection);

   //Iterate through blocks:
   uint64_t bufferBlock=0;
   for(uint i=0;i<localCells;i++){
     //Go through all spatial cells
     uint cell=fileCells[localCellStartOffset+i];
     for (uint blockIndex=0;blockIndex<nBlocks[localCellStartOffset+i];blockIndex++){
        //Get the block id:
        uint64_t blockId;
        //Note: blockid's datatype is uint
        const short unsigned int cellsInBlocksPerDirection = 4;
        if( blockIdByteSize == sizeof(unsigned int) ) {
           unsigned int * blockIds = reinterpret_cast<unsigned int*>(blockIdBuffer_char);
           blockId = blockIds[bufferBlock];
        } else {
           //blockIds_dataSize == sizeof(uint64_t)
           uint64_t * blockIds = reinterpret_cast<uint64_t*>(blockIdBuffer_char);
           blockId = blockIds[bufferBlock];
        }
        //Get the block's coordinates (min coordinates)
        boost::array<Real, 3> blockCoordinates;
        getVelocityBlockCoordinates( blockId, blockCoordinates );

        // set    volume average of distrib. function for each cell in the block.
        for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
           creal vx_cell_center = blockCoordinates[0] + (ic+convert<Real>(0.5))*dvx_blockCell;
           creal vy_cell_center = blockCoordinates[1] + (jc+convert<Real>(0.5))*dvy_blockCell;
           creal vz_cell_center = blockCoordinates[2] + (kc+convert<Real>(0.5))*dvz_blockCell;
           //TODO: use faster set_value
           mpiGrid[cell]->set_value(vx_cell_center,vy_cell_center,vz_cell_center,avgBuffer[bufferBlock*avgVectorSize+cellIndex(ic,jc,kc)]);
        }
        bufferBlock++; 
     }
   }


   delete(avgBuffer);
   delete[] blockIdBuffer_char;
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
template <typename fileReal, class U>
bool readCellParamsVariable(U & file,
			    const vector<uint64_t>& fileCells,
                            const uint64_t localCellStartOffset,
			    const uint64_t localCells,
			    const string& variableName,
                            const size_t cellParamsIndex,
                            const size_t expectedVectorSize,
                            dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid){
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
      logFile << "(RESTARTBUILDER)  ERROR: Failed to read " << endl << write;
      return false;
   }

   if(vectorSize!=expectedVectorSize){
      logFile << "(RESTARTBUILDER)  vectorsize wrong " << endl << write;
      return false;
   }
   
   buffer=new fileReal[vectorSize*localCells];
   if(file.readArray("VARIABLE",attribs,localCellStartOffset,localCells,(char *)buffer) == false ) {
      logFile << "(RESTARTBUILDER)  ERROR: Failed to read " << variableName << endl << write;
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



/*! A function for reading parameters e.g. 'timestep'
 \param file Some vlsv parallel reader with a file open
 \param name Name of the parameter
 \param value Variable in which to store the scalar variable (double, float, int .. )
 \param masterRank The master process' id (Vlasiator uses 0 so this should equal 0 by default)
 \param comm MPI comm (MPI_COMM_WORLD should be the default)
 \return Returns true if the operation is successful
 */
template <typename T>
bool readScalarParameter(VLSVParReader & file, string name,T& value, int masterRank,MPI_Comm comm){
   bool success=true;
   int myRank;
   uint64_t arraySize;
   uint64_t vectorSize;
   VLSV::datatype dataType;
   uint64_t byteSize;
   MPI_Comm_rank(comm,&myRank);
   if(myRank==masterRank){
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",name));

      if (file.getArrayInfoMaster("PARAMETERS",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
         logFile << "(RESTART)  ERROR: Failed to read info for parameter"<< name << endl << write;
         success=false;
      }
      
      if(vectorSize!=1 || arraySize!=1){
         logFile << "(RESTART) Parameter not scalar" << endl << write;
         success=false;
      }
      
      if(dataType == VLSV::INT && byteSize == 4) {
         int32_t buffer;         
         if(file.readArrayMaster("PARAMETERS",attribs,0,1,(char *)&buffer) == false ) {
            logFile << "(RESTART)  ERROR: Failed to read value for parameter"<< name << endl << write;
            success=false;
         }
         value=buffer;
      }

      else if(dataType == VLSV::UINT && byteSize == 4) {
         uint32_t buffer;         
         if(file.readArrayMaster("PARAMETERS",attribs,0,1,(char *)&buffer) == false ) {
            logFile << "(RESTART)  ERROR: Failed to read value for parameter"<< name << endl << write;
            success=false;
         }
         value=buffer;
      }
      else if(dataType == VLSV::UINT && byteSize == 8) {
         uint64_t buffer;         
         if(file.readArrayMaster("PARAMETERS",attribs,0,1,(char *)&buffer) == false ) {
            logFile << "(RESTART)  ERROR: Failed to read value for parameter"<< name << endl << write;
            success=false;
         }
         value=buffer;
      }
      else  if(dataType == VLSV::FLOAT && byteSize == 4) {
         float buffer;         
         if(file.readArrayMaster("PARAMETERS",attribs,0,1,(char *)&buffer) == false ) {
            logFile << "(RESTART)  ERROR: Failed to read value for parameter"<< name << endl << write;
            success=false;
         }
         value=buffer;
      }
      else  if(dataType == VLSV::FLOAT && byteSize == 8) {
         double buffer;         
         if(file.readArrayMaster("PARAMETERS",attribs,0,1,(char *)&buffer) == false ) {
            logFile << "(RESTART)  ERROR: Failed to read value for parameter"<< name << endl << write;
            success=false;
         }
         value=buffer;
      }
      else {
         logFile << "(RESTART) Unsupported parameter type"<< name << endl << write;
         success=false;
      }
   }
   
   
   MPI_Bcast(&value,sizeof(T),MPI_BYTE,masterRank,comm);
   return success;
}

/*! A function for reading parameters e.g. 'timestep'
 \param file Some vlsv parallel reader with a file open
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
template <class U, typename T>
bool checkScalarParameter(U & file, const string & name, T correctValue, int masterRank,MPI_Comm comm){
   T value;
   readScalarParameter(file,name,value,masterRank,comm);
   if(value!=correctValue){
      std::ostringstream s;
      s << "(RESTART) Parameter " << name << " has missmatching value.";
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

/*! A function for checking the version of the file
 \param vlsvReader Some vlsv reader with a file open
 \return Returns the version number
 */
float checkVersion( Reader & vlsvReader ) {
   string versionTag = "version";
   float version;
   if( vlsvReader.readParameter( versionTag, version ) == false ) {
      return 0;
   }
   if( version == 1.00 ) {
      return version;
   } else {
      cerr << "Invalid version!" << endl;
      exit(1);
      return 0;
   }
}

/*! This function is used to read the restart file. The template class T can be ParallelReader or VLSVParReader
 \param mpiGrid Vlasiator's grid
 \param name Name of restart file
 \return Returns true if the operation was successful
 \sa readGrid
 */
template <class T>
bool exec_readGrid(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
                   const std::string& name) {
   vector<uint64_t> fileCells; /*< CellIds for all cells in file*/
   vector<uint> nBlocks;/*< Number of blocks for all cells in file*/
   bool success=true;
   int myRank,processes;
   
   
   // Attempt to open VLSV file for reading:
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   MPI_Comm_size(MPI_COMM_WORLD,&processes);
   
   
   phiprof::start("readGrid");

   T file;
   MPI_Info mpiInfo = MPI_INFO_NULL;

   const short int masterRank = 0;
   if (file.open(name,MPI_COMM_WORLD,masterRank,mpiInfo) == false) {
      success=false;
   }
   exitOnError(success,"(RESTART) Could not open file",MPI_COMM_WORLD);

   if(readScalarParameter(file,"t",P::t,0,MPI_COMM_WORLD) ==false) success=false;//CONT
   //if( file.readParameter( "t", P::t ) == false ) success = false;
   P::t_min=P::t;

   //FIXME: If we use the dt we read in then the restarted simulation
   //has much greater deviation from original trajectory-> do we have
   //a latent bug, is there something we do not read in?
   //         if(readScalarParameter(file,"dt",P::dt,0,MPI_COMM_WORLD) ==false) success=false;
   
   if( readScalarParameter(file,"tstep",P::tstep,0,MPI_COMM_WORLD) ==false) success=false;
   //if( file.readParameter( "tstep", P::tstep ) == false ) success = false;
   P::tstep_min=P::tstep;
   if(readScalarParameter(file,"dt",P::dt,0,MPI_COMM_WORLD) ==false) success=false;
   //if( file.readParameter( "dt", P::dt ) == false ) success = false;
   checkScalarParameter(file,"xmin",P::xmin,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"ymin",P::ymin,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"zmin",P::zmin,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"xmax",P::xmax,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"ymax",P::ymax,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"zmax",P::zmax,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"xcells_ini",P::xcells_ini,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"ycells_ini",P::ycells_ini,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"zcells_ini",P::zcells_ini,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"vxmin",P::vxmin,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"vymin",P::vymin,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"vzmin",P::vzmin,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"vxmax",P::vxmax,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"vymax",P::vymax,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"vzmax",P::vzmax,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"vxblocks_ini",P::vxblocks_ini,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"vyblocks_ini",P::vyblocks_ini,masterRank,MPI_COMM_WORLD);
   checkScalarParameter(file,"vzblocks_ini",P::vzblocks_ini,masterRank,MPI_COMM_WORLD);

   phiprof::start("readDatalayout");
   if(success) { success=readCellIds(file,fileCells,masterRank,MPI_COMM_WORLD); }

   //check that the cellID lists are identical in file and grid
   if(myRank==0){
      vector<uint64_t> allGridCells=mpiGrid.get_all_cells();
      if(fileCells.size() != allGridCells.size()){
         success=false;
      }
   }
   
   exitOnError(success,"(RESTART) Wrong number of cells in restartfile",MPI_COMM_WORLD);
   if(success) {
      success = readNBlocks(file,nBlocks,masterRank,MPI_COMM_WORLD);
   }
   //make sure all cells are empty, we will anyway overwrite everything and in that case moving cells is easier...
   vector<uint64_t> gridCells = mpiGrid.get_cells();
   for(uint i=0;i<gridCells.size();i++){
      mpiGrid[gridCells[i]]->clear();
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
   //get new list of local gridcells
   gridCells = mpiGrid.get_cells();
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
      mpiGrid[gridCells[i]]->parameters[CellParams::XCRD] = mpiGrid.get_cell_x_min(gridCells[i]);
      mpiGrid[gridCells[i]]->parameters[CellParams::YCRD] = mpiGrid.get_cell_y_min(gridCells[i]);
      mpiGrid[gridCells[i]]->parameters[CellParams::ZCRD] = mpiGrid.get_cell_z_min(gridCells[i]);
      mpiGrid[gridCells[i]]->parameters[CellParams::DX  ] = mpiGrid.get_cell_length_x(gridCells[i]);
      mpiGrid[gridCells[i]]->parameters[CellParams::DY  ] = mpiGrid.get_cell_length_y(gridCells[i]);
      mpiGrid[gridCells[i]]->parameters[CellParams::DZ  ] = mpiGrid.get_cell_length_z(gridCells[i]);
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


   if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"perturbed_B",CellParams::PERBX,3,mpiGrid); }
// This has to be set anyway, there are also the derivatives tahat should be written/read if we want to only read in bancground field
//   if(success)
//     success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"background_B",CellParams::BGBX,3,mpiGrid);
   if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"moments",CellParams::RHO,4,mpiGrid); }
   if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"moments_dt2",CellParams::RHO_DT2,4,mpiGrid); }
   if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"moments_r",CellParams::RHO_R,4,mpiGrid); }
   if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"moments_v",CellParams::RHO_V,4,mpiGrid); }
   if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"LB_weight",CellParams::LBWEIGHTCOUNTER,1,mpiGrid); }
   if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"max_v_dt",CellParams::MAXVDT,1,mpiGrid); }
   if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"max_r_dt",CellParams::MAXRDT,1,mpiGrid); }
   if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"max_fields_dt",CellParams::MAXFDT,1,mpiGrid); }
   if( typeid(T) == typeid(ParallelReader) ) {
      // Read rho losses Note: vector size = 1 (In the older versions the rho loss wasn't recorded)
      if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"rho_loss_adjust",CellParams::RHOLOSSADJUST,1,mpiGrid); }
      if(success) { success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"rho_loss_velocity_boundary",CellParams::RHOLOSSVELBOUNDARY,1,mpiGrid); }
   }
   
   phiprof::stop("readCellParameters");
   phiprof::start("readBlockData");
   if(success) { success=readBlockData<double>(file,fileCells,localCellStartOffset,localCells,nBlocks,localBlockStartOffset,localBlocks,mpiGrid); }
   phiprof::stop("readBlockData");
   if(success) { success=file.close(); }
   phiprof::stop("readGrid");

   exitOnError(success,"(RESTART) Other failure",MPI_COMM_WORLD);
   return success;
}

//FIXME, readGrid has no support for checking or converting endianness
/*!
\brief Read in state from a vlsv file in order to restart simulations
\param mpiGrid Vlasiator's grid
\param name Name of the restart file e.g. "restart.00052.vlsv"
*/
bool readGrid(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
              const std::string& name){
   Reader vlsvCheck;
   vlsvCheck.open( name );
   //Check the vlsv version from the file:
   const bool newLib = (checkVersion( vlsvCheck ) == 1.00);
   vlsvCheck.close();
   if( newLib ) {
      return exec_readGrid<ParallelReader>(mpiGrid, name);
   } else {
      return exec_readGrid<VLSVParReader>(mpiGrid, name);
   }
   return false;
}


