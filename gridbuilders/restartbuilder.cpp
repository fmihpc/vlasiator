#include <cstdlib>
#include <iostream>
#include <stdint.h>

#include "../mpilogger.h"
#include "../parameters.h"
#include "restartbuilder.h"

using namespace std;

extern MPILogger mpilogger;

static MPI_File filePtr; /**< MPI file pointer.*/

RestartBuilder::RestartBuilder(): GridBuilder(),initialized(false) { }

RestartBuilder::~RestartBuilder() { }

bool RestartBuilder::addCellBlockDataRequests(VirtualCell::ID& totalCells,VirtualCell::ID& blockOffset,VirtualCell::ID* cellIDs,
					      uint* blocksPerCell,Real** avgsBuffer,Real** blockParamsBuffer,uint** nbrsVelBuffer) {
   bool success = true;
   
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name","avgs"));
   attribs.push_back(make_pair("mesh","SpatialGrid"));

   // Read values of distribution function to all velocity blocks:
   if (vlsvReader.multiReadStart("BLOCKVARIABLE",attribs) == false) {
      success = false;
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to get file info for array containing distrib. function values!" << endl << write;
   }   
   if (success == true) {
      for (VirtualCell::ID i=0; i<totalCells; ++i) {
	 if (vlsvReader.multiReadAddUnit(blocksPerCell[i],reinterpret_cast<char*>(avgsBuffer[i])) == false) {
	    success = false;
	    mpilogger << "(RESTARTBUILDER) ERROR: Failed to add multiread units for distrib. function!" << endl << write;
	 }
      }
   }
   if (success == true) {
      if (vlsvReader.multiReadEnd(blockOffset) == false) {
	 success = false;
	 mpilogger << "(RESTARTBUILDER) ERROR: Failed to read distrib. function values!" << endl << write;
      }
   }   
   // Read velocity block parameters:
   attribs.clear();
   attribs.push_back(make_pair("name","SpatialGrid"));
   if (vlsvReader.multiReadStart("BLOCKCOORDINATES",attribs) == false) {
      success = false;
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to get file info for array containing vel. block coordinates!" << endl << write;
   }   
   if (success == true) {
      for (VirtualCell::ID i=0; i<totalCells; ++i) {
	 if (vlsvReader.multiReadAddUnit(blocksPerCell[i],reinterpret_cast<char*>(blockParamsBuffer[i])) == false) {
	    success = false;
	    mpilogger << "(RESTARTBUILDER) ERROR: Failed to add multiread units for block coordinates!" << endl << write;
	 }
      }
   }
   if (success == true) {
      if (vlsvReader.multiReadEnd(blockOffset) == false) {
	 success = false;
	 mpilogger << "(RESTARTBUILDER) ERROR: Failed to read vel. block coordinates!" << endl << write;
      }
   }   
   // Read velocity block neighbour lists:
   attribs.clear();
   attribs.push_back(make_pair("name","SpatialGrid"));
   attribs.push_back(make_pair("mesh","SpatialGrid"));
   if (vlsvReader.multiReadStart("BLOCKNBRS",attribs) == false) {
      success = false;
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to get file info for array containing vel. block neighbours!" << endl << write;
   }
   if (success == true) {
      for (VirtualCell::ID i=0; i<totalCells; ++i) {
	 if (vlsvReader.multiReadAddUnit(blocksPerCell[i],reinterpret_cast<char*>(nbrsVelBuffer[i])) == false) {
	    success = false;
	    mpilogger << "(RESTARTBUILDER) ERROR: Failed to add multiread units for vel. block neighbours!" << endl << write;
	 }
      }
   }
   if (success == true) {
      if (vlsvReader.multiReadEnd(blockOffset) == false) {
	 success = false;
	 mpilogger << "(RESTARTBUILDER) ERROR: Failed to read vel. block neighbours!" << endl << write;
      }
   }
   return success;
}

bool RestartBuilder::addCellBlockNumberRequests(VirtualCell::ID& totalCells,VirtualCell::ID& cellOffset,VirtualCell::ID* cellIDs,uint* N_blocks) { 
   bool success = true;
   bool convertEndianness = false;
   
   bool deleteArray = false;
   char* buffer = NULL;
   uint64_t arraySize;
   uint64_t vectorSize;
   uint64_t dataSize;
   VLSV::datatype dataType;
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name","SpatialGrid"));

   // Read block number array info and create a temporary buffer for data (if necessary):
   if (vlsvReader.getArrayInfo("NBLOCKS",attribs,arraySize,vectorSize,dataType,dataSize) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read block number array info!" << endl << write;
      success = false;
   }
   if (dataSize != sizeof(uint)) {
      buffer = new char[totalCells*vectorSize*dataSize];
      deleteArray = true;
   } else {
      buffer = reinterpret_cast<char*>(N_blocks);
   }
   // Read block number array data and convert into requested datatype (if necessary):
   if (vlsvReader.readArray("NBLOCKS",attribs,cellOffset,totalCells,buffer) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read block number array data!" << endl << write;
      success = false;
   }
   if (deleteArray == true) { // convert
      VlsHeader::UInt (*converter)(const unsigned char* const,const bool&);
      if (dataSize == 1) converter = convUInt1;
      else if (dataSize == 2) converter = convUInt2;
      else if (dataSize == 4) converter = convUInt4;
      else if (dataSize == 8) converter = convUInt8;
      else {cerr << "Converter not implemented in restartbuilder!" << endl; exit(1);}
      
      const unsigned char* ptr = reinterpret_cast<const unsigned char*>(buffer);
      for (VirtualCell::ID i=0; i<totalCells; ++i) N_blocks[i] = converter(ptr + i*dataSize,convertEndianness);
      delete buffer;
      buffer = NULL;
   }
   
   return success;
}

bool RestartBuilder::addCellNbrRequests(VirtualCell::ID& totalCells,VirtualCell::ID& totalNbrs,VirtualCell::ID& cellOffset,
					VirtualCell::ID& nbrOffset,VirtualCell::ID* cellIDs,
					uchar* nbrsPerCell,Real* coords,VirtualCell::ID* nbrIDs,uchar* nbrTypes) {
   bool success = true;
   bool convertEndianness = false;

   bool deleteArray = false;
   char* buffer = NULL;
   uint64_t arraySize;
   uint64_t vectorSize;
   uint64_t dataSize;
   VLSV::datatype dataType;
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name","SpatialGrid"));
   
   // Read coordinate array info and create a temporary buffer for data (if necessary):
   if (vlsvReader.getArrayInfo("COORDS",attribs,arraySize,vectorSize,dataType,dataSize) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read coordinate array info!" << endl << write;
      success = false;
   }
   if (dataSize != sizeof(Real)) {
      buffer = new char[totalCells*vectorSize*dataSize];
      deleteArray = true;
   } else {
      buffer = reinterpret_cast<char*>(coords);
   }
   // Read coordinate array data and convert into requested data type (if necessary):
   if (vlsvReader.readArray("COORDS",attribs,cellOffset,totalCells,buffer) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read coordinate array data!" << endl << write;
      success = false;
   }
   if (deleteArray == true) { // Convert
      VlsHeader::Real (*converter)(const unsigned char* const,const bool&);
      if (dataSize == 4) converter = convReal4;
      else if (dataSize == 8) converter = convReal8;
      else if (dataSize == 16) converter = convReal16;
      else {cerr << "Converter not implemented in restartbuilder!" << endl; exit(1);}
      
      const unsigned char* ptr = reinterpret_cast<const unsigned char*>(buffer);
      for (VirtualCell::ID i=0; i<totalCells; ++i) coords[i] = converter(ptr + i*dataSize,false);      
      delete buffer;
      deleteArray = false;
   }

   // Read neighbour ID array info and create a temporary buffer for data (if necessary):
   if (vlsvReader.getArrayInfo("NBRIDS",attribs,arraySize,vectorSize,dataType,dataSize) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read neighbour ID array info!" << endl << write;
      success = false;
   }
   if (dataSize != sizeof(VirtualCell::ID)) {
      buffer = new char[totalNbrs*vectorSize*dataSize];
      deleteArray = true;
   } else {
      buffer = reinterpret_cast<char*>(nbrIDs);
   }
   // Read neighbour ID array data and convert into requested data type (if necessary):
   if (vlsvReader.readArray("NBRIDS",attribs,nbrOffset,totalNbrs,buffer) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read neighbour ID array data!" << endl << write;
      success = false;
   }
   if (deleteArray == true) { // convert
      VlsHeader::UInt (*converter)(const unsigned char* const,const bool&);
      if (dataSize == 4) converter = convUInt4;
      else if (dataSize == 8) converter = convUInt8;
      else {cerr << "Converter not implemented in restartbuilder!" << endl; exit(1);}

      const unsigned char* ptr = reinterpret_cast<const unsigned char*>(buffer);
      for (VirtualCell::ID i=0; i<totalNbrs; ++i) nbrIDs[i] = converter(ptr+i*dataSize,convertEndianness);      
      delete buffer;
      deleteArray = false;
   }
   
   // Read neighbour type array and create a temporary buffer (if necessary):
   if (vlsvReader.getArrayInfo("NBRTYPES",attribs,arraySize,vectorSize,dataType,dataSize) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read neighbour type array info!" << endl << write;
      success = false;
   }
   if (dataSize != sizeof(uchar)) {
      buffer = new char[totalNbrs*vectorSize*dataSize];
      deleteArray = true;
   } else {
      buffer = reinterpret_cast<char*>(nbrTypes);
   }
   // Read neighbour type array data and convert into requested data type (if necessary):
   if (vlsvReader.readArray("NBRTYPES",attribs,nbrOffset,totalNbrs,buffer) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read neighbour type array data!" << endl << write;
      success = false;
   }
   if (deleteArray == true) { // convert
      VlsHeader::UInt (*converter)(const unsigned char* const,const bool&);
      if (dataSize == 1) converter = convUInt1;
      else if (dataSize == 2) converter = convUInt2;
      else if (dataSize == 4) converter = convUInt4;
      else if (dataSize == 8) converter = convUInt8;
      else {cerr << "Converter not implemented in restartbuilder!" << endl; exit(1);}
      
      const unsigned char* ptr = reinterpret_cast<const unsigned char*>(buffer);
      for (VirtualCell::ID i=0; i<totalNbrs; ++i) nbrTypes[i] = converter(ptr+i*dataSize,convertEndianness);
      delete buffer;
      deleteArray = false;
   }
   return success;
}

bool RestartBuilder::addCellParamsRequests(VirtualCell::ID& totalCells,VirtualCell::ID& cellOffset,VirtualCell::ID* cellIDs,Real* cellParams) { 
   bool success = true;
   
   bool deleteArray = false;
   char* buffer = NULL;
   uint64_t arraySize;
   uint64_t vectorSize;
   uint64_t dataSize;
   VLSV::datatype dataType;
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name","SpatialGrid"));
   
   // Read cell parameter array info and create a temporary buffer (if necessary):
   if (vlsvReader.getArrayInfo("CELLPARAMS",attribs,arraySize,vectorSize,dataType,dataSize) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read cell parameter array info!" << endl << write;
      success = false;
   }
   if (dataSize != sizeof(Real)) {
      buffer = new char[totalCells*vectorSize*dataSize];
      deleteArray = true;
   } else {
      buffer = reinterpret_cast<char*>(cellParams);
   }
   // Read cell parameter array data and convert into requested data type (if necessary):
   if (vlsvReader.readArray("CELLPARAMS",attribs,cellOffset,totalCells,buffer) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read cell parameter array info!" << endl << write;
      success = false;
   }
   if (deleteArray == true) {
      VlsHeader::Real (*converter)(const unsigned char* const,const bool&);
      if (dataSize == 4) converter = convReal4;
      if (dataSize == 4) converter = convReal8;
      if (dataSize == 4) converter = convReal16;
      else {cerr << "Converter not implemented in restartbuilder!" << endl; exit(1);}
      
      const unsigned char* ptr = reinterpret_cast<const unsigned char*>(buffer);
      for (VirtualCell::ID i=0; i<totalCells; ++i) cellParams[i] = converter(ptr+i*dataSize,false);
      delete buffer;
      deleteArray = false;
   }   
   return success;
}

bool RestartBuilder::doInitialLoadBalance() {return false;}

bool RestartBuilder::finalize() {
   return initialized;
}

bool RestartBuilder::getCellIDs(std::vector<VirtualCell::ID>& cellIDs,std::vector<uchar>& N_nbrs) {
   if (myRank != masterRank) {
      mpilogger << "(RESTARTBUILDER) ERROR: It is prohibited to call getCellIDs on process #" << myRank << endl << write;
      return false;
   }
   bool success = true;
   
   // Get info on array containing cell IDs:
   uint64_t arraySize;
   uint64_t vectorSize;
   VLSV::datatype dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name","SpatialGrid"));
   if (vlsvReader.getArrayInfoMaster("MESH",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read cell ID array info!" << endl << write;
      return false;
   }
   
   // Read cell IDs:
   char* IDbuffer = new char[arraySize*vectorSize*byteSize];
   if (vlsvReader.readArrayMaster("MESH",attribs,0,arraySize,IDbuffer) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read cell IDs!" << endl << write;
      success = false;
   }
   
   // Convert global IDs into VirtualCell::ID type:
   cellIDs.resize(arraySize);
   this->N_cells = arraySize;
   N_nbrs.resize(arraySize);
   if (dataType == VLSV::UINT && byteSize == 4) {
      uint32_t* ptr = reinterpret_cast<uint32_t*>(IDbuffer);
      for (uint64_t i=0; i<arraySize; ++i) cellIDs[i] = ptr[i];
   } else if (dataType == VLSV::UINT && byteSize == 8) {
      uint64_t* ptr = reinterpret_cast<uint64_t*>(IDbuffer);
      for (uint64_t i=0; i<arraySize; ++i) cellIDs[i] = ptr[i];
   } else {
      mpilogger << "(RESTARTBUILDER) ERROR: VLSVParReader returned an unsupported datatype for cell IDs!" << endl << write;
      success = false;
   }
   delete IDbuffer;
   
   // Get info on array containing the number of spatial neighbours per cell:
   if (vlsvReader.getArrayInfoMaster("NBRSUM",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read neighbour sum array info!" << endl << write;
      return false;
   }
   
   // Read number of neighbours per cell:
   IDbuffer = new char[arraySize*vectorSize*byteSize];
   cellNbrs = new uchar[this->N_cells];
   if (vlsvReader.readArrayMaster("NBRSUM",attribs,0,arraySize,IDbuffer) == false) {
      mpilogger << "(RESTARTBUILDER) ERROR: Failed to read neighbour sum array data!" << endl << write;
      success = false;
   }
   
   // Convert number of neighbours into uchars:
   if (dataType == VLSV::UINT && byteSize == 1) {
      uint8_t* ptr = reinterpret_cast<uint8_t*>(IDbuffer);
      for (uint64_t i=0; i<arraySize; ++i) N_nbrs[i] = ptr[i];
      for (uint64_t i=0; i<arraySize; ++i) cellNbrs[i] = ptr[i];
   } else {
      mpilogger << "(RESTARTBUILDER) ERROR: VLSVParReader returned an unsupported datatype for neighbour sum array!" << endl << write;
      success = false;
   }
   delete IDbuffer;
   return success;
}

bool RestartBuilder::getParameter(const std::string& parameterName,std::string& value) {
   map<string,string>::const_iterator it = options.find(parameterName);
   if (it == options.end()) return false;
   value = it->second;
   return true;
}

bool RestartBuilder::getTotalNumberOfCells(VirtualCell::ID& N_cells) {
   if (myRank != masterRank) {
      mpilogger << "(RESTARTBUILDER) ERROR: it is prohibited to call getTotalNumberOfCells on process #" << myRank << endl << write;
      return false;
   }
   
   bool success = true;
   uint64_t arraySize;
   uint64_t vectorSize;
   VLSV::datatype dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name","SpatialGrid"));
   if (vlsvReader.getArrayInfoMaster("MESH",attribs,arraySize,vectorSize,dataType,byteSize) == false) success=false;
   N_cells = arraySize;
   return success;
}

bool RestartBuilder::initialize(MPI_Comm comm,const int& masterRank) {
   typedef Parameters P;
   typedef Readparameters RP;
   
   this->comm = comm;
   this->masterRank = masterRank;
   initialized = true;
   
   // Define parameters read from parameter file:
   RP::add("gridbuilder.q","Charge of simulated particle species, in Coulombs.",numeric_limits<Real>::max());
   RP::add("gridbuilder.m","Mass of simulated particle species, in kilograms.",numeric_limits<Real>::max());
   RP::add("gridbuilder.dt","Timestep in seconds.",numeric_limits<Real>::max());
   RP::add("gridbuilder.t_min","Simulation time at timestep 0, in seconds.",numeric_limits<Real>::max());
   RP::add("gridbuilder.timestep","Timestep when grid is loaded. Defaults to value zero.",0);
   RP::add("gridbuilder.max_timesteps","Max. value for timesteps. Defaults to value zero.",0);
   
   RP::add("restartbuilder.filename","Name of the restart file","");
   
   // Fetch parameter values and scatter to all processes. 
   // Note that q,m,dt,t_min,timestep,max_timesteps are defined in class gridbuilder.
   RP::parse();
   RP::get("gridbuilder.q",q);
   RP::get("gridbuilder.m",m);
   RP::get("gridbuilder.dt",dt);
   RP::get("gridbuilder.t_min",t_min);
   RP::get("gridbuilder.timestep",timestep);
   RP::get("gridbuilder.max_timesteps",max_timesteps);
   
   RP::get("restartbuilder.filename",fileName);
   
   // Attempt to open VLSV file for reading:
   MPI_Comm_rank(comm,&myRank);
   MPI_Comm_size(comm,&processes);
   
   MPI_Info mpiInfo = MPI_INFO_NULL;
   if (vlsvReader.open(fileName,comm,masterRank,mpiInfo) == false) {
      mpilogger << "(RESTARTBUILDER) VLSVParReader failed to open restart file '" << fileName << "' for reading!" << endl << write;
      initialized = false;
   }
   return initialized;
}

bool RestartBuilder::processCellBlockDataRequests() {return true;}

bool RestartBuilder::processCellBlockNumberRequests() {return true;}

bool RestartBuilder::processCellNbrRequests() {return true;}

bool RestartBuilder::processCellParamsRequests() {return true;}

bool RestartBuilder::waitCellBlockDataRequests() {return true;}

bool RestartBuilder::waitCellBlockNumberRequests() {return true;}

bool RestartBuilder::waitCellNbrRequests() {return true;}

bool RestartBuilder::waitCellParamsRequests() {return true;}

// Register RestartBuilder
GridBuilder* newBuilder() {return new RestartBuilder;}

class Dummy {
 public:
   Dummy() {
      GridBuilderFactory::registerBuilder(newBuilder);
   }
   ~Dummy() { }
};

static Dummy dummy;
