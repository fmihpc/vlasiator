#include <cstdlib>
#include <iostream>
#include <iomanip> // for setprecision()
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>
#include "fileio.h"
#include "phiprof.hpp"
#include "logger.h"
#include "vlsvwriter2.h"
#include "vlsvreader2.h"
#include "vlasovmover.h"

using namespace std;
using namespace phiprof;

extern Logger logFile, diagnostic;


bool exitOnError(bool success,string message,MPI_Comm comm){
   int successInt;
   if(success)
      successInt=1;
   else
      successInt=0;
   
   MPI_Bcast(&successInt,1,MPI_INT,0,MPI_COMM_WORLD);
   if(successInt==1) {
      return true;
   }
   else{
      logFile << message << endl<<write ;
      exit(1);
   }
}


/*!
 \brief Read cell ID's
 Read in cell ID's from file
*/

bool readCellIds(VLSVParReader & file,
                 vector<uint64_t>& fileCells){
   // Get info on array containing cell Ids:
   uint64_t arraySize;
   uint64_t vectorSize;
   VLSV::datatype dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   bool success=true;
   
   //Get array info for 
   attribs.push_back(make_pair("name","SpatialGrid"));
   if (file.getArrayInfo("MESH",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      logFile << "(RESTARTBUILDER) ERROR: Failed to read cell ID array info!" << endl << write;
      return false;
   }
   
   // Read cell Ids:
   char* IDbuffer = new char[arraySize*vectorSize*byteSize];
   if (file.readArray("MESH",attribs,0,arraySize,IDbuffer) == false) {
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
   return success;
}


/*!
 \brief Read number of blocks per cell

*/


bool readNBlocks(VLSVParReader & file,
                 vector<unsigned int>& nBlocks){
   // Get info on array containing cell Ids:
   uint64_t arraySize;
   uint64_t vectorSize;
   VLSV::datatype dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   bool success=true;
   
   //Get array info for 
   attribs.push_back(make_pair("name","Blocks"));
   attribs.push_back(make_pair("mesh","SpatialGrid"));
   if (file.getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      logFile << "(RESTARTBUILDER) ERROR: Failed to read number of blocks" << endl << write;
      return false;
   }

   
   nBlocks.resize(vectorSize*arraySize);
   if (file.readArray("VARIABLE",attribs,0,arraySize,(char*)&(nBlocks[0])) == false) {
      logFile << "(RESTARTBUILDER) ERROR: Failed to read number of blocks!" << endl << write;
      success = false;
   }

   
   
   return success;
}


/*   
     This reads in data one cell at a time. It is not the most efficient way but has the following benefits
     - For large datasets (distribution function), we avoid any problem with having to store all distribution functions twice in memory
     - Machinery in readvlsv does not at the moment support setting fileviews, this should be improved.
     
*/
template <typename fileReal>
bool readBlockData(VLSVParReader & file,
                   const vector<uint64_t>& fileCells,
		   const uint localCellStartOffset,
		   const uint localCells,
		   const vector<uint>& nBlocks,
		   const uint localBlockStartOffset,
		   const uint localBlocks,
                   dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid){
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
  int maxNumberOfBlocks;
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
           mpiGrid[cell]->set_value(vx_cell_center,vy_cell_center,vz_cell_center,avgBuffer[bufferBlock*avgVectorSize+cellIndex(ic,jc,kc)]);
        }
        bufferBlock++; 
     }
   }


   delete(avgBuffer);
   delete(coordBuffer);
   return success;
}



template <typename fileReal>
bool readCellParamsVariable(VLSVParReader & file,
			    const vector<uint64_t>& fileCells,
                            const uint localCellStartOffset,
			    const uint localCells,
			    const string& variableName,
                            const size_t cellParamsIndex,
                            const size_t expectedVectorSize,
                            dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid){
   uint64_t arraySize;
   uint64_t vectorSize;
   VLSV::datatype dataType;
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
   
   delete(buffer);
   return success;
}

bool readGrid(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
              const std::string& name){
   vector<uint64_t> fileCells; /*< CellIds for all cells in file*/
   vector<uint> nBlocks;/*< Number of blocks for all cells in file*/
   bool success=true;
   int myRank,processes;
   
   
   // Attempt to open VLSV file for reading:
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   MPI_Comm_size(MPI_COMM_WORLD,&processes);
   
   
   phiprof::start("readGrid");

   VLSVParReader file;
   MPI_Info mpiInfo = MPI_INFO_NULL;

   
   if (file.open(name,MPI_COMM_WORLD,0,mpiInfo) == false) {
      success=false;
   }
   exitOnError(success,"(RESTART) Could not open file",MPI_COMM_WORLD);


   phiprof::start("readDatalayout");
   if(success) 
      success=readCellIds(file,fileCells);

   //check that the cellID lists are identical in file and grid
   if(myRank==0){
      vector<uint64_t> allGridCells=mpiGrid.get_all_cells();
      if(fileCells.size() != allGridCells.size()){
         success=false;
      }
   }
   
   exitOnError(success,"(RESTART) Wrong number of cells in restartfile",MPI_COMM_WORLD);
   
   if(success) 
      success=readNBlocks(file,nBlocks);
   
   //make sure all cells are empty, we will anyway overwrite everything and in that case moving cells is easier...
   vector<uint64_t> gridCells = mpiGrid.get_cells();
   for(uint i=0;i<gridCells.size();i++){
      mpiGrid[gridCells[i]]->clear();
   }
   

   //prepare to migrate cells so that each process has its cells contiguously in the file
   //First processA processes get cellsPerProcessA cells, the next cellsPerProcessB cells
   uint cellsPerProcessA=fileCells.size()/processes;
   uint cellsPerProcessB=cellsPerProcessA+1;
   uint processesB=fileCells.size()%processes;
   uint processesA=processes-processesB;
   
   //pin local cells to remote processes
   for(uint i=0;i<fileCells.size();i++){
      if(mpiGrid.is_local(fileCells[i])){
         uint newCellProcess;
         if(i>cellsPerProcessA*processesA)
            newCellProcess=processesA + (i-cellsPerProcessA*processesA)/cellsPerProcessB;
         else
            newCellProcess=i/cellsPerProcessA;
         mpiGrid.pin(fileCells[i],newCellProcess);
      }
   }
   
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_DATA);
   mpiGrid.migrate_cells();
   //get new list of local gridcells
   gridCells = mpiGrid.get_cells();
   
   //this is where local cells start in file-list after migration
   uint localCellStartOffset;
   uint localCells;
   if(myRank < processesA) {
      localCells=cellsPerProcessA;
      localCellStartOffset=cellsPerProcessA*myRank;
   }
   else {
      localCells=cellsPerProcessB;
      localCellStartOffset=cellsPerProcessA*processesA +
         cellsPerProcessB*(myRank-processesA);
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
      mpiGrid[gridCells[i]]->parameters[CellParams::DX  ] = mpiGrid.get_cell_x_size(gridCells[i]);
      mpiGrid[gridCells[i]]->parameters[CellParams::DY  ] = mpiGrid.get_cell_y_size(gridCells[i]);
      mpiGrid[gridCells[i]]->parameters[CellParams::DZ  ] = mpiGrid.get_cell_z_size(gridCells[i]);
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

   if(success)
     success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"B",CellParams::BX,3,mpiGrid);
   if(success)
     success=readCellParamsVariable<double>(file,fileCells,localCellStartOffset,localCells,"B0",CellParams::BX0,3,mpiGrid);

   phiprof::stop("readCellParameters");
   phiprof::start("readBlockData");
   if(success) 
      success=readBlockData<double>(file,fileCells,localCellStartOffset,localCells,nBlocks,localBlockStartOffset,localBlocks,mpiGrid);
   phiprof::stop("readBlockData");
   if(success)
      success=file.close();
   phiprof::stop("readGrid");

   if(success)
      calculateVelocityMoments(mpiGrid);

   exitOnError(success,"(RESTART) Other failure",MPI_COMM_WORLD);
   return success;
}



//-------------------------------------
//write functions

bool writeDataReducer(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                      DataReducer& dataReducer,
                      int dataReducerIndex,
                      VLSVWriter& vlsvWriter){
   map<string,string> attribs;                      
   string variableName,dataType;
   bool success=true;
   
   uint dataSize,vectorSize;
   attribs["mesh"] = "SpatialGrid";
   variableName = dataReducer.getName(dataReducerIndex);
   if (dataReducer.getDataVectorInfo(dataReducerIndex,dataType,dataSize,vectorSize) == false) {
      cerr << "ERROR when requesting info from DRO " << dataReducerIndex << endl;
      return false;
   }
   vector<uint64_t> cells = mpiGrid.get_cells();
   uint64_t arraySize = cells.size()*vectorSize*dataSize;
   
   //Request DataReductionOperator to calculate the reduced data for all local cells:
   char* varBuffer = new char[arraySize];
   for (uint64_t cell=0; cell<cells.size(); ++cell) {
      if (dataReducer.reduceData(mpiGrid[cells[cell]],dataReducerIndex,varBuffer + cell*vectorSize*dataSize) == false){
         success = false;
      }
      
      if (success == false){
         logFile << "(MAIN) writeGrid: ERROR datareductionoperator '" << dataReducer.getName(dataReducerIndex) <<
            "' returned false!" << endl << writeVerbose;
      }
   }
   // Write  reduced data to file if DROP was successful:
   if(success){
      if (vlsvWriter.writeArray("VARIABLE",variableName,attribs,cells.size(),vectorSize,dataType,dataSize,varBuffer) == false)
         success = false;
      if (success == false){
         logFile << "(MAIN) writeGrid: ERROR failed to write datareductionoperator data to file!" << endl << writeVerbose;
      }
   }
   delete[] varBuffer;
   varBuffer = NULL;
   return success;
}

bool writeGrid(const dccrg::Dccrg<SpatialCell>& mpiGrid,
               DataReducer& dataReducer,
               const string& name,
               const uint& index,
               const bool& writeRestart) {
    double allStart = MPI_Wtime();
    bool success = true;
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
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
   
   // Get all local cell Ids and write to file:
   map<string,string> attribs;
   vector<uint64_t> cells = mpiGrid.get_cells();

   if (vlsvWriter.writeArray("MESH","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) {
      cerr << "Proc #" << myRank << " failed to write cell Ids!" << endl;
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
      cerr << "Proc #" << myRank << " failed to write cell coords!" << endl;
   }
   delete[] buffer;   

   
   if(writeRestart == false ) {
   
      // Write variables calculate d by DataReductionOperators (DRO). We do not know how many 
      // numbers each DRO calculates, so a buffer has to be re-allocated for each DRO:
      for (uint i=0; i<dataReducer.size(); ++i) {
         writeDataReducer(mpiGrid,dataReducer,i,vlsvWriter);
      }
      
      phiprof::initializeTimer("Barrier","MPI","Barrier");
      phiprof::start("Barrier");
      MPI_Barrier(MPI_COMM_WORLD);
      phiprof::stop("Barrier");
      vlsvWriter.close();
      phiprof::stop("writeGrid-reduced");
   }
   else {
      //write restart
      uint64_t totalBlocks = 0;  
      for(size_t cell=0;cell<cells.size();++cell){
         totalBlocks+=mpiGrid[cells[cell]]->size();
      }
      //write out DROs we need for restarts
      DataReducer restartReducer;
      restartReducer.addOperator(new DRO::VariableB);
      restartReducer.addOperator(new DRO::VariableB0);
      restartReducer.addOperator(new DRO::Blocks);

      for (uint i=0; i<restartReducer.size(); ++i) {
         writeDataReducer(mpiGrid,restartReducer,i,vlsvWriter);
      }
      
      // Write velocity blocks and related data. 
      // In restart we just write velocity grids for all cells.
      // First write global Ids of those cells which write velocity blocks (here: all cells):
      if (vlsvWriter.writeArray("CELLSWITHBLOCKS","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) success = false;
      if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;
      //Write velocity block coordinates.
      std::vector<Real> velocityBlockParameters;
      velocityBlockParameters.reserve(totalBlocks*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      
      //gather data for writing
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
      
      attribs.clear();
      if (vlsvWriter.writeArray("BLOCKCOORDINATES","SpatialGrid",attribs,totalBlocks,BlockParams::N_VELOCITY_BLOCK_PARAMS,&(velocityBlockParameters[0])) == false) success = false;
      if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write BLOCKCOORDINATES to file!" << endl << writeVerbose;
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
      
      attribs.clear();
      attribs["mesh"] = "SpatialGrid";
      if (vlsvWriter.writeArray("BLOCKVARIABLE","avgs",attribs,totalBlocks,SIZE_VELBLOCK,&(velocityBlockData[0])) == false) success=false;
      if (success ==false)      logFile << "(MAIN) writeGrid: ERROR occurred when writing BLOCKVARIABLE avgs" << endl << writeVerbose;
      velocityBlockData.clear();
      vlsvWriter.close();

      phiprof::stop("writeGrid-restart");//,1.0e-6*bytesWritten,"MB");

   }

   return success;
}
   



        





bool writeDiagnostic(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                     DataReducer& dataReducer,
                     luint tstep){
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   
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
   
   if (printDiagnosticHeader == true && myRank == 0) {
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
      
      if (success == false) logFile << "(MAIN) computeDiagnostic: ERROR datareductionoperator '" << dataReducer.getName(i) << "' returned false!" << endl << writeVerbose;
   }
   
   MPI_Reduce(&localMin[0], &globalMin[0], nOps, MPI_Type<Real>(), MPI_MIN, 0, MPI_COMM_WORLD);
   MPI_Reduce(&localMax[0], &globalMax[0], nOps, MPI_Type<Real>(), MPI_MAX, 0, MPI_COMM_WORLD);
   MPI_Reduce(&localSum[0], &globalSum[0], nOps + 1, MPI_Type<Real>(), MPI_SUM, 0, MPI_COMM_WORLD);

   diagnostic << setprecision(12); 
   diagnostic << tstep << "\t";
   
for (uint i=0; i<nOps; ++i) {
      if (globalSum[0] != 0.0) globalAvg[i] = globalSum[i+1] / globalSum[0];
      else globalAvg[i] = globalSum[i+1];
      if (myRank == 0) {
	 diagnostic << globalMin[i] << "\t" <<
	 globalMax[i] << "\t" <<
	 globalSum[i+1] << "\t" <<
	 globalAvg[i] << "\t";
      }
   }
   if (myRank == 0) diagnostic << endl << write;
   return true;
}
