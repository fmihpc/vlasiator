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



using namespace std;
using namespace phiprof;

extern Logger logFile, diagnostic;



bool readGrid(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
              const std::string& name){
   
   vector<uint64_t> cellIDs;
   double allStart = MPI_Wtime();
   bool success=true;
   int myRank,processes;

   
   // Attempt to open VLSV file for reading:
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   MPI_Comm_size(MPI_COMM_WORLD,&processes);
   

   phiprof::start("readGrid");
   VLSVParReader file;
   MPI_Info mpiInfo = MPI_INFO_NULL;

   if (file.open(name,MPI_COMM_WORLD,0,mpiInfo) == false) {
      logFile << "(RESTART) VLSVParReader failed to open restart file '" << name << "' for reading!" << endl << write;
      success=false;
   }

   // Get info on array containing cell IDs:
   uint64_t arraySize;
   uint64_t vectorSize;
   VLSV::datatype dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name","SpatialGrid"));
   if (file.getArrayInfo("MESH",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      logFile << "(RESTARTBUILDER) ERROR: Failed to read cell ID array info!" << endl << write;
      return false;
   }
   
   // Read cell IDs:
   char* IDbuffer = new char[arraySize*vectorSize*byteSize];
   if (file.readArrayMaster("MESH",attribs,0,arraySize,IDbuffer) == false) {
      logFile << "(RESTARTBUILDER) ERROR: Failed to read cell IDs!" << endl << write;
      success = false;
   }
   
   // Convert global IDs into VirtualCell::ID type:
   cellIDs.resize(arraySize);
   int N_cells = arraySize;
   if (dataType == VLSV::UINT && byteSize == 4) {
      uint32_t* ptr = reinterpret_cast<uint32_t*>(IDbuffer);
      for (uint64_t i=0; i<arraySize; ++i) cellIDs[i] = ptr[i];
   } else if (dataType == VLSV::UINT && byteSize == 8) {
      uint64_t* ptr = reinterpret_cast<uint64_t*>(IDbuffer);
      for (uint64_t i=0; i<arraySize; ++i) cellIDs[i] = ptr[i];
   } else {
      logFile << "(RESTARTBUILDER) ERROR: VLSVParReader returned an unsupported datatype for cell IDs!" << endl << write;
      success = false;
   }



   


   
   if(myRank==0)
      for(uint i=0;i<cellIDs.size();i++)
         cout << i <<" " <<cellIDs[i]<<endl;
   
   


      
   phiprof::stop("readGrid");
   

   return success;

}


bool writeGrid(const dccrg::Dccrg<SpatialCell>& mpiGrid,
               DataReducer& dataReducer,
               const string& name,
               const uint& index,
               const bool& writeRestart) {
    double allStart = MPI_Wtime();
    bool success = true;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
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
   
   // Get all local cell IDs and write to file:
   map<string,string> attribs;
   vector<uint64_t> cells = mpiGrid.get_cells();

   if (vlsvWriter.writeArray("MESH","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) {
      cerr << "Proc #" << myrank << " failed to write cell IDs!" << endl;
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
      cerr << "Proc #" << myrank << " failed to write cell coords!" << endl;
   }
   delete[] buffer;

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
      if (success == false) logFile << "(MAIN) writeGrid: ERROR datareductionoperator '" << dataReducer.getName(i) << "' returned false!" << endl << writeVerbose;
      
      // Write reduced data to file:
      if (vlsvWriter.writeArray("VARIABLE",variableName,attribs,cells.size(),vectorSize,dataType,dataSize,varBuffer) == false) success = false;
      if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write datareductionoperator data to file!" << endl << writeVerbose;
      delete[] varBuffer;
      varBuffer = NULL;
   }

   // If restart data is not written, exit here:
   if (writeRestart == false) {
      phiprof::initializeTimer("Barrier","MPI","Barrier");
      phiprof::start("Barrier");
      MPI_Barrier(MPI_COMM_WORLD);
      phiprof::stop("Barrier");
      vlsvWriter.close();
      phiprof::stop("writeGrid-reduced");
      return success;
   }

   attribs.clear();
//START TO WRITE RESTART
   
   // Write spatial cell parameters:
   Real* paramsBuffer = new Real[cells.size()*CellParams::N_SPATIAL_CELL_PARAMS];

   for (size_t i = 0; i < cells.size(); ++i)
      for (uint j = 0; j < CellParams::N_SPATIAL_CELL_PARAMS; ++j) {
      paramsBuffer[i*CellParams::N_SPATIAL_CELL_PARAMS+j] = mpiGrid[cells[i]]->parameters[j];
   }
   
   if (vlsvWriter.writeArray("CELLPARAMS","SpatialGrid",attribs,cells.size(),CellParams::N_SPATIAL_CELL_PARAMS,paramsBuffer) == false) {
      logFile << "(MAIN) writeGrid: ERROR failed to write spatial cell parameters!" << endl << writeVerbose;
      success = false;
   }
   delete[] paramsBuffer;
   
   // Write velocity blocks and related data. Which cells write velocity grids 
   // should be requested from a function, but for now we just write velocity grids for all cells.
   // First write global IDs of those cells which write velocity blocks (here: all cells):
   if (vlsvWriter.writeArray("CELLSWITHBLOCKS","SpatialGrid",attribs,cells.size(),1,&(cells[0])) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;
   
   // Write the number of velocity blocks in each spatial cell. Again a temporary buffer is used:
   uint* N_blocks = new uint[cells.size()];
   uint64_t totalBlocks = 0;
   for (size_t cell=0; cell<cells.size(); ++cell) {
      N_blocks[cell] = mpiGrid[cells[cell]]->size();
      totalBlocks += mpiGrid[cells[cell]]->size();
   }
   if (vlsvWriter.writeArray("NBLOCKS","SpatialGrid",attribs,cells.size(),1,N_blocks) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write NBLOCKS to file!" << endl << writeVerbose;
   delete[] N_blocks;
   
   double start = MPI_Wtime();

   // Write velocity block coordinates.// TODO: add support for MPI_Datatype in startMultiwrite... or use normal writeArray as all data is collected already in one place

   std::vector<Real> velocityBlockParameters;
   velocityBlockParameters.reserve(totalBlocks*BlockParams::N_VELOCITY_BLOCK_PARAMS);

   // gather data for writing
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

   attribs["mesh"] = "SpatialGrid";
   if (vlsvWriter.writeArray("BLOCKVARIABLE","avgs",attribs,totalBlocks,SIZE_VELBLOCK,&(velocityBlockData[0])) == false) success=false;
   if (success ==false)      logFile << "(MAIN) writeGrid: ERROR occurred when writing BLOCKVARIABLE avgs" << endl << writeVerbose;
   velocityBlockData.clear();
   
   double end = MPI_Wtime();


   vlsvWriter.close();

   double allEnd = MPI_Wtime();
   
   //double bytesWritten = totalBlocks*((SIZE_VELBLOCK+SIZE_BLOCKPARAMS)*sizeof(Real)+(SIZE_NBRS_VEL)*sizeof(uint));
   double secs = end-start;
   //FIXME, should be global info
   //logFile << "Wrote " << bytesWritten/1.0e6 << " MB of data in " << secs << " seconds, datarate is " << bytesWritten/secs/1.0e9 << " GB/s" << endl << writeVerbose;
   
   double allSecs = allEnd-allStart;

   phiprof::stop("writeGrid-restart");//,1.0e-6*bytesWritten,"MB");
   return success;
}











bool writeDiagnostic(const dccrg::Dccrg<SpatialCell>& mpiGrid,
		       DataReducer& dataReducer,
		       luint tstep)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   
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
   
   if (printDiagnosticHeader == true && myrank == 0) {
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
      if (myrank == 0) {
	 diagnostic << globalMin[i] << "\t" <<
	 globalMax[i] << "\t" <<
	 globalSum[i+1] << "\t" <<
	 globalAvg[i] << "\t";
      }
   }
   if (myrank == 0) diagnostic << endl << write;
   return true;
}
