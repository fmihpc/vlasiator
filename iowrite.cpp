#include <cstdlib>
#include <iostream>
#include <iomanip> // for setprecision()
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>
#include "iowrite.h"
#include "phiprof.hpp"
#include "parameters.h"
#include "logger.h"
#include "vlsvwriter2.h"
#include "vlsvreader2.h"
#include "vlasovmover.h"

using namespace std;
using namespace phiprof;

extern Logger logFile, diagnostic;

typedef Parameters P;


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

bool writeDataReducer(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                      const vector<uint64_t>& cells,
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

template <typename T>
bool writeScalarParameter(string name,T value,VLSVWriter& vlsvWriter,int masterRank,MPI_Comm comm){
   int myRank;
   MPI_Comm_rank(comm,&myRank);
   if(myRank==masterRank){
      map<string,string> attribs;
      std::ostringstream s;
      s << value;
      attribs["value"]=s.str();
      vlsvWriter.writeArrayMaster("PARAMETERS",name,attribs,1,1,&value);
   }
   return true;
}



/*write variables in grid needed by both restart and normal files*/


bool writeCommonGridData(
   VLSVWriter& vlsvWriter,
   const dccrg::Dccrg<SpatialCell>& mpiGrid,
   vector<uint64_t> &cells,
   const uint& index,
   MPI_Comm comm
) {
   int myRank;
   MPI_Comm_rank(comm,&myRank);
   map<string,string> attribs;
   attribs.clear();
   if (vlsvWriter.writeArray("MESH","SpatialGrid",attribs,cells.size(),1,cells.data()) == false) {
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
   attribs.clear();
   if (vlsvWriter.writeArray("COORDS","SpatialGrid",attribs,cells.size(),6,buffer) == false) {
      cerr << "Proc #" << myRank << " failed to write cell coords!" << endl;
   }
   delete[] buffer;   
   writeScalarParameter("t",P::t,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("dt",P::dt,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("tstep",P::tstep,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("fileIndex",index,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("xmin",P::xmin,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("xmax",P::xmax,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("ymin",P::ymin,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("ymax",P::ymax,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("zmin",P::zmin,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("zmax",P::zmax,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("xcells_ini",P::xcells_ini,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("ycells_ini",P::ycells_ini,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("zcells_ini",P::zcells_ini,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("vxmin",P::vxmin,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("vxmax",P::vxmax,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("vymin",P::vymin,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("vymax",P::vymax,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("vzmin",P::vzmin,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("vzmax",P::vzmax,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("vxblocks_ini",P::vxblocks_ini,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("vyblocks_ini",P::vyblocks_ini,vlsvWriter,0,MPI_COMM_WORLD);
   writeScalarParameter("vzblocks_ini",P::vzblocks_ini,vlsvWriter,0,MPI_COMM_WORLD);
   return true; //to make compiler happy,no real errorchecking done
}

bool writeVelocityDistributionData(
   VLSVWriter& vlsvWriter,
   const dccrg::Dccrg<SpatialCell>& mpiGrid,
   vector<uint64_t> &cells,
   MPI_Comm comm) {
    // Write velocity blocks and related data. 
   // In restart we just write velocity grids for all cells.
   //   First write global Ids of those cells which write velocity blocks (here: all cells):
   map<string,string> attribs;
   bool success=true;

   
   //Compute totalBlocks
   uint64_t totalBlocks = 0;  
   vector<uint> blocksPerCell;   
   for(size_t cell=0;cell<cells.size();++cell){
      totalBlocks+=mpiGrid[cells[cell]]->number_of_blocks;
      blocksPerCell.push_back(mpiGrid[cells[cell]]->number_of_blocks);
   }
   
   if (vlsvWriter.writeArray("CELLSWITHBLOCKS","SpatialGrid",attribs,cells.size(),1,cells.data()) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;
   //Write blocks per cell, this has to be in the same order as cellswitblocks so that extracting works
   if(vlsvWriter.writeArray("BLOCKSPERCELL","SpatialGrid",attribs,blocksPerCell.size(),1,blocksPerCell.data()) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;

   //Write velocity block coordinates.
   std::vector<Real> velocityBlockParameters;
   try {
      velocityBlockParameters.reserve(totalBlocks*BlockParams::N_VELOCITY_BLOCK_PARAMS);
      
      //gather data for writing
      for (size_t cell=0; cell<cells.size(); ++cell) {
         SpatialCell* SC = mpiGrid[cells[cell]];
         for (unsigned int block_i=0;block_i < SC->number_of_blocks;block_i++){
            unsigned int block = SC->velocity_block_list[block_i];
            Velocity_Block* block_data = SC->at(block);
            for(unsigned int p=0;p<BlockParams::N_VELOCITY_BLOCK_PARAMS;++p){
               velocityBlockParameters.push_back(block_data->parameters[p]);
            }
         }
      }
      }
   catch (...) {
      success=false;
      }
   
   if( globalSuccess(success,"(MAIN) writeGrid: ERROR: Failed to fill temporary array velocityBlockParameters",MPI_COMM_WORLD) == false) {
      vlsvWriter.close();
      return false;
   }
   
   
   attribs.clear();
   if (vlsvWriter.writeArray("BLOCKCOORDINATES","SpatialGrid",attribs,totalBlocks,BlockParams::N_VELOCITY_BLOCK_PARAMS,velocityBlockParameters.data()) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write BLOCKCOORDINATES to file!" << endl << writeVerbose;
   velocityBlockParameters.clear();
   
   
   // Write values of distribution function:
   std::vector<Real> velocityBlockData;
   try {
      velocityBlockData.reserve(totalBlocks*SIZE_VELBLOCK);
      for (size_t cell=0; cell<cells.size(); ++cell) {
         SpatialCell* SC = mpiGrid[cells[cell]];
         for (unsigned int block_i=0;block_i < SC->number_of_blocks;block_i++){
            unsigned int block = SC->velocity_block_list[block_i];
            Velocity_Block* block_data = SC->at(block);
            for(unsigned int vc=0;vc<SIZE_VELBLOCK;++vc){
               velocityBlockData.push_back(block_data->data[vc]);
            }
         }
         }
   }
   catch (...) {
      success=false;
   }
   
   if( globalSuccess(success,"(MAIN) writeGrid: ERROR: Failed to fill temporary velocityBlockData array",MPI_COMM_WORLD) == false) {
      vlsvWriter.close();
      return false;
   }

   
   attribs.clear();
   attribs["mesh"] = "SpatialGrid";
   if (vlsvWriter.writeArray("BLOCKVARIABLE","avgs",attribs,totalBlocks,SIZE_VELBLOCK,velocityBlockData.data()) == false) success=false;
   if (success ==false)      logFile << "(MAIN) writeGrid: ERROR occurred when writing BLOCKVARIABLE f" << endl << writeVerbose;
   velocityBlockData.clear();
    
   return success;
}

bool writeGrid(
   const dccrg::Dccrg<SpatialCell>& mpiGrid,
   DataReducer& dataReducer,
   const uint& index
) {
   double allStart = MPI_Wtime();
   bool success = true;
   int myRank;

   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   phiprof::start("writeGrid-reduced");
    
   // Create a name for the output file and open it with VLSVWriter:
   stringstream fname;
   fname << P::systemWriteName[index] <<".";
   fname.width(7);
   fname.fill('0');
   fname << P::systemWrites[index] << ".vlsv";
   
   VLSVWriter vlsvWriter;
   vlsvWriter.open(fname.str(),MPI_COMM_WORLD,0);
   
   // Get all local cell Ids 
   map<string,string> attribs;
   vector<uint64_t> cells = mpiGrid.get_cells();
   //no order assumed so let's order cells here
   std::sort(cells.begin(),cells.end());

   //write basic description of grid
   writeCommonGridData(vlsvWriter,mpiGrid,cells,P::systemWrites[index],MPI_COMM_WORLD);
   
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
      lineX =  cells[i] % P::xcells_ini;
      lineY = (cells[i] / P::xcells_ini) % P::ycells_ini;
      lineZ = (cells[i] /(P::xcells_ini *  P::ycells_ini)) % P::zcells_ini;
      // Check that indices are in correct intersection at least in one plane
      if ((P::systemWriteDistributionWriteXlineStride[index] > 0 &&
           P::systemWriteDistributionWriteYlineStride[index] > 0 &&
           lineX % P::systemWriteDistributionWriteXlineStride[index] == 0 &&
           lineY % P::systemWriteDistributionWriteYlineStride[index] == 0)
          ||
          (P::systemWriteDistributionWriteYlineStride[index] > 0 &&
           P::systemWriteDistributionWriteZlineStride[index] > 0 &&
           lineY % P::systemWriteDistributionWriteYlineStride[index] == 0 &&
           lineZ % P::systemWriteDistributionWriteZlineStride[index] == 0)
          ||
          (P::systemWriteDistributionWriteZlineStride[index] > 0 &&
           P::systemWriteDistributionWriteXlineStride[index] > 0 &&
           lineZ % P::systemWriteDistributionWriteZlineStride[index] == 0 &&
           lineX % P::systemWriteDistributionWriteXlineStride[index] == 0)
      ) {
         velSpaceCells.push_back(cells[i]);
         mpiGrid[cells[i]]->parameters[CellParams::ISCELLSAVINGF] = 1.0;
      }
   }
   
   // Write variables calculate d by DataReductionOperators (DRO). We do not know how many 
   // numbers each DRO calculates, so a buffer has to be re-allocated for each DRO:
   for (uint i = 0; i < dataReducer.size(); ++i) {
      writeDataReducer(mpiGrid, cells, dataReducer, i, vlsvWriter);
   }

   uint64_t numVelSpaceCells;
   uint64_t localNumVelSpaceCells;
   localNumVelSpaceCells=velSpaceCells.size();
   MPI_Allreduce(&localNumVelSpaceCells,&numVelSpaceCells,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
   if(numVelSpaceCells>0) {
      //write out velocity space data, if there are cells with this data
      writeVelocityDistributionData(vlsvWriter,mpiGrid,velSpaceCells,MPI_COMM_WORLD);
   }

   
   phiprof::initializeTimer("Barrier","MPI","Barrier");
   phiprof::start("Barrier");
   MPI_Barrier(MPI_COMM_WORLD);
   phiprof::stop("Barrier");
   vlsvWriter.close();
   phiprof::stop("writeGrid-reduced");

   return success;
}
   



bool writeRestart(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                  DataReducer& dataReducer,
                  const string& name,
                  const uint& index) {
   double allStart = MPI_Wtime();
   bool success = true;
   int myRank;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   phiprof::start("writeGrid-restart");

    
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
   //no order assumed so let's order cells here
   std::sort(cells.begin(),cells.end());


   //write basic description of grid
   writeCommonGridData(vlsvWriter,mpiGrid,cells,index,MPI_COMM_WORLD);
   
 
   //write out DROs we need for restarts
   DataReducer restartReducer;
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("background_B",CellParams::BGBX,3));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("perturbed_B",CellParams::PERBX,3));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("moments",CellParams::RHO,4));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("moments_dt2",CellParams::RHO_DT2,4));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("moments_r",CellParams::RHO_R,4));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("moments_v",CellParams::RHO_V,4));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("LB_weight",CellParams::LBWEIGHTCOUNTER,1));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("max_v_dt",CellParams::MAXVDT,1));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("max_r_dt",CellParams::MAXRDT,1));
   restartReducer.addOperator(new DRO::DataReductionOperatorCellParams("max_fields_dt",CellParams::MAXFDT,1));
   
   restartReducer.addOperator(new DRO::MPIrank);
   restartReducer.addOperator(new DRO::BoundaryType);
   restartReducer.addOperator(new DRO::BoundaryLayer);
   restartReducer.addOperator(new DRO::VelocitySubSteps);
   
   for (uint i=0; i<restartReducer.size(); ++i) {
      writeDataReducer(mpiGrid,cells,restartReducer,i,vlsvWriter);
   }
   
   writeVelocityDistributionData(vlsvWriter,mpiGrid,cells,MPI_COMM_WORLD);
   
   vlsvWriter.close();
   
   phiprof::stop("writeGrid-restart");//,1.0e-6*bytesWritten,"MB");

   return success;
}
   



bool writeDiagnostic(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                     DataReducer& dataReducer)
{
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



