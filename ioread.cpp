/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip> // for setprecision()
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>
#include <array>
#include <sys/types.h>
#include <sys/stat.h>

#include "ioread.h"
#include "phiprof.hpp"
#include "parameters.h"
#include "logger.h"
#include "vlsv_reader_parallel.h"
#include "vlasovmover.h"
#include "object_wrapper.h"
#include "grid.h"

using namespace std;
using namespace phiprof;

extern Logger logFile, diagnostic;

typedef Parameters P;

/*!
 * \brief Checks for command files written to the local directory.
 * If a file STOP was written and is readable, then a bailout with restart writing is initiated.
 * If a file KILL was written and is readable, then a bailout without a restart is initiated.
 * If a file SAVE was written and is readable, then restart writing without a bailout is initiated.
 * If a file DOLB was written and is readable, then a new load balancing is initiated.
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
   if(stat("DOLB", &tempStat) == 0) {
      cerr << "Received an external DOLB command. Balancing load." << endl;
      globalflags::balanceLoad = true;
      char newName[80];
      // Get the current time.
      const time_t rawTime = time(NULL);
      const struct tm * timeInfo = localtime(&rawTime);
      strftime(newName, 80, "DOLB_%F_%H-%M-%S", timeInfo);
      rename("DOLB", newName);
      return;
   }
   if(stat("DOMR", &tempStat) == 0) {
      cerr << "Received an external DOMR command. Refining grid." << endl;
      globalflags::doRefine = true;
      char newName[80];
      // Get the current time.
      const time_t rawTime = time(NULL);
      const struct tm * timeInfo = localtime(&rawTime);
      strftime(newName, 80, "DOMR_%F_%H-%M-%S", timeInfo);
      rename("DOMR", newName);
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
bool exitOnError(bool success, const string& message, MPI_Comm comm) {
   int successInt;
   int globalSuccessInt;
   if(success)
      successInt=1;
   else
      successInt=0;

   MPI_Allreduce(&successInt,&globalSuccessInt,1,MPI_INT,MPI_MIN,comm);

   if(globalSuccessInt==1) {
      return true;
   } else {
      logFile << message << endl<<write ;
      exit(1);
   }
}

/*!
 \brief Read cell ID's
 Read in cell ID's from file. Note: Uses the newer version of vlsv parallel reader
 \param file Some vlsv reader with a file open
 \param fileCells Vector in whic to store the cell ids
 \param masterRank The simulation's master rank id (Vlasiator uses 0, which should be the default)
 \param comm MPI comm (MPI_COMM_WORLD should be the default)
*/
bool readCellIds(vlsv::ParallelReader & file, vector<CellID>& fileCells, const int masterRank,MPI_Comm comm)
{
   // Get info on array containing cell Ids:
   uint64_t arraySize = 0;
   uint64_t vectorSize;
   vlsv::datatype::type dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   bool success=true;
   int rank;
   MPI_Comm_rank(comm,&rank);
   if (rank==masterRank) {
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
      const uint64_t& numberOfCells = arraySize;
      fileCells.resize(numberOfCells);
      if (dataType == vlsv::datatype::type::UINT && byteSize == 4) {
         uint32_t* ptr = reinterpret_cast<uint32_t*>(IDbuffer);
         //Input cell ids
         for (uint64_t i=0; i<numberOfCells; ++i) {
            const CellID cellID = ptr[i];
            fileCells[i] = cellID;
         }
      } else if (dataType == vlsv::datatype::type::UINT && byteSize == 8) {
         uint64_t* ptr = reinterpret_cast<uint64_t*>(IDbuffer);
         for (uint64_t i=0; i<numberOfCells; ++i) {
            const CellID cellID = ptr[i];
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

/* Read the total number of velocity blocks per spatial cell in the spatial mesh.
 * The returned value for each cell is a sum of the velocity blocks associated in each particle species. 
 * The value is used to calculate an initial load balance after restart.
 * @param file Some vlsv reader with a file open (can be old or new vlsv reader)
 * @param nBlocks Vector for holding information on cells and the number of blocks in them -- this function saves data here
 * @param masterRank The master rank of this process (Vlasiator uses masterRank = 0 and so it should be the default)
 * @param comm MPI comm
 * @return Returns true if the operation was successful
 @ @see exec_readGrid
*/
bool readNBlocks(vlsv::ParallelReader& file,const std::string& meshName,
                 std::vector<size_t>& nBlocks,int masterRank,MPI_Comm comm) {
   bool success = true;

   // Get info on array containing cell IDs:
   uint64_t arraySize;
   uint64_t vectorSize;
   vlsv::datatype::type dataType;
   uint64_t byteSize;

   // Read mesh bounding box to all processes, the info in bbox contains 
   // the number of spatial cells in the mesh.
   // (This is *not* the physical coordinate bounding box.)
   list<pair<string,string> > attribsIn;
   map<string,string> attribsOut;
   attribsIn.push_back(make_pair("mesh",meshName));

   // Read number of domains and domain sizes
   uint64_t N_domains;
   file.getArrayAttributes("MESH_DOMAIN_SIZES",attribsIn,attribsOut);
   auto it = attribsOut.find("arraysize");
   if (it == attribsOut.end()) {
      cerr << "VLSV\t\t ERROR: Array 'MESH_DOMAIN_SIZES' XML tag does not have attribute 'arraysize'" << endl;
      return false;
   } else {
      N_domains = atoi(it->second.c_str());
   }

   uint64_t N_spatialCells = 0;

	int64_t* domainInfo = NULL;
	if (file.read("MESH_DOMAIN_SIZES",attribsIn,0,N_domains,domainInfo) == false) return false;

	for (uint i_domain = 0; i_domain < N_domains; ++i_domain) {
		
		N_spatialCells += domainInfo[2*i_domain];

	}

   nBlocks.resize(N_spatialCells);


   #pragma omp parallel for
   for (size_t i=0; i<nBlocks.size(); ++i) nBlocks[i] = 0;

   // Note: the input file contains N particle species, which are also 
   // defined in the configuration file. We need to read the BLOCKSPERCELL 
   // array for each species, and sum the values for each spatial cell.
   set<string> speciesNames;
   if (file.getUniqueAttributeValues("BLOCKSPERCELL","name",speciesNames) == false) return false;

   // Iterate over all particle species and read in BLOCKSPERCELL array 
   // to all processes, and add the values to nBlocks
   uint64_t* buffer = new uint64_t[N_spatialCells];
   for (set<string>::const_iterator s=speciesNames.begin(); s!=speciesNames.end(); ++s) {      
      attribsIn.clear();
      attribsIn.push_back(make_pair("mesh",meshName));
      attribsIn.push_back(make_pair("name",*s));
      if (file.getArrayInfo("BLOCKSPERCELL",attribsIn,arraySize,vectorSize,dataType,byteSize) == false) return false;

      if (file.read("BLOCKSPERCELL",attribsIn,0,arraySize,buffer) == false) {
         delete [] buffer; buffer = NULL;
         return false;
      }

      #pragma omp parallel for
      for (size_t i=0; i<N_spatialCells; ++i) {
         nBlocks[i] += buffer[i];
      }
   }
   delete [] buffer; buffer = NULL;
   return success;
}

/*! A function for reading parameters, e.g., 'timestep'.
 \param file VLSV parallel reader with a file open.
 \param name Name of the parameter.
 \param value Variable in which to store the scalar variable (double, float, int .. ).
 \param masterRank The master process' id (Vlasiator uses 0 so this should equal 0 by default).
 \param comm MPI comm (MPI_COMM_WORLD should be the default).
 \return Returns true if the operation is successful. */
template <typename T>
bool readScalarParameter(vlsv::ParallelReader& file,string name,T& value,int masterRank,MPI_Comm comm) {
   if (file.readParameter(name,value) == false) {
      logFile << "(RESTART) ERROR: Failed to read parameter '" << name << "' value in ";
      logFile << __FILE__ << ":" << __LINE__ << endl << write;
      return false;
   }
   return true;
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
bool checkScalarParameter(vlsv::ParallelReader& file,const string& name,T correctValue,int masterRank,MPI_Comm comm) {
   T value;
   if (readScalarParameter(file,name,value,masterRank,comm) == false) {
      ostringstream s;
      s << "(RESTART) ERROR: Failed to read parameter '" << name << "' value in " << __FILE__ << ":" << __LINE__ << endl;
      exitOnError(false, s.str(), MPI_COMM_WORLD);
      return false;
   }
   if (value != correctValue){
      ostringstream s;
      s << "(RESTART) Parameter " << name << " has mismatching value.";
      s << " CFG value = " << correctValue;
      s << " Restart file value = " << value;
      exitOnError(false,s.str(),MPI_COMM_WORLD);
      return false;
   } else {
      return true;
   }
}

/** Read velocity block mesh data and distribution function data belonging to this process 
 * for the given particle species. This function must be called simultaneously by all processes.
 * @param file VLSV reader with input file open.
 * @param spatMeshName Name of the spatial mesh.
 * @param fileCells List of all spatial cell IDs.
 * @param localCellStartOffset The offset from which to start reading cells.
 * @param localCells How many spatial cells after the offset to read.
 * @param blocksPerCell Number of velocity blocks for this particle species in each spatial cell belonging to this process.
 * @param localBlockStartOffset Offset into velocity block data arrays from which to start reading data.
 * @param localBlocks Number of velocity blocks for this species assigned to this process.
 * @param mpiGrid Parallel grid library.
 * @param popID ID of the particle species who's data is to be read.
 * @return If true, velocity block data was read successfully.*/
template <typename fileReal>
bool _readBlockData(
   vlsv::ParallelReader & file,
   const std::string& spatMeshName,
   const std::vector<uint64_t>& fileCells,
   const uint64_t localCellStartOffset,
   const uint64_t localCells,
   const vmesh::LocalID* blocksPerCell,
   const uint64_t localBlockStartOffset,
   const uint64_t localBlocks,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   std::function<vmesh::GlobalID(vmesh::GlobalID)> blockIDremapper,
   const uint popID
) {   
   uint64_t arraySize;
   uint64_t avgVectorSize;
   vlsv::datatype::type dataType;
   uint64_t byteSize;
   list<pair<string,string> > avgAttribs;
   bool success=true;
   const string popName = getObjectWrapper().particleSpecies[popID].name;
   const string tagName = "BLOCKIDS";
   
   avgAttribs.push_back(make_pair("mesh",spatMeshName));
   avgAttribs.push_back(make_pair("name",popName));
   
    //Get block id array info and store them into blockIdAttribs, lockIdByteSize, blockIdDataType, blockIdVectorSize
  list<pair<string,string> > blockIdAttribs;
  uint64_t blockIdVectorSize, blockIdByteSize;
  vlsv::datatype::type blockIdDataType;
  blockIdAttribs.push_back( make_pair("mesh", spatMeshName));
  blockIdAttribs.push_back( make_pair("name", popName));
  if (file.getArrayInfo("BLOCKIDS",blockIdAttribs,arraySize,blockIdVectorSize,blockIdDataType,blockIdByteSize) == false ){
    logFile << "(RESTART) ERROR: Failed to read BLOCKCOORDINATES array info " << endl << write;
    return false;
  }
  if(file.getArrayInfo("BLOCKVARIABLE",avgAttribs,arraySize,avgVectorSize,dataType,byteSize) == false ){
    logFile << "(RESTART) ERROR: Failed to read BLOCKVARIABLE array info " << endl << write;
    return false;
  }

   //Some routine error checks:
   if( avgVectorSize!=WID3 ){
      logFile << "(RESTART) ERROR: Blocksize does not match in restart file " << endl << write;
      return false;
   }
   if( byteSize != sizeof(fileReal) ) {
      logFile << "(RESTART) ERROR: Bad avgs bytesize at " << __FILE__ << " " << __LINE__ << endl << write;
      return false;
   }
   
   if( blockIdByteSize != sizeof(vmesh::GlobalID)) {
      logFile << "(RESTART) ERROR: BlockID data size does not match " << __FILE__ << " " << __LINE__ << endl << write;
      return false;
   }

   fileReal* avgBuffer = new fileReal[avgVectorSize * localBlocks]; //avgs data for all cells
   vmesh::GlobalID * blockIdBuffer = new vmesh::GlobalID[blockIdVectorSize * localBlocks]; //blockids of all cells

   //Read block ids and data
   if (file.readArray("BLOCKIDS", blockIdAttribs, localBlockStartOffset, localBlocks, (char*)blockIdBuffer ) == false) {
      cerr << "ERROR, failed to read BLOCKIDS in " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   if (file.readArray("BLOCKVARIABLE", avgAttribs, localBlockStartOffset, localBlocks, (char*)avgBuffer) == false) {
      cerr << "ERROR, failed to read BLOCKVARIABLE in " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   
   uint64_t blockBufferOffset=0;
   //Go through all spatial cells     
   vector<vmesh::GlobalID> blockIdsInCell; //blockIds in a particular cell, temporary usage
   for(uint64_t i=0; i<localCells; i++) {
      CellID cell = fileCells[localCellStartOffset + i]; //spatial cell id 
      vmesh::LocalID nBlocksInCell = blocksPerCell[i];
      //copy blocks in this cell to vector blockIdsInCell, size of read in data has been checked earlier
      blockIdsInCell.reserve(nBlocksInCell);
      blockIdsInCell.assign(blockIdBuffer + blockBufferOffset, blockIdBuffer + blockBufferOffset + nBlocksInCell);
      for(auto& id : blockIdsInCell) {
         id = blockIDremapper(id);
      }
      mpiGrid[cell]->add_velocity_blocks(blockIdsInCell,popID); //allocate space for all blocks and create them
      //copy avgs data, here a conversion may happen between float and double
      Realf *cellBlockData=mpiGrid[cell]->get_data(popID);
      for(uint64_t i = 0; i< WID3 * nBlocksInCell ; i++){
         cellBlockData[i] =  avgBuffer[blockBufferOffset*WID3 + i];
      }
      blockBufferOffset += nBlocksInCell; //jump to location of next local cell
   }

   delete[] avgBuffer;
   delete[] blockIdBuffer;
   return success;
}

/** Read velocity block data of all existing particle species.
 * @param file VLSV reader.
 * @param meshName Name of the spatial mesh.
 * @param fileCells Vector containing spatial cell IDs.
 * @param localCellStartOffset Offset into fileCells, determines where the cells belonging 
 * to this process start.
 * @param localCells Number of spatial cells assigned to this process.
 * @param mpiGrid Parallel grid library.
 * @return If true, velocity block data was read successfully.*/
bool readBlockData(
        vlsv::ParallelReader& file,
        const string& meshName,
        const vector<CellID>& fileCells,
        const uint64_t localCellStartOffset,
        const uint64_t localCells,
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
   ) {
   bool success = true;

   const uint64_t bytesReadStart = file.getBytesRead();
   int N_processes;
   MPI_Comm_size(MPI_COMM_WORLD,&N_processes);

   uint64_t arraySize;
   uint64_t vectorSize;
   vlsv::datatype::type dataType;
   uint64_t byteSize;
   uint64_t* offsetArray = new uint64_t[N_processes];

   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      const string& popName = getObjectWrapper().particleSpecies[popID].name;

      // Create a cellID remapping lambda that can renumber our velocity space, should it's size have changed.
      // By default, this is a no-op that keeps the blockIDs untouched.
      std::function<vmesh::GlobalID(vmesh::GlobalID)> blockIDremapper = [](vmesh::GlobalID oldID) -> vmesh::GlobalID {return oldID;};

      // Check that velocity space extents and DV matches the grids we have created
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("mesh",popName));
      std::array<unsigned int, 6> fileMeshBBox;
      unsigned int* bufferpointer = &fileMeshBBox[0];
      if (file.read("MESH_BBOX",attribs,0,6,bufferpointer,false) == false) {
         logFile << "(RESTART) ERROR: Failed to read MESH_BBOX at " << __FILE__ << ":" << __LINE__ << endl << write;
         success = false;
      }

      const size_t meshID = getObjectWrapper().particleSpecies[popID].velocityMesh;
      const vmesh::MeshParameters& ourMeshParams = getObjectWrapper().velocityMeshes[meshID];
      if(fileMeshBBox[0] != ourMeshParams.gridLength[0] ||
            fileMeshBBox[1] != ourMeshParams.gridLength[1] ||
            fileMeshBBox[2] != ourMeshParams.gridLength[2]) {

         logFile << "(RESTART) INFO: velocity mesh sizes don't match:" << endl
                 << "    restart file has " << fileMeshBBox[0] << " x " << fileMeshBBox[1] << " x " << fileMeshBBox[2] << "," << endl
                 << "    config specifies " << ourMeshParams.gridLength[0] << " x " <<  ourMeshParams.gridLength[1] << " x " <<  ourMeshParams.gridLength[2] << endl << write;

         if(ourMeshParams.gridLength[0] < fileMeshBBox[0] ||
               ourMeshParams.gridLength[1] < fileMeshBBox[1] ||
               ourMeshParams.gridLength[2] < fileMeshBBox[2]) {
            logFile << "(RESTART) ERROR: trying to shrink velocity space." << endl << write;
            abort();
         }

         // If we are mismatched, we have to iterate through the velocity coords to see if we have a
         // chance at renumbering.
         std::vector<Real> fileVelCoordsX(fileMeshBBox[0]*fileMeshBBox[3]+1);
         std::vector<Real> fileVelCoordsY(fileMeshBBox[1]*fileMeshBBox[4]+1);
         std::vector<Real> fileVelCoordsZ(fileMeshBBox[2]*fileMeshBBox[5]+1);

         Real* tempPointer = fileVelCoordsX.data();
         if (file.read("MESH_NODE_CRDS_X",attribs,0,fileMeshBBox[0]*fileMeshBBox[3]+1,tempPointer,false) == false) {
            logFile << "(RESTART) ERROR: Failed to read MESH_NODE_CRDS_X at " << __FILE__ << ":" << __LINE__ << endl << write;
            success = false;
         }
         tempPointer = fileVelCoordsY.data();
         if (file.read("MESH_NODE_CRDS_Y",attribs,0,fileMeshBBox[1]*fileMeshBBox[4]+1,tempPointer,false) == false) {
            logFile << "(RESTART) ERROR: Failed to read MESH_NODE_CRDS_X at " << __FILE__ << ":" << __LINE__ << endl << write;
            success = false;
         }
         tempPointer = fileVelCoordsZ.data();
         if (file.read("MESH_NODE_CRDS_Z",attribs,0,fileMeshBBox[2]*fileMeshBBox[5]+1,tempPointer,false) == false) {
            logFile << "(RESTART) ERROR: Failed to read MESH_NODE_CRDS_X at " << __FILE__ << ":" << __LINE__ << endl << write;
            success = false;
         }

         const Real dVx = getObjectWrapper().velocityMeshes[meshID].cellSize[0];
         for(const auto& c : fileVelCoordsX) {
            Real cellindex = (c - getObjectWrapper().velocityMeshes[meshID].meshMinLimits[0]) / dVx;
            if(fabs(nearbyint(cellindex) - cellindex) > 1./10000.) {
               logFile << "(RESTART) ERROR: Can't resize velocity space as cell coordinates don't match." << endl
                  << "          (X coordinate " << c << " = " << cellindex <<" * " << dVx << " + " << getObjectWrapper().velocityMeshes[meshID].meshMinLimits[0] << endl
                  << "           coordinate  = cellindex *   dV  +  meshMinLimits)" << endl << write;
               abort();
            }
         }

         const Real dVy = getObjectWrapper().velocityMeshes[meshID].cellSize[1];
         for(const auto& c : fileVelCoordsY) {
            Real cellindex = (c - getObjectWrapper().velocityMeshes[meshID].meshMinLimits[1]) / dVy;
            if(fabs(nearbyint(cellindex) - cellindex) > 1./10000.) {
               logFile << "(RESTART) ERROR: Can't resize velocity space as cell coordinates don't match." << endl
                  << "           (Y coordinate " << c << " = " << cellindex <<" * " << dVy << " + " << getObjectWrapper().velocityMeshes[meshID].meshMinLimits[1] << endl
                  << "           coordinate  = cellindex *   dV  +  meshMinLimits)" << endl << write;
               abort();
            }
         }

         const Real dVz = getObjectWrapper().velocityMeshes[meshID].cellSize[2];
         for(const auto& c : fileVelCoordsY) {
            Real cellindex = (c - getObjectWrapper().velocityMeshes[meshID].meshMinLimits[2]) / dVz;
            if(fabs(nearbyint(cellindex) - cellindex) > 1./10000.) {
               logFile << "(RESTART) ERROR: Can't resize velocity space as cell coordinates don't match." << endl
                  << "           (Z coordinate " << c << " = " << cellindex <<" * " << dVz << " + " << getObjectWrapper().velocityMeshes[meshID].meshMinLimits[2] << endl
                  << "           coordinate  = cellindex *   dV  +  meshMinLimits)" << endl << write;
               abort();
            }
         }

         // If we haven't aborted above, we can apparently renumber our
         // cellIDs. Build an approprita blockIDremapper lambda for this purpose.
         std::array<int, 3> velGridOffset;
         velGridOffset[0] = (fileVelCoordsX[0] - getObjectWrapper().velocityMeshes[meshID].meshMinLimits[0]) / dVx;
         velGridOffset[1] = (fileVelCoordsY[0] - getObjectWrapper().velocityMeshes[meshID].meshMinLimits[1]) / dVy;
         velGridOffset[2] = (fileVelCoordsZ[0] - getObjectWrapper().velocityMeshes[meshID].meshMinLimits[2]) / dVz;

         if((velGridOffset[0] % ourMeshParams.blockLength[0] != 0) ||
               (velGridOffset[1] % ourMeshParams.blockLength[1] != 0) ||
               (velGridOffset[2] % ourMeshParams.blockLength[2] != 0)) {
            logFile << "(RESTART) ERROR: resizing velocity space on restart must end up with the old velocity space" << endl
                    << "                 at a block boundary of the new space!" << endl
                    << "                 (It now starts at cell [" << velGridOffset[0] << ", " << velGridOffset[1] << "," << velGridOffset[2] << "])" << endl << write;
            abort();
         }

         velGridOffset[0] /= ourMeshParams.blockLength[0];
         velGridOffset[1] /= ourMeshParams.blockLength[1];
         velGridOffset[2] /= ourMeshParams.blockLength[2];

         blockIDremapper = [fileMeshBBox,velGridOffset,ourMeshParams](vmesh::GlobalID oldID) -> vmesh::GlobalID {
            unsigned int x,y,z;
            x = oldID % fileMeshBBox[0];
            y = (oldID / fileMeshBBox[0]) % fileMeshBBox[1];
            z = oldID / (fileMeshBBox[0] * fileMeshBBox[1]);

            x += velGridOffset[0];
            y += velGridOffset[1];
            z += velGridOffset[2];

            //logFile << " Remapping " << oldID << "(" << x << "," << y << "," << z << ") to " << x + y * ourMeshParams.gridLength[0] + z* ourMeshParams.gridLength[0] * ourMeshParams.gridLength[1] << endl << write;
            return x + y * ourMeshParams.gridLength[0] + z* ourMeshParams.gridLength[0] * ourMeshParams.gridLength[1];
         };

         logFile << "    => Resizing velocity space by renumbering GlobalIDs." << endl << endl << write;
      }

      // In restart files each spatial cell has an entry in CELLSWITHBLOCKS. 
      // Each process calculates how many velocity blocks it has for this species.
      attribs.clear();
      attribs.push_back(make_pair("mesh",meshName));
      attribs.push_back(make_pair("name",popName));
      vmesh::LocalID* blocksPerCell = NULL;
      
      if (file.read("BLOCKSPERCELL",attribs,localCellStartOffset,localCells,blocksPerCell,true) == false) {
         logFile << "(RESTART) ERROR: Failed to read BLOCKSPERCELL at " << __FILE__ << ":" << __LINE__ << endl << write;
         success = false;
      }

      // Count how many velocity blocks this process gets
      uint64_t blockSum = 0;
      for (uint64_t i=0; i<localCells; ++i){
         blockSum += blocksPerCell[i];
      }
      
      // Gather all block sums to master process who will them broadcast 
      // the values to everyone
      MPI_Allgather(&blockSum,1,MPI_Type<uint64_t>(),offsetArray,1,MPI_Type<uint64_t>(),MPI_COMM_WORLD);      
      
      // Calculate the offset from which this process starts reading block data
      uint64_t myOffset = 0;
      for (int64_t i=0; i<mpiGrid.get_rank(); ++i) myOffset += offsetArray[i];
      
      if (file.getArrayInfo("BLOCKVARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
         logFile << "(RESTART)  ERROR: Failed to read BLOCKVARIABLE INFO" << endl << write;
         return false;
      }

      // Call _readBlockData
      if (dataType == vlsv::datatype::type::FLOAT) {
         switch (byteSize) {
            case sizeof(double):
               if (_readBlockData<double>(file,meshName,fileCells,localCellStartOffset,localCells,blocksPerCell,
                                          myOffset,blockSum,mpiGrid,blockIDremapper,popID) == false) success = false;
               break;
            case sizeof(float):
               if (_readBlockData<float>(file,meshName,fileCells,localCellStartOffset,localCells,blocksPerCell,
                                         myOffset,blockSum,mpiGrid,blockIDremapper,popID) == false) success = false;
               break;
         }
      } else if (dataType == vlsv::datatype::type::UINT) {
         switch (byteSize) {
            case sizeof(uint32_t):
               if (_readBlockData<uint32_t>(file,meshName,fileCells,localCellStartOffset,localCells,blocksPerCell,
                                            myOffset,blockSum,mpiGrid,blockIDremapper,popID) == false) success = false;
               break;
            case sizeof(uint64_t):
               if (_readBlockData<uint64_t>(file,meshName,fileCells,localCellStartOffset,localCells,blocksPerCell,
                                            myOffset,blockSum,mpiGrid,blockIDremapper,popID) == false) success = false;
               break;
         }
      } else if (dataType == vlsv::datatype::type::INT) {
         switch (byteSize) {
            case sizeof(int32_t):
               if (_readBlockData<int32_t>(file,meshName,fileCells,localCellStartOffset,localCells,blocksPerCell,
                                           myOffset,blockSum,mpiGrid,blockIDremapper,popID) == false) success = false;
               break;
            case sizeof(int64_t):
               if (_readBlockData<int64_t>(file,meshName,fileCells,localCellStartOffset,localCells,blocksPerCell,
                                           myOffset,blockSum,mpiGrid,blockIDremapper,popID) == false) success = false;
               break;
         }
      } else {
         logFile << "(RESTART) ERROR: Failed to read data type at readCellParamsVariable" << endl << write;
         success = false;
      }
      delete [] blocksPerCell; blocksPerCell = NULL;
   } // for-loop over particle species

   delete [] offsetArray; offsetArray = NULL;
   
   const uint64_t bytesReadEnd = file.getBytesRead() - bytesReadStart;
   logFile << "Velocity meshes and data read, approximate data rate is ";
   logFile << vlsv::printDataRate(bytesReadEnd,file.getReadTime()) << endl << write;

   return success;
}

/*! Reads cell parameters from the file and saves them in the right place in mpiGrid
 \param file Some parallel vlsv reader with a file open
 \param fileCells List of all cell ids
 \param localCellStartOffset Offset in the fileCells list for this process ( calculated so that the amount of blocks is distributed somewhat evenly between processes)
 \param localCells The amount of cells to read in this process after localCellStartOffset
 \param cellParamsIndex The parameter of the cell index e.g. CellParams::RHOM
 \param expectedVectorSize The amount of elements in the parameter (parameter can be a scalar or a vector of size N)
 \param mpiGrid Vlasiator's grid (the parameters are saved here)
 \return Returns true if the operation is successful
 */
template <typename fileReal>
static bool _readCellParamsVariable(
                                    vlsv::ParallelReader& file,
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
   vlsv::datatype::type dataType;
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
   if(file.readArray("VARIABLE", attribs, localCellStartOffset, localCells, (char*) buffer) == false ) {
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
 \param cellParamsIndex The parameter of the cell index e.g. CellParams::RHOM
 \param expectedVectorSize The amount of elements in the parameter (parameter can be a scalar or a vector of size N)
 \param mpiGrid Vlasiator's grid (the parameters are saved here)
 \return Returns true if the operation is successful
 */
bool readCellParamsVariable(
   vlsv::ParallelReader& file,
   const vector<CellID>& fileCells,
   const uint64_t localCellStartOffset,
   const uint64_t localCells,
   const string& variableName,
   const size_t cellParamsIndex,
   const size_t expectedVectorSize,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
) {
   uint64_t arraySize;
   uint64_t vectorSize;
   vlsv::datatype::type dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   
   attribs.push_back(make_pair("name",variableName));
   attribs.push_back(make_pair("mesh","SpatialGrid"));
   
   if (file.getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      logFile << "(RESTART)  ERROR: Failed to read " << endl << write;
      return false;
   }

   // Call _readCellParamsVariable
   if( dataType == vlsv::datatype::type::FLOAT ) {
      switch (byteSize) {
         case sizeof(double):
            return _readCellParamsVariable<double>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
         case sizeof(float):
            return _readCellParamsVariable<float>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
      }
   } else if( dataType == vlsv::datatype::type::UINT ) {
      switch (byteSize) {

         case sizeof(uint32_t):
            return _readCellParamsVariable<uint32_t>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
         case sizeof(uint64_t):
            return _readCellParamsVariable<uint64_t>( file, fileCells, localCellStartOffset, localCells, variableName, cellParamsIndex, expectedVectorSize, mpiGrid );
            break;
      }
   } else if( dataType == vlsv::datatype::type::INT ) {
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
   // For compiler purposes
   return false;
}

/*! Read a fsgrid variable (consinting of N real values) from the given vlsv file.
 * \param file VLSV parallel reader with a file open.
 * \param variableName Name of the variable in the file
 * \param numWritingRanks Number of mpi ranks that were used to write this file (used for reconstruction of the spatial order)
 * \param targetGrid target location where the data will be stored.
 */
template<unsigned long int N> bool readFsGridVariable(
   vlsv::ParallelReader& file, const string& variableName, int numWritingRanks, FsGrid<std::array<Real, N>,FS_STENCIL_WIDTH> & targetGrid) {

   phiprof::Timer preparations {"preparations"};

   uint64_t arraySize;
   uint64_t vectorSize;
   vlsv::datatype::type dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   bool convertFloatType = false;
   
   attribs.push_back(make_pair("name",variableName));
   attribs.push_back(make_pair("mesh","fsgrid"));

   phiprof::Timer getArrayInfo {"getArrayInfo"};
   if (file.getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      logFile << "(RESTART)  ERROR: Failed to read " << endl << write;
      return false;
   }
   if(! (dataType == vlsv::datatype::type::FLOAT && byteSize == sizeof(Real))) {
      logFile << "(RESTART) Converting floating point format of fsgrid variable " << variableName << " from " << byteSize * 8 << " bits to " << sizeof(Real) * 8 << " bits." << endl << write;
      convertFloatType = true;
      // Note: this implicitly assumes that Real is of type double, and we either read a double in directly, or read a float and convert it to double.
   }
   getArrayInfo.stop();

   // Are we restarting from the same number of tasks, or a different number?
   int size, myRank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   std::array<FsGridTools::FsIndex_t,3>& localSize = targetGrid.getLocalSize();
   std::array<FsGridTools::FsIndex_t,3>& localStart = targetGrid.getLocalStart();
   std::array<FsGridTools::FsSize_t,3>& globalSize = targetGrid.getGlobalSize();

   // Determine our tasks storage size
   size_t storageSize = localSize[0]*localSize[1]*localSize[2];

   preparations.stop();
   std::array<FsGridTools::Task_t,3> fileDecomposition={0,0,0};
   // No override given (zero array)
   if(P::overrideReadFsGridDecomposition == fileDecomposition){
      // Try and read the decomposition from file
      if (readFsgridDecomposition(file, fileDecomposition) == false) {
         exitOnError(false, "(RESTART) Failed to read Fsgrid decomposition", MPI_COMM_WORLD);
      }
   } else { 
      // Override
      logFile << "(RESTART) Using manual override for FsGrid MESH_DECOMPOSITION." << endl << write;
      fileDecomposition = P::overrideReadFsGridDecomposition;
      int fsgridInputRanks=0;
      // Read numWritingRanks from file, that should exist and be sane
      if(readScalarParameter(file,"numWritingRanks",fsgridInputRanks, MASTER_RANK, MPI_COMM_WORLD) == false) {
         exitOnError(false, "(RESTART) FSGrid writing rank number not found in restart file", MPI_COMM_WORLD);
      }
      // Check that the override is sane wrt. numWritingRanks
      if(fileDecomposition[0]*fileDecomposition[1]*fileDecomposition[2] != fsgridInputRanks){
         exitOnError(false, "(RESTART) Trying to use a manual FsGrid decomposition for a file with a differing number of input ranks.", MPI_COMM_WORLD);
      }
   }

   std::array<FsGridTools::Task_t,3> decomposition = targetGrid.getDecomposition();
   // targetGrid.computeDomainDecomposition(globalSize, size, decomposition);

   if(decomposition == fileDecomposition) {
      // Easy case: same decomposition => slurp it in.
      //

      // Determine offset in file by summing up all the previous tasks' sizes.
      size_t localStartOffset = 0;
      for(FsGridTools::Task_t task = 0; task < myRank; task++) {
         std::array<FsGridTools::FsIndex_t, 3> thatTasksSize;
         thatTasksSize[0] = targetGrid.calcLocalSize(globalSize[0], decomposition[0], task/decomposition[2]/decomposition[1]);
         thatTasksSize[1] = targetGrid.calcLocalSize(globalSize[1], decomposition[1], (task/decomposition[2])%decomposition[1]);
         thatTasksSize[2] = targetGrid.calcLocalSize(globalSize[2], decomposition[2], task%decomposition[2]);
         localStartOffset += thatTasksSize[0] * thatTasksSize[1] * thatTasksSize[2];
      }
      
      // Read into buffer
      std::vector<Real> buffer(storageSize*N);

      if(file.readArray("VARIABLE",attribs, localStartOffset, storageSize, buffer.data()) == false) {
         logFile << "(RESTART)  ERROR: Failed to read fsgrid variable " << variableName << endl << write;
         return false;
      }
      
      // Assign buffer into fsgrid
      int index=0;
      for(FsGridTools::FsIndex_t z=0; z<localSize[2]; z++) {
         for(FsGridTools::FsIndex_t y=0; y<localSize[1]; y++) {
            for(FsGridTools::FsIndex_t x=0; x<localSize[0]; x++) {
               memcpy(targetGrid.get(x,y,z), &buffer[index], N*sizeof(Real));
               index += N;
            }
         }
      }
      
   } else {

      // More difficult case: different number of tasks.
      // In this case, our own fsgrid domain overlaps (potentially many) domains in the file.
      // We read the whole source rank into a temporary buffer, and transfer the overlapping
      // part.
      //
      // +------------+----------------+
      // |            |                |
      // |    . . . . . . . . . . . .  |
      // |    .<----->|<----------->.  |
      // |    .<----->|<----------->.  |
      // |    .<----->|<----------->.  |
      // +----+-------+-------------+--|
      // |    .<----->|<----------->.  |
      // |    .<----->|<----------->.  |
      // |    .<----->|<----------->.  |
      // |    . . . . . . . . . . . .  |
      // |            |                |
      // +------------+----------------+

      // Iterate through tasks and find their overlap with our domain.
      size_t fileOffset = 0;
      for(int task = 0; task < numWritingRanks; task++) {

         phiprof::Timer taskArithmetics1 {"task overlap arithmetics 1"};

         std::array<FsGridTools::FsIndex_t,3> thatTasksSize;
         std::array<FsGridTools::FsIndex_t,3> thatTasksStart;
         thatTasksSize[0] = targetGrid.calcLocalSize(globalSize[0], fileDecomposition[0], task/fileDecomposition[2]/fileDecomposition[1]);
         thatTasksSize[1] = targetGrid.calcLocalSize(globalSize[1], fileDecomposition[1], (task/fileDecomposition[2])%fileDecomposition[1]);
         thatTasksSize[2] = targetGrid.calcLocalSize(globalSize[2], fileDecomposition[2], task%fileDecomposition[2]);

         thatTasksStart[0] = targetGrid.calcLocalStart(globalSize[0], fileDecomposition[0], task/fileDecomposition[2]/fileDecomposition[1]);
         thatTasksStart[1] = targetGrid.calcLocalStart(globalSize[1], fileDecomposition[1], (task/fileDecomposition[2])%fileDecomposition[1]);
         thatTasksStart[2] = targetGrid.calcLocalStart(globalSize[2], fileDecomposition[2], task%fileDecomposition[2]);

         // Iterate through overlap area
         std::array<FsGridTools::FsIndex_t,3> overlapStart,overlapEnd,overlapSize;
         overlapStart[0] = max(localStart[0],thatTasksStart[0]);
         overlapStart[1] = max(localStart[1],thatTasksStart[1]);
         overlapStart[2] = max(localStart[2],thatTasksStart[2]);

         overlapEnd[0] = min(localStart[0]+localSize[0], thatTasksStart[0]+thatTasksSize[0]);
         overlapEnd[1] = min(localStart[1]+localSize[1], thatTasksStart[1]+thatTasksSize[1]);
         overlapEnd[2] = min(localStart[2]+localSize[2], thatTasksStart[2]+thatTasksSize[2]);

         overlapSize[0] = max(overlapEnd[0]-overlapStart[0],(FsGridTools::FsIndex_t)0);
         overlapSize[1] = max(overlapEnd[1]-overlapStart[1],(FsGridTools::FsIndex_t)0);
         overlapSize[2] = max(overlapEnd[2]-overlapStart[2],(FsGridTools::FsIndex_t)0);

         taskArithmetics1.stop();

         // Read into buffer
         std::vector<Real> buffer(thatTasksSize[0]*thatTasksSize[1]*thatTasksSize[2]*N);

         phiprof::Timer multiRead {"multiRead"};
         file.startMultiread("VARIABLE", attribs);
         // Read every source rank that we have an overlap with.
         if(overlapSize[0]*overlapSize[1]*overlapSize[2] > 0) {


            if(!convertFloatType) {
               if(file.addMultireadUnit((char*)buffer.data(), thatTasksSize[0]*thatTasksSize[1]*thatTasksSize[2])==false) {
                  logFile << "(RESTART)  ERROR: Failed to read fsgrid variable " << variableName << endl << write;
                  return false;
               }
               file.endMultiread(fileOffset);
            } else {
               std::vector<float> readBuffer(thatTasksSize[0]*thatTasksSize[1]*thatTasksSize[2]*N);
               if(file.addMultireadUnit((char*)readBuffer.data(), thatTasksSize[0]*thatTasksSize[1]*thatTasksSize[2])==false) {
                  logFile << "(RESTART)  ERROR: Failed to read fsgrid variable " << variableName << endl << write;
                  return false;
               }
               file.endMultiread(fileOffset);

               for(uint64_t i=0; i< thatTasksSize[0]*thatTasksSize[1]*thatTasksSize[2]*N; i++) {
                  buffer[i]=readBuffer[i];
               }
            }

            // Copy continuous stripes in x direction.
            for(FsGridTools::FsIndex_t z=overlapStart[2]; z<overlapEnd[2]; z++) {
               for(FsGridTools::FsIndex_t y=overlapStart[1]; y<overlapEnd[1]; y++) {
                  for(FsGridTools::FsIndex_t x=overlapStart[0]; x<overlapEnd[0]; x++) {
                     FsGridTools::FsIndex_t index = (z - thatTasksStart[2]) * thatTasksSize[0]*thatTasksSize[1]
                        + (y - thatTasksStart[1]) * thatTasksSize[0]
                        + (x - thatTasksStart[0]);

                     memcpy(targetGrid.get(x - localStart[0], y - localStart[1], z - localStart[2]), &buffer[index*N], N*sizeof(Real));
                  }
               }
            }
         } else {
            // If we don't overlap, just perform a dummy read.
            file.endMultiread(fileOffset);
         }
         fileOffset += thatTasksSize[0] * thatTasksSize[1] * thatTasksSize[2];
         multiRead.stop();
      }
   }
   phiprof::Timer updateGhostsTimer {"updateGhostCells"};
   targetGrid.updateGhostCells();
   updateGhostsTimer.stop();
   return true;
}

/*! Read an ionosphere variable from the given vlsv file.
 * Note that only singular floating point values (no vectors) can be read at this time.
 * \param file VLSV parallel reader with a file open.
 * \param variableName Name of the variable in the file
 * \param grid the ionosphere grid that data will be deposited into
 * \param index index into the nodes' parameters array, where the data will end up.
 */
bool readIonosphereNodeVariable(
   vlsv::ParallelReader& file, const string& variableName, SBC::SphericalTriGrid& grid, ionosphereParameters index) {

   uint64_t arraySize;
   uint64_t vectorSize;
   vlsv::datatype::type dataType;
   uint64_t byteSize;
   list<pair<string,string> > attribs;
   
   attribs.push_back(make_pair("name",variableName));
   attribs.push_back(make_pair("mesh","ionosphere"));

   // If we don't have an ionosphere (zero nodes), we simply skip trying to read any restart data for this.
   if(grid.nodes.size() == 0) {
      return true;
   }

   if (file.getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
      logFile << "(RESTART)  ERROR: Failed to read array info for " << variableName << endl << write;
      return false;
   }

   // Verify that this is a scalar variable
   if(vectorSize != 1) {
      logFile << "(RESTART) ERROR: Trying to read vector valued (" << vectorSize << " components) ionosphere parameter from restart file. Only scalars are supported." << endl << write;
      return false;
   }

   // Verify that the size matches our constructed ionosphere object
   if(grid.nodes.size() != arraySize) {
      logFile << "(RESTART) ERROR: Ionosphere restart size mismatch: trying to read variable " << variableName << " with " << arraySize << " values into a ionosphere grid with " << grid.nodes.size() << " nodes!" << endl << write;
      return false;
   }

   std::vector<Real> buffer(arraySize);
   if(file.readArray("VARIABLE", attribs, 0, arraySize, buffer.data()) == false) {
      logFile << "(RESTART) ERROR: Failed to read ionosphere variable " << variableName << endl << write;
   }

   for(uint i=0; i<grid.nodes.size(); i++) {
      grid.nodes[i].parameters[index] = buffer[i];
   }

   return true;
}

/*!
\brief Read in state from a vlsv file in order to restart simulations
\param mpiGrid Vlasiator's grid
\param name Name of the restart file e.g. "restart.00052.vlsv"
 \return Returns true if the operation was successful
 \sa readGrid
 */
bool exec_readGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
                   const std::string& name) {
   vector<CellID> fileCells; /*< CellIds for all cells in file*/
   vector<size_t> nBlocks;/*< Number of blocks for all cells in file*/
   bool success=true;
   int myRank,processes;

   // Note: Spatial grid name hard-coded here.
   // But so are the other mesh names below.
   const string meshName = "SpatialGrid";
   
   // Attempt to open VLSV file for reading:
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   MPI_Comm_size(MPI_COMM_WORLD,&processes);

   phiprof::Timer readGridTimer {"readGrid"};

   phiprof::Timer readScalarsTimer {"readScalars"};

   vlsv::ParallelReader file;

   MPI_Info MPIinfo;
   if (P::restartReadHints.size() == 0) {
      MPIinfo = MPI_INFO_NULL;
   } else {
      MPI_Info_create(&MPIinfo);
      
      for (std::vector<std::pair<std::string,std::string>>::const_iterator it = P::restartReadHints.begin();
           it != P::restartReadHints.end();
           it++)
      {
         MPI_Info_set(MPIinfo, it->first.c_str(), it->second.c_str());
      }
   }

   if (file.open(name,MPI_COMM_WORLD,MASTER_RANK,MPIinfo) == false) {
      success=false;
   }
   exitOnError(success,"(RESTART) Could not open file",MPI_COMM_WORLD);

   // Around May 2015 time was renamed from "t" to "time", we try to read both, 
   // new way is read first
   if (readScalarParameter(file,"time",P::t,MASTER_RANK,MPI_COMM_WORLD) == false)
     if (readScalarParameter(file,"t", P::t,MASTER_RANK,MPI_COMM_WORLD) == false) 
       success=false;
   P::t_min=P::t;

   // Around May 2015 timestep was renamed from "tstep" to "timestep", we to read
   // both, new way is read first
   if (readScalarParameter(file,"timestep",P::tstep,MASTER_RANK,MPI_COMM_WORLD) == false)
     if (readScalarParameter(file,"tstep", P::tstep,MASTER_RANK,MPI_COMM_WORLD) ==false) 
       success = false;
   P::tstep_min=P::tstep;

   if(readScalarParameter(file,"dt",P::dt,MASTER_RANK,MPI_COMM_WORLD) ==false) success=false;

   if(readScalarParameter(file,"fieldSolverSubcycles",P::fieldSolverSubcycles,MASTER_RANK,MPI_COMM_WORLD) ==false) {
      // Legacy restarts do not have this field, it "should" be safe for one or two steps...
      P::fieldSolverSubcycles = 1.0;
      cout << " No P::fieldSolverSubcycles found in restart, setting 1." << endl;
   }
   MPI_Bcast(&(P::fieldSolverSubcycles),1,MPI_Type<uint>(),MASTER_RANK,MPI_COMM_WORLD);
   



   checkScalarParameter(file,"xmin",P::xmin,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"ymin",P::ymin,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"zmin",P::zmin,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"xmax",P::xmax,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"ymax",P::ymax,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"zmax",P::zmax,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"xcells_ini",P::xcells_ini,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"ycells_ini",P::ycells_ini,MASTER_RANK,MPI_COMM_WORLD);
   checkScalarParameter(file,"zcells_ini",P::zcells_ini,MASTER_RANK,MPI_COMM_WORLD);

   readScalarsTimer.stop();

   phiprof::Timer readLayoutimer {"readDatalayout"};
   if (success) {
      success = readCellIds(file,fileCells,MASTER_RANK,MPI_COMM_WORLD);
   }

   // Check that the cellID lists are identical in file and grid
   if (myRank==0) {
      vector<CellID> allGridCells = mpiGrid.get_all_cells();
      if (fileCells.size() != allGridCells.size()) {
         std::cout << "File has " << fileCells.size() << " cells, got " << allGridCells.size() << " cells!" << std::endl;
         success=false;
      }
   }
   
   exitOnError(success,"(RESTART) Wrong number of cells in restart file",MPI_COMM_WORLD);

   // Read the total number of velocity blocks in each spatial cell.
   // Note that this is a sum over all existing particle species.
   if (success == true) {
      success = readNBlocks(file,meshName,nBlocks,MASTER_RANK,MPI_COMM_WORLD);
   }

   //make sure all cells are empty, we will anyway overwrite everything and 
   // in that case moving cells is easier...
     {
        const vector<CellID>& gridCells = getLocalCells();
        for (size_t i=0; i<gridCells.size(); i++) {
           for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
             mpiGrid[gridCells[i]]->clear(popID);
        }
     }

   uint64_t totalNumberOfBlocks=0;
   unsigned int numberOfBlocksPerProcess;
   for(uint i=0; i<nBlocks.size(); ++i){
      totalNumberOfBlocks += nBlocks[i];
   }
   numberOfBlocksPerProcess= 1 + totalNumberOfBlocks/processes;

   uint64_t localCellStartOffset=0; // This is where local cells start in file-list after migration.
   uint64_t localCells=0;
   uint64_t numberOfBlocksCount=0;
   
   // Pin local cells to remote processes, we try to balance number of blocks so that 
   // each process has the same amount of blocks, more or less.
   for (size_t i=0; i<fileCells.size(); ++i) {
      numberOfBlocksCount += nBlocks[i];
      int newCellProcess = numberOfBlocksCount/numberOfBlocksPerProcess;
      if (newCellProcess == myRank) {
         if (localCells == 0)
            localCellStartOffset=i; //here local cells start
         ++localCells;
      }
      if (mpiGrid.is_local(fileCells[i])) {
         mpiGrid.pin(fileCells[i],newCellProcess);
      }
   }

   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);

   //Do initial load balance based on pins. Need to transfer at least sysboundaryflags
   mpiGrid.balance_load(false);

   //update list of local gridcells
   recalculateLocalCellsCache(mpiGrid);

   //get new list of local gridcells
   const vector<CellID>& gridCells = getLocalCells();

   // Unpin cells, otherwise we will never change this initial bad balance
   for (size_t i=0; i<gridCells.size(); ++i) {
      mpiGrid.unpin(gridCells[i]);
   }

   // Check for errors, has migration succeeded
   if (localCells != gridCells.size() ) {
      success=false;
   } 

   if (success == true) {
      for (uint64_t i=localCellStartOffset; i<localCellStartOffset+localCells; ++i) {
         if(mpiGrid.is_local(fileCells[i]) == false) {
            success = false;
         }
      }
   }

   exitOnError(success,"(RESTART) Cell migration failed",MPI_COMM_WORLD);

   // Set cell coordinates based on cfg (mpigrid) information
   for (size_t i=0; i<gridCells.size(); ++i) {
      array<double, 3> cell_min = mpiGrid.geometry.get_min(gridCells[i]);
      array<double, 3> cell_length = mpiGrid.geometry.get_length(gridCells[i]);

      mpiGrid[gridCells[i]]->parameters[CellParams::XCRD] = cell_min[0];
      mpiGrid[gridCells[i]]->parameters[CellParams::YCRD] = cell_min[1];
      mpiGrid[gridCells[i]]->parameters[CellParams::ZCRD] = cell_min[2];
      mpiGrid[gridCells[i]]->parameters[CellParams::DX  ] = cell_length[0];
      mpiGrid[gridCells[i]]->parameters[CellParams::DY  ] = cell_length[1];
      mpiGrid[gridCells[i]]->parameters[CellParams::DZ  ] = cell_length[2];
   }

   // Where local data start in the blocklists
   //uint64_t localBlocks=0;
   //for(uint64_t i=localCellStartOffset; i<localCellStartOffset+localCells; ++i) {
   //  localBlocks += nBlocks[i];
   //}
   readLayoutimer.stop();

   //todo, check file datatype, and do not just use double
   phiprof::Timer readParametersTimer {"readCellParameters"};
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"moments",CellParams::RHOM,5,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"moments_dt2",CellParams::RHOM_DT2,5,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"moments_r",CellParams::RHOM_R,5,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"moments_v",CellParams::RHOM_V,5,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"pressure",CellParams::P_11,3,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"pressure_dt2",CellParams::P_11_DT2,3,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"pressure_r",CellParams::P_11_R,3,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"pressure_v",CellParams::P_11_V,3,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"LB_weight",CellParams::LBWEIGHTCOUNTER,1,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"max_v_dt",CellParams::MAXVDT,1,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"max_r_dt",CellParams::MAXRDT,1,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"max_fields_dt",CellParams::MAXFDT,1,mpiGrid); }
   if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"vg_drift",CellParams::BULKV_FORCING_X,3,mpiGrid); }
   if (P::refineOnRestart) {
      // Refinement indices alpha_1 and alpha_2
      if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"vg_amr_alpha1",CellParams::AMR_ALPHA1,1,mpiGrid); }
      if(success) { success=readCellParamsVariable(file,fileCells,localCellStartOffset,localCells,"vg_amr_alpha2",CellParams::AMR_ALPHA2,1,mpiGrid); }
   }

   // Backround B has to be set, there are also the derivatives that should be written/read if we wanted to only read in background field
   readParametersTimer.stop();

   phiprof::Timer readBlocksTimer {"readBlockData"};
   if (success == true) {
      success = readBlockData(file,meshName,fileCells,localCellStartOffset,localCells,mpiGrid); 
   }
   readBlocksTimer.stop();

   phiprof::Timer updateNeighborsTimer {"updateMpiGridNeighbors"};
   mpiGrid.update_copies_of_remote_neighbors(FULL_NEIGHBORHOOD_ID);
   updateNeighborsTimer.stop();
   
   phiprof::Timer readfsTimer {"readFsGrid"};
   // Read fsgrid data back in
   int fsgridInputRanks=0;
   phiprof::Timer tReadScalarParameter {"readScalarParameter"};
   if(readScalarParameter(file,"numWritingRanks",fsgridInputRanks, MASTER_RANK, MPI_COMM_WORLD) == false) {
      exitOnError(false, "(RESTART) FSGrid writing rank number not found in restart file", MPI_COMM_WORLD);
   }
   tReadScalarParameter.stop();

   if (success) { success = readFsGridVariable(file, "fg_PERB", fsgridInputRanks, perBGrid); }
   if (success) { success = readFsGridVariable(file, "fg_E", fsgridInputRanks, EGrid); }
   exitOnError(success,"(RESTART) Failure reading fsgrid restart variables",MPI_COMM_WORLD);
   readfsTimer.stop();
   
   phiprof::Timer readIonosphereTimer {"readIonosphere"};
   bool ionosphereSuccess=true;
   ionosphereSuccess = readIonosphereNodeVariable(file, "ig_fac", SBC::ionosphereGrid, ionosphereParameters::SOURCE);
   // Reconstruct source term by multiplying the fac density with the element area
   for(uint i = 0; i<SBC::ionosphereGrid.nodes.size(); i++) {
      Real area = 0;
      for(uint e=0; e< SBC::ionosphereGrid.nodes[i].numTouchingElements; e++) {
         area += SBC::ionosphereGrid.elementArea(SBC::ionosphereGrid.nodes[i].touchingElements[e]);
      }
      area /= 3.; // As every element has 3 corners, don't double-count areas
      SBC::ionosphereGrid.nodes[i].parameters[ionosphereParameters::SOURCE] *= area;
   }
   ionosphereSuccess &= readIonosphereNodeVariable(file, "ig_rhon", SBC::ionosphereGrid, ionosphereParameters::RHON);
   ionosphereSuccess &= readIonosphereNodeVariable(file, "ig_electrontemp", SBC::ionosphereGrid, ionosphereParameters::TEMPERATURE);
   ionosphereSuccess &= readIonosphereNodeVariable(file, "ig_potential", SBC::ionosphereGrid, ionosphereParameters::SOLUTION);
   if(!ionosphereSuccess) {
      logFile << "(RESTART) Reading ionosphere variables failed. Continuing anyway. Variables will be zero, assuming this is an ionosphere cold start?" << std::endl;
   }

   // Read additional variables that are not formally required for solving the
   // ionosphere, but help making the first output consistent if ionosphere
   // timestep is very large.
   // If these are missing from the restart file, we are fine continuing with
   // zeros.
   bool ionosphereOptionalSuccess = readIonosphereNodeVariable(file, "ig_sigmah", SBC::ionosphereGrid, ionosphereParameters::SIGMAH);
   ionosphereOptionalSuccess &= readIonosphereNodeVariable(file, "ig_sigmap", SBC::ionosphereGrid, ionosphereParameters::SIGMAP);
   ionosphereOptionalSuccess &= readIonosphereNodeVariable(file, "ig_sigmaparallel", SBC::ionosphereGrid, ionosphereParameters::SIGMAPARALLEL);
   ionosphereOptionalSuccess &= readIonosphereNodeVariable(file, "ig_precipitation", SBC::ionosphereGrid, ionosphereParameters::PRECIP);
   if(ionosphereSuccess && !ionosphereOptionalSuccess) {
      logFile << "(RESTART) Restart file contains no ionosphere conductivity data. Ionosphere will run fine, but first output bulk file might have bogus conductivities." << std::endl;
   }
   readIonosphereTimer.stop();

   success = file.close();

   exitOnError(success,"(RESTART) Other failure",MPI_COMM_WORLD);
   return success;
}

//FIXME, readGrid has no support for checking or converting endianness
/*!
\brief Read in state from a vlsv file in order to restart simulations
\param mpiGrid Vlasiator's grid
\param name Name of the restart file e.g. "restart.00052.vlsv"
*/
bool readGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
              const std::string& name){
   //Check the vlsv version from the file:
   return exec_readGrid(mpiGrid,perBGrid,EGrid,technicalGrid,name);
}

/*!
\brief Refine the grid to be identical to the file's
\param mpiGrid Vlasiator's grid
\param name Name of the restart file e.g. "restart.00052.vlsv"
*/
bool readFileCells(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const std::string& name)
{
   phiprof::Timer readCellIdsTimer {"Restart read File cellIDs"};
   vector<CellID> fileCells; /*< CellIds for all cells in file*/
   bool success = true;
   vlsv::ParallelReader file;
   MPI_Info mpiInfo = MPI_INFO_NULL;

   // Not sure if this success business is useful at all...
   success = file.open(name,MPI_COMM_WORLD,MASTER_RANK,mpiInfo);
   exitOnError(success,"(READ_FILE_CELLS) Could not open file",MPI_COMM_WORLD);

   readCellIds(file,fileCells,MASTER_RANK,MPI_COMM_WORLD);
   phiprof::Timer loadCellsTimer {"load CellIDs into grid"};
   success = mpiGrid.load_cells(fileCells);
   exitOnError(success,"(READ_FILE_CELLS) Failed to refine grid",MPI_COMM_WORLD);
   loadCellsTimer.stop();

   success = file.close();
   exitOnError(success,"(READ_FILE_CELLS) Other error",MPI_COMM_WORLD);
   return success;
}


bool readFsgridDecomposition(vlsv::ParallelReader& file, std::array<FsGridTools::Task_t,3>& decomposition){
   list<pair<string,string> > attribs;

   int myRank;   
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   phiprof::Timer readFsGridDecomposition {"readFsGridDecomposition"};

   attribs.push_back(make_pair("mesh","fsgrid"));

   std::array<FsGridTools::FsSize_t,3> gridSize;
   FsGridTools::FsSize_t* gridSizePtr = &gridSize[0];
   bool success = file.read("MESH_BBOX",attribs, 0, 3, gridSizePtr, false);
   if(success == false){
      exitOnError(false, "(RESTART) FSGrid gridsize not found in file.", MPI_COMM_WORLD);
      return false;
   }

   std::array<FsGridTools::Task_t,3> fsGridDecomposition={0,0,0}; 
   FsGridTools::Task_t* ptr = &fsGridDecomposition[0];

   success = file.read("MESH_DECOMPOSITION",attribs, 0, 3, ptr, false);
   if (success == false) {
      if (myRank == MASTER_RANK){
         std::cout << "Could not read MESH_DECOMPOSITION, attempting to calculate it from MESH." << endl;
      }
      int fsgridInputRanks=0;
      if(file.readParameter("numWritingRanks",fsgridInputRanks) == false) {
         exitOnError(false, "(RESTART) FSGrid writing rank number not found in restart file.", MPI_COMM_WORLD);
         return false;
      }

      int64_t* domainInfo = NULL;
      success = file.read("MESH_DOMAIN_SIZES",attribs, 0, fsgridInputRanks, domainInfo);
      if(success == false){
         if (myRank == MASTER_RANK){
            std::cerr << "Could not read MESH_DOMAIN_SIZES from file" << endl;
         }
         return false;
      }
      std::vector<uint64_t> mesh_domain_sizes;
      for (int i = 0; i < 2*fsgridInputRanks; i+=2){
         mesh_domain_sizes.push_back(domainInfo[i]);
      }
      list<pair<string,string> > mesh_attribs;
      mesh_attribs.push_back(make_pair("name","fsgrid"));
      std::vector<FsGridTools::FsSize_t> rank_first_ids(fsgridInputRanks);
      FsGridTools::FsSize_t* ids_ptr = &rank_first_ids[0];

      std::set<FsGridTools::FsIndex_t> x_corners, y_corners, z_corners;
      
      int64_t begin_rank = 0;
      for(auto rank_size : mesh_domain_sizes){
         if(myRank == MASTER_RANK){
            if(file.read("MESH", mesh_attribs, begin_rank, 1, ids_ptr, false) == false){
               if (myRank == 0){
                  std::cerr << "Reading MESH failed.\n";
               }
               return false;
            }
            std::array<FsGridTools::FsIndex_t,3> inds = FsGridTools::globalIDtoCellCoord(*ids_ptr, gridSize);
            x_corners.insert(inds[0]);
            y_corners.insert(inds[1]);
            z_corners.insert(inds[2]);
            ++ids_ptr;
            begin_rank += rank_size;
         } else {
            file.read("MESH", mesh_attribs, begin_rank, 0, ids_ptr, false);
         }
      }

      decomposition[0] = x_corners.size();
      decomposition[1] = y_corners.size();
      decomposition[2] = z_corners.size();
      MPI_Bcast(&decomposition, 3, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

      if(decomposition[0]*decomposition[1]*decomposition[2] == fsgridInputRanks){
         if (myRank == MASTER_RANK){
            std::cout << "Fsgrid decomposition computed from MESH to be " << decomposition[0] << " " << decomposition[1] << " " <<decomposition[2] << endl;
         }
         return true;
      } else {
         if (myRank == MASTER_RANK){
            std::cout << "Fsgrid decomposition computed from MESH to be " << decomposition[0] << " " << decomposition[1] << " " <<decomposition[2] << ", which is not compatible with numWritingRanks ("<<fsgridInputRanks << ")" << endl;
         }
         return false;
      }

   } else {
      decomposition[0] = fsGridDecomposition[0];
      decomposition[1] = fsGridDecomposition[1];
      decomposition[2] = fsGridDecomposition[2];
      // logFile << "(RESTART) Fsgrid decomposition read as " << decomposition[0] << " " << decomposition[1] << " " <<decomposition[2] << "\n";
      return true;
   }
   return false;
}
