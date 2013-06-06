#include <cstdlib>
#include <iostream>
#include <iomanip> // for setprecision()
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>
#include <array>
#include "iowrite.h"
#include "phiprof.hpp"
#include "parameters.h"
#include "logger.h"
#include "vlsvwriter2.h"
#include "vlsvreader2.h"
#include "vlsv_writer.h"
#include "vlasovmover.h"

using namespace std;
using namespace phiprof;
using namespace vlsv;

extern Logger logFile, diagnostic;

typedef Parameters P;


//O: Should Zone be moved to another file?
class Zone {
   private:
      //Cell id of the zone -- Note: zone is basically a cell
      uint64_t cellId;
      //local id of the zone
      uint64_t localId;
      //MPI process rank owning the zone:
      //O: NOTE: I actually don't know what this means
      int myRank;
      //Has the zone's indices
      array<unsigned int, 3> indices;
      //Check initialization function for input parameters of int
      template<typename T>
      bool checkInitialization( const T check ) {
         if( check == numeric_limits<T>::max() ) {
            //Not initialized
            cerr << "Called uninitialized cell variable!" << endl;
            return false;
         } else {
            //everything ok
            return true;
         }
      }
      template<typename T>
      bool checkInitialization( const T check ) const {
         if( check == numeric_limits<T>::max() ) {
            //Not initialized
            cerr << "Called uninitialized cell variable!" << endl;
            return false;
         } else {
            //everything ok
            return true;
         }
      }
   public:
      uint64_t getCellId() const {
         checkInitialization( cellId );
         return cellId;
      }
      uint64_t getCellId() {
         checkInitialization( cellId );
         return cellId;
      }
      void setGeometry( const uint64_t _cellId, const array<uint64_t, 3> & cell_bounds ) {
         cellId = _cellId;
         
         //Set cell indices:
         indices[0] = cellId % cell_bounds[0];
         indices[1] = ((cellId - indices[0]) / cell_bounds[0]) % cell_bounds[1];
         indices[2] = ((cellId - cell_bounds[0]*indices[1]) / (cell_bounds[0]*cell_bounds[1]));
      }
      uint64_t getLocalId() const {
         checkInitialization( localId );
         return localId;
      }
      uint64_t getLocalId() {
         checkInitialization( localId );
         return localId;
      }
      void setLocalId( const unsigned int _localId ) {
         localId = _localId;
      }
      int getRank() {
         checkInitialization( myRank );
         return myRank;
      }
      int getRank() const {
         checkInitialization( myRank );
         return myRank;
      }
      void setRank( const int rank ) {
         myRank = rank;
      }
      const array<unsigned int, 3> & getIndices() {
         return indices;
      }
      const array<unsigned int, 3> & getIndices() const {
         return indices;
      }
      Zone() {
         //Put on default values
         cellId = numeric_limits<uint64_t>::max();
         localId = numeric_limits<uint64_t>::max();
         myRank = numeric_limits<int>::max();
         indices[0] = 0;
         indices[1] = 0;
         indices[2] = 0;
      }
      ~Zone() {
      }
};

// array<unsigned int, 3> & cell_bounds

bool createZone(  const dccrg::Dccrg<SpatialCell> & mpiGrid,
                  const vector<uint64_t> & local_cells,
                  const vector<uint64_t> & ghost_cells,
                  vector<Zone> & local_zones,
                  vector<Zone> & ghost_zones ) {
   if( local_zones.empty() || ghost_zones.empty() ) {
      cerr << "Error at: ";
      cerr << __FILE__ << " " << __LINE__;
      cerr << ", zone empty!" << endl;
      return false;
   }
   //Get the cell boundaries (Tells how many cells there are in x, y, z direction)
   //Note: P::xcells_ini stands for the number of cells in x direction
   const array<uint64_t, 3> cell_bounds { { P::xcells_ini, P::ycells_ini, P::zcells_ini } };

   //Declare an iterator for iterating though the cell ids
   vector<uint64_t>::const_iterator it;
   //Local ids for the process start from 0 (this is used in the iteration)
   uint64_t localId = 0;
   //Iterate through local cells and save the data into zones
   for( it = local_cells.begin(); it != local_cells.end(); ++it ) {
      //Declare Zone -- this will be appended to vector<Zone> & zone
      Zone _zone;
      const uint64_t cellId = (*it);
      //Set the cell id and indices (cell coordinates):
      _zone.setGeometry( cellId, cell_bounds );
      //Set the cell's MPI process:
      _zone.setRank( mpiGrid.get_process( cellId ) );

      //O: FOR DEBUGGING: (REMOVE)
      cout << _zone.getRank() << endl;

      //O: Set the local id:
      //O: Hmm, I'm not sure whether MPI is needed for the local ids if every process has their own local ids
      //O: Let's do it here anyway in case it turns out it's needed
      /*
      SpatialCell * cell = mpiGrid[cellId];
      cell->ioLocalCellId = localId;
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_IOLOCALCELLID);
      mpiGrid.update_remote_neighbor_data(NEAREST_NEIGHBORHOOD_ID);
      */

      //Set the local id:
      _zone.setLocalId( localId );

      //Append _zone to the list of zones:
      local_zones.push_back( _zone );
      //Increment the local id
      localId++;
   }
   //Do the same for ghost_zones:
   for( it = ghost_cells.begin(); it != ghost_cells.end(); ++it ) {
      //Declare Zone -- this will be appended to vector<Zone> & zone
      Zone _zone;
      const uint64_t cellId = (*it);
      //Set the cell id and indices (cell coordinates):
      _zone.setGeometry( cellId, cell_bounds );
      //Set the cell's MPI process:
      _zone.setRank( mpiGrid.get_process( cellId ) );

      //O: FOR DEBUGGING: (REMOVE)
      cout << _zone.getRank() << endl;

      //O: Set the local id:
      //O: Hmm, I'm not sure whether MPI is needed for the local ids if every process has their own local ids
      //O: Let's do it here anyway in case it turns out it's needed
      /*
      SpatialCell * cell = mpiGrid[cellId];
      cell->ioLocalCellId = localId;
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_IOLOCALCELLID);
      mpiGrid.update_remote_neighbor_data(NEAREST_NEIGHBORHOOD_ID);
      */

      //Set the local id:
      _zone.setLocalId( localId );

      //Append _zone to the list of zones:
      ghost_zones.push_back( _zone );
      //Increment the local id
      localId++;
   }
   return true;
}

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

bool o_writeVelocityDistributionData(
                                    Writer& vlsvWriter,
                                    const dccrg::Dccrg<SpatialCell>& mpiGrid,
                                    const vector<const Zone*> & zones,
                                    MPI_Comm comm
                                    ) {
   // O: NOTE: I'm not entirely sure this works with the VisIt plugin
   // Write velocity blocks and related data. 
   // In restart we just write velocity grids for all cells.
   // First write global Ids of those cells which write velocity blocks (here: all cells):
   map<string,string> attribs;
   bool success=true;

   //Make sure zones is not empty:
   if( zones.empty() ) {
      cerr << "Error, passed empty zones vector at: " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   //O: I guess the problem with writeArray and Zone class implementation shows up here
   //I need to create vector called cells based on the Zones class
   vector<const Zone*>::const_iterator it;
   vector<uint64_t> cells;
   //Iterate through zones and put data into std::vector cells:
   for( it = zones.begin(); it != zones.end(); ++it ) {
      //Get the zone that is currently being iterated:
      const Zone * zone = *it;
      //Make sure zone is not empty:
      if( !zone ) {
         cerr << "Error, passed NULL pointer at: " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      cells.push_back( zone->getCellId() );
   }
   //Make sure everything went smoothly:
   if( cells.empty() ) {
      cerr << "Error, passed empty zones vector at: " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   //Compute totalBlocks
   uint64_t totalBlocks = 0;  
   vector<uint> blocksPerCell;   
   for(size_t cell=0;cell<cells.size();++cell){
      totalBlocks+=mpiGrid[cells[cell]]->number_of_blocks;
      blocksPerCell.push_back(mpiGrid[cells[cell]]->number_of_blocks);
   }
   
   //O: NOT SURE IF THIS IS CORRECT
   attribs["mesh"] = "SpatialGrid";
   if (vlsvWriter.writeArray("CELLSWITHBLOCKS",attribs,cells.size(),1,cells.data()) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;
   //Write blocks per cell, this has to be in the same order as cellswitblocks so that extracting works
   if(vlsvWriter.writeArray("BLOCKSPERCELL",attribs,blocksPerCell.size(),1,blocksPerCell.data()) == false) success = false;
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
   } catch (...) {
      success=false;
   }
   
   if( globalSuccess(success,"(MAIN) writeGrid: ERROR: Failed to fill temporary array velocityBlockParameters",MPI_COMM_WORLD) == false) {
      vlsvWriter.close();
      return false;
   }
   
   
   attribs.clear();
   if (vlsvWriter.writeArray("BLOCKCOORDINATES",attribs,totalBlocks,BlockParams::N_VELOCITY_BLOCK_PARAMS,velocityBlockParameters.data()) == false) success = false;
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
   //O: MAKE SURE THIS IS CORRECT!
   attribs["mesh"] = "SpatialGrid";
   attribs["name"] = "avgs";
   if (vlsvWriter.writeArray("BLOCKVARIABLE",attribs,totalBlocks,SIZE_VELBLOCK,velocityBlockData.data()) == false) success=false;
   if (success ==false)      logFile << "(MAIN) writeGrid: ERROR occurred when writing BLOCKVARIABLE f" << endl << writeVerbose;
   velocityBlockData.clear();
    
   return success;
}


bool o_writeDataReducer(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                      const vector<Zone> & local_zones,
                      DataReducer& dataReducer,
                      int dataReducerIndex,
                      Writer& vlsvWriter){
   //Make sure local_zones is not empty:
   if ( local_zones.empty() ) {
      cerr << "ERROR, passed an empty local_zones vector at: " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   map<string,string> attribs;
   string variableName, dataType;
   bool success=true;
   
   uint dataSize,vectorSize;
   //Input the variable name:
   variableName = dataReducer.getName(dataReducerIndex);
   attribs["name"] = variableName;
   //Input the mesh name:
   attribs["mesh"] = "SpatialGrid";
   if (dataReducer.getDataVectorInfo(dataReducerIndex,dataType,dataSize,vectorSize) == false) {
      cerr << "ERROR when requesting info from DRO " << dataReducerIndex << endl;
      return false;
   }

   //Get arraySize:
   uint64_t arraySize = local_zones.size()*vectorSize*dataSize;
   
   //Request DataReductionOperator to calculate the reduced data for all local cells:
   //O: NOTE: Might need to be cast into proper data type
   char* varBuffer = new char[arraySize];
   //Iterate through zones:
   vector<Zone>::const_iterator zone = local_zones.begin();
   for ( uint i = 0; zone != local_zones.end(); ++zone, ++i ) {
      //Get the cell id of the zone:
      const uint64_t cellId = zone->getCellId();
      if (dataReducer.reduceData(mpiGrid[cellId],dataReducerIndex,varBuffer + i*vectorSize*dataSize) == false){
         success = false;
      }
      
      if (success == false){
         logFile << "(MAIN) writeGrid: ERROR datareductionoperator '" << dataReducer.getName(dataReducerIndex) <<
            "' returned false!" << endl << writeVerbose;
      }
   }
   // Write  reduced data to file if DROP was successful:
   if(success){
      if (vlsvWriter.writeArray("VARIABLE", attribs, dataType, arraySize, vectorSize, dataSize, varBuffer) == false)
         success = false;
      if (success == false){
         logFile << "(MAIN) writeGrid: ERROR failed to write datareductionoperator data to file!" << endl << writeVerbose;
      }
   }
   delete[] varBuffer;
   varBuffer = NULL;
   return success;
}


bool o_writeVariables( const dccrg::Dccrg<SpatialCell>& mpiGrid,
                       DataReducer& dataReducer,
                       Writer & vlsvWriter, 
                       const uint & index,
                       const string & meshName, 
                       const vector<Zone> & local_zones ) {
   //Make sure local_zones is not empty:
   if ( local_zones.empty() ) {
      cerr << "ERROR, passed an empty local_zones vector at: " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //O: NOTE: '
   //Compute which zones will write out their velocity space
   vector<const Zone*> velSpaceZones;
   int lineX, lineY, lineZ;
   //Iterate through local_zones
   //O: NOTE: Might need to remove the const correctness unfortunately
   vector<Zone>::const_iterator it;
   for( it = local_zones.begin(); it != local_zones.end(); ++it ) {
      //Get the zone's cell id:
      const uint64_t cellId = it->getCellId();
      mpiGrid[cellId]->parameters[CellParams::ISCELLSAVINGF] = 0.0;
      // CellID stride selection
      if (P::systemWriteDistributionWriteStride[index] > 0 &&
          cellId % P::systemWriteDistributionWriteStride[index] == 0) {
         //This zone writes out its velocity space:
         //O: FIX THIS!
         velSpaceZones.push_back( &(*it) );
         mpiGrid[cellId]->parameters[CellParams::ISCELLSAVINGF] = 1.0;
         continue; // Avoid double entries in case the cell also matches following conditions.
      }
      // Cell lines selection
      // Determine cellID's 3D indices
      lineX =  cellId % P::xcells_ini;
      lineY = (cellId / P::xcells_ini) % P::ycells_ini;
      lineZ = (cellId /(P::xcells_ini *  P::ycells_ini)) % P::zcells_ini;
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
         //This zone writes out its velocity space so add it to the std::vector
         //O: FIX THIS!
         velSpaceZones.push_back( &(*it) );
         mpiGrid[cellId]->parameters[CellParams::ISCELLSAVINGF] = 1.0;
      }
   }

   // Write variables calculate d by DataReductionOperators (DRO). We do not know how many 
   // numbers each DRO calculates, so a buffer has to be re-allocated for each DRO:
   for (uint i = 0; i < dataReducer.size(); ++i) {
      o_writeDataReducer(mpiGrid, local_zones, dataReducer, i, vlsvWriter);
   }

   uint64_t numVelSpaceZones;
   uint64_t localNumVelSpaceZones;
   localNumVelSpaceZones=velSpaceZones.size();
   MPI_Allreduce(&localNumVelSpaceZones,&numVelSpaceZones,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
   if(numVelSpaceZones>0) {
      //write out velocity space data, if there are cells with this data
      o_writeVelocityDistributionData(vlsvWriter,mpiGrid,velSpaceZones,MPI_COMM_WORLD);
   }


   return true;
}

template <typename T>
bool o_writeScalarParameter(string name,T value,Writer& vlsvWriter,int masterRank,MPI_Comm comm){
//bool Writer::writeParameter(const std::string& parameterName,const T* const array)
   int myRank;
   MPI_Comm_rank(comm,&myRank);
   if(myRank==masterRank){
      map<string,string> attributes;
      std::ostringstream s;
      s << value;
      attributes["value"]=s.str();
      vlsvWriter.writeArray("PARAMETER", attributes, 1, 1, &value);
   }
   return true;
}

bool o_writeCommonGridData(
   Writer& vlsvWriter,
   const dccrg::Dccrg<SpatialCell>& mpiGrid,
   vector<uint64_t> &cells,
   const uint& index,
   MPI_Comm comm
) {
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   const int masterProcessId = 0;
   o_writeScalarParameter("t",P::t,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("dt",P::dt,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("tstep",P::tstep,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("fileIndex",index,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("xmin",P::xmin,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("xmax",P::xmax,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("ymin",P::ymin,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("ymax",P::ymax,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("zmin",P::zmin,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("zmax",P::zmax,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("xcells_ini",P::xcells_ini,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("ycells_ini",P::ycells_ini,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("zcells_ini",P::zcells_ini,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("vxmin",P::vxmin,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("vxmax",P::vxmax,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("vymin",P::vymin,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("vymax",P::vymax,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("vzmin",P::vzmin,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("vzmax",P::vzmax,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("vxblocks_ini",P::vxblocks_ini,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("vyblocks_ini",P::vyblocks_ini,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   o_writeScalarParameter("vzblocks_ini",P::vzblocks_ini,vlsvWriter,masterProcessId,MPI_COMM_WORLD);
   return true; //to make compiler happy,no real errorchecking done
}


//O: Ok this function name is ridiculously long. I just lack the imagination for a proper name.
//O: Not done

bool o_writeGhostZoneDomainAndLocalIdNumbers( Writer & vlsvWriter,
                                              const string & meshName,
                                              const vector<Zone> & ghost_zones ) {
   //Make sure ghost_zones is not empty:
   if( ghost_zones.empty() ) {
      cerr << "Passed an empty vector at: " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //We need to keep track of how many ghost domain ids and ghost local ids there are:
   unsigned int numberOfGhosts = 0;  

   //Declare vectors for storing data
   vector<unsigned int> ghostDomainIds;
   vector<unsigned int> ghostLocalIds;

   //Iterate through all ghost zones:
   vector<Zone>::const_iterator it;
   for( it = ghost_zones.begin(); it != ghost_zones.end(); ++it ) {
      //Domain id is the MPI process rank owning the ghost zone
      //O: Find out about the mpi process rank
      const int domainId = it->getRank();
      //get the local id
      //O: How to define?
      const int localId = it->getLocalId();
      //Append to the vectors
      ghostDomainIds.push_back( domainId );
      ghostLocalIds.push_back( localId );

      //Increment the number of ghosts as ghostDomainIds has been pushed back
      //O: Couldn't this be simply fetched with ghostDomainIds.size() as well?
      numberOfGhosts++;
   }

   //Write the array:
   map<string, string> xmlAttributes;
   //Note: should be "SpatialGrid"
   xmlAttributes["mesh"] = meshName;
   if( vlsvWriter.writeArray( "MESH_GHOST_DOMAINS", xmlAttributes, numberOfGhosts, 1, &(ghostDomainIds[0]) ) == false ) {
      return false;
   }
   if( vlsvWriter.writeArray( "MESH_DOMAIN_SIZES", xmlAttributes, numberOfGhosts, 1, &(ghostLocalIds[0]) ) == false ) {
      return false;
   }
   //Everything good
   return true;
}



bool o_writeDomainSizes( Writer & vlsvWriter,
                         const string & meshName,
                         const unsigned int numberOfLocalZones,
                         const unsigned int numberOfGhostZones ) {
   //Declare domainSize. There are two types of domain sizes -- ghost and local
   const int numberOfDomainTypes = 2;
   int domainSize[numberOfDomainTypes];
   domainSize[0] = numberOfLocalZones;
   domainSize[1] = numberOfGhostZones;

   //Write the array:
   map<string, string> xmlAttributes;
   //Put the meshName
   xmlAttributes["mesh"] = meshName;
   //Write (writeArray does the writing) Note: Here the important part is "MESH_DOMAIN_SIZES" -- visit plugin needs this
   if( vlsvWriter.writeArray( "MESH_DOMAIN_SIZES", xmlAttributes, 1, 2, domainSize ) == false ) {
      return false;
   } else {
      //write successful, return true
      return true;
   }
}



//Writes the global id numbers. Note: vector< array<unsigned int, 3> > & local_zones is a vector which has coordinates for local zones
//NOTE: local_zones and ghost_zones are edited!
bool o_writeZoneGlobalIdNumbers( Writer & vlsvWriter,
                                 const string & meshName,
                                 vector<Zone> & local_zones,
                                 vector<Zone> & ghost_zones ) {
   //Get the cells in x, y, z direction right off the bat (for the sake of clarity):
   const unsigned int xCells = P::xcells_ini;
   const unsigned int yCells = P::ycells_ini;
   const unsigned int zCells = P::zcells_ini;

   vector<uint64_t> globalIds;

   //Iterate through local_zones and store the values into globalIDs
   //Note: globalID is defined as follows: global ID = z*yCells*xCells + y*xCells + x
   vector<Zone>::iterator it;
   for( it = local_zones.begin(); it != local_zones.end(); ++it ) {
      /* O: THIS IS NOT NEEDED
      const array<unsigned int, 3> & indices = it->getIndices();
      //Get the coordinates for the sake of clarity:
      //Note: These are cell coordinates (integers)
      const unsigned int x = indices[0];
      const unsigned int y = indices[1];
      const unsigned int z = indices[2];
      */

      //Get the global id:
      const uint64_t globalId = it->getCellId();

      //Add the global id:
      globalIds.push_back( globalId );
   }
   //Do the same for ghost zones: (Append to the end of the list of global ids)
   for( it = ghost_zones.begin(); it != ghost_zones.end(); ++it ) {
      const array<unsigned int, 3> & indices = it->getIndices();
      //Get the coordinates for the sake of clarity:
      //Note: These are cell coordinates (integers)
      const unsigned int x = indices[0];
      const unsigned int y = indices[1];
      const unsigned int z = indices[2];

      //Calculate the global id:
      const unsigned int globalId = z * yCells * xCells + y * xCells + x;

      //Add the global id:
      globalIds.push_back( globalId );
   }

   //Get the total number of zones:
   const unsigned int numberOfZones = globalIds.size();

   //Write the array:
   map<string, string> xmlAttributes;
   //The name of the mesh (user input -- should be "SpatialGrid")
   xmlAttributes["name"] = meshName;
   //A mandatory 'type' -- just something visit hopefully understands, because I dont :)
   xmlAttributes["type"] = "multi_ucd";
   //Write:
   //O: NOTE: FIND OUT WHY ITS globalIds[0]!!!!!!!!!!!!!!!!!!!!!!!!!!
   if( vlsvWriter.writeArray( "MESH", xmlAttributes, numberOfZones, 1, &(globalIds[0]) ) == false ) {
      cout << "Unsuccessful writing of MESH at: " << __FILE__ << " " << __LINE__ << endl;
      return false;
   } else {
      //Successfully wrote the array
      return true;
   }
}

//Writes the cell coordinates
//Note: meshName should probably be spatialGrid
//_cells contains number of cells in x, y and z direction
bool o_writeBoundingBoxNodeCoordinates ( Writer & vlsvWriter,
                                         const string & meshName,
                                         const int masterRank,
                                         MPI_Comm comm ) {
   //Create variables xCells, yCells, zCells which tell the number of zones in the given direction
   //Note: This is for the sake of clarity.
   const unsigned int xCells = P::xcells_ini;
   const unsigned int yCells = P::ycells_ini;
   const unsigned int zCells = P::zcells_ini;

   //Create variables xmin, ymin, zmin for calculations
   const float xmin = (float)P::xmin;
   const float ymin = (float)P::ymin;
   const float zmin = (float)P::zmin;

   //Create variables for cell lengths in x, y, z directions for calculations
   const float xCellLength = (float)P::dx_ini;
   const float yCellLength = (float)P::dy_ini;
   const float zCellLength = (float)P::dz_ini;
   

   //Create node coordinates:
   //These are the coordinates for any given node in x y or z direction
   //Note: Why is it xCells + 1 and not xCells?
   float * xNodeCoordinates = new float[xCells + 1];
   float * yNodeCoordinates = new float[yCells + 1];
   float * zNodeCoordinates = new float[zCells + 1];

   //Input the coordinates for the nodes:
   //O: NOTE! I am not actually sure whether these are supposed to be the cell coordinates(integers) or the 'real' coordinates(doubles)
   for( unsigned int i = 0; i < xCells + 1; ++i ) {
      //The x coordinate of the first node should be xmin, the second xmin + xCellLength and so on
      xNodeCoordinates[i] = xmin + xCellLength * i;
   }
   for( unsigned int i = 0; i < yCells + 1; ++i ) {
      yNodeCoordinates[i] = ymin + yCellLength * i;
   }
   for( unsigned int i = 0; i < zCells + 1; ++i ) {
      zNodeCoordinates[i] = zmin + zCellLength * i;
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
   //Check the rank:
   if( myRank == masterRank ) {
      //Save with the correct name "MESH_NODE_CRDS_X" -- writeArray returns false if something goes wrong
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_X", xmlAttributes, xCells + 1, 1, xNodeCoordinates) == false ) success = false;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Y", xmlAttributes, yCells + 1, 1, yNodeCoordinates) == false ) success = false;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Z", xmlAttributes, zCells + 1, 1, zNodeCoordinates) == false ) success = false;
   } else {
      //O: Wait what's happening here? 0 here and xCells in the other?
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_X", xmlAttributes, 0, 1, xNodeCoordinates) == false ) success = false;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Y", xmlAttributes, 0, 1, yNodeCoordinates) == false ) success = false;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Z", xmlAttributes, 0, 1, zNodeCoordinates) == false ) success = false;
   }
   return success;
}


//Writes the boundaries for the grid. box_values should be a vector of size 6 and should contain:
//num of blocks in x, y, z direction, number of cells in x, y, z direction
//O: meshName should probably be "SpatialGrid"
//O: Should this do anything if it isn't a master process?
bool o_writeMeshBoundingBox( Writer & vlsvWriter, 
                             const string & meshName, 
                             const array<long unsigned int, 6> & box_values,
                             const int masterRank,
                             MPI_Comm comm ) {
   const unsigned int box_size = 6;

   //Check for correct input
   if( box_values.size() != box_size ) {
      cerr << "Box_values must be of size 6 at: " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }

   //Get my rank from the MPI_Comm
   int myRank;
   MPI_Comm_rank(comm, &myRank);

   //Declare boundaryBox (writeArray expects it like this but it will have the same values as box_values)
   long unsigned int boundaryBox[box_size];
   //Save the box_values into boundaryBox:
   array<long unsigned int, 6>::const_iterator it = box_values.begin();
   for( unsigned int i = 0; i < box_size; ++it, ++i ) {
      //Save the value to bbox:
      boundaryBox[i] = *it;
   }

   //Write the array:
   //Declare attributes
   map<string, string> xmlAttributes;
   //We received mesh name as a parameter: MOST LIKELY THIS IS SpatialGrid!
   xmlAttributes["mesh"] = meshName;

   //Write an array (NOTE: success will be returned and writeArray will return true or false depending on whether or not the write is successful)
   bool success;
   if( myRank == masterRank ) {
      //The visit plugin expects MESH_BBOX as a keyword
      //NOTE: writeArray writes boundaryBox
      success = vlsvWriter.writeArray("MESH_BBOX", xmlAttributes, 6, 1, boundaryBox);
   } else {
      success = vlsvWriter.writeArray("MESH_BBOX", xmlAttributes, 0, 1, boundaryBox);
   }
   return success;
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

//O: Added bool new here which is false by default. bool new determines whether we want to use the new vlsv library or the old one
bool writeGrid(
               const dccrg::Dccrg<SpatialCell>& mpiGrid,
               DataReducer& dataReducer,
               const uint& index,
               bool newLib ) {
   if( newLib ) {
      //Go with the new vlsv library:
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

      //Open the file with vlsvWriter:
      Writer vlsvWriter;
      const int masterProcessId = 0;
      vlsvWriter.open( fname.str(), MPI_COMM_WORLD, masterProcessId );

      //Declare local zones:
      vector<Zone> local_zones;
      // Get all local cell Ids 
      vector<uint64_t> local_cells = mpiGrid.get_cells();
      //no order assumed so let's order cells here
      //O: NOTE: Should this be done? It might be better to input the local cell ids first -- or not. Probably better to sort first
      std::sort(local_cells.begin(), local_cells.end());

      //Create the local zones (inputs proper values into local_zones):
      //O: THIS IS INCORRECT! createZone should create both local zone and ghost zone at the same time
      //if( createZone( mpiGrid, local_cells, local_zones ) == false ) return false;

      //Declare ghost zones:
      vector<Zone> ghost_zones;
      // Get all ghost cell Ids
      vector<uint64_t> ghost_cells = mpiGrid.get_remote_cells_on_process_boundary( NEAREST_NEIGHBORHOOD_ID );
      //Sort the cells:
      std::sort(ghost_cells.begin(), ghost_cells.end());
      //Create the zones (inputs proper values into local_zones and ghost_zones):
      if( createZone( mpiGrid, local_cells, ghost_cells, local_zones, ghost_zones ) == false ) return false;

      //The mesh name is "SpatialGrid"
      const string meshName = "SpatialGrid";

      //Write mesh boundaries: NOTE: master process only
      //Note: box_values[0] is the number of cells in x direction, [1] is in y direction, [2] is in z direction
      //[3], [4], [5] are 1 because we are not using a block-based mesh
      const long unsigned int notBlockBased = 1;
      array<long unsigned int, 6> box_values{ { P::xcells_ini, P::ycells_ini, P::zcells_ini,
                                                 notBlockBased, notBlockBased, notBlockBased } };
      if( o_writeMeshBoundingBox( vlsvWriter, meshName, box_values, masterProcessId, MPI_COMM_WORLD ) == false ) return false;

      //Write the node coordinates: NOTE: master process only
      if( o_writeBoundingBoxNodeCoordinates( vlsvWriter, meshName, masterProcessId, MPI_COMM_WORLD ) == false ) return false;

      //Write basic grid parameters: NOTE: master process only ( I think )
      //O: Not sure if the VisIt plugin understands this
      if( o_writeCommonGridData(vlsvWriter, mpiGrid, local_cells, P::systemWrites[index], MPI_COMM_WORLD) == false ) return false;

      //Write zone global id numbers:
      if( o_writeZoneGlobalIdNumbers( vlsvWriter, meshName, local_zones, ghost_zones ) == false ) return false;

      //Write domain sizes:
      if( o_writeDomainSizes( vlsvWriter, meshName, local_zones.size(), ghost_zones.size() ) == false ) return false;

      //Write ghost zone domain and local id numbers
      if( o_writeGhostZoneDomainAndLocalIdNumbers( vlsvWriter, meshName, ghost_zones ) == false ) return false;

      //Write necessary variables:
      if( o_writeVariables( mpiGrid, dataReducer, vlsvWriter, index, meshName, local_zones ) == false ) return false;
      
      phiprof::initializeTimer("Barrier","MPI","Barrier");
      phiprof::start("Barrier");
      MPI_Barrier(MPI_COMM_WORLD);
      phiprof::stop("Barrier");
      vlsvWriter.close();
      phiprof::stop("writeGrid-reduced");

      return success;

   } else {
      //Go with the old vlsv library
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
      vlsvWriter.open(fname.str(),MPI_COMM_WORLD,0,0);
   
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
   return false;
}
   



bool writeRestart(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                  DataReducer& dataReducer,
                  const string& name,
                  const uint& index,
                  const int& stripe,
                  bool newLib) {
   if( newLib ) {
      //Go with the new vlsv library
   } else {
      //Go with the old vlsv library
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
      vlsvWriter.open(fname.str(),MPI_COMM_WORLD,0,stripe); 
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
   return false;
}
   



bool writeDiagnostic(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                     DataReducer& dataReducer,
                     bool newLib)
{
   if( newLib ) {
      //Go with the new vlsv library
   } else {
      //Go with the old vlsv library
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
   return false;
}



