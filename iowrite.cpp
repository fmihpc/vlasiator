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


class Zone_base {
   protected:
      //Points to some given mpiGrid
      const dccrg::Dccrg<SpatialCell> * mpiGrid;
      //Points to some given cell in the mpiGrid
      const SpatialCell * cellInfo;
      //Cell id of the zone
      //NOTE: Zone is basically equivalent to a cell
      //NOTE: This cell id is following the logic used for grids in visit
      //so it's starting from 0 as opposed to 1 as is the case in vlasiator code
      uint64_t cellId;
   public:
      const uint64_t & getCellId() const {
         return cellId;
      }
      const uint64_t & getLocalId() const {
         //Return the zone's local id:
         return cellInfo->ioLocalCellId;
      }
      const SpatialCell * getSpatialCell() const {
         return cellInfo;
      }
      inline int getRank() const {
         //Return the zone's rank:
         return mpiGrid->get_process( cellId + 1 );
      }
      Zone_base( const dccrg::Dccrg<SpatialCell> & _mpiGrid, const uint64_t & _cellId ) {
         //Point the mpiGrid in the right direction:
         mpiGrid = &_mpiGrid;
         //Get the cell info
         cellInfo = _mpiGrid[_cellId];
         //Note: Visit is following a grid logic where cell ids start from number 0 and move up
         //Vlasiator's grid starts from cell id 1. What we want is the visit's logic
         cellId = _cellId - 1;
      }
      virtual ~Zone_base() {
      }
};

//A version of Zone_base with error checks:
class Zone_safe : public Zone_base {
   private:
      //A few checks to make sure everything is ok.
      //There might be a chance the cellInfo is pointing to null or the given cell id is out of bounds
      //So let the user decide if he wants to check for it:
      bool cellIdIsOk() const {
         //Make sure the cell id is not out of bounds:
         const uint64_t & xCells = P::xcells_ini; //Number of cells in x-direction in our grid
         const uint64_t & yCells = P::ycells_ini;
         const uint64_t & zCells = P::zcells_ini;
         const uint64_t maxXCoordinate = xCells - 1; //The largest coordinate number for x in our grid
         const uint64_t maxYCoordinate = yCells - 1;
         const uint64_t maxZCoordinate = zCells - 1;

         //Calculate the largest possible cell id:
         const uint64_t maxCellId =  maxZCoordinate * yCells * xCells + maxYCoordinate * xCells + maxXCoordinate;
         if( cellId > maxCellId ) return false; //If the cell id is out of bounds, return false
         return true; //Everything ok
      }
      bool cellInfoIsOk() const {
         //Checks if cellInfo is pointing to anything:
         if( cellInfo == NULL ) return false;
         return true; //Everything ok
      }
      bool cellGridIsOk() const {
         //Checks if mpiGrid is pointing to anything:
         if( mpiGrid == NULL ) return false;
         return true; //Everything ok
      }
      bool isOk() const {
         //Check to make sure the cell id is ok:
         if( !cellIdIsOk() ) return false;
         //Check to make sure cellInfo is pointing to something
         if( !cellInfoIsOk() ) return false;
         //Check to make sure mpiGrid is pointing to something
         if( !cellGridIsOk() ) return false;
         return true; //Everything ok
      }
   public:
      const uint64_t & getLocalId() const {
         //Do error checking
         if( !cellInfoIsOk() ) {
            cerr << "ERROR at: " << __FILE__ << " " << __LINE__ << ", invalid cell id!" << endl;
         }
         //Return the same as Zone_base would return:
         return Zone_base::getLocalId();
      }
      inline int getRank() const {
         if( !cellGridIsOk() ) {
            cerr << "ERROR at: " << __FILE__ << " " << __LINE__ << ", NULL GRID!" << endl;
         }
         return Zone_base::getRank();
      }
      const SpatialCell * getSpatialCell() const {
         if( !cellInfoIsOk() ) {
            cerr << "ERROR at: " << __FILE__ << " " << __LINE__ << ", invalid cell id!" << endl;
         }
         return Zone_base::getSpatialCell();
      }
      Zone_safe( const dccrg::Dccrg<SpatialCell> & _mpiGrid, const uint64_t & _cellId )
               : Zone_base( _mpiGrid, _cellId ) {
         if( !cellIdIsOk() ) {
            cerr << "ERROR at: " << __FILE__ << " " << __LINE__ << ", INVALID CELL ID " << cellId << endl;
         }
         if( !cellInfoIsOk() ) {
            cerr << "ERROR at: " << __FILE__ << " " << __LINE__ << ", INVALID CELL INFO" << endl;
         }
         if( !cellGridIsOk() ) {
            cerr << "ERROR at: " << __FILE__ << " " << __LINE__ << ", INVALID CELL GRID" << endl;
         }
         if( _cellId == 0 ) {
            cerr << "ERROR at: " << __FILE__ << " " << __LINE__ << ", INVALID CELL ID INPUT!" << endl;
         }
      }
      ~Zone_safe() {
      }
};

//We want to use Zone_base:
//typedef Zone_base Zone;
typedef Zone_safe Zone;

bool createZone(  dccrg::Dccrg<SpatialCell> & mpiGrid,
                  const vector<uint64_t> & local_cells,
                  const vector<uint64_t> & ghost_cells,
                  vector<Zone> & local_zones,
                  vector<Zone> & ghost_zones,
                  MPI_Comm comm ) {
   if( local_cells.empty() ) {
      if( !ghost_cells.empty() ) {
         //Local cells empty but ghost cells not empty -- something very wrong
         cerr << "ERROR! LOCAL CELLS EMPTY BUT GHOST CELLS NOT AT: " << __FILE__ << " " << __LINE__ << endl;
         exit( 1 );
      }
      cerr << "Warning, inputting empty local cells at: " << __FILE__ << " " << __LINE__ << endl;;
   }
   int myRank;
   MPI_Comm_rank(comm,&myRank);

   //Get the cell boundaries (Tells how many cells there are in x, y, z direction)
   //Note: P::xcells_ini stands for the number of cells in x direction
   const array<uint64_t, 3> cell_bounds { { P::xcells_ini, P::ycells_ini, P::zcells_ini } };

   //Declare an iterator for iterating though the cell ids
   vector<uint64_t>::const_iterator it;
   //Local ids for the process start from 0 (this is used in the iteration)
   uint64_t thisProcessLocalId = 0;
   //Iterate through local cells and save the data into zones
   for( it = local_cells.begin(); it != local_cells.end(); ++it ) {
      //NOTE: (*it) = cellId
      //Set the local id
      mpiGrid[(*it)]->ioLocalCellId = thisProcessLocalId;
      //Append _zone to the list of zones: (Note: *it = cell id)
      local_zones.push_back( Zone( mpiGrid, (*it) ) );
      //Increment the local id
      thisProcessLocalId++;
   }
   //Update the local ids:
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_IOLOCALCELLID);
   mpiGrid.update_remote_neighbor_data(NEAREST_NEIGHBORHOOD_ID);

   //Do the same for ghost_zones:
   for( it = ghost_cells.begin(); it != ghost_cells.end(); ++it ) {
      //Append zone to the list of zones:
      ghost_zones.push_back( Zone( mpiGrid, (*it) ) );
   }
   if( local_zones.empty() ) {
      if( !ghost_zones.empty() ) {
         //Local cells empty but ghost cells not empty -- something very wrong
         cerr << "ERROR! LOCAL ZONES EMPTY BUT GHOST ZONES NOT AT: " << __FILE__ << " " << __LINE__ << endl;
         exit(1);
      }
      cerr << "Warning, outputting empty local zones at: " << __FILE__ << " " << __LINE__ << endl;;
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
                                    dccrg::Dccrg<SpatialCell>& mpiGrid,
                                    const vector<uint64_t> & cells,
                                    MPI_Comm comm
                                    ) {
   // O: NOTE: I'm not entirely sure this works with the VisIt plugin
   // Write velocity blocks and related data. 
   // In restart we just write velocity grids for all cells.
   // First write global Ids of those cells which write velocity blocks (here: all cells):
   map<string,string> attribs;
   bool success=true;

   //Make sure zones is not empty:
   if( cells.empty() ) {
      cerr << "Warning, passed empty zones vector at: " << __FILE__ << " " << __LINE__ << endl;
      //return false
   }

   //Compute totalBlocks
   uint64_t totalBlocks = 0;  
   vector<uint> blocksPerCell;   
   for(size_t cell=0;cell<cells.size();++cell){
      totalBlocks+=mpiGrid[cells[cell]]->number_of_blocks;
      blocksPerCell.push_back(mpiGrid[cells[cell]]->number_of_blocks);
   }
   
   //The name of the mesh is "SpatialGrid"
   attribs["mesh"] = "SpatialGrid";
   const unsigned int vectorSize = 1;
   //Write the array:
   if (vlsvWriter.writeArray("CELLSWITHBLOCKS",attribs,cells.size(),vectorSize,cells.data()) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;
   //Write blocks per cell, this has to be in the same order as cellswitblocks so that extracting works
   if(vlsvWriter.writeArray("BLOCKSPERCELL",attribs,blocksPerCell.size(),vectorSize,blocksPerCell.data()) == false) success = false;
   if (success == false) logFile << "(MAIN) writeGrid: ERROR failed to write CELLSWITHBLOCKS to file!" << endl << writeVerbose;
   cerr << "WROTE BLOCKSPERCELL" << endl;

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
      cerr << "FAILED TO WRITE VELOCITY BLOCK COORDINATES AT: " << __FILE__ << " " << __LINE__ << endl;
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
   //O: FIX! MAKE SURE THIS IS CORRECT!
   attribs["mesh"] = "SpatialGrid";
   attribs["name"] = "avgs";
   //O: note: This could be done with &(velocityBlockData[0]), too
   if (vlsvWriter.writeArray("BLOCKVARIABLE",attribs,totalBlocks,SIZE_VELBLOCK,velocityBlockData.data()) == false) success=false;
   if (success ==false)      logFile << "(MAIN) writeGrid: ERROR occurred when writing BLOCKVARIABLE f" << endl << writeVerbose;
   velocityBlockData.clear();
    
   return success;
}

bool o_writeDataReducer(const dccrg::Dccrg<SpatialCell>& mpiGrid,
                      const vector<uint64_t>& cells,
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
  
   
   uint64_t arraySize = cells.size()*vectorSize*dataSize;
   
   //Request DataReductionOperator to calculate the reduced data for all local cells:
   char* varBuffer = NULL;
   try {
      varBuffer = new char[arraySize];
   } catch( bad_alloc& ) {
      cerr << "ERROR, FAILED TO ALLOCATE MEMORY AT: " << __FILE__ << " " << __LINE__ << endl;
      exit( 1 );
   }
   for (uint64_t cell=0; cell<cells.size(); ++cell) {
      //Reduce data ( return false if the operation fails )
      if (dataReducer.reduceData(mpiGrid[cells[cell]],dataReducerIndex,varBuffer + cell*vectorSize*dataSize) == false){
         success = false;
         logFile << "(MAIN) writeGrid: ERROR datareductionoperator '" << dataReducer.getName(dataReducerIndex) <<
            "' returned false!" << endl << writeVerbose;
      }
   }
   // Write  reduced data to file if DROP was successful:
   if(success){
      //Write the array:
      if (vlsvWriter.writeArray("VARIABLE",attribs, dataType, cells.size(), vectorSize, dataSize, varBuffer) == false) {
         success = false;
         logFile << "(MAIN) writeGrid: ERROR failed to write datareductionoperator data to file!" << endl << writeVerbose;
      }
   }
   delete[] varBuffer;
   varBuffer = NULL;
   return success;
}



bool o_writeVariables( dccrg::Dccrg<SpatialCell>& mpiGrid,
                       DataReducer& dataReducer,
                       Writer & vlsvWriter, 
                       const vector<uint64_t> & cells ) {
   //Make sure local_zones is not empty, or at least warn if it is:
   if ( cells.empty() ) {
      cerr << "WARNING: passed an empty cells vector at: " << __FILE__ << " " << __LINE__ << endl;
   }

   // Write variables calculate d by DataReductionOperators (DRO). We do not know how many 
   // numbers each DRO calculates, so a buffer has to be re-allocated for each DRO:
   for (uint i = 0; i < dataReducer.size(); ++i) {
      o_writeDataReducer(mpiGrid, cells, dataReducer, i, vlsvWriter);
   }
   return true;
}


//template <typename T>
//bool o_writeScalarParameter(string name,T value,Writer& vlsvWriter,int masterRank,MPI_Comm comm){
//   int myRank;
//   MPI_Comm_rank(comm,&myRank);
//   unsigned int arraySize;
//   unsigned int vectorSize;
//   map<string,string> xmlAttributes;
//   if( myRank==masterRank ) {
//      //Only master rank writes out the data:
//      arraySize = 1;
//      vectorSize = 1;
//      std::ostringstream s;
//      s << value;
//      xmlAttributes["value"]=s.str();
//   } else {
//      //Write dummy data:
//      arraySize = 0;
//      vectorSize = 1;
//      xmlAttributes["value"] = "";
//   }
//   xmlAttributes["name"] = name;
//   //Write the data:
//   vlsvWriter.writeArray("PARAMETER", xmlAttributes, arraySize, vectorSize, &value);
//   return true;
//}


bool o_writeCommonGridData(
   Writer& vlsvWriter,
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   const uint& index,
   MPI_Comm comm
) {
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   const int masterProcessId = 0;
   vlsvWriter.writeParameter("t", &P::t);
   vlsvWriter.writeParameter("dt", &P::dt);
   vlsvWriter.writeParameter("tstep", &P::tstep);
   vlsvWriter.writeParameter("xmin", &P::xmin);
   vlsvWriter.writeParameter("xmax", &P::xmax);
   vlsvWriter.writeParameter("ymin", &P::ymin);
   vlsvWriter.writeParameter("ymax", &P::ymax);
   vlsvWriter.writeParameter("zmin", &P::zmin);
   vlsvWriter.writeParameter("zmax", &P::zmax);
   vlsvWriter.writeParameter("xcells_ini", &P::xcells_ini);
   vlsvWriter.writeParameter("ycells_ini", &P::ycells_ini);
   vlsvWriter.writeParameter("zcells_ini", &P::zcells_ini);
   vlsvWriter.writeParameter("vxmin", &P::vxmin);
   vlsvWriter.writeParameter("vxmax", &P::vxmax);
   vlsvWriter.writeParameter("vymin", &P::vymin);
   vlsvWriter.writeParameter("vymax", &P::vymax);
   vlsvWriter.writeParameter("vzmin", &P::vzmin);
   vlsvWriter.writeParameter("vzmax", &P::vzmax);
   vlsvWriter.writeParameter("vxblocks_ini", &P::vxblocks_ini);
   vlsvWriter.writeParameter("vyblocks_ini", &P::vyblocks_ini);
   vlsvWriter.writeParameter("vzblocks_ini", &P::vzblocks_ini);
   return true; //to make compiler happy,no real errorchecking done
}


bool o_writeGhostZoneDomainAndLocalIdNumbers( dccrg::Dccrg<SpatialCell>& mpiGrid,
                                              Writer & vlsvWriter,
                                              const string & meshName,
                                              const vector<Zone> & ghost_zones ) {
   //Declare vectors for storing data
   vector<uint64_t> ghostDomainIds;
   vector<uint64_t> ghostLocalIds;

   //Iterate through all ghost zones:
   vector<Zone>::const_iterator it;
   for( it = ghost_zones.begin(); it != ghost_zones.end(); ++it ) {
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

      //Append to the vectors
      ghostDomainIds.push_back( it->getRank() );
      ghostLocalIds.push_back( it->getLocalId() );
   }

   //We need the number of ghost zones for vlsvWriter:
   uint64_t numberOfGhosts = ghost_zones.size();

   //Write:
   map<string, string> xmlAttributes; //Used for writing in info
   //Note: should be "SpatialGrid"
   xmlAttributes["mesh"] = meshName;
   const unsigned int vectorSize = 1;
   //Write the in the number of ghost domains: (Returns false if writing fails)
   if( vlsvWriter.writeArray( "MESH_GHOST_DOMAINS", xmlAttributes, numberOfGhosts, vectorSize, &(ghostDomainIds[0]) ) == false ) {
      cerr << "Error, failed to write MEST_GHOST_DOMAINS at: " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Write the in the number of ghost local ids: (Returns false if writing fails)
   if( vlsvWriter.writeArray( "MESH_GHOST_LOCALIDS", xmlAttributes, numberOfGhosts, vectorSize, &(ghostLocalIds[0]) ) == false ) {
      cerr << "Error, failed to write MEST_GHOST_LOCALIDS at: " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Everything good
   return true;
}



bool o_writeDomainSizes( Writer & vlsvWriter,
                         const string & meshName,
                         const unsigned int & numberOfLocalZones,
                         const unsigned int & numberOfGhostZones ) {
   //Declare domainSize. There are two types of domain sizes -- ghost and local
   const unsigned int numberOfDomainTypes = 2;
   unsigned int domainSize[numberOfDomainTypes];
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
      return false;
   }
   return true;
}



//Writes the global id numbers. Note: vector< array<unsigned int, 3> > & local_zones is a vector which has coordinates for local zones
bool o_writeZoneGlobalIdNumbers( const dccrg::Dccrg<SpatialCell>& mpiGrid,
                                 Writer & vlsvWriter,
                                 const string & meshName,
                                 const vector<Zone> & local_zones,
                                 const vector<Zone> & ghost_zones ) {
   if( local_zones.empty() ) {
      if( !ghost_zones.empty() ) {
         //Something very wrong -- local zones should always have members when ghost zones has members
         cerr << "ERROR, LOCAL ZONES EMPTY BUT GHOST ZONES NOT AT " << __FILE__ << __LINE__ << endl;
         exit( 1 );
      }
      cerr << "Warning: Inputting empty local zones at: " << __FILE__ << " " << __LINE__ << endl;
   }

   //Get the cells in x, y, z direction right off the bat (for the sake of clarity):
   const unsigned int xCells = P::xcells_ini;
   const unsigned int yCells = P::ycells_ini;
   const unsigned int zCells = P::zcells_ini;

   vector<uint64_t> globalIds;

   //Iterate through local_zones and store the values into globalIDs
   //Note: globalID is defined as follows: global ID = z*yCells*xCells + y*xCells + x
   vector<Zone>::const_iterator it;
   for( it = local_zones.begin(); it != local_zones.end(); ++it ) {
      //Add the global id:
      globalIds.push_back( it->getCellId() );
   }
   //Do the same for ghost zones: (Append to the end of the list of global ids)
   for( it = ghost_zones.begin(); it != ghost_zones.end(); ++it ) {
      //Add the global id:
      globalIds.push_back( it->getCellId() );
   }

   //Get the total number of zones:
   const uint64_t numberOfZones = globalIds.size();

   //Write the array:
   map<string, string> xmlAttributes;
   //The name of the mesh (user input -- should be "SpatialGrid")
   xmlAttributes["name"] = meshName;
   //A mandatory 'type' -- just something visit hopefully understands, because I dont :)
   xmlAttributes["type"] = "multi_ucd";
   //Set periodicity:
   if( mpiGrid.is_periodic( 0 ) ) xmlAttributes["xperiodic"] = "yes"; else xmlAttributes["xperiodic"] = "no";
   if( mpiGrid.is_periodic( 1 ) ) xmlAttributes["yperiodic"] = "yes"; else xmlAttributes["yperiodic"] = "no";
   if( mpiGrid.is_periodic( 2 ) ) xmlAttributes["zperiodic"] = "yes"; else xmlAttributes["zperiodic"] = "no";
   //Write:
   if( numberOfZones == 0 ) {
      const uint64_t dummy_data = 0;
      const unsigned int dummy_array = 0;
      if( vlsvWriter.writeArray( "MESH", xmlAttributes, dummy_array, 1, &dummy_data ) == false ) {
         cerr << "Unsuccessful writing of MESH at: " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
   } else {
      if( vlsvWriter.writeArray( "MESH", xmlAttributes, numberOfZones, 1, &(globalIds[0]) ) == false ) {
         cerr << "Unsuccessful writing of MESH at: " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
   }
   //Successfully wrote the array
   return true;
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
   Real * xNodeCoordinates = NULL;
   Real * yNodeCoordinates = NULL;
   Real * zNodeCoordinates = NULL;
   try{
      xNodeCoordinates = new Real[xCells + 1];
      yNodeCoordinates = new Real[yCells + 1];
      zNodeCoordinates = new Real[zCells + 1];
   } catch( bad_alloc& ) {
      cerr << "ERROR, FAILED TO ALLOCATE MEMORY AT: " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }
   if( yNodeCoordinates == NULL || yNodeCoordinates == NULL || zNodeCoordinates == NULL ) {
      cerr << "ERROR, NULL POINTER AT: " << __FILE__ << " " << __LINE__ << endl;
      exit(1);
   }

   //Input the coordinates for the nodes:
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
   //Check the rank and write the arrays:
   const unsigned int vectorSize = 1;
   unsigned int arraySize;
   if( myRank == masterRank ) {
      //Save with the correct name "MESH_NODE_CRDS_X" -- writeArray returns false if something goes wrong
      arraySize = xCells + 1;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_X", xmlAttributes, arraySize, vectorSize, xNodeCoordinates) == false ) success = false;
      arraySize = yCells + 1;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Y", xmlAttributes, arraySize, vectorSize, yNodeCoordinates) == false ) success = false;
      arraySize = zCells + 1;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Z", xmlAttributes, arraySize, vectorSize, zNodeCoordinates) == false ) success = false;
   } else {
      //Not a master process, so write empty:
      arraySize = 0;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_X", xmlAttributes, arraySize, vectorSize, xNodeCoordinates) == false ) success = false;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Y", xmlAttributes, arraySize, vectorSize, yNodeCoordinates) == false ) success = false;
      if( vlsvWriter.writeArray("MESH_NODE_CRDS_Z", xmlAttributes, arraySize, vectorSize, zNodeCoordinates) == false ) success = false;
   }
   delete[] xNodeCoordinates;
   delete[] yNodeCoordinates;
   delete[] zNodeCoordinates;
   return success;
}


//Writes the boundaries for the grid. box_values should be a vector of size 6 and should contain:
//num of blocks in x, y, z direction, number of cells in x, y, z direction
//O: meshName should probably be "SpatialGrid"
//O: Should this do anything if it isn't a master process?
bool o_writeMeshBoundingBox( Writer & vlsvWriter, 
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
   long unsigned int boundaryBox[box_size] = { numberOfXCells, numberOfYCells, numberOfZCells, 
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
         return true;
      }
   }
   return false;
}


bool o_writeVelocitySpace( dccrg::Dccrg<SpatialCell>& mpiGrid,
                           Writer & vlsvWriter,
                           int index,
                           vector<uint64_t> & cells ) {
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

      uint64_t numVelSpaceCells;
      uint64_t localNumVelSpaceCells;
      localNumVelSpaceCells=velSpaceCells.size();
      MPI_Allreduce(&localNumVelSpaceCells,&numVelSpaceCells,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
      if(numVelSpaceCells>0) {
         //write out velocity space data, if there are cells with this data
         if( o_writeVelocityDistributionData( vlsvWriter, mpiGrid, velSpaceCells, MPI_COMM_WORLD ) == false ) {
            cerr << "ERROR, FAILED TO WRITE VELOCITY DISTRIBUTION DATA AT " << __FILE__ << " " << __LINE__ << endl;
         }
      }
      return true;
}


bool writeGrid(
               dccrg::Dccrg<SpatialCell>& mpiGrid,
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
      //std::sort(local_cells.begin(), local_cells.end());

      //Declare ghost zones:
      vector<Zone> ghost_zones;
      // Get all ghost cell Ids (NOTE: this works slightly differently depending on whether the grid is periodic or not)
      vector<uint64_t> ghost_cells = mpiGrid.get_remote_cells_on_process_boundary( NEAREST_NEIGHBORHOOD_ID );
      //Sort the cells:
      //std::sort(ghost_cells.begin(), ghost_cells.end());

      //Make sure the local cells and ghost cells are fetched properly
      if( local_cells.empty() ) {
         if( !ghost_cells.empty() ) {
            //Local cells empty but ghost cells not empty -- something very wrong
            cerr << "ERROR! LOCAL CELLS EMPTY BUT GHOST CELLS NOT AT: " << __FILE__ << " " << __LINE__ << endl;
            exit( 1 );
         }
         cerr << "Warning, inputting empty local cells at: " << __FILE__ << " " << __LINE__ << endl;;
      }

      //Create the zones (inputs proper values into local_zones and ghost_zones):
      //Note: there's some MPI communication between processes done in createZone
      if( createZone( mpiGrid, local_cells, ghost_cells, local_zones, ghost_zones, MPI_COMM_WORLD ) == false ) return false;     

      //The mesh name is "SpatialGrid" (This is used for writing in data)
      const string meshName = "SpatialGrid";

      //Write mesh boundaries: NOTE: master process only
      //Visit plugin needs to know the boundaries of the mesh so the number of cells in x, y, z direction
      if( o_writeMeshBoundingBox( vlsvWriter, meshName, masterProcessId, MPI_COMM_WORLD ) == false ) return false;

      //Write the node coordinates: NOTE: master process only
      if( o_writeBoundingBoxNodeCoordinates( vlsvWriter, meshName, masterProcessId, MPI_COMM_WORLD ) == false ) return false;

      //Write basic grid variables: NOTE: master process only
      if( o_writeCommonGridData(vlsvWriter, mpiGrid, P::systemWrites[index], MPI_COMM_WORLD) == false ) return false;

      //Write zone global id numbers:
      if( o_writeZoneGlobalIdNumbers( mpiGrid, vlsvWriter, meshName, local_zones, ghost_zones ) == false ) return false;

      //Write domain sizes:
      if( o_writeDomainSizes( vlsvWriter, meshName, local_zones.size(), ghost_zones.size() ) == false ) return false;

      //Write ghost zone domain and local id numbers ( VisIt plugin needs this for MPI )
      if( o_writeGhostZoneDomainAndLocalIdNumbers( mpiGrid, vlsvWriter, meshName, ghost_zones ) == false ) return false;

      //Write necessary variables:
      if( o_writeVariables( mpiGrid, dataReducer, vlsvWriter, local_cells ) == false ) return false;

      if( o_writeVelocitySpace( mpiGrid, vlsvWriter, index, local_cells ) == false ) return false;

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

      if( cells.empty() ) {
         cerr << "EMPTY CELLS ON PROCESS: " << myRank << endl;
      }

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
   



bool writeRestart(dccrg::Dccrg<SpatialCell>& mpiGrid,
                  DataReducer& dataReducer,
                  const string& name,
                  const uint& index,
                  const int& stripe,
                  bool newLib) {
   if( newLib ) {
      //Go with the new vlsv library
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

      //Declare ghost zones:
      vector<Zone> ghost_zones;
      // Get all ghost cell Ids
      //O: Make sure this actually IS getting the ghost cells!
      vector<uint64_t> ghost_cells = mpiGrid.get_remote_cells_on_process_boundary( NEAREST_NEIGHBORHOOD_ID );
      //Sort the cells:
      std::sort(ghost_cells.begin(), ghost_cells.end());

      //Create the zones (inputs proper values into local_zones and ghost_zones):
      //Note: there's some MPI communication between processes done in createZone
      if( createZone( mpiGrid, local_cells, ghost_cells, local_zones, ghost_zones, MPI_COMM_WORLD ) == false ) return false;

      //The mesh name is "SpatialGrid"
      const string meshName = "SpatialGrid";

      //Write mesh boundaries: NOTE: master process only
      //Visit plugin needs to know the boundaries of the mesh so the number of cells in x, y, z direction
      if( o_writeMeshBoundingBox( vlsvWriter, meshName, masterProcessId, MPI_COMM_WORLD ) == false ) return false;

      //Write the node coordinates: NOTE: master process only
      if( o_writeBoundingBoxNodeCoordinates( vlsvWriter, meshName, masterProcessId, MPI_COMM_WORLD ) == false ) return false;

      //Write basic grid parameters: NOTE: master process only ( I think )
      //O: Not sure if the VisIt plugin understands this
      if( o_writeCommonGridData(vlsvWriter, mpiGrid, P::systemWrites[index], MPI_COMM_WORLD) == false ) return false;

      //Write zone global id numbers:
      if( o_writeZoneGlobalIdNumbers( mpiGrid, vlsvWriter, meshName, local_zones, ghost_zones ) == false ) return false;

      //Write domain sizes:
      if( o_writeDomainSizes( vlsvWriter, meshName, local_zones.size(), ghost_zones.size() ) == false ) return false;

      //Write ghost zone domain and local id numbers
      if( o_writeGhostZoneDomainAndLocalIdNumbers( mpiGrid, vlsvWriter, meshName, ghost_zones ) == false ) return false;

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

      //Write necessary variables:
      //REMOVE: if( o_writeVariables( mpiGrid, dataReducer, vlsvWriter, index, meshName, local_zones ) == false ) return false;
      //O: NOTE: Looking at the original write restart function I think we only need o_writeDataReducer and
      //o_writeVelocityDistributionData
      for (uint i=0; i<restartReducer.size(); ++i) {
         o_writeDataReducer(mpiGrid, local_cells, restartReducer, i, vlsvWriter);
      }

      //write the velocity distribution data -- note: it's expecting a vector of pointers:
      o_writeVelocityDistributionData(vlsvWriter, mpiGrid, local_cells, MPI_COMM_WORLD);

      phiprof::stop("writeGrid-restart");//,1.0e-6*bytesWritten,"MB");

      return success;


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
      //O: I just copypasted this from the old -section. I don't think vlsvWriter is being used here.
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



