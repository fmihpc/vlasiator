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
#include <stdint.h>
#include <cmath>
#include <list>
#include <silo.h>
#include <sstream>
#include <dirent.h>
#include <mpi.h>
#include <array> //std::array from here
//#include <typeinfo>

#include "vlsv_reader.h"
#include "definitions.h"
#include "vlsvreaderinterface.h"

using namespace std;


//A struct for holding info on cell structure (the grid)
struct CellStructure {
   //The number of cells in x, y, z direction (initialized somewhere in read parameters)
   uint64_t cell_bounds[3];
   //Length of a cell in x, y, z direction
   Real cell_length[3];
   //x_min, y_min, z_min are stored here 
   Real min_coordinates[3];

   //The number of cells in x, y, z direction (initialized somewhere in read parameters)
   uint64_t vcell_bounds[3];
   //Length of a cell in x, y, z direction
   Real vblock_length[3];
   //vx_min, vy_min, vz_min are stored here
   Real min_vcoordinates[3];
};


using namespace vlsv;

//Calculates the cell coordinates and outputs into coordinates
//Input:
//[0] CellStructure cellStruct -- A struct for holding cell information. Has the cell length in x,y,z direction, for example
//[1] uint64_t cellId -- Some given cell id
//Output:
//[0] Real * coordinates -- Some coordinates x, y, z (NOTE: the vector size should be 3!)
void getCellCoordinates( const CellStructure & cellStruct, const uint64_t _cellId, array<Real, 3> & coordinates ) {
   //In vlasiator the cell ids start from 1 but for fetching coordinates it's more logical to start from 0
   const uint64_t cellId = _cellId - 1;
   //Calculate the cell coordinates in block coordinates (so in the cell grid where the coordinates are integers)
   uint64_t currentCellCoordinate[3];
   //Note: cell_bounds is a variable that tells the length of a cell in x, y or z direction (depending on the index)
   //cellStruct is a struct that holds info on the cell structure used in simulation (such as the length of the cell and the mininum
   //value of x within the cell grid)
   currentCellCoordinate[0] = cellId % cellStruct.cell_bounds[0];
   currentCellCoordinate[1] = ((cellId - currentCellCoordinate[0]) / cellStruct.cell_bounds[0]) % cellStruct.cell_bounds[1];
   currentCellCoordinate[2] = ((cellId - cellStruct.cell_bounds[0]*currentCellCoordinate[1]) / (cellStruct.cell_bounds[0]*cellStruct.cell_bounds[1]));
   //Get the coordinates of the cell. These are the coordinates in actual space (not cell coordinates, which are integers from 1 up to some number)
   coordinates[0] = cellStruct.min_coordinates[0] + currentCellCoordinate[0] * cellStruct.cell_length[0];
   coordinates[1] = cellStruct.min_coordinates[1] + currentCellCoordinate[1] * cellStruct.cell_length[1];
   coordinates[2] = cellStruct.min_coordinates[2] + currentCellCoordinate[2] * cellStruct.cell_length[2];
   //all done
   return;
}

//Initalizes cellStruct
//Input:
//[0] vlsv::Reader vlsvReader -- some reader with a file open (used for loading parameters)
//Output:
//[0] CellStructure cellStruct -- Holds info on cellStruct. The members are given the correct values here (Note: CellStructure could be made into a class
//instead of a struct with this as the constructor but since a geometry class has already been coded before, it would be a waste)
void setCellVariables( vlsvinterface::Reader & vlsvReader, CellStructure & cellStruct ) {
   //Get x_min, x_max, y_min, y_max, etc so that we know where the given cell id is in (loadParameter returns char*, hence the cast)
   //O: Note: Not actually sure if these are Real valued or not
   Real x_min, x_max, y_min, y_max, z_min, z_max, vx_min, vx_max, vy_min, vy_max, vz_min, vz_max;
   //Read in the parameter:
   if( vlsvReader.readParameter( "xmin", x_min ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "xmax", x_max ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "ymin", y_min ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "ymax", y_max ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "zmin", z_min ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "zmax", z_max ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;

   if( vlsvReader.readParameter( "vxmin", vx_min ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "vxmax", vx_max ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "vymin", vy_min ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "vymax", vy_max ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "vzmin", vz_min ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   if( vlsvReader.readParameter( "vzmax", vz_max ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;

   //Number of cells in x, y, z directions (used later for calculating where in the cell coordinates the given
   //coordinates are) (Done in getCellCoordinates)
   //There's x, y and z coordinates so the number of different coordinates is 3:
   const short int NumberOfCoordinates = 3;
   uint64_t cell_bounds[NumberOfCoordinates];
   uint64_t vcell_bounds[NumberOfCoordinates];
   //Get the number of velocity blocks in x,y,z direction from the file:
   //x-direction
   if( vlsvReader.readParameter( "vxblocks_ini", vcell_bounds[0] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   //y-direction
   if( vlsvReader.readParameter( "vyblocks_ini", vcell_bounds[1] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   //z-direction
   if( vlsvReader.readParameter( "vzblocks_ini", vcell_bounds[2] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   //Get the number of cells in x,y,z direction from the file:
   //x-direction
   if( vlsvReader.readParameter( "xcells_ini", cell_bounds[0] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   //y-direction
   if( vlsvReader.readParameter( "ycells_ini", cell_bounds[1] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   //z-direction
   if( vlsvReader.readParameter( "zcells_ini", cell_bounds[2] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
   //Now we have the needed variables, so let's calculate how much in one block equals in length:
   //Total length of x, y, z:
   Real x_length = x_max - x_min;
   Real y_length = y_max - y_min;
   Real z_length = z_max - z_min;

   Real vx_length = vx_max - vx_min;
   Real vy_length = vy_max - vy_min;
   Real vz_length = vz_max - vz_min;
   //Set the cell structure properly:
   for( int i = 0; i < NumberOfCoordinates; ++i ) {
      cellStruct.cell_bounds[i] = cell_bounds[i];
      cellStruct.vcell_bounds[i] = vcell_bounds[i];
   }
   //Calculate the cell length
   cellStruct.cell_length[0] = ( x_length / (Real)(cell_bounds[0]) );
   cellStruct.cell_length[1] = ( y_length / (Real)(cell_bounds[1]) );
   cellStruct.cell_length[2] = ( z_length / (Real)(cell_bounds[2]) );
   //Calculate the velocity cell length
   cellStruct.vblock_length[0] = ( vx_length / (Real)(vcell_bounds[0]) );
   cellStruct.vblock_length[1] = ( vy_length / (Real)(vcell_bounds[1]) );
   cellStruct.vblock_length[2] = ( vz_length / (Real)(vcell_bounds[2]) );
   //Calculate the minimum coordinates
   cellStruct.min_coordinates[0] = x_min;
   cellStruct.min_coordinates[1] = y_min;
   cellStruct.min_coordinates[2] = z_min;
   //Calculate the minimum coordinates for velocity cells
   cellStruct.min_vcoordinates[0] = vx_min;
   cellStruct.min_vcoordinates[1] = vy_min;
   cellStruct.min_vcoordinates[2] = vz_min;


   for( int i = 0; i < 3; ++i ) {
      if( cellStruct.cell_length[i] == 0 || cellStruct.cell_bounds[i] == 0 || cellStruct.vblock_length[i] == 0 || cellStruct.vcell_bounds[i] == 0 ) {
         cerr << "ERROR, ZERO CELL LENGTH OR CELL_BOUNDS AT " << __FILE__ << " " << __LINE__ << endl;
         exit(1);
      }
   }
   return;
}
//


static DBfile* fileptr = NULL; // Pointer to file opened by SILO


bool isDataTypeUint( vlsv::datatype::type& dataType ) {
   if( dataType == vlsv::datatype::type::UINT ) {
      return true;
   } else {
      return false;
   }
}

template <typename T>
uint64_t convUInt(const char* ptr,const T& dataType,const uint64_t& dataSize) {
   if ( isDataTypeUint(dataType) == false ) {
      cerr << "Erroneous datatype given to convUInt" << endl;
      exit(1);
   }
   
   switch (dataSize) {
    case 1:
      return *reinterpret_cast<const unsigned char*>(ptr);
      break;
    case 2:
      return *reinterpret_cast<const unsigned short int*>(ptr);
      break;
    case 4:
      return *reinterpret_cast<const unsigned int*>(ptr);
      break;
    case 8:
      return *reinterpret_cast<const unsigned long int*>(ptr);
      break;
   }
   return 0;
}


int SiloType(const datatype::type & dataType, const uint64_t & dataSize) {
   switch (dataType) {
      case datatype::type::INT:
         if (dataSize == 2) return DB_SHORT;
         else if (dataSize == 4) return DB_INT;
         else if (dataSize == 8) return DB_LONG;
         else return -1;
         break;
      case datatype::type::UINT:
         if (dataSize == 2) return DB_SHORT;
         else if (dataSize == 4) return DB_INT;
         else if (dataSize == 8) return DB_LONG;
         else return -1;
         break;
      case datatype::type::FLOAT:
         if (dataSize == 4) return DB_FLOAT;
         else if (dataSize == 8) return DB_DOUBLE;
         else return -1;
         break;
      case datatype::type::UNKNOWN:
         cerr << "INVALID DATATYPE AT " << __FILE__ << " " << __LINE__ << endl;
         exit(1);
   }
   return -1;
}


template<typename REAL> struct NodeCrd {
   static REAL EPS;
   REAL x;
   REAL y;
   REAL z;
   NodeCrd(const REAL& x,const REAL& y,const REAL& z): x(x),y(y),z(z) { }
   
   bool comp(const NodeCrd<REAL>& n) const {
      REAL EPS1,EPS2,EPS;
      EPS1 = 1.0e-6 * fabs(x);
      EPS2 = 1.0e-6 * fabs(n.x);
      if (x == 0.0) EPS1 = 1.0e-7;
      if (n.x == 0.0) EPS2 = 1.0e-7;
      EPS = max(EPS1,EPS2);
      if (fabs(x - n.x) > EPS) return false;
      
      EPS1 = 1.0e-6 * fabs(y);
      EPS2 = 1.0e-6 * fabs(n.y);
      if (y == 0.0) EPS1 = 1.0e-7;
      if (n.y == 0.0) EPS2 = 1.0e-7;
      EPS = max(EPS1,EPS2);
      if (fabs(y - n.y) > EPS) return false;
      
      EPS1 = 1.0e-6 * fabs(z);
      EPS2 = 1.0e-6 * fabs(n.z);
      if (z == 0.0) EPS1 = 1.0e-7;
      if (n.z == 0.0) EPS2 = 1.0e-7;
      EPS = max(EPS1,EPS2);
      if (fabs(z - n.z) > EPS) return false;
      return true;
   }
};

struct NodeComp {
    bool operator()(const NodeCrd<double>& a,const NodeCrd<double>& b) const {
      if (a.comp(b) == true) return false;
      double EPS = 0.5e-5 * (fabs(a.z) + fabs(b.z));      
      if (a.z > b.z + EPS) return false;
      if (a.z < b.z - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.y) + fabs(b.y));      
      if (a.y > b.y + EPS) return false;
      if (a.y < b.y - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.x) + fabs(b.x));      
      if (a.x > b.x + EPS) return false;
      if (a.x < b.x - EPS) return true;
      //cerr << "ERROR" << endl;
      return false;
   }
    
   bool operator()(const NodeCrd<float>& a,const NodeCrd<float>& b) const {
      if (a.comp(b) == true) return false;
      float EPS = 0.5e-5 * (fabs(a.z) + fabs(b.z));      
      if (a.z > b.z + EPS) return false;
      if (a.z < b.z - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.y) + fabs(b.y));      
      if (a.y > b.y + EPS) return false;
      if (a.y < b.y - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.x) + fabs(b.x));      
      if (a.x > b.x + EPS) return false;
      if (a.x < b.x - EPS) return true;
      //cerr << "ERROR" << endl;
      return false;
   }
};


//Function for converting a mesh variable (Saves the variable into an open SILO file)
//Input:
//[0] vlsvReader -- some vlsv reader with a file open
//[1] meshName -- name of the mesh, e.g. "SpatialGrid"
//[2] varName -- Name of the variable
template <class T>
bool convertMeshVariable(T & vlsvReader,const string& meshName,const string& varName) {
   bool success = true;

   // Writing a unstructured grid variable is a rather straightforward process. The 
   // only compilation here is that some of the variables are actually vectors, i.e. 
   // vectorSize > 1 (vectorSize == 1 for scalars). Format in which vectors are stored in VLSV 
   // differ from format in which they are written to SILO files.
   vlsv::datatype::type dataType;
   uint64_t arraySize,vectorSize,dataSize;
   list<pair<string, string> > xmlAttributes;
   xmlAttributes.push_back(make_pair("name", varName));
   xmlAttributes.push_back(make_pair("mesh", meshName));
   if (vlsvReader.getArrayInfo("VARIABLE", xmlAttributes, arraySize, vectorSize, dataType, dataSize) == false) return false;
   // Read variable data. Note that we do not actually need to care if 
   // the data is given as floats or doubles.
   char* buffer = new char[arraySize*vectorSize*dataSize];
   const short unsigned int startingIndex = 0;
   if (vlsvReader.readArray("VARIABLE", xmlAttributes, startingIndex, arraySize, buffer ) == false) success = false;
   if (success == false) {
      cerr << "FAILED TO READ VARIABLE AT " << __FILE__ << " " << __LINE__ << endl;
      delete buffer;
      return success;
   }
   
   // Vector variables need to be copied to temporary arrays before 
   // writing to SILO file:
   char** components = new char*[vectorSize];
   for (uint64_t i=0; i<vectorSize; ++i) {
      components[i] = new char[arraySize*dataSize];
      for (uint64_t j=0; j<arraySize; ++j) for (uint64_t k=0; k<dataSize; ++k) 
         components[i][j*dataSize+k] = buffer[j*vectorSize*dataSize + i*dataSize + k];
   }

   // SILO requires one variable name per (vector) component, but we only have one. 
   // That is, for electric field SILO would like to get "Ex","Ey", and "Ez", but we 
   // only have "E" in VLSV file. Use the VLSV variable name for all components.
   vector<string> varNames(vectorSize);
   vector<char*> varNamePtrs(vectorSize);
   for (uint64_t i=0; i<vectorSize; ++i) {
      stringstream ss;
      ss << varName << (i+1);
      varNames[i] = ss.str();
      varNamePtrs[i] = const_cast<char*>(varNames[i].c_str());
   }
   
   // Write the unstructured mesh variable to SILO.
   if (DBPutUcdvar(fileptr,varName.c_str(),meshName.c_str(),vectorSize,&(varNamePtrs[0]),components,arraySize,NULL,0,SiloType(dataType,dataSize),DB_ZONECENT,NULL) < 0) success = false;

   for (uint64_t i=0; i<vectorSize; ++i) {
      delete [] components[i];
   }
   delete [] components;
   delete [] buffer;   
   return success;
}




bool convertMesh(vlsvinterface::Reader & vlsvReader,const string& meshName) {
   bool success = true;
   const float EPS = 1.0e-7;
   
   // First task is to push all unique node coordinates into a map.
   // This is not too difficult for unrefined grids, since each spatial cell stores
   // its bottom lower left corner coordinate and size. For refined grid the situation
   // is more complex, as there are more unique nodes than the lower left corners:
   map<NodeCrd<Real>,uint64_t,NodeComp> nodes;
 
   //Read in all cell ids:
   vector<uint64_t> cellIds;
   if( vlsvReader.getCellIds( cellIds ) == false ) {
      cerr << "Failed to read cell ids at "  << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   //Read in cell structure (Used in calculating coordinates from cell ids)
   CellStructure cellStruct;
   setCellVariables( vlsvReader, cellStruct );

   // Read the coordinate array one node (of a spatial cell) at a time 
   // and create a map which only contains each existing node once.
   //char* coordsBuffer = new char[vectorSize*dataSize];
   //Real* ptr = reinterpret_cast<Real*>(coordsBuffer);
   //for (uint64_t i=0; i<arraySize; ++i) {
   int i = 0;
   for( vector<uint64_t>::const_iterator it = cellIds.begin(); it != cellIds.end(); ++it, ++i ) {
      //if (vlsvReader.readArray("COORDS",meshName,i,1,coordsBuffer) == false) {success = false; break;}
      
      // Insert all eight nodes of a cell into map nodes.
      // NOTE: map is a unique associative container - given a suitable comparator, map 
      // will filter out duplicate nodes.
      //Get cell coordinates:
      const uint64_t cellId = *it;
      array<Real, 3> coordinates;
      //Store the coordinates in 'coordinates' std::array
      getCellCoordinates( cellStruct, cellId, coordinates );

      //Note: cellStruct.cell_length[i] = a cell's length in i direction
      Real X0 = coordinates[0];
      Real X1 = coordinates[0]+cellStruct.cell_length[0];
      Real Y0 = coordinates[1];
      Real Y1 = coordinates[1]+cellStruct.cell_length[1];
      Real Z0 = coordinates[2];
      Real Z1 = coordinates[2]+cellStruct.cell_length[2];
      
      // Flush very small coordinate values to zero:
      if (fabs(X0) < EPS) X0 = 0.0;
      if (fabs(X1) < EPS) X1 = 0.0;
      if (fabs(Y0) < EPS) Y0 = 0.0;
      if (fabs(Y1) < EPS) Y1 = 0.0;
      if (fabs(Z0) < EPS) Z0 = 0.0;
      if (fabs(Z1) < EPS) Z1 = 0.0;

      nodes.insert(make_pair(NodeCrd<Real>(X0,Y0,Z0),0));
      nodes.insert(make_pair(NodeCrd<Real>(X1,Y0,Z0),0));
      nodes.insert(make_pair(NodeCrd<Real>(X1,Y1,Z0),0));
      nodes.insert(make_pair(NodeCrd<Real>(X0,Y1,Z0),0));
      nodes.insert(make_pair(NodeCrd<Real>(X0,Y0,Z1),0));
      nodes.insert(make_pair(NodeCrd<Real>(X1,Y0,Z1),0));
      nodes.insert(make_pair(NodeCrd<Real>(X1,Y1,Z1),0));
      nodes.insert(make_pair(NodeCrd<Real>(X0,Y1,Z1),0));
   }
   if (success == false) {
      cerr << "ERROR reading array COORDS" << endl;
   }

   // Copy unique node x,y,z coordinates into separate arrays, 
   // which will be passed to silo writer:
   uint64_t counter = 0;
   Real* xcrds = new Real[nodes.size()];
   Real* ycrds = new Real[nodes.size()];
   Real* zcrds = new Real[nodes.size()];
   for (map<NodeCrd<Real>,uint64_t>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
      it->second = counter;
      xcrds[counter] = it->first.x;
      ycrds[counter] = it->first.y;
      zcrds[counter] = it->first.z;
      ++counter;
   }

   // Read through the coordinate array again and create a node list. Each 3D spatial cell is 
   // associated with 8 nodes, and most of these nodes are shared with neighbouring cells. In 
   // order to get VisIt display the data correctly, the duplicate nodes should not be used. 
   // Here we create a list of indices into xcrds,ycrds,zcrds arrays, with eight entries per cell:
   int* nodelist = new int[8*cellIds.size()];
   i = 0;
   for( vector<uint64_t>::const_iterator it = cellIds.begin(); it != cellIds.end(); ++it, ++i ) {
      // Read the bottom lower left corner coordinates of a cell and its sizes. Note 
      // that zones will end up in SILO file in the same order as they are in VLSV file.

      //Get cell coordinates:
      const uint64_t cellId = *it;
      array<Real, 3> coordinates;
      //Store the coordinates in 'coordinates' std::array
      getCellCoordinates( cellStruct, cellId, coordinates );

      //Note: cellStruct.cell_length[i] = a cell's length in i direction
      Real X0 = coordinates[0];
      Real X1 = coordinates[0]+cellStruct.cell_length[0];
      Real Y0 = coordinates[1];
      Real Y1 = coordinates[1]+cellStruct.cell_length[1];
      Real Z0 = coordinates[2];
      Real Z1 = coordinates[2]+cellStruct.cell_length[2];
      
      // Flush very small coordinate values to zero:
      if (fabs(X0) < EPS) X0 = 0.0;
      if (fabs(X1) < EPS) X1 = 0.0;
      if (fabs(Y0) < EPS) Y0 = 0.0;
      if (fabs(Y1) < EPS) Y1 = 0.0;
      if (fabs(Z0) < EPS) Z0 = 0.0;
      if (fabs(Z1) < EPS) Z1 = 0.0;
      
      // Search the cell's nodes from the map created above. For each node in nodelist 
      // store an index into an array which only contains the unique nodes:
      map<NodeCrd<Real>,uint64_t,NodeComp>::const_iterator it2;
      it2 = nodes.find(NodeCrd<Real>(X0,Y0,Z0)); if (it2 == nodes.end()) success = false; nodelist[i*8+0] = it2->second;
      it2 = nodes.find(NodeCrd<Real>(X1,Y0,Z0)); if (it2 == nodes.end()) success = false; nodelist[i*8+1] = it2->second;
      it2 = nodes.find(NodeCrd<Real>(X1,Y1,Z0)); if (it2 == nodes.end()) success = false; nodelist[i*8+2] = it2->second;
      it2 = nodes.find(NodeCrd<Real>(X0,Y1,Z0)); if (it2 == nodes.end()) success = false; nodelist[i*8+3] = it2->second;
      it2 = nodes.find(NodeCrd<Real>(X0,Y0,Z1)); if (it2 == nodes.end()) success = false; nodelist[i*8+4] = it2->second;
      it2 = nodes.find(NodeCrd<Real>(X1,Y0,Z1)); if (it2 == nodes.end()) success = false; nodelist[i*8+5] = it2->second;
      it2 = nodes.find(NodeCrd<Real>(X1,Y1,Z1)); if (it2 == nodes.end()) success = false; nodelist[i*8+6] = it2->second;
      it2 = nodes.find(NodeCrd<Real>(X0,Y1,Z1)); if (it2 == nodes.end()) success = false; nodelist[i*8+7] = it2->second;
   }
   //O: REMOVE THIS
   //delete coordsBuffer;
   if (success == false) {
      cerr << "Failed to find node(s)" << endl;
   }
   
   // Write the unstructured mesh to SILO file:
   const int N_dims  = 3;                      // Number of dimensions
   const int N_nodes = nodes.size();           // Total number of nodes
   const int N_zones = cellIds.size();              // Total number of zones (=spatial cells)
   int shapeTypes[] = {DB_ZONETYPE_HEX};       // Hexahedrons only
   int shapeSizes[] = {8};                     // Each hexahedron has 8 nodes
   int shapeCnt[] = {N_zones};                 // Only 1 shape type (hexahedron)
   const int N_shapes = 1;                     //  -- "" --
   
   void* coords[3];                            // Pointers to coordinate arrays
   coords[0] = xcrds;
   coords[1] = ycrds;
   coords[2] = zcrds;
   
   // Write zone list into silo file:
   const string zoneListName = meshName + "Zones";
   if (DBPutZonelist2(fileptr,zoneListName.c_str(),N_zones,N_dims,nodelist,8*cellIds.size(),0,0,0,shapeTypes,shapeSizes,shapeCnt,N_shapes,NULL) < 0) success = false;
   
   // Write grid into silo file:
   if (DBPutUcdmesh(fileptr,meshName.c_str(),N_dims,NULL,coords,N_nodes,N_zones,zoneListName.c_str(),NULL,SiloType(datatype::type::FLOAT,sizeof(Real)),NULL) < 0) success = false;

   nodes.clear();
   delete nodelist;
   delete xcrds;
   delete ycrds;
   delete zcrds;

   // Write the cell IDs as a variable:
   string cellIDlabel = "Cell ID";
   DBoptlist* optList = DBMakeOptlist(1);   
   DBAddOption(optList,DBOPT_LABEL,const_cast<char*>(cellIDlabel.c_str()));
   //Note: cellIds is of type uint
   if (DBPutUcdvar1(fileptr,"CellID",meshName.c_str(),reinterpret_cast<char*>(cellIds.data()),cellIds.size(),NULL,0,SiloType(datatype::type::UINT,sizeof(uint64_t)),DB_ZONECENT,NULL) < 0) success = false;
//   delete buffer;
   DBFreeOptlist(optList);
   
   // Write all variables of this mesh into silo file:
   list<string> variables;
   vlsvReader.getVariableNames(meshName,variables);
   for (list<string>::const_iterator it=variables.begin(); it!=variables.end(); ++it) {
      if (convertMeshVariable(vlsvReader,meshName,*it) == false) success = false;
   }
   return success;
   
}


template <class T>
bool convertSILO(const string& fname) {
   bool success = true;

   
   // Open VLSV file for reading:
   T vlsvReader;
   if (vlsvReader.open(fname) == false) {
      cerr << "Failed to open '" << fname << "'" << endl;
      return false;
   }

   
   // Open SILO file for writing:
   size_t found=fname.find_last_of("/\\");
   //remove path from vlsvfile name
   string fileout =  fname.substr(found+1);
   size_t pos = fileout.rfind(".vlsv");
   if (pos != string::npos) fileout.replace(pos,5,".silo");
   
   fileptr = DBCreate(fileout.c_str(),DB_CLOBBER,DB_LOCAL,"Vlasov data file",DB_PDB);
   if (fileptr == NULL) return false;
   
   // Get the names of all meshes in vlsv file, and write into silo file:
   list<string> meshNames;
   if (vlsvReader.getMeshNames(meshNames) == false) {
      DBClose(fileptr);
      return false;
   }
   for (list<string>::const_iterator it=meshNames.begin(); it!=meshNames.end(); ++it) {
      if (convertMesh(vlsvReader,*it) == false) {
         DBClose(fileptr);
         return false;
      }
   }
   vlsvReader.close();
   DBClose(fileptr);
   return success;
}

int main(int argn,char* args[]) {
   int ntasks, rank;
   MPI_Init(&argn, &args);
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   if (rank == 0 && argn < 2) {
      cout << endl;
      cout << "USAGE: ./vlsv2vtk <input file mask(s)>" << endl;
      cout << "Each VLSV in the current directory is compared against the given file mask(s)," << endl;
      cout << "and if match is found, that file is converted into SILO format." << endl;
      cout << endl;
      return 1;
   }

   // Convert file masks into strings:
   vector<string> masks, fileList;
   for (int i=1; i<argn; ++i) masks.push_back(args[i]);

   // Compare directory contents against each mask:

   const string suffix = ".vlsv";
   int filesFound = 0, filesConverted = 0;
   for (size_t mask=0; mask<masks.size(); ++mask) {
      size_t found=masks[mask].find_last_of("/\\");
      string directory=".";
      if(found != string::npos)
         directory = masks[mask].substr(0,found);
      const string maskName =  masks[mask].substr(found+1);

      if(rank == 0) {cout << "Comparing mask '" << maskName << "' in folder '" << directory <<"'" << endl;}
      DIR* dir = opendir(directory.c_str());
      if (dir == NULL) continue;
      
      struct dirent* entry = readdir(dir);
      while (entry != NULL) {
         const string entryName = entry->d_name;
         // Compare entry name against given mask and file suffix ".vlsv":
         if (entryName.find(maskName) == string::npos || entryName.find(suffix) == string::npos) {
            entry = readdir(dir);
            continue;
         }
         
         fileList.push_back(directory);
         fileList.back().append("/");
         fileList.back().append(entryName);
         filesFound++;
         entry = readdir(dir);
      }
      closedir(dir);
      if (rank == 0 && filesFound == 0) cout << "\t no matches found" << endl;
   }
   
   for(size_t entryName = 0; entryName < fileList.size(); entryName++) {
      if(entryName%ntasks == (uint)rank) {
         cout << "\tProc " << rank << " converting '" << fileList[entryName] << "'" << endl;
         convertSILO<vlsvinterface::Reader>(fileList[entryName]);
         filesConverted++;
      }
   }
   
   int totalFilesConverted =0;   
   MPI_Reduce(&filesConverted, &totalFilesConverted, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   if (rank == 0 && totalFilesConverted == 0) cout << "\t no files converted" << endl;
   
   MPI_Finalize();
   return 0;
}
