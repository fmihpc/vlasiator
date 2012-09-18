/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>

#include <limits>
#include <stdint.h>
#include <cmath>
#include <list>
#include <silo.h>
#include <sstream>
#include <dirent.h>
#include <stdio.h>

#include "vlsvreader2.h"
#include "definitions.h"

using namespace std;

static DBfile* fileptr = NULL; // Pointer to file opened by SILO

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
     double EPS = 0.5e-3 * (fabs(a.z) + fabs(b.z));
     if (a.z > b.z + EPS) return false;
     if (a.z < b.z - EPS) return true;
     
     EPS = 0.5e-3 * (fabs(a.y) + fabs(b.y));
     if (a.y > b.y + EPS) return false;
     if (a.y < b.y - EPS) return true;
     
     EPS = 0.5e-3 * (fabs(a.x) + fabs(b.x));
     if (a.x > b.x + EPS) return false;
     if (a.x < b.x - EPS) return true;
     return false;
   }
   
   bool operator()(const NodeCrd<float>& a,const NodeCrd<float>& b) const {
      float EPS = 0.5e-3 * (fabs(a.z) + fabs(b.z));
      if (a.z > b.z + EPS) return false;
      if (a.z < b.z - EPS) return true;
      
      EPS = 0.5e-3 * (fabs(a.y) + fabs(b.y));
      if (a.y > b.y + EPS) return false;
      if (a.y < b.y - EPS) return true;
      
      EPS = 0.5e-3 * (fabs(a.x) + fabs(b.x));
      if (a.x > b.x + EPS) return false;
      if (a.x < b.x - EPS) return true;
      return false;
   }
};    

int SiloType(const VLSV::datatype& dataType,const uint64_t& dataSize) {
   switch (dataType) {
    case VLSV::INT:
      if (dataSize == 2) return DB_SHORT;
      else if (dataSize == 4) return DB_INT;
      else if (dataSize == 8) return DB_LONG;
      else return -1;
      break;
    case VLSV::UINT:
      if (dataSize == 2) return DB_SHORT;
      else if (dataSize == 4) return DB_INT;
      else if (dataSize == 8) return DB_LONG;
      else return -1;
      break;
    case VLSV::FLOAT:
      if (dataSize == 4) return DB_FLOAT;
      else if (dataSize == 8) return DB_DOUBLE;
      else return -1;
      break;
   }
   return -1;
}

uint64_t convUInt(const char* ptr,const VLSV::datatype& dataType,const uint64_t& dataSize) {
   if (dataType != VLSV::UINT) {
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

bool convertVelocityBlockVariable(VLSVReader& vlsvReader,const string& spatGridName,const string& velGridName,const uint64_t& cellID,
				  const uint64_t& N_blocks,const uint64_t& blockOffset,const string& varName) {
   bool success = true;
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name",varName));
   attribs.push_back(make_pair("mesh",spatGridName));
   
   VLSV::datatype dataType;
   uint64_t arraySize,vectorSize,dataSize;
   if (vlsvReader.getArrayInfo("BLOCKVARIABLE",attribs,arraySize,vectorSize,dataType,dataSize) == false) {
      cerr << "Could not read BLOCKVARIABLE array info" << endl;
      return false;
   }
   
   attribs.clear();
   attribs.push_back(make_pair("mesh",spatGridName));
   char* buffer = new char[N_blocks*vectorSize*dataSize];
   if (vlsvReader.readArray("BLOCKVARIABLE",varName,attribs,blockOffset,N_blocks,buffer) == false) success = false;
   if (success == false) {
      cerr << "ERROR could not read block variable" << endl;
      delete buffer;
      return success;
   }
   
   stringstream ss;
   ss << varName << cellID;
   string varName2 = ss.str();
   
   string label = "Distrib.function";
   string unit = "1/m^3 (m/s)^3";
   int conserved = 1;
   int extensive = 1;
   DBoptlist* optList = DBMakeOptlist(4);
   DBAddOption(optList,DBOPT_LABEL,const_cast<char*>(label.c_str()));
   DBAddOption(optList,DBOPT_EXTENTS_SIZE,const_cast<char*>(unit.c_str()));
   DBAddOption(optList,DBOPT_CONSERVED,&conserved);
   DBAddOption(optList,DBOPT_EXTENSIVE,&extensive);
   DBPutUcdvar1(fileptr,varName2.c_str(),velGridName.c_str(),buffer,N_blocks*vectorSize,NULL,0,SiloType(dataType,dataSize),DB_ZONECENT,optList);
   
   DBFreeOptlist(optList);
   delete buffer;
   return success;
}

bool convertVelocityBlocks2(VLSVReader& vlsvReader,const string& meshName,const uint64_t& cellID) {
   //return true;
   bool success = true;
   
   // Get some required info from VLSV file:
   // "cwb" = "cells with blocks" IDs of spatial cells which wrote their velocity grid
   // "nb"  = "number of blocks"  Number of velocity blocks in each velocity grid
   // "bc"  = "block coordinates" Coordinates of each block of each velocity grid
   list<pair<string,string> > attribs;

   
   VLSV::datatype cwb_dataType,nb_dataType,bc_dataType;
   uint64_t cwb_arraySize,cwb_vectorSize,cwb_dataSize;
   uint64_t nb_arraySize,nb_vectorSize,nb_dataSize;
   uint64_t bc_arraySize,bc_vectorSize,bc_dataSize;

   //first read in number of blocks
   attribs.clear();
   attribs.push_back(make_pair("name","Blocks"));
   attribs.push_back(make_pair("mesh",meshName));   
   if (vlsvReader.getArrayInfo("VARIABLE",attribs,nb_arraySize,nb_vectorSize,nb_dataType,nb_dataSize) == false) {
      //cerr << "Could not find array NBLOCKS" << endl;
      return false;
   }

   char* nb_buffer = new char[nb_arraySize*nb_vectorSize*nb_dataSize];
   if (vlsvReader.readArray("VARIABLE",attribs,0,nb_arraySize,nb_buffer) == false) success = false;
   if (success == false) {
      cerr << "Failed to read number of block  for mesh '" << meshName << "'" << endl;
      delete nb_buffer;
   }

   //read in other metadata
   attribs.clear();
   attribs.push_back(make_pair("name",meshName));   
   if (vlsvReader.getArrayInfo("CELLSWITHBLOCKS",attribs,cwb_arraySize,cwb_vectorSize,cwb_dataType,cwb_dataSize) == false) {
      //cerr << "Could not find array CELLSWITHBLOCKS" << endl;
      return false;
   }

   // Create buffers for cwb,nb and read data:
   char* cwb_buffer = new char[cwb_arraySize*cwb_vectorSize*cwb_dataSize];
   if (vlsvReader.readArray("CELLSWITHBLOCKS",meshName,0,cwb_arraySize,cwb_buffer) == false) success = false;
   if (success == false) {
      cerr << "Failed to read block metadata for mesh '" << meshName << "'" << endl;
      delete cwb_buffer;
      return success;
   }

   
   if (vlsvReader.getArrayInfo("BLOCKCOORDINATES",attribs,bc_arraySize,bc_vectorSize,bc_dataType,bc_dataSize) == false) {
      //cerr << "Could not find array BLOCKCOORDINATES" << endl;
      return false;
   }


   // Search for the given cellID:
   uint64_t blockOffset = 0;
   uint64_t cellIndex = numeric_limits<uint64_t>::max();
   uint64_t N_blocks;
   for (uint64_t cell=0; cell<cwb_arraySize; ++cell) {
      const uint64_t readCellID = convUInt(cwb_buffer+cell*cwb_dataSize,cwb_dataType,cwb_dataSize);
      N_blocks   = convUInt(nb_buffer+cell*nb_dataSize,nb_dataType,nb_dataSize);
      if (cellID == readCellID) {
	 cellIndex = cell;
	 break;
      }
      blockOffset += N_blocks;
   }
   if (cellIndex == numeric_limits<uint64_t>::max()) {
      //cerr << "Spatial cell #" << cellID << " not found!" << endl;
      return false;
   } else {
      //cout << "Spatial cell #" << cellID << " has offset " << blockOffset << endl;
   }
   
   map<NodeCrd<Real>,uint64_t,NodeComp> nodes;
   
   // Read all block coordinates of the velocity grid:
   char* bc_buffer = new char[N_blocks*bc_vectorSize*bc_dataSize];
   vlsvReader.readArray("BLOCKCOORDINATES",meshName,blockOffset,N_blocks,bc_buffer);
   for (uint64_t b=0; b<N_blocks; ++b) {
      Real vx_min,vy_min,vz_min,dvx,dvy,dvz;
      if (bc_dataSize == 4) {
	 vx_min = *reinterpret_cast<float*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 0*bc_dataSize );
	 vy_min = *reinterpret_cast<float*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 1*bc_dataSize );
	 vz_min = *reinterpret_cast<float*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 2*bc_dataSize );
	 dvx    = *reinterpret_cast<float*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 3*bc_dataSize );
	 dvy    = *reinterpret_cast<float*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 4*bc_dataSize );
	 dvz    = *reinterpret_cast<float*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 5*bc_dataSize );
      } else {
	 vx_min = *reinterpret_cast<double*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 0*bc_dataSize );
	 vy_min = *reinterpret_cast<double*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 1*bc_dataSize );
	 vz_min = *reinterpret_cast<double*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 2*bc_dataSize );
	 dvx    = *reinterpret_cast<double*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 3*bc_dataSize );
	 dvy    = *reinterpret_cast<double*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 4*bc_dataSize );
	 dvz    = *reinterpret_cast<double*>(bc_buffer +b*bc_vectorSize*bc_dataSize + 5*bc_dataSize );
      }
      
      creal EPS = 1.0e-7;
      for (int kv=0; kv<5; ++kv) {
	 Real VZ = vz_min + kv*dvz;
	 if (fabs(VZ) < EPS) VZ = 0.0;
	 for (int jv=0; jv<5; ++jv) {
	    Real VY = vy_min + jv*dvy;
	    if (fabs(VY) < EPS) VY = 0.0;
	    for (int iv=0; iv<5; ++iv) {
	       Real VX = vx_min + iv*dvx;
	       if (fabs(VX) < EPS) VX = 0.0;
	       nodes.insert(make_pair(NodeCrd<Real>(VX,VY,VZ),(uint64_t)0));
	    }
	 }
      }
   }
      
   Real* vx_crds = new Real[nodes.size()];
   Real* vy_crds = new Real[nodes.size()];
   Real* vz_crds = new Real[nodes.size()];
   uint64_t counter = 0;
   for (map<NodeCrd<Real>,uint64_t>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
      it->second = counter;
      vx_crds[counter] = it->first.x;
      vy_crds[counter] = it->first.y;
      vz_crds[counter] = it->first.z;
      ++counter;
   }

   const uint64_t nodeListSize = N_blocks*64*8;
   int* nodeList = new int[nodeListSize];
   for (uint64_t b=0; b<N_blocks; ++b) {
      Real vx_min,vy_min,vz_min,dvx,dvy,dvz;
      if (bc_dataSize == 4) {
	 // floats
	 vx_min = *reinterpret_cast<float*>(bc_buffer+b*bc_vectorSize*bc_dataSize+0*sizeof(float));
	 vy_min = *reinterpret_cast<float*>(bc_buffer+b*bc_vectorSize*bc_dataSize+1*sizeof(float));
	 vz_min = *reinterpret_cast<float*>(bc_buffer+b*bc_vectorSize*bc_dataSize+2*sizeof(float));
	 dvx    = *reinterpret_cast<float*>(bc_buffer+b*bc_vectorSize*bc_dataSize+3*sizeof(float));
	 dvy    = *reinterpret_cast<float*>(bc_buffer+b*bc_vectorSize*bc_dataSize+4*sizeof(float));
	 dvz    = *reinterpret_cast<float*>(bc_buffer+b*bc_vectorSize*bc_dataSize+5*sizeof(float));
      } else {
	 // doubles
	 vx_min = *reinterpret_cast<double*>(bc_buffer+b*bc_vectorSize*bc_dataSize+0*sizeof(double));
	 vy_min = *reinterpret_cast<double*>(bc_buffer+b*bc_vectorSize*bc_dataSize+1*sizeof(double));
	 vz_min = *reinterpret_cast<double*>(bc_buffer+b*bc_vectorSize*bc_dataSize+2*sizeof(double));
	 dvx    = *reinterpret_cast<double*>(bc_buffer+b*bc_vectorSize*bc_dataSize+3*sizeof(double));
	 dvy    = *reinterpret_cast<double*>(bc_buffer+b*bc_vectorSize*bc_dataSize+4*sizeof(double));
	 dvz    = *reinterpret_cast<double*>(bc_buffer+b*bc_vectorSize*bc_dataSize+5*sizeof(double));
      }
      Real VX,VY,VZ;
      creal EPS = 1.0e-7;
      map<NodeCrd<Real>,uint64_t>::const_iterator it;
      for (int kv=0; kv<4; ++kv) {
	 for (int jv=0; jv<4; ++jv) {
	    for (int iv=0; iv<4; ++iv) {
	       const unsigned int cellInd = kv*16+jv*4+iv;
	       VX = vx_min + iv*dvx; if (fabs(VX) < EPS) VX = 0.0;
	       VY = vy_min + jv*dvy; if (fabs(VY) < EPS) VY = 0.0;
	       VZ = vz_min + kv*dvz; if (fabs(VZ) < EPS) VZ = 0.0;
	       it = nodes.find(NodeCrd<Real>(VX,VY,VZ));
	       if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
	       nodeList[b*64*8+cellInd*8+0] = it->second;
	 
	       VX = vx_min + (iv+1)*dvx; if (fabs(VX) < EPS) VX = 0.0;
	       it = nodes.find(NodeCrd<Real>(VX,VY,VZ));
	       if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
	       nodeList[b*64*8+cellInd*8+1] = it->second;
	       
	       VY = vy_min + (jv+1)*dvy; if (fabs(VY) < EPS) VY = 0.0;
	       it = nodes.find(NodeCrd<Real>(VX,VY,VZ));
	       if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
	       nodeList[b*64*8+cellInd*8+2] = it->second;
	       
	       VX = vx_min + iv*dvx; if (fabs(VX) < EPS) VX = 0.0;
	       it = nodes.find(NodeCrd<Real>(VX,VY,VZ));
	       if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
	       nodeList[b*64*8+cellInd*8+3] = it->second;
	 
	       VY = vy_min + jv*dvy; if (fabs(VY) < EPS) VY = 0.0;
	       VZ = vz_min + (kv+1)*dvz; if (fabs(VZ) < EPS) VZ = 0.0;
	       it = nodes.find(NodeCrd<Real>(VX,VY,VZ));
	       if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
	       nodeList[b*64*8+cellInd*8+4] = it->second;
	       
	       VX = vx_min + (iv+1)*dvx; if (fabs(VX) < EPS) VX = 0.0;
	       it = nodes.find(NodeCrd<Real>(VX,VY,VZ));
	       if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
	       nodeList[b*64*8+cellInd*8+5] = it->second;
	       
	       VY = vy_min + (jv+1)*dvy; if (fabs(VY) < EPS) VY = 0.0;
	       it = nodes.find(NodeCrd<Real>(VX,VY,VZ));
	       if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY << ' ' << VZ << endl;
	       nodeList[b*64*8+cellInd*8+6] = it->second;
	       
	       VX = vx_min + iv*dvx; if (fabs(VX) < EPS) VX = 0.0;
	       it = nodes.find(NodeCrd<Real>(VX,VY,VZ));
	       if (it == nodes.end()) cerr << "Could not find node " << VX << ' ' << VY+dvy << ' ' << VZ+dvz << endl;
	       nodeList[b*64*8+cellInd*8+7] = it->second;
	    }
	 }
      }
   }

   const int N_dims  = 3;                      // Number of dimensions
   const int N_nodes = nodes.size();           // Total number of nodes
   const int N_zones = N_blocks*64;            // Total number of zones (=spatial cells)
   int shapeTypes[] = {DB_ZONETYPE_HEX};       // Hexahedrons only
   int shapeSizes[] = {8};                     // Each hexahedron has 8 nodes
   int shapeCnt[] = {N_zones};                 // Only 1 shape type (hexahedron)
   const int N_shapes = 1;                     //  -- "" --
   
   void* coords[3];                            // Pointers to coordinate arrays
   coords[0] = vx_crds;
   coords[1] = vy_crds;
   coords[2] = vz_crds;
   
   // Write zone list into silo file:
   stringstream ss;
   ss << "VelGrid" << cellID;
   const string zoneListName = ss.str()+"Zones";
   if (DBPutZonelist2(fileptr,zoneListName.c_str(),N_zones,N_dims,nodeList,nodeListSize,0,0,0,shapeTypes,shapeSizes,shapeCnt,N_shapes,NULL) < 0) success = false;
   
   // Write grid into silo file:
   const string gridName = ss.str();
   if (DBPutUcdmesh(fileptr,gridName.c_str(),N_dims,NULL,coords,N_nodes,N_zones,zoneListName.c_str(),NULL,SiloType(bc_dataType,bc_dataSize),NULL) < 0) success = false;
   
   delete nodeList;
   delete vx_crds;
   delete vy_crds;
   delete vz_crds;
   delete bc_buffer;
   
   list<string> blockVarNames;
   if (vlsvReader.getBlockVariableNames(meshName,blockVarNames) == false) {
      cerr << "Failed to read block variable names!" << endl;
      success = false;
   }
   if (success == true) {
      for (list<string>::iterator it=blockVarNames.begin(); it!=blockVarNames.end(); ++it) {
	 if (convertVelocityBlockVariable(vlsvReader,meshName,gridName,cellID,N_blocks,blockOffset,*it) == false) success = false;
      }
   }
   
   return success;
}

int main(int argn,char* args[]) {
   int ntasks, rank;
   MPI_Init(&argn, &args);
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   if (rank == 0 && argn != 3) {
      cout << endl;
      cout << "USAGE: ./vlsvextract <file name mask> <cell ID>" << endl;
      cout << endl;
      cout << "Each VLSV file in the currect directory is compared against the mask," << endl;
      cout << "and if the file name matches the mask, the given velocity grid is " << endl;
      cout << "written to a SILO file." << endl;
      cout << endl;
      cout << "Cell ID is the ID of the spatial cell whose velocity grid is to be extracted." << endl;
      cout << endl;
      return 1;
   }
   //const string fname = args[1];
   const string mask = args[1];
   const uint64_t cellID = atoi(args[2]);
   
   const string directory = ".";
   const string suffix = ".vlsv";
   DIR* dir = opendir(directory.c_str());
   if (dir == NULL) {
      cerr << "ERROR in reading directory contents!" << endl;
      closedir(dir);
      return 1;
   }
   
   VLSVReader vlsvReader;
   int filesFound = 0, entryCounter = 0;
   vector<string> fileList;
   struct dirent* entry = readdir(dir);
   while (entry != NULL) {
      const string entryName = entry->d_name;
      if (entryName.find(mask) == string::npos || entryName.find(suffix) == string::npos) {
	 entry = readdir(dir);
	 continue;
      }
      fileList.push_back(entryName);
      filesFound++;
      entry = readdir(dir);
   }
   if (rank == 0 && filesFound == 0) cout << "\t no matches found" << endl;
   closedir(dir);
   
   for(size_t entryName = 0; entryName < fileList.size(); entryName++) {
      if(entryCounter++%ntasks == rank) {
	 // Open VLSV file and read mesh names:
	 vlsvReader.open(fileList[entryName]);
	 list<string> meshNames;
	 if (vlsvReader.getMeshNames(meshNames) == false) {
	    cout << "\t file '" << fileList[entryName] << "' not compatible" << endl;
	    vlsvReader.close();
	    continue;
	 }
	 
	 // Create a new file suffix for the output file:
	 stringstream ss1;
	 ss1 << ".silo";
	 string newSuffix;
	 ss1 >> newSuffix;

	 // Create a new file prefix for the output file:
	 stringstream ss2;
	 ss2 << "velgrid" << '.' << cellID;
	 string newPrefix;
	 ss2 >> newPrefix; 
	 
	 // Replace .vlsv with the new suffix:
	 string fileout = fileList[entryName];
	 size_t pos = fileout.rfind(".vlsv");
	 if (pos != string::npos) fileout.replace(pos,5,newSuffix);
	 
	 pos = fileout.find(".");
	 if (pos != string::npos) fileout.replace(0,pos,newPrefix);
	 
	 // Create a SILO file for writing:
	 fileptr = DBCreate(fileout.c_str(),DB_CLOBBER,DB_LOCAL,"Vlasov data file",DB_PDB);
	 if (fileptr == NULL) {
	    cerr << "\t failed to create output SILO file for input file '" << fileList[entryName] << "'" << endl;
	    DBClose(fileptr);
	    vlsvReader.close();
	    continue;
	 }
	 
	 // Extract velocity grid from VLSV file, if possible, and convert into SILO format:
	 bool velGridExtracted = true;
	 for (list<string>::const_iterator it=meshNames.begin(); it!=meshNames.end(); ++it) {
	    if (convertVelocityBlocks2(vlsvReader,*it,cellID) == false) {
	       velGridExtracted = false;
	    } else {
	       cout << "\t extracted from '" << fileList[entryName] << "'" << endl;
	    }
	 }
	 DBClose(fileptr);
	 
	 // If velocity grid was not extracted, delete the SILO file:
	 if (velGridExtracted == false) {
	    if (remove(fileout.c_str()) != 0) {
	       cerr << "\t ERROR: failed to remote dummy output file!" << endl;
	    }
	 }
	 
	 vlsvReader.close();
      }
   }
   
   MPI_Finalize();
   return 0;
}


