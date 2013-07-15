/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












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

#include "vlsvreader2.h"
#include "definitions.h"

using namespace std;





static DBfile* fileptr = NULL; // Pointer to file opened by SILO

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



   
bool convertMeshVariable(VLSVReader& vlsvReader,const string& meshName,const string& varName) {
   bool success = true;

   // Writing a unstructured grid variable is a rather straightforward process. The 
   // only compilation here is that some of the variables are actually vectors, i.e. 
   // vectorSize > 1 (vectorSize == 1 for scalars). Format in which vectors are stored in VLSV 
   // differ from format in which they are written to SILO files.
   VLSV::datatype dataType;
   uint64_t arraySize,vectorSize,dataSize;
   if (vlsvReader.getArrayInfo("VARIABLE",varName,meshName,arraySize,vectorSize,dataType,dataSize) == false) return false;
   
   // Read variable data. Note that we do not actually need to care if 
   // the data is given as floats or doubles.
   char* buffer = new char[arraySize*vectorSize*dataSize];
   if (vlsvReader.readArray("VARIABLE",varName,0,arraySize,buffer) == false) success = false;
   if (success == false) {
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
   if (DBPutUcdvar(fileptr,varName.c_str(),meshName.c_str(),vectorSize,&(varNamePtrs[0]),components,arraySize,NULL,0,
		   SiloType(dataType,dataSize),DB_ZONECENT,NULL) < 0) success = false;

   for (uint64_t i=0; i<vectorSize; ++i) {
      delete [] components[i];
   }
   delete [] components;
   delete [] buffer;   
   return success;
}

bool convertMesh(VLSVReader& vlsvReader,const string& meshName) {
   bool success = true;
   const float EPS = 1.0e-7;
   
   // First task is to push all unique node coordinates into a map.
   // This is not too difficult for unrefined grids, since each spatial cell stores
   // its bottom lower left corner coordinate and size. For refined grid the situation
   // is more complex, as there are more unique nodes than the lower left corners:
   map<NodeCrd<Real>,uint64_t,NodeComp> nodes;
   
   VLSV::datatype dataType;
   uint64_t arraySize,vectorSize,dataSize;
   if (vlsvReader.getArrayInfo("COORDS",meshName,arraySize,vectorSize,dataType,dataSize) == false) return false;
   
   // Read the coordinate array one node (of a spatial cell) at a time 
   // and create a map which only contains each existing node once.
   char* coordsBuffer = new char[vectorSize*dataSize];
   Real* ptr = reinterpret_cast<Real*>(coordsBuffer);
   for (uint64_t i=0; i<arraySize; ++i) {
      if (vlsvReader.readArray("COORDS",meshName,i,1,coordsBuffer) == false) {success = false; break;}
      
      // Insert all eight nodes of a cell into map nodes.
      // NOTE: map is a unique associative container - given a suitable comparator, map 
      // will filter out duplicate nodes.
      creal x  = ptr[0];
      creal y  = ptr[1];
      creal z  = ptr[2];
      creal dx = ptr[3];
      creal dy = ptr[4];
      creal dz = ptr[5];

      Real X0 = x;
      Real X1 = x+dx;
      Real Y0 = y;
      Real Y1 = y+dy;
      Real Z0 = z;
      Real Z1 = z+dz;
      
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
   int* nodelist = new int[8*arraySize];
   for (uint64_t i=0; i<arraySize; ++i) {
      // Read the bottom lower left corner coordinates of a cell and its sizes. Note 
      // that zones will end up in SILO file in the same order as they are in VLSV file.
      if (vlsvReader.readArray("COORDS",meshName,i,1,coordsBuffer) == false) {
	 cerr << "Failed to read array coords" << endl;
	 success = false; 
	 break;
      }
      creal x  = ptr[0];
      creal y  = ptr[1];
      creal z  = ptr[2];
      creal dx = ptr[3];
      creal dy = ptr[4];
      creal dz = ptr[5];

      // Calculate x,y,z coordinates of the eight nodes of the cell:
      Real X0 = x;
      Real X1 = x+dx;
      Real Y0 = y;
      Real Y1 = y+dy;
      Real Z0 = z;
      Real Z1 = z+dz;
      
      // Flush very small coordinate values to zero:
      if (fabs(X0) < EPS) X0 = 0.0;
      if (fabs(X1) < EPS) X1 = 0.0;
      if (fabs(Y0) < EPS) Y0 = 0.0;
      if (fabs(Y1) < EPS) Y1 = 0.0;
      if (fabs(Z0) < EPS) Z0 = 0.0;
      if (fabs(Z1) < EPS) Z1 = 0.0;
      
      // Search the cell's nodes from the map created above. For each node in nodelist 
      // store an index into an array which only contains the unique nodes:
      map<NodeCrd<Real>,uint64_t,NodeComp>::const_iterator it;
      it = nodes.find(NodeCrd<Real>(X0,Y0,Z0)); if (it == nodes.end()) success = false; nodelist[i*8+0] = it->second;
      it = nodes.find(NodeCrd<Real>(X1,Y0,Z0)); if (it == nodes.end()) success = false; nodelist[i*8+1] = it->second;
      it = nodes.find(NodeCrd<Real>(X1,Y1,Z0)); if (it == nodes.end()) success = false; nodelist[i*8+2] = it->second;
      it = nodes.find(NodeCrd<Real>(X0,Y1,Z0)); if (it == nodes.end()) success = false; nodelist[i*8+3] = it->second;
      it = nodes.find(NodeCrd<Real>(X0,Y0,Z1)); if (it == nodes.end()) success = false; nodelist[i*8+4] = it->second;
      it = nodes.find(NodeCrd<Real>(X1,Y0,Z1)); if (it == nodes.end()) success = false; nodelist[i*8+5] = it->second;
      it = nodes.find(NodeCrd<Real>(X1,Y1,Z1)); if (it == nodes.end()) success = false; nodelist[i*8+6] = it->second;
      it = nodes.find(NodeCrd<Real>(X0,Y1,Z1)); if (it == nodes.end()) success = false; nodelist[i*8+7] = it->second;
   }
   delete coordsBuffer;
   if (success == false) {
      cerr << "Failed to find node(s)" << endl;
   }
   
   // Write the unstructured mesh to SILO file:
   const int N_dims  = 3;                      // Number of dimensions
   const int N_nodes = nodes.size();           // Total number of nodes
   const int N_zones = arraySize;              // Total number of zones (=spatial cells)
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
   if (DBPutZonelist2(fileptr,zoneListName.c_str(),N_zones,N_dims,nodelist,8*arraySize,0,0,0,shapeTypes,shapeSizes,shapeCnt,N_shapes,NULL) < 0) success = false;
   
   // Write grid into silo file:
   if (DBPutUcdmesh(fileptr,meshName.c_str(),N_dims,NULL,coords,N_nodes,N_zones,zoneListName.c_str(),NULL,SiloType(dataType,dataSize),NULL) < 0) success = false;

   nodes.clear();
   delete nodelist;
   delete xcrds;
   delete ycrds;
   delete zcrds;

   // Write the cell IDs as a variable:
   if (vlsvReader.getArrayInfo("MESH",meshName,arraySize,vectorSize,dataType,dataSize) == false) return false;
   char* buffer = new char[arraySize*vectorSize*dataSize];
   if (vlsvReader.readArray("MESH",meshName,0,arraySize,buffer) == false) success = false;
   string cellIDlabel = "Cell ID";
   DBoptlist* optList = DBMakeOptlist(1);   
   DBAddOption(optList,DBOPT_LABEL,const_cast<char*>(cellIDlabel.c_str()));   
   if (DBPutUcdvar1(fileptr,"CellID",meshName.c_str(),buffer,arraySize,NULL,0,SiloType(dataType,dataSize),DB_ZONECENT,NULL) < 0) success = false;
   delete buffer;
   DBFreeOptlist(optList);
   
   // Write all variables of this mesh into silo file:
   list<string> variables;
   vlsvReader.getVariableNames(meshName,variables);
   for (list<string>::const_iterator it=variables.begin(); it!=variables.end(); ++it) {
      if (convertMeshVariable(vlsvReader,meshName,*it) == false) success = false;
   }
   return success;
  
}

bool convertSILO(const string& fname) {
   bool success = true;

   
   // Open VLSV file for reading:
   VLSVReader vlsvReader;
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
         convertSILO(fileList[entryName]);
	 filesConverted++;
      }
   }
   
   int totalFilesConverted =0;   
   MPI_Reduce(&filesConverted, &totalFilesConverted, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   if (rank == 0 && totalFilesConverted == 0) cout << "\t no files converted" << endl;
   
   MPI_Finalize();
   return 0;
}
