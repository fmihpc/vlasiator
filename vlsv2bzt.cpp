/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

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
#include <stdint.h>
#include <cmath>
#include <list>
#include <set>
#include <sstream>
#include <dirent.h>

#include <mpi.h>

#include "vlsvreader2.h"
#include "definitions.h"

using namespace std;

bool convertMesh(VLSVReader& vlsvReader,
		 const string& meshName,
		 const unsigned long int entryName) {
   bool coordSuccess = true;
   bool variableSuccess = true;
   
   // Creating a map from the x coordinate to the variable sought
   // This gives the spatial ordering needed (default is a<b)
   map<Real, Real> orderedData;
   
   VLSV::datatype coordDataType;
   VLSV::datatype variableDataType;
   uint64_t coordArraySize, coordVectorSize, coordDataSize;
   uint64_t variableArraySize, variableVectorSize, variableDataSize;
   
   if (vlsvReader.getArrayInfo("COORDS", meshName, coordArraySize, coordVectorSize, coordDataType, coordDataSize) == false) return false;
   if (vlsvReader.getArrayInfo("VARIABLE", "B", meshName, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false) return false;
   if (coordArraySize != variableArraySize) {
      cerr << "ERROR array size mismatch" << endl;
   }
   
   // Read the coordinate array one node (of a spatial cell) at a time 
   // and create a map which contains each cell's x and Bz
   char* coordBuffer = new char[coordVectorSize*coordDataSize];
   Real* coordPtr = reinterpret_cast<Real*>(coordBuffer);
   char* variableBuffer = new char[variableVectorSize*variableDataSize];
   Real* variablePtr = reinterpret_cast<Real*>(variableBuffer);
   for (uint64_t i=0; i<coordArraySize; ++i) {
      if (vlsvReader.readArray("COORDS",meshName,i,1,coordBuffer) == false) {coordSuccess = false; break;}
      if (vlsvReader.readArray("VARIABLE", "B", i, 1, variableBuffer) == false) {variableSuccess = false; break;}
      
      // Get the x coordinate
      creal x  = coordPtr[0];
      // Get the Bz value
      creal Bz = variablePtr[2];
      // Put those into the map
      orderedData.insert(pair<Real, Real>(x, Bz));
   }
   if (coordSuccess == false) {
      cerr << "ERROR reading array COORDS" << endl;
   }
   if (variableSuccess == false) {
      cerr << "ERROR reading array VARIABLE (B)" << endl;
   }
   unsigned long int numberOfValuesToSend = orderedData.size();
   Real valuesToSend[numberOfValuesToSend];
   unsigned long int j;
   map<Real, Real>::const_iterator it;
   for(it = orderedData.begin(), j = 0;
       it != orderedData.end();
       it++, j++) {
      valuesToSend[j] = orderedData[it->first];
      }
   static int rank = -1;
   if(rank == -1) {MPI_Comm_rank(MPI_COMM_WORLD, &rank);}
   if(rank == 1) {
      static bool sendSize = true;
      if(sendSize)
      {
	 sendSize == false;
	 MPI_Send(&numberOfValuesToSend, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      }
   }
   MPI_Ssend(&valuesToSend, numberOfValuesToSend, MPI_DOUBLE, 0, entryName + 1, MPI_COMM_WORLD);
   return coordSuccess && variableSuccess;
}

bool convertSILO(const string& fname, const unsigned long int entryName) {
   bool success = true;
   
   // Open VLSV file for reading:
   VLSVReader vlsvReader;
   if (vlsvReader.open(fname) == false) {
      cerr << "Failed to open '" << fname << "'" << endl;
      return false;
   }
   
   // Get the names of all meshes in vlsv file, and write into silo file:
   list<string> meshNames;
   if (vlsvReader.getMeshNames(meshNames) == false) {
      return false;
   }
   for (list<string>::const_iterator it=meshNames.begin(); it!=meshNames.end(); ++it) {
      if (convertMesh(vlsvReader, *it, entryName) == false) {
	 return false;
      }
   }
   
   vlsvReader.close();
   return success;
}

int main(int argn,char* args[]) {
   int ntasks, rank;
   MPI_Init(&argn, &args);
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   if (argn < 3) {
      cout << endl;
      cout << "USAGE: ./vlsv2bzt <Output file name> <input file mask(s)>" << endl;
      cout << "Each VLSV in the current directory is compared against the given file mask(s)," << endl;
      cout << "and if match is found, Bz is extracted and ordered before being output to an ascii file *with the name of the 1st argument*." << endl;
      cout << endl;
      return 1;
   }
   
   // Convert file masks into strings:
   vector<string> masks;
   set<string> fileList;
   for (int i=2; i<argn; ++i) masks.push_back(args[i]);

   // Compare directory contents against each mask:
   const string directory = ".";
   const string suffix = ".vlsv";
   unsigned long int filesFound = 0, filesConverted = 0;
   for (size_t mask=0; mask<masks.size(); ++mask) {
      if(rank == 0) {cout << "Comparing mask '" << masks[mask] << "'" << endl;}
      DIR* dir = opendir(directory.c_str());
      if (dir == NULL) continue;
      
      struct dirent* entry = readdir(dir);
      while (entry != NULL) {
	 const string entryName = entry->d_name;
	 // Compare entry name against given mask and file suffix ".vlsv":
	 if (entryName.find(masks[mask]) == string::npos || entryName.find(suffix) == string::npos) {
	    entry = readdir(dir);
	    continue;
	 }
	 fileList.insert(entryName);
	 filesFound++;
	 entry = readdir(dir);
      }
      closedir(dir);
      if (rank == 0 && filesFound == 0) cout << "\t no matches found" << endl;
   }
   
   if(rank != 0) {
      unsigned long int entryName;
      set<string>::iterator it;
      for(entryName = 0, it=fileList.begin(); entryName < fileList.size(); entryName++, it++)
      {
	 if(entryName%(ntasks - 1) + 1 == rank) {
	    cout << "\tProc " << rank << " converting '" << *it << "'" << endl;
	    convertSILO(*it, entryName);
	    filesConverted++;
	 }
      }
      if (filesConverted == 0) cout << "\t no files converted" << endl;
   } else {
      // YK Create/open ascii file to append to
      FILE * outputFile;
      outputFile = fopen(args[1], "w");
      
      int lineLength;
      MPI_Recv(&lineLength, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      Real recvBuffer[lineLength];
      for(unsigned long int entryName = 0; entryName < fileList.size(); entryName++)
      {
	 MPI_Recv(&recvBuffer, lineLength, MPI_DOUBLE, entryName%(ntasks - 1) + 1, entryName + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 for(int i = 0; i < lineLength; i++)
	 {
	    fprintf(outputFile, "%e ", recvBuffer[i]);
	 }
	 fprintf(outputFile, "\n");
      }
      fclose(outputFile);
   }
   MPI_Finalize();
   return 0;
}
