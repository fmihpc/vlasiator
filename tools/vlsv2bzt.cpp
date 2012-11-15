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

/*! \file vlsv2bzt.cpp
 *  \brief Extraction utility to get an (x-t) dataset for further analysis (e.g. Fourier transform).
 * 
 * This tool is used to extract the spatial profile of a given variable from a 1D dataset at each time step and save this into an ascii or binary file to be imported into MATLAB/Octave/Python... for further analysis, e.g. 2D Fourier transformation to produce dispersion plots.
 * 
 * The calling pattern is:
 * "$ ./vlsv2bzt <Variable name> <component> <Output file name> <Output file format> <input file mask(s)>"
 */

#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <cmath>
#include <list>
#include <set>
#include <sstream>
#include <dirent.h>
#include <fcntl.h>

#include <mpi.h>

#include "vlsvreader2.h"
#include "definitions.h"

using namespace std;

/*! Extracts the dataset from the VLSV file opened by convertSILO and sends it to the master rank for writing.
 * \param vlsvReader VLSVReader class object used to access the VLSV file
 * \param entryName Identification number of the file, used in MPI communication
 * \param meshName Address of the string containing the name of the mesh to be extracted
 * \param varToExtract Pointer to the char array containing the name of the variable to extract
 * \param compToExtract Unsigned int designating the component to extract (0 for scalars)
 */
bool convertMesh(VLSVReader& vlsvReader,
   const string& meshName,
   const unsigned long int entryName,
   const char * varToExtract,
   const int compToExtract
) {
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
   if (vlsvReader.getArrayInfo("VARIABLE", varToExtract, meshName, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false) return false;
   if (coordArraySize != variableArraySize) {
      cerr << "ERROR array size mismatch" << endl;
   }
   
   // Read the coordinate array one node (of a spatial cell) at a time 
   // and create a map which contains each cell's x and variable to be extracted
   char* coordBuffer = new char[coordVectorSize*coordDataSize];
   Real* coordPtr = reinterpret_cast<Real*>(coordBuffer);
   char* variableBuffer = new char[variableVectorSize*variableDataSize];
   Real* variablePtr = reinterpret_cast<Real*>(variableBuffer);
   for (uint64_t i=0; i<coordArraySize; ++i) {
      if (vlsvReader.readArray("COORDS",meshName,i,1,coordBuffer) == false) {coordSuccess = false; break;}
      if (vlsvReader.readArray("VARIABLE", varToExtract, i, 1, variableBuffer) == false) {variableSuccess = false; break;}
      
      // Get the x coordinate
      creal x  = coordPtr[0]; /*!< Changing this allows to extract datasets along x (0), y (1) or z (2). */
      // Get the variable value
      creal extract = variablePtr[compToExtract];
      // Put those into the map
      orderedData.insert(pair<Real, Real>(x, extract));
   }
   if (coordSuccess == false) {
      cerr << "ERROR reading array COORDS" << endl;
   }
   if (variableSuccess == false) {
      cerr << "ERROR reading array VARIABLE " << varToExtract << endl;
   }
   unsigned long int numberOfValuesToSend = orderedData.size();
   vector<Real> valuesToSend(numberOfValuesToSend);
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
   MPI_Ssend(&valuesToSend[0], numberOfValuesToSend, MPI_DOUBLE, 0, entryName + 1, MPI_COMM_WORLD);
   return coordSuccess && variableSuccess;
}

/*! Opens the VLSV file and extracts the mesh names. Sends for processing to convertMesh.
 * \param fileName String containing the name of the file to be processed
 * \param entryName Identification number of the file, used in MPI communication
 * \param varToExtract Pointer to the char array containing the name of the variable to extract
 * \param compToExtract Unsigned int designating the component to extract (0 for scalars)
 * \sa convertMesh
 */
bool convertSILO(const string& fileName,
   const unsigned long int entryName,
   const char * varToExtract,
   const int compToExtract
) {
   bool success = true;
   
   // Open VLSV file for reading:
   VLSVReader vlsvReader;
   if (vlsvReader.open(fileName) == false) {
      cerr << "Failed to open '" << fileName << "'" << endl;
      return false;
   }
   
   // Get the names of all meshes in vlsv file
   list<string> meshNames;
   if (vlsvReader.getMeshNames(meshNames) == false) {
      return false;
   }
   for (list<string>::const_iterator it=meshNames.begin(); it!=meshNames.end(); ++it) {
      if (convertMesh(vlsvReader, *it, entryName, varToExtract, compToExtract) == false) {
         return false;
      }
   }
   
   vlsvReader.close();
   return success;
}

/*! Main function, reads the arguments, processes the file masks to get the file list and distributes the work.
 * 
 * The parallelisation follows the master-slave pattern, process 0 distributes tasks and receives the result, which it writes to file.
 * 
 * Using binary is ~3 times faster and produces ~0.6 times the volume of data in ascii.
 */
int main(int argn,char* args[]) {
   int ntasks, rank;
   MPI_Init(&argn, &args);
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   if ((argn < 6) && (rank == 0)) {
      cout << endl;
      cout << "USAGE: ./vlsv2bzt <Variable name> <component> <Output file name> <Output file format> <input file mask(s)>" << endl;
      cout << "Each VLSV in the current directory is compared against the given file mask(s)," << endl;
      cout << "and if match is found, the variable's component is extracted and ordered before being output to a file *with the name of the 3rd argument* and the *format given by the 4th argument* (a: ascii, b: binary)." << endl;
      cout << "Note that due to the parallelisation algorithm, there should be at least 2 processes used." << endl;
      cout << endl;
      return 1;
   }
   // 1st arg is variable name
   char * varToExtract = args[1];
   // 2nd arg is its component, 0 for scalars, 2 for z component etc
   int compToExtract = atoi(args[2]);
   // 3rd arg is file name
   char * fileName = args[3];
   // 4th arg is format
   char * fileFormat = args[4];
   
   // Convert file masks into strings:
   vector<string> masks;
   set<string> fileList;
   for (int i=5; i<argn; ++i) masks.push_back(args[i]);

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
            convertSILO(*it, entryName, varToExtract, compToExtract);
            filesConverted++;
         }
      }
      if (filesConverted == 0) cout << "\t no files converted" << endl;
   } else {
      // YK Create/open file to append to
      FILE* outputFile = NULL;
      char endOfLine = '\n';
      if(*fileFormat == 'a') outputFile = fopen(fileName, "w");
      if(*fileFormat == 'b') outputFile = fopen(fileName, "wb");
      
      int lineLength;
      MPI_Recv(&lineLength, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      vector<Real> recvBuffer(lineLength);
      for(unsigned long int entryName = 0; entryName < fileList.size(); entryName++) {
         MPI_Recv(&recvBuffer[0], lineLength, MPI_DOUBLE, entryName%(ntasks - 1) + 1, entryName + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         if(*fileFormat == 'a') {
            for(int i = 0; i < lineLength; i++) {
               fprintf(outputFile, "%e ", recvBuffer[i]);
            }
            fprintf(outputFile, "\n");
         }
         if(*fileFormat == 'b') {
            fwrite(&recvBuffer[0], sizeof(recvBuffer[0]), lineLength, outputFile);
//            fwrite(&endOfLine, sizeof(endOfLine), 1, outputFile);
         }
      }
      if(*fileFormat == 'b') {
         creal format[2] = {fileList.size()*1.0, lineLength*1.0};
         fwrite(&format[0], sizeof(format[0]), 2, outputFile);
      }
      fclose(outputFile);
   }
   MPI_Finalize();
   return 0;
}
