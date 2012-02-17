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

# include <cmath>

#include "vlsvreader2.h"
#include "definitions.h"

using namespace std;

bool convertMesh(VLSVReader& vlsvReader,
		 const string& meshName,
		 const char * varToExtract,
		 const int compToExtract,
		 map<uint, Real> * orderedData)
{
   bool meshSuccess = true;
   bool variableSuccess = true;
   
   VLSV::datatype meshDataType;
   VLSV::datatype variableDataType;
   uint64_t meshArraySize, meshVectorSize, meshDataSize;
   uint64_t variableArraySize, variableVectorSize, variableDataSize;
   
   if (vlsvReader.getArrayInfo("MESH", meshName, meshArraySize, meshVectorSize, meshDataType, meshDataSize) == false) return false;
   if (vlsvReader.getArrayInfo("VARIABLE", varToExtract, meshName, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false) return false;
   if (meshArraySize != variableArraySize) {
      cerr << "ERROR array size mismatch" << endl;
   }
   
   // Read the mesh array one node (of a spatial cell) at a time 
   // and create a map which contains each cell's CellID and variable to be extracted
   char* meshBuffer = new char[meshVectorSize*meshDataSize];
   uint* meshPtr = reinterpret_cast<uint*>(meshBuffer);
   char* variableBuffer = new char[variableVectorSize*variableDataSize];
   Real* variablePtr = reinterpret_cast<Real*>(variableBuffer);
   for (uint64_t i=0; i<meshArraySize; ++i) {
      if (vlsvReader.readArray("MESH", meshName, i, 1, meshBuffer) == false) {meshSuccess = false; break;}
      if (vlsvReader.readArray("VARIABLE", varToExtract, i, 1, variableBuffer) == false) {variableSuccess = false; break;}
      // Get the CellID
      uint CellID  = meshPtr[0];
      // Get the variable value
      Real extract = variablePtr[compToExtract];
      // Put those into the map
      orderedData->insert(pair<uint, Real>(CellID, extract));
   }
   
   if (meshSuccess == false) {
      cerr << "ERROR reading array MESH" << endl;
   }
   if (variableSuccess == false) {
      cerr << "ERROR reading array VARIABLE " << varToExtract << endl;
   }
   return meshSuccess && variableSuccess;
}

bool convertSILO(const string fileName,
		 const char * varToExtract,
		 const int compToExtract,
		 map<uint, Real> * orderedData)
{
   bool success = true;
   
   // Open VLSV file for reading:
   VLSVReader vlsvReader;
   
   if (vlsvReader.open(fileName) == false)
   {
      cerr << "Failed to open '" << fileName << "'" << endl;
      return false;
   }
   
   // Get the names of all meshes in vlsv file
   list<string> meshNames;
   if (vlsvReader.getMeshNames(meshNames) == false)
   {
      return false;
   }
   
   for (list<string>::const_iterator it=meshNames.begin(); it!=meshNames.end(); ++it)
   {
      if (convertMesh(vlsvReader, *it, varToExtract, compToExtract, orderedData) == false)
      {
	 return false;
      }
      
   }
   vlsvReader.close();
   return success;
}

bool pDistance(map<uint, Real> * orderedData1, map<uint, Real> * orderedData2, creal p)
{
   // Computes the p-distance between two datasets X(x) provided in the maps
   /*
     p-distance defined as:
     ||X_1 - X_2||_p = [\sum_i |X_1(i) - X_2(i)|^p]^{1/p}
   */
   
   Real result;
   map<uint, Real>::const_iterator it;
   
   for(it=orderedData1->begin(); it != orderedData1->end() ; it++)
   {
      result += pow(abs(orderedData1->at(it->first) - orderedData2->at(it->first)), p);
   }
   result = pow(result, 1.0 / p);
   
   cout << "The " << p << "-distance between both datasets is " << result << "." << endl;
   return 0;
}

int main(int argn,char* args[])
{
   if (argn < 5)
   {
      cout << endl;
      cout << "USAGE: ./vlsvdiff <file1> <file2> <Variable> <component>" << endl;
      cout << "Performs various comparison operations between the two files given, for the variable and component given" << endl;
      cout << endl;
      return 1;
   }
   
   // 1st arg is file1 name
   const string fileName1 = args[1];
   // 2nd arg is file2 name
   const string fileName2 = args[2];
   // 3rd arg is variable name
   char * varToExtract = args[3];
   // 4th arg is its component, 0 for scalars, 2 for z component etc
   int compToExtract = atoi(args[4]);
   
   map<uint, Real> orderedData1;
   map<uint, Real> orderedData2;
   
   if(convertSILO(fileName1, varToExtract, compToExtract, &orderedData1) != true)
   {
      cerr << "Data import error with " << fileName1 << endl;
      return 1;
   }
   
   if(convertSILO(fileName2, varToExtract, compToExtract, &orderedData2) != true)
   {
      cerr << "Data import error with " << fileName2 << endl;
      return 1;
   }
   
   // Basic consistency check
   if(orderedData1.size() != orderedData2.size())
   {
      cerr << "Datasets have differing size" << endl;
      return 1;
   }
   
   pDistance(&orderedData1, &orderedData2, 1);
   pDistance(&orderedData1, &orderedData2, 2);
   
   
   
   return 0;
}
