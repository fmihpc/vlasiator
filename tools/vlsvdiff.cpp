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
#include <stdint.h>
#include <cmath>
#include <limits> // YK
#include <list>
#include <set>
#include <sstream>
#include <dirent.h>

#include "vlsvreader2.h"
#include "definitions.h"

using namespace std;

bool convertMesh(VLSVReader& vlsvReader,
		 const string& meshName,
		 const char * varToExtract,
		 const uint compToExtract,
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
   if (compToExtract + 1 > variableVectorSize) {
      cerr << "ERROR invalid component, this variable has size " << variableVectorSize << endl;
      abort();
   }
   
   // Read the mesh array one node (of a spatial cell) at a time 
   // and create a map which contains each cell's CellID and variable to be extracted
   char* meshBuffer = new char[meshVectorSize*meshDataSize];
   uint* meshPtr = reinterpret_cast<uint*>(meshBuffer);
   char* variableBuffer = new char[variableVectorSize*variableDataSize];
   Real* variablePtrReal = reinterpret_cast<Real*>(variableBuffer);
   uint* variablePtrUint = reinterpret_cast<uint*>(variableBuffer);
   int* variablePtrInt = reinterpret_cast<int*>(variableBuffer);
   
   for (uint64_t i=0; i<meshArraySize; ++i) {
      if (vlsvReader.readArray("MESH", meshName, i, 1, meshBuffer) == false) {meshSuccess = false; break;}
      if (vlsvReader.readArray("VARIABLE", varToExtract, i, 1, variableBuffer) == false) {variableSuccess = false; break;}
      // Get the CellID
      uint CellID  = meshPtr[0];
      // Get the variable value
      Real extract = NAN;
      switch (variableDataType)
      {
	 case VLSV::FLOAT:
	    extract = variablePtrReal[compToExtract];
	    break;
	 case VLSV::UINT:
	    extract = (Real)(variablePtrUint[compToExtract]);
	    break;
	 case VLSV::INT:
	    extract = (Real)(variablePtrInt[compToExtract]);
	    break;
      }
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
		 const uint compToExtract,
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

bool shiftAverage(const map<uint, Real> * const orderedData1,
		  const map<uint, Real> * const orderedData2,
		  map<uint, Real> * shiftedData2
		 )
{
   map<uint, Real>::const_iterator it1, it2;
   Real avg1 = 0.0;
   Real avg2 = 0.0;
   
   for(it1=orderedData1->begin(), it2=orderedData2->begin();
       it1 != orderedData1->end(), it2 != orderedData2->end();
       it1++, it2++)
   {
      avg1 += orderedData1->at(it1->first);
      avg2 += orderedData2->at(it2->first);
   }
   avg1 /= orderedData1->size();
   avg2 /= orderedData1->size();
   
   for(it2=orderedData2->begin(); it2 != orderedData2->end(); it2++)
   {
      shiftedData2->insert(pair<uint, Real>(it2->first, it2->second - avg2 + avg1));
   }
   
   return 0;
}

bool pDistance(map<uint, Real> * const orderedData1,
	       map<uint, Real> * const orderedData2,
	       creal p,
	       Real * absolute,
	       Real * relative,
	       bool doShiftAverage
	      )
{
   // Computes the relative and absolute p-distance between two datasets X(x) provided in the maps
   /*   
    *   for p != 0.0
    *   absolute p-distance defined as:
    *   ||X_1 - X_2||_p = [\sum_i |X_1(i) - X_2(i)|^p]^{1/p}
    *   relative p-distance defined as:
    *   ||X_1 - X_2||_p = [\sum_i |X_1(i) - X_2(i)|^p]^{1/p} / ||X_1||_p
    * 
    *   for p == 0.0 it is the infinity distance
    *   absolute infinity-distance defined as:
    *    ||X_1 - X_2||_\infinity = max_i(|X_1(i) - X_2(i)|)
    *    relative infinity-distance defined as:
    *    ||X_1 - X_2||_\infinity = max_i(|X_1(i) - X_2(i)|) / ||X_1||_\infinity
    */
   map<uint, Real> shiftedData2;
   map<uint, Real> * data2 = orderedData2;
   
   if(doShiftAverage == true)
   {
      shiftAverage(orderedData1, orderedData2, &shiftedData2);
      data2 = &shiftedData2;
   }
   
   *absolute = 0.0;
   Real length = 0.0;
   *relative = 0.0;
   map<uint, Real>::const_iterator it;
   if(p != 0.0)
   {
      for(it=orderedData1->begin(); it != orderedData1->end() ; it++)
      {
	 *absolute += pow(abs(orderedData1->at(it->first) - data2->at(it->first)), p);
	 length += pow(abs(orderedData1->at(it->first)), p);
      }
      *absolute = pow(*absolute, 1.0 / p);
      length = pow(length, 1.0 / p);
   }
   else
   {
      for(it=orderedData1->begin(); it != orderedData1->end() ; it++)
      {
	 *absolute = max(*absolute, abs(orderedData1->at(it->first) - data2->at(it->first)));
	 length = max(length, abs(orderedData1->at(it->first)));
      }
   }
   
   if(length != 0.0)
   {
      *relative = *absolute / length;
   }
   else
   {
      cout << "WARNING (pDistance) : length of reference is 0.0, cannot divide to give relative distance." << endl;
      *relative = -1;
   }
   return 0;
}

bool singleStatistics(map<uint, Real> * orderedData,
                      Real * size,
		      Real * mini,
		      Real * maxi,
		      Real * avg,
		      Real * stdev
)
{
   /*
    * Returns basic statistics on the map passed to it.
    * 
    * 
    */
   map<uint, Real>::const_iterator it;
   
   *size = orderedData->size();
   *mini = numeric_limits<Real>::max();
   *maxi = numeric_limits<Real>::min();
   *avg = 0.0;
   *stdev = 0.0;
   
   for(it=orderedData->begin(); it != orderedData->end() ; it++)
   {
      *mini = min(*mini, orderedData->at(it->first));
      *maxi = max(*maxi, orderedData->at(it->first));
      *avg += orderedData->at(it->first);
   }
   *avg /= *size;
   for(it=orderedData->begin(); it != orderedData->end() ; it++)
   {
      *stdev += pow(orderedData->at(it->first) - *avg, 2.0);
   }
   *stdev = sqrt(*stdev);
   *stdev /= (*size - 1);
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
   uint compToExtract = atoi(args[4]);
   
   map<uint, Real> orderedData1;
   map<uint, Real> orderedData2;
   Real absolute, relative, mini, maxi, size, avg, stdev;
   
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
   
   singleStatistics(&orderedData1, &size, &mini, &maxi, &avg, &stdev);
   cout << "Statistics on first file: size " << size
        << " min = " << mini
        << " max = " << maxi
        << " average = " << avg
        << " standard deviation " << stdev
        << endl;
   
   singleStatistics(&orderedData2, &size, &mini, &maxi, &avg, &stdev);
   cout << "Statistics on second file: size " << size
	<< " min = " << mini
	<< " max = " << maxi
	<< " average = " << avg
	<< " standard deviation " << stdev
	<< endl;
   
   pDistance(&orderedData1, &orderedData2, 1, &absolute, &relative, false);
   cout << "The absolute 1-distance between both datasets is " << absolute << "." << endl;
   cout << "The relative 1-distance between both datasets is " << relative << "." << endl;
   
   pDistance(&orderedData1, &orderedData2, 1, &absolute, &relative, true);
   cout << "The average-shifted absolute 1-distance between both datasets is " << absolute << "." << endl;
   cout << "The average-shifted relative 1-distance between both datasets is " << relative << "." << endl;
   
   pDistance(&orderedData1, &orderedData2, 2, &absolute, &relative, false);
   cout << "The absolute 2-distance between both datasets is " << absolute << "." << endl;
   cout << "The relative 2-distance between both datasets is " << relative << "." << endl;
   pDistance(&orderedData1, &orderedData2, 2, &absolute, &relative, true);
   cout << "The average-shifted absolute 2-distance between both datasets is " << absolute << "." << endl;
   cout << "The average-shifted relative 2-distance between both datasets is " << relative << "." << endl;
   
   pDistance(&orderedData1, &orderedData2, 0, &absolute, &relative, false);
   cout << "The absolute infinity-distance between both datasets is " << absolute << "." << endl;
   cout << "The relative infinity-distance between both datasets is " << relative << "." << endl;
   
   pDistance(&orderedData1, &orderedData2, 0, &absolute, &relative, true);
   cout << "The average-shifted absolute infinity-distance between both datasets is " << absolute << "." << endl;
   cout << "The average-shifted relative infinity-distance between both datasets is " << relative << "." << endl;
   
   return 0;
}
