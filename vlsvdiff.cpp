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

bool infinityDistance(map<uint, Real> * orderedData1,
		      map<uint, Real> * orderedData2,
		      Real * absolute,
		      Real * relative
		     )
{
   /*
    * Computes the relative and absolute infinity-distance between two datasets X(x) provided in the maps
    *    absolute infinity-distance defined as:
    *    ||X_1 - X_2||_\infinity = max_i(|X_1(i) - X_2(i)|)
    * 
    *    relative infinity-distance defined as:
    *    ||X_1 - X_2||_\infinity = max_i(|(X_1(i) - X_2(i)) / X_1(i)|)
    */
   
   *absolute = 0.0;
   *relative = 0.0;
   map<uint, Real>::const_iterator it;
   
   for(it=orderedData1->begin(); it != orderedData1->end() ; it++)
   {
      *absolute = max(*absolute, abs(orderedData1->at(it->first) - orderedData2->at(it->first)));
      if(orderedData1->at(it->first) != 0.0)
      {
	 *relative = *absolute / abs(orderedData1->at(it->first));
      }
      else if (orderedData2->at(it->first) != 0.0)
      {
	 cout << "WARNING (infinityDistance) : 1st cell has 0.0, dividing by 2nd." << endl;
	 *relative = *absolute / abs(orderedData2->at(it->first));
      }
      else
      {
	 cout << "WARNING (infinityDistance) : both cells have 0.0, cell skipped." << endl;
      }
   }
   
   //cout << "The absolute infinity-distance between both datasets is " << absolute << "." << endl;
   //cout << "The relative infinity-distance between both datasets is " << relative << "." << endl;
   return 0;
}

bool infinityDistanceAvgSubtracted(map<uint, Real> * orderedData1,
				   map<uint, Real> * orderedData2,
				   Real * absolute,
				   Real * relative
)
{
   /*
    * Subtracts the average of the dataset from itself in order to compute
    * the relative and absolute infinity-distance between the two datasets X(x) provided in the maps
    * using infinityDistance().
    */
   
   map<uint, Real>::const_iterator it1, it2;
   map <uint, Real> subtractedData1, subtractedData2;
   Real avg1 = 0.0;
   Real avg2 = 0.0;
   
   for(it1=orderedData1->begin(); it1 != orderedData1->end() ; it1++)
   {
      avg1 += orderedData1->at(it1->first);
      avg2 += orderedData1->at(it1->first);
   }
   avg1 /= orderedData1->size();
   avg2 /= orderedData1->size();
   
   for(it1=orderedData1->begin(), it2=orderedData2->begin();
       it1 != orderedData1->end(), it2 != orderedData2->end();
       it1++, it2++)
   {
      subtractedData1.insert(pair<uint, Real>(it1->first, it1->second - avg1));
      subtractedData2.insert(pair<uint, Real>(it2->first, it2->second - avg2));
   }
   
   infinityDistance(&subtractedData1, &subtractedData2, absolute, relative);
   
   //cout << "The average-subtracted absolute infinity-distance between both datasets is " << absolute << "." << endl;
   //cout << "The average-subtracted relative infinity-distance between both datasets is " << relative << "." << endl;
   return 0;
}

bool pDistance(map<uint, Real> * orderedData1,
	       map<uint, Real> * orderedData2,
	       creal p,
	       Real * absolute,
	       Real * relative
	      )
{
   // Computes the relative and absolute p-distance between two datasets X(x) provided in the maps
   /*
    *   absolute p-distance defined as:
    *   ||X_1 - X_2||_p = [\sum_i |X_1(i) - X_2(i)|^p]^{1/p}
    *   relative p-distance defined as:
    *   ||X_1 - X_2||_p = [\sum_i |(X_1(i) - X_2(i)) / X_1(i)|^p]^{1/p}
    */
   
   *absolute = 0.0;
   *relative = 0.0;
   map<uint, Real>::const_iterator it;
   
   for(it=orderedData1->begin(); it != orderedData1->end() ; it++)
   {
      *absolute += pow(abs(orderedData1->at(it->first) - orderedData2->at(it->first)), p);
      if(orderedData1->at(it->first) != 0.0)
      {
	 *relative = *absolute / pow(abs(orderedData1->at(it->first)), p);
      }
      else if (orderedData2->at(it->first) != 0.0)
      {
	 cout << "WARNING (pDistance) : 1st cell has 0.0, dividing by 2nd." << endl;
	 *relative = *absolute / pow(abs(orderedData2->at(it->first)), p);
      }
      else
      {
	 cout << "WARNING (infinityDistance) : both cells have 0.0, cell skipped." << endl;
      }
   }
   *absolute = pow(*absolute, 1.0 / p);
   *relative = pow(*relative, 1.0 / p);
   //cout << "The absolute " << p << "-distance between both datasets is " << absolute << "." << endl;
   //cout << "The relative " << p << "-distance between both datasets is " << relative << "." << endl;
   return 0;
}

bool pDistanceAvgSubtracted(map<uint, Real> * orderedData1,
			    map<uint, Real> * orderedData2,
			    creal p,
			    Real * absolute,
			    Real * relative
)
{
   /*
    * Subtracts the average of the dataset from itself in order to compute
    * the relative and absolute p-distance between the two datasets X(x) provided in the maps
    * using pDistance().
    */
   
   map<uint, Real>::const_iterator it1, it2;
   map <uint, Real> subtractedData1, subtractedData2;
   Real avg1 = 0.0;
   Real avg2 = 0.0;
   
   for(it1=orderedData1->begin(); it1 != orderedData1->end() ; it1++)
   {
      avg1 += orderedData1->at(it1->first);
      avg2 += orderedData1->at(it1->first);
   }
   avg1 /= orderedData1->size();
   avg2 /= orderedData1->size();
   
   for(it1=orderedData1->begin(), it2=orderedData2->begin();
       it1 != orderedData1->end(), it2 != orderedData2->end();
   it1++, it2++)
       {
	  subtractedData1.insert(pair<uint, Real>(it1->first, it1->second - avg1));
	  subtractedData2.insert(pair<uint, Real>(it2->first, it2->second - avg2));
       }
       
       pDistance(&subtractedData1, &subtractedData2, p, absolute, relative);
       
       //cout << "The average-subtracted absolute p-distance between both datasets is " << absolute << "." << endl;
       //cout << "The average-subtracted relative p-distance between both datasets is " << relative << "." << endl;
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
   Real absolute, relative;
   
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
   
   pDistance(&orderedData1, &orderedData2, 1, &absolute, &relative);
   cout << "The absolute 1-distance between both datasets is " << absolute << "." << endl;
   cout << "The relative 1-distance between both datasets is " << relative << "." << endl;
   
   pDistanceAvgSubtracted(&orderedData1, &orderedData2, 1, &absolute, &relative);
   cout << "The average-subtracted absolute 1-distance between both datasets is " << absolute << "." << endl;
   cout << "The average-subtracted relative 1-distance between both datasets is " << relative << "." << endl;
   
   pDistance(&orderedData1, &orderedData2, 2, &absolute, &relative);
   cout << "The absolute 2-distance between both datasets is " << absolute << "." << endl;
   cout << "The relative 2-distance between both datasets is " << relative << "." << endl;
   pDistanceAvgSubtracted(&orderedData1, &orderedData2, 2, &absolute, &relative);
   cout << "The average-subtracted absolute 2-distance between both datasets is " << absolute << "." << endl;
   
   cout << "The average-subtracted relative 2-distance between both datasets is " << relative << "." << endl;
   infinityDistance(&orderedData1, &orderedData2, &absolute, &relative);
   cout << "The absolute infinity-distance between both datasets is " << absolute << "." << endl;
   cout << "The relative infinity-distance between both datasets is " << relative << "." << endl;
   
   infinityDistanceAvgSubtracted(&orderedData1, &orderedData2, &absolute, &relative);
   cout << "The average-subtracted absolute infinity-distance between both datasets is " << absolute << "." << endl;
   cout << "The average-subtracted relative infinity-distance between both datasets is " << relative << "." << endl;
   
   
   return 0;
}
