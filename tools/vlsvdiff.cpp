/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
*/

/*! \file vlsvdiff.cpp
 * \brief Tool to compare VLSV files.
 * 
 * Tool to compare two VLSV files, two folders with the same number of VLSV files or a folder to a reference file.
 * The tool assumes the files have a name syntax 'grid\.[0-9]+\.vlsv'. It takes four arguments.
 * 
 * Calling patterns are:
 * 
 * "$ vlsvdiff <file1> <file2> <Variable> <component>": Gives single-file statistics and distances between the two files given, for the variable and component given
 * 
 * "$ vlsvdiff <folder1> <folder2> <Variable> <component>": Gives single-file statistics and distances between pairs of files grid*.vlsv taken in alphanumeric order in the two folders given, for the variable and component given
 * 
 * "$ vlsvdiff <file1> <folder2> <Variable> <component>" or "$ vlsvdiff <folder1> <file2> <Variable> <component>": Gives single-file statistics and distances between a file, and files grid*.vlsv taken in alphanumeric order in the given folder, for the variable and component given
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
#include <typeinfo>
#include <algorithm>
#include <cstring>

#include "vlsvreader2.h"
#include "definitions.h"
#include "vlsv_reader.h"
#include "vlsvreaderinterface.h"
#include "vlsv_writer.h"

using namespace std;
using namespace vlsv;

static uint64_t convUInt(const char* ptr, const VLSV::datatype& dataType, const uint64_t& dataSize) {
   if (dataType != VLSV::UINT) {
      cerr << "Erroneous datatype given to convUInt" << endl;
      exit(1);
   }

   switch (dataSize) {
      case 1:
         return *reinterpret_cast<const unsigned char*> (ptr);
         break;
      case 2:
         return *reinterpret_cast<const unsigned short int*> (ptr);
         break;
      case 4:
         return *reinterpret_cast<const unsigned int*> (ptr);
         break;
      case 8:
         return *reinterpret_cast<const unsigned long int*> (ptr);
         break;
   }
   return 0;
}

static uint64_t convUInt(const char* ptr, const vlsv::datatype::type & dataType, const uint64_t& dataSize) {
   if (dataType != vlsv::datatype::type::UINT) {
      cerr << "Erroneous datatype given to convUInt" << endl;
      exit(1);
   }

   switch (dataSize) {
      case 1:
         return *reinterpret_cast<const unsigned char*> (ptr);
         break;
      case 2:
         return *reinterpret_cast<const unsigned short int*> (ptr);
         break;
      case 4:
         return *reinterpret_cast<const unsigned int*> (ptr);
         break;
      case 8:
         return *reinterpret_cast<const unsigned long int*> (ptr);
         break;
   }
   return 0;
}

/** Read given array data from input file, and byte-copy it to the output file.
 * @param input Input file reader.
 * @param output Output file reader.
 * @param tagName Name of the copied array.
 * @param inputAttributes XML attributes for the copied array.
 * @return If true, the array was copied successfully.*/
bool copyArray(vlsv::Reader& input,vlsv::Writer& output,
	       const std::string& tagName,
	       const list<pair<string,string> >& inputAttribs) {
   bool success = true;

   // Read input array attributes
   map<string,string> outputAttribs;
   if (input.getArrayAttributes(tagName,inputAttribs,outputAttribs) == false) {
      cerr << "ERROR: Failed to read array '" << tagName << "' attributes in " << __FILE__ << ":" << __LINE__ << endl;
      cerr << "Input attributes are:" << endl;
      for (list<pair<string,string> >::const_iterator it=inputAttribs.begin(); it!=inputAttribs.end(); ++it) {
	 cerr << "\t '" << it->first << "' = '" << it->second << "'" << endl;
      }
      return false;
   }

   // Figure out arraysize, vectorsize, datasize, and datatype of the copied array
   map<string,string>::const_iterator it;
   it = outputAttribs.find("arraysize"); if (it == outputAttribs.end()) return false;
   const uint64_t arraysize = atol(it->second.c_str());
   it = outputAttribs.find("vectorsize"); if (it == outputAttribs.end()) return false;
   const uint64_t vectorsize = atol(it->second.c_str());
   it = outputAttribs.find("datasize"); if (it == outputAttribs.end()) return false;
   const uint64_t datasize = atol(it->second.c_str());
   it = outputAttribs.find("datatype"); if (it == outputAttribs.end()) return false;
   const string datatype = it->second;

   const size_t bytes = arraysize*vectorsize*datasize;

   // Read values from input file
   char* ptr = new char[bytes];
   if (input.readArray(tagName,inputAttribs,0,arraysize,ptr) == false) {
      cerr << "ERROR: Failed to clone array '" << tagName << "' in " << __FILE__ << ":" << __LINE__ << endl;
      delete [] ptr; return false;
   }

   // Write array to output file
   if (output.writeArray(tagName,outputAttribs,datatype,arraysize,vectorsize,datasize,ptr) == false) {
      cerr << "ERROR: Failed to write array '" << tagName << "' in " << __FILE__ << ":" << __LINE__ << endl;
      success = false;
   }
   
   delete [] ptr; ptr = NULL;
   return success;
}

/** Copy the spatial mesh from input to output.
 * @param inputFileName Name of the input file where the mesh is copied from.
 * @param output VLSV reader for the file where the cloned mesh is written.
 * @param meshName Name of the mesh.
 * @return If true, the mesh was successfully cloned.*/
bool cloneMesh(const string& inputFileName,vlsv::Writer& output,const string& meshName) {
   bool success = true;
   
   vlsv::Reader input;
   if (input.open(inputFileName) == false) {
      cerr << "ERROR failed to open input file '" << inputFileName << "' in " << __FILE__ << ":" << __LINE__ << endl;
      return false;      
   }
   
   list<pair<string,string> > inputAttribs;
   inputAttribs.push_back(make_pair("name",meshName));   
   if (copyArray(input,output,"MESH",inputAttribs) == false) success = false;
   
   inputAttribs.clear();
   inputAttribs.push_back(make_pair("mesh",meshName));
   if (copyArray(input,output,"MESH_BBOX",inputAttribs) == false) success = false;
   if (copyArray(input,output,"MESH_DOMAIN_SIZES",inputAttribs) == false) success = false;
   if (copyArray(input,output,"MESH_NODE_CRDS_X",inputAttribs) == false) success = false;
   if (copyArray(input,output,"MESH_NODE_CRDS_Y",inputAttribs) == false) success = false;
   if (copyArray(input,output,"MESH_NODE_CRDS_Z",inputAttribs) == false) success = false;
   if (copyArray(input,output,"MESH_GHOST_LOCALIDS",inputAttribs) == false) success = false;
   if (copyArray(input,output,"MESH_GHOST_DOMAINS",inputAttribs) == false) success = false;

   input.close();
   return success;
}

/*! Extracts the dataset from the VLSV file opened by convertSILO.
 * \param vlsvReader oldVlsv::Reader class object used to access the VLSV file
 * \param meshName Address of the string containing the name of the mesh to be extracted
 * \param varToExtract Pointer to the char array containing the name of the variable to extract
 * \param compToExtract Unsigned int designating the component to extract (0 for scalars)
 * \param orderedData Pointer to the return argument map which will get the extracted dataset
 */
bool convertMesh(newVlsv::Reader& vlsvReader,
                 const string& meshName,
                 const char * varToExtract,
                 const uint compToExtract,
                 map<uint, Real> * orderedData,
		 unordered_map<size_t,size_t>& cellOrder,
		 const bool& storeCellOrder) {
   //Check for null pointer:
   if( !varToExtract || !orderedData ) {
      cerr << "ERROR, PASSED A NULL POINTER AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   bool meshSuccess = true;
   bool variableSuccess = true;
   
   datatype::type meshDataType;
   datatype::type variableDataType;
   uint64_t meshArraySize, meshVectorSize, meshDataSize;
   uint64_t variableArraySize, variableVectorSize, variableDataSize;

   //Get local cell ids:
   vector<uint64_t> local_cells;
   if( vlsvReader.getCellIds( local_cells ) == false ) {
      cerr << "Failed to read cell ids at "  << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   list<pair<string, string> > variableAttributes;
   const string _varToExtract( varToExtract );
   variableAttributes.push_back( make_pair("mesh", meshName) );
   variableAttributes.push_back( make_pair("name", _varToExtract) );
   //Read in array size, vector size, data type and data size of the array "VARIABLE" in the vlsv file (Needed in reading the array)
   if (vlsvReader.getArrayInfo("VARIABLE", variableAttributes, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false) {
      cerr << "ERROR, failed to read variable at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   //Check for correct output:
   if (local_cells.size() != variableArraySize) {
      cerr << "ERROR array size mismatch: " << local_cells.size() << " " << variableArraySize << endl;
   }
   if (compToExtract + 1 > variableVectorSize) {
      cerr << "ERROR invalid component, this variable has size " << variableVectorSize << endl;
      abort();
   }
   
   // Read the mesh array one node (of a spatial cell) at a time 
   // and create a map which contains each cell's CellID and variable to be extracted
   char* variableBuffer = new char[variableVectorSize*variableDataSize];
   float* variablePtrFloat = reinterpret_cast<float*>(variableBuffer);
   double* variablePtrDouble = reinterpret_cast<double*>(variableBuffer);
   uint* variablePtrUint = reinterpret_cast<uint*>(variableBuffer);
   int* variablePtrInt = reinterpret_cast<int*>(variableBuffer);

   if (storeCellOrder == true) {
      cellOrder.clear();
   }

   orderedData->clear();

   for (uint64_t i=0; i<local_cells.size(); ++i) {
      const short int amountToReadIn = 1;
      const uint64_t & startingReadIndex = i;
      if (vlsvReader.readArray("VARIABLE", variableAttributes, startingReadIndex, amountToReadIn, variableBuffer) == false) {
         cerr << "ERROR, failed to read variable at " << __FILE__ << " " << __LINE__ << endl;
         variableSuccess = false; 
         break;
      }
      // Get the CellID
      uint64_t & CellID = local_cells[i];
      
      // Get the variable value
      Real extract = NAN;

      switch (variableDataType) {
         case datatype::type::FLOAT:
            if(variableDataSize == sizeof(float)) extract = (Real)(variablePtrFloat[compToExtract]);
            if(variableDataSize == sizeof(double)) extract = (Real)(variablePtrDouble[compToExtract]);
            break;
         case datatype::type::UINT:
            extract = (Real)(variablePtrUint[compToExtract]);
            break;
         case datatype::type::INT:
            extract = (Real)(variablePtrInt[compToExtract]);
            break;
         case datatype::type::UNKNOWN:
            cerr << "ERROR, BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << endl;
            break;
      }
      // Put those into the map
      orderedData->insert(pair<uint64_t, Real>(CellID, extract));
      if (storeCellOrder == true) {
	 cellOrder[CellID] = i;
      }
   }

   if (meshSuccess == false) {
      cerr << "ERROR reading array MESH" << endl;
   }
   if (variableSuccess == false) {
      cerr << "ERROR reading array VARIABLE " << varToExtract << endl;
   }
   return meshSuccess && variableSuccess;
}

/*! Opens the VLSV file and extracts the mesh names. Sends for processing to convertMesh.
 * \param fileName String containing the name of the file to be processed
 * \param varToExtract Pointer to the char array containing the name of the variable to extract
 * \param compToExtract Unsigned int designating the component to extract (0 for scalars)
 * \param orderedData Pointer to the return argument map which will get the extracted dataset
 * \sa convertMesh
 */
template <class T>
bool convertSILO(const string fileName,
                 const char * varToExtract,
                 const uint compToExtract,
                 map<uint, Real> * orderedData,
		 unordered_map<size_t,size_t>& cellOrder,
		 const bool& storeCellOrder=false) 
{
   bool success = true;

   // Open VLSV file for reading:
   T vlsvReader;
   if (vlsvReader.open(fileName) == false) {
      cerr << "Failed to open '" << fileName << "'" << endl;
      return false;
   }

   // Get the names of all meshes in vlsv file
   list<string> meshNames;
   if (vlsvReader.getMeshNames(meshNames) == false) {
      cerr << "Failed to read mesh names" << endl;
      exit(1);
   }
   
   // Clear old data
   orderedData->clear();

   for (list<string>::const_iterator it=meshNames.begin(); it!=meshNames.end(); ++it) {
      if (convertMesh(vlsvReader, *it, varToExtract, compToExtract, orderedData, cellOrder, storeCellOrder) == false) {
         return false;
      }
   }
   vlsvReader.close();
   return success;
}

/*! Shift the second file to the average of the first
 * \param orderedData1 Pointer to the reference file's data
 * \param orderedData2 Pointer to the data to be shifted
 * \param shiftedData2 Pointer to where the shifted data of the second file will be put
 */
bool shiftAverage(const map<uint, Real> * const orderedData1,
                  const map<uint, Real> * const orderedData2,
                  map<uint, Real> * shiftedData2
                  ) {
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

/*! Compute the absolute and relative \f$ p \f$-distance between two datasets X(x) provided in the maps orderedData1 and orderedData2. Note that the dataset passed in orderedData1 will be taken as the reference dataset both when shifting averages and when computing relative distances.
 * 
 * For \f$ p \neq 0 \f$:
 * 
 * absolute \f$ p \f$-distance defined as:
 * 
 * \f$ \|X_1 - X_2\|_p = \left[\sum_i |X_1(i) - X_2(i)|^p\right]^{1/p}\f$,
 * 
 * relative \f$ p \f$-distance defined as:
 * 
 * \f$ \|X_1 - X_2\|_p = \left[\sum_i |X_1(i) - X_2(i)|^p\right]^{1/p} / \|X_1\|_p \f$.
 * 
 * For \f$ p = 0 \f$ it is the \f$ \infty \f$-distance:
 * 
 * absolute \f$ \infty \f$-distance defined as:
 * 
 * \f$ \|X_1 - X_2\|_\infty = \max_i\left(|X_1(i) - X_2(i)|\right)\f$
 * 
 * relative \f$ \infty \f$-distance defined as:
 * 
 * \f$ \|X_1 - X_2\|_\infty = \max_i\left(|X_1(i) - X_2(i)|\right) / \|X_1\|_\infty \f$
 * 
 * \param orderedData1 Pointer to the first file's data map
 * \param orderedData2 Pointer to the second file's data map
 * \param p Parameter of the distance formula
 * \param absolute Return argument pointer, absolute value
 * \param relative Return argument pointer, relative value
 * \param doShiftAverage Boolean argument to determine whether to shift the second file's data
 * \sa shiftAverage
 */
bool pDistance(const map<uint, Real>& orderedData1,
               const map<uint, Real>& orderedData2,
               creal p,
               Real * absolute,
               Real * relative,
               const bool doShiftAverage,
	       const unordered_map<size_t,size_t>& cellOrder,
	       vlsv::Writer& outputFile,
	       const std::string& meshName,
	       const std::string& varName
               ) {
   #warning shift average not re-implemented
   
   // Reset old values
   *absolute = 0.0;
   *relative = 0.0;

   vector<Real> array(orderedData1.size());
   for (size_t i=0; i<array.size(); ++i) array[i] = -1.0;

   Real length = 0.0;
   if (p == 0) {
      for (map<uint,Real>::const_iterator it1=orderedData1.begin(); it1!=orderedData1.end(); ++it1) {
	 map<uint,Real>::const_iterator it2 = orderedData2.find(it1->first);
	 Real value = 0.0;
	 if (it2 != orderedData2.end()) {
	    value = abs(it1->second - it2->second);
	    *absolute = max(*absolute, value);
	    length    = max(length, abs(it1->second));
	 }
	 array[cellOrder.at(it1->first)] = value;
      }
   } else {
      for (map<uint,Real>::const_iterator it1=orderedData1.begin(); it1!=orderedData1.end(); ++it1) {
	 map<uint,Real>::const_iterator it2 = orderedData2.find(it1->first);
	 Real value = 0.0;
	 if (it2 != orderedData2.end()) {
	    value = pow(abs(it1->second - it2->second), p);
	    *absolute += value;
	    length    += pow(abs(it1->second), p);
	 }
	 *absolute = pow(*absolute, 1.0 / p);
	 length = pow(length, 1.0 / p);
	 array[cellOrder.at(it1->first)] = value;
      }
   }

   if (length != 0.0) *relative = *absolute / length;
   else {
      cout << "WARNING (pDistance) : length of reference is 0.0, cannot divide to give relative distance." << endl;
      *relative = -1;
   }

   map<string,string> attributes;
   attributes["mesh"] = meshName;
   attributes["name"] = varName;
   if (outputFile.writeArray("VARIABLE",attributes,array.size(),1,&(array[0])) == false) {
      cerr << "ERROR failed to write variable '" << varName << "' to output file in " << __FILE__ << ":" << __LINE__ << endl;
      return 1;
   }
   
   return 0;
}

/*! In verbose mode print the distance, in non-verbose store them for later output when lastCall is true
 * \param p Parameter of the distance
 * \param absolute Absolute value pointer
 * \param relative Relative value pointer
 * \param shiftedAverage Boolean parameter telling whether the dataset is average-shifted
 * \param verboseOutput Boolean parameter telling whether the output is verbose or compact
 * \param lastCall Boolean parameter telling whether this is the last call to the function
 * \sa shiftAverage pDistance
 */
bool outputDistance(const Real p,
                    const Real * absolute,
                    const Real * relative,
                    const bool shiftedAverage,
                    const bool verboseOutput,
                    const bool lastCall
)
{
   if(verboseOutput == true) {
      if(shiftedAverage == false) {
         cout << "The absolute " << p << "-distance between both datasets is " << *absolute  << endl;
         cout << "The relative " << p << "-distance between both datasets is " << *relative  << endl;
      } else {
         cout << "The average-shifted absolute " << p << "-distance between both datasets is " << *absolute  << endl;
         cout << "The average-shifted relative " << p << "-distance between both datasets is " << *relative  << endl;
      }
   } else {
      static vector<Real> fileOutputData;
      static uint fileNumber = 0;
      
      if(lastCall == true) {
         vector<Real>::const_iterator it;
         for(it = fileOutputData.begin(); it != fileOutputData.end(); it++) {
            cout << *it << "\t";
         }
         fileOutputData.clear();
         return 0;
      }
      
      fileOutputData.push_back(*absolute);
      fileOutputData.push_back(*relative);
   }
   return 0;
}

/*! Compute statistics on a single file
 * \param size Return argument pointer, dataset size
 * \param mini Return argument pointer, dataset minimum
 * \param maxi Return argument pointer, dataset maximum
 * \param avg Return argument pointer, dataset average
 * \param stdev Return argument pointer, dataset standard deviation
 */
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

/*! In verbose mode print the statistics, in non-verbose store them for later output when lastCall is true
 * \param size Pointer to dataset size
 * \param mini Pointer to dataset minimum
 * \param maxi Pointer to dataset maximum
 * \param avg Pointer to dataset average
 * \param stdev Pointer to dataset standard deviation
 * \param verboseOutput Boolean parameter telling whether the output is verbose or compact
 * \param lastCall Boolean parameter telling whether this is the last call to the function
 * \sa singleStatistics
 */
bool outputStats(const Real * size,
                 const Real * mini,
                 const Real * maxi,
                 const Real * avg,
                 const Real * stdev,
                 const bool verboseOutput,
                 const bool lastCall
                 ) {
   if(verboseOutput == true)
   {
      cout << "Statistics on file: size " << *size
      << " min = " << *mini
      << " max = " << *maxi
      << " average = " << *avg
      << " standard deviation " << *stdev
      << endl;
   }
   else
   {
      static uint fileNumber = 0;
      static vector<Real> pairStats;
      
      if(lastCall == true)
      {
         vector<Real>::const_iterator it;
         for(it = pairStats.begin(); it != pairStats.end(); it++)
         {
            cout << *it << "\t";
         }
         pairStats.clear();
         return 0;
      }
      
      if(fileNumber%2 == 0)
      {
         pairStats.push_back(fileNumber / 2 + 1);
      }
      pairStats.push_back(*size);
      pairStats.push_back(*mini);
      pairStats.push_back(*maxi);
      pairStats.push_back(*avg);
      pairStats.push_back(*stdev);
      fileNumber++;
   }
   return 0;
}

/*! In folder-processing, non-verbose mode the data are stored during the processing and output at the end to have the data sorted properly
 * \sa outputStats outputDistance
 */
bool printNonVerboseData()
{
   static bool header = true;
   if(header == true)
   {
      // Key to contents
      cout << "#1   File number in folder\n" <<
              "#2   File 1 size\n" <<
              "#3   File 1 min\n" <<
              "#4   File 1 max\n" <<
              "#5   File 1 average\n" <<
              "#6   File 1 standard deviation\n" <<
              "#7   File 2 size\n" <<
              "#8   File 2 min\n" <<
              "#9   File 2 max\n" <<
              "#10  File 2 average\n" <<
              "#11  File 2 standard deviation\n" <<
              "#12  absolute infinity-distance\n" <<
              "#13  relative infinity-distance\n" <<
              "#14  absolute average-shifted infinity-distance\n" <<
              "#15  relative average-shifted infinity-distance\n" <<
              "#16  absolute 1-distance\n" <<
              "#17  relative 1-distance\n" <<
              "#18  absolute average-shifted 1-distance\n" <<
              "#19  relative average-shifted 1-distance\n" <<
              "#20  absolute 2-distance\n" <<
              "#21  relative 2-distance\n" <<
              "#22  absolute average-shifted 2-distance\n" <<
              "#23  relative average-shifted 2-distance\n" <<
              endl;
      header = false;
   }
   
   // Data
   // last argument (lastCall) is true to get the output of the whole stored dataset
   outputStats(NULL, NULL, NULL, NULL, NULL, false, true);
   outputDistance(0, NULL, NULL, false, false, true);
   
   return 0;
}

bool getBlockIds( newVlsv::Reader & vlsvReader,
                  const unordered_map<uint64_t, pair<uint64_t, uint32_t>> & cellsWithBlocksLocations,
                  const uint64_t & cellId,
                  vector<uint32_t> & blockIds ) {
   // Read the block ids:
   //Check if the cell id can be found:
   unordered_map<uint64_t, pair<uint64_t, uint32_t>>::const_iterator it = cellsWithBlocksLocations.find( cellId );
   if( it == cellsWithBlocksLocations.end() ) {
      cerr << "COULDNT FIND CELL ID " << cellId << " AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Get offset and number of blocks:
   pair<uint64_t, uint32_t> offsetAndBlocks = it->second;
   const uint64_t blockOffset = get<0>(offsetAndBlocks);
   const uint32_t N_blocks = get<1>(offsetAndBlocks);

   // Get some required info from VLSV file:
   list<pair<string, string> > attribs;
   attribs.push_back(make_pair("mesh", "SpatialGrid"));

   //READ BLOCK IDS:
   uint64_t blockIds_arraySize, blockIds_vectorSize, blockIds_dataSize;
   vlsv::datatype::type blockIds_dataType;
   //Input blockIds_arraySize, blockIds_vectorSize, blockIds_dataSize blockIds_dataType: (Returns false if fails)
   if (vlsvReader.getArrayInfo("BLOCKIDS", attribs, blockIds_arraySize, blockIds_vectorSize, blockIds_dataType, blockIds_dataSize) == false) {
      cerr << "ERROR, COULD NOT FIND BLOCKIDS AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Make sure blockid's datatype is correct:
   if( blockIds_dataType != vlsv::datatype::type::UINT ) {
      cerr << "ERROR, bad datatype at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Create buffer for reading in data:  (Note: arraySize, vectorSize, etc were fetched from getArrayInfo)
   char * blockIds_buffer = new char[N_blocks*blockIds_vectorSize*blockIds_dataSize];
   //Read the data into the buffer:
   if( vlsvReader.readArray( "BLOCKIDS", attribs, blockOffset, N_blocks, blockIds_buffer ) == false ) {
      cerr << "ERROR, FAILED TO READ BLOCKIDS AT " << __FILE__ << " " << __LINE__ << endl;
      delete[] blockIds_buffer;
      return false;
   }
   //Input the block ids:
   blockIds.reserve(N_blocks);
   for (uint64_t i = 0; i < N_blocks; ++i) {
      const uint64_t blockId = convUInt(blockIds_buffer + i*blockIds_dataSize, blockIds_dataType, blockIds_dataSize);
      blockIds.push_back( (uint32_t)(blockId) );
   }
   delete[] blockIds_buffer;
   return true;

}

uint32_t getBlockId( const double vx,
                     const double vy,
                     const double vz,
                     const double dvx,
                     const double dvy,
                     const double dvz,
                     const double vx_min,
                     const double vy_min,
                     const double vz_min,
                     const double vx_length,
                     const double vy_length,
                     const double vz_length ) {

   const array<unsigned int, 3> indices{ { (unsigned int) floor((vx - vx_min) / (double)(dvx*4)),
                                     (unsigned int) floor((vy - vy_min) / (double)(dvy*4)),
                                     (unsigned int) floor((vz - vz_min) / (double)(dvz*4)) } };
   const uint32_t blockId = indices[0]
                + indices[1] * vx_length
                + indices[2] * vx_length * vy_length;

    return blockId;
}



//Loads a parameter from a file
//usage: Real x = loadParameter( vlsvReader, nameOfParameter );
//Note: this is used in getCellIdFromCoords
//Input:
//[0] vlsvReader -- some VLSVReader which has a file opened
//[1] name -- name of the parameter, e.g. "xmin"
//Output:
//[0] Parameter -- Saves the parameter into the parameter variable
template <typename T>
bool loadParameter( VLSVReader& vlsvReader, const string& name, T & parameter ) {
   //Declare dataType, arraySize, vectorSize, dataSize so we know how much data we want to store
   VLSV::datatype dataType;
   uint64_t arraySize, vectorSize, dataSize; //vectorSize should be 1
   //Write into dataType, arraySize, etc with getArrayInfo -- if fails to read, give error message
   if( vlsvReader.getArrayInfo( "PARAMETERS", name, arraySize, vectorSize, dataType, dataSize ) == false ) {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      return false;
   }
   //Declare a buffer to write the parameter's data in (arraySize, etc was received from getArrayInfo)
   char * buffer = new char[arraySize * vectorSize * dataSize];
  
   //Read data into the buffer and return error if something went wrong
   if( vlsvReader.readArray( "PARAMETERS", name, 0, vectorSize, buffer ) == false ) {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      return false;
   }
   //SHOULD be a vector of size 1 and since I am going to want to assume that, making a check here
   if( vectorSize != 1 ) {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      return false;
   }
   //Input the parameter
   if( typeid(T) == typeid(double) ) {
      if( dataSize == 8 ) {
         parameter = *reinterpret_cast<double*>(buffer);
      } else if( dataSize == 4 ) {
         parameter = *reinterpret_cast<float*>(buffer);
      } else {
         cerr << "Error, bad datasize while reading parameters at " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
   } else if( typeid(T) == typeid(uint64_t) ) {
      if( dataSize == 8 ) {
         parameter = *reinterpret_cast<uint64_t*>(buffer);
      } else if( dataSize == 4 ) {
         parameter = *reinterpret_cast<uint32_t*>(buffer);
      } else {
         cerr << "Error, bad datasize while reading parameters at " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
   } else {
      cerr << "Error, could not read parameter '" << name << "' at: " << __FILE__ << " " << __LINE__; //FIX
      cerr << " Error message: invalid type in loadParameters" << endl;
      return false;
   }
   return true;
}



bool getBlockIds( oldVlsv::Reader & vlsvReader,
                  const unordered_map<uint64_t, pair<uint64_t, uint32_t>> & cellsWithBlocksLocations,
                  const uint64_t & cellId,
                  vector<uint32_t> & blockIds ) {
   // Read the block ids:
   //Check if the cell id can be found:
   unordered_map<uint64_t, pair<uint64_t, uint32_t>>::const_iterator it = cellsWithBlocksLocations.find( cellId );
   if( it == cellsWithBlocksLocations.end() ) {
      cerr << "COULDNT FIND CELL ID " << cellId << " AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Get offset and number of blocks:
   pair<uint64_t, uint32_t> offsetAndBlocks = it->second;
   const uint64_t blockOffset = get<0>(offsetAndBlocks);
   const uint32_t N_blocks = get<1>(offsetAndBlocks);

   // Get some required info from VLSV file:
   list<pair<string, string> > attribs;
   attribs.push_back(make_pair("name", "SpatialGrid"));

   //READ BLOCK COORDINATES:
   uint64_t blockCoordinates_arraySize, blockCoordinates_vectorSize, blockCoordinates_dataSize;
   vlsv::datatype::type blockCoordinates_dataType;
   //Input blockCoordinates_arraySize, blockCoordinates_vectorSize, blockCoordinates_dataSize blockCoordinates_dataType: (Returns false if fails)
   if (vlsvReader.getArrayInfo("BLOCKCOORDINATES", attribs, blockCoordinates_arraySize, blockCoordinates_vectorSize, blockCoordinates_dataType, blockCoordinates_dataSize) == false) {
      cerr << "ERROR, COULD NOT FIND BLOCKCOORDINATES AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Make sure blockid's datatype is correct:
   if( blockCoordinates_dataType != vlsv::datatype::type::FLOAT ) {
      cerr << "ERROR, bad datatype at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Create buffer for reading in data:  (Note: arraySize, vectorSize, etc were fetched from getArrayInfo)
   char * blockCoordinates_buffer = new char[N_blocks*blockCoordinates_vectorSize*blockCoordinates_dataSize];
   //Read the data into the buffer:
   if( vlsvReader.readArray( "BLOCKCOORDINATES", attribs, blockOffset, N_blocks, blockCoordinates_buffer ) == false ) {
      cerr << "ERROR, FAILED TO READ BLOCKCOORDINATES AT " << __FILE__ << " " << __LINE__ << endl;
      delete[] blockCoordinates_buffer;
      return false;
   }
   //Input the block ids:
   // First read in vxmin, vymin, vzmin etc
   double vx_min, vy_min, vz_min;
   uint64_t vx_length, vy_length, vz_length;
   loadParameter( vlsvReader, "vxmin", vx_min );
   loadParameter( vlsvReader, "vymin", vy_min );
   loadParameter( vlsvReader, "vzmin", vz_min );
   loadParameter( vlsvReader, "vxblocks_ini", vx_length );
   loadParameter( vlsvReader, "vyblocks_ini", vy_length );
   loadParameter( vlsvReader, "vzblocks_ini", vz_length );

   double vx, vy, vz, dvx, dvy, dvz;

   blockIds.reserve(N_blocks);
   for (uint64_t b = 0; b < N_blocks; ++b) {
      if (blockCoordinates_dataSize == 4) {
         vx = *reinterpret_cast<float*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 0 * blockCoordinates_dataSize);
         vy = *reinterpret_cast<float*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 1 * blockCoordinates_dataSize);
         vz = *reinterpret_cast<float*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 2 * blockCoordinates_dataSize);
         dvx = *reinterpret_cast<float*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 3 * blockCoordinates_dataSize);
         dvy = *reinterpret_cast<float*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 4 * blockCoordinates_dataSize);
         dvz = *reinterpret_cast<float*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 5 * blockCoordinates_dataSize);
      } else {
         vx = *reinterpret_cast<double*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 0 * blockCoordinates_dataSize);
         vy = *reinterpret_cast<double*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 1 * blockCoordinates_dataSize);
         vz = *reinterpret_cast<double*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 2 * blockCoordinates_dataSize);
         dvx = *reinterpret_cast<double*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 3 * blockCoordinates_dataSize);
         dvy = *reinterpret_cast<double*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 4 * blockCoordinates_dataSize);
         dvz = *reinterpret_cast<double*> (blockCoordinates_buffer + b * blockCoordinates_vectorSize * blockCoordinates_dataSize + 5 * blockCoordinates_dataSize);
      }
      const uint32_t blockId = getBlockId( vx, vy, vz, dvx, dvy, dvz, vx_min, vy_min, vz_min, vx_length, vy_length, vz_length);
      blockIds.push_back( (uint32_t)(blockId) );
   }
   delete[] blockCoordinates_buffer;
   return true;

}


// Reads avgs values of some given cell id
// Input:
// [0] vlsvReader -- Some vlsv reader with a file open
// [1] cellId -- The spatial cell's ID
// Output:
// [2] avgs -- Saves the output into an unordered map with block id as the key and an array of avgs as the value
// return false or true depending on whether the operation was successful
template <class T>
bool readAvgs( T & vlsvReader,
               const unordered_map<uint64_t, pair<uint64_t, uint32_t>> & cellsWithBlocksLocations,
               const uint64_t & cellId, 
               unordered_map<uint32_t, array<double, 64> > & avgs ) {
   // Get the block ids:
   vector<uint32_t> blockIds;
   if( getBlockIds( vlsvReader, cellsWithBlocksLocations, cellId, blockIds ) == false ) { return false; }
   // Read avgs:
   list<pair<string, string> > attribs;
   attribs.push_back(make_pair("name", "avgs"));
   attribs.push_back(make_pair("mesh", "SpatialGrid"));

   datatype::type dataType;
   uint64_t arraySize, vectorSize, dataSize;
   if (vlsvReader.getArrayInfo("BLOCKVARIABLE", attribs, arraySize, vectorSize, dataType, dataSize) == false) {
      cerr << "ERROR READING BLOCKVARIABLE AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   // Make a routine error checks:
   if( vectorSize != 64 ) {
      cerr << "ERROR, BAD AVGS VECTOR SIZE AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   unordered_map<uint64_t, pair<uint64_t, uint32_t>>::const_iterator it = cellsWithBlocksLocations.find( cellId );
   if( it == cellsWithBlocksLocations.end() ) {
      cerr << "COULDNT FIND CELL ID " << cellId << " AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   //Get offset and number of blocks:
   pair<uint64_t, uint32_t> offsetAndBlocks = it->second;
   const uint64_t blockOffset = get<0>(offsetAndBlocks);
   const uint32_t N_blocks = get<1>(offsetAndBlocks);

   if( N_blocks != blockIds.size() ) {
      cerr << "ERROR, BAD AVGS ARRAY SIZE AT " << __FILE__ << " " << __LINE__ << endl;
      cerr << "AVGS SIZE: " << N_blocks << endl;
      cerr << "BLOCKIDS SIZE: " << blockIds.size() << endl;
      return false;
   }

   char* buffer = new char[N_blocks * vectorSize * dataSize];
   if (vlsvReader.readArray("BLOCKVARIABLE", attribs, blockOffset, N_blocks, buffer) == false) {
      cerr << "ERROR could not read block variable at " << __FILE__ << " " << __LINE__ << endl;
      delete[] buffer;
      return false;
   }
   // Input avgs values:
   array<double, 64> avgs_temp;
   if( dataSize == 4 ) {
      float * buffer_float = reinterpret_cast<float*>( buffer );
      for( uint b = 0; b < blockIds.size(); ++b ) {
         const uint32_t & blockId = blockIds[b];
         for( uint i = 0; i < vectorSize; ++i ) {
            avgs_temp[i] = buffer_float[vectorSize * b + i];
         }
         avgs.insert(make_pair(blockId, avgs_temp));
      }
   } else if( dataSize == 8 ) {
      double * buffer_double = reinterpret_cast<double*>( buffer );
      for( uint b = 0; b < blockIds.size(); ++b ) {
         const uint32_t & blockId = blockIds[b];
         for( uint i = 0; i < vectorSize; ++i ) {
            avgs_temp[i] = buffer_double[vectorSize * b + i];
         }
         avgs.insert(make_pair(blockId, avgs_temp));
      }
   } else {
      cerr << "ERROR, BAD AVGS DATASIZE AT " << __FILE__ << " " << __LINE__ << endl;
      delete[] buffer;
      return false;
   }
   delete[] buffer;
   return true;
}

template <class T>
bool getCellsWithBlocksLocations( T & vlsvReader, 
                                  unordered_map<uint64_t, pair<uint64_t, uint32_t>> & cellsWithBlocksLocations ) {
   if(cellsWithBlocksLocations.empty() == false) {
      cellsWithBlocksLocations.clear();
   }
   const string meshName = "SpatialGrid";
   vlsv::datatype::type cwb_dataType;
   uint64_t cwb_arraySize, cwb_vectorSize, cwb_dataSize;
   list<pair<string, string> > attribs;

   //Get the mesh name for reading in data from the correct place
   if( typeid(vlsvReader) == typeid(oldVlsv::Reader) ) {
      attribs.push_back( make_pair("name", meshName) );
   } else {
      attribs.push_back( make_pair("mesh", meshName) );
   }

   //Get array info
   if (vlsvReader.getArrayInfo("CELLSWITHBLOCKS", attribs, cwb_arraySize, cwb_vectorSize, cwb_dataType, cwb_dataSize) == false) {
      cerr << "ERROR, COULD NOT FIND ARRAY CELLSWITHBLOCKS AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   //Make sure the data format is correct:
   if( cwb_vectorSize != 1 ) {
      cerr << "ERROR, BAD VECTORSIZE AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   if( cwb_dataType != vlsv::datatype::type::UINT ) {
      cerr << "ERROR, BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   if( cwb_dataSize != sizeof(uint64_t) ) {
      cerr << "ERROR, BAD DATASIZE AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   // Create buffer and read data:
   const uint64_t cwb_amountToReadIn = cwb_arraySize * cwb_vectorSize * cwb_dataSize;
   const uint16_t cwb_startingPoint = 0;
   char * cwb_buffer = new char[cwb_amountToReadIn];
   if (vlsvReader.readArray("CELLSWITHBLOCKS", attribs, cwb_startingPoint, cwb_arraySize, cwb_buffer) == false) {
      cerr << "Failed to read block metadata for mesh '" << meshName << "'" << endl;
      delete[] cwb_buffer;
      return false;
   }

   vlsv::datatype::type nb_dataType;
   uint64_t nb_arraySize, nb_vectorSize, nb_dataSize;  

   //Get the mesh name for reading in data from the correct place
   //Read array info -- stores output in nb_arraySize, nb_vectorSize, nb_dataType, nb_dataSize
   if (vlsvReader.getArrayInfo("BLOCKSPERCELL", attribs, nb_arraySize, nb_vectorSize, nb_dataType, nb_dataSize) == false) {
      cerr << "ERROR, COULD NOT FIND ARRAY BLOCKSPERCELL AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   // Create buffers for  number of blocks (nb) and read data:
   const short int startingPoint = 0; //Read the array from 0 (the beginning)
   char* nb_buffer = new char[nb_arraySize * nb_vectorSize * nb_dataSize];
   if (vlsvReader.readArray("BLOCKSPERCELL", attribs, startingPoint, nb_arraySize, nb_buffer) == false) {
      cerr << "Failed to read number of blocks for mesh '" << meshName << "'" << endl;
      delete[] nb_buffer;
      delete[] cwb_buffer;
      return false;
   }

   // Input cellswithblock locations:
   uint64_t blockOffset = 0;
   uint64_t N_blocks;
   for (uint64_t cell = 0; cell < cwb_arraySize; ++cell) {
      const uint64_t readCellID = convUInt(cwb_buffer + cell*cwb_dataSize, cwb_dataType, cwb_dataSize);
      N_blocks = convUInt(nb_buffer + cell*nb_dataSize, nb_dataType, nb_dataSize);
      const pair<uint64_t, uint32_t> input = make_pair( blockOffset, N_blocks );
      //Insert the location and number of blocks into the map
      cellsWithBlocksLocations.insert( make_pair(readCellID, input) );
      blockOffset += N_blocks;
   }

   delete[] cwb_buffer;
   delete[] nb_buffer;
   return true;
}

template <class T, class U>
bool compareAvgs( const string fileName1,
                  const string fileName2,
                  const bool verboseOutput,
                  vector<uint64_t> & cellIds1,
                  vector<uint64_t> & cellIds2
                ) {
   if( cellIds1.empty() == true || cellIds2.empty() == true ) {
      cerr << "ERROR, CELL IDS EMPTY IN COMPARE AVGS" << endl;
      return false;
   }
   // Declare map for locating velocity spaces within cell ids
   // Note: Key = cell id, value->first = blockOffset, value->second = numberOfBlocksToRead
   unordered_map<uint64_t, pair<uint64_t, uint32_t>> cellsWithBlocksLocations1;
   unordered_map<uint64_t, pair<uint64_t, uint32_t>> cellsWithBlocksLocations2;
   // Open the files for reading:
   T vlsvReader1;
   if( vlsvReader1.open(fileName1) == false ) {
      cerr << "Error opening file name " << fileName1 << " at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   U vlsvReader2;
   if( vlsvReader2.open(fileName2) == false ) {
      cerr << "Error opening file name " << fileName2 << " at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   if( getCellsWithBlocksLocations( vlsvReader1, cellsWithBlocksLocations1 ) == false ) {
      cerr << "ERROR AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   if( getCellsWithBlocksLocations( vlsvReader2, cellsWithBlocksLocations2 ) == false ) {
      cerr << "ERROR AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   // Consistency check:
   if( cellsWithBlocksLocations2.size() != cellsWithBlocksLocations1.size() ) {
      cerr << "BAD CELLS WITH BLOCKS SIZE AT "  << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   // Create a few variables for the cell id loop:
   vector< double > avgsDiffs;
   double totalAbsAvgs = 0;
   double totalAbsDiff = 0;
   double totalAbsLog10Diff = 0;
   double threshold=1e-16;
   uint64_t numOfRelevantCells = 0;
   uint64_t numOfIdenticalBlocks = 0;
   uint64_t numOfNonIdenticalBlocks = 0;
   if( cellIds1[0] == 0 || cellIds2[0] == 0 ) {
      // User input 0 as the cell id -- compare all cell ids
      cellIds1.clear();
      cellIds2.clear();
      for( unordered_map<uint64_t, pair<uint64_t, uint32_t>>::const_iterator it = cellsWithBlocksLocations1.begin(); it != cellsWithBlocksLocations1.end(); ++it ) {
         cellIds1.push_back(it->first);
         cellIds2.push_back(it->first);
      }
   }

   if( cellIds1.size() != cellIds2.size() ) {
      cerr << "ERROR, BAD CELL ID SIZES AT " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }
   // Go through cell ids:
   for( uint cellIndex = 0; cellIndex < cellIds2.size(); cellIndex++ ) {
      const uint64_t & cellId1 = cellIds1[cellIndex];
      const uint64_t & cellId2 = cellIds2[cellIndex];
      // Get the avgs in a hash map (The velocity block id is the key and avgs is the value):
      const uint velocityCellsPerBlock = 64;
      unordered_map<uint32_t, array<double, velocityCellsPerBlock> > avgs1;
      unordered_map<uint32_t, array<double, velocityCellsPerBlock> > avgs2;
      // Store the avgs in avgs1 and 2:
      if( readAvgs( vlsvReader1, cellsWithBlocksLocations1, cellId1, avgs1 ) == false ) {
         cerr << "ERROR, FAILED TO READ AVGS AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
      if( readAvgs( vlsvReader2, cellsWithBlocksLocations2, cellId2, avgs2 ) == false ) {
         cerr << "ERROR, FAILED TO READ AVGS AT " << __FILE__ << " " << __LINE__ << endl;
         return false;
      }
   
      //Compare the avgs values:
      // First make a check on how many of the block ids are identical:
      const size_t sizeOfAvgs1 = avgs1.size();
      const size_t sizeOfAvgs2 = avgs2.size();
      // Vector of block ids that are the same
      vector<uint32_t> blockIds1;
      vector<uint32_t> blockIds2;
      blockIds1.reserve(sizeOfAvgs1);
      blockIds2.reserve(sizeOfAvgs2);
      // Input block ids:
      for( unordered_map<uint32_t, array<double, velocityCellsPerBlock> >::const_iterator it = avgs1.begin(); it != avgs1.end(); ++it ) {
         blockIds1.push_back(it->first);
      }
      for( unordered_map<uint32_t, array<double, velocityCellsPerBlock> >::const_iterator it = avgs2.begin(); it != avgs2.end(); ++it ) {
         blockIds2.push_back(it->first);
      }
      // Compare block ids:
      // Sort
      sort( blockIds1.begin(), blockIds1.end() );
      sort( blockIds2.begin(), blockIds2.end() );
      // Create iterators
      vector<uint32_t>::const_iterator it1 = blockIds1.begin();
      vector<uint32_t>::const_iterator it2 = blockIds2.begin();
      // Separate block ids into two categories -- the ones that blockids1 and blockids2 share and ones that only one of them shares
      vector<uint32_t> identicalBlockIds;
      vector<uint32_t> nonIdenticalBlockIds;
   
      while( true ) {
         if( it1 == blockIds1.end() || it2 == blockIds2.end() ) {
            // Reach end of block ids
            break;
         }
         if( *it1 == *it2 ) {
            // Identical block id
            identicalBlockIds.push_back(*it1);
            it1++; it2++;
         } else if( *it1 < *it2 ) {
            // Non identical block id
            // The block ids are sorted so to get identical block ids one must increment the lower value
            nonIdenticalBlockIds.push_back(*it1);
            it1++;
         } else if( *it2 < *it1 ) {
            // Non identical block id
            // The block ids are sorted so to get identical block ids one must increment the lower value
            nonIdenticalBlockIds.push_back(*it2);
            it2++;
         }
      }
      // Get the rest of the non identical block ids (If there are any)
      // Note: This is only needed if for example it1 hit the end of the iteration and it2 still isn't at the end
      for( ; it1 != blockIds1.end(); ++it1 ) {
         nonIdenticalBlockIds.push_back(*it1);
      }
      for( ; it2 != blockIds2.end(); ++it2 ) {
         nonIdenticalBlockIds.push_back(*it2);
      }
      // Compare block ids:
      const uint64_t totalNumberOfBlocks = identicalBlockIds.size() + nonIdenticalBlockIds.size();
      const double percentageOfIdenticalBlocks = (double)(totalNumberOfBlocks) / (double)(identicalBlockIds.size());
      // Compare the avgs values of the identical blocks:
      avgsDiffs.reserve(avgsDiffs.size() + identicalBlockIds.size() * velocityCellsPerBlock);
      for( vector<uint32_t>::const_iterator it = identicalBlockIds.begin(); it != identicalBlockIds.end(); ++it ) {
         // Get the block id
         const uint32_t blockId = *it;
         // Get avgs values:
         const array<double, velocityCellsPerBlock> & avgsValues1 = avgs1.at(blockId);
         const array<double, velocityCellsPerBlock> & avgsValues2 = avgs2.at(blockId);
         // Get the diff:
         for( uint i = 0; i < velocityCellsPerBlock; ++i ) {
            double val1=avgsValues1[i]>threshold?avgsValues1[i]:threshold;
            double val2=avgsValues2[i]>threshold?avgsValues2[i]:threshold;
            if(avgsValues1[i]>threshold || avgsValues2[i]>threshold)
               numOfRelevantCells++;
            
            avgsDiffs.push_back( abs(val1 - val2) );
            totalAbsAvgs += (abs(val1) + abs(val2));
            totalAbsDiff +=  abs(val1 - val2);
            totalAbsLog10Diff += abs(log10(val1) - log10(val2));
            
         }
      }
      // Compare the avgs values of nonidentical blocks:
      array<double, velocityCellsPerBlock> zeroAvgs;
      for( uint i = 0; i < velocityCellsPerBlock; ++i ) {
         zeroAvgs[i] = 0;
      }
      for( vector<uint32_t>::const_iterator it = nonIdenticalBlockIds.begin(); it != nonIdenticalBlockIds.end(); ++it ) {
         // Get the block id
         const uint32_t blockId = *it;
         // Get avgs values: 

         const array<double, velocityCellsPerBlock> * avgsValues1;
         const array<double, velocityCellsPerBlock> * avgsValues2;

         unordered_map<uint32_t, array<double, velocityCellsPerBlock> >::const_iterator it2 = avgs1.find( blockId );
         if( it2 == avgs1.end() ) {
            avgsValues1 = &zeroAvgs;
         } else {
            avgsValues1 = &(it2->second);
         }

         it2 = avgs2.find( blockId );
         if( it2 == avgs2.end() ) {
            avgsValues2 = &zeroAvgs;
         } else {
            avgsValues2 = &(it2->second);
         }
         // Get the diff:
         for( uint i = 0; i < velocityCellsPerBlock; ++i ) {
            double val1=avgsValues1->operator[](i)>threshold?avgsValues1->operator[](i):threshold;
            double val2=avgsValues2->operator[](i)>threshold?avgsValues2->operator[](i):threshold;
            if( avgsValues1->operator[](i)>threshold || avgsValues2->operator[](i)>threshold)
               numOfRelevantCells++;
            
            avgsDiffs.push_back( abs(val1 - val2) );
            totalAbsAvgs += (abs(val1) + abs(val2));
            totalAbsDiff +=  abs(val1 - val2);
            totalAbsLog10Diff += abs(log10(val1) - log10(val2));
         }
      }
      numOfIdenticalBlocks += identicalBlockIds.size();
      numOfNonIdenticalBlocks += nonIdenticalBlockIds.size();
   }
   // Get the max and min diff, and the sum of the diff
   double maxDiff = 0;
   double minDiff = numeric_limits<Real>::max();
   double sumDiff = 0;
   for( vector<double>::const_iterator it = avgsDiffs.begin(); it != avgsDiffs.end(); ++it ) {
      sumDiff += *it;
      if( maxDiff < *it ) {
         maxDiff = *it;
      }
      if( minDiff > *it ) {
         minDiff = *it;
      }
   }
   
   const double relativeSumDiff = sumDiff / totalAbsAvgs;
   cout << "File names: " << fileName1 << " & " << fileName2 << endl <<
      "NonIdenticalBlocks:      " << numOfNonIdenticalBlocks << endl <<
      "IdenticalBlocks:         " << numOfIdenticalBlocks <<  endl <<
      "Absolute_Error:          " << totalAbsDiff  << endl <<
      "Mean-Absolute-Error:     " << totalAbsDiff / numOfRelevantCells << endl <<
      "Max-Absolute-Error:      " << maxDiff << endl <<
      "Absolute-log-Error:      " << totalAbsLog10Diff << endl <<
      "Mean-Absolute-log-Error: " << totalAbsLog10Diff / numOfRelevantCells << endl;
   

   return true;
}

/*! Read in the contents of the variable component in both files passed in strings fileName1 and fileName2, and compute statistics and distances as wished
 * \param fileName1 String argument giving the location of the first file to process
 * \param fileName2 String argument giving the location of the second file to process
 * \param varToExtract Pointer to the char array containing the name of the variable to extract
 * \param compToExtract Unsigned int designating the component to extract (0 for scalars)
 * \param verboseOutput Boolean parameter telling whether the output will be verbose or compact
 * \sa convertSILO singleStatistics outputStats pDistance outputDistance printNonVerboseData
 */
bool process2Files(const string fileName1,
                   const string fileName2,
                   const char * varToExtract,
                   const uint compToExtract,
                   const bool verboseOutput,
                   const uint compToExtract2 = 0
                  ) {
   //Check whether the file(s) use new or old vlsv library:
   const bool file1UsesNewVlsvLib = (checkVersion( fileName1 ) == 1.00);
   const bool file2UsesNewVlsvLib = (checkVersion( fileName2 ) == 1.00);

   map<uint, Real> orderedData1;
   map<uint, Real> orderedData2;
   Real absolute, relative, mini, maxi, size, avg, stdev;

   // If the user wants to check avgs, call the avgs check function and return it. Otherwise move on to compare variables:
   if( strcmp(varToExtract, "avgs") == 0 ) {
      vector<uint64_t> cellIds1;
      vector<uint64_t> cellIds2;
      cellIds1.reserve(1);
      cellIds2.reserve(1);
      cellIds1.push_back(compToExtract);
      cellIds2.push_back(compToExtract2);
      // Compare files:
      if( file1UsesNewVlsvLib && file2UsesNewVlsvLib ) {
         if( compareAvgs<newVlsv::Reader, newVlsv::Reader>(fileName1, fileName2, verboseOutput, cellIds1, cellIds2) == false ) { return false; }
      } else if( file1UsesNewVlsvLib && !file2UsesNewVlsvLib ) {
         if( compareAvgs<newVlsv::Reader, oldVlsv::Reader>(fileName1, fileName2, verboseOutput, cellIds1, cellIds2) == false ) { return false; }
      } else if( !file1UsesNewVlsvLib && file2UsesNewVlsvLib ) {
         if( compareAvgs<oldVlsv::Reader, newVlsv::Reader>(fileName1, fileName2, verboseOutput, cellIds1, cellIds2) == false ) { return false; }
      } else {
         // Both are old vlsv format
         if( compareAvgs<oldVlsv::Reader, oldVlsv::Reader>(fileName1, fileName2, verboseOutput, cellIds1, cellIds2) == false ) { return false; }
      }
   } else {
   
      unordered_map<size_t,size_t> cellOrder;

      bool success = true;
      if( file1UsesNewVlsvLib ) {
         success = convertSILO<newVlsv::Reader>(fileName1, varToExtract, compToExtract, &orderedData1, cellOrder, true);
      } else {
         //success = convertSILO<oldVlsv::Reader>(fileName1, varToExtract, compToExtract, &orderedData1, cellOrder, true);
      }
      if( success == false ) {
         cerr << "ERROR Data import error with " << fileName1 << endl;
         return 1;
      }
   
      if( file2UsesNewVlsvLib ) {
         success = convertSILO<newVlsv::Reader>(fileName2, varToExtract, compToExtract, &orderedData2, cellOrder, false);
      } else {
         //success = convertSILO<oldVlsv::Reader>(fileName2, varToExtract, compToExtract, &orderedData2, cellOrder, false);
      }
      if( success == false ) {
         cerr << "ERROR Data import error with " << fileName2 << endl;
         return 1;
      }

      // Basic consistency check
      if(orderedData1.size() != orderedData2.size()) {
         cerr << "ERROR Datasets have different size." << endl;
         return 1;
      }

      // Open VLSV file where the diffence in the chosen variable is written
      const string prefix = fileName1.substr(0,fileName1.find_last_of('.'));
      const string suffix = fileName1.substr(fileName1.find_last_of('.'),fileName1.size());
      const string outputFileName = prefix + ".diff." + varToExtract + suffix;

      vlsv::Writer outputFile;
      if (outputFile.open(outputFileName,MPI_COMM_SELF,0) == false) {
	 cerr << "ERROR failed to open output file '" << outputFileName << "' in " << __FILE__ << ":" << __LINE__ << endl;
	 return false;
      }

      // Clone mesh from input file to diff file
      if (cloneMesh(fileName1,outputFile,"SpatialGrid") == false) return false;
      const string varName = varToExtract;

      singleStatistics(&orderedData1, &size, &mini, &maxi, &avg, &stdev); //CONTINUE
      outputStats(&size, &mini, &maxi, &avg, &stdev, verboseOutput, false);

      singleStatistics(&orderedData2, &size, &mini, &maxi, &avg, &stdev);
      outputStats(&size, &mini, &maxi, &avg, &stdev, verboseOutput, false);

      pDistance(orderedData1, orderedData2, 0, &absolute, &relative, false, cellOrder,outputFile,"SpatialGrid","d0_"+varName);
      outputDistance(0, &absolute, &relative, false, verboseOutput, false);
      pDistance(orderedData1, orderedData2, 0, &absolute, &relative, true, cellOrder,outputFile,"SpatialGrid","d0_sft_"+varName);
      outputDistance(0, &absolute, &relative, true, verboseOutput, false);

      pDistance(orderedData1, orderedData2, 1, &absolute, &relative, false, cellOrder,outputFile,"SpatialGrid","d1_"+varName);
      outputDistance(1, &absolute, &relative, false, verboseOutput, false);
      pDistance(orderedData1, orderedData2, 1, &absolute, &relative, true, cellOrder,outputFile,"SpatialGrid","d1_sft_"+varName);
      outputDistance(1, &absolute, &relative, true, verboseOutput, false);

      pDistance(orderedData1, orderedData2, 2, &absolute, &relative, false, cellOrder,outputFile,"SpatialGrid","d2_"+varName);
      outputDistance(2, &absolute, &relative, false, verboseOutput, false);
      pDistance(orderedData1, orderedData2, 2, &absolute, &relative, true, cellOrder,outputFile,"SpatialGrid","d2_sft_"+varName);
      outputDistance(2, &absolute, &relative, true, verboseOutput, false);

      outputFile.close();
   }
   
   if(verboseOutput == false)
   {
      printNonVerboseData();
      cout << endl;
   }
   
   return 0;
}

/*! Creates the list of grid*.vlsv files present in the folder passed
 * \param dir DIR type pointer to the directory entry to process
 * \param fileList Pointer to a set of strings, return argument for the produced file list
 */
bool processDirectory(DIR* dir, set<string> * fileList) {
   int filesFound = 0, entryCounter = 0;
   
   //const string mask = "grid";
   const string mask = "distributions";
   const string suffix = ".vlsv";
   
   struct dirent* entry = readdir(dir);
   while (entry != NULL) {
      const string entryName = entry->d_name;
      if (entryName.find(mask) == string::npos || entryName.find(suffix) == string::npos) {
         entry = readdir(dir);
         continue;
      }
      fileList->insert(entryName);
      filesFound++;
      entry = readdir(dir);
   }
   if (filesFound == 0) cout << "INFO no matches found" << endl;
   
   return 0;
}

/*! Main function, detects which calling pattern is used and sends to the corresponding processing functions.
 * 
 * \sa process2Files processDirectory
 */
int main(int argn,char* args[])
{
   MPI_Init(&argn,&args);
   
   if (argn < 5)
   {
      cout << endl;
      cout << "USAGE 1: ./vlsvdiff <file1> <file2> <Variable> <component>" << endl;
      cout << "Gives single-file statistics and distances between the two files given, for the variable and component given" << endl;
      cout << "USAGE 2: ./vlsvdiff <folder1> <folder2> <Variable> <component>" << endl;
      cout << "Gives single-file statistics and distances between pairs of files grid*.vlsv taken in alphanumeric order in the two folders given, for the variable and component given" << endl;
      cout << "USAGE 3: ./vlsvdiff <file1> <folder2> <Variable> <component>" << endl;
      cout << "         ./vlsvdiff <folder1> <file2> <Variable> <component>" << endl;
      cout << "Gives single-file statistics and distances between a file, and files grid*.vlsv taken in alphanumeric order in the given folder, for the variable and component given" << endl;
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
   // 5h arg if there is one:
   uint compToExtract2;
   if( argn > 5 ) {
      compToExtract2 = atoi(args[5]);
   } else {
      compToExtract2 = compToExtract;
   }
   
   DIR* dir1 = opendir(fileName1.c_str());
   DIR* dir2 = opendir(fileName2.c_str());
   
   if (dir1 == NULL && dir2 == NULL) {
      cout << "INFO Reading in two files." << endl;
      
      // Process two files with verbose output (last argument true)
      process2Files(fileName1, fileName2, varToExtract, compToExtract, true, compToExtract2);
      //CONTINUE
      
      closedir(dir1);
      closedir(dir2);
   }
   else if (dir1 == NULL || dir2 == NULL)
   {
      // Mixed file and directory
      cout << "#INFO Reading in one file and one directory." << endl;
      set<string> fileList;
      set<string>::iterator it;

      if(dir1 == NULL){
         //file in 1, directory in 2
         processDirectory(dir2, &fileList);
         for(it = fileList.begin(); it != fileList.end();++it){
            // Process two files with non-verbose output (last argument false), give full path to the file processor
            process2Files(fileName1,fileName2 + "/" + *it, varToExtract, compToExtract, false, compToExtract2);
         }
      }

      if(dir2 == NULL){
         //directory in 1, file in 2
         processDirectory(dir1, &fileList);
         for(it = fileList.begin(); it != fileList.end();++it){
            // Process two files with non-verbose output (last argument false), give full path to the file processor
            process2Files(fileName1+"/"+*it,fileName2, varToExtract, compToExtract, false, compToExtract2);
         }
      }

      closedir(dir1);
      closedir(dir2);
      return 1;
   }
   else if (dir1 != NULL && dir2 != NULL)
   {
      // Process two folders, files of the same rank compared, first folder is reference in relative distances
      cout << "#INFO Reading in two directories." << endl;
      set<string> fileList1, fileList2;
      
      // Produce a sorted file list
      processDirectory(dir1, &fileList1);
      processDirectory(dir2, &fileList2);
      
      // Basic consistency check
      if(fileList1.size() != fileList2.size())
      {
         cerr << "ERROR Folders have different number of files." << endl;
         return 1;
      }
      
      set<string>::iterator it1, it2;
      for(it1 = fileList1.begin(), it2 = fileList2.begin();
          it1 != fileList2.end(), it2 != fileList2.end();
          it1++, it2++)
      {
      // Process two files with non-verbose output (last argument false), give full path to the file processor
      process2Files(fileName1 + "/" + *it1,
                    fileName2 + "/" + *it2, varToExtract, compToExtract, false, compToExtract2);
      }
      
      closedir(dir1);
      closedir(dir2);
   }

   MPI_Finalize();
   return 0;
}
