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

/*! \file vlsvdiff.cpp
 * \brief Tool to compare VLSV files.
 * 
 * Tool to compare two VLSV files, two folders with the same number of VLSV files or a folder to a reference file.
 * The tool assumes the files have a name syntax 'grid\.[0-9]+\.vlsv'. It takes four arguments.
 * 
 * Calling patterns are:
 * 
 * "$ vlsvdiff <file1> <file2> <Variable> <component>": Gives single-file statistics and distances between the two files given, for the variable and component given
 
 * "$ vlsvdiff --diff --meshname=<Meshname> <file1> <file2> <Variable> <component>": Gives single-file statistics and distances between the two files given, for the variable and component given
 * 
 * "$ vlsvdiff <folder1> <folder2> <Variable> <component>": Gives single-file statistics and distances between pairs of files grid*.vlsv taken in alphanumeric order in the two folders given, for the variable and component given
 * 
 * "$ vlsvdiff <file1> <folder2> <Variable> <component>" or "$ vlsvdiff <folder1> <file2> <Variable> <component>": Gives single-file statistics and distances between a file, and files grid*.vlsv taken in alphanumeric order in the given folder, for the variable and component given
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits> // YK
#include <list>
#include <set>
#include <sstream>
#include <stdint.h>
#include <string>
#include <typeinfo>

#include "definitions.h"
#include <vlsv_reader.h>
#include "vlsvreaderinterface.h"
#include <vlsv_writer.h>

// #include "../ioread.h" //getFsGridDomainDecomposition
#include <fsgrid.hpp> // computeDomainDecomposition

using namespace std;
using namespace vlsv;

// Command line option,value pairs are parsed and stored to map attributes.
// The key is the option name, and the value is the value. For example, 
// "vlsvdiff --meshname=plaa" would cause 'attributes["meshname"]' to be 
// equal to 'plaa'.
static map<string,string> attributes;

//Global enum and variable
static int gridName; 
enum gridType{
   SpatialGrid,
   fsgrid,
   ionosphere
};


static uint64_t convUInt(const char* ptr, const vlsv::datatype::type& dataType, const uint64_t& dataSize) {
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
 * @param optional If true, this parameter is OK to be missing.
 * @return If true, the array was copied successfully.*/
bool copyArray(vlsv::Reader& input,vlsv::Writer& output,
               const std::string& tagName,
               const list<pair<string,string> >& inputAttribs,
               bool optional=false) {
   bool success = true;

   // Read input array attributes
   map<string,string> outputAttribs;
   if (input.getArrayAttributes(tagName,inputAttribs,outputAttribs) == false) {

      if(!optional) {
         cerr << "ERROR: Failed to read array '" << tagName << "' attributes in " << __FILE__ << ":" << __LINE__ << endl;
         cerr << "Input attributes are:" << endl;
         for (list<pair<string,string> >::const_iterator it=inputAttribs.begin(); it!=inputAttribs.end(); ++it) {
            cerr << "\t '" << it->first << "' = '" << it->second << "'" << endl;
         }
         return false;
      } else {
         // This was an optional parameter, so whatever.
         return true;
      }
   }

   // Figure out arraysize, vectorsize, datasize, and datatype of the copied array
   map<string,string>::const_iterator it;
   map<string,string>::iterator it2;
   it = outputAttribs.find("arraysize"); if (it == outputAttribs.end()) return false;
   uint64_t arraysize = atol(it->second.c_str());
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


/* Small function that overrides how fsgrid diff files are written*/
bool HandleFsGrid(const string& inputFileName,
                  vlsv::Writer& output,
                  std::map<uint, Real> orderedData)
{
   

   //Open input file
  vlsv::Reader input;
   if (input.open(inputFileName) == false) {
      cerr << "ERROR failed to open input file '" << inputFileName << "' in " << __FILE__ << ":" << __LINE__ << endl;
      return false;
   }

   //Read Mesh Attributes
   std::string tagName="MESH";
   list<pair<string,string> > inputAttribs;
   inputAttribs.push_back(make_pair("name","fsgrid"));
   map<string,string> outputAttribs;

   if (input.getArrayAttributes(tagName,inputAttribs,outputAttribs) == false) {
      cerr << "ERROR: Failed to read array '" << tagName << "' attributes in " << __FILE__ << ":" << __LINE__ << endl;
      cerr << "Input attributes are:" << endl;
      for (list<pair<string,string> >::const_iterator it=inputAttribs.begin(); it!=inputAttribs.end(); ++it) {
         cerr << "\t '" << it->first << "' = '" << it->second << "'" << endl;
      }
      return false;
   }

   //Collect needed attributes to a map named patch
   map<string,string>::const_iterator it;
   it = outputAttribs.find("arraysize"); if (it == outputAttribs.end()) return false;
   uint64_t arraysize = atol(it->second.c_str());
   it = outputAttribs.find("vectorsize"); if (it == outputAttribs.end()) return false;
   const uint64_t vectorsize = atol(it->second.c_str());
   it = outputAttribs.find("datasize"); if (it == outputAttribs.end()) return false;
   const uint64_t datasize = atol(it->second.c_str());
   it = outputAttribs.find("datatype"); if (it == outputAttribs.end()) return false;
   const string datatype = it->second;
   it = outputAttribs.find("xperiodic"); if (it == outputAttribs.end()) return false;
   const string xperiodic = it->second;
   it = outputAttribs.find("yperiodic"); if (it == outputAttribs.end()) return false;
   const string yperiodic = it->second;
   it = outputAttribs.find("zperiodic"); if (it == outputAttribs.end()) return false;
   const string zperiodic = it->second;
   it = outputAttribs.find("type"); if (it == outputAttribs.end()) return false;
   const string type = it->second;

   map<string,string>patch;
   patch["arraysize"]=std::to_string(arraysize);
   patch["datasize"]=std::to_string(datasize);
   patch["datatype"]=datatype;
   patch["name"]="fsgrid";
   patch["type"]=type;
   patch["vectorsize"]=std::to_string(vectorsize);
   patch["xperiodic"]=xperiodic;
   patch["yperiodic"]=yperiodic;
   patch["zperiodic"]=zperiodic;


   //Get the global IDs in a vector
   std::vector<uint64_t> globalIds;
   for (const auto iter : orderedData){
      globalIds.push_back( iter.first   );
   }
   
   //Write to file
   output.writeArray("MESH",patch,arraysize,1,&globalIds[0]);

   std::array<int,1> numWritingRanks = {1};
   output.writeParameter("numWritingRanks", &numWritingRanks[0]);

   // Save the FSgrid decomposition
   std::map<std::string, std::string> xmlAttributes;
   const std::string meshName="fsgrid";
   xmlAttributes["mesh"] = meshName;
   std::array<FsGridTools::Task_t, 3> decom = {1,1,1};
   output.writeArray("MESH_DECOMPOSITION", outputAttribs, 3u, 1u, &decom[0]);
   
   //Now for MESH_DOMAIN_SIZES
   inputAttribs.clear();
   inputAttribs.push_back(make_pair("mesh","fsgrid"));
   tagName="MESH_DOMAIN_SIZES";

   if (input.getArrayAttributes(tagName,inputAttribs,outputAttribs) == false) {
      cerr << "ERROR: Failed to read array '" << tagName << "' attributes in " << __FILE__ << ":" << __LINE__ << endl;
      cerr << "Input attributes are:" << endl;
      for (list<pair<string,string> >::const_iterator it=inputAttribs.begin(); it!=inputAttribs.end(); ++it) {
         cerr << "\t '" << it->first << "' = '" << it->second << "'" << endl;
      }
      return false;
   }
   
   
   //Read some attributes we need and parse to our map
   it = outputAttribs.find("datasize"); if (it == outputAttribs.end()) return false;
   const uint64_t  datasize2 = atol(it->second.c_str());
   it = outputAttribs.find("datatype"); if (it == outputAttribs.end()) return false;
   const string  datatype2 = it->second;
   it = outputAttribs.find("vectorsize"); if (it == outputAttribs.end()) return false;
   const uint64_t vectorsize2 = atol(it->second.c_str());
   
   patch.clear();
   patch["arraysize"]="1";
   patch["datasize"]=to_string(datasize2);
   patch["datatype"]=datatype2;
   patch["mesh"]="fsgrid";
   patch["vectorsize"]=to_string(vectorsize2);
   
   //Override MESH_DOMAIN_SIZES
   std::array<uint64_t,2> meshDomainSize({globalIds.size(), 0});
   output.writeArray("MESH_DOMAIN_SIZES",patch ,1,vectorsize2, &meshDomainSize[0]);


   //Close the file
   input.close();


   return true;


}

bool getFsgridDecomposition(vlsvinterface::Reader& file, std::array<int,3>& decomposition){
   uint64_t arraySize;
   uint64_t vectorSize;
   vlsv::datatype::type dataType;
   uint64_t byteSize;
   
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("mesh","fsgrid"));


   std::array<FsGridTools::Task_t,3> fsGridDecomposition={0,0,0}; 
   int* ptr = &fsGridDecomposition[0];

   // Check if array exists:
   bool success = file.getArrayInfo("MESH_DECOMPOSITION",attribs,arraySize,vectorSize,dataType,byteSize);
   if (success == false) {
      // std::cout << "Could not read MESH_DECOMPOSITION" << endl;
      // std::cerr << "ptr " << fsGridDecomposition[0] <<" "<<  fsGridDecomposition[1] << " " <<  fsGridDecomposition[2]<<"\n";
      // std::cerr << "No decomposition found in restart file. Computing fsgrid decomposition for ioread, check results!" <<std::endl;

      int fsgridInputRanks=0;
      if(file.readParameter("numWritingRanks",fsgridInputRanks) == false) {
         std::cerr << "FSGrid writing rank number not found in restart file" << endl;
         exit(1);
      }
      std::array<FsGridTools::FsSize_t,3> gridSize;
      FsGridTools::FsSize_t* gridSizePtr = &gridSize[0];
      success = file.read("MESH_BBOX",attribs, 0, 3, gridSizePtr, false);
      if(success == false){
         std::cerr << "Could not read MESH_BBOX from file" << endl;
         exit(1);
      }
      int64_t* domainInfo = NULL;
      success = file.read("MESH_DOMAIN_SIZES",attribs, 0, fsgridInputRanks, domainInfo);
      if(success == false){
         std::cerr << "Could not read MESH_DOMAIN_SIZES from file" << endl;
         exit(1);
      }
      std::vector<uint64_t> mesh_domain_sizes;
      for (int i = 0; i < 2*fsgridInputRanks; i+=2){
         mesh_domain_sizes.push_back(domainInfo[i]);
      }
      list<pair<string,string> > mesh_attribs;
      mesh_attribs.push_back(make_pair("name","fsgrid"));
      std::vector<FsGridTools::FsSize_t> rank_first_ids(fsgridInputRanks);
      FsGridTools::FsSize_t* ids_ptr = &rank_first_ids[0];

      std::set<FsGridTools::FsIndex_t> x_corners, y_corners, z_corners;
      
      int64_t begin_rank = 0;
      int i = 0;
      for(auto rank_size : mesh_domain_sizes){
         if(file.read("MESH", mesh_attribs, begin_rank, 1, ids_ptr, false) == false){
            std::cerr << "Reading MESH failed.\n";
            exit(1);
         }
         std::array<FsGridTools::FsIndex_t,3> inds = FsGridTools::globalIDtoCellCoord(*ids_ptr, gridSize);
         x_corners.insert(inds[0]);
         y_corners.insert(inds[1]);
         z_corners.insert(inds[2]);
         ++ids_ptr;
         begin_rank += rank_size;
      }

      decomposition[0] = x_corners.size();
      decomposition[1] = y_corners.size();
      decomposition[2] = z_corners.size();
      std::cout << "Fsgrid decomposition computed from MESH to be " << decomposition[0] << " " << decomposition[1] << " " <<decomposition[2] << endl;

      return true;   
   } else {
      // data exists, now read it
      success = file.read("MESH_DECOMPOSITION",attribs, 0, 3, ptr, false);
      decomposition[0] = fsGridDecomposition[0];
      decomposition[1] = fsGridDecomposition[1];
      decomposition[2] = fsGridDecomposition[2];
      std::cout << "Fsgrid decomposition read as " << decomposition[0] << " " << decomposition[1] << " " <<decomposition[2] << endl;
      return true;
   }

   return false;
}

/** Copy the spatial mesh from input to output.
 * @param inputFileName Name of the input file where the mesh is copied from.
 * @param output VLSV reader for the file where the cloned mesh is written.
 * @param meshName Name of the mesh.
 * @return If true, the mesh was successfully cloned.*/
bool cloneMesh(const string& inputFileName,vlsv::Writer& output,const string& meshName, std::map<uint, Real> orderedData) {
   bool success = true;
            
   vlsv::Reader input;
   if (input.open(inputFileName) == false) {
      cerr << "ERROR failed to open input file '" << inputFileName << "' in " << __FILE__ << ":" << __LINE__ << endl;
      return false;
   }
   
   list<pair<string,string> > inputAttribs;
   inputAttribs.push_back(make_pair("name",meshName));
   inputAttribs.clear();
   inputAttribs.push_back(make_pair("mesh",meshName));
   if (copyArray(input,output,"MESH_BBOX",inputAttribs) == false) success = false;

   // Mesh have either individual coordinate arrays (for cartesian geometries)...
   if (copyArray(input,output,"MESH_NODE_CRDS_X",inputAttribs, meshName == "ionosphere") == false) success = false;
   if (copyArray(input,output,"MESH_NODE_CRDS_Y",inputAttribs, meshName == "ionosphere") == false) success = false;
   if (copyArray(input,output,"MESH_NODE_CRDS_Z",inputAttribs, meshName == "ionosphere") == false) success = false;
   
   // Or they have per-node coordinate arrays (for unstructured meshes)
   if (copyArray(input,output,"MESH_NODE_CRDS",inputAttribs, meshName != "ionosphere") == false) success = false;
   if (copyArray(input,output,"MESH_OFFSETS",inputAttribs, meshName != "ionosphere") == false) success = false;

   if (copyArray(input,output,"MESH_GHOST_LOCALIDS",inputAttribs, meshName == "ionosphere") == false) success = false;
   if (copyArray(input,output,"MESH_GHOST_DOMAINS",inputAttribs, meshName == "ionosphere") == false) success = false;
   
   //Only do this if we diff SpatialGrid data
   if (gridName==gridType::SpatialGrid || gridName==gridType::ionosphere){
      if (copyArray(input,output,"MESH_DOMAIN_SIZES",inputAttribs) == false) success = false;

      inputAttribs.clear();
      inputAttribs.push_back(make_pair("name",meshName));
      if (copyArray(input,output,"MESH",inputAttribs) == false) success = false;
   }else{
      HandleFsGrid(inputFileName,output,orderedData);
   }

   input.close();
   return success;
}

/*! Extracts the dataset from the VLSV file opened by convertSILO.
 * \param vlsvReader vlsvinterface::Reader class object used to access the VLSV file
 * \param meshName Address of the string containing the name of the mesh to be extracted
 * \param varToExtract Pointer to the char array containing the name of the variable to extract
 * \param compToExtract Unsigned int designating the component to extract (0 for scalars)
 * \param orderedData Pointer to the return argument map which will get the extracted dataset
 */
bool convertMesh(vlsvinterface::Reader& vlsvReader,
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

   list<pair<string, string> > variableAttributes;
   const string _varToExtract( varToExtract );
   variableAttributes.push_back( make_pair("mesh", meshName) );
   variableAttributes.push_back( make_pair("name", _varToExtract) );
   //Read in array size, vector size, data type and data size of the array "VARIABLE" in the vlsv file (Needed in reading the array)
   if (vlsvReader.getArrayInfo("VARIABLE", variableAttributes, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false) {
      cerr << "ERROR, failed to get array info for '" << _varToExtract << "' on mesh '" << meshName << "' at " << __FILE__ << " " << __LINE__ << endl;
      return false;
   }

   switch(gridName) {
      case gridType::SpatialGrid:
         {
            std::vector<char> variableBuffer(variableVectorSize * variableDataSize);
            float *variablePtrFloat = reinterpret_cast<float *>(variableBuffer.data());
            double *variablePtrDouble = reinterpret_cast<double *>(variableBuffer.data());
            uint *variablePtrUint = reinterpret_cast<uint *>(variableBuffer.data());
            int *variablePtrInt = reinterpret_cast<int *>(variableBuffer.data());

            // Read the mesh array one node (of a spatial cell) at a time
            // and create a map which contains each cell's CellID and variable to be extracted
            //Get local cell ids:
            vector<uint64_t> local_cells;
            if ( vlsvReader.getCellIds( local_cells, meshName) == false ) {
               cerr << "Failed to read cell ids at "  << __FILE__ << " " << __LINE__ << endl;
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

            if (storeCellOrder == true) {
               cellOrder.clear();
            }

            orderedData->clear();

            for (uint64_t i=0; i<local_cells.size(); ++i) {
               const short int amountToReadIn = 1;
               const uint64_t & startingReadIndex = i;
               if (vlsvReader.readArray("VARIABLE", variableAttributes, startingReadIndex, amountToReadIn, variableBuffer.data()) == false) {
                  cerr << "ERROR, failed to read variable '" << _varToExtract << "' at " << __FILE__ << " " << __LINE__ << endl;
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
         }
         break;
 
      case gridType::fsgrid:

         {
            // Get Spatial Grid's  max refinement Level
            int maxRefLevel=0;
            list<pair<string, string>> meshAttributesIn;
            meshAttributesIn.push_back(make_pair("name", "SpatialGrid"));
            map<string,string> meshAttributesOut;
            if (vlsvReader.getArrayAttributes("MESH", meshAttributesIn,meshAttributesOut) == false)
            {
               cerr << "ERROR, failed to get array info for '" << _varToExtract << "' at " << __FILE__ << " " << __LINE__ << endl;
               return false;
            }

            std::map<string, string>::iterator attributesOutIt;
            attributesOutIt = meshAttributesOut.find("max_refinement_level");
            if (attributesOutIt != meshAttributesOut.end())
            {
               maxRefLevel = stoi(attributesOutIt->second);
            }
            int numtasks;
            int xcells,ycells,zcells;
            vlsvReader.readParameter("numWritingRanks",numtasks);
            vlsvReader.readParameter("xcells_ini",xcells);
            vlsvReader.readParameter("ycells_ini",ycells);
            vlsvReader.readParameter("zcells_ini",zcells);
            xcells*=pow(2,maxRefLevel);
            ycells*=pow(2,maxRefLevel);
            zcells*=pow(2,maxRefLevel);
            std::array<int,3> GlobalBox={xcells,ycells,zcells};
            std::array<int,3> thisDomainDecomp;

            //Compute Domain Decomposition Scheme for this vlsv file
            //FsGridTools::computeDomainDecomposition(GlobalBox,numtasks,thisDomainDecomp);
            getFsgridDecomposition(vlsvReader, thisDomainDecomp);


            std::array<int32_t,3> taskSize,taskStart;
            std::array<int32_t,3> taskEnd;
            int readOffset=0;
            int index,my_x,my_y,my_z;
            orderedData->clear();

            for (int task=0; task<numtasks; task++){

               my_x=task/thisDomainDecomp[2]/thisDomainDecomp[1];
               my_y=(task/thisDomainDecomp[2])%thisDomainDecomp[1];
               my_z=task%thisDomainDecomp[2];


               taskStart[0] = FsGridTools::calcLocalStart(GlobalBox[0], thisDomainDecomp[0], my_x);
               taskStart[1] = FsGridTools::calcLocalStart(GlobalBox[1], thisDomainDecomp[1], my_y);
               taskStart[2] = FsGridTools::calcLocalStart(GlobalBox[2], thisDomainDecomp[2], my_z);

               taskSize[0] = FsGridTools::calcLocalSize(GlobalBox[0], thisDomainDecomp[0], my_x);
               taskSize[1] = FsGridTools::calcLocalSize(GlobalBox[1], thisDomainDecomp[1], my_y);
               taskSize[2] = FsGridTools::calcLocalSize(GlobalBox[2], thisDomainDecomp[2], my_z);

               taskEnd[0]= taskStart[0]+taskSize[0];
               taskEnd[1]= taskStart[1]+taskSize[1];
               taskEnd[2]= taskStart[2]+taskSize[2];

               int64_t readSize=  taskSize[0] * taskSize[1] * taskSize[2] ;
               //Allocate vector for reading
               std::vector<Real> buffer(readSize*variableVectorSize);

               if ( variableDataSize==sizeof(Real)){
                  if (vlsvReader.readArray("VARIABLE", variableAttributes, readOffset, readSize,  (char*)buffer.data()) == false) {
                     cerr << "ERROR, failed to read variable '" << _varToExtract << "' at " << __FILE__ << " " << __LINE__ << endl;
                     variableSuccess = false; 
                     break;
                  }
               }else{
                  std::vector<float> tmpbuffer(readSize * variableVectorSize);
                  if (vlsvReader.readArray("VARIABLE", variableAttributes, readOffset, readSize, (char *)tmpbuffer.data()) == false){
                     cerr << "ERROR, failed to read variable '" << _varToExtract << "' at " << __FILE__ << " " << __LINE__ << endl;
                     variableSuccess = false;
                     break;
                  }
                  for (unsigned int i = 0; i < readSize * variableVectorSize; i++){
                     buffer[i] = tmpbuffer[i];
                  }
               }

               uint64_t globalindex,counter=0;;
               for (int z=taskStart[2]; z<taskEnd[2]; z++){
                  for (int y=taskStart[1]; y< taskEnd[1]; y++){
                     for (int x=taskStart[0]; x<taskEnd[0]; x++){
                        globalindex= x + y*xcells + z*xcells*ycells;
                        Real data;
                        switch (variableDataType){
                           case datatype::type::FLOAT:
                              if (variableDataSize == sizeof(float))
                                 memcpy(&data, &buffer[counter + compToExtract], sizeof(float));
                              if (variableDataSize == sizeof(double))
                                 memcpy(&data, &buffer[counter + compToExtract], sizeof(double));
                              break;
                           case datatype::type::UINT:
                              memcpy(&data, &buffer[counter + compToExtract], sizeof(uint));
                              break;
                           case datatype::type::INT:
                              memcpy(&data, &buffer[counter + compToExtract], sizeof(int));
                              break;
                           case datatype::type::UNKNOWN:
                              cerr << "ERROR, BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << endl;
                              break;
                        }
                        //Add to map
                        orderedData->insert(pair<uint64_t, Real>(globalindex, data));
                        counter+=variableVectorSize;
                     }
                  }
               }
               readOffset+=readSize;

            }
         }
         break;

      case gridType::ionosphere:

         if(compToExtract >= variableVectorSize) {
            cerr << "ERROR invalid component, this variable has size " << variableVectorSize << endl;
            abort();
         }
         orderedData->clear();
         
         switch(variableDataType) {
            case datatype::type::FLOAT: 
               {
                  if(variableDataSize == sizeof(double)) { 
                     std::vector<double> buffer(variableVectorSize * variableArraySize);
                     // The mesh is simply one big blob that can be read in one go.
                     if(vlsvReader.readArray("VARIABLE", variableAttributes, 0, variableArraySize, (char*)buffer.data()) == false) {
                        cerr << "ERROR, failed to read variable '" << _varToExtract << "' at " << __FILE__ << " " << __LINE__ << endl;
                        variableSuccess = false; 
                        break;
                     }

                     for(unsigned int i=0; i<variableArraySize; i++) {
                        orderedData->insert(pair<uint64_t, Real>(i, buffer[i*variableVectorSize + compToExtract]));
                     }
                  } else if(variableDataSize == sizeof(float)) {
                     std::vector<double> buffer(variableVectorSize * variableArraySize);
                     // The mesh is simply one big blob that can be read in one go.
                     if(vlsvReader.readArray("VARIABLE", variableAttributes, 0, variableArraySize, (char*)buffer.data()) == false) {
                        cerr << "ERROR, failed to read variable '" << _varToExtract << "' at " << __FILE__ << " " << __LINE__ << endl;
                        variableSuccess = false; 
                        break;
                     }

                     for(unsigned int i=0; i<variableArraySize; i++) {
                        orderedData->insert(pair<uint64_t, Real>(i, buffer[i*variableVectorSize + compToExtract]));
                     }
                  }
               }
               break;
            default:
               cerr << "Error: No support for ionosphere parameters that are not float-valued implemented, at " << __FILE__ << " " << __LINE__ << endl;
               break;
         }

         break;
      default:
         cerr<<"meshName not recognized\t" << __FILE__ << " " << __LINE__ <<endl;
         abort();
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
                 Real& time,
                 const bool& storeCellOrder=false) {
   bool success = true;

   // Open VLSV file for reading:
   T vlsvReader;   

   if (vlsvReader.open(fileName) == false) {
      cerr << "Failed to open '" << fileName << "'" << endl;
      cerr << "VLSV error " << vlsvReader.getErrorString() << endl;
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
      if (*it != attributes["--meshname"]) continue;

      if (convertMesh(vlsvReader, *it, varToExtract, compToExtract, orderedData, cellOrder, storeCellOrder) == false) {
         return false;
      }      
   }

   vlsvReader.readParameter("time", time);

   vlsvReader.close();
   return success;
}

/*! Shift the second file to the average of the first
 * \param orderedData1 Pointer to the reference file's data
 * \param orderedData2 Pointer to the data to be shifted
 * \param shiftedData2 Pointer to where the shifted data of the second file will be put
 */
bool shiftAverage(const map<uint, Real>* const orderedData1,
                  const map<uint, Real>* const orderedData2,
                  map<uint,Real>* shiftedData2
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
   map<uint,Real> shiftedData2;
   map<uint,Real>* data2 = const_cast< map<uint,Real>* >(&orderedData2);   

   if (doShiftAverage == true) {
      shiftAverage(&orderedData1,&orderedData2,&shiftedData2);
      data2 = &shiftedData2;
   }

   // Reset old values
   *absolute = 0.0;
   *relative = 0.0;

   vector<Real> array(orderedData1.size());
   for (size_t i=0; i<array.size(); ++i) array[i] = -1.0;

   Real length = 0.0;
   if (p == 0) {
      for (map<uint,Real>::const_iterator it1=orderedData1.begin(); it1!=orderedData1.end(); ++it1) {
         map<uint,Real>::const_iterator it2 = data2->find(it1->first);
         Real value = 0.0;
         if (it2 != data2->end()) {
            value = abs(it1->second - it2->second);
            *absolute = max(*absolute, value);
            length    = max(length, abs(it1->second));
         
            }
         if (gridName==gridType::SpatialGrid){  
            array[cellOrder.at(it1->first)] = value;
         }else if (gridName==gridType::fsgrid || gridName==gridType::ionosphere) {   
            array.at(it1->first)=value;
         }  
      }
   } else if (p == 1) {
      for (map<uint,Real>::const_iterator it1=orderedData1.begin(); it1!=orderedData1.end(); ++it1) {
         map<uint,Real>::const_iterator it2 = data2->find(it1->first);
         Real value = 0.0;
         if (it2 != data2->end()) {
            value = abs(it1->second - it2->second);
            *absolute += value;
            length    += abs(it1->second);
         
            }
         if (gridName==gridType::SpatialGrid){  
            array[cellOrder.at(it1->first)] = value;
         }else if (gridName==gridType::fsgrid || gridName==gridType::ionosphere) {   
            array[it1->first]=value;
         }  
      }
   } else {
      for (map<uint,Real>::const_iterator it1=orderedData1.begin(); it1!=orderedData1.end(); ++it1) {
         map<uint,Real>::const_iterator it2 = data2->find(it1->first);
         Real value = 0.0;
         if (it2 != data2->end()) {
            value = pow(abs(it1->second - it2->second), p);
            *absolute += value;
            length    += pow(abs(it1->second), p);
         
            }
         if (gridName==gridType::SpatialGrid){  
            array[cellOrder.at(it1->first)] = pow(value,1.0/p);
         }else if (gridName==gridType::fsgrid || gridName==gridType::ionosphere) {   
            array[it1->first]=pow(value,1.0/p);
         }  
      }
      *absolute = pow(*absolute, 1.0 / p);
      length = pow(length, 1.0 / p);
   }

   if (length != 0.0) *relative = *absolute / length;
   else {
      cout << "WARNING (pDistance) : length of reference is 0.0, cannot divide to give relative distance." << endl;
      *relative = -1;
   }

   // Write out the difference (if requested):
   if (attributes.find("--diff") != attributes.end()) {
      map<string,string> attributes;
      attributes["mesh"] = meshName;
      attributes["name"] = varName;
      if(meshName == "ionosphere") {
         attributes["centering"] = "node";
      }


      if (outputFile.writeArray("VARIABLE",attributes,array.size(),1,&(array[0])) == false) {
         cerr << "ERROR failed to write variable '" << varName << "' to output file in " << __FILE__ << ":" << __LINE__ << endl;
         return 1;
      }
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
         cout << "The absolute " << p << "-distance between both datasets is " << setprecision(3) << *absolute  << endl;
         cout << "The relative " << p << "-distance between both datasets is " << setprecision(3) << *relative  << endl;
      } else {
         cout << "The average-shifted absolute " << p << "-distance between both datasets is " << setprecision(3) << *absolute  << endl;
         cout << "The average-shifted relative " << p << "-distance between both datasets is " << setprecision(3) << *relative  << endl;
      }
   } else {
      static vector<Real> fileOutputData;
      static uint fileNumber = 0;
      
      if(lastCall == true) {
         vector<Real>::const_iterator it;
         for(it = fileOutputData.begin(); it != fileOutputData.end(); it++) {
            cout << setprecision(3) << *it << "\t";
         }
         fileOutputData.clear();
         return 0;
      }
      
      fileOutputData.push_back(*absolute);
      fileOutputData.push_back(*relative);
   }
   return 0;
}

/*! In verbose mode print delta t, in non-verbose store them for later output when lastCall is true
 * \param dt delta t
 * \param verboseOutput Boolean parameter telling whether the output is verbose or compact
 * \param lastCall Boolean parameter telling whether this is the last call to the function
 */
bool outputDt(
   const Real dt,
   const bool verboseOutput,
   const bool lastCall
) {
   if(verboseOutput == true) {
      cout << "The delta t between both datasets is " << dt << endl;
   } else {
      static vector<Real> fileOutputData;
      static uint fileNumber = 0;
      
      if(lastCall == true) {
         vector<Real>::const_iterator it;
         for(auto f : fileOutputData) {
            cout << f << "\t";
         }
         fileOutputData.clear();
         return 0;
      }
      
      fileOutputData.push_back(dt);
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
   outputDt(0, false, true);
   
   return 0;
}

bool getBlockIds(vlsvinterface::Reader& vlsvReader,
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
   attribs.push_back(make_pair("mesh", attributes["--meshname"]));

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

// Reads avgs values of some given cell id
// Input:
// [0] vlsvReader -- Some vlsv reader with a file open
// [1] cellId -- The spatial cell's ID
// Output:
// [2] avgs -- Saves the output into an unordered map with block id as the key and an array of avgs as the value
// return false or true depending on whether the operation was successful
template <class T>
bool readAvgs( T & vlsvReader,
               string name,
               const unordered_map<uint64_t, pair<uint64_t, uint32_t>> & cellsWithBlocksLocations,
               const uint64_t & cellId, 
               unordered_map<uint32_t, array<double, 64> > & avgs ) {
   // Get the block ids:
   vector<uint32_t> blockIds;
   if( getBlockIds( vlsvReader, cellsWithBlocksLocations, cellId, blockIds ) == false ) { return false; }
   // Read avgs:
   list<pair<string, string> > attribs;
   attribs.push_back(make_pair("name", name));
   attribs.push_back(make_pair("mesh", attributes["--meshname"]));

   datatype::type dataType;
   uint64_t arraySize, vectorSize, dataSize;
   if (vlsvReader.getArrayInfo("BLOCKVARIABLE", attribs, arraySize, vectorSize, dataType, dataSize) == false) {
      //no 
//      cerr << "ERROR READING BLOCKVARIABLE AT " << __FILE__ << " " << __LINE__ << endl;
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
   const string meshName = attributes["--meshname"];
   vlsv::datatype::type cwb_dataType;
   uint64_t cwb_arraySize, cwb_vectorSize, cwb_dataSize;
   list<pair<string, string> > attribs;

   //Get the mesh name for reading in data from the correct place
   attribs.push_back( make_pair("mesh", meshName) );

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
      if( readAvgs( vlsvReader1, "proton", cellsWithBlocksLocations1, cellId1, avgs1 ) == false ) {
         if( readAvgs( vlsvReader1, "avgs", cellsWithBlocksLocations1, cellId1, avgs1 ) == false ) {
            cerr << "ERROR, FAILED TO READ AVGS AT " << __FILE__ << " " << __LINE__ << endl;
            return false;
         }
      }
      
      if( readAvgs( vlsvReader2, "proton", cellsWithBlocksLocations2, cellId2, avgs2 ) == false ) {
         if( readAvgs( vlsvReader2, "avgs", cellsWithBlocksLocations2, cellId2, avgs2 ) == false ) {
            cerr << "ERROR, FAILED TO READ AVGS AT " << __FILE__ << " " << __LINE__ << endl;
            return false;
         }
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

   Real time1 {0.0};
   Real time2 {0.0};
   vlsvReader1.readParameter("time", time1);
   vlsvReader2.readParameter("time", time2);

   const double relativeSumDiff = sumDiff / totalAbsAvgs;
   cout << "File names: " << fileName1 << " & " << fileName2 << endl <<
      setprecision(3) <<
      "NonIdenticalBlocks:      " << numOfNonIdenticalBlocks << endl <<
      "IdenticalBlocks:         " << numOfIdenticalBlocks <<  endl <<
      "Absolute_Error:          " << totalAbsDiff  << endl <<
      "Mean-Absolute-Error:     " << totalAbsDiff / numOfRelevantCells << endl <<
      "Max-Absolute-Error:      " << maxDiff << endl <<
      "Absolute-log-Error:      " << totalAbsLog10Diff << endl <<
      "Mean-Absolute-log-Error: " << totalAbsLog10Diff / numOfRelevantCells << endl <<
      "Delta-t: " << time2 - time1 << endl;

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
   map<uint, Real> orderedData1;
   map<uint, Real> orderedData2;
   Real absolute, relative, mini, maxi, size, avg, stdev;

   // If the user wants to check avgs, call the avgs check function and return it. Otherwise move on to compare variables:
   if( strcmp(varToExtract, "proton") == 0 && attributes.find("--no-distrib") == attributes.end()) {
      vector<uint64_t> cellIds1;
      vector<uint64_t> cellIds2;
      cellIds1.reserve(1);
      cellIds2.reserve(1);
      cellIds1.push_back(compToExtract);
      cellIds2.push_back(compToExtract2);
      // Compare files:
      if (compareAvgs<vlsvinterface::Reader, vlsvinterface::Reader>(fileName1, fileName2, verboseOutput, cellIds1, cellIds2) == false) { 
         return false; 
      }
   } else {
      unordered_map<size_t,size_t> cellOrder;
   
      bool success = true;
      Real time1 {0.0};
      success = convertSILO<vlsvinterface::Reader>(fileName1, varToExtract, compToExtract, &orderedData1, cellOrder, time1, true);

      if( success == false ) {
         cerr << "ERROR Data import error with " << fileName1 << endl;
         return 1;
      }

      Real time2 {0.0};
      success = convertSILO<vlsvinterface::Reader>(fileName2, varToExtract, compToExtract, &orderedData2, cellOrder, time2, false);

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
      string outputFileName = prefix + ".diff." + varToExtract + suffix;
      const string varName = varToExtract;
      vlsv::Writer outputFile;
      if (attributes.find("--diff") != attributes.end()) {
         if (outputFileName[0] == '.' && outputFileName[1] == '/') {
            outputFileName = outputFileName.substr(2,string::npos);
         }
         
         for (size_t s=0; s<outputFileName.size(); ++s)
           if (outputFileName[s] == '/') outputFileName[s] = '_';

         if (outputFile.open(outputFileName,MPI_COMM_SELF,0) == false) {
            cerr << "ERROR failed to open output file '" << outputFileName << "' in " << __FILE__ << ":" << __LINE__ << endl;
            return false;
         }


         map<string,string>::const_iterator it = attributes.find("--meshname");
         if (cloneMesh(fileName1,outputFile,it->second,orderedData1) == false) {
            std::cerr<<"Failed"<<std::endl;
            return false;
         }
      }

      singleStatistics(&orderedData1, &size, &mini, &maxi, &avg, &stdev); //CONTINUE
      // Clone mesh from input file to diff file
      outputStats(&size, &mini, &maxi, &avg, &stdev, verboseOutput, false);

      singleStatistics(&orderedData2, &size, &mini, &maxi, &avg, &stdev);
      outputStats(&size, &mini, &maxi, &avg, &stdev, verboseOutput, false);

      pDistance(orderedData1, orderedData2, 0, &absolute, &relative, false, cellOrder,outputFile,attributes["--meshname"],"d0_"+varName);
      outputDistance(0, &absolute, &relative, false, verboseOutput, false);
      pDistance(orderedData1, orderedData2, 0, &absolute, &relative, true, cellOrder,outputFile,attributes["--meshname"],"d0_sft_"+varName);
      outputDistance(0, &absolute, &relative, true, verboseOutput, false);

      pDistance(orderedData1, orderedData2, 1, &absolute, &relative, false, cellOrder,outputFile,attributes["--meshname"],"d1_"+varName);
      outputDistance(1, &absolute, &relative, false, verboseOutput, false);
      pDistance(orderedData1, orderedData2, 1, &absolute, &relative, true, cellOrder,outputFile,attributes["--meshname"],"d1_sft_"+varName);
      outputDistance(1, &absolute, &relative, true, verboseOutput, false);

      pDistance(orderedData1, orderedData2, 2, &absolute, &relative, false, cellOrder,outputFile,attributes["--meshname"],"d2_"+varName);
      outputDistance(2, &absolute, &relative, false, verboseOutput, false);
      pDistance(orderedData1, orderedData2, 2, &absolute, &relative, true, cellOrder,outputFile,attributes["--meshname"],"d2_sft_"+varName);
      outputDistance(2, &absolute, &relative, true, verboseOutput, false);

      outputDt(time2 - time1, verboseOutput, false);

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
bool processDirectory(DIR* dir, set<string>* fileList) {
   int filesFound = 0, entryCounter = 0;
   
   const string mask = attributes["--filemask"];
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

void printHelp(const map<string,string>& defAttribs,const map<string,string>& descriptions) {
   cout << endl;
   cout << "VLSVDIFF command line attributes are given as option=value pairs, value can be empty." << endl;
   cout << "If the default value is 'unset', then giving the option in the command line turns it on." << endl;
   cout << "For example, \"vlsvdiff --help\" displays this message and the option '--help' does not have a value." << endl << endl;
   
   cout << "Known attributes and default values are:" << endl;
   for (map<string,string>::const_iterator it=defAttribs.begin(); it!=defAttribs.end(); ++it) {
      cout << endl;
      const size_t optionWidth   = 30;
      const size_t descrMaxWidth = 120;
      
      // Print the option,value pair so that the field width is always 30 characters
      string option = it->first;
      if (it->second.size() > 0) option = option + "=" + it->second;
      else option = option + " (unset)";
      
      if (option.size() < optionWidth) {
         size_t padding = optionWidth-option.size();
         for (size_t i=0; i<padding; ++i) option = option + ' ';
      }
      cout << option;
      
      // Print the description, possibly on multiple lines.
      map<string,string>::const_iterator descr = descriptions.find(it->first);
      if (descr == descriptions.end()) {
         cout << "(no description given)" << endl;
         continue;
      }
      
      // If the description fits in the first line, print it and continue
      if (descr->second.size() <= descrMaxWidth-optionWidth) {
         cout << descr->second << endl;
         continue;
      }

      // Print the description on multiple lines. First parse the description 
      // string and store each word to a vector.
      vector<string> text;
      size_t i=0;
      while (i < descr->second.size()) {
         size_t i_space = descr->second.find_first_of(' ',i);
         if (i_space == string::npos) i_space = descr->second.size();
         text.push_back(descr->second.substr(i,i_space-i));
         i = i_space+1;
      }

      // Write out the words in vector 'text' so that the length of any line 
      // does not exceed descrMaxWidth characters.
      i = optionWidth;
      for (size_t s=0; s<text.size(); ++s) {
         if (i+text[s].size() <= descrMaxWidth) {
            cout << text[s] << ' ';
            i += text[s].size()+1;
         } else {
            cout << endl;
            for (unsigned int j=0; j<optionWidth; ++j) cout << ' ';
            i = optionWidth;
            
            cout << text[s] << ' ';
            i += text[s].size()+1;
         }
      }
      cout << endl;
   }
   cout << endl;
}

/*! Main function, detects which calling pattern is used and sends to the corresponding processing functions.
 * 
 * \sa process2Files processDirectory
 */
int main(int argn,char* args[]) {
   MPI_Init(&argn,&args);

   // Create default attributes
   map<string,string> defAttribs;
   map<string,string> descriptions;
   defAttribs.insert(make_pair("--meshname","SpatialGrid"));
   defAttribs.insert(make_pair("--filemask","bulk"));
   defAttribs.insert(make_pair("--help",""));
   defAttribs.insert(make_pair("--no-distrib",""));
   defAttribs.insert(make_pair("--diff",""));

   descriptions["--meshname"] = "Name of the spatial mesh that is used in diff.";
   descriptions["--filemask"] = "File mask used in directory comparison mode. For example, if you want to compare files starting with 'fullf', set '--filemask=fullf'.";
   descriptions["--help"]     = "Print this help message.";
   descriptions["--diff"]     = "If set, difference file(s) are written.";
   descriptions["--no-distrib"] = "If set, velocity block data are not compared even if the given variable corresponds to velocity block data.";


   // Create default attributes
   for (map<string,string>::const_iterator it=defAttribs.begin(); it!=defAttribs.end(); ++it) {
      if (it->second.size() == 0) continue;
      attributes.insert(make_pair(it->first,it->second));
   }

   vector<string> argsVector;

   // Parse attributes,value pairs from command line
   int i=0;
   while (i < argn) {
      if (args[i][1] == '\0') {argsVector.push_back(args[i]); ++i; continue;}
      if (args[i][0] == '-' && args[i][1] == '-') {
         string s = args[i];
         if (argn > i) {
            if (s.find("=") == string::npos) {
               attributes.insert(make_pair(string(args[i]),""));
            } else {
               size_t pos = s.find("=");
               string arg = s.substr(0,s.find('='));
               string val = s.substr(s.find('=')+1,string::npos);
               attributes[arg] = val;
            }
            ++i;
            continue;
         } else {
            if (s.find("=") == string::npos) {
               attributes.insert(make_pair(string(args[i]),""));
            } else {
               size_t pos = s.find("=");
               string arg = s.substr(0,s.find('='));
               string val = s.substr(s.find('=')+1,string::npos);
               attributes[arg] = val;
            }
            ++i;
            break;
         }
      } else {
         argsVector.push_back(args[i]);
      }
      ++i;
   }

   if (attributes.find("--help") != attributes.end()) {
      printHelp(defAttribs,descriptions);
      return 0;
   }



   if (argsVector.size() < 5) {
      cout << endl;
      cout << "USAGE 1: ./vlsvdiff <file1> <file2> <Variable> <component>" << endl;
      cout << "Gives single-file statistics and distances between the two files given, for the variable and component given" << endl;
      cout << "USAGE 2: ./vlsvdiff <folder1> <folder2> <Variable> <component>" << endl;
      cout << "Gives single-file statistics and distances between pairs of files grid*.vlsv taken in alphanumeric order in the two folders given, for the variable and component given" << endl;
      cout << "USAGE 3: ./vlsvdiff <file1> <folder2> <Variable> <component>" << endl;
      cout << "         ./vlsvdiff <folder1> <file2> <Variable> <component>" << endl;
      cout << "Gives single-file statistics and distances between a file, and files grid*.vlsv taken in alphanumeric order in the given folder, for the variable and component given" << endl;
      cout << endl;
      cout << "Type ./vlsvdiff --help for more info" << endl;
      cout << endl;
      return 1;
   }

   // 1st arg is file1 name
   const string fileName1 = argsVector[1];
   // 2nd arg is file2 name
   const string fileName2 = argsVector[2];
   // 3rd arg is variable name
   const char* varToExtract = argsVector[3].c_str();

   // 4th arg is its component, 0 for scalars, 2 for z component etc
   uint compToExtract = atoi(argsVector[4].c_str());
   // 5h arg if there is one:
   uint compToExtract2;
   if(argsVector.size() > 5 ) {
      compToExtract2 = atoi(argsVector[5].c_str());
   } else {
      compToExtract2 = compToExtract;
   }
   

   //Figure out Meshname
   if (attributes["--meshname"] == "SpatialGrid") { 
      gridName=gridType::SpatialGrid ;
   }else if (attributes["--meshname"]=="fsgrid"){
      gridName=gridType::fsgrid ;
   }else if (attributes["--meshname"]=="ionosphere"){
      gridName=gridType::ionosphere ;
   }else{
      std::cout<<attributes["--meshname"]<<std::endl;
      std::cerr<<"Wrong grid type"<<std::endl;
      abort();
   }




   DIR* dir1 = opendir(fileName1.c_str());
   DIR* dir2 = opendir(fileName2.c_str());

   if (dir1 == nullptr && dir2 == nullptr) {
      cout << "INFO Reading in two files." << endl;
      
      // Process two files with verbose output (last argument true)
      process2Files(fileName1, fileName2, varToExtract, compToExtract, true, compToExtract2);
   } else if (dir1 == nullptr || dir2 == nullptr) {
      // Mixed file and directory
      cout << "#INFO Reading in one file and one directory." << endl;
      set<string> fileList;

      if(dir1 == nullptr){
         //file in 1, directory in 2
         processDirectory(dir2, &fileList);
         for(auto f : fileList) {
            // Process two files with non-verbose output (last argument false), give full path to the file processor
            process2Files(fileName1,fileName2 + "/" + f, varToExtract, compToExtract, false, compToExtract2);
         }
         closedir(dir2);
      }

      if(dir2 == nullptr){
         //directory in 1, file in 2
         processDirectory(dir1, &fileList);
         for(auto f : fileList) {
            // Process two files with non-verbose output (last argument false), give full path to the file processor
            process2Files(fileName1 + "/" + f,fileName2, varToExtract, compToExtract, false, compToExtract2);
         }
         closedir(dir1);
      }
   } else if (dir1 && dir2) {
      // Process two folders, files of the same rank compared, first folder is reference in relative distances
      cout << "#INFO Reading in two directories." << endl;
      set<string> fileList1, fileList2;
      
      // Produce a sorted file list
      processDirectory(dir1, &fileList1);
      processDirectory(dir2, &fileList2);
      
      // Basic consistency check
      if(fileList1.size() != fileList2.size()) {
         cerr << "ERROR Folders have different number of files." << endl;
         return 1;
      }
      
      // TODO zip these once we're using C++23
      for (auto it1 = fileList1.begin(), it2 = fileList2.begin(); it1 != fileList2.end(), it2 != fileList2.end(); it1++, it2++) {
         // Process two files with non-verbose output (last argument false), give full path to the file processor
         process2Files(fileName1 + "/" + *it1, fileName2 + "/" + *it2, varToExtract, compToExtract, false, compToExtract2);
      }
      
      closedir(dir1);
      closedir(dir2);
   }

   MPI_Finalize();
   return 0;
}
