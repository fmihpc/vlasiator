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

#include "definitions.h"
#include "mpiconversion.h"
#include "muxml.h"
#include "vlscommon.h"
#include "vlsvreader2.h"

using namespace std;

VLSVReader::VLSVReader() {
   endiannessReader = detectEndianness();
   fileOpen = false;
   swapIntEndianness = false;
}

VLSVReader::~VLSVReader() {
   filein.close();   
}

bool VLSVReader::close() {
   filein.close();
   xmlReader.clear();
   return true;
}

bool VLSVReader::getArrayInfo(const std::string& tagName,const std::string& arrayName,uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const {   
   if (fileOpen == false) return false;
   list<pair<string,string> > attribs;
   XMLNode* node = xmlReader.find(tagName,attribs);
   if (node == NULL) return false;

   arraySize = atoi(node->attributes["arraysize"].c_str());
   vectorSize = atoi(node->attributes["vectorsize"].c_str());
   byteSize = atoi(node->attributes["datasize"].c_str());
   if (node->attributes["datatype"] == "int") dataType = VLSV::INT;
   else if (node->attributes["datatype"] == "uint") dataType = VLSV::UINT;
   else if (node->attributes["datatype"] == "float") dataType = VLSV::FLOAT;
   else {
      cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
      return false;
   }
   return true;
}

bool VLSVReader::getArrayInfo(const std::string& tagName,const std::string& arrayName,const std::string& meshName,uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const {
   if (fileOpen == false) return false;
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name",arrayName));
   attribs.push_back(make_pair("mesh",meshName));
   XMLNode* node = xmlReader.find(tagName,attribs);
   if (node == NULL) return false;

   arraySize = atoi(node->attributes["arraysize"].c_str());
   vectorSize = atoi(node->attributes["vectorsize"].c_str());
   byteSize = atoi(node->attributes["datasize"].c_str());
   if (node->attributes["datatype"] == "int") dataType = VLSV::INT;
   else if (node->attributes["datatype"] == "uint") dataType = VLSV::UINT;
   else if (node->attributes["datatype"] == "float") dataType = VLSV::FLOAT;
   else {
      cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
      return false;
   }
   return true;
}

bool VLSVReader::getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& dataSize) const {
   if (fileOpen == false) return false;
   XMLNode* node = xmlReader.find(tagName,attribs);
   if (node == NULL) return false;
   
   arraySize = atoi(node->attributes["arraysize"].c_str());
   vectorSize = atoi(node->attributes["vectorsize"].c_str());
   dataSize = atoi(node->attributes["datasize"].c_str());
   if (node->attributes["datatype"] == "int") dataType = VLSV::INT;
   else if (node->attributes["datatype"] == "uint") dataType = VLSV::UINT;
   else if (node->attributes["datatype"] == "float") dataType = VLSV::FLOAT;
   else {
      cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
      return false;
   }
   return true;
}

bool VLSVReader::getMeshNames(list<string>& meshNames) const {
   meshNames.clear();
   if (fileOpen == false) return false;
   
   XMLNode* node = xmlReader.find("VLSV");
   for (multimap<string,XMLNode*>::const_iterator it=node->children.lower_bound("MESH"); it!=node->children.upper_bound("MESH"); ++it) {
      map<string,string>::const_iterator tmp = it->second->attributes.find("name");
      if (tmp == it->second->attributes.end()) {
	 cerr << "VLSVReader ERROR: XML tag does not contain attribute name!" << endl;
	 return false;
      }
      meshNames.push_back(tmp->second);
   }   
   return true;
}

bool VLSVReader::getBlockVariableNames(const std::string& meshName,std::list<std::string>& varNames) const {
   varNames.clear();
   if (fileOpen == false) return false;
   
   XMLNode* node = xmlReader.find("VLSV");
   if (node == NULL) return false;
   
   for (multimap<string,XMLNode*>::const_iterator it=node->children.lower_bound("BLOCKVARIABLE"); it!=node->children.upper_bound("BLOCKVARIABLE"); ++it) {
      map<string,string>::const_iterator tmp = it->second->attributes.find("mesh");
      if (tmp == it->second->attributes.end()) {
	 cerr << "VLSVReader ERROR: XML tag does not contain attribute mesh!" << endl;
	 return false;
      }
      
      if (tmp->second == meshName) {
	 tmp = it->second->attributes.find("name");
	 if (tmp == it->second->attributes.end()) {
	    cerr << "VLSVReader ERROR: XML tag does not contain attribute name!" << endl;
	 } else {
	    varNames.push_back(tmp->second);
	 }
      }
   }
   return true;
}
   

bool VLSVReader::getVariableNames(const std::string& meshName,std::list<std::string>& varNames) const {
   varNames.clear();
   if (fileOpen == false) return false;
   
   XMLNode* node = xmlReader.find("VLSV");
   if (node == NULL) return false;
   
   for (multimap<string,XMLNode*>::const_iterator it=node->children.lower_bound("VARIABLE"); it!=node->children.upper_bound("VARIABLE"); ++it) {
      map<string,string>::const_iterator tmp = it->second->attributes.find("mesh");
      if (tmp == it->second->attributes.end()) {
	 cerr << "VLSVReader ERROR: XML tag does not contain attribute mesh!" << endl;
	 return false;
      }
      
      if (tmp->second == meshName) {
	 tmp = it->second->attributes.find("name");
	 if (tmp == it->second->attributes.end()) {
	    cerr << "VLSVReader ERROR: XML tag does not contain attribute name!" << endl;
	 } else {
	    varNames.push_back(tmp->second);
	 }
      }
   }
   return true;
}

bool VLSVReader::loadArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs) {
   if (fileOpen == false) return false;
   
   // Find tag corresponding to given array:
   XMLNode* node = xmlReader.find(tagName,attribs);
   if (node == NULL) return false;
   
   // Find array name:
   list<pair<string,string> >::const_iterator arrayName = attribs.end();
   for (list<pair<string,string> >::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) {
      if (it->first == "name") arrayName = it;
   }
   if (arrayName == attribs.end()) {
      cerr << "VLSVReader ERROR: attribs does not contain array name!" << endl;
      return false;
   }
   
   // Copy array information from tag:
   arrayOpen.offset = atoi(node->value.c_str());
   arrayOpen.tagName = tagName;
   arrayOpen.arrayName = arrayName->second;
   arrayOpen.arraySize = atoi(node->attributes["arraysize"].c_str());
   arrayOpen.vectorSize = atoi(node->attributes["vectorsize"].c_str());
   arrayOpen.dataSize = atoi(node->attributes["datasize"].c_str());
   if (node->attributes["datatype"] == "int") arrayOpen.dataType = VLSV::INT;
   else if (node->attributes["datatype"] == "uint") arrayOpen.dataType = VLSV::UINT;
   else if (node->attributes["datatype"] == "float") arrayOpen.dataType = VLSV::FLOAT;
   else {
      cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
      return false;
   }   
   if (arrayOpen.arraySize == 0) return false;
   if (arrayOpen.vectorSize == 0) return false;
   if (arrayOpen.dataSize == 0) return false;
   
   return true;
}

bool VLSVReader::open(const std::string& fname) {
   bool success = true;
   filein.open(fname.c_str(), fstream::in);
   if (filein.good() == true) {
      fileName = fname;
      fileOpen = true;
   } else {
      filein.close();
      success = false;
   }
   if (success == false) return success;
   
   // Detect file endianness:
   char* ptr = reinterpret_cast<char*>(&endiannessFile);
   filein.read(ptr,1);
   if (endiannessFile != endiannessReader) swapIntEndianness = true;

   // Read footer offset:
   uint64_t footerOffset;
   char buffer[16];
   filein.seekg(8);
   filein.read(buffer,8);
   footerOffset = convUInt64(buffer,swapIntEndianness);
   
   // Read footer XML tree:
   filein.seekg(footerOffset);
   xmlReader.read(filein);
   filein.clear();
   filein.seekg(16);
   
   //xmlReader.print(cout);

   
   
   return success;
}

bool VLSVReader::readArray(
   const std::string& tagName,
   const std::list<std::pair<std::string,std::string> >& attribs,
   const uint64_t& begin,
   const uint64_t& amount,
   char* buffer
) {
   if (fileOpen == false) {
      cerr << "VLSVReader ERROR: readArray called but a file is not open!" << endl;
      return false;
   }
   
   // If array info has not been loaded, find it from XML tree:
   list<pair<string,string> >::const_iterator arrayName = attribs.end();
   for (list<pair<string,string> >::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) {
      if (it->first == "name") arrayName = it;
   }
   if (arrayName == attribs.end()) {
      cerr << "VLSVReader ERROR: attribs does not contain array name!" << endl;
      return false;
   }
   
   if (arrayOpen.tagName != tagName || arrayOpen.arrayName != arrayName->second) {
      // Find tag corresponding to given array:
      XMLNode* node = xmlReader.find(tagName,attribs);
      if (node == NULL) {
	 cerr << "VLSVReader ERROR: Failed to find tag='" << tagName << "' attribs:" << endl;
	 for (list<pair<string,string> >::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) {
	    cerr << '\t' << it->first << " = '" << it->second << "'" << endl;
	 }
	 return false;
      }
      
      // Copy array information from tag:
      arrayOpen.offset = atoi(node->value.c_str());
      arrayOpen.tagName = tagName;
      arrayOpen.arrayName = arrayName->second;
      arrayOpen.arraySize = atoi(node->attributes["arraysize"].c_str());
      arrayOpen.vectorSize = atoi(node->attributes["vectorsize"].c_str());
      arrayOpen.dataSize = atoi(node->attributes["datasize"].c_str());
      if (node->attributes["datatype"] == "int") arrayOpen.dataType = VLSV::INT;
      else if (node->attributes["datatype"] == "uint") arrayOpen.dataType = VLSV::UINT;
      else if (node->attributes["datatype"] == "float") arrayOpen.dataType = VLSV::FLOAT;
      else {
	 cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
	 return false;
      }
      
      if (arrayOpen.arraySize == 0) return false;
      if (arrayOpen.vectorSize == 0) return false;
      if (arrayOpen.dataSize == 0) return false;
   }
   // Sanity check on values:
   if (begin + amount > arrayOpen.arraySize) {
      cerr << "VLSVReader ERROR: Requested read exceeds array size. begin: " << begin << " amount: " << amount << " size: " << arrayOpen.arraySize << endl;
      return false;
   }
   
   streamoff start = arrayOpen.offset + begin*arrayOpen.vectorSize*arrayOpen.dataSize;
   streamsize readBytes = amount*arrayOpen.vectorSize*arrayOpen.dataSize;
   filein.clear();
   filein.seekg(start);
   filein.read(buffer,readBytes);
   if (filein.gcount() != readBytes) {
      cerr << "VLSVReader ERROR: Failed to read requested amount of bytes!" << endl;
      return false;
   }
   return true;
}
   
bool VLSVReader::readArray(
   const std::string& tagName,
   const std::string& arrayName,
   const std::list<std::pair<std::string,std::string> >& attribs,
   const uint64_t& begin,
   const uint64_t& amount,
   char* buffer
) {
   if (fileOpen == false) {
      std::cerr << __FILE__ << ":" << __LINE__ << " File not open" << std::endl;
      return false;
   }

   // If array info has not been loaded, find it from XML tree:
   if (arrayOpen.tagName != tagName || arrayOpen.arrayName != arrayName) {
      // Find tag corresponding to given array:
      list<pair<string,string> > attribs2 = attribs;
      attribs2.push_back(make_pair("name",arrayName));
      XMLNode* node = xmlReader.find(tagName,attribs2);
      if (node == NULL) {
         std::cerr << __FILE__ << ":" << __LINE__ << " node == NULL" << std::endl;
         return false;
      }
      
      // Copy array information from tag:
      arrayOpen.offset = atoi(node->value.c_str());
      arrayOpen.tagName = tagName;
      arrayOpen.arrayName = arrayName;
      arrayOpen.arraySize = atoi(node->attributes["arraysize"].c_str());
      arrayOpen.vectorSize = atoi(node->attributes["vectorsize"].c_str());
      arrayOpen.dataSize = atoi(node->attributes["datasize"].c_str());
      if (node->attributes["datatype"] == "int") arrayOpen.dataType = VLSV::INT;
      else if (node->attributes["datatype"] == "uint") arrayOpen.dataType = VLSV::UINT;
      else if (node->attributes["datatype"] == "float") arrayOpen.dataType = VLSV::FLOAT;
      else {
	 cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
	 return false;
      }
      if (arrayOpen.arraySize == 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << " Array size == 0" << std::endl;
         return false;
      }
      if (arrayOpen.vectorSize == 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << " Vector size == 0" << std::endl;
         return false;
      }
      if (arrayOpen.dataSize == 0) {
         std::cerr << __FILE__ << ":" << __LINE__ << " Data size == 0" << std::endl;
         return false;
      }
   }
   
   // Sanity check on values:
   if (begin + amount > arrayOpen.arraySize) {
      std::cerr << __FILE__ << ":" << __LINE__ << " Too much data requested" << std::endl;
      return false;
   }
   
   streamoff start = arrayOpen.offset + begin*arrayOpen.vectorSize*arrayOpen.dataSize;
   streamsize readBytes = amount*arrayOpen.vectorSize*arrayOpen.dataSize;
   filein.clear();
   filein.seekg(start);
   filein.read(buffer,readBytes);
   if (filein.gcount() != readBytes) {
      std::cerr << __FILE__ << ":" << __LINE__ << " Reading of file failed" << std::endl;
      return false;
   }
   
   // Swap endianness, if necessary:
   
   return true;
}
   
bool VLSVReader::readArray(
   const std::string& tagName,
   const std::string& arrayName,
   const uint64_t& begin,
   const uint64_t& amount,
   char* buffer
) {
   if (fileOpen == false) return false;

   // If array info has not been loaded, find it from XML tree:
   if (arrayOpen.tagName != tagName || arrayOpen.arrayName != arrayName) {
      // Find tag corresponding to given array:
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",arrayName));
      XMLNode* node = xmlReader.find(tagName,attribs);
      if (node == NULL) return false;

      // Copy array information from tag:
      arrayOpen.offset = atoi(node->value.c_str());
      arrayOpen.tagName = tagName;
      arrayOpen.arrayName = arrayName;	
      arrayOpen.arraySize = atoi(node->attributes["arraysize"].c_str());
      arrayOpen.vectorSize = atoi(node->attributes["vectorsize"].c_str());
      arrayOpen.dataSize = atoi(node->attributes["datasize"].c_str());
      if (node->attributes["datatype"] == "int") arrayOpen.dataType = VLSV::INT;
      else if (node->attributes["datatype"] == "uint") arrayOpen.dataType = VLSV::UINT;
      else if (node->attributes["datatype"] == "float") arrayOpen.dataType = VLSV::FLOAT;
      else {
	 cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
	 return false;
      }
      if (arrayOpen.arraySize == 0) return false;
      if (arrayOpen.vectorSize == 0) return false;
      if (arrayOpen.dataSize == 0) return false;
   }
   
   // Sanity check on values:
   if (begin + amount > arrayOpen.arraySize) return false;
   
   streamoff start = arrayOpen.offset + begin*arrayOpen.vectorSize*arrayOpen.dataSize;
   streamsize readBytes = amount*arrayOpen.vectorSize*arrayOpen.dataSize;
   filein.clear();
   filein.seekg(start);
   filein.read(buffer,readBytes);
   if (filein.gcount() != readBytes) return false;
   
   // Swap endianness, if necessary:
   
   return true;
}


// ********************************
// ***** VLSV PARALLEL READER *****
// ********************************

VLSVParReader::VLSVParReader(): VLSVReader() { }

VLSVParReader::~VLSVParReader() {
   close();
}

bool VLSVParReader::close() {
   if (fileOpen == false) return true;
   MPI_File_close(&filePtr);
   
   if (myRank == masterRank) filein.close();
   
   return true;
}

bool VLSVParReader::getArrayInfoMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
				       uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& dataSize) {
   if (myRank != masterRank) {
      cerr << "(VLSVPARREADER): getArrayInfoMaster called on process #" << myRank << endl;
      exit(1);
   }
   return VLSVReader::getArrayInfo(tagName,attribs,arraySize,vectorSize,dataType,dataSize);
}

bool VLSVParReader::getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs) {
   bool success = true;
   if (myRank == masterRank) {
      //success = getArrayInfoMaster(tagName,attribs,arrayOpen.arraySize,arrayOpen.vectorSize,arrayOpen.dataType,arrayOpen.dataSize);
      success = VLSVReader::loadArray(tagName,attribs);
   }
   // Check that master read array info correctly:
   int globalSuccess = 0;
   if (success == true) globalSuccess = 1;
   MPI_Bcast(&globalSuccess,1,MPI_Type<int>(),masterRank,comm);
   if (globalSuccess == 0) success = false;
   if (success == false) return false;
   
   // Master broadcasts array info to all processes:
   MPI_Bcast(&arrayOpen.offset,    1,MPI_Type<streamoff>(),masterRank,comm);
   MPI_Bcast(&arrayOpen.arraySize, 1,MPI_Type<uint64_t>(), masterRank,comm);
   MPI_Bcast(&arrayOpen.vectorSize,1,MPI_Type<uint64_t>(), masterRank,comm);
   MPI_Bcast(&arrayOpen.dataType,  1,MPI_Type<int>(),      masterRank,comm);
   MPI_Bcast(&arrayOpen.dataSize,  1,MPI_Type<uint64_t>(), masterRank,comm);
   return success;
}

bool VLSVParReader::getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
				 uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) {
   if (getArrayInfo(tagName,attribs) == false) return false;

   // Copy values to output variables:
   arraySize  = arrayOpen.arraySize;
   vectorSize = arrayOpen.vectorSize;
   dataType   = arrayOpen.dataType;
   byteSize   = arrayOpen.dataSize;
   return true;
}

bool VLSVParReader::multiReadAddUnit(const uint64_t& amount,char* buffer) {
   bool success = true;
   multiReadUnits.insert(make_pair(buffer,amount*arrayOpen.vectorSize*arrayOpen.dataSize));
   return success;
}

bool VLSVParReader::multiReadEnd(const uint64_t& offset) {
   bool success = true;
   
   // Get the byte size of multiread vector MPI datatype:
   MPI_Aint vectorLowerBound,vectorExtent;
   MPI_Type_get_extent(multiReadVectorType,&vectorLowerBound,&vectorExtent);
   
   // Create an MPI datatype for reading all data at once:
   const uint64_t N_reads = multiReadUnits.size();
   int* blockLengths  = new int[N_reads];
   int* displacements = new int[N_reads];
   
   displacements[0] = 0;
   blockLengths[0] = multiReadUnits.begin()->second;
   map<char*,uint64_t>::const_iterator it = multiReadUnits.begin();
   ++it;
   /*
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   if (myrank == 0) cerr << "read vector size = " << vectorExtent << endl;
   if (myrank == 0) cerr << "0" << '\t' << multiReadUnits.begin()->first << '\t' << displacements[0] << '\t' << blockLengths[0] << endl;
   cerr << "offset: " << offset << " arrayOpen.offset: " << arrayOpen.offset << endl;
   */
   uint64_t counter = 1;
   while (it != multiReadUnits.end()) {
      // Note: MPI_Type_indexed takes displacements in ints
      //const int baseOffset    = reinterpret_cast<size_t>(multiReadUnits.begin()->first);
      //const size_t offset     = reinterpret_cast<size_t>(it->first);
      
      blockLengths[counter]  = it->second;
      displacements[counter] = it->first - multiReadUnits.begin()->first;
      /*
      if (myrank == 0) {
	 cerr << counter << '\t' << it->first << '\t' << displacements[counter] << '\t' << blockLengths[counter] << endl;
      }
      */
      ++it;
      ++counter;
   }
   
   MPI_Datatype readType;
   MPI_Type_indexed(N_reads,blockLengths,displacements,MPI_Type<char>(),&readType);
   
   // Commit datatype and read everything in parallel:
   MPI_Type_commit(&readType);
   const uint64_t byteOffset = arrayOpen.offset + offset*arrayOpen.vectorSize*arrayOpen.dataSize;
   //cerr << "Proc #" << myrank << " reading from byte offset #" << byteOffset << endl;
   if (MPI_File_read_at_all(filePtr,byteOffset,multiReadUnits.begin()->first,1,readType,MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      //cerr << "ERROR in read!" << endl;
   }
   MPI_Type_free(&readType);
   
   delete blockLengths;
   delete displacements;   
   multiReadUnits.clear();
   return success;
}

bool VLSVParReader::multiReadStart(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs) {
   if (fileOpen == false) return false;
   bool success = true;
   multiReadUnits.clear();
   if (getArrayInfo(tagName,attribs) == false) return false;
   
   if (MPI_Type_contiguous(arrayOpen.vectorSize*arrayOpen.dataSize,MPI_Type<char>(),&multiReadVectorType) != MPI_SUCCESS) success = false;
   return success;
}

bool VLSVParReader::open(const std::string& fname,MPI_Comm comm,const int& masterRank,MPI_Info mpiInfo) {
   bool success = true;
   this->comm = comm;
   this->masterRank = masterRank;
   MPI_Comm_rank(comm,&myRank);
   MPI_Comm_size(comm,&processes);
   
   // Attempt to open the given input file using MPI:
   fileName = fname;
   int accessMode = MPI_MODE_RDONLY;
   if (MPI_File_open(comm,const_cast<char*>(fileName.c_str()),accessMode,mpiInfo,&filePtr) != MPI_SUCCESS) success = false;
   else fileOpen = true;
   
   // Only master process reads file footer and endianness. This is done using 
   // VLSVReader open() member function:
   if (myRank == masterRank) {
      VLSVReader::open(fname);
   }
   
   // Check that all processes have opened the file successfully:
   uchar globalSuccess = 0;
   if (success == true) globalSuccess = 1;
   uchar* results = new uchar[processes];
   MPI_Allgather(&globalSuccess,1,MPI_Type<uchar>(),results,1,MPI_Type<uchar>(),comm);
   for (int i=0; i<processes; ++i) if (results[i] == 0) success = false;
   delete results;
   results = NULL;
   if (success == false) return success;
   
   // Broadcast file endianness to all processes:
   MPI_Bcast(&endiannessFile,1,MPI_Type<uchar>(),masterRank,comm);
   
   return success;
}

bool VLSVParReader::readArrayMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
				    const uint64_t& begin,const uint64_t& amount,char* buffer) {
   if (myRank != masterRank) {
      cerr << "(VLSVPARREADER) readArrayMaster erroneously called on process #" << myRank << endl;
      exit(1);
   }
   // VLSVreadArray reads offset from XML tree into master only
   return VLSVReader::readArray(tagName,attribs,begin,amount,buffer);
}

bool VLSVParReader::readArray(
   const std::string& tagName,
   const std::list<std::pair<std::string,std::string> >& attribs,
   const uint64_t& begin,
   const uint64_t& amount,
   char* buffer
) {
   if (fileOpen == false) return false;
   bool success = true;

   // Fetch array info to all processes:
   if (getArrayInfo(tagName,attribs) == false) return false;
   const MPI_Offset start = arrayOpen.offset + begin*arrayOpen.vectorSize*arrayOpen.dataSize;
   const int readBytes    = amount*arrayOpen.vectorSize*arrayOpen.dataSize;
   
   // Read data on all processes in parallel:
   if (MPI_File_read_at_all(filePtr,start,buffer,readBytes,MPI_Type<char>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      cerr << "(VLSVPARREADER) MPI_File_read_at_all failed!" << endl;
      success = false;
   }
   return success;
}









