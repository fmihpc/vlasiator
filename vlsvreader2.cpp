#include <cstdlib>
#include <iostream>

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
   return true;
}

bool VLSVReader::getArrayInfo(const std::string& tagName,const std::string& arrayName,uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const {   
   if (fileOpen == false) return false;
   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name",arrayName));
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

bool VLSVReader::readArray(const std::string& tagName,const std::string& arrayName,const std::list<std::pair<std::string,std::string> >& attribs,const uint64_t& begin,const uint64_t& amount,char* buffer) {
   if (fileOpen == false) return false;
   
   // If array info has not been loaded, find it from XML tree:
   if (arrayOpen.tagName != tagName || arrayOpen.arrayName != arrayName) {
      // Find tag corresponding to given array:
      list<pair<string,string> > attribs2 = attribs;
      attribs2.push_back(make_pair("name",arrayName));
      XMLNode* node = xmlReader.find(tagName,attribs2);
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
   
bool VLSVReader::readArray(const std::string& tagName,const std::string& arrayName,const uint64_t& begin,const uint64_t& amount,char* buffer) {
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







