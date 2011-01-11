#include <cstdlib>
#include <iostream>

#include "datareducer.h"

using namespace std;

static unsigned int dataReducers = 0;
//std::vector<DataReductionOperator*> DataReducer::operators;

DataReducer::DataReducer() { 
   ++dataReducers;
}

DataReducer::~DataReducer() {
   --dataReducers;
   if (dataReducers != 0) return;
   
   // Call delete for each DataReductionOperator:
   for (vector<DRO::DataReductionOperator*>::iterator it=operators.begin(); it!=operators.end(); ++it) {
      delete *it;
      *it = NULL;
   }
}

bool DataReducer::addOperator(DRO::DataReductionOperator* op) {
   operators.push_back(op);
   return true;
}

bool DataReducer::appendReducedData(const SpatialCell& cell,unsigned char* const byteArray) {
   unsigned char* ptr = byteArray;
   for (vector<DRO::DataReductionOperator*>::iterator it=operators.begin(); it!=operators.end(); ++it) {
      if ((*it)->appendReducedData(ptr) == false) return false;
      ptr += (*it)->getOutputByteSize();
   }
   return true;
}

unsigned int DataReducer::getByteSize() const {
   unsigned int size = 0;
   for (vector<DRO::DataReductionOperator*>::const_iterator it=operators.begin(); it!=operators.end(); ++it) {
      size += (*it)->getOutputByteSize();
   }
   return size;
}

bool DataReducer::getDescription(unsigned char*& byteArray,unsigned int& arraySize) {
   //cerr << "DataReducer::getDescription called" << endl;
   
   vector<unsigned char> array; // Temp array 
   
   unsigned int nameByteSize;
   string varName;
   unsigned char varType;
   unsigned char elementByteSize;
   const unsigned char* ptr;
   
   // Request necessary information from each DataReductionOperator and append the data 
   // into array as unsigned chars:
   for (vector<DRO::DataReductionOperator*>::const_iterator it=operators.begin(); it!=operators.end(); ++it) {
      varName = (*it)->getName();
      nameByteSize = varName.size() + 1; // +1 is for end-of-character.
      varType = (*it)->getVariableType();
      elementByteSize = (*it)->getElementByteSize();
      
      ptr = reinterpret_cast<unsigned char*>(&nameByteSize);
      for (size_t i=0; i<sizeof(nameByteSize); ++i) array.push_back(ptr[i]);
      
      for (size_t i=0; i<varName.size(); ++i) array.push_back(static_cast<unsigned char>(varName[i]));
      array.push_back(static_cast<unsigned char>('\0'));
      array.push_back(varType);
      array.push_back(elementByteSize);
   }
   // nameByteSize=0 signals the end of variable descriptions:
   nameByteSize = 0;
   ptr = reinterpret_cast<unsigned char*>(&nameByteSize);
   for (size_t i=0; i<sizeof(nameByteSize); ++i) array.push_back(ptr[i]);
   /*
   for (size_t i=0; i<array.size(); ++i) {
      cerr << i << ": '" << (int)array[i] << "'" << endl;
   }
   */
   // Copy contents of temp array to byteArray:
   delete byteArray;
   byteArray = new unsigned char[array.size()];
   for (size_t i=0; i<array.size(); ++i) byteArray[i] = array[i];
   arraySize = array.size();
   
   return true;
}

unsigned char DataReducer::getNameByteSize() {return sizeof(unsigned int);}

bool DataReducer::reduceData(const SpatialCell& cell) {
   // Tell each reduction operator which spatial cell we are dealing with:
   for (vector<DRO::DataReductionOperator*>::iterator it=operators.begin(); it!=operators.end(); ++it) {
      if ((*it)->setSpatialCell(cell) == false) return false;
   }
   
   // Go through all velocity blocks in the spatial cell:
   for (unsigned int i=0; i<cell.N_blocks; ++i) {
      // Send each velocity block to data reduction operator:
      for (vector<DRO::DataReductionOperator*>::iterator it=operators.begin(); it!=operators.end(); ++it) {
	 if ((*it)->reduceData(cell.cpu_avgs+i*SIZE_VELBLOCK,cell.cpu_blockParams+i*SIZE_BLOCKPARAMS) == false) return false;
      }
   }
   return true;
}
