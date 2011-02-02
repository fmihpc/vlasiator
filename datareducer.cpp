#include <cstdlib>
#include <iostream>

#include "datareducer.h"

using namespace std;

static unsigned int dataReducers = 0;

/** Constructor for class DataReducer. Increases the value of DataReducer::dataReducers by one.
 */
DataReducer::DataReducer() { 
   ++dataReducers;
}

/** Destructor for class DataReducer. Reduces the value of DataReducer::dataReducers by one, 
 * and if after the reduction DataReducer::dataReducers equals zero all stored DataReductionOperators 
 * are deleted.
 */
DataReducer::~DataReducer() {
   --dataReducers;
   if (dataReducers != 0) return;
   
   // Call delete for each DataReductionOperator:
   for (vector<DRO::DataReductionOperator*>::iterator it=operators.begin(); it!=operators.end(); ++it) {
      delete *it;
      *it = NULL;
   }
}

/** Add a new DRO::DataReductionOperator which has been created with new operation. 
 * DataReducer will take care of deleting it.
 * @return If true, the given DRO::DataReductionOperator was added successfully.
 */
bool DataReducer::addOperator(DRO::DataReductionOperator* op) {
   operators.push_back(op);
   return true;
}

/** Asks each DRO::DataReductionOperator to write its data into the given byte array.
 * @param cell SpatialCell whose distribution function was reduced.
 * @param byteArray Byte array where reduced data should be inserted.
 * @return If true, all data was appended successfully.
 * @see DataReducer::reduceData.
 */
bool DataReducer::appendReducedData(const SpatialCell& cell,unsigned char* const byteArray) {
   unsigned char* ptr = byteArray;
   for (vector<DRO::DataReductionOperator*>::iterator it=operators.begin(); it!=operators.end(); ++it) {
      if ((*it)->appendReducedData(ptr) == false) return false;
      ptr += (*it)->getOutputByteSize();
   }
   return true;
}

bool DataReducer::appendReducedData(const SpatialCell& cell,std::vector<unsigned char>& byteArray) {
   for (vector<DRO::DataReductionOperator*>::iterator it=operators.begin(); it!=operators.end(); ++it) {
      //if ((*it)->appendReducedData(ptr) == false) return false;
      if ((*it)->appendReducedData(byteArray) == false) {
	 cerr << "DRO named '" << (*it)->getName() << "' failed" << endl;
	 return false;
      }
      //ptr += (*it)->getOutputByteSize();
   }
   return true;
}

/** Get the total size of the output data from all DRO::DataReductionOperators, in bytes.
 * @return The output size of reduced data for one SpatialCell.
 */
unsigned int DataReducer::getByteSize() const {
   unsigned int size = 0;
   for (vector<DRO::DataReductionOperator*>::const_iterator it=operators.begin(); it!=operators.end(); ++it) {
      size += (*it)->getOutputByteSize();
   }
   return size;
}

/** Get description of reduced data from each stored DRO::DataReductionOperator, and write 
 * the description and its size (in bytes) to the given variables. Array byteArray 
 * is allocated here, so the user has to take care of deleting it elsewhere.
 * The byteArray should be written into a .vlsv file as it is.
 * 
 * @param byteArray A pointer to a byte array which should be allocated, and 
 * where the descriptions should be written to.
 * @param arraySize The size of allocated byteArray.
 * @return If true, the descriptions were written successfully. 
 * It is always safe to delete byteArray.
 * @see VlsWriter for more information on vlsv file format.
 */
bool DataReducer::getDescription(unsigned char*& byteArray,unsigned int& arraySize) {
   byteArray = NULL;
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

   // Copy contents of temp array to byteArray:
   delete byteArray;
   byteArray = new unsigned char[array.size()];
   for (size_t i=0; i<array.size(); ++i) byteArray[i] = array[i];
   arraySize = array.size();
   
   return true;
}

/** Get the byte size of an entry, which gives the byte size of a variable name 
 * (which is a string). 
 * @return The size of a field which contains the size of a variable name.
 */
unsigned char DataReducer::getNameSizeEntryByteSize() {return sizeof(unsigned int);}

/** Reduce the given SpatialCell data. This function 
 * calls DRO::DataReductionOperator::reduceData member function of each
 * DRO::DataReductionOperator stored in DataReducer for the given SpatialCell.
 * @param cell SpatialCell whose distribution function should be reduced.
 * @return If true, all reduced data was calculated successfully.
 */
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
