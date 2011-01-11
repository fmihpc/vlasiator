#ifndef DATAREDUCER_H
#define DATAREDUCER_H

#include <vector>

#include "datareductionoperator.h"

class DataReducer {
 public:
   DataReducer();
   ~DataReducer();
   
   bool addOperator(DRO::DataReductionOperator* op);
   bool appendReducedData(const SpatialCell& cell,unsigned char* const byteArray);
   unsigned int getByteSize() const;
   bool getDescription(unsigned char*& byteArray,unsigned int& arraySize);
   static unsigned char getNameByteSize();
   bool reduceData(const SpatialCell& cell);
   
 private:
   DataReducer(const DataReducer& dr);
   std::vector<DRO::DataReductionOperator*> operators;
};

#endif
