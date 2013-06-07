/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












*/

#ifndef DATAREDUCER_H
#define DATAREDUCER_H

#include <vector>

#include "../spatial_cell.hpp"
#include "datareductionoperator.h"

/** The purpose of DataReducer is to contain DRO::DataReductionOperators, and apply 
 * them to simulation data when writing output files. Files containing full 
 * distribution functions of every spatial cell require so much disk space 
 * that they cannot be written out so often as user would want. Thus, derived 
 * quantities need to be calculated for every spatial cell, which are then 
 * written to data files. This process is here called data reduction.
 */
class DataReducer {
 public:
   DataReducer();
   ~DataReducer();
   
   bool addOperator(DRO::DataReductionOperator* op);
   bool getDataVectorInfo(const unsigned int& operatorID,std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
   std::string getName(const unsigned int& operatorID) const;
   bool reduceData(const SpatialCell* cell,const unsigned int& operatorID,char* buffer);
   bool reduceData(const SpatialCell* cell,const unsigned int& operatorID,Real * result);
   unsigned int size() const;
   
 private:
   /** Private copy-constructor to prevent copying the class.
    */
   DataReducer(const DataReducer& dr);
   
   std::vector<DRO::DataReductionOperator*> operators;
   /**< A container for all DRO::DataReductionOperators stored in DataReducer.*/
};

void initializeDataReducers(DataReducer * outputReducer, DataReducer * diagnosticReducer);

#endif
