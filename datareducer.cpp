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

#include "datareducer.h"
#include "readparameters.h"

using namespace std;

void initializeDataReducers(DataReducer * outputReducer, DataReducer * diagnosticReducer)
{
   typedef Readparameters RP;
   typedef Parameters P;
   
   vector<string>::const_iterator it;
   for (it = P::outputVariableList.begin();
	it != P::outputVariableList.end();
   it++) {
      if(*it == "B")
	 outputReducer->addOperator(new DRO::VariableB);
      if(*it == "E")
	 outputReducer->addOperator(new DRO::VariableE);
      if(*it == "Rho")
	 outputReducer->addOperator(new DRO::VariableRho);
      if(*it == "RhoV")
	 outputReducer->addOperator(new DRO::VariableRhoV);
      if(*it == "RhoLossAdjust")
	 outputReducer->addOperator(new DRO::VariableRhoLossAdjust);
      if(*it == "RhoLossVelBoundary")
	 outputReducer->addOperator(new DRO::VariableRhoLossVelBoundary);
      if(*it == "MPIrank")
	 outputReducer->addOperator(new DRO::MPIrank);
      if(*it == "Blocks")
	 outputReducer->addOperator(new DRO::Blocks);
      if(*it == "VolE")
	 outputReducer->addOperator(new DRO::VariableVolE);
      if(*it == "VolB")
	 outputReducer->addOperator(new DRO::VariableVolB);
      if(*it == "Pressure")
	 outputReducer->addOperator(new DRO::VariablePressure);
      if(*it == "PTensor") {
	 outputReducer->addOperator(new DRO::VariablePTensorDiagonal);
	 outputReducer->addOperator(new DRO::VariablePTensorOffDiagonal);
      }      
      if(*it == "dBxdz")
	 outputReducer->addOperator(new DRO::VariabledBxdz);
   }
   
   for (it = P::diagnosticVariableList.begin();
	it != P::diagnosticVariableList.end();
   it++) {
      if(*it == "FluxB")
	 diagnosticReducer->addOperator(new DRO::DiagnosticFluxB);
      if(*it == "Blocks")
	 diagnosticReducer->addOperator(new DRO::Blocks);
      if(*it == "Rho")
	 diagnosticReducer->addOperator(new DRO::VariableRho);
      if(*it == "RhoLossAdjust")
	 diagnosticReducer->addOperator(new DRO::VariableRhoLossAdjust);
      if(*it == "RhoLossVelBoundary")
	 diagnosticReducer->addOperator(new DRO::VariableRhoLossVelBoundary);
      if(*it == "MaxVi")
	 diagnosticReducer->addOperator(new DRO::MaxVi);
   }
}

// ************************************************************
// ***** DEFINITIONS FOR DATAREDUCER CLASS *****
// ************************************************************

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

/** Get the name of a DataReductionOperator.
 * @param operatorID ID number of the operator whose name is requested.
 * @return Name of the operator.
 */
std::string DataReducer::getName(const unsigned int& operatorID) const {
   if (operatorID >= operators.size()) return "";
   return operators[operatorID]->getName();
}

/** Get info on the type of data calculated by the given DataReductionOperator.
 * A DataReductionOperator writes an array on disk. Each element of the array is a vector with n elements. Finally, each
 * vector element has a byte size, as given by the sizeof function.
 * @param operatorID ID number of the DataReductionOperator whose output data info is requested.
 * @param dataType Basic datatype, must be int, uint, or float.
 * @param dataSize Byte size of written datatype, for example double-precision floating points
 * have byte size of sizeof(double).
 * @param vectorSize How many elements are in the vector returned by the DataReductionOperator.
 * @return If true, DataReductionOperator was found and it returned sensible values.
 */
bool DataReducer::getDataVectorInfo(const unsigned int& operatorID,std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
   if (operatorID >= operators.size()) return false;
   return operators[operatorID]->getDataVectorInfo(dataType,dataSize,vectorSize);
}

/** Request a DataReductionOperator to calculate its output data and to write it to the given buffer.
 * @param cell Pointer to spatial cell whose data is to be reduced.
 * @param operatorID ID number of the applied DataReductionOperator.
 * @param buffer Buffer in which DataReductionOperator should write its data.
 * @return If true, DataReductionOperator calculated and wrote data successfully.
 */
bool DataReducer::reduceData(const SpatialCell* cell,const unsigned int& operatorID,char* buffer) {
   // Tell the chosen operator which spatial cell we are counting:
   if (operatorID >= operators.size()) return false;
   if (operators[operatorID]->setSpatialCell(cell) == false) return false;

   if (operators[operatorID]->reduceData(cell,buffer) == false) return false;
   return true;
}

/** Request a DataReductionOperator to calculate its output data and to write it to the given variable.
 * @param cell Pointer to spatial cell whose data is to be reduced.
 * @param operatorID ID number of the applied DataReductionOperator.
 * @param result Real variable in which DataReductionOperator should write its result.
 * @return If true, DataReductionOperator calculated and wrote data successfully.
 */
bool DataReducer::reduceData(const SpatialCell* cell,const unsigned int& operatorID,Real * result) {
   // Tell the chosen operator which spatial cell we are counting:
   if (operatorID >= operators.size()) return false;
   if (operators[operatorID]->setSpatialCell(cell) == false) return false;
   
   if (operators[operatorID]->reduceData(cell,result) == false) return false;
   return true;
}

/** Get the number of DataReductionOperators stored in DataReducer.
 * @return Number of DataReductionOperators stored in DataReducer.
 */
unsigned int DataReducer::size() const {return operators.size();}

