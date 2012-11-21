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
#include "../common.h"
using namespace std;

void initializeDataReducers(DataReducer * outputReducer, DataReducer * diagnosticReducer)
{
   typedef Parameters P;
         
   vector<string>::const_iterator it;
   for (it = P::outputVariableList.begin();
        it != P::outputVariableList.end();
        it++) {
      if(*it == "B")
         outputReducer->addOperator(new DRO::VariableB);
      if(*it == "BackgroundB")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("background_B",CellParams::BGBX,3));
      if(*it == "PerturbedB")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("perturbed_B",CellParams::PERBX,3));
      if(*it == "E")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("E",CellParams::EX,3));
      if(*it == "Rho")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho",CellParams::RHO,1));
      if(*it == "RhoV")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho_v",CellParams::RHOVX,3));
      if(*it == "RhoLossAdjust")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho_loss_adjust",CellParams::RHOLOSSADJUST,1));
      if(*it == "RhoLossVelBoundary")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho_loss_velocity_boundary",CellParams::RHOLOSSVELBOUNDARY,1));
      if(*it == "LBweight")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("LB_weight",CellParams::LBWEIGHTCOUNTER,1));
      if(*it == "MaxVdt")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_v_dt",CellParams::MAXVDT,1));
      if(*it == "MaxRdt")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_r_dt",CellParams::MAXRDT,1));
      if(*it == "MaxFieldsdt")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_fields_dt",CellParams::MAXFDT,1));
      if(*it == "MPIrank")
         outputReducer->addOperator(new DRO::MPIrank);
      if(*it == "BoundaryType")
         outputReducer->addOperator(new DRO::BoundaryType);
      if(*it == "BoundaryLayer")
         outputReducer->addOperator(new DRO::BoundaryLayer);
      if(*it == "Blocks")
         outputReducer->addOperator(new DRO::Blocks);
      if(*it == "VolE")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("E_vol",CellParams::EXVOL,3));
      if(*it == "VolB")
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("B_vol",CellParams::BXVOL,3));
      if(*it == "Pressure")
         outputReducer->addOperator(new DRO::VariablePressure);
      if(*it == "PTensor") {
         outputReducer->addOperator(new DRO::VariablePTensorDiagonal);
         outputReducer->addOperator(new DRO::VariablePTensorOffDiagonal);
      }
      if(*it == "derivs") {
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("drhodx",fieldsolver::drhodx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("drhody",fieldsolver::drhody,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("drhodz",fieldsolver::drhodz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBxdy",fieldsolver::dBxdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBxdz",fieldsolver::dBxdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBydx",fieldsolver::dBydx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBydz",fieldsolver::dBydz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBzdx",fieldsolver::dBzdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBzdy",fieldsolver::dBzdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVxdx",fieldsolver::dVxdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVxdy",fieldsolver::dVxdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVxdz",fieldsolver::dVxdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVydx",fieldsolver::dVydx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVydy",fieldsolver::dVydy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVydz",fieldsolver::dVydz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVzdx",fieldsolver::dVzdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVzdy",fieldsolver::dVzdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVzdz",fieldsolver::dVzdz,1));
      }
      if(*it == "BVOLderivs") {
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBXVOLdy",bvolderivatives::dBXVOLdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBXVOLdz",bvolderivatives::dBXVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBYVOLdx",bvolderivatives::dBYVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBYVOLdz",bvolderivatives::dBYVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBZVOLdx",bvolderivatives::dBZVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBZVOLdy",bvolderivatives::dBZVOLdy,1));
      }
   }
   
   for (it = P::diagnosticVariableList.begin();
        it != P::diagnosticVariableList.end();
        it++) {
      if(*it == "FluxB")
         diagnosticReducer->addOperator(new DRO::DiagnosticFluxB);
      if(*it == "FluxE")
         diagnosticReducer->addOperator(new DRO::DiagnosticFluxE);
      if(*it == "Blocks")
         diagnosticReducer->addOperator(new DRO::Blocks);
      if(*it == "Rho")
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho",CellParams::RHO,1));
      if(*it == "RhoLossAdjust")
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho_loss_adjust",CellParams::RHOLOSSADJUST,1));
      if(*it == "RhoLossVelBoundary")
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho_loss_velocity_boundary",CellParams::RHOLOSSVELBOUNDARY,1));
      if(*it == "LBweight")
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("LB_weight",CellParams::LBWEIGHTCOUNTER,1));
      if(*it == "MaxVdt")
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_v_dt",CellParams::MAXVDT,1));
      if(*it == "MaxRdt")
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_r_dt",CellParams::MAXRDT,1));
      if(*it == "MaxFieldsdt")
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_fields_dt",CellParams::MAXFDT,1));
      if(*it == "MaxDistributionFunction")
         diagnosticReducer->addOperator(new DRO::MaxDistributionFunction);
      if(*it == "MinDistributionFunction")
         diagnosticReducer->addOperator(new DRO::MinDistributionFunction);
      if(*it == "BoundaryType")
         diagnosticReducer->addOperator(new DRO::BoundaryType);
      if(*it == "BoundaryLayer")
         diagnosticReducer->addOperator(new DRO::BoundaryLayer);
   }
}

// ************************************************************
// ***** DEFINITIONS FOR DATAREDUCER CLASS *****
// ************************************************************

/** Constructor for class DataReducer.
 */
DataReducer::DataReducer() { }

/** Destructor for class DataReducer. All stored DataReductionOperators 
 * are deleted.
 */
DataReducer::~DataReducer() {
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

