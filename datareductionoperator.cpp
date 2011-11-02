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
#include <mpi.h>
#include <iostream>
#include <limits>

#include "datareductionoperator.h"
#include "vlscommon.h"
#include "cpu/cpu_common.h"


using namespace std;

namespace DRO {

   // ************************************************************
   // ***** DEFINITIONS FOR DATAREDUCTIONOPERATOR BASE CLASS *****
   // ************************************************************
   
   /** DataReductionOperator base class constructor. The constructor is empty.*/
   DataReductionOperator::DataReductionOperator() { }

   /** DataReductionOperator base class virtual destructor. The destructor is empty.*/
   DataReductionOperator::~DataReductionOperator() { }

   /** Get info on the data the DataReductionOperator writes on disk. A DataReductionOperator writes 
    * an array on disk. Each element of the array is a vector with n elements. Finally, each 
    * vector element has a byte size, as given by the sizeof function.
    * @param dataType Basic datatype, must be int, uint, float
    * @param dataSize Byte size of written datatype, for example double-precision floating points
    * have byte size of sizeof(double).
    * @param vectorSize How many elements are in the vector returned by the DataReductionOperator.
    * @return If true, DataReductionOperator returned sensible values.
    */
   bool DataReductionOperator::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      cerr << "ERROR: DataReductionOperator::getDataVectorInfo called insted of derived class function!" << endl;
      return false;
   }
   
   /** Get the name of the reduced data variable. The name is written to the disk as-is 
    * and is used in visualization.
    * @return The name of the data. The base class function returns an empty string.
    */
   std::string DataReductionOperator::getName() const {
      cerr << "ERROR: DataReductionOperator::getName called instead of derived class function!" << endl;
      return string("");
   }
   
   /** Reduce the data and write the data vector to the given buffer.
    * @param N_blocks Number of velocity blocks in array avgs.
    * @param avgs Array containing distribution function values for each velocity block.
    * @param blockParams Array containing the parameters of each velocity block.
    * @param buffer Buffer in which the reduced data is written.
    * @return If true, DataReductionOperator reduced data successfully.
    */
   bool DataReductionOperator::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      cerr << "ERROR: DataReductionOperator::reduceData called instead of derived class function!" << endl;
      return false;
   }
   
   /** Set the SpatialCell whose data is going to be reduced by subsequent calls to 
    * DRO::DataReductionOperator::reduceData. This function is provided so that 
    * variables stored per SpatialCell can be accessed.
    * 
    * Spatial cell variables are stored in array SpatialCell::cpu_cellParams. 
    * The contents of array elements are stored in namespace CellParams. For example, 
    * cell.cpu_cellParams[%CellParams::EX] contains the electric field.
    * @param cell The SpatialCell whose data is to be reduced next.
    * @return If true, the SpatialCell was set correctly.
    */
   bool DataReductionOperator::setSpatialCell(const SpatialCell& cell) {
      cerr << "ERROR: DataReductionOperator::setSpatialCell called instead of derived class function!" << endl;
      return false;
   }
   
   VariableE::VariableE(): DataReductionOperator() { }
   VariableE::~VariableE() { }
   
   bool VariableE::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize = sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   std::string VariableE::getName() const {return "E";}
   
   bool VariableE::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(E);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableE::setSpatialCell(const SpatialCell& cell) {
      E = &(cell.cpu_cellParams[CellParams::EX]);
      Ex = cell.cpu_cellParams[CellParams::EX];
      Ey = cell.cpu_cellParams[CellParams::EY];
      Ez = cell.cpu_cellParams[CellParams::EZ];
      return true;
   }

   VariableVolE::VariableVolE(): DataReductionOperator() { }
   VariableVolE::~VariableVolE() { }
   
   bool VariableVolE::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
         
   std::string VariableVolE::getName() const {return "E_vol";}
      
      bool VariableVolE::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(E);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableVolE::setSpatialCell(const SpatialCell& cell) {
      E = &(cell.cpu_cellParams[CellParams::EXVOL]);
      return true;
   }
   
   VariableB::VariableB(): DataReductionOperator() { }
   VariableB::~VariableB() { }

   bool VariableB::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   std::string VariableB::getName() const {return "B";}
   
   bool VariableB::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(B);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableB::setSpatialCell(const SpatialCell& cell) {
      B  = &(cell.cpu_cellParams[CellParams::BX]);
      Bx = cell.cpu_cellParams[CellParams::BX];
      By = cell.cpu_cellParams[CellParams::BY];
      Bz = cell.cpu_cellParams[CellParams::BZ];
      return true;
   }

   VariableVolB::VariableVolB(): DataReductionOperator() { }
   VariableVolB::~VariableVolB() { }
   
   bool VariableVolB::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   std::string VariableVolB::getName() const {return "B_vol";}
   
   bool VariableVolB::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(B);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableVolB::setSpatialCell(const SpatialCell& cell) {
      B  = &(cell.cpu_cellParams[CellParams::BXVOL]);
      return true;
   }
   
   VariableRho::VariableRho(): DataReductionOperator() { }
   VariableRho::~VariableRho() { }
   
   std::string VariableRho::getName() const {return "rho";}
   
   bool VariableRho::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   bool VariableRho::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&rho);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRho::setSpatialCell(const SpatialCell& cell) {
      rho = cell.cpu_cellParams[CellParams::RHO];
      return true;
   }
   
   
   MPIrank::MPIrank(): DataReductionOperator() { }
   MPIrank::~MPIrank() { }
   
   bool MPIrank::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = 4;
      vectorSize = 1;
      return true;
   }
   
   std::string MPIrank::getName() const {return "MPI_rank";}
   
   bool MPIrank::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&mpiRank);
      for (uint i=0; i<sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool MPIrank::setSpatialCell(const SpatialCell& cell) {
      int intRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&intRank);
      rank = 1.0*intRank;
      mpiRank = intRank;
      return true;
   }

   VariableRhoV::VariableRhoV(): DataReductionOperator() { }
   VariableRhoV::~VariableRhoV() { }
   
   bool VariableRhoV::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   std::string VariableRhoV::getName() const {return "rho_v";}
   
   bool VariableRhoV::reduceData(const unsigned int& N_blocks,const Real* const avgs,const Real* const blockParams,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(rhov);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoV::setSpatialCell(const SpatialCell& cell) {
      rhov  = &(cell.cpu_cellParams[CellParams::RHOVX]);
      rhovx = cell.cpu_cellParams[CellParams::RHOVX];
      rhovy = cell.cpu_cellParams[CellParams::RHOVY];
      rhovz = cell.cpu_cellParams[CellParams::RHOVZ];
      return true;
   }
   
   // Added by YK
   VariablePressure::VariablePressure(): DataReductionOperator() { }
   VariablePressure::~VariablePressure() { }
   
   std::string VariablePressure::getName() const {return "Pressure";}
   
   bool VariablePressure::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
     dataType = "float";
     dataSize =  sizeof(Real);
     vectorSize = 1;
     return true;
   }
   
   // YK Adding pressure calculations to Vlasiator.
   // p = m/3 * integral((v - <V>)^2 * f(r,v) dV), doing the sum of the x, y and z components.
   bool VariablePressure::reduceData(const unsigned int& N_blocks, const Real* const avgs, const Real* const blockParams, char* buffer) {
      const Real HALF = 0.5;
      const Real THIRD = 1.0/3.0;
     
      Real nvx2_sum = 0.0;
      Real nvy2_sum = 0.0;
      Real nvz2_sum = 0.0;
      for (uint n=0; n<N_blocks; ++n)
	 for (uint k=0; k<WID; ++k)
	    for (uint j=0; j<WID; ++j)
	       for (uint i=0; i<WID; ++i) {
		  const Real VX = blockParams[n*SIZE_BLOCKPARAMS + BlockParams::VXCRD] + (i+HALF) * blockParams[n*SIZE_BLOCKPARAMS + BlockParams::DVX];
		  const Real VY = blockParams[n*SIZE_BLOCKPARAMS + BlockParams::VYCRD] + (j+HALF) * blockParams[n*SIZE_BLOCKPARAMS + BlockParams::DVY];
		  const Real VZ = blockParams[n*SIZE_BLOCKPARAMS + BlockParams::VZCRD] + (k+HALF) * blockParams[n*SIZE_BLOCKPARAMS + BlockParams::DVZ];
	 
		  nvx2_sum += avgs[n*SIZE_VELBLOCK + cellIndex(i,j,k)]*(VX - averageVX)*(VX - averageVX);
		  nvy2_sum += avgs[n*SIZE_VELBLOCK + cellIndex(i,j,k)]*(VY - averageVY)*(VY - averageVY);
		  nvz2_sum += avgs[n*SIZE_VELBLOCK + cellIndex(i,j,k)]*(VZ - averageVZ)*(VZ - averageVZ);
     }
     
     // Accumulate contributions coming from this velocity block to the 
     // spatial cell velocity moments. If multithreading / OpenMP is used, 
     // these updates need to be atomic:
     const Real DV3 = blockParams[BlockParams::DVX] * blockParams[BlockParams::DVY] * blockParams[BlockParams::DVZ];
     Pressure += physicalconstants::MASS_PROTON * THIRD * (nvx2_sum + nvy2_sum + nvz2_sum) * DV3;
     
     const char* ptr = reinterpret_cast<const char*>(&Pressure);
     for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
     return true;
   }
   
   bool VariablePressure::setSpatialCell(const SpatialCell& cell) {
      averageVX = cell.cpu_cellParams[CellParams::RHOVX] / cell.cpu_cellParams[CellParams::RHO];
      averageVY = cell.cpu_cellParams[CellParams::RHOVY] / cell.cpu_cellParams[CellParams::RHO];
      averageVZ = cell.cpu_cellParams[CellParams::RHOVZ] / cell.cpu_cellParams[CellParams::RHO];
      Pressure = 0.0;
      return true;
   }
} // namespace DRO
