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
#include <mpi.h>
#include <iostream>
#include <limits>

#include "datareductionoperator.h"
#include "vlscommon.h"
#include "vlasovsolver/cpu_common.h"


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
   bool DataReductionOperator::reduceData(const SpatialCell* cell,char* buffer) {
      cerr << "ERROR: DataReductionOperator::reduceData called instead of derived class function!" << endl;
      return false;
   }
   
   /** Reduce the data and write the data vector to the given variable.
    * @param N_blocks Number of velocity blocks in array avgs.
    * @param avgs Array containing distribution function values for each velocity block.
    * @param blockParams Array containing the parameters of each velocity block.
    * @param result Real variable in which the reduced data is written.
    * @return If true, DataReductionOperator reduced data successfully.
    */
   bool DataReductionOperator::reduceData(const SpatialCell* cell,Real* result) {
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
   bool DataReductionOperator::setSpatialCell(const SpatialCell* cell) {
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
   
   bool VariableE::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(E);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableE::setSpatialCell(const SpatialCell* cell) {
      E = &(cell->parameters[CellParams::EX]);
      Ex = cell->parameters[CellParams::EX];
      Ey = cell->parameters[CellParams::EY];
      Ez = cell->parameters[CellParams::EZ];
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
      
      bool VariableVolE::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(E);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableVolE::setSpatialCell(const SpatialCell* cell) {
      E = &(cell->parameters[CellParams::EXVOL]);
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
   
   bool VariableB::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(B);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableB::setSpatialCell(const SpatialCell* cell) {
      B  = &(cell->parameters[CellParams::BX]);
      Bx = cell->parameters[CellParams::BX];
      By = cell->parameters[CellParams::BY];
      Bz = cell->parameters[CellParams::BZ];
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
   
   bool VariableVolB::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(B);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableVolB::setSpatialCell(const SpatialCell* cell) {
      B  = &(cell->parameters[CellParams::BXVOL]);
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
   
   bool VariableRho::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&rho);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }

   bool VariableRho::reduceData(const SpatialCell* cell,Real* buffer){
      *buffer=rho;
      return true;
   }
   
   bool VariableRho::setSpatialCell(const SpatialCell* cell) {
      rho = cell->parameters[CellParams::RHO];
      return true;
   }


//RHOLOSSADJUST

   VariableRhoLossAdjust::VariableRhoLossAdjust(): DataReductionOperator() { }
   VariableRhoLossAdjust::~VariableRhoLossAdjust() { }

   std::string VariableRhoLossAdjust::getName() const {return "rho_loss_adjust";}

   bool VariableRhoLossAdjust::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   bool VariableRhoLossAdjust::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&rhoLoss);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoLossAdjust::reduceData(const SpatialCell* cell,Real* buffer){
      *buffer=rhoLoss;
      return true;
   }
   
   bool VariableRhoLossAdjust::setSpatialCell(const SpatialCell* cell) {
      rhoLoss = cell->parameters[CellParams::RHOLOSSADJUST];
      return true;
   }


   //RHOLOSSVELBOUNDARY

   VariableRhoLossVelBoundary::VariableRhoLossVelBoundary(): DataReductionOperator() { }
   VariableRhoLossVelBoundary::~VariableRhoLossVelBoundary() { }

   std::string VariableRhoLossVelBoundary::getName() const {return "rho_loss_velocity_boundary";}

   bool VariableRhoLossVelBoundary::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   bool VariableRhoLossVelBoundary::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&rhoLoss);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoLossVelBoundary::reduceData(const SpatialCell* cell,Real* buffer){
      *buffer=rhoLoss;
      return true;
   }
   
   bool VariableRhoLossVelBoundary::setSpatialCell(const SpatialCell* cell) {
      rhoLoss = cell->parameters[CellParams::RHOLOSSVELBOUNDARY];
      return true;
   }


//MPI rank
      
   
   MPIrank::MPIrank(): DataReductionOperator() { }
   MPIrank::~MPIrank() { }
   
   bool MPIrank::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = 4;
      vectorSize = 1;
      return true;
   }
   
   std::string MPIrank::getName() const {return "MPI_rank";}
   
   bool MPIrank::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&mpiRank);
      for (uint i=0; i<sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool MPIrank::setSpatialCell(const SpatialCell* cell) {
      int intRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&intRank);
      rank = 1.0*intRank;
      mpiRank = intRank;
      return true;
   }




   Blocks::Blocks(): DataReductionOperator() { }
   Blocks::~Blocks() { }
   
   bool Blocks::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = 4;
      vectorSize = 1;
      return true;
   }
   
   std::string Blocks::getName() const {return "Blocks";}
   
   bool Blocks::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&nBlocks);
      for (uint i=0; i<sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool Blocks::reduceData(const SpatialCell* cell,Real* buffer) {
      *buffer = 1.0 * nBlocks;
      return true;
   }
   
   bool Blocks::setSpatialCell(const SpatialCell* cell) {
      nBlocks = cell->number_of_blocks;
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
   
   bool VariableRhoV::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(rhov);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoV::setSpatialCell(const SpatialCell* cell) {
      rhov  = &(cell->parameters[CellParams::RHOVX]);
      rhovx = cell->parameters[CellParams::RHOVX];
      rhovy = cell->parameters[CellParams::RHOVY];
      rhovz = cell->parameters[CellParams::RHOVZ];
      return true;
   }
   
   // Scalar pressure added by YK
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
   bool VariablePressure::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      const Real THIRD = 1.0/3.0;
     
      Real nvx2_sum = 0.0;
      Real nvy2_sum = 0.0;
      Real nvz2_sum = 0.0;

      for(uint n=0; n<cell->number_of_blocks; n++) {
	 unsigned int blockId = cell->velocity_block_list[n];
	 const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
	 for (uint k=0; k<WID; ++k)
	    for (uint j=0; j<WID; ++j)
	       for (uint i=0; i<WID; ++i) {
		  const Real VX = block-> parameters[BlockParams::VXCRD] + (i+HALF) * block-> parameters[BlockParams::DVX];
		  const Real VY = block-> parameters[BlockParams::VYCRD] + (j+HALF) * block-> parameters[BlockParams::DVY];
		  const Real VZ = block-> parameters[BlockParams::VZCRD] + (k+HALF) * block-> parameters[BlockParams::DVZ];
		  
		  const Real DV3 = block-> parameters[BlockParams::DVX] * block-> parameters[BlockParams::DVY] * block-> parameters[BlockParams::DVZ];
		  
		  nvx2_sum += block-> data[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX) * DV3;
		  nvy2_sum += block-> data[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY) * DV3;
		  nvz2_sum += block-> data[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ) * DV3;
	       }
     }
      
     // Accumulate contributions coming from this velocity block to the 
     // spatial cell velocity moments. If multithreading / OpenMP is used, 
     // these updates need to be atomic:
     
     Pressure += physicalconstants::MASS_PROTON * THIRD * (nvx2_sum + nvy2_sum + nvz2_sum);
     
     const char* ptr = reinterpret_cast<const char*>(&Pressure);
     for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
     return true;
   }
   
   bool VariablePressure::setSpatialCell(const SpatialCell* cell) {
      if(cell-> parameters[CellParams::RHO] != 0.0) {
	 averageVX = cell-> parameters[CellParams::RHOVX] / cell-> parameters[CellParams::RHO];
	 averageVY = cell-> parameters[CellParams::RHOVY] / cell-> parameters[CellParams::RHO];
	 averageVZ = cell-> parameters[CellParams::RHOVZ] / cell-> parameters[CellParams::RHO];
      } else {
	 averageVX = 0.0;
	 averageVY = 0.0;
	 averageVZ = 0.0;
      }
      Pressure = 0.0;
      return true;
   }
   
   
   // YK Adding pressure calculations to Vlasiator.
   // p_ij = m/3 * integral((v - <V>)_i(v - <V>)_j * f(r,v) dV)
   
   // Pressure tensor 6 components (11, 22, 33, 23, 13, 12) added by YK
   // Split into VariablePTensorDiagonal (11, 22, 33)
   // and VariablePTensorOffDiagonal (23, 13, 12)
   VariablePTensorDiagonal::VariablePTensorDiagonal(): DataReductionOperator() { }
   VariablePTensorDiagonal::~VariablePTensorDiagonal() { }

   std::string VariablePTensorDiagonal::getName() const {return "PTensorDiagonal";}

   bool VariablePTensorDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   bool VariablePTensorDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      const Real THIRD = 1.0/3.0;
      
      Real nvxvx_sum = 0.0;
      Real nvyvy_sum = 0.0;
      Real nvzvz_sum = 0.0;
      
      for(uint n=0; n<cell->number_of_blocks; n++) {
	 unsigned int blockId = cell->velocity_block_list[n];
	 const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
	 for (uint k=0; k<WID; ++k)
	    for (uint j=0; j<WID; ++j)
	       for (uint i=0; i<WID; ++i) {
		  const Real VX = block-> parameters[BlockParams::VXCRD] + (i+HALF) * block-> parameters[BlockParams::DVX];
		  const Real VY = block-> parameters[BlockParams::VYCRD] + (j+HALF) * block-> parameters[BlockParams::DVY];
		  const Real VZ = block-> parameters[BlockParams::VZCRD] + (k+HALF) * block-> parameters[BlockParams::DVZ];
		  
		  const Real DV3 = block-> parameters[BlockParams::DVX] * block-> parameters[BlockParams::DVY] * block-> parameters[BlockParams::DVZ];
		  
		  nvxvx_sum += block-> data[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX) * DV3;
		  nvyvy_sum += block-> data[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY) * DV3;
		  nvzvz_sum += block-> data[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ) * DV3;
	       }
      }
      
      // Accumulate contributions coming from this velocity block to the 
      // spatial cell velocity moments. If multithreading / OpenMP is used, 
      // these updates need to be atomic:
      
      PTensor[0] += physicalconstants::MASS_PROTON * THIRD * nvxvx_sum;
      PTensor[1] += physicalconstants::MASS_PROTON * THIRD * nvyvy_sum;
      PTensor[2] += physicalconstants::MASS_PROTON * THIRD * nvzvz_sum;
      
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }

   bool VariablePTensorDiagonal::setSpatialCell(const SpatialCell* cell) {
      if(cell-> parameters[CellParams::RHO] != 0.0) {
	 averageVX = cell-> parameters[CellParams::RHOVX] / cell-> parameters[CellParams::RHO];
	 averageVY = cell-> parameters[CellParams::RHOVY] / cell-> parameters[CellParams::RHO];
	 averageVZ = cell-> parameters[CellParams::RHOVZ] / cell-> parameters[CellParams::RHO];
      } else {
	 averageVX = 0.0;
	 averageVY = 0.0;
	 averageVZ = 0.0;
      }
      for(int i = 0; i < 3; i++) PTensor[i] = 0.0;
      return true;
   }
   
   VariablePTensorOffDiagonal::VariablePTensorOffDiagonal(): DataReductionOperator() { }
   VariablePTensorOffDiagonal::~VariablePTensorOffDiagonal() { }
   
   std::string VariablePTensorOffDiagonal::getName() const {return "PTensorOffDiagonal";}
   
   bool VariablePTensorOffDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   bool VariablePTensorOffDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      const Real THIRD = 1.0/3.0;
      
      Real nvxvy_sum = 0.0;
      Real nvzvx_sum = 0.0;
      Real nvyvz_sum = 0.0;
      
      for(uint n=0; n<cell->number_of_blocks; n++) {
	 unsigned int blockId = cell->velocity_block_list[n];
	 const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
	 for (uint k=0; k<WID; ++k)
	    for (uint j=0; j<WID; ++j)
	       for (uint i=0; i<WID; ++i) {
		  const Real VX = block-> parameters[BlockParams::VXCRD] + (i+HALF) * block-> parameters[BlockParams::DVX];
		  const Real VY = block-> parameters[BlockParams::VYCRD] + (j+HALF) * block-> parameters[BlockParams::DVY];
		  const Real VZ = block-> parameters[BlockParams::VZCRD] + (k+HALF) * block-> parameters[BlockParams::DVZ];
		  
		  const Real DV3 = block-> parameters[BlockParams::DVX] * block-> parameters[BlockParams::DVY] * block-> parameters[BlockParams::DVZ];
		  
		  nvxvy_sum += block-> data[cellIndex(i,j,k)] * (VX - averageVX) * (VY - averageVY) * DV3;
		  nvzvx_sum += block-> data[cellIndex(i,j,k)] * (VZ - averageVZ) * (VX - averageVX) * DV3;
		  nvyvz_sum += block-> data[cellIndex(i,j,k)] * (VY - averageVY) * (VZ - averageVZ) * DV3;
	       }
      }
      
      // Accumulate contributions coming from this velocity block to the 
      // spatial cell velocity moments. If multithreading / OpenMP is used, 
      // these updates need to be atomic:
      
      PTensor[0] += physicalconstants::MASS_PROTON * THIRD * nvyvz_sum;
      PTensor[1] += physicalconstants::MASS_PROTON * THIRD * nvzvx_sum;
      PTensor[2] += physicalconstants::MASS_PROTON * THIRD * nvxvy_sum;
      
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorOffDiagonal::setSpatialCell(const SpatialCell* cell) {
      if(cell-> parameters[CellParams::RHO] != 0.0) {
	 averageVX = cell-> parameters[CellParams::RHOVX] / cell-> parameters[CellParams::RHO];
	 averageVY = cell-> parameters[CellParams::RHOVY] / cell-> parameters[CellParams::RHO];
	 averageVZ = cell-> parameters[CellParams::RHOVZ] / cell-> parameters[CellParams::RHO];
      } else {
	 averageVX = 0.0;
	 averageVY = 0.0;
	 averageVZ = 0.0;
      }
      for(int i = 0; i < 3; i++) PTensor[i] = 0.0;
      return true;
   }



   //maximum absolute velocity in x,y, or z direction
   MaxVi::MaxVi(): DataReductionOperator() { }
   MaxVi::~MaxVi() { }
   
   std::string MaxVi::getName() const {return "MaximumVelocityComponent";}
   
   bool MaxVi::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
     dataType = "float";
     dataSize =  sizeof(Real);
     vectorSize = 1;
     return true;
   }
   

   bool MaxVi::reduceData(const SpatialCell* cell,Real* buffer) {
      const Real HALF = 0.5;
      const Real THIRD = 1.0/3.0;
     
      Real maxV=0;

      for(uint n=0; n<cell->number_of_blocks;n++) {
	 unsigned int blockId = cell->velocity_block_list[n];
	 const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
	 for (uint k=0; k<WID; ++k)
	    for (uint j=0; j<WID; ++j)
	       for (uint i=0; i<WID; ++i) {
                  const int celli=k*WID*WID+j*WID+i;
                  
		  const Real VX = block-> parameters[BlockParams::VXCRD] + (i+HALF) * block-> parameters[BlockParams::DVX];
		  const Real VY = block-> parameters[BlockParams::VYCRD] + (j+HALF) * block-> parameters[BlockParams::DVY];
		  const Real VZ = block-> parameters[BlockParams::VZCRD] + (k+HALF) * block-> parameters[BlockParams::DVZ];
		  
		  const Real DV3 = block-> parameters[BlockParams::DVX] * block-> parameters[BlockParams::DVY] * block-> parameters[BlockParams::DVZ];
		  
                  if(block->data[celli]!=0.0){
                     if(fabs(VX)>maxV) maxV=fabs(VX);
                     if(fabs(VY)>maxV) maxV=fabs(VY);
                     if(fabs(VZ)>maxV) maxV=fabs(VZ);
                  }
	       }
      }
      *buffer=maxV;
      return true;
   }
   bool MaxVi::reduceData(const SpatialCell* cell,char* buffer) {
      Real maxV;
      reduceData(cell,&maxV);
      const char* ptr = reinterpret_cast<const char*>(&maxV);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool MaxVi::setSpatialCell(const SpatialCell* cell) {
      return true;
   }

   // YK Integrated divergence of magnetic field
   // Integral of div B over the simulation volume =
   // Integral of flux of B on simulation volume surface
   DiagnosticFluxB::DiagnosticFluxB(): DataReductionOperator() { }
   DiagnosticFluxB::~DiagnosticFluxB() { }
   
   std::string DiagnosticFluxB::getName() const {return "FluxB";}
   
   bool DiagnosticFluxB::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   bool DiagnosticFluxB::reduceData(const SpatialCell* cell,Real * result) {
      creal x = cell->parameters[CellParams::XCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal y = cell->parameters[CellParams::YCRD];
      creal dy = cell->parameters[CellParams::DY];
      creal z = cell->parameters[CellParams::ZCRD];
      creal dz = cell->parameters[CellParams::DZ];
      creal cx = x + 0.5 * dx;
      creal cy = y + 0.5 * dy;
      creal cz = z + 0.5 * dz;
      
      Real value = 0.0;
      if(cx > Parameters::xmax - 2.0 * dx && cx < Parameters::xmax - dx) {
	 value += cell->parameters[CellParams::BX];
      } else if (cx < Parameters::xmin + 2.0 * dx && cx > Parameters::xmin + dx) {
	 value += -1.0*cell->parameters[CellParams::BX];
      }
      if(cy > Parameters::ymax - 2.0 * dy && cy < Parameters::ymax - dy) {
	 value += cell->parameters[CellParams::BY];
      } else if (cy < Parameters::ymin + 2.0 * dy && cy > Parameters::ymin + dy) {
	 value += -1.0*cell->parameters[CellParams::BY];
      }
      if(cz > Parameters::zmax - 2.0 * dz && cz < Parameters::zmax - dz) {
	 value += cell->parameters[CellParams::BZ];
      } else if (cz < Parameters::zmin + 2.0 * dz && cz > Parameters::zmin + dz) {
	 value += -1.0*cell->parameters[CellParams::BZ];
      }
      *result = value;
      
      return true;
   }
   
   bool DiagnosticFluxB::setSpatialCell(const SpatialCell* cell) {return true;}
   
   VariabledBxdz::VariabledBxdz(): DataReductionOperator() { }
   VariabledBxdz::~VariabledBxdz() { }
   
   bool VariabledBxdz::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize = sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   std::string VariabledBxdz::getName() const {return "dBxdz";}
   
   bool VariabledBxdz::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&value);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariabledBxdz::setSpatialCell(const SpatialCell* cell) {
      value = cell->derivatives[fieldsolver::dBxdz];
      return true;
   }
   
   
} // namespace DRO
