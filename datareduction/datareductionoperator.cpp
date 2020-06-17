/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <limits>
#include <array>
#include "datareductionoperator.h"
#include "../object_wrapper.h"

using namespace std;

typedef Parameters P;

namespace DRO {
   
   // ************************************************************
   // ***** DEFINITIONS FOR DATAREDUCTIONOPERATOR BASE CLASS *****
   // ************************************************************
   
   /** DataReductionOperator base class constructor. The constructor is empty.*/
   DataReductionOperator::DataReductionOperator() { }
   
   /** DataReductionOperator base class virtual destructor. The destructor is empty.*/
   DataReductionOperator::~DataReductionOperator() { }

   /** Reduce the data and write the data vector to the given vlsv output buffer.
    * @param cell the SpatialCell to reduce data out of
    * @param buffer Buffer in which the reduced data is written.
    * @return If true, DataReductionOperator reduced data successfully.
    */
   bool DataReductionOperator::reduceData(const SpatialCell* cell,char* buffer) {
      cerr << "ERROR: DataReductionOperator::reduceData called instead of derived class function! (variable" <<
               getName() << ")" << endl;
      cerr << "       Did you use a diagnostic reducer for writing bulk data?" << endl;
      return false;
   }
   
   /** Reduce the data and write the data vector to the given variable.
    * If the vector length is larger than one, memory gets corrupted.
    * Note that this function is only used for writing into diagnostic files.
    * @param cell the SpatialCell to reduce data out of
    * @param buffer Buffer in which the reduced data is written.
    * @return If true, DataReductionOperator reduced data successfully.
    */
   bool DataReductionOperator::reduceDiagnostic(const SpatialCell* cell,Real* result) {
      cerr << "ERROR: DataReductionOperator::reduceData called instead of derived class function! (variable " <<
              getName() << ")" << endl;
      cerr << "       Did you use a bulk reducer for writing diagnostic data?" << endl;
      return false;
   }

   DataReductionOperatorCellParams::DataReductionOperatorCellParams(const std::string& name,const unsigned int parameterIndex,const unsigned int _vectorSize):
   DataReductionOperator() {
      vectorSize=_vectorSize;
      variableName=name;
      _parameterIndex=parameterIndex;
   }
   DataReductionOperatorCellParams::~DataReductionOperatorCellParams() { }
   
   bool DataReductionOperatorCellParams::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& _vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      _vectorSize = vectorSize;
      return true;
   }

   std::string DataReductionOperatorCellParams::getName() const {return variableName;}
   
   bool DataReductionOperatorCellParams::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(data);
      for (uint i = 0; i < vectorSize*sizeof(Real); ++i){
         buffer[i] = ptr[i];
      }
      return true;
   }
   
   bool DataReductionOperatorCellParams::reduceDiagnostic(const SpatialCell* cell,Real* buffer){
      //If vectorSize is >1 it still works, we just give the first value and no other ones..
      *buffer=data[0];
      return true;
   }
   bool DataReductionOperatorCellParams::setSpatialCell(const SpatialCell* cell) {
      for (uint i=0; i<vectorSize; i++) {
         if(std::isinf(cell->parameters[_parameterIndex+i]) || std::isnan(cell->parameters[_parameterIndex+i])) {
            string message = "The DataReductionOperator " + this->getName() + " returned a nan or an inf in its " + std::to_string(i) + "-component.";
            bailout(true, message, __FILE__, __LINE__);
         }
      }
      data  = &(cell->parameters[_parameterIndex]);
      return true;
   }

   std::string DataReductionOperatorFsGrid::getName() const {return variableName;}
   bool DataReductionOperatorFsGrid::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      // These are only set to dmmy values, as this reducer writes its own vlsv dataset anyway
      dataType = "float";
      dataSize = sizeof(double);
      vectorSize = 1;
      return true;
   }
   bool DataReductionOperatorFsGrid::reduceData(const SpatialCell* cell,char* buffer) {
      // This returns false, since it will handle writing itself in writeFsGridData below.
      return false;
   }
   bool DataReductionOperatorFsGrid::reduceDiagnostic(const SpatialCell* cell,Real * result) {
      return false;
   }
   bool DataReductionOperatorFsGrid::setSpatialCell(const SpatialCell* cell) {
      return true;
   }
   
   bool DataReductionOperatorFsGrid::writeFsGridData(
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid,
                      const std::string& meshName, vlsv::Writer& vlsvWriter,
                      const bool writeAsFloat) {

      std::map<std::string,std::string> attribs;
      attribs["mesh"]=meshName;
      attribs["name"]=variableName;
      attribs["unit"]=unit;
      attribs["unitLaTeX"]=unitLaTeX;
      attribs["unitConversion"]=unitConversion;
      attribs["variableLaTeX"]=variableLaTeX;

      std::vector<double> varBuffer =
         lambda(perBGrid,EGrid,EHallGrid,EGradPeGrid,momentsGrid,dPerBGrid,dMomentsGrid,BgBGrid,volGrid,technicalGrid);

      std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
      int vectorSize = varBuffer.size() / (gridSize[0]*gridSize[1]*gridSize[2]);

      if(writeAsFloat) {
         // Convert down to 32bit floats to save output space
         std::vector<float> varBufferFloat(varBuffer.size());
         for(int i=0; i<varBuffer.size(); i++) {
            varBufferFloat[i] = (float)varBuffer[i];
         }
         if(vlsvWriter.writeArray("VARIABLE",attribs, "float", gridSize[0]*gridSize[1]*gridSize[2], vectorSize, sizeof(float), reinterpret_cast<const char*>(varBufferFloat.data())) == false) {
            string message = "The DataReductionOperator " + this->getName() + " failed to write its data.";
            bailout(true, message, __FILE__, __LINE__);
         }

      } else {
         if(vlsvWriter.writeArray("VARIABLE",attribs, "float", gridSize[0]*gridSize[1]*gridSize[2], vectorSize, sizeof(double), reinterpret_cast<const char*>(varBuffer.data())) == false) {
            string message = "The DataReductionOperator " + this->getName() + " failed to write its data.";
            bailout(true, message, __FILE__, __LINE__);
         }
      }

      return true;
   }

   DataReductionOperatorBVOLDerivatives::DataReductionOperatorBVOLDerivatives(const std::string& name,const unsigned int parameterIndex,const unsigned int vectorSize):
   DataReductionOperatorCellParams(name,parameterIndex,vectorSize) {
      
   }
   //a version with derivatives, this is the only function that is different
   bool DataReductionOperatorBVOLDerivatives::setSpatialCell(const SpatialCell* cell) {
      data  = &(cell->derivativesBVOL[_parameterIndex]);
      return true;
   }
   
   
   
   //------------------ total BVOL --------------------------------------- 
   VariableBVol::VariableBVol(): DataReductionOperator() { }
   VariableBVol::~VariableBVol() { }
   
   bool VariableBVol::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   std::string VariableBVol::getName() const {return "vg_b_vol";}
   
   bool VariableBVol::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(B);
      for (uint i = 0; i < 3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableBVol::setSpatialCell(const SpatialCell* cell) {
      B[0] = cell->parameters[CellParams::PERBXVOL] +  cell->parameters[CellParams::BGBXVOL];
      B[1] = cell->parameters[CellParams::PERBYVOL] +  cell->parameters[CellParams::BGBYVOL];
      B[2] = cell->parameters[CellParams::PERBZVOL] +  cell->parameters[CellParams::BGBZVOL];
      if(std::isinf(B[0]) || std::isnan(B[0]) ||
         std::isinf(B[1]) || std::isnan(B[1]) ||
         std::isinf(B[2]) || std::isnan(B[2])
      ) {
         string message = "The DataReductionOperator " + this->getName() + " returned a nan or an inf.";
         bailout(true, message, __FILE__, __LINE__);
      }
      return true;
   }

   //MPI rank
   MPIrank::MPIrank(): DataReductionOperator() { }
   MPIrank::~MPIrank() { }
   
   bool MPIrank::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = sizeof(int);
      vectorSize = 1;
      return true;
   }
   
   std::string MPIrank::getName() const {return "vg_rank";}
   
   bool MPIrank::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&mpiRank);
      for (uint i = 0; i < sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool MPIrank::setSpatialCell(const SpatialCell* cell) {
      int intRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&intRank);
      rank = 1.0*intRank;
      mpiRank = intRank;
      return true;
   }
   
   // BoundaryType
   BoundaryType::BoundaryType(): DataReductionOperator() { }
   BoundaryType::~BoundaryType() { }
   
   bool BoundaryType::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = sizeof(int);
      vectorSize = 1;
      return true;
   }
   
   std::string BoundaryType::getName() const {return "vg_boundarytype";}
   
   bool BoundaryType::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&boundaryType);
      for (uint i = 0; i < sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool BoundaryType::setSpatialCell(const SpatialCell* cell) {
      boundaryType = (int)cell->sysBoundaryFlag;
      return true;
   }


      // BoundaryLayer
   BoundaryLayer::BoundaryLayer(): DataReductionOperator() { }
   BoundaryLayer::~BoundaryLayer() { }
   
   bool BoundaryLayer::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = sizeof(int);
      vectorSize = 1;
      return true;
   }
   
   std::string BoundaryLayer::getName() const {return "vg_boundarylayer";}
   
   bool BoundaryLayer::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&boundaryLayer);
      for (uint i = 0; i < sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool BoundaryLayer::setSpatialCell(const SpatialCell* cell) {
      boundaryLayer = (int)cell->sysBoundaryLayer;
      return true;
   }
   
   // Blocks
   Blocks::Blocks(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName=getObjectWrapper().particleSpecies[popID].name;
   }
   Blocks::~Blocks() { }
   
   bool Blocks::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "uint";
      dataSize = sizeof(int);
      vectorSize = 1;
      return true;
   }
   
   std::string Blocks::getName() const {return popName + "/vg_blocks";}
   
   bool Blocks::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&nBlocks);
      for (uint i = 0; i < sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool Blocks::reduceDiagnostic(const SpatialCell* cell,Real* buffer) {
      *buffer = 1.0 * nBlocks;
      return true;
   }
  
   bool Blocks::setSpatialCell(const SpatialCell* cell) {
      nBlocks = cell->get_number_of_velocity_blocks(popID);
      return true;
   }
   
   // Scalar pressure from the stored values which were calculated to be used by the solvers
   VariablePressureSolver::VariablePressureSolver(): DataReductionOperator() { }
   VariablePressureSolver::~VariablePressureSolver() { }
   
   std::string VariablePressureSolver::getName() const {return "vg_pressure";}
   
   bool VariablePressureSolver::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   bool VariablePressureSolver::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&Pressure);
      for (uint i = 0; i < sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePressureSolver::setSpatialCell(const SpatialCell* cell) {
      Pressure = 1.0/3.0 * (cell->parameters[CellParams::P_11] + cell->parameters[CellParams::P_22] + cell->parameters[CellParams::P_33]);
      return true;
   }
   
   // YK Adding pressure calculations to Vlasiator.
   // p_ij = m/3 * integral((v - <V>)_i(v - <V>)_j * f(r,v) dV)
   
   // Pressure tensor 6 components (11, 22, 33, 23, 13, 12) added by YK
   // Split into VariablePTensorDiagonal (11, 22, 33)
   // and VariablePTensorOffDiagonal (23, 13, 12)
   VariablePTensorDiagonal::VariablePTensorDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
   }
   VariablePTensorDiagonal::~VariablePTensorDiagonal() { }
   
   std::string VariablePTensorDiagonal::getName() const {return popName + "/vg_ptensor_diagonal";}
   
   bool VariablePTensorDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   bool VariablePTensorDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      # pragma omp parallel
      {
         Real thread_nvxvx_sum = 0.0;
         Real thread_nvyvy_sum = 0.0;
         Real thread_nvzvz_sum = 0.0;
         
         const Real* parameters  = cell->get_block_parameters(popID);
         const Realf* block_data = cell->get_data(popID);
         
         # pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); n++) {
	    for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
	       const Real VX 
		 =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
		 + (i + HALF)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	       const Real VY 
		 =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
		 + (j + HALF)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
	       const Real VZ 
		 =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] 
		 + (k + HALF)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
	       const Real DV3 
		 = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]
		 * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY] 
		 * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
                     
	       thread_nvxvx_sum += block_data[n * SIZE_VELBLOCK+cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX) * DV3;
	       thread_nvyvy_sum += block_data[n * SIZE_VELBLOCK+cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY) * DV3;
	       thread_nvzvz_sum += block_data[n * SIZE_VELBLOCK+cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ) * DV3;
            }
         }
         thread_nvxvx_sum *= getObjectWrapper().particleSpecies[popID].mass;
         thread_nvyvy_sum *= getObjectWrapper().particleSpecies[popID].mass;
         thread_nvzvz_sum *= getObjectWrapper().particleSpecies[popID].mass;

         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            PTensor[0] += thread_nvxvx_sum;
            PTensor[1] += thread_nvyvy_sum;
            PTensor[2] += thread_nvzvz_sum;
         }
      }
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < 3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorDiagonal::setSpatialCell(const SpatialCell* cell) {
      averageVX = cell-> parameters[CellParams::VX];
      averageVY = cell-> parameters[CellParams::VY];
      averageVZ = cell-> parameters[CellParams::VZ];
      for(int i = 0; i < 3; i++) PTensor[i] = 0.0;
      return true;
   }
   
   VariablePTensorOffDiagonal::VariablePTensorOffDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
   }
   VariablePTensorOffDiagonal::~VariablePTensorOffDiagonal() { }
   
   std::string VariablePTensorOffDiagonal::getName() const {return popName + "/vg_ptensor_offdiagonal";}
   
   bool VariablePTensorOffDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   bool VariablePTensorOffDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      # pragma omp parallel
      {
         Real thread_nvxvy_sum = 0.0;
         Real thread_nvzvx_sum = 0.0;
         Real thread_nvyvz_sum = 0.0;
         
         const Real* parameters = cell->get_block_parameters(popID);
         const Realf* block_data = cell->get_data(popID);
         
         # pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); n++) {               
	    for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
	       const Real VX 
		 =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
		 + (i + HALF)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	       const Real VY 
		 =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
		 + (j + HALF)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
	       const Real VZ 
		 =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] 
		 + (k + HALF)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
	       const Real DV3 
		 = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]
		 * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY] 
		 * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
	       
	       thread_nvxvy_sum += block_data[n * SIZE_VELBLOCK+cellIndex(i,j,k)] * (VX - averageVX) * (VY - averageVY) * DV3;
	       thread_nvzvx_sum += block_data[n * SIZE_VELBLOCK+cellIndex(i,j,k)] * (VZ - averageVZ) * (VX - averageVX) * DV3;
	       thread_nvyvz_sum += block_data[n * SIZE_VELBLOCK+cellIndex(i,j,k)] * (VY - averageVY) * (VZ - averageVZ) * DV3;
            }
         }
         thread_nvxvy_sum *= getObjectWrapper().particleSpecies[popID].mass;
         thread_nvzvx_sum *= getObjectWrapper().particleSpecies[popID].mass;
         thread_nvyvz_sum *= getObjectWrapper().particleSpecies[popID].mass;
         
         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            PTensor[0] += thread_nvyvz_sum;
            PTensor[1] += thread_nvzvx_sum;
            PTensor[2] += thread_nvxvy_sum;
         }
      }
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < 3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorOffDiagonal::setSpatialCell(const SpatialCell* cell) {
      averageVX = cell-> parameters[CellParams::VX];
      averageVY = cell-> parameters[CellParams::VY];
      averageVZ = cell-> parameters[CellParams::VZ];
      for(int i = 0; i < 3; i++) PTensor[i] = 0.0;
      return true;
   }   
   
   // YK maximum value of the distribution function (diagnostic)
   MaxDistributionFunction::MaxDistributionFunction(cuint _popID): DataReductionOperator(),popID(_popID) {
     popName=getObjectWrapper().particleSpecies[popID].name;
   }
   MaxDistributionFunction::~MaxDistributionFunction() { }
   
   std::string MaxDistributionFunction::getName() const {return popName + "/maximumdistributionfunctionvalue";}
   
   bool MaxDistributionFunction::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }   
   
   bool MaxDistributionFunction::reduceDiagnostic(const SpatialCell* cell,Real* buffer) {
      maxF = std::numeric_limits<Real>::min();
      
      #pragma omp parallel 
      {
         Real threadMax = std::numeric_limits<Real>::min();
         
         const Realf* block_data = cell->get_data(popID);
         
         #pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); ++n) {
	    for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
	       threadMax = max((Real)(block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)]), threadMax);
            }
         }

         #pragma omp critical
         {
            maxF = max(threadMax, maxF);
         }
      }
      
      *buffer = maxF;
      return true;
   }
   
   bool MaxDistributionFunction::reduceData(const SpatialCell* cell,char* buffer) {
      Real dummy;
      reduceDiagnostic(cell,&dummy);
      const char* ptr = reinterpret_cast<const char*>(&dummy);
      for (uint i = 0; i < sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool MaxDistributionFunction::setSpatialCell(const SpatialCell* cell) {
      return true;
   }
   
   
   // YK minimum value of the distribution function (diagnostic)
   MinDistributionFunction::MinDistributionFunction(cuint _popID): DataReductionOperator(),popID(_popID) {
     popName=getObjectWrapper().particleSpecies[popID].name;
   }
   MinDistributionFunction::~MinDistributionFunction() { }
   
   std::string MinDistributionFunction::getName() const {return popName + "/minimumdistributionfunctionvalue";}
   
   bool MinDistributionFunction::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }   
   
   bool MinDistributionFunction::reduceDiagnostic(const SpatialCell* cell,Real* buffer) {
      minF =  std::numeric_limits<Real>::max();

      #pragma omp parallel 
      {
         Real threadMin = std::numeric_limits<Real>::max();
         
         const Realf* block_data = cell->get_data(popID);

         #pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); ++n) {
	    for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
	       threadMin = min((Real)(block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)]), threadMin);
            }
         }
         
         #pragma omp critical
         {
            minF = min(threadMin, minF);
         }
      }
      
      *buffer = minF;
      return true;
   }
   
   bool MinDistributionFunction::reduceData(const SpatialCell* cell,char* buffer) {
      Real dummy;
      reduceDiagnostic(cell,&dummy);
      const char* ptr = reinterpret_cast<const char*>(&dummy);
      for (uint i = 0; i < sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool MinDistributionFunction::setSpatialCell(const SpatialCell* cell) {
      return true;
   }

  /*******
	  Helper functions for finding the velocity cell indices or IDs within a single velocity block
	  either belonging to the thermal or the non-thermal population. 
	  There is some code duplication here, but as these helper functions are called within threads for
	  block separately, it's preferable to have them fast even at the cost of code repetition.
  ********/

   //Helper function for getting the velocity cell ids that are a part of the nonthermal population:
   static void getNonthermalVelocityCells(
      const Real* block_parameters,
      vector<uint64_t> & vCellIds,
      cuint popID
   ) {
      creal HALF = 0.5;
      const std::array<Real, 3> thermalV = getObjectWrapper().particleSpecies[popID].thermalV;
      creal thermalRadius = getObjectWrapper().particleSpecies[popID].thermalRadius;
      // Go through every velocity cell (i, j, k are indices)
      for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
         // Get the vx, vy, vz coordinates of the velocity cell
         const Real VX = block_parameters[BlockParams::VXCRD] + (i + HALF) * block_parameters[BlockParams::DVX];
         const Real VY = block_parameters[BlockParams::VYCRD] + (j + HALF) * block_parameters[BlockParams::DVY];
         const Real VZ = block_parameters[BlockParams::VZCRD] + (k + HALF) * block_parameters[BlockParams::DVZ];
         // Compare the distance of the velocity cell from the center of the maxwellian distribution to the radius of the maxwellian distribution
         if( ( (thermalV[0] - VX) * (thermalV[0] - VX)
             + (thermalV[1] - VY) * (thermalV[1] - VY)
             + (thermalV[2] - VZ) * (thermalV[2] - VZ) )
             >
             thermalRadius*thermalRadius ) {
             //The velocity cell is a part of the nonthermal population:
             vCellIds.push_back(cellIndex(i,j,k));
          }
      }
   }
   //Helper function for getting the velocity cell ids that are a part of the nonthermal population:
   static void getThermalVelocityCells(
      const Real* block_parameters,
      vector<uint64_t> & vCellIds,
      cuint popID
   ) {
      creal HALF = 0.5;
      const std::array<Real, 3> thermalV = getObjectWrapper().particleSpecies[popID].thermalV;
      creal thermalRadius = getObjectWrapper().particleSpecies[popID].thermalRadius;
      for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
         const Real VX = block_parameters[BlockParams::VXCRD] + (i + HALF) * block_parameters[BlockParams::DVX];
         const Real VY = block_parameters[BlockParams::VYCRD] + (j + HALF) * block_parameters[BlockParams::DVY];
         const Real VZ = block_parameters[BlockParams::VZCRD] + (k + HALF) * block_parameters[BlockParams::DVZ];
         if( ( (thermalV[0] - VX) * (thermalV[0] - VX)
             + (thermalV[1] - VY) * (thermalV[1] - VY)
             + (thermalV[2] - VZ) * (thermalV[2] - VZ) )
             <=
             thermalRadius*thermalRadius ) {
             //The velocity cell is not a part of the nonthermal population:
             vCellIds.push_back(cellIndex(i,j,k));
          }
      }
   }
   //Helper function for getting the velocity cell indices that are a part of the nonthermal population:
   static void getNonthermalVelocityCellIndices(
      const Real* block_parameters,
      vector<array<uint, 3>> & vCellIndices,
      cuint popID
   ) {
      creal HALF = 0.5;
      const std::array<Real, 3> thermalV = getObjectWrapper().particleSpecies[popID].thermalV;
      creal thermalRadius = getObjectWrapper().particleSpecies[popID].thermalRadius;
      // Go through a block's every velocity cell
      for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
         // Get the coordinates of the velocity cell (e.g. VX = block_vx_min_coordinates + (velocity_cell_indice_x+0.5)*length_of_velocity_cell_in_x_direction
         const Real VX = block_parameters[BlockParams::VXCRD] + (i + HALF) * block_parameters[BlockParams::DVX];
         const Real VY = block_parameters[BlockParams::VYCRD] + (j + HALF) * block_parameters[BlockParams::DVY];
         const Real VZ = block_parameters[BlockParams::VZCRD] + (k + HALF) * block_parameters[BlockParams::DVZ];
         // Calculate the distance of the velocity cell from the center of the maxwellian distribution and compare it to the approximate radius of the maxwellian distribution
         if( ( (thermalV[0] - VX) * (thermalV[0] - VX)
             + (thermalV[1] - VY) * (thermalV[1] - VY)
             + (thermalV[2] - VZ) * (thermalV[2] - VZ) )
             >
             thermalRadius*thermalRadius ) {
             //The velocity cell is a part of the nonthermal population because it is not within the radius:
             const array<uint, 3> indices{{i, j, k}};
             vCellIndices.push_back( indices );
          }
      }
   }
   //Helper function for getting the velocity cell indices that are not a part of the nonthermal population:
   static void getThermalVelocityCellIndices(
      const Real* block_parameters,
      vector<array<uint, 3>> & vCellIndices,
      cuint popID
   ) {
      creal HALF = 0.5;
      const std::array<Real, 3> thermalV = getObjectWrapper().particleSpecies[popID].thermalV;
      creal thermalRadius = getObjectWrapper().particleSpecies[popID].thermalRadius;
      // Go through a block's every velocity cell
      for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
         // Get the coordinates of the velocity cell (e.g. VX = block_vx_min_coordinates + (velocity_cell_indice_x+0.5)*length_of_velocity_cell_in_x_direction
         const Real VX = block_parameters[BlockParams::VXCRD] + (i + HALF) * block_parameters[BlockParams::DVX];
         const Real VY = block_parameters[BlockParams::VYCRD] + (j + HALF) * block_parameters[BlockParams::DVY];
         const Real VZ = block_parameters[BlockParams::VZCRD] + (k + HALF) * block_parameters[BlockParams::DVZ];
         // Calculate the distance of the velocity cell from the center of the maxwellian distribution and compare it to the approximate radius of the maxwellian distribution
         if( ( (thermalV[0] - VX) * (thermalV[0] - VX)
             + (thermalV[1] - VY) * (thermalV[1] - VY)
             + (thermalV[2] - VZ) * (thermalV[2] - VZ) )
             <=
             thermalRadius*thermalRadius ) {
             //The velocity cell is part of the thermal population because it is within the radius:
             const array<uint, 3> indices{{i, j, k}};
             vCellIndices.push_back( indices );
          }
      }
   }

  /********
	   Next level of helper functions - these include threading and calculate zeroth or first velocity moments or the
	   diagonal / off-diagonal pressure tensor components for
	   thermal or non-thermal populations  ********/

   //Calculates rho thermal or rho non-thermal
   static void rhoNonthermalCalculation( const SpatialCell * cell, const bool calculateNonthermal, cuint popID, Real & rho ) {
      const Real HALF = 0.5;
      # pragma omp parallel
      {
         Real thread_n_sum = 0.0;
         
         const Real* parameters = cell->get_block_parameters(popID);
         const Realf* block_data = cell->get_data(popID);
         
         # pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); ++n) {
            const Real DV3
            = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]
            * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]
            * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
            vector< uint64_t > vCells; //Velocity cell ids
            vCells.clear();
            if ( calculateNonthermal == true ) {
               getNonthermalVelocityCells(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCells, popID);
            } else {
               getThermalVelocityCells(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCells, popID);
            }
            for( vector< uint64_t >::const_iterator it = vCells.begin(); it != vCells.end(); ++it ) {
               //velocity cell id = *it
               thread_n_sum += block_data[n * SIZE_VELBLOCK + (*it)] * DV3;
            }
         }
         // Accumulate contributions coming from this velocity block
         // If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         // todo: use omp reduction
         # pragma omp critical
         {
            rho += thread_n_sum;
         }
      }
      return;
   }

   static void VNonthermalCalculation( const SpatialCell * cell, const bool calculateNonthermal, cuint popID, Real * V ) {
      const Real HALF = 0.5;
      // Make sure the V is initialized
      V[0] = 0;
      V[1] = 0;
      V[2] = 0;
      Real n_sum = 0;
      # pragma omp parallel
      {
         Real thread_nvx_sum = 0.0;
         Real thread_nvy_sum = 0.0;
         Real thread_nvz_sum = 0.0;
         Real thread_n_sum = 0.0;

         const Real* parameters = cell->get_block_parameters(popID);
         const Realf* block_data = cell->get_data(popID);
         
         # pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); ++n) {
            // Get the volume of a velocity cell
            const Real DV3
            = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]
            * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]
            * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
            // Get the velocity cell indices of the cells that are a part of the nonthermal population
            vector< array<uint, 3> > vCellIndices;
            vCellIndices.clear();
            // Save indices to the std::vector
            if( calculateNonthermal == true ) {
               getNonthermalVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            } else {
               getThermalVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            }
            // We have now fetched all of the needed velocity cell indices, so now go through them:
            for( vector< array<uint, 3> >::const_iterator it = vCellIndices.begin(); it != vCellIndices.end(); ++it ) {
               // Get the indices of the current iterated velocity cell
               const array<uint, 3> indices = *it;
               const uint i = indices[0];
               const uint j = indices[1];
               const uint k = indices[2];
               // Get the coordinates of the velocity cell (e.g. VX = block_vx_min_coordinates + (velocity_cell_indice_x+0.5)*length_of_velocity_cell_in_x_direction)
               const Real VX = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] + (i + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
               const Real VY = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] + (j + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
               const Real VZ = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] + (k + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
               // Add the value of the coordinates and multiply by the AVGS value of the velocity cell and the volume of the velocity cell
               thread_nvx_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)]*VX*DV3;
               thread_nvy_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)]*VY*DV3;
               thread_nvz_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)]*VZ*DV3;
               thread_n_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)]*DV3;
            }
         } // for-loop over velocity blocks

         // Accumulate contributions coming from this velocity block.
         // If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            V[0] += thread_nvx_sum;
            V[1] += thread_nvy_sum;
            V[2] += thread_nvz_sum;
            n_sum += thread_n_sum;
         }
      }

      // Finally, divide n_sum*V by V.
      V[0]/=n_sum;
      V[1]/=n_sum;
      V[2]/=n_sum;
      return;
   }

   static void PTensorDiagonalNonthermalCalculations( const SpatialCell * cell,
                                                      const bool calculateNonthermal,
                                                      const Real averageVX,
                                                      const Real averageVY,
                                                      const Real averageVZ,
                                                      cuint popID,
                                                      Real * PTensor ) {
      const Real HALF = 0.5;
      # pragma omp parallel
      {
         Real thread_nvxvx_sum = 0.0;
         Real thread_nvyvy_sum = 0.0;
         Real thread_nvzvz_sum = 0.0;
         
         const Real* parameters = cell->get_block_parameters(popID);
         const Realf* block_data = cell->get_data(popID);
      
         # pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); ++n) {
            const Real DV3
            = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]
            * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]
            * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
            vector< array<uint, 3> > vCellIndices;
            vCellIndices.clear();
            if( calculateNonthermal == true ) {
               getNonthermalVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            } else {
               getThermalVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            }
            for( vector< array<uint, 3> >::const_iterator it = vCellIndices.begin(); it != vCellIndices.end(); ++it ) {
               //Go through every velocity cell:
               const array<uint, 3> indices = *it;
               const uint i = indices[0];
               const uint j = indices[1];
               const uint k = indices[2];
               const Real VX = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] + (i + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
               const Real VY = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] + (j + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
               const Real VZ = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] + (k + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
               thread_nvxvx_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX) * DV3;
               thread_nvyvy_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY) * DV3;
               thread_nvzvz_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ) * DV3;
            }
         }
         thread_nvxvx_sum *= getObjectWrapper().particleSpecies[popID].mass;
         thread_nvyvy_sum *= getObjectWrapper().particleSpecies[popID].mass;
         thread_nvzvz_sum *= getObjectWrapper().particleSpecies[popID].mass;

         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            PTensor[0] += thread_nvxvx_sum;
            PTensor[1] += thread_nvyvy_sum;
            PTensor[2] += thread_nvzvz_sum;
         }
      }
      return;
   }

   static void PTensorOffDiagonalNonthermalCalculations( const SpatialCell * cell,
                                                         const bool calculateNonthermal,
                                                         const Real averageVX,
                                                         const Real averageVY,
                                                         const Real averageVZ,
                                                         cuint popID,
                                                         Real * PTensor ) {
      const Real HALF = 0.5;
      # pragma omp parallel
      {
         Real thread_nvxvy_sum = 0.0;
         Real thread_nvzvx_sum = 0.0;
         Real thread_nvyvz_sum = 0.0;
         
         const Real* parameters = cell->get_block_parameters(popID);
         const Realf* block_data = cell->get_data(popID);
      
         # pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); ++n) {
            const Real DV3
            = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]
            * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]
            * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
            vector< array<uint, 3> > vCellIndices;
            if( calculateNonthermal == true ) {
               getNonthermalVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            } else {
               getThermalVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            }
            for( vector< array<uint, 3> >::const_iterator it = vCellIndices.begin(); it != vCellIndices.end(); ++it ) {
               //Go through every velocity cell:
               const array<uint, 3> indices = *it;
               const uint i = indices[0];
               const uint j = indices[1];
               const uint k = indices[2];
               const Real VX = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] + (i + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
               const Real VY = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] + (j + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
               const Real VZ = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] + (k + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
               thread_nvxvy_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)] * (VX - averageVX) * (VY - averageVY) * DV3;
               thread_nvzvx_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)] * (VZ - averageVZ) * (VX - averageVX) * DV3;
               thread_nvyvz_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)] * (VY - averageVY) * (VZ - averageVZ) * DV3;
            }
         }
         thread_nvxvy_sum *= getObjectWrapper().particleSpecies[popID].mass;
         thread_nvzvx_sum *= getObjectWrapper().particleSpecies[popID].mass;
         thread_nvyvz_sum *= getObjectWrapper().particleSpecies[popID].mass;
         
         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            PTensor[0] += thread_nvyvz_sum;
            PTensor[1] += thread_nvzvx_sum;
            PTensor[2] += thread_nvxvy_sum;
         }
      }
   }

  /********* 
	     End velocity moment / thermal/non-thermal helper functions
  *********/


   VariableMeshData::VariableMeshData(): DataReductionOperatorHandlesWriting() { }
   VariableMeshData::~VariableMeshData() { }
   
   std::string VariableMeshData::getName() const {return "vg_meshdata";}
   
   bool VariableMeshData::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      return true;
   }
   
   bool VariableMeshData::setSpatialCell(const SpatialCell* cell) {return true;}
   
   bool VariableMeshData::writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                    const std::vector<CellID>& cells,const std::string& meshName,
                                    vlsv::Writer& vlsvWriter) {
      bool success = true;
      for (size_t i = 0; i < getObjectWrapper().meshData.size(); ++i) {
         const string dataName = getObjectWrapper().meshData.getName(i);
         
         // If dataName equals "" then something is wrong, skip array
         if (dataName.size() == 0) continue;
         
         size_t dataSize = getObjectWrapper().meshData.getDataSize(i);
         const std::string dataType = getObjectWrapper().meshData.getDataType(i);
         size_t vectorSize = getObjectWrapper().meshData.getVectorSize(i);
         size_t arraySize = getObjectWrapper().meshData.getMeshSize();
         char* pointer = getObjectWrapper().meshData.getData<char>(i);

         if (vectorSize == 0 || vectorSize > 3) continue;
         
         map<string,string> attribs;
         attribs["mesh"] = meshName;
         attribs["name"] = dataName;
         
         if (vlsvWriter.writeArray("VARIABLE",attribs,dataType,arraySize,vectorSize,dataSize,pointer) == false) {
            cerr << "write failed!" << endl;
            success = false;
         }
      }
      return success;
   }
   
   // Rho nonthermal:
   VariableRhoNonthermal::VariableRhoNonthermal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].thermalRadius == 0.0) ? true : false;
   }
   VariableRhoNonthermal::~VariableRhoNonthermal() { }
   
   std::string VariableRhoNonthermal::getName() const {return popName + "/vg_rho_nonthermal";}
   
   bool VariableRhoNonthermal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 1;
      return true;
   }
   
   bool VariableRhoNonthermal::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateNonthermal = true;
      rhoNonthermalCalculation( cell, calculateNonthermal, popID, RhoNonthermal );
      const char* ptr = reinterpret_cast<const char*>(&RhoNonthermal);
      for (uint i = 0; i < sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoNonthermal::setSpatialCell(const SpatialCell* cell) {
      RhoNonthermal = 0.0;
      return true;
   }

   // Rho thermal:
   VariableRhoThermal::VariableRhoThermal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].thermalRadius == 0.0) ? true : false;
   }
   VariableRhoThermal::~VariableRhoThermal() { }
   
   std::string VariableRhoThermal::getName() const {return popName + "/vg_rho_thermal";}
   
   bool VariableRhoThermal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 1;
      return true;
   }
   
   bool VariableRhoThermal::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateNonthermal = false; //We don't want nonthermal
      rhoNonthermalCalculation( cell, calculateNonthermal, popID, RhoThermal );
      const char* ptr = reinterpret_cast<const char*>(&RhoThermal);
      for (uint i = 0; i < sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoThermal::setSpatialCell(const SpatialCell* cell) {
      RhoThermal = 0.0;
      return true;
   }

   // v nonthermal:
   VariableVNonthermal::VariableVNonthermal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].thermalRadius == 0.0) ? true : false;
   }
   VariableVNonthermal::~VariableVNonthermal() { }
   
   std::string VariableVNonthermal::getName() const {return popName + "/vg_v_nonthermal";}
   
   bool VariableVNonthermal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }

   bool VariableVNonthermal::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateNonthermal = true;
      //Calculate v nonthermal
      VNonthermalCalculation( cell, calculateNonthermal, popID, VNonthermal );
      const uint VNonthermalSize = 3;
      const char* ptr = reinterpret_cast<const char*>(&VNonthermal);
      for (uint i = 0; i < VNonthermalSize*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableVNonthermal::setSpatialCell(const SpatialCell* cell) {
      // Initialize values
      for( uint i = 0; i < 3; ++i ) {
         VNonthermal[i] = 0.0;
      }
      return true;
   }

   //v thermal:
   VariableVThermal::VariableVThermal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].thermalRadius == 0.0) ? true : false;
   }
   VariableVThermal::~VariableVThermal() { }
   
   std::string VariableVThermal::getName() const {return popName + "/vg_v_thermal";}
   
   bool VariableVThermal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }

   bool VariableVThermal::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateNonthermal = false;
      //Calculate v nonthermal
      VNonthermalCalculation( cell, calculateNonthermal, popID, VThermal );
      const uint vectorSize = 3;
      const char* ptr = reinterpret_cast<const char*>(&VThermal);
      for (uint i = 0; i < vectorSize*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableVThermal::setSpatialCell(const SpatialCell* cell) {
      // Initialize values
      for( uint i = 0; i < 3; ++i ) {
         VThermal[i] = 0.0;
      }
      return true;
   }

   // Adding pressure calculations for nonthermal population to Vlasiator.
   // p_ij = m/3 * integral((v - <V>)_i(v - <V>)_j * f(r,v) dV)
   
   // Pressure tensor 6 components (11, 22, 33, 23, 13, 12) added by YK
   // Split into VariablePTensorNonthermalDiagonal (11, 22, 33)
   // and VariablePTensorNonthermalOffDiagonal (23, 13, 12)
   VariablePTensorNonthermalDiagonal::VariablePTensorNonthermalDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].thermalRadius == 0.0) ? true : false;
   }
   VariablePTensorNonthermalDiagonal::~VariablePTensorNonthermalDiagonal() { }
   
   std::string VariablePTensorNonthermalDiagonal::getName() const {return popName + "/vg_ptensor_nonthermal_diagonal";}
   
   bool VariablePTensorNonthermalDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }
   
   bool VariablePTensorNonthermalDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateNonthermal = true;
      //Calculate PTensor and save it in PTensorArray:
      PTensorDiagonalNonthermalCalculations( cell, calculateNonthermal, averageVX, averageVY, averageVZ, popID, PTensor );
      const uint vectorSize = 3;
      //Save the data into buffer:
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < vectorSize*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorNonthermalDiagonal::setSpatialCell(const SpatialCell* cell) {
      //Get v of the nonthermal:
      Real V[3] = {0};
      const bool calculateNonthermal = true; //We are calculating nonthermal
      VNonthermalCalculation( cell, calculateNonthermal, popID, V );
      //Set the average velocities:
      averageVX = V[0];
      averageVY = V[1];
      averageVZ = V[2];
      const uint vectorSize = 3;
      for(uint i = 0; i < vectorSize; i++) PTensor[i] = 0.0;
      return true;
   }

   // Adding pressure calculations for thermal population to Vlasiator.
   // p_ij = m/3 * integral((v - <V>)_i(v - <V>)_j * f(r,v) dV)
   
   // Pressure tensor 6 components (11, 22, 33, 23, 13, 12) added by YK
   // Split into VariablePTensorThermalDiagonal (11, 22, 33)
   // and VariablePTensorThermalOffDiagonal (23, 13, 12)
   VariablePTensorThermalDiagonal::VariablePTensorThermalDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].thermalRadius == 0.0) ? true : false;
   }
   VariablePTensorThermalDiagonal::~VariablePTensorThermalDiagonal() { }
   
   std::string VariablePTensorThermalDiagonal::getName() const {return popName + "/vg_ptensor_thermal_diagonal";}
   
   bool VariablePTensorThermalDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }
   
   bool VariablePTensorThermalDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateNonthermal = false;
      //Calculate PTensor and save it in PTensorArray:
      PTensorDiagonalNonthermalCalculations( cell, calculateNonthermal, averageVX, averageVY, averageVZ, popID, PTensor );
      const uint vectorSize = 3;
      //Save the data into buffer:
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < vectorSize*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorThermalDiagonal::setSpatialCell(const SpatialCell* cell) {
      //Get v of the thermal:
      Real V[3] = {0};
      const bool calculateNonthermal = false; //We are not calculating nonthermal
      VNonthermalCalculation( cell, calculateNonthermal, popID, V );
      //Set the average velocities:
      averageVX = V[0];
      averageVY = V[1];
      averageVZ = V[2];
      const uint vectorSize = 3;
      for(uint i = 0; i < vectorSize; i++) PTensor[i] = 0.0;
      return true;
   }

   VariablePTensorNonthermalOffDiagonal::VariablePTensorNonthermalOffDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].thermalRadius == 0.0) ? true : false;
   }
   VariablePTensorNonthermalOffDiagonal::~VariablePTensorNonthermalOffDiagonal() { }
   
   std::string VariablePTensorNonthermalOffDiagonal::getName() const {return popName + "/vg_ptensor_nonthermal_offdiagonal";}
   
   bool VariablePTensorNonthermalOffDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }
   
   bool VariablePTensorNonthermalOffDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      //Calculate PTensor for PTensorArray:
      const bool calculateNonthermal = true;
      //Calculate and save:
      PTensorOffDiagonalNonthermalCalculations( cell, calculateNonthermal, averageVX, averageVY, averageVZ, popID, PTensor );
      const uint vectorSize = 3;
      //Input data into buffer
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < 3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorNonthermalOffDiagonal::setSpatialCell(const SpatialCell* cell) {
      //Get v of the nonthermal:
      Real V[3] = {0};
      const bool calculateNonthermal = true; //We are calculating nonthermal
      VNonthermalCalculation( cell, calculateNonthermal, popID, V );
      //Set the average velocities:
      averageVX = V[0];
      averageVY = V[1];
      averageVZ = V[2];
      for(int i = 0; i < 3; i++) PTensor[i] = 0.0;
      return true;
   }

   VariablePTensorThermalOffDiagonal::VariablePTensorThermalOffDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].thermalRadius == 0.0) ? true : false;
   }
   VariablePTensorThermalOffDiagonal::~VariablePTensorThermalOffDiagonal() { }
   
   std::string VariablePTensorThermalOffDiagonal::getName() const {return popName + "/vg_ptensor_thermal_offdiagonal";}
   
   bool VariablePTensorThermalOffDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }
   
   bool VariablePTensorThermalOffDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      //Calculate PTensor for PTensorArray:
      const bool calculateNonthermal = false;
      //Calculate and save:
      PTensorOffDiagonalNonthermalCalculations( cell, calculateNonthermal, averageVX, averageVY, averageVZ, popID, PTensor );
      const uint vectorSize = 3;
      //Input data into buffer
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < 3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorThermalOffDiagonal::setSpatialCell(const SpatialCell* cell) {
      //Get v of the nonthermal:
      Real V[3] = {0};
      const bool calculateNonthermal = false; //We are not calculating nonthermal
      VNonthermalCalculation( cell, calculateNonthermal, popID, V );
      //Set the average velocities:
      averageVX = V[0];
      averageVY = V[1];
      averageVZ = V[2];
      for(int i = 0; i < 3; i++) PTensor[i] = 0.0;
      return true;
   }


   VariableEffectiveSparsityThreshold::VariableEffectiveSparsityThreshold(cuint _popID): DataReductionOperator(),popID(_popID) { 
     popName=getObjectWrapper().particleSpecies[popID].name;
   }
   VariableEffectiveSparsityThreshold::~VariableEffectiveSparsityThreshold() { }

   bool VariableEffectiveSparsityThreshold::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }

   std::string VariableEffectiveSparsityThreshold::getName() const {return popName + "/vg_effectivesparsitythreshold";}
   
   bool VariableEffectiveSparsityThreshold::reduceData(const spatial_cell::SpatialCell* cell,char* buffer) {
      Real dummy;
      reduceDiagnostic(cell,&dummy);
      const char* ptr = reinterpret_cast<const char*>(&dummy);
      for (uint i = 0; i < sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }

   bool VariableEffectiveSparsityThreshold::reduceDiagnostic(const spatial_cell::SpatialCell* cell,Real* result) {
      *result = cell->getVelocityBlockMinValue(popID);
      return true;
   }

   bool VariableEffectiveSparsityThreshold::setSpatialCell(const spatial_cell::SpatialCell* cell) {
      return true;
   }

   /*! \brief Precipitation directional differential number flux
    * Evaluation of the precipitating differential flux (per population).
    * In a selected number (default: 16) of logarithmically spaced energy bins, the average of
    *      V*V/mass
    * is calculated within the loss cone of fixed angular opening (default: 10 deg).
    * The differential flux is converted in part. / cm^2 / s / sr / eV (unit used by observers).
    * Parameters that can be set in cfg file under [{species}_precipitation]: nChannels, emin [eV], emax [eV], lossConeAngle [deg]
    * The energy channels are saved in bulk files as PrecipitationCentreEnergy{channel_number}.
    */
   VariablePrecipitationDiffFlux::VariablePrecipitationDiffFlux(cuint _popID): DataReductionOperatorHasParameters(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      lossConeAngle = getObjectWrapper().particleSpecies[popID].precipitationLossConeAngle; // deg
      emin = getObjectWrapper().particleSpecies[popID].precipitationEmin;    // already converted to SI
      emax = getObjectWrapper().particleSpecies[popID].precipitationEmax;    // already converted to SI
      nChannels = getObjectWrapper().particleSpecies[popID].precipitationNChannels; // number of energy channels, logarithmically spaced between emin and emax
      for (int i=0; i<nChannels; i++){
         channels.push_back(emin * pow(emax/emin,float(i)/(nChannels-1)));
      }
   }
   VariablePrecipitationDiffFlux::~VariablePrecipitationDiffFlux() { }
   
   std::string VariablePrecipitationDiffFlux::getName() const {return popName + "/vg_precipitationdifferentialflux";}
   
   bool VariablePrecipitationDiffFlux::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = nChannels; //Number of energy channels
      return true;
   }
   
   bool VariablePrecipitationDiffFlux::reduceData(const SpatialCell* cell,char* buffer) {

      dataDiffFlux.assign(nChannels,0.0);

      std::vector<Real> sumWeights(nChannels,0.0);

      std::array<Real,3> B;
      B[0] = cell->parameters[CellParams::PERBXVOL] +  cell->parameters[CellParams::BGBXVOL];
      B[1] = cell->parameters[CellParams::PERBYVOL] +  cell->parameters[CellParams::BGBYVOL];
      B[2] = cell->parameters[CellParams::PERBZVOL] +  cell->parameters[CellParams::BGBZVOL];

      Real cosAngle = cos(lossConeAngle*M_PI/180.0); // cosine of fixed loss cone angle

      // Unit B-field direction
      creal normB = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
      std::array<Real,3> b_unit;
      for (uint i=0; i<3; i++){
         B[i] /= normB;
      }

      // If southern hemisphere, loss cone is around -B
      if (cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ] < 0.0){
         for (uint i=0; i<3; i++){
            B[i] = -B[i];
         }
      }

      # pragma omp parallel
      {
         std::vector<Real> thread_lossCone_sum(nChannels,0.0);
         std::vector<Real> thread_count(nChannels,0.0);
         
         const Real* parameters  = cell->get_block_parameters(popID);
         const Realf* block_data = cell->get_data(popID);
         
         # pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); n++) {
            for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
               const Real VX 
                  =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                  + (i + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
               const Real VY 
                  =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
                  + (j + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
               const Real VZ 
                  =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] 
                  + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

               const Real DV3 
                  = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]
                  * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY] 
                  * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

               const Real normV = sqrt(VX*VX + VY*VY + VZ*VZ);
               const Real VdotB_norm = (B[0]*VX + B[1]*VY + B[2]*VZ)/normV;
               Real countAndGate = floor(VdotB_norm/cosAngle);  // gate function: 0 outside loss cone, 1 inside
               countAndGate = max(0.,countAndGate);
               const Real energy = 0.5 * getObjectWrapper().particleSpecies[popID].mass * normV*normV; // in SI

               // Find the correct energy bin number to update
               int binNumber = round((log(energy) - log(emin)) / log(emax/emin) * (nChannels-1));
               binNumber = max(binNumber,0); // anything < emin goes to the lowest channel
               binNumber = min(binNumber,nChannels-1); // anything > emax goes to the highest channel

               thread_lossCone_sum[binNumber] += block_data[n * SIZE_VELBLOCK + cellIndex(i,j,k)] * countAndGate * normV*normV * DV3;
               thread_count[binNumber] += countAndGate * DV3;
            }
         }

         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            for (int i=0; i<nChannels; i++) {
               dataDiffFlux[i] += thread_lossCone_sum[i];
               sumWeights[i] += thread_count[i];
            }
         }
      }

      // Averaging within each bin and conversion to unit of part. cm-2 s-1 sr-1 ev-1
      for (int i=0; i<nChannels; i++) {
         if (sumWeights[i] != 0) {
            dataDiffFlux[i] *= 1.0 / (getObjectWrapper().particleSpecies[popID].mass * sumWeights[i]) * physicalconstants::CHARGE * 1.0e-4;
         }
      }

      const char* ptr = reinterpret_cast<const char*>(dataDiffFlux.data());
      for (uint i = 0; i < nChannels*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePrecipitationDiffFlux::setSpatialCell(const SpatialCell* cell) {
      return true;
   }

   bool VariablePrecipitationDiffFlux::writeParameters(vlsv::Writer& vlsvWriter) {
      for (int i=0; i<nChannels; i++) {
         const Real channelev = channels[i]/physicalconstants::CHARGE; // in eV
         if( vlsvWriter.writeParameter(popName+"_PrecipitationCentreEnergy"+std::to_string(i), &channelev) == false ) { return false; }
      }
      if( vlsvWriter.writeParameter(popName+"_LossConeAngle", &lossConeAngle) == false ) { return false; }
      return true;
   }

   /*! \brief Energy density
    * Calculates the energy density of particles in three bins: total energy density, above E1limit*solar wind energy, and above E2limit*solar wind energy
    * Energy densities are given in eV/cm^3.
    * Parameters that can be set in cfg file under [{species}_energydensity]: 
    *    - solarwindspeed [m/s], 
    *    - solarwindenergy [eV], 
    *    - limit1 [scalar, default: 5.], 
    *    - limit2 [scalar, default: 10.].
    * The energy thresholds are saved in bulk files as parameters: 
    *    - EnergyDensityESW (in eV),
    *    - EnergyDensityELimit1 (as scalar multiplier of EnergyDensityESW),
    *    - EnergyDensityELimit2 (as scalar multiplier of EnergyDensityESW).
    */
   VariableEnergyDensity::VariableEnergyDensity(cuint _popID): DataReductionOperatorHasParameters(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      // Store internally in SI units
      solarwindenergy = getObjectWrapper().particleSpecies[popID].SolarWindEnergy;
      E1limit = solarwindenergy * getObjectWrapper().particleSpecies[popID].EnergyDensityLimit1;
      E2limit = solarwindenergy * getObjectWrapper().particleSpecies[popID].EnergyDensityLimit2;
   }
   VariableEnergyDensity::~VariableEnergyDensity() { }
   
   std::string VariableEnergyDensity::getName() const {return popName + "/vg_energydensity";}
   
   bool VariableEnergyDensity::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3; // This is not components, but rather total energy density, density over E1, and density over E2
      return true;
   }
   
   bool VariableEnergyDensity::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;

      for(int i = 0; i < 3; i++) {
         EDensity[i] = 0.0;
      }

      # pragma omp parallel
      {
         Real thread_E0_sum = 0.0;
         Real thread_E1_sum = 0.0;
         Real thread_E2_sum = 0.0;
         
         const Real* parameters  = cell->get_block_parameters(popID);
         const Realf* block_data = cell->get_data(popID);
        
         # pragma omp for
         for (vmesh::LocalID n=0; n<cell->get_number_of_velocity_blocks(popID); n++) {
            const Real DV3 
               = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]
               * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY] 
               * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

            for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
               const Real VX 
                  =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                  + (i + HALF)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
               const Real VY 
                  =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
                  + (j + HALF)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
               const Real VZ 
                  =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] 
                  + (k + HALF)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
                     
               const Real ENERGY = (VX*VX + VY*VY + VZ*VZ) * HALF * getObjectWrapper().particleSpecies[popID].mass;
               thread_E0_sum += block_data[n * SIZE_VELBLOCK+cellIndex(i,j,k)] * ENERGY * DV3;
               if (ENERGY > E1limit) thread_E1_sum += block_data[n * SIZE_VELBLOCK+cellIndex(i,j,k)] * ENERGY * DV3;
               if (ENERGY > E2limit) thread_E2_sum += block_data[n * SIZE_VELBLOCK+cellIndex(i,j,k)] * ENERGY * DV3;
            }
         }

         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            EDensity[0] += thread_E0_sum;
            EDensity[1] += thread_E1_sum;
            EDensity[2] += thread_E2_sum;
         }

      }
      // Output energy density in units eV/cm^3 instead of Joules per m^3
      EDensity[0] *= (1.0e-6)/physicalconstants::CHARGE;
      EDensity[1] *= (1.0e-6)/physicalconstants::CHARGE;
      EDensity[2] *= (1.0e-6)/physicalconstants::CHARGE;

      const char* ptr = reinterpret_cast<const char*>(&EDensity);
      for (uint i = 0; i < 3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableEnergyDensity::setSpatialCell(const SpatialCell* cell) {
      return true;
   }
   
   bool VariableEnergyDensity::writeParameters(vlsv::Writer& vlsvWriter) {
      // Output solar wind energy in eV
      Real swe = solarwindenergy/physicalconstants::CHARGE;
      // Output other bin limits as multipliers
      Real e1l = getObjectWrapper().particleSpecies[popID].EnergyDensityLimit1;
      Real e2l = getObjectWrapper().particleSpecies[popID].EnergyDensityLimit2;

      if( vlsvWriter.writeParameter(popName+"_EnergyDensityESW", &swe) == false ) { return false; }
      if( vlsvWriter.writeParameter(popName+"_EnergyDensityELimit1", &e1l) == false ) { return false; }
      if( vlsvWriter.writeParameter(popName+"_EnergyDensityELimit2", &e2l) == false ) { return false; }
      return true;
   }

} // namespace DRO
