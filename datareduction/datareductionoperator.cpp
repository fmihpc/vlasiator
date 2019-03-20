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




   
   DataReductionOperatorDerivatives::DataReductionOperatorDerivatives(const std::string& name,const unsigned int parameterIndex,const unsigned int vectorSize):
   DataReductionOperatorCellParams(name,parameterIndex,vectorSize) {

   }
   //a version with derivatives, this is the only function that is different
   bool DataReductionOperatorDerivatives::setSpatialCell(const SpatialCell* cell) {
      data  = &(cell->derivatives[_parameterIndex]);
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
   
   std::string VariableBVol::getName() const {return "B_vol";}
   
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




   //------------------ total B --------------------------------------- 
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
      for (uint i = 0; i < 3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableB::setSpatialCell(const SpatialCell* cell) {
      B[0] = cell->parameters[CellParams::PERBX] +  cell->parameters[CellParams::BGBX];
      B[1] = cell->parameters[CellParams::PERBY] +  cell->parameters[CellParams::BGBY];
      B[2] = cell->parameters[CellParams::PERBZ] +  cell->parameters[CellParams::BGBZ];
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
   
   std::string MPIrank::getName() const {return "MPI_rank";}
   
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
   
   //FsGrid cartcomm mpi rank
   FsGridRank::FsGridRank(): DataReductionOperator() { }
   FsGridRank::~FsGridRank() { }
   
   bool FsGridRank::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = sizeof(int);
      vectorSize = 1;
      return true;
   }
   
   std::string FsGridRank::getName() const {return "FSgrid_rank";}
   
   bool FsGridRank::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&fsgridRank);
      for (uint i=0; i<sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool FsGridRank::setSpatialCell(const SpatialCell* cell) {
      fsgridRank = cell->get_cell_parameters()[CellParams::FSGRID_RANK];
      return true;
   }

   //FsGrids idea of what the boundaryType ist
   FsGridBoundaryType::FsGridBoundaryType(): DataReductionOperator() { }
   FsGridBoundaryType::~FsGridBoundaryType() { }
   
   bool FsGridBoundaryType::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = sizeof(int);
      vectorSize = 1;
      return true;
   }
   
   std::string FsGridBoundaryType::getName() const {return "FSgrid_boundaryType";}
   
   bool FsGridBoundaryType::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&fsgridBoundaryType);
      for (uint i=0; i<sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool FsGridBoundaryType::setSpatialCell(const SpatialCell* cell) {
      fsgridBoundaryType = cell->get_cell_parameters()[CellParams::FSGRID_BOUNDARYTYPE];
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
   
   std::string BoundaryType::getName() const {return "Boundary_type";}
   
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
   
   std::string BoundaryLayer::getName() const {return "Boundary_layer";}
   
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
   
   std::string Blocks::getName() const {return popName + "/Blocks";}
   
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
   
   std::string VariablePressureSolver::getName() const {return "Pressure";}
   
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
   
   std::string VariablePTensorDiagonal::getName() const {return popName + "/PTensorDiagonal";}
   
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
   
   std::string VariablePTensorOffDiagonal::getName() const {return popName + "/PTensorOffDiagonal";}
   
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
   
   // Integrated divergence of magnetic field
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
   
   bool DiagnosticFluxB::reduceDiagnostic(const SpatialCell* cell,Real * result) {
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
         value += cell->parameters[CellParams::PERBX];
      } else if (cx < Parameters::xmin + 2.0 * dx && cx > Parameters::xmin + dx) {
         value += -1.0*cell->parameters[CellParams::PERBX];
      }
      if(cy > Parameters::ymax - 2.0 * dy && cy < Parameters::ymax - dy) {
         value += cell->parameters[CellParams::PERBY];
      } else if (cy < Parameters::ymin + 2.0 * dy && cy > Parameters::ymin + dy) {
         value += -1.0*cell->parameters[CellParams::PERBY];
      }
      if(cz > Parameters::zmax - 2.0 * dz && cz < Parameters::zmax - dz) {
         value += cell->parameters[CellParams::PERBZ];
      } else if (cz < Parameters::zmin + 2.0 * dz && cz > Parameters::zmin + dz) {
         value += -1.0*cell->parameters[CellParams::PERBZ];
      }
      *result = value;
      
      return true;
   }
   
   bool DiagnosticFluxB::setSpatialCell(const SpatialCell* cell) {return true;}
   
   
   
   // YK Integrated divergence of electric field
   // Integral of div E over the simulation volume =
   // Integral of flux of E on simulation volume surface
   DiagnosticFluxE::DiagnosticFluxE(): DataReductionOperator() { }
   DiagnosticFluxE::~DiagnosticFluxE() { }
   
   std::string DiagnosticFluxE::getName() const {return "FluxE";}
   
   bool DiagnosticFluxE::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   bool DiagnosticFluxE::reduceDiagnostic(const SpatialCell* cell,Real * result) {
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
         value += cell->parameters[CellParams::EX];
      } else if (cx < Parameters::xmin + 2.0 * dx && cx > Parameters::xmin + dx) {
         value += -1.0*cell->parameters[CellParams::EX];
      }
      if(cy > Parameters::ymax - 2.0 * dy && cy < Parameters::ymax - dy) {
         value += cell->parameters[CellParams::EY];
      } else if (cy < Parameters::ymin + 2.0 * dy && cy > Parameters::ymin + dy) {
         value += -1.0*cell->parameters[CellParams::EY];
      }
      if(cz > Parameters::zmax - 2.0 * dz && cz < Parameters::zmax - dz) {
         value += cell->parameters[CellParams::EZ];
      } else if (cz < Parameters::zmin + 2.0 * dz && cz > Parameters::zmin + dz) {
         value += -1.0*cell->parameters[CellParams::EZ];
      }
      *result = value;
      
      return true;
   }
   
   bool DiagnosticFluxE::setSpatialCell(const SpatialCell* cell) {return true;}
   
   // YK maximum value of the distribution function
   MaxDistributionFunction::MaxDistributionFunction(cuint _popID): DataReductionOperator(),popID(_popID) {
     popName=getObjectWrapper().particleSpecies[popID].name;
   }
   MaxDistributionFunction::~MaxDistributionFunction() { }
   
   std::string MaxDistributionFunction::getName() const {return popName + "/MaximumDistributionFunctionValue";}
   
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
   
   
   // YK minimum value of the distribution function
   MinDistributionFunction::MinDistributionFunction(cuint _popID): DataReductionOperator(),popID(_popID) {
     popName=getObjectWrapper().particleSpecies[popID].name;
   }
   MinDistributionFunction::~MinDistributionFunction() { }
   
   std::string MinDistributionFunction::getName() const {return popName + "/MinimumDistributionFunctionValue";}
   
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
	  either belonging to the backstreaming or the non-backstreaming population. 
	  There is some code duplication here, but as these helper functions are called within threads for
	  block separately, it's preferable to have them fast even at the cost of code repetition.
  ********/

   //Helper function for getting the velocity cell ids that are a part of the backstream population:
   static void getBackstreamVelocityCells(
      const Real* block_parameters,
      vector<uint64_t> & vCellIds,
      cuint popID
   ) {
      creal HALF = 0.5;
      const std::array<Real, 3> backstreamV = getObjectWrapper().particleSpecies[popID].backstreamV;
      creal backstreamRadius = getObjectWrapper().particleSpecies[popID].backstreamRadius;
      // Go through every velocity cell (i, j, k are indices)
      for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
         // Get the vx, vy, vz coordinates of the velocity cell
         const Real VX = block_parameters[BlockParams::VXCRD] + (i + HALF) * block_parameters[BlockParams::DVX];
         const Real VY = block_parameters[BlockParams::VYCRD] + (j + HALF) * block_parameters[BlockParams::DVY];
         const Real VZ = block_parameters[BlockParams::VZCRD] + (k + HALF) * block_parameters[BlockParams::DVZ];
         // Compare the distance of the velocity cell from the center of the maxwellian distribution to the radius of the maxwellian distribution
         if( ( (backstreamV[0] - VX) * (backstreamV[0] - VX)
             + (backstreamV[1] - VY) * (backstreamV[1] - VY)
             + (backstreamV[2] - VZ) * (backstreamV[2] - VZ) )
             >
             backstreamRadius*backstreamRadius ) {
             //The velocity cell is a part of the backstream population:
             vCellIds.push_back(cellIndex(i,j,k));
          }
      }
   }
   //Helper function for getting the velocity cell ids that are a part of the backstream population:
   static void getNonBackstreamVelocityCells(
      const Real* block_parameters,
      vector<uint64_t> & vCellIds,
      cuint popID
   ) {
      creal HALF = 0.5;
      const std::array<Real, 3> backstreamV = getObjectWrapper().particleSpecies[popID].backstreamV;
      creal backstreamRadius = getObjectWrapper().particleSpecies[popID].backstreamRadius;
      for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
         const Real VX = block_parameters[BlockParams::VXCRD] + (i + HALF) * block_parameters[BlockParams::DVX];
         const Real VY = block_parameters[BlockParams::VYCRD] + (j + HALF) * block_parameters[BlockParams::DVY];
         const Real VZ = block_parameters[BlockParams::VZCRD] + (k + HALF) * block_parameters[BlockParams::DVZ];
         if( ( (backstreamV[0] - VX) * (backstreamV[0] - VX)
             + (backstreamV[1] - VY) * (backstreamV[1] - VY)
             + (backstreamV[2] - VZ) * (backstreamV[2] - VZ) )
             <=
             backstreamRadius*backstreamRadius ) {
             //The velocity cell is not a part of the backstream population:
             vCellIds.push_back(cellIndex(i,j,k));
          }
      }
   }
   //Helper function for getting the velocity cell indices that are a part of the backstream population:
   static void getBackstreamVelocityCellIndices(
      const Real* block_parameters,
      vector<array<uint, 3>> & vCellIndices,
      cuint popID
   ) {
      creal HALF = 0.5;
      const std::array<Real, 3> backstreamV = getObjectWrapper().particleSpecies[popID].backstreamV;
      creal backstreamRadius = getObjectWrapper().particleSpecies[popID].backstreamRadius;
      // Go through a block's every velocity cell
      for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
         // Get the coordinates of the velocity cell (e.g. VX = block_vx_min_coordinates + (velocity_cell_indice_x+0.5)*length_of_velocity_cell_in_x_direction
         const Real VX = block_parameters[BlockParams::VXCRD] + (i + HALF) * block_parameters[BlockParams::DVX];
         const Real VY = block_parameters[BlockParams::VYCRD] + (j + HALF) * block_parameters[BlockParams::DVY];
         const Real VZ = block_parameters[BlockParams::VZCRD] + (k + HALF) * block_parameters[BlockParams::DVZ];
         // Calculate the distance of the velocity cell from the center of the maxwellian distribution and compare it to the approximate radius of the maxwellian distribution
         if( ( (backstreamV[0] - VX) * (backstreamV[0] - VX)
             + (backstreamV[1] - VY) * (backstreamV[1] - VY)
             + (backstreamV[2] - VZ) * (backstreamV[2] - VZ) )
             >
             backstreamRadius*backstreamRadius ) {
             //The velocity cell is a part of the backstream population because it is not within the radius:
             const array<uint, 3> indices{{i, j, k}};
             vCellIndices.push_back( indices );
          }
      }
   }
   //Helper function for getting the velocity cell indices that are not a part of the backstream population:
   static void getNonBackstreamVelocityCellIndices(
      const Real* block_parameters,
      vector<array<uint, 3>> & vCellIndices,
      cuint popID
   ) {
      creal HALF = 0.5;
      const std::array<Real, 3> backstreamV = getObjectWrapper().particleSpecies[popID].backstreamV;
      creal backstreamRadius = getObjectWrapper().particleSpecies[popID].backstreamRadius;
      // Go through a block's every velocity cell
      for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
         // Get the coordinates of the velocity cell (e.g. VX = block_vx_min_coordinates + (velocity_cell_indice_x+0.5)*length_of_velocity_cell_in_x_direction
         const Real VX = block_parameters[BlockParams::VXCRD] + (i + HALF) * block_parameters[BlockParams::DVX];
         const Real VY = block_parameters[BlockParams::VYCRD] + (j + HALF) * block_parameters[BlockParams::DVY];
         const Real VZ = block_parameters[BlockParams::VZCRD] + (k + HALF) * block_parameters[BlockParams::DVZ];
         // Calculate the distance of the velocity cell from the center of the maxwellian distribution and compare it to the approximate radius of the maxwellian distribution
         if( ( (backstreamV[0] - VX) * (backstreamV[0] - VX)
             + (backstreamV[1] - VY) * (backstreamV[1] - VY)
             + (backstreamV[2] - VZ) * (backstreamV[2] - VZ) )
             <=
             backstreamRadius*backstreamRadius ) {
             //The velocity cell is not a part of the backstream population because it is within the radius:
             const array<uint, 3> indices{{i, j, k}};
             vCellIndices.push_back( indices );
          }
      }
   }

  /********
	   Next level of helper functions - these include threading and calculate zeroth or first velocity moments or the
	   diagonal / off-diagonal pressure tensor components for
	   backstreaming or non-backstreaming populations  ********/

   //Calculates rho backstream or rho non backstream
   static void rhoBackstreamCalculation( const SpatialCell * cell, const bool calculateBackstream, cuint popID, Real & rho ) {
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
            if ( calculateBackstream == true ) {
               getBackstreamVelocityCells(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCells, popID);
            } else {
               getNonBackstreamVelocityCells(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCells, popID);
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

   static void VBackstreamCalculation( const SpatialCell * cell, const bool calculateBackstream, cuint popID, Real * V ) {
      const Real HALF = 0.5;
      // Make sure the V is initialized
      V[0] = 0;
      V[1] = 0;
      V[2] = 0;
      Real n = 0;
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
            // Get the velocity cell indices of the cells that are a part of the backstream population
            vector< array<uint, 3> > vCellIndices;
            vCellIndices.clear();
            // Save indices to the std::vector
            if( calculateBackstream == true ) {
               getBackstreamVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            } else {
               getNonBackstreamVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            }
            // We have now fethced all of the needed velocity cell indices, so now go through them:
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
            n += thread_n_sum;
         }
      }

      // Finally, divide n*V by V.
      V[0]/=n;
      V[1]/=n;
      V[2]/=n;
      return;
   }

   static void PTensorDiagonalBackstreamCalculations( const SpatialCell * cell,
                                                      const bool calculateBackstream,
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
            if( calculateBackstream == true ) {
               getBackstreamVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            } else {
               getNonBackstreamVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
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

   static void PTensorOffDiagonalBackstreamCalculations( const SpatialCell * cell,
                                                         const bool calculateBackstream,
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
            if( calculateBackstream == true ) {
               getBackstreamVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
            } else {
               getNonBackstreamVelocityCellIndices(&parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS], vCellIndices, popID);
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
	     End velocity moment / backstream/non-backstreamn helper functions
  *********/


   VariableMeshData::VariableMeshData(): DataReductionOperatorHandlesWriting() { }
   VariableMeshData::~VariableMeshData() { }
   
   std::string VariableMeshData::getName() const {return "MeshData";}
   
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
   
   // Rho backstream:
   VariableRhoBackstream::VariableRhoBackstream(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].backstreamRadius == 0.0) ? true : false;
   }
   VariableRhoBackstream::~VariableRhoBackstream() { }
   
   std::string VariableRhoBackstream::getName() const {return popName + "/RhoBackstream";}
   
   bool VariableRhoBackstream::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 1;
      return true;
   }
   
   // Adding rho backstream calculations to Vlasiator.
   bool VariableRhoBackstream::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateBackstream = true;
      rhoBackstreamCalculation( cell, calculateBackstream, popID, RhoBackstream );
      const char* ptr = reinterpret_cast<const char*>(&RhoBackstream);
      for (uint i = 0; i < sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoBackstream::setSpatialCell(const SpatialCell* cell) {
      RhoBackstream = 0.0;
      return true;
   }


   // Rho non backstream:
   VariableRhoNonBackstream::VariableRhoNonBackstream(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].backstreamRadius == 0.0) ? true : false;
   }
   VariableRhoNonBackstream::~VariableRhoNonBackstream() { }
   
   std::string VariableRhoNonBackstream::getName() const {return popName + "/RhoNonBackstream";}
   
   bool VariableRhoNonBackstream::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 1;
      return true;
   }
   
   // Rho non backstream calculation.
   bool VariableRhoNonBackstream::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateBackstream = false; //We don't want backstream
      rhoBackstreamCalculation( cell, calculateBackstream, popID, Rho );
      const char* ptr = reinterpret_cast<const char*>(&Rho);
      for (uint i = 0; i < sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoNonBackstream::setSpatialCell(const SpatialCell* cell) {
      Rho = 0.0;
      return true;
   }

   // v backstream:
   VariableVBackstream::VariableVBackstream(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].backstreamRadius == 0.0) ? true : false;
   }
   VariableVBackstream::~VariableVBackstream() { }
   
   std::string VariableVBackstream::getName() const {return popName + "/VBackstream";}
   
   bool VariableVBackstream::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }

   // Adding v backstream calculations to Vlasiator.
   bool VariableVBackstream::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateBackstream = true;
      //Calculate v backstream
      VBackstreamCalculation( cell, calculateBackstream, popID, VBackstream );
      const uint VBackstreamSize = 3;
      const char* ptr = reinterpret_cast<const char*>(&VBackstream);
      for (uint i = 0; i < VBackstreamSize*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableVBackstream::setSpatialCell(const SpatialCell* cell) {
      // Initialize values
      for( uint i = 0; i < 3; ++i ) {
         VBackstream[i] = 0.0;
      }
      return true;
   }

   //v non backstream:
   VariableVNonBackstream::VariableVNonBackstream(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].backstreamRadius == 0.0) ? true : false;
   }
   VariableVNonBackstream::~VariableVNonBackstream() { }
   
   std::string VariableVNonBackstream::getName() const {return popName + "/VNonBackstream";}
   
   bool VariableVNonBackstream::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }

   // Adding v non backstream calculations to Vlasiator.
   bool VariableVNonBackstream::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateBackstream = false;
      //Calculate v backstream
      VBackstreamCalculation( cell, calculateBackstream, popID, V );
      const uint vectorSize = 3;
      const char* ptr = reinterpret_cast<const char*>(&V);
      for (uint i = 0; i < vectorSize*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableVNonBackstream::setSpatialCell(const SpatialCell* cell) {
      // Initialize values
      for( uint i = 0; i < 3; ++i ) {
         V[i] = 0.0;
      }
      return true;
   }

   // Adding pressure calculations for backstream population to Vlasiator.
   // p_ij = m/3 * integral((v - <V>)_i(v - <V>)_j * f(r,v) dV)
   
   // Pressure tensor 6 components (11, 22, 33, 23, 13, 12) added by YK
   // Split into VariablePTensorBackstreamDiagonal (11, 22, 33)
   // and VariablePTensorOffDiagonal (23, 13, 12)
   VariablePTensorBackstreamDiagonal::VariablePTensorBackstreamDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].backstreamRadius == 0.0) ? true : false;
   }
   VariablePTensorBackstreamDiagonal::~VariablePTensorBackstreamDiagonal() { }
   
   std::string VariablePTensorBackstreamDiagonal::getName() const {return popName + "/PTensorBackstreamDiagonal";}
   
   bool VariablePTensorBackstreamDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }
   
   bool VariablePTensorBackstreamDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateBackstream = true;
      //Calculate PTensor and save it in PTensorArray:
      PTensorDiagonalBackstreamCalculations( cell, calculateBackstream, averageVX, averageVY, averageVZ, popID, PTensor );
      const uint vectorSize = 3;
      //Save the data into buffer:
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < vectorSize*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorBackstreamDiagonal::setSpatialCell(const SpatialCell* cell) {
      //Get v of the backstream:
      Real V[3] = {0};
      const bool calculateBackstream = true; //We are calculating backstream
      VBackstreamCalculation( cell, calculateBackstream, popID, V );
      //Set the average velocities:
      averageVX = V[0];
      averageVY = V[1];
      averageVZ = V[2];
      const uint vectorSize = 3;
      for(uint i = 0; i < vectorSize; i++) PTensor[i] = 0.0;
      return true;
   }

   // Adding pressure calculations for backstream population to Vlasiator.
   // p_ij = m/3 * integral((v - <V>)_i(v - <V>)_j * f(r,v) dV)
   
   // Pressure tensor 6 components (11, 22, 33, 23, 13, 12) added by YK
   // Split into VariablePTensorNonBackstreamDiagonal (11, 22, 33)
   // and VariablePTensorOffDiagonal (23, 13, 12)
   VariablePTensorNonBackstreamDiagonal::VariablePTensorNonBackstreamDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].backstreamRadius == 0.0) ? true : false;
   }
   VariablePTensorNonBackstreamDiagonal::~VariablePTensorNonBackstreamDiagonal() { }
   
   std::string VariablePTensorNonBackstreamDiagonal::getName() const {return popName + "/PTensorNonBackstreamDiagonal";}
   
   bool VariablePTensorNonBackstreamDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }
   
   bool VariablePTensorNonBackstreamDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      const bool calculateBackstream = false;
      //Calculate PTensor and save it in PTensorArray:
      PTensorDiagonalBackstreamCalculations( cell, calculateBackstream, averageVX, averageVY, averageVZ, popID, PTensor );
      const uint vectorSize = 3;
      //Save the data into buffer:
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < vectorSize*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorNonBackstreamDiagonal::setSpatialCell(const SpatialCell* cell) {
      //Get v of the backstream:
      Real V[3] = {0};
      const bool calculateBackstream = false; //We are not calculating backstream
      VBackstreamCalculation( cell, calculateBackstream, popID, V );
      //Set the average velocities:
      averageVX = V[0];
      averageVY = V[1];
      averageVZ = V[2];
      const uint vectorSize = 3;
      for(uint i = 0; i < vectorSize; i++) PTensor[i] = 0.0;
      return true;
   }

   VariablePTensorBackstreamOffDiagonal::VariablePTensorBackstreamOffDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].backstreamRadius == 0.0) ? true : false;
   }
   VariablePTensorBackstreamOffDiagonal::~VariablePTensorBackstreamOffDiagonal() { }
   
   std::string VariablePTensorBackstreamOffDiagonal::getName() const {return popName + "/PTensorBackstreamOffDiagonal";}
   
   bool VariablePTensorBackstreamOffDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }
   
   bool VariablePTensorBackstreamOffDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      //Calculate PTensor for PTensorArray:
      const bool calculateBackstream = true;
      //Calculate and save:
      PTensorOffDiagonalBackstreamCalculations( cell, calculateBackstream, averageVX, averageVY, averageVZ, popID, PTensor );
      const uint vectorSize = 3;
      //Input data into buffer
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < 3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorBackstreamOffDiagonal::setSpatialCell(const SpatialCell* cell) {
      //Get v of the backstream:
      Real V[3] = {0};
      const bool calculateBackstream = true; //We are calculating backstream
      VBackstreamCalculation( cell, calculateBackstream, popID, V );
      //Set the average velocities:
      averageVX = V[0];
      averageVY = V[1];
      averageVZ = V[2];
      for(int i = 0; i < 3; i++) PTensor[i] = 0.0;
      return true;
   }

   VariablePTensorNonBackstreamOffDiagonal::VariablePTensorNonBackstreamOffDiagonal(cuint _popID): DataReductionOperator(),popID(_popID) {
      popName = getObjectWrapper().particleSpecies[popID].name;
      doSkip = (getObjectWrapper().particleSpecies[popID].backstreamRadius == 0.0) ? true : false;
   }
   VariablePTensorNonBackstreamOffDiagonal::~VariablePTensorNonBackstreamOffDiagonal() { }
   
   std::string VariablePTensorNonBackstreamOffDiagonal::getName() const {return popName + "/PTensorNonBackstreamOffDiagonal";}
   
   bool VariablePTensorNonBackstreamOffDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = (doSkip == true) ? 0 : 3;
      return true;
   }
   
   bool VariablePTensorNonBackstreamOffDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      //Calculate PTensor for PTensorArray:
      const bool calculateBackstream = false;
      //Calculate and save:
      PTensorOffDiagonalBackstreamCalculations( cell, calculateBackstream, averageVX, averageVY, averageVZ, popID, PTensor );
      const uint vectorSize = 3;
      //Input data into buffer
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i = 0; i < 3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorNonBackstreamOffDiagonal::setSpatialCell(const SpatialCell* cell) {
      //Get v of the backstream:
      Real V[3] = {0};
      const bool calculateBackstream = false; //We are not calculating backstream
      VBackstreamCalculation( cell, calculateBackstream, popID, V );
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

   std::string VariableEffectiveSparsityThreshold::getName() const {return popName + "/EffectiveSparsityThreshold";}
   
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
} // namespace DRO
