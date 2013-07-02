/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */

#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <limits>

#include "datareductionoperator.h"
#include "../vlscommon.h"
//#include "../parameters.h"

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
   
   // TODO update this documentation snippet.
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
   
   // TODO update this documentation snippet.
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
   
   
   
   
   
   DataReductionOperatorCellParams::DataReductionOperatorCellParams(const std::string& name,const unsigned int parameterIndex,const unsigned int vectorSize):
   DataReductionOperator() {
      _vectorSize=vectorSize;
      _name=name;
      _parameterIndex=parameterIndex;
   }
   DataReductionOperatorCellParams::~DataReductionOperatorCellParams() { }
   
   bool DataReductionOperatorCellParams::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = _vectorSize;
      return true;
   }
   
   std::string DataReductionOperatorCellParams::getName() const {return _name;}
   
   bool DataReductionOperatorCellParams::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(_data);
      for (uint i=0; i<_vectorSize*sizeof(Real); ++i){
         buffer[i] = ptr[i];
      }
      return true;
   }
   
   bool DataReductionOperatorCellParams::reduceData(const SpatialCell* cell,Real* buffer){
      //If _vectorSize is >1 it still works, we just give the first value and no other ones..
      *buffer=_data[0];
      return true;
   }
   bool DataReductionOperatorCellParams::setSpatialCell(const SpatialCell* cell) {
      _data  = &(cell->parameters[_parameterIndex]);
      return true;
   }




   
   DataReductionOperatorDerivatives::DataReductionOperatorDerivatives(const std::string& name,const unsigned int parameterIndex,const unsigned int vectorSize):
   DataReductionOperatorCellParams(name,parameterIndex,vectorSize) {

   }
   //a version with derivatives, this is the only function that is different
   bool DataReductionOperatorDerivatives::setSpatialCell(const SpatialCell* cell) {
      _data  = &(cell->derivatives[_parameterIndex]);
      return true;
   }


   DataReductionOperatorBVOLDerivatives::DataReductionOperatorBVOLDerivatives(const std::string& name,const unsigned int parameterIndex,const unsigned int vectorSize):
   DataReductionOperatorCellParams(name,parameterIndex,vectorSize) {
      
   }
   //a version with derivatives, this is the only function that is different
   bool DataReductionOperatorBVOLDerivatives::setSpatialCell(const SpatialCell* cell) {
      _data  = &(cell->derivativesBVOL[_parameterIndex]);
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
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableBVol::setSpatialCell(const SpatialCell* cell) {
      B[0] = cell->parameters[CellParams::PERBXVOL] +  cell->parameters[CellParams::BGBXVOL];
      B[1] = cell->parameters[CellParams::PERBYVOL] +  cell->parameters[CellParams::BGBYVOL];
      B[2] = cell->parameters[CellParams::PERBZVOL] +  cell->parameters[CellParams::BGBZVOL];
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
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableB::setSpatialCell(const SpatialCell* cell) {
      B[0] = cell->parameters[CellParams::PERBX] +  cell->parameters[CellParams::BGBX];
      B[1] = cell->parameters[CellParams::PERBY] +  cell->parameters[CellParams::BGBY];
      B[2] = cell->parameters[CellParams::PERBZ] +  cell->parameters[CellParams::BGBZ];
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
      for (uint i=0; i<sizeof(int); ++i) buffer[i] = ptr[i];
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
      for (uint i=0; i<sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool BoundaryLayer::setSpatialCell(const SpatialCell* cell) {
      boundaryLayer = (int)cell->sysBoundaryLayer;
      return true;
   }

   // VelocitySubSteps
   VelocitySubSteps::VelocitySubSteps(): DataReductionOperator() { }
   VelocitySubSteps::~VelocitySubSteps() { }
   
   bool VelocitySubSteps::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "int";
      dataSize = sizeof(int);
      vectorSize = 1;
      return true;
   }
   
   std::string VelocitySubSteps::getName() const {return "Velocity_substeps";}
   
   bool VelocitySubSteps::reduceData(const SpatialCell* cell,char* buffer) {
      const char* ptr = reinterpret_cast<const char*>(&substeps);
      for (uint i=0; i<sizeof(int); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VelocitySubSteps::setSpatialCell(const SpatialCell* cell) {
      substeps = (int)cell->subStepsAcceleration;
      return true;
   }

   
   // Blocks
   Blocks::Blocks(): DataReductionOperator() { }
   Blocks::~Blocks() { }
   
   bool Blocks::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "uint";
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
   
   
   // Scalar pressure 
   VariablePressure::VariablePressure(): DataReductionOperator() { }
   VariablePressure::~VariablePressure() { }
   
   std::string VariablePressure::getName() const {return "Pressure";}
   
   bool VariablePressure::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   // Adding pressure calculations to Vlasiator.
   // p = m/3 * integral((v - <V>)^2 * f(r,v) dV), doing the sum of the x, y and z components.
   bool VariablePressure::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      const Real THIRD = 1.0/3.0;
      # pragma omp parallel
      {
         Real thread_nvx2_sum = 0.0;
         Real thread_nvy2_sum = 0.0;
         Real thread_nvz2_sum = 0.0;
         # pragma omp for
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
                     
                     thread_nvx2_sum += block-> data[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX) * DV3;
                     thread_nvy2_sum += block-> data[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY) * DV3;
                     thread_nvz2_sum += block-> data[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ) * DV3;
            }
         }
         
         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            Pressure += physicalconstants::MASS_PROTON * THIRD * (thread_nvx2_sum + thread_nvy2_sum + thread_nvz2_sum);
         }
      }
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
      # pragma omp parallel
      {
         Real thread_nvxvx_sum = 0.0;
         Real thread_nvyvy_sum = 0.0;
         Real thread_nvzvz_sum = 0.0;
         # pragma omp for
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
                     
                     thread_nvxvx_sum += block-> data[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX) * DV3;
                     thread_nvyvy_sum += block-> data[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY) * DV3;
                     thread_nvzvz_sum += block-> data[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ) * DV3;
            }
         }
         
         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            PTensor[0] += physicalconstants::MASS_PROTON * THIRD * thread_nvxvx_sum;
            PTensor[1] += physicalconstants::MASS_PROTON * THIRD * thread_nvyvy_sum;
            PTensor[2] += physicalconstants::MASS_PROTON * THIRD * thread_nvzvz_sum;
         }
      }
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
      # pragma omp parallel
      {
         Real thread_nvxvy_sum = 0.0;
         Real thread_nvzvx_sum = 0.0;
         Real thread_nvyvz_sum = 0.0;
         # pragma omp for
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
                     
                     thread_nvxvy_sum += block-> data[cellIndex(i,j,k)] * (VX - averageVX) * (VY - averageVY) * DV3;
                     thread_nvzvx_sum += block-> data[cellIndex(i,j,k)] * (VZ - averageVZ) * (VX - averageVX) * DV3;
                     thread_nvyvz_sum += block-> data[cellIndex(i,j,k)] * (VY - averageVY) * (VZ - averageVZ) * DV3;
            }
         }
         
         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            PTensor[0] += physicalconstants::MASS_PROTON * THIRD * thread_nvyvz_sum;
            PTensor[1] += physicalconstants::MASS_PROTON * THIRD * thread_nvzvx_sum;
            PTensor[2] += physicalconstants::MASS_PROTON * THIRD * thread_nvxvy_sum;
         }
      }
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
   
   bool DiagnosticFluxE::reduceData(const SpatialCell* cell,Real * result) {
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
   MaxDistributionFunction::MaxDistributionFunction(): DataReductionOperator() { }
   MaxDistributionFunction::~MaxDistributionFunction() { }
   
   std::string MaxDistributionFunction::getName() const {return "MaximumDistributionFunctionValue";}
   
   bool MaxDistributionFunction::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   
   bool MaxDistributionFunction::reduceData(const SpatialCell* cell,Real* buffer) {
      const Real HALF = 0.5;
      
      #pragma omp parallel 
      {
         Real threadMax = std::numeric_limits<Real>::min();
         #pragma omp for
         for(uint n=0; n<cell->number_of_blocks; n++) {
            unsigned int blockId = cell->velocity_block_list[n];
            const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
            for (uint k=0; k<WID; ++k)
               for (uint j=0; j<WID; ++j)
                  for (uint i=0; i<WID; ++i) {
                     const int celli=k*WID*WID+j*WID+i;
                     threadMax = max(block->data[celli], threadMax);
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
      reduceData(cell,&dummy);
      const char* ptr = reinterpret_cast<const char*>(&dummy);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool MaxDistributionFunction::setSpatialCell(const SpatialCell* cell) {
      return true;
   }
   
   
   // YK minimum value of the distribution function
   MinDistributionFunction::MinDistributionFunction(): DataReductionOperator() { }
   MinDistributionFunction::~MinDistributionFunction() { }
   
   std::string MinDistributionFunction::getName() const {return "MinimumDistributionFunctionValue";}
   
   bool MinDistributionFunction::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   
   bool MinDistributionFunction::reduceData(const SpatialCell* cell,Real* buffer) {
      const Real HALF = 0.5;
      
      #pragma omp parallel 
      {
         Real threadMin = std::numeric_limits<Real>::max();
         #pragma omp for
         for(uint n=0; n<cell->number_of_blocks; n++) {
            unsigned int blockId = cell->velocity_block_list[n];
            const Velocity_Block* block = cell->at(blockId); //returns a reference to block
            for (uint k=0; k<WID; ++k)
               for (uint j=0; j<WID; ++j)
                  for (uint i=0; i<WID; ++i) {
                     const int celli=k*WID*WID+j*WID+i;
                     threadMin = min(block->data[celli], threadMin);
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
      reduceData(cell,&dummy);
      const char* ptr = reinterpret_cast<const char*>(&dummy);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool MinDistributionFunction::setSpatialCell(const SpatialCell* cell) {
      return true;
   }

   //Helper function for getting the velocity cell ids that are a part of the backstream population:
   static void getBackstreamVelocityCells( const Velocity_Block * block, vector<uint64_t> & vCellIds ) {
      const Real HALF = 0.5;
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
         const Real VX = block-> parameters[BlockParams::VXCRD] + (i+HALF) * block-> parameters[BlockParams::DVX];
         const Real VY = block-> parameters[BlockParams::VYCRD] + (j+HALF) * block-> parameters[BlockParams::DVY];
         const Real VZ = block-> parameters[BlockParams::VZCRD] + (k+HALF) * block-> parameters[BlockParams::DVZ];
         if( ( (P::backstreamvx - VX)*(P::backstreamvx - VX)
             + (P::backstreamvy - VY)*(P::backstreamvy - VY)
             + (P::backstreamvz - VZ)*(P::backstreamvz - VZ) )
             >
             P::backstreamradius*P::backstreamradius ) {
             //The velocity cell is a part of the backstream population:
             vCellIds.push_back(cellIndex(i,j,k));
          }
      }
   }
   //Helper function for getting the velocity cell ids that are a part of the backstream population as well as their coordinates:
   static void getBackstreamVelocityCellIndices( const Velocity_Block * block, 
                                                 vector<array<uint, 3>> & vCellIndices ) {
      const Real HALF = 0.5;
      for (uint k=0; k<WID; ++k) for (uint j=0; j<WID; ++j) for (uint i=0; i<WID; ++i) {
         const Real VX = block-> parameters[BlockParams::VXCRD] + (i+HALF) * block-> parameters[BlockParams::DVX];
         const Real VY = block-> parameters[BlockParams::VYCRD] + (j+HALF) * block-> parameters[BlockParams::DVY];
         const Real VZ = block-> parameters[BlockParams::VZCRD] + (k+HALF) * block-> parameters[BlockParams::DVZ];
         if( ( (P::backstreamvx - VX)*(P::backstreamvx - VX)
             + (P::backstreamvy - VY)*(P::backstreamvy - VY)
             + (P::backstreamvz - VZ)*(P::backstreamvz - VZ) )
             >
             P::backstreamradius*P::backstreamradius ) {
             //The velocity cell is a part of the backstream population:
             const array<uint, 3> indices{{i, j, k}};
             vCellIndices.push_back( indices );
          }
      }
   }

   // Rho backstream:
   VariableRhoBackstream::VariableRhoBackstream(): DataReductionOperator() { }
   VariableRhoBackstream::~VariableRhoBackstream() { }
   
   std::string VariableRhoBackstream::getName() const {return "RhoBackstream";}
   
   bool VariableRhoBackstream::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   // Adding rho backstream calculations to Vlasiator.
   bool VariableRhoBackstream::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      # pragma omp parallel
      {
         Real thread_n_sum = 0.0;
         # pragma omp for
         for(uint n=0; n<cell->number_of_blocks; n++) {
            const unsigned int blockId = cell->velocity_block_list[n];
            const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
            const Real DV3 = block-> parameters[BlockParams::DVX] * block-> parameters[BlockParams::DVY] * block-> parameters[BlockParams::DVZ];
            vector< uint64_t > vCells; //Velocity cell ids
            getBackstreamVelocityCells(block, vCells);
            for( vector< uint64_t >::const_iterator it = vCells.begin(); it != vCells.end(); ++it ) {
               //velocity cell id = *it
               thread_n_sum += block-> data[(*it)];
            }
            //Multiply by the volume:
            thread_n_sum = thread_n_sum * DV3;
         }

         // Accumulate contributions coming from this velocity block
         // If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            RhoBackstream += thread_n_sum;
         }
      }
      const char* ptr = reinterpret_cast<const char*>(&RhoBackstream);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoBackstream::setSpatialCell(const SpatialCell* cell) {
      RhoBackstream = 0.0;
      return true;
   }


   // Rho non backstream:
   VariableRhoNonBackstream::VariableRhoNonBackstream(): DataReductionOperator() { }
   VariableRhoNonBackstream::~VariableRhoNonBackstream() { }
   
   std::string VariableRhoNonBackstream::getName() const {return "RhoNonBackstream";}
   
   bool VariableRhoNonBackstream::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   // Adding rho non backstream calculations to Vlasiator.
   bool VariableRhoNonBackstream::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      # pragma omp parallel
      {
         Real thread_n_sum = 0.0;
         # pragma omp for
         for(uint n=0; n<cell->number_of_blocks; n++) {
            const unsigned int blockId = cell->velocity_block_list[n];
            const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
            const Real DV3 = block-> parameters[BlockParams::DVX] * block-> parameters[BlockParams::DVY] * block-> parameters[BlockParams::DVZ];
            vector< uint64_t > vCells; //Velocity cell ids
            getBackstreamVelocityCells(block, vCells);
            for( vector< uint64_t >::const_iterator it = vCells.begin(); it != vCells.end(); ++it ) {
               //velocity cell id = *it
               thread_n_sum += block-> data[(*it)];
            }
            //Multiply by the volume:
            thread_n_sum = thread_n_sum * DV3;
         }

         // Accumulate contributions coming from this velocity block
         // If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            Rho += thread_n_sum;
         }
      }
      const char* ptr = reinterpret_cast<const char*>(&Rho);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoNonBackstream::setSpatialCell(const SpatialCell* cell) {
      Rho = 0.0;
      return true;
   }

   //Rho v backstream:
   VariableRhoVBackstream::VariableRhoVBackstream(): DataReductionOperator() { }
   VariableRhoVBackstream::~VariableRhoVBackstream() { }
   
   std::string VariableRhoVBackstream::getName() const {return "RhoVBackstream";}
   
   bool VariableRhoVBackstream::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }



   // Adding rho v backstream calculations to Vlasiator.
   bool VariableRhoVBackstream::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      # pragma omp parallel
      {
         Real thread_nvx_sum = 0.0;
         Real thread_nvy_sum = 0.0;
         Real thread_nvz_sum = 0.0;
         # pragma omp for
         for(uint n=0; n<cell->number_of_blocks; n++) {
            const unsigned int blockId = cell->velocity_block_list[n];
            const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
            const Real DV3 = block-> parameters[BlockParams::DVX] * block-> parameters[BlockParams::DVY] * block-> parameters[BlockParams::DVZ];
            vector< array<uint, 3> > vCellIndices;
            getBackstreamVelocityCellIndices(block, vCellIndices);
            for( vector< array<uint, 3> >::const_iterator it = vCellIndices.begin(); it != vCellIndices.end(); ++it ) {
               const array<uint, 3> indices = *it;
               const uint i = indices[0];
               const uint j = indices[1];
               const uint k = indices[2];
               const Real VX = block-> parameters[BlockParams::VXCRD] + (i+HALF) * block-> parameters[BlockParams::DVX];
               const Real VY = block-> parameters[BlockParams::VYCRD] + (j+HALF) * block-> parameters[BlockParams::DVY];
               const Real VZ = block-> parameters[BlockParams::VZCRD] + (k+HALF) * block-> parameters[BlockParams::DVZ];
               thread_nvx_sum += block->data[cellIndex(i,j,k)]*VX;
               thread_nvy_sum += block->data[cellIndex(i,j,k)]*VY;
               thread_nvz_sum += block->data[cellIndex(i,j,k)]*VZ;
            }
            thread_nvx_sum = thread_nvx_sum * DV3;
            thread_nvy_sum = thread_nvy_sum * DV3;
            thread_nvz_sum = thread_nvz_sum * DV3;
         }

         // Accumulate contributions coming from this velocity block
         // If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            RhoVBackstream[0] += thread_nvx_sum;
            RhoVBackstream[1] += thread_nvy_sum;
            RhoVBackstream[2] += thread_nvz_sum;
         }
      }
      const uint RhoVBackstreamSize = 3;
      const char* ptr = reinterpret_cast<const char*>(&RhoVBackstream);
      for (uint i=0; i<RhoVBackstreamSize*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariableRhoVBackstream::setSpatialCell(const SpatialCell* cell) {
      for( uint i = 0; i < 3; ++i ) {
         RhoVBackstream[i] = 0.0;
      }
      return true;
   }

   // Scalar pressure of backstream
   VariablePressureBackstream::VariablePressureBackstream(): DataReductionOperator() { }
   VariablePressureBackstream::~VariablePressureBackstream() { }
   
   std::string VariablePressureBackstream::getName() const {return "PressureBackstream";}
   
   bool VariablePressureBackstream::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   // Adding pressure backstream calculations to Vlasiator.
   // p = m/3 * integral((v - <V>)^2 * f(r,v) dV), doing the sum of the x, y and z components.
   bool VariablePressureBackstream::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      const Real THIRD = 1.0/3.0;
      # pragma omp parallel
      {
         Real thread_nvx2_sum = 0.0;
         Real thread_nvy2_sum = 0.0;
         Real thread_nvz2_sum = 0.0;
         # pragma omp for
         for(uint n=0; n<cell->number_of_blocks; n++) {
            unsigned int blockId = cell->velocity_block_list[n];
            const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
            const Real DV3 = block-> parameters[BlockParams::DVX] * block-> parameters[BlockParams::DVY] * block-> parameters[BlockParams::DVZ];
            vector< array<uint, 3> > vCellIndices;
            getBackstreamVelocityCellIndices(block, vCellIndices);
            for( vector< array<uint, 3> >::const_iterator it = vCellIndices.begin(); it != vCellIndices.end(); ++it ) {
               //Go through every velocity cell:
               const array<uint, 3> indices = *it;
               const uint i = indices[0];
               const uint j = indices[1];
               const uint k = indices[2];
               const Real VX = block-> parameters[BlockParams::VXCRD] + (i+HALF) * block-> parameters[BlockParams::DVX];
               const Real VY = block-> parameters[BlockParams::VYCRD] + (j+HALF) * block-> parameters[BlockParams::DVY];
               const Real VZ = block-> parameters[BlockParams::VZCRD] + (k+HALF) * block-> parameters[BlockParams::DVZ];
               thread_nvx2_sum += block-> data[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX);
               thread_nvy2_sum += block-> data[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY);
               thread_nvz2_sum += block-> data[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ);
            }
            //Multiply by volume:
            thread_nvx2_sum = thread_nvx2_sum * DV3;
            thread_nvy2_sum = thread_nvy2_sum * DV3;
            thread_nvz2_sum = thread_nvz2_sum * DV3;
         }
         
         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            Pressure += physicalconstants::MASS_PROTON * THIRD * (thread_nvx2_sum + thread_nvy2_sum + thread_nvz2_sum);
         }
      }
      const char* ptr = reinterpret_cast<const char*>(&Pressure);
      for (uint i=0; i<sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePressureBackstream::setSpatialCell(const SpatialCell* cell) {
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

   // Adding pressure calculations for backstream population to Vlasiator.
   // p_ij = m/3 * integral((v - <V>)_i(v - <V>)_j * f(r,v) dV)
   
   // Pressure tensor 6 components (11, 22, 33, 23, 13, 12) added by YK
   // Split into VariablePTensorBackstreamDiagonal (11, 22, 33)
   // and VariablePTensorOffDiagonal (23, 13, 12)
   VariablePTensorBackstreamDiagonal::VariablePTensorBackstreamDiagonal(): DataReductionOperator() { }
   VariablePTensorBackstreamDiagonal::~VariablePTensorBackstreamDiagonal() { }
   
   std::string VariablePTensorBackstreamDiagonal::getName() const {return "PTensorDiagonal";}
   
   bool VariablePTensorBackstreamDiagonal::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 3;
      return true;
   }
   
   bool VariablePTensorBackstreamDiagonal::reduceData(const SpatialCell* cell,char* buffer) {
      const Real HALF = 0.5;
      const Real THIRD = 1.0/3.0;
      # pragma omp parallel
      {
         Real thread_nvxvx_sum = 0.0;
         Real thread_nvyvy_sum = 0.0;
         Real thread_nvzvz_sum = 0.0;
         # pragma omp for
         for(uint n=0; n<cell->number_of_blocks; n++) {
            unsigned int blockId = cell->velocity_block_list[n];
            const Velocity_Block* block = cell->at(blockId); //returns a reference to block   
            const Real DV3 = block-> parameters[BlockParams::DVX] * block-> parameters[BlockParams::DVY] * block-> parameters[BlockParams::DVZ];
            vector< array<uint, 3> > vCellIndices;
            getBackstreamVelocityCellIndices(block, vCellIndices);
            for( vector< array<uint, 3> >::const_iterator it = vCellIndices.begin(); it != vCellIndices.end(); ++it ) {
               //Go through every velocity cell:
               const array<uint, 3> indices = *it;
               const uint i = indices[0];
               const uint j = indices[1];
               const uint k = indices[2];
               const Real VX = block-> parameters[BlockParams::VXCRD] + (i+HALF) * block-> parameters[BlockParams::DVX];
               const Real VY = block-> parameters[BlockParams::VYCRD] + (j+HALF) * block-> parameters[BlockParams::DVY];
               const Real VZ = block-> parameters[BlockParams::VZCRD] + (k+HALF) * block-> parameters[BlockParams::DVZ];
               thread_nvxvx_sum += block-> data[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX);
               thread_nvyvy_sum += block-> data[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY);
               thread_nvzvz_sum += block-> data[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ);
            }
            //Multiply by volume:
            thread_nvxvx_sum = thread_nvxvx_sum * DV3;
            thread_nvyvy_sum = thread_nvyvy_sum * DV3;
            thread_nvzvz_sum = thread_nvzvz_sum * DV3;
         }
         
         // Accumulate contributions coming from this velocity block to the 
         // spatial cell velocity moments. If multithreading / OpenMP is used, 
         // these updates need to be atomic:
         # pragma omp critical
         {
            PTensor[0] += physicalconstants::MASS_PROTON * THIRD * thread_nvxvx_sum;
            PTensor[1] += physicalconstants::MASS_PROTON * THIRD * thread_nvyvy_sum;
            PTensor[2] += physicalconstants::MASS_PROTON * THIRD * thread_nvzvz_sum;
         }
      }
      const char* ptr = reinterpret_cast<const char*>(&PTensor);
      for (uint i=0; i<3*sizeof(Real); ++i) buffer[i] = ptr[i];
      return true;
   }
   
   bool VariablePTensorBackstreamDiagonal::setSpatialCell(const SpatialCell* cell) {
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
   
} // namespace DRO
