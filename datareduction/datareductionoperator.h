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

#ifndef DATAREDUCTIONOPERATOR_H
#define DATAREDUCTIONOPERATOR_H

#include <vector>
#include "../definitions.h"
#include "../spatial_cell.hpp"
using namespace spatial_cell;

namespace DRO {

   /** DRO::DataReductionOperator defines a base class for reducing simulation data
    * (six-dimensional distribution function) into more compact variables, e.g. 
    * scalar fields, which can be written into file(s) and visualized.
    * 
    * The intention is that each DRO::DataReductionOperator stores the reduced data 
    * into internal variables, whose values are written into a byte array when 
    * DRO::DataReductionOperator::appendReducedData is called.
    * 
    * If needed, a user can write his or her own DRO::DataReductionOperators, which 
    * are loaded when the simulation initializes.
    */
   class DataReductionOperator {
    public:
      DataReductionOperator();
      virtual ~DataReductionOperator();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceData(const SpatialCell* cell,Real * result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
    protected:
	
   };
   
   class MPIrank: public DataReductionOperator {
    public:
      MPIrank();
      ~MPIrank();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      Real rank;
      int mpiRank;
   };
   
   class BoundaryType: public DataReductionOperator {
   public:
      BoundaryType();
      ~BoundaryType();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      int boundaryType;
   };

   class Blocks: public DataReductionOperator {
    public:
      Blocks();
      ~Blocks();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool reduceData(const SpatialCell* cell,Real* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      int nBlocks;
   };

   
   class VariableB: public DataReductionOperator {
    public:
      VariableB();
      ~VariableB();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      Real Bx;
      Real By;
      Real Bz;
      const Real* B;
   };

   class VariableVolB: public DataReductionOperator {
    public:
      VariableVolB();
      ~VariableVolB();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
            bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      const Real* B;
   };
   
   class VariableE: public DataReductionOperator {
    public:
      VariableE();
      ~VariableE();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      Real Ex;
      Real Ey;
      Real Ez;
      const Real* E;
   };

   class VariableVolE: public DataReductionOperator {
    public:
      VariableVolE();
      ~VariableVolE();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      const Real* E;
   };
   
   class VariableRho: public DataReductionOperator {
    public:
      VariableRho();
      ~VariableRho();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool reduceData(const SpatialCell* cell,Real* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      Real rho;
   };

   class VariableRhoLossAdjust: public DataReductionOperator {
    public:
      VariableRhoLossAdjust();
      ~VariableRhoLossAdjust();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool reduceData(const SpatialCell* cell,Real* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      Real rhoLoss;
   };


   class VariableRhoLossVelBoundary: public DataReductionOperator {
    public:
      VariableRhoLossVelBoundary();
      ~VariableRhoLossVelBoundary();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool reduceData(const SpatialCell* cell,Real* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      Real rhoLoss;
   };

   
   class VariableRhoV: public DataReductionOperator {
    public:
      VariableRhoV();
      ~VariableRhoV();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:
      Real rhovx;
      Real rhovy;
      Real rhovz;
      const Real* rhov;
   };

   class MaxVi: public DataReductionOperator {
   public:
      MaxVi();
      ~MaxVi();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool reduceData(const SpatialCell* cell,Real *buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
    protected:

   };
   
  // Added by YK
  class VariablePressure: public DataReductionOperator {
     public:
        VariablePressure();
        ~VariablePressure();
     
        bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
        std::string getName() const;
        bool reduceData(const SpatialCell* cell,char* buffer);
        bool setSpatialCell(const SpatialCell* cell);
     
     protected:
        Real averageVX, averageVY, averageVZ;
	Real Pressure;
   };
   
   class VariablePTensorDiagonal: public DataReductionOperator {
   public:
      VariablePTensorDiagonal();
      ~VariablePTensorDiagonal();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
   };
   
   class VariablePTensorOffDiagonal: public DataReductionOperator {
   public:
      VariablePTensorOffDiagonal();
      ~VariablePTensorOffDiagonal();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
   };
   
   class DiagnosticFluxB: public DataReductionOperator {
   public:
      DiagnosticFluxB();
      ~DiagnosticFluxB();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,Real* result);
      bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      
   };
   
   class VariabledBxdz: public DataReductionOperator {
   public:
      VariabledBxdz();
      ~VariabledBxdz();
      
      bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      std::string getName() const;
      bool reduceData(const SpatialCell* cell,char* buffer);
      bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real value;
   };
   
} // namespace DRO

#endif

