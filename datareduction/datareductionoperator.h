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

#ifndef DATAREDUCTIONOPERATOR_H
#define DATAREDUCTIONOPERATOR_H

#include <vector>

#include <vlsv_writer.h>
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include "../definitions.h"
#include "../spatial_cell.hpp"
#include "../parameters.h"
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
    *
    * Datareduction oeprators are not thread-safe, some of the more intensive ones are threaded within. 
    */

   class DataReductionOperator {
   public:
      DataReductionOperator();
      virtual ~DataReductionOperator();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const = 0;
      virtual std::string getName() const = 0;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceDiagnostic(const SpatialCell* cell,Real * result);
      virtual bool setSpatialCell(const SpatialCell* cell) = 0;
      
   protected:
   
   };

   class DataReductionOperatorHandlesWriting: public DataReductionOperator {
   public:
      DataReductionOperatorHandlesWriting() : DataReductionOperator() {};
      virtual bool writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,const std::string& meshName,
                             vlsv::Writer& vlsvWriter) = 0;
   };

   class DataReductionOperatorCellParams: public DataReductionOperator {
   public:
      DataReductionOperatorCellParams(const std::string& name,const unsigned int parameterIndex,const unsigned int vectorSize);
      virtual ~DataReductionOperatorCellParams();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceDiagnostic(const SpatialCell* cell,Real * result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      uint _parameterIndex;
      uint vectorSize;
      std::string variableName;
      const Real *data;
   };

   class DataReductionOperatorDerivatives: public DataReductionOperatorCellParams {
   public:
      DataReductionOperatorDerivatives(const std::string& name,const unsigned int parameterIndex,const unsigned int vectorSize);
      virtual bool setSpatialCell(const SpatialCell* cell);
   };
   
   class DataReductionOperatorBVOLDerivatives: public DataReductionOperatorCellParams {
   public:
      DataReductionOperatorBVOLDerivatives(const std::string& name,const unsigned int parameterIndex,const unsigned int vectorSize);
      virtual bool setSpatialCell(const SpatialCell* cell);
   };
   
   class MPIrank: public DataReductionOperator {
   public:
      MPIrank();
      virtual ~MPIrank();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real rank;
      int mpiRank;
   };

   class FsGridRank: public DataReductionOperator {
   public:
      FsGridRank();
      virtual ~FsGridRank();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real rank;
      int fsgridRank;
   };

   class FsGridBoundaryType: public DataReductionOperator {
   public:
      FsGridBoundaryType();
      virtual ~FsGridBoundaryType();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real rank;
      int fsgridBoundaryType;
   };
   
   class BoundaryType: public DataReductionOperator {
   public:
      BoundaryType();
      virtual ~BoundaryType();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      int boundaryType;
   };

   class BoundaryLayer: public DataReductionOperator {
   public:
      BoundaryLayer();
      virtual ~BoundaryLayer();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      int boundaryLayer;
   };

   class Blocks: public DataReductionOperator {
   public:
      Blocks(cuint popID);
      virtual ~Blocks();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceDiagnostic(const SpatialCell* cell,Real* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      uint nBlocks;
      uint popID;
      std::string popName;
   };
   
   class VariableB: public DataReductionOperator {
   public:
      VariableB();
      virtual ~VariableB();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real B[3];
   };

      
   class VariableBVol: public DataReductionOperator {
   public:
      VariableBVol();
      virtual ~VariableBVol();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real B[3];
   };

   
   class VariablePressureSolver: public DataReductionOperator {
   public:
      VariablePressureSolver();
      virtual ~VariablePressureSolver();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real Pressure;
   };
   
   class VariablePTensorDiagonal: public DataReductionOperator {
   public:
      VariablePTensorDiagonal(cuint popID);
      virtual ~VariablePTensorDiagonal();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
      uint popID;
      std::string popName;
   };
   
   class VariablePTensorOffDiagonal: public DataReductionOperator {
   public:
      VariablePTensorOffDiagonal(cuint popID);
      virtual ~VariablePTensorOffDiagonal();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
      uint popID;
      std::string popName;
   };
   
   class DiagnosticFluxB: public DataReductionOperator {
   public:
      DiagnosticFluxB();
      virtual  ~DiagnosticFluxB();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceDiagnostic(const SpatialCell* cell,Real* result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      
   };
   
   class DiagnosticFluxE: public DataReductionOperator {
   public:
      DiagnosticFluxE();
      virtual  ~DiagnosticFluxE();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceDiagnostic(const SpatialCell* cell,Real* result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      
   };
   
   class MaxDistributionFunction: public DataReductionOperator {
   public:
      MaxDistributionFunction(cuint popID);
      virtual ~MaxDistributionFunction();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceDiagnostic(const SpatialCell* cell,Real *buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real maxF;
      uint popID;
      std::string popName;
   };
   
   class MinDistributionFunction: public DataReductionOperator {
   public:
      MinDistributionFunction(cuint popID);
      virtual ~MinDistributionFunction();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceDiagnostic(const SpatialCell* cell,Real *buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real minF;
      uint popID;
      std::string popName;
   };

   /** This class writes all scalar and two- or three-component vector data 
    * that is stored to MeshDataContainer to output file.*/
   class VariableMeshData: public DataReductionOperatorHandlesWriting {
   public:
      VariableMeshData();
      virtual ~VariableMeshData();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool setSpatialCell(const SpatialCell* cell);
      virtual bool writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,const std::string& meshName,
                             vlsv::Writer& vlsvWriter);
      
   private:
      
   };
   
   class VariableRhoBackstream: public DataReductionOperator {
   public:
      VariableRhoBackstream(cuint popID);
      virtual ~VariableRhoBackstream();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
      Real RhoBackstream;
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariableRhoNonBackstream: public DataReductionOperator {
   public:
      VariableRhoNonBackstream(cuint popID);
      virtual ~VariableRhoNonBackstream();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
      Real Rho;
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariableVBackstream: public DataReductionOperator {
   public:
      VariableVBackstream(cuint popID);
      virtual ~VariableVBackstream();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
      Real VBackstream[3];
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariableVNonBackstream: public DataReductionOperator {
   public:
      VariableVNonBackstream(cuint popID);
      virtual ~VariableVNonBackstream();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real V[3];
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariablePTensorBackstreamDiagonal: public DataReductionOperator {
   public:
      VariablePTensorBackstreamDiagonal(cuint popID);
      virtual ~VariablePTensorBackstreamDiagonal();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariablePTensorNonBackstreamDiagonal: public DataReductionOperator {
   public:
      VariablePTensorNonBackstreamDiagonal(cuint popID);
      virtual ~VariablePTensorNonBackstreamDiagonal();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariablePTensorBackstreamOffDiagonal: public DataReductionOperator {
   public:
      VariablePTensorBackstreamOffDiagonal(cuint popID);
      virtual ~VariablePTensorBackstreamOffDiagonal();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariablePTensorNonBackstreamOffDiagonal: public DataReductionOperator {
   public:
      VariablePTensorNonBackstreamOffDiagonal(cuint popID);
      virtual ~VariablePTensorNonBackstreamOffDiagonal();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
      uint popID;
      std::string popName;
      bool doSkip;
   };
   
   class VariableEffectiveSparsityThreshold: public DataReductionOperator {
   public:
      VariableEffectiveSparsityThreshold(cuint popID);
      virtual ~VariableEffectiveSparsityThreshold();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceDiagnostic(const spatial_cell::SpatialCell* cell,Real* result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      uint popID;
      std::string popName;
   };
   
} // namespace DRO

#endif

