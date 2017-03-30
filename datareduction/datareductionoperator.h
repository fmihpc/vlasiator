/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool handlesWriting() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceData(const SpatialCell* cell,Real * result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      virtual bool writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,const std::string& meshName,
                             vlsv::Writer& vlsvWriter);
      
   protected:
   
   };

   class DataReductionOperatorCellParams: public DataReductionOperator {
   public:
      DataReductionOperatorCellParams(const std::string& name,const unsigned int parameterIndex,const unsigned int vectorSize);
      virtual ~DataReductionOperatorCellParams();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceData(const SpatialCell* cell,Real * result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      uint _parameterIndex;
      uint _vectorSize;
      std::string _name;
      const Real *_data;
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

   class BoundaryLayerNew: public DataReductionOperator {
   public:
      BoundaryLayerNew();
      virtual ~BoundaryLayerNew();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      int boundaryLayer;
   };

   class VelocitySubSteps: public DataReductionOperator {
   public:
      VelocitySubSteps();
      virtual ~VelocitySubSteps();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      int substeps;
   };
   

   class Blocks: public DataReductionOperator {
   public:
      Blocks();
      virtual ~Blocks();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceData(const SpatialCell* cell,Real* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      int nBlocks;
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

   
  // Added by YK
   class VariablePressure: public DataReductionOperator {
   public:
      VariablePressure();
      virtual ~VariablePressure();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceData(const SpatialCell* cell,Real * result);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
      Real averageVX, averageVY, averageVZ;
      Real Pressure;
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
      VariablePTensorDiagonal();
      virtual ~VariablePTensorDiagonal();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
   };
   
   class VariablePTensorOffDiagonal: public DataReductionOperator {
   public:
      VariablePTensorOffDiagonal();
      virtual ~VariablePTensorOffDiagonal();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
   };
   
   class DiagnosticFluxB: public DataReductionOperator {
   public:
      DiagnosticFluxB();
      virtual  ~DiagnosticFluxB();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,Real* result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      
   };
   
   class DiagnosticFluxE: public DataReductionOperator {
   public:
      DiagnosticFluxE();
      virtual  ~DiagnosticFluxE();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,Real* result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      
   };
   
   class MaxDistributionFunction: public DataReductionOperator {
   public:
      MaxDistributionFunction();
      virtual ~MaxDistributionFunction();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceData(const SpatialCell* cell,Real *buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real maxF;
   };
   
   class MinDistributionFunction: public DataReductionOperator {
   public:
      MinDistributionFunction();
      virtual ~MinDistributionFunction();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceData(const SpatialCell* cell,Real *buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real minF;
   };

   /** This class writes all scalar and two- or three-component vector data 
    * that is stored to MeshDataContainer to output file.*/
   class VariableMeshData: public DataReductionOperator {
   public:
      VariableMeshData();
      virtual ~VariableMeshData();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool handlesWriting() const;
      virtual bool setSpatialCell(const SpatialCell* cell);
      virtual bool writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,const std::string& meshName,
                             vlsv::Writer& vlsvWriter);
      
   private:
      
   };
   
   class VariableRhomBackstream: public DataReductionOperator {
   public:
      VariableRhomBackstream();
      virtual ~VariableRhomBackstream();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
   Real RhomBackstream;
   };

   class VariableRhomNonBackstream: public DataReductionOperator {
   public:
      VariableRhomNonBackstream();
      virtual ~VariableRhomNonBackstream();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
   Real Rhom;
   };


   class VariableRhomVBackstream: public DataReductionOperator {
   public:
      VariableRhomVBackstream();
      virtual ~VariableRhomVBackstream();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
      Real RhomVBackstream[3];
   };

   class VariableRhomVNonBackstream: public DataReductionOperator {
   public:
      VariableRhomVNonBackstream();
      virtual ~VariableRhomVNonBackstream();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real RhomV[3];
   };

   class VariablePressureBackstream: public DataReductionOperator {
   public:
      VariablePressureBackstream();
      virtual ~VariablePressureBackstream();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real averageVX, averageVY, averageVZ;
      Real Pressure;
   };

   class VariablePressureNonBackstream: public DataReductionOperator {
   public:
      VariablePressureNonBackstream();
      virtual ~VariablePressureNonBackstream();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real averageVX, averageVY, averageVZ;
      Real Pressure;
   };

   class VariablePTensorBackstreamDiagonal: public DataReductionOperator {
   public:
      VariablePTensorBackstreamDiagonal();
      virtual ~VariablePTensorBackstreamDiagonal();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      
   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
   };

   class VariablePTensorNonBackstreamDiagonal: public DataReductionOperator {
   public:
      VariablePTensorNonBackstreamDiagonal();
      virtual ~VariablePTensorNonBackstreamDiagonal();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
   };

   class VariablePTensorBackstreamOffDiagonal: public DataReductionOperator {
   public:
      VariablePTensorBackstreamOffDiagonal();
      virtual ~VariablePTensorBackstreamOffDiagonal();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
   };

   class VariablePTensorNonBackstreamOffDiagonal: public DataReductionOperator {
   public:
      VariablePTensorNonBackstreamOffDiagonal();
      virtual ~VariablePTensorNonBackstreamOffDiagonal();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real averageVX, averageVY, averageVZ;
      Real PTensor[3];
   };
   
   class VariableMinValue: public DataReductionOperator {
   public:
      VariableMinValue();
      virtual ~VariableMinValue();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool handlesWriting() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceData(const spatial_cell::SpatialCell* cell,Real* result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      virtual bool writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,const std::string& meshName,vlsv::Writer& vlsvWriter);
      
   protected:
      
   };
   
} // namespace DRO

#endif

