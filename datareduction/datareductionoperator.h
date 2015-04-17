/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute


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

   class VariableRhoBackstream: public DataReductionOperator {
   public:
      VariableRhoBackstream();
      virtual ~VariableRhoBackstream();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
   Real RhoBackstream;
   };

   class VariableRhoNonBackstream: public DataReductionOperator {
   public:
      VariableRhoNonBackstream();
      virtual ~VariableRhoNonBackstream();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
   Real Rho;
   };


   class VariableRhoVBackstream: public DataReductionOperator {
   public:
      VariableRhoVBackstream();
      virtual ~VariableRhoVBackstream();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
      Real RhoVBackstream[3];
   };

   class VariableRhoVNonBackstream: public DataReductionOperator {
   public:
      VariableRhoVNonBackstream();
      virtual ~VariableRhoVNonBackstream();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real RhoV[3];
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
      virtual bool handlesWriting();
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceData(const spatial_cell::SpatialCell* cell,Real* result);
      virtual bool setSpatialCell(const SpatialCell* cell);
      virtual bool writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,const std::string& meshName,vlsv::Writer& vlsvWriter);
      
   protected:
      
   };
   
} // namespace DRO

#endif

