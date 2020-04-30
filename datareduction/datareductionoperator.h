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

#include "fsgrid.hpp"
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
      virtual bool getUnitMetadata(std::string& _unit,std::string& _unitLaTeX,std::string& _variableLaTeX,std::string& _unitConversion) {
	_unit=unit;
	_unitLaTeX=unitLaTeX;
	_unitConversion=unitConversion;
	_variableLaTeX=variableLaTeX;
	return true;
      };
      virtual bool setUnitMetadata(std::string& _unit,std::string& _unitLaTeX,std::string& _variableLaTeX,std::string& _unitConversion) {
	unit = _unit;
	unitLaTeX = _unitLaTeX;
	unitConversion = _unitConversion;
	variableLaTeX = _variableLaTeX;
	return true;
      }

      virtual std::string getName() const = 0;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool reduceDiagnostic(const SpatialCell* cell,Real * result);
      virtual bool setSpatialCell(const SpatialCell* cell) = 0;
      
   protected:
      std::string unit;
      std::string unitLaTeX;
      std::string variableLaTeX;
      std::string unitConversion;
      
   };

   class DataReductionOperatorHandlesWriting: public DataReductionOperator {
   public:
      DataReductionOperatorHandlesWriting() : DataReductionOperator() {};
      virtual bool writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,const std::string& meshName,
                             vlsv::Writer& vlsvWriter) = 0;
   };

   class DataReductionOperatorHasParameters: public DataReductionOperator {
   public:
      DataReductionOperatorHasParameters() : DataReductionOperator() {};
      virtual bool writeParameters(vlsv::Writer& vlsvWriter) = 0;
   };

   class DataReductionOperatorFsGrid : public DataReductionOperator {

      public:
        typedef std::function<std::vector<double>(
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)> ReductionLambda;
      private:
         ReductionLambda lambda;
         std::string variableName;

      public:
         DataReductionOperatorFsGrid(const std::string& name, ReductionLambda l) : DataReductionOperator(),lambda(l),variableName(name) {};
	 virtual std::string getName() const;
	 virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
	 virtual bool setSpatialCell(const SpatialCell* cell);
	 virtual bool reduceData(const SpatialCell* cell,char* buffer);
	 virtual bool reduceDiagnostic(const SpatialCell* cell,Real * result);
         virtual bool writeFsGridData(
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
                      const bool writeAsFloat=false);
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
   
   class VariableRhoThermal: public DataReductionOperator {
   public:
      VariableRhoThermal(cuint popID);
      virtual ~VariableRhoThermal();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
      Real RhoThermal;
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariableRhoNonthermal: public DataReductionOperator {
   public:
      VariableRhoNonthermal(cuint popID);
      virtual ~VariableRhoNonthermal();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
      Real RhoNonthermal;
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariableVThermal: public DataReductionOperator {
   public:
      VariableVThermal(cuint popID);
      virtual ~VariableVThermal();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
     
   protected:
      Real VThermal[3];
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariableVNonthermal: public DataReductionOperator {
   public:
      VariableVNonthermal(cuint popID);
      virtual ~VariableVNonthermal();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);

   protected:
      Real VNonthermal[3];
      uint popID;
      std::string popName;
      bool doSkip;
   };

   class VariablePTensorThermalDiagonal: public DataReductionOperator {
   public:
      VariablePTensorThermalDiagonal(cuint popID);
      virtual ~VariablePTensorThermalDiagonal();
      
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

   class VariablePTensorNonthermalDiagonal: public DataReductionOperator {
   public:
      VariablePTensorNonthermalDiagonal(cuint popID);
      virtual ~VariablePTensorNonthermalDiagonal();

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

   class VariablePTensorThermalOffDiagonal: public DataReductionOperator {
   public:
      VariablePTensorThermalOffDiagonal(cuint popID);
      virtual ~VariablePTensorThermalOffDiagonal();

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

   class VariablePTensorNonthermalOffDiagonal: public DataReductionOperator {
   public:
      VariablePTensorNonthermalOffDiagonal(cuint popID);
      virtual ~VariablePTensorNonthermalOffDiagonal();

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

   class VariableEnergyDensity: public DataReductionOperatorHasParameters {
   public:
      VariableEnergyDensity(cuint popID);
      virtual ~VariableEnergyDensity();

      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      virtual bool writeParameters(vlsv::Writer& vlsvWriter);

   protected:
      uint popID;
      std::string popName;
      Real EDensity[3];
      Real solarwindenergy;
      Real E1limit;
      Real E2limit;
   };
   
   // Precipitation directional differential number flux
   class VariablePrecipitationDiffFlux: public DataReductionOperatorHasParameters {
   public:
      VariablePrecipitationDiffFlux(cuint popID);
      virtual ~VariablePrecipitationDiffFlux();
      
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool reduceData(const SpatialCell* cell,char* buffer);
      virtual bool setSpatialCell(const SpatialCell* cell);
      virtual bool writeParameters(vlsv::Writer& vlsvWriter);
      
   protected:
      uint popID;
      std::string popName;
      int nChannels;
      Real emin, emax;
      Real lossConeAngle;
      std::vector<Real> channels, dataDiffFlux;
   };
} // namespace DRO

#endif

