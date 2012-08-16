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

#ifndef SOLARWIND_H
#define SOLARWIND_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

# define DATA_LENGTH 9

namespace SBC {
   class SolarWind: public SysBoundaryCondition {
   public:
      SolarWind();
      ~SolarWind();
      
      void getParameters();
      
      bool initSysBoundary(creal& t);
      int assignSysBoundary(creal* cellParams);
      virtual bool applyInitialState(dccrg::Dccrg<SpatialCell>& mpiGrid);
      bool loadInputData();
      std::vector<std::vector<Real> > loadFile(const char* file);
      bool generateTemplateCells(creal& t);
      void generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t);
      void interpolate(const int inputDataIndex, creal t, Real* rho, Real* T, Real* vx, Real* vy, Real* vz, Real* Bx, Real* By, Real* Bz);
      void determineFace(creal x, creal y, creal z, creal dx, creal dy, creal dz);
      
      std::string getName() const;
      virtual uint getIndex() const;
      virtual uint getPrecedence() const;
      virtual bool isDynamic() const;
      
   protected:
      uint faces : 6;
      uint isThisCellOnAFace : 6;
      std::vector<std::vector<Real> > inputData[6];
      spatial_cell::SpatialCell templateCells[6];
      std::vector<std::string> faceList; /*!< List of faces on which solar wind boundary conditions are to be applied ([+-][xyz]). */
      std::string files[6]; /*!< Input files for the solar wind boundary conditions. */
      bool isThisDynamic; /*!< Is the solar wind inflow dynamic in time or not. */
      uint precedence; /*! Precedence value of the solar wind system boundary condition. */
   };
}

#endif
