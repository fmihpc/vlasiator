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
#include "../spatial_cell.hpp"

# define DATA_LENGTH 9

namespace SBC {
   class SolarWind: public SysBoundaryCondition {
   public:
      SolarWind();
      ~SolarWind();
      
      bool initSysBoundary(creal& t);
      int assignSysBoundary(creal* cellParams);
      Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,
                                 creal& dx,creal& dy,creal& dz,
                                 creal& vx,creal& vy,creal& vz,
                                 creal& dvx,creal& dvy,creal& dvz);
      void calcCellParameters(Real* cellParams, creal& t);
      bool loadInputData();
      std::vector<std::vector<Real> > loadFile(const char* file);
      bool generateTemplateCells(creal& t);
      void generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t);
      void interpolate(const int inputDataIndex, creal t, Real* rho, Real* T, Real* vx, Real* vy, Real* vz, Real* Bx, Real* By, Real* Bz);
      void determineFace(creal x, creal y, creal z, creal dx, creal dy, creal dz);
      std::string getName() const;
      virtual uint getIndex() const;
      virtual uint getPrecedence() const;
      
   protected:
      uint faces : 6;
      uint isThisCellOnAFace : 6;
      std::vector<std::vector<Real> > inputData[6];
//      int faceToInputData[6];
      spatial_cell::SpatialCell templateCells[6];
   };
}

#endif
