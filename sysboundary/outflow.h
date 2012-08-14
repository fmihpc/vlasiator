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

#ifndef OUTFLOW_H
#define OUTFLOW_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {
   class Outflow: public SysBoundaryCondition {
   public:
      Outflow();
      ~Outflow();
      
      void getParameters();
      
      bool initSysBoundary(creal& t);
      int assignSysBoundary(creal* cellParams);
      Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,
                                 creal& dx,creal& dy,creal& dz,
                                 creal& vx,creal& vy,creal& vz,
                                 creal& dvx,creal& dvy,creal& dvz);
      void calcCellParameters(Real* cellParams, creal& t);
      std::string getName() const;
      virtual uint getIndex() const;
      virtual uint getPrecedence() const;
      virtual bool isDynamic() const;
      
   protected:
      uint faces : 6;
      uint isThisCellOnAFace : 6;
      std::vector<std::string> faceList; /*!< List of faces on which outflow boundary conditions are to be applied ([+-][xyz]). */
      uint precedence; /*!< Precedence value of the outflow system boundary condition. */
   };
}

#endif
