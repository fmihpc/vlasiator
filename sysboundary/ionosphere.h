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

#ifndef IONOSPHERE_H
#define IONOSPHERE_H

#include <vector>
#include "../definitions.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {
   class Ionosphere: public SysBoundaryCondition {
   public:
      Ionosphere();
      ~Ionosphere();
      
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
      Real center[3]; /*!< Coordinates of the centre of the ionosphere. */
      Real radius; /*!< Radius of the ionosphere. */
      uint precedence; /*! Precedence value of the ionosphere system boundary condition. */
   };
}

#endif
