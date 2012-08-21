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
#include "../project.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {
   /*!\brief Ionosphere is a class applying ionospheric boundary conditions.
    * 
    * Ionosphere is a class handling cells tagged as sysboundarytype::IONOSPHERE by this
    * system boundary condition. It applies ionospheric boundary conditions.
    */
   class Ionosphere: public SysBoundaryCondition {
   public:
      Ionosphere();
      virtual ~Ionosphere();
      
      static void addParameters();
      virtual void getParameters();
      
      virtual bool initSysBoundary(creal& t);
      virtual int assignSysBoundary(creal* cellParams);
      virtual bool applyInitialState(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid);
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      Real center[3]; /*!< Coordinates of the centre of the ionosphere. */
      Real radius; /*!< Radius of the ionosphere. */
   };
}

#endif
