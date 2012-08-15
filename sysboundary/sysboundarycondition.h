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

#ifndef SYSBOUNDARYCONDITION_H
#define SYSBOUNDARYCONDITION_H

#include <dccrg.hpp>
#include <vector>
#include "../definitions.h"
#include "../spatial_cell.hpp"

using namespace spatial_cell;

namespace SBC {
   /*! SBC::SysBoundaryCondition defines a base class for applying boundary conditions.
    * 
    * TODO detail what functions are included
    * 
    * If needed, a user can write his or her own SBC::SysBoundaryConditions, which 
    * are loaded when the simulation initializes.
    */
   class SysBoundaryCondition {
      public:
         SysBoundaryCondition();
         virtual ~SysBoundaryCondition();
         
         void getParameters();
         
         virtual bool initSysBoundary(creal& t);
         virtual int assignSysBoundary(creal* cellParams);
         virtual bool applyInitialState(dccrg::Dccrg<SpatialCell>& mpiGrid);
         virtual std::string getName() const;
         virtual uint getIndex() const;
         virtual uint getPrecedence() const;
         virtual bool isDynamic() const;
      
      protected:
         
   };
} // namespace SBC

#endif
