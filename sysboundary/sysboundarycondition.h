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

#include <vector>
#include "../definitions.h"
#include "../spatial_cell.hpp"

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
         
         virtual bool initSysBoundary();
         virtual int assignSysBoundary(creal* cellParams);
         virtual Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,
                                            creal& dx,creal& dy,creal& dz,
                                            creal& vx,creal& vy,creal& vz,
                                            creal& dvx,creal& dvy,creal& dvz);
         virtual void calcCellParameters(Real* cellParams, creal& t);
         virtual std::string getName() const;
         virtual uint getIndex() const;
         virtual uint getPrecedence() const;
         
      protected:
         
   };
} // namespace SBC

#endif
