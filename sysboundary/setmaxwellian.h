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

#ifndef SETMAXWELLIAN_H
#define SETMAXWELLIAN_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"
#include "setbyuser.h"

namespace SBC {
   /*!\brief SetMaxwellian is a class applying fixed Maxwellian conditions according to parameters read from an input file.
    * 
    * Maxwellian is a class handling cells tagged as sysboundarytype::MAXWELLIAN by this
    * system boundary condition.
    * 
    * It applies fixed Maxwellian settings to the system boundary cells, the parameters of
    * which are being read from an input file.
    * 
    * The class inherits most of its machinery from
    * SysBoundaryCondition::SetByUser. The parameters are more general than for Maxwellian
    * and could be put in SysBoundaryCondition::SetByUser but this way they can have a
    * specific prefix which is needed if several inheriting classes are needed.
    * 
    */
   class SetMaxwellian: public SetByUser {
   public:
      SetMaxwellian();
      virtual ~SetMaxwellian();
      
      static void addParameters();
      virtual void getParameters();
      
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      void generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t);
      
      uint nSpaceSamples;
      uint nVelocitySamples;
   };
}

#endif
