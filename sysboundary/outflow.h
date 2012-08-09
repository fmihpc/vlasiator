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
#include "../spatial_cell.hpp"

namespace SBC {
   class Outflow: public SysBoundaryCondition {
   public:
      Outflow();
      ~Outflow();
      
      bool initSysBoundary();
      int assignSysBoundary(creal* cellParams);
      std::string getName() const;
      virtual int getIndex() const;
      
   protected:
      uint faces : 6;
      uint isThisCellOnAFace : 6;
   };
}

#endif
