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
#include "../project.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {
   /*!\brief Outflow is a class applying outflow boundary conditions.
    * 
    * Outflow is a class handling cells tagged as sysboundarytype::OUTFLOW by this system
    * boundary condition. It applies outflow boundary conditions.
    */
   class Outflow: public SysBoundaryCondition {
   public:
      Outflow();
      virtual ~Outflow();
      
      static void addParameters();
      virtual void getParameters();
      
      virtual bool initSysBoundary(creal& t);
      virtual int assignSysBoundary(creal* cellParams);
      virtual bool applyInitialState(dccrg::Dccrg<SpatialCell>& mpiGrid);
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      /*! Array of bool telling which faces are going to be processed by the system boundary condition.*/
      bool facesToProcess[6];
      /*! Array of bool used to tell on which face(s) (if any) a given cell is. \sa determineFace */
      bool isThisCellOnAFace[6];
      /*! List of faces on which outflow boundary conditions are to be applied ([xyz][+-]). */
      std::vector<std::string> faceList;
   };
}

#endif
