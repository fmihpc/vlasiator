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

#ifndef SYSBOUNDARY_H
#define SYSBOUNDARY_H

#include <vector>

#include <dccrg.hpp>
#include "../spatial_cell.hpp"
using namespace spatial_cell;
#include "sysboundarycondition.h"
#include "ionosphere.h"
#include "outflow.h"
#include "solarwind.h"


/*! The purpose of SysBoundary is to contain SBC::SysBoundaryConditions, and apply 
 * them to the simulation volume cells if the cells pertain to a specific boundary type.
 * If the simulation domain is not fully periodic then the behaviour at the edges or boundaries of the volume has to be properly defined.
 * TODO detail the workings when coded
 * 
 * If needed, a user can write his or her own SBC::SysBoundaryConditions, which 
 * are loaded when the simulation initializes.
 */
class SysBoundary {
   public:
      SysBoundary();
      ~SysBoundary();
      
      bool addSysBoundary(SBC::SysBoundaryCondition* sbc);
      std::string getName(const unsigned int& sysBoundaryID) const;
//    bool (const unsigned int& sysBoundaryID) const;
      unsigned int size() const;
      std::vector<SBC::SysBoundaryCondition*> getSysBoundariesList() const;
   
 private:
   /*! Private copy-constructor to prevent copying the class. */
   SysBoundary(const SysBoundary& bc);
   
   std::vector<SBC::SysBoundaryCondition*> sysBoundaries;
   /*!< A container for all SBC::SysBoundaryConditions stored in SysBoundary.*/
};

bool initializeSysBoundaries(SysBoundary * sbc);
bool assignSysBoundaryType(SysBoundary* sbc, dccrg::Dccrg<SpatialCell>& mpiGrid);

#endif
