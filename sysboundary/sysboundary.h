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
      
      void addParameters();
      void getParameters();
      
      bool addSysBoundary(SBC::SysBoundaryCondition* sbc, creal& t);
      bool initializeSysBoundaries(creal& t);
      bool assignSysBoundaryType(SpatialCell& cell);
      unsigned int size() const;
      SBC::SysBoundaryCondition* getSysBoundary(uint sysBoundaryType) const;
      bool isDynamic() const;
      bool isBoundaryPeriodic(uint direction) const;
   
   private:
      /*! Private copy-constructor to prevent copying the class. */
      SysBoundary(const SysBoundary& bc);
   
      /*! A container for all SBC::SysBoundaryConditions stored in SysBoundary.*/
      std::list<SBC::SysBoundaryCondition*> sysBoundaries;
      /*! A map from the system boundary types to the corresponding class member. */
      std::map<uint, SBC::SysBoundaryCondition*> indexToSysBoundary;
      std::vector<std::string> sysBoundaryCondList; /*!< List of system boundary conditions (SBC) to be used. */
      bool isThisDynamic;
      bool isPeriodic[3];
      
   //Add getParameters, readParameters to read in parameters from cfg,
   //this is called from main() where we also add project
   //parameters. getParameters and readParameters again calls get and
   //read functions in each boundary conditions.

   //Add functions below into class


};


Real calcSysBoundaryPhaseSpaceDensity(SysBoundary& sysBoundaries,
                                      uint sysBoundaryFlag,
                                      creal& x, creal& y, creal& z,
                                      creal& dx, creal& dy, creal& dz,
                                      creal& vx, creal& vy, creal& vz,
                                      creal& dvx, creal& dvy, creal& dvz);
void calcSysBoundaryCellParameters(SysBoundary& sysBoundaries,
                                   uint sysBoundaryFlag,
                                   Real* cellParams,
                                   creal& t);

bool precedenceSort(const SBC::SysBoundaryCondition* first, 
                    const SBC::SysBoundaryCondition* second);

#endif
