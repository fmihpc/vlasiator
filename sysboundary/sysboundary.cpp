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

#include <cstdlib>
#include <iostream>

#include "sysboundary.h"

using namespace std;

bool precedenceSort(const SBC::SysBoundaryCondition* first,
                    const SBC::SysBoundaryCondition* second) {
   if(first->getPrecedence() < second->getPrecedence()) return true;
   else return false;
}

// ************************************************************
// ***** DEFINITIONS FOR BOUNDARY CLASS *****
// ************************************************************

/*! Constructor for class SysBoundary. Increases the value of SysBoundary::nSysBoundaries by one. */
SysBoundary::SysBoundary() { }

/*! Destructor for class SysBoundary. Reduces the value of SysBoundary::nSysBoundaries by one,
 * and if after the destruction SysBoundary::nSysBoundaries equals zero all stored SysBoundaries are deleted.
 */
SysBoundary::~SysBoundary() {
   // Call delete for each SysBoundaryCondition:
   for (list<SBC::SysBoundaryCondition*>::iterator it=sysBoundaries.begin();
        it!=sysBoundaries.end();
        ++it) {
      delete *it;
      *it = NULL;
   }
}

void SysBoundary::addParameters() {
   Readparameters::addComposing("boundaries.boundary", "List of boundary condition (BC) types to be used. Each boundary condition to be used has to be on a new line boundary = YYY. Available (20120807) are outflow ionosphere solarwind.");
   Readparameters::add("boundaries.periodic_x","If 'yes' the grid is periodic in x-direction. Defaults to 'no'.","no");
   Readparameters::add("boundaries.periodic_y","If 'yes' the grid is periodic in y-direction. Defaults to 'no'.","no");
   Readparameters::add("boundaries.periodic_z","If 'yes' the grid is periodic in z-direction. Defaults to 'no'.","no");
   
   //call static addParameter functions in all bc's
   SBC::DoNotCompute::addParameters();
   SBC::Ionosphere::addParameters();
   SBC::Outflow::addParameters();
   SBC::SolarWind::addParameters();
}

void SysBoundary::getParameters() {
   Readparameters::get("boundaries.boundary", sysBoundaryCondList);
   std::string periodic_x,periodic_y,periodic_z;
   Readparameters::get("gridbuilder.periodic_x",periodic_x);
   Readparameters::get("gridbuilder.periodic_y",periodic_y);
   Readparameters::get("gridbuilder.periodic_z",periodic_z);
   isPeriodic[0] = false;
   isPeriodic[1] = false;
   isPeriodic[2] = false;
   if (periodic_x == "yes") isPeriodic[0] = true;
   if (periodic_y == "yes") isPeriodic[1] = true;
   if (periodic_z == "yes") isPeriodic[2] = true;
   //getParameters for each system boundary condition is called initialization.
}

/*! Add a new SBC::SysBoundaryCondition which has been created with new sysBoundary. 
 * SysBoundary will take care of deleting it.
 * @return If true, the given SBC::SysBoundaryCondition was added successfully.
 */
bool SysBoundary::addSysBoundary(SBC::SysBoundaryCondition* bc, creal& t) {
   sysBoundaries.push_back(bc);
   if(sysBoundaries.size() > 1) {
      sysBoundaries.sort(precedenceSort);
   }
   
   bool success = true;
   if(bc->initSysBoundary(t) == false) {
      success = false;
   }
   
   // This assumes that only one instance of each type is created.
   indexToSysBoundary[bc->getIndex()] = bc;
   
   return success;
}

bool SysBoundary::initSysBoundaries(creal& t) {
   bool success = true;
   vector<string>::const_iterator it;
   for (it = sysBoundaryCondList.begin();
        it != sysBoundaryCondList.end();
   it++) {
      if(*it == "Outflow") {
         if(this->addSysBoundary(new SBC::Outflow, t) == false) {
            cerr << "Error in adding Outflow boundary." << endl;
            success = false;
         }
         isThisDynamic = isThisDynamic|this->getSysBoundary(sysboundarytype::OUTFLOW)->isDynamic();
      }
      if(*it == "Ionosphere") {
         if(this->addSysBoundary(new SBC::Ionosphere, t) == false) {
            cerr << "Error in adding Ionosphere boundary." << endl;
            success = false;
         }
         if(this->addSysBoundary(new SBC::DoNotCompute, t) == false) {
            cerr << "Error in adding DoNotCompute boundary (for Ionosphere)." << endl;
            success = false;
         }
         isThisDynamic = isThisDynamic|
         this->getSysBoundary(sysboundarytype::IONOSPHERE)->isDynamic();
      }
      if(*it == "SolarWind") {
         if(this->addSysBoundary(new SBC::SolarWind, t) == false) {
            cerr << "Error in adding SolarWind boundary." << endl;
            success = false;
         }
         isThisDynamic = isThisDynamic|
         this->getSysBoundary(sysboundarytype::SW)->isDynamic();
      }
   }
   
   return success;
}

bool SysBoundary::classifyCells(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   using namespace sysboundarytype;
   uint indexToAssign, tmpType;
   
   vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i=0; i<cells.size(); ++i) {
      indexToAssign = NOT_SYSBOUNDARY;
      list<SBC::SysBoundaryCondition*>::iterator it;
      for (it = sysBoundaries.begin();
         it != sysBoundaries.end();
         it++) {
         tmpType=(*it)->assignSysBoundary(&(mpiGrid[cells[i]]->parameters[0]));
         
         if(tmpType == DO_NOT_COMPUTE) {
            indexToAssign = tmpType;
            break; 
         } else if (tmpType != NOT_SYSBOUNDARY) {
            indexToAssign = tmpType;
         }
      }
      mpiGrid[cells[i]]->sysBoundaryFlag = indexToAssign;
   }
   
   return true;
}

bool SysBoundary::applyInitialState(dccrg::Dccrg<SpatialCell>& mpiGrid) {
   bool success = true;
   using namespace sysboundarytype;
   
   list<SBC::SysBoundaryCondition*>::iterator it;
   for (it = sysBoundaries.begin();
        it != sysBoundaries.end();
   it++) {
      if((*it)->applyInitialState(mpiGrid) == false) {
         cerr << "ERROR: " << (*it)->getName() << " system boundary condition not applied correctly." << endl;
         success = false;
      }
   }
   
   return success;
}

/*! Get a pointer to the SysBoundaryCondition of given index.
 * \retval Pointer to the instance of the SysBoundaryCondition.
 */
SBC::SysBoundaryCondition* SysBoundary::getSysBoundary(uint sysBoundaryType) const {
   return indexToSysBoundary.find(sysBoundaryType)->second;
}

/*! Get the number of SysBoundaryConditions stored in SysBoundary.
 * \retval Number of SysBoundaryConditions stored in SysBoundary.
 */
unsigned int SysBoundary::size() const {return sysBoundaries.size();}
bool SysBoundary::isDynamic() const {return isThisDynamic;}
bool SysBoundary::isBoundaryPeriodic(uint direction) const {return isPeriodic[direction];}
