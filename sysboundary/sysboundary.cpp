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

#include "../parameters.h"
#include "sysboundary.h"
#include "ionosphere.h"
#include "outflow.h"
#include "solarwind.h"

using namespace std;

bool initializeSysBoundaries(SysBoundary * bc, creal& t)
{
   bool isSysBoundaryCondDynamic = false;
   vector<string>::const_iterator it;
   for (it = Parameters::sysBoundaryCondList.begin();
        it != Parameters::sysBoundaryCondList.end();
        it++) {
      if(*it == "Outflow") {
         bc->addSysBoundary(new SBC::Outflow, t);
      }
      if(*it == "Ionosphere") {
         bc->addSysBoundary(new SBC::Ionosphere, t);
      }
      if(*it == "SolarWind") {
         bc->addSysBoundary(new SBC::SolarWind, t);
         isSysBoundaryCondDynamic = Parameters::isSolarWindDynamic;
      }
   }
   return isSysBoundaryCondDynamic;
}

bool assignSysBoundaryType(SysBoundary * bc, SpatialCell& cell)
{
   using namespace sysboundarytype;
   
   uint indexToAssign, tmpType;
   vector<SBC::SysBoundaryCondition*> sysBoundariesList = bc->getSysBoundariesList();
   map<uint, uint> precedenceMap = bc->getPrecedenceMap();
   
   indexToAssign = NOT_SYSBOUNDARY;
   for (uint j = 0;
         j < sysBoundariesList.size();
         j++) {
      if((tmpType=sysBoundariesList[j]->assignSysBoundary(&(cell.parameters[0]))) == DO_NOT_COMPUTE) {
         indexToAssign = tmpType;
         break;
      } else if (precedenceMap.find(tmpType)->second > precedenceMap.find(indexToAssign)->second) {
         indexToAssign = tmpType;
      }
   }
   cell.sysBoundaryFlag = 1.0*indexToAssign;
   
   return true;
}

Real calcSysBoundaryPhaseSpaceDensity(SysBoundary* sysBoundaries,
                                      uint sysBoundaryFlag,
                                      creal& x,creal& y,creal& z,
                                      creal& dx,creal& dy,creal& dz,
                                      creal& vx,creal& vy,creal& vz,
                                      creal& dvx,creal& dvy,creal& dvz
) {
   Real average;
   
   using namespace sysboundarytype;
   
   switch(sysBoundaryFlag) {
      case DO_NOT_COMPUTE:
      case NOT_SYSBOUNDARY:
         cerr << "ERROR: called SysBoundary::calcPhaseSpaceDensity for a cell flagged DO_NOT_COMPUTE or NOT_SYSBOUNDARY!!" << endl;
         exit(1);
         break;
      case OUTFLOW:
         average = sysBoundaries->getSysBoundary(OUTFLOW)->calcPhaseSpaceDensity(x,y,z,
                                                                                 dx,dy,dz,
                                                                                 vx,vy,vz,
                                                                                 dvx,dvy,dvz);
         break;
      case IONOSPHERE:
         average = sysBoundaries->getSysBoundary(IONOSPHERE)->calcPhaseSpaceDensity(x,y,z,
                                                                                    dx,dy,dz,
                                                                                    vx,vy,vz,
                                                                                    dvx,dvy,dvz);
         break;
      case SW:
         average = sysBoundaries->getSysBoundary(SW)->calcPhaseSpaceDensity(x,y,z,
                                                                            dx,dy,dz,
                                                                            vx,vy,vz,
                                                                            dvx,dvy,dvz);
         break;
      default:
         average = 0.0;
   }
   
   return average;
}

void calcSysBoundaryCellParameters(SysBoundary* sysBoundaries,
                                   uint sysBoundaryFlag,
                                   Real* cellParams,
                                   creal& t) {
   using namespace sysboundarytype;
   
   switch(sysBoundaryFlag) {
      case DO_NOT_COMPUTE:
      case NOT_SYSBOUNDARY:
         cerr << "ERROR: called SysBoundary::calcCellParameters for a cell flagged DO_NOT_COMPUTE or NOT_SYSBOUNDARY!!" << endl;
         exit(1);
         break;
      case OUTFLOW:
         sysBoundaries->getSysBoundary(OUTFLOW)->calcCellParameters(cellParams, t);
         break;
      case IONOSPHERE:
         sysBoundaries->getSysBoundary(IONOSPHERE)->calcCellParameters(cellParams, t);
         break;
      case SW:
         sysBoundaries->getSysBoundary(SW)->calcCellParameters(cellParams, t);
         break;
      default:
         break;
   }
}

// ************************************************************
// ***** DEFINITIONS FOR BOUNDARY CLASS *****
// ************************************************************

static unsigned int nSysBoundaries = 0;

/*! Constructor for class SysBoundary. Increases the value of SysBoundary::nSysBoundaries by one. */
SysBoundary::SysBoundary() { 
   ++nSysBoundaries;
}

/*! Destructor for class SysBoundary. Reduces the value of SysBoundary::nSysBoundaries by one,
 * and if after the destruction SysBoundary::nSysBoundaries equals zero all stored SysBoundaries are deleted.
 */
SysBoundary::~SysBoundary() {
   --nSysBoundaries;
   if (nSysBoundaries != 0) return;
   
   // Call delete for each SysBoundaryCondition:
   for (vector<SBC::SysBoundaryCondition*>::iterator it=sysBoundaries.begin();
        it!=sysBoundaries.end();
        ++it) {
      delete *it;
      *it = NULL;
   }
}

/*! Add a new SBC::SysBoundaryCondition which has been created with new sysBoundary. 
 * SysBoundary will take care of deleting it.
 * @return If true, the given SBC::SysBoundaryCondition was added successfully.
 */
bool SysBoundary::addSysBoundary(SBC::SysBoundaryCondition* bc, creal& t) {
   sysBoundaries.push_back(bc);
   
   bool success = true;
   if(bc->initSysBoundary(t) == false) {
      cerr << "Error in setting system boundary condition " << bc->getName() << endl;
      success = false;
   }
   
   indexToSysBoundary[bc->getIndex()] = bc;
   
   indexToPrecedence[bc->getIndex()] = bc->getPrecedence();
   
   return success;
}

/*! Get the name of a SysBoundaryCondition.
 * @param sysBoundaryID ID number of the system boundary whose name is requested.
 * @return Name of the system boundary.
 */
string SysBoundary::getName(const unsigned int& sysBoundaryID) const {
   if (sysBoundaryID >= sysBoundaries.size()) return "";
   return sysBoundaries[sysBoundaryID]->getName();
}

/*! Get the list of system boundary conditions contained in SysBoundary.
 * @return List of the system boundary conditions.
 */
vector<SBC::SysBoundaryCondition*> SysBoundary::getSysBoundariesList() const {return sysBoundaries;}

/*! Get the map of system boundary condition precedences contained in SysBoundary.
 * @return Map of the system boundary condition precedences.
 */
std::map<uint, uint> SysBoundary::getPrecedenceMap() const {return indexToPrecedence;}

/*! Get a pointer to the SysBoundaryCondition of given index.
 * @return Pointer to the instance of the SysBoundaryCondition.
 */
SBC::SysBoundaryCondition* SysBoundary::getSysBoundary(uint sysBoundaryType) const {
   return indexToSysBoundary.find(sysBoundaryType)->second;
}

/*! Get the number of SysBoundaryConditions stored in SysBoundary.
 * @return Number of SysBoundaryConditions stored in SysBoundary.
 */
unsigned int SysBoundary::size() const {return sysBoundaries.size();}

