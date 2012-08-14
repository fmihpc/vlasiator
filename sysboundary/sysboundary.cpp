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

/*
//Instead of calcSysBoundaryPhaseSpaceDensity below, have higher-level function
void setState(mpiGrid){
   loop over bcond
      boundarycond.setstate(mpigrid);
}
*/

Real calcSysBoundaryPhaseSpaceDensity(SysBoundary& sysBoundaries,
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
         average = sysBoundaries.getSysBoundary(OUTFLOW)->calcPhaseSpaceDensity(x,y,z,
                                                                                 dx,dy,dz,
                                                                                 vx,vy,vz,
                                                                                 dvx,dvy,dvz);
         break;
      case IONOSPHERE:
         average = sysBoundaries.getSysBoundary(IONOSPHERE)->calcPhaseSpaceDensity(x,y,z,
                                                                                    dx,dy,dz,
                                                                                    vx,vy,vz,
                                                                                    dvx,dvy,dvz);
         break;
      case SW:
         average = sysBoundaries.getSysBoundary(SW)->calcPhaseSpaceDensity(x,y,z,
                                                                            dx,dy,dz,
                                                                            vx,vy,vz,
                                                                            dvx,dvy,dvz);
         break;
      default:
         average = 0.0;
   }
   
   return average;
}

void calcSysBoundaryCellParameters(SysBoundary& sysBoundaries,
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
         sysBoundaries.getSysBoundary(OUTFLOW)->calcCellParameters(cellParams, t);
         break;
      case IONOSPHERE:
         sysBoundaries.getSysBoundary(IONOSPHERE)->calcCellParameters(cellParams, t);
         break;
      case SW:
         sysBoundaries.getSysBoundary(SW)->calcCellParameters(cellParams, t);
         break;
      default:
         break;
   }
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
   for (vector<SBC::SysBoundaryCondition*>::iterator it=sysBoundaries.begin();
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
   
   
   Readparameters::addComposing("outflow.face", "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
   
   Readparameters::addComposing("solarwind.face", "List of faces on which solar wind boundary conditions are to be applied ([xyz][+-]).");
   Readparameters::add("solarwind.file_x+", "Input files for the solar wind inflow parameters on face x+.", "");
   Readparameters::add("solarwind.file_x-", "Input files for the solar wind inflow parameters on face x-.", "");
   Readparameters::add("solarwind.file_y+", "Input files for the solar wind inflow parameters on face y+.", "");
   Readparameters::add("solarwind.file_y-", "Input files for the solar wind inflow parameters on face y-.", "");
   Readparameters::add("solarwind.file_z+", "Input files for the solar wind inflow parameters on face z+.", "");
   Readparameters::add("solarwind.file_z-", "Input files for the solar wind inflow parameters on face z-.", "");
   Readparameters::add("solarwind.dynamic", "Boolean value, is the solar wind inflow dynamic in time or not.", 0);
   
   Readparameters::add("ionosphere.centerX", "X coordinate of ionosphere center (m)", 0.0);
   Readparameters::add("ionosphere.centerY", "Y coordinate of ionosphere center (m)", 0.0);
   Readparameters::add("ionosphere.centerZ", "Z coordinate of ionosphere center (m)", 0.0);
   Readparameters::add("ionosphere.radius", "Radius of ionosphere (m).", 1.0e7);
   
   Readparameters::add("outflow.precedence", "Precedence value of the outflow system boundary condition (integer), the higher the stronger.", 3);
   Readparameters::add("solarwind.precedence", "Precedence value of the solar wind system boundary condition (integer), the higher the stronger.", 2);
   Readparameters::add("ionosphere.precedence", "Precedence value of the ionosphere system boundary condition (integer), the higher the stronger.", 1);
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

bool SysBoundary::initializeSysBoundaries(creal& t)
{
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

bool SysBoundary::assignSysBoundaryType(SpatialCell& cell)
{
   using namespace sysboundarytype;
   
   uint indexToAssign, tmpType;
   
   indexToAssign = NOT_SYSBOUNDARY;
   for (uint j = 0;
        j < sysBoundaries.size();
        j++) {
      tmpType=sysBoundaries[j]->assignSysBoundary(&(cell.parameters[0]));
      
      if(tmpType == DO_NOT_COMPUTE) {
         indexToAssign = tmpType;
         break;
      } else if (indexToPrecedence.find(tmpType)->second >
                 indexToPrecedence.find(indexToAssign)->second) {
         indexToAssign = tmpType;
      }
   }
   cell.sysBoundaryFlag = indexToAssign;
   
   return true;
}

/*! Get the name of a SysBoundaryCondition.
 * @param sysBoundaryID ID number of the system boundary whose name is requested.
 * @return Name of the system boundary.
 */
string SysBoundary::getName(const unsigned int& sysBoundaryID) const {
   if (sysBoundaryID >= sysBoundaries.size()) return "";
   return sysBoundaries[sysBoundaryID]->getName();
}

// /*! Get the list of system boundary conditions contained in SysBoundary.
//  * @return List of the system boundary conditions.
//  */
// vector<SBC::SysBoundaryCondition*> SysBoundary::getSysBoundariesList() const {return sysBoundaries;}

// /*! Get the map of system boundary condition precedences contained in SysBoundary.
//  * @return Map of the system boundary condition precedences.
//  */
// std::map<uint, uint> SysBoundary::getPrecedenceMap() const {return indexToPrecedence;}

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
bool SysBoundary::isDynamic() const {return isThisDynamic;}
bool SysBoundary::isBoundaryPeriodic(uint direction) const {return isPeriodic[direction];}
