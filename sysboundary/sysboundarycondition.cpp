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
#include <mpi.h>
#include <iostream>
#include <limits>

#include "../parameters.h"
#include "sysboundarycondition.h"
#include "ionosphere.h"
#include "outflow.h"
#include "solarwind.h"

using namespace std;

namespace SBC {
   // ************************************************************
   // ***** DEFINITIONS FOR BOUNDARYCONDITION BASE CLASS *****
   // ************************************************************
   
   /*! SysBoundaryCondition base class constructor. The constructor is empty.*/
   SysBoundaryCondition::SysBoundaryCondition() { }
   
   /*! SysBoundaryCondition base class virtual destructor. The destructor is empty.*/
   SysBoundaryCondition::~SysBoundaryCondition() { }
   
   void SysBoundaryCondition::addParameters() {
      cerr << "ERROR: SysBoundaryCondition::addParameters called instead of derived class function!" << endl;
   }
   
   void SysBoundaryCondition::getParameters() {
      cerr << "ERROR: SysBoundaryCondition::getParameters called instead of derived class function!" << endl;
   }
   
   /*! Function called at initialisation to set the system boundary condition's parameters.
    */
   bool SysBoundaryCondition::initSysBoundary(creal& t) {
      cerr << "ERROR: SysBoundaryCondition::initSysBoundary called instead of derived class function!" << endl;
      return false;
   }
   
   /*! Function used to assign the system boundary condition type to a cell.
    * \param 
    * \return The system boundary condition type's index
    */
   int SysBoundaryCondition::assignSysBoundary(creal* cellParams) {
      cerr << "ERROR: SysBoundaryCondition::assignSysBoundary called instead of derived class function!" << endl;
      return false;
   }
   
   /*! Function used to apply the system boundary condition initial state to a cell. */
   bool SysBoundaryCondition::applyInitialState(dccrg::Dccrg<SpatialCell>& mpiGrid) {
      cerr << "ERROR: SysBoundaryCondition::applyInitialState called instead of derived class function!" << endl;
      return false;
   }
      
   /*! Get the name of the system boundary condition.
    * @return The name of the data. The base class function returns an empty string.
    */
   std::string SysBoundaryCondition::getName() const {
      cerr << "ERROR: SysBoundaryCondition::getName called instead of derived class function!" << endl;
      return string("");
   }
   
   /*! Get the enum index of the system boundary condition.
    * @return The index of the system boundary condition as enumerated in namespace sysboundarytype. The base class function returns 0.
    */
   uint SysBoundaryCondition::getIndex() const {
      cerr << "ERROR: SysBoundaryCondition::getIndex called instead of derived class function!" << endl;
      return 0;
   }
   
   /*! Get the precedence value of the system boundary condition.
    * @return The precedence value of the system boundary condition as set by parameter. The base class function returns 0.
    */
   uint SysBoundaryCondition::getPrecedence() const {
      cerr << "ERROR: SysBoundaryCondition::getPrecedence called instead of derived class function!" << endl;
      return 0;
   }
   
   /*! Returns whether this boundary condition is dynamic in time.
    * @return Boolean value.
    */
   bool SysBoundaryCondition::isDynamic() const {
      cerr << "ERROR: SysBoundaryCondition::isDynamic called instead of derived class function!" << endl;
      return false;
   }
   
} // namespace SBC
