/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012 Finnish Meteorological Institute
 * 
 * Vlasiator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 * 
 * Vlasiator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/*!\file donotcompute.cpp
 * \brief Implementation of the class SysBoundaryCondition::DoNotCompute to handle cells classified as sysboundarytype::DO_NOT_COMPUTE.
 */

#include <cstdlib>
#include <iostream>

#include "donotcompute.h"

using namespace std;

namespace SBC {
   DoNotCompute::DoNotCompute(): SysBoundaryCondition() { }
   DoNotCompute::~DoNotCompute() { }
   
   void DoNotCompute::addParameters() { }
   void DoNotCompute::getParameters() { }
   
   bool DoNotCompute::initSysBoundary(creal& t) {
      precedence = 0;
      isThisDynamic = false;
      return true;
   }
   
   int DoNotCompute::assignSysBoundary(creal* cellParams) {
      return sysboundarytype::NOT_SYSBOUNDARY;
   }
   
   bool DoNotCompute::applyInitialState(dccrg::Dccrg<SpatialCell>& mpiGrid) {
      vector<uint64_t> cells = mpiGrid.get_cells();
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if(cell->sysBoundaryFlag != this->getIndex()) continue;
         
         cell->parameters[CellParams::BX] = 0.0;
         cell->parameters[CellParams::BY] = 0.0;
         cell->parameters[CellParams::BZ] = 0.0;
         cell->parameters[CellParams::EX] = 0.0;
         cell->parameters[CellParams::EY] = 0.0;
         cell->parameters[CellParams::EZ] = 0.0;
         cell->parameters[CellParams::RHO] = 0.0;
         cell->parameters[CellParams::RHOVX] =0.0;
         cell->parameters[CellParams::RHOVY] = 0.0;
         cell->parameters[CellParams::RHOVZ] =0.0;
         cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
         cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
         
         //let's get rid of blocks not fulfilling the criteria here to save
         //memory.
         cell->adjustSingleCellVelocityBlocks();
      }
      
      return true;
   }
   
   std::string DoNotCompute::getName() const {return "DoNotCompute";}
   
   uint DoNotCompute::getIndex() const {return sysboundarytype::DO_NOT_COMPUTE;}
}
