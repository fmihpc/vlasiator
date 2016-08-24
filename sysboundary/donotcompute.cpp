/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */

/*!\file donotcompute.cpp
 * \brief Implementation of the class SysBoundaryCondition::DoNotCompute to handle cells classified as sysboundarytype::DO_NOT_COMPUTE.
 */

#include <cstdlib>
#include <iostream>

#include "donotcompute.h"
#include "../object_wrapper.h"

using namespace std;

namespace SBC {
   DoNotCompute::DoNotCompute(): SysBoundaryCondition() { }
   DoNotCompute::~DoNotCompute() { }
   
   void DoNotCompute::addParameters() { }
   void DoNotCompute::getParameters() { }
   
   bool DoNotCompute::initSysBoundary(
      creal& t,
      Project &project
   ) {
      precedence = 0;
      isThisDynamic = false;
      return true;
   }
   
   bool DoNotCompute::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& ) {
      return true;
   }
   
   bool DoNotCompute::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      Project&
   ) {
      vector<CellID> cells = mpiGrid.get_cells();
#pragma omp parallel for
      for (size_t i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if(cell->sysBoundaryFlag != this->getIndex()) continue;
         
         cell->parameters[CellParams::PERBX] = 0.0;
         cell->parameters[CellParams::PERBY] = 0.0;
         cell->parameters[CellParams::PERBZ] = 0.0;
         cell->parameters[CellParams::PERBX_DT2] = 0.0;
         cell->parameters[CellParams::PERBY_DT2] = 0.0;
         cell->parameters[CellParams::PERBZ_DT2] = 0.0;
         cell->parameters[CellParams::EX] = 0.0;
         cell->parameters[CellParams::EY] = 0.0;
         cell->parameters[CellParams::EZ] = 0.0;
         cell->parameters[CellParams::RHO] = 0.0;
         cell->parameters[CellParams::RHOVX] = 0.0;
         cell->parameters[CellParams::RHOVY] = 0.0;
         cell->parameters[CellParams::RHOVZ] = 0.0;
         cell->parameters[CellParams::RHO_DT2] = 0.0;
         cell->parameters[CellParams::RHOVX_DT2] = 0.0;
         cell->parameters[CellParams::RHOVY_DT2] = 0.0;
         cell->parameters[CellParams::RHOVZ_DT2] = 0.0;
         cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
         cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
         
         //let's get rid of blocks not fulfilling the criteria here to save
         //memory.
         for (unsigned int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
            cell->adjustSingleCellVelocityBlocks(popID);
      }
      
      return true;
   }
   
   std::string DoNotCompute::getName() const {return "DoNotCompute";}
   
   uint DoNotCompute::getIndex() const {return sysboundarytype::DO_NOT_COMPUTE;}
}
