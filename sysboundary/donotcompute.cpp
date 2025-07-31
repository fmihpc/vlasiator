/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
   
   void DoNotCompute::initSysBoundary(
      creal& t,
      Project &project
   ) {
      precedence = 0;
      dynamic = false;
   }
   
   void DoNotCompute::assignSysBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>&,
                                        std::span< fsgrids::technical> technical, FieldSolverGrid &fsgrid) {
      // Does nothing.
   }

   void DoNotCompute::applyInitialState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                        std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid,
                                        std::span<array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                        fsgrids::bgbspan bgb, Project&) {
      const vector<CellID>& cells = getLocalCells();
#pragma omp parallel for
      for (size_t i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if(cell->sysBoundaryFlag != this->getIndex()) continue;

         //TODO: Set fields on B grid to 0         
         cell->parameters[CellParams::RHOM] = 0.0;
         cell->parameters[CellParams::VX] = 0.0;
         cell->parameters[CellParams::VY] = 0.0;
         cell->parameters[CellParams::VZ] = 0.0;
         cell->parameters[CellParams::RHOQ] = 0.0;
         cell->parameters[CellParams::RHOM_DT2] = 0.0;
         cell->parameters[CellParams::VX_DT2] = 0.0;
         cell->parameters[CellParams::VY_DT2] = 0.0;
         cell->parameters[CellParams::VZ_DT2] = 0.0;
         cell->parameters[CellParams::RHOQ_DT2] = 0.0;
      }
   }

   void DoNotCompute::updateState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                  std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid,
                                  std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                  fsgrids::bgbspan bgb, creal t) {}

   void DoNotCompute::getFaces(bool *faces) {}

   string DoNotCompute::getName() const {return "DoNotCompute";}
   
   uint DoNotCompute::getIndex() const {return sysboundarytype::DO_NOT_COMPUTE;}
}
