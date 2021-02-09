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

/*!\file nothing.cpp
 * \brief Implementation of the class BoundaryCondition::Nothing to handle
 * cells classified as boundarytype::NOTHING.
 */

#include <cstdlib>
#include <iostream>

#include "../object_wrapper.h"
#include "nothing.h"

using namespace std;

namespace BC {
Nothing::Nothing() : BoundaryCondition() {}
Nothing::~Nothing() {}

void Nothing::addParameters() {}
void Nothing::getParameters() {}

void Nothing::initBoundary(creal t, Project &project) {
   precedence = 0;
   dynamic = false;
}

void Nothing::assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &,
                               FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid) {}

void Nothing::applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                  FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                                  Project &project) {
   vector<CellID> cells = mpiGrid.get_cells();
#pragma omp parallel for
   for (size_t i = 0; i < cells.size(); ++i) {
      SpatialCell *cell = mpiGrid[cells[i]];
      if (cell->boundaryFlag != this->getIndex())
         continue;

      // TODO: Set fields on B grid to 0
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

      // let's get rid of blocks not fulfilling the criteria here to save
      // memory.
      for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID)
         cell->adjustSingleCellVelocityBlocks(popID);
   }
}

void Nothing::updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                            FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid, creal t) {}

std::string Nothing::getName() const { return "Nothing"; }

uint Nothing::getIndex() const { return boundarytype::NOTHING; }
} // namespace BC
