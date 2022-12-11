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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef IOREAD_H
#define IOREAD_H
#include "mpi.h"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <string>

#include "vlsv_reader_parallel.h"
#include "definitions.h"
#include "spatial_cell.hpp"
#include "datareduction/datareducer.h"


/*!

\brief Read in state from a vlsv file in order to restart simulations
\param mpiGrid Vlasiator's grid
\param name Name of the restart file e.g. "restart.00052.vlsv"
*/
bool readGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
              const std::string& name);

/*!

\brief Refine the grid to be identical to the file's
\param mpiGrid Vlasiator's grid
\param name Name of the restart file e.g. "restart.00052.vlsv"
*/
bool readFileCells(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const std::string& name);

/*!
 * \brief Check in local directory for external commands passed to the simulation. Only executed by MASTER_RANK
 */
void checkExternalCommands();

bool readIonosphereNodeVariable(
   vlsv::ParallelReader& file, const string& variableName, SBC::SphericalTriGrid& grid, ionosphereParameters index);

#endif
