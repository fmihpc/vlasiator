/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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
#ifndef IOWRITE_H
#define IOWRITE_H
#include "mpi.h"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <string>
#include <vector>
#include <vlsv_writer.h>

#include "spatial_cell.hpp"
#include "datareduction/datareducer.h"

/*!

\brief Write out system into a vlsv file

\param mpiGrid     The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute data that is added into file
\param index       Index to call the correct member of the various parameter vectors
\param writeGhosts Write ghost zones
*/
bool writeGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
               DataReducer* dataReducer,
               const uint& index,
               const bool writeGhosts = true
);

/*!

\brief Write out a restart of the simulation into a vlsv file. All block data in remote cells will be reset.

\param mpiGrid   The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute data that is added into file
\param name       File name prefix, file will be called "name.index.vlsv"
\param fileIndex  File index, file will be called "name.index.vlsv"
*/
bool writeRestart(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  DataReducer& dataReducer,
                  const std::string& name,
                  const uint& fileIndex,
                  const int& stripe);

/*!

\brief Write out simulation diagnostics into diagnostic.txt

@param mpiGrid   The DCCRG grid with spatial cells
@param dataReducer Contains datareductionoperators that are used to compute diagnostic data
*/
bool writeDiagnostic(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,DataReducer& dataReducer);

bool writeVelocitySpace(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                        vlsv::Writer& vlsvWriter,int index,const std::vector<uint64_t>& cells);

bool writeVelocityDistributionData(vlsv::Writer& vlsvWriter,dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                   const std::vector<uint64_t>& cells,MPI_Comm comm);

#endif
