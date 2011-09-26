/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

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

#ifndef WRITEVARS_H
#define WRITEVARS_H

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "definitions.h"
#include "cell_spatial.h"

#ifdef PARGRID
   #include "pargrid.h"
#endif

// Variable saving needs to be completely rewritten

void openOutputFile(const std::string& fileName);
void closeOutputFile(const std::string& fileName);

#ifdef PARGRID
   template<class C> void writeCellDistribution(const C& mpiGrid);
#endif
void writeVelocityBlockGrid(const std::string& gridName,cuint& BLOCKS,Real* blockParams);
void writeVelocityBlockScalar(const std::string& varName,const std::string& gridName,cuint& BLOCKS,Real* array);


#ifdef PARGRID
template<> void writeCellDistribution(const ParGrid<SpatialCell>& mpiGrid) {
   if (mpiGrid.rank() != 0) return;
   
   // Get all cells and their hosts:
   std::map<ID::type,int> hosts;
   mpiGrid.getCellDistribution(hosts);
   
   // Open output file:
   std::stringstream fname;
   fname << "loadbalance." << mpiGrid.rank() << '.';
   fname.width(5);
   fname.fill('0');
   fname << Parameters::tstep << ".silo";
   openOutputFile(fname.str(),"spatial_cells");
   reserveSpatialCells(hosts.size());
   
   // Write each cell and its host:
   SpatialCell cell;
   for (std::map<ID::type,int>::const_iterator it=hosts.begin(); it!=hosts.end(); ++it) {
      cell.cpu_cellParams[CellParams::DX] = mpiGrid.get_cell_x_size(it->first);
      cell.cpu_cellParams[CellParams::DY] = mpiGrid.get_cell_y_size(it->first);
      cell.cpu_cellParams[CellParams::DZ] = mpiGrid.get_cell_z_size(it->first);
      cell.cpu_cellParams[CellParams::XCRD] = mpiGrid.get_cell_x(it->first) - 0.5*cell.cpu_cellParams[CellParams::DX];
      cell.cpu_cellParams[CellParams::YCRD] = mpiGrid.get_cell_y(it->first) - 0.5*cell.cpu_cellParams[CellParams::DY];
      cell.cpu_cellParams[CellParams::ZCRD] = mpiGrid.get_cell_z(it->first) - 0.5*cell.cpu_cellParams[CellParams::DZ];
      addSpatialCell(cell.cpu_cellParams);
   }
   
   // Write output file and release resources:
   writeSpatialCells("spatcells");
   closeOutputFile();
   freeCells();
}
#endif

#endif
