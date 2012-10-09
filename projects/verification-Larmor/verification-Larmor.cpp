/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"
#include "vlasovmover.h"

typedef larmorParameters LarmP;
Real LarmP::BX0 = NAN;
Real LarmP::BY0 = NAN;
Real LarmP::BZ0 = NAN;
Real LarmP::VX0 = NAN;
Real LarmP::VY0 = NAN;
Real LarmP::VZ0 = NAN;
Real LarmP::X0 = NAN;
Real LarmP::Y0 = NAN;
Real LarmP::Z0 = NAN;
Real LarmP::DENSITY = NAN;


bool initializeProject(void) {return true;}

bool addProjectParameters() {
   typedef Readparameters RP;
   RP::add("Larmor.BX0", "Background field value (T)", 0.0);
   RP::add("Larmor.BY0", "Background field value (T)", 0.0);
   RP::add("Larmor.BZ0", "Background field value (T)", 0.0);
   RP::add("Larmor.VX0", "Bulk velocity in x", 0.0);
   RP::add("Larmor.VY0", "Bulk velocity in y", 0.0);
   RP::add("Larmor.VZ0", "Bulk velocity in z", 0.0);
   RP::add("Larmor.X0", "Initial Position", 0.0);
   RP::add("Larmor.Y0", "Initial Position", 0.0);
   RP::add("Larmor.Z0", "Initial Position", 0.0);
   RP::add("Larmor.rho", "Number density (m^-3)", 1.0e7);
   return true;
}

bool getProjectParameters() {
   typedef Readparameters RP;
   RP::get("Larmor.BX0", LarmP::BX0);
   RP::get("Larmor.BY0", LarmP::BY0);
   RP::get("Larmor.BZ0", LarmP::BZ0);
   RP::get("Larmor.VX0", LarmP::VX0);
   RP::get("Larmor.VY0", LarmP::VY0);
   RP::get("Larmor.VZ0", LarmP::VZ0);
   RP::get("Larmor.X0", LarmP::X0);
   RP::get("Larmor.Y0", LarmP::Y0);
   RP::get("Larmor.Z0", LarmP::Z0);
   RP::get("Larmor.rho", LarmP::DENSITY);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

void setProjectCell(SpatialCell* cell) {
   // Set up cell parameters:
   calcCellParameters(&((*cell).parameters[0]), 0.0);
   
   cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
   cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
   
   // Go through each velocity block in the velocity phase space grid.
   // Set the initial state and block parameters:
   creal dvx_block = SpatialCell::block_dvx; // Size of a block in vx-direction
   creal dvy_block = SpatialCell::block_dvy; //                    vy
   creal dvz_block = SpatialCell::block_dvz; //                    vz
   creal dvx_blockCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
   creal dvy_blockCell = SpatialCell::cell_dvy; //                                vy
   creal dvz_blockCell = SpatialCell::cell_dvz; //                                vz
   
   for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
      for (uint jv=0; jv<P::vyblocks_ini; ++jv)
         for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
            creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
            creal vy_block = P::vymin + jv*dvy_block; // vy-
            creal vz_block = P::vzmin + kv*dvz_block; // vz-
            
            // Calculate volume average of distrib. function for each cell in the block.
            for (uint kc=0; kc<WID; ++kc) 
               for (uint jc=0; jc<WID; ++jc) 
                  for (uint ic=0; ic<WID; ++ic) {
                     creal vx_cell = vx_block + ic*dvx_blockCell;
                     creal vy_cell = vy_block + jc*dvy_blockCell;
                     creal vz_cell = vz_block + kc*dvz_blockCell;
                     Real average = 
                     calcPhaseSpaceDensity(cell->parameters[CellParams::XCRD],
                                           cell->parameters[CellParams::YCRD],
                                           cell->parameters[CellParams::ZCRD],
                                           cell->parameters[CellParams::DX],
                                           cell->parameters[CellParams::DY],
                                           cell->parameters[CellParams::DZ],
                                           vx_cell,vy_cell,vz_cell,
                                           dvx_blockCell,dvy_blockCell,dvz_blockCell);
                     
                     if(average!=0.0){
                        creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
                        creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
                        creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
                        cell->set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
                     }
                  }
         }
         calculateCellVelocityMoments(cell);
         
         //let's get rid of blocks not fulfilling the criteria here to save memory.
         cell->adjustSingleCellVelocityBlocks();
}


Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {

   static bool isSet=false;
   //static variables should be threadprivate
#pragma omp threadprivate(isSet)

   if(vx < Parameters::vxmin + 0.5 * dvx ||
      vy < Parameters::vymin + 0.5 * dvy ||
      vz < Parameters::vzmin + 0.5 * dvz ||
      vx > Parameters::vxmax - 1.5 * dvx ||
      vy > Parameters::vymax - 1.5 * dvy ||
      vz > Parameters::vzmax - 1.5 * dvz
   ) return 0.0;

   if(isSet)
      return 0.0; //exactly one value to be set

   
   creal mass = Parameters::m;
   creal q = Parameters::q;
   creal k = 1.3806505e-23; // Boltzmann
   creal mu0 = 1.25663706144e-6; // mu_0

   if( fabs(vx-LarmP::VX0)<dvx &&
       fabs(vy-LarmP::VY0)<dvy &&
       fabs(vz-LarmP::VZ0)<dvz &&
       fabs(x-LarmP::X0)<dx &&
       fabs(y-LarmP::Y0)<dy &&
       fabs(z-LarmP::Z0)<dz){
      isSet=true;
      return LarmP::DENSITY/(dvx*dvy*dvz);
   }

   return 0.0;
   
}
      
void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD];
   creal dy = cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD];
   creal dz = cellParams[CellParams::DZ];
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::PERBX   ] = 0.0;
   cellParams[CellParams::PERBY   ] = 0.0;
   cellParams[CellParams::PERBZ   ] = 0.0;
   cellParams[CellParams::BGBX   ] = LarmP::BX0;
   cellParams[CellParams::BGBY   ] = LarmP::BY0;
   cellParams[CellParams::BGBZ   ] = LarmP::BZ0;

}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   const std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

