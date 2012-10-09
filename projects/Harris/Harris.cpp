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

using namespace std;

typedef harrisParameters HP;
Real HP::SCA_LAMBDA = NAN;
Real HP::B0 = NAN;
Real HP::TEMPERATURE = NAN;
Real HP::DENSITY = NAN;

bool initializeProject(void) {
   return true;
}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("Harris.Scale_size", "Harris sheet scale size (m)", 150000.0);
   RP::add("Harris.B0", "Magnetic field at infinity (T)", 8.33061003094e-8);
   RP::add("Harris.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Harris.rho", "Number density at infinity (m^-3)", 1.0e7);
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("Harris.Scale_size", HP::SCA_LAMBDA);
   RP::get("Harris.B0", HP::B0);
   RP::get("Harris.Temperature", HP::TEMPERATURE);
   RP::get("Harris.rho", HP::DENSITY);
   return true;
}

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

bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
  
  creal mass = 1.67262171e-27; // m_p in kg
  creal k = 1.3806505e-23; // Boltzmann
  creal mu0 = 1.25663706144e-6; // mu_0
  creal q = 1.60217653e-19; // q_i
  
  creal Vy0 = 0.0;
  
  return HP::DENSITY * pow(mass / (2.0 * M_PI * k * HP::TEMPERATURE), 1.5) * (
     5.0 / pow(cosh((x + 0.5 * dx) / (HP::SCA_LAMBDA)), 2.0) * exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy - Vy0, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * k * HP::TEMPERATURE))
    +
    exp(- mass * (pow(vx + 0.5 * dvx, 2.0) + pow(vy + 0.5 * dvy, 2.0) + pow(vz + 0.5 * dvz, 2.0)) / (2.0 * k * HP::TEMPERATURE)));
}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
  creal x = cellParams[CellParams::XCRD];
  creal dx = cellParams[CellParams::DX];
  
  cellParams[CellParams::EX   ] = 0.0;
  cellParams[CellParams::EY   ] = 0.0;
  cellParams[CellParams::EZ   ] = 0.0;
  cellParams[CellParams::PERBX   ] = 0.0;
  cellParams[CellParams::PERBY   ] = 0.0;
  cellParams[CellParams::PERBZ   ] = 0.0;
  cellParams[CellParams::BGBX   ] = 0.0;
  cellParams[CellParams::BGBY   ] = 0.0;
  cellParams[CellParams::BGBZ   ] = HP::B0 * tanh((x + 0.5 * dx) / HP::SCA_LAMBDA);
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

