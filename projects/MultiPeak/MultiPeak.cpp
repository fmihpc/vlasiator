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
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"
#include "vlasovmover.h"

using namespace std;

typedef multiPeakParameters MPP;
Real MPP::rho[2] = {NAN};
Real MPP::T[2] = {NAN};
Real MPP::Vx[2] = {NAN};
Real MPP::Vy[2] = {NAN};
Real MPP::Vz[2] = {NAN};
Real MPP::Bx = NAN;
Real MPP::By = NAN;
Real MPP::Bz = NAN;
uint MPP::nVelocitySamples = 0;

bool initializeProject(void) {return true;}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("MultiPeak.rho1", "Number density, first peak (m^-3)", 0.0);
   RP::add("MultiPeak.rho2", "Number density, second peak (m^-3)", 0.0);
   RP::add("MultiPeak.T1", "Temperature, first peak (K)", 0.0);
   RP::add("MultiPeak.T2", "Temperature, second peak (K)", 0.0);
   RP::add("MultiPeak.Vx1", "Bulk velocity x component, first peak (m/s)", 0.0);
   RP::add("MultiPeak.Vx2", "Bulk velocity x component, second peak (m/s)", 0.0);
   RP::add("MultiPeak.Vy1", "Bulk velocity y component, first peak (m/s)", 0.0);
   RP::add("MultiPeak.Vy2", "Bulk velocity y component, second peak (m/s)", 0.0);
   RP::add("MultiPeak.Vz1", "Bulk velocity z component, first peak (m/s)", 0.0);
   RP::add("MultiPeak.Vz2", "Bulk velocity z component, second peak (m/s)", 0.0);
   RP::add("MultiPeak.Bx", "Magnetic field x component (T)", 0.0);
   RP::add("MultiPeak.By", "Magnetic field y component (T)", 0.0);
   RP::add("MultiPeak.Bz", "Magnetic field z component (T)", 0.0);
   RP::add("MultiPeak.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("MultiPeak.rho1", MPP::rho[1]);
   RP::get("MultiPeak.rho2", MPP::rho[2]);
   RP::get("MultiPeak.T1", MPP::T[1]);
   RP::get("MultiPeak.T2", MPP::T[2]);
   RP::get("MultiPeak.Vx1", MPP::Vx[1]);
   RP::get("MultiPeak.Vx2", MPP::Vx[2]);
   RP::get("MultiPeak.Vy1", MPP::Vy[1]);
   RP::get("MultiPeak.Vy2", MPP::Vy[2]);
   RP::get("MultiPeak.Vz1", MPP::Vz[1]);
   RP::get("MultiPeak.Vz2", MPP::Vz[2]);
   RP::get("MultiPeak.Bx", MPP::Bx);
   RP::get("MultiPeak.By", MPP::By);
   RP::get("MultiPeak.Bz", MPP::Bz);
   RP::get("MultiPeak.nVelocitySamples", MPP::nVelocitySamples);
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

Real getDistribValue(creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal mass = 1.67262171e-27; // m_p in kg
   creal k = 1.3806505e-23; // Boltzmann
   //  creal mu0 = 1.25663706144e-6; // mu_0
   //  creal q = 1.60217653e-19; // q_i
   //  creal gamma = 5./3.;
   
   return
   MPP::rho[1] * pow(mass / (2.0 * M_PI * k * MPP::T[1]), 1.5) *
   exp(- mass * (pow(vx - MPP::Vx[1], 2.0) + pow(vy - MPP::Vy[1], 2.0) + pow(vz - MPP::Vz[1], 2.0)) / (2.0 * k * MPP::T[1])) +
   MPP::rho[2] * pow(mass / (2.0 * M_PI * k * MPP::T[2]), 1.5) *
   exp(- mass * (pow(vx - MPP::Vx[2], 2.0) + pow(vy - MPP::Vy[2], 2.0) + pow(vz - MPP::Vz[2], 2.0)) / (2.0 * k * MPP::T[2]));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
   creal d_vx = dvx / (MPP::nVelocitySamples-1);
   creal d_vy = dvy / (MPP::nVelocitySamples-1);
   creal d_vz = dvz / (MPP::nVelocitySamples-1);
   Real avg = 0.0;
//#pragma omp parallel for collapse(6) reduction(+:avg)
   for (uint vi=0; vi<MPP::nVelocitySamples; ++vi)
      for (uint vj=0; vj<MPP::nVelocitySamples; ++vj)
	 for (uint vk=0; vk<MPP::nVelocitySamples; ++vk)
	    {
	       avg += getDistribValue(vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
	    }
   return avg / pow(MPP::nVelocitySamples, 3.0);
}

bool cellParametersChanged(creal& t) {return false;}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = MPP::Bx;
   cellParams[CellParams::BY   ] = MPP::By;
   cellParams[CellParams::BZ   ] = MPP::Bz;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

