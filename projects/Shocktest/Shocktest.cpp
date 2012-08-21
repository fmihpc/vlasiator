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

typedef riemannParameters RiP;
Real RiP::rho[2] = {NAN};
Real RiP::T[2] = {NAN};
Real RiP::Vx[2] = {NAN};
Real RiP::Vy[2] = {NAN};
Real RiP::Vz[2] = {NAN};
Real RiP::Bx[2] = {NAN};
Real RiP::By[2] = {NAN};
Real RiP::Bz[2] = {NAN};
uint RiP::nSpaceSamples = 0;
uint RiP::nVelocitySamples = 0;

bool initializeProject(void) {return true;}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("Riemann.rho1", "Number density, left state (m^-3)", 0.0);
   RP::add("Riemann.rho2", "Number density, right state (m^-3)", 0.0);
   RP::add("Riemann.T1", "Temperature, left state (K)", 0.0);
   RP::add("Riemann.T2", "Temperature, right state (K)", 0.0);
   RP::add("Riemann.Vx1", "Bulk velocity x component, left state (m/s)", 0.0);
   RP::add("Riemann.Vx2", "Bulk velocity x component, right state (m/s)", 0.0);
   RP::add("Riemann.Vy1", "Bulk velocity y component, left state (m/s)", 0.0);
   RP::add("Riemann.Vy2", "Bulk velocity y component, right state (m/s)", 0.0);
   RP::add("Riemann.Vz1", "Bulk velocity z component, left state (m/s)", 0.0);
   RP::add("Riemann.Vz2", "Bulk velocity z component, right state (m/s)", 0.0);
   RP::add("Riemann.Bx1", "Magnetic field x component, left state (T)", 0.0);
   RP::add("Riemann.Bx2", "Magnetic field x component, right state (T)", 0.0);
   RP::add("Riemann.By1", "Magnetic field y component, left state (T)", 0.0);
   RP::add("Riemann.By2", "Magnetic field y component, right state (T)", 0.0);
   RP::add("Riemann.Bz1", "Magnetic field z component, left state (T)", 0.0);
   RP::add("Riemann.Bz2", "Magnetic field z component, right state (T)", 0.0);
   RP::add("Riemann.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Riemann.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("Riemann.rho1", RiP::rho[RiP::LEFT]);
   RP::get("Riemann.rho2", RiP::rho[RiP::RIGHT]);
   RP::get("Riemann.T1", RiP::T[RiP::LEFT]);
   RP::get("Riemann.T2", RiP::T[RiP::RIGHT]);
   RP::get("Riemann.Vx1", RiP::Vx[RiP::LEFT]);
   RP::get("Riemann.Vx2", RiP::Vx[RiP::RIGHT]);
   RP::get("Riemann.Vy1", RiP::Vy[RiP::LEFT]);
   RP::get("Riemann.Vy2", RiP::Vy[RiP::RIGHT]);
   RP::get("Riemann.Vz1", RiP::Vz[RiP::LEFT]);
   RP::get("Riemann.Vz2", RiP::Vz[RiP::RIGHT]);
   RP::get("Riemann.Bx1", RiP::Bx[RiP::LEFT]);
   RP::get("Riemann.Bx2", RiP::Bx[RiP::RIGHT]);
   RP::get("Riemann.By1", RiP::By[RiP::LEFT]);
   RP::get("Riemann.By2", RiP::By[RiP::RIGHT]);
   RP::get("Riemann.Bz1", RiP::Bz[RiP::LEFT]);
   RP::get("Riemann.Bz2", RiP::Bz[RiP::RIGHT]);
   RP::get("Riemann.nSpaceSamples", RiP::nSpaceSamples);
   RP::get("Riemann.nVelocitySamples", RiP::nVelocitySamples);
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

Real getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal mass = 1.67262171e-27; // m_p in kg
   creal k = 1.3806505e-23; // Boltzmann
   //  creal mu0 = 1.25663706144e-6; // mu_0
   //  creal q = 1.60217653e-19; // q_i
   //  creal gamma = 5./3.;
   
   cint side = (x < 0.0) ? RiP::LEFT : RiP::RIGHT;
   
   return RiP::rho[side] * pow(mass / (2.0 * M_PI * k * RiP::T[side]), 1.5) *
   exp(- mass * (pow(vx - RiP::Vx[side], 2.0) + pow(vy - RiP::Vy[side], 2.0) + pow(vz - RiP::Vz[side], 2.0)) / (2.0 * k * RiP::T[side]));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
   creal d_x = dx / (RiP::nSpaceSamples-1);
   creal d_y = dy / (RiP::nSpaceSamples-1);
   creal d_z = dz / (RiP::nSpaceSamples-1);
   creal d_vx = dvx / (RiP::nVelocitySamples-1);
   creal d_vy = dvy / (RiP::nVelocitySamples-1);
   creal d_vz = dvz / (RiP::nVelocitySamples-1);
   Real avg = 0.0;
   for (uint i=0; i<RiP::nSpaceSamples; ++i)
      for (uint j=0; j<RiP::nSpaceSamples; ++j)
	 for (uint k=0; k<RiP::nSpaceSamples; ++k)
	    for (uint vi=0; vi<RiP::nVelocitySamples; ++vi)
	       for (uint vj=0; vj<RiP::nVelocitySamples; ++vj)
		  for (uint vk=0; vk<RiP::nVelocitySamples; ++vk)
		     {
			avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
		     }
   return avg / pow(RiP::nSpaceSamples, 3.0) / pow(RiP::nVelocitySamples, 3.0);
}

bool cellParametersChanged(creal& t) {return false;}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   
   Real Bxavg, Byavg, Bzavg;
   Bxavg = Byavg = Bzavg = 0.0;
   Real d_x = dx / (RiP::nSpaceSamples - 1);
   
   for (uint i=0; i<RiP::nSpaceSamples; ++i)
      for (uint j=0; j<RiP::nSpaceSamples; ++j)
	 for (uint k=0; k<RiP::nSpaceSamples; ++k) {
	    Bxavg += ((x + i * d_x) < 0.0) ? RiP::Bx[RiP::LEFT] : RiP::Bx[RiP::RIGHT];
	    Byavg += ((x + i * d_x) < 0.0) ? RiP::By[RiP::LEFT] : RiP::By[RiP::RIGHT];
	    Bzavg += ((x + i * d_x) < 0.0) ? RiP::Bz[RiP::LEFT] : RiP::Bz[RiP::RIGHT];
	 }
   cuint nPts = pow(RiP::nSpaceSamples, 3.0);
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = Bxavg / nPts;
   cellParams[CellParams::BY   ] = Byavg / nPts;
   cellParams[CellParams::BZ   ] = Bzavg / nPts;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

