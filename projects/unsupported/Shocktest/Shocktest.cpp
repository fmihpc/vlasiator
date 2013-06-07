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

typedef shocktestParameters stP;
Real stP::rho[2] = {NAN};
Real stP::T[2] = {NAN};
Real stP::Vx[2] = {NAN};
Real stP::Vy[2] = {NAN};
Real stP::Vz[2] = {NAN};
Real stP::Bx[2] = {NAN};
Real stP::By[2] = {NAN};
Real stP::Bz[2] = {NAN};
uint stP::nSpaceSamples = 0;
uint stP::nVelocitySamples = 0;

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
   RP::get("Riemann.rho1", stP::rho[stP::LEFT]);
   RP::get("Riemann.rho2", stP::rho[stP::RIGHT]);
   RP::get("Riemann.T1", stP::T[stP::LEFT]);
   RP::get("Riemann.T2", stP::T[stP::RIGHT]);
   RP::get("Riemann.Vx1", stP::Vx[stP::LEFT]);
   RP::get("Riemann.Vx2", stP::Vx[stP::RIGHT]);
   RP::get("Riemann.Vy1", stP::Vy[stP::LEFT]);
   RP::get("Riemann.Vy2", stP::Vy[stP::RIGHT]);
   RP::get("Riemann.Vz1", stP::Vz[stP::LEFT]);
   RP::get("Riemann.Vz2", stP::Vz[stP::RIGHT]);
   RP::get("Riemann.Bx1", stP::Bx[stP::LEFT]);
   RP::get("Riemann.Bx2", stP::Bx[stP::RIGHT]);
   RP::get("Riemann.By1", stP::By[stP::LEFT]);
   RP::get("Riemann.By2", stP::By[stP::RIGHT]);
   RP::get("Riemann.Bz1", stP::Bz[stP::LEFT]);
   RP::get("Riemann.Bz2", stP::Bz[stP::RIGHT]);
   RP::get("Riemann.nSpaceSamples", stP::nSpaceSamples);
   RP::get("Riemann.nVelocitySamples", stP::nVelocitySamples);
   return true;
}

Real getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal mass = 1.67262171e-27; // m_p in kg
   creal k = 1.3806505e-23; // Boltzmann
   //  creal mu0 = 1.25663706144e-6; // mu_0
   //  creal q = 1.60217653e-19; // q_i
   //  creal gamma = 5./3.;
   
   cint side = (x < 0.0) ? stP::LEFT : stP::RIGHT;
   
   return stP::rho[side] * pow(mass / (2.0 * M_PI * k * stP::T[side]), 1.5) *
   exp(- mass * (pow(vx - stP::Vx[side], 2.0) + pow(vy - stP::Vy[side], 2.0) + pow(vz - stP::Vz[side], 2.0)) / (2.0 * k * stP::T[side]));
}

/** Integrate the distribution function over the given six-dimensional phase-space cell.
 * @param x Starting value of the x-coordinate of the cell.
 * @param y Starting value of the y-coordinate of the cell.
 * @param z Starting value of the z-coordinate of the cell.
 * @param dx The size of the cell in x-direction.
 * @param dy The size of the cell in y-direction.
 * @param dz The size of the cell in z-direction.
 * @param vx Starting value of the vx-coordinate of the cell.
 * @param vy Starting value of the vy-coordinate of the cell.
 * @param vz Starting value of the vz-coordinate of the cell.
 * @param dvx The size of the cell in vx-direction.
 * @param dvy The size of the cell in vy-direction.
 * @param dvz The size of the cell in vz-direction.
 * @return The volume average of the distribution function in the given phase space cell.
 * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
 */
Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
   creal d_x = dx / (stP::nSpaceSamples-1);
   creal d_y = dy / (stP::nSpaceSamples-1);
   creal d_z = dz / (stP::nSpaceSamples-1);
   creal d_vx = dvx / (stP::nVelocitySamples-1);
   creal d_vy = dvy / (stP::nVelocitySamples-1);
   creal d_vz = dvz / (stP::nVelocitySamples-1);
   Real avg = 0.0;
   for (uint i=0; i<stP::nSpaceSamples; ++i)
      for (uint j=0; j<stP::nSpaceSamples; ++j)
	 for (uint k=0; k<stP::nSpaceSamples; ++k)
	    for (uint vi=0; vi<stP::nVelocitySamples; ++vi)
	       for (uint vj=0; vj<stP::nVelocitySamples; ++vj)
		  for (uint vk=0; vk<stP::nVelocitySamples; ++vk)
		     {
			avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
		     }
   return avg / pow(stP::nSpaceSamples, 3.0) / pow(stP::nVelocitySamples, 3.0);
}

/** Calculate parameters for the given spatial cell at the given time.
 * Here you need to set values for the following array indices:
 * CellParams::EX, CellParams::EY, CellParams::EZ, CellParams::BX, CellParams::BY, and CellParams::BZ.
 * 
 * The following array indices contain the coordinates of the "lower left corner" of the cell: 
 * CellParams::XCRD, CellParams::YCRD, and CellParams::ZCRD.
 * The cell size is given in the following array indices: CellParams::DX, CellParams::DY, and CellParams::DZ.
 * @param cellParams Array containing cell parameters.
 * @param t The current value of time. This is passed as a convenience. If you need more detailed information 
 * of the state of the simulation, you can read it from Parameters.
 */
void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   
   Real Bxavg, Byavg, Bzavg;
   Bxavg = Byavg = Bzavg = 0.0;
   Real d_x = dx / (stP::nSpaceSamples - 1);
   
   for (uint i=0; i<stP::nSpaceSamples; ++i)
      for (uint j=0; j<stP::nSpaceSamples; ++j)
	 for (uint k=0; k<stP::nSpaceSamples; ++k) {
	    Bxavg += ((x + i * d_x) < 0.0) ? stP::Bx[stP::LEFT] : stP::Bx[stP::RIGHT];
	    Byavg += ((x + i * d_x) < 0.0) ? stP::By[stP::LEFT] : stP::By[stP::RIGHT];
	    Bzavg += ((x + i * d_x) < 0.0) ? stP::Bz[stP::LEFT] : stP::Bz[stP::RIGHT];
	 }
   cuint nPts = pow(stP::nSpaceSamples, 3.0);
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::PERBX   ] = 0.0;
   cellParams[CellParams::PERBY   ] = 0.0;
   cellParams[CellParams::PERBZ   ] = 0.0;
   cellParams[CellParams::BGBX   ] = Bxavg / nPts;
   cellParams[CellParams::BGBY   ] = Byavg / nPts;
   cellParams[CellParams::BGBZ   ] = Bzavg / nPts;
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
