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
   RP::add("Riemann.rhoL", "Number density, left state (m^-3)", 0.0);
   RP::add("Riemann.rhoR", "Number density, right state (m^-3)", 0.0);
   RP::add("Riemann.TL", "Temperature, left state (K)", 0.0);
   RP::add("Riemann.TR", "Temperature, right state (K)", 0.0);
   RP::add("Riemann.VxL", "Bulk velocity x component, left state (m/s)", 0.0);
   RP::add("Riemann.VxR", "Bulk velocity x component, right state (m/s)", 0.0);
   RP::add("Riemann.VyL", "Bulk velocity y component, left state (m/s)", 0.0);
   RP::add("Riemann.VyR", "Bulk velocity y component, right state (m/s)", 0.0);
   RP::add("Riemann.VzL", "Bulk velocity z component, left state (m/s)", 0.0);
   RP::add("Riemann.VzR", "Bulk velocity z component, right state (m/s)", 0.0);
   RP::add("Riemann.BxL", "Magnetic field x component, left state (T)", 0.0);
   RP::add("Riemann.BxR", "Magnetic field x component, right state (T)", 0.0);
   RP::add("Riemann.ByL", "Magnetic field y component, left state (T)", 0.0);
   RP::add("Riemann.ByR", "Magnetic field y component, right state (T)", 0.0);
   RP::add("Riemann.BzL", "Magnetic field z component, left state (T)", 0.0);
   RP::add("Riemann.BzR", "Magnetic field z component, right state (T)", 0.0);
   RP::add("Riemann.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Riemann.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("Riemann.rhoL", RiP::rho[RiP::LEFT]);
   RP::get("Riemann.rhoR", RiP::rho[RiP::RIGHT]);
   RP::get("Riemann.TL", RiP::T[RiP::LEFT]);
   RP::get("Riemann.TR", RiP::T[RiP::RIGHT]);
   RP::get("Riemann.VxL", RiP::Vx[RiP::LEFT]);
   RP::get("Riemann.VxR", RiP::Vx[RiP::RIGHT]);
   RP::get("Riemann.VyL", RiP::Vy[RiP::LEFT]);
   RP::get("Riemann.VyR", RiP::Vy[RiP::RIGHT]);
   RP::get("Riemann.VzL", RiP::Vz[RiP::LEFT]);
   RP::get("Riemann.VzR", RiP::Vz[RiP::RIGHT]);
   RP::get("Riemann.BxL", RiP::Bx[RiP::LEFT]);
   RP::get("Riemann.BxR", RiP::Bx[RiP::RIGHT]);
   RP::get("Riemann.ByL", RiP::By[RiP::LEFT]);
   RP::get("Riemann.ByR", RiP::By[RiP::RIGHT]);
   RP::get("Riemann.BzL", RiP::Bz[RiP::LEFT]);
   RP::get("Riemann.BzR", RiP::Bz[RiP::RIGHT]);
   RP::get("Riemann.nSpaceSamples", RiP::nSpaceSamples);
   RP::get("Riemann.nVelocitySamples", RiP::nVelocitySamples);
   return true;
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

