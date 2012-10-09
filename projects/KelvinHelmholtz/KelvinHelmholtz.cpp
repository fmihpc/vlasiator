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

typedef kelvinHelmholtzParameters KHP;
Real KHP::rho[2] = {NAN};
Real KHP::T[2] = {NAN};
Real KHP::Vx[2] = {NAN};
Real KHP::Vy[2] = {NAN};
Real KHP::Vz[2] = {NAN};
Real KHP::Bx[2] = {NAN};
Real KHP::By[2] = {NAN};
Real KHP::Bz[2] = {NAN};
Real KHP::lambda = 0;
Real KHP::amp = NAN;
Real KHP::offset = 0;
Real KHP::transitionWidth = 0;
uint KHP::nSpaceSamples = 0;
uint KHP::nVelocitySamples = 0;

bool initializeProject(void) {return true;}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("KelvinHelmholtz.rho1", "Number density, KHP::TOP state (m^-3)", 0.0);
   RP::add("KelvinHelmholtz.rho2", "Number density, KHP::BOTTOM state (m^-3)", 0.0);
   RP::add("KelvinHelmholtz.T1", "Temperature, KHP::TOP state (K)", 0.0);
   RP::add("KelvinHelmholtz.T2", "Temperature, KHP::BOTTOM state (K)", 0.0);
   RP::add("KelvinHelmholtz.Vx1", "Bulk velocity x component, KHP::TOP state (m/s)", 0.0);
   RP::add("KelvinHelmholtz.Vx2", "Bulk velocity x component, KHP::BOTTOM state (m/s)", 0.0);
   RP::add("KelvinHelmholtz.Vy1", "Bulk velocity y component, KHP::TOP state (m/s)", 0.0);
   RP::add("KelvinHelmholtz.Vy2", "Bulk velocity y component, KHP::BOTTOM state (m/s)", 0.0);
   RP::add("KelvinHelmholtz.Vz1", "Bulk velocity z component, KHP::TOP state (m/s)", 0.0);
   RP::add("KelvinHelmholtz.Vz2", "Bulk velocity z component, KHP::BOTTOM state (m/s)", 0.0);
   RP::add("KelvinHelmholtz.Bx1", "Magnetic field x component, KHP::TOP state (T)", 0.0);
   RP::add("KelvinHelmholtz.Bx2", "Magnetic field x component, KHP::BOTTOM state (T)", 0.0);
   RP::add("KelvinHelmholtz.By1", "Magnetic field y component, KHP::TOP state (T)", 0.0);
   RP::add("KelvinHelmholtz.By2", "Magnetic field y component, KHP::BOTTOM state (T)", 0.0);
   RP::add("KelvinHelmholtz.Bz1", "Magnetic field z component, KHP::TOP state (T)", 0.0);
   RP::add("KelvinHelmholtz.Bz2", "Magnetic field z component, KHP::BOTTOM state (T)", 0.0);
   RP::add("KelvinHelmholtz.lambda", "Initial perturbation wavelength (m)", 0.0);
   RP::add("KelvinHelmholtz.amp", "Initial perturbation amplitude (m)", 0.0);
   RP::add("KelvinHelmholtz.offset", "Boundaries offset from 0 (m)", 0.0);
   RP::add("KelvinHelmholtz.transitionWidth", "Width of tanh transition for all changing values", 0.0);
   RP::add("KelvinHelmholtz.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("KelvinHelmholtz.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("KelvinHelmholtz.rho1", KHP::rho[KHP::TOP]);
   RP::get("KelvinHelmholtz.rho2", KHP::rho[KHP::BOTTOM]);
   RP::get("KelvinHelmholtz.T1", KHP::T[KHP::TOP]);
   RP::get("KelvinHelmholtz.T2", KHP::T[KHP::BOTTOM]);
   RP::get("KelvinHelmholtz.Vx1", KHP::Vx[KHP::TOP]);
   RP::get("KelvinHelmholtz.Vx2", KHP::Vx[KHP::BOTTOM]);
   RP::get("KelvinHelmholtz.Vy1", KHP::Vy[KHP::TOP]);
   RP::get("KelvinHelmholtz.Vy2", KHP::Vy[KHP::BOTTOM]);
   RP::get("KelvinHelmholtz.Vz1", KHP::Vz[KHP::TOP]);
   RP::get("KelvinHelmholtz.Vz2", KHP::Vz[KHP::BOTTOM]);
   RP::get("KelvinHelmholtz.Bx1", KHP::Bx[KHP::TOP]);
   RP::get("KelvinHelmholtz.Bx2", KHP::Bx[KHP::BOTTOM]);
   RP::get("KelvinHelmholtz.By1", KHP::By[KHP::TOP]);
   RP::get("KelvinHelmholtz.By2", KHP::By[KHP::BOTTOM]);
   RP::get("KelvinHelmholtz.Bz1", KHP::Bz[KHP::TOP]);
   RP::get("KelvinHelmholtz.Bz2", KHP::Bz[KHP::BOTTOM]);
   RP::get("KelvinHelmholtz.lambda", KHP::lambda);
   RP::get("KelvinHelmholtz.amp", KHP::amp);
   RP::get("KelvinHelmholtz.offset", KHP::offset);
   RP::get("KelvinHelmholtz.transitionWidth", KHP::transitionWidth);
   RP::get("KelvinHelmholtz.nSpaceSamples", KHP::nSpaceSamples);
   RP::get("KelvinHelmholtz.nVelocitySamples", KHP::nVelocitySamples);
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

Real profile(creal top, creal bottom, creal x, creal z) {
   if(top == bottom) {
      return top;
   }
   if(KHP::offset != 0.0) {
      return 0.5 * ((top-bottom) * (
	 tanh((z + KHP::offset + KHP::amp * cos(2.0*M_PI*x/KHP::lambda))/KHP::transitionWidth) -
	 tanh((z-(KHP::offset + KHP::amp * cos(2.0*M_PI*x/KHP::lambda)))/KHP::transitionWidth) -1) + top+bottom);
   } else {
      return 0.5 * ((top-bottom) * tanh(z/KHP::transitionWidth) + top+bottom);
   }
}

Real getDistribValue(creal& x, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal mass = 1.67262171e-27; // m_p in kg
   creal k = 1.3806505e-23; // Boltzmann
   //  creal mu0 = 1.25663706144e-6; // mu_0
   //  creal q = 1.60217653e-19; // q_i
   //  creal gamma = 5./3.;
   Real rho = profile(KHP::rho[KHP::TOP], KHP::rho[KHP::TOP], x, z);
   Real T = profile(KHP::T[KHP::BOTTOM], KHP::T[KHP::TOP], x, z);
   Real Vx = profile(KHP::Vx[KHP::BOTTOM], KHP::Vx[KHP::TOP], x, z);
   Real Vy = profile(KHP::Vy[KHP::BOTTOM], KHP::Vy[KHP::TOP], x, z);
   Real Vz = profile(KHP::Vz[KHP::BOTTOM], KHP::Vz[KHP::TOP], x, z);
   
   return rho * pow(mass / (2.0 * M_PI * k * T), 1.5) *
   exp(- mass * (pow(vx - Vx, 2.0) + pow(vy - Vy, 2.0) + pow(vz - Vz, 2.0)) / (2.0 * k * T));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
   creal d_x = dx / (KHP::nSpaceSamples-1);
   creal d_z = dz / (KHP::nSpaceSamples-1);
   creal d_vx = dvx / (KHP::nVelocitySamples-1);
   creal d_vy = dvy / (KHP::nVelocitySamples-1);
   creal d_vz = dvz / (KHP::nVelocitySamples-1);
   Real avg = 0.0;
//#pragma omp parallel for collapse(6) reduction(+:avg)
   for (uint i=0; i<KHP::nSpaceSamples; ++i)
      for (uint k=0; k<KHP::nSpaceSamples; ++k)
	 for (uint vi=0; vi<KHP::nVelocitySamples; ++vi)
	    for (uint vj=0; vj<KHP::nVelocitySamples; ++vj)
	       for (uint vk=0; vk<KHP::nVelocitySamples; ++vk)
		  {
		     avg += getDistribValue(x+i*d_x, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
		  }
   return avg / pow(KHP::nSpaceSamples, 2.0) / pow(KHP::nVelocitySamples, 3.0);
}

bool cellParametersChanged(creal& t) {return false;}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   creal z = cellParams[CellParams::ZCRD];
   creal dz = cellParams[CellParams::DZ];
   
   Real Bxavg, Byavg, Bzavg;
   Bxavg = Byavg = Bzavg = 0.0;
   Real d_x = dx / (KHP::nSpaceSamples - 1);
   Real d_z = dz / (KHP::nSpaceSamples - 1);
   for (uint i=0; i<KHP::nSpaceSamples; ++i)
      for (uint k=0; k<KHP::nSpaceSamples; ++k) {
	 Bxavg += profile(KHP::Bx[KHP::BOTTOM], KHP::Bx[KHP::TOP], x+i*d_x, z+k*d_z);
	 Byavg += profile(KHP::By[KHP::BOTTOM], KHP::By[KHP::TOP], x+i*d_x, z+k*d_z);
	 Bzavg += profile(KHP::Bz[KHP::BOTTOM], KHP::Bz[KHP::TOP], x+i*d_x, z+k*d_z);
      }
   cuint nPts = pow(KHP::nSpaceSamples, 2.0);
   
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

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

