/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute
*/

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "Shocktest.h"
#include "../../spatial_cell.hpp"
#include "../../common.h"
#include "../project.h"
#include "../../parameters.h"
#include "../../readparameters.h"
#include "../../vlasovmover.h"

using namespace std;
namespace projects {

   Shocktest::Shocktest() : IsotropicMaxwellian::IsotropicMaxwellian() {} // Constructor
   Shocktest::~Shocktest() {} // Destructor

   
   bool Shocktest::initialize(void) { return true; }
   
   void Shocktest::addParameters(){
      typedef Readparameters RP;
      RP::add("Shocktest.rho1", "Number density, left state (m^-3)", 0.0);
      RP::add("Shocktest.rho2", "Number density, right state (m^-3)", 0.0);
      RP::add("Shocktest.T1", "Temperature, left state (K)", 0.0);
      RP::add("Shocktest.T2", "Temperature, right state (K)", 0.0);
      RP::add("Shocktest.Vx1", "Bulk velocity x component, left state (m/s)", 0.0);
      RP::add("Shocktest.Vx2", "Bulk velocity x component, right state (m/s)", 0.0);
      RP::add("Shocktest.Vy1", "Bulk velocity y component, left state (m/s)", 0.0);
      RP::add("Shocktest.Vy2", "Bulk velocity y component, right state (m/s)", 0.0);
      RP::add("Shocktest.Vz1", "Bulk velocity z component, left state (m/s)", 0.0);
      RP::add("Shocktest.Vz2", "Bulk velocity z component, right state (m/s)", 0.0);
      RP::add("Shocktest.Bx1", "Magnetic field x component, left state (T)", 0.0);
      RP::add("Shocktest.Bx2", "Magnetic field x component, right state (T)", 0.0);
      RP::add("Shocktest.By1", "Magnetic field y component, left state (T)", 0.0);
      RP::add("Shocktest.By2", "Magnetic field y component, right state (T)", 0.0);
      RP::add("Shocktest.Bz1", "Magnetic field z component, left state (T)", 0.0);
      RP::add("Shocktest.Bz2", "Magnetic field z component, right state (T)", 0.0);
      RP::add("Shocktest.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Shocktest.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   }
   
   void Shocktest::getParameters(){
      this->rho[2] = {NAN};
      this->T[2] = {NAN};
      this->Vx[2] = {NAN};
      this->Vy[2] = {NAN};
      this->Vz[2] = {NAN};
      this->Bx[2] = {NAN};
      this->By[2] = {NAN};
      this->Bz[2] = {NAN};
      this->nSpaceSamples = 0;
      this->nVelocitySamples = 0;

      typedef Readparameters RP;
      RP::get("Shocktest.rho1", this->rho[this->LEFT]);
      RP::get("Shocktest.rho2", this->rho[this->RIGHT]);
      RP::get("Shocktest.T1", this->T[this->LEFT]);
      RP::get("Shocktest.T2", this->T[this->RIGHT]);
      RP::get("Shocktest.Vx1", this->Vx[this->LEFT]);
      RP::get("Shocktest.Vx2", this->Vx[this->RIGHT]);
      RP::get("Shocktest.Vy1", this->Vy[this->LEFT]);
      RP::get("Shocktest.Vy2", this->Vy[this->RIGHT]);
      RP::get("Shocktest.Vz1", this->Vz[this->LEFT]);
      RP::get("Shocktest.Vz2", this->Vz[this->RIGHT]);
      RP::get("Shocktest.Bx1", this->Bx[this->LEFT]);
      RP::get("Shocktest.Bx2", this->Bx[this->RIGHT]);
      RP::get("Shocktest.By1", this->By[this->LEFT]);
      RP::get("Shocktest.By2", this->By[this->RIGHT]);
      RP::get("Shocktest.Bz1", this->Bz[this->LEFT]);
      RP::get("Shocktest.Bz2", this->Bz[this->RIGHT]);
      RP::get("Shocktest.nSpaceSamples", this->nSpaceSamples);
      RP::get("Shocktest.nVelocitySamples", this->nVelocitySamples);
   }
   
   Real Shocktest::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      creal mass = 1.67262171e-27; // m_p in kg
      creal k = 1.3806505e-23; // Boltzmann
      //  creal mu0 = 1.25663706144e-6; // mu_0
      //  creal q = 1.60217653e-19; // q_i
      //  creal gamma = 5./3.;
      
      cint side = (x < 0.0) ? this->LEFT : this->RIGHT;

      // Disable compiler warnings: (unused variables but the function is inherited)
      (void)y; (void)z; (void)dvx; (void)dvy; (void)dvz;
      
      return this->rho[side] * pow(mass / (2.0 * M_PI * k * this->T[side]), 1.5) *
      exp(- mass * (pow(vx - this->Vx[side], 2.0) + pow(vy - this->Vy[side], 2.0) + pow(vz - this->Vz[side], 2.0)) / (2.0 * k * this->T[side]));
   }


   /** Returns the center coordinates of the maxwellian distribution
   @ param x The x coordinate of the given spatial cell
   @ param y The x coordinate of the given spatial cell
   @ param z The x coordinate of the given spatial cell
   @ param component The component to be returned
   */

   Real Shocktest::getV0( 
      creal x,
      creal y,
      creal z,
      cuint component
   ) {
      Real V0;
      if( component == 0 ) {
         V0 = this->Vx[this->LEFT];
      } else if( component == 1 ) {
         V0 = this->Vy[this->LEFT];
      } else if( component == 2 ) {
         V0 = this->Vz[this->LEFT];
      }
      return V0;
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
   Real Shocktest::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
      creal d_x = dx / (this->nSpaceSamples-1);
      creal d_y = dy / (this->nSpaceSamples-1);
      creal d_z = dz / (this->nSpaceSamples-1);
      creal d_vx = dvx / (this->nVelocitySamples-1);
      creal d_vy = dvy / (this->nVelocitySamples-1);
      creal d_vz = dvz / (this->nVelocitySamples-1);
      Real avg = 0.0;
      for (uint i=0; i<this->nSpaceSamples; ++i)
         for (uint j=0; j<this->nSpaceSamples; ++j)
   	 for (uint k=0; k<this->nSpaceSamples; ++k)
   	    for (uint vi=0; vi<this->nVelocitySamples; ++vi)
   	       for (uint vj=0; vj<this->nVelocitySamples; ++vj)
   		  for (uint vk=0; vk<this->nVelocitySamples; ++vk)
   		     {
   			avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
   		     }
      return avg / pow(this->nSpaceSamples, 3.0) / pow(this->nVelocitySamples, 3.0);
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
   void Shocktest::calcCellParameters(Real* cellParams,creal& t) {
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      
      Real Bxavg, Byavg, Bzavg;
      Bxavg = Byavg = Bzavg = 0.0;
      Real d_x = dx / (this->nSpaceSamples - 1);
      
      for (uint i=0; i<this->nSpaceSamples; ++i)
         for (uint j=0; j<this->nSpaceSamples; ++j)
   	 for (uint k=0; k<this->nSpaceSamples; ++k) {
   	    Bxavg += ((x + i * d_x) < 0.0) ? this->Bx[this->LEFT] : this->Bx[this->RIGHT];
   	    Byavg += ((x + i * d_x) < 0.0) ? this->By[this->LEFT] : this->By[this->RIGHT];
   	    Bzavg += ((x + i * d_x) < 0.0) ? this->Bz[this->LEFT] : this->Bz[this->RIGHT];
   	 }
      cuint nPts = pow(this->nSpaceSamples, 3.0);
      
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
   
   void Shocktest::setCell(SpatialCell* cell) {
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
} // Namespace projects
