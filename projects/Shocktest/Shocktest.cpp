/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute
*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"


#include "Shocktest.h"
#include "../../spatial_cell.hpp"
#include "../../common.h"
#include "../project.h"
#include "../../parameters.h"
#include "../../readparameters.h"
#include "../../vlasovmover.h"

using namespace std;
using namespace spatial_cell;

namespace projects {

   Shocktest::Shocktest() : TriAxisSearch() {} // Constructor
   Shocktest::~Shocktest() {} // Destructor

   
   bool Shocktest::initialize(void) { return Project::initialize(); }
   
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
      Project::getParameters();
      this->rho[this->LEFT] = {NAN};
      this->T[this->LEFT] = {NAN};
      this->Vx[this->LEFT] = {NAN};
      this->Vy[this->LEFT] = {NAN};
      this->Vz[this->LEFT] = {NAN};
      this->Bx[this->LEFT] = {NAN};
      this->By[this->LEFT] = {NAN};
      this->Bz[this->LEFT] = {NAN};
      this->rho[this->RIGHT] = {NAN};
      this->T[this->RIGHT] = {NAN};
      this->Vx[this->RIGHT] = {NAN};
      this->Vy[this->RIGHT] = {NAN};
      this->Vz[this->RIGHT] = {NAN};
      this->Bx[this->RIGHT] = {NAN};
      this->By[this->RIGHT] = {NAN};
      this->Bz[this->RIGHT] = {NAN};
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
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      
      cint side = (x < 0.0) ? this->LEFT : this->RIGHT;

      // Disable compiler warnings: (unused variables but the function is inherited)
      (void)y; (void)z; (void)dvx; (void)dvy; (void)dvz;
      
      return this->rho[side] * pow(mass / (2.0 * M_PI * kb * this->T[side]), 1.5) *
      exp(- mass * (pow(vx - this->Vx[side], 2.0) + pow(vy - this->Vy[side], 2.0) + pow(vz - this->Vz[side], 2.0)) / (2.0 * kb * this->T[side]));
   }


   /** Returns the center coordinates of the maxwellian distribution
   @ param x The x coordinate of the given spatial cell
   @ param y The x coordinate of the given spatial cell
   @ param z The x coordinate of the given spatial cell
   */

   vector<std::array<Real, 3>> Shocktest::getV0(
      creal x,
      creal y,
      creal z
   ) {
      vector<std::array<Real, 3>> centerPoints;
      cint side = (x < 0.0) ? this->LEFT : this->RIGHT;
      std::array<Real, 3> V0 {{this->Vx[side], this->Vy[side], this->Vz[side]}};
      centerPoints.push_back(V0);
      return centerPoints;
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
   Real Shocktest::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const int& popID) {   
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
   void Shocktest::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
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
      cellParams[CellParams::PERBX   ] = Bxavg / nPts;
      cellParams[CellParams::PERBY   ] = Byavg / nPts;
      cellParams[CellParams::PERBZ   ] = Bzavg / nPts;
   }
   

   /*! Base class sets zero background field */
   void Shocktest::setCellBackgroundField(SpatialCell* cell) {
      ConstantField bgField;
      bgField.initialize(0,0,0); //bg bx, by,bz
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
   }

} // Namespace projects
