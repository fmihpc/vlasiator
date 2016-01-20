/*
This file is part of Vlasiator.

Copyright 2011, 2012, 2015 Finnish Meteorological Institute
 * 
*/

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../object_wrapper.h"

#include "KHB.h"

namespace projects {
   using namespace std;
   KHB::KHB(): Project() { }
   KHB::~KHB() { }
   
   bool KHB::initialize(void) {return Project::initialize();}
   
   void KHB::addParameters() {
      typedef Readparameters RP;
      RP::add("KHB.rho1", "Number density, this->TOP state (m^-3)", 0.0);
      RP::add("KHB.rho2", "Number density, this->BOTTOM state (m^-3)", 0.0);
      RP::add("KHB.T1", "Temperature, this->TOP state (K)", 0.0);
      RP::add("KHB.T2", "Temperature, this->BOTTOM state (K)", 0.0);
      RP::add("KHB.Vx1", "Bulk velocity x component, this->TOP state (m/s)", 0.0);
      RP::add("KHB.Vx2", "Bulk velocity x component, this->BOTTOM state (m/s)", 0.0);
      RP::add("KHB.Vy1", "Bulk velocity y component, this->TOP state (m/s)", 0.0);
      RP::add("KHB.Vy2", "Bulk velocity y component, this->BOTTOM state (m/s)", 0.0);
      RP::add("KHB.Vz1", "Bulk velocity z component, this->TOP state (m/s)", 0.0);
      RP::add("KHB.Vz2", "Bulk velocity z component, this->BOTTOM state (m/s)", 0.0);
      RP::add("KHB.Bx1", "Magnetic field x component, this->TOP state (T)", 0.0);
      RP::add("KHB.Bx2", "Magnetic field x component, this->BOTTOM state (T)", 0.0);
      RP::add("KHB.By1", "Magnetic field y component, this->TOP state (T)", 0.0);
      RP::add("KHB.By2", "Magnetic field y component, this->BOTTOM state (T)", 0.0);
      RP::add("KHB.Bz1", "Magnetic field z component, this->TOP state (T)", 0.0);
      RP::add("KHB.Bz2", "Magnetic field z component, this->BOTTOM state (T)", 0.0);
      RP::add("KHB.lambda", "Initial perturbation wavelength (m)", 0.0);
      RP::add("KHB.amp", "Initial perturbation amplitude (m)", 0.0);
      RP::add("KHB.offset", "Boundaries offset from 0 (m)", 0.0);
      RP::add("KHB.transitionWidth", "Width of tanh transition for all changing values", 0.0);
      RP::add("KHB.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("KHB.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   }

   void KHB::getParameters() {
      Project::getParameters();
      typedef Readparameters RP;
      RP::get("KHB.rho1", this->rho[this->TOP]);
      RP::get("KHB.rho2", this->rho[this->BOTTOM]);
      RP::get("KHB.T1", this->T[this->TOP]);
      RP::get("KHB.T2", this->T[this->BOTTOM]);
      RP::get("KHB.Vx1", this->Vx[this->TOP]);
      RP::get("KHB.Vx2", this->Vx[this->BOTTOM]);
      RP::get("KHB.Vy1", this->Vy[this->TOP]);
      RP::get("KHB.Vy2", this->Vy[this->BOTTOM]);
      RP::get("KHB.Vz1", this->Vz[this->TOP]);
      RP::get("KHB.Vz2", this->Vz[this->BOTTOM]);
      RP::get("KHB.Bx1", this->Bx[this->TOP]);
      RP::get("KHB.Bx2", this->Bx[this->BOTTOM]);
      RP::get("KHB.By1", this->By[this->TOP]);
      RP::get("KHB.By2", this->By[this->BOTTOM]);
      RP::get("KHB.Bz1", this->Bz[this->TOP]);
      RP::get("KHB.Bz2", this->Bz[this->BOTTOM]);
      RP::get("KHB.lambda", this->lambda);
      RP::get("KHB.amp", this->amp);
      RP::get("KHB.offset", this->offset);
      RP::get("KHB.transitionWidth", this->transitionWidth);
      RP::get("KHB.nSpaceSamples", this->nSpaceSamples);
      RP::get("KHB.nVelocitySamples", this->nVelocitySamples);
   }
   
   
   Real KHB::profile(creal top, creal bottom, creal x, creal z) {
      if(top == bottom) {
         return top;
      }
      if(this->offset != 0.0) {
         return 0.5 * ((top-bottom) * (
         tanh((x + this->offset + this->amp * cos(2.0*M_PI*z/this->lambda))/this->transitionWidth) -
         tanh((x-(this->offset + this->amp * cos(2.0*M_PI*z/this->lambda)))/this->transitionWidth) -1) + top+bottom);
      } else {
         return 0.5 * ((top-bottom) * tanh(x/this->transitionWidth) + top+bottom);
      }
   }
   
   Real KHB::getDistribValue(creal& x, creal& z, creal& vx, creal& vy, creal& vz){
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      Real rho = profile(this->rho[this->BOTTOM], this->rho[this->TOP], x, z);
      Real T = profile(this->T[this->BOTTOM], this->T[this->TOP], x, z);
      Real Vx = profile(this->Vx[this->BOTTOM], this->Vx[this->TOP], x, z);
      Real Vy = profile(this->Vy[this->BOTTOM], this->Vy[this->TOP], x, z);
      Real Vz = profile(this->Vz[this->BOTTOM], this->Vz[this->TOP], x, z);
      
      return rho * pow(mass / (2.0 * M_PI * kb * T), 1.5) *
      exp(- mass * (pow(vx - Vx, 2.0) + pow(vy - Vy, 2.0) + pow(vz - Vz, 2.0)) / (2.0 * kb * T));
   }

   Real KHB::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const int& popID) {   
      creal d_x = dx / (this->nSpaceSamples-1);
      creal d_z = dz / (this->nSpaceSamples-1);
      creal d_vx = dvx / (this->nVelocitySamples-1);
      creal d_vy = dvy / (this->nVelocitySamples-1);
      creal d_vz = dvz / (this->nVelocitySamples-1);
      Real avg = 0.0;
      uint samples=0;

      Real middleValue=getDistribValue(x+0.5*dx, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz);
      if (middleValue < 0.000001*getObjectWrapper().particleSpecies[popID].sparseMinValue) {
         return middleValue; //abort, this will not be accepted anyway
      }
      
   //#pragma omp parallel for collapse(6) reduction(+:avg)
      for (uint i=0; i<this->nSpaceSamples; ++i)
         for (uint k=0; k<this->nSpaceSamples; ++k)
            for (uint vi=0; vi<this->nVelocitySamples; ++vi)
               for (uint vj=0; vj<this->nVelocitySamples; ++vj)
                  for (uint vk=0; vk<this->nVelocitySamples; ++vk){
                     avg +=getDistribValue(x+i*d_x, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz);
                     samples++;
                  }
      return avg / samples;
   }
   

   void KHB::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal z = cellParams[CellParams::ZCRD];
      creal dz = cellParams[CellParams::DZ];
      
      Real Bxavg, Byavg, Bzavg;
      Bxavg = Byavg = Bzavg = 0.0;
      Real d_x = dx / (this->nSpaceSamples - 1);
      Real d_z = dz / (this->nSpaceSamples - 1);
      for (uint i=0; i<this->nSpaceSamples; ++i)
         for (uint k=0; k<this->nSpaceSamples; ++k) {
            Bxavg += profile(this->Bx[this->BOTTOM], this->Bx[this->TOP], x+i*d_x, z+k*d_z);
            Byavg += profile(this->By[this->BOTTOM], this->By[this->TOP], x+i*d_x, z+k*d_z);
            Bzavg += profile(this->Bz[this->BOTTOM], this->Bz[this->TOP], x+i*d_x, z+k*d_z);
         }
      cuint nPts = pow(this->nSpaceSamples, 2.0);
      
      cellParams[CellParams::PERBX   ] = Bxavg / nPts;
      cellParams[CellParams::PERBY   ] = Byavg / nPts;
      cellParams[CellParams::PERBZ   ] = Bzavg / nPts;
   }

} // namespace projects
