/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../object_wrapper.h"

#include "Alfven.h"

using namespace std;

namespace projects {
   Alfven::Alfven(): Project() { }
   Alfven::~Alfven() { }
   
   bool Alfven::initialize(void) {
      bool success = Project::initialize();

      Real norm = sqrt(this->Bx_guiding*this->Bx_guiding + this->By_guiding*this->By_guiding + this->Bz_guiding*this->Bz_guiding);
      this->Bx_guiding /= norm;
      this->By_guiding /= norm;
      this->By_guiding /= norm;
      this->ALPHA = atan(this->By_guiding/this->Bx_guiding);
      
      return success;
   } 

   void Alfven::addParameters() {
      typedef Readparameters RP;
      RP::add("Alfven.B0", "Guiding field value (T)", 1.0e-10);
      RP::add("Alfven.Bx_guiding", "Guiding field x component", 1);
      RP::add("Alfven.By_guiding", "Guiding field y component", 0);
      RP::add("Alfven.Bz_guiding", "Guiding field z component", 0);
      RP::add("Alfven.rho", "Number density (m^-3)", 1.0e8);
      RP::add("Alfven.Wavelength", "Wavelength (m)", 100000.0);
      RP::add("Alfven.Temperature", "Temperature (K)", 0.86456498092);
      RP::add("Alfven.A_mag", "Amplitude of the magnetic perturbation", 0.1);
      RP::add("Alfven.A_vel", "Amplitude of the velocity perturbation", 0.1);
      RP::add("Alfven.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Alfven.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   }

   void Alfven::getParameters(){
      Project::getParameters();
      
      typedef Readparameters RP;
      RP::get("Alfven.B0", this->B0);
      RP::get("Alfven.Bx_guiding", this->Bx_guiding);
      RP::get("Alfven.By_guiding", this->By_guiding);
      RP::get("Alfven.Bz_guiding", this->Bz_guiding);
      RP::get("Alfven.rho", this->DENSITY);
      RP::get("Alfven.Wavelength", this->WAVELENGTH);
      RP::get("Alfven.Temperature", this->TEMPERATURE);
      RP::get("Alfven.A_mag", this->A_MAG);
      RP::get("Alfven.A_vel", this->A_VEL);
      RP::get("Alfven.nSpaceSamples", this->nSpaceSamples);
      RP::get("Alfven.nVelocitySamples", this->nVelocitySamples);
   }

   /*Real calcPhaseSpaceDensity(creal& z,creal& x,creal& y,creal& dz,creal& dx,creal& dy,
               creal& vz,creal& vx,creal& vy,creal& dvz,creal& dvx,creal& dvy) {*/
   Real Alfven::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const int& popID) {
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;
      creal mu0 = physicalconstants::MU_0;
      creal ALFVEN_VEL = this->B0 / sqrt(mu0 * this->DENSITY * mass);

      creal ksi = (x * cos(this->ALPHA) + y * sin(this->ALPHA)) / this->WAVELENGTH;
      creal Vx = this->A_VEL * ALFVEN_VEL * sin(this->ALPHA) * sin(2.0 * M_PI * ksi);
      creal Vy = - this->A_VEL * ALFVEN_VEL * cos(this->ALPHA) * sin(2.0 * M_PI * ksi);
      creal Vz = - this->A_VEL * ALFVEN_VEL * cos(2.0 * M_PI * ksi);
   
      creal den = this->DENSITY * pow(mass / (2.0 * M_PI * kb * this->TEMPERATURE), 1.5) *
      exp(- mass * (pow(vx - Vx, 2.0) + pow(vy - Vy, 2.0) + pow(vz - Vz, 2.0)) / (2.0 * kb * this->TEMPERATURE));
   return den;
   }
   
   Real Alfven::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const int& popID) {
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
            avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz,popID);
         }
      return avg / pow(this->nSpaceSamples, 3.0) / pow(this->nVelocitySamples, 3.0);
   }
   
   void Alfven::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      
      Real dBxavg, dByavg, dBzavg;
      dBxavg = dByavg = dBzavg = 0.0;
      Real d_x = dx / (this->nSpaceSamples - 1);
      Real d_y = dy / (this->nSpaceSamples - 1);

      for (uint i=0; i<this->nSpaceSamples; ++i)
         for (uint j=0; j<this->nSpaceSamples; ++j)
      for (uint k=0; k<this->nSpaceSamples; ++k) {
         Real ksi = ((x + i * d_x)  * cos(this->ALPHA) + (y + j * d_y) * sin(this->ALPHA)) / this->WAVELENGTH;
         dBxavg += sin(2.0 * M_PI * ksi);
         dByavg += sin(2.0 * M_PI * ksi);
         dBzavg += cos(2.0 * M_PI * ksi);
      }
      cuint nPts = pow(this->nSpaceSamples, 3.0);
      
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      //Field below could laso be set as background field
      cellParams[CellParams::PERBX   ] = this->B0 * cos(this->ALPHA) - this->A_MAG * this->B0 * sin(this->ALPHA) * dBxavg / nPts;
      cellParams[CellParams::PERBY   ] = this->B0 * sin(this->ALPHA) + this->A_MAG * this->B0 * cos(this->ALPHA) * dByavg / nPts;
      cellParams[CellParams::PERBZ   ] = this->B0 * this->A_MAG * dBzavg / nPts;
   }
} // namespace projects
