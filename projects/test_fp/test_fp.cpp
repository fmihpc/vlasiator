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

#include "test_fp.h"

enum cases {BXCASE,BYCASE,BZCASE};

using namespace std;

namespace projects {
   test_fp::test_fp(): Project() { }
   test_fp::~test_fp() { }


   /*typedef test_fpParameters tfP;
   Real this->B0 = NAN;
   Real this->DENSITY = NAN;
   Real this->TEMPERATURE = NAN;
   Real this->ALPHA = NAN;
   int  this->CASE = 5;
   bool this->shear = false;
   */

   bool test_fp::initialize(void) {
      this->ALPHA *= M_PI / 4.0;
      return true;
   }

   void test_fp::addParameters(void){
      typedef Readparameters RP;
      RP::add("test_fp.B0", "Magnetic field value in the non-zero patch (T)", 1.0e-9);
      RP::add("test_fp.rho", "Number density (m^-3)", 1.0e7);
      RP::add("test_fp.Temperature", "Temperature (K)", 1.0e-6);
      RP::add("test_fp.angle", "Orientation of the propagation expressed in pi/4", 0.0);
      RP::add("test_fp.Bdirection", "Direction of the magnetic field (0:x, 1:y, 2:z)", 0);
      RP::add("test_fp.shear", "Add a shear (if false, V=0.5 everywhere).", true);
   }

   void test_fp::getParameters(void){
      typedef Readparameters RP;
      RP::get("test_fp.B0", this->B0);
      RP::get("test_fp.rho", this->DENSITY);
      RP::get("test_fp.Temperature", this->TEMPERATURE);
      RP::get("test_fp.angle", this->ALPHA);
      RP::get("test_fp.Bdirection", this->CASE);
      RP::get("test_fp.shear", this->shear);
   }

   Real test_fp::sign(creal value)
   {
      if(abs(value) < 1e-5) return 0.0;
      else return value / abs(value);
   }


   Real test_fp::calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
      Real VX,VY,VZ;
      if (this->shear == true)
      {
         Real ksi,eta;
         switch (this->CASE) {
            case BXCASE:
               ksi = ((y + 0.5 * dy)  * cos(this->ALPHA) + (z + 0.5 * dz) * sin(this->ALPHA)) / (2.0 * sqrt(2.0));
               eta = (-(y + 0.5 * dy)  * sin(this->ALPHA) + (z + 0.5 * dz) * cos(this->ALPHA)) / (2.0 * sqrt(2.0));
               VX = 0.0;
               VY = sign(cos(this->ALPHA)) * 0.5 + 0.1*cos(this->ALPHA) * sin(2.0 * M_PI * eta);
               VZ = sign(sin(this->ALPHA)) * 0.5 + 0.1*sin(this->ALPHA) * sin(2.0 * M_PI * eta);
               break;
            case BYCASE:
               ksi = ((z + 0.5 * dz)  * cos(this->ALPHA) + (x + 0.5 * dx) * sin(this->ALPHA)) / (2.0 * sqrt(2.0));
               eta = (-(z + 0.5 * dz)  * sin(this->ALPHA) + (x + 0.5 * dx) * cos(this->ALPHA)) / (2.0 * sqrt(2.0));
               VX = sign(sin(this->ALPHA)) * 0.5 + 0.1*sin(this->ALPHA) * sin(2.0 * M_PI * eta);
               VY = 0.0;
               VZ = sign(cos(this->ALPHA)) * 0.5 + 0.1*cos(this->ALPHA) * sin(2.0 * M_PI * eta);
               break;
               case BZCASE:
               ksi = ((x + 0.5 * dx)  * cos(this->ALPHA) + (y + 0.5 * dy) * sin(this->ALPHA)) / (2.0 * sqrt(2.0));
               eta = (-(x + 0.5 * dx)  * sin(this->ALPHA) + (y + 0.5 * dy) * cos(this->ALPHA)) / (2.0 * sqrt(2.0));
               VX = sign(cos(this->ALPHA)) * 0.5 + 0.1*cos(this->ALPHA) * sin(2.0 * M_PI * eta);
               VY = sign(sin(this->ALPHA)) * 0.5 + 0.1*sin(this->ALPHA) * sin(2.0 * M_PI * eta);
               VZ = 0.0;
               break;
         }
      } else {
         switch (this->CASE) {
            case BXCASE:
               VX = 0.0;
               VY = cos(this->ALPHA) * 0.5;
               VZ = sin(this->ALPHA) * 0.5; 
               break;
            case BYCASE:
               VX = sin(this->ALPHA) * 0.5;
               VY = 0.0;
               VZ = cos(this->ALPHA) * 0.5;
               break;
            case BZCASE:
               VX = cos(this->ALPHA) * 0.5;
               VY = sin(this->ALPHA) * 0.5;
               VZ = 0.0;
               break;
         }
      }
      
      creal VX2 = (vx+0.5*dvx-VX)*(vx+0.5*dvx-VX);
      creal VY2 = (vy+0.5*dvy-VY)*(vy+0.5*dvy-VY);
      creal VZ2 = (vz+0.5*dvz-VZ)*(vz+0.5*dvz-VZ);
      
      creal CONST = physicalconstants::MASS_PROTON / 2.0 / physicalconstants::K_B / this->TEMPERATURE;
      Real NORM = (physicalconstants::MASS_PROTON / 2.0 / M_PI / physicalconstants::K_B / this->TEMPERATURE);
      NORM = this->DENSITY * pow(NORM,1.5);
      
      return NORM*exp(-CONST*(VX2+VY2+VZ2));
   }

   void test_fp::calcCellParameters(Real* cellParams,creal& t) {
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      cellParams[CellParams::PERBX   ] = 0.0;
      cellParams[CellParams::PERBY   ] = 0.0;
      cellParams[CellParams::PERBZ   ] = 0.0;
      
      typedef Parameters P;
      creal x = cellParams[CellParams::XCRD] + 0.5 * cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD] + 0.5 * cellParams[CellParams::DY];
      creal z = cellParams[CellParams::ZCRD] + 0.5 * cellParams[CellParams::DZ];
      
      switch (this->CASE) {
      case BXCASE:
         if (y >= -0.2 && y <= 0.2)
      if (z >= -0.2 && z <= 0.2)
         cellParams[CellParams::PERBX] = this->B0;
         break;
      case BYCASE:
         if (x >= -0.2 && x <= 0.2)
      if (z >= -0.2 && z <= 0.2)
         cellParams[CellParams::PERBY] = this->B0;
         break;
      case BZCASE:
         if (x >= -0.2 && x <= 0.2)
      if (y >= -0.2 && y <= 0.2)
         cellParams[CellParams::PERBZ] = this->B0;
         break;
      }
   }

}// namespace projects