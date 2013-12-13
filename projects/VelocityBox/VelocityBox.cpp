/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"

#include "VelocityBox.h"

using namespace std;

namespace projects {
   VelocityBox::VelocityBox(): Project() { }
   VelocityBox::~VelocityBox() { }


   bool VelocityBox::initialize(void) {return true;}

   void VelocityBox::addParameters(){
      typedef Readparameters RP;
      RP::add("VelocityBox.rho", "Number density in full 6 dimensions (m^-6 s^3)", 0.0);
      RP::add("VelocityBox.Vx1", "Box min x (m/s)", 0.0);
      RP::add("VelocityBox.Vx2", "Box max x (m/s)", 0.0);
      RP::add("VelocityBox.Vy1", "Box min y (m/s)", 0.0);
      RP::add("VelocityBox.Vy2", "Box max y (m/s)", 0.0);
      RP::add("VelocityBox.Vz1", "Box min z (m/s)", 0.0);
      RP::add("VelocityBox.Vz2", "Box max z (m/s)", 0.0);
      RP::add("VelocityBox.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("VelocityBox.By", "Magnetic field y component (T)", 0.0);
      RP::add("VelocityBox.Bz", "Magnetic field z component (T)", 0.0);
   }

   void VelocityBox::getParameters(){
      typedef Readparameters RP;
      RP::get("VelocityBox.rho", this->rho);
      RP::get("VelocityBox.Vx1", this->Vx[0]);
      RP::get("VelocityBox.Vx2", this->Vx[1]);
      RP::get("VelocityBox.Vy1", this->Vy[0]);
      RP::get("VelocityBox.Vy2", this->Vy[1]);
      RP::get("VelocityBox.Vz1", this->Vz[0]);
      RP::get("VelocityBox.Vz2", this->Vz[1]);
      RP::get("VelocityBox.Bx", this->Bx);
      RP::get("VelocityBox.By", this->By);
      RP::get("VelocityBox.Bz", this->Bz);
   }

  Real VelocityBox::getDistribValue(creal& vx, creal& vy, creal& vz){
     if (vx >= this->Vx[0] && vx <= this->Vx[1] &&
	 vy >= this->Vy[0] && vy <= this->Vy[1] &&
	 vz >= this->Vz[0] && vz <= this->Vz[1])
       return this->rho;
   }



   Real VelocityBox::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
     return getDistribValue(vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz);
   }


  
   void VelocityBox::calcCellParameters(Real* cellParams,creal& t) {
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      cellParams[CellParams::PERBX   ] = 0.0;
      cellParams[CellParams::PERBY   ] = 0.0;
      cellParams[CellParams::PERBZ   ] = 0.0;
   }

   void VelocityBox::setCellBackgroundField(SpatialCell* cell) {
     ConstantField bgField;
      bgField.initialize(this->Bx,
                         this->By,
                         this->Bz);
      
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
   }
   
}// namespace projects
