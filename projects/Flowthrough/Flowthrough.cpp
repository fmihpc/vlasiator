/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"

#include "Flowthrough.h"

using namespace std;


/** Enumerates spatial density models Flowthrough project supports. 
 * In most cases you want to use 'Maxwellian'. However, test package 
 * uses 'SheetMaxwellian'.*/

enum DensityModel {
   Maxwellian,
   SheetMaxwellian
};

static DensityModel densityModel;

namespace projects {
   Flowthrough::Flowthrough(): TriAxisSearch() { }
   Flowthrough::~Flowthrough() { }
   
   bool Flowthrough::initialize(void) {
      return Project::initialize();
   }

   void Flowthrough::addParameters(){
      typedef Readparameters RP;
      RP::add("Flowthrough.emptyBox","Is the simulation domain empty initially?",false);
      RP::add("Flowthrough.rho", "Number density (m^-3)", 0.0);
      RP::add("Flowthrough.T", "Temperature (K)", 0.0);
      RP::add("Flowthrough.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("Flowthrough.By", "Magnetic field y component (T)", 0.0);
      RP::add("Flowthrough.Bz", "Magnetic field z component (T)", 0.0);
      RP::add("Flowthrough.VX0", "Initial bulk velocity in x-direction", 0.0);
      RP::add("Flowthrough.VY0", "Initial bulk velocity in y-direction", 0.0);
      RP::add("Flowthrough.VZ0", "Initial bulk velocity in z-direction", 0.0);
      RP::add("Flowthrough.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Flowthrough.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
      RP::add("Flowthrough.densityModel","Plasma density model, 'Maxwellian' or 'SheetMaxwellian'",string("Maxwellian"));
   }
   
   void Flowthrough::getParameters(){
      Project::getParameters();
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      typedef Readparameters RP;
      if (!RP::get("Flowthrough.emptyBox",emptyBox)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.rho", this->rho)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.T", this->T)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.Bx", this->Bx)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.By", this->By)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.Bz", this->Bz)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.VX0", this->V0[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.VY0", this->V0[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.VZ0", this->V0[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.nSpaceSamples", this->nSpaceSamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Flowthrough.nVelocitySamples", this->nVelocitySamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      string densityModelString;
      if (!RP::get("Flowthrough.densityModel",densityModelString)) {
         if (myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if (densityModelString == "Maxwellian") densityModel = Maxwellian;
      else if (densityModelString == "SheetMaxwellian") densityModel = SheetMaxwellian;
      else {
         if (myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: Unknown option value!" << endl;
         exit(1);
      }
   }

   Real Flowthrough::getDistribValue(creal& x,creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) const {
      Real rvalue = 0;
      switch (densityModel) {
       case Maxwellian:
         rvalue = rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5)
           * exp(- physicalconstants::MASS_PROTON * (  (vx-this->V0[0])*(vx-this->V0[0])
                                                       + (vy-this->V0[1])*(vy-this->V0[1])
                                                       + (vz-this->V0[2])*(vz-this->V0[2])
                                                    ) / (2.0 * physicalconstants::K_B * this->T));
         break;
       case SheetMaxwellian:
         rvalue = sqrt(x*x + y*y + z*z);
         if (rvalue <= +3e7) {
            rvalue = 4*rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5)
              * exp(- physicalconstants::MASS_PROTON * ((  vx-this->V0[0])*(vx-this->V0[0]) + (vy-this->V0[1])*(vy-this->V0[1])
                                                        + (vz-this->V0[2])*(vz-this->V0[2])) / (2.0 * physicalconstants::K_B * this->T));
         } else {
            rvalue = 0;
         }
         break;
      }
      
      return rvalue;
   }

   Real Flowthrough::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const int& popID) const {
      if (emptyBox == true) return 0.0;
      if((this->nSpaceSamples > 1) && (this->nVelocitySamples > 1)) {
         creal d_x = dx / (nSpaceSamples-1);
         creal d_y = dy / (nSpaceSamples-1);
         creal d_z = dz / (nSpaceSamples-1);
         creal d_vx = dvx / (nVelocitySamples-1);
         creal d_vy = dvy / (nVelocitySamples-1);
         creal d_vz = dvz / (nVelocitySamples-1);

         Real avg = 0.0;
         for (uint i=0; i<nSpaceSamples; ++i) for (uint j=0; j<nSpaceSamples; ++j) for (uint k=0; k<nSpaceSamples; ++k) {
            for (uint vi=0; vi<nVelocitySamples; ++vi) for (uint vj=0; vj<nVelocitySamples; ++vj) for (uint vk=0; vk<nVelocitySamples; ++vk) {
               avg += getDistribValue(x+i*d_x,y+j*d_y,z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
            }
         }
         return avg / (nSpaceSamples*nSpaceSamples*nSpaceSamples*nVelocitySamples*nVelocitySamples*nVelocitySamples);
      } else {
         return getDistribValue(x+0.5*dx,y+0.5*dy,z+0.5*dz,vx+0.5*dvx,vy+0.5*dvy,vz+0.5*dvz,dvx,dvy,dvz);
         
      }
   }

   void Flowthrough::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      cellParams[CellParams::PERBX] = this->Bx;
      cellParams[CellParams::PERBY] = this->By;
      cellParams[CellParams::PERBZ] = this->Bz;
   }

   void Flowthrough::setCellBackgroundField(spatial_cell::SpatialCell* cell) const {
      if (Parameters::propagateField == true) {
         ConstantField bgField;
         bgField.initialize(0,0,0); //bg bx, by,bz
         setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
      } else {
         ConstantField bgField;
         bgField.initialize(Bx,By,Bz); //bg bx, by,bz
         setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
      }
   }
   
   std::vector<std::array<Real, 3> > Flowthrough::getV0(
      creal x,
      creal y,
      creal z
   ) const {
      vector<std::array<Real, 3>> centerPoints;
      std::array<Real, 3> point {{this->V0[0], this->V0[1], this->V0[2]}};
      centerPoints.push_back(point);
      return centerPoints;
   }
   
} //namespace projects
