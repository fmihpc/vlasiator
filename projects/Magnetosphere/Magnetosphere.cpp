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

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"

#include "Magnetosphere.h"

using namespace std;

namespace projects {
   Magnetosphere::Magnetosphere(): IsotropicMaxwellian() { }
   Magnetosphere::~Magnetosphere() { }
   
   void Magnetosphere::addParameters() {
      typedef Readparameters RP;
      RP::add("Magnetosphere.rho", "Number density (m^-3)", 0.0);
      RP::add("Magnetosphere.T", "Temperature (K)", 0.0);
      RP::add("Magnetosphere.VX0", "Initial bulk velocity in x-direction", 0.0);
      RP::add("Magnetosphere.VY0", "Initial bulk velocity in y-direction", 0.0);
      RP::add("Magnetosphere.VZ0", "Initial bulk velocity in z-direction", 0.0);
      RP::add("Magnetosphere.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Magnetosphere.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   }
   
   void Magnetosphere::getParameters(){
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      typedef Readparameters RP;
      if(!RP::get("Magnetosphere.rho", this->rho)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.T", this->T)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.VX0", this->V0[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.VY0", this->V0[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.VZ0", this->V0[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.nSpaceSamples", this->nSpaceSamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.nVelocitySamples", this->nVelocitySamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
   }
   
   bool Magnetosphere::initialize() {
      return true;
   }

   Real Magnetosphere::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      if((this->nSpaceSamples > 1) && (this->nVelocitySamples > 1)) {
         creal d_x = dx / (this->nSpaceSamples-1);
         creal d_y = dy / (this->nSpaceSamples-1);
         creal d_z = dz / (this->nSpaceSamples-1);
         creal d_vx = dvx / (this->nVelocitySamples-1);
         creal d_vy = dvy / (this->nVelocitySamples-1);
         creal d_vz = dvz / (this->nVelocitySamples-1);
         
         Real avg = 0.0;
         // #pragma omp parallel for collapse(6) reduction(+:avg)
         // WARNING No threading here if calling functions are already threaded
         for (uint i=0; i<this->nSpaceSamples; ++i)
            for (uint j=0; j<this->nSpaceSamples; ++j)
               for (uint k=0; k<this->nSpaceSamples; ++k)
                  for (uint vi=0; vi<this->nVelocitySamples; ++vi)
                     for (uint vj=0; vj<this->nVelocitySamples; ++vj)
                        for (uint vk=0; vk<this->nVelocitySamples; ++vk) {
                           avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
                        }
                        return avg /
                        (this->nSpaceSamples*this->nSpaceSamples*this->nSpaceSamples) /
                        (this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
      } else {
         return getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, dvx, dvy, dvz);
      }
      
      //    CellID cellID = 1 + round((x - Parameters::xmin) / dx + 
      //    (y - Parameters::ymin) / dy * Parameters::xcells_ini +
      //    (z - Parameters::zmin) / dz * Parameters::ycells_ini * Parameters::xcells_ini);
      
      //    return cellID * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5) *
      //    exp(- physicalconstants::MASS_PROTON * (vx*vx + vy*vy + vz*vz) / (2.0 * physicalconstants::K_B * this->T));
   }
   
   void Magnetosphere::calcCellParameters(Real* cellParams,creal& t) {
      set_dipole(cellParams);
   }

   Real Magnetosphere::getDistribValue(creal& x,creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      return this->rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5) *
      exp(- physicalconstants::MASS_PROTON * ((vx-this->V0[0])*(vx-this->V0[0]) + (vy-this->V0[1])*(vy-this->V0[1]) + (vz-this->V0[2])*(vz-this->V0[2])) / (2.0 * physicalconstants::K_B * this->T));
   }
   
} // namespace projects

