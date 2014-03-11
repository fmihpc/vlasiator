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

#include "Flowthrough.h"

using namespace std;

namespace projects {
   Flowthrough::Flowthrough(): TriAxisSearch() { }
   Flowthrough::~Flowthrough() { }
   
   bool Flowthrough::initialize(void) {return true;}
   
   void Flowthrough::addParameters(){
      typedef Readparameters RP;
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
   }
   
   void Flowthrough::getParameters(){
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      typedef Readparameters RP;
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
   }

   Real Flowthrough::getDistribValue(creal& x,creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      return this->rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5) *
      exp(- physicalconstants::MASS_PROTON * ((vx-this->V0[0])*(vx-this->V0[0]) + (vy-this->V0[1])*(vy-this->V0[1]) + (vz-this->V0[2])*(vz-this->V0[2])) / (2.0 * physicalconstants::K_B * this->T));
   }

   Real Flowthrough::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
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
                     return avg / (this->nSpaceSamples*this->nSpaceSamples*this->nSpaceSamples*this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
      
   //    CellID cellID = 1 + round((x - Parameters::xmin) / dx + 
   //    (y - Parameters::ymin) / dy * Parameters::xcells_ini +
   //    (z - Parameters::zmin) / dz * Parameters::ycells_ini * Parameters::xcells_ini);
      
   //    return cellID * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5) *
   //    exp(- physicalconstants::MASS_PROTON * (vx*vx + vy*vy + vz*vz) / (2.0 * physicalconstants::K_B * this->T));
   }

   void Flowthrough::calcCellParameters(Real* cellParams,creal& t) {
      cellParams[CellParams::PERBX] = this->Bx;
      cellParams[CellParams::PERBY] = this->By;
      cellParams[CellParams::PERBZ] = this->Bz;
   }
   
   vector<std::array<Real, 3>> Flowthrough::getV0(
      creal x,
      creal y,
      creal z
   ) {
      vector<std::array<Real, 3>> centerPoints;
      std::array<Real, 3> point {{this->V0[0], this->V0[1], this->V0[2]}};
      centerPoints.push_back(point);
      return centerPoints;
   }
   
} //namespace projects
