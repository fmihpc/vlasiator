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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

%
Interplanetary shock project by Markus Battarbee (markus.battarbee@gmail.com)
Based on SilvaShock project by Urs Ganse
Previous development version name was UtuShock
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <vector>
#include "vectorclass.h"
#include "vector3d.h"

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../object_wrapper.h"

#include "IPShock.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
  IPShock::IPShock(): TriAxisSearch() { }
  IPShock::~IPShock() { }
  
  bool IPShock::initialize() {
    return Project::initialize();
  }
  
  void IPShock::addParameters() {
    typedef Readparameters RP;
    // Common (field / etc.) parameters
    RP::add("IPShock.BX0u", "Upstream mag. field value (T)", 1.0e-9);
    RP::add("IPShock.BY0u", "Upstream mag. field value (T)", 2.0e-9);
    RP::add("IPShock.BZ0u", "Upstream mag. field value (T)", 3.0e-9);
    RP::add("IPShock.BX0d", "Downstream mag. field value (T)", 1.0e-9);
    RP::add("IPShock.BY0d", "Downstream mag. field value (T)", 2.0e-9);
    RP::add("IPShock.BZ0d", "Downstream mag. field value (T)", 3.0e-9);
    RP::add("IPShock.Width", "Shock Width (m)", 50000);

    // Per-population parameters
    for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
       const std::string& pop = getObjectWrapper().particleSpecies[i].name;
       RP::add(pop + "_IPShock.VX0u", "Upstream Bulk velocity in x", 0.0);
       RP::add(pop + "_IPShock.VY0u", "Upstream Bulk velocity in y", 0.0);
       RP::add(pop + "_IPShock.VZ0u", "Upstream Bulk velocuty in z", 0.0);
       RP::add(pop + "_IPShock.rhou", "Upstream Number density (m^-3)", 1.0e7);
       RP::add(pop + "_IPShock.Temperatureu", "Upstream Temperature (K)", 2.0e6);

       RP::add(pop + "_IPShock.VX0d", "Downstream Bulk velocity in x", 0.0);
       RP::add(pop + "_IPShock.VY0d", "Downstream Bulk velocity in y", 0.0);
       RP::add(pop + "_IPShock.VZ0d", "Downstream Bulk velocuty in z", 0.0);
       RP::add(pop + "_IPShock.rhod", "Downstream Number density (m^-3)", 1.0e7);
       RP::add(pop + "_IPShock.Temperatured", "Downstream Temperature (K)", 2.0e6);

       RP::add(pop + "_IPShock.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
       RP::add(pop + "_IPShock.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
       RP::add(pop + "_IPShock.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
    }

  }

  void IPShock::getParameters() {
    Project::getParameters();
    int myRank;

    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    typedef Readparameters RP;
    if(!RP::get("IPShock.BX0u", this->B0u[0])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.BY0u", this->B0u[1])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.BZ0u", this->B0u[2])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.BX0d", this->B0d[0])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.BY0d", this->B0d[1])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.BZ0d", this->B0d[2])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.Width", this->Shockwidth)) {
       if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
       exit(1);
    }

    // Per-population parameters
    for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
       const std::string& pop = getObjectWrapper().particleSpecies[i].name;
       IPShockSpeciesParameters sP;

       if(!RP::get(pop + "_IPShock.VX0u", sP.V0u[0])) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.VY0u", sP.V0u[1])) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.VZ0u", sP.V0u[2])) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.rhou", sP.DENSITYu)) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.Temperatureu", sP.TEMPERATUREu)) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }

       if(!RP::get(pop + "_IPShock.VX0d", sP.V0d[0])) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.VY0d", sP.V0d[1])) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.VZ0d", sP.V0d[2])) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.rhod", sP.DENSITYd)) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.Temperatured", sP.TEMPERATUREd)) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }

       if(!RP::get(pop + "_IPShock.nSpaceSamples", sP.nSpaceSamples)) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.nVelocitySamples", sP.nVelocitySamples)) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }
       if(!RP::get(pop + "_IPShock.maxwCutoff", sP.maxwCutoff)) {
          if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
          exit(1);
       }

       speciesParams.push_back(sP);
    }
    
    if(myRank == MASTER_RANK) {
      std::cerr << "B0x u = " << this->B0u[0] << std::endl;
      std::cerr << "B0y u = " << this->B0u[1] << std::endl;
      std::cerr << "B0z u = " << this->B0u[2] << std::endl;
      std::cerr << "B0x d = " << this->B0d[0] << std::endl;
      std::cerr << "B0y d = " << this->B0d[1] << std::endl;
      std::cerr << "B0z d = " << this->B0d[2] << std::endl;
      //std::cerr << "V0x u = " << this->V0u[0] << std::endl;
      //std::cerr << "V0y u = " << this->V0u[1] << std::endl;
      //std::cerr << "V0z u = " << this->V0u[2] << std::endl;
      //std::cerr << "V0x d = " << this->V0d[0] << std::endl;
      //std::cerr << "V0y d = " << this->V0d[1] << std::endl;
      //std::cerr << "V0z d = " << this->V0d[2] << std::endl;

      //std::cerr << "rhou = " << this->DENSITYu << std::endl;
      //std::cerr << "rhod = " << this->DENSITYd << std::endl;
      //std::cerr << "tempu = " << this->TEMPERATUREu << std::endl;
      //std::cerr << "tempd = " << this->TEMPERATUREd << std::endl;

      //std::cerr << "nSpaceSamples = " << this->nSpaceSamples << std::endl;
      //std::cerr << "nVelocitySamples = " << this->nVelocitySamples << std::endl;
      //std::cerr << "maxwCutoff = " << this->maxwCutoff << std::endl;
      //std::cerr << "Width = " << this->Shockwidth << std::endl;
    }

    /* 
       Now allows flow and field both in z and y -directions. As assuming we're
       in the dHT frame, all flow and magnetic field should be in a single plane.
    */

    /* Magnitude of tangential B-field and flow components */
    this->B0utangential = sqrt(this->B0u[1]*this->B0u[1] + this->B0u[2]*this->B0u[2]);
    this->B0dtangential = sqrt(this->B0d[1]*this->B0d[1] + this->B0d[2]*this->B0d[2]);
    for(auto& sP : speciesParams) {
       sP.V0utangential = sqrt(sP.V0u[1]*sP.V0u[1] + sP.V0u[2]*sP.V0u[2]);
       sP.V0dtangential = sqrt(sP.V0d[1]*sP.V0d[1] + sP.V0d[2]*sP.V0d[2]);
    }

    /* Check direction of upstream and downstream flows and fields 
       Define y-z-directional angle phi so that 
       By = cos(phi_B)*B_tang
       Bz = sin(phi_B)*B_tang
       Vy = cos(phi_V)*V_tang
       Vz = sin(phi_V)*V_tang
       If we're in the dHT frame, phi_B and phi_V should be the same, and also the same
       both in the upstream and in the downstream.
    */
    this->Bucosphi = abs(this->B0u[1])/this->B0utangential;   
    this->Bdcosphi = abs(this->B0d[1])/this->B0dtangential;   
    for(auto& sP : speciesParams) {
       sP.Vucosphi = abs(sP.V0u[1])/sP.V0utangential;
       sP.Vdcosphi = abs(sP.V0d[1])/sP.V0dtangential;
    }

    /* Save signs as well for reconstruction during interpolation.
       For both components of B and V, upstream and downstream signs should be the same. */
    this->Byusign=0;
    if (this->B0u[1] < 0) this->Byusign=-1;
    if (this->B0u[1] > 0) this->Byusign=+1;
    this->Bzusign=0;
    if (this->B0u[2] < 0) this->Bzusign=-1;
    if (this->B0u[2] > 0) this->Bzusign=+1;
    this->Bydsign=0;
    if (this->B0d[1] < 0) this->Bydsign=-1;
    if (this->B0d[1] > 0) this->Bydsign=+1;
    this->Bzdsign=0;
    if (this->B0d[2] < 0) this->Bzdsign=-1;
    if (this->B0d[2] > 0) this->Bzdsign=+1;

    for(auto& sP : speciesParams) {
       sP.Vyusign=0;
       if (sP.V0u[1] < 0) sP.Vyusign=-1;
       if (sP.V0u[1] > 0) sP.Vyusign=+1;
       sP.Vzusign=0;
       if (sP.V0u[2] < 0) sP.Vzusign=-1;
       if (sP.V0u[2] > 0) sP.Vzusign=+1;
       sP.Vydsign=0;
       if (sP.V0d[1] < 0) sP.Vydsign=-1;
       if (sP.V0d[1] > 0) sP.Vydsign=+1;
       sP.Vzdsign=0;
       if (sP.V0d[2] < 0) sP.Vzdsign=-1;
       if (sP.V0d[2] > 0) sP.Vzdsign=+1;

       /* Check that upstream and downstream values both are separately parallel */
       if ( (abs(this->Bucosphi)-abs(sP.Vucosphi) > 1e-10) || (this->Byusign*this->Bzusign != sP.Vyusign*sP.Vzusign) )
       {
          if(myRank == MASTER_RANK) {
             std::cout<<" Warning: Upstream B and V not parallel"<<std::endl;
             std::cout<<" Bucosphi "<<Bucosphi<<" Vucosphi "<<sP.Vucosphi<<" Byusign "<<Byusign<<" Bzusign "<<Bzusign<<" Vyusign "<<sP.Vyusign<<" Vzusign "<<sP.Vzusign<<std::endl;
          }
       }
       if ( (abs(this->Bdcosphi)-abs(sP.Vdcosphi) > 1e-10) || (this->Bydsign*this->Bzdsign != sP.Vydsign*sP.Vzdsign) )
       {
          if(myRank == MASTER_RANK) {
             std::cout<<" Warning: Downstream B and V not parallel"<<std::endl;
             std::cout<<" Bdcosphi "<<Bdcosphi<<" Vdcosphi "<<sP.Vdcosphi<<" Bydsign "<<Bydsign<<" Bzdsign "<<Bzdsign<<" Vydsign "<<sP.Vydsign<<" Vzdsign "<<sP.Vzdsign<<std::endl;
          }
       }
       /* Verify that upstream and downstream flows are in a plane */
       if ( (abs(this->Bdcosphi)-abs(this->Bucosphi) > 1e-10) && (this->Bydsign*this->Bzdsign != this->Byusign*this->Bzusign) )
       {
          if(myRank == MASTER_RANK) {
             std::cout<<" Warning: Upstream and downstream B_tangentials not in same plane"<<std::endl;
             std::cout<<" Bdcosphi "<<Bdcosphi<<" Bucosphi "<<Bucosphi<<" Bydsign "<<Bydsign<<" Bzdsign "<<Bzdsign<<" Byusign "<<Byusign<<" Bzusign "<<Bzusign<<std::endl;
          }
       }
       if ( (abs(sP.Vdcosphi)-abs(sP.Vucosphi) > 1e-10) && (sP.Vydsign*sP.Vzdsign != sP.Vyusign*sP.Vzusign) )
       {
          if(myRank == MASTER_RANK) {
             std::cout<<" Warning: Upstream and downstream V_tangentials not in same plane"<<std::endl;
             std::cout<<" Vdcosphi "<<sP.Vdcosphi<<" Vucosphi "<<sP.Vucosphi<<" Vydsign "<<sP.Vydsign<<" Vzdsign "<<sP.Vzdsign<<" Vyusign "<<sP.Vyusign<<" Vzusign "<<sP.Vzusign<<std::endl;
          }
       }
    }
    
  }


  std::vector<std::array<Real, 3>> IPShock::getV0(creal x, creal y, creal z, const uint popID) const {
    Real mass = getObjectWrapper().particleSpecies[popID].mass;
    Real mu0 = physicalconstants::MU_0;
    const IPShockSpeciesParameters& sP = this->speciesParams[popID];

    // Interpolate density between upstream and downstream
    // All other values are calculated from jump conditions
    Real DENSITY = interpolate(sP.DENSITYu,sP.DENSITYd, x);
    if (DENSITY < 1e-20) {
      std::cout<<"density too low! "<<DENSITY<<" x "<<x<<" y "<<y<<" z "<<z<<std::endl;
    }
    
    // Solve tangential components for B and V
    Real VX = sP.DENSITYu * sP.V0u[0] / DENSITY;
    Real BX = this->B0u[0];
    Real MAsq = std::pow((sP.V0u[0]/this->B0u[0]), 2) * sP.DENSITYu * mass * mu0;
    Real Btang = this->B0utangential * (MAsq - 1.0)/(MAsq*VX/sP.V0u[0] -1.0);
    Real Vtang = VX * Btang / BX;

    /* Reconstruct Y and Z components using cos(phi) values and signs. Tangential variables are always positive. */
    //Real BY = Btang * this->Bucosphi * this->Byusign;
    //Real BZ = Btang * sqrt(1. - this->Bucosphi * this->Bucosphi) * this->Bzusign;
    Real VY = abs(Vtang) * sP.Vucosphi * sP.Vyusign;
    Real VZ = abs(Vtang) * sqrt(1. - sP.Vucosphi * sP.Vucosphi) * sP.Vzusign;

    // Disable compiler warnings: (unused variables but the function is inherited)
    (void)y;
    (void)z;
    
    std::array<Real, 3> V0 {{VX, VY, VZ}};
    std::vector<std::array<Real, 3>> retval;
    retval.push_back(V0);

    return retval;
  }

  Real IPShock::getDistribValue(creal& x, creal& y, creal& z, 
				creal& vx, creal& vy, creal& vz, 
				creal& dvx, creal& dvy, creal& dvz,
        const uint popID) const {

    Real mass = getObjectWrapper().particleSpecies[popID].mass;
    Real KB = physicalconstants::K_B;
    Real mu0 = physicalconstants::MU_0;
    Real adiab = 5./3.;
    const IPShockSpeciesParameters& sP = this->speciesParams[popID];

    // Interpolate density between upstream and downstream
    // All other values are calculated from jump conditions
    Real DENSITY = interpolate(sP.DENSITYu,sP.DENSITYd, x);
    if (DENSITY < 1e-20) {
      std::cout<<"density too low! "<<DENSITY<<" x "<<x<<" y "<<y<<" z "<<z<<std::endl;
    }
    
    // Solve tangential components for B and V
    Real hereVX = sP.DENSITYu * sP.V0u[0] / DENSITY;
    Real hereBX = this->B0u[0];
    Real MAsq = std::pow((sP.V0u[0]/this->B0u[0]), 2) * sP.DENSITYu * mass * mu0;
    Real hereBtang = this->B0u[2] * (MAsq - 1.0)/(MAsq*hereVX/sP.V0u[0] -1.0);
    Real hereVtang = hereVX * hereBtang / hereBX;

    /* Reconstruct Y and Z components using cos(phi) values and signs. Tangential variables are always positive. */
    //Real hereBY = hereBtang * this->Bucosphi * this->Byusign;
    //Real hereBZ = hereBtang * sqrt(1. - this->Bucosphi * this->Bucosphi) * this->Bzusign;
    Real hereVY = abs(hereVtang) * sP.Vucosphi * sP.Vyusign;
    Real hereVZ = abs(hereVtang) * sqrt(1. - sP.Vucosphi * sP.Vucosphi) * sP.Vzusign;

    // Old incorrect temperature - just interpolate for now
    //Real TEMPERATURE = this->TEMPERATUREu + (mass*(adiab-1.0)/(2.0*KB*adiab)) * 
    //  ( std::pow(this->V0u[0],2) + std::pow(this->V0u[2],2) - std::pow(hereVX,2) - std::pow(hereVZ,2) );
    Real TEMPERATURE = interpolate(sP.TEMPERATUREu,sP.TEMPERATUREd, x);

    std::array<Real, 3> pertV0 {{hereVX, hereVY, hereVZ}};

    Real result = 0.0;

    result = DENSITY * std::pow(mass / (2.0 * M_PI * KB * TEMPERATURE), 1.5) *
      exp(- mass * ((vx-pertV0[0])*(vx-pertV0[0]) + (vy-pertV0[1])*(vy-pertV0[1]) + (vz-pertV0[2])*(vz-pertV0[2])) / (2.0 * KB * TEMPERATURE));

    return result;
  }

  Real IPShock::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz, const uint popID) const {

    const IPShockSpeciesParameters& sP = this->speciesParams[popID];
    Real result = 0.0;
    if((sP.nSpaceSamples > 1) && (sP.nVelocitySamples > 1)) {
      creal d_x = dx / (sP.nSpaceSamples-1);
      creal d_y = dy / (sP.nSpaceSamples-1);
      creal d_z = dz / (sP.nSpaceSamples-1);
      creal d_vx = dvx / (sP.nVelocitySamples-1);
      creal d_vy = dvy / (sP.nVelocitySamples-1);
      creal d_vz = dvz / (sP.nVelocitySamples-1);

      Real avg = 0.0;
      
      for (uint i=0; i<sP.nSpaceSamples; ++i)
         for (uint j=0; j<sP.nSpaceSamples; ++j)
            for (uint k=0; k<sP.nSpaceSamples; ++k)      
               for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
                  for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
                     for (uint vk=0; vk<sP.nVelocitySamples; ++vk)
                     {
                        avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz, popID);
                     }
      
      result = avg /
	(sP.nSpaceSamples*sP.nSpaceSamples*sP.nSpaceSamples) / 
	(sP.nVelocitySamples*sP.nVelocitySamples*sP.nVelocitySamples);
    } else {
      result = getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, dvx, dvy, dvz, popID);
    }               

    if(result < sP.maxwCutoff) {
      return 0.0;
    } else {
      return result;
    }
  }
  
  void IPShock::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {
    // Disable compiler warnings: (unused variables but the function is inherited)
    (void)t;

    /* Maintain all values in BPERT for simplicity */
    Real* cellParams = cell->get_cell_parameters();

    creal x = cellParams[CellParams::XCRD];
    creal dx = cellParams[CellParams::DX];
    creal y = cellParams[CellParams::YCRD];
    creal dy = cellParams[CellParams::DY];
    creal z = cellParams[CellParams::ZCRD];
    creal dz = cellParams[CellParams::DZ];

    cellParams[CellParams::EX   ] = 0.0;
    cellParams[CellParams::EY   ] = 0.0;
    cellParams[CellParams::EZ   ] = 0.0;

    Real KB = physicalconstants::K_B;
    Real mu0 = physicalconstants::MU_0;
    Real adiab = 5./3.;

    // Interpolate density between upstream and downstream
    // All other values are calculated from jump conditions
    Real MassDensity = 0.;
    Real MassDensityU = 0.;
    Real EffectiveVu0 = 0.;
    for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
       const IPShockSpeciesParameters& sP = speciesParams[i];
       Real mass = getObjectWrapper().particleSpecies[i].mass;

       MassDensity += mass * interpolate(sP.DENSITYu,sP.DENSITYd, x);
       MassDensityU += mass * sP.DENSITYu;
       EffectiveVu0 += sP.V0u[0] * mass * sP.DENSITYu;
    }
    EffectiveVu0 /= MassDensityU;

    // Solve tangential components for B and V
    Real VX = MassDensityU * EffectiveVu0 / MassDensity;
    Real BX = this->B0u[0];
    Real MAsq = std::pow((EffectiveVu0/this->B0u[0]), 2) * MassDensityU * mu0;
    Real Btang = this->B0utangential * (MAsq - 1.0)/(MAsq*VX/EffectiveVu0 -1.0);
    Real Vtang = VX * Btang / BX;

    /* Reconstruct Y and Z components using cos(phi) values and signs. Tangential variables are always positive. */
    Real BY = abs(Btang) * this->Bucosphi * this->Byusign;
    Real BZ = abs(Btang) * sqrt(1. - this->Bucosphi * this->Bucosphi) * this->Bzusign;
    //Real VY = Vtang * this->Vucosphi * this->Vyusign;
    //Real VZ = Vtang * sqrt(1. - this->Vucosphi * this->Vucosphi) * this->Vzusign;

    cellParams[CellParams::PERBX   ] = BX;
    cellParams[CellParams::PERBY   ] = BY;
    cellParams[CellParams::PERBZ   ] = BZ;

  }

  Real IPShock::interpolate(Real upstream, Real downstream, Real x) const {
    Real coord = 0.5 + x/this->Shockwidth; //Now shock will be from 0 to 1
    //x /= 0.5 * this->Shockwidth;
    Real a = 0.0;
    if (coord <= 0.0) a = downstream;
    if (coord >= 1.0) a = upstream;
    if ((coord > 0.0) && (coord < 1.0)) {
      // Ken Perlin Smootherstep
      Real interpolation = ( 6.0 * coord * coord - 15.0 * coord +10. ) * coord * coord * coord;
      a = upstream * interpolation + downstream * (1. - interpolation);
    }
    return a;
  }

  void IPShock::setCellBackgroundField(spatial_cell::SpatialCell* cell) const {
    setBackgroundFieldToZero(cell->parameters, cell->derivatives,cell->derivativesBVOL);
  }

}//namespace projects
