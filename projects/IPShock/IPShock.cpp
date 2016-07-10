/*
This file is part of Vlasiator.

Copyright 2016 Finnish Meteorological Institute
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

#include "IPShock.h"
#include "noise.h"


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
    RP::add("IPShock.BX0u", "Upstream mag. field value (T)", 1.0e-9);
    RP::add("IPShock.BY0u", "Upstream mag. field value (T)", 2.0e-9);
    RP::add("IPShock.BZ0u", "Upstream mag. field value (T)", 3.0e-9);
    RP::add("IPShock.VX0u", "Upstream Bulk velocity in x", 0.0);
    RP::add("IPShock.VY0u", "Upstream Bulk velocity in y", 0.0);
    RP::add("IPShock.VZ0u", "Upstream Bulk velocuty in z", 0.0);
    RP::add("IPShock.rhou", "Upstream Number density (m^-3)", 1.0e7);
    RP::add("IPShock.Temperatureu", "Upstream Temperature (K)", 2.0e6);

    RP::add("IPShock.BX0d", "Downstream mag. field value (T)", 1.0e-9);
    RP::add("IPShock.BY0d", "Downstream mag. field value (T)", 2.0e-9);
    RP::add("IPShock.BZ0d", "Downstream mag. field value (T)", 3.0e-9);
    RP::add("IPShock.VX0d", "Downstream Bulk velocity in x", 0.0);
    RP::add("IPShock.VY0d", "Downstream Bulk velocity in y", 0.0);
    RP::add("IPShock.VZ0d", "Downstream Bulk velocuty in z", 0.0);
    RP::add("IPShock.rhod", "Downstream Number density (m^-3)", 1.0e7);
    RP::add("IPShock.Temperatured", "Downstream Temperature (K)", 2.0e6);

    RP::add("IPShock.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
    RP::add("IPShock.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
    RP::add("IPShock.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
    RP::add("IPShock.Width", "Shock Width (m)", 50000);

    RP::add("IPShock.BPertAmp", "Amplitude of magnetic perturbation", 0.0);
    RP::add("IPShock.BPertScale", "Spatial scale of magnetic perturbation", 1.0);
    RP::add("IPShock.BPertOctaves", "Octaves of magnetic perturbation", 1);
  }

  void IPShock::getParameters() {
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
    if(!RP::get("IPShock.VX0u", this->V0u[0])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.VY0u", this->V0u[1])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.VZ0u", this->V0u[2])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.rhou", this->DENSITYu)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.Temperatureu", this->TEMPERATUREu)) {
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
    if(!RP::get("IPShock.VX0d", this->V0d[0])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.VY0d", this->V0d[1])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.VZ0d", this->V0d[2])) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.rhod", this->DENSITYd)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.Temperatured", this->TEMPERATUREd)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }

    if(!RP::get("IPShock.nSpaceSamples", this->nSpaceSamples)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.nVelocitySamples", this->nVelocitySamples)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.maxwCutoff", this->maxwCutoff)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.Width", this->Shockwidth)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }

    if(!RP::get("IPShock.BPertAmp", this->BPerturbationAmp)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.BPertScale", this->BPerturbationScale)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    if(!RP::get("IPShock.BPertOctaves", this->BPerturbationOctaves)) {
      if(myRank == MASTER_RANK) std::cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << std::endl;
      exit(1);
    }
    
    if(myRank == MASTER_RANK) {
      std::cerr << "B0x u = " << this->B0u[0] << std::endl;
      std::cerr << "B0y u = " << this->B0u[1] << std::endl;
      std::cerr << "B0z u = " << this->B0u[2] << std::endl;
      std::cerr << "B0x d = " << this->B0d[0] << std::endl;
      std::cerr << "B0y d = " << this->B0d[1] << std::endl;
      std::cerr << "B0z d = " << this->B0d[2] << std::endl;
      std::cerr << "V0x u = " << this->V0u[0] << std::endl;
      std::cerr << "V0y u = " << this->V0u[1] << std::endl;
      std::cerr << "V0z u = " << this->V0u[2] << std::endl;
      std::cerr << "V0x d = " << this->V0d[0] << std::endl;
      std::cerr << "V0y d = " << this->V0d[1] << std::endl;
      std::cerr << "V0z d = " << this->V0d[2] << std::endl;

      std::cerr << "rhou = " << this->DENSITYu << std::endl;
      std::cerr << "rhod = " << this->DENSITYd << std::endl;
      std::cerr << "tempu = " << this->TEMPERATUREu << std::endl;
      std::cerr << "tempd = " << this->TEMPERATUREd << std::endl;

      std::cerr << "nSpaceSamples = " << this->nSpaceSamples << std::endl;
      std::cerr << "nVelocitySamples = " << this->nVelocitySamples << std::endl;
      std::cerr << "maxwCutoff = " << this->maxwCutoff << std::endl;
      std::cerr << "Width = " << this->Shockwidth << std::endl;

      std::cerr << "BPertAmp = " << this->BPerturbationAmp << std::endl;
      std::cerr << "BPertScale = " << this->BPerturbationScale << std::endl;
      std::cerr << "BPertOctaves = " << this->BPerturbationOctaves << std::endl;

      if ((abs(this->V0u[1]) > 1e-10)||(abs(this->B0u[1]) > 1e-10)||(abs(this->V0d[1]) > 1e-10)||(abs(this->B0d[1]) > 1e-10))
	std::cout<<" Enforcing nil y-directional flow and magnetic flux "<<std::endl;
    }

    if ((abs(this->V0u[1]) > 1e-10)||(abs(this->B0u[1]) > 1e-10)||(abs(this->V0d[1]) > 1e-10)||(abs(this->B0d[1]) > 1e-10))
      {
	this->V0u[1] = 0.0;	this->B0u[1] = 0.0;	this->V0d[1] = 0.0;	this->B0d[1] = 0.0;
      }
    
  }


  std::vector<std::array<Real, 3>> IPShock::getV0(creal x, creal y, creal z) const {
    Real mass = physicalconstants::MASS_PROTON;
    Real mu0 = physicalconstants::MU_0;

    // Interpolate density between upstream and downstream
    // All other values are calculated from jump conditions
    Real DENSITY = interpolate(this->DENSITYu,this->DENSITYd, x);
    if (DENSITY < 1e-20) {
      std::cout<<"density too low! "<<DENSITY<<" x "<<x<<" y "<<y<<" z "<<z<<std::endl;
    }
    
    // Assume all velocities and magnetic fields are in x-z-plane
    // Alternatively, could set it so a rotation is performed for solving calculations and then returned to initial frame
    Real VX = this->DENSITYu * this->V0u[0] / DENSITY;
    Real BX = this->B0u[0];
    Real MAsq = std::pow((this->V0u[0]/this->B0u[0]), 2) * this->DENSITYu * mass * mu0;
    Real BZ = this->B0u[2] * (MAsq - 1.0)/(MAsq*VX/this->V0u[0] -1.0);
    Real VZ = VX * BZ / BX;

    // Disable compiler warnings: (unused variables but the function is inherited)
    (void)y; (void)z;
    
    std::array<Real, 3> V0 {{VX, 0.0, VZ}};
    std::vector<std::array<Real, 3>> retval;
    retval.push_back(V0);

    return retval;
  }

  Real IPShock::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) const {

    Real mass = physicalconstants::MASS_PROTON;
    Real KB = physicalconstants::K_B;
    Real mu0 = physicalconstants::MU_0;
    Real adiab = 5./3.;

    // Interpolate density between upstream and downstream
    // All other values are calculated from jump conditions
    Real DENSITY = interpolate(this->DENSITYu,this->DENSITYd, x);
    if (DENSITY < 1e-20) {
      std::cout<<"density too low! "<<DENSITY<<" x "<<x<<" y "<<y<<" z "<<z<<std::endl;
    }
    
    // Assume all velocities and magnetic fields are in x-z-plane
    // Alternatively, could set it so a rotation is performed for solving calculations and then returned to initial frame

    Real hereVX = this->DENSITYu * this->V0u[0] / DENSITY;
    Real hereBX = this->B0u[0];
    Real MAsq = std::pow((this->V0u[0]/this->B0u[0]), 2) * this->DENSITYu * mass * mu0;
    Real hereBZ = this->B0u[2] * (MAsq - 1.0)/(MAsq*hereVX/this->V0u[0] -1.0);
    Real hereVZ = hereVX * hereBZ / hereBX;
    // Old incorrect temperature - just interpolate for now
    //Real TEMPERATURE = this->TEMPERATUREu + (mass*(adiab-1.0)/(2.0*KB*adiab)) * 
    //  ( std::pow(this->V0u[0],2) + std::pow(this->V0u[2],2) - std::pow(hereVX,2) - std::pow(hereVZ,2) );
    Real TEMPERATURE = interpolate(this->TEMPERATUREu,this->TEMPERATUREd, x);

    std::array<Real, 3> pertV0 {{hereVX, 0.0, hereVZ}};

    Real result = 0.0;

    result = DENSITY * std::pow(mass / (2.0 * M_PI * KB * TEMPERATURE), 1.5) *
      exp(- mass * ((vx-pertV0[0])*(vx-pertV0[0]) + (vy-pertV0[1])*(vy-pertV0[1]) + (vz-pertV0[2])*(vz-pertV0[2])) / (2.0 * KB * TEMPERATURE));

    return result;
  }

  Real IPShock::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz, const int& popID) const {
    Real result = 0.0;
    if((this->nSpaceSamples > 1) && (this->nVelocitySamples > 1)) {
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
      
      result = avg /
	(this->nSpaceSamples*this->nSpaceSamples*this->nSpaceSamples) / 
	(this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
    } else {
      result = getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, dvx, dvy, dvz);
    }               

    if(result < this->maxwCutoff) {
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

    Vec3d position = Vec3d(x,y,z);
    Vec3d noise_scale = Vec3d(BPerturbationScale);

    Vec3d Bpert = BPerturbationAmp * divergence_free_noise(position, noise_scale, BPerturbationOctaves);

    Real mass = physicalconstants::MASS_PROTON;
    Real KB = physicalconstants::K_B;
    Real mu0 = physicalconstants::MU_0;
    Real adiab = 5./3.;

    // Interpolate density between upstream and downstream
    // All other values are calculated from jump conditions
    Real DENSITY = interpolate(this->DENSITYu,this->DENSITYd, x);
    
    // Assume all velocities and magnetic fields are in x-z-plane
    // Alternatively, could set it so a rotation is performed for solving calculations and then returned to initial frame
    Real VX = this->DENSITYu * this->V0u[0] / DENSITY;
    Real BX = this->B0u[0];
    Real MAsq = std::pow((this->V0u[0]/this->B0u[0]), 2) * this->DENSITYu * mass * mu0;
    Real BZ = this->B0u[2] * (MAsq - 1.0)/(MAsq*VX/this->V0u[0] -1.0);

    cellParams[CellParams::PERBX   ] = BX + Bpert[0];
    cellParams[CellParams::PERBY   ] = 0.0+ Bpert[1];
    cellParams[CellParams::PERBZ   ] = BZ + Bpert[2];

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
      //a = .5 * (1. + tanh(x * 2. * M_PI));
    }
    return a;
  }

  void IPShock::setCellBackgroundField(spatial_cell::SpatialCell* cell) const {
    setBackgroundFieldToZero(cell->parameters, cell->derivatives,cell->derivativesBVOL);
  }

}//namespace projects
