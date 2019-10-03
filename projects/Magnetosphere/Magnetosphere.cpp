/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
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
#include <array>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../backgroundfield/dipole.hpp"
#include "../../backgroundfield/linedipole.hpp"
#include "../../backgroundfield/vectordipole.hpp"
#include "../../object_wrapper.h"

#include "Magnetosphere.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   Magnetosphere::Magnetosphere(): TriAxisSearch() { }
   Magnetosphere::~Magnetosphere() { }
   
   void Magnetosphere::addParameters() {
      typedef Readparameters RP;
      // Common (field / etc.) parameters
      RP::add("Magnetosphere.constBgBX", "Constant flat Bx component in the whole simulation box. Default is none.", 0.0);
      RP::add("Magnetosphere.constBgBY", "Constant flat By component in the whole simulation box. Default is none.", 0.0);
      RP::add("Magnetosphere.constBgBZ", "Constant flat Bz component in the whole simulation box. Default is none.", 0.0);
      RP::add("Magnetosphere.noDipoleInSW", "If set to 1, the dipole magnetic field is not set in the solar wind inflow cells. Default 0.", 0.0);
      RP::add("Magnetosphere.dipoleScalingFactor","Scales the field strength of the magnetic dipole compared to Earths.", 1.0);
      RP::add("Magnetosphere.dipoleType","0: Normal 3D dipole, 1: line-dipole for 2D polar simulations, 2: line-dipole with mirror, 3: 3D dipole with mirror", 0);
      RP::add("Magnetosphere.dipoleMirrorLocationX","x-coordinate of dipole Mirror", -1.0);

      RP::add("Magnetosphere.refine_L3radius","Radius of L3-refined sphere", 6.371e7); // 10 RE
      RP::add("Magnetosphere.refine_L3nosexmin","Low x-value of nose L3-refined box", 5.0e7); //
      RP::add("Magnetosphere.refine_L3tailheight","Height in +-z of tail L3-refined box", 1.0e7); //
      RP::add("Magnetosphere.refine_L3tailwidth","Width in +-y of tail L3-refined box", 5.0e7); // 10 RE
      RP::add("Magnetosphere.refine_L3tailxmin","Low x-value of tail L3-refined box", -20.0e7); // 10 RE
      RP::add("Magnetosphere.refine_L3tailxmax","High x-value of tail L3-refined box", -5.0e7); // 10 RE
      
      RP::add("Magnetosphere.refine_L2radius","Radius of L2-refined sphere", 9.5565e7); // 15 RE
      RP::add("Magnetosphere.refine_L2tailthick","Thickness of L2-refined tail region", 3.1855e7); // 5 RE
      RP::add("Magnetosphere.refine_L1radius","Radius of L1-refined sphere", 1.59275e8); // 25 RE
      RP::add("Magnetosphere.refine_L1tailthick","Thickness of L1-refined tail region", 6.371e7); // 10 RE

      RP::add("Magnetosphere.dipoleTiltPhi","Magnitude of dipole tilt, in degrees", 0.0);
      RP::add("Magnetosphere.dipoleTiltTheta","Direction of dipole tilt from Sun-Earth-line, in degrees", 0.0);
      RP::add("Magnetosphere.dipoleXFull","X-coordinate up to which dipole is at full strength, in metres", 9.5565e7); // 15 RE
      RP::add("Magnetosphere.dipoleXZero","X-coordinate after which dipole is at zero strength, in metres", 1.9113e8); // 30 RE
      RP::add("Magnetosphere.dipoleInflowBX","Inflow magnetic field Bx component to which the vector potential dipole converges. Default is none.", 0.0);
      RP::add("Magnetosphere.dipoleInflowBY","Inflow magnetic field By component to which the vector potential dipole converges. Default is none.", 0.0);
      RP::add("Magnetosphere.dipoleInflowBZ","Inflow magnetic field Bz component to which the vector potential dipole converges. Default is none.", 0.0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_Magnetosphere.rho", "Tail region number density (m^-3)", 0.0);
         RP::add(pop + "_Magnetosphere.T", "Temperature (K)", 0.0);
         RP::add(pop + "_Magnetosphere.VX0", "Initial bulk velocity in x-direction", 0.0);
         RP::add(pop + "_Magnetosphere.VY0", "Initial bulk velocity in y-direction", 0.0);
         RP::add(pop + "_Magnetosphere.VZ0", "Initial bulk velocity in z-direction", 0.0);
         RP::add(pop + "_Magnetosphere.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
         RP::add(pop + "_Magnetosphere.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
      }
   }
   
   void Magnetosphere::getParameters(){
      Project::getParameters();
      
      int myRank;
      Real dummy;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      typedef Readparameters RP;
      if(!RP::get("Magnetosphere.constBgBX", this->constBgB[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.constBgBY", this->constBgB[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.constBgBZ", this->constBgB[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.noDipoleInSW", dummy)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      this->noDipoleInSW = dummy == 1 ? true:false;
      if(!RP::get("Magnetosphere.dipoleScalingFactor", this->dipoleScalingFactor)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }

      if(!RP::get("Magnetosphere.dipoleMirrorLocationX", this->dipoleMirrorLocationX)) {
           if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }

      if(!RP::get("Magnetosphere.dipoleType", this->dipoleType)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);       
      }
      if(!RP::get("ionosphere.radius", this->ionosphereRadius)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.centerX", this->center[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.centerY", this->center[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.centerZ", this->center[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.geometry", this->ionosphereGeometry)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }


      if(!Readparameters::get("Magnetosphere.refine_L3radius", this->refine_L3radius)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.refine_L3nosexmin", this->refine_L3nosexmin)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.refine_L3tailwidth", this->refine_L3tailwidth)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.refine_L3tailheight", this->refine_L3tailheight)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.refine_L3tailxmin", this->refine_L3tailxmin)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.refine_L3tailxmax", this->refine_L3tailxmax)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }

      if(!Readparameters::get("Magnetosphere.refine_L2radius", this->refine_L2radius)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.refine_L2tailthick", this->refine_L2tailthick)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.refine_L1radius", this->refine_L1radius)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.refine_L1tailthick", this->refine_L1tailthick)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }

      if(!Readparameters::get("Magnetosphere.dipoleTiltPhi", this->dipoleTiltPhi)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.dipoleTiltTheta", this->dipoleTiltTheta)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.dipoleXFull", this->dipoleXFull)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.dipoleXZero", this->dipoleXZero)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.dipoleInflowBX", this->dipoleInflowB[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.dipoleInflowBY", this->dipoleInflowB[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.dipoleInflowBZ", this->dipoleInflowB[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         MagnetosphereSpeciesParameters sP;

         if(!RP::get(pop + "_Magnetosphere.rho", sP.rho)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_Magnetosphere.T", sP.T)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_Magnetosphere.VX0", sP.V0[0])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_Magnetosphere.VY0", sP.V0[1])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_Magnetosphere.VZ0", sP.V0[2])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }

         if(!RP::get(pop + "_Magnetosphere.nSpaceSamples", sP.nSpaceSamples)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_Magnetosphere.nVelocitySamples", sP.nVelocitySamples)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }

         if(!RP::get(pop + "_ionosphere.rho", sP.ionosphereRho)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_ionosphere.VX0", sP.ionosphereV0[0])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_ionosphere.VY0", sP.ionosphereV0[1])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_ionosphere.VZ0", sP.ionosphereV0[2])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!RP::get(pop + "_ionosphere.taperRadius", sP.ionosphereTaperRadius)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
            exit(1);
         }

        speciesParams.push_back(sP);
      }

   }
   
   bool Magnetosphere::initialize() {
      return Project::initialize();
   }

   Real Magnetosphere::calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
                                             creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,
                                             creal& dvz,const uint popID) const {

      const MagnetosphereSpeciesParameters& sP = this->speciesParams[popID];

      if((sP.nSpaceSamples > 1) && (sP.nVelocitySamples > 1)) {
         creal d_x = dx / (sP.nSpaceSamples-1);
         creal d_y = dy / (sP.nSpaceSamples-1);
         creal d_z = dz / (sP.nSpaceSamples-1);
         creal d_vx = dvx / (sP.nVelocitySamples-1);
         creal d_vy = dvy / (sP.nVelocitySamples-1);
         creal d_vz = dvz / (sP.nVelocitySamples-1);
         
         Real avg = 0.0;
         // #pragma omp parallel for collapse(6) reduction(+:avg)
         // WARNING No threading here if calling functions are already threaded
         for (uint i=0; i<sP.nSpaceSamples; ++i)
            for (uint j=0; j<sP.nSpaceSamples; ++j)
               for (uint k=0; k<sP.nSpaceSamples; ++k)
                  for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
                     for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
                        for (uint vk=0; vk<sP.nVelocitySamples; ++vk) {
                           avg += getDistribValue(x+i*d_x,y+j*d_y,z+k*d_z,vx+vi*d_vx,vy+vj*d_vy,vz+vk*d_vz,dvx,dvy,dvz,popID);
                        }
         return avg /
         (sP.nSpaceSamples*sP.nSpaceSamples*sP.nSpaceSamples) /
         (sP.nVelocitySamples*sP.nVelocitySamples*sP.nVelocitySamples);
      } else {
         return getDistribValue(x+0.5*dx,y+0.5*dy,z+0.5*dz,vx+0.5*dvx,vy+0.5*dvy,vz+0.5*dvz,dvx,dvy,dvz,popID);
      }
   }
   
   /*! Magnetosphere does not set any extra perturbed B. */
   void Magnetosphere::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   /* set 0-centered dipole */
   void Magnetosphere::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
      FsGrid< fsgrids::technical, 2>& technicalGrid
   ) {
      Dipole bgFieldDipole;
      LineDipole bgFieldLineDipole;
      VectorDipole bgVectorDipole;

      // The hardcoded constants of dipole and line dipole moments are obtained
      // from Daldorff et al (2014), see
      // https://github.com/fmihpc/vlasiator/issues/20 for a derivation of the
      // values used here.
      switch(this->dipoleType) {
            case 0:
               bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
               setBackgroundField(bgFieldDipole, BgBGrid);
               break;
            case 1:
               bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
               setBackgroundField(bgFieldLineDipole, BgBGrid);
               break;
            case 2:
               bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
               setBackgroundField(bgFieldLineDipole, BgBGrid);
               //Append mirror dipole
               bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, this->dipoleMirrorLocationX, 0.0, 0.0 );
               setBackgroundField(bgFieldLineDipole, BgBGrid, true);
               break;
            case 3:
               bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
               setBackgroundField(bgFieldDipole, BgBGrid);
               //Append mirror dipole                
               bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, this->dipoleMirrorLocationX, 0.0, 0.0, 0.0 );//mirror
               setBackgroundField(bgFieldDipole, BgBGrid, true);
               break; 
            case 4:  // Vector potential dipole, vanishes or optionally scales to static inflow value after a given x-coordinate
	       // What we in fact do is we place the regular dipole in the background field, and the
	       // corrective terms in the perturbed field. This maintains the BGB as curl-free.
	       bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
               setBackgroundField(bgFieldDipole, BgBGrid);
	       // Difference into perBgrid, only if not restarting
	       if (P::isRestart == false) {
		  bgFieldDipole.initialize(-8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );
		  setPerturbedField(bgFieldDipole, perBGrid);
		  bgVectorDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, this->dipoleTiltPhi*3.14159/180., this->dipoleTiltTheta*3.14159/180., this->dipoleXFull, this->dipoleXZero, this->dipoleInflowB[0], this->dipoleInflowB[1], this->dipoleInflowB[2]);
		  setPerturbedField(bgVectorDipole, perBGrid, true);
	       }
               break;              
            default:
               setBackgroundFieldToZero(BgBGrid);
      }
      
      const auto localSize = BgBGrid.getLocalSize();
      
#pragma omp parallel
      {
         //Force field to zero in the perpendicular direction for 2D (1D) simulations. Otherwise we have unphysical components.
         if(P::xcells_ini==1) {
#pragma omp for collapse(3)
            for (int x = 0; x < localSize[0]; ++x) {
               for (int y = 0; y < localSize[1]; ++y) {
                  for (int z = 0; z < localSize[2]; ++z) {
                     std::array<Real, fsgrids::bgbfield::N_BGB>* cell = BgBGrid.get(x, y, z);
                     cell->at(fsgrids::bgbfield::BGBX)=0;
                     cell->at(fsgrids::bgbfield::BGBXVOL)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBzdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBZVOLdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdz)=0.0;
                  }
               }
            }
         }
         if(P::ycells_ini==1) {
            /*2D simulation in x and z. Set By and derivatives along Y, and derivatives of By to zero*/
#pragma omp for collapse(3)
            for (int x = 0; x < localSize[0]; ++x) {
               for (int y = 0; y < localSize[1]; ++y) {
                  for (int z = 0; z < localSize[2]; ++z) {
                     std::array<Real, fsgrids::bgbfield::N_BGB>* cell = BgBGrid.get(x, y, z);
                     cell->at(fsgrids::bgbfield::BGBY)=0.0;
                     cell->at(fsgrids::bgbfield::BGBYVOL)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBzdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBZVOLdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdz)=0.0;
                  }
               }
            }
         }
         if(P::zcells_ini==1) {
#pragma omp for collapse(3)
            for (int x = 0; x < localSize[0]; ++x) {
               for (int y = 0; y < localSize[1]; ++y) {
                  for (int z = 0; z < localSize[2]; ++z) {
                     std::array<Real, fsgrids::bgbfield::N_BGB>* cell = BgBGrid.get(x, y, z);
                     cell->at(fsgrids::bgbfield::BGBX)=0;
                     cell->at(fsgrids::bgbfield::BGBY)=0;
                     cell->at(fsgrids::bgbfield::BGBYVOL)=0.0;
                     cell->at(fsgrids::bgbfield::BGBXVOL)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBxdz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBydz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdy)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBXVOLdz)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdx)=0.0;
                     cell->at(fsgrids::bgbfield::dBGBYVOLdz)=0.0;
                  }
               }
            }
         }
         
         // Remove dipole from inflow cells if this is requested
         if(this->noDipoleInSW) {
#pragma omp for collapse(3)
            for (int x = 0; x < localSize[0]; ++x) {
               for (int y = 0; y < localSize[1]; ++y) {
                  for (int z = 0; z < localSize[2]; ++z) {
                     if(technicalGrid.get(x, y, z)->sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN ) {
                        for (int i = 0; i < fsgrids::bgbfield::N_BGB; ++i) {
                           BgBGrid.get(x,y,z)->at(i) = 0;
                        }
                     }
                  }
               }
            }
         }
      } // end of omp parallel region
      // Superimpose constant background field if needed
      if(this->constBgB[0] != 0.0 || this->constBgB[1] != 0.0 || this->constBgB[2] != 0.0) {
         ConstantField bgConstantField;
         bgConstantField.initialize(this->constBgB[0], this->constBgB[1], this->constBgB[2]);
         setBackgroundField(bgConstantField, BgBGrid, true);
      }
   }
   
   
   Real Magnetosphere::getDistribValue(
           creal& x,creal& y,creal& z,
           creal& vx,creal& vy,creal& vz,
           creal& dvx,creal& dvy,creal& dvz,
           const uint popID) const
   {
      const MagnetosphereSpeciesParameters& sP = this->speciesParams[popID];
      Real initRho = sP.rho;
      std::array<Real, 3> initV0 = this->getV0(x, y, z, popID)[0];
      
      Real radius;
      
      switch(this->ionosphereGeometry) {
         case 0:
            // infinity-norm, result is a diamond/square with diagonals aligned on the axes in 2D
            radius = fabs(x-center[0]) + fabs(y-center[1]) + fabs(z-center[2]);
            break;
         case 1:
            // 1-norm, result is is a grid-aligned square in 2D
            radius = max(max(fabs(x-center[0]), fabs(y-center[1])), fabs(z-center[2]));
            break;
         case 2:
            // 2-norm (Cartesian), result is a circle in 2D
            radius = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
            break;
         case 3:
            // cylinder aligned with y-axis, use with polar plane/line dipole
            radius = sqrt((x-center[0])*(x-center[0]) + (z-center[2])*(z-center[2]));
            break;
         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3." << std::endl;
            abort();
      }
      
      if(radius < sP.ionosphereTaperRadius) {
         // linear tapering
         //initRho = this->ionosphereRho - (ionosphereRho-tailRho)*(radius-this->ionosphereRadius) / (this->ionosphereTaperRadius-this->ionosphereRadius);
         
         // sine tapering
         initRho = sP.rho - (sP.rho-sP.ionosphereRho)*0.5*(1.0+sin(M_PI*(radius-this->ionosphereRadius)/(sP.ionosphereTaperRadius-this->ionosphereRadius)+0.5*M_PI));
         if(radius < this->ionosphereRadius) {
            // Just to be safe, there are observed cases where this failed.
            initRho = sP.ionosphereRho;
         }
      }

      Real mass = getObjectWrapper().particleSpecies[popID].mass;

      return initRho * pow(mass / (2.0 * M_PI * physicalconstants::K_B * sP.T), 1.5) *
      exp(- mass * ((vx-initV0[0])*(vx-initV0[0]) + (vy-initV0[1])*(vy-initV0[1]) + (vz-initV0[2])*(vz-initV0[2])) / (2.0 * physicalconstants::K_B * sP.T));
   }

   vector<std::array<Real, 3> > Magnetosphere::getV0(
      creal x,
      creal y,
      creal z,
      const uint popID
   ) const {
      const MagnetosphereSpeciesParameters& sP = this->speciesParams[popID];

      vector<std::array<Real, 3> > centerPoints;
      std::array<Real, 3> V0 {{sP.V0[0], sP.V0[1], sP.V0[2]}};
      std::array<Real, 3> ionosphereV0 = {{sP.ionosphereV0[0], sP.ionosphereV0[1], sP.ionosphereV0[2]}};
      
      Real radius;
      
      switch(this->ionosphereGeometry) {
         case 0:
            // infinity-norm, result is a diamond/square with diagonals aligned on the axes in 2D
            radius = fabs(x-center[0]) + fabs(y-center[1]) + fabs(z-center[2]);
            break;
         case 1:
            // 1-norm, result is is a grid-aligned square in 2D
            radius = max(max(fabs(x-center[0]), fabs(y-center[1])), fabs(z-center[2]));
            break;
         case 2:
            // 2-norm (Cartesian), result is a circle in 2D
            radius = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
            break;
         case 3:
            // cylinder aligned with y-axis, use with polar plane/line dipole
            radius = sqrt((x-center[0])*(x-center[0]) + (z-center[2])*(z-center[2]));
            break;
         default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3." << std::endl;
            abort();
      }
      
      if(radius < sP.ionosphereTaperRadius) {
         // linear tapering
         //initV0[i] *= (radius-this->ionosphereRadius) / (this->ionosphereTaperRadius-this->ionosphereRadius);
         
         // sine tapering
         Real q=0.5*(1.0-sin(M_PI*(radius-this->ionosphereRadius)/(sP.ionosphereTaperRadius-this->ionosphereRadius)+0.5*M_PI));
         
         for(uint i=0; i<3; i++) {
            V0[i]=q*(V0[i]-ionosphereV0[i])+ionosphereV0[i];
            if(radius < this->ionosphereRadius) {
               // Just to be safe, there are observed cases where this failed.
               V0[i] = ionosphereV0[i];
            }
         }
      }
      
      centerPoints.push_back(V0);
      return centerPoints;
   }

   bool Magnetosphere::refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const {
 
     int myRank;       
     MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

     // mpiGrid.set_maximum_refinement_level(std::min(this->maxSpatialRefinementLevel, mpiGrid.mapping.get_maximum_refinement_level()));

     std::vector<CellID> refinedCells;

      // cout << "I am at line " << __LINE__ << " of " << __FILE__ <<  endl;
     if(myRank == MASTER_RANK) std::cout << "Maximum refinement level is " << mpiGrid.mapping.get_maximum_refinement_level() << std::endl;
      
     // Leave boundary cells and a bit of safety margin
     const int bw = 2* VLASOV_STENCIL_WIDTH;
     const int bw2 = 2*(bw + VLASOV_STENCIL_WIDTH);
     const int bw3 = 2*(bw2 + VLASOV_STENCIL_WIDTH);

     // Calculate regions for refinement
     if (P::amrMaxSpatialRefLevel > 0) {

	// L1 refinement.
	for (uint i = bw; i < P::xcells_ini-bw; ++i) {
	   for (uint j = bw; j < P::ycells_ini-bw; ++j) {
	      for (uint k = bw; k < P::zcells_ini-bw; ++k) {
		 
		 std::array<double,3> xyz;
		 xyz[0] = P::xmin + (i+0.5)*P::dx_ini;
		 xyz[1] = P::ymin + (j+0.5)*P::dy_ini;
		 xyz[2] = P::zmin + (k+0.5)*P::dz_ini;
                 
		 Real radius2 = (xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
		 // Check if cell is within L1 sphere, or within L1 tail slice
		 if ((radius2 < refine_L1radius*refine_L1radius) ||
                     ((xyz[0] < 0) && (std::abs(xyz[1]) < refine_L1radius) && 
		      (std::abs(xyz[2])<refine_L1tailthick)))
		    {
		       CellID myCell = mpiGrid.get_existing_cell(xyz);
		       mpiGrid.refine_completely(myCell);
		    }
	      }
	   }
	}
	refinedCells = mpiGrid.stop_refining(true);      
	if(myRank == MASTER_RANK) std::cout << "Finished first level of refinement" << endl;
#ifndef NDEBUG
	if(refinedCells.size() > 0) {
	   std::cout << "Rank " << myRank << " refined " << refinedCells.size() << " cells. " << std::endl;
	}
#endif
	mpiGrid.balance_load();
     }
     
     if (P::amrMaxSpatialRefLevel > 1) {
	
	// L2 refinement.
	for (uint i = bw2; i < 2*P::xcells_ini-bw2; ++i) {
	   for (uint j = bw2; j < 2*P::ycells_ini-bw2; ++j) {
	      for (uint k = bw2; k < 2*P::zcells_ini-bw2; ++k) {
		 
		 std::array<double,3> xyz;
		 xyz[0] = P::xmin + (i+0.5)*0.5*P::dx_ini;
		 xyz[1] = P::ymin + (j+0.5)*0.5*P::dy_ini;
		 xyz[2] = P::zmin + (k+0.5)*0.5*P::dz_ini;
                 
		 Real radius2 = (xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
		 // Check if cell is within L1 sphere, or within L1 tail slice
		 if ((radius2 < refine_L2radius*refine_L2radius) ||
		     ((xyz[0] < 0) && (std::abs(xyz[1]) < refine_L2radius) && 
		      (std::abs(xyz[2])<refine_L2tailthick)))
		    {
		       CellID myCell = mpiGrid.get_existing_cell(xyz);
		       // Check if the cell is tagged as do not compute
		       mpiGrid.refine_completely(myCell);
		    }
	      }
	   }
	}
	refinedCells = mpiGrid.stop_refining(true);
	if(myRank == MASTER_RANK) std::cout << "Finished second level of refinement" << endl;
#ifndef NDEBUG
	if(refinedCells.size() > 0) {
	   std::cout << "Rank " << myRank << " refined " << refinedCells.size() << " cells. " << std::endl;
	}
#endif
	
	mpiGrid.balance_load();
     }
     
     if (P::amrMaxSpatialRefLevel > 2) {
	
	//	if (refine_L3radius < refine_L2radius && refine_L3radius > ionosphereRadius) {
	// L3 refinement.
	   for (uint i = bw3; i < 4*P::xcells_ini-bw3; ++i) {
	      for (uint j = bw3; j < 4*P::ycells_ini-bw3; ++j) {
		 for (uint k = bw3; k < 4*P::zcells_ini-bw3; ++k) {
		    
		    std::array<double,3> xyz;
		    xyz[0] = P::xmin + (i+0.5)*0.25*P::dx_ini;
		    xyz[1] = P::ymin + (j+0.5)*0.25*P::dy_ini;
		    xyz[2] = P::zmin + (k+0.5)*0.25*P::dz_ini;
                    
 		    Real radius2 = (xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
// 		    // Check if cell is within L1 sphere, or within L1 tail slice
// 		    if (radius2 < refine_L3radius*refine_L3radius)
// 		       {
// 			  CellID myCell = mpiGrid.get_existing_cell(xyz);
// 			  // Check if the cell is tagged as do not compute
// 			  mpiGrid.refine_completely(myCell);
// 		       }

		    // Check if cell is within the nose cap
		    if ((xyz[0]>refine_L3nosexmin) && (radius2<refine_L3radius*refine_L3radius))
		       {
			  CellID myCell = mpiGrid.get_existing_cell(xyz);
			  // Check if the cell is tagged as do not compute
			  mpiGrid.refine_completely(myCell);			  
		       }

		    // Check if cell is within the tail box
		    if ((xyz[0]>refine_L3tailxmin) && (xyz[0]<refine_L3tailxmax) &&
			(abs(xyz[1])<refine_L3tailwidth) && (abs(xyz[2])<refine_L3tailheight))
		       {
			  CellID myCell = mpiGrid.get_existing_cell(xyz);
			  // Check if the cell is tagged as do not compute
			  mpiGrid.refine_completely(myCell);
		       }

 		 }
	      }
	   }
	   refinedCells = mpiGrid.stop_refining(true);
	   if(myRank == MASTER_RANK) std::cout << "Finished third level of refinement" << endl;
#ifndef NDEBUG
	   if(refinedCells.size() > 0) {
	      std::cout << "Rank " << myRank << " refined " << refinedCells.size() << " cells. " << std::endl;
	   }
#endif
	   
	   mpiGrid.balance_load();
        } else {
	   std::cout << "Skipping third level of refinement because the radius is larger than the 2nd level radius or smaller than the ionosphere radius." << std::endl;
        }
     //}
     
     return true;
   }
   
} // namespace projects

