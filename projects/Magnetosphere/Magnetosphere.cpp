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
#include "../../sysboundary/ionosphere.h"

#include "Magnetosphere.h"
#include "../../fieldsolver/derivatives.hpp"

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

      RP::add("Magnetosphere.refine_L4radius","Radius of L3-refined sphere or cap", 6.0e7);
      RP::add("Magnetosphere.refine_L4nosexmin","Low x-value of nose L3-refined box", 5.5e7);

      RP::add("Magnetosphere.refine_L3radius","Radius of L3-refined sphere or cap", 6.371e7); // 10 RE
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
      //New Parameter for zeroing out derivativeNew Parameter for zeroing out derivativess
      RP::add("Magnetosphere.zeroOutDerivativesX","Zero Out Perpendicular components", 1.0);
      RP::add("Magnetosphere.zeroOutDerivativesY","Zero Out Perpendicular components", 1.0);
      RP::add("Magnetosphere.zeroOutDerivativesZ","Zero Out Perpendicular components", 1.0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_Magnetosphere.rho", "Tail region number density (m^-3)", 0.0);
         RP::add(pop + "_Magnetosphere.T", "Temperature (K)", 0.0);
         RP::add(pop + "_Magnetosphere.VX0", "Initial bulk velocity in x-direction", 0.0);
         RP::add(pop + "_Magnetosphere.VY0", "Initial bulk velocity in y-direction", 0.0);
         RP::add(pop + "_Magnetosphere.VZ0", "Initial bulk velocity in z-direction", 0.0);
         RP::add(pop + "_Magnetosphere.taperInnerRadius", "Inner radius of the zone with a density tapering from the ionospheric value to the background (m)", 0.0);
         RP::add(pop + "_Magnetosphere.taperOuterRadius", "Outer radius of the zone with a density tapering from the ionospheric value to the background (m)", 0.0);
      }
   }
   
   void Magnetosphere::getParameters(){
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

      Project::getParameters();
      SysBoundary& sysBoundaryContainer = getObjectWrapper().sysBoundaryContainer;

      Real dummy;
      typedef Readparameters RP;
      RP::get("Magnetosphere.constBgBX", this->constBgB[0]);
      RP::get("Magnetosphere.constBgBY", this->constBgB[1]);
      RP::get("Magnetosphere.constBgBZ", this->constBgB[2]);
      RP::get("Magnetosphere.noDipoleInSW", dummy);
      this->noDipoleInSW = dummy == 1 ? true:false;
      RP::get("Magnetosphere.dipoleScalingFactor", this->dipoleScalingFactor);

      RP::get("Magnetosphere.dipoleMirrorLocationX", this->dipoleMirrorLocationX);

      RP::get("Magnetosphere.dipoleType", this->dipoleType);

      /* Enforce "dipole" (incl. correction terms) in solar wind with dipole type 4. */
      if ((this->dipoleType == 4) && (this->noDipoleInSW)) {
         if(myRank == MASTER_RANK) {
            std::cerr<<"Note: Initializing Magnetosphere with dipole type 4, which requires the dipole + vector potential "
            <<"correction terms in the solar wind. Thus overriding the config and setting Magnetosphere.noDipoleInSW=0."<<std::endl;
         }
         this->noDipoleInSW = false;
      }

      /** Read inner boundary parameters from either ionospheric or copysphere sysboundary condition */
      if (sysBoundaryContainer.existSysBoundary("Copysphere")) {
         RP::get("copysphere.radius", this->ionosphereRadius);
         RP::get("copysphere.centerX", this->center[0]);
         RP::get("copysphere.centerY", this->center[1]);
         RP::get("copysphere.centerZ", this->center[2]);
         RP::get("copysphere.geometry", this->ionosphereGeometry);
      } else if (sysBoundaryContainer.existSysBoundary("Ionosphere")) {
         RP::get("ionosphere.radius", this->ionosphereRadius);
         RP::get("ionosphere.centerX", this->center[0]);
         RP::get("ionosphere.centerY", this->center[1]);
         RP::get("ionosphere.centerZ", this->center[2]);
         RP::get("ionosphere.geometry", this->ionosphereGeometry);
      } else {
         if(myRank == MASTER_RANK) {
            std::cerr<<"Warning in initializing Magnetosphere: Could not find inner boundary (ionosphere or copysphere)!"<<std::endl;
         }
      }
      if(ionosphereRadius < 1000.) {
         // For really small ionospheric radius values, assume R_E units
         ionosphereRadius *= physicalconstants::R_E;
      }

      RP::get("Magnetosphere.refine_L4radius", this->refine_L4radius);
      RP::get("Magnetosphere.refine_L4nosexmin", this->refine_L4nosexmin);

      RP::get("Magnetosphere.refine_L3radius", this->refine_L3radius);
      RP::get("Magnetosphere.refine_L3nosexmin", this->refine_L3nosexmin);
      RP::get("Magnetosphere.refine_L3tailwidth", this->refine_L3tailwidth);
      RP::get("Magnetosphere.refine_L3tailheight", this->refine_L3tailheight);
      RP::get("Magnetosphere.refine_L3tailxmin", this->refine_L3tailxmin);
      RP::get("Magnetosphere.refine_L3tailxmax", this->refine_L3tailxmax);

      RP::get("Magnetosphere.refine_L2radius", this->refine_L2radius);
      RP::get("Magnetosphere.refine_L2tailthick", this->refine_L2tailthick);
      RP::get("Magnetosphere.refine_L1radius", this->refine_L1radius);
      RP::get("Magnetosphere.refine_L1tailthick", this->refine_L1tailthick);

      RP::get("Magnetosphere.dipoleTiltPhi", this->dipoleTiltPhi);
      RP::get("Magnetosphere.dipoleTiltTheta", this->dipoleTiltTheta);
      RP::get("Magnetosphere.dipoleXFull", this->dipoleXFull);
      RP::get("Magnetosphere.dipoleXZero", this->dipoleXZero);
      RP::get("Magnetosphere.dipoleInflowBX", this->dipoleInflowB[0]);
      RP::get("Magnetosphere.dipoleInflowBY", this->dipoleInflowB[1]);
      RP::get("Magnetosphere.dipoleInflowBZ", this->dipoleInflowB[2]);

      RP::get("Magnetosphere.zeroOutDerivativesX", this->zeroOutComponents[0]);
      RP::get("Magnetosphere.zeroOutDerivativesY", this->zeroOutComponents[1]);
      RP::get("Magnetosphere.zeroOutDerivativesZ", this->zeroOutComponents[2]);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         MagnetosphereSpeciesParameters sP;

         RP::get(pop + "_Magnetosphere.rho", sP.rho);
         RP::get(pop + "_Magnetosphere.T", sP.T);
         RP::get(pop + "_Magnetosphere.VX0", sP.V0[0]);
         RP::get(pop + "_Magnetosphere.VY0", sP.V0[1]);
         RP::get(pop + "_Magnetosphere.VZ0", sP.V0[2]);

         /** Read inner boundary parameters from either ionospheric or copysphere sysboundary condition */
         if (sysBoundaryContainer.existSysBoundary("Copysphere")) {
            RP::get(pop + "_copysphere.rho", sP.ionosphereRho);
            RP::get(pop + "_copysphere.T", sP.ionosphereT);
            RP::get(pop + "_copysphere.VX0", sP.ionosphereV0[0]);
            RP::get(pop + "_copysphere.VY0", sP.ionosphereV0[1]);
            RP::get(pop + "_copysphere.VZ0", sP.ionosphereV0[2]);
         } else if (sysBoundaryContainer.existSysBoundary("Ionosphere")) {
            RP::get(pop + "_ionosphere.rho", sP.ionosphereRho);
            RP::get(pop + "_ionosphere.T", sP.ionosphereT);
            RP::get(pop + "_ionosphere.VX0", sP.ionosphereV0[0]);
            RP::get(pop + "_ionosphere.VY0", sP.ionosphereV0[1]);
            RP::get(pop + "_ionosphere.VZ0", sP.ionosphereV0[2]);
         }
         RP::get(pop + "_Magnetosphere.taperInnerRadius", sP.taperInnerRadius);
         RP::get(pop + "_Magnetosphere.taperOuterRadius", sP.taperOuterRadius);
         // Backward-compatibility: cfgs from before Sep 2021 setting pop_ionosphere.taperRadius will fail with the unknown option.
         // Some fail-safety checks
         if(sP.taperInnerRadius < 0 || sP.taperOuterRadius < 0) {
            if(myRank == MASTER_RANK) {
               cerr << "Error: " << pop << "_Magnetosphere.taperInnerRadius and tapeOuterRadius should be >= 0! Aborting." << endl;
            }
            abort();
         }
         if(sP.taperInnerRadius > sP.taperOuterRadius) {
            if(myRank == MASTER_RANK) {
               cerr << "Error: " << pop << "_Magnetosphere.taperInnerRadius should be <= taperOuterRadius! Aborting." << endl;
            }
            abort();
         }
         if(sP.taperOuterRadius > 0 && sP.taperOuterRadius <= this->ionosphereRadius) {
            if(myRank == MASTER_RANK) {
               cerr << "Error: " << pop << "_Magnetosphere.taperOuterRadius is non-zero yet smaller than ionosphere.radius / copysphere.radius! Aborting." << endl;
            }
            abort();
         }
         if(sP.taperInnerRadius == 0 && sP.taperOuterRadius > 0) {
            if(myRank == MASTER_RANK) {
               cerr << "Warning: " << pop << "_Magnetosphere.taperInnerRadius is zero (default), now setting this to the same value as ionosphere.radius / copysphere.radius, that is " << this->ionosphereRadius << ". Set/change " << pop << "_Magnetosphere.taperInnerRadius if this is not the expected behavior." << endl;
            }
            sP.taperInnerRadius = this->ionosphereRadius;
         }
         if(sP.ionosphereT == 0) {
            if(myRank == MASTER_RANK) {
               if (sysBoundaryContainer.existSysBoundary("Copysphere")) {
                  cerr << "Warning: " << pop << "_copysphere.T is zero (default), now setting to the same value as " << pop << "_Magnetosphere.T, that is " << sP.T << ". Set/change " << pop << "_copysphere.T if this is not the expected behavior." << endl;
               } else if (sysBoundaryContainer.existSysBoundary("Ionosphere")) {
                  cerr << "Warning: " << pop << "_ionosphere.T is zero (default), now setting to the same value as " << pop << "_Magnetosphere.T, that is " << sP.T << ". Set/change " << pop << "_ionosphere.T if this is not the expected behavior." << endl;
               }
            }
            sP.ionosphereT = sP.T;
         }
         if(sP.ionosphereRho == 0) {
            if(myRank == MASTER_RANK) {
               if (sysBoundaryContainer.existSysBoundary("Copysphere")) {
                  cerr << "Warning: " << pop << "_copysphere.rho is zero (default), now setting to the same value as " << pop << "_Magnetosphere.rho, that is " << sP.rho << ". Set/change " << pop << "_copysphere.rho if this is not the expected behavior." << endl;
               } else if (sysBoundaryContainer.existSysBoundary("Ionosphere")) {
                  cerr << "Warning: " << pop << "_ionosphere.rho is zero (default), now setting to the same value as " << pop << "_Magnetosphere.rho, that is " << sP.rho << ". Set/change " << pop << "_ionosphere.rho if this is not the expected behavior." << endl;
               }
            }
            sP.ionosphereRho = sP.rho;
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

      return getDistribValue(x+0.5*dx,y+0.5*dy,z+0.5*dz,vx+0.5*dvx,vy+0.5*dvy,vz+0.5*dvz,dvx,dvy,dvz,popID);
   }
   
   /*! Magnetosphere does not set any extra perturbed B. */
   void Magnetosphere::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   /* set 0-centered dipole */
   void Magnetosphere::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      Dipole bgFieldDipole;
      LineDipole bgFieldLineDipole;
      VectorDipole bgVectorDipole;

      phiprof::Timer switchDipoleTypeTimer {"switch-dipoleType"};
      // The hardcoded constants of dipole and line dipole moments are obtained
      // from Daldorff et al (2014), see
      // https://github.com/fmihpc/vlasiator/issues/20 for a derivation of the
      // values used here.
      switch(this->dipoleType) {
            case 0:
               bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
               setBackgroundField(bgFieldDipole, BgBGrid);
               SBC::ionosphereGrid.setDipoleField(bgFieldDipole);
               break;
            case 1:
               bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
               setBackgroundField(bgFieldLineDipole, BgBGrid);
               SBC::ionosphereGrid.setDipoleField(bgFieldLineDipole);
               break;
            case 2:
               bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
               setBackgroundField(bgFieldLineDipole, BgBGrid);
               //Append mirror dipole
               bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, this->dipoleMirrorLocationX, 0.0, 0.0 );
               setBackgroundField(bgFieldLineDipole, BgBGrid, true);
               SBC::ionosphereGrid.setDipoleField(bgFieldLineDipole);
               break;
            case 3:
               bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
               setBackgroundField(bgFieldDipole, BgBGrid);
               SBC::ionosphereGrid.setDipoleField(bgFieldDipole);
               //Append mirror dipole                
               bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, this->dipoleMirrorLocationX, 0.0, 0.0, 0.0 );//mirror
               setBackgroundField(bgFieldDipole, BgBGrid, true);
               break; 
            case 4:  // Vector potential dipole, vanishes or optionally scales to static inflow value after a given x-coordinate
               // What we in fact do is we place the regular dipole in the background field, and the
               // corrective terms in the perturbed field. This maintains the BGB as curl-free.
               bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 ); //set dipole moment
               setBackgroundField(bgFieldDipole, BgBGrid);
               SBC::ionosphereGrid.setDipoleField(bgFieldDipole);
               // Now we calculate the difference required to scale the dipole to zero as we approach the inflow,
               // and store it inside the BgBGrid object for use by e.g. boundary conditions.
               bgFieldDipole.initialize(-8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );
               setPerturbedField(bgFieldDipole, BgBGrid, fsgrids::bgbfield::BGBXVDCORR);
               bgVectorDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, this->dipoleTiltPhi*M_PI/180., this->dipoleTiltTheta*M_PI/180., this->dipoleXFull, this->dipoleXZero, this->dipoleInflowB[0], this->dipoleInflowB[1], this->dipoleInflowB[2]);
               setPerturbedField(bgVectorDipole, BgBGrid, fsgrids::bgbfield::BGBXVDCORR, true);
               if (P::isRestart == false) {
                  // If we are starting a new simulation, we also copy this data into perB.
                  const auto localSize = BgBGrid.getLocalSize().data();
                  #pragma omp parallel for collapse(2)
                  for (int z = 0; z < localSize[2]; ++z) {
                     for (int y = 0; y < localSize[1]; ++y) {
                        for (int x = 0; x < localSize[0]; ++x) {
                           std::array<Real, fsgrids::bgbfield::N_BGB>* BGBcell = BgBGrid.get(x, y, z);
                           std::array<Real, fsgrids::bfield::N_BFIELD>* PERBcell = perBGrid.get(x, y, z);
                           PERBcell->at(fsgrids::bfield::PERBX) = BGBcell->at(fsgrids::bgbfield::BGBXVDCORR);
                           PERBcell->at(fsgrids::bfield::PERBY) = BGBcell->at(fsgrids::bgbfield::BGBYVDCORR);
                           PERBcell->at(fsgrids::bfield::PERBZ) = BGBcell->at(fsgrids::bgbfield::BGBZVDCORR);
                        }
                     }
                  }
               }
               break;
            default:
               setBackgroundFieldToZero(BgBGrid);
      }
      switchDipoleTypeTimer.stop();

      const auto localSize = BgBGrid.getLocalSize().data();
      
      phiprof::Timer zeroingTimer {"zeroing-out"};

#pragma omp parallel
      {
         bool doZeroOut;
         //Force field to zero in the perpendicular direction for 2D (1D) simulations. Otherwise we have unphysical components.
         doZeroOut = P::xcells_ini ==1 && this->zeroOutComponents[0]==1;
      
         if(doZeroOut) {
#pragma omp for collapse(2)
            for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
               for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
                  for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
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

          doZeroOut = P::ycells_ini ==1 && this->zeroOutComponents[1]==1;
          if(doZeroOut) {
             /*2D simulation in x and z. Set By and derivatives along Y, and derivatives of By to zero*/
#pragma omp for collapse(2)
             for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
                   for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
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

         doZeroOut = P::zcells_ini ==1 && this->zeroOutComponents[2]==1;
         if(doZeroOut) {
#pragma omp for collapse(2)
            for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
               for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
                  for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
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
#pragma omp for collapse(2)
            for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
               for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
                  for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
                     if(technicalGrid.get(x, y, z)->sysBoundaryFlag == sysboundarytype::MAXWELLIAN ) {
                        for (int i = 0; i < fsgrids::bgbfield::N_BGB; ++i) {
                           BgBGrid.get(x,y,z)->at(i) = 0;
                        }
                        if ( (this->dipoleType==4) && (P::isRestart == false) ) {
                           // If we set BGB to zero here, then we should also set perB in new runs to zero.
                           for (int i = 0; i < fsgrids::bfield::N_BFIELD; ++i) {
                              perBGrid.get(x,y,z)->at(i) = 0;
                           }
                        }
                     }
                  }
               }
            }
         }
      } // end of omp parallel region

      zeroingTimer.stop();

      phiprof::Timer addConstantTimer {"add-constant-field"};
      // Superimpose constant background field if needed
      if(this->constBgB[0] != 0.0 || this->constBgB[1] != 0.0 || this->constBgB[2] != 0.0) {
         ConstantField bgConstantField;
         bgConstantField.initialize(this->constBgB[0], this->constBgB[1], this->constBgB[2]);
         setBackgroundField(bgConstantField, BgBGrid, true);
         SBC::ionosphereGrid.setConstantBackgroundField(this->constBgB);
      }
      addConstantTimer.stop();
      phiprof::Timer storeNodeTimer {"ionosphereGrid.storeNodeB"};
      SBC::ionosphereGrid.storeNodeB();
      storeNodeTimer.stop();
   }
   
   
   Real Magnetosphere::getDistribValue(
           creal& x,creal& y,creal& z,
           creal& vx,creal& vy,creal& vz,
           creal& dvx,creal& dvy,creal& dvz,
           const uint popID) const
   {
      const MagnetosphereSpeciesParameters& sP = this->speciesParams[popID];
      Real initRho = sP.rho;
      Real initT = sP.T;
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
      
      if(radius < sP.taperOuterRadius) {
         // sine tapering
         initRho = sP.rho - (sP.rho-sP.ionosphereRho)*0.5*(1.0+sin(M_PI*(radius-sP.taperInnerRadius)/(sP.taperOuterRadius-sP.taperInnerRadius)+0.5*M_PI));
         initT = sP.T - (sP.T-sP.ionosphereT)*0.5*(1.0+sin(M_PI*(radius-sP.taperInnerRadius)/(sP.taperOuterRadius-sP.taperInnerRadius)+0.5*M_PI));
         if(radius <= sP.taperInnerRadius) {
            initRho = sP.ionosphereRho;
            initT = sP.ionosphereT;
         }
      }

      Real mass = getObjectWrapper().particleSpecies[popID].mass;

      return initRho * pow(mass / (2.0 * M_PI * physicalconstants::K_B * initT), 1.5) *
      exp(- mass * ((vx-initV0[0])*(vx-initV0[0]) + (vy-initV0[1])*(vy-initV0[1]) + (vz-initV0[2])*(vz-initV0[2])) / (2.0 * physicalconstants::K_B * initT));
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
      
      if(radius < sP.taperOuterRadius) {
         // sine tapering
         Real q=0.5*(1.0-sin(M_PI*(radius-sP.taperInnerRadius)/(sP.taperOuterRadius-sP.taperInnerRadius)+0.5*M_PI));
         
         for(uint i=0; i<3; i++) {
            V0[i]=q*(V0[i]-ionosphereV0[i])+ionosphereV0[i];
            if(radius <= sP.taperInnerRadius) {
               V0[i] = ionosphereV0[i];
            }
         }
      }
      
      centerPoints.push_back(V0);
      return centerPoints;
   }

   bool Magnetosphere::refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const {
      phiprof::Timer refineSCTimer {"Magnetosphere: refine spatial cells"};
   
      int myRank;       
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

      if(myRank == MASTER_RANK) {
         std::cout << "Maximum refinement level is " << mpiGrid.mapping.get_maximum_refinement_level() << std::endl;
      }

      std::vector<CellID> cells = getLocalCells();

      // L1 refinement.
      if (P::amrMaxSpatialRefLevel > 0 && P::amrMaxAllowedSpatialRefLevel > 0) {
         //#pragma omp parallel for
         for (uint i = 0; i < cells.size(); ++i) {
            CellID id = cells[i];
            std::array<double,3> xyz = mpiGrid.get_center(id);
                     
            Real radius2 = pow(xyz[0], 2) + pow(xyz[1], 2) + pow(xyz[2], 2);
            bool inSphere = radius2 < refine_L1radius*refine_L1radius;
            bool inTail = xyz[0] < 0 && fabs(xyz[1]) < refine_L1radius && fabs(xyz[2]) < refine_L1tailthick;
            if ((inSphere || inTail) && radius2 < P::refineRadius * P::refineRadius) {
               //#pragma omp critical
               mpiGrid.refine_completely(id);
            }
         }

         cells = mpiGrid.stop_refining();      
         if (myRank == MASTER_RANK) {
            std::cout << "Finished first level of refinement" << endl;
         }
         #ifndef NDEBUG
         if (cells.size() > 0) {
            std::cout << "Rank " << myRank << " refined " << cells.size() << " cells to level 1" << std::endl;
         }
         #endif //NDEBUG
      }
      
      // L2 refinement.
      if (P::amrMaxSpatialRefLevel > 1 && P::amrMaxAllowedSpatialRefLevel > 1) {
         //#pragma omp parallel for
         for (uint i = 0; i < cells.size(); ++i) {
            CellID id = cells[i];
            std::array<double,3> xyz = mpiGrid.get_center(id);
                     
            Real radius2 = pow(xyz[0], 2) + pow(xyz[1], 2) + pow(xyz[2], 2);
            bool inSphere = radius2 < pow(refine_L2radius, 2);
            bool inTail = xyz[0] < 0 && fabs(xyz[1]) < refine_L2radius && fabs(xyz[2])<refine_L2tailthick;
            if ((inSphere || inTail) && radius2 < P::refineRadius * P ::refineRadius) {
               //#pragma omp critical
               mpiGrid.refine_completely(id);
            }
         }
         cells = mpiGrid.stop_refining();
         if(myRank == MASTER_RANK) {
            std::cout << "Finished second level of refinement" << endl;
         }
         #ifndef NDEBUG
         if (cells.size() > 0) {
            std::cout << "Rank " << myRank << " refined " << cells.size() << " cells to level 2" << std::endl;
         }
         #endif //NDEBUG

      }
      
      // L3 refinement.
      if (P::amrMaxSpatialRefLevel > 2 && P::amrMaxAllowedSpatialRefLevel > 2) {
         //#pragma omp parallel for
         for (uint i = 0; i < cells.size(); ++i) {
            CellID id = cells[i];
            std::array<double,3> xyz = mpiGrid.get_center(id);
                     
            Real radius2 = pow(xyz[0], 2) + pow(xyz[1], 2) + pow(xyz[2], 2);
            bool inNoseCap = (xyz[0]>refine_L3nosexmin) && (radius2<refine_L3radius*refine_L3radius);
            bool inTail = (xyz[0]>refine_L3tailxmin) && (xyz[0]<refine_L3tailxmax) && (fabs(xyz[1])<refine_L3tailwidth) && (fabs(xyz[2])<refine_L3tailheight);
            if ((inNoseCap || inTail) && radius2 < P::refineRadius * P::refineRadius) {
               //#pragma omp critical
               mpiGrid.refine_completely(id);			  
            }
         }
         cells = mpiGrid.stop_refining();
         if (myRank == MASTER_RANK) {
            std::cout << "Finished third level of refinement" << endl;
         }
         #ifndef NDEBUG
         if (cells.size() > 0) {
            std::cout << "Rank " << myRank << " refined " << cells.size() << " cells to level 3" << std::endl;
         }
         #endif //NDEBUG
      }

      // L4 refinement.
      if (P::amrMaxSpatialRefLevel > 3 && P::amrMaxAllowedSpatialRefLevel > 3) {
         //#pragma omp parallel for
         for (uint i = 0; i < cells.size(); ++i) {
            CellID id = cells[i];
            std::array<double,3> xyz = mpiGrid.get_center(id);
                     
            Real radius2 = (xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);

            // Check if cell is within the nose cap
            bool inNose = refine_L4nosexmin && radius2<refine_L4radius*refine_L4radius;
            if (inNose && radius2 < P::refineRadius * P::refineRadius) {
               //#pragma omp critical
               mpiGrid.refine_completely(id);			  
            }
         }

         cells = mpiGrid.stop_refining();
         if (myRank == MASTER_RANK) {
            std::cout << "Finished fourth level of refinement" << endl;
         }
         #ifndef NDEBUG
         if (cells.size() > 0) {
            std::cout << "Rank " << myRank << " refined " << cells.size() << " cells to level 4" << std::endl;
         }
         #endif //NDEBUG
      }

      return true;
   }

   bool Magnetosphere::forceRefinement( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, int n ) const {
   
      int myRank;       
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

      if(myRank == MASTER_RANK) {
         std::cout << "Maximum refinement level is " << mpiGrid.mapping.get_maximum_refinement_level() << std::endl;
      }

      for (CellID id : getLocalCells()) {
         std::array<double,3> xyz {mpiGrid.get_center(id)};
         Real radius2 {pow(xyz[0], 2) + pow(xyz[1], 2) + pow(xyz[2], 2)};
         int refLevel {mpiGrid.get_refinement_level(id)};
         int refineTarget {0};

         if (P::amrMaxSpatialRefLevel > 0 && P::amrMaxAllowedSpatialRefLevel > 0) {
            bool inSphere = radius2 < refine_L1radius*refine_L1radius;
            bool inTail = xyz[0] < 0 && fabs(xyz[1]) < refine_L1radius && fabs(xyz[2]) < refine_L1tailthick;
            if ((inSphere || inTail) && radius2 < P::refineRadius * P ::refineRadius)
               ++refineTarget;
         }
         if (P::amrMaxSpatialRefLevel > 1 && P::amrMaxAllowedSpatialRefLevel > 1) {
            bool inSphere = radius2 < pow(refine_L2radius, 2);
            bool inTail = xyz[0] < 0 && fabs(xyz[1]) < refine_L2radius && fabs(xyz[2])<refine_L2tailthick;
            if ((inSphere || inTail) && radius2 < P::refineRadius * P ::refineRadius)
               ++refineTarget;
         }
         if (P::amrMaxSpatialRefLevel > 2 && P::amrMaxAllowedSpatialRefLevel > 2) {
            bool inNoseCap = (xyz[0]>refine_L3nosexmin) && (radius2<refine_L3radius*refine_L3radius);
            bool inTail = (xyz[0]>refine_L3tailxmin) && (xyz[0]<refine_L3tailxmax) && (fabs(xyz[1])<refine_L3tailwidth) && (fabs(xyz[2])<refine_L3tailheight);
            if ((inNoseCap || inTail) && radius2 < P::refineRadius * P ::refineRadius)
               ++refineTarget;
         }
         if (P::amrMaxSpatialRefLevel > 3 && P::amrMaxAllowedSpatialRefLevel > 3) {
            bool inNose = refine_L4nosexmin && radius2<refine_L4radius*refine_L4radius;
            if (inNose && radius2 < P::refineRadius * P ::refineRadius)
               ++refineTarget;
         }

         if (!canRefine(mpiGrid[id])) {
            mpiGrid.dont_refine(id);
            mpiGrid.dont_unrefine(id);
         } else if (refLevel <= n && refLevel < refineTarget) {
            mpiGrid.refine_completely(id);
         } else if (refLevel >= mpiGrid.mapping.get_maximum_refinement_level() - n && refLevel > refineTarget) {
            mpiGrid.unrefine_completely(id);
         } else {
            mpiGrid.dont_unrefine(id);
         }
      }

      return true;
   }
   
} // namespace projects

