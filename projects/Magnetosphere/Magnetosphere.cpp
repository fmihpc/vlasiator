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
#include <array>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/dipole.hpp"
#include "../../backgroundfield/linedipole.hpp"


#include "Magnetosphere.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   Magnetosphere::Magnetosphere(): TriAxisSearch() { }
   Magnetosphere::~Magnetosphere() { }
   
   void Magnetosphere::addParameters() {
      typedef Readparameters RP;
      RP::add("Magnetosphere.rho", "Tail region number density (m^-3)", 0.0);
      RP::add("Magnetosphere.T", "Temperature (K)", 0.0);
      RP::add("Magnetosphere.VX0", "Initial bulk velocity in x-direction", 0.0);
      RP::add("Magnetosphere.VY0", "Initial bulk velocity in y-direction", 0.0);
      RP::add("Magnetosphere.VZ0", "Initial bulk velocity in z-direction", 0.0);
      RP::add("Magnetosphere.constBgBX", "Constant flat Bx component in the whole simulation box. Default is none.", 0.0);
      RP::add("Magnetosphere.constBgBY", "Constant flat By component in the whole simulation box. Default is none.", 0.0);
      RP::add("Magnetosphere.constBgBZ", "Constant flat Bz component in the whole simulation box. Default is none.", 0.0);
      RP::add("Magnetosphere.noDipoleInSW", "If set to 1, the dipole magnetic field is not set in the solar wind inflow cells. Default 0.", 0.0);
      RP::add("Magnetosphere.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Magnetosphere.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
      RP::add("Magnetosphere.dipoleScalingFactor","Scales the field strength of the magnetic dipole compared to Earths.", 1.0);
      RP::add("Magnetosphere.dipoleType","0: Normal 3D dipole, 1: line-dipole for 2D polar simulations, 2: line-dipole with mirror, 3: 3D dipole with mirror", 0);
      RP::add("Magnetosphere.dipoleMirrorLocationX","x-coordinate of dipole Mirror", -1.0);
   }
   
   void Magnetosphere::getParameters(){
      Project::getParameters();
      
      int myRank;
      Real dummy;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      typedef Readparameters RP;
      if(!RP::get("Magnetosphere.rho", this->tailRho)) {
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

      if(!RP::get("Magnetosphere.nSpaceSamples", this->nSpaceSamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("Magnetosphere.nVelocitySamples", this->nVelocitySamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }

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
      
      if(!RP::get("ionosphere.rho", this->ionosphereRho)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.radius", this->ionosphereRadius)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.taperRadius", this->ionosphereTaperRadius)) {
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
      if(!RP::get("ionosphere.VX0", this->ionosphereV0[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.VY0", this->ionosphereV0[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ionosphere.VZ0", this->ionosphereV0[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }

   }
   
   bool Magnetosphere::initialize() {
      return Project::initialize();
   }

   Real Magnetosphere::calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
                                             creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,
                                             creal& dvz,const int& popID) const {
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
                           avg += getDistribValue(x+i*d_x,y+j*d_y,z+k*d_z,vx+vi*d_vx,vy+vj*d_vy,vz+vk*d_vz,dvx,dvy,dvz);
                        }
                        return avg /
                        (this->nSpaceSamples*this->nSpaceSamples*this->nSpaceSamples) /
                        (this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
      } else {
         return getDistribValue(x+0.5*dx,y+0.5*dy,z+0.5*dz,vx+0.5*dvx,vy+0.5*dvy,vz+0.5*dvz,dvx,dvy,dvz);
      }
   }
   
   /*! Magnetosphere does not set any extra perturbed B. */
   void Magnetosphere::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      cellParams[CellParams::PERBX] = 0.0;
      cellParams[CellParams::PERBY] = 0.0;
      cellParams[CellParams::PERBZ] = 0.0;
   }

   /* set 0-centered dipole */
   void Magnetosphere::setCellBackgroundField(SpatialCell *cell) const {
      if(cell->sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN && this->noDipoleInSW) {
         setBackgroundFieldToZero(cell->parameters, cell->derivatives,cell->derivativesBVOL);
      }
      else {
         Dipole bgFieldDipole;
         LineDipole bgFieldLineDipole;

         switch(this->dipoleType) {
             case 0:
                bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
                setBackgroundField(bgFieldDipole,cell->parameters, cell->derivatives,cell->derivativesBVOL);
                break;
             case 1:
                bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
                setBackgroundField(bgFieldLineDipole,cell->parameters, cell->derivatives,cell->derivativesBVOL);
                break;
             case 2:
                bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
                setBackgroundField(bgFieldLineDipole,cell->parameters, cell->derivatives,cell->derivativesBVOL);
                //Append mirror dipole
                bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, this->dipoleMirrorLocationX, 0.0, 0.0 );
                setBackgroundField(bgFieldLineDipole,cell->parameters, cell->derivatives,cell->derivativesBVOL, true);
                break;
             case 3:
                bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
                setBackgroundField(bgFieldDipole,cell->parameters, cell->derivatives,cell->derivativesBVOL);
                //Append mirror dipole                
                bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, this->dipoleMirrorLocationX, 0.0, 0.0, 0.0 );//mirror
                setBackgroundField(bgFieldDipole,cell->parameters, cell->derivatives,cell->derivativesBVOL, true);
                break;
                
             default:
                setBackgroundFieldToZero(cell->parameters, cell->derivatives,cell->derivativesBVOL);
                
         }
      }
      

      //Force field to zero in the perpendicular direction for 2D (1D) simulations. Otherwise we have unphysical components.
      if(P::xcells_ini==1) {
         cell->parameters[CellParams::BGBX]=0;
         cell->parameters[CellParams::BGBXVOL]=0.0;
         cell->derivatives[fieldsolver::dBGBydx]=0.0;
         cell->derivatives[fieldsolver::dBGBzdx]=0.0;
         cell->derivatives[fieldsolver::dBGBxdy]=0.0;
         cell->derivatives[fieldsolver::dBGBxdz]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBYVOLdx]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBZVOLdx]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBXVOLdy]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBXVOLdz]=0.0;
      }
      
      if(P::ycells_ini==1) {
         /*2D simulation in x and z. Set By and derivatives along Y, and derivatives of By to zero*/
         cell->parameters[CellParams::BGBY]=0.0;
         cell->parameters[CellParams::BGBYVOL]=0.0;
         cell->derivatives[fieldsolver::dBGBxdy]=0.0;
         cell->derivatives[fieldsolver::dBGBzdy]=0.0;
         cell->derivatives[fieldsolver::dBGBydx]=0.0;
         cell->derivatives[fieldsolver::dBGBydz]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBXVOLdy]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBZVOLdy]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBYVOLdx]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBYVOLdz]=0.0;
      }
      if(P::zcells_ini==1) {
         cell->parameters[CellParams::BGBX]=0;
         cell->parameters[CellParams::BGBY]=0;
         cell->parameters[CellParams::BGBYVOL]=0.0;
         cell->parameters[CellParams::BGBXVOL]=0.0;
         cell->derivatives[fieldsolver::dBGBxdy]=0.0;
         cell->derivatives[fieldsolver::dBGBxdz]=0.0;
         cell->derivatives[fieldsolver::dBGBydx]=0.0;
         cell->derivatives[fieldsolver::dBGBydz]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBXVOLdy]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBXVOLdz]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBYVOLdx]=0.0;
         cell->derivativesBVOL[bvolderivatives::dBGBYVOLdz]=0.0;
      }
      for(uint component=0; component<3; component++) {
         if(this->constBgB[component] != 0.0) {
            cell->parameters[CellParams::BGBX+component] += this->constBgB[component];
            cell->parameters[CellParams::BGBXVOL+component] += this->constBgB[component];
         }
      }
      
//       // FIXME TESTING HACK to be used when one wants to get the "zero" Hall field from the dipole
//       cell->parameters[CellParams::PERBX] = cell->parameters[CellParams::BGBX];
//       cell->parameters[CellParams::PERBXVOL] = cell->parameters[CellParams::BGBXVOL];
//       cell->parameters[CellParams::BGBX] = 0.0;
//       cell->parameters[CellParams::BGBXVOL] = 0.0;
//       cell->parameters[CellParams::PERBY] = cell->parameters[CellParams::BGBY];
//       cell->parameters[CellParams::PERBYVOL] = cell->parameters[CellParams::BGBYVOL];
//       cell->parameters[CellParams::BGBY] = 0.0;
//       cell->parameters[CellParams::BGBYVOL] = 0.0;
//       cell->parameters[CellParams::PERBZ] = cell->parameters[CellParams::BGBY];
//       cell->parameters[CellParams::PERBZVOL] = cell->parameters[CellParams::BGBZVOL];
//       cell->parameters[CellParams::BGBZ] = 0.0;
//       cell->parameters[CellParams::BGBZVOL] = 0.0;
//       // END OF TESTING HACK
      
   }
      
      
   Real Magnetosphere::getDistribValue(
           creal& x,creal& y,creal& z,
           creal& vx,creal& vy,creal& vz,
           creal& dvx,creal& dvy,creal& dvz) const 
   {
      Real initRho = this->tailRho;
      std::array<Real, 3> initV0 = this->getV0(x, y, z)[0];
      
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
      
      if(radius < this->ionosphereTaperRadius) {
         // linear tapering
         //initRho = this->ionosphereRho - (ionosphereRho-tailRho)*(radius-this->ionosphereRadius) / (this->ionosphereTaperRadius-this->ionosphereRadius);
         
         // sine tapering
         initRho = this->tailRho - (this->tailRho-this->ionosphereRho)*0.5*(1.0+sin(M_PI*(radius-this->ionosphereRadius)/(this->ionosphereTaperRadius-this->ionosphereRadius)+0.5*M_PI));
         if(radius < this->ionosphereRadius) {
            // Just to be safe, there are observed cases where this failed.
            initRho = this->ionosphereRho;
         }
      }
      
      return initRho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5) *
      exp(- physicalconstants::MASS_PROTON * ((vx-initV0[0])*(vx-initV0[0]) + (vy-initV0[1])*(vy-initV0[1]) + (vz-initV0[2])*(vz-initV0[2])) / (2.0 * physicalconstants::K_B * this->T));
   }
   
   vector<std::array<Real, 3> > Magnetosphere::getV0(
      creal x,
      creal y,
      creal z
   ) const {
      vector<std::array<Real, 3> > centerPoints;
      std::array<Real, 3> V0 {{this->V0[0], this->V0[1], this->V0[2]}};
      std::array<Real, 3> ionosphereV0 = {{this->ionosphereV0[0], this->ionosphereV0[1], this->ionosphereV0[2]}};
      
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
      
      if(radius < this->ionosphereTaperRadius) {
         // linear tapering
         //initV0[i] *= (radius-this->ionosphereRadius) / (this->ionosphereTaperRadius-this->ionosphereRadius);
         
         // sine tapering
         Real q=0.5*(1.0-sin(M_PI*(radius-this->ionosphereRadius)/(this->ionosphereTaperRadius-this->ionosphereRadius)+0.5*M_PI));
         
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
   
} // namespace projects

