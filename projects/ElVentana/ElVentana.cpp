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
#include "../../object_wrapper.h"
#include "../../ioread.h"
#include "../../memoryallocation.h"
#include "../../fieldsolver/gridGlue.hpp"
#include "../../tools/vlsvreaderinterface.h"

#include "ElVentana.h"
#include "../../grid.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   ElVentana::ElVentana(): TriAxisSearch() { }
   ElVentana::~ElVentana() { }

   void ElVentana::addParameters() {
      typedef Readparameters RP;
      
      RP::add("ElVentana.constBgBX", "Constant flat Bx component in the whole simulation box. Default is none.", 0.0);
      RP::add("ElVentana.constBgBY", "Constant flat By component in the whole simulation box. Default is none.", 0.0);
      RP::add("ElVentana.constBgBZ", "Constant flat Bz component in the whole simulation box. Default is none.", 0.0);
      RP::add("ElVentana.StartFile", "Restart file to be used to read the fields and population parameters.", "restart.0000000.vlsv");
      RP::add("ElVentana.WindowX_min", "Section boundary (X-direction) of the domain contained in the StartFile to be selected (meters).",6.0e7);
      RP::add("ElVentana.WindowX_max", "Section boundary (X-direction) of the domain contained in the StartFile to be selected (meters).",7.8e7);
      RP::add("ElVentana.WindowY_min", "Section boundary (Y-direction) of the domain contained in the StartFile to be selected (meters).",-.1e6);
      RP::add("ElVentana.WindowY_max", "Section boundary (Y-direction) of the domain contained in the StartFile to be selected (meters).",.1e6);
      RP::add("ElVentana.WindowZ_min", "Section boundary (Z-direction) of the domain contained in the StartFile to be selected (meters).",-1.8e7);
      RP::add("ElVentana.WindowZ_max", "Section boundary (Z-direction) of the domain contained in the StartFile to be selected (meters).",1.8e7);
      RP::add("ElVentana.dipoleScalingFactor","Scales the field strength of the magnetic dipole compared to Earths.", 1.0);
      RP::add("ElVentana.dipoleType","0: Normal 3D dipole, 1: line-dipole for 2D polar simulations", 0);
      RP::add("ElVentana.noDipoleInSW", "If set to 1, the dipole magnetic field is not set in the solar wind inflow cells. Default 0.", 0.0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         RP::add(pop + "_ElVentana.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
         RP::add(pop + "_ElVentana.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
         RP::add(pop + "_ElVentana.Temperatureratio", "Scale temperature from input values. Default 1. (Set to 1/4 for electrons)", 1.0);
      }
   }
   
   void ElVentana::getParameters(){
      Project::getParameters();
      
      int myRank;
      Real dummy;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      typedef Readparameters RP;
      if(!RP::get("ElVentana.constBgBX", this->constBgB[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.constBgBY", this->constBgB[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.constBgBZ", this->constBgB[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.StartFile", this->StartFile)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.WindowX_min", this->WindowX[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.WindowX_max", this->WindowX[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!RP::get("ElVentana.WindowY_min", this->WindowY[0])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
            exit(1);
         }
         if(!RP::get("ElVentana.WindowY_max", this->WindowY[1])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
            exit(1);
         }
         if(!RP::get("ElVentana.WindowZ_min", this->WindowZ[0])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
            exit(1);
         }
         if(!RP::get("ElVentana.WindowZ_max", this->WindowZ[1])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
            exit(1);
         }
         if(!RP::get("ElVentana.dipoleScalingFactor", this->dipoleScalingFactor)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
            exit(1);
         }
         if(!RP::get("ElVentana.dipoleType", this->dipoleType)) {
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
         if(!RP::get("ElVentana.noDipoleInSW", dummy)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
            exit(1);
         }
         this->noDipoleInSW = dummy == 1 ? true:false;


         // Per-population parameters
         for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
            const std::string& pop = getObjectWrapper().particleSpecies[i].name;
            ElVentanaSpeciesParameters sP;

            if(!RP::get(pop + "_ElVentana.nSpaceSamples", sP.nSpaceSamples)) {
               if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
               exit(1);
            }
            if(!RP::get(pop + "_ElVentana.nVelocitySamples", sP.nVelocitySamples)) {
               if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
               exit(1);
            }
            if(!RP::get(pop + "_ElVentana.Temperatureratio", sP.Temperatureratio)) {
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

            speciesParams.push_back(sP);
         }

      }
      
      bool ElVentana::initialize() {
         bool success = Project::initialize();

         return success;
      }

      bool ElVentana::rescalesDensity(const uint popID) const {
         // Rescale all population densities to ensure conservation of mass for accurate current
         return true;
      }

      Real ElVentana::getCorrectNumberDensity(spatial_cell::SpatialCell* cell,const uint popID) const {
         // This is the stored and read number density, not mass density!
         return cell->parameters[CellParams::RHOM];
      }


      Real ElVentana::getDistribValue(
         creal& x,creal& y, creal& z,
         creal& vx, creal& vy, creal& vz,
         creal& dvx, creal& dvy, creal& dvz,
         const uint popID
      ) const {

         const ElVentanaSpeciesParameters& sP = speciesParams[popID]; //Use parameters from cfg file in some cases?
         Real mass = getObjectWrapper().particleSpecies[popID].mass;
         Real temperature, initRho;
         Real radius;
         Real distvalue;
         CellID cellID = findCellIDXYZ(x, y, z);
         SpatialCell *cell = (*newmpiGrid)[cellID];

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
         
         // This is the stored and read number density, not mass density!
         initRho = cell->parameters[CellParams::RHOM];
         if(radius < this->ionosphereRadius) {
      // Just to be safe, there are observed cases where this failed.
      initRho = sP.ionosphereRho;
         }

         int myRank;
         MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

         // Defines the distribution function as a drifting maxwellian based on read values.
         // Assumes thus that the input pressure, bulk velocity, and density are from
         // a proton-only run and can be converted to electrons (by comparing with curl of B)
         if (initRho > 0.) {
            std::array<Real, 3> pressure_T;
            for(uint j=0;j<this->vecsizepressure;j++){
               pressure_T[j] = cell->parameters[CellParams::P_11+j];
            }
            temperature = (pressure_T[0] + pressure_T[1] + pressure_T[2]) / (3.*physicalconstants::K_B * initRho);
      // Scale temperatures from input values. For electrons, this should be about 1/4, for protons, 1
      temperature = temperature * sP.Temperatureratio;
            const std::array<Real, 3> v0 = this->getV0(x, y, z, popID)[0];
      if (mass > 0.5*physicalconstants::MASS_PROTON) {
         // protons, regular Maxwell-Boltzmann distribution
         distvalue = initRho * pow(mass / (2.0 * M_PI * physicalconstants::K_B * temperature), 1.5) *
         exp(- mass * ( pow(vx - v0[0], 2.0) + pow(vy - v0[1], 2.0) + pow(vz - v0[2], 2.0) ) /
            (2.0 * physicalconstants::K_B * temperature));
         return distvalue;
      } else {
         // electrons: assume that we have only protons and electrons as active populations
         std::array<Real, 3> Ji;
         Ji[0] = v0[0] * initRho * physicalconstants::CHARGE;
         Ji[1] = v0[1] * initRho * physicalconstants::CHARGE;
         Ji[2] = v0[2] * initRho * physicalconstants::CHARGE;
         // Now account for current requirement from curl of B
         const Real dBXdy = cell->derivativesBVOL[bvolderivatives::dPERBXVOLdy];
         const Real dBXdz = cell->derivativesBVOL[bvolderivatives::dPERBXVOLdz];
         const Real dBYdx = cell->derivativesBVOL[bvolderivatives::dPERBYVOLdx];
         const Real dBYdz = cell->derivativesBVOL[bvolderivatives::dPERBYVOLdz];
         const Real dBZdx = cell->derivativesBVOL[bvolderivatives::dPERBZVOLdx];
         const Real dBZdy = cell->derivativesBVOL[bvolderivatives::dPERBZVOLdy];
         std::array<Real, 3> Jreq;
         Jreq[0] = (dBYdz - dBZdy)/physicalconstants::MU_0;
         Jreq[1] = (dBZdx - dBXdz)/physicalconstants::MU_0;
         Jreq[2] = (dBXdy - dBYdx)/physicalconstants::MU_0;
         // Total required current density Jreq = Ji + Je
         // Electron velocity is Je/(density*charge)
         std::array<Real, 3> ve;
         ve[0] = (Jreq[0] - Ji[0])/(initRho*getObjectWrapper().particleSpecies[popID].charge);
         ve[1] = (Jreq[1] - Ji[1])/(initRho*getObjectWrapper().particleSpecies[popID].charge);
         ve[2] = (Jreq[2] - Ji[2])/(initRho*getObjectWrapper().particleSpecies[popID].charge);

         distvalue = initRho * pow(mass / (2.0 * M_PI * physicalconstants::K_B * temperature), 1.5) *
         exp(- mass * ( pow(vx - ve[0], 2.0) + pow(vy - ve[1], 2.0) + pow(vz - ve[2], 2.0) ) /
            (2.0 * physicalconstants::K_B * temperature));
         
         /* Try Maxwell-Jüttner relativistic distribution */
         /*
         Real theta = physicalconstants::K_B * temperature /(mass*physicalconstants::light*physicalconstants::light);
         ve[0] -= vx;
         ve[1] -= vy;
         ve[2] -= vz;
         Real vpc2 = (ve[0]*ve[0]+ve[1]*ve[1]+ve[2]*ve[2])/(physicalconstants::light*physicalconstants::light);
         Real gamma2 = 1.0/(1.-vpc2);
         Real gamma = sqrt(gamma2);
         Real beta = sqrt(vpc2);
         Real vv = sqrt(ve[0]*ve[0]+ve[1]*ve[1]+ve[2]*ve[2]);
         //Real K2 = boost::math::cyl_bessel_k((Real)2.0, (Real)(1./theta)); // Apparently need to use K2 i.e. cylindrical
         Real K2 = std::cyl_bessel_k(2.0, 1./theta); // Apparently need to use K2 i.e. cylindrical
         Real mom = gamma*mass*vv;
         Real mom_p = mass*(vv+0.5*dvx)/sqrt(1.-pow((vv+0.5*dvx)/physicalconstants::light,2));
         Real mom_m = mass*(vv-0.5*dvx)/sqrt(1.-pow((vv-0.5*dvx)/physicalconstants::light,2));
         //distvalue = initRho * (gamma2*beta/(theta*K2))*exp(-gamma/theta);
         distvalue = (initRho / (4.0*M_PI*pow(mass*physicalconstants::light,3.0)*theta*K2) )*exp(-gamma/theta);
         // convert from per momentum to per velocity
         distvalue *= (M_PI*4./3.)*(pow(mom_p,3)-pow(mom_m,3))/(pow(vv+0.5*dvx,3)-pow(vv-0.5*dvx,3));
         //if(myRank == MASTER_RANK) std::cerr << "Edist " << distvalue << std::endl;
         //if((cellID>1450)&&(cellID<1460)) std::cerr << "Edist " << distvalue << std::endl;
         */
         return distvalue;
      }
         } else {
            return 0;
         }

      }


      Real ElVentana::calcPhaseSpaceDensity(
         creal& x,creal& y,creal& z,
         creal& dx,creal& dy,creal& dz,
         creal& vx,creal& vy,creal& vz,
         creal& dvx,creal& dvy,creal& dvz,const uint popID
      ) const {
         const ElVentanaSpeciesParameters& sP = speciesParams[popID];
         creal d_x = dx / (2*sP.nSpaceSamples);
         creal d_y = dy / (2*sP.nSpaceSamples);
         creal d_z = dz / (2*sP.nSpaceSamples);
         creal d_vx = dvx / (2*sP.nVelocitySamples);
         creal d_vy = dvy / (2*sP.nVelocitySamples);
         creal d_vz = dvz / (2*sP.nVelocitySamples);
         Real avg = 0.0;         
         // #pragma omp parallel for collapse(6) reduction(+:avg)
         // WARNING No threading here if calling functions are already threaded
         for (uint i=0; i<sP.nSpaceSamples; ++i)
      for (uint j=0; j<sP.nSpaceSamples; ++j)
         for (uint k=0; k<sP.nSpaceSamples; ++k)
            for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
         for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
            for (uint vk=0; vk<sP.nVelocitySamples; ++vk) {
            avg += getDistribValue(x+(1+2*i)*d_x, y+(1+2*j)*d_y, z+(1+2*k)*d_z, vx+(1+2*vi)*d_vx, vy+(1+2*vj)*d_vy, vz+(1+2*vk)*d_vz, 2*d_vx, 2*d_vy, 2*d_vz, popID);
            }
         return avg /
            (sP.nSpaceSamples*sP.nSpaceSamples*sP.nSpaceSamples) /
            (sP.nVelocitySamples*sP.nVelocitySamples*sP.nVelocitySamples);
      }    

      vector<std::array<Real, 3>> ElVentana::getV0(
         creal x,
         creal y,
         creal z,
         const uint popID
      ) const {

         // Return the values calculated from the start value in each cell
         const ElVentanaSpeciesParameters& sP = this->speciesParams[popID];
         vector<std::array<Real, 3>> V0;
         std::array<Real, 3> v = {{0.0, 0.0, 0.0}};
         CellID cellID;
         
         cellID = findCellIDXYZ(x, y, z);
         SpatialCell *cell = (*newmpiGrid)[cellID];
         v[0] = cell->parameters[CellParams::VX];
         v[1] = cell->parameters[CellParams::VY];
         v[2] = cell->parameters[CellParams::VZ];
         
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
         
         if(radius < this->ionosphereRadius) {
      // Just to be safe, there are observed cases where this failed.
      for (uint64_t i=0; i<3; i++) {
         v[i] = ionosphereV0[i];
      }
         }
         
         // Check if velocity is within the velocity space boundaries
         if ( v[0] < getObjectWrapper().velocityMeshes[popID].meshMinLimits[0] ||
            v[0] > getObjectWrapper().velocityMeshes[popID].meshMaxLimits[0] ||
         v[1] < getObjectWrapper().velocityMeshes[popID].meshMinLimits[1] ||
            v[1] > getObjectWrapper().velocityMeshes[popID].meshMaxLimits[1] ||
            v[2] < getObjectWrapper().velocityMeshes[popID].meshMinLimits[2] ||
            v[2] > getObjectWrapper().velocityMeshes[popID].meshMaxLimits[2] ) {
            cerr << "ABORTING!!! Bulk velocity read from StartFile is outside the velocity space boundaries. " << endl;
            exit(1);
         }  

         V0.push_back(v);
         return V0;
      }

      /* Function to read relevant variables and store them to be read when each cell is being setup */
      void ElVentana::setupBeforeSetCell(const std::vector<CellID>& cells, 
                     dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                     bool& needCurl) {
         vector<CellID> fileCellsID; /*< CellIds for all cells in file*/
         int myRank,processes;
         const string filename = this->StartFile;
         int isbulk = 0;
         Real *buffer;
         uint64_t fileOffset, vecsize;
         vlsv::datatype::type dataType;
         //int k;

         // Read densities, velocities, pressures and magnetic fields from VLSV file
         // Attempt to open VLSV file for reading:
         MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
         MPI_Comm_size(MPI_COMM_WORLD,&processes);
         MPI_Info mpiInfo = MPI_INFO_NULL;
         std::array<double, 3> fileMin, fileMax, fileDx;
         std::array<uint64_t, 3> fileCells;

         if (this->vlsvSerialReader.open(filename) == false) {
            if(myRank == MASTER_RANK) 
               cout << "Could not open file: " << filename << endl;
            exit(1);
         }

         //if (this->vlsvParaReader.open(filename,MPI_COMM_WORLD,MASTER_RANK,mpiInfo) == false) {
         //   if (myRank == MASTER_RANK) 
         //      cout << "Could not open file: " << filename << endl;
         //   exit(1);
         //}

         if (readCellIds(this->vlsvParaReader,fileCellsID,MASTER_RANK,MPI_COMM_WORLD) == false) {
            if (myRank == MASTER_RANK)
               cout << "Could not read cell IDs." << endl;
            exit(1);
         }

         MPI_Bcast(&(fileCellsID[0]),fileCellsID.size(),MPI_UINT64_T,MASTER_RANK,MPI_COMM_WORLD);

         readGridSize(fileMin, fileMax, fileCells, fileDx);
         //this->vlsvParaReader.close();
         //MPI_Barrier(MPI_COMM_WORLD);
         // Closed later
         // :D

         // Check if cell size from file is the same as from new grid!! 

         // Check file type... couldn't get the VLSVreader interface to work yet
         /*      vlsvinterface::Reader interfacer;
            interfacer.open(filename);
            std::list<std::string> variableNames;
            std::string gridname("SpatialGrid");
            interfacer.getVariableNames(gridname,variableNames);
            interfacer.close();
            if(find(variableNames.begin(), variableNames.end(), std::string("moments"))!=variableNames.end()) {
            // Moments were found, i.e. it's a restart file
            isbulk = 0;
            } else {
            // It is a bulk file
            isbulk = 1;
            } */

         // TODO: pointless now that we can determine variable names?
         if (filename.find("bulk.") != string::npos) { 
            // It is a bulk file
            isbulk = 1; // Hard setting for now...
            // Could be determined by checking if a "moments" variable exists in the vlsv file?
         }

         std::string varname = "";
         
         for (uint64_t i=0; i<cells.size(); i++) {
            SpatialCell* cell = mpiGrid[cells[i]];
            // Calculate cellID in old grid
            CellID oldCellID = getOldCellID(cells[i], mpiGrid, fileCells, fileMin, fileDx);

            // Calculate fileoffset corresponding to old cellID
            auto cellIt = std::find(fileCellsID.begin(), fileCellsID.end(), oldCellID);
            if (cellIt == fileCellsID.end()) {
               cerr << "Could not find cell " << cells[i] << " old ID " << oldCellID << " Reflevel " << cell->parameters[CellParams::REFINEMENT_LEVEL] << " at " << cell->parameters[CellParams::XCRD] << " " << cell->parameters[CellParams::YCRD] << " " << cell->parameters[CellParams::ZCRD] << std::endl;
               exit(1);
            } else {
               fileOffset = cellIt - fileCellsID.begin();
            }

            // NOTE: This section assumes that the magnetic field values are saved on the spatial
            // (vlasov) grid, instead of FSgrid. FSgrid input reading isn't supported yet.
            // Also assumes the perturbed_B values have been saved to the bulk files. If
            // Perturbed_B is missing, it could perhaps be reconstructed from total B and
            // the analytically calculated background fields.
            //
            // The values are face-averages, not cell-averages, but are temporarily read into PERBXVOL anyway.
            varname = pickVarName("SpatialGrid", {"perturbed_B", "B"});
            if (varname == "") {
               perBSet = false;
            } else {
               buffer = readVar(varname, fileOffset, vecsize);
               for (uint j=0; j < vecsize; j++) {
                  mpiGrid[cells[i]]->parameters[CellParams::PERBXVOL+j] = buffer[j];
               }
               perBSet = true;

               totalBRead = (varname == "B");

               delete[] buffer;
            }

            // NOTE: This section assumes that the electric field values are saved on the spatial
            // (vlasov) grid, instead of FSgrid. FSgrid input reading isn't supported yet.
            // Also assumes the E values have been saved to the bulk files.
            //
            // The values are edge-averages, not cell-averages, but are temporarily read into EXVOL anyway.
            //attribs.push_back(make_pair("mesh","SpatialGrid"));
            //attribs.push_back(make_pair("name","E"));
            //if (this->vlsvSerialReader.getArrayInfo("VARIABLE",attribs,arraySize,this->vecsizeE,dataType,byteSize) == false) {
            //   if(myRank == MASTER_RANK) logFile << "(START)  ERROR: Failed to read E array info" << endl << write;
            //   exit(1);
            //}
            //buffer=new Real[this->vecsizeE];
            //if (this->vlsvSerialReader.readArray("VARIABLE", attribs, fileOffset, 1, (char *)buffer) == false ) {
            //   if(myRank == MASTER_RANK) logFile << "(START)  ERROR: Failed to read E"  << endl << write;
            //   exit(1);
            //}
            //for (uint j=0; j<vecsizeE; j++) {
            //   mpiGrid[cells[i]]->parameters[CellParams::EXVOL+j] = buffer[j];
            //}
            //delete[] buffer;
            //attribs.pop_back();
            //attribs.pop_back();
               
            // The background fields are initialized directly on FSgrid and are not read in.
            // The commented out block literally reads them in. What the fuck?

            //attribs.push_back(make_pair("mesh","SpatialGrid"));
            //varname = pickVarName({"background_B", "B"});
            //if (varname == "") {
            //   if (myRank == MASTER_RANK) 
            //      logFile << "(START)  ERROR: No background B variable found!" << endl << write;
            //   exit(1);
            //}
            //attribs.push_back(make_pair("name", varname));

            //if (this->vlsvSerialReader.getArrayInfo("VARIABLE",attribs,arraySize,this->vecsizebackground_B,dataType,byteSize) == false) {
            //   logFile << "(START)  ERROR: Failed to read " << varname << " array info" << endl << write; 
            //   exit(1);
            //}
            //buffer=new Real[this->vecsizebackground_B];
            //if (this->vlsvSerialReader.readArray("VARIABLE", attribs, fileOffset, 1, (char *)buffer) == false ) {
            //   logFile << "(START)  ERROR: Failed to read " << varname << endl << write;
            //   exit(1);
            //}
            //if (varname == "B") {
            //   for (uint j=0; j<vecsizebackground_B; j++) {
            //      mpiGrid[cells[i]]->parameters[CellParams::BGBXVOL+j] = buffer[j] - mpiGrid[cells[i]]->parameters[CellParams::PERBXVOL+j]; 
            //   }
            //} else {
            //   for (uint j=0; j<vecsizebackground_B; j++) {
            //      mpiGrid[cells[i]]->parameters[CellParams::BGBXVOL+j] = buffer[j];
            //   }
            //}
            //delete[] buffer;
            //attribs.pop_back();
            //attribs.pop_back();

            // The following sections are not multipop-safe. Multipop-handling can probably be added via the
            // variableNames list (see detection of bulk / restart files)
            // if (isbulk == 1) {
            //    attribs.push_back(make_pair("name","vg_rho"));
            // } else {
            //    attribs.push_back(make_pair("name","moments"));
            // }

            varname = pickVarName("SpatialGrid", {"proton/vg_rho", "rho", "moments"});
            if (varname == "") {
               if (myRank == MASTER_RANK) 
                  logFile << "(START)  ERROR: No rho variable found!" << endl << write;
               exit(1);
            }

            buffer = readVar(varname, fileOffset, vecsize);

            // Parse zeroth and first moments data
            if (vecsize == 5) {
               // Reading the new format (multipop) restart file (rhom, massVx, massVy, massVz, rhoq)
               mpiGrid[cells[i]]->parameters[CellParams::RHOM] = buffer[4] / physicalconstants::CHARGE;
               mpiGrid[cells[i]]->parameters[CellParams::VX] = buffer[1]/buffer[0];
               mpiGrid[cells[i]]->parameters[CellParams::VY] = buffer[2]/buffer[0];
               mpiGrid[cells[i]]->parameters[CellParams::VZ] = buffer[3]/buffer[0];
            } else if (vecsize == 4) {
               // Reading the old format restart file (rho, Vx, Vy, Vz)
               for (uint j=0; j < vecsize; j++) {
                  mpiGrid[cells[i]]->parameters[CellParams::RHOM+j] = buffer[j];
               }	   
            } else if (vecsize == 1)  {
               // Reading a bulk file. First store number density.
               mpiGrid[cells[i]]->parameters[CellParams::RHOM] = buffer[0];

               // Now read rho_v in a separate read call	    
               varname = pickVarName("SpatialGrid", {"proton/vg_v", "rho_v"});
               if (varname == "") {
                  if (myRank == MASTER_RANK) 
                     logFile << "(START)  ERROR: No v variable found!" << endl << write;
                  exit(1);
               }

               buffer = readVar(varname, fileOffset, vecsize);

               if(varname=="rho_v") {
                  for (uint j=0; j < vecsize; j++) {
                     mpiGrid[cells[i]]->parameters[CellParams::VX+j] = buffer[j]/mpiGrid[cells[i]]->parameters[CellParams::RHOM];
                  }
               } else {
                  for (uint j=0; j < vecsize; j++) {
                     mpiGrid[cells[i]]->parameters[CellParams::VX+j] = buffer[j];
                  }
               }
            } else {
               if(myRank == MASTER_RANK) cout << "Could not identify restart file format for file: " << this->StartFile << endl;
               exit(1);
            }
            delete[] buffer;

            // Read second velocity moments data (pressure)
            // if (isbulk == 1) {
            //    attribs.push_back(make_pair("name","vg_ptensor_diagonal"));
            // } else {
            //    attribs.push_back(make_pair("name","pressure"));
            // }
            varname = pickVarName("SpatialGrid", {"proton/vg_ptensor_diagonal", "PTensorDiagonal", "pressure"});
            if (varname == "") {
               if (myRank == MASTER_RANK) 
                  logFile << "(START)  ERROR: No PTensorDiagonal variable found!" << endl << write;
               exit(1);
            }

            buffer = readVar(varname, fileOffset, this->vecsizepressure);
            for (uint j=0; j < this->vecsizepressure; j++) {
               mpiGrid[cells[i]]->parameters[CellParams::P_11+j] = buffer[j];
            }
            delete[] buffer;
         }

         newmpiGrid = &mpiGrid;
         this->vlsvSerialReader.close();

         // Let initializeGrid know that this project needs the Curl of B
         needCurl = true;
      }

      void ElVentana::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
         // EJE initialization to zero shouldn't be tied to this project so isn't done here anymore.
      }


      /* set 0-centered dipole */
      void ElVentana::setProjectBField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
         FsGrid< fsgrids::technical, 2>& technicalGrid
   ) {

     Dipole bgFieldDipole;
     LineDipole bgFieldLineDipole;

     // The hardcoded constants of dipole and line dipole moments are obtained
     // from Daldorff et al (2014), see
     // https://github.com/fmihpc/vlasiator/issues/20 for a derivation of the
     // values used here.
     switch(this->dipoleType) {
      case 0:
         bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
         setBackgroundField(bgFieldDipole,BgBGrid);
         break;
      case 1:
         bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
         setBackgroundField(bgFieldLineDipole,BgBGrid);
         break;
      case 2:
         bgFieldLineDipole.initialize(126.2e6 *this->dipoleScalingFactor, 0.0, 0.0, 0.0 );//set dipole moment     
         setBackgroundField(bgFieldLineDipole,BgBGrid);
         // Ignore mirror dipole for ElVentana runs. Could be added here if needed.
         break;
      case 3:
      case 4: // Vector potential dipole stored in perturbed field, which is read later
         bgFieldDipole.initialize(8e15 *this->dipoleScalingFactor, 0.0, 0.0, 0.0, 0.0 );//set dipole moment
         setBackgroundField(bgFieldDipole,BgBGrid);
         // Ignore mirror dipole for ElVentana runs. Could be added here if needed.
         break;
      default:
	    setBackgroundFieldToZero(BgBGrid);
	    
     }

      const auto localSize = BgBGrid.getLocalSize().data();
      
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

      // Try reading perturbed B-field here if we're not in restart
      if (!P::isRestart) {
         std::string varName = pickVarName("fsgrid", {"fg_b_perturbed", "fg_b"});
         if (varName == "" || !readFsGridVariable(varName, perBGrid)) {            
            if (!perBSet) {
               logFile << "(START)  ERROR: No B field in file!" << endl << write;
               exit(1);
            } else {
               logFile << "(START)  No B field in FsGrid, using volumetric averages" << endl << write;
            }
         } else {
            totalBRead = (varName == "fg_b");
            std::cerr << "B Read!" << std::endl;
         }

         if (totalBRead) {
            logFile << "(START)  No b_perturbed in FsGrid, calculating from b!" << endl << write;
            #pragma omp for collapse(3)
            for (int x = 0; x < localSize[0]; ++x) {
               for (int y = 0; y < localSize[1]; ++y) {
                  for (int z = 0; z < localSize[2]; ++z) {
                     std::array<Real, fsgrids::bfield::N_BFIELD>* bcell = perBGrid.get(x, y, z);
                     std::array<Real, fsgrids::bgbfield::N_BGB>* bgcell = BgBGrid.get(x, y, z);
                     bcell->at(fsgrids::bfield::PERBX) -= bgcell->at(fsgrids::bgbfield::BGBX);
                     bcell->at(fsgrids::bfield::PERBY) -= bgcell->at(fsgrids::bgbfield::BGBY);
                     bcell->at(fsgrids::bfield::PERBZ) -= bgcell->at(fsgrids::bgbfield::BGBZ);
                  }
               }
            }
         }
      }
   }

   CellID ElVentana::findCellIDXYZ(creal x, creal y, creal z) const {
      return newmpiGrid->get_existing_cell({x, y, z});
   }

   // Reads physical minima and maxima, amount of cells and their dimensions
   bool ElVentana::readGridSize(std::array<double, 3> &fileMin, std::array<double, 3> &fileMax, std::array<uint64_t, 3> &fileCells, std::array<double, 3> &fileDx) {
      int myRank = 0;
      //MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      double filexmin, fileymin, filezmin, filexmax, fileymax, filezmax;
      uint filexcells, fileycells, filezcells;

      if (myRank == MASTER_RANK) 
         cout << "Trying to read parameters from file... " << endl;

      // TODO looks like hot garbage, rewrite
      if (this->vlsvParaReader.readParameter("xmin", filexmin) == false) {
         if (myRank == MASTER_RANK)
            cout << " Could not read parameter xmin. " << endl;
         return false;
      }
      if (this->vlsvParaReader.readParameter("ymin", fileymin) == false) {
         if (myRank == MASTER_RANK)
            cout << " Could not read parameter ymin. " << endl;
         return false;
      }
      if (this->vlsvParaReader.readParameter("zmin", filezmin) == false) {
         if (myRank == MASTER_RANK) 
            cout << " Could not read parameter zmin. " << endl;
         return false;
      }
      if (this->vlsvParaReader.readParameter("xmax", filexmax) == false) {
         if (myRank == MASTER_RANK)
            cout << " Could not read parameter xmax. " << endl;
         return false;
      }
      if (this->vlsvParaReader.readParameter("ymax", fileymax) == false) {
         if (myRank == MASTER_RANK)
            cout << " Could not read parameter ymax. " << endl;
         return false;
      }
      if (this->vlsvParaReader.readParameter("zmax", filezmax) == false) {
         if (myRank == MASTER_RANK) 
            cout << " Could not read parameter zmax. " << endl;
         return false;
      }
      if (this->vlsvParaReader.readParameter("xcells_ini", filexcells) == false) {
         if (myRank == MASTER_RANK)
            cout << " Could not read parameter xcells_ini. " << endl;
         return false;
      }
      if (this->vlsvParaReader.readParameter("ycells_ini", fileycells) == false) {
         if (myRank == MASTER_RANK)
            cout << " Could not read parameter ycells_ini. " << endl;
         return false;
      }
      if (this->vlsvParaReader.readParameter("zcells_ini", filezcells) == false) {
         if (myRank == MASTER_RANK) 
            cout << " Could not read parameter zcells_ini. " << endl;
         return false;
      }
      if (myRank == MASTER_RANK) 
         cout << "All parameters read." << endl;

      fileMin = {filexmin, fileymin, filezmin};
      fileMax = {filexmax, fileymax, filezmax};
      fileCells = {filexcells, fileycells, filezcells};
      for (int i = 0; i < 3; ++i) {
         fileDx[i] = (fileMax[i] - fileMin[i])/fileCells[i];
      }

      return true;
   }

   // Read variable from FsGrid to the ElVentana window
   template<unsigned long int N> bool ElVentana::readFsGridVariable(const string& variableName, FsGrid<std::array<Real, N>,2>& targetGrid) {

      uint64_t arraySize;
      uint64_t vectorSize;
      vlsv::datatype::type dataType;
      uint64_t byteSize;
      list<pair<string,string> > attribs;
      const string filename = this->StartFile;

      int numWritingRanks = 0;

      // Are we restarting from the same number of tasks, or a different number?
      int size, myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      MPI_Comm_size(MPI_COMM_WORLD,&size);
      MPI_Info mpiInfo = MPI_INFO_NULL;
      if (myRank == MASTER_RANK)
         std::cerr << "Reading " << variableName << "from FsGrid!" << std::endl;

      //if (this->vlsvParaReader.open(filename,MPI_COMM_WORLD,MASTER_RANK,mpiInfo) == false) {
      //   if (myRank == MASTER_RANK) 
      //      cout << "Could not open file: " << filename << endl;
      //   return false;
      //}

      if(this->vlsvParaReader.readParameter("numWritingRanks",numWritingRanks) == false) {
         std::cerr << "FSGrid writing rank number not found";
         return false;
      }
      
      attribs.push_back(make_pair("name",variableName));
      attribs.push_back(make_pair("mesh","fsgrid"));

      if (this->vlsvParaReader.getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,dataType,byteSize) == false) {
         logFile << "(RESTART)  ERROR: Failed to read array info for " << variableName << endl << write;
         return false;
      }

      std::array<double, 3> fileMin, fileMax, fileDx;
      std::array<uint64_t, 3> fileCells;
      if(!readGridSize(fileMin, fileMax, fileCells, fileDx)) {
         return false;
      }

      for (int i = 0; i < P::amrMaxSpatialRefLevel; ++i) {
         for (auto it = fileDx.begin(); it < fileDx.end(); ++it)
            *it /= 2;
      }

      std::array<int32_t,3>& localSize = targetGrid.getLocalSize();
      std::array<int32_t,3>& localStart = targetGrid.getLocalStart();
      std::array<int32_t,3> localOffset;
      for (int i = 0; i < 3; ++i) {
         localOffset[i] = (targetGrid.physicalGlobalStart[i] - fileMin[i]) / fileDx[i];
      }

      // Determine our tasks storage size
      size_t storageSize = localSize[0]*localSize[1]*localSize[2];

      // More difficult case: different number of tasks.
      // In this case, our own fsgrid domain overlaps (potentially many) domains in the file.
      // We read the whole source rank into a temporary buffer, and transfer the overlapping
      // part.
      //
      // +------------+----------------+
      // |            |                |
      // |    . . . . . . . . . . . .  |
      // |    .<----->|<----------->.  |
      // |    .<----->|<----------->.  |
      // |    .<----->|<----------->.  |
      // +----+-------+-------------+--|
      // |    .<----->|<----------->.  |
      // |    .<----->|<----------->.  |
      // |    .<----->|<----------->.  |
      // |    . . . . . . . . . . . .  |
      // |            |                |
      // +------------+----------------+

      // Determine the decomposition in the file and the one in RAM for our restart
      std::array<int,3> fileDecomposition;
      // WHY does this want int instead of uint_64t (CellID)?
      std::array<int, 3> fileCellsCopy;
      for (int i = 0; i < 3; ++i) {
         fileCells[i] *= pow(2, P::amrMaxSpatialRefLevel);
         fileCellsCopy[i] = fileCells[i];
      }
 
      targetGrid.computeDomainDecomposition(fileCellsCopy, numWritingRanks, fileDecomposition);

      // Iterate through tasks and find their overlap with our domain.
      size_t fileOffset = 0;
      for(int task = 0; task < numWritingRanks; task++) {
         std::array<int32_t,3> thatTasksSize;
         std::array<int32_t,3> thatTasksStart;
         thatTasksSize[0] = targetGrid.calcLocalSize(fileCells[0], fileDecomposition[0], task/fileDecomposition[2]/fileDecomposition[1]);
         thatTasksSize[1] = targetGrid.calcLocalSize(fileCells[1], fileDecomposition[1], (task/fileDecomposition[2])%fileDecomposition[1]);
         thatTasksSize[2] = targetGrid.calcLocalSize(fileCells[2], fileDecomposition[2], task%fileDecomposition[2]);

         thatTasksStart[0] = targetGrid.calcLocalStart(fileCells[0], fileDecomposition[0], task/fileDecomposition[2]/fileDecomposition[1]);
         thatTasksStart[1] = targetGrid.calcLocalStart(fileCells[1], fileDecomposition[1], (task/fileDecomposition[2])%fileDecomposition[1]);
         thatTasksStart[2] = targetGrid.calcLocalStart(fileCells[2], fileDecomposition[2], task%fileDecomposition[2]);

         // Iterate through overlap area
         std::array<int,3> overlapStart,overlapEnd,overlapSize;
         for (int i = 0; i < 3; ++i) {
            overlapStart[i] = max(localStart[i] + localOffset[i], thatTasksStart[i]);
            overlapEnd[i] = min(localStart[i] + localOffset[i] + localSize[i], thatTasksStart[i] + thatTasksSize[i]);
            overlapSize[i] = max(overlapEnd[i] - overlapStart[i], 0);
         }

         // Read into buffer
         std::vector<Real> buffer(thatTasksSize[0]*thatTasksSize[1]*thatTasksSize[2]*N);

         phiprof::start("readArray");
         // TODO: Should these be multireads instead? And/or can this be parallelized?
         if(this->vlsvParaReader.readArray("VARIABLE",attribs, fileOffset, thatTasksSize[0]*thatTasksSize[1]*thatTasksSize[2], buffer.data()) == false) {
            logFile << "(START)  ERROR: Failed to read fsgrid variable " << variableName << endl << write;
            return false;
         }
         phiprof::stop("readArray");

         phiprof::start("memcpy");

         // Copy continuous stripes in x direction.
         for(int z=overlapStart[2]; z<overlapEnd[2]; z++) {
            for(int y=overlapStart[1]; y<overlapEnd[1]; y++) {
               for(int x=overlapStart[0]; x<overlapEnd[0]; x++) {
                  int index = (z - thatTasksStart[2]) * thatTasksSize[0]*thatTasksSize[1]
                     + (y - thatTasksStart[1]) * thatTasksSize[0]
                     + (x - thatTasksStart[0]);

                  memcpy(targetGrid.get(x - localStart[0] - localOffset[0], y - localStart[1] - localOffset[1], z - localStart[2] - localOffset[2]), &buffer[index*N], N*sizeof(Real));
               }
            }
         }
         phiprof::stop("memcpy");
         fileOffset += thatTasksSize[0] * thatTasksSize[1] * thatTasksSize[2];
      } 
      phiprof::start("updateGhostCells");
      targetGrid.updateGhostCells();
      phiprof::stop("updateGhostCells");

      this->vlsvParaReader.close();

      return true;
   }

   // Pick first variable found in file from list of names given
   std::string ElVentana::pickVarName(const std::string &grid, const std::list<std::string> &varNames) {
      vlsvinterface::Reader r;
      r.open(StartFile);

      std::list<std::string> fileVarNames;
      r.getVariableNames(grid, fileVarNames);

      std::string varName = "";
      for (std::string s : varNames) {
         if (find(fileVarNames.begin(), fileVarNames.end(), s) != fileVarNames.end()) {
            varName = s;
            break;
         }
      }

      r.close();
      return varName;
   }

   // Read variable in single cell with given offset. Returns pointer to buffer size vecsize containing the variable
   Real* ElVentana::readVar(string varname, CellID fileOffset, uint64_t &vecsize) {
      list<pair<string,string>> attribs;
      uint64_t arraySize, byteSize;
      vlsv::datatype::type dataType;

      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

      attribs.push_back(make_pair("mesh","SpatialGrid"));
      attribs.push_back(make_pair("name", varname));
      if (this->vlsvSerialReader.getArrayInfo("VARIABLE",attribs,arraySize,vecsize,dataType,byteSize) == false) {
         if(myRank == MASTER_RANK) logFile << "(START)  ERROR: Failed to read " << varname << " array info" << endl << write;
         exit(1);
      }

      Real* buffer = new Real[vecsize];

      if (this->vlsvSerialReader.readArray("VARIABLE", attribs, fileOffset, 1, buffer) == false ) {
         if(myRank == MASTER_RANK) logFile << "(START)  ERROR: Failed to read " << varname << endl << write;
         exit(1);
      }

      return buffer;
   }

   // Returrns CellID in file grid corresponding to newID in the window
   CellID ElVentana::getOldCellID(CellID newID, dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, std::array<CellID, 3> fileCells, std::array<double, 3> &fileMin, std::array<double, 3> fileDx) {
      int refLevel = mpiGrid.get_refinement_level(newID);
      dccrg::Mapping fileMapping;
      if (!(fileMapping.set_length(fileCells) && fileMapping.set_maximum_refinement_level(P::amrMaxSpatialRefLevel))) {
         std::cerr << "Could not set file mapping!" << std::endl;
      }

      std::array<double,3> xyz = mpiGrid.get_center(newID);      

      std::array<uint64_t, 3> indices;
      for (int i = 0; i < 3; ++i) {
         indices[i] = (xyz[i] - fileMin[i]) / (fileDx[i] / pow(2, P::amrMaxSpatialRefLevel));
      }

      return fileMapping.get_cell_from_indices(indices, refLevel);
   }

   // Refine spatial cells to same level as file
   bool ElVentana::refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      vector<CellID> fileCellsID; /*< CellIds for all cells in file*/
      std::array<double, 3> fileMin, fileMax, fileDx;
      std::array<CellID, 3> fileCells;
      const string filename = this->StartFile;
      std::vector<CellID> cells = getLocalCells();

      int myRank, processes;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      MPI_Comm_size(MPI_COMM_WORLD,&processes);
      MPI_Info mpiInfo = MPI_INFO_NULL;

      if (this->vlsvParaReader.open(filename,MPI_COMM_WORLD,MASTER_RANK,mpiInfo) == false) {
         if (myRank == MASTER_RANK)
            cout << "Could not open file: " << filename << endl;
         exit(1);
      }

      if (readCellIds(this->vlsvParaReader,fileCellsID,MASTER_RANK,MPI_COMM_WORLD) == false) {
         if (myRank == MASTER_RANK)
            cout << "Could not read cell IDs." << endl;
         exit(1);
      }

      readGridSize(fileMin, fileMax, fileCells, fileDx);

      // Refine all cells that aren't found in file.
      for (uint64_t refLevel=0; refLevel < P::amrMaxSpatialRefLevel; ++refLevel) {
         int oldCells = cells.size();
         for (int i = 0; i < cells.size(); ++i) {
            CellID id = cells[i];
            //if (myRank == MASTER_RANK)
            //   std::cout << "ID: " << id;
            CellID oldCellID = getOldCellID(id, mpiGrid, fileCells, fileMin, fileDx);
            if (std::find(fileCellsID.begin(), fileCellsID.end(), oldCellID) == fileCellsID.end()) {
               //std::cerr << "Refined " << id << " old ID " << oldCellID << std::endl;
               mpiGrid.refine_completely(id);
            } else {
               //cerr << "Found " << id << " old ID " << oldCellID << " Reflevel " << refLevel << std::endl;
               //std::cerr << "Found " << id << " old ID " << oldCellID << std::endl;
            }
         }
         //bool done = false;
         //done = mpiGrid.stop_refining().empty();
         cells = mpiGrid.stop_refining(true);
         std::cerr << "Refined " << cells.size() << " cells out of " << oldCells << std::endl;
         //mpiGrid.balance_load();
         //recalculateLocalCellsCache();
         //initSpatialCellCoordinates(mpiGrid);
         //setFaceNeighborRanks(mpiGrid);
         //if (done)
         //   break;
      }

      //this->vlsvParaReader.close();
      return true;
   }

} // namespace projects

