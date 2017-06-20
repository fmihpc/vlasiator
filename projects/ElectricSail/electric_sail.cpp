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
 * 
 * File:   electric_sail.cpp
 * Author: sandroos
 *
 * Created on March 3, 2015
 */

#include <cstdlib>
#include <iostream>

#include "../../readparameters.h"
#include "../../logger.h"
#include "../../object_wrapper.h"
#include "../../poisson_solver/poisson_solver.h"
#include "../read_gaussian_population.h"

#include "electric_sail.h"

using namespace std;
using namespace spatial_cell;

extern Logger logFile;

namespace projects {
   
   static Real radius = 0;
   
   Population::Population(const Real& rho,const Real& Tx,const Real& Ty,const Real& Tz,
                          const Real& Vx,const Real& Vy,const Real& Vz) {
      this->rho = rho;
      T[0] = Tx;
      T[1] = Ty;
      T[2] = Tz;
      V[0] = Vx;
      V[1] = Vy;
      V[2] = Vz;
   }
   
   ElectricSail::ElectricSail(): TriAxisSearch() { }
   
   ElectricSail::~ElectricSail() { }
   
   void ElectricSail::addParameters() {
      typedef Readparameters RP;
      RP::add("ElectricSail.solver","Name of the Poisson solver",string("SOR"));
      RP::add("ElectricSail.radius","Radius where charge density is non-zero",(Real)15e3);
      RP::add("ElectricSail.max_iterations","Maximum number of iterations",(uint)1000);
      RP::add("ElectricSail.min_relative_change","Potential is iterated until it the relative change is less than this value",(Real)1e-5);
      RP::add("ElectricSail.clear_potential","Clear potential each timestep before solving Poisson?",true);
      RP::add("ElectricSail.is_2D","If true then system is two-dimensional in xy-plane",true);
      RP::add("ElectricSail.tether_x","Electric sail tether x-position",(Real)0.0);
      RP::add("ElectricSail.tether_y","Electric sail tether y-position",(Real)0.0);
      RP::add("ElectricSail.max_absolute_error","Maximum absolute error allowed in Poisson solution",(Real)1e-4);
      RP::add("ElectricSail.add_particle_cloud","If true, add charge neutralizing particle cloud around tethet (bool)",false);
      RP::add("ElectricSail.tetherCharge","Tether charge per meter in elementary charges",(Real)200e9);
      RP::add("ElectricSail.timeDependentCharge","If true, tether charge is time dependent (bool)",false);
      RP::add("ElectricSail.tetherChargeRiseTime","Time when tether charge reaches its maximum value",0.0);
      RP::add("ElectricSail.useBackgroundField","Use background field instead of external charge density?",false);

      projects::ReadGaussianPopulation rgp;
      rgp.addParameters("ElectricSail");
   }

   Real ElectricSail::getCorrectNumberDensity(spatial_cell::SpatialCell* cell,const uint popID) const {
      if (addParticleCloud == false) return populations[popID].rho;
      if (getObjectWrapper().particleSpecies[popID].name != "Electron") return populations[popID].rho;

      if (tetherUnitCharge < 0) {
         cerr << "negative tether not implemented in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }

      const Real* parameters = cell->get_cell_parameters();

      Real pos[3];
      pos[0] = parameters[CellParams::XCRD] + 0.5*parameters[CellParams::DX];
      pos[1] = parameters[CellParams::YCRD] + 0.5*parameters[CellParams::DY];
      pos[2] = parameters[CellParams::ZCRD] + 0.5*parameters[CellParams::DZ];
      Real radius2 = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];

      if (radius2 > particleCloudRadius*particleCloudRadius) return populations[popID].rho;

      const Real charge = getObjectWrapper().particleSpecies[popID].charge;
      const Real DZ = parameters[CellParams::DZ];
      Real cloudDens = -tetherUnitCharge / (M_PI*particleCloudRadius*particleCloudRadius*charge);
      return populations[popID].rho + cloudDens;
   }

   Real ElectricSail::getDistribValue(creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz,const uint popID) const {
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;
      const Population& pop = populations[popID];

      Real value 
         = pop.rho
         * mass / (2.0 * M_PI * kb )
         * exp(-mass * (  pow(vx-pop.V[0],2.0) / (2.0*kb*pop.T[0]) 
                        + pow(vy-pop.V[1],2.0) / (2.0*kb*pop.T[1])
                       )
              );

      return value;
   }

   void ElectricSail::getParameters() {
      bool success = true;
      
      Project::getParameters();
      typedef Readparameters RP;
      RP::get("ElectricSail.solver",poisson::Poisson::solverName);
      RP::get("ElectricSail.radius",radius);
      RP::get("ElectricSail.max_iterations",poisson::Poisson::maxIterations);
      RP::get("ElectricSail.min_relative_change",poisson::Poisson::minRelativePotentialChange);
      RP::get("ElectricSail.is_2D",poisson::Poisson::is2D);
      RP::get("ElectricSail.clear_potential",poisson::Poisson::clearPotential);
      RP::get("ElectricSail.tether_x",tether_x);
      RP::get("ElectricSail.tether_y",tether_y);
      RP::get("ElectricSail.max_absolute_error",poisson::Poisson::maxAbsoluteError);
      RP::get("ElectricSail.add_particle_cloud",addParticleCloud);
      RP::get("ElectricSail.tetherCharge",tetherUnitCharge);
      RP::get("ElectricSail.timeDependentCharge",timeDependentCharge);
      RP::get("ElectricSail.tetherChargeRiseTime",tetherChargeRiseTime);
      RP::get("ElectricSail.useBackgroundField",useBackgroundField);
      
      projects::ReadGaussianPopulation rgp;
      projects::GaussianPopulation gaussPops;
      if (rgp.getParameters("ElectricSail",gaussPops) == false) success = false;
      
      if (success == false) {
         stringstream ss;
         ss << "(ElectricSail) ERROR: Could not find parameters for all species.";
         ss << "\t CHECK YOUR CONFIG FILE!" << endl;
         cerr << ss.str();
         exit(1);
      }

      for (size_t i=0; i<gaussPops.rho.size(); ++i) {
         populations.push_back(projects::Population(gaussPops.rho[i],gaussPops.Tx[i],gaussPops.Ty[i],
                                                    gaussPops.Tz[i],gaussPops.Vx[i],gaussPops.Vy[i],gaussPops.Vz[i]));
      }
      
      particleCloudRadius = 50.0;
      tetherUnitCharge *= physicalconstants::CHARGE*Parameters::dz_ini;
      if (timeDependentCharge == true) poisson::Poisson::timeDependentBackground = true;
   }

   bool ElectricSail::initialize() {
      bool success = Project::initialize();
      if (populations.size() < getObjectWrapper().particleSpecies.size()) {
         cerr << "(ElectricSail) ERROR: you have not defined parameters for all populations in ";
         cerr << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      
      logFile << "(ElectricSail) Population parameters are:" << endl;
      for (size_t i=0; i<populations.size(); ++i) {
         logFile << "\t Population: " << i << endl;
         logFile << "\t rho       : " << populations[i].rho << endl;
         logFile << "\t T         : " << populations[i].T[0] << '\t' << populations[i].T[1] << '\t' << populations[i].T[2] << endl;
         logFile << "\t V         : " << populations[i].V[0] << '\t' << populations[i].V[1] << '\t' << populations[i].V[2] << endl;
      }
      logFile << endl;
      
      logFile << "(ElectricSail) Tether parameters are:" << endl;
      logFile << "\t add prtcl cloud?: ";
      if (addParticleCloud == true) logFile << "Yes" << endl;
      else logFile << "No" << endl;
      logFile << "\t charge per meter: " << tetherUnitCharge/Parameters::dz_ini << endl;
      logFile << "\t charge (total)  : " << tetherUnitCharge << endl;
      logFile << "\t charge time-dep?: ";
      if (timeDependentCharge == true) {
         logFile << "Yes" << endl;
         logFile << "\t charge rise time: " << tetherChargeRiseTime << endl;
      } else logFile << "No" << endl;
      logFile << writeVerbose;

      return success;
   }

   bool ElectricSail::rescalesDensity(const uint popID) const {
      return true;
   }

   /**
    * 
    * NOTE: This is only called in grid.cpp:initializeGrid.
    * NOTE: This function must be thread-safe.
    */
   void ElectricSail::setCellBackgroundField(SpatialCell* cell) const {
      Real X  = cell->parameters[CellParams::XCRD];
      Real Y  = cell->parameters[CellParams::YCRD];
      Real Z  = cell->parameters[CellParams::ZCRD];
      Real DX = cell->parameters[CellParams::DX];
      Real DY = cell->parameters[CellParams::DY];
      Real DZ = cell->parameters[CellParams::DZ];

      cell->parameters[CellParams::RHOQ_EXT] = 0;
      Real pos[3];
      pos[0] = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      pos[1] = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      pos[2] = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      Real factor = 1.0;
      if (timeDependentCharge == true) {
         factor = max((Real)0.0,(Real)1.0+(Parameters::t-tetherChargeRiseTime)/tetherChargeRiseTime);
         factor = min((Real)1.0,factor);
      }

      if (useBackgroundField == false) {
         Real rad = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
         Real D3 = cell->parameters[CellParams::DX]*cell->parameters[CellParams::DY];
         if (rad <= 5) cell->parameters[CellParams::RHOQ_EXT] = 0.25*factor*tetherUnitCharge/D3/physicalconstants::EPS_0;
         
         cell->parameters[CellParams::BGEXVOL] = 0;
         cell->parameters[CellParams::BGEYVOL] = 0;
         cell->parameters[CellParams::BGEZVOL] = 0;
         return;
      }
       
      cell->parameters[CellParams::RHOQ_EXT] = 0;

      const Real EPSILON = 1e-30;
      uint N = 1;
      int N3_sum = 0;
      Real E_vol[3] = {0,0,0};
      
      bool ok = false;
      do {
         Real E_current[3] = {0,0,0};

         const Real DX_N = DX / N;
         const Real DY_N = DY / N;
         const Real DZ_N = DZ / N;

         // Sample E using N points
         Real E_dummy[3] = {0,0,0};
         for (uint k=0; k<N; ++k) for (uint j=0; j<N; ++j) for (uint i=0; i<N; ++i) {
            Real x[3];
            x[0] = X + 0.5*DX_N + i*DX_N;
            x[1] = Y + 0.5*DY_N + j*DY_N;
            x[2] = Z + 0.5*DZ_N + k*DZ_N;

            tetherElectricField(x,E_dummy);
            E_current[0] += E_dummy[0];
            E_current[1] += E_dummy[1];
            E_current[1] += E_dummy[2];
         }

         // Check if the current estimate of volume-averaged E is good enough
         Real delta = 0;
         delta = max(delta,(E_current[0]-E_vol[0])/(E_current[0]+EPSILON));
         delta = max(delta,(E_current[1]-E_vol[1])/(E_current[1]+EPSILON));
         delta = max(delta,(E_current[2]-E_vol[2])/(E_current[2]+EPSILON));
         if (delta < poisson::Poisson::minRelativePotentialChange) ok = true;
         if (N >= poisson::Poisson::maxIterations) ok = true;

         // Add new E values to accumulated sums
         for (int i=0; i<3; ++i) E_vol[i] += E_current[i];
         N3_sum += N*N*N;
         ++N;
      } while (ok == false);
      
      // Store the computed volume-average
      cell->parameters[CellParams::BGEXVOL] = E_vol[0] / N3_sum;
      cell->parameters[CellParams::BGEYVOL] = E_vol[1] / N3_sum;
      cell->parameters[CellParams::BGEZVOL] = E_vol[2] / N3_sum;

   }

   void ElectricSail::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      Real dx = cellParams[CellParams::DX];
      Real dy = cellParams[CellParams::DY];
      Real dz = cellParams[CellParams::DZ];
      Real x = cellParams[CellParams::XCRD] + 0.5*dx;
      Real y = cellParams[CellParams::YCRD] + 0.5*dy;
      Real z = cellParams[CellParams::ZCRD] + 0.5*dz;
      Real R = sqrt(x*x + y*y + z*z);
      
      cellParams[CellParams::RHOQ_TOT] = 0;
      
      if (R > radius) return;
      
      const Real volume = dx*dy*dz;
      cellParams[CellParams::RHOQ_TOT] = physicalconstants::CHARGE
              / volume / physicalconstants::EPS_0;
      
      if (Parameters::isRestart == true) return;
      cellParams[CellParams::PHI] = 0;
   }

   Real ElectricSail::calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
      
      // Iterative sampling of the distribution function. Keep track of the 
      // accumulated volume average over the iterations. When the next 
      // iteration improves the average by less than 1%, return the value.
      Real avgTotal = 0.0;
      bool ok = false;
      uint N = 2;                // Start by using a single velocity sample
      int N3_sum = 0;           // Sum of sampling points used so far
      do {
         Real avg = 0.0;        // Volume average obtained during this sampling
         creal DVX = dvx / N; 
         creal DVY = dvy / N;
         //creal DVZ = dvz / N;
         creal DVZ = dvz / 1;

         // Sample the distribution using N*N*N points
         for (uint vi=0; vi<N; ++vi) {
            for (uint vj=0; vj<N; ++vj) {
               //for (uint vk=0; vk<N; ++vk) {
               for (uint vk=0; vk<1; ++vk) {
                  creal VX = vx + 0.5*DVX + vi*DVX;
                  creal VY = vy + 0.5*DVY + vj*DVY;
                  creal VZ = vz + 0.5*DVZ + vk*DVZ;
                  avg += getDistribValue(VX,VY,VZ,DVX,DVY,DVZ,popID);
               }
            }
         }
         
         // Compare the current and accumulated volume averages:
         Real eps = max(numeric_limits<creal>::min(),avg * static_cast<Real>(1e-6));
         Real avgAccum   = avgTotal / (avg + N3_sum);
         //Real avgCurrent = avg / (N*N*N);
         Real avgCurrent = avg / (N*N);
         if (fabs(avgCurrent-avgAccum)/(avgAccum+eps) < 0.01) ok = true;
         else if (avg < getObjectWrapper().particleSpecies[popID].sparseMinValue*0.01) ok = true;
         else if (N > 10) {
            ok = true;
         }

         avgTotal += avg;
         //N3_sum += N*N*N;
         N3_sum += N*N;
         ++N;
      } while (ok == false);

      return avgTotal / N3_sum;
   }

   std::vector<std::array<Real, 3>> ElectricSail::getV0(creal x,creal y,creal z, const uint popID) const {
      vector<std::array<Real, 3>> centerPoints;
      for(uint i=0; i<populations.size(); ++i) {
         std::array<Real,3> point {{populations[i].V[0],populations[i].V[1],populations[i].V[2]}};
         centerPoints.push_back(point);
      }
      return centerPoints;
   }

   void ElectricSail::tetherElectricField(Real* x,Real* E) const {
      const Real minRadius2 = 5.0*5.0;
      Real constant = tetherUnitCharge / (2*M_PI*physicalconstants::EPS_0);
      if (timeDependentCharge == true) {
         Real factor = max((Real)0.0,(Real)1.0+(Parameters::t-tetherChargeRiseTime)/tetherChargeRiseTime);
         factor = min((Real)1.0,factor);
         constant *= factor;
      }

      E[0] = constant * (x[0] - tether_x);
      E[1] = constant * (x[1] - tether_y);
      Real radius2 = (x[0]-tether_x)*(x[0]-tether_x) + (x[1]-tether_y)*(x[1]-tether_y);
      radius2 = max(minRadius2,radius2);
      E[0] /= radius2;
      E[1] /= radius2;
   }

} // namespace projects
