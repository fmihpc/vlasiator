/* This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute
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
      RP::add("ElectricSail.tether_voltage","Electric sail tether voltage",(Real)-10000.0);
      
      //RP::addComposing("ElectricSail.rho", "Number density (m^-3)");
      //RP::addComposing("ElectricSail.Tx", "Temperature (K)");
      //RP::addComposing("ElectricSail.Ty", "Temperature");
      //RP::addComposing("ElectricSail.Tz", "Temperature");
      //RP::addComposing("ElectricSail.Vx", "Bulk velocity x component (m/s)");
      //RP::addComposing("ElectricSail.Vz", "Bulk velocity z component (m/s)");
      
      projects::ReadGaussianPopulation rgp;
      rgp.addParameters("ElectricSail");
   }

   Real ElectricSail::getCorrectNumberDensity(const int& popID) const {
      return populations[popID].rho;
   }
   
   Real ElectricSail::getDistribValue(creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz,const int& popID) const {
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;
      const Population& pop = populations[popID];

      Real value 
         = pop.rho
         * pow(mass / (2.0 * M_PI * kb ), 1.5)
         * 1.0 / sqrt(pop.T[0]*pop.T[1]*pop.T[2])
         * exp(-mass * (  pow(vx-pop.V[0],2.0) / (2.0*kb*pop.T[0]) 
                        + pow(vy-pop.V[1],2.0) / (2.0*kb*pop.T[1]) 
                        + pow(vz-pop.V[2],2.0) / (2.0*kb*pop.T[2])
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
      RP::get("ElectricSail.tether_voltage",tetherVoltage);
      
      // Get parameters for each population
      /*
      vector<Real> rho;
      vector<Real> Tx;
      vector<Real> Ty;
      vector<Real> Tz;
      vector<Real> Vx;
      vector<Real> Vy;
      vector<Real> Vz;
      RP::get("ElectricSail.rho", rho);
      RP::get("ElectricSail.Tx", Tx);
      RP::get("ElectricSail.Ty", Ty);
      RP::get("ElectricSail.Tz", Tz);
      RP::get("ElectricSail.Vx", Vx);
      RP::get("ElectricSail.Vy", Vy);
      RP::get("ElectricSail.Vz", Vz);
      
      if (Tx.size() != rho.size()) success = false;
      if (Ty.size() != rho.size()) success = false;
      if (Tz.size() != rho.size()) success = false;
      if (Vx.size() != rho.size()) success = false;
      if (Vy.size() != rho.size()) success = false;
      if (Vz.size() != rho.size()) success = false;
      */
      
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

      //for (size_t i=0; i<rho.size(); ++i) {
      //   populations.push_back(projects::Population(rho[i],Tx[i],Ty[i],Tz[i],Vx[i],Vy[i],Vz[i]));
      //}

      for (size_t i=0; i<gaussPops.rho.size(); ++i) {
         populations.push_back(projects::Population(gaussPops.rho[i],gaussPops.Tx[i],gaussPops.Ty[i],
                                                    gaussPops.Tz[i],gaussPops.Vx[i],gaussPops.Vy[i],gaussPops.Vz[i]));
      }
      
      tether_y = 2*Parameters::dy_ini;
      
#warning TESTING FIXME
      tetherUnitCharge = 100e9 * 1.602e-19;
      //tetherUnitCharge = 0.0;
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
      logFile << endl << writeVerbose;

      return success;
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

      const Real EPSILON = 1e-30;
      int N = 1;
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
         for (int k=0; k<N; ++k) for (int j=0; j<N; ++j) for (int i=0; i<N; ++i) {
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

   void ElectricSail::calcCellParameters(Real* cellParams,creal& t) {
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
            creal& dvx, creal& dvy, creal& dvz,const int& popID) const {
      //this->popID = popID;
      
// #warning TESTING remove me
//      if (y < Parameters::dy_ini) return 0.0;
//      if (y > 2*Parameters::dy_ini) return 0.0;
//      if (x < -2e5) return 0.0;
//      if (x+dx > 2e5) return 0.0;
      
      // Iterative sampling of the distribution function. Keep track of the 
      // accumulated volume average over the iterations. When the next 
      // iteration improves the average by less than 1%, return the value.
      Real avgTotal = 0.0;
      bool ok = false;
      int N = 2;                // Start by using a single velocity sample
      int N3_sum = 0;           // Sum of sampling points used so far
      do {
         Real avg = 0.0;        // Volume average obtained during this sampling
         creal DVX = dvx / N; 
         creal DVY = dvy / N;
         creal DVZ = dvz / N;

         // Sample the distribution using N*N*N points
         for (uint vi=0; vi<N; ++vi) {
            for (uint vj=0; vj<N; ++vj) {
               for (uint vk=0; vk<N; ++vk) {
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
         Real avgCurrent = avg / (N*N*N);
         if (fabs(avgCurrent-avgAccum)/(avgAccum+eps) < 0.01) ok = true;
         else if (avg < getObjectWrapper().particleSpecies[popID].sparseMinValue*0.01) ok = true;
         else if (N > 10) {
            ok = true;
         }

         avgTotal += avg;
         N3_sum += N*N*N;
         ++N;
      } while (ok == false);

      #warning TESTING remove me
      Real rad = sqrt(x*x + (y-1e5)*(y-1e5) + z*z);
      if (rad < 2e4) avgTotal *= 2;
      
      return avgTotal / N3_sum;
   }

   vector<std::array<Real, 3>> ElectricSail::getV0(creal x,creal y,creal z) {
      vector<std::array<Real, 3>> centerPoints;
      for(uint i=0; i<populations.size(); ++i) {
         std::array<Real,3> point {{populations[i].V[0],populations[i].V[1],populations[i].V[2]}};
         centerPoints.push_back(point);
      }
      return centerPoints;
   }

   void ElectricSail::tetherElectricField(Real* x,Real* E) const {
      const Real constant = tetherUnitCharge / (2*M_PI*physicalconstants::EPS_0);
      E[0] = constant * (x[0] - tether_x);
      E[1] = constant * (x[1] - tether_y);

      Real radius2 = (x[0]-tether_x)*(x[0]-tether_x) + (x[1]-tether_y)*(x[1]-tether_y);

      radius2 = max(1e-12,radius2);
      E[0] /= radius2;
      E[1] /= radius2;
   }

} // namespace projects
