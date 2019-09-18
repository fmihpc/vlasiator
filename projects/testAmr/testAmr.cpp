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
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../object_wrapper.h"

#include "testAmr.h"

using namespace std;
using namespace spatial_cell;

Real projects::testAmr::rhoRnd;

namespace projects {
   testAmr::testAmr(): TriAxisSearch() { }
   
   testAmr::~testAmr() { }

   bool testAmr::initialize(void) {
      return Project::initialize();
   }

   void testAmr::addParameters(){
      typedef Readparameters RP;

      RP::add("testAmr.Bx", "Magnetic field x component (T)", 0.0);
      RP::add("testAmr.By", "Magnetic field y component (T)", 0.0);
      RP::add("testAmr.Bz", "Magnetic field z component (T)", 0.0);
      RP::add("testAmr.dBx", "Magnetic field x component cosine perturbation amplitude (T)", 0.0);
      RP::add("testAmr.dBy", "Magnetic field y component cosine perturbation amplitude (T)", 0.0);
      RP::add("testAmr.dBz", "Magnetic field z component cosine perturbation amplitude (T)", 0.0);
      RP::add("testAmr.magXPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along x (T)", 1.0e-9);
      RP::add("testAmr.magYPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along y (T)", 1.0e-9);
      RP::add("testAmr.magZPertAbsAmp", "Absolute amplitude of the random magnetic perturbation along z (T)", 1.0e-9);
      RP::add("testAmr.lambda", "B cosine perturbation wavelength (m)", 1.0);
      RP::add("testAmr.nVelocitySamples", "Number of sampling points per velocity dimension", 2);
      RP::add("testAmr.densityModel","Which spatial density model is used?",string("uniform"));
      RP::add("testAmr.maxSpatialRefinementLevel", "Maximum level for spatial refinement", 1.0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         RP::add(pop+"_testAmr.n", "Number of peaks to create", 0);
         RP::addComposing(pop+"_testAmr.rho", "Number density (m^-3)");
         RP::addComposing(pop+"_testAmr.Tx", "Temperature (K)");
         RP::addComposing(pop+"_testAmr.Ty", "Temperature");
         RP::addComposing(pop+"_testAmr.Tz", "Temperature");
         RP::addComposing(pop+"_testAmr.Vx", "Bulk velocity x component (m/s)");
         RP::addComposing(pop+"_testAmr.Vy", "Bulk velocity y component (m/s)");
         RP::addComposing(pop+"_testAmr.Vz", "Bulk velocity z component (m/s)");
         RP::addComposing(pop+"_testAmr.rhoPertAbsAmp", "Absolute amplitude of the density perturbation");
      }
   }

   void testAmr::getParameters(){

      typedef Readparameters RP;
      Project::getParameters();
      RP::get("testAmr.Bx", this->Bx);
      RP::get("testAmr.By", this->By);
      RP::get("testAmr.Bz", this->Bz);
      RP::get("testAmr.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("testAmr.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("testAmr.magZPertAbsAmp", this->magZPertAbsAmp);
      RP::get("testAmr.dBx", this->dBx);
      RP::get("testAmr.dBy", this->dBy);
      RP::get("testAmr.dBz", this->dBz);
      RP::get("testAmr.lambda", this->lambda);
      RP::get("testAmr.maxSpatialRefinementLevel", this->maxSpatialRefinementLevel);
      RP::get("testAmr.nVelocitySamples", this->nVelocitySamples);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         testAmrSpeciesParameters sP;
         RP::get(pop + "_testAmr.n",  sP.numberOfPeaks);
         RP::get(pop + "_testAmr.rho",sP.rho);
         RP::get(pop + "_testAmr.Tx", sP.Tx);
         RP::get(pop + "_testAmr.Ty", sP.Ty);
         RP::get(pop + "_testAmr.Tz", sP.Tz);
         RP::get(pop + "_testAmr.Vx", sP.Vx);
         RP::get(pop + "_testAmr.Vy", sP.Vy);
         RP::get(pop + "_testAmr.Vz", sP.Vz);

         RP::get(pop + "_testAmr.rhoPertAbsAmp", sP.rhoPertAbsAmp);
         
         if(!sP.isConsistent()) {
            cerr << "You should define all parameters (testAmr.rho, testAmr.Tx, testAmr.Ty, testAmr.Tz, testAmr.Vx, testAmr.Vy, testAmr.Vz, testAmr.rhoPertAbsAmp) for all " << sP.numberOfPeaks << " peaks of population " << pop << "." << endl;
            abort();
         }

         speciesParams.push_back(sP);
      }
      
      string densModelString;
      RP::get("testAmr.densityModel",densModelString);
      
      if (densModelString == "uniform") densityModel = Uniform;
      else if (densModelString == "testcase") densityModel = TestCase;
   }

   Real testAmr::getDistribValue(creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
      const testAmrSpeciesParameters& sP = speciesParams[popID];
      creal mass = getObjectWrapper().particleSpecies[popID].mass;
      creal kb = physicalconstants::K_B;

      Real value = 0.0;

      for (uint i=0; i<sP.numberOfPeaks; ++i) {
         value += (sP.rho[i] + sP.rhoPertAbsAmp[i] * rhoRnd)
               * pow(mass / (2.0 * M_PI * kb ), 1.5) * 1.0
               / sqrt(sP.Tx[i]*sP.Ty[i]*sP.Tz[i]) 
               * exp(- mass * (pow(vx - sP.Vx[i], 2.0) / (2.0 * kb * sP.Tx[i]) 
                             + pow(vy - sP.Vy[i], 2.0) / (2.0 * kb * sP.Ty[i]) 
                             + pow(vz - sP.Vz[i], 2.0) / (2.0 * kb * sP.Tz[i])));
      }
      return value;
   }

   Real testAmr::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, 
                                       creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,
                                       const uint popID) const {
      // Iterative sampling of the distribution function. Keep track of the 
      // accumulated volume average over the iterations. When the next 
      // iteration improves the average by less than 1%, return the value.
      Real avgTotal = 0.0;
      bool ok = false;
      uint N = nVelocitySamples; // Start by using nVelocitySamples
      int N3_sum = 0;           // Sum of sampling points used so far

      const testAmrSpeciesParameters& sP = speciesParams[popID];

      const Real avgLimit = 0.01*getObjectWrapper().particleSpecies[popID].sparseMinValue;
      do {
         Real avg = 0.0;        // Volume average obtained during this sampling
         creal DVX = dvx / N; 
         creal DVY = dvy / N;
         creal DVZ = dvz / N;

         Real rhoFactor = 1.0;
         switch (densityModel) {
            case Uniform:
               rhoFactor = 1.0;
               break;
            case TestCase:
               rhoFactor = 1.0;
               if (x < P::xmin + 0.31 * (P::xmax - P::xmin) &&
                   y < P::ymin + 0.31 * (P::ymax - P::ymin)) {
                  rhoFactor = 3.0;
               }
               break;
            default:
               rhoFactor = 1.0;
               break;
         }

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
         avg *= rhoFactor;
         
         // Compare the current and accumulated volume averages:
         Real eps = max(numeric_limits<creal>::min(),avg * static_cast<Real>(1e-6));
         Real avgAccum   = avgTotal / (avg + N3_sum);
         Real avgCurrent = avg / (N*N*N);
         if (fabs(avgCurrent-avgAccum)/(avgAccum+eps) < 0.01) ok = true;
         else if (avg < avgLimit) ok = true;
         else if (N > 10) {
            ok = true;
         }

         avgTotal += avg;
         N3_sum += N*N*N;
         ++N;
      } while (ok == false);

      return avgTotal / N3_sum;
   }

   void testAmr::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      setRandomCellSeed(cell);
      rhoRnd = 0.5 - getRandomNumber();
   }

   void testAmr::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
      FsGrid< fsgrids::technical, 2>& technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(this->Bx,
                         this->By,
                         this->Bz);
      
      setBackgroundField(bgField, BgBGrid);
      
      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();
         
         #pragma omp parallel for collapse(3)
         for (int x = 0; x < localSize[0]; ++x) {
            for (int y = 0; y < localSize[1]; ++y) {
               for (int z = 0; z < localSize[2]; ++z) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  
                  const int64_t cellid = perBGrid.GlobalIDForCoords(x, y, z);
                  
                  setRandomSeed(cellid);
                  
                  if (this->lambda != 0.0) {
                     cell->at(fsgrids::bfield::PERBX) = this->dBx*cos(2.0 * M_PI * xyz[0] / this->lambda);
                     cell->at(fsgrids::bfield::PERBY) = this->dBy*sin(2.0 * M_PI * xyz[0] / this->lambda);
                     cell->at(fsgrids::bfield::PERBZ) = this->dBz*cos(2.0 * M_PI * xyz[0] / this->lambda);
                  }
                  
                  cell->at(fsgrids::bfield::PERBX) += this->magXPertAbsAmp * (0.5 - getRandomNumber());
                  cell->at(fsgrids::bfield::PERBY) += this->magYPertAbsAmp * (0.5 - getRandomNumber());
                  cell->at(fsgrids::bfield::PERBZ) += this->magZPertAbsAmp * (0.5 - getRandomNumber());
               }
            }
         }
      }
   }
   
   std::vector<std::array<Real, 3> > testAmr::getV0(
                                                creal x,
                                                creal y,
                                                creal z,
                                                const uint popID
                                               ) const {
      const testAmrSpeciesParameters& sP = speciesParams[popID];
      vector<std::array<Real, 3> > centerPoints;
      for(uint i=0; i<sP.numberOfPeaks; i++) {
         array<Real, 3> point {{sP.Vx[i], sP.Vy[i], sP.Vz[i]}};
         centerPoints.push_back(point);
      }
      return centerPoints;
   }

   bool testAmr::refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const {

     int myRank;       
     MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

     if(myRank == MASTER_RANK) std::cout << "Maximum refinement level is " << mpiGrid.mapping.get_maximum_refinement_level() << std::endl;
      
      std::vector<bool> refineSuccess;

      for (int i = 0; i < 2 * P::amrBoxHalfWidthX; ++i) {
         for (int j = 0; j < 2 * P::amrBoxHalfWidthY; ++j) {
            for (int k = 0; k < 2 * P::amrBoxHalfWidthZ; ++k) {
     
               std::array<double,3> xyz;
               xyz[0] = P::amrBoxCenterX + (0.5 + i - P::amrBoxHalfWidthX) * P::dx_ini;
               xyz[1] = P::amrBoxCenterY + (0.5 + j - P::amrBoxHalfWidthY) * P::dy_ini;
               xyz[2] = P::amrBoxCenterZ + (0.5 + k - P::amrBoxHalfWidthZ) * P::dz_ini;
               
               CellID myCell = mpiGrid.get_existing_cell(xyz);
               if (mpiGrid.refine_completely_at(xyz)) {
#ifndef NDEBUG
                  std::cout << "Rank " << myRank << " is refining cell " << myCell << std::endl;
#endif
               }
            }
         }
      }
      std::vector<CellID> refinedCells = mpiGrid.stop_refining(true);      
      if(myRank == MASTER_RANK) std::cout << "Finished first level of refinement" << endl;
#ifndef NDEBUG
      if(refinedCells.size() > 0) {
	std::cout << "Refined cells produced by rank " << myRank << " are: ";
	for (auto cellid : refinedCells) {
	  std::cout << cellid << " ";
	}
	std::cout << endl;
      }
#endif
                  
      mpiGrid.balance_load();

      if(mpiGrid.get_maximum_refinement_level() > 1) {

         for (int i = 0; i < 2 * P::amrBoxHalfWidthX; ++i) {
            for (int j = 0; j < 2 * P::amrBoxHalfWidthY; ++j) {
               for (int k = 0; k < 2 * P::amrBoxHalfWidthZ; ++k) {
                  
                  std::array<double,3> xyz;
                  xyz[0] = P::amrBoxCenterX + 0.5 * (0.5 + i - P::amrBoxHalfWidthX) * P::dx_ini;
                  xyz[1] = P::amrBoxCenterY + 0.5 * (0.5 + j - P::amrBoxHalfWidthY) * P::dy_ini;
                  xyz[2] = P::amrBoxCenterZ + 0.5 * (0.5 + k - P::amrBoxHalfWidthZ) * P::dz_ini;
                  
                  CellID myCell = mpiGrid.get_existing_cell(xyz);
                  if (mpiGrid.refine_completely_at(xyz)) {
#ifndef NDEBUG
                     std::cout << "Rank " << myRank << " is refining cell " << myCell << std::endl;
#endif
                  }
               }
            }
         }
         
         std::vector<CellID> refinedCells = mpiGrid.stop_refining(true);      
         if(myRank == MASTER_RANK) std::cout << "Finished second level of refinement" << endl;
#ifndef NDEBUG
         if(refinedCells.size() > 0) {
            std::cout << "Refined cells produced by rank " << myRank << " are: ";
            for (auto cellid : refinedCells) {
               std::cout << cellid << " ";
            }
            std::cout << endl;
         }
#endif              
         mpiGrid.balance_load();
      }
      
      return true;
   }
   
}// namespace projects
