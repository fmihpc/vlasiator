
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../object_wrapper.h"

#include "KH_instability.h"

using namespace std;

namespace projects {
   KH_instability::KH_instability(): Project() { }
   KH_instability::~KH_instability() { }

   bool KH_instability::initialize(void) {

      // generate random numbers for each cell, used for bulk velocity noise
      const uint length = Parameters::xcells_ini * Parameters::ycells_ini * Parameters::zcells_ini + 1;
      this->random_x.resize(length);
      this->random_y.resize(length);
      this->random_z.resize(length);
      setRandomSeed(1337);
      for (uint i=0; i<length; i++) {
         this->random_x[i] = getRandomNumber();
         this->random_y[i] = getRandomNumber();
         this->random_z[i] = getRandomNumber();
      }

      return Project::initialize();
   }

   void KH_instability::addParameters(){
      typedef Readparameters RP;
      RP::add("KH_instability.Vx1", "Initial bulk velocity in x-direction, LEFT (m/s)", 0.0);
      RP::add("KH_instability.Vx2", "Initial bulk velocity in x-direction, RIGHT (m/s)", 0.0);
      RP::add("KH_instability.Vy1", "Initial bulk velocity in y-direction, LEFT (m/s)", 0.0);
      RP::add("KH_instability.Vy2", "Initial bulk velocity in y-direction, RIGHT (m/s)", 0.0);
      RP::add("KH_instability.Vz1", "Initial bulk velocity in z-direction, LEFT (m/s)", 0.0);
      RP::add("KH_instability.Vz2", "Initial bulk velocity in z-direction, RIGHT (m/s)", 0.0);
      RP::add("KH_instability.Bx1", "Magnetic field x component, LEFT (T)", 0.0);
      RP::add("KH_instability.Bx2", "Magnetic field x component, RIGHT (T)", 0.0);
      RP::add("KH_instability.By1", "Magnetic field y component, LEFT (T)", 0.0);
      RP::add("KH_instability.By2", "Magnetic field y component, RIGHT (T)", 0.0);
      RP::add("KH_instability.Bz1", "Magnetic field z component, LEFT (T)", 0.0);
      RP::add("KH_instability.Bz2", "Magnetic field z component, RIGHT (T)", 0.0);    
      RP::add("KH_instability.transitionWidth", "Width of the tanh-shaped velocity/magnetic shear layer (m)", 1.0);
      RP::add("KH_instability.perturbationAmplitude", "Amplitude of the initial velocity perturbation (fraction of velocity shear)", 0.0);
      RP::add("KH_instability.perturbationWaveNumber", "Wave number of the initial velocity perturbation (dimensionless, actual wavenumber multiplied by 2*transitionWidth)", 0.0);
      RP::add("KH_instability.noiseAmplitude", "Amplitude of noise added to the initial velocity perturbation (fraction of velocity shear)", 0.0); 

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         RP::add(pop + "_KH_instability.rho1", "Number density, LEFT (m^-3)", 0.0);
         RP::add(pop + "_KH_instability.rho2", "Number density, RIGHT (m^-3)", 0.0);
	     RP::add(pop + "_KH_instability.P", "Thermal pressure (Pa)", 0.0);
         RP::add(pop + "_KH_instability.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
         RP::add(pop + "_KH_instability.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
      }  
   }
   
   void KH_instability::getParameters(){
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

      Project::getParameters();

      typedef Readparameters RP;
      RP::get("KH_instability.Vx1", this->Vx[LEFT]);
      RP::get("KH_instability.Vx2", this->Vx[RIGHT]);
      RP::get("KH_instability.Vy1", this->Vy[LEFT]);
      RP::get("KH_instability.Vy2", this->Vy[RIGHT]);
      RP::get("KH_instability.Vz1", this->Vz[LEFT]);
      RP::get("KH_instability.Vz2", this->Vz[RIGHT]);
      RP::get("KH_instability.Bx1", this->Bx[LEFT]);
      RP::get("KH_instability.Bx2", this->Bx[RIGHT]);
      RP::get("KH_instability.By1", this->By[LEFT]);
      RP::get("KH_instability.By2", this->By[RIGHT]);
      RP::get("KH_instability.Bz1", this->Bz[LEFT]);
      RP::get("KH_instability.Bz2", this->Bz[RIGHT]);
      RP::get("KH_instability.transitionWidth", this->transitionWidth);
      RP::get("KH_instability.perturbationAmplitude", this->perturbationAmplitude);
      RP::get("KH_instability.perturbationWaveNumber", this->perturbationWaveNumber);
      RP::get("KH_instability.noiseAmplitude", this->noiseAmplitude);

      // Per-population parameters
      for (uint i=0; i<getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        KH_instabilitySpeciesParameters sP;
        RP::get(pop + "_KH_instability.rho1", sP.rho[LEFT]);
        RP::get(pop + "_KH_instability.rho2", sP.rho[RIGHT]);
        RP::get(pop + "_KH_instability.P", sP.P);
        RP::get(pop + "_KH_instability.nSpaceSamples", sP.nSpaceSamples);
        RP::get(pop + "_KH_instability.nVelocitySamples", sP.nVelocitySamples);

        speciesParams.push_back(sP);
      }
    
   }

   Real KH_instability::profile(creal left, creal right, creal x, creal y, creal z) const {
      return 0.5 * ((right - left) * tanh(x / this->transitionWidth) + right + left);
   }

   Real KH_instability::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& rand_x, creal& rand_y, creal& rand_z, const uint popID) const {

      Real mass = getObjectWrapper().particleSpecies[popID].mass;
      const KH_instabilitySpeciesParameters& sP = speciesParams[popID];
      Real delta_v = this->Vy[RIGHT] - this->Vy[LEFT];
      
      // single maxwellian, pressure balance implementation
      Real rho = profile(sP.rho[LEFT], sP.rho[RIGHT], x, y, z);

      Real V0x = perturbationAmplitude * delta_v * sin(0.5 * perturbationWaveNumber * y / this->transitionWidth) * exp(-pow(x / this->transitionWidth,2))
	       + 0.001 * perturbationAmplitude * delta_v * sin(0.25 * perturbationWaveNumber * y / this->transitionWidth) * exp(-pow(x / this->transitionWidth,2))
               + this->noiseAmplitude * delta_v * 2.0 * (0.5 - rand_x);

      Real V0y = profile(this->Vy[LEFT], this->Vy[RIGHT], x, y, z)
               + this->noiseAmplitude * delta_v * 2.0 * (0.5 - rand_y);

      Real V0z = this->noiseAmplitude * delta_v * 2.0 * (0.5 - rand_z);

      Real kBT = sP.P / rho;

      Real rvalue = 0;
      // Maxwellian
      rvalue = rho * pow(mass / (2.0 * M_PI * kBT), 1.5)
                   * exp(-mass * (  (vx-V0x)*(vx-V0x)
                                  + (vy-V0y)*(vy-V0y)
                                  + (vz-V0z)*(vz-V0z)
                                 ) / (2.0 * kBT));

      return rvalue;    
   }

   Real KH_instability::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {

      const KH_instabilitySpeciesParameters& sP = speciesParams[popID];

      // fetch pre-calculated random values
      int cellid = (int) ((x + 0.5 * dx - Parameters::xmin) / dx) +
                   (int) ((y + 0.5 * dy - Parameters::ymin) / dy) * Parameters::xcells_ini +
                   (int) ((z + 0.5 * dz - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
      creal rand_x = this->random_x[cellid];
      creal rand_y = this->random_y[cellid];
      creal rand_z = this->random_z[cellid];

      if((sP.nSpaceSamples > 1) && (sP.nVelocitySamples > 1)) {
         creal d_x = dx / (sP.nSpaceSamples-1);
         creal d_y = dy / (sP.nSpaceSamples-1);
         creal d_z = dz / (sP.nSpaceSamples-1);
         creal d_vx = dvx / (sP.nVelocitySamples-1);
         creal d_vy = dvy / (sP.nVelocitySamples-1);
         creal d_vz = dvz / (sP.nVelocitySamples-1);

         Real avg = 0.0;
         for (uint i=0; i<sP.nSpaceSamples; ++i) for (uint j=0; j<sP.nSpaceSamples; ++j) for (uint k=0; k<sP.nSpaceSamples; ++k) {
            for (uint vi=0; vi<sP.nVelocitySamples; ++vi) for (uint vj=0; vj<sP.nVelocitySamples; ++vj) for (uint vk=0; vk<sP.nVelocitySamples; ++vk) {
               avg += getDistribValue(x+i*d_x,y+j*d_y,z+k*d_z,vx+vi*d_vx,vy+vj*d_vy,vz+vk*d_vz,rand_x,rand_y,rand_z,popID);
            }
         }
         return avg / (sP.nSpaceSamples*sP.nSpaceSamples*sP.nSpaceSamples*sP.nVelocitySamples*sP.nVelocitySamples*sP.nVelocitySamples);
      } else {
         return getDistribValue(x+0.5*dx,y+0.5*dy,z+0.5*dz,vx+0.5*dvx,vy+0.5*dvy,vz+0.5*dvz,rand_x,rand_y,rand_z,popID);
      }
   }

   void KH_instability::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }

   void KH_instability::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      setBackgroundFieldToZero(BgBGrid);

      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();

#pragma omp parallel for collapse(3)
         for (int x = 0; x < localSize[0]; ++x) {
            for (int y = 0; y < localSize[1]; ++y) {
               for (int z = 0; z < localSize[2]; ++z) {
                   const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                   std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);

                   // XXX does this have pressure balance?
                   cell->at(fsgrids::bfield::PERBX) = profile(this->Bx[LEFT], this->Bx[RIGHT], xyz[0]+0.5*perBGrid.DX, xyz[1]+0.5*perBGrid.DY, xyz[2]+0.5*perBGrid.DZ);
                   cell->at(fsgrids::bfield::PERBY) = profile(this->By[LEFT], this->By[RIGHT], xyz[0]+0.5*perBGrid.DX, xyz[1]+0.5*perBGrid.DY, xyz[2]+0.5*perBGrid.DZ);
                   cell->at(fsgrids::bfield::PERBZ) = profile(this->Bz[LEFT], this->Bz[RIGHT], xyz[0]+0.5*perBGrid.DX, xyz[1]+0.5*perBGrid.DY, xyz[2]+0.5*perBGrid.DZ);
               }
            }
         }
      }
     }

} //namespace projects
