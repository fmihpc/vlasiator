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
#include "../../object_wrapper.h"
#include "../../velocity_mesh_parameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"

#include "Dispersion.h"

Real projects::Dispersion::rndRho, projects::Dispersion::rndVel[3];

using namespace std;
using namespace spatial_cell;

namespace projects {
   Dispersion::Dispersion(): Project() { }
   Dispersion::~Dispersion() { }
   
   bool Dispersion::initialize(void) {return Project::initialize();}
   
   void Dispersion::addParameters() {
      typedef Readparameters RP;
      RP::add("Dispersion.B0", "Guide magnetic field strength (T)", 1.0e-9);
      RP::add("Dispersion.magXPertAbsAmp", "Absolute amplitude of the magnetic perturbation along x (T)", 1.0e-9);
      RP::add("Dispersion.magYPertAbsAmp", "Absolute amplitude of the magnetic perturbation along y (T)", 1.0e-9);
      RP::add("Dispersion.magZPertAbsAmp", "Absolute amplitude of the magnetic perturbation along z (T)", 1.0e-9);
      RP::add("Dispersion.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
      RP::add("Dispersion.angleXY", "Orientation of the guide magnetic field with respect to the x-axis in x-y plane (rad)", 0.001);
      RP::add("Dispersion.angleXZ", "Orientation of the guide magnetic field with respect to the x-axis in x-z plane (rad)", 0.001);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        RP::add(pop + "_Dispersion.VX0", "Bulk velocity (m/s)", 0.0);
        RP::add(pop + "_Dispersion.VY0", "Bulk velocity (m/s)", 0.0);
        RP::add(pop + "_Dispersion.VZ0", "Bulk velocity (m/s)", 0.0);
        RP::add(pop + "_Dispersion.rho", "Number density (m^-3)", 1.0e7);
        RP::add(pop + "_Dispersion.Temperature", "Temperature (K)", 2.0e6);
        RP::add(pop + "_Dispersion.densityPertRelAmp", "Relative amplitude of the density perturbation", 0.1);
        RP::add(pop + "_Dispersion.velocityPertAbsAmp", "Absolute amplitude of the velocity perturbation", 1.0e6);
      }
   }
   
   void Dispersion::getParameters() {
      Project::getParameters();
      typedef Readparameters RP;
      Project::getParameters();
      RP::get("Dispersion.B0", this->B0);
      RP::get("Dispersion.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("Dispersion.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("Dispersion.magZPertAbsAmp", this->magZPertAbsAmp);
      RP::get("Dispersion.maxwCutoff", this->maxwCutoff);
      RP::get("Dispersion.angleXY", this->angleXY);
      RP::get("Dispersion.angleXZ", this->angleXZ);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        DispersionSpeciesParameters sP;
        RP::get(pop + "_Dispersion.VX0", sP.VX0);
        RP::get(pop + "_Dispersion.VY0", sP.VY0);
        RP::get(pop + "_Dispersion.VZ0", sP.VZ0);
        RP::get(pop + "_Dispersion.rho", sP.DENSITY);
        RP::get(pop + "_Dispersion.Temperature", sP.TEMPERATURE);
        RP::get(pop + "_Dispersion.densityPertRelAmp", sP.densityPertRelAmp);
        RP::get(pop + "_Dispersion.velocityPertAbsAmp", sP.velocityPertAbsAmp);

         speciesParams.push_back(sP);
      }
   }
   
   void Dispersion::hook(
      cuint& stage,
      const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid
   ) const {
      if(hook::END_OF_TIME_STEP == stage) {
         int myRank;
         MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

         vector<Real> localRhom(P::xcells_ini, 0.0), outputRhom(P::xcells_ini, 0.0);

         const vector<CellID>& cells = getLocalCells();

         for(uint i=0; i<cells.size(); i++) {
            if(cells[i] <= P::xcells_ini) {
               localRhom[cells[i] - 1] = mpiGrid[cells[i]]->parameters[CellParams::RHOM];
            }
         }
         
         MPI_Reduce(&(localRhom[0]), &(outputRhom[0]), P::xcells_ini, MPI_DOUBLE, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD);

         vector<Real> localPerBx(P::xcells_ini, 0.0);
         vector<Real> localPerBy(P::xcells_ini, 0.0);
         vector<Real> localPerBz(P::xcells_ini, 0.0);
         vector<Real> outputPerBx(P::xcells_ini, 0.0);
         vector<Real> outputPerBy(P::xcells_ini, 0.0);
         vector<Real> outputPerBz(P::xcells_ini, 0.0);
         
         const std::array<FsGridTools::FsIndex_t, 3> localSize = perBGrid.getLocalSize();
         const std::array<FsGridTools::FsIndex_t, 3> localStart = perBGrid.getLocalStart();
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            localPerBx[x + localStart[0]] = perBGrid.get(x, 0, 0)->at(fsgrids::bfield::PERBX);
            localPerBy[x + localStart[0]] = perBGrid.get(x, 0, 0)->at(fsgrids::bfield::PERBY);
            localPerBz[x + localStart[0]] = perBGrid.get(x, 0, 0)->at(fsgrids::bfield::PERBZ);
         }
         
         MPI_Reduce(&(localPerBx[0]), &(outputPerBx[0]), P::xcells_ini, MPI_DOUBLE, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD);
         MPI_Reduce(&(localPerBy[0]), &(outputPerBy[0]), P::xcells_ini, MPI_DOUBLE, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD);
         MPI_Reduce(&(localPerBz[0]), &(outputPerBz[0]), P::xcells_ini, MPI_DOUBLE, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD);
         
         if(myRank == MASTER_RANK) {
            FILE* outputFile = fopen("perBxt.bin", "ab");
            fwrite(&(outputPerBx[0]), sizeof(outputPerBx[0]), P::xcells_ini, outputFile);
            fclose(outputFile);
            outputFile = fopen("perByt.bin", "ab");
            fwrite(&(outputPerBy[0]), sizeof(outputPerBy[0]), P::xcells_ini, outputFile);
            fclose(outputFile);
            outputFile = fopen("perBzt.bin", "ab");
            fwrite(&(outputPerBz[0]), sizeof(outputPerBz[0]), P::xcells_ini, outputFile);
            fclose(outputFile);
            outputFile = fopen("rhomt.bin", "ab");
            fwrite(&(outputRhom[0]), sizeof(outputRhom[0]), P::xcells_ini, outputFile);
            fclose(outputFile);
         }
      }
   }
   
   Realf Dispersion::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      const DispersionSpeciesParameters& sP = speciesParams[popID];
      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      Real initT = sP.TEMPERATURE;
      Real initRho = sP.DENSITY * (1.0 + sP.densityPertRelAmp * (0.5 - this->rndRho));
      const Real initV0X = sP.VX0 + sP.velocityPertAbsAmp * (0.5 - this->rndVel[0]);
      const Real initV0Y = sP.VY0 + sP.velocityPertAbsAmp * (0.5 - this->rndVel[1]);
      const Real initV0Z = sP.VZ0 + sP.velocityPertAbsAmp * (0.5 - this->rndVel[2]);

      #ifdef USE_GPU
      vmesh::VelocityMesh *vmesh = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh *vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      #endif
      // Loop over blocks
      Realf rhosum = 0;
      arch::parallel_reduce<arch::null>(
         {WID, WID, WID, nRequested},
         ARCH_LOOP_LAMBDA (const uint i, const uint j, const uint k, const uint initIndex, Realf *lsum ) {
            vmesh::GlobalID *GIDlist = vmesh->getGrid()->data();
            Realf* bufferData = VBC->getData();
            const vmesh::GlobalID blockGID = GIDlist[initIndex];
            // Calculate parameters for new block
            Real blockCoords[6];
            vmesh->getBlockInfo(blockGID,&blockCoords[0]);
            creal vxBlock = blockCoords[0];
            creal vyBlock = blockCoords[1];
            creal vzBlock = blockCoords[2];
            creal dvxCell = blockCoords[3];
            creal dvyCell = blockCoords[4];
            creal dvzCell = blockCoords[5];
            ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
               creal vx = vxBlock + (i+0.5)*dvxCell - initV0X;
               creal vy = vyBlock + (j+0.5)*dvyCell - initV0Y;
               creal vz = vzBlock + (k+0.5)*dvzCell - initV0Z;
               Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
               bufferData[initIndex*WID3 + k*WID2 + j*WID + i] = value;
               //lsum[0] += value;
            };
         }, rhosum);
      return rhosum;
   }

   void Dispersion::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      std::default_random_engine rndState;
      setRandomCellSeed(cell,rndState);
      
      this->rndRho=getRandomNumber(rndState);
      
      this->rndVel[0]=getRandomNumber(rndState);
      this->rndVel[1]=getRandomNumber(rndState);
      this->rndVel[2]=getRandomNumber(rndState);
   }
   
   void Dispersion::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
   ) {
      ConstantField bgField;
      bgField.initialize(this->B0 * cos(this->angleXY) * cos(this->angleXZ),
                         this->B0 * sin(this->angleXY) * cos(this->angleXZ),
                         this->B0 * sin(this->angleXZ));
                         
      setBackgroundField(bgField, BgBGrid);
      
      if(!P::isRestart) {
         const auto localSize = BgBGrid.getLocalSize().data();
         
#pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  const int64_t cellid = perBGrid.GlobalIDForCoords(x, y, z);
                  
                  std::default_random_engine rndState;
                  setRandomSeed(cellid,rndState);
                  
                  Real rndBuffer[3];
                  rndBuffer[0]=getRandomNumber(rndState);
                  rndBuffer[1]=getRandomNumber(rndState);
                  rndBuffer[2]=getRandomNumber(rndState);
                  
                  cell->at(fsgrids::bfield::PERBX) = this->magXPertAbsAmp * (0.5 - rndBuffer[0]);
                  cell->at(fsgrids::bfield::PERBY) = this->magYPertAbsAmp * (0.5 - rndBuffer[1]);
                  cell->at(fsgrids::bfield::PERBZ) = this->magZPertAbsAmp * (0.5 - rndBuffer[2]);
               }
            }
         }
      }
   }
} // namespace projects
