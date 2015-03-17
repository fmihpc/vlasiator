/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2011, 2012, 2013, 2014 Finnish Meteorological Institute
 * 
 */

/* NOTE
 * The following piece of code has to be pasted into the main loop in vlasiator.cpp
 * to run the Dispersion project and get proper bin output files.
 * This should be kept here for future reference and reuse!!
      phiprof::stop("Propagate",computedCells,"Cells");
      
      if(P::projectName == "Dispersion") {
         vector<Real> localRho(P::xcells_ini, 0.0),
                      outputRho(P::xcells_ini, 0.0),
                      localPerBy(P::xcells_ini, 0.0),
                      outputPerBy(P::xcells_ini, 0.0);
         for(uint i=0; i<cells.size(); i++) {
            if(cells[i] <= P::xcells_ini) {
               localPerBy[cells[i] - 1] = mpiGrid[cells[i]]->parameters[CellParams::PERBY];
               localRho[cells[i] - 1] = mpiGrid[cells[i]]->parameters[CellParams::RHO];
            }
         }
         
         MPI_Reduce(&(localPerBy[0]), &(outputPerBy[0]), P::xcells_ini, MPI_DOUBLE, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD);
         MPI_Reduce(&(localRho[0]), &(outputRho[0]), P::xcells_ini, MPI_DOUBLE, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD);
         
         if(myRank == MASTER_RANK) {
            FILE* outputFile = fopen("perByt.bin", "ab");
            fwrite(&(outputPerBy[0]), sizeof(outputPerBy[0]), P::xcells_ini, outputFile);
            fclose(outputFile);
            outputFile = fopen("rhot.bin", "ab");
            fwrite(&(outputRho[0]), sizeof(outputRho[0]), P::xcells_ini, outputFile);
            fclose(outputFile);
         }
      }
      
      //Move forward in time      
      ++P::tstep;
      P::t += P::dt;
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"

#include "Dispersion.h"

Real projects::Dispersion::rndRho, projects::Dispersion::rndVel[3];

using namespace std;
using namespace spatial_cell;

namespace projects {
   Dispersion::Dispersion(): Project() { }
   Dispersion::~Dispersion() { }
   
   bool Dispersion::initialize(void) {return true;}
   
   void Dispersion::addParameters() {
      typedef Readparameters RP;
      RP::add("Dispersion.B0", "Guide magnetic field strength (T)", 1.0e-9);
      RP::add("Dispersion.VX0", "Bulk velocity (m/s)", 0.0);
      RP::add("Dispersion.VY0", "Bulk velocity (m/s)", 0.0);
      RP::add("Dispersion.VZ0", "Bulk velocity (m/s)", 0.0);
      RP::add("Dispersion.angleXY", "Orientation of the guide magnetic field with respect to the x-axis in x-y plane (rad)", 0.001);
      RP::add("Dispersion.angleXZ", "Orientation of the guide magnetic field with respect to the x-axis in x-z plane (rad)", 0.001);
      RP::add("Dispersion.rho", "Number density (m^-3)", 1.0e7);
      RP::add("Dispersion.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Dispersion.magXPertAbsAmp", "Absolute amplitude of the magnetic perturbation along x (T)", 1.0e-9);
      RP::add("Dispersion.magYPertAbsAmp", "Absolute amplitude of the magnetic perturbation along y (T)", 1.0e-9);
      RP::add("Dispersion.magZPertAbsAmp", "Absolute amplitude of the magnetic perturbation along z (T)", 1.0e-9);
      RP::add("Dispersion.densityPertRelAmp", "Relative amplitude of the density perturbation", 0.1);
      RP::add("Dispersion.velocityPertAbsAmp", "Absolute amplitude of the velocity perturbation", 1.0e6);
      RP::add("Dispersion.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("Dispersion.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
      RP::add("Dispersion.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
   }
   
   void Dispersion::getParameters() {
      typedef Readparameters RP;
      Project::getParameters();
      RP::get("Dispersion.B0", this->B0);
      RP::get("Dispersion.VX0", this->VX0);
      RP::get("Dispersion.VY0", this->VY0);
      RP::get("Dispersion.VZ0", this->VZ0);
      RP::get("Dispersion.angleXY", this->angleXY);
      RP::get("Dispersion.angleXZ", this->angleXZ);
      RP::get("Dispersion.rho", this->DENSITY);
      RP::get("Dispersion.Temperature", this->TEMPERATURE);
      RP::get("Dispersion.magXPertAbsAmp", this->magXPertAbsAmp);
      RP::get("Dispersion.magYPertAbsAmp", this->magYPertAbsAmp);
      RP::get("Dispersion.magZPertAbsAmp", this->magZPertAbsAmp);
      RP::get("Dispersion.densityPertRelAmp", this->densityPertRelAmp);
      RP::get("Dispersion.velocityPertAbsAmp", this->velocityPertAbsAmp);
      RP::get("Dispersion.nSpaceSamples", this->nSpaceSamples);
      RP::get("Dispersion.nVelocitySamples", this->nVelocitySamples);
      RP::get("Dispersion.maxwCutoff", this->maxwCutoff);
   }
   
   Real Dispersion::getDistribValue(creal& vx,creal& vy, creal& vz) {
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      return exp(- mass * ((vx-this->VX0)*(vx-this->VX0) + (vy-this->VY0)*(vy-this->VY0) + (vz-this->VZ0)*(vz-this->VZ0)) / (2.0 * kb * this->TEMPERATURE));
   }
   
   Real Dispersion::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      if(vx < Parameters::vxmin + 0.5 * dvx ||
         vy < Parameters::vymin + 0.5 * dvy ||
         vz < Parameters::vzmin + 0.5 * dvz ||
         vx > Parameters::vxmax - 1.5 * dvx ||
         vy > Parameters::vymax - 1.5 * dvy ||
         vz > Parameters::vzmax - 1.5 * dvz
      ) return 0.0;
      
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      
      creal d_vx = dvx / (this->nVelocitySamples-1);
      creal d_vy = dvy / (this->nVelocitySamples-1);
      creal d_vz = dvz / (this->nVelocitySamples-1);
      Real avg = 0.0;
      
      for (uint vi=0; vi<this->nVelocitySamples; ++vi)
         for (uint vj=0; vj<this->nVelocitySamples; ++vj)
            for (uint vk=0; vk<this->nVelocitySamples; ++vk)
            {
               avg += getDistribValue(
                  vx+vi*d_vx - this->velocityPertAbsAmp * (0.5 - this->rndVel[0]),
                  vy+vj*d_vy - this->velocityPertAbsAmp * (0.5 - this->rndVel[1]),
                  vz+vk*d_vz - this->velocityPertAbsAmp * (0.5 - this->rndVel[2])
               );
            }
            
            creal result = avg *
            this->DENSITY * (1.0 + this->densityPertRelAmp * (0.5 - this->rndRho)) *
            pow(mass / (2.0 * M_PI * kb * this->TEMPERATURE), 1.5) /
            //            (Parameters::vzmax - Parameters::vzmin) / 
            (this->nVelocitySamples*this->nVelocitySamples*this->nVelocitySamples);
            if(result < this->maxwCutoff) {
               return 0.0;
            } else {
               return result;
            }
   }
   
   void Dispersion::calcCellParameters(Real* cellParams,creal& t) {
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      creal z = cellParams[CellParams::ZCRD];
      creal dz = cellParams[CellParams::DZ];
      
      CellID cellID = (int) ((x - Parameters::xmin) / dx) +
         (int) ((y - Parameters::ymin) / dy) * Parameters::xcells_ini +
         (int) ((z - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
      
      setRandomSeed(cellID);
      
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      
      this->rndRho=getRandomNumber();
      
      this->rndVel[0]=getRandomNumber();
      this->rndVel[1]=getRandomNumber();
      this->rndVel[2]=getRandomNumber();
      
      Real rndBuffer[3];
      rndBuffer[0]=getRandomNumber();
      rndBuffer[1]=getRandomNumber();
      rndBuffer[2]=getRandomNumber();

      cellParams[CellParams::PERBX] = this->magXPertAbsAmp * (0.5 - rndBuffer[0]);
      cellParams[CellParams::PERBY] = this->magYPertAbsAmp * (0.5 - rndBuffer[1]);
      cellParams[CellParams::PERBZ] = this->magZPertAbsAmp * (0.5 - rndBuffer[2]);

   }
   
   void Dispersion::setCellBackgroundField(SpatialCell* cell) {
      ConstantField bgField;
      bgField.initialize(this->B0 * cos(this->angleXY) * cos(this->angleXZ),
                         this->B0 * sin(this->angleXY) * cos(this->angleXZ),
                         this->B0 * sin(this->angleXZ));
                         
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
   }
} // namespace projects
