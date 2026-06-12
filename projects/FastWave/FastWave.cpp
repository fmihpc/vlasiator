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

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "FastWave.h"

using namespace spatial_cell;

namespace projects {
FastWave::FastWave() : TriAxisSearch() {}
FastWave::~FastWave() {}

bool FastWave::initialize(void) {
   bool success = Project::initialize();

   creal m = physicalconstants::MASS_PROTON;
   creal e = physicalconstants::CHARGE;
   creal kB = physicalconstants::K_B;
   creal gamma = 5.0 / 3.0;
   creal mu0 = physicalconstants::MU_0;

   // Three reference values: BRef, nRef, mRef=m
   // All other reference values can be derived.
   lRef = sqrt(m / (mu0 * nRef)) / e;
   uRef = BRef / sqrt(mu0 * m * nRef);
   tRef = lRef / uRef;
   pRef = m * nRef * sqr(uRef);
   TRef = pRef / (kB * nRef);

   // Defined in dimensionless units
   rho1 = A * rho0;
   u1 = A * sqrt((sqr(B0) + 2 * p0) / rho0);
   ppar1 = A * p0;
   pperp1 = A * 0.5 * (3 * gamma - 1) * p0;

   cosalpha = cos(alpha);
   sinalpha = sin(alpha);
   // Convert to SI units
   rho0 *= nRef * m;
   p0 *= pRef;
   B0 *= BRef;

   rho1 *= nRef * m;
   u1 *= uRef;
   p1 *= pRef;
   ppar1 *= pRef;
   pperp1 *= pRef;
   // wave vector
   kwave = 2 * M_PI / (lambda * lRef);

   if (verbose) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK) {
         std::cout << "wavelength = " << lambda << "\n";
         std::cout << "amplitude = " << A << "\n";
         std::cout << "angle = " << alpha / 3.1415926 * 180 << " [degree]\n";
         std::cout << "-------------------------------------------------\n";
         std::cout << "Unit conversion factors from dimensionless to SI:\n";
         std::cout << "lRef: " << lRef << " [m]\n";
         std::cout << "tRef: " << tRef << " [s]\n";
         std::cout << "vRef: " << uRef << " [m/s]\n";
         std::cout << "nRef: " << nRef << " [#]\n";
         std::cout << "TRef: " << TRef << " [K]\n";
         std::cout << "pRef: " << pRef << " [Pa]\n";
         std::cout << "BRef: " << BRef << " [T]\n";
         std::cout << "-------------------------------------------------\n";
         creal Lx = P::xmax - P::xmin;
         std::cout << "wavelength / di = " << Lx / lRef << "\n";
         // Useful for setting the Vspace ~ 10*Vth width
         Real Vth = sqrt(gamma * p0 / rho0);
         std::cout << "Vth: " << Vth << " [m/s]\n";
         const vmesh::MeshParameters& mp = vmesh::getMeshWrapper()->velocityMeshesCreation->at(0);
         Real dvx = (mp.meshLimits[1] - mp.meshLimits[0]) / (mp.gridLength[0] * mp.blockLength[0]);
         std::cout << "Vmax:" << mp.meshLimits[1] << "\n";
         std::cout << "Vmin:" << mp.meshLimits[0] << "\n";
         std::cout << "nvcells:" << mp.gridLength[0] * mp.blockLength[0] << "\n";
         std::cout << "dvx = " << dvx << " [m/s]\n";
         std::cout << "Vth/dvx = " << Vth / dvx << "\n";
         std::cout << "Vmax/Vth = " << mp.meshLimits[1] / Vth << "\n";
         if (10 * dvx > Vth) {
            std::cout << "vmesh maybe too coarse!\n";
         }
         if (5 * Vth > mp.meshLimits[1]) {
            std::cout << "vmesh extent maybe too small!\n";
         }
         std::cout << "------------------------------------------------\n";
      }
   }

   return success;
}

void FastWave::addParameters() {
   typedef Readparameters RP;
   // Dimensionless inputs
   RP::add("FastWave.B0", "Background magnetic field", 0.04);
   RP::add("FastWave.lambda", "Dimensionless wavelength", 32.0);
   RP::add("FastWave.amplitude", "Amplitude of perturbations", 1e-2);
   RP::add("FastWave.alpha", "In-plane magnetic field tilted angle in radian w.r.t. +x", 1.5707963267948966);
   RP::add("FastWave.isBiMaxwellian", "Use BiMaxwellian distribution", 1);
   RP::add("FastWave.verbose", "Turn on/off detailed information", 0);
   // SI units inputs
   RP::add("FastWave.BRef", "Reference magnetic field strength", 1e-8);
   RP::add("FastWave.nRef", "Reference number density", 1e6);
   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;
      RP::add(pop + "_FastWave.rho0", "Background mass density", 1.0);
      RP::add(pop + "_FastWave.p0", "Background pressure", 4.5e-4);
   }
}

void FastWave::getParameters() {
   Project::getParameters();

   typedef Readparameters RP;
   RP::get("FastWave.B0", B0);
   RP::get("FastWave.lambda", lambda);
   RP::get("FastWave.amplitude", A);
   RP::get("FastWave.alpha", alpha);
   RP::get("FastWave.verbose", verbose);
   RP::get("FastWave.isBiMaxwellian", isBiMaxwellian);
   RP::get("FastWave.BRef", BRef);
   RP::get("FastWave.nRef", nRef);

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;

      if (pop != "proton") {
         int myRank;
         MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
         if (myRank == MASTER_RANK) {
            std::cerr << "Linear wave test only supports proton!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
      }
      RP::get(pop + "_FastWave.rho0", rho0);
      RP::get(pop + "_FastWave.p0", p0);
   }
}

Real FastWave::getMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy,
                             creal& dvz, const uint popID) const {
   creal m = getObjectWrapper().particleSpecies[popID].mass;
   creal kB = physicalconstants::K_B;

   Real n, ux, uy, T;

   Real xkpar = x * sinalpha - y * cosalpha;
   Real du = u1 * sin(kwave * xkpar); // assume originally along +x

   creal perturb = sin(kwave * x);
   creal rho = rho0 + perturb * rho1;

   n = rho / m;
   // Assume background velocity 0
   ux = du * cosalpha;
   uy = du * sinalpha;

   creal p = p0 + perturb * p1;
   T = p / (n * kB);

   creal coef = m / (2 * M_PI * kB * T);
   // Maxwellian f(v) = f(v;n,u,T)
   creal f = n * sqrt(coef) * coef * exp(-coef * M_PI * (sqr(vx - ux) + sqr(vy - uy) + sqr(vz)));

   return f;
}

Real FastWave::getBiMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy,
                               creal& dvz, const uint popID) const {
   creal m = getObjectWrapper().particleSpecies[popID].mass;
   creal kB = physicalconstants::K_B;

   Real n, ux, uy, Tpar, Tperp;

   creal xkpar = x * sinalpha - y * cosalpha;
   creal du = u1 * sin(kwave * xkpar); // assume originally along +x
   creal perturb = sin(kwave * xkpar);

   creal rho = rho0 + perturb * rho1;

   n = rho / m;
   // Assume background velocity 0
   ux = du * sinalpha;
   uy = -du * cosalpha;

   Tpar = (p0 + perturb * ppar1) / (n * kB);
   Tperp = (p0 + perturb * pperp1) / (n * kB);

   creal coef = m / (2 * M_PI * kB);
   // BiMaxwellian f(v) = f(v;n,u,Tpar,Tperp)
   creal vpar = vx * cosalpha + vy * sinalpha;
   creal vperpx = vx * sinalpha * sinalpha - vy * sinalpha * cosalpha;
   creal vperpy = -vx * cosalpha * sinalpha + vy * cosalpha * cosalpha;

   creal f = n * sqrt(coef) * coef / (Tperp * sqrt(Tpar)) *
             exp(-coef * M_PI / Tperp * (sqr(vperpx - ux) + sqr(vperpy - uy) + sqr(vz)) - coef * M_PI / Tpar * sqr(vpar));
   // BiMaxwellian f(v) = f(v;n,u,Tpar, Tperp) with B0 in the y-direction
   //creal f = n * sqrt(coef) * coef / (Tperp * sqrt(Tpar)) *
   //          exp(-coef * M_PI / Tperp * (sqr(vx - ux) + sqr(vz)) - coef * M_PI / Tpar * sqr(vy - uy));
   return f;
}

Real FastWave::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx,
                                     creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz, const uint popID) const {
   Real f = 0.0;
   if (!isBiMaxwellian) {
      f = getMaxwellian(x + 0.5 * dx, y + 0.5 * dy, z + 0.5 * dz, vx + 0.5 * dvx, vy + 0.5 * dvy, vz + 0.5 * dvz, dvx,
                        dvy, dvz, popID);
   } else {
      f = getBiMaxwellian(x + 0.5 * dx, y + 0.5 * dy, z + 0.5 * dz, vx + 0.5 * dvx, vy + 0.5 * dvy, vz + 0.5 * dvz, dvx,
                          dvy, dvz, popID);
   }

   return f;
}

   Realf FastWave::probePhaseSpace(spatial_cell::SpatialCell *cell,
                                        const uint popID,
                                        Real vx_in, Real vy_in, Real vz_in
      ) const {
      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];
      creal initRho = 1.0;
      creal initT = 1.0;
      const std::array<Real, 3> V0 = this->getV0(x, y, z, popID)[0];
      creal initV0X = V0[0];
      creal initV0Y = V0[1];
      creal initV0Z = V0[2];
      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;
      const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
      return value;
   }

std::vector<std::array<Real, 3>> FastWave::getV0(creal x, creal y, creal z, const uint popID) const {
   std::vector<std::array<Real, 3>> V0;
   std::array<Real, 3> v = {{0.0, 0.0, 0.0}};
   V0.push_back(v);
   return V0;
}

   Realf FastWave::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      //const speciesParameters& sP = this->speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
      creal initRho = 1.0;
      creal initT = 1.0;
      const std::array<Real, 3> V0 = this->getV0(x, y, z, popID)[0];
      creal initV0X = V0[0];
      creal initV0Y = V0[1];
      creal initV0Z = V0[2];

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
               const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
               bufferData[initIndex*WID3 + k*WID2 + j*WID + i] = value;
               //lsum[0] += value;
            };
         }, rhosum);
      return rhosum;
   }

void FastWave::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

void FastWave::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   // Background field
   ConstantField bgField;
   bgField.initialize(B0 * cosalpha, B0 * sinalpha, 0.0);
   setBackgroundField(bgField, BgBGrid);

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();

#pragma omp parallel for collapse(3)
      for (int i = 0; i < localSize[0]; ++i) {
         for (int j = 0; j < localSize[1]; ++j) {
            for (int k = 0; k < localSize[2]; ++k) {
               const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
               std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);
               Real xkpar = x[0] * sinalpha - x[1] * cosalpha;
               Real Bpar = A * B0 * sin(kwave * xkpar);
               cell->at(fsgrids::bfield::PERBX) = Bpar * cosalpha;
               cell->at(fsgrids::bfield::PERBY) = Bpar * sinalpha;
               cell->at(fsgrids::bfield::PERBZ) = 0.0;
            }
         }
      }
   }
}

} // namespace projects