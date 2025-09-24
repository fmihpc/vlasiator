/*
 * This file is part of Vlasiator.
 * Copyright 2010-2025 Finnish Meteorological Institute and University of Helsinki
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

#include "../object_wrapper.h"
#include "../parameters.h"
#include <math.h>
// #include <cmath> // NaN Inf checks
#include "common_pitch_angle_diffusion.hpp"
#include "vec.h"
#include <Eigen/Geometry>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>

#define MUSPACE(var, v_ind, mu_ind) var.at((mu_ind) * nbins_v + (v_ind))

using namespace spatial_cell;
using namespace Eigen;

template <typename Lambda> inline static void loop_over_block(Lambda loop_body) {

   for (int k = 0; k < WID; ++k) {
      for (int j = 0; j < WID; j += VECL / WID) { // Iterate through coordinates (z,y)

         // create vectors with the i and j indices in the vector position on the plane.
#if VECL == 4 && WID == 4
         const Veci i_indices = Veci({0, 1, 2, 3});
         const Veci j_indices = Veci({j, j, j, j});
#elif VECL == 4 && WID == 8
#error "__FILE__ : __LINE__ : VECL == 4 && WID == 8 cannot work!"
#elif VECL == 8 && WID == 4
         const Veci i_indices = Veci({0, 1, 2, 3, 0, 1, 2, 3});
         const Veci j_indices = Veci({j, j, j, j, j + 1, j + 1, j + 1, j + 1});
#elif VECL == 8 && WID == 8
         const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7});
         const Veci j_indices = Veci({j, j, j, j, j, j, j, j});
#elif VECL == 16 && WID == 4
         const Veci i_indices = Veci({0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3});
         const Veci j_indices = Veci({j, j, j, j, j + 1, j + 1, j + 1, j + 1, j + 2, j + 2, j + 2, j + 2, j + 3, j + 3, j + 3, j + 3});
#elif VECL == 16 && WID == 8
         const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7});
         const Veci j_indices = Veci({j, j, j, j, j, j, j, j, j + 1, j + 1, j + 1, j + 1, j + 1, j + 1, j + 1, j + 1});
#elif VECL == 16 && WID == 16
         const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15});
         const Veci j_indices = Veci({j, j, j, j, j, j, j, j, j, j, j, j, j, j, j, j});
#elif VECL == 32 && WID == 4
#error "__FILE__ : __LINE__ : VECL == 32 && WID == 4 cannot work, too long vector for one plane!"
#elif VECL == 32 && WID == 8
         const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7});
         const Veci j_indices = Veci({j,     j,     j,     j,     j,     j,     j,     j,     j + 1, j + 1, j + 1, j + 1, j + 1, j + 1, j + 1, j + 1,
                                      j + 2, j + 2, j + 2, j + 2, j + 2, j + 2, j + 2, j + 2, j + 3, j + 3, j + 3, j + 3, j + 3, j + 3, j + 3, j + 3});
#elif VECL == 64 && WID == 4
#error "__FILE__ : __LINE__ : VECL == 64 && WID == 4 cannot work, too long vector for one plane!"
#elif VECL == 64 && WID == 8
         const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
                                      0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7});
         const Veci j_indices = Veci({j,     j,     j,     j,     j,     j,     j,     j,     j + 1, j + 1, j + 1, j + 1, j + 1, j + 1, j + 1, j + 1, j + 2, j + 2, j + 2, j + 2, j + 2, j + 2,
                                      j + 2, j + 2, j + 3, j + 3, j + 3, j + 3, j + 3, j + 3, j + 3, j + 3, j + 4, j + 4, j + 4, j + 4, j + 4, j + 4, j + 4, j + 4, j + 5, j + 5, j + 5, j + 5,
                                      j + 5, j + 5, j + 5, j + 5, j + 6, j + 6, j + 6, j + 6, j + 6, j + 6, j + 6, j + 6, j + 7, j + 7, j + 7, j + 7, j + 7, j + 7, j + 7, j + 7});
#else
#error "This VECL && This WID cannot work!"
#define xstr(s) str(s)
#define str(s) #s
#pragma message "VECL =" xstr(VECL)
#pragma message "WID = "xstr(WID)
#endif

         loop_body(i_indices, j_indices, k);
      }
   }
}

void pitchAngleDiffusion(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, const uint popID) {

   // Ensure nu0 dat file is read, if requested
   if (P::PADcoefficient < 0) {
      readNuArrayFromFile();
   }

   int nbins_v = Parameters::PADvbins;
   int nbins_mu = Parameters::PADmubins;
   const Real dmubins = 2.0 / nbins_mu;

   // resonance gap filling coefficient, not needed assuming even number of bins in mu-space
   const Real epsilon = 0.0;

   phiprof::Timer diffusionTimer{"pitch-angle-diffusion"};

   const auto LocalCells = getLocalCells();
   #pragma omp parallel
   {
      std::vector<int> fcount(nbins_v * nbins_mu, 0);                    // Array to count number of f stored
      std::vector<Realf> fmu(nbins_v * nbins_mu, 0);                     // Array to store f(v,mu)
      std::vector<Realf> dfdmu(nbins_v * nbins_mu, 0);                   // Array to store dfdmu
      std::vector<Realf> dfdmu2(nbins_v * nbins_mu, 0);                  // Array to store dfdmumu
      std::vector<Realf> dfdt_mu(nbins_v * nbins_mu, 0);                 // Array to store dfdt_mu
      #pragma omp for
      for (size_t CellIdx = 0; CellIdx < LocalCells.size(); CellIdx++) { // Iterate over all spatial cells

         const auto CellID = LocalCells[CellIdx];
         SpatialCell& cell = *mpiGrid[CellID];
         const Real* parameters = cell.get_block_parameters(popID);
         const vmesh::LocalID* nBlocks = cell.get_velocity_grid_length(popID);
         const size_t meshID = getObjectWrapper().particleSpecies[popID].velocityMesh;
         const vmesh::MeshParameters& vMesh = vmesh::getMeshWrapper()->velocityMeshes->at(meshID);

         // Ensure mass conservation
         Realf density_pre_adjust = 0.0;
         Realf density_post_adjust = 0.0;
         if (getObjectWrapper().particleSpecies[popID].sparse_conserve_mass) {
            Vec vectorSum{0};
            Vec vectorAdd{0};
            for (size_t i = 0; i < cell.get_number_of_velocity_blocks(popID) * WID3 / VECL; ++i) {
               vectorAdd.load(&cell.get_data(popID)[i * VECL]);
               vectorSum += vectorAdd;
               // density_pre_adjust += cell.get_data(popID)[i];
            }
            density_pre_adjust = horizontal_add(vectorSum);
         }

         Real dtTotalDiff = 0.0; // Diffusion time elapsed

         const Real Vmax = 2 * sqrt(3) * vMesh.meshLimits[1];
         const Real dVbins = Vmax / nbins_v;

         const Real bulkVX = cell.parameters[CellParams::VX];
         const Real bulkVY = cell.parameters[CellParams::VY];
         const Real bulkVZ = cell.parameters[CellParams::VZ];

         bool currentSpatialLoopComplete;
         Realf Sparsity;
         std::array<Real, 3> b;
         Real nu0;

         // Compute parameters
         computePitchAngleDiffusionParameters(cell, popID, CellIdx, currentSpatialLoopComplete, Sparsity, b, nu0);

         // Enable nu0 disk output; skip cells where diffusion is not required (or diffusion coefficient is very small).
         cell.parameters[CellParams::NU0] = nu0;
         if (nu0 <= 0.001) {
            continue;
         }

         while (dtTotalDiff < Parameters::dt) { // Substep loop

            const Real RemainT = Parameters::dt - dtTotalDiff; // Remaining time before reaching simulation time step
            Real checkCFL = std::numeric_limits<Real>::max();

            // Initialised at each substep
            std::fill(fmu.begin(), fmu.end(), 0.0);
            std::fill(fcount.begin(), fcount.end(), 0);

            // Build 2d array of f(v,mu)
            for (vmesh::LocalID n = 0; n < cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks

               loop_over_block([&](Veci i_indices, Veci j_indices, int k) -> void { // Lambda function processor
                  // Get velocity space coordinates
                  const Vec VX(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] +
                               (to_realf(i_indices) + 0.5) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]);
                  const Vec VY(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] +
                               (to_realf(j_indices) + 0.5) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]);
                  const Vec VZ(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] + (k + 0.5) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ]);

                  const Vec VplasmaX = VX - bulkVX;
                  const Vec VplasmaY = VY - bulkVY;
                  const Vec VplasmaZ = VZ - bulkVZ;

                  const Vec normV = sqrt(VplasmaX * VplasmaX + VplasmaY * VplasmaY + VplasmaZ * VplasmaZ);
                  const Vec Vpara = VplasmaX * b[0] + VplasmaY * b[1] + VplasmaZ * b[2];
                  const Vec mu = Vpara / (normV + std::numeric_limits<Real>::min()); // + min value to avoid division by 0.

                  const Veci Vindex = roundi(floor((normV) / dVbins));
                  const Vec Vmu = dVbins * (to_realf(Vindex) + 0.5); // Take value at the center of the mu cell
                  Veci muindex = roundi(floor((mu + 1.0) / dmubins));

                  Vec CellValue;
                  CellValue.load(&cell.get_data(n, popID)[WID2 * k + WID * j_indices[0] + i_indices[0]]);
                  const Vec increment = 2.0 * M_PI * Vmu * Vmu * CellValue;
                  for (uint i = 0; i < VECL; i++) {
                     // Safety check to handle edge case where mu = exactly 1.0
                     const int mui = std::max(0, std::min((int)muindex[i], nbins_mu - 1));
                     const int vi = std::max(0, std::min((int)Vindex[i], nbins_v - 1));
                     MUSPACE(fmu, vi, mui) += increment[i];
                     MUSPACE(fcount, vi, mui) += 1;
                  }
               }); // End of Lambda
            } // End blocks

            int cRight;
            int cLeft;

            // Compute space/time derivatives (take first non-zero neighbours) & CFL & Ddt
            for (int indv = 0; indv < nbins_v; indv++) {
               const Real Vmu = dVbins * (float(indv) + 0.5);

               // Divide f by count (independent of v but needs to be computed for all mu before derivatives)
               for (int indmu = 0; indmu < nbins_mu; indmu++) {
                  if (MUSPACE(fcount, indv, indmu) == 0 || MUSPACE(fmu, indv, indmu) <= 0.0) {
                     MUSPACE(fmu, indv, indmu) = 0;
                  } else {
                     MUSPACE(fmu, indv, indmu) = MUSPACE(fmu, indv, indmu) / MUSPACE(fcount, indv, indmu);
                  }
               }

               // Search limits for how many cells in mu-direction should be max evaluated when searching for a near neighbour?
               // Assuming some oversampling; changing these values may result in method breaking at very small plasma frame velocities.
               const int rlimit = nbins_mu - 1;
               const int llimit = 0;

               for (int indmu = 0; indmu < nbins_mu; indmu++) {
                  // Compute spatial derivatives
                  if (indmu == 0) {
                     cLeft = 0;
                     cRight = 1;
                     while ((MUSPACE(fcount, indv, indmu + cRight) == 0) && (indmu + cRight < rlimit)) {
                        cRight += 1;
                     }
                     if ((MUSPACE(fcount, indv, indmu + cRight) == 0) && (indmu + cRight == rlimit)) {
                        cRight = 0;
                     }
                  } else if (indmu == nbins_mu - 1) {
                     cLeft = 1;
                     cRight = 0;
                     while ((MUSPACE(fcount, indv, indmu - cLeft) == 0) && (indmu - cLeft > llimit)) {
                        cLeft += 1;
                     }
                     if ((MUSPACE(fcount, indv, indmu - cLeft) == 0) && (indmu - cLeft == llimit)) {
                        cLeft = 0;
                     }
                  } else {
                     cLeft = 1;
                     cRight = 1;
                     while ((MUSPACE(fcount, indv, indmu + cRight) == 0) && (indmu + cRight < rlimit)) {
                        cRight += 1;
                     }
                     if ((MUSPACE(fcount, indv, indmu + cRight) == 0) && (indmu + cRight == rlimit)) {
                        cRight = 0;
                     }
                     while ((MUSPACE(fcount, indv, indmu - cLeft) == 0) && (indmu - cLeft > llimit)) {
                        cLeft += 1;
                     }
                     if ((MUSPACE(fcount, indv, indmu - cLeft) == 0) && (indmu - cLeft == llimit)) {
                        cLeft = 0;
                     }
                  }
                  if ((cRight == 0) && (cLeft != 0)) {
                     MUSPACE(dfdmu, indv, indmu) = (MUSPACE(fmu, indv, indmu + cRight) - MUSPACE(fmu, indv, indmu - cLeft)) / ((cRight + cLeft) * dmubins);
                     MUSPACE(dfdmu2, indv, indmu) = 0.0;
                  } else if ((cLeft == 0) && (cRight != 0)) {
                     MUSPACE(dfdmu, indv, indmu) = (MUSPACE(fmu, indv, indmu + cRight) - MUSPACE(fmu, indv, indmu - cLeft)) / ((cRight + cLeft) * dmubins);
                     MUSPACE(dfdmu2, indv, indmu) = 0.0;
                  } else if ((cLeft == 0) && (cRight == 0)) {
                     MUSPACE(dfdmu, indv, indmu) = 0.0;
                     MUSPACE(dfdmu2, indv, indmu) = 0.0;
                  } else {
                     MUSPACE(dfdmu, indv, indmu) = (MUSPACE(fmu, indv, indmu + cRight) - MUSPACE(fmu, indv, indmu - cLeft)) / ((cRight + cLeft) * dmubins);
                     MUSPACE(dfdmu2, indv, indmu) =
                         ((MUSPACE(fmu, indv, indmu + cRight) - MUSPACE(fmu, indv, indmu)) / (cRight * dmubins) - (MUSPACE(fmu, indv, indmu) - MUSPACE(fmu, indv, indmu - cLeft)) / (cLeft * dmubins)) /
                         (0.5 * dmubins * (cRight + cLeft));
                  }

                  // Compute time derivative
                  const Realf mu = (indmu + 0.5) * dmubins - 1.0;
                  const Realf Dmumu = nu0 / 2.0 * (abs(mu) / (1.0 + abs(mu)) + epsilon) * (1.0 - mu * mu);
                  const Realf dDmu = nu0 / 2.0 * ((mu / abs(mu)) * ((1.0 - mu * mu) / ((1.0 + abs(mu)) * (1.0 + abs(mu)))) - 2.0 * mu * (abs(mu) / (1.0 + abs(mu)) + epsilon));
                  // We divide dfdt_mu by the normalization factor 2pi*v^2 already here.
                  const Realf dfdt_mu_val = (dDmu * MUSPACE(dfdmu, indv, indmu) + Dmumu * MUSPACE(dfdmu2, indv, indmu)) / (2.0 * M_PI * Vmu * Vmu);
                  MUSPACE(dfdt_mu, indv, indmu) = dfdt_mu_val;

                  // Only consider CFL for non-negative phase-space cells above the sparsity threshold
                  const Realf CellValue = MUSPACE(fmu, indv, indmu) / (2.0 * M_PI * Vmu * Vmu);
                  const Realf absdfdt = abs(MUSPACE(dfdt_mu, indv, indmu)); // Already scaled
                  if (absdfdt > 0.0 && CellValue > Sparsity) {
                     checkCFL = std::min(CellValue * Parameters::PADCFL * (1.0 / absdfdt), checkCFL);
                  }
               } // End mu loop
            } // End v loop

            // Compute Ddt
            Real Ddt = checkCFL;
            if (Ddt > RemainT) {
               Ddt = RemainT;
            }
            dtTotalDiff = dtTotalDiff + Ddt;

            for (vmesh::LocalID n = 0; n < cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks

               loop_over_block([&](Veci i_indices, Veci j_indices, int k) -> void { // Lambda function processor
                  // Get velocity space coordinates
                  const Vec VX(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] +
                               (to_realf(i_indices) + 0.5) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]);
                  const Vec VY(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] +
                               (to_realf(j_indices) + 0.5) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]);
                  const Vec VZ(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] + (k + 0.5) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ]);

                  const Vec VplasmaX = VX - bulkVX;
                  const Vec VplasmaY = VY - bulkVY;
                  const Vec VplasmaZ = VZ - bulkVZ;

                  const Vec normV = sqrt(VplasmaX * VplasmaX + VplasmaY * VplasmaY + VplasmaZ * VplasmaZ);
                  const Vec Vpara = VplasmaX * b[0] + VplasmaY * b[1] + VplasmaZ * b[2];
                  const Vec mu = Vpara / (normV + std::numeric_limits<Real>::min()); // + min value to avoid division by 0.

                  const Veci Vindex = roundi(floor((normV) / dVbins));
                  const Vec Vmu = dVbins * (to_realf(Vindex) + 0.5); // Take value at the center of the mu cell
                  Veci muindex = roundi(floor((mu + 1.0) / dmubins));

                  // Compute dfdt
                  std::array<Realf, VECL> dfdt = {0};
                  for (uint i = 0; i < VECL; i++) {
                     // Safety check to handle edge case where mu = exactly 1.0
                     const int mui = std::max(0, std::min((int)muindex[i], nbins_mu - 1));
                     const int vi = std::max(0, std::min((int)Vindex[i], nbins_v - 1));
                     dfdt[i] = MUSPACE(dfdt_mu, vi, mui); // dfdt_mu was scaled back down by 2pi*v^2 on creation
                  }
                  Vec dfdtUpdate;
                  dfdtUpdate.load(&dfdt[0]);

                  // Update cell value, ensuring result is non-negative
                  Vec CellValue;
                  CellValue.load(&cell.get_data(n, popID)[WID2 * k + WID * j_indices[0] + i_indices[0]]);
                  Vec NewCellValue = CellValue + dfdtUpdate * Ddt;
                  const Vecb lessZero = NewCellValue < 0.0;
                  NewCellValue = select(lessZero, 0.0, NewCellValue);
                  NewCellValue.store(&cell.get_data(n, popID)[WID2 * k + WID * j_indices[0] + i_indices[0]]);
               }); // End of Lambda
            } // End Blocks

         } // End Time loop

         // Ensure mass conservation
         if (getObjectWrapper().particleSpecies[popID].sparse_conserve_mass) {
            Vec vectorSum{0};
            Vec vectorAdd{0};
            for (size_t i = 0; i < cell.get_number_of_velocity_blocks(popID) * WID3 / VECL; ++i) {
               vectorAdd.load(&cell.get_data(popID)[i * VECL]);
               vectorSum += vectorAdd;
            }
            density_post_adjust = horizontal_add(vectorSum);

            if (density_post_adjust != 0.0 && density_pre_adjust != density_post_adjust) {
               const Vec adjustRatio = density_pre_adjust / density_post_adjust;
               Vec vectorAdjust;
               for (size_t i = 0; i < cell.get_number_of_velocity_blocks(popID) * WID3 / VECL; ++i) {
                  vectorAdjust.load(&cell.get_data(popID)[i * VECL]);
                  vectorAdjust *= adjustRatio;
                  vectorAdjust.store(&cell.get_data(popID)[i * VECL]);
               }
            }
         }
      } // End spatial cell loop
   } // End parallel workshare region
} // End function
