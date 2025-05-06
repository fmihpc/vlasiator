/*
 * This file is part of Vlasiator.
 * Copyright 2024-2025 University of Helsinki, CSC
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

#ifndef ARCH_MOMENTS_H
#define ARCH_MOMENTS_H

#include <vector>
#include <limits>
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include "../definitions.h"
#include "../common.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"

using namespace spatial_cell;

#define nMom1 4
#define nMom2 6

// ***** FUNCTION DECLARATIONS ***** //

template<typename REAL, uint SIZE>
void blockVelocityFirstMoments(vmesh::VelocityBlockContainer *blockContainer,
                               REAL (&array)[SIZE],
                               uint nBlocks);

template<typename REAL, uint SIZE>
void blockVelocitySecondMoments(vmesh::VelocityBlockContainer *blockContainer,
                                const REAL averageVX,
                                const REAL averageVY,
                                const REAL averageVZ,
                                REAL (&array)[SIZE],
                                uint nBlocks);

void calculateMoments_R(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells,
   const bool& computeSecond,
   const bool initialCompute=false
);

void calculateMoments_V(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID>& cells,
   const bool& computeSecond,
   const bool initialCompute=false
);


// ***** TEMPLATE FUNCTION DEFINITIONS ***** //

/** Calculate the zeroth and first velocity moments for the given
 * velocity block and add results to 'array', which must have at
 * least size four. After this function returns, the contents of
 * 'array' are as follows: array[0]=n; array[1]=n*Vx; array[2]=nVy;
 * array[3]=nVz; Here n is the scaled number density, i.e., number density
 * times population mass / proton mass. This function is AMR safe.
 * @param data Distribution functions
 * @param blockParameters Parameters for the velocity blocks
 * @param array Array where the calculated moments are added
 * @param nBlocks The number of blocks */
template<typename REAL, uint SIZE> inline
void blockVelocityFirstMoments(
   vmesh::VelocityBlockContainer *blockContainer,
   REAL (&array)[SIZE],
   uint nBlocks) {

   arch::parallel_reduce<arch::sum>({WID, WID, WID, nBlocks},
     ARCH_LOOP_LAMBDA (const uint i, const uint j, const uint k, const uint blockLID, Real *lsum ) {

       Realf *data = blockContainer->getData();
       Real *blockParameters = blockContainer->getParameters();
       const Realf* avgs = &data[blockLID*WID3];
       const Real* blockParamsZ = &blockParameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS];
       const Real DV3 = blockParamsZ[BlockParams::DVX]*blockParamsZ[BlockParams::DVY]*blockParamsZ[BlockParams::DVZ];
       const Real HALF = 0.5;

       ARCH_INNER_BODY(i, j, k, blockLID, lsum) {
         const Real VX = blockParamsZ[BlockParams::VXCRD] + (i+HALF)*blockParamsZ[BlockParams::DVX];
         const Real VY = blockParamsZ[BlockParams::VYCRD] + (j+HALF)*blockParamsZ[BlockParams::DVY];
         const Real VZ = blockParamsZ[BlockParams::VZCRD] + (k+HALF)*blockParamsZ[BlockParams::DVZ];
         lsum[0] += avgs[cellIndex(i,j,k)] * DV3;
         lsum[1] += avgs[cellIndex(i,j,k)]*VX * DV3;
         lsum[2] += avgs[cellIndex(i,j,k)]*VY * DV3;
         lsum[3] += avgs[cellIndex(i,j,k)]*VZ * DV3;
       };
     }, array);
}

/** Calculate the second velocity moments for the velocity blocks, and add
 * results to 'array', which must have at least size three. After this function
 * returns, the contents of 'array' are as follows: array[0]=n(Vx-Vx0);
 * array[1]=n(Vy-Vy0); array[2]=n(Vz-Vz0); Here Vx0,Vy0,Vz0 are the components
 * of the bulk velocity (calculated over all species). This function is AMR safe.
 * @param data Distribution functions
 * @param blockParameters Parameters for the velocity blocks
 * @param averageVX Bulk velocity x
 * @param averageVY Bulk velocity y
 * @param averageVZ Bulk velocity z
 * @param array Array where the calculated moments are added
 * @param nBlocks The number of blocks */
template<typename REAL, uint SIZE> inline
void blockVelocitySecondMoments(
   vmesh::VelocityBlockContainer *blockContainer,
   const REAL averageVX,
   const REAL averageVY,
   const REAL averageVZ,
   REAL (&array)[SIZE],
   uint nBlocks) {

   arch::parallel_reduce<arch::sum>({WID, WID, WID, nBlocks},
     ARCH_LOOP_LAMBDA (const uint i, const uint j, const uint k, const uint blockLID, Real *lsum ) {

       Realf *data = blockContainer->getData();
       Real *blockParameters = blockContainer->getParameters();
       const Realf* avgs = &data[blockLID*WID3];
       const Real* blockParams = &blockParameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS];
       const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
       const Real HALF = 0.5;

       ARCH_INNER_BODY(i, j, k, blockLID, lsum) {
         const Real VX = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX];
         const Real VY = blockParams[BlockParams::VYCRD] + (j+HALF)*blockParams[BlockParams::DVY];
         const Real VZ = blockParams[BlockParams::VZCRD] + (k+HALF)*blockParams[BlockParams::DVZ];
         lsum[0] += avgs[cellIndex(i,j,k)] * (VX - averageVX) * (VX - averageVX) * DV3;
         lsum[1] += avgs[cellIndex(i,j,k)] * (VY - averageVY) * (VY - averageVY) * DV3;
         lsum[2] += avgs[cellIndex(i,j,k)] * (VZ - averageVZ) * (VZ - averageVZ) * DV3;
         lsum[3] += avgs[cellIndex(i,j,k)] * (VY - averageVY) * (VZ - averageVZ) * DV3;
         lsum[4] += avgs[cellIndex(i,j,k)] * (VX - averageVX) * (VZ - averageVZ) * DV3;
         lsum[5] += avgs[cellIndex(i,j,k)] * (VX - averageVX) * (VY - averageVY) * DV3;
       };
     }, array);
}

#endif
