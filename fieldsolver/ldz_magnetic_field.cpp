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
 */

#ifdef _OPENMP
   #include <omp.h>
#endif

#include "fs_cache.h"
#include "ldz_magnetic_field.hpp"

using namespace std;

/*! \brief Low-level magnetic field propagation function.
 * 
 * Propagates the cell's face-averaged magnetic field components by
 * using Faraday's law on the face edges. Depending on the time order
 * of accuracy it is done in one stage or in two stages using the
 * intermediate E1 components for the first stage of the second-order
 * Runge-Kutta method and E for the other cases.
 * 
 * \param cellCache Field solver cell cache
 * \param cells Vector of cells to process
 * \param dt Length of the time step
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param doX If true, compute the x component (default true).
 * \param doY If true, compute the y component (default true).
 * \param doZ If true, compute the z component (default true).
 */
void propagateMagneticField(
   const std::vector<fs_cache::CellCache>& cellCache,
   const std::vector<uint16_t>& cells,
   creal& dt,
   cint& RKCase,
   const bool doX, //=true (default)
   const bool doY, //=true (default)
   const bool doZ  //=true (default)
) {

   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const uint16_t localID = cells[c];

      cuint existingCellsFlag = cellCache[localID].existingCellsFlags;
      Real* cp0 = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
      creal dx = cp0[CellParams::DX];
      creal dy = cp0[CellParams::DY];
      creal dz = cp0[CellParams::DZ];

      if ((existingCellsFlag & PROPAGATE_BX) == PROPAGATE_BX && doX == true) {
         const Real* cp1 = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->parameters;
         const Real* cp2 = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->parameters;
         switch (RKCase) {
          case RK_ORDER1:
            cp0[CellParams::PERBX] += dt/dz*(cp2[CellParams::EY] - cp0[CellParams::EY]) + dt/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ]);
            break;
          case RK_ORDER2_STEP1:
            cp0[CellParams::PERBX_DT2] =
              cp0[CellParams::PERBX] + 0.5*dt*(1.0/dz*(cp2[CellParams::EY] - cp0[CellParams::EY]) +
                                               1.0/dy*(cp0[CellParams::EZ] - cp1[CellParams::EZ]));
            break;
          case RK_ORDER2_STEP2:
            cp0[CellParams::PERBX] += dt * (1.0/dz*(cp2[CellParams::EY_DT2] - cp0[CellParams::EY_DT2]) +
                                            1.0/dy*(cp0[CellParams::EZ_DT2] - cp1[CellParams::EZ_DT2]));
            break;
          default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "Invalid RK case." << std::endl;
            abort();
         }
      }

      if ((existingCellsFlag & PROPAGATE_BY) == PROPAGATE_BY && doY == true) {
         const Real* cp1 = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->parameters;
         const Real* cp2 = cellCache[localID].cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->parameters;

         switch (RKCase) {
          case RK_ORDER1:
            cp0[CellParams::PERBY] += dt/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) + dt/dz*(cp0[CellParams::EX] - cp1[CellParams::EX]);
            break;
          case RK_ORDER2_STEP1:
            cp0[CellParams::PERBY_DT2] =
              cp0[CellParams::PERBY] + 0.5*dt*(1.0/dx*(cp2[CellParams::EZ] - cp0[CellParams::EZ]) +
                                               1.0/dz*(cp0[CellParams::EX] - cp1[CellParams::EX]));
            break;
          case RK_ORDER2_STEP2:
                  cp0[CellParams::PERBY] += dt * (1.0/dx*(cp2[CellParams::EZ_DT2] - cp0[CellParams::EZ_DT2]) +
                                                  1.0/dz*(cp0[CellParams::EX_DT2] - cp1[CellParams::EX_DT2]));
            break;
          default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "Invalid RK case." << std::endl;
            abort();
         }
      }

      if ((existingCellsFlag & PROPAGATE_BZ) == PROPAGATE_BZ && doZ == true) {
         const Real* cp1 = cellCache[localID].cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->parameters;
         const Real* cp2 = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->parameters;

         switch (RKCase) {
          case RK_ORDER1:
            cp0[CellParams::PERBZ] += dt/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) + dt/dx*(cp0[CellParams::EY] - cp1[CellParams::EY]);
            break;
          case RK_ORDER2_STEP1:
            cp0[CellParams::PERBZ_DT2] =
              cp0[CellParams::PERBZ] + 0.5*dt*(1.0/dy*(cp2[CellParams::EX] - cp0[CellParams::EX]) +
                                               1.0/dx*(cp0[CellParams::EY] - cp1[CellParams::EY]));
            break;
          case RK_ORDER2_STEP2:
            cp0[CellParams::PERBZ] += dt  * (1.0/dy*(cp2[CellParams::EX_DT2] - cp0[CellParams::EX_DT2]) +
                                             1.0/dx*(cp0[CellParams::EY_DT2] - cp1[CellParams::EY_DT2]));
            break;
          default:
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << "Invalid RK case." << std::endl;
            abort();
         }
      }
   }
}

/*! \brief High-level magnetic field propagation function.
 * 
 * Propagates the magnetic field and applies the field boundary conditions defined in project.h where needed.
 * 
 * \param mpiGrid Grid
 * \param sysBoundaries System boundary conditions existing
 * \param dt Length of the time step
 * \param localCells Vector of local cells to process
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa propagateMagneticField propagateSysBoundaryMagneticField
 */
void propagateMagneticFieldSimple(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   const std::vector<CellID>& localCells,
   cint& RKCase
) {

   phiprof::start("Propagate magnetic field");
   
   fs_cache::CacheContainer& cacheContainer = fs_cache::getCache();

   int timer=phiprof::initializeTimer("Compute system inner cells");
   phiprof::start(timer);

   // Propagate B on all local cells:
   propagateMagneticField(cacheContainer.localCellsCache,cacheContainer.local_NOT_SYSBOUND_DO_NOT_COMPUTE,dt,RKCase);

   //phiprof::stop("propagate not sysbound",localCells.size(),"Spatial Cells");
   phiprof::stop(timer,cacheContainer.local_NOT_SYSBOUND_DO_NOT_COMPUTE.size(),"Spatial Cells");

   //This communication is needed for boundary conditions, in practice almost all
   //of the communication is going to be redone in calculateDerivativesSimple
   //TODO: do not transfer if there are no field boundaryconditions
   phiprof::start("MPI");
   if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
      // Exchange PERBX,PERBY,PERBZ with neighbours
      spatial_cell::SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERB,true);
   } else { // RKCase == RK_ORDER2_STEP1
      // Exchange PERBX_DT2,PERBY_DT2,PERBZ_DT2 with neighbours
      spatial_cell::SpatialCell::set_mpi_transfer_type(Transfer::CELL_PERBDT2,true);
   }

   mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID);
   phiprof::stop("MPI");

   // Propagate B on system boundary/process inner cells
   timer=phiprof::initializeTimer("Compute system boundary/process inner cells");
   phiprof::start(timer);
   #pragma omp parallel for
   for (size_t c=0; c<cacheContainer.boundaryCellsWithLocalNeighbours.size(); ++c) {
      const uint16_t localID = cacheContainer.boundaryCellsWithLocalNeighbours[c];
      propagateSysBoundaryMagneticField(mpiGrid, cacheContainer.localCellsCache, localID, sysBoundaries, dt, RKCase);
   }
   phiprof::stop(timer,cacheContainer.boundaryCellsWithLocalNeighbours.size(),"Spatial Cells");

   // Propagate B on system boundary/process boundary cells
   timer=phiprof::initializeTimer("Compute system boundary/process boundary cells");
   phiprof::start(timer);
   #pragma omp parallel for
   for (size_t c=0; c<cacheContainer.boundaryCellsWithRemoteNeighbours.size(); ++c) {
      const uint16_t localID = cacheContainer.boundaryCellsWithRemoteNeighbours[c];
      propagateSysBoundaryMagneticField(mpiGrid, cacheContainer.localCellsCache, localID, sysBoundaries, dt, RKCase);
   }
   phiprof::stop(timer,cacheContainer.boundaryCellsWithRemoteNeighbours.size(),"Spatial Cells");

   const size_t N_cells 
     = cacheContainer.boundaryCellsWithRemoteNeighbours.size()
     + cacheContainer.boundaryCellsWithLocalNeighbours.size()
     + cacheContainer.local_NOT_SYSBOUND_DO_NOT_COMPUTE.size();
   phiprof::stop("Propagate magnetic field",N_cells,"Spatial Cells");
}

/*! \brief Low-level magnetic field propagation function.
 * 
 * Propagates the magnetic field according to the system boundary conditions.
 * 
 * \param mpiGrid Grid
 * \param cellCache Field solver cell cache
 * \param localID Field solver cache local cell ID
 * \param sysBoundaries System boundary conditions existing
 * \param dt Length of the time step
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa propagateMagneticFieldSimple propagateMagneticField
 */
void propagateSysBoundaryMagneticField(
   const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<fs_cache::CellCache>& cellCache,
   const uint16_t& localID,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& RKCase
) {
   Real* cp0 = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
   cuint sysBoundaryFlag = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->sysBoundaryFlag;
   int offset = 0;
   for (uint component = 0; component < 3; component++) {
      if (RKCase == RK_ORDER2_STEP1) {
         offset = CellParams::PERBX_DT2 - CellParams::PERBX;
      }
      cp0[CellParams::PERBX + offset + component] =
        sysBoundaries.getSysBoundary(sysBoundaryFlag)->
        fieldSolverBoundaryCondMagneticField(mpiGrid, cellCache, localID, dt, RKCase, offset, component);
   }
}
