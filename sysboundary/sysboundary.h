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

#ifndef SYSBOUNDARY_H
#define SYSBOUNDARY_H

#include <map>
#include <list>
#include <vector>
#include <mpi.h>
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "../definitions.h"
#include "../parameters.h"
#include "../readparameters.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"

#include "sysboundarycondition.h"

/*! \brief SysBoundary contains the SysBoundaryConditions used in the simulation.
 *
 * The purpose of SysBoundary is to contain SBC::SysBoundaryConditions, and apply
 * them to the simulation volume cells if the cells pertain to a specific boundary type.
 * If the simulation domain is not fully periodic then the behaviour at the edges or boundaries of the volume has to be properly defined.
 *
 * initSysBoundaries creates the instances of SBC::SysBoundaryConditions that are needed.
 * They are then initialised, which means the internals are prepared for the system
 * boundary condition to be applied (import and process parameters, generate template cells
 * etc.). When the whole simulation domain is initialised, the boundary conditions are
 * applied to the cells they have by calling applyInitialState.
 *
 * If needed, a user can write his or her own SBC::SysBoundaryConditions, which
 * are loaded when the simulation initializes.
 */
class SysBoundary {
 public:
   SysBoundary();
   ~SysBoundary();

   void addParameters();
   void getParameters();

   void addSysBoundary(
                       SBC::SysBoundaryCondition* sbc,
                       Project& project,
                       creal& t
                      );
   void initSysBoundaries(
                          Project& project,
                          creal& t
                         );
   bool existSysBoundary(std::string name);
   void checkRefinement(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
   void classifyCells(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                     FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid);
   void applyInitialState(
                          dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
                          FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                          FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                          Project& project
                         );
   void updateState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                    FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
                    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                    FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                    creal t);
   void applySysBoundaryVlasovConditions(dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, creal& t, const bool calculate_V_moments);
   void setupL2OutflowAtRestart(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid);

   unsigned int size() const;
   SBC::SysBoundaryCondition* getSysBoundary(cuint sysBoundaryType) const;
   bool isAnyDynamic() const;
   bool isPeriodic(uint direction) const;
   void updateSysBoundariesAfterLoadBalance(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid);
   void clear() { // Clears all conts of SBC (destructing template cells for GPU branch)
      sysBoundaries.clear();
   }
   private:
      /*! Private copy-constructor to prevent copying the class. */
      SysBoundary(const SysBoundary& bc);

      //std::set<SBC::SysBoundaryCondition*,SBC::Comparator> sysBoundaries;

      /*! A container for all SBC::SysBoundaryConditions stored in SysBoundary.*/
      std::list<SBC::SysBoundaryCondition*> sysBoundaries;
      /*! A map from the system boundary types to the corresponding class member. */
      std::map<uint, SBC::SysBoundaryCondition*> indexToSysBoundary;
      /*! List of system boundary conditions (SBC) to be used. */
      std::vector<std::string> sysBoundaryCondList;
      /*! bool telling whether any system boundary condition is dynamic in time (and thus needs updating). */
      bool anyDynamic;

      /*! Array of bool telling whether the system is periodic in any direction. */
      bool periodic[3];
};

bool precedenceSort(const SBC::SysBoundaryCondition* first,
                    const SBC::SysBoundaryCondition* second);



/*
   Input a vector of cellIDs (cellList) and compute a new vector with only those cells which are on a sysboundary and are to be computed
*/
void getBoundaryCellList(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                         const std::vector<CellID>& cellList,
                         std::vector<CellID>& boundaryCellList);

#endif
