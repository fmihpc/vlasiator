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

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "../parameters.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <list>
#include <map>
#include <mpi.h>
#include <vector>

#include "boundarycondition.h"

/*! \brief Boundary contains the BoundaryConditions used in the simulation.
 *
 * The purpose of Boundary is to contain BC::BoundaryConditions, and
 * apply them to the simulation volume cells if the cells pertain to a specific
 * boundary type. If the simulation domain is not fully periodic then the
 * behaviour at the edges or boundaries of the volume has to be properly
 * defined.
 *
 * initBoundaries creates the instances of BC::BoundaryConditions that
 * are needed. They are then initialised, which means the internals are prepared
 * for the boundary condition to be applied (import and process
 * parameters, generate template cells etc.). When the whole simulation domain
 * is initialised, the boundary conditions are applied to the cells they have by
 * calling applyInitialState.
 *
 * If needed, a user can write his or her own BC::BoundaryConditions, which
 * are loaded when the simulation initializes.
 */
class Boundary {
public:
   Boundary();
   ~Boundary();

   void addParameters();
   void getParameters();

   void addBoundary(BC::BoundaryCondition *sbc, Project &project, creal t);
   void initBoundaries(Project &project, creal t);
   void checkRefinement(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid);
   void classifyCells(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                      FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid);
   void applyInitialState(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                          FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                          Project &project);
   void updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid, creal t);
   void applyBoundaryVlasovConditions(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                      creal t, const bool doCalcMomentsV);
   unsigned int size() const;
   BC::BoundaryCondition *getBoundary(cuint boundaryType) const;
   bool isDynamic() const;
   bool isPeriodic(uint direction) const;
   void updateBoundariesAfterLoadBalance(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid);

private:
   /*! Private copy-constructor to prevent copying the class. */
   Boundary(const Boundary &bc);

   // std::set<BC::BoundaryCondition*,BC::Comparator> boundaries;

   /*! A container for all BC::BoundaryConditions stored in Boundary.*/
   std::list<BC::BoundaryCondition *> boundaries;
   /*! A map from the boundary types to the corresponding class member. */
   std::map<uint, BC::BoundaryCondition *> indexToBoundary;
   /*! List of boundary conditions (BC) to be used. */
   std::vector<std::string> boundaryCondList;
   /*! bool telling whether any boundary condition is dynamic in time (and thus needs updating). */
   bool anyDynamic;

   /*! Array of bool telling whether the system is periodic in any direction. */
   bool periodic[3];
};

bool precedenceSort(const BC::BoundaryCondition *first, const BC::BoundaryCondition *second);

/*
   Input a vector of cellIDs (cellList) and compute a new vector with only those cells which are on a boundary and are
   to be computed
*/
void getBoundaryCellList(const dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                         const std::vector<CellID> &cellList, std::vector<CellID> &boundaryCellList);

#endif
