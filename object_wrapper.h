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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef OBJECT_WRAPPER_H
#define OBJECT_WRAPPER_H

#include <vector>
#include <fsgrid.hpp>

#include "definitions.h"
#include "item_storage.h"
#include "object_factory.h"
#include "vamr_refinement_criteria.h"
#include "particle_species.h"
#include "projects/project.h"
#include "velocity_mesh_parameters.h"
#include "sysboundary/sysboundarycondition.h"
#include "sysboundary/sysboundary.h"
#include "common.h"


struct ObjectWrapper {
   ObjectWrapper() { }

   ObjectFactory<vamr_ref_criteria::Base> vamrVelRefCriteria; /**< Factory for all known VAMR refinement criteria.*/
   std::vector<species::Species> particleSpecies;           /**< Parameters for all particle species.*/
   projects::Project*                    project;           /**< Simulated project.*/
   std::vector<vmesh::MeshParameters> velocityMeshes;       /**< Parameters for velocity mesh(es).*/
   SysBoundary sysBoundaryContainer;                        /**< Container for sysboundaries.*/

   bool addParameters();                                    /**< Add config file parameters for objects held in this wrapper */
   bool addPopulationParameters();                          /**< After parsing the names of populations, create parameters for each of them */
   bool getParameters();                                    /**< Use parsed config file parameters for objects held in this wrapper */

 private:
   ObjectWrapper(const ObjectWrapper& ow);
   ObjectWrapper& operator=(const ObjectWrapper& ow);
};

// Currently defined in vlasiator.cpp
ObjectWrapper& getObjectWrapper();

// All fsgrids together in a crunchy wrapper
struct FsGridWrapper {

      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> perBGrid;
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> perBDt2Grid;
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> EGrid;
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> EDt2Grid;
      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> EHallGrid;
      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> EGradPeGrid;
      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> EGradPeDt2Grid;
      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> momentsGrid;
      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> momentsDt2Grid;
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> dPerBGrid;
      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> dMomentsGrid;
      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> dMomentsDt2Grid;
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> BgBGrid;
      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> volGrid;
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> technicalGrid;

      FsGridWrapper(std::array<FsGridTools::FsSize_t,3> fsGridDimensions, 
            MPI_Comm communicator, std::array<bool, 3> periodicity,
            std::array<Real, 3> DX, std::array<Real, 3> xmin) :
      perBGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      perBDt2Grid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      EGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      EDt2Grid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      EHallGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      EGradPeGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      EGradPeDt2Grid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      momentsGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      momentsDt2Grid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      dPerBGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      dMomentsGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      dMomentsDt2Grid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      BgBGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      volGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition),
      technicalGrid(fsGridDimensions, communicator, periodicity,P::manualFsGridDecomposition) {

         // Set DX, DY and DZ
         // TODO: This is currently just taking the values from cell 1, and assuming them to be
         // constant throughout the simulation.
         perBGrid.DX = perBDt2Grid.DX = EGrid.DX = EDt2Grid.DX = EHallGrid.DX = EGradPeGrid.DX = EGradPeDt2Grid.DX = momentsGrid.DX
            = momentsDt2Grid.DX = dPerBGrid.DX = dMomentsGrid.DX = dMomentsDt2Grid.DX = BgBGrid.DX = volGrid.DX = technicalGrid.DX
            = DX[0];

         perBGrid.DY = perBDt2Grid.DY = EGrid.DY = EDt2Grid.DY = EHallGrid.DY = EGradPeGrid.DY = EGradPeDt2Grid.DY = momentsGrid.DY
            = momentsDt2Grid.DY = dPerBGrid.DY = dMomentsGrid.DY = dMomentsDt2Grid.DY = BgBGrid.DY = volGrid.DY = technicalGrid.DY
            = DX[1];

         perBGrid.DZ = perBDt2Grid.DZ = EGrid.DZ = EDt2Grid.DZ = EHallGrid.DZ = EGradPeGrid.DZ = EGradPeDt2Grid.DZ = momentsGrid.DZ
            = momentsDt2Grid.DZ = dPerBGrid.DZ = dMomentsGrid.DZ = dMomentsDt2Grid.DZ = BgBGrid.DZ = volGrid.DZ = technicalGrid.DZ
            = DX[2];

         // Set the physical start (lower left corner) X, Y, Z
         perBGrid.physicalGlobalStart = perBDt2Grid.physicalGlobalStart = EGrid.physicalGlobalStart = EDt2Grid.physicalGlobalStart
            = EHallGrid.physicalGlobalStart = EGradPeGrid.physicalGlobalStart = EGradPeDt2Grid.physicalGlobalStart = momentsGrid.physicalGlobalStart
            = momentsDt2Grid.physicalGlobalStart = dPerBGrid.physicalGlobalStart = dMomentsGrid.physicalGlobalStart = dMomentsDt2Grid.physicalGlobalStart
            = BgBGrid.physicalGlobalStart = volGrid.physicalGlobalStart = technicalGrid.physicalGlobalStart
            = xmin;

      }
};

#endif
