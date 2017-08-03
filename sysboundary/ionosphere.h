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

#ifndef IONOSPHERE_H
#define IONOSPHERE_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

using namespace projects;
using namespace std;

namespace SBC {

   struct IonosphereSpeciesParameters {
      Real rho;
      Real V0[3];
      Real T;
      uint nSpaceSamples;
      uint nVelocitySamples;
   };

   /*!\brief Ionosphere is a class applying ionospheric boundary conditions.
    * 
    * Ionosphere is a class handling cells tagged as sysboundarytype::IONOSPHERE by this system boundary condition. It applies ionospheric boundary conditions.
    * 
    * These consist in:
    * - Do nothing for the distribution (keep the initial state constant in time);
    * - Keep only the normal perturbed B component and null out the other perturbed components (perfect conductor behavior);
    * - Null out the electric fields.
    */
   class Ionosphere: public SysBoundaryCondition {
   public:
      Ionosphere();
      virtual ~Ionosphere();
      
      static void addParameters();
      virtual void getParameters();
      
      virtual bool initSysBoundary(
         creal& t,
         Project &project
      );
      virtual bool assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
      virtual bool applyInitialState(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         Project &project
      );
      virtual Real fieldSolverBoundaryCondMagneticField(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const std::vector<fs_cache::CellCache>& cellCache,
         const uint16_t& localID,
         creal& dt,
         cuint& RKCase,
         cint& offset,
         cuint& component
      );
      virtual void fieldSolverBoundaryCondElectricField(
         dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         cuint RKCase,
         cuint component
      );
      virtual void fieldSolverBoundaryCondHallElectricField(
         fs_cache::CellCache& cache,
         cuint RKCase,
         cuint component
      );
      virtual void fieldSolverBoundaryCondGradPeElectricField(
         fs_cache::CellCache& cache,
         cuint RKCase,
         cuint component
      );
      virtual void fieldSolverBoundaryCondDerivatives(
         dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         cuint& RKCase,
         cuint& component
      );
      virtual void fieldSolverBoundaryCondBVOLDerivatives(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         cuint& component
      );
      virtual void vlasovBoundaryCondition(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         const uint popID
      );
      
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      void generateTemplateCell(Project &project);
      void setCellFromTemplate(SpatialCell* cell,const uint popID);
      
      Real shiftedMaxwellianDistribution(const uint popID,creal& vx, creal& vy, creal& vz);
      
      vector<vmesh::GlobalID> findBlocksToInitialize(
         SpatialCell& cell,const uint popID
      );
      
      std::array<Real, 3> fieldSolverGetNormalDirection(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID
      );
      
      Real center[3]; /*!< Coordinates of the centre of the ionosphere. */
      Real radius; /*!< Radius of the ionosphere. */
      uint geometry; /*!< Geometry of the ionosphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: polar-plane cylinder with line dipole. */

      std::vector<IonosphereSpeciesParameters> speciesParams;
      Real T;
      Real rho;
      Real VX0;
      Real VY0;
      Real VZ0;
      
      uint nSpaceSamples;
      uint nVelocitySamples;
      
      spatial_cell::SpatialCell templateCell;
   };
}

#endif
