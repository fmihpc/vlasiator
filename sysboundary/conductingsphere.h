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

#ifndef CONDUCTINGSPHERE_H
#define CONDUCTINGSPHERE_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"
#include "conductingsphereFieldBoundary.h"

using namespace projects;
using namespace std;

namespace SBC {

   struct ConductingsphereSpeciesParameters {
      Real rho;
      Real V0[3];
      Real T;
      Real fluffiness;
      uint nSpaceSamples;
      uint nVelocitySamples;
   };

   /*!\brief Conductingsphere is a class applying an ionosphere-ish boundary conditions.
    * 
    * Conductingsphere is a class handling cells tagged as sysboundarytype::CONDUCTINGSPHERE by this system boundary condition. It applies perfectly conducting boundary conditions.
    * 
    * These consist in:
    * - Do nothing for the distribution (keep the initial state constant in time);
    * - Keep only the normal perturbed B component and null out the other perturbed components (perfect conductor behavior);
    * - Null out the electric fields.
    *
    * For 3D magnetospheric simulations, you might be interesting in trying the ionosphere boundary instead!
    */
   class Conductingsphere: public SysBoundaryCondition {
   public:
      Conductingsphere();
      virtual ~Conductingsphere();
      
      static void addParameters();
      virtual void getParameters();
      
      ConductingSphereFieldBoundary* getFieldBoundary() {return fieldBoundary;} 
      bool initFieldBoundary(); 
      
      bool initSysBoundary(
         creal& t,
         Project &project
      );
      virtual bool assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                     FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid);
      virtual bool applyInitialState(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
         FsGrid< Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> & perBGrid,
         Project &project
      );
      ARCH_HOSTDEV Real fieldSolverBoundaryCondMagneticField(
         const arch::buf<FsGrid< Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
         const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
         cint i,
         cint j,
         cint k,
         creal& dt,
         cuint& component
      ) {
         return fieldBoundary->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondMagneticFieldProjection(
         const arch::buf<FsGrid< Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
         const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
         cint i,
         cint j,
         cint k
      ) {
         fieldBoundary->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondElectricField(
         const arch::buf<FsGrid< Real, fsgrids::efield::N_EFIELD, FS_STENCIL_WIDTH>> & EGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) {
         fieldBoundary->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondHallElectricField(
         const arch::buf<FsGrid< Real, fsgrids::ehall::N_EHALL, FS_STENCIL_WIDTH>> & EHallGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) {
         fieldBoundary->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondGradPeElectricField(
         const arch::buf<FsGrid< Real, fsgrids::egradpe::N_EGRADPE, FS_STENCIL_WIDTH>> & EGradPeGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) {
         fieldBoundary->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondDerivatives(
         const arch::buf<FsGrid< Real, fsgrids::dperb::N_DPERB, FS_STENCIL_WIDTH>> & dPerBGrid,
         const arch::buf<FsGrid< Real, fsgrids::dmoments::N_DMOMENTS, FS_STENCIL_WIDTH>> & dMomentsGrid,
         cint i,
         cint j,
         cint k,
         cuint& RKCase,
         cuint& component
      ) {
         fieldBoundary->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondBVOLDerivatives(
         const arch::buf<FsGrid< Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH>> & volGrid,
         cint i,
         cint j,
         cint k,
         cuint& component
      ) {
         fieldBoundary->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
      }
      virtual void vlasovBoundaryCondition(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         const uint popID,
         const bool calculate_V_moments
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
         FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
         cint i,
         cint j,
         cint k
      );
      
      Real center[3]; /*!< Coordinates of the centre of the conducting sphere. */
      Real radius; /*!< Radius of the conducting sphere. */
      uint geometry; /*!< Geometry of the conducting sphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: polar-plane cylinder with line dipole. */

      std::vector<ConductingsphereSpeciesParameters> speciesParams;
      Real T;
      Real rho;
      Real VX0;
      Real VY0;
      Real VZ0;
      
      uint nSpaceSamples;
      uint nVelocitySamples;
      
      spatial_cell::SpatialCell templateCell;

      ConductingSphereFieldBoundary* fieldBoundary;

   };
}

#endif
