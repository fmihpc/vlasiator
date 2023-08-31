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

#ifndef DONOTCOMPUTE_H
#define DONOTCOMPUTE_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

using namespace projects;

namespace SBC {
   /*!\brief DoNotCompute is a class handling cells not to be computed.
    * 
    * DoNotCompute is a class handling cells tagged as sysboundarytype::DO_NOT_COMPUTE by a system boundary condition (e.g. SysBoundaryCondition::Ionosphere).
    */
   class DoNotCompute: public SysBoundaryCondition {
   public:
      DoNotCompute();
      virtual ~DoNotCompute();
      
      static void addParameters();
      virtual void getParameters();
      
      virtual bool initSysBoundary(
         creal& t,
         Project &project
      );
      virtual bool assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                     FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid);
      virtual bool applyInitialState(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
         FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> & perBGrid,
         Project &project
      );
      virtual std::string getName() const;
      virtual uint getIndex() const;

      bool initFieldBoundary() {return false;}

      // Explicit warning functions to inform the user if a doNotCompute cell gets computed
      ARCH_HOSTDEV virtual Real fieldSolverBoundaryCondMagneticField(
         const arch::buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & perBGrid,
         const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
         cint i,
         cint j,
         cint k,
         creal& dt,
         cuint& component
      ) { 
         #ifndef __CUDA_ARCH__
         std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondMagneticField called!" << std::endl;
         #endif 
         return 0.;
      }
      ARCH_HOSTDEV virtual void fieldSolverBoundaryCondMagneticFieldProjection(
         const arch::buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & perBGrid,
         const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
         cint i,
         cint j,
         cint k
      ) {
         #ifndef __CUDA_ARCH__
         std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondMagneticFieldProjection called!" << std::endl;
         #endif
      }
      ARCH_HOSTDEV virtual void fieldSolverBoundaryCondElectricField(
         const arch::buf<FsGrid<Real, fsgrids::efield::N_EFIELD, FS_STENCIL_WIDTH>> & EGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) {
         #ifndef __CUDA_ARCH__
         std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondElectricField called!" << std::endl;
         #endif
      }
      ARCH_HOSTDEV virtual void fieldSolverBoundaryCondHallElectricField(
         const arch::buf<FsGrid<Real, fsgrids::ehall::N_EHALL, FS_STENCIL_WIDTH>> & EHallGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) {
         #ifndef __CUDA_ARCH__
         std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondHallElectricField called!" << std::endl;
         #endif
      }
      ARCH_HOSTDEV virtual void fieldSolverBoundaryCondGradPeElectricField(
         const arch::buf<FsGrid<Real, fsgrids::egradpe::N_EGRADPE, FS_STENCIL_WIDTH>> & EGradPeGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) {
         #ifndef __CUDA_ARCH__ 
         std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondGradPeElectricField called!" << std::endl;
         #endif
      }
      ARCH_HOSTDEV virtual void fieldSolverBoundaryCondDerivatives(
         const arch::buf<FsGrid<Real, fsgrids::dperb::N_DPERB, FS_STENCIL_WIDTH>> & dPerBGrid,
         const arch::buf<FsGrid<Real, fsgrids::dmoments::N_DMOMENTS, FS_STENCIL_WIDTH>> & dMomentsGrid,
         cint i,
         cint j,
         cint k,
         cuint& RKCase,
         cuint& component
      ) {
         #ifndef __CUDA_ARCH__
         std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondDerivatives called!" << std::endl;
         #endif
      }
      ARCH_HOSTDEV virtual void fieldSolverBoundaryCondBVOLDerivatives(
         const arch::buf<FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH>> & volGrid,
         cint i,
         cint j,
         cint k,
         cuint& component
      ) {
         #ifndef __CUDA_ARCH__
         std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondBVOLDerivatives called!" << std::endl;
         #endif
      }
      virtual void vlasovBoundaryCondition(
          const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
          const CellID& cellID,
          const uint popID,
          const bool calculate_V_moments
      ) { std::cerr << "ERROR: DoNotCompute::vlasovBoundaryCondition called!" << std::endl;}
   };
}

#endif
