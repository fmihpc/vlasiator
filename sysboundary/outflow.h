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

#ifndef OUTFLOW_H
#define OUTFLOW_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {

	struct OutflowSpeciesParameters {
      /*! Array of bool telling which faces are going to be skipped by the Vlasov system boundary condition.*/
			std::array<bool, 6> facesToSkipVlasov;
      /*! List of schemes to use for the Vlasov outflow boundary conditions on each face ([xyz][+-]). */
      std::array<uint, 6> faceVlasovScheme;

      /*! Factor by which to quench the inflowing parts of the velocity distribution function.*/
      Real quenchFactor;
	};

   /*!\brief Outflow is a class applying copy/outflow boundary conditions.
    * 
    * Outflow is a class handling cells tagged as sysboundarytype::OUTFLOW by this system boundary condition. It applies copy/outflow boundary conditions.
    * 
    * These consist in:
    * - Copy the distribution and moments from the nearest NOT_SYSBOUNDARY cell;
    * - Copy the perturbed B components from the nearest NOT_SYSBOUNDARY cell. EXCEPTION: the face components adjacent to the simulation domain at the +x/+y/+z faces are propagated still.
    */
   class Outflow: public SysBoundaryCondition {
   public:
      Outflow();
      virtual ~Outflow();
      
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
      
      virtual void getFaces(bool* faces);
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      /*! Array of bool telling which faces are going to be processed by the system boundary condition.*/
      bool facesToProcess[6];
      /*! Array of bool telling which faces are going to be processed by the fields system boundary condition.*/
      bool facesToSkipFields[6];
      /*! List of faces on which outflow boundary conditions are to be applied ([xyz][+-]). */
      std::vector<std::string> faceList;
      /*! List of faces on which no fields outflow boundary conditions are to be applied ([xyz][+-]). */
      std::vector<std::string> faceNoFieldsList;
      Real fieldBoundaryCopyFromExistingFaceNbrMagneticField(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         cuint& component
      );
			std::vector<OutflowSpeciesParameters> speciesParams;
      
      enum vlasovscheme {
         NONE,
         COPY,
         LIMIT,
         N_SCHEMES
      };
      
   }; // class Outflow
} // namespace SBC

#endif
