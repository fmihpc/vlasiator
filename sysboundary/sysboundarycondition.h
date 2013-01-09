/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SYSBOUNDARYCONDITION_H
#define SYSBOUNDARYCONDITION_H

#include <dccrg.hpp>
#include <vector>
#include "../definitions.h"
#include "../spatial_cell.hpp"
#include "../projects/project.h"

using namespace spatial_cell;
using namespace projects;

namespace SBC {
   /*!\brief SBC::SysBoundaryCondition is the base class for system boundary conditions.
    * 
    * SBC::SysBoundaryCondition defines a base class for applying boundary conditions.
    * Specific system boundary conditions inherit from this base class, that's why most
    * functions defined here are not meant to be called and contain a corresponding error
    * message. The functions to be called are the inherited class members.
    * 
    * The initSysBoundary function is used to initialise the internal workings needed by the
    * system boundary condition to run (e.g. importing parameters, initialising class
    * members). assignSysBoundary is used to determine whether a given cell is within the
    * domain of system boundary condition. applyInitialState is called to initialise a system
    * boundary cell's parameters and velocity space.
    * 
    * If needed, a user can write his or her own SBC::SysBoundaryConditions, which 
    * are loaded when the simulation initializes.
    */
   class SysBoundaryCondition {
      public:
         SysBoundaryCondition();
         virtual ~SysBoundaryCondition();
         
         static void addParameters();
         virtual void getParameters();
         
         virtual bool initSysBoundary(
            creal& t,
            Project &project
         );
         virtual bool assignSysBoundary(dccrg::Dccrg<SpatialCell>& mpiGrid);
         virtual bool applyInitialState(
            const dccrg::Dccrg<SpatialCell>& mpiGrid,
            Project &project
         );
//          virtual bool applySysBoundaryCondition(
//             const dccrg::Dccrg<SpatialCell>& mpiGrid,
//             creal& t
//          );
         virtual Real fieldSolverBoundaryCondMagneticField(
            const dccrg::Dccrg<SpatialCell>& mpiGrid,
            const CellID& cellID,
            creal& dt,
            cuint& component
         );
         virtual void fieldSolverBoundaryCondElectricField(
            dccrg::Dccrg<SpatialCell>& mpiGrid,
            const CellID& cellID,
            cuint RKCase,
            cuint component
         );
         virtual void fieldSolverBoundaryCondDerivatives(
            dccrg::Dccrg<SpatialCell>& mpiGrid,
            const CellID& cellID,
            cuint& RKCase,
            cuint& component
         );
         virtual void fieldSolverBoundaryCondBVOLDerivatives(
            const dccrg::Dccrg<SpatialCell>& mpiGrid,
            const CellID& cellID,
            cuint& component
         );
         static void setCellDerivativesToZero(
            const dccrg::Dccrg<SpatialCell>& mpiGrid,
            const CellID& cellID,
            cuint& component
         );
         static void setCellBVOLDerivativesToZero(
            const dccrg::Dccrg<SpatialCell>& mpiGrid,
            const CellID& cellID,
            cuint& component
         );
         virtual void vlasovBoundaryCondition(
            const dccrg::Dccrg<SpatialCell>& mpiGrid,
            const CellID& cellID
         );
         
         virtual void getFaces(bool* faces);
         virtual std::string getName() const;
         virtual uint getIndex() const;
         uint getPrecedence() const;
         bool isDynamic() const;
      
      protected:
         void determineFace(
            bool* isThisCellOnAFace,
            creal x, creal y, creal z,
            creal dx, creal dy, creal dz
         );
         void copyCellData(SpatialCell *from, SpatialCell *to);
//          void zeroCellData( SpatialCell *to);
         
         CellID getClosestNonsysboundaryCell(
            const dccrg::Dccrg<SpatialCell>& mpiGrid,
            const CellID& cellID
         );
         
         /*! Precedence value of the system boundary condition. */
         uint precedence;
         /*! Is the boundary condition dynamic in time or not. */
         bool isThisDynamic;
   };
} // namespace SBC

#endif
