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
   /*!\brief Ionosphere is a class applying ionospheric boundary conditions.
    * 
    * Ionosphere is a class handling cells tagged as sysboundarytype::IONOSPHERE by this
    * system boundary condition. It applies ionospheric boundary conditions.
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
      virtual bool assignSysBoundary(dccrg::Dccrg<SpatialCell>& mpiGrid);
      virtual bool applyInitialState(
         const dccrg::Dccrg<SpatialCell>& mpiGrid,
         Project &project
      );
//       virtual bool applySysBoundaryCondition(
//          const dccrg::Dccrg<SpatialCell>& mpiGrid,
//          creal& t
//       );
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
      virtual void vlasovBoundaryCondition(
         const dccrg::Dccrg<SpatialCell>& mpiGrid,
         const CellID& cellID
      );
      
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      void generateTemplateCell(Project &project);
      void setCellFromTemplate(SpatialCell *cell);
      
      Real maxwellianDistribution(
         creal& vx, creal& vy, creal& vz
      );
      
      vector<uint> findBlocksToInitialize(
         SpatialCell& cell
      );
      
      Real center[3]; /*!< Coordinates of the centre of the ionosphere. */
      Real radius; /*!< Radius of the ionosphere. */
      
      Real T;
      Real rho;
      
      uint nSpaceSamples;
      uint nVelocitySamples;
      
      spatial_cell::SpatialCell templateCell;
   };
}

#endif
