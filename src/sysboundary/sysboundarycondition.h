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

#ifndef SYSBOUNDARYCONDITION_H
#define SYSBOUNDARYCONDITION_H

#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <fsgrid.hpp>

#include <vector>
#include "../common.h"
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
         virtual void getParameters()=0;
         
         virtual bool initSysBoundary(
            creal& t,
            Project &project
         )=0;
         virtual bool assignSysBoundary(dccrg::Dccrg<SpatialCell,
                                        dccrg::Cartesian_Geometry>& mpiGrid,
                                        FsGrid< fsgrids::technical, 2> & technicalGrid)=0;
         virtual bool applyInitialState(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
            Project &project
         )=0;
         virtual Real fieldSolverBoundaryCondMagneticField(
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & bGrid,
            FsGrid< fsgrids::technical, 2> & technicalGrid,
            cint i,
            cint j,
            cint k,
            creal& dt,
            cuint& component
         )=0;
         virtual void fieldSolverBoundaryCondElectricField(
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
            cint i,
            cint j,
            cint k,
            cuint component
         )=0;
         virtual void fieldSolverBoundaryCondHallElectricField(
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
            cint i,
            cint j,
            cint k,
            cuint component
         )=0;
         virtual void fieldSolverBoundaryCondGradPeElectricField(
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
            cint i,
            cint j,
            cint k,
            cuint component
         )=0;
         virtual void fieldSolverBoundaryCondDerivatives(
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
            cint i,
            cint j,
            cint k,
            cuint& RKCase,
            cuint& component
         )=0;
         virtual void fieldSolverBoundaryCondBVOLDerivatives(
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
            cint i,
            cint j,
            cint k,
            cuint& component
         )=0;
         static void setCellDerivativesToZero(
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
            cint i,
            cint j,
            cint k,
            cuint& component
         );
         static void setCellBVOLDerivativesToZero(
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
            cint i,
            cint j,
            cint k,
            cuint& component
         );
        
         /** This function computes the Vlasov (distribution function) 
          * boundary condition for the given particle species only. 
          * It is not! allowed to change block structure in cell.
          * @param mpiGrid Parallel grid.
          * @param cellID Spatial cell ID.
          * @param popID Particle species ID.*/
        virtual void vlasovBoundaryCondition(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const uint popID,
            const bool calculate_V_moments
        )=0;

         virtual void getFaces(bool* faces);
         virtual std::string getName() const=0;
         virtual uint getIndex() const=0;
         uint getPrecedence() const;
         bool isDynamic() const;
      
         bool updateSysBoundaryConditionsAfterLoadBalance(
            dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const std::vector<CellID> & local_cells_on_boundary
         );
      bool doApplyUponRestart() const;
      void setPeriodicity(
         bool isFacePeriodic[3]
      );
      protected:
         void determineFace(
            bool* isThisCellOnAFace,
            creal x, creal y, creal z,
            creal dx, creal dy, creal dz,
            const bool excludeSlicesAndPeriodicDimensions = false
         );
         void copyCellData(
            SpatialCell *from,
            SpatialCell *to,
            const bool copyMomentsOnly,
            const uint popID,
            const bool calculate_V_moments
         );
         void averageCellData(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            std::vector<CellID> cellList,
            SpatialCell *to,
            const uint popID,
            const bool calculate_V_moments,
            creal fluffiness = 0
         );
         std::array<SpatialCell*,27> & getFlowtoCells(
               const CellID& cellID
         );

         std::array<Realf*,27> getFlowtoCellsBlock(
               const std::array<SpatialCell*,27> flowtoCells,
               const vmesh::GlobalID blockGID,
               const uint popID
         );
      

      /*! Helper function to get the index of a neighboring cell in the arrays in allFlowtoCells.
       * \param i Offset in x direction (-1, 0 or 1)
       * \param j Offset in y direction (-1, 0 or 1)
       * \param k Offset in z direction (-1, 0 or 1)
       * \retval int Index in the flowto cell array (0 to 26, indexed from - to + x, y, z.
       */
      inline int nbrID(const int i, const int j, const int k){
         return (k+1)*9 + (j+1)*3 + i + 1;
      }
      
         void vlasovBoundaryCopyFromTheClosestNbr(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const bool& copyMomentsOnly,
            const uint popID,
            const bool calculate_V_moments
         );
         void vlasovBoundaryCopyFromTheClosestNbrAndLimit(
               const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
               const CellID& cellID,
               const uint popID
         );
         void vlasovBoundaryCopyFromAllClosestNbrs(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const uint popID,
            const bool calculate_V_moments
         );
         void vlasovBoundaryFluffyCopyFromAllCloseNbrs(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const uint popID,
            const bool calculate_V_moments,
            creal fluffiness
         );
         void vlasovBoundaryReflect(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            creal& nx,
            creal& ny,
            creal& nz,
            const uint popID
         );
         void vlasovBoundaryAbsorb(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            creal& nx,
            creal& ny,
            creal& nz,
            creal& quenchingFactor,
            const uint popID
         );
         std::array<int, 3> getTheClosestNonsysboundaryCell(
            FsGrid< fsgrids::technical, 2> & technicalGrid,
            cint i,
            cint j,
            cint k
         );
         std::vector< std::array<int, 3> > getAllClosestNonsysboundaryCells(
            FsGrid< fsgrids::technical, 2> & technicalGrid,
            cint i,
            cint j,
            cint k
         );
         CellID & getTheClosestNonsysboundaryCell(
            const CellID& cellID
         );
         std::vector<CellID> & getAllClosestNonsysboundaryCells(
            const CellID& cellID
         );
         std::vector<CellID> & getAllCloseNonsysboundaryCells(
            const CellID& cellID
         );
         Real fieldBoundaryCopyFromSolvingNbrMagneticField(
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & bGrid,
            FsGrid< fsgrids::technical, 2> & technicalGrid,
            cint i,
            cint j,
            cint k,
            cuint component,
            cuint mask
         );
         
         /*! Precedence value of the system boundary condition. */
         uint precedence;
         /*! Is the boundary condition dynamic in time or not. */
         bool isThisDynamic;
         /*! Array of bool telling whether the system is periodic in any direction. */
         bool isPeriodic[3];
         /*! Map of closest nonsysboundarycells. Used in getAllClosestNonsysboundaryCells. */
         std::unordered_map<CellID, std::vector<CellID>> allClosestNonsysboundaryCells;
         /*! Map of close nonsysboundarycells. Used in getAllCloseNonsysboundaryCells. */
         std::unordered_map<CellID, std::vector<CellID>> allCloseNonsysboundaryCells;
      
         /*! Array of cells into which the distribution function can flow. Used in getAllFlowtoCells. Cells into which one cannot flow are set to INVALID_CELLID. */
         std::unordered_map<CellID, std::array<SpatialCell*, 27>> allFlowtoCells;
         /*! bool telling whether to call again applyInitialState upon restarting the simulation. */
         bool applyUponRestart;
   };
   


} // namespace SBC

#endif
