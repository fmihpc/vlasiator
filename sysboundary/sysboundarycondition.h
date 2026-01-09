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
#include "../spatial_cells/spatial_cell_wrapper.hpp"
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
         
         virtual void initSysBoundary(
            creal& t,
            Project &project
         )=0;
         virtual void assignSysBoundary(dccrg::Dccrg<SpatialCell,
                                        dccrg::Cartesian_Geometry>& mpiGrid,
                                        FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)=0;
         virtual void applyInitialState(
            dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
            Project &project
         )=0;
         virtual void updateState(
            dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
            FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
            FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
            creal t
         )=0;
         virtual Real fieldSolverBoundaryCondMagneticField(
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & bgbGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
            cint i,
            cint j,
            cint k,
            creal dt,
            cuint component
         )=0;
         virtual void fieldSolverBoundaryCondElectricField(
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            cint i,
            cint j,
            cint k,
            cuint component
         )=0;
         virtual void fieldSolverBoundaryCondHallElectricField(
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            cint i,
            cint j,
            cint k,
            cuint component
         )=0;
         virtual void fieldSolverBoundaryCondGradPeElectricField(
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            cint i,
            cint j,
            cint k,
            cuint component
         )=0;
         virtual void fieldSolverBoundaryCondDerivatives(
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            cint i,
            cint j,
            cint k,
            cuint RKCase,
            cuint component
         )=0;
         virtual void fieldSolverBoundaryCondBVOLDerivatives(
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            cint i,
            cint j,
            cint k,
            cuint component
         )=0;
         static void setCellDerivativesToZero(
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            cint i,
            cint j,
            cint k,
            cuint component
         );
         static void setCellBVOLDerivativesToZero(
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            cint i,
            cint j,
            cint k,
            cuint component
         );
        
         virtual void mapCellPotentialAndGetEXBDrift(
            std::array<Real, CellParams::N_SPATIAL_CELL_PARAMS>& cellParams
         );

         /** This function computes the Vlasov (distribution function) 
          * boundary condition for the given particle species only. 
          * It is not! allowed to change block structure in cell.
          * @param mpiGrid Parallel grid.
          * @param cellID Spatial cell ID.
          * @param popID Particle species ID.*/
        virtual void vlasovBoundaryCondition(
            dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const uint popID,
            const bool calculate_V_moments
        )=0;

        virtual void setupL2OutflowAtRestart(
            dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid
        ) {
            std::cerr << "ERROR: base class SysBoundaryCondition::setupL2OutflowAtRestart called!" << std::endl;
        }



         /*! Function used to know which faces the boundary condition is applied to.
          * @param faces Pointer to array of 6 bool in which the values are returned whether the corresponding face is of that
          * type. Order: 0 x+; 1 x-; 2 y+; 3 y-; 4 z+; 5 z-
          */
         virtual void getFaces(bool *faces) = 0;
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
         void determineFace(
            std::array<bool, 6> &isThisCellOnAFace,
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            CellID id,
            const bool excludeSlicesAndPeriodicDimensions = false
         );
         void copyCellData(
            SpatialCell *from,
            SpatialCell *to,
            const bool copyMomentsOnly,
            const uint popID,
            const bool copy_V_moments
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
            dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const bool& copyMomentsOnly,
            const uint popID,
            const bool calculate_V_moments
         );
         void vlasovBoundaryCopyFromTheClosestL1OutflowNbr(
            dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const bool& copyMomentsOnly,
            const uint popID,
            const bool calculate_V_moments
         );
         void vlasovBoundaryCopyFromAllClosestNbrs(
            dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const uint popID,
            const bool calculate_V_moments
         );
         void vlasovBoundaryFluffyCopyFromAllCloseNbrs(
            dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const uint popID,
            const bool calculate_V_moments,
            creal fluffiness
         );
         std::array<int, 3> getTheClosestNonsysboundaryCell(
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
            cint i,
            cint j,
            cint k
         );
         std::vector< std::array<int, 3> > getAllClosestNonsysboundaryCells(
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
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
         CellID & getTheClosestL1OutflowCell(
            const CellID& cellID
         );
         std::vector<CellID> & getAllClosestL1OutflowCells(
            const CellID& cellID
         );
         Real fieldBoundaryCopyFromSolvingNbrMagneticField(
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
            cint i,
            cint j,
            cint k,
            cuint component,
            cuint mask
         );
         
         /*! Precedence value of the system boundary condition. */
         uint precedence;
         /*! Is the boundary condition dynamic in time or not. */
         bool dynamic;
         /*! Array of bool telling whether the system is periodic in any direction. */
         bool periodic[3];
         /*! Map of closest nonsysboundarycells. Used in getAllClosestNonsysboundaryCells. */
         std::unordered_map<CellID, std::vector<CellID>> allClosestNonsysboundaryCells;
         /*! Map of close nonsysboundarycells. Used in getAllCloseNonsysboundaryCells. */
         std::unordered_map<CellID, std::vector<CellID>> allCloseNonsysboundaryCells;
         /*! Map of closest Outflow L1 cells. Used in getAllClosestL1OutflowCells. */
         std::unordered_map<CellID, std::vector<CellID>> allClosestL1OutflowCells;
         /*! Map of close Outflow L1 cells. Used in getAllCloseL1OutflowCells. */
         std::unordered_map<CellID, std::vector<CellID>> allCloseL1OutflowCells;
      
         /*! bool telling whether to call again applyInitialState upon restarting the simulation. */
         bool applyUponRestart;
   };

   class OuterBoundaryCondition: public SysBoundaryCondition {
      public:
         virtual void assignSysBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid);
      protected:
         /*! Array of bool telling which faces are going to be processed by the system boundary condition.*/
         bool facesToProcess[6];
   };
   
   // Moved outside the class since it's a helper function that doesn't require member access
   void averageCellData (
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      std::vector<CellID> cellList,
      SpatialCell *to,
      const uint popID,
      creal fluffiness = 0
   );

   /*!\brief SBC::findMaxwellianBlocksToInitialize returns a list of blocks to construct the VDF with.
    * 
    *  Here the while loop iterates  from the centre of the maxwellian in blocksize (4*dvx) increments,
    *  and looks at the centre of the first velocity cell in the block (+0.5dvx), checking if the
    *  phase-space density there is large enough to be included due to sparsity threshold.
    *  That results in a "blocks radius"  vRadiusSquared from the centre of the Maxwellian distribution.
    *  Then we iterate through the actual blocks and calculate their radius R2 based on their velocity coordinates
    *  and the plasma bulk velocity. Blocks that fullfil R2<vRadiusSquared are included to blocksToInitialize.
    */
   vmesh::LocalID findMaxwellianBlocksToInitialize(
      const uint popID,
      spatial_cell::SpatialCell& cell,
      creal& rho,
      creal& T,
      creal& VX0,
      creal& VY0,
      creal& VZ0);

} // namespace SBC

#endif
