/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute


*/

#ifndef SYSBOUNDARYCONDITION_H
#define SYSBOUNDARYCONDITION_H

#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include <vector>
#include "../definitions.h"
#include "../spatial_cell.hpp"
#include "../projects/project.h"
#include "../fieldsolver/fs_cache.h"

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
         static void setCellDerivativesToZero(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            cuint& component
         );
         static void setCellBVOLDerivativesToZero(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
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
            const int& popID
        );
        virtual void fieldSolverBoundaryCondGradPeElectricField(
           fs_cache::CellCache& cache,
           cuint RKCase,
           cuint component
           );

         virtual void getFaces(bool* faces);
         virtual std::string getName() const;
         virtual uint getIndex() const;
         uint getPrecedence() const;
         bool isDynamic() const;

         bool updateSysBoundaryConditionsAfterLoadBalance(
            dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const std::vector<CellID> & local_cells_on_boundary
         );
         bool doApplyUponRestart() const;
      protected:
         void determineFace(
            bool* isThisCellOnAFace,
            creal x, creal y, creal z,
            creal dx, creal dy, creal dz,
            const bool excludeSlices = false
         );
         void copyCellData(
            SpatialCell *from,
            SpatialCell *to,
            bool allowBlockAdjustment,
            const int& popID
         );
         void averageCellData(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            std::vector<CellID> cellList,
            SpatialCell *to,
            const int& popID
         );
         CellID & getTheClosestNonsysboundaryCell(
            const CellID& cellID
         );
         std::vector<CellID> & getAllClosestNonsysboundaryCells(
            const CellID& cellID
         );
         void vlasovBoundaryCopyFromTheClosestNbr(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const int& popID
         );
         void vlasovBoundaryCopyFromAllClosestNbrs(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            const int& popID
         );
         void vlasovBoundaryReflect(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            creal& nx,
            creal& ny,
            creal& nz,
            const int& popID
         );
         void vlasovBoundaryAbsorb(
            const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
            const CellID& cellID,
            creal& nx,
            creal& ny,
            creal& nz,
            const int& popID
         );
         
         /*! Precedence value of the system boundary condition. */
         uint precedence;
         /*! Is the boundary condition dynamic in time or not. */
         bool isThisDynamic;
         /*! Map of closest nonsysboundarycells. Used in getAllClosestNonsysboundaryCells. */
         std::unordered_map<CellID, std::vector<CellID>> allClosestNonsysboundaryCells;
         /*! bool telling whether to call again applyInitialState upon restarting the simulation. */
         bool applyUponRestart;
   };
} // namespace SBC

#endif
