/*
 This file is part of Vlasiator.
 
 Copyright 2015 Finnish Meteorological Institute
  
 */

#ifndef PROJECT_BOUNDARY_H
#define PROJECT_BOUNDARY_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {
   /*!\brief Base class for system boundary conditions with user-set settings and parameters read from file.
    * 
    * ProjectBoundary uses the simulated project to set the boundary conditions.
    * 
    * This class handles the import and interpolation in time of the input parameters read
    * from file as well as the assignment of the state from the template cells.
    * 
    * The daughter classes have then to handle parameters and generate the template cells as
    * wished from the data returned. 
    */
   class ProjectBoundary: public SysBoundaryCondition {
   public:
      ProjectBoundary();
      virtual ~ProjectBoundary();
      
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
//       virtual bool applySysBoundaryCondition(
//          const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
//          creal& t
//       );
      virtual Real fieldSolverBoundaryCondMagneticField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
         FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
         cint i,
         cint j,
         cint k,
         creal& dt,
         cuint& RKCase,
         cuint& component
      );
      virtual void fieldSolverBoundaryCondElectricField(
         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondHallElectricField(
         FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 3, 2> & EHallGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondGradPeElectricField(
         FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondDerivatives(
         FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
         FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
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
         const CellID& cellID,const int& popID
      );
      
      virtual void getFaces(bool* faces);
      
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:

      bool generateTemplateCell();

      /*! Array of bool telling which faces are going to be processed by the system boundary condition.*/
      bool facesToProcess[6];

      Project* project;
      spatial_cell::SpatialCell templateCell;

      /*! List of faces on which user-set boundary conditions are to be applied ([xyz][+-]). */
      std::vector<std::string> faceList;
   };
}

#endif
