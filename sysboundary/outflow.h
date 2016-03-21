/*
 This file is part of Vlasiator.
 
 Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
 */

#ifndef OUTFLOW_H
#define OUTFLOW_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {
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
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBDt2Grid,
         FsGrid< fsgrids::technical, 2> & technicalGrid,
         cint i,
         cint j,
         cint k,
         creal& dt,
         cuint& RKCase,
         cuint& component
      );
      virtual void fieldSolverBoundaryCondElectricField(
         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondHallElectricField(
         FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondGradPeElectricField(
         FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondDerivatives(
         FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
         FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
         cint i,
         cint j,
         cint k,
         cuint& RKCase,
         cuint& component
      );
      virtual void fieldSolverBoundaryCondBVOLDerivatives(
         FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
         cint i,
         cint j,
         cint k,
         cuint& component
      );
      virtual void vlasovBoundaryCondition(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         const int& popID
      );
      
      virtual void getFaces(bool* faces);
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      /*! Array of bool telling which faces are going to be processed by the system boundary condition.*/
      bool facesToProcess[6];
      /*! List of faces on which outflow boundary conditions are to be applied ([xyz][+-]). */
      std::vector<std::string> faceList;
   }; // class Outflow
} // namespace SBC

#endif
