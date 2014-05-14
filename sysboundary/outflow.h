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
   /*!\brief Outflow is a class applying outflow boundary conditions.
    * 
    * Outflow is a class handling cells tagged as sysboundarytype::OUTFLOW by this system
    * boundary condition. It applies outflow boundary conditions.
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
//       virtual bool applySysBoundaryCondition(
//          const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
//          creal& t
//       );
      virtual Real fieldSolverBoundaryCondMagneticField(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         creal& dt,
         cuint& component
      );
      virtual void fieldSolverBoundaryCondElectricField(
         dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
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
         const CellID& cellID
      );
      
      virtual void getFaces(bool* faces);
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      /*! Array of bool telling which faces are going to be processed by the system boundary condition.*/
      bool facesToProcess[6];
      /*! List of faces on which outflow boundary conditions are to be applied ([xyz][+-]). */
      std::vector<std::string> faceList;
      Real fieldBoundaryCopyFromExistingFaceNbrMagneticField(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         cuint& component
      );
      void vlasovBoundaryCopyFromExistingFaceNbr(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID
      );
   }; // class Outflow
} // namespace SBC

#endif
