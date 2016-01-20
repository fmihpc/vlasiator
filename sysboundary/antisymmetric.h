/*
 This file is part of Vlasiator.
 
 Copyright 2015 Finnish Meteorological Institute
  
 */

#ifndef ANTISYMMETRIC_H
#define ANTISYMMETRIC_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {
   /*!\brief Antisymmetric is a class applying antisymmetric boundary conditions.
    * Normal field components are antisymmetric.
    * Tangential field components are symmetric.
    * 
    * Distribution function is mirrored in the same way.
    * 
    * Antisymmetric is a class handling cells tagged as sysboundarytype::ANTISYMMETRIC by this system
    * boundary condition. It applies antisymmetric boundary conditions.
    */
   class Antisymmetric: public SysBoundaryCondition {
   public:
      Antisymmetric();
      virtual ~Antisymmetric();
      
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
      Real fieldBoundaryCopyFromExistingFaceNbrMagneticField(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         cuint& component
      );
   }; // class Antisymmetric
} // namespace SBC

#endif
