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
         dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
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
      
      bool generateTemplateCells(creal& t);
      virtual void generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t);
      bool setCellsFromTemplate(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
      
      /*! Array of bool telling which faces are going to be processed by the system boundary condition.*/
      bool facesToProcess[6];
      /*! Vector containing a vector for each face which has the current boundary condition. Each of these vectors has one line per input data line (time point). The length of the lines is nParams.*/
      std::vector<std::vector<Real> > inputData[6];
      /*! Array of template spatial cells replicated over the corresponding simulation volume face. Only the template for an active face is actually being touched at all by the code. */
      spatial_cell::SpatialCell templateCells[6];
      /*! List of faces on which user-set boundary conditions are to be applied ([xyz][+-]). */
      std::vector<std::string> faceList;
      /*! Input files for the user-set boundary conditions. */
      std::string files[6];
      /*! Number of parameters per input file line. */
      uint nParams;
   };
}

#endif
