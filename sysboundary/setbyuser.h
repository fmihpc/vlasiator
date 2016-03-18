/*
 This file is part of Vlasiator.
 
 Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
 */

#ifndef SETBYUSER_H
#define SETBYUSER_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {
   /*!\brief Base class for system boundary conditions with user-set settings and parameters read from file.
    * 
    * SetByUser is a base class for e.g. SysBoundaryConditon::SetMaxwellian.
    * It defines the managing functions to set boundary conditions on the faces of the
    * simulation domain.
    * 
    * This class handles the import and interpolation in time of the input parameters read
    * from file as well as the assignment of the state from the template cells.
    * 
    * The daughter classes have then to handle parameters and generate the template cells as
    * wished from the data returned.
    */
   class SetByUser: public SysBoundaryCondition {
   public:
      SetByUser();
      virtual ~SetByUser();
      
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
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
         const std::vector<fs_cache::CellCache>& cellCache,
         const uint16_t& localID,
         creal& dt,
         cuint& RKCase,
         cint& offset,
         cuint& component
      );
      virtual void fieldSolverBoundaryCondElectricField(
         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
         cuint i,
         cuint j,
         cuint k,
         cuint component
      );
      virtual void fieldSolverBoundaryCondHallElectricField(
         FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 3, 2> & EHallGrid,
         fs_cache::CellCache& cache,
         cuint RKCase,
         cuint component
      );
      virtual void fieldSolverBoundaryCondGradPeElectricField(
         FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
         cuint i,
         cuint j,
         cuint k,
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
         const CellID& cellID,
         const int& popID
      );
      
      virtual void getFaces(bool* faces);
      
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      bool loadInputData();
      std::vector<std::vector<Real> > loadFile(const char* file);
      void interpolate(const int inputDataIndex, creal t, Real* outputData);
      
      bool generateTemplateCells(creal& t);
      virtual void generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t);
      bool setCellsFromTemplate(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const int& popID);
      
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
