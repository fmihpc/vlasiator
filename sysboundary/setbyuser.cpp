/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010-2013,2015 Finnish Meteorological Institute
 * 
 */

/*!\file setbyuser.cpp
 * \brief Implementation of the class SysBoundaryCondition::SetByUser. 
 * This serves as the base class for further classes like SysBoundaryCondition::SetMaxwellian.
 */

#include <cstdlib>
#include <iostream>

#include "setbyuser.h"
#include "../vlasovmover.h"
#include "../fieldsolver/fs_common.h"
#include "../object_wrapper.h"

#ifndef NDEBUG
   #define DEBUG_SETBYUSER
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_SETBYUSER
#endif

using namespace std;

namespace SBC {
   SetByUser::SetByUser(): SysBoundaryCondition() { }
   SetByUser::~SetByUser() { }
   
   void SetByUser::addParameters() {
      cerr << "Base class SetByUser::addParameters() called instead of derived class function!" << endl;
   }
   
   void SetByUser::getParameters() {
      cerr << "Base class SetByUser::getParameters() called instead of derived class function!" << endl;
   }
   
   bool SetByUser::initSysBoundary(
      creal& t,
      Project &project
   ) {
      /* The array of bool describes which of the x+, x-, y+, y-, z+, z- faces are to have user-set system boundary conditions.
       * A true indicates the corresponding face will have user-set system boundary conditions.
       * The 6 elements correspond to x+, x-, y+, y-, z+, z- respectively.
       */
      bool success = true;
      for(uint i=0; i<6; i++) facesToProcess[i] = false;
      
      this->getParameters();
      
      vector<string>::const_iterator it;
      for (it = faceList.begin(); it != faceList.end(); ++it) {
         if(*it == "x+") facesToProcess[0] = true;
         if(*it == "x-") facesToProcess[1] = true;
         if(*it == "y+") facesToProcess[2] = true;
         if(*it == "y-") facesToProcess[3] = true;
         if(*it == "z+") facesToProcess[4] = true;
         if(*it == "z-") facesToProcess[5] = true;
      }
      
      success = loadInputData();
      success = success & generateTemplateCells(t);
      
      return success;
   }
   
   bool SetByUser::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      vector<CellID> cells = mpiGrid.get_cells();
      for(uint i = 0; i < cells.size(); i++) {
         if(mpiGrid[cells[i]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         creal* const cellParams = &(mpiGrid[cells[i]]->parameters[0]);
         creal dx = cellParams[CellParams::DX];
         creal dy = cellParams[CellParams::DY];
         creal dz = cellParams[CellParams::DZ];
         creal x = cellParams[CellParams::XCRD] + 0.5*dx;
         creal y = cellParams[CellParams::YCRD] + 0.5*dy;
         creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
         
         bool isThisCellOnAFace[6];
         determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);
         // Comparison of the array defining which faces to use and the array telling on which faces this cell is
         bool doAssign = false;
         for(uint j=0; j<6; j++) doAssign = doAssign || (facesToProcess[j] && isThisCellOnAFace[j]);
         if(doAssign) {
            mpiGrid[cells[i]]->sysBoundaryFlag = this->getIndex();
         }
      }
      return true;
   }
   
   bool SetByUser::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      Project &project
   ) {
      bool success = true;
      for (unsigned int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         if (setCellsFromTemplate(mpiGrid, popID) == false) success = false;
      }
      
      return success;
   }
   
   Real SetByUser::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<fs_cache::CellCache>& cellCache,
      const uint16_t& localID,
      creal& dt,
      cuint& RKCase,
      cint& offset,
      cuint& component
   ) {
      Real result = 0.0;
      creal* cp0 = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
      creal dx = cp0[CellParams::DX];
      creal dy = cp0[CellParams::DY];
      creal dz = cp0[CellParams::DZ];
      creal x = cp0[CellParams::XCRD] + 0.5*dx;
      creal y = cp0[CellParams::YCRD] + 0.5*dy;
      creal z = cp0[CellParams::ZCRD] + 0.5*dz;
      
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);

      for (uint i=0; i<6; i++) {
         if (isThisCellOnAFace[i]) {
            result = templateCells[i].parameters[CellParams::PERBX + offset + component];
            break; // This effectively sets the precedence of faces through the order of faces.
         }
      }
      return result;
   }

   void SetByUser::fieldSolverBoundaryCondElectricField(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint RKCase,
      cuint component
   ) {
      if((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
         mpiGrid[cellID]->parameters[CellParams::EX+component] = 0.0;
      } else {// RKCase == RK_ORDER2_STEP1
         mpiGrid[cellID]->parameters[CellParams::EX_DT2+component] = 0.0;
      }
   }

   void SetByUser::fieldSolverBoundaryCondHallElectricField(
      fs_cache::CellCache& cache,
      cuint RKCase,
      cuint component
   ) {

      Real* cp = cache.cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
      
      switch (component) {
         case 0:
            cp[CellParams::EXHALL_000_100] = 0.0;
            cp[CellParams::EXHALL_010_110] = 0.0;
            cp[CellParams::EXHALL_001_101] = 0.0;
            cp[CellParams::EXHALL_011_111] = 0.0;
            break;
         case 1:
            cp[CellParams::EYHALL_000_010] = 0.0;
            cp[CellParams::EYHALL_100_110] = 0.0;
            cp[CellParams::EYHALL_001_011] = 0.0;
            cp[CellParams::EYHALL_101_111] = 0.0;
            break;
         case 2:
            cp[CellParams::EZHALL_000_001] = 0.0;
            cp[CellParams::EZHALL_100_101] = 0.0;
            cp[CellParams::EZHALL_010_011] = 0.0;
            cp[CellParams::EZHALL_110_111] = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   void SetByUser::fieldSolverBoundaryCondGradPeElectricField(
      fs_cache::CellCache& cache,
      cuint RKCase,
      cuint component
   ) {
      
      Real* cp = cache.cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
      
      switch (component) {
         case 0:
            //             cp[CellParams::EXGRADPE_000_100] = 0.0;
            //             cp[CellParams::EXGRADPE_010_110] = 0.0;
            //             cp[CellParams::EXGRADPE_001_101] = 0.0;
            //             cp[CellParams::EXGRADPE_011_111] = 0.0;
            cp[CellParams::EXGRADPE] = 0.0;
            break;
         case 1:
            //             cp[CellParams::EYGRADPE_000_010] = 0.0;
            //             cp[CellParams::EYGRADPE_100_110] = 0.0;
            //             cp[CellParams::EYGRADPE_001_011] = 0.0;
            //             cp[CellParams::EYGRADPE_101_111] = 0.0;
            cp[CellParams::EYGRADPE] = 0.0;
            break;
         case 2:
            //             cp[CellParams::EZGRADPE_000_001] = 0.0;
            //             cp[CellParams::EZGRADPE_100_101] = 0.0;
            //             cp[CellParams::EZGRADPE_010_011] = 0.0;
            //             cp[CellParams::EZGRADPE_110_111] = 0.0;
            cp[CellParams::EZGRADPE] = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   void SetByUser::fieldSolverBoundaryCondDerivatives(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& RKCase,
      cuint& component
   ) {
      this->setCellDerivativesToZero(mpiGrid, cellID, component);
   }

   void SetByUser::fieldSolverBoundaryCondBVOLDerivatives(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      this->setCellBVOLDerivativesToZero(mpiGrid, cellID, component);
   }

   void SetByUser::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const int& popID
   ) {
      // No need to do anything in this function, as the propagators do not touch the distribution function   
   }
   
   bool SetByUser::setCellsFromTemplate(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const int& popID) {
      vector<CellID> cells = mpiGrid.get_cells();
      #pragma omp parallel for
      for (size_t i=0; i<cells.size(); i++) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if(cell->sysBoundaryFlag != this->getIndex()) continue;
         
         creal dx = cell->parameters[CellParams::DX];
         creal dy = cell->parameters[CellParams::DY];
         creal dz = cell->parameters[CellParams::DZ];
         creal x = cell->parameters[CellParams::XCRD] + 0.5*dx;
         creal y = cell->parameters[CellParams::YCRD] + 0.5*dy;
         creal z = cell->parameters[CellParams::ZCRD] + 0.5*dz;
         
         bool isThisCellOnAFace[6];
         determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);
         
         for(uint i=0; i<6; i++) {
            if(facesToProcess[i] && isThisCellOnAFace[i]) {
               if (popID == 0) {
                  cell->parameters[CellParams::PERBX] = templateCells[i].parameters[CellParams::PERBX];
                  cell->parameters[CellParams::PERBY] = templateCells[i].parameters[CellParams::PERBY];
                  cell->parameters[CellParams::PERBZ] = templateCells[i].parameters[CellParams::PERBZ];
               
                  cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
                  cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
               }

               copyCellData(&templateCells[i], cell,true,false,popID);
               break; // This effectively sets the precedence of faces through the order of faces.
            }
         }
      }
      return true;
   }
   
   void SetByUser::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) faces[i] = facesToProcess[i];
   }
   
   bool SetByUser::loadInputData() {
      for(uint i=0; i<6; i++) {
         if(facesToProcess[i]) {
            inputData[i] = loadFile(&(files[i][0]));
         } else {
            vector<Real> tmp1;
            vector<vector<Real> > tmp2;
            for(uint j=0; j<nParams; j++) {
               tmp1.push_back(-1.0);
            }
            tmp2.push_back(tmp1);
            inputData[i] = tmp2;
         }
      }
      return true;
   }
   
   /*! Load user-set boundary data from given file.
    * The first entry of each line is assumed to be the time.
    * The number of entries per line is given by nParams which is defined as a parameter
    * from the configuration file/command line.
    * 
    * Function adapted from GUMICS-5.
    * 
    * \param fn Name of the file to be opened.
    * \retval dataset Vector of Real vectors. Each line of length nParams is put into a vector. Each of these is then put into the vector returned here.
    */
   vector<vector<Real> > SetByUser::loadFile(const char *fn) {
      vector<vector<Real> > dataset;
 
   
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      
      // Count lines with data
      FILE *fp;
      fp = fopen(fn,"r");
      if (fp == NULL) {
         cerr << "Couldn't open parameter file " << fn << endl;
         exit(1);
      }
      uint nlines = 0;
      int ret = nParams;

      // Make sure the type id of Real is correct
      phiprof_assert( typeid( Real ) == typeid(float) || typeid( Real ) == typeid(double) );

      while (!feof(fp) && ret == (int)nParams) {
         Real readParam;
         ret = 0;
         if ( typeid( readParam ) == typeid(double) ) {
            for(uint i=0; i<nParams; i++) ret += fscanf(fp, "%lf", &readParam);
         } else if( typeid( readParam ) == typeid(float) ) {
            for(uint i=0; i<nParams; i++) ret += fscanf(fp, "%f", &readParam);
         } else {
            phiprof_assert( typeid( readParam ) == typeid(float) || typeid( readParam ) == typeid(double) );
            
         }
         nlines++;
      }
      nlines--;
      
      fclose(fp);
      
      if (nlines < 1) {
         cerr << "Parameter file must have at least one value (t, n, T...)" << endl;
         exit(1);
      }
      
      if (myRank == 0) cout << "Parameter data file (" << fn << ") has " << nlines << " values"<< endl;
      
      fp = fopen(fn,"r");
      for (uint line=0; line<nlines; line++) {
         vector<Real> tempData;
         for (uint i=0; i<nParams; i++) {
            Real readParam;
            int ret;
            if ( typeid( readParam ) == typeid(double) ) {
               ret = fscanf(fp,"%lf",&readParam);
            } else if( typeid( readParam ) == typeid(float) ) {
               ret = fscanf(fp,"%f",&readParam);
            } else {
               phiprof_assert( typeid( readParam ) == typeid(float) || typeid( readParam ) == typeid(double) ); 
            }
            if (ret != 1) {
               cerr << "Couldn't read a number from parameter file " << *fn << " for line value " << line << endl;
            }
            tempData.push_back(readParam);
         }
         dataset.push_back(tempData);
      }
      
      // check that sw data is in ascending temporal order
      for (uint line = 1; line < nlines; line++) {
         if (dataset[line][0] < dataset[line - 1][0]) {
            cerr << "Parameter data must be in ascending temporal order" << endl;
            exit(1);
         }
      }
      
      fclose(fp);
      
      return dataset;
   }
   
   /*! Loops through the array of template cells and generates the ones needed. The function
    * generateTemplateCell is defined in the inheriting class such as to have the specific
    * condition needed.
    * \param t Simulation time.
    * \sa generateTemplateCell
    */
   bool SetByUser::generateTemplateCells(creal& t) {
      #pragma omp parallel for
      for(uint i=0; i<6; i++) {
         int index;
         if(facesToProcess[i]) {
            generateTemplateCell(templateCells[i], i, t);
         }
      }
      return true;
   }
   
   /*!Interpolate the input data to the given time.
    * The first entry of each line is assumed to be the time.
    * \param inputDataIndex Index used to get the correct face's input data.
    * \param t Current simulation time.
    * \param outputData Pointer to the location where to write the result. Make sure from the calling side that nParams Real values can be written there!
    */
   void SetByUser::interpolate(
      const int inputDataIndex,
      creal t,
      Real* outputData
   ) {
      // Find first data[0] value which is >= t
      int i1=0,i2=0;
      bool found = false;
      Real s;      // 0 <= s < 1
      
      // use first value of sw data if interpolating for time before sw data starts
      if (t < inputData[inputDataIndex][0][0]) {
         i1 = i2 = 0;
         s = 0;
      } else {
         for (uint i=0; i<inputData[inputDataIndex].size(); i++) {
            if (inputData[inputDataIndex][i][0] >= t) {
               found = true;
               i2 = (int)i;
               break;
            }
         }
         if (found) {
            // i2 is now "ceil(t)"
            i1 = i2 - 1;
            if (i1 < 0) {
               i1 = i2 = 0;
               s = 0.0;
            } else {
               // normal case, now both i1 and i2 are >= 0 and < nlines, and i1 = i2-1
               s = (t - inputData[inputDataIndex][i1][0])/(inputData[inputDataIndex][i2][0] - inputData[inputDataIndex][i1][0]);
            }
         } else {
            i1 = i2 = inputData[inputDataIndex].size()-1;
            s = 0.0;
         }
      }
      
      creal s1 = 1 - s;
      
      for(uint i=0; i<nParams-1; i++) {
         outputData[i] = s1*inputData[inputDataIndex][i1][i+1] +
                           s*inputData[inputDataIndex][i2][i+1];
      }
   }

   void SetByUser::generateTemplateCell(
      spatial_cell::SpatialCell& templateCell,
      int inputDataIndex,
      creal& t
   ) {
      cerr << "Base class SetByUser::generateTemplateCell() called instead of derived class function!" << endl;
   }
   
   string SetByUser::getName() const {
      cerr << "Base class SetByUser::getName() called instead of derived class function!" << endl;
      return "SetByUser";
   }
   
   uint SetByUser::getIndex() const {
      cerr << "Base class SetByUser::getIndex() called instead of derived class function!" << endl;
      return sysboundarytype::N_SYSBOUNDARY_CONDITIONS;
   }
}
