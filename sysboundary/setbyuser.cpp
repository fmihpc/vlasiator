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

/*!\file setbyuser.cpp
 * \brief Implementation of the class SysBoundaryCondition::SetByUser. 
 * This serves as the base class for further classes like SysBoundaryCondition::SetMaxwellian.
 */

#include <cstdlib>
#include <iostream>

#include <assert.h>
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
   SetByUser::SetByUser(): OuterBoundaryCondition() { }
   SetByUser::~SetByUser() { }
   
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
      
      for(unsigned int i=0; i<speciesParams.size(); i++) {
         success = loadInputData(i);
      }
      success = success & generateTemplateCells(t);
      
      return success;
   }
   
   bool SetByUser::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      Project &project
   ) {
      bool success = true;
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         if (!setCellsFromTemplate(mpiGrid, popID)) success = false;
      }
      if (!setBFromTemplate(technicalGrid, perBGrid)) success = false;
      
      return success;
   }
   
   Real SetByUser::fieldSolverBoundaryCondMagneticField(
      FsGrid< array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k,
      creal& dt,
      cuint& component
   ) {
      Real result = 0.0;
      const array<int, 3> globalIndices = technicalGrid.getGlobalIndices(i,j,k);

      creal x = (convert<Real>(globalIndices[0])+0.5)*technicalGrid.DX + Parameters::xmin;
      creal y = (convert<Real>(globalIndices[1])+0.5)*technicalGrid.DY + Parameters::ymin;
      creal z = (convert<Real>(globalIndices[2])+0.5)*technicalGrid.DZ + Parameters::zmin;
      int refLevel = technicalGrid.get(i, j, k)->refLevel;

      // if refLevel isn't 0, assume neighbour might be on a lower refinement level
      if (refLevel > 0) {
         --refLevel;
      }

      creal dx = Parameters::dx_ini * pow(2, -refLevel);
      creal dy = Parameters::dy_ini * pow(2, -refLevel);
      creal dz = Parameters::dz_ini * pow(2, -refLevel);
      
      bool isThisCellOnAFace[6];
      determineFace(isThisCellOnAFace, x, y, z, dx, dy, dz, true);

      for (uint i=0; i<6; i++) {
         if (isThisCellOnAFace[i]) {
            result = templateB[i][component];
            break; // This effectively sets the precedence of faces through the order of faces.
         }
      }
      return result;
   }

   void SetByUser::fieldSolverBoundaryCondElectricField(
      FsGrid< array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGrid.get(i,j,k)->at(fsgrids::efield::EX+component) = 0.0;
   }

   void SetByUser::fieldSolverBoundaryCondHallElectricField(
      FsGrid< array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      array<Real, fsgrids::ehall::N_EHALL> * cp = EHallGrid.get(i,j,k);
      switch (component) {
         case 0:
            cp->at(fsgrids::ehall::EXHALL_000_100) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_010_110) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_001_101) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_011_111) = 0.0;
            break;
         case 1:
            cp->at(fsgrids::ehall::EYHALL_000_010) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_100_110) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_001_011) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_101_111) = 0.0;
            break;
         case 2:
            cp->at(fsgrids::ehall::EZHALL_000_001) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_100_101) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_010_011) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_110_111) = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   void SetByUser::fieldSolverBoundaryCondGradPeElectricField(
      FsGrid< array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
         EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE+component) = 0.0;
   }
   
   void SetByUser::fieldSolverBoundaryCondDerivatives(
      FsGrid< array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      FsGrid< array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint& RKCase,
      cuint& component
   ) {
      this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
   }

   void SetByUser::fieldSolverBoundaryCondBVOLDerivatives(
      FsGrid< array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {
      this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   }

   void SetByUser::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const uint popID,
      const bool calculate_V_moments
   ) {
      // No need to do anything in this function, as the propagators do not touch the distribution function   
   }
   
   bool SetByUser::setBFromTemplate(FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid, FsGrid< array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid) {

      array<bool,6> isThisCellOnAFace;
      const array<int, 3> gridDims(perBGrid.getLocalSize());

      for (int k=0; k<gridDims[2]; k++) {
         for (int j=0; j<gridDims[1]; j++) {
            for (int i=0; i<gridDims[0]; i++) {
               if (technicalGrid.get(i, j, k)->sysBoundaryFlag != this->getIndex())
                  continue;

               const auto coords = technicalGrid.getPhysicalCoords(i,j,k);
               
               // TODO: This code up to determineFace() should be in a separate function, it gets called in a lot of places.
               // Shift to the center of the fsgrid cell
               auto cellCenterCoords = coords;
               cellCenterCoords[0] += 0.5 * perBGrid.DX;
               cellCenterCoords[1] += 0.5 * perBGrid.DY;
               cellCenterCoords[2] += 0.5 * perBGrid.DZ;

               int refLevel = technicalGrid.get(i, j, k)->refLevel;

               if(refLevel == -1) {
                  cerr << "Error, could not get refinement level of remote DCCRG cell " << __FILE__ << " " << __LINE__ << endl;
                  return false;
               }

               // if refLevel isn't 0, assume neighbour might be on a lower refinement level
               if (refLevel > 0)
                  --refLevel;

               creal dx = P::dx_ini / pow(2, refLevel);
               creal dy = P::dy_ini / pow(2, refLevel);
               creal dz = P::dz_ini / pow(2, refLevel);
               
               isThisCellOnAFace.fill(false);

               determineFace(isThisCellOnAFace.data(), cellCenterCoords[0], cellCenterCoords[1], cellCenterCoords[2], dx, dy, dz);

               for(uint iface=0; iface < 6; iface++) {
                  if(facesToProcess[iface] && isThisCellOnAFace[iface]) {
                     perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBX) = templateB[iface][0];
                     perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBY) = templateB[iface][1];
                     perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBZ) = templateB[iface][2];
                     break;
                  }
               }
            }
         }
      }
      return true;
   }


   bool SetByUser::setCellsFromTemplate(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const uint popID) {

      //#pragma omp parallel for
      for (const auto& id : getLocalCells()) {
         SpatialCell* cell = mpiGrid[id];
         if(cell->sysBoundaryFlag != this->getIndex()) 
            continue;
         
         std::array<bool, 6> isThisCellOnAFace;
         determineFace(isThisCellOnAFace, mpiGrid, id);

         for (uint i=0; i < 6; i++) {
            if (facesToProcess[i] && isThisCellOnAFace[i]) {
               copyCellData(&templateCells[i], cell ,false,popID,true); // copy also vdf, _V
               copyCellData(&templateCells[i], cell ,true,popID,false); // don't copy vdf again but copy _R now
               break; // This effectively sets the precedence of faces through the order of faces.
            }
         }
      }

      return true;
   }
   
   void SetByUser::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) faces[i] = facesToProcess[i];
   }
   
   bool SetByUser::loadInputData(const uint popID) {
      UserSpeciesParameters& sP = speciesParams[popID];

      for(uint i=0; i<6; i++) {
         if(facesToProcess[i]) {
            sP.inputData[i] = loadFile(sP.files[i].c_str(), sP.nParams);
         } else {
            vector<Real> tmp1;
            vector<vector<Real> > tmp2;
            for(uint j=0; j<sP.nParams; j++) {
               tmp1.push_back(-1.0);
            }
            tmp2.push_back(tmp1);
            sP.inputData[i] = tmp2;
         }
      }
      return true;
   }
   
// Suppress fscanf warnings since we check type
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
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
   vector<vector<Real> > SetByUser::loadFile(const char *fn, unsigned int nParams) {
      vector<vector<Real> > dataset;
 
   
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      
      // Count lines with data
      FILE *fp;
      fp = fopen(fn,"r");
      if (fp == NULL) {
         cerr << "Couldn't open parameter file " << fn << endl;
         exit(ExitCodes::FAILURE);
      }
      uint nlines = 0;
      int ret = nParams;

      // Make sure the type id of Real is correct
      assert( typeid( Real ) == typeid(float) || typeid( Real ) == typeid(double) );

      while (!feof(fp) && ret == (int)nParams) {
         Real readParam;
         ret = 0;
         if ( typeid( readParam ) == typeid(double) ) {
            for(uint i=0; i<nParams; i++) {
               ret += fscanf(fp, "%lf", &readParam);
            }
         } else if( typeid( readParam ) == typeid(float) ) {
            for(uint i=0; i<nParams; i++) {
               ret += fscanf(fp, "%f", &readParam);
            }
         } else {
            assert( typeid( readParam ) == typeid(float) || typeid( readParam ) == typeid(double) );
         }
         nlines++;
      }
      nlines--;
      
      fclose(fp);
      
      if (nlines < 1) {
         cerr << "Parameter file must have at least one value (t, n, T...)" << endl;
         exit(ExitCodes::FAILURE);
      }
      
      if (myRank == 0) cout << "Parameter data file (" << fn << ") has " << nlines << " values"<< endl;
      
      fp = fopen(fn,"r");
      for (uint line=0; line<nlines; line++) {
         vector<Real> tempData;
         for (uint i=0; i<nParams; i++) {
            Real readParam;
            int ret = 0;
            if ( typeid( readParam ) == typeid(double) ) {
               ret = fscanf(fp,"%lf",&readParam);
            } else if( typeid( readParam ) == typeid(float) ) {
               ret = fscanf(fp,"%f",&readParam);
            } else {
               assert( typeid( readParam ) == typeid(float) || typeid( readParam ) == typeid(double) ); 
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
            exit(ExitCodes::FAILURE);
         }
      }
      
      fclose(fp);
      
      return dataset;
   }
#pragma GCC diagnostic pop
   
   /*! Loops through the array of template cells and generates the ones needed. The function
    * generateTemplateCell is defined in the inheriting class such as to have the specific
    * condition needed.
    * \param t Simulation time.
    * \sa generateTemplateCell
    */
   bool SetByUser::generateTemplateCells(creal& t) {
      #pragma omp parallel for
      for(uint i=0; i<6; i++) {
         if(facesToProcess[i]) {
            generateTemplateCell(templateCells[i], templateB[i], i, t);
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
      const int inputDataIndex, const uint popID,
      creal t,
      Real* outputData
   ) {

      UserSpeciesParameters& sP = speciesParams[popID];

      // Find first data[0] value which is >= t
      int i1=0,i2=0;
      bool found = false;
      Real s;      // 0 <= s < 1
      
      // use first value of sw data if interpolating for time before sw data starts
      if (t < sP.inputData[inputDataIndex][0][0]) {
         i1 = i2 = 0;
         s = 0;
      } else {
         for (uint i=0; i<sP.inputData[inputDataIndex].size(); i++) {
            if (sP.inputData[inputDataIndex][i][0] >= t) {
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
               s = (t - sP.inputData[inputDataIndex][i1][0])/(sP.inputData[inputDataIndex][i2][0] - sP.inputData[inputDataIndex][i1][0]);
            }
         } else {
            i1 = i2 = sP.inputData[inputDataIndex].size()-1;
            s = 0.0;
         }
      }
      
      creal s1 = 1 - s;
      
      for(uint i=0; i<sP.nParams-1; i++) {
         outputData[i] = s1*sP.inputData[inputDataIndex][i1][i+1] +
                           s*sP.inputData[inputDataIndex][i2][i+1];
      }
   }
}
