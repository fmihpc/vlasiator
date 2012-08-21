/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012 Finnish Meteorological Institute
 * 
 * Vlasiator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 * 
 * Vlasiator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/*!\file setbyuser.cpp
 * \brief Implementation of the class SysBoundaryCondition::SetByUser. This serves as the base class for further classes like SysBoundaryCondition::SetMaxwellian.
 */

#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <limits>

#include "setbyuser.h"
#include "../vlasovmover.h"


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
   
   bool SetByUser::initSysBoundary(creal& t) {
      /* The array of bool describes which of the x+, x-, y+, y-, z+, z- faces are to have user-set system boundary conditions.
       * A true indicates the corresponding face will have user-set system boundary conditions.
       * The 6 elements correspond to x+, x-, y+, y-, z+, z- respectively.
       */
      bool success = true;
      for(uint i=0; i<6; i++) facesToProcess[i] = false;
      
      this->getParameters();
      
      vector<string>::const_iterator it;
      for (it = faceList.begin();
           it != faceList.end();
      it++) {
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
   
   int SetByUser::assignSysBoundary(creal* cellParams) {
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5*dx;
      creal y = cellParams[CellParams::YCRD] + 0.5*dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
      
      int typeToAssign = sysboundarytype::NOT_SYSBOUNDARY;
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);
      // Comparison of the array defining which faces to use and the array telling on which faces this cell is
      bool doAssign = false;
      for(uint i=0; i<6; i++) doAssign = doAssign || (facesToProcess[i] && isThisCellOnAFace[i]);
      if(doAssign) typeToAssign = this->getIndex();
      
      return typeToAssign;
   }
   
   bool SetByUser::applyInitialState(dccrg::Dccrg<SpatialCell>& mpiGrid) {
      vector<uint64_t> cells = mpiGrid.get_cells();
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if(cell->sysBoundaryFlag != this->getIndex()) continue;
         
         creal dx = cell->parameters[CellParams::DX];
         creal dy = cell->parameters[CellParams::DY];
         creal dz = cell->parameters[CellParams::DZ];
         creal x = cell->parameters[CellParams::XCRD] + 0.5*dx;
         creal y = cell->parameters[CellParams::YCRD] + 0.5*dy;
         creal z = cell->parameters[CellParams::ZCRD] + 0.5*dz;
         
         determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);
         
         for(uint i=0; i<6; i++) {
            if(facesToProcess[i] && isThisCellOnAFace[i]) {
               cell->parameters[CellParams::BX] = templateCells[i].parameters[CellParams::BX];
               cell->parameters[CellParams::BY] = templateCells[i].parameters[CellParams::BY];
               cell->parameters[CellParams::BZ] = templateCells[i].parameters[CellParams::BZ];
               
               cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
               cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
               
               // Go through each velocity block in the velocity phase space grid.
               // Set the initial state and block parameters:
               creal dvx_block = SpatialCell::block_dvx; // Size of a block in vx-direction
               creal dvy_block = SpatialCell::block_dvy; //                    vy
               creal dvz_block = SpatialCell::block_dvz; //                    vz
               creal dvx_blockCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
               creal dvy_blockCell = SpatialCell::cell_dvy; //                                vy
               creal dvz_blockCell = SpatialCell::cell_dvz; //                                vz
               
               for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
                  for (uint jv=0; jv<P::vyblocks_ini; ++jv)
                     for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
                        creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
                        creal vy_block = P::vymin + jv*dvy_block; // vy-
                        creal vz_block = P::vzmin + kv*dvz_block; // vz-
                        
                        // Calculate volume average of distrib. function for each cell in the block.
                        for (uint kc=0; kc<WID; ++kc) 
                           for (uint jc=0; jc<WID; ++jc) 
                              for (uint ic=0; ic<WID; ++ic) {
                                 creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
                                 creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
                                 creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
                                 Real average = templateCells[i].get_value(vx_cell_center,
                                                                           vy_cell_center,
                                                                           vz_cell_center);
                                 
                                 if(average!=0.0){
                                    cell->set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
                                 }
                              }
                     }
               calculateCellVelocityMoments(cell);
               
               //let's get rid of blocks not fulfilling the criteria here to save memory.
               cell->adjustSingleCellVelocityBlocks();
            }
         }
      }
      return true;
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
    * \retval dataset Vector of Real vectors. Each line of length nParams is put into
    * a vector. Each of these is then put into the vector returned here.
    */
   vector<vector<Real> > SetByUser::loadFile(const char *fn) {
      vector<vector<Real> > dataset;
      
      int myrank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
      
      // Count lines with data
      FILE *fp;
      fp = fopen(fn,"r");
      if (fp == NULL) {
         cerr << "Couldn't open parameter file " << fn << endl;
         exit(1);
      }
      uint nlines = 0;
      int ret = nParams;
      while (!feof(fp) && ret == (int)nParams) {
         Real x;
         ret = 0;
         for(uint i=0; i<nParams; i++) ret += fscanf(fp, "%lf", &x);
         nlines++;
      }
      nlines--;
      
      fclose(fp);
      
      if (nlines < 1) {
         cerr << "Parameter file must have at least one value (t, n, T...)" << endl;
         exit(1);
      }
      
      if (myrank == 0) cout << "Parameter data file (" << fn << ") has " << nlines << " values"<< endl;
      
      fp = fopen(fn,"r");
      for (uint line=0; line<nlines; line++) {
         vector<Real> tempData;
         for (uint i=0; i<nParams; i++) {
            Real x;
            int ret = fscanf(fp,"%lf",&x);
            if (ret != 1) {
               cerr << "Couldn't read a number from parameter file " << *fn << " for line value " << line << endl;
            }
            tempData.push_back(x);
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
   void SetByUser::interpolate(const int inputDataIndex, creal t, Real* outputData) {
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

   void SetByUser::generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t) {
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
