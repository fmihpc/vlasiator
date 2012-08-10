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

#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <limits>

#include "sysboundarycondition.h"
#include "solarwind.h"
#include "../parameters.h"

using namespace std;

namespace SBC {
   SolarWind::SolarWind(): SysBoundaryCondition() { }
   SolarWind::~SolarWind() { }
   
   bool SolarWind::initSysBoundary(creal& t) {
      /* The bit field encodes which of the +x, -x, +y, -y, +z, -z faces are to have solar wind system boundary conditions.
       * The bit field has length 6, a bit raised to 1 indicates the corresponding face will have outflow.
       * The 6 bits left-to-right correspond to +x, -x, +y, -y, +z, -z respectively.
       */
      bool success = true;
      numberOfFaces = Parameters::solarWindFaceList.size();
      faces = 0;
      vector<string>::const_iterator it;
      for (it = Parameters::solarWindFaceList.begin();
           it != Parameters::solarWindFaceList.end();
      it++) {
         if(*it == "x+") faces = faces|(1<<5);
         if(*it == "x-") faces = faces|(1<<4);
         if(*it == "y+") faces = faces|(1<<3);
         if(*it == "y-") faces = faces|(1<<2);
         if(*it == "z+") faces = faces|(1<<1);
         if(*it == "z-") faces = faces|1;
      }
      
      success = loadInputData();
      success = success & generateTemplateCells(t);
      
      return success;
   }
   
   int SolarWind::assignSysBoundary(creal* cellParams) {
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5*dx;
      creal y = cellParams[CellParams::YCRD] + 0.5*dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
      
      int typeToAssign = sysboundarytype::NOT_SYSBOUNDARY;
      determineFace(x, y, z, dx, dy, dz);
      // Bitwise comparison of the field defining which faces to use and the field telling on which faces this cell is
      if((faces & isThisCellOnAFace) != 0) typeToAssign = getIndex();
      
      return typeToAssign;
   }
   
   Real SolarWind::calcPhaseSpaceDensity(creal& x,creal& y,creal& z,
                                         creal& dx,creal& dy,creal& dz,
                                         creal& vx,creal& vy,creal& vz,
                                         creal& dvx,creal& dvy,creal& dvz
   ) {
      Real average = 0.0;
      determineFace(x, y, z, dx, dy, dz);
      for(uint i=0; i<6; i++) {
         if((isThisCellOnAFace&(1<<(5-i)))==(1<<(5-i))) {
            average = templateCells[faceToInputData.find(i)->second].get_value(vx, vy, vz);
            break;
         }
      }
      return average;
   }
   
   void SolarWind::calcCellParameters(Real* cellParams, creal& t) {
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5*dx;
      creal y = cellParams[CellParams::YCRD] + 0.5*dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
      
      determineFace(x, y, z, dx, dy, dz);
      
      for(uint i=0; i<6; i++) {
         if((isThisCellOnAFace&(1<<(5-i)))==(1<<(5-i))) {
            cellParams[CellParams::BX] = templateCells[i].parameters[CellParams::BX];
            cellParams[CellParams::BY] = templateCells[i].parameters[CellParams::BY];
            cellParams[CellParams::BZ] = templateCells[i].parameters[CellParams::BZ];
            break;
         }
      }
   }
   
   bool SolarWind::loadInputData() {
      for(uint i=0; i<6; i++) {
         if((faces&(1<<(5-i))) == true) {
            inputData.push_back(loadFile(&(Parameters::solarWindFiles[i][0])));
            faceToInputData[i] = inputData.size();
         } else {
            faceToInputData[i] = -1;
         }
      }
      return true;
   }
   
   /*! Load solar wind boundary data from given file with the following format (in SI units):
    * t n T vx vy vz Bx By Bz
    * where t is the time of data and n the number of protons per m^3.
    * 
    * Function adapted from GUMICS-5.
    */
   vector<Real*> SolarWind::loadFile(const char *fn)
   {
      vector<Real*> dataset;
      
      int myrank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
      
      // Count lines with data
      FILE *fp;
      fp = fopen(fn,"r");
      if (fp == NULL) {
         std::cerr << "Couldn't open solar wind file " << fn << std::endl;
         exit(1);
      }
      int nlines = 0;
      int ret = 9;
      while (!feof(fp) && ret == 9) {
         Real x;
         ret = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &x, &x, &x, &x, &x, &x, &x, &x, &x);
         nlines++;
      }
      nlines--;
      fclose(fp);
      
      if (nlines < 1) {
         std::cerr << "Solar wind file must have at least one value (t, n, T...)" << std::endl;
         exit(1);
      }
      
      if (myrank == 0) std::cout << "Solar wind data file (" << fn << ") has " << nlines << " values"<< std::endl;
      
      Real tempData[9];
      
      fp = fopen(fn,"r");
      for (int line=0; line<nlines; line++) {
         for (int i=0; i<9; i++) {
            Real x;
            int ret = fscanf(fp,"%lf",&x);
            if (ret != 1) {
               std::cerr << "Couldn't read a number from solar wind file " << *fn << " for solar wind value " << line << std::endl;
            }
            tempData[i] = x;
         }
         dataset.push_back(tempData);
      }
      
      // check that sw data is in ascending temporal order
      for (int line = 1; line < nlines; line++) {
         if (dataset[line][0] < dataset[line - 1][0]) {
            std::cerr << "Solar wind data must be in ascending temporal order" << std::endl;
            exit(1);
         }
      }
      
      fclose(fp);
      
      return dataset;
   }
   
   bool SolarWind::generateTemplateCells(creal& t) {
      for(uint i=0; i<numberOfFaces; i++) {
         int index;
         if((index=faceToInputData.find(i)->second)!=-1) {
            generateTemplateCell(templateCells[i], index, t);
         }
      }
      return true;
   }
   
   void SolarWind::generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t) {
      Real rho, T, Vx, Vy, Vz, Bx, By, Bz;
      interpolate(inputDataIndex, t, &rho, &T, &Vx, &Vy, &Vz, &Bx, &By, &Bz);
      templateCell.parameters[CellParams::XCRD] = 0.0;
      templateCell.parameters[CellParams::YCRD] = 0.0;
      templateCell.parameters[CellParams::ZCRD] = 0.0;
      templateCell.parameters[CellParams::DX] = 1;
      templateCell.parameters[CellParams::DY] = 1;
      templateCell.parameters[CellParams::DZ] = 1;
      templateCell.parameters[CellParams::BX] = Bx;
      templateCell.parameters[CellParams::BY] = By;
      templateCell.parameters[CellParams::BZ] = Bz;
      
      templateCell.parameters[CellParams::RHO  ] = 0.0;
      templateCell.parameters[CellParams::RHOVX] = 0.0;
      templateCell.parameters[CellParams::RHOVY] = 0.0;
      templateCell.parameters[CellParams::RHOVZ] = 0.0;
      templateCell.parameters[CellParams::RHOLOSSADJUST] = 0.0;
      templateCell.parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
      
      // Go through each velocity block in the velocity phase space grid.
      // Set the initial state and block parameters:
      creal dvx_block = spatial_cell::SpatialCell::block_dvx; // Size of a block in vx-direction
      creal dvy_block = spatial_cell::SpatialCell::block_dvy; //                    vy
      creal dvz_block = spatial_cell::SpatialCell::block_dvz; //                    vz
      creal dvx_blockCell = spatial_cell::SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
      creal dvy_blockCell = spatial_cell::SpatialCell::cell_dvy; //                                vy
      creal dvz_blockCell = spatial_cell::SpatialCell::cell_dvz; //                                vz
      
      for (uint kv=0; kv<P::vzblocks_ini; ++kv) for (uint jv=0; jv<P::vyblocks_ini; ++jv) for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
         creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
         creal vy_block = P::vymin + jv*dvy_block; // vy-
         creal vz_block = P::vzmin + kv*dvz_block; // vz-
         
         // Calculate volume average of distrib. function for each cell in the block.
         for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic)
         {
            creal vx_cell = vx_block + ic*dvx_blockCell;
            creal vy_cell = vy_block + jc*dvy_blockCell;
            creal vz_cell = vz_block + kc*dvz_blockCell;
            creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
            creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
            creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
            creal average = rho * pow(physicalconstants::MASS_PROTON / 
                            (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
                            exp(-physicalconstants::MASS_PROTON *
                            (pow(vx_cell_center - Vx, 2.0) + 
                             pow(vy_cell_center - Vy, 2.0) + 
                             pow(vz_cell_center - Vz, 2.0)) / 
                             (2.0 * physicalconstants::K_B * T));
            if(average!=0.0){
               templateCell.set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
               // Add contributions to spatial cell velocity moments:
               creal dV = dvx_blockCell*dvy_blockCell*dvz_blockCell;  // Volume of one cell in a block      
               templateCell.parameters[CellParams::RHO  ] += average*dV;
               templateCell.parameters[CellParams::RHOVX] += average*vx_cell_center*dV;
               templateCell.parameters[CellParams::RHOVY] += average*vy_cell_center*dV;
               templateCell.parameters[CellParams::RHOVZ] += average*vz_cell_center*dV;
            }
         }
      }
   }
   
   void SolarWind::interpolate(const int inputDataIndex, creal t, Real* rho, Real* T, Real* vx, Real* vy, Real* vz, Real* Bx, Real* By, Real* Bz) {
      // Find first data[0] value which is >= t
      uint i;
      int i1=0,i2=0;
      bool found = false;
      Real s;      // 0 <= s < 1
      
      // use first value of sw data if interpolating for time before sw data starts
      if (t < inputData[inputDataIndex][0][0]) {
         i1 = i2 = 0;
         s = 0;
      } else {
         for (i=0; i<inputData[inputDataIndex].size(); i++) {
            if (inputData[inputDataIndex][i][0] >= t) {
               found = true;
               i2 = i;
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
      *rho  = s1*inputData[inputDataIndex][i1][1] + s*inputData[inputDataIndex][i2][1];
      *T  = s1*inputData[inputDataIndex][i1][2] + s*inputData[inputDataIndex][i2][2];
      *vx = s1*inputData[inputDataIndex][i1][3] + s*inputData[inputDataIndex][i2][3];
      *vy = s1*inputData[inputDataIndex][i1][4] + s*inputData[inputDataIndex][i2][4];
      *vz = s1*inputData[inputDataIndex][i1][5] + s*inputData[inputDataIndex][i2][5];
      *Bx = s1*inputData[inputDataIndex][i1][6] + s*inputData[inputDataIndex][i2][6];
      *By = s1*inputData[inputDataIndex][i1][7] + s*inputData[inputDataIndex][i2][7];
      *Bz = s1*inputData[inputDataIndex][i1][8] + s*inputData[inputDataIndex][i2][8];
   }
   
   std::string SolarWind::getName() const {return "SolarWind";}
   
   void SolarWind::determineFace(creal x, creal y, creal z, creal dx, creal dy, creal dz) {
      isThisCellOnAFace = 0;
      if(x > Parameters::xmax - dx) isThisCellOnAFace = isThisCellOnAFace|(1<<5);
      if(x < Parameters::xmin + dx) isThisCellOnAFace = isThisCellOnAFace|(1<<4);
      if(y > Parameters::ymax - dy) isThisCellOnAFace = isThisCellOnAFace|(1<<3);
      if(y < Parameters::ymin + dy) isThisCellOnAFace = isThisCellOnAFace|(1<<2);
      if(z > Parameters::zmax - dz) isThisCellOnAFace = isThisCellOnAFace|(1<<1);
      if(z < Parameters::zmin + dz) isThisCellOnAFace = isThisCellOnAFace|1;
   }
   
   uint SolarWind::getIndex() const {return sysboundarytype::SW;}
   uint SolarWind::getPrecedence() const {return Parameters::solarWindPrecedence;}
}
