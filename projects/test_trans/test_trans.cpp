/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"

#include "test_trans.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
   test_trans::test_trans(): Project() { }
   test_trans::~test_trans() { }
      
  // Real this->cellPosition = 0;

   bool test_trans::initialize(void) {
      return Project::initialize();
   }

   void test_trans::addParameters(){
      typedef Readparameters RP;
      RP::add("test_trans.cellPosition", "Position of the centre of the cells initiated (same used in velocity and space).",(Real)1.5);
      RP::add("test_trans.peakValue","Value of the distribution function",(Real)1.0);
   }

   void test_trans::getParameters() {
      Project::getParameters();
      
      typedef Readparameters RP;
      RP::get("test_trans.cellPosition", this->cellPosition);
      RP::get("test_trans.peakValue" ,peakValue);
   }


   Real test_trans::calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
                                          creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz,
                                          const uint popID) const {
      //Please use even number of cells in velocity and real space
      Real xyz[3];
      Real vxyz[3];
      
   //location of this cell
      vxyz[0]=(vx+0.5*dvx)/dvx;
      vxyz[1]=(vy+0.5*dvy)/dvy;
      vxyz[2]=(vz+0.5*dvz)/dvz;

         
      xyz[0]=(x+0.5*dx)/dx;
      xyz[1]=(y+0.5*dy)/dy;
      xyz[2]=(z+0.5*dz)/dz;

      creal pos = this->cellPosition;
      //real space coordinates of boxes
      //Assume an even number of spatial cells per grid dimension
      const Real box_real[8][3] = { { pos, pos, pos},
                                    {-pos, pos, pos},
                                    { pos,-pos, pos},
                                    { pos, pos,-pos},
                                    {-pos,-pos, pos},
                                    {-pos, pos,-pos},
                                    { pos,-pos,-pos},
                                    {-pos,-pos,-pos}};
      
      //velocity space coordinates of boxes in reduced units
      //there is always an even amount of velocity cells per dimension (assuming WID is even) 
      const Real box_vel[8][3] = { { pos, pos, pos},
                                 {-pos, pos, pos},
                                 { pos,-pos, pos},
                                 { pos, pos,-pos},
                                 {-pos,-pos, pos},
                                 {-pos, pos,-pos},
                                 { pos,-pos,-pos},
                                 {-pos,-pos,-pos}};
      
      
      for(int box=0;box<8;box++){
         bool outsideBox=false;
         for(int i=0;i<3;i++){
            if(xyz[i]<(box_real[box][i]-0.1) ||
               xyz[i]>(box_real[box][i]+0.1) ||
               vxyz[i]<(box_vel[box][i]-0.1) ||
               vxyz[i]>(box_vel[box][i]+0.1)){
               outsideBox=true;
               break;
            }
         }
         
         if (!outsideBox) {
            return peakValue;
         }
      }
      return 0.0;
   }

   void test_trans::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      Real* cellParams = cell->get_cell_parameters();
      cellParams[CellParams::EX   ] = 0.0;
      cellParams[CellParams::EY   ] = 0.0;
      cellParams[CellParams::EZ   ] = 0.0;
      cellParams[CellParams::PERBX   ] = 0.0;
      cellParams[CellParams::PERBY   ] = 0.0;
      cellParams[CellParams::PERBZ   ] = 0.0;
   }
   
   void test_trans::setCellBackgroundField(SpatialCell* cell) const {
      ConstantField bgField;
      bgField.initialize(0.0,0.0,1e-9);
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
   }

}// namespace projects

