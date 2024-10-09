#pragma once
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
#include <vector>
#include <array>
#include "boundaries.h"
#include "particleparameters.h"

// A 3D cartesian vector field with suitable interpolation properties for
// particle pushing
struct Field
{
   // Time at which this field is "valid"
   double time;

   // Mesh spacing
   double dx[3];

   // Information about spatial dimensions (like extent, boundaries, etc)
   Boundary* dimension[3];

   // The actual field data
   std::vector<double> data;

   // Constructor (primarily here to make sure boundaries are properly initialized as zero)
   Field() {
      for(int i=0; i<3; i++) {
         dimension[i] = nullptr;
      }
   }

   double* getCellRef(int x, int y, int z) {

      if(dimension[2]->cells == 1) {
         // Equatorial plane
         return &(data[4*(y*dimension[0]->cells+x)]);
      } else {
         // General 3d case
         return &(data[4*(z*dimension[0]->cells*dimension[1]->cells + y*dimension[0]->cells + x)]);
      }
   }

   std::array<double,3> getCell(int x, int y, int z) {

      // Map these cell coordinates using the boundaries
      x = dimension[0]->cellCoordinate(x);
      y = dimension[1]->cellCoordinate(y);
      z = dimension[2]->cellCoordinate(z);

      double* cell = getCellRef(x,y,z);
      return {cell[0],cell[1],cell[2]};
   }

   // Round-Brace indexing: indexing by physical location, with interpolation
   virtual std::array<double,3> operator()(std::array<double,3> v) {
      double min[3] = { min[0] = dimension[0]->min, dimension[1]->min, dimension[2]->min};

      int32_t index[3];
      double fract[3];

      for(int i=0; i<3; i++) {
         v[i] -= min[i];
         v[i] /= dx[i];
         index[i] = (int32_t)std::trunc(v[i]);
         fract[i] = v[i] - (double)std::trunc(v[i]);
      }

      if(dimension[2]->cells <= 1) {
         // Equatorial plane
         std::vector<double> interp(12);
         std::array<double,3> tmp;
         tmp = getCell(index[0],index[1],index[2]);
         interp[0] = tmp[0]; // former interp 0
         interp[1] = tmp[1];
         interp[2] = tmp[2];
         tmp = getCell(index[0]+1,index[1],index[2]);
         interp[3] = tmp[0]; // former interp 1
         interp[4] = tmp[1];
         interp[5] = tmp[2];
         tmp = getCell(index[0],index[1]+1,index[2]);
         interp[6] = tmp[0]; // former interp 2
         interp[7] = tmp[1];
         interp[8] = tmp[2];
         tmp = getCell(index[0]+1,index[1]+1,index[2]);
         interp[9] = tmp[0]; // former interp 3
         interp[10] = tmp[1];
         interp[11] = tmp[2];

         tmp[0] = fract[0]*(fract[1]*interp[9]+(1.-fract[1])*interp[3])
            + (1.-fract[0])*(fract[1]*interp[6]+(1.-fract[1])*interp[0]);
         tmp[1] = fract[0]*(fract[1]*interp[10]+(1.-fract[1])*interp[4])
            + (1.-fract[0])*(fract[1]*interp[7]+(1.-fract[1])*interp[1]);
         tmp[2] = fract[0]*(fract[1]*interp[11]+(1.-fract[1])*interp[5])
            + (1.-fract[0])*(fract[1]*interp[8]+(1.-fract[1])*interp[2]);


         return tmp;
      } else if (dimension[1]->cells <= 1) {
         // Polar plane
         std::vector<double> interp(12);
         std::array<double,3> tmp;

         tmp = getCell(index[0],index[1],index[2]);
         interp[0] = tmp[0]; // former interp 0
         interp[1] = tmp[1];
         interp[2] = tmp[2];
         tmp = getCell(index[0]+1,index[1],index[2]);
         interp[3] = tmp[0]; // former interp 1
         interp[4] = tmp[1];
         interp[5] = tmp[2];
         tmp = getCell(index[0],index[1],index[2]+1);
         interp[6] = tmp[0]; // former interp 2
         interp[7] = tmp[1];
         interp[8] = tmp[2];
         tmp = getCell(index[0]+1,index[1],index[2]+1);
         interp[9] = tmp[0]; // former interp 3
         interp[10] = tmp[1];
         interp[11] = tmp[2];

         tmp[0] = fract[0]*(fract[2]*interp[9]+(1.-fract[2])*interp[3])
            + (1.-fract[0])*(fract[2]*interp[6]+(1.-fract[2])*interp[0]);
         tmp[1] = fract[0]*(fract[2]*interp[10]+(1.-fract[2])*interp[4])
            + (1.-fract[0])*(fract[2]*interp[7]+(1.-fract[2])*interp[1]);
         tmp[2] = fract[0]*(fract[2]*interp[11]+(1.-fract[2])*interp[5])
            + (1.-fract[0])*(fract[2]*interp[8]+(1.-fract[2])*interp[2]);
         return tmp;
      } else {
         // Proper 3D
         std::vector<double> interp(24);
         std::array<double,3> tmp;
         tmp = getCell(index[0],index[1],index[2]);
         interp[0] = tmp[0]; // former interp 0
         interp[1] = tmp[1];
         interp[2] = tmp[2];
         tmp = getCell(index[0]+1,index[1],index[2]);
         interp[3] = tmp[0]; // former interp 1
         interp[4] = tmp[1];
         interp[5] = tmp[2];
         tmp = getCell(index[0],index[1]+1,index[2]);
         interp[6] = tmp[0]; // former interp 2
         interp[7] = tmp[1];
         interp[8] = tmp[2];
         tmp = getCell(index[0]+1,index[1]+1,index[2]);
         interp[9] = tmp[0]; // former interp 3
         interp[10] = tmp[1];
         interp[11] = tmp[2];
         tmp = getCell(index[0],index[1],index[2]+1);
         interp[12] = tmp[0]; // former interp 4
         interp[13] = tmp[1];
         interp[14] = tmp[2];
         tmp = getCell(index[0]+1,index[1],index[2]+1);
         interp[15] = tmp[0]; // former interp 5
         interp[16] = tmp[1];
         interp[17] = tmp[2];
         tmp = getCell(index[0],index[1]+1,index[2]+1);
         interp[18] = tmp[0]; // former interp 6
         interp[19] = tmp[1];
         interp[20] = tmp[2];
         tmp = getCell(index[0]+1,index[1]+1,index[2]+1);
         interp[21] = tmp[0]; // former interp 7
         interp[22] = tmp[1];
         interp[23] = tmp[2];



         tmp[0] = fract[2] * (
               fract[0]*(fract[1]*interp[9]+(1.-fract[1])*interp[3])
               + (1.-fract[0])*(fract[1]*interp[6]+(1.-fract[1])*interp[0]))
            + (1.-fract[2]) * (
                  fract[0]*(fract[1]*interp[21]+(1.-fract[1])*interp[15])
                  + (1.-fract[0])*(fract[1]*interp[18]+(1.-fract[1])*interp[12]));
         tmp[1] = fract[2] * (
               fract[0]*(fract[1]*interp[10]+(1.-fract[1])*interp[4])
               + (1.-fract[0])*(fract[1]*interp[7]+(1.-fract[1])*interp[1]))
            + (1.-fract[2]) * (
                  fract[0]*(fract[1]*interp[22]+(1.-fract[1])*interp[16])
                  + (1.-fract[0])*(fract[1]*interp[19]+(1.-fract[1])*interp[13]));
         tmp[2] = fract[2] * (
               fract[0]*(fract[1]*interp[11]+(1.-fract[1])*interp[5])
               + (1.-fract[0])*(fract[1]*interp[8]+(1.-fract[1])*interp[2]))
            + (1.-fract[2]) * (
                  fract[0]*(fract[1]*interp[23]+(1.-fract[1])*interp[17])
                  + (1.-fract[0])*(fract[1]*interp[20]+(1.-fract[1])*interp[14]));
         return tmp;
      }

   }
   virtual std::array<double,3> operator()(double x, double y, double z) {
      std::array<double,3> v = {x,y,z};
      return operator()(v);
   }

};

// Linear Temporal interpolation between two input fields
struct Interpolated_Field : Field{
   Field &a, &b;
   double t;

   /* Constructor:
    * Inputs are the two fields to interpolate between
    * and the current time.
    */
   Interpolated_Field(Field& _a, Field& _b, float _t) : a(_a),b(_b),t(_t) {
   }

   virtual std::array<double,3> operator()(std::array<double,3> v) {
      std::array<double,3> aval=a(v);
      std::array<double,3> bval=b(v);

      double fract = (t - a.time)/(b.time-a.time);
      std::array<double,3> tmp;
      tmp[0] = fract*bval[0] + (1.-fract)*aval[0];
      tmp[1] = fract*bval[1] + (1.-fract)*aval[1];
      tmp[2] = fract*bval[2] + (1.-fract)*aval[2];
      return tmp;
   }
};
