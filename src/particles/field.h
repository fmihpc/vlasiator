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
#include "vectorclass.h"
#include "vector3d.h"
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

   Vec3d getCell(int x, int y, int z) {

      // Map these cell coordinates using the boundaries
      x = dimension[0]->cellCoordinate(x);
      y = dimension[1]->cellCoordinate(y);
      z = dimension[2]->cellCoordinate(z);

      double* cell = getCellRef(x,y,z);
      Vec3d retval;
      retval.load_partial(3,cell);
      return retval;
   }

   // Round-Brace indexing: indexing by physical location, with interpolation
   virtual Vec3d operator()(Vec3d v) {
      Vec3d vmin,vdx;
      double min[3] = { min[0] = dimension[0]->min, dimension[1]->min, dimension[2]->min};
      vmin.load(min);
      vdx.load(dx);

      v -= vmin;
      v /= vdx;

      int index[3];
      double fract[3];
      truncate_to_int(v).store(index);
      (v-Vec3d(truncate(v))).store(fract);

      if(dimension[2]->cells <= 1) {
         // Equatorial plane
         Vec3d interp[4];
         interp[0] = getCell(index[0],index[1],index[2]);
         interp[1] = getCell(index[0]+1,index[1],index[2]);
         interp[2] = getCell(index[0],index[1]+1,index[2]);
         interp[3] = getCell(index[0]+1,index[1]+1,index[2]);

         return fract[0]*(fract[1]*interp[3]+(1.-fract[1])*interp[1])
            + (1.-fract[0])*(fract[1]*interp[2]+(1.-fract[1])*interp[0]);
      } else if (dimension[1]->cells <= 1) {
         // Polar plane
         Vec3d interp[4];

         interp[0] = getCell(index[0],index[1],index[2]);
         interp[1] = getCell(index[0]+1,index[1],index[2]);
         interp[2] = getCell(index[0],index[1],index[2]+1);
         interp[3] = getCell(index[0]+1,index[1],index[2]+1);

         return fract[0]*(fract[2]*interp[3]+(1.-fract[2])*interp[1])
            + (1.-fract[0])*(fract[2]*interp[2]+(1.-fract[2])*interp[0]);
      } else {
         // Proper 3D
         Vec3d interp[8];

         interp[0] = getCell(index[0],index[1],index[2]);
         interp[1] = getCell(index[0]+1,index[1],index[2]);
         interp[2] = getCell(index[0],index[1]+1,index[2]);
         interp[3] = getCell(index[0]+1,index[1]+1,index[2]);
         interp[4] = getCell(index[0],index[1],index[2]+1);
         interp[5] = getCell(index[0]+1,index[1],index[2]+1);
         interp[6] = getCell(index[0],index[1]+1,index[2]+1);
         interp[7] = getCell(index[0]+1,index[1]+1,index[2]+1);

         return fract[2] * (
               fract[0]*(fract[1]*interp[3]+(1.-fract[1])*interp[1])
               + (1.-fract[0])*(fract[1]*interp[2]+(1.-fract[1])*interp[0]))
            + (1.-fract[2]) * (
                  fract[0]*(fract[1]*interp[7]+(1.-fract[1])*interp[5])
                  + (1.-fract[0])*(fract[1]*interp[6]+(1.-fract[1])*interp[4]));
      }

   }
   virtual Vec3d operator()(double x, double y, double z) {
      Vec3d v(x,y,z);
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

   virtual Vec3d operator()(Vec3d v) {
      Vec3d aval=a(v);
      Vec3d bval=b(v);

      double fract = (t - a.time)/(b.time-a.time);
      return fract*bval + (1.-fract)*aval;
   }
};
