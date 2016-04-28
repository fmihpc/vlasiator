#pragma once
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

   // Coordinate boundaries
   double min[3];
   double max[3];

   // Mesh spacing
   double dx[3];

   // Mesh cells
   int cells[3];

   bool periodic[3];

   // The actual field data
   std::vector<double> data;

   double* getCellRef(int x, int y, int z) {

      if(cells[2] == 1) {
         return &(data[4*(y*cells[0]+x)]);
      } else {
         return &(data[4*(z*cells[0]*cells[1] + y*cells[0] + x)]);
      }
   }

   Vec3d getCell(int x, int y, int z) {

      // Map these cell coordinates using the boundaries
      x = ParticleParameters::boundary_behaviour_x->cellCoordinate(x);
      y = ParticleParameters::boundary_behaviour_y->cellCoordinate(y);
      z = ParticleParameters::boundary_behaviour_z->cellCoordinate(z);

      double* cell = getCellRef(x,y,z);
      Vec3d retval;
      retval.load_partial(3,cell);
      return retval;
   }

   // Round-Brace indexing: indexing by physical location, with interpolation
   Vec3d operator()(Vec3d v) {
      Vec3d vmin,vdx;
      vmin.load(min);
      vdx.load(dx);

      v -= vmin;
      v /= vdx;

      int index[3];
      double fract[3];
      truncate_to_int(v).store(index);
      (v-truncate(v)).store(fract);

      if(cells[2] <= 1) {
         // Equatorial plane
         Vec3d interp[4];
         interp[0] = getCell(index[0],index[1],index[2]);
         interp[1] = getCell(index[0]+1,index[1],index[2]);
         interp[2] = getCell(index[0],index[1]+1,index[2]);
         interp[3] = getCell(index[0]+1,index[1]+1,index[2]);

         return fract[0]*(fract[1]*interp[3]+(1.-fract[1])*interp[1])
            + (1.-fract[0])*(fract[1]*interp[2]+(1.-fract[1])*interp[0]);
      } else if (cells[1] <= 1) {
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
   Field& a,b;
   double t;

   /* Constructor:
    * Inputs are the two fields to interpolate between
    * and the current time.
    */
   Interpolated_Field(Field& _a, Field& _b, float _t) : a(_a),b(_b),t(_t) {
   }

   Vec3d operator()(Vec3d v) {
      Vec3d aval=a(v);
      Vec3d bval=b(v);

      double fract = (t - a.time)/(b.time-a.time);
      return fract*bval + (1.-fract)*aval;
   }
};
