/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011, 2012 Finnish Meteorological Institute
*/

#include <stdlib.h>
#include <math.h>
#include "dipole.hpp"
#include "../common.h"


void Dipole::initialize(const double moment)
{
   this->initialized = true;
   q_x=0.0;
   q_y=0.0;
   q_z=-moment;
   this->center_x = 0;
   this->center_y = 0;
   this->center_z = 0;
}


double Dipole::call( double x, double y, double z) const
{
   const double minimumR=1e-3*physicalconstants::R_E; //The dipole field is defined to be outside of Earth, and units are in meters     
   if(this->initialized==false)
      return 0.0;
   x-= center_x;
   y-= center_y;
   z-= center_z;
   double r2 = x*x + y*y + z*z;
   if(r2<minimumR*minimumR)
      r2=minimumR*minimumR;
   const double r = sqrt(r2);
   const double invr5 = 1.0/(r2*r2*r);

   
   switch (_fComponent) {
       case X:
          //Dipole Bx terms
          return (-(r2*q_x) + 3*q_x*x*x + 3*q_y*x*y
                  + 3*q_z*x*z)*invr5;
          break;
       case Y:
          //Dipole By terms   
          return (-(r2*q_y) + 3*q_y*y*y + 3*q_z*y*z
                  + 3*q_x*y*x)*invr5;
          break;
       case Z:
          //Dipole Bz terms   
          return (-(r2*q_z) + 3*q_z*z*z + 3*q_x*z*x
                  + 3*q_y*z*y)*invr5;          
          break;
   }
   return 0;	// dummy, but prevents gcc from yelling
}







