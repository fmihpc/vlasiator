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
/*
Background magnetic field class of Vlasiator.
*/

#include <stdlib.h>
#include <math.h>
#include "dipole.hpp"
#include "../common.h"

//tilt_angle is agains the z axis in the x-z plane. In radians*/
void Dipole::initialize(const double moment,const double center_x, const double center_y, const double center_z, const double tilt_angle=0){
   this->initialized = true;
   q[0]=-sin(tilt_angle)*moment;
   q[1]=0.0;
   q[2]=-cos(tilt_angle)*moment;
   center[0]=center_x;
   center[1]=center_y;
   center[2]=center_z;
}



double Dipole::call( double x, double y, double z) const
{
   const double minimumR=1e-3*physicalconstants::R_E; //The dipole field is defined to be outside of Earth, and units are in meters     
   if(this->initialized==false)
      return 0.0;
   double r[3];
   
   r[0]= x-center[0];
   r[1]= y-center[1];
   r[2]= z-center[2];
   
   double r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
   
   if(r2<minimumR*minimumR)
      //  r2=minimumR*minimumR;
      return 0.0; //set zero field inside dipole
   
   const double r5 = (r2*r2*sqrt(r2));
   const double rdotq=q[0]*r[0] + q[1]*r[1] +q[2]*r[2];
   
   const double B=( 3*r[_fComponent]*rdotq-q[_fComponent]*r2)/r5;
   
   if(_derivative == 0) {
      //Value of B
      return B;
   }
   else if(_derivative == 1) {
      //first derivatives       
      unsigned int sameComponent;
      if(_dComponent==_fComponent)
         sameComponent=1;
      else
         sameComponent=0;
      
      return -5*B*r[_dComponent]/r2+
         (3*q[_dComponent]*r[_fComponent] -
          2*q[_fComponent]*r[_dComponent] +
          3*rdotq*sameComponent)/r5;
   }
   return 0; // dummy, but prevents gcc from yelling
}






