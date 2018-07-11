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
#include "linedipole.hpp"
#include "../common.h"


void LineDipole::initialize(const double moment,const double center_x, const double center_y, const double center_z){
   this->initialized = true;
   q[0]=0.0;
   q[1]=0.0;
   q[2]=moment;
   center[0]=center_x;
   center[1]=center_y;
   center[2]=center_z;
}



double LineDipole::call( double x, double y, double z) const
{
   const double minimumR=1e-3*physicalconstants::R_E; //The dipole field is defined to be outside of Earth, and units are in meters     
   if(this->initialized==false)
      return 0.0;
   double r[3];
   
   r[0]= x-center[0];
   r[1]= y-center[1];
   r[2]= z-center[2];
   
   double r2 = r[0]*r[0]+r[2]*r[2]; // r[1] not necessary in this case, removed to enable proper cylindrical ionosphere (ionosphere.geometry = 3)
   
   if(r2<minimumR*minimumR)
      //  r2=minimumR*minimumR;
      return 0.0; //set zero field inside dipole
   
   const double r6 = (r2*r2*r2);
//    const double rdotq=q[0]*r[0] + q[1]*r[1] +q[2]*r[2];
   const double D = -q[2]; 
   
   const double DerivativeSameComponent=D*( 2*r[2]*(r[2]*r[2]-3*r[0]*r[0]))/r6;
   const double DerivativeDiffComponent=D*( 2*r[0]*(r[0]*r[0]-3*r[2]*r[2]))/r6;
   //const double B;
   //const double der;
   
   if(_derivative == 0) {
      if(_fComponent == 0)
         return D*2*r[0]*r[2]/(r2*r2);
      if(_fComponent == 2)
         return D*(r[2]*r[2]-r[0]*r[0])/(r2*r2); 
      if(_fComponent == 1)
         return 0;
   }
   else if(_derivative == 1) {
      //first derivatives
      if(_dComponent== 1 || _fComponent==1) {
         return 0;
      }
      else if(_dComponent==_fComponent) {
         if(_fComponent == 0) {
            return DerivativeSameComponent;
         }
         else if(_fComponent == 2) {
            return -DerivativeSameComponent;
         }
      }
      else { 
         return DerivativeDiffComponent;
      }
 
   }
   return 0;   // dummy, but prevents gcc from yelling
}






