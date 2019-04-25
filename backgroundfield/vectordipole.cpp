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
#include "vectordipole.hpp"
#include "../common.h"

// tilt_angle_phi is from the z-axis in radians
// tilt_angle_theta is from the Sun-Earth-line in radians
void VectorDipole::initialize(const double moment,const double center_x, const double center_y, const double center_z, const double tilt_angle_phi=0, const double tilt_angle_theta=0, const double radius_f, const double radius_z){
   this->initialized = true;

   q[0]=-sin(tilt_angle_phi)*cos(tilt_angle_theta)*moment;
   q[1]=-sin(tilt_angle_phi)*sin(tilt_angle_theta)*moment;
   q[2]=-cos(tilt_angle_phi)*moment;

   center[0]=center_x;
   center[1]=center_y;
   center[2]=center_z;

   radius[0]=radius_f;
   radius[1]=radius_z;
}



double VectorDipole::call( double x, double y, double z) const
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

   if(r2>=radius[1]*radius[1])
     return 0.0; //set zero field and derivatives outside zero radius

   /* This function is called from within other calls, one component at a time.
      The component in question is defined using the _fComponent index. */

   const double r1 = sqrt(r2);
   const double r5 = (r2*r2*r1);
   const double rdotq=q[0]*r[0] + q[1]*r[1] +q[2]*r[2];   
   const double B=( 3*r[_fComponent]*rdotq-q[_fComponent]*r2)/r5;

   if(_derivative == 0) && (r1 <= radius[0]) 
     // Full dipole field within full radius
     return B;

   if(_derivative == 1)  && (r1 <= radius[0]){
      //first derivatives of full field
      unsigned int sameComponent;
      if(_dComponent==_fComponent)
         sameComponent=1;
      else
         sameComponent=0;
      
      // TODO: verify that this doesn't assume dipole aligned with z
      return -5*B*r[_dComponent]/r2+
         (3*q[_dComponent]*r[_fComponent] -
          2*q[_fComponent]*r[_dComponent] +
          3*rdotq*sameComponent)/r5;
   }

   // Calculate vector potential within transition range 
   double qcrossr[3];
   qcrossr[0] = q[1]*r[2]-q[2]*r[1]; 
   qcrossr[1] = q[2]*r[0]-q[0]*r[2]; 
   qcrossr[2] = q[0]*r[1]-q[1]*r[0]; 
   const double A = qcrossr / (r2*r1);
   // Coordinate within smootherstep function
   const double Sx = -(r1-radius[1])/(radius[1]-radius[0]);
   const double Sx2 = Sx*Sx;
   // Smootherstep and its radial derivative
   const double S2 = 6.*Sx2*Sx2*Sx - 15.*Sx2*Sx2 + 10.*Sx2*Sx;
   const double dS2dr = -(30.*Sx2*Sx2 - 60.*Sx2*Sx + 30.*Sx2)/(radius[1]-radius[0]);

   // Radial unit vector at that location in cartesian components
   double er[3];
   er[0]=r[0]/r1;
   er[1]=r[1]/r1;
   er[2]=r[2]/r1;

   // Cartesian derivatives of S2
   double dS2cart;
   dS2cart[0] = er[0]*dS2dr;
   dS2cart[1] = er[1]*dS2dr;
   dS2cart[2] = er[2]*dS2dr;      

   if(_derivative == 0) && (r1 > radius[0]) {
     /* Within transition range (between radius[0] and radius[1]) we
	multiply the magnetic field with the S2 smootherstep function
	and an additional corrective term to remove divergence. This
	is based on using the dipole field vector potential and scaling
	it using the smootherstep function S2.

	Notation:
	 m = dipole moment (vector)
	 r = position vector
	 R = position distance

	The regular dipole field vector potential
	A(r) = (mu0/4 pi R^3) * (q cross r)

	The smootherstep function
	         ( 0,                  x<=0
	S2(Sx) = ( 6x^5 -15x^4 +10x^3, 0<=x<=1
                 ( 1,                  x>=1

	Radial distance scaling for S2
        Sx = -(R-radius[1])/(radius[1]-radius[0])

	The scaled vector potential is A'(r) = A(r)*S2(Sx)

	The scaled magnetic field is
	del cross A'(r)
            = S2(Sx) del cross A(r) + del S2(Sx) cross A(r)
            = S2(sx) B(r)           + del S2(Sx) cross A(r)

     */
     double correctionterm[3];
     correctionterm[0] = dS2cart[1]*A[2] - dS2cart[2]*A[1];
     correctionterm[1] = dS2cart[2]*A[0] - dS2cart[0]*A[2];
     correctionterm[2] = dS2cart[0]*A[1] - dS2cart[1]*A[0];

     return B*S2 + correctionterm[_fComponent];
   }

   else if(_derivative == 1) && (r1 > radius[0]) {
      // first derivatives of field calculated from diminishing vector potential

      // TODO: calculate derivatives and implement
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






