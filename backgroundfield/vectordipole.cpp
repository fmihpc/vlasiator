/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * Copyright 2017-2019 University of Helsinki
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
     return 0.0; //set zero field and derivatives outside "zero radius"

   /* This function is called from within other calls, one component at a time.
      The component in question is defined using the _fComponent index. If a derivative
      is requested, the direction of the derivative is defined using _dComponent. */

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
      
      /* Confirmed Battarbee 26.04.2019: This is the correct
	 3D dipole derivative.  */
      return -5*B*r[_dComponent]/r2+
         (3*q[_dComponent]*r[_fComponent] -
          2*q[_fComponent]*r[_dComponent] +
          3*rdotq*sameComponent)/r5;
   }

   /* Within transition range (between "full radius" and "zero radius"), use
      a vector potential scaled with the smootherstep function. Calculated
      and coded by Markus Battarbee, 30.04.2019 */

   // Calculate vector potential within transition range 
   double A[3];
   A[0] = (q[1]*r[2]-q[2]*r[1]) / (r2*r1); 
   A[1] = (q[2]*r[0]-q[0]*r[2]) / (r2*r1); 
   A[2] = (q[0]*r[1]-q[1]*r[0]) / (r2*r1); 
   // Coordinate within smootherstep function
   const double Sx = -(r1-radius[1])/(radius[1]-radius[0]);
   const double Sx2 = Sx*Sx;
   // Smootherstep and its radial derivative
   const double S2 = 6.*Sx2*Sx2*Sx - 15.*Sx2*Sx2 + 10.*Sx2*Sx;
   const double dS2dr = -(30.*Sx2*Sx2 - 60.*Sx2*Sx + 30.*Sx2)/(radius[1]-radius[0]);

   // Cartesian derivatives of S2
   double dS2cart[3];
   dS2cart[0] = (r[0]/r1)*dS2dr;
   dS2cart[1] = (r[1]/r1)*dS2dr;
   dS2cart[2] = (r[2]/r1)*dS2dr;      

   if(_derivative == 0) && (r1 > radius[0]) {
     /* Within transition range (between radius[0] and radius[1]) we
	multiply the magnetic field with the S2 smootherstep function
	and add an additional corrective term to remove divergence. This
	is based on using the dipole field vector potential and scaling
	it using the smootherstep function S2.

	Notation:
	 q = dipole moment (vector)
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
	B'(r) = del cross A'(r)
              =(NRL)= S2(Sx) del cross A(r) + del S2(Sx) cross A(r)
                    = S2(Sx) B(r)           + del S2(Sx) cross A(r)

     */
       double delS2crossA[3];
       delS2crossA[0] = dS2cart[1]*A[2] - dS2cart[2]*A[1];
       delS2crossA[1] = dS2cart[2]*A[0] - dS2cart[0]*A[2];
       delS2crossA[2] = dS2cart[0]*A[1] - dS2cart[1]*A[0];
       
       return S2*B + delS2crossA[_fComponent];
   }

   else if(_derivative == 1) && (r1 > radius[0]) {
       /* first derivatives of field calculated from diminishing vector potential

	  del B'(r) = S2(Sx) del B(r) + B(r) del S2(Sx) + del (del S2(Sx) cross A(r))
	  
	  component-wise:

	  del Bx = S2(Sx) del Bx + del S2(Sx) Bx + del(del S2(Sx) cross A)@i=x
	  del By = S2(Sx) del By + del S2(Sx) By + del(del S2(Sx) cross A)@i=y
	  del Bz = S2(Sx) del Bz + del S2(Sx) Bz + del(del S2(Sx) cross A)@i=z

	  where

	  del(del S2(Sx) cross A)@i=x = del (dS2/dy Az - dS/dz Ay)
	         = del(dS/dy) Az + dS/dy del Az - del(DS/dz) Ay - dS/dz del Ay

	  del(del S2(Sx) cross A)@i=y = del (dS2/dz Ax - dS/dx Az)
	         = del(dS/dz) Ax + dS/dz del Ax - del(DS/dx) Az - dS/dx del Az

	  del(del S2(Sx) cross A)@i=z = del (dS2/dx Ay - dS/dy Ax)
	         = del(dS/dx) Ay + dS/dx del Ay - del(DS/dy) Ax - dS/dy del Ax
       **********/

      unsigned int sameComponent;
      if(_dComponent==_fComponent)
         sameComponent=1;
      else
         sameComponent=0;

      // Regular derivative of B
      const double delB = -5*B*r[_dComponent]/r2+
         (3*q[_dComponent]*r[_fComponent] -
          2*q[_fComponent]*r[_dComponent] +
          3*rdotq*sameComponent)/r5;

       // Calculate del Ax, del Ay, del Az
       double delAx[3];
       double delAy[3];
       double delAz[3];
       delAx[0] = (-3./(r2*r2*r1))*(q[1]*r[2]-q[2]*r[1])*r[0];
       delAx[1] = (-3./(r2*r2*r1))*(q[1]*r[2]-q[2]*r[1])*r[1] -q[2]/(r2*r1);
       delAx[2] = (-3./(r2*r2*r1))*(q[1]*r[2]-q[2]*r[1])*r[2] +q[1]/(r2*r1);
       delAy[0] = (-3./(r2*r2*r1))*(q[2]*r[0]-q[0]*r[2])*r[0] +q[2]/(r2*r1);
       delAy[1] = (-3./(r2*r2*r1))*(q[2]*r[0]-q[0]*r[2])*r[1];
       delAy[2] = (-3./(r2*r2*r1))*(q[2]*r[0]-q[0]*r[2])*r[2] -q[0]/(r2*r1);
       delAz[0] = (-3./(r2*r2*r1))*(q[0]*r[1]-q[1]*r[0])*r[0] -q[1]/(r2*r1);
       delAz[1] = (-3./(r2*r2*r1))*(q[0]*r[1]-q[1]*r[0])*r[1] +q[0]/(r2*r1);
       delAz[2] = (-3./(r2*r2*r1))*(q[0]*r[1]-q[1]*r[0])*r[2];

       // Calculate del (dS2/dx), del (dS2/dy), del (dS2/dz)
       double deldS2dx[3];
       double deldS2dy[3];
       double deldS2dz[3];
       deldS2dx[0] = (-r[0]/(r2*r2*r1))*dS2dr*r[0] + dS2dr/r1;
       deldS2dx[1] = (-r[0]/(r2*r2*r1))*dS2dr*r[1];
       deldS2dx[2] = (-r[0]/(r2*r2*r1))*dS2dr*r[2];
       deldS2dy[0] = (-r[1]/(r2*r2*r1))*dS2dr*r[0];
       deldS2dy[1] = (-r[1]/(r2*r2*r1))*dS2dr*r[1] + dS2dr/r1;
       deldS2dy[2] = (-r[1]/(r2*r2*r1))*dS2dr*r[2];
       deldS2dz[0] = (-r[2]/(r2*r2*r1))*dS2dr*r[0];
       deldS2dz[1] = (-r[2]/(r2*r2*r1))*dS2dr*r[1];
       deldS2dz[2] = (-r[2]/(r2*r2*r1))*dS2dr*r[2] + dS2dr/r1;

       // Calculate del(del S2(Sx) cross A)@i=x, del(del S2(Sx) cross A)@i=y, del(del S2(Sx) cross A)@i=z
       double ddS2crossA[3][3];
       // derivatives of X-directional field
       ddS2crossA[0][0] = deldS2dy[0]*A[2] + dS2cart[1]*delAz[0] - deldS2dz[0]*A[1] - dS2cart[2]*delAy[0];
       ddS2crossA[0][1] = deldS2dy[1]*A[2] + dS2cart[1]*delAz[1] - deldS2dz[1]*A[1] - dS2cart[2]*delAy[1];
       ddS2crossA[0][2] = deldS2dy[2]*A[2] + dS2cart[1]*delAz[2] - deldS2dz[2]*A[1] - dS2cart[2]*delAy[2];
       // derivatives of Y-directional field
       ddS2crossA[1][0] = deldS2dz[0]*A[0] + dS2cart[2]*delAx[0] - deldS2dx[0]*A[2] - dS2cart[0]*delAz[0];
       ddS2crossA[1][1] = deldS2dz[1]*A[0] + dS2cart[2]*delAx[1] - deldS2dx[1]*A[2] - dS2cart[0]*delAz[1];
       ddS2crossA[1][2] = deldS2dz[2]*A[0] + dS2cart[2]*delAx[2] - deldS2dx[2]*A[2] - dS2cart[0]*delAz[2];
       // derivatives of Z-directional field
       ddS2crossA[2][0] = deldS2dx[0]*A[1] + dS2cart[0]*delAy[0] - deldS2dy[0]*A[0] - dS2cart[1]*delAx[0];
       ddS2crossA[2][1] = deldS2dx[1]*A[1] + dS2cart[0]*delAy[1] - deldS2dy[1]*A[0] - dS2cart[1]*delAx[1];
       ddS2crossA[2][2] = deldS2dx[2]*A[1] + dS2cart[0]*delAy[2] - deldS2dy[2]*A[0] - dS2cart[1]*delAx[2];

       return S2*delB + dS2cart[_dComponent]*B + ddS2crossA[_fComponent][_dComponent];
   }

   return 0; // dummy, but prevents gcc from yelling
}






