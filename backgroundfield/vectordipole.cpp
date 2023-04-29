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
void VectorDipole::initialize(const double moment,const double center_x, const double center_y, const double center_z, const double tilt_angle_phi, const double tilt_angle_theta, const double xlimit_f, const double xlimit_z, const double IMF_Bx, const double IMF_By, const double IMF_Bz){
   this->initialized = true;

   q[0]=-sin(tilt_angle_phi)*cos(tilt_angle_theta)*moment;
   q[1]=-sin(tilt_angle_phi)*sin(tilt_angle_theta)*moment;
   q[2]=-cos(tilt_angle_phi)*moment;

   center[0]=center_x;
   center[1]=center_y;
   center[2]=center_z;

   // Scale dipole as a function of x-coordinate
   xlimit[0]=xlimit_f; // Full dipole when x < xlimit_f
   xlimit[1]=xlimit_z; // Zero dipole when x > xlimit_z

   // Going from xlimit_f to xlimit_z, scale in IMF B-field
   IMF[0]=IMF_Bx;
   IMF[1]=IMF_By;
   IMF[2]=IMF_Bz;

   // TODO: If values for xlimit are zero, instead place them as 15 RE and Xmax-2*cellsize?
}



double VectorDipole::operator()( double x, double y, double z, coordinate component, unsigned int derivative, coordinate dcomponent) const {
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
   
   if(r[0]>=xlimit[1]){
      //set zero or IMF field and derivatives outside "zero x limit"
      if(derivative == 0) {
         return IMF[component]; 
      } else {
         return 0.0;
      }
   }
   /* This function is called from within other calls, one component at a time.
      The component in question is defined using the component index. If a derivative
      is requested, the direction of the derivative is defined using dcomponent. */

   const double r1 = sqrt(r2);
   const double r5 = (r2*r2*r1);
   const double rdotq=q[0]*r[0] + q[1]*r[1] +q[2]*r[2];   
   const double B=( 3*r[component]*rdotq-q[component]*r2)/r5;

   if((derivative == 0) && (r[0] <= xlimit[0])) {
      // Full dipole field within full xlimit
      return B;
   }

   if((derivative == 1) && (r[0] <= xlimit[0])){
      //first derivatives of full field
      unsigned int sameComponent;
      if(dcomponent==component) {
         sameComponent=1;
      } else {
         sameComponent=0;
      }
      
      /* Confirmed Battarbee 26.04.2019: This is the correct
	 3D dipole derivative.  */
      return -5*B*r[dcomponent]/r2+
         (3*q[dcomponent]*r[component] -
          2*q[component]*r[dcomponent] +
          3*rdotq*sameComponent)/r5;
   }

   /* Within transition range (between "full x limit" and "zero x limit"), use
      a vector potential scaled with the smootherstep function. Calculated
      and coded by Markus Battarbee, 08.05.2019 */

   // Calculate vector potential within transition range 
   double A[3];
   A[0] = (q[1]*r[2]-q[2]*r[1]) / (r2*r1); 
   A[1] = (q[2]*r[0]-q[0]*r[2]) / (r2*r1); 
   A[2] = (q[0]*r[1]-q[1]*r[0]) / (r2*r1); 
   // Calculate vector potential for IMF scaling
   double IMFA[3];
   IMFA[0] = 0.5*(IMF[1]*r[2] - IMF[2]*r[1]);
   IMFA[1] = 0.5*(IMF[2]*r[0] - IMF[0]*r[2]);
   IMFA[2] = 0.5*(IMF[0]*r[1] - IMF[1]*r[0]);
   const double IMFB = IMF[component];
   
   // Coordinate within smootherstep function (x-coordinate only)
   const double s = -(r[0]-xlimit[1])/(xlimit[1]-xlimit[0]);
   const double ss = s*s;
   // Smootherstep and its x-directional derivative
   const double S2 = 6.*ss*ss*s - 15.*ss*ss + 10.*ss*s;
   const double dS2dx = -(30.*ss*ss - 60.*ss*s + 30.*ss)/(xlimit[1]-xlimit[0]);

   // Smootherstep for IMF
   const double IMFs = (r[0]-xlimit[0])/(xlimit[1]-xlimit[0]);
   const double IMFss = IMFs*IMFs;
   // Smootherstep and its x-directional derivative
   const double IMFS2 = 6.*IMFss*IMFss*IMFs - 15.*IMFss*IMFss + 10.*IMFss*IMFs;
   const double IMFdS2dx = (30.*IMFss*IMFss - 60.*IMFss*IMFs + 30.*IMFss)/(xlimit[1]-xlimit[0]);

   // Cartesian derivatives of S2
   double dS2cart[3];
   dS2cart[0] = dS2dx; //(r[0]/r1)*dS2dr;
   dS2cart[1] = 0;     //(r[1]/r1)*dS2dr;
   dS2cart[2] = 0;     //(r[2]/r1)*dS2dr;      

   // Cartesian derivatives of S2
   double IMFdS2cart[3];
   IMFdS2cart[0] = IMFdS2dx; //(r[0]/r1)*dS2dr;
   IMFdS2cart[1] = 0;     //(r[1]/r1)*dS2dr;
   IMFdS2cart[2] = 0;     //(r[2]/r1)*dS2dr;      

   if(derivative == 0) {
     /* Within transition range (between xlimit[0] and xlimit[1]) we
        multiply the magnetic field with the S2 smootherstep function
        and add an additional corrective term to remove divergence. This
        is based on using the dipole field vector potential and scaling
        it using the smootherstep function S2.

        Notation:
        q = dipole moment (vector)
        r = position vector
        x = x-coordinate r[0]
        R = position distance

        The regular dipole field vector potential
        A(r) = (mu0/4 pi R^3) * (q cross r)

        The smootherstep function
                ( 0,                  s<=0
        S2(s) = ( 6s^5 -15s^4 +10s^3, 0<=s<=1
                ( 1,                  s>=1

        Radial distance scaling for S2
        s = -(x-xlimit[1])/(xlimit[1]-xlimit[0])
        ds = -dx/(xlimit[1]-xlimit[0])

        The scaled vector potential is A'(r) = A(r)*S2(s)

        The scaled magnetic field is
        B'(r) = del cross A'(r)
              =(NRL)= S2(s) del cross A(r) + del S2(s) cross A(r)
                    = S2(s) B(r)           + del S2(s) cross A(r)

     */
       double delS2crossA[3];
       //delS2crossA[0] = dS2cart[1]*A[2] - dS2cart[2]*A[1]; 
       //delS2crossA[1] = dS2cart[2]*A[0] - dS2cart[0]*A[2]; 
       //delS2crossA[2] = dS2cart[0]*A[1] - dS2cart[1]*A[0]; 
       // Don't calculate zero terms
       delS2crossA[0] = 0;
       delS2crossA[1] = -dS2cart[0]*A[2];
       delS2crossA[2] = dS2cart[0]*A[1];

       double IMFdelS2crossA[3];
       // Don't calculate zero terms
       IMFdelS2crossA[0] = 0;
       IMFdelS2crossA[1] = -IMFdS2cart[0]*IMFA[2];
       IMFdelS2crossA[2] = IMFdS2cart[0]*IMFA[1];

       //return S2*B + delS2crossA[component];
       return S2*B + delS2crossA[component] + IMFS2*IMFB + IMFdelS2crossA[component];
   }

   else if(derivative == 1) {
       /* first derivatives of field calculated from diminishing vector potential

          del B'(r) = S2(s) del B(r) + B(r) del S2(s) + del (del S2(s) cross A(r))

          component-wise:

          del Bx = S2(s) del Bx + del S2(s) Bx + del(del S2(s) cross A)@i=x
          del By = S2(s) del By + del S2(s) By + del(del S2(s) cross A)@i=y
          del Bz = S2(s) del Bz + del S2(s) Bz + del(del S2(s) cross A)@i=z

          where

          del(del S2(s) cross A)@i=x = del (dS2/dy Az - dS2/dz Ay)
             = del(dS2/dy) Az + dS2/dy del Az - del(DS/dz) Ay - dS2/dz del Ay

          del(del S2(s) cross A)@i=y = del (dS2/dz Ax - dS2/dx Az)
             = del(dS2/dz) Ax + dS2/dz del Ax - del(DS/dx) Az - dS2/dx del Az

          del(del S2(s) cross A)@i=z = del (dS2/dx Ay - dS2/dy Ax)
             = del(dS2/dx) Ay + dS2/dx del Ay - del(DS/dy) Ax - dS2/dy del Ax

          note that dS2/dy == dS2/dz == 0
        **********/

      unsigned int sameComponent;
      if(dcomponent==component) {
         sameComponent=1;
      } else {
         sameComponent=0;
      }

      // Regular derivative of B
      const double delB = -5*B*r[dcomponent]/r2+
         (3*q[dcomponent]*r[component] -
          2*q[component]*r[dcomponent] +
          3*rdotq*sameComponent)/r5;

      // IMF field is constant
      const double IMFdelB = 0.;

       // Calculate del Ax, del Ay, del Az
      double delAy[3];
      double delAz[3];
      delAy[0] = (-3./(r2*r2*r1))*(q[2]*r[0]-q[0]*r[2])*r[0] +q[2]/(r2*r1);
      delAy[1] = (-3./(r2*r2*r1))*(q[2]*r[0]-q[0]*r[2])*r[1];
      delAy[2] = (-3./(r2*r2*r1))*(q[2]*r[0]-q[0]*r[2])*r[2] -q[0]/(r2*r1);
      delAz[0] = (-3./(r2*r2*r1))*(q[0]*r[1]-q[1]*r[0])*r[0] -q[1]/(r2*r1);
      delAz[1] = (-3./(r2*r2*r1))*(q[0]*r[1]-q[1]*r[0])*r[1] +q[0]/(r2*r1);
      delAz[2] = (-3./(r2*r2*r1))*(q[0]*r[1]-q[1]*r[0])*r[2];
      //  derivatives of x-directional component of A are not needed here
      //double delAx[3];
      //delAx[0] = (-3./(r2*r2*r1))*(q[1]*r[2]-q[2]*r[1])*r[0];
      //delAx[1] = (-3./(r2*r2*r1))*(q[1]*r[2]-q[2]*r[1])*r[1] -q[2]/(r2*r1);
      //delAx[2] = (-3./(r2*r2*r1))*(q[1]*r[2]-q[2]*r[1])*r[2] +q[1]/(r2*r1);

      // Calculate del IMFAx, del IMFAy, del IMFAz
      //double IMFdelAx[3]; // Unused, since always zero
      //IMFdelAx[0] = 0.;
      //IMFdelAx[1] = -0.5*IMF[2];
      //IMFdelAx[2] =  0.5*IMF[1];
      double IMFdelAy[3];
      IMFdelAy[0] =  0.5*IMF[2];
      IMFdelAy[1] = 0.0;
      IMFdelAy[2] = -0.5*IMF[0];
      double IMFdelAz[3];
      IMFdelAz[0] = -0.5*IMF[1];
      IMFdelAz[1] =  0.5*IMF[0];
      IMFdelAz[2] = 0.0;

      // Calculate del (dS2/dx), del (dS2/dy), del (dS2/dz)
      // Of course now only del (dS2/dx) is non-zero
      const double ddidS2dx = 60.*(2.*ss*s - 3.*ss + s)/((xlimit[1]-xlimit[0])*(xlimit[1]-xlimit[0]));
      // This is the same for IMF field scaling as well

      double deldS2dx[3];
      //double deldS2dy[3];
      //double deldS2dz[3];
      deldS2dx[0] = ddidS2dx;
      deldS2dx[1] = 0;
      deldS2dx[2] = 0;
      /*
      ddidS2dr = 60.*(2.*ss*s - 3.*ss + s)/(r2*(xlimit[1]-xlimit[0])*(xlimit[1]-xlimit[0]));
      deldS2dx[0] = ddidS2dr*r[0]*r[0] -(r[0]/(r2*r1))*dS2dr*r[0] + dS2dr/r1;
      deldS2dx[1] = ddidS2dr*r[0]*r[1] -(r[0]/(r2*r1))*dS2dr*r[1];
      deldS2dx[2] = ddidS2dr*r[0]*r[2] -(r[0]/(r2*r1))*dS2dr*r[2];
      deldS2dy[0] = ddidS2dr*r[1]*r[0] -(r[1]/(r2*r1))*dS2dr*r[0];
      deldS2dy[1] = ddidS2dr*r[1]*r[1] -(r[1]/(r2*r1))*dS2dr*r[1] + dS2dr/r1;
      deldS2dy[2] = ddidS2dr*r[1]*r[2] -(r[1]/(r2*r1))*dS2dr*r[2];
      deldS2dz[0] = ddidS2dr*r[2]*r[0] -(r[2]/(r2*r1))*dS2dr*r[0];
      deldS2dz[1] = ddidS2dr*r[2]*r[1] -(r[2]/(r2*r1))*dS2dr*r[1];
      deldS2dz[2] = ddidS2dr*r[2]*r[2] -(r[2]/(r2*r1))*dS2dr*r[2] + dS2dr/r1;

      // Calculate del(del S2(s) cross A)@i=x, del(del S2(s) cross A)@i=y, del(del S2(s) cross A)@i=z
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
      */

      // Only include components which are nonzero
      double ddS2crossA[3][3];
      // derivatives of X-directional field
      ddS2crossA[0][0] = 0;
      ddS2crossA[0][1] = 0;
      ddS2crossA[0][2] = 0;
      // derivatives of Y-directional field
      ddS2crossA[1][0] = - deldS2dx[0]*A[2] - dS2cart[0]*delAz[0];
      ddS2crossA[1][1] = - deldS2dx[1]*A[2] - dS2cart[0]*delAz[1];
      ddS2crossA[1][2] = - deldS2dx[2]*A[2] - dS2cart[0]*delAz[2];
      // derivatives of Z-directional field
      ddS2crossA[2][0] = deldS2dx[0]*A[1] + dS2cart[0]*delAy[0];
      ddS2crossA[2][1] = deldS2dx[1]*A[1] + dS2cart[0]*delAy[1];
      ddS2crossA[2][2] = deldS2dx[2]*A[1] + dS2cart[0]*delAy[2];

      // Now for IMF portion
      // Only include components which are nonzero
      double IMFddS2crossA[3][3];
      // derivatives of X-directional field
      IMFddS2crossA[0][0] = 0;
      IMFddS2crossA[0][1] = 0;
      IMFddS2crossA[0][2] = 0;
      // derivatives of Y-directional field
      IMFddS2crossA[1][0] = - deldS2dx[0]*IMFA[2] - IMFdS2cart[0]*IMFdelAz[0];
      IMFddS2crossA[1][1] = - deldS2dx[1]*IMFA[2] - IMFdS2cart[0]*IMFdelAz[1];
      IMFddS2crossA[1][2] = - deldS2dx[2]*IMFA[2] - IMFdS2cart[0]*IMFdelAz[2];
      // derivatives of Z-directional field
      IMFddS2crossA[2][0] = deldS2dx[0]*IMFA[1] + IMFdS2cart[0]*IMFdelAy[0];
      IMFddS2crossA[2][1] = deldS2dx[1]*IMFA[1] + IMFdS2cart[0]*IMFdelAy[1];
      IMFddS2crossA[2][2] = deldS2dx[2]*IMFA[1] + IMFdS2cart[0]*IMFdelAy[2];

      //return S2*delB + dS2cart[dcomponent]*B + ddS2crossA[component][dcomponent];
      return S2*delB + dS2cart[dcomponent]*B + ddS2crossA[component][dcomponent] + 
         IMFS2*IMFdelB + IMFdS2cart[dcomponent]*IMFB + IMFddS2crossA[component][dcomponent];
   }

   return 0; // dummy, but prevents gcc from yelling
}






