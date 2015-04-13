/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011, 2012 Finnish Meteorological Institute
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






