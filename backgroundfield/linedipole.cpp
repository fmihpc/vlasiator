/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011, 2012 Finnish Meteorological Institute
*/

#include <stdlib.h>
#include <math.h>
#include "linedipole.hpp"
#include "../common.h"

//tilt_angle is agains the z axis in the x-z plane. In radians*/
void LineDipole::initialize(const double moment,const double tilt_angle=0)
{
   
   this->initialized = true;
//    q[0]=-sin(tilt_angle)*moment;
//    q[1]=0.0;
//    q[2]=-cos(tilt_angle)*moment;
   q[2]=moment;
   center[0]=0.0;
   center[1]=0.0;
   center[2]=0.0;
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
   
   double r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
   
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






