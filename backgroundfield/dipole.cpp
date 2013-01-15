/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011, 2012 Finnish Meteorological Institute
*/

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "../common.h"
#include "B0.hpp"

void Dipole::initialize(const double moment)
{
   this->initialized = true;
   q[0]=0.0;
   q[1]=0.0;
   q[2]=moment;
   this->center_x = 0;
   this->center_y = 0;
   this->center_z = 0;
}





inline double Dipole::DipoleBxTerms(const double x[3], double r, double invr, double invr2) const
{
   const double invr5 = invr*sqr(invr2);
   const double r2 = sqr(r);
   //recflops(17);
   return (-(r2*q[0]) + 3*q[0]*sqr(x[0]) + 3*q[1]*x[0]*x[1]
           + 3*q[2]*x[0]*x[2])*invr5;
}

inline double Dipole::DipoleByTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr5 = invr*sqr(invr2);
	const double r2 = sqr(r);
	//recflops(17);
	return (-(r2*q[1]) + 3*q[1]*sqr(x[1]) + 3*q[2]*x[1]*x[2]
			+ 3*q[0]*x[1]*x[0])*invr5;
}

inline double Dipole::DipoleBzTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr5 = invr*sqr(invr2);
	const double r2 = sqr(r);
	//recflops(17);
	return (-(r2*q[2]) + 3*q[2]*sqr(x[2]) + 3*q[0]*x[2]*x[0]
			+ 3*q[1]*x[2]*x[1])*invr5;
}

double Dipole::value(unsigned int fComponent, double x, double y, double z) const
{
   if(this->initialized==false)
      return 0.0;
   
   x-= center_x;
   y-= center_y;
   z-= center_z;
   const double r2 = max(sqr(minimum_r), sqr(x) + sqr(y) + sqr(z));
   const double r = sqrt(r2);
   const double invr2 = 1.0/r2;
   const double invr = 1.0/r;
   const double R[3] = {x,y,z};
   
   switch (fComponent) {
       case 0: return DipoleBxTerms(R,r,invr,invr2);
       case 1: return DipoleByTerms(R,r,invr,invr2);
       case 2: return DipoleBzTerms(R,r,invr,invr2);
          
   }
   return 0;	// dummy, but prevents gcc from yelling
}



double Dipole::derivative(unsigned int fComponent, unsigned int dComponent,  double x, double y, double z) const
{
   if(this->initialized==false)
      return 0.0;
   
   x-= center_x;
   y-= center_y;
   z-= center_z;
   const double r2 = max(sqr(minimum_r), sqr(x) + sqr(y) + sqr(z));
   const double r = sqrt(r2);
   const double invr2 = 1.0/r2;
   const double invr = 1.0/r;
   const double R[3] = {x,y,z};
   return 0;	// dummy, but prevents gcc from yelling
}








