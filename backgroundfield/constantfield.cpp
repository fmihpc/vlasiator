/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011, 2012 Finnish Meteorological Institute
*/

#include <stdlib.h>
#include <math.h>
#include "constantfield.hpp"
#include "../common.h"


void ConstantField::initialize(const double Bx,const double By, const double Bz){
   _B[0]=Bx;
   _B[1]=By;
   _B[2]=Bz;
   _initialized=true;
}



double ConstantField::call( double , double , double ) const
{
   if(_derivative == 0) {
      //Value of B
      return _B[_fComponent];
   }
   else if(_derivative > 0) {
      //all derivatives are zero
      return 0.0;
   }
   return 0; // dummy, but prevents gcc from yelling
}






