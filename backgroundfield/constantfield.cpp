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






