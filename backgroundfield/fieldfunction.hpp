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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*
Background magnetic field class of Vlasiator.
*/

#ifndef FIELDFUNCTION_HPP
#define FIELDFUNCTION_HPP
#include "functions.hpp"
#include <iostream>
#include <cstdlib>

class FieldFunction: public T3DFunction {
private:
protected:
   coordinate _fComponent;
   coordinate _dComponent;
   unsigned int _derivative;
public:
   FieldFunction(){
      //set sane initial values (x component of field)
      setComponent(X);
      setDerivComponent(X);
      setDerivative(0);
   }
   inline void setComponent(coordinate fComponent){ _fComponent=fComponent; }
   inline void setDerivComponent(coordinate dComponent){ _dComponent=dComponent; }
   inline void setDerivative(unsigned int derivative){
      _derivative=derivative;

      if( !(derivative==0 || derivative==1) ) {
         std::cerr<< "Only first derivatives supported!! "<< derivative<< std::endl;
         std::exit(1);
      } 
   }
};
#endif

