/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef FIELDFUNCTION_HPP
#define FIELDFUNCTION_HPP
#include "functions.hpp"
#include <iostream>


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

