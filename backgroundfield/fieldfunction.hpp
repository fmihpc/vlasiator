/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef FIELDFUNCTION_HPP
#define FIELDFUNCTION_HPP
#include "functions.hpp"



class FieldFunction: public T3DFunction {
private:
protected:
  coordinate _fComponent;
public:
   inline void setComponent(coordinate fComponent){ _fComponent=fComponent; } 
};

#endif

