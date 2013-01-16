/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef DIPOLE_HPP
#define DIPOLE_HPP
#include "functions.hpp"



class Dipole: public T3DFunction {
private:
  bool initialized;
  double q_x,q_y,q_z;                  // Dipole moment; set to (0,0,moment)
  double center_x, center_y, center_z;	// Coordinates where the dipole sits; set to (0,0,0)
  coordinate _fComponent;
public:
  
  Dipole(){
    this->initialized = false;
  }
  void initialize(const double moment);
  inline void setComponent(coordinate fComponent){ _fComponent=fComponent; } 
  
  virtual double call(double x, double y, double z) const;
  
   virtual ~Dipole() {}
};

#endif

