/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef DIPOLE_HPP
#define DIPOLE_HPP
#include "fieldfunction.hpp"



class Dipole: public FieldFunction {
private:
   bool initialized;
   double q[3];                  // Dipole moment; set to (0,0,moment)
   double center[3];	// Coordinates where the dipole sits; set to (0,0,0)
public:
  
  Dipole(){
     this->initialized = false;
  }
   void initialize(const double moment,const double center_x, const double center_y, const double center_z, const double tilt_angle);
   virtual double call(double x, double y, double z) const;  
   virtual ~Dipole() {}
};

#endif

