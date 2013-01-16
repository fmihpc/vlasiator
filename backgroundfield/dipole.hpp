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
   double q_x,q_y,q_z;                  // Dipole moment; set to (0,0,moment)
   double center_x, center_y, center_z;	// Coordinates where the dipole sits; set to (0,0,0)

   double bx(double x,double y,double z, double r2, double invr5) const;
   double by(double x,double y,double z, double r2, double invr5) const;
   double bz(double x,double y,double z, double r2, double invr5) const;

      
public:
  
  Dipole(){
     this->initialized = false;
  }

   void initialize(const double moment);
  
   virtual double call(double x, double y, double z) const;
  
   virtual ~Dipole() {}
};

#endif

