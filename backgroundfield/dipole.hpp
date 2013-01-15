/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef DIPOLE_HPP
#define DIPOLE_HPP

class Dipole: public T3DFunction {
private:
   bool initialized;
   double q[3];
   double DipoleBxTerms(const double x[3], double r, double invr, double invr2) const;
   double DipoleByTerms(const double x[3], double r, double invr, double invr2) const;
   double DipoleBzTerms(const double x[3], double r, double invr, double invr2) const;

   double dipmom;
   double center_x, center_y, center_z;	// coordinates where the dipole (quadrupole...) sits; default (0,0,0)

public:
   
   Dipole(const double moment){
      this->initialized = false;
   }
   void initialize(const double moment);
   
 
   virtual double value(unsigned int fComponent,double x, double y, double z) const;
   virtual double derivative(unsigned int fComponent,unsigned int dComponent,double x, double y, double z) const;

   virtual ~Dipole() {}
};

#endif

