/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef CONSTANTFIELD_HPP
#define CONSTANTFIELD_HPP
#include "fieldfunction.hpp"



class ConstantField: public FieldFunction {
private:
   bool _initialized;
   double _B[3]; // constant backgroundfield
public:
  
   ConstantField(){
      this->_initialized = false;
   }
   virtual ~ConstantField() {};

   
   void initialize(const double Bx,const double By, const double Bz);
   virtual double call(double x, double y, double z) const;
};

#endif

