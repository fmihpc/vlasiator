/*
Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef FIELDFUNCTION_HPP
#define FIELDFUNCTION_HPP
class FieldFunction {
public:
   virtual double value(unsigned int fComponent,double x,double y,double z) const =0;
   virtual double derivative(unsigned int fComponent,unsigned int dComponent,double x,double y,double z) const =0;
   virtual ~T3DFunction() {}
};

#endif
