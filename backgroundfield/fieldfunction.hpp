/*
Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef FIELDFUNCTION_HPP
#define FIELDFUNCTION_HPP

enum component { X, Y, Z };

class FieldFunction {
public:
  inline bool setFComponent(component fComponent){ _fComponent=fComponent; }
  inline bool setDComponent(component dComponent){ _dComponent=dComponent; }
  
  virtual double value(double x,double y,double z) const =0;
  virtual double derivative(double x,double y,double z) const =0;
  virtual ~FieldFunction() {}
  
protected:
  component _fComponent;
  component _dComponent;
};

#endif
