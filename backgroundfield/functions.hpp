/*
Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

enum coordinate { X, Y, Z };


class T1DFunction {public: virtual double call(double) const =0; virtual ~T1DFunction() {}};
class T2DFunction {public: virtual double call(double,double) const =0; virtual ~T2DFunction() {}};
class T3DFunction {public: virtual double call(double,double,double) const =0; virtual ~T3DFunction() {}};

// T2D_fix1, T2D_fix2: Fixing 1st or 2nd arg of a 2D function, thus making a 1D function

class T2D_fix1 : public T1DFunction {
private:
	const T2DFunction& f;
	double x;
public:
	T2D_fix1(const T2DFunction& f1, double x1) : f(f1),x(x1) {}
	virtual double call(double y) const {return f.call(x,y);}
	virtual ~T2D_fix1() {}
};

class T2D_fix2 : public T1DFunction {
private:
	const T2DFunction& f;
	double y;
public:
	T2D_fix2(const T2DFunction& f1, double y1) : f(f1),y(y1) {}
	virtual double call(double x) const {return f.call(x,y);}
	virtual ~T2D_fix2() {}
};

// T3D_fix1, T3D_fix2, T3D_fix3: Fixing 1st, 2nd or 3rd arg of a 3D function, thus making a 2D function

class T3D_fix1 : public T2DFunction {
private:
	const T3DFunction& f;
	double x;
public:
	T3D_fix1(const T3DFunction& f1, double x1) : f(f1),x(x1) {}
	virtual double call(double y, double z) const {return f.call(x,y,z);}
	virtual ~T3D_fix1() {}
};

class T3D_fix2 : public T2DFunction {
private:
	const T3DFunction& f;
	double y;
public:
	T3D_fix2(const T3DFunction& f1, double y1) : f(f1),y(y1) {}
	virtual double call(double x, double z) const {return f.call(x,y,z);}
	virtual ~T3D_fix2() {}
};

class T3D_fix3 : public T2DFunction {
private:
	const T3DFunction& f;
	double z;
public:
	T3D_fix3(const T3DFunction& f1, double z1) : f(f1),z(z1) {}
	virtual double call(double x, double y) const {return f.call(x,y,z);}
	virtual ~T3D_fix3() {}
};

// T3D_fix12, T3D_fix13, T3D_fix23: Fixing arguments 1,2, 1,3, or 2,3 of a 3D function, thus making a 1D function

class T3D_fix12 : public T1DFunction {
private:
	const T3DFunction& f;
	double x,y;
public:
	T3D_fix12(const T3DFunction& f1, double x1, double y1) : f(f1),x(x1),y(y1) {}
	virtual double call(double z) const {return f.call(x,y,z);}
	virtual ~T3D_fix12() {}
};

class T3D_fix13 : public T1DFunction {
private:
	const T3DFunction& f;
	double x,z;
public:
	T3D_fix13(const T3DFunction& f1, double x1, double z1) : f(f1),x(x1),z(z1) {}
	virtual double call(double y) const {return f.call(x,y,z);}
	virtual ~T3D_fix13() {}
};

class T3D_fix23 : public T1DFunction {
private:
	const T3DFunction& f;
	double y,z;
public:
	T3D_fix23(const T3DFunction& f1, double y1, double z1) : f(f1),y(y1),z(z1) {}
	virtual double call(double x) const {return f.call(x,y,z);}
	virtual ~T3D_fix23() {}
};

#endif
