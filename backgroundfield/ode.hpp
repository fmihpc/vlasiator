/*
Ordinary differential equation stuff of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2011 Finnish Meteorological Institute
*/

#ifndef ODE_HPP
#define ODE_HPP

/*!
Type Ordinary Differential Equation Right-Hand-Side
*/
class TODERHS {
public:
	virtual void derivs(double x, const double y[], double yout[]) = 0;
	virtual ~TODERHS() {}
};

extern
void BulirschStoerStep(
	double y[],	// nv-vector, the dynamic variable vector on input and output
	const double dydx[],	// must be initialized to derivs(x,y,dydx) on input
	int nv,	// length of y,dydx,yscal vectors
	double &x,	// the argument variable (input/output)
	double htry,	// the first step size to try
	double eps,	// relative tolerance
	const double yscal[],	// 'typical value' vector for dynamic variables
	double& hdid,	// output: the step size it doublely did. If hdid < htry, there was some problem.
	double& hnext,	// output: suggestion for next step size
	TODERHS& rhs);	// object that computes the right-hand sides of the equation (member function derivs)

#endif
