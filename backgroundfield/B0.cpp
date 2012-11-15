/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011, 2012 Finnish Meteorological Institute
*/

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "../common.h"
#include "B0.hpp"

static int PowerOfThree(int n)		// return 3^n
{
	int i, result = 1;
	for (i=0; i<n; i++) result*= 3;
	return result;
}

static int FlatIndex(const int ind[], int n)
{
	int n1, result = 0;
	for (n1=0; n1<n; n1++) {
		result = 3*result + ind[n1];
	}
	return result;
}

void TB0::initialize(int maxorder1)
{
	using namespace physicalconstants;

	this->initialized = true;
	this->maxorder = maxorder1;
	if (this->maxorder < 0 || this->maxorder > MAX_MULTIPOLE_ORDER) {
		cerr << "*** TB0: bad maxorder=" << maxorder1 << "\n";
		exit(1);
	}
	int n;
	for (n=0; n<=MAX_MULTIPOLE_ORDER; n++) q[n] = 0;
	// Allocate q memory pools
	for (n = 0; n <= this->maxorder; n++) {
		const int N = PowerOfThree(n);
		q[n] = new double [N];
		memset(q[n],0,N*sizeof(double));
	}

	this->setcoeff(0, 0);
	this->setcoeff(1, 0);
	this->setcoeff(2, -this->dipmom);
	/*minimum_r = 1e-10;
	Baccuracy = 1e-10;*/
	this->minimum_r = 1e-3 * R_E;
	this->Baccuracy = 1e-17;
	this->center_x = this->center_y = this->center_z = 0;
}

void TB0::dealloc()
{
	if (!this->initialized) {
		return;
	}

	int n;
	for (n=0; n<=MAX_MULTIPOLE_ORDER; n++) {
		delete [] q[n];
	}
}

void TB0::setcoeff(const int ind[], int n1, double q1)	// general case
{
	if (n1 < 0 || n1 > this->maxorder) {cerr << "*** TB0::setcoeff bad n=" << n1 << "\n"; exit(1);}
	q[n1][FlatIndex(ind,n1)] = q1;
	if (q1 != 0) all_qs_are_zero = 0;
}

void TB0::setcoeff(double q1) {q[0][0] = q1;if (q1 != 0) all_qs_are_zero = 0;}

void TB0::setcoeff(int i, double q1) {q[1][i] = q1; if (q1 != 0) all_qs_are_zero = 0;}

void TB0::setcoeff(int i, int j, double q1) {q[2][i*3+j] = q1; if (q1 != 0) all_qs_are_zero = 0;}

double TB0::EvaluateMultipoleTerm(int n, const double q1[], const double r[3]) const
/*
   n=0:   q1[0]
   n=1:   q1.r
   n=2:   sum_{ij} q1_{ij} x_i x_j
   n=3:   sum_{ijk} q1_{ijk} x_i x_j x_k
     ...
*/
{
	if (n == 0) {
		// monopole term
		return q1[0];
	} else if (n == 1) {
		// dipole term
		//recflops(5);
		return q1[0]*r[0] + q1[1]*r[1] + q1[2]*r[2];
	} else {
		// higher term, use recursion
		//recflops(8);
		const int N = PowerOfThree(n-1);
		return EvaluateMultipoleTerm(n-1,q1,r)*r[0]
			+  EvaluateMultipoleTerm(n-1,q1+N,r)*r[1]
			+  EvaluateMultipoleTerm(n-1,q1+2*N,r)*r[2];
	}
}

double TB0::EvaluatePhi(const double r[3]) const
/*
          q    q.r            xi xj              xi xj xk
   phi = --- + ---- + sum Q   -----  +  sum Q    -------- + ...
          r    r^3     ij  ij  r^5      ijk  ijk   r^7
*/
{
	const double x = r[0] - center_x;
	const double y = r[1] - center_y;
	const double z = r[2] - center_z;
	
	const double invr2 = 1.0/(max(sqr(minimum_r), sqr(x) + sqr(y) + sqr(z)));
	const double invr = sqrt(invr2);
	const double rshifted[3] = {x,y,z};
	double phi;
	if (this->maxorder == 1 && q[0][0] == 0.) {
		// Fast branch for dipole potential only
		phi = (invr*invr2)*(q[1][0]*rshifted[0] + q[1][1]*rshifted[1] + q[1][2]*rshifted[2]);
	} else {
		// General branch
		int n;
		phi = 0;
		double Rcoeff = invr;
		for (n = 0; n <= this->maxorder; n++) {
			const double mpoleterm = EvaluateMultipoleTerm(n,q[n],rshifted);
			phi+= Rcoeff*mpoleterm;
			Rcoeff*= invr2;
		}
	}
	return phi;
}

double TB0::MonopoleBxTerms(const double x[3], double /*r*/, double invr, double invr2) const
{
	const double invr3 = invr*invr2;
	//recflops(3);
	return (invr3*q[0][0])*x[0];
}

double TB0::MonopoleByTerms(const double x[3], double /*r*/, double invr, double invr2) const
{
	const double invr3 = invr*invr2;
	//recflops(3);
	return (invr3*q[0][0])*x[1];
}

double TB0::MonopoleBzTerms(const double x[3], double /*r*/, double invr, double invr2) const
{
	const double invr3 = invr*invr2;
	//recflops(3);
	return (invr3*q[0][0])*x[2];
}

#define qref(i) q[1][i]

inline double TB0::DipoleBxTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr5 = invr*sqr(invr2);
	const double r2 = sqr(r);
	//recflops(17);
	return (-(r2*qref(0)) + 3*qref(0)*sqr(x[0]) + 3*qref(1)*x[0]*x[1]
			+ 3*qref(2)*x[0]*x[2])*invr5;
}

inline double TB0::DipoleByTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr5 = invr*sqr(invr2);
	const double r2 = sqr(r);
	//recflops(17);
	return (-(r2*qref(1)) + 3*qref(1)*sqr(x[1]) + 3*qref(2)*x[1]*x[2]
			+ 3*qref(0)*x[1]*x[0])*invr5;
}

inline double TB0::DipoleBzTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr5 = invr*sqr(invr2);
	const double r2 = sqr(r);
	//recflops(17);
	return (-(r2*qref(2)) + 3*qref(2)*sqr(x[2]) + 3*qref(0)*x[2]*x[0]
			+ 3*qref(1)*x[2]*x[1])*invr5;
}

#undef qref

#define qref(i,j) q[2][(i-1)*3+j-1]

double TB0::QuadrupoleBxTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr7 = invr*invr2*sqr(invr2);
	const double r2 = sqr(r);
	//recflops(63);
	return (-2*r2*qref(1, 1)*x[0] + 5*qref(1, 1)*pow3(x[0]) - 
			r2*qref(1, 2)*x[1] - r2*qref(2, 1)*x[1] + 
			5*qref(1, 2)*sqr(x[0])*x[1] + 5*qref(2, 1)*sqr(x[0])*x[1] + 
			5*qref(2, 2)*x[0]*sqr(x[1]) - r2*qref(1, 3)*x[2] - 
			r2*qref(3, 1)*x[2] + 5*qref(1, 3)*sqr(x[0])*x[2] + 
			5*qref(3, 1)*sqr(x[0])*x[2] + 5*qref(2, 3)*x[0]*x[1]*x[2] + 
			5*qref(3, 2)*x[0]*x[1]*x[2] + 5*qref(3, 3)*x[0]*sqr(x[2]))*invr7;
}

double TB0::QuadrupoleByTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr7 = invr*invr2*sqr(invr2);
	const double r2 = sqr(r);
	//recflops(63);
	return (-(r2*qref(1, 2)*x[0]) - r2*qref(2, 1)*x[0] - 
			2*r2*qref(2, 2)*x[1] + 5*qref(1, 1)*sqr(x[0])*x[1] + 
			5*qref(1, 2)*x[0]*sqr(x[1]) + 5*qref(2, 1)*x[0]*sqr(x[1]) + 
			5*qref(2, 2)*pow3(x[1]) - r2*qref(2, 3)*x[2] - r2*qref(3, 2)*x[2] + 
			5*qref(1, 3)*x[0]*x[1]*x[2] + 5*qref(3, 1)*x[0]*x[1]*x[2] + 
			5*qref(2, 3)*sqr(x[1])*x[2] + 5*qref(3, 2)*sqr(x[1])*x[2] + 
			5*qref(3, 3)*x[1]*sqr(x[2]))*invr7;
}

double TB0::QuadrupoleBzTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr7 = invr*invr2*sqr(invr2);
	const double r2 = sqr(r);
	//recflops(63);
	return (-(r2*qref(1, 3)*x[0]) - r2*qref(3, 1)*x[0] - 
			r2*qref(2, 3)*x[1] - r2*qref(3, 2)*x[1] - 2*r2*qref(3, 3)*x[2] + 
			5*qref(1, 1)*sqr(x[0])*x[2] + 5*qref(1, 2)*x[0]*x[1]*x[2] + 
			5*qref(2, 1)*x[0]*x[1]*x[2] + 5*qref(2, 2)*sqr(x[1])*x[2] + 
			5*qref(1, 3)*x[0]*sqr(x[2]) + 5*qref(3, 1)*x[0]*sqr(x[2]) + 
			5*qref(2, 3)*x[1]*sqr(x[2]) + 5*qref(3, 2)*x[1]*sqr(x[2]) + 
			5*qref(3, 3)*pow3(x[3]))*invr7;
}

#undef qref

#define qref(i,j,k) q[3][((i-1)*3+j-1)*3+k-1]

double TB0::OctupoleBxTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr9 = invr*sqr(invr2)*sqr(invr2);
	const double r2 = sqr(r);
	const double x02 = sqr(x[0]);
	const double x03 = x[0]*x02;
	const double x04 = sqr(x02);
	const double x12 = sqr(x[1]);
	const double x22 = sqr(x[2]);
	const double x23 = x[2]*x22;
	const double x13 = x[1]*x12;

	//recflops(95);
	
	return (-3*r2*qref(1, 1, 1)*x02
			- 2*r2*((qref(1, 1, 2) + qref(1, 2, 1) + qref(2, 1, 1))*x[0]*x[1]
					+ (qref(1, 1, 3) + qref(1, 3, 1) + qref(3, 1, 1))*x[0]*x[2])
			- r2*((qref(1, 2, 2) + qref(2, 1, 2) + qref(2, 2, 1))*x12
				  + (qref(1, 2, 3) + qref(1, 3, 2) + qref(2, 1, 3) + qref(2, 3, 1) + qref(3, 1, 2) + qref(3, 2, 1))*x[1]*x[2]
				  + (qref(1, 3, 3) + qref(3, 1, 3) + qref(3, 3, 1))*x22)
			+ 7*(qref(1, 1, 1)*x04
				 + (qref(1, 1, 2) + qref(1, 2, 1) + qref(2, 1, 1))*x03*x[1]
				 + (qref(1, 2, 2) + qref(2, 1, 2) + qref(2, 2, 1))*x02*x12
				 + qref(2, 2, 2)*x[0]*x13
				 + (qref(1, 1, 3) + qref(1, 3, 1) + qref(3, 1, 1))*x03*x[2]
				 + (qref(1, 2, 3) + qref(1, 3, 2) + qref(2, 1, 3) + qref(2, 3, 1) + qref(3, 1, 2) + qref(3, 2, 1))*x02*x[1]*x[2]
				 + (qref(2, 2, 3) + qref(2, 3, 2) + qref(3, 2, 2))*x[0]*x12*x[2]
				 + (qref(1, 3, 3) + qref(3, 1, 3) + qref(3, 3, 1))*x02*x22
				 + (qref(2, 3, 3) + qref(3, 2, 3) + qref(3, 3, 2))*x[0]*x[1]*x22
				 + qref(3, 3, 3)*x[0]*x23
				)
		)*invr9;
}

double TB0::OctupoleByTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr9 = invr*sqr(invr2)*sqr(invr2);
	const double r2 = sqr(r);
	const double x02 = sqr(x[0]);
	const double x03 = x[0]*x02;
//	const double x04 = sqr(x02);
	const double x12 = sqr(x[1]);
	const double x22 = sqr(x[2]);
	const double x23 = x[2]*x22;
	const double x13 = x[1]*x12;
	const double x14 = sqr(x12);

	//recflops(95);

	return (
			- 3*r2*qref(2, 2, 2)*x12
			- 2*r2*((qref(1, 2, 2) + qref(2, 1, 2) + qref(2, 2, 1))*x[0]*x[1]
					+ (qref(2, 2, 3) + qref(2, 3, 2) + qref(3, 2, 2))*x[1]*x[2])
		    - r2*((qref(1, 1, 2) + qref(1, 2, 1) + qref(2, 1, 1))*x02
				  + (qref(2, 3, 3) + qref(3, 2, 3) + qref(3, 3, 2))*x22
				  + (qref(1, 2, 3) + qref(1, 3, 2) + qref(2, 1, 3)
					 + qref(2, 3, 1) + qref(3, 1, 2) + qref(3, 2, 1))*x[0]*x[2])
			+ 7*(qref(1, 1, 1)*x03*x[1]
				 + qref(2, 2, 2)*x14
				 + qref(3, 3, 3)*x[1]*x23
				 + (qref(1, 1, 2) + qref(1, 2, 1) + qref(2, 1, 1))*x02*x12
				 + (qref(1, 2, 2) + qref(2, 1, 2) + qref(2, 2, 1))*x[0]*x13
				 + (qref(1, 1, 3) + qref(1, 3, 1) + qref(3, 1, 1))*x02*x[1]*x[2]
				 + (qref(1, 2, 3) + qref(1, 3, 2) + qref(2, 1, 3))*x[0]*x12*x[2]
				 + (qref(2, 3, 1) + qref(3, 1, 2) + qref(3, 2, 1))*x[0]*x12*x[2]
				 + (qref(2, 2, 3) + qref(2, 3, 2) + qref(3, 2, 2))*x13*x[2]
				 + (qref(1, 3, 3) + qref(3, 1, 3) + qref(3, 3, 1))*x[0]*x[1]*x22
				 + (qref(2, 3, 3) + qref(3, 2, 3) + qref(3, 3, 2))*x12*x22
				)
		)*invr9;
}

double TB0::OctupoleBzTerms(const double x[3], double r, double invr, double invr2) const
{
	const double invr9 = invr*sqr(invr2)*sqr(invr2);
	const double r2 = sqr(r);
	const double x02 = sqr(x[0]);
	const double x03 = x[0]*x02;
//	const double x04 = sqr(x02);
	const double x12 = sqr(x[1]);
	const double x22 = sqr(x[2]);
	const double x23 = x[2]*x22;
	const double x24 = sqr(x22);
	const double x13 = x[1]*x12;

	//recflops(95);

	return (-(r2*qref(1, 1, 3)*x02) - r2*qref(1, 3, 1)*x02 - 
     r2*qref(3, 1, 1)*x02 - r2*qref(1, 2, 3)*x[0]*x[1] - 
     r2*qref(1, 3, 2)*x[0]*x[1] - r2*qref(2, 1, 3)*x[0]*x[1] - 
     r2*qref(2, 3, 1)*x[0]*x[1] - r2*qref(3, 1, 2)*x[0]*x[1] - 
     r2*qref(3, 2, 1)*x[0]*x[1] - r2*qref(2, 2, 3)*x12 - 
     r2*qref(2, 3, 2)*x12 - r2*qref(3, 2, 2)*x12 - 
     2*r2*qref(1, 3, 3)*x[0]*x[2] - 2*r2*qref(3, 1, 3)*x[0]*x[2] - 
     2*r2*qref(3, 3, 1)*x[0]*x[2] + 7*qref(1, 1, 1)*x03*x[2] - 
     2*r2*qref(2, 3, 3)*x[1]*x[2] - 2*r2*qref(3, 2, 3)*x[1]*x[2] - 
     2*r2*qref(3, 3, 2)*x[1]*x[2] + 7*qref(1, 1, 2)*x02*x[1]*x[2] + 
     7*qref(1, 2, 1)*x02*x[1]*x[2] + 
     7*qref(2, 1, 1)*x02*x[1]*x[2] + 
     7*qref(1, 2, 2)*x[0]*x12*x[2] + 
     7*qref(2, 1, 2)*x[0]*x12*x[2] +
     7*qref(2, 2, 1)*x[0]*x12*x[2] + 7*qref(2, 2, 2)*x13*x[2] - 
     3*r2*qref(3, 3, 3)*x22 + 7*qref(1, 1, 3)*x02*x22 + 
     7*qref(1, 3, 1)*x02*x22 + 7*qref(3, 1, 1)*x02*x22 + 
     7*qref(1, 2, 3)*x[0]*x[1]*x22 + 
     7*qref(1, 3, 2)*x[0]*x[1]*x22 + 
     7*qref(2, 1, 3)*x[0]*x[1]*x22 + 
     7*qref(2, 3, 1)*x[0]*x[1]*x22 + 
     7*qref(3, 1, 2)*x[0]*x[1]*x22 + 
     7*qref(3, 2, 1)*x[0]*x[1]*x22 + 7*qref(2, 2, 3)*x12*x22 + 
     7*qref(2, 3, 2)*x12*x22 + 7*qref(3, 2, 2)*x12*x22 + 
     7*qref(1, 3, 3)*x[0]*x23 + 7*qref(3, 1, 3)*x[0]*x23 + 
     7*qref(3, 3, 1)*x[0]*x23 + 7*qref(2, 3, 3)*x[1]*x23 + 
     7*qref(3, 2, 3)*x[1]*x23 + 7*qref(3, 3, 2)*x[1]*x23 + 
     7*qref(3, 3, 3)*x24)*invr9;
	
}

#undef qref



double TB0::EvalBx(double x, double y, double z) const
{
	double Bx = 0;
	x-= center_x;
	y-= center_y;
	z-= center_z;
	const double r2 = max(sqr(minimum_r), sqr(x) + sqr(y) + sqr(z));
	const double r = sqrt(r2);
	const double invr2 = 1.0/r2;
	const double invr = 1.0/r;
	const double R[3] = {x,y,z};
	//recflops(9+flops_sqrt+2*flops_div+1);
	if (q[0][0] != 0)
		Bx+= MonopoleBxTerms(R,r,invr,invr2);
	if (this->maxorder >= 1) {
		Bx+= DipoleBxTerms(R,r,invr,invr2);
		if (this->maxorder >= 2) {
			Bx+= QuadrupoleBxTerms(R,r,invr,invr2);
			if (this->maxorder >= 3)
				Bx+= OctupoleBxTerms(R,r,invr,invr2);
		}
	}
	return Bx;
}

double TB0::EvalBy(double x, double y, double z) const
{
	double By = 0;
	x-= center_x;
	y-= center_y;
	z-= center_z;
	const double r2 = max(sqr(minimum_r), sqr(x) + sqr(y) + sqr(z));
	const double r = sqrt(r2);
	const double invr2 = 1.0/r2;
	const double invr = 1.0/r;
	const double R[3] = {x,y,z};
	//recflops(9+flops_sqrt+2*flops_div+1);
	if (q[0][0] != 0)
		By+= MonopoleByTerms(R,r,invr,invr2);
	if (this->maxorder >= 1) {
		By+= DipoleByTerms(R,r,invr,invr2);
		if (this->maxorder >= 2) {
			By+= QuadrupoleByTerms(R,r,invr,invr2);
			if (this->maxorder >= 3)
				By+= OctupoleByTerms(R,r,invr,invr2);
		}
	}
	return By;
}

double TB0::EvalBz(double x, double y, double z) const
{
	double Bz = 0;
	x-= center_x;
	y-= center_y;
	z-= center_z;
	const double r2 = max(sqr(minimum_r), sqr(x) + sqr(y) + sqr(z));
	const double r = sqrt(r2);
	const double invr2 = 1.0/r2;
	const double invr = 1.0/r;
	const double R[3] = {x,y,z};
	//recflops(9+flops_sqrt+2*flops_div+1);
	if (q[0][0] != 0)
		Bz+= MonopoleBzTerms(R,r,invr,invr2);
	if (this->maxorder >= 1) {
		Bz+= DipoleBzTerms(R,r,invr,invr2);
		if (this->maxorder >= 2) {
			Bz+= QuadrupoleBzTerms(R,r,invr,invr2);
			if (this->maxorder >= 3)
				Bz+= OctupoleBzTerms(R,r,invr,invr2);
		}
	}
	return Bz;
}

double TB0::call(double x, double y, double z) const
{
	switch (this->eval_comp) {
	case EVAL_BX: return EvalBx(x,y,z);
	case EVAL_BY: return EvalBy(x,y,z);
	case EVAL_BZ: return EvalBz(x,y,z);
	case EVAL_PHI: {double r[3] = {x,y,z}; return EvaluatePhi(r);}
	}
	return 0;	// dummy, but prevents gcc from yelling
}

void TB0::BackgroundEvaluate(const double r[3], double& Bx, double& By, double& Bz) const
{
	if (all_qs_are_zero) {
		Bx = By = Bz = 0;
	} else {
		Bx = EvalBx(r[0],r[1],r[2]);
		By = EvalBy(r[0],r[1],r[2]);
		Bz = EvalBz(r[0],r[1],r[2]);
	}
}


void TB0::BackgroundLineAverage(const double r1[3], unsigned short d, double L, double& Bx, double& By, double& Bz)
{
	if (all_qs_are_zero) {
		Bx = By = Bz = 0;
		return;
	}
	const double norm = 1/L;
	const double acc = Baccuracy*L;
	const double a = r1[d];
	const double b = r1[d] + L;
	T1DFunction *f=0;
	switch (d) {
	case 0: f = new T3D_fix23(*this,r1[1],r1[2]); break;
	case 1: f = new T3D_fix13(*this,r1[0],r1[2]); break;
	case 2: f = new T3D_fix12(*this,r1[0],r1[1]); break;
	default:
		cerr << "*** TB0::BackgroundLineAverage: d=" << d << " is bad\n";
		exit(1);
	}
	this->eval_comp = EVAL_BX; Bx = Romberg(*f,a,b,acc)*norm;
	this->eval_comp = EVAL_BY; By = Romberg(*f,a,b,acc)*norm;
	this->eval_comp = EVAL_BZ; Bz = Romberg(*f,a,b,acc)*norm;
	delete f;

	Bx += this->const_Bx0;
	By += this->const_By0;
	Bz += this->const_Bz0;
}


void TB0::BackgroundSurfaceAverage(
	const double r1[3],
	unsigned short d,
	double L1,
	double L2,
	double& Bx, double& By, double& Bz
) {
	if (all_qs_are_zero) {
		Bx = By = Bz = 0;
		return;
	}
	const double acc = Baccuracy*L1*L2;
	const double norm = 1/(L1*L2);
	switch (d) {
	case 0:
	{
		T3D_fix1 f(*this,r1[0]);
		this->eval_comp = EVAL_BX; Bx = Romberg(f, r1[1],r1[1]+L1, r1[2],r1[2]+L2, acc)*norm;
		this->eval_comp = EVAL_BY; By = Romberg(f, r1[1],r1[1]+L1, r1[2],r1[2]+L2, acc)*norm;
		this->eval_comp = EVAL_BZ; Bz = Romberg(f, r1[1],r1[1]+L1, r1[2],r1[2]+L2, acc)*norm;
	}
		break;
	case 1:
	{
		T3D_fix2 f(*this,r1[1]);
		this->eval_comp = EVAL_BX; Bx = Romberg(f, r1[0],r1[0]+L1, r1[2],r1[2]+L2, acc)*norm;
		this->eval_comp = EVAL_BY; By = Romberg(f, r1[0],r1[0]+L1, r1[2],r1[2]+L2, acc)*norm;
		this->eval_comp = EVAL_BZ; Bz = Romberg(f, r1[0],r1[0]+L1, r1[2],r1[2]+L2, acc)*norm;
	}
		break;
	case 2:
	{
		T3D_fix3 f(*this,r1[2]);
		this->eval_comp = EVAL_BX; Bx = Romberg(f, r1[0],r1[0]+L1, r1[1],r1[1]+L2, acc)*norm;
		this->eval_comp = EVAL_BY; By = Romberg(f, r1[0],r1[0]+L1, r1[1],r1[1]+L2, acc)*norm;
		this->eval_comp = EVAL_BZ; Bz = Romberg(f, r1[0],r1[0]+L1, r1[1],r1[1]+L2, acc)*norm;
	}
		break;
	default:
		cerr << "*** TB0::BackgroundSurfaceAverage: d=" << d << " is bad\n";
		exit(1);
	}

	Bx += this->const_Bx0;
	By += this->const_By0;
	Bz += this->const_Bz0;
}


void TB0::BackgroundVolumeAverageFast(const double r1[3], const double r2[3],
									  double& Bx, double& By, double& Bz)
{
	if (all_qs_are_zero) {
		Bx = By = Bz = 0;
		return;
	}
	const double L = r2[0] - r1[0];
	// Since B ~ phi/r, phi accuracy is r*Baccuracy. We adopt rmin as representing r here.
	const double rr1 = sqrt(sqr(r1[0]-center_x) + sqr(r1[1]-center_y) + sqr(r1[2]-center_z));
	const double rr2 = sqrt(sqr(r2[0]-center_x) + sqr(r2[1]-center_y) + sqr(r2[2]-center_z));
	const double rmin = max(minimum_r,min(rr1,rr2));
	const double acc = rmin*Baccuracy*sqr(L);
	const double norm = 1/(L*sqr(L));		// 1/V

	this->eval_comp = EVAL_PHI;

	const double Bxm = Romberg(T3D_fix1(*this,r1[0]), r1[1],r2[1], r1[2],r2[2], acc)*norm;
	const double Bxp = Romberg(T3D_fix1(*this,r2[0]), r1[1],r2[1], r1[2],r2[2], acc)*norm;
	Bx = Bxm - Bxp;

	const double Bym = Romberg(T3D_fix2(*this,r1[1]), r1[0],r2[0], r1[2],r2[2], acc)*norm;
	const double Byp = Romberg(T3D_fix2(*this,r2[1]), r1[0],r2[0], r1[2],r2[2], acc)*norm;
	By = Bym - Byp;

	const double Bzm = Romberg(T3D_fix3(*this,r1[2]), r1[0],r2[0], r1[1],r2[1], acc)*norm;
	const double Bzp = Romberg(T3D_fix3(*this,r2[2]), r1[0],r2[0], r1[1],r2[1], acc)*norm;
	Bz = Bzm - Bzp;

	Bx += this->const_Bx0;
	By += this->const_By0;
	Bz += this->const_Bz0;
}

void TB0::BackgroundVolumeAverage(const double r1[3], const double r2[3],
								  double& Bx, double& By, double& Bz)
{
	if (all_qs_are_zero) {
		Bx = By = Bz = 0;
		return;
	}
	const double acc = Baccuracy*(r2[0]-r1[0])*(r2[1]-r1[1])*(r2[2]-r1[2]);
	const double norm = 1.0/((r2[0]-r1[0])*(r2[1]-r1[1])*(r2[2]-r1[2]));
	this->eval_comp = EVAL_BX; Bx = Romberg(*this, r1[0],r2[0], r1[1],r2[1], r1[2],r2[2], acc)*norm;
	this->eval_comp = EVAL_BY; By = Romberg(*this, r1[0],r2[0], r1[1],r2[1], r1[2],r2[2], acc)*norm;
	this->eval_comp = EVAL_BZ; Bz = Romberg(*this, r1[0],r2[0], r1[1],r2[1], r1[2],r2[2], acc)*norm;

	Bx += this->const_Bx0;
	By += this->const_By0;
	Bz += this->const_Bz0;
}

void TB0::derivs(double, const double y[], double dydx[])
{
	double Bx,By,Bz, bx,by,bz;
	BackgroundEvaluate(y, Bx,By,Bz);
	const double B2 = sqr(Bx) + sqr(By) + sqr(Bz);
	const double invB = 1.0/sqrt(B2);
	bx = Bx*invB; by = By*invB; bz = Bz*invB;
	dydx[0] = integration_direction*bx;
	dydx[1] = integration_direction*by;
	dydx[2] = integration_direction*bz;
	const double norm2 = sqr(dydx[0]) + sqr(dydx[1]) + sqr(dydx[2]);
	maxnorm2 = (norm2 > maxnorm2) ? norm2 : maxnorm2;
}

bool TB0::TraceFieldLine(const double r1[3], double r1len, double sgn, double rs, double condsgn, double result[3])
// Trace field line from point r1 in direction sgn. r1len must be equal to |r1|.
// Stop when |r|*condsgn >= rs*condsgn. Iterate back until |r| == rs as closely as possible.
// If |r|*condsgn < |r1|*condsgn, stop and return false. On success return true.
{
	const double htry_vs_rlen = 0.9;
	const double eps = 1e-9;
	const double hepsilon = 0.01;
	const double stepreduction = 0.5;
	int d;
	integration_direction = sgn;
	double r[3],drds[3],rlen=r1len,rold[3];
	for (d=0; d<3; d++) r[d] = rold[d] = r1[d];
	TB0::derivs(0,r,drds);
	double s = 0;
	const double yscal[3] = {r1len,r1len,r1len};
	double htry = htry_vs_rlen*(r1len < rs ? r1len : rs);
	const double hthreshold = hepsilon*htry;
	const double rs2 = sqr(rs);
	double hdid,hnext;
	maxnorm2 = 0;
beginning:
	while (1) {
		BulirschStoerStep(r,drds,3,s, htry,eps, yscal, hdid,hnext, *this);
		if (hdid < htry) cerr << "TB0::TraceFieldLine: Bulirsch-Stoer step hdid < htry\n";
		//cout << "BulirschStoerStep(htry=" << htry << ", hdid=" << hdid << ", hnext=" << hnext << ")\n";
		if (condsgn > 0 && maxnorm2 >= rs2) break;
		rlen = sqrt( sqr(r[0]) + sqr(r[1]) + sqr(r[2]) );
		if (rlen*condsgn >= rs*condsgn) break;
		if (rlen*condsgn < r1len*condsgn) return false;
		for (d=0; d<3; d++) rold[d] = r[d];
	}
	// now rold is inside, r is outside
	if (hdid > hthreshold) {
		htry = stepreduction*hdid;		// reduce leap size and restart from rold
		for (d=0; d<3; d++) r[d] = rold[d];
		goto beginning;
	}
	const double rlenold = sqrt( sqr(rold[0]) + sqr(rold[1]) + sqr(rold[2]) );
	if (rlen == rlenold)
		for (d=0; d<3; d++) result[d] = r[d];
	else {
		assert(rlenold*condsgn < rs*condsgn);
		assert(rlen*condsgn >= rs*condsgn);
		// (1-t)*rlenold + t*rlen == rs
		// rlenold + t*(rlen-rlenold) == rs
		// t == (rs - rlenold)/(rlen - rlenold)
		const double t = (rs - rlenold)/(rlen - rlenold);
		for (d=0; d<3; d++) result[d] = (1-t)*rold[d] + t*r[d];
	}
	return true;
}

bool TB0::TraceFLToShell(const double r1[3], double rs, double result[3])
{
	const double r1len = sqrt( sqr(r1[0]) + sqr(r1[1]) + sqr(r1[2]) );
	double Bx1,By1,Bz1;
	BackgroundEvaluate(r1, Bx1,By1,Bz1);
	const double B_dot_er = Bx1*r1[0] + By1*r1[1] + Bz1*r1[2];
	double sgn,condsgn;
	if (rs > r1len) {
		// Integrate towards outer shell
		sgn = (B_dot_er >= 0) ? +1 : -1;
		condsgn = +1;
	} else if (rs < r1len) {
		// Integrate towards inner shell
		sgn = (B_dot_er <= 0) ? +1 : -1;
		condsgn = -1;
	} else {
		cerr << "*** TB0::TraceFLToShell: rs == |r1|\n";
		return false;
	}
	return TraceFieldLine(r1,r1len,sgn,rs,condsgn, result);
}

