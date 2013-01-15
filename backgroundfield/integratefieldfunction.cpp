/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011, 2012 Finnish Meteorological Institute
*/

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "../common.h"
#include "integratefieldfunction.hpp"

static int PowerOfThree(int n)		// return 3^n
{
	int i, result = 1;
	for (i=0; i<n; i++) result*= 3;
	return result;
}


void IntegrateFieldFunction::initialize(int maxorder1)
{
   this->initialized = true;
   this->accuracy = 1e-17;
   
}



void IntegrateFieldFunction::LineAverage(const double r1[3], unsigned short d, double L, double& Bx, double& By, double& Bz)
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
		cerr << "*** IntegrateFieldFunction::LineAverage: d=" << d << " is bad\n";
		exit(1);
	}
	this->eval_comp = EVAL_BX; Bx = Romberg(*f,a,b,acc)*norm;
	this->eval_comp = EVAL_BY; By = Romberg(*f,a,b,acc)*norm;
	this->eval_comp = EVAL_BZ; Bz = Romberg(*f,a,b,acc)*norm;
	delete f;
}


void IntegrateFieldFunction::SurfaceAverage(
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
		cerr << "*** IntegrateFieldFunction::SurfaceAverage: d=" << d << " is bad\n";
		exit(1);
	}

}


void IntegrateFieldFunction::VolumeAverage(const double r1[3], const double r2[3],
                                           double& Bx, double& By, double& Bz)
{

	const double acc = accuracy*(r2[0]-r1[0])*(r2[1]-r1[1])*(r2[2]-r1[2]);
	const double norm = 1.0/((r2[0]-r1[0])*(r2[1]-r1[1])*(r2[2]-r1[2]));
	this->eval_comp = EVAL_BX; Bx = Romberg(*this, r1[0],r2[0], r1[1],r2[1], r1[2],r2[2], acc)*norm;
	this->eval_comp = EVAL_BY; By = Romberg(*this, r1[0],r2[0], r1[1],r2[1], r1[2],r2[2], acc)*norm;
	this->eval_comp = EVAL_BZ; Bz = Romberg(*this, r1[0],r2[0], r1[1],r2[1], r1[2],r2[2], acc)*norm;

}

