/*
This file is part of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2011 Finnish Meteorological Institute
*/

#include "cassert"
#include "cmath"
#include "cstdlib"
#include "iostream"

#include "ode.hpp"

static void RationalExtrapolation(int iest, double xest, const double yest[], double yz[], double dy[], int nv, int nuse)
{
	const int imax = 11;
	const int nmax = 10;
	const int ncol = 7;
	assert(0 <= iest && iest < imax);
	assert(1 <= nuse && nuse <= ncol);
	assert(1 <= nv && nv < nmax);
	int j,k;
	double V,B,B1,C,yy,ddy=0;
	static double x[imax],D[nmax][ncol],fx[ncol];
	x[iest] = xest;
	if (iest == 0) {
		for (j=0; j<nv; j++) {
			yz[j] = yest[j];
			D[j][0] = yest[j];
			dy[j] = yest[j];
		}
	} else {
		const int m1 = (iest < nuse) ? iest : nuse;		// use at most nuse previous members
		for (k=0; k<m1-1; k++)
			fx[k+1] = x[iest-k]/xest;
		for (j=0; j<nv; j++) {	// evaluate next diagonal in tableau
			yy = yest[j];
			V = D[j][0];
			C = yy;
			D[j][0] = yy;
			for (k=1; k<m1; k++) {
				B1 = fx[k]*V;
				B = B1 - C;
				if (B != 0) {
					B = (C-V)/B;
					ddy = C*B;
					C = B1*B;
				} else
					ddy = V;
				if (k != m1) V = D[j][k];
				D[j][k] = ddy;
				yy+= ddy;
			}
			dy[j] = ddy;
			yz[j] = yy;
		}
	}
}

static void ModifiedMidpointMethod
    (double y[],
	 const double dydx[],
	 int nvar,
	 double xs,
	 double htot,
	 int nstep,
	 double yout[],
	 TODERHS& rhs)
{
	const int nmax = 10;
	double ym[nmax],yn[nmax];
	int i,n;
	double h,h2,x,swap;
	h = htot/nstep;
	for (i=0; i<nvar; i++) {
		ym[i] = y[i];
		yn[i] = y[i] + h*dydx[i];		// first step
	}
	x = xs + h;
	rhs.derivs(x,yn,yout);		// will use yout for temporary storage of derivatives
	h2 = 2*h;
	for (n=2; n<=nstep; n++) {		// general step
		for (i=0; i<nvar; i++) {
			swap = ym[i] + h2*yout[i];
			ym[i] = yn[i];
			yn[i] = swap;
		}
		x+= h;
		rhs.derivs(x,yn,yout);
	}
	for (i=0; i<nvar; i++)
		yout[i] = 0.5*(ym[i] + yn[i] + h*yout[i]);
}

void BulirschStoerStep(
	double y[],
	const double dydx[],
	int nv,
	double &x,
	double htry,
	double eps,
	const double yscal[],
	double& hdid,
	double& hnext,
	TODERHS& rhs)
{
	const int nmax = 10;
	const int imax = 11;		// don't change, please, nseq vector hard-coded below
	const int nuse = 7;
	const double shrink = 0.95;
	const double grow = 1.2;
	assert(1 <= nv && nv < nmax);
	double yerr[nmax], ysav[nmax],dysav[nmax],yseq[nmax];
	static const int nseq[imax] = {2,4,6,8,12,16,24,32,48,64,96};
	int i,j;
	double h,xsav,xest,errmax,errmax1;
	h = htry;
	xsav = x;
	for (i=0; i<nv; i++) {	// save the starting values
		ysav[i] = y[i];
		dysav[i] = dydx[i];
	}
beginning:
	for (i=0; i<imax; i++) {		// evaluate the sequence of modified midpoint integrations
		ModifiedMidpointMethod(ysav,dysav,nv,xsav,h,nseq[i],yseq,rhs);
		xest = h/nseq[i];
		xest = xest*xest;
		RationalExtrapolation(i,xest,yseq,y,yerr,nv,nuse);
		errmax = 0;
		for (j=0; j<nv; j++) {
			errmax1 = fabs(yerr[j]/yscal[j]);
			errmax = (errmax1 > errmax) ? errmax1 : errmax;
		}
		errmax/= eps;		// scale accuracy relative to tolerance

		if (errmax < 1) {		// step converged
			x+= h;
			hdid = h;
			if (i == nuse-1)
				hnext = h*shrink;
			else if (i == nuse-2)
				hnext = h*grow;
			else
				hnext = (h*nseq[nuse-2])/nseq[i];
			return;		// normal return
		}
	}
	// if here, then step failed. We reduce the stepsize and try again.
	h = 0.25*h/(1 << ((imax-nuse)/2));
	if (x+h == x) {
		std::cerr << "*** ode.C:BSstep: Step size underflow\n";
		exit(1);
	}
	goto beginning;
}
