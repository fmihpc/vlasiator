/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef B0_HPP
#define B0_HPP

#include "boost/program_options.hpp"

#include "ode.hpp"
#include "quadr.hpp"

#define MAX_MULTIPOLE_ORDER 3	// 1:dipole, 2:quadrupole, 3:octupole,...
// The expressions have been programmed only up to octupole order, so do not change
// this definition unless you have added more members to TB0.

/*

  The magnetic potential expansion is as follows.
  The computed magnetic field B = -grad(Phi).
  The q,Q coefficients are the parameters set by setcoeff() below.
  
          q    q.r            xi xj              xi xj xk
   Phi = --- + ---- + sum Q   -----  +  sum Q    -------- + ...
          r    r^3     ij  ij  r^5      ijk  ijk   r^7

*/

enum Eval_Comp {
	EVAL_BX,
	EVAL_BY,
	EVAL_BZ,
	EVAL_PHI
};

class TB0: public TODERHS, public T3DFunction {
private:
	bool initialized;
	int maxorder;
	double const_Bx0, const_By0, const_Bz0;
	double *q[MAX_MULTIPOLE_ORDER+1];
	double EvaluateMultipoleTerm(int n, const double q[], const double r[3]) const;
	double MonopoleBxTerms(const double x[3], double r, double invr, double invr2) const;
	double MonopoleByTerms(const double x[3], double r, double invr, double invr2) const;
	double MonopoleBzTerms(const double x[3], double r, double invr, double invr2) const;
	double DipoleBxTerms(const double x[3], double r, double invr, double invr2) const;
	double DipoleByTerms(const double x[3], double r, double invr, double invr2) const;
	double DipoleBzTerms(const double x[3], double r, double invr, double invr2) const;
	double QuadrupoleBxTerms(const double x[3], double r, double invr, double invr2) const;
	double QuadrupoleByTerms(const double x[3], double r, double invr, double invr2) const;
	double QuadrupoleBzTerms(const double x[3], double r, double invr, double invr2) const;
	double OctupoleBxTerms(const double x[3], double r, double invr, double invr2) const;
	double OctupoleByTerms(const double x[3], double r, double invr, double invr2) const;
	double OctupoleBzTerms(const double x[3], double r, double invr, double invr2) const;
	void dealloc();
	// The magnetic potential Phi (B = -grad(Phi)) can be evaluated also
	double EvaluatePhi(const double r[3]) const;
	int all_qs_are_zero;
	// Parameters needed in derivs:
	double integration_direction;		// either +1 or -1
	double maxnorm2;					// records maximum norm^2 of y[] where derivs() is called
	// Our private field line tracer. Cannot be const since calls derivs (see below).
	bool TraceFieldLine(const double r1[3], double r1len, double sgn, double rs, double condsgn, double result[3]);
	// Switch for call(): will call either EvalBx, EvalBy, or EvalBz
	Eval_Comp eval_comp;
	double dipmom;
public:

	/*!
	Maxorder determines how many terms are included. The default is 1, which means
	that dipole (and monopole) terms are included, but not quadrupole or octupole.
	To include quadrupole, set maxorder to 2. To include also octupole, set it to 3.
	Alternatively, you can set it later by calling set_maxorder (but you must call
	set_maxorder before calling setcoeff or computing anything).
	Do NOT set maxorder higher than necessary, it will unnecessarily slow down
	all computations here.
	*/
	void initialize(int maxorder1=1);

	TB0()
	{
		this->initialized = false;
	}

	TB0(boost::program_options::options_description* options)
	{
		this->initialized = false;
		options->add_options()
			("constBx0",
				boost::program_options::value<double>(&(this->const_Bx0))->default_value(0),
				"Add a constant x component arg to the background magnetic field")
			("constBy0",
				boost::program_options::value<double>(&(this->const_By0))->default_value(0),
				"Add a constant y component arg to the background magnetic field")
			("constBz0",
				boost::program_options::value<double>(&(this->const_Bz0))->default_value(0),
				"Add a constant z component arg to the background magnetic field")
			("dipmom",
				boost::program_options::value<double>(&(this->dipmom))->default_value(8e15),
				"Magnitude of Earth's dipole moment (Tm^3, positive number)");
	}

	void set_dipole_moment(const double given_moment)
	{
		this->dipmom = given_moment;
	}

	void set_constant_Bx(const double given_Bx0)
	{
		this->const_Bx0 = given_Bx0;
	}

	void set_constant_By(const double given_By0)
	{
		this->const_By0 = given_By0;
	}

	void set_constant_Bz(const double given_Bz0)
	{
		this->const_Bz0 = given_Bz0;
	}

	void set_maxorder(int maxorder1) {dealloc(); this->initialize(maxorder1);}

	double minimum_r;	// singularity avoider in the origin; default 1e-10, but you can reassign it
	double Baccuracy;	// absolute accuracy to which to compute integrals etc., in Tesla, default 1e-10
	double center_x, center_y, center_z;	// coordinates where the dipole (quadrupole...) sits; default (0,0,0)

	void setcoeff(double q1);	// monopole coefficient, usually 0 (the default)
	void setcoeff(int i, double q1);	// setting of the three dipole coefficients (i=0,1,2)
	void setcoeff(int i, int j, double q1);	// quadrupole coefficients setting (nine of them)
	void setcoeff(const int ind[], int n1, double q1);	// general case; n1 is 1 for dipole, 2 for quadrupole etc. this form must be used to set the octupole coefficients; it can be used to set also monopole/dipole/quadrupole if wanted.

	// You can evaluate B components anywhere:
	double EvalBx(double x, double y, double z) const;
	double EvalBy(double x, double y, double z) const;
	double EvalBz(double x, double y, double z) const;

	// call() was virtually inherited from T3DFunction
	virtual double call(double x, double y, double z) const;

	/*!
	Or you can evaluate all components by this call
	*/
	void BackgroundEvaluate(const double r[3],
							double& Bx, double& By, double& Bz) const;

	/*!
	Average of B (all components) along a coordinate-aligned line starting from r1,
	having length L (can be negative) and proceeding to d'th coordinate
	*/
	void BackgroundLineAverage(const double r1[3], unsigned short d, double L,
							   double& Bx, double& By, double& Bz);

	/*!
	Average of B (all components) along a rectangular coordinate-aligned surface
	which is orthogonal to d'th coordinate, has lower left corner at r1,
	and surface side lengths (positive) equal to L1,L2 (either yz, xz or xy,
	depending on d).
	*/
	void BackgroundSurfaceAverage(const double r1[3], unsigned short d, double L1, double L2,
								  double& Bx, double& By, double& Bz);

	/*!
	Average of B (all components) over a rectangular coordinate-aligned volume
	having lower left corner at r1 and upper right corner at r2.
	*/
	void BackgroundVolumeAverage(const double r1[3], const double r2[3],
								 double& Bx, double& By, double& Bz);
	/*!
	Fast version of BackgroundVolumeAverage (uses surface integrals of magnetic potential Phi
	instead of directly computing volume integrals as the previous call does).
	*/
	void BackgroundVolumeAverageFast(const double r1[3], const double r2[3],
									 double& Bx, double& By, double& Bz);

	/*!
	derivs is passed to ODE integrator, it is virtually inherited from TODERHS
	*/
	void derivs(double x, const double y[], double dyds[]);

	/*!
	Starting from point r1, follow field line until it hits shell r=rs.
	Return the point in result. Chooses the direction in which to follow
	the field line automatically based on the sign of the radial component
	of B0 at r1. The shell rs can be either smaller or larger than |r1|,
	the routine detects this and takes the appropriate branch.
	In case of successful integration, return true.
	If the field line never reaches rs (this could happen only if rs > |r1|),
	return false. The result vector is not set in this case.
	This case is detected when the point returns to a distance which is
	smaller than |r1|.
	TraceFLToShell cannot be const as it must call derivs().
	*/
	bool TraceFLToShell(const double r1[3], double rs, double result[3]);
	
	virtual ~TB0() {dealloc();}
};

#endif

