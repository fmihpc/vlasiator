#pragma once
#include "physconst.h"
#include "vectorclass.h"
#include "vector3d.h"

static double gamma(Vec3d v) {
	double uq = dot_product(v,v);
	return sqrt(1. + (uq / (PhysicalConstantsSI::c * PhysicalConstantsSI::c)));
}

/* Lorentz-boost */
static Vec3d lorentz_boost(Vec3d& u, Vec3d& v) {
	/* Transform velocities into betas */
	Vec3d uprime = u / PhysicalConstantsSI::c;
	Vec3d vprime = v / PhysicalConstantsSI::c;

	double udotv = dot_product(uprime,vprime);
	double vsqr = dot_product(vprime,vprime);

	if(vsqr == 0) { return u; }

	Vec3d uparallel = udotv / vsqr * vprime;
	Vec3d uperp = uprime - uparallel;

	return PhysicalConstantsSI::c * (vprime + uparallel + sqrt(1 - vsqr) * uperp) / (1 + udotv);
}

/* Lorentz transform of E and B fields */
static Vec3d lorentz_transformed_E(Vec3d& E, Vec3d& B, Vec3d& v) {

	Vec3d nv = normalize_vector(v);
	double g = gamma(v);

	return g * (E + cross_product(v,B)/PhysicalConstantsSI::c) - (g-1.)*dot_product(E,nv)*nv;
}
static Vec3d lorentz_transformed_B(Vec3d& E, Vec3d& B, Vec3d& v) {
	Vec3d nv = normalize_vector(v);
	double g = gamma(v);

	return g * (B - cross_product(v,E)/PhysicalConstantsSI::c) - (g-1.)*dot_product(B,nv)*nv;
}
