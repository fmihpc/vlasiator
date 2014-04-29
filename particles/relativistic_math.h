#pragma once
#include <glm/glm.hpp>
#include "physconst.h"

static double gamma(glm::dvec3 v) {
	double uq = glm::dot(v,v);
	return sqrt(1. + (uq / (PhysicalConstantsSI::c * PhysicalConstantsSI::c)));
}

/* Lorentz-boost */
static glm::dvec3 lorentz_boost(glm::dvec3& u, glm::dvec3& v) {
	/* Transform velocities into betas */
	glm::dvec3 uprime = u / PhysicalConstantsSI::c;
	glm::dvec3 vprime = v / PhysicalConstantsSI::c;

	double udotv = glm::dot(uprime,vprime);
	double vsqr = glm::dot(vprime,vprime);

	if(vsqr == 0) { return u; }

	glm::dvec3 uparallel = udotv / vsqr * vprime;
	glm::dvec3 uperp = uprime - uparallel;

	return PhysicalConstantsSI::c * (vprime + uparallel + sqrt(1 - vsqr) * uperp) / (1 + udotv);
}

/* Lorentz transform of E and B fields */
static glm::dvec3 lorentz_transformed_E(glm::dvec3& E, glm::dvec3& B, glm::dvec3& v) {

	glm::dvec3 nv = glm::normalize(v);
	double g = gamma(v);

	return g * (E + glm::cross(v,B)/PhysicalConstantsSI::c) - (g-1.)*glm::dot(E,nv)*nv;
}
static glm::dvec3 lorentz_transformed_B(glm::dvec3& E, glm::dvec3& B, glm::dvec3& v) {
	glm::dvec3 nv = glm::normalize(v);
	double g = gamma(v);

	return g * (B - glm::cross(v,E)/PhysicalConstantsSI::c) - (g-1.)*glm::dot(B,nv)*nv;
}
