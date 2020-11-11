#pragma once
/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#include "physconst.h"
#include "vectorclass.h"
#include "vector3d.h"

static double gamma(Vec3d v) {
   double uq = dot_product(v,v);
   return sqrt(1. + (uq / (PhysicalConstantsSI::c * PhysicalConstantsSI::c)));
}

/* Lorentz-boost */
static Vec3d lorentzBoost(Vec3d& u, Vec3d& v) {
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
