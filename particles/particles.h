#pragma once
/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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
#include <vector>
#include "vectorclass.h"
#include "vector3d.h"
#include "../definitions.h"

struct Particle {
      Real m;
      Real q;
      Vec3d x;
      Vec3d v;

      Particle(Real mass, Real charge, Vec3d _x, Vec3d _v) :
         m(mass),q(charge),x(_x),v(_v) {}

      /* Particle propagation given E- and B-Field at the particle location
       * with the Boris-Method */
      void push(Vec3d& B, Vec3d& E, double dt);
};


void writeParticles(std::vector<Particle>& p, const char* filename);
