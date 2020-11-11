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
#include <vector>
#include "vectorclass.h"
#include "vector3d.h"
#include "../definitions.h"
#include "../memoryallocation.h"

struct Particle {
      Vec3d x;
      Vec3d v;
      Real m;
      Real q;
      char padding[128-sizeof(Vec3d)*2-sizeof(Real)*2];

      Particle(Real mass, Real charge, const Vec3d& _x, const Vec3d& _v) :
         x(_x),v(_v),m(mass),q(charge) {}

      /* Particle propagation given E- and B-Field at the particle location
       * with the Boris-Method */
      void push(Vec3d& B, Vec3d& E, double dt);
};

typedef std::vector<Particle, aligned_allocator<Particle, 32>> ParticleContainer;

void writeParticles(ParticleContainer& p, const char* filename);

