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
#include "particles.h"
#include "physconst.h"
#include "relativistic_math.h"
#include "vectorclass.h"
#include "vector3d.h"
#include "vlsv_writer.h"


/* Particle propagation given E- and B-Field at the particle location
 * with the Boris-Method */
void Particle::push(Vec3d& B, Vec3d& E, double dt) {

   Vec3d uminus = v + (q * E * dt)/(2. * m);
   Vec3d h = (q * B * dt)/(2. * m * gamma(uminus));
   Vec3d uprime = uminus + cross_product(uminus, h);
   h = (2.* h)/(1. + dot_product(h,h));
   Vec3d uplus = uminus + cross_product(uprime, h);

   v = uplus + (q * E * dt)/(2. * m);
   x += dt * v;
}

void writeParticles(ParticleContainer& p,const char* filename) {

   vlsv::Writer vlsvWriter;
   vlsvWriter.open(filename,MPI_COMM_WORLD,0);

   std::vector<double> writebuf(p.size() * 3);

   /* First, store particle positions */
   uint writable_particles=0;
   for(unsigned int i=0; i < p.size(); i++) {
      if(vector_length(p[i].x) == 0) {
        continue;
      }

      p[i].x.store(&(writebuf[3*writable_particles]));
      writable_particles++;
   }

   std::map<std::string,std::string> attribs;
   attribs["name"] = "proton_position";
   attribs["type"] = vlsv::mesh::STRING_POINT;
   if (vlsvWriter.writeArray("MESH",attribs,writable_particles,3,writebuf.data()) == false) {
      std::cerr << "\t ERROR failed to write particle positions!" << std::endl;
   }

   /* Then, velocities */
   writable_particles=0;
   for(unsigned int i=0; i < p.size(); i++) {
      if(vector_length(p[i].x) == 0) {
        continue;
      }
      p[i].v.store(&(writebuf[3*writable_particles]));
      writable_particles++;
   }

   attribs["name"] = "proton_velocity";
   if (vlsvWriter.writeArray("MESH",attribs,writable_particles,3,writebuf.data()) == false) {
      std::cerr << "\t ERROR failed to write particle velocities!" << std::endl;
   }
   vlsvWriter.close();
}
