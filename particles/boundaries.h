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

#include "particles.h"

struct Boundary
{
   // Handle a particle's boundary behaviour.
   // returns "true" if the particle is still part of the simulation
   // afterwards, or "false" if it is to be removed.
   virtual bool handleParticle(Particle& p) = 0;

   // Handle cell coordinate in the spatial dimension this boundary
   // object cares about (for example: wrap in a periodic direction)
   virtual int cellCoordinate(int c) = 0;

   Boundary(int _dimension) : dimension(_dimension) {};
   virtual void setExtent(double _min, double _max, int _cells) {
      min=_min;
      max=_max;
      cells=_cells;
   }

   virtual ~Boundary(){};

   // Which spatial dimension to handle
   int dimension;

   // Minimum and maximum spatial extents in this dimension
   double min, max;

   // Number of cells in this dimension
   int cells;
};

// Boundary for a spatial dimension that is only 1 cell thick (pseudo-periodic)
struct CompactSpatialDimension : public Boundary
{
   virtual bool handleParticle(Particle& p) {
      // This boundary does not affect particles
      return true;
   }
   virtual int cellCoordinate(int c) {
      // Actual cell coordinates in this direction are
      // always mapped to 0.
      return 0;
   }
   CompactSpatialDimension(int _dimension) : Boundary(_dimension){};
};

// Open boundary, which removes particles if they fly out
struct OpenBoundary : public Boundary
{
   virtual bool handleParticle(Particle& p) {

      // Delete particles that are outside our boundaries.
      if(p.x[dimension] <= min || p.x[dimension] >= max) {
         return false;
      } else {
         // Leave all others be.
         return true;
      }
   }

   virtual int cellCoordinate(int c) {
      // Cell coordinates are clamped
      // TODO: Should this print warnings?
      if(c < 0) {
         return 0;
      } else if(c >= cells) {
         return cells-1;
      } else {
         return c;
      }
   }

   virtual void setExtent(double _min, double _max, int _cells) {
      double dx = (_max-_min)/((double)_cells);
      min=_min+2*dx; // 2 Cells border.
      max=_max-2*dx;
      cells=_cells;
   }
   OpenBoundary(int _dimension) : Boundary(_dimension){};
};

struct ReflectBoundary : public Boundary
{
   virtual bool handleParticle(Particle& p) {
      // Particles outside of bounds get their velocities flipped
      if(p.x[dimension] <= min || p.x[dimension] >= max) {
         p.v *= flip_v;
      }
      return true;
   }

   virtual int cellCoordinate(int c) {
      // Cell coordinates are clamped
      // TODO: Should this print warnings?
      if(c < 0) {
         return 0;
      } else if(c >= cells) {
         return cells-1;
      } else {
         return c;
      }
   }

   // Constructor
   ReflectBoundary(int _dimension) : Boundary(_dimension) {
      double flip[3] = {1.,1.,1.};
      flip[dimension] = -1.;
      flip_v.load(flip);
   }
   virtual void setExtent(double _min, double _max, int _cells) {
      double dx = (_max-_min)/((double)_cells);
      min=_min+2*dx; // 2 Cells border.
      max=_max-2*dx;
      cells=_cells;
   }

   // Vector to multiply with in order to flip velocity
   // vectors for our dimension
   Vec3d flip_v;

};

struct PeriodicBoundary : public Boundary
{
   virtual bool handleParticle(Particle& p) {
      if(p.x[dimension] < min) {
         p.x += offset_p;
      } else if(p.x[dimension] >= max) {
         p.x -= offset_p;
      }
      return true;
   }

   virtual int cellCoordinate(int c) {
      return c % cells;
   }

   // Constructor
   PeriodicBoundary(int _dimension) : Boundary(_dimension) {
   }
   virtual void setExtent(double _min, double _max, int _cells) {
      min=_min;
      max=_max;
      cells=_cells;

      double offset[3] = {0.,0.,0.};
      offset[dimension] = max-min;
      offset_p.load(offset);
   }

   // Vector to offset particle positions that leave through
   // one boundary with, to come out the other end
   Vec3d offset_p;

};

template<typename T> Boundary* createBoundary(int dimension)
{
   return new T(dimension);
}
