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
#include <math.h>

// The default is (Gaussian) CGS units
class PhysicalConstantsCGS
{
   public:
      PhysicalConstantsCGS(){};

      static const double me;      // Electron mass (g)
      static const double mp;      // Proton mass (g)
      static const double mpi0;    // Mass of uncharged pion (g)
      static const double e;       // Elementary charge (statcoulomb)
      static const double c;       // Speed of light (cm/s)
      static const double eps0;    // Permitivity of vacuum
      static const double mu0;     // Permeability of vacuum
      static const double thomcs;  // Thomson cross section (cm^2)
      static const double h;       // Planck constant (erg sec)
      static const double hbar;    // reduced Planck constant (erg sec)
      static const double k;       // Boltzmann constant (erg / K)
      static const double r0;      // Classical electron radius / Thomson scattering length (cm)
      static const double G;       // Gravitational constant cm^3/(g sec^2)
};

// Constants in SI
class PhysicalConstantsSI
{
   public:
      PhysicalConstantsSI(){};

      static const double me;      // Electron mass (kg)
      static const double mp;      // Proton mass (kg)
      static const double mpi0;    // Mass of uncharged pion (kg)
      static const double e;       // Elementary charge (C)
      static const double c;       // Speed of light (m/s)
      static const double eps0;    // Permitivity of vacuum (F / m)
      static const double mu0;     // Permeability of vacuum (H / m)
      static const double thomcs;  // Thomson cross section (m^2)
      static const double h;       // Planck constant (J s)
      static const double hbar;    // reduced Planck constant (J s)
      static const double k;       // Boltzmann constant (J / K)
      static const double r0;      // Classical electron radius / Thomson scattering length (m)
      static const double G;       // Gravitational constant m^3/(kg s^2)
};

// Natural Units with e = c = k = me = 1 (Stoney units)
class PhysicalConstantsnorm
{
   public:
      PhysicalConstantsnorm(){};

      static const double me;      // Electron mass
      static const double mp;      // Proton mass
      static const double mpi0;    // Mass of uncharged pion
      static const double e;       // Elementary charge
      static const double c;       // Speed of light
      static const double eps0;    // Permitivity of vacuum
      static const double mu0;     // Permeability of vacuum
      static const double thomcs;  // Thomson cross section
      static const double h;       // Planck constant
      static const double hbar;    // reduced Planck constant
      static const double k;       // Boltzmann constant
      static const double r0;      // Classical electron radius
      static const double G;       // Gravitational constant

      // l = e^2 / c^2 me
      // t = e^2 / c^3 me
      // T = = me c^2 / k
};

