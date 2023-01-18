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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*
Background magnetic field class of Vlasiator.
*/

#ifndef LINEDIPOLE_HPP
#define LINEDIPOLE_HPP
#include "fieldfunction.hpp"



class LineDipole {
private:
   bool initialized = false;
   double q[3];                  // Dipole moment; set to (0,0,moment)
   double center[3]; // Coordinates where the dipole sits; set to (0,0,0)
public:
  
  LineDipole() {}

   void initialize(const double moment, const double center_x, const double center_y, const double center_z);
  
   double operator()(double x, double y, double z, coordinate component, unsigned int derivative=0, coordinate dcomponent=X) const;
};

#endif

