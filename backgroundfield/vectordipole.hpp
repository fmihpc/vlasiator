/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * Copyright 2017-2019 University of Helsinki
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

#ifndef VECTORDIPOLE_HPP
#define VECTORDIPOLE_HPP
#include "fieldfunction.hpp"



class VectorDipole: public FieldFunction {
private:
   bool initialized;
   double q[3];      // Dipole moment; set to (0,0,moment) for z-aligned
   double center[3]; // Coordinates where the dipole sits; set to (0,0,0)
   double xlimit[2]; // X-coodrinate extents of full and zero dipole
   double IMF[3];    // IMF value to scale to, starting at xlimit[0] and finishing at xlimit[1]
public:
   
   VectorDipole(){
      this->initialized = false;
   }
   void initialize(const double moment,const double center_x, const double center_y, const double center_z, const double tilt_angle_phi, const double tilt_angle_theta, const double xlimit_f, const double xlimit_z, const double IMF_Bx, const double IMF_By, const double IMF_Bz);
   virtual double call(double x, double y, double z) const;  
   virtual ~VectorDipole() {}
};

#endif

