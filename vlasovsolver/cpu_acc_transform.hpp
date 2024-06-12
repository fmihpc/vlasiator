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

#ifndef CPU_ACC_TRANSFORM_H
#define CPU_ACC_TRANSFORM_H

#include <Eigen/Geometry>
#include <Eigen/Core>

#include "../common.h"
#include "../spatial_cell_wrapper.hpp"

Eigen::Transform<Real,3,Eigen::Affine> compute_acceleration_transformation(
   spatial_cell::SpatialCell* spatial_cell,const uint popID,const Real& dt);

void updateAccelerationMaxdt(spatial_cell::SpatialCell* spatial_cell, 
                             const uint popID);



#endif
