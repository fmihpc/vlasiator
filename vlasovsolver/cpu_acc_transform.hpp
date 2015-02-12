/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010-2015 Finnish Meteorological Institute
 * 
 * */

#ifndef CPU_ACC_TRANSFORM_H
#define CPU_ACC_TRANSFORM_H

#include <Eigen/Geometry>
#include <Eigen/Core>

#include "../common.h"
#include "../spatial_cell.hpp"

Eigen::Transform<Real,3,Eigen::Affine> compute_acceleration_transformation(
        spatial_cell::SpatialCell* spatial_cell,const int& popID,const Real& dt);

#endif
