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
#ifndef CPU_ACC_INTERSECTIONS_H
#define CPU_ACC_INTERSECTIONS_H

#include <Eigen/Core>

#include "../definitions.h"
#include "../spatial_cell_wrapper.hpp"

Eigen::Matrix<Real,3,1> line_plane_intersection(const Eigen::Matrix<Real,3,1>& l_point,
                                                const Eigen::Matrix<Real,3,1>& l_direction,
                                                const Eigen::Matrix<Real,3,1>& p_point,
                                                const Eigen::Matrix<Real,3,1>& p_normal);

void compute_intersections_1st(
        const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        const Eigen::Transform<Real,3,Eigen::Affine>& bwd_transform,
        const Eigen::Transform<Real,3,Eigen::Affine>& fwd_transform,
        uint dimension,const uint8_t& refLevel,
        Real& intersection,Real& intersection_di,
        Real& intersection_dj,Real& intersection_dk);

void compute_intersections_2nd(
        const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        const Eigen::Transform<Real,3,Eigen::Affine>& bwd_transform,
        const Eigen::Transform<Real,3,Eigen::Affine>& fwd_transform,
        uint dimension,const uint8_t& refLevel,
        Real& intersection,Real& intersection_di,
        Real& intersection_dj,Real& intersection_dk);

void compute_intersections_3rd(
        const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        const Eigen::Transform<Real,3,Eigen::Affine>& bwd_transform,
        const Eigen::Transform<Real,3,Eigen::Affine>& fwd_transform,
        uint dimension,const uint8_t& refLevel,
        Real& intersection,Real& intersection_di,
        Real& intersection_dj,Real& intersection_dk);

#endif
