#ifndef CPU_ACC_INTERSECTIONS_H
#define CPU_ACC_INTERSECTIONS_H

#include <Eigen/Core>

#include "../definitions.h"
#include "../spatial_cell.hpp"

Eigen::Matrix<Real,3,1> line_plane_intersection(const Eigen::Matrix<Real,3,1>& l_point,
                                                const Eigen::Matrix<Real,3,1>& l_direction,
                                                const Eigen::Matrix<Real,3,1>& p_point,
                                                const Eigen::Matrix<Real,3,1>& p_normal);

void compute_intersections_1st(
        const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        const Eigen::Transform<Real,3,Eigen::Affine>& bwd_transform,
        const Eigen::Transform<Real,3,Eigen::Affine>& fwd_transform,
        uint dimension,
        Real& intersection,Real& intersection_di,
        Real& intersection_dj,Real& intersection_dk);

void compute_intersections_2nd(
        const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        const Eigen::Transform<Real,3,Eigen::Affine>& bwd_transform,
        const Eigen::Transform<Real,3,Eigen::Affine>& fwd_transform,
        uint dimension,
        Real& intersection,Real& intersection_di,
        Real& intersection_dj,Real& intersection_dk);

void compute_intersections_3rd(
        const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        const Eigen::Transform<Real,3,Eigen::Affine>& bwd_transform,
        const Eigen::Transform<Real,3,Eigen::Affine>& fwd_transform,
        uint dimension,
        Real& intersection,Real& intersection_di,
        Real& intersection_dj,Real& intersection_dk);

#endif
