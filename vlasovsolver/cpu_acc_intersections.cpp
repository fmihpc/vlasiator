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

#include <algorithm>
#include <cmath>
#include <utility>

/*TODO - replace with standard library c++11 functions as soon as possible*/
//#include <boost/array.hpp>
//#include <boost/unordered_map.hpp>
//#include <boost/unordered_set.hpp>

#include <Eigen/Geometry>

#include "../common.h"
#include "../spatial_cell_wrapper.hpp"
#include "cpu_acc_intersections.hpp"

using namespace std;
using namespace Eigen;

/**
  Computes the intersection point of a plane and a line

  \param l_point Point on the line
  \param l_direction Vector in the direction of the line
  \param p_point Point on plane
  \param p_normal Normal vector to plane
  \param intersection The function will set this to the intersection point.*/
Eigen::Matrix<Real,3,1> line_plane_intersection(const Eigen::Matrix<Real,3,1>& l_point,const Eigen::Matrix<Real,3,1>& l_direction,
                                                const Eigen::Matrix<Real,3,1>& p_point,const Eigen::Matrix<Real,3,1>& p_normal){
  const Real nom=p_normal.dot(p_point-l_point);
  const Real dem=p_normal.dot(l_direction);
  return l_point+(nom/dem)*l_direction;
}

/*!
Computes the first intersection data; this is z~ in section 2.4 in Zerroukat et al (2012). We assume all velocity cells have the same dimensions. 
Intersection z coordinate for (i,j,k) is: intersection + i * intersection_di + j * intersection_dj + k * intersection_dk 
\param spatial_cell spatial cell that is accelerated
\param fwd_transform Transform that describes acceleration forward in time
\param bwd_transform Transform that describes acceleration backward in time, used to compute the lagrangian departure gri
\param dimension Along which dimension is this intersection/mapping computation done. It is assumed the three mappings are in order 012, 120 or 201 
\param refLevel Refinement level at which intersections are computed
\param intersection Intersection z coordinate at i,j,k=0
\param intersection_di Change in z-coordinate for a change in i index of 1
\param intersection_dj Change in z-coordinate for a change in j index of 1
\param intersection_dk Change in z-coordinate for a change in k index of 1 */
void compute_intersections_1st(
        const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        const Transform<Real,3,Affine>& bwd_transform,const Transform<Real,3,Affine>& fwd_transform,
        uint dimension,const uint8_t& refLevel,
        Real& intersection,Real& intersection_di,Real& intersection_dj,Real& intersection_dk) {
    
    if (dimension == 0) { //Prepare intersections for mapping along X first (mapping order X-Y-Z)
        // Normal of Lagrangian planes
        const Eigen::Matrix<Real,3,1> plane_normal = bwd_transform.linear()*Eigen::Matrix<Real,3,1>(1.0, 0.0, 0.0);
        // Point on lowest potential Lagrangian plane
        const Eigen::Matrix<Real,3,1> plane_point
            = bwd_transform*Eigen::Matrix<Real,3,1>(vmesh.getMeshMinLimits()[0], 0.0, 0.0);
        // line along Euclidian x direction, unit vector
        const Eigen::Matrix<Real,3,1> line_direction = Eigen::Matrix<Real,3,1>(1.0, 0.0, 0.0);
        const Eigen::Matrix<Real,3,1> line_point(
            0.0,
            0.5*vmesh.getCellSize(refLevel)[1]+vmesh.getMeshMinLimits()[1],
            0.5*vmesh.getCellSize(refLevel)[2]+vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> lagrangian_di = bwd_transform.linear()
            * Eigen::Matrix<Real,3,1>(vmesh.getCellSize(refLevel)[0],0,0.0);
        const Eigen::Matrix<Real,3,1> euclidian_dj = Eigen::Matrix<Real,3,1>(0,vmesh.getCellSize(refLevel)[1],0.0);
        const Eigen::Matrix<Real,3,1> euclidian_dk = Eigen::Matrix<Real,3,1>(0.0,0.0,vmesh.getCellSize(refLevel)[2]);
          
        // compute intersections, varying lines and plane in i,j,k
        const Eigen::Matrix<Real,3,1> intersection_0_0_0 = line_plane_intersection(line_point, line_direction, plane_point, plane_normal);
        const Eigen::Matrix<Real,3,1> intersection_1_0_0 = line_plane_intersection(line_point, line_direction, plane_point + lagrangian_di, plane_normal);
        const Eigen::Matrix<Real,3,1> intersection_0_1_0 = line_plane_intersection(line_point + euclidian_dj, line_direction, plane_point, plane_normal);
        const Eigen::Matrix<Real,3,1> intersection_0_0_1 = line_plane_intersection(line_point + euclidian_dk, line_direction, plane_point, plane_normal);
        intersection=intersection_0_0_0[dimension];
        intersection_di=intersection_1_0_0[dimension]-intersection_0_0_0[dimension];
        intersection_dj=intersection_0_1_0[dimension]-intersection_0_0_0[dimension];
        intersection_dk=intersection_0_0_1[dimension]-intersection_0_0_0[dimension];
    }
    if (dimension == 1) { //Prepare intersections for mapping along Y first (mapping order Y-Z-X)
        // Normal of Lagrangian planes
        const Eigen::Matrix<Real,3,1> plane_normal = bwd_transform.linear()*Eigen::Matrix<Real,3,1>(0.0, 1.0, 0.0); 
        // Point on lowest possible Lagrangian plane
        const Eigen::Matrix<Real,3,1> plane_point
            = bwd_transform*Eigen::Matrix<Real,3,1>(0.0, vmesh.getMeshMinLimits()[1], 0.0); 
        // line along Euclidian y direction, unit vector
        const Eigen::Matrix<Real,3,1> line_direction = Eigen::Matrix<Real,3,1>(0.0, 1.0, 0.0);
        const Eigen::Matrix<Real,3,1> line_point(
            0.5*vmesh.getCellSize(refLevel)[0]+vmesh.getMeshMinLimits()[0],
            0.0,
            0.5*vmesh.getCellSize(refLevel)[2]+vmesh.getMeshMinLimits()[2]);  
        const Eigen::Matrix<Real,3,1> euclidian_di
            = Eigen::Matrix<Real,3,1>(vmesh.getCellSize(refLevel)[0], 0.0, 0.0); 
        const Eigen::Matrix<Real,3,1> lagrangian_dj 
            = bwd_transform.linear()*Eigen::Matrix<Real,3,1>(0.0 ,vmesh.getCellSize(refLevel)[1], 0.0); 
        const Eigen::Matrix<Real,3,1> euclidian_dk
            = Eigen::Matrix<Real,3,1>(0.0 , 0.0 ,vmesh.getCellSize(refLevel)[2]); 
        
        // compute intersections, varying lines and plane in i,j,k
        const Eigen::Matrix<Real,3,1> intersection_0_0_0 = line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
        const Eigen::Matrix<Real,3,1> intersection_1_0_0 = line_plane_intersection(line_point + euclidian_di, line_direction, plane_point, plane_normal);
        const Eigen::Matrix<Real,3,1> intersection_0_1_0 = line_plane_intersection(line_point, line_direction, plane_point + lagrangian_dj, plane_normal);
        const Eigen::Matrix<Real,3,1> intersection_0_0_1 = line_plane_intersection(line_point + euclidian_dk, line_direction, plane_point, plane_normal);

        intersection=intersection_0_0_0[dimension];
        intersection_di=intersection_1_0_0[dimension]-intersection_0_0_0[dimension];
        intersection_dj=intersection_0_1_0[dimension]-intersection_0_0_0[dimension];
        intersection_dk=intersection_0_0_1[dimension]-intersection_0_0_0[dimension];
    }

    if (dimension == 2) {
        // This is the  case presented in the Slice 3D article
        // Prepare intersections for mapping along Z first (mapping order Z-X-Y)
        
        //Normal of Lagrangian planes
        const Eigen::Matrix<Real,3,1> plane_normal
            = bwd_transform.linear()*Eigen::Matrix<Real,3,1>(0,0,1.0); 
      
        // Point on lowest possible Lagrangian plane
        const Eigen::Matrix<Real,3,1> plane_point
            = bwd_transform*Eigen::Matrix<Real,3,1>(0.0,0.0,vmesh.getMeshMinLimits()[2]);

        // line along Euclidian z direction, unit vector
        const Eigen::Matrix<Real,3,1> line_direction = Eigen::Matrix<Real,3,1>(0,0,1.0); 
        const Eigen::Matrix<Real,3,1> line_point(
            0.5*vmesh.getCellSize(refLevel)[0]+vmesh.getMeshMinLimits()[0],
            0.5*vmesh.getCellSize(refLevel)[1]+vmesh.getMeshMinLimits()[1],
            0.0);
        const Eigen::Matrix<Real,3,1> euclidian_di
            = Eigen::Matrix<Real,3,1>(vmesh.getCellSize(refLevel)[0],0,0.0);
        const Eigen::Matrix<Real,3,1> euclidian_dj
            = Eigen::Matrix<Real,3,1>(0,vmesh.getCellSize(refLevel)[1],0.0); 
        const Eigen::Matrix<Real,3,1> lagrangian_dk
            = bwd_transform.linear()*Eigen::Matrix<Real,3,1>(0.0,0.0,vmesh.getCellSize(refLevel)[2]);

      // compute intersections, varying lines and plane in i,j,k
      const Eigen::Matrix<Real,3,1> intersection_0_0_0 = line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
      const Eigen::Matrix<Real,3,1> intersection_1_0_0 = line_plane_intersection(line_point + euclidian_di, line_direction, plane_point, plane_normal);
      const Eigen::Matrix<Real,3,1> intersection_0_1_0 = line_plane_intersection(line_point + euclidian_dj, line_direction, plane_point, plane_normal);
      const Eigen::Matrix<Real,3,1> intersection_0_0_1 = line_plane_intersection(line_point, line_direction, plane_point + lagrangian_dk, plane_normal);
      intersection=intersection_0_0_0[dimension];
      intersection_di=intersection_1_0_0[dimension]-intersection_0_0_0[dimension];
      intersection_dj=intersection_0_1_0[dimension]-intersection_0_0_0[dimension];
      intersection_dk=intersection_0_0_1[dimension]-intersection_0_0_0[dimension];
   }
}

/*!
  Computes the second intersection data; this is x~ in section 2.4 in Zerroukat et al (2012). We assume all velocity cells have the same dimensions.
  Intersection x coordinate for (i,j,k) is: intersection + i * intersection_di + j * intersection_dj + k * intersection_dk 
  \param spatial_cell spatial cell that is accelerated
  \param fwd_transform Transform that describes acceleration forward in time
  \param bwd_transform Transform that describes acceleration backward in time, used to compute the lagrangian departure grid
  \param dimension Along which dimension is this intersection/mapping computation done. It is assumed the three mappings are in order 012, 120 or 201 
  \param refLevel Refinement level at which the intersections are computed
  \param intersection Intersection x-coordinate at i,j,k=0
  \param intersection_di Change in x-coordinate for a change in i index of 1
  \param intersection_dj Change in x-coordinate for a change in j index of 1
  \param intersection_dk Change in x-coordinate for a change in k index of 1*/
void compute_intersections_2nd(
        const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        const Transform<Real,3,Affine>& bwd_transform,const Transform<Real,3,Affine>& fwd_transform,
        uint dimension,const uint8_t& refLevel,
        Real& intersection,Real& intersection_di,Real& intersection_dj,Real& intersection_dk){
       
    if (dimension == 0) { // Prepare intersections for mapping along X second (mapping order Z-X-Y)
        // This is the case presented in the Slice 3D article, 
        // data along z has been moved to Lagrangian coordinates.
        
        // Normal of Euclidian y-plane
        const Eigen::Matrix<Real,3,1> plane_normal = Eigen::Matrix<Real,3,1>(0.0, 1.0, 0.0);
        
        //Point on lowest Euclidian y-plane through middle of cells
        Eigen::Matrix<Real,3,1> plane_point
            = Eigen::Matrix<Real,3,1>(0,vmesh.getMeshMinLimits()[1]+vmesh.getCellSize(refLevel)[1]*0.5,0); 
        const Eigen::Matrix<Real,3,1> lagrangian_di
            = bwd_transform.linear()*Eigen::Matrix<Real,3,1>(vmesh.getCellSize(refLevel)[0],0,0.0);
        
        // Distance between Euclidian planes
        const Eigen::Matrix<Real,3,1> euclidian_dj
            = Eigen::Matrix<Real,3,1>(0,vmesh.getCellSize(refLevel)[1],0.0);
        const Eigen::Matrix<Real,3,1> lagrangian_dk
            = bwd_transform.linear()*Eigen::Matrix<Real,3,1>(0.0,0.0,vmesh.getCellSize(refLevel)[2]);
      
        // line along Lagrangian y line, unit vector. Only rotation here, not translation
        const Eigen::Matrix<Real,3,1> line_direction = bwd_transform.linear() * Eigen::Matrix<Real,3,1>(0,1.0,0.0); 
        const Eigen::Matrix<Real,3,1> line_point = bwd_transform * Eigen::Matrix<Real,3,1>(
            vmesh.getMeshMinLimits()[0],
            0.5*vmesh.getCellSize(refLevel)[1]+vmesh.getMeshMinLimits()[1],
            0.5*vmesh.getCellSize(refLevel)[2]+vmesh.getMeshMinLimits()[2]);

        // Compute two intersections between Lagrangian line (absolute position 
        // does not matter so set to 0,0,0, and two Euclidian planes.
        Eigen::Matrix<Real,3,1> intersect_0_0_0 = line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
        Eigen::Matrix<Real,3,1> intersect_1_0_0 = line_plane_intersection(line_point + lagrangian_di, line_direction, plane_point, plane_normal);
        Eigen::Matrix<Real,3,1> intersect_0_1_0 = line_plane_intersection(line_point, line_direction, plane_point + euclidian_dj, plane_normal);
        Eigen::Matrix<Real,3,1> intersect_0_0_1 = line_plane_intersection(line_point + lagrangian_dk, line_direction, plane_point, plane_normal);

        intersection=intersect_0_0_0[dimension];
        intersection_di = intersect_1_0_0[dimension] - intersect_0_0_0[dimension];
        intersection_dj = intersect_0_1_0[dimension] - intersect_0_0_0[dimension];
        intersection_dk = intersect_0_0_1[dimension] - intersect_0_0_0[dimension];
    }
    if (dimension == 1) { //Prepare intersections for mapping along Y second (mapping order X-Y-Z)
        // Normal of Euclidian z-plane
        const Eigen::Matrix<Real,3,1> plane_normal = Eigen::Matrix<Real,3,1>(0.0, 0.0, 1.0);
        
        // Point on lowest Euclidian z-plane through middle of cells
        Eigen::Matrix<Real,3,1> plane_point
            = Eigen::Matrix<Real,3,1>(0.0, 0.0,vmesh.getMeshMinLimits()[2]+vmesh.getCellSize(refLevel)[2] * 0.5);

        const Eigen::Matrix<Real,3,1> lagrangian_di
            = bwd_transform.linear() * Eigen::Matrix<Real,3,1>(vmesh.getCellSize(refLevel)[0], 0.0,  0.0); 
        const Eigen::Matrix<Real,3,1> lagrangian_dj
            = bwd_transform.linear() * Eigen::Matrix<Real,3,1>(0.0, vmesh.getCellSize(refLevel)[1], 0.0);
        // Distance between Euclidian planes
        const Eigen::Matrix<Real,3,1> euclidian_dk
            = Eigen::Matrix<Real,3,1>(0.0, 0.0, vmesh.getCellSize(refLevel)[2]); 
  
        // line along Lagrangian z line, unit vector. Only rotation here, not translation
        const Eigen::Matrix<Real,3,1> line_direction 
            = bwd_transform.linear() * Eigen::Matrix<Real,3,1>(0.0, 0.0, 1.0); 
        const Eigen::Matrix<Real,3,1> line_point = bwd_transform * Eigen::Matrix<Real,3,1>(
            0.5*vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
            vmesh.getMeshMinLimits()[1],
            0.5*vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        
        // Compute two intersections between Lagrangian line (absolute position 
        // does not matter so set to 0,0,0, and two Euclidian planes.
        Eigen::Matrix<Real,3,1> intersect_0_0_0 = line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
        Eigen::Matrix<Real,3,1> intersect_1_0_0 = line_plane_intersection(line_point + lagrangian_di, line_direction, plane_point, plane_normal);
        Eigen::Matrix<Real,3,1> intersect_0_1_0 = line_plane_intersection(line_point + lagrangian_dj, line_direction, plane_point, plane_normal);
        Eigen::Matrix<Real,3,1> intersect_0_0_1 = line_plane_intersection(line_point, line_direction, plane_point + euclidian_dk, plane_normal);

        intersection=intersect_0_0_0[dimension];
        intersection_di = intersect_1_0_0[dimension] - intersect_0_0_0[dimension];
        intersection_dj = intersect_0_1_0[dimension] - intersect_0_0_0[dimension];
        intersection_dk = intersect_0_0_1[dimension] - intersect_0_0_0[dimension];
   
    }    
    if (dimension == 2) { //Prepare intersections for mapping along Z second (mapping order Y-Z-X)
        // Normal of Euclidian x-plane
        const Eigen::Matrix<Real,3,1> plane_normal = Eigen::Matrix<Real,3,1>(1.0, 0.0, 0.0);
        //Point on lowest Euclidian x-plane through middle of cells
        Eigen::Matrix<Real,3,1> plane_point 
            = Eigen::Matrix<Real,3,1>(vmesh.getMeshMinLimits()[0]+vmesh.getCellSize(refLevel)[0]*0.5, 0.0, 0.0);
        // Distance between Euclidian planes
        const Eigen::Matrix<Real,3,1> euclidian_di = Eigen::Matrix<Real,3,1>(vmesh.getCellSize(refLevel)[0], 0.0, 0.0);
        const Eigen::Matrix<Real,3,1> lagrangian_dj 
            = bwd_transform.linear() * Eigen::Matrix<Real,3,1>(0.0, vmesh.getCellSize(refLevel)[1], 0.0); 
        const Eigen::Matrix<Real,3,1> lagrangian_dk = bwd_transform.linear() * Eigen::Matrix<Real,3,1>(0.0, 0.0, vmesh.getCellSize(refLevel)[2]); 
  
        // line along Lagrangian x line, unit vector. Only rotation here, not translation
        const Eigen::Matrix<Real,3,1> line_direction = bwd_transform.linear() * Eigen::Matrix<Real,3,1>(1.0, 0.0, 0.0);
        const Eigen::Matrix<Real,3,1> line_point = bwd_transform * Eigen::Matrix<Real,3,1>(
            0.5 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
            0.5 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
            vmesh.getMeshMinLimits()[2]);  
        
        // Compute two intersections between Lagrangian line (absolute position 
        // does not matter so set to 0,0,0, and two Euclidian planes.
        Eigen::Matrix<Real,3,1> intersect_0_0_0 = line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
        Eigen::Matrix<Real,3,1> intersect_1_0_0 = line_plane_intersection(line_point, line_direction, plane_point + euclidian_di, plane_normal);
        Eigen::Matrix<Real,3,1> intersect_0_1_0 = line_plane_intersection(line_point + lagrangian_dj, line_direction, plane_point, plane_normal);
        Eigen::Matrix<Real,3,1> intersect_0_0_1 = line_plane_intersection(line_point + lagrangian_dk, line_direction, plane_point, plane_normal);

        intersection=intersect_0_0_0[dimension];
        intersection_di = intersect_1_0_0[dimension] - intersect_0_0_0[dimension];
        intersection_dj = intersect_0_1_0[dimension] - intersect_0_0_0[dimension];
        intersection_dk = intersect_0_0_1[dimension] - intersect_0_0_0[dimension];
    }
}

/*!
  Computes the third intersection data; this is y intersections in Zerroukat et al (2012). We assume all velocity cells have the same dimensions.
  Intersection y-coordinate for (i,j,k) is: intersection + i * intersection_di + j * intersection_dj + k * intersection_dk 
  \param spatial_cell spatial cell that is accelerated
  \param fwd_transform Transform that describes acceleration forward in time
  \param bwd_transform Transform that describes acceleration backward in time, used to compute the lagrangian departure grid
  \param dimension Along which dimension is this intersection/mapping computation done. It is assumed the three mappings are in order 012, 120 or 201 
  \param refLevel Refinement level at which intersections are computed
  \param intersection Intersection y-coordinate at i,j,k=0
  \param intersection_di Change in y-coordinate for a change in i index of 1
  \param intersection_dj Change in y-coordinate for a change in j index of 1
  \param intersection_dk Change in y-coordinate for a change in k index of 1


  euclidian y goes from vy_min to vy_max, this is mapped to wherever y plane is in lagrangian.*/
void compute_intersections_3rd(
    const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
    const Transform<Real,3,Affine>& bwd_transform,const Transform<Real,3,Affine>& fwd_transform,
    uint dimension,const uint8_t& refLevel,
    Real& intersection,Real& intersection_di,Real& intersection_dj,Real& intersection_dk) {
    
    if (dimension == 0) { //Prepare intersections for mapping along X third (mapping order Y-Z-X)
        const Eigen::Matrix<Real,3,1> point_0_0_0 = bwd_transform
            * Eigen::Matrix<Real,3,1>(0.0 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      0.5 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
                                      0.5 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> point_1_0_0 = bwd_transform
            * Eigen::Matrix<Real,3,1>(1.0 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      0.5 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
                                      0.5 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> point_0_1_0 = bwd_transform
            * Eigen::Matrix<Real,3,1>(0.0 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      1.5 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
                                      0.5 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> point_0_0_1 = bwd_transform
            * Eigen::Matrix<Real,3,1>(0.0 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      0.5 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
                                      1.5 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        intersection = point_0_0_0[dimension];
        intersection_di = point_1_0_0[dimension]-point_0_0_0[dimension];
        intersection_dj = point_0_1_0[dimension]-point_0_0_0[dimension];
        intersection_dk = point_0_0_1[dimension]-point_0_0_0[dimension];
   }
   if (dimension == 1) { //Prepare intersections for mapping along Y third (mapping order Z-X-Y)
        // This is the case presented in the Slice 3D article, 
        // data along z has been moved to Lagrangian coordinates
        const Eigen::Matrix<Real,3,1> point_0_0_0 = bwd_transform
            * Eigen::Matrix<Real,3,1>(0.5 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      vmesh.getMeshMinLimits()[1],
                                      0.5 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> point_1_0_0 = bwd_transform
            * Eigen::Matrix<Real,3,1>(1.5 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      vmesh.getMeshMinLimits()[1],
                                      0.5 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> point_0_1_0 = bwd_transform
            * Eigen::Matrix<Real,3,1>(0.5 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      1.0 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
                                      0.5 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> point_0_0_1 = bwd_transform
            * Eigen::Matrix<Real,3,1>(0.5 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      vmesh.getMeshMinLimits()[1],
                                      1.5 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        intersection = point_0_0_0[dimension];
        intersection_di = point_1_0_0[dimension]-point_0_0_0[dimension];
        intersection_dj = point_0_1_0[dimension]-point_0_0_0[dimension];
        intersection_dk = point_0_0_1[dimension]-point_0_0_0[dimension];
    }
    if (dimension == 2) { //Prepare intersections for mapping along Z third (mapping order X-Y-Z)
        const Eigen::Matrix<Real,3,1> point_0_0_0 = bwd_transform
            * Eigen::Matrix<Real,3,1>(0.5 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      0.5 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
                                      0.0 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> point_1_0_0 = bwd_transform
            * Eigen::Matrix<Real,3,1>(1.5 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      0.5 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
                                      0.0 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> point_0_1_0 = bwd_transform
            * Eigen::Matrix<Real,3,1>(0.5 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      1.5 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
                                      0.0 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        const Eigen::Matrix<Real,3,1> point_0_0_1 = bwd_transform
            * Eigen::Matrix<Real,3,1>(0.5 * vmesh.getCellSize(refLevel)[0] + vmesh.getMeshMinLimits()[0],
                                      0.5 * vmesh.getCellSize(refLevel)[1] + vmesh.getMeshMinLimits()[1],
                                      1.0 * vmesh.getCellSize(refLevel)[2] + vmesh.getMeshMinLimits()[2]);
        intersection = point_0_0_0[dimension];
        intersection_di = point_1_0_0[dimension]-point_0_0_0[dimension];
        intersection_dj = point_0_1_0[dimension]-point_0_0_0[dimension];
        intersection_dk = point_0_0_1[dimension]-point_0_0_0[dimension];
    }
}
