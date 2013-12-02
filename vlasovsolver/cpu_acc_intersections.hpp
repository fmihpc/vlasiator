#ifndef CPU_ACC_INTERSECTIONS_H
#define CPU_ACC_INTERSECTIONS_H

#include "algorithm"
#include "cmath"
#include "utility"

/*TODO - replace with standard library c++11 functions as soon as possible*/
#include "boost/array.hpp"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "common.h"
#include "spatial_cell.hpp"
#include <Eigen/Geometry>
#include <Eigen/Core>

using namespace std;
using namespace spatial_cell;
using namespace Eigen;


/*!
  Computes the intersection point of a plane and a line

  \param l_point Point on the line
  \param l_direction Vector in the direction of the line
  \param p_point Point on plane
  \param p_normal Normal vector to plane
  \param intersection The function will set this to the intersection point
*/


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
\param intersection Intersection z coordinate at i,j,k=0
\param intersection_di Change in z-coordinate for a change in i index of 1
\param intersection_dj Change in z-coordinate for a change in j index of 1
\param intersection_dk Change in z-coordinate for a change in k index of 1
*/


void compute_intersections_z(const SpatialCell* spatial_cell,
			     const Transform<Real,3,Affine>& bwd_transform,const Transform<Real,3,Affine>& fwd_transform,
			     Real& intersection,Real& intersection_di,Real& intersection_dj,Real& intersection_dk){
  
  const Eigen::Matrix<Real,3,1> plane_normal = bwd_transform*Eigen::Matrix<Real,3,1>(0,0,1.0); //Normal of lagrangian planes
  const Eigen::Matrix<Real,3,1> plane_point =  bwd_transform*Eigen::Matrix<Real,3,1>(0.0,0.0,SpatialCell::vz_min); /*<Point on lowest potential lagrangian plane */
  const Eigen::Matrix<Real,3,1> plane_delta =  bwd_transform*Eigen::Matrix<Real,3,1>(0,0,SpatialCell::cell_dvz); //vector between two lagrangian planes
  const Eigen::Matrix<Real,3,1> line_direction = Eigen::Matrix<Real,3,1>(0,0,1.0); //line along euclidian z direction, unit vector
  const Eigen::Matrix<Real,3,1> line_point(0.5*SpatialCell::cell_dvx+SpatialCell::vx_min,
					   0.5*SpatialCell::cell_dvy+SpatialCell::vy_min,
					   0.0);

  const Eigen::Matrix<Real,3,1> euclidian_di  = Eigen::Matrix<Real,3,1>(SpatialCell::cell_dvx,0,0.0); 
  const Eigen::Matrix<Real,3,1> euclidian_dj  = Eigen::Matrix<Real,3,1>(0,SpatialCell::cell_dvy,0.0); 
  const Eigen::Matrix<Real,3,1> lagrangian_dk = bwd_transform*Eigen::Matrix<Real,3,1>(0.0,0.0,SpatialCell::cell_dvz); 
  
  
  
  /*compute intersections, varying lines and plane in i,j,k*/
  const Eigen::Matrix<Real,3,1> intersection_0_0_0 = line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
  const Eigen::Matrix<Real,3,1> intersection_1_0_0 = line_plane_intersection(line_point + euclidian_di, line_direction, plane_point, plane_normal);
  const Eigen::Matrix<Real,3,1> intersection_0_1_0 = line_plane_intersection(line_point + euclidian_dj, line_direction, plane_point, plane_normal);
  const Eigen::Matrix<Real,3,1> intersection_0_0_1 = line_plane_intersection(line_point, line_direction, plane_point + lagrangian_dk, plane_normal);

  intersection=intersection_0_0_0[2];
  intersection_di=intersection_1_0_0[2]-intersection_0_0_0[2];
  intersection_dj=intersection_0_1_0[2]-intersection_0_0_0[2];
  intersection_dk=intersection_0_0_1[2]-intersection_0_0_0[2];
}


/*!
  Computes the second intersection data; this is x~ in section 2.4 in Zerroukat et al (2012). We assume all velocity cells have the same dimensions.
  Intersection x coordinate for (i,j,k) is: intersection + i * intersection_di + j * intersection_dj + k * intersection_dk 
  \param spatial_cell spatial cell that is accelerated
  \param fwd_transform Transform that describes acceleration forward in time
  \param bwd_transform Transform that describes acceleration backward in time, used to compute the lagrangian departure grid
  \param intersection Intersection x-coordinate at i,j,k=0
  \param intersection_di Change in x-coordinate for a change in i index of 1
  \param intersection_dj Change in x-coordinate for a change in j index of 1
  \param intersection_dk Change in x-coordinate for a change in k index of 1
*/
void compute_intersections_x(const SpatialCell* spatial_cell,
			     const Transform<Real,3,Affine>& bwd_transform,const Transform<Real,3,Affine>& fwd_transform,
			     Real& intersection,Real& intersection_di,Real& intersection_dj,Real& intersection_dk){

  const Eigen::Matrix<Real,3,1> plane_normal = Eigen::Matrix<Real,3,1>(0.0, 1.0, 0.0); //Normal of Euclidian y-plane
  Eigen::Matrix<Real,3,1> plane_point =  Eigen::Matrix<Real,3,1>(0,SpatialCell::vy_min+SpatialCell::cell_dvy*0.5,0); //Point on lowest euclidian y-plane through middle of cells
  const Eigen::Matrix<Real,3,1> plane_delta = Eigen::Matrix<Real,3,1>(0,SpatialCell::cell_dvy,0.0); //Distance between euclidian planes

  const Eigen::Matrix<Real,3,1> lagrangian_di = bwd_transform * Eigen::Matrix<Real,3,1>(SpatialCell::cell_dvx,0,0.0); 
  const Eigen::Matrix<Real,3,1> euclidian_dj  = Eigen::Matrix<Real,3,1>(0,SpatialCell::cell_dvy,0.0); //Distance between euclidian planes
  const Eigen::Matrix<Real,3,1> lagrangian_dk = bwd_transform * Eigen::Matrix<Real,3,1>(0.0,0.0,SpatialCell::cell_dvz); 
  
  
  const Eigen::Matrix<Real,3,1> line_direction = bwd_transform * Eigen::Matrix<Real,3,1>(0,1.0,0.0); //line along lagrangian y line, unit vector
  const Eigen::Matrix<Real,3,1> line_point = bwd_transform * Eigen::Matrix<Real,3,1>(SpatialCell::vx_min,
										     0.5*SpatialCell::cell_dvy+SpatialCell::vy_min,
										     0.5*SpatialCell::cell_dvz+SpatialCell::vz_min);  
  
  /*Compute two intersections between lagrangian line (absolute position does not matter so set to 0,0,0, and two euclidian planes*/
  Eigen::Matrix<Real,3,1> intersect_0_0_0 = line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
  Eigen::Matrix<Real,3,1> intersect_1_0_0 = line_plane_intersection(line_point + lagrangian_di, line_direction, plane_point, plane_normal);
  Eigen::Matrix<Real,3,1> intersect_0_1_0 = line_plane_intersection(line_point, line_direction, plane_point + euclidian_dj, plane_normal);
  Eigen::Matrix<Real,3,1> intersect_0_0_1 = line_plane_intersection(line_point + lagrangian_dk, line_direction, plane_point, plane_normal);
  
  intersection=intersect_0_0_0[0];
  intersection_di = intersect_0_0_0[0] - intersect_1_0_0[0];
  intersection_dj = intersect_0_0_0[0] - intersect_0_1_0[0];
  intersection_dk = intersect_0_0_0[0] - intersect_0_0_1[0];

}

#endif
