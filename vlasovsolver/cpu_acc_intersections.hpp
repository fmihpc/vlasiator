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
Computes the first intersection data; this is z~ in section 2.4 in Zerroukat et al (2012). We assume all velocity cells have the same dimensions

\param bwd_transform Transform that is used to compute the lagrangian departure grid
\param intersections map keys are indices of line (with regards to coord1 and coord2), values are value of coordinate at first intersection, number of intersections.
\param intersection_distance Distance between intersections
*/


void compute_intersections_z(SpatialCell* spatial_cell, std::vector<unsigned int> downstream_blocks,
		     const Transform<Real,3,Affine>& bwd_transform,const Transform<Real,3,Affine>& fwd_transform,
		     boost::unordered_map< boost::array<int,2> , std::vector<Real> >& intersections, 
		     Real& intersection_distance){
  
  cout <<"number of downstream blocks " <<downstream_blocks.size()<< "total blocks "<<spatial_cell->number_of_blocks <<endl;
  Real max_z=spatial_cell->at(downstream_blocks[0])->parameters[BlockParams::VZCRD];
  Real min_z=max_z;
  Real dvz=spatial_cell->cell_dvz;
   
  /*compute the i,j indices of all blocks, use set to get rid of dupicates
    Also compute max and min z of blocks*/
   boost::unordered_set< boost::array<int,2> > block_ij;
   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
     const unsigned int block = spatial_cell->velocity_block_list[block_i];
     Velocity_Block* block_ptr = spatial_cell->at(block);
     const Real block_start_vz=block_ptr->parameters[BlockParams::VZCRD];
     if(block_start_vz<min_z) min_z= block_start_vz;
     if(block_start_vz+WID*dvz>max_z) max_z=block_start_vz+WID*dvz;
     block_ij.insert( {{ block%WID , (block/WID)%WID }} );
   }

   const Eigen::Matrix<Real,3,1> plane_normal=bwd_transform*Eigen::Matrix<Real,3,1>(0,0,1.0); //Normal of lagrangian planes
   const Eigen::Matrix<Real,3,1> plane_point=bwd_transform*Eigen::Matrix<Real,3,1>(0,0,min_z); //Point on lowest lagrangian plane (not all blocks exist of course)
   const Eigen::Matrix<Real,3,1> plane_delta=bwd_transform*Eigen::Matrix<Real,3,1>(0,0,dvz); //vector between two lagrangian planes
   const Eigen::Matrix<Real,3,1> line_direction=Eigen::Matrix<Real,3,1>(0,0,1.0); //line along euclidian z direction, unit vector
   
   cout << "Angle of rotation "<< acos( plane_normal.dot(Eigen::Matrix<Real,3,1>(0,0,1.0))) * 180.0/(M_PI) <<endl;
   cout << "dvz "<< dvz;
   /*Distance between lagrangian planes along line direction in Euclidian coordinates,assuming line direction has unit length*/
   intersection_distance = plane_delta.dot(plane_delta)/plane_delta.dot(line_direction); 
   
   /*loop over instersections and add intersection position to map*/
   for (const auto& ij: block_ij){
    for (uint cell_xi=0;cell_xi<WID;cell_xi++){
      for (uint cell_yi=0;cell_yi<WID;cell_yi++){

	const Eigen::Matrix<Real,3,1> line_point((ij[0]*WID+cell_xi+0.5)*SpatialCell::block_dvx+SpatialCell::vx_min,
						 (ij[1]*WID+cell_yi+0.5)*SpatialCell::block_dvy+SpatialCell::vy_min,
						 0.0);
	const boost::array<int,2> cell_ij = {{ ij[0] * WID + cell_xi , ij[1] * WID + cell_yi }}; /* cell indices in x and y*/
	//	intersections.insert(std::make_pair(cell_ij,std::vector<Real>()));
	
	/*compute intersection between the line along z, with the lowest possible plane*/
	Eigen::Matrix<Real,3,1> intersection=line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
	
	/*loop over z coordinate of all possible intersections*/
	while(intersection[2]<max_z) {

	  /* Transform velocity coordinate of lagrangian cell to the actual one it has in the propagated state*/
	  const Eigen::Matrix<Real,3,1> lagrangian_cell_velocity=fwd_transform*(intersection + plane_delta);
	  
	  /*check if block exists*/
	  const unsigned int block=spatial_cell->get_velocity_block(lagrangian_cell_velocity[0],
								    lagrangian_cell_velocity[1],
								    lagrangian_cell_velocity[2]);

	  /* check if block exists, add z value to intersections list if it does*/
	  if(! spatial_cell->is_null_block(spatial_cell->at(block))){
	    intersections[cell_ij].push_back(intersection[2]);
	  }
	  /*go to next possible intersection*/
	  intersection[2]+=intersection_distance;
	}
      }
    }
   }
}

#endif
