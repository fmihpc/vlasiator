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


void compute_intersections_z(SpatialCell* spatial_cell,
			     const Transform<Real,3,Affine>& bwd_transform,const Transform<Real,3,Affine>& fwd_transform,
			     boost::unordered_map< boost::array<int,3> , Real >& intersections, 
			     Real& intersection_distance){
  
  if(spatial_cell->number_of_blocks==0)
    return; //not much to do
  
  uint min_z_bindex,max_z_bindex
;  /*gather unique list of the i,j indices of all blocks, use set to get rid of dupicates
    Also compute max and min z of blocks*/
  boost::unordered_set< boost::array<int,2> > block_ij;
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
     const unsigned int block = spatial_cell->velocity_block_list[block_i];
     Velocity_Block* block_ptr = spatial_cell->at(block);
     const velocity_block_indices_t block_ijk= spatial_cell->get_velocity_block_indices(block);

     if(block_i==0) {
       /*initialize values on first iteration, will not cost much as
	 the branch predictor will optimize this if away*/
       min_z_bindex = block_ijk[2]; 
       max_z_bindex = block_ijk[2]; 
     }
     if( min_z_bindex > block_ijk[2] )
       min_z_bindex = block_ijk[2]; //this is the lowest block index in z-direction

     if( max_z_bindex < block_ijk[2] )
       max_z_bindex = block_ijk[2]; //this is the highest block index in z-direction

     block_ij.insert( {{ block_ijk[0] , block_ijk[1] }} );
  }
  
  cout << "min/mz v index "<< min_z_bindex<< " " <<max_z_bindex <<endl;
  const Eigen::Matrix<Real,3,1> plane_normal = bwd_transform*Eigen::Matrix<Real,3,1>(0,0,1.0); //Normal of lagrangian planes
  const Eigen::Matrix<Real,3,1> plane_point =  bwd_transform*Eigen::Matrix<Real,3,1>(0,0,SpatialCell::vz_min+SpatialCell::block_dvz*min_z_bindex); //Point on lowest lagrangian plane (not all blocks exist of course)
  const Eigen::Matrix<Real,3,1> plane_delta =  bwd_transform*Eigen::Matrix<Real,3,1>(0,0,SpatialCell::cell_dvz); //vector between two lagrangian planes
  const Eigen::Matrix<Real,3,1> line_direction = Eigen::Matrix<Real,3,1>(0,0,1.0); //line along euclidian z direction, unit vector
  
  cout << "Angle of rotation "<< acos( plane_normal.dot(Eigen::Matrix<Real,3,1>(0,0,1.0))) * 180.0/(M_PI) <<endl;

  /*Distance between lagrangian planes along line direction in Euclidian coordinates,assuming line direction has unit length*/
  intersection_distance = plane_delta.dot(plane_delta)/plane_delta.dot(line_direction); 
  
  /*loop over instersections and add intersection position to map*/
  /*first loop over ij in euclidian space*/
  for (const auto& ij: block_ij){    
    for (uint cell_xi=0;cell_xi<WID;cell_xi++){
      for (uint cell_yi=0;cell_yi<WID;cell_yi++){
	
	const Eigen::Matrix<Real,3,1> line_point((ij[0]*WID+cell_xi+0.5)*SpatialCell::cell_dvx+SpatialCell::vx_min,
						 (ij[1]*WID+cell_yi+0.5)*SpatialCell::cell_dvy+SpatialCell::vy_min,
						 0.0);
	boost::array<int,3> cell_ijk = {{ ij[0] * WID + cell_xi , ij[1] * WID + cell_yi , 0 }}; /* cell indices in x and y, z has a dummy value*/
	/*compute intersection between the line along z, with the lowest possible plane corresponding to z-block min_z_bindex*/
	Eigen::Matrix<Real,3,1> intersection=line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
	
	
	/*loop over z of all possible intersections, indices in lagrangian space*/
	for(uint  block_zi= min_z_bindex;block_zi<=max_z_bindex;block_zi++){
	  for(uint cell_zi=0;cell_zi<WID;cell_zi++) {
	    /* Transform velocity coordinate of lagrangian cell to the actual one it has in the propagated state*/
	    const Eigen::Matrix<Real,3,1> lagrangian_cell_velocity=fwd_transform * ( intersection + 0.5*plane_delta );
	    
	    /*get block corresponding to velocity*/
	    const unsigned int block=spatial_cell->get_velocity_block(lagrangian_cell_velocity[0],
								      lagrangian_cell_velocity[1],
								      lagrangian_cell_velocity[2]);
	    
	    /* check if block exists, add z value to intersections list if it does*/
	    /*TODO, this check is probably not foolproof, would need an additional cell to make sure nothing is missed?*/
	    if(! spatial_cell->is_null_block(spatial_cell->at(block))){
	      cell_ijk[2]=block_zi*WID+cell_zi;
	      intersections[cell_ijk]=intersection[2];
	    }
	    /*go to next possible intersection*/
	    intersection[2]+=intersection_distance;
	  }
	}
      }
    }
  }
}










/*!
Computes the second intersection data; this is x~ in section 2.4 in Zerroukat et al (2012). We assume all velocity cells have the same dimensions

\param bwd_transform Transform that is used to compute the lagrangian departure grid
\param intersections map keys are indices of the new cells to which mass in mapped, x and y in euclidian space, z in lagrangian. Value is z-coordinate of cell
\param intersection_distance Distance between intersections
*/


void compute_intersections_x(SpatialCell* spatial_cell,
			     const Transform<Real,3,Affine>& bwd_transform,const Transform<Real,3,Affine>& fwd_transform,
			     boost::unordered_map< boost::array<int,3> , Real >& intersections, 
			     Real& intersection_distance){
  
  if(spatial_cell->number_of_blocks==0)
    return; //not much to do
  
  /*TODO, in z-intersecttions we have essentially the same computation, should be combined...*/
  
  uint min_y_bindex,max_y_bindex;
  /*gather unique list of the i,k indices of all blocks, use set to get rid of dupicates
    Also compute max and min z of blocks*/
  boost::unordered_set< boost::array<int,2> > block_ik;  
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
    const unsigned int block = spatial_cell->velocity_block_list[block_i];
    Velocity_Block* block_ptr = spatial_cell->at(block);
    const velocity_block_indices_t block_ijk= spatial_cell->get_velocity_block_indices(block);
    
     if(block_i==0) {
       /*initialize values on first iteration, will not cost much as
	 the branch predictor will optimize this if away*/
       min_y_bindex = block_ijk[1]; 
       max_y_bindex = block_ijk[1]; 
     }
     if( min_y_bindex > block_ijk[1] )
       min_y_bindex = block_ijk[1]; //this is the lowest block index in z-direction
     
     if( max_y_bindex < block_ijk[1] )
       max_y_bindex = block_ijk[1]; //this is the highest block index in z-direction
     
     block_ik.insert( {{ block_ijk[0] , block_ijk[2] }} );
  }
  
  
  
  const Eigen::Matrix<Real,3,1> plane_normal = Eigen::Matrix<Real,3,1>(0,1.0,0.0); //Normal of Euclidian y-plane
  Eigen::Matrix<Real,3,1> plane_point =  Eigen::Matrix<Real,3,1>(0,SpatialCell::vy_min+SpatialCell::block_dvy*min_y_bindex+SpatialCell::cell_dvy*0.5,0); //Point on lowest euclidian y-plane through middle of cells (not all blocks exist of course)
  const Eigen::Matrix<Real,3,1> plane_delta = Eigen::Matrix<Real,3,1>(0,SpatialCell::cell_dvy,0.0); //Distance between euclidian planes
  const Eigen::Matrix<Real,3,1> line_direction = bwd_transform*Eigen::Matrix<Real,3,1>(0,1.0,0.0); //line along lagrangian y line, unit vector
  const Eigen::Matrix<Real,3,1> x_l = bwd_transform*Eigen::Matrix<Real,3,1>(1.0,0.0,0.0); //line along lagrangian x line, unit vector
  

  /*Distance between x-intersections of line_direction (lagrangian y) and the euclidian y-planes when going along euclidian y*/
  intersection_distance = SpatialCell::cell_dvx/x_l.dot(Eigen::Matrix<Real,3,1>(1.0,0.0,0.0));
  cout<< "x-dist" << intersection_distance << " dvx "<<SpatialCell::cell_dvx<<endl;

  /*Compute two intersections between lagrangian line (absolute position does not matter so set to 0,0,0, and two euclidian planes
  Eigen::Matrix<Real,3,1> intersect_1=line_plane_intersection(Eigen::Matrix<Real,3,1>(0,0,0),line_direction,plane_point,plane_normal);
  Eigen::Matrix<Real,3,1> intersect_2=line_plane_intersection(Eigen::Matrix<Real,3,1>(0,0,0),line_direction,plane_point+plane_delta,plane_normal);
  Real intersection_distance_L_x =intersection_2[0]-intersection_1[0];
  

  /*loop over instersections and add intersection position to map*/
  /*first loop over ij in euclidian space*/


  for (const auto& ik: block_ik){    
    for (uint cell_xi=0;cell_xi<WID;cell_xi++){
      for (uint cell_zi=0;cell_zi<WID;cell_zi++){
	
	const Eigen::Matrix<Real,3,1> line_point = bwd_transform*Eigen::Matrix<Real,3,1>((ik[0]*WID+cell_xi+0.5)*SpatialCell::cell_dvx+SpatialCell::vx_min,
										       0.0,
										       (ik[1]*WID+cell_zi+0.5)*SpatialCell::cell_dvz+SpatialCell::vz_min);
	
	boost::array<int,3> cell_ijk = {{ ik[0] * WID + cell_xi,0 , ik[1] * WID + cell_zi  }}; /* cell indices in x and z, y has a dummy value*/

	/*compute intersection between the line along lagrangian y, with the lowest possible plane corresponding to y-block min_y_bindex*/
	
	
	/*loop over y index of euclidian planes to get all intersection*/
	plane_point =  Eigen::Matrix<Real,3,1>(0,SpatialCell::vy_min+SpatialCell::block_dvy*min_y_bindex+SpatialCell::cell_dvy*0.5,0); //Point on lowest euclidian y-plane through middle of cells (not all blocks exist of course)
	for(uint  block_yi= min_y_bindex;block_yi<=max_y_bindex;block_yi++){
	  for(uint cell_yi=0;cell_yi<WID;cell_yi++) {
	    
	    Eigen::Matrix<Real,3,1> intersection=line_plane_intersection(line_point,line_direction,plane_point,plane_normal);
	    plane_point[1]+=SpatialCell::cell_dvy; //increase y coordinated
	    
	    /* Transform velocity coordinate of lagrangian cell to the actual one it has in the propagated state*/
	    const Eigen::Matrix<Real,3,1> lagrangian_cell_velocity=fwd_transform * ( intersection + Eigen::Matrix<Real,3,1>(0.01,0,0) ); //intersection + a small delta (0.01 is small compared to normal velocities)
	    
	    /*check if block exists*/
	    const unsigned int block=spatial_cell->get_velocity_block(lagrangian_cell_velocity[0],
								      lagrangian_cell_velocity[1],
								      lagrangian_cell_velocity[2]);
	    
	    /* check if block exists, add z value to intersections list if it does*/
	    if(! spatial_cell->is_null_block(spatial_cell->at(block))){
	      cell_ijk[1]=block_yi*WID+cell_yi;
	      intersections[cell_ijk]=intersection[0];
	    }

	  }
	}
      }
    }
  }
}



#endif
