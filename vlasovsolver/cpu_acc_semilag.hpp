/*
This file is part of Vlasiator.

Copyright 2012 Finnish Meteorological Institute

*/

#ifndef CPU_ACC_SEMILAG_H
#define CPU_ACC_SEMILAG_H

#include "algorithm"
#include "cmath"
#include "utility"

/*TODO - replace with standard library c++11 functions as soon as possible*/
#include "boost/array.hpp"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "boost/lexical_cast.hpp"

#include "common.h"
#include "spatial_cell.hpp"

#include <Eigen/Geometry>
#include "vlasovsolver/cpu_compute_downstream_blocks.hpp"
using namespace std;
using namespace spatial_cell;
using namespace Eigen;



/*Compute transform during on timestep, and update the bulk velocity of the cell*/

Transform<Real,3,Affine> compute_acceleration_transformation( SpatialCell* spatial_cell, const Real dt) {
   /*total field*/
   const Real Bx = spatial_cell->parameters[CellParams::BGBXVOL]+spatial_cell->parameters[CellParams::PERBXVOL];
   const Real By = spatial_cell->parameters[CellParams::BGBYVOL]+spatial_cell->parameters[CellParams::PERBYVOL];
   const Real Bz = spatial_cell->parameters[CellParams::BGBZVOL]+spatial_cell->parameters[CellParams::PERBZVOL];
   /*perturbed field*/
   const Real perBx = spatial_cell->parameters[CellParams::PERBXVOL];
   const Real perBy = spatial_cell->parameters[CellParams::PERBYVOL];
   const Real perBz = spatial_cell->parameters[CellParams::PERBZVOL];   
   //read in derivatives need for curl of B (only pertrubed, curl of background field is always 0!)
   const Real dBXdy = spatial_cell->derivativesBVOL[bvolderivatives::dPERBXVOLdy]/spatial_cell->parameters[CellParams::DY];
   const Real dBXdz = spatial_cell->derivativesBVOL[bvolderivatives::dPERBXVOLdz]/spatial_cell->parameters[CellParams::DZ];
   const Real dBYdx = spatial_cell->derivativesBVOL[bvolderivatives::dPERBYVOLdx]/spatial_cell->parameters[CellParams::DX];

   const Real dBYdz = spatial_cell->derivativesBVOL[bvolderivatives::dPERBYVOLdz]/spatial_cell->parameters[CellParams::DZ];
   const Real dBZdx = spatial_cell->derivativesBVOL[bvolderivatives::dPERBZVOLdx]/spatial_cell->parameters[CellParams::DX];
   const Real dBZdy = spatial_cell->derivativesBVOL[bvolderivatives::dPERBZVOLdy]/spatial_cell->parameters[CellParams::DY];

   
   const Eigen::Matrix<Real,3,1> B(Bx,By,Bz);
   const Eigen::Matrix<Real,3,1> unit_B(B.normalized());
   const Real gyro_period = 2 * M_PI * Parameters::m  / (fabs(Parameters::q) * B.norm());
   //Set maximum timestep limit for this cell, based on a  maximum allowed rotation angle
   //TODO, max angle could be read in from cfg
   spatial_cell->parameters[CellParams::MAXVDT]=gyro_period*(10.0/360.0);
   
  //compute initial moments, based on actual distribution function
   spatial_cell->parameters[CellParams::RHO_V  ] = 0.0;
   spatial_cell->parameters[CellParams::RHOVX_V] = 0.0;
   spatial_cell->parameters[CellParams::RHOVY_V] = 0.0;
   spatial_cell->parameters[CellParams::RHOVZ_V] = 0.0;
   
   for(unsigned int block_i=0; block_i< spatial_cell->number_of_blocks;block_i++){
      unsigned int block = spatial_cell->velocity_block_list[block_i];         
      cpu_calcVelocityMoments(spatial_cell,block,CellParams::RHO_V,CellParams::RHOVX_V,CellParams::RHOVY_V,CellParams::RHOVZ_V);   
   }

   
   
   const Real rho=spatial_cell->parameters[CellParams::RHO_V];
   //scale rho for hall term, if user requests
   const Real hallRho =  (rho <= Parameters::lorentzHallMinimumRho ) ? Parameters::lorentzHallMinimumRho : rho ;
   const Real hallPrefactor = 1.0 / (physicalconstants::MU_0 * hallRho * Parameters::q );

   Eigen::Matrix<Real,3,1> bulk_velocity(spatial_cell->parameters[CellParams::RHOVX_V]/rho,
                                 spatial_cell->parameters[CellParams::RHOVY_V]/rho,
                                 spatial_cell->parameters[CellParams::RHOVZ_V]/rho);   

   /*compute total transformation*/
   Transform<Real,3,Affine> total_transform(Matrix4d::Identity());
      
   unsigned int bulk_velocity_substeps; /*!<in this many substeps we iterate forward bulk velocity when the complete transformation is computed (0.1 deg per substep*/

   if(Parameters::lorentzHallTerm)
      bulk_velocity_substeps=dt/(gyro_period*(0.1/360.0)); 
   else
      bulk_velocity_substeps=1;
   
   /*note, we assume q is positive (pretty good assumption though)*/
   const Real substeps_radians=-(2.0*M_PI*dt/gyro_period)/bulk_velocity_substeps; /*!< how many radians each substep is*/
   for(uint i=0;i<bulk_velocity_substeps;i++){
   
      /*rotation origin is the point through which we place our rotation axis (direction of which is unitB)*/
      /*first add bulk velocity (using the total transform computed this far*/
      Eigen::Matrix<Real,3,1> rotation_pivot(total_transform*bulk_velocity);
      
      if(Parameters::lorentzHallTerm) {
         //inlude lorentzHallTerm (we should include, always)      
         rotation_pivot[0]-=hallPrefactor*(dBZdy - dBYdz);
         rotation_pivot[1]-=hallPrefactor*(dBXdz - dBZdx);
         rotation_pivot[2]-=hallPrefactor*(dBYdx - dBXdy);
      }
      
      /*add to transform matrix the small rotation around  pivot
        when added like thism, and not using *= operator, the transformations
        are in the correct order
       */
      total_transform=Translation<Real,3>(-rotation_pivot)*total_transform;
      total_transform=AngleAxis<Real>(substeps_radians,unit_B)*total_transform;
      total_transform=Translation<Real,3>(rotation_pivot)*total_transform;
   }

   return total_transform;
}

/*!
  Computes the intersection point of a plane and a line

  \param l_point Point on the line
  \param l_direction Vector in the direction of the line
  \param p_point Point on plane
  \param p_normal Normal vector to plane
  \param intersection The function will set this to the intersection point
*/


Eigen::Matrix<Real,3,1> line_plane_intersection(const Eigen::Matrix<Real,3,1>& l_point,const Eigen::Matrix<Real,3,1>& l_direction,
						const Eigen::Matrix<Real,3,1>& p_point,const Eigen::Matrix<Real,3,1>& p_normal,
						Eigen::Matrix<Real,3,1>& intersection){
  const Real nom=p_normal.dot(p_point-l_point);
  const Real dem=p_normal.dot(l_direction);
  return l_point+(nom/dem)*l_direction;
}


/*!
Computes the first intersection data; this is z~ in section 2.4 in Zerroukat et al (2012). We assume all velocity cells have the same dimensions

\param bwd_transform Transform that is used to compute the lagrangian departure grid
\param intersections map keys are indices of line (with regards to coord1 and coord2), values are value of coordinate at first intersection, index of  first intersection and index of last intersection.
\param intersection_distance Distance between intersections
*/


void intersections_1(SpatialCell* spatial_cell, std::vector<unsigned int> downstream_blocks,
		     const Transform<Real,3,Affine>& bwd_transform,
		     boost::unordered_map<std::tuple<uint64_t,uint64_t>, std::tuple<Real,int,int > >& intersections, 
		     Real& intersection_distance){
  
   Real max_z=spatial_cell->at(downstream_blocks[0])->parameters[BlockParams::VZCRD];
   Real min_z=max_z;
   Real dvz;


  /*compute the i,j indices of all blocks, use set to get rid of dupicates
   Also compute max and min z of blocks*/
  boost::unordered_set< boost::array<int,2> > block_ij;
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
    const unsigned int block = spatial_cell->velocity_block_list[block_i];
    Velocity_Block* block_ptr = spatial_cell->at(block);
    dvz=block_ptr->parameters[BlockParams::DVZ];
    const Real block_start_vz=block_ptr->parameters[BlockParams::VZCRD];
    if(block_start_vz<min_z) min_z= block_start_vz;
    if(block_start_vz+WID*dvz>max_z) max_z=block_start_vz+WID*dvz;
    block_ij.insert( {{ block%WID , (block/WID)%WID }} );

  }
  
  

  /*compute xdyd plane normal*/
  const Eigen::Matrix<Real,3,1> plane_normal=bwd_transform*Eigen::Matrix<Real,3,1>(0,0,1.0);
  const Eigen::Matrix<Real,3,1> plane_point=bwd_transform*Eigen::Matrix<Real,3,1>(0,0,min_z);
  const Eigen::Matrix<Real,3,1> plane_delta=bwd_transform*Eigen::Matrix<Real,3,1>(0,0,dvz); //vector between two lagrangian planes
  const Eigen::Matrix<Real,3,1> line_direction=Eigen::Matrix<Real,3,1>(0,0,1.0);
  

  
  Real dvx_L = plane_delta.norm2()/plane_delta.dot(line_direction); /*<Distance between lagrangian planes along line direction in Euclidian coordinates,assuming line direction has unit length!*/
  
  /*loop over block_ij*/
  for (const auto& ij: block_ij){
    for (uint cell_xi=0;cell_xi<WID;cell_xi++){
      for (uint cell_yi=0;cell_yi<WID;cell_yi++){
	Eigen::Matrix<Real,3,1> line_point((ij[0]*WID+cell_xi+0.5)*SpatialCell::block_dvx+SpatialCell::vx_min,
					   (ij[1]*WID+cell_yi+0.5)*SpatialCell::block_dvy+SpatialCell::vy_min,
					   0.0);
	
	
      }
    }
  }
}
      
/*!  

Propagates the distribution function in velocity space of given real
space cell.

Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
three‐dimensional monotone and conservative semi‐Lagrangian scheme
(SLICE‐3D) for transport problems." Quarterly Journal of the Royal
Meteorological Society 138.667 (2012): 1640-1651.

*/

void cpu_accelerate_cell(SpatialCell* spatial_cell,const Real dt) {

 
  // Make a copy of the blocklist, these are the current Eulerian cells
  std::vector<unsigned int> upstream_blocks;
  for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
    upstream_blocks.push_back(spatial_cell->velocity_block_list[block_i]);
  }
  
   /* Compute masses based on densities, and store it in fx. data is
      cleared to make way for comutation of downstream blocks, and
      finally also the actual accelerated values.
   */
   for (unsigned int block_i = 0; block_i < spatial_cell->number_of_blocks; block_i++) {
      const unsigned int block = spatial_cell->velocity_block_list[block_i];
      Velocity_Block* block_ptr = spatial_cell->at(block);
      const Real volume=block_ptr->parameters[BlockParams::DVX]*
         block_ptr->parameters[BlockParams::DVY]*
         block_ptr->parameters[BlockParams::DVZ];
      for (unsigned int cell = 0; cell < VELOCITY_BLOCK_LENGTH; cell++) {
         block_ptr->fx[cell] = block_ptr->data[cell]*volume;
         block_ptr->data[cell] = 0.0;
      }
   }
   
   /*compute transform, forward in time and backward in time*/
   phiprof::start("compute-transform");
   //compute the transform performed in this acceleration
   Transform<Real,3,Affine> fwd_transform= compute_acceleration_transformation(spatial_cell,dt);
   Transform<Real,3,Affine> bwd_transform= fwd_transform.inverse();
   phiprof::stop("compute-transform");
   
   /*compute all downstream blocks, blocks to which the distribution
     flows during this timestep, this will also create new downstream
     blocks
    */   
   std::vector<unsigned int> downstream_blocks;
   compute_downstream_blocks(spatial_cell,fwd_transform,downstream_blocks);
   
   



      
}


   

#endif

