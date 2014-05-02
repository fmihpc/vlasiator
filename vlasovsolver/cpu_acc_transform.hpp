#ifndef CPU_ACC_TRANSFORM_H
#define CPU_ACC_TRANSFORM_H

#include "common.h"
#include "spatial_cell.hpp"

#include <Eigen/Geometry>
#include <Eigen/Core>

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
   spatial_cell->parameters[CellParams::MAXVDT]=gyro_period*(P::maxSlAccelerationRotation/360.0);
   
  //compute initial moments, based on actual distribution function
   spatial_cell->parameters[CellParams::RHO_V  ] = 0.0;
   spatial_cell->parameters[CellParams::RHOVX_V] = 0.0;
   spatial_cell->parameters[CellParams::RHOVY_V] = 0.0;
   spatial_cell->parameters[CellParams::RHOVZ_V] = 0.0;
   
   for(unsigned int block_i=0; block_i< spatial_cell->number_of_blocks;block_i++){
      unsigned int block = spatial_cell->velocity_block_list[block_i];         
      cpu_calcVelocityFirstMoments(spatial_cell,block,CellParams::RHO_V,CellParams::RHOVX_V,CellParams::RHOVY_V,CellParams::RHOVZ_V);   
   }

   
   
   const Real rho=spatial_cell->parameters[CellParams::RHO_V];
   //scale rho for hall term, if user requests
   const Real hallRho =  (rho <= Parameters::lorentzHallMinimumRho ) ? Parameters::lorentzHallMinimumRho : rho ;
   const Real hallPrefactor = 1.0 / (physicalconstants::MU_0 * hallRho * Parameters::q );

   
   Eigen::Matrix<Real,3,1> bulk_velocity(spatial_cell->parameters[CellParams::RHOVX_V]/rho,
                                 spatial_cell->parameters[CellParams::RHOVY_V]/rho,
                                 spatial_cell->parameters[CellParams::RHOVZ_V]/rho);   
   /*compute total transformation*/
   Transform<Real,3,Affine> total_transform(Matrix<Real, 4, 4>::Identity()); //CONTINUE

   unsigned int bulk_velocity_substeps; /*!<in this many substeps we iterate forward bulk velocity when the complete transformation is computed (0.1 deg per substep*/
   bulk_velocity_substeps=dt/(gyro_period*(0.1/360.0)); 
   if(bulk_velocity_substeps<1)
      bulk_velocity_substeps=1;
   
   /*note, we assume q is positive (pretty good assumption though)*/
   const Real substeps_radians=-(2.0*M_PI*dt/gyro_period)/bulk_velocity_substeps; /*!< how many radians each substep is*/
   for(uint i=0;i<bulk_velocity_substeps;i++){
   
      /*rotation origin is the point through which we place our rotation axis (direction of which is unitB)*/
      /*first add bulk velocity (using the total transform computed this far*/
      Eigen::Matrix<Real,3,1> rotation_pivot(total_transform*bulk_velocity);
      
      //inlude lorentzHallTerm (we should include, always)      
      rotation_pivot[0]-=hallPrefactor*(dBZdy - dBYdz);
      rotation_pivot[1]-=hallPrefactor*(dBXdz - dBZdx);
      rotation_pivot[2]-=hallPrefactor*(dBYdx - dBXdy);

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



#endif
