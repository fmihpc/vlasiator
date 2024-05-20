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

#include "common.h"
#include "spatial_cell_wrapper.hpp"

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
   const Real dBXdy = spatial_cell->derivativesBVOL[bvolderivatives::dPERBXVOLdy];
   const Real dBXdz = spatial_cell->derivativesBVOL[bvolderivatives::dPERBXVOLdz];
   const Real dBYdx = spatial_cell->derivativesBVOL[bvolderivatives::dPERBYVOLdx];

   const Real dBYdz = spatial_cell->derivativesBVOL[bvolderivatives::dPERBYVOLdz];
   const Real dBZdx = spatial_cell->derivativesBVOL[bvolderivatives::dPERBZVOLdx];
   const Real dBZdy = spatial_cell->derivativesBVOL[bvolderivatives::dPERBZVOLdy];

   
   const Eigen::Matrix<Real,3,1> B(Bx,By,Bz);
   const Eigen::Matrix<Real,3,1> unit_B(B.normalized());

   #warning This is wrong for multipop 
   const Real gyro_period = 2 * M_PI * physicalconstants::MASS_PROTON  / (fabs(physicalconstants::CHARGE) * B.norm());
   
   //Set maximum timestep limit for this cell, based on a  maximum allowed rotation angle
   spatial_cell->parameters[CellParams::MAXVDT]=gyro_period*(P::maxSlAccelerationRotation/360.0);
   
  //compute initial moments, based on actual distribution function
   spatial_cell->parameters[CellParams::RHOM_V  ] = 0.0;
   spatial_cell->parameters[CellParams::VX_V] = 0.0;
   spatial_cell->parameters[CellParams::VY_V] = 0.0;
   spatial_cell->parameters[CellParams::VZ_V] = 0.0;
   spatial_cell->parameters[CellParams::RHOQ_V] = 0.0;
   
   for (vmesh::LocalID block_i=0; block_i<spatial_cell->get_number_of_velocity_blocks(); ++block_i) {
      cpu_calcVelocityFirstMoments(spatial_cell,block_i,CellParams::RHOM_V,CellParams::VX_V,CellParams::VY_V,CellParams::VZ_V);
   }
   
   const Real rhoq=spatial_cell->parameters[CellParams::RHOQ_V];
   //scale rho for hall term, if user requests
   const Real hallRhoq =  (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
   const Real hallPrefactor = 1.0 / (physicalconstants::MU_0 * hallRhoq );

   
   Eigen::Matrix<Real,3,1> bulk_velocity(spatial_cell->parameters[CellParams::VX_V],
                                 spatial_cell->parameters[CellParams::VY_V],
                                 spatial_cell->parameters[CellParams::VZ_V]);
   /*compute total transformation*/
   Transform<Real,3,Affine> total_transform(Matrix<Real, 4, 4>::Identity()); //CONTINUE

   unsigned int bulk_velocity_substeps; /*!<in this many substeps we iterate forward bulk velocity when the complete transformation is computed (0.1 deg per substep*/
   bulk_velocity_substeps=fabs(dt)/(gyro_period*(0.1/360.0)); 
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
