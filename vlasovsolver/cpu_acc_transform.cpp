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

#include "../object_wrapper.h"
#include "cpu_moments.h"
#include "cpu_acc_transform.hpp"

using namespace std;
using namespace spatial_cell;
using namespace Eigen;


/*!
  Compute max timestep for vlasov acceleration for the particular population
  in one spatial cell.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
*/
void updateAccelerationMaxdt(
   SpatialCell* spatial_cell,
   const uint popID) 
{
   const Real EPSILON = 1e2 * numeric_limits<Real>::min();
   const Real Bx = spatial_cell->parameters[CellParams::BGBXVOL]+spatial_cell->parameters[CellParams::PERBXVOL];
   const Real By = spatial_cell->parameters[CellParams::BGBYVOL]+spatial_cell->parameters[CellParams::PERBYVOL];
   const Real Bz = spatial_cell->parameters[CellParams::BGBZVOL]+spatial_cell->parameters[CellParams::PERBZVOL];
   const Eigen::Matrix<Real,3,1> B(Bx,By,Bz);
   const Real B_mag = EPSILON > B.norm() ? EPSILON : B.norm();
   const Real gyro_period = fabs(2 * M_PI * getObjectWrapper().particleSpecies[popID].mass
				 / (getObjectWrapper().particleSpecies[popID].charge * B_mag));

   // Set maximum timestep limit for this cell, based on a maximum allowed rotation angle
   spatial_cell->set_max_v_dt(popID,gyro_period*(P::maxSlAccelerationRotation/360.0));
   
   // Constrain Vlasov solver with plasma frequency?
   if (P::ResolvePlasmaPeriod) {
     Real rho = EPSILON > spatial_cell->get_population(popID).RHO_V ? EPSILON : spatial_cell->get_population(popID).RHO_V;
     const Real plasma_period
       = fabs(2 * M_PI * sqrt(physicalconstants::EPS_0 * getObjectWrapper().particleSpecies[popID].mass / 
			      rho)/getObjectWrapper().particleSpecies[popID].charge); 
     Real smallest = gyro_period < plasma_period ? gyro_period : plasma_period;
     spatial_cell->set_max_v_dt(popID,smallest*(P::maxSlAccelerationRotation/360.0));
     }
}


/*!
 Compute transform during on timestep, and update the bulk velocity of the
 cell
 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
 * @param dt Time step of one subcycle.
*/

Eigen::Transform<Real,3,Eigen::Affine> compute_acceleration_transformation(
        SpatialCell* spatial_cell,
        const uint popID,
        const Real& dt) {

   const Real EPSILON = 1e2 * numeric_limits<Real>::min();

   // total field is BGB + perturbed B
   const Real Bx = spatial_cell->parameters[CellParams::BGBXVOL]+spatial_cell->parameters[CellParams::PERBXVOL];
   const Real By = spatial_cell->parameters[CellParams::BGBYVOL]+spatial_cell->parameters[CellParams::PERBYVOL];
   const Real Bz = spatial_cell->parameters[CellParams::BGBZVOL]+spatial_cell->parameters[CellParams::PERBZVOL];

   // read in derivatives need for curl of B (only perturbed, curl of background field is always 0!)
   const Real dBXdy = spatial_cell->derivativesBVOL[bvolderivatives::dPERBXVOLdy];
   const Real dBXdz = spatial_cell->derivativesBVOL[bvolderivatives::dPERBXVOLdz];
   const Real dBYdx = spatial_cell->derivativesBVOL[bvolderivatives::dPERBYVOLdx];

   const Real dBYdz = spatial_cell->derivativesBVOL[bvolderivatives::dPERBYVOLdz];
   const Real dBZdx = spatial_cell->derivativesBVOL[bvolderivatives::dPERBZVOLdx];
   const Real dBZdy = spatial_cell->derivativesBVOL[bvolderivatives::dPERBZVOLdy];

   const Eigen::Matrix<Real,3,1> B(Bx,By,Bz);
   Eigen::Matrix<Real,3,1> unit_B(B.normalized());

   // If B equals zero then gyro_period and unit_B are NAN.
   // Guard against that with epsilons:
   const Real B_mag = EPSILON > B.norm() ? EPSILON : B.norm();
   if (B_mag < 2*EPSILON) {
      unit_B(0,0) = 0; unit_B(1,0) = 0; unit_B(2,0) = 1;
   }

   const Real gyro_period
     = 2 * M_PI * getObjectWrapper().particleSpecies[popID].mass
     / (getObjectWrapper().particleSpecies[popID].charge * B_mag);
   const Real plasma_period
     = fabs(2 * M_PI * sqrt(physicalconstants::EPS_0 * getObjectWrapper().particleSpecies[popID].mass / 
			    spatial_cell->get_population(popID).RHO_V)/getObjectWrapper().particleSpecies[popID].charge);

   // scale rho for hall term, if user requests
   const Real rhoq = EPSILON > spatial_cell->parameters[CellParams::RHOQ_V] ? EPSILON : spatial_cell->parameters[CellParams::RHOQ_V];
   const Real hallRhoq =  (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
   const Real hallPrefactor = 1.0 / (physicalconstants::MU_0 * hallRhoq );

   // Bulk velocity is used to transform to a frame where the motional E-field vanishes
   Eigen::Matrix<Real,3,1> bulk_velocity(spatial_cell->parameters[CellParams::VX_V],
                                         spatial_cell->parameters[CellParams::VY_V],
                                         spatial_cell->parameters[CellParams::VZ_V]);

   // Use electron solvers for anything below half a proton mass
   bool smallparticle = false;
   Eigen::Matrix<Real,3,1> electronV(0.,0.,0.);
   if (getObjectWrapper().particleSpecies[popID].mass < 0.5*physicalconstants::MASS_PROTON) {
      smallparticle = true;
      // Store the original electron bulk velocity
      bulk_velocity(0,0) = electronV[0] = spatial_cell->get_population(popID).V_V[0];
      bulk_velocity(1,0) = electronV[1] = spatial_cell->get_population(popID).V_V[1];
      bulk_velocity(2,0) = electronV[2] = spatial_cell->get_population(popID).V_V[2];
   }  

    // compute total transformation
   Transform<Real,3,Affine> total_transform(Matrix<Real, 4, 4>::Identity());

   // in this many substeps we iterate forward when the complete transformation is computed (0.1 deg per substep).
   unsigned int transformation_substeps; 
   transformation_substeps = fabs(dt) / fabs(gyro_period*(0.1/360.0));
   if (smallparticle) {
      unsigned int transformation_substeps_2; 
      transformation_substeps_2 = fabs(dt) / fabs(plasma_period*(0.1/360.0));
      transformation_substeps = transformation_substeps_2 > transformation_substeps ? transformation_substeps_2 : transformation_substeps;

      // Account for extra substeps when plasma period ~ gyro period
      Real logFactor = 1.0/fabs(log( fabs(gyro_period)/fabs(plasma_period) ));
      if(std::isnan(logFactor) == false){
         //constraint: max n times more subcycles
         logFactor = min(logFactor, (Real)P::maxResonantSubcycleFactor);
         //constraint: do not lower subcycle count
         logFactor = max(logFactor, 1.0);
         transformation_substeps = int(transformation_substeps * logFactor);
      }
   }
   if ((transformation_substeps < 1) && (fabs(dt)>0)) transformation_substeps=1;
      
   const Real substeps_radians = -(2.0*M_PI*dt/gyro_period)/transformation_substeps; // how many radians each substep is.
   const Real substeps_dt=dt/(Real)transformation_substeps; /*!< how many s each substep is*/
 
   Eigen::Matrix<Real,3,1> EgradPe(
      spatial_cell->parameters[CellParams::EXGRADPE],
      spatial_cell->parameters[CellParams::EYGRADPE],
      spatial_cell->parameters[CellParams::EZGRADPE]);

   Eigen::Matrix<Real,3,1> EfromJe(
      spatial_cell->parameters[CellParams::EXJE],
      spatial_cell->parameters[CellParams::EYJE],
      spatial_cell->parameters[CellParams::EZJE]);

   // Calculate E from charge density imbalance?
   /* Eigen::Matrix<Real,3,1> Efromrq(
    spatial_cell->parameters[CellParams::ERHOQX],
    spatial_cell->parameters[CellParams::ERHOQY],
    spatial_cell->parameters[CellParams::ERHOQZ]);*/
   Eigen::Matrix<Real,3,1> Efromrq(0.0, 0.0, 0.0);
      
   const Real q = getObjectWrapper().particleSpecies[popID].charge;
   const Real mass = getObjectWrapper().particleSpecies[popID].mass;
   const Real rho = spatial_cell->get_population(popID).RHO_V;
   const Real h = substeps_dt;
   
   // Gather ion current density for electron calculations
   Eigen::Matrix<Real,3,1> Ji(0.,0.,0.);
   for (uint popID_EJE=0; popID_EJE<getObjectWrapper().particleSpecies.size(); ++popID_EJE) {
     if (getObjectWrapper().particleSpecies[popID_EJE].mass > 0.5*physicalconstants::MASS_PROTON) {
       Ji[0] += getObjectWrapper().particleSpecies[popID_EJE].charge * spatial_cell->get_population(popID_EJE).RHO_V 
	      * spatial_cell->get_population(popID_EJE).V_V[0];
       Ji[1] += getObjectWrapper().particleSpecies[popID_EJE].charge * spatial_cell->get_population(popID_EJE).RHO_V 
	 * spatial_cell->get_population(popID_EJE).V_V[1];
       Ji[2] += getObjectWrapper().particleSpecies[popID_EJE].charge * spatial_cell->get_population(popID_EJE).RHO_V 
	 * spatial_cell->get_population(popID_EJE).V_V[2];
     }
   }
   // Now account for current requirement from curl of B
   Ji[0] += (dBZdy - dBYdz)/physicalconstants::MU_0;
   Ji[1] += (dBXdz - dBZdx)/physicalconstants::MU_0;
   Ji[2] += (dBYdx - dBXdy)/physicalconstants::MU_0;

   for (uint i=0; i<transformation_substeps; ++i) {
      Eigen::Matrix<Real,3,1> dEJEt(0.,0.,0.);
      Eigen::Matrix<Real,3,1> Je(0.,0.,0.);
      Eigen::Matrix<Real,3,1> k11, k12, k21, k22;
      Eigen::Matrix<Real,3,1> k31, k32, k41, k42;
      Eigen::Matrix<Real,3,1> deltaV; 

      // rotation origin is the point through which we place our rotation axis (direction of which is unitB).
      // first add bulk velocity (using the total transform computed this far.
      Eigen::Matrix<Real,3,1> rotation_pivot(total_transform*bulk_velocity);
      
      /* include HallTerm       
	 This performs a transformation into a frame where the newly generated motional
	 electric field cancels out the Hall electric field. This is identical to the frame
	 in which electrons are, assuming that the local current density corresponds with 
	 the local curl of B.
	 If we are propagating electrons, we are already in the electron frame and the hall term
	 isn't needed.
      */
      if (!smallparticle) {
         rotation_pivot[0]-= hallPrefactor*(dBZdy - dBYdz);
         rotation_pivot[1]-= hallPrefactor*(dBXdz - dBZdx);
         rotation_pivot[2]-= hallPrefactor*(dBYdx - dBXdy);
      }
      
      // Calculate EJE only for the electron population
      if ((smallparticle) && (fabs(substeps_dt) > EPSILON)) {
	 // First find the current electron moments, this results in leapfrog-like propagation of EJE
	 Eigen::Matrix<Real,3,1> electronVcurr(total_transform*electronV);
	    
         // This is a traditional RK4 integrator 
         // In effect, it runs two RK4 integrators in parallel, one for velocity, one for electric field
         const Eigen::Matrix<Real,3,1> beta  = -q * Ji / (mass * physicalconstants::EPS_0);
         const Real alpha = -pow(q, 2.) * rho / (mass * physicalconstants::EPS_0);
         // derivative estimates for acceleration and field changes at start of step
         k11 = h * q / mass * (EfromJe+Efromrq);  // h * dv/dt
         k12 = h * (beta + alpha * electronVcurr); // h * (-q/m eps) * J_tot  ==  h * d^2 v / dt^2 == (q/m) dE/dt

         k21 = h * (q / mass * (EfromJe+Efromrq) + k12/2); // estimate acceleration using k12 field estimate (at half interval)
         k22 = h * (beta + alpha * (electronVcurr + k11/2)); // estimate field change using k11 current estimate (at half interval)

         k31 = h * (q / mass * (EfromJe+Efromrq) + k22/2); // estimate acceleration using k22 field estimate (at half interval)
         k32 = h * (beta + alpha * (electronVcurr + k21/2)); // estimate field change using k21 current estimate (at half interval)

         k41 = h * (q / mass * (EfromJe+Efromrq) + k32); // estimate acceleration using k32 field estimate (at full interval)
         k42 = h * (beta + alpha * (electronVcurr + k31)); // estimate field change using k31 current estimate (at full interval)
	    
         deltaV = (k11 + 2*k21 + 2*k31 + k41) / 6.; // Finally update velocity based on weighted acceleration estimate
         EfromJe += mass / q * (k12 + 2*k22 + 2*k32 + k42) / 6.; // And update fields based on weighted velocity (current) estimate
      } // end if (smallparticle==true) and dt>0

      if (smallparticle) {	  
         total_transform=Translation<Real,3>(deltaV) * total_transform;
      }

      /* Evaluate electron pressure gradient term. This is treated as a simple
	 nudge so we do half before and half after the rotation */
      if (Parameters::ohmGradPeTerm > 0) {
	 total_transform=Translation<Real,3>( (fabs(getObjectWrapper().particleSpecies[popID].charge)
	 				      /getObjectWrapper().particleSpecies[popID].mass) *
	 				     EgradPe * 0.5 * substeps_dt) * total_transform;
      }
      // add to transform matrix the small rotation around  pivot
      // when added like this, and not using *= operator, the transformations
      // are in the correct order
      total_transform = Translation<Real,3>(-rotation_pivot)*total_transform;
      total_transform = AngleAxis<Real>(substeps_radians,unit_B)*total_transform;
      total_transform = Translation<Real,3>(rotation_pivot)*total_transform;
      /* Evaluate second half of electron pressure gradient term. */
      if (Parameters::ohmGradPeTerm > 0) {
	 total_transform=Translation<Real,3>( (fabs(getObjectWrapper().particleSpecies[popID].charge)
	 				      /getObjectWrapper().particleSpecies[popID].mass) *
	 				     EgradPe * 0.5 * substeps_dt) * total_transform;
      }

   }
   
   // Update EJE in CELLPARAMS
   spatial_cell->parameters[CellParams::EXJE] = EfromJe[0];
   spatial_cell->parameters[CellParams::EYJE] = EfromJe[1];
   spatial_cell->parameters[CellParams::EZJE] = EfromJe[2];

   
   return total_transform;
}
