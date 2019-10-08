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
#include <fstream>
//#include <cstdlib>

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
   const Real Bx = spatial_cell->parameters[CellParams::BGBXVOL]+spatial_cell->parameters[CellParams::PERBXVOL];
   const Real By = spatial_cell->parameters[CellParams::BGBYVOL]+spatial_cell->parameters[CellParams::PERBYVOL];
   const Real Bz = spatial_cell->parameters[CellParams::BGBZVOL]+spatial_cell->parameters[CellParams::PERBZVOL];
   const Eigen::Matrix<Real,3,1> B(Bx,By,Bz);
   const Real B_mag = B.norm() + 1e-30;      
   const Real gyro_period = 2 * M_PI * getObjectWrapper().particleSpecies[popID].mass
      / (getObjectWrapper().particleSpecies[popID].charge * B_mag);

   // Set maximum timestep limit for this cell, based on a maximum allowed rotation angle
   spatial_cell->set_max_v_dt(popID,fabs(gyro_period)*(P::maxSlAccelerationRotation/360.0));
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
   // total field
   const Real Bx = spatial_cell->parameters[CellParams::BGBXVOL]+spatial_cell->parameters[CellParams::PERBXVOL];
   const Real By = spatial_cell->parameters[CellParams::BGBYVOL]+spatial_cell->parameters[CellParams::PERBYVOL];
   const Real Bz = spatial_cell->parameters[CellParams::BGBZVOL]+spatial_cell->parameters[CellParams::PERBZVOL];

   // perturbed field
   //const Real perBx = spatial_cell->parameters[CellParams::PERBXVOL];
   //const Real perBy = spatial_cell->parameters[CellParams::PERBYVOL];
   //const Real perBz = spatial_cell->parameters[CellParams::PERBZVOL];   

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
   // Guard against that by adding epsilons:
   const Real B_mag = B.norm() + 1e-30;
   if (B_mag < 1e-28) {
      unit_B(0,0) = 0; unit_B(1,0) = 0; unit_B(2,0) = 1;
   }

   const Real gyro_period
     = 2 * M_PI * getObjectWrapper().particleSpecies[popID].mass
     / (getObjectWrapper().particleSpecies[popID].charge * B_mag);
   const Real plasma_period
     = 2 * M_PI * sqrt(physicalconstants::EPS_0 * getObjectWrapper().particleSpecies[popID].mass / 
     getObjectWrapper().particleSpecies[popID].charge /
     getObjectWrapper().particleSpecies[popID].charge / spatial_cell->get_population(popID).RHO ); 

   // scale rho for hall term, if user requests
   const Real EPSILON = 1e10 * numeric_limits<Real>::min();
   const Real rhoq = spatial_cell->parameters[CellParams::RHOQ_V] + EPSILON;
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
      bulk_velocity(0,0) = spatial_cell->get_population(popID).V_V[0];
      bulk_velocity(1,0) = spatial_cell->get_population(popID).V_V[1];                                                     
      bulk_velocity(2,0) = spatial_cell->get_population(popID).V_V[2];                                                     
      // Store the original electron bulk velocity
      electronV[0] = spatial_cell->get_population(popID).V_V[0];
      electronV[1] = spatial_cell->get_population(popID).V_V[1];
      electronV[2] = spatial_cell->get_population(popID).V_V[2];
   }  

    // compute total transformation
   Transform<Real,3,Affine> total_transform(Matrix<Real, 4, 4>::Identity()); //CONTINUE

   // in this many substeps we iterate forward bulk velocity when the complete transformation is computed (0.1 deg per substep).
   unsigned int bulk_velocity_substeps; 
   bulk_velocity_substeps = fabs(dt) / fabs(gyro_period*(0.1/360.0));
   if (smallparticle) {
      unsigned int bulk_velocity_substeps_2; 
      bulk_velocity_substeps_2 = fabs(dt) / fabs(plasma_period/10.);
      bulk_velocity_substeps =  bulk_velocity_substeps_2 > bulk_velocity_substeps ? bulk_velocity_substeps_2 : bulk_velocity_substeps;
   }
   if (bulk_velocity_substeps < 1) bulk_velocity_substeps=1;
      
   const Real substeps_radians = -(2.0*M_PI*dt/gyro_period)/bulk_velocity_substeps; // how many radians each substep is.
   const Real substeps_dt=dt/bulk_velocity_substeps; /*!< how many s each substep is*/
   Eigen::Matrix<Real,3,1> EgradPe(
      spatial_cell->parameters[CellParams::EXGRADPE],
      spatial_cell->parameters[CellParams::EYGRADPE],
      spatial_cell->parameters[CellParams::EZGRADPE]);
   Eigen::Matrix<Real,3,1> EfromJe(
      spatial_cell->parameters[CellParams::EXJE],
      spatial_cell->parameters[CellParams::EYJE],
      spatial_cell->parameters[CellParams::EZJE]);

  /* 
     ofstream substepFile;
     substepFile.open ("substep_test.txt", ios::app);
     if (dt > 0.001*fabs(gyro_period)) {
     substepFile << " bulk_velocity_substeps = " << bulk_velocity_substeps << endl;
     substepFile << " dt: " << dt << endl;
     substepFile << " substeps_radians: " << substeps_radians << endl;
     substepFile << " gyro and plasma periods: " << fabs(gyro_period) << "\t" << plasma_period << endl;
     substepFile << " CellID: " << spatial_cell->parameters[CellParams::CELLID] << endl;
     substepFile << " N: " << floor(dt/fabs(gyro_period)) << endl;
     substepFile << endl;
   }
   */
      
   const Real q = getObjectWrapper().particleSpecies[popID].charge;
   const Real mass = getObjectWrapper().particleSpecies[popID].mass;
   const Real rho = spatial_cell->get_population(popID).RHO;
   const Real h = substeps_dt;

   bool RKN = true; 
   for (uint i=0; i<bulk_velocity_substeps; ++i) {
      Eigen::Matrix<Real,3,1> dEJEt(0.,0.,0.);
      Eigen::Matrix<Real,3,1> Je(0.,0.,0.);
      Eigen::Matrix<Real,3,1> Ji(0.,0.,0.);
      //Eigen::Matrix<Real,3,1> k1(0.,0.,0.);
      //Eigen::Matrix<Real,3,1> k3(0.,0.,0.);
      Eigen::Matrix<Real,3,1> k11, k12, k21, k22;
      Eigen::Matrix<Real,3,1> k31, k32, k41, k42;
      Eigen::Matrix<Real,3,1> deltaV; 
      
      // rotation origin is the point through which we place our rotation axis (direction of which is unitB).
      // first add bulk velocity (using the total transform computed this far.
      Eigen::Matrix<Real,3,1> rotation_pivot(total_transform*bulk_velocity);
      
      /* include lorentzHallTerm (we should include, always)      
	 This performs a transformation into a frame where the newly generated motional
	 electric field cancels out the Hall electric field  */
      if (!smallparticle) { //(for large particles i.e. not electrons or positrons)      
	 rotation_pivot[0]-= hallPrefactor*(dBZdy - dBYdz);
	 rotation_pivot[1]-= hallPrefactor*(dBXdz - dBZdx);
	 rotation_pivot[2]-= hallPrefactor*(dBYdx - dBXdy);
      } // For the electron run, since we use the electron bulk velocity, the Hall term should not be used
      
      // Calculate EJE only for the electron population
      if ((smallparticle) && (substeps_dt > 0)) {
	 // First find the current electron moments, this results in leapfrog-like propagation of EJE
	 Eigen::Matrix<Real,3,1> electronVcurr(total_transform*electronV);
	 if (RKN) { // Use second order solver or...
            for (uint popID_EJE=0; popID_EJE<getObjectWrapper().particleSpecies.size(); ++popID_EJE) {
               if (getObjectWrapper().particleSpecies[popID_EJE].charge > 0) {
                  Ji[0] += getObjectWrapper().particleSpecies[popID_EJE].charge * spatial_cell->get_population(popID_EJE).RHO 
		     * spatial_cell->get_population(popID_EJE).V_V[0];
                  Ji[1] += getObjectWrapper().particleSpecies[popID_EJE].charge * spatial_cell->get_population(popID_EJE).RHO 
		     * spatial_cell->get_population(popID_EJE).V_V[1];
                  Ji[2] += getObjectWrapper().particleSpecies[popID_EJE].charge * spatial_cell->get_population(popID_EJE).RHO 
		     * spatial_cell->get_population(popID_EJE).V_V[2];
               } else { // Here we assume there is no more than one negatively charged popoulation (FIXME)
                  Je[0] += getObjectWrapper().particleSpecies[popID_EJE].charge * spatial_cell->get_population(popID_EJE).RHO
		     * electronVcurr[0];
                  Je[1] += getObjectWrapper().particleSpecies[popID_EJE].charge * spatial_cell->get_population(popID_EJE).RHO 
		     * electronVcurr[1];
                  Je[2] += getObjectWrapper().particleSpecies[popID_EJE].charge * spatial_cell->get_population(popID_EJE).RHO 
		     * electronVcurr[2];
               }         
            }
            // This is a traditional RK4 integrator 
            const Eigen::Matrix<Real,3,1> beta  = -q / mass / physicalconstants::EPS_0 * Ji;
            const Real alpha = pow(q, 2.)  / mass / physicalconstants::EPS_0 * rho;
            k11 = h * q / mass * EfromJe;
            k12 = h * (beta - alpha * electronVcurr);
            k21 = h * (q / mass * EfromJe + k12/2);
            k22 = h * (beta - alpha * (electronVcurr + k11/2)); 
            k31 = h * (q / mass * EfromJe + k22/2);
            k32 = h * (beta - alpha * (electronVcurr + k21/2)); 
            k41 = h * (q / mass * EfromJe + k32);
            k42 = h * (beta - alpha * (electronVcurr + k31));
            deltaV = (k11 + 2*k21 + 2*k31 + k41) / 6.; 
            EfromJe += mass / q * (k12 + 2*k22 + 2*k32 + k42) / 6.; 
            
            /* This RKN solver did not work
	       k1 = -getObjectWrapper().particleSpecies[popID].charge / getObjectWrapper().particleSpecies[popID].mass 
	       / physicalconstants::EPS_0 * (Ji + Je);
	       k3 = -getObjectWrapper().particleSpecies[popID].charge / getObjectWrapper().particleSpecies[popID].mass 
	       / physicalconstants::EPS_0 * 
	       ( Ji + getObjectWrapper().particleSpecies[popID].charge * spatial_cell->get_population(popID).RHO
	       * (  electronVcurr + 2./3. * substeps_dt * EfromJe * getObjectWrapper().particleSpecies[popID].charge 
	       / getObjectWrapper().particleSpecies[popID].mass + 2./9. * pow(substeps_dt,2.) * k1 ) );
	       dEJEt = getObjectWrapper().particleSpecies[popID].mass / getObjectWrapper().particleSpecies[popID].charge
	       * substeps_dt / 4. * (k1 + 3. * k3); 
	       */
	    
         } else {  
	    // Use simple Eulerian solver
	    /* Calculate electrostatic field derivative via current
	       using the substep-transformed electron bulkV and exising other population bulkVs */
	    for (uint popID_EJE=0; popID_EJE<getObjectWrapper().particleSpecies.size(); ++popID_EJE) {
	       if (getObjectWrapper().particleSpecies[popID_EJE].mass < 0.5*physicalconstants::MASS_PROTON) {
		  dEJEt[0] += -getObjectWrapper().particleSpecies[popID_EJE].charge *spatial_cell->get_population(popID_EJE).RHO
		     * electronVcurr[0] / physicalconstants::EPS_0;
		  dEJEt[1] += -getObjectWrapper().particleSpecies[popID_EJE].charge *spatial_cell->get_population(popID_EJE).RHO
		     * electronVcurr[1] / physicalconstants::EPS_0;
		  dEJEt[2] += -getObjectWrapper().particleSpecies[popID_EJE].charge *spatial_cell->get_population(popID_EJE).RHO
		     * electronVcurr[2] / physicalconstants::EPS_0;
	       } else {
		  dEJEt[0] += -getObjectWrapper().particleSpecies[popID_EJE].charge *spatial_cell->get_population(popID_EJE).RHO
		     * spatial_cell->get_population(popID_EJE).V_V[0] / physicalconstants::EPS_0;
		  dEJEt[1] += -getObjectWrapper().particleSpecies[popID_EJE].charge *spatial_cell->get_population(popID_EJE).RHO
		     * spatial_cell->get_population(popID_EJE).V_V[1] / physicalconstants::EPS_0;
		  dEJEt[2] += -getObjectWrapper().particleSpecies[popID_EJE].charge *spatial_cell->get_population(popID_EJE).RHO
		     * spatial_cell->get_population(popID_EJE).V_V[2] / physicalconstants::EPS_0;
	       }         
	    }
	    
	    // Increment EfromJe with derivative times half of substep to get representative field throughout integration step
	    EfromJe += dEJEt*0.5*substeps_dt;
	    // Find B-perpendicular and B-parallel components of EfromJe
	    Eigen::Matrix<Real,3,1> EfromJe_parallel(EfromJe.dot(unit_B)*unit_B);
	    Eigen::Matrix<Real,3,1> EfromJe_perpendicular(EfromJe-EfromJe_parallel);
	    Eigen::Matrix<Real,3,1> unit_EJEperp(EfromJe_perpendicular.normalized());
	    
	    // Add pivot transformation to negate component of EfromJe perpendicular to B
	    // Vnorm = Bnorm cross Enorm
	    Real EJEperpperB = EfromJe_perpendicular.norm() / B.norm();
	    rotation_pivot[0]-= EJEperpperB * (unit_B[1]*unit_EJEperp[2] - unit_B[2]*unit_EJEperp[1]);
	    rotation_pivot[1]-= EJEperpperB * (unit_B[2]*unit_EJEperp[0] - unit_B[0]*unit_EJEperp[2]);
	    rotation_pivot[2]-= EJEperpperB * (unit_B[0]*unit_EJEperp[1] - unit_B[1]*unit_EJEperp[0]);
	    
         }
	 
         /* if (getObjectWrapper().particleSpecies[popID].charge < 0 && int(spatial_cell->parameters[CellParams::CELLID]) == 1 &&
	    i % 10 == 0) {
	    substepFile << " EJEperpperB :\n" << EJEperpperB << endl;
	    substepFile << " EfromJe: \n" << EfromJe << endl;
	    substepFile << " EfromJe_parallel: \n" << EfromJe_parallel << endl;
	    substepFile << " EfromJe_perpendicular: \n" << EfromJe_perpendicular << endl;
	    substepFile << " unit_EJEperp: \n" << unit_EJEperp << endl;
	    substepFile << " electronVcurr: \n" << electronVcurr << endl;
            substepFile << " substep_i: \n" << i << endl;
            substepFile << " dEJEt: \n" << dEJEt << endl;
            substepFile << " V0: \n" << electronVcurr << endl;
            substepFile << " k11: \n" << k11 << endl;
            substepFile << " k12: \n" << k12 << endl;
            substepFile << " Ji: \n" << Ji << endl;
            substepFile << " Je: \n" << Je << endl;
            substepFile << " deltaV: \n" << deltaV << endl;
            substepFile << " rotation_pivot matrix at 2 :\n" << rotation_pivot << endl;
	    substepFile << " CellID: " << spatial_cell->parameters[CellParams::CELLID] << endl;
	    substepFile << " electronV: \n" << electronV << endl;
	    substepFile << " substeps_dt, dt: " << substeps_dt << "\t" << dt << endl;
	    substepFile << " substeps_radians: \n" << substeps_radians << endl;
            substepFile << " " << endl;
	    } */
      } // end if (smallparticle==true) and dt>0
      
      // add to transform matrix the small rotation around  pivot
      // when added like this, and not using *= operator, the transformations
      // are in the correct order
      total_transform = Translation<Real,3>(-rotation_pivot)*total_transform;
      total_transform = AngleAxis<Real>(substeps_radians,unit_B)*total_transform;
      total_transform = Translation<Real,3>(rotation_pivot)*total_transform;
      
      if (smallparticle) {
	 Eigen::Matrix<Real,3,1> EfromJe_parallel(EfromJe.dot(unit_B)*unit_B);
	  
	 if (RKN) {
	    total_transform=Translation<Real,3>( (getObjectWrapper().particleSpecies[popID].charge/
						  getObjectWrapper().particleSpecies[popID].mass) * 
						 EfromJe * substeps_dt ) * total_transform;
	    //total_transform=Translation<Real,3>( (getObjectWrapper().particleSpecies[popID].charge/
	    //   getObjectWrapper().particleSpecies[popID].mass) * 
	    //   EfromJe * substeps_dt + (k1 + k3) * substeps_dt / 4.) * total_transform;

	 } else {
	    // If using the Eulerian scheme, then the rotation algorithm is used
	    // and only the parallel nudge is required:
	    // Perform B-parallel acceleration from EJE field
	    total_transform=Translation<Real,3>( (getObjectWrapper().particleSpecies[popID].charge/
						  getObjectWrapper().particleSpecies[popID].mass) * 
						 EfromJe_parallel * substeps_dt) * total_transform;

	    /* The alternative to decomposing the EJE field into parallel and perpendicular components is to
	       treat it a simple acceleration term. This acceleration was found by integrating over
	       the time-varying electric field. */
	    //total_transform=Translation<Real,3>( (getObjectWrapper().particleSpecies[popID].charge/
	    //   getObjectWrapper().particleSpecies[popID].mass) * 
	    //   EfromJe * substeps_dt) * total_transform;
	    EfromJe += dEJEt*0.5*substeps_dt;
	 }
      }
      
      // if (getObjectWrapper().particleSpecies[popID].charge < 0 && 
      //     int(spatial_cell->parameters[CellParams::CELLID]) == 1 && i % 10 == 0) {
      //   substepFile << " total_transform: \n" << total_transform.matrix() << endl;
      //   substepFile << " " << endl;
      // }
   
      // Electron pressure gradient term (this is still untested and might also need to be decomposed into perp and parallel portions)
      if((!smallparticle) && (Parameters::ohmGradPeTerm > 0)) {
	 total_transform=Translation<Real,3>( (fabs(getObjectWrapper().particleSpecies[popID].charge)
					       /getObjectWrapper().particleSpecies[popID].mass) *
					      EgradPe * substeps_dt) * total_transform;
      }
   }
   //substepFile.close();
   
   // Update EJE in CELLPARAMS
   spatial_cell->parameters[CellParams::EXJE] = EfromJe[0];
   spatial_cell->parameters[CellParams::EYJE] = EfromJe[1];
   spatial_cell->parameters[CellParams::EZJE] = EfromJe[2];

   
   return total_transform;
}
