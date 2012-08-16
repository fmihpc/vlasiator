/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef RIEMANN_H
#define RIEMANN_H

#include "definitions.h"
#include "spatial_cell.hpp"
#include "projects/projects_common.h"
#include "projects/projects_fieldboundary.h"
#include "projects/projects_vlasov_acceleration.h"
#include "projects/projects_vlasov_boundary.h"
#include "fieldsolver.h"

#include "dccrg.hpp"


struct kelvinHelmholtzParameters {
   enum {
      TOP,
      BOTTOM
   };
   static Real rho[2];
   static Real T[2];
   static Real Vx[2];
   static Real Vy[2];
   static Real Vz[2];
   static Real Bx[2];
   static Real By[2];
   static Real Bz[2];
   static Real lambda;
   static Real amp;
   static Real offset;
   static Real transitionWidth;
   static uint nSpaceSamples;
   static uint nVelocitySamples;
};

typedef kelvinHelmholtzParameters KHP;

/**
 * Initialize project. Can be used, e.g., to read in parameters from the input file
 */
bool initializeProject(void);

/** Register parameters that should be read in
 */
bool addProjectParameters(void);
/** Get the value that was read in
 */
bool getProjectParameters(void);

/** Query if spatial cell parameters (of any cell) have changed and need to be 
 * recalculated. If you have a completely static case, then you can always return 
 * false here. Otherwise you need to return true so that function calcCellParameters
 * gets called for each spatial cell.
 * @param t The current value of time.
 * @return If true, the spatial cell parameters need to be recalculated.
 */
bool cellParametersChanged(creal& t);

/** Calculate parameters for the given velocity grid block.
 * Here you need to set a value for array indices: BlockParams::Q_PER_M.
 * @param blockParams Array containing the block parameters
 */
void calcBlockParameters(Real* blockParams);

/** Calculate parameters for the given spatial cell at the given time.
 * Here you need to set values for the following array indices:
 * CellParams::EX, CellParams::EY, CellParams::EZ, CellParams::BX, CellParams::BY, and CellParams::BZ.
 * 
 * The following array indices contain the coordinates of the "lower left corner" of the cell: 
 * CellParams::XCRD, CellParams::YCRD, and CellParams::ZCRD.
 * The cell size is given in the following array indices: CellParams::DX, CellParams::DY, and CellParams::DZ.
 * @param cellParams Array containing cell parameters.
 * @param t The current value of time. This is passed as a convenience. If you need more detailed information 
 * of the state of the simulation, you can read it from Parameters.
 */
void calcCellParameters(Real* cellParams,creal& t);

void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& dt);

/** Integrate the distribution function over the given six-dimensional phase-space cell.
 * @param x Starting value of the x-coordinate of the cell.
 * @param y Starting value of the y-coordinate of the cell.
 * @param z Starting value of the z-coordinate of the cell.
 * @param dx The size of the cell in x-direction.
 * @param dy The size of the cell in y-direction.
 * @param dz The size of the cell in z-direction.
 * @param vx Starting value of the vx-coordinate of the cell.
 * @param vy Starting value of the vy-coordinate of the cell.
 * @param vz Starting value of the vz-coordinate of the cell.
 * @param dvx The size of the cell in vx-direction.
 * @param dvy The size of the cell in vy-direction.
 * @param dvz The size of the cell in vz-direction.
 * @return The volume average of the distribution function in the given phase space cell.
 * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
 */
Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz);

/* Split out because it is used in system boundary conditions which are initialised using the same values. */
//let's make it into a project function -> move it there
template<class SC>
void setProjectCell(SC* cell) {
   // Set up cell parameters:
   calcCellParameters(&((*cell).parameters[0]), 0.0);
   cell->parameters[CellParams::RHO  ] = 0.0;
   cell->parameters[CellParams::RHOVX] = 0.0;
   cell->parameters[CellParams::RHOVY] = 0.0;
   cell->parameters[CellParams::RHOVZ] = 0.0;

   cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
   cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
   
   // Go through each velocity block in the velocity phase space grid.
   // Set the initial state and block parameters:
   creal dvx_block = SpatialCell::block_dvx; // Size of a block in vx-direction
   creal dvy_block = SpatialCell::block_dvy; //                    vy
   creal dvz_block = SpatialCell::block_dvz; //                    vz
   creal dvx_blockCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
   creal dvy_blockCell = SpatialCell::cell_dvy; //                                vy
   creal dvz_blockCell = SpatialCell::cell_dvz; //                                vz
   
   for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
      for (uint jv=0; jv<P::vyblocks_ini; ++jv)
         for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
            creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
            creal vy_block = P::vymin + jv*dvy_block; // vy-
            creal vz_block = P::vzmin + kv*dvz_block; // vz-
            
            // Calculate volume average of distrib. function for each cell in the block.
            for (uint kc=0; kc<WID; ++kc) 
               for (uint jc=0; jc<WID; ++jc) 
                  for (uint ic=0; ic<WID; ++ic) {
                     creal vx_cell = vx_block + ic*dvx_blockCell;
                     creal vy_cell = vy_block + jc*dvy_blockCell;
                     creal vz_cell = vz_block + kc*dvz_blockCell;
                     Real average = 
                     calcPhaseSpaceDensity(cell->parameters[CellParams::XCRD],
                                           cell->parameters[CellParams::YCRD],
                                           cell->parameters[CellParams::ZCRD],
                                           cell->parameters[CellParams::DX],
                                           cell->parameters[CellParams::DY],
                                           cell->parameters[CellParams::DZ],
                                           vx_cell,vy_cell,vz_cell,
                                           dvx_blockCell,dvy_blockCell,dvz_blockCell);
                     
                     if(average!=0.0){
                        creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
                        creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
                        creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
                        cell->set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
                        // Add contributions to spatial cell velocity moments:
                        creal dV = dvx_blockCell*dvy_blockCell*dvz_blockCell;  // Volume of one cell in a block      
                        cell->parameters[CellParams::RHO  ] += average*dV;
                        cell->parameters[CellParams::RHOVX] += average*vx_cell_center*dV;
                        cell->parameters[CellParams::RHOVY] += average*vy_cell_center*dV;
                        cell->parameters[CellParams::RHOVZ] += average*vz_cell_center*dV;
                     }
                  }
         }
         creal spatialVolume = cell->parameters[CellParams::DX]*
         cell->parameters[CellParams::DY]*
         cell->parameters[CellParams::DZ];
         cell->parameters[CellParams::RHO  ] /= spatialVolume;
         cell->parameters[CellParams::RHOVX] /= spatialVolume;
         cell->parameters[CellParams::RHOVY] /= spatialVolume;
         cell->parameters[CellParams::RHOVZ] /= spatialVolume;
         
         //lets get rid of blocks not fulfilling the criteria here to save
         //memory. neighbor_ptrs is empty as we do not have any consistent
         //data in neighbours yet, adjustments done only based on velocity
         //space.
         //SVA: move the next three lines into a separate function in spatial cell class
         std::vector<SC*> neighbor_ptrs;
         cell->update_all_block_has_content();
         cell->adjust_velocity_blocks(neighbor_ptrs);
}

/** Calculate the boundary value of volume average of distribution function. This function 
 * should calculate the value of distribution function on the other side of given phase space 
 * cell face. The coordinate direction is given with parameter crd, and positive or negative 
 * face is given with parameter negSide. For example, if crd == 0 and negSide == true, then 
 * boundary value is requested at -x face of the given phase space cell.
 * 
 * The spatial coordinate values can be calculated using cellParams. This array contains the parameters 
 * of the spatial cell which is inside the simulation domain, and you should assume that the 
 * ghost cell has the same (dx,dy,dz) values. For example, if crd == 0 and negSide == true, then 
 * the spatial coordinates of the lower left corner of the ghost cell are 
 * (cellParams[XCRD]-cellParams[DX], cellParams[YCRD], cellParams[ZCRD]).
 * 
 * The velocity coordinates can be calculated using blockParams. This array contains the parameters 
 * of the velocity grid block which is just inside the simulation domain, and you should assume that 
 * the corresponding velocity grid block of the ghost cell has the same values. For example, the 
 * velocity coordinates of the lower left corner of the velocity block are 
 * (blockParams[VXCRD]+iv*blockParams[DVX], blockParams[VYCRD]+jv*blockParams[DVY], 
 * blockParams[VZCRD]+kv*blockParams[DVZ]).
 * 
 * Note that this function does not need to be a template.
 * 
 * @param iv The vx-index of the cell in velocity block.
 * @param jv The vy-index of the cell in velocity block.
 * @param kv The vz-index of the cell in velocity block.
 * @param cellParams Array containing the spatial cell parameters.
 * @param blockParams Array containing the velocity block parameters.
 * @param avg Volume average of distribution function in the velocity block cell just inside the 
 * simulation domain. 
 * @param crd The spatial coordinate direction (0 = x, 1 = y, 2 = z).
 * @param negSide If true, then the boundary value at the negative coordinate side is requested.
 * @return Volume average of distribution function at the given boundary. The physical 
 * unit of this quantity is 1 / (m^3 (m/s)^3).
 */
template<typename T>
T calcBoundVolAvg(cuint& iv,cuint& jv,cuint& kv,const T* const cellParams,
		  const T* const blockParams,const T& avg,const int& crd,const bool& negSide) {
   return avg;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundaryCondDerivX(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,creal* const derivatives,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
   array[fieldsolver::drhodx] = 0.0;
   array[fieldsolver::dBydx]  = 0.0;
   array[fieldsolver::dBzdx]  = 0.0;
   array[fieldsolver::dVxdx]  = 0.0;
   array[fieldsolver::dVydx]  = 0.0;
   array[fieldsolver::dVzdx]  = 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundaryCondDerivY(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,creal* const derivatives,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
   array[fieldsolver::drhody] = 0.0;
   array[fieldsolver::dBxdy]  = 0.0;
   array[fieldsolver::dBzdy]  = 0.0;
   array[fieldsolver::dVxdy]  = 0.0;
   array[fieldsolver::dVydy]  = 0.0;
   array[fieldsolver::dVzdy]  = 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundaryCondDerivZ(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,creal* const derivatives,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
   array[fieldsolver::drhodz] = 0.0;
   array[fieldsolver::dBxdz]  = 0.0;
   array[fieldsolver::dBydz]  = 0.0;
   array[fieldsolver::dVxdz]  = 0.0;
   array[fieldsolver::dVydz]  = 0.0;
   array[fieldsolver::dVzdz]  = 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldSolverBoundaryCondBx(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
   Real Bx;
   creal x = mpiGrid[cellID]->parameters[CellParams::XCRD];
   creal z = mpiGrid[cellID]->parameters[CellParams::ZCRD];
   creal d_z = mpiGrid[cellID]->parameters[CellParams::DZ];
   if(KHP::offset != 0.0) {
      Bx = ((z + 0.5 * d_z) < KHP::offset + KHP::amp * cos(2.0 * M_PI * x / KHP::lambda)) &&
      ((z + 0.5 * d_z) > -1.0 * KHP::offset + KHP::amp * cos(2.0 * M_PI * x / KHP::lambda)) ?
      KHP::Bx[KHP::TOP] : KHP::Bx[KHP::BOTTOM];
      
   } else {
      Bx = ((z + 0.5 * d_z) > 0.0) ? KHP::Bx[KHP::TOP] : KHP::Bx[KHP::BOTTOM];
   }
   return Bx;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldSolverBoundaryCondBy(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
   Real By;
   creal x = mpiGrid[cellID]->parameters[CellParams::XCRD];
   creal z = mpiGrid[cellID]->parameters[CellParams::ZCRD];
   creal d_z = mpiGrid[cellID]->parameters[CellParams::DZ];
   if(KHP::offset != 0.0) {
      By = ((z + 0.5 * d_z) < KHP::offset + KHP::amp * cos(2.0 * M_PI * x / KHP::lambda)) &&
      ((z + 0.5 * d_z) > -1.0 * KHP::offset + KHP::amp * cos(2.0 * M_PI * x / KHP::lambda)) ?
      KHP::By[KHP::TOP] : KHP::By[KHP::BOTTOM];
      
   } else {
      By = ((z + 0.5 * d_z) > 0.0) ? KHP::By[KHP::TOP] : KHP::By[KHP::BOTTOM];
   }
   return By;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldSolverBoundaryCondBz(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
   Real Bz;
   creal x = mpiGrid[cellID]->parameters[CellParams::XCRD];
   creal z = mpiGrid[cellID]->parameters[CellParams::ZCRD];
   creal d_z = mpiGrid[cellID]->parameters[CellParams::DZ];
   if(KHP::offset != 0.0) {
      Bz = ((z + 0.5 * d_z) < KHP::offset + KHP::amp * cos(2.0 * M_PI * x / KHP::lambda)) &&
      ((z + 0.5 * d_z) > -1.0 * KHP::offset + KHP::amp * cos(2.0 * M_PI * x / KHP::lambda)) ?
      KHP::Bz[KHP::TOP] : KHP::Bz[KHP::BOTTOM];
      
   } else {
      Bz = ((z + 0.5 * d_z) > 0.0) ? KHP::Bz[KHP::TOP] : KHP::Bz[KHP::BOTTOM];
   }
   return Bz;
}

template<typename CELLID,typename UINT> 
void vlasovBoundaryCondition(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
   return;
}

template<typename UINT,typename REAL> void calcAccFaceX(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,const REAL* const cellParams,const REAL* const blockParams) {
   lorentzForceFaceX(ax,ay,az,I,J,K,cellParams,blockParams);
}
   
template<typename UINT,typename REAL> void calcAccFaceY(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,const REAL* const cellParams,const REAL* const blockParams) {
   lorentzForceFaceY(ax,ay,az,I,J,K,cellParams,blockParams);
}

template<typename UINT,typename REAL> void calcAccFaceZ(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,const REAL* const cellParams,const REAL* const blockParams) {
   lorentzForceFaceZ(ax,ay,az,I,J,K,cellParams,blockParams);
}

#endif

