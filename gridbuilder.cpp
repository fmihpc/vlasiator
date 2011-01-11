#ifndef GRIDBUILDER_H
#define GRIDBUILDER_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <map>

#include "definitions.h"
#include "common.h"
#include "parameters.h"
#include "cell_spatial.h"
#include "logger.h"
#include "project.h"

using namespace std;

extern Logger logger;

inline uint velblock(cuint& iv,cuint& jv,cuint& kv) {
   typedef Parameters P;
   return kv*P::vyblocks_ini*P::vxblocks_ini + jv*P::vxblocks_ini + iv;
}

uint cellIndex(cuint& i,cuint& j,cuint& k) {
   typedef Parameters P;
   return k*P::ycells_ini*P::xcells_ini + j*P::xcells_ini + i;
}

/** Set up a spatial cell.
 * @param cell The spatial cell which is to be initialized.
 * @param xmin x-coordinate of the lower left corner of the cell.
 * @param ymin y-coordinate of the lower left corner of the cell.
 * @param zmin z-coordinate of the lower left corner of the cell.
 * @param dx Size of the cell in x-direction.
 * @param dy Size of the cell in y-direction.
 * @param dz Size of the cell in z-direction.
 * @rvalue If true, the cell was initialized successfully. Otherwise an error has 
 * occurred and the simulation should be aborted.
 */
bool buildSpatialCell(SpatialCell& cell,creal& xmin,creal& ymin,
		      creal& zmin,creal& dx,creal& dy,creal& dz,
		     const bool& isRemote) {
   typedef Parameters P;

   cuint VELBLOCKS = P::vxblocks_ini*P::vyblocks_ini*P::vzblocks_ini;
   // Set up cell parameters:
   cell.cpu_cellParams[CellParams::XCRD] = xmin;
   cell.cpu_cellParams[CellParams::YCRD] = ymin;
   cell.cpu_cellParams[CellParams::ZCRD] = zmin;
   cell.cpu_cellParams[CellParams::DX  ] = dx;
   cell.cpu_cellParams[CellParams::DY  ] = dy;
   cell.cpu_cellParams[CellParams::DZ  ] = dz;
   calcCellParameters(cell.cpu_cellParams,0.0);
   cell.cpu_cellParams[CellParams::RHO  ] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVX] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVY] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVZ] = 0.0;
   
   // Go through each velocity block in the velocity phase space grid.
   // Set the initial state and block parameters:
   creal dvx_block = (P::vxmax-P::vxmin)/P::vxblocks_ini; // Size of a block in vx-direction
   creal dvy_block = (P::vymax-P::vymin)/P::vyblocks_ini; //                    vy
   creal dvz_block = (P::vzmax-P::vzmin)/P::vzblocks_ini; //                    vz
   creal dvx_blockCell = dvx_block / WID;                 // Size of one cell in a block in vx-direction
   creal dvy_blockCell = dvy_block / WID;                 //                                vy
   creal dvz_blockCell = dvz_block / WID;                 //                                vz
   creal dV = dvx_blockCell*dvy_blockCell*dvz_blockCell;  // Volume of one cell in a block
   Real* const blockParams = cell.cpu_blockParams;
   Real* const avgs = cell.cpu_avgs;
   uint* const nbrsVel = cell.cpu_nbrsVel;

   for (uint kv=0; kv<P::vzblocks_ini; ++kv) for (uint jv=0; jv<P::vyblocks_ini; ++jv) for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
      cuint velIndex = velblock(iv, jv, kv);
      
      creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
      creal vy_block = P::vymin + jv*dvy_block; // vy-
      creal vz_block = P::vzmin + kv*dvz_block; // vz-
      
      calcBlockParameters(blockParams + velIndex*SIZE_BLOCKPARAMS);
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::VXCRD] = vx_block;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::VYCRD] = vy_block;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::VZCRD] = vz_block;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::DVX  ] = dvx_blockCell;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::DVY  ] = dvy_blockCell;
      blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::DVZ  ] = dvz_blockCell;

      if (isRemote == true) continue;
      // Calculate volume average of distrib. function for each cell in the block.
      for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
	 creal vx_cell = vx_block + ic*dvx_blockCell;
	 creal vy_cell = vy_block + jc*dvy_blockCell;
	 creal vz_cell = vz_block + kc*dvz_blockCell;
	 avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic] =
	   calcPhaseSpaceDensity(xmin,ymin,zmin,dx,dy,dz,vx_cell,vy_cell,vz_cell,dvx_blockCell,dvy_blockCell,dvz_blockCell);
	 
	 // Add contributions to spatial cell velocity moments:
	 cell.cpu_cellParams[CellParams::RHO  ] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*dV;
	 cell.cpu_cellParams[CellParams::RHOVX] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*(vx_block + (ic+convert<Real>(0.5))*dvx_blockCell)*dV;
	 cell.cpu_cellParams[CellParams::RHOVY] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*(vy_block + (jc+convert<Real>(0.5))*dvy_blockCell)*dV;
	 cell.cpu_cellParams[CellParams::RHOVZ] += avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic]*(vz_block + (kc+convert<Real>(0.5))*dvz_blockCell)*dV;
      }
      
      // Since the whole velocity space is inside the spatial cell and the initial 
      // velocity grid contains just unrefined blocks,
      // the velocity neighbour lists can be constructed now:
      uint vxneg_nbr = iv-1;
      uint vxpos_nbr = iv+1;
      uint vyneg_nbr = jv-1;
      uint vypos_nbr = jv+1;
      uint vzneg_nbr = kv-1;
      uint vzpos_nbr = kv+1;
      uint state = 0;
      if (iv == 0) {              // -vx boundary
	 state = state | NbrsVel::VX_NEG_BND;
	 vxneg_nbr = P::vxblocks_ini-1;
      }
      if (iv == P::vxblocks_ini-1) { // +vx boundary
	 state = state | NbrsVel::VX_POS_BND;
	 vxpos_nbr = 0;
      }
      if (jv == 0) {              // -vy boundary
	 state = state | NbrsVel::VY_NEG_BND;
	 vyneg_nbr = P::vyblocks_ini-1;
      }
      if (jv == P::vyblocks_ini-1) { // +vy boundary
	 state = state | NbrsVel::VY_POS_BND;
	 vypos_nbr = 0;
      }
      if (kv == 0) {              // -vz boundary
	 state = state | NbrsVel::VZ_NEG_BND;
	 vzneg_nbr = P::vzblocks_ini-1;
      }
      if (kv == P::vzblocks_ini-1) { // +vz boundary
	 state = state | NbrsVel::VZ_POS_BND;
	 vzpos_nbr = 0;
      }
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE] = state;
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::MYIND] = velIndex;
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VXNEG] = velblock(vxneg_nbr,jv,kv);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VXPOS] = velblock(vxpos_nbr,jv,kv);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VYNEG] = velblock(iv,vyneg_nbr,kv);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VYPOS] = velblock(iv,vypos_nbr,kv);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VZNEG] = velblock(iv,jv,vzneg_nbr);
      nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VZPOS] = velblock(iv,jv,vzpos_nbr);
   }
   return true;
}
   
#endif
