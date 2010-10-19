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
#include "grid.h"
#include "devicegrid.h"
#include "logger.h"
#include "project.h"

using namespace std;

extern Logger logger;

inline uint velblock(cuint& iv,cuint& jv,cuint& kv) {
   typedef Parameters P;
   return kv*P::vyblocks_ini*P::vxblocks_ini + jv*P::vxblocks_ini + iv;
}

inline uint index(cuint& i,cuint& j,cuint& k) {
   typedef Parameters P;
   return k*P::ycells_ini*P::xcells_ini + j*P::xcells_ini + i;
}

bool buildBaseGrid(Grid& grid,DeviceGrid& deviceGrid) {
   typedef Parameters P;
   
   real x,y,z;
   creal dx = (P::xmax-P::xmin)/P::xcells_ini;
   creal dy = (P::ymax-P::ymin)/P::ycells_ini;
   creal dz = (P::zmax-P::zmin)/P::zcells_ini;

   cuint VELBLOCKS = P::vxblocks_ini*P::vyblocks_ini*P::vzblocks_ini;
   uint nbr_xneg,nbr_xpos,nbr_yneg,nbr_ypos,nbr_zneg,nbr_zpos;
   uint spaCellIndex;
   map<uint,uint> cellIndices;
   
   // Go through each spatial cell and set its parameters:
   for (uint k=0; k<P::zcells_ini; ++k) for (uint j=0; j<P::ycells_ini; ++j) for (uint i=0; i<P::xcells_ini; ++i) {
      x = P::xmin + i*dx;
      y = P::ymin + j*dy;
      z = P::zmin + k*dz;

      // Initialize the spatial cell:
      spaCellIndex = grid.getSpatialCell(VELBLOCKS);
      if (spaCellIndex == numeric_limits<uint>::max()) {
	 logger << "(GRIDBUILDER): Failed to get a spatial cell from grid!" << endl;
	 return false;
      }
      deviceGrid.initSpatialCell(*(grid[spaCellIndex]));
      cellIndices[index(i,j,k)] = spaCellIndex;
      
      grid[spaCellIndex]->i_ind = i;
      grid[spaCellIndex]->j_ind = j;
      grid[spaCellIndex]->k_ind = k;
      grid[spaCellIndex]->cellType = Cell::INNER;
      
      // Calculate parameters for the spatial cell:
      grid[spaCellIndex]->cpu_cellParams[CellParams::XCRD ] = x;
      grid[spaCellIndex]->cpu_cellParams[CellParams::YCRD ] = y;
      grid[spaCellIndex]->cpu_cellParams[CellParams::ZCRD ] = z;
      grid[spaCellIndex]->cpu_cellParams[CellParams::DX   ] = dx;
      grid[spaCellIndex]->cpu_cellParams[CellParams::DY   ] = dy;
      grid[spaCellIndex]->cpu_cellParams[CellParams::DZ   ] = dz;
      calcCellParameters(grid[spaCellIndex]->cpu_cellParams);
      grid[spaCellIndex]->cpu_cellParams[CellParams::RHO  ] = 0.0;
      grid[spaCellIndex]->cpu_cellParams[CellParams::RHOVX] = 0.0;
      grid[spaCellIndex]->cpu_cellParams[CellParams::RHOVY] = 0.0;
      grid[spaCellIndex]->cpu_cellParams[CellParams::RHOVZ] = 0.0;
   }
   
   // Now that the spatial cells exist, go through each one again and set up spatial neighbour lists 
   // and initial state for each velocity grid block:
   uint iniVelBlock = 0;
   for (uint k=0; k<P::zcells_ini; ++k) for (uint j=0; j<P::ycells_ini; ++j) for (uint i=0; i<P::xcells_ini; ++i) {      
      // Set spatial neighbour list: (periodic boundaries for now)
      nbr_xneg = i-1;
      nbr_xpos = i+1;
      nbr_yneg = j-1;
      nbr_ypos = j+1;
      nbr_zneg = k-1;
      nbr_zpos = k+1;
      if (nbr_xneg >= P::xcells_ini) nbr_xneg = P::xcells_ini-1;
      if (nbr_xpos == P::xcells_ini) nbr_xpos = 0;
      if (nbr_yneg >= P::ycells_ini) nbr_yneg = P::ycells_ini-1;
      if (nbr_ypos == P::ycells_ini) nbr_ypos = 0;
      if (nbr_zneg >= P::zcells_ini) nbr_zneg = P::zcells_ini-1;
      if (nbr_zpos == P::zcells_ini) nbr_zpos = 0;
      
      spaCellIndex = cellIndices[index(i,j,k)];
      grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::STATE] = NbrsSpa::INNER;
      grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::MYIND] = spaCellIndex;
      grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::XNEG] = cellIndices[index(nbr_xneg,j,k)];
      grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::XPOS] = cellIndices[index(nbr_xpos,j,k)];
      grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::YNEG] = cellIndices[index(i,nbr_yneg,k)];
      grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::YPOS] = cellIndices[index(i,nbr_ypos,k)];
      grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::ZNEG] = cellIndices[index(i,j,nbr_zneg)];
      grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::ZPOS] = cellIndices[index(i,j,nbr_zpos)];
      /*
      cout << "Spatial Cell " << i << ' ' << j << ' ' << k << endl;
      cout << "\t" << nbr_xneg << ' ' << nbr_xpos << ' ' << nbr_yneg << ' ' << nbr_ypos << ' ' << nbr_zneg << ' ' << nbr_zpos << endl;
      cout << "\t" << cellIndices[index(nbr_xneg,j,k)] << ' ' << cellIndices[index(nbr_xpos,j,k)] << endl;
      cout << "\t" << cellIndices[index(i,nbr_yneg,k)] << ' ' << cellIndices[index(i,nbr_ypos,k)] << endl;
      cout << "\t" << cellIndices[index(i,j,nbr_zneg)] << ' ' << cellIndices[index(i,j,nbr_zpos)] << endl;
      */
      // Set initial state for the velocity grid blocks
      //buildVelocityGrid(grid,spaCellIndex,i,j,k);

      /*
      // Request an index from the device. We need this gpu index to calculate correct values 
      // for the neighbour indices.
      Device::BlockType blockType = Device::INNER;
      if (grid[index(i,j,k)]->cpu_nbrsSpa[NbrsSpa::STATE] != NbrsSpa::INNER) blockType = Device::BOUNDARY;
      grid[index(i,j,k)]->gpuIndex = deviceGrid.getBlocks(VELBLOCKS,blockType);
      */
      
      creal x_cell = grid[spaCellIndex]->cpu_cellParams[CellParams::XCRD];
      creal y_cell = grid[spaCellIndex]->cpu_cellParams[CellParams::YCRD];
      creal z_cell = grid[spaCellIndex]->cpu_cellParams[CellParams::ZCRD];
      creal dx_cell = grid[spaCellIndex]->cpu_cellParams[CellParams::DX];
      creal dy_cell = grid[spaCellIndex]->cpu_cellParams[CellParams::DY];
      creal dz_cell = grid[spaCellIndex]->cpu_cellParams[CellParams::DZ];
      creal dvx_block = (P::vxmax-P::vxmin)/P::vxblocks_ini;
      creal dvy_block = (P::vymax-P::vymin)/P::vyblocks_ini;
      creal dvz_block = (P::vzmax-P::vzmin)/P::vzblocks_ini;
      creal dvx_blockCell = dvx_block / WID;
      creal dvy_blockCell = dvy_block / WID;
      creal dvz_blockCell = dvz_block / WID;
      creal dV = dvx_blockCell*dvy_blockCell*dvz_blockCell;
      
      for (uint kv=0; kv<P::vzblocks_ini; ++kv) for (uint jv=0; jv<P::vyblocks_ini; ++jv) for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
	 cuint velIndex = velblock(iv,jv,kv);
	 
	 creal vx_block = P::vxmin + iv*dvx_block;
	 creal vy_block = P::vymin + jv*dvy_block;
	 creal vz_block = P::vzmin + kv*dvz_block;
	 
	 // Set parameters of the velocity block:
	 // NOTE: THIS MIGHT NEED IMPROVEMENT
	 calcBlockParameters(grid[spaCellIndex]->cpu_blockParams + velIndex*SIZE_BLOCKPARAMS);
	 //grid[spaCellIndex]->cpu_blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::Q_PER_M] = 1.0;
	 grid[spaCellIndex]->cpu_blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::VXCRD  ] = vx_block;
	 grid[spaCellIndex]->cpu_blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::VYCRD  ] = vy_block;
	 grid[spaCellIndex]->cpu_blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::VZCRD  ] = vz_block;
	 grid[spaCellIndex]->cpu_blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::DVX    ] = dvx_blockCell;
	 grid[spaCellIndex]->cpu_blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::DVY    ] = dvy_blockCell;
	 grid[spaCellIndex]->cpu_blockParams[velIndex*SIZE_BLOCKPARAMS + BlockParams::DVZ    ] = dvz_blockCell;
	 
	 // Calculate volume average for each cell in the block:
	 // NOTE: THIS NEEDS TO BE IMPROVED
	 for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
	    creal vx_blockCell = vx_block + (ic+0.5)*dvx_blockCell;
	    creal vy_blockCell = vy_block + (jc+0.5)*dvy_blockCell;
	    creal vz_blockCell = vz_block + (kc+0.5)*dvz_blockCell;
	    
	    grid[spaCellIndex]->cpu_avgs[velIndex*SIZE_VELBLOCK + kc*WID2+jc*WID+ic] =
	      calcPhaseSpaceDensity(x_cell+0.5*dx_cell,y_cell+0.5*dy_cell,z_cell+0.5*dz_cell,vx_blockCell,vy_blockCell,vz_blockCell)*dV;
	 }
	 
	 // Calculate velocity neighbour list entries for the block:
	 uint vxneg_nbr = iv-1;
	 uint vxpos_nbr = iv+1;
	 uint vyneg_nbr = jv-1;
	 uint vypos_nbr = jv+1;
	 uint vzneg_nbr = kv-1;
	 uint vzpos_nbr = kv+1;
	 
	 // --------- If the velocity block is a boundary block, raise the appropriate boundary flags.
	 // --------- YOU ALSO NEED TO SET SPATIAL BOUNDARY FLAGS HERE !!!
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE] = 0;
	 if (iv == 0) {              // -vx boundary
	    grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE]
	      = grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE] | NbrsVel::VX_NEG_BND;
	    vxneg_nbr = P::vxblocks_ini-1;
	 }
	 if (iv == P::vxblocks_ini-1) { // +vx boundary
	    grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE]
	      = grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE] | NbrsVel::VX_POS_BND;
	    vxpos_nbr = 0;
	 }
	 if (jv == 0) {              // -vy boundary
	    grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE]
	      = grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE] | NbrsVel::VY_NEG_BND;
	    vyneg_nbr = P::vyblocks_ini-1;
	 }
	 if (jv == P::vyblocks_ini-1) { // +vy boundary
	    grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE]
	      = grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE] | NbrsVel::VY_POS_BND;
	    vypos_nbr = 0;
	 }
	 if (kv == 0) {              // -vz boundary
	    grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE]
	      = grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE] | NbrsVel::VZ_NEG_BND;
	    vzneg_nbr = P::vzblocks_ini-1;
	 }
	 if (kv == P::vzblocks_ini-1) { // +vz boundary
	    grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE]
	      = grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::STATE] | NbrsVel::VZ_POS_BND;
	    vzpos_nbr = 0;
	 }
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::MYIND] = iniVelBlock + velIndex;
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VXNEG] = iniVelBlock + velblock(vxneg_nbr,jv,kv);
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VXPOS] = iniVelBlock + velblock(vxpos_nbr,jv,kv);
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VYNEG] = iniVelBlock + velblock(iv,vyneg_nbr,kv);
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VYPOS] = iniVelBlock + velblock(iv,vypos_nbr,kv);
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VZNEG] = iniVelBlock + velblock(iv,jv,vzneg_nbr);
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::VZPOS] = iniVelBlock + velblock(iv,jv,vzpos_nbr);
	 
	 // Calculate spatial grid neighbour list entries for the block:
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::XNEG] = grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::XNEG]*VELBLOCKS + velIndex;
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::XPOS] = grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::XPOS]*VELBLOCKS + velIndex;
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::YNEG] = grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::YNEG]*VELBLOCKS + velIndex;
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::YPOS] = grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::YPOS]*VELBLOCKS + velIndex;
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::ZNEG] = grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::ZNEG]*VELBLOCKS + velIndex;
	 grid[spaCellIndex]->cpu_nbrsVel[velIndex*SIZE_NBRS_VEL + NbrsVel::ZPOS] = grid[spaCellIndex]->cpu_nbrsSpa[NbrsSpa::ZPOS]*VELBLOCKS + velIndex;
      }
      iniVelBlock += grid[spaCellIndex]->N_blocks;
   }
   
   return true;
}

#endif
