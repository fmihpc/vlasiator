/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
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

#include <vector>
#include "../definitions.h"
#include "../spatial_cell_wrapper.hpp"
#include "../object_wrapper.h"
#include "../arch/gpu_base.hpp"
#include "gpu_trans_map_amr.hpp" // for loaning of allVmeshPointer
//#include <stdint.h>
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

//using namespace std;
using namespace spatial_cell;

/* Mini-kernel for evalutaing all blocks in all velocity meshes
 * finding the low and high corner velocities
 * comparing with the spatial cell size
 * and storing the largest allowed spatial dt for each cell
 *
 * @param allVmeshPointer Vector of pointers to velocitymeshes, used for gathering active blocks
 * @param dev_max_dt Buffer to store max allowed dt into (siz of nAllCells)
 * @param dev_dxdydz Buffer of cell spatial extents (size of 3*nAllCells)
 * @param nAllCells count of cells to read from allVmeshPointer
 */
__global__ void reduce_v_dt_kernel(
   const split::SplitVector<vmesh::VelocityMesh*> *allVmeshPointer,
   Real* dev_max_dt,
   const Real* dev_dxdydz,
   const uint nAllCells)
{
   const uint ti = threadIdx.x; // [0,GPUTHREADS*WARPSPERBLOCK)
   const uint blockSize = blockDim.x;
   const uint cellIndex = blockIdx.x; // userd for pointer to cell (or population)

   __shared__ Real smallest[GPUTHREADS*WARPSPERBLOCK]; //==blockSize
   smallest[ti] = numeric_limits<Real>::max();

   vmesh::VelocityMesh* thisVmesh = allVmeshPointer->at(cellIndex);
   uint thisVmeshSize = thisVmesh->size();
   Real blockInfo[6];
   const Real dx = dev_dxdydz[3*cellIndex + 0];
   const Real dy = dev_dxdydz[3*cellIndex + 1];
   const Real dz = dev_dxdydz[3*cellIndex + 2];
   const Real EPS = numeric_limits<Real>::min() * 1000;
   const Real HALF = 0.5;

   for (uint blockIndex = ti/2; blockIndex < thisVmeshSize; blockIndex += blockSize/2) {
      if (blockIndex < thisVmeshSize) {
         const vmesh::GlobalID GID = thisVmesh->getGlobalID(blockIndex);
         thisVmesh->getBlockInfo(GID,blockInfo); //This now calculates instead of reading from stored arrays
         // Indices 0-2 contain coordinates of the lower left corner.
         // Indices 3-5 contain the cell size.
         int i = (ti % 2) * (WID-1);
         // low and high corners, i.e., i == 0, i == WID - 1
         const Real Vx = blockInfo[0] + (i + HALF) * blockInfo[3] + EPS;
         const Real Vy = blockInfo[1] + (i + HALF) * blockInfo[4] + EPS;
         const Real Vz = blockInfo[2] + (i + HALF) * blockInfo[5] + EPS;
         smallest[ti] = min({dx / fabs(Vx), dy / fabs(Vy), dz / fabs(Vz), smallest[ti]});
      }
   }
   __syncthreads();
   // Now reduce for cell
   for (unsigned int s=blockSize/2; s>0; s>>=1) {
      if (ti < s) {
         smallest[ti] = min(smallest[ti],smallest[ti + s]);
      }
      __syncthreads();
   }
   // Kostis' suggestion: use two-stage warp votes to reduce
   // for (int offset = GPUTHREADS/2; offset > 0; offset >>=1){
   //    val = min(val,__shfl_down_sync(FULL_MASK, val, offset));
   // }
   // __syncthreads();
   if (ti==0) {
      dev_max_dt[cellIndex] = smallest[0];
   }   
}


void reduce_vlasov_dt(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                      const vector<CellID>& cells,
                      Real (&dtMaxLocal)[3]) {

   phiprof::Timer computeGpuTimestepTimer {"compute-vlasov-gpu-timestep"};
   // Does not use streams
   const uint nAllCells = cells.size();
   const uint nPOP = getObjectWrapper().particleSpecies.size();

   // Resize allVmeshPointer, one for each cell and each pop
   gpu_trans_allocate(nAllCells*nPOP,0,0,0);

   Real* host_max_dt;
   Real* host_dxdydz;
   Real* dev_max_dt;
   Real* dev_dxdydz;

   // Host memory will be pinned
   CHK_ERR( gpuMallocHost((void**)&host_max_dt, nAllCells*nPOP*sizeof(Real)) );
   CHK_ERR( gpuMallocHost((void**)&host_dxdydz, nAllCells*nPOP*3*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_max_dt, nAllCells*nPOP*sizeof(Real)) );
   CHK_ERR( gpuMalloc((void**)&dev_dxdydz, nAllCells*nPOP*3*sizeof(Real)) );

   // Gather vmeshes
   #pragma omp parallel for schedule(static)
   for(uint celli = 0; celli < nAllCells; celli++){
      SpatialCell* cell = mpiGrid[cells[celli]];
      cell->parameters[CellParams::MAXRDT] = numeric_limits<Real>::max();
      //cell->parameters[CellParams::MAXRDT] = numeric_limits<Real>::max();
      for (uint popID = 0; popID < nPOP; ++popID) {
         host_dxdydz[3*celli*nPOP + 3*popID + 0] = cell->parameters[CellParams::DX];
         host_dxdydz[3*celli*nPOP + 3*popID + 1] = cell->parameters[CellParams::DY];
         host_dxdydz[3*celli*nPOP + 3*popID + 2] = cell->parameters[CellParams::DZ];
         allVmeshPointer->at(celli*nPOP + popID) = cell->dev_get_velocity_mesh(popID); // GPU-side vmesh
      }
   }
   CHK_ERR( gpuMemcpy(dev_dxdydz, host_dxdydz, nAllCells*nPOP*3*sizeof(Real), gpuMemcpyHostToDevice) );
   allVmeshPointer->optimizeGPU(); // no stream given so blocking

   // Launch kernel gathering largest allowed dt for velocity
   reduce_v_dt_kernel<<<nAllCells, GPUTHREADS*WARPSPERBLOCK, 0, 0>>> (
      allVmeshPointer,
      dev_max_dt,
      dev_dxdydz,
      nAllCells*nPOP
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuMemcpy(host_max_dt, dev_max_dt, nAllCells*nPOP*sizeof(Real), gpuMemcpyDeviceToHost) );
   // CHK_ERR( gpuStreamSynchronize(bgStream) );

   // Real min_v_dt = numeric_limits<Real>::max();
   // #pragma omp parallel
   // {
   //    Real thread_min_v_dt = numeric_limits<Real>::max();
   //    #pragma omp for schedule(static)
   //    for(uint celli = 0; celli < nAllCells; celli++){
   //       SpatialCell* cell = mpiGrid[cells[celli]];
   //       cell->parameters[CellParams::MAXRDT] = host_max_dt[celli];
   //       thread_min_v_dt = thread_min_v_dt < host_max_dt[celli] ? thread_min_v_dt : host_max_dt[celli];
   //    }
   //    #pragma omp critical
   //    {
   //       min_v_dt = min_v_dt < thread_min_v_dt ? min_v_dt : thread_min_v_dt;
   //   }
   // }
   // dtMaxLocal[0]
   #pragma omp parallel for schedule(static)
   for(uint celli = 0; celli < nAllCells; celli++){
      SpatialCell* cell = mpiGrid[cells[celli]];
      for (uint popID = 0; popID < nPOP; ++popID) {
         cell->set_max_r_dt(popID, host_max_dt[celli*nPOP + popID]);
         cell->parameters[CellParams::MAXRDT] = min(cell->get_max_r_dt(popID), cell->parameters[CellParams::MAXRDT]);
      }
   }
   computeGpuTimestepTimer.stop();

   CHK_ERR( gpuFreeHost(host_max_dt) );
   CHK_ERR( gpuFreeHost(host_dxdydz) );
   CHK_ERR( gpuFree(dev_max_dt) );
   CHK_ERR( gpuFree(dev_dxdydz) );

   // GPUTODO thread this?
   phiprof::Timer computeRestTimestepTimer {"compute-vlasov-rest-timestep"};
   for (vector<CellID>::const_iterator cell_id = cells.begin(); cell_id != cells.end(); ++cell_id) {
      SpatialCell* cell = mpiGrid[*cell_id];

      if (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
          (cell->sysBoundaryLayer == 1 && cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY)) {
         // spatial fluxes computed also for L1 boundary cells
         dtMaxLocal[0] = min(dtMaxLocal[0], cell->parameters[CellParams::MAXRDT]);
      }

      if (cell->parameters[CellParams::MAXVDT] != 0 &&
          (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
           (P::vlasovAccelerateMaxwellianBoundaries && cell->sysBoundaryFlag == sysboundarytype::MAXWELLIAN))) {
         // acceleration only done on non-boundary cells
         dtMaxLocal[1] = min(dtMaxLocal[1], cell->parameters[CellParams::MAXVDT]);
      }
   }
   computeRestTimestepTimer.stop();
}
