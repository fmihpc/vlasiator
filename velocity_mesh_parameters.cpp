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

#include "velocity_mesh_parameters.h"
#include <iostream>
#include <cstdlib>

#ifdef USE_GPU
   #include "include/splitvector/splitvec.h"
   #include "arch/gpu_base.hpp"
#endif

// Pointers to MeshWrapper objects
static vmesh::MeshWrapper *meshWrapper;

#ifdef USE_GPU
vmesh::MeshWrapper* MWdev;
std::array<vmesh::MeshParameters,MAX_VMESH_PARAMETERS_COUNT> *velocityMeshes_upload;

__global__ void debug_kernel(const uint popID) {
   vmesh::printVelocityMesh(0);
}
#endif

void vmesh::allocateMeshWrapper() {
   meshWrapper = new vmesh::MeshWrapper();
}

vmesh::MeshWrapper* vmesh::host_getMeshWrapper() {
   return meshWrapper;
}

#ifdef USE_GPU
//#pragma hd_warning_disable // only applies to next function
#pragma nv_diag_suppress=20091
// We record the set of constant memory symbols to static arrays. We are avoiding using objects with non-trivial Ctors
// to not risk they to be constructed after the static object that uses them. There will be an instance of the static
// objects per compilation unit, so the ordering is hard to control.
static vmesh::MeshWrapper** meshWrapperDevRegister[128] = {0};
vmesh::meshWrapperDevRegistor::meshWrapperDevRegistor(vmesh::MeshWrapper*& v) {
   for (size_t InstanceIdx = 0; InstanceIdx < sizeof(meshWrapperDevRegister) / sizeof(vmesh::MeshWrapper**);
        ++InstanceIdx) {
      if (auto*& slot = meshWrapperDevRegister[InstanceIdx]; !slot) {
         slot = &v;
         //printf("Got instance of device mesh wrapper handler %p (index %ld)\n", &v, InstanceIdx);
         return;
      }
   }
   assert(false && "Not enough slots to register mesh wrapper slots.");
}

void vmesh::MeshWrapper::uploadMeshWrapper() {
   // Store address to velocityMeshes array
   std::array<vmesh::MeshParameters,MAX_VMESH_PARAMETERS_COUNT> * temp = meshWrapper->velocityMeshes;
   // gpu-Malloc space on device, copy array contents
   CHK_ERR( gpuMalloc((void **)&velocityMeshes_upload, sizeof(std::array<vmesh::MeshParameters,MAX_VMESH_PARAMETERS_COUNT>)) );
   CHK_ERR( gpuMemcpy(velocityMeshes_upload, meshWrapper->velocityMeshes, sizeof(std::array<vmesh::MeshParameters,MAX_VMESH_PARAMETERS_COUNT>),gpuMemcpyHostToDevice) );
   // Make wrapper point to device-side array
   meshWrapper->velocityMeshes = velocityMeshes_upload;
   // Allocate and copy meshwrapper on device
   CHK_ERR( gpuMalloc((void **)&MWdev, sizeof(vmesh::MeshWrapper)) );
   CHK_ERR( gpuMemcpy(MWdev, meshWrapper, sizeof(vmesh::MeshWrapper),gpuMemcpyHostToDevice) );
   // Set the global symbol of meshWrapper
   int count=0;
   for (size_t InstanceIdx = 0; InstanceIdx < sizeof(meshWrapperDevRegister) / sizeof(vmesh::MeshWrapper**);
        ++InstanceIdx) {
      if (auto* slot = meshWrapperDevRegister[InstanceIdx]; slot) {
         //printf("Setting device mesh wrapper handler %p (index %ld)\n", slot, InstanceIdx);
         CHK_ERR( gpuMemcpyToSymbol(*slot, &MWdev, sizeof(vmesh::MeshWrapper*)) );
         count++;
      } else {
         break;
      }
   }
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   if(myRank == MASTER_RANK) {
      printf("Done setting all %d instances of device mesh wrapper handler!\n",count);
   }

   // Copy host-side address back
   meshWrapper->velocityMeshes = temp;
   // And sync
   CHK_ERR( gpuDeviceSynchronize() );
}
void vmesh::deallocateMeshWrapper() {
   CHK_ERR( gpuFree(velocityMeshes_upload) );
   CHK_ERR( gpuFree(MWdev) );
   CHK_ERR( gpuFree(meshWrapperDevInstance) );
   // And sync
   CHK_ERR( gpuDeviceSynchronize() );
}
#endif

void vmesh::MeshWrapper::initVelocityMeshes(const uint nMeshes) {
   // Verify lengths match?
   if (meshWrapper->velocityMeshesCreation->size() != nMeshes) {
      printf("Error! Initialized only %d velocity meshes out of %d created ones.\n",nMeshes,
             (int)meshWrapper->velocityMeshesCreation->size());
      abort();
   }
   // Create pointer to array of sufficient length
   meshWrapper->velocityMeshes = new std::array<vmesh::MeshParameters,MAX_VMESH_PARAMETERS_COUNT>;

   // Copy data in, also set auxiliary values
   for (uint i=0; i<nMeshes; ++i) {
      vmesh::MeshParameters* vMesh = &(meshWrapper->velocityMeshes->at(i));
      vmesh::MeshParameters* vMeshIn = &(meshWrapper->velocityMeshesCreation->at(i));

      // Limits
      vMesh->meshLimits[0] = vMeshIn->meshLimits[0];
      vMesh->meshLimits[1] = vMeshIn->meshLimits[1];
      vMesh->meshLimits[2] = vMeshIn->meshLimits[2];
      vMesh->meshLimits[3] = vMeshIn->meshLimits[3];
      vMesh->meshLimits[4] = vMeshIn->meshLimits[4];
      vMesh->meshLimits[5] = vMeshIn->meshLimits[5];
      // Grid length
      vMesh->gridLength[0] = vMeshIn->gridLength[0];
      vMesh->gridLength[1] = vMeshIn->gridLength[1];
      vMesh->gridLength[2] = vMeshIn->gridLength[2];
      // Block length
      vMesh->blockLength[0] = vMeshIn->blockLength[0];
      vMesh->blockLength[1] = vMeshIn->blockLength[1];
      vMesh->blockLength[2] = vMeshIn->blockLength[2];

      // Calculate derived mesh parameters:
      vMesh->meshMinLimits[0] = vMesh->meshLimits[0];
      vMesh->meshMinLimits[1] = vMesh->meshLimits[2];
      vMesh->meshMinLimits[2] = vMesh->meshLimits[4];
      vMesh->meshMaxLimits[0] = vMesh->meshLimits[1];
      vMesh->meshMaxLimits[1] = vMesh->meshLimits[3];
      vMesh->meshMaxLimits[2] = vMesh->meshLimits[5];

      vMesh->gridSize[0] = vMesh->meshMaxLimits[0] - vMesh->meshMinLimits[0];
      vMesh->gridSize[1] = vMesh->meshMaxLimits[1] - vMesh->meshMinLimits[1];
      vMesh->gridSize[2] = vMesh->meshMaxLimits[2] - vMesh->meshMinLimits[2];

      vMesh->blockSize[0] = vMesh->gridSize[0] / vMesh->gridLength[0];
      vMesh->blockSize[1] = vMesh->gridSize[1] / vMesh->gridLength[1];
      vMesh->blockSize[2] = vMesh->gridSize[2] / vMesh->gridLength[2];

      vMesh->cellSize[0] = vMesh->blockSize[0] / vMesh->blockLength[0];
      vMesh->cellSize[1] = vMesh->blockSize[1] / vMesh->blockLength[1];
      vMesh->cellSize[2] = vMesh->blockSize[2] / vMesh->blockLength[2];

      vMesh->max_velocity_blocks
         = vMeshIn->gridLength[0]
         * vMeshIn->gridLength[1]
         * vMeshIn->gridLength[2];
      vMesh->initialized = true;
   }
#ifdef USE_GPU
   // Now all velocity meshes have been initialized on host, into
   // the array. Now we need to upload a copy onto GPU.
   vmesh::MeshWrapper::uploadMeshWrapper();
   // printf("Host printout\n");
   // vmesh::printVelocityMesh(0);
   // printf("Device printout\n");
   // debug_kernel<<<1, 1, 0, 0>>> (0);
   // CHK_ERR( gpuDeviceSynchronize() );
#endif

   return;
}
