#include "velocity_mesh_parameters.h"
#include <iostream>
#include <cstdlib>

#ifdef USE_CUDA
   #include "include/splitvector/splitvec.h"
   #include "cuda_runtime.h"
   #include "cuda_context.cuh"
#endif

// Pointers to MeshWrapper objects
static vmesh::MeshWrapper *meshWrapper;
#ifdef USE_CUDA
__device__ vmesh::MeshWrapper *meshWrapperDev;
#endif

__global__ void debug_kernel(
   const uint popID
   ) {
   vmesh::printVelocityMesh(0);
}

void vmesh::allocMeshWrapper() {
   // This is now allocated in unified memory
   meshWrapper = new vmesh::MeshWrapper();
}

vmesh::MeshWrapper* vmesh::host_getMeshWrapper() {
   return meshWrapper;
}

// This needs to be CUDA_HOSTDEV for compilation although it's called only from device side
#ifdef USE_CUDA
//#pragma hd_warning_disable // only applies to next function
#pragma nv_diag_suppress=20091
CUDA_HOSTDEV vmesh::MeshWrapper* vmesh::dev_getMeshWrapper() {
   return meshWrapperDev;
}
void vmesh::MeshWrapper::uploadMeshWrapper() {
   // Makes a copy of the meshWrapper (both the original and the new one now reside in unified memory).
   // The copy can then reside in device memory without causing page faults between host and device
   vmesh::MeshWrapper* MWdev = new vmesh::MeshWrapper(*meshWrapper);
   // Set the global symbol of meshWrapper (points to the copy)
   HANDLE_ERROR( cudaMemcpyToSymbol(meshWrapperDev, &MWdev, sizeof(vmesh::MeshWrapper*)) );
   MWdev->prefetchDevice();
   HANDLE_ERROR( cudaDeviceSynchronize() );
}
#endif

void vmesh::MeshWrapper::initVelocityMeshes(const uint nMeshes) {
   for (uint i=0; i<nMeshes; ++i) {
      vmesh::MeshParameters* vMesh = &(meshWrapper->velocityMeshes->at(i));
      // Duplicate definition of limits
      vMesh->meshMinLimits[0] = vMesh->meshLimits[0];
      vMesh->meshMinLimits[1] = vMesh->meshLimits[2];
      vMesh->meshMinLimits[2] = vMesh->meshLimits[4];
      vMesh->meshMaxLimits[0] = vMesh->meshLimits[1];
      vMesh->meshMaxLimits[1] = vMesh->meshLimits[3];
      vMesh->meshMaxLimits[2] = vMesh->meshLimits[5];

      // Calculate derived mesh parameters:
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
         = vMesh->gridLength[0]
         * vMesh->gridLength[1]
         * vMesh->gridLength[2];
      vMesh->initialized = true;
   }
#ifdef USE_CUDA
   // Now all velocity meshes have been initialized on host
   vmesh::MeshWrapper::uploadMeshWrapper();
#endif

   // vmesh::printVelocityMesh(0);
   // dim3 block(1,1,1);
   // debug_kernel<<<1, block, 0, 0>>> (0);
   // HANDLE_ERROR( cudaDeviceSynchronize() );
   return;
}
