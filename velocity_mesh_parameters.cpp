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

void vmesh::allocMeshWrapper() {
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
   // Create (copy) splitvector of meshparams in unified memory
   split::SplitVector<vmesh::MeshParameters> *velocityMeshesDev =
      new split::SplitVector<vmesh::MeshParameters>(*(meshWrapper->velocityMeshes));
   // Upload it to GPU
   split::SplitVector<vmesh::MeshParameters> *velocityMeshesDevUp = velocityMeshesDev->upload();
   // Create temporary copy of meshWrapper
   vmesh::MeshWrapper meshWrapperTemp(*meshWrapper);
   // Make it point to splitvector in unified (device) memory
   meshWrapperTemp.velocityMeshes = velocityMeshesDevUp;
   // Allocate device-side room for meshWrapper
   vmesh::MeshWrapper *meshWrapperDevUpload;
   HANDLE_ERROR( cudaMalloc((void **)&meshWrapperDevUpload, sizeof(vmesh::MeshWrapper)) );
   // Copy meshWrapper to device
   HANDLE_ERROR( cudaMemcpy(meshWrapperDevUpload, &meshWrapperTemp, sizeof(vmesh::MeshWrapper),cudaMemcpyHostToDevice) );
   // Make __device__ pointer (symbol) point to this MeshWrapper
   HANDLE_ERROR( cudaMemcpyToSymbol(meshWrapperDev, meshWrapperDevUpload, sizeof(meshWrapperDev)) );
   // Don't bother freeing this small amount of device memory
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
   return;
}

CUDA_HOSTDEV void vmesh::MeshWrapper::printVelocityMesh(const uint meshIndex) {
   //vmesh::MeshParameters *vMesh = &(meshWrapper->velocityMeshes->at(meshIndex));
   //std::cerr<<"printing mesh "<<meshIndex<<" at "<<vMesh<<" for wrapper "<<meshWrapper<<std::endl;
   vmesh::MeshParameters *vMesh = &(getMeshWrapper()->velocityMeshes->at(meshIndex));

   printf("\nPrintout of velocity mesh %d \n",meshIndex);
   printf("Mesh size\n");
   printf(" %d %d %d \n",vMesh->gridLength[0],vMesh->gridLength[1],vMesh->gridLength[2]);
   printf("Block size (max reflevel %d)\n",vMesh->refLevelMaxAllowed);
   printf(" %d %d %d \n",vMesh->blockLength[0],vMesh->blockLength[1],vMesh->blockLength[2]);
   printf("Mesh limits \n");
   printf(" %f %f %f %f \n",vMesh->meshMinLimits[0],vMesh->meshLimits[0],vMesh->meshMaxLimits[0],vMesh->meshLimits[1]);
   printf(" %f %f %f %f \n",vMesh->meshMinLimits[1],vMesh->meshLimits[2],vMesh->meshMaxLimits[1],vMesh->meshLimits[3]);
   printf(" %f %f %f %f \n",vMesh->meshMinLimits[2],vMesh->meshLimits[4],vMesh->meshMaxLimits[2],vMesh->meshLimits[5]);
   printf("Derived mesh paramters \n");
   printf(" gridSize %f %f %f \n",vMesh->gridSize[0],vMesh->gridSize[1],vMesh->gridSize[2]);
   printf(" blockSize %f %f %f \n",vMesh->blockSize[0],vMesh->blockSize[1],vMesh->blockSize[2]);
   printf(" cellSize %f %f %f \n",vMesh->cellSize[0],vMesh->cellSize[1],vMesh->cellSize[2]);
   printf(" max velocity blocks %d \n\n",vMesh->max_velocity_blocks);
}
