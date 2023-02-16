#include "velocity_mesh_parameters.h"
#include <iostream>
#include <cstdlib>


#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#define CUDA_DEV __device__
#else
#define CUDA_HOSTDEV
#define CUDA_DEV
#endif

#ifdef USE_CUDA
   #include "include/splitvector/splitvec.h"
   #include "cuda_runtime.h"
   #include "cuda_context.cuh"
#endif

static vmesh::MeshWrapper *meshWrapper;
#ifdef USE_CUDA
__device__ vmesh::MeshWrapper *meshWrapperDev;
#endif

void vmesh::allocMeshWrapper() {
   meshWrapper = new vmesh::MeshWrapper();
   std::cerr<<"Allocated wrapper "<<meshWrapper<<std::endl;
}

vmesh::MeshWrapper* vmesh::getMeshWrapper() {
   return meshWrapper;
}

CUDA_DEV vmesh::MeshWrapper* vmesh::dev_getMeshWrapper() {
   return meshWrapperDev;
}

#ifdef USE_CUDA
void vmesh::MeshWrapper::uploadMeshWrapper() {
   // Create (copy) splitvector of meshes
   split::SplitVector<vmesh::MeshParameters> *velocityMeshesDev =
      new split::SplitVector<vmesh::MeshParameters>(*(meshWrapper->velocityMeshes));
   // Upload it to GPU
   split::SplitVector<vmesh::MeshParameters> *velocityMeshesDevUp = velocityMeshesDev->upload();
   // Create copy of meshWrapper
   vmesh::MeshWrapper meshWrapperTemp(*meshWrapper);
   // Make it point to splitvector in device memory
   printf("changing splitvector to point from 0x%llx ",meshWrapperTemp.velocityMeshes);
   meshWrapperTemp.velocityMeshes = velocityMeshesDevUp;
   printf("to 0x%llx \n",meshWrapperTemp.velocityMeshes);
   // Allocate device-side room for meshWrapper
   vmesh::MeshWrapper *meshWrapperDevUpload;
   HANDLE_ERROR( cudaMalloc((void **)&meshWrapperDevUpload, sizeof(vmesh::MeshWrapper)) );
   // Copy meshWrapper to device
   HANDLE_ERROR( cudaMemcpy(meshWrapperDevUpload, &meshWrapperTemp, sizeof(vmesh::MeshWrapper),cudaMemcpyHostToDevice) );
   //printf("Uploaded velocity mesh from 0x%llx to 0x%llx\n",&meshWrapperTemp,*meshWrapperDevUpload);

   HANDLE_ERROR( cudaMemcpyToSymbol(meshWrapperDev, meshWrapperDevUpload, sizeof(meshWrapperDev)) );
   // Don't bother freeing this small amount of device memory on exit.
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
      std::cerr<<"inited mesh "<<i<<" at "<<vMesh<<" for wrapper "<<meshWrapper<<std::endl;

      printf("Printout of velocity mesh %d\n",i);
   printf("Mesh limits \n");
   printf(" %f %f %f %f \n",vMesh->meshMinLimits[0],vMesh->meshLimits[0],vMesh->meshMaxLimits[0],vMesh->meshLimits[1]);
   printf(" %f %f %f %f \n",vMesh->meshMinLimits[1],vMesh->meshLimits[2],vMesh->meshMaxLimits[1],vMesh->meshLimits[3]);
   printf(" %f %f %f %f \n",vMesh->meshMinLimits[2],vMesh->meshLimits[4],vMesh->meshMaxLimits[2],vMesh->meshLimits[5]);
   printf("Derived mesh paramters \n");
   printf(" gridSize %f %f %f \n",vMesh->gridSize[0],vMesh->gridSize[1],vMesh->gridSize[2]);
   printf(" blockSize %f %f %f \n",vMesh->blockSize[0],vMesh->blockSize[1],vMesh->blockSize[2]);
   printf(" cellSize %f %f %f \n",vMesh->cellSize[0],vMesh->cellSize[1],vMesh->cellSize[2]);
   printf(" max velocity blocks %d \n",vMesh->max_velocity_blocks);
      
   }
#ifdef USE_CUDA
   vmesh::MeshWrapper::uploadMeshWrapper();
#endif
   return;
}

void vmesh::MeshWrapper::printVelocityMesh(const uint meshIndex) {
   vmesh::MeshParameters *vMesh = &(meshWrapper->velocityMeshes->at(meshIndex));
   std::cerr<<"printing mesh "<<meshIndex<<" at "<<vMesh<<" for wrapper "<<meshWrapper<<std::endl;

   printf("Printout of velocity mesh %d \n",meshIndex);
   printf("Mesh limits \n");
   printf(" %f %f %f %f \n",vMesh->meshMinLimits[0],vMesh->meshLimits[0],vMesh->meshMaxLimits[0],vMesh->meshLimits[1]);
   printf(" %f %f %f %f \n",vMesh->meshMinLimits[1],vMesh->meshLimits[2],vMesh->meshMaxLimits[1],vMesh->meshLimits[3]);
   printf(" %f %f %f %f \n",vMesh->meshMinLimits[2],vMesh->meshLimits[4],vMesh->meshMaxLimits[2],vMesh->meshLimits[5]);
   printf("Derived mesh paramters \n");
   printf(" gridSize %f %f %f \n",vMesh->gridSize[0],vMesh->gridSize[1],vMesh->gridSize[2]);
   printf(" blockSize %f %f %f \n",vMesh->blockSize[0],vMesh->blockSize[1],vMesh->blockSize[2]);
   printf(" cellSize %f %f %f \n",vMesh->cellSize[0],vMesh->cellSize[1],vMesh->cellSize[2]);
   printf(" max velocity blocks %d \n",vMesh->max_velocity_blocks);
}
