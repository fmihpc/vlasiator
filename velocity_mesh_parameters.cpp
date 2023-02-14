#include "velocity_mesh_parameters.h"
#include <iostream>
#include <cstdlib>


#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#ifdef USE_CUDA
   #include "include/splitvector/splitvec.h"
   #include "cuda_runtime.h"
   #include "cuda_context.cuh"
   static bool meshWrapperIsOnDevice;
#endif

static vmesh::MeshWrapper *meshWrapper;
static vmesh::MeshWrapper *meshWrapperDev;

void vmesh::allocMeshWrapper() {
   meshWrapper = new vmesh::MeshWrapper();
   std::cerr<<"Allocated wrapper "<<meshWrapper<<std::endl;
   //return (meshWrapper!=NULL);
}

vmesh::MeshWrapper* vmesh::getMeshWrapper() {
   return meshWrapper;
}

CUDA_HOSTDEV vmesh::MeshWrapper* vmesh::dev_getMeshWrapper() {
   // Will this even work?
   return meshWrapperDev;
}

#ifdef USE_CUDA
void vmesh::MeshWrapper::uploadMeshWrapper() {
   // Create splitvector of meshes and upload it
   split::SplitVector<vmesh::MeshParameters> *velocityMeshesDev =
      new split::SplitVector<vmesh::MeshParameters>(*(meshWrapper->velocityMeshes));
   split::SplitVector<vmesh::MeshParameters> *velocityMeshesDevUp;
   velocityMeshesDevUp = velocityMeshesDev->upload();
   // Copy of meshWrapper which points to splitvector in device memory
   vmesh::MeshWrapper *meshWrapperTemp = new vmesh::MeshWrapper(*meshWrapper);
   //meshWrapperTemp = meshWrapper;
   meshWrapperTemp->velocityMeshes = velocityMeshesDevUp;
   //  Upload device-side wrapper

   HANDLE_ERROR( cudaMalloc((void **)&meshWrapperDev, sizeof(vmesh::MeshWrapper)) );
   HANDLE_ERROR( cudaMemcpyAsync(meshWrapperDev, &meshWrapperTemp, sizeof(vmesh::MeshWrapper),cudaMemcpyHostToDevice,cuda_getStream()) );
   // store flag indicating that required information exists on device
   meshWrapperIsOnDevice = true;
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
   meshWrapperIsOnDevice = false;
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

   //printf(" BLOCK LENGTH %d \n",(meshWrapper->velocityMeshes)[meshIndex].blockLength[0]);
}
