#include <cstdlib>
#include <cuda_runtime_api.h>

#include "cudafuncs.h"
#include "logger.h"

using namespace std;

extern Logger logger;

bool deviceCreateArray(Real*& arrptr,const size_t& bytes) {
   cudaError_t error = cudaHostAlloc(reinterpret_cast<void**>(&arrptr),bytes,cudaHostAllocDefault);
   if (error == cudaSuccess) return true;
   logger << "\t CUDA error '" << cudaGetErrorString(error) << endl;
   return false;
}

bool deviceCreateArray(uint*& arrptr,const size_t& bytes) {
   cudaError_t error = cudaHostAlloc(reinterpret_cast<void**>(&arrptr),bytes,cudaHostAllocDefault);
   if (error == cudaSuccess) return true;
   logger << "\t CUDA error '" << cudaGetErrorString(error) << endl;
   return false;
}

bool deviceDeleteArray(Real*& arrptr) {
   cudaError_t error = cudaFreeHost(arrptr);
   if (error == cudaSuccess) {
      arrptr = NULL;
      return true;
   }
   return false;
}

bool deviceDeleteArray(uint*& arrptr) {
   cudaError_t error = cudaFreeHost(arrptr);
   if (error == cudaSuccess) {
      arrptr = NULL;
      return true;
   }
   return false;
}

void* gpuCreateArray(const std::string& name,const size_t& bytes) {
   void* ptr;
   logger << "Creating array '" << name << "', result = " << cudaGetErrorString( cudaMalloc(&ptr,bytes) ) << endl;
   return ptr;
}

void gpuDeleteArray(const std::string& name,void* ptr) {
   logger << "Deleting array '" << name << "', result = " << cudaGetErrorString( cudaFree(ptr) ) << endl;
}

void gpuCopyArray(const std::string& name,const size_t& bytes,void* cpuPtr,void* gpuPtr,const bool& cpuToGpu) {
   void* src;                // Pointer to source array
   void* dest;               // Pointer to destination array
   cudaMemcpyKind direction; // Copy direction
   
   if (cpuToGpu == true) { // Copying from host to device
      src = cpuPtr;
      dest = gpuPtr;
      direction = cudaMemcpyHostToDevice;
   } else { // Copying from device to host
      src = gpuPtr;
      dest = cpuPtr;
      direction = cudaMemcpyDeviceToHost;
   }
   
   logger << "Copying array '" << name << "', result = " << cudaGetErrorString( cudaMemcpy(dest,src,bytes,direction) ) << endl;
}

