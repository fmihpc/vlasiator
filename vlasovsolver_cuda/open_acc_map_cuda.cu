#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "open_acc_map_h.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h"

/*
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ));
static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}
*/
__global__ void cudaFunction(int *b)
{
  int index = threadIdx.x + blockIdx.x*blockDim.x;
  printf("CUDA [%d]: \n", index);
  if(index<CUDASIZE)
  {
    b[index] = b[index]+index;
  }
}

void wrapper(int c)
{
  printf("STAGE 3\n");
  printf("c: %d\n", c);
  int b[CUDASIZE];
  int *dev_b;
  //HANDLE_ERROR( cudaMalloc((void**)&dev_b, CUDASIZE * sizeof(int)) );
  cudaMalloc((void**)&dev_b, CUDASIZE * sizeof(int));
  for(int a_c=0; a_c<CUDASIZE; a_c++)
  {
    b[a_c] = c-a_c;
  }
  printf("before: b: %d\n", b[0]);
	cudaMemcpy(dev_b, b, CUDASIZE*sizeof(int), cudaMemcpyHostToDevice);
  cudaFunction<<<BLOCKS, THREADS>>>(dev_b);
  cudaMemcpy(b, dev_b, CUDASIZE*sizeof(int), cudaMemcpyHostToDevice);
  printf("after: b: %d\n", b[CUDASIZE-1]);
  cudaFree(dev_b);
}
