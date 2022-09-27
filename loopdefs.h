
namespace devices{

#ifndef USE_CUDA
#define CUDA_HOSTDEV
template <typename Lambda, typename T, uint nDim, uint nReductions>
inline void parallel_reduce(uint (&limits)[nDim], Lambda loop_body, T (&sum)[nReductions]) {

  if(nDim == 1){
    uint idx[1];
           
    # pragma omp for
    for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
      loop_body(idx, sum);
  }
  else if(nDim == 2){
    uint idx[2];
           
    # pragma omp for collapse(2)
    for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
      for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
        loop_body(idx, sum);
  }
  else if(nDim == 3){
    uint idx[3];
           
    # pragma omp for collapse(3)
    for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
      for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
        for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
          loop_body(idx, sum);  
  }
  else if(nDim == 4){
    uint idx[4];
           
    # pragma omp for collapse(4)
    for (idx[3] = 0; idx[3] < limits[3]; ++idx[3]) 
      for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
        for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
          for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
            loop_body(idx, sum);  
  }
  else{
    printf("ERROR in %s at line %d: This loop nest dimension is not supported!\n", __FILE__, __LINE__);
    exit(1);
  }
}

#else

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#define BLOCKSIZE_R 64

// The macros and functions that can be compiled with a C compiler
#define CUDA_ERR(err) (cuda_error(err, __FILE__, __LINE__))
  inline static void cuda_error(cudaError_t err, const char *file, int line) {
  	if (err != cudaSuccess) {
  		printf("\n\n%s in %s at line %d\n", cudaGetErrorString(err), file, line);
  		// exit(1);
  	}
}

template <int nDim, int nReductions, typename LambdaFun, typename T>
__global__ static void 
__launch_bounds__(BLOCKSIZE_R)
reduction_kernel(LambdaFun loop_fun, T * __restrict__ rslt, uint* lims, uint n_total)
{
    // Specialize BlockReduce for a 1D block of BLOCKSIZE_R * 1 * 1 threads on type T
#ifdef __CUDA_ARCH__
    typedef cub::BlockReduce<T, BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1, __CUDA_ARCH__> BlockReduce;
#else
    typedef cub::BlockReduce<T, BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1> BlockReduce;
#endif

    __shared__ typename BlockReduce::TempStorage temp_storage[nReductions];

    const uint idxGlob = ((blockIdx.x*blockDim.x)+threadIdx.x);

    if (idxGlob >= n_total) return;
    
    uint idx[nDim];
    switch (nDim)
    {
      case 4:
        idx[3] = (idxGlob / (lims[0] * lims[1] * lims[2])) % lims[3];
      case 3:
        idx[2] = (idxGlob / (lims[0] * lims[1])) % lims[2];
      case 2:
        idx[1] = (idxGlob / lims[0]) % lims[1];
      case 1:
        idx[0] = idxGlob % lims[0];  
    }

    T thread_data[nReductions];
    for(uint i = 0; i < nReductions; i++)
      thread_data[i] = 0;

    // Evaluate the loop body
    loop_fun(idx, thread_data);

    // Perform reductions
    // Compute the block-wide sum for thread0
    T aggregate[nReductions];
    for(uint i = 0; i < nReductions; i++)
      aggregate[i] = BlockReduce(temp_storage[i]).Sum(thread_data[i]); 

   //  // Store aggregate
    if(threadIdx.x == 0) 
      for(uint i = 0; i < nReductions; i++)
        atomicAdd(&rslt[i], aggregate[i]);
}


template <typename Lambda, typename T, uint nDim, uint nReductions>
__forceinline__ static void parallel_reduce(uint (&limits)[nDim], Lambda loop_body, T (&sum)[nReductions]) {

  uint n_total = 1;
  for(uint i = 0; i < nDim; i++)
    n_total *= limits[i];

  const uint blocksize = 64;
  const uint gridsize = (n_total - 1 + blocksize) / blocksize;

  T* d_buf;
  CUDA_ERR(cudaMalloc(&d_buf, nReductions*sizeof(T)));
  CUDA_ERR(cudaMemcpy(d_buf, sum, nReductions*sizeof(T), cudaMemcpyHostToDevice));
  
  uint* d_limits;
  CUDA_ERR(cudaMalloc(&d_limits, nDim*sizeof(uint)));
  CUDA_ERR(cudaMemcpy(d_limits, limits, nDim*sizeof(uint), cudaMemcpyHostToDevice));

  reduction_kernel<nDim,nReductions><<<gridsize, blocksize>>>(loop_body, d_buf, d_limits, n_total);
  
  CUDA_ERR(cudaStreamSynchronize(0));
  CUDA_ERR(cudaMemcpy(sum, d_buf, nReductions*sizeof(T), cudaMemcpyDeviceToHost));
  CUDA_ERR(cudaFree(d_buf));
  CUDA_ERR(cudaFree(d_limits));
}
#endif
}
