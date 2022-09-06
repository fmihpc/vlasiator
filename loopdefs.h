
namespace devices{

#ifndef USE_CUDA
#define CUDA_HOSTDEV
template <typename Lambda, typename T, int nDim, int nReductions>
inline void parallel_reduce(int (&size)[nDim], Lambda loop_body, T (&sum)[nReductions]) {

  int idx[4];
         
  # pragma omp for
  for (idx[3] = 0; idx[3] < size[3]; ++idx[3]) 
    for (idx[2] = 0; idx[2] < size[2]; ++idx[2]) 
      for (idx[1] = 0; idx[1] < size[1]; ++idx[1]) 
        for (idx[0] = 0; idx[0] < size[0]; ++idx[0])
          loop_body(idx, sum);
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

template<int nDim>
struct Reduction
{
   int nSize[nDim];
   int nTot = 1;
};

template <int nDim, int nReductions, typename LambdaFun, typename T>
__global__ static void 
__launch_bounds__(BLOCKSIZE_R)
reduction_kernel(LambdaFun loop_fun, T * __restrict__ rslt, Reduction<nDim> lims)
{
    // Specialize BlockReduce for a 1D block of BLOCKSIZE_R * 1 * 1 threads on type T
#ifdef __CUDA_ARCH__
    typedef cub::BlockReduce<T, BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1, __CUDA_ARCH__> BlockReduce;
#else
    typedef cub::BlockReduce<T, BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1> BlockReduce;
#endif

    __shared__ typename BlockReduce::TempStorage temp_storage[nReductions];

    const int idxGlob = ((blockIdx.x*blockDim.x)+threadIdx.x);

    if (idxGlob >= lims.nTot) return;
    
    int idx[nDim];
    switch (nDim)
    {
      case 4:
        idx[3] = (idxGlob / (lims.nSize[0] * lims.nSize[1] * lims.nSize[2])) % lims.nSize[3];
      case 3:
        idx[2] = (idxGlob / (lims.nSize[0] * lims.nSize[1])) % lims.nSize[2];
      case 2:
        idx[1] = (idxGlob / lims.nSize[0]) % lims.nSize[1];
      case 1:
        idx[0] = idxGlob % lims.nSize[0];  
    }

    T thread_data[nReductions];
    for(int i = 0; i < nReductions; i++)
      thread_data[i] = 0;

    // Evaluate the loop body
    loop_fun(idx, thread_data);

    // Perform reductions
    // Compute the block-wide sum for thread0
    T aggregate[nReductions];
    for(int i = 0; i < nReductions; i++)
      aggregate[i] = BlockReduce(temp_storage[i]).Sum(thread_data[i]); 

   //  // Store aggregate
    if(threadIdx.x == 0) 
      for(int i = 0; i < nReductions; i++)
        atomicAdd(&rslt[i], aggregate[i]);
}


template <typename Lambda, typename T, int nDim, int nReductions>
__forceinline__ static void parallel_reduce(int (&size)[nDim], Lambda loop_body, T (&sum)[nReductions]) {

  Reduction<nDim> lims;
  for(int i = 0; i < nDim; i++){
    lims.nSize[i] = size[i];
    lims.nTot *= size[i];
  }
  
  const int blocksize = 64;
  const int gridsize = (lims.nTot- 1 + blocksize) / blocksize;

  T* buf;
  CUDA_ERR(cudaMalloc(&buf, nReductions*sizeof(T)));
  CUDA_ERR(cudaMemcpy(buf, sum, nReductions*sizeof(T), cudaMemcpyHostToDevice));

  reduction_kernel<nDim,nReductions><<<gridsize, blocksize>>>(loop_body, buf, lims);
  
  CUDA_ERR(cudaStreamSynchronize(0));
  CUDA_ERR(cudaMemcpy(sum, buf, nReductions*sizeof(T), cudaMemcpyDeviceToHost));
  CUDA_ERR(cudaFree(buf));
}
#endif
}
