#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

namespace arch{

enum reduceOp { max, min, sum, prod };

template <typename Lambda, typename T>
CUDA_HOSTDEV inline static void lambda_eval(const uint (&idx)[1], T *thread_data, Lambda loop_body) { loop_body(idx[0], thread_data); }

template <typename Lambda, typename T>
CUDA_HOSTDEV inline static void lambda_eval(const uint (&idx)[2], T *thread_data, Lambda loop_body) { loop_body(idx[0], idx[1], thread_data); }

template <typename Lambda, typename T>
CUDA_HOSTDEV inline static void lambda_eval(const uint (&idx)[3], T *thread_data, Lambda loop_body) { loop_body(idx[0], idx[1], idx[2], thread_data); }

template <typename Lambda, typename T>
CUDA_HOSTDEV inline static void lambda_eval(const uint (&idx)[4], T *thread_data, Lambda loop_body) { loop_body(idx[0], idx[1], idx[2], idx[3], thread_data); }

#ifndef USE_CUDA
template <reduceOp op, uint nDim, uint nReductions, typename Lambda, typename T>
inline void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T (&sum)[nReductions]) {

  if(nDim == 1){
    uint idx[1];
           
    # pragma omp for
    for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
      lambda_eval(idx, sum, loop_body);
  }
  else if(nDim == 2){
    uint idx[2];
           
    # pragma omp for collapse(2)
    for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
      for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
        lambda_eval(idx, sum, loop_body);
  }
  else if(nDim == 3){
    uint idx[3];
           
    # pragma omp for collapse(3)
    for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
      for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
        for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
          lambda_eval(idx, sum, loop_body);
  }
  else if(nDim == 4){
    uint idx[4];
           
    # pragma omp for collapse(4)
    for (idx[3] = 0; idx[3] < limits[3]; ++idx[3]) 
      for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
        for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
          for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
            lambda_eval(idx, sum, loop_body);
  }
  else{
    printf("ERROR in %s at line %d: This loop nest dimension is not supported!\n", __FILE__, __LINE__);
    exit(1);
  }
}

#else

#define BLOCKSIZE_R 64

// The macros and functions that can be compiled with a C compiler
#define CUDA_ERR(err) (cuda_error(err, __FILE__, __LINE__))
  inline static void cuda_error(cudaError_t err, const char *file, int line) {
  	if (err != cudaSuccess) {
  		printf("\n\n%s in %s at line %d\n", cudaGetErrorString(err), file, line);
  		// exit(1);
  	}
}

__device__ __forceinline__ static void atomicMax(double *address, double val2) {
  unsigned long long ret = __double_as_longlong(*address);
  while(val2 > __longlong_as_double(ret)) {
    unsigned long long old = ret;
    if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old)
    break;
  }
}

__device__ __forceinline__ static void atomicMin(double *address, double val2) {
  unsigned long long ret = __double_as_longlong(*address);
  while(val2 < __longlong_as_double(ret)) {
    unsigned long long old = ret;
    if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old)
    break;
  }
}

template <reduceOp op, uint nDim, uint nReductions, typename Lambda, typename T>
__global__ static void 
__launch_bounds__(BLOCKSIZE_R)
reduction_kernel(Lambda loop_body, const T * __restrict__ init_val, T * __restrict__ rslt, const uint * __restrict__ lims, const uint n_total)
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
      thread_data[i] = init_val[i];

    // Evaluate the loop body
    lambda_eval(idx, &thread_data[0], loop_body);

    // Perform reductions
    T aggregate[nReductions];
    if(op == reduceOp::sum){
      // Compute the block-wide sum for thread0
      for(uint i = 0; i < nReductions; i++)
        aggregate[i] = BlockReduce(temp_storage[i]).Sum(thread_data[i]); 
  
      // Store aggregate
      if(threadIdx.x == 0) 
        for(uint i = 0; i < nReductions; i++)
          atomicAdd(&rslt[i], aggregate[i]);
    }
    else if(op == reduceOp::max) {
      // Compute the block-wide sum for thread0
      for(uint i = 0; i < nReductions; i++)
        aggregate[i] = BlockReduce(temp_storage[i]).Reduce(thread_data[i], cub::Max()); 
  
      // Store aggregate
      if(threadIdx.x == 0) 
        for(uint i = 0; i < nReductions; i++)
          atomicMax(&rslt[i], aggregate[i]);
    }
    else if(op == reduceOp::min) {
      // Compute the block-wide sum for thread0
      for(uint i = 0; i < nReductions; i++)
        aggregate[i] = BlockReduce(temp_storage[i]).Reduce(thread_data[i], cub::Min()); 
  
      // Store aggregate
      if(threadIdx.x == 0) 
        for(uint i = 0; i < nReductions; i++)
          atomicMin(&rslt[i], aggregate[i]);
    }
    else {
      printf("ERROR at %s:%d: Invalid reduction identifier \"op\".", __FILE__, __LINE__);
    }
}

template <reduceOp op, uint nDim, uint nReductions, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T (&sum)[nReductions]) {

  uint n_total = 1;
  for(uint i = 0; i < nDim; i++)
    n_total *= limits[i];

  const uint blocksize = 64;
  const uint gridsize = (n_total - 1 + blocksize) / blocksize;

  T* d_buf;
  CUDA_ERR(cudaMalloc(&d_buf, nReductions*sizeof(T)));
  CUDA_ERR(cudaMemcpy(d_buf, sum, nReductions*sizeof(T), cudaMemcpyHostToDevice));
  
  T* d_const_buf;
  CUDA_ERR(cudaMalloc(&d_const_buf, nReductions*sizeof(T)));
  CUDA_ERR(cudaMemcpy(d_const_buf, d_buf, nReductions*sizeof(T), cudaMemcpyDeviceToDevice));

  uint* d_limits;
  CUDA_ERR(cudaMalloc(&d_limits, nDim*sizeof(uint)));
  CUDA_ERR(cudaMemcpy(d_limits, limits, nDim*sizeof(uint), cudaMemcpyHostToDevice));

  reduction_kernel<op, nDim, nReductions><<<gridsize, blocksize>>>(loop_body, d_const_buf, d_buf, d_limits, n_total);
  
  CUDA_ERR(cudaStreamSynchronize(0));
  CUDA_ERR(cudaMemcpy(sum, d_buf, nReductions*sizeof(T), cudaMemcpyDeviceToHost));
  CUDA_ERR(cudaFree(d_buf));
  CUDA_ERR(cudaFree(d_const_buf));
  CUDA_ERR(cudaFree(d_limits));
}

#endif

template <reduceOp op, uint nDim, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T &sum) {
  T sum_array[1] = { sum };
  arch::parallel_reduce<op>(limits, loop_body, sum_array);
}

template <reduceOp op, uint nDim, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, std::unique_ptr<T[]> &sum) {
  T sum_array[1] = { sum };
  arch::parallel_reduce<op>(limits, loop_body, sum_array);
}

}
