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

template <reduceOp op, uint nDim, uint nReduStatic, typename Lambda, typename T>
__global__ static void 
__launch_bounds__(BLOCKSIZE_R)
reduction_kernel(Lambda loop_body, const T * __restrict__ init_val, T * __restrict__ rslt, const uint * __restrict__ lims, const uint n_total, const uint nReduDynamic)
{
    typedef cub::BlockReduce<T, BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1> BlockReduce;

    extern __shared__ typename BlockReduce::TempStorage temp_storage_dynamic[];
    __shared__ typename BlockReduce::TempStorage temp_storage_static[nReduStatic ? 2 * nReduStatic : 1];

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

    T data_static[nReduStatic ? 2 * nReduStatic : 1];
    T *aggregate_static = &data_static[0];

    T *aggregate_ptr;
    T *thread_data_ptr;

    typename BlockReduce::TempStorage *temp_storage_ptr;

    uint nReductions;

    if(nReduStatic == 0){
      aggregate_ptr = (T*)malloc((2 * nReduDynamic - 1) * sizeof(T));
      thread_data_ptr = &thread_data_ptr[nReduDynamic - 1];
      temp_storage_ptr = &temp_storage_dynamic[0];
      nReductions = nReduDynamic;
    }
    else{
      thread_data_ptr = &data_static[nReduDynamic - 1];
      aggregate_ptr = &aggregate_static[1];
      temp_storage_ptr = &temp_storage_static[1];
      nReductions = nReduStatic;
    }

    for(uint i = 0; i < nReductions; i++)
      thread_data_ptr[i] = init_val[i];

    // Evaluate the loop body
    lambda_eval(idx, &thread_data_ptr[0], loop_body);

    auto cub_call = [] __device__ (T &aggregate, typename BlockReduce::TempStorage &temp_storage, T &thread_data ) { 
      if(op == reduceOp::sum)
        aggregate = BlockReduce(temp_storage).Sum(thread_data); 
      else if(op == reduceOp::max)
        aggregate = BlockReduce(temp_storage).Reduce(thread_data, cub::Max()); 
      else if(op == reduceOp::min)
        aggregate = BlockReduce(temp_storage).Reduce(thread_data, cub::Min()); 
      else {}
    };

    auto store_call = [] __device__ (T &aggregate, T &rslt) { 
      if(op == reduceOp::sum)
        atomicAdd(&rslt, aggregate);
      else if(op == reduceOp::max)
        atomicMax(&rslt, aggregate);
      else if(op == reduceOp::min)
        atomicMin(&rslt, aggregate); 
      else
        printf("ERROR at %s:%d: Invalid reduction identifier \"op\".", __FILE__, __LINE__);
    };

    // Perform reductions
    // Compute the block-wide sum for thread0
    cub_call(aggregate_static[0], temp_storage_static[0], thread_data_ptr[0]);
    for(uint i = 0; i < nReductions - 1; i++)
      cub_call(aggregate_ptr[i], temp_storage_ptr[i], thread_data_ptr[i + 1]);

    // Store aggregate
    if(threadIdx.x == 0) 
      store_call(aggregate_static[0], rslt[0]);
      for(uint i = 0; i < nReductions - 1; i++)
        store_call(aggregate_ptr[i], rslt[i + 1]);
    
    if(nReduStatic == 0)
      free(aggregate_ptr);
}

template <reduceOp op, uint nReduStatic, uint nDim, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T *sum, const uint nReduDynamic) {

  uint nReductions;
  if(nReduStatic == 0)
    nReductions = nReduDynamic;
  else
    nReductions = nReduStatic;

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


  if(nReduStatic == 0) {
    uint cubTempStorageTypeSize = sizeof(cub::BlockReduce<T, BLOCKSIZE_R,cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1>);

    reduction_kernel<op, nDim, 0><<<gridsize, blocksize, (nReductions - 1) * cubTempStorageTypeSize>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, nReductions);
  }
  else{
    reduction_kernel<op, nDim, nReduStatic><<<gridsize, blocksize>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, nReductions);
  }
  
  CUDA_ERR(cudaStreamSynchronize(0));
  CUDA_ERR(cudaMemcpy(sum, d_buf, nReductions*sizeof(T), cudaMemcpyDeviceToHost));
  CUDA_ERR(cudaFree(d_buf));
  CUDA_ERR(cudaFree(d_const_buf));
  CUDA_ERR(cudaFree(d_limits));
}

#endif

template <reduceOp op, uint nDim, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T &sum) {
  constexpr uint nReductions = 1;
  arch::parallel_reduce<op, nReductions>(limits, loop_body, &sum, nReductions);
}

template <reduceOp op, uint nDim, uint nReductions, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T (&sum)[nReductions]) { 
  arch::parallel_reduce<op, nReductions>(limits, loop_body, &sum[0], nReductions);
}

template <reduceOp op, uint nDim, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, std::vector<T> &sum) {
  arch::parallel_reduce<op, 0>(limits, loop_body, sum.data(), sum.size());
}

}
