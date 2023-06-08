#pragma once
#include "ionosphere_gpu_solver.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>
namespace ionogpu {

constexpr size_t warp_size = 32;

struct CudaStream {

   cudaStream_t stream;

   CudaStream() {
      cudaStreamCreate(&stream);
   }

   ~CudaStream() {
      cudaStreamDestroy(stream);
   }

   CudaStream(const CudaStream& other) = delete;
   CudaStream(CudaStream&& other) = delete;
   CudaStream& operator=(const CudaStream& other) = delete;
   CudaStream& operator=(CudaStream&& other) = delete;
   


};


template <typename T, size_t N = 0> class IonogpuMemoryAllocator {
   size_t size_;
   T* data_;

   size_t size_in_bytes_() const { return sizeof(T) * size_; }

   std::array<T*, N> pointers_to_data_;

public:
   struct PaddedVec {
      const std::vector<T>& vec;
      size_t minimun_size;
      bool zero_out_rest = false;

      size_t size() const noexcept { return std::max(vec.size(), minimun_size); }
   };

   IonogpuMemoryAllocator(const size_t size) : size_{size} { cudaMalloc((void**)&data_, size_in_bytes_()); }

   IonogpuMemoryAllocator(std::initializer_list<PaddedVec> padded_vecs)
       : size_{std::accumulate(
             padded_vecs.begin(), padded_vecs.end(), size_t{0},
             [](const auto s, const auto& padded_vec) noexcept -> size_t { return s + padded_vec.size(); })} {

      cudaMalloc((void**)&data_, size_in_bytes_());
      auto pointers_to_data = std::vector<T*>{};
      size_t current_place = 0;
      for (const auto& padded_vec : padded_vecs) {
         pointers_to_data.push_back(copy_vector_to_device_p(padded_vec.vec, current_place));
         if (padded_vec.zero_out_rest) {
            zero_out_memory(current_place + padded_vec.vec.size(), padded_vec.size() - padded_vec.vec.size());
         }
         current_place += std::max(padded_vec.vec.size(), padded_vec.minimun_size);
      }

      assert(N == padded_vecs.size());
      std::copy(pointers_to_data.begin(), pointers_to_data.end(), pointers_to_data_.begin());
   }

   ~IonogpuMemoryAllocator() {
      if (data_) {
         cudaFree(data_);
      }
   }

   // We want copying to be explicit
   // Copy constructor
   IonogpuMemoryAllocator(const IonogpuMemoryAllocator<T>& other) = delete;
   // Copy assignment
   IonogpuMemoryAllocator& operator=(const IonogpuMemoryAllocator<T>& other) = delete;

   // Move constructor
   IonogpuMemoryAllocator(IonogpuMemoryAllocator<T>&& other) noexcept {
      size_ = other.size_;
      data_ = other.data_;
      other.data_ = nullptr;
      pointers_to_data_ = std::move(other.pointers_to_data_);
   }
   // Move assignment
   IonogpuMemoryAllocator& operator=(IonogpuMemoryAllocator&& other) noexcept {
      size_ = other.size_;
      data_ = other.data_;
      other.data_ = nullptr;
      pointers_to_data_ = std::move(other.pointers_to_data_);
      return *this;
   }

   size_t size() const noexcept { return size_; }

   T* data() noexcept { return data_; }

   void copy_data_to_device(const T* data) { cudaMemcpy(data_, data, size_in_bytes_(), cudaMemcpyHostToDevice); }

   /**
    * Returns pointer to copied data
    */

   T* copy_data_to_device_p(const T* data, const size_t offset_in_elements, const size_t size_in_elements) {
      cudaMemcpy(data_ + offset_in_elements, data, size_in_elements * sizeof(T), cudaMemcpyHostToDevice);
      return data_ + offset_in_elements;
   }

   void copy_data_to_device(const T* data, const size_t offset_in_elements, const size_t size_in_elements) {
      cudaMemcpy(data_ + offset_in_elements, data, size_in_elements * sizeof(T), cudaMemcpyHostToDevice);
   }

   void copy_data_to_device(T const* const data, T* const place, const size_t size_in_elements) {
      copy_data_to_device(data, static_cast<size_t>(place - data_), size_in_elements);
   }

   void copy_vector_to_device(const std::vector<T>& v, const size_t offset_in_elements) {
      copy_data_to_device(v.data(), offset_in_elements, v.size());
   }

   void copy_vector_to_device(const std::vector<T>& v, T* const place) {
      copy_data_to_device(v.data(), static_cast<size_t>(place - data_), v.size());
   }

   T* copy_vector_to_device_p(const std::vector<T>& v, const size_t offset_in_elements) {
      return copy_data_to_device_p(v.data(), offset_in_elements, v.size());
   }

   void zero_out_memory() { cudaMemset(data_, 0, size_in_bytes_()); }

   void zero_out_memory(const size_t offset_in_elements, const size_t size_in_elements) {
      cudaMemset(data_ + offset_in_elements, 0, size_in_elements * sizeof(T));
   }

   void zero_out_memory(T* const place, const size_t size_in_elements) {
      zero_out_memory(static_cast<size_t>(place - data_), size_in_elements);
   }

   std::vector<T> copy_data_to_host_vector() const {
      auto temp_vector = std::vector<T>(size_);
      cudaMemcpy(temp_vector.data(), data_, size_in_bytes_(), cudaMemcpyDeviceToHost);
      return temp_vector;
   }

   std::vector<T> copy_data_to_host_vector(const size_t offset_in_elements, const size_t size_in_elements) const {
      auto temp_vector = std::vector<T>(size_in_elements);
      cudaMemcpy(temp_vector.data(), data_ + offset_in_elements, size_in_elements * sizeof(T), cudaMemcpyDeviceToHost);
      return temp_vector;
   }
   std::vector<T> copy_data_to_host_vector(const T* place, const size_t size_in_elements) const {
      return copy_data_to_host_vector(static_cast<size_t>(place - data_), size_in_elements);
   }
   std::array<T*, N> get_pointers_to_data() const noexcept { return pointers_to_data_; }

   /**
    * Constructor to create IonogpuMemoryAllocator from vector and order/multiply elements based on SparseMatrix
    * indecies So if we would calculate A * x we would get IonogpuMemoryAllocator with elements: flatten({
    *    {"elements of x for multiplying first row of A"},
    *    {"elements of x for multiplying second row of A"},
    *    ...,
    *    {"elements of x for multiplying n:th row of A"} })
    *
    */
};

template <typename T>
__global__ void matrixVectorProduct(T const* const M, T const* const v, T* const Mv, const size_t n) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   Mv[i] = T{0};
   for (size_t j = 0; j < n; ++j) {
      Mv[i] += M[i * n + j] * v[j];
   }
}

template <typename T> std::vector<T> matrixVectorProduct(const std::vector<T>& M, const std::vector<T>& v) {
   assert(M.size() == v.size() * v.size());
   const auto height = ((v.size() / warp_size) + 1) * warp_size;
   const auto width = v.size();

   const auto M_v_and_Mv_on_device =
       IonogpuMemoryAllocator<T, 3>{{M, width * height}, {v, width, true}, {{}, height, true}};

   const auto [M_device_p, v_device_p, Mv_device_p] = M_v_and_Mv_on_device.get_pointers_to_data();

   const auto blocks = height / warp_size;
   const auto threads_per_block = warp_size;

   matrixVectorProduct<T><<<blocks, threads_per_block>>>(M_device_p, v_device_p, Mv_device_p, width);
   return M_v_and_Mv_on_device.copy_data_to_host_vector(Mv_device_p, v.size());
}

template <typename T>
__global__ void preSparseMatrixVectorProduct(const size_t m, T const* const sparse_M, T const* const pre_b,
                                             T* const x) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   x[i] = T{0};
   for (size_t j = 0; j < m; ++j) {
      x[i] += sparse_M[i * m + j] * pre_b[i * m + j];
   }
}

template <typename T>
std::vector<T> preSparseMatrixVectorProduct(const size_t n, const size_t m, const std::vector<T>& sparse_M,
                                            const std::vector<size_t>& indecies, const std::vector<T>& b) {
   assert(sparse_M.size() == n * m);
   assert(indecies.size() == n * m);
   assert(b.size() == n);

   const auto height = ((n / warp_size) + 1) * warp_size;

   const auto pre_b = [=, &b] {
      auto temp = std::vector<T>(n * m);
      std::transform(indecies.begin(), indecies.end(), temp.begin(), [&b](const auto i) noexcept { return b[i]; });
      return temp;
   }();

   const auto sparse_M_pre_b_and_x_on_device =
       IonogpuMemoryAllocator<T, 3>{{sparse_M, height * m, true}, {pre_b, height * m}, {{}, height, true}};

   const auto [sparse_A_device_p, pre_b_device_p, x_device_p] = sparse_M_pre_b_and_x_on_device.get_pointers_to_data();

   const auto blocks = height / warp_size;
   const auto threads_per_block = warp_size;

   preSparseMatrixVectorProduct<T><<<blocks, threads_per_block>>>(m, sparse_A_device_p, pre_b_device_p, x_device_p);
   return sparse_M_pre_b_and_x_on_device.copy_data_to_host_vector(x_device_p, n);
};

template <typename T>
__global__ void sparseMatrixVectorProduct(const size_t m, T const* const sparse_M, size_t const* const indecies,
                                          T const* const b, T* x) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   T sum = T{0};
   for (size_t j = 0; j < m; ++j) {
      sum += sparse_M[i * m + j] * b[indecies[i * m + j]];
   }
   x[i] = sum;
}

template <typename T>
std::vector<T> sparseMatrixVectorProduct(const size_t n, const size_t m, const std::vector<T>& sparse_M,
                                         const std::vector<size_t>& indecies, const std::vector<T>& b) {
   assert(sparse_M.size() == n * m);
   assert(indecies.size() == n * m);
   assert(b.size() == n);

   const auto height = ((n / warp_size) + 1) * warp_size;

   const auto sparse_M_b_and_x_on_device =
       IonogpuMemoryAllocator<T, 3>{{sparse_M, height * m, true}, {b, height}, {{}, height}};

   const auto [sparse_A_device_p, b_device_p, x_device_p] = sparse_M_b_and_x_on_device.get_pointers_to_data();

   const auto indecies_on_device = IonogpuMemoryAllocator<size_t, 1>{{indecies, height * m}};

   const auto [indecies_device_p] = indecies_on_device.get_pointers_to_data();

   const auto blocks = height / warp_size;
   const auto threads_per_block = warp_size;

   sparseMatrixVectorProduct<T>
       <<<blocks, threads_per_block>>>(m, sparse_A_device_p, indecies_device_p, b_device_p, x_device_p);
   return sparse_M_b_and_x_on_device.copy_data_to_host_vector(x_device_p, n);
};

template <typename T> __global__ void vectorAddition(T const* const a, T const* const b, T* const result) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   result[i] = a[i] + b[i];
}

template <typename T> __global__ void vectorSubtraction(T const* const a, T const* const b, T* const result) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   result[i] = a[i] - b[i];
}
template <typename T>
__global__ void vectorElementwiseMultiplication(T const* const a, T const* const b, T* const result) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   result[i] = a[i] * b[i];
}

template <typename T> __global__ void vectorElementwiseDivision(T const* const a, T const* const b, T* const result) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   result[i] = a[i] / b[i];
}

template <typename T> std::vector<T> vectorAddition(const std::vector<T>& a, const std::vector<T>& b) {
   assert(a.size() == b.size());
   const auto height = ((a.size() / warp_size) + 1) * warp_size;
   const auto data_on_device = IonogpuMemoryAllocator<T, 3>{{a, height}, {b, height}, {{}, height}};
   const auto [a_device_p, b_device_p, result_device_p] = data_on_device.get_pointers_to_data();
   vectorAddition<<<height / warp_size, warp_size>>>(a_device_p, b_device_p, result_device_p);
   return data_on_device.copy_data_to_host_vector(result_device_p, a.size());
}

template <typename T> std::vector<T> vectorSubtraction(const std::vector<T>& a, const std::vector<T>& b) {
   assert(a.size() == b.size());
   const auto height = ((a.size() / warp_size) + 1) * warp_size;
   const auto data_on_device = IonogpuMemoryAllocator<T, 3>{{a, height}, {b, height}, {{}, height}};
   const auto [a_device_p, b_device_p, result_device_p] = data_on_device.get_pointers_to_data();
   vectorSubtraction<<<height / warp_size, warp_size>>>(a_device_p, b_device_p, result_device_p);
   return data_on_device.copy_data_to_host_vector(result_device_p, a.size());
}

template <typename T> std::vector<T> vectorElementwiseMultiplication(const std::vector<T>& a, const std::vector<T>& b) {
   assert(a.size() == b.size());
   const auto height = ((a.size() / warp_size) + 1) * warp_size;
   const auto data_on_device = IonogpuMemoryAllocator<T, 3>{{a, height}, {b, height}, {{}, height}};
   const auto [a_device_p, b_device_p, result_device_p] = data_on_device.get_pointers_to_data();
   vectorElementwiseMultiplication<<<height / warp_size, warp_size>>>(a_device_p, b_device_p, result_device_p);
   return data_on_device.copy_data_to_host_vector(result_device_p, a.size());
}

template <typename T> std::vector<T> vectorElementwiseDivision(const std::vector<T>& a, const std::vector<T>& b) {
   assert(a.size() == b.size());
   const auto height = ((a.size() / warp_size) + 1) * warp_size;
   const auto data_on_device = IonogpuMemoryAllocator<T, 3>{{a, height}, {b, height}, {{}, height}};
   const auto [a_device_p, b_device_p, result_device_p] = data_on_device.get_pointers_to_data();
   vectorElementwiseDivision<<<height / warp_size, warp_size>>>(a_device_p, b_device_p, result_device_p);
   return data_on_device.copy_data_to_host_vector(result_device_p, a.size());
}

// https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
// Modified to calculate dotproduct
template <typename T, unsigned int ThreadsPerBlock> __device__ void warpReduce(volatile T* sdata, unsigned int tid) {
   if constexpr (ThreadsPerBlock >= 64)
      sdata[tid] += sdata[tid + 32];
   if constexpr (ThreadsPerBlock >= 32)
      sdata[tid] += sdata[tid + 16];
   if constexpr (ThreadsPerBlock >= 16)
      sdata[tid] += sdata[tid + 8];
   if constexpr (ThreadsPerBlock >= 8)
      sdata[tid] += sdata[tid + 4];
   if constexpr (ThreadsPerBlock >= 4)
      sdata[tid] += sdata[tid + 2];
   if constexpr (ThreadsPerBlock >= 2)
      sdata[tid] += sdata[tid + 1];
}

// We have to have overloads for floats and doubles because of limitation of CUDA shared memory in template
// nvcc wont allow to have same dynamic shared memory array with the same name but different type (or alignment)
// so we can create overload that has different named shared dynamic memory (sdata_float and sdata_double)
// More info: https://stackoverflow.com/questions/27570552/templated-cuda-kernel-with-dynamic-shared-memory

// Each block will reduce 2 * ThreadsPerBlock amount of elements
template <unsigned int ThreadsPerBlock>
__global__ void reduce6([[maybe_unused]] float const* const g_idata1, [[maybe_unused]] float const* const g_idata2,
                        float* const partial_sums) {
   extern __shared__ float sdata_float[];
   const auto tid = threadIdx.x;
   constexpr auto elements_per_block = 2 * ThreadsPerBlock;
   const auto i = blockIdx.x * elements_per_block + tid;
   // const auto gridSize = elements_per_block * gridDim.x;
   // while (i < n) { sdata[tid] += g_idata1[i] * g_idata2[i] + g_idata1[i+ThreadsPerBlock] *
   // g_idata2[i+ThreadsPerBlock]; i += gridSize; }
   sdata_float[tid] = g_idata1[i] * g_idata2[i] + g_idata1[i + ThreadsPerBlock] * g_idata2[i + ThreadsPerBlock];
   __syncthreads();
   if constexpr (ThreadsPerBlock >= 512) {
      if (tid < 256) {
         sdata_float[tid] += sdata_float[tid + 256];
      }
      __syncthreads();
   }
   if constexpr (ThreadsPerBlock >= 256) {
      if (tid < 128) {
         sdata_float[tid] += sdata_float[tid + 128];
      }
      __syncthreads();
   }
   if constexpr (ThreadsPerBlock >= 128) {
      if (tid < 64) {
         sdata_float[tid] += sdata_float[tid + 64];
      }
      __syncthreads();
   }
   // At this point we don't have to sync threads beacause last 32 threads are from same warp so they are already synced
   if (tid < 32)
      warpReduce<float, ThreadsPerBlock>(sdata_float, tid);

   if (tid == 0)
      partial_sums[blockIdx.x] = sdata_float[0];
}

// Each block will reduce 2 * ThreadsPerBlock amount of elements
template <unsigned int ThreadsPerBlock>
__global__ void reduce6([[maybe_unused]] double const* const g_idata1, [[maybe_unused]] double const* const g_idata2,
                        double* const partial_sums) {
   extern __shared__ double sdata_double[];
   const auto tid = threadIdx.x;
   constexpr auto elements_per_block = 2 * ThreadsPerBlock;
   const auto i = blockIdx.x * elements_per_block + tid;
   // const auto gridSize = elements_per_block * gridDim.x;
   // while (i < n) { sdata[tid] += g_idata1[i] * g_idata2[i] + g_idata1[i+ThreadsPerBlock] *
   // g_idata2[i+ThreadsPerBlock]; i += gridSize; }
   sdata_double[tid] = g_idata1[i] * g_idata2[i] + g_idata1[i + ThreadsPerBlock] * g_idata2[i + ThreadsPerBlock];
   __syncthreads();
   if constexpr (ThreadsPerBlock >= 512) {
      if (tid < 256) {
         sdata_double[tid] += sdata_double[tid + 256];
      }
      __syncthreads();
   }
   if constexpr (ThreadsPerBlock >= 256) {
      if (tid < 128) {
         sdata_double[tid] += sdata_double[tid + 128];
      }
      __syncthreads();
   }
   if constexpr (ThreadsPerBlock >= 128) {
      if (tid < 64) {
         sdata_double[tid] += sdata_double[tid + 64];
      }
      __syncthreads();
   }
   // At this point we don't have to sync threads beacause last 32 threads are from same warp so they are already synced
   if (tid < 32)
      warpReduce<double, ThreadsPerBlock>(sdata_double, tid);

   if (tid == 0)
      partial_sums[blockIdx.x] = sdata_double[0];
}

template <typename T> __global__ void naive_sum(T* const p_device, const size_t length) {
   for (size_t i = 1; i < length; ++i) {
      p_device[0] += p_device[i];
   }
}

template <typename T> __global__ void zeroOutData(T* const data, const int n) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   if (i < n) {
      data[i] = 0;
   }
}

namespace dotProductConfig {
constexpr auto elements_per_block = 256;
constexpr auto threads_per_block = elements_per_block / 2;
} // namespace dotProductConfig

// when this function it is assumed that v_device_p and w_device_p have zero padding after n elements until padded_size
// [partial_sums_device_p, partial_sums_device_p + blocks] must be allocated
template <typename T>
T dotProduct(T const* const v_device_p, T const* const w_device_p, const size_t n, T* const partial_sums_device_p, const CudaStream& stream = CudaStream()) {
   const auto blocks = (n / dotProductConfig::elements_per_block) + 1;

   reduce6<dotProductConfig::threads_per_block>
       <<<blocks, dotProductConfig::threads_per_block, sizeof(T) * dotProductConfig::threads_per_block, stream.stream>>>(
           v_device_p, w_device_p, partial_sums_device_p);

   /*
       const auto partial_sums = [&] {
           auto temp = std::vector<T>(blocks);
           cudaMemcpy(temp.data(), partial_sums_device_p, blocks * sizeof(T), cudaMemcpyDeviceToHost);
           return temp;
       }(); */
   naive_sum<T><<<1, 1, 0, stream.stream>>>(partial_sums_device_p, blocks);
   T sum;
   cudaMemcpy(&sum, partial_sums_device_p, sizeof(T), cudaMemcpyDeviceToHost);
   return sum;
}

template <typename T> T dotProduct(const std::vector<T>& v, const std::vector<T>& w) {
   assert(v.size() == w.size());
   constexpr auto elements_per_block = dotProductConfig::elements_per_block;
   const auto blocks = ((v.size() / elements_per_block) + 1);
   const auto padded_size = blocks * elements_per_block;

   const auto data_device = IonogpuMemoryAllocator<T, 3>{{v, padded_size, true}, {w, padded_size, true}, {{}, blocks}};

   const auto [v_device_p, w_device_p, partial_sums_device_p] = data_device.get_pointers_to_data();

   return dotProduct<T>(v_device_p, w_device_p, v.size(), partial_sums_device_p);
}

// when this function it is assumed that v_device_p and w_device_p have zero padding after n elements until padded_size
// Look at dotProduct()
template <typename T> T vectorNormSquared(T const* const v_device_p, const size_t n, T* partial_sums_device_p, const CudaStream& stream = CudaStream()) {
   return dotProduct<T>(v_device_p, v_device_p, n, partial_sums_device_p, stream);
}

template <typename T> T vectorNormSquared(const std::vector<T>& v) {
   constexpr auto elements_per_block = dotProductConfig::elements_per_block;
   const auto blocks = ((v.size() / elements_per_block) + 1);
   const auto padded_size = blocks * elements_per_block;

   const auto data_device = IonogpuMemoryAllocator<T, 2>{{v, padded_size, true}, {{}, blocks}};

   const auto [v_device_p, partial_sums_device_p] = data_device.get_pointers_to_data();

   return dotProduct<T>(v_device_p, v_device_p, v.size(), partial_sums_device_p);
}

// scalar * v + w
template <typename T>
__global__ void multiplyVectorWithScalarAndAddItToAnotherVector(const T scalar, T const* const v_device_p,
                                                                T const* const w_device_p, T* const result_device_p) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   result_device_p[i] = scalar * v_device_p[i] + w_device_p[i];
}

// scalar * v + w
template <typename T>
std::vector<T> multiplyVectorWithScalarAndAddItToAnotherVector(const T scalar, const std::vector<T>& v,
                                                               const std::vector<T>& w) {
   assert(v.size() == w.size());
   const auto padded_size = (((v.size() / warp_size) + 1) * warp_size);

   const auto data_device = IonogpuMemoryAllocator<T, 3>{{v, padded_size}, {w, padded_size}, {{}, padded_size}};

   const auto [v_device_p, w_device_p, result_device_p] = data_device.get_pointers_to_data();

   multiplyVectorWithScalarAndAddItToAnotherVector<<<padded_size / warp_size, warp_size>>>(scalar, v_device_p,
                                                                                           w_device_p, result_device_p);
   return data_device.copy_data_to_host_vector(result_device_p, v.size());
}

// Modified Asolve from page 86 (110 pdf) of
// https://www.grad.hr/nastava/gs/prg/NumericalRecipesinC.pdf
// We assume that element (i, i) of sparse_A is stored at sparse_A[i * m]
template <typename T>
__global__ void Asolve_diagonal(size_t n, size_t m, T const* const sparse_A, T const* const b, T* const x) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;

   // This might be abel to be optimized away by setting padded part of sparse A to non_zero
   if (i < n) {
      x[i] = b[i] / sparse_A[i * m];
   } else {
      x[i] = 0;
   }
}

template <typename T> __global__ void copyData(T const* const source_device_p, T* const destination_device_p) {
   const auto i = blockDim.x * blockIdx.x + threadIdx.x;
   destination_device_p[i] = source_device_p[i];
}

template <typename T>
ReturnOfSparseBiCGCUDA<T> sparseBiCGCUDA(const size_t n, const size_t m, const std::vector<T>& sparse_A,
                                         const std::vector<T>& sparse_A_transposed, const std::vector<size_t>& indecies,
                                         const std::vector<T>& b,
                                         const ConfigurationForIonosphereGPUSolver<T>& config) {

   assert(sparse_A.size() == n * m);
   assert(sparse_A_transposed.size() == n * m);
   assert(indecies.size() == n * m);
   assert(b.size() == n);

   const auto height = ((n / warp_size) + 1) * warp_size;

#ifdef IONOGPU_STATIC_MEMORY_ALLOCATION
   static
#endif
       const auto indecies_on_device = IonogpuMemoryAllocator<size_t, 1>{{indecies, height * m, true}};
   const auto [indecies_device_p] = indecies_on_device.get_pointers_to_data();

   const auto blocks_for_dot_product = n / dotProductConfig::elements_per_block + 1;
   // Dot product needs more padding
   [[maybe_unused]] const auto padded_size_for_dot_product =
       blocks_for_dot_product * dotProductConfig::elements_per_block;

   const auto space_for_mask_gauge = (config.gauge == ionogpu::Gauge::mask) ? height : 0;

#ifdef IONOGPU_STATIC_MEMORY_ALLOCATION
   static
#endif
       auto data_on_device =
           IonogpuMemoryAllocator<T, 14>{{sparse_A, height * m, true},
                                         {sparse_A_transposed, height * m, true},
                                         {b, padded_size_for_dot_product, true},
                                         {/* x */ {}, height, true}, // Initial guess of x is zero
                                         {/* r */ {}, padded_size_for_dot_product, true},
                                         {/* rr */ {}, padded_size_for_dot_product, true},
                                         {/* z */ {}, padded_size_for_dot_product, true},
                                         {/* zz */ {}, height},
                                         {/* p */ {}, height},
                                         {/* pp */ {}, padded_size_for_dot_product, true},
                                         {/* best_solution */ {}, height},
                                         {/* partial sums for dot product */ {}, blocks_for_dot_product},
                                         {/* temp */ {}, height},
                                         {/* mask_gauge */ config.mask_gauge, space_for_mask_gauge, true}};

   const auto [sparse_A_device_p, sparse_A_transposed_device_p, b_device_p, x_device_p, r_device_p, rr_device_p,
               z_device_p, zz_device_p, p_device_p, pp_device_p, best_solution_device_p,
               partial_sums_for_dot_product_device_p, temp_device_p, mask_gauge_device_p] =
       data_on_device.get_pointers_to_data();

#ifdef IONOSPHERE_GPU_STATIC_MEMORY_ALLOCATION
   static bool initialized = false;
   if (initialized) {
      indices_on_device.zero_out_memory();
      data_on_device.zero_out_memory();
      indices_on_device.copy_vector_to_device(indecies, indecies_device_p);
      data_on_device.copy_vector_to_device(sparse_A, sparse_A_device_p);
      data_on_device.copy_vector_to_device(sparse_A_transposed, sparse_A_transposed_device_p);
      data_on_device.copy_vector_to_device(b, b_device_p);
      data_on_device.copy_vector_to_device(config.mask_gauge, mask_gauge_device_p);
   } else {
      initialized = true;
   }
#endif

   auto number_of_restarts = int{0};
   auto min_error = std::numeric_limits<T>::max();
   auto iteration = int{0};
   const auto stream = CudaStream();
   // This is part of the restart mechanism
   const auto bnrm = std::sqrt(vectorNormSquared<T>(b_device_p, n, partial_sums_for_dot_product_device_p, stream));


   do {
      copyData<<<height / warp_size, warp_size, 0, stream.stream>>>(best_solution_device_p, x_device_p);

      // *******************************************************************
      // From this point onwards comments will refere numecial recipies in C
      // *******************************************************************
      // Atime(n, x, r, 0);
      sparseMatrixVectorProduct<T><<<height / warp_size, warp_size, 0, stream.stream>>>(
          m, sparse_A_device_p, indecies_device_p, x_device_p, r_device_p);

      // r[j]=b[j]-r[j];
      vectorSubtraction<T><<<height / warp_size, warp_size, 0, stream.stream>>>(b_device_p, r_device_p, r_device_p);

      if (config.use_minimum_residual_variant) {
         // atimes(n,r,rr,0);
         sparseMatrixVectorProduct<T><<<height / warp_size, warp_size, 0, stream.stream>>>(
             m, sparse_A_device_p, indecies_device_p, r_device_p, rr_device_p);
      } else {
         // rr[j]=r[j];
         copyData<<<height / warp_size, warp_size, 0, stream.stream>>>(r_device_p, rr_device_p);
      }

      // Asumed we use itol == 1 convergence test
      // if (itol == 1) {
      //    bnrm=snrm(n,b,itol);

      auto Asolve = [&, sparse_A_device_p = sparse_A_device_p,
                     m](T const* const local_b_device_p, T* const local_x_device_p,
                        [[maybe_unused]] const bool transposed = false) -> void {
         switch (config.precondition) {
         case Precondition::none: {
            copyData<<<height / warp_size, warp_size, 0, stream.stream>>>(local_b_device_p, local_x_device_p);
         } break;
         case Precondition::diagonal: {
            Asolve_diagonal<T><<<height / warp_size, warp_size, 0, stream.stream>>>(n, m, sparse_A_device_p, local_b_device_p,
                                                                             local_x_device_p);
         } break;
         }
      };

      //    asolve(n,r,z,0);
      Asolve(r_device_p, z_device_p);

      // }
      // else if ...

      // while (*iter <= itmax) {
      // (*iter == 1)

      // This gets initialized in the loop. (First loop should get refactored out)
      T bkden;
      auto failcount = int{0};
      auto first_iteration = true;
      for (; iteration < config.max_iterations;) {
         // somehow iteration can not be in above for statement when linkin this with mpicc
         iteration += 1;

         // asolve(n,rr,zz,1);
         Asolve(rr_device_p, zz_device_p, true);
         // for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
         const auto bknum = dotProduct<T>(z_device_p, rr_device_p, n, partial_sums_for_dot_product_device_p, stream);

         if (first_iteration) {
            first_iteration = false;
            // p[j]=z[j];
            copyData<<<height / warp_size, warp_size, 0, stream.stream>>>(z_device_p, p_device_p);
            // pp[j]=zz[j];
            copyData<<<height / warp_size, warp_size, 0, stream.stream>>>(zz_device_p, pp_device_p);
         } else {
            const auto bk = bknum / bkden;
            // p[j]=bk*p[j]+z[j];
            // pp[j]=bk*pp[j]+zz[j];
            multiplyVectorWithScalarAndAddItToAnotherVector<T>
                <<<height / warp_size, warp_size, 0, stream.stream>>>(bk, p_device_p, z_device_p, p_device_p);
            multiplyVectorWithScalarAndAddItToAnotherVector<T>
                <<<height / warp_size, warp_size, 0, stream.stream>>>(bk, pp_device_p, zz_device_p, pp_device_p);
         }
         bkden = bknum;
         // atimes(n,p,z,0);
         sparseMatrixVectorProduct<T><<<height / warp_size, warp_size, 0, stream.stream>>>(
             m, sparse_A_device_p, indecies_device_p, p_device_p, z_device_p);

         // for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
         const auto akden = dotProduct<T>(z_device_p, pp_device_p, n, partial_sums_for_dot_product_device_p, stream);
         const auto ak = bknum / akden;

         // for (j=1;j<=n;j++) {
         //     x[j] += ak*p[j];
         //     r[j] -= ak*z[j];
         //     rr[j] -= ak*zz[j];
         // }
         multiplyVectorWithScalarAndAddItToAnotherVector<T>
             <<<height / warp_size, warp_size, 0, stream.stream>>>(ak, p_device_p, x_device_p, x_device_p);

         if (config.gauge == ionogpu::Gauge::mask) {
            vectorElementwiseMultiplication<T>
                <<<height / warp_size, warp_size, 0, stream.stream>>>(x_device_p, mask_gauge_device_p, x_device_p);
         }
         // *******************************************************************************
         // From this point onwards comments will not refere numecial recipies in C anymore
         // *******************************************************************************
         // Copying from orginal implementation instead of this:
         // multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size, 0, stream.stream>>>(-ak,
         // z_device_p, r_device_p, r_device_p); multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height /
         // warp_size, warp_size, 0, stream.stream>>>(-ak, zz_device_p, rr_device_p, rr_device_p);

         // We do this:

         sparseMatrixVectorProduct<T><<<height / warp_size, warp_size, 0, stream.stream>>>(
             m, sparse_A_device_p, indecies_device_p, x_device_p, temp_device_p);
         vectorSubtraction<T><<<height / warp_size, warp_size, 0, stream.stream>>>(b_device_p, temp_device_p, r_device_p);

         sparseMatrixVectorProduct<T><<<height / warp_size, warp_size, 0, stream.stream>>>(
             m, sparse_A_transposed_device_p, indecies_device_p, x_device_p, temp_device_p);

         vectorSubtraction<T><<<height / warp_size, warp_size, 0, stream.stream>>>(b_device_p, temp_device_p, rr_device_p);

         if (config.gauge == ionogpu::Gauge::mask) {
            vectorElementwiseMultiplication<T>
                <<<height / warp_size, warp_size, 0, stream.stream>>>(r_device_p, mask_gauge_device_p, r_device_p);
            vectorElementwiseMultiplication<T>
                <<<height / warp_size, warp_size, 0, stream.stream>>>(rr_device_p, mask_gauge_device_p, rr_device_p);
         }

         // We only support pole gauge at the moment:
         // if (config.gauge == Gauge::pole) {
         //     data_on_device.zero_out_memory(x_device_p, 1);
         // }

         // asolve(n,r,z,0);
         Asolve(r_device_p, z_device_p);

         // I think olderr in vlasiator solver is unused
         const auto error =
             std::sqrt(vectorNormSquared<T>(r_device_p, n, partial_sums_for_dot_product_device_p, stream)) / bnrm;

         if (error < min_error) {
            copyData<<<height / warp_size, warp_size, 0, stream.stream>>>(x_device_p, best_solution_device_p);
            min_error = error;
            failcount = 0;
         } else {
            copyData<<<height / warp_size, warp_size, 0, stream.stream>>>(best_solution_device_p, x_device_p);
            ++failcount;
         }

         if (min_error < config.relative_L2_convergence_threshold) {
            break;
         }
         if (failcount > config.max_failure_count || error > config.max_error_growth_factor * min_error) {
            ++number_of_restarts;
            break;
         }
      }
      // Restart mechanism
   } while (min_error > config.relative_L2_convergence_threshold && iteration < config.max_iterations);

   return ReturnOfSparseBiCGCUDA<T>{iteration, number_of_restarts, min_error,
                                    data_on_device.copy_data_to_host_vector(best_solution_device_p, n)};
};
} // namespace ionogpu
