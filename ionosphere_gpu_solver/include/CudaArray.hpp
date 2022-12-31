#pragma once
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <vector>
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <cmath>
#include "timer.hpp"
namespace ionogpu {

constexpr size_t warp_size = 32;
    
template <typename T, size_t N = 0>
class CudaArray {
    size_t size_;
    T* data_;

    size_t size_in_bytes_() const {
        return sizeof(T) * size_;
    }

    std::array<T*, N> pointers_to_data_;
    
public:

    struct PaddedVec {
        const std::vector<T>& vec;
        size_t minimun_size;        
        bool zero_out_rest = false;

        size_t size() const noexcept {
            return std::max(vec.size(), minimun_size);
        }
    };

    CudaArray(const size_t size) : size_{size} {
        cudaMalloc((void**)&data_, size_in_bytes_());
    }

    CudaArray(std::initializer_list<PaddedVec> padded_vecs)
     : size_ {std::accumulate(padded_vecs.begin(), padded_vecs.end(), size_t{ 0 },
        [] (const auto s, const auto& padded_vec) noexcept -> size_t {
            return s + padded_vec.size();
        }
        )} {

        cudaMalloc((void**)&data_, size_in_bytes_());
        auto pointers_to_data = std::vector<T*> {};
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

    ~CudaArray() {
        if (data_) {
            cudaFree(data_);
        }
    }

    // We want copying to be explicit
    // Copy constructor
    CudaArray(const CudaArray<T>& other) = delete;
    // Copy assignment
    CudaArray& operator=(const CudaArray<T>& other) = delete;
    
    // Move constructor
    CudaArray(CudaArray<T>&& other) noexcept {
        size_ = other.size_; 
        data_ = other.data_;
        other.data_ = nullptr;
        pointers_to_data_ = std::move(other.pointers_to_data_);
        
    }
    // Move assignment
    CudaArray& operator=(CudaArray&& other) noexcept {
        size_ = other.size_; 
        data_ = other.data_;
        other.data_ = nullptr;
        pointers_to_data_ = std::move(other.pointers_to_data_);
        return *this;
    }
    
    size_t size() const noexcept {
        return size_;
    }

    T* data() noexcept {
        return data_;
    }
    
    void copy_data_to_device(const T* data) {
        cudaMemcpy(data_, data, size_in_bytes_(), cudaMemcpyHostToDevice); 
    }

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

    void copy_data_to_device(T const * const data, T* const place, const size_t size_in_elements) {
        copy_data_to_device(data, static_cast<size_t>(place - data_), size_in_elements);
    }

    void copy_vector_to_device(const std::vector<T>& v, const size_t offset_in_elements) {
        copy_data_to_device(v.data(), offset_in_elements, v.size());
    }

    T* copy_vector_to_device_p(const std::vector<T>& v, const size_t offset_in_elements) {
        return copy_data_to_device_p(v.data(), offset_in_elements, v.size());
    }

    void zero_out_memory() {
        cudaMemset(data_, 0, size_in_bytes_());
    }

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
    std::array<T*, N> get_pointers_to_data() const noexcept {
        return pointers_to_data_;
    }

    /**
    * Constructor to create CudaArray from vector and order/multiply elements based on SparseMatrix indecies
    * So if we would calculate A * x we would get CudaArray with elements:
    * flatten({
    *    {"elements of x for multiplying first row of A"},
    *    {"elements of x for multiplying second row of A"},
    *    ...,
    *    {"elements of x for multiplying n:th row of A"} })
    *
    */
};



template <typename T>
__global__ void matrixVectorProduct(T const * const  M, T const * const v, T* const Mv, const size_t n) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    Mv[i] = T{ 0 };
    for (size_t j = 0; j < n; ++j) {
        Mv[i] += M[i * n + j] * v[j];
    }
}

template <typename T>
std::vector<T> matrixVectorProduct(const std::vector<T>& M, const std::vector<T>& v) { 
    assert(M.size() == v.size() * v.size());
    const auto height = ((v.size() / warp_size) + 1) * warp_size;
    const auto width = v.size();

    const auto M_v_and_Mv_on_device = CudaArray<T, 3> {
        {M, width * height},
        {v, width, true},
        {{}, height, true}
    };

    const auto [M_device_p, v_device_p, Mv_device_p] = M_v_and_Mv_on_device.get_pointers_to_data();

    const auto blocks = height / warp_size; 
    const auto threads_per_block = warp_size;
 
    matrixVectorProduct<T><<<blocks, threads_per_block>>>(M_device_p, v_device_p, Mv_device_p, width); 
    return M_v_and_Mv_on_device.copy_data_to_host_vector(Mv_device_p, v.size());  
    
}


template <typename T>
__global__ void preSparseMatrixVectorProduct(
    const size_t m,
    T const * const sparse_M,
    T const * const pre_b,
    T* const x
) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    x[i] = T{ 0 };
    for (size_t j = 0; j < m; ++j) {
        x[i] += sparse_M[i * m + j] * pre_b[i * m + j];
    }
}

template<typename T>
std::vector<T> preSparseMatrixVectorProduct(
    const size_t n, const size_t m,    
    const std::vector<T>& sparse_M,
    const std::vector<size_t>& indecies,
    const std::vector<T>& b
) {
    assert(sparse_M.size() == n * m);
    assert(indecies.size() == n * m);
    assert(b.size() == n);

    const auto height = ((n / warp_size) + 1) * warp_size;

    const auto pre_b = [=, &b] {
        auto temp  = std::vector<T>(n * m);
        std::transform(indecies.begin(), indecies.end(), temp.begin(), [&b](const auto i) noexcept {return b[i];});
        return temp;
    }();

    
    const auto sparse_M_pre_b_and_x_on_device = CudaArray<T, 3> {
        {sparse_M, height * m, true},
        {pre_b, height * m},
        {{}, height, true}
    };

    const auto [sparse_A_device_p, pre_b_device_p, x_device_p] = sparse_M_pre_b_and_x_on_device.get_pointers_to_data();

    const auto blocks = height / warp_size;
    const auto threads_per_block = warp_size;

    preSparseMatrixVectorProduct<T><<<blocks, threads_per_block>>>(m, sparse_A_device_p, pre_b_device_p, x_device_p);
    return sparse_M_pre_b_and_x_on_device.copy_data_to_host_vector(x_device_p, n);
};


template <typename T>
__global__ void sparseMatrixVectorProduct(
    const size_t m,
    T const * const sparse_M,
    size_t const * const indecies,
    T const * const b,
    T* x
) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    x[i] = T{ 0 };
    for (size_t j = 0; j < m; ++j) {
        x[i] += sparse_M[i * m + j] * b[indecies[i * m + j]];
    }
}

template<typename T>
std::vector<T> sparseMatrixVectorProduct(
    const size_t n, const size_t m,    
    const std::vector<T>& sparse_M,
    const std::vector<size_t>& indecies,
    const std::vector<T>& b
) {
    assert(sparse_M.size() == n * m);
    assert(indecies.size() == n * m);
    assert(b.size() == n);

    const auto height = ((n / warp_size) + 1) * warp_size;


    const auto sparse_M_b_and_x_on_device = CudaArray<T, 3> {
        {sparse_M, height * m, true},
        {b, height},
        {{}, height}
    };

    const auto [sparse_A_device_p, b_device_p, x_device_p] = sparse_M_b_and_x_on_device.get_pointers_to_data();

    const auto indecies_on_device = CudaArray<size_t, 1> {
        {indecies, height * m}
    };

    const auto [indecies_device_p] = indecies_on_device.get_pointers_to_data();

    const auto blocks = height / warp_size;
    const auto threads_per_block = warp_size;

    sparseMatrixVectorProduct<T><<<blocks, threads_per_block>>>(m, sparse_A_device_p, indecies_device_p, b_device_p, x_device_p);
    return sparse_M_b_and_x_on_device.copy_data_to_host_vector(x_device_p, n);
};

template<typename T>
__global__ void vectorAddition(T const * const a, T const * const b, T * const result) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    result[i] = a[i] + b[i];
}

template<typename T> 
__global__ void vectorSubtraction(T const * const a, T const * const b, T * const result) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    result[i] = a[i] - b[i];
}
template<typename T>
__global__ void vectorElementwiseMultiplication(T const * const a, T const * const b, T * const result) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    result[i] = a[i] * b[i];
}

template<typename T>
__global__ void vectorElementwiseDivision(T const * const a, T const * const b, T * const result) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    result[i] = a[i] / b[i];
}

template<typename T>
std::vector<T> vectorAddition(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    const auto height = ((a.size() / warp_size) + 1) * warp_size;
    const auto data_on_device = CudaArray<T, 3> {
        {a, height},
        {b, height},
        {{}, height}
    };
    const auto [a_device_p, b_device_p, result_device_p] = data_on_device.get_pointers_to_data();
    vectorAddition<<<height / warp_size, warp_size>>>(a_device_p, b_device_p, result_device_p);
    return data_on_device.copy_data_to_host_vector(result_device_p, a.size());
}

template<typename T>
std::vector<T> vectorSubtraction(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    const auto height = ((a.size() / warp_size) + 1) * warp_size;
    const auto data_on_device = CudaArray<T, 3> {
        {a, height},
        {b, height},
        {{}, height}
    };
    const auto [a_device_p, b_device_p, result_device_p] = data_on_device.get_pointers_to_data();
    vectorSubtraction<<<height / warp_size, warp_size>>>(a_device_p, b_device_p, result_device_p);
    return data_on_device.copy_data_to_host_vector(result_device_p, a.size());
}

template<typename T>
std::vector<T> vectorElementwiseMultiplication(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    const auto height = ((a.size() / warp_size) + 1) * warp_size;
    const auto data_on_device = CudaArray<T, 3> {
        {a, height},
        {b, height},
        {{}, height}
    };
    const auto [a_device_p, b_device_p, result_device_p] = data_on_device.get_pointers_to_data();
    vectorElementwiseMultiplication<<<height / warp_size, warp_size>>>(a_device_p, b_device_p, result_device_p);
    return data_on_device.copy_data_to_host_vector(result_device_p, a.size());
}

template<typename T>
std::vector<T> vectorElementwiseDivision(const std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());
    const auto height = ((a.size() / warp_size) + 1) * warp_size;
    const auto data_on_device = CudaArray<T, 3> {
        {a, height},
        {b, height},
        {{}, height}
    };
    const auto [a_device_p, b_device_p, result_device_p] = data_on_device.get_pointers_to_data();
    vectorElementwiseDivision<<<height / warp_size, warp_size>>>(a_device_p, b_device_p, result_device_p);
    return data_on_device.copy_data_to_host_vector(result_device_p, a.size());
}

// https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
// Modified to calculate dotproduct
template <typename T, unsigned int ThreadsPerBlock>
__device__ void warpReduce(volatile T *sdata, unsigned int tid) {
    if constexpr (ThreadsPerBlock >= 64) sdata[tid] += sdata[tid + 32];
    if constexpr (ThreadsPerBlock >= 32) sdata[tid] += sdata[tid + 16];
    if constexpr (ThreadsPerBlock >= 16) sdata[tid] += sdata[tid + 8];
    if constexpr (ThreadsPerBlock >= 8) sdata[tid] += sdata[tid + 4];
    if constexpr (ThreadsPerBlock >= 4) sdata[tid] += sdata[tid + 2];
    if constexpr (ThreadsPerBlock >= 2) sdata[tid] += sdata[tid + 1];
}

// Each block will reduce 2 * ThreadsPerBlock amount of elements 
template <typename T, unsigned int ThreadsPerBlock>
__global__ void reduce6([[maybe_unused]]T const * const g_idata1, [[maybe_unused]]T const * const g_idata2, T * const partial_sums) {
    extern __shared__ T sdata[];
    const auto tid = threadIdx.x;
    constexpr auto elements_per_block = 2 * ThreadsPerBlock;
    const auto i = blockIdx.x * elements_per_block + tid;
    //const auto gridSize = elements_per_block * gridDim.x;
    //while (i < n) { sdata[tid] += g_idata1[i] * g_idata2[i] + g_idata1[i+ThreadsPerBlock] * g_idata2[i+ThreadsPerBlock]; i += gridSize; }
    sdata[tid] = g_idata1[i] * g_idata2[i] + g_idata1[i + ThreadsPerBlock] * g_idata2[i + ThreadsPerBlock];
    __syncthreads();
    if constexpr (ThreadsPerBlock >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if constexpr (ThreadsPerBlock >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if constexpr (ThreadsPerBlock >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
    // At this point we don't have to sync threads beacause last 32 threads are from same warp so they are already synced
    if (tid < 32) warpReduce<T, ThreadsPerBlock>(sdata, tid);

    if (tid == 0) partial_sums[blockIdx.x] = sdata[0]; 
}

template <typename T>
__global__ void zeroOutData(T * const data, const int n) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) {
        data[i] = 0;
    }
    
}

namespace dotProductConfig {
    constexpr auto elements_per_block = 256;
    constexpr auto threads_per_block = elements_per_block / 2;
}



// when this function it is assumed that v_device_p and w_device_p have zero padding after n elements until padded_size
// [partial_sums_device_p, partial_sums_device_p + blocks] must be allocated
template<typename T>
T dotProduct(T const * const v_device_p, T const * const w_device_p, const size_t n, T * const partial_sums_device_p) {
    const auto blocks = (n / dotProductConfig::elements_per_block) + 1;
    
    reduce6<T, dotProductConfig::threads_per_block>
        <<<blocks, dotProductConfig::threads_per_block, sizeof(T) * dotProductConfig::threads_per_block>>>(v_device_p, w_device_p, partial_sums_device_p);
    

    const auto partial_sums = [&] {
        auto temp = std::vector<T>(blocks);
        cudaMemcpy(temp.data(), partial_sums_device_p, blocks * sizeof(T), cudaMemcpyDeviceToHost);
        return temp;
    }();

    return std::accumulate(partial_sums.begin(), partial_sums.end(), T{ 0 });
}

template<typename T>
T dotProduct(const std::vector<T>& v, const std::vector<T>& w) {
    assert(v.size() == w.size());
    constexpr auto elements_per_block = dotProductConfig::elements_per_block; 
    const auto blocks = ((v.size() / elements_per_block) + 1);
    const auto padded_size = blocks * elements_per_block;

    const auto data_device = CudaArray<T, 3> {
        {v, padded_size, true},
        {w, padded_size, true},
        {{}, blocks}
    };

    const auto [v_device_p, w_device_p, partial_sums_device_p] = data_device.get_pointers_to_data();

    return dotProduct<T>(v_device_p, w_device_p, v.size(), partial_sums_device_p);
}

// when this function it is assumed that v_device_p and w_device_p have zero padding after n elements until padded_size
// Look at dotProduct()
template<typename T>
T vectorNormSquared(T const * const v_device_p, const size_t n, T * partial_sums_device_p) {
    return dotProduct<T>(v_device_p, v_device_p, n, partial_sums_device_p);
}

template<typename T>
T vectorNormSquared(const std::vector<T>& v) {
    constexpr auto elements_per_block = dotProductConfig::elements_per_block; 
    const auto blocks = ((v.size() / elements_per_block) + 1);
    const auto padded_size = blocks * elements_per_block;

    const auto data_device = CudaArray<T, 2> {
        {v, padded_size, true},
        {{}, blocks}
    };

    const auto [v_device_p, partial_sums_device_p] = data_device.get_pointers_to_data();

    return dotProduct<T>(v_device_p, v_device_p, v.size(), partial_sums_device_p);
}

// scalar * v + w
template<typename T>
__global__ void multiplyVectorWithScalarAndAddItToAnotherVector(const T scalar, T const * const v_device_p, T const * const w_device_p, T * const result_device_p) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    result_device_p[i] = scalar * v_device_p[i] + w_device_p[i];
}

// scalar * v + w
template<typename T>
std::vector<T> multiplyVectorWithScalarAndAddItToAnotherVector(const T scalar, const std::vector<T>& v, const std::vector<T>& w) {
    assert(v.size() == w.size());
    const auto padded_size = (((v.size() / warp_size) + 1) * warp_size);

    const auto data_device = CudaArray<T, 3> {
        {v, padded_size},
        {w, padded_size},
        {{}, padded_size}
    };

    const auto [v_device_p, w_device_p, result_device_p] = data_device.get_pointers_to_data();


    multiplyVectorWithScalarAndAddItToAnotherVector<<<padded_size / warp_size, warp_size>>>(scalar, v_device_p, w_device_p, result_device_p);
    return data_device.copy_data_to_host_vector(result_device_p, v.size());
}

// Modified Asolve from page 86 (110 pdf) of
// https://www.grad.hr/nastava/gs/prg/NumericalRecipesinC.pdf
// We assume that element (i, i) of sparse_A is stored at sparse_A[i * m]
template<typename T>
__global__ void Asolve_diagonal(size_t n, size_t m, T const * const sparse_A, T const * const b, T * const x) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;

    // This might be abel to be optimized away by setting padded part of sparse A to non_zero
    if (i < n) {
        x[i] = b[i] / sparse_A[i * m];
    } else {
        x[i] = 0;
    }
}


template<typename T>
__global__ void copyData(T const * const source_device_p, T * const destination_device_p) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;
    destination_device_p[i] = source_device_p[i];
}

 

template <typename T>
ReturnOfSparseBiCGCUDA<T> sparseBiCGCUDA(
        const size_t n, const size_t m,    
        const std::vector<T>& sparse_A,
        const std::vector<T>& sparse_A_transposed,
        const std::vector<size_t>& indecies,
        const std::vector<T>& b,
        const ConfigurationForIonosphereGPUSolver<T>& config) {

    assert(sparse_A.size() == n * m);
    assert(sparse_A_transposed.size() == n * m);
    assert(indecies.size() == n * m);
    assert(b.size() == n);

    const auto height = ((n / warp_size) + 1) * warp_size;

    const auto indecies_on_device = CudaArray<size_t, 1> {
        {indecies, height * m, true}
    };
    const auto [indecies_device_p] = indecies_on_device.get_pointers_to_data();

    const auto blocks_for_dot_product = n / dotProductConfig::elements_per_block + 1;
    // Dot product needs more padding
    [[maybe_unused]] const auto padded_size_for_dot_product = blocks_for_dot_product * dotProductConfig::elements_per_block;
    
    auto data_on_device = CudaArray<T, 13> {
        {sparse_A, height * m, true},
        {sparse_A_transposed, height * m, true},
        {b, padded_size_for_dot_product, true},
        {/* x */{}, height, true}, // Initial guess of x is zero
        {/* r */{}, padded_size_for_dot_product, true},
        {/* rr */{}, padded_size_for_dot_product, true},
        {/* z */{}, padded_size_for_dot_product, true},
        {/* zz */{}, height},
        {/* p */{}, height},
        {/* pp */{}, padded_size_for_dot_product, true},
        {/* best_solution */{}, height},
        {/* partial sums for dot product */{}, blocks_for_dot_product},
        {/* temp */{}, height}
    };
    
    const auto [
        sparse_A_device_p,
        sparse_A_transposed_device_p,
        b_device_p,
        x_device_p,
        r_device_p,
        rr_device_p,
        z_device_p,
        zz_device_p,
        p_device_p,
        pp_device_p,
        best_solution_device_p,
        partial_sums_for_dot_product_device_p,
        temp_device_p
    ] = data_on_device.get_pointers_to_data();
    
    
    auto number_of_restarts = int { 0 };
    auto min_error = std::numeric_limits<T>::max();
    auto iteration = int { 0 };
    // This is part of the restart mechanism
   const auto bnrm = std::sqrt(vectorNormSquared<T>(b_device_p, n, partial_sums_for_dot_product_device_p));
    do {
        copyData<<<height / warp_size, warp_size>>>(best_solution_device_p, x_device_p);

       // *******************************************************************
        // From this point onwards comments will refere numecial recipies in C
        // *******************************************************************
        // Atime(n, x, r, 0);
        sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
            m, sparse_A_device_p, indecies_device_p, x_device_p, r_device_p
        );

        // r[j]=b[j]-r[j];
        vectorSubtraction<T><<<height / warp_size, warp_size>>>(b_device_p, r_device_p, r_device_p);

        if (config.use_minimum_residual_variant) {
            // atimes(n,r,rr,0);
            sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
                m, sparse_A_device_p, indecies_device_p, r_device_p, rr_device_p
            );
        } else {
            // rr[j]=r[j];
            copyData<<<height / warp_size, warp_size>>>(r_device_p, rr_device_p);
        }

        // Asumed we use itol == 1 convergence test
        // if (itol == 1) {
        //    bnrm=snrm(n,b,itol);

        auto Asolve = [&, sparse_A_device_p = sparse_A_device_p, m]
        (T const * const local_b_device_p, T * const local_x_device_p, [[maybe_unused]] const bool transposed = false) -> void {
            switch (config.precondition) {
                case Precondition::none: {
                    copyData<<<height / warp_size, warp_size>>>(local_b_device_p, local_x_device_p);
                } break;
                case Precondition::diagonal: {
                    Asolve_diagonal<T><<<height / warp_size, warp_size>>>(n, m, sparse_A_device_p, local_b_device_p, local_x_device_p);
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
        auto failcount = int{ 0 };
        auto first_iteration = true;
        for (;iteration < config.max_iterations;) {
            // somehow iteration can not be in above for statement when linkin this with mpicc
            iteration += 1;

            // asolve(n,rr,zz,1);
            Asolve(rr_device_p, zz_device_p, true);
            // for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
            const auto bknum = dotProduct<T>(z_device_p, rr_device_p, n, partial_sums_for_dot_product_device_p);

            if (first_iteration) {
                first_iteration = false;
                // p[j]=z[j];
                copyData<<<height / warp_size, warp_size>>>(z_device_p, p_device_p);
                // pp[j]=zz[j];
                copyData<<<height / warp_size, warp_size>>>(zz_device_p, pp_device_p);
            } else {
                const auto bk = bknum / bkden;
                // p[j]=bk*p[j]+z[j];
                // pp[j]=bk*pp[j]+zz[j];
                multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(bk, p_device_p, z_device_p, p_device_p);
                multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(bk, pp_device_p, zz_device_p, pp_device_p);


            }
            bkden = bknum;
            // atimes(n,p,z,0);
            sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
                m, sparse_A_device_p, indecies_device_p, p_device_p, z_device_p
            );

            // for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
            const auto akden = dotProduct<T>(z_device_p, pp_device_p, n, partial_sums_for_dot_product_device_p);
            const auto ak = bknum / akden;

            // for (j=1;j<=n;j++) {
            //     x[j] += ak*p[j];
            //     r[j] -= ak*z[j];
            //     rr[j] -= ak*zz[j];
            // }   
            multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(ak, p_device_p, x_device_p, x_device_p);
            // *******************************************************************************
            // From this point onwards comments will not refere numecial recipies in C anymore
            // *******************************************************************************
            // Copying from orginal implementation instead of this:
            //multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(-ak, z_device_p, r_device_p, r_device_p);
            //multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(-ak, zz_device_p, rr_device_p, rr_device_p);

            // We do this:
            
            sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
                m, sparse_A_device_p, indecies_device_p, x_device_p, temp_device_p
            );
            vectorSubtraction<T><<<height / warp_size, warp_size>>>(b_device_p, temp_device_p, r_device_p);

            sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
                m, sparse_A_transposed_device_p, indecies_device_p, x_device_p, temp_device_p
            );

            vectorSubtraction<T><<<height / warp_size, warp_size>>>(b_device_p, temp_device_p, rr_device_p); 
            // We only support pole gauge at the moment:
            // if (config.gauge == Gauge::pole) {
            //     data_on_device.zero_out_memory(x_device_p, 1);
            // }
            

            // asolve(n,r,z,0);
            Asolve(r_device_p, z_device_p);

            
            // I think olderr in vlasiator solver is unused
            const auto error = std::sqrt(vectorNormSquared<T>(r_device_p, n, partial_sums_for_dot_product_device_p)) / bnrm;
            
            if (error < min_error) {
                copyData<<<height / warp_size, warp_size>>>(x_device_p, best_solution_device_p);
                min_error = error;
                failcount = 0;
            } else {
                copyData<<<height / warp_size, warp_size>>>(best_solution_device_p, x_device_p);
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
    
    return ReturnOfSparseBiCGCUDA<T>{
        iteration,
        number_of_restarts,
        min_error,
        data_on_device.copy_data_to_host_vector(best_solution_device_p, n)
    };
    
};

template <typename T>
ReturnOfSparseBiCGCUDA<T> sparseBiCGSTABCUDA(
        const size_t n, const size_t m,    
        const std::vector<T>& sparse_A,
        const std::vector<size_t>& indecies,
        const std::vector<T>& b,
        const ConfigurationForIonosphereGPUSolver<T>& config) {
    assert(sparse_A.size() == n * m);
    assert(indecies.size() == n * m);
    assert(b.size() == n);

    timer::time("sparseBiCGSTABCUDA::init");
    const auto height = ((n / warp_size) + 1) * warp_size;
    assert(height % warp_size == 0);

    const auto indecies_on_device = CudaArray<size_t, 1> {
        {indecies, height * m, true}
    };
    const auto [indecies_device_p] = indecies_on_device.get_pointers_to_data();

    const auto blocks_for_dot_product = n / dotProductConfig::elements_per_block + 1;
    // Dot product needs more padding
    [[maybe_unused]] const auto padded_size_for_dot_product = blocks_for_dot_product * dotProductConfig::elements_per_block;
    
    //const auto size_of_padding_for_dot_product = padded_size_for_dot_product - n;

    auto data_on_device = CudaArray<T, 15> {
        {sparse_A, height * m, true},
        {b, padded_size_for_dot_product, true},
        {/* x */{}, padded_size_for_dot_product},
        {/* r */{}, padded_size_for_dot_product, true},
        {/* r_hat */{}, padded_size_for_dot_product, true},
        {/* p */{}, height, true},
        {/* v */{}, padded_size_for_dot_product, true},
        {/* h */{}, height},
        {/* s */{}, padded_size_for_dot_product, true},
        {/* t */{}, padded_size_for_dot_product, true},
        {/* y */{}, padded_size_for_dot_product, true},
        {/* z */{}, padded_size_for_dot_product, true},
        {/* temp */{}, padded_size_for_dot_product},
        {/* best_solution */{}, height, true},
        {/* partial sums for dot product */{}, blocks_for_dot_product}
    };

    const auto [
        sparse_A_device_p,
        b_device_p,
        x_device_p,
        r_device_p,
        r_hat_device_p,
        p_device_p,
        v_device_p,
        h_device_p,
        s_device_p,
        t_device_p,
        y_device_p,
        z_device_p,
        temp_device_p,
        best_solution_device_p,
        partial_sums_for_dot_product_device_p
    ] = data_on_device.get_pointers_to_data();

    
    timer::time("sparseBiCGSTABCUDA::init");
    //zeroOutData<<<(size_of_padding_for_dot_product / warp_size) + 1, warp_size>>>(b_device_p + n, size_of_padding_for_dot_product);
    const auto b_norm = std::sqrt(vectorNormSquared<T>(b_device_p, n, partial_sums_for_dot_product_device_p));
    
    
    auto number_of_restarts = int { 0 };
    auto min_error = std::numeric_limits<T>::max();
    auto iteration = int { 0 };
    // This is part of the restart mechanism
    do {
        // *******************************************************************************
        // Numbers with star (eg. 3* ) will refere preconditionless BiCGSTAB steps Wikipedia
        // *******************************************************************************
        // 1*
        timer::time("sparseBiCGSTABCUDA::1*");
        copyData<T><<<height / warp_size, warp_size>>>(best_solution_device_p, x_device_p);
        sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
            m, sparse_A_device_p, indecies_device_p, x_device_p, temp_device_p
        );
        vectorSubtraction<T><<<height / warp_size, warp_size>>>(b_device_p, temp_device_p, r_device_p);
    
        timer::time("sparseBiCGSTABCUDA::1*", "sparseBiCGSTABCUDA::2*");


    
        // 2*
        copyData<<<height / warp_size, warp_size>>>(r_device_p, r_hat_device_p);

        timer::time("sparseBiCGSTABCUDA::2*", "sparseBiCGSTABCUDA::3*");
        // 3*
        T rho_im1 { 1 };
        T alpha { 1 };
        T omega { 1 };

        timer::time("sparseBiCGSTABCUDA::3*", "sparseBiCGSTABCUDA::4*");
        // 4*
        zeroOutData<<<padded_size_for_dot_product / warp_size, warp_size>>>(v_device_p, padded_size_for_dot_product); 
        zeroOutData<<<height / warp_size, warp_size>>>(p_device_p, height); 


        timer::time("sparseBiCGSTABCUDA::4*", "sparseBiCGSTABCUDA::5*");
        // 5* (After this numbers with two ** (ed. 3**) will refere the steps inside the loop)
        auto failcount = int{ 0 };
        for (;iteration < config.max_iterations;) {
            // somehow iteration can not be in above for statement when linkin this with mpicc
            iteration += 1;
            timer::time("sparseBiCGSTABCUDA::1**");
            // 1**

            T rho_i = dotProduct(r_hat_device_p, r_device_p, n, partial_sums_for_dot_product_device_p);
            timer::time("sparseBiCGSTABCUDA::1**", "sparseBiCGSTABCUDA::2**");
            // 2**
            T beta = (rho_i / alpha) * (rho_im1 / omega);

            // Above is last use of rho_im1 so we can store current rho for next iteration
            rho_im1 = rho_i;

            timer::time("sparseBiCGSTABCUDA::2**", "sparseBiCGSTABCUDA::3**");
            // 3**
            multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(-omega, v_device_p, p_device_p, temp_device_p);
            multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(beta, temp_device_p, r_device_p, p_device_p);
         

            timer::time("sparseBiCGSTABCUDA::3**", "sparseBiCGSTABCUDA::4**");
            // 4**
            switch (config.precondition) {
                
                case Precondition::diagonal: {
                    Asolve_diagonal<T><<<height / warp_size, warp_size>>>(n, m, sparse_A_device_p, p_device_p, y_device_p);
                } break;
                case Precondition::none: {
                    copyData<<<height / warp_size, warp_size>>>(p_device_p, y_device_p);
                } break;
            }

            sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
                m, sparse_A_device_p, indecies_device_p, y_device_p, v_device_p);
            
           timer::time("sparseBiCGSTABCUDA::4**", "sparseBiCGSTABCUDA::5**");
            // 5**
            alpha = rho_i / dotProduct(r_hat_device_p, y_device_p, n, partial_sums_for_dot_product_device_p);

            timer::time("sparseBiCGSTABCUDA::5**", "sparseBiCGSTABCUDA::6**");
            // 6**
            multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(alpha, y_device_p, x_device_p, h_device_p);
           timer::time("sparseBiCGSTABCUDA::6**", "sparseBiCGSTABCUDA::7**");
            // 7**
            sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
                m, sparse_A_device_p, indecies_device_p, h_device_p, temp_device_p  
            );
            // residue in temp
            vectorSubtraction<T><<<height / warp_size, warp_size>>>(b_device_p, temp_device_p, temp_device_p);
            const auto error_of_h = std::sqrt(vectorNormSquared<T>(temp_device_p, n, partial_sums_for_dot_product_device_p)) / b_norm;

            if (error_of_h < min_error) {
                copyData<<<height / warp_size, warp_size>>>(h_device_p, best_solution_device_p);
                min_error = error_of_h;
                failcount = 0;
            }
            if (min_error < config.relative_L2_convergence_threshold) {
                break;
            }

            timer::time("sparseBiCGSTABCUDA::7**", "sparseBiCGSTABCUDA::8**");
            // 8** 

            multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(-alpha, v_device_p, r_device_p, s_device_p);

            timer::time("sparseBiCGSTABCUDA::8**", "sparseBiCGSTABCUDA::9**");
            // 9**

            switch (config.precondition) {
                
                case Precondition::diagonal: {
                    // K_2^(-1)s in z
                    Asolve_diagonal<T><<<height / warp_size, warp_size>>>(n, m, sparse_A_device_p, s_device_p, z_device_p);
                } break;
                case Precondition::none: {
                    // s in z
                    copyData<<<height / warp_size, warp_size>>>(s_device_p, z_device_p);
                } break;
            }


            sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
                m, sparse_A_device_p, indecies_device_p, z_device_p, t_device_p
            );

        
            timer::time("sparseBiCGSTABCUDA::9**", "sparseBiCGSTABCUDA::10**");
            // 10**

            omega = dotProduct(t_device_p, s_device_p, n, partial_sums_for_dot_product_device_p) / vectorNormSquared(t_device_p, n, partial_sums_for_dot_product_device_p);

            timer::time("sparseBiCGSTABCUDA::10**", "sparseBiCGSTABCUDA::11**");
            // 11**

            multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(omega, z_device_p, h_device_p, x_device_p);

            timer::time("sparseBiCGSTABCUDA::11**", "sparseBiCGSTABCUDA::12**");
            // 12** 
    

            sparseMatrixVectorProduct<T><<<height / warp_size, warp_size>>>(
                m, sparse_A_device_p, indecies_device_p, x_device_p, temp_device_p
            );
 
            vectorSubtraction<T><<<height / warp_size, warp_size>>>(b_device_p, temp_device_p, temp_device_p);

            const auto error = std::sqrt(vectorNormSquared<T>(temp_device_p, n, partial_sums_for_dot_product_device_p)) / b_norm;

            if (error < min_error) {
                copyData<<<height / warp_size, warp_size>>>(x_device_p, best_solution_device_p);
                min_error = error;
                failcount = 0;
            } else {
                copyData<<<height / warp_size, warp_size>>>(best_solution_device_p, x_device_p);
                ++failcount;
            }

            if (min_error < config.relative_L2_convergence_threshold) {
                break;
            }
            if (failcount > config.max_failure_count || error > config.max_error_growth_factor * min_error) {
                ++number_of_restarts;
                break;
            }
            timer::time("sparseBiCGSTABCUDA::12**", "sparseBiCGSTABCUDA::13**");
            // 13**
            multiplyVectorWithScalarAndAddItToAnotherVector<T><<<height / warp_size, warp_size>>>(-omega, t_device_p, s_device_p, r_device_p);
            timer::time("sparseBiCGSTABCUDA::13**");

        }
        timer::time("sparseBiCGSTABCUDA::5*");
    // Restart mechanism
    } while (min_error > config.relative_L2_convergence_threshold && iteration < config.max_iterations);

    return ReturnOfSparseBiCGCUDA<T>{
        iteration,
        number_of_restarts,
        min_error,
        data_on_device.copy_data_to_host_vector(best_solution_device_p, n)
    };
    
};
 
}
