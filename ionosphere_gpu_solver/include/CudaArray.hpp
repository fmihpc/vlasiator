#pragma once
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <vector>
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <initializer_list>
#include <numeric>
#include "ionosphere_gpu_solver.hpp"
namespace ionogpu {
    
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
                zeroOutMemory(current_place + padded_vec.vec.size(), padded_vec.size() - padded_vec.vec.size());
            }
            current_place += std::max(padded_vec.vec.size(), padded_vec.minimun_size);
        }

        assert(N == padded_vecs.size());
        std::copy(pointers_to_data.begin(), pointers_to_data.end(), pointers_to_data_.begin());
        
    }

    ~CudaArray() {
        cudaFree(data_);
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
    }
    // Move assignment
    CudaArray& operator=(CudaArray&& other) noexcept {
        size_ = other.size_; 
        data_ = other.data_;
        other.data_ = nullptr;
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

    void copy_vector_to_device(const std::vector<T>& v, const size_t offset_in_elements) {
        copy_data_to_device(v.data(), offset_in_elements, v.size());
    }

    T* copy_vector_to_device_p(const std::vector<T>& v, const size_t offset_in_elements) {
        return copy_data_to_device_p(v.data(), offset_in_elements, v.size());
    }

    void zeroOutMemory() {
        cudaMemset(data_, 0, size_in_bytes_());
    }

    void zeroOutMemory(const size_t offset_in_elements, const size_t size_in_elements) {
        cudaMemset(data_ + offset_in_elements, 0, size_in_elements * sizeof(T));
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

    for (size_t j = 0; j < n; ++j) {
        Mv[i] += M[i * n + j] * v[j];
    }
}

template <typename T>
std::vector<T> matrixVectorProduct(const std::vector<T>& M, const std::vector<T>& v) { 
    assert(M.size() == v.size() * v.size());
    const auto height = ((v.size() / 32) + 1) * 32;
    const auto width = v.size();

    const auto M_v_and_Mv_on_device = CudaArray<T, 3> {
        {M, width * height},
        {v, width, true},
        {{}, height, true}
    };

    const auto [M_device_p, v_device_p, Mv_device_p] = M_v_and_Mv_on_device.get_pointers_to_data();

    const auto blocks = height / 32; 
    const auto threads_per_block = 32;
 
    matrixVectorProduct<double><<<blocks, threads_per_block>>>(M_device_p, v_device_p, Mv_device_p, width); 
    return M_v_and_Mv_on_device.copy_data_to_host_vector(Mv_device_p, v.size());  
    
}


template <typename T>
__global__ void preSparseMatrixVectorProduct(
    const size_t m,
    T const * const sparse_M,
    T const * const pre_b,
    T* x
) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;

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

    const auto height = ((n / 32) + 1) * 32;

    const auto pre_b = [=, &b] {
        auto temp  = std::vector<T>(n * m);
        std::transform(indecies.begin(), indecies.end(), temp.begin(), [&b](const auto i) noexcept {return b[i];});
        return temp;
    }();

    
    const auto sparse_M_pre_b_and_x_on_device = CudaArray<double, 3> {
        {sparse_M, height * m, true},
        {pre_b, height * m},
        {{}, n, true}
    };

    const auto [sparse_A_device_p, pre_b_device_p, x_device_p] = sparse_M_pre_b_and_x_on_device.get_pointers_to_data();

    const auto blocks = height / 32;
    const auto threads_per_block = 32;

    preSparseMatrixVectorProduct<double><<<blocks, threads_per_block>>>(m, sparse_A_device_p, pre_b_device_p, x_device_p);
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

    const auto height = ((n / 32) + 1) * 32;


    const auto sparse_M_b_and_x_on_device = CudaArray<double, 3> {
        {sparse_M, height * m, true},
        {b, height},
        {{}, height}
    };

    const auto [sparse_A_device_p, b_device_p, x_device_p] = sparse_M_b_and_x_on_device.get_pointers_to_data();

    const auto indecies_on_device = CudaArray<size_t, 1> {
        {indecies, height * m}
    };

    const auto [indecies_device_p] = indecies_on_device.get_pointers_to_data();

    const auto blocks = height / 32;
    const auto threads_per_block = 32;

    sparseMatrixVectorProduct<double><<<blocks, threads_per_block>>>(m, sparse_A_device_p, indecies_device_p, b_device_p, x_device_p);
    return sparse_M_b_and_x_on_device.copy_data_to_host_vector(x_device_p, n);
};


 
}
