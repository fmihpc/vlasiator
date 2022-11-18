#pragma once
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include "ionosphere_gpu_solver.hpp"
namespace ionogpu {
    
template <typename T>
class CudaArray {
    size_t size_;
    T* data_;

    size_t size_in_bytes_() const {
        return sizeof(T) * size_;
    }
    public:
    CudaArray(const size_t size) : size_{size} {
        cudaMalloc((void**)&data_, size_in_bytes_());
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

    auto M_v_and_Mv_on_device = CudaArray<T>(height * width + width + height);
    
    M_v_and_Mv_on_device.copy_data_to_device(M.data(), 0, M.size());
    M_v_and_Mv_on_device.copy_data_to_device(v.data(), width * height, v.size());
    // We zero out Mv but also the badded part of v
    M_v_and_Mv_on_device.zeroOutMemory(width * height + width, height);
   
    const auto blocks = height / 32; 
    const auto threads_per_block = 32;
    
    auto M_device_p = M_v_and_Mv_on_device.data();
    auto v_device_p = M_device_p + width * height;
    auto Mv_device_p = v_device_p + width;
    matrixVectorProduct<double><<<blocks, threads_per_block>>>(M_device_p, v_device_p, Mv_device_p, width); 
    return M_v_and_Mv_on_device.copy_data_to_host_vector(width * height + width, v.size());  
    
}


template <typename T>
__global__ void sparseMatrixVectorProduct(
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

    const auto pre_b = [&] () {
        auto temp  = std::vector<T>(n * m);
        std::transform(indecies.begin(), indecies.end(), temp.begin(), [&](const auto i) noexcept {return b[i];});
        return temp;
    }();


    auto sparse_A_pre_b_and_x_on_device = CudaArray<double>(height * m + height * m + height);
    const T* sparse_A_device_p = sparse_A_pre_b_and_x_on_device.copy_vector_to_device_p(sparse_M, 0);
    sparse_A_pre_b_and_x_on_device.zeroOutMemory(sparse_M.size(), n * height - sparse_M.size());
    const T* pre_b_device_p = sparse_A_pre_b_and_x_on_device.copy_vector_to_device_p(pre_b, height * m);
    sparse_A_pre_b_and_x_on_device.zeroOutMemory(height * m + height * m, n);
    auto x_device_p = const_cast<T*>(pre_b_device_p + height * m);

    const auto blocks = height / 32;
    const auto threads_per_block = 32;

    sparseMatrixVectorProduct<double><<<blocks, threads_per_block>>>(m, sparse_A_device_p, pre_b_device_p, x_device_p);
    return sparse_A_pre_b_and_x_on_device.copy_data_to_host_vector(height * m + height * m, n);
};
/* 
std::vector<double> Atimes(
    const size_t n, const size_t m,    
    const std::vector<double>& A,
    const std::vector<size_t>& indecies,
    const std::vector<double>& b
);
 */
}
