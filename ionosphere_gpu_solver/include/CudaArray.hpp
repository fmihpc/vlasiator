#pragma once
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <vector>
#include <cassert>
#include <iostream>
namespace ionogpu__ {
    
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

    void copy_data_to_device(const T* data, const size_t offset_in_elements, const size_t size_in_elements) {
        cudaMemcpy(data_ + offset_in_elements, data, size_in_elements * sizeof(T), cudaMemcpyHostToDevice); 
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
};



template <typename T>
__global__ void MatrixVectorProduct(T const * const  M, T const * const v, T* const Mv, const size_t n) {
    const auto i = blockDim.x * blockIdx.x + threadIdx.x;

    for (size_t j = 0; j < n; ++j) {
        Mv[i] += M[i * n + j] * v[j];
    }
}

template <typename T>
std::vector<T> MatrixVectorProduct(const std::vector<T>& M, const std::vector<T>& v) { 
    assert(M.size() == v.size() * v.size());
    const auto heigth = ((v.size() / 32) + 1) * 32;
    const auto width = v.size();

    auto M_v_and_Mv_on_device = CudaArray<T>(heigth * width + width + heigth);
    
    M_v_and_Mv_on_device.copy_data_to_device(M.data(), 0, M.size());
    M_v_and_Mv_on_device.copy_data_to_device(v.data(), width * heigth, v.size());
    // We zero out Mv but also the badded part of v
    M_v_and_Mv_on_device.zeroOutMemory(width * heigth + width, heigth);
   
    const auto blocks = heigth / 32; 
    const auto threads_per_block = 32;
    
    auto M_device_p = M_v_and_Mv_on_device.data();
    auto v_device_p = M_device_p + width * heigth;
    auto Mv_device_p = v_device_p + width;
    MatrixVectorProduct<double><<<blocks, threads_per_block>>>(M_device_p, v_device_p, Mv_device_p, width); 
    return M_v_and_Mv_on_device.copy_data_to_host_vector(width * heigth + width, v.size());  
    
}

}
