#include "../include/ionosphere_gpu_solver.hpp"
#include "../include/CudaArray.hpp"
#include <iostream>

template
std::vector<double> ionogpu::matrixVectorProduct<double>(const std::vector<double>& M, const std::vector<double>& v);

template
std::vector<double> ionogpu::sparseMatrixVectorProduct<double>(
    const size_t n, const size_t m,    
    const std::vector<double>& A,
    const std::vector<size_t>& indecies,
    const std::vector<double>& b);


std::vector<double> ionogpu::sparsebcgCUDA(
        const size_t n, const size_t m,    
        const std::vector<double>& A,
        const std::vector<size_t>& indecies,
        const std::vector<double>& b) {
    
            
    
    return {};
};
 