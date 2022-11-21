#include "../include/ionosphere_gpu_solver.hpp"
#include "../include/CudaArray.hpp"
#include <iostream>

template
std::vector<double> ionogpu::matrixVectorProduct<double>(const std::vector<double>& M, const std::vector<double>& v);

template
std::vector<double> ionogpu::preSparseMatrixVectorProduct<double>(
    const size_t n, const size_t m,    
    const std::vector<double>& A,
    const std::vector<size_t>& indecies,
    const std::vector<double>& b);

template
std::vector<double> ionogpu::sparseMatrixVectorProduct<double>(
    const size_t n, const size_t m,    
    const std::vector<double>& A,
    const std::vector<size_t>& indecies,
    const std::vector<double>& b);

template 
std::vector<double> ionogpu::sparseBiCGCUDA<double>(
        const size_t n, const size_t m,    
        const std::vector<double>& sparse_A,
        const std::vector<double>& sparse_A_trasposed,
        const std::vector<size_t>& indecies,
        const std::vector<double>& b,
        const ConfigurationForSparsebcgCUDA<double>& config);
 

template
std::vector<double> ionogpu::vectorAddition<double>(const std::vector<double>& a, const std::vector<double>& b);

template
std::vector<double> ionogpu::vectorSubtraction<double>(const std::vector<double>& a, const std::vector<double>& b);

template
double ionogpu::vectorNorm<double>(const std::vector<double>& v);