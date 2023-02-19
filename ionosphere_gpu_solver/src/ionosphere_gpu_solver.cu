#include "../include/ionosphere_gpu_solver.hpp"
#include "../include/IonogpuMemoryAllocator.hpp"

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
ionogpu::ReturnOfSparseBiCGCUDA<double> ionogpu::sparseBiCGCUDA<double>(
        const size_t n, const size_t m,    
        const std::vector<double>& sparse_A,
        const std::vector<double>& sparse_A_trasposed,
        const std::vector<size_t>& indecies,
        const std::vector<double>& b,
        const ionogpu::ConfigurationForIonosphereGPUSolver<double>& config);

template 
ionogpu::ReturnOfSparseBiCGCUDA<double> ionogpu::sparseBiCGSTABCUDA<double>(
        const size_t n, const size_t m,    
        const std::vector<double>& sparse_A,
        const std::vector<size_t>& indecies,
        const std::vector<double>& b,
        const ionogpu::ConfigurationForIonosphereGPUSolver<double>& config);

        
 

template
std::vector<double> ionogpu::vectorAddition<double>(const std::vector<double>& a, const std::vector<double>& b);

template
std::vector<double> ionogpu::vectorSubtraction<double>(const std::vector<double>& a, const std::vector<double>& b);

template
double ionogpu::vectorNormSquared<double>(const std::vector<double>& v);

template
std::vector<double> ionogpu::vectorElementwiseMultiplication<double>(const std::vector<double>& a, const std::vector<double>& b);

template
std::vector<double> ionogpu::vectorElementwiseDivision<double>(const std::vector<double>& a, const std::vector<double>& b);


template
double ionogpu::dotProduct<double>(const std::vector<double>& v, const std::vector<double>& w);

// scalar * v + w
template
std::vector<double> ionogpu::multiplyVectorWithScalarAndAddItToAnotherVector<double>(const double scalar, const std::vector<double>& v, const std::vector<double>& w);


template
std::vector<float> ionogpu::matrixVectorProduct<float>(const std::vector<float>& M, const std::vector<float>& v);

template
std::vector<float> ionogpu::preSparseMatrixVectorProduct<float>(
    const size_t n, const size_t m,    
    const std::vector<float>& A,
    const std::vector<size_t>& indecies,
    const std::vector<float>& b);

template
std::vector<float> ionogpu::sparseMatrixVectorProduct<float>(
    const size_t n, const size_t m,    
    const std::vector<float>& A,
    const std::vector<size_t>& indecies,
    const std::vector<float>& b);

template 
ionogpu::ReturnOfSparseBiCGCUDA<float> ionogpu::sparseBiCGCUDA<float>(
        const size_t n, const size_t m,    
        const std::vector<float>& sparse_A,
        const std::vector<float>& sparse_A_trasposed,
        const std::vector<size_t>& indecies,
        const std::vector<float>& b,
        const ionogpu::ConfigurationForIonosphereGPUSolver<float>& config);

template 
ionogpu::ReturnOfSparseBiCGCUDA<float> ionogpu::sparseBiCGSTABCUDA<float>(
        const size_t n, const size_t m,    
        const std::vector<float>& sparse_A,
        const std::vector<size_t>& indecies,
        const std::vector<float>& b,
        const ionogpu::ConfigurationForIonosphereGPUSolver<float>& config);

        
 

template
std::vector<float> ionogpu::vectorAddition<float>(const std::vector<float>& a, const std::vector<float>& b);

template
std::vector<float> ionogpu::vectorSubtraction<float>(const std::vector<float>& a, const std::vector<float>& b);

template
float ionogpu::vectorNormSquared<float>(const std::vector<float>& v);

template
std::vector<float> ionogpu::vectorElementwiseMultiplication<float>(const std::vector<float>& a, const std::vector<float>& b);

template
std::vector<float> ionogpu::vectorElementwiseDivision<float>(const std::vector<float>& a, const std::vector<float>& b);


template
float ionogpu::dotProduct<float>(const std::vector<float>& v, const std::vector<float>& w);

// scalar * v + w
template
std::vector<float> ionogpu::multiplyVectorWithScalarAndAddItToAnotherVector<float>(const float scalar, const std::vector<float>& v, const std::vector<float>& w);