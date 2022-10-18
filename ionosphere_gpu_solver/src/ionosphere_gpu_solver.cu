#include "../include/CudaArray.hpp"
#include "../include/ionosphere_gpu_solver.hpp"







template<>
std::vector<double> ionogpu::MatrixVectorProduct<double>(const std::vector<double>& M, const std::vector<double>& v) {
return ionogpu__::MatrixVectorProduct<double>(M, v);
}
