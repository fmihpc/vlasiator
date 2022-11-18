#include "ionosphere_gpu_solver.hpp"
#include "tools.hpp"
#include <cassert>
#include <vector>


auto main() -> int {
    
    const auto M = std::vector<double>{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    const auto v = std::vector<double>{ 1, 2, 3 }; 

    const auto Mv = ionogpu::matrixVectorProduct<double>(M, v);
    const auto Mv_correct = std::vector<double>{14.0, 32.0, 50.0}; 

    const auto [absolute_error, relative_error] = ionogpu::testing::calculate_absolute_and_relative_error_of_range(Mv, Mv_correct);
    
    assert(absolute_error < 0.00001);
    assert(relative_error < 0.0001);

    return 0;
}
