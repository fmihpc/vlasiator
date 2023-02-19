#include "ionosphere_gpu_solver.hpp"
#include "tools.hpp"
#include <cassert>
#include <vector>

using test_type = TEST_TYPE_PROTOTYPE;

auto main() -> int {
    
    const auto M = std::vector<test_type>{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    const auto v = std::vector<test_type>{ 1, 2, 3 }; 

    const auto Mv = ionogpu::matrixVectorProduct<test_type>(M, v);
    const auto Mv_correct = std::vector<test_type>{14.0, 32.0, 50.0}; 

    [[maybe_unused]] const auto [absolute_error, relative_error] = ionogpu::testing::calculate_absolute_and_relative_error_of_range(Mv, Mv_correct);
    
    assert(absolute_error < test_type{ 0.00001 });
    assert(relative_error < test_type{ 0.0001 });

    return 0;
}
