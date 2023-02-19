#include "ionosphere_gpu_solver.hpp"
#include "tools.hpp"
#include <cassert>
#include <iterator>
#include <vector>

using test_type = TEST_TYPE_PROTOTYPE;

constexpr size_t n = 1000;
static const auto M = [] {
    auto temp = std::vector<test_type>(n * n, 0);
    for (size_t i = 0; i < n; ++i) {
        temp[i * n + i] = 102;
        temp[(i * n + i + 2) % (n * n)] = 1;
        temp[(i * n + i + 3) % (n * n)] = 0.4123;
        temp[(i * n + i - 4) % (n * n)] = 0.4123;
        temp[(i * n + i - 1) % (n * n)] = 0.4123;
    }
    return temp;
}(); 
static const auto Mx = [] {
    auto temp = std::vector<test_type>(n);
    for (size_t i = 0; i < temp.size(); i += 1) {
        temp[i] = 1;
    }
    return temp;
}();

auto main() -> int {

    const auto max_number_of_nonzero_elements_on_each_row = 32;
    const auto [indecies, sparse_M] = ionogpu::testing::create_sparse_matrix_from_dense_matrix(M, n, max_number_of_nonzero_elements_on_each_row);
    
    const auto M_T = ionogpu::testing::transpose_matrix(M, n);
    [[maybe_unused]] const auto [sink, sparse_M_T] = ionogpu::testing::create_sparse_matrix_from_dense_matrix(M_T, n, max_number_of_nonzero_elements_on_each_row);
        

    
    for (size_t i = 0; i < n; ++i) {
        assert(M[i * n + i] != 0);
        assert(M_T[i * n + i] != 0);
        assert(sparse_M[max_number_of_nonzero_elements_on_each_row * i] != 0);
        assert(sparse_M_T[max_number_of_nonzero_elements_on_each_row * i] != 0);
    }


{

    const auto config = ionogpu::ConfigurationForIonosphereGPUSolver<test_type> {
            .max_iterations = 2000,
            .max_failure_count = 10,
            .max_error_growth_factor = 100,
            .relative_L2_convergence_threshold = 1e-6,
            .precondition = ionogpu::Precondition::diagonal,
            .use_minimum_residual_variant = true,
            .gauge = ionogpu::Gauge::none
    };

    const auto [number_of_iterations, number_of_restarts, min_error, x] = ionogpu::sparseBiCGCUDA(
        n, max_number_of_nonzero_elements_on_each_row,
        sparse_M,
        sparse_M_T,
        indecies,
        Mx,
        config
    );
    assert(number_of_iterations <= config.max_iterations);
    assert(min_error < config.relative_L2_convergence_threshold);
}
{ 
    const auto config = ionogpu::ConfigurationForIonosphereGPUSolver<test_type> {
            .max_iterations = 1000,
            .max_failure_count = 10,
            .max_error_growth_factor = 100,
            .relative_L2_convergence_threshold = 1e-6,
            .precondition = ionogpu::Precondition::diagonal,
            .use_minimum_residual_variant = false,
            .gauge = ionogpu::Gauge::none
    };

    const auto [number_of_iterations, number_of_restarts, min_error, x] = ionogpu::sparseBiCGCUDA(
        n, max_number_of_nonzero_elements_on_each_row,
        sparse_M,
        sparse_M_T,
        indecies,
        Mx,
        config
    );
    
    assert(number_of_iterations <= config.max_iterations);
    assert(min_error < config.relative_L2_convergence_threshold);
} 

    return 0;

}