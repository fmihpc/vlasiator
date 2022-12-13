#include "ionosphere_gpu_solver.hpp"
#include "tools.hpp"
#include <cassert>
#include <iterator>
#include <vector>

using test_type = double;

constexpr size_t n = 10000;
static const auto M = [] {
    auto temp = std::vector<test_type>(n * n, 0);
    for (size_t i = 0; i < n; ++i) {
        temp[i * n + i] = 30;
        temp[(i * n - 7 + i) % n * n] = 2 * 20;
        temp[(i * n - 9 + i) % n * n] = 3 * 20;
        temp[(i * n + 7 + i) % n * n] = 2 * 20;
        temp[(i * n + 9 + i) % n * n] = 3 * 20;
    }
    return temp;
}(); 
static const auto x_correct = [] {
    auto temp = std::vector<test_type>(n);
    for (size_t i = 0; i < temp.size(); i += 1) {
        temp[i] = 1 + (i % 2);
    }
    return temp;
}();

auto main() -> int {

    const auto max_number_of_nonzero_elements_on_each_row = 32;
    const auto [indecies, sparse_M] = ionogpu::testing::create_sparse_matrix_from_dense_matrix(M, n, max_number_of_nonzero_elements_on_each_row);
        
    for (size_t i = 0; i < n; ++i) {
        assert(M[i * n + i] != 0);
        assert(sparse_M[max_number_of_nonzero_elements_on_each_row * i] != 0);
    }


    const auto Mx = [&] {
        auto temp = std::vector<test_type>(n, 0);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                temp[i] += M[i * n + j] * x_correct[j];
            }
        }
        return temp;
    }();

    const auto config = ionogpu::ConfigurationForIonosphereGPUSolver<test_type> {
            .max_iterations = 1000,
            .max_failure_count = 20,
            .max_error_growth_factor = 100,
            .relative_L2_convergence_threshold = 0.000001,
            .precondition = ionogpu::Precondition::diagonal,
            .use_minimum_residual_variant = false,
            .gauge = ionogpu::Gauge::none
    };

    const auto [number_of_iterations, number_of_restarts, min_error, x] = ionogpu::sparseBiCGSTABCUDA(
        n, max_number_of_nonzero_elements_on_each_row,
        sparse_M,
        indecies,
        Mx,
        config
    );
    std::cout << "Number of iterations: " << number_of_iterations << "\nNumber of restarts: " << number_of_restarts << "\n"; 
    std::cout << "Min error: " << min_error << "\n";
/* 
    std::cout  << "\nSolved x:\n";
    for (const auto e : x) {
        std::cout << e << " ";
    }
    std::cout << "\n";
    std::cout << "\nMx_gpu:\n";
    const auto Mx_gpu = [&] {
        auto temp = std::vector<test_type>(n, 0);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                temp[i] += M[i * n + j] * x[j];
            }
        }
        return temp;
    }();
    for (const auto e : Mx_gpu) {
        std::cout << e << " ";
    }
    std::cout << "\n\n Mx_correct:\n";

    for (const auto e : Mx) {
        std::cout << e << " ";
    }
    std::cout << "\n";  */

    assert(number_of_iterations <= config.max_iterations);
    
    [[maybe_unused]] const auto [absolute_error, relative_error] = ionogpu::testing::calculate_absolute_and_relative_error_of_range(x, x_correct);

    assert(absolute_error < 0.0001);
    assert(relative_error < 0.0001);

    return 0;

}