#pragma once
#include "timer.hpp"
#include <concepts>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <utility>
#include <chrono> 
#include <algorithm>
#include <ranges>
#include <cmath>
#include <cassert>
#include <tuple>
namespace ionogpu {
namespace testing {

template <typename... Args, std::invocable<Args...> F>
void benchmark_functions_with_parameters(
        const std::integral auto runtimes,
        const std::initializer_list<std::pair<F, const char *>> functions_with_names,
        const Args... arguments) {
   
    std::cout << "Average, min and max runtimes:\n";
    for (const auto [f, name] : functions_with_names) {
        auto total_time = std::chrono::duration<double>{ 0 }; 
        auto times = std::vector<std::chrono::duration<double>>(runtimes);
        for (auto i = decltype(runtimes){ 0 }; i < runtimes; ++i) {
            const auto begin = std::chrono::steady_clock::now();
            std::invoke(f, arguments...);
            const auto end = std::chrono::steady_clock::now();
            const auto runtime = end - begin;
            total_time += runtime;
            times[i] = runtime;
        }
        const auto average_runtime = total_time / static_cast<double>(runtimes);
        const auto min_runtime = std::ranges::min(times); 
        const auto max_runtime = std::ranges::max(times); 
        std::cout << name << ": " << average_runtime.count() << " "
            << min_runtime.count() << " " << max_runtime.count() <<  "\n";
    }
}
/**
 *  Intermediate values are stored in doubles.
 *  Returns tuple of absolute and relative error
 */
template <std::ranges::range R>
std::tuple<double, double> calculate_absolute_and_relative_error_of_range(const R& v, const R& v_correct) {
    auto total_absolute_error = double{ 0.0 };
    auto total_relative_error = double{ 0.0 };
    
    assert(v.size() == v_correct.size());
    
    
    // Should be as follow in c++23  
    //for (const auto [x, x_correct] : std::ranges::views::zip(v, v_correct)) {
    for (size_t i = 0; i < v.size(); ++i) {
        const auto absolute_error = static_cast<double>(std::fabs(v[i] - v_correct[i]));
        const auto relative_error = [=]() -> double{
            if (v_correct[i] == 0) {
                return static_cast<bool>(v[i]);
            } else {
                return absolute_error / static_cast<double>(v_correct[i]);
            }
        }();

        total_absolute_error += absolute_error;
        total_relative_error += relative_error;
    }

    const auto average_relative_error = total_relative_error / static_cast<double>(v.size());
    
    return std::tuple{total_absolute_error, average_relative_error};
}

template <typename T>
auto create_sparse_matrix_from_dense_matrix(
    const std::vector<T>& M,
    const size_t n,
    const size_t max_num_of_nonzero_elements_on_each_row
) -> std::tuple<std::vector<size_t>, std::vector<T>> {
    auto sparse_M = std::vector<T>(n * max_num_of_nonzero_elements_on_each_row, 0);
    auto indecies = std::vector<size_t>(n * max_num_of_nonzero_elements_on_each_row, 0);
    for (size_t i = 0; i < n; ++i) {
        auto nonzero_elements_on_this_row { 0 };
        for (size_t j = 0; j < n; ++j) {
            if (M[i * n + j] != 0) {
                sparse_M[i * max_num_of_nonzero_elements_on_each_row + nonzero_elements_on_this_row] = M[i * n + j];
                indecies[i * max_num_of_nonzero_elements_on_each_row + nonzero_elements_on_this_row] = j;
                ++nonzero_elements_on_this_row;
            }
        }
    }
    return std::tuple {std::move(indecies), std::move(sparse_M)};
}

template <typename T>
auto transpose_matrix(
    const std::vector<T>& M,
    const size_t n
) {
    assert(M.size() == n * n);
    
    auto M_T = std::vector<T>(M.size());

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            M_T[i * n + j] = M[j * n + i];
        }
    }

    return M_T;
}

}
}
