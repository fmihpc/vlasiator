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
namespace iono_gpu {
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
 */
template <std::ranges::range R>
void calculate_numerical_error_of_range(const R& v, const R& v_correct) {
    auto total_absolute_error = double{ 0.0 };
    auto total_relative_error = double{ 0.0 };
    
    if (v.size() != v_correct.size()) {
        std::cout << "v (" << v.size() << ") and v_correct ("<< v_correct.size() << ") contain different amount of elements!";
        return;
    }
    
    // Should be as follow in c++23  
    //for (const auto [x, x_correct] : std::ranges::views::zip(v, v_correct)) {
    for (size_t i = 0; i < v.size(); ++i) {
        const auto absolute_error = static_cast<double>(std::fabs(v[i] - v_correct[i]));
        const auto relative_error = absolute_error / static_cast<double>(v_correct[i]);

        total_absolute_error += absolute_error;
        total_relative_error += relative_error;
    }

    const auto average_relative_error = total_relative_error / static_cast<double>(v.size());
    
    std::cout << "Total absolute error: " << total_absolute_error << " Average relative error per element: " << average_relative_error << "\n";
}


}
}
