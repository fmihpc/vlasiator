#pragma once
#include "timer.hpp"
#include <concepts>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <utility>
#include <chrono> 
#include <algorithm>
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


}
}
