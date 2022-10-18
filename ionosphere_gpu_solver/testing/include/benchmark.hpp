#pragma once
#include "timer.hpp"
#include <concepts>
#include <functional>
#include <initializer_list>
#include <iostream>
namespace iono_gpu {
namespace testing {

template <typename... Args, std::invocable<Args...> F>
void benchmark_functions_with_parameters(const std::initializer_list<F> functions, const Args... arguments) {
    for (const auto f : functions) {
        std::cout << std::invoke(f, arguments...) << " ";
    }
    std::cout << "\n";
}


}
}
