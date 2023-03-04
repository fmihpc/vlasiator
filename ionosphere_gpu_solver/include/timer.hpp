#pragma once

#include <chrono>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <optional>
#include <iomanip>
#include <ios>


namespace ionogpu {

namespace timer {
/* 
    We assume that every function in form of

        constexpr bool condition = ...;

        void foo() {
            if constexpr (condition) {

            }
        }

    will get optimized away ie. function is never called if condition is false.

    In compiler explorer this was true for nvcc with atleast -O2

 */

constexpr bool enable_benchmark = 0;
namespace timer_internal {
    
    using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;
    struct TimeInfo {
         std::chrono::microseconds total_time;
        std::optional<TimePoint> last_time_point;
        size_t number_of_measurements;
    };

    struct TimePrinter {
        std::map<std::string, TimeInfo> times;

        void key_found(TimeInfo& time_info) {

            if (time_info.last_time_point.has_value()) {
                const auto now = std::chrono::steady_clock::now();
                time_info.total_time += std::chrono::duration_cast<decltype(TimeInfo::total_time)>(now - time_info.last_time_point.value());
                time_info.number_of_measurements += 1;
                time_info.last_time_point.reset();
            } else {
                time_info.last_time_point = std::chrono::steady_clock::now();
            }
        }

        TimePrinter() = default; 
        ~TimePrinter() {
            std::cout <<
            "│  Name                      │ Cumulative time[μs]  │  Number of measurements  │ Average time[μs] │\n";
            const auto sorted_times_in_vec = [&] {
                auto temp = std::vector<std::pair<std::string, TimeInfo>>();
                temp.reserve(times.size());
                for (const auto& x: times) {
                    temp.push_back(x);
                }
                std::sort(temp.begin(), temp.end(), [](const auto& a, const auto& b) {return std::get<1>(a).total_time > std::get<1>(b).total_time; });
                return temp;
            }();

            for (const auto& [name, time_info] : sorted_times_in_vec) {
                std::cout << "│ " << std::setw(26) << std::left << name << " │ " 
                << std::setw(20) << std::right << time_info.total_time.count() << " │ " 
                << std::setw(24) << std::right << time_info.number_of_measurements << " │ " 
                << std::setw(16) << std::right << static_cast<double>(time_info.total_time.count()) / static_cast<double>(time_info.number_of_measurements) << " │\n";
            }
        }


        TimePrinter(const TimePrinter& other) = delete;
        TimePrinter(TimePrinter&& other) = delete;
        TimePrinter& operator=(const TimePrinter& other) = delete;
        TimePrinter& operator=(TimePrinter&& other) = delete;
    };
}

// First time this is called with spesific key it will start timer for that key
// Second time this is called with spesific key it will stop the timer for that key and then adds to total time for that key
// When program terminates the destructor of TimePrinter will print all of the total times for all of the keys
inline void time(const std::string& key) {
    if constexpr (enable_benchmark) {
        static auto time_printer = timer_internal::TimePrinter{};

        if (auto time_info_iterator = time_printer.times.find(key); time_info_iterator == time_printer.times.end()) {
            time_printer.times[key] = timer_internal::TimeInfo{std::chrono::milliseconds{0}, std::chrono::steady_clock::now(), 0};
        } else {
            time_printer.key_found(std::get<1>(*time_info_iterator));
        }
    }
}

// Overload of the above that calls it twice so
//
// timer::time("fucntion 1");
// function1();              
// timer::time("function 1", "function 2");  
// function2();              
// timer::time("function 2");
//
inline void time(const std::string& key1, const std::string& key2) {
    time(key1);
    time(key2);
}


} // namespace timer
} // namespace ionogpu
