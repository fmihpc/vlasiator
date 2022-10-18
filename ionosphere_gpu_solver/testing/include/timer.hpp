#pragma once

#include <chrono>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace timer {
class Manager {
  private:
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;
    typedef std::pair<time_point, std::string> time_point_with_info;

    std::vector<time_point_with_info> time_points {};

  public:
    Manager(std::string info) { push_time(info); };

    void push_time(std::string info = "") {
        time_points.push_back(std::make_pair(std::chrono::high_resolution_clock::now(), info));
    };


    template <bool print_pairs>
    void print_times() {
        if (time_points.size() % 2 != 0 && print_pairs) {
            throw std::runtime_error("Tried to print pairs with odd number of time_points");
        }

        if (time_points.size() < 2) {
            throw std::runtime_error("Timer needs at least 2 time_points");
        }

        std::cout << "\n Times:\n";

        size_t increment;
        if (print_pairs) {
            increment = 2;
        } else {
            increment = 1;
        }
        for (size_t i = 0; i < time_points.size() - 1; i += increment) {
            const auto duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(time_points[i + 1].first - time_points[i].first)
                    .count();

            std::cout << time_points[i].second << " --> " << time_points[i + 1].second << ": " << duration << "ms\n";
        }
        std::cout << "\n";
    }
};

} // namespace timer
