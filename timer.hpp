#pragma once

#include <chrono>
#include <iostream>
#include <string>

class timer {
  std::string m_function_name;
  std::chrono::time_point<std::chrono::system_clock> m_start;

public:
  explicit timer(std::string start_string)
      : m_function_name(std::move(start_string)) {
    m_start = std::chrono::system_clock::now();
  }
  void reset() { m_start = std::chrono::system_clock::now(); }
  ~timer() {
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - m_start;

    std::cout << "Finished " << m_function_name << " function in  "
              << elapsed_seconds.count() << "s\n";
  }
};
