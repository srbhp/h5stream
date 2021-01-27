#pragma once
#include <ctime>
#include <functional>
#include <random>

  std::random_device rd;
  std::mt19937 gen(rd());

double unirandom(double a, double b) {
  std::uniform_real_distribution<double> distribution(a, b);
  return distribution(gen);
}
