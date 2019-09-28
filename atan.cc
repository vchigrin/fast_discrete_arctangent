#include <algorithm>
#include <iomanip>
#include <array>
#include <random>
#include <cmath>
#include <iostream>
#include <chrono>
#include "fast_discrete_arctangent.h"

namespace {
const int kNumComputations = 10000000;
std::pair<float, float> kTestTable[kNumComputations];

void InitTestTable() {
  std::mt19937 generator;
 // std::uniform_real_distribution<float> dist(-5.f, 5.f);
  std::uniform_real_distribution<float> dist(0.f, 5.f);
  for (auto& point : kTestTable) {
    do {
      point.first = dist(generator);
      point.second = dist(generator);
    // Ensure we'll not attempt to compute sector of (0,0) point.
    } while (point.first == 0.f && point.second == 0.f);
  }
}

}  // namespace

template<typename T>
void PrintResult(const T& calc, double x, double y) {
  const size_t n = calc.SectorNumer(x, y);
  std::cout << "x=" << x << ", y=" << y << " sector " << n << std::endl;
}


template<typename T>
void TestPerf(const T& calc, const std::string& name) {

  const auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < kNumComputations; ++i) {
    calc.SectorNumer(kTestTable[i].first, kTestTable[i].second);
  }
  const auto end = std::chrono::high_resolution_clock::now();
  size_t ns_count = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
  std::cout << name << " took " << ns_count << " nanoseconds, " << static_cast<float>(ns_count) / kNumComputations << " ns/computation" << std::endl;
}

template<typename T>
void DoIt(const T& calc,  const std::string& name) {
  std::cout << "\n==== Using " << name << "=========" << std::endl;

  PrintResult(calc, 1.f, 1.f);
  PrintResult(calc, 2.f, 1.f);
  PrintResult(calc, 3.f, 1.f);
  PrintResult(calc, 4.f, 1.f);
  PrintResult(calc, 5.f, 1.f);
  std::cout << std::endl;

  PrintResult(calc, 1.f, 0.f);
  PrintResult(calc, 1.f, 1.f);
  PrintResult(calc, 0.f, 1.f);
  PrintResult(calc, -1.f, 1.f);
  PrintResult(calc, -1.f, 0.f);
  PrintResult(calc, -1.f, -1.f);
  PrintResult(calc, 0.f, -1.f);
  PrintResult(calc, 1.f, -1.f);


  TestPerf(calc, name);
}

template<typename T1, typename T2>
void Compare(const T1& expected_calc, const T2& tested_calc) {
  int num_errors = 0;
  int num_by_one_errors = 0;
  for (int i = 0; i < kNumComputations; ++i) {
    const auto& point = kTestTable[i];
    size_t expected_sector = expected_calc.SectorNumer(
        point.first, point.second);
    size_t tested_sector = tested_calc.SectorNumer(
        point.first, point.second);
    const int diff = std::abs(static_cast<int>(expected_sector) -
        static_cast<int>(tested_sector));
    if (diff == 1) {
      num_by_one_errors++;
      ////
//      std::cout << "Point " << std::setprecision(11) << point.first << ";" << point.second
//                << " expected sector " << expected_sector
//                << " tested sector " << tested_sector
//                << std::endl;
//      break;
      ////
    } else if (diff > 1) {
      std::cerr << "Error on point " << point.first
                << ";" << point.second << std::endl;
      return;
      num_errors++;
    }
  }
  std::cout << "Testing finished "  << num_errors
            << " big errors and "
            << num_by_one_errors << " by one errors found among "
            << kNumComputations << " points"  << std::endl;
}

int main(int, char*[]) {
  InitTestTable();
  static const int kNumSectors = 2200;

  const DiscreteAtanSimple<kNumSectors, double> slow_atan;
  DoIt(slow_atan, "Slow computing:");

  const DiscreteAtanTableBased<kNumSectors, double> fast_atan;
  DoIt(fast_atan, "Fast computing:");

  Compare(slow_atan, fast_atan);
 // std::cout << "Avg dist " << ((float)total_dist) /  total_calls << std::endl;
  /*
  Linear search, with [x,4/3 range]
  Fast computing: took 371566842 nanoseconds, 37.1567 ns/computation
  Testing finished 0 big errors and 96 by one errors found among 10000000 points
  Avg dist 17.1494


  Binary search, with [x,4/3 range]
  Fast computing: took 453782378 nanoseconds, 45.3782 ns/computation
  Testing finished 0 big errors and 96 by one errors found among 10000000 points
  Avg dist 35.6263

  Separate index table of size  700
Fast computing: took 167781529 nanoseconds, 16.7782 ns/computation
Testing finished 0 big errors and 96 by one errors found among 10000000 points
Avg dist 1.20894


Fast computing: took 175751151 nanoseconds, 17.5751 ns/computation
Testing finished 0 big errors and 96 by one errors found among 10000000 points
Avg dist 1.25375
  */

  return 0;
}
