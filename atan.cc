#include <algorithm>
#include <iomanip>
#include <array>
#include <random>
#include <cmath>
#include <iostream>
#include <chrono>
#include "fast_discrete_arctangent.h"

namespace {
const int kNumComputations = 1e7;
using FloatT=float;

std::pair<FloatT, FloatT> kTestTable[kNumComputations];

void InitTestTable() {
  std::mt19937 generator;
  std::uniform_real_distribution<FloatT> dist(-5.f, 5.f);
  for (auto& point : kTestTable) {
    do {
      point.first = dist(generator);
      point.second = dist(generator);
    // Ensure we'll not attempt to compute sector of (0,0) point.
    } while (point.first == 0.f && point.second == 0.f);
  }
}

}  // namespace

template<typename T, typename FloatT>
void PrintResult(const T& calc, FloatT x, FloatT y) {
  const int n = calc.SectorNumer(x, y);
  std::cout << "x=" << x << ", y=" << y << " sector " << n << std::endl;
}


template<typename T>
int TestPerf(const T& calc, const std::string& name) {
  // Unused result to prevent compiler optimize out all SectorNumer() calls.
  int unused = 0;
  const auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < kNumComputations; ++i) {
    unused += calc.SectorNumer(kTestTable[i].first, kTestTable[i].second);
  }
  const auto end = std::chrono::high_resolution_clock::now();
  size_t ns_count = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
  std::cout << name << " took " << ns_count << " nanoseconds, " << static_cast<float>(ns_count) / kNumComputations << " ns/computation" << std::endl;
  return unused;
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
    int expected_sector = expected_calc.SectorNumer(
        point.first, point.second);
    int tested_sector = tested_calc.SectorNumer(
        point.first, point.second);
    const int diff = std::abs(expected_sector - tested_sector);
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

  const DiscreteAtanSimple<FloatT> slow_atan(kNumSectors);
  DoIt(slow_atan, "Slow computing:");

  const DiscreteAtanTableBased<FloatT> fast_atan(kNumSectors);
  DoIt(fast_atan, "Fast computing:");

  Compare(slow_atan, fast_atan);
  return 0;
}
