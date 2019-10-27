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
    int diff = std::abs(expected_sector - tested_sector);
    diff = std::min(diff, expected_calc.sector_count() - diff);
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
                << ";" << point.second  << "."
                << " expected " << expected_sector
                << " tested " << tested_sector
                << std::endl;
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
  TestPerf(slow_atan, "Slow computing:");

  const DiscreteAtanTableBased<FloatT> fast_atan(kNumSectors);
  TestPerf(fast_atan, "Fast computing:");

  std::cout << "Comparing Slow computing (double) and fast computing "
            << typeid(FloatT).name() << std::endl;
  Compare(DiscreteAtanSimple<double>(kNumSectors), fast_atan);
  std::cout << "Comparing Slow computing (double) and Slow computing (float) "
            << std::endl;
  Compare(DiscreteAtanSimple<double>(kNumSectors),
          DiscreteAtanSimple<float>(kNumSectors));
  return 0;
}
