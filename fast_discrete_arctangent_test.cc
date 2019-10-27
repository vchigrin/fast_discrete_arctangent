#include <gtest/gtest.h>
#include <iostream>
#include <chrono>
#include <random>
#include "fast_discrete_arctangent.h"

namespace {

int get_sector_diff(int a, int b, int sector_count) {
  int diff = std::abs(a - b);
  return std::min(diff, sector_count - diff);
}

}  // namespace

TEST(FastDiscreteArctangentTest, TestCorrectness) {
  static const int kNumSectors = 2000;
  static const int kNumComputations = 2e7;

  std::mt19937 generator;
  std::uniform_real_distribution<double> dist(-5.f, 5.f);

  // The most precise results expected during computations in double.
  const DiscreteAtanSimple<double> expected_calc(kNumSectors);
  const DiscreteAtanTableBased<double> tested_calc_double(kNumSectors);
  const DiscreteAtanTableBased<float> tested_calc_float(kNumSectors);
  const DiscreteAtanSimple<float> simple_calc_float(kNumSectors);

  int num_tested_double_rounding_errors = 0;
  int num_tested_float_rounding_errors = 0;
  int num_tested_simple_float_rounding_errors = 0;
  for (int i = 0; i < kNumComputations; ++i) {
    double x = dist(generator);
    double y = dist(generator);
    const int expected_sector = expected_calc.SectorNumer(x, y);
    const int tested_sector_double = tested_calc_double.SectorNumer(x, y);
    const int tested_sector_float = tested_calc_float.SectorNumer(x, y);
    const int tested_sector_simple_float = simple_calc_float.SectorNumer(x, y);

    {
      const int diff = get_sector_diff(
          expected_sector, tested_sector_double, kNumSectors);
      // Error greater then 1 is definitely not a rounding error and must
      // be considered as bug.
      ASSERT_LE(diff, 1);
      num_tested_double_rounding_errors += (diff > 0);
    }

    {
      const int diff = get_sector_diff(
          expected_sector, tested_sector_float, kNumSectors);
      ASSERT_LE(diff, 1);
      num_tested_float_rounding_errors += (diff > 0);
    }

    {
      const int diff = get_sector_diff(
          expected_sector, tested_sector_simple_float, kNumSectors);
      ASSERT_LE(diff, 1);
      num_tested_simple_float_rounding_errors += (diff > 0);
    }
  }
  // Expect less then 1  rounding error per 1 million
  // measurements in double mode.
  EXPECT_LE(
      num_tested_double_rounding_errors,
      kNumSectors / 1e6);
  // Expect calculations in table-based form in float mode is not worse then
  // simple calculations using floats.
  EXPECT_LE(num_tested_float_rounding_errors,
            num_tested_simple_float_rounding_errors);
  std::cout << "Performed " << kNumComputations
            << " measurements" << std::endl;
  std::cout << "Rounding error count in double mode "
            << num_tested_double_rounding_errors << std::endl;
  std::cout << "Rounding error count in float mode "
            << num_tested_float_rounding_errors << std::endl;
  std::cout << "Rounding error count in simple float mode "
            << num_tested_simple_float_rounding_errors << std::endl;
}


template<typename FloatT>
void TestSpeed() {
  static const int kNumSectors = 2000;
  static const int kNumComputations = 2e7;
  std::vector<std::pair<FloatT, FloatT>> test_table(kNumComputations);
  std::mt19937 generator;
  std::uniform_real_distribution<FloatT> dist(-5.f, 5.f);
  for (int i = 0; i < kNumComputations; ++i) {
    auto& point = test_table[i];
    do {
      point.first = dist(generator);
      point.second = dist(generator);
    // Ensure we'll not attempt to compute sector of (0,0) point.
    } while (point.first == 0 && point.second == 0);
  }
  DiscreteAtanSimple<FloatT> simple_calc(kNumSectors);
  DiscreteAtanTableBased<FloatT> fast_calc(kNumSectors);

  int64_t sum_simple = 0;
  std::chrono::nanoseconds duration_simple;
  {
    const auto start = std::chrono::steady_clock::now();
    for (const auto& point : test_table) {
      // Use unused sum to prevent compiler optimize out this call or
      // some it parts.
      sum_simple += simple_calc.SectorNumer(point.first, point.second);
    }
    const auto end = std::chrono::steady_clock::now();
    duration_simple = std::chrono::duration_cast<std::chrono::nanoseconds>
        (end - start);
  }

  int64_t sum_fast = 0;
  std::chrono::nanoseconds duration_fast;
  {
    const auto start = std::chrono::steady_clock::now();
    for (const auto& point : test_table) {
      // Use unused sum to prevent compiler optimize out this call or
      // some it parts.
      sum_fast += fast_calc.SectorNumer(point.first, point.second);
    }
    const auto end = std::chrono::steady_clock::now();
    duration_fast = std::chrono::duration_cast<std::chrono::nanoseconds>
        (end - start);
  }
  EXPECT_LE(std::abs(sum_fast - sum_simple), kNumComputations);
  const double speed_up_percent = (
      (duration_simple - duration_fast).count() * 100.) /
          duration_simple.count();
  std::cout << "Table based implementation is " << speed_up_percent
            << "% faster then simple" << std::endl;
  std::cout << "Simple duration "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                duration_simple).count() << " ms. fast "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                duration_fast).count() << " ms."
            << " for " << kNumComputations << " computations" << std::endl;
  // Expect that fast implementation at least 25% faster.
  // Typical numbers are ~50% (float) and ~70% (double).
  ASSERT_LE(duration_fast.count(), duration_simple.count() * 0.75);
}

TEST(FastDiscreteArctangentTest, TestSpeedDouble) {
  TestSpeed<double>();
}

TEST(FastDiscreteArctangentTest, TestSpeedFloat) {
  TestSpeed<float>();

}
