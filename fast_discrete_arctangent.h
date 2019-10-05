#pragma once
#include <array>
#include <cmath>


// May be useful with large SectorCount - in that case in is faster
// then DiscreteAtanTableBased.
template<size_t SectorCount, typename FloatT = float>
class DiscreteAtanSimple {
  static constexpr FloatT kInvSectorNumber =
      SectorCount / static_cast<FloatT>(2. * M_PI);
 public:
  size_t SectorNumer(FloatT x, FloatT y) const {
    // -pi ... pi
    FloatT alpha = std::atan2(y, x);
    if (alpha < 0.f) {
      // 0 ... 2. * pi
      alpha += static_cast<FloatT>(2. * M_PI);
    }
    return static_cast<FloatT>(alpha * kInvSectorNumber);
  }
};

template<size_t SectorCount, typename FloatT = float>
class DiscreteAtanTableBased {
 public:
  static_assert(SectorCount > 4, "Too few sectors");
  static_assert(SectorCount % 8 == 0, "Sector count must be divisible by 8");

  DiscreteAtanTableBased() {
    InitTables();
  }

  size_t SectorNumer(FloatT x, FloatT y) const {
    // Compute result as base + multiplier * SectorNumer0_45(x, y);
    int multiplier = 1;
    int base = 0;
    if (y < 0) {
      // return SectorCount - SectorNumer_0_180(x, -y) - 1;
      base = SectorCount - 1;
      multiplier = -multiplier;
      y = -y;
    }
    // Now x,y in 0...pi range.
    if (x < 0) {
      // return SectorCount / 2 - SectorNumer0_90(-x, y) - 1;
      base += (SectorCount / 2 - 1) * multiplier;
      multiplier = -multiplier;
      x = -x;
    }
    // Now x,y in 0...pi / 2 range.
    if (x <= y) {
      // -1 since SectorNumer0_45 rounds angle to zero, we want preserve
      // this  behavior here.
      // return SectorCount / 4 - SectorNumer0_45(y, x) - 1;
      base += (SectorCount / 4 - 1) * multiplier;
      multiplier = -multiplier;
      std::swap(x, y);
    }
    // Now x,y in 0...pi / 4 range.
    return base + multiplier * SectorNumer0_45(x, y);
  }

 private:
  static constexpr size_t kIndicesCount = SectorCount / 4;

  // Computes sector number for angles [0.. pi/4] (that is x >= y, x>=0, y>=0).
  inline size_t SectorNumer0_45(FloatT x, FloatT y) const {
    const FloatT tan_angle = y / x;
    size_t upper_bound = upper_bound_indices_[
        static_cast<size_t>(tan_angle * kIndicesCount)];
    /*
    Since we round down (tan_angle * kIndicesCount), we get upper bound
    not greater then required.
    We can not start with "too big" position here.
    PROOF:
    By definition how we build upper_bound_indices_, and because entries
    in table_pi_4_ monotonically increase, (as does tan function on [0...pi/4],
    we have for each integral i
    table_pi_4_[upper_bound_indices_[i]] >= i / kIdxCount
    table_pi_4_[upper_bound_indices_[i] - 1] < i / kIdxCount

    Here we take i = floor(tan_angle * kIndicesCount). Fromm here and
    from last inequality we have

    table_pi_4_[upper_bound_indices_[i] - 1] < floor(tan * kIdxCount) / kIdxCount
    and
    floor(tan * kIdxCount) / kIdxCount <= tan * kIdxCount / kIdxCount = tan

    So we've proved that previous entry in table_pi_4_ - k
    table_pi_4_[upper_bound_indices_[i] - 1] < tan
    that we required.
    */

    // Find exact upper_bound for current angle (it may be less then required).
    const FloatT* p = table_pi_4_.begin() + upper_bound;
    while (*p <= tan_angle)
      ++p;

    return p - table_pi_4_.begin() - 1;
  }

  void InitTables() {
    for (size_t k = 0; k < SectorCount / 8; ++k) {
      const FloatT alpha = static_cast<FloatT>((2. * M_PI * k) / SectorCount);
      table_pi_4_[k] = std::tan(alpha);
    }
    // Last, out-of-range element to make search code simpler.
    table_pi_4_[SectorCount / 8] = table_pi_4_[SectorCount / 8 - 1] + 1;
    size_t start = 0;
    for (size_t k = 0; k < kIndicesCount; ++k) {
      const FloatT value = static_cast<FloatT>(k) / kIndicesCount;
      while (table_pi_4_[start] < value) {
        start++;
      }
      upper_bound_indices_[k] = start;
    }
    upper_bound_indices_[kIndicesCount] =
      upper_bound_indices_[kIndicesCount - 1];
  }
  // Computed for sectors in range [0..pi/4]. All other can be obtained by
  // mirroring.
  // i-th element of the table contain lower boundary of the tan(alpha),
  // that go to the sector with number i.
  std::array<FloatT, SectorCount / 8 + 1> table_pi_4_;
  // i-th element contain lowest index in table_pi_4_, that value greater
  // then or equal to  i / kIndicesCount
  std::array<int, kIndicesCount + 1> upper_bound_indices_;
};
