#pragma once
#include <array>
#include <cassert>
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

template<typename FloatT = float>
class DiscreteAtanTableBased {
 public:
  DiscreteAtanTableBased(int sector_count)
      : sector_count_(sector_count),
        indices_count_(static_cast<int>(
            std::ceil(1. / (std::tan((2 * M_PI) / sector_count))))) {
    assert(sector_count > 8);
    assert(sector_count % 8 == 0);
    InitTables();
  }

  size_t SectorNumer(FloatT x, FloatT y) const {
    // Compute result as offset + SectorNumer0_45(x, y);
    int offset = 0;
    if (y < 0) {
      offset += sector_count_ / 2;
      x = -x;
      y = -y;
    }
    // Now x,y in 0...pi range.
    if (x < 0) {
      offset += sector_count_ / 4;
      const FloatT new_y = -x;
      x = y;
      y = new_y;
    }
    // Now x,y in 0...pi / 2 range.
    if (x <= y) {
      offset += sector_count_ / 8;
      const FloatT new_x = x + y;
      const FloatT new_y = y - x;
      x = new_x;
      y = new_y;
    }
    // Now x,y in 0...pi / 4 range.
    return offset + SectorNumer0_45(x, y);
  }

 private:
  // Computes sector number for angles [0.. pi/4] (that is x >= y, x>=0, y>=0).
  inline size_t SectorNumer0_45(FloatT x, FloatT y) const {
    const FloatT tan_angle = y / x;
    const auto& entry = table_[static_cast<int>(tan_angle * indices_count_)];
    int r = entry.j_value;
    return r + (entry.t_value <= tan_angle);
  }

  void InitTables() {
    std::vector<FloatT> table_pi_4;
    table_pi_4.resize(sector_count_ / 8 + 1);
    for (size_t k = 0; k < sector_count_ / 8; ++k) {
      const FloatT alpha = static_cast<FloatT>((2. * M_PI * k) / sector_count_);
      table_pi_4[k] = std::tan(alpha);
    }
    // Last, out-of-range element to make search code simpler.
    table_pi_4[sector_count_ / 8] = table_pi_4[sector_count_ / 8 - 1] + 1;
    size_t start = 0;
    table_.resize(indices_count_ + 1);
    for (size_t k = 0; k <= indices_count_; ++k) {
      const FloatT value = static_cast<FloatT>(k) / indices_count_;
      while (table_pi_4[start] < value) {
        start++;
      }
      table_[k].j_value = start - 1;
      table_[k].t_value = table_pi_4[start];
    }
  }
  const int sector_count_;
  const int indices_count_;

  struct Entry {
    FloatT t_value;
    // Value of J[k] - 1 from the document.
    int j_value;
  };

  std::vector<Entry> table_;
};
