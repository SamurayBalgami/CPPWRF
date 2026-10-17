#pragma once
// Arrays indexed like WRF's Fortran arrays: arbitrary lower bounds, (i,k,j) order with i fastest.
// real.exe is ported statement by statement, so the C++ keeps WRF's indices and loop bounds.
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cppwrf::init_real {

template <typename T>
class F1 {
 public:
  F1() = default;
  F1(int k0, int k1, T v = T(0)) : k0_(k0), k1_(k1), d_(std::size_t(k1 - k0 + 1), v) {}
  T& operator()(int k) { return d_[std::size_t(k - k0_)]; }
  const T& operator()(int k) const { return d_[std::size_t(k - k0_)]; }
  int lo() const { return k0_; }
  int hi() const { return k1_; }
  std::vector<T>& data() { return d_; }
  const std::vector<T>& data() const { return d_; }

 private:
  int k0_ = 1, k1_ = 0;
  std::vector<T> d_;
};

template <typename T>
class F2 {
 public:
  F2() = default;
  F2(int i0, int i1, int j0, int j1, T v = T(0))
      : i0_(i0), i1_(i1), j0_(j0), j1_(j1), ni_(i1 - i0 + 1), d_(std::size_t(i1 - i0 + 1) * (j1 - j0 + 1), v) {}
  T& operator()(int i, int j) { return d_[idx(i, j)]; }
  const T& operator()(int i, int j) const { return d_[idx(i, j)]; }
  int ilo() const { return i0_; }
  int ihi() const { return i1_; }
  int jlo() const { return j0_; }
  int jhi() const { return j1_; }
  bool empty() const { return d_.empty(); }
  std::vector<T>& data() { return d_; }
  const std::vector<T>& data() const { return d_; }

 private:
  std::size_t idx(int i, int j) const { return std::size_t(i - i0_) + std::size_t(ni_) * std::size_t(j - j0_); }
  int i0_ = 1, i1_ = 0, j0_ = 1, j1_ = 0, ni_ = 0;
  std::vector<T> d_;
};

template <typename T>
class F3 {
 public:
  F3() = default;
  F3(int i0, int i1, int k0, int k1, int j0, int j1, T v = T(0))
      : i0_(i0), i1_(i1), k0_(k0), k1_(k1), j0_(j0), j1_(j1), ni_(i1 - i0 + 1), nk_(k1 - k0 + 1),
        d_(std::size_t(i1 - i0 + 1) * (k1 - k0 + 1) * (j1 - j0 + 1), v) {}
  T& operator()(int i, int k, int j) { return d_[idx(i, k, j)]; }
  const T& operator()(int i, int k, int j) const { return d_[idx(i, k, j)]; }
  int ilo() const { return i0_; }
  int ihi() const { return i1_; }
  int klo() const { return k0_; }
  int khi() const { return k1_; }
  int jlo() const { return j0_; }
  int jhi() const { return j1_; }
  bool empty() const { return d_.empty(); }
  std::vector<T>& data() { return d_; }
  const std::vector<T>& data() const { return d_; }

 private:
  std::size_t idx(int i, int k, int j) const {
    return std::size_t(i - i0_) + std::size_t(ni_) * (std::size_t(k - k0_) + std::size_t(nk_) * std::size_t(j - j0_));
  }
  int i0_ = 1, i1_ = 0, k0_ = 1, k1_ = 0, j0_ = 1, j1_ = 0, ni_ = 0, nk_ = 0;
  std::vector<T> d_;
};

}  // namespace cppwrf::init_real
