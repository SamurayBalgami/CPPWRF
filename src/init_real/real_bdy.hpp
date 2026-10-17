#pragma once
// Lateral boundary arrays of real.exe (main/real_em.F assemble_output): couple
// (dyn_em/module_big_step_utilities_em.F, with calc_mu_staggered), stuff_bdy_new and
// stuff_bdytend_new (share/module_bc.F).
#include <algorithm>

#include "init_real/farray.hpp"
#include "init_real/real_core.hpp"

namespace cppwrf::init_real {

// Boundary arrays: xs/xe (j, k, w), ys/ye (i, k, w); w = 1..spec_bdy_width
template <typename T>
struct BdyVar {
  F3<T> xs, xe, ys, ye;
  void alloc(const Dims& d, int nk_lo, int nk_hi, int w) {
    xs = F3<T>(d.jds, d.jde, nk_lo, nk_hi, 1, w);
    xe = F3<T>(d.jds, d.jde, nk_lo, nk_hi, 1, w);
    ys = F3<T>(d.ids, d.ide, nk_lo, nk_hi, 1, w);
    ye = F3<T>(d.ids, d.ide, nk_lo, nk_hi, 1, w);
  }
};

template <typename T>
void calc_mu_staggered(const F2<T>& mu, const F2<T>& mub, F2<T>& muu, F2<T>& muv, const Dims& d) {
  const int its = d.its(), ite = d.ite(), jts = d.jts(), jte = d.jte();
  int itf = ite, jtf = std::min(jte, d.jde - 1);
  // its = ids and ite = ide (serial)
  for (int j = jts; j <= jtf; ++j)
    for (int i = its + 1; i <= itf - 1; ++i)
      muu(i, j) = T(0.5) * (mu(i, j) + mu(i - 1, j) + mub(i, j) + mub(i - 1, j));
  for (int j = jts; j <= jtf; ++j) muu(its, j) = mu(its, j) + mub(its, j);
  for (int j = jts; j <= jtf; ++j) muu(ite, j) = mu(ite - 1, j) + mub(ite - 1, j);
  itf = std::min(ite, d.ide - 1);
  jtf = jte;
  for (int j = jts + 1; j <= jtf - 1; ++j)
    for (int i = its; i <= itf; ++i)
      muv(i, j) = T(0.5) * (mu(i, j) + mu(i, j - 1) + mub(i, j) + mub(i, j - 1));
  for (int i = its; i <= itf; ++i) muv(i, jts) = mu(i, jts) + mub(i, jts);
  for (int i = its; i <= itf; ++i) muv(i, jte) = mu(i, jte - 1) + mub(i, jte - 1);
}

// couple: name 'u', 'v', 'h' or 't' (anything else: mass point, k = kts..kde-1)
template <typename T>
void couple(const F2<T>& mu, const F2<T>& mub, F3<T>& rfield, const F3<T>& field, char name, const F2<T>& msf,
            const F1<T>& c1h, const F1<T>& c2h, const F1<T>& c1, const F1<T>& c2, const Dims& d) {
  const int its = d.its(), ite = d.ite(), jts = d.jts(), jte = d.jte(), kts = d.kts(), kte = d.kte();
  const int ktf = std::min(kte, d.kde - 1);
  if (name == 'u' || name == 'v') {
    F2<T> muu(d.ids, d.ide, d.jds, d.jde), muv(d.ids, d.ide, d.jds, d.jde);
    calc_mu_staggered(mu, mub, muu, muv, d);
    if (name == 'u') {
      for (int j = jts; j <= std::min(jte, d.jde - 1); ++j)
        for (int k = kts; k <= ktf; ++k)
          for (int i = its; i <= ite; ++i) rfield(i, k, j) = field(i, k, j) * (c1h(k) * muu(i, j) + c2h(k)) / msf(i, j);
    } else {
      for (int j = jts; j <= jte; ++j)
        for (int k = kts; k <= ktf; ++k)
          for (int i = its; i <= std::min(ite, d.ide - 1); ++i)
            rfield(i, k, j) = field(i, k, j) * (c1h(k) * muv(i, j) + c2h(k)) / msf(i, j);
    }
  } else if (name == 'h') {
    for (int j = jts; j <= std::min(jte, d.jde - 1); ++j)
      for (int k = kts; k <= kte; ++k)
        for (int i = its; i <= std::min(ite, d.ide - 1); ++i)
          rfield(i, k, j) = field(i, k, j) * ((c1(k) * mu(i, j)) + (c1(k) * mub(i, j) + c2(k)));
  } else {
    for (int j = jts; j <= std::min(jte, d.jde - 1); ++j)
      for (int k = kts; k <= ktf; ++k)
        for (int i = its; i <= std::min(ite, d.ide - 1); ++i)
          rfield(i, k, j) = field(i, k, j) * ((c1(k) * mu(i, j)) + (c1(k) * mub(i, j) + c2(k)));
  }
}

// stuff_bdy_new / stuff_bdytend_new: val(i, k, j) gives the value to store. kd0..kd1 = kds..kde of the field.
template <typename T, typename F>
void stuff(BdyVar<T>& b, char st, int w, int kd0, int kd1, const Dims& dd, F val) {
  const int ids = dd.ids, ide = dd.ide, jds = dd.jds, jde = dd.jde;
  const int its = ids, ite = ide, jts = jds, jte = jde;
  const bool full_k = (st == 'W' || st == 'M');
  const int k1 = full_k ? kd1 : kd1 - 1;
  // xs
  {
    const int j1 = st == 'V' ? std::min(jde, jte) : std::min(jde - 1, jte);
    for (int j = std::max(jds, jts); j <= j1; ++j)
      for (int k = kd0; k <= k1; ++k)
        for (int i = std::max(ids, its); i <= std::min(ids + w - 1, ite); ++i) b.xs(j, k, i) = val(i, k, j);
  }
  // xe
  if (st == 'U') {
    for (int j = std::max(jds, jts); j <= std::min(jde - 1, jte); ++j)
      for (int k = kd0; k <= kd1 - 1; ++k)
        for (int i = std::min(ide, ite); i >= std::max(ide - w + 1, its); --i) b.xe(j, k, ide - i + 1) = val(i, k, j);
  } else {
    const int j1 = st == 'V' ? std::min(jde, jte) : std::min(jde - 1, jte);
    for (int j = std::max(jds, jts); j <= j1; ++j)
      for (int k = kd0; k <= k1; ++k)
        for (int i = std::min(ide - 1, ite); i >= std::max(ide - w, its); --i) b.xe(j, k, ide - i) = val(i, k, j);
  }
  // ys
  {
    const int i1 = st == 'U' ? std::min(ide, ite) : std::min(ide - 1, ite);
    for (int j = std::max(jds, jts); j <= std::min(jds + w - 1, jte); ++j)
      for (int k = kd0; k <= k1; ++k)
        for (int i = std::max(ids, its); i <= i1; ++i) b.ys(i, k, j) = val(i, k, j);
  }
  // ye
  if (st == 'V') {
    for (int j = std::min(jde, jte); j >= std::max(jde - w + 1, jts); --j)
      for (int k = kd0; k <= kd1 - 1; ++k)
        for (int i = std::max(ids, its); i <= std::min(ide - 1, ite); ++i) b.ye(i, k, jde - j + 1) = val(i, k, j);
  } else {
    const int i1 = st == 'U' ? std::min(ide, ite) : std::min(ide - 1, ite);
    for (int j = std::min(jde - 1, jte); j >= std::max(jde - w, jts); --j)
      for (int k = kd0; k <= k1; ++k)
        for (int i = std::max(ids, its); i <= i1; ++i) b.ye(i, k, jde - j) = val(i, k, j);
  }
}

}  // namespace cppwrf::init_real
