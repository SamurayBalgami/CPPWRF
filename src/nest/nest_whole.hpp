// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Nests in a decomposed run (docs/MPI.md): a worker's patch of the parent and its patch of the nest cover different
// areas, so a kernel that reads one domain and writes the other reads a copy of the whole source field
// (tile::WholeCopy, gathered from the workers) and writes the points of the destination's patch. Without a
// decomposition the copy is the field itself.
#include <functional>

#include "core/field.hpp"
#include "core/tiles.hpp"

namespace cppwrf::nest {

// A real field (FieldT, a non-uniform BaseFieldT) of a domain as a whole on this worker: read-only, with the
// members the nest kernels' launchers use (data, layout, cview). The constructor is collective in a decomposed run.
template <int R>
class WholeField {
 public:
  template <class F>
  explicit WholeField(const F& f) : data_(f.data()), L_(f.layout()) {
    w_.get(f.base(), f.layout());
    if (w_.on()) {
      data_ = w_.template data<real>();
      L_ = w_.layout();
    }
  }
  // only the box of the field each worker reads (tile::WholeCopy::get with a box)
  template <class F>
  WholeField(const F& f, const std::function<tile::PatchBox(int)>& box) : data_(f.data()), L_(f.layout()) {
    w_.get(f.base(), f.layout(), box);
    if (w_.on()) {
      data_ = w_.template data<real>();
      L_ = w_.layout();
    }
  }
  real* data() const { return data_; }
  const Layout& layout() const { return L_; }
  auto cview() const {
    if constexpr (R == 3) return View3<const real>{data_, int(L_.sj), int(L_.sk)};
    else return View2<const real>{data_, int(L_.sj)};
  }

 private:
  tile::WholeCopy w_;
  real* data_;
  Layout L_;
};

// The halo exchange of a field whose own points a nest kernel has just set (no-op without a decomposition)
template <class F>
void exchange_field(const F& f) {
  if (tile::tiled() && f.layout().si != 0) tile::exchange_base(f.base());
}

}  // namespace cppwrf::nest
