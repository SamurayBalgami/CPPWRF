// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Lateral halo fill (periodic or symmetric) for x and, if ny>1, y. Halo width 3.
// Symmetric (rigid wall, mirror at the boundary face, WRF symmetric_xs/xe):
//   normal velocity (U in x, V in y): odd about the face, whose value is left as is;  everything else: even.
// Fill order x then y, so corners are consistent.
#include <initializer_list>
#include <vector>

#include "core/field.hpp"

namespace cppwrf {

// OddFace: odd about the face, whose value changes sign too (WRF set_physical_bc2d 'u' in x: dat(ids-i) =
// -dat(ids+i) from i = 0, i.e. muu at the walls; the 3-D 'u' and the 'v' loops start at i = 1)
enum class Parity { Even, Odd, OddFace };

void fill_halo(Field3& f, const Grid& g, Parity px, Parity py);
void fill_halo(Field2& f, const Grid& g, Parity px, Parity py);
// Default parity from staggering: U odd in x, V odd in y, else even.
void fill_halo(Field3& f, const Grid& g);
// 2D fields (mu, muu, ...) are always even.
void fill_halo(Field2& f, const Grid& g);
// Default-parity fill of several fields with one launch per direction (same result as one by one).
void fill_halos(std::initializer_list<Field3*> fs, const Grid& g);
void fill_halos(std::initializer_list<Field2*> fs, const Grid& g);
void fill_halos(const std::vector<Field3*>& fs, const Grid& g);
// One field of a batched fill: 3-D or 2-D, the parities of fill_halo(f, g) or the given ones.
struct HaloField {
  real* p;
  Layout L;
  Stagger st;
  Parity px, py;
  HaloField(Field3& f, Parity x, Parity y) : p(f.data()), L(f.layout()), st(f.stagger()), px(x), py(y) {}
  HaloField(Field2& f, Parity x, Parity y) : p(f.data()), L(f.layout()), st(f.stagger()), px(x), py(y) {}
  HaloField(Field3& f)  // NOLINT: implicit, as the field itself in a list
      : HaloField(f, f.stagger() == Stagger::U ? Parity::Odd : Parity::Even,
                  f.stagger() == Stagger::V ? Parity::Odd : Parity::Even) {}
  HaloField(Field2& f) : HaloField(f, Parity::Even, Parity::Even) {}  // NOLINT
};
// Several fields of either rank and any parities, one launch per direction for up to 8 fields: the same values
// as fill_halo of each field in turn (each field's x pass before its y pass; the fields must be distinct).
void fill_halos(std::initializer_list<HaloField> fs, const Grid& g);
void fill_halos(const std::vector<HaloField>& fs, const Grid& g);

}  // namespace cppwrf
