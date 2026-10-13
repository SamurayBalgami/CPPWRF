#pragma once
// Lateral halo fill (periodic or symmetric) for x and, if ny>1, y. Halo width 3.
// Symmetric (rigid wall, mirror at the boundary face, WRF symmetric_xs/xe):
//   normal velocity (U in x, V in y): odd, face value = 0;  everything else: even.
// Fill order x then y, so corners are consistent.
#include <initializer_list>

#include "core/field.hpp"

namespace cppwrf {

enum class Parity { Even, Odd };

void fill_halo(Field3& f, const Grid& g, Parity px, Parity py);
void fill_halo(Field2& f, const Grid& g, Parity px, Parity py);
// Default parity from staggering: U odd in x, V odd in y, else even.
void fill_halo(Field3& f, const Grid& g);
// 2D fields (mu, muu, ...) are always even.
void fill_halo(Field2& f, const Grid& g);
// Default-parity fill of several fields with one launch per direction (same result as one by one).
void fill_halos(std::initializer_list<Field3*> fs, const Grid& g);
void fill_halos(std::initializer_list<Field2*> fs, const Grid& g);

}  // namespace cppwrf
