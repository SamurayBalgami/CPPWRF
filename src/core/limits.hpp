// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Host-side guards of the compile-time column sizes (kMaxLevels, soil and snow layers): an error in every build
// type (an assert would vanish with NDEBUG in Release builds and the kernels would overrun their column arrays).
#include <stdexcept>
#include <string>

namespace cppwrf {

inline void check_limit(const char* scheme, const char* what, int n, int max) {
  if (n > max)
    throw std::runtime_error(std::string(scheme) + ": " + what + " = " + std::to_string(n) + " exceeds the column limit " +
                             std::to_string(max));
}

}  // namespace cppwrf
