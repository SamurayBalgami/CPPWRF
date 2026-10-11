#pragma once
// Working precision. CPPWRF_DOUBLE is set by CMake option CPPWRF_DOUBLE.
namespace cppwrf {
#ifdef CPPWRF_DOUBLE
using real = double;
inline constexpr const char* precision_name = "float64";
#else
using real = float;
inline constexpr const char* precision_name = "float32";
#endif
using real_t = real;
}  // namespace cppwrf
