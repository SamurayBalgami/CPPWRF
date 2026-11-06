// SPDX-License-Identifier: GPL-3.0-or-later
// glibc libmvec's 2-lane vector log and log10 (_ZGVbN2v_log, _ZGVbN2v_log10) as they are on an Intel CPU, for the
// Fortran drivers of the differential tests (tests/phys/wrf_fortran.sh links this object before -lmvec, so the
// vectorised loops of WRF's objects call these instead of libmvec's). glibc's two functions take an approximate
// reciprocal with rcpps, whose result the hardware defines: on an AMD CPU libmvec's log / log10 differ from Intel's
// in the last bit for some inputs, and WRF's side of a test with them. Each lane here is cppwrf::mvec::log / log10
// (src/core/libmvec_x86.hpp: libmvec's SSE4 code statement by statement with Intel's rcpps table; on Intel equal to
// libmvec for every input, tools/colab/colab_mvec.sh), so WRF's code gives on every x86-64 CPU the bits it gives on
// Intel. On an AMD EPYC, libmvec's log made KF, Noah-MP and the MYNN surface layer differ, its log10 the WSM6 and
// Thompson reflectivity. Of glibc's other SSE4 vector functions in WRF's objects, acos, asin, tan and the float
// acosf, asinf, tanf, powf also take rsqrtps / rcpps (exp, pow, sin, cos, atan ... none); the tests do not depend on
// them (every differential test gives Intel's bits on AMD with only these two replaced), a new test that does
// would need them here too. Compiled without FMA contraction and without -march (IEEE double operations in
// libmvec's order). -DCPPWRF_MVEC_INTEL_NO_LOG / -DCPPWRF_MVEC_INTEL_NO_LOG10 leave one of them out (to find which one a
// difference comes from, e.g. as an LD_PRELOAD library).
#include <emmintrin.h>

#include "core/libmvec_x86.hpp"

extern "C" {
#if !defined(CPPWRF_MVEC_INTEL_NO_LOG)
__m128d _ZGVbN2v_log(__m128d x) {
  alignas(16) double v[2];
  _mm_store_pd(v, x);
  v[0] = cppwrf::mvec::log(v[0]);
  v[1] = cppwrf::mvec::log(v[1]);
  return _mm_load_pd(v);
}
#endif
#if !defined(CPPWRF_MVEC_INTEL_NO_LOG10)
__m128d _ZGVbN2v_log10(__m128d x) {
  alignas(16) double v[2];
  _mm_store_pd(v, x);
  v[0] = cppwrf::mvec::log10(v[0]);
  v[1] = cppwrf::mvec::log10(v[1]);
  return _mm_load_pd(v);
}
#endif
}
