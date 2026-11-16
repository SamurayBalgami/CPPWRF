#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Step 5: RRTM LW (ra_lw_physics = 1) in the FP32 build. rrtm_lw.cpp keeps double (as in the FP64 build); the
driver widens the columns it hands over and rounds the results once, instead of the stops of steps 3 and 4."""
import sys
from pathlib import Path

p = Path(sys.argv[1]) / "src/phys/driver/radiation.cpp"
s = p.read_text()

def stop_block(s, indent):
    """the '#else ... stop ... #endif' block of steps 3 / 4 whose fprintf line has this indent"""
    a = s.index("#else\n" + indent + 'std::fprintf(stderr, "EXPERIMENT: RRTM LW')
    return s[a:s.index("#endif\n", a) + len("#endif\n")]


stop_batch = stop_block(s, " " * 6)
new_batch = """#else
      // FP32 build: RRTM LW computes in double (rrtm_lw.cpp is not converted); the staged columns are widened and
      // the tendencies and fluxes rounded once
      const std::size_t n3 = std::size_t(ncol) * nz, n3w = std::size_t(ncol) * (nz + 1), nc = std::size_t(ncol);
      auto host = [](const rr* d, std::size_t n) {
        std::vector<rr> h(n);
        copy_to_host(h.data(), d, n * sizeof(rr));
        return std::vector<double>(h.begin(), h.end());
      };
      auto dev = [](rr* d, const std::vector<double>& h) {
        const std::vector<rr> f(h.begin(), h.end());
        copy_to_dev(d, f.data(), f.size() * sizeof(rr));
      };
      const std::vector<double> hp8w = host(q.p8w, n3w), ht8w = host(q.t8w, n3w), hp = host(q.p3d, n3),
                                ht3 = host(q.t3d, n3), hdz = host(q.dz, n3), hpi = host(q.pi3d, n3),
                                hqv = host(q.qv, n3), hcf = host(q.cf, n3), htsk = host(q.tsk, nc),
                                hem = host(q.emiss, nc);
      const std::vector<double> hqc = in.qc ? host(q.qc, n3) : std::vector<double>(),
                                hqr = in.qr ? host(q.qr, n3) : std::vector<double>(),
                                hqi = in.qi ? host(q.qi, n3) : std::vector<double>(),
                                hqs = in.qs ? host(q.qs, n3) : std::vector<double>();
      const RrtmLwBatchInput bi{ncol, nz, hp8w.data(), ht8w.data(), hp.data(), ht3.data(), hdz.data(), hpi.data(),
                                hqv.data(), in.qc ? hqc.data() : nullptr, in.qr ? hqr.data() : nullptr,
                                in.qi ? hqi.data() : nullptr, in.qs ? hqs.data() : nullptr, hcf.data(), htsk.data(),
                                hem.data(), in.qi != nullptr, lc.co2, lc.n2o, lc.ch4, cst::r_d, cst::g, o_.icloud,
                                o_.warm_rain, li.col_begin, li.col_end};
      std::vector<double> hr(n3), hrc(n3), hg(nc), ho(nc);
      rrtm_lw_batch(*q.rrtm, bi, q.rrtm_nlayers, hr.data(), hrc.data(), hg.data(), ho.data());
      dev(q.rlw, hr);
      dev(q.rlwc, hrc);
      dev(q.glw, hg);
      dev(q.olr, ho);
      const std::vector<rr> zero(8 * nc, rr(0.0));
      copy_to_dev(q.lwd, zero.data(), 8 * nc * sizeof(rr));
      copy_to_dev(q.lwcf, zero.data(), nc * sizeof(rr));
#endif
"""
stop_col = stop_block(s, " " * 8)
new_col = """#else
        // FP32 build: RRTM LW computes in double (rrtm_lw.cpp is not converted); the column is widened and the
        // results rounded once
        auto dv = [](const std::vector<rr>& v) { return std::vector<double>(v.begin(), v.end()); };
        const std::vector<double> dp8w = dv(cp8w), dt8w = dv(ct8w), dp = dv(cp), dt = dv(ct), ddz = dv(cdz),
                                  dpi = dv(cpi), dqv = dv(cqv), dqc = dv(cqc), dqr = dv(cqr), dqi = dv(cqi),
                                  dqs = dv(cqs), dcf = dv(ccf);
        const RrtmLwColumn rc{dp8w.data(), dt8w.data(), dp.data(), dt.data(), ddz.data(), dpi.data(), dqv.data(),
                              in.qc ? dqc.data() : nullptr, in.qr ? dqr.data() : nullptr,
                              in.qi ? dqi.data() : nullptr, in.qs ? dqs.data() : nullptr, dcf.data(),
                              in.qi != nullptr, tsk[q2], ems[q2], lc.co2, lc.n2o, lc.ch4, cst::r_d,
                              cst::g, o_.icloud, o_.warm_rain};
        std::vector<double> drlw(nz), drlwc(nz);
        double dglw = 0, dolr = 0;
        rrtm_lw_column(*p_->rrtm, rc, *rws, dglw, dolr, drlw.data(), drlwc.data(), rmsg);
        g_lw = rr(dglw);
        olr = rr(dolr);
        for (int k = 0; k < nz; ++k) {
          rlw[k] = rr(drlw[k]);
          rlwc[k] = rr(drlwc[k]);
        }
#endif
"""
assert s.count(stop_batch) == 1 and s.count(stop_col) == 1, (s.count(stop_batch), s.count(stop_col))
s = s.replace(stop_batch, new_batch).replace(stop_col, new_col)
p.write_text(s)
print("ok")

# rr: float only in CUDA FP32 builds. The gain is a GPU one; on the CPU float was not faster (N3, 16 threads), so
# CPU builds keep double and CPU FP32 is the same arithmetic as before this change.
root = Path(sys.argv[1]) / "src/phys/ra_rrtmg_common"
(root / "rrtmg_real.hpp").write_text("""// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// RRTMG's arithmetic type (src/phys/ra_rrtmg_lw, ra_rrtmg_sw, ra_rrtmg_common, the radiation driver, Dudhia SW).
// CUDA FP32 builds: float, as RRTMG's rb in WRF real*4; the FP64 rate of most GPUs is 1/32 to 1/64 of the FP32 rate,
// and RRTMG in double dominated the FP32 step with radt = 1 (docs/RRTMG_FP32_EXPERIMENT.md). CPU builds and FP64
// builds: double, the same code as with double written out (rr(x) is x).
#include "core/precision.hpp"
namespace cppwrf::phys {
#if defined(CPPWRF_USE_CUDA)
using rr = ::cppwrf::real;
#else
using rr = double;
#endif
}  // namespace cppwrf::phys
""")
t = (root / "fortran_read.hpp").read_text()
a = "// EXPERIMENT: float targets (FP32 build): read as double, rounded once"
assert a in t
t = t.replace(a, "// float targets (CUDA FP32 builds, RRTMG's tables): read as double, rounded once")
(root / "fortran_read.hpp").write_text(t)
print("rr ok")

# the suite hands the latitudes / longitudes to the driver in its arithmetic type (step 3 used cppwrf::real)
p = Path(sys.argv[1]) / "src/phys/driver/suite.cpp"
t = p.read_text()
assert t.count("std::vector<::cppwrf::real>(xl") == 2
t = t.replace("std::vector<::cppwrf::real>(xl", "std::vector<::cppwrf::phys::rr>(xl")
if '#include "phys/ra_rrtmg_common/rrtmg_real.hpp"' not in t:
    i = t.index("#include")
    t = t[:i] + '#include "phys/ra_rrtmg_common/rrtmg_real.hpp"\n' + t[i:]
p.write_text(t)
print("suite ok")
