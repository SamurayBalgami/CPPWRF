#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Step 3: the two call sites outside RRTMG."""
import sys
from pathlib import Path

root = Path(sys.argv[1])

# the old RRTM LW scheme keeps double: not used by the measured cases, so in the FP32 experiment build it stops
p = root / "src/phys/driver/radiation.cpp"
t = p.read_text()
a = t.index("        const RrtmLwColumn rc{")
b = t.index("rrtm_lw_column(*p_->rrtm, rc, *rws, g_lw, olr, rlw.data(), rlwc.data(), rmsg);")
b = t.index("\n", b) + 1
t = (t[:a] + "#ifdef CPPWRF_DOUBLE\n" + t[a:b] +
     "#else\n        std::fprintf(stderr, \"EXPERIMENT: RRTM LW (ra_lw_physics 1) is not in this FP32 build\\n\");\n"
     "        std::abort();\n#endif\n" + t[b:])
if "#include <cstdlib>" not in t:
    t = t.replace("#include <cstdio>", "#include <cstdio>\n#include <cstdlib>", 1)
p.write_text(t)

# the suite hands the latitudes / longitudes as double: convert to the driver's type
p = root / "src/phys/driver/suite.cpp"
t = p.read_text()
old = "rad_ = std::make_unique<RadiationDriver>(g_, o_.rad, xlat, xlong);"
assert old in t
t = t.replace(old, "rad_ = std::make_unique<RadiationDriver>(g_, o_.rad, std::vector<::cppwrf::real>(xlat.begin(), xlat.end()),\n"
                   "                                             std::vector<::cppwrf::real>(xlong.begin(), xlong.end()));")
p.write_text(t)
print("ok")
