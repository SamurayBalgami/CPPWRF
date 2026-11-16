#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Step 4: the GPU build's RRTM LW host path (radiation.cpp, batch staging) keeps double, as in step 3."""
import sys
from pathlib import Path

p = Path(sys.argv[1]) / "src/phys/driver/radiation.cpp"
s = p.read_text()
a = "    if (rrtm) {\n      // RRTMLWRAD on the host from the staged columns"
b = "      copy_to_dev(q.lwcf, zero.data(), nc * sizeof(rr));\n    } else {\n      q.lwb->run("
assert s.count(a) == 1 and s.count(b) == 1
s = s.replace(a, "    if (rrtm) {\n#ifdef CPPWRF_DOUBLE\n      // RRTMLWRAD on the host from the staged columns")
s = s.replace(b, "      copy_to_dev(q.lwcf, zero.data(), nc * sizeof(rr));\n#else\n"
                 "      std::fprintf(stderr, \"EXPERIMENT: RRTM LW (ra_lw_physics 1) is not in this FP32 build" + chr(92) + "n\"); "
                 "std::abort();\n#endif\n    } else {\n      q.lwb->run(")
p.write_text(s)
print("ok")
