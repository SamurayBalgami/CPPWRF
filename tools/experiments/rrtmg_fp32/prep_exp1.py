#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Step 1: RRTMG arithmetic in the working precision.

In the exp clone: every `double` in src/phys/ra_rrtmg_{lw,sw,common} becomes `rr`
(= cppwrf::real: double in FP64 builds, so the FP64 build is the same code; float in FP32 builds).
The binary tables stay double on disk: they are read into a double copy of the structs
(namespace f64) and converted member by member.
"""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
phys = root / "src/phys"
INC = '#include "phys/ra_rrtmg_common/rrtmg_real.hpp"\n'
SKIP = {"fortran_read.hpp"}

(phys / "ra_rrtmg_common/rrtmg_real.hpp").write_text(
    "// SPDX-License-Identifier: GPL-3.0-or-later\n#pragma once\n"
    "// EXPERIMENT: RRTMG's arithmetic type = the build's working precision.\n"
    '#include "core/precision.hpp"\n'
    "namespace cppwrf::phys {\nusing rr = ::cppwrf::real;\n}\n")

# 1) keep the double definition of the table structs (file layout) in namespace f64
tables = {
    "lw": (phys / "ra_rrtmg_lw/rrtmg_lw_tables.hpp", "RrtmgLwTables", "load_rrtmg_lw_tables"),
    "sw": (phys / "ra_rrtmg_sw/rrtmg_sw_tables.hpp", "RrtmgSwTables", "load_rrtmg_sw_tables"),
}
struct_members = {}
for key, (hpp, top, loader) in tables.items():
    src = hpp.read_text()
    f64 = src.replace("namespace cppwrf::phys {", "namespace cppwrf::phys::f64 {")
    f64 = re.sub(r"bool %s\(.*?\);\n" % loader, "", f64)
    (hpp.parent / (hpp.stem + "_f64.hpp")).write_text(f64)
    # parse structs: name -> [(type, member)]
    structs = {}
    for m in re.finditer(r"struct (\w+) \{(.*?)\n\};", src, re.S):
        mem = []
        for line in m.group(2).splitlines():
            line = line.split("//")[0].strip()
            if not line:
                continue
            mm = re.fullmatch(r"(\w+)\s+(\w+)(\[[^\]]+\])?;", line)
            assert mm, (hpp, line)
            mem.append((mm.group(1), mm.group(2), mm.group(3) is not None))
        structs[m.group(1)] = mem
    struct_members[key] = (structs, top)

# 2) double -> rr in the RRTMG sources
files = [p for d in ("ra_rrtmg_lw", "ra_rrtmg_sw", "ra_rrtmg_common") for p in sorted((phys / d).iterdir())
         if p.suffix in (".hpp", ".cpp") and p.name not in SKIP and p.name != "rrtmg_real.hpp"
         and not p.name.endswith("_f64.hpp")]
n = 0
for p in files:
    t = p.read_text()
    t2, k = re.subn(r"\bdouble\b", "rr", t)
    if k == 0:
        continue
    n += k
    if "#pragma once\n" in t2:
        t2 = t2.replace("#pragma once\n", "#pragma once\n" + INC, 1)
    else:
        i = t2.index("#include")
        t2 = t2[:i] + INC + t2[i:]
    p.write_text(t2)
    print(f"{p.relative_to(root)}: {k}")
print("total double -> rr:", n)

# 3) table headers: no size assert for rr structs; converter + loader through the f64 copy
for key, (hpp, top, loader) in tables.items():
    t = hpp.read_text()
    t = re.sub(r"static_assert\(sizeof\(%s\).*?\n" % top, "", t)
    hpp.write_text(t)
    structs, _ = struct_members[key]
    conv = []
    for name, mem in structs.items():
        body = []
        for typ, mname, arr in mem:
            if typ in structs:
                body.append(f"  convert_tables(s.{mname}, d.{mname});")
            elif arr:
                body.append(f"  for (std::size_t i = 0; i < std::size(s.{mname}); ++i) d.{mname}[i] = "
                            f"static_cast<std::remove_reference_t<decltype(d.{mname}[0])>>(s.{mname}[i]);")
            else:
                body.append(f"  d.{mname} = static_cast<decltype(d.{mname})>(s.{mname});")
        conv.append(f"[[maybe_unused]] static void convert_tables(const f64::{name}& s, {name}& d) {{\n"
                    + "\n".join(body) + "\n}\n")
    cpp = hpp.with_suffix(".cpp")
    c = cpp.read_text()
    c = c.replace(f'#include "phys/ra_rrtmg_{key}/rrtmg_{key}_tables.hpp"\n',
                  f'#include "phys/ra_rrtmg_{key}/rrtmg_{key}_tables.hpp"\n'
                  f'#include "phys/ra_rrtmg_{key}/rrtmg_{key}_tables_f64.hpp"\n#include <iterator>\n#include <memory>\n'
                  f'#include <type_traits>\n')
    c = c.replace("namespace cppwrf::phys {", "namespace cppwrf::phys {\n" + "\n".join(conv), 1)
    old_read = f"file.read(reinterpret_cast<char*>(&tables), sizeof({top}));"
    old_chk = f"if (file.gcount() == static_cast<std::streamsize>(sizeof({top}))) return true;"
    assert old_read in c and old_chk in c, cpp
    c = c.replace(old_read, f"auto f = std::make_unique<f64::{top}>();\n"
                            f"    file.read(reinterpret_cast<char*>(f.get()), sizeof(f64::{top}));")
    c = c.replace(old_chk, f"if (file.gcount() == static_cast<std::streamsize>(sizeof(f64::{top}))) {{\n"
                           f"      convert_tables(*f, tables);\n      return true;\n    }}")
    cpp.write_text(c)
print("tables done")
