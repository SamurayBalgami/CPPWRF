#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Step 2.

- The radiation driver and Dudhia (which share RRTMG's strided views) also use rr.
- Floating literals in all these files become rr(literal): in C++ an unsuffixed literal is double and
  would pull float arithmetic back to double (Fortran's default-real literals of WRF real*4 are single).
  FP64 build: rr = double, rr(x) is the same value and type, so that build is unchanged.
- fortran_read.hpp: float overloads (read as double, then rounded) for the FP32 build.
"""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
phys = root / "src/phys"
INC = '#include "phys/ra_rrtmg_common/rrtmg_real.hpp"\n'

extra = [phys / "driver/radiation.cpp", phys / "driver/radiation.hpp"] + \
        sorted(p for p in (phys / "ra_sw_dudhia").iterdir() if p.suffix in (".cpp", ".hpp"))
for p in extra:
    t = p.read_text()
    t2, k = re.subn(r"\bdouble\b", "rr", t)
    if k and INC not in t2:
        if "#pragma once\n" in t2:
            t2 = t2.replace("#pragma once\n", "#pragma once\n" + INC, 1)
        else:
            i = t2.index("#include")
            t2 = t2[:i] + INC + t2[i:]
    p.write_text(t2)
    print(f"{p.relative_to(root)}: double -> rr {k}")

LIT = re.compile(r"(?<![\w.])((?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?|\d+[eE][+-]?\d+)(?![\w.])")
targets = [p for d in ("ra_rrtmg_lw", "ra_rrtmg_sw", "ra_rrtmg_common") for p in sorted((phys / d).iterdir())
           if p.suffix in (".hpp", ".cpp") and p.name not in ("fortran_read.hpp", "rrtmg_real.hpp")
           and not p.name.endswith("_f64.hpp") and not p.name.endswith("_tables.cpp")] + extra


def wrap_code(line):
    # leave comments, strings and preprocessor lines alone
    if line.lstrip().startswith("#"):
        return line, 0
    out, n, i = [], 0, 0
    code = line
    cpos = code.find("//")
    tail = ""
    if cpos >= 0:
        code, tail = code[:cpos], code[cpos:]
    parts = re.split(r'("(?:\\.|[^"\\])*")', code)
    for j, seg in enumerate(parts):
        if j % 2 == 1:
            out.append(seg)
            continue
        seg2, k = LIT.subn(r"rr(\1)", seg)
        n += k
        out.append(seg2)
    return "".join(out) + tail, n


total = 0
for p in targets:
    lines = p.read_text().split("\n")
    in_block = False
    cnt = 0
    for idx, line in enumerate(lines):
        if in_block:
            if "*/" in line:
                in_block = False
            continue
        if line.lstrip().startswith("/*") and "*/" not in line:
            in_block = True
            continue
        lines[idx], k = wrap_code(line)
        cnt += k
    p.write_text("\n".join(lines))
    total += cnt
    print(f"{p.relative_to(root)}: literals {cnt}")
print("total literals:", total)

fr = phys / "ra_rrtmg_common/fortran_read.hpp"
t = fr.read_text()
add = """
#ifndef CPPWRF_DOUBLE
// EXPERIMENT: float targets (FP32 build): read as double, rounded once
inline bool fortran_real(const std::string& s, int implied_d, float& v) {
  double t = 0; bool ok = fortran_real(s, implied_d, t); v = float(t); return ok; }
inline bool fortran_f_field(const std::string& rec, std::size_t a, std::size_t w, int d, float& v) {
  double t = 0; bool ok = fortran_f_field(rec, a, w, d, t); v = float(t); return ok; }
inline int fortran_list_real(const std::string& rec, float& v) {
  double t = 0; int st = fortran_list_real(rec, t); v = float(t); return st; }
#endif
"""
k = t.rfind("}  // namespace")
assert k > 0
fr.write_text(t[:k] + add + t[k:])
print("fortran_read overloads added")
