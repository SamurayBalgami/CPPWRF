#!/usr/bin/env python3
"""DATA statements of WRF's phys/module_ra_rrtm.F (RRTM longwave, ra_lw_physics = 1) as a C++ header.

  python3 tools/rrtm_lw_data.py WRF/phys/module_ra_rrtm.F > src/phys/ra_rrtm_lw/rrtm_lw_data.hpp

Every DATA statement of the module (its specification part) and of its subroutines becomes a constexpr array
(or scalar) in Fortran element order (column-major, the first index fastest); a subroutine's DATA is named
SUBROUTINE_NAME. The forms the module uses: NAME /list/, NAME(i) /v/, (NAME(IDATA[,j]),IDATA=lo,hi) /list/,
several of them in one statement, and repeat counts n*v. The values keep their Fortran spelling, which in the
real*8 build (-fdefault-real-8) are double constants: the C++ compiler rounds them as gfortran does. A module
array must be filled completely (checked against its declaration).
"""
import re
import sys

PARAMS = {"MG": 16, "NBANDS": 16, "NGPT": 140, "MAXINPX": 35, "MAXXSEC": 4}
INTS = {"NGM", "NGC", "NGS", "NGN", "NGB", "NXMOL", "IXINDX", "NG", "NSPA", "NSPB"}


def statements(path):
    """(scope, statement) for every free-form statement, comments removed and continuations joined."""
    scope, buf = "MODULE", ""
    for raw in open(path, encoding="latin-1"):
        line = raw.split("!", 1)[0].rstrip()
        if not line.strip():
            continue
        if buf:
            line = line.lstrip()
            if line.startswith("&"):
                line = line[1:]
        if line.endswith("&"):
            buf += line[:-1]
            continue
        st, buf = buf + line, ""
        m = re.match(r"\s*SUBROUTINE\s+(\w+)", st, re.I)
        if m:
            scope = m.group(1).upper()
        yield scope, st.strip()


def ev(s):
    s = s.strip().upper()
    for k, v in PARAMS.items():
        s = re.sub(r"\b%s\b" % k, str(v), s)
    return int(eval(s))  # small integer expressions: 5*(59-13+1)


def module_dims(path):
    """Declared bounds [(lo, hi), ...] of the module's arrays."""
    dims = {}
    for scope, st in statements(path):
        if scope != "MODULE":
            continue
        m = re.match(r"INTEGER\s*,\s*PARAMETER\s*::\s*(\w+)\s*=\s*(\d+)$", st, re.I)
        if m:
            PARAMS[m.group(1).upper()] = int(m.group(2))
            continue
        if re.match(r"(REAL|INTEGER)\b.*::", st, re.I):
            decl = st.split("::", 1)[1]
            for name, d in re.findall(r"(\w+)\s*\(([^()]*)\)", decl):
                b = []
                for part in d.split(","):
                    lo, hi = part.split(":") if ":" in part else ("1", part)
                    b.append((ev(lo), ev(hi)))
                dims[name.upper()] = b
    return dims


def values(lst):
    out = []
    for tok in lst.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "*" in tok:
            n, v = tok.split("*")
            out += [v.strip()] * int(n)
        else:
            out.append(tok)
    return out


def data_items(st):
    """[(target, [values])] of one DATA statement."""
    body = st[4:].strip()
    items = []
    while body:
        a = body.index("/")
        b = body.index("/", a + 1)
        items.append((body[:a].strip().strip(",").strip(), values(body[a + 1:b])))
        body = body[b + 1:].strip().lstrip(",").strip()
    return items


def cnum(v, is_int):
    v = v.strip()
    if is_int:
        assert re.fullmatch(r"[-+]?\d+", v), v
        return v
    assert re.fullmatch(r"[-+]?(\d+\.?\d*|\.\d+)([EeDd][-+]?\d+)?", v), v
    return v.replace("D", "E").replace("d", "e")


def main():
    path = sys.argv[1]
    dims = module_dims(path)
    arrays = {}  # name -> (bounds, {flat index: value}) for the module; locals: name -> [values]
    order = []
    for scope, st in statements(path):
        if not re.match(r"DATA\b", st, re.I):
            continue
        for target, vals in data_items(st):
            t = target.upper().replace(" ", "")
            m = re.fullmatch(r"\((\w+)\(IDATA(?:,(\d+))?\),IDATA=(\d+),(\d+)\)", t)
            if m:
                name, j, lo, hi = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
                idx = [(i, int(j)) if j else (i,) for i in range(lo, hi + 1)]
            else:
                m = re.fullmatch(r"(\w+)(?:\((\d+)\))?", t)
                assert m, target
                name = m.group(1)
                idx = [(int(m.group(2)),)] if m.group(2) else None
            key = name if scope == "MODULE" else scope + "_" + name
            if key not in arrays:
                arrays[key] = {}
                order.append((key, scope, name))
            if idx is None:
                if scope == "MODULE" and name in dims:
                    b = dims[name]
                    n = 1
                    for lo, hi in b:
                        n *= hi - lo + 1
                    assert len(vals) == n, (name, len(vals), n)
                for k, v in enumerate(vals):
                    arrays[key][k] = v
            else:
                b = dims[name]
                assert len(idx) == len(vals), (target, len(idx), len(vals))
                for ix, v in zip(idx, vals):
                    flat, mul = 0, 1
                    for (lo, hi), i in zip(b, ix):
                        assert lo <= i <= hi
                        flat += (i - lo) * mul
                        mul *= hi - lo + 1
                    assert flat not in arrays[key], (name, ix)
                    arrays[key][flat] = v
    print("#pragma once")
    print("// Generated by tools/rrtm_lw_data.py from WRF 4.8.0 phys/module_ra_rrtm.F: its DATA statements, in Fortran")
    print("// element order (first index fastest); a subroutine's DATA as SUBROUTINE_NAME. Do not edit.")
    print()
    print("namespace cppwrf::phys::rrtm_data {")
    print()
    for key, scope, name in order:
        a = arrays[key]
        is_int = scope == "MODULE" and name in INTS
        ty = "int" if is_int else "double"
        if scope == "MODULE" and name in dims:
            b = dims[name]
            n = 1
            for lo, hi in b:
                n *= hi - lo + 1
            assert sorted(a) == list(range(n)), (name, len(a), n)
            decl = ", ".join("%d:%d" % lb if lb[0] != 1 else str(lb[1]) for lb in b)
            comment = "  // %s(%s)" % (name, decl)
        else:
            n = len(a)
            assert sorted(a) == list(range(n)), key
            comment = ""
        vals = [cnum(a[k], is_int) for k in range(n)]
        if n == 1 and not (scope == "MODULE" and name in dims):
            print("inline constexpr %s %s = %s;" % (ty, key, vals[0]))
            continue
        print("inline constexpr %s %s[%d] = {%s" % (ty, key, n, comment))
        for k in range(0, n, 8):
            print("    " + ", ".join(vals[k:k + 8]) + ("," if k + 8 < n else ""))
        print("};")
    print()
    print("}  // namespace cppwrf::phys::rrtm_data")


if __name__ == "__main__":
    main()
