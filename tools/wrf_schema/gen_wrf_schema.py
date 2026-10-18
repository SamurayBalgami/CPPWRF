#!/usr/bin/env python3
"""Extract the WRF I/O schema of a compiled WRF tree into a tab-separated table.

WRF's Registry is compiled into Fortran includes (inc/allocs_*.f90: the state variables with their
I/O metadata; inc/in_use_for_config_*.inc: the package conditions; inc/scalar_indices.inc: the
members of the 4-D arrays such as moist and scalar). This script reads them and writes one line per
state variable and per 4-D array member, so that cppwrf_real writes wrfinput/wrfbdy files with the
same variables, dimensions and attributes as WRF's real.exe.

usage: gen_wrf_schema.py <WRF build dir> <output .tsv>
"""
import glob
import os
import re
import sys


def parse_conditions(inc_dir):
    """in_use_for_config: var -> list of (option, op, value); in use if any holds (no entry: always)."""
    cond = {}
    pat_name = re.compile(r"IF\(TRIM\(vname\)\.EQ\.'([^']+)'\)THEN")
    pat_expr = re.compile(r"model_config_rec%(\w+)(?:\(id\))?\.(EQ|NE|GT|GE|LT|LE)\.([^\s)]+)")
    for f in sorted(glob.glob(os.path.join(inc_dir, "in_use_for_config_*.inc"))):
        cur = None
        for line in open(f):
            m = pat_name.search(line)
            if m:
                cur = m.group(1)
                continue
            if cur and "in_use = " in line and "in_use.OR." not in line:
                e = pat_expr.search(line)
                if not e:
                    raise SystemExit("unparsed condition for %s: %s" % (cur, line.strip()))
                cond.setdefault(cur, []).append((e.group(1), e.group(2), e.group(3)))
    return cond


def parse_allocs(inc_dir):
    files = glob.glob(os.path.join(inc_dir, "allocs_*.f90")) or glob.glob(os.path.join(inc_dir, "allocs_*.F"))
    files.sort(key=lambda p: int(re.search(r"allocs_(\d+)", p).group(1)))
    out = []
    pat = re.compile(r"grid%tail_statevars%(\w+(?:\(\d\))?)\s*=\s*(.*)$")
    alloc_cond = re.compile(r"in_use_for_config\(id,'([^']+)'\)")
    for f in files:
        cur_alloc = None
        e = None
        for line in open(f):
            m = alloc_cond.search(line)
            if m and line.lstrip().startswith("IF("):
                cur_alloc = m.group(1)
            m = pat.search(line)
            if not m:
                continue
            key, val = m.group(1), m.group(2).strip()
            if key == "VarName":
                e = {"alloc": cur_alloc}
                out.append(e)
            if e is None or key in ("next",):
                continue
            if val.startswith("'"):
                val = val[1 : val.rindex("'")]
            elif "!" in val:  # allocs_*.F: trailing "! hex" comments
                val = val[: val.index("!")].strip()
            elif val.startswith("=>") or key.endswith("table") or key.startswith("rfield") or key.startswith("ifield"):
                pass
            e[key] = val
    return out


def parse_scalar_tables(inc_dir):
    """4-D array members: (array, member, package option, value, dname, desc, units, streams1, streams2, bdy)."""
    f = os.path.join(inc_dir, "scalar_indices.inc")
    out = []
    cond = None
    pkg = re.compile(r"IF\s*\(\s*model_config_rec%(\w+)\(idomain\)\s*==\s*([^)\s]+)\s*\)\s*THEN")
    tab = re.compile(r"(\w+)_(dname|desc|units|boundary)_table\(\s*idomain,\s*P_(\w+)\s*\)\s*=\s*(.*)$")
    st = re.compile(r"(\w+)_streams_table\(\s*idomain,\s*P_(\w+)\s*\)%stream\((\d)\)\s*=\s*(\d+)")
    cur = {}
    for line in open(f):
        m = pkg.search(line)
        if m:
            cond = (m.group(1), m.group(2))
            continue
        m = tab.search(line)
        if m:
            arr, what, mem, val = m.group(1), m.group(2), m.group(3), m.group(4).strip()
            if val.startswith("'"):
                val = val[1 : val.rindex("'")]
            elif "!" in val:  # allocs_*.F: trailing "! hex" comments
                val = val[: val.index("!")].strip()
            k = (arr, mem, cond)
            cur.setdefault(k, {})[what] = val
            continue
        m = st.search(line)
        if m:
            k = (m.group(1), m.group(2), cond)
            cur.setdefault(k, {})["streams" + m.group(3)] = m.group(4)
    for (arr, mem, c), d in cur.items():
        out.append((arr, mem, c, d))
    return out


def rconfig_defaults(reg_dir):
    """namelist option -> (type, group, default) from the Registry rconfig lines (Registry.EM and its
    includes; the first definition wins)."""
    d = {}
    files = [os.path.join(reg_dir, "Registry.EM_COMMON")] + sorted(glob.glob(os.path.join(reg_dir, "registry.*")))
    for f in files:
        if not os.path.exists(f):
            continue
        for line in open(f, errors="replace"):
            t = line.split()
            if len(t) >= 6 and t[0] == "rconfig":
                grp = t[3].split(",")[1] if "," in t[3] else t[3]
                d.setdefault(t[2].lower(), (t[1], grp, t[5]))
    return d


def nest_specs(reg_dir):
    """state variable (lower case) -> (down, up): the nest interpolation (d) and feedback (u) of the
    Registry io column, "" when absent, "interp_fcn" / "copy_fcn" for the default functions."""
    out = {}
    files = [os.path.join(reg_dir, "Registry.EM_COMMON")] + sorted(glob.glob(os.path.join(reg_dir, "registry.*")))
    for f in files:
        if not os.path.exists(f):
            continue
        for line in open(f, errors="replace"):
            t = line.split()
            if len(t) < 8 or t[0].lower() != "state":
                continue
            io = t[7]
            spec = {}
            i = 0
            while i < len(io):
                ch = io[i]
                if ch in "duf" and i + 1 < len(io) and io[i + 1] == "=":
                    j = io.index(")", i)
                    spec[ch] = io[i + 3:j]
                    i = j + 1
                    continue
                if ch in "duf":
                    spec.setdefault(ch, {"d": "interp_fcn", "u": "copy_fcn", "f": "bdy_interp"}[ch])
                if ch == "{":
                    i = io.index("}", i) + 1
                    continue
                i += 1
            out.setdefault(t[2].lower(), (spec.get("d", ""), spec.get("u", "")))
    return out


def main():
    wrf, outp = sys.argv[1], sys.argv[2]
    inc = os.path.join(wrf, "inc")
    cond = parse_conditions(inc)
    allocs = parse_allocs(inc)
    scal = parse_scalar_tables(inc)
    defaults = rconfig_defaults(os.path.join(wrf, "Registry"))
    nests = nest_specs(os.path.join(wrf, "Registry"))

    def nest_of(v):  # time-level variables (u_2) carry the Registry name (u)
        v = v.lower()
        if v in nests:
            return nests[v]
        if re.match(r".*_[12]$", v) and v[:-2] in nests:
            return nests[v[:-2]]
        return ("", "")
    used_opts = sorted({c[0] for v in cond.values() for c in v} | {c[2][0] for c in scal if c[2]})
    with open(outp, "w") as o:
        o.write("# WRF I/O schema extracted by tools/wrf_schema/gen_wrf_schema.py from %s\n" % os.path.basename(os.path.normpath(wrf)))
        o.write("# RCONFIG\tname\ttype\tgroup\tdefault\n")
        for op in sorted(defaults):
            t, grp, dv = defaults[op]
            o.write("RCONFIG\t%s\t%s\t%s\t%s\n" % (op, t, grp, dv))
        missing = [op for op in used_opts if op.lower() not in defaults]
        for op in missing:
            o.write("RCONFIG\t%s\tinteger\tderived\t0\n" % op.lower())
        o.write("# VAR\tvarname\tdataname\ttype\tndim\tmemorder\tstagger\tntl\tstreams1\tstreams2\tdims(sd:ed;...)\tdimnames\tdesc\tunits\tconds\tnest_down\tnest_up\n")
        for e in allocs:
            v = e.get("VarName", "")
            if e.get("scalar_array", "").upper() == ".TRUE.":
                kind = "ARRAY4D"
            else:
                kind = "VAR"
            dims = ";".join("%s:%s" % (e.get("sd%d" % i, ""), e.get("ed%d" % i, "")) for i in (1, 2, 3, 4) if ("sd%d" % i) in e)
            dimnames = ";".join(e.get("dimname%d" % i, "") for i in (1, 2, 3, 4) if ("dimname%d" % i) in e)
            cs = cond.get(v, [])
            c = "|".join("%s.%s.%s" % x for x in cs) if cs else "-"
            o.write("\t".join([kind, v, e.get("DataName", ""), e.get("Type", ""), e.get("Ndim", ""),
                               e.get("MemoryOrder", ""), e.get("Stagger", ""), e.get("Ntl", "0"),
                               e.get("streams(1)", "0"), e.get("streams(2)", "0"), dims, dimnames,
                               e.get("Description", ""), e.get("Units", ""), c,
                               nest_of(v)[0], nest_of(v)[1]]) + "\n")
        o.write("# MEMBER\tarray\tmember\tcond\tdname\tdesc\tunits\tstreams1\tstreams2\tboundary\n")
        for arr, mem, c, d in scal:
            cs = "%s.EQ.%s" % c if c else "-"
            o.write("\t".join(["MEMBER", arr, mem, cs, d.get("dname", ""), d.get("desc", ""), d.get("units", ""),
                               d.get("streams1", "0"), d.get("streams2", "0"), d.get("boundary", "")]) + "\n")


if __name__ == "__main__":
    main()
