#!/usr/bin/env python3
"""The data tars each suite case reads (DATA in tests/verify/cases/*.py), derived from the tars' file lists.

  case_data.py LISTDIR [--write]

LISTDIR holds TAR.txt (`tar tvf TAR`, tools/colab/colab_tarlist.sh) for the tars of tests/verify/vm_data.txt (and of
the slow cases' tars, if listed). For every case of tests/verify/wrf_suite.py (and its special checks) the files it
reads are: the reference wrfout files; next to them the run's namelist.input and myoutfields.txt, for a tolerance-0
case the auxhist3 files and, restarted in the middle (suite.restart_at), WRF's wrfrst files at the restart; and
every input file (not wrfout / wrfrst / auxhist / rsl / log files) directly in each directory and each file its
configuration and --set list name under ${CPPWRF_WRF_RUNS} or ${CPPWRF_WRF_REF} (real.dir, phys.table_dir,
case.sounding_file). The tars are unpacked in vm_data.txt's order, a later tar's file replacing an earlier one's,
so a file is the last tar's version; an earlier tar that holds the same file (same size and time) serves as well.
The case needs the fewest tars (by size) that give each of its files as the full suite has it. A reference wrfout
file in no tar is an error. Prints each case whose recorded DATA differs from this and exits 1 if one does;
--write records the derived tars in the group files instead (DATA = {...} at the end of each file).
"""
import datetime
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wrf_suite  # noqa: E402
from cases import DATA, GROUPS  # noqa: E402
import cases.restarts as restarts  # noqa: E402

REPO = wrf_suite.REPO
OUTPUT = re.compile(r"(wrfout|wrfrst|auxhist\d+)_d\d\d_.*|rsl\..*|.*\.log")


def order():
    """vm_data.txt's tars in its order."""
    return [x.strip() for x in open(os.path.join(HERE, "vm_data.txt")) if x.strip() and not x.startswith("#")]


def listing(path):
    """{path in the tar: (size, date time)} of a `tar tvf` list (links: their target's entry is the file)."""
    out = {}
    for line in open(path):
        f = line.rstrip("\n").split(None, 5)
        if len(f) < 6 or not re.match(r"[-hl]", f[0]):
            continue
        p = f[5].split(" -> ")[0].split(" link to ")[0]
        p = p[2:] if p.startswith("./") else p
        out[p] = (int(f[2]), f[3] + " " + f[4])
    return out


def rel(path):
    """A path under ${CPPWRF_WRF_RUNS} / ${CPPWRF_WRF_REF} as it is in the tars (wrf_runs/..., wrf_ref/...)."""
    path = path.replace("${CPPWRF_WRF_RUNS}", wrf_suite.WR).replace("${CPPWRF_WRF_REF}", wrf_suite.WREF)
    for root, name in ((wrf_suite.WR, "wrf_runs"), (wrf_suite.WREF, "wrf_ref")):
        if path.startswith(root + "/"):
            return os.path.normpath(name + path[len(root):])
    return None


def config_paths(cfg, sets):
    """The data paths a configuration and its --set list name (section.key -> path; --set replaces the file's)."""
    conf, sec = {}, None
    for line in open(os.path.join(REPO, cfg)):
        line = line.split("#")[0].strip()
        if re.fullmatch(r"\[.+\]", line):
            sec = line[1:-1]
        elif sec and "=" in line:
            conf[f"{sec}.{line.split('=', 1)[0].strip()}"] = line.split("=", 1)[1].strip()
    for x in sets:
        if "=" in x:
            conf[x.split("=", 1)[0]] = x.split("=", 1)[1]
    return [r for r in (rel(v) for v in conf.values()) if r]


def needs(name):
    """(files the check needs: {path: required}, data paths that may be directories) of a case or special check."""
    files, dirs = {}, []
    case = next((c for c in wrf_suite.CASES if c[0] == name), None)
    if case:
        _, _, cfg, sets, wo, t, tol = case
        dirs += config_paths(cfg, [s for s in sets if not s.startswith("suite.")])
        mid = wrf_suite.restart_at(sets)
        for w in (wo.values() if isinstance(wo, dict) else [wo]):
            r = rel(w)
            files[r] = True
            d, dom = os.path.dirname(r), re.search(r"wrfout_(d\d\d)_", r).group(1)
            for f in ("namelist.input", "myoutfields.txt"):
                files.setdefault(os.path.join(d, f), False)
            if all(x == 0.0 for x in tol.values()):
                files.setdefault(os.path.join(d, f"auxhist3_{dom}_*"), False)
                if mid:  # the reference's start (its file name) + mid
                    t0 = datetime.datetime.strptime(r[-19:], "%Y-%m-%d_%H:%M:%S")
                    at = (t0 + datetime.timedelta(seconds=mid)).strftime("%Y-%m-%d_%H:%M:%S")
                    files[os.path.join(d, f"wrfrst_{dom}_{at}")] = True
        return files, dirs
    if name.startswith("straka_wrf"):
        files[rel(f"{wrf_suite.WREF}/run_straka_d1/wrfout_d01_0001-01-01_00:00:00")] = True
        return files, config_paths("configs/straka_wrf.ini", [])
    m = re.fullmatch(r"restart_(.+?)_(fp64|fp64_cuda|fp32_cuda)", name)
    tag, cfg, rd = next(r for r in wrf_suite.RESTARTS if r[0] == m.group(1))
    return files, config_paths(cfg, [])


def derive(name, lists, tars):
    """The tars `name` needs, in vm_data.txt's order, and the problems (a required file in no tar)."""
    files, dirs = needs(name)
    want, problems = {}, []  # path -> required
    for p, req in files.items():
        if "*" in p:  # the auxhist3 files
            for t in tars:
                for q in lists[t]:
                    if re.fullmatch(re.escape(p).replace(r"\*", ".*"), q):
                        want.setdefault(q, False)
        else:
            want[p] = want.get(p, False) or req
    for d in dirs:
        hit = False
        for t in tars:
            for q in lists[t]:
                if os.path.dirname(q) == d and not OUTPUT.fullmatch(os.path.basename(q)):
                    want.setdefault(q, False)
                    hit = True
                elif q == d:
                    want.setdefault(q, False)
                    hit = True
        if not hit:
            problems.append(f"{d}: in no tar")
    # each file: the tars that hold it as the last one does (the full suite's version)
    good = {}
    for p, req in want.items():
        have = [t for t in tars if p in lists[t]]
        if not have:
            if req:
                problems.append(f"{p}: in no tar")
            continue
        good[p] = [t for t in have if lists[t][p] == lists[have[-1]][p]]
    size = {t: sum(s for s, _ in lists[t].values()) for t in tars}
    chosen = set()
    for p in sorted(good, key=lambda p: len(good[p])):  # the files only one tar gives first
        if not chosen & set(good[p]):
            chosen.add(min(good[p], key=lambda t: size[t]))
    while True:  # a later chosen tar's other version of a file wins: then the last one's tar too
        fix = [p for p in good if [t for t in tars if t in chosen and p in lists[t]][-1] not in good[p]]
        if not fix:
            break
        chosen.add([t for t in tars if fix[0] in lists[t]][-1])
    return [t for t in tars if t in chosen], problems


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        sys.exit(__doc__)
    tars = order() + [t for v in DATA.values() for t in v.split(",") if t and t not in order()]
    tars = list(dict.fromkeys(tars))
    lists = {t: listing(os.path.join(args[0], t + ".txt")) for t in tars if os.path.exists(os.path.join(args[0], t + ".txt"))}
    unlisted = [t for t in tars if t not in lists]
    names = [c[0] for c in wrf_suite.CASES] + wrf_suite.SPECIAL
    derived, bad = {}, 0
    for n in names:
        rec = DATA.get(n)
        if rec and any(t in unlisted for t in rec.split(",")):
            print(f"unchecked {n}: {rec} (no list of {', '.join(t for t in rec.split(',') if t in unlisted)})")
            derived[n] = rec
            continue
        got, problems = derive(n, lists, [t for t in tars if t in lists])
        derived[n] = ",".join(got)
        for x in problems:
            print(f"PROBLEM {n}: {x}")
        if rec != derived[n]:
            bad += 1
            print(f"{'DIFFERS' if rec else 'MISSING'} {n}: recorded {rec} derived {derived[n]}")
    if "--write" in sys.argv:
        for g in GROUPS + [restarts]:
            mine = [c[0] for c in getattr(g, "CASES", [])] if g is not restarts else wrf_suite.SPECIAL
            path = g.__file__
            s = open(path).read()
            k = s.find("\n# the data tars (tests/verify/vm_data.txt)")
            s = s if k < 0 else s[:k + 1]
            s = s.rstrip("\n") + "\n\n" + (
                "# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py\n"
                "# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)\n"
                "DATA = {\n" + "".join(f"    {n!r}: {derived[n]!r},\n" for n in mine) + "}\n")
            open(path, "w").write(s)
        print(f"DATA written for {len(names)} checks")
        return
    print(f"{len(names)} checks, {bad} differ from the tars' lists; unlisted tars: {', '.join(unlisted) or 'none'}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
