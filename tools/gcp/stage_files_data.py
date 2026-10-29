#!/usr/bin/env python3
"""Pack the references of tests/verify/run_wrf_files_test.sh and tests/verify/ideal_files_test.py whole (every frame
and variable) for the VM suite (tools/gcp/remote_suite.sh runs both tests after wrf_suite.py).

  stage_files_data.py [--list] [--links-only] OUT.tar|-   (-: the tar to stdout, e.g. | gcloud storage cp - gs://BUCKET/X.tar)

From ${CPPWRF_WRF_RUNS}: for each case of run_wrf_files_test.sh the input directory of its configuration (real.dir;
inputs and tables as tools/colab/stage_data.py keeps them), its wrfout reference directory and its restart (rst_) and
resumed (rsm_) directories (their wrfout / wrfrst files, namelist.input and iofields file; the resumed run's link to
the restart file it started from is left out); the steps/ references of wrf_suite.py's tolerance-0 cases whole
(stage_data.py --no-steps leaves them out of its tar); from ${CPPWRF_WRF_REF}: the reference directories of the ideal
cases (wrfout, namelist.input, input_sounding). Links are followed. --list prints the files and their size instead.
--links-only packs only the files that are links (their content): v05_files_links.tar, after v05_files_data.tar in
tests/verify/vm_data.txt, replaces the links that tar packed as links before the dereference fix.
"""
import os
import re
import sys
import tarfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "tools", "colab"))
sys.path.insert(0, os.path.join(REPO, "tests", "verify"))
import stage_data  # noqa: E402
import wrf_suite  # noqa: E402

WR, WREF = stage_data.WR, stage_data.WREF
OUTPUTS = re.compile(r"(wrfout|wrfrst)_d\d\d_.*|namelist\.input|myoutfields\.txt$")


def files_test_cases():
    """(config, wrfout reference dir, restart dir or None, extra --set) of run_wrf_files_test.sh's CASES."""
    text = open(os.path.join(REPO, "tests", "verify", "run_wrf_files_test.sh")).read()
    for line in text.split('CASES="', 1)[1].split('\n"', 1)[0].strip().splitlines():
        _, cfg, ref, rst, extra = line.split("|")
        yield cfg, ref.replace("$R", WR), None if rst == "-" else rst.replace("$R", WR), extra.replace("$R", WR)


def main():
    args = sys.argv[1:]
    listing = "--list" in args
    links_only = "--links-only" in args
    out = [a for a in args if a not in ("--list", "--links-only")][0]
    files = {}  # path in the tar -> file

    def add(path, root, arc):
        files.setdefault(os.path.join(arc, os.path.relpath(path, root)), path)

    def add_outputs(d, links=True):
        for f in sorted(os.listdir(d)):
            p = os.path.join(d, f)
            if OUTPUTS.match(f) and os.path.isfile(p) and (links or not os.path.islink(p)):
                add(p, WR, "wrf_runs")

    for cfg, ref, rst, extra in files_test_cases():
        m = re.search(r"real\.dir=(\S+)", extra)
        rd = m.group(1) if m else stage_data.config_value(cfg, "real", "dir").replace("${CPPWRF_WRF_RUNS}", WR)
        thompson = (stage_data.config_value(cfg, "phys", "mp_physics") or "0").strip() == "8"
        for f in sorted(os.listdir(rd)):
            p = os.path.join(rd, f)
            if os.path.isdir(p) or (f in stage_data.THOMPSON and not thompson):
                continue
            if f in stage_data.THOMPSON or not any(stage_data.fnmatch.fnmatch(f, s) for s in stage_data.SKIP):
                add(p, WR, "wrf_runs")
        add_outputs(ref)
        if rst:
            add_outputs(rst)
            rsm = rst.replace("/rst_", "/rsm_")
            if os.path.isdir(rsm):
                add_outputs(rsm, links=False)
    for name, key, cfg, sets, wo, t, tol in wrf_suite.CASES:
        refs = list(wo.values()) if isinstance(wo, dict) else [wo]
        if not all(x == 0.0 for x in tol.values()):
            continue
        for r in refs:
            if r.startswith(WREF) and key == "cpu_dp" and os.path.exists(r):  # ideal_files_test.py
                d = os.path.dirname(r)
                for f in ("namelist.input", "input_sounding", os.path.basename(r)):
                    if os.path.exists(os.path.join(d, f)):
                        add(os.path.join(d, f), WREF, "wrf_ref")
            elif "/steps/" in r and os.path.exists(r):
                add_outputs(os.path.dirname(r))
    if links_only:
        files = {arc: p for arc, p in files.items() if os.path.islink(p)}
    if listing:
        tot = 0
        for arc, p in sorted(files.items()):
            s = os.path.getsize(p)
            tot += s
            print(f"{s:12d} {arc}")
        print(f"{len(files)} files, {tot / 2**30:.2f} GiB")
        return
    tf = tarfile.open(fileobj=sys.stdout.buffer, mode="w|", dereference=True) if out == "-" else \
        tarfile.open(out, "w", dereference=True)  # a link's file content, not the link
    with tf:
        for arc, p in sorted(files.items()):
            tf.add(p, arcname=arc, recursive=False)
    print(f"{out}: {len(files)} files", file=sys.stderr)


if __name__ == "__main__":
    main()
