#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The drop-in wrf.exe and the namelist checker (tools/cppwrf-wrf, tools/cppwrf-check; python/cppwrf/dropin.py) on a
tiny em_quarter_ss run (10 x 10 x 10 points, 5 steps of 12 s, a generated sounding): WRF's file names in the run
directory (history_outname, rst_outname), the rsl logs with WRF's timing lines and SUCCESS COMPLETE WRF, the exit
status and the FATAL message of an option without a port and of a failing model, a run resumed from its wrfrst file
with restart = .true., and the checker's exit status.

  CPPWRF_BIN=<build>/cppwrf CPPWRF_REAL_BIN=<build>/src/init_real/cppwrf_real CPPWRF_TEST_DIR=<dir> \
      python3 -m unittest tests/tools/test_dropin.py
Without CPPWRF_BIN (CPPWRF_REAL_BIN) the model (cppwrf_real) runs are skipped. cppwrf-realexe is checked on its
failures only (no met_em files here): the FATAL message, the exit status, no SUCCESS line."""
import os
import shutil
import subprocess
import sys
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
WRF = os.path.join(REPO, "tools", "cppwrf-wrf")
CHECK = os.path.join(REPO, "tools", "cppwrf-check")
REALEXE = os.path.join(REPO, "tools", "cppwrf-realexe")
BIN = os.environ.get("CPPWRF_BIN", "")
REAL_BIN = os.environ.get("CPPWRF_REAL_BIN", "")
WORK = os.environ.get("CPPWRF_TEST_DIR", os.path.join(REPO, "build", "dropin_test"))

NAMELIST = """&time_control
 run_minutes = 1, start_year = 0001, start_month = 01, start_day = 01, start_hour = 00, start_minute = 00,
 start_second = 00, end_year = 0001, end_month = 01, end_day = 01, end_hour = 01,
 history_interval_s = 24, frames_per_outfile = 2, restart = {restart}, restart_interval = 1,
 io_form_history = 2, io_form_restart = 2, {extra_tc}
/
&domains
 time_step = 12, max_dom = 1, e_we = 11, e_sn = 11, e_vert = 11, dx = 2000, dy = 2000, ztop = 20000,
/
&physics
 mp_physics = {mp}, ra_lw_physics = 0, ra_sw_physics = 0, sf_sfclay_physics = 0, sf_surface_physics = 0,
 bl_pbl_physics = 0, cu_physics = 0,
/
&dynamics
 hybrid_opt = 0, diff_opt = 1, km_opt = 1, damp_opt = 3, zdamp = 5000., dampcoef = 0.2, khdif = 500, kvdif = 500,
 time_step_sound = 6, non_hydrostatic = .true., moist_adv_opt = 1, scalar_adv_opt = 1, mix_full_fields = .true.,
/
&bdy_control
 periodic_x = .true., periodic_y = .true.,
/
&ideal
 ideal_case = 2,
/
"""


def sounding() -> str:
    """A conditionally unstable sounding in WRF's input_sounding format (surface p, theta, qv; then z, theta, qv, u, v)."""
    rows = ["1000.0 300.0 14.0"]
    for k in range(1, 41):
        z = 500.0 * k
        th = 300.0 + (3.0 * z / 1000.0 if z > 12000 else 1.0 * z / 1000.0)
        qv = max(14.0 * (1.0 - z / 12000.0), 0.0)
        rows.append(f"{z:.1f} {th:.3f} {qv:.4f} {min(z / 1000.0, 10.0) - 5.0:.3f} 0.0")
    return "\n".join(rows) + "\n"


def make_dir(name, restart=".false.", mp=1, extra_tc="", with_sounding=True):
    d = os.path.join(WORK, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    with open(os.path.join(d, "namelist.input"), "w") as f:
        f.write(NAMELIST.format(restart=restart, mp=mp, extra_tc=extra_tc))
    if with_sounding:
        with open(os.path.join(d, "input_sounding"), "w") as f:
            f.write(sounding())
    return d


def run(tool, d, *args):
    env = dict(os.environ, OMP_NUM_THREADS="2", OMP_WAIT_POLICY="passive")
    return subprocess.run([sys.executable, tool, *args], cwd=d, env=env, capture_output=True, text=True)


def text(path):
    with open(path) as f:
        return f.read()


class Checker(unittest.TestCase):
    def test_supported(self):
        d = make_dir("check_ok")
        r = run(CHECK, d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("CPPWRF runs this namelist", r.stdout)

    def test_not_ported(self):
        d = make_dir("check_bad", mp=24)
        r = run(CHECK, d, "namelist.input")
        self.assertEqual(r.returncode, 1)
        self.assertIn("NOT PORTED  &physics.mp_physics", r.stdout)
        self.assertIn("suggestion: use one of 0, 1, 3, 4, 6, 8, 10, 28", r.stdout)

    def test_json_and_missing(self):
        d = make_dir("check_json")
        r = run(CHECK, d, "--json")
        self.assertEqual(r.returncode, 0)
        self.assertIn('"status": "ported"', r.stdout)
        self.assertEqual(run(CHECK, d, "nonexistent.input").returncode, 2)


@unittest.skipUnless(BIN and os.path.exists(BIN), "CPPWRF_BIN (a cppwrf binary) not set")
class DropIn(unittest.TestCase):
    def test_run_names_logs(self):
        d = make_dir("run")
        r = run(WRF, d, "--exe", BIN)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr + text(os.path.join(d, "rsl.error.0000")))
        self.assertIn("d01 0001-01-01_00:01:00 wrf: SUCCESS COMPLETE WRF", r.stdout)
        files = sorted(os.listdir(d))
        # frames_per_outfile = 2 with a frame every 24 s: two history files, as WRF names them
        for f in ("wrfout_d01_0001-01-01_00:00:00", "wrfout_d01_0001-01-01_00:00:48", "wrfrst_d01_0001-01-01_00:01:00",
                  "rsl.out.0000", "rsl.error.0000"):
            self.assertIn(f, files)
        out = text(os.path.join(d, "rsl.out.0000"))
        self.assertIn("Timing for main: time 0001-01-01_00:00:12 on domain   1:", out)
        self.assertIn("Timing for Writing wrfout_d01_0001-01-01_00:00:48 for domain        1:", out)
        self.assertIn("Timing for Writing restart for domain        1:", out)
        self.assertEqual(out.rstrip("\n").splitlines()[-1], "d01 0001-01-01_00:01:00 wrf: SUCCESS COMPLETE WRF")
        self.assertIn("SUCCESS COMPLETE WRF", text(os.path.join(d, "rsl.error.0000")))
        # CPPWRF's own files: the configuration and run.json in cppwrf_work, no restart_<step>.nc by default
        self.assertTrue(os.path.exists(os.path.join(d, "cppwrf_work", "run.ini")))
        for root, _, fs in os.walk(d):
            self.assertFalse([f for f in fs if f.startswith("restart_")], root)
        self.assertFalse([f for f in files if f.startswith("d0")])

        # the same run resumed from its wrfrst file at minute 1... of a 1-minute run: WRF's restart = .true. with the
        # start date of the file and one more minute
        d2 = make_dir("resumed", restart=".true.")
        nl = text(os.path.join(d2, "namelist.input")).replace("start_minute = 00", "start_minute = 01")
        with open(os.path.join(d2, "namelist.input"), "w") as f:
            f.write(nl)
        os.symlink(os.path.join(d, "wrfrst_d01_0001-01-01_00:01:00"), os.path.join(d2, "wrfrst_d01_0001-01-01_00:01:00"))
        r = run(WRF, d2, "--exe", BIN)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr + text(os.path.join(d2, "rsl.error.0000")))
        self.assertIn("d01 0001-01-01_00:02:00 wrf: SUCCESS COMPLETE WRF", r.stdout)
        self.assertIn("wrfout_d01_0001-01-01_00:01:12", os.listdir(d2))  # WRF: no frame at the restart time

    def test_restart_timers_refused(self):
        # WRF's override_restart_timers / write_hist_at_0h_rst change a resumed run's alarms: not ported, refused
        for k in ("override_restart_timers", "write_hist_at_0h_rst"):
            d = make_dir(f"rst_{k}", restart=".true.", extra_tc=f"{k} = .true.,")
            r = run(WRF, d, "--exe", BIN)
            self.assertEqual(r.returncode, 1, k)
            err = text(os.path.join(d, "rsl.error.0000"))
            self.assertIn("FATAL CALLED", err)
            self.assertIn(f"{k} = .true. with restart = .true. is not supported", err)

    def test_outnames(self):
        d = make_dir("outnames", extra_tc="history_outname = 'hist_d<domain>_<date>.nc', "
                                          "rst_outname = 'rst_d<domain>_<date>',")
        r = run(WRF, d, "--exe", BIN)
        self.assertEqual(r.returncode, 0, text(os.path.join(d, "rsl.error.0000")))
        files = os.listdir(d)
        self.assertIn("hist_d01_0001-01-01_00:00:00.nc", files)
        self.assertIn("rst_d01_0001-01-01_00:01:00", files)

    def test_not_ported_option(self):
        d = make_dir("bad_option", mp=24)
        r = run(WRF, d, "--exe", BIN)
        self.assertEqual(r.returncode, 1)
        err = text(os.path.join(d, "rsl.error.0000"))
        self.assertIn("FATAL CALLED", err)
        self.assertIn("mp_physics = 24 (WSM7", err)
        self.assertNotIn("SUCCESS", err + r.stdout)
        self.assertFalse([f for f in os.listdir(d) if f.startswith("wrfout")])

    def test_model_failure(self):
        d = make_dir("no_sounding", with_sounding=False)
        r = run(WRF, d, "--exe", BIN)
        self.assertNotEqual(r.returncode, 0)
        err = text(os.path.join(d, "rsl.error.0000"))
        self.assertIn("FATAL CALLED", err)
        self.assertIn("input_sounding", err)
        self.assertNotIn("SUCCESS COMPLETE WRF", err)

    def test_wrong_mode(self):
        d = make_dir("wrong_mode")
        other = "gpu-fp32" if "float64" in subprocess.run([BIN, "--version"], capture_output=True, text=True).stdout \
            else "cpu-fp64"
        r = run(WRF, d, "--exe", BIN, "--mode", other)
        self.assertEqual(r.returncode, 1)
        self.assertIn(f"mode {other}", text(os.path.join(d, "rsl.error.0000")))


@unittest.skipUnless(REAL_BIN and os.path.exists(REAL_BIN), "CPPWRF_REAL_BIN (cppwrf_real) not set")
class RealExe(unittest.TestCase):
    def test_no_met_em(self):
        d = make_dir("realexe")
        r = run(REALEXE, d, "--exe", REAL_BIN)
        self.assertEqual(r.returncode, 1)
        err = text(os.path.join(d, "rsl.error.0000"))
        self.assertIn("FATAL CALLED", err)
        self.assertIn("met_em.d01.0001-01-01_00:00:00.nc not found", err)
        self.assertNotIn("SUCCESS COMPLETE REAL_EM INIT", err + r.stdout)

    def test_cppwrf_real_failure(self):
        # a met_em file that is not NetCDF: cppwrf_real's own error, passed on with its exit status
        d = make_dir("realexe_bad")
        open(os.path.join(d, "met_em.d01.0001-01-01_00:00:00.nc"), "w").write("not netcdf\n")
        r = run(REALEXE, d, "--exe", REAL_BIN)
        self.assertNotEqual(r.returncode, 0)
        err = text(os.path.join(d, "rsl.error.0000"))
        self.assertIn("FATAL CALLED", err)
        self.assertIn("cppwrf_real stopped with exit status", err)


if __name__ == "__main__":
    unittest.main()
