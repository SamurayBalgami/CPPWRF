#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Verification tiers (data/option_tiers.txt; src/io/option_tiers.cpp, tools/wrf2cppwrf.py) with a dummy row that
puts Kessler (mp_physics = 1) below verified, through $CPPWRF_OPTION_TIERS_EXTRA:
- the translator warns and passes &cppwrf allow_unverified on as [cppwrf] allow_unverified;
- the model (tests/tools/test_dropin.py's tiny em_quarter_ss through cppwrf-wrf) stops without the switch, naming the
  option, the value and the tier; runs with it, the wrfout and wrfrst files carrying CPPWRF_UNVERIFIED_OPTIONS;
- without a row, the switch changes no byte of the wrfout files;
- cppwrf_real (cppwrf-realexe) stops on the row without &cppwrf allow_unverified.

  CPPWRF_BIN=<build>/cppwrf CPPWRF_REAL_BIN=<build>/src/init_real/cppwrf_real CPPWRF_TEST_DIR=<dir> \\
      python3 -m unittest tests/tools/test_option_tiers.py"""
import filecmp
import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "tools"))
import test_dropin as td  # noqa: E402
import wrf2cppwrf as w2c  # noqa: E402

DUMMY = "phys.mp_physics  physics.mp_physics  1  column_verified  dummy row of the tier test\n"
ALLOW = "&cppwrf\n allow_unverified = .true.,\n/\n"


def extra_file() -> str:
    os.makedirs(td.WORK, exist_ok=True)
    p = os.path.join(td.WORK, "tiers_dummy.txt")
    with open(p, "w") as f:
        f.write(DUMMY)
    return p


def make(name, allow=False):
    d = td.make_dir(name)
    if allow:
        with open(os.path.join(d, "namelist.input"), "a") as f:
            f.write(ALLOW)
    return d


def run(tool, d, rows, *args):
    env = dict(os.environ, OMP_NUM_THREADS="2", OMP_WAIT_POLICY="passive")
    env.pop("CPPWRF_OPTION_TIERS_EXTRA", None)
    if rows:
        env["CPPWRF_OPTION_TIERS_EXTRA"] = extra_file()
    return subprocess.run([sys.executable, tool, *args], cwd=d, env=env, capture_output=True, text=True)


def global_att(path, name):
    import netCDF4
    with netCDF4.Dataset(path) as nc:
        return nc.getncattr(name) if name in nc.ncattrs() else None


class Table(unittest.TestCase):
    def test_builtin_table(self):
        # the shipped table parses; each row (a port merged before its whole-model proof) names a WRF option and its
        # scheme and column driver, or why there is none
        with open(os.path.join(REPO, "data", "option_tiers.txt")) as f:
            for row in w2c.parse_option_tiers(f.read(), "data/option_tiers.txt"):
                self.assertTrue(row["wrf_key"] and row["note"], row)
                self.assertIn(row["tier"], w2c.TIERS)

    def test_parse_and_match(self):
        rows = w2c.parse_option_tiers("phys.mp_physics physics.mp_physics 4,20:23 column_verified WSM5 (x)\n"
                                      "- physics.sst_skin true unverified\n", "t")
        self.assertEqual(rows[0]["note"], "WSM5 (x)")
        self.assertIsNone(rows[1]["config_key"])
        self.assertTrue(w2c.tier_value_matches(rows[0]["values"], 21))
        self.assertTrue(w2c.tier_value_matches(rows[0]["values"], 23))
        self.assertTrue(w2c.tier_value_matches(rows[0]["values"], 20))
        self.assertFalse(w2c.tier_value_matches(rows[0]["values"], 19))
        self.assertTrue(w2c.tier_value_matches(rows[0]["values"], 4.0))
        self.assertFalse(w2c.tier_value_matches(rows[0]["values"], 5))
        self.assertTrue(w2c.tier_value_matches(rows[1]["values"], True))
        self.assertFalse(w2c.tier_value_matches(rows[1]["values"], False))
        for bad in ("a b c", "a physics.x 1 verified", "a physicsx 1 unverified", "a physics.x 3:2 unverified"):
            with self.assertRaises(ValueError):
                w2c.parse_option_tiers(bad, "t")


class Translator(unittest.TestCase):
    def translate(self, text, rows):
        saved = w2c.OPTION_TIERS
        w2c.OPTION_TIERS = w2c.parse_option_tiers(rows, "dummy")
        try:
            tr = w2c.WrfTranslator(w2c.parse_namelist(text), namelist_path="/x/qss/namelist.input")
            sections = tr.translate()
            return tr, sections, tr.format_ini(sections)
        finally:
            w2c.OPTION_TIERS = saved

    def test_warns_without_switch(self):
        tr, sections, ini = self.translate(td.NAMELIST.format(restart=".false.", mp=1, extra_tc=""), DUMMY)
        self.assertEqual(tr.errors, [])
        w = [x for x in tr.warnings if "verification tier" in x]
        self.assertEqual(len(w), 1, tr.warnings)
        self.assertIn("&physics.mp_physics = 1: verification tier column_verified", w[0])
        self.assertIn("stops on it unless &cppwrf allow_unverified = .true.", w[0])
        self.assertNotIn("cppwrf", sections)
        self.assertNotIn("[cppwrf]", ini)

    def test_switch_passed_through(self):
        text = td.NAMELIST.format(restart=".false.", mp=1, extra_tc="") + ALLOW
        tr, sections, ini = self.translate(text, DUMMY)
        self.assertEqual(tr.errors, [])
        self.assertIn("runs because &cppwrf allow_unverified = .true.", " ".join(tr.warnings))
        self.assertIn("[cppwrf]\nallow_unverified = true\n", ini)
        # no row: no warning, the switch still passed on; another value of the option: no warning
        tr, sections, ini = self.translate(text, "")
        self.assertFalse([x for x in tr.warnings if "verification tier" in x or "&cppwrf" in x], tr.warnings)
        self.assertEqual(sections["cppwrf"], {"allow_unverified": True})
        tr, _, _ = self.translate(td.NAMELIST.format(restart=".false.", mp=6, extra_tc=""), DUMMY)
        self.assertFalse([x for x in tr.warnings if "verification tier" in x])

    def test_bad_cppwrf_group(self):
        base = td.NAMELIST.format(restart=".false.", mp=1, extra_tc="")
        tr, _, _ = self.translate(base + "&cppwrf\n allow_unverified = 3,\n/\n", "")
        self.assertTrue([e for e in tr.errors if "allow_unverified = 3" in e], tr.errors)
        tr, _, _ = self.translate(base + "&cppwrf\n allow_unverifed = .true.,\n/\n", "")
        self.assertTrue([e for e in tr.errors if "&cppwrf.allow_unverifed is not a CPPWRF option" in e], tr.errors)

    def test_checker_reports_ported(self):
        # cppwrf-check: the option is ported (below verified), not "not ported" or "ignored"
        r = run(td.CHECK, make("tier_check"), True, "--verbose")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertNotIn("NOT PORTED", r.stdout)
        self.assertIn("verification tier column_verified", r.stdout)
        self.assertRegex(r.stdout, r"PORTED +&physics\.mp_physics +=")
        self.assertIn(" 0 ignored, 0 not ported", r.stdout)


@unittest.skipUnless(td.BIN and os.path.exists(td.BIN), "CPPWRF_BIN (a cppwrf binary) not set")
class Model(unittest.TestCase):
    def test_rejected_then_allowed(self):
        d = make("tier_rejected")
        r = run(td.WRF, d, True, "--exe", td.BIN)
        err = td.text(os.path.join(d, "rsl.error.0000"))
        self.assertNotEqual(r.returncode, 0, r.stdout + err)
        self.assertIn("phys.mp_physics = 1 is column_verified (dummy row of the tier test)", err)
        self.assertIn("set cppwrf.allow_unverified = true", err)
        self.assertFalse([f for f in os.listdir(d) if f.startswith("wrfout")])

        d = make("tier_allowed", allow=True)
        r = run(td.WRF, d, True, "--exe", td.BIN)
        err = td.text(os.path.join(d, "rsl.error.0000"))
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr + err)
        self.assertIn("cppwrf: warning: phys.mp_physics = 1 is column_verified", err)
        for f in ("wrfout_d01_0001-01-01_00:00:00", "wrfout_d01_0001-01-01_00:00:48", "wrfrst_d01_0001-01-01_00:01:00"):
            self.assertEqual(global_att(os.path.join(d, f), "CPPWRF_UNVERIFIED_OPTIONS"),
                             "mp_physics=1:column_verified", f)

    def test_switch_alone_changes_no_byte(self):
        plain, allowed = make("tier_plain"), make("tier_switch_only", allow=True)
        for d in (plain, allowed):
            r = run(td.WRF, d, False, "--exe", td.BIN)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr + td.text(os.path.join(d, "rsl.error.0000")))
        outs = sorted(f for f in os.listdir(plain) if f.startswith(("wrfout", "wrfrst")))
        self.assertEqual(len(outs), 3)
        for f in outs:
            self.assertIsNone(global_att(os.path.join(allowed, f), "CPPWRF_UNVERIFIED_OPTIONS"))
            self.assertTrue(filecmp.cmp(os.path.join(plain, f), os.path.join(allowed, f), shallow=False), f)


@unittest.skipUnless(td.REAL_BIN and os.path.exists(td.REAL_BIN), "CPPWRF_REAL_BIN (cppwrf_real) not set")
class RealExe(unittest.TestCase):
    def test_gate(self):
        # a met_em file that is not NetCDF: the gate comes first, with the switch cppwrf_real goes on and fails on it
        for name, allow in (("tier_realexe", False), ("tier_realexe_allowed", True)):
            d = make(name, allow=allow)
            with open(os.path.join(d, "met_em.d01.0001-01-01_00:00:00.nc"), "w") as f:
                f.write("not netcdf\n")
            r = run(td.REALEXE, d, True, "--exe", td.REAL_BIN)
            self.assertNotEqual(r.returncode, 0)
            out = r.stdout + r.stderr + td.text(os.path.join(d, "rsl.error.0000"))
            if allow:
                self.assertNotIn("verified tier", out)
                self.assertIn("warning: real: d01: &physics.mp_physics = 1 is column_verified", out)
            else:
                self.assertIn("real: d01: options below the verified tier", out)
                self.assertIn("&physics.mp_physics = 1 is column_verified", out)


if __name__ == "__main__":
    unittest.main()
