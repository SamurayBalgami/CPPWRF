#!/usr/bin/env python3
"""Unit tests for tools/wrf2cppwrf.py namelist translator.

Tests:
1. Pure Python Fortran namelist tokenizer and parser.
2. Translation and key-by-key comparison against CPPWRF reference configs.
3. Multiple domains detection and warning emission.
4. Error detection on unsupported physics/dynamics options.
5. 2-step model execution under tools/heavy.sh.


"""

import os
import subprocess
import sys
import tempfile
import unittest

from typing import Any, Dict, List, Optional

# Add workspace to path
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO_ROOT)

from tools.wrf2cppwrf import parse_namelist, WrfTranslator, NamelistParseError


def load_ini_dict(path_or_text: str, is_text: bool = False) -> dict:
    """Parse INI file or text into a dict of {section: {key: value}}."""
    sections = {}
    cur_sec = None
    if is_text:
        lines = path_or_text.splitlines()
    else:
        with open(path_or_text) as f:
            lines = f.readlines()
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            cur_sec = line[1:-1]
            sections[cur_sec] = {}
        elif "=" in line and cur_sec:
            k, v = line.split("=", 1)
            sections[cur_sec][k.strip()] = v.strip()
    return sections


class TestNamelistParser(unittest.TestCase):
    """Test suite for the f90nml-free Fortran namelist parser."""

    def test_basic_parsing(self):
        text = """
        &time_control
         run_days = 0,
         run_hours = 1,
         restart = .false.,
         debug_level = 0
        /
        &domains
         time_step = 60,
         dx = 15000,
         dy = 15000.0,
        /
        """
        nml = parse_namelist(text)
        self.assertIn("time_control", nml)
        self.assertIn("domains", nml)
        self.assertEqual(nml["time_control"]["run_days"], [0])
        self.assertEqual(nml["time_control"]["run_hours"], [1])
        self.assertEqual(nml["time_control"]["restart"], [False])
        self.assertEqual(nml["domains"]["time_step"], [60])
        self.assertEqual(nml["domains"]["dx"], [15000])
        self.assertEqual(nml["domains"]["dy"], [15000.0])

    def test_comments_and_whitespace(self):
        text = """
        ! Top-level comment
        &time_control ! inline comment
         run_days = 2, ! days
         ! full line comment in group
         run_hours = 12
        /
        """
        nml = parse_namelist(text)
        self.assertEqual(nml["time_control"]["run_days"], [2])
        self.assertEqual(nml["time_control"]["run_hours"], [12])

    def test_string_literals(self):
        text = '''
        &io
         filename = 'my_file.nc',
         msg = "Hello, world!",
         nested_quotes = 'It''s fine'
        /
        '''
        nml = parse_namelist(text)
        self.assertEqual(nml["io"]["filename"], ["my_file.nc"])
        self.assertEqual(nml["io"]["msg"], ["Hello, world!"])
        self.assertEqual(nml["io"]["nested_quotes"], ["It's fine"])

    def test_booleans(self):
        text = """
        &flags
         b1 = .true.,
         b2 = .FALSE.,
         b3 = .t.,
         b4 = .f.,
         b5 = T,
         b6 = F
        /
        """
        nml = parse_namelist(text)
        grp = nml["flags"]
        self.assertEqual(grp["b1"], [True])
        self.assertEqual(grp["b2"], [False])
        self.assertEqual(grp["b3"], [True])
        self.assertEqual(grp["b4"], [False])
        self.assertEqual(grp["b5"], [True])
        self.assertEqual(grp["b6"], [False])

    def test_repetition_syntax(self):
        text = """
        &dynamics
         diff_opt = 3*0,
         factors = 2*0.12, 1.0
        /
        """
        nml = parse_namelist(text)
        self.assertEqual(nml["dynamics"]["diff_opt"], [0, 0, 0])
        self.assertEqual(nml["dynamics"]["factors"], [0.12, 0.12, 1.0])

    def test_scientific_notation_and_fortran_d(self):
        text = """
        &floats
         a = 1.0d-4,
         b = 2.5D+02,
         c = 5000.,
         d = .12,
         e = -1.5e-3
        /
        """
        nml = parse_namelist(text)
        grp = nml["floats"]
        self.assertAlmostEqual(grp["a"][0], 1.0e-4)
        self.assertAlmostEqual(grp["b"][0], 250.0)
        self.assertAlmostEqual(grp["c"][0], 5000.0)
        self.assertAlmostEqual(grp["d"][0], 0.12)
        self.assertAlmostEqual(grp["e"][0], -1.5e-3)

    def test_multiple_domains_warning(self):
        text = """
        &domains
         max_dom = 2,
         e_we = 120, 240,
         dx = 15000, 5000
        /
        """
        nml = parse_namelist(text)
        translator = WrfTranslator(nml)
        translator.translate()
        warnings_joined = " ".join(translator.warnings)
        self.assertIn("multiple domains", warnings_joined.lower())


class TestRealCasesReproduction(unittest.TestCase):
    """Test suite comparing generated INI against configs/real_*.ini key by key."""

    def _compare_configs(self, nml_rel_path: str, ref_ini_rel_path: str,
                          ignore_keys: Optional[dict] = None):
        nml_path = os.path.expanduser(nml_rel_path)
        ref_path = os.path.join(REPO_ROOT, ref_ini_rel_path)
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")
        if not os.path.exists(ref_path):
            self.skipTest(f"Reference INI not found: {ref_path}")

        with open(nml_path) as f:
            nml = parse_namelist(f.read())
        translator = WrfTranslator(nml, namelist_path=nml_path)
        gen = translator.translate()
        ref = load_ini_dict(ref_path)

        all_secs = sorted(set(list(gen.keys()) + list(ref.keys())))
        for s in all_secs:
            self.assertIn(s, gen, f"[{s}] missing in generated config")
            self.assertIn(s, ref, f"[{s}] extra in generated config")

            all_keys = sorted(set(list(gen[s].keys()) + list(ref[s].keys())))
            for k in all_keys:
                if ignore_keys and s in ignore_keys and k in ignore_keys[s]:
                    continue
                self.assertIn(k, gen[s], f"[{s}].{k} missing in generated config")
                self.assertIn(k, ref[s], f"[{s}].{k} extra in generated config")

                # Compare normalized string values
                gv = str(gen[s][k]).replace(" ", "")
                rv = str(ref[s][k]).replace(" ", "")
                self.assertEqual(gv, rv, f"[{s}].{k} mismatch: generated '{gen[s][k]}' vs ref '{ref[s][k]}'")

    def test_real_r0_key_by_key(self):
        self._compare_configs("~/wrf_real/run_R0_r8/namelist.input", "configs/real_r0.ini")

    def test_real_r1_key_by_key(self):
        self._compare_configs("~/wrf_real/run_R1_r8/namelist.input", "configs/real_r1.ini")

    def test_real_r2_key_by_key(self):
        self._compare_configs("~/wrf_real/run_R2/namelist.input", "configs/real_r2.ini")

    def test_real_r3_key_by_key(self):
        self._compare_configs("~/wrf_real/run_R3/namelist.input", "configs/real_r3.ini")

    def test_real_mp6_key_by_key(self):
        self._compare_configs("~/wrf_real/run_mp6/namelist.input", "configs/real_mp6.ini")

    def test_real_p1_key_by_key(self):
        self._compare_configs("~/wrf_real/run_P1/namelist.input", "configs/real_p1.ini")

    def test_real_p2_key_by_key(self):
        self._compare_configs("~/wrf_real/run_P2/namelist.input", "configs/real_p2.ini")

    def test_real_p3_key_by_key(self):
        self._compare_configs("~/wrf_real/run_P3/namelist.input", "configs/real_p3.ini")

    def test_real_p4_diff_opt_mapped(self):
        # run_P4 namelist specifies diff_opt = 1, km_opt = 4
        nml_path = os.path.expanduser("~/wrf_real/run_P4/namelist.input")
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")
        with open(nml_path) as f:
            nml = parse_namelist(f.read())
        t = WrfTranslator(nml, namelist_path=nml_path)
        gen = t.translate()
        self.assertEqual(gen["dyn"].get("diff_opt"), 1)
        self.assertEqual(gen["dyn"].get("km_opt"), 4)
        self.assertEqual(gen["phys"].get("mp_physics"), 6)


class TestIdealCases(unittest.TestCase):
    """Test suite for ideal cases translation."""

    def test_hill2d_mapping(self):
        nml_path = os.path.expanduser("~/wrf_ref/run_hill_p0/namelist.input")
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")
        with open(nml_path) as f:
            nml = parse_namelist(f.read())
        t = WrfTranslator(nml, namelist_path=nml_path)
        gen = t.translate()
        self.assertEqual(gen["run"]["case"], "hill2d")
        self.assertEqual(gen["domain"]["nx"], 201)
        self.assertEqual(gen["domain"]["ny"], 1)
        self.assertEqual(gen["domain"]["nz"], 40)
        self.assertEqual(gen["domain"]["dx"], 2000)
        self.assertEqual(gen["time"]["dt"], 20)
        self.assertEqual(gen["time"]["time_step_sound"], 6)

    def test_quarter_ss_mapping(self):
        nml_path = os.path.expanduser("~/wrf_ref/run_qss/namelist.input")
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")
        with open(nml_path) as f:
            nml = parse_namelist(f.read())
        t = WrfTranslator(nml, namelist_path=nml_path)
        gen = t.translate()
        self.assertEqual(gen["run"]["case"], "quarter_ss")
        self.assertEqual(gen["domain"]["nx"], 41)
        self.assertEqual(gen["domain"]["ny"], 41)
        self.assertEqual(gen["domain"]["nz"], 40)
        self.assertEqual(gen["time"]["dt"], 12)
        self.assertEqual(gen["time"]["time_step_sound"], 6)

    def test_straka_mapping(self):
        nml_path = os.path.expanduser("~/wrf_ref/run_straka_d1/namelist.input")
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")
        with open(nml_path) as f:
            nml = parse_namelist(f.read())
        t = WrfTranslator(nml, namelist_path=nml_path)
        gen = t.translate()
        self.assertEqual(gen["run"]["case"], "straka")
        self.assertEqual(gen["domain"]["nx"], 512)
        self.assertEqual(gen["domain"]["ny"], 1)
        self.assertEqual(gen["domain"]["nz"], 64)
        self.assertEqual(gen["time"]["dt"], 1)
        self.assertEqual(gen["time"]["time_step_sound"], 6)


class TestValidationAndWarnings(unittest.TestCase):
    """Test error and warning detection on unsupported/deviating options."""

    def test_unsupported_mp_physics(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &physics mp_physics = 2 /
        &dynamics /
        """
        nml = parse_namelist(text)
        t = WrfTranslator(nml)
        t.translate()
        self.assertTrue(any("mp_physics = 2" in w for w in t.warnings))

    def test_unsupported_diff_opt(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &dynamics diff_opt = 2 /
        """
        nml = parse_namelist(text)
        t = WrfTranslator(nml)
        t.translate()
        self.assertTrue(any("diff_opt = 2" in w for w in t.warnings))

    def test_hydrostatic_error(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &dynamics non_hydrostatic = .false. /
        """
        nml = parse_namelist(text)
        t = WrfTranslator(nml)
        t.translate()
        self.assertTrue(any("non_hydrostatic = .false." in e for e in t.errors))


class TestExecution2Steps(unittest.TestCase):
    """Verify that generated INI executes for 2 steps with build-a0-cpu/cppwrf."""

    def test_run_r0_2steps(self):
        exe = os.path.join(REPO_ROOT, "build-a0-cpu", "cppwrf")
        if not os.path.exists(exe):
            self.skipTest("build-a0-cpu/cppwrf binary not found")

        nml_path = os.path.expanduser("~/wrf_real/run_R0_r8/namelist.input")
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")

        with tempfile.NamedTemporaryFile(suffix=".ini", mode="w", delete=False) as tf:
            ini_path = tf.name

        try:
            # Generate INI using translator CLI
            cmd_trans = [sys.executable, os.path.join(REPO_ROOT, "tools", "wrf2cppwrf.py"),
                         nml_path, "-o", ini_path]
            ret_trans = subprocess.run(cmd_trans, capture_output=True, text=True)
            self.assertEqual(ret_trans.returncode, 0, f"Translator failed: {ret_trans.stderr}")

            # Run 2 steps with heavy.sh
            cmd_run = [os.path.join(REPO_ROOT, "tools", "heavy.sh"), exe, ini_path,
                       "--set", "run.run_time=120"]
            ret_run = subprocess.run(cmd_run, capture_output=True, text=True, cwd=REPO_ROOT)
            self.assertEqual(ret_run.returncode, 0, f"Run failed: {ret_run.stderr}\n{ret_run.stdout}")
            self.assertIn("done, 2 steps", ret_run.stdout)
        finally:
            if os.path.exists(ini_path):
                os.remove(ini_path)


if __name__ == "__main__":
    unittest.main()
