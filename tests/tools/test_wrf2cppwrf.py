#!/usr/bin/env python3
"""Unit tests for tools/wrf2cppwrf.py namelist translator.

Tests:
1. Pure Python Fortran namelist tokenizer and parser.
2. Translation and key-by-key comparison against CPPWRF reference configs.
3. Multiple domains detection and warning emission.
4. Error detection on unsupported physics/dynamics options.
5. 2-step model execution under tools/heavy.sh.

The tests of WRF reference runs read ${CPPWRF_WRF_RUNS} and ${CPPWRF_WRF_REF} (defaults as tests/verify/wrf_suite.py:
~/wrf_real480, ~/wrf_ref480) and skip when the files are missing; the model runs also need $CPPWRF_BIN.
"""

import os
import subprocess
import sys
import tempfile
import unittest

# set in the environment: the configurations and the model read ${CPPWRF_WRF_RUNS}
os.environ.setdefault("CPPWRF_WRF_RUNS", os.path.expanduser("~/wrf_real480"))
os.environ.setdefault("CPPWRF_WRF_REF", os.path.expanduser("~/wrf_ref480"))


def wrf_runs_path(rel: str) -> str:
    """Path of a WRF reference run file under ${CPPWRF_WRF_RUNS}."""
    return os.path.join(os.environ["CPPWRF_WRF_RUNS"], rel)


def wrf_ref_path(rel: str) -> str:
    """Path of a WRF ideal-case run file under ${CPPWRF_WRF_REF}."""
    return os.path.join(os.environ["CPPWRF_WRF_REF"], rel)

from typing import Any, Dict, List, Optional

# Add workspace to path
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO_ROOT)

from tools.wrf2cppwrf import parse_namelist, WrfTranslator, NamelistParseError


def real_case(text: str) -> str:
    """A real-case namelist: domain 1 needs specified = .true. (Registry default .false.; with no condition on a
    side wrf.exe stops in share/module_bc.F boundary_condition_check, and the translator asks for --case)."""
    return text + "\n&bdy_control specified = .true. /\n"


# real-case keys of the translation that wrf.exe checks against the input file: (namelist group, key, Registry
# default; None: written only when the namelist has it)
REAL_INPUT_CHECKS = {
    ("real", "p_top_requested"): ("domains", "p_top_requested", 5000.0),
    ("real", "iso_temp"): ("dynamics", "iso_temp", 200.0),
    ("real", "dx"): ("domains", "dx", None),
    ("real", "dy"): ("domains", "dy", None),
    ("dyn", "hybrid_opt"): ("dynamics", "hybrid_opt", 2),
}


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

    def test_multiple_domains_translation(self):
        text = """
        &time_control
         run_hours = 1,
         history_interval = 10, 10, 10,
         input_from_file = .true., .true., .true.
        /
        &domains
         time_step = 60,
         max_dom = 3,
         e_we = 120, 121, 91,
         e_sn = 90, 91, 61,
         e_vert = 45, 45, 45,
         dx = 15000, 5000, 1000,
         dy = 15000, 5000, 1000,
         parent_id = 0, 1, 2,
         i_parent_start = 1, 40, 30,
         j_parent_start = 1, 30, 20,
         parent_grid_ratio = 1, 3, 5,
         parent_time_step_ratio = 1, 3, 5,
         feedback = 1,
         smooth_option = 2
        /
        &dynamics
         khdif = 0, 0, 10,
         kvdif = 0, 0, 5
        /
        &bdy_control
         specified = .true., .false., .false.,
         nested = .false., .true., .true.
        /
        """
        nml = parse_namelist(text)
        translator = WrfTranslator(nml)
        gen = translator.translate()
        self.assertEqual(translator.errors, [])
        self.assertEqual(translator.warnings, [])
        self.assertEqual(gen["domains"]["max_dom"], 3)
        self.assertEqual(gen["domains"]["feedback"], 1)
        self.assertEqual(gen["domains"]["smooth_option"], 2)
        self.assertEqual(gen["domains"]["i_parent_start"], 40)
        self.assertEqual(gen["domains"]["j_parent_start"], 30)
        self.assertEqual(gen["domains"]["parent_grid_ratio"], 3)
        self.assertEqual(gen["domains"]["parent_time_step_ratio"], 3)
        self.assertEqual(gen["d02"]["wrfinput"], "wrfinput_d02")
        # no [d0N] dt: the model takes the parent's dt / parent_time_step_ratio (20 s, 4 s) as WRF does and
        # refuses another value (src/cases/nest_driver.cpp); the history intervals below are in those steps
        self.assertNotIn("dt", gen["d02"])
        self.assertEqual(gen["d02"]["out_interval"], 30)
        self.assertEqual(gen["d02"]["diag_interval"], 30)
        self.assertEqual(gen["d03"]["parent_id"], 2)
        self.assertEqual(gen["d03"]["i_parent_start"], 30)
        self.assertEqual(gen["d03"]["j_parent_start"], 20)
        self.assertEqual(gen["d03"]["parent_grid_ratio"], 5)
        self.assertEqual(gen["d03"]["parent_time_step_ratio"], 5)
        self.assertEqual(gen["d03"]["wrfinput"], "wrfinput_d03")
        self.assertNotIn("dt", gen["d03"])
        self.assertEqual(gen["d03"]["out_interval"], 150)
        self.assertEqual(gen["d03"]["diag_interval"], 150)
        self.assertEqual(gen["d03"]["khdif"], 10)
        self.assertEqual(gen["d03"]["kvdif"], 5)


class TestRealCasesReproduction(unittest.TestCase):
    """Test suite comparing generated INI against configs/real_*.ini key by key."""

    def _compare_configs(self, nml_rel_path: str, ref_ini_rel_path: str,
                          ignore_keys: Optional[dict] = None, sections: Optional[tuple] = None):
        nml_path = wrf_runs_path(nml_rel_path)
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
        self.assertEqual(translator.errors, [])

        ignore_keys = {s: set(keys) for s, keys in (ignore_keys or {}).items()}
        if gen["run"]["case"] == "real":
            # namelist values the model checks against wrfinput_d01 as wrf.exe does (start_em.F, input_wrf.F stop
            # when they differ): the namelist's value or the Registry default; the hand-written configs leave them
            # to the input file
            for (sec, key), (grp, nkey, default) in REAL_INPUT_CHECKS.items():
                nvals = nml.get(grp, {}).get(nkey)
                if not nvals and default is None:
                    self.assertNotIn(key, gen.get(sec, {}), f"[{sec}].{key}")
                    continue
                self.assertEqual(float(gen[sec][key]), float(nvals[0] if nvals else default), f"[{sec}].{key}")
                if key not in ref.get(sec, {}):
                    ignore_keys.setdefault(sec, set()).add(key)
        # WRF writes wrfrst files every restart_interval minutes whatever restart is; a hand-written config of a run
        # that ends before the first restart time has no [restart]
        rst = nml.get("time_control", {}).get("restart_interval", [0])[0]
        if rst > 0 and "restart" not in ref:
            self.assertLess(gen["run"]["run_time"], rst * 60)
            self.assertEqual(gen.pop("restart"), {"interval": rst * 60})

        all_secs = sorted(set(list(gen.keys()) + list(ref.keys())))
        if sections:  # only these sections (the hand-written configs add run/output settings)
            all_secs = [x for x in all_secs if x in sections]
        for s in all_secs:
            self.assertIn(s, gen, f"[{s}] missing in generated config")
            self.assertIn(s, ref, f"[{s}] extra in generated config")

            all_keys = sorted(set(list(gen[s].keys()) + list(ref[s].keys())))
            for k in all_keys:
                if ignore_keys and s in ignore_keys and k in ignore_keys[s]:
                    continue
                if sections and (k not in ref[s] or k not in gen[s]):
                    continue  # default values written by one side only
                self.assertIn(k, gen[s], f"[{s}].{k} missing in generated config")
                self.assertIn(k, ref[s], f"[{s}].{k} extra in generated config")

                # Compare normalized string values (treating '20' and '20.0' as equal)
                gv = str(gen[s][k]).replace(" ", "")
                rv = str(ref[s][k]).replace(" ", "")
                if gv != rv and (gv.endswith(".0") or rv.endswith(".0")):
                    try:
                        if float(gv) == float(rv):
                            gv = rv
                    except ValueError:
                        pass
                self.assertEqual(gv, rv, f"[{s}].{k} mismatch: generated '{gen[s][k]}' vs ref '{ref[s][k]}'")

    def test_real_r0_key_by_key(self):
        self._compare_configs("run_R0_r8/namelist.input", "configs/real_r0.ini")

    def test_real_r1_key_by_key(self):
        self._compare_configs("run_R1_r8/namelist.input", "configs/real_r1.ini")

    def test_real_r2_key_by_key(self):
        self._compare_configs("run_R2/namelist.input", "configs/real_r2.ini")

    def test_real_r3_key_by_key(self):
        self._compare_configs("run_R3/namelist.input", "configs/real_r3.ini")

    def test_real_mp6_key_by_key(self):
        self._compare_configs("run_mp6/namelist.input", "configs/real_mp6.ini")

    def test_real_p1_key_by_key(self):
        self._compare_configs("run_P1_r8/steps/namelist.input", "configs/real_p1_r8.ini", sections=("phys", "dyn"))

    def test_real_p2_key_by_key(self):
        self._compare_configs("run_P2_r8/steps/namelist.input", "configs/real_p2_r8.ini", sections=("phys", "dyn"))

    def test_real_p3_key_by_key(self):
        self._compare_configs("run_P3_r8/steps/namelist.input", "configs/real_p3_r8.ini", sections=("phys", "dyn"))

    def test_real_p8_key_by_key(self):
        self._compare_configs("run_P8_r8/steps/namelist.input", "configs/real_p8_r8.ini", sections=("phys", "dyn"))

    def test_real_n1_key_by_key(self):
        self._compare_configs("run_N1/namelist.input", "configs/real_n1.ini", sections=("domains", "d02"))

    def test_real_n2_key_by_key(self):
        self._compare_configs("run_N2/namelist.input", "configs/real_n2.ini", sections=("domains", "d02"))

    def test_real_n1_steps_key_by_key(self):
        self._compare_configs("run_N1_r8/steps/namelist.input", "configs/real_n1_steps.ini", sections=("domains", "d02"))

    def test_real_n2_steps_key_by_key(self):
        self._compare_configs("run_N2_r8/steps/namelist.input", "configs/real_n2_steps.ini", sections=("domains", "d02"))

    def test_real_p4_diff_opt_mapped(self):
        # run_P4 namelist specifies diff_opt = 1, km_opt = 4
        nml_path = wrf_runs_path("run_P4/namelist.input")
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
        nml_path = wrf_ref_path("run_hill_p0/namelist.input")
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")
        with open(nml_path) as f:
            nml = parse_namelist(f.read())
        t = WrfTranslator(nml, namelist_path=nml_path)
        gen = t.translate()
        self.assertEqual(gen["run"]["case"], "em_hill2d_x")
        self.assertEqual(gen["domain"]["nx"], 201)
        self.assertEqual(gen["domain"]["ny"], 1)
        self.assertEqual(gen["domain"]["nz"], 40)
        self.assertEqual(gen["domain"]["dx"], 2000)
        self.assertEqual(gen["time"]["dt"], 20)
        self.assertEqual(gen["time"]["time_step_sound"], 6)

    def test_quarter_ss_mapping(self):
        nml_path = wrf_ref_path("run_qss/namelist.input")
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")
        with open(nml_path) as f:
            nml = parse_namelist(f.read())
        t = WrfTranslator(nml, namelist_path=nml_path)
        gen = t.translate()
        self.assertEqual(gen["run"]["case"], "em_quarter_ss")
        self.assertEqual(gen["domain"]["nx"], 41)
        self.assertEqual(gen["domain"]["ny"], 41)
        self.assertEqual(gen["domain"]["nz"], 40)
        self.assertEqual(gen["time"]["dt"], 12)
        self.assertEqual(gen["time"]["time_step_sound"], 6)

    def test_straka_mapping(self):
        nml_path = wrf_ref_path("run_straka_d1/namelist.input")
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")
        with open(nml_path) as f:
            nml = parse_namelist(f.read())
        t = WrfTranslator(nml, namelist_path=nml_path)
        gen = t.translate()
        self.assertEqual(gen["run"]["case"], "em_grav2d_x")
        self.assertEqual(gen["domain"]["nx"], 512)
        self.assertEqual(gen["domain"]["ny"], 1)
        self.assertEqual(gen["domain"]["nz"], 64)
        self.assertEqual(gen["time"]["dt"], 1)
        self.assertEqual(gen["time"]["time_step_sound"], 6)

    def test_wrf_ideal_case_names_and_keys(self):
        # Every WRF ideal_case that CPPWRF ports maps to the src/cases/wrf_ideal.cpp case name; the
        # sounding is the input_sounding next to the namelist; hybrid_opt = 0 is kept (the model's
        # default is 2); open_xs/open_xe give bc_x = open.
        tmp = tempfile.mkdtemp()
        nml_path = os.path.join(tmp, "namelist.input")
        expected = {1: "em_hill2d_x", 2: "em_quarter_ss", 4: "em_squall2d_x", 5: "em_squall2d_y", 6: "em_grav2d_x"}
        for ideal_case, name in expected.items():
            text = f"""
            &time_control run_minutes = 10 /
            &domains time_step = 3, e_we = 202, e_sn = 3, e_vert = 81, dx = 250, dy = 250, ztop = 20000 /
            &dynamics hybrid_opt = 0 /
            &bdy_control periodic_x = .false., open_xs = .true., open_xe = .true., periodic_y = .true. /
            &ideal ideal_case = {ideal_case} /
            """
            t = WrfTranslator(parse_namelist(text), namelist_path=nml_path)
            gen = t.translate()
            self.assertEqual(t.errors, [], msg=f"ideal_case = {ideal_case}")
            self.assertEqual(gen["run"]["case"], name)
            self.assertEqual(gen["case"]["sounding_file"], os.path.join(tmp, "input_sounding"))
            self.assertEqual(gen["dyn"]["hybrid_opt"], 0)
            self.assertEqual(gen["domain"]["bc_x"], "open")
            self.assertEqual(gen["domain"]["dy"], 250)
            # no damp_opt in the namelist: WRF's Registry default 3 (the model's own default is 0)
            self.assertEqual(gen["dyn"]["damp_opt"], 3)


class TestValidationAndWarnings(unittest.TestCase):
    """Test error and warning detection on unsupported/deviating options."""

    def test_unsupported_mp_physics(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &physics mp_physics = 2 /
        &dynamics /
        """
        nml = parse_namelist(real_case(text))
        t = WrfTranslator(nml)
        t.translate()
        self.assertTrue(any("mp_physics = 2" in e for e in t.errors), t.errors)

    def test_unsupported_diff_opt(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &dynamics diff_opt = 3 /
        """
        nml = parse_namelist(real_case(text))
        t = WrfTranslator(nml)
        t.translate()
        self.assertTrue(any("diff_opt = 3" in e for e in t.errors), t.errors)

    def test_hydrostatic_error(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &dynamics non_hydrostatic = .false. /
        """
        nml = parse_namelist(real_case(text))
        t = WrfTranslator(nml)
        t.translate()
        self.assertTrue(any("non_hydrostatic = .false." in e for e in t.errors))

    def test_noahmp_mapping(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &physics sf_sfclay_physics = 1, sf_surface_physics = 4, bl_pbl_physics = 1 /
        &noah_mp dveg = 4, opt_runsrf = 3, opt_runsub = 3 /
        &dynamics /
        """
        nml = parse_namelist(real_case(text))
        t = WrfTranslator(nml)
        gen = t.translate()
        self.assertEqual(t.warnings, [])
        self.assertEqual(t.errors, [])
        self.assertEqual(gen["phys"]["sf_surface_physics"], 4)
        self.assertEqual(gen["phys"]["num_soil_layers"], 4)
        self.assertEqual(gen["phys"]["dveg"], 4)
        self.assertEqual(gen["phys"]["opt_runsrf"], 3)
        self.assertEqual(gen["phys"]["opt_runsub"], 3)

    def test_extended_physics_options_mapping(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &physics
         mp_physics = 6,
         mp_tend_lim = 5.0,
         hail_opt = 1,
         use_mp_re = 0,
         sf_sfclay_physics = 1,
         sf_surface_physics = 2,
         num_soil_layers = 4,
         bl_pbl_physics = 1,
         ysu_topdown_pblmix = 0,
         cu_physics = 3,
         ishallow = 1,
         clos_choice = 2,
         kfeta_trigger = 3,
         convtrans_avglen_m = 45.0,
         ra_lw_physics = 4,
         ra_sw_physics = 4,
         radt = 10,
         icloud = 0,
         cldovrlp = 4,
         idcor = 1,
         o3input = 0,
         ghg_input = 0,
         aer_opt = 0,
         rdmaxalb = .false.,
         opt_thcnd = 2,
         kim_tofd = .true.,
         tofd_factor = 0.005
        /
        &dynamics /
        """
        nml = parse_namelist(real_case(text))
        t = WrfTranslator(nml)
        gen = t.translate()
        self.assertEqual(t.errors, [])
        p = gen["phys"]
        self.assertEqual(p["mp_tend_lim"], 5.0)
        self.assertEqual(p["hail_opt"], 1)
        self.assertEqual(p["use_mp_re"], 0)
        self.assertEqual(p["num_soil_layers"], 4)
        self.assertEqual(p["ysu_topdown_pblmix"], 0)
        self.assertEqual(p["ishallow"], 1)
        self.assertEqual(p["clos_choice"], 2)
        self.assertEqual(p["kfeta_trigger"], 3)
        self.assertEqual(p["convtrans_avglen_m"], 45.0)
        self.assertEqual(p["icloud"], 0)
        self.assertEqual(p["cldovrlp"], 4)
        self.assertEqual(p["idcor"], 1)
        self.assertEqual(p["o3input"], 0)
        self.assertEqual(p["ghg_input"], 0)
        self.assertEqual(p["rdmaxalb"], False)
        self.assertEqual(p["opt_thcnd"], 2)
        self.assertEqual(p["kim_tofd"], True)
        self.assertEqual(p["tofd_factor"], 0.005)

    def test_force_read_thompson(self):
        def phys(mp, extra):
            t = WrfTranslator(parse_namelist(real_case(f"""
            &time_control /
            &domains time_step = 60 /
            &physics mp_physics = {mp}, {extra} /
            &dynamics /
            """)))
            gen = t.translate()
            self.assertEqual(t.errors, [])
            self.assertFalse([w for w in t.warnings if "thompson" in w])
            return gen["phys"]
        self.assertTrue(phys(8, "force_read_thompson = .true.")["force_read_thompson"])
        self.assertNotIn("force_read_thompson", phys(8, "force_read_thompson = .false."))
        self.assertNotIn("force_read_thompson", phys(6, "force_read_thompson = .true."))  # read only by Thompson
        self.assertNotIn("write_thompson_tables", phys(8, "write_thompson_tables = .false."))

    def test_opt_thcnd(self):
        base = """
        &time_control /
        &domains time_step = 60 /
        &physics sf_sfclay_physics = 1, sf_surface_physics = 2, num_soil_layers = 4, bl_pbl_physics = 1,
         opt_thcnd = {} /
        &dynamics /
        """
        t = WrfTranslator(parse_namelist(real_case(base.format(2))))
        gen = t.translate()
        self.assertEqual(t.errors, [])
        self.assertEqual(gen["phys"]["opt_thcnd"], 2)
        t = WrfTranslator(parse_namelist(real_case(base.format(1))))
        gen = t.translate()
        self.assertEqual(t.errors, [])
        self.assertNotIn("opt_thcnd", gen["phys"])
        t = WrfTranslator(parse_namelist(real_case(base.format(3))))
        t.translate()
        self.assertTrue(any("opt_thcnd = 3" in e for e in t.errors))

    def test_unported_dynamics_switches_rejected(self):
        for key, val in (("diff_6th_slopeopt", "1"), ("sfs_opt", "1"), ("m_opt", "1"), ("tracer_mix6_off", ".true.")):
            text = f"&time_control /\n&domains time_step = 60 /\n&physics /\n&dynamics {key} = {val} /\n"
            t = WrfTranslator(parse_namelist(real_case(text)))
            t.translate()
            self.assertTrue(any(key in e for e in t.errors), key)
        t = WrfTranslator(parse_namelist(real_case("&time_control /\n&domains time_step = 60 /\n&physics /\n&dynamics sfs_opt = 0 /\n")))
        t.translate()
        self.assertEqual(t.errors, [])

    def test_seaice_parameters_rejected(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &physics sf_sfclay_physics = 1, sf_surface_physics = 2, num_soil_layers = 4, bl_pbl_physics = 1,
         seaice_albedo_default = {} /
        &dynamics /
        """
        t = WrfTranslator(parse_namelist(real_case(text.format(0.65))))
        t.translate()
        self.assertEqual(t.errors, [])
        t = WrfTranslator(parse_namelist(real_case(text.format(0.8))))
        t.translate()
        self.assertTrue(any("seaice_albedo_default = 0.8" in e for e in t.errors))

    def test_extended_dynamics_options_mapping(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &dynamics
         diff_opt = 2,
         km_opt = 2,
         c_s = 0.18,
         c_k = 0.12,
         mix_isotropic = 1,
         mix_upper_bound = 0.2,
         mix_full_fields = .true.,
         tke_drag_coefficient = 0.01,
         tke_heat_flux = 0.05,
         tke_upper_bound = 500.0,
         moist_mix2_off = .true.,
         scalar_mix2_off = .true.,
         tke_mix2_off = .true.,
         moist_mix6_off = .true.,
         scalar_mix6_off = .true.,
         tke_mix6_off = .true.,
         w_damping = 1,
         w_crit_cfl = 1.2,
         scalar_adv_opt = 0,
         tke_adv_opt = 0,
         epssm = 0.2,
         smdiv = 0.15,
         emdiv = 0.02,
         etac = 0.25
        /
        """
        nml = parse_namelist(real_case(text))
        t = WrfTranslator(nml)
        gen = t.translate()
        self.assertEqual(t.errors, [])
        d = gen["dyn"]
        self.assertEqual(d["c_s"], 0.18)
        self.assertEqual(d["c_k"], 0.12)
        self.assertEqual(d["mix_isotropic"], 1)
        self.assertEqual(d["mix_upper_bound"], 0.2)
        self.assertEqual(d["mix_full_fields"], True)
        self.assertEqual(d["tke_drag_coefficient"], 0.01)
        self.assertEqual(d["tke_heat_flux"], 0.05)
        self.assertEqual(d["tke_upper_bound"], 500.0)
        self.assertEqual(d["moist_mix2_off"], True)
        self.assertEqual(d["scalar_mix2_off"], True)
        self.assertEqual(d["tke_mix2_off"], True)
        self.assertEqual(d["moist_mix6_off"], True)
        self.assertEqual(d["scalar_mix6_off"], True)
        self.assertEqual(d["tke_mix6_off"], True)
        self.assertEqual(d["w_crit_cfl"], 1.2)
        self.assertEqual(d["scalar_adv_opt"], 0)
        self.assertEqual(d["tke_adv_opt"], 0)
        self.assertEqual(d["epssm"], 0.2)
        self.assertEqual(d["smdiv"], 0.15)
        self.assertEqual(d["emdiv"], 0.02)
        self.assertEqual(d["etac"], 0.25)

    def test_physics_suite_conus(self):
        """WRF's default test/em_real namelist: physics_suite CONUS, -1 options, gwd_opt and fractional_seaice."""
        text = ("&time_control /\n&domains time_step = 90, max_dom = 1 /\n"
                "&physics physics_suite = 'CONUS', mp_physics = -1, cu_physics = -1, ra_lw_physics = -1,\n"
                " ra_sw_physics = -1, bl_pbl_physics = -1, sf_sfclay_physics = -1, sf_surface_physics = -1,\n"
                " radt = 15, fractional_seaice = 1 /\n&dynamics gwd_opt = 1, diff_opt = 2, km_opt = 4 /\n")
        t = WrfTranslator(parse_namelist(real_case(text)))
        sec = t.translate()
        self.assertEqual(t.errors, [])
        expect = {"mp_physics": 8, "cu_physics": 6, "ra_lw_physics": 4, "ra_sw_physics": 4, "bl_pbl_physics": 2,
                  "sf_sfclay_physics": 2, "sf_surface_physics": 2, "fractional_seaice": 1}
        for k, v in expect.items():
            self.assertEqual(sec["phys"][k], v, k)
        self.assertEqual(sec["dyn"]["gwd_opt"], 1)

    def test_sst_update(self):
        """sst_update = 1: the auxinput4 interval and a custom wrflowinp name."""
        text = ("&time_control auxinput4_inname = 'lowinp_d<domain>', auxinput4_interval = 180,\n"
                " io_form_auxinput4 = 2 /\n&domains time_step = 60 /\n&physics sst_update = 1, mp_physics = 6, sf_sfclay_physics = 1,\n"
                " sf_surface_physics = 2, bl_pbl_physics = 1 /\n")
        t = WrfTranslator(parse_namelist(real_case(text)))
        sec = t.translate()
        self.assertEqual(t.errors, [])
        self.assertEqual(sec["phys"]["sst_update"], 1)
        self.assertEqual(sec["phys"]["auxinput4_interval"], 180)
        self.assertEqual(sec["real"]["wrflowinp"], "lowinp_d01")

    def test_physics_suite_tropical_keeps_explicit_values(self):
        text = ("&time_control /\n&domains time_step = 60 /\n"
                "&physics physics_suite = 'tropical', cu_physics = 1 /\n")
        t = WrfTranslator(parse_namelist(real_case(text)))
        sec = t.translate()
        self.assertEqual(t.errors, [])
        self.assertEqual(sec["phys"]["mp_physics"], 6)
        self.assertEqual(sec["phys"]["cu_physics"], 1)
        self.assertEqual(sec["phys"]["bl_pbl_physics"], 1)
        self.assertEqual(sec["phys"]["sf_sfclay_physics"], 1)

    def test_physics_suite_errors(self):
        for fragment, expected_err in [
                ("&physics mp_physics = -1 /", "needs a physics_suite"),
                ("&physics physics_suite = 'arctic' /", "not recognized"),
                ("&dynamics gwd_opt = 3 /", "gwd_opt = 3"),
                ("&physics sf_sfclay_physics = 1, sf_surface_physics = 2, fractional_seaice = 2 /",
                 "fractional_seaice = 2"),
                ("&physics sf_sfclay_physics = 1, sf_surface_physics = 2, bl_pbl_physics = 2 /",
                 "use myjsfc scheme"),
                ("&physics sf_sfclay_physics = 2, sf_surface_physics = 2, bl_pbl_physics = 1 /",
                 "use sfclay scheme"),
                ("&physics sf_sfclay_physics = 2, sf_surface_physics = 2, bl_pbl_physics = 2, isfflx = 0 /",
                 "uninitialised memory")]:
            t = WrfTranslator(parse_namelist(real_case(f"&time_control /\n&domains time_step = 60 /\n{fragment}\n")))
            t.translate()
            self.assertIn(expected_err, " ".join(t.errors), fragment)

    def test_nest_gwd_opt_equal_or_zero(self):
        """check_a_mundo: a nest's gwd_opt is domain 1's or 0."""
        doms = ("&domains time_step = 60, max_dom = 2, e_we = 100, 61, e_sn = 100, 61, e_vert = 45, 45,"
                " parent_id = 0, 1, i_parent_start = 1, 30, j_parent_start = 1, 30,"
                " parent_grid_ratio = 1, 3, parent_time_step_ratio = 1, 3 /\n")
        for gwd, bad in (("0, 1", True), ("1, 0", False), ("1", False), ("1, 1", False)):
            t = WrfTranslator(parse_namelist(f"&time_control /\n{doms}&dynamics gwd_opt = {gwd} /\n"))
            t.translate()
            self.assertEqual(any("gwd_opt must be equal for all domains" in e for e in t.errors), bad, gwd)

    def test_unsupported_physics_and_dynamics_errors(self):
        cases = [
            # share/module_check_a_mundo.F: "bl_mynn_edmf > 0 requires both shcu_physics=0 & ishallow=0"
            ("&physics mp_physics = 6, sf_sfclay_physics = 5, bl_pbl_physics = 5, ishallow = 1 /",
             "requires shcu_physics = 0 and ishallow = 0"),
            ("&physics cu_physics = 3, clos_choice = 99 /", "clos_choice = 99"),
            ("&physics cu_physics = 1, kfeta_trigger = 2 /", "kfeta_trigger = 2"),
            ("&physics sf_surface_physics = 2, num_soil_layers = 5 /", "num_soil_layers = 5"),
            ("&physics sf_surface_physics = 4, opt_crop = 1 /", "opt_crop = 1"),
            ("&physics ra_lw_physics = 4, ra_sw_physics = 0 /", "requires both ra_lw_physics = 4 and ra_sw_physics = 4"),
            ("&physics ra_lw_physics = 4, ra_sw_physics = 4, aer_opt = 1 /", "aer_opt = 1"),
            ("&physics ra_lw_physics = 4, ra_sw_physics = 4, cldovrlp = 9 /", "cldovrlp = 9"),
            ("&physics cu_rad_feedback = .true. /", "cu_rad_feedback"),
            ("&physics shcu_physics = 2 /", "shcu_physics"),
            ("&physics topo_wind = 1 /", "topo_wind"),
            ("&dynamics moist_adv_opt = 2 /", "moist_adv_opt = 2"),
            ("&dynamics diff_6th_opt = 3 /", "diff_6th_opt = 3"),
            ("&dynamics h_mom_adv_order = 5, h_sca_adv_order = 6 /", "h_sca_adv_order = 6"),
            ("&dynamics hybrid_opt = 9 /", "hybrid_opt = 9"),
        ]
        for fragment, expected_err in cases:
            text = f"&time_control /\n&domains time_step = 60 /\n{fragment}\n"
            nml = parse_namelist(real_case(text))
            t = WrfTranslator(nml)
            t.translate()
            joined = " ".join(t.errors)
            self.assertIn(expected_err, joined, f"Expected '{expected_err}' in errors for {fragment}, got {t.errors}")

    def test_real_case_needs_specified(self):
        # specified defaults to .false. (Registry); a real case without it has no lateral boundary condition on
        # domain 1 and wrf.exe stops in share/module_bc.F boundary_condition_check
        text = "&time_control /\n&domains time_step = 60 /\n"
        t = WrfTranslator(parse_namelist(text))
        t.translate()
        self.assertEqual(len(t.errors), 1, t.errors)
        self.assertIn("pass --case", t.errors[0])
        t = WrfTranslator(parse_namelist(text), override_case="real")
        t.translate()
        self.assertEqual(len(t.errors), 1, t.errors)
        self.assertIn("specified = .false. in a real case", t.errors[0])
        t = WrfTranslator(parse_namelist(real_case(text)), override_case="real")
        self.assertEqual(t.translate()["run"]["case"], "real")
        self.assertEqual(t.errors, [])

    def test_ishallow_without_grell_freitas(self):
        # WRF passes ishallow to the G3 and GF schemes only (phys/module_cumulus_driver.F) and check_a_mundo has no
        # rule for it with other cumulus schemes: Kain-Fritsch with ishallow = 1 runs as without it
        t = WrfTranslator(parse_namelist(real_case("&time_control /\n&domains time_step = 60 /\n"
                                                   "&physics cu_physics = 1, ishallow = 1 /\n")))
        t.translate()
        self.assertEqual(t.errors, [])

    def test_invalid_nest_configurations_errors(self):
        base = """
        &time_control input_from_file = .true., .true. /
        &domains
         time_step = 60, max_dom = 2,
         e_we = 120, {e_we}, e_sn = 90, 91, e_vert = 45, {e_vert},
         parent_id = 0, 1, i_parent_start = 1, {i_start}, j_parent_start = 1, 30,
         parent_grid_ratio = 1, 3, parent_time_step_ratio = 1, 3{extra}
        /
        &bdy_control specified = .true., .false. /
        """

        def errors(e_we=121, e_vert=45, i_start=40, extra=""):
            t = WrfTranslator(parse_namelist(base.format(e_we=e_we, e_vert=e_vert, i_start=i_start, extra=extra)))
            t.translate()
            return t.errors

        self.assertEqual(errors(), [])
        # 1. (e_we - s_we) not a multiple of parent_grid_ratio: "Nested dimensions are illegal" (dyn_em/start_em.F)
        e = errors(e_we=120)
        self.assertEqual(len(e), 1, e)
        self.assertIn("not divisible by parent_grid_ratio", e[0])

        # 2. the nest (40 parent cells from i_parent_start) must end 2 parent cells inside d01 (nx = 119): SINT's
        # stencil reaches 2 cells beyond the nest; 78 + 40 - 1 = 117 is the last allowed end
        self.assertEqual(errors(i_start=78), [])
        for i_start in (79, 100):
            e = errors(i_start=i_start)
            self.assertEqual(len(e), 1, e)
            self.assertIn("reaches within 2 cells of parent domain d01's edge", e[0])

        # 3. Mismatched e_vert and invalid smooth_option
        e = errors(e_vert=60, extra=",\n         smooth_option = 5")
        self.assertTrue(any("e_vert" in x for x in e), e)
        self.assertTrue(any("smooth_option = 5" in x for x in e), e)


class TestExecution2Steps(unittest.TestCase):
    """Verify that generated INI executes for 2 steps (binary: $CPPWRF_BIN, default build/cppwrf)."""

    def test_run_r0_2steps(self):
        exe = os.environ.get("CPPWRF_BIN", os.path.join(REPO_ROOT, "build", "cppwrf"))
        if not os.path.exists(exe):
            self.skipTest(f"{exe} not found (set CPPWRF_BIN)")

        nml_path = wrf_runs_path("run_R0_r8/namelist.input")
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

    def test_run_n1_2steps(self):
        exe = os.environ.get("CPPWRF_BIN", os.path.join(REPO_ROOT, "build", "cppwrf"))
        if not os.path.exists(exe):
            self.skipTest(f"{exe} not found (set CPPWRF_BIN)")

        nml_path = wrf_runs_path("run_N1_r8/steps/namelist.input")
        if not os.path.exists(nml_path):
            self.skipTest(f"Namelist not found: {nml_path}")

        with tempfile.NamedTemporaryFile(suffix=".ini", mode="w", delete=False) as tf:
            ini_path = tf.name

        try:
            cmd_trans = [sys.executable, os.path.join(REPO_ROOT, "tools", "wrf2cppwrf.py"),
                         nml_path, "-o", ini_path]
            ret_trans = subprocess.run(cmd_trans, capture_output=True, text=True)
            self.assertEqual(ret_trans.returncode, 0, f"Translator failed: {ret_trans.stderr}")

            cmd_run = [os.path.join(REPO_ROOT, "tools", "heavy.sh"), exe, ini_path,
                       "--set", "run.run_time=120"]
            ret_run = subprocess.run(cmd_run, capture_output=True, text=True, cwd=REPO_ROOT)
            self.assertEqual(ret_run.returncode, 0, f"Nested run failed: {ret_run.stderr}\n{ret_run.stdout}")
        finally:
            if os.path.exists(ini_path):
                os.remove(ini_path)


if __name__ == "__main__":
    unittest.main()
