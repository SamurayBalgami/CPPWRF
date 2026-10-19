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


def wrf_runs_path(rel: str) -> str:
    """Path of a WRF reference run file under ${CPPWRF_WRF_RUNS} ('' when the variable is unset)."""
    runs = os.environ.get("CPPWRF_WRF_RUNS")
    return os.path.join(runs, rel) if runs else ""

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
        self.assertEqual(gen["d02"]["dt"], 20)
        self.assertEqual(gen["d02"]["out_interval"], 30)
        self.assertEqual(gen["d02"]["diag_interval"], 30)
        self.assertEqual(gen["d03"]["parent_id"], 2)
        self.assertEqual(gen["d03"]["i_parent_start"], 30)
        self.assertEqual(gen["d03"]["j_parent_start"], 20)
        self.assertEqual(gen["d03"]["parent_grid_ratio"], 5)
        self.assertEqual(gen["d03"]["parent_time_step_ratio"], 5)
        self.assertEqual(gen["d03"]["wrfinput"], "wrfinput_d03")
        self.assertEqual(gen["d03"]["dt"], 4)
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
                if gv != rv and gv.endswith(".0") or rv.endswith(".0"):
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
        nml_path = os.path.join(os.environ.get("CPPWRF_WRF_REF", ""), "run_hill_p0/namelist.input")
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
        nml_path = os.path.join(os.environ.get("CPPWRF_WRF_REF", ""), "run_qss/namelist.input")
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
        nml_path = os.path.join(os.environ.get("CPPWRF_WRF_REF", ""), "run_straka_d1/namelist.input")
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
        &dynamics diff_opt = 3 /
        """
        nml = parse_namelist(text)
        t = WrfTranslator(nml)
        t.translate()
        self.assertTrue(any("diff_opt = 3" in w for w in t.warnings))

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

    def test_noahmp_mapping(self):
        text = """
        &time_control /
        &domains time_step = 60 /
        &physics sf_sfclay_physics = 1, sf_surface_physics = 4, bl_pbl_physics = 1 /
        &noah_mp dveg = 4, opt_runsrf = 3, opt_runsub = 3 /
        &dynamics /
        """
        nml = parse_namelist(text)
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
         kim_tofd = .true.,
         tofd_factor = 0.005
        /
        &dynamics /
        """
        nml = parse_namelist(text)
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
        self.assertEqual(p["kim_tofd"], True)
        self.assertEqual(p["tofd_factor"], 0.005)

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
        nml = parse_namelist(text)
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
        t = WrfTranslator(parse_namelist(text))
        sec = t.translate()
        self.assertEqual(t.errors, [])
        expect = {"mp_physics": 8, "cu_physics": 6, "ra_lw_physics": 4, "ra_sw_physics": 4, "bl_pbl_physics": 2,
                  "sf_sfclay_physics": 2, "sf_surface_physics": 2, "fractional_seaice": 1}
        for k, v in expect.items():
            self.assertEqual(sec["phys"][k], v, k)
        self.assertEqual(sec["dyn"]["gwd_opt"], 1)

    def test_physics_suite_tropical_keeps_explicit_values(self):
        text = ("&time_control /\n&domains time_step = 60 /\n"
                "&physics physics_suite = 'tropical', cu_physics = 1 /\n")
        t = WrfTranslator(parse_namelist(text))
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
                 "fractional_seaice = 2")]:
            t = WrfTranslator(parse_namelist(f"&time_control /\n&domains time_step = 60 /\n{fragment}\n"))
            t.translate()
            self.assertIn(expected_err, " ".join(t.errors), fragment)

    def test_unsupported_physics_and_dynamics_errors(self):
        cases = [
            ("&physics cu_physics = 1, ishallow = 1 /", "ishallow = 1 requires cu_physics = 3"),
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
            nml = parse_namelist(text)
            t = WrfTranslator(nml)
            t.translate()
            joined = " ".join(t.errors)
            self.assertIn(expected_err, joined, f"Expected '{expected_err}' in errors for {fragment}, got {t.errors}")

    def test_invalid_nest_configurations_errors(self):
        # 1. Non-divisible child extent (e_we - 1) % parent_grid_ratio != 0
        text1 = """
        &domains
         time_step = 60, max_dom = 2,
         e_we = 120, 120, e_sn = 90, 91, e_vert = 45, 45,
         parent_id = 0, 1, i_parent_start = 1, 40, j_parent_start = 1, 30,
         parent_grid_ratio = 1, 3, parent_time_step_ratio = 1, 3
        /
        &bdy_control specified = .true., .false. /
        """
        t1 = WrfTranslator(parse_namelist(text1))
        t1.translate()
        self.assertTrue(any("not divisible by parent_grid_ratio" in e for e in t1.errors))

        # 2. Child extent exceeds parent domain
        text2 = """
        &domains
         time_step = 60, max_dom = 2,
         e_we = 120, 121, e_sn = 90, 91, e_vert = 45, 45,
         parent_id = 0, 1, i_parent_start = 1, 100, j_parent_start = 1, 30,
         parent_grid_ratio = 1, 3, parent_time_step_ratio = 1, 3
        /
        &bdy_control specified = .true., .false. /
        """
        t2 = WrfTranslator(parse_namelist(text2))
        t2.translate()
        self.assertTrue(any("exceeds parent domain" in e for e in t2.errors))

        # 3. Mismatched e_vert and invalid smooth_option
        text3 = """
        &domains
         time_step = 60, max_dom = 2,
         e_we = 120, 121, e_sn = 90, 91, e_vert = 45, 60,
         parent_id = 0, 1, i_parent_start = 1, 40, j_parent_start = 1, 30,
         parent_grid_ratio = 1, 3, parent_time_step_ratio = 1, 3,
         smooth_option = 5
        /
        &bdy_control specified = .true., .false. /
        """
        t3 = WrfTranslator(parse_namelist(text3))
        t3.translate()
        self.assertTrue(any("e_vert" in e for e in t3.errors))
        self.assertTrue(any("smooth_option = 5" in e for e in t3.errors))


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
