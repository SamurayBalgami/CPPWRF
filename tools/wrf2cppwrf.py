#!/usr/bin/env python3
"""WRF namelist.input to CPPWRF configuration (.ini) translator.

Converts WRF namelist.input into CPPWRF-compatible INI configuration files.
Implements a pure Python (f90nml-free) Fortran namelist parser.
Supports single-domain execution (maps domain 1, warns on multiple domains).
Categorizes and warns/errors on unsupported or non-standard WRF options.

"""

import argparse
import datetime
import math
import os
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple


# -----------------------------------------------------------------------------
# Fortran Namelist Tokenizer and Parser (f90nml-free)
# -----------------------------------------------------------------------------

class NamelistParseError(Exception):
    """Raised when namelist parsing fails."""
    pass


def _tokenize(text: str) -> List[Tuple[str, str, int]]:
    """Tokenize Fortran namelist text.
    
    Returns a list of (token_type, token_value, line_number) tuples.
    Token types:
      - 'GROUP_START': '&' followed by group name
      - 'GROUP_END': '/' or '&end'
      - 'EQUALS': '='
      - 'COMMA': ','
      - 'VALUE': unquoted literal or quoted string
    """
    tokens = []
    lines = text.splitlines()
    
    for line_idx, line in enumerate(lines, 1):
        i = 0
        n = len(line)
        while i < n:
            c = line[i]
            
            # Whitespace
            if c.isspace():
                i += 1
                continue
                
            # Comment: '!' to end of line
            if c == '!':
                break
                
            # Group start: &name
            if c == '&':
                match = re.match(r'&([A-Za-z0-9_]+)', line[i:])
                if match:
                    gname = match.group(1)
                    if gname.lower() == 'end':
                        tokens.append(('GROUP_END', '&end', line_idx))
                    else:
                        tokens.append(('GROUP_START', gname.lower(), line_idx))
                    i += match.end()
                    continue
                else:
                    raise NamelistParseError(f"Line {line_idx}: invalid group syntax near '{line[i:]}'")
                    
            # Group end: '/'
            if c == '/':
                tokens.append(('GROUP_END', '/', line_idx))
                i += 1
                continue
                
            # Equals: '='
            if c == '=':
                tokens.append(('EQUALS', '=', line_idx))
                i += 1
                continue
                
            # Comma: ','
            if c == ',':
                tokens.append(('COMMA', ',', line_idx))
                i += 1
                continue
                
            # Quoted string ('...' or "...")
            if c in ("'", '"'):
                quote = c
                start = i
                i += 1
                s_chars = []
                while i < n:
                    if line[i] == quote:
                        if i + 1 < n and line[i + 1] == quote:
                            # Escaped quote: '' or ""
                            s_chars.append(quote)
                            i += 2
                        else:
                            i += 1
                            break
                    else:
                        s_chars.append(line[i])
                        i += 1
                else:
                    raise NamelistParseError(f"Line {line_idx}: unterminated string literal starting with {quote}")
                tokens.append(('VALUE', ''.join(s_chars), line_idx))
                continue
                
            # Word / Number / Expression / Identifier
            start = i
            while i < n and not line[i].isspace() and line[i] not in (',', '=', '/', '!', '&'):
                i += 1
            raw_val = line[start:i]
            tokens.append(('VALUE', raw_val, line_idx))
            
    return tokens


def _parse_scalar_value(s: str) -> Any:
    """Parse a single token string into a Python bool, int, float, or str."""
    low = s.strip().lower()
    
    # Booleans
    if low in ('.true.', '.t.', 't'):
        return True
    if low in ('.false.', '.f.', 'f'):
        return False
        
    # Integers
    if re.fullmatch(r'[+-]?\d+', low):
        try:
            return int(low)
        except ValueError:
            pass
            
    # Floats (including Fortran exponents d/D and trailing dots like 5000.)
    f_str = re.sub(r'[dD]', 'e', low)
    if re.fullmatch(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?', f_str):
        try:
            return float(f_str)
        except ValueError:
            pass
            
    # Fallback to string as-is
    return s.strip()


def parse_namelist(text: str) -> Dict[str, Dict[str, List[Any]]]:
    """Parse namelist text into a dictionary of groups and key -> list of values.
    
    Structure: { 'group_name': { 'key_name': [val1, val2, ...] } }
    Keys and group names are normalized to lowercase.
    Handles Fortran repetition syntax (e.g. 3*0 -> [0, 0, 0]).
    """
    tokens = _tokenize(text)
    data: Dict[str, Dict[str, List[Any]]] = {}
    current_group: Optional[str] = None
    
    idx = 0
    num_tokens = len(tokens)
    
    while idx < num_tokens:
        tok_type, tok_val, line_num = tokens[idx]
        
        if tok_type == 'GROUP_START':
            current_group = tok_val
            if current_group not in data:
                data[current_group] = {}
            idx += 1
            continue
            
        elif tok_type == 'GROUP_END':
            current_group = None
            idx += 1
            continue
            
        elif tok_type == 'VALUE':
            if current_group is None:
                idx += 1
                continue
                
            key = tok_val.lower()
            if '(' in key and key.endswith(')'):
                key = key[:key.find('(')]
                
            idx += 1
            if idx >= num_tokens:
                break
                
            next_type, next_val, _ = tokens[idx]
            if next_type != 'EQUALS':
                continue
            idx += 1  # Skip '='
            
            values: List[Any] = []
            while idx < num_tokens:
                t_type, t_val, _ = tokens[idx]
                if t_type == 'GROUP_END' or t_type == 'GROUP_START':
                    break
                if t_type == 'EQUALS':
                    if values:
                        last_val = values.pop()
                        idx -= 1
                    break
                if t_type == 'COMMA':
                    idx += 1
                    continue
                if t_type == 'VALUE':
                    repeat_match = re.fullmatch(r'(\d+)\*(.*)', t_val)
                    if repeat_match:
                        count = int(repeat_match.group(1))
                        val_str = repeat_match.group(2)
                        parsed_val = _parse_scalar_value(val_str)
                        values.extend([parsed_val] * count)
                    else:
                        values.append(_parse_scalar_value(t_val))
                    idx += 1
                    
            if current_group in data:
                if key in data[current_group]:
                    data[current_group][key].extend(values)
                else:
                    data[current_group][key] = values
                    
        else:
            idx += 1
            
    return data


# -----------------------------------------------------------------------------
# Options Mapping Specification & Validation
# -----------------------------------------------------------------------------

SUPPORTED_OPTIONS = {
    "mp_physics": {
        0: "None / off",
        1: "Kessler warm-rain",
        6: "WSM6 6-class single-moment",
        8: "Thompson",
    },
    "sf_sfclay_physics": {
        0: "None / off",
        1: "Revised MM5 surface layer (sfclayrev)",
        2: "Eta similarity surface layer (MYJ)",
    },
    "sf_surface_physics": {
        0: "None / off",
        2: "Unified Noah land surface model (4 soil layers)",
    },
    "bl_pbl_physics": {
        0: "None / off",
        1: "YSU (Yonsei University) planetary boundary layer",
        2: "MYJ (Mellor-Yamada-Janjic) planetary boundary layer",
    },
    "cu_physics": {
        0: "None / off",
        1: "Kain-Fritsch cumulus parameterization",
        3: "Grell-Freitas cumulus parameterization",
        16: "New Tiedtke cumulus parameterization",
    },
    "ra_lw_physics": {
        0: "None / off",
        4: "RRTMG longwave radiation",
    },
    "ra_sw_physics": {
        0: "None / off",
        4: "RRTMG shortwave radiation",
    },
    "diff_opt": {
        0: "None / off",
        1: "Cartesian or coordinate-surface diffusion",
        2: "Physical-space (full metric) diffusion",
    },
    "km_opt": {
        1: "Constant eddy viscosity (khdif / kvdif)",
        2: "1.5-order prognostic TKE",
        4: "2D horizontal Smagorinsky",
    },
    "damp_opt": {
        0: "None / off",
        1: "Upper diffusive damping layer",
        2: "Rayleigh relaxation layer",
        3: "Implicit Rayleigh damping of w in upper layer",
    },
    "h_adv_order": {5: "5th-order upwind", 6: "6th-order centered"},
    "v_adv_order": {3: "3rd-order upwind", 4: "4th-order centered"},
}

IGNORED_WRF_OPTIONS = {
    # &time_control
    "frames_per_outfile": "Frame count per NetCDF file handled by output format",
    "debug_level": "WRF internal logging verbosity",
    "input_from_file": "Real case initialization always reads wrfinput",
    "iofields_filename": "Custom NetCDF field selection",
    "io_form_history": "History output format specified in [output] format",
    "io_form_restart": "Restart format handled via NetCDF",
    "io_form_input": "Input format handled via NetCDF",
    "io_form_boundary": "Boundary format handled via NetCDF",
    "start_year": "Start date component",
    "start_month": "Start date component",
    "start_day": "Start date component",
    "start_hour": "Start date component",
    "start_minute": "Start date component",
    "start_second": "Start date component",
    "end_year": "End date component",
    "end_month": "End date component",
    "end_day": "End date component",
    "end_hour": "End date component",
    "end_minute": "End date component",
    "end_second": "End date component",
    "interval_seconds": "Lateral boundary interval (read directly from wrfbdy_d01 NetCDF in real case)",
    # &domains
    "time_step_fract_num": "Combined into time.dt if non-zero",
    "time_step_fract_den": "Combined into time.dt if non-zero",
    "grid_id": "Domain grid hierarchy identifier (domain 1 only)",
    "parent_id": "Domain nesting parent identifier",
    "i_parent_start": "Domain nesting starting index",
    "j_parent_start": "Domain nesting starting index",
    "parent_grid_ratio": "Nesting grid refinement ratio",
    "parent_time_step_ratio": "Nesting timestep ratio",
    "feedback": "Two-way nesting feedback",
    "smooth_option": "Nesting boundary smoothing",
    "num_metgrid_levels": "Metgrid vertical levels used by real.exe",
    "num_metgrid_soil_levels": "Metgrid soil levels used by real.exe",
    "e_we": "Real case grid dimension read directly from wrfinput_d01",
    "e_sn": "Real case grid dimension read directly from wrfinput_d01",
    "e_vert": "Real case grid dimension read directly from wrfinput_d01",
    "p_top_requested": "Real case p_top read directly from wrfinput_d01",
    "dx": "Real case grid spacing read directly from wrfinput_d01",
    "dy": "Real case grid spacing read directly from wrfinput_d01",
    # &physics
    "isfflx": "Surface heat/moisture fluxes flag",
    "ifsnow": "Snow effects flag",
    "icloud": "Cloud fractions for radiation",
    "num_soil_layers": "Fixed at 4 for Noah LSM in CPPWRF",
    "sf_urban_physics": "Urban canopy model",
    # &dynamics
    "base_temp": "Base potential temperature (read from sounding or wrfinput)",
    "mix_full_fields": "Diffusion on full fields",
    "hybrid_opt": "Real case hybrid coordinate flag read directly from wrfinput_d01 global attribute",
    "h_sca_adv_order": "Scalar horizontal advection order (coupled with h_mom_adv_order in CPPWRF)",
    "v_sca_adv_order": "Scalar vertical advection order (coupled with v_mom_adv_order in CPPWRF)",
    "epssm": "Acoustic divergence damping parameter (default 0.1)",
    "smdiv": "Acoustic divergence damping parameter (default 0.1)",
    "emdiv": "External mode filter coefficient (default 0.01)",
    "diff_6th_factor": "Sixth-order diffusion factor (default 0.12)",
    # &bdy_control
    "spec_bdy_width": "Real case lateral boundary width (must be 5; read from wrfbdy_d01)",
    # &fdda, &grib2, &namelist_quilt
    "nio_tasks_per_group": "WRF quilt I/O server",
    "nio_groups": "WRF quilt I/O server",
}


def _format_num(val: Any) -> str:
    """Format integer or float cleanly for INI."""
    if isinstance(val, bool):
        return "true" if val else "false"
    if isinstance(val, int):
        return str(val)
    if isinstance(val, float):
        if val.is_integer():
            return str(int(val))
        return f"{val:g}"
    return str(val)


# -----------------------------------------------------------------------------
# Translation Engine
# -----------------------------------------------------------------------------

class WrfTranslator:
    """Translates parsed WRF namelist into CPPWRF INI configuration."""

    def __init__(self, namelist: Dict[str, Dict[str, List[Any]]],
                 namelist_path: str = "",
                 override_case: Optional[str] = None,
                 real_dir: Optional[str] = None,
                 table_dir: Optional[str] = None,
                 output_format: str = "native",
                 explicit_files: bool = False):
        self.nml = namelist
        self.namelist_path = os.path.abspath(namelist_path) if namelist_path else ""
        self.override_case = override_case
        self.real_dir = real_dir
        self.table_dir = table_dir
        self.output_format = output_format
        self.explicit_files = explicit_files

        self.warnings: List[str] = []
        self.errors: List[str] = []
        self.notices: List[str] = []
        self.used_keys: Set[Tuple[str, str]] = set()

    def warn(self, msg: str):
        self.warnings.append(msg)

    def error(self, msg: str):
        self.errors.append(msg)

    def notice(self, msg: str):
        self.notices.append(msg)

    def get_opt(self, group: str, key: str, default: Any = None) -> Any:
        grp = self.nml.get(group, {})
        self.used_keys.add((group, key))
        vals = grp.get(key, [])
        if not vals:
            return default
        if len(vals) > 1:
            self.warn(f"Option &{group}.{key} has {len(vals)} values {vals}; selecting domain 1 value ({vals[0]}).")
        return vals[0]

    def has_opt(self, group: str, key: str) -> bool:
        grp = self.nml.get(group, {})
        return key in grp and len(grp[key]) > 0

    def translate(self) -> Dict[str, Dict[str, Any]]:
        """Perform the translation and return INI sections dict."""
        sections: Dict[str, Dict[str, Any]] = {
            "run": {},
        }

        # 1. Check for multiple domains in &domains
        max_dom = self.get_opt("domains", "max_dom", 1)
        if isinstance(max_dom, int) and max_dom > 1:
            self.warn(f"Multiple domains specified in namelist (max_dom = {max_dom}). CPPWRF supports domain 1 only.")

        # 2. Determine Case
        ideal_grp = self.nml.get("ideal", {})
        ideal_case = self.get_opt("ideal", "ideal_case", None)
        specified = self.get_opt("bdy_control", "specified", False)

        case_name = "real"
        if self.override_case:
            case_name = self.override_case
        elif ideal_case is not None:
            ideal_map = {
                0: "real",
                1: "hill2d",
                2: "quarter_ss",
                4: "squall2d",
                6: "straka",
            }
            if ideal_case in ideal_map:
                case_name = ideal_map[ideal_case]
            else:
                self.error(f"WRF ideal_case = {ideal_case} has no direct CPPWRF equivalent.")
                case_name = f"ideal_{ideal_case}"
        elif specified:
            case_name = "real"
        else:
            dir_name = os.path.basename(os.path.dirname(self.namelist_path)).lower() if self.namelist_path else ""
            if "hill" in dir_name:
                case_name = "hill2d"
            elif "qss" in dir_name or "quarter" in dir_name:
                case_name = "quarter_ss"
            elif "squall" in dir_name:
                case_name = "squall2d"
            elif "straka" in dir_name or "grav" in dir_name:
                case_name = "straka"
            else:
                case_name = "real"

        sections["run"]["case"] = case_name

        # 3. Determine run_time
        run_days = self.get_opt("time_control", "run_days", 0)
        run_hours = self.get_opt("time_control", "run_hours", 0)
        run_minutes = self.get_opt("time_control", "run_minutes", 0)
        run_seconds = self.get_opt("time_control", "run_seconds", 0)
        total_s = (run_days * 86400) + (run_hours * 3600) + (run_minutes * 60) + run_seconds

        if total_s <= 0:
            s_yr = self.get_opt("time_control", "start_year", 2000)
            s_mo = self.get_opt("time_control", "start_month", 1)
            s_dy = self.get_opt("time_control", "start_day", 1)
            s_hr = self.get_opt("time_control", "start_hour", 0)
            s_mi = self.get_opt("time_control", "start_minute", 0)
            s_se = self.get_opt("time_control", "start_second", 0)

            e_yr = self.get_opt("time_control", "end_year", s_yr)
            e_mo = self.get_opt("time_control", "end_month", s_mo)
            e_dy = self.get_opt("time_control", "end_day", s_dy)
            e_hr = self.get_opt("time_control", "end_hour", s_hr)
            e_mi = self.get_opt("time_control", "end_minute", s_mi)
            e_se = self.get_opt("time_control", "end_second", s_se)

            try:
                t0 = datetime.datetime(s_yr, s_mo, s_dy, s_hr, s_mi, s_se)
                t1 = datetime.datetime(e_yr, e_mo, e_dy, e_hr, e_mi, e_se)
                total_s = int((t1 - t0).total_seconds())
            except Exception as e:
                self.warn(f"Could not compute run duration from start/end dates: {e}")
                total_s = 3600

        sections["run"]["run_time"] = total_s

        # 4. Physics Options (evaluated early to determine physics suite activity)
        mp_physics = self.get_opt("physics", "mp_physics", 0)
        sf_sfclay = self.get_opt("physics", "sf_sfclay_physics", 0)
        sf_surface = self.get_opt("physics", "sf_surface_physics", 0)
        bl_pbl = self.get_opt("physics", "bl_pbl_physics", 0)
        cu_physics = self.get_opt("physics", "cu_physics", 0)
        ra_lw = self.get_opt("physics", "ra_lw_physics", 0)
        ra_sw = self.get_opt("physics", "ra_sw_physics", 0)
        radt = self.get_opt("physics", "radt", 0)
        bldt = self.get_opt("physics", "bldt", 0)
        cudt = self.get_opt("physics", "cudt", 0)

        has_physics_suite = (sf_sfclay != 0 or sf_surface != 0 or bl_pbl != 0 or
                             cu_physics != 0 or ra_lw != 0 or ra_sw != 0)

        # 5. Real Data Section [real]
        if case_name == "real":
            sections["real"] = {}
            r_dir = self.real_dir
            if not r_dir:
                r_dir = os.path.dirname(self.namelist_path) if self.namelist_path else "."
            r_dir_abs = os.path.abspath(os.path.expanduser(r_dir))
            wrf_runs = os.environ.get("CPPWRF_WRF_RUNS")
            if wrf_runs:
                wrf_runs_abs = os.path.abspath(os.path.expanduser(wrf_runs))
                if r_dir_abs == wrf_runs_abs:
                    sections["real"]["dir"] = "${CPPWRF_WRF_RUNS}"
                elif r_dir_abs.startswith(wrf_runs_abs + os.sep):
                    rel = os.path.relpath(r_dir_abs, wrf_runs_abs)
                    sections["real"]["dir"] = f"${{CPPWRF_WRF_RUNS}}/{rel}"
                else:
                    sections["real"]["dir"] = r_dir_abs
            else:
                sections["real"]["dir"] = r_dir_abs

            # In configs/real_p1..p4, wrfinput/wrfbdy were written explicitly
            if self.explicit_files or has_physics_suite:
                sections["real"]["wrfinput"] = "wrfinput_d01"
                sections["real"]["wrfbdy"] = "wrfbdy_d01"

        # 6. Domain Dimensions & Boundary Conditions [domain] (ideal cases)
        if case_name != "real":
            sections["domain"] = {}
            e_we = self.get_opt("domains", "e_we", 100)
            e_sn = self.get_opt("domains", "e_sn", 3)
            e_vert = self.get_opt("domains", "e_vert", 41)
            s_we = self.get_opt("domains", "s_we", 1)
            s_sn = self.get_opt("domains", "s_sn", 1)
            s_vert = self.get_opt("domains", "s_vert", 1)

            nx = e_we - s_we
            nz = e_vert - s_vert
            ny = 1 if (e_sn - s_sn <= 2) else (e_sn - s_sn)

            sections["domain"]["nx"] = nx
            sections["domain"]["ny"] = ny
            sections["domain"]["nz"] = nz

            dx = self.get_opt("domains", "dx", 1000.0)
            dy = self.get_opt("domains", "dy", dx)
            sections["domain"]["dx"] = dx
            if ny > 1 or case_name in ("hill2d", "quarter_ss"):
                sections["domain"]["dy"] = dy

            ztop = self.get_opt("domains", "ztop", None)
            if ztop is not None:
                sections["domain"]["ztop"] = ztop

            p_top = self.get_opt("domains", "p_top_requested", self.get_opt("domains", "p_top", None))
            if p_top is not None:
                sections["domain"]["p_top"] = p_top

            per_x = self.get_opt("bdy_control", "periodic_x", True)
            sym_xs = self.get_opt("bdy_control", "symmetric_xs", False)
            open_xs = self.get_opt("bdy_control", "open_xs", False)

            if per_x:
                sections["domain"]["bc_x"] = "periodic"
            elif sym_xs:
                sections["domain"]["bc_x"] = "symmetric"
            elif open_xs:
                sections["domain"]["bc_x"] = "specified"

            per_y = self.get_opt("bdy_control", "periodic_y", True)
            sym_ys = self.get_opt("bdy_control", "symmetric_ys", False)
            open_ys = self.get_opt("bdy_control", "open_ys", False)

            if ny > 1 or "bc_x" in sections["domain"]:
                if per_y:
                    sections["domain"]["bc_y"] = "periodic"
                elif sym_ys:
                    sections["domain"]["bc_y"] = "symmetric"
                elif open_ys:
                    sections["domain"]["bc_y"] = "specified"

        # 7. Ideal Case Def Section [case]
        if case_name == "hill2d":
            sections["case"] = {
                "sounding_file": "configs/input_sounding_hill2d",
                "hm": 100.0,
                "xa": 5.0,
                "u0": 10.0,
            }
        elif case_name == "quarter_ss":
            sections["case"] = {
                "sounding_file": "configs/input_sounding_qss",
                "delt": 3.0,
            }
        elif case_name == "straka":
            sections["case"] = {
                "xc": 25550,
            }

        # 8. Time Section [time]
        sections["time"] = {}
        dt = float(self.get_opt("domains", "time_step", 60.0))
        num = self.get_opt("domains", "time_step_fract_num", 0)
        den = self.get_opt("domains", "time_step_fract_den", 1)
        if den and den > 0 and num:
            dt += float(num) / float(den)
        sections["time"]["dt"] = int(dt) if dt.is_integer() else dt

        tss = self.get_opt("dynamics", "time_step_sound", 0)
        if tss and tss > 0:
            sections["time"]["time_step_sound"] = tss

        # 9. Dynamics Section [dyn]
        sections["dyn"] = {}
        khdif = self.get_opt("dynamics", "khdif", 0.0)
        kvdif = self.get_opt("dynamics", "kvdif", 0.0)
        sections["dyn"]["khdif"] = int(khdif) if isinstance(khdif, (int, float)) and float(khdif).is_integer() else khdif
        sections["dyn"]["kvdif"] = int(kvdif) if isinstance(kvdif, (int, float)) and float(kvdif).is_integer() else kvdif

        w_damping = self.get_opt("dynamics", "w_damping", 0)
        if w_damping != 0:
            sections["dyn"]["w_damping"] = w_damping

        damp_opt = self.get_opt("dynamics", "damp_opt", 0)
        if damp_opt != 0:
            sections["dyn"]["damp_opt"] = damp_opt
            if damp_opt not in SUPPORTED_OPTIONS["damp_opt"]:
                self.warn(f"WRF damp_opt = {damp_opt} is not supported; CPPWRF supports damp_opt 0, 1, 2 and 3.")
            zdamp = self.get_opt("dynamics", "zdamp", 5000.0)
            dampcoef = self.get_opt("dynamics", "dampcoef", 0.2)
            sections["dyn"]["zdamp"] = int(zdamp) if float(zdamp).is_integer() else zdamp
            sections["dyn"]["dampcoef"] = dampcoef
        elif case_name in ("hill2d",):
            sections["dyn"]["damp_opt"] = 0

        diff_opt = self.get_opt("dynamics", "diff_opt", 0)
        if diff_opt != 0:
            sections["dyn"]["diff_opt"] = diff_opt
            if diff_opt not in SUPPORTED_OPTIONS["diff_opt"]:
                self.warn(f"WRF diff_opt = {diff_opt} is not directly supported; CPPWRF supports diff_opt 0, 1 and 2.")

        km_opt = self.get_opt("dynamics", "km_opt", 1)
        if (diff_opt != 0 or km_opt != 1) and km_opt is not None:
            if km_opt not in (0, 1):
                sections["dyn"]["km_opt"] = km_opt
            if km_opt not in SUPPORTED_OPTIONS["km_opt"] and km_opt != 0:
                self.warn(f"WRF km_opt = {km_opt} is not supported; CPPWRF supports km_opt 1, 2 and 4.")

        c_s = self.get_opt("dynamics", "c_s", None)
        if km_opt == 4 and c_s is not None:
            sections["dyn"]["c_s"] = c_s

        diff_6th_opt = self.get_opt("dynamics", "diff_6th_opt", 0)
        if diff_6th_opt != 0:
            sections["dyn"]["diff_6th_opt"] = diff_6th_opt
            diff_6th_factor = self.get_opt("dynamics", "diff_6th_factor", 0.12)
            sections["dyn"]["diff_6th_factor"] = diff_6th_factor

        moist_adv_opt = self.get_opt("dynamics", "moist_adv_opt", 1)
        if moist_adv_opt is not None and case_name != "hill2d":
            sections["dyn"]["moist_adv_opt"] = moist_adv_opt

        # scalar_adv_opt: included when physics suite is enabled, or when != 0 for non-dry runs
        scalar_adv_opt = self.get_opt("dynamics", "scalar_adv_opt", 0)
        if has_physics_suite or (scalar_adv_opt != 0 and case_name in ("quarter_ss",)):
            sections["dyn"]["scalar_adv_opt"] = scalar_adv_opt if scalar_adv_opt != 0 else 1

        prandtl = self.get_opt("dynamics", "prandtl", None)
        if prandtl is not None and prandtl > 0:
            pr_inv = 1.0 / prandtl
            if abs(pr_inv - 3.0) > 1e-4 or case_name in ("straka", "quarter_ss"):
                sections["dyn"]["pr_inv"] = pr_inv
        elif case_name in ("straka", "quarter_ss"):
            sections["dyn"]["pr_inv"] = 3.0

        top_lid = self.get_opt("dynamics", "top_lid", False)
        if top_lid or case_name in ("quarter_ss",):
            sections["dyn"]["top_lid"] = top_lid

        h_adv = self.get_opt("dynamics", "h_mom_adv_order", 5)
        if h_adv != 5 or case_name in ("quarter_ss",):
            sections["dyn"]["h_adv_order"] = h_adv
            if h_adv not in SUPPORTED_OPTIONS["h_adv_order"]:
                self.error(f"WRF h_mom_adv_order = {h_adv} is not supported (5 or 6 only).")

        v_adv = self.get_opt("dynamics", "v_mom_adv_order", 3)
        if v_adv != 3 or case_name in ("quarter_ss",):
            sections["dyn"]["v_adv_order"] = v_adv
            if v_adv not in SUPPORTED_OPTIONS["v_adv_order"]:
                self.error(f"WRF v_mom_adv_order = {v_adv} is not supported (3 or 4 only).")

        # Validate core numerical assumptions
        non_hydro = self.get_opt("dynamics", "non_hydrostatic", True)
        if not non_hydro:
            self.error("WRF non_hydrostatic = .false. is unsupported; CPPWRF is strictly non-hydrostatic.")

        rk_ord = self.get_opt("dynamics", "rk_ord", 3)
        if rk_ord != 3:
            self.error(f"WRF rk_ord = {rk_ord} is unsupported; CPPWRF strictly uses 3rd-order Runge-Kutta (RK3).")

        use_thm = self.get_opt("dynamics", "use_theta_m", 1)
        if use_thm == 0:
            self.warn("WRF use_theta_m = 0 is unsupported; CPPWRF always uses moist potential temperature (theta_m).")

        gwd_opt = self.get_opt("dynamics", "gwd_opt", 0)
        if gwd_opt != 0:
            self.warn(f"WRF gravity wave drag (gwd_opt = {gwd_opt}) is not implemented in CPPWRF.")

        chem_adv = self.get_opt("dynamics", "chem_adv_opt", 0)
        if chem_adv != 0:
            self.warn("WRF chemistry advection (chem_adv_opt) is not supported.")

        tke_adv = self.get_opt("dynamics", "tke_adv_opt", None)
        if tke_adv is not None and km_opt == 2:
            sections["dyn"]["tke_adv_opt"] = tke_adv

        # 10. Moisture Species [moist]
        if case_name == "real":
            sections["moist"] = {}
            if mp_physics == 6:
                sections["moist"]["species"] = "qv, qc, qr, qi, qs, qg" if not has_physics_suite else "qv,qc,qr,qi,qs,qg"
            elif mp_physics == 1:
                sections["moist"]["species"] = "qv, qc, qr"
            else:
                sections["moist"]["species"] = "qv"
        elif mp_physics > 0 or case_name in ("quarter_ss", "squall2d"):
            sections["moist"] = {}
            if mp_physics == 6:
                sections["moist"]["species"] = "qv, qc, qr, qi, qs, qg"
            else:
                sections["moist"]["species"] = "qv, qc, qr"

        # 11. Physics Section [phys]
        if mp_physics != 0 or has_physics_suite:
            sections["phys"] = {}
            if mp_physics != 0:
                sections["phys"]["mp_physics"] = mp_physics
                if mp_physics not in SUPPORTED_OPTIONS["mp_physics"]:
                    self.warn(f"WRF mp_physics = {mp_physics} is not supported; CPPWRF supports 0, 1 (Kessler), 6 (WSM6).")

            if has_physics_suite:
                sections["phys"]["sf_sfclay_physics"] = sf_sfclay
                sections["phys"]["sf_surface_physics"] = sf_surface
                sections["phys"]["bl_pbl_physics"] = bl_pbl
                sections["phys"]["cu_physics"] = cu_physics
                sections["phys"]["ra_lw_physics"] = ra_lw
                sections["phys"]["ra_sw_physics"] = ra_sw

                if ra_lw != 0 or ra_sw != 0 or radt != 0:
                    sections["phys"]["radt"] = 0 if radt <= 1 else radt
                if bldt > 0:
                    sections["phys"]["bldt"] = bldt
                if cudt > 0:
                    sections["phys"]["cudt"] = cudt

                tbl_dir = self.table_dir
                if not tbl_dir:
                    wrf_runs = os.environ.get("CPPWRF_WRF_RUNS")
                    candidates = []
                    if wrf_runs:
                        candidates.append(os.path.join(wrf_runs, "wrf_run"))
                    if self.namelist_path:
                        candidates.append(os.path.join(os.path.dirname(self.namelist_path), "wrf_run"))
                        candidates.append(os.path.dirname(self.namelist_path))
                    for cand in candidates:
                        if cand and os.path.exists(os.path.join(cand, "LANDUSE.TBL")):
                            tbl_dir = cand
                            break
                    if not tbl_dir:
                        tbl_dir = os.path.join(wrf_runs, "wrf_run") if wrf_runs else "${CPPWRF_WRF_RUNS}/wrf_run"
                sections["phys"]["table_dir"] = tbl_dir

                if sf_sfclay not in SUPPORTED_OPTIONS["sf_sfclay_physics"]:
                    self.warn(f"WRF sf_sfclay_physics = {sf_sfclay} not supported (0 or 1 only).")
                if sf_surface not in SUPPORTED_OPTIONS["sf_surface_physics"]:
                    self.warn(f"WRF sf_surface_physics = {sf_surface} not supported (0 or 2 only).")
                if bl_pbl not in SUPPORTED_OPTIONS["bl_pbl_physics"]:
                    self.warn(f"WRF bl_pbl_physics = {bl_pbl} not supported (0 or 1 only).")
                if cu_physics not in SUPPORTED_OPTIONS["cu_physics"]:
                    self.warn(f"WRF cu_physics = {cu_physics} not supported (0 or 1 only).")
                if ra_lw not in SUPPORTED_OPTIONS["ra_lw_physics"]:
                    self.warn(f"WRF ra_lw_physics = {ra_lw} not supported (0 or 4 only).")
                if ra_sw not in SUPPORTED_OPTIONS["ra_sw_physics"]:
                    self.warn(f"WRF ra_sw_physics = {ra_sw} not supported (0 or 4 only).")

        # 12. Output Section [output]
        sections["output"] = {}
        hist_s = self.get_opt("time_control", "history_interval_s", 0)
        hist_m = self.get_opt("time_control", "history_interval", 0)

        dt_val = float(sections["time"]["dt"])
        if hist_s and hist_s > 0:
            out_interval = max(1, round(float(hist_s) / dt_val))
        elif hist_m and hist_m > 0:
            out_interval = max(1, round(float(hist_m * 60) / dt_val))
        else:
            out_interval = max(1, round(float(total_s) / dt_val))

        diag_interval = 1 if (out_interval == 1 or case_name in ("hill2d",)) else min(10, out_interval)
        sections["output"]["diag_interval"] = diag_interval
        sections["output"]["out_interval"] = out_interval

        if self.output_format and (self.output_format != "native" or case_name in ("hill2d", "quarter_ss")):
            sections["output"]["format"] = self.output_format

        # 13. Restart Section [restart]
        restart_on = self.get_opt("time_control", "restart", False)
        restart_m = self.get_opt("time_control", "restart_interval", 0)
        if restart_on and restart_m and restart_m > 0:
            sections["restart"] = {
                "interval": int(restart_m * 60)
            }

        # 14. Check for unused / ignored options
        for grp_name, grp_dict in self.nml.items():
            for key_name, vals in grp_dict.items():
                if (grp_name, key_name) not in self.used_keys:
                    if key_name in IGNORED_WRF_OPTIONS:
                        self.notice(f"Option &{grp_name}.{key_name} = {vals[0] if vals else None}: {IGNORED_WRF_OPTIONS[key_name]}")
                    else:
                        self.warn(f"Option &{grp_name}.{key_name} = {vals} is not recognized or supported by CPPWRF.")

        return sections

    def format_ini(self, sections: Dict[str, Dict[str, Any]]) -> str:
        """Format translated sections into standard CPPWRF INI string."""
        lines = []
        lines.append("# Translated from WRF namelist.input by tools/wrf2cppwrf.py")
        if self.namelist_path:
            lines.append(f"# Source: {self.namelist_path}")
        lines.append(f"# Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        if self.warnings:
            lines.append("#")
            lines.append("# Translation Warnings:")
            for w in self.warnings:
                lines.append(f"#   - [WARN] {w}")

        if self.notices:
            lines.append("#")
            lines.append("# Translation Notices (unmapped/ignored WRF infrastructure):")
            for n in self.notices:
                lines.append(f"#   - [INFO] {n}")

        lines.append("")

        section_order = ["run", "real", "domain", "case", "time", "dyn", "moist", "phys", "restart", "output", "tracers"]
        remaining = [s for s in sections if s not in section_order]

        for sec in section_order + remaining:
            if sec not in sections or not sections[sec]:
                continue
            lines.append(f"[{sec}]")
            for k, v in sections[sec].items():
                lines.append(f"{k} = {_format_num(v)}")
            lines.append("")

        return "\n".join(lines).rstrip() + "\n"


# -----------------------------------------------------------------------------
# CLI Entry Point
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Translate WRF namelist.input to CPPWRF .ini configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 tools/wrf2cppwrf.py ~/wrf_real/run_R0_r8/namelist.input -o run_r0.ini
  python3 tools/wrf2cppwrf.py ~/wrf_ref/run_hill_p0/namelist.input -o hill.ini
  python3 tools/wrf2cppwrf.py ~/wrf_real/run_P1/namelist.input --format wrfout
"""
    )
    parser.add_argument("namelist", help="Path to WRF namelist.input")
    parser.add_argument("-o", "--output", help="Output .ini file path (default: stdout)")
    parser.add_argument("--case", help="Explicitly specify CPPWRF case (e.g. real, hill2d, quarter_ss, straka)")
    parser.add_argument("--real-dir", help="Path to directory containing wrfinput_d01 / wrfbdy_d01 (for real case)")
    parser.add_argument("--table-dir", help="Path to physics tables directory (for Noah LSM / radiation)")
    parser.add_argument("--format", choices=["native", "wrfout", "both"], default="native",
                        help="Output format (native binary, WRF NetCDF wrfout, or both)")
    parser.add_argument("--explicit-files", action="store_true",
                        help="Explicitly write wrfinput = wrfinput_d01 and wrfbdy = wrfbdy_d01 in [real]")
    parser.add_argument("--strict", action="store_true",
                        help="Treat warnings as errors and exit with non-zero code")

    args = parser.parse_args()

    if not os.path.exists(args.namelist):
        sys.stderr.write(f"Error: namelist file not found: {args.namelist}\n")
        sys.exit(2)

    try:
        with open(args.namelist, "r", encoding="utf-8", errors="replace") as f:
            text = f.read()
    except Exception as e:
        sys.stderr.write(f"Error reading namelist: {e}\n")
        sys.exit(2)

    try:
        nml = parse_namelist(text)
    except Exception as e:
        sys.stderr.write(f"Error parsing Fortran namelist: {e}\n")
        sys.exit(2)

    translator = WrfTranslator(
        nml,
        namelist_path=args.namelist,
        override_case=args.case,
        real_dir=args.real_dir,
        table_dir=args.table_dir,
        output_format=args.format,
        explicit_files=args.explicit_files
    )

    sections = translator.translate()
    ini_content = translator.format_ini(sections)

    for notice in translator.notices:
        sys.stderr.write(f"[INFO] {notice}\n")
    for warn in translator.warnings:
        sys.stderr.write(f"[WARNING] {warn}\n")
    for err in translator.errors:
        sys.stderr.write(f"[ERROR] {err}\n")

    if args.output:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(ini_content)
    else:
        sys.stdout.write(ini_content)

    if translator.errors:
        sys.exit(1)
    if args.strict and translator.warnings:
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
