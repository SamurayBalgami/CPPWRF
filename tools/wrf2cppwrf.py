#!/usr/bin/env python3
"""WRF namelist.input to CPPWRF configuration (.ini) translator.

Converts WRF namelist.input into CPPWRF-compatible INI configuration files.
Implements a pure Python (f90nml-free) Fortran namelist parser.
Supports single-domain and nested multi-domain (max_dom > 1) execution.
Categorizes and warns/errors on unsupported or non-standard WRF options.

"""

import argparse
import datetime
import math
import os
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wrf_max_domains import MAX_DOMAINS_DEFAULTS  # noqa: E402
from wrf_namelist_options import WRF_NAMELIST_OPTIONS  # noqa: E402


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

# physics options a nest may set differently from its parent (verified bit-identical with WRF nests)
PER_DOMAIN_PHYS = {"cu_physics", "cudt", "radt", "bldt"}  # max_domains entries of the WRF Registry
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
        5: "MYNN surface layer",
    },
    "sf_surface_physics": {
        0: "None / off",
        2: "Unified Noah land surface model (4 soil layers)",
        4: "Noah-MP land surface model (v5.1.2, 4 soil layers)",
    },
    "bl_pbl_physics": {
        0: "None / off",
        1: "YSU (Yonsei University) planetary boundary layer",
        2: "MYJ (Mellor-Yamada-Janjic) planetary boundary layer",
        5: "MYNN-EDMF planetary boundary layer (WRF's default bl_mynn_* options)",
    },
    "cu_physics": {
        0: "None / off",
        1: "Kain-Fritsch cumulus parameterization",
        3: "Grell-Freitas cumulus parameterization",
        6: "Tiedtke cumulus parameterization",
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
        3: "3D Smagorinsky",
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

# share/module_check_a_mundo.F setup_physics_suite: the -1 options of each suite
PHYSICS_SUITES = {
    "conus": {"mp_physics": 8, "cu_physics": 6, "ra_lw_physics": 4, "ra_sw_physics": 4,
              "bl_pbl_physics": 2, "sf_sfclay_physics": 2, "sf_surface_physics": 2},
    "tropical": {"mp_physics": 6, "cu_physics": 16, "ra_lw_physics": 4, "ra_sw_physics": 4,
                 "bl_pbl_physics": 1, "sf_sfclay_physics": 1, "sf_surface_physics": 2},
}

SUPPORTED_NOAHMP_OPTIONS = {
    "dveg": (1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    "opt_crs": (1, 2),
    "opt_btr": (1, 2, 3),
    "opt_runsrf": (1, 2, 3, 4),
    "opt_runsub": (1, 2, 3, 4),
    "opt_sfc": (1, 2),
    "opt_frz": (1, 2),
    "opt_inf": (1, 2),
    "opt_rad": (1, 2, 3),
    "opt_alb": (1, 2),
    "opt_snf": (1, 2, 3, 4),
    "opt_tbot": (1, 2),
    "opt_stc": (1, 2, 3),
    "opt_gla": (1, 2),
    "opt_rsf": (1, 2, 3, 4),
    "opt_soil": (1,),
    "opt_pedo": (1,),
    "opt_crop": (0,),
    "opt_irr": (0,),
    "opt_irrm": (0,),
    "opt_infdv": (0,),
    "opt_tdrn": (0,),
    "opt_tksno": (1, 2, 3, 4, 5),
    "opt_compact": (1, 2),
    "opt_scf": (1, 2, 3, 4, 5),
    "opt_wetland": (0, 1, 2),
}

# &dynamics options that WRF applies but CPPWRF does not port: any non-default value is an error
UNSUPPORTED_DYNAMICS_SWITCHES = {
    "diff_6th_slopeopt": "6th-order diffusion slope limiter",
    "sfs_opt": "nonlinear backscatter / anisotropy (NBA) subfilter stresses",
    "m_opt": "subfilter stress output",
    "tracer_mix2_off": "tracer 2nd-order diffusion switch",
    "tracer_mix6_off": "tracer 6th-order diffusion switch",
    "chem_mix2_off": "chemistry 2nd-order diffusion switch",
    "chem_mix6_off": "chemistry 6th-order diffusion switch",
    "zadvect_implicit": "implicit-explicit vertical advection (IEVA)",
    "pert_coriolis": "perturbation Coriolis",
    "tracer_opt": "passive tracers",
}

# dynamics options whose value must stay at the WRF default (other values are not ported)
DEFAULT_ONLY_OPTIONS = {
    ("dynamics", "momentum_adv_opt"): (1, "WENO momentum advection"),
    ("dynamics", "phi_adv_z"): (1, "geopotential advection without double-staggered omega"),
    ("domains", "use_adaptive_time_step"): (False, "adaptive time step"),
    ("physics", "no_mp_heating"): (0, "microphysics without latent heating"),
    ("physics", "mp_zero_out"): (0, "microphysics zero-out of small moist values"),
    ("bdy_control", "polar"): (False, "polar filters and boundaries"),
    ("bdy_control", "have_bcs_moist"): (False, "moist species from the boundary file"),
    ("bdy_control", "have_bcs_scalar"): (False, "scalar species from the boundary file"),
    ("dynamics", "use_baseparam_fr_nml"): (False, "base-state parameters from the namelist"),
    ("dynamics", "use_input_w"): (False, "vertical velocity from the input file"),
    ("domains", "rebalance"): (0, "rebalancing of the input state"),
    ("time_control", "cycling"): (False, "cycling of the input state"),
    ("dfi_control", "dfi_opt"): (0, "digital filter initialisation"),
    ("domains", "allowed"): (True, "domains that are not integrated (grid_allowed)"),
    # physics and surface options of WRF's drivers that CPPWRF does not port (they change results or the output)
    ("physics", "ua_phys"): (False, "UA Noah snow physics"),
    ("physics", "sf_surface_mosaic"): (0, "Noah mosaic land use"),
    ("physics", "tmn_update"): (0, "deep soil temperature update"),
    ("physics", "sst_skin"): (0, "skin SST"),
    ("physics", "sf_ocean_physics"): (0, "ocean mixed layer / 3-D ocean model"),
    ("physics", "omlcall"): (0, "ocean mixed layer"),
    ("physics", "bucket_mm"): (-1, "rain accumulation buckets (I_RAINC, I_RAINNC)"),
    ("physics", "bucket_j"): (-1, "radiation accumulation buckets (I_ACSWUPT, ...)"),
    ("physics", "prec_acc_dt"): (0, "precipitation accumulation outputs (PREC_ACC_C, PREC_ACC_NC, SNOW_ACC_NC)"),
    ("physics", "aercu_opt"): (0, "aerosol-aware cumulus"),
    ("physics", "windfarm_opt"): (0, "wind farm parameterisation"),
    ("physics", "sf_surf_irr_scheme"): (0, "irrigation"),
    ("physics", "distributed_ahe_opt"): (0, "anthropogenic heat (AHE of the urban canopy model)"),
    ("physics", "do_radar_ref"): (0, "radar reflectivity diagnostic (REFL_10CM)"),
    ("physics", "kf_edrates"): (0, "KF entrainment / detrainment rate outputs (UDR_KF, DDR_KF, UER_KF, DER_KF, TIMEC_KF)"),
    ("physics", "ra_sw_eclipse"): (0, "solar eclipse in the shortwave radiation (ECOBSC, ECMASK)"),
    ("diags", "solar_diagnostics"): (0, "solar diagnostics (radiation_driver's Solar_step every step, SZA, CLRNIDX ...)"),
    ("fdda", "grid_fdda"): (0, "grid nudging"),
    ("fdda", "grid_sfdda"): (0, "surface grid nudging"),
    ("fdda", "obs_nudge_opt"): (0, "observation nudging"),
    ("fire", "ifire"): (0, "wildland fire model"),
    ("stoch", "multi_perturb"): (0, "stochastic multi-physics perturbations"),
    ("stoch", "sppt"): (0, "stochastically perturbed physics tendencies"),
    ("stoch", "skebs"): (0, "stochastic kinetic-energy backscatter"),
    ("stoch", "rand_perturb"): (0, "random perturbations"),
    ("stoch", "spp"): (0, "stochastically perturbed parameters"),
    ("stoch", "spp_pbl"): (0, "stochastically perturbed PBL parameters"),
    ("stoch", "spp_conv"): (0, "stochastically perturbed convection parameters"),
    ("stoch", "spp_lsm"): (0, "stochastically perturbed land-surface parameters"),
    ("stoch", "perturb_bdy"): (0, "perturbed lateral boundary tendencies (spec_bdy_dry_perturb)"),
    ("stoch", "perturb_chem_bdy"): (0, "perturbed chemistry boundaries"),
    ("bdy_control", "constant_bc"): (False, "constant lateral boundaries (zero_bdytend)"),
    ("bdy_control", "multi_bdy_files"): (False, "one boundary file per interval"),
}

# max_domains options that are single values in CPPWRF: every domain must have domain 1's value
SAME_ON_ALL_DOMAINS = (("physics", "kim_tofd"), ("physics", "tofd_factor"), ("dynamics", "gwd_dx_factor"),
                       ("dynamics", "gwd_if_nonhyd"))

UNSUPPORTED_PHYSICS_SWITCHES = {
    "cu_rad_feedback": "Cumulus-radiation feedback",
    "shcu_physics": "Separate shallow cumulus scheme",
    "sf_lake_physics": "Lake physics model",
    "topo_wind": "Topographic wind correction",
    "slope_rad": "Slope-dependent radiation",
    "topo_shading": "Topographic shading",
    "swint_opt": "Shortwave radiation interpolation between radt calls",
    "seaice_thickness_opt": "Variable sea-ice thickness option",
    "seaice_snowdepth_opt": "Sea-ice snow depth option",
    "seaice_albedo_opt": "Sea-ice albedo option",
    "grav_settling": "Fog/droplet gravitational settling",
    "ra_call_offset": "Radiation call offset",
    "isftcflx": "Alternative tropical cyclone surface flux formulation",
    "iz0tlnd": "Alternative thermal roughness length over land",
    "shalwater_z0": "Shallow-water roughness formulation",
    "scm_force_flux": "Single-column model prescribed surface flux",
    "rdlai2d": "2D LAI read from input (Noah)",
    "usemonalb": "Monthly surface albedo from geogrid",
}

IGNORED_WRF_OPTIONS = {
    # &time_control
    "frames_per_outfile": "Frame count per NetCDF file handled by output format",
    "debug_level": "WRF internal logging verbosity",
    "input_from_file": "Real case initialization always reads wrfinput",
    "iofields_filename": "Custom NetCDF field selection (handled by WrfHistory via output.wrf_namelist)",
    "io_form_history": "History output format specified in [output] format",
    "io_form_restart": "No restart files (restart_interval = 0)",
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
    "grid_id": "Domain grid hierarchy identifier",
    "parent_id": "Domain nesting parent identifier",
    "i_parent_start": "Domain nesting starting index",
    "j_parent_start": "Domain nesting starting index",
    "parent_grid_ratio": "Nesting grid refinement ratio",
    "parent_time_step_ratio": "Nesting timestep ratio",
    "feedback": "Two-way nesting feedback",
    "smooth_option": "Nesting boundary smoothing",
    "num_metgrid_levels": "Metgrid vertical levels used by real.exe",
    "num_metgrid_soil_levels": "Metgrid soil levels used by real.exe",
    "e_we": "Real case grid dimension read directly from wrfinput",
    "e_sn": "Real case grid dimension read directly from wrfinput",
    "e_vert": "Real case grid dimension read directly from wrfinput",
    "s_we": "Domain starting index in X (always 1)",
    "s_sn": "Domain starting index in Y (always 1)",
    "s_vert": "Domain starting index in Z (always 1)",
    "p_top_requested": "Real case p_top read directly from wrfinput_d01",
    "dzbot": "real.exe vertical level generation (eta levels read from wrfinput_d01)",
    "dzstretch_s": "real.exe vertical level generation (eta levels read from wrfinput_d01)",
    "dzstretch_u": "real.exe vertical level generation (eta levels read from wrfinput_d01)",
    "num_land_cat": "Land-use category count (land use read from wrfinput_d01)",
    "io_form_auxinput4": "wrflowinp format (NetCDF); used with sst_update = 1",
    "auxinput4_inname": "wrflowinp file name; used with sst_update = 1",
    "auxinput4_interval": "wrflowinp interval; used with sst_update = 1",
    "dx": "Real case grid spacing read directly from wrfinput",
    "dy": "Real case grid spacing read directly from wrfinput",
    # &physics
    "ifsnow": "Snow effects flag (used by the slab LSM only)",
    "num_soil_layers": "4 soil layers (Noah, Noah-MP)",
    # &dynamics
    "base_temp": "Base potential temperature (read from sounding or wrfinput)",
    "hybrid_opt": "Real case hybrid coordinate flag read directly from wrfinput_d01 global attribute",
    "h_sca_adv_order": "Scalar horizontal advection order (coupled with h_mom_adv_order in CPPWRF)",
    "v_sca_adv_order": "Scalar vertical advection order (coupled with v_mom_adv_order in CPPWRF)",
    "epssm": "Acoustic divergence damping parameter (default 0.1)",
    "smdiv": "Acoustic divergence damping parameter (default 0.1)",
    "emdiv": "External mode filter coefficient (default 0.01)",
    "diff_6th_factor": "Sixth-order diffusion factor (default 0.12)",
    # &bdy_control
    "spec_bdy_width": "Real case lateral boundary width (must be 5; read from wrfbdy_d01)",
    "nested": "Child domain nested lateral boundary flag (handled automatically for [d02]..[d0N])",
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
        return repr(val)  # shortest text of the same double: WRF reads the namelist real in full
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
        self._max_dom: int = 1

    def warn(self, msg: str):
        self.warnings.append(msg)

    def error(self, msg: str):
        self.errors.append(msg)

    def check_supported(self, option: str, value: Any):
        """Error when a WRF option value has no CPPWRF implementation (SUPPORTED_OPTIONS)."""
        table = SUPPORTED_OPTIONS[option]
        if value not in table:
            listed = ", ".join(f"{k} ({v})" for k, v in table.items() if k != 0)
            self.error(f"WRF {option} = {value} is not supported; CPPWRF supports 0, {listed}.")

    def notice(self, msg: str):
        self.notices.append(msg)

    def get_opt_list(self, group: str, key: str) -> List[Any]:
        grp = self.nml.get(group, {})
        self.used_keys.add((group, key))
        return grp.get(key, [])

    def get_dom_opt(self, group: str, key: str, dom_idx: int, default: Any = None,
                    fallback_d0: bool = True) -> Any:
        vals = self.get_opt_list(group, key)
        if not vals:
            return default
        if dom_idx < len(vals):
            return vals[dom_idx]
        if fallback_d0:
            return vals[0]
        return default

    def get_opt(self, group: str, key: str, default: Any = None) -> Any:
        grp = self.nml.get(group, {})
        self.used_keys.add((group, key))
        vals = grp.get(key, [])
        if not vals:
            return default
        if len(vals) > 1 and self._max_dom == 1:
            self.warn(f"Option &{group}.{key} has {len(vals)} values {vals}; selecting domain 1 value ({vals[0]}).")
        return vals[0]

    def has_opt(self, group: str, key: str) -> bool:
        grp = self.nml.get(group, {})
        return key in grp and len(grp[key]) > 0

    def pad_max_domains(self):
        """A max_domains option with fewer entries than max_dom: WRF's namelist read leaves the missing domains at the
        Registry default (not domain 1's value); check_a_mundo copies domain 1's km_opt and diff_opt into -1 entries.
        Only the groups whose options change results (the output and file options are handled where they are read)."""
        for group in ("physics", "dynamics", "fdda", "stoch", "fire", "noah_mp"):
            for key, vals in self.nml.get(group, {}).items():
                if (group, key) not in MAX_DOMAINS_DEFAULTS or not 0 < len(vals) < self._max_dom:
                    continue
                default = MAX_DOMAINS_DEFAULTS[(group, key)]
                if key in ("km_opt", "diff_opt") and default == -1:
                    default = vals[0]
                vals.extend([default] * (self._max_dom - len(vals)))

    def apply_physics_suite(self):
        """WRF setup_physics_suite: the -1 physics options take the values of physics_suite (CONUS, tropical)."""
        suite = str(self.get_opt("physics", "physics_suite", "none")).strip().lower()
        phys = self.nml.get("physics", {})
        if suite == "none":
            for key in PHYSICS_SUITES["conus"]:
                if -1 in phys.get(key, []):
                    self.error(f"WRF {key} = -1 needs a physics_suite (CONUS or tropical).")
            return
        if suite not in PHYSICS_SUITES:
            self.error(f"WRF physics_suite = '{suite}' is not recognized (CONUS, tropical or none).")
            return
        for key, val in PHYSICS_SUITES[suite].items():
            vals = phys.get(key)
            if not vals:  # Registry default -1
                phys[key] = [val] * self._max_dom
            else:
                phys[key] = [val if v == -1 else v for v in vals]
        self.notice(f"physics_suite = '{suite}': " + ", ".join(f"{k} = {phys[k][0]}" for k in PHYSICS_SUITES[suite]))

    def translate(self) -> Dict[str, Dict[str, Any]]:
        """Perform the translation and return INI sections dict."""
        sections: Dict[str, Dict[str, Any]] = {
            "run": {},
        }

        # 1. Determine max_dom in &domains
        max_dom_raw = self.get_dom_opt("domains", "max_dom", 0, 1)
        if not isinstance(max_dom_raw, int) or max_dom_raw < 1:
            self.error(f"WRF max_dom = {max_dom_raw} is invalid; max_dom must be >= 1.")
            max_dom = 1
        else:
            max_dom = max_dom_raw
        self._max_dom = max_dom
        self.pad_max_domains()
        self.apply_physics_suite()

        # 2. Determine Case
        ideal_case = self.get_opt("ideal", "ideal_case", None)
        specified = self.get_opt("bdy_control", "specified", False)

        case_name = "real"
        if self.override_case:
            case_name = self.override_case
        elif ideal_case is not None:
            # WRF ideal_case -> the ports of WRF's ideal.exe (src/cases/wrf_ideal.cpp)
            ideal_map = {
                0: "real",
                1: "em_hill2d_x",
                2: "em_quarter_ss",
                4: "em_squall2d_x",
                5: "em_squall2d_y",
                6: "em_grav2d_x",
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
                case_name = "em_hill2d_x"
            elif "qss" in dir_name or "quarter" in dir_name:
                case_name = "em_quarter_ss"
            elif "squall" in dir_name:
                case_name = "em_squall2d_y" if dir_name.endswith("_y") else "em_squall2d_x"
            elif "straka" in dir_name or "grav" in dir_name:
                case_name = "em_grav2d_x"
            else:  # a real case has specified = .true. (WRF's lateral boundaries need one condition per side)
                self.error("WRF case not identified: no ideal_case, specified = .false. and an unknown directory "
                           "name; pass --case.")
                case_name = None

        if case_name is None:
            case_name = "real"
        elif case_name == "real" and not specified:
            # domain 1 then has no condition on any side (open and periodic sides are errors below): wrf.exe stops
            # with "Error in boundary condition specification" (share/module_bc.F boundary_condition_check)
            self.error("WRF specified = .false. in a real case is not supported: WRF stops without a lateral "
                       "boundary condition on domain 1.")
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
                self.error(f"Could not compute run duration from start/end dates: {e}")
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
        # the column physics suite (surface layer, LSM, PBL, cumulus, radiation) runs in real-data cases only
        if case_name != "real" and has_physics_suite:
            self.error("WRF sf_sfclay_physics / sf_surface_physics / bl_pbl_physics / cu_physics / ra_lw_physics / "
                       "ra_sw_physics in an ideal case are not supported (CPPWRF runs them in real-data cases only).")

        # WRF tiles (frame/module_tiles.F set_tiles2) that split the i rows change WRF's results: row-coupled schemes
        # (Tiedtke's icum test) and the libmvec lanes of vectorised loops run over its..ite of a tile. CPPWRF gives
        # WRF's single-tile result, so tilings with more than one tile in x are errors; tiles in y are equivalent.
        numtiles = self.get_opt("domains", "numtiles", 1)
        numtiles_inc = max(self.get_opt("domains", "numtiles_inc", 0), 1)
        tile_sz_x = self.get_opt("domains", "tile_sz_x", 0)
        tile_sz_y = self.get_opt("domains", "tile_sz_y", 0)
        tile_strategy = self.get_opt("domains", "tile_strategy", 0)
        for d in range(self._max_dom):
            e_we = self.get_dom_opt("domains", "e_we", d, 0)
            e_sn = self.get_dom_opt("domains", "e_sn", d, 0)
            nt, ntx, nty = numtiles, self.get_opt("domains", "numtiles_x", 0), self.get_opt("domains", "numtiles_y", 0)
            strategy = tile_strategy
            if strategy == 0:  # TILE_NONE: 1-D in y unless there are more tiles than rows
                strategy = 2
                if e_sn > 0 and nt > e_sn and ntx == 0 and nty == 0:
                    ntx, nty = 1, e_sn
                    while ntx * numtiles_inc * nty < nt:
                        ntx += 1
                    ntx, strategy = ntx * numtiles_inc, 3
            if tile_sz_x >= 1 and tile_sz_y >= 1:
                ntx = -(-e_we // tile_sz_x) if e_we > 0 else 2
            elif ntx >= 1 or nty >= 1:
                if ntx < 1:
                    ntx = nt // nty
            elif strategy == 1:  # TILE_X
                ntx = nt
            elif strategy == 3:  # TILE_XY: least_aspect
                ntx, mini = nt, 2 * nt
                for y in range(1, nt + 1):
                    if nt % y == 0 and abs(y - nt // y) < mini:
                        mini, ntx = abs(y - nt // y), nt // y
            else:
                ntx = 1
            if ntx > 1:
                self.error(f"WRF tiling with {ntx} tiles in x on domain {d + 1} (numtiles / numtiles_x / tile_sz_x / "
                           "tile_strategy) is not supported: WRF's results depend on the tiles' i ranges, CPPWRF "
                           "matches WRF with one tile in x.")
                break

        # 5. Real Data Section [real]
        if case_name == "real":
            # domain 1 of a real case has specified boundaries on all sides; WRF's channel setup (specified +
            # periodic_x) is not ported, and open or symmetric sides with specified are fatal in WRF (module_bc)
            for key in ("periodic_x", "periodic_y", "open_xs", "open_xe", "open_ys", "open_ye", "symmetric_xs",
                        "symmetric_xe", "symmetric_ys", "symmetric_ye"):
                if self.get_opt("bdy_control", key, False):
                    self.error(f"WRF {key} = .true. in a real case is not supported (specified boundaries only).")
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

            # namelist values that wrf.exe checks against the input file (start_em.F, input_wrf.F: fatal
            # when they differ); the model makes the same checks
            sections["real"]["p_top_requested"] = self.get_opt("domains", "p_top_requested", 5000.0)
            sections["real"]["iso_temp"] = self.get_opt("dynamics", "iso_temp", 200.0)
            for k in ("dx", "dy"):
                if self.has_opt("domains", k):
                    sections["real"][k] = self.get_opt("domains", k, None)

            # In configs/real_p1..p4 and nested configs, wrfinput/wrfbdy are written explicitly
            if self.explicit_files or has_physics_suite or max_dom > 1:
                sections["real"]["wrfinput"] = "wrfinput_d01"
                sections["real"]["wrfbdy"] = "wrfbdy_d01"

        # 6. Domain Dimensions & Boundary Conditions [domain] (ideal cases)
        if case_name != "real":
            sections["domain"] = {}
            em = case_name.startswith("em_")
            # missing keys: WRF's Registry defaults for the em_* cases (e_we = e_sn = 32, e_vert = 31,
            # dx = dy = 200 m, ztop = 15 km); the analytic cases keep CPPWRF's own
            e_we = self.get_opt("domains", "e_we", 32 if em else 100)
            e_sn = self.get_opt("domains", "e_sn", 32 if em else 3)
            e_vert = self.get_opt("domains", "e_vert", 31 if em else 41)
            s_we = self.get_opt("domains", "s_we", 1)
            s_sn = self.get_opt("domains", "s_sn", 1)
            s_vert = self.get_opt("domains", "s_vert", 1)

            nx = e_we - s_we
            nz = e_vert - s_vert
            ny = 1 if (e_sn - s_sn <= 2) else (e_sn - s_sn)

            sections["domain"]["nx"] = nx
            sections["domain"]["ny"] = ny
            sections["domain"]["nz"] = nz

            dx = self.get_opt("domains", "dx", 200.0 if em else 1000.0)
            dy = self.get_opt("domains", "dy", 200.0 if em else dx)
            sections["domain"]["dx"] = dx
            if ny > 1 or em:
                sections["domain"]["dy"] = dy

            ztop = self.get_opt("domains", "ztop", 15000.0 if em else None)
            if ztop is not None:
                sections["domain"]["ztop"] = ztop

            p_top = self.get_opt("domains", "p_top_requested", self.get_opt("domains", "p_top", None))
            if p_top is not None:
                sections["domain"]["p_top"] = p_top

            def lateral_bc(d):
                """share/module_bc.F: exactly one condition per side (all .false. by default, else fatal);
                CPPWRF has one condition per direction (periodic or open)."""
                if not em:  # CPPWRF's analytic cases: periodic unless set
                    if self.get_opt("bdy_control", f"periodic_{d}", True):
                        return "periodic"
                    if self.get_opt("bdy_control", f"symmetric_{d}s", False):
                        return "symmetric"
                    return "open" if self.get_opt("bdy_control", f"open_{d}s", False) else None
                kinds = []
                for side in ("s", "e"):
                    on = [k for k, key in (("periodic", f"periodic_{d}"), ("symmetric", f"symmetric_{d}{side}"),
                                           ("open", f"open_{d}{side}"), ("specified", "specified"),
                                           ("nested", "nested"))
                          if self.get_opt("bdy_control", key, False)]
                    if len(on) != 1:
                        self.error(f"WRF needs exactly one lateral boundary condition on the {d}{side} side "
                                   f"(got {', '.join(on) or 'none'}).")
                        return None
                    kinds.append(on[0])
                if kinds[0] != kinds[1]:
                    self.error(f"WRF lateral boundaries {kinds[0]}_{d}s / {kinds[1]}_{d}e are not supported "
                               f"(CPPWRF has one condition per direction).")
                    return None
                if kinds[0] in ("specified", "nested"):
                    self.error(f"WRF {kinds[0]} lateral boundaries are not supported in ideal cases.")
                    return None
                return kinds[0]

            bcx = lateral_bc("x")
            if bcx is not None:
                sections["domain"]["bc_x"] = bcx  # open: BC::Open, symmetric: BC::Symmetric (src/cases/model.cpp)

            bcy = lateral_bc("y")
            if (ny > 1 or "bc_x" in sections["domain"]) and bcy is not None:
                sections["domain"]["bc_y"] = bcy

        # 7. Ideal Case Def Section [case]
        if case_name.startswith("em_"):
            # ideal.exe reads input_sounding from the run directory, i.e. next to namelist.input
            nml_dir = os.path.dirname(self.namelist_path) if self.namelist_path else "."
            snd = os.path.join(os.path.abspath(nml_dir), "input_sounding")
            sections["case"] = {"sounding_file": snd}
            # symmetric walls need zero normal wind there: with a wind across them WRF keeps the wall value of
            # the normal velocity while muu changes sign at the wall, and its solution breaks down (the model
            # checks the interpolated wall winds; here any nonzero u / v of the sounding)
            dom = sections.get("domain", {})
            cols = [c for c, bc in ((3, dom.get("bc_x")), (4, dom.get("bc_y") if dom.get("ny", 1) > 1 else None))
                    if bc == "symmetric"]
            if cols and os.path.exists(snd):
                rows = [l.split() for l in open(snd).read().splitlines()[1:]]
                for c in cols:
                    if any(len(r) > c and float(r[c]) != 0.0 for r in rows):
                        self.error(f"WRF symmetric_{'x' if c == 3 else 'y'}s/{'x' if c == 3 else 'y'}e with a "
                                   f"nonzero input_sounding {'u' if c == 3 else 'v'}: a wind across the walls is "
                                   f"not supported (WRF's solution breaks down).")
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
            self.error("WRF time_step_fract_num/den (fractional time step) is not supported: WRF's boundary "
                       "time dtbc is a running sum of dt, not verified for a non-integer dt.")
            dt += float(num) / float(den)
        sections["time"]["dt"] = int(dt) if dt.is_integer() else dt
        # start_em.F: d01 of a real-data run (map_proj /= 0) with dt / dx [s/km] above reasonable_time_step_ratio
        # stops with "Time step too large"
        if case_name == "real":
            dx_km = min(float(self.get_opt("domains", "dx", 200.0)), float(self.get_opt("domains", "dy", 200.0))) / 1000.
            ratio = float(self.get_opt("domains", "reasonable_time_step_ratio", 6.0))
            if dt / dx_km > ratio:
                self.error(f"WRF stops with 'Time step too large': time_step / dx = {dt / dx_km:g} s/km is above "
                           f"reasonable_time_step_ratio = {ratio:g}.")

        tss =self.get_opt("dynamics", "time_step_sound", 0)
        if tss and tss > 0:
            sections["time"]["time_step_sound"] = tss

        # 9. Dynamics Section [dyn]
        sections["dyn"] = {}
        khdif = self.get_opt("dynamics", "khdif", 0.0)
        kvdif = self.get_opt("dynamics", "kvdif", 0.0)
        sections["dyn"]["khdif"] = int(khdif) if isinstance(khdif, (int, float)) and float(khdif).is_integer() else khdif
        sections["dyn"]["kvdif"] = int(kvdif) if isinstance(kvdif, (int, float)) and float(kvdif).is_integer() else kvdif
        if self.get_opt("dynamics", "gwd_opt", 0) == 1:  # GWDO in the PBL driver (needs CON, OA1-4, OL1-4, ELVMAX in wrfinput)
            sections["dyn"]["gwd_opt"] = 1
            gdx = self.get_opt("dynamics", "gwd_dx_factor", 2.0)
            if gdx != 2.0:
                sections["dyn"]["gwd_dx_factor"] = gdx
            if not self.get_opt("dynamics", "gwd_if_nonhyd", True):
                sections["dyn"]["gwd_if_nonhyd"] = False

        w_damping = self.get_opt("dynamics", "w_damping", 0)
        if w_damping not in (0, 1):
            self.error(f"WRF w_damping = {w_damping} is not supported; CPPWRF supports 0 and 1.")
        elif w_damping != 0:
            sections["dyn"]["w_damping"] = w_damping

        if self.has_opt("dynamics", "w_crit_cfl"):
            w_crit_cfl = self.get_opt("dynamics", "w_crit_cfl", 1.0)
            if float(w_crit_cfl) <= 0.0:
                self.error(f"WRF w_crit_cfl = {w_crit_cfl} must be > 0.")
            elif float(w_crit_cfl) != 1.0:
                sections["dyn"]["w_crit_cfl"] = w_crit_cfl

        # WRF v4.8.0 Registry default damp_opt = 3 (w Rayleigh damping, zdamp 5000, dampcoef 0.2); the
        # model's default is 0, so a namelist without damp_opt must still produce damp_opt = 3
        damp_opt = self.get_opt("dynamics", "damp_opt", 3)
        if damp_opt != 0:
            sections["dyn"]["damp_opt"] = damp_opt
            if damp_opt not in SUPPORTED_OPTIONS["damp_opt"]:
                self.error(f"WRF damp_opt = {damp_opt} is not supported; CPPWRF supports damp_opt 0, 1, 2 and 3.")
            zdamp = self.get_opt("dynamics", "zdamp", 5000.0)
            dampcoef = self.get_opt("dynamics", "dampcoef", 0.2)
            sections["dyn"]["zdamp"] = int(zdamp) if float(zdamp).is_integer() else zdamp
            ztop_nl = self.get_opt("domains", "ztop", None)
            if ztop_nl is not None and float(zdamp) >= float(ztop_nl):
                # WRF's damping layer then reaches below the model bottom (rk_rayleigh_damp reads out of
                # bounds; wrf.exe crashes)
                self.error(f"WRF zdamp = {zdamp} >= ztop = {ztop_nl} is not supported (WRF fails).")
            sections["dyn"]["dampcoef"] = dampcoef
        elif case_name in ("hill2d",):
            sections["dyn"]["damp_opt"] = 0

        # WRF Registry default diff_opt = km_opt = -1: wrf.exe stops ("Both km_opt and diff_opt need to
        # be set in the namelist.input file", share/module_check_a_mundo.F)
        for opt_key in ("diff_opt", "km_opt"):
            if not self.has_opt("dynamics", opt_key):
                self.notice(f"WRF {opt_key} is not set: wrf.exe v4.8.0 refuses this namelist (both km_opt and "
                            "diff_opt must be set); CPPWRF assumes diff_opt = 0, km_opt = 1.")
        diff_opt = self.get_opt("dynamics", "diff_opt", 0)
        if diff_opt == 0 and (float(khdif) != 0.0 or float(kvdif) != 0.0):
            sections["dyn"]["diff_opt"] = 0  # WRF: no diffusion; CPPWRF's default with khdif/kvdif is 1
        if diff_opt != 0:
            sections["dyn"]["diff_opt"] = diff_opt
            if diff_opt not in SUPPORTED_OPTIONS["diff_opt"]:
                self.error(f"WRF diff_opt = {diff_opt} is not supported; CPPWRF supports diff_opt 0, 1 and 2.")

        km_opt = self.get_opt("dynamics", "km_opt", 1)
        if (diff_opt != 0 or km_opt != 1) and km_opt is not None:
            if km_opt not in (0, 1):
                sections["dyn"]["km_opt"] = km_opt
            if km_opt not in SUPPORTED_OPTIONS["km_opt"] and km_opt != 0:
                self.error(f"WRF km_opt = {km_opt} is not supported; CPPWRF supports km_opt "
                           + ", ".join(str(k) for k in sorted(SUPPORTED_OPTIONS["km_opt"])) + ".")

        c_s = self.get_opt("dynamics", "c_s", None)
        if c_s is not None and (km_opt in (3, 4) or float(c_s) != 0.25):
            sections["dyn"]["c_s"] = c_s

        c_k = self.get_opt("dynamics", "c_k", None)
        if c_k is not None and (km_opt == 2 or float(c_k) != 0.15):
            sections["dyn"]["c_k"] = c_k

        for key in ("mix_isotropic", "mix_upper_bound"):
            if self.has_opt("dynamics", key):
                val = self.get_opt("dynamics", key, None)
                if key == "mix_isotropic" and val not in (0, 1):
                    self.error(f"WRF mix_isotropic = {val} is not supported; CPPWRF supports 0 and 1.")
                elif key == "mix_upper_bound" and float(val) <= 0.0:
                    self.error(f"WRF mix_upper_bound = {val} must be > 0.")
                elif diff_opt != 0 or (key == "mix_isotropic" and val != 0) or (key == "mix_upper_bound" and float(val) != 0.1):
                    sections["dyn"][key] = val

        if self.has_opt("dynamics", "mix_full_fields"):
            mix_ff = self.get_opt("dynamics", "mix_full_fields", False)
            if bool(mix_ff):
                sections["dyn"]["mix_full_fields"] = True

        for tke_key, def_v in (("tke_drag_coefficient", 0.0), ("tke_heat_flux", 0.0), ("tke_upper_bound", 1000.0)):
            if self.has_opt("dynamics", tke_key):
                val = self.get_opt("dynamics", tke_key, def_v)
                if tke_key == "tke_upper_bound" and float(val) <= 0.0:
                    self.error(f"WRF tke_upper_bound = {val} must be > 0.")
                elif float(val) != def_v:
                    sections["dyn"][tke_key] = val

        for mix_off_key in ("moist_mix2_off", "scalar_mix2_off", "tke_mix2_off",
                            "moist_mix6_off", "scalar_mix6_off", "tke_mix6_off"):
            if self.has_opt("dynamics", mix_off_key):
                val = self.get_opt("dynamics", mix_off_key, False)
                if bool(val):
                    sections["dyn"][mix_off_key] = True

        diff_6th_opt = self.get_opt("dynamics", "diff_6th_opt", 0)
        if diff_6th_opt not in (0, 1, 2):
            self.error(f"WRF diff_6th_opt = {diff_6th_opt} is not supported; CPPWRF supports 0, 1, and 2.")
        elif diff_6th_opt != 0:
            sections["dyn"]["diff_6th_opt"] = diff_6th_opt
            diff_6th_factor = self.get_opt("dynamics", "diff_6th_factor", 0.12)
            if float(diff_6th_factor) <= 0.0:
                self.error(f"WRF diff_6th_factor = {diff_6th_factor} must be > 0.")
            sections["dyn"]["diff_6th_factor"] = diff_6th_factor

        moist_adv_opt = self.get_opt("dynamics", "moist_adv_opt", 1)
        if moist_adv_opt not in (0, 1):
            self.error(f"WRF moist_adv_opt = {moist_adv_opt} is not supported; CPPWRF supports 0 (simple) and 1 (positive-definite).")
        elif moist_adv_opt is not None and case_name != "hill2d":
            sections["dyn"]["moist_adv_opt"] = moist_adv_opt

        # scalar_adv_opt: included when physics suite is enabled, when explicitly 0, or when != 0 for non-dry ideal runs
        scalar_adv_opt = self.get_opt("dynamics", "scalar_adv_opt", 0)
        if scalar_adv_opt not in (0, 1):
            self.error(f"WRF scalar_adv_opt = {scalar_adv_opt} is not supported; CPPWRF supports 0 (simple) and 1 (positive-definite).")
        elif self.has_opt("dynamics", "scalar_adv_opt") and scalar_adv_opt == 0:
            sections["dyn"]["scalar_adv_opt"] = 0
        elif has_physics_suite or (scalar_adv_opt != 0 and case_name in ("quarter_ss",)):
            sections["dyn"]["scalar_adv_opt"] = scalar_adv_opt if scalar_adv_opt != 0 else 1

        # WRF's prandtl is the constant 1./3. of module_model_constants.F, not a namelist option (an
        # &dynamics prandtl is rejected with the other options WRF does not have, below)
        if case_name in ("straka", "quarter_ss"):
            sections["dyn"]["pr_inv"] = 3.0

        top_lid = self.get_opt("dynamics", "top_lid", False)
        if top_lid or case_name in ("quarter_ss",):
            sections["dyn"]["top_lid"] = top_lid

        h_adv = self.get_opt("dynamics", "h_mom_adv_order", 5)
        if h_adv != 5 or case_name in ("quarter_ss",):
            sections["dyn"]["h_adv_order"] = h_adv
            if h_adv not in SUPPORTED_OPTIONS["h_adv_order"]:
                self.error(f"WRF h_mom_adv_order = {h_adv} is not supported (5 or 6 only).")
        # Registry: h_sca_adv_order defaults to 5 whatever h_mom_adv_order is
        h_sca = self.get_opt("dynamics", "h_sca_adv_order", 5)
        if h_sca != h_adv:
            self.error(f"WRF h_sca_adv_order = {h_sca} differs from h_mom_adv_order = {h_adv}; "
                       "CPPWRF uses a single horizontal advection order (h_adv_order).")

        v_adv = self.get_opt("dynamics", "v_mom_adv_order", 3)
        if v_adv != 3 or case_name in ("quarter_ss",):
            sections["dyn"]["v_adv_order"] = v_adv
            if v_adv not in SUPPORTED_OPTIONS["v_adv_order"]:
                self.error(f"WRF v_mom_adv_order = {v_adv} is not supported (3 or 4 only).")
        # Registry: v_sca_adv_order defaults to 3 whatever v_mom_adv_order is
        v_sca = self.get_opt("dynamics", "v_sca_adv_order", 3)
        if v_sca != v_adv:
            self.error(f"WRF v_sca_adv_order = {v_sca} differs from v_mom_adv_order = {v_adv}; "
                       "CPPWRF uses a single vertical advection order (v_adv_order).")

        for ac_key, def_v in (("epssm", 0.1), ("smdiv", 0.1), ("emdiv", 0.01)):
            if self.has_opt("dynamics", ac_key):
                val = self.get_opt("dynamics", ac_key, def_v)
                if float(val) != def_v:
                    sections["dyn"][ac_key] = val

        if self.has_opt("dynamics", "etac"):
            etac = self.get_opt("dynamics", "etac", 0.2)
            if float(etac) <= 0.0 or float(etac) >= 1.0:
                self.error(f"WRF etac = {etac} must satisfy 0 < etac < 1.")
            elif float(etac) != 0.2:
                sections["dyn"]["etac"] = etac

        if case_name == "real":
            # Registry default 2; input_wrf.F stops when it differs from the input's HYBRID_OPT
            hybrid_opt = self.get_opt("dynamics", "hybrid_opt", 2)
            if hybrid_opt not in (0, 1, 2):
                self.error(f"WRF hybrid_opt = {hybrid_opt} is not supported; CPPWRF supports 0, 1, and 2.")
            sections["dyn"]["hybrid_opt"] = hybrid_opt
        elif self.has_opt("dynamics", "hybrid_opt"):
            hybrid_opt = self.get_opt("dynamics", "hybrid_opt", 0)
            if hybrid_opt not in (0, 1, 2):
                self.error(f"WRF hybrid_opt = {hybrid_opt} is not supported; CPPWRF supports 0, 1, and 2.")
            elif case_name != "real":
                # always written: the ideal cases default to hybrid_opt = 2 (WRF Registry), so a
                # namelist hybrid_opt = 0 must reach the model
                sections["dyn"]["hybrid_opt"] = hybrid_opt

        # hypsometric_opt is a &domains option (Registry.EM_COMMON, default 2). wrf.exe uses the namelist value
        # (the input's HYPSOMETRIC_OPT attribute only resets it to 1 when missing); the real case's default is 2
        if self.has_opt("domains", "hypsometric_opt"):
            hypso = self.get_opt("domains", "hypsometric_opt", 2)
            if hypso not in (1, 2):
                self.error(f"WRF hypsometric_opt = {hypso} is not supported; CPPWRF supports 1 and 2.")
            elif (case_name != "real" and hypso != 1) or (case_name == "real" and hypso != 2):
                sections["dyn"]["hypsometric_opt"] = hypso

        # WRF 4.8 options without a CPPWRF implementation
        if self.get_opt("time_control", "iau", 0):
            self.error("WRF iau (incremental analysis update) is not supported.")
        if self.get_opt("physics", "dyn_lightning_option", 0):
            self.error("WRF dyn_lightning_option is not supported.")

        # Validate core numerical assumptions
        non_hydro = self.get_opt("dynamics", "non_hydrostatic", True)
        if not non_hydro:
            self.error("WRF non_hydrostatic = .false. is unsupported; CPPWRF is strictly non-hydrostatic.")

        rk_ord = self.get_opt("dynamics", "rk_ord", 3)
        if rk_ord != 3:
            self.error(f"WRF rk_ord = {rk_ord} is unsupported; CPPWRF strictly uses 3rd-order Runge-Kutta (RK3).")

        use_thm = self.get_opt("dynamics", "use_theta_m", 1)
        if use_thm == 0:
            self.error("WRF use_theta_m = 0 is not supported; CPPWRF requires use_theta_m = 1.")

        gwd_opt = self.get_opt("dynamics", "gwd_opt", 0)
        if gwd_opt not in (0, 1):
            self.error(f"WRF gwd_opt = {gwd_opt} is not supported; CPPWRF supports 0 and 1 (GWDO).")
        for d in range(1, self._max_dom):  # module_check_a_mundo.F: a nest has domain 1's gwd_opt or 0
            gwd_d = self.get_dom_opt("dynamics", "gwd_opt", d, 0)
            if gwd_d != gwd_opt and gwd_d != 0:
                self.error(f"WRF gwd_opt must be equal for all domains (or = zero): domain {d + 1} has gwd_opt = "
                           f"{gwd_d}, domain 1 {gwd_opt}.")

        chem_adv = self.get_opt("dynamics", "chem_adv_opt", 0)
        if chem_adv != 0:
            self.warn("WRF chemistry advection (chem_adv_opt) is not supported.")

        tke_adv = self.get_opt("dynamics", "tke_adv_opt", None)
        if tke_adv is not None:
            if tke_adv not in (0, 1):
                self.error(f"WRF tke_adv_opt = {tke_adv} is not supported; CPPWRF supports 0 and 1.")
            elif km_opt == 2 or tke_adv != 1:
                sections["dyn"]["tke_adv_opt"] = tke_adv

        # 10. Moisture Species [moist]
        if case_name == "real":
            sections["moist"] = {}
            if mp_physics == 8:  # Thompson: two-moment ice and rain (QNICE, QNRAIN)
                sections["moist"]["species"] = "qv,qc,qr,qi,qs,qg,qni,qnr"
            elif mp_physics == 6:
                sections["moist"]["species"] = "qv, qc, qr, qi, qs, qg" if not has_physics_suite else "qv,qc,qr,qi,qs,qg"
            elif mp_physics == 1:
                sections["moist"]["species"] = "qv, qc, qr"
            else:
                sections["moist"]["species"] = "qv"
        elif mp_physics > 0 or case_name in ("quarter_ss", "squall2d"):
            sections["moist"] = {}
            if mp_physics == 8:
                sections["moist"]["species"] = "qv,qc,qr,qi,qs,qg,qni,qnr"
            elif mp_physics == 6:
                sections["moist"]["species"] = "qv, qc, qr, qi, qs, qg"
            else:
                sections["moist"]["species"] = "qv, qc, qr"

        # 11. Physics Section [phys]
        if self.get_opt("physics", "sst_update", 0) != 0 and not (mp_physics != 0 or has_physics_suite):
            self.error("WRF sst_update = 1 needs the physics suite (mp_physics /= 0) in CPPWRF.")
        if mp_physics != 0 or has_physics_suite:
            sections["phys"] = {}
            if mp_physics != 0:
                sections["phys"]["mp_physics"] = mp_physics
                self.check_supported("mp_physics", mp_physics)

            if has_physics_suite:
                sections["phys"]["sf_sfclay_physics"] = sf_sfclay
                sections["phys"]["sf_surface_physics"] = sf_surface
                sections["phys"]["bl_pbl_physics"] = bl_pbl
                sections["phys"]["cu_physics"] = cu_physics
                sections["phys"]["ra_lw_physics"] = ra_lw
                sections["phys"]["ra_sw_physics"] = ra_sw

                if (ra_lw == 4 or ra_sw == 4) and (ra_lw != 4 or ra_sw != 4):
                    self.error(f"CPPWRF RRTMG requires both ra_lw_physics = 4 and ra_sw_physics = 4 "
                               f"(got ra_lw_physics = {ra_lw}, ra_sw_physics = {ra_sw}).")

                if ra_lw != 0 or ra_sw != 0 or radt != 0:
                    # as given: stepra = nint(radt*60/dt) and calc_coszen at xtime + radt/2 use it
                    sections["phys"]["radt"] = radt
                if bldt > 0:
                    sections["phys"]["bldt"] = bldt
                if cudt > 0:
                    sections["phys"]["cudt"] = cudt

                if sf_surface != 0:
                    # check_a_mundo uses 4 layers for Noah and Noah-MP, but the history output keeps the
                    # namelist's soil dimension (P1 with num_soil_layers = 5: same values, ZS/TSLB/SMOIS shapes differ)
                    num_soil = self.get_opt("physics", "num_soil_layers", 4)
                    if num_soil != 4:
                        self.error(f"WRF num_soil_layers = {num_soil} is not supported with "
                                   f"sf_surface_physics = {sf_surface}; CPPWRF requires 4 soil layers.")
                    else:
                        sections["phys"]["num_soil_layers"] = 4
                    opt_thcnd = self.get_opt("physics", "opt_thcnd", 1)
                    if opt_thcnd not in (1, 2):
                        self.error(f"WRF opt_thcnd = {opt_thcnd} is not supported; CPPWRF supports 1 and 2.")
                    elif opt_thcnd == 2 and sf_surface == 2:  # Noah only
                        sections["phys"]["opt_thcnd"] = 2

                tbl_dir = self.table_dir
                if not tbl_dir:
                    wrf_runs = os.environ.get("CPPWRF_WRF_RUNS")
                    # the WRF run directory first: wrf.exe reads its tables there (including the Thompson
                    # .dat tables, which CPPWRF otherwise computes)
                    candidates = []
                    if self.namelist_path:
                        candidates.append(os.path.dirname(os.path.abspath(self.namelist_path)))
                        candidates.append(os.path.join(os.path.dirname(self.namelist_path), "wrf_run"))
                    if wrf_runs:
                        candidates.append(os.path.join(wrf_runs, "wrf_run"))
                    req_tbls = ["LANDUSE.TBL"]
                    if sf_surface == 4:
                        req_tbls.append("NoahmpTable.TBL")
                    for cand in candidates:
                        if cand and all(os.path.exists(os.path.join(cand, t)) for t in req_tbls):
                            tbl_dir = cand
                            break
                    if not tbl_dir:
                        tbl_dir = os.path.join(wrf_runs, "wrf_run") if wrf_runs else "${CPPWRF_WRF_RUNS}/wrf_run"
                sections["phys"]["table_dir"] = tbl_dir

                # WRF 4.8 turbulent orographic form drag of sfclayrev (needs VAR in wrfinput)
                sstu = self.get_opt("physics", "sst_update", 0)
                if sstu not in (0, 1):
                    self.error(f"WRF sst_update = {sstu} is not supported; CPPWRF supports 0 and 1.")
                elif sstu == 1:
                    if self.get_opt("time_control", "io_form_auxinput4", 0) == 0:
                        self.error("WRF sst_update = 1 needs io_form_auxinput4 /= 0 (as WRF).")
                    for k in ("auxinput4_interval_s", "auxinput4_interval_h", "auxinput4_interval_d", "auxinput4_interval_y",
                              "auxinput4_begin_s", "auxinput4_begin_m", "auxinput4_begin_h", "auxinput4_begin_d",
                              "auxinput4_begin_y", "auxinput4_interval_m"):
                        if self.get_opt("time_control", k, 0) not in (0, None):  # only the minutes interval
                            self.error(f"WRF {k} is not supported; CPPWRF reads auxinput4_interval (minutes).")
                    iv = self.get_opt("time_control", "auxinput4_interval", 0)
                    if not iv:
                        self.error("WRF sst_update = 1 needs auxinput4_interval (minutes).")
                    sections["phys"]["sst_update"] = 1
                    sections["phys"]["auxinput4_interval"] = iv
                    inname = str(self.get_opt("time_control", "auxinput4_inname", "wrflowinp_d<domain>")).strip("'\"")
                    if inname != "wrflowinp_d<domain>":
                        sections["real"]["wrflowinp"] = inname.replace("<domain>", "01")
                # sea-ice parameters: CPPWRF runs with WRF's defaults only (no reference run with others)
                for si_key, si_def in (("seaice_albedo_default", 0.65), ("seaice_thickness_default", 3.0),
                                       ("seaice_snowdepth_max", 1.e10), ("seaice_snowdepth_min", 0.001)):
                    si_val = self.get_opt("physics", si_key, si_def)
                    if abs(float(si_val) - si_def) > 1e-12 * si_def:
                        self.error(f"WRF {si_key} = {si_val} is not supported; CPPWRF uses WRF's default {si_def}.")
                fsi = self.get_opt("physics", "fractional_seaice", 0)
                if fsi not in (0, 1):
                    self.error(f"WRF fractional_seaice = {fsi} is not supported; CPPWRF supports 0 and 1.")
                elif fsi == 1:
                    if sf_surface not in (2, 4) or sf_sfclay not in (1, 2):
                        self.error("WRF fractional_seaice = 1 needs sf_surface_physics 2 or 4 and sf_sfclay_physics 1 or 2.")
                    sections["phys"]["fractional_seaice"] = 1
                    if self.get_opt("physics", "tice2tsk_if2cold", False):
                        sections["phys"]["tice2tsk_if2cold"] = True
                if self.get_opt("physics", "kim_tofd", False):
                    sections["phys"]["kim_tofd"] = True
                    tofd_factor = self.get_opt("physics", "tofd_factor", 0.003)
                    if tofd_factor != 0.003:
                        sections["phys"]["tofd_factor"] = tofd_factor
                elif self.has_opt("physics", "tofd_factor"):
                    tofd_factor = self.get_opt("physics", "tofd_factor", 0.003)
                    if tofd_factor != 0.003:
                        sections["phys"]["tofd_factor"] = tofd_factor

                if bl_pbl == 5 or sf_sfclay == 5:
                    # MYNN: only WRF's default options are ported (Registry.EM_COMMON)
                    mynn_defaults = (("bl_mynn_tkeadvect", False), ("bl_mynn_cloudpdf", 2), ("bl_mynn_mixlength", 1),
                                     ("bl_mynn_edmf", 1), ("bl_mynn_edmf_dd", 0), ("bl_mynn_edmf_mom", 1),
                                     ("bl_mynn_edmf_tke", 0), ("bl_mynn_mixscalars", 0), ("bl_mynn_mixaerosols", 1),
                                     ("bl_mynn_mixnumcon", 0), ("bl_mynn_output", 0), ("bl_mynn_cloudmix", 1),
                                     ("bl_mynn_mixqt", 0), ("bl_mynn_ess", 1), ("tke_budget", 0),
                                     ("sf_mynn_sfcflux_land", 0), ("sf_mynn_sfcflux_water", 1), ("pert_mynn", False),
                                     ("spp_pbl", 0), ("mynn_chem_vertmx", False))
                    for mk, mdef in mynn_defaults:
                        if self.has_opt("physics", mk) and self.get_opt("physics", mk, mdef) != mdef:
                            self.error(f"WRF {mk} = {self.get_opt('physics', mk, mdef)} is not supported; "
                                       f"CPPWRF's MYNN uses WRF's default {mk} = {mdef}.")
                    closure = self.get_opt("physics", "bl_mynn_closure", 2.6)
                    if abs(float(closure) - 2.6) > 1e-9:
                        self.error(f"WRF bl_mynn_closure = {closure} is not supported; CPPWRF's MYNN uses 2.6.")
                    if bl_pbl == 5 and sf_sfclay != 5:
                        self.error("WRF bl_pbl_physics = 5 is supported with sf_sfclay_physics = 5 only.")
                    if bl_pbl == 5 and mp_physics not in (6, 8):
                        self.error("WRF bl_pbl_physics = 5 needs a microphysics with qc, qi and qs (mp_physics 6 or 8): "
                                   "WRF passes the MYNN-EDMF an uninitialised snow content otherwise.")
                    if bl_pbl == 5 and (self.get_opt("physics", "ishallow", 0) > 0 or self.get_opt("physics", "shcu_physics", 0) > 0):
                        self.error("WRF bl_mynn_edmf > 0 requires shcu_physics = 0 and ishallow = 0.")
                # module_physics_init.F pbl_select: YSU needs isfc = 1 (sfclay), MYJ isfc = 2 (myjsfc)
                if bl_pbl == 1 and sf_sfclay != 1:
                    self.error("WRF bl_pbl_physics = 1 (YSU) needs sf_sfclay_physics = 1 "
                               "(module_physics_init: use sfclay scheme for this pbl option).")
                if bl_pbl == 2 and sf_sfclay != 2:
                    self.error("WRF bl_pbl_physics = 2 (MYJ) needs sf_sfclay_physics = 2 "
                               "(module_physics_init: use myjsfc scheme for this pbl option).")
                if bl_pbl == 5:
                    icb = self.get_opt("physics", "icloud_bl", 1)
                    if icb not in (0, 1):
                        self.error(f"WRF icloud_bl = {icb} is not supported (0 or 1).")
                    elif icb != 1:
                        sections["phys"]["icloud_bl"] = icb
                if bl_pbl == 1 and sf_sfclay != 1:
                    self.error("WRF bl_pbl_physics = 1 (YSU) needs sf_sfclay_physics = 1 (module_physics_init: "
                               "'Use sf_sfclay_physics= 1 or 91 for this pbl option').")
                # module_pbl_driver: scalar_pblmix > 0 mixes the 4-D scalars (Thompson's QNI, QNR) with exch_h in
                # diff4d after the PBL scheme (not ported; Kessler and WSM6 have no scalars); MYNN turns it off (check_a_mundo)
                scalar_pblmix = self.get_opt_list("physics", "scalar_pblmix")
                if bl_pbl not in (0, 5) and mp_physics not in (0, 1, 6) and any(v not in (0, None) for v in scalar_pblmix):
                    self.error("WRF scalar_pblmix = 1 (PBL mixing of the number concentrations with diff4d) is not "
                               "supported by CPPWRF.")
                self.check_supported("sf_sfclay_physics", sf_sfclay)
                self.check_supported("sf_surface_physics", sf_surface)
                self.check_supported("bl_pbl_physics", bl_pbl)
                self.check_supported("cu_physics", cu_physics)
                self.check_supported("ra_lw_physics", ra_lw)
                self.check_supported("ra_sw_physics", ra_sw)

                if sf_surface == 4:
                    noahmp_keys = (
                        "dveg", "opt_crs", "opt_btr", "opt_runsrf", "opt_runsub",
                        "opt_sfc", "opt_frz", "opt_inf", "opt_rad", "opt_alb",
                        "opt_snf", "opt_tbot", "opt_stc", "opt_gla", "opt_rsf",
                        "opt_soil", "opt_pedo", "opt_crop", "opt_irr", "opt_irrm",
                        "opt_infdv", "opt_tdrn", "opt_tksno", "opt_compact",
                        "opt_scf", "opt_wetland", "soiltstep",
                    )
                    for nk in noahmp_keys:
                        if self.has_opt("noah_mp", nk):
                            sections["phys"][nk] = self.get_opt("noah_mp", nk)
                        elif self.has_opt("physics", nk):
                            sections["phys"][nk] = self.get_opt("physics", nk)
                    for grp_nm in ("noah_mp", "physics"):
                        if self.has_opt(grp_nm, "opt_run"):
                            opt_run_val = self.get_opt(grp_nm, "opt_run")
                            sections["phys"].setdefault("opt_runsrf", opt_run_val)
                            sections["phys"].setdefault("opt_runsub", opt_run_val)
                    for nk, allowed in SUPPORTED_NOAHMP_OPTIONS.items():
                        if nk in sections["phys"] and sections["phys"][nk] not in allowed:
                            self.error(f"WRF Noah-MP {nk} = {sections['phys'][nk]} is not supported; "
                                       f"CPPWRF supports {list(allowed)}.")
                    if "soiltstep" in sections["phys"] and float(sections["phys"]["soiltstep"]) < 0.0:
                        self.error(f"WRF Noah-MP soiltstep = {sections['phys']['soiltstep']} must be >= 0.")

        # Additional microphysics / cumulus / PBL / surface / radiation options in [phys]
        if self.has_opt("physics", "mp_tend_lim"):
            mp_tend_lim = self.get_opt("physics", "mp_tend_lim", 10.0)
            if float(mp_tend_lim) <= 0.0:
                self.error(f"WRF mp_tend_lim = {mp_tend_lim} must be > 0.")
            elif float(mp_tend_lim) != 10.0:
                sections.setdefault("phys", {})["mp_tend_lim"] = mp_tend_lim

        if self.has_opt("physics", "hail_opt"):
            hail_opt = self.get_opt("physics", "hail_opt", 0)
            if hail_opt not in (0, 1):
                self.error(f"WRF hail_opt = {hail_opt} is not supported; CPPWRF supports 0 and 1.")
            elif hail_opt != 0:
                sections.setdefault("phys", {})["hail_opt"] = hail_opt

        if self.has_opt("physics", "use_mp_re"):
            use_mp_re = self.get_opt("physics", "use_mp_re", 1)
            if use_mp_re not in (0, 1):
                self.error(f"WRF use_mp_re = {use_mp_re} is not supported; CPPWRF supports 0 and 1.")
            elif use_mp_re != 1:
                sections.setdefault("phys", {})["use_mp_re"] = use_mp_re

        if self.has_opt("physics", "ishallow"):
            # module_cu_gfl_wrf: the shallow plume runs for ishallow == 1 only; the other ported schemes ignore it
            if self.get_opt("physics", "ishallow", 0) == 1:
                sections.setdefault("phys", {})["ishallow"] = 1

        if self.has_opt("physics", "clos_choice"):
            clos_choice = self.get_opt("physics", "clos_choice", 0)
            if not isinstance(clos_choice, int) or clos_choice > 16:  # > 16 reads past xf_ens (also in WRF)
                self.error(f"WRF clos_choice = {clos_choice} is not supported; CPPWRF supports up to 16.")
            elif clos_choice != 0:
                sections.setdefault("phys", {})["clos_choice"] = clos_choice

        if self.has_opt("physics", "kfeta_trigger"):
            kfeta_trigger = self.get_opt("physics", "kfeta_trigger", 1)
            if kfeta_trigger not in (1, 3):
                self.error(f"WRF kfeta_trigger = {kfeta_trigger} is not supported; CPPWRF supports 1 and 3.")
            elif kfeta_trigger != 1:
                sections.setdefault("phys", {})["kfeta_trigger"] = kfeta_trigger

        if self.has_opt("physics", "convtrans_avglen_m"):
            convtrans = self.get_opt("physics", "convtrans_avglen_m", 30.0)
            if float(convtrans) <= 0.0:
                self.error(f"WRF convtrans_avglen_m = {convtrans} must be > 0.")
            elif float(convtrans) != 30.0:
                sections.setdefault("phys", {})["convtrans_avglen_m"] = convtrans
            # convtrans_prep (cu_diag = 1: Grell-Freitas): stepave = convtrans_avglen_m*60./dt, then MOD(ktau, stepave)
            # divides by zero when stepave is 0 (WRF stops with a floating-point exception)
            dts = [dt]
            for d in range(1, max(1, self._max_dom)):
                pid = self.get_dom_opt("domains", "parent_id", d, d, fallback_d0=False)
                ratio = self.get_dom_opt("domains", "parent_time_step_ratio", d, 1, fallback_d0=False)
                dts.append(dts[pid - 1] / float(ratio) if isinstance(pid, int) and 1 <= pid <= d else dt)
            for d, dt_d in enumerate(dts):
                if self.get_dom_opt("physics", "cu_physics", d, 0) == 3 and int(float(convtrans) * 60.0 / dt_d) < 1:
                    self.error(f"WRF convtrans_avglen_m = {convtrans} is shorter than domain {d + 1}'s time step "
                               f"({dt_d} s): WRF's convtrans_prep divides by zero (MOD(ktau, 0)).")

        if self.has_opt("physics", "ysu_topdown_pblmix"):
            ysu_td = self.get_opt("physics", "ysu_topdown_pblmix", 1)
            if ysu_td not in (0, 1):
                self.error(f"WRF ysu_topdown_pblmix = {ysu_td} is not supported; CPPWRF supports 0 and 1.")
            elif ysu_td != 1:
                sections.setdefault("phys", {})["ysu_topdown_pblmix"] = ysu_td

        if self.has_opt("physics", "rdmaxalb"):
            rdmaxalb = self.get_opt("physics", "rdmaxalb", True)
            if not bool(rdmaxalb):
                sections.setdefault("phys", {})["rdmaxalb"] = False

        # surface fluxes (surface layer, LSM and the diffusion's surface flux), cloud fractions in RRTMG
        isfflx = self.get_opt("physics", "isfflx", 1)
        if isfflx not in (0, 1):
            self.error(f"WRF isfflx = {isfflx} is not supported; CPPWRF supports 0 and 1.")
        elif isfflx != 1 and (has_physics_suite or diff_opt == 2):  # used by these only
            sections.setdefault("phys", {})["isfflx"] = isfflx

        icloud = self.get_opt("physics", "icloud", 1)
        if icloud not in (0, 1):
            self.error(f"WRF icloud = {icloud} is not supported; CPPWRF supports icloud 0 and 1.")
        elif (ra_lw != 0 or ra_sw != 0) and icloud != 1:
            sections.setdefault("phys", {})["icloud"] = icloud

        if self.has_opt("physics", "cldovrlp"):
            cldovrlp = self.get_opt("physics", "cldovrlp", 2)
            if cldovrlp not in (1, 2, 3, 4, 5):
                self.error(f"WRF cldovrlp = {cldovrlp} is not supported; CPPWRF supports 1..5.")
            elif cldovrlp != 2:
                sections.setdefault("phys", {})["cldovrlp"] = cldovrlp

        if self.has_opt("physics", "idcor"):
            idcor = self.get_opt("physics", "idcor", 0)
            if idcor not in (0, 1):
                self.error(f"WRF idcor = {idcor} is not supported; CPPWRF supports 0 and 1.")
            elif idcor != 0:
                sections.setdefault("phys", {})["idcor"] = idcor

        if self.has_opt("physics", "o3input"):
            o3input = self.get_opt("physics", "o3input", 2)
            if o3input not in (0, 2):
                self.error(f"WRF o3input = {o3input} is not supported; CPPWRF supports 0 and 2.")
            elif o3input != 2:
                sections.setdefault("phys", {})["o3input"] = o3input

        if self.has_opt("physics", "ghg_input"):
            ghg_input = self.get_opt("physics", "ghg_input", 1)
            if ghg_input not in (0, 1):
                self.error(f"WRF ghg_input = {ghg_input} is not supported; CPPWRF supports 0 and 1.")
            elif ghg_input != 1:
                sections.setdefault("phys", {})["ghg_input"] = ghg_input

        if self.has_opt("physics", "aer_opt"):
            aer_opt = self.get_opt("physics", "aer_opt", 0)
            if aer_opt != 0:
                self.error(f"WRF aer_opt = {aer_opt} is not supported; CPPWRF supports aer_opt = 0 only.")

        sf_urban = self.get_opt("physics", "sf_urban_physics", 0)
        if sf_urban != 0:
            self.error(f"WRF sf_urban_physics = {sf_urban} (urban canopy) is not supported by CPPWRF.")

        for unsup_key, desc in UNSUPPORTED_DYNAMICS_SWITCHES.items():
            if self.has_opt("dynamics", unsup_key):
                val = self.get_opt("dynamics", unsup_key, 0)
                if val not in (0, False):
                    self.error(f"WRF {unsup_key} = {val} ({desc}) is not supported by CPPWRF.")

        # nest options with per-domain values that are not ported (every domain's value is checked)
        for group, key, default, desc in (("domains", "interp_method_type", 2, "nest interpolation other than SINT"),
                                          ("domains", "nest_interp_coord", 0, "isobaric nest interpolation"),
                                          ("time_control", "input_from_hires", False, "high-resolution nest terrain")):
            for val in self.get_opt_list(group, key):
                if val != default:
                    self.error(f"WRF {key} = {val} ({desc}) is not supported by CPPWRF.")
        # nests open at t = 0 and run to the end of domain 1 (no delayed or early-ending nests)
        for key in ("start_year", "start_month", "start_day", "start_hour", "start_minute", "start_second",
                    "end_year", "end_month", "end_day", "end_hour", "end_minute", "end_second"):
            vals = self.get_opt_list("time_control", key)[:max(1, self._max_dom)]
            if len(set(vals)) > 1:
                self.error(f"WRF {key} = {vals}: nests with their own start or end time are not supported.")

        for (group, key), (default, desc) in DEFAULT_ONLY_OPTIONS.items():
            for val in self.get_opt_list(group, key)[:max(1, self._max_dom)]:
                if val != default:
                    self.error(f"WRF {key} = {val} ({desc}) is not supported by CPPWRF.")
                    break
        # module_diagnostics_driver: noahmp_output = 3 adds the AC* Noah-MP accumulations to the history output
        if self.get_opt("noah_mp", "noahmp_output", 1) == 3 and 4 in self.get_opt_list("physics", "sf_surface_physics"):
            self.error("WRF noahmp_output = 3 (Noah-MP AC* accumulation outputs) is not supported by CPPWRF.")
        for group, key in SAME_ON_ALL_DOMAINS:
            vals = self.get_opt_list(group, key)[:max(1, self._max_dom)]
            if len(set(vals)) > 1:
                self.error(f"WRF {key} = {vals}: CPPWRF uses one value for all domains.")

        # every domain's value: sf_lake_physics, topo_wind, shalwater_z0, ... are max_domains options in the Registry
        for unsup_key, desc in UNSUPPORTED_PHYSICS_SWITCHES.items():
            for val in self.get_opt_list("physics", unsup_key)[:max(1, self._max_dom)]:
                if val not in (0, False):
                    self.error(f"WRF {unsup_key} = {val} ({desc}) is not supported by CPPWRF.")
                    break

        if self.has_opt("bdy_control", "spec_bdy_width"):
            spec_bdy = self.get_opt("bdy_control", "spec_bdy_width", 5)
            if case_name == "real" and max_dom == 1 and spec_bdy != 5:
                self.error(f"WRF spec_bdy_width = {spec_bdy} is not supported; CPPWRF requires spec_bdy_width = 5.")
        # lbc_fcx_gcx (dyn_em/module_bc_em.F:1320): CPPWRF's relaxation weights assume spec_exp = 0
        spec_exp = self.get_opt("bdy_control", "spec_exp", 0.0)
        if float(spec_exp) != 0.0:
            self.error(f"WRF spec_exp = {spec_exp} is not supported; CPPWRF requires spec_exp = 0.")
        # with specified = .true. WRF sets spec_zone = 1, relax_zone = spec_bdy_width - 1
        # (share/module_check_a_mundo.F:248-250); otherwise the namelist values are used
        spec_zone = self.get_opt("bdy_control", "spec_zone", 1)
        relax_zone = self.get_opt("bdy_control", "relax_zone", 4)
        if not specified and (spec_zone != 1 or relax_zone != 4):
            self.error(f"WRF spec_zone = {spec_zone}, relax_zone = {relax_zone} are not supported; "
                       "CPPWRF requires spec_zone = 1 and relax_zone = 4.")

        # 12. Output Section [output]
        sections["output"] = {}
        # maybe_remove_colons (share/module_io_domain.F:468): ':' -> '_' in the wrfout / wrfrst file names
        if self.get_opt("time_control", "nocolons", False):
            self.error("WRF nocolons = .true. is not supported (CPPWRF's file names keep the colons).")
        if case_name.startswith("em_") and self.namelist_path:
            # ideal cases: the history / restart files take WRF's namelist values (radt, cudt, iofields, ...)
            sections["output"]["wrf_namelist"] = os.path.abspath(self.namelist_path)
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
        if restart_on:
            # WRF counts run_* from the restart time; CPPWRF resumes with run.restart_file and run.run_time
            # counted from the simulation start
            self.error("WRF restart = .true. is not translated: resume CPPWRF with --set run.restart_file=<wrfrst file> "
                       "and run.run_time counted from the simulation start.")
        # restart files are written whenever restart_interval > 0 (Registry default 0), whatever restart is;
        # a non-zero restart_interval_m takes precedence (share/set_timekeeping.F:400)
        restart_m = self.get_opt("time_control", "restart_interval_m", 0) or \
            self.get_opt("time_control", "restart_interval", 0)
        for unit_key in ("restart_interval_s", "restart_interval_h", "restart_interval_d", "restart_begin",
                         "restart_begin_s", "restart_begin_m", "restart_begin_h", "restart_begin_d",
                         "restart_begin_y"):
            if self.get_opt("time_control", unit_key, 0):
                self.error(f"WRF {unit_key} is not supported; use restart_interval (minutes).")
        # main/module_wrf_top.F:416: WRF also writes a restart file of the initial state
        if self.get_opt("time_control", "write_restart_at_0h", False):
            self.error("WRF write_restart_at_0h = .true. is not supported.")
        if restart_m and restart_m > 0:
            # the restart file of every domain is wrfrst_d0N_<date> in NetCDF (med_restart_out, output_restart)
            io_form_rst = self.get_opt("time_control", "io_form_restart", 2)
            if io_form_rst != 2:
                self.error(f"WRF io_form_restart = {io_form_rst} is not supported; CPPWRF writes NetCDF (2).")
            rst_outname = self.get_opt("time_control", "rst_outname", "wrfrst_d<domain>_<date>")
            if rst_outname != "wrfrst_d<domain>_<date>":
                self.error(f"WRF rst_outname = '{rst_outname}' is not supported; CPPWRF writes wrfrst_d<domain>_<date>.")
            # nests: med_restart_out writes every domain when d01's restart alarm rings, each domain's own alarm
            # (DomainRunner) rings at the same times only when the interval is a multiple of d01's time step
            steps = restart_m * 60 / dt_val
            if max_dom > 1 and abs(steps - round(steps)) > 1e-9 * steps:
                self.error(f"WRF restart_interval = {restart_m} min is not a multiple of time_step = {dt_val:g} s "
                           "in a nested run: not supported.")
            sections["restart"] = {
                "interval": int(restart_m * 60)
            }

        # 14. Nested Domains [domains] and [d02]..[d0N]
        for mv_key in ("vortex_interval", "max_vortex_speed", "corral_dist", "track_level",
                       "time_to_move", "move_id", "move_interval", "move_cd_x", "move_cd_y"):
            if self.has_opt("domains", mv_key):
                vals = self.get_opt_list("domains", mv_key)
                if any(v not in (0, False) for v in vals):
                    self.error(f"WRF moving nest option &domains.{mv_key} = {vals} is not supported by CPPWRF.")

        for vr_key in ("vert_refine_method", "vert_refine_fact"):
            if self.has_opt("domains", vr_key):
                vals = self.get_opt_list("domains", vr_key)
                if any((v not in (0, 1) if vr_key == "vert_refine_fact" else v != 0) for v in vals):
                    self.error(f"WRF vertical nest refinement (&domains.{vr_key} = {vals}) is not supported by CPPWRF.")

        if max_dom > 1:
            if case_name != "real":
                self.error(f"Nested domains (max_dom = {max_dom}) are only supported for real-data cases (got '{case_name}').")

            feedback = self.get_dom_opt("domains", "feedback", 0, 1)
            if feedback not in (0, 1):
                self.error(f"WRF feedback = {feedback} is not supported; CPPWRF supports 0 (one-way) and 1 (two-way).")

            smooth_option = self.get_dom_opt("domains", "smooth_option", 0, 2)  # Registry default 2
            if smooth_option not in (0, 1, 2):
                self.error(f"WRF smooth_option = {smooth_option} is not supported; CPPWRF supports 0, 1, and 2.")

            sections["domains"] = {
                "max_dom": max_dom,
                "feedback": feedback,
                "smooth_option": smooth_option,
            }

            if self.has_opt("domains", "save_topo_from_real"):  # a wrfinput variable, not a namelist option
                self.error("save_topo_from_real is not a WRF namelist option (WRF reads it from wrfinput_d0N).")

            if self.has_opt("domains", "blend_width"):
                blend_w = self.get_dom_opt("domains", "blend_width", 0, 5)
                if not isinstance(blend_w, int) or blend_w < 0:
                    self.error(f"WRF blend_width = {blend_w} must be an integer >= 0.")
                elif blend_w != 5:
                    sections["domains"]["blend_width"] = blend_w

            if self.has_opt("bdy_control", "spec_bdy_width"):
                spec_bdy = self.get_dom_opt("bdy_control", "spec_bdy_width", 0, 5)
                if not isinstance(spec_bdy, int) or spec_bdy < 1:
                    self.error(f"WRF spec_bdy_width = {spec_bdy} must be an integer >= 1.")
                elif spec_bdy != 5:
                    sections["domains"]["spec_bdy_width"] = spec_bdy

            if self.has_opt("bdy_control", "nested"):
                self.get_opt_list("bdy_control", "nested")

            dom_dt: List[float] = [dt_val] * max_dom
            dom_out_iv: List[int] = [out_interval] * max_dom
            dom_diag_iv: List[int] = [diag_interval] * max_dom

            for d in range(max_dom):
                dom_id = d + 1
                tag = f"d{dom_id:02d}"
                # a missing entry is the Registry default .false. (the nest interpolated from its parent)
                if d > 0 and not self.get_dom_opt("time_control", "input_from_file", d, False, fallback_d0=False):
                    self.error(f"Domain {tag}: input_from_file = .false. (also the default of a missing entry) is "
                               f"not supported; CPPWRF requires wrfinput_{tag}.")

                if d == 0:
                    continue

                parent_id = self.get_dom_opt("domains", "parent_id", d, dom_id - 1, fallback_d0=False)
                if not isinstance(parent_id, int) or parent_id < 1 or parent_id >= dom_id:
                    self.error(f"Domain {tag}: invalid parent_id = {parent_id} (must satisfy 1 <= parent_id < {dom_id}).")
                    p_idx = 0
                else:
                    p_idx = parent_id - 1

                i_start = self.get_dom_opt("domains", "i_parent_start", d, None, fallback_d0=False)
                j_start = self.get_dom_opt("domains", "j_parent_start", d, None, fallback_d0=False)
                grid_ratio = self.get_dom_opt("domains", "parent_grid_ratio", d, None, fallback_d0=False)
                ts_ratio = self.get_dom_opt("domains", "parent_time_step_ratio", d, 1, fallback_d0=False)  # Registry default 1

                # SINT's stencil reaches 2 parent cells beyond the nest; nearer the parent's edge WRF reads its
                # intermediate grid there (not verified): nests start at parent cell 3 and end 2 cells inside
                if i_start is None or not isinstance(i_start, int) or i_start < 3:
                    self.error(f"Domain {tag}: i_parent_start = {i_start} is not supported (must be an integer >= 3).")
                    if not isinstance(i_start, int) or i_start < 1:
                        i_start = 2
                if j_start is None or not isinstance(j_start, int) or j_start < 3:
                    self.error(f"Domain {tag}: j_parent_start = {j_start} is not supported (must be an integer >= 3).")
                    if not isinstance(j_start, int) or j_start < 1:
                        j_start = 2
                if grid_ratio is None or not isinstance(grid_ratio, int) or grid_ratio < 2:
                    self.error(f"Domain {tag}: parent_grid_ratio = {grid_ratio} is invalid (must be an integer >= 2).")
                    if not isinstance(grid_ratio, int) or grid_ratio < 1:
                        grid_ratio = 3
                if ts_ratio is None or not isinstance(ts_ratio, int) or ts_ratio < 1:
                    self.error(f"Domain {tag}: parent_time_step_ratio = {ts_ratio} is invalid (must be an integer >= 1).")
                    if not isinstance(ts_ratio, int) or ts_ratio < 1:
                        ts_ratio = grid_ratio

                # Validate child geometry against parent domain when dimensions are provided in &domains
                if self.has_opt("domains", "e_we"):
                    c_e_we = self.get_dom_opt("domains", "e_we", d, None, fallback_d0=False)
                    p_e_we = self.get_dom_opt("domains", "e_we", p_idx, None, fallback_d0=False)
                    c_s_we = self.get_dom_opt("domains", "s_we", d, 1)
                    p_s_we = self.get_dom_opt("domains", "s_we", p_idx, 1)
                    if c_e_we is not None and p_e_we is not None:
                        c_nx = c_e_we - c_s_we
                        p_nx = p_e_we - p_s_we
                        if c_nx <= 0:
                            self.error(f"Domain {tag}: invalid e_we = {c_e_we} (must be > {c_s_we}).")
                        elif c_nx % grid_ratio != 0:
                            self.error(f"Domain {tag}: (e_we - s_we) = {c_nx} is not divisible by parent_grid_ratio = {grid_ratio}.")
                        elif i_start + c_nx // grid_ratio - 1 > p_nx - 2:
                            self.error(f"Domain {tag}: nest X extent (i_parent_start={i_start}, nx={c_nx}, ratio={grid_ratio}) "
                                       f"reaches within 2 cells of parent domain d{parent_id:02d}'s edge (nx={p_nx}).")

                if self.has_opt("domains", "e_sn"):
                    c_e_sn = self.get_dom_opt("domains", "e_sn", d, None, fallback_d0=False)
                    p_e_sn = self.get_dom_opt("domains", "e_sn", p_idx, None, fallback_d0=False)
                    c_s_sn = self.get_dom_opt("domains", "s_sn", d, 1)
                    p_s_sn = self.get_dom_opt("domains", "s_sn", p_idx, 1)
                    if c_e_sn is not None and p_e_sn is not None:
                        c_ny = c_e_sn - c_s_sn
                        p_ny = p_e_sn - p_s_sn
                        if c_ny <= 0:
                            self.error(f"Domain {tag}: invalid e_sn = {c_e_sn} (must be > {c_s_sn}).")
                        elif c_ny % grid_ratio != 0:
                            self.error(f"Domain {tag}: (e_sn - s_sn) = {c_ny} is not divisible by parent_grid_ratio = {grid_ratio}.")
                        elif j_start + c_ny // grid_ratio - 1 > p_ny - 2:
                            self.error(f"Domain {tag}: nest Y extent (j_parent_start={j_start}, ny={c_ny}, ratio={grid_ratio}) "
                                       f"reaches within 2 cells of parent domain d{parent_id:02d}'s edge (ny={p_ny}).")

                if self.has_opt("domains", "e_vert"):
                    c_ev = self.get_dom_opt("domains", "e_vert", d, None, fallback_d0=False)
                    p_ev = self.get_dom_opt("domains", "e_vert", p_idx, None, fallback_d0=False)
                    if c_ev is not None and p_ev is not None and c_ev != p_ev:
                        self.error(f"Domain {tag}: e_vert = {c_ev} must match parent domain e_vert = {p_ev}.")

                # a nest's dx / dy entries are not used: check_a_mundo sets dx(1) / dy(1) over the accumulated
                # parent_grid_ratio (the model passes that to the nest's wrfinput check, src/cases/nest_driver.cpp)

                # Validate child physics options and ensure consistency with root [phys]
                phys_over = {}
                for pk in ("mp_physics", "sf_sfclay_physics", "sf_surface_physics", "bl_pbl_physics",
                           "cu_physics", "ra_lw_physics", "ra_sw_physics", "radt", "bldt", "cudt",
                           "ishallow", "clos_choice", "kfeta_trigger", "sf_urban_physics"):
                    if self.has_opt("physics", pk):
                        v_d = self.get_dom_opt("physics", pk, d, None)
                        v_0 = self.get_dom_opt("physics", pk, 0, None)
                        if pk == "sf_urban_physics" and v_d not in (0, None):
                            self.error(f"Domain {tag}: WRF sf_urban_physics = {v_d} is not supported by CPPWRF.")
                        elif pk in SUPPORTED_OPTIONS and v_d is not None and v_d not in SUPPORTED_OPTIONS[pk]:
                            self.error(f"Domain {tag}: WRF {pk} = {v_d} is not supported by CPPWRF.")
                        elif v_d != v_0 and pk in PER_DOMAIN_PHYS:
                            phys_over[pk] = v_d  # [d0N] phys_<key>
                        elif v_d != v_0:
                            self.error(f"Domain {tag}: per-domain physics override &physics.{pk} = {v_d} "
                                       f"(domain 1 has {v_0}) is not supported; only "
                                       f"{', '.join(sorted(PER_DOMAIN_PHYS))} may differ between domains.")

                if dom_id == 2:
                    sections["domains"]["i_parent_start"] = i_start
                    sections["domains"]["j_parent_start"] = j_start
                    sections["domains"]["parent_grid_ratio"] = grid_ratio
                    sections["domains"]["parent_time_step_ratio"] = ts_ratio

                sections[tag] = {}
                if dom_id >= 3 or parent_id != dom_id - 1:
                    sections[tag]["parent_id"] = parent_id
                    sections[tag]["i_parent_start"] = i_start
                    sections[tag]["j_parent_start"] = j_start
                    sections[tag]["parent_grid_ratio"] = grid_ratio
                    sections[tag]["parent_time_step_ratio"] = ts_ratio

                sections[tag]["wrfinput"] = f"wrfinput_{tag}"
                for pk, v in phys_over.items():
                    sections[tag]["phys_" + pk] = v
                child_dt = dom_dt[p_idx] / float(ts_ratio)  # the model derives it the same way (no [d0N] dt)
                dom_dt[d] = child_dt

                hist_s_d = self.get_dom_opt("time_control", "history_interval_s", d, 0, fallback_d0=False)
                hist_m_d = self.get_dom_opt("time_control", "history_interval", d, 0, fallback_d0=False)
                if hist_s_d and hist_s_d > 0:
                    child_out_iv = max(1, round(float(hist_s_d) / child_dt))
                elif hist_m_d and hist_m_d > 0:
                    child_out_iv = max(1, round(float(hist_m_d * 60) / child_dt))
                else:
                    child_out_iv = dom_out_iv[p_idx] * ts_ratio
                child_diag_iv = 1 if child_out_iv == 1 else min(child_out_iv, dom_diag_iv[p_idx] * ts_ratio)
                dom_out_iv[d] = child_out_iv
                dom_diag_iv[d] = child_diag_iv
                sections[tag]["out_interval"] = child_out_iv
                sections[tag]["diag_interval"] = child_diag_iv

                if self.has_opt("dynamics", "non_hydrostatic"):
                    if not self.get_dom_opt("dynamics", "non_hydrostatic", d, True):
                        self.error(f"Domain {tag}: non_hydrostatic = .false. is unsupported; CPPWRF is strictly non-hydrostatic.")

                dyn_dom_map = (
                    ("diff_opt", "diff_opt"),
                    ("km_opt", "km_opt"),
                    ("diff_6th_opt", "diff_6th_opt"),
                    ("diff_6th_factor", "diff_6th_factor"),
                    ("zdamp", "zdamp"),
                    ("dampcoef", "dampcoef"),
                    ("khdif", "khdif"),
                    ("kvdif", "kvdif"),
                    ("moist_adv_opt", "moist_adv_opt"),
                    ("scalar_adv_opt", "scalar_adv_opt"),
                    ("tke_adv_opt", "tke_adv_opt"),
                    ("h_mom_adv_order", "h_adv_order"),
                    ("v_mom_adv_order", "v_adv_order"),
                    ("c_s", "c_s"),
                    ("c_k", "c_k"),
                    ("gwd_opt", "gwd_opt"),
                    ("epssm", "epssm"),
                    ("smdiv", "smdiv"),
                    ("emdiv", "emdiv"),
                    ("mix_isotropic", "mix_isotropic"),
                    ("mix_upper_bound", "mix_upper_bound"),
                    ("mix_full_fields", "mix_full_fields"),
                    ("tke_drag_coefficient", "tke_drag_coefficient"),
                    ("tke_heat_flux", "tke_heat_flux"),
                    ("tke_upper_bound", "tke_upper_bound"),
                    ("moist_mix2_off", "moist_mix2_off"),
                    ("scalar_mix2_off", "scalar_mix2_off"),
                    ("tke_mix2_off", "tke_mix2_off"),
                    ("moist_mix6_off", "moist_mix6_off"),
                    ("scalar_mix6_off", "scalar_mix6_off"),
                    ("tke_mix6_off", "tke_mix6_off"),
                )
                for nml_k, ini_k in dyn_dom_map:
                    if self.has_opt("dynamics", nml_k):
                        v_d = self.get_dom_opt("dynamics", nml_k, d, None)
                        v_0 = self.get_dom_opt("dynamics", nml_k, 0, None)
                        if ini_k in SUPPORTED_OPTIONS and v_d is not None and v_d not in SUPPORTED_OPTIONS[ini_k] and not (ini_k == "km_opt" and v_d == 0):
                            self.error(f"Domain {tag}: WRF {nml_k} = {v_d} is not supported.")
                        elif ini_k in ("moist_adv_opt", "scalar_adv_opt", "tke_adv_opt", "mix_isotropic") and v_d not in (0, 1):
                            self.error(f"Domain {tag}: WRF {nml_k} = {v_d} is not supported (0 or 1 only).")
                        elif ini_k == "diff_6th_opt" and v_d not in (0, 1, 2):
                            self.error(f"Domain {tag}: WRF diff_6th_opt = {v_d} is not supported (0, 1, or 2 only).")
                        elif ini_k in ("diff_6th_factor", "mix_upper_bound", "tke_upper_bound") and float(v_d) <= 0.0:
                            self.error(f"Domain {tag}: WRF {nml_k} = {v_d} must be > 0.")
                        if v_d != v_0 and v_d is not None:
                            sections[tag][ini_k] = int(v_d) if isinstance(v_d, (int, float)) and float(v_d).is_integer() and ini_k in ("khdif", "kvdif", "zdamp") else v_d
                            if ini_k in ("khdif", "kvdif") and self.get_dom_opt("dynamics", "diff_opt", d, 0) == 0:
                                sections[tag]["diff_opt"] = 0  # WRF: no diffusion; CPPWRF's default with khdif/kvdif is 1
                # the other &dynamics max_domains options the translator reads have one value in CPPWRF; WRF runs each
                # domain with its own entry (a missing entry is the Registry default, pad_max_domains)
                for nml_k in ("h_sca_adv_order", "v_sca_adv_order", "top_lid", "chem_adv_opt", "momentum_adv_opt",
                              "phi_adv_z", *UNSUPPORTED_DYNAMICS_SWITCHES):
                    if self.has_opt("dynamics", nml_k):
                        v_d = self.get_dom_opt("dynamics", nml_k, d, None)
                        v_0 = self.get_dom_opt("dynamics", nml_k, 0, None)
                        if v_d != v_0:
                            self.error(f"Domain {tag}: WRF {nml_k} = {v_d} differs from domain 1's {v_0}; "
                                       "CPPWRF uses domain 1's value on every domain.")

                if self.has_opt("dynamics", "time_step_sound"):
                    tss_d = self.get_dom_opt("dynamics", "time_step_sound", d, 0)
                    tss_0 = self.get_dom_opt("dynamics", "time_step_sound", 0, 0)
                    if tss_d != tss_0:
                        sections[tag]["time_step_sound"] = tss_d

        # 15. Check for unused / ignored options
        wrf_groups = {g for g, _ in WRF_NAMELIST_OPTIONS}  # a group WRF does not read is not checked (WRF skips it)
        for grp_name, grp_dict in self.nml.items():
            for key_name, vals in grp_dict.items():
                if grp_name in wrf_groups and (grp_name, key_name) not in WRF_NAMELIST_OPTIONS:
                    # module_io: a namelist read error ("ERROR reading namelist ...") is fatal in wrf.exe
                    self.error(f"&{grp_name}.{key_name} is not a WRF 4.8.0 namelist option (wrf.exe stops reading "
                               f"the namelist).")
                elif (grp_name, key_name) not in self.used_keys:
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

        section_order = ["run", "real", "domain", "case", "time", "dyn", "moist", "phys", "restart", "output", "tracers", "domains"]
        remaining = [s for s in sections if s not in section_order]

        for sec in section_order + remaining:
            if sec not in sections or not sections[sec]:
                continue
            lines.append(f"[{sec}]")
            for k, v in sections[sec].items():
                text = _format_num(v)
                if "#" in text or ";" in text:  # the INI reader cuts a value at a comment character
                    self.error(f"{sec}.{k} = '{text}': '#' and ';' cannot appear in a CPPWRF configuration value")
                lines.append(f"{k} = {text}")
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
    parser.add_argument("--case", help="Explicitly specify CPPWRF case (e.g. real, em_hill2d_x, em_quarter_ss, straka)")
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
