"""MYNN-EDMF PBL and the MYNN surface layer: MYNN, MYNNSF1 / MYNNSF2, MYNNX1 / MYNNX2.

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
    # MYNN = P6 with the MYNN surface layer (sf_sfclay_physics = 5), the MYNN-EDMF PBL (bl_pbl_physics = 5,
    # icloud_bl = 1) and RRTMG every 10 min; reference WRF-4.8.0-r8d8. Bit-identical for 3 steps and 1 h.
    ("real_mynn_r8_fp64_cpu", "cpu_dp", "configs/real_mynn_r8.ini", [],
     f"{WR}/run_MYNN_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynn_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynn_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_MYNN_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNN_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynn_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynn_r8.ini",  # L4: U 2.6e-4 W 5.8e-3 THM 5.5e-6 QV 1.1e-4
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_MYNN_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNN_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # MYNNSF1 / MYNNSF2 = MYNN with the MYNN-EDMF PBL on the revised MM5 (sf_sfclay_physics = 1) / Eta (2) surface
    # layer instead of the MYNN one, from the P6 inputs; reference WRF-4.8.0-r8d8. Bit-identical for 3 steps and 1 h.
    ("real_mynnsf1_r8_fp64_cpu", "cpu_dp", "configs/real_mynnsf1_r8.ini", [],
     f"{WR}/run_MYNNSF1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnsf1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynnsf1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNSF1_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnsf1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynnsf1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNSF1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),  # L4: U 2.5e-4 W 5.3e-3 THM 5.3e-6 QV 2.9e-4, G4: U 2.7e-4
    ("real_mynnsf2_r8_fp64_cpu", "cpu_dp", "configs/real_mynnsf2_r8.ini", [],
     f"{WR}/run_MYNNSF2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnsf2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynnsf2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNSF2_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnsf2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynnsf2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNSF2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),  # L4: U 2.6e-4 W 5.4e-3 THM 5.6e-6 QV 1.1e-4
    # MYNNX1 / MYNNX2 = MYNN (P6 inputs, MYNN-EDMF PBL and MYNN surface layer, RRTMG radt 10) with the non-default
    # MYNN options: X1 bl_mynn_tkeadvect, cloudpdf 1, mixlength 2, closure 3.0, edmf_tke 1, mixqt 1, output 1,
    # tke_budget 1, sf_mynn_sfcflux_land 2 / _water 3, grav_settling 2; X2 bl_mynn_edmf 0, cloudmix 0, cloudpdf 0,
    # mixlength 0, closure 2.5, ess 0, mixscalars 1, sf_mynn_sfcflux_land 3 / _water 4, grav_settling 1.
    # Reference WRF-4.8.0-r8d8; bit-identical for 3 steps and 1 h. CUDA FP64 L4 / G4: X1 U 3.2e-4 / 3.1e-4,
    # W 8.6e-3 / 5.9e-3, THM 9.0e-6 / 5.4e-6, QV 1.4e-4 / 1.3e-4; X2 U 2.1e-4, W 6.3e-3 / 7.3e-3, THM 4.4e-6 / 4.8e-6,
    # QV 6.4e-5 / 6.7e-5.
    ("real_mynnx1_r8_fp64_cpu", "cpu_dp", "configs/real_mynnx1_r8.ini", [],
     f"{WR}/run_MYNNX1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnx1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynnx1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNX1_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnx1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynnx1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNX1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_mynnx2_r8_fp64_cpu", "cpu_dp", "configs/real_mynnx2_r8.ini", [],
     f"{WR}/run_MYNNX2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnx2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynnx2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNX2_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnx2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynnx2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNX2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('MYNN', 'P6 with the MYNN surface layer and the MYNN-EDMF PBL (icloud_bl 1) + RRTMG (radt 10 min)', '1 h'),
    ('MYNNSF1 / MYNNSF2', 'MYNN with the revised MM5 / Eta surface layer instead of the MYNN one', '1 h'),
    ('MYNNX1 / MYNNX2',
     'MYNN with non-default MYNN options (X1: tkeadvect, cloudpdf 1, mixlength 2, closure 3.0, edmf_tke 1,'
     ' mixqt 1, output 1, tke_budget 1, sfcflux land 2 / water 3, grav_settling 2; X2: edmf 0, cloudmix 0,'
     ' cloudpdf 0, mixlength 0, closure 2.5, ess 0, mixscalars 1, sfcflux land 3 / water 4, grav_settling '
     '1)',
     '1 h'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {
    'real_mynn_r8_fp64_cpu': 'mynn_data.tar,v05_files_data.tar,v05_files_links.tar',
    'real_mynn_r8_fp64_cpu_1h': 'mynn_data.tar,v05_suite_refs.tar',
    'real_mynn_r8_fp64_cuda_1h': 'mynn_data.tar,v05_suite_refs.tar',
    'real_mynnsf1_r8_fp64_cpu': 'port_a06_data.tar,port_c01_data.tar,v05_iofields.tar,v05_files_data.tar,v05_files_links.tar',
    'real_mynnsf1_r8_fp64_cpu_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs.tar,v05_iofields.tar',
    'real_mynnsf1_r8_fp64_cuda_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs.tar,v05_iofields.tar',
    'real_mynnsf2_r8_fp64_cpu': 'port_a06_data.tar,port_c01_data.tar,v05_iofields.tar,v05_files_data.tar,v05_files_links.tar',
    'real_mynnsf2_r8_fp64_cpu_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs.tar,v05_iofields.tar',
    'real_mynnsf2_r8_fp64_cuda_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs.tar,v05_iofields.tar',
    'real_mynnx1_r8_fp64_cpu': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs2.tar,v05_iofields.tar',
    'real_mynnx1_r8_fp64_cpu_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs2.tar,v05_iofields.tar',
    'real_mynnx1_r8_fp64_cuda_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs2.tar,v05_iofields.tar',
    'real_mynnx2_r8_fp64_cpu': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs2.tar,v05_iofields.tar',
    'real_mynnx2_r8_fp64_cpu_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs2.tar,v05_iofields.tar',
    'real_mynnx2_r8_fp64_cuda_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs2.tar,v05_iofields.tar',
}
