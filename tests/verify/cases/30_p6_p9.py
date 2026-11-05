"""P6 (Thompson), P7 (New Tiedtke), P8 (Grell-Freitas-Li), P9 (MYJ + Eta).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
    # P6 = P5 with Thompson (mp_physics = 8; WRF r8 built with -fdefault-double-8, see PHYSICS_SUITE.md
    # section 11): bit-identical for 3 steps and 1 h
    ("real_p6_r8_fp64_cpu", "cpu_dp", "configs/real_p6_r8.ini", [], f"{WR}/run_P6_r8/steps/wrfout_d01_2026-10-12_00:00:00",
     180, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0}),
    ("real_p6_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p6_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P6_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0, "QNICE": 0.0,
      "QNRAIN": 0.0}),
    ("real_p6_fp32_1h", "cuda", "configs/real_p6_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P6/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 5e-4, "W": 1e-2, "THM": 1e-5, "QV": 5e-4}),
    ("real_p6_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p6_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P6_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # P7 = P2 with New Tiedtke (cu_physics = 16) instead of KF; reference WRF-4.8.0-r8d8 (r8 + double-8).
    # The scheme has thresholds that roundoff flips (one column at step 4 jumps 1.3e-4 K under 1e-15 input noise).
    ("real_p7_r8_fp64_cpu", "cpu_dp", "configs/real_p7_r8.ini", [],
     f"{WR}/run_P7_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p7_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p7_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P7_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P7_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p7_fp32_1h", "cuda", "configs/real_p7_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P7", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P7/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-3, "W": 5e-2, "THM": 3e-5, "QV": 8e-4}),
    # P8 = P2 with Grell-Freitas-Li (cu_physics = 3, GFLSCHEME, ishallow = 1); reference WRF-4.8.0-r8d8.
    # Column replay of cu_grell_freitas_li through WRF's objects is bit-identical (tests/phys/cu_gfl).
    ("real_p8_r8_fp64_cpu_10min", "cpu_dp", "configs/real_p8_r8.ini", ["run.run_time=600"],
     f"{WR}/run_P8_r8/steps/wrfout_d01_2026-10-12_00:00:00", 600,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    # P8 CUDA: Grell-Freitas thresholds amplify the GPU's last-bit differences, which depend on the architecture
    # (deterministic per GPU): QV 5.4e-9 on L4 / RTX 2060, 1.5e-7 on G4 (Blackwell: nvcc's sm_120 code contracts
    # other products of the surface layer and Noah into FMAs; its compute_89 PTX or --fmad=false gives the L4 result)
    ("real_p8_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p8_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P8_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P8_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p8_fp32_1h", "cuda", "configs/real_p8_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-3, "W": 5e-2, "THM": 3e-5, "QV": 8e-4}),
    # P9 = P1 with MYJ PBL (bl_pbl_physics = 2) + Eta surface layer (sf_sfclay_physics = 2);
    # reference WRF real*8. Bit-identical on CPU FP64 for both 3 steps and 1 h.
    ("real_p9_r8_fp64_cpu", "cpu_dp", "configs/real_p9_r8.ini", [],
     f"{WR}/run_P9_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p9_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p9_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P9_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P9_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p9_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p9_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P9_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P9_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p9_fp32_1h", "cuda", "configs/real_p9_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P9", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P9/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 6e-4, "W": 1e-2, "THM": 1e-5, "QV": 2e-3}),
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('P6', 'P5 with Thompson', '1 h'),
    ('P7', 'P2 with New Tiedtke', '1 h'),
    ('P8', 'P2 with Grell-Freitas-Li (deep and shallow)', '1 h'),
    ('P9', 'P1 with MYJ PBL and the Eta surface layer', '1 h'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {
    'real_p6_r8_fp64_cpu': 'port_a06_data.tar,port_c01_data.tar,v05_files_data.tar,v05_files_links.tar',
    'real_p6_r8_fp64_cpu_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs.tar',
    'real_p6_fp32_1h': 'suite_data.tar',
    'real_p6_r8_fp64_cuda_1h': 'port_a06_data.tar,port_c01_data.tar,v05_suite_refs.tar',
    'real_p7_r8_fp64_cpu': 'v05_files_data.tar,v05_files_links.tar',
    'real_p7_r8_fp64_cuda_1h': 'suite_data.tar,v05_suite_refs.tar',
    'real_p7_fp32_1h': 'suite_data.tar',
    'real_p8_r8_fp64_cpu_10min': 'v05_files_data.tar,v05_files_links.tar',
    'real_p8_r8_fp64_cuda_1h': 'suite_data.tar,v05_suite_refs.tar',
    'real_p8_fp32_1h': 'suite_data.tar',
    'real_p9_r8_fp64_cpu': 'v05_files_data.tar,v05_files_links.tar',
    'real_p9_r8_fp64_cpu_1h': 'suite_data.tar,v05_suite_refs.tar',
    'real_p9_r8_fp64_cuda_1h': 'suite_data.tar,v05_suite_refs.tar',
    'real_p9_fp32_1h': 'suite_data.tar',
}
