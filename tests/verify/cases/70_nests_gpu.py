"""Nests on the GPU: N3 and N2 in CUDA FP64, N2 in FP32.

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
    # N3 CUDA (needs > 6 GB; measured on L4 and G4): the P4 physics' roundoff growth, so P4's tolerances
    ("real_n3_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n3.ini", [],
     {d: f"{WR}/run_N3h_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {}),
    ("real_n2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n2.ini", ["real.dir=" + f"{WR}/run_N2_r8"],
     {d: f"{WR}/run_N2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 0.0, "d01.THM": 0.0, "d02.U": 0.0, "d02.THM": 0.0}),
    ("real_n2_fp32_1h", "cuda", "configs/real_n2.ini", [],
     {d: f"{WR}/run_N2/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 5e-4, "d01.THM": 3e-6, "d02.U": 6e-4, "d02.THM": 3e-6}),
]
