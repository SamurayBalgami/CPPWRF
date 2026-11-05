"""restart_<tag>_fp64: 3 steps straight vs restarted after step 1 (wrf_suite.py's restart_p1); restart_<tag>_fp64_cuda:
the same in CUDA FP64 (the GPU run is deterministic: restarted = straight bit for bit). (tag, config, input run
directory under ${CPPWRF_WRF_RUNS})."""

RESTARTS = (("p1", "configs/real_p1_r8.ini", "run_P1_r8"), ("a01", "configs/real_a01_r8.ini", "run_P1_r8"),
            ("p7", "configs/real_p7_r8.ini", "run_P7_r8"),
            ("p8", "configs/real_p8_r8.ini", "run_P8_r8"), ("p9", "configs/real_p9_r8.ini", "run_P9_r8"),
            ("nmp1", "configs/real_nmp1_r8.ini", "run_NMP1_r8"), ("mynn", "configs/real_mynn_r8.ini", "run_MYNN_r8"),
            ("mynnsf1", "configs/real_mynnsf1_r8.ini", "run_P6_r8"), ("mynnsf2", "configs/real_mynnsf2_r8.ini", "run_P6_r8"),
            ("mynnx1", "configs/real_mynnx1_r8.ini", "run_P6_r8"), ("mynnx2", "configs/real_mynnx2_r8.ini", "run_P6_r8"),
            ("t1", "configs/real_t1.ini", "run_R1_r8"),
            ("a05", "configs/real_a05_r8.ini", "run_A05_r8"), ("a05g", "configs/real_a05g_r8.ini", "run_A05G_r8"),
            ("slab1", "configs/real_slab1_r8.ini", "run_SLAB1_r8"))
