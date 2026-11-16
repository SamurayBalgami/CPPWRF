# RRTMG in the working precision: experiment scripts

Results and reasoning: [docs/RRTMG_FP32_EXPERIMENT.md](../../../docs/RRTMG_FP32_EXPERIMENT.md).

**Colab:** [rrtmg_fp32_colab.ipynb](https://colab.research.google.com/github/SamurayBalgami/CPPWRF/blob/exp/rrtmg-fp32/tools/experiments/rrtmg_fp32/rrtmg_fp32_colab.ipynb)
builds upstream c9c95c5 and this branch on a Colab GPU, runs N3 in both and prints the speed-up and the accuracy
against GPU FP64.

| File | What it does |
|---|---|
| `prep_exp1.py` ... `prep_exp5.py` | the source change, applied in this order to a clean upstream checkout (`python3 prep_expN.py CHECKOUT`); they give the change's `src/` exactly |
| `local_cuda.sh` | CUDA 12.9 with glibc >= 2.41: a user-local nvcc whose `crt/math_functions.h` declares sinpi / cospi / rsqrt noexcept |
| `build_gpu.sh` | GPU build of a checkout (`build_gpu.sh CHECKOUT OFF build_dir [-DCMAKE_CUDA_FLAGS=--fmad=false]`) |
| `run_case.sh` | one timed run (`run_case.sh CHECKOUT build_dir configs/real_n3.ini OUT --set ...`) |
| `sass_f64.sh` | FP64 SASS instructions in the RRTMG kernels of a binary |
| `compare_6h.sh` | RMSE tables of the 6 h runs against the FP64 run (`OUT_ROOT=...`) |
| `rrtmg_fp32_colab.ipynb` | the whole comparison on Colab |
