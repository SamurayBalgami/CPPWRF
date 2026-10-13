# CPPWRF

A from-scratch C++20 / CUDA re-implementation of the WRF-ARW dynamical core and physics, with
one source for CPU (OpenMP) and GPU (CUDA), FP32 or FP64. Every formula follows WRF v4.6.1 and is
verified against the original: physics schemes column by column against the Fortran compiled
with gfortran, the dynamics against real WRF runs.

Status and results: [docs/STATUS.md](docs/STATUS.md).

## Build

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=OFF
cmake --build build -j
build/cppwrf configs/straka.ini
```

## License

GPL-3.0-or-later, see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md) (includes the WRF notice).
