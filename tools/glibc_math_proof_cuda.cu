// cppwrf::glibc (src/core/glibc_math.hpp) on the GPU against the same functions on the host, bit for bit: the device
// computes what the host port computes (and tools/glibc_math_proof.cpp shows the host port equals glibc). Inputs: the
// edge values and random sets of tools/glibc_math_proof.cpp (shared generator, tools/glibc_math_inputs.hpp).
//
//   nvcc -std=c++20 -O3 --fmad=false -arch=native -DCPPWRF_GLIBC_MATH -Isrc -Itools -Xcompiler -fopenmp,-ffp-contract=off
//        tools/glibc_math_proof_cuda.cu -o glibc_math_proof_cuda
//   ./glibc_math_proof_cuda [N] [FUNCTION,...]
// Exit 1 when any result differs.
#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "core/glibc_math.hpp"
#include "glibc_math_inputs.hpp"

using namespace glibc_math_inputs;

template <int F>
__host__ __device__ double eval(double x, double y) {
  return dispatch<F>(x, y);
}

template <int F>
__global__ void kern(const double* x, const double* y, double* r, long long n) {
  const long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
  if (i < n) r[i] = eval<F>(x[i], y[i]);
}

template <int F>
long long run(long long n_random, const std::vector<double>& edges, int& shown) {
  const long long batch = 1 << 22;
  std::vector<double> hx(batch), hy(batch), hr(batch);
  double *dx, *dy, *dr;
  cudaMalloc(&dx, batch * 8);
  cudaMalloc(&dy, batch * 8);
  cudaMalloc(&dr, batch * 8);
  long long bad = 0, total = 0;
  std::uint64_t s = 0x5eed0000ull + F;
  const long long ne = (long long)edges.size();
  const long long n_edges = kArity[F] == 1 ? ne : ne * 64;
  for (long long done = 0; done < n_edges + n_random; done += batch) {
    const long long m = std::min(batch, n_edges + n_random - done);
    for (long long t = 0; t < m; ++t) {
      const long long u = done + t;
      if (u < n_edges) edge_input(F, edges, u, s, hx[t], hy[t]);
      else gen(F, s, int(u % 4), hx[t], hy[t]);
    }
    cudaMemcpy(dx, hx.data(), m * 8, cudaMemcpyHostToDevice);
    cudaMemcpy(dy, hy.data(), m * 8, cudaMemcpyHostToDevice);
    kern<F><<<(m + 255) / 256, 256>>>(dx, dy, dr, m);
    cudaMemcpy(hr.data(), dr, m * 8, cudaMemcpyDeviceToHost);
#pragma omp parallel for reduction(+ : bad) schedule(static)
    for (long long t = 0; t < m; ++t) {
      const double h = eval<F>(hx[t], hy[t]);
      if (!same(h, hr[t])) {
        ++bad;
#pragma omp critical
        if (shown++ < 10)
          std::printf("  %s(%016llx, %016llx): host %016llx, device %016llx\n", kNames[F], (unsigned long long)B(hx[t]),
                      (unsigned long long)B(hy[t]), (unsigned long long)B(h), (unsigned long long)B(hr[t]));
      }
    }
    total += m;
  }
  cudaFree(dx);
  cudaFree(dy);
  cudaFree(dr);
  std::printf("%-10s %13lld inputs, %lld differ\n", kNames[F], total, bad);
  std::fflush(stdout);
  return bad;
}

template <int... F>
long long run_all(long long n, const std::string& only, std::integer_sequence<int, F...>) {
  const std::vector<double> edges = edge_values();
  long long bad = 0;
  int shown = 0;
  ((only.empty() || ("," + only + ",").find("," + std::string(kNames[F]) + ",") != std::string::npos
        ? (void)(bad += run<F>(n, edges, shown))
        : (void)0),
   ...);
  return bad;
}

int main(int argc, char** argv) {
  const long long n = argc > 1 ? std::atoll(argv[1]) : 100000000ll;
  const std::string only = argc > 2 ? argv[2] : "";
  cudaDeviceProp p;
  cudaGetDeviceProperties(&p, 0);
  std::printf("glibc_math_proof_cuda: %s (sm_%d%d), %lld random inputs per function\n", p.name, p.major, p.minor, n);
  const long long bad = run_all(n, only, std::make_integer_sequence<int, kNumFuncs>{});
  std::printf(bad ? "glibc_math_proof_cuda: %lld DIFFER\n" : "glibc_math_proof_cuda: identical\n", bad);
  return bad ? 1 : 0;
}
