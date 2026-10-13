// driver_diffusion.cpp
// C++ differential testing harness for WRF v4.6.1 Diffusion and Damping kernels in CPPWRF.
// Supports both CPU (OpenMP) and CUDA execution with single-source semantics.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "dyn/diffusion/diffusion.hpp"

using namespace cppwrf;
using namespace cppwrf::dyn::diffusion;

template <typename T>
struct DevField {
  std::vector<T> host;
  T* dev = nullptr;
  std::size_t n = 0;

  DevField() = default;
  explicit DevField(std::size_t count) { init(count); }
  ~DevField() {
    if (dev) dev_free(dev);
  }
  DevField(const DevField&) = delete;
  DevField& operator=(const DevField&) = delete;
  DevField(DevField&& o) noexcept : host(std::move(o.host)), dev(o.dev), n(o.n) {
    o.dev = nullptr;
    o.n = 0;
  }
  DevField& operator=(DevField&& o) noexcept {
    if (this != &o) {
      if (dev) dev_free(dev);
      host = std::move(o.host);
      dev = o.dev;
      n = o.n;
      o.dev = nullptr;
      o.n = 0;
    }
    return *this;
  }

  void init(std::size_t count) {
    if (dev) dev_free(dev);
    n = count;
    host.assign(n, T(0));
    dev = static_cast<T*>(dev_alloc(n * sizeof(T)));
  }
  void upload() {
    copy_to_dev(dev, host.data(), n * sizeof(T));
  }
  void download() {
    copy_to_host(host.data(), dev, n * sizeof(T));
  }
  void fill(T val = T(0)) {
    std::fill(host.begin(), host.end(), val);
    upload();
  }
  T* data() { return dev; }
  const T* data() const { return dev; }
  T* host_data() { return host.data(); }
  const T* host_data() const { return host.data(); }
  std::size_t size() const { return n; }
  T& operator[](std::size_t i) { return host[i]; }
  const T& operator[](std::size_t i) const { return host[i]; }
};

static void write_chunk(std::ofstream& out, const std::string& name, const double* data, int count) {
  char tag[32] = {0};
  std::strncpy(tag, name.c_str(), 31);
  out.write(tag, 32);
  out.write(reinterpret_cast<const char*>(&count), sizeof(int));
  out.write(reinterpret_cast<const char*>(data), count * sizeof(double));
}

template <typename T>
static void write_field(std::ofstream& out, const std::string& name, DevField<T>& f) {
  f.download();
  write_chunk(out, name, f.host_data(), static_cast<int>(f.size()));
}

template <typename T>
static void read_field(std::ifstream& in, DevField<T>& f) {
  in.read(reinterpret_cast<char*>(f.host_data()), f.size() * sizeof(T));
  f.upload();
}

inline View3<double> make_view3(double* p, int nx, int nz) {
  return View3<double>{p, nx * nz, nx};
}
inline View3<const double> make_const_view3(const double* p, int nx, int nz) {
  return View3<const double>{p, nx * nz, nx};
}
inline View2<const double> make_view2(const double* p, int nx) {
  return View2<const double>{p, nx};
}
inline BView2<const double> make_bview2(const double* p, int nx) {
  return BView2<const double>{p, 1, nx};
}
inline BView3<const double> make_bview3(const double* p, int nx, int nz) {
  return BView3<const double>{p, 1, nx * nz, nx};
}
inline View1<const double> make_view1(const double* p) {
  return View1<const double>{p};
}

int main(int argc, char* argv[]) {
  std::string in_path = (argc > 1) ? argv[1] : "diffusion_input.bin";
  std::string out_path = (argc > 2) ? argv[2] : "diffusion_output_cpp.bin";

  std::ifstream in(in_path, std::ios::binary);
  if (!in) {
    std::cerr << "Failed to open input file: " << in_path << std::endl;
    return 1;
  }

  int nx, ny, nz;
  double dx, dy, dt;
  in.read(reinterpret_cast<char*>(&nx), sizeof(int));
  in.read(reinterpret_cast<char*>(&ny), sizeof(int));
  in.read(reinterpret_cast<char*>(&nz), sizeof(int));
  in.read(reinterpret_cast<char*>(&dx), sizeof(double));
  in.read(reinterpret_cast<char*>(&dy), sizeof(double));
  in.read(reinterpret_cast<char*>(&dt), sizeof(double));

  // 1D arrays
  DevField<double> dn(nz), rdn(nz), dnw(nz), rdnw(nz), fnm(nz), fnp(nz);
  double cf1, cf2, cf3;
  DevField<double> c1(nz), c2(nz), c3(nz), c4(nz), c1f(nz + 1), c2f(nz + 1);
  DevField<double> u_base(nz), v_base(nz), t_base(nz), qv_base(nz);

  read_field(in, dn);
  read_field(in, rdn);
  read_field(in, dnw);
  read_field(in, rdnw);
  read_field(in, fnm);
  read_field(in, fnp);
  in.read(reinterpret_cast<char*>(&cf1), sizeof(double));
  in.read(reinterpret_cast<char*>(&cf2), sizeof(double));
  in.read(reinterpret_cast<char*>(&cf3), sizeof(double));
  read_field(in, c1);
  read_field(in, c2);
  read_field(in, c3);
  read_field(in, c4);
  read_field(in, c1f);
  read_field(in, c2f);
  read_field(in, u_base);
  read_field(in, v_base);
  read_field(in, t_base);
  read_field(in, qv_base);

  // 2D fields
  DevField<double> mut(nx * ny), muu((nx + 1) * ny), muv(nx * (ny + 1));
  DevField<double> msftx(nx * ny), msfty(nx * ny);
  DevField<double> msfux((nx + 1) * ny), msfuy((nx + 1) * ny);
  DevField<double> msfvx(nx * (ny + 1)), msfvy(nx * (ny + 1));
  DevField<double> ust(nx * ny), hfx(nx * ny), qfx(nx * ny);

  read_field(in, mut);
  read_field(in, muu);
  read_field(in, muv);
  read_field(in, msftx);
  read_field(in, msfty);
  read_field(in, msfux);
  read_field(in, msfuy);
  read_field(in, msfvx);
  read_field(in, msfvy);
  read_field(in, ust);
  read_field(in, hfx);
  read_field(in, qfx);

  DevField<double> msfvx_inv(nx * (ny + 1));
  for (size_t idx = 0; idx < msfvx.size(); ++idx) {
    msfvx_inv[idx] = 1.0 / msfvx[idx];
  }
  msfvx_inv.upload();

  // 3D fields
  DevField<double> phb(nx * (nz + 1) * ny), ph(nx * (nz + 1) * ny);
  DevField<double> u((nx + 1) * nz * ny), v(nx * nz * (ny + 1));
  DevField<double> w(nx * (nz + 1) * ny), ww(nx * (nz + 1) * ny);
  DevField<double> thp(nx * nz * ny), theta(nx * nz * ny), t(nx * nz * ny);
  DevField<double> pb(nx * nz * ny), p(nx * nz * ny);
  DevField<double> alt(nx * nz * ny), rho(nx * nz * ny);
  DevField<double> moist(nx * nz * ny), tke(nx * (nz + 1) * ny);

  read_field(in, phb);
  read_field(in, ph);
  read_field(in, u);
  read_field(in, v);
  read_field(in, w);
  read_field(in, ww);
  read_field(in, thp);
  read_field(in, theta);
  read_field(in, t);
  read_field(in, pb);
  read_field(in, p);
  read_field(in, alt);
  read_field(in, rho);
  read_field(in, moist);
  read_field(in, tke);
  in.close();

  // Outputs
  DevField<double> z(nx * nz * ny), rdz(nx * nz * ny), rdzw(nx * nz * ny);
  DevField<double> zx((nx + 1) * (nz + 1) * ny), zy(nx * (nz + 1) * (ny + 1));
  DevField<double> defor11(nx * nz * ny), defor22(nx * nz * ny), defor33(nx * nz * ny);
  DevField<double> defor12((nx + 1) * nz * (ny + 1));
  DevField<double> defor13((nx + 1) * (nz + 1) * ny);
  DevField<double> defor23(nx * (nz + 1) * (ny + 1));
  DevField<double> div(nx * nz * ny);
  DevField<double> xkmh(nx * nz * ny), xkmv(nx * nz * ny);
  DevField<double> xkhh(nx * nz * ny), xkhv(nx * nz * ny);
  DevField<double> xkmh_k1(nx * nz * ny), xkmv_k1(nx * nz * ny);
  DevField<double> xkhh_k1(nx * nz * ny), xkhv_k1(nx * nz * ny);
  DevField<double> titau11(nx * nz * ny), titau22(nx * nz * ny), titau33(nx * nz * ny);
  DevField<double> titau12((nx + 1) * nz * (ny + 1));
  DevField<double> titau13((nx + 1) * (nz + 1) * ny);
  DevField<double> titau23(nx * (nz + 1) * (ny + 1));
  DevField<double> ru_tend((nx + 1) * nz * ny), rv_tend(nx * nz * (ny + 1));
  DevField<double> rw_tend(nx * (nz + 1) * ny), rt_tend(nx * nz * ny);
  DevField<double> moist_tend(nx * nz * ny);
  DevField<double> dummy3d(nx * nz * ny);
  dummy3d.fill(0.0);
  dummy3d.upload();

  std::ofstream out(out_path, std::ios::binary);

  // ============================================================================
  // TEST 1: compute_diff_metrics
  // Case A: Periodic
  // ============================================================================
  {
    DiffusionParams p;
    p.periodic_x = true;
    p.periodic_y = true;
    p.specified = false;

    compute_diff_metrics(p,
                         make_const_view3(ph.data(), nx, nz + 1),
                         make_bview3(phb.data(), nx, nz + 1),
                         make_view3(z.data(), nx, nz),
                         make_view3(rdz.data(), nx, nz),
                         make_view3(rdzw.data(), nx, nz),
                         make_view3(zx.data(), nx + 1, nz + 1),
                         make_view3(zy.data(), nx, nz + 1),
                         1.0 / dx, 1.0 / dy, nx, ny, nz);

    zx.download();
    zy.download();
    std::vector<double> zx_out((nx + 1) * nz * ny);
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        for (int i = 0; i <= nx; ++i) {
          zx_out[i + k * (nx + 1) + j * (nx + 1) * nz] = zx[i + k * (nx + 1) + j * (nx + 1) * (nz + 1)];
        }
      }
    }
    std::vector<double> zy_out(nx * nz * (ny + 1));
    for (int j = 0; j <= ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
          zy_out[i + k * nx + j * nx * nz] = zy[i + k * nx + j * nx * (nz + 1)];
        }
      }
    }

    write_field(out, "metrics_z_periodic", z);
    write_field(out, "metrics_rdz_periodic", rdz);
    write_field(out, "metrics_rdzw_periodic", rdzw);
    write_chunk(out, "metrics_zx_periodic", zx_out.data(), static_cast<int>(zx_out.size()));
    write_chunk(out, "metrics_zy_periodic", zy_out.data(), static_cast<int>(zy_out.size()));
  }

  // Case B: Specified
  {
    DiffusionParams p;
    p.periodic_x = false;
    p.periodic_y = false;
    p.specified = true;

    compute_diff_metrics(p,
                         make_const_view3(ph.data(), nx, nz + 1),
                         make_bview3(phb.data(), nx, nz + 1),
                         make_view3(z.data(), nx, nz),
                         make_view3(rdz.data(), nx, nz),
                         make_view3(rdzw.data(), nx, nz),
                         make_view3(zx.data(), nx + 1, nz + 1),
                         make_view3(zy.data(), nx, nz + 1),
                         1.0 / dx, 1.0 / dy, nx, ny, nz);

    zx.download();
    zy.download();
    std::vector<double> zx_out((nx + 1) * nz * ny);
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        for (int i = 0; i <= nx; ++i) {
          zx_out[i + k * (nx + 1) + j * (nx + 1) * nz] = zx[i + k * (nx + 1) + j * (nx + 1) * (nz + 1)];
        }
      }
    }
    std::vector<double> zy_out(nx * nz * (ny + 1));
    for (int j = 0; j <= ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
          zy_out[i + k * nx + j * nx * nz] = zy[i + k * nx + j * nx * (nz + 1)];
        }
      }
    }

    write_field(out, "metrics_z_specified", z);
    write_field(out, "metrics_rdz_specified", rdz);
    write_field(out, "metrics_rdzw_specified", rdzw);
    write_chunk(out, "metrics_zx_specified", zx_out.data(), static_cast<int>(zx_out.size()));
    write_chunk(out, "metrics_zy_specified", zy_out.data(), static_cast<int>(zy_out.size()));
  }

  // ============================================================================
  // TEST 2: cal_deform_and_div
  // Case A: Periodic
  // ============================================================================
  {
    DiffusionParams p;
    p.periodic_x = true;
    p.periodic_y = true;
    p.specified = false;

    cal_deform_and_div(p,
                       make_const_view3(u.data(), nx + 1, nz),
                       make_const_view3(v.data(), nx, nz),
                       make_const_view3(w.data(), nx, nz + 1),
                       make_view1(u_base.data()),
                       make_view1(v_base.data()),
                       make_const_view3(zx.data(), nx + 1, nz + 1),
                       make_const_view3(zy.data(), nx, nz + 1),
                       make_const_view3(rdz.data(), nx, nz),
                       make_const_view3(rdzw.data(), nx, nz),
                       make_bview2(msftx.data(), nx),
                       make_bview2(msfty.data(), nx),
                       make_bview2(msfux.data(), nx + 1),
                       make_bview2(msfuy.data(), nx + 1),
                       make_bview2(msfvx.data(), nx),
                       make_bview2(msfvy.data(), nx),
                       make_view1(fnm.data()),
                       make_view1(fnp.data()),
                       make_view1(dn.data()),
                       make_view1(dnw.data()),
                       cf1, cf2, cf3, 1.0 / dx, 1.0 / dy,
                       make_view3(defor11.data(), nx, nz),
                       make_view3(defor22.data(), nx, nz),
                       make_view3(defor33.data(), nx, nz),
                       make_view3(defor12.data(), nx + 1, nz),
                       make_view3(defor13.data(), nx + 1, nz + 1),
                       make_view3(defor23.data(), nx, nz + 1),
                       make_view3(div.data(), nx, nz),
                       nx, ny, nz);

    write_field(out, "defor11_periodic", defor11);
    write_field(out, "defor22_periodic", defor22);
    write_field(out, "defor33_periodic", defor33);
    write_field(out, "defor12_periodic", defor12);
    write_field(out, "defor13_periodic", defor13);
    write_field(out, "defor23_periodic", defor23);
    write_field(out, "div_periodic", div);
  }

  // Case B: Specified
  {
    DiffusionParams p;
    p.periodic_x = false;
    p.periodic_y = false;
    p.specified = true;

    defor11.fill(0.0);
    defor22.fill(0.0);
    defor33.fill(0.0);
    defor12.fill(0.0);
    defor13.fill(0.0);
    defor23.fill(0.0);
    div.fill(0.0);

    cal_deform_and_div(p,
                       make_const_view3(u.data(), nx + 1, nz),
                       make_const_view3(v.data(), nx, nz),
                       make_const_view3(w.data(), nx, nz + 1),
                       make_view1(u_base.data()),
                       make_view1(v_base.data()),
                       make_const_view3(zx.data(), nx + 1, nz + 1),
                       make_const_view3(zy.data(), nx, nz + 1),
                       make_const_view3(rdz.data(), nx, nz),
                       make_const_view3(rdzw.data(), nx, nz),
                       make_bview2(msftx.data(), nx),
                       make_bview2(msfty.data(), nx),
                       make_bview2(msfux.data(), nx + 1),
                       make_bview2(msfuy.data(), nx + 1),
                       make_bview2(msfvx.data(), nx),
                       make_bview2(msfvy.data(), nx),
                       make_view1(fnm.data()),
                       make_view1(fnp.data()),
                       make_view1(dn.data()),
                       make_view1(dnw.data()),
                       cf1, cf2, cf3, 1.0 / dx, 1.0 / dy,
                       make_view3(defor11.data(), nx, nz),
                       make_view3(defor22.data(), nx, nz),
                       make_view3(defor33.data(), nx, nz),
                       make_view3(defor12.data(), nx + 1, nz),
                       make_view3(defor13.data(), nx + 1, nz + 1),
                       make_view3(defor23.data(), nx, nz + 1),
                       make_view3(div.data(), nx, nz),
                       nx, ny, nz);

    write_field(out, "defor11_specified", defor11);
    write_field(out, "defor22_specified", defor22);
    write_field(out, "defor33_specified", defor33);
    write_field(out, "defor12_specified", defor12);
    write_field(out, "defor13_specified", defor13);
    write_field(out, "defor23_specified", defor23);
    write_field(out, "div_specified", div);
  }

  // ============================================================================
  // TEST 3: calculate_km_kh
  // ============================================================================
  // Case A: km_opt = 1
  {
    DiffusionParams p;
    p.km_opt = 1;
    p.khdif = 12.0;
    p.kvdif = 6.0;

    calculate_km_kh(p, dx, dy,
                    make_const_view3(defor11.data(), nx, nz),
                    make_const_view3(defor22.data(), nx, nz),
                    make_const_view3(defor12.data(), nx + 1, nz),
                    make_const_view3(rdzw.data(), nx, nz),
                    make_const_view3(zx.data(), nx + 1, nz + 1),
                    make_const_view3(zy.data(), nx, nz + 1),
                    make_bview2(msftx.data(), nx),
                    make_bview2(msfty.data(), nx),
                    make_view3(xkmh.data(), nx, nz),
                    make_view3(xkmv.data(), nx, nz),
                    make_view3(xkhh.data(), nx, nz),
                    make_view3(xkhv.data(), nx, nz),
                    nx, ny, nz);

    write_field(out, "xkmh_kmopt1", xkmh);
    write_field(out, "xkmv_kmopt1", xkmv);
    write_field(out, "xkhh_kmopt1", xkhh);
    write_field(out, "xkhv_kmopt1", xkhv);

    copy_dev(xkmh_k1.data(), xkmh.data(), xkmh.size() * sizeof(double));
    copy_dev(xkmv_k1.data(), xkmv.data(), xkmv.size() * sizeof(double));
    copy_dev(xkhh_k1.data(), xkhh.data(), xkhh.size() * sizeof(double));
    copy_dev(xkhv_k1.data(), xkhv.data(), xkhv.size() * sizeof(double));
  }

  // Case B: km_opt = 4 (with diff_opt = 2 slope reduction)
  {
    DiffusionParams p;
    p.km_opt = 4;
    p.diff_opt = 2;
    p.c_s = 0.25;
    p.specified = true;

    calculate_km_kh(p, dx, dy,
                    make_const_view3(defor11.data(), nx, nz),
                    make_const_view3(defor22.data(), nx, nz),
                    make_const_view3(defor12.data(), nx + 1, nz),
                    make_const_view3(rdzw.data(), nx, nz),
                    make_const_view3(zx.data(), nx + 1, nz + 1),
                    make_const_view3(zy.data(), nx, nz + 1),
                    make_bview2(msftx.data(), nx),
                    make_bview2(msfty.data(), nx),
                    make_view3(xkmh.data(), nx, nz),
                    make_view3(xkmv.data(), nx, nz),
                    make_view3(xkhh.data(), nx, nz),
                    make_view3(xkhv.data(), nx, nz),
                    nx, ny, nz);

    write_field(out, "xkmh_kmopt4", xkmh);
    write_field(out, "xkmv_kmopt4", xkmv);
    write_field(out, "xkhh_kmopt4", xkhh);
    write_field(out, "xkhv_kmopt4", xkhv);
  }

  // ============================================================================
  // TEST 4: diff_opt = 1 (coordinate surface diffusion)
  // ============================================================================
  {
    DiffusionParams p;
    p.diff_opt = 1;
    p.specified = true;

    ru_tend.fill(0.0);
    rv_tend.fill(0.0);
    rw_tend.fill(0.0);
    rt_tend.fill(0.0);
    moist_tend.fill(0.0);

    horizontal_diffusion('u', p,
                         make_const_view3(u.data(), nx + 1, nz),
                         make_view3(ru_tend.data(), nx + 1, nz),
                         make_view2(mut.data(), nx),
                         make_bview2(msftx.data(), nx),
                         make_bview2(msfty.data(), nx),
                         make_bview2(msfux.data(), nx + 1),
                         make_bview2(msfuy.data(), nx + 1),
                         make_bview2(msfvx.data(), nx),
                         make_bview2(msfvy.data(), nx),
                         make_bview2(msfvx_inv.data(), nx),
                         make_view1(c1.data()), make_view1(c2.data()),
                         make_const_view3(xkmh.data(), nx, nz),
                         1.0 / dx, 1.0 / dy,
                         nx, ny, nz);

    horizontal_diffusion('v', p,
                         make_const_view3(v.data(), nx, nz),
                         make_view3(rv_tend.data(), nx, nz),
                         make_view2(mut.data(), nx),
                         make_bview2(msftx.data(), nx),
                         make_bview2(msfty.data(), nx),
                         make_bview2(msfux.data(), nx + 1),
                         make_bview2(msfuy.data(), nx + 1),
                         make_bview2(msfvx.data(), nx),
                         make_bview2(msfvy.data(), nx),
                         make_bview2(msfvx_inv.data(), nx),
                         make_view1(c1.data()), make_view1(c2.data()),
                         make_const_view3(xkmh.data(), nx, nz),
                         1.0 / dx, 1.0 / dy,
                         nx, ny, nz);

    horizontal_diffusion('w', p,
                         make_const_view3(w.data(), nx, nz + 1),
                         make_view3(rw_tend.data(), nx, nz + 1),
                         make_view2(mut.data(), nx),
                         make_bview2(msftx.data(), nx),
                         make_bview2(msfty.data(), nx),
                         make_bview2(msfux.data(), nx + 1),
                         make_bview2(msfuy.data(), nx + 1),
                         make_bview2(msfvx.data(), nx),
                         make_bview2(msfvy.data(), nx),
                         make_bview2(msfvx_inv.data(), nx),
                         make_view1(c1f.data()), make_view1(c2f.data()),
                         make_const_view3(xkmh.data(), nx, nz),
                         1.0 / dx, 1.0 / dy,
                         nx, ny, nz);

    horizontal_diffusion('m', p,
                         make_const_view3(theta.data(), nx, nz),
                         make_view3(rt_tend.data(), nx, nz),
                         make_view2(mut.data(), nx),
                         make_bview2(msftx.data(), nx),
                         make_bview2(msfty.data(), nx),
                         make_bview2(msfux.data(), nx + 1),
                         make_bview2(msfuy.data(), nx + 1),
                         make_bview2(msfvx.data(), nx),
                         make_bview2(msfvy.data(), nx),
                         make_bview2(msfvx_inv.data(), nx),
                         make_view1(c1.data()), make_view1(c2.data()),
                         make_const_view3(xkmh.data(), nx, nz),
                         1.0 / dx, 1.0 / dy,
                         nx, ny, nz);

    horizontal_diffusion_3dmp(p,
                              make_const_view3(moist.data(), nx, nz),
                              make_const_view3(dummy3d.data(), nx, nz),
                              make_view3(moist_tend.data(), nx, nz),
                              make_view2(mut.data(), nx),
                              make_bview2(msftx.data(), nx),
                              make_bview2(msfty.data(), nx),
                              make_bview2(msfux.data(), nx + 1),
                              make_bview2(msfuy.data(), nx + 1),
                              make_bview2(msfvy.data(), nx),
                              make_bview2(msfvx_inv.data(), nx),
                              make_view1(c1.data()), make_view1(c2.data()),
                              make_const_view3(xkmh.data(), nx, nz),
                              1.0 / dx, 1.0 / dy,
                              nx, ny, nz);

    write_field(out, "u_hdiff1", ru_tend);
    write_field(out, "v_hdiff1", rv_tend);
    write_field(out, "w_hdiff1", rw_tend);
    write_field(out, "t_hdiff1", rt_tend);
    write_field(out, "m_hdiff1", moist_tend);

    // Vertical diffusion
    ru_tend.fill(0.0);
    rv_tend.fill(0.0);
    rw_tend.fill(0.0);
    rt_tend.fill(0.0);
    moist_tend.fill(0.0);

    vertical_diffusion_u(p,
                         make_const_view3(u.data(), nx + 1, nz),
                         make_view3(ru_tend.data(), nx + 1, nz),
                         make_view1(u_base.data()),
                         make_view2(muu.data(), nx + 1),
                         make_view1(c1.data()), make_view1(c2.data()),
                         make_const_view3(alt.data(), nx, nz),
                         make_view1(rdn.data()), make_view1(rdnw.data()),
                         5.0,
                         nx, ny, nz);

    vertical_diffusion_v(p,
                         make_const_view3(v.data(), nx, nz),
                         make_view3(rv_tend.data(), nx, nz),
                         make_view1(v_base.data()),
                         make_view2(muv.data(), nx),
                         make_view1(c1.data()), make_view1(c2.data()),
                         make_const_view3(alt.data(), nx, nz),
                         make_view1(rdn.data()), make_view1(rdnw.data()),
                         5.0,
                         nx, ny, nz);

    vertical_diffusion('w',
                       make_const_view3(w.data(), nx, nz + 1),
                       make_view3(rw_tend.data(), nx, nz + 1),
                       make_view2(mut.data(), nx),
                       make_view1(c1f.data()), make_view1(c2f.data()),
                       make_const_view3(alt.data(), nx, nz),
                       make_view1(rdn.data()), make_view1(rdnw.data()),
                       5.0,
                       nx, ny, nz);

    vertical_diffusion('m',
                       make_const_view3(theta.data(), nx, nz),
                       make_view3(rt_tend.data(), nx, nz),
                       make_view2(mut.data(), nx),
                       make_view1(c1.data()), make_view1(c2.data()),
                       make_const_view3(alt.data(), nx, nz),
                       make_view1(rdn.data()), make_view1(rdnw.data()),
                       5.0,
                       nx, ny, nz);

    vertical_diffusion_mp(make_const_view3(moist.data(), nx, nz),
                          make_view3(moist_tend.data(), nx, nz),
                          make_const_view3(dummy3d.data(), nx, nz),
                          make_view2(mut.data(), nx),
                          make_view1(c1.data()), make_view1(c2.data()),
                          make_const_view3(alt.data(), nx, nz),
                          make_view1(rdn.data()), make_view1(rdnw.data()),
                          5.0,
                          nx, ny, nz);

    write_field(out, "u_vdiff1", ru_tend);
    write_field(out, "v_vdiff1", rv_tend);
    write_field(out, "w_vdiff1", rw_tend);
    write_field(out, "t_vdiff1", rt_tend);
    write_field(out, "m_vdiff1", moist_tend);
  }

  // ============================================================================
  // TEST 5: diff_opt = 2 (physical diffusion with stress tensor)
  // ============================================================================
  {
    DiffusionParams p;
    p.diff_opt = 2;
    p.km_opt = 4;
    p.specified = true;

    ru_tend.fill(0.0);
    rv_tend.fill(0.0);
    rw_tend.fill(0.0);
    rt_tend.fill(0.0);
    moist_tend.fill(0.0);

    horizontal_diffusion_2(p,
                           make_const_view3(defor11.data(), nx, nz),
                           make_const_view3(defor22.data(), nx, nz),
                           make_const_view3(defor33.data(), nx, nz),
                           make_const_view3(defor12.data(), nx + 1, nz),
                           make_const_view3(defor13.data(), nx + 1, nz + 1),
                           make_const_view3(defor23.data(), nx, nz + 1),
                           make_const_view3(rho.data(), nx, nz),
                           make_const_view3(xkmh.data(), nx, nz),
                           make_const_view3(xkmv.data(), nx, nz),
                           make_const_view3(xkhh.data(), nx, nz),
                           make_const_view3(zx.data(), nx + 1, nz + 1),
                           make_const_view3(zy.data(), nx, nz + 1),
                           make_const_view3(rdz.data(), nx, nz),
                           make_const_view3(rdzw.data(), nx, nz),
                           make_bview2(msftx.data(), nx),
                           make_bview2(msfty.data(), nx),
                           make_bview2(msfux.data(), nx + 1),
                           make_bview2(msfuy.data(), nx + 1),
                           make_bview2(msfvx.data(), nx),
                           make_bview2(msfvy.data(), nx),
                           make_view1(fnm.data()), make_view1(fnp.data()),
                           make_view1(dn.data()), make_view1(dnw.data()),
                           cf1, cf2, cf3, 1.0 / dx, 1.0 / dy,
                           make_view3(titau11.data(), nx, nz),
                           make_view3(titau22.data(), nx, nz),
                           make_view3(titau33.data(), nx, nz),
                           make_view3(titau12.data(), nx + 1, nz),
                           make_view3(titau13.data(), nx + 1, nz + 1),
                           make_view3(titau23.data(), nx, nz + 1),
                           make_view3(ru_tend.data(), nx + 1, nz),
                           make_view3(rv_tend.data(), nx, nz),
                           make_view3(rw_tend.data(), nx, nz + 1),
                           make_const_view3(thp.data(), nx, nz),
                           make_view3(rt_tend.data(), nx, nz),
                           make_const_view3(moist.data(), nx, nz),
                           make_view3(moist_tend.data(), nx, nz),
                           true, nx, ny, nz);

    write_field(out, "u_hdiff2", ru_tend);
    write_field(out, "v_hdiff2", rv_tend);
    write_field(out, "w_hdiff2", rw_tend);
    write_field(out, "t_hdiff2", rt_tend);
    write_field(out, "m_hdiff2", moist_tend);

    // Vertical diffusion 2 with km_opt = 4
    ru_tend.fill(0.0);
    rv_tend.fill(0.0);
    rw_tend.fill(0.0);
    rt_tend.fill(0.0);
    moist_tend.fill(0.0);

    vertical_diffusion_2(p,
                         make_const_view3(titau13.data(), nx + 1, nz + 1),
                         make_const_view3(titau23.data(), nx, nz + 1),
                         make_const_view3(titau33.data(), nx, nz),
                         make_const_view3(thp.data(), nx, nz),
                         make_view1(t_base.data()),
                         make_const_view3(xkhv.data(), nx, nz),
                         make_const_view3(rho.data(), nx, nz),
                         make_const_view3(rdz.data(), nx, nz),
                         make_view1(fnm.data()), make_view1(fnp.data()),
                         make_view1(dn.data()), make_view1(dnw.data()),
                         make_view3(ru_tend.data(), nx + 1, nz),
                         make_view3(rv_tend.data(), nx, nz),
                         make_view3(rw_tend.data(), nx, nz + 1),
                         make_view3(rt_tend.data(), nx, nz),
                         make_const_view3(moist.data(), nx, nz),
                         make_view3(moist_tend.data(), nx, nz),
                         true, nx, ny, nz);

    write_field(out, "u_vdiff2", ru_tend);
    write_field(out, "v_vdiff2", rv_tend);
    write_field(out, "w_vdiff2", rw_tend);
    write_field(out, "t_vdiff2", rt_tend);
    write_field(out, "m_vdiff2", moist_tend);

    // Case C: Vertical diffusion 2 with km_opt = 1 (nonzero xkhv to exercise theta/scalar vertical terms)
    p.km_opt = 1;
    ru_tend.fill(0.0);
    rv_tend.fill(0.0);
    rw_tend.fill(0.0);
    rt_tend.fill(0.0);
    moist_tend.fill(0.0);

    cal_titau_vert(make_const_view3(rho.data(), nx, nz),
                   make_const_view3(xkmh_k1.data(), nx, nz),
                   make_const_view3(xkmv_k1.data(), nx, nz),
                   make_const_view3(defor13.data(), nx + 1, nz + 1),
                   make_const_view3(defor23.data(), nx, nz + 1),
                   make_const_view3(defor33.data(), nx, nz),
                   make_view1(fnm.data()), make_view1(fnp.data()),
                   make_view3(titau13.data(), nx + 1, nz + 1),
                   make_view3(titau23.data(), nx, nz + 1),
                   make_view3(titau33.data(), nx, nz),
                   true, nx, ny, nz);

    vertical_diffusion_2(p,
                         make_const_view3(titau13.data(), nx + 1, nz + 1),
                         make_const_view3(titau23.data(), nx, nz + 1),
                         make_const_view3(titau33.data(), nx, nz),
                         make_const_view3(thp.data(), nx, nz),
                         make_view1(t_base.data()),
                         make_const_view3(xkhv_k1.data(), nx, nz),
                         make_const_view3(rho.data(), nx, nz),
                         make_const_view3(rdz.data(), nx, nz),
                         make_view1(fnm.data()), make_view1(fnp.data()),
                         make_view1(dn.data()), make_view1(dnw.data()),
                         make_view3(ru_tend.data(), nx + 1, nz),
                         make_view3(rv_tend.data(), nx, nz),
                         make_view3(rw_tend.data(), nx, nz + 1),
                         make_view3(rt_tend.data(), nx, nz),
                         make_const_view3(moist.data(), nx, nz),
                         make_view3(moist_tend.data(), nx, nz),
                         true, nx, ny, nz);

    write_field(out, "u_vdiff2_kmopt1", ru_tend);
    write_field(out, "v_vdiff2_kmopt1", rv_tend);
    write_field(out, "w_vdiff2_kmopt1", rw_tend);
    write_field(out, "t_vdiff2_kmopt1", rt_tend);
    write_field(out, "m_vdiff2_kmopt1", moist_tend);
  }

  // ============================================================================
  // TEST 6: sixth_order_diffusion
  // ============================================================================
  // Case A: diff_6th_opt = 1
  ru_tend.fill(0.0);
  rv_tend.fill(0.0);
  rw_tend.fill(0.0);
  rt_tend.fill(0.0);
  moist_tend.fill(0.0);

  {
    DiffusionParams p;
    p.diff_6th_opt = 1;
    p.diff_6th_factor = 0.12;
    p.specified = true;

    sixth_order_diffusion('u', p, make_const_view3(u.data(), nx + 1, nz),
                          make_view3(ru_tend.data(), nx + 1, nz),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1.data()), make_view1(c2.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    sixth_order_diffusion('v', p, make_const_view3(v.data(), nx, nz),
                          make_view3(rv_tend.data(), nx, nz),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1.data()), make_view1(c2.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    sixth_order_diffusion('w', p, make_const_view3(w.data(), nx, nz + 1),
                          make_view3(rw_tend.data(), nx, nz + 1),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1f.data()), make_view1(c2f.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    sixth_order_diffusion('m', p, make_const_view3(theta.data(), nx, nz),
                          make_view3(rt_tend.data(), nx, nz),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1.data()), make_view1(c2.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    sixth_order_diffusion('m', p, make_const_view3(moist.data(), nx, nz),
                          make_view3(moist_tend.data(), nx, nz),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1.data()), make_view1(c2.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    write_field(out, "u_6th_opt1", ru_tend);
    write_field(out, "v_6th_opt1", rv_tend);
    write_field(out, "w_6th_opt1", rw_tend);
    write_field(out, "t_6th_opt1", rt_tend);
    write_field(out, "m_6th_opt1", moist_tend);
  }

  // Case B: diff_6th_opt = 2
  ru_tend.fill(0.0);
  rv_tend.fill(0.0);
  rw_tend.fill(0.0);
  rt_tend.fill(0.0);
  moist_tend.fill(0.0);

  {
    DiffusionParams p;
    p.diff_6th_opt = 2;
    p.diff_6th_factor = 0.12;
    p.specified = true;

    sixth_order_diffusion('u', p, make_const_view3(u.data(), nx + 1, nz),
                          make_view3(ru_tend.data(), nx + 1, nz),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1.data()), make_view1(c2.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    sixth_order_diffusion('v', p, make_const_view3(v.data(), nx, nz),
                          make_view3(rv_tend.data(), nx, nz),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1.data()), make_view1(c2.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    sixth_order_diffusion('w', p, make_const_view3(w.data(), nx, nz + 1),
                          make_view3(rw_tend.data(), nx, nz + 1),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1f.data()), make_view1(c2f.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    sixth_order_diffusion('m', p, make_const_view3(theta.data(), nx, nz),
                          make_view3(rt_tend.data(), nx, nz),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1.data()), make_view1(c2.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    sixth_order_diffusion('m', p, make_const_view3(moist.data(), nx, nz),
                          make_view3(moist_tend.data(), nx, nz),
                          make_view2(mut.data(), nx),
                          dt,
                          make_view1(c1.data()), make_view1(c2.data()),
                          make_bview2(msftx.data(), nx),
                          make_bview2(msfty.data(), nx),
                          make_bview2(msfux.data(), nx + 1),
                          make_bview2(msfuy.data(), nx + 1),
                          make_bview2(msfvx.data(), nx),
                          make_bview2(msfvy.data(), nx),
                          nx, ny, nz);

    write_field(out, "u_6th_opt2", ru_tend);
    write_field(out, "v_6th_opt2", rv_tend);
    write_field(out, "w_6th_opt2", rw_tend);
    write_field(out, "t_6th_opt2", rt_tend);
    write_field(out, "m_6th_opt2", moist_tend);
  }

  // ============================================================================
  // TEST 7: w_damp
  // ============================================================================
  rw_tend.fill(0.0);
  {
    DiffusionParams p;
    p.w_damping = 1;
    p.w_crit_cfl = 1.0;
    p.zadvect_implicit = 0;
    p.polar = false;

    double max_vert_cfl = 0.0;
    double max_horiz_cfl = 0.0;

    w_damp(make_view3(rw_tend.data(), nx, nz + 1),
           max_vert_cfl, max_horiz_cfl,
           make_const_view3(u.data(), nx + 1, nz),
           make_const_view3(v.data(), nx, nz),
           make_const_view3(ww.data(), nx, nz + 1),
           make_const_view3(w.data(), nx, nz + 1),
           make_view2(mut.data(), nx),
           make_view1(c1f.data()), make_view1(c2f.data()),
           make_view1(rdnw.data()),
           1.0 / dx, 1.0 / dy,
           make_bview2(msfux.data(), nx + 1),
           make_bview2(msfuy.data(), nx + 1),
           make_bview2(msfvx.data(), nx),
           make_bview2(msfvy.data(), nx),
           dt, p,
           nx, ny, nz);

    write_field(out, "rw_tend_damp", rw_tend);
    write_chunk(out, "max_vert_cfl", &max_vert_cfl, 1);
    write_chunk(out, "max_horiz_cfl", &max_horiz_cfl, 1);
  }

  out.close();
  std::cout << "C++ driver finished successfully." << std::endl;
  return 0;
}
