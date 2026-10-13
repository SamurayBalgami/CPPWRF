#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/grid.hpp"
#include "core/halo.hpp"
#include "dyn/state.hpp"
#include "phys/driver/phy_prep.hpp"
#include "phys/driver/tendency_coupling.hpp"

using namespace cppwrf;
using namespace cppwrf::dyn;
using namespace cppwrf::phys::driver;

template <class T>
void read_flat(FILE* f, T* buf, std::size_t n) {
  if (std::fread(buf, sizeof(T), n, f) != n) {
    std::fprintf(stderr, "Unexpected EOF reading %zu elements\n", n);
    std::exit(1);
  }
}

template <class T>
void write_flat(FILE* f, const T* buf, std::size_t n) {
  if (std::fwrite(buf, sizeof(T), n, f) != n) {
    std::fprintf(stderr, "Failed writing %zu elements\n", n);
    std::exit(1);
  }
}

// Reads 3D field stored in (nz, ny, nx) row-major order into Field3
void read_field3_m(FILE* f, Field3& fld, int nx, int ny, int nz) {
  std::vector<real> grid(nz * ny * nx);
  std::vector<double> buf(nx);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      read_flat(f, buf.data(), nx);
      for (int i = 0; i < nx; ++i) grid[k * ny * nx + j * nx + i] = real(buf[i]);
    }
  }
  fld.fill_host([&](int i, int j, int k) {
    return grid[k * ny * nx + j * nx + i];
  });
}

// Writes Field3 in (nz, ny, nx) order
void write_field3_m(FILE* f, const Field3& fld, int nx, int ny, int nz) {
  std::vector<real> h;
  fld.download(h);
  const auto& L = fld.layout();
  std::vector<double> buf(nx);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) buf[i] = double(h[L.index(i, j, k)]);
      write_flat(f, buf.data(), nx);
    }
  }
}

// Reads U-staggered 3D field (nz, ny, nx+1)
void read_field3_u(FILE* f, Field3& fld, int nx, int ny, int nz) {
  const int nx_u = nx + 1;
  std::vector<real> grid(nz * ny * nx_u);
  std::vector<double> buf(nx_u);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      read_flat(f, buf.data(), nx_u);
      for (int i = 0; i < nx_u; ++i) grid[k * ny * nx_u + j * nx_u + i] = real(buf[i]);
    }
  }
  fld.fill_host([&](int i, int j, int k) {
    return grid[k * ny * nx_u + j * nx_u + i];
  });
}

void write_field3_u(FILE* f, const Field3& fld, int nx, int ny, int nz) {
  std::vector<real> h;
  fld.download(h);
  const auto& L = fld.layout();
  const int nx_u = nx + 1;
  std::vector<double> buf(nx_u);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx_u; ++i) buf[i] = double(h[L.index(i, j, k)]);
      write_flat(f, buf.data(), nx_u);
    }
  }
}

// Reads V-staggered 3D field (nz, ny+1, nx)
void read_field3_v(FILE* f, Field3& fld, int nx, int ny, int nz) {
  const int ny_v = ny + 1;
  std::vector<real> grid(nz * ny_v * nx);
  std::vector<double> buf(nx);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny_v; ++j) {
      read_flat(f, buf.data(), nx);
      for (int i = 0; i < nx; ++i) grid[k * ny_v * nx + j * nx + i] = real(buf[i]);
    }
  }
  fld.fill_host([&](int i, int j, int k) {
    return grid[k * ny_v * nx + j * nx + i];
  });
}

void write_field3_v(FILE* f, const Field3& fld, int nx, int ny, int nz) {
  std::vector<real> h;
  fld.download(h);
  const auto& L = fld.layout();
  const int ny_v = ny + 1;
  std::vector<double> buf(nx);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny_v; ++j) {
      for (int i = 0; i < nx; ++i) buf[i] = double(h[L.index(i, j, k)]);
      write_flat(f, buf.data(), nx);
    }
  }
}

// Reads 2D field (ny, nx)
void read_field2_m(FILE* f, Field2& fld, int nx, int ny) {
  std::vector<real> grid(ny * nx);
  std::vector<double> buf(nx);
  for (int j = 0; j < ny; ++j) {
    read_flat(f, buf.data(), nx);
    for (int i = 0; i < nx; ++i) grid[j * nx + i] = real(buf[i]);
  }
  fld.fill_host([&](int i, int j, int) {
    return grid[j * nx + i];
  });
}

void read_field2_u(FILE* f, Field2& fld, int nx, int ny) {
  const int nx_u = nx + 1;
  std::vector<real> grid(ny * nx_u);
  std::vector<double> buf(nx_u);
  for (int j = 0; j < ny; ++j) {
    read_flat(f, buf.data(), nx_u);
    for (int i = 0; i < nx_u; ++i) grid[j * nx_u + i] = real(buf[i]);
  }
  fld.fill_host([&](int i, int j, int) {
    return grid[j * nx_u + i];
  });
}

void read_field2_v(FILE* f, Field2& fld, int nx, int ny) {
  const int ny_v = ny + 1;
  std::vector<real> grid(ny_v * nx);
  std::vector<double> buf(nx);
  for (int j = 0; j < ny_v; ++j) {
    read_flat(f, buf.data(), nx);
    for (int i = 0; i < nx; ++i) grid[j * nx + i] = real(buf[i]);
  }
  fld.fill_host([&](int i, int j, int) {
    return grid[j * nx + i];
  });
}

// Reads 2D BaseField2 (ny, nx) with halo
void read_basefield2_m(FILE* f, BaseField2& fld, int nx, int ny) {
  std::vector<real> grid(ny * nx);
  std::vector<double> buf(nx);
  for (int j = 0; j < ny; ++j) {
    read_flat(f, buf.data(), nx);
    for (int i = 0; i < nx; ++i) grid[j * nx + i] = real(buf[i]);
  }
  fld.fill_host_halo([&](int i, int j, int) {
    int ii = ((i % nx) + nx) % nx;
    int jj = ((j % ny) + ny) % ny;
    return grid[jj * nx + ii];
  });
}

void read_basefield2_u(FILE* f, BaseField2& fld, int nx, int ny) {
  const int nx_u = nx + 1;
  std::vector<real> grid(ny * nx_u);
  std::vector<double> buf(nx_u);
  for (int j = 0; j < ny; ++j) {
    read_flat(f, buf.data(), nx_u);
    for (int i = 0; i < nx_u; ++i) grid[j * nx_u + i] = real(buf[i]);
  }
  fld.fill_host_halo([&](int i, int j, int) {
    int ii = ((i % nx) + nx) % nx;
    if (i == nx) ii = nx;
    int jj = ((j % ny) + ny) % ny;
    return grid[jj * nx_u + ii];
  });
}

void read_basefield2_v(FILE* f, BaseField2& fld, int nx, int ny) {
  const int ny_v = ny + 1;
  std::vector<real> grid(ny_v * nx);
  std::vector<double> buf(nx);
  for (int j = 0; j < ny_v; ++j) {
    read_flat(f, buf.data(), nx);
    for (int i = 0; i < nx; ++i) grid[j * nx + i] = real(buf[i]);
  }
  fld.fill_host_halo([&](int i, int j, int) {
    int ii = ((i % nx) + nx) % nx;
    int jj = ((j % ny) + ny) % ny;
    if (j == ny) jj = ny;
    return grid[jj * nx + ii];
  });
}

int main(int argc, char** argv) {
  std::string infile = "driver_input.bin";
  std::string outfile = "driver_output_cpp.bin";
  if (argc >= 2) infile = argv[1];
  if (argc >= 3) outfile = argv[2];

  FILE* f = std::fopen(infile.c_str(), "rb");
  if (!f) {
    std::fprintf(stderr, "Cannot open %s\n", infile.c_str());
    return 1;
  }

  int32_t nx = 0, ny = 0, nz = 0, n_moist = 0, use_theta_m = 1;
  double p_top = 0.0;
  read_flat(f, &nx, 1);
  read_flat(f, &ny, 1);
  read_flat(f, &nz, 1);
  read_flat(f, &n_moist, 1);
  read_flat(f, &p_top, 1);
  read_flat(f, &use_theta_m, 1);

  Grid g{nx, ny, nz};
  g.p_top = p_top;
  g.bc_x = BC::Periodic;
  g.bc_y = BC::Periodic;

  // BaseState: non-flat (terrain present), mapped (map factors present)
  BaseState b{g, false, true};
  b.p_top = p_top;

  // 1D vertical arrays
  std::vector<double> znw(nz + 1), znu(nz), dnw(nz), dn(nz + 1), fnm(nz + 1), fnp(nz + 1);
  std::vector<double> c1h(nz), c2h(nz), c1f(nz + 1), c2f(nz + 1);
  read_flat(f, znw.data(), nz + 1);
  read_flat(f, znu.data(), nz);
  read_flat(f, dnw.data(), nz);
  read_flat(f, dn.data(), nz + 1);
  read_flat(f, fnm.data(), nz + 1);
  read_flat(f, fnp.data(), nz + 1);
  read_flat(f, c1h.data(), nz);
  read_flat(f, c2h.data(), nz);
  read_flat(f, c1f.data(), nz + 1);
  read_flat(f, c2f.data(), nz + 1);

  for (int k = 0; k <= nz; ++k) {
    b.znw[k] = znw[k];
    b.dn[k] = dn[k];
    b.fnm[k] = fnm[k];
    b.fnp[k] = fnp[k];
    b.c1f[k] = c1f[k];
    b.c2f[k] = c2f[k];
  }
  for (int k = 0; k < nz; ++k) {
    b.znu[k] = znu[k];
    b.dnw[k] = dnw[k];
    b.c1h[k] = c1h[k];
    b.c2h[k] = c2h[k];
  }
  b.upload_levels();

  // State
  State s{g};
  read_field2_m(f, s.mut, nx, ny);
  read_field2_u(f, s.muu, nx, ny);
  read_field2_v(f, s.muv, nx, ny);

  // 2D BaseState map factors and terrain
  read_basefield2_m(f, b.msftx, nx, ny);
  read_basefield2_m(f, b.msfty, nx, ny);
  read_basefield2_u(f, b.msfux, nx, ny);
  read_basefield2_u(f, b.msfuy, nx, ny);
  read_basefield2_v(f, b.msfvx, nx, ny);
  read_basefield2_v(f, b.msfvy, nx, ny);
  read_basefield2_v(f, b.msfvx_inv, nx, ny);
  read_basefield2_m(f, b.ht, nx, ny);

  // 3D prognostics and diagnostics
  read_field3_u(f, s.u_2, nx, ny, nz);
  read_field3_v(f, s.v_2, nx, ny, nz);
  read_field3_m(f, s.t_2, nx, ny, nz);
  read_field3_m(f, s.p, nx, ny, nz);

  // pb
  {
    std::vector<real> pb_grid(nz * ny * nx);
    std::vector<double> buf(nx);
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        read_flat(f, buf.data(), nx);
        for (int i = 0; i < nx; ++i) pb_grid[k * ny * nx + j * nx + i] = real(buf[i]);
      }
    }
    b.pb.fill_host_halo([&](int i, int j, int k) {
      int ii = ((i % nx) + nx) % nx;
      int jj = ((j % ny) + ny) % ny;
      return pb_grid[k * ny * nx + jj * nx + ii];
    });
  }

  read_field3_m(f, s.alt, nx, ny, nz);
  read_field3_m(f, s.ph_2, nx, ny, nz + 1);

  // phb
  {
    std::vector<real> phb_grid((nz + 1) * ny * nx);
    std::vector<double> buf(nx);
    for (int k = 0; k < nz + 1; ++k) {
      for (int j = 0; j < ny; ++j) {
        read_flat(f, buf.data(), nx);
        for (int i = 0; i < nx; ++i) phb_grid[k * ny * nx + j * nx + i] = real(buf[i]);
      }
    }
    b.phb.fill_host_halo([&](int i, int j, int k) {
      int ii = ((i % nx) + nx) % nx;
      int jj = ((j % ny) + ny) % ny;
      return phb_grid[k * ny * nx + jj * nx + ii];
    });
  }

  // Register moist species
  s.add_scalar(g, "qv", true);
  s.add_scalar(g, "qc", true);
  s.add_scalar(g, "qr", true);
  s.add_scalar(g, "qi", true);
  s.add_scalar(g, "qs", true);
  s.add_scalar(g, "qg", true);
  s.P_QV = 0; s.P_QC = 1; s.P_QR = 2; s.P_QI = 3; s.P_QS = 4; s.P_QG = 5;

  for (int m = 0; m < n_moist; ++m) {
    read_field3_m(f, s.q_2[m], nx, ny, nz);
    fill_halo(s.q_2[m], g);
  }

  // qtot = sum of moist species
  {
    std::vector<real> qtot_grid(nz * ny * nx, real(0));
    for (int m = 0; m < n_moist; ++m) {
      std::vector<real> q_h;
      s.q_2[m].download(q_h);
      const auto& L = s.q_2[m].layout();
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
          for (int i = 0; i < nx; ++i)
            qtot_grid[k * ny * nx + j * nx + i] += q_h[L.index(i, j, k)];
    }
    s.qtot.fill_host([&](int i, int j, int k) {
      return qtot_grid[k * ny * nx + j * nx + i];
    });
    fill_halo(s.qtot, g);
  }

  fill_halo(s.u_2, g);
  fill_halo(s.v_2, g);
  fill_halo(s.t_2, g);
  fill_halo(s.p, g);
  fill_halo(s.alt, g);
  fill_halo(s.ph_2, g);
  fill_halo(s.mut, g);
  fill_halo(s.muu, g);
  fill_halo(s.muv, g);

  // Read physics tendencies
  PhysicsTendencies pt{g};
  read_field3_m(f, pt.rthraten, nx, ny, nz);
  read_field3_m(f, pt.rthblten, nx, ny, nz);
  read_field3_m(f, pt.rublten, nx, ny, nz);
  read_field3_m(f, pt.rvblten, nx, ny, nz);
  read_field3_m(f, pt.rqvblten, nx, ny, nz);
  read_field3_m(f, pt.rqcblten, nx, ny, nz);
  read_field3_m(f, pt.rqiblten, nx, ny, nz);
  read_field3_m(f, pt.rthcuten, nx, ny, nz);
  read_field3_m(f, pt.rucuten, nx, ny, nz);
  read_field3_m(f, pt.rvcuten, nx, ny, nz);
  read_field3_m(f, pt.rqvcuten, nx, ny, nz);
  read_field3_m(f, pt.rqccuten, nx, ny, nz);
  read_field3_m(f, pt.rqrcuten, nx, ny, nz);
  read_field3_m(f, pt.rqicuten, nx, ny, nz);
  read_field3_m(f, pt.rqscuten, nx, ny, nz);
  read_field3_m(f, pt.rthshten, nx, ny, nz);
  read_field3_m(f, pt.rushten, nx, ny, nz);
  read_field3_m(f, pt.rvshten, nx, ny, nz);
  read_field3_m(f, pt.rqvshten, nx, ny, nz);
  read_field3_m(f, pt.rqcshten, nx, ny, nz);
  read_field3_m(f, pt.rqrshten, nx, ny, nz);
  read_field3_m(f, pt.rqishten, nx, ny, nz);
  read_field3_m(f, pt.rqsshten, nx, ny, nz);
  read_field3_m(f, pt.rqgshten, nx, ny, nz);

  std::fclose(f);

  // Tendencies container
  Tendencies t{g};
  t.add_scalar(g, "qv");
  t.add_scalar(g, "qc");
  t.add_scalar(g, "qr");
  t.add_scalar(g, "qi");
  t.add_scalar(g, "qs");
  t.add_scalar(g, "qg");

  // 1. Run phy_prep
  PhyPrepFields pp{g};
  run_phy_prep(s, b, pp, PhyPrepParams{use_theta_m});

  // 2. Run Tendency Coupling (calculate_phy_tend, update_phy_ten, apply_tendf_map_factors)
  couple_physics_tendencies(t, pt, s, b, g, true);

  // Write C++ output binary
  FILE* fo = std::fopen(outfile.c_str(), "wb");
  if (!fo) {
    std::fprintf(stderr, "Cannot open %s for writing\n", outfile.c_str());
    return 1;
  }

  // phy_prep outputs
  write_field3_m(fo, pp.th_phy, nx, ny, nz);
  write_field3_m(fo, pp.th_phy_m_t0, nx, ny, nz);
  write_field3_m(fo, pp.p_phy, nx, ny, nz);
  write_field3_m(fo, pp.pi_phy, nx, ny, nz);
  write_field3_m(fo, pp.t_phy, nx, ny, nz);
  write_field3_m(fo, pp.rho, nx, ny, nz);
  write_field3_m(fo, pp.u_phy, nx, ny, nz);
  write_field3_m(fo, pp.v_phy, nx, ny, nz);
  write_field3_m(fo, pp.z, nx, ny, nz);
  write_field3_m(fo, pp.dz8w, nx, ny, nz + 1);
  write_field3_m(fo, pp.z_at_w, nx, ny, nz + 1);
  write_field3_m(fo, pp.p8w, nx, ny, nz + 1);
  write_field3_m(fo, pp.t8w, nx, ny, nz + 1);
  write_field3_m(fo, pp.p_hyd, nx, ny, nz);
  write_field3_m(fo, pp.p_hyd_w, nx, ny, nz + 1);

  // Tendency coupling outputs
  write_field3_u(fo, t.ru_tendf, nx, ny, nz);
  write_field3_v(fo, t.rv_tendf, nx, ny, nz);
  write_field3_m(fo, t.t_tendf, nx, ny, nz);
  for (int m = 0; m < n_moist; ++m) {
    write_field3_m(fo, t.q_tendf[m], nx, ny, nz);
  }

  write_field3_u(fo, t.ru_tend, nx, ny, nz);
  write_field3_v(fo, t.rv_tend, nx, ny, nz);
  write_field3_m(fo, t.t_tend, nx, ny, nz);

  std::fclose(fo);

  return 0;
}
