#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdint>
#include <cassert>

namespace cppwrf::nest::test {

struct CoarseDataset {
  int nx = 0, ny = 0, nz = 0;
  double p_top = 0.0;

  std::vector<double> u;      // (nx+1, nz, ny)
  std::vector<double> v;      // (nx, nz, ny+1)
  std::vector<double> w;      // (nx, nz+1, ny)
  std::vector<double> ph;     // (nx, nz+1, ny)
  std::vector<double> t;      // (nx, nz, ny)
  std::vector<double> mu;     // (nx, ny)
  std::vector<double> qvapor; // (nx, nz, ny)
  std::vector<double> hgt;    // (nx, ny)
  std::vector<double> mub;    // (nx, ny)
  std::vector<double> phb;    // (nx, nz+1, ny)

  std::vector<double> msfuy;  // (nx+1, ny)
  std::vector<double> msfvx;  // (nx, ny+1)
  std::vector<double> msfty;  // (nx, ny)

  std::vector<double> c1h, c2h, c3h, c4h; // nz
  std::vector<double> c1f, c2f;           // nz+1

  bool load(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
      std::cerr << "Error: cannot open binary dataset at " << path << std::endl;
      return false;
    }

    f.read(reinterpret_cast<char*>(&nx), sizeof(int));
    f.read(reinterpret_cast<char*>(&ny), sizeof(int));
    f.read(reinterpret_cast<char*>(&nz), sizeof(int));
    f.read(reinterpret_cast<char*>(&p_top), sizeof(double));

    auto read_vec = [&](std::vector<double>& v, std::size_t sz) {
      v.resize(sz);
      f.read(reinterpret_cast<char*>(v.data()), sz * sizeof(double));
    };

    read_vec(u, static_cast<std::size_t>((nx + 1) * nz * ny));
    read_vec(v, static_cast<std::size_t>(nx * nz * (ny + 1)));
    read_vec(w, static_cast<std::size_t>(nx * (nz + 1) * ny));
    read_vec(ph, static_cast<std::size_t>(nx * (nz + 1) * ny));
    read_vec(t, static_cast<std::size_t>(nx * nz * ny));
    read_vec(mu, static_cast<std::size_t>(nx * ny));
    read_vec(qvapor, static_cast<std::size_t>(nx * nz * ny));
    read_vec(hgt, static_cast<std::size_t>(nx * ny));
    read_vec(mub, static_cast<std::size_t>(nx * ny));
    read_vec(phb, static_cast<std::size_t>(nx * (nz + 1) * ny));

    read_vec(msfuy, static_cast<std::size_t>((nx + 1) * ny));
    read_vec(msfvx, static_cast<std::size_t>(nx * (ny + 1)));
    read_vec(msfty, static_cast<std::size_t>(nx * ny));

    read_vec(c1h, static_cast<std::size_t>(nz));
    read_vec(c2h, static_cast<std::size_t>(nz));
    read_vec(c3h, static_cast<std::size_t>(nz));
    read_vec(c4h, static_cast<std::size_t>(nz));
    read_vec(c1f, static_cast<std::size_t>(nz + 1));
    read_vec(c2f, static_cast<std::size_t>(nz + 1));

    return true;
  }
};

}  // namespace cppwrf::nest::test
