#pragma once
// wrfinput_d01 / wrfbdy_d01 writers (WRF names, dimensions and per-variable attributes).
#include <map>
#include <string>

#include "init_real/met_em.hpp"
#include "init_real/nc_out.hpp"
#include "init_real/real_bdy.hpp"
#include "init_real/real_domain.hpp"

namespace cppwrf::init_real {

template <typename T>
void write_wrfinput(const RealDomain<T>& dm, const MetEm& m, const std::string& path, const std::string& date) {
  (void)m;
  const int nx = dm.d.ide - 1, ny = dm.d.jde - 1, nz = dm.d.kde - 1;
  const char ft = sizeof(T) == 8 ? 'd' : 'f';
  const int fieldtype = sizeof(T) == 8 ? 105 : 104;
  NcOut nc(path);
  nc.dim("Time", 0);
  nc.dim("DateStrLen", 19);
  nc.dim("west_east", std::size_t(nx));
  nc.dim("south_north", std::size_t(ny));
  nc.dim("bottom_top", std::size_t(nz));
  nc.dim("bottom_top_stag", std::size_t(nz + 1));
  nc.dim("soil_layers_stag", std::size_t(dm.o.num_soil_layers));
  nc.dim("west_east_stag", std::size_t(nx + 1));
  nc.dim("south_north_stag", std::size_t(ny + 1));
  nc.var("Times", 'c', {"Time", "DateStrLen"}, {}, 0);
  struct V3 { const char* name; const F3<T>* f; char st; };
  struct V2 { const char* name; const F2<T>* f; char st; };
  struct V1 { const char* name; const F1<T>* f; bool stag; };
  const V3 v3[] = {{"U", &dm.u_2, 'X'}, {"V", &dm.v_2, 'Y'}, {"PH", &dm.ph_2, 'Z'}, {"PHB", &dm.phb, 'Z'},
                   {"T", &dm.th_phy_m_t0, ' '}, {"THM", &dm.t_2, ' '}, {"T_INIT", &dm.t_init, ' '}, {"P", &dm.p, ' '},
                   {"AL", &dm.al, ' '}, {"ALB", &dm.alb, ' '}, {"PB", &dm.pb, ' '}, {"P_HYD", &dm.p_hyd, ' '},
                   {"QVAPOR", &dm.qv, ' '}};
  const V2 v2[] = {{"MU", &dm.mu_2, ' '}, {"MUB", &dm.mub, ' '}, {"PSFC", &dm.psfc, ' '}, {"U10", &dm.u10, ' '},
                   {"V10", &dm.v10, ' '}, {"T2", &dm.t2, ' '}, {"TH2", &dm.th2, ' '}, {"Q2", &dm.q2, ' '},
                   {"HGT", &dm.ht, ' '}, {"TSK", &dm.tsk, ' '}, {"TMN", &dm.tmn, ' '}, {"LANDMASK", &dm.landmask, ' '},
                   {"XLAND", &dm.xland, ' '}, {"LU_INDEX", &dm.lu_index, ' '}, {"VEGFRA", &dm.vegfra, ' '},
                   {"SHDMAX", &dm.shdmax, ' '}, {"SHDMIN", &dm.shdmin, ' '}, {"SHDAVG", &dm.shdavg, ' '},
                   {"ALBBCK", &dm.albbck, ' '}, {"SNOALB", &dm.snoalb, ' '}, {"LAI", &dm.lai, ' '},
                   {"SNOW", &dm.snow, ' '}, {"SNOWH", &dm.snowh, ' '}, {"SST", &dm.sst, ' '},
                   {"SEAICE", &dm.xice, ' '}, {"LAKEMASK", &dm.lakemask, ' '}, {"XLAT", &dm.xlat, ' '},
                   {"XLONG", &dm.xlong, ' '}, {"MF_VX_INV", &dm.msfvx_inv, 'Y'}, {"WATER_DEPTH", &dm.water_depth, ' '}};
  const V1 v1[] = {{"ZNW", &dm.znw, true}, {"ZNU", &dm.znu, false}, {"DNW", &dm.dnw, false}, {"RDNW", &dm.rdnw, false},
                   {"DN", &dm.dn, false}, {"RDN", &dm.rdn, false}, {"FNM", &dm.fnm, false}, {"FNP", &dm.fnp, false},
                   {"C1H", &dm.c1h, false}, {"C2H", &dm.c2h, false}, {"C3H", &dm.c3h, false}, {"C4H", &dm.c4h, false},
                   {"C1F", &dm.c1f, true}, {"C2F", &dm.c2f, true}, {"C3F", &dm.c3f, true}, {"C4F", &dm.c4f, true}};
  for (const V3& v : v3) {
    std::string zd = v.st == 'Z' ? "bottom_top_stag" : "bottom_top";
    std::string yd = v.st == 'Y' ? "south_north_stag" : "south_north";
    std::string xd = v.st == 'X' ? "west_east_stag" : "west_east";
    nc.var(v.name, ft, {"Time", zd, yd, xd}, {{"MemoryOrder", "XYZ"}}, fieldtype);
  }
  for (const V2& v : v2) {
    std::string yd = v.st == 'Y' ? "south_north_stag" : "south_north";
    nc.var(v.name, ft, {"Time", yd, "west_east"}, {{"MemoryOrder", "XY "}}, fieldtype);
  }
  for (const V1& v : v1) nc.var(v.name, ft, {"Time", v.stag ? "bottom_top_stag" : "bottom_top"}, {{"MemoryOrder", "Z  "}}, fieldtype);
  nc.var("SMOIS", ft, {"Time", "soil_layers_stag", "south_north", "west_east"}, {{"MemoryOrder", "XYZ"}}, fieldtype);
  nc.var("TSLB", ft, {"Time", "soil_layers_stag", "south_north", "west_east"}, {{"MemoryOrder", "XYZ"}}, fieldtype);
  nc.var("SH2O", ft, {"Time", "soil_layers_stag", "south_north", "west_east"}, {{"MemoryOrder", "XYZ"}}, fieldtype);
  nc.var("ZS", ft, {"Time", "soil_layers_stag"}, {{"MemoryOrder", "Z  "}}, fieldtype);
  nc.var("DZS", ft, {"Time", "soil_layers_stag"}, {{"MemoryOrder", "Z  "}}, fieldtype);
  nc.var("IVGTYP", 'i', {"Time", "south_north", "west_east"}, {{"MemoryOrder", "XY "}}, 106);
  nc.var("ISLTYP", 'i', {"Time", "south_north", "west_east"}, {{"MemoryOrder", "XY "}}, 106);
  for (const char* s : {"P_TOP", "CF1", "CF2", "CF3", "CFN", "CFN1", "RDX", "RDY", "T00", "P00", "TLP", "TISO", "TLP_STRAT",
                        "P_STRAT"})
    nc.var(s, ft, {"Time"}, {{"MemoryOrder", "0  "}}, fieldtype);
  nc.enddef();
  nc.put_text("Times", date);
  for (const V3& v : v3) {
    int ni = v.st == 'X' ? nx + 1 : nx, nj = v.st == 'Y' ? ny + 1 : ny, nk = v.st == 'Z' ? nz + 1 : nz;
    nc.put(v.name, nc_order(*v.f, ni, nk, nj));
  }
  for (const V2& v : v2) nc.put(v.name, nc_order(*v.f, nx, v.st == 'Y' ? ny + 1 : ny));
  for (const V1& v : v1) nc.put(v.name, nc_order(*v.f, v.stag ? nz + 1 : nz));
  nc.put("SMOIS", nc_order(dm.smois, nx, dm.o.num_soil_layers, ny));
  nc.put("TSLB", nc_order(dm.tslb, nx, dm.o.num_soil_layers, ny));
  nc.put("SH2O", nc_order(dm.sh2o, nx, dm.o.num_soil_layers, ny));
  nc.put("ZS", nc_order(dm.zs, dm.o.num_soil_layers));
  nc.put("DZS", nc_order(dm.dzs, dm.o.num_soil_layers));
  {
    std::vector<int> a(static_cast<std::size_t>(nx) * ny), b(static_cast<std::size_t>(nx) * ny);
    for (int j = 1; j <= ny; ++j)
      for (int i = 1; i <= nx; ++i) {
        a[std::size_t(j - 1) * nx + (i - 1)] = dm.ivgtyp(i, j);
        b[std::size_t(j - 1) * nx + (i - 1)] = dm.isltyp(i, j);
      }
    nc.put_int("IVGTYP", a);
    nc.put_int("ISLTYP", b);
  }
  auto s = [&](const char* n, T v) { nc.put(n, {double(v)}); };
  s("P_TOP", dm.p_top); s("CF1", dm.cf1); s("CF2", dm.cf2); s("CF3", dm.cf3); s("CFN", dm.cfn); s("CFN1", dm.cfn1);
  s("RDX", dm.rdx); s("RDY", dm.rdy); s("T00", dm.t00); s("P00", dm.p00); s("TLP", dm.a); s("TISO", dm.tiso);
  s("TLP_STRAT", dm.a_strat); s("P_STRAT", dm.p_strat);
}

}  // namespace cppwrf::init_real

namespace cppwrf::init_real {

// wrfbdy_d01: one record per boundary interval (values at the start time of the interval and
// tendencies over it).
template <typename T>
class WrfBdyWriter {
 public:
  struct Vars {  // values and tendencies of one field
    const BdyVar<T>* val;
    const BdyVar<T>* tend;
  };
  WrfBdyWriter(const std::string& path, const Dims& d, int w) : nc_(path), d_(d), w_(w) {
    const int nx = d.ide - 1, ny = d.jde - 1, nz = d.kde - 1;
    nc_.dim("Time", 0);
    nc_.dim("DateStrLen", 19);
    nc_.dim("south_north", std::size_t(ny));
    nc_.dim("bottom_top", std::size_t(nz));
    nc_.dim("bdy_width", std::size_t(w));
    nc_.dim("west_east_stag", std::size_t(nx + 1));
    nc_.dim("south_north_stag", std::size_t(ny + 1));
    nc_.dim("west_east", std::size_t(nx));
    nc_.dim("bottom_top_stag", std::size_t(nz + 1));
    nc_.var("Times", 'c', {"Time", "DateStrLen"}, {}, 0);
    nc_.var("md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", 'c', {"Time", "DateStrLen"}, {}, 0);
    nc_.var("md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", 'c', {"Time", "DateStrLen"}, {}, 0);
    const char ft = sizeof(T) == 8 ? 'd' : 'f';
    const int fieldtype = sizeof(T) == 8 ? 105 : 104;
    for (const Spec& s : specs()) {
      for (const char* side : {"BXS", "BXE", "BYS", "BYE", "BTXS", "BTXE", "BTYS", "BTYE"}) {
        const bool x = side[1] == 'X' || side[2] == 'X';
        std::vector<std::string> dims = {"Time", "bdy_width"};
        if (s.kdim) dims.push_back(s.kdim);
        dims.push_back(x ? s.xdim : s.ydim);
        nc_.var(std::string(s.name) + "_" + side, ft, dims, {}, fieldtype);
      }
    }
    nc_.enddef();
  }

  void write(std::size_t rec, const std::string& this_date, const std::string& next_date,
             const std::map<std::string, Vars>& v) {
    nc_.put_text("Times", this_date, rec);
    nc_.put_text("md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", this_date, rec);
    nc_.put_text("md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", next_date, rec);
    const int nx = d_.ide - 1, ny = d_.jde - 1, nz = d_.kde - 1;
    for (const Spec& s : specs()) {
      const int nk = s.kdim == nullptr ? 1 : (std::string(s.kdim) == "bottom_top_stag" ? nz + 1 : nz);
      const int nxs = std::string(s.xdim) == "south_north_stag" ? ny + 1 : ny;  // length along a west/east side
      const int nys = std::string(s.ydim) == "west_east_stag" ? nx + 1 : nx;
      auto it = v.find(s.name);
      for (int t = 0; t < 2; ++t)
        for (int sd = 0; sd < 4; ++sd) {
          static const char* names[2][4] = {{"BXS", "BXE", "BYS", "BYE"}, {"BTXS", "BTXE", "BTYS", "BTYE"}};
          const int n = sd < 2 ? nxs : nys;
          std::vector<double> out(static_cast<std::size_t>(w_) * nk * n, 0.0);
          if (it != v.end()) {
            const BdyVar<T>* b = t == 0 ? it->second.val : it->second.tend;
            const F3<T>& a = sd == 0 ? b->xs : sd == 1 ? b->xe : sd == 2 ? b->ys : b->ye;
            const int k0 = a.klo();
            for (int ww = 1; ww <= w_; ++ww)
              for (int k = 0; k < nk; ++k)
                for (int p = 1; p <= n; ++p)
                  out[(std::size_t(ww - 1) * nk + k) * n + (p - 1)] = double(a(p, k0 + k, ww));
          }
          nc_.put(std::string(s.name) + "_" + names[t][sd], out, rec);
        }
    }
  }

 private:
  struct Spec {
    const char* name;
    const char* kdim;  // nullptr: 2-D
    const char* xdim;  // dimension along the west/east sides
    const char* ydim;  // dimension along the south/north sides
  };
  static const std::vector<Spec>& specs() {
    static const std::vector<Spec> s = {
        {"U", "bottom_top", "south_north", "west_east_stag"}, {"V", "bottom_top", "south_north_stag", "west_east"},
        {"W", "bottom_top_stag", "south_north", "west_east"}, {"PH", "bottom_top_stag", "south_north", "west_east"},
        {"T", "bottom_top", "south_north", "west_east"},       {"MU", nullptr, "south_north", "west_east"},
        {"QVAPOR", "bottom_top", "south_north", "west_east"},  {"HT_SHAD", nullptr, "south_north", "west_east"},
        {"PC", nullptr, "south_north", "west_east"}};
    return s;
  }
  NcOut nc_;
  Dims d_;
  int w_;
};

}  // namespace cppwrf::init_real
