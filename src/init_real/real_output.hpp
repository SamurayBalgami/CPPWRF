#pragma once
// wrfinput_d01 / wrfbdy_d01 writers (WRF names, dimensions and per-variable attributes).
#include <algorithm>
#include <cstdio>
#include <functional>
#include <map>
#include <string>

#include "init_real/met_em.hpp"
#include "init_real/nc_out.hpp"
#include "init_real/real_bdy.hpp"
#include "init_real/real_domain.hpp"
#include "init_real/real_values.hpp"
#include "init_real/wrf_schema.hpp"

namespace cppwrf::init_real {

// Global attributes of WRF's output_wrf (switch input_only / boundary_only), in WRF's order.
template <typename T>
void write_global_atts(NcOut& nc, const WrfConfig& c, const MetEm& m, const RealDomain<T>& dm,
                       const std::string& start_date, bool input) {
  auto r = [&](const std::string& n, double v) {
    if (sizeof(T) == 8) nc.att_d(n, v); else nc.att_f(n, float(v));
  };
  auto i = [&](const std::string& n, int v) { nc.att(n, v); };
  const int ide = dm.d.ide, jde = dm.d.jde, kde = dm.d.kde;
  nc.att("TITLE", std::string(" OUTPUT FROM CPPWRF REAL V4.8.0 PREPROCESSOR"));
  nc.att("START_DATE", start_date);
  if (input) nc.att("SIMULATION_START_DATE", start_date);
  i("WEST-EAST_GRID_DIMENSION", ide);
  i("SOUTH-NORTH_GRID_DIMENSION", jde);
  i("BOTTOM-TOP_GRID_DIMENSION", kde);
  r("DX", c.num("dx")); r("DY", c.num("dy"));
  i("AERCU_OPT", c.integer("aercu_opt")); r("AERCU_FCT", c.num("aercu_fct"));
  i("IDEAL_CASE", c.integer("ideal_case")); i("DIFF_6TH_SLOPEOPT", c.integer("diff_6th_slopeopt"));
  i("AUTO_LEVELS_OPT", c.integer("auto_levels_opt")); r("DIFF_6TH_THRESH", c.num("diff_6th_thresh"));
  r("DZBOT", c.num("dzbot")); r("DZSTRETCH_S", c.num("dzstretch_s")); r("DZSTRETCH_U", c.num("dzstretch_u"));
  nc.att("GRIDTYPE", std::string("C"));
  i("DIFF_OPT", c.integer("diff_opt")); i("KM_OPT", c.integer("km_opt")); i("DAMP_OPT", c.integer("damp_opt"));
  r("DAMPCOEF", c.num("dampcoef")); r("KHDIF", c.num("khdif")); r("KVDIF", c.num("kvdif"));
  i("MP_PHYSICS", c.integer("mp_physics")); i("RA_LW_PHYSICS", c.integer("ra_lw_physics"));
  i("RA_SW_PHYSICS", c.integer("ra_sw_physics")); i("SF_SFCLAY_PHYSICS", c.integer("sf_sfclay_physics"));
  i("SF_SURFACE_PHYSICS", c.integer("sf_surface_physics")); i("BL_PBL_PHYSICS", c.integer("bl_pbl_physics"));
  i("CU_PHYSICS", c.integer("cu_physics")); i("SF_LAKE_PHYSICS", c.integer("sf_lake_physics"));
  i("SURFACE_INPUT_SOURCE", c.integer("surface_input_source")); i("SST_UPDATE", c.integer("sst_update"));
  i("GHG_INPUT", c.integer("ghg_input")); i("GRID_FDDA", c.integer("grid_fdda"));
  i("GFDDA_INTERVAL_M", c.integer("auxinput10_interval_m")); i("GFDDA_END_H", c.integer("auxinput10_end_h"));
  i("GRID_SFDDA", c.integer("grid_sfdda")); i("SGFDDA_INTERVAL_M", c.integer("auxinput9_interval_m"));
  i("SGFDDA_END_H", c.integer("auxinput9_end_h")); i("HYPSOMETRIC_OPT", c.integer("hypsometric_opt"));
  i("USE_THETA_M", c.integer("use_theta_m"));
  if (input) { i("USE_MAXW_LEVEL", c.integer("use_maxw_level")); i("USE_TROP_LEVEL", c.integer("use_trop_level")); }
  i("GWD_OPT", c.integer("gwd_opt")); i("SF_URBAN_PHYSICS", c.integer("sf_urban_physics"));
  i("SF_SURFACE_MOSAIC", c.integer("sf_surface_mosaic")); i("SF_OCEAN_PHYSICS", c.integer("sf_ocean_physics"));
  i("SLUCM_DISTRIBUTED_DRAG", c.logical("slucm_distributed_drag") ? 1 : 0);
  i("DISTRIBUTED_AHE_OPT", c.integer("distributed_ahe_opt"));
  if (input) nc.att("SIMULATION_INITIALIZATION_TYPE", std::string("REAL-DATA CASE"));
  i("WEST-EAST_PATCH_START_UNSTAG", 1); i("WEST-EAST_PATCH_END_UNSTAG", ide - 1);
  i("WEST-EAST_PATCH_START_STAG", 1); i("WEST-EAST_PATCH_END_STAG", ide);
  i("SOUTH-NORTH_PATCH_START_UNSTAG", 1); i("SOUTH-NORTH_PATCH_END_UNSTAG", jde - 1);
  i("SOUTH-NORTH_PATCH_START_STAG", 1); i("SOUTH-NORTH_PATCH_END_STAG", jde);
  i("BOTTOM-TOP_PATCH_START_UNSTAG", 1); i("BOTTOM-TOP_PATCH_END_UNSTAG", kde - 1);
  i("BOTTOM-TOP_PATCH_START_STAG", 1); i("BOTTOM-TOP_PATCH_END_STAG", kde);
  i("GRID_ID", dm.o.dom); i("PARENT_ID", dm.o.parent_id); i("I_PARENT_START", dm.o.i_parent_start);
  i("J_PARENT_START", dm.o.j_parent_start); i("PARENT_GRID_RATIO", dm.o.parent_grid_ratio);
  {  // grid%dt: time_step divided by the parent time-step ratios down to this domain
    double dt = c.num("time_step", 1) + c.num("time_step_fract_num", 1) / std::max(1.0, c.num("time_step_fract_den", 1));
    for (int id = dm.o.dom; id > 1; id = c.integer("parent_id", id)) dt = dt / c.num("parent_time_step_ratio", id);
    r("DT", dt);
  }
  // the projection comes from the metgrid input (float attributes there)
  for (const char* a : {"CEN_LAT", "CEN_LON", "TRUELAT1", "TRUELAT2", "MOAD_CEN_LAT", "STAND_LON", "POLE_LAT", "POLE_LON"})
    r(a, double(float(m.att_real(a))));
  if (input) {
    int jy, jd;
    double gm;
    get_julgmt(start_date, jy, jd, gm);
    r("GMT", gm); i("JULYR", jy); i("JULDAY", jd);
  }
  const int mp = m.att_int("MAP_PROJ");
  i("MAP_PROJ", mp);
  nc.att("MAP_PROJ_CHAR", std::string(mp == 1 ? "Lambert Conformal" : mp == 2 ? "Polar Stereographic"
                                     : mp == 3 ? "Mercator" : mp == 6 ? "Cylindrical Equidistant" : "unknown"));
  nc.att("MMINLU", dm.mminlu);
  i("NUM_LAND_CAT", m.att_int("NUM_LAND_CAT")); i("ISWATER", dm.iswater); i("ISLAKE", dm.islake);
  i("ISICE", dm.isice); i("ISURBAN", dm.isurban); i("ISOILWATER", dm.isoilwater);
  i("HYBRID_OPT", c.integer("hybrid_opt")); r("ETAC", c.integer("hybrid_opt") == 0 ? 0.0 : c.num("etac"));
}

// wrfinput_d01: every variable of WRF's input stream for this configuration (tools/wrf_schema), with the
// fields real.exe computes and WRF's initial value 0 elsewhere.
template <typename T>
void write_wrfinput(const RealDomain<T>& dm, const MetEm& m, const WrfSchema& sc, const WrfConfig& c,
                    const std::string& path, const std::string& date) {
  const int ide = dm.d.ide, jde = dm.d.jde, kde = dm.d.kde;
  std::vector<OutVar> vars = sc.input_vars(c, ide, jde, kde, int(sizeof(T)));
  NcOut nc(path, !c.logical("use_netcdf_classic"));
  nc.dim("Time", 0);
  nc.dim("DateStrLen", 19);
  {
    // dimensions in the order WRF defines them (Fortran order of each variable); an unnamed dimension
    // reuses an existing one of the same length or becomes DIM<nnnn> (external/io_netcdf)
    std::vector<std::pair<std::string, std::size_t>> defined = {{"Time", 0}, {"DateStrLen", 19}};
    for (OutVar& v : vars)
      for (int k = int(v.dimnames.size()) - 1; k >= 0; --k) {
        std::string& n = v.dimnames[std::size_t(k)];
        const std::size_t len = v.dimlens[std::size_t(k)];
        if (n.empty()) {
          for (auto& [dn, dl] : defined)
            if (dl == len && dn != "Time") { n = dn; break; }
          if (n.empty()) {
            char b[16];
            std::snprintf(b, sizeof b, "DIM%04d", int(defined.size()));
            n = b;
          }
        }
        bool found = false;
        for (auto& [dn, dl] : defined) found = found || dn == n;
        if (!found) {
          nc.dim(n, len);
          defined.push_back({n, len});
        }
      }
  }
  nc.var("Times", 'c', {"Time", "DateStrLen"}, {}, 0);
  const int ftr = sizeof(T) == 8 ? 105 : 104;
  for (const OutVar& v : vars) {
    std::vector<std::string> dims = {"Time"};
    for (const std::string& n : v.dimnames) dims.push_back(n);
    std::map<std::string, std::string> at = {{"MemoryOrder", v.memorder_att}, {"description", v.desc},
                                             {"units", v.units}, {"stagger", v.stagger}};
    if (!v.coordinates.empty()) at["coordinates"] = v.coordinates;
    const char t = v.type == 'l' ? 'i' : v.type;
    const int ft = t == 'f' ? 104 : t == 'd' ? 105 : t == 'i' ? 106 : ftr;
    nc.var(v.name, t, dims, at, ft);
  }
  write_global_atts(nc, c, m, dm, date, true);
  nc.enddef();
  nc.put_text("Times", date);
  const auto values = wrfinput_values(dm, m);
  for (const OutVar& v : vars) {
    std::size_t n = 1;
    for (std::size_t l : v.dimlens) n *= l;
    std::vector<double> data(n, 0.0);
    auto it = values.find(v.name);
    if (it != values.end()) {
      data = it->second(v);
      if (data.size() != n) throw std::runtime_error("wrfinput: size mismatch for " + v.name);
    }
    if (v.type == 'i' || v.type == 'l') nc.put_int(v.name, std::vector<int>(data.begin(), data.end()));
    else nc.put(v.name, data);
  }
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
  WrfBdyWriter(const std::string& path, const RealDomain<T>& dm, const MetEm& m, const WrfSchema& sc,
               const WrfConfig& c, const std::string& start_date, int w)
      : nc_(path, !c.logical("use_netcdf_classic")), d_(dm.d), w_(w) {
    const int nx = d_.ide - 1, ny = d_.jde - 1, nz = d_.kde - 1;
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
    // share/wrf_bdyout: U V W PH T MU, the moist and scalar members with boundary data, HT_SHAD, PC
    auto state = [&](const char* name, const char* varname, const char* kdim, const char* xdim, const char* ydim) {
      const SchemaVar* v = sc.find(varname);
      if (!v) throw std::runtime_error(std::string("wrfbdy: no schema entry for ") + varname);
      specs_.push_back({name, kdim, xdim, ydim, v->desc_full, v->units, v->stagger, ""});
    };
    state("U", "u_2", "bottom_top", "south_north", "west_east_stag");
    state("V", "v_2", "bottom_top", "south_north_stag", "west_east");
    state("W", "w_2", "bottom_top_stag", "south_north", "west_east");
    state("PH", "ph_2", "bottom_top_stag", "south_north", "west_east");
    state("T", "t_2", "bottom_top", "south_north", "west_east");
    state("MU", "mu_2", nullptr, "south_north", "west_east");
    for (const char* array : {"moist", "scalar"})
      for (const SchemaMember* mb : sc.boundary_members(array, c))
        specs_.push_back({mb->dname, "bottom_top", "south_north", "west_east", mb->desc, mb->units, "", array});
    state("HT_SHAD", "ht_shad", nullptr, "south_north", "west_east");
    state("PC", "pc_2", nullptr, "south_north", "west_east");
    const char ft = sizeof(T) == 8 ? 'd' : 'f';
    const int fieldtype = sizeof(T) == 8 ? 105 : 104;
    // write order: a state variable's 8 arrays together; for a 4-D array all members' values, then
    // all members' tendencies
    for (std::size_t k = 0; k < specs_.size();) {
      if (specs_[k].state()) {
        for (int sd = 0; sd < 8; ++sd) order_.push_back({k, sd});
        ++k;
      } else {
        std::size_t e = k;
        while (e < specs_.size() && specs_[e].array == specs_[k].array) ++e;
        for (int t = 0; t < 2; ++t)
          for (std::size_t q = k; q < e; ++q)
            for (int sd = 0; sd < 4; ++sd) order_.push_back({q, t * 4 + sd});
        k = e;
      }
    }
    static const char* sides[8] = {"BXS", "BXE", "BYS", "BYE", "BTXS", "BTXE", "BTYS", "BTYE"};
    for (auto [k, sdi] : order_) {
      const Spec& s = specs_[k];
      {
        const char* side = sides[sdi];
        const bool tend = side[1] == 'T';
        const char* sd = side + (tend ? 2 : 1);  // XS, XE, YS, YE
        std::vector<std::string> dims = {"Time", "bdy_width"};
        if (s.kdim) dims.push_back(s.kdim);
        dims.push_back(sd[0] == 'X' ? s.xdim : s.ydim);
        std::string mo = std::string(sd) + (s.kdim ? "Z" : " ");
        std::string desc = s.state() ? (tend ? "bdy tend " : "bdy ") + s.desc : s.desc;
        std::string units = s.state() && tend ? "(" + s.units + ")/dt" : s.units;
        nc_.var(s.name + "_" + side, ft, dims,
                {{"MemoryOrder", mo}, {"description", desc}, {"units", units}, {"stagger", s.stagger}}, fieldtype);
      }
    }
    write_global_atts(nc_, c, m, dm, start_date, false);
    nc_.enddef();
  }

  void write(std::size_t rec, const std::string& this_date, const std::string& next_date,
             const std::map<std::string, Vars>& v) {
    nc_.put_text("Times", this_date, rec);
    nc_.put_text("md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", this_date, rec);
    nc_.put_text("md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_", next_date, rec);
    const int nx = d_.ide - 1, ny = d_.jde - 1, nz = d_.kde - 1;
    for (auto [k, sdi] : order_) {
      const Spec& s = specs_[k];
      const int nk = s.kdim == nullptr ? 1 : (std::string(s.kdim) == "bottom_top_stag" ? nz + 1 : nz);
      const int nxs = std::string(s.xdim) == "south_north_stag" ? ny + 1 : ny;  // length along a west/east side
      const int nys = std::string(s.ydim) == "west_east_stag" ? nx + 1 : nx;
      auto it = v.find(s.name);
      {
        const int t = sdi / 4, sd = sdi % 4;
        {
          static const char* names[2][4] = {{"BXS", "BXE", "BYS", "BYE"}, {"BTXS", "BTXE", "BTYS", "BTYE"}};
          const int n = sd < 2 ? nxs : nys;
          std::vector<double> out(static_cast<std::size_t>(w_) * nk * n, 0.0);
          if (it != v.end()) {
            const BdyVar<T>* b = t == 0 ? it->second.val : it->second.tend;
            const F3<T>& a = sd == 0 ? b->xs : sd == 1 ? b->xe : sd == 2 ? b->ys : b->ye;
            const int k0 = a.klo();
            for (int ww = 1; ww <= w_; ++ww)
              for (int kk = 0; kk < nk; ++kk)
                for (int p = 1; p <= n; ++p)
                  out[(std::size_t(ww - 1) * nk + kk) * n + (p - 1)] = double(a(p, k0 + kk, ww));
          }
          nc_.put(s.name + "_" + names[t][sd], out, rec);
        }
      }
    }
  }

 private:
  struct Spec {
    std::string name;
    const char* kdim;  // nullptr: 2-D
    const char* xdim;  // dimension along the west/east sides
    const char* ydim;  // dimension along the south/north sides
    std::string desc, units, stagger;
    std::string array;  // 4-D array of a member; empty for a state variable (bdy/bdy tend descriptions)
    bool state() const { return array.empty(); }
  };
  std::vector<Spec> specs_;
  std::vector<std::pair<std::size_t, int>> order_;  // (spec, side 0..7)
  NcOut nc_;
  Dims d_;
  int w_;
};

}  // namespace cppwrf::init_real
