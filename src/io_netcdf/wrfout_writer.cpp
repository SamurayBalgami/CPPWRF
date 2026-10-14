#include "io_netcdf/wrfout_writer.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <map>
#include <set>
#include <stdexcept>

namespace cppwrf::io_netcdf {

namespace {

std::string to_upper(std::string_view s) {
  std::string r(s);
  for (char& c : r) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return r;
}

const VarView* find_var(const std::vector<VarView>& vars, std::string_view name) {
  // 1. Exact case-sensitive match
  for (const auto& v : vars) {
    if (v.name == name) return &v;
  }
  std::string target = to_upper(name);
  // 2. Specific perturbation / temperature aliases (strict, no accidental full-field matching)
  if (target == "P") {
    for (const auto& v : vars) {
      if (to_upper(v.name) == "P_PERT") return &v;
    }
    return nullptr;
  }
  if (target == "PH") {
    for (const auto& v : vars) {
      if (to_upper(v.name) == "PH_PERT") return &v;
    }
    return nullptr;
  }
  if (target == "MU") {
    for (const auto& v : vars) {
      if (to_upper(v.name) == "MU_PERT") return &v;
    }
    return nullptr;
  }
  if (target == "T") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "T_PERT" || n == "TH_300" || n == "THETA_300") return &v;
    }
    return nullptr;
  }
  // 3. Exact case-insensitive match for other variables
  for (const auto& v : vars) {
    if (to_upper(v.name) == target) return &v;
  }
  if (target == "QVAPOR") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "QV") return &v;
    }
  }
  if (target == "QCLOUD") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "QC") return &v;
    }
  }
  if (target == "QRAIN") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "QR") return &v;
    }
  }
  if (target == "QICE") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "QI") return &v;
    }
  }
  if (target == "QSNOW") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "QS") return &v;
    }
  }
  if (target == "QGRAUP") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "QG") return &v;
    }
  }
  return nullptr;
}

// Packs a variable into a contiguous float buffer of shape [nz, ny, nx]
void pack_to_float(const VarView& v, int expected_nx, int expected_ny, int expected_nz,
                   std::vector<float>& out) {
  out.resize(size_t(expected_nx) * expected_ny * expected_nz);
  const int nx = v.nx > 0 ? v.nx : 1;
  const int ny = v.ny > 0 ? v.ny : 1;
  const int nz = v.nz > 0 ? v.nz : 1;
  const std::ptrdiff_t sj = v.sj != 0 ? v.sj : nx;
  const std::ptrdiff_t sk = v.sk != 0 ? v.sk : (std::ptrdiff_t)nx * ny;

  for (int k = 0; k < expected_nz; ++k) {
    int src_k = std::min(k, nz - 1);
    for (int j = 0; j < expected_ny; ++j) {
      int src_j = std::min(j, ny - 1);
      for (int i = 0; i < expected_nx; ++i) {
        int src_i = std::min(i, nx - 1);
        std::ptrdiff_t src_idx = src_i + src_j * sj + src_k * sk;
        size_t dst_idx = (size_t(k) * expected_ny + j) * expected_nx + i;
        if (v.dtype == DType::F64) {
          const double* p = static_cast<const double*>(v.data);
          out[dst_idx] = static_cast<float>(p[src_idx]);
        } else {
          const float* p = static_cast<const float*>(v.data);
          out[dst_idx] = p[src_idx];
        }
      }
    }
  }
}

void put_att_text(int ncid, int varid, const char* name, const std::string& val) {
  NC_CHECK(nc_put_att_text(ncid, varid, name, val.size(), val.data()));
}

void put_att_int(int ncid, int varid, const char* name, int val) {
  NC_CHECK(nc_put_att_int(ncid, varid, name, NC_INT, 1, &val));
}

void put_att_float(int ncid, int varid, const char* name, float val) {
  NC_CHECK(nc_put_att_float(ncid, varid, name, NC_FLOAT, 1, &val));
}

struct VarDef {
  std::string name;
  std::vector<int> dimids;
  std::string memory_order;
  std::string stagger;
  std::string desc;
  std::string units;
  std::string coords;
  int varid = -1;
};

}  // namespace

std::string format_wrf_date(const std::string& base_date, double time_sec) {
  int Y = 2000, M = 1, D = 1, hr = 0, mn = 0, sc = 0;
  if (base_date.size() >= 19) {
    std::sscanf(base_date.c_str(), "%4d-%2d-%2d_%2d:%2d:%2d", &Y, &M, &D, &hr, &mn, &sc);
  }
  // Convert base date to civil days from 1970-01-01 (proleptic Gregorian)
  const int y = Y - (M <= 2);
  const int era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(y - era * 400);
  const unsigned doy = (153 * (M + (M > 2 ? -3 : 9)) + 2) / 5 + D - 1;
  const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  const long base_days = long(era) * 146097 + doe - 719468;

  long total_sec = base_days * 86400L + hr * 3600L + mn * 60L + sc + std::lround(time_sec);
  long new_days = total_sec / 86400L;
  long sec_day = total_sec % 86400L;
  if (sec_day < 0) {
    sec_day += 86400L;
    new_days -= 1;
  }
  hr = int(sec_day / 3600);
  mn = int((sec_day % 3600) / 60);
  sc = int(sec_day % 60);

  // Convert civil days back to date (Hinnant's algorithm)
  const long z = new_days + 719468;
  const long res_era = (z >= 0 ? z : z - 146096) / 146097;
  const unsigned res_doe = static_cast<unsigned>(z - res_era * 146097);
  const unsigned res_yoe = (res_doe - res_doe / 1460 + res_doe / 36524 - res_doe / 146096) / 365;
  const long res_y = static_cast<long>(res_yoe) + res_era * 400;
  const unsigned res_doy = res_doe - (365 * res_yoe + res_yoe / 4 - res_yoe / 100);
  const unsigned mp = (5 * res_doy + 2) / 153;
  const unsigned res_d = res_doy - (153 * mp + 2) / 5 + 1;
  const unsigned res_m = mp < 10 ? mp + 3 : mp - 9;
  const long final_y = res_y + (res_m <= 2);

  char buf[32];
  std::snprintf(buf, sizeof(buf), "%04ld-%02u-%02u_%02d:%02d:%02d", final_y, res_m, res_d, hr, mn, sc);
  return std::string(buf);
}

class WrfoutWriter::Impl {
 public:
  Impl(std::string out_dir, WrfoutMeta meta) : out_dir_(std::move(out_dir)), meta_(std::move(meta)) {}

  void write_frame(long step, double time_sec, const std::vector<VarView>& vars,
                   const std::string& target_file, std::string& out_written_file) {
    (void)step;
    std::string date_str = format_wrf_date(meta_.start_date, time_sec);
    std::string file_path = target_file;
    if (file_path.empty()) {
      std::filesystem::create_directories(out_dir_);
      file_path = out_dir_ + "/wrfout_d01_" + meta_.start_date;
    } else {
      auto p = std::filesystem::path(file_path).parent_path();
      if (!p.empty()) std::filesystem::create_directories(p);
    }
    out_written_file = file_path;

    if (!is_initialized_ || file_path != current_file_) {
      ncid_.close();
      current_file_ = file_path;
      defs_.clear();
      generic_.clear();
      time_idx_ = 0;

      bool file_exists = std::filesystem::exists(file_path);
      // W5: recreate with NC_CLOBBER at run start, append only on restart
      bool do_open_existing = meta_.is_restart && file_exists;

      int dim_time = -1, dim_datestrlen = -1, dim_we = -1, dim_sn = -1, dim_bt = -1;
      int dim_we_stag = -1, dim_sn_stag = -1, dim_bt_stag = -1;

      if (!do_open_existing) {
        int nid = -1;
        NC_CHECK(nc_create(file_path.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &nid));
        ncid_ = ScopedNcFile(nid);

        // 1. Dimensions
        NC_CHECK(nc_def_dim(ncid_.get(), "Time", NC_UNLIMITED, &dim_time));
        NC_CHECK(nc_def_dim(ncid_.get(), "DateStrLen", 19, &dim_datestrlen));
        NC_CHECK(nc_def_dim(ncid_.get(), "west_east", meta_.nx, &dim_we));
        NC_CHECK(nc_def_dim(ncid_.get(), "south_north", meta_.ny, &dim_sn));
        NC_CHECK(nc_def_dim(ncid_.get(), "bottom_top", meta_.nz, &dim_bt));
        NC_CHECK(nc_def_dim(ncid_.get(), "bottom_top_stag", meta_.nz + 1, &dim_bt_stag));
        NC_CHECK(nc_def_dim(ncid_.get(), "west_east_stag", meta_.nx + 1, &dim_we_stag));
        NC_CHECK(nc_def_dim(ncid_.get(), "south_north_stag", meta_.ny + 1, &dim_sn_stag));

        // 2. Global attributes
        write_global_attributes(ncid_.get());

        // 3. Define standard variables
        // Times
        defs_["Times"] = {"Times", {dim_time, dim_datestrlen}, "", "", "WRF Output Times", "", ""};

        // 3D Staggered U, V, W, PH, PHB
        defs_["U"] = {"U", {dim_time, dim_bt, dim_sn, dim_we_stag}, "XYZ", "X", "x-wind component", "m s-1", "XLONG_U XLAT_U XTIME"};
        defs_["V"] = {"V", {dim_time, dim_bt, dim_sn_stag, dim_we}, "XYZ", "Y", "y-wind component", "m s-1", "XLONG_V XLAT_V XTIME"};
        defs_["W"] = {"W", {dim_time, dim_bt_stag, dim_sn, dim_we}, "XYZ", "Z", "z-wind component", "m s-1", "XLONG XLAT XTIME"};
        defs_["PH"] = {"PH", {dim_time, dim_bt_stag, dim_sn, dim_we}, "XYZ", "Z", "perturbation geopotential", "m2 s-2", "XLONG XLAT XTIME"};
        defs_["PHB"] = {"PHB", {dim_time, dim_bt_stag, dim_sn, dim_we}, "XYZ", "Z", "base-state geopotential", "m2 s-2", "XLONG XLAT XTIME"};

        // 3D Unstaggered T, P, PB
        defs_["T"] = {"T", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "perturbation potential temperature (theta-t0)", "K", "XLONG XLAT XTIME"};
        defs_["P"] = {"P", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "perturbation pressure", "Pa", "XLONG XLAT XTIME"};
        defs_["PB"] = {"PB", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "BASE STATE PRESSURE", "Pa", "XLONG XLAT XTIME"};

        // Optional moist species and surface diagnostics
        if (find_var(vars, "QVAPOR") || find_var(vars, "QV")) {
          defs_["QVAPOR"] = {"QVAPOR", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "Water vapor mixing ratio", "kg kg-1", "XLONG XLAT XTIME"};
        }
        if (find_var(vars, "QCLOUD") || find_var(vars, "QC")) {
          defs_["QCLOUD"] = {"QCLOUD", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "Cloud water mixing ratio", "kg kg-1", "XLONG XLAT XTIME"};
        }
        if (find_var(vars, "QRAIN") || find_var(vars, "QR")) {
          defs_["QRAIN"] = {"QRAIN", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "Rain water mixing ratio", "kg kg-1", "XLONG XLAT XTIME"};
        }
        if (find_var(vars, "QICE") || find_var(vars, "QI")) {
          defs_["QICE"] = {"QICE", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "Ice mixing ratio", "kg kg-1", "XLONG XLAT XTIME"};
        }
        if (find_var(vars, "QSNOW") || find_var(vars, "QS")) {
          defs_["QSNOW"] = {"QSNOW", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "Snow mixing ratio", "kg kg-1", "XLONG XLAT XTIME"};
        }
        if (find_var(vars, "QGRAUP") || find_var(vars, "QG")) {
          defs_["QGRAUP"] = {"QGRAUP", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "Graupel mixing ratio", "kg kg-1", "XLONG XLAT XTIME"};
        }
        if (find_var(vars, "RAINNC")) {
          defs_["RAINNC"] = {"RAINNC", {dim_time, dim_sn, dim_we}, "XY ", "", "ACCUMULATED TOTAL GRID SCALE PRECIPITATION", "mm", "XLONG XLAT XTIME"};
        }
        if (find_var(vars, "SNOWNC")) {
          defs_["SNOWNC"] = {"SNOWNC", {dim_time, dim_sn, dim_we}, "XY ", "", "ACCUMULATED TOTAL GRID SCALE SNOW", "mm", "XLONG XLAT XTIME"};
        }
        if (find_var(vars, "GRAUPELNC")) {
          defs_["GRAUPELNC"] = {"GRAUPELNC", {dim_time, dim_sn, dim_we}, "XY ", "", "ACCUMULATED TOTAL GRID SCALE GRAUPEL", "mm", "XLONG XLAT XTIME"};
        }
        if (find_var(vars, "SR")) {
          defs_["SR"] = {"SR", {dim_time, dim_sn, dim_we}, "XY ", "", "fraction of frozen precipitation", "-", "XLONG XLAT XTIME"};
        }

        // 2D MU, MUB, XLAT, XLONG, HGT
        defs_["MU"] = {"MU", {dim_time, dim_sn, dim_we}, "XY ", "", "perturbation dry air mass in column", "Pa", "XLONG XLAT XTIME"};
        defs_["MUB"] = {"MUB", {dim_time, dim_sn, dim_we}, "XY ", "", "base state dry air mass in column", "Pa", "XLONG XLAT XTIME"};
        defs_["XLAT"] = {"XLAT", {dim_time, dim_sn, dim_we}, "XY ", "", "LATITUDE, SOUTH IS NEGATIVE", "degree_north", ""};
        defs_["XLONG"] = {"XLONG", {dim_time, dim_sn, dim_we}, "XY ", "", "LONGITUDE, WEST IS NEGATIVE", "degree_east", ""};
        defs_["HGT"] = {"HGT", {dim_time, dim_sn, dim_we}, "XY ", "", "Terrain Height", "m", "XLONG XLAT XTIME"};

        // 1D Vertical coordinates
        defs_["ZNU"] = {"ZNU", {dim_time, dim_bt}, "Z  ", "", "eta values on half (mass) levels", "", ""};
        defs_["ZNW"] = {"ZNW", {dim_time, dim_bt_stag}, "Z  ", "Z", "eta values on full (w) levels", "", ""};

        // Hybrid vertical coordinate terms
        defs_["C1H"] = {"C1H", {dim_time, dim_bt}, "Z  ", "", "half levels, c1h = d bf / d eta, using znw", "Dimensionless", ""};
        defs_["C2H"] = {"C2H", {dim_time, dim_bt}, "Z  ", "", "half levels, c2h = (1-c1h)*(p0-pt)", "Pa", ""};
        defs_["C1F"] = {"C1F", {dim_time, dim_bt_stag}, "Z  ", "Z", "full levels, c1f = d bf / d eta, using znu", "Dimensionless", ""};
        defs_["C2F"] = {"C2F", {dim_time, dim_bt_stag}, "Z  ", "Z", "full levels, c2f = (1-c1f)*(p0-pt)", "Pa", ""};
        defs_["C3H"] = {"C3H", {dim_time, dim_bt}, "Z  ", "", "half levels, c3h = bh", "Dimensionless", ""};
        defs_["C4H"] = {"C4H", {dim_time, dim_bt}, "Z  ", "", "half levels, c4h = (eta-bh)*(p0-pt), using znu", "Pa", ""};
        defs_["C3F"] = {"C3F", {dim_time, dim_bt_stag}, "Z  ", "Z", "full levels, c3f = bf", "Dimensionless", ""};
        defs_["C4F"] = {"C4F", {dim_time, dim_bt_stag}, "Z  ", "Z", "full levels, c4f = (eta-bf)*(p0-pt), using znw", "Pa", ""};

        // 0D Scalars
        defs_["RDX"] = {"RDX", {dim_time}, "0  ", "", "INVERSE X GRID LENGTH", "m-1", ""};
        defs_["RDY"] = {"RDY", {dim_time}, "0  ", "", "INVERSE Y GRID LENGTH", "m-1", ""};
        defs_["P_TOP"] = {"P_TOP", {dim_time}, "0  ", "", "PRESSURE TOP OF THE MODEL", "Pa", ""};
        defs_["XTIME"] = {"XTIME", {dim_time}, "0  ", "", "minutes since simulation start", "minutes", ""};

        // Any other mass-point field of the frame (physics diagnostics, R*BLTEN, soil state):
        // 2D -> (Time, south_north, west_east), nz levels -> bottom_top, other level counts ->
        // soil_layers_stag (WRF's soil dimension).
        int dim_soil = -1;
        for (const VarView& v : vars) {
          if (v.stagger != 'M' || defs_.count(v.name) || v.nx != meta_.nx || v.ny != meta_.ny) continue;
          if (v.nz == 1) {
            defs_[v.name] = {v.name, {dim_time, dim_sn, dim_we}, "XY ", "", v.name, "", "XLONG XLAT XTIME"};
          } else if (v.nz == meta_.nz) {
            defs_[v.name] = {v.name, {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", v.name, "", "XLONG XLAT XTIME"};
          } else {
            if (dim_soil < 0) NC_CHECK(nc_def_dim(ncid_.get(), "soil_layers_stag", v.nz, &dim_soil));
            defs_[v.name] = {v.name, {dim_time, dim_soil, dim_sn, dim_we}, "XYZ", "Z", v.name, "", "XLONG XLAT XTIME"};
          }
          generic_.insert(v.name);
        }

        for (auto& [name, d] : defs_) {
          nc_type xtype = (name == "Times") ? NC_CHAR : NC_FLOAT;
          NC_CHECK(nc_def_var(ncid_.get(), name.c_str(), xtype, int(d.dimids.size()), d.dimids.data(), &d.varid));
          if (name != "Times") {
            put_att_int(ncid_.get(), d.varid, "FieldType", 104);
            put_att_text(ncid_.get(), d.varid, "MemoryOrder", d.memory_order);
            put_att_text(ncid_.get(), d.varid, "description", d.desc);
            put_att_text(ncid_.get(), d.varid, "units", d.units);
            put_att_text(ncid_.get(), d.varid, "stagger", d.stagger);
            if (!d.coords.empty()) {
              put_att_text(ncid_.get(), d.varid, "coordinates", d.coords);
            }
          }
        }

        NC_CHECK(nc_enddef(ncid_.get()));
        time_idx_ = 0;
      } else {
        int nid = -1;
        NC_CHECK(nc_open(file_path.c_str(), NC_WRITE, &nid));
        ncid_ = ScopedNcFile(nid);

        NC_CHECK(nc_inq_dimid(ncid_.get(), "Time", &dim_time));
        size_t tlen = 0;
        NC_CHECK(nc_inq_dimlen(ncid_.get(), dim_time, &tlen));
        time_idx_ = tlen;

        // Populate varids from existing file
        int nvars = 0;
        NC_CHECK(nc_inq_nvars(ncid_.get(), &nvars));
        for (int i = 0; i < nvars; ++i) {
          char vname[NC_MAX_NAME + 1];
          NC_CHECK(nc_inq_varname(ncid_.get(), i, vname));
          defs_[vname] = VarDef{vname, {}, "", "", "", "", "", i};
        }
        for (const VarView& v : vars)  // restart append: the extra fields already in the file
          if (v.stagger == 'M' && defs_.count(v.name) && v.nx == meta_.nx && v.ny == meta_.ny) generic_.insert(v.name);
      }
      is_initialized_ = true;
    }

    // Write data at time_idx_
    write_frame_data(ncid_.get(), time_idx_, time_sec, date_str, vars, defs_);
    NC_CHECK(nc_sync(ncid_.get()));
    time_idx_++;
  }

 private:
  void write_global_attributes(int ncid) {
    put_att_text(ncid, NC_GLOBAL, "TITLE", "OUTPUT FROM CPPWRF");
    put_att_text(ncid, NC_GLOBAL, "START_DATE", meta_.start_date);
    put_att_text(ncid, NC_GLOBAL, "SIMULATION_START_DATE", meta_.simulation_start_date);

    put_att_int(ncid, NC_GLOBAL, "WEST-EAST_GRID_DIMENSION", meta_.nx + 1);
    put_att_int(ncid, NC_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION", meta_.ny + 1);
    put_att_int(ncid, NC_GLOBAL, "BOTTOM-TOP_GRID_DIMENSION", meta_.nz + 1);

    put_att_float(ncid, NC_GLOBAL, "DX", float(meta_.dx));
    put_att_float(ncid, NC_GLOBAL, "DY", float(meta_.dy));
    put_att_float(ncid, NC_GLOBAL, "DT", float(meta_.dt));

    put_att_text(ncid, NC_GLOBAL, "GRIDTYPE", "C");
    put_att_int(ncid, NC_GLOBAL, "DIFF_OPT", meta_.diff_opt);
    put_att_int(ncid, NC_GLOBAL, "KM_OPT", meta_.km_opt);
    put_att_int(ncid, NC_GLOBAL, "DAMP_OPT", meta_.damp_opt);
    put_att_float(ncid, NC_GLOBAL, "DAMPCOEF", float(meta_.dampcoef));
    put_att_float(ncid, NC_GLOBAL, "KHDIF", float(meta_.khdif));
    put_att_float(ncid, NC_GLOBAL, "KVDIF", float(meta_.kvdif));

    put_att_int(ncid, NC_GLOBAL, "MP_PHYSICS", meta_.mp_physics);
    put_att_int(ncid, NC_GLOBAL, "RA_LW_PHYSICS", meta_.ra_lw_physics);
    put_att_int(ncid, NC_GLOBAL, "RA_SW_PHYSICS", meta_.ra_sw_physics);
    put_att_int(ncid, NC_GLOBAL, "SF_SFCLAY_PHYSICS", meta_.sf_sfclay_physics);
    put_att_int(ncid, NC_GLOBAL, "SF_SURFACE_PHYSICS", meta_.sf_surface_physics);
    put_att_int(ncid, NC_GLOBAL, "BL_PBL_PHYSICS", meta_.bl_pbl_physics);
    put_att_int(ncid, NC_GLOBAL, "CU_PHYSICS", meta_.cu_physics);
    put_att_int(ncid, NC_GLOBAL, "SF_LAKE_PHYSICS", 0);

    put_att_int(ncid, NC_GLOBAL, "SURFACE_INPUT_SOURCE", 1);
    put_att_int(ncid, NC_GLOBAL, "SST_UPDATE", 0);
    put_att_int(ncid, NC_GLOBAL, "GRID_FDDA", 0);
    put_att_int(ncid, NC_GLOBAL, "GFDDA_INTERVAL_M", 0);
    put_att_int(ncid, NC_GLOBAL, "GFDDA_END_H", 0);
    put_att_int(ncid, NC_GLOBAL, "GRID_SFDDA", 0);
    put_att_int(ncid, NC_GLOBAL, "SGFDDA_INTERVAL_M", 0);
    put_att_int(ncid, NC_GLOBAL, "SGFDDA_END_H", 0);
    put_att_int(ncid, NC_GLOBAL, "HYPSOMETRIC_OPT", meta_.hypsometric_opt);
    put_att_int(ncid, NC_GLOBAL, "USE_THETA_M", meta_.use_theta_m);
    put_att_int(ncid, NC_GLOBAL, "GWD_OPT", 0);
    put_att_int(ncid, NC_GLOBAL, "SF_URBAN_PHYSICS", 0);
    put_att_int(ncid, NC_GLOBAL, "SF_SURFACE_MOSAIC", 0);
    put_att_int(ncid, NC_GLOBAL, "SF_OCEAN_PHYSICS", 0);

    put_att_int(ncid, NC_GLOBAL, "WEST-EAST_PATCH_START_UNSTAG", 1);
    put_att_int(ncid, NC_GLOBAL, "WEST-EAST_PATCH_END_UNSTAG", meta_.nx);
    put_att_int(ncid, NC_GLOBAL, "WEST-EAST_PATCH_START_STAG", 1);
    put_att_int(ncid, NC_GLOBAL, "WEST-EAST_PATCH_END_STAG", meta_.nx + 1);

    put_att_int(ncid, NC_GLOBAL, "SOUTH-NORTH_PATCH_START_UNSTAG", 1);
    put_att_int(ncid, NC_GLOBAL, "SOUTH-NORTH_PATCH_END_UNSTAG", meta_.ny);
    put_att_int(ncid, NC_GLOBAL, "SOUTH-NORTH_PATCH_START_STAG", 1);
    put_att_int(ncid, NC_GLOBAL, "SOUTH-NORTH_PATCH_END_STAG", meta_.ny + 1);

    put_att_int(ncid, NC_GLOBAL, "BOTTOM-TOP_PATCH_START_UNSTAG", 1);
    put_att_int(ncid, NC_GLOBAL, "BOTTOM-TOP_PATCH_END_UNSTAG", meta_.nz);
    put_att_int(ncid, NC_GLOBAL, "BOTTOM-TOP_PATCH_START_STAG", 1);
    put_att_int(ncid, NC_GLOBAL, "BOTTOM-TOP_PATCH_END_STAG", meta_.nz + 1);

    put_att_int(ncid, NC_GLOBAL, "GRID_ID", 1);
    put_att_int(ncid, NC_GLOBAL, "PARENT_ID", 0);
    put_att_int(ncid, NC_GLOBAL, "I_PARENT_START", 1);
    put_att_int(ncid, NC_GLOBAL, "J_PARENT_START", 1);
    put_att_int(ncid, NC_GLOBAL, "PARENT_GRID_RATIO", 1);

    put_att_float(ncid, NC_GLOBAL, "CEN_LAT", float(meta_.cen_lat));
    put_att_float(ncid, NC_GLOBAL, "CEN_LON", float(meta_.cen_lon));
    put_att_float(ncid, NC_GLOBAL, "TRUELAT1", float(meta_.truelat1));
    put_att_float(ncid, NC_GLOBAL, "TRUELAT2", float(meta_.truelat2));
    put_att_float(ncid, NC_GLOBAL, "MOAD_CEN_LAT", float(meta_.moad_cen_lat));
    put_att_float(ncid, NC_GLOBAL, "STAND_LON", float(meta_.stand_lon));
    put_att_float(ncid, NC_GLOBAL, "POLE_LAT", float(meta_.pole_lat));
    put_att_float(ncid, NC_GLOBAL, "POLE_LON", float(meta_.pole_lon));

    put_att_float(ncid, NC_GLOBAL, "GMT", 0.0f);
    put_att_int(ncid, NC_GLOBAL, "JULYR", 1);
    put_att_int(ncid, NC_GLOBAL, "JULDAY", 1);

    put_att_int(ncid, NC_GLOBAL, "MAP_PROJ", meta_.map_proj);
    std::string proj_str = (meta_.map_proj == 0) ? "Cartesian" : "Lambert Conformal";
    put_att_text(ncid, NC_GLOBAL, "MAP_PROJ_CHAR", proj_str);

    put_att_text(ncid, NC_GLOBAL, "MMINLU", "    ");
    put_att_int(ncid, NC_GLOBAL, "NUM_LAND_CAT", 21);
    put_att_int(ncid, NC_GLOBAL, "ISWATER", 16);
    put_att_int(ncid, NC_GLOBAL, "ISLAKE", -1);
    put_att_int(ncid, NC_GLOBAL, "ISICE", 24);
    put_att_int(ncid, NC_GLOBAL, "ISURBAN", 1);
    put_att_int(ncid, NC_GLOBAL, "ISOILWATER", 14);

    put_att_int(ncid, NC_GLOBAL, "HYBRID_OPT", meta_.hybrid_opt);
    put_att_float(ncid, NC_GLOBAL, "ETAC", float(meta_.etac));
  }

  void write_frame_data(int ncid, size_t time_idx, double time_sec, const std::string& date_str,
                        const std::vector<VarView>& vars, const std::map<std::string, VarDef>& defs) {
    const int nx = meta_.nx;
    const int ny = meta_.ny;
    const int nz = meta_.nz;

    // 1. Times
    if (defs.count("Times")) {
      size_t start[2] = {time_idx, 0};
      size_t count[2] = {1, 19};
      char dstr[19];
      std::memset(dstr, ' ', 19);
      std::memcpy(dstr, date_str.data(), std::min(size_t(19), date_str.size()));
      NC_CHECK(nc_put_vara_text(ncid, defs.at("Times").varid, start, count, dstr));
    }

    // Helper lambda to write 3D field
    auto write_3d = [&](std::string_view name, int enx, int eny, int enz) {
      if (!defs.count(std::string(name))) return;
      int varid = defs.at(std::string(name)).varid;
      const VarView* v = find_var(vars, name);
      std::vector<float> buf;
      if (v) {
        pack_to_float(*v, enx, eny, enz, buf);
      } else {
        buf.assign(size_t(enx) * eny * enz, 0.0f);
      }
      size_t start[4] = {time_idx, 0, 0, 0};
      size_t count[4] = {1, size_t(enz), size_t(eny), size_t(enx)};
      NC_CHECK(nc_put_vara_float(ncid, varid, start, count, buf.data()));
    };

    // Helper lambda to write 2D field
    auto write_2d = [&](std::string_view name, int enx, int eny) {
      if (!defs.count(std::string(name))) return;
      int varid = defs.at(std::string(name)).varid;
      const VarView* v = find_var(vars, name);
      std::vector<float> buf;
      if (v) {
        pack_to_float(*v, enx, eny, 1, buf);
      } else {
        buf.assign(size_t(enx) * eny, 0.0f);
      }
      size_t start[3] = {time_idx, 0, 0};
      size_t count[3] = {1, size_t(eny), size_t(enx)};
      NC_CHECK(nc_put_vara_float(ncid, varid, start, count, buf.data()));
    };

    // Helper lambda to write 1D vertical field
    auto write_1d = [&](std::string_view name, int enz, float def_val = 0.0f) {
      if (!defs.count(std::string(name))) return;
      int varid = defs.at(std::string(name)).varid;
      const VarView* v = find_var(vars, name);
      std::vector<float> buf;
      if (v) {
        pack_to_float(*v, 1, 1, enz, buf);
      } else {
        buf.assign(enz, def_val);
      }
      size_t start[2] = {time_idx, 0};
      size_t count[2] = {1, size_t(enz)};
      NC_CHECK(nc_put_vara_float(ncid, varid, start, count, buf.data()));
    };

    // Helper lambda to write 0D scalar field
    auto write_scalar = [&](std::string_view name, float def_val) {
      if (!defs.count(std::string(name))) return;
      int varid = defs.at(std::string(name)).varid;
      const VarView* v = find_var(vars, name);
      float val = def_val;
      if (v && v->data) {
        if (v->dtype == DType::F64) {
          val = static_cast<float>(*static_cast<const double*>(v->data));
        } else {
          val = *static_cast<const float*>(v->data);
        }
      }
      size_t start[1] = {time_idx};
      size_t count[1] = {1};
      NC_CHECK(nc_put_vara_float(ncid, varid, start, count, &val));
    };

    // 3D Staggered
    write_3d("U", nx + 1, ny, nz);
    write_3d("V", nx, ny + 1, nz);
    write_3d("W", nx, ny, nz + 1);
    write_3d("PH", nx, ny, nz + 1);
    write_3d("PHB", nx, ny, nz + 1);

    // 3D Unstaggered
    write_3d("T", nx, ny, nz);
    write_3d("P", nx, ny, nz);
    write_3d("PB", nx, ny, nz);
    if (defs.count("QVAPOR")) {
      write_3d("QVAPOR", nx, ny, nz);
    }
    if (defs.count("QCLOUD")) {
      write_3d("QCLOUD", nx, ny, nz);
    }
    if (defs.count("QRAIN")) {
      write_3d("QRAIN", nx, ny, nz);
    }
    if (defs.count("QICE")) {
      write_3d("QICE", nx, ny, nz);
    }
    if (defs.count("QSNOW")) {
      write_3d("QSNOW", nx, ny, nz);
    }
    if (defs.count("QGRAUP")) {
      write_3d("QGRAUP", nx, ny, nz);
    }

    // 2D Fields
    write_2d("MU", nx, ny);
    write_2d("MUB", nx, ny);
    if (defs.count("RAINNC")) {
      write_2d("RAINNC", nx, ny);
    }
    if (defs.count("SNOWNC")) {
      write_2d("SNOWNC", nx, ny);
    }
    if (defs.count("GRAUPELNC")) {
      write_2d("GRAUPELNC", nx, ny);
    }
    if (defs.count("SR")) {
      write_2d("SR", nx, ny);
    }
    for (const VarView& v : vars) {  // the extra mass-point fields (see the definitions)
      if (!generic_.count(v.name) || !defs.count(v.name)) continue;
      if (v.nz == 1) write_2d(v.name, nx, ny);
      else write_3d(v.name, nx, ny, v.nz);
    }
    write_2d("XLAT", nx, ny);
    write_2d("XLONG", nx, ny);
    write_2d("HGT", nx, ny);

    // 1D Vertical fields
    // Default eta if ZNU/ZNW not provided
    std::vector<float> znw_def(nz + 1), znu_def(nz);
    for (int k = 0; k <= nz; ++k) znw_def[k] = 1.0f - float(k) / float(nz);
    for (int k = 0; k < nz; ++k) znu_def[k] = 0.5f * (znw_def[k] + znw_def[k + 1]);

    const VarView* v_znu = find_var(vars, "ZNU");
    const VarView* v_znw = find_var(vars, "ZNW");
    if (v_znu) write_1d("ZNU", nz);
    else {
      size_t start[2] = {time_idx, 0};
      size_t count[2] = {1, size_t(nz)};
      NC_CHECK(nc_put_vara_float(ncid, defs.at("ZNU").varid, start, count, znu_def.data()));
    }
    if (v_znw) write_1d("ZNW", nz + 1);
    else {
      size_t start[2] = {time_idx, 0};
      size_t count[2] = {1, size_t(nz + 1)};
      NC_CHECK(nc_put_vara_float(ncid, defs.at("ZNW").varid, start, count, znw_def.data()));
    }

    // Hybrid vertical coordinate terms
    write_1d("C1H", nz, 1.0f);
    write_1d("C2H", nz, 0.0f);
    write_1d("C1F", nz + 1, 1.0f);
    write_1d("C2F", nz + 1, 0.0f);
    write_1d("C3H", nz, 0.5f);
    write_1d("C4H", nz, 0.0f);
    write_1d("C3F", nz + 1, 0.5f);
    write_1d("C4F", nz + 1, 0.0f);

    // 0D Scalars
    float rdx_val = float(meta_.dx > 0 ? 1.0 / meta_.dx : 0.0);
    float rdy_val = float(meta_.dy > 0 ? 1.0 / meta_.dy : 0.0);
    write_scalar("RDX", rdx_val);
    write_scalar("RDY", rdy_val);
    write_scalar("P_TOP", float(meta_.p_top));
    write_scalar("XTIME", float(time_sec / 60.0));
  }

  std::string out_dir_;
  WrfoutMeta meta_;
  ScopedNcFile ncid_;
  std::string current_file_;
  size_t time_idx_ = 0;
  std::map<std::string, VarDef> defs_;
  std::set<std::string> generic_;  // extra mass-point fields defined from the frame
  bool is_initialized_ = false;
};

WrfoutWriter::WrfoutWriter(std::string out_dir, WrfoutMeta meta)
    : impl_(std::make_unique<Impl>(std::move(out_dir), std::move(meta))) {}

WrfoutWriter::~WrfoutWriter() = default;
WrfoutWriter::WrfoutWriter(WrfoutWriter&&) noexcept = default;
WrfoutWriter& WrfoutWriter::operator=(WrfoutWriter&&) noexcept = default;

void WrfoutWriter::write_frame(long step, double time_sec, const std::vector<VarView>& vars,
                              const std::string& target_file) {
  impl_->write_frame(step, time_sec, vars, target_file, last_file_);
}

}  // namespace cppwrf::io_netcdf
