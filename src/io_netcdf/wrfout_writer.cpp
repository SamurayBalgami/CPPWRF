#include "io_netcdf/wrfout_writer.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <map>
#include <stdexcept>

namespace cppwrf::io_netcdf {

namespace {

std::string to_upper(std::string_view s) {
  std::string r(s);
  for (char& c : r) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return r;
}

const VarView* find_var(const std::vector<VarView>& vars, std::string_view name) {
  std::string target = to_upper(name);
  for (const auto& v : vars) {
    if (to_upper(v.name) == target) return &v;
  }
  // Check known aliases
  if (target == "T") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "TH_PERT" || n == "T_PERT") return &v;
    }
  }
  if (target == "P") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "P_PERT") return &v;
    }
  }
  if (target == "MU") {
    for (const auto& v : vars) {
      std::string n = to_upper(v.name);
      if (n == "MU_PERT") return &v;
    }
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
  int y = 1, m = 1, d = 1, hr = 0, mn = 0, sc = 0;
  if (base_date.size() >= 19) {
    std::sscanf(base_date.c_str(), "%4d-%2d-%2d_%2d:%2d:%2d", &y, &m, &d, &hr, &mn, &sc);
  }
  long total_sec = std::lround(time_sec) + sc + mn * 60 + hr * 3600;
  long days_add = total_sec / 86400;
  total_sec %= 86400;
  if (total_sec < 0) {
    total_sec += 86400;
    days_add -= 1;
  }
  hr = int(total_sec / 3600);
  mn = int((total_sec % 3600) / 60);
  sc = int(total_sec % 60);
  d += int(days_add);

  // Simplified calendar for simulation timestamps
  while (d > 30) {
    d -= 30;
    m++;
    if (m > 12) {
      m = 1;
      y++;
    }
  }
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d_%02d:%02d:%02d", y, m, d, hr, mn, sc);
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

    bool file_exists = std::filesystem::exists(file_path);
    int ncid = -1;
    size_t time_idx = 0;

    int dim_time = -1, dim_datestrlen = -1, dim_we = -1, dim_sn = -1, dim_bt = -1;
    int dim_we_stag = -1, dim_sn_stag = -1, dim_bt_stag = -1;

    std::map<std::string, VarDef> defs;

    if (!file_exists) {
      NC_CHECK(nc_create(file_path.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &ncid));

      // 1. Dimensions
      NC_CHECK(nc_def_dim(ncid, "Time", NC_UNLIMITED, &dim_time));
      NC_CHECK(nc_def_dim(ncid, "DateStrLen", 19, &dim_datestrlen));
      NC_CHECK(nc_def_dim(ncid, "west_east", meta_.nx, &dim_we));
      NC_CHECK(nc_def_dim(ncid, "south_north", meta_.ny, &dim_sn));
      NC_CHECK(nc_def_dim(ncid, "bottom_top", meta_.nz, &dim_bt));
      NC_CHECK(nc_def_dim(ncid, "bottom_top_stag", meta_.nz + 1, &dim_bt_stag));
      NC_CHECK(nc_def_dim(ncid, "west_east_stag", meta_.nx + 1, &dim_we_stag));
      NC_CHECK(nc_def_dim(ncid, "south_north_stag", meta_.ny + 1, &dim_sn_stag));

      // 2. Global attributes
      write_global_attributes(ncid);

      // 3. Define standard variables
      // Times
      defs["Times"] = {"Times", {dim_time, dim_datestrlen}, "", "", "WRF Output Times", "", ""};

      // 3D Staggered U, V, W, PH, PHB
      defs["U"] = {"U", {dim_time, dim_bt, dim_sn, dim_we_stag}, "XYZ", "X", "x-wind component", "m s-1", "XLONG_U XLAT_U XTIME"};
      defs["V"] = {"V", {dim_time, dim_bt, dim_sn_stag, dim_we}, "XYZ", "Y", "y-wind component", "m s-1", "XLONG_V XLAT_V XTIME"};
      defs["W"] = {"W", {dim_time, dim_bt_stag, dim_sn, dim_we}, "XYZ", "Z", "z-wind component", "m s-1", "XLONG XLAT XTIME"};
      defs["PH"] = {"PH", {dim_time, dim_bt_stag, dim_sn, dim_we}, "XYZ", "Z", "perturbation geopotential", "m2 s-2", "XLONG XLAT XTIME"};
      defs["PHB"] = {"PHB", {dim_time, dim_bt_stag, dim_sn, dim_we}, "XYZ", "Z", "base-state geopotential", "m2 s-2", "XLONG XLAT XTIME"};

      // 3D Unstaggered T, P, PB
      defs["T"] = {"T", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "perturbation potential temperature (theta-t0)", "K", "XLONG XLAT XTIME"};
      defs["P"] = {"P", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "perturbation pressure", "Pa", "XLONG XLAT XTIME"};
      defs["PB"] = {"PB", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "BASE STATE PRESSURE", "Pa", "XLONG XLAT XTIME"};

      // Optional moist species and surface diagnostics
      if (find_var(vars, "QVAPOR") || find_var(vars, "QV")) {
        defs["QVAPOR"] = {"QVAPOR", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "Water vapor mixing ratio", "kg kg-1", "XLONG XLAT XTIME"};
      }
      if (find_var(vars, "QCLOUD") || find_var(vars, "QC")) {
        defs["QCLOUD"] = {"QCLOUD", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "Cloud water mixing ratio", "kg kg-1", "XLONG XLAT XTIME"};
      }
      if (find_var(vars, "QRAIN") || find_var(vars, "QR")) {
        defs["QRAIN"] = {"QRAIN", {dim_time, dim_bt, dim_sn, dim_we}, "XYZ", "", "Rain water mixing ratio", "kg kg-1", "XLONG XLAT XTIME"};
      }
      if (find_var(vars, "RAINNC")) {
        defs["RAINNC"] = {"RAINNC", {dim_time, dim_sn, dim_we}, "XY ", "", "ACCUMULATED TOTAL GRID SCALE PRECIPITATION", "mm", "XLONG XLAT XTIME"};
      }

      // 2D MU, MUB, XLAT, XLONG, HGT
      defs["MU"] = {"MU", {dim_time, dim_sn, dim_we}, "XY ", "", "perturbation dry air mass in column", "Pa", "XLONG XLAT XTIME"};
      defs["MUB"] = {"MUB", {dim_time, dim_sn, dim_we}, "XY ", "", "base state dry air mass in column", "Pa", "XLONG XLAT XTIME"};
      defs["XLAT"] = {"XLAT", {dim_time, dim_sn, dim_we}, "XY ", "", "LATITUDE, SOUTH IS NEGATIVE", "degree_north", ""};
      defs["XLONG"] = {"XLONG", {dim_time, dim_sn, dim_we}, "XY ", "", "LONGITUDE, WEST IS NEGATIVE", "degree_east", ""};
      defs["HGT"] = {"HGT", {dim_time, dim_sn, dim_we}, "XY ", "", "Terrain Height", "m", "XLONG XLAT XTIME"};

      // 1D Vertical coordinates
      defs["ZNU"] = {"ZNU", {dim_time, dim_bt}, "Z  ", "", "eta values on half (mass) levels", "", ""};
      defs["ZNW"] = {"ZNW", {dim_time, dim_bt_stag}, "Z  ", "Z", "eta values on full (w) levels", "", ""};

      // Hybrid vertical coordinate terms
      defs["C1H"] = {"C1H", {dim_time, dim_bt}, "Z  ", "", "half levels, c1h = d bf / d eta, using znw", "Dimensionless", ""};
      defs["C2H"] = {"C2H", {dim_time, dim_bt}, "Z  ", "", "half levels, c2h = (1-c1h)*(p0-pt)", "Pa", ""};
      defs["C1F"] = {"C1F", {dim_time, dim_bt_stag}, "Z  ", "Z", "full levels, c1f = d bf / d eta, using znu", "Dimensionless", ""};
      defs["C2F"] = {"C2F", {dim_time, dim_bt_stag}, "Z  ", "Z", "full levels, c2f = (1-c1f)*(p0-pt)", "Pa", ""};
      defs["C3H"] = {"C3H", {dim_time, dim_bt}, "Z  ", "", "half levels, c3h = bh", "Dimensionless", ""};
      defs["C4H"] = {"C4H", {dim_time, dim_bt}, "Z  ", "", "half levels, c4h = (eta-bh)*(p0-pt), using znu", "Pa", ""};
      defs["C3F"] = {"C3F", {dim_time, dim_bt_stag}, "Z  ", "Z", "full levels, c3f = bf", "Dimensionless", ""};
      defs["C4F"] = {"C4F", {dim_time, dim_bt_stag}, "Z  ", "Z", "full levels, c4f = (eta-bf)*(p0-pt), using znw", "Pa", ""};

      // 0D Scalars
      defs["RDX"] = {"RDX", {dim_time}, "0  ", "", "INVERSE X GRID LENGTH", "m-1", ""};
      defs["RDY"] = {"RDY", {dim_time}, "0  ", "", "INVERSE Y GRID LENGTH", "m-1", ""};
      defs["P_TOP"] = {"P_TOP", {dim_time}, "0  ", "", "PRESSURE TOP OF THE MODEL", "Pa", ""};
      defs["XTIME"] = {"XTIME", {dim_time}, "0  ", "", "minutes since simulation start", "minutes", ""};

      for (auto& [name, d] : defs) {
        nc_type xtype = (name == "Times") ? NC_CHAR : NC_FLOAT;
        NC_CHECK(nc_def_var(ncid, name.c_str(), xtype, int(d.dimids.size()), d.dimids.data(), &d.varid));
        if (name != "Times") {
          put_att_int(ncid, d.varid, "FieldType", 104);
          put_att_text(ncid, d.varid, "MemoryOrder", d.memory_order);
          put_att_text(ncid, d.varid, "description", d.desc);
          put_att_text(ncid, d.varid, "units", d.units);
          put_att_text(ncid, d.varid, "stagger", d.stagger);
          if (!d.coords.empty()) {
            put_att_text(ncid, d.varid, "coordinates", d.coords);
          }
        }
      }

      NC_CHECK(nc_enddef(ncid));
      time_idx = 0;
    } else {
      NC_CHECK(nc_open(file_path.c_str(), NC_WRITE, &ncid));
      NC_CHECK(nc_inq_dimid(ncid, "Time", &dim_time));
      size_t tlen = 0;
      NC_CHECK(nc_inq_dimlen(ncid, dim_time, &tlen));
      time_idx = tlen;

      // Populate varids from existing file
      int nvars = 0;
      NC_CHECK(nc_inq_nvars(ncid, &nvars));
      for (int i = 0; i < nvars; ++i) {
        char vname[NC_MAX_NAME + 1];
        NC_CHECK(nc_inq_varname(ncid, i, vname));
        defs[vname] = VarDef{vname, {}, "", "", "", "", "", i};
      }
    }

    // Write data at time_idx
    write_frame_data(ncid, time_idx, time_sec, date_str, vars, defs);
    NC_CHECK(nc_close(ncid));
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
    put_att_int(ncid, NC_GLOBAL, "USE_THETA_M", 0);
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

    // 2D Fields
    write_2d("MU", nx, ny);
    write_2d("MUB", nx, ny);
    if (defs.count("RAINNC")) {
      write_2d("RAINNC", nx, ny);
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
