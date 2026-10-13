#include "io_netcdf/wrf_readers.hpp"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "io_netcdf/netcdf_common.hpp"

namespace cppwrf::io_netcdf {

namespace {

bool read_dim(int ncid, const char* name, int& len, bool required = false) {
  int dimid = -1;
  if (nc_inq_dimid(ncid, name, &dimid) != NC_NOERR) {
    if (required) throw std::runtime_error("wrf_readers: missing required dimension '" + std::string(name) + "'");
    return false;
  }
  size_t dlen = 0;
  NC_CHECK(nc_inq_dimlen(ncid, dimid, &dlen));
  len = static_cast<int>(dlen);
  if (required && len <= 0) {
    throw std::runtime_error("wrf_readers: dimension '" + std::string(name) + "' has invalid length " + std::to_string(len));
  }
  return true;
}

std::vector<std::string> read_time_strings(int ncid, const char* varname) {
  std::vector<std::string> res;
  int varid = -1;
  if (nc_inq_varid(ncid, varname, &varid) != NC_NOERR) return res;
  int ndims = 0;
  nc_inq_varndims(ncid, varid, &ndims);
  if (ndims < 2) return res;
  std::vector<int> dimids(ndims);
  nc_inq_vardimid(ncid, varid, dimids.data());
  size_t n_times = 0, str_len = 0;
  nc_inq_dimlen(ncid, dimids[0], &n_times);
  nc_inq_dimlen(ncid, dimids[1], &str_len);
  std::vector<char> buf(n_times * str_len);
  size_t start[2] = {0, 0};
  size_t count[2] = {n_times, str_len};
  if (nc_get_vara_text(ncid, varid, start, count, buf.data()) == NC_NOERR) {
    for (size_t t = 0; t < n_times; ++t) {
      std::string s(buf.data() + t * str_len, str_len);
      while (!s.empty() && (s.back() == ' ' || s.back() == '\0')) s.pop_back();
      res.push_back(s);
    }
  }
  return res;
}

}  // namespace

WrfInputData read_wrfinput(const std::string& path) {
  ScopedNcFile ncid;
  int status = nc_open(path.c_str(), NC_NOWRITE, ncid.ptr());
  if (status != NC_NOERR) {
    throw std::runtime_error("read_wrfinput: failed to open '" + path + "': " + nc_strerror(status));
  }

  WrfInputData data;

  // Inquire dimensions
  read_dim(ncid, "west_east", data.nx, true);
  read_dim(ncid, "south_north", data.ny, true);
  read_dim(ncid, "bottom_top", data.nz, true);
  read_dim(ncid, "west_east_stag", data.nx_stag, false);
  read_dim(ncid, "south_north_stag", data.ny_stag, false);
  read_dim(ncid, "bottom_top_stag", data.nz_stag, false);

  if (data.nx_stag == 0 && data.nx > 0) data.nx_stag = data.nx + 1;
  if (data.ny_stag == 0 && data.ny > 0) data.ny_stag = data.ny + 1;
  if (data.nz_stag == 0 && data.nz > 0) data.nz_stag = data.nz + 1;

  // Global attributes
  int natts = 0;
  NC_CHECK(nc_inq_natts(ncid, &natts));
  for (int i = 0; i < natts; ++i) {
    char name[NC_MAX_NAME + 1];
    nc_type type;
    size_t len = 0;
    NC_CHECK(nc_inq_attname(ncid, NC_GLOBAL, i, name));
    NC_CHECK(nc_inq_att(ncid, NC_GLOBAL, name, &type, &len));
    if (type == NC_CHAR) {
      std::vector<char> buf(len + 1, '\0');
      nc_get_att_text(ncid, NC_GLOBAL, name, buf.data());
      data.attrs_str[name] = std::string(buf.data(), len);
    } else if (type == NC_INT || type == NC_SHORT) {
      int val = 0;
      nc_get_att_int(ncid, NC_GLOBAL, name, &val);
      data.attrs_int[name] = val;
    } else if (type == NC_FLOAT || type == NC_DOUBLE) {
      double val = 0.0;
      nc_get_att_double(ncid, NC_GLOBAL, name, &val);
      data.attrs_real[name] = val;
    }
  }

  // Populate common metadata
  if (data.attrs_real.count("DX")) data.dx = data.attrs_real["DX"];
  if (data.attrs_real.count("DY")) data.dy = data.attrs_real["DY"];
  if (data.attrs_real.count("DT")) data.dt = data.attrs_real["DT"];
  if (data.attrs_int.count("HYBRID_OPT")) data.hybrid_opt = data.attrs_int["HYBRID_OPT"];
  if (data.attrs_real.count("ETAC")) data.etac = data.attrs_real["ETAC"];
  if (data.attrs_int.count("MAP_PROJ")) data.map_proj = data.attrs_int["MAP_PROJ"];
  if (data.attrs_real.count("CEN_LAT")) data.cen_lat = data.attrs_real["CEN_LAT"];
  if (data.attrs_real.count("CEN_LON")) data.cen_lon = data.attrs_real["CEN_LON"];
  if (data.attrs_real.count("TRUELAT1")) data.truelat1 = data.attrs_real["TRUELAT1"];
  if (data.attrs_real.count("TRUELAT2")) data.truelat2 = data.attrs_real["TRUELAT2"];
  if (data.attrs_real.count("MOAD_CEN_LAT")) data.moad_cen_lat = data.attrs_real["MOAD_CEN_LAT"];
  if (data.attrs_real.count("STAND_LON")) data.stand_lon = data.attrs_real["STAND_LON"];
  if (data.attrs_str.count("START_DATE")) data.start_date = data.attrs_str["START_DATE"];
  if (data.attrs_str.count("SIMULATION_START_DATE")) data.simulation_start_date = data.attrs_str["SIMULATION_START_DATE"];

  // Helper lambda to read a variable (Time frame 0)
  auto read_var = [&](const char* name, std::vector<double>& out) -> bool {
    int varid = -1;
    if (nc_inq_varid(ncid, name, &varid) != NC_NOERR) return false;
    int ndims = 0;
    nc_inq_varndims(ncid, varid, &ndims);
    std::vector<int> dimids(ndims);
    nc_inq_vardimid(ncid, varid, dimids.data());
    std::vector<size_t> start(ndims, 0);
    std::vector<size_t> count(ndims, 1);
    size_t total = 1;
    for (int d = 0; d < ndims; ++d) {
      size_t dlen = 0;
      nc_inq_dimlen(ncid, dimids[d], &dlen);
      if (d == 0) {
        char dname[NC_MAX_NAME + 1];
        nc_inq_dimname(ncid, dimids[d], dname);
        if (std::string(dname) == "Time") {
          start[d] = 0;
          count[d] = 1;
          continue;
        }
      }
      count[d] = dlen;
      total *= dlen;
    }
    out.resize(total);
    if (nc_get_vara_double(ncid, varid, start.data(), count.data(), out.data()) != NC_NOERR) {
      out.clear();
      return false;
    }
    return true;
  };

  // Helper lambda to read a required variable with exact size verification
  auto read_required_var = [&](const char* name, std::vector<double>& out, size_t expected_size) {
    if (!read_var(name, out) || out.empty()) {
      throw std::runtime_error("read_wrfinput: required variable '" + std::string(name) + "' missing in '" + path + "'");
    }
    if (out.size() != expected_size) {
      throw std::runtime_error("read_wrfinput: variable '" + std::string(name) + "' size mismatch: got " +
                               std::to_string(out.size()) + ", expected " + std::to_string(expected_size));
    }
  };

  const size_t s_1d = size_t(data.nz);
  const size_t s_1d_stag = size_t(data.nz_stag);
  const size_t s_2d = size_t(data.ny) * data.nx;
  const size_t s_2d_u = size_t(data.ny) * data.nx_stag;
  const size_t s_2d_v = size_t(data.ny_stag) * data.nx;
  const size_t s_3d = size_t(data.nz) * data.ny * data.nx;
  const size_t s_3d_u = size_t(data.nz) * data.ny * data.nx_stag;
  const size_t s_3d_v = size_t(data.nz) * data.ny_stag * data.nx;
  const size_t s_3d_w = size_t(data.nz_stag) * data.ny * data.nx;

  // Scalars / 1D constants
  std::vector<double> tmp;
  if (read_var("CF1", tmp) && !tmp.empty()) data.cf1 = tmp[0];
  if (read_var("CF2", tmp) && !tmp.empty()) data.cf2 = tmp[0];
  if (read_var("CF3", tmp) && !tmp.empty()) data.cf3 = tmp[0];
  if (read_var("CFN", tmp) && !tmp.empty()) data.cfn = tmp[0];
  if (read_var("CFN1", tmp) && !tmp.empty()) data.cfn1 = tmp[0];
  if (read_var("P_TOP", tmp) && !tmp.empty()) data.p_top = tmp[0];

  // 1D Vertical fields
  read_required_var("ZNW", data.znw, s_1d_stag);
  read_required_var("ZNU", data.znu, s_1d);
  read_required_var("DN", data.dn, s_1d);
  read_required_var("DNW", data.dnw, s_1d);
  read_required_var("RDN", data.rdn, s_1d);
  read_required_var("RDNW", data.rdnw, s_1d);
  read_required_var("FNM", data.fnm, s_1d);
  read_required_var("FNP", data.fnp, s_1d);

  // Hybrid coordinate weights
  read_required_var("C1H", data.c1h, s_1d);
  read_required_var("C2H", data.c2h, s_1d);
  read_required_var("C3H", data.c3h, s_1d);
  read_required_var("C4H", data.c4h, s_1d);
  read_required_var("C1F", data.c1f, s_1d_stag);
  read_required_var("C2F", data.c2f, s_1d_stag);
  read_required_var("C3F", data.c3f, s_1d_stag);
  read_required_var("C4F", data.c4f, s_1d_stag);

  // 2D Mass fields
  read_required_var("MUB", data.mub, s_2d);
  read_required_var("MU", data.mu, s_2d);
  read_required_var("HGT", data.hgt, s_2d);
  read_required_var("XLAT", data.xlat, s_2d);
  read_required_var("XLONG", data.xlong, s_2d);
  read_required_var("F", data.f, s_2d);
  read_required_var("E", data.e, s_2d);
  read_required_var("SINALPHA", data.sinalpha, s_2d);
  read_required_var("COSALPHA", data.cosalpha, s_2d);

  // Map factors
  if (!read_var("MAPFAC_MX", data.mapfac_mx) || data.mapfac_mx.size() != s_2d) {
    read_required_var("MAPFAC_M", data.mapfac_mx, s_2d);
  }
  if (!read_var("MAPFAC_MY", data.mapfac_my) || data.mapfac_my.size() != s_2d) {
    read_required_var("MAPFAC_M", data.mapfac_my, s_2d);
  }
  if (!read_var("MAPFAC_UX", data.mapfac_ux) || data.mapfac_ux.size() != s_2d_u) {
    read_required_var("MAPFAC_U", data.mapfac_ux, s_2d_u);
  }
  if (!read_var("MAPFAC_UY", data.mapfac_uy) || data.mapfac_uy.size() != s_2d_u) {
    read_required_var("MAPFAC_U", data.mapfac_uy, s_2d_u);
  }
  if (!read_var("MAPFAC_VX", data.mapfac_vx) || data.mapfac_vx.size() != s_2d_v) {
    read_required_var("MAPFAC_V", data.mapfac_vx, s_2d_v);
  }
  if (!read_var("MAPFAC_VY", data.mapfac_vy) || data.mapfac_vy.size() != s_2d_v) {
    read_required_var("MAPFAC_V", data.mapfac_vy, s_2d_v);
  }

  // 3D Mass fields
  read_required_var("PB", data.pb, s_3d);
  read_required_var("P", data.p, s_3d);
  read_required_var("T", data.t, s_3d);
  read_required_var("T_INIT", data.t_init, s_3d);
  read_required_var("ALB", data.alb, s_3d);
  read_required_var("AL", data.al, s_3d);
  read_required_var("QVAPOR", data.qvapor, s_3d);

  // 3D Staggered fields
  read_required_var("U", data.u, s_3d_u);
  read_required_var("V", data.v, s_3d_v);
  read_required_var("W", data.w, s_3d_w);
  read_required_var("PHB", data.phb, s_3d_w);
  read_required_var("PH", data.ph, s_3d_w);

  // Moist species (optional, not duplicated)
  read_var("QCLOUD", data.qcloud);
  read_var("QRAIN", data.qrain);
  read_var("QICE", data.qice);
  read_var("QSNOW", data.qsnow);
  read_var("QGRAUP", data.qgraup);

  // Inquire all 3D moisture variables starting with 'Q'
  int nvars = 0;
  nc_inq_nvars(ncid, &nvars);
  for (int v = 0; v < nvars; ++v) {
    char vname[NC_MAX_NAME + 1];
    nc_inq_varname(ncid, v, vname);
    if (vname[0] == 'Q' && !data.moist_species.count(vname)) {
      int nd = 0;
      nc_inq_varndims(ncid, v, &nd);
      if (nd == 4) {
        std::vector<double> buf;
        if (read_var(vname, buf) && !buf.empty()) {
          data.moist_species[vname] = std::move(buf);
        }
      }
    }
  }
  // Soil dimensions and fields
  read_dim(ncid, "soil_layers_stag", data.num_soil_layers);
  read_var("ZS", data.zs);
  read_var("DZS", data.dzs);
  read_var("TSK", data.tsk);
  read_var("TSLB", data.tslb);
  read_var("SMOIS", data.smois);
  read_var("SH2O", data.sh2o);
  read_var("SNOW", data.snow);
  read_var("SNOWH", data.snowh);
  read_var("SNOWC", data.snowc);
  read_var("CANWAT", data.canwat);
  read_var("ALBEDO", data.albedo);
  read_var("ALBBCK", data.albbck);
  read_var("EMISS", data.emiss);
  read_var("EMBCK", data.embck);
  read_var("SNOALB", data.snoalb);
  read_var("TMN", data.tmn);
  read_var("XLAND", data.xland);
  read_var("SEAICE", data.seaice);
  if (data.seaice.empty()) read_var("XICE", data.seaice);
  read_var("LAKEMASK", data.lakemask);
  read_var("WATER_DEPTH", data.water_depth);
  read_var("VEGFRA", data.vegfra);
  read_var("SHDMIN", data.shdmin);
  read_var("SHDMAX", data.shdmax);
  read_var("LAI", data.lai);
  read_var("LU_INDEX", data.lu_index);
  read_var("IVGTYP", data.ivgtyp);
  read_var("ISLTYP", data.isltyp);
  read_var("U10", data.u10);
  read_var("V10", data.v10);
  read_var("TH2", data.th2);
  read_var("T2", data.t2);
  read_var("Q2", data.q2);
  read_var("UST", data.ust);
  read_var("PBLH", data.pblh);
  read_var("HFX", data.hfx);
  read_var("QFX", data.qfx);
  read_var("LH", data.lh);
  read_var("ZNT", data.znt);
  read_var("QSFC", data.qsfc);

  return data;
}

WrfBdyData read_wrfbdy(const std::string& path) {
  ScopedNcFile ncid;
  int status = nc_open(path.c_str(), NC_NOWRITE, ncid.ptr());
  if (status != NC_NOERR) {
    throw std::runtime_error("read_wrfbdy: failed to open '" + path + "': " + nc_strerror(status));
  }

  WrfBdyData data;

  read_dim(ncid, "Time", data.n_times, true);
  read_dim(ncid, "bdy_width", data.bdy_width, true);
  read_dim(ncid, "west_east", data.nx, true);
  read_dim(ncid, "south_north", data.ny, true);
  read_dim(ncid, "bottom_top", data.nz, true);
  read_dim(ncid, "west_east_stag", data.nx_stag, false);
  read_dim(ncid, "south_north_stag", data.ny_stag, false);
  read_dim(ncid, "bottom_top_stag", data.nz_stag, false);

  if (data.nx_stag == 0 && data.nx > 0) data.nx_stag = data.nx + 1;
  if (data.ny_stag == 0 && data.ny > 0) data.ny_stag = data.ny + 1;
  if (data.nz_stag == 0 && data.nz > 0) data.nz_stag = data.nz + 1;

  data.times = read_time_strings(ncid, "Times");
  data.this_bdy_times = read_time_strings(ncid, "md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_");
  data.next_bdy_times = read_time_strings(ncid, "md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_");

  auto read_bdy_var = [&](const char* name, std::vector<double>& out) -> bool {
    int varid = -1;
    if (nc_inq_varid(ncid, name, &varid) != NC_NOERR) return false;
    int ndims = 0;
    nc_inq_varndims(ncid, varid, &ndims);
    std::vector<int> dimids(ndims);
    nc_inq_vardimid(ncid, varid, dimids.data());
    std::vector<size_t> start(ndims, 0);
    std::vector<size_t> count(ndims, 1);
    size_t total = 1;
    for (int d = 0; d < ndims; ++d) {
      size_t dlen = 0;
      nc_inq_dimlen(ncid, dimids[d], &dlen);
      count[d] = dlen;
      total *= dlen;
    }
    out.resize(total);
    if (nc_get_vara_double(ncid, varid, start.data(), count.data(), out.data()) != NC_NOERR) {
      out.clear();
      return false;
    }
    return true;
  };

  auto read_required_bdy_var = [&](const char* name, std::vector<double>& out, size_t expected_size) {
    if (!read_bdy_var(name, out) || out.empty()) {
      throw std::runtime_error("read_wrfbdy: required variable '" + std::string(name) + "' missing in '" + path + "'");
    }
    if (out.size() != expected_size) {
      throw std::runtime_error("read_wrfbdy: variable '" + std::string(name) + "' size mismatch: got " +
                               std::to_string(out.size()) + ", expected " + std::to_string(expected_size));
    }
  };

  const size_t b_x = size_t(data.n_times) * data.bdy_width * data.nz * data.ny;
  const size_t b_y = size_t(data.n_times) * data.bdy_width * data.nz * data.nx;
  const size_t b_x_u = b_x;
  const size_t b_y_u = size_t(data.n_times) * data.bdy_width * data.nz * data.nx_stag;
  const size_t b_x_v = size_t(data.n_times) * data.bdy_width * data.nz * data.ny_stag;
  const size_t b_y_v = b_y;
  const size_t b_w_x = size_t(data.n_times) * data.bdy_width * data.nz_stag * data.ny;
  const size_t b_w_y = size_t(data.n_times) * data.bdy_width * data.nz_stag * data.nx;
  const size_t b_mu_x = size_t(data.n_times) * data.bdy_width * data.ny;
  const size_t b_mu_y = size_t(data.n_times) * data.bdy_width * data.nx;

  // U
  read_required_bdy_var("U_BXS", data.u_bxs, b_x_u);
  read_required_bdy_var("U_BXE", data.u_bxe, b_x_u);
  read_required_bdy_var("U_BYS", data.u_bys, b_y_u);
  read_required_bdy_var("U_BYE", data.u_bye, b_y_u);
  read_required_bdy_var("U_BTXS", data.u_btxs, b_x_u);
  read_required_bdy_var("U_BTXE", data.u_btxe, b_x_u);
  read_required_bdy_var("U_BTYS", data.u_btys, b_y_u);
  read_required_bdy_var("U_BTYE", data.u_btye, b_y_u);

  // V
  read_required_bdy_var("V_BXS", data.v_bxs, b_x_v);
  read_required_bdy_var("V_BXE", data.v_bxe, b_x_v);
  read_required_bdy_var("V_BYS", data.v_bys, b_y_v);
  read_required_bdy_var("V_BYE", data.v_bye, b_y_v);
  read_required_bdy_var("V_BTXS", data.v_btxs, b_x_v);
  read_required_bdy_var("V_BTXE", data.v_btxe, b_x_v);
  read_required_bdy_var("V_BTYS", data.v_btys, b_y_v);
  read_required_bdy_var("V_BTYE", data.v_btye, b_y_v);

  // W
  read_required_bdy_var("W_BXS", data.w_bxs, b_w_x);
  read_required_bdy_var("W_BXE", data.w_bxe, b_w_x);
  read_required_bdy_var("W_BYS", data.w_bys, b_w_y);
  read_required_bdy_var("W_BYE", data.w_bye, b_w_y);
  read_required_bdy_var("W_BTXS", data.w_btxs, b_w_x);
  read_required_bdy_var("W_BTXE", data.w_btxe, b_w_x);
  read_required_bdy_var("W_BTYS", data.w_btys, b_w_y);
  read_required_bdy_var("W_BTYE", data.w_btye, b_w_y);

  // PH
  read_required_bdy_var("PH_BXS", data.ph_bxs, b_w_x);
  read_required_bdy_var("PH_BXE", data.ph_bxe, b_w_x);
  read_required_bdy_var("PH_BYS", data.ph_bys, b_w_y);
  read_required_bdy_var("PH_BYE", data.ph_bye, b_w_y);
  read_required_bdy_var("PH_BTXS", data.ph_btxs, b_w_x);
  read_required_bdy_var("PH_BTXE", data.ph_btxe, b_w_x);
  read_required_bdy_var("PH_BTYS", data.ph_btys, b_w_y);
  read_required_bdy_var("PH_BTYE", data.ph_btye, b_w_y);

  // T
  read_required_bdy_var("T_BXS", data.t_bxs, b_x);
  read_required_bdy_var("T_BXE", data.t_bxe, b_x);
  read_required_bdy_var("T_BYS", data.t_bys, b_y);
  read_required_bdy_var("T_BYE", data.t_bye, b_y);
  read_required_bdy_var("T_BTXS", data.t_btxs, b_x);
  read_required_bdy_var("T_BTXE", data.t_btxe, b_x);
  read_required_bdy_var("T_BTYS", data.t_btys, b_y);
  read_required_bdy_var("T_BTYE", data.t_btye, b_y);

  // MU
  read_required_bdy_var("MU_BXS", data.mu_bxs, b_mu_x);
  read_required_bdy_var("MU_BXE", data.mu_bxe, b_mu_x);
  read_required_bdy_var("MU_BYS", data.mu_bys, b_mu_y);
  read_required_bdy_var("MU_BYE", data.mu_bye, b_mu_y);
  read_required_bdy_var("MU_BTXS", data.mu_btxs, b_mu_x);
  read_required_bdy_var("MU_BTXE", data.mu_btxe, b_mu_x);
  read_required_bdy_var("MU_BTYS", data.mu_btys, b_mu_y);
  read_required_bdy_var("MU_BTYE", data.mu_btye, b_mu_y);

  // QVAPOR
  read_required_bdy_var("QVAPOR_BXS", data.qvapor_bxs, b_x);
  read_required_bdy_var("QVAPOR_BXE", data.qvapor_bxe, b_x);
  read_required_bdy_var("QVAPOR_BYS", data.qvapor_bys, b_y);
  read_required_bdy_var("QVAPOR_BYE", data.qvapor_bye, b_y);
  read_required_bdy_var("QVAPOR_BTXS", data.qvapor_btxs, b_x);
  read_required_bdy_var("QVAPOR_BTXE", data.qvapor_btxe, b_x);
  read_required_bdy_var("QVAPOR_BTYS", data.qvapor_btys, b_y);
  read_required_bdy_var("QVAPOR_BTYE", data.qvapor_btye, b_y);

  // QCLOUD (optional)
  read_bdy_var("QCLOUD_BXS", data.qcloud_bxs);
  read_bdy_var("QCLOUD_BXE", data.qcloud_bxe);
  read_bdy_var("QCLOUD_BYS", data.qcloud_bys);
  read_bdy_var("QCLOUD_BYE", data.qcloud_bye);
  read_bdy_var("QCLOUD_BTXS", data.qcloud_btxs);
  read_bdy_var("QCLOUD_BTXE", data.qcloud_btxe);
  read_bdy_var("QCLOUD_BTYS", data.qcloud_btys);
  read_bdy_var("QCLOUD_BTYE", data.qcloud_btye);

  // QRAIN (optional)
  read_bdy_var("QRAIN_BXS", data.qrain_bxs);
  read_bdy_var("QRAIN_BXE", data.qrain_bxe);
  read_bdy_var("QRAIN_BYS", data.qrain_bys);
  read_bdy_var("QRAIN_BYE", data.qrain_bye);
  read_bdy_var("QRAIN_BTXS", data.qrain_btxs);
  read_bdy_var("QRAIN_BTXE", data.qrain_btxe);
  read_bdy_var("QRAIN_BTYS", data.qrain_btys);
  read_bdy_var("QRAIN_BTYE", data.qrain_btye);

  // QICE (optional)
  read_bdy_var("QICE_BXS", data.qice_bxs);
  read_bdy_var("QICE_BXE", data.qice_bxe);
  read_bdy_var("QICE_BYS", data.qice_bys);
  read_bdy_var("QICE_BYE", data.qice_bye);
  read_bdy_var("QICE_BTXS", data.qice_btxs);
  read_bdy_var("QICE_BTXE", data.qice_btxe);
  read_bdy_var("QICE_BTYS", data.qice_btys);
  read_bdy_var("QICE_BTYE", data.qice_btye);

  return data;
}

}  // namespace cppwrf::io_netcdf
