#include "io_netcdf/wrf_readers.hpp"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "io_netcdf/netcdf_common.hpp"

namespace cppwrf::io_netcdf {

namespace {

bool read_dim(int ncid, const char* name, int& len) {
  int dimid = -1;
  if (nc_inq_dimid(ncid, name, &dimid) != NC_NOERR) return false;
  size_t dlen = 0;
  if (nc_inq_dimlen(ncid, dimid, &dlen) != NC_NOERR) return false;
  len = static_cast<int>(dlen);
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
  int ncid = -1;
  int status = nc_open(path.c_str(), NC_NOWRITE, &ncid);
  if (status != NC_NOERR) {
    throw std::runtime_error("read_wrfinput: failed to open '" + path + "': " + nc_strerror(status));
  }

  WrfInputData data;

  // Inquire dimensions
  read_dim(ncid, "west_east", data.nx);
  read_dim(ncid, "south_north", data.ny);
  read_dim(ncid, "bottom_top", data.nz);
  read_dim(ncid, "west_east_stag", data.nx_stag);
  read_dim(ncid, "south_north_stag", data.ny_stag);
  read_dim(ncid, "bottom_top_stag", data.nz_stag);

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

  // Scalars / 1D constants
  std::vector<double> tmp;
  if (read_var("CF1", tmp) && !tmp.empty()) data.cf1 = tmp[0];
  if (read_var("CF2", tmp) && !tmp.empty()) data.cf2 = tmp[0];
  if (read_var("CF3", tmp) && !tmp.empty()) data.cf3 = tmp[0];
  if (read_var("CFN", tmp) && !tmp.empty()) data.cfn = tmp[0];
  if (read_var("CFN1", tmp) && !tmp.empty()) data.cfn1 = tmp[0];
  if (read_var("P_TOP", tmp) && !tmp.empty()) data.p_top = tmp[0];

  // 1D Vertical fields
  read_var("ZNW", data.znw);
  read_var("ZNU", data.znu);
  read_var("DN", data.dn);
  read_var("DNW", data.dnw);
  read_var("RDN", data.rdn);
  read_var("RDNW", data.rdnw);
  read_var("FNM", data.fnm);
  read_var("FNP", data.fnp);

  // Hybrid coordinate weights
  read_var("C1H", data.c1h);
  read_var("C2H", data.c2h);
  read_var("C3H", data.c3h);
  read_var("C4H", data.c4h);
  read_var("C1F", data.c1f);
  read_var("C2F", data.c2f);
  read_var("C3F", data.c3f);
  read_var("C4F", data.c4f);

  // 2D Mass fields
  read_var("MUB", data.mub);
  read_var("MU", data.mu);
  read_var("HGT", data.hgt);
  read_var("XLAT", data.xlat);
  read_var("XLONG", data.xlong);
  read_var("F", data.f);
  read_var("E", data.e);
  read_var("SINALPHA", data.sinalpha);
  read_var("COSALPHA", data.cosalpha);

  // Map factors
  if (!read_var("MAPFAC_MX", data.mapfac_mx)) read_var("MAPFAC_M", data.mapfac_mx);
  if (!read_var("MAPFAC_MY", data.mapfac_my)) read_var("MAPFAC_M", data.mapfac_my);
  if (!read_var("MAPFAC_UX", data.mapfac_ux)) read_var("MAPFAC_U", data.mapfac_ux);
  if (!read_var("MAPFAC_UY", data.mapfac_uy)) read_var("MAPFAC_U", data.mapfac_uy);
  if (!read_var("MAPFAC_VX", data.mapfac_vx)) read_var("MAPFAC_V", data.mapfac_vx);
  if (!read_var("MAPFAC_VY", data.mapfac_vy)) read_var("MAPFAC_V", data.mapfac_vy);

  // 3D Mass fields
  read_var("PB", data.pb);
  read_var("P", data.p);
  read_var("T", data.t);
  read_var("T_INIT", data.t_init);
  read_var("ALB", data.alb);
  read_var("AL", data.al);
  read_var("QVAPOR", data.qvapor);

  // 3D Staggered fields
  read_var("U", data.u);
  read_var("V", data.v);
  read_var("W", data.w);
  read_var("PHB", data.phb);
  read_var("PH", data.ph);

  // Moist species
  read_var("QCLOUD", data.qcloud);
  read_var("QRAIN", data.qrain);
  read_var("QICE", data.qice);
  read_var("QSNOW", data.qsnow);
  read_var("QGRAUP", data.qgraup);

  if (!data.qcloud.empty()) data.moist_species["QCLOUD"] = data.qcloud;
  if (!data.qrain.empty()) data.moist_species["QRAIN"] = data.qrain;
  if (!data.qice.empty()) data.moist_species["QICE"] = data.qice;
  if (!data.qsnow.empty()) data.moist_species["QSNOW"] = data.qsnow;
  if (!data.qgraup.empty()) data.moist_species["QGRAUP"] = data.qgraup;

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

  nc_close(ncid);
  return data;
}

WrfBdyData read_wrfbdy(const std::string& path) {
  int ncid = -1;
  int status = nc_open(path.c_str(), NC_NOWRITE, &ncid);
  if (status != NC_NOERR) {
    throw std::runtime_error("read_wrfbdy: failed to open '" + path + "': " + nc_strerror(status));
  }

  WrfBdyData data;

  read_dim(ncid, "Time", data.n_times);
  read_dim(ncid, "bdy_width", data.bdy_width);
  read_dim(ncid, "west_east", data.nx);
  read_dim(ncid, "south_north", data.ny);
  read_dim(ncid, "bottom_top", data.nz);
  read_dim(ncid, "west_east_stag", data.nx_stag);
  read_dim(ncid, "south_north_stag", data.ny_stag);
  read_dim(ncid, "bottom_top_stag", data.nz_stag);

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
    data.raw_fields[name] = out;
    return true;
  };

  // U
  read_bdy_var("U_BXS", data.u_bxs);
  read_bdy_var("U_BXE", data.u_bxe);
  read_bdy_var("U_BYS", data.u_bys);
  read_bdy_var("U_BYE", data.u_bye);
  read_bdy_var("U_BTXS", data.u_btxs);
  read_bdy_var("U_BTXE", data.u_btxe);
  read_bdy_var("U_BTYS", data.u_btys);
  read_bdy_var("U_BTYE", data.u_btye);

  // V
  read_bdy_var("V_BXS", data.v_bxs);
  read_bdy_var("V_BXE", data.v_bxe);
  read_bdy_var("V_BYS", data.v_bys);
  read_bdy_var("V_BYE", data.v_bye);
  read_bdy_var("V_BTXS", data.v_btxs);
  read_bdy_var("V_BTXE", data.v_btxe);
  read_bdy_var("V_BTYS", data.v_btys);
  read_bdy_var("V_BTYE", data.v_btye);

  // W
  read_bdy_var("W_BXS", data.w_bxs);
  read_bdy_var("W_BXE", data.w_bxe);
  read_bdy_var("W_BYS", data.w_bys);
  read_bdy_var("W_BYE", data.w_bye);
  read_bdy_var("W_BTXS", data.w_btxs);
  read_bdy_var("W_BTXE", data.w_btxe);
  read_bdy_var("W_BTYS", data.w_btys);
  read_bdy_var("W_BTYE", data.w_btye);

  // PH
  read_bdy_var("PH_BXS", data.ph_bxs);
  read_bdy_var("PH_BXE", data.ph_bxe);
  read_bdy_var("PH_BYS", data.ph_bys);
  read_bdy_var("PH_BYE", data.ph_bye);
  read_bdy_var("PH_BTXS", data.ph_btxs);
  read_bdy_var("PH_BTXE", data.ph_btxe);
  read_bdy_var("PH_BTYS", data.ph_btys);
  read_bdy_var("PH_BTYE", data.ph_btye);

  // T
  read_bdy_var("T_BXS", data.t_bxs);
  read_bdy_var("T_BXE", data.t_bxe);
  read_bdy_var("T_BYS", data.t_bys);
  read_bdy_var("T_BYE", data.t_bye);
  read_bdy_var("T_BTXS", data.t_btxs);
  read_bdy_var("T_BTXE", data.t_btxe);
  read_bdy_var("T_BTYS", data.t_btys);
  read_bdy_var("T_BTYE", data.t_btye);

  // MU
  read_bdy_var("MU_BXS", data.mu_bxs);
  read_bdy_var("MU_BXE", data.mu_bxe);
  read_bdy_var("MU_BYS", data.mu_bys);
  read_bdy_var("MU_BYE", data.mu_bye);
  read_bdy_var("MU_BTXS", data.mu_btxs);
  read_bdy_var("MU_BTXE", data.mu_btxe);
  read_bdy_var("MU_BTYS", data.mu_btys);
  read_bdy_var("MU_BTYE", data.mu_btye);

  // QVAPOR
  read_bdy_var("QVAPOR_BXS", data.qvapor_bxs);
  read_bdy_var("QVAPOR_BXE", data.qvapor_bxe);
  read_bdy_var("QVAPOR_BYS", data.qvapor_bys);
  read_bdy_var("QVAPOR_BYE", data.qvapor_bye);
  read_bdy_var("QVAPOR_BTXS", data.qvapor_btxs);
  read_bdy_var("QVAPOR_BTXE", data.qvapor_btxe);
  read_bdy_var("QVAPOR_BTYS", data.qvapor_btys);
  read_bdy_var("QVAPOR_BTYE", data.qvapor_btye);

  // QCLOUD
  read_bdy_var("QCLOUD_BXS", data.qcloud_bxs);
  read_bdy_var("QCLOUD_BXE", data.qcloud_bxe);
  read_bdy_var("QCLOUD_BYS", data.qcloud_bys);
  read_bdy_var("QCLOUD_BYE", data.qcloud_bye);
  read_bdy_var("QCLOUD_BTXS", data.qcloud_btxs);
  read_bdy_var("QCLOUD_BTXE", data.qcloud_btxe);
  read_bdy_var("QCLOUD_BTYS", data.qcloud_btys);
  read_bdy_var("QCLOUD_BTYE", data.qcloud_btye);

  // QRAIN
  read_bdy_var("QRAIN_BXS", data.qrain_bxs);
  read_bdy_var("QRAIN_BXE", data.qrain_bxe);
  read_bdy_var("QRAIN_BYS", data.qrain_bys);
  read_bdy_var("QRAIN_BYE", data.qrain_bye);
  read_bdy_var("QRAIN_BTXS", data.qrain_btxs);
  read_bdy_var("QRAIN_BTXE", data.qrain_btxe);
  read_bdy_var("QRAIN_BTYS", data.qrain_btys);
  read_bdy_var("QRAIN_BTYE", data.qrain_btye);

  // QICE
  read_bdy_var("QICE_BXS", data.qice_bxs);
  read_bdy_var("QICE_BXE", data.qice_bxe);
  read_bdy_var("QICE_BYS", data.qice_bys);
  read_bdy_var("QICE_BYE", data.qice_bye);
  read_bdy_var("QICE_BTXS", data.qice_btxs);
  read_bdy_var("QICE_BTXE", data.qice_btxe);
  read_bdy_var("QICE_BTYS", data.qice_btys);
  read_bdy_var("QICE_BTYE", data.qice_btye);

  // Generic sweep for any other boundary variable in the file
  int nvars = 0;
  nc_inq_nvars(ncid, &nvars);
  for (int v = 0; v < nvars; ++v) {
    char vname[NC_MAX_NAME + 1];
    nc_inq_varname(ncid, v, vname);
    std::string sname(vname);
    if (!data.raw_fields.count(sname)) {
      if (sname.find("_BX") != std::string::npos || sname.find("_BY") != std::string::npos ||
          sname.find("_BTX") != std::string::npos || sname.find("_BTY") != std::string::npos) {
        std::vector<double> buf;
        read_bdy_var(vname, buf);
      }
    }
  }

  nc_close(ncid);
  return data;
}

}  // namespace cppwrf::io_netcdf
