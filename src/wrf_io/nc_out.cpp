#include "wrf_io/nc_out.hpp"

#include <netcdf.h>

#include <cctype>
#include <stdexcept>

namespace cppwrf::wrf_io {

namespace {
void ck(int st, const std::string& what) {
  if (st != NC_NOERR) throw std::runtime_error("netcdf: " + what + ": " + nc_strerror(st));
}
}  // namespace

NcOut::NcOut(const std::string& path, bool netcdf4) : nc4_(netcdf4) {
  ck(nc_create(path.c_str(), netcdf4 ? (NC_CLOBBER | NC_NETCDF4) : (NC_CLOBBER | NC_64BIT_OFFSET), &nc_), "create " + path);
}

NcOut::~NcOut() {
  if (nc_ >= 0) nc_close(nc_);
}

void NcOut::dim(const std::string& name, std::size_t len) {
  int id;
  ck(nc_def_dim(nc_, name.c_str(), len == 0 ? NC_UNLIMITED : len, &id), "dim " + name);
  dims_[name] = id;
  dimlen_[name] = len == 0 ? 1 : len;
}

void NcOut::att(const std::string& name, const std::string& v) {
  ck(nc_put_att_text(nc_, NC_GLOBAL, name.c_str(), v.size(), v.c_str()), "att " + name);
}
void NcOut::att(const std::string& name, int v) { ck(nc_put_att_int(nc_, NC_GLOBAL, name.c_str(), NC_INT, 1, &v), "att " + name); }
void NcOut::att_f(const std::string& name, float v) {
  ck(nc_put_att_float(nc_, NC_GLOBAL, name.c_str(), NC_FLOAT, 1, &v), "att " + name);
}
void NcOut::att_d(const std::string& name, double v) {
  ck(nc_put_att_double(nc_, NC_GLOBAL, name.c_str(), NC_DOUBLE, 1, &v), "att " + name);
}
void NcOut::att_f(const std::string& name, const std::vector<float>& v) {
  ck(nc_put_att_float(nc_, NC_GLOBAL, name.c_str(), NC_FLOAT, v.size(), v.data()), "att " + name);
}

void NcOut::var(const std::string& name, char type, const std::vector<std::string>& dims,
                const std::map<std::string, std::string>& satts, int field_type) {
  std::vector<int> ids;
  for (const std::string& dn : dims) {
    auto it = dims_.find(dn);
    if (it == dims_.end()) throw std::runtime_error("netcdf: undefined dimension " + dn + " for " + name);
    ids.push_back(it->second);
  }
  nc_type t = type == 'f' ? NC_FLOAT : type == 'd' ? NC_DOUBLE : type == 'i' ? NC_INT : NC_CHAR;
  int id;
  ck(nc_def_var(nc_, name.c_str(), t, int(ids.size()), ids.data(), &id), "var " + name);
  // external/io_netcdf ext_ncd_write_field: chunk the first two (Fortran-order) dimensions in halves,
  // shuffle + deflate level 2, for memory orders xyz/xy/xsz/xez/ysz/yez
  auto mo = satts.find("MemoryOrder");
  std::string m = mo == satts.end() ? "" : mo->second;
  while (!m.empty() && m.back() == ' ') m.pop_back();
  for (char& ch : m) ch = char(std::tolower(static_cast<unsigned char>(ch)));
  const bool chunk = m == "xyz" || m == "xzy" || m == "yxz" || m == "yzx" || m == "zxy" || m == "zyx" || m == "xsz" ||
                     m == "xez" || m == "ysz" || m == "yez" || m == "xy" || m == "yx";
  if (nc4_ && chunk) {
    const std::size_t nd = dims.size();  // includes Time first
    std::vector<std::size_t> len(nd);
    for (std::size_t k = 0; k < nd; ++k) len[k] = dimlen_[dims[k]];
    // Fortran order: L(1) = fastest = dims[nd-1], ..., time = dims[0]
    std::vector<std::size_t> ch(nd);
    for (std::size_t k = 0; k < nd; ++k) ch[k] = len[k];
    ch[0] = 1;
    std::size_t& c1 = ch[nd - 1];
    std::size_t& c2 = ch[nd - 2];
    c1 = (len[nd - 1] + 1) / 2;
    c2 = (len[nd - 2] + 1) / 2;
    auto block = [&] {
      std::size_t b = 1;
      for (std::size_t k = 1; k < nd; ++k) b *= ch[k];
      return b;
    };
    while (block() > 32000000) {
      c1 = (c1 + 1) / 2;
      c2 = (c2 + 1) / 2;
    }
    ck(nc_def_var_chunking(nc_, id, NC_CHUNKED, ch.data()), "chunking " + name);
    ck(nc_def_var_deflate(nc_, id, 1, 1, 2), "deflate " + name);
  }
  if (field_type != 0) ck(nc_put_att_int(nc_, id, "FieldType", NC_INT, 1, &field_type), "FieldType");
  for (const char* key : {"MemoryOrder", "description", "units", "stagger", "coordinates"}) {
    auto it = satts.find(key);
    if (it != satts.end()) ck(nc_put_att_text(nc_, id, key, it->second.size(), it->second.c_str()), key);
  }
  vars_[name] = id;
  types_[name] = type;
}

void NcOut::enddef() { ck(nc_enddef(nc_), "enddef"); }

void NcOut::sync() { ck(nc_sync(nc_), "sync"); }

void NcOut::put(const std::string& name, const std::vector<double>& v, std::size_t rec) {
  int id = vars_.at(name);
  if (types_[name] == 'f') {
    std::vector<float> f(v.begin(), v.end());
    std::vector<std::size_t> st, ct;
    int nd;
    ck(nc_inq_varndims(nc_, id, &nd), "ndims");
    std::vector<int> dimids(static_cast<std::size_t>(nd));
    ck(nc_inq_vardimid(nc_, id, dimids.data()), "dimids");
    for (int d = 0; d < nd; ++d) {
      std::size_t len;
      ck(nc_inq_dimlen(nc_, dimids[std::size_t(d)], &len), "dimlen");
      st.push_back(d == 0 ? rec : 0);
      ct.push_back(d == 0 ? 1 : len);
    }
    ck(nc_put_vara_float(nc_, id, st.data(), ct.data(), f.data()), "put " + name);
  } else {
    std::vector<std::size_t> st, ct;
    int nd;
    ck(nc_inq_varndims(nc_, id, &nd), "ndims");
    std::vector<int> dimids(static_cast<std::size_t>(nd));
    ck(nc_inq_vardimid(nc_, id, dimids.data()), "dimids");
    for (int d = 0; d < nd; ++d) {
      std::size_t len;
      ck(nc_inq_dimlen(nc_, dimids[std::size_t(d)], &len), "dimlen");
      st.push_back(d == 0 ? rec : 0);
      ct.push_back(d == 0 ? 1 : len);
    }
    ck(nc_put_vara_double(nc_, id, st.data(), ct.data(), v.data()), "put " + name);
  }
}

void NcOut::put_int(const std::string& name, const std::vector<int>& v, std::size_t rec) {
  int id = vars_.at(name);
  int nd;
  ck(nc_inq_varndims(nc_, id, &nd), "ndims");
  std::vector<int> dimids(static_cast<std::size_t>(nd));
  ck(nc_inq_vardimid(nc_, id, dimids.data()), "dimids");
  std::vector<std::size_t> st, ct;
  for (int d = 0; d < nd; ++d) {
    std::size_t len;
    ck(nc_inq_dimlen(nc_, dimids[std::size_t(d)], &len), "dimlen");
    st.push_back(d == 0 ? rec : 0);
    ct.push_back(d == 0 ? 1 : len);
  }
  ck(nc_put_vara_int(nc_, id, st.data(), ct.data(), v.data()), "put " + name);
}

void NcOut::put_text(const std::string& name, const std::string& s, std::size_t rec) {
  int id = vars_.at(name);
  std::size_t st[2] = {rec, 0}, ct[2] = {1, s.size()};
  ck(nc_put_vara_text(nc_, id, st, ct, s.c_str()), "put " + name);
}

}  // namespace cppwrf::wrf_io
