#include "init_real/met_em.hpp"

#include <netcdf.h>

#include <stdexcept>

namespace cppwrf::init_real {

namespace {
void ck(int st, const std::string& what) {
  if (st != NC_NOERR) throw std::runtime_error("met_em: " + what + ": " + nc_strerror(st));
}
}  // namespace

MetEm MetEm::read(const std::string& path) {
  MetEm m;
  m.path = path;
  int nc = -1;
  ck(nc_open(path.c_str(), NC_NOWRITE, &nc), "open " + path);
  int natt = 0, nvar = 0;
  ck(nc_inq_natts(nc, &natt), "natts");
  for (int a = 0; a < natt; ++a) {
    char name[NC_MAX_NAME + 1];
    ck(nc_inq_attname(nc, NC_GLOBAL, a, name), "attname");
    nc_type t;
    std::size_t len = 0;
    ck(nc_inq_att(nc, NC_GLOBAL, name, &t, &len), "att");
    if (t == NC_CHAR) {
      std::string s(len, '\0');
      ck(nc_get_att_text(nc, NC_GLOBAL, name, s.data()), "att text");
      m.satt_[name] = s;
    } else if (t == NC_INT || t == NC_SHORT || t == NC_BYTE) {
      std::vector<int> v(len);
      ck(nc_get_att_int(nc, NC_GLOBAL, name, v.data()), "att int");
      m.iatt_[name] = v;
    } else {
      std::vector<double> v(len);
      ck(nc_get_att_double(nc, NC_GLOBAL, name, v.data()), "att real");
      m.ratt_[name] = v;
    }
  }
  ck(nc_inq_nvars(nc, &nvar), "nvars");
  for (int id = 0; id < nvar; ++id) {
    char name[NC_MAX_NAME + 1];
    nc_type t;
    int nd = 0, dimids[NC_MAX_VAR_DIMS], na = 0;
    ck(nc_inq_var(nc, id, name, &t, &nd, dimids, &na), "var");
    std::vector<std::size_t> start(static_cast<std::size_t>(nd), 0), count(static_cast<std::size_t>(nd), 1);
    MetVar mv;
    std::size_t total = 1;
    for (int d = 0; d < nd; ++d) {
      char dn[NC_MAX_NAME + 1];
      std::size_t len = 0;
      ck(nc_inq_dim(nc, dimids[d], dn, &len), "dim");
      if (d == 0 && std::string(dn) == "Time") continue;  // first time only
      count[std::size_t(d)] = len;
      mv.dims.emplace_back(dn);
      mv.shape.push_back(len);
      total *= len;
    }
    if (std::string(name) == "Times") {
      std::string s(total, '\0');
      ck(nc_get_vara_text(nc, id, start.data(), count.data(), s.data()), "Times");
      m.time = s.substr(0, 19);
      continue;
    }
    if (t == NC_CHAR) continue;
    mv.v.resize(total);
    ck(nc_get_vara_double(nc, id, start.data(), count.data(), mv.v.data()), std::string("read ") + name);
    m.vars_[name] = std::move(mv);
  }
  nc_close(nc);
  return m;
}

const MetVar& MetEm::var(const std::string& name) const {
  auto it = vars_.find(name);
  if (it == vars_.end()) throw std::runtime_error("met_em: no variable " + name + " in " + path);
  return it->second;
}

int MetEm::att_int(const std::string& name, int def) const {
  auto it = iatt_.find(name);
  if (it != iatt_.end() && !it->second.empty()) return it->second[0];
  auto r = ratt_.find(name);
  if (r != ratt_.end() && !r->second.empty()) return int(r->second[0]);
  return def;
}

double MetEm::att_real(const std::string& name, double def) const {
  auto r = ratt_.find(name);
  if (r != ratt_.end() && !r->second.empty()) return r->second[0];
  auto it = iatt_.find(name);
  if (it != iatt_.end() && !it->second.empty()) return it->second[0];
  return def;
}

std::string MetEm::att_str(const std::string& name, const std::string& def) const {
  auto it = satt_.find(name);
  return it == satt_.end() ? def : it->second;
}

}  // namespace cppwrf::init_real
