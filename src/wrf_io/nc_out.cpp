#include "wrf_io/nc_out.hpp"

#include "wrf_io/nc_deflate.hpp"

#include <netcdf.h>

#include <cctype>
#include <cstring>
#include <functional>
#include <stdexcept>

namespace cppwrf::wrf_io {

namespace {
void ck(int st, const std::string& what) {
  if (st != NC_NOERR) throw std::runtime_error("netcdf: " + what + ": " + nc_strerror(st));
}
}  // namespace

NcOut::NcOut(const std::string& path, bool netcdf4) : nc4_(netcdf4) {
  if (netcdf4) deflate::install();
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
  VarInfo info;
  info.id = id;
  info.type = type;
  for (const std::string& dn : dims) info.len.push_back(dimlen_[dn]);
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
    info.chunk = ch;
    info.shuffle = true;
    info.level = 2;
  }
  if (field_type != 0) ck(nc_put_att_int(nc_, id, "FieldType", NC_INT, 1, &field_type), "FieldType");
  for (const char* key : {"MemoryOrder", "description", "units", "stagger", "coordinates"}) {
    auto it = satts.find(key);
    if (it != satts.end()) ck(nc_put_att_text(nc_, id, key, it->second.size(), it->second.c_str()), key);
  }
  vars_[name] = std::move(info);
}

void NcOut::enddef() { ck(nc_enddef(nc_), "enddef"); }

void NcOut::sync() { ck(nc_sync(nc_), "sync"); }

void NcOut::put(const std::string& name, const std::vector<double>& v, std::size_t rec) {
  stage(name, v, rec);
  flush();
}

void NcOut::put_int(const std::string& name, const std::vector<int>& v, std::size_t rec) {
  stage_int(name, v, rec);
  flush();
}

void NcOut::stage(const std::string& name, const std::vector<double>& v, std::size_t rec) {
  if (vars_.at(name).type == 'f')  // converted here, as the float put did
    stage_bytes(name, 'f', v.size(), sizeof(float), rec, [&](unsigned char* p) {
      float* f = reinterpret_cast<float*>(p);
      for (std::size_t i = 0; i < v.size(); ++i) f[i] = float(v[i]);
    });
  else
    stage_bytes(name, 'd', v.size(), sizeof(double), rec,
                [&](unsigned char* p) { std::memcpy(p, v.data(), v.size() * sizeof(double)); });
}

void NcOut::stage_int(const std::string& name, const std::vector<int>& v, std::size_t rec) {
  stage_bytes(name, 'i', v.size(), sizeof(int), rec, [&](unsigned char* p) { std::memcpy(p, v.data(), v.size() * sizeof(int)); });
}

void NcOut::stage_bytes(const std::string& name, char mem, std::size_t n, std::size_t esize, std::size_t rec,
                        const std::function<void(unsigned char*)>& fill_in) {
  const VarInfo& vi = vars_.at(name);
  Staged s;
  s.v = &vi;
  s.name = name;
  s.mem = mem;
  s.rec = rec;
  if (staged_.size() < spare_.size()) s.bytes = std::move(spare_[staged_.size()]);  // the previous record's buffer
  s.bytes.resize(n * esize);
  if (n) fill_in(s.bytes.data());
  // the chunks HDF5 will compress: one record, memory type = file type, every value given
  std::size_t total = 1;
  for (std::size_t d = 1; d < vi.len.size(); ++d) total *= vi.len[d];
  if (!vi.chunk.empty() && vi.len.size() >= 2 && vi.chunk[0] == 1 && mem == vi.type && total == n) {
    deflate::Record r;
    r.data = s.bytes.data();
    r.len.assign(vi.len.begin() + 1, vi.len.end());
    r.chunk.assign(vi.chunk.begin() + 1, vi.chunk.end());
    r.elsize = esize;
    r.shuffle = vi.shuffle;
    r.level = vi.level;
    r.fill.resize(esize);
    // HDF5 fills the rest of an edge chunk with the dataset's fill value: NetCDF's default fill of the type
    if (mem == 'f') { const float f = NC_FILL_FLOAT; std::memcpy(r.fill.data(), &f, esize); }
    else if (mem == 'd') { const double d = NC_FILL_DOUBLE; std::memcpy(r.fill.data(), &d, esize); }
    else { const int i = NC_FILL_INT; std::memcpy(r.fill.data(), &i, esize); }
    deflate::submit(r);
  }
  staged_.push_back(std::move(s));
}

void NcOut::flush() {
  if (staged_.empty()) return;
  deflate::wait();
  for (Staged& s : staged_) {
    const VarInfo& vi = *s.v;
    const std::size_t nd = vi.len.size();
    std::vector<std::size_t> st(nd), ct(nd);
    for (std::size_t d = 0; d < nd; ++d) {
      st[d] = d == 0 ? s.rec : 0;
      ct[d] = d == 0 ? 1 : vi.len[d];
    }
    if (s.mem == 'f')
      ck(nc_put_vara_float(nc_, vi.id, st.data(), ct.data(), reinterpret_cast<const float*>(s.bytes.data())), "put " + s.name);
    else if (s.mem == 'd')
      ck(nc_put_vara_double(nc_, vi.id, st.data(), ct.data(), reinterpret_cast<const double*>(s.bytes.data())), "put " + s.name);
    else
      ck(nc_put_vara_int(nc_, vi.id, st.data(), ct.data(), reinterpret_cast<const int*>(s.bytes.data())), "put " + s.name);
  }
  // the buffers stay for the next record (the same variables in the same order: no new pages)
  if (spare_.size() < staged_.size()) spare_.resize(staged_.size());
  for (std::size_t k = 0; k < staged_.size(); ++k) spare_[k] = std::move(staged_[k].bytes);
  staged_.clear();
}

void NcOut::put_text(const std::string& name, const std::string& s, std::size_t rec) {
  int id = vars_.at(name).id;
  std::size_t st[2] = {rec, 0}, ct[2] = {1, s.size()};
  ck(nc_put_vara_text(nc_, id, st, ct, s.c_str()), "put " + name);
}

}  // namespace cppwrf::wrf_io
