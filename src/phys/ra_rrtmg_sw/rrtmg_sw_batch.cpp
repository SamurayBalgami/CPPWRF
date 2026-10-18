#include "phys/ra_rrtmg_sw/rrtmg_sw_batch.hpp"

#include <algorithm>
#include <stdexcept>

#include "core/backend.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"

namespace cppwrf::phys {

namespace {

// chunk column ci on threadIdx.x (CUDA) or as the outer collapsed loop (OpenMP)
#if defined(CPPWRF_USE_CUDA)
#define SW_CI(a, c) (a)
#define SW_X(b, c) (c)
inline Range3D sw_range(int n, int nx) { return Range3D{0, n, 0, 1, 0, nx}; }
#else
#define SW_CI(a, c) (c)
#define SW_X(b, c) (b)
inline Range3D sw_range(int n, int nx) { return Range3D{0, 1, 0, nx, 0, n}; }
#endif

struct SwChunk {
  SwLayout lay;
  double* d;
  int* ip;
  int n;
  int col0;
  CPPWRF_FORCEINLINE SwViews views(int ci) const { return lay.views(d, ip, ci, n); }
};

CPPWRF_FORCEINLINE SD col_view(const double* a, int col, int ncol) { return SD{const_cast<double*>(a) + col, ncol}; }

struct KSwPrep {
  SwChunk ch;
  RrtmgSwCallConfig c;
  SwBatchInput in;
  const RrtmgSwTables* t;
  const RrtmgLwConsts* k;
  CPPWRF_HD void operator()(int a, int, int cc) const {
    const int ci = SW_CI(a, cc);
    const int col = ch.col0 + ci, nc = in.ncol;
    SwColumnIn x{};
    x.coszen = in.coszen[col]; x.albedo = in.albedo[col]; x.xland = in.xland[col]; x.xice = in.xice[col];
    x.snow = in.snow[col];
    x.p8w = col_view(in.p8w, col, nc); x.t8w = col_view(in.t8w, col, nc); x.p3d = col_view(in.p3d, col, nc);
    x.pi3d = col_view(in.pi3d, col, nc); x.t3d = col_view(in.t3d, col, nc); x.qv3d = col_view(in.qv, col, nc);
    x.cldfra3d = col_view(in.cldfra, col, nc); x.o33d = col_view(in.o3, col, nc);
    x.has_qc = in.qc != nullptr; x.has_qi = in.qi != nullptr; x.has_qs = in.qs != nullptr;
    x.qc3d = x.has_qc ? col_view(in.qc, col, nc) : SD{nullptr, 0};
    x.qi3d = x.has_qi ? col_view(in.qi, col, nc) : SD{nullptr, 0};
    x.qs3d = x.has_qs ? col_view(in.qs, col, nc) : SD{nullptr, 0};
    x.re_cloud = in.re_cloud ? col_view(in.re_cloud, col, nc) : SD{nullptr, 0};
    x.re_ice = in.re_ice ? col_view(in.re_ice, col, nc) : SD{nullptr, 0};
    x.re_snow = in.re_snow ? col_view(in.re_snow, col, nc) : SD{nullptr, 0};
    rrtmg_sw_prep(c, x, ch.views(ci), *t, *k);
  }
};

template <int ML>
struct KSwGpoint {
  SwChunk ch;
  RrtmgSwCallConfig c;
  const RrtmgSwTables* t;
  CPPWRF_HD void operator()(int a, int b, int cc) const {
    rrtmg_sw_gpoint<ML>(c, SW_X(b, cc), ch.views(SW_CI(a, cc)), *t);
  }
};

struct KSwLevel {
  SwChunk ch;
  const RrtmgSwTables* t;
  CPPWRF_HD void operator()(int a, int b, int cc) const { rrtmg_sw_level(SW_X(b, cc), ch.views(SW_CI(a, cc)), *t); }
};

struct KSwFinish {
  SwChunk ch;
  RrtmgSwCallConfig c;
  SwBatchInput in;
  SwBatchOutput out;
  CPPWRF_HD void operator()(int a, int, int cc) const {
    const int ci = SW_CI(a, cc);
    const int col = ch.col0 + ci, nc = in.ncol;
    const SwColumnOut o =
        rrtmg_sw_finish(c, ch.views(ci), col_view(in.pi3d, col, nc), SD{out.rthraten + col, nc},
                        out.rthratenc ? SD{out.rthratenc + col, nc} : SD{nullptr, 0});
    out.gsw[col] = o.gsw;
    if (out.cf) out.cf[col] = o.swcf;
    if (out.diag) {
      const double d[8] = {o.swupt, o.swuptc, o.swdnt, o.swdntc, o.swupb, o.swupbc, o.swdnb, o.swdnbc};
      for (int q = 0; q < 8; ++q) out.diag[std::size_t(q) * nc + col] = d[q];
    }
  }
};

}  // namespace

RrtmgSwBatch::RrtmgSwBatch(const RrtmgSwTables& tables, int nlayers, int max_chunk)
    : lay_(nlayers), chunk_(std::max(1, max_chunk)) {
  if (nlayers > RRTMG_SW_COL_MAX) throw std::runtime_error("RRTMG SW: too many layers");
  d_tables_ = static_cast<RrtmgSwTables*>(dev_alloc(sizeof(RrtmgSwTables)));
  copy_to_dev(d_tables_, &tables, sizeof(RrtmgSwTables));
  d_consts_ = static_cast<RrtmgLwConsts*>(dev_alloc(sizeof(RrtmgLwConsts)));
  copy_to_dev(d_consts_, &rrtmg_lw_consts(), sizeof(RrtmgLwConsts));
  d_ = static_cast<double*>(dev_alloc(lay_.nd * std::size_t(chunk_) * sizeof(double)));
  i_ = static_cast<int*>(dev_alloc(lay_.ni * std::size_t(chunk_) * sizeof(int)));
}

RrtmgSwBatch::~RrtmgSwBatch() {
  dev_free(d_tables_);
  dev_free(d_consts_);
  dev_free(d_);
  dev_free(i_);
}

void RrtmgSwBatch::run(const RrtmgSwCallConfig& c, const SwBatchInput& in, const SwBatchOutput& out) {
  if (c.nlayers != lay_.nl) throw std::runtime_error("RrtmgSwBatch: nlayers differs from the constructor's");
  for (int col0 = 0; col0 < in.ncol; col0 += chunk_) {
    const int n = std::min(chunk_, in.ncol - col0);
    const SwChunk ch{lay_, d_, i_, n, col0};
    parallel_for(sw_range(n, 1), KSwPrep{ch, c, in, d_tables_, d_consts_});
    if (c.nlayers <= 48)
      parallel_for(sw_range(n, 112), KSwGpoint<48>{ch, c, d_tables_});
    else
      parallel_for(sw_range(n, 112), KSwGpoint<RRTMG_SW_COL_MAX>{ch, c, d_tables_});
    parallel_for(sw_range(n, c.nlayers + 1), KSwLevel{ch, d_tables_});
    parallel_for(sw_range(n, 1), KSwFinish{ch, c, in, out});
  }
}

}  // namespace cppwrf::phys
