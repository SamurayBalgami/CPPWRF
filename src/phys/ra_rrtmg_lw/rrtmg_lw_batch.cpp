#include "phys/ra_rrtmg_lw/rrtmg_lw_batch.hpp"

#include <algorithm>
#include <stdexcept>

#include "core/backend.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"

namespace cppwrf::phys {

namespace {

// Kernel index mapping: chunk column ci and a second index x. CUDA: ci on threadIdx.x (coalesced
// column-fastest arrays); OpenMP: ci and x are the two collapsed parallel loops.
#if defined(CPPWRF_USE_CUDA)
#define LW_CI(a, b, c) (a)
#define LW_X(a, b, c) (c)
inline Range3D lw_range(int n, int nx) { return Range3D{0, n, 0, 1, 0, nx}; }
#else
#define LW_CI(a, b, c) (c)
#define LW_X(a, b, c) (b)
inline Range3D lw_range(int n, int nx) { return Range3D{0, 1, 0, nx, 0, n}; }
#endif

// common part of the kernels: the views of chunk column ci (global column col0 + ci)
struct LwChunk {
  LwLayout lay;
  double* d;
  int* ip;
  int n;     // columns in this chunk
  int col0;  // first global column
  CPPWRF_FORCEINLINE LwViews views(int ci) const { return lay.views(d, ip, ci, n); }
};

struct KLwPrep {
  LwChunk ch;
  RrtmgLwCallConfig c;
  LwBatchInput in;
  const RrtmgLwTables* t;
  const RrtmgLwConsts* k;
  CPPWRF_HD void operator()(int a, int b, int cc) const {
    const int ci = LW_CI(a, b, cc);
    const int col = ch.col0 + ci, nc = in.ncol;
    LwColumnIn x{};
    x.tsk = in.tsk[col]; x.emiss = in.emiss[col]; x.xland = in.xland[col]; x.xice = in.xice[col];
    x.snow = in.snow[col]; x.xlat = 0.0;
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
    rrtmg_lw_prep(c, x, ch.views(ci), *t, *k);
  }
};

struct KLwLayer {
  LwChunk ch;
  int nlayers;
  const RrtmgLwTables* t;
  CPPWRF_HD void operator()(int a, int b, int cc) const {
    rrtmg_lw_layer(nlayers, LW_X(a, b, cc), ch.views(LW_CI(a, b, cc)), *t);
  }
};

template <int MAX_LAY>
struct KLwBand {
  LwChunk ch;
  int nlayers;
  const RrtmgLwTables* t;
  CPPWRF_HD void operator()(int a, int b, int cc) const {
    rrtmg_lw_bandstage<MAX_LAY>(nlayers, LW_X(a, b, cc), ch.views(LW_CI(a, b, cc)), *t);
  }
};

struct KLwLevel {
  LwChunk ch;
  const RrtmgLwTables* t;
  CPPWRF_HD void operator()(int a, int b, int cc) const {
    rrtmg_lw_level(LW_X(a, b, cc), ch.views(LW_CI(a, b, cc)), *t);
  }
};

struct KLwFinish {
  LwChunk ch;
  RrtmgLwCallConfig c;
  LwBatchInput in;
  LwBatchOutput out;
  CPPWRF_HD void operator()(int a, int b, int cc) const {
    const int ci = LW_CI(a, b, cc);
    const int col = ch.col0 + ci, nc = in.ncol;
    const LwColumnOut o = rrtmg_lw_finish(c, ch.views(ci), col_view(in.pi3d, col, nc),
                                          SD{out.rthraten + col, nc},
                                          out.rthratenc ? SD{out.rthratenc + col, nc} : SD{nullptr, 0});
    out.glw[col] = o.glw;
    if (out.cf) out.cf[col] = o.lwcf;
    if (out.olr) out.olr[col] = o.olr;
    if (out.diag) {
      const double d[8] = {o.lwupt, o.lwuptc, o.lwdnt, o.lwdntc, o.lwupb, o.lwupbc, o.lwdnb, o.lwdnbc};
      for (int q = 0; q < 8; ++q) out.diag[std::size_t(q) * nc + col] = d[q];
    }
  }
};

}  // namespace

RrtmgLwBatch::RrtmgLwBatch(const RrtmgLwTables& tables, int nlayers, int max_chunk)
    : lay_(nlayers), chunk_(std::max(1, max_chunk)) {
  if (nlayers > RRTMG_LW_SWEEP_MAX) throw std::runtime_error("RRTMG LW: too many layers");
  d_tables_ = static_cast<RrtmgLwTables*>(dev_alloc(sizeof(RrtmgLwTables)));
  copy_to_dev(d_tables_, &tables, sizeof(RrtmgLwTables));
  d_consts_ = static_cast<RrtmgLwConsts*>(dev_alloc(sizeof(RrtmgLwConsts)));
  copy_to_dev(d_consts_, &rrtmg_lw_consts(), sizeof(RrtmgLwConsts));
  d_ = static_cast<double*>(dev_alloc(lay_.nd * std::size_t(chunk_) * sizeof(double)));
  i_ = static_cast<int*>(dev_alloc(lay_.ni * std::size_t(chunk_) * sizeof(int)));
}

RrtmgLwBatch::~RrtmgLwBatch() {
  dev_free(d_tables_);
  dev_free(d_consts_);
  dev_free(d_);
  dev_free(i_);
}

void RrtmgLwBatch::run(const RrtmgLwCallConfig& c, const LwBatchInput& in, const LwBatchOutput& out) {
  if (c.nlayers != lay_.nl) throw std::runtime_error("RrtmgLwBatch: nlayers differs from the constructor's");
  for (int col0 = 0; col0 < in.ncol; col0 += chunk_) {
    const int n = std::min(chunk_, in.ncol - col0);
    const LwChunk ch{lay_, d_, i_, n, col0};
    parallel_for(lw_range(n, 1), KLwPrep{ch, c, in, d_tables_, d_consts_});
    parallel_for(lw_range(n, c.nlayers), KLwLayer{ch, c.nlayers, d_tables_});
    if (c.nlayers <= 60)
      parallel_for(lw_range(n, 16), KLwBand<60>{ch, c.nlayers, d_tables_});
    else
      parallel_for(lw_range(n, 16), KLwBand<RRTMG_LW_SWEEP_MAX>{ch, c.nlayers, d_tables_});
    parallel_for(lw_range(n, c.nlayers + 1), KLwLevel{ch, d_tables_});
    parallel_for(lw_range(n, 1), KLwFinish{ch, c, in, out});
  }
}

}  // namespace cppwrf::phys
