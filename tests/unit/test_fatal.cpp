// WRF's fatal stops inside physics kernels (raise_fatal in the kernel, fatal_check on the host):
//  (a) Noah land SFLX, glacial ice and sea ice with SNEQV / SNOWH > 1 stop with
//      "Physical snow depth is less than snow water equiv.";
//      Noah land with a category beyond the tables or NROOT > NSOIL stops as WRF's REDPRM;
//  (b) the Kain-Fritsch codes (QG(1) < 0, ISTOP = 1) raised in a kernel stop with WRF's messages;
//  (c) a physical snowpack raises nothing, and fatal_check clears the flag after a stop.
// The kernels run on the backend of the build (CPU or CUDA device flag).
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "core/precision.hpp"
#include "phys/driver/noah_ice.hpp"
#include "phys/lsm_noah/noah.hpp"

using namespace cppwrf;
static int fails = 0;
#define CHECK(c, ...)                                                                                       \
  do {                                                                                                      \
    if (!(c)) { std::printf("FAIL %s:%d %s ", __FILE__, __LINE__, #c); std::printf(__VA_ARGS__); std::printf("\n"); ++fails; } \
  } while (0)

// Runs f, returns the FatalError code (0: no stop).
template <class F> static int stop_code(F&& f) {
  try {
    f();
    fatal_check();
  } catch (const FatalError& e) {
    return e.code;
  }
  return 0;
}

static const Range2D one{0, 1, 0, 1};

struct KGlacial {
  real sneqv, snowh;
  CPPWRF_HD void operator()(int, int) const {
    real sldpth[4] = {real(0.1), real(0.3), real(0.6), real(1.0)};
    real stc[4] = {real(260.), real(261.), real(262.), real(263.)};
    real emissi = real(0.95), t1 = real(260.), sh = snowh, sw = sneqv, sncovr = real(1.), snotime1 = real(0.),
         ribb = real(0.);
    phys::noah_ice::glacial::Out o{};
    phys::noah_ice::glacial::sflx_glacial(real(0.), real(60.), 4, sldpth, real(250.), real(100.), real(1.e5),
                                          real(0.), real(260.), real(1.e-3), real(260.), real(1.2e-3),
                                          real(1.e-4), real(0.8), real(0.8), real(263.), real(0.001), emissi,
                                          real(0.95), t1, stc, sh, sw, real(0.01), sncovr, snotime1, ribb,
                                          real(2.834e6), o);
  }
};

struct KSeaice {
  real sneqv, snowh;
  CPPWRF_HD void operator()(int, int) const {
    real stc[4] = {real(265.), real(266.), real(268.), real(270.)};
    real z0brd = real(0.001), emissi = real(0.98), t1 = real(265.), sh = snowh, sw = sneqv, sncovr = real(1.),
         ribb = real(0.);
    phys::noah_ice::seaice::Opts op{};
    phys::noah_ice::seaice::Out o{};
    phys::noah_ice::seaice::sflx_seaice(op, real(0.), real(60.), 4, real(3.), real(250.), real(100.), real(1.e5),
                                        real(0.), real(265.), real(1.e-3), real(265.), real(1.2e-3), real(1.e-4),
                                        real(0.8), real(271.4), z0brd, emissi, t1, stc, sh, sw, real(0.01),
                                        sncovr, ribb, o);
  }
};

struct KRaise {
  int code;
  CPPWRF_HD void operator()(int, int) const { raise_fatal(code); }
};

// Noah land SFLX on one column: every real view points into a zeroed scratch block except SNEQV and
// SNOWH; VEGTYP = SOILTYP = SLOPETYP = ACTIVE = 1; one category in each table (cats: SLCATS, LUCATS,
// SLPCATS; nroot: NROTBL(1)).
static int noah_land(real sneqv, real snowh, int cats = 1, int nroot = 1) {
  const int nscr = 64;
  real* scr = static_cast<real*>(dev_alloc(sizeof(real) * nscr));
  real* snow = static_cast<real*>(dev_alloc(sizeof(real) * 2));
  int* ones = static_cast<int*>(dev_alloc(sizeof(int)));
  const real hs[2] = {sneqv, snowh};
  const int one_i = 1;
  copy_to_dev(snow, hs, sizeof hs);
  copy_to_dev(ones, &one_i, sizeof one_i);
  phys::noah::NoahKernel k{};
  View2<real> s{scr, 1};
  View2<const real> c{scr, 1};
  View3<real> s3{scr, 1, 1};
  View2<const int> n{ones, 1};
  k.lwdn = k.soldn = k.solnet = k.sfcprs = k.prcp = k.sfctmp = k.th2 = k.q2 = k.sfcspd = k.q2sat = k.dqsdt2 =
      k.ffrozp = k.shdfac = k.shdmin = k.shdmax = k.alb = k.snoalb = k.tbot = c;
  k.vegtyp = k.soiltyp = k.slopetyp = k.active = n;
  k.z0brd = k.embrd = k.cmc = k.t1 = k.albedo = k.ch = k.cm = k.emissi = k.z0 = k.xlai = k.snotime1 = k.ribb = s;
  k.stc = k.smc = k.sh2o = k.et = k.smav = s3;
  k.eta = k.sheat = k.eta_kinematic = k.fdown = k.ec = k.edir = k.ett = k.esnow = k.drip = k.dew = k.beta = k.etp =
      k.ssoil = k.flx1 = k.flx2 = k.flx3 = k.snomlt = k.sncovr = k.runoff1 = k.runoff2 = k.runoff3 = k.rc = k.pc =
          k.rsmin = k.soilw = k.soilm = k.q1 = s;
  k.sneqv = View2<real>{snow, 1};
  k.snowh = View2<real>{snow + 1, 1};
  k.dt = real(60.);
  k.zlvl = real(10.);
  k.nsoil = 4;
  const real sl[4] = {real(0.1), real(0.3), real(0.6), real(1.0)};
  for (int z = 0; z < 4; ++z) k.sldpth[z] = sl[z];
  k.isurban = 1;
  k.tables.slcats = k.tables.lucats = k.tables.slpcats = cats;
  k.tables.nrotbl[0] = nroot;
  const int code = stop_code([&] { phys::noah::run_noah(k, one); });
  dev_free(scr);
  dev_free(snow);
  dev_free(ones);
  return code;
}

int main() {
  // (a) snow density > 1: 1 m of water in 0.5 m of snow
  CHECK(noah_land(real(1.0), real(0.5)) == kFatalSnowDensity, "noah land");
  CHECK(stop_code([] { column_for(one, KGlacial{real(1.0), real(0.5)}); }) == kFatalSnowDensity, "glacial");
  CHECK(stop_code([] { column_for(one, KSeaice{real(1.0), real(0.5)}); }) == kFatalSnowDensity, "sea ice");
  // Noah REDPRM: a category beyond the tables, more root layers than soil layers
  CHECK(noah_land(real(0.2), real(0.8), 0) == kFatalNoahSoilTypes, "noah soil types");
  CHECK(noah_land(real(0.2), real(0.8), 1, 5) == kFatalNoahRootLayers, "noah root layers");
  // (b) Kain-Fritsch codes and the messages
  CHECK(stop_code([] { column_for(one, KRaise{kFatalKfQgSurface}); }) == kFatalKfQgSurface, "kf qg");
  CHECK(stop_code([] { column_for(one, KRaise{kFatalKfIstop}); }) == kFatalKfIstop, "kf istop");
  CHECK(std::strcmp(fatal_message(kFatalSnowDensity), "Physical snow depth is less than snow water equiv.") == 0, "msg");
  CHECK(std::strcmp(fatal_message(kFatalKfQgSurface), "QG, QG(NK).LT.0") == 0, "msg");
  CHECK(std::strcmp(fatal_message(kFatalKfIstop), "KAIN-FRITSCH, istop=1, diags") == 0, "msg");
  // (c) physical snow (density 0.25) raises nothing; the flag is clear after each stop above
  CHECK(stop_code([] { column_for(one, KGlacial{real(0.2), real(0.8)}); }) == 0, "glacial ok");
  CHECK(stop_code([] { column_for(one, KSeaice{real(0.2), real(0.8)}); }) == 0, "sea ice ok");
  CHECK(stop_code([] {}) == 0, "flag cleared");
  std::printf("%s (%d failures)\n", fails ? "FAILED" : "OK", fails);
  return fails ? 1 : 0;
}
