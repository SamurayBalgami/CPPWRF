// Column driver for noah_ice::glacial::sflx_glacial and noah_ice::seaice::sflx_seaice; reads the
// input of generate_inputs.py and writes the same output layout as driver_noah_ice.f90.
#include <cstdint>
#include <cstdio>
#include <vector>

#include "phys/driver/noah_ice.hpp"

using namespace cppwrf;
using namespace cppwrf::phys::noah_ice;

namespace {
constexpr int kNin = 24, kNg = 27, kNs = 24;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s input.bin output.bin\n", argv[0]);
    return 1;
  }
  std::FILE* f = std::fopen(argv[1], "rb");
  if (!f) return 1;
  int32_t hdr[2];
  double dl[2];
  if (std::fread(hdr, sizeof(int32_t), 2, f) != 2 || std::fread(dl, sizeof(double), 2, f) != 2) return 1;
  const int ncols = hdr[0], nsoil = hdr[1];
  if (nsoil > kMaxSoilLevels) return 1;
  std::vector<double> sl(nsoil), vin(size_t(ncols) * kNin), stcin(size_t(ncols) * nsoil);
  if (std::fread(sl.data(), sizeof(double), sl.size(), f) != sl.size() ||
      std::fread(vin.data(), sizeof(double), vin.size(), f) != vin.size() ||
      std::fread(stcin.data(), sizeof(double), stcin.size(), f) != stcin.size())
    return 1;
  std::fclose(f);
  const real dt = real(dl[0]), lvcoef = real(dl[1]);
  real sldpth[kMaxSoilLevels];
  for (int k = 0; k < nsoil; ++k) sldpth[k] = real(sl[k]);

  std::vector<double> og(size_t(ncols) * (kNg + nsoil)), os(size_t(ncols) * (kNs + nsoil));
  auto v = [&](int i, int n) { return real(vin[size_t(n) * ncols + i]); };
  for (int i = 0; i < ncols; ++i) {
    const real ffrozp = v(i, 0), lwdn = v(i, 1), solnet = v(i, 2), sfcprs = v(i, 3), prcp = v(i, 4),
               sfctmp = v(i, 5), q2 = v(i, 6), th2 = v(i, 7), q2sat = v(i, 8), dqsdt2 = v(i, 9), alb = v(i, 10),
               snoalb = v(i, 11), tbot = v(i, 12), embrd = v(i, 15), ch = v(i, 19), sithick = v(i, 23);
    real stc[kMaxSoilLevels];
    {
      real z0brd = v(i, 13), emissi = v(i, 14), t1 = v(i, 16), snowh = v(i, 17), sneqv = v(i, 18),
           sncovr = v(i, 20), snotime1 = v(i, 21), ribb = v(i, 22);
      for (int k = 0; k < nsoil; ++k) stc[k] = real(stcin[size_t(k) * ncols + i]);
      glacial::Out o{};
      glacial::sflx_glacial(ffrozp, dt, nsoil, sldpth, lwdn, solnet, sfcprs, prcp, sfctmp, q2, th2, q2sat, dqsdt2,
                            alb, snoalb, tbot, z0brd, emissi, embrd, t1, stc, snowh, sneqv, ch, sncovr, snotime1,
                            ribb, lvcoef, o);
      const real out[kNg] = {t1,     snowh,  sneqv,  sncovr, snotime1, ribb,      emissi,  o.albedo, o.z0,
                             o.eta,  o.sheat, o.eta_kinematic, o.fdown, o.esnow, o.dew, o.etp,  o.ssoil,
                             o.flx1, o.flx2, o.flx3, o.snomlt, o.runoff1, o.q1, ch, alb, z0brd, embrd};
      for (int n = 0; n < kNg; ++n) og[size_t(n) * ncols + i] = double(out[n]);
      for (int k = 0; k < nsoil; ++k) og[size_t(kNg + k) * ncols + i] = double(stc[k]);
    }
    {
      real z0brd = v(i, 13), emissi = v(i, 14), t1 = v(i, 16), snowh = v(i, 17), sneqv = v(i, 18),
           sncovr = v(i, 20), ribb = v(i, 22);
      for (int k = 0; k < nsoil; ++k) stc[k] = real(stcin[size_t(k) * ncols + i]);
      seaice::Opts op{};
      seaice::Out o{};
      if (!seaice::sflx_seaice(op, ffrozp, dt, nsoil, sithick, lwdn, solnet, sfcprs, prcp, sfctmp, q2, th2,
                               q2sat, dqsdt2, snoalb, tbot, z0brd, emissi, t1, stc, snowh, sneqv, ch, sncovr, ribb,
                               o)) {
        std::fprintf(stderr, "column %d: snow density > 1 (WRF stops)\n", i + 1);
        return 1;
      }
      const real out[kNs] = {t1,      snowh,  sneqv,  sncovr,          ribb,    emissi,  o.albedo, o.z0,
                             z0brd,   o.eta,  o.sheat, o.eta_kinematic, o.fdown, o.esnow, o.dew,   o.etp,
                             o.ssoil, o.flx1, o.flx2, o.flx3,          o.snomlt, o.runoff1, o.q1, ch};
      for (int n = 0; n < kNs; ++n) os[size_t(n) * ncols + i] = double(out[n]);
      for (int k = 0; k < nsoil; ++k) os[size_t(kNs + k) * ncols + i] = double(stc[k]);
    }
  }

  std::FILE* g = std::fopen(argv[2], "wb");
  if (!g) return 1;
  std::fwrite(hdr, sizeof(int32_t), 2, g);
  std::fwrite(og.data(), sizeof(double), og.size(), g);
  std::fwrite(os.data(), sizeof(double), os.size(), g);
  std::fclose(g);
  return 0;
}
