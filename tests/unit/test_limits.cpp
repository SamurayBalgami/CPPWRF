// SPDX-License-Identifier: GPL-3.0-or-later
// The column-size guards of the physics runners are errors in every build type (Release: NDEBUG): a column above
// kMaxLevels (or more soil / snow layers than the kernels hold) throws before any kernel is launched, with the
// scheme, the size and the limit in the message (src/core/limits.hpp).
#include <cstdio>
#include <stdexcept>
#include <string>

#include "core/backend.hpp"
#include "core/precision.hpp"
#include "phys/bl_ysu/ysu.hpp"
#include "phys/cu_kf/kf.hpp"
#include "phys/lsm_noah/noah.hpp"
#include "phys/lsm_noahmp/noahmp.hpp"
#include "phys/mp_thompson/thompson.hpp"

using namespace cppwrf;
static int fails = 0;

// f() must throw std::runtime_error whose message contains want
template <class F> static void expect(const char* what, const std::string& want, F&& f) {
  try {
    f();
  } catch (const std::runtime_error& e) {
    if (std::string(e.what()).find(want) != std::string::npos) {
      std::printf("ok   %s: %s\n", what, e.what());
      return;
    }
    std::printf("FAIL %s: message '%s' lacks '%s'\n", what, e.what(), want.c_str());
    ++fails;
    return;
  }
  std::printf("FAIL %s: no error\n", what);
  ++fails;
}

int main() {
  const Range2D cols{0, 1, 0, 1};
  {
    namespace y = phys::ysu;
    const int nz = y::kMaxLevels + 1;
    expect("YSU", "YSU: nz = " + std::to_string(nz) + " exceeds the column limit " + std::to_string(y::kMaxLevels),
           [&] { y::run_ysu({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, real(60), cols, nz); });
  }
  {
    namespace kf = phys::cu_kf;
    kf::KFKernel<real> k{};
    k.nz = kf::kMaxLevels + 1;
    expect("Kain-Fritsch", "Kain-Fritsch: nz = " + std::to_string(k.nz) + " exceeds the column limit",
           [&] { kf::run_kf(k, cols); });
  }
  {
    namespace th = phys::mp_thompson;
    const int nz = th::kMaxLevels + 1;
    const th::ThompsonTablesView<real> tables{};
    expect("Thompson", "Thompson: nz = " + std::to_string(nz) + " exceeds the column limit",
           [&] { th::run_thompson<real>({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, real(60), cols, nz, tables); });
  }
  {
    namespace no = phys::noah;
    no::NoahKernel k{};
    k.nsoil = no::kMaxSoilLevels + 1;
    expect("Noah", "Noah: num_soil_layers = " + std::to_string(k.nsoil) + " exceeds the column limit",
           [&] { no::run_noah(k, cols); });
  }
  {
    namespace mp = phys::noahmp;
    mp::NoahmpKernel k{};
    k.nsoil = mp::kMaxSoilLayers + 1;
    expect("Noah-MP soil", "Noah-MP: num_soil_layers", [&] { mp::run_noahmp(k, cols); });
    k.nsoil = mp::kMaxSoilLayers;
    k.nsnow = mp::kMaxSnowLayers + 1;
    expect("Noah-MP snow", "Noah-MP: snow layers", [&] { mp::run_noahmp(k, cols); });
    mp::NoahmpInitKernel ki{};
    ki.nsoil = mp::kMaxSoilLayers + 1;
    expect("Noah-MP init", "Noah-MP: num_soil_layers", [&] { mp::run_noahmp_init(ki, cols); });
  }
  std::printf("%s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
  return fails ? 1 : 0;
}
