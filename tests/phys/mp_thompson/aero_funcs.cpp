// The C++ aerosol functions of the Thompson kernel (activ_ncloud, eff_aero, ice_demott, ice_koop) on the samples
// of aero_funcs.py, against fortran/aero_funcs.f90 (WRF's own functions); CCN_ACTIVATE.BIN read from ".".
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "phys/mp_thompson/thompson.hpp"
#include "phys/mp_thompson/thompson_tables.hpp"

using namespace cppwrf;
namespace mt = cppwrf::phys::mp_thompson;

int main(int argc, char** argv) {
    if (argc < 3) return 2;
    FILE* f = std::fopen(argv[1], "rb");
    if (!f) return 1;
    int32_t n = 0;
    if (std::fread(&n, 4, 1, f) != 1) return 1;
    std::vector<std::vector<double>> in(17, std::vector<double>(std::size_t(n)));
    for (auto& v : in)
        if (std::fread(v.data(), 8, std::size_t(n), f) != std::size_t(n)) return 1;
    std::fclose(f);
    mt::ThompsonHostTables<double> h;
    mt::init_thompson_tables(h, ".", false, true);
    const mt::ThompsonTablesView<double> tb = h.view();
    const char species[3] = {'r', 's', 'g'};
    std::vector<double> out(std::size_t(n) * 4);
    for (int i = 0; i < n; ++i) {
        auto a = [&](int c) { return in[std::size_t(c)][std::size_t(i)]; };
        out[std::size_t(i)] = mt::activ_ncloud(a(0), a(1), a(2), tb);
        out[std::size_t(n) + i] = mt::eff_aero(a(3), a(4), a(5), a(6), a(7), species[int(std::lround(a(8)))]);
        out[2 * std::size_t(n) + i] = mt::ice_demott(a(9), a(10), a(11));
        out[3 * std::size_t(n) + i] = mt::ice_koop(a(12), a(13), a(14), a(15), a(16));
    }
    FILE* o = std::fopen(argv[2], "wb");
    if (!o) return 1;
    std::fwrite(out.data(), 8, out.size(), o);
    std::fclose(o);
    return 0;
}
