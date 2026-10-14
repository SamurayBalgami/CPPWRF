#include "phys/mp_thompson/thompson_tables.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>

using namespace cppwrf::phys::mp_thompson;

struct TableStat {
    std::string name;
    size_t count;
    double max_abs_diff;
    double max_rel_diff;
    bool pass;
};

template <typename T>
bool compare_buffer(const std::string& name, const T* cpp_data, const double* f_data, size_t count, double tol, TableStat& stat) {
    double max_abs = 0.0;
    double max_rel = 0.0;
    bool pass = true;

    for (size_t i = 0; i < count; ++i) {
        double c = static_cast<double>(cpp_data[i]);
        double f = f_data[i];
        double diff = std::abs(c - f);
        if (diff > max_abs) max_abs = diff;

        double denom = std::max(1.0e-12, std::abs(f));
        double rel = diff / denom;
        if (rel > max_rel) max_rel = rel;

        if (rel > tol) {
            pass = false;
        }
    }

    stat = {name, count, max_abs, max_rel, pass};
    return pass;
}

int main(int argc, char** argv) {
    std::string table_path = "tables_fortran.bin";
    if (argc > 1) {
        table_path = argv[1];
    }

    std::ifstream in(table_path, std::ios::binary);
    if (!in) {
        std::cerr << "Failed to open " << table_path << "\n";
        return 1;
    }

    std::cout << "Initializing C++ Thompson tables...\n";
    ThompsonHostTables<double> tables_cpp{};
    init_thompson_tables(tables_cpp);
    std::cout << "C++ Thompson tables initialized.\n\n";

    std::vector<TableStat> stats;
    bool overall_pass = true;
    constexpr double TOL = 1.0e-12;

    auto test_one = [&](const std::string& name, const double* cpp_data, size_t count) {
        std::vector<double> f_data(count);
        in.read(reinterpret_cast<char*>(f_data.data()), count * sizeof(double));
        if (!in) {
            std::cerr << "Error reading " << name << " from " << table_path << "\n";
            overall_pass = false;
            return;
        }
        TableStat st;
        bool ok = compare_buffer(name, cpp_data, f_data.data(), count, TOL, st);
        stats.push_back(st);
        if (!ok) overall_pass = false;
    };

    // 1. Graupel-rain collection (5 tables)
    test_one("tcg_racg", tables_cpp.tcg_racg.data(), N_RACG);
    test_one("tmr_racg", tables_cpp.tmr_racg.data(), N_RACG);
    test_one("tcr_gacr", tables_cpp.tcr_gacr.data(), N_RACG);
    test_one("tnr_racg", tables_cpp.tnr_racg.data(), N_RACG);
    test_one("tnr_gacr", tables_cpp.tnr_gacr.data(), N_RACG);

    // 2. Snow-rain collection (12 tables)
    test_one("tcs_racs1", tables_cpp.tcs_racs1.data(), N_RACS);
    test_one("tmr_racs1", tables_cpp.tmr_racs1.data(), N_RACS);
    test_one("tcs_racs2", tables_cpp.tcs_racs2.data(), N_RACS);
    test_one("tmr_racs2", tables_cpp.tmr_racs2.data(), N_RACS);
    test_one("tcr_sacr1", tables_cpp.tcr_sacr1.data(), N_RACS);
    test_one("tms_sacr1", tables_cpp.tms_sacr1.data(), N_RACS);
    test_one("tcr_sacr2", tables_cpp.tcr_sacr2.data(), N_RACS);
    test_one("tms_sacr2", tables_cpp.tms_sacr2.data(), N_RACS);
    test_one("tnr_racs1", tables_cpp.tnr_racs1.data(), N_RACS);
    test_one("tnr_racs2", tables_cpp.tnr_racs2.data(), N_RACS);
    test_one("tnr_sacr1", tables_cpp.tnr_sacr1.data(), N_RACS);
    test_one("tnr_sacr2", tables_cpp.tnr_sacr2.data(), N_RACS);

    // 3. Rain freezing (4 tables)
    test_one("tpi_qrfz", tables_cpp.tpi_qrfz.data(), N_QRFZ);
    test_one("tpg_qrfz", tables_cpp.tpg_qrfz.data(), N_QRFZ);
    test_one("tni_qrfz", tables_cpp.tni_qrfz.data(), N_QRFZ);
    test_one("tnr_qrfz", tables_cpp.tnr_qrfz.data(), N_QRFZ);

    // 4. Cloud water freezing (2 tables)
    test_one("tpi_qcfz", tables_cpp.tpi_qcfz.data(), N_QCFZ);
    test_one("tni_qcfz", tables_cpp.tni_qcfz.data(), N_QCFZ);

    // 5. Ice autoconversion to snow (3 tables)
    test_one("tps_iaus", tables_cpp.tps_iaus.data(), N_IAUS);
    test_one("tni_iaus", tables_cpp.tni_iaus.data(), N_IAUS);
    test_one("tpi_ide",  tables_cpp.tpi_ide.data(),  N_IAUS);

    // 6. Collision efficiencies (2 tables)
    test_one("t_Efrw", tables_cpp.t_Efrw.data(), N_EFRW);
    test_one("t_Efsw", tables_cpp.t_Efsw.data(), N_EFSW);

    // 7. Droplet evaporation (2 tables)
    test_one("tpc_wev", tables_cpp.tpc_wev.data(), N_WEV);
    test_one("tnc_wev", tables_cpp.tnc_wev.data(), N_WEV);

    std::cout << std::left << std::setw(15) << "Table Name"
              << std::right << std::setw(12) << "Elements"
              << std::setw(18) << "Max Abs Diff"
              << std::setw(18) << "Max Rel Diff"
              << std::setw(10) << "Status" << "\n";
    std::cout << std::string(73, '-') << "\n";

    std::cout << std::scientific << std::setprecision(6);
    for (const auto& s : stats) {
        std::cout << std::left << std::setw(15) << s.name
                  << std::right << std::setw(12) << s.count
                  << std::setw(18) << s.max_abs_diff
                  << std::setw(18) << s.max_rel_diff
                  << std::setw(10) << (s.pass ? "PASS" : "FAIL") << "\n";
    }
    std::cout << std::string(73, '-') << "\n";

    if (overall_pass) {
        std::cout << "All Thompson lookup tables verified against Fortran! (PASS <= 1e-12)\n";
        return 0;
    } else {
        std::cerr << "Table verification FAILED for one or more tables!\n";
        return 1;
    }
}
