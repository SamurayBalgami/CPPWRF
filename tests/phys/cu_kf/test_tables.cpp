#include "src/phys/cu_kf/kf_tables.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace cppwrf::phys::cu_kf;

int main() {
    KFTables<double> tables_cpp{};
    init_kf_tables(tables_cpp);

    std::ifstream in("tables_fortran.bin", std::ios::binary);
    if (!in) {
        std::cerr << "Failed to open tables_fortran.bin\n";
        return 1;
    }

    double rdpr_f = 0, rdthk_f = 0, plutop_f = 0;
    in.read(reinterpret_cast<char*>(&rdpr_f), sizeof(double));
    in.read(reinterpret_cast<char*>(&rdthk_f), sizeof(double));
    in.read(reinterpret_cast<char*>(&plutop_f), sizeof(double));

    double the0k_f[KFNP];
    in.read(reinterpret_cast<char*>(the0k_f), sizeof(the0k_f));

    double alu_f[KFNALU];
    in.read(reinterpret_cast<char*>(alu_f), sizeof(alu_f));

    // Fortran ttab(250, 220) stream read: 250*220 doubles
    double ttab_f[KFNP][KFNT];
    in.read(reinterpret_cast<char*>(ttab_f), sizeof(ttab_f));

    double qstab_f[KFNP][KFNT];
    in.read(reinterpret_cast<char*>(qstab_f), sizeof(qstab_f));

    std::cout << std::scientific << std::setprecision(10);

    auto check_diff = [](const char* name, double c, double f, double tol = 1e-12) {
        double diff = std::abs(c - f);
        double denom = std::max(1.0, std::max(std::abs(c), std::abs(f)));
        double rel = diff / denom;
        if (rel > tol) {
            std::cerr << "FAIL " << name << ": C++=" << c << " Fortran=" << f << " rel=" << rel << "\n";
            return false;
        }
        return true;
    };

    bool pass = true;
    pass &= check_diff("rdpr", tables_cpp.rdpr, rdpr_f);
    pass &= check_diff("rdthk", tables_cpp.rdthk, rdthk_f);
    pass &= check_diff("plutop", tables_cpp.plutop, plutop_f);

    double max_rel_the0k = 0;
    for (int i = 0; i < KFNP; ++i) {
        double diff = std::abs(tables_cpp.the0k[i] - the0k_f[i]);
        double rel = diff / the0k_f[i];
        if (rel > max_rel_the0k) max_rel_the0k = rel;
        pass &= check_diff("the0k", tables_cpp.the0k[i], the0k_f[i]);
    }

    double max_rel_alu = 0;
    for (int i = 0; i < KFNALU; ++i) {
        double diff = std::abs(tables_cpp.alu[i] - alu_f[i]);
        double rel = diff / std::max(1.0, std::abs(alu_f[i]));
        if (rel > max_rel_alu) max_rel_alu = rel;
        pass &= check_diff("alu", tables_cpp.alu[i], alu_f[i]);
    }

    double max_rel_ttab = 0;
    for (int kp = 0; kp < KFNP; ++kp) {
        for (int it = 0; it < KFNT; ++it) {
            double diff = std::abs(tables_cpp.ttab[kp][it] - ttab_f[kp][it]);
            double rel = diff / ttab_f[kp][it];
            if (rel > max_rel_ttab) max_rel_ttab = rel;
            if (rel > 1e-12) {
                std::cerr << "FAIL ttab[" << kp << "][" << it << "]: C=" << tables_cpp.ttab[kp][it] << " F=" << ttab_f[kp][it] << " rel=" << rel << "\n";
                pass = false;
                break;
            }
        }
    }

    double max_rel_qstab = 0;
    for (int kp = 0; kp < KFNP; ++kp) {
        for (int it = 0; it < KFNT; ++it) {
            double diff = std::abs(tables_cpp.qstab[kp][it] - qstab_f[kp][it]);
            double denom = std::max(1e-10, std::abs(qstab_f[kp][it]));
            double rel = diff / denom;
            if (rel > max_rel_qstab) max_rel_qstab = rel;
            if (rel > 1e-12) {
                std::cerr << "FAIL qstab[" << kp << "][" << it << "]: C=" << tables_cpp.qstab[kp][it] << " F=" << qstab_f[kp][it] << " rel=" << rel << "\n";
                pass = false;
                break;
            }
        }
    }

    std::cout << "Max rel diff the0k: " << max_rel_the0k << "\n";
    std::cout << "Max rel diff alu:   " << max_rel_alu << "\n";
    std::cout << "Max rel diff ttab:  " << max_rel_ttab << "\n";
    std::cout << "Max rel diff qstab: " << max_rel_qstab << "\n";

    if (pass) {
        std::cout << "All Kain-Fritsch tables verified successfully! PASS <= 1e-12\n";
        return 0;
    } else {
        std::cerr << "Table verification FAILED!\n";
        return 1;
    }
}
