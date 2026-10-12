#include "init_real.hpp"
#include <netcdf.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>

using namespace cppwrf::init_real;
namespace fs = std::filesystem;

#define NC_CHECK(cmd) do { \
    int _e = (cmd); \
    if (_e != NC_NOERR) { \
        std::cerr << "NetCDF error: " << nc_strerror(_e) << " at line " << __LINE__ << std::endl; \
        return false; \
    } \
} while(0)

struct DiffStats {
    std::string varname;
    double max_abs = 0.0;
    double max_rel = 0.0;
    double rms = 0.0;
    size_t count = 0;
    std::string worst_loc;
    double worst_ref = 0.0;
    double worst_cpp = 0.0;
};

static DiffStats compare_buffers(const std::string& name,
                                 const float* ref,
                                 const float* cpp,
                                 size_t total_size,
                                 const std::vector<int>& dims) {
    DiffStats stat;
    stat.varname = name;
    stat.count = total_size;
    double sum_sq = 0.0;

    for (size_t idx = 0; idx < total_size; ++idx) {
        double r = ref[idx];
        double c = cpp[idx];
        double diff = std::abs(c - r);
        double rel = (std::abs(r) > 1.0e-7) ? (diff / std::abs(r)) : diff;

        sum_sq += diff * diff;

        if (diff > stat.max_abs) {
            stat.max_abs = diff;
            stat.max_rel = rel;
            stat.worst_ref = r;
            stat.worst_cpp = c;

            // Compute indices
            if (dims.size() == 1) {
                stat.worst_loc = "[" + std::to_string(idx) + "]";
            } else if (dims.size() == 2) {
                int j = idx / dims[1];
                int i = idx % dims[1];
                stat.worst_loc = "[" + std::to_string(j) + "," + std::to_string(i) + "]";
            } else if (dims.size() == 3) {
                int k = idx / (dims[1] * dims[2]);
                int rem = idx % (dims[1] * dims[2]);
                int j = rem / dims[2];
                int i = rem % dims[2];
                stat.worst_loc = "[" + std::to_string(k) + "," + std::to_string(j) + "," + std::to_string(i) + "]";
            }
        }
    }

    if (total_size > 0) {
        stat.rms = std::sqrt(sum_sq / static_cast<double>(total_size));
    }
    return stat;
}

static bool read_nc_var(int ncid, const std::string& name, std::vector<float>& buf, std::vector<size_t>& shape) {
    int varid;
    if (nc_inq_varid(ncid, name.c_str(), &varid) != NC_NOERR) {
        return false;
    }
    int ndims;
    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
    std::vector<int> dimids(ndims);
    NC_CHECK(nc_inq_vardimid(ncid, varid, dimids.data()));

    shape.resize(ndims);
    size_t total = 1;
    for (int d = 0; d < ndims; ++d) {
        NC_CHECK(nc_inq_dimlen(ncid, dimids[d], &shape[d]));
        total *= shape[d];
    }
    buf.resize(total);
    NC_CHECK(nc_get_var_float(ncid, varid, buf.data()));
    return true;
}

static void print_table_header() {
    std::cout << "| Variable | Max Abs Diff | Max Rel Diff | RMS Diff | Worst Point | Ref Value | C++ Value |\n";
    std::cout << "|:---------|:-------------|:-------------|:---------|:------------|:----------|:----------|\n";
}

static void print_stat(const DiffStats& s) {
    std::cout << "| " << std::left << std::setw(8) << s.varname
              << " | " << std::scientific << std::setprecision(3) << s.max_abs
              << " | " << std::scientific << std::setprecision(3) << s.max_rel
              << " | " << std::scientific << std::setprecision(3) << s.rms
              << " | " << std::left << std::setw(11) << s.worst_loc
              << " | " << std::scientific << std::setprecision(4) << s.worst_ref
              << " | " << std::scientific << std::setprecision(4) << s.worst_cpp
              << " |\n";
}

int main(int argc, char* argv[]) {
    std::string run_dir = "/home/ostriquetrum/wrf_real/wrf_run";
    if (argc > 1) run_dir = argv[1];

    std::vector<std::string> met_files = {
        run_dir + "/met_em.d01.2026-10-12_00:00:00.nc",
        run_dir + "/met_em.d01.2026-10-12_06:00:00.nc",
        run_dir + "/met_em.d01.2026-10-12_12:00:00.nc",
        run_dir + "/met_em.d01.2026-10-12_18:00:00.nc",
        run_dir + "/met_em.d01.2026-10-13_00:00:00.nc"
    };

    for (const auto& f : met_files) {
        if (!fs::exists(f)) {
            std::cout << "Reference met_em file missing: " << f << ". Skipping test.\n";
            return 0; // CTest skip
        }
    }

    std::cout << "=================================================================\n";
    std::cout << " running differential tests against WRF 4.6.1 reference\n";
    std::cout << "=================================================================\n\n";

    // Test 1: hybrid_opt = 2
    std::string ref_input_h2 = run_dir + "/wrfinput_d01_hybrid2";
    if (fs::exists(ref_input_h2)) {
        std::cout << "### Test 1: wrfinput_d01 Comparison (hybrid_opt = 2)\n\n";
        auto res = run_real_init(met_files, 2);
        const auto& s = res.initial_state;

        int ncid;
        if (nc_open(ref_input_h2.c_str(), NC_NOWRITE, &ncid) == NC_NOERR) {
            print_table_header();
            auto compare_1d = [&](const std::string& name, const Array1D& arr) {
                std::vector<float> ref;
                std::vector<size_t> shp;
                if (read_nc_var(ncid, name, ref, shp)) {
                    print_stat(compare_buffers(name, ref.data(), arr.data_ptr(), arr.size(), {static_cast<int>(arr.size())}));
                }
            };
            auto compare_2d = [&](const std::string& name, const Array2D& arr) {
                std::vector<float> ref;
                std::vector<size_t> shp;
                if (read_nc_var(ncid, name, ref, shp)) {
                    print_stat(compare_buffers(name, ref.data(), arr.data_ptr(), arr.size(), {arr.d1, arr.d2}));
                }
            };
            auto compare_3d = [&](const std::string& name, const Array3D& arr) {
                std::vector<float> ref;
                std::vector<size_t> shp;
                if (read_nc_var(ncid, name, ref, shp)) {
                    print_stat(compare_buffers(name, ref.data(), arr.data_ptr(), arr.size(), {arr.d1, arr.d2, arr.d3}));
                }
            };

            compare_1d("ZNW", s.znw);
            compare_1d("ZNU", s.znu);
            compare_1d("DNW", s.dnw);
            compare_1d("DN", s.dn);
            compare_1d("RDNW", s.rdnw);
            compare_1d("RDN", s.rdn);
            compare_1d("C1F", s.c1f);
            compare_1d("C2F", s.c2f);
            compare_1d("C3F", s.c3f);
            compare_1d("C4F", s.c4f);
            compare_1d("C1H", s.c1h);
            compare_1d("C2H", s.c2h);
            compare_1d("C3H", s.c3h);
            compare_1d("C4H", s.c4h);

            compare_2d("MUB", s.mub);
            compare_2d("MU", s.mu);
            compare_2d("MAPFAC_M", s.msft);
            compare_2d("MAPFAC_U", s.msfu);
            compare_2d("MAPFAC_V", s.msfv);
            compare_2d("F", s.f);
            compare_2d("E", s.e);
            compare_2d("SINALPHA", s.sina);
            compare_2d("COSALPHA", s.cosa);
            compare_2d("HGT", s.hgt);

            compare_3d("U", s.u);
            compare_3d("V", s.v);
            compare_3d("W", s.w);
            compare_3d("PH", s.ph);
            compare_3d("PHB", s.phb);
            compare_3d("T", s.t);
            compare_3d("T_INIT", s.t_init);
            compare_3d("P", s.p);
            compare_3d("PB", s.pb);
            compare_3d("AL", s.al);
            compare_3d("ALB", s.alb);
            compare_3d("QVAPOR", s.qvapor);

            nc_close(ncid);
            std::cout << "\n";
        }
    }

    // Test 2: hybrid_opt = 0
    std::string ref_input_h0 = run_dir + "/wrfinput_d01_hybrid0";
    if (fs::exists(ref_input_h0)) {
        std::cout << "### Test 2: wrfinput_d01 Comparison (hybrid_opt = 0)\n\n";
        auto res = run_real_init({met_files[0]}, 0);
        const auto& s = res.initial_state;

        int ncid;
        if (nc_open(ref_input_h0.c_str(), NC_NOWRITE, &ncid) == NC_NOERR) {
            print_table_header();
            auto compare_1d = [&](const std::string& name, const Array1D& arr) {
                std::vector<float> ref;
                std::vector<size_t> shp;
                if (read_nc_var(ncid, name, ref, shp)) {
                    print_stat(compare_buffers(name, ref.data(), arr.data_ptr(), arr.size(), {static_cast<int>(arr.size())}));
                }
            };
            auto compare_2d = [&](const std::string& name, const Array2D& arr) {
                std::vector<float> ref;
                std::vector<size_t> shp;
                if (read_nc_var(ncid, name, ref, shp)) {
                    print_stat(compare_buffers(name, ref.data(), arr.data_ptr(), arr.size(), {arr.d1, arr.d2}));
                }
            };
            auto compare_3d = [&](const std::string& name, const Array3D& arr) {
                std::vector<float> ref;
                std::vector<size_t> shp;
                if (read_nc_var(ncid, name, ref, shp)) {
                    print_stat(compare_buffers(name, ref.data(), arr.data_ptr(), arr.size(), {arr.d1, arr.d2, arr.d3}));
                }
            };

            compare_1d("C1F", s.c1f);
            compare_1d("C2F", s.c2f);
            compare_1d("C3F", s.c3f);
            compare_1d("C4F", s.c4f);
            compare_1d("C1H", s.c1h);
            compare_1d("C2H", s.c2h);
            compare_1d("C3H", s.c3h);
            compare_1d("C4H", s.c4h);

            compare_2d("MUB", s.mub);
            compare_2d("MU", s.mu);

            compare_3d("U", s.u);
            compare_3d("V", s.v);
            compare_3d("T", s.t);
            compare_3d("P", s.p);
            compare_3d("PB", s.pb);
            compare_3d("PH", s.ph);
            compare_3d("PHB", s.phb);

            nc_close(ncid);
            std::cout << "\n";
        }
    }

    // Test 3: wrfbdy_d01
    std::string ref_bdy = run_dir + "/wrfbdy_d01";
    if (fs::exists(ref_bdy)) {
        std::cout << "### Test 3: wrfbdy_d01 Lateral Boundary Comparison (hybrid_opt = 2)\n\n";
        auto res = run_real_init(met_files, 2);
        const auto& bdy = res.bdy_data;

        int ncid;
        if (nc_open(ref_bdy.c_str(), NC_NOWRITE, &ncid) == NC_NOERR) {
            print_table_header();
            auto compare_bdy_3d = [&](const std::string& name, const Array3D& arr, int interval_idx) {
                int varid;
                if (nc_inq_varid(ncid, name.c_str(), &varid) == NC_NOERR) {
                    size_t start[4] = {static_cast<size_t>(interval_idx), 0, 0, 0};
                    size_t count[4] = {1, static_cast<size_t>(arr.d1), static_cast<size_t>(arr.d2), static_cast<size_t>(arr.d3)};
                    std::vector<float> ref(arr.size());
                    if (nc_get_vara_float(ncid, varid, start, count, ref.data()) == NC_NOERR) {
                        print_stat(compare_buffers(name + "[t0]", ref.data(), arr.data_ptr(), arr.size(), {arr.d1, arr.d2, arr.d3}));
                    }
                }
            };
            auto compare_bdy_2d = [&](const std::string& name, const Array3D& arr, int interval_idx) {
                int varid;
                if (nc_inq_varid(ncid, name.c_str(), &varid) == NC_NOERR) {
                    size_t start[3] = {static_cast<size_t>(interval_idx), 0, 0};
                    size_t count[3] = {1, static_cast<size_t>(arr.d1), static_cast<size_t>(arr.d3)};
                    std::vector<float> ref(arr.d1 * arr.d3);
                    if (nc_get_vara_float(ncid, varid, start, count, ref.data()) == NC_NOERR) {
                        print_stat(compare_buffers(name + "[t0]", ref.data(), arr.data_ptr(), ref.size(), {arr.d1, arr.d3}));
                    }
                }
            };

            if (!bdy.intervals.empty()) {
                const auto& iv = bdy.intervals[0];
                compare_bdy_3d("U_BXS", iv.u_bxs, 0);
                compare_bdy_3d("U_BXE", iv.u_bxe, 0);
                compare_bdy_3d("U_BYS", iv.u_bys, 0);
                compare_bdy_3d("U_BYE", iv.u_bye, 0);
                compare_bdy_3d("U_BTXS", iv.u_btxs, 0);
                compare_bdy_3d("U_BTXE", iv.u_btxe, 0);

                compare_bdy_3d("V_BXS", iv.v_bxs, 0);
                compare_bdy_3d("V_BXE", iv.v_bxe, 0);
                compare_bdy_3d("V_BYS", iv.v_bys, 0);
                compare_bdy_3d("V_BYE", iv.v_bye, 0);
                compare_bdy_3d("V_BTXS", iv.v_btxs, 0);
                compare_bdy_3d("V_BTXE", iv.v_btxe, 0);

                compare_bdy_3d("T_BXS", iv.t_bxs, 0);
                compare_bdy_3d("T_BXE", iv.t_bxe, 0);
                compare_bdy_3d("T_BTXS", iv.t_btxs, 0);
                compare_bdy_3d("T_BTXE", iv.t_btxe, 0);

                compare_bdy_3d("PH_BXS", iv.ph_bxs, 0);
                compare_bdy_3d("PH_BXE", iv.ph_bxe, 0);
                compare_bdy_3d("PH_BTXS", iv.ph_btxs, 0);
                compare_bdy_3d("PH_BTXE", iv.ph_btxe, 0);

                compare_bdy_2d("MU_BXS", iv.mu_bxs, 0);
                compare_bdy_2d("MU_BXE", iv.mu_bxe, 0);
                compare_bdy_2d("MU_BYS", iv.mu_bys, 0);
                compare_bdy_2d("MU_BYE", iv.mu_bye, 0);
                compare_bdy_2d("MU_BTXS", iv.mu_btxs, 0);
                compare_bdy_2d("MU_BTXE", iv.mu_btxe, 0);
            }
            nc_close(ncid);
            std::cout << "\n";
        }
    }

    return 0;
}
