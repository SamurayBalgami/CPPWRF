#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/cu_kf/kf.hpp"
#include "phys/cu_kf/kf_tables.hpp"

using namespace cppwrf;

template <typename T>
struct DevBuffer {
    T* p = nullptr;
    std::size_t n = 0;
    DevBuffer() = default;
    explicit DevBuffer(std::size_t count)
        : n(count), p(static_cast<T*>(dev_alloc(count * sizeof(T)))) {}
    ~DevBuffer() {
        if (p) dev_free(p);
    }
    DevBuffer(DevBuffer&& o) noexcept : p(o.p), n(o.n) {
        o.p = nullptr;
        o.n = 0;
    }
    DevBuffer& operator=(DevBuffer&& o) noexcept {
        if (this != &o) {
            if (p) dev_free(p);
            p = o.p;
            n = o.n;
            o.p = nullptr;
            o.n = 0;
        }
        return *this;
    }
    DevBuffer(const DevBuffer&) = delete;
    DevBuffer& operator=(const DevBuffer&) = delete;
    T* get() const { return p; }
};

template <typename D, typename S>
void h2d(D* dst, const S* src, std::size_t count) {
    if constexpr (std::is_same_v<D, S>) {
        copy_to_dev(dst, src, count * sizeof(D));
    } else {
        std::vector<D> tmp(count);
        for (std::size_t i = 0; i < count; ++i) tmp[i] = static_cast<D>(src[i]);
        copy_to_dev(dst, tmp.data(), count * sizeof(D));
    }
}

template <typename D, typename S>
void d2h(D* dst, const S* src, std::size_t count) {
    if constexpr (std::is_same_v<D, S>) {
        copy_to_host(dst, src, count * sizeof(D));
    } else {
        std::vector<S> tmp(count);
        copy_to_host(tmp.data(), src, count * sizeof(S));
        for (std::size_t i = 0; i < count; ++i) dst[i] = static_cast<D>(tmp[i]);
    }
}

int main(int argc, char** argv) {
    std::string infile = "input_kf.bin";
    std::string outfile = "output_cpp.bin";
    bool bench_mode = false;
    bool bench_100k = false;
    int bench_iters = 10;

    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg == "--bench") {
            bench_mode = true;
            if (a + 1 < argc && argv[a + 1][0] != '-') {
                bench_iters = std::atoi(argv[++a]);
            }
        } else if (arg == "--bench-100k") {
            bench_mode = true;
            bench_100k = true;
            if (a + 1 < argc && argv[a + 1][0] != '-') {
                bench_iters = std::atoi(argv[++a]);
            }
        } else if (a == 1 && arg[0] != '-') {
            infile = arg;
        } else if (a == 2 && arg[0] != '-') {
            outfile = arg;
        }
    }

    FILE* f = std::fopen(infile.c_str(), "rb");
    if (!f) {
        std::fprintf(stderr, "Cannot open input file: %s\n", infile.c_str());
        return 1;
    }

    int32_t ncols = 0, nlevels = 0;
    double dt_in = 0, dx_in = 0, cudt_in = 0;
    int32_t stepcu = 0, ktau = 0, trigger = 0, warm_rain_int = 0;
    int32_t flag_qi_int = 0, flag_qs_int = 0, flag_qr_int = 0;

    if (std::fread(&ncols, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&nlevels, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&dt_in, sizeof(double), 1, f) != 1 ||
        std::fread(&dx_in, sizeof(double), 1, f) != 1 ||
        std::fread(&cudt_in, sizeof(double), 1, f) != 1 ||
        std::fread(&stepcu, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&ktau, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&trigger, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&warm_rain_int, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&flag_qi_int, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&flag_qs_int, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&flag_qr_int, sizeof(int32_t), 1, f) != 1) {
        std::fprintf(stderr, "Failed reading header from %s\n", infile.c_str());
        std::fclose(f);
        return 1;
    }

    std::size_t n2d = std::size_t(ncols);
    std::size_t n3d = std::size_t(ncols) * std::size_t(nlevels);
    std::size_t n3d_w = std::size_t(ncols) * std::size_t(nlevels + 1);

    // Host buffers for 3D inputs
    std::vector<double> h_u(n3d), h_v(n3d), h_th(n3d), h_t(n3d), h_qv(n3d);
    std::vector<double> h_dz8w(n3d), h_pcps(n3d), h_rho(n3d), h_pi(n3d), h_w0avg(n3d);
    std::vector<double> h_w(n3d_w);

    // Read 3D inputs
    auto read_3d = [&](std::vector<double>& v, std::size_t count) {
        if (std::fread(v.data(), sizeof(double), count, f) != count) {
            std::fprintf(stderr, "Failed reading 3D array from %s\n", infile.c_str());
            return false;
        }
        return true;
    };

    if (!read_3d(h_u, n3d) || !read_3d(h_v, n3d) || !read_3d(h_th, n3d) ||
        !read_3d(h_t, n3d) || !read_3d(h_qv, n3d) || !read_3d(h_dz8w, n3d) ||
        !read_3d(h_pcps, n3d) || !read_3d(h_rho, n3d) || !read_3d(h_pi, n3d) ||
        !read_3d(h_w0avg, n3d) || !read_3d(h_w, n3d_w)) {
        std::fclose(f);
        return 1;
    }

    // Host buffers for 2D inputs
    std::vector<double> h_raincv(n2d), h_pratec(n2d), h_nca(n2d);
    std::vector<double> h_cutop(n2d), h_cubot(n2d);

    auto read_2d = [&](std::vector<double>& v) {
        if (std::fread(v.data(), sizeof(double), n2d, f) != n2d) {
            std::fprintf(stderr, "Failed reading 2D array from %s\n", infile.c_str());
            return false;
        }
        return true;
    };

    if (!read_2d(h_raincv) || !read_2d(h_pratec) || !read_2d(h_nca) ||
        !read_2d(h_cutop) || !read_2d(h_cubot)) {
        std::fclose(f);
        return 1;
    }
    std::fclose(f);
 
    if (bench_100k) {
        int n_orig = ncols;
        int n_target = 100000;
        auto expand_3d = [&](std::vector<double>& v, int lvls) {
            std::vector<double> expanded(std::size_t(n_target) * std::size_t(lvls));
            for (int k = 0; k < lvls; ++k) {
                for (int c = 0; c < n_target; ++c) {
                    expanded[k * n_target + c] = v[k * n_orig + (c % n_orig)];
                }
            }
            v = std::move(expanded);
        };
        auto expand_2d = [&](std::vector<double>& v) {
            std::vector<double> expanded(n_target);
            for (int c = 0; c < n_target; ++c) {
                expanded[c] = v[c % n_orig];
            }
            v = std::move(expanded);
        };
        expand_3d(h_u, nlevels);
        expand_3d(h_v, nlevels);
        expand_3d(h_th, nlevels);
        expand_3d(h_t, nlevels);
        expand_3d(h_qv, nlevels);
        expand_3d(h_dz8w, nlevels);
        expand_3d(h_pcps, nlevels);
        expand_3d(h_rho, nlevels);
        expand_3d(h_pi, nlevels);
        expand_3d(h_w0avg, nlevels);
        expand_3d(h_w, nlevels + 1);

        expand_2d(h_raincv);
        expand_2d(h_pratec);
        expand_2d(h_nca);
        expand_2d(h_cutop);
        expand_2d(h_cubot);

        ncols = n_target;
        n2d = std::size_t(ncols);
        n3d = std::size_t(ncols) * std::size_t(nlevels);
        n3d_w = std::size_t(ncols) * std::size_t(nlevels + 1);
    }

    // Initialize lookup tables
    phys::cu_kf::KFTables<real> h_tables;
    phys::cu_kf::init_kf_tables(h_tables);

    DevBuffer<phys::cu_kf::KFTables<real>> d_tables(1);
    h2d(d_tables.get(), &h_tables, 1);

    // Allocate device buffers
    DevBuffer<real> d_u(n3d), d_v(n3d), d_th(n3d), d_t(n3d), d_qv(n3d);
    DevBuffer<real> d_dz8w(n3d), d_pcps(n3d), d_rho(n3d), d_pi(n3d), d_w0avg(n3d);
    DevBuffer<real> d_w(n3d_w);

    DevBuffer<real> d_raincv(n2d), d_pratec(n2d), d_nca(n2d);
    DevBuffer<real> d_cutop(n2d), d_cubot(n2d);

    DevBuffer<real> d_rthcuten(n3d), d_rqvcuten(n3d), d_rqccuten(n3d);
    DevBuffer<real> d_rqrcuten(n3d), d_rqicuten(n3d), d_rqscuten(n3d);
    DevBuffer<real> d_cldfra_dp_kf(n3d), d_cldfra_sh_kf(n3d), d_qc_kf(n3d), d_qi_kf(n3d);

    auto copy_to_dev = [](const std::vector<double>& h, DevBuffer<real>& d, std::size_t count) {
        h2d(d.get(), h.data(), count);
    };

    copy_to_dev(h_u, d_u, n3d);
    copy_to_dev(h_v, d_v, n3d);
    copy_to_dev(h_th, d_th, n3d);
    copy_to_dev(h_t, d_t, n3d);
    copy_to_dev(h_qv, d_qv, n3d);
    copy_to_dev(h_dz8w, d_dz8w, n3d);
    copy_to_dev(h_pcps, d_pcps, n3d);
    copy_to_dev(h_rho, d_rho, n3d);
    copy_to_dev(h_pi, d_pi, n3d);
    copy_to_dev(h_w0avg, d_w0avg, n3d);
    copy_to_dev(h_w, d_w, n3d_w);

    copy_to_dev(h_raincv, d_raincv, n2d);
    copy_to_dev(h_pratec, d_pratec, n2d);
    copy_to_dev(h_nca, d_nca, n2d);
    copy_to_dev(h_cutop, d_cutop, n2d);
    copy_to_dev(h_cubot, d_cubot, n2d);

    // Setup 3D and 2D views
    int sj3d = ncols * nlevels;
    int sk3d = ncols;
    int sj3d_w = ncols * (nlevels + 1);
    int sk3d_w = ncols;
    int sj2d = ncols;

    auto v2c = [&](const DevBuffer<real>& b) { return View2<const real>{b.get(), sj2d}; };
    auto v2 = [&](DevBuffer<real>& b) { return View2<real>{b.get(), sj2d}; };
    auto v3c = [&](const DevBuffer<real>& b) { return View3<const real>{b.get(), sj3d, sk3d}; };
    auto v3 = [&](DevBuffer<real>& b) { return View3<real>{b.get(), sj3d, sk3d}; };
    auto v3c_w = [&](const DevBuffer<real>& b) { return View3<const real>{b.get(), sj3d_w, sk3d_w}; };

    phys::cu_kf::KFKernel<real> kernel;
    kernel.nx = ncols;
    kernel.ny = 1;
    kernel.nz = nlevels;
    kernel.dt = static_cast<real>(dt_in);
    kernel.dx = static_cast<real>(dx_in);
    kernel.cudt = static_cast<real>(cudt_in);
    kernel.stepcu = stepcu;
    kernel.ktau = ktau;
    kernel.trigger = trigger;
    kernel.warm_rain = (warm_rain_int != 0);
    kernel.flag_qi = (flag_qi_int != 0);
    kernel.flag_qs = (flag_qs_int != 0);
    kernel.flag_qr = (flag_qr_int != 0);

    kernel.xlv0 = static_cast<real>(3.15e6);
    kernel.xlv1 = static_cast<real>(2370.0);
    kernel.xls0 = static_cast<real>(2.905e6);
    kernel.xls1 = static_cast<real>(259.532);
    kernel.cp = static_cast<real>(1004.5);
    kernel.r = static_cast<real>(287.0);
    kernel.g = static_cast<real>(9.81);
    kernel.ep1 = static_cast<real>(461.6 / 287.0 - 1.0);
    kernel.ep2 = static_cast<real>(287.0 / 461.6);
    kernel.svp1 = static_cast<real>(0.6112);
    kernel.svp2 = static_cast<real>(17.67);
    kernel.svp3 = static_cast<real>(29.65);
    kernel.svpt0 = static_cast<real>(273.15);

    kernel.tables = d_tables.get();

    kernel.u = v3c(d_u);
    kernel.v = v3c(d_v);
    kernel.w = v3c_w(d_w);
    kernel.th = v3c(d_th);
    kernel.t = v3c(d_t);
    kernel.qv = v3c(d_qv);
    kernel.dz8w = v3c(d_dz8w);
    kernel.pcps = v3c(d_pcps);
    kernel.rho = v3c(d_rho);
    kernel.pi = v3c(d_pi);

    kernel.w0avg = v3(d_w0avg);
    kernel.raincv = v2(d_raincv);
    kernel.pratec = v2(d_pratec);
    kernel.nca = v2(d_nca);
    kernel.cutop = v2(d_cutop);
    kernel.cubot = v2(d_cubot);

    kernel.rthcuten = v3(d_rthcuten);
    kernel.rqvcuten = v3(d_rqvcuten);
    kernel.rqccuten = v3(d_rqccuten);
    kernel.rqrcuten = v3(d_rqrcuten);
    kernel.rqicuten = v3(d_rqicuten);
    kernel.rqscuten = v3(d_rqscuten);

    kernel.cldfra_dp_kf = v3(d_cldfra_dp_kf);
    kernel.cldfra_sh_kf = v3(d_cldfra_sh_kf);
    kernel.qc_kf = v3(d_qc_kf);
    kernel.qi_kf = v3(d_qi_kf);

    Range2D range{0, ncols, 0, 1};

    // Warm-up run
    phys::cu_kf::run_kf(kernel, range);
    dev_sync();

    if (bench_mode) {
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < bench_iters; ++iter) {
            phys::cu_kf::run_kf(kernel, range);
        }
        dev_sync();
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed_sec = std::chrono::duration<double>(t1 - t0).count();
        double ms_per_iter = (elapsed_sec / bench_iters) * 1000.0;
        double ns_per_col = (elapsed_sec / (bench_iters * double(ncols))) * 1.0e9;

        std::cout << "Benchmark results (" << bench_iters << " iterations, "
                  << ncols << " columns x " << nlevels << " levels):\n";
        std::cout << "  Total time:    " << elapsed_sec * 1000.0 << " ms\n";
        std::cout << "  Time/iter:     " << ms_per_iter << " ms\n";
        std::cout << "  ns/column:     " << ns_per_col << " ns/col\n";
    }

    if (bench_100k) {
        return 0;
    }

    // Copy outputs back to host
    std::vector<double> out_raincv(n2d), out_pratec(n2d), out_nca(n2d);
    std::vector<double> out_cutop(n2d), out_cubot(n2d);
    std::vector<double> out_w0avg(n3d), out_rthcuten(n3d), out_rqvcuten(n3d);
    std::vector<double> out_rqccuten(n3d), out_rqrcuten(n3d), out_rqicuten(n3d), out_rqscuten(n3d);
    std::vector<double> out_cldfra_dp(n3d), out_cldfra_sh(n3d), out_qc(n3d), out_qi(n3d);

    auto copy_from_dev = [](const DevBuffer<real>& d, std::vector<double>& h, std::size_t count) {
        d2h(h.data(), d.get(), count);
    };

    copy_from_dev(d_raincv, out_raincv, n2d);
    copy_from_dev(d_pratec, out_pratec, n2d);
    copy_from_dev(d_nca, out_nca, n2d);
    copy_from_dev(d_cutop, out_cutop, n2d);
    copy_from_dev(d_cubot, out_cubot, n2d);

    copy_from_dev(d_w0avg, out_w0avg, n3d);
    copy_from_dev(d_rthcuten, out_rthcuten, n3d);
    copy_from_dev(d_rqvcuten, out_rqvcuten, n3d);
    copy_from_dev(d_rqccuten, out_rqccuten, n3d);
    copy_from_dev(d_rqrcuten, out_rqrcuten, n3d);
    copy_from_dev(d_rqicuten, out_rqicuten, n3d);
    copy_from_dev(d_rqscuten, out_rqscuten, n3d);
    copy_from_dev(d_cldfra_dp_kf, out_cldfra_dp, n3d);
    copy_from_dev(d_cldfra_sh_kf, out_cldfra_sh, n3d);
    copy_from_dev(d_qc_kf, out_qc, n3d);
    copy_from_dev(d_qi_kf, out_qi, n3d);

    FILE* fo = std::fopen(outfile.c_str(), "wb");
    if (!fo) {
        std::fprintf(stderr, "Cannot open output file: %s\n", outfile.c_str());
        return 1;
    }

    auto write_buf = [&](const std::vector<double>& v) {
        std::fwrite(v.data(), sizeof(double), v.size(), fo);
    };

    // 2D fields
    write_buf(out_raincv);
    write_buf(out_pratec);
    write_buf(out_nca);
    write_buf(out_cutop);
    write_buf(out_cubot);

    // 3D fields
    write_buf(out_w0avg);
    write_buf(out_rthcuten);
    write_buf(out_rqvcuten);
    write_buf(out_rqccuten);
    write_buf(out_rqrcuten);
    write_buf(out_rqicuten);
    write_buf(out_rqscuten);
    write_buf(out_cldfra_dp);
    write_buf(out_cldfra_sh);
    write_buf(out_qc);
    write_buf(out_qi);

    std::fclose(fo);
    std::cout << "Wrote " << outfile << " successfully for " << ncols << " columns.\n";
    return 0;
}
