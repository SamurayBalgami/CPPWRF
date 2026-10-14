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
#include "phys/mp_thompson/thompson.hpp"
#include "phys/mp_thompson/thompson_tables.hpp"

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
    std::string infile = "input_thompson.bin";
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
    double dt_in = 0.0;

    if (std::fread(&ncols, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&nlevels, sizeof(int32_t), 1, f) != 1 ||
        std::fread(&dt_in, sizeof(double), 1, f) != 1) {
        std::fprintf(stderr, "Failed reading header from %s\n", infile.c_str());
        std::fclose(f);
        return 1;
    }

    std::size_t n2d = std::size_t(ncols);
    std::size_t n3d = std::size_t(ncols) * std::size_t(nlevels);

    std::vector<double> h_th(n3d), h_qv(n3d), h_qc(n3d), h_qr(n3d);
    std::vector<double> h_qi(n3d), h_qs(n3d), h_qg(n3d), h_ni(n3d), h_nr(n3d);
    std::vector<double> h_pii(n3d), h_p(n3d), h_w(n3d), h_dz(n3d);

    auto read_3d = [&](std::vector<double>& v) {
        if (std::fread(v.data(), sizeof(double), n3d, f) != n3d) {
            std::fprintf(stderr, "Failed reading 3D array from %s\n", infile.c_str());
            return false;
        }
        return true;
    };

    if (!read_3d(h_th) || !read_3d(h_qv) || !read_3d(h_qc) || !read_3d(h_qr) ||
        !read_3d(h_qi) || !read_3d(h_qs) || !read_3d(h_qg) || !read_3d(h_ni) || !read_3d(h_nr) ||
        !read_3d(h_pii) || !read_3d(h_p) || !read_3d(h_w) || !read_3d(h_dz)) {
        std::fclose(f);
        return 1;
    }

    std::vector<double> h_rainnc(n2d), h_rainncv(n2d), h_snownc(n2d), h_snowncv(n2d);
    std::vector<double> h_graupelnc(n2d), h_graupelncv(n2d), h_sr(n2d);

    auto read_2d = [&](std::vector<double>& v) {
        if (std::fread(v.data(), sizeof(double), n2d, f) != n2d) {
            std::fprintf(stderr, "Failed reading 2D array from %s\n", infile.c_str());
            return false;
        }
        return true;
    };

    if (!read_2d(h_rainnc) || !read_2d(h_rainncv) || !read_2d(h_snownc) || !read_2d(h_snowncv) ||
        !read_2d(h_graupelnc) || !read_2d(h_graupelncv) || !read_2d(h_sr)) {
        std::fclose(f);
        return 1;
    }
    std::fclose(f);

    if (bench_100k) {
        int n_orig = ncols;
        int n_target = 100000;
        auto expand_3d = [&](std::vector<double>& v) {
            std::vector<double> expanded(std::size_t(n_target) * std::size_t(nlevels));
            for (int k = 0; k < nlevels; ++k) {
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

        expand_3d(h_th); expand_3d(h_qv); expand_3d(h_qc); expand_3d(h_qr);
        expand_3d(h_qi); expand_3d(h_qs); expand_3d(h_qg); expand_3d(h_ni); expand_3d(h_nr);
        expand_3d(h_pii); expand_3d(h_p); expand_3d(h_w); expand_3d(h_dz);

        expand_2d(h_rainnc); expand_2d(h_rainncv);
        expand_2d(h_snownc); expand_2d(h_snowncv);
        expand_2d(h_graupelnc); expand_2d(h_graupelncv);
        expand_2d(h_sr);

        ncols = n_target;
        n2d = std::size_t(ncols);
        n3d = std::size_t(ncols) * std::size_t(nlevels);
    }

    // Initialize/load tables
    phys::mp_thompson::ThompsonHostTables<real> h_tables;
    if (!phys::mp_thompson::load_thompson_tables(h_tables, "tables_fortran.bin")) {
        std::cout << "Generating Thompson tables in C++...\n";
        phys::mp_thompson::init_thompson_tables(h_tables);
    }

    phys::mp_thompson::ThompsonDeviceTables<real> dev_tables;
    dev_tables.upload(h_tables);
    phys::mp_thompson::ThompsonTablesView<real> d_tables_view = dev_tables.view();

    // Device buffers for 3D state
    DevBuffer<real> d_th(n3d), d_qv(n3d), d_qc(n3d), d_qr(n3d);
    DevBuffer<real> d_qi(n3d), d_qs(n3d), d_qg(n3d), d_ni(n3d), d_nr(n3d);
    DevBuffer<real> d_pii(n3d), d_p(n3d), d_w(n3d), d_dz(n3d);

    // Device buffers for 2D state
    DevBuffer<real> d_rainnc(n2d), d_rainncv(n2d), d_snownc(n2d), d_snowncv(n2d);
    DevBuffer<real> d_graupelnc(n2d), d_graupelncv(n2d), d_sr(n2d);

    auto copy_to_d = [](const std::vector<double>& h, DevBuffer<real>& d, std::size_t count) {
        h2d(d.get(), h.data(), count);
    };

    copy_to_d(h_th, d_th, n3d);
    copy_to_d(h_qv, d_qv, n3d);
    copy_to_d(h_qc, d_qc, n3d);
    copy_to_d(h_qr, d_qr, n3d);
    copy_to_d(h_qi, d_qi, n3d);
    copy_to_d(h_qs, d_qs, n3d);
    copy_to_d(h_qg, d_qg, n3d);
    copy_to_d(h_ni, d_ni, n3d);
    copy_to_d(h_nr, d_nr, n3d);
    copy_to_d(h_pii, d_pii, n3d);
    copy_to_d(h_p, d_p, n3d);
    copy_to_d(h_w, d_w, n3d);
    copy_to_d(h_dz, d_dz, n3d);

    copy_to_d(h_rainnc, d_rainnc, n2d);
    copy_to_d(h_rainncv, d_rainncv, n2d);
    copy_to_d(h_snownc, d_snownc, n2d);
    copy_to_d(h_snowncv, d_snowncv, n2d);
    copy_to_d(h_graupelnc, d_graupelnc, n2d);
    copy_to_d(h_graupelncv, d_graupelncv, n2d);
    copy_to_d(h_sr, d_sr, n2d);

    int sj3d = ncols * nlevels;
    int sk3d = ncols;
    int sj2d = ncols;

    auto v3 = [&](DevBuffer<real>& b) { return View3<real>{b.get(), sj3d, sk3d}; };
    auto v3c = [&](const DevBuffer<real>& b) { return View3<const real>{b.get(), sj3d, sk3d}; };
    auto v2 = [&](DevBuffer<real>& b) { return View2<real>{b.get(), sj2d}; };

    Range2D range{0, ncols, 0, 1};

    auto launch_kernel = [&]() {
        phys::mp_thompson::run_thompson<real>(
            v3(d_th), v3(d_qv), v3(d_qc), v3(d_qr),
            v3(d_qi), v3(d_qs), v3(d_qg), v3(d_ni), v3(d_nr),
            v3c(d_pii), v3c(d_p), v3c(d_w), v3c(d_dz),
            v2(d_rainnc), v2(d_rainncv), v2(d_snownc), v2(d_snowncv),
            v2(d_graupelnc), v2(d_graupelncv), v2(d_sr),
            static_cast<real>(dt_in), range, nlevels, d_tables_view
        );
    };

    // Warm-up run
    launch_kernel();
    dev_sync();

    if (bench_mode) {
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < bench_iters; ++iter) {
            launch_kernel();
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
    std::vector<double> out_qv(n3d), out_qc(n3d), out_qr(n3d), out_qi(n3d);
    std::vector<double> out_qs(n3d), out_qg(n3d), out_ni(n3d), out_nr(n3d), out_th(n3d);
    std::vector<double> out_rainnc(n2d), out_rainncv(n2d), out_snownc(n2d), out_snowncv(n2d);
    std::vector<double> out_graupelnc(n2d), out_graupelncv(n2d), out_sr(n2d);

    auto copy_from_d = [](const DevBuffer<real>& d, std::vector<double>& h, std::size_t count) {
        d2h(h.data(), d.get(), count);
    };

    copy_from_d(d_qv, out_qv, n3d);
    copy_from_d(d_qc, out_qc, n3d);
    copy_from_d(d_qr, out_qr, n3d);
    copy_from_d(d_qi, out_qi, n3d);
    copy_from_d(d_qs, out_qs, n3d);
    copy_from_d(d_qg, out_qg, n3d);
    copy_from_d(d_ni, out_ni, n3d);
    copy_from_d(d_nr, out_nr, n3d);
    copy_from_d(d_th, out_th, n3d);

    copy_from_d(d_rainnc, out_rainnc, n2d);
    copy_from_d(d_rainncv, out_rainncv, n2d);
    copy_from_d(d_snownc, out_snownc, n2d);
    copy_from_d(d_snowncv, out_snowncv, n2d);
    copy_from_d(d_graupelnc, out_graupelnc, n2d);
    copy_from_d(d_graupelncv, out_graupelncv, n2d);
    copy_from_d(d_sr, out_sr, n2d);

    FILE* fo = std::fopen(outfile.c_str(), "wb");
    if (!fo) {
        std::fprintf(stderr, "Cannot open output file: %s\n", outfile.c_str());
        return 1;
    }

    auto write_buf = [&](const std::vector<double>& v) {
        std::fwrite(v.data(), sizeof(double), v.size(), fo);
    };

    // 3D outputs: qv, qc, qr, qi, qs, qg, ni, nr, th
    write_buf(out_qv);
    write_buf(out_qc);
    write_buf(out_qr);
    write_buf(out_qi);
    write_buf(out_qs);
    write_buf(out_qg);
    write_buf(out_ni);
    write_buf(out_nr);
    write_buf(out_th);

    // 2D outputs: rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, sr
    write_buf(out_rainnc);
    write_buf(out_rainncv);
    write_buf(out_snownc);
    write_buf(out_snowncv);
    write_buf(out_graupelnc);
    write_buf(out_graupelncv);
    write_buf(out_sr);

    std::fclose(fo);
    std::cout << "Wrote " << outfile << " successfully for " << ncols << " columns.\n";
    return 0;
}
