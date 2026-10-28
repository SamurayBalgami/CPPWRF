#pragma once

#include <vector>
#include <cmath>
#include <cstdint>
#include <string>
#include "core/backend.hpp"

namespace cppwrf::phys::mp_thompson {

// Table dimensions
constexpr int nbins = 100;
constexpr int nbc = nbins;
constexpr int nbi = nbins;
constexpr int nbr = nbins;
constexpr int nbs = nbins;
constexpr int nbg = nbins;

constexpr int ntb_c = 37;
constexpr int ntb_i = 64;
constexpr int ntb_r = 37;
constexpr int ntb_s = 37;
constexpr int ntb_g = 37;
constexpr int ntb_g1 = 37;
constexpr int ntb_r1 = 37;
constexpr int ntb_i1 = 55;
constexpr int ntb_t = 9;
constexpr int ntb_IN = 55;
constexpr int dimNRHG1 = 1;
constexpr int ntb_freeze_temp = 45;

// Total elements per table
constexpr size_t N_RACG = static_cast<size_t>(ntb_g1) * ntb_g * dimNRHG1 * ntb_r1 * ntb_r; // 1,874,161
constexpr size_t N_RACS = static_cast<size_t>(ntb_s) * ntb_t * ntb_r1 * ntb_r;             // 455,949
constexpr size_t N_QRFZ = static_cast<size_t>(ntb_r) * ntb_r1 * ntb_freeze_temp * ntb_IN;  // 3,388,275
constexpr size_t N_QCFZ = static_cast<size_t>(ntb_c) * nbc * ntb_freeze_temp * ntb_IN;    // 9,157,500
constexpr size_t N_IAUS = static_cast<size_t>(ntb_i) * ntb_i1;                             // 3,520
constexpr size_t N_EFRW = static_cast<size_t>(nbr) * nbc;                                  // 10,000
constexpr size_t N_EFSW = static_cast<size_t>(nbs) * nbc;                                  // 10,000
constexpr size_t N_WEV  = static_cast<size_t>(nbc) * ntb_c * nbc;                          // 370,000

// Trivially copyable view holding pointers and scalar constants for CPU/CUDA kernels
template <typename Real>
struct ThompsonTablesView {
    // 1. Graupel-rain collection (5 tables, 37x37x1x37x37)
    // Fortran indices: (i=ntb_g1, j=ntb_g, n3=dimNRHG, k=ntb_r1, m=ntb_r)
    const Real* tcg_racg{nullptr};
    const Real* tmr_racg{nullptr};
    const Real* tcr_gacr{nullptr};
    const Real* tnr_racg{nullptr};
    const Real* tnr_gacr{nullptr};

    // 2. Snow-rain collection (12 tables, 37x9x37x37)
    // Fortran indices: (i=ntb_s, j=ntb_t, k=ntb_r1, m=ntb_r)
    const Real* tcs_racs1{nullptr};
    const Real* tmr_racs1{nullptr};
    const Real* tcs_racs2{nullptr};
    const Real* tmr_racs2{nullptr};
    const Real* tcr_sacr1{nullptr};
    const Real* tms_sacr1{nullptr};
    const Real* tcr_sacr2{nullptr};
    const Real* tms_sacr2{nullptr};
    const Real* tnr_racs1{nullptr};
    const Real* tnr_racs2{nullptr};
    const Real* tnr_sacr1{nullptr};
    const Real* tnr_sacr2{nullptr};

    // 3. Rain freezing (4 tables, 37x37x45x55)
    // Fortran indices: (i=ntb_r, j=ntb_r1, k=45, m=ntb_IN)
    const Real* tpi_qrfz{nullptr};
    const Real* tpg_qrfz{nullptr};
    const Real* tni_qrfz{nullptr};
    const Real* tnr_qrfz{nullptr};

    // 4. Cloud water freezing (2 tables, 37x100x45x55)
    // Fortran indices: (i=ntb_c, j=nbc, k=45, m=ntb_IN)
    const Real* tpi_qcfz{nullptr};
    const Real* tni_qcfz{nullptr};

    // 5. Ice autoconversion to snow (3 tables, 64x55)
    // Fortran indices: (i=ntb_i, j=ntb_i1)
    const Real* tps_iaus{nullptr};
    const Real* tni_iaus{nullptr};
    const Real* tpi_ide{nullptr};

    // 6. Collision efficiencies (2 tables, 100x100)
    // Fortran indices: (i=nbr/nbs, j=nbc)
    const Real* t_Efrw{nullptr};
    const Real* t_Efsw{nullptr};

    // 7. Cloud droplet evaporation (2 tables, 100x37x100)
    // Fortran indices: (i=nbc, j=ntb_c, k=nbc)
    const Real* tpc_wev{nullptr};
    const Real* tnc_wev{nullptr};

    // Precomputed physical constants
    Real t1_qr_qc;
    Real t1_qr_qi;
    Real t2_qr_qi;
    Real t1_qs_qc;
    Real t1_qs_qi;
    Real t1_qr_ev;
    Real t2_qr_ev;
    Real t1_qs_sd;
    Real t2_qs_sd;
    Real t1_qs_me;
    Real t2_qs_me;
    Real t1_qg_sd;
    Real t2_qg_sd;
    Real t1_qg_me;
    Real t2_qg_me;

    // Indexing constants
    int nic1, nic2, nii2, nii3, nir2, nir3, nis2, nig2, nig3, niIN2;

    // Table lower-bound coordinates (index 0)
    Real r_c[1];
    Real r_i[1];
    Real r_r[1];
    Real r_g[1];
    Real r_s[1];
    Real Nt_i[1];
    Real Nt_IN[1];

    // Gamma and exponent parameters
    Real cce[5][15];
    Real ccg[5][15];
    Real ocg1[15];
    Real ocg2[15];

    Real cie[7];
    Real cig[7];
    Real oig1;
    Real oig2;
    Real obmi;

    Real cre[13];
    Real crg[13];
    Real ore1;
    Real org1;
    Real org2;
    Real org3;
    Real obmr;

    Real cse[17];
    Real csg[17];
    Real oams;
    Real obms;
    Real ocms;

    Real cge[12];
    Real cgg[12];
    Real oamg;
    Real obmg;
    Real ocmg;
    Real oge1;
    Real ogg1;
    Real ogg2;
    Real ogg3;

    Real sa[10];
    Real sb[10];

    Real Dr_0;
    Real log_Dr;
    Real inv_log_Dr;
    Real Ds_0;
    Real log_Ds;
    Real inv_log_Ds;
    Real t_Nc_0;
    Real inv_nic1;

    Real t1_qg_qc;
    Real D0i;
    Real xm0s;
    Real xm0g;

    CPPWRF_HD inline int get_table_index(Real val, int ni_ref, int ntb_max) const {
        int ni_val = static_cast<int>(std::round(std::log10(val)));
        int n = ni_val;
        for (int nn = ni_val - 1; nn <= ni_val + 1; ++nn) {
            n = nn;
            Real p10 = powi(Real(10.0), nn);  // WRF: 10.**nn, INTEGER nn (__powidf2)
            Real scaled = val / p10;
            if (scaled >= Real(1.0) && scaled < Real(10.0)) {
                break;
            }
        }
        Real p10 = powi(Real(10.0), n);
        int idx = static_cast<int>(val / p10) + 10 * (n - ni_ref) - (n - ni_ref);
        if (idx < 1) idx = 1;
        if (idx > ntb_max) idx = ntb_max;
        return idx;
    }

    // Helper functions for table indexing matching Fortran 1:1
    CPPWRF_HD inline size_t idx_racg(int i, int j, int n3, int k, int m) const {
        return static_cast<size_t>(i) + static_cast<size_t>(ntb_g1) * (
               static_cast<size_t>(j) + static_cast<size_t>(ntb_g) * (
               static_cast<size_t>(n3) + static_cast<size_t>(dimNRHG1) * (
               static_cast<size_t>(k) + static_cast<size_t>(ntb_r1) * static_cast<size_t>(m))));
    }

    CPPWRF_HD inline size_t idx_racs(int i, int j, int k, int m) const {
        return static_cast<size_t>(i) + static_cast<size_t>(ntb_s) * (
               static_cast<size_t>(j) + static_cast<size_t>(ntb_t) * (
               static_cast<size_t>(k) + static_cast<size_t>(ntb_r1) * static_cast<size_t>(m)));
    }

    CPPWRF_HD inline size_t idx_qrfz(int i, int j, int k, int m) const {
        return static_cast<size_t>(i) + static_cast<size_t>(ntb_r) * (
               static_cast<size_t>(j) + static_cast<size_t>(ntb_r1) * (
               static_cast<size_t>(k) + static_cast<size_t>(ntb_freeze_temp) * static_cast<size_t>(m)));
    }

    CPPWRF_HD inline size_t idx_qcfz(int i, int j, int k, int m) const {
        return static_cast<size_t>(i) + static_cast<size_t>(ntb_c) * (
               static_cast<size_t>(j) + static_cast<size_t>(nbc) * (
               static_cast<size_t>(k) + static_cast<size_t>(ntb_freeze_temp) * static_cast<size_t>(m)));
    }

    CPPWRF_HD inline size_t idx_iaus(int i, int j) const {
        return static_cast<size_t>(i) + static_cast<size_t>(ntb_i) * static_cast<size_t>(j);
    }

    CPPWRF_HD inline size_t idx_efrw(int i, int j) const {
        return static_cast<size_t>(i) + static_cast<size_t>(nbr) * static_cast<size_t>(j);
    }

    CPPWRF_HD inline size_t idx_efsw(int i, int j) const {
        return static_cast<size_t>(i) + static_cast<size_t>(nbs) * static_cast<size_t>(j);
    }

    CPPWRF_HD inline size_t idx_wev(int i, int j, int k) const {
        return static_cast<size_t>(i) + static_cast<size_t>(nbc) * (
               static_cast<size_t>(j) + static_cast<size_t>(ntb_c) * static_cast<size_t>(k));
    }
};

// Host table container owning allocated memory
template <typename Real>
struct ThompsonHostTables {
    // 1. Graupel-rain collection (5 tables)
    std::vector<Real> tcg_racg{std::vector<Real>(N_RACG, Real(0))};
    std::vector<Real> tmr_racg{std::vector<Real>(N_RACG, Real(0))};
    std::vector<Real> tcr_gacr{std::vector<Real>(N_RACG, Real(0))};
    std::vector<Real> tnr_racg{std::vector<Real>(N_RACG, Real(0))};
    std::vector<Real> tnr_gacr{std::vector<Real>(N_RACG, Real(0))};

    // 2. Snow-rain collection (12 tables)
    std::vector<Real> tcs_racs1{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tmr_racs1{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tcs_racs2{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tmr_racs2{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tcr_sacr1{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tms_sacr1{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tcr_sacr2{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tms_sacr2{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tnr_racs1{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tnr_racs2{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tnr_sacr1{std::vector<Real>(N_RACS, Real(0))};
    std::vector<Real> tnr_sacr2{std::vector<Real>(N_RACS, Real(0))};

    // 3. Rain freezing (4 tables)
    std::vector<Real> tpi_qrfz{std::vector<Real>(N_QRFZ, Real(0))};
    std::vector<Real> tpg_qrfz{std::vector<Real>(N_QRFZ, Real(0))};
    std::vector<Real> tni_qrfz{std::vector<Real>(N_QRFZ, Real(0))};
    std::vector<Real> tnr_qrfz{std::vector<Real>(N_QRFZ, Real(0))};

    // 4. Cloud water freezing (2 tables)
    std::vector<Real> tpi_qcfz{std::vector<Real>(N_QCFZ, Real(0))};
    std::vector<Real> tni_qcfz{std::vector<Real>(N_QCFZ, Real(0))};

    // 5. Ice autoconversion to snow (3 tables)
    std::vector<Real> tps_iaus{std::vector<Real>(N_IAUS, Real(0))};
    std::vector<Real> tni_iaus{std::vector<Real>(N_IAUS, Real(0))};
    std::vector<Real> tpi_ide{std::vector<Real>(N_IAUS, Real(0))};

    // 6. Collision efficiencies (2 tables)
    std::vector<Real> t_Efrw{std::vector<Real>(N_EFRW, Real(0))};
    std::vector<Real> t_Efsw{std::vector<Real>(N_EFSW, Real(0))};

    // 7. Cloud droplet evaporation (2 tables)
    std::vector<Real> tpc_wev{std::vector<Real>(N_WEV, Real(0))};
    std::vector<Real> tnc_wev{std::vector<Real>(N_WEV, Real(0))};

    // Precomputed constants
    Real t1_qr_qc{0};
    Real t1_qr_qi{0};
    Real t2_qr_qi{0};
    Real t1_qs_qc{0};
    Real t1_qs_qi{0};
    Real t1_qr_ev{0};
    Real t2_qr_ev{0};
    Real t1_qs_sd{0};
    Real t2_qs_sd{0};
    Real t1_qs_me{0};
    Real t2_qs_me{0};
    Real t1_qg_sd{0};
    Real t2_qg_sd{0};
    Real t1_qg_me{0};
    Real t2_qg_me{0};

    int nic1{0}, nic2{0}, nii2{0}, nii3{0}, nir2{0}, nir3{0}, nis2{0}, nig2{0}, nig3{0}, niIN2{0};

    Real r_c[ntb_c];
    Real r_i[ntb_i];
    Real r_r[ntb_r];
    Real r_g[ntb_g];
    Real r_s[ntb_s];
    Real N0r_exp[ntb_r1];
    Real N0g_exp[ntb_g1];
    Real Nt_i[ntb_i1];
    Real Nt_IN[ntb_IN];
    Real Tc[ntb_t];

    // Gamma and exponent parameters
    Real cce[5][15];
    Real ccg[5][15];
    Real ocg1[15];
    Real ocg2[15];

    Real cie[7];
    Real cig[7];
    Real oig1{0};
    Real oig2{0};
    Real obmi{0};

    Real cre[13];
    Real crg[13];
    Real ore1{0};
    Real org1{0};
    Real org2{0};
    Real org3{0};
    Real obmr{0};

    Real cse[17];
    Real csg[17];
    Real oams{0};
    Real obms{0};
    Real ocms{0};

    Real cge[12];
    Real cgg[12];
    Real oamg{0};
    Real obmg{0};
    Real ocmg{0};
    Real oge1{0};
    Real ogg1{0};
    Real ogg2{0};
    Real ogg3{0};

    Real sa[10];
    Real sb[10];

    Real Dr_0{0};
    Real log_Dr{0};
    Real inv_log_Dr{0};
    Real Ds_0{0};
    Real log_Ds{0};
    Real inv_log_Ds{0};
    Real t_Nc_0{0};
    Real inv_nic1{0};

    Real t1_qg_qc{0};
    Real D0i{0};
    Real xm0s{0};
    Real xm0g{0};

    // Export view
    ThompsonTablesView<Real> view() const {
        ThompsonTablesView<Real> v{};
        v.tcg_racg = tcg_racg.data();
        v.tmr_racg = tmr_racg.data();
        v.tcr_gacr = tcr_gacr.data();
        v.tnr_racg = tnr_racg.data();
        v.tnr_gacr = tnr_gacr.data();

        v.tcs_racs1 = tcs_racs1.data();
        v.tmr_racs1 = tmr_racs1.data();
        v.tcs_racs2 = tcs_racs2.data();
        v.tmr_racs2 = tmr_racs2.data();
        v.tcr_sacr1 = tcr_sacr1.data();
        v.tms_sacr1 = tms_sacr1.data();
        v.tcr_sacr2 = tcr_sacr2.data();
        v.tms_sacr2 = tms_sacr2.data();
        v.tnr_racs1 = tnr_racs1.data();
        v.tnr_racs2 = tnr_racs2.data();
        v.tnr_sacr1 = tnr_sacr1.data();
        v.tnr_sacr2 = tnr_sacr2.data();

        v.tpi_qrfz = tpi_qrfz.data();
        v.tpg_qrfz = tpg_qrfz.data();
        v.tni_qrfz = tni_qrfz.data();
        v.tnr_qrfz = tnr_qrfz.data();

        v.tpi_qcfz = tpi_qcfz.data();
        v.tni_qcfz = tni_qcfz.data();

        v.tps_iaus = tps_iaus.data();
        v.tni_iaus = tni_iaus.data();
        v.tpi_ide = tpi_ide.data();

        v.t_Efrw = t_Efrw.data();
        v.t_Efsw = t_Efsw.data();

        v.tpc_wev = tpc_wev.data();
        v.tnc_wev = tnc_wev.data();

        v.t1_qr_qc = t1_qr_qc;
        v.t1_qr_qi = t1_qr_qi;
        v.t2_qr_qi = t2_qr_qi;
        v.t1_qs_qc = t1_qs_qc;
        v.t1_qs_qi = t1_qs_qi;
        v.t1_qr_ev = t1_qr_ev;
        v.t2_qr_ev = t2_qr_ev;
        v.t1_qs_sd = t1_qs_sd;
        v.t2_qs_sd = t2_qs_sd;
        v.t1_qs_me = t1_qs_me;
        v.t2_qs_me = t2_qs_me;
        v.t1_qg_sd = t1_qg_sd;
        v.t2_qg_sd = t2_qg_sd;
        v.t1_qg_me = t1_qg_me;
        v.t2_qg_me = t2_qg_me;

        v.nic1 = nic1; v.nic2 = nic2; v.nii2 = nii2; v.nii3 = nii3;
        v.nir2 = nir2; v.nir3 = nir3; v.nis2 = nis2; v.nig2 = nig2;
        v.nig3 = nig3; v.niIN2 = niIN2;

        v.r_c[0] = r_c[0];
        v.r_i[0] = r_i[0];
        v.r_r[0] = r_r[0];
        v.r_g[0] = r_g[0];
        v.r_s[0] = r_s[0];
        v.Nt_i[0] = Nt_i[0];
        v.Nt_IN[0] = Nt_IN[0];

        for (int k = 0; k < 5; ++k) {
            for (int n = 0; n < 15; ++n) {
                v.cce[k][n] = cce[k][n];
                v.ccg[k][n] = ccg[k][n];
            }
        }
        for (int n = 0; n < 15; ++n) {
            v.ocg1[n] = ocg1[n];
            v.ocg2[n] = ocg2[n];
        }
        for (int i = 0; i < 7; ++i) {
            v.cie[i] = cie[i];
            v.cig[i] = cig[i];
        }
        v.oig1 = oig1; v.oig2 = oig2; v.obmi = obmi;

        for (int i = 0; i < 13; ++i) {
            v.cre[i] = cre[i];
            v.crg[i] = crg[i];
        }
        v.ore1 = ore1; v.org1 = org1; v.org2 = org2; v.org3 = org3; v.obmr = obmr;

        for (int i = 0; i < 17; ++i) {
            v.cse[i] = cse[i];
            v.csg[i] = csg[i];
        }
        v.oams = oams; v.obms = obms; v.ocms = ocms;

        for (int i = 0; i < 12; ++i) {
            v.cge[i] = cge[i];
            v.cgg[i] = cgg[i];
        }
        v.oamg = oamg; v.obmg = obmg; v.ocmg = ocmg;
        v.oge1 = oge1; v.ogg1 = ogg1; v.ogg2 = ogg2; v.ogg3 = ogg3;

        for (int i = 0; i < 10; ++i) {
            v.sa[i] = sa[i];
            v.sb[i] = sb[i];
        }

        v.Dr_0 = Dr_0;
        v.log_Dr = log_Dr;
        v.inv_log_Dr = inv_log_Dr;
        v.Ds_0 = Ds_0;
        v.log_Ds = log_Ds;
        v.inv_log_Ds = inv_log_Ds;
        v.t_Nc_0 = t_Nc_0;
        v.inv_nic1 = inv_nic1;

        v.t1_qg_qc = t1_qg_qc;
        v.D0i = D0i;
        v.xm0s = xm0s;
        v.xm0g = xm0g;

        return v;
    }
};

// Device-resident table container owning a single contiguous device allocation for all 28 tables
template <typename Real>
struct ThompsonDeviceTables {
    Real* buf_{nullptr};
    ThompsonTablesView<Real> view_{};

    ThompsonDeviceTables() = default;
    ~ThompsonDeviceTables() {
        if (buf_) dev_free(buf_);
    }
    ThompsonDeviceTables(const ThompsonDeviceTables&) = delete;
    ThompsonDeviceTables& operator=(const ThompsonDeviceTables&) = delete;
    ThompsonDeviceTables(ThompsonDeviceTables&& o) noexcept : buf_(o.buf_), view_(o.view_) {
        o.buf_ = nullptr;
    }
    ThompsonDeviceTables& operator=(ThompsonDeviceTables&& o) noexcept {
        if (this != &o) {
            if (buf_) dev_free(buf_);
            buf_ = o.buf_;
            view_ = o.view_;
            o.buf_ = nullptr;
        }
        return *this;
    }

    void upload(const ThompsonHostTables<Real>& h) {
        if (buf_) dev_free(buf_);
        constexpr size_t total = 5 * N_RACG + 12 * N_RACS + 4 * N_QRFZ + 2 * N_QCFZ +
                                 3 * N_IAUS + N_EFRW + N_EFSW + 2 * N_WEV;
        buf_ = static_cast<Real*>(dev_alloc(total * sizeof(Real)));
        view_ = h.view();
        size_t off = 0;
        auto up = [&](const std::vector<Real>& src) -> const Real* {
            Real* dst = buf_ + off;
            copy_to_dev(dst, src.data(), src.size() * sizeof(Real));
            off += src.size();
            return dst;
        };
        view_.tcg_racg = up(h.tcg_racg);
        view_.tmr_racg = up(h.tmr_racg);
        view_.tcr_gacr = up(h.tcr_gacr);
        view_.tnr_racg = up(h.tnr_racg);
        view_.tnr_gacr = up(h.tnr_gacr);

        view_.tcs_racs1 = up(h.tcs_racs1);
        view_.tmr_racs1 = up(h.tmr_racs1);
        view_.tcs_racs2 = up(h.tcs_racs2);
        view_.tmr_racs2 = up(h.tmr_racs2);
        view_.tcr_sacr1 = up(h.tcr_sacr1);
        view_.tms_sacr1 = up(h.tms_sacr1);
        view_.tcr_sacr2 = up(h.tcr_sacr2);
        view_.tms_sacr2 = up(h.tms_sacr2);
        view_.tnr_racs1 = up(h.tnr_racs1);
        view_.tnr_racs2 = up(h.tnr_racs2);
        view_.tnr_sacr1 = up(h.tnr_sacr1);
        view_.tnr_sacr2 = up(h.tnr_sacr2);

        view_.tpi_qrfz = up(h.tpi_qrfz);
        view_.tpg_qrfz = up(h.tpg_qrfz);
        view_.tni_qrfz = up(h.tni_qrfz);
        view_.tnr_qrfz = up(h.tnr_qrfz);

        view_.tpi_qcfz = up(h.tpi_qcfz);
        view_.tni_qcfz = up(h.tni_qcfz);

        view_.tps_iaus = up(h.tps_iaus);
        view_.tni_iaus = up(h.tni_iaus);
        view_.tpi_ide = up(h.tpi_ide);

        view_.t_Efrw = up(h.t_Efrw);
        view_.t_Efsw = up(h.t_Efsw);

        view_.tpc_wev = up(h.tpc_wev);
        view_.tnc_wev = up(h.tnc_wev);
    }

    const ThompsonTablesView<Real>& view() const { return view_; }
};

// WRF's table files when present, else the tables computed as thompson_init does; force_read (WRF
// force_read_thompson): a missing or unreadable table file is an error, as in qr_acr_qg / qr_acr_qs / freezeH2O
template <typename Real>
void init_thompson_tables(ThompsonHostTables<Real>& tables, const std::string& wrf_dat_dir = "",
                          bool force_read = false);

// Explicit instantiations declared
extern template void init_thompson_tables<double>(ThompsonHostTables<double>& tables, const std::string&, bool);
extern template void init_thompson_tables<float>(ThompsonHostTables<float>& tables, const std::string&, bool);

} // namespace cppwrf::phys::mp_thompson
