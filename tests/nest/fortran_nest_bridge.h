#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void wrf_bridge_init_constants();

void wrf_bridge_sint(double* xf, int ims, int ime, int jms, int jme, const int* icmask_int,
                     int its, int ite, int jts, int jte, int nf, int xstag_int, int ystag_int);

void wrf_bridge_sintb(const double* xf1, double* xf, int ims, int ime, int jms, int jme,
                      const int* icmask_int, int its, int ite, int jts, int jte, int nf,
                      int xstag_int, int ystag_int);

void wrf_bridge_interp_fcn_sint(const double* cfld, int cids, int cide, int ckds, int ckde, int cjds, int cjde,
                                int cims, int cime, int ckms, int ckme, int cjms, int cjme,
                                int cits, int cite, int ckts, int ckte, int cjts, int cjte,
                                double* nfld, int nids, int nide, int nkds, int nkde, int njds, int njde,
                                int nims, int nime, int nkms, int nkme, int njms, int njme,
                                int nits, int nite, int nkts, int nkte, int njts, int njte,
                                int shw, const int* imask, int xstag_int, int ystag_int,
                                int ipos, int jpos, int nri, int nrj);

void wrf_bridge_bdy_interp1(const double* cfld, int cids, int cide, int ckds, int ckde, int cjds, int cjde,
                            int cims, int cime, int ckms, int ckme, int cjms, int cjme,
                            int cits, int cite, int ckts, int ckte, int cjts, int cjte,
                            const double* nfld, int spec_bdy_width,
                            int nids, int nide, int nkds, int nkde, int njds, int njde,
                            int nims, int nime, int nkms, int nkme, int njms, int njme,
                            int nits, int nite, int nkts, int nkte, int njts, int njte,
                            int shw1, const int* imask, int xstag_int, int ystag_int,
                            int ipos, int jpos, int nri, int nrj,
                            double* bdy_xs, double* bdy_txs, double* bdy_xe, double* bdy_txe,
                            double* bdy_ys, double* bdy_tys, double* bdy_ye, double* bdy_tye,
                            double cdt, double ndt);

void wrf_bridge_copy_fcn(double* cfld, int cids, int cide, int ckds, int ckde, int cjds, int cjde,
                         int cims, int cime, int ckms, int ckme, int cjms, int cjme,
                         int cits, int cite, int ckts, int ckte, int cjts, int cjte,
                         const double* nfld, int nids, int nide, int nkds, int nkde, int njds, int njde,
                         int nims, int nime, int nkms, int nkme, int njms, int njme,
                         int nits, int nite, int nkts, int nkte, int njts, int njte,
                         int shw, const int* imask, int xstag_int, int ystag_int,
                         int ipos, int jpos, int nri, int nrj);

void wrf_bridge_sm121(double* cfld, int cids, int cide, int ckds, int ckde, int cjds, int cjde,
                      int cims, int cime, int ckms, int ckme, int cjms, int cjme,
                      int cits, int cite, int ckts, int ckte, int cjts, int cjte,
                      int xstag_int, int ystag_int,
                      int nids, int nide, int nkds, int nkde, int njds, int njde,
                      int nims, int nime, int nkms, int nkme, int njms, int njme,
                      int nits, int nite, int nkts, int nkte, int njts, int njte,
                      int nri, int nrj, int ipos, int jpos);

void wrf_bridge_smdsm(double* cfld, int cids, int cide, int ckds, int ckde, int cjds, int cjde,
                      int cims, int cime, int ckms, int ckme, int cjms, int cjme,
                      int cits, int cite, int ckts, int ckte, int cjts, int cjte,
                      int xstag_int, int ystag_int,
                      int nids, int nide, int nkds, int nkde, int njds, int njde,
                      int nims, int nime, int nkms, int nkme, int njms, int njme,
                      int nits, int nite, int nkts, int nkte, int njts, int njte,
                      int nri, int nrj, int ipos, int jpos);

void wrf_bridge_blend_terrain(const double* ter_interpolated, double* ter_input,
                              int ids, int ide, int jds, int jde, int kds, int kde,
                              int ims, int ime, int jms, int jme, int kms, int kme,
                              int ips, int ipe, int jps, int jpe, int kps, int kpe);

void wrf_bridge_copy_3d_field(double* ter_interpolated, const double* ter_input,
                              int ids, int ide, int jds, int jde, int kds, int kde,
                              int ims, int ime, int jms, int jme, int kms, int kme,
                              int ips, int ipe, int jps, int jpe, int kps, int kpe);

void wrf_bridge_adjust_tempqv(const double* mub, const double* save_mub,
                              const double* c3, const double* c4, const double* znw, double p_top,
                              double* th, const double* pp, double* qv, int use_theta_m,
                              int ids, int ide, int jds, int jde, int kds, int kde,
                              int ims, int ime, int jms, int jme, int kms, int kme,
                              int ips, int ipe, int jps, int jpe, int kps, int kpe);

void wrf_bridge_couple_or_uncouple_em(double* u_2, double* v_2, double* w_2, double* ph_2, double* t_2, double* qv,
                                      const double* mub, const double* mu_2,
                                      const double* msfuy, const double* msfvx, const double* msfty,
                                      const double* c1f, const double* c2f,
                                      const double* c1h, const double* c2h,
                                      int nx, int ny, int nz, int couple_int);

#ifdef __cplusplus
}
#endif
