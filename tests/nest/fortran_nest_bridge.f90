! SPDX-License-Identifier: GPL-3.0-or-later
! tests/nest/fortran_nest_bridge.f90
! Fortran-to-C bridge for WRF 4.8.0 nesting routines
MODULE wrf_nest_bridge
  USE ISO_C_BINDING
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE nl_set_spec_zone ( id, val )
      INTEGER, INTENT(IN) :: id, val
    END SUBROUTINE nl_set_spec_zone
    SUBROUTINE nl_set_relax_zone ( id, val )
      INTEGER, INTENT(IN) :: id, val
    END SUBROUTINE nl_set_relax_zone
    SUBROUTINE nl_set_spec_bdy_width ( id, val )
      INTEGER, INTENT(IN) :: id, val
    END SUBROUTINE nl_set_spec_bdy_width
    SUBROUTINE nl_set_blend_width ( id, val )
      INTEGER, INTENT(IN) :: id, val
    END SUBROUTINE nl_set_blend_width

    SUBROUTINE sint ( xf, ims, ime, jms, jme, icmask, its, ite, jts, jte, nf, xstag, ystag )
      USE ISO_C_BINDING
      INTEGER :: ims, ime, jms, jme, its, ite, jts, jte, nf
      LOGICAL :: xstag, ystag
      LOGICAL :: icmask(ims:ime, jms:jme)
      REAL(C_DOUBLE) :: xf(ims:ime, jms:jme, nf)
    END SUBROUTINE sint

    SUBROUTINE sintb ( xf1, xf, ims, ime, jms, jme, icmask, its, ite, jts, jte, nf, xstag, ystag )
      USE ISO_C_BINDING
      INTEGER :: ims, ime, jms, jme, its, ite, jts, jte, nf
      LOGICAL :: xstag, ystag
      LOGICAL :: icmask(ims:ime, jms:jme)
      REAL(C_DOUBLE) :: xf1(ims:ime, jms:jme, nf), xf(ims:ime, jms:jme, nf)
    END SUBROUTINE sintb

    SUBROUTINE interp_fcn_sint ( cfld, cids, cide, ckds, ckde, cjds, cjde, &
                                 cims, cime, ckms, ckme, cjms, cjme, &
                                 cits, cite, ckts, ckte, cjts, cjte, &
                                 nfld, nids, nide, nkds, nkde, njds, njde, &
                                 nims, nime, nkms, nkme, njms, njme, &
                                 nits, nite, nkts, nkte, njts, njte, &
                                 shw, imask, xstag, ystag, ipos, jpos, nri, nrj )
      USE ISO_C_BINDING
      INTEGER :: cids, cide, ckds, ckde, cjds, cjde, &
                 cims, cime, ckms, ckme, cjms, cjme, &
                 cits, cite, ckts, ckte, cjts, cjte, &
                 nids, nide, nkds, nkde, njds, njde, &
                 nims, nime, nkms, nkme, njms, njme, &
                 nits, nite, nkts, nkte, njts, njte, &
                 shw, ipos, jpos, nri, nrj
      LOGICAL :: xstag, ystag
      REAL(C_DOUBLE) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
      REAL(C_DOUBLE) :: nfld(nims:nime, nkms:nkme, njms:njme)
      INTEGER :: imask(nims:nime, njms:njme)
    END SUBROUTINE interp_fcn_sint

    SUBROUTINE bdy_interp1 ( cfld, cids, cide, ckds, ckde, cjds, cjde, &
                            cims, cime, ckms, ckme, cjms, cjme, &
                            cits, cite, ckts, ckte, cjts, cjte, &
                            nfld, nijds, nijde, spec_bdy_width, &
                            nids, nide, nkds, nkde, njds, njde, &
                            nims, nime, nkms, nkme, njms, njme, &
                            nits, nite, nkts, nkte, njts, njte, &
                            shw1, imask, xstag, ystag, ipos, jpos, nri, nrj, &
                            cbdy_xs, bdy_xs, cbdy_xe, bdy_xe, &
                            cbdy_ys, bdy_ys, cbdy_ye, bdy_ye, &
                            cbdy_txs, bdy_txs, cbdy_txe, bdy_txe, &
                            cbdy_tys, bdy_tys, cbdy_tye, bdy_tye, &
                            cdt, ndt )
      USE ISO_C_BINDING
      INTEGER :: cids, cide, ckds, ckde, cjds, cjde, &
                 cims, cime, ckms, ckme, cjms, cjme, &
                 cits, cite, ckts, ckte, cjts, cjte, &
                 nijds, nijde, spec_bdy_width, &
                 nids, nide, nkds, nkde, njds, njde, &
                 nims, nime, nkms, nkme, njms, njme, &
                 nits, nite, nkts, nkte, njts, njte, &
                 shw1, ipos, jpos, nri, nrj
      LOGICAL :: xstag, ystag
      REAL(C_DOUBLE) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
      REAL(C_DOUBLE) :: nfld(nims:nime, nkms:nkme, njms:njme)
      INTEGER :: imask(nims:nime, njms:njme)
      REAL(C_DOUBLE) :: cbdy_xs(*), cbdy_txs(*)
      REAL(C_DOUBLE) :: cbdy_xe(*), cbdy_txe(*)
      REAL(C_DOUBLE) :: cbdy_ys(*), cbdy_tys(*)
      REAL(C_DOUBLE) :: cbdy_ye(*), cbdy_tye(*)
      REAL(C_DOUBLE) :: bdy_xs(njms:njme, nkms:nkme, spec_bdy_width), bdy_txs(njms:njme, nkms:nkme, spec_bdy_width)
      REAL(C_DOUBLE) :: bdy_xe(njms:njme, nkms:nkme, spec_bdy_width), bdy_txe(njms:njme, nkms:nkme, spec_bdy_width)
      REAL(C_DOUBLE) :: bdy_ys(nims:nime, nkms:nkme, spec_bdy_width), bdy_tys(nims:nime, nkms:nkme, spec_bdy_width)
      REAL(C_DOUBLE) :: bdy_ye(nims:nime, nkms:nkme, spec_bdy_width), bdy_tye(nims:nime, nkms:nkme, spec_bdy_width)
      REAL(C_DOUBLE) :: cdt, ndt
    END SUBROUTINE bdy_interp1

    SUBROUTINE copy_fcn ( cfld, cids, cide, ckds, ckde, cjds, cjde, &
                          cims, cime, ckms, ckme, cjms, cjme, &
                          cits, cite, ckts, ckte, cjts, cjte, &
                          nfld, nids, nide, nkds, nkde, njds, njde, &
                          nims, nime, nkms, nkme, njms, njme, &
                          nits, nite, nkts, nkte, njts, njte, &
                          shw, imask, xstag, ystag, ipos, jpos, nri, nrj )
      USE ISO_C_BINDING
      INTEGER :: cids, cide, ckds, ckde, cjds, cjde, &
                 cims, cime, ckms, ckme, cjms, cjme, &
                 cits, cite, ckts, ckte, cjts, cjte, &
                 nids, nide, nkds, nkde, njds, njde, &
                 nims, nime, nkms, nkme, njms, njme, &
                 nits, nite, nkts, nkte, njts, njte, &
                 shw, ipos, jpos, nri, nrj
      LOGICAL :: xstag, ystag
      REAL(C_DOUBLE) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
      REAL(C_DOUBLE) :: nfld(nims:nime, nkms:nkme, njms:njme)
      INTEGER :: imask(nims:nime, njms:njme)
    END SUBROUTINE copy_fcn

    SUBROUTINE sm121 ( cfld, cids, cide, ckds, ckde, cjds, cjde, &
                       cims, cime, ckms, ckme, cjms, cjme, &
                       cits, cite, ckts, ckte, cjts, cjte, &
                       xstag, ystag, &
                       nids, nide, nkds, nkde, njds, njde, &
                       nims, nime, nkms, nkme, njms, njme, &
                       nits, nite, nkts, nkte, njts, njte, &
                       nri, nrj, ipos, jpos )
      USE ISO_C_BINDING
      INTEGER :: cids, cide, ckds, ckde, cjds, cjde, &
                 cims, cime, ckms, ckme, cjms, cjme, &
                 cits, cite, ckts, ckte, cjts, cjte, &
                 nids, nide, nkds, nkde, njds, njde, &
                 nims, nime, nkms, nkme, njms, njme, &
                 nits, nite, nkts, nkte, njts, njte, &
                 nri, nrj, ipos, jpos
      LOGICAL :: xstag, ystag
      REAL(C_DOUBLE) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
    END SUBROUTINE sm121

    SUBROUTINE smdsm ( cfld, cids, cide, ckds, ckde, cjds, cjde, &
                       cims, cime, ckms, ckme, cjms, cjme, &
                       cits, cite, ckts, ckte, cjts, cjte, &
                       xstag, ystag, &
                       nids, nide, nkds, nkde, njds, njde, &
                       nims, nime, nkms, nkme, njms, njme, &
                       nits, nite, nkts, nkte, njts, njte, &
                       nri, nrj, ipos, jpos )
      USE ISO_C_BINDING
      INTEGER :: cids, cide, ckds, ckde, cjds, cjde, &
                 cims, cime, ckms, ckme, cjms, cjme, &
                 cits, cite, ckts, ckte, cjts, cjte, &
                 nids, nide, nkds, nkde, njds, njde, &
                 nims, nime, nkms, nkme, njms, njme, &
                 nits, nite, nkts, nkte, njts, njte, &
                 nri, nrj, ipos, jpos
      LOGICAL :: xstag, ystag
      REAL(C_DOUBLE) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
    END SUBROUTINE smdsm

    SUBROUTINE blend_terrain ( ter_interpolated, ter_input, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               ips, ipe, jps, jpe, kps, kpe )
      USE ISO_C_BINDING
      INTEGER :: ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 ips, ipe, jps, jpe, kps, kpe
      REAL(C_DOUBLE) :: ter_interpolated(ims:ime, kms:kme, jms:jme)
      REAL(C_DOUBLE) :: ter_input(ims:ime, kms:kme, jms:jme)
    END SUBROUTINE blend_terrain

    SUBROUTINE copy_3d_field ( ter_interpolated, ter_input, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               ips, ipe, jps, jpe, kps, kpe )
      USE ISO_C_BINDING
      INTEGER :: ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 ips, ipe, jps, jpe, kps, kpe
      REAL(C_DOUBLE) :: ter_interpolated(ims:ime, kms:kme, jms:jme)
      REAL(C_DOUBLE) :: ter_input(ims:ime, kms:kme, jms:jme)
    END SUBROUTINE copy_3d_field

    SUBROUTINE adjust_tempqv ( mub, save_mub, c3, c4, znw, p_top, &
                               th, pp, qv, use_theta_m, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               ips, ipe, jps, jpe, kps, kpe )
      USE ISO_C_BINDING
      INTEGER :: ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 ips, ipe, jps, jpe, kps, kpe, use_theta_m
      REAL(C_DOUBLE) :: mub(ims:ime, jms:jme), save_mub(ims:ime, jms:jme)
      REAL(C_DOUBLE) :: c3(kms:kme), c4(kms:kme), znw(kms:kme), p_top
      REAL(C_DOUBLE) :: th(ims:ime, kms:kme, jms:jme)
      REAL(C_DOUBLE) :: pp(ims:ime, kms:kme, jms:jme)
      REAL(C_DOUBLE) :: qv(ims:ime, kms:kme, jms:jme)
    END SUBROUTINE adjust_tempqv
  END INTERFACE

CONTAINS

  SUBROUTINE wrf_bridge_init_constants() BIND(C, name="wrf_bridge_init_constants")
    CALL nl_set_spec_zone(1, 1)
    CALL nl_set_relax_zone(1, 4)
    CALL nl_set_spec_bdy_width(1, 5)
    CALL nl_set_blend_width(1, 5)
  END SUBROUTINE wrf_bridge_init_constants

  SUBROUTINE wrf_bridge_sint(xf, ims, ime, jms, jme, icmask_int, its, ite, jts, jte, nf, xstag_int, ystag_int) &
      BIND(C, name="wrf_bridge_sint")
    INTEGER(C_INT), VALUE :: ims, ime, jms, jme, its, ite, jts, jte, nf, xstag_int, ystag_int
    REAL(C_DOUBLE), INTENT(INOUT) :: xf(ims:ime, jms:jme, nf)
    INTEGER(C_INT), INTENT(IN) :: icmask_int(ims:ime, jms:jme)
    LOGICAL :: l_icmask(ims:ime, jms:jme)
    LOGICAL :: xstag, ystag
    INTEGER :: i, j

    DO j = jms, jme
      DO i = ims, ime
        l_icmask(i, j) = (icmask_int(i, j) /= 0)
      END DO
    END DO
    xstag = (xstag_int /= 0)
    ystag = (ystag_int /= 0)

    CALL sint(xf, ims, ime, jms, jme, l_icmask, its, ite, jts, jte, nf, xstag, ystag)
  END SUBROUTINE wrf_bridge_sint

  SUBROUTINE wrf_bridge_sintb(xf1, xf, ims, ime, jms, jme, icmask_int, its, ite, jts, jte, nf, xstag_int, ystag_int) &
      BIND(C, name="wrf_bridge_sintb")
    INTEGER(C_INT), VALUE :: ims, ime, jms, jme, its, ite, jts, jte, nf, xstag_int, ystag_int
    REAL(C_DOUBLE), INTENT(IN) :: xf1(ims:ime, jms:jme, nf)
    REAL(C_DOUBLE), INTENT(OUT) :: xf(ims:ime, jms:jme, nf)
    INTEGER(C_INT), INTENT(IN) :: icmask_int(ims:ime, jms:jme)
    LOGICAL :: l_icmask(ims:ime, jms:jme)
    LOGICAL :: xstag, ystag
    INTEGER :: i, j

    DO j = jms, jme
      DO i = ims, ime
        l_icmask(i, j) = (icmask_int(i, j) /= 0)
      END DO
    END DO
    xstag = (xstag_int /= 0)
    ystag = (ystag_int /= 0)

    CALL sintb(xf1, xf, ims, ime, jms, jme, l_icmask, its, ite, jts, jte, nf, xstag, ystag)
  END SUBROUTINE wrf_bridge_sintb

  SUBROUTINE wrf_bridge_interp_fcn_sint(cfld, cids, cide, ckds, ckde, cjds, cjde, &
                                        cims, cime, ckms, ckme, cjms, cjme, &
                                        cits, cite, ckts, ckte, cjts, cjte, &
                                        nfld, nids, nide, nkds, nkde, njds, njde, &
                                        nims, nime, nkms, nkme, njms, njme, &
                                        nits, nite, nkts, nkte, njts, njte, &
                                        shw, imask, xstag_int, ystag_int, ipos, jpos, nri, nrj) &
      BIND(C, name="wrf_bridge_interp_fcn_sint")
    INTEGER(C_INT), VALUE :: cids, cide, ckds, ckde, cjds, cjde, &
                             cims, cime, ckms, ckme, cjms, cjme, &
                             cits, cite, ckts, ckte, cjts, cjte, &
                             nids, nide, nkds, nkde, njds, njde, &
                             nims, nime, nkms, nkme, njms, njme, &
                             nits, nite, nkts, nkte, njts, njte, &
                             shw, xstag_int, ystag_int, ipos, jpos, nri, nrj
    REAL(C_DOUBLE), INTENT(IN) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
    REAL(C_DOUBLE), INTENT(INOUT) :: nfld(nims:nime, nkms:nkme, njms:njme)
    INTEGER(C_INT), INTENT(IN) :: imask(nims:nime, njms:njme)
    LOGICAL :: xstag, ystag

    xstag = (xstag_int /= 0)
    ystag = (ystag_int /= 0)

    CALL interp_fcn_sint(cfld, cids, cide, ckds, ckde, cjds, cjde, &
                         cims, cime, ckms, ckme, cjms, cjme, &
                         cits, cite, ckts, ckte, cjts, cjte, &
                         nfld, nids, nide, nkds, nkde, njds, njde, &
                         nims, nime, nkms, nkme, njms, njme, &
                         nits, nite, nkts, nkte, njts, njte, &
                         shw, imask, xstag, ystag, ipos, jpos, nri, nrj)
  END SUBROUTINE wrf_bridge_interp_fcn_sint

  SUBROUTINE wrf_bridge_bdy_interp1(cfld, cids, cide, ckds, ckde, cjds, cjde, &
                                    cims, cime, ckms, ckme, cjms, cjme, &
                                    cits, cite, ckts, ckte, cjts, cjte, &
                                    nfld, spec_bdy_width, &
                                    nids, nide, nkds, nkde, njds, njde, &
                                    nims, nime, nkms, nkme, njms, njme, &
                                    nits, nite, nkts, nkte, njts, njte, &
                                    shw1, imask, xstag_int, ystag_int, ipos, jpos, nri, nrj, &
                                    bdy_xs, bdy_txs, bdy_xe, bdy_txe, &
                                    bdy_ys, bdy_tys, bdy_ye, bdy_tye, &
                                    cdt, ndt) BIND(C, name="wrf_bridge_bdy_interp1")
    INTEGER(C_INT), VALUE :: cids, cide, ckds, ckde, cjds, cjde, &
                             cims, cime, ckms, ckme, cjms, cjme, &
                             cits, cite, ckts, ckte, cjts, cjte, &
                             spec_bdy_width, &
                             nids, nide, nkds, nkde, njds, njde, &
                             nims, nime, nkms, nkme, njms, njme, &
                             nits, nite, nkts, nkte, njts, njte, &
                             shw1, xstag_int, ystag_int, ipos, jpos, nri, nrj
    REAL(C_DOUBLE), INTENT(IN) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
    REAL(C_DOUBLE), INTENT(IN) :: nfld(nims:nime, nkms:nkme, njms:njme)
    INTEGER(C_INT), INTENT(IN) :: imask(nims:nime, njms:njme)
    REAL(C_DOUBLE), INTENT(OUT) :: bdy_xs(njms:njme, nkms:nkme, spec_bdy_width), bdy_txs(njms:njme, nkms:nkme, spec_bdy_width)
    REAL(C_DOUBLE), INTENT(OUT) :: bdy_xe(njms:njme, nkms:nkme, spec_bdy_width), bdy_txe(njms:njme, nkms:nkme, spec_bdy_width)
    REAL(C_DOUBLE), INTENT(OUT) :: bdy_ys(nims:nime, nkms:nkme, spec_bdy_width), bdy_tys(nims:nime, nkms:nkme, spec_bdy_width)
    REAL(C_DOUBLE), INTENT(OUT) :: bdy_ye(nims:nime, nkms:nkme, spec_bdy_width), bdy_tye(nims:nime, nkms:nkme, spec_bdy_width)
    REAL(C_DOUBLE), VALUE :: cdt, ndt

    LOGICAL :: xstag, ystag
    INTEGER :: nijds, nijde
    REAL(C_DOUBLE) :: dummy(1)

    REAL(C_DOUBLE), ALLOCATABLE :: t_xs(:,:,:), t_txs(:,:,:)
    REAL(C_DOUBLE), ALLOCATABLE :: t_xe(:,:,:), t_txe(:,:,:)
    REAL(C_DOUBLE), ALLOCATABLE :: t_ys(:,:,:), t_tys(:,:,:)
    REAL(C_DOUBLE), ALLOCATABLE :: t_ye(:,:,:), t_tye(:,:,:)
    REAL(C_DOUBLE), ALLOCATABLE :: t_nfld(:,:,:)
    INTEGER(C_INT), ALLOCATABLE :: t_imask(:,:)

    xstag = (xstag_int /= 0)
    ystag = (ystag_int /= 0)
    nijds = min(nids, njds)
    nijde = max(nide, njde)

    ALLOCATE(t_xs(njms:njme+2, nkms:nkme, 1:spec_bdy_width))
    ALLOCATE(t_txs(njms:njme+2, nkms:nkme, 1:spec_bdy_width))
    ALLOCATE(t_xe(njms:njme+2, nkms:nkme, 1:spec_bdy_width))
    ALLOCATE(t_txe(njms:njme+2, nkms:nkme, 1:spec_bdy_width))
    ALLOCATE(t_ys(nims:nime+2, nkms:nkme, 1:spec_bdy_width))
    ALLOCATE(t_tys(nims:nime+2, nkms:nkme, 1:spec_bdy_width))
    ALLOCATE(t_ye(nims:nime+2, nkms:nkme, 1:spec_bdy_width))
    ALLOCATE(t_tye(nims:nime+2, nkms:nkme, 1:spec_bdy_width))
    ALLOCATE(t_nfld(nims:nime+2, nkms:nkme, njms:njme+2))
    ALLOCATE(t_imask(nims:nime+2, njms:njme+2))

    t_xs = 0.0_C_DOUBLE
    t_txs = 0.0_C_DOUBLE
    t_xe = 0.0_C_DOUBLE
    t_txe = 0.0_C_DOUBLE
    t_ys = 0.0_C_DOUBLE
    t_tys = 0.0_C_DOUBLE
    t_ye = 0.0_C_DOUBLE
    t_tye = 0.0_C_DOUBLE
    t_nfld = 0.0_C_DOUBLE
    t_imask = 0
    t_nfld(nims:nime, nkms:nkme, njms:njme) = nfld(nims:nime, nkms:nkme, njms:njme)
    t_imask(nims:nime, njms:njme) = imask(nims:nime, njms:njme)

    CALL bdy_interp1(cfld, cids, cide, ckds, ckde, cjds, cjde, &
                     cims, cime, ckms, ckme, cjms, cjme, &
                     cits, cite, ckts, ckte, cjts, cjte, &
                     t_nfld, nijds, nijde, spec_bdy_width, &
                     nids, nide, nkds, nkde, njds, njde, &
                     nims, nime+2, nkms, nkme, njms, njme+2, &
                     nits, nite, nkts, nkte, njts, njte, &
                     shw1, t_imask, xstag, ystag, ipos, jpos, nri, nrj, &
                     dummy, t_xs, dummy, t_xe, dummy, t_ys, dummy, t_ye, &
                     dummy, t_txs, dummy, t_txe, dummy, t_tys, dummy, t_tye, &
                     cdt, ndt)

    bdy_xs(njms:njme, nkms:nkme, 1:spec_bdy_width) = t_xs(njms:njme, nkms:nkme, 1:spec_bdy_width)
    bdy_txs(njms:njme, nkms:nkme, 1:spec_bdy_width) = t_txs(njms:njme, nkms:nkme, 1:spec_bdy_width)
    bdy_xe(njms:njme, nkms:nkme, 1:spec_bdy_width) = t_xe(njms:njme, nkms:nkme, 1:spec_bdy_width)
    bdy_txe(njms:njme, nkms:nkme, 1:spec_bdy_width) = t_txe(njms:njme, nkms:nkme, 1:spec_bdy_width)
    bdy_ys(nims:nime, nkms:nkme, 1:spec_bdy_width) = t_ys(nims:nime, nkms:nkme, 1:spec_bdy_width)
    bdy_tys(nims:nime, nkms:nkme, 1:spec_bdy_width) = t_tys(nims:nime, nkms:nkme, 1:spec_bdy_width)
    bdy_ye(nims:nime, nkms:nkme, 1:spec_bdy_width) = t_ye(nims:nime, nkms:nkme, 1:spec_bdy_width)
    bdy_tye(nims:nime, nkms:nkme, 1:spec_bdy_width) = t_tye(nims:nime, nkms:nkme, 1:spec_bdy_width)

    DEALLOCATE(t_xs, t_txs, t_xe, t_txe, t_ys, t_tys, t_ye, t_tye, t_nfld, t_imask)
  END SUBROUTINE wrf_bridge_bdy_interp1

  SUBROUTINE wrf_bridge_copy_fcn(cfld, cids, cide, ckds, ckde, cjds, cjde, &
                                 cims, cime, ckms, ckme, cjms, cjme, &
                                 cits, cite, ckts, ckte, cjts, cjte, &
                                 nfld, nids, nide, nkds, nkde, njds, njde, &
                                 nims, nime, nkms, nkme, njms, njme, &
                                 nits, nite, nkts, nkte, njts, njte, &
                                 shw, imask, xstag_int, ystag_int, ipos, jpos, nri, nrj) &
      BIND(C, name="wrf_bridge_copy_fcn")
    INTEGER(C_INT), VALUE :: cids, cide, ckds, ckde, cjds, cjde, &
                             cims, cime, ckms, ckme, cjms, cjme, &
                             cits, cite, ckts, ckte, cjts, cjte, &
                             nids, nide, nkds, nkde, njds, njde, &
                             nims, nime, nkms, nkme, njms, njme, &
                             nits, nite, nkts, nkte, njts, njte, &
                             shw, xstag_int, ystag_int, ipos, jpos, nri, nrj
    REAL(C_DOUBLE), INTENT(OUT) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
    REAL(C_DOUBLE), INTENT(IN) :: nfld(nims:nime, nkms:nkme, njms:njme)
    INTEGER(C_INT), INTENT(IN) :: imask(nims:nime, njms:njme)
    LOGICAL :: xstag, ystag

    xstag = (xstag_int /= 0)
    ystag = (ystag_int /= 0)

    CALL copy_fcn(cfld, cids, cide, ckds, ckde, cjds, cjde, &
                  cims, cime, ckms, ckme, cjms, cjme, &
                  cits, cite, ckts, ckte, cjts, cjte, &
                  nfld, nids, nide, nkds, nkde, njds, njde, &
                  nims, nime, nkms, nkme, njms, njme, &
                  nits, nite, nkts, nkte, njts, njte, &
                  shw, imask, xstag, ystag, ipos, jpos, nri, nrj)
  END SUBROUTINE wrf_bridge_copy_fcn

  SUBROUTINE wrf_bridge_sm121(cfld, cids, cide, ckds, ckde, cjds, cjde, &
                              cims, cime, ckms, ckme, cjms, cjme, &
                              cits, cite, ckts, ckte, cjts, cjte, &
                              xstag_int, ystag_int, &
                              nids, nide, nkds, nkde, njds, njde, &
                              nims, nime, nkms, nkme, njms, njme, &
                              nits, nite, nkts, nkte, njts, njte, &
                              nri, nrj, ipos, jpos) BIND(C, name="wrf_bridge_sm121")
    INTEGER(C_INT), VALUE :: cids, cide, ckds, ckde, cjds, cjde, &
                             cims, cime, ckms, ckme, cjms, cjme, &
                             cits, cite, ckts, ckte, cjts, cjte, &
                             xstag_int, ystag_int, &
                             nids, nide, nkds, nkde, njds, njde, &
                             nims, nime, nkms, nkme, njms, njme, &
                             nits, nite, nkts, nkte, njts, njte, &
                             nri, nrj, ipos, jpos
    REAL(C_DOUBLE), INTENT(INOUT) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
    LOGICAL :: xstag, ystag

    xstag = (xstag_int /= 0)
    ystag = (ystag_int /= 0)

    CALL sm121(cfld, cids, cide, ckds, ckde, cjds, cjde, &
               cims, cime, ckms, ckme, cjms, cjme, &
               cits, cite, ckts, ckte, cjts, cjte, &
               xstag, ystag, &
               nids, nide, nkds, nkde, njds, njde, &
               nims, nime, nkms, nkme, njms, njme, &
               nits, nite, nkts, nkte, njts, njte, &
               nri, nrj, ipos, jpos)
  END SUBROUTINE wrf_bridge_sm121

  SUBROUTINE wrf_bridge_smdsm(cfld, cids, cide, ckds, ckde, cjds, cjde, &
                              cims, cime, ckms, ckme, cjms, cjme, &
                              cits, cite, ckts, ckte, cjts, cjte, &
                              xstag_int, ystag_int, &
                              nids, nide, nkds, nkde, njds, njde, &
                              nims, nime, nkms, nkme, njms, njme, &
                              nits, nite, nkts, nkte, njts, njte, &
                              nri, nrj, ipos, jpos) BIND(C, name="wrf_bridge_smdsm")
    INTEGER(C_INT), VALUE :: cids, cide, ckds, ckde, cjds, cjde, &
                             cims, cime, ckms, ckme, cjms, cjme, &
                             cits, cite, ckts, ckte, cjts, cjte, &
                             xstag_int, ystag_int, &
                             nids, nide, nkds, nkde, njds, njde, &
                             nims, nime, nkms, nkme, njms, njme, &
                             nits, nite, nkts, nkte, njts, njte, &
                             nri, nrj, ipos, jpos
    REAL(C_DOUBLE), INTENT(INOUT) :: cfld(cims:cime, ckms:ckme, cjms:cjme)
    LOGICAL :: xstag, ystag

    xstag = (xstag_int /= 0)
    ystag = (ystag_int /= 0)

    CALL smdsm(cfld, cids, cide, ckds, ckde, cjds, cjde, &
               cims, cime, ckms, ckme, cjms, cjme, &
               cits, cite, ckts, ckte, cjts, cjte, &
               xstag, ystag, &
               nids, nide, nkds, nkde, njds, njde, &
               nims, nime, nkms, nkme, njms, njme, &
               nits, nite, nkts, nkte, njts, njte, &
               nri, nrj, ipos, jpos)
  END SUBROUTINE wrf_bridge_smdsm

  SUBROUTINE wrf_bridge_blend_terrain(ter_interpolated, ter_input, &
                                      ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      ips, ipe, jps, jpe, kps, kpe) &
      BIND(C, name="wrf_bridge_blend_terrain")
    INTEGER(C_INT), VALUE :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             ips, ipe, jps, jpe, kps, kpe
    REAL(C_DOUBLE), INTENT(IN) :: ter_interpolated(ims:ime, kms:kme, jms:jme)
    REAL(C_DOUBLE), INTENT(INOUT) :: ter_input(ims:ime, kms:kme, jms:jme)

    CALL blend_terrain(ter_interpolated, ter_input, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       ips, ipe, jps, jpe, kps, kpe)
  END SUBROUTINE wrf_bridge_blend_terrain

  SUBROUTINE wrf_bridge_copy_3d_field(ter_interpolated, ter_input, &
                                      ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      ips, ipe, jps, jpe, kps, kpe) &
      BIND(C, name="wrf_bridge_copy_3d_field")
    INTEGER(C_INT), VALUE :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             ips, ipe, jps, jpe, kps, kpe
    REAL(C_DOUBLE), INTENT(OUT) :: ter_interpolated(ims:ime, kms:kme, jms:jme)
    REAL(C_DOUBLE), INTENT(IN) :: ter_input(ims:ime, kms:kme, jms:jme)

    CALL copy_3d_field(ter_interpolated, ter_input, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       ips, ipe, jps, jpe, kps, kpe)
  END SUBROUTINE wrf_bridge_copy_3d_field

  SUBROUTINE wrf_bridge_adjust_tempqv(mub, save_mub, c3, c4, znw, p_top, &
                                      th, pp, qv, use_theta_m, &
                                      ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      ips, ipe, jps, jpe, kps, kpe) &
      BIND(C, name="wrf_bridge_adjust_tempqv")
    INTEGER(C_INT), VALUE :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             ips, ipe, jps, jpe, kps, kpe, use_theta_m
    REAL(C_DOUBLE), INTENT(IN) :: mub(ims:ime, jms:jme), save_mub(ims:ime, jms:jme)
    REAL(C_DOUBLE), INTENT(IN) :: c3(kms:kme), c4(kms:kme), znw(kms:kme)
    REAL(C_DOUBLE), VALUE :: p_top
    REAL(C_DOUBLE), INTENT(INOUT) :: th(ims:ime, kms:kme, jms:jme)
    REAL(C_DOUBLE), INTENT(IN) :: pp(ims:ime, kms:kme, jms:jme)
    REAL(C_DOUBLE), INTENT(INOUT) :: qv(ims:ime, kms:kme, jms:jme)

    CALL adjust_tempqv(mub, save_mub, c3, c4, znw, p_top, &
                       th, pp, qv, use_theta_m, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       ips, ipe, jps, jpe, kps, kpe)
  END SUBROUTINE wrf_bridge_adjust_tempqv

END MODULE wrf_nest_bridge
