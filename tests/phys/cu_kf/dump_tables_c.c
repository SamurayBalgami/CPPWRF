// SPDX-License-Identifier: GPL-3.0-or-later
/* Writes the Kain-Fritsch lookup tables of WRF's module_cu_kfeta.o after kf_lutab (dump_tables_f.f90) to
   tables_fortran.bin. The tables are PRIVATE module variables; gfortran still gives them global symbols
   (__module_cu_kfeta_MOD_<name>), which this file reads. REAL is real*8 in WRF's -r8 build. */
#include <stdio.h>

extern double __module_cu_kfeta_MOD_rdpr, __module_cu_kfeta_MOD_rdthk, __module_cu_kfeta_MOD_plutop;
extern double __module_cu_kfeta_MOD_the0k[220];
extern double __module_cu_kfeta_MOD_alu[200];
extern double __module_cu_kfeta_MOD_ttab[220][250]; /* TTAB(KFNT=250, KFNP=220) */
extern double __module_cu_kfeta_MOD_qstab[220][250];

void dump_kf_tables_(void) {
  FILE* f = fopen("tables_fortran.bin", "wb");
  if (!f) {
    perror("tables_fortran.bin");
    return;
  }
  fwrite(&__module_cu_kfeta_MOD_rdpr, sizeof(double), 1, f);
  fwrite(&__module_cu_kfeta_MOD_rdthk, sizeof(double), 1, f);
  fwrite(&__module_cu_kfeta_MOD_plutop, sizeof(double), 1, f);
  fwrite(__module_cu_kfeta_MOD_the0k, sizeof(double), 220, f);
  fwrite(__module_cu_kfeta_MOD_alu, sizeof(double), 200, f);
  fwrite(__module_cu_kfeta_MOD_ttab, sizeof(double), 220 * 250, f);
  fwrite(__module_cu_kfeta_MOD_qstab, sizeof(double), 220 * 250, f);
  fclose(f);
}
