# compare WRF debug dycore checkpoints (vm_*.bin, Fortran memory arrays) with CPPWRF's (logical arrays)
import sys, numpy as np
wdir, cdir = sys.argv[1], sys.argv[2]
specs = {'vm_tend_rk': (['ru_tend', 'rv_tend', 'rw_tend', 'ph_tend', 't_tend', '2:mu_tend', 'ru', 'rv', 'php', 'alt', 'ww', '2:muu', '2:mut'],
                        ['ru_tend', 'rv_tend', 'rw_tend', 'ph_tend', 't_tend', 'ru', 'rv', 'php', 'alt', 'ww', 'muu', 'mut']),
         'vm_state_rk': (['u', 'v', 'w', 't', 'ph', '2:mu_2', '2:muts', 'qv', 'ww', 'ru_m', 'rv_m', 'ww_m'], ['u', 'v', 'w', 't', 'ph', 'qv', 'ww', 'ru_m', 'rv_m', 'ww_m', 'mu_2', 'muts']),
         'vm_p_rk': (['p', 'al'], ['p', 'al']),
         'vm_s0_rk': (['p', 'al', 't_save', 'u_save', 'v_save', 'ph', 'php', 'alt', 'cqu', 'pb', '2:muu', '2:mu_2'],
                      ['p', 'al', 't_save', 'u_save', 'v_save', 'ph', 'php', 'alt', 'cqu', 'muu', 'mu_2', 'pb']),
         'vm_s1_rk': (['u', 'v'], ['u', 'v']),
         'vm_s2_rk': (['t', 'ww', '2:mu_2', '2:muts'], ['t', 'ww', 'mu_2', 'muts']),
         'vm_s3_rk': (['w', 'ph', 'p', 'al'], ['w', 'ph', 'p', 'al'])}
order = ['vm_tend_rk', 'vm_s0_rk', 'vm_s1_rk', 'vm_s2_rk', 'vm_s3_rk', 'vm_state_rk', 'vm_p_rk']
def read_wrf(fn, names):
    raw = open(fn, 'rb').read()
    h = np.frombuffer(raw[:24], '>i4').astype(int); ims, ime, kms, kme, jms, jme = h
    ni, nk, nj = ime - ims + 1, kme - kms + 1, jme - jms + 1
    off = 24; out = {}
    for n in names:
        if n.startswith('2:'):
            a = np.frombuffer(raw[off:off + 8 * ni * nj], '>f8').reshape(nj, ni); off += 8 * ni * nj
            out[n[2:]] = (a, ims, jms, True)
        else:
            a = np.frombuffer(raw[off:off + 8 * ni * nk * nj], '>f8').reshape(nj, nk, ni); off += 8 * ni * nk * nj
            out[n] = (a, ims, jms, False)
    return out
def read_c(fn, names):
    raw = open(fn, 'rb').read(); off = 0; out = {}
    for n in names:
        nx, ny, nz = np.frombuffer(raw[off:off + 12], np.int32); off += 12
        a = np.frombuffer(raw[off:off + 8 * nx * ny * nz], np.float64).reshape(ny, nz, nx); off += 8 * nx * ny * nz
        out[n] = a
    return out
for rk in (1, 2, 3):
    for base in order:
        wn, cn = specs[base]
        W = read_wrf(f'{wdir}/{base}{rk}.bin', wn); C = read_c(f'{cdir}/{base}{rk}.bin', cn)
        for n in cn:
            a, ims, jms, two = W[n]; c = C[n]
            ny, nz, nx = c.shape
            w = a[1 - jms:1 - jms + ny, 1 - ims:1 - ims + nx][:, None, :] if two else a[1 - jms:1 - jms + ny, 0:nz, 1 - ims:1 - ims + nx]
            d = np.abs(w - c); m = np.abs(w).max()
            frac = (w != c).mean()
            print(f'rk{rk} {base[3:-3]:6s} {n:8s} maxabs {d.max():9.2e} maxrel {d.max()/max(m,1e-300):9.2e} frac {frac:6.3f}')
