"""Tolerance-0 comparison of model output, shared by the harness (wrf_suite.py, compare_real.py,
tests/init_real/compare_real.py, tools/gcp/matrix_variant.py).

Two values are identical when they are the same number: NaN equals NaN, NaN differs from any number, -0 differs
from +0 (WRF writes -0, e.g. in W, C1H, C1F), infinities equal themselves. netCDF variables are read as stored
(no masking of fill values, no scaling)."""
import numpy as np


def raw(ds):
    """The netCDF4 dataset `ds` with masking and scaling off (values as stored)."""
    ds.set_auto_maskandscale(False)
    return ds


def differ(x, y):
    """Mask of the points where x and y (same shape) are not the same value."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    nx, ny = np.isnan(x), np.isnan(y)
    d = (x != y) & ~(nx & ny)
    return d | ((x == 0) & (y == 0) & (np.signbit(x) != np.signbit(y)))


def describe(x, y):
    """None when x and y are identical, else 'n/size max|d| (first: a vs b at index)' with the NaN-vs-number and
    signed-zero point counts."""
    x, y = np.atleast_1d(np.asarray(x, dtype=np.float64)), np.atleast_1d(np.asarray(y, dtype=np.float64))
    if x.shape != y.shape:
        return f"shape {x.shape} vs {y.shape}"
    d = differ(x, y)
    n = int(d.sum())
    if not n:
        return None
    i = int(np.flatnonzero(d)[0])
    nan = int((np.isnan(x) != np.isnan(y)).sum())
    sz = int(((x == 0) & (y == 0) & (np.signbit(x) != np.signbit(y))).sum())
    with np.errstate(invalid="ignore"):
        m = np.abs(x - y)[d & ~np.isnan(x) & ~np.isnan(y)]
    return (f"{n}/{x.size} max|d| {m.max() if m.size else 0.0:.3e} (first {float(x.flat[i])!r} vs {float(y.flat[i])!r} at {i})"
            + (f" nan-vs-number {nan}" if nan else "") + (f" signed-zero {sz}" if sz else ""))


def text(v, t=None):
    """Bytes of a character variable (all of it, or frame t)."""
    a = np.asarray(v[:] if t is None else v[t])
    return b"".join(a.ravel().tolist()) if a.dtype.kind == "S" else a.tobytes()


def metadata(a, b, ignore=("TITLE",), only=None):
    """Differences of two netCDF datasets other than values: variable list and order, dimensions, each variable's
    dimensions, type and attributes, the global attributes (names in `ignore` aside) and their order.
    only=[variables of b]: b holds these variables of a full file (their order in a, the dimensions they use)."""
    out = []
    va, vb = list(a.variables), list(b.variables)
    if only is not None:
        va = [x for x in va if x in only]
    if va != vb:
        out.append(f"variables: missing {[x for x in vb if x not in va]} extra {[x for x in va if x not in vb]}"
                   + ("" if set(va) != set(vb) else " (order)"))
    dims = lambda f: [(d, "U" if f.dimensions[d].isunlimited() else len(f.dimensions[d])) for d in f.dimensions]  # noqa: E731
    da, db = dims(a), dims(b)
    if only is not None:
        used = {d for n in vb for d in b.variables[n].dimensions}
        da, db = [d for d in da if d[0] in used], [d for d in db if d[0] in used]
    if da != db:
        out.append(f"dimensions {da} vs {db}")
    for n in vb:
        if n not in va:
            continue
        x, y = a.variables[n], b.variables[n]
        if x.dimensions != y.dimensions or x.dtype != y.dtype:
            out.append(f"{n}: dims/type {x.dimensions} {x.dtype} vs {y.dimensions} {y.dtype}")
        ax = [(k, repr(x.getncattr(k))) for k in x.ncattrs()]
        ay = [(k, repr(y.getncattr(k))) for k in y.ncattrs()]
        if ax != ay:
            out.append(f"{n}: attributes {ax} vs {ay}")
    ga = [(k, repr(a.getncattr(k)), type(a.getncattr(k)).__name__) for k in a.ncattrs() if k not in ignore]
    gb = [(k, repr(b.getncattr(k)), type(b.getncattr(k)).__name__) for k in b.ncattrs() if k not in ignore]
    if ga != gb:
        ma, mb = {k: v for k, *v in ga}, {k: v for k, *v in gb}
        dif = {k: f"{ma.get(k)} vs {mb.get(k)}" for k in sorted(set(ma) | set(mb)) if ma.get(k) != mb.get(k)}
        out.append(f"global attributes {dif}" if dif else "global attribute order differs")
    return out


def frames(files):
    """{Times string: (dataset, frame index)} over the files of one stream (raw datasets)."""
    import netCDF4
    m = {}
    for f in sorted(files):
        ds = raw(netCDF4.Dataset(f))
        for t in range(len(ds["Times"])):
            m[text(ds["Times"], t).decode()] = (ds, t)
    return m


def compare_stream(cfiles, wfiles, subset=False, until=None, global_atts=True):
    """Differences of CPPWRF's files of one stream and domain (wrfout_d0N or wrfrst_d0N, any frames per file) from
    WRF's: ({label: difference}, number of common frames). Every variable of every common frame (by Times) at
    tolerance 0, types, and the metadata of the files holding the first common frame; times or variables only one
    model wrote are differences ('(times)', '(wrf-only-vars)', '(cppwrf-only-vars)', '(metadata)').
    subset=True: WRF's files hold a part of the output (tools/colab/stage_data.py keeps some frames and variables):
    CPPWRF's other variables and frames are no difference, WRF's frames up to `until` (Times string) must all be
    there and the metadata covers WRF's variables. global_atts=False: the global attributes are not compared."""
    A, B = frames(cfiles), frames(wfiles)
    bad, nframes = {}, {}
    common = sorted(set(A) & set(B))
    only_c = [] if subset else sorted(set(A) - set(B))
    only_w = sorted(t for t in set(B) - set(A) if until is None or t <= until)
    if only_c or only_w:  # times written by one model only
        bad["(times)"] = f"only-cppwrf {only_c} only-wrf {only_w}"
    wrf_only, cppwrf_only = set(), set()
    for k in common:
        (a, ta), (b, tb) = A[k], B[k]
        cppwrf_only |= {v for v in a.variables if v not in b.variables}
        wrf_only |= {v for v in b.variables if v not in a.variables}
        for v in b.variables:
            if v not in a.variables or "Time" not in b[v].dimensions:
                continue
            if b[v].dtype.kind == "S":
                if v != "Times" and text(a[v], ta) != text(b[v], tb):
                    bad.setdefault(v, f"{k} text differs")
                continue
            x, y = np.asarray(a[v][ta]), np.asarray(b[v][tb])
            if x.shape != y.shape:
                if v not in bad:  # WRF's values per leading index (e.g. soil layers)
                    yl = y.reshape(y.shape[0], -1) if y.ndim > 1 else y.reshape(-1, 1)
                    bad[v] = (f"shape cppwrf {x.shape} wrf {y.shape} wrf min/max per level "
                              + " ".join(f"{r.min():.4g}/{r.max():.4g}" for r in yl))
                continue
            if a[v].dtype != b[v].dtype:
                bad.setdefault(v, f"type cppwrf {a[v].dtype} wrf {b[v].dtype}")
            dif = describe(x, y)
            if dif:
                nframes[v] = nframes.get(v, 0) + 1
                bad.setdefault(v, f"{k} {dif}")
    if wrf_only:  # variables WRF writes and CPPWRF does not
        bad["(wrf-only-vars)"] = sorted(wrf_only)
    if cppwrf_only and not subset:  # variables CPPWRF writes and WRF does not
        bad["(cppwrf-only-vars)"] = sorted(cppwrf_only)
    for v in nframes:
        if nframes[v] > 1:
            bad[v] += f" [{nframes[v]} frames]"
    if common:  # dimensions, types, variable and global attributes of the files holding the first common frame
        a, b = A[common[0]][0], B[common[0]][0]
        ignore = ("TITLE",) if global_atts else set(a.ncattrs()) | set(b.ncattrs())
        meta = metadata(a, b, ignore, only=list(b.variables) if subset else None)
        if meta:
            bad["(metadata)"] = meta
    if not common:
        bad["(frames)"] = "no common frame"
    return bad, len(common)
