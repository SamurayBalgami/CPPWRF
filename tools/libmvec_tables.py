#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Tables of glibc 2.41 libmvec's SSE4 vector log / log10 with Intel's rcpps as a C++ header.

  python3 tools/libmvec_tables.py GLIBC-2.41-SOURCE > src/core/libmvec_tables.hpp     (on an Intel x86-64 CPU)
  python3 tools/libmvec_tables.py --check src/core/libmvec_tables.hpp    (any x86-64 CPU: exit 1 if its rcpps differs)
  python3 tools/libmvec_tables.py --hash        (Intel: the hashes of libmvec's results that tests/unit/test_libmvec.cpp
                                                 expects)

WRF's gfortran build calls _ZGVbN2v_log / _ZGVbN2v_log10 in vectorised loops; on any CPU with SSE4.1 glibc
selects their SSE4 variants (sysdeps/x86_64/fpu/multiarch/svml_d_log2_core_sse4.S, svml_d_log102_core_sse4.S).
Both take an approximate reciprocal of the mantissa with `rcpps` and round it to an integer: the only operation of
the algorithms whose result the IEEE standard does not define. Intel and AMD implement rcpps differently, so the
vector log differs between them; WRF's references are made on Intel. src/core/libmvec_x86.hpp re-implements
the two variants with that integer taken from the table written here (Intel's rcpps, dumped on this machine):

  m = mantissa of x with the exponent of 2^-10 (double), f = (float)m in [2^-10, 2^-9]  (cvtpd2ps)
  k = roundeven((double)rcpps(f))                                       (cvtps2pd + roundpd $0 / + 1.5*2^52)

The dump covers every input the algorithms give rcpps: the 2^23 floats with f's biased exponent 117 and
f = 2^-9 (cvtpd2ps rounding up). It checks that rcpps depends only on the top 11 bits of f's mantissa (so
t = (bits(f) >> 12) - (117 << 11), 0 ... 2048, indexes every input class, 2048 being f = 2^-9), that it scales
exactly with the exponent (every biased exponent 1 ... 252, all mantissas: the table is independent of the
exponent the algorithm chose) and that k(t) falls monotonically from 1024 to 512 taking every integer value;
the header stores the 512 steps (the first t where k drops below 1024 - n), the minimal exact form. The glibc
data (the -log(Rcp) tables, polynomial coefficients and constants) are copied bit for bit from the source.
The rcpps loop is machine code in an executable page (no compiler needed):
  loop: movups (%rdi),%xmm0; rcpps %xmm0,%xmm0; movups %xmm0,(%rdi); add $16,%rdi; sub $4,%rsi; ja loop; ret
"""
import ctypes
import mmap
import re
import struct
import sys

import numpy as np

EXP = 117  # biased exponent of f = (float)m, m in [2^-10, 2^-9)


def cpu_model():
    vendor = model = "?"
    for line in open("/proc/cpuinfo"):
        if line.startswith("vendor_id") and vendor == "?":
            vendor = line.split(":", 1)[1].strip()
        if line.startswith("model name") and model == "?":
            model = line.split(":", 1)[1].strip()
    return vendor, model


class Rcpps:
    CODE = bytes.fromhex("0f1007" "0f53c0" "0f1107" "4883c710" "4883ee04" "77ed" "c3")

    def __init__(self):
        self.buf = mmap.mmap(-1, 4096, prot=mmap.PROT_READ | mmap.PROT_WRITE | mmap.PROT_EXEC)
        self.buf.write(self.CODE)
        addr = ctypes.addressof(ctypes.c_char.from_buffer(self.buf))
        self.fn = ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_size_t)(addr)

    def __call__(self, bits):
        """rcpps of float32 values given by their bit patterns (uint32) -> result bit patterns."""
        n = len(bits)
        b = np.zeros((n + 3) // 4 * 4, dtype=np.uint32)
        b[:n] = bits
        b[n:] = 0x3F800000
        self.fn(b.ctypes.data, b.size)
        return b[:n]


def dump_k(rcpps):
    mant = np.arange(1 << 23, dtype=np.uint32)
    r = rcpps((EXP << 23) | mant)
    groups = r.reshape(2048, 4096)
    if not np.all(groups == groups[:, :1]):
        sys.exit("rcpps depends on more than the top 11 mantissa bits")
    r_mant = r & 0x7FFFFF
    r_exp = (r >> 23).astype(np.int64)
    for e in range(1, 253):
        re_ = rcpps((np.uint32(e) << 23) | mant)
        if not (np.array_equal(re_ & 0x7FFFFF, r_mant) and np.array_equal((re_ >> 23).astype(np.int64),
                                                                          r_exp + (EXP - e))):
            sys.exit(f"rcpps does not scale with the exponent at biased exponent {e}")
    top = rcpps(np.array([(EXP + 1) << 23], dtype=np.uint32))  # f = 2^-9
    rt = np.concatenate([groups[:, 0], top]).view(np.float32).astype(np.float64)
    k = np.rint(rt)  # round half to even, as roundpd $0 and the 1.5*2^52 addition
    if k[0] != 1024 or k[-1] != 512 or np.any(np.diff(k) > 0) or len(np.unique(k)) != 513:
        sys.exit("k(t) is not monotone from 1024 to 512 with every value")
    steps = [int(np.argmax(k <= 1024 - n)) for n in range(1, 513)]
    expand = np.array([1024 - sum(1 for s in steps if s <= t) for t in range(2049)])
    assert np.array_equal(expand, k)
    return steps


def check(header):
    """This CPU's roundeven(rcpps(f)) for every input of the algorithms against the header's table (any vendor)."""
    text = open(header).read()
    body = text[text.index("kRcpSteps[512] = {"):]
    steps = [int(v) for v in re.findall(r"\d+", body[body.index("{") + 1:body.index("}")])]
    assert len(steps) == 512
    table = np.array([1024 - sum(1 for s in steps if s <= t) for t in range(2049)], dtype=np.float64)
    rcpps = Rcpps()
    mant = np.arange(1 << 23, dtype=np.uint32)
    f = np.concatenate([(EXP << 23) | mant, np.array([(EXP + 1) << 23], dtype=np.uint32)])
    k = np.rint(rcpps(f).view(np.float32).astype(np.float64))
    differ = int(np.sum(k != table[(f >> 12) - (EXP << 11)]))
    print(f"{cpu_model()[1]}: roundeven(rcpps(f)) differs from {header} for {differ} of {len(f)} inputs f")
    return differ


def test_inputs():
    """The input set of tests/unit/test_libmvec.cpp (bit patterns, in its order)."""
    M = (1 << 64) - 1
    st = [0]

    def sm():  # splitmix64
        st[0] = (st[0] + 0x9E3779B97F4A7C15) & M
        z = st[0]
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M
        return z ^ (z >> 31)
    x = [0x0, 0x8000000000000000, 0x1, 0x000FFFFFFFFFFFFF, 0x8000000000000001, 0x0010000000000000, 0x0010000000000001,
         0x7FEFFFFFFFFFFFFF, 0x7FF0000000000000, 0xFFF0000000000000, 0x7FF8000000000000, 0xFFF8000000000000,
         0x7FF0000000000001, 0x7FF4000000000000, 0xBFF0000000000000, 0x3FF0000000000000, 0x3FEFFFFFFFFFFFFF,
         0x3FF0000000000001, 0x4024000000000000]
    low_mask, tie = (1 << 41) - 1, (0xFFF << 29) | (1 << 28)
    for e in (1, 2, 512, 1000, 1021, 1022, 1023, 1024, 1025, 1100, 1536, 2045, 2046):
        for t in range(2048):
            for low in (0, low_mask, tie, tie - 1, sm() & low_mask):
                x.append(e << 52 | t << 41 | low)
    x += [sm() for _ in range(65536)]
    x += [0x0010000000000000 + sm() % 0x7FE0000000000000 for _ in range(65536)]
    x += [struct.unpack("<Q", struct.pack("<d", 0.5 + 1.5 * float(sm() >> 11) * 2.0 ** -53))[0] for _ in range(65536)]
    return x


def test_hashes():
    """FNV-1a (64-bit words) of glibc libmvec's log and log10 over test_inputs() on this (Intel) CPU."""
    mv = ctypes.CDLL("libmvec.so.1")
    out = []
    for name in ("_ZGVbN2v_log", "_ZGVbN2v_log10"):
        fn = getattr(mv, name)
        fn.restype = ctypes.c_double
        fn.argtypes = [ctypes.c_double]
        h = 0xCBF29CE484222325
        for b in test_inputs():
            y = struct.unpack("<Q", struct.pack("<d", fn(struct.unpack("<d", struct.pack("<Q", b))[0])))[0]
            h = ((h ^ y) * 0x100000001B3) & ((1 << 64) - 1)
        out.append(h)
    return out


def quads(text, start, n):
    i = text.index(start)
    q = re.findall(r"\.quad\s+(0x[0-9a-fA-F]+)", text[i:])[:n]
    assert len(q) == n
    return [int(x, 16) for x in q]


def log_data(glibc):
    text = open(f"{glibc}/sysdeps/x86_64/fpu/svml_d_log_data.S").read()
    la = quads(text, "!= _Log_LA_table", 513)
    consts = {}
    for name in ("poly_coeff_1", "poly_coeff_2", "poly_coeff_3", "poly_coeff_4", "ExpMask", "Two10", "MinNorm",
                 "MaxNorm", "HalfMask", "One", "Threshold", "Bias", "Bias1", "L2"):
        consts[name] = int(re.search(rf"double_vector _{name} (0x[0-9a-fA-F]+)", text).group(1), 16)
    return la, consts


def log10_data(glibc):
    path = f"{glibc}/sysdeps/x86_64/fpu/multiarch/svml_d_log102_core_sse4.S"
    text = open(path).read()
    data = text[text.index("__svml_dlog10_data_internal:"):]
    la = quads(data, "/* Log_LA_table */", 513)
    poly = quads(data, "/* poly_coeff[5] */", 5)  # one .quad line (both lanes) per coefficient: the first lane
    consts = {f"poly_coeff_{i}": v for i, v in enumerate(poly)}  # +0 (coeff5) ... +64 (coeff1)
    for name in ("ExpMask", "Two10", "MinNorm", "MaxNorm", "HalfMask", "One", "Threshold", "Bias", "Bias1", "L2"):
        consts[name] = quads(data, f"/* {name} */", 1)[0]
    flt12 = re.search(r"\.FLT_12:\s*\.long\s+0x([0-9a-fA-F]+),\s*0x([0-9a-fA-F]+)", text)
    consts["FLT_12"] = (int(flt12.group(2), 16) << 32) | int(flt12.group(1), 16)
    return la, consts


def table(name, ctype, values, fmt, per_line):
    out = [f"inline constexpr {ctype} {name}[{len(values)}] = {{"]
    for i in range(0, len(values), per_line):
        out.append("    " + ", ".join(fmt(v) for v in values[i:i + per_line]) + ",")
    out.append("};")
    return out


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--hash":
        if cpu_model()[0] != "GenuineIntel":
            sys.exit("not an Intel CPU")
        hl, hl10 = test_hashes()
        print(f"constexpr std::uint64_t kHashLog = 0x{hl:016x}ull, kHashLog10 = 0x{hl10:016x}ull;")
        return
    if len(sys.argv) == 3 and sys.argv[1] == "--check":
        sys.exit(1 if check(sys.argv[2]) else 0)
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    glibc = sys.argv[1]
    vendor, model = cpu_model()
    if vendor != "GenuineIntel":
        sys.exit(f"not an Intel CPU ({vendor}): the table must be Intel's rcpps")
    steps = dump_k(Rcpps())
    la, lc = log_data(glibc)
    la10, l10c = log10_data(glibc)
    if l10c["FLT_12"] != 0x4338000000000000:
        sys.exit("log10: .FLT_12 is not 1.5*2^52 (rounding to an integer as roundpd $0)")
    hx = lambda v: f"0x{v:016x}ull"
    out = [
        "// SPDX-License-Identifier: GPL-3.0-or-later",
        "#pragma once",
        "// GENERATED by tools/libmvec_tables.py from glibc 2.41's sources and Intel's rcpps; do not edit.",
        f"// rcpps dumped on: {model}",
        "// glibc: sysdeps/x86_64/fpu/svml_d_log_data.S (__svml_dlog_data), multiarch/svml_d_log102_core_sse4.S",
        "// (__svml_dlog10_data_internal). Values as IEEE bit patterns.",
        "#include <cstdint>",
        "",
        "namespace cppwrf::mvec_tables {",
        "",
        "// k(t) = roundeven(rcpps(f)) for t = (bits(f) >> 12) - (117 << 11), t = 0 ... 2048 (Intel): k(t) = 1024 -",
        "// (the number of entries <= t)",
        *table("kRcpSteps", "std::uint16_t", steps, str, 16),
        "",
        "// _ZGVbN2v_log_sse4: _Log_LA_table (-log(mRcp), 9-bit index k - 512) and the constants",
        *table("kLogLA", "std::uint64_t", la, hx, 4),
    ]
    for name, v in lc.items():
        out.append(f"inline constexpr std::uint64_t kLog_{name} = {hx(v)};")
    out += ["", "// _ZGVbN2v_log10_sse4: Log_LA_table and the constants (poly_coeff_0 ... 4 at poly_coeff + 0 ... 64)",
            *table("kLog10LA", "std::uint64_t", la10, hx, 4)]
    for name, v in l10c.items():
        if name != "FLT_12":
            out.append(f"inline constexpr std::uint64_t kLog10_{name} = {hx(v)};")
    out += ["", "}  // namespace cppwrf::mvec_tables", ""]
    sys.stdout.write("\n".join(out))


if __name__ == "__main__":
    main()
