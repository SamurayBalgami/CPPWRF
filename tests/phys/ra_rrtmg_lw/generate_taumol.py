#!/usr/bin/env python3
"""
Translates taugb1..taugb16 from tests/phys/fortran/module_ra_rrtmg_lw.F into
src/phys/ra_rrtmg_lw/rrtmg_lw_taumol.hpp.
"""
import re

ngc = [10, 12, 16, 14, 16, 8, 12, 8, 12, 6, 8, 8, 4, 2, 2, 2]
ngs = [0, 10, 22, 38, 52, 68, 76, 88, 96, 108, 114, 122, 130, 134, 136, 138]
nspa = [1, 1, 9, 9, 9, 1, 9, 1, 9, 1, 1, 9, 9, 1, 9, 9]
nspb = [1, 1, 5, 5, 5, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1]

ngs_map = {
    "ngs1": 10, "ngs2": 22, "ngs3": 38, "ngs4": 52,
    "ngs5": 68, "ngs6": 76, "ngs7": 88, "ngs8": 96,
    "ngs9": 108, "ngs10": 114, "ngs11": 122, "ngs12": 130,
    "ngs13": 134, "ngs14": 136, "ngs15": 138
}

def join_continuation_lines(raw_lines):
    joined = []
    curr = ""
    for line in raw_lines:
        if "!" in line:
            line = line[:line.index("!")]
        s = line.strip()
        if not s:
            continue
        if curr:
            if s.startswith("&"):
                s = s[1:].strip()
            curr = curr + " " + s
        else:
            curr = s

        if curr.endswith("&"):
            curr = curr[:-1].strip()
        else:
            joined.append(curr)
            curr = ""
    if curr:
        joined.append(curr)
    return joined

def build_converter():
    with open("tests/phys/fortran/module_ra_rrtmg_lw.F", "r", errors="ignore") as f:
        lines = f.readlines()

    out = []
    out.append('#pragma once\n\n')
    out.append('#include "core/backend.hpp"\n')
    out.append('#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"\n')
    out.append('#include "phys/ra_rrtmg_lw/rrtmg_lw_setcoef.hpp"\n')
    out.append('#include <algorithm>\n#include <cmath>\n\n')
    out.append('namespace cppwrf::phys {\n\n')

    for b in range(1, 17):
        b_name = f"taugb{b}"
        start_idx = -1
        end_idx = -1
        for idx, line in enumerate(lines):
            l = line.strip().lower()
            if l.startswith(f"subroutine {b_name}"):
                start_idx = idx
            elif start_idx != -1 and l.startswith(f"end subroutine {b_name}"):
                end_idx = idx
                break

        assert start_idx != -1 and end_idx != -1, f"Could not find subroutine {b_name}"
        band_raw = lines[start_idx+1:end_idx]
        band_joined = join_continuation_lines(band_raw)
        cpp_band = translate_band(b, band_joined)
        out.append(cpp_band)
        out.append('\n')

    # Write wrapper taumol function
    out.append('''CPPWRF_HD inline void rrtmg_lw_taumol(
    int nlayers,
    const double* pavel,
    const double* const* wx,
    const double* coldry,
    int laytrop,
    const int* jp,
    const int* jt,
    const int* jt1,
    const double* colh2o,
    const double* colco2,
    const double* colo3,
    const double* coln2o,
    const double* colco,
    const double* colch4,
    const double* colo2,
    const double* colbrd,
    const double* fac00,
    const double* fac01,
    const double* fac10,
    const double* fac11,
    const double* rat_h2oco2,
    const double* rat_h2oco2_1,
    const double* rat_h2oo3,
    const double* rat_h2oo3_1,
    const double* rat_h2on2o,
    const double* rat_h2on2o_1,
    const double* rat_h2och4,
    const double* rat_h2och4_1,
    const double* rat_n2oco2,
    const double* rat_n2oco2_1,
    const double* rat_o3co2,
    const double* rat_o3co2_1,
    const double* selffac,
    const double* selffrac,
    const int* indself,
    const double* forfac,
    const double* forfrac,
    const int* indfor,
    const double* minorfrac,
    const double* scaleminor,
    const double* scaleminorn2,
    const int* indminor,
    double* const* fracs,  // [nlayers][140]
    double* const* taug,   // [nlayers][140]
    const RrtmgLwTables& tables
) {
''')
    for b in range(1, 17):
        out.append(f'  rrtmg_lw_taugb{b:02d}(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables);\n')
    out.append('}\n\n}  // namespace cppwrf::phys\n')

    with open("src/phys/ra_rrtmg_lw/rrtmg_lw_taumol.hpp", "w") as f:
        f.writelines(out)

    print("Successfully generated src/phys/ra_rrtmg_lw/rrtmg_lw_taumol.hpp!")

def translate_band(b, code_lines):
    ng = ngc[b - 1]
    ng_offset = ngs[b - 1]
    na = nspa[b - 1]
    nb = nspb[b - 1]
    b_str = f"{b:02d}"

    header = f'''CPPWRF_HD inline void rrtmg_lw_taugb{b_str}(
    int nlayers,
    const double* pavel,
    const double* const* wx,
    const double* coldry,
    int laytrop,
    const int* jp,
    const int* jt,
    const int* jt1,
    const double* colh2o,
    const double* colco2,
    const double* colo3,
    const double* coln2o,
    const double* colco,
    const double* colch4,
    const double* colo2,
    const double* colbrd,
    const double* fac00,
    const double* fac01,
    const double* fac10,
    const double* fac11,
    const double* rat_h2oco2,
    const double* rat_h2oco2_1,
    const double* rat_h2oo3,
    const double* rat_h2oo3_1,
    const double* rat_h2on2o,
    const double* rat_h2on2o_1,
    const double* rat_h2och4,
    const double* rat_h2och4_1,
    const double* rat_n2oco2,
    const double* rat_n2oco2_1,
    const double* rat_o3co2,
    const double* rat_o3co2_1,
    const double* selffac,
    const double* selffrac,
    const int* indself,
    const double* forfac,
    const double* forfrac,
    const int* indfor,
    const double* minorfrac,
    const double* scaleminor,
    const double* scaleminorn2,
    const int* indminor,
    double* const* fracs,
    double* const* taug,
    const RrtmgLwTables& tables
) {{
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b{b_str};
'''
    res = [header]

    for line in code_lines:
        l = line.strip().lower()
        if l.startswith("use ") or l.startswith("implicit ") or l.startswith("save"):
            continue
        if l.startswith("integer") or l.startswith("real"):
            m = re.match(r"(?:integer|real)\s*(?:\([^)]*\))?\s*::\s*(.*)", line, re.I)
            if m:
                vars_str = m.group(1).strip()
                t = "int" if l.startswith("integer") else "double"
                res.append(f"  {t} {vars_str};\n")
            continue

        stmt = translate_statement(line, b, ng, ng_offset, na, nb)
        res.append(stmt + "\n")

    res.append("}\n")
    return "".join(res)

def translate_statement(line, b, ng, ng_offset, na, nb):
    l = line.strip().lower()

    # Do loops
    m_do_lay1 = re.match(r"do\s+lay\s*=\s*1\s*,\s*laytrop", l)
    if m_do_lay1:
        return "  for (int lay = 0; lay < laytrop; ++lay) {"

    m_do_lay2 = re.match(r"do\s+lay\s*=\s*laytrop\s*\+\s*1\s*,\s*nlayers", l)
    if m_do_lay2:
        return "  for (int lay = laytrop; lay < nlayers; ++lay) {"

    m_do_ig = re.match(r"do\s+ig\s*=\s*1\s*,\s*ng\d+", l)
    if m_do_ig:
        return f"  for (int ig = 1; ig <= {ng}; ++ig) {{"

    if l == "enddo" or l == "end do":
        return "  }"

    # If statements
    if l.startswith("if (") and l.endswith(") then"):
        cond = line.strip()[4:-6].strip()
        cond_cpp = translate_expr(cond, b, ng, ng_offset, na, nb)
        return f"  if ({cond_cpp}) {{"
    elif l.startswith("else if (") and l.endswith(") then"):
        cond = line.strip()[9:-6].strip()
        cond_cpp = translate_expr(cond, b, ng, ng_offset, na, nb)
        return f"  }} else if ({cond_cpp}) {{"
    elif l == "else":
        return "  } else {"
    elif l == "endif" or l == "end if":
        return "  }"
    elif l.startswith("if (") and "then" not in l:
        m_if = re.match(r"if\s*\((.*)\)\s*(.*)", line.strip(), re.I)
        if m_if:
            cond = translate_expr(m_if.group(1), b, ng, ng_offset, na, nb)
            stmt = translate_expr(m_if.group(2), b, ng, ng_offset, na, nb)
            return f"  if ({cond}) {stmt};"

    expr = translate_expr(line.strip(), b, ng, ng_offset, na, nb)
    if not expr.endswith(";"):
        expr += ";"
    return "  " + expr

def translate_expr(text, b, ng, ng_offset, na, nb):
    s = text

    # Remove Fortran real literal suffixes: _rb, _im, _8
    s = re.sub(r"_rb\b", "", s, flags=re.I)
    s = re.sub(r"_im\b", "", s, flags=re.I)
    s = re.sub(r"_8\b", "", s, flags=re.I)

    # Substitute ngsN with its numerical value
    for k, v in ngs_map.items():
        s = re.sub(rf"\b{k}\b", str(v), s, flags=re.I)

    # General powers: (expr)**power or var**power
    s = re.sub(r"\bp\s*\*\*\s*4\b", "(p * p * p * p)", s)
    s = re.sub(r"\b(\w+)\s*\*\*\s*2\b", r"(\1 * \1)", s)
    s = re.sub(r"(\([^)]+\)|\w+)\s*\*\*\s*([0-9.]+)", r"std::pow(\1, \2)", s)

    # Logical operators
    s = re.sub(r"\.eq\.", " == ", s, flags=re.I)
    s = re.sub(r"\.ne\.", " != ", s, flags=re.I)
    s = re.sub(r"\.lt\.", " < ", s, flags=re.I)
    s = re.sub(r"\.le\.", " <= ", s, flags=re.I)
    s = re.sub(r"\.gt\.", " > ", s, flags=re.I)
    s = re.sub(r"\.ge\.", " >= ", s, flags=re.I)
    s = re.sub(r"\.and\.", " && ", s, flags=re.I)
    s = re.sub(r"\.or\.", " || ", s, flags=re.I)

    # Intrinsics
    s = re.sub(r"\bexp\(", "std::exp(", s)
    s = re.sub(r"\bmin\(", "std::min(", s)
    s = re.sub(r"\bmax\(", "std::max(", s)
    s = re.sub(r"\bint\(", "static_cast<int>(", s)
    s = re.sub(r"\bfloat\(", "static_cast<double>(", s)
    s = re.sub(r"\bmod\(([^,]+),\s*1\.0\)", r"std::fmod(\1, 1.0)", s)

    # 1D array parameters indexed by lay
    arrays_lay = [
        "pavel", "colh2o", "colco2", "colo3", "coln2o", "colco", "colch4", "colo2", "colbrd", "coldry",
        "fac00", "fac01", "fac10", "fac11", "jp", "jt", "jt1",
        "rat_h2oco2", "rat_h2oco2_1", "rat_h2oo3", "rat_h2oo3_1", "rat_h2on2o", "rat_h2on2o_1", "rat_h2och4", "rat_h2och4_1",
        "rat_n2oco2", "rat_n2oco2_1", "rat_o3co2", "rat_o3co2_1",
        "selffac", "selffrac", "indself", "forfac", "forfrac", "indfor",
        "minorfrac", "scaleminor", "scaleminorn2", "indminor"
    ]
    for arr in arrays_lay:
        s = re.sub(rf"\b{arr}\(lay\)", f"{arr}[lay]", s)

    # wx(ix, lay)
    s = re.sub(r"\bwx\(([^,]+),\s*lay\)", r"wx[\1 - 1][lay]", s)

    # nspa and nspb references
    s = re.sub(rf"\bnspa\({b}\)", str(na), s)
    s = re.sub(rf"\bnspb\({b}\)", str(nb), s)

    # chi_mls(ispec, ilev)
    s = re.sub(r"\bchi_mls\(([^,]+),\s*([^)]+)\)", r"get_chi_mls(tables, \1, \2)", s)

    # taug and fracs output arrays: taug(lay, X) -> taug[lay][X - 1]
    # In band 1, it might be taug(lay, ig) -> taug[lay][ig - 1]
    s = re.sub(r"\btaug\(lay,\s*([^)]+)\)", r"taug[lay][\1 - 1]", s)
    s = re.sub(r"\bfracs\(lay,\s*([^)]+)\)", r"fracs[lay][\1 - 1]", s)

    # Table lookups for band b
    na_65 = na * 65
    nb_235 = nb * 235

    s = re.sub(r"\babsa\(([^,]+),\s*ig\)", rf"b.ka[(ig - 1) * {na_65} + (\1 - 1)]", s)
    s = re.sub(r"\babsb\(([^,]+),\s*ig\)", rf"b.kb[(ig - 1) * {nb_235} + (\1 - 1)]", s)
    s = re.sub(r"\bselfref\(([^,]+),\s*ig\)", r"b.selfref[(ig - 1) * 10 + (\1 - 1)]", s)
    s = re.sub(r"\bforref\(([^,]+),\s*ig\)", r"b.forref[(ig - 1) * 4 + (\1 - 1)]", s)

    # 1D/2D fracrefa / fracrefb
    if na == 1:
        s = re.sub(r"\bfracrefa\(ig\)", r"b.fracrefa[ig - 1]", s)
    else:
        s = re.sub(r"\bfracrefa\(ig,\s*([^)]+)\)", rf"b.fracrefa[(\1 - 1) * {ng} + (ig - 1)]", s)

    if nb == 1 or nb == 0:
        s = re.sub(r"\bfracrefb\(ig\)", r"b.fracrefb[ig - 1]", s)
    else:
        s = re.sub(r"\bfracrefb\(ig,\s*([^)]+)\)", rf"b.fracrefb[(\1 - 1) * {ng} + (ig - 1)]", s)

    # Minor gas arrays per band
    if b == 1:
        s = re.sub(r"\bka_mn2\(([^,]+),\s*ig\)", r"b.ka_mn2[(ig - 1) * 19 + (\1 - 1)]", s)
        s = re.sub(r"\bkb_mn2\(([^,]+),\s*ig\)", r"b.kb_mn2[(ig - 1) * 19 + (\1 - 1)]", s)
    elif b == 3:
        s = re.sub(r"\bka_mn2o\(([^,]+),\s*([^,]+),\s*ig\)", r"b.ka_mn2o[(ig - 1) * (19 * 9) + (\2 - 1) * 9 + (\1 - 1)]", s)
        s = re.sub(r"\bkb_mn2o\(([^,]+),\s*([^,]+),\s*ig\)", r"b.kb_mn2o[(ig - 1) * (19 * 5) + (\2 - 1) * 5 + (\1 - 1)]", s)
    elif b == 5:
        s = re.sub(r"\bka_mo3\(([^,]+),\s*([^,]+),\s*ig\)", r"b.ka_mo3[(ig - 1) * (19 * 9) + (\2 - 1) * 9 + (\1 - 1)]", s)
        s = re.sub(r"\bccl4\(ig\)", r"b.ccl4[ig - 1]", s)
    elif b == 6:
        s = re.sub(r"\bka_mco2\(([^,]+),\s*ig\)", r"b.ka_mco2[(ig - 1) * 19 + (\1 - 1)]", s)
        s = re.sub(r"\bcfc11adj\(ig\)", r"b.cfc11adj[ig - 1]", s)
        s = re.sub(r"\bcfc12\(ig\)", r"b.cfc12[ig - 1]", s)
    elif b == 7:
        s = re.sub(r"\bka_mco2\(([^,]+),\s*([^,]+),\s*ig\)", r"b.ka_mco2[(ig - 1) * (19 * 9) + (\2 - 1) * 9 + (\1 - 1)]", s)
        s = re.sub(r"\bkb_mco2\(([^,]+),\s*ig\)", r"b.kb_mco2[(ig - 1) * 19 + (\1 - 1)]", s)
    elif b == 8:
        s = re.sub(r"\bka_mco2\(([^,]+),\s*ig\)", r"b.ka_mco2[(ig - 1) * 19 + (\1 - 1)]", s)
        s = re.sub(r"\bkb_mco2\(([^,]+),\s*ig\)", r"b.kb_mco2[(ig - 1) * 19 + (\1 - 1)]", s)
        s = re.sub(r"\bka_mn2o\(([^,]+),\s*ig\)", r"b.ka_mn2o[(ig - 1) * 19 + (\1 - 1)]", s)
        s = re.sub(r"\bkb_mn2o\(([^,]+),\s*ig\)", r"b.kb_mn2o[(ig - 1) * 19 + (\1 - 1)]", s)
        s = re.sub(r"\bka_mo3\(([^,]+),\s*ig\)", r"b.ka_mo3[(ig - 1) * 19 + (\1 - 1)]", s)
        s = re.sub(r"\bcfc12\(ig\)", r"b.cfc12[ig - 1]", s)
        s = re.sub(r"\bcfc22adj\(ig\)", r"b.cfc22adj[ig - 1]", s)
    elif b == 9:
        s = re.sub(r"\bka_mn2o\(([^,]+),\s*([^,]+),\s*ig\)", r"b.ka_mn2o[(ig - 1) * (19 * 9) + (\2 - 1) * 9 + (\1 - 1)]", s)
        s = re.sub(r"\bkb_mn2o\(([^,]+),\s*ig\)", r"b.kb_mn2o[(ig - 1) * 19 + (\1 - 1)]", s)
    elif b == 11:
        s = re.sub(r"\bka_mo2\(([^,]+),\s*ig\)", r"b.ka_mo2[(ig - 1) * 19 + (\1 - 1)]", s)
        s = re.sub(r"\bkb_mo2\(([^,]+),\s*ig\)", r"b.kb_mo2[(ig - 1) * 19 + (\1 - 1)]", s)
    elif b == 13:
        s = re.sub(r"\bka_mco2\(([^,]+),\s*([^,]+),\s*ig\)", r"b.ka_mco2[(ig - 1) * (19 * 9) + (\2 - 1) * 9 + (\1 - 1)]", s)
        s = re.sub(r"\bka_mco\(([^,]+),\s*([^,]+),\s*ig\)", r"b.ka_mco[(ig - 1) * (19 * 9) + (\2 - 1) * 9 + (\1 - 1)]", s)
        s = re.sub(r"\bkb_mo3\(([^,]+),\s*ig\)", r"b.kb_mo3[(ig - 1) * 19 + (\1 - 1)]", s)
    elif b == 15:
        s = re.sub(r"\bka_mn2\(([^,]+),\s*([^,]+),\s*ig\)", r"b.ka_mn2[(ig - 1) * (19 * 9) + (\2 - 1) * 9 + (\1 - 1)]", s)

    return s

if __name__ == "__main__":
    build_converter()
