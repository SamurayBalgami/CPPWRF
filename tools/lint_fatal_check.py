#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Lint: every launch of a kernel that can report a WRF fatal stop is followed by fatal_check().

A kernel reports WRF's wrf_error_fatal with raise_fatal(code) (src/core/backend.hpp); the host stops the run only in
fatal_check(), and on the GPU the flag is per translation unit, so the check must follow the launch in the same
function. This script finds, in src/:
  1. the raising code: every struct / class or function whose body calls raise_fatal, and (repeatedly) every struct or
     function without a launch of its own that names one of them (a kernel functor calling a raising device function,
     a runner wrapping a raising kernel, ...; a name defined in several files only within its own file);
  2. every launch (parallel_for, column_for, column_for_work, column_for_threads, ring_for) whose arguments name raising
     code, directly or through a variable or parameter of a raising type;
and fails when such a launch is not followed by fatal_check() later in the same function (lambdas belong to the
function that defines them; a launch inside a lambda needs the check inside that lambda). It prints the raising launches it found with their checks.

  tools/lint_fatal_check.py [SRC_DIR]     exit 0: clean; 1: a launch without fatal_check (or no raising launch found)
"""
import os
import re
import sys

LAUNCH = re.compile(r"\b(parallel_for|column_for|column_for_work|column_for_threads|ring_for)\s*\(")
CHECK = re.compile(r"\bfatal_check\s*\(\s*\)")
CONTROL = {"if", "for", "while", "switch", "catch", "return", "sizeof", "decltype", "alignof", "static_assert",
           "noexcept", "requires", "defined"}
SKIP_FILES = {os.path.join("core", "backend.hpp"), os.path.join("core", "backend.cpp")}


def strip_code(text):
    """Comments, string / char literals and preprocessor lines blanked (newlines kept, so offsets keep their lines)."""
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append(re.sub(r"[^\n]", " ", text[i:j]))
            i = j
        elif c in "\"'":
            if c == "'" and i > 0 and (text[i - 1].isalnum()):  # digit separator 1'000
                out.append(" ")
                i += 1
                continue
            j = i + 1
            while j < n and text[j] != c and text[j] != "\n":
                j += 2 if text[j] == "\\" else 1
            j = min(j + 1, n)
            out.append(c + " " * (j - i - 2) + c if j - i >= 2 else " " * (j - i))
            i = j
        else:
            out.append(c)
            i += 1
    s = "".join(out)
    # preprocessor directives (with continuation lines)
    return re.sub(r"(?m)^[ \t]*#(?:[^\n]*\\\n)*[^\n]*", lambda m: re.sub(r"[^\n]", " ", m.group(0)), s)


def matching(s, i, open_ch, close_ch):
    d = 0
    for j in range(i, len(s)):
        if s[j] == open_ch:
            d += 1
        elif s[j] == close_ch:
            d -= 1
            if d == 0:
                return j
    return len(s) - 1


def classify(header):
    """('struct', name) | ('function', name) | ('namespace', '') | ('lambda', '') | ('block', '')"""
    h = " ".join(header.split())
    if re.search(r"\bnamespace\b", h) or re.match(r'^extern\s*"', h):
        return "namespace", ""
    if re.search(r"\]\s*(\([^()]*\))?\s*(mutable|constexpr|noexcept|->[^{]*|\s)*$", h) and "[" in h:
        return "lambda", ""
    h2 = re.sub(r"\btemplate\s*<[^<>]*(<[^<>]*>[^<>]*)*>", " ", h)
    m = re.search(r"\b(struct|class|union)\s+(?:CPPWRF_\w+\s+|alignas\s*\([^)]*\)\s*)*(\w+)", h2)
    if m and not re.search(r"\benum\b", h2) and "(" not in h2[: m.start()]:
        return "struct", m.group(2)
    # a function: NAME ( ... ) [qualifiers] [: init list]
    for m in re.finditer(r"([A-Za-z_][\w:]*(?:operator\s*\(\s*\)|operator\s*\S+)?)\s*\(", h2):
        name = m.group(1)
        if name.split("::")[-1] in CONTROL:
            return "block", ""
        if re.match(r"^\s*(if|for|while|switch|catch|do|else|try)\b", h2):
            return "block", ""
        if "=" in h2[: m.start()] and "operator" not in name:  # an initialiser: x = f(...) {
            return "block", ""
        return "function", name
    return "block", ""


class Entity:
    def __init__(self, kind, name, path, hstart, start, end, text):
        self.kind, self.name, self.path = kind, name, path
        self.hstart, self.start, self.end, self.text = hstart, start, end, text

    @property
    def body(self):
        return self.text[self.start:self.end]

    @property
    def full(self):
        return self.text[self.hstart:self.end]


def entities(path, text):
    """The outermost struct / function definitions of a file (members and lambdas belong to them)."""
    res, stack, last = [], [], 0  # stack: (kind, name, hstart, brace position)
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in ";}":
            if c == "}" and stack:
                kind, name, hstart, b = stack.pop()
                if kind in ("struct", "function") and not any(k in ("struct", "function") for k, *_ in stack):
                    qual = name.split("::")
                    ename = qual[-2] if kind == "function" and len(qual) > 1 else qual[-1]
                    res.append(Entity(kind, ename, path, hstart, b + 1, i, text))
            last = i + 1
        elif c == "{":
            header = text[last:i]
            kind, name = classify(header)
            stack.append((kind, name, last, i))
            last = i + 1
        i += 1
    return res


def lambda_spans(s):
    """(start, end) of the lambda bodies in s"""
    spans, stack, last = [], [], 0
    for i, c in enumerate(s):
        if c in ";}":
            if c == "}" and stack:
                kind, b = stack.pop()
                if kind == "lambda":
                    spans.append((b, i))
            last = i + 1
        elif c == "{":
            stack.append((classify(s[last:i])[0], i))
            last = i + 1
    return spans


def idents(s):
    return set(re.findall(r"\b[A-Za-z_]\w*\b", s))


def line_of(text, pos):
    return text.count("\n", 0, pos) + 1


def main():
    src = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "src"))
    ents = []
    for root, _, files in os.walk(src):
        for f in sorted(files):
            if not f.endswith((".hpp", ".cpp", ".cu", ".cuh", ".h")):
                continue
            p = os.path.join(root, f)
            if os.path.relpath(p, src) in SKIP_FILES:
                continue
            ents += entities(os.path.relpath(p, src), strip_code(open(p, encoding="utf-8", errors="replace").read()))
    for e in ents:
        e.launch = bool(LAUNCH.search(e.body))
        e.ids = idents(e.full)
    # a name defined in several files (run, compute, ...) is matched only within the file that defines it as raising
    files = {}
    for e in ents:
        files.setdefault(e.name, set()).add(e.path)
    rfiles = {}  # raising name -> the files where it is raising code

    def add(e):
        rfiles.setdefault(e.name, set()).add(e.path)

    def refs(e):
        return {n for n in e.ids if n in rfiles and (len(files.get(n, ())) <= 1 or e.path in rfiles[n])}

    for e in ents:
        if re.search(r"\braise_fatal\s*\(", e.body):
            add(e)
    origin = set(rfiles)
    changed = True
    while changed:
        changed = False
        for e in ents:
            if not e.launch and e.path not in rfiles.get(e.name, ()) and refs(e):
                add(e)
                changed = True
    bad, found = [], []
    for e in ents:
        if not e.launch:
            continue
        full, off = e.full, e.start - e.hstart
        # variables and parameters of a raising type: T var, T<...> var, auto var = T...
        raising = refs(e)
        rvars = set()
        for m in re.finditer(r"\b([A-Za-z_]\w*)\s*(?:<[^;{}()]*>)?\s*[&*]?\s*\b([A-Za-z_]\w*)\s*(?=[{;=(),])", full):
            if m.group(1) in raising:
                rvars.add(m.group(2))
        for m in re.finditer(r"\bauto\s*&?\s*([A-Za-z_]\w*)\s*=\s*([A-Za-z_]\w*)", full):
            if m.group(2) in raising:
                rvars.add(m.group(1))
        lams = lambda_spans(full)
        for m in LAUNCH.finditer(full, off):
            close = matching(full, m.end() - 1, "(", ")")
            # in a lambda the check must be in that lambda (it may be called alone)
            end = min([b for a, b in lams if a < m.start() < b] or [len(full)])
            args = idents(full[m.end():close])
            hit = sorted(args & (raising | rvars))
            if not hit:
                continue
            where = f"{e.path}:{line_of(e.text, e.hstart + m.start())}"
            chk = CHECK.search(full, close, end)
            if chk:
                found.append(f"{where}: {m.group(1)}({', '.join(hit)}) -> fatal_check() at line "
                             f"{line_of(e.text, e.hstart + chk.start())}")
            else:
                bad.append(f"{where}: {m.group(1)} of raising code ({', '.join(hit)}) in {e.name} is not followed by "
                           f"fatal_check() in the same function")
    print(f"raise_fatal in: {', '.join(sorted(origin))}")
    print(f"raising code: {len(rfiles)} structs / functions; launches of it: {len(found) + len(bad)}")
    for x in found:
        print("  ok   " + x)
    for x in bad:
        print("  FAIL " + x)
    if not found and not bad:
        print("FAIL no launch of raising code found (the lint no longer recognises the code)")
        return 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
