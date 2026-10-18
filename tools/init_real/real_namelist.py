#!/usr/bin/env python3
"""Copy a wrf.exe namelist.input for real.exe: real.exe needs the end time at least one
interval_seconds after the start time, so a shorter run gets end = start + interval_seconds."""
import datetime as dt
import re
import sys

src, dst = sys.argv[1], sys.argv[2]
text = open(src).read()
KEYS = ['year', 'month', 'day', 'hour', 'minute', 'second']


def values(key):
    m = re.search(r'^\s*%s\s*=\s*([^\n]*)$' % key, text, re.M | re.I)
    return [int(x) for x in re.findall(r'-?\d+', m.group(1))] if m else None


def stamp(prefix, i):
    v = [values(prefix + k) for k in KEYS]
    return dt.datetime(*[(x[i] if x and len(x) > i else (x[0] if x else 0)) for x in v])


interval = (values('interval_seconds') or [21600])[0]
ndom = (values('max_dom') or [1])[0]
for i in range(ndom):
    s, e = stamp('start_', i), stamp('end_', i)
    if (e - s).total_seconds() < interval:
        e = s + dt.timedelta(seconds=interval)
        for k in KEYS:
            def sub(m, k=k, e=e):
                vals = re.findall(r'-?\d+', m.group(2))
                vals[i] = str(getattr(e, k))
                return m.group(1) + ', '.join(vals) + ','
            text = re.sub(r'^(\s*end_%s\s*=\s*)([^\n]*)$' % k, sub, text, flags=re.M | re.I)
        text = re.sub(r'^(\s*run_\w+\s*=\s*)[^\n]*\n', '', text, flags=re.M | re.I)
open(dst, 'w').write(text)
