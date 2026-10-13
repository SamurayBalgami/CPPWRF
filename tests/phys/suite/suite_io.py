"""Reader for the suite record files (suite_input.bin, suite_ref.bin, suite_cpp.bin)."""
import struct

import numpy as np


def read_records(path):
    recs = {}
    with open(path, "rb") as f:
        while True:
            hdr = f.read(32)
            if len(hdr) < 32:
                break
            name = hdr.decode().strip()
            nd, n1, n2, n3 = struct.unpack("<4i", f.read(16))
            a = np.frombuffer(f.read(8 * n1 * n2 * n3), dtype="<f8")
            recs[name] = a.reshape((n1, n2, n3), order="F")[(slice(None),) * nd + (0,) * (3 - nd)]
    return recs
