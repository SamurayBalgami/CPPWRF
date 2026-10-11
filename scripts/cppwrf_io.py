"""Python reader for CPPWRF outputs.

Supported functions:
- read_run(out_dir) -> dict
- read_diag(out_dir) -> dict[str, np.ndarray]
- list_snapshots(out_dir) -> list[int]
- read_snapshot(out_dir, step) -> dict[str, np.ndarray]
"""

import csv
import json
import os
import re
from typing import Any, Union
import numpy as np

PathLike = Union[str, os.PathLike]


def read_run(out_dir: PathLike) -> dict[str, Any]:
    """Read run.json configuration / metadata from out_dir."""
    path = os.path.join(out_dir, "run.json")
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def read_diag(out_dir: PathLike) -> dict[str, np.ndarray]:
    """Read diag.csv from out_dir into a mapping from column name to 1D numpy array."""
    path = os.path.join(out_dir, "diag.csv")
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            return {}
        cols = [c.strip() for c in header if c.strip()]
        data: dict[str, list[float]] = {col: [] for col in cols}
        for row in reader:
            if not row or not any(field.strip() for field in row):
                continue
            for col, val in zip(cols, row):
                data[col].append(float(val.strip()))
    return {col: np.array(vals, dtype=np.float64) for col, vals in data.items()}


def list_snapshots(out_dir: PathLike) -> list[int]:
    """List available snapshot step numbers in out_dir, sorted in ascending order."""
    if not os.path.exists(out_dir):
        return []
    pattern = re.compile(r"^snap_(\d+)\.json$")
    steps: list[int] = []
    for entry in os.listdir(out_dir):
        match = pattern.match(entry)
        if match:
            step = int(match.group(1))
            bin_name = f"snap_{step:06d}.bin"
            if os.path.exists(os.path.join(out_dir, bin_name)):
                steps.append(step)
    steps.sort()
    return steps


def read_snapshot(out_dir: PathLike, step: int) -> dict[str, np.ndarray]:
    """Read snapshot step into a dictionary of arrays shaped [nz, ny, nx] with dtype as stored."""
    json_path = os.path.join(out_dir, f"snap_{step:06d}.json")
    bin_path = os.path.join(out_dir, f"snap_{step:06d}.bin")

    with open(json_path, "r", encoding="utf-8") as f:
        meta = json.load(f)
    with open(bin_path, "rb") as f:
        bin_data = f.read()

    result: dict[str, np.ndarray] = {}
    for var in meta.get("variables", []):
        name = var["name"]
        shape = tuple(var["shape"])
        dtype_str = var["dtype"]
        if dtype_str == "float32":
            dtype = np.float32
        elif dtype_str == "float64":
            dtype = np.float64
        else:
            raise ValueError(f"Unknown dtype '{dtype_str}' for variable '{name}'")

        offset = var["offset"]
        nbytes = var["nbytes"]
        raw = bin_data[offset : offset + nbytes]
        arr = np.frombuffer(raw, dtype=dtype).reshape(shape)
        # Make a copy so that it's a writable array detached from byte buffer
        result[name] = arr.copy()
    return result
