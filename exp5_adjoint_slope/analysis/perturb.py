"""Adjoint-guided perturbations for V2/V3 (notes/pd0_preregistro.md): the pattern is the GA
sensitivity normalised to its maximum |value| and scaled to the PD0 amplitude (0.01 N m-2 for
tau_x, constant in time; 0.1 degC for T0; 1e-5 m2 s-1 for kappa_v). The xx_<name> files are read
by pkg/ctrl when doInitXX = .FALSE. (scripts/run_adjoint.py install_xx)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import write_bin

AMPLITUDE = {"fu": 0.01, "theta": 0.1, "diffkr": 1e-5}      # PD0 (propuesta), units of u


def pattern(sens, amp: float) -> np.ndarray:
    """amp * sens / max|sens|."""
    s = np.asarray(sens, dtype=float)
    return amp * s / np.abs(s).max()


def predicted_dj(sens, du) -> float:
    """First-order dJ = sum(dJ/du * du) (Losch & Heimbach 2007)."""
    return float((np.asarray(sens, dtype=float) * np.asarray(du, dtype=float)).sum())


def write_xx(out_dir: Path, name: str, arr) -> Path:
    """xx_<name>.0000000000.data/.meta, one record, big-endian float64, (ny, nx) or (nr, ny, nx)
    in the format pkg/ctrl writes (MDS meta as in GSF/xx_fu.0000000000.meta)."""
    a = np.asarray(arr, dtype=float)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    data, meta = out / f"xx_{name}.0000000000.data", out / f"xx_{name}.0000000000.meta"
    write_bin(data, a)
    dims = list(reversed(a.shape))                       # x first
    lines = ",\n".join(f"   {n:3d},    1,   {n:3d}" for n in dims)
    meta.write_text(
        f" nDims = [   {a.ndim} ];\n dimList = [\n{lines}\n ];\n dataprec = [ 'float64' ];\n"
        f" nrecords = [          1 ];\n timeStepNumber = [          0 ];\n")
    return data
