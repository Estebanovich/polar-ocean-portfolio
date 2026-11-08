"""Realised pycnocline oscillation inside the (c) fjord in MITgcm (informative; same definition as
criterion 4 of the pyOCEAN cross-check, notes/pd1_c_fjord3d.md).

Pycnocline = depth of the T = 1.25 degC isotherm (eval_pyo_c.pycnocline_depth), y-mean over wet rows of a
column, from the 3-h mean THETA (mean3h record 1; a 3-h boxcar attenuates a 6-d sine by < 0.1 %).
Sine fit with period 6 d over days 6-30 (eval_pyo_c.sine_fit). Columns: i = 100 (x = 20 km, mid-fjord,
the pyOCEAN I_MID) and i = 25 (x = 5 km, near the ice face). Nominal forcing: 30 m (Jackson et al. 2018).

Usage: c_shelf_amplitude.py [--runs-root DIR]  -> figures/c_shelf_amplitude.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pfcommon.paths import EXP1_RUNS  # noqa: E402
from exp1_fjord_plumes.analysis import c_runs  # noqa: E402
from exp1_fjord_plumes.crosscheck.pyocean_c.eval_pyo_c import pycnocline_depth, sine_fit  # noqa: E402
from pfcommon.mitgcm_io import read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402

DAY = 86400.0
COLUMNS = {"x20km_mid": 100, "x5km_near_face": 25}
FIT_DAYS = (6.0, 30.0)
NOMINAL_M = 30.0
OUT = Path(__file__).resolve().parents[1] / "figures" / "c_shelf_amplitude.json"


def column_series(run: Path, i: int) -> tuple[np.ndarray, np.ndarray]:
    p = c_runs.params(run)["params"]
    dep = c_runs.depth(run)
    zc = -(np.arange(p["nz"]) + 0.5) * p["dz"]
    rows = np.flatnonzero(dep[:, i] > 0)
    t, d = [], []
    for it, dd in c_runs.mds_iters(c_runs.case_parts(run), "mean3h"):
        T = read_mds(dd / "mean3h", it)[0][1]                 # (nz, ny, nx)
        t.append(it * p["deltaT"] - c_runs.MEAN_S / 2)
        d.append(np.nanmean([pycnocline_depth(T[zc >= -dep[j, i], j, i], zc[zc >= -dep[j, i]]) for j in rows]))
    return np.array(t), np.array(d)


def main(runs_root: Path) -> dict:
    out = {"definition": __doc__.split("\n\n")[1], "nominal_m": NOMINAL_M, "cases": {}}
    for case in ("c_fjord_shelf", "c_fjord_all"):
        run = Path(runs_root) / case
        for part in c_runs.case_parts(run):
            require_ok(part)
        res = {}
        for name, i in COLUMNS.items():
            t, d = column_series(run, i)
            sel = (t >= FIT_DAYS[0] * DAY) & (t <= FIT_DAYS[1] * DAY) & np.isfinite(d)
            f = sine_fit(t[sel], d[sel])
            f["amplitude_over_nominal"] = f["amplitude"] / NOMINAL_M
            f["n_samples"] = int(sel.sum())
            res[name] = f
        out["cases"][case] = res
    OUT.write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-root", type=Path, default=EXP1_RUNS)
    r = main(ap.parse_args().runs_root)
    for c, v in r["cases"].items():
        for k, f in v.items():
            print(f"{c:14s} {k:15s} A = {f['amplitude']:5.1f} m ({f['amplitude_over_nominal']:.2f} of 30)  "
                  f"lag {f['lag_days']:.2f} d  mean {f['mean']:.1f} m  rms resid {f['rms_residual']:.1f} m  n {f['n_samples']}")
