"""Post hoc V4 forcing screening (notes/v4_revision_posthoc.md, registered before running).

For a 10-year screening run (S0, SW, SF, SWF), from the 12 monthly means of the last stage:
  - shelf_theta: volume-mean theta over y < 80 km (degC);
  - J: onshore heat transport across the J line (TW), heat_transport.J_from_means with the
    bolus (GM_vbT) and Redi (DFyE_TH) fluxes; also the three parts;
  - passes: shelf_theta < -1.0 degC and J > 0.
Writes <run>/screen_v4.json.

Usage: python -m exp5_adjoint_slope.analysis.screen_v4 <run_dir> [<run_dir> ...]"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pfcommon.campaign import nml_float
from pfcommon.mitgcm_io import read_bin, read_mds

from . import heat_transport as ht

SHELF_Y_MAX = 80e3          # m
THETA_MAX = -1.0            # degC


def shelf_mean(theta, hfac_c, drf, rac, yc, y_max: float = SHELF_Y_MAX) -> float:
    """Volume-weighted mean of theta (nr, ny, nx) over wet cells with y < y_max."""
    vol = np.asarray(hfac_c, float) * np.asarray(drf, float)[:, None, None] * np.asarray(rac, float)
    sel = (np.asarray(yc, float) < y_max)[None, :, None] & (vol > 0)
    w = np.where(sel, vol, 0.0)
    return float((np.asarray(theta, float) * w).sum() / w.sum())


def shelf_bottom_mean(theta, hfac_c, yc, y_max: float = SHELF_Y_MAX) -> float:
    """Mean theta of the bottom wet cell over the wet columns with y < y_max."""
    wet = np.asarray(hfac_c, float) > 0
    n = wet.sum(axis=0)
    kb = np.maximum(n - 1, 0)
    tb = np.take_along_axis(np.asarray(theta, float), kb[None], axis=0)[0]
    sel = (np.asarray(yc, float) < y_max)[:, None] & (n > 0)
    return float(tb[sel].mean())


def passes(shelf_theta: float, J: float) -> bool:
    return bool(shelf_theta < THETA_MAX and J > 0)


def screen(run: Path) -> dict:
    run = Path(run)
    case = json.loads((run / "input" / "case.json").read_text())
    last = sorted((d for d in run.glob("stage*") if d.is_dir()),
                  key=lambda d: int(d.name[5:]))[-1]
    its = sorted(int(f.name.split(".")[1]) for f in last.glob("state3d.*.data"))[-12:]
    recs = np.array([np.asarray(read_mds(last / "state3d", it)[0]) for it in its])
    mean = recs.mean(axis=0)   # THETA, SALT, UVEL, VVEL, VVELTH, DFyE_TH, GM_vbT
    nr, ny, nx = mean.shape[1:]
    hfc = read_bin(run / "hFacC.data", (nr, ny, nx))
    hfs = read_bin(run / "hFacS.data", (nr, ny, nx))
    drf = read_bin(run / "DRF.data", (nr,))
    rac = read_bin(run / "RAC.data", (ny, nx))
    yc = read_bin(run / "YC.data", (ny, nx))[:, 0]
    j = case["j_qb"]
    dxg = read_bin(run / "DXG.data", (ny, nx))[j]
    data = (last / "data").read_text() if (last / "data").exists() else (run / "data").read_text()
    rho0, cp = nml_float(data, "rhoConst"), nml_float(data, "HeatCapacity_Cp")
    kw = dict(theta_f=case["theta_f"], rho0=rho0, cp=cp)
    J = ht.J_from_means(mean[4], mean[3], hfs, drf, dxg, j, dfy=mean[5], vbt=mean[6], **kw)
    Jbol = ht.J_from_means(mean[4], mean[3], hfs, drf, dxg, j, vbt=mean[6], **kw)
    Jres = ht.J_from_means(mean[4], mean[3], hfs, drf, dxg, j, **kw)
    th = shelf_mean(mean[0], hfc, drf, rac, yc)
    out = {"case": run.name, "stage": last.name, "iters": [its[0], its[-1]],
           "shelf_theta_degC": th,
           "shelf_bottom_theta_degC": shelf_bottom_mean(mean[0], hfc, yc),
           "theta_range_degC": [float(mean[0][hfc > 0].min()), float(mean[0][hfc > 0].max())],
           "shelf_salt": shelf_mean(mean[1], hfc, drf, rac, yc),
           "J_TW": J, "J_resolved_TW": Jres, "J_bolus_TW": Jbol - Jres, "J_redi_TW": J - Jbol,
           "passes": passes(th, J)}
    (run / "screen_v4.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("runs", type=Path, nargs="+")
    for r in ap.parse_args().runs:
        print(json.dumps(screen(r.resolve()), indent=1))
