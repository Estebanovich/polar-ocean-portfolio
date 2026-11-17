"""J anchor (plan Task 3): the Fortran cost (fc of a GJ run) against analysis.heat_transport.J
computed from the every-step snapshots of the same run. Writes <run>/anchor.json.

Usage: python -m exp5_adjoint_slope.analysis.anchor <run_dir>"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pfcommon.campaign import nml_float, nml_int
from pfcommon.mitgcm_io import read_bin, read_mds

from . import heat_transport as ht


def window_iters(n_total: int, last_steps: int) -> list[int]:
    """Iterations whose end-of-step state enters J: n_total - last_steps + 1 .. n_total."""
    return list(range(n_total - last_steps + 1, n_total + 1))


def rel_diff_extrapolated(series, fc: float, n_window: int) -> float:
    """|mean - fc| / |fc| with the missing final records (n_window - len(series)) extrapolated
    linearly from the last two. Post hoc diagnostic (2026-11-16): it explains rel_diff, which
    stays the reported figure."""
    s = list(np.asarray(series, dtype=float))
    while len(s) < n_window:
        s.append(2 * s[-1] - s[-2])
    m = float(np.mean(s))
    return abs(m - fc) / abs(fc)


def anchor(run: Path) -> dict:
    run = Path(run)
    case = json.loads((run / "input" / "case.json").read_text())
    data = (run / "data").read_text()
    n, dt = nml_int(data, "nTimeSteps"), nml_float(data, "deltaT")
    L = nml_float((run / "data.cost").read_text(), "lastinterval")
    # Diagnostics snapshots exist for iterations 1 .. n-1 (none is written after the last
    # step), so the anchor uses the window records that exist and reports how many.
    its = [i for i in window_iters(n, int(round(L / dt))) if (run / f"snap3d.{i:010d}.data").exists()]
    recs = []
    for it in its:
        a, _ = read_mds(run / "snap3d", it)
        recs.append(np.asarray(a))   # VVEL, THETA, GM_PsiY, GM_Kvy, GM_Kvz, GM_vbT, DFyE_TH
    recs = np.array(recs)
    v, th, psi, kvy, kvz, vbt, dfy = (recs[:, n] for n in range(7))
    nr, ny, nx = v.shape[1:]
    j = case["j_qb"]
    hfs = read_bin(run / "hFacS.data", (nr, ny, nx))          # grid files carry no iteration
    hfc = read_bin(run / "hFacC.data", (nr, ny, nx))
    drf = read_bin(run / "DRF.data", (nr,))
    drc = read_bin(run / "DRC.data", (nr + 1,))[:nr]          # drC(k), k = 1..Nr
    dxg = read_bin(run / "DXG.data", (ny, nx))[j]
    dyc = read_bin(run / "DYC.data", (ny, nx))[j]
    rho0 = nml_float(data, "rhoConst")
    cp = nml_float(data, "HeatCapacity_Cp")
    kw = dict(theta_f=case["theta_f"], rho0=rho0, cp=cp)
    # Advective form: GM_ExtraDiag is true (GM_isopycK /= 0, gmredi_readparms.F:242-247), so the
    # Redi flux has the Kvz term too.
    gm = ht.gm_yflux(th, kvy, hfs, drf, dxg, dyc, j, kvz=kvz, drc=drc, mask_c=(hfc > 0))
    jpy = ht.J(v, th, hfs, drf, dxg, j, gm, psi, **kw)
    series = ht.J_series(v, th, hfs, drf, dxg, j, gm, psi, **kw)
    jdiag = ht.J(v, th, hfs, drf, dxg, j, dfy[:, :, j, :] + vbt[:, :, j, :], None, **kw)
    jmean = ht.J(v, th, hfs, drf, dxg, j, None, None, **kw)
    jbol = ht.J(v, th, hfs, drf, dxg, j, None, psi, **kw) - jmean
    fc = json.loads((run / "fc.json").read_text())["fc"]
    out = {"fc_fortran_TW": fc, "J_python_TW": jpy, "J_python_from_model_fluxes_TW": jdiag,
           "J_resolved_TW": jmean, "J_bolus_TW": jbol, "J_redi_TW": jpy - jmean - jbol,
           "n_records": len(its), "n_window_steps": int(round(L / dt)),
           "rel_diff": abs(jpy - fc) / abs(fc) if fc else None,
           "rel_diff_last_extrapolated": rel_diff_extrapolated(series, fc, int(round(L / dt)))}
    (run / "anchor.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run", type=Path)
    print(json.dumps(anchor(ap.parse_args().run.resolve()), indent=1))
