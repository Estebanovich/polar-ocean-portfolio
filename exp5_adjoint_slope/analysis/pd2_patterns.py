"""PD2 perturbation patterns (notes/pd0_preregistro.md, "Enmienda PD2"), written as pkg/ctrl
control files for the forward ensembles:

  xx_v2d/<u>_<p|m>  : +/- the GM30 pattern of u (tau_x, T0, kappa_v), PD0 amplitudes; pred.json
  xx_ens/<n>_<p|m>  : v3e = EA-mean tau_x(y) pattern (max 0.01 N m-2, also the FP pattern),
                      h5ai / h5aii = H5a uniform and 37.5-km dipole (RMS 0.01 over y < 275 km);
                      patterns.json (profiles, per-member predictions), fp_pattern.json (H5b).

All perturbations are time-constant and x-uniform; the runs use one gentim2d record
(ctrl_period_days 0) with periodic forcing reloaded every step (scripts/run_adjoint.py).

Usage: python -m exp5_adjoint_slope.analysis.pd2_patterns {v2d|ens}"""
from __future__ import annotations

import json
import sys

import numpy as np

from pfcommon.paths import EXP5_RUNS

from . import adjsens, h5a, perturb

YEARS = range(21, 31)
NR, NY, NX = 35, 45, 40
YC = (np.arange(NY) + 0.5) * 10e3


def _write_set(out, name: str, arr) -> None:
    z3, z2 = np.zeros((NR, NY, NX)), np.zeros((NY, NX))
    for m in ("fu", "fv", "theta", "salt", "diffkr"):
        perturb.write_xx(out, m, arr if m == name else (z2 if m in ("fu", "fv") else z3))


def v2d() -> dict:
    d = EXP5_RUNS / "GM30"
    sens = {"fu": adjsens.load_adxx(d, "fu").sum(0), "theta": adjsens.load_adxx(d, "theta")[0],
            "diffkr": adjsens.load_adxx(d, "diffkr")[0]}
    pred = {}
    for n, s in sens.items():
        du = perturb.pattern(s, perturb.AMPLITUDE[n])
        for sg, lab in ((1, "p"), (-1, "m")):
            _write_set(EXP5_RUNS / "xx_v2d" / f"{n}_{lab}", n, sg * du)
        pred[n] = perturb.predicted_dj(s, du)
    (EXP5_RUNS / "xx_v2d" / "pred.json").write_text(json.dumps(pred) + "\n")
    return pred


def ens() -> dict:
    S = np.array([adjsens.tau_x_sensitivity(adjsens.load_adxx(EXP5_RUNS / f"EA{y}", "fu"), None)
                  for y in YEARS]).sum(1)                       # (members, ny), summed over records
    pats = {"v3e": perturb.pattern(S.mean(0), 0.01), "h5ai": h5a.uniform(YC),
            "h5aii": h5a.dipole(YC, 37.5e3)}
    out = {"y_m": YC.tolist()}
    for n, pf in pats.items():
        for sg, lab in ((1, "p"), (-1, "m")):
            _write_set(EXP5_RUNS / "xx_ens" / f"{n}_{lab}", "fu", sg * np.repeat(pf[:, None], NX, 1))
        out[n] = {"profile": pf.tolist(), "pred_members": [float(np.dot(s, pf)) for s in S]}
    (EXP5_RUNS / "xx_ens" / "patterns.json").write_text(json.dumps(out) + "\n")
    (EXP5_RUNS / "xx_ens" / "fp_pattern.json").write_text(
        json.dumps({"y_m": YC.tolist(), "dtau": pats["v3e"].tolist()}) + "\n")
    return {n: float(np.mean(out[n]["pred_members"])) for n in pats}


if __name__ == "__main__":
    print(json.dumps({"v2d": v2d, "ens": ens}[sys.argv[1]](), indent=1))
