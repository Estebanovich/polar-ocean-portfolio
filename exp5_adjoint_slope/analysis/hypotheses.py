"""H5a and H5c (plan Task 6), with the PD1b scales (y < 275 km, break y0 = 125 km, half-widths
37.5 / 18.75 / 75 km) and the PD2 inputs (EA ensemble-adjoint sensitivity, leads <= 60 d).
Thin interface over analysis.h5a and analysis.ensemble.

Usage: python -m exp5_adjoint_slope.analysis.hypotheses   (writes figures/hypotheses.json from
figures/validation.json and figures/sensitivity.json)"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import ensemble, h5a as _h5a

Y_MAX, Y0, WIDTHS = _h5a.Y_MAX, _h5a.Y0, _h5a.WIDTHS
FIG = Path(__file__).resolve().parents[1] / "figures"


def pattern_uniform(y_m) -> np.ndarray:
    """(i) shape: 1 for y < 275 km, 0 beyond (stronger easterlies when multiplied by -amplitude)."""
    return np.where(np.asarray(y_m, float) < Y_MAX, 1.0, 0.0)


def pattern_dipole(y_m, half_width_m: float, y0_m: float = Y0) -> np.ndarray:
    """(ii) shape: derivative of a Gaussian centred at the break, x exp(-x^2), x = (y - y0)/w."""
    x = (np.asarray(y_m, float) - y0_m) / half_width_m
    return x * np.exp(-x ** 2)


def normalize_rms(p, y_m, y_max_m: float = Y_MAX) -> np.ndarray:
    p = np.asarray(p, float)
    m = np.asarray(y_m, float) < y_max_m
    return p / np.sqrt((p[m] ** 2).mean())


def h5a(sens_y, y_m, widths_m=WIDTHS) -> dict:
    """Frozen verdict with one sensitivity profile (or members x ny): pass at the first width."""
    r = _h5a.evaluate(np.atleast_2d(sens_y), y_m)
    first = f"{widths_m[0] / 1e3:g}"
    return {**r, "pass": r["widths"][first]["pass"]}


def h5c(sens_month, lead_days) -> dict:
    lead = np.asarray(lead_days, float)
    o = np.argsort(lead)
    r = ensemble.h5c(lead[o], np.asarray(sens_month, float)[o])
    return {**r, "pass": r["verdict"] == "pass"}


def main() -> None:
    v = json.loads((FIG / "validation.json").read_text())["B"]
    out = {"note": "post hoc (PD2): EA21-EA30 ensemble adjoint, 90-day windows",
           "H5a": v["H5a"], "H5c": v["H5c"]}
    (FIG / "hypotheses.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
