"""H6c and V6d for the 10-km model (spec Sec. 6-7).

H6c: dS = mean_m S_m(2060) - mean_m S_m(1900) over the members; detectable if |dS| > 2 sigma,
sigma = sqrt(s_1900^2/N + s_2060^2/N) (sample std, ddof 1), as sigma_delta of Exp. 5
(exp5_adjoint_slope/analysis/crosscheck.py); same sign as the coarse dS; ratio in [0.5, 2].
With fewer than 3 members the result is only "indicative".
V6d: eke_equilibrated of Exp. 5 (last 3 months within 10 % of the 3 before) on monthly EKE,
with at least 3 years of spin-up."""
from __future__ import annotations

import numpy as np

from exp5_adjoint_slope.analysis.crosscheck import eke_equilibrated, sigma_delta


def h6c(ds_fine_members: dict[int, list[float]], ds_coarse: float) -> dict:
    a, b = np.asarray(ds_fine_members[1900], float), np.asarray(ds_fine_members[2060], float)
    ds = float(b.mean() - a.mean())
    sig = sigma_delta(a, b)
    det = bool(abs(ds) > 2 * sig)
    sign = bool(np.sign(ds) == np.sign(ds_coarse) and ds != 0)
    ratio = ds / ds_coarse if ds_coarse else float("nan")
    within = bool(sign and 0.5 <= ratio <= 2.0)
    return {"dS": ds, "sigma": float(sig), "detectable": det, "sign_ok": sign, "ratio": float(ratio),
            "pass": bool(det and within), "indicative": bool(min(a.size, b.size) < 3), "n": int(min(a.size, b.size))}


def v6d(eke_monthly, min_years: float = 3.0) -> dict:
    s = list(eke_monthly)
    r = eke_equilibrated(s) if len(s) >= 6 else {"pass": False}
    enough = len(s) >= 12 * min_years
    out = dict(r, months=len(s), min_years_ok=bool(enough))
    out["pass"] = bool(r["pass"] and enough)
    return out


def eke_series(run) -> list[float]:
    """Monthly volume-mean EKE (m2 s-2) of a fine run from eke_mon (UVEL, VVEL, UVELSQ, VVELSQ),
    eke of Exp. 5 on the wet cells (hFacC of the run or of its pickup_from ancestor)."""
    from pathlib import Path

    from exp5_adjoint_slope.analysis.crosscheck import eke
    from . import io
    from pfcommon.mitgcm_io import read_bin
    run = Path(run)
    c = io.case(run)
    wet = read_bin(io.grid_dir(run) / "hFacC.data", (c["nr"], c["ny"], c["nx"])) > 0
    out = []
    for it in io.iters(run, "eke_mon"):
        f = io.fields(run, "eke_mon", it)
        out.append(eke(f["UVEL"], f["VVEL"], f["UVELSQ"], f["VVELSQ"], wet))
    return out


def f0_next_years(months_done: int, v6d_pass: bool, min_years: int = 3, max_years: int = 8) -> int:
    """Years to add to an F0 spin-up: up to the 3-year minimum first, then 1 at a time until V6d
    passes or the 8-year cap is reached (0 = stop)."""
    years = int(months_done) // 12
    if years < min_years:
        return min_years - years
    return 0 if (v6d_pass or years >= max_years) else 1
