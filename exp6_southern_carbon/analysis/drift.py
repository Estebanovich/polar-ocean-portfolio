"""V6a, the spin-up drift rule (spec Sec. 6; notes/pd0_preregistro.md)."""
from __future__ import annotations

import numpy as np

FLUX_TOL = 0.01                 # mol C m-2 yr-1 per century
DEEP_TOL = 1e-6 * 1035.0        # 1 umol kg-1 per century, in mol m-3 (rho0 = 1035)
MIN_YEARS = 1000.0


def _century_change(years, series) -> float:
    y, s = np.asarray(years, float), np.asarray(series, float)
    last = y >= y[-1] - 100.0
    return float(abs(s[last][-1] - s[last][0]))


def v6a(years, flux_decadal, deep_dic_decadal, run_years: float | None = None) -> dict:
    """years: tracer years of decadal means; flux in mol C m-2 yr-1; deep DIC (below 1000 m) in mol m-3.
    run_years: length of the run in tracer years (the ">= 1000 tracer years" of the rule); without
    it, the span of the decadal means is used, which is one decade shorter than the run."""
    span = float(run_years) if run_years is not None else float(np.asarray(years, float)[-1] - np.asarray(years, float)[0])
    df, dd = _century_change(years, flux_decadal), _century_change(years, deep_dic_decadal)
    ok_f, ok_d, ok_n = df < FLUX_TOL, dd < DEEP_TOL, span >= MIN_YEARS
    reason = [] if ok_n else [f"span {span:.0f} < {MIN_YEARS:.0f} years"]
    return {"d_flux_century": df, "d_deep_dic_century": dd, "span_years": span,
            "pass": bool(ok_f and ok_d and ok_n), "reason": "; ".join(reason)}
