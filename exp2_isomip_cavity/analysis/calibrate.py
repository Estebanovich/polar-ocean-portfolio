#!/opt/anaconda3/envs/pyocean/bin/python
r"""K1 calibration of the ISOMIP+ heat-exchange coefficient Gamma_T (Exp. 2, Task 4).

Usage::

    calibrate.py --res 2km|4km [--target 30] [--out figures/calibration_<res>.json] RUN_DIR...

The binding rule is pd0 K1 (``notes/pd0_preregistro.md``), after AD16 sect. 3.2.1, Eq. 37:

* quantity: area-weighted mean melt [m a^-1] of Ocean0 over the deep shelf (ice base
  < -300 m, AD16 Eq. 37) in the last 6 months of a 2-year run (monthly records 19-24),
  computed with ``metrics.mean_melt(rd, calendar_window(rd, 'last6'), 'deep300')``;
* samples: Gamma_T in {0.015, 0.03, 0.06}, Gamma_S = Gamma_T / 35 (AD16 Table 4);
* interpolation: ln<m> against ln Gamma_T, piecewise linear between neighbouring samples,
  to <m> = 30 m a^-1. **No extrapolation**: if 30 is not inside [min <m>, max <m>] the
  result is NOT_BRACKETED and the calibration goes to PD2 (exit code 3);
* verification: one run at Gamma_T,cal; K1 passes if 28 <= <m> <= 32 (``k1_verdict``). If
  it fails, pd0 allows exactly one more iteration with the 4 points (the same function).

Every run is read through ``metrics.load_run``, which calls ``pfcommon.runctl.require_ok``
first: a run that is not OK, did not reach nTimeSteps, or has NaN in its log (pd0 N3)
raises and stops the calibration (no partial calibration from the remaining runs).

The printed (last stdout line) and stored ``gammaT_cal`` is rounded to 5 significant digits,
the value passed to the verification run with ``--set gammaT=``; ``gammaT_cal_raw`` is the
unrounded interpolant (relative difference <= 5e-5, i.e. <= ~5e-5 in melt for slope ~1).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402

TARGET_M_PER_A = 30.0          # pd0 K1 / AD16 Eq. 37: deep Ocean0 melt, m a^-1
K1_TOL_M_PER_A = 2.0           # pd0 K1: 30 +/- 2 m a^-1 (closed interval 28-32)
GAMMA_S_RATIO = 35.0           # pd0 K1 / AD16 Table 4: Gamma_S = Gamma_T / 35
SIG_DIGITS = 5


class NotBracketed(ValueError):
    """The target is not inside the sampled melt range (or the samples are not monotonic):
    pd0 K1 forbids extrapolation -> report and go to PD2."""


# ---- one run ----------------------------------------------------------------------------

def calibration_point(run_dir) -> dict:
    """K1 quantity of one calibration run, plus what is needed to check it (N1, ids)."""
    rd = mt.load_run(run_dir)                     # require_ok (N3) before anything is read
    p = rd.params
    win = mt.calendar_window(rd, "last6")         # asserts 24 monthly records
    return {
        "run_dir": str(rd.run_dir),
        "case": p.get("case"),
        "gammaT": float(p["gammaT"]),
        "gammaS": float(p["gammaS"]),
        "dx_km": p.get("dx_km"),
        "deltaT": p.get("deltaT"),
        "init_profile": p.get("init_profile"),
        "restore_profile": p.get("restore_profile"),
        "s": p.get("s", rd.actual.get("s")),
        "deep_melt_last6_m_per_a": mt.mean_melt(rd, win, "deep300"),
        "shelf_melt_last6_m_per_a": mt.mean_melt(rd, win, "shelf"),
        "N1_drift": mt.equilibrium_drift(rd, "2yr"),
        "N3": "PASS",                             # load_run would have raised otherwise
    }


def _check_consistent(points: list) -> None:
    """All runs: Ocean0 (init and restoring profile ocean0, s = 1), Gamma_S = Gamma_T / 35,
    one resolution and one dt."""
    for pt in points:
        if pt["init_profile"] != "ocean0" or pt["restore_profile"] != "ocean0" \
                or float(pt["s"] if pt["s"] is not None else 1.0) != 1.0:
            raise ValueError(f"{pt['run_dir']}: not an Ocean0 run (init "
                             f"{pt['init_profile']}, restore {pt['restore_profile']}, "
                             f"s {pt['s']})")
        if not math.isclose(pt["gammaS"], pt["gammaT"] / GAMMA_S_RATIO, rel_tol=1e-12):
            raise ValueError(f"{pt['run_dir']}: gammaS {pt['gammaS']} != gammaT / 35")
    for key, what in (("dx_km", "resolution"), ("deltaT", "time step")):
        vals = {pt[key] for pt in points}
        if len(vals) > 1:
            raise ValueError(f"calibration runs mix {what}s: {sorted(vals)}")


# ---- interpolation ----------------------------------------------------------------------

def interpolate_loglog(gammas, melts, target: float = TARGET_M_PER_A):
    """Gamma_T at which the piecewise-linear interpolant of ln<m>(ln Gamma_T) through the
    samples equals ``target`` (pd0 K1). Returns (gamma_cal, (gamma_lo, gamma_hi)).

    Raises ValueError for repeated Gamma_T or non-positive values; NotBracketed if the
    target is outside the sampled range (no extrapolation) or the samples are not strictly
    monotonic in Gamma_T (the root would be ambiguous)."""
    if len(gammas) != len(melts) or len(gammas) < 2:
        raise ValueError("need at least 2 (gammaT, melt) pairs of equal length")
    if any(g <= 0 for g in gammas) or any(m <= 0 for m in melts):
        raise ValueError(f"gammaT and melt must be > 0 for ln: {gammas}, {melts}")
    pairs = sorted(zip(map(float, gammas), map(float, melts)))
    g = [a for a, _ in pairs]
    m = [b for _, b in pairs]
    if len(set(g)) != len(g):
        raise ValueError(f"repeated gammaT in {g}")
    d = [m[i + 1] - m[i] for i in range(len(m) - 1)]
    if not (all(x > 0 for x in d) or all(x < 0 for x in d)):
        raise NotBracketed(f"melt is not strictly monotonic in gammaT ({g} -> {m}); "
                           "K1 interpolation is ambiguous: report and go to PD2")
    lo, hi = min(m), max(m)
    if not lo <= target <= hi:
        raise NotBracketed(f"target {target} m/yr not bracketed by the sampled deep melt "
                           f"[{lo:.4g}, {hi:.4g}] (gammaT {g[0]}-{g[-1]}); no extrapolation "
                           "(pd0 K1): report and go to PD2")
    lt = math.log(target)
    for i in range(len(g) - 1):
        a, b = m[i], m[i + 1]
        if min(a, b) <= target <= max(a, b):
            if target == a:
                return g[i], (g[i], g[i + 1])
            if target == b:
                return g[i + 1], (g[i], g[i + 1])
            # ln Gamma = ln g_i + (ln t - ln m_i) (ln g_i+1 - ln g_i) / (ln m_i+1 - ln m_i)
            x = math.log(g[i]) + (lt - math.log(a)) * (math.log(g[i + 1]) - math.log(g[i])) \
                / (math.log(b) - math.log(a))
            return math.exp(x), (g[i], g[i + 1])
    raise AssertionError("unreachable: bracketed target without a segment")


def k1_verdict(deep_melt: float, target: float = TARGET_M_PER_A,
               tol: float = K1_TOL_M_PER_A) -> dict:
    """pd0 K1 for a verification run: pass if target - tol <= <m>_deep <= target + tol."""
    return {"deep_melt_last6_m_per_a": float(deep_melt), "range": [target - tol, target + tol],
            "pass": bool(target - tol <= deep_melt <= target + tol)}


def round_sig(x: float, n: int = SIG_DIGITS) -> float:
    return float(f"{x:.{n}g}")


# ---- driver -----------------------------------------------------------------------------

def calibrate(run_dirs, target: float = TARGET_M_PER_A) -> dict:
    """Read every run (any N3 failure raises), check consistency, interpolate. Returns a
    report dict with status 'OK' or 'NOT_BRACKETED' (gammaT_cal None)."""
    points = [calibration_point(r) for r in run_dirs]
    _check_consistent(points)
    points.sort(key=lambda p: p["gammaT"])
    g = [p["gammaT"] for p in points]
    m = [p["deep_melt_last6_m_per_a"] for p in points]
    rep = {"rule": "pd0 K1 (AD16 sect. 3.2.1, Eq. 37): deep (< -300 m) Ocean0 melt, last 6 "
                   "months (records 19-24), piecewise-linear ln<m> vs ln GammaT, no "
                   "extrapolation",
           "target_m_per_a": target, "k1_range": [target - K1_TOL_M_PER_A,
                                                   target + K1_TOL_M_PER_A],
           "dx_km": points[0]["dx_km"], "deltaT": points[0]["deltaT"], "points": points}
    try:
        gc, seg = interpolate_loglog(g, m, target)
    except NotBracketed as e:
        rep.update(status="NOT_BRACKETED", gammaT_cal=None, gammaT_cal_raw=None,
                   segment=None, message=str(e))
        return rep
    rep.update(status="OK", gammaT_cal_raw=gc, gammaT_cal=round_sig(gc),
               gammaS_cal=round_sig(gc) / GAMMA_S_RATIO, segment=list(seg),
               message="GammaT_cal from the bracketing segment; verify with one run (K1)")
    return rep


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="pd0 K1 GammaT calibration from Ocean0 runs")
    ap.add_argument("run_dirs", nargs="+", type=Path)
    ap.add_argument("--res", required=True, help="resolution label, e.g. 2km or 4km")
    ap.add_argument("--target", type=float, default=TARGET_M_PER_A)
    ap.add_argument("--out", type=Path, default=None,
                    help="JSON report (default figures/calibration_<res>.json)")
    a = ap.parse_args(argv)
    rep = {"resolution": a.res, **calibrate(a.run_dirs, a.target)}
    out = a.out or REPO / "exp2_isomip_cavity" / "figures" / f"calibration_{a.res}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=2, default=float) + "\n")
    for p in rep["points"]:
        print(f"# gammaT {p['gammaT']:<8g} deep melt {p['deep_melt_last6_m_per_a']:.4f} m/yr "
              f"(shelf {p['shelf_melt_last6_m_per_a']:.4f}, N1 {p['N1_drift']:.4f})",
              file=sys.stderr)
    if rep["status"] != "OK":
        print(f"# {rep['status']}: {rep['message']}", file=sys.stderr)
        return 3
    print(f"{rep['gammaT_cal']:.{SIG_DIGITS}g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
