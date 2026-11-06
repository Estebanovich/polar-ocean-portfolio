#!/usr/bin/env python
"""Evaluator of the pyOCEAN cross-check of the (c) open boundary.

Applies, mechanically, the four criteria pre-registered in ``exp1_fjord_plumes/notes/pd1_c_fjord3d.md``
(last section, 2026-11-06) to the runs written by ``run_pyo_c.py``:

1. reaches 30 d with finite fields;
2. volume balance: Q = 0: |eta_mean| < 0.01 m at all times; Q = 100: the net outgoing OB transport averaged over
   days 20-30 is Q +- 2 %, and the drift of eta_mean over days 20-30 is < 0.01 m;
3. at each daily output, max|u| on the 3 columns next to the OB <= 1.5 x max|u| in the interior (columns 10-189);
4. informative: amplitude and phase of the realised mid-fjord pycnocline (6-day-period sine fit, days 6-30)
   against the nominal 30 m.

Operational choices, fixed before any run was evaluated (also in README.md):
* "at all times" (2, Q = 0) = every model step (the run keeps the hourly max of |eta_mean| over steps);
* "drift" (2, Q = 100) = max over t in [20, 30] d of |eta_mean(t) - eta_mean(20 d)| (hourly records): the
  strictest reading;
* the day-20-30 mean transport = (cumulative outflow volume at 30 d - at 20 d)/(10 d), from per-step sums;
* (3) daily outputs are days 1..30 (day 0 is the state at rest, u = 0); faces nx-2 .. nx (the last is the OB
  face), interior u-faces 10..189, every layer and row;
* (4) the pycnocline is the depth of the isotherm T_MID = (T_PW + T_AW)/2 = 1.25 degC, the midpoint of the
  Sermilik profile (it sits at z_c exactly: profiles.continuous_sermilik), linearly interpolated, y-mean over
  the wet rows of column nx/2, hourly.  Fit: z = a - A sin(omega t - phi), omega = 2 pi/6 d; phi > 0 = lag
  behind the exterior z_ext = 200 - 30 sin(omega t).

Writes ``exp1_fjord_plumes/figures/crosscheck_pyo_c.json``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

T_MID = 0.5 * (-1.5 + 4.0)          # [degC] (T_PW + T_AW)/2, Sciascia et al. 2013 Table 2 (gendata/profiles.py)
DAY = 86400.0
PERIOD = 6.0 * DAY                  # [s] shelf period (cases.yaml shelf_period_days)
A_NOMINAL = 30.0                    # [m] shelf_amp
ETA_TOL = 0.01                      # [m] criteria 2
Q_TOL = 0.02                        # relative, criterion 2
RATIO_MAX = 1.5                     # criterion 3
FIT_WINDOW = (6.0, 30.0)            # [d] criterion 4
BAL_WINDOW = (20.0, 30.0)           # [d] criterion 2
PKG = Path(__file__).resolve().parents[2]
RUNS = Path("/Volumes/LaCie/portfolio_runs/exp1/crosscheck_pyo")


def pycnocline_depth(T: np.ndarray, z: np.ndarray, T_mid: float = T_MID) -> float:
    """Depth [m, positive] of the first T = T_mid crossing going down from the surface, linear in z between
    the two layer centres around it.  T and z (layer centres, m, negative down) in any vertical order.
    NaN if there is no crossing."""
    o = np.argsort(-np.asarray(z, float))                   # surface first
    zz, TT = np.asarray(z, float)[o], np.asarray(T, float)[o] - T_mid
    s = np.flatnonzero((TT[:-1] <= 0.0) & (TT[1:] > 0.0) | (TT[:-1] > 0.0) & (TT[1:] <= 0.0))
    if s.size == 0:
        return float("nan")
    k = s[0]
    w = TT[k] / (TT[k] - TT[k + 1])
    return float(-(zz[k] + w * (zz[k + 1] - zz[k])))


def sine_fit(t: np.ndarray, z: np.ndarray, period: float = PERIOD) -> dict:
    """Least squares z = a + b sin(w t) + c cos(w t); with z = a - A sin(w t - phi): A = hypot(b, c),
    phi = atan2(c, -b) [rad] (phi > 0: lag behind -sin(w t))."""
    w = 2 * np.pi / period
    G = np.column_stack([np.ones_like(t), np.sin(w * t), np.cos(w * t)])
    (a, b, c), *_ = np.linalg.lstsq(G, z, rcond=None)
    A, phi = float(np.hypot(b, c)), float(np.arctan2(c, -b))
    resid = z - G @ np.array([a, b, c])
    return dict(mean=float(a), amplitude=A, phase_rad=phi, phase_deg=float(np.degrees(phi)),
                lag_days=float(phi / w / DAY), rms_residual=float(np.sqrt(np.mean(resid ** 2))))


def criterion1(status: dict, days: float = 30.0) -> dict:
    ok = bool(status.get("completed")) and bool(status.get("finite")) and status.get("t_s", 0) >= days * DAY - 1
    return dict(passed=ok, completed=bool(status.get("completed")), finite=bool(status.get("finite")),
                t_end_days=float(status.get("t_s", 0)) / DAY, error=status.get("error"))


def criterion2(t: np.ndarray, eta_mean: np.ndarray, eta_absmax: np.ndarray, cum_out: np.ndarray, Q: float) -> dict:
    """t [s] hourly; eta_mean [m] at t; eta_absmax [m] max over the steps of each hour; cum_out [m^3]
    cumulative outward OB volume at t."""
    if Q == 0:
        mx = float(np.nanmax(eta_absmax)) if eta_absmax.size else float("nan")
        return dict(passed=bool(np.isfinite(mx) and mx < ETA_TOL), max_abs_eta_mean_m=mx, Q=0.0)
    t0, t1 = BAL_WINDOW[0] * DAY, BAL_WINDOW[1] * DAY
    if t.size == 0 or t[-1] < t1 - 1:
        return dict(passed=False, reason="series does not reach day 30", Q=Q)
    tt = np.concatenate([[0.0], t])
    cc = np.concatenate([[0.0], cum_out])
    q_mean = float((np.interp(t1, tt, cc) - np.interp(t0, tt, cc)) / (t1 - t0))
    e20 = float(np.interp(t0, t, eta_mean))
    sel = (t >= t0 - 1) & (t <= t1 + 1)
    drift = float(np.max(np.abs(eta_mean[sel] - e20)))
    rel = q_mean / Q - 1.0
    return dict(passed=bool(abs(rel) <= Q_TOL and drift < ETA_TOL), Q=Q, Q_out_mean_d20_30=q_mean,
                rel_error=rel, q_passed=bool(abs(rel) <= Q_TOL), eta_drift_d20_30_m=drift,
                drift_passed=bool(drift < ETA_TOL),
                eta_mean_d20_m=e20, eta_mean_d30_m=float(np.interp(t1, t, eta_mean)))


def criterion3(days: np.ndarray, umax_ob3: np.ndarray, umax_int: np.ndarray) -> dict:
    r = np.where(umax_int > 0, umax_ob3 / np.where(umax_int > 0, umax_int, 1.0), np.inf)
    k = int(np.argmax(r)) if r.size else 0
    return dict(passed=bool(r.size > 0 and np.all(r <= RATIO_MAX)), max_ratio=float(r.max()) if r.size else None,
                day_of_max=int(days[k]) if r.size else None, n_days=int(r.size),
                n_days_failed=int(np.sum(r > RATIO_MAX)), ratio_by_day=[float(x) for x in r],
                umax_ob3_by_day=[float(x) for x in umax_ob3], umax_int_by_day=[float(x) for x in umax_int])


def criterion4(t: np.ndarray, pyc: np.ndarray) -> dict:
    sel = (t >= FIT_WINDOW[0] * DAY - 1) & (t <= FIT_WINDOW[1] * DAY + 1) & np.isfinite(pyc)
    if sel.sum() < 10:
        return dict(informative=True, reason="not enough records in days 6-30")
    f = sine_fit(t[sel], pyc[sel])
    f.update(informative=True, nominal_amplitude=A_NOMINAL, amplitude_ratio=f["amplitude"] / A_NOMINAL,
             n_records=int(sel.sum()))
    return f


def daily_umax(run: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """max|u| next to the OB (u-faces nx-2..nx) and in the interior (10..189) of each day_XX.npz, XX >= 1."""
    days, ob3, inn = [], [], []
    for f in sorted(run.glob("day_*.npz")):
        d = int(f.stem.split("_")[1])
        if d == 0:
            continue
        u = np.abs(np.nan_to_num(np.load(f)["u"]))
        days.append(d)
        ob3.append(u[:, :, -3:].max())
        inn.append(u[:, :, 10:190].max())
    return np.array(days), np.array(ob3), np.array(inn)


def evaluate_run(run: Path, days: float = 30.0) -> dict:
    status = json.loads((run / "status.json").read_text())
    Q = float(status["config"]["Q_runoff"])
    s = np.load(run / "series.npz")
    cols = [str(c) for c in s["columns"]]
    data = s["data"] if s["data"].size else np.zeros((0, len(cols)))
    col = {c: data[:, i] for i, c in enumerate(cols)}
    dd, ob3, inn = daily_umax(run)
    return dict(case=status["config"]["case"], Q=Q,
                c1=criterion1(status, days),
                c2=criterion2(col["t"], col["eta_mean"], col["eta_mean_absmax_hour"], col["cum_ob_out"], Q),
                c3=criterion3(dd, ob3, inn),
                c4=criterion4(col["t"], col["pyc_mid"]),
                cost=dict(wall_s=status.get("wall_s"), s_per_step=status.get("s_per_step")))


def print_report(res: dict) -> None:
    for r in res["cases"]:
        print(f"== {r['case']} (Q = {r['Q']:g} m3/s)")
        c1, c2, c3, c4 = r["c1"], r["c2"], r["c3"], r["c4"]
        print(f"  C1 {'PASS' if c1['passed'] else 'FAIL'}: t_end {c1['t_end_days']:.3f} d, finite {c1['finite']}, "
              f"completed {c1['completed']}" + (f", error {c1['error']}" if c1.get("error") else ""))
        if r["Q"] == 0:
            print(f"  C2 {'PASS' if c2['passed'] else 'FAIL'}: max|eta_mean| = {c2['max_abs_eta_mean_m']:.3e} m")
        elif "Q_out_mean_d20_30" in c2:
            print(f"  C2 {'PASS' if c2['passed'] else 'FAIL'}: Q_out(d20-30) = {c2['Q_out_mean_d20_30']:.3f} m3/s "
                  f"(rel {c2['rel_error']:+.4f}), eta drift d20-30 = {c2['eta_drift_d20_30_m']:.3e} m")
        else:
            print(f"  C2 FAIL: {c2.get('reason')}")
        print(f"  C3 {'PASS' if c3['passed'] else 'FAIL'}: max ratio {c3['max_ratio']} (day {c3['day_of_max']}), "
              f"{c3['n_days_failed']}/{c3['n_days']} days > {RATIO_MAX}")
        if "amplitude" in c4:
            print(f"  C4 (informative): A = {c4['amplitude']:.2f} m ({c4['amplitude_ratio']:.3f} of 30 m), "
                  f"phase lag {c4['phase_deg']:.1f} deg ({c4['lag_days']:.2f} d), mean {c4['mean']:.1f} m")
        else:
            print(f"  C4 (informative): {c4.get('reason')}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=Path, default=RUNS)
    ap.add_argument("--cases", nargs="+", default=["pyo_shelf", "pyo_shelf_runoff"])
    ap.add_argument("--days", type=float, default=30.0)
    ap.add_argument("--json", type=Path, default=PKG / "figures" / "crosscheck_pyo_c.json")
    a = ap.parse_args()
    res = dict(criteria_source="exp1_fjord_plumes/notes/pd1_c_fjord3d.md (last section, 2026-11-06)",
               cases=[evaluate_run(a.runs / c, a.days) for c in a.cases])
    print_report(res)
    a.json.parent.mkdir(parents=True, exist_ok=True)
    a.json.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {a.json}")


if __name__ == "__main__":
    main()
