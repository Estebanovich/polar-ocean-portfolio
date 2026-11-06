"""Figure 2 (Task 6): total face melt of two 3-D point plumes against their separation; C4, C4'.

Usage: fig2_two_plume_melt.py [--runs-root DIR] [--window T0 T1] [--out-dir DIR]

Looks for ``<runs-root>/b_3d_<strat>_<tag>`` (strat = twolayer|sermilik, tag = single|x0_<m>).
Definitions are pre-registered (notes/pd1_b_twoplume3d.md, last section) and are not changed here:
  window 2-3 h (informational: 1-2 h and the relative change); M = face_total_melt;
  normalised M(x0) / (2 M_single) against x0_actual / x0c (x0c from case.json);
  C4  |z_detach(model, single) - z_detach(half-cone theory)| <= 0.15 H, BOTH stratifications;
  C4' M(x0 = 0) / (2 M_single) in [0.4, 0.8], BOTH stratifications.
Choices the pre-registration leaves open (declared in fig2_metrics.json): the "far from the face"
ambient is the y-mean (wet cells) window-mean T/S of the column nearest x = AMBIENT_X_M; sensitivity to
this x is reported. A case whose pre-registered window contains NaN (run diverged) is NOT part of
any verdict; it is drawn as an open marker using the 1-2 h window (post-hoc, labelled).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from exp1_fjord_plumes.analysis.plume_diag import (  # noqa: E402
    ambient_far_from_face, face_total_melt, face_w_profile)
from exp1_fjord_plumes.gendata.make_case import _profile  # noqa: E402
from exp1_fjord_plumes.theory.melt import MeltParams  # noqa: E402
from exp1_fjord_plumes.theory.plume import point_plume  # noqa: E402
from pfcommon import style  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402

DEFAULT_ROOT = Path("/Volumes/LaCie/portfolio_runs/exp1")
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
DATA = Path(__file__).resolve().parent / "data" / "cenedese_gatto_2016a_fig2.csv"
STRATS = ("twolayer", "sermilik")
X0S = (0, 20, 40, 80, 160, 320)
WIN_MAIN, WIN_INFO = (2.0, 3.0), (1.0, 2.0)
AMBIENT_X_M = 855.0  # m: outermost column outside the 100 m OB sponge (declared choice)
AMBIENT_SENS = (500.0, 700.0, 855.0)
C4_TOL, C4P_RANGE = 0.15, (0.4, 0.8)
LAB_ALPHA = 0.1  # entrainment constant of the lab theory [verificar]: not given numerically in the paper
LAB_H2_CM = 20.5  # cm, bottom-layer thickness H20 (Cenedese & Gatto 2016a, section 2)


def theory_zdetach(params: dict, ambient=None) -> float:
    """Half-cone point-plume detachment elevation [m] (Slater et al. 2016 eqs. 1a-1d; mode "mitgcm"
    as R13, MeltParams(ice_depth=H)). ``ambient`` = (z centres surface-first, T, S) or None = t = 0
    profile (make_case._profile)."""
    H = params["H"]
    z_grid = np.linspace(-H, 0.0, int(H) + 1)
    if ambient is None:
        Ta, Sa = _profile(params["strat"], z_grid, params)
    else:
        zc, Tm, Sm = (np.asarray(a, float)[::-1] for a in ambient)
        Ta, Sa = np.interp(z_grid, zc, Tm), np.interp(z_grid, zc, Sm)
    sol = point_plume(z_grid, Ta, Sa, params["Q0"], params["T0_src"], params["S0_src"],
                      MeltParams(ice_depth=H), "mitgcm", alpha=params.get("alpha", 0.1),
                      geometry="half_cone")
    return float(sol.z_detach)


def _case(run: Path):
    return json.loads((run / "case.json").read_text())


def c4_for(run_single: Path, window) -> dict:
    """C4 numbers for one stratification from its `_single` run."""
    meta = _case(run_single)
    p = meta["params"]
    H = p["H"]
    assert abs(meta["Q0_total"] - p["Q0"]) < 1e-9, "single run must carry Q0 only"
    zm = face_w_profile(run_single, window)["z_detach"]
    amb = ambient_far_from_face(run_single, window, AMBIENT_X_M)
    zt = theory_zdetach(p, (amb["z"], amb["T"], amb["S"]))
    zt0 = theory_zdetach(p)
    sens = {}
    for x in AMBIENT_SENS:
        a = ambient_far_from_face(run_single, window, x)
        sens[f"x={a['x']:.0f}m"] = theory_zdetach(p, (a["z"], a["T"], a["S"]))
    out = {"z_detach_model": zm, "z_detach_theory": zt, "diff_over_H": (zm - zt) / H,
           "pass": bool(abs(zm - zt) / H <= C4_TOL), "ambient_x_m": amb["x"],
           "ambient_sensitivity_z_detach_theory": sens,
           "ambient_sensitivity_diff_over_H": {k: (zm - v) / H for k, v in sens.items()},
           "post_hoc": {"label": "post-hoc: theory with the initial (t = 0) profile",
                        "z_detach_theory_initial": zt0, "diff_over_H_initial": (zm - zt0) / H,
                        "pass_initial": bool(abs(zm - zt0) / H <= C4_TOL)}}
    return out


def collect(root: Path, window=WIN_MAIN):
    """Per-case M for both windows and everything needed for the figure."""
    cases = {}
    for s in STRATS:
        for tag in ("single",) + tuple(f"x0_{x}" for x in X0S):
            run = root / f"b_3d_{s}_{tag}"
            if not (run / "run_status.json").exists():
                continue
            require_ok(run)
            meta = _case(run)
            m_main = face_total_melt(run, window)
            m_info = face_total_melt(run, WIN_INFO)
            cases[(s, tag)] = {
                "run": run.name, "x0_actual": meta.get("x0_actual", 0.0), "x0c": meta["params"]["x0c"],
                "Q0_total": meta["Q0_total"], "M_main": m_main, "M_info": m_info,
                "valid_main": bool(np.isfinite(m_main)),
                "rel_change_main_vs_info": (m_main - m_info) / m_info if np.isfinite(m_main) else float("nan")}
    return cases


def lab_points():
    rows = [l.strip().split(",") for l in DATA.read_text().splitlines() if l.strip() and not l.startswith("#")][1:]
    x0 = np.array([float(r[0]) for r in rows])
    smr = np.array([float(r[1]) for r in rows])
    single_half = smr[-1]  # open circle: single source of Qsg/2; dashed line = 2x this
    x0c = LAB_H2_CM * LAB_ALPHA / 0.35  # cm; Cenedese & Gatto 2016a eq. in sec. 2/conclusions: x0 < H2 alpha / 0.35
    sel = slice(0, len(rows) - 1)
    return x0[sel] / x0c, smr[sel] / (2 * single_half), x0c


def make_figure(cases, stem: Path):
    import matplotlib.pyplot as plt
    style.apply()
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    cols = dict(zip(STRATS, style.OKABE_ITO))
    for s in STRATS:
        ms = cases.get((s, "single"))
        if ms is None:
            continue
        for tag in [f"x0_{x}" for x in X0S]:
            c = cases.get((s, tag))
            if c is None:
                continue
            xs = c["x0_actual"] / c["x0c"]
            if c["valid_main"]:
                ax.plot(xs, c["M_main"] / (2 * ms["M_main"]), "o", color=cols[s], ms=5.5)
            else:
                ax.plot(xs, c["M_info"] / (2 * ms["M_info"]), "o", mfc="none", color=cols[s], ms=6)
        pts = [(cases[(s, f"x0_{x}")]["x0_actual"] / cases[(s, f"x0_{x}")]["x0c"],
                cases[(s, f"x0_{x}")]["M_main"] / (2 * ms["M_main"])) for x in X0S
               if (s, f"x0_{x}") in cases and cases[(s, f"x0_{x}")]["valid_main"]]
        if pts:
            ax.plot(*zip(*pts), "-", color=cols[s], lw=1, label=f"MITgcm, {s}")
    lx, ly, x0c_lab = lab_points()
    ax.plot(lx, ly, "s", color="k", mfc="none", ms=6,
            label="Cenedese & Gatto 2016a (hand-digitised)")
    ax.axhspan(*C4P_RANGE, color="0.9", zorder=0, label="C4' range 0.4-0.8")
    ax.axhline(1, color="0.5", ls="--", lw=0.8)
    ax.set_xlabel(r"$x_0^{actual}/x_{0c}$")
    ax.set_ylabel(r"$M(x_0)\,/\,(2\,M_{single})$")
    ax.set_xlim(-0.1, 6.0)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.text(0.01, 0.005,
             "Window 2-3 h. Open MITgcm markers: run diverged (NaN) inside 2-3 h; 1-2 h window, post-hoc.\n"
             "x0c = 57.14 m (nominal for sermilik: no interface). Lab x0c = H2*alpha/0.35, alpha = 0.1 assumed.",
             fontsize=6, va="bottom")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    return style.save(fig, stem)


def main(runs_root: Path, window, out_dir: Path):
    cases = collect(runs_root, tuple(window))
    out_dir.mkdir(parents=True, exist_ok=True)
    res = {"window_h": list(window), "info_window_h": list(WIN_INFO),
           "unit": "M in m3 s-1 of ice", "cases": {}, "C4": {}, "C4p": {}}
    for (s, tag), c in cases.items():
        res["cases"][c["run"]] = c
    for s in STRATS:
        ms, m0 = cases.get((s, "single")), cases.get((s, "x0_0"))
        if ms is None or m0 is None or not (ms["valid_main"] and m0["valid_main"]):
            res["C4"][s] = res["C4p"][s] = {"error": "single or x0_0 missing/invalid"}
            continue
        res["C4"][s] = c4_for(runs_root / ms["run"], tuple(window))
        r = m0["M_main"] / (2 * ms["M_main"])
        res["C4p"][s] = {"ratio": r, "pass": bool(C4P_RANGE[0] <= r <= C4P_RANGE[1])}
        res["C4p"][s]["ratio_info_window"] = m0["M_info"] / (2 * ms["M_info"])
        res["C4p"][s]["normalised"] = {t: cases[(s, t)]["M_main"] / (2 * ms["M_main"]) for (ss, t) in cases
                                         if ss == s and t != "single"}
    ok4 = all(isinstance(v, dict) and v.get("pass") for v in res["C4"].values())
    ok4p = all(isinstance(v, dict) and v.get("pass") for v in res["C4p"].values())
    res["verdict"] = {"C4": "PASS" if ok4 else "FAIL", "C4p": "PASS" if ok4p else "FAIL",
                      "rule": "both stratifications must pass"}
    res["invalid_cases_nan_in_window"] = [c["run"] for c in cases.values() if not c["valid_main"]]
    res["ambient_choice"] = (f"y-mean window-mean T/S of the column nearest x = {AMBIENT_X_M} m; "
                             "the pre-registration says only 'far from the face'")
    res["lab"] = {"file": str(DATA.name), "x0c_lab_cm": LAB_H2_CM * LAB_ALPHA / 0.35,
                  "note": "alpha = 0.1 assumed [verificar]; denominator = 2 x open circle (single 1.13 cm3/s)"}
    make_figure(cases, out_dir / "fig2_two_plume_melt")
    (out_dir / "fig2_metrics.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--window", type=float, nargs=2, default=WIN_MAIN)
    ap.add_argument("--out-dir", type=Path, default=FIG_DIR)
    a = ap.parse_args()
    r = main(a.runs_root, tuple(a.window), a.out_dir)
    print(json.dumps({k: r[k] for k in ("C4", "C4p", "verdict", "invalid_cases_nan_in_window")}, indent=1, default=float))
