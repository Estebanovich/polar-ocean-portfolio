"""Figure 1 (Task 5): MITgcm (a) line plume against line-plume theory; C2, C2', C2''.

Usage: fig1_plume_validation.py [--runs-root DIR] [--window T0 T1] [--out-dir DIR]
                                [--section-case NAME]

Looks for ``<runs-root>/a_line2d_<strat>_q<q0>`` (strat = twolayer|sermilik,
q0 = 0.01|0.03|0.1); cases that are missing are skipped. Every run must be OK
(`plume_profile` calls `require_ok`). Writes ``fig1_plume_validation.{png,pdf}``
and ``fig1_metrics.json`` in ``--out-dir`` (default exp1_fjord_plumes/figures).

Criteria (global-constraints.md):
  C2   |z_detach(MITgcm) - z_detach(theory)| <= 0.10 H
  C2'  mean plume velocity in the lower half of the plume within +-25 %
       (MITgcm: width-averaged w; theory: top-hat u)
  C2'' mean face melt within a factor 1.5, theory with mode="mitgcm" and
       MeltParams(ice_depth=H) (MITgcm's own transfer law and linear ice T).
Ruling R13: all verdicts (C2, C2', C2'') use the theory in mode "mitgcm" with
MeltParams(ice_depth=H); the "jenkins"-mode numbers are context only.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from exp1_fjord_plumes.analysis.plume_diag import plume_profile  # noqa: E402
from exp1_fjord_plumes.gendata.make_case import _profile  # noqa: E402
from exp1_fjord_plumes.theory.melt import MeltParams  # noqa: E402
from exp1_fjord_plumes.theory.plume import line_plume  # noqa: E402
from pfcommon import style  # noqa: E402
from pfcommon.mitgcm_io import read_mds  # noqa: E402

DEFAULT_ROOT = Path("/Volumes/LaCie/portfolio_runs/exp1")
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
STRATS = ("twolayer", "sermilik")
QS = ("0.01", "0.03", "0.1")
XMAX = 1500.0  # m shown in the section

KNOWN_DIFFERENCES = [
    {"class": "c", "item": "MITgcm drives gammaT with |wVEL| of the 5 m face cell, the theory with the top-hat u"},
    {"class": "c", "item": "SW_TEMP theta -> in-situ temperature conversion in MITgcm (not in the theory)"},
    {"class": "c", "item": "no-slip wall in MITgcm vs Cd u^2 drag in the theory"},
    {"class": "b/c (small)",
     "item": "MITgcm adds melt as a virtual salt flux with no volume: the m_w term in dQ/dz is missing"},
]


def theory_solution(params: dict, mode: str):
    H = params["H"]
    z_grid = np.linspace(-H, 0.0, int(H) + 1)
    Ta, Sa = _profile(params["strat"], z_grid, params)
    p = MeltParams(ice_depth=H)
    return line_plume(z_grid, Ta, Sa, params["q0"], params["T0_src"], params["S0_src"], p, mode)


def _mean(a: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    return float(a.mean()) if a.size else float("nan")


def evaluate(prof: dict, params: dict) -> dict:
    """C2, C2', C2'' for one case from `plume_profile` output and the theory."""
    H = params["H"]
    z = prof["z"]
    out = {}
    for mode in ("jenkins", "mitgcm"):
        th = theory_solution(params, mode)
        zd_m, zd_t = prof["z_detach"], th.z_detach
        zlow = (-H + zd_t) / 2.0  # lower half of the (theory) plume, base to mid-height
        u_th = _mean(np.interp(z[(z <= zlow)], th.z, th.u))
        u_mo = _mean(prof["w_mean"][z <= zlow])
        m_mo = _mean(np.where(z < zd_m, prof["melt"], np.nan))
        m_th = _mean(th.melt)
        out[mode] = {
            "z_detach_model": zd_m, "z_detach_theory": zd_t,
            "z_detach_diff_over_H": (zd_m - zd_t) / H,
            "u_lower_model": u_mo, "u_lower_theory": u_th, "u_ratio": u_mo / u_th,
            "melt_model": m_mo, "melt_theory": m_th, "melt_ratio": m_mo / m_th,
            "theory_reached_top": bool(th.reached_top),
        }
    m = out["mitgcm"]  # R13: verdicts use the mitgcm-mode theory
    verdict = {
        "C2": bool(abs(m["z_detach_diff_over_H"]) <= 0.10),
        "C2p": bool(abs(m["u_ratio"] - 1.0) <= 0.25),
        "C2pp": bool(1 / 1.5 <= m["melt_ratio"] <= 1.5),
    }
    return {"modes": out, "pass": verdict, "melt_mean_whole_face_model": _mean(prof["melt"])}


def find_cases(root: Path) -> dict[tuple[str, str], Path]:
    return {(s, q): root / f"a_line2d_{s}_q{q}" for s in STRATS for q in QS
            if (root / f"a_line2d_{s}_q{q}" / "run_status.json").exists()}


def _final_section(run: Path):
    T, mt = read_mds(run / "T")
    W, _ = read_mds(run / "W")
    return T[0], W[0]


def make_figure(cases, profs, params, section_case, stem: Path):
    import matplotlib.pyplot as plt

    style.apply()
    fig = plt.figure(figsize=(7.5, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.5, 1, 1])
    ax_s = fig.add_subplot(gs[:, 0])
    axw = [fig.add_subplot(gs[r, 1]) for r in range(2)]
    axm = [fig.add_subplot(gs[r, 2]) for r in range(2)]

    key = section_case if section_case in cases else next(iter(cases))
    run = cases[key]
    prof = profs[key]
    T, W = _final_section(run)
    xf, zf, z = prof["xf"], prof["zf"], prof["z"]
    nz, nx = T.shape
    n = int(np.searchsorted(xf, XMAX))
    Wc = 0.5 * (W + np.vstack([W[1:], np.zeros((1, nx))]))
    pc = ax_s.pcolormesh(xf[: n + 1], zf, T[:, :n], cmap=style.CMAPS["temperature"], shading="flat")
    fig.colorbar(pc, ax=ax_s, label="T (°C)", shrink=0.8)
    cs = ax_s.contour(prof["x"][:n], z, Wc[:, :n], levels=[-0.1, -0.02, 0.02, 0.1, 0.3], colors="k", linewidths=0.6)
    ax_s.clabel(cs, fmt="%.2f", fontsize=6)
    ax_s.set_xlabel("distance from ice face (m)")
    ax_s.set_ylabel("z (m)")
    ax_s.set_title(f"{run.name}\nT and w (m/s)", fontsize=8)

    colors = dict(zip(QS, style.OKABE_ITO))
    for r, strat in enumerate(STRATS):
        for q in QS:
            if (strat, q) not in profs:
                continue
            pr, pa = profs[(strat, q)], params[(strat, q)]
            th = theory_solution(pa, "mitgcm")
            c = colors[q]
            axw[r].plot(pr["w_mean"], pr["z"], color=c, label=f"q0={q}")
            axw[r].plot(th.u, th.z, color=c, ls="--")
            axm[r].plot(pr["melt"] * 86400, pr["z"], color=c)
            axm[r].plot(th.melt * 86400, th.z, color=c, ls="--")
        axw[r].set_ylabel(f"{strat}\nz (m)")
        axw[r].set_ylim(-params[next(iter(params))]["H"], 0)
        axm[r].set_ylim(axw[r].get_ylim())
        axm[r].set_yticklabels([])
    axw[0].set_title("plume velocity", fontsize=8)
    axm[0].set_title("face melt", fontsize=8)
    axw[1].set_xlabel("w, u (m/s)")
    axm[1].set_xlabel("melt (m/day)")
    axw[0].plot([], [], "k-", label="MITgcm")
    axw[0].plot([], [], "k--", label="theory")
    axw[0].legend(loc="lower right", frameon=False)
    return style.save(fig, stem)


def main(runs_root: Path, window: tuple[float, float], out_dir: Path, section_case: str | None = None):
    cases = find_cases(runs_root)
    if not cases:
        raise SystemExit(f"no a_line2d_<strat>_q<q0> run found under {runs_root}")
    profs, params, metrics = {}, {}, {}
    for key, run in cases.items():
        profs[key] = plume_profile(run, window)  # RuntimeError if the run is not OK
        params[key] = json.loads((run / "case.json").read_text())["params"]
        metrics[run.name] = evaluate(profs[key], params[key])
    sec = next((k for k, r in cases.items() if r.name == section_case), None)
    make_figure(cases, profs, params, sec if sec else ("twolayer", "0.03"), out_dir / "fig1_plume_validation")
    allpass = {c: all(m["pass"][c] for m in metrics.values()) for c in ("C2", "C2p", "C2pp")}
    (out_dir / "fig1_metrics.json").write_text(json.dumps({
        "window_h": list(window), "cases": metrics, "all_pass": allpass,
        "known_differences": KNOWN_DIFFERENCES,
    }, indent=2, default=float))
    return metrics


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--window", type=float, nargs=2, default=(36.0, 48.0))
    ap.add_argument("--out-dir", type=Path, default=FIG_DIR)
    ap.add_argument("--section-case", default="a_line2d_twolayer_q0.03")
    a = ap.parse_args()
    main(a.runs_root, tuple(a.window), a.out_dir, a.section_case)
