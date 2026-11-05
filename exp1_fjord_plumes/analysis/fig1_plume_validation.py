"""Figure 1 (Task 5): MITgcm (a) line plume against line-plume theory; C2, C2', C2''.

Usage: fig1_plume_validation.py [--runs-root DIR] [--window T0 T1] [--out-dir DIR]

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
XMAX = 500.0  # m shown in the section

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
    # Post-hoc metrics (R20): NOT part of the pre-registered verdicts.
    th = theory_solution(params, "mitgcm")
    zlow = (-H + th.z_detach) / 2.0
    low = z <= zlow
    u_th = _mean(np.interp(z[low], th.z, th.u))
    core = _mean(prof["w_face"][low])
    wmax = _mean(prof["w_max"][low])
    q_low = prof["Q"][low]
    mq = _mean(prof["M"][low] / np.where(q_low > 0, q_low, np.nan))
    post_hoc = {
        "label": "post-hoc, defined after seeing C2' fail; not a verdict",
        "u_theory_lower_half": u_th,
        "w_face_core_lower_half": core, "w_face_ratio": core / u_th,
        "w_max_lower_half": wmax, "w_max_ratio": wmax / u_th,
        "M_over_Q_lower_half": mq, "M_over_Q_ratio": mq / u_th,
        "width_10pct_lower_half": _mean(prof["width"][low]),
        "theory_b_lower_half": _mean(np.interp(z[low], th.z, th.b)),
        "Q_model_lower_half": _mean(q_low), "Q_theory_lower_half": _mean(np.interp(z[low], th.z, th.Q)),
    }
    return {"modes": out, "pass": verdict, "post_hoc": post_hoc,
            "melt_mean_whole_face_model": _mean(prof["melt"])}


def find_cases(root: Path) -> dict[tuple[str, str], Path]:
    return {(s, q): root / f"a_line2d_{s}_q{q}" for s in STRATS for q in QS
            if (root / f"a_line2d_{s}_q{q}" / "run_status.json").exists()}


def _final_section(run: Path):
    T, mt = read_mds(run / "T")
    W, _ = read_mds(run / "W")
    return T[0], W[0]


def make_figure(cases, profs, params, section_q, stem: Path):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    style.apply()
    fig = plt.figure(figsize=(7.6, 6.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1, 1])
    axs = [fig.add_subplot(gs[r, 0]) for r in range(2)]
    axw = [fig.add_subplot(gs[r, 1]) for r in range(2)]
    axm = [fig.add_subplot(gs[r, 2]) for r in range(2)]
    H = params[next(iter(params))]["H"]
    tags = "abcdef"

    # left: x-z sections of T with w contours (one case per stratification)
    for r, strat in enumerate(STRATS):
        key = (strat, section_q if (strat, section_q) in cases else next(q for s_, q in cases if s_ == strat))
        prof = profs[key]
        T, W = _final_section(cases[key])
        xf, zf, z = prof["xf"], prof["zf"], prof["z"]
        nz, nx = T.shape
        n = int(np.searchsorted(xf, XMAX))
        Wc = 0.5 * (W + np.vstack([W[1:], np.zeros((1, nx))]))
        ax = axs[r]
        pc = ax.pcolormesh(xf[: n + 1], zf, T[:, :n], cmap=style.CMAPS["temperature"],
                           shading="flat", vmin=-1.6, vmax=4.1)
        cs = ax.contour(prof["x"][:n], z, Wc[:, :n], levels=[0.1, 0.3, 0.5], colors="k", linewidths=0.7)
        ax.clabel(cs, fmt="%.2f", fontsize=6, inline_spacing=1)
        ax.set_ylabel("z (m)")
        ax.set_title(f"({tags[r]}) {strat}, q0 = {key[1]} m$^2$ s$^{{-1}}$\nT (colour), w = 0.1, 0.3, 0.5 m s$^{{-1}}$ (contours)",
                     fontsize=8, loc="left")
    axs[1].set_xlabel("distance from ice face (m)")
    cb = fig.colorbar(pc, ax=axs, location="bottom", shrink=0.8, aspect=30, pad=0.02)
    cb.set_label("potential temperature (°C)")

    colors = dict(zip(QS, style.OKABE_ITO))
    for r, strat in enumerate(STRATS):
        for q in QS:
            if (strat, q) not in profs:
                continue
            pr, pa = profs[(strat, q)], params[(strat, q)]
            th = theory_solution(pa, "mitgcm")
            c = colors[q]
            axw[r].plot(pr["w_mean"], pr["z"], color=c, lw=1.6)
            axw[r].plot(th.u, th.z, color=c, ls="--", lw=1.3)
            axw[r].plot(pr["w_face"], pr["z"], color=c, ls=":", lw=1.0)
            axm[r].plot(pr["melt"] * 86400, pr["z"], color=c, lw=1.6)
            axm[r].plot(th.melt * 86400, th.z, color=c, ls="--", lw=1.3)
        axw[r].set_ylim(-H, 0)
        axm[r].set_ylim(-H, 0)
        axw[r].set_ylabel("z (m)")
        axm[r].set_yticklabels([])
        axw[r].set_title(f"({tags[2 + r]}) {strat}: plume velocity", fontsize=8, loc="left")
        axm[r].set_title(f"({tags[4 + r]}) {strat}: face melt", fontsize=8, loc="left")
        axw[r].set_xlim(left=0)
        axm[r].set_xlim(left=0)
    axw[1].set_xlabel("plume velocity (m s$^{-1}$)")
    axm[1].set_xlabel("ice melt rate (m day$^{-1}$)")
    h = [Line2D([], [], color=colors[q], lw=1.6, label=f"q0 = {q} m$^2$ s$^{{-1}}$") for q in QS]
    h += [Line2D([], [], color="k", lw=1.6, label="MITgcm (width-mean w)"),
          Line2D([], [], color="k", ls=":", lw=1.0, label="MITgcm, face cell (post-hoc)"),
          Line2D([], [], color="k", ls="--", lw=1.3, label="theory (mitgcm mode)")]
    fig.legend(handles=h, loc="outside lower right", ncol=2, fontsize=7, frameon=False)
    return style.save(fig, stem)


def main(runs_root: Path, window: tuple[float, float], out_dir: Path):
    cases = find_cases(runs_root)
    if not cases:
        raise SystemExit(f"no a_line2d_<strat>_q<q0> run found under {runs_root}")
    profs, params, metrics = {}, {}, {}
    for key, run in cases.items():
        profs[key] = plume_profile(run, window)  # RuntimeError if the run is not OK
        params[key] = json.loads((run / "case.json").read_text())["params"]
        metrics[run.name] = evaluate(profs[key], params[key])
    make_figure(cases, profs, params, "0.03", out_dir / "fig1_plume_validation")
    allpass = {c: all(m["pass"][c] for m in metrics.values()) for c in ("C2", "C2p", "C2pp")}
    (out_dir / "fig1_metrics.json").write_text(json.dumps({
        "window_h": list(window), "cases": metrics, "all_pass": allpass,
        "known_differences": KNOWN_DIFFERENCES,
        "pending_convergence_study": "exp1_fjord_plumes/notes/pd2_convergence.md",
        "note": "C2' is reported as FAILED under its pre-registered definition (width-mean "
                "velocity, 10 % of w_max, lower half); post_hoc entries are not verdicts (R20).",
    }, indent=2, default=float))
    return metrics


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--window", type=float, nargs=2, default=(36.0, 48.0))
    ap.add_argument("--out-dir", type=Path, default=FIG_DIR)
    a = ap.parse_args()
    main(a.runs_root, tuple(a.window), a.out_dir)
