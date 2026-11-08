"""Mechanical evaluator of the pre-registered PD2 convergence/sensitivity study (R20, R21).

Implements exactly the criteria of exp1_fjord_plumes/notes/pd2_convergence.md
(window, stationarity, sigma, thresholds, validity checks, Richardson with r_eff,
outcome classes V/G/S/P1/P2/W/D/F/M evaluated in that order, separate sensitivity
results). Nothing here is tuned to data: constants are the note's.

Metrics reuse the C2' code unchanged: `plume_profile` (R22 interpolated width) +
`fig1_plume_validation.evaluate` (theory in mode "mitgcm", t = 0 ambient):
  R_Q = u_ratio, R_w = post_hoc w_max_ratio, b = post_hoc width_10pct_lower_half.

Usage: pd2_convergence_eval.py [--runs-root DIR] [--out-dir DIR]
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pfcommon.paths import EXP1_RUNS  # noqa: E402
from exp1_fjord_plumes.analysis.fig1_plume_validation import evaluate  # noqa: E402
from exp1_fjord_plumes.analysis.plume_diag import plume_profile  # noqa: E402
from pfcommon.mitgcm_io import read_bin  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402

DEFAULT_ROOT = EXP1_RUNS
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"

WINDOW = (12.0, 24.0)
STAT_TOL = 0.02
RQ_FLOOR, B_FLOOR = 0.03, 2.5  # m (finest dx0)
SIGMA_K = 2.0
R21, R32 = 1.881, 1.912  # r_eff of the note (dx5->dx2p5, dx10->dx5)
G_TOL = 0.2
S_RANGE = (0.75, 1.25)
P_RANGE = (0.5, 3.0)
CG3D_CAP = 2000
CFL_MAX = 0.5
REF = "a_line2d_twolayer_q0.03"
CONV = ("conv_dx10", "conv_dx5", "conv_dx2p5")  # coarse -> fine
DX0 = {"conv_dx10": 10.0, "conv_dx5": 5.0, "conv_dx2p5": 2.5}


# ----------------------------------------------------------------- validity
def parse_monitor(text: str) -> dict:
    """cg3d iteration counts and advective CFLs from MITgcm output.txt monitor lines."""
    return {
        "cg3d": [int(v) for v in re.findall(r"cg3d_iters \(last\)\s*=\s*(\d+)", text)],
        "cfl_u": [float(v) for v in re.findall(r"advcfl_uvel_max\s*=\s*([-+0-9.Ee]+)", text)],
        "cfl_w": [float(v) for v in re.findall(r"advcfl_wvel_max\s*=\s*([-+0-9.Ee]+)", text)],
    }


def check_monitor(mon: dict, n_timesteps: int, cap: int = CG3D_CAP, cfl_max: float = CFL_MAX,
                  require_per_step: bool = True) -> dict:
    """Validity items 2 and 3: exactly nTimeSteps cg3d lines, all < cap; CFL <= 0.5.

    With ``require_per_step`` False (production reference, monitored every 200 steps)
    the lines that exist are checked and ``per_step_verified`` is False.
    """
    reasons = []
    n = len(mon["cg3d"])
    per_step = n == n_timesteps
    if require_per_step and not per_step:
        reasons.append(f"cg3d lines: found {n}, expected nTimeSteps = {n_timesteps}")
    if mon["cg3d"] and max(mon["cg3d"]) >= cap:
        reasons.append(f"cg3d reached the cap: max {max(mon['cg3d'])} >= {cap}")
    cmax = max(mon["cfl_u"] + mon["cfl_w"], default=float("nan"))
    if not mon["cfl_u"] or not mon["cfl_w"]:
        reasons.append("no CFL monitor lines")
    elif cmax > cfl_max:
        reasons.append(f"CFL {cmax:.3f} > {cfl_max}")
    note = None if per_step else (
        f"sampled monitor lines ({n} of {n_timesteps} steps), not per-step verified")
    return {"valid": not reasons, "reasons": reasons, "per_step_verified": per_step, "note": note,
            "n_cg3d_lines": n, "cg3d_max": max(mon["cg3d"], default=None),
            "cg3d_mean": float(np.mean(mon["cg3d"])) if mon["cg3d"] else None, "cfl_max": cmax}


def _cap_from_data(run_dir: Path) -> int:
    f = Path(run_dir) / "data"
    if f.exists():
        m = re.search(r"cg3dMaxIters\s*=\s*(\d+)", f.read_text())
        if m:
            return int(m.group(1))
    return CG3D_CAP


# -------------------------------------------------------- thresholds / trend
def threshold(sigma_a: float, sigma_b: float, floor: float) -> float:
    return max(SIGMA_K * math.hypot(sigma_a, sigma_b), floor)


def trend(coarse: float, fine: float, thr: float) -> str:
    """Label of one refinement step (fine minus coarse): sube / baja / plano."""
    d = fine - coarse
    return "sube" if d > thr else "baja" if d < -thr else "plano"


def is_stationary(r_first: float, r_second: float, r_all: float, tol: float = STAT_TOL) -> bool:
    return bool(abs(r_second - r_first) / abs(r_all) < tol)


# ------------------------------------------------------------------ Richardson
def richardson(phi1: float, phi2: float, phi3: float, r21: float = R21, r32: float = R32,
               max_iter: int = 200, tol: float = 1e-10) -> dict:
    """Variable-ratio Richardson extrapolation (Celik et al. 2008 [verificar]).

    phi1 fine, phi2 medium, phi3 coarse. Inconclusive (``conclusive`` False, no
    R_ext) if s < 0, an eps is zero, the iteration does not converge, or p is
    outside [0.5, 3]. ``p_r2`` (r = 2) is a diagnostic only.
    """
    e21, e32 = phi2 - phi1, phi3 - phi2
    out = {"conclusive": False, "p": None, "R_ext": None, "reason": None, "p_r2": None}
    if e21 == 0 or e32 == 0:
        out["reason"] = "epsilon = 0"
        return out
    ratio = e32 / e21
    out["p_r2"] = abs(math.log(abs(ratio))) / math.log(2.0) if ratio != 0 else None
    s = 1.0 if ratio > 0 else -1.0
    if s < 0:
        out["reason"] = "non-monotone (s < 0)"
        return out
    if abs(ratio) <= 1.0:  # R27: increments do not shrink under refinement
        out["reason"] = "monotone divergence (|eps32/eps21| <= 1)"
        return out
    q, p = 0.0, None
    for _ in range(max_iter):
        p_new = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if p_new <= 0:
            p = None
            break
        try:
            a, b = r21**p_new - s, r32**p_new - s
            q_new = math.log(a / b)
        except (ValueError, ZeroDivisionError, OverflowError):
            p = None
            break
        if abs(p_new - (p if p is not None else 1e9)) < tol:
            p = p_new
            break
        p, q = p_new, q_new
    else:
        p = None
    if p is None:
        out["reason"] = "iteration did not converge"
        return out
    out["p"] = p
    if not (P_RANGE[0] <= p <= P_RANGE[1]):
        out["reason"] = f"p = {p:.3f} outside [{P_RANGE[0]}, {P_RANGE[1]}]"
        return out
    out["conclusive"] = True
    out["R_ext"] = (r21**p * phi1 - phi2) / (r21**p - 1.0)
    return out


# --------------------------------------------------------------- classification
def classify(valid: bool, rq_trend: tuple[str, str], b_trend: tuple[str, str],
             width_fixed: bool, rich: dict) -> str:
    """Outcome class of the convergence study, evaluated in the note's order.

    valid        all conv runs valid and stationary (V otherwise)
    rq_trend     labels of (dx10->dx5, dx5->dx2p5) for R_Q (sube/baja/plano)
    b_trend      same for b
    width_fixed  G: |(b_coarse/b_fine)/r_eff - 1| < 0.2 in BOTH steps
    rich         `richardson` result (uses ``conclusive`` and ``R_ext``)
    """
    if not valid:
        return "V"
    if width_fixed:
        return "G"
    rq_up = all(t == "sube" for t in rq_trend)
    b_down = all(t == "baja" for t in b_trend)
    conclusive = bool(rich.get("conclusive"))
    if rq_up and b_down:
        if conclusive and S_RANGE[0] <= rich["R_ext"] <= S_RANGE[1]:
            return "S"
        if conclusive and rich["R_ext"] < S_RANGE[0]:
            return "P1"
        return "P2"  # inconclusive Richardson or R_ext > 1.25
    if rq_up:
        return "W"
    if all(t == "baja" for t in rq_trend):
        return "D"
    if all(t == "plano" for t in rq_trend):
        return "F"
    return "M"


CLASS_TEXT = {
    "V": "Inconclusive (validity)",
    "G": "Width fixed by the grid: plume not resolved; inconclusive for R_Q",
    "S": "Class (c) attribution to resolution: SUPPORTED",
    "P1": "Partial support: numerical part quantified, remainder unexplained (candidate class b)",
    "P2": "Trend in favour, no extrapolated number",
    "W": "Not supported: R_Q changes without the plume narrowing (entrainment or 2D physics, class b)",
    "D": "Refutes class (c): points to our error (a) or missing physics (b)",
    "F": "Refutes class (c): deficit independent of resolution; evaluate (b)",
    "M": "Inconclusive (non-monotone): no attribution",
}


def width_fixed_flag(b10: float, b5: float, b25: float) -> tuple[bool, tuple[float, float]]:
    """G test: |(b_coarse/b_fine)/r_eff - 1| < 0.2 in both steps.

    A zero, negative or non-finite width gives (False, (nan, nan)); `assemble`
    treats that as an invalid run (class V), never as a width ratio."""
    if not all(math.isfinite(b) and b > 0 for b in (b10, b5, b25)):
        return False, (float("nan"), float("nan"))
    v = ((b10 / b5) / R32, (b5 / b25) / R21)
    return bool(all(abs(x - 1.0) < G_TOL for x in v)), v


# ------------------------------------------------------------ per-run metrics
def _metrics_at(run_dir: Path, window: tuple[float, float], params: dict) -> dict:
    prof = plume_profile(run_dir, window)
    ev = evaluate(prof, params)
    m = ev["modes"]["mitgcm"]
    return {"R_Q": m["u_ratio"], "R_w": ev["post_hoc"]["w_max_ratio"],
            "b": ev["post_hoc"]["width_10pct_lower_half"], "ev": ev}


def metric_windows(t0: float = 12.0) -> dict:
    """Averaging windows (hours) of the note's protocol, starting at ``t0``.

    The protocol (notes/pd2_convergence.md) uses t0 = 12 h: X over 12-24 h,
    stationarity halves 12-18 / 18-24 h, sigma from hourly windows 12-13 ... 23-24 h.
    t0 = 24 h is the pre-registered shift for the extended 36 h `sens_xu2012` run
    (note section "Extensión de sens_xu2012"); the rules themselves are unchanged."""
    return {"X": (t0, t0 + 12.0), "first": (t0, t0 + 6.0), "second": (t0 + 6.0, t0 + 12.0),
            "samples": [(t0 + k, t0 + k + 1.0) for k in range(12)]}


def run_metrics(run_dir: Path, require_per_step: bool = True, t0: float = 12.0) -> dict:
    """Validity, X(t0..t0+12 h), sigma, stationarity, C2/C2'' of one run. Invalid runs
    (failed status, monitor check) return ``valid`` False with reasons."""
    run_dir = Path(run_dir)
    out = {"run": run_dir.name, "valid": True, "reasons": []}
    try:
        require_ok(run_dir)
    except RuntimeError as e:
        return {**out, "valid": False, "reasons": [str(e)]}
    case = json.loads((run_dir / "case.json").read_text())
    params = case["params"]
    n = int(params["nTimeSteps"])
    mon = check_monitor(parse_monitor((run_dir / "output.txt").read_text()), n,
                        cap=_cap_from_data(run_dir), require_per_step=require_per_step)
    out["monitor"] = mon
    if not mon["valid"]:
        out["valid"] = False
        out["reasons"] += mon["reasons"]
    w = metric_windows(t0)
    x = _metrics_at(run_dir, w["X"], params)
    first = _metrics_at(run_dir, w["first"], params)["R_Q"]
    second = _metrics_at(run_dir, w["second"], params)["R_Q"]
    samples = [_metrics_at(run_dir, win, params) for win in w["samples"]]
    out.update({
        "R_Q": x["R_Q"], "R_w": x["R_w"], "b": x["b"],
        "sigma_R_Q": float(np.std([s["R_Q"] for s in samples], ddof=1)),
        "sigma_b": float(np.std([s["b"] for s in samples], ddof=1)),
        "R_Q_12_18": first, "R_Q_18_24": second,
        "stationary": is_stationary(first, second, x["R_Q"]),
        "C2_pass": x["ev"]["pass"]["C2"], "C2pp_pass": x["ev"]["pass"]["C2pp"],
        "C2p_pass_info": x["ev"]["pass"]["C2p"],
        "z_detach_diff_over_H": x["ev"]["modes"]["mitgcm"]["z_detach_diff_over_H"],
        "melt_ratio": x["ev"]["modes"]["mitgcm"]["melt_ratio"],
        "theory_b_lower_half": x["ev"]["post_hoc"]["theory_b_lower_half"],
    })
    if not out["stationary"]:
        out["valid"] = False
        out["reasons"].append(f"not stationary (R_Q second half vs first half of {w['X']} h >= 2 %)")
    # r_eff check: mean dx of cell centres within 0-25 m of the face
    nx = int(params["nx"])
    dx = read_bin(run_dir / "delX.bin", (nx,))
    xc = np.cumsum(dx) - 0.5 * dx
    out["dx_mean_0_25m"] = float(dx[xc <= 25.0].mean())
    return out


# --------------------------------------------------------------------- study
def _diff(a: dict, b: dict) -> dict:
    """Change from `a` (baseline) to `b` with the note's thresholds."""
    tR = threshold(a["sigma_R_Q"], b["sigma_R_Q"], RQ_FLOOR)
    tb = threshold(a["sigma_b"], b["sigma_b"], B_FLOOR)
    return {"dR_Q": b["R_Q"] - a["R_Q"], "thr_R_Q": tR, "R_Q_trend": trend(a["R_Q"], b["R_Q"], tR),
            "db": b["b"] - a["b"], "thr_b": tb, "b_trend": trend(a["b"], b["b"], tb)}


def assemble(metrics: dict[str, dict]) -> dict:
    """Pure assembly of the study result from per-run metrics (see `run_metrics`)."""
    res: dict = {"runs": metrics, "window_h": list(WINDOW)}
    conv = {k: metrics.get(k) for k in CONV}
    missing = [k for k, v in conv.items() if v is None]
    if missing:
        res["convergence"] = {"class": "INCOMPLETE", "reason": f"missing runs: {missing}"}
    else:
        bad = {k: v["reasons"] for k, v in conv.items() if not v["valid"]}
        if bad:  # class V: no steps, Richardson or G (invalid runs may carry no metrics)
            res["convergence"] = {"class": "V", "meaning": CLASS_TEXT["V"], "invalid_runs": bad}
            return _finish(res, metrics)
        c10, c5, c25 = (conv[k] for k in CONV)
        fixed, gv = width_fixed_flag(c10["b"], c5["b"], c25["b"])
        if not all(math.isfinite(v) for v in gv):
            bad = {k: ["zero or non-finite plume width"] for k in CONV}
            res["convergence"] = {"class": "V", "meaning": CLASS_TEXT["V"], "invalid_runs": bad}
            return _finish(res, metrics)
        steps = [_diff(c10, c5), _diff(c5, c25)]
        rich = richardson(c25["R_Q"], c5["R_Q"], c10["R_Q"])
        cls = classify(True, tuple(s["R_Q_trend"] for s in steps),
                       tuple(s["b_trend"] for s in steps), fixed, rich)
        conv_res = {"class": cls, "meaning": CLASS_TEXT[cls], "invalid_runs": bad,
                    "steps": {"dx10_to_dx5": steps[0], "dx5_to_dx2p5": steps[1]},
                    "richardson": rich, "width_fixed_ratios": list(gv), "r_eff": {"r21": R21, "r32": R32}}
        if cls == "P1":
            conv_res["numerical_part_R_ext_minus_RQ_dx5"] = rich["R_ext"] - c5["R_Q"]
            conv_res["numerical_part_caveat"] = (
                "R_Q(dx5) is conv_dx5, which also has dz = 5 m (production 10 m): the numerical part "
                "mixes dx and dz refinement")
        res["convergence"] = conv_res
    return _finish(res, metrics)


def _finish(res: dict, metrics: dict) -> dict:
    """Sensitivities and per-run tolerance flags (shared tail of `assemble`)."""
    ref = metrics.get(REF)
    sens = {}
    if ref is not None and ref["valid"]:
        for name, key in (("dz_sensitivity", "conv_dx5"), ("wall_free_slip", "sens_freeslip"),
                          ("xu2012_closure_joint_with_wall_drag", "sens_xu2012")):
            m = metrics.get(key)
            if m is None or not m["valid"]:
                sens[name] = {"status": "unavailable or invalid", "run": key}
                continue
            d = _diff(ref, m)
            d["run"] = key
            if name == "wall_free_slip":
                d["wall_drag_contributes"] = d["R_Q_trend"] == "sube"
            sens[name] = d
    else:
        sens["status"] = "reference missing or invalid"
    res["sensitivity"] = sens
    res["per_run_tolerances"] = {k: {"C2_pass": v.get("C2_pass"), "C2pp_pass": v.get("C2pp_pass")}
                                 for k, v in metrics.items() if v.get("valid") or "C2_pass" in v}
    return res


XU_EXT = "sens_xu2012_36h"


def evaluate_study(run_dirs: dict[str, Path], out_dir: Path | None = None) -> dict:
    """Evaluate the study; ``run_dirs`` maps names (conv_dx10, conv_dx5, conv_dx2p5,
    sens_freeslip, sens_xu2012, a_line2d_twolayer_q0.03) to run directories.
    Writes ``pd2_convergence.json`` (default exp1_fjord_plumes/figures)."""
    ext = run_dirs.get(XU_EXT)
    metrics = {k: run_metrics(Path(v), require_per_step=(k != REF))
               for k, v in run_dirs.items() if k != XU_EXT}
    if ext is not None and not metrics.get("sens_xu2012", {}).get("valid", False):
        # pre-registered extension (note "Extensión de sens_xu2012"): the 24 h run was not
        # stationary; the identical 36 h run is evaluated over 24-36 h and replaces it.
        m = run_metrics(Path(ext), t0=24.0)
        m["replaces"] = "sens_xu2012 (24 h run not stationary); window 24-36 h"
        metrics["sens_xu2012"] = m
    res = assemble(metrics)
    out_dir = Path(out_dir or FIG_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "pd2_convergence.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    return res


def make_figure(res: dict, stem: Path):
    import matplotlib.pyplot as plt

    from pfcommon import style

    style.apply()
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.0), constrained_layout=True)
    runs = [k for k in CONV if k in res["runs"] and "R_Q" in res["runs"][k]]
    dx = [DX0[k] for k in runs]
    for ax, key, sk, lab in ((axs[0], "R_Q", "sigma_R_Q", "R_Q = mean(Q/b) / theory u"),
                             (axs[1], "b", "sigma_b", "plume width b (m)")):
        ax.errorbar(dx, [res["runs"][k][key] for k in runs], yerr=[2 * res["runs"][k][sk] for k in runs],
                    fmt="o-", capsize=3, label="MITgcm (12-24 h, ±2σ)")
        ax.set_xscale("log")
        ax.set_xlabel("dx0 at the ice face (m)")
        ax.set_xticks(sorted(DX0.values()), [f"{v:g}" for v in sorted(DX0.values())])
        ax.minorticks_off()
        ax.set_ylabel(lab)
        ax.invert_xaxis()
    conv = res.get("convergence", {})
    rich = conv.get("richardson", {})
    if rich.get("conclusive"):
        axs[0].plot([DX0["conv_dx2p5"] / 2], [rich["R_ext"]], "s", color="C1", label="Richardson")
    axs[0].axhspan(0.75, 1.25, color="0.9", zorder=0, label="C2' tolerance")
    if runs and "theory_b_lower_half" in res["runs"][runs[-1]]:
        axs[1].axhline(res["runs"][runs[-1]]["theory_b_lower_half"], color="k", ls="--", label="theory top-hat b")
    axs[0].set_title(f"class {conv.get('class', '?')}", fontsize=8)
    axs[0].legend(fontsize=6.5, frameon=False)
    axs[1].legend(fontsize=6.5, frameon=False)
    sens = res.get("sensitivity", {})
    names = [n for n in ("dz_sensitivity", "wall_free_slip", "xu2012_closure_joint_with_wall_drag")
             if "dR_Q" in sens.get(n, {})]
    if names:
        labels = {"dz_sensitivity": "dz 10 → 5 m", "wall_free_slip": "free slip",
                  "xu2012_closure_joint_with_wall_drag": "Xu-style closure\n(joint with wall drag)"}
        axs[2].barh(range(len(names)), [sens[n]["dR_Q"] for n in names], color="C0")
        for i, n in enumerate(names):  # detection threshold of the note, drawn as a band
            axs[2].plot([-sens[n]["thr_R_Q"], sens[n]["thr_R_Q"]], [i, i], "k|-", lw=1, ms=8,
                        label="± detection threshold" if i == 0 else None)
        axs[2].legend(fontsize=6.5, frameon=False, loc="lower right")
        axs[2].set_yticks(range(len(names)), [labels[n] for n in names], fontsize=7)
        axs[2].axvline(0, color="k", lw=0.6)
    axs[2].set_xlabel("change in R_Q vs production")
    return style.save(fig, stem)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--out-dir", type=Path, default=FIG_DIR)
    ap.add_argument("--single", type=Path, default=None,
                    help="evaluate one run only (e.g. the extended sens_xu2012) and print metrics")
    ap.add_argument("--t0", type=float, default=12.0, help="window start (h) for --single")
    a = ap.parse_args()
    if a.single is not None:
        m = run_metrics(a.single, t0=a.t0)
        print(json.dumps({k: v for k, v in m.items() if k != "monitor"}, indent=1, default=str))
        raise SystemExit(0)
    dirs = {k: a.runs_root / k for k in (*CONV, "sens_freeslip", "sens_xu2012", XU_EXT, REF)
            if (a.runs_root / k).exists()}
    r = evaluate_study(dirs, a.out_dir)
    make_figure(r, a.out_dir / "pd2_convergence")
    print(json.dumps({"class": r["convergence"].get("class")}))
