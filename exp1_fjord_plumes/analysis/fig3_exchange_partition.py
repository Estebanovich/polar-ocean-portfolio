"""Figure 3 (Task 7): fjord-shelf exchange partition for the four (c) cases.

Usage: fig3_exchange_partition.py [--runs-root DIR] [--suffix S] [--out-dir DIR]
                                  [--i-section 100]

Cases ``<runs-root>/c_fjord_{plume,wind,shelf,all}<suffix>``; missing ones are
skipped (listed in ``missing``). Every case must be OK: two-stage cases (wind,
all) need ``run_status_final.json`` OK (``c_runs.case_parts``).

* Q_ex(t) = 1/2 sum |u| dA (``analysis/exchange.py``; Jackson et al. 2018
  Q_in/Q_out form) on the cross-fjord section at the U points of column
  ``--i-section`` (default 100: x = 20 km, mid-fjord, outside the OBCS sponge
  i >= 188 and away from the plume column i = 0), from HOURLY snapshots of u
  (ruling R31: |u| is nonlinear, 3-h means bias it low).
* Heat toward the front H = rho0 c_p sum [-ADVx_TH + T_ref mean(u) dA]
  from the 3-h mean of MITgcm's ADVx_TH (degC m^3/s; exact time mean of the
  scheme's own advective flux, incl. u'T'), i.e. rho0 c_p <sum (-u)(T - T_ref) dA>,
  with T_ref = -1.9 degC, the surface freezing point at S = 34.6 psu
  (lam1 S + lam2, MeltParams; a class-(c) choice: with the plume's net
  down-fjord transport Q0 the value depends on T_ref). ``c_runs.exchange_series``.
* ``cross_check`` per case: Q_ex from the 3-h-mean u and heat from the
  products of 3-h-mean u and T, with relative differences to the headline
  values (sub-3-h variability they miss).
* Panels: (a) Q_ex time series; (b) time-mean Q_ex; (c) time-mean heat
  toward the front. W/Rd (``c_runs.w_over_rd``) is annotated.
* C5' per case from the monitor (``c_runs.c5p_check``).
* Control (``c_fjord_control<suffix>``, no forcing; notes/pd1_c_fjord3d.md, pre-registered
  2026-11-06): same metric on the same section gives the spurious exchange floor. Each case
  reports ``Q_ex_over_control``; a case with mean Q_ex < 3x the control mean is flagged
  ``interpretable: false`` (pre-registered rule). The control is drawn grey in (a) and the
  3x level as a dashed line in (b). If the control run is missing, the flag is None.

Writes ``fig3_exchange_partition.{png,pdf}`` and ``fig3_metrics.json`` in
``--out-dir`` (default exp1_fjord_plumes/figures).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from exp1_fjord_plumes.analysis import c_runs  # noqa: E402
from exp1_fjord_plumes.theory.melt import MeltParams  # noqa: E402
from pfcommon import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

DEFAULT_ROOT = Path("/Volumes/LaCie/portfolio_runs/exp1")
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
CASES = ("c_fjord_plume", "c_fjord_wind", "c_fjord_shelf", "c_fjord_all")
LABELS = {"c_fjord_plume": "plume", "c_fjord_wind": "wind", "c_fjord_shelf": "shelf",
          "c_fjord_all": "all"}
CONTROL = "c_fjord_control"
CONTROL_FACTOR = 3.0  # pre-registered: signal < 3x control is not interpreted (pd1_c_fjord3d.md)
I_SECTION = 100
_mp = MeltParams()
T_REF = _mp.lam1 * 34.6 + _mp.lam2      # -1.899 degC


def main(runs_root: Path = DEFAULT_ROOT, out_dir: Path = FIG_DIR, suffix: str = "",
         i_section: int = I_SECTION) -> dict:
    runs_root, out_dir = Path(runs_root), Path(out_dir)
    found = {c: runs_root / f"{c}{suffix}" for c in CASES if (runs_root / f"{c}{suffix}").is_dir()}
    if not found:
        raise FileNotFoundError(f"no c_fjord_* runs with suffix {suffix!r} in {runs_root}")
    series, metrics = {}, {}
    for c, run in found.items():
        ts = c_runs.exchange_series(run, i_section, T_REF)
        series[c] = ts
        qm, hm = float(ts["Q_ex"].mean()), float(ts["heat_W"].mean())
        q3, h3 = float(ts["Q_ex_3h"].mean()), float(ts["heat_W_3h_product"].mean())
        metrics[c] = {"run": str(run), "n_hourly": int(ts["t_hourly_s"].size),
                      "n_3h": int(ts["t_3h_s"].size),
                      "t_first_h": float(ts["t_hourly_s"][0] / 3600),
                      "t_last_h": float(ts["t_hourly_s"][-1] / 3600),
                      "Q_ex_mean_m3_s": qm, "Q_ex_max_m3_s": float(ts["Q_ex"].max()),
                      "heat_toward_front_mean_W": hm,
                      "C5p": c_runs.c5p_check(run),
                      "cross_check": {
                          "Q_ex_mean_from_3h_means_m3_s": q3,
                          "Q_ex_rel_diff": (q3 - qm) / qm if qm else None,
                          "heat_from_mean_products_W": h3,
                          "heat_rel_diff": (h3 - hm) / abs(hm) if hm else None}}
    ctrl_run = runs_root / f"{CONTROL}{suffix}"
    ctrl = None
    if ctrl_run.is_dir():
        cts = c_runs.exchange_series(ctrl_run, i_section, T_REF)
        series[CONTROL] = cts
        ctrl = {"run": str(ctrl_run), "Q_ex_mean_m3_s": float(cts["Q_ex"].mean()),
                "Q_ex_max_m3_s": float(cts["Q_ex"].max()),
                "heat_toward_front_mean_W": float(cts["heat_W"].mean())}
    for c in found:
        if ctrl is None:
            metrics[c]["Q_ex_over_control"], metrics[c]["interpretable"] = None, None
        else:
            r = metrics[c]["Q_ex_mean_m3_s"] / ctrl["Q_ex_mean_m3_s"] if ctrl["Q_ex_mean_m3_s"] else float("inf")
            metrics[c]["Q_ex_over_control"] = r
            metrics[c]["interpretable"] = bool(r >= CONTROL_FACTOR)
    first = next(iter(found.values()))
    wrd = c_runs.w_over_rd(first)
    dx = c_runs.params(first)["params"]["dx"]

    style.apply()
    fig, axs = plt.subplots(1, 3, figsize=(11.0, 3.8),
                            gridspec_kw={"width_ratios": [2.2, 1, 1], "wspace": 0.45})
    names = list(found)
    colors = {c: style.OKABE_ITO[k] for k, c in enumerate(CASES)}
    for c in names:
        axs[0].plot(series[c]["t_hourly_s"] / 86400, series[c]["Q_ex"] / 1e3, color=colors[c],
                    label=LABELS[c], lw=1.4)
    if ctrl is not None:
        axs[0].plot(series[CONTROL]["t_hourly_s"] / 86400, series[CONTROL]["Q_ex"] / 1e3,
                    color="0.55", ls="--", lw=1.0, label="control (no forcing)")
    axs[0].set(xlabel="time (days)", ylabel=r"$Q_{ex}$ (10$^3$ m$^3$ s$^{-1}$)",
               title=f"(a) hourly exchange flux at x = {i_section * dx / 1e3:.0f} km")
    axs[0].legend(frameon=False, ncol=2)
    xb = np.arange(len(names))
    axs[1].bar(xb, [metrics[c]["Q_ex_mean_m3_s"] / 1e3 for c in names],
               color=[colors[c] for c in names])
    axs[1].set(xticks=xb, xticklabels=[LABELS[c] for c in names],
               ylabel=r"mean $Q_{ex}$ (10$^3$ m$^3$ s$^{-1}$)", title="(b) time-mean")
    if ctrl is not None:
        axs[1].axhline(CONTROL_FACTOR * ctrl["Q_ex_mean_m3_s"] / 1e3, color="0.55", ls="--", lw=1.0)
        axs[1].text(xb[-1] + 0.4, CONTROL_FACTOR * ctrl["Q_ex_mean_m3_s"] / 1e3, "3x control",
                    ha="right", va="bottom", fontsize=7, color="0.4")
    axs[2].bar(xb, [metrics[c]["heat_toward_front_mean_W"] / 1e9 for c in names],
               color=[colors[c] for c in names])
    axs[2].set(xticks=xb, xticklabels=[LABELS[c] for c in names],
               ylabel="heat toward front (GW)", title=r"(c) $T_{ref}$ = %.2f °C" % T_REF)
    axs[0].text(0.02, 0.95, f"W/Rd = {wrd['W_over_Rd']:.2f}", transform=axs[0].transAxes,
                va="top", fontsize=9)
    stem = out_dir / "fig3_exchange_partition"
    style.save(fig, stem)
    plt.close(fig)
    out = {"i_section": i_section, "T_ref_degC": T_REF,
           "Q_ex_source": "hourly UVEL snapshots (snap1h), R31",
           "heat_source": "ADVx_TH 3-h means (mean3h), R31: rho0 cp sum(-ADVx + T_ref u dA)", "W_over_Rd": wrd["W_over_Rd"],
           "Rd": wrd, "cases": metrics, "control": ctrl, "control_factor": CONTROL_FACTOR, "missing": [c for c in CASES if c not in found],
           "C5p_all_pass": all(m["C5p"]["pass"] for m in metrics.values())}
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "fig3_metrics.json").write_text(json.dumps(out, indent=2))
    return json.loads(json.dumps(out))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--out-dir", type=Path, default=FIG_DIR)
    ap.add_argument("--i-section", type=int, default=I_SECTION)
    a = ap.parse_args()
    print(json.dumps(main(a.runs_root, a.out_dir, a.suffix, a.i_section), indent=2))
