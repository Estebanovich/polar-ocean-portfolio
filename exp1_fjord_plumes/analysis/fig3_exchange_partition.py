"""Figure 3 (Task 7): fjord-shelf exchange partition for the four (c) cases.

Usage: fig3_exchange_partition.py [--runs-root DIR] [--suffix S] [--out-dir DIR]
                                  [--i-section 100]

Cases ``<runs-root>/c_fjord_{plume,wind,shelf,all}<suffix>``; missing ones are
skipped (listed in ``missing``). Every case must be OK: two-stage cases (wind,
all) need ``run_status_final.json`` OK (``c_runs.case_parts``).

* Q_ex(t) = 1/2 sum |u| dA (``analysis/exchange.py``; Jackson et al. 2018
  Q_in/Q_out form) on the cross-fjord section at the U points of column
  ``--i-section`` (default 100: x = 20 km, mid-fjord, outside the OBCS sponge
  i >= 188 and away from the plume column i = 0), from the 3-h means.
* Heat toward the front H = rho0 c_p sum (-u)(T - T_ref) dA with
  T_ref = -1.9 degC, the surface freezing point at S = 34.6 psu
  (lam1 S + lam2, MeltParams; a class-(c) choice: with the plume's net
  down-fjord transport Q0 the value depends on T_ref).
* Panels: (a) Q_ex time series; (b) time-mean Q_ex; (c) time-mean heat
  toward the front. W/Rd (``c_runs.w_over_rd``) is annotated.
* C5' per case from the monitor (``c_runs.c5p_check``).

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
        metrics[c] = {"run": str(run), "n_records": int(ts["t_s"].size),
                      "t_first_h": float(ts["t_s"][0] / 3600), "t_last_h": float(ts["t_s"][-1] / 3600),
                      "Q_ex_mean_m3_s": float(ts["Q_ex"].mean()),
                      "Q_ex_max_m3_s": float(ts["Q_ex"].max()),
                      "heat_toward_front_mean_W": float(ts["heat_W"].mean()),
                      "C5p": c_runs.c5p_check(run)}
    first = next(iter(found.values()))
    wrd = c_runs.w_over_rd(first)
    dx = c_runs.params(first)["params"]["dx"]

    style.apply()
    fig, axs = plt.subplots(1, 3, figsize=(11.0, 3.8),
                            gridspec_kw={"width_ratios": [2.2, 1, 1], "wspace": 0.45})
    names = list(found)
    colors = {c: style.OKABE_ITO[k] for k, c in enumerate(CASES)}
    for c in names:
        axs[0].plot(series[c]["t_s"] / 86400, series[c]["Q_ex"] / 1e3, color=colors[c],
                    label=LABELS[c], lw=1.4)
    axs[0].set(xlabel="time (days)", ylabel=r"$Q_{ex}$ (10$^3$ m$^3$ s$^{-1}$)",
               title=f"(a) exchange flux at x = {i_section * dx / 1e3:.0f} km")
    axs[0].legend(frameon=False, ncol=2)
    xb = np.arange(len(names))
    axs[1].bar(xb, [metrics[c]["Q_ex_mean_m3_s"] / 1e3 for c in names],
               color=[colors[c] for c in names])
    axs[1].set(xticks=xb, xticklabels=[LABELS[c] for c in names],
               ylabel=r"mean $Q_{ex}$ (10$^3$ m$^3$ s$^{-1}$)", title="(b) time-mean")
    axs[2].bar(xb, [metrics[c]["heat_toward_front_mean_W"] / 1e9 for c in names],
               color=[colors[c] for c in names])
    axs[2].set(xticks=xb, xticklabels=[LABELS[c] for c in names],
               ylabel="heat toward front (GW)", title=r"(c) $T_{ref}$ = %.2f °C" % T_REF)
    axs[0].text(0.02, 0.95, f"W/Rd = {wrd['W_over_Rd']:.2f}", transform=axs[0].transAxes,
                va="top", fontsize=9)
    stem = out_dir / "fig3_exchange_partition"
    style.save(fig, stem)
    plt.close(fig)
    out = {"i_section": i_section, "T_ref_degC": T_REF, "W_over_Rd": wrd["W_over_Rd"],
           "Rd": wrd, "cases": metrics, "missing": [c for c in CASES if c not in found],
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
