#!/opt/anaconda3/envs/pyocean/bin/python
r"""Figures 2-4 of Exp. 4 (Task 7, Step 3), in the pfcommon.style of the Exp. 2-3 figures.

Usage::

    figures.py [--fig-dir exp4_arctic_bayiw/figures] [--cache EXP4_RUNS/analysis_cache]

* Figure 2 ``fig2_spectra_partition``: (a) surface horizontal KE spectra (one Hann window,
  partition.py station_band_energy) of the fjord (F1: south-coast mean and shelf station) against
  the thesis bay (R1k, bay interior station 'Seiche') and the bridge (R2 bay interior, R3 shelf
  station in front of the absent bay), with f of each case; (b) band shares of the surface KE at
  those stations; (c) signed band shares F_b / sum_5 |F_b| of the net outward Kang & Fringer
  baroclinic pressure work F_bc (PD0 Enmienda 3; straits for the fjord, mouth for the bay). The
  bridge panel is descriptive (spec 7).
* Figure 3 ``fig3_strait_flux_budget``: (a) the KF baroclinic pressure work F_bc into the inner
  fjord through each strait in F1 against time (PD0 Enmienda 3); (b) the total for F1, F1c, F3 and F2; (c, d) the cumulative baroclinic
  budget of F1 and F2 (P_bc, C, F = F_bc + A_bc in, dE_bc and the residual D; the budget needs
  the advection, so F here is pressure work plus advection).
* Figure 4 ``fig4_ice_maps``: depth-integrated time-mean baroclinic KE per resolved band in the
  inner fjord, F2 (landfast ice) above F3 (open water, same winter profile), one log colour scale
  per band.

Every figure writes PNG + PDF through ``pfcommon.style.save`` and its numbers to a JSON of the
same stem. Figure functions take the JSON dicts (and, for Fig. 4, the maps read from the
partition caches by ``main``).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pfcommon import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

FIG_DIR = REPO / "exp4_arctic_bayiw" / "figures"
BANDS = ("low", "kelvin_subinertial", "near_inertial", "superinertial", "seiche")
RESOLVED = BANDS[1:]
BAND_LABELS = {"low": "low (unresolved)", "kelvin_subinertial": "Kelvin/subinertial",
               "near_inertial": "near-inertial", "superinertial": "superinertial", "seiche": "seiche"}
BAND_COLOURS = {"low": "#BBBBBB", "kelvin_subinertial": style.OKABE_ITO[0],
                "near_inertial": style.OKABE_ITO[1], "superinertial": style.OKABE_ITO[2],
                "seiche": style.OKABE_ITO[3]}
RUN_COLOURS = {"F1": style.OKABE_ITO[0], "F1c": style.OKABE_ITO[3], "F2": style.OKABE_ITO[4],
               "F3": style.OKABE_ITO[1], "R1": "#888888", "R1k": style.OKABE_ITO[7],
               "R2": style.OKABE_ITO[2], "R3": style.OKABE_ITO[5]}
F_CPH = {"fjord": 0.0817353, "thesis": 0.0399352}      # theory.json f (78 N; thesis f0)


def _write_json(rec: dict, stem: Path) -> None:
    Path(stem).with_suffix(".json").write_text(json.dumps(rec, indent=1) + "\n")


# ---------------------------------------------------------------------------- Figure 2
SPECTRA = (("F1", "coast_south_mean", "F1 fjord, south coast (mean)", "-"),
           ("F1", "shelf", "F1 fjord, shelf station", "--"),
           ("R1k", "Seiche", "R1k thesis bay, interior", "-"),
           ("R2", "Seiche", "R2 bay at 78 N, interior", "-"),
           ("R3", "Omega", "R3 no bay, shelf (Omega)", "--"))
FLUX_RUNS = (("F1", "F1 straits"), ("F1c", "F1c straits"), ("F3", "F3 straits"),
             ("F2", "F2 straits"), ("R1k", "R1k mouth"), ("R2", "R2 mouth"))


def fig2_spectra_partition(part: dict, out_dir: Path) -> dict:
    style.apply()
    P = part["runs"]
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(10.0, 3.8), gridspec_kw={"width_ratios": [1.3, 1, 1]})
    rec = {"spectra_stations": [], "station_share": {}, "flux_share": {}}
    for run, st, lab, ls in SPECTRA:
        if run not in P or st not in P[run].get("spectra", {}):
            continue
        s = P[run]["spectra"][st]
        f, p = np.asarray(s["freq_cph"]), np.asarray(s["psd_m2_s-2_cph-1"])
        a.loglog(f[1:], p[1:], ls, color=RUN_COLOURS[run], lw=0.9, label=lab)
        rec["spectra_stations"].append([run, st])
    for k, v in F_CPH.items():
        a.axvline(v, color="k", lw=0.6, ls=":" if k == "thesis" else "--")
    a.text(F_CPH["fjord"] * 1.05, 0.98, "f (78°N)", transform=a.get_xaxis_transform(), fontsize=7, va="top")
    a.text(F_CPH["thesis"] * 1.05, 0.90, "f thesis", transform=a.get_xaxis_transform(), fontsize=7, va="top")
    a.set_xlabel("frequency (cph)")
    a.set_ylabel("surface HKE PSD (m$^2$ s$^{-2}$ cph$^{-1}$)")
    a.set_title("(a) surface KE spectra")
    a.legend(fontsize=6.5, loc="lower left", frameon=False)
    # (b) station band shares
    labels = []
    for q, (run, st, lab, _) in enumerate(SPECTRA):
        if run not in P or st not in P[run].get("stations", {}):
            continue
        fr = P[run]["stations"][st]["fraction"]
        left = 0.0
        for bnd in BANDS:
            v = fr.get(bnd) or 0.0
            b.barh(len(labels), v, left=left, color=BAND_COLOURS[bnd], height=0.7,
                   label=BAND_LABELS[bnd] if not labels else None)
            left += v
        rec["station_share"][f"{run}:{st}"] = fr
        labels.append(lab.replace(", ", "\n"))
    b.set_yticks(range(len(labels)), labels, fontsize=6.5)
    b.invert_yaxis()
    b.set_xlim(0, 1)
    b.set_xlabel("share of surface KE")
    b.set_title("(b) surface KE by band")
    b.legend(fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False)
    # (c) flux band shares (signed)
    flabels = []
    for run, lab in FLUX_RUNS:
        if run not in P or not P[run].get("straits"):
            continue
        fb = P[run]["straits"]["total"]["F_bc_out_J"]
        tot = sum(abs(fb[x]) for x in BANDS)
        sh = {x: fb[x] / tot for x in BANDS}
        y = len(flabels)
        for q, bnd in enumerate(BANDS):
            c.barh(y + (q - 2) * 0.15, sh[bnd], height=0.15, color=BAND_COLOURS[bnd])
        rec["flux_share"][run] = {"share_signed": sh, "sum_abs_J": tot}
        flabels.append(lab)
    c.axvline(0, color="k", lw=0.6)
    c.set_yticks(range(len(flabels)), flabels, fontsize=7)
    c.invert_yaxis()
    c.set_xlim(-1, 1)
    c.set_xlabel("$F_b$ / $\\Sigma_5|F_b|$ (+ out)")
    c.set_title("(c) outward $F_{bc}$ (pressure work) by band")
    fig.tight_layout()
    style.save(fig, Path(out_dir) / "fig2_spectra_partition")
    plt.close(fig)
    _write_json(rec, Path(out_dir) / "fig2_spectra_partition")
    return rec


# ---------------------------------------------------------------------------- Figure 3
def fig3_strait_flux_budget(bud: dict, out_dir: Path) -> dict:
    style.apply()
    B = bud["runs"]
    fig, ax = plt.subplots(2, 2, figsize=(9.0, 6.0), sharex=True)
    rec = {"final_cumulative_bc_J": {}, "strait_series_runs": []}
    a = ax[0, 0]
    if "F1" in B:
        s = B["F1"]["series"]
        for q, (sec, v) in enumerate(s["F_bc_in_W"].items()):
            a.plot(s["t_h"], np.asarray(v) / 1e3, color=style.OKABE_ITO[q], lw=0.8, label=sec)
        a.legend(fontsize=7, frameon=False)
    a.axhline(0, color="k", lw=0.5)
    a.set_ylabel("$F_{bc}$ into the fjord (kW)")
    a.set_title("(a) F1, per strait")
    a = ax[0, 1]
    for run in ("F1", "F1c", "F3", "F2"):
        if run not in B:
            continue
        s = B[run]["series"]
        tot = np.sum([np.asarray(v) for v in s["F_bc_in_W"].values()], axis=0)
        a.plot(s["t_h"], tot / 1e3, color=RUN_COLOURS[run], lw=0.8, label=run)
        rec["strait_series_runs"].append(run)
    a.axhline(0, color="k", lw=0.5)
    a.legend(fontsize=7, frameon=False)
    a.set_title("(b) both straits")
    for a, run in ((ax[1, 0], "F1"), (ax[1, 1], "F2")):
        if run not in B:
            continue
        s = B[run]["series"]
        lab = {"P": "$P_{bc}$ (wind)", "C": "$C$ (conversion)", "F": "$F_{bc}+A_{bc}$ in",
               "dE": "$\\Delta E_{bc}$", "D": "$D$ (residual)"}
        for q, k in enumerate(("P", "C", "F", "dE", "D")):
            a.plot(s["t_h"], np.asarray(s["bc_cumulative_J"][k]) / 1e9, color=style.OKABE_ITO[q],
                   lw=0.9, ls="--" if k == "D" else "-", label=lab[k])
        a.axhline(0, color="k", lw=0.5)
        a.set_title(f"({'c' if run == 'F1' else 'd'}) {run}: cumulative baroclinic budget")
        a.set_xlabel("time (h)")
        rec["final_cumulative_bc_J"][run] = B[run]["closure_J"]["bc"]
    ax[1, 0].set_ylabel("energy (GJ)")
    ax[1, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    style.save(fig, Path(out_dir) / "fig3_strait_flux_budget")
    plt.close(fig)
    _write_json(rec, Path(out_dir) / "fig3_strait_flux_budget")
    return rec


# ---------------------------------------------------------------------------- Figure 4
def fig4_ice_maps(maps: dict, out_dir: Path) -> dict:
    style.apply()
    names = list(maps["band_names"])
    runs = ("F2", "F3")
    fig, ax = plt.subplots(2, len(RESOLVED), figsize=(10.0, 3.4), sharex=True, sharey=True,
                           layout="constrained")
    rec = {"box_total_ke_bc_J": {}, "units": "J m^-2 (depth-integrated time-mean KE_bc)"}
    dx = maps["dx_m"]
    allj = np.concatenate([maps[r]["jj"] for r in runs])
    alli = np.concatenate([maps[r]["ii"] for r in runs])
    j0, j1, i0, i1 = allj.min(), allj.max() + 1, alli.min(), alli.max() + 1
    ext = [(maps["x0_m"] + i0 * dx) / 1e3, (maps["x0_m"] + i1 * dx) / 1e3,
           (maps["y0_m"] + j0 * dx) / 1e3, (maps["y0_m"] + j1 * dx) / 1e3]
    for c, bnd in enumerate(RESOLVED):
        q = names.index(bnd)
        vals = np.concatenate([maps[r]["ke_bc_band"][:, q] for r in runs])
        pos = vals[vals > 0]
        hi = np.percentile(pos, 99.5) if pos.size else 1.0
        norm = LogNorm(vmin=hi * 1e-3, vmax=hi)
        for rr, r in enumerate(runs):
            img = np.full((j1 - j0, i1 - i0), np.nan)
            img[maps[r]["jj"] - j0, maps[r]["ii"] - i0] = np.clip(maps[r]["ke_bc_band"][:, q], hi * 1e-3, None)
            im = ax[rr, c].imshow(img, origin="lower", extent=ext, norm=norm, cmap="cmo.amp")
            ax[rr, c].set_aspect("equal")
            if rr == 0:
                ax[rr, c].set_title(BAND_LABELS[bnd])
        fig.colorbar(im, ax=ax[:, c], orientation="horizontal", fraction=0.08, pad=0.02,
                     label="J m$^{-2}$")
    for rr, r in enumerate(runs):
        ax[rr, 0].set_ylabel(f"{r}\n{'landfast ice' if r == 'F2' else 'open water'}\nnorthing (km)")
        area = maps[r].get("area_m2")
        A = np.asarray(area) if area is not None else dx * dx * np.ones(maps[r]["jj"].size)
        rec["box_total_ke_bc_J"][r] = {b: float((maps[r]["ke_bc_band"][:, names.index(b)] * A).sum())
                                       for b in names}
    for c in range(len(RESOLVED)):
        ax[1, c].set_xlabel("easting (km)")
    style.save(fig, Path(out_dir) / "fig4_ice_maps")
    plt.close(fig)
    _write_json(rec, Path(out_dir) / "fig4_ice_maps")
    return rec


def load_maps(cache: Path) -> dict:
    """The F2 and F3 box-energy caches (partition.cache_file 'energy') as Fig. 4 maps."""
    import yaml
    from exp4_arctic_bayiw.analysis import partition as pa
    spec = yaml.safe_load((REPO / "exp4_arctic_bayiw" / "cases.yaml").read_text())
    bb = spec["configs"]["fjord_vmf"]["bbox_utm"]
    m = {"x0_m": float(bb[0]), "y0_m": float(bb[2]), "dx_m": float(spec["configs"]["fjord_vmf"]["dx_m"])}
    for r in ("F2", "F3"):
        e = pa._load(pa.cache_file(cache, "energy", r))
        m["band_names"] = [str(x) for x in e["band_names"]]
        m[r] = {"jj": e["jj"], "ii": e["ii"], "ke_bc_band": e["ke_bc_band"], "area_m2": e["area_m2"]}
    return m


def main(argv=None) -> int:
    from pfcommon.paths import EXP4_RUNS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fig-dir", type=Path, default=FIG_DIR)
    ap.add_argument("--cache", type=Path, default=EXP4_RUNS / "analysis_cache")
    a = ap.parse_args(argv)
    part = json.loads((a.fig_dir / "partition.json").read_text())
    bud = json.loads((a.fig_dir / "budget.json").read_text())
    fig2_spectra_partition(part, a.fig_dir)
    fig3_strait_flux_budget(bud, a.fig_dir)
    fig4_ice_maps(load_maps(a.cache), a.fig_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
