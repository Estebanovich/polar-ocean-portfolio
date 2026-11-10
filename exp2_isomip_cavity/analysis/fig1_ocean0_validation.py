#!/opt/anaconda3/envs/pyocean/bin/python
r"""Figure 1 of Exp. 2 and the pd0 Ocean1 verdicts (Task 4, phase B).

Usage::

    fig1_ocean0_validation.py [--run /Volumes/LaCie/portfolio_runs/exp2/ocean1_2km_r21]
                              [--out-dir exp2_isomip_cavity/figures]

Mechanical against ``notes/pd0_preregistro.md`` (Ocean1 at 2 km, year 20 = monthly records
229-240):

* V2a whole-shelf mean melt in [4, 12] m a^-1; V2b max |psi(x, z)| in the cavity
  (x < 640 km) in [0.20, 0.25] Sv; V2c cavity max |psi_bt| (y > 50 km, x < 640 km) in
  [0.09, 0.43] Sv; V2d max |psi_bt| for 630 < x < 650 km in [0.14, 0.68] Sv (R13: |psi|).
  Closed intervals, time-mean fields of year 20 (``metrics.isomip_summary``).
* V3: the year-20 mean (m_bar, C_bar) of the monthly whole-shelf melt and circulation
  strength (max - min psi, x < 640 km), per kind (barotropic 'bt', overturning 'moc'), must
  lie in the COM band (without COCO) of Y26 Fig. 13 (read from figure, approx.), linear in
  m_bar between the anchors m = 5 and m = 10. NOT_EVALUABLE if m_bar is outside [4, 12].
  pd0 makes [4, 12] evaluable but gives anchors only at 5 and 10; outside [5, 10] the band
  used for the verdict is the linear extension of the anchor lines (``extend=True``). The
  variant held at the anchor values (``extend=False``) is reported next to it; if the two
  disagree the JSON says so. Informative: the 240-month least-squares fit against the
  published one (bt y = 0.024 x + 0.071, r = 0.89; moc y = 0.02 x + 0.019, r = 0.97).
* N1 Ocean1 form (year 20 vs year 19, < 3 %), N2 (< 1 %), N3 (require_ok, through
  ``metrics.load_run``, which also checks the R21 melt_law and the contiguity of the 240
  monthly records over the 4 restart stages).
* R17 CFL maxima of year 20 (stage-4 monitor blocks with t >= 19 years) and of the whole run;
  R15 supercooling at the end of year 20 and in the year-20 mean state. Not criteria.

Figure panels: (a) Ocean1 year-20 melt map, grounding line and ice front; (b) 20-year series
of the whole-shelf mean melt; (c) V2 metrics against the published ISOMIP+ ranges; (d) the
melt-circulation diagram over the V3 band; (e) Ocean0 deep melt against GammaT at 2 and 4 km
(R21 calibration runs only; AD16 item 17).

Outputs: fig1_ocean.{png,pdf} and fig1_ocean.json (verdicts, numbers, caption).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402
from exp2_isomip_cavity.analysis import pd1_calibration_eval as ce  # noqa: E402
from pfcommon import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

FIG_DIR = REPO / "exp2_isomip_cavity" / "figures"
DEFAULT_RUN = Path("/Volumes/LaCie/portfolio_runs/exp2/ocean1_2km_r21")
YEAR_S = mt.SECONDS_PER_YEAR
CFL_KEYS = ce.CFL_KEYS
CFL_MAX = ce.CFL_MAX

# pd0 V2 (Y26 sect. 4.2-4.3); closed intervals
V2_RANGES = {
    "V2a": {"key": "mean_melt_shelf", "range": (4.0, 12.0), "units": "m a-1",
            "label": "whole-shelf mean melt", "source": "Y26 p. 2068, L1260-1261"},
    "V2b": {"key": "overturning_max", "range": (0.20, 0.25), "units": "Sv",
            "label": "max overturning", "source": "Y26 p. 2066, L1100-1102"},
    "V2c": {"key": "psi_bt_max_cavity", "range": (0.09, 0.43), "units": "Sv",
            "label": "cavity max |psi_bt|", "source": "Y26 p. 2066, L1118 (0.26 +- 0.17)"},
    "V2d": {"key": "psi_bt_max_front", "range": (0.14, 0.68), "units": "Sv",
            "label": "near-front max |psi_bt|", "source": "Y26 p. 2067, L1155 (0.41 +- 0.27)"},
}
# pd0 V3 anchors (Y26 Fig. 13, read from figure, approx.): (m, lo, hi)
V3_ANCHORS = {"bt": ((5.0, 0.17, 0.30), (10.0, 0.25, 0.45)),
              "moc": ((5.0, 0.10, 0.18), (10.0, 0.19, 0.28))}
V3_MELT_RANGE = (4.0, 12.0)
V3_FIT_PUBLISHED = {"bt": {"slope": 0.024, "intercept": 0.071, "r": 0.89},
                    "moc": {"slope": 0.02, "intercept": 0.019, "r": 0.97}}
V3_BAND_LABEL = "read from Yung et al. 2026 Fig. 13 (approx.)"
# PD2 (Esteban, 2026-11-10): V2b stays the pre-registered FAIL. The pre-registration took the
# x < 640 km restriction from Y26 Fig. 13, but the 0.20-0.25 Sv of Y26 L1099-1109 (and the
# Fig. 8 caption, L1182) is the full-domain maximum, including the peak near the northern
# boundary: our error (a) in operationalising the criterion. Reported as a POST HOC row only.
V2B_POSTHOC_LABEL = ("V2b with the Yung 2026 definition (full-domain year-20 max of "
                     "|psi(x,z)|), post hoc")
KIND_LABEL = {"bt": "barotropic", "moc": "overturning"}
KIND_COLOUR = {"bt": style.OKABE_ITO[0], "moc": style.OKABE_ITO[1]}
RES_COLOUR = {"2km": style.OKABE_ITO[0], "4km": style.OKABE_ITO[1]}


# ---- pd0 rules (pure) ----------------------------------------------------------------------

def range_verdict(value: float, rng) -> str:
    return "PASS" if rng[0] <= value <= rng[1] else "FAIL"


def v3_band(m: float, kind: str, extend: bool = True):
    """(lo, hi) of the V3 band [Sv] at melt m [m a^-1]: linear in m through the anchors
    m = 5 and 10 (pd0). ``extend``: the anchor lines continue to [4, 12]; otherwise the band
    is held at the anchor values outside [5, 10]. None if m is outside [4, 12]."""
    if kind not in V3_ANCHORS:
        raise ValueError(f"kind must be 'bt' or 'moc', not {kind!r}")
    if not V3_MELT_RANGE[0] <= m <= V3_MELT_RANGE[1]:
        return None
    (m0, lo0, hi0), (m1, lo1, hi1) = V3_ANCHORS[kind]
    mm = m if extend else min(max(m, m0), m1)
    f = (mm - m0) / (m1 - m0)
    return (lo0 + f * (lo1 - lo0), hi0 + f * (hi1 - hi0))


def _margin(c: float, band) -> dict:
    """Distance of c outside the band [Sv] and relative to the violated edge (0 inside)."""
    lo, hi = band
    if c < lo:
        return {"abs_Sv": lo - c, "rel": (lo - c) / lo, "side": "below"}
    if c > hi:
        return {"abs_Sv": c - hi, "rel": (c - hi) / hi, "side": "above"}
    return {"abs_Sv": 0.0, "rel": 0.0, "side": "inside"}


def v2b_full_domain_posthoc(psi, xg, rf, jpl_value: float | None = None) -> dict:
    """POST HOC (PD2, 2026-11-10): max |psi(x, z)| over the whole domain (Y26 L1099-1109,
    Fig. 8), with its location, against the V2b range. Not the pd0 verdict."""
    psi = np.asarray(psi, float)
    k, i = np.unravel_index(np.nanargmax(np.abs(psi)), psi.shape)
    val = float(abs(psi[k, i]))
    rng = V2_RANGES["V2b"]["range"]
    return {"label": V2B_POSTHOC_LABEL, "post_hoc": True, "value": val,
            "x_km": float(xg[i]) / 1e3, "z_m": float(rf[k]), "range": list(rng),
            "verdict_against_range": range_verdict(val, rng), "jpl_value": jpl_value,
            "note": "pd0 V2b verdict unchanged (FAIL, cavity x < 640 km); the restriction was "
                    "taken from Y26 Fig. 13 by mistake, class (a) in the operationalisation"}


def v3_verdict(m_bar: float, c_bar: float, kind: str) -> dict:
    if kind not in V3_ANCHORS:
        raise ValueError(f"kind must be 'bt' or 'moc', not {kind!r}")
    band = v3_band(m_bar, kind)
    if band is None:
        return {"verdict": "NOT_EVALUABLE", "m_bar": m_bar, "C_bar": c_bar, "band": None,
                "reason": f"m_bar outside {list(V3_MELT_RANGE)} (pd0)"}
    strict = v3_band(m_bar, kind, extend=False)
    v = range_verdict(c_bar, band)
    return {"verdict": v, "m_bar": m_bar, "C_bar": c_bar, "band": list(band),
            "band_held_at_anchors": list(strict),
            "verdict_band_held_at_anchors": range_verdict(c_bar, strict),
            "variants_agree": v == range_verdict(c_bar, strict),
            "margin_extended": _margin(c_bar, band),
            "margin_held_at_anchors": _margin(c_bar, strict),
            "m_bar_in_extrapolated_band": not (V3_ANCHORS[kind][0][0] <= m_bar
                                               <= V3_ANCHORS[kind][1][0])}


def linear_fit(x, y) -> dict:
    x, y = np.asarray(x, float), np.asarray(y, float)
    slope, intercept = np.polyfit(x, y, 1)
    return {"slope": float(slope), "intercept": float(intercept),
            "r": float(np.corrcoef(x, y)[0, 1]), "n": int(x.size)}


def ocean1_verdicts(summary: dict) -> dict:
    """pd0 V2a-d, V3, N1, N2 from ``metrics.isomip_summary`` of the 240-record Ocean1 run."""
    out = {}
    for name, spec in V2_RANGES.items():
        val = float(summary[spec["key"]]["value"])
        out[name] = {"value": val, "range": list(spec["range"]), "units": spec["units"],
                     "verdict": range_verdict(val, spec["range"]), "source": spec["source"],
                     "quantity": spec["label"]}
    cs = summary["circulation_strength"]
    mon = cs["monthly"]
    out["V3"] = {}
    for kind in ("bt", "moc"):
        v = v3_verdict(float(cs["window_mean"]["melt_shelf"]), float(cs["window_mean"][kind]),
                       kind)
        v["fit_240_months"] = linear_fit(mon["melt_shelf"], mon[kind])
        v["fit_published"] = V3_FIT_PUBLISHED[kind]
        v["band_source"] = V3_BAND_LABEL
        out["V3"][kind] = v
    out["N1"] = {"value": summary["N1"]["value"], "pass": bool(summary["N1"]["pass"]),
                 "form": "ocean1: |<m>(year 20) - <m>(year 19)| / <m>(year 20) < 3 %"}
    out["N2"] = {"value": summary["N2"]["value"], "pass": bool(summary["N2"]["pass"])}
    return out


def cfl_maxima_window(log_text: str, t_min: float) -> dict:
    """R17 CFL maxima over the %MON blocks whose time_secondsf >= t_min."""
    blocks = re.split(r"(?=%MON time_secondsf)", log_text)
    vals = {k: [] for k in CFL_KEYS}
    n = 0
    for b in blocks:
        mt_ = re.match(r"%MON time_secondsf\s*=\s*(\S+)", b)
        if not mt_ or float(mt_.group(1).replace("D", "E")) < t_min:
            continue
        n += 1
        for k in CFL_KEYS:
            vals[k] += [abs(float(v.replace("D", "E")))
                        for v in re.findall(rf"%MON {k}\s*=\s*(\S+)", b)]
    if not n or any(not v for v in vals.values()):
        raise ValueError(f"no complete %MON blocks with t >= {t_min:g} s")
    out = {k: max(v) for k, v in vals.items()}
    out["max"] = max(out[k] for k in CFL_KEYS)
    out["n_monitor_blocks"] = n
    out["above_0.5"] = bool(out["max"] > CFL_MAX)
    return out


# ---- figure data (pure) --------------------------------------------------------------------

def figure_data(summary: dict, verdicts: dict, melt_map, shelf, wet, x, y,
                cal: dict, ver: dict) -> dict:
    """Plain data of every panel (also written to the JSON, and what the tests check)."""
    melt_map = np.asarray(melt_map, float)
    a = {"melt": np.where(shelf, melt_map, np.nan).tolist(),
         "shelf": np.asarray(shelf, bool).tolist(), "wet": np.asarray(wet, bool).tolist(),
         "x_km": (np.asarray(x) / 1e3).tolist(), "y_km": (np.asarray(y) / 1e3).tolist()}
    mon = summary["circulation_strength"]["monthly"]
    n = len(mon["melt_shelf"])
    b = {"years": [(k + 0.5) / 12 for k in range(n)], "melt": list(mon["melt_shelf"]),
         "year20": [19.0, 20.0], "year19": [18.0, 19.0], "N1": verdicts["N1"]}
    c = []
    for name, spec in V2_RANGES.items():
        v = verdicts[name]
        lo, hi = spec["range"]
        c.append({"name": name, "label": spec["label"], "value": v["value"], "range": [lo, hi],
                  "units": spec["units"], "verdict": v["verdict"],
                  "position": (v["value"] - lo) / (hi - lo)})
    d = {"band_label": V3_BAND_LABEL}
    m_band = list(np.linspace(*V3_MELT_RANGE, 33))
    for kind in ("bt", "moc"):
        bands = [v3_band(m, kind) for m in m_band]
        vv = verdicts["V3"][kind]
        d[kind] = {"monthly_melt": list(mon["melt_shelf"]), "monthly_C": list(mon[kind]),
                   "year20_idx": [n - 12, n], "m_bar": vv["m_bar"], "C_bar": vv["C_bar"],
                   "verdict": vv["verdict"], "band_m": m_band,
                   "band_extrapolated": [not (V3_ANCHORS[kind][0][0] <= m
                                              <= V3_ANCHORS[kind][1][0]) for m in m_band],
                   "m_bar_in_extrapolated_band": vv["m_bar_in_extrapolated_band"],
                   "band_lo": [bb[0] for bb in bands], "band_hi": [bb[1] for bb in bands],
                   "fit_published": V3_FIT_PUBLISHED[kind], "fit_ours": vv["fit_240_months"]}
    e = {"k1_band": [28.0, 32.0], "target": 30.0}
    for res in ("2km", "4km"):
        pts = sorted(cal[res]["points"], key=lambda p: p["gammaT"])
        vp = ver[res]["point"]
        e[res] = {"gammaT": [p["gammaT"] for p in pts],
                  "deep_melt": [p["deep_melt_last6_m_per_a"] for p in pts],
                  "gammaT_cal": cal[res]["gammaT_cal"],
                  "verification": [vp["gammaT"], vp["deep_melt_last6_m_per_a"]]}
    c_post = verdicts.get("V2b_posthoc_full_domain")
    return {"a": a, "b": b, "c": c, "c_posthoc": c_post, "d": d, "e": e}


# ---- caption and figure --------------------------------------------------------------------

def caption(verdicts: dict, resolution: dict | None) -> str:
    v = verdicts
    v3 = v["V3"]
    res = resolution or {}
    psi = res.get("psi_bt_max_cavity", {}).get("rel_diff_4_vs_2")
    r = res.get("melt_map_correlation", {}).get("r")
    res_txt = ("Ocean0 2 vs 4 km (descriptive, pd0): the deep and whole-shelf melt agree "
               "largely BY CONSTRUCTION, because GammaT is calibrated separately at each "
               "resolution to the same Ocean0 deep melt (pd0 K1), so that agreement is not "
               "independent evidence of convergence")
    if psi is not None:
        res_txt += f"; the cavity |psi_bt| maximum differs by {100 * psi:+.0f} % (4 vs 2 km)"
    if r is not None:
        res_txt += f" and the melt maps correlate at r = {r:.2f}"
    return (
        "Exp. 2 ISOMIP+ validation at 2 km. Ocean1 (COLD start, WARM restoring, 20 years, "
        "GammaT_cal,2km = 0.030959): (a) melt rate averaged over year 20 (records 229-240), "
        "grounding line (black) and ice-shelf edge (grey dashed); (b) whole-shelf mean melt, "
        "monthly, 20 years (year 20 shaded; N1 year 20 vs 19); (c) pd0 V2 metrics of year 20 "
        "against the published ISOMIP+ ranges (Yung et al. 2026), as position within each "
        f"range: V2a {v['V2a']['verdict']}, V2b {v['V2b']['verdict']} (pd0, cavity x < 640 km"
        + (f"; POST HOC with the Yung 2026 full-domain definition, open diamond: "
           f"{v['V2b_posthoc_full_domain']['value']:.3f} Sv, inside 0.20-0.25"
           if v.get("V2b_posthoc_full_domain") else "") + "), V2c "
        f"{v['V2c']['verdict']}, V2d {v['V2d']['verdict']}; (d) monthly melt vs circulation "
        "strength (max - min psi, x < 640 km) over the multi-model COM band "
        f"({V3_BAND_LABEL}; hatched where extrapolated beyond the 5 and 10 m/a anchors); "
        f"V3 barotropic {v3['bt']['verdict']}, overturning {v3['moc']['verdict']}"
        + (f" (marginal, unresolved within the figure-reading uncertainty: "
           f"{100 * v3['moc']['margin_extended']['rel']:.1f} % below the extended band, "
           f"{100 * v3['moc']['margin_held_at_anchors']['rel']:.1f} % below the band held at "
           f"the 10 m/a anchor)" if v3["moc"]["verdict"] == "FAIL" else "")
        + "; dashed: published fits. (e) Ocean0 deep melt (ice base "
        "< -300 m, months 19-24) against GammaT, R21 calibration runs (circles) and the "
        "verification run at GammaT_cal (star); grey band: pd0 K1 28-32 m/a. "
        + res_txt + ". All verdicts mechanical against the pre-registration (pd0); FAILs "
        "shown as such.")


def make_figure(data: dict, verdicts: dict, out_stem: Path, cap: str) -> list:
    style.apply()
    fig = plt.figure(figsize=(7.0, 9.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.62, 1, 1], hspace=0.5, wspace=0.32,
                          left=0.09, right=0.97, top=0.96, bottom=0.17)
    axA = fig.add_subplot(gs[0, :])
    axB, axC = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    axD, axE = fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])

    # (a) melt map
    a = data["a"]
    x, y = np.array(a["x_km"]), np.array(a["y_km"])
    melt = np.ma.masked_invalid(np.array(a["melt"], float))
    pm = axA.pcolormesh(x, y, melt, cmap="cmo.amp", shading="nearest", vmin=0)
    fig.colorbar(pm, ax=axA, label="melt (m a$^{-1}$)", fraction=0.03, pad=0.01)
    wet, shelf = np.array(a["wet"], float), np.array(a["shelf"], float)
    # grounding line = edge of the wet domain (grounded ice is dry); ice-shelf edge = edge of
    # (shelf or grounded), i.e. where floating ice meets ice-free open water
    axA.contour(x, y, wet, levels=[0.5], colors="k", linewidths=0.9)
    axA.contour(x, y, np.maximum(shelf, 1.0 - wet), levels=[0.5], colors="0.45",
                linewidths=0.8, linestyles="--")
    axA.set(xlabel="x (km)", ylabel="y (km)", xlim=(420.0, x.max()), aspect="equal",
            title="(a) Ocean1 year-20 melt; grounding line (black), ice-shelf edge (dashed)")

    # (b) series
    b = data["b"]
    axB.plot(b["years"], b["melt"], "-", color="k", lw=0.9)
    axB.axvspan(*b["year20"], color=style.OKABE_ITO[3], alpha=0.25, lw=0, label="year 20")
    axB.axvspan(*b["year19"], color="0.85", lw=0, label="year 19")
    n1 = b["N1"]
    axB.set(xlabel="time (years)", ylabel=r"$\langle m \rangle$ (m a$^{-1}$)", xlim=(0, 20),
            title=f"(b) whole-shelf mean melt; N1 = {100 * n1['value']:.2f} % "
                  f"({'PASS' if n1['pass'] else 'FAIL'})")
    axB.legend(fontsize=7, frameon=False, loc="lower right")

    # (c) V2
    rows = data["c"]
    for i, r in enumerate(rows):
        yy = len(rows) - 1 - i
        axC.plot([0, 1], [yy, yy], "-", color="0.75", lw=6, solid_capstyle="butt")
        col = style.OKABE_ITO[2] if r["verdict"] == "PASS" else style.OKABE_ITO[1]
        axC.plot(r["position"], yy, "o" if r["verdict"] == "PASS" else "X", color=col, ms=8)
        axC.text(-0.08, yy + 0.22, f"{r['name']} {r['label']}", fontsize=6.5, ha="left")
        axC.text(1.08, yy, f"{r['value']:.3g} {r['units']}\n[{r['range'][0]:g}, "
                 f"{r['range'][1]:g}] {r['verdict']}", fontsize=6.3, va="center", color=col)
    ph = data.get("c_posthoc")
    if ph:
        lo, hi = ph["range"]
        yy = len(rows) - 1 - [r["name"] for r in rows].index("V2b")
        axC.plot((ph["value"] - lo) / (hi - lo), yy - 0.28, "D", mfc="w",
                 mec=style.OKABE_ITO[4], mew=1.3, ms=6)
        txt = f"POST HOC, Y26 full-domain\ndef.: {ph['value']:.3f} Sv"
        if ph.get("jpl_value"):
            txt += f" (JPL {ph['jpl_value']:.3f})"
        axC.text(1.08, yy - 0.42, txt, fontsize=5.6, va="center", color=style.OKABE_ITO[4])
    axC.set(xlim=(-0.6, 1.9), ylim=(-0.7, len(rows) - 0.3), yticks=[],
            xticks=[0, 1], xticklabels=["low", "high"],
            title="(c) V2, year 20: position in published range")
    axC.axvline(0, color="0.6", lw=0.5)
    axC.axvline(1, color="0.6", lw=0.5)

    # (d) melt-circulation
    d = data["d"]
    for kind in ("bt", "moc"):
        k = d[kind]
        col = KIND_COLOUR[kind]
        bm, ex = np.array(k["band_m"]), np.array(k["band_extrapolated"])
        blo, bhi = np.array(k["band_lo"]), np.array(k["band_hi"])
        inner = ~ex
        axD.fill_between(bm[inner], blo[inner], bhi[inner], color=col, alpha=0.18, lw=0)
        for seg in (bm <= 5.0, bm >= 10.0):          # extrapolated parts (edges shared)
            axD.fill_between(bm[seg], blo[seg], bhi[seg], facecolor="none", edgecolor=col,
                             hatch="////", alpha=0.45, lw=0)
        i0, i1 = k["year20_idx"]
        axD.plot(k["monthly_melt"][:i0], k["monthly_C"][:i0], ".", color=col, ms=2, alpha=0.35)
        axD.plot(k["monthly_melt"][i0:i1], k["monthly_C"][i0:i1], ".", color=col, ms=4)
        mm = np.array(k["band_m"])
        fp = k["fit_published"]
        axD.plot(mm, fp["slope"] * mm + fp["intercept"], "--", color=col, lw=0.8)
        axD.plot(k["m_bar"], k["C_bar"], "o" if k["verdict"] == "PASS" else "X", mfc="w",
                 mec=col, mew=1.5, ms=8,
                 label=f"{KIND_LABEL[kind]} year-20 mean: V3 {k['verdict']}")
    axD.plot([], [], "s", color="0.6", alpha=0.4, label=f"COM band, {d['band_label']}")
    axD.fill_between([], [], [], facecolor="none", edgecolor="0.4", hatch="////",
                     label="band extrapolated beyond the 5 and 10 m/a anchors")
    axD.set_ylim(top=0.66)                        # room for the legend above the band
    axD.set(xlabel=r"melt (m a$^{-1}$)", ylabel="circulation strength (Sv)",
            title="(d) melt vs circulation (V3)")
    axD.legend(fontsize=5.8, frameon=False, loc="upper left")

    # (e) GammaT
    e = data["e"]
    axE.axhspan(*e["k1_band"], color="0.88", lw=0)
    for res in ("2km", "4km"):
        r = e[res]
        col = RES_COLOUR[res]
        axE.plot(r["gammaT"], r["deep_melt"], "o-", color=col, lw=0.9, ms=4,
                 label=f"{res}: GammaT_cal = {r['gammaT_cal']:.5g}")
        axE.plot(*r["verification"], "*", color=col, ms=10, mec="k", mew=0.4)
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    axE.set(xscale="log", yscale="log")
    axE.xaxis.set_major_locator(FixedLocator([0.015, 0.03, 0.06]))
    axE.yaxis.set_major_locator(FixedLocator([10, 15, 20, 30, 40, 60]))
    for ax_ in (axE.xaxis, axE.yaxis):
        ax_.set_minor_locator(NullLocator())
        ax_.set_minor_formatter(NullFormatter())
    axE.set_xticklabels(["0.015", "0.03", "0.06"])
    axE.set_yticklabels(["10", "15", "20", "30", "40", "60"])
    axE.set(xlabel=r"$\Gamma_T$",
            ylabel="Ocean0 deep melt (m a$^{-1}$)",
            title="(e) Ocean0 melt vs GammaT (R21 calibration)")
    axE.legend(fontsize=6.5, frameon=False, loc="upper left")
    fig.text(0.09, 0.005, cap, ha="left", va="bottom", fontsize=6.0, wrap=True, color="0.25")
    paths = style.save(fig, out_stem)
    plt.close(fig)
    return [str(p) for p in paths]


# ---- run -----------------------------------------------------------------------------------

def seam_check(rd) -> dict:
    """Records at each restart seam: the last record of stage k ends where the first of
    stage k + 1 starts (load_run already refuses gaps; this reports them)."""
    seams = []
    for a, b in zip(rd.records[:-1], rd.records[1:]):
        if a.stage_dir != b.stage_dir:
            seams.append({"record_1based": rd.records.index(b) + 1, "t_end_prev_s": a.t1,
                          "t_start_next_s": b.t0, "contiguous": abs(a.t1 - b.t0) < 1e-6,
                          "stage": str(b.stage_dir.name)})
    return {"n_records": rd.n_records, "seams": seams,
            "all_contiguous": all(s["contiguous"] for s in seams)}


def evaluate(run_dir: Path) -> tuple[dict, dict, object]:
    rd = mt.load_run(run_dir)                     # require_ok (N3) + melt_law + contiguity
    summary = mt.isomip_summary(rd, write=False)
    v = ocean1_verdicts(summary)
    v["N3"] = {"pass": True, "source": "metrics.load_run -> require_ok (all 4 stages)"}
    v["records"] = seam_check(rd)
    win = list(mt.calendar_window(rd, "year20"))
    t20 = rd.t_start + 19 * YEAR_S
    logs = [ce.run_log(sd).read_text(errors="replace") for sd in rd.stage_dirs]
    v["R17_cfl_year20"] = cfl_maxima_window(logs[-1], t20)
    v["R17_cfl_whole_run"] = ce.cfl_maxima("\n".join(logs))
    nyw = int(rd.params["ny_wet_domain"])
    itr = int(rd.params["nIter0"]) + int(rd.params["nTimeSteps"])
    T = mt._read_state(rd.stage_dirs[-1], "T", itr)
    S = mt._read_state(rd.stage_dirs[-1], "S", itr)
    v["R15_supercooling_end_state"] = ce.supercooling_stats(T, S, rd.hfac, rd.rc,
                                                            wall_rows=(0, nyw - 1))
    v["R15_supercooling_year20_mean"] = ce.supercooling_stats(
        rd.mean_field("dynDiag", "THETA", win), rd.mean_field("dynDiag", "SALT", win),
        rd.hfac, rd.rc, wall_rows=(0, nyw - 1))
    v["ustar_check_last_record (informative)"] = summary["ustar_check_last_record (informative)"]
    v["mean_melt_deep300_year20"] = summary["mean_melt_deep300"]["value"]
    jpl_val = None
    jf = FIG_DIR / "jpl_comparison.json"
    if jf.exists():
        jpl_val = (json.loads(jf.read_text()).get("ocean1", {}).get("jpl", {})
                   .get("circulation", {}).get("psi_moc_max_abs_full_domain_Sv"))
    v["V2b_posthoc_full_domain"] = v2b_full_domain_posthoc(
        mt.overturning_streamfunction(rd, win), rd.xg[0], rd.rf, jpl_value=jpl_val)
    return v, summary, rd


def run(run_dir: Path, out_dir: Path) -> dict:
    out_dir = Path(out_dir)
    v, summary, rd = evaluate(run_dir)
    win = list(mt.calendar_window(rd, "year20"))
    melt_map = mt.melt_rate_m_per_a(rd.mean_field(mt._SHI, "SHIfwFlx", win))
    nyw = int(rd.params["ny_wet_domain"])
    cal = {r: json.loads((FIG_DIR / f"calibration_{r}_r21.json").read_text())
           for r in ("2km", "4km")}
    ver = {r: json.loads((FIG_DIR / f"verify_{r}_r21.json").read_text()) for r in ("2km", "4km")}
    res = json.loads((FIG_DIR / "pd1_calibration_eval_r21.json").read_text())["resolution_test"]
    data = figure_data(summary, v, melt_map[:nyw], rd.shelf[:nyw], rd.wet[:nyw],
                       rd.xc[0], rd.yc[:nyw, 0], cal, ver)
    cap = caption(v, res)
    paths = make_figure(data, v, out_dir / "fig1_ocean", cap)
    out = {"run": str(run_dir), "window": "Ocean1 year 20, monthly records 229-240 (pd0)",
           "gammaT": float(rd.params["gammaT"]), "verdicts": v,
           "resolution_test_ocean0": res, "caption": cap, "figure": paths,
           "panel_data_summary": {"d": {k: {kk: data["d"][k][kk] for kk in
                                             ("m_bar", "C_bar", "verdict", "fit_ours")}
                                        for k in ("bt", "moc")}, "e": data["e"],
                                  "c": data["c"]}}
    (out_dir / "fig1_ocean.json").write_text(json.dumps(out, indent=2, default=float) + "\n")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Exp. 2 Figure 1 and pd0 Ocean1 verdicts")
    ap.add_argument("--run", type=Path, default=DEFAULT_RUN)
    ap.add_argument("--out-dir", type=Path, default=FIG_DIR)
    a = ap.parse_args(argv)
    out = run(a.run, a.out_dir)
    v = out["verdicts"]
    print(json.dumps({k: v[k]["verdict"] for k in V2_RANGES}
                     | {f"V3_{k}": v["V3"][k]["verdict"] for k in ("bt", "moc")}
                     | {"N1": v["N1"]["pass"], "N2": v["N2"]["pass"], "N3": v["N3"]["pass"]},
                     indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
