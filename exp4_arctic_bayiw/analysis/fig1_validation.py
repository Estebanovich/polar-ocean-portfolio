"""Figure 1 of Exp. 4: validation V2-V4 against the frozen PD0 thresholds, from figures/validation.json.

Panels: (a) V2a horizontal kinetic energy spectrum at the shelf station with the near-inertial
band [0.85 f, 1.20 f) and f; (b) V2b front arrival time against along-coast distance with the
fitted speed and the frozen range [0.8 c1(H mean), 1.2 c1(H max)]; (c) V3 M2 pycnocline
amplitude against distance to land with the fitted decay scale (1-8 km); (d) V3 M2 phase along the
south coast with the fitted phase speed (0.21-0.84 m s^-1); (e) V4 flux fraction per resolved band,
F1 (200 m) and V4 (400 m), with V4a. Each panel title gives PASS or FAIL as evaluated by validate.py;
the figure draws numbers, it decides nothing. Band labels in (e): Kelvin = Kelvin/subinertial,
NI = near-inertial, SI = super-inertial (PD0 section 3).

Usage: fig1_validation.py [--validation figures/validation.json] [--out figures/fig1_validation]
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

PKG = Path(__file__).resolve().parents[1]
BANDS = ("kelvin_subinertial", "near_inertial", "superinertial", "seiche")
LABEL = {"kelvin_subinertial": "Kelvin", "near_inertial": "NI", "superinertial": "SI",
         "seiche": "seiche"}


def _tag(ok) -> str:
    return "not evaluated" if ok is None else ("PASS" if ok else "FAIL")


def _missing(ax, name):
    ax.text(0.5, 0.5, f"{name}: not evaluated", ha="center", va="center", transform=ax.transAxes)
    ax.set_xticks([]); ax.set_yticks([])


def figure(val: dict):
    import matplotlib.pyplot as plt
    style.apply()
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 5.2), constrained_layout=True)
    a, b, c, d, e, f = axs.ravel()
    C = style.OKABE_ITO
    # (a) V2a
    v = val.get("V2a")
    if v:
        fr = np.asarray(v["spectrum"]["freq_hz"]) * 3600.0              # cph
        P = np.asarray(v["spectrum"]["hke"]) / 3600.0                    # per cph
        fc = v["f_hz"] * 3600.0
        a.semilogy(fr[1:], P[1:], color=C[0], lw=0.8)
        lo, hi = v["band_f"]
        a.axvspan(lo * fc, hi * fc, color=C[4], alpha=0.25, lw=0)
        a.axvline(fc, color="k", lw=0.7, ls="--")
        a.plot(v["f_peak_hz"] * 3600.0, np.interp(v["f_peak_hz"] * 3600.0, fr, P), "o", color=C[1], ms=4)
        a.set_xlim(0, 0.3)
        a.set_xlabel("frequency (cph)"); a.set_ylabel("HKE (m$^2$ s$^{-2}$ cph$^{-1}$)")
        a.set_title(f"(a) V2a {v['offset_df']:.2f}$\\Delta f$ {_tag(v['pass'])}", loc="left")
    else:
        _missing(a, "V2a")
    # (b) V2b
    v = val.get("V2b")
    if v:
        s = np.asarray(v["s_m"]) / 1e3
        ta = np.asarray(v["arrival_s"], float) / 3600.0
        b.plot(s, ta, "o", color=C[0], ms=4)
        if v.get("c_m_s") not in (None, float("inf")):
            ok = np.isfinite(ta)
            t0 = np.nanmean(ta - s * 1e3 / v["c_m_s"] / 3600.0)
            ss = np.linspace(s.min(), s.max(), 20)
            b.plot(ss, t0 + ss * 1e3 / v["c_m_s"] / 3600.0, color=C[1], lw=1,
                   label=f"c = {v['c_m_s']:.2f} m s$^{{-1}}$")
            for cr, ls in zip(v["range_m_s"], (":", "--")):
                b.plot(ss, np.nanmean(ta[ok]) + (ss - s[ok].mean()) * 1e3 / cr / 3600.0, color="0.5",
                       lw=0.7, ls=ls, label=f"{cr:.2f} m s$^{{-1}}$")
            b.legend(frameon=False, fontsize=7)
        b.set_xlabel("along-coast distance (km)"); b.set_ylabel("front arrival (h)")
        b.set_title(f"(b) V2b {_tag(v['pass'])}", loc="left")
    else:
        _missing(b, "V2b")
    v = val.get("V3")
    if v:
        # (c) V3 presence
        p = v["presence"]
        dd = np.asarray(p["d_m"]) / 1e3
        am = np.asarray(p["amplitude"])
        c.semilogy(dd, am, "o", color=C[0], ms=4)
        L = p["decay_scale_km"]
        if np.isfinite(L) and am.size:
            fit = np.exp(np.polyfit(dd, np.log(am), 1)[1]) * np.exp(-dd / L)
            c.semilogy(dd, fit, color=C[1], lw=1, label=f"L = {L:.1f} km")
            c.legend(frameon=False, fontsize=7)
        c.set_xlabel("distance to land (km)"); c.set_ylabel("M2 amplitude of $\\zeta$ (m)")
        c.set_title(f"(c) V3 trapping {_tag(p['pass'])}", loc="left")
        # (d) V3 sense and speed
        q = v["speed"]
        s = np.asarray(q["s_m"]) / 1e3
        ph = np.rad2deg(np.asarray(q["phase_rad"]))
        d.plot(s, ph, "o", color=C[0], ms=4)
        k = v["sense"]["dphase_ds_rad_per_m"]
        d.plot(s, ph.mean() + np.rad2deg(k) * 1e3 * (s - s.mean()), color=C[1], lw=1,
               label=f"c = {q['c_m_s']:.2f} m s$^{{-1}}$ ({v['sense']['direction']}ward)")
        d.legend(frameon=False, fontsize=7, loc="upper left")
        d.set_xlabel("along the south coast (km)"); d.set_ylabel("M2 phase lag (deg)")
        d.set_title(f"(d) V3 sense {_tag(v['sense']['pass'])}, c {_tag(q['pass'])}", loc="left")
    else:
        _missing(c, "V3"); _missing(d, "V3")
    # (e) V4
    v = val.get("V4")
    if v:
        x = np.arange(len(BANDS))
        p1 = [v["V4b"]["phi_F1"].get(bb, np.nan) for bb in BANDS]
        p4 = [v["V4b"]["phi_V4"].get(bb, np.nan) for bb in BANDS]
        e.bar(x - 0.2, p1, 0.4, color=C[0], label="F1 200 m")
        e.bar(x + 0.2, p4, 0.4, color=C[1], label="V4 400 m")
        e.set_xticks(x, [LABEL[bb] for bb in BANDS], fontsize=8)
        e.set_ylabel("fraction of outward flux")
        e.legend(frameon=False, fontsize=7)
        e.set_title(f"(e) V4a {_tag(v['V4a']['pass'])}, V4b {_tag(v['V4b']['pass'])}", loc="left")
        e.text(0.98, 0.70, f"V4a: {100 * v['V4a']['rel_diff']:.0f} %", ha="right", fontsize=7,
               transform=e.transAxes)
    else:
        _missing(e, "V4")
    # (f) the frozen thresholds
    T = val.get("thresholds", {})
    f.axis("off")
    lines = ["Frozen PD0 thresholds",
             f"V2a: peak within {T.get('V2a', {}).get('max_peak_offset_df')} $\\Delta f$ of f",
             f"V2b: c in {T.get('V2b', {}).get('range_m_s')} m s$^{{-1}}$",
             f"V3: decay {T.get('V3', {}).get('decay_scale_km')} km; cyclonic;",
             f"     c in {T.get('V3', {}).get('c_range_m_s')} m s$^{{-1}}$",
             f"V4a: $\\leq$ {T.get('V4a', {}).get('max_rel_diff')} relative",
             f"V4b: $\\leq$ {T.get('V4b', {}).get('max_abs_diff')} per band"]
    f.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=8, transform=f.transAxes)
    return fig


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--validation", type=Path, default=PKG / "figures" / "validation.json")
    ap.add_argument("--out", type=Path, default=PKG / "figures" / "fig1_validation")
    a = ap.parse_args(argv)
    val = json.loads(a.validation.read_text())
    fig = figure(val)
    style.save(fig, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
