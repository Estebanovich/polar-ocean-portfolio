"""Figures 0-5 of the Exp. 6 README (plan Task 10; spec Sec. 8), from the committed figures/*.json only:
figdata.json (analysis/figdata.py), hypotheses.json (analysis/hypotheses.py) and crosscheck.json
(fine model, plan Task 9). A figure whose JSON is missing is skipped, so the coarse README can be
rendered before the fine runs exist. Flat article style (pfcommon.style).

Units: flux and S in mol C m-2 yr-1 (+ into the ocean); H6b terms in Gmol C yr-1 (1e9 mol C yr-1);
DIC in umol kg-1; G200 in umol kg-1 m-1.
Usage: python -m exp6_southern_carbon.analysis.figures"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pfcommon import style

style.apply()
import matplotlib.pyplot as plt  # noqa: E402

FIG = Path(__file__).resolve().parents[1] / "figures"
C = style.OKABE_ITO
YEARS = ("1900", "1980", "2020", "2060")
GMOL = 1e-9


def _arr(x):
    return np.array([np.nan if v is None else v for v in np.ravel(x)], float).reshape(np.shape(x))


def _levels(*fields, n: int = 15):
    """Contour levels spanning every field (NaN ignored), so no section falls outside them."""
    a = np.concatenate([_arr(f).ravel() for f in fields])
    return np.linspace(np.nanmin(a), np.nanmax(a), n)


def fig0_schematic(d, stem):
    """Geometry and forcing of the channel (tutorial_reentrant_channel + pkg/dic, spec Sec. 3)."""
    g = d["figdata"]["geometry"]
    y, H = _arr(g["y_km"]), _arr(g["depth_m"])
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(7.0, 4.6), sharex=True, gridspec_kw={"height_ratios": [1, 2.2]})
    a0.plot(y, g["tau_x"], color=C[0])
    a0.set_ylabel(r"$\tau_x$ (N m$^{-2}$)", color=C[0])
    b0 = a0.twinx()
    sst = _arr(g["sst_mean"])
    b0.fill_between(y, sst - g["sst_amp"], sst + g["sst_amp"], color=C[1], alpha=0.2, lw=0)
    b0.plot(y, sst, color=C[1])
    b0.set_ylabel(r"SST$_{relax}$ (°C)", color=C[1])
    a0.set_title("(a) wind stress (fixed) and SST restoring with its seasonal range", loc="left")
    a1.fill_between(y, -H, -H.max() * 1.02, color="0.82", lw=0)
    a1.plot(y, -H, color="0.2", lw=1)
    a1.axvspan(y[0] - (y[1] - y[0]) / 2, g["R_edge_km"], color=C[2], alpha=0.15, lw=0)
    a1.text(g["R_edge_km"] / 2, -0.12 * H.max(), "region R\n(upwelling)", ha="center", va="center", fontsize=8)
    if g["sponge_y_km"]:
        a1.axvspan(*g["sponge_y_km"], color="0.6", alpha=0.3, lw=0)
        a1.text(np.mean(g["sponge_y_km"]), -0.5 * H.max(), "sponge (pre-industrial)", ha="center", va="center",
                fontsize=8, rotation=90)
    a1.set_ylim(-H.max() * 1.02, 0)
    a1.set_xlim(0, y[-1] + (y[1] - y[0]) / 2)
    a1.set_xlabel("y (km)")
    a1.set_ylabel("depth (m)")
    a1.set_title("(b) cross-channel section (periodic in x); air-sea CO$_2$ exchange everywhere", loc="left")
    return style.save(fig, stem)


def fig1_sections(d, stem):
    """x-mean DIC of GT in the decades centred on 1900 and 2060, and G200 by row."""
    s = d["figdata"]["sections"]
    y, z = _arr(s["y_km"]), _arr(s["z_m"])
    fig = plt.figure(figsize=(7.0, 6.0), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.6, 1])
    top = []
    lv = _levels(s["dic"]["1900"], s["dic"]["2060"])
    for k, yr in enumerate(("1900", "2060")):
        a = fig.add_subplot(gs[0, k], sharey=top[0] if top else None)
        top.append(a)
        cf = a.contourf(y, -z, _arr(s["dic"][yr]), levels=lv, cmap="viridis")
        a.axvline(d["figdata"]["geometry"]["R_edge_km"], color="w", lw=0.8, ls="--")
        a.set_title(f"({'ab'[k]}) DIC, GT {yr} decade", loc="left")
        a.set_xlabel("y (km)")
        if k == 0:
            a.set_ylabel("depth (m)")
        else:
            a.tick_params(labelleft=False)
    fig.colorbar(cf, ax=top, shrink=0.9, label=r"DIC ($\mu$mol kg$^{-1}$)", format="%.0f")
    a = fig.add_subplot(gs[1, :])
    for k, yr in enumerate(("1900", "2060")):
        a.plot(y, _arr(s["g200_y"][yr]), color=C[k], label=f"{yr}: $G_{{200}}$(R) = {s['g200_R'][yr]:.3f}")
    a.axvspan(0, d["figdata"]["geometry"]["R_edge_km"], color=C[2], alpha=0.12, lw=0)
    a.set_xlabel("y (km)")
    a.set_ylabel(r"$G_{200}$ ($\mu$mol kg$^{-1}$ m$^{-1}$)")
    a.set_title("(c) upper-200 m DIC gradient by row; shading = R", loc="left")
    a.legend(frameon=False, fontsize=8)
    return style.save(fig, stem)


def fig2_series(d, stem):
    """Annual channel-mean flux of GC, GT and the wind branches."""
    ser = d["figdata"]["series"]
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for k, p in enumerate(("GC", "GT")):
        if p in ser:
            ax.plot(ser[p]["years"], ser[p]["flux"], color=C[k], lw=1.2, label=p)
        for yr in YEARS:
            b = ser.get(f"{p}_w{yr}")
            if b:
                ax.plot(b["years"], b["flux"], color=C[k], lw=0.9, ls="--")
    ax.plot([], [], color="0.3", ls="--", label=r"wind $\times$1.2 branches")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel("year")
    ax.set_ylabel(r"$F$ (mol C m$^{-2}$ yr$^{-1}$, + into ocean)")
    ax.legend(frameon=False, fontsize=8)
    return style.save(fig, stem)


def fig3_central(d, stem):
    """S_GT, S_nat and S_ant of each branch against the parent G200 of GT (H6a)."""
    h = d["hypotheses"]
    g = np.array([h["terms"][y]["parent"]["g200"] for y in YEARS])
    sg = np.array([h["s_gt"][y] for y in YEARS])
    sn = np.array([h["h6a"]["s_nat"][y] for y in YEARS])
    sa = np.array([h["h6a"]["s_ant"][y] for y in YEARS])
    fig, ax = plt.subplots(figsize=(6.0, 4.3))
    ax.plot(g, sg, "o-", color=C[1], label=r"$S_{GT}$")
    ax.plot(g, sn, "s--", color=C[0], mfc="none", label=r"$S_{nat} = S_{GC}$")
    ax.plot(g, sa, "^:", color=C[2], label=r"$S_{ant} = S_{GT} - S_{GC}$")
    for x, s_, yr in zip(g, sg, YEARS):
        ax.annotate(yr, (x, s_), textcoords="offset points", xytext=(4, 4), fontsize=8)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.invert_xaxis()                              # G200 erodes with time: later branches to the right
    ax.set_xlabel(r"$G_{200}$ of GT in R ($\mu$mol kg$^{-1}$ m$^{-1}$; reversed axis)")
    ax.set_ylabel(r"$S$ (mol C m$^{-2}$ yr$^{-1}$)")
    ax.set_title(f"H6a {'passes' if h['h6a']['pass'] else 'fails'}: "
                 rf"$\Delta S$ = {h['h6a']['dS']:+.3f}", loc="left")
    ax.legend(frameon=False, fontsize=8)
    return style.save(fig, stem)


def fig4_h6b(d, stem):
    """Offline terms (Jarnikova et al. 2025 eqs. 1, 4, PO4 bio) against the measured budget change (H6b)."""
    h = d["hypotheses"]
    keys = (("mld", "eq. 1 (MLD)"), ("moc", "eq. 4 (MOC)"), ("bio", "bio (PO4)"),
            ("sum_circ_bio", "offline sum"))
    x = np.arange(len(YEARS))
    w = 0.14
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for k, (key, lab) in enumerate(keys):
        ax.bar(x + (k - 2.5) * w, [h["terms"][y][key] * GMOL for y in YEARS], w, color=C[k], label=lab)
    ax.bar(x + 1.5 * w, [h["meas"][y] * GMOL for y in YEARS], w, color="0.15", label="measured (budget)")
    ax.bar(x + 2.5 * w, [h["terms"][y].get("temp", np.nan) * GMOL for y in YEARS], w, color="w",
           edgecolor="0.3", hatch="///", label="temperature (not in sum)")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xticks(x, [f"GT {y}" for y in YEARS])
    ax.set_ylabel(r"branch $-$ parent (Gmol C yr$^{-1}$)")
    ax.set_title(f"H6b {'passes' if h['h6b']['pass'] else 'fails'} "
                 f"({h['h6b']['n_within']}/4 branches within ±25 %)", loc="left")
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    return style.save(fig, stem)


def fig5_fine(d, stem):
    """Fine (10 km) members against the coarse model at 1900 and 2060 (H6c)."""
    x, h = d["crosscheck"], d["hypotheses"]
    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    for k, yr in enumerate(("1900", "2060")):
        m = _arr(x["S_fine_members"][yr])
        ax.plot(np.full(m.size, k) + np.linspace(-0.05, 0.05, m.size), m, "o", color=C[1], mfc="none")
        ax.plot(k, m.mean(), "o", color=C[1], ms=8, label="10 km members / mean" if k == 0 else None)
        ax.plot(k + 0.2, h["s_gt"][yr], "s", color=C[0], ms=7, label="50 km (GT)" if k == 0 else None)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xticks([0.1, 1.1], ["1900", "2060"])
    ax.set_xlim(-0.5, 1.6)
    ax.set_ylabel(r"$S$ (mol C m$^{-2}$ yr$^{-1}$)")
    c = x["H6c"]
    ax.set_title(f"H6c: ratio fine/coarse = {c['ratio']:.2f}" + (" (indicative)" if c["indicative"] else ""),
                 loc="left")
    ax.legend(frameon=False, fontsize=8)
    return style.save(fig, stem)


FIGS = {"fig0_schematic": fig0_schematic, "fig1_sections": fig1_sections, "fig2_series": fig2_series,
        "fig3_central": fig3_central, "fig4_h6b": fig4_h6b, "fig5_fine": fig5_fine}
NEEDS = {"fig0_schematic": ("figdata",), "fig1_sections": ("figdata",), "fig2_series": ("figdata",),
         "fig3_central": ("hypotheses",), "fig4_h6b": ("hypotheses",), "fig5_fine": ("crosscheck", "hypotheses")}


def main() -> list[str]:
    data = {n: json.loads((FIG / f"{n}.json").read_text())
            for n in ("figdata", "hypotheses", "crosscheck") if (FIG / f"{n}.json").exists()}
    made = []
    for name, fn in FIGS.items():
        if all(n in data for n in NEEDS[name]):
            print(fn(data, FIG / name))
            made.append(name)
        else:
            print(f"{name}: skipped (missing {', '.join(n for n in NEEDS[name] if n not in data)}.json)")
        plt.close("all")
    return made


if __name__ == "__main__":
    main()
