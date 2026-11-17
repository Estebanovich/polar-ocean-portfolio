"""Figures 0-5 of the Exp. 5 README (plan Task 7), from the committed figures/*.json only:
mean_state.json and diagnostics.json (analysis/figdata.py), sensitivity.json and validation.json
(analysis/pd2_eval.py), crosscheck.json (analysis/crosscheck.py). Flat article style (pfcommon.style).

Usage: python -m exp5_adjoint_slope.analysis.figures"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pfcommon import style

from . import h5a

style.apply()
import matplotlib.pyplot as plt  # noqa: E402

FIG = Path(__file__).resolve().parents[1] / "figures"
C = style.OKABE_ITO
KM = 1e-3
GW = 1e3          # TW -> GW


def _arr(x):
    return np.array([np.nan if v is None else v for v in np.ravel(x)], float).reshape(np.shape(x))


def fig0_schematic(d, stem):
    """Geometry and forcing of the channel (Dettling et al. 2024, Table 1)."""
    ms = d["mean_state"]
    y, H, tau = np.array(ms["y_m"]) * KM, np.array(ms["depth_m"]), np.array(ms["tau_x"])
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(7.0, 4.3), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 2.4]})
    a0.plot(y, tau, color=C[0])
    a0.axhline(0, color="0.6", lw=0.6)
    a0.set_ylabel(r"$\tau_x$ (N m$^{-2}$)")
    a0.set_title("(a) easterly wind stress", loc="left")
    a1.fill_between(y, -H, -3100, color="0.82", lw=0)
    a1.plot(y, -H, color="0.2", lw=1)
    a1.axvspan(0, 50, color=C[4], alpha=0.25, lw=0)
    a1.text(25, -250, "polynya\nsalt flux", ha="center", va="center", fontsize=8)
    a1.axvspan(400, 450, color="0.6", alpha=0.25, lw=0)
    a1.text(425, -1500, "sponge", ha="center", va="center", fontsize=8, rotation=90)
    yq = ms["y_qb_m"] * KM
    a1.axvline(yq, color=C[1], lw=1.2)
    a1.text(yq + 4, -2900, "J line", color=C[1], fontsize=8)
    a1.axvline(125, color="0.4", lw=0.8, ls=":")
    a1.text(127, -2600, "shelf break", color="0.3", fontsize=8)
    a1.set_ylim(-3100, 0)
    a1.set_xlim(0, 450)
    a1.set_xlabel("y (km)")
    a1.set_ylabel("depth (m)")
    a1.set_title("(b) bathymetry, uniform in x (periodic, 400 km)", loc="left")
    return style.save(fig, stem)


def fig1_mean_state(d, stem):
    ms = d["mean_state"]
    y, z = np.array(ms["y_m"]) * KM, np.array(ms["z_m"])
    th, sg = _arr(ms["theta"]), _arr(ms["sigma_anom"])
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw={"width_ratios": [2.6, 1]})
    pc = a0.pcolormesh(y, z, th, cmap=style.CMAPS["temperature"], shading="nearest")
    a0.contour(y, z, np.where(np.isnan(sg), np.nan, sg), 10, colors="k", linewidths=0.4)
    a0.axvline(ms["y_qb_m"] * KM, color=C[1], lw=1)
    a0.set_xlabel("y (km)")
    a0.set_ylabel("z (m)")
    a0.set_title(r"(a) G0 year 30: $\theta$ (colour) and $\sigma$ (lines)", loc="left")
    fig.colorbar(pc, ax=a0, label=r"$\theta$ (°C)")
    J = ms["J"]
    parts = [("resolved", J["J_resolved_TW"]), ("GM bolus", J["J_bolus_TW"]),
             ("Redi", J["J_redi_TW"]), ("total", J["J_TW"])]
    a1.bar(range(4), [v * GW for _, v in parts], color=[C[0], C[2], C[3], "0.3"])
    a1.axhline(0, color="0.5", lw=0.6)
    a1.set_xticks(range(4), [n for n, _ in parts], rotation=30)
    a1.set_ylabel("onshore heat transport (GW)")
    a1.set_title("(b) J and its parts", loc="left")
    fig.tight_layout()
    return style.save(fig, stem)


def fig2_sensitivity(d, stem):
    s = d["sensitivity"]
    y = np.array(s["y_m"]) * KM
    m, se = np.array(s["dJdtaux_rows_mean"]), np.array(s["dJdtaux_rows_se"])
    fig, ax = plt.subplots(figsize=(7.0, 3.4))
    for k, lead in enumerate(s["leads_days"]):
        ax.plot(y, m[k] * GW, color=C[k], label=f"lead {lead:.0f} d")
        ax.fill_between(y, (m[k] - 2 * se[k]) * GW, (m[k] + 2 * se[k]) * GW, color=C[k], alpha=0.2, lw=0)
    ax.axhline(0, color="0.5", lw=0.6)
    ax.axvline(125, color="0.4", lw=0.8, ls=":")
    ax.axvline(170, color=C[1], lw=0.8, ls="--")
    ax.set_xlim(0, 450)
    ax.set_xlabel("y (km)")
    ax.set_ylabel(r"$\partial J/\partial\tau_x$ per row (GW per N m$^{-2}$)")
    ax.set_title("Ensemble-adjoint sensitivity of J to monthly zonal wind stress "
                 "(mean of 10 windows, ±2 s.e.)", loc="left")
    ax.legend(frameon=False)
    return style.save(fig, stem)


def fig3_validation(d, stem):
    v, g = d["validation"], d["diagnostics"]
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.6))
    a = ax[0, 0]
    pts = v["A"]["V1d"]["points"]
    adj = np.abs([p["grad_adj"] for p in pts])
    a.loglog(adj, np.abs([p["fd_grad"] for p in pts]), "o", color=C[0], label="finite difference")
    tl = [abs(p["grad_tlm"]) if p.get("grad_tlm") is not None else np.nan for p in pts]
    a.loglog(adj, tl, "x", color=C[1], label="tangent linear")
    lim = [adj.min() / 3, adj.max() * 3]
    a.plot(lim, lim, color="0.5", lw=0.6)
    a.set_xlabel("|adjoint gradient|")
    a.set_ylabel("|reference|")
    a.set_title("(a) V1: gradient check, 10 points", loc="left")
    a.legend(frameon=False, fontsize=7)
    a = ax[0, 1]
    w = sorted(g["adjoint_growth"], key=float)
    a.semilogy([float(k) for k in w], [g["adjoint_growth"][k]["max_abs_dJdtheta0"] for k in w], "o-",
               color=C[3])
    a.set_xlabel("adjoint window (days)")
    a.set_ylabel(r"max $|\partial J/\partial\theta_0|$ (TW °C$^{-1}$)")
    a.set_title("(b) backward growth of the adjoint", loc="left")
    a = ax[1, 0]
    for k, (n, js) in enumerate(g["ensemble_J_TW"].items()):
        a.plot(np.arange(1, len(js) + 1), np.array(js) * GW, color=C[k % len(C)], lw=1, label=n)
    a.set_xlabel("month")
    a.set_ylabel("monthly J (GW)")
    a.set_title(r"(c) forward ensemble, $\theta_0$ noise $\leq 10^{-6}$ °C", loc="left")
    a.legend(frameon=False, fontsize=7, ncol=2)
    a = ax[1, 1]
    pr, fw = np.array(v["B"]["V3e"]["pred_members"]) * GW, np.array(v["B"]["V3e"]["fwd_members"]) * GW
    a.plot(pr, fw, "o", color=C[2])
    lim = [0, max(pr.max(), fw.max()) * 1.1]
    a.plot(lim, lim, color="0.5", lw=0.6)
    a.set_xlabel(r"adjoint $\Delta J$ (GW)")
    a.set_ylabel(r"forward $\Delta J$ (GW)")
    a.set_title("(d) V3: adjoint against forward, 10 windows", loc="left")
    fig.tight_layout()
    return style.save(fig, stem)


def fig4_hypotheses(d, stem):
    v = d["validation"]["B"]
    y = np.array(d["sensitivity"]["y_m"])
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.8))
    a = ax[0]
    a.plot(y * KM, h5a.uniform(y), color="0.3", label="(i) uniform")
    for k, w in enumerate(h5a.WIDTHS):
        a.plot(y * KM, h5a.dipole(y, w), color=C[k], label=f"(ii) {w * KM:g} km")
    a.set_xlim(0, 450)
    a.set_xlabel("y (km)")
    a.set_ylabel(r"$\delta\tau_x$ (N m$^{-2}$)")
    a.set_title("(a) H5a patterns", loc="left")
    a.legend(frameon=False, fontsize=6)
    a = ax[1]
    ws = v["H5a"]["widths"]
    keys = list(ws)
    x = np.arange(len(keys))
    a.bar(x - 0.2, [abs(ws[k]["dj_uniform"]) * GW for k in keys], 0.4, color="0.6", label="(i)")
    a.bar(x + 0.2, [abs(ws[k]["dj_dipole"]) * GW for k in keys], 0.4, color=C[0], label="(ii)")
    dr = v["H5a"].get("direct_37.5km_years26-30")
    if dr:
        a.errorbar([-0.2, 0.2], [abs(dr["i"]["fwd"][0]) * GW, abs(dr["ii"]["fwd"][0]) * GW],
                   yerr=[2 * dr["i"]["fwd"][1] * GW, 2 * dr["ii"]["fwd"][1] * GW], fmt="k_", ms=10,
                   label="forward")
    a.set_xticks(x, keys)
    a.set_xlabel("dipole half-width (km)")
    a.set_ylabel(r"$|\Delta J|$ (GW)")
    a.set_title("(b) H5a", loc="left")
    a.legend(frameon=False, fontsize=6)
    a = ax[2]
    c = v["H5c"]
    a.errorbar(c["leads_days"], np.array(c["sens_mean"]) * GW, yerr=2 * np.array(c["sens_se"]) * GW,
               fmt="o-", color=C[1])
    a.set_xlabel("lead (days)")
    a.set_ylabel(r"$|\sum_y \partial J/\partial\tau_x|$ (GW per N m$^{-2}$)")
    a.set_title("(c) H5c", loc="left")
    fig.tight_layout()
    return style.save(fig, stem)


def fig5_crosscheck(d, stem):
    x = d["crosscheck"]
    h = x["H5b"]
    jc = np.array(list(x["J_control_TW"].values())) * GW
    jp = np.array(list(x["J_pert_TW"].values())) * GW
    fig, ax = plt.subplots(1, 2, figsize=(5.6, 2.8))
    a = ax[0]
    for k in range(min(jc.size, jp.size)):
        a.plot([0, 1], [jc[k], jp[k]], color="0.7", lw=0.8, zorder=1)
    a.scatter(np.zeros(jc.size), jc, color="0.3", zorder=2, label="control")
    a.scatter(np.ones(jp.size), jp, color=C[0], zorder=2, label="perturbed")
    a.set_xticks([0, 1], ["control", "perturbed"])
    a.set_xlim(-0.5, 1.5)
    a.set_ylabel("J, last 30 days (GW)")
    a.set_title("(a) fine-model members", loc="left")
    a = ax[1]
    p = h["dj_pred"] * GW
    a.axhspan(p / 2, 2 * p, color=C[1], alpha=0.15, lw=0, label="factor 2 of prediction")
    a.errorbar([0], [p], yerr=[2 * x["dj_pred_se_TW"] * GW], fmt="s", color=C[1], label="coarse prediction")
    a.errorbar([1], [h["dj"] * GW], yerr=[2 * h["sigma_delta"] * GW], fmt="o", color=C[0], label="fine model")
    a.axhline(0, color="0.5", lw=0.6)
    a.set_xticks([0, 1], ["coarse", "fine"])
    a.set_xlim(-0.6, 1.6)
    a.set_ylabel(r"$\Delta J$ (GW), $\pm 2\sigma$")
    a.set_title("(b) H5b", loc="left")
    a.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    return style.save(fig, stem)


FIGS = {"fig0_schematic": fig0_schematic, "fig1_mean_state": fig1_mean_state,
        "fig2_sensitivity": fig2_sensitivity, "fig3_validation": fig3_validation,
        "fig4_hypotheses": fig4_hypotheses, "fig5_crosscheck": fig5_crosscheck}


def main() -> None:
    data = {n: json.loads((FIG / f"{n}.json").read_text())
            for n in ("mean_state", "sensitivity", "validation", "diagnostics", "crosscheck")}
    for name, fn in FIGS.items():
        print(fn(data, FIG / name))


if __name__ == "__main__":
    main()
