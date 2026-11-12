#!/usr/bin/env python3
"""
figura_energia_balance_kf.py
============================
Lámina de la sección "Energía de las ondas internas": balance de energía de la
bahía con las ecuaciones de Kang y Fringer (2012), integradas de 0 a t:

    barotrópico:  ΔE_bt + F_bt = P_bt − C − D_bt
    baroclínico:  ΔE_bc + F_bc = P_bc + C − D_bc

P trabajo del viento, F energía que sale por la boca (trabajo de presión más
advección, Kang y Fringer ecs. 15–16), C conversión
barotrópica → baroclínica, D residuo (disipación y términos no lineales).
Datos de balance_kang_fringer.py (npz_energy/{esc}_balance_kf_bahia.npz).

Columnas: escenarios. Filas: barotrópico, baroclínico y, a la escala de las
ondas, la boca: el trabajo de presión F_bc (el flujo de las ondas) que sale de la
bahía frente a la conversión en el escalón
600 → 163 m de la boca (fila 323, fuera de la caja; C_escalon de
balance_kang_fringer.py). Las ondas generadas en el escalón radian al océano
sin cruzar la cara de la boca, así que la conversión despreciable DENTRO de la
bahía no implica que la bahía no genere ondas por conversión.

Uso:  /opt/anaconda3/envs/pyocean/bin/python figura_energia_balance_kf.py
Salida: figs/tesis_energia_balance_kf.{pdf,png}
"""
from __future__ import annotations

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import estructura_vertical as ev

NPZ = ev.BASE / "npz_energy"
ESCENARIOS = ["2layer", "linear", "realistic"]
LINEAS = [("P", "P viento", "k", "-"), ("dE", "ΔE", "0.55", "-"),
          ("F", "F boca (presión + advección)", "#2166ac", "-"), ("C", "C", "#b2182b", "-"),
          ("D", "D", "0.3", (0, (3, 1.5)))]


def _acumular(t_h, x):
    return np.concatenate(([0.0], np.cumsum(0.5 * (x[1:] + x[:-1]) * np.diff(t_h) * 3600.0)))


def terminos(r, ape_campo=None):
    """Términos acumulados [J] de los dos balances (dict bt / bc). La APE es la
    de la interfaz (½ρ₀g′ζ_i²) si el npz la trae con valores (dos capas) y la
    exacta en los demás; ape_campo='ape' fuerza la exacta."""
    t = r["t_horas"]
    if ape_campo is None:
        ai = r["ape_interfaz"] if "ape_interfaz" in r else None
        ape_campo = "ape_interfaz" if ai is not None and np.isfinite(ai).all() else "ape"
    E = {"bt": r["ke_bt"] + r["pe_bt"], "bc": r["ke_bc"] + r[ape_campo]}
    C = _acumular(t, r["C"])
    out = {}
    for comp, signo in (("bt", -1.0), ("bc", 1.0)):
        Fp = _acumular(t, r[f"F_{comp}"])
        a = r.get(f"A_{comp}", np.zeros_like(t))
        if comp == "bc" and ape_campo == "ape_interfaz" and "A_ape_interfaz" in r:
            a = a - r["A_ape"] + r["A_ape_interfaz"]        # se advecta la APE de la interfaz
        A = _acumular(t, a)
        x = dict(P=_acumular(t, r[f"P_{comp}"]), Fp=Fp, A=A, F=Fp + A,
                 C=C, dE=E[comp] - E[comp][0])
        x["D"] = x["P"] + signo * C - x["dE"] - x["F"]
        out[comp] = x
    if "C_escalon" in r:                    # fuera de la caja: no entra al balance
        out["escalon"] = _acumular(t, r["C_escalon"])
    return out


def _panel(ax, t, x, comp):
    for k, _, color, ls in LINEAS:
        y = x[k] if not (k == "C" and comp == "bt") else -x[k]
        ax.plot(t, y / 1e9, color=color, ls=ls, lw=1.0)
    ax.axhline(0, color="0.6", lw=0.4)
    ax.text(0.97, 0.42 if comp == "bt" else 0.45,
            f"F = {x['F'][-1] / 1e9:.2f} GJ\nC = {x['C'][-1] / 1e9:.3f} GJ",
            transform=ax.transAxes,
            ha="right", va="center", fontsize=6.5,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.6))


def _boca(ax, t, T):
    """F_bc que sale de la bahía frente a la conversión en el escalón de la boca."""
    ax.plot(t, T["bc"]["Fp"] / 1e9, color="#2166ac", lw=1.0)
    ax.plot(t, T["escalon"] / 1e9, color="#b2182b", lw=1.0, ls=(0, (3, 1.5)))
    ax.axhline(0, color="0.6", lw=0.4)
    ax.text(0.03, 0.80, f"F$_{{presión}}$ = {T['bc']['Fp'][-1] / 1e9:.2f} GJ\n"
            f"C$_{{esc}}$ = {T['escalon'][-1] / 1e9:.2f} GJ", transform=ax.transAxes,
            ha="left", va="top", fontsize=6.5,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.6))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--formatos", nargs="+", default=["pdf", "png"])
    args = ap.parse_args()
    ev.estilo()
    fig, axs = plt.subplots(3, 3, figsize=(190 * ev.MM, 150 * ev.MM), sharex=True,
                            gridspec_kw=dict(left=0.10, right=0.985, top=0.91,
                                             bottom=0.08, wspace=0.10, hspace=0.12))
    letras = iter("adgbehcfi")
    for j, esc in enumerate(ESCENARIOS):
        r = np.load(NPZ / f"{esc}_balance_kf_bahia.npz")
        T = terminos(r)
        for i, comp in enumerate(("bt", "bc")):
            ax = axs[i, j]
            _panel(ax, r["t_horas"], T[comp], comp)
            ax.text(0.03, 0.96, f"({next(letras)})", transform=ax.transAxes, va="top",
                    fontweight="bold", bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.6))
            print(f"{esc:9s} {comp}: " + "  ".join(
                f"{k}={T[comp][k][-1] / 1e9:8.3f}" for k, *_ in LINEAS) + " GJ")
        _boca(axs[2, j], r["t_horas"], T)
        axs[2, j].text(0.03, 0.96, f"({next(letras)})", transform=axs[2, j].transAxes,
                       va="top", fontweight="bold",
                       bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.6))
        print(f"{esc:9s} escalón: C_esc = {T['escalon'][-1] / 1e9:.3f} GJ")
        axs[0, j].set_title(ev.SCEN_LABEL[esc], fontsize=9)
    for i, (nombre, lim) in enumerate((("barotrópico", (-0.7, 2.1)),
                                       ("baroclínico", (-20, 430)),
                                       ("boca", (-0.45, 1.25)))):
        for ax in axs[i]:
            ax.set_ylim(*lim)
        for ax in axs[i, 1:]:
            ax.tick_params(labelleft=False)
        axs[i, 0].set_ylabel("energía acumulada (GJ)", labelpad=1)
        axs[i, 0].text(-0.30, 0.5, nombre, transform=axs[i, 0].transAxes, rotation=90,
                       ha="center", va="center", fontsize=9)
    for ax in axs[2]:
        ax.set_xlim(0, 120)
        ax.set_xticks(range(0, 121, 20))
        ax.set_xlabel("tiempo (h)", labelpad=1)
    h = [plt.Line2D([], [], color=c, ls=ls, lw=1.2, label=n) for _, n, c, ls in LINEAS]
    h.append(plt.Line2D([], [], color="#b2182b", ls=(0, (3, 1.5)), lw=1.2,
                        label="C escalón de la boca (fuera de la bahía)"))
    fig.legend(handles=h, loc="upper center", ncol=6, frameon=False,
               bbox_to_anchor=(0.54, 1.0), fontsize=6.3)
    for f in args.formatos:
        p = ev.BASE / "figs" / f"tesis_energia_balance_kf.{f}"
        fig.savefig(p, dpi=args.dpi, facecolor="white")
        print(f"  → {p}")


if __name__ == "__main__":
    main()
