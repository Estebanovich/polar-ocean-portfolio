#!/usr/bin/env python3
"""
prototipo_hovmoller_cpcg.py
===========================
PROTOTIPO — Hovmöller DISTANCIA–TIEMPO para visualizar la propagación de ondas
internas y estimar velocidad de fase (c_p) y de grupo (c_g), con las rectas
modales teóricas c_n superpuestas.

Fuente: secciones yz `sec{W}Anom_yz_Gamma.npz` (anom = t × z × coord), donde
`coord` es el eje N–S de propagación desde la bahía. A una profundidad fija se
extrae A(t, coord) → Hovmöller.

Física del caso IWBay: las ondas internas están enganchadas al SEICHE barotrópico
(T=14430 s ≈ 0.25 cph), no al viento de 1200 s. Se demodula/filtra en la banda
del seiche para aislar la señal. c_p = pendiente de las franjas; c_g = pendiente
del frente de energía (mejor visible en la etapa LIBRE, t>60 h).

Paneles:
  (1) Hovmöller crudo + rectas modales c_n (desde el origen de generación)
  (2) Hovmöller filtrado en banda del seiche (franjas limpias → c_p) + envolvente
  (3) espectro f–k con línea f_seiche y rayos c_n (c_p = ω/k)

Uso (entorno base con matplotlib/cmocean):
    /opt/anaconda3/bin/python prototipo_hovmoller_cpcg.py --scen realistic
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_BASE = Path(__file__).resolve().parent
NPZ_DIR = _BASE / "npz_out"
MODES_DIR = _BASE / "modes_npz"
FIG_DIR = _BASE / "figs"

SEC_TAG = {"Temp": "secThetaAnom", "W": "secWAnom", "U": "secUAnom", "V": "secVAnom"}
SCEN_TITLE = {"2layer": "Dos capas", "linear": "Lineal", "realistic": "Realista"}

T_SEICHE_S = 14430.0                 # s (0.25 cph) — periodo del seiche barotrópico
F_SEICHE_CPH = 3600.0 / T_SEICHE_S   # ≈ 0.2495 cph
SEICHE_BAND_CPH = (0.235, 0.290)     # banda del seiche (mapas_psd_dominio.py)
T_FREE_H = 60.0                      # inicio de la etapa libre (viento OFF)


def load_section(scen: str, var: str) -> dict:
    with np.load(NPZ_DIR / f"{scen}_{SEC_TAG[var]}_yz_Gamma.npz", allow_pickle=True) as d:
        return dict(
            t_h=np.asarray(d["t_horas"], float),
            RC=np.asarray(d["RC"], float),
            coord=np.asarray(d["coord"], float),
            anom=np.asarray(d["anom"], np.float32),
            mask=np.asarray(d["mask2d"]).astype(bool),      # True = agua
            labels=[str(x) for x in d["station_labels"]],
            scoord=np.asarray(d["station_coords"], float),
        )


def load_cn(scen: str) -> np.ndarray:
    """c_n del océano profundo (donde radia la onda). Devuelve al menos [c1].
    2layer sólo tiene un modo baroclínico → array de longitud 1."""
    for name in (f"{scen}_modes_deep_t0.npz", f"{scen}_mode_deep_t0.npz",
                 f"{scen}_mode_2layer_CR2005.npz",
                 f"{scen}_modes_verticales_tRef.npz"):
        p = MODES_DIR / name
        if p.exists():
            with np.load(p, allow_pickle=True) as d:
                return np.atleast_1d(np.asarray(d["c"], float))
    raise FileNotFoundError(f"No hay modos para {scen}")


def bandpass_time(A: np.ndarray, dt_s: float, band_cph) -> np.ndarray:
    """Filtro pasa-banda en tiempo (eje 0) por FFT: conserva band_cph."""
    nt = A.shape[0]
    f_cph = np.fft.rfftfreq(nt, d=dt_s) * 3600.0
    F = np.fft.rfft(A - A.mean(axis=0, keepdims=True), axis=0)
    keep = (f_cph >= band_cph[0]) & (f_cph <= band_cph[1])
    F[~keep] = 0.0
    return np.fft.irfft(F, n=nt, axis=0)


def envelope_time(A_bp: np.ndarray) -> np.ndarray:
    """Envolvente por transformada de Hilbert analítica (eje tiempo)."""
    nt = A_bp.shape[0]
    Fa = np.fft.fft(A_bp, axis=0)
    h = np.zeros(nt)
    h[0] = 1
    if nt % 2 == 0:
        h[nt // 2] = 1; h[1:nt // 2] = 2
    else:
        h[1:(nt + 1) // 2] = 2
    analytic = np.fft.ifft(Fa * h[:, None], axis=0)
    return np.abs(analytic)


def make_figure(scen: str, var: str, depth_m: float | None, dpi: int) -> Path:
    S = load_section(scen, var)
    cn = load_cn(scen)
    t_h = S["t_h"]
    dt_s = float((t_h[1] - t_h[0]) * 3600.0)
    ck = S["coord"] / 1e3                          # km (eje de propagación)
    anom = S["anom"]                               # (t, z, coord)

    # profundidad: la más energética en la banda media (núcleo del modo interno);
    # se excluye la capa superficial, donde ΔW es señal del seiche, no de la onda.
    water_col = S["mask"].any(axis=1)              # niveles con agua
    rms_z = np.sqrt(np.nanmean(anom ** 2, axis=(0, 2)))
    if depth_m is None:
        midband = water_col & (S["RC"] < -50) & (S["RC"] > -400)
        rms_masked = np.where(midband, rms_z, -1)
        kz = int(np.argmax(rms_masked))
    else:
        kz = int(np.argmin(np.abs(S["RC"] - (-abs(depth_m)))))
    z_sel = float(S["RC"][kz])

    # recorte a la región activa (donde hay señal) + margen
    var_x_full = np.nanvar(anom[:, kz, :], axis=0)
    active = np.where(var_x_full > 0.03 * var_x_full.max())[0]
    if active.size:
        lo = max(0, active.min() - 8); hi = min(len(ck), active.max() + 8)
        csl = slice(lo, hi)
    else:
        csl = slice(None)
    ck = ck[csl]

    A = anom[:, kz, csl].astype(float)             # (t, coord) recortado
    # normaliza por su máximo absoluto para color simétrico
    A_n = A / (np.nanmax(np.abs(A)) or 1.0)

    # banda del seiche → franjas limpias
    A_bp = bandpass_time(A, dt_s, SEICHE_BAND_CPH)
    A_bp_n = A_bp / (np.nanmax(np.abs(A_bp)) or 1.0)
    env = envelope_time(A_bp)                       # envolvente (para c_g)

    # origen de generación: coord de máxima varianza filtrada (boca de bahía)
    var_x = np.nanvar(A_bp, axis=0)
    x0_km = float(ck[np.nanargmax(var_x)])

    # ── figura: 3 paneles ────────────────────────────────────────────────────
    fig, axs = plt.subplots(1, 3, figsize=(16.5, 6.2))

    try:
        import cmocean
        cmap = cmocean.cm.balance
    except Exception:
        cmap = "RdBu_r"

    def draw_hov(ax, field, title):
        pc = ax.pcolormesh(ck, t_h, field, cmap=cmap, vmin=-1, vmax=1,
                           shading="auto", rasterized=True)
        ax.axhline(T_FREE_H, color="0.2", lw=1.0, ls="--")
        ax.text(ck.min(), T_FREE_H + 1, " viento OFF (etapa libre)", fontsize=8,
                color="0.2", va="bottom")
        ax.set_xlabel("distancia N–S  (km)")
        ax.set_ylabel("tiempo  (h)")
        ax.set_title(title, fontsize=11)
        return pc

    # panel 1: crudo + rectas modales c_n desde (x0, 0)
    pc = draw_hov(axs[0], A_n, f"(1) Δ{var} crudo · z={z_sel:.0f} m")
    tt = np.array([0.0, t_h.max()])
    for n, c in enumerate(cn, start=1):
        # ola que se aleja hacia -N–S (ajusta signo si radia al otro lado)
        for sgn in (+1, -1):
            xline = x0_km + sgn * c * (tt * 3600.0) / 1e3
            axs[0].plot(xline, tt, lw=1.3, ls="-",
                        color=plt.cm.viridis(n / (len(cn) + 1)),
                        label=f"c$_{n}$={c:.2f} m/s" if sgn == +1 else None)
    axs[0].plot([x0_km], [0], "k*", ms=9)
    axs[0].legend(fontsize=8, loc="upper right", framealpha=0.85)
    axs[0].set_xlim(ck.min(), ck.max()); axs[0].set_ylim(0, t_h.max())

    # panel 2: filtrado en banda seiche (franjas = c_p) + contorno de envolvente
    draw_hov(axs[1], A_bp_n,
             f"(2) filtrado banda seiche {SEICHE_BAND_CPH[0]:.2f}–{SEICHE_BAND_CPH[1]:.2f} cph")
    env_n = env / (np.nanmax(env) or 1.0)
    axs[1].contour(ck, t_h, env_n, levels=[0.3, 0.6], colors="k",
                   linewidths=0.7, alpha=0.7)
    axs[1].set_xlim(ck.min(), ck.max()); axs[1].set_ylim(0, t_h.max())

    # panel 3: espectro f–k
    win = np.hanning(A.shape[0])[:, None] * np.hanning(A.shape[1])[None, :]
    P = np.abs(np.fft.fftshift(np.fft.fft2(A_bp * win))) ** 2
    dk = float(np.median(np.abs(np.diff(ck))))     # km
    k_cpkm = np.fft.fftshift(np.fft.fftfreq(A.shape[1], d=dk))   # ciclos/km
    f_cph = np.fft.fftshift(np.fft.fftfreq(A.shape[0], d=dt_s)) * 3600.0
    Pn = P / P.max()
    ax = axs[2]
    pcf = ax.pcolormesh(k_cpkm, f_cph, 10 * np.log10(Pn + 1e-6),
                        cmap="magma", vmin=-40, vmax=0, shading="auto",
                        rasterized=True)
    ax.axhline(F_SEICHE_CPH, color="cyan", lw=1.2, ls="--")
    ax.text(k_cpkm.min(), F_SEICHE_CPH, " f seiche", color="cyan", fontsize=8,
            va="bottom")
    # rayos c_n: f_cph = c_n * k_cpkm  (ω=c k → f[cph]=c[m/s]*k[cyc/km]*3600/1000)
    kk = np.array([0, k_cpkm.max()])
    for n, c in enumerate(cn, start=1):
        ax.plot(kk, c * kk * 3600.0 / 1000.0, lw=1.1,
                color=plt.cm.viridis(n / (len(cn) + 1)), label=f"c$_{n}$")
    ax.set_xlim(0, k_cpkm.max()); ax.set_ylim(0, min(1.0, f_cph.max()))
    ax.set_xlabel("número de onda k  (ciclos/km)")
    ax.set_ylabel("frecuencia  (cph)")
    ax.set_title("(3) espectro f–k  (banda seiche)")
    ax.legend(fontsize=8, loc="upper right")

    fig.colorbar(pc, ax=axs[:2], shrink=0.7, pad=0.01,
                 label=f"Δ{var} / máx|Δ{var}|")
    fig.colorbar(pcf, ax=axs[2], shrink=0.7, pad=0.02, label="PSD (dB)")
    fig.suptitle(
        f"Hovmöller distancia–tiempo · {SCEN_TITLE.get(scen, scen)} · "
        f"ondas internas @ seiche (T={T_SEICHE_S/60:.0f} min) · "
        f"c_p = pendiente franjas, c_g = frente de energía",
        fontsize=12)
    fig.subplots_adjust(left=0.05, right=0.97, top=0.90, bottom=0.10, wspace=0.28)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = FIG_DIR / f"proto_hovmoller_cpcg_{scen}_{var}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[{scen}/{var}] z={z_sel:.0f} m · x0={x0_km:.0f} km · "
          f"c_n={np.round(cn,2)} → {out}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scen", default="realistic",
                    choices=["2layer", "linear", "realistic"])
    ap.add_argument("--var", default="W", choices=list(SEC_TAG))
    ap.add_argument("--depth", type=float, default=None,
                    help="profundidad fija (m, positiva); def. la más energética")
    ap.add_argument("--dpi", type=int, default=160)
    args = ap.parse_args()
    make_figure(args.scen, args.var, args.depth, args.dpi)


if __name__ == "__main__":
    main()
