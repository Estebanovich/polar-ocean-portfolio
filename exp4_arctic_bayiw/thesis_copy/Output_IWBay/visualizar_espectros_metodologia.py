"""
visualizar_espectros_metodologia.py
Series de tiempo + PSD para Δη, ΔT, Δw, Δu/Δv
Variables y estaciones definidas en metodología de tesis.

Figuras producidas
──────────────────
Por estación, comparando escenarios:
  fig1_eta_{stn}.png              Δη  serie + PSD
  fig2_temp_{stn}_{d}m.png        ΔT  serie + PSD  (110 y 250 m)
  fig3_w_{stn}_{d}m.png           Δw  baroclinica serie + PSD
  fig4_rotary_{stn}_{d}m.png      Δu+iΔv  espectros rotatorios CW/CCW

Por escenario, todas las estaciones juntas:
  fig5_eta_{scen}.png             Δη  todas estaciones
  fig6_temp_{scen}_{d}m.png       ΔT  todas estaciones
  fig7_w_{scen}_{d}m.png          Δw  todas estaciones
  fig8_rotary_{scen}_{d}m.png     rotatorio CW/CCW  todas estaciones
"""
from __future__ import annotations
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend
from pathlib import Path

BASE_DIR = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay")
sys.path.insert(0, str(BASE_DIR))

# ══ CONFIGURATION ══════════════════════════════════════════════════════════════

NPZ_OUT_DIR  = BASE_DIR / "npz_out"
NPZ_BARO_DIR = BASE_DIR / "npz_baro_full"
FIGS_DIR     = BASE_DIR / "figs" / "espectros_metodologia"
FIGS_DIR.mkdir(parents=True, exist_ok=True)

SCENARIOS = ["linear", "2layer", "realistic"]
STATIONS  = ["Seiche", "A", "B", "C", "Omega"]
DEPTHS    = [110.0, 250.0]

SEICHE_CPH      = 0.259          # pico observado en los espectros
SEICHE_TEO_CPH  = 1.0 / 3.29      # = 0.304 cph; Merian 1/4-onda: T0=4L/√(gh), L=119km, h=164m
CORIOLIS_F      = 0.0322
XLIM_CPH        = (0.0, 0.5)

# Colores vivos para que las líneas verticales se distingan del grid (negro punteado).
VLINE_STYLE = {
    "f0":         dict(color="#e91e63", ls="-",  lw=1.4),   # magenta
    "seiche_obs": dict(color="#ff8c00", ls="--", lw=1.4),   # naranja
    "seiche_teo": dict(color="#00bcd4", ls="-.", lw=1.4),   # cian
}

# Bandas de frecuencia usadas en los mapas 2D de PSD (mapas_psd_dominio.py)
BANDS = {
    "NI":     (0.035, 0.065),
    #"SubNI":  (0.065, 0.085),
    "Seiche": (0.223, 0.290),
    "modo 1": (0.15, 0.2125),
}
BAND_COLORS = {
    "NI":     "#1f77b4",   # azul
    #"SubNI":  "#2ca02c",   # verde
    "Seiche": "#d62728",   # rojo
    "modo 1": "#9c27b0",   # morado
}
BAND_ALPHA = 0.12

SCENARIO_COLORS = {"linear": "#1f77b4", "2layer": "#ff7f0e", "realistic": "#2ca02c"}
SCENARIO_LS     = {"linear": "-",        "2layer": "--",      "realistic": "-."}

STATION_COLORS  = {"Seiche": "#9d02d7", "A": "#ea5f94", "B": "#0000ff", "C": "#fa8775",
                   "Omega": "#0aa3a3"}
STATION_LS      = {"Seiche": "-",       "A": "--",      "B": "-.",      "C": ":",
                   "Omega": (0, (3, 1, 1, 1))}

# ── NPZ path helpers ───────────────────────────────────────────────────────────

def _depth_tag(d: float, prefix: str) -> str:
    return f"{prefix}{int(round(d))}m"

def _eta_series_path(scen):  return NPZ_OUT_DIR / f"{scen}_series.npz"
def _eta_espec_path(scen):   return NPZ_OUT_DIR / f"{scen}_espectros.npz"

def _temp_series_path(scen, d):
    return NPZ_OUT_DIR / f"{scen}_Temp_{_depth_tag(d,'T')}_series.npz"
def _temp_espec_path(scen, d):
    return NPZ_OUT_DIR / f"{scen}_Temp_{_depth_tag(d,'T')}_espectros.npz"

def _w_series_path(scen, d):
    return NPZ_OUT_DIR / f"{scen}_W_{_depth_tag(d,'W')}_series.npz"
def _w_espec_path(scen, d):
    return NPZ_OUT_DIR / f"{scen}_W_{_depth_tag(d,'W')}_espectros.npz"

def _uv_series_path(scen, d):
    return NPZ_OUT_DIR / f"{scen}_UV_{_depth_tag(d,'UV')}_series.npz"
def _uv_espec_path(scen, d):
    return NPZ_OUT_DIR / f"{scen}_UV_{_depth_tag(d,'UV')}_espectros.npz"

def _baro_path(scen, var, label):
    return NPZ_BARO_DIR / scen / f"{var}_{label}_baroFULL.npz"


# ── utilidades ────────────────────────────────────────────────────────────────

def _valid(y: np.ndarray, thr: float = 1e-12) -> bool:
    """True si la señal tiene variación real (no celda seca ni toda ceros)."""
    return bool(np.max(np.abs(y)) > thr)

def _vlines(ax):
    ax.axvline(CORIOLIS_F,     **VLINE_STYLE["f0"],
               label=f"f₀ = {CORIOLIS_F:.4f} cph (T = {1/CORIOLIS_F:.2f} h)")
    ax.axvline(SEICHE_CPH,     **VLINE_STYLE["seiche_obs"],
               label=f"f seiche obs = {SEICHE_CPH:.4f} cph (T = {1/SEICHE_CPH:.2f} h)")
    ax.axvline(SEICHE_TEO_CPH, **VLINE_STYLE["seiche_teo"],
               label=f"f seiche teor = {SEICHE_TEO_CPH:.4f} cph (T = {1/SEICHE_TEO_CPH:.2f} h)")

def _bands(ax):
    """Sombrea las 3 bandas de frecuencia usadas en los mapas 2D de PSD."""
    for tag, (f0, f1) in BANDS.items():
        ax.axvspan(f0, f1, color=BAND_COLORS[tag], alpha=BAND_ALPHA, lw=0,
                   label=f"banda {tag} ({f0:g}–{f1:g} cph)", zorder=0)

def _finalize_spec(ax, ylabel: str):
    from matplotlib.ticker import AutoMinorLocator
    _bands(ax)
    ax.set_xlim(*XLIM_CPH)
    ax.set_xlabel("Frecuencia (cph)")
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=7, ncols=2)
    # Grid horizontal x4: 4 minor ticks entre cada major → 4× más densidad
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.grid(True, which="major", ls=":", lw=0.7, color="0.5")
    ax.grid(True, which="minor", axis="x", ls=":", lw=0.4, color="0.7")

def _finalize_ts(ax, ylabel: str, title: str):
    ax.set_xlabel("Tiempo (h)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=7)
    ax.grid(ls=":")


# ══ FIGURA 1 — Δη series + PSD (una figura por estación, escenarios solapados) ═

def fig1_eta():
    """Fig. 1: Δη serie temporal + PSD — por estación, comparando escenarios."""
    for stn in STATIONS:
        fig, (ax_ts, ax_sp) = plt.subplots(2, 1, figsize=(10, 7), sharex=False)
        has_data = False

        for scen in SCENARIOS:
            ps = _eta_series_path(scen)
            pe = _eta_espec_path(scen)
            if not ps.exists() or not pe.exists():
                continue

            ds = np.load(ps, allow_pickle=False)
            de = np.load(pe, allow_pickle=False)

            t_key   = f"{stn}_t_horas"
            eta_key = f"{stn}_eta_diff"
            f_key   = f"{stn}_freq_cph"
            p_key   = f"{stn}_PSD"

            if t_key not in ds.files or eta_key not in ds.files:
                continue

            t_h = np.asarray(ds[t_key],   dtype=float)
            eta = np.asarray(ds[eta_key],  dtype=float)
            if not _valid(eta):
                continue

            f_cph = np.asarray(de[f_key], dtype=float)
            PSD   = np.asarray(de[p_key], dtype=float) * 3600.0  # Hz⁻¹ → cph⁻¹

            col = SCENARIO_COLORS[scen]
            ls  = SCENARIO_LS[scen]
            ax_ts.plot(t_h, eta,   color=col, ls=ls, lw=1.3, label=scen)
            ax_sp.plot(f_cph, PSD, color=col, ls=ls, lw=1.3, label=scen)
            has_data = True

        if not has_data:
            plt.close(fig)
            continue

        _finalize_ts(ax_ts, "Δη (m)", f"Δη — estación {stn}")
        _vlines(ax_sp)
        if any(line.get_ydata().max() > 0 for line in ax_sp.get_lines()):
            ax_sp.set_yscale("log")
        _finalize_spec(ax_sp, "PSD(Δη) [m² cph⁻¹]")
        fig.tight_layout()
        fig.savefig(FIGS_DIR / f"fig1_eta_{stn}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


# ══ FIGURA 2 — ΔT series + PSD ═════════════════════════════════════════════════

def fig2_temperatura():
    """Fig. 2: ΔT serie temporal + PSD — por estación y profundidad."""
    for stn in STATIONS:
        for depth in DEPTHS:
            tag  = _depth_tag(depth, "T")
            fig, (ax_ts, ax_sp) = plt.subplots(2, 1, figsize=(10, 7), sharex=False)
            has_data = False

            for scen in SCENARIOS:
                ps = _temp_series_path(scen, depth)
                pe = _temp_espec_path(scen, depth)
                if not ps.exists() or not pe.exists():
                    continue

                ds = np.load(ps, allow_pickle=False)
                de = np.load(pe, allow_pickle=False)

                t_key = f"{stn}_t_horas"
                y_key = f"{stn}_{tag}"
                f_key = f"{stn}_freq_cph"
                p_key = f"{stn}_PSD_{tag}"

                if t_key not in ds.files or y_key not in ds.files:
                    continue

                y = np.asarray(ds[y_key], dtype=float)
                if not _valid(y):
                    continue

                t_h   = np.asarray(ds[t_key], dtype=float)
                f_cph = np.asarray(de[f_key], dtype=float)
                PSD   = np.asarray(de[p_key], dtype=float)

                col = SCENARIO_COLORS[scen]
                ls  = SCENARIO_LS[scen]
                ax_ts.plot(t_h, y,     color=col, ls=ls, lw=1.3, label=scen)
                ax_sp.plot(f_cph, PSD, color=col, ls=ls, lw=1.3, label=scen)
                has_data = True

            if not has_data:
                plt.close(fig)
                continue

            _finalize_ts(ax_ts, "ΔT (°C)", f"ΔT — estación {stn} — z≈{depth:.0f} m")
            _vlines(ax_sp)
            if any(line.get_ydata().max() > 0 for line in ax_sp.get_lines()):
                ax_sp.set_yscale("log")
            _finalize_spec(ax_sp, "PSD(ΔT) [°C² cph⁻¹]")
            fig.tight_layout()
            fig.savefig(FIGS_DIR / f"fig2_temp_{stn}_{int(depth)}m.png", dpi=300, bbox_inches="tight")
            plt.close(fig)


# ══ FIGURA 3 — Δw baroclinica (baroFULL) ═══════════════════════════════════════

def fig3_w_baroclinico():
    """Fig. 3: Δw baroclinica serie + PSD — por estación y profundidad (baroFULL)."""
    from scipy.signal import get_window

    for stn in STATIONS:
        for depth in DEPTHS:
            fig, (ax_ts, ax_sp) = plt.subplots(2, 1, figsize=(10, 7), sharex=False)
            has_data = False

            for scen in SCENARIOS:
                p = _baro_path(scen, "W", stn)
                if not p.exists():
                    continue

                d   = np.load(p, allow_pickle=True)
                t_s = np.asarray(d["t"],           dtype=float)
                RC  = np.asarray(d["RC"],          dtype=float)
                bc  = np.asarray(d["baroclinica"], dtype=float)

                kz   = int(np.argmin(np.abs(RC + abs(depth))))
                y    = bc[:, kz]
                if not _valid(y):
                    continue

                zref = float(RC[kz])
                t_h  = t_s / 3600.0
                y    = detrend(y, type="linear")
                dt_s = float(np.median(np.diff(t_s)))

                N   = len(y)
                win = get_window("hann", N, fftbins=True)
                U_w = float((win ** 2).sum())
                Y   = np.fft.rfft(y * win)
                PSD = (np.abs(Y) ** 2) / (U_w * (1.0 / dt_s)) * 3600.0
                f_cph = np.fft.rfftfreq(N, d=dt_s) * 3600.0

                col = SCENARIO_COLORS[scen]
                ls  = SCENARIO_LS[scen]
                ax_ts.plot(t_h, y,     color=col, ls=ls, lw=1.3,
                           label=f"{scen} (z≈{abs(zref):.0f}m)")
                ax_sp.plot(f_cph, PSD, color=col, ls=ls, lw=1.3, label=scen)
                has_data = True

            if not has_data:
                plt.close(fig)
                continue

            _finalize_ts(ax_ts, "Δw baro (m s⁻¹)",
                         f"Δw baroclinico — estación {stn} — z≈{depth:.0f} m")
            _vlines(ax_sp)
            if any(line.get_ydata().max() > 0 for line in ax_sp.get_lines()):
                ax_sp.set_yscale("log")
            _finalize_spec(ax_sp, "PSD(Δw) [m² s⁻² cph⁻¹]")
            fig.tight_layout()
            fig.savefig(FIGS_DIR / f"fig3_w_{stn}_{int(depth)}m.png", dpi=300, bbox_inches="tight")
            plt.close(fig)


# ══ FIGURA 4 — Espectros rotatorios CW/CCW (pre-calculados) ════════════════════

def fig4_rotary_uv():
    """Fig. 4: Espectros rotatorios CW/CCW — por estación, usando espectros pre-calculados."""
    for stn in STATIONS:
        for depth in DEPTHS:
            tag_uv = _depth_tag(depth, "UV")
            fig, ax = plt.subplots(figsize=(9, 5))
            has_data = False

            for scen in SCENARIOS:
                pe = _uv_espec_path(scen, depth)
                ps = _uv_series_path(scen, depth)
                if not pe.exists():
                    continue

                de = np.load(pe, allow_pickle=False)
                f_key   = f"{stn}_freq_cph"
                cw_key  = f"{stn}_PSD_CW_{tag_uv}"
                ccw_key = f"{stn}_PSD_CCW_{tag_uv}"

                if f_key not in de.files or cw_key not in de.files:
                    continue

                # verificar señal válida desde series
                if ps.exists():
                    ds = np.load(ps, allow_pickle=False)
                    u_key = f"{stn}_U{int(round(depth))}m"
                    if u_key in ds.files and not _valid(ds[u_key]):
                        continue

                f_cph = np.asarray(de[f_key],  dtype=float)
                p_cw  = np.asarray(de[cw_key], dtype=float)
                p_ccw = np.asarray(de[ccw_key], dtype=float)

                col = SCENARIO_COLORS[scen]
                ax.plot(f_cph, p_cw,  color=col, ls="-",  lw=1.4, label=f"{scen} CW")
                ax.plot(f_cph, p_ccw, color=col, ls="--", lw=1.0, alpha=0.75,
                        label=f"{scen} CCW")
                has_data = True

            if not has_data:
                plt.close(fig)
                continue

            _vlines(ax)
            if any(line.get_ydata().max() > 0 for line in ax.get_lines()):
                ax.set_yscale("log")
            ax.set_title(f"Espectros rotatorios — estación {stn} — z≈{depth:.0f} m")
            _finalize_spec(ax, "PSD(Δu+iΔv) [m² s⁻² cph⁻¹]")
            fig.tight_layout()
            fig.savefig(FIGS_DIR / f"fig4_rotary_{stn}_{int(depth)}m.png", dpi=300,
                        bbox_inches="tight")
            plt.close(fig)


# ══ FIGURAS 5-8 — Por escenario, todas las estaciones ══════════════════════════

def _plot_all_stations_ts_psd(
    t_list: list[tuple],   # [(stn_label, t_h, y), ...]
    f_list: list[tuple],   # [(stn_label, f_cph, PSD), ...]
    ts_ylabel: str,
    sp_ylabel: str,
    title: str,
) -> plt.Figure:
    """2-panel: series temporales (arriba) + PSD (abajo), todas las estaciones."""
    fig, (ax_ts, ax_sp) = plt.subplots(2, 1, figsize=(10, 7), sharex=False)

    for stn, t_h, y in t_list:
        ax_ts.plot(t_h, y, color=STATION_COLORS.get(stn, "k"),
                   ls=STATION_LS.get(stn, "-"), lw=1.3, label=stn)

    for stn, f_cph, PSD in f_list:
        ax_sp.plot(f_cph, PSD, color=STATION_COLORS.get(stn, "k"),
                   ls=STATION_LS.get(stn, "-"), lw=1.3, label=stn)

    _finalize_ts(ax_ts, ts_ylabel, title)
    _vlines(ax_sp)
    if any(line.get_ydata().max() > 0 for line in ax_sp.get_lines()):
        ax_sp.set_yscale("log")
    _finalize_spec(ax_sp, sp_ylabel)
    fig.tight_layout()
    return fig


def fig5_todas_estaciones_eta():
    """Fig. 5: Δη todas las estaciones — una figura por escenario."""
    for scen in SCENARIOS:
        ps = _eta_series_path(scen)
        pe = _eta_espec_path(scen)
        if not ps.exists() or not pe.exists():
            continue

        ds = np.load(ps, allow_pickle=False)
        de = np.load(pe, allow_pickle=False)

        t_list, f_list = [], []
        for stn in STATIONS:
            t_key = f"{stn}_t_horas";  y_key = f"{stn}_eta_diff"
            f_key = f"{stn}_freq_cph"; p_key = f"{stn}_PSD"
            if t_key not in ds.files:
                continue
            y = np.asarray(ds[y_key], dtype=float)
            if not _valid(y):
                continue
            t_list.append((stn, np.asarray(ds[t_key], dtype=float), y))
            PSD = np.asarray(de[p_key], dtype=float) * 3600.0
            f_list.append((stn, np.asarray(de[f_key], dtype=float), PSD))

        if not t_list:
            continue
        fig = _plot_all_stations_ts_psd(
            t_list, f_list,
            "Δη (m)", "PSD(Δη) [m² cph⁻¹]",
            f"Δη — todas las estaciones — {scen}",
        )
        fig.savefig(FIGS_DIR / f"fig5_eta_{scen}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def fig6_todas_estaciones_temp():
    """Fig. 6: ΔT todas las estaciones — por escenario y profundidad."""
    for scen in SCENARIOS:
        for depth in DEPTHS:
            tag = _depth_tag(depth, "T")
            ps = _temp_series_path(scen, depth)
            pe = _temp_espec_path(scen, depth)
            if not ps.exists() or not pe.exists():
                continue

            ds = np.load(ps, allow_pickle=False)
            de = np.load(pe, allow_pickle=False)

            t_list, f_list = [], []
            for stn in STATIONS:
                t_key = f"{stn}_t_horas"; y_key = f"{stn}_{tag}"
                f_key = f"{stn}_freq_cph"; p_key = f"{stn}_PSD_{tag}"
                if t_key not in ds.files:
                    continue
                y = np.asarray(ds[y_key], dtype=float)
                if not _valid(y):
                    continue
                t_list.append((stn, np.asarray(ds[t_key], dtype=float), y))
                f_list.append((stn, np.asarray(de[f_key], dtype=float),
                               np.asarray(de[p_key], dtype=float)))

            if not t_list:
                continue
            fig = _plot_all_stations_ts_psd(
                t_list, f_list,
                "ΔT (°C)", "PSD(ΔT) [°C² cph⁻¹]",
                f"ΔT — todas las estaciones — {scen} — z≈{depth:.0f} m",
            )
            fig.savefig(FIGS_DIR / f"fig6_temp_{scen}_{int(depth)}m.png", dpi=300,
                        bbox_inches="tight")
            plt.close(fig)


def fig7_todas_estaciones_w():
    """Fig. 7: Δw baroclinica todas las estaciones — por escenario y profundidad."""
    for scen in SCENARIOS:
        for depth in DEPTHS:
            tag = _depth_tag(depth, "W")
            ps = _w_series_path(scen, depth)
            pe = _w_espec_path(scen, depth)
            if not ps.exists() or not pe.exists():
                continue

            ds = np.load(ps, allow_pickle=False)
            de = np.load(pe, allow_pickle=False)

            t_list, f_list = [], []
            for stn in STATIONS:
                t_key = f"{stn}_t_horas"; y_key = f"{stn}_{tag}"
                f_key = f"{stn}_freq_cph"; p_key = f"{stn}_PSD_{tag}"
                if t_key not in ds.files or y_key not in ds.files:
                    continue
                y = np.asarray(ds[y_key], dtype=float)
                if not _valid(y):
                    continue
                t_list.append((stn, np.asarray(ds[t_key], dtype=float), y))
                f_list.append((stn, np.asarray(de[f_key], dtype=float),
                               np.asarray(de[p_key], dtype=float)))

            if not t_list:
                continue
            fig = _plot_all_stations_ts_psd(
                t_list, f_list,
                "Δw (m s⁻¹)", "PSD(Δw) [m² s⁻² cph⁻¹]",
                f"Δw — todas las estaciones — {scen} — z≈{depth:.0f} m",
            )
            fig.savefig(FIGS_DIR / f"fig7_w_{scen}_{int(depth)}m.png", dpi=300,
                        bbox_inches="tight")
            plt.close(fig)


def fig8_todas_estaciones_rotary():
    """Fig. 8: Espectros rotatorios CW/CCW todas las estaciones — por escenario."""
    for scen in SCENARIOS:
        for depth in DEPTHS:
            tag_uv = _depth_tag(depth, "UV")
            pe = _uv_espec_path(scen, depth)
            ps = _uv_series_path(scen, depth)
            if not pe.exists():
                continue

            de = np.load(pe, allow_pickle=False)
            fig, ax = plt.subplots(figsize=(9, 5))
            has_data = False

            for stn in STATIONS:
                f_key   = f"{stn}_freq_cph"
                cw_key  = f"{stn}_PSD_CW_{tag_uv}"
                ccw_key = f"{stn}_PSD_CCW_{tag_uv}"
                if f_key not in de.files or cw_key not in de.files:
                    continue

                # filtro de señal válida
                if ps.exists():
                    ds = np.load(ps, allow_pickle=False)
                    u_key = f"{stn}_U{int(round(depth))}m"
                    if u_key in ds.files and not _valid(ds[u_key]):
                        continue

                f_cph = np.asarray(de[f_key],   dtype=float)
                p_cw  = np.asarray(de[cw_key],  dtype=float)
                p_ccw = np.asarray(de[ccw_key], dtype=float)

                col = STATION_COLORS.get(stn, "k")
                ax.plot(f_cph, p_cw,  color=col, ls="-",  lw=1.4, label=f"{stn} CW")
                ax.plot(f_cph, p_ccw, color=col, ls="--", lw=1.0, alpha=0.75,
                        label=f"{stn} CCW")
                has_data = True

            if not has_data:
                plt.close(fig)
                continue

            _vlines(ax)
            if any(line.get_ydata().max() > 0 for line in ax.get_lines()):
                ax.set_yscale("log")
            ax.set_title(f"Espectros rotatorios — todas estaciones — {scen} — z≈{depth:.0f} m")
            _finalize_spec(ax, "PSD(Δu+iΔv) [m² s⁻² cph⁻¹]")
            fig.tight_layout()
            fig.savefig(FIGS_DIR / f"fig8_rotary_{scen}_{int(depth)}m.png", dpi=300,
                        bbox_inches="tight")
            plt.close(fig)


# ══ MAIN ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")

    print("fig1: Δη serie + PSD (por estación)")
    fig1_eta()
    print("fig2: ΔT serie + PSD (por estación)")
    fig2_temperatura()
    print("fig3: Δw baroclinica serie + PSD (por estación)")
    fig3_w_baroclinico()
    print("fig4: rotatorios CW/CCW (por estación)")
    fig4_rotary_uv()
    print("fig5: Δη todas las estaciones (por escenario)")
    fig5_todas_estaciones_eta()
    print("fig6: ΔT todas las estaciones (por escenario)")
    fig6_todas_estaciones_temp()
    print("fig7: Δw todas las estaciones (por escenario)")
    fig7_todas_estaciones_w()
    print("fig8: rotatorios todas las estaciones (por escenario)")
    fig8_todas_estaciones_rotary()
    print("Figuras guardadas en:", FIGS_DIR)
