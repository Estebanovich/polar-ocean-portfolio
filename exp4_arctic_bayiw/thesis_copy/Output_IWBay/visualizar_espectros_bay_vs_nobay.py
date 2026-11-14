"""
visualizar_espectros_bay_vs_nobay.py
Series de tiempo + PSD comparando el caso CON bahía (run_expand) y SIN bahía
(run_expand_nobay, control) — graficados por separado, NO su diferencia.

A diferencia de `visualizar_espectros_metodologia.py` (que grafica Δ = bahía − sin bahía),
aquí cada caso se superpone en los mismos ejes:
    con bahía  → línea sólida   (rojo)
    sin bahía  → línea punteada (azul)

Requiere que `regenerar_espectros_psd.py` haya guardado las claves separadas
`_bay` / `_nobay` en los NPZ de npz_out/.

Figuras (una por escenario × estación):
  fig_eta_{scen}_{stn}.png            η   serie + PSD
  fig_temp_{scen}_{stn}_{d}m.png      T   serie + PSD  (110 y 250 m)
  fig_w_{scen}_{stn}_{d}m.png         w   serie + PSD  (W cruda a profundidad)
  fig_rotary_{scen}_{stn}_{d}m.png    u+iv  espectros rotatorios CW/CCW
"""
from __future__ import annotations
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay")
sys.path.insert(0, str(BASE_DIR))

# Reutiliza configuración y helpers de estilo del script de diferencias.
# (importar el módulo NO ejecuta su main(): está bajo `if __name__ == "__main__"`)
from visualizar_espectros_metodologia import (
    SCENARIOS, STATIONS, DEPTHS,
    _eta_series_path, _eta_espec_path,
    _temp_series_path, _temp_espec_path,
    _w_series_path, _w_espec_path,
    _uv_series_path, _uv_espec_path,
    _depth_tag, _valid, _vlines, _finalize_spec, _finalize_ts,
)

# ══ CONFIGURATION ══════════════════════════════════════════════════════════════

FIGS_DIR = BASE_DIR / "figs" / "espectros_bay_vs_nobay"
FIGS_DIR.mkdir(parents=True, exist_ok=True)

# Estilo por caso: color distingue bahía/sin-bahía, line-style refuerza.
CASE_STYLE = {
    "bay":   dict(color="#d62728", ls="-",  lw=1.4),   # rojo, sólido
    "nobay": dict(color="#1f77b4", ls="--", lw=1.4),   # azul, punteado
}
CASE_LABEL = {"bay": "con bahía", "nobay": "sin bahía"}


def _log_if_positive(ax):
    if any(line.get_ydata().size and line.get_ydata().max() > 0 for line in ax.get_lines()):
        ax.set_yscale("log")


def _plot_two_panel(t_curves, f_curves, ts_ylabel, sp_ylabel, title, outname):
    """2 paneles (serie arriba, PSD abajo) con curvas bahía/sin-bahía superpuestas.

    t_curves : [(case, t_h, y), ...]
    f_curves : [(case, f_cph, PSD), ...]
    """
    if not t_curves:
        return
    fig, (ax_ts, ax_sp) = plt.subplots(2, 1, figsize=(10, 7), sharex=False)
    for case, t_h, y in t_curves:
        ax_ts.plot(t_h, y, label=CASE_LABEL[case], **CASE_STYLE[case])
    for case, f_cph, PSD in f_curves:
        ax_sp.plot(f_cph, PSD, label=CASE_LABEL[case], **CASE_STYLE[case])

    _finalize_ts(ax_ts, ts_ylabel, title)
    _vlines(ax_sp)
    _log_if_positive(ax_sp)
    _finalize_spec(ax_sp, sp_ylabel)
    fig.tight_layout()
    fig.savefig(FIGS_DIR / outname, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ══ η — serie + PSD ════════════════════════════════════════════════════════════

def fig_eta():
    for scen in SCENARIOS:
        ps, pe = _eta_series_path(scen), _eta_espec_path(scen)
        if not ps.exists() or not pe.exists():
            continue
        ds, de = np.load(ps, allow_pickle=False), np.load(pe, allow_pickle=False)

        for stn in STATIONS:
            t_key, f_key = f"{stn}_t_horas", f"{stn}_freq_cph"
            if t_key not in ds.files or f_key not in de.files:
                continue
            t_h   = np.asarray(ds[t_key], dtype=float)
            f_cph = np.asarray(de[f_key], dtype=float)

            t_curves, f_curves = [], []
            for case in ("bay", "nobay"):
                y_key, p_key = f"{stn}_eta_{case}", f"{stn}_PSD_{case}"
                if y_key not in ds.files or p_key not in de.files:
                    continue
                y = np.asarray(ds[y_key], dtype=float)
                if not _valid(y):
                    continue
                t_curves.append((case, t_h, y))
                PSD = np.asarray(de[p_key], dtype=float) * 3600.0  # Hz⁻¹ → cph⁻¹
                f_curves.append((case, f_cph, PSD))

            _plot_two_panel(
                t_curves, f_curves,
                "η (m)", "PSD(η) [m² cph⁻¹]",
                f"η — bahía vs sin bahía — {scen} — estación {stn}",
                f"fig_eta_{scen}_{stn}.png",
            )


# ══ T — serie + PSD ════════════════════════════════════════════════════════════

def fig_temp():
    for scen in SCENARIOS:
        for depth in DEPTHS:
            tag = _depth_tag(depth, "T")
            ps, pe = _temp_series_path(scen, depth), _temp_espec_path(scen, depth)
            if not ps.exists() or not pe.exists():
                continue
            ds, de = np.load(ps, allow_pickle=False), np.load(pe, allow_pickle=False)

            for stn in STATIONS:
                t_key, f_key = f"{stn}_t_horas", f"{stn}_freq_cph"
                if t_key not in ds.files or f_key not in de.files:
                    continue
                t_h   = np.asarray(ds[t_key], dtype=float)
                f_cph = np.asarray(de[f_key], dtype=float)

                t_curves, f_curves = [], []
                for case in ("bay", "nobay"):
                    y_key, p_key = f"{stn}_{tag}_{case}", f"{stn}_PSD_{tag}_{case}"
                    if y_key not in ds.files or p_key not in de.files:
                        continue
                    y = np.asarray(ds[y_key], dtype=float)
                    if not _valid(y):
                        continue
                    t_curves.append((case, t_h, y))
                    f_curves.append((case, f_cph, np.asarray(de[p_key], dtype=float)))

                _plot_two_panel(
                    t_curves, f_curves,
                    "T (°C)", "PSD(T) [°C² cph⁻¹]",
                    f"T — bahía vs sin bahía — {scen} — estación {stn} — z≈{depth:.0f} m",
                    f"fig_temp_{scen}_{stn}_{int(depth)}m.png",
                )


# ══ w — serie + PSD (W cruda a profundidad) ════════════════════════════════════

def fig_w():
    for scen in SCENARIOS:
        for depth in DEPTHS:
            tag = _depth_tag(depth, "W")
            ps, pe = _w_series_path(scen, depth), _w_espec_path(scen, depth)
            if not ps.exists() or not pe.exists():
                continue
            ds, de = np.load(ps, allow_pickle=False), np.load(pe, allow_pickle=False)

            for stn in STATIONS:
                t_key, f_key = f"{stn}_t_horas", f"{stn}_freq_cph"
                if t_key not in ds.files or f_key not in de.files:
                    continue
                t_h   = np.asarray(ds[t_key], dtype=float)
                f_cph = np.asarray(de[f_key], dtype=float)

                t_curves, f_curves = [], []
                for case in ("bay", "nobay"):
                    y_key, p_key = f"{stn}_{tag}_{case}", f"{stn}_PSD_{tag}_{case}"
                    if y_key not in ds.files or p_key not in de.files:
                        continue
                    y = np.asarray(ds[y_key], dtype=float)
                    if not _valid(y):
                        continue
                    t_curves.append((case, t_h, y))
                    f_curves.append((case, f_cph, np.asarray(de[p_key], dtype=float)))

                _plot_two_panel(
                    t_curves, f_curves,
                    "w (m s⁻¹)", "PSD(w) [m² s⁻² cph⁻¹]",
                    f"w — bahía vs sin bahía — {scen} — estación {stn} — z≈{depth:.0f} m",
                    f"fig_w_{scen}_{stn}_{int(depth)}m.png",
                )


# ══ rotatorio CW/CCW (1 panel PSD) ═════════════════════════════════════════════

def fig_rotary():
    for scen in SCENARIOS:
        for depth in DEPTHS:
            tag_uv = _depth_tag(depth, "UV")
            pe = _uv_espec_path(scen, depth)
            ps = _uv_series_path(scen, depth)
            if not pe.exists():
                continue
            de = np.load(pe, allow_pickle=False)
            ds = np.load(ps, allow_pickle=False) if ps.exists() else None

            for stn in STATIONS:
                f_key = f"{stn}_freq_cph"
                if f_key not in de.files:
                    continue
                f_cph = np.asarray(de[f_key], dtype=float)

                fig, ax = plt.subplots(figsize=(9, 5))
                has_data = False
                for case in ("bay", "nobay"):
                    cw_key  = f"{stn}_PSD_CW_{tag_uv}_{case}"
                    ccw_key = f"{stn}_PSD_CCW_{tag_uv}_{case}"
                    if cw_key not in de.files or ccw_key not in de.files:
                        continue
                    # filtro de señal válida desde la serie U del caso
                    if ds is not None:
                        u_key = f"{stn}_U{int(round(depth))}m_{case}"
                        if u_key in ds.files and not _valid(ds[u_key]):
                            continue
                    col = CASE_STYLE[case]["color"]
                    ax.plot(f_cph, np.asarray(de[cw_key], dtype=float),
                            color=col, ls="-",  lw=1.4,
                            label=f"{CASE_LABEL[case]} CW")
                    ax.plot(f_cph, np.asarray(de[ccw_key], dtype=float),
                            color=col, ls="--", lw=1.0, alpha=0.75,
                            label=f"{CASE_LABEL[case]} CCW")
                    has_data = True

                if not has_data:
                    plt.close(fig)
                    continue

                _vlines(ax)
                _log_if_positive(ax)
                ax.set_title(f"Espectros rotatorios — bahía vs sin bahía — "
                             f"{scen} — estación {stn} — z≈{depth:.0f} m")
                _finalize_spec(ax, "PSD(u+iv) [m² s⁻² cph⁻¹]")
                fig.tight_layout()
                fig.savefig(FIGS_DIR / f"fig_rotary_{scen}_{stn}_{int(depth)}m.png",
                            dpi=300, bbox_inches="tight")
                plt.close(fig)


# ══ MAIN ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")

    print("fig_eta: η serie + PSD (bahía vs sin bahía)")
    fig_eta()
    print("fig_temp: T serie + PSD (bahía vs sin bahía)")
    fig_temp()
    print("fig_w: w serie + PSD (bahía vs sin bahía)")
    fig_w()
    print("fig_rotary: rotatorios CW/CCW (bahía vs sin bahía)")
    fig_rotary()
    print("Figuras guardadas en:", FIGS_DIR)
