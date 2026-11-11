"""mapas_psd_dominio_bay_vs_nobay.py
=====================================
Igual que `mapas_psd_dominio.py` pero para visualizar DÓNDE viven las
frecuencias (bandas NI / modo 1 / Seiche) en el dominio para CADA caso por
separado: CON bahía (run_expand) y SIN bahía (run_expand_nobay) — no la
diferencia.

Lee los cubos zoom separados que produce `regenerar_zoom_bay_nobay.py`
(`{scen}_{var}_zoom_{keyz}_{bay|nobay}.npz`) y genera, por variable, DOS
figuras (una por caso), filas = escenarios, columnas = bandas:

    mapas_psd_Eta_sfc_{bay,nobay}.png
    mapas_psd_W_z110m_{bay,nobay}.png
    mapas_psd_Temp_z110m_{bay,nobay}.png
    mapas_psd_UV_CCW_z110m_{bay,nobay}.png
    mapas_psd_UV_CW_z110m_{bay,nobay}.png

Nota: la máscara de tierra proviene del grid CON bahía; en el caso sin-bahía
la región de la ensenada aparece como energía ~0 (sin señal), no coloreada
como tierra.

Uso:
    python mapas_psd_dominio_bay_vs_nobay.py
    python mapas_psd_dominio_bay_vs_nobay.py --vars eta w --source bay
"""
import argparse
import sys
from pathlib import Path

_BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(_BASE))

# Reutiliza toda la maquinaria de PSD/bandas/plot del script original.
from mapas_psd_dominio import (
    SCENARIOS, NPZ_DIR, FIGS_DIR, BANDS,
    maps_scalar, maps_rotary, plot_maps, _band_color_limits,
)

CASE_TITLE = {"bay": "con bahía", "nobay": "sin bahía"}


def _shared_band_lim(maps_by_case):
    """Límites de color por banda compartidos entre casos para comparación justa.
    maps_by_case: {case: {scen: {band: 2D}}}. Devuelve {band: (vmin, vmax)}."""
    merged, order = {}, []
    for case, m_all in maps_by_case.items():
        for s, m in m_all.items():
            key = f"{s}|{case}"
            merged[key] = m
            order.append(key)
    return {t: _band_color_limits(merged, order, t) for t in BANDS}


def _scalar_var(sources, varname, keyz, var_label, units, cmap):
    """Mapas escalares de ambos casos con escala de color compartida."""
    maps_by_case, xy_by_case, order_by_case = {}, {}, {}
    for case in sources:
        m_all, xy_all, order = {}, {}, []
        for s in SCENARIOS:
            p = NPZ_DIR / f"{s}_{varname}_zoom_{keyz}_{case}.npz"
            if not p.exists():
                print(f"  [skip] {p.name} no existe")
                continue
            m, X, Y, st = maps_scalar(p, s)
            m_all[s] = m
            xy_all[s] = (X / 1e3, Y / 1e3, st)
            order.append(s)
        if order:
            maps_by_case[case] = m_all
            xy_by_case[case]   = xy_all
            order_by_case[case] = order
    if not maps_by_case:
        return
    band_lim = _shared_band_lim(maps_by_case)
    for case, m_all in maps_by_case.items():
        cl = CASE_TITLE[case]
        plot_maps(m_all, xy_by_case[case], order_by_case[case],
                  f"{var_label} ({cl})", units,
                  FIGS_DIR / f"mapas_psd_{varname}_{keyz}_{case}.png",
                  cmap=cmap, title_extra=cl, band_lim=band_lim)


def _rotary_var(sources, keyz, units):
    """Mapas rotatorios CW/CCW de ambos casos con escala compartida (por sentido)."""
    ccw_by_case, cw_by_case, xy_by_case, order_by_case = {}, {}, {}, {}
    for case in sources:
        ccw_all, cw_all, xy_all, order = {}, {}, {}, []
        for s in SCENARIOS:
            pu = NPZ_DIR / f"{s}_U_zoom_{keyz}_{case}.npz"
            pv = NPZ_DIR / f"{s}_V_zoom_{keyz}_{case}.npz"
            if not pu.exists() or not pv.exists():
                print(f"  [skip] {pu.name}/{pv.name} no existen")
                continue
            ccw, cw, X, Y, st = maps_rotary(pu, pv, s)
            ccw_all[s] = ccw; cw_all[s] = cw
            xy_all[s]  = (X / 1e3, Y / 1e3, st)
            order.append(s)
        if order:
            ccw_by_case[case] = ccw_all; cw_by_case[case] = cw_all
            xy_by_case[case]  = xy_all;  order_by_case[case] = order
    if not ccw_by_case:
        return
    lim_ccw = _shared_band_lim(ccw_by_case)
    lim_cw  = _shared_band_lim(cw_by_case)
    for case in ccw_by_case:
        cl = CASE_TITLE[case]
        plot_maps(ccw_by_case[case], xy_by_case[case], order_by_case[case],
                  f"UV rotario CCW ({cl})", units,
                  FIGS_DIR / f"mapas_psd_UV_CCW_{keyz}_{case}.png",
                  cmap="plasma", title_extra=cl, band_lim=lim_ccw)
        plot_maps(cw_by_case[case], xy_by_case[case], order_by_case[case],
                  f"UV rotario CW ({cl})", units,
                  FIGS_DIR / f"mapas_psd_UV_CW_{keyz}_{case}.png",
                  cmap="plasma", title_extra=cl, band_lim=lim_cw)


def run(vars_to_run, sources):
    if "eta" in vars_to_run:
        _scalar_var(sources, "Eta", "sfc", "η", "m²·cph", "viridis")
    if "w" in vars_to_run:
        _scalar_var(sources, "W", "z110m", "W", "m²·s⁻²·cph", "plasma")
    if "temp" in vars_to_run:
        _scalar_var(sources, "Temp", "z110m", "Temp", "°C²·cph", "inferno")
    if "uv" in vars_to_run:
        _rotary_var(sources, "z110m", "m²·s⁻²·cph")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vars", nargs="+", default=["eta", "w", "temp", "uv"],
                    choices=["eta", "w", "temp", "uv"])
    ap.add_argument("--source", nargs="+", default=["bay", "nobay"],
                    choices=["bay", "nobay"])
    args = ap.parse_args()
    run(set(args.vars), args.source)
    print("Figuras guardadas en:", FIGS_DIR)
