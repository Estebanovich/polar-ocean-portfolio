"""mapas_psd_dominio.py
========================
Mapas 2D de energía espectral por banda sobre el dominio recortado, para:

    Eta (superficie)         → mapas_psd_Eta_sfc.png
    W   (z ≈ 110 m)          → mapas_psd_W_z110m.png
    Temp (z ≈ 110 m)         → mapas_psd_Temp_z110m.png
    UV  (z ≈ 110 m, rotario) → mapas_psd_UV_CCW_z110m.png
                                mapas_psd_UV_CW_z110m.png

Metodología (idéntica a regenerar_espectros_psd.py + Espectros.ipynb Sec. 11):
    - Recorte temporal t ≥ 27 h (descarta transitorio del impulso de viento).
    - Detrend lineal + ventana de Hann por píxel.
    - FFT (rfft para escalares; fft compleja u+iv para UV → CW/CCW).
    - Normalización Parseval: PSD = |Y|² / (U_w · fs) · 3600 → var²·cph⁻¹.
    - Integración en cada banda de frecuencia.

Bandas (basadas en análisis de los *_espectros.npz por estación):
    NI     0.035–0.065 cph   T≈10.5–28.6 h  cuasi-inercial + sub-inercial bajo
    modo 1 0.155–0.223 cph   T≈4.5–6.5 h    modo baroclínico/seiche secundario
    Seiche 0.235–0.290 cph   T≈3.4–4.3 h    seiche barotrópico principal

CLI:
    python mapas_psd_dominio.py                  # genera todo
    python mapas_psd_dominio.py --vars eta w     # solo Eta y W
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.signal import detrend as sp_detrend

_BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(_BASE))

from coords_iwbay import offsets_recentrado

from MIT_DatProcesW import (
    _as_2d_from_grid,
    _close_mnc,
    _is_index_axis,
    _open_mnc_or_raise,
    get_iwbay_paths,
    make_stations_from_grid,
)

# ── constantes ──────────────────────────────────────────────────────────────
NPZ_DIR    = _BASE / "npz_out"
FIGS_DIR   = _BASE / "figs"
TMIN_HORAS = 27.0
# Los ejes del grid.nc NO tienen el origen en la bahía: xgOrigin/ygOrigin son
# números redondos que no coinciden con lo que acumula la malla telescópica, así
# que la costa cae en y ≈ +76.5 km y el eje de la ensenada en x ≈ +76.1 km. Se
# corrige al dibujar (ver coords_iwbay.py); poner False para los ejes crudos.
RECENTRAR  = True
DT_S       = 0.25 * 3600.0     # 15 min
J0, J1     = 200, 352
I0, I1     = 90,  210
J_CENTRO   = 315
RADIO      = 50

BANDS = {
    "NI":     (0.035, 0.065),
    "modo 1": (0.155, 0.223),
    "Seiche": (0.235, 0.290),
}

# Estaciones a NO mostrar en las figuras de mapas PSD.
EXCLUDED_STATIONS = {"Gamma", "Seiche2"}

SCENARIOS = ["2layer", "linear", "realistic"]


# ── IO helpers ──────────────────────────────────────────────────────────────
def _load_zoom(npz_path: Path):
    d = np.load(npz_path, allow_pickle=True)
    t_h  = np.asarray(d["t_horas"], dtype=float)
    Xz   = np.asarray(d["X_zoom"],  dtype=float)
    Yz   = np.asarray(d["Y_zoom"],  dtype=float)
    Mz   = np.asarray(d["mask_zoom"], dtype=bool)
    A    = np.asarray(d["data"],    dtype=np.float32)
    return t_h, Xz, Yz, Mz, A


def _remap_axes_and_stations(Xz, Yz, scen):
    """Si los ejes están en índices, los remapea a metros via grid.
    Devuelve Xz, Yz (metros) y lista de estaciones (xs_m, ys_m, label, color)."""
    paths    = get_iwbay_paths()
    grid_dir = Path(paths["scenarios"][scen]) / paths["bay_dir"]
    gr       = _open_mnc_or_raise(grid_dir, "grid.t*.nc")
    try:
        depth = np.array(gr.variables["Depth"][:], dtype=float)
        ny, nx = depth.shape
        XC = _as_2d_from_grid(gr, ("XC","xC","X","x","XG","xG"), axis="x", ny=ny, nx=nx)
        YC = _as_2d_from_grid(gr, ("YC","yC","Y","y","YG","yG"), axis="y", ny=ny, nx=nx)
        stations_idx = make_stations_from_grid(gr, j_centro=J_CENTRO, radio=RADIO)

        if _is_index_axis(Xz, nx) and _is_index_axis(Yz, ny):
            ix = np.clip(np.round(Xz).astype(int), 0, nx - 1)
            iy = np.clip(np.round(Yz).astype(int), 0, ny - 1)
            Xz = XC[iy[0], ix].astype(float)
            Yz = YC[iy, ix[0]].astype(float)

        st_phys = []
        for (i, j, lab, col) in stations_idx:
            if str(lab) in EXCLUDED_STATIONS:
                continue
            i, j = int(i), int(j)
            st_phys.append((float(XC[j, i]), float(YC[j, i]),
                            str(lab), str(col or "k")))
    finally:
        _close_mnc(gr)
    return Xz, Yz, st_phys


# ── cómputo PSD por píxel ──────────────────────────────────────────────────
def _psd_scalar(A_cut):
    """A_cut: (nt, ny, nx). Devuelve (f_cph, PSD shape (nf, ny*nx))."""
    nt, nyc, nxc = A_cut.shape
    A2d = sp_detrend(A_cut.reshape(nt, nyc * nxc).astype(np.float64),
                     axis=0, type="linear")
    win = np.hanning(nt)
    U2  = float((win ** 2).sum())
    A2d *= win[:, np.newaxis]
    Af    = np.fft.rfft(A2d, axis=0)
    f_cph = np.fft.rfftfreq(nt, d=DT_S) * 3600.0
    PSD   = (np.abs(Af) ** 2) / (U2 * (1.0 / DT_S)) * 3600.0
    return f_cph, PSD


def _psd_rotary(U_cut, V_cut):
    """Espectro rotario por píxel. Devuelve (f_pos_cph, PSD_CCW, PSD_CW)
    cada uno shape (nf_pos, ny*nx)."""
    nt, nyc, nxc = U_cut.shape
    Ud = sp_detrend(U_cut.reshape(nt, nyc * nxc).astype(np.float64),
                    axis=0, type="linear")
    Vd = sp_detrend(V_cut.reshape(nt, nyc * nxc).astype(np.float64),
                    axis=0, type="linear")
    win = np.hanning(nt)
    U2  = float((win ** 2).sum())
    Z   = (Ud + 1j * Vd) * win[:, np.newaxis]
    Zf  = np.fft.fft(Z, axis=0)
    f   = np.fft.fftfreq(nt, d=DT_S) * 3600.0
    PSD = (np.abs(Zf) ** 2) / (U2 * (1.0 / DT_S)) * 3600.0
    pos = np.where(f > 0)[0]
    neg = np.where(f < 0)[0][::-1]   # |f| creciente
    return f[pos], PSD[pos], PSD[neg]


def _integrate_bands(f_cph, PSD, mask_2d, shape_yx):
    """PSD: (nf, ny*nx). Devuelve dict[band] -> 2D (ny, nx) con NaN en tierra."""
    if len(f_cph) < 2:
        raise ValueError("Frecuencia insuficiente para integrar.")
    df = float(f_cph[1] - f_cph[0])
    nyc, nxc = shape_yx
    out = {}
    for tag, (f0, f1) in BANDS.items():
        bm = (f_cph >= f0) & (f_cph <= f1)
        E  = (PSD[bm] * df).sum(axis=0).reshape(nyc, nxc)
        out[tag] = np.where(~mask_2d, E, np.nan)
    return out


# ── pipelines por variable ─────────────────────────────────────────────────
def maps_scalar(npz_path, scen):
    t_h, Xz, Yz, Mz, A = _load_zoom(npz_path)
    Xz, Yz, st_phys    = _remap_axes_and_stations(Xz, Yz, scen)
    mt = t_h >= TMIN_HORAS
    Ac = A[mt][:, J0:J1, I0:I1]
    mc = Mz[J0:J1, I0:I1]
    Xc = Xz[I0:I1]; Yc = Yz[J0:J1]
    f, PSD = _psd_scalar(Ac)
    return _integrate_bands(f, PSD, mc, Ac.shape[1:]), Xc, Yc, st_phys


def maps_rotary(u_npz, v_npz, scen):
    t_h, Xz, Yz, Mz, U = _load_zoom(u_npz)
    _,   _,  _,  _,  V = _load_zoom(v_npz)
    Xz, Yz, st_phys    = _remap_axes_and_stations(Xz, Yz, scen)
    mt = t_h >= TMIN_HORAS
    Uc = U[mt][:, J0:J1, I0:I1]
    Vc = V[mt][:, J0:J1, I0:I1]
    mc = Mz[J0:J1, I0:I1]
    Xc = Xz[I0:I1]; Yc = Yz[J0:J1]
    f, P_ccw, P_cw = _psd_rotary(Uc, Vc)
    shape_yx = Uc.shape[1:]
    return (_integrate_bands(f, P_ccw, mc, shape_yx),
            _integrate_bands(f, P_cw,  mc, shape_yx),
            Xc, Yc, st_phys)


# ── plotting ────────────────────────────────────────────────────────────────
def _band_color_limits(all_maps, scen_order, tag):
    vals = []
    for s in scen_order:
        v = all_maps[s][tag]
        v = v[np.isfinite(v) & (v > 0)]
        if v.size: vals.append(v)
    if not vals:
        return 1e-12, 1e-3
    cat  = np.concatenate(vals)
    vmax = float(np.percentile(cat, 99.5))
    vmin = max(float(np.percentile(cat, 5)), vmax * 1e-5)
    return vmin, vmax


def plot_maps(all_maps, all_xy, scen_order, var_label, units, out_path,
              cmap="plasma", title_extra="", land_color="#ece2c6",
              band_lim=None):
    bands = list(BANDS.keys())
    nrows = len(scen_order); ncols = len(bands)

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.6 * ncols, 4.0 * nrows),
        gridspec_kw={"wspace": 0.08, "hspace": 0.12},
    )
    axes = np.atleast_2d(axes)

    # band_lim externo (p.ej. compartido entre bahía/sin-bahía) o auto por figura
    if band_lim is None:
        band_lim = {t: _band_color_limits(all_maps, scen_order, t) for t in bands}

    # colormap con celdas NaN (tierra) coloreadas tipo arena
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(land_color)

    x0, y0 = offsets_recentrado() if RECENTRAR else (0.0, 0.0)
    for r, scen in enumerate(scen_order):
        Xkm, Ykm, st_phys = all_xy[scen]
        Xkm, Ykm = np.asarray(Xkm, float) - x0, np.asarray(Ykm, float) - y0
        for c, tag in enumerate(bands):
            ax = axes[r, c]
            ax.set_facecolor(land_color)
            E  = all_maps[scen][tag]
            vmin, vmax = band_lim[tag]
            pcm = ax.pcolormesh(
                Xkm, Ykm, E,
                norm=mcolors.LogNorm(vmin=vmin, vmax=vmax),
                cmap=cmap_obj, shading="auto", rasterized=True,
            )
            fig.colorbar(pcm, ax=ax, label=f"E [{units}]",
                         pad=0.01, fraction=0.04)

            for (xs, ys, lab, col) in st_phys:
                xs_km, ys_km = xs / 1e3 - x0, ys / 1e3 - y0
                if not (Xkm.min() <= xs_km <= Xkm.max()
                        and Ykm.min() <= ys_km <= Ykm.max()):
                    continue
                ax.plot(xs_km, ys_km, marker="o", color=col,
                        ms=6, mec="k", mew=0.5, zorder=5)
                if c == 0:
                    ax.annotate(lab, (xs_km, ys_km), fontsize=7, color="w",
                                xytext=(4, 2), textcoords="offset points")

            ax.set_xlim(Xkm.min(), Xkm.max())
            ax.set_ylim(Ykm.min(), Ykm.max())
            ax.set_aspect("equal")
            ax.tick_params(labelsize=8)

            if r == 0:
                f0, f1 = BANDS[tag]
                T0, T1 = 1.0 / f1, 1.0 / f0
                ax.set_title(f"{tag}  {f0}–{f1} cph  (T≈{T0:.1f}–{T1:.1f} h)",
                             fontsize=10, fontweight="bold")
            if r == nrows - 1:
                ax.set_xlabel("X (km)", fontsize=9)
            else:
                ax.tick_params(labelbottom=False)
            if c == 0:
                ax.set_ylabel(f"{scen}\nY (km)", fontsize=9)
            else:
                ax.tick_params(labelleft=False)

    origen = ("  •  0 = eje de la ensenada (X) y línea de costa (Y)"
              if RECENTRAR else "  •  ejes crudos del grid.nc")
    fig.suptitle(
        f"Energía espectral 2D — {var_label}  {title_extra}\n"
        f"t > {TMIN_HORAS:g} h  •  Welch/Parseval + Hann{origen}",
        fontsize=12, y=1.0,
    )
    FIGS_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


# ── orquestador ─────────────────────────────────────────────────────────────
def run(vars_to_run):
    if "eta" in vars_to_run:
        m_all, xy_all = {}, {}
        for s in SCENARIOS:
            m, X, Y, st = maps_scalar(NPZ_DIR / f"{s}_Eta_zoom_sfc.npz", s)
            m_all[s] = m; xy_all[s] = (X / 1e3, Y / 1e3, st)
            print(f"  [Eta {s}] max por banda: " +
                  " ".join(f"{t}={np.nanmax(v):.2e}" for t, v in m.items()))
        plot_maps(m_all, xy_all, SCENARIOS, "Δη", "m²·cph",
                  FIGS_DIR / "mapas_psd_Eta_sfc.png",
                  cmap="viridis", title_extra="@ superficie")

    if "w" in vars_to_run:
        m_all, xy_all = {}, {}
        for s in SCENARIOS:
            m, X, Y, st = maps_scalar(NPZ_DIR / f"{s}_W_zoom_z110m.npz", s)
            m_all[s] = m; xy_all[s] = (X / 1e3, Y / 1e3, st)
            print(f"  [W {s}] max por banda: " +
                  " ".join(f"{t}={np.nanmax(v):.2e}" for t, v in m.items()))
        plot_maps(m_all, xy_all, SCENARIOS, "ΔW", "m²·s⁻²·cph",
                  FIGS_DIR / "mapas_psd_W_z110m.png",
                  cmap="plasma", title_extra="@ z=110 m")

    if "temp" in vars_to_run:
        m_all, xy_all = {}, {}
        for s in SCENARIOS:
            m, X, Y, st = maps_scalar(NPZ_DIR / f"{s}_Temp_zoom_z110m.npz", s)
            m_all[s] = m; xy_all[s] = (X / 1e3, Y / 1e3, st)
            print(f"  [Temp {s}] max por banda: " +
                  " ".join(f"{t}={np.nanmax(v):.2e}" for t, v in m.items()))
        plot_maps(m_all, xy_all, SCENARIOS, "ΔTemp", "°C²·cph",
                  FIGS_DIR / "mapas_psd_Temp_z110m.png",
                  cmap="inferno", title_extra="@ z=110 m")

    if "uv" in vars_to_run:
        ccw_all, cw_all, xy_all = {}, {}, {}
        for s in SCENARIOS:
            ccw, cw, X, Y, st = maps_rotary(
                NPZ_DIR / f"{s}_U_zoom_z110m.npz",
                NPZ_DIR / f"{s}_V_zoom_z110m.npz", s,
            )
            ccw_all[s] = ccw; cw_all[s] = cw
            xy_all[s] = (X / 1e3, Y / 1e3, st)
            print(f"  [UV-CCW {s}] " +
                  " ".join(f"{t}={np.nanmax(v):.2e}" for t, v in ccw.items()))
            print(f"  [UV-CW  {s}] " +
                  " ".join(f"{t}={np.nanmax(v):.2e}" for t, v in cw.items()))
        plot_maps(ccw_all, xy_all, SCENARIOS,
                  "UV rotario CCW (antihorario)", "m²·s⁻²·cph",
                  FIGS_DIR / "mapas_psd_UV_CCW_z110m.png",
                  cmap="plasma", title_extra="@ z=110 m")
        plot_maps(cw_all, xy_all, SCENARIOS,
                  "UV rotario CW (horario)", "m²·s⁻²·cph",
                  FIGS_DIR / "mapas_psd_UV_CW_z110m.png",
                  cmap="plasma", title_extra="@ z=110 m")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vars", nargs="+",
                    default=["eta", "w", "temp", "uv"],
                    choices=["eta", "w", "temp", "uv"],
                    help="Variables a procesar (default: todas).")
    args = ap.parse_args()
    run(set(args.vars))
