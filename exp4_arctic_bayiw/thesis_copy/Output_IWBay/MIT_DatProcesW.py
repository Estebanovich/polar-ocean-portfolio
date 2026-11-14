"""
MIT_DatProcesW — Post-processing utilities for the IWBay MITgcm experiments.

Provides helpers for loading MNC/NetCDF model output, extracting transects and
station time-series, performing barotropic/baroclinic decompositions, and
generating figures and animations for the internal-wave bay thesis.
"""
from collections import namedtuple
from pathlib import Path
import glob
import json
import math
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from scipy.signal import detrend, get_window
import cmocean as cmo
import MITgcmutils as MIT
from netCDF4 import Dataset

# Default station-geometry parameters used throughout the module
_J_CENTRO_DEFAULT: int = 315
_RADIO_DEFAULT: int = 50

# Named tuple for observation stations (grid indices + metadata)
Station = namedtuple('Station', ['i', 'j', 'label', 'color'])

def get_iwbay_paths():
    """
    Devuelve rutas base y escenarios estándar para IWBay.

    Retorna
    -------
    dict con llaves:
      - base_linear: Path
      - base_karina: Path
      - scenarios: dict[str, Path]
      - bay_dir: str
      - nobay_dir: str
    """
    base_linear = Path("/Volumes/Esteban_SSD/Maestría/MITgcm/verification")
    base_karina = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification")

    scenarios = {
        "realistic": base_karina / "BayIW_Rectan",
        "linear": base_linear / "BayIW_Rectan_linear",
        "2layer": base_karina / "BayIW_Rectan_2layer",
        "linear_2xwind": base_karina / "BayIW_Rectan_linear_2xwind",
    }

    return {
        "base_linear": base_linear,
        "base_karina": base_karina,
        "scenarios": scenarios,
        "bay_dir": "run_expand/mnc_00*",
        "nobay_dir": "run_expand_nobay/mnc_00*",
    }

def set_plot_style_publication(
    context: str = "talk",
    style: str = "whitegrid",
    font_scale: float = 0.9,
    despine: bool = False
):
    """
    Aplica un estilo consistente para figuras de notebook/publicación.
    - Usa seaborn si está disponible.
    - Ajusta rcParams de Matplotlib para paneles múltiples (hovmöller/modos).
    """
    try:
        import seaborn as sns
        sns.set_theme(context=context, style=style, font_scale=font_scale)
        if despine:
            sns.despine()
    except Exception:
        # Fallback limpio si seaborn no está instalado.
        plt.style.use("default")

    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": ":",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

# ==== Helpers de concatenación tiempo/datos (si no existen) ====

def _get_time_from_nc(nc):
    if 'T' in nc.variables:
        return np.array(nc.variables['T'][:], dtype=float)
    if 'time' in nc.variables:
        return np.array(nc.variables['time'][:], dtype=float)
    raise KeyError("No encontré variable de tiempo ('T' o 'time').")

def _concat_time_only(times_list, tol=1e-9):
    """
    Concatena una lista de vectores de tiempo [t1, t2, ...] eliminando el
    primer elemento del siguiente bloque si duplica al último del acumulado.
    """
    out = []
    for k, tk in enumerate(times_list):
        if tk is None or np.size(tk) == 0:
            continue
        tk = np.asarray(tk, dtype=float)
        if k > 0 and out and abs(tk[0] - out[-1][-1]) <= tol:
            tk = tk[1:]
        out.append(tk)
    return np.concatenate(out) if out else np.array([])

def _concat_time_and_data(times_list, data_list, tol=1e-9):
    """
    Concatena tiempos [t1, t2, ...] y matrices/arrays [X1, X2, ...] a lo largo
    del eje de tiempo (0), eliminando el primer registro del siguiente bloque
    si duplica al último del acumulado (|Δt|<=tol).

    times_list: list de arrays shape (nt,)
    data_list : list de arrays shape (nt, ...), por ej. (nt, nz) o (nt, nz, ny, nx)
    Devuelve (T_concat, X_concat).
    """
    t_all, X_all = [], []
    for k, (t, X) in enumerate(zip(times_list, data_list)):
        if t is None or np.size(t) == 0 or X is None or np.size(X) == 0:
            continue
        t = np.asarray(t, dtype=float)
        X = np.asarray(X)
        if k > 0 and t_all and abs(t[0] - t_all[-1][-1]) <= tol:
            t = t[1:]
            X = X[1:, ...]
        t_all.append(t)
        X_all.append(X)
    T = np.concatenate(t_all) if t_all else np.array([])
    X = np.concatenate(X_all, axis=0) if X_all else np.array([])
    return T, X
# ==== fin helpers ====


def _get_RC(gr):
    """Return the 1-D vertical cell-centre coordinate array (negative metres) from a grid file."""
    if 'RC' in gr.variables:
        RC = np.array(gr.variables['RC'][:], dtype=float).squeeze()
    elif 'Z' in gr.variables:
        RC = np.array(gr.variables['Z'][:], dtype=float).squeeze()
    else:
        raise KeyError("No encontré coordenadas verticales ('RC' o 'Z') en el grid.")
    if RC.ndim > 1:
        RC = RC[..., 0, 0]
    return RC


def _to_1d_x(gr, nx, varname_primary='XC', varname_fallback='XG'):
    names = [varname_primary, varname_fallback, "xC", "X", "x", "xG"]
    seen = set()
    for vn in names:
        if vn in seen:
            continue
        seen.add(vn)
        if vn in gr.variables:
            arr = np.array(gr.variables[vn][:], dtype=float)
            if arr.size > 0:
                if arr.ndim == 2:
                    return arr[0, :]
                if arr.ndim == 1:
                    if arr.size == nx + 1:
                        return arr[:-1]
                    return arr
    return np.arange(nx, dtype=float)

def _to_1d_y(gr, ny, varname_primary='YC', varname_fallback='YG'):
    names = [varname_primary, varname_fallback, "yC", "Y", "y", "yG"]
    seen = set()
    for vn in names:
        if vn in seen:
            continue
        seen.add(vn)
        if vn in gr.variables:
            arr = np.array(gr.variables[vn][:], dtype=float)
            if arr.size > 0:
                if arr.ndim == 2:
                    return arr[:, 0]
                if arr.ndim == 1:
                    if arr.size == ny + 1:
                        return arr[:-1]
                    return arr
    return np.arange(ny, dtype=float)

def _as_2d_from_grid(gr, names, axis, ny, nx):
    arr = None
    for nm in names:
        if nm in gr.variables:
            arr = np.array(gr.variables[nm][:], dtype=float)
            break
    if arr is None:
        return None
    if arr.ndim == 2:
        return arr
    if arr.ndim == 1:
        if axis == "x":
            if arr.size == nx + 1:
                arr = arr[:-1]
            if arr.size == nx:
                return np.tile(arr[None, :], (ny, 1))
        else:
            if arr.size == ny + 1:
                arr = arr[:-1]
            if arr.size == ny:
                return np.tile(arr[:, None], (1, nx))
    return None

def _is_index_axis(v, n_expected):
    vv = np.asarray(v, dtype=float).ravel()
    if vv.size < 2 or not np.all(np.isfinite(vv)):
        return False
    int_like = float(np.max(np.abs(vv - np.rint(vv)))) < 1e-6
    in_range = (float(np.min(vv)) >= -1.0) and (float(np.max(vv)) <= float(n_expected) + 1.0)
    mono = np.all(np.diff(vv) >= -1e-9)
    return bool(int_like and in_range and mono)


def _open_mnc_or_raise(dirpath: Path, pattern: str):
    """
    Abre archivos MNC que coincidan con ``dirpath / pattern``.

    Si no encuentra nada directamente, reintenta dentro de los
    subdirectorios ``mnc_00*/`` (salida estándar de MITgcm con MNC).
    """
    search = str(dirpath / pattern)
    if not glob.glob(search):
        # Intentar dentro de subdirectorios MNC (mnc_0001/, mnc_0002/, ...)
        search = str(dirpath / "mnc_00*" / pattern)
        if not glob.glob(search):
            raise FileNotFoundError(
                f"No encontré '{pattern}' en:\n"
                f"  {dirpath}\n"
                f"  {dirpath}/mnc_00*/"
            )
    return MIT.mnc_files(search)


def discover_stage_patterns(outdir, var_prefix: str = "state"):
    """
    Descubre automáticamente los prefijos de iteración de los archivos MNC
    ``<var_prefix>.<iter>.t*.nc`` en ``outdir`` (o ``outdir/mnc_00*/``) y devuelve
    la tupla de patrones ordenada cronológicamente por número de iteración:
    ``("<var_prefix>.0000000000.t*.nc", "<var_prefix>.0000001260.t*.nc", ...)``.

    Funciona tanto con salida monolítica (p.ej. 2 etapas: 0000000000 y 0000007200)
    como con salida troceada en muchos chunks (p.ej. cada 1260 iteraciones), por lo
    que reemplaza con seguridad al antiguo ``STAGE_PATTERNS`` fijo de 2 elementos.
    """
    outdir = Path(outdir)
    search = str(outdir / f"{var_prefix}.*.t*.nc")
    files = glob.glob(search)
    if not files:
        search = str(outdir / "mnc_00*" / f"{var_prefix}.*.t*.nc")
        files = glob.glob(search)
    pat = re.compile(rf"^{re.escape(var_prefix)}\.(\d+)\.t")
    iters = set()
    for f in files:
        m = pat.match(Path(f).name)
        if m:
            iters.add(int(m.group(1)))
    if not iters:
        raise FileNotFoundError(
            f"No encontré archivos '{var_prefix}.<iter>.t*.nc' en:\n"
            f"  {outdir}\n  {outdir}/mnc_00*/"
        )
    return tuple(f"{var_prefix}.{it:010d}.t*.nc" for it in sorted(iters))


def cargar_datos_mitgcm_etapas(outdir_bay: Path,
                               outdir_nobay: Path,
                               stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
                               grid_pattern   = "grid.t*.nc"):
    gr_bay = _open_mnc_or_raise(outdir_bay, grid_pattern)
    nc_bay_list, nc_nby_list = [], []
    tb_list, tn_list = [], []

    for patt in stage_patterns:
        nc_b = _open_mnc_or_raise(outdir_bay, patt)
        nc_n = _open_mnc_or_raise(outdir_nobay, patt)
        nc_bay_list.append(nc_b)
        nc_nby_list.append(nc_n)
        tb_list.append(_get_time_from_nc(nc_b))
        tn_list.append(_get_time_from_nc(nc_n))

    time_concat_b = _concat_time_only(tb_list)
    time_concat_n = _concat_time_only(tn_list)
    if time_concat_b.size and time_concat_n.size:
        if time_concat_b.shape != time_concat_n.shape or np.any(np.abs(time_concat_b - time_concat_n) > 1e-7):
            print("[AVISO] Ejes de tiempo difieren entre bay/nobay; se usará el de 'bay' para recortar.")

    return {"nc_bay_list": nc_bay_list,
            "nc_nobay_list": nc_nby_list,
            "gr_bay": gr_bay,
            "time": time_concat_b}

def _close_mnc(mnc_obj):
    # cierra todos los archivos abiertos dentro del objeto MNC
    for fh in getattr(mnc_obj, "nc", []):
        try:
            fh.close()
        except Exception:
            pass

# ─────────────────────────────────────────────────────────────
# HELPER: encontrar la profundidad de la isoterma en un perfil
# ─────────────────────────────────────────────────────────────
def _isoterma_depth_1col(Tz, z, T_iso):
    """
    Encuentra la profundidad (positiva hacia abajo) donde T(z)=T_iso.
    Si no se detecta cruce exacto, usa el punto de temperatura más cercano.
    Si el perfil no permite estimar profundidad, devuelve np.nan.

    Parámetros
    ----------
    Tz : array (nz,)
        Perfil vertical de temperatura [°C]
    z : array (nz,)
        Coordenadas verticales del MITgcm (negativas hacia abajo)
    T_iso : float
        Temperatura objetivo de la isoterma [°C]

    Retorna
    -------
    depth : float
        Profundidad positiva en metros (np.nan si no se estima).
    """
    nz = Tz.size
    # Convertir a profundidad positiva
    depth = -z  # z<0 → profundidad>0

    # Buscar cruce entre niveles
    for k in range(nz - 1):
        Tk, Tk1 = Tz[k], Tz[k+1]

        # Cambio de signo → cruce
        if (Tk - T_iso) * (Tk1 - T_iso) <= 0:
            # Interpolación lineal
            if Tk1 != Tk:
                alpha = (T_iso - Tk) / (Tk1 - Tk)
                return depth[k] + alpha * (depth[k+1] - depth[k])
            else:
                # Ambos niveles son exactamente la isoterma
                return depth[k]

    # Si llegamos aquí, no hubo cruce exacto.
    # Buscar el nivel más cercano a T_iso.
    idx = np.argmin(np.abs(Tz - T_iso))
    if np.isfinite(Tz[idx]):
        return depth[idx]

    return np.nan



# ──────────────────────────
# Selección de nivel por profundidad
# ──────────────────────────
def _select_k_for_depth(gr_bay, depth_m: float, prefer_RF=False):
    """
    depth_m > 0 (profundidad bajo la superficie).
    Usa RC (centros) por defecto. Si prefer_RF=True, usa RF (caras).
    Devuelve (k_index, z_ref) donde z_ref es la profundidad (negativa) de ese nivel.
    """
    depth_neg = -abs(depth_m)
    if prefer_RF and 'RF' in gr_bay.variables:
        z = np.array(gr_bay.variables['RF'][:], dtype=float).squeeze()
    elif 'RC' in gr_bay.variables:
        z = np.array(gr_bay.variables['RC'][:], dtype=float).squeeze()
    elif 'Z' in gr_bay.variables:
        z = np.array(gr_bay.variables['Z'][:], dtype=float).squeeze()
    else:
        raise KeyError("No encontré coordenadas verticales ('RC'/'RF'/'Z') en grid.")

    # Si z es 1D (nz,), ok. Si es 3D, tomamos columna 1D
    if z.ndim > 1:
        z = z[..., 0, 0] if z.shape[-2:] != () else z.squeeze()

    k = int(np.argmin(np.abs(z - depth_neg)))
    return k, float(z[k])

# ──────────────────────────
# Lectura y concatenación de W en (k,j,i)
# ──────────────────────────
def _read_var_concat_3d(nc_list, varname: str, k: int, j: int, i: int):
    times, series = [], []
    for nc in nc_list:
        t = _get_time_from_nc(nc)
        arr = np.array(nc.variables[varname][:, k, j, i], dtype=float)
        if t.size == 0 or arr.size == 0:
            continue
        if arr.shape[0] != t.shape[0]:
            raise ValueError(f"Inconsistencia {varname}: longitudes no coinciden en una etapa.")
        times.append(t); series.append(arr)
    if not times:
        return np.array([]), np.array([])
    return _concat_time_and_data(times, series)
#def _read_var_concat_

# ──────────────────────────
# Guardar ΔW a una lista de profundidades
# ──────────────────────────
def guardar_W_series_y_espectros(
    nc_bay_list,        # lista de MNC por etapas (de cargar_datos_mitgcm_etapas)
    nc_nobay_list,      # lista de MNC por etapas
    gr_bay,
    time_concat,        # tiempo concatenado (de bay)
    stations,           # [(i,j,label)] o [(i,j,label,color)]
    depths_m=(110.0, 250.0),
    tmin_horas: float = 27.0,
    save_path: Path = Path("./npz_out"),
    save_prefix_base: str = "escenario_W",
    extra_meta: dict = None,
    prefer_RF: bool = False,     # True si prefieres usar RF para W
):
    """
    Para cada profundidad en 'depths_m', guarda:
      - <save_prefix_base>_W<depth>m_series.npz
      - <save_prefix_base>_W<depth>m_espectros.npz
    Claves:  '<label>_t_horas', '<label>_W<depth>m'  y  '<label>_freq_cph', '<label>_PSD_W<depth>m'
    """
    save_path.mkdir(parents=True, exist_ok=True)

    # recorte temporal
    t = np.array(time_concat, dtype=float)
    tmin_seg = tmin_horas * 3600.0
    idx_start = np.where(t >= tmin_seg)[0][0] if np.any(t >= tmin_seg) else 0
    t_cut = t[idx_start:]
    t_cut_h = t_cut / 3600.0

    for d in depths_m:
        # Nivel k más cercano a la profundidad
        k, z_ref = _select_k_for_depth(gr_bay, d, prefer_RF=prefer_RF)

        series_npz, espec_npz, tabla_est = {}, {}, []
        for st in stations:
            i_st, j_st, label, color = st

            # W en bay/nobay concatenado por etapas
            tb, Wb = _read_var_concat_3d(nc_bay_list,  "W", k, j_st, i_st)
            tn, Wn = _read_var_concat_3d(nc_nobay_list, "W", k, j_st, i_st)
            if tb.size == 0 or tn.size == 0:
                raise ValueError(f"No se pudo leer W para estación {label} a {d} m.")

            # recortar según índice de 't' (tiempo de bay)
            Wb_cut = Wb[idx_start:]
            Wn_cut = Wn[idx_start:]
            dW_cut = Wb_cut - Wn_cut

            # guardar serie
            keyW = f"W{int(round(d))}m"
            series_npz[f"{label}_t_horas"]   = t_cut_h
            series_npz[f"{label}_{keyW}"]    = dW_cut

            # espectro Welch/Parseval con ventana Hann
            if len(t_cut) > 1:
                dt  = t_cut[1] - t_cut[0]           # segundos
                N   = len(dW_cut)
                win = get_window('hann', N, fftbins=True)
                U_w = float((win**2).sum())
                eta_detr = detrend(dW_cut) * win
                Y    = np.fft.rfft(eta_detr)
                freqs = np.fft.rfftfreq(N, d=dt)
                PSD  = (np.abs(Y)**2) / (U_w * (1.0 / dt)) * 3600.0
                freq_cph = freqs * 3600.0
            else:
                PSD = np.array([np.nan]); freq_cph = np.array([np.nan])

            espec_npz[f"{label}_freq_cph"]     = freq_cph
            espec_npz[f"{label}_PSD_{keyW}"]   = PSD

            tabla_est.append({"label": str(label), "i": int(i_st), "j": int(j_st), "color": color})

        # metadatos
        meta = {
            "var": "W",
            "depth_m": float(d),
            "z_ref": float(z_ref),    # profundidad negativa del nivel elegido
            "k_index": int(k),
            "tmin_horas": float(tmin_horas),
            "idx_start": int(idx_start),
            "dt_seconds": float(t_cut[1]-t_cut[0]) if len(t_cut) > 1 else None,
            "stations": tabla_est,
            "prefer_RF": bool(prefer_RF)
        }
        if extra_meta: meta.update(extra_meta)
        meta_json = json.dumps(meta, ensure_ascii=False)

        series_npz["meta_json"] = np.array(meta_json)
        espec_npz["meta_json"]  = np.array(meta_json)

        keyW = f"W{int(round(d))}m"
        out_series    = save_path / f"{save_prefix_base}_{keyW}_series.npz"
        out_espectros = save_path / f"{save_prefix_base}_{keyW}_espectros.npz"
        np.savez_compressed(out_series, **series_npz)
        np.savez_compressed(out_espectros, **espec_npz)
        print(f"[OK] Guardados ({keyW}):\n  - {out_series}\n  - {out_espectros}")

# ──────────────────────────
# Envoltura por escenario
# ──────────────────────────
def make_stations_from_grid(gr_bay, j_centro: int = _J_CENTRO_DEFAULT, radio: int = _RADIO_DEFAULT):
    """
    Construye la lista de estaciones a partir del grid.

    Busca una variable 2D en el grid para obtener (ny, nx). 
    Intenta primero 'XC', luego 'xC', y si no, toma la primera var 2D que encuentre.
    """
    # --- DEBUG opcional: ver qué variables tiene este archivo ---
    # print("DEBUG make_stations_from_grid vars:", list(gr_bay.variables.keys()))

    # Buscar una variable 2D para obtener dimensiones
    if "XC" in gr_bay.variables:
        base = gr_bay.variables["XC"]
    elif "xC" in gr_bay.variables:
        base = gr_bay.variables["xC"]
    else:
        # fallback: primera variable 2D que aparezca
        base = None
        for vname, v in gr_bay.variables.items():
            if hasattr(v, "shape") and len(v.shape) == 2:
                base = v
                # print(f"DEBUG usando {vname} como base para (ny,nx)")
                break
        if base is None:
            raise KeyError(
                "No encontré 'XC', 'xC' ni ninguna variable 2D para inferir (ny,nx)."
            )

    ny, nx = base.shape

    i_centro = nx // 2
    clamp = lambda v, lo, hi: max(lo, min(hi, v))

    j_centro = clamp(j_centro, 0, ny - 1)

    iS1, jS1 = i_centro, clamp(j_centro + 30, 0, ny - 1)
    iS2, jS2 = i_centro, clamp(j_centro + 10, 0, ny - 1)
    iA, jA  = clamp(i_centro - radio, 0, nx - 1), j_centro
    iB, jB  = clamp(i_centro + radio, 0, nx - 1), j_centro
    iC, jC  = i_centro, clamp(j_centro - radio, 0, ny - 1)
    iG, jG  = i_centro, clamp(j_centro - 2 * radio, 0, ny - 1)
    iO, jO  = i_centro, j_centro

    return [
        Station(iS1, jS1, "Seiche",  "#9d02d7"),
        Station(iS2, jS2, "Seiche2", "#ffd700"),
        Station(iA,  jA,  "A",       "#ea5f94"),
        Station(iB,  jB,  "B",       "#0000ff"),
        Station(iC,  jC,  "C",       "#fa8775"),
        Station(iG,  jG,  "Gamma",   "limegreen"),
        Station(iO,  jO,  "Omega",   "cyan"),
    ]


def guardar_W_series_y_espectros_escenario(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    tmin_horas: float = 27.0,
    save_dir: Path = Path("./npz_out"),
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    depths_m=(110.0, 250.0),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern   = "grid.t*.nc",
    prefer_RF: bool = False
):
    data = cargar_datos_mitgcm_etapas(outdir_bay, outdir_nobay,
                                      stage_patterns=stage_patterns,
                                      grid_pattern=grid_pattern)
    stations = make_stations_from_grid(data["gr_bay"], j_centro=j_centro, radio=radio)
    extra_meta = {
        "scenario": scenario_label,
        "paths": {"bay": str(outdir_bay), "nobay": str(outdir_nobay)},
        "stage_patterns": list(stage_patterns)
    }
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    try:
        guardar_W_series_y_espectros(
            nc_bay_list  = data["nc_bay_list"],
            nc_nobay_list= data["nc_nobay_list"],
            gr_bay       = data["gr_bay"],
            time_concat  = data["time"],
            stations     = stations,
            depths_m     = depths_m,
            tmin_horas   = tmin_horas,
            save_path    = save_dir,
            save_prefix_base = f"{scenario_label}_W",
            extra_meta   = extra_meta,
            prefer_RF    = prefer_RF
        )
    finally:
        _close_mnc(data["gr_bay"])
        for nc in data["nc_bay_list"]:
            _close_mnc(nc)
        for nc in data["nc_nobay_list"]:
            _close_mnc(nc)

def procesar_varios_escenarios_W(
    scenarios_dict: dict,
    subdir_bay: str = "run_expand",
    subdir_nobay: str = "run_expand_nobay",
    tmin_horas: float = 27.0,
    save_dir: Path = Path("./npz_out"),
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    depths_m=(110.0, 250.0),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    prefer_RF: bool = False,
) -> dict:
    """
    Procesa ΔW para todos los escenarios en ``scenarios_dict``.

    Retorna
    -------
    dict[str, list[Path]]
        Claves = etiquetas de escenario; valores = lista de rutas NPZ generadas.
    """
    save_dir = Path(save_dir)
    results: dict = {}
    for scen_label, base_path in scenarios_dict.items():
        outdir_bay   = Path(base_path) / subdir_bay
        outdir_nobay = Path(base_path) / subdir_nobay
        for d, name in [(outdir_bay, "bay"), (outdir_nobay, "nobay")]:
            if not d.exists():
                raise FileNotFoundError(
                    f"[{scen_label}] Directorio {name} no encontrado: {d}"
                )
        print(f"\n=== Escenario (W): {scen_label} ===")
        out = guardar_W_series_y_espectros_escenario(
            scenario_label = scen_label,
            outdir_bay     = outdir_bay,
            outdir_nobay   = outdir_nobay,
            tmin_horas     = tmin_horas,
            save_dir       = save_dir,
            j_centro       = j_centro,
            radio          = radio,
            depths_m       = depths_m,
            stage_patterns = stage_patterns,
            prefer_RF      = prefer_RF,
        )
        results[scen_label] = out if isinstance(out, list) else []
    print("\n[FIN] W para todos los escenarios.")
    return results

# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────
def _load_npz_with_meta(series_npz_path, espectros_npz_path):
    dseries = np.load(series_npz_path, allow_pickle=True)
    despec  = np.load(espectros_npz_path, allow_pickle=True)
    meta = json.loads(str(dseries["meta_json"]))
    return dseries, despec, meta

def _keyW_from_meta(meta):
    d = int(round(meta.get("depth_m", 0)))
    return f"W{d}m"

def _auto_scenario_colors(escenarios_npz: dict, scenario_colors: dict | None):
    keys = list(escenarios_npz.keys())
    if scenario_colors:
        return {k: scenario_colors.get(k) for k in keys}
    cmap = plt.get_cmap("tab10")
    return {k: cmap(i % 10) for i, k in enumerate(keys)}

def _detectar_estaciones(d0, keyW: str) -> list:
    """Infiere etiquetas de estación a partir de las claves del NPZ."""
    labels: set = set()
    pat1 = re.compile(r"^(?P<label>.+?)_t_horas$")
    pat2 = re.compile(rf"^(?P<label>.+?)_{re.escape(keyW)}$")
    for k in d0.files:
        m = pat1.match(k) or pat2.match(k)
        if m:
            labels.add(m.group("label"))
    ordered = sorted(labels)
    if "Seiche" in ordered:
        ordered.insert(0, ordered.pop(ordered.index("Seiche")))
    return ordered

def _build_color_map(stations_meta: list, estaciones: list) -> dict:
    """Construye color por estación desde meta; rellena con tab10 si falta."""
    color_map: dict = {}
    for st in stations_meta:
        lab = st.get("label")
        col = st.get("color")
        if lab and col not in (None, "", "None"):
            color_map[lab] = col
    cmap = plt.get_cmap("tab10")
    for i, lab in enumerate(sorted(estaciones)):
        if color_map.get(lab) in (None, "", "None", {}):
            color_map[lab] = cmap(i % 10)
    return color_map

def _verificar_consistencia_keyW(escenarios_npz: dict, keyW: str) -> None:
    """Lanza ValueError si algún NPZ tiene profundidad distinta a keyW."""
    for scen, paths in escenarios_npz.items():
        dm = np.load(paths["series"], allow_pickle=True)
        meta = json.loads(str(dm["meta_json"]))
        kw = _keyW_from_meta(meta)
        if kw != keyW:
            raise ValueError(
                f"Profundidad inconsistente en '{scen}': "
                f"{paths['series']} tiene '{kw}', esperado '{keyW}'."
            )

def _cargar_datos_W(
    escenarios_npz: dict,
    estaciones: list,
    keyW: str,
) -> tuple:
    """
    Carga series y espectros de todos los escenarios en una sola pasada.

    Retorna
    -------
    data_series : list[(estacion, scen_name, t_h, dW)]
    data_spec   : list[(estacion, scen_name, f_cph, PSD)]
    scenarios_order : list[str]   — orden de iteración
    """
    data_series: list = []
    data_spec:   list = []
    scenarios_order: list = []

    for scen, paths in escenarios_npz.items():
        dser = np.load(paths["series"],    allow_pickle=True)
        desp = np.load(paths["espectros"], allow_pickle=True)
        meta = json.loads(str(dser["meta_json"]))
        scen_name = meta.get("scenario", scen)
        scenarios_order.append(scen_name)

        for lab in estaciones:
            try:
                t_h = dser[f"{lab}_t_horas"]
                dW  = dser[f"{lab}_{keyW}"]
            except KeyError:
                raise KeyError(
                    f"{paths['series']}: sin claves para estación '{lab}' ({keyW})."
                )
            try:
                f_cph = desp[f"{lab}_freq_cph"]
                PSD   = desp[f"{lab}_PSD_{keyW}"]
            except KeyError:
                raise KeyError(
                    f"{paths['espectros']}: sin claves espectrales para '{lab}' ({keyW})."
                )
            data_series.append((lab, scen_name, t_h, dW))
            data_spec.append((lab, scen_name, f_cph, PSD))

    return data_series, data_spec, scenarios_order

def _make_legend_handles(
    scenarios_order: list,
    estaciones: list,
    color_map: dict,
    scenario_linestyles: dict,
    linewidth: float,
) -> tuple:
    """Construye los handles de leyenda para escenarios y estaciones."""
    scen_handles = [
        Line2D([0], [0], color="k",
               linestyle=scenario_linestyles.get(s, "-"),
               linewidth=linewidth, label=s)
        for s in scenarios_order
    ]
    stn_handles = [
        Line2D([0], [0], color=color_map.get(lab, "k"),
               linestyle="-", linewidth=linewidth, label=lab)
        for lab in estaciones
    ]
    return scen_handles, stn_handles

# ─────────────────────────────────────────────────────────────
# 1) Estaciones en un escenario (W a una profundidad)
# ─────────────────────────────────────────────────────────────
def graficar_W_estaciones_en_escenario(
    series_npz_path: str,
    espectros_npz_path: str,
    estaciones_orden: list | None = None,
    seiche_freq_cph: float = 0.303356,
    coriolis_f: float = 0.04,
    xlim_cph: tuple = (0, 0.8),
    ylim_psd: tuple = (1e-13, 1e-7),
    tlim: tuple | None = None,
    which: tuple = ("series", "espectro"),
    titulo_suffix: str | None = None,
    figsize: tuple = (7, 5),
    linewidth: float = 1.4,
    alpha: float = 1.0,
    save_fig: bool = False,
    show_fig: bool = True,
    save_dir: str | Path = "./figs",
    dpi: int = 150,
) -> dict:
    """
    Grafica ΔW (bahía − sin bahía) para todas las estaciones de un escenario.

    Parámetros
    ----------
    which       : subconjunto de figuras a generar: ``("series",)``,
                  ``("espectro",)`` o ambas ``("series", "espectro")``.
    tlim        : (t_min, t_max) en horas para el eje temporal; None = auto.
    save_fig    : si True guarda las figuras en ``save_dir``.
    show_fig    : si True llama plt.show().

    Retorna
    -------
    dict con claves ``"series"`` y/o ``"espectro"`` apuntando a los objetos Figure.
    """
    dseries, despec, meta = _load_npz_with_meta(series_npz_path, espectros_npz_path)
    scenario = meta.get("scenario", "escenario")
    keyW  = _keyW_from_meta(meta)
    depth = meta.get("depth_m", None)

    stations_meta = meta.get("stations", [])
    if estaciones_orden is None:
        estaciones_orden = _detectar_estaciones(dseries, keyW)
    color_map = _build_color_map(stations_meta, estaciones_orden)

    depth_tag = f"{int(round(depth))}m" if depth is not None else "Xm"
    save_dir  = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    out: dict = {}

    # ── Serie temporal ───────────────────────────────────────────────────
    if "series" in which:
        fig, ax = plt.subplots(figsize=figsize)
        for lab in estaciones_orden:
            t_h = dseries[f"{lab}_t_horas"]
            dW  = dseries[f"{lab}_{keyW}"]
            ax.plot(t_h, dW, label=f"Est. {lab}",
                    color=color_map.get(lab), linewidth=linewidth, alpha=alpha)
        ax.set_xlabel("Tiempo (horas)")
        ax.set_ylabel(f"Δ{keyW} (m s$^{{-1}}$)")
        title = f"ΔW (bahía − sin bahía) • {scenario}"
        if depth is not None:
            title += f" • {depth_tag}"
        if titulo_suffix:
            title += f"\n{titulo_suffix}"
        ax.set_title(title)
        if tlim is not None:
            ax.set_xlim(*tlim)
        ax.grid(True)
        ax.legend()
        fig.tight_layout()
        if save_fig:
            p = save_dir / f"{scenario}_W_{depth_tag}_series_estaciones.png"
            fig.savefig(p, dpi=dpi)
            print(f"[OK] {p}")
        if show_fig:
            plt.show()
        out["series"] = fig

    # ── Espectro ─────────────────────────────────────────────────────────
    if "espectro" in which:
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_yscale("log", base=10)
        for lab in estaciones_orden:
            f_cph = despec[f"{lab}_freq_cph"]
            PSD   = despec[f"{lab}_PSD_{keyW}"]
            ax.plot(f_cph, PSD, label=f"Est. {lab}",
                    color=color_map.get(lab), linewidth=linewidth, alpha=alpha)
        ax.axvline(seiche_freq_cph, color="k",    ls="--", lw=0.8, label="f seiche")
        ax.axvline(coriolis_f,      color="gray",  ls="--", lw=0.8, label="f Coriolis")
        ax.set_xlim(*xlim_cph)
        ax.set_ylim(*ylim_psd)
        ax.set_xlabel("Frecuencia (ciclos/hora)")
        ax.set_ylabel(f"PSD($\\Delta${keyW}) [m$^2$ s$^{-2}$ cph$^{-1}$]")
        ax.set_title(f"Espectro Δ{keyW} • {scenario}" +
                     (f" • {depth_tag}" if depth is not None else ""))
        ax.grid(True, which="both", ls=":")
        ax.legend()
        fig.tight_layout()
        if save_fig:
            p = save_dir / f"{scenario}_W_{depth_tag}_espectro_estaciones.png"
            fig.savefig(p, dpi=dpi)
            print(f"[OK] {p}")
        if show_fig:
            plt.show()
        out["espectro"] = fig

    return out

# ─────────────────────────────────────────────────────────────
# 2) Comparar escenarios en una estación (W a una profundidad)
# ─────────────────────────────────────────────────────────────
def comparar_escenarios_W_en_estacion(
    escenarios_npz: dict,
    # dict: {
    #   "linear":    {"series": ".../linear_W_W110m_series.npz", "espectros":".../linear_W_W110m_espectros.npz"},
    #   "realistic": {"series": "...",                            "espectros":"..."},
    #   ...
    # }
    estacion_label: str = "A",
    seiche_freq_cph: float = 0.303356,
    coriolis_f: float = 0.04,   # cph for lat ~ 41.5°N
    xlim_cph=(0, 0.8),
    ylim_psd=(1e-13, 1e-7),
    scenario_colors: dict | None = None,   # {"linear":"tab:orange", "realistic":"tab:green", ...}
    titulo_suffix: str | None = None
):
    """
    Compara ΔW (bahía − sin bahía) a una profundidad (implícita en los archivos)
    entre varios escenarios, para UNA estación.
    """
    scen_color = _auto_scenario_colors(escenarios_npz, scenario_colors)

    # Detecta la profundidad esperada a partir del primer archivo
    first_key = next(iter(escenarios_npz))
    d0 = np.load(escenarios_npz[first_key]["series"], allow_pickle=False)
    meta0 = json.loads(str(d0["meta_json"]))
    keyW = _keyW_from_meta(meta0)
    depth = meta0.get("depth_m", None)

    series_data = []
    espec_data  = []

    for scen, paths in escenarios_npz.items():
        dser = np.load(paths["series"], allow_pickle=False)
        desp = np.load(paths["espectros"], allow_pickle=False)
        meta = json.loads(str(dser["meta_json"]))
        scen_name = meta.get("scenario", scen)

        # Verifica que todos correspondan a la misma profundidad
        keyW_here = _keyW_from_meta(meta)
        if keyW_here != keyW:
            raise ValueError(f"Profundidad inconsistente: {paths['series']} es {keyW_here}, esperado {keyW}.")

        # Extrae la estación
        try:
            t_h = dser[f"{estacion_label}_t_horas"]
            dW  = dser[f"{estacion_label}_{keyW}"]
        except KeyError:
            raise KeyError(f"{paths['series']} no tiene claves para la estación '{estacion_label}' ({keyW}).")

        try:
            f_cph = desp[f"{estacion_label}_freq_cph"]
            PSD   = desp[f"{estacion_label}_PSD_{keyW}"]
        except KeyError:
            raise KeyError(f"{paths['espectros']} no tiene claves espectrales para '{estacion_label}' ({keyW}).")

        series_data.append((scen_name, t_h, dW))
        espec_data.append((scen_name, f_cph, PSD))

    # Serie temporal
    plt.figure(figsize=(7,5))
    for scen_name, t_h, dW in series_data:
        plt.plot(t_h, dW, label=scen_name, color=scen_color.get(scen_name))

    plt.xlabel("Tiempo (horas)")
    plt.ylabel(f"Δ {keyW} (m s$^{{-1}}$)")
    title = f"Comparación de escenarios • Estación {estacion_label}"
    if depth is not None:
        title += f" • {int(round(depth))} m"
    if titulo_suffix:
        title += f"\n{titulo_suffix}"
    plt.title(title)
    plt.grid(True)
    plt.legend(title="Escenario")
    plt.tight_layout()
    plt.show()

    # Espectro
    plt.figure(figsize=(7,5))
    plt.yscale("log", base=10)
    for scen_name, f_cph, PSD in espec_data:
        plt.plot(f_cph, PSD, label=scen_name, color=scen_color.get(scen_name))

    plt.axvline(seiche_freq_cph, color='k', linestyle='--', linewidth=0.7, label='f seiche')
    plt.axvline(coriolis_f, color='gray', linestyle='--', linewidth=0.7, label='f Coriolis')
    plt.xlim(*xlim_cph)
    plt.ylim(*ylim_psd)
    plt.xlabel("Frecuencia (ciclos/hora)")
    plt.ylabel(f"PSD($\\Delta${keyW}) [m$^2$ s$^{{-2}}$ cph$^{{-1}}$]")
    plt.title(f"Comparación de espectros • Estación {estacion_label}")
    plt.grid(True, which="both", ls=':')
    plt.legend(title="Escenario")
    plt.tight_layout()
    plt.show()

# ─────────────────────────────────────────────────────────────
# 2b) Comparar escenarios en TODAS las estaciones (W a una profundidad)
#     • Color = estación  •  Estilo de línea = escenario
# ─────────────────────────────────────────────────────────────
def comparar_escenarios_W_en_todas_estaciones(
    escenarios_npz: dict,
    estaciones_orden: list | None = None,
    seiche_freq_cph: float = 0.303356,
    coriolis_f: float = 0.04,
    xlim_cph: tuple = (0, 0.8),
    ylim_psd: tuple = (1e-13, 1e-7),
    tlim: tuple | None = None,
    which: tuple = ("series", "espectro"),
    scenario_linestyles: dict | None = None,
    titulo_suffix: str | None = None,
    figsize: tuple = (7, 5),
    figsize_panel: tuple = (7, 3),
    linewidth: float = 1.4,
    alpha: float = 1.0,
    panel_layout: bool = False,
    save_fig: bool = False,
    show_fig: bool = True,
    save_dir: str | Path = "./figs",
    dpi: int = 150,
) -> dict:
    """
    Compara ΔW entre escenarios para TODAS las estaciones.

    Convención visual: color = estación, estilo de línea = escenario.

    Modos
    -----
    ``panel_layout=False`` (default)
        Una figura separada por estación para cada tipo (serie / espectro).
        Retorna ``{estacion: {"series": Fig, "espectro": Fig}}``.

    ``panel_layout=True``
        Una única figura con N filas (una por estación) y 2 columnas:
        columna izquierda = series temporales, columna derecha = espectros.
        Retorna ``{"panel": Fig}``.
    """
    _ls = scenario_linestyles or {
        "linear": "-", "2layer": "--", "realistic": ":", "barotropic": "-.",
    }

    # ── Cargar referencia y validar ──────────────────────────────────────
    first_key = next(iter(escenarios_npz))
    d0    = np.load(escenarios_npz[first_key]["series"], allow_pickle=True)
    meta0 = json.loads(str(d0["meta_json"]))
    keyW  = _keyW_from_meta(meta0)
    depth = meta0.get("depth_m", None)
    _verificar_consistencia_keyW(escenarios_npz, keyW)

    estaciones = estaciones_orden if estaciones_orden else _detectar_estaciones(d0, keyW)
    color_map  = _build_color_map(meta0.get("stations", []), estaciones)

    # ── Cargar todos los datos en una sola pasada ────────────────────────
    data_series, data_spec, scenarios_order = _cargar_datos_W(
        escenarios_npz, estaciones, keyW
    )

    depth_tag = f"{int(round(depth))}m" if depth is not None else "Xm"
    psd_ylabel = f"PSD($\\Delta${keyW}) [m$^2$ s$^{{-2}}$ cph$^{{-1}}$]"
    save_dir   = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    out: dict  = {}

    # ════════════════════════════════════════════════════════════════════
    # MODO PANEL: una figura, N filas × 2 columnas
    # ════════════════════════════════════════════════════════════════════
    if panel_layout:
        nrows  = len(estaciones)
        fw     = figsize_panel[0] * 2 + 0.5
        fh     = figsize_panel[1] * nrows
        fig, axes = plt.subplots(
            nrows=nrows, ncols=2,
            figsize=(fw, fh),
            gridspec_kw={"wspace": 0.35, "hspace": 0.18},
        )
        if nrows == 1:
            axes = axes[np.newaxis, :]   # asegurar 2D

        for row, estacion_label in enumerate(estaciones):
            col  = color_map.get(estacion_label)
            s_data = [(sn, t, dw) for lab, sn, t, dw in data_series if lab == estacion_label]
            e_data = [(sn, f, p)  for lab, sn, f, p  in data_spec   if lab == estacion_label]

            ax_s = axes[row, 0]
            ax_e = axes[row, 1]

            # ── Serie temporal ────────────────────────────────────────
            for scen_name, t_h, dW in s_data:
                ax_s.plot(t_h, dW, label=scen_name,
                          linestyle=_ls.get(scen_name, "-"),
                          color=col, linewidth=linewidth, alpha=alpha)
            if tlim is not None:
                ax_s.set_xlim(*tlim)
            ax_s.grid(True, ls=":")
            ax_s.set_ylabel(
                f"Est. {estacion_label}\n"
                f"$\\Delta${keyW} (m s$^{{-1}}$)",
                color=col, fontsize=9,
            )
            ax_s.tick_params(axis="y", labelcolor=col)
            if row < nrows - 1:
                ax_s.tick_params(labelbottom=False)
            else:
                ax_s.set_xlabel("Tiempo (horas)")

            # ── Espectro ──────────────────────────────────────────────
            ax_e.set_yscale("log", base=10)
            for scen_name, f_cph, PSD in e_data:
                ax_e.plot(f_cph, PSD, label=scen_name,
                          linestyle=_ls.get(scen_name, "-"),
                          color=col, linewidth=linewidth, alpha=alpha)
            ax_e.axvline(seiche_freq_cph, color="k",   ls="--", lw=0.7)
            ax_e.axvline(coriolis_f,      color="gray", ls="--", lw=0.7)
            ax_e.set_xlim(*xlim_cph)
            ax_e.set_ylim(*ylim_psd)
            ax_e.grid(True, which="both", ls=":")
            ax_e.set_ylabel(psd_ylabel, fontsize=8)
            if row < nrows - 1:
                ax_e.tick_params(labelbottom=False)
            else:
                ax_e.set_xlabel("Frecuencia (ciclos/hora)")

        # ── Títulos de columna ────────────────────────────────────────
        axes[0, 0].set_title(
            f"ΔW(t) — comparación de escenarios • {depth_tag}" +
            (f"\n{titulo_suffix}" if titulo_suffix else ""),
            fontsize=10,
        )
        axes[0, 1].set_title(
            f"PSD(ΔW) — comparación de escenarios • {depth_tag}" +
            (f"\n{titulo_suffix}" if titulo_suffix else ""),
            fontsize=10,
        )

        # ── Leyenda de escenarios (una por columna) ───────────────────
        scen_handles, _ = _make_legend_handles(
            scenarios_order, estaciones, color_map, _ls, linewidth
        )
        for ax_leg in (axes[0, 0], axes[0, 1]):
            ax_leg.legend(handles=scen_handles, title="Escenario",
                          fontsize=8, title_fontsize=8,
                          loc="upper right", framealpha=0.8)

        fig.tight_layout()
        if save_fig:
            p = save_dir / f"comp_W_{depth_tag}_panel_layout.png"
            fig.savefig(p, dpi=dpi)
            print(f"[OK] {p}")
        if show_fig:
            plt.show()
        out["panel"] = fig
        return out

    # ════════════════════════════════════════════════════════════════════
    # MODO CLÁSICO: figura separada por estación
    # ════════════════════════════════════════════════════════════════════
    for estacion_label in estaciones:
        s_data = [(sn, t, dw) for lab, sn, t, dw  in data_series if lab == estacion_label]
        e_data = [(sn, f, p)  for lab, sn, f, p   in data_spec   if lab == estacion_label]
        out[estacion_label] = {}

        if "series" in which:
            fig, ax = plt.subplots(figsize=figsize)
            for scen_name, t_h, dW in s_data:
                ax.plot(t_h, dW, label=scen_name,
                        linestyle=_ls.get(scen_name, "-"),
                        color=color_map.get(estacion_label),
                        linewidth=linewidth, alpha=alpha)
            ax.set_xlabel("Tiempo (horas)")
            ax.set_ylabel(f"$\\Delta${keyW} (m s$^{{-1}}$)")
            title = f"Comparación de escenarios • Est. {estacion_label}"
            if depth is not None:
                title += f" • {depth_tag}"
            if titulo_suffix:
                title += f"\n{titulo_suffix}"
            ax.set_title(title)
            if tlim is not None:
                ax.set_xlim(*tlim)
            ax.grid(True)
            ax.legend(title="Escenario")
            fig.tight_layout()
            if save_fig:
                p = save_dir / f"comp_W_{depth_tag}_series_est{estacion_label}.png"
                fig.savefig(p, dpi=dpi)
                print(f"[OK] {p}")
            if show_fig:
                plt.show()
            out[estacion_label]["series"] = fig

        if "espectro" in which:
            fig, ax = plt.subplots(figsize=figsize)
            ax.set_yscale("log", base=10)
            for scen_name, f_cph, PSD in e_data:
                ax.plot(f_cph, PSD, label=scen_name,
                        linestyle=_ls.get(scen_name, "-"),
                        color=color_map.get(estacion_label),
                        linewidth=linewidth, alpha=alpha)
            ax.axvline(seiche_freq_cph, color="k",   ls="--", lw=0.8, label="f seiche")
            ax.axvline(coriolis_f,      color="gray", ls="--", lw=0.8, label="f Coriolis")
            ax.set_xlim(*xlim_cph)
            ax.set_ylim(*ylim_psd)
            ax.set_xlabel("Frecuencia (ciclos/hora)")
            ax.set_ylabel(psd_ylabel)
            ax.set_title(f"Comparación de espectros • Est. {estacion_label}" +
                         (f" • {depth_tag}" if depth is not None else ""))
            ax.grid(True, which="both", ls=":")
            ax.legend(title="Escenario")
            fig.tight_layout()
            if save_fig:
                p = save_dir / f"comp_W_{depth_tag}_espectro_est{estacion_label}.png"
                fig.savefig(p, dpi=dpi)
                print(f"[OK] {p}")
            if show_fig:
                plt.show()
            out[estacion_label]["espectro"] = fig

    return out


# ─────────────────────────────────────────────────────────────
# 2c) Escenarios combinados en una sola figura por tipo (todas las estaciones)
#     — Colores = estaciones;  estilo de línea = escenarios
# ─────────────────────────────────────────────────────────────
def comparar_escenarios_W_en_todas_estaciones_combinado(
    escenarios_npz: dict,
    estaciones_orden: list | None = None,
    seiche_freq_cph: float = 0.303356,
    coriolis_f: float = 0.04,
    xlim_cph: tuple = (0, 0.8),
    ylim_psd: tuple = (1e-13, 1e-7),
    tlim: tuple | None = None,
    which: tuple = ("series", "espectro"),
    scenario_linestyles: dict | None = None,
    titulo_suffix: str | None = None,
    alpha: float = 0.95,
    linewidth: float = 1.6,
    figsize: tuple = (9, 6),
    save_fig: bool = False,
    show_fig: bool = True,
    save_dir: str | Path = "./figs",
    dpi: int = 150,
) -> dict:
    """
    Genera hasta dos figuras combinadas (una por tipo):

    * ``"series"``  — ΔW(t) de todas las estaciones y escenarios.
    * ``"espectro"`` — PSD(f) ídem.

    Convención visual: color = estación, estilo de línea = escenario.
    Dos leyendas por figura: una para escenarios (estilos) y otra para
    estaciones (colores).

    Retorna
    -------
    dict con claves ``"series"`` y/o ``"espectro"`` apuntando a los objetos Figure.
    """
    _ls = scenario_linestyles or {
        "linear": "-", "2layer": "--", "realistic": ":", "barotropic": "-.",
    }

    # ── Referencia y validación ──────────────────────────────────────────
    first_key = next(iter(escenarios_npz))
    d0    = np.load(escenarios_npz[first_key]["series"], allow_pickle=True)
    meta0 = json.loads(str(d0["meta_json"]))
    keyW  = _keyW_from_meta(meta0)
    depth = meta0.get("depth_m", None)
    _verificar_consistencia_keyW(escenarios_npz, keyW)

    estaciones = estaciones_orden if estaciones_orden else _detectar_estaciones(d0, keyW)
    color_map  = _build_color_map(meta0.get("stations", []), estaciones)

    # ── Carga única de datos ─────────────────────────────────────────────
    data_series, data_spec, scenarios_order = _cargar_datos_W(
        escenarios_npz, estaciones, keyW
    )

    depth_tag = f"{int(round(depth))}m" if depth is not None else "Xm"
    save_dir  = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    out: dict = {}

    # ── Helper interno para añadir las dos leyendas ──────────────────────
    def _add_legends(ax):
        sh, sth = _make_legend_handles(
            scenarios_order, estaciones, color_map, _ls, linewidth
        )
        leg_sc  = ax.legend(handles=sh,  title="Escenarios", loc="upper left",  frameon=True)
        ax.legend(handles=sth, title="Estaciones", loc="upper right", frameon=True)
        ax.add_artist(leg_sc)  # conservar ambas

    # ============== FIGURA 1: Series temporales ==========================
    if "series" in which:
        fig1, ax1 = plt.subplots(figsize=figsize)
        for lab, scen_name, t_h, dW in data_series:
            ax1.plot(t_h, dW,
                     color=color_map.get(lab),
                     linestyle=_ls.get(scen_name, "-"),
                     linewidth=linewidth, alpha=alpha)
        ax1.set_xlabel("Tiempo (horas)")
        ax1.set_ylabel(f"Δ{keyW} (m s$^{{-1}}$)")
        title1 = f"Δ{keyW}(t) — Todas las estaciones"
        if depth is not None:
            title1 += f" • {depth_tag}"
        if titulo_suffix:
            title1 += f"\n{titulo_suffix}"
        ax1.set_title(title1)
        if tlim is not None:
            ax1.set_xlim(*tlim)
        ax1.grid(True, which="both", ls=":")
        _add_legends(ax1)
        fig1.tight_layout()
        if save_fig:
            p = save_dir / f"comp_W_{depth_tag}_series_combinado.png"
            fig1.savefig(p, dpi=dpi)
            print(f"[OK] {p}")
        if show_fig:
            plt.show()
        out["series"] = fig1

    # ============== FIGURA 2: Espectros ==================================
    if "espectro" in which:
        fig2, ax2 = plt.subplots(figsize=figsize)
        ax2.set_yscale("log", base=10)
        for lab, scen_name, f_cph, PSD in data_spec:
            ax2.plot(f_cph, PSD,
                     color=color_map.get(lab),
                     linestyle=_ls.get(scen_name, "-"),
                     linewidth=linewidth, alpha=alpha)
        ax2.axvline(seiche_freq_cph, color="k",   ls="--", lw=0.9, label="f seiche")
        ax2.axvline(coriolis_f,      color="gray", ls="--", lw=0.9, label="f Coriolis")
        ax2.set_xlim(*xlim_cph)
        ax2.set_ylim(*ylim_psd)
        ax2.set_xlabel("Frecuencia (ciclos/hora)")
        ax2.set_ylabel(f"PSD($\\Delta${keyW}) [m$^2$ s$^{{-2}}$ cph$^{{-1}}$]")
        title2 = f"PSD(Δ{keyW}) — Todas las estaciones"
        if depth is not None:
            title2 += f" • {depth_tag}"
        if titulo_suffix:
            title2 += f"\n{titulo_suffix}"
        ax2.set_title(title2)
        ax2.grid(True, which="both", ls=":")
        _add_legends(ax2)
        fig2.tight_layout()
        if save_fig:
            p = save_dir / f"comp_W_{depth_tag}_espectro_combinado.png"
            fig2.savefig(p, dpi=dpi)
            print(f"[OK] {p}")
        if show_fig:
            plt.show()
        out["espectro"] = fig2

    return out


# Asumo que ya tienes estas funciones en tu módulo:
# - cargar_datos_mitgcm_etapas(outdir_bay, outdir_nobay, stage_patterns, grid_pattern)
# - make_stations_from_grid(gr_bay, j_centro=315, radio=50)
# - _get_time_from_nc(nc), _concat_time_and_data(times_list, data_list, tol)
# Si alguna te falta, dímelo y te la pego completa aquí.

def guardar_hovmoller_npz_escenario_unico(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    varname: str = 'W',                 # variable 3D (nt, nz, y, x)
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern   = "grid.t*.nc"
):
    """
    Construye y guarda EN UN SOLO NPZ el Hovmöller (tiempo × profundidad × estación)
    de (bay − nobay) para todas las estaciones del escenario.
    Salida: <save_dir>/<scenario_label>_hov_ALL_<varname>.npz
    Contiene:
      - t_horas: (nt,)
      - RC: (nz,)
      - diff: (nt, nz, nstations)
      - labels: (nstations,)   (dtype=object)
      - ij: (nstations, 2)     int32 con (i, j) por estación
      - colors: (nstations,)   (dtype=object, puede contener None)
      - idx_100, idx_400: int
      - meta_json: json con info de escenario y parámetros
    """
    save_dir.mkdir(parents=True, exist_ok=True)

    # A) Carga multi-etapas (usa tus helpers existentes)
    data = cargar_datos_mitgcm_etapas(
        outdir_bay, outdir_nobay,
        stage_patterns=stage_patterns,
        grid_pattern=grid_pattern
    )
    nc_bay_list   = data["nc_bay_list"]
    nc_nobay_list = data["nc_nobay_list"]
    gr_bay        = data["gr_bay"]
    time_concat   = data["time"]  # segundos

    # B) Estaciones y eje vertical
    stations = make_stations_from_grid(gr_bay, j_centro=j_centro, radio=radio)  # [(i,j,label,color), ...]
    RC = _get_RC(gr_bay)
    nz = RC.shape[0]

    # C) Recorte temporal
    t_all = np.array(time_concat, dtype=float)
    tmin_seg = tmin_horas * 3600.0
    idx0 = np.where(t_all >= tmin_seg)[0][0] if np.any(t_all >= tmin_seg) else 0
    t_cut = t_all[idx0:]
    t_h   = t_cut / 3600.0
    nt = t_h.shape[0]

    # D) Índices guía (≈100 y ≈400 m)
    idx_100 = int(np.argmin(np.abs(RC + 100.0)))
    idx_400 = int(np.argmin(np.abs(RC + 400.0)))

    # E) Construcción del cubo diff (nt × nz × nstations)
    nst = len(stations)
    diff_all = np.empty((nt, nz, nst), dtype=np.float32)
    labels = []
    colors = []
    ij     = np.empty((nst, 2), dtype=np.int32)

    for k, st in enumerate(stations):
        i, j, label, color = st

        # Recolecta por etapas (tiempos y bloques nt×nz)
        times_b, blocks_b = [], []
        times_n, blocks_n = [], []
        for nc_b, nc_n in zip(nc_bay_list, nc_nobay_list):
            tb = _get_time_from_nc(nc_b)
            tn = _get_time_from_nc(nc_n)
            Xb = np.array(nc_b.variables[varname][:, :, j, i], dtype=float)  # (nt, nz)
            Xn = np.array(nc_n.variables[varname][:, :, j, i], dtype=float)
            if Xb.shape[0] != tb.shape[0] or Xn.shape[0] != tn.shape[0]:
                raise ValueError("Longitudes tiempo/serie no coinciden en alguna etapa.")
            times_b.append(tb); blocks_b.append(Xb)
            times_n.append(tn); blocks_n.append(Xn)

        t_b, Xb_all = _concat_time_and_data(times_b, blocks_b)
        t_n, Xn_all = _concat_time_and_data(times_n, blocks_n)

        # Diferencia y recorte temporal (usando idx0 de t_all)
        Xdiff = Xb_all - Xn_all        # (nt_full, nz)
        Xdiff_cut = Xdiff[idx0:, :]    # (nt, nz)
        if Xdiff_cut.shape[0] != nt or Xdiff_cut.shape[1] != nz:
            raise ValueError(f"Forma inesperada tras recorte en estación {label}: {Xdiff_cut.shape}, esperado {(nt, nz)}.")

        diff_all[:, :, k] = Xdiff_cut
        labels.append(label)
        colors.append(color)
        ij[k, 0] = int(i)
        ij[k, 1] = int(j)

    # F) Metadatos
    meta = {
        "scenario": scenario_label,
        "var": varname,
        "paths": {"bay": str(outdir_bay), "nobay": str(outdir_nobay)},
        "stage_patterns": list(stage_patterns),
        "tmin_horas": float(tmin_horas),
        "idx0": int(idx0),
        "dt_seconds": float(t_cut[1] - t_cut[0]) if len(t_cut) > 1 else None,
        "stations": [{"label": labels[k], "i": int(ij[k,0]), "j": int(ij[k,1]), "color": colors[k]} for k in range(nst)]
    }
    meta_json = json.dumps(meta, ensure_ascii=False)

    # G) Guardado
    out_path = save_dir / f"{scenario_label}_hov_ALL_{varname}.npz"
    np.savez_compressed(
        out_path,
        t_horas=t_h.astype(np.float64),
        RC=RC.astype(np.float64),
        diff=diff_all,                               # (nt, nz, nstations)
        labels=np.array(labels, dtype=object),
        colors=np.array(colors, dtype=object),
        ij=ij.astype(np.int32),
        idx_100=np.array(idx_100, dtype=np.int32),
        idx_400=np.array(idx_400, dtype=np.int32),
        meta_json=np.array(meta_json)
    )
    print(f"[OK] Bundle Hovmöller guardado: {out_path}")
    return out_path

def guardar_hovmoller_uv_npz_escenario_unico(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern   = "grid.t*.nc",
    source: str = "diff",               # "diff" (bay-nobay), "bay", "nobay"
    tol_time: float = 1e-9
):
    """
    Construye y guarda EN UN SOLO NPZ los Hovmöller (tiempo × profundidad × estación)
    de U y V interpoladas a la malla de centros (T,Z,Y,X):
      - U(T,Z,Y,Xp1) -> Uc(T,Z,Y,X) con 0.5*(U[..., :-1] + U[..., 1:])
      - V(T,Z,Yp1,X) -> Vc(T,Z,Y,X) con 0.5*(V[..., :-1, :] + V[..., 1:, :])

    Salida: <save_dir>/<scenario_label>_hov_ALL_UV_<source>.npz
    Contiene:
      - t_horas: (nt,)
      - RC: (nz,)
      - U_hov: (nt, nz, nstations)
      - V_hov: (nt, nz, nstations)
      - labels: (nstations,)   (dtype=object)
      - ij: (nstations, 2)     int32 con (i, j) por estación
      - idx_100, idx_400: int
      - meta_json: json con info de escenario y parámetros
    """
    save_dir.mkdir(parents=True, exist_ok=True)

    def _u_center_col(varU, j: int, i: int):
        """
        Extrae Uc(:, :, j, i) sin construir todo Uc 4D.
        U(T,Z,Y,Xp1) -> Uc(T,Z,Y,X) por promedio en x.
        """
        nxp1 = int(varU.shape[-1])
        if nxp1 <= 1:
            return np.array(varU[:, :, j, 0], dtype=float)
        if i > 0:
            return 0.5 * (
                np.array(varU[:, :, j, i], dtype=float) +
                np.array(varU[:, :, j, i - 1], dtype=float)
            )
        # borde izquierdo (i=0)
        return 0.5 * (
            np.array(varU[:, :, j, 0], dtype=float) +
            np.array(varU[:, :, j, 1], dtype=float)
        )

    def _v_center_col(varV, j: int, i: int):
        """
        Extrae Vc(:, :, j, i) sin construir todo Vc 4D.
        V(T,Z,Yp1,X) -> Vc(T,Z,Y,X) por promedio en y.
        """
        nyp1 = int(varV.shape[-2])
        if nyp1 <= 1:
            return np.array(varV[:, :, 0, i], dtype=float)
        if j > 0:
            return 0.5 * (
                np.array(varV[:, :, j, i], dtype=float) +
                np.array(varV[:, :, j - 1, i], dtype=float)
            )
        # borde sur (j=0)
        return 0.5 * (
            np.array(varV[:, :, 0, i], dtype=float) +
            np.array(varV[:, :, 1, i], dtype=float)
        )

    data = cargar_datos_mitgcm_etapas(
        outdir_bay, outdir_nobay,
        stage_patterns=stage_patterns,
        grid_pattern=grid_pattern
    )
    nc_bay_list   = data["nc_bay_list"]
    nc_nobay_list = data["nc_nobay_list"]
    gr_bay        = data["gr_bay"]
    try:
        stations = make_stations_from_grid(gr_bay, j_centro=j_centro, radio=radio)

        RC = _get_RC(gr_bay)
        nz = int(RC.shape[0])

        nst = len(stations)
        t_lists = [[] for _ in range(nst)]
        U_lists = [[] for _ in range(nst)]
        labels, colors = [], []
        ij = np.empty((nst, 2), dtype=np.int32)

        for k, st in enumerate(stations):
            i, j, label, color = st
            labels.append(label)
            colors.append(color)
            ij[k, 0] = int(i)
            ij[k, 1] = int(j)

        # ============================================================
        # Fase 1: construir Hovmöller de U (antes de tocar V)
        # ============================================================
        for nc_b, nc_n in zip(nc_bay_list, nc_nobay_list):
            try:
                tb = _get_time_from_nc(nc_b)
                tn = _get_time_from_nc(nc_n)
                Ub = nc_b.variables['U']
                Un = nc_n.variables['U']

                if Ub.shape[0] != tb.shape[0] or Un.shape[0] != tn.shape[0]:
                    raise ValueError("Longitudes tiempo/U no coinciden en una etapa.")
                if tb.shape != tn.shape or np.any(np.abs(tb - tn) > 1e-7):
                    raise ValueError("Los tiempos de bay y nobay difieren en una etapa.")

                # Evita duplicar el primer timestep cuando se concatenan etapas.
                trim0 = False
                if t_lists[0]:
                    trim0 = abs(tb[0] - t_lists[0][-1][-1]) <= tol_time
                sl = slice(1, None) if trim0 else slice(None)
                t_use = tb[sl]

                ny_u = int(Ub.shape[-2])
                nx_u = int(Ub.shape[-1]) - 1  # U en xp1 -> centros en x
                for sidx, (i, j) in enumerate(ij):
                    if not (0 <= j < ny_u and 0 <= i < nx_u):
                        raise IndexError(
                            f"Estación fuera de rango para U centrada: "
                            f"(i,j)=({i},{j}) vs (nx,ny)=({nx_u},{ny_u})"
                        )

                    u_b = _u_center_col(Ub, int(j), int(i))
                    u_n = _u_center_col(Un, int(j), int(i))
                    if source == "diff":
                        u_k = u_b - u_n
                    elif source == "bay":
                        u_k = u_b
                    elif source == "nobay":
                        u_k = u_n
                    else:
                        raise ValueError("source debe ser 'diff' | 'bay' | 'nobay'.")

                    t_lists[sidx].append(t_use)
                    U_lists[sidx].append(np.array(u_k[sl, :], dtype=float))
            finally:
                # Cierra la etapa tan pronto deja de usarse para reducir memoria/handles abiertos.
                _close_mnc(nc_b)
                _close_mnc(nc_n)

        t_full0 = np.concatenate(t_lists[0]) if t_lists[0] else np.array([])
        if t_full0.size == 0:
            raise RuntimeError("No se pudo construir el eje temporal concatenado para U/V.")

        tmin_seg = tmin_horas * 3600.0
        idx0 = np.where(t_full0 >= tmin_seg)[0][0] if np.any(t_full0 >= tmin_seg) else 0
        t_cut = t_full0[idx0:]
        t_h = t_cut / 3600.0
        nt = int(t_h.shape[0])

        U_hov = np.full((nt, nz, nst), np.nan, dtype=np.float32)

        for sidx in range(nst):
            U_full = np.concatenate(U_lists[sidx], axis=0) if U_lists[sidx] else np.empty((0, nz))
            U_cut = U_full[idx0:, :]
            if U_cut.shape != (nt, nz):
                raise ValueError(
                    f"Forma inesperada en estación {labels[sidx]}: "
                    f"U={U_cut.shape}, esperado {(nt, nz)}"
                )
            U_hov[:, :, sidx] = U_cut.astype(np.float32)

        # Liberar listas de U antes de iniciar V (baja pico de memoria).
        U_lists = None

        # ============================================================
        # Fase 2: construir Hovmöller de V (después de U)
        # ============================================================
        t_lists_v = [[] for _ in range(nst)]
        V_lists = [[] for _ in range(nst)]

        # Reabrimos etapas para V en streaming.
        for patt in stage_patterns:
            try:
                nc_b = _open_mnc_or_raise(outdir_bay, patt)
                nc_n = _open_mnc_or_raise(outdir_nobay, patt)
            except FileNotFoundError:
                continue
            try:
                tb = _get_time_from_nc(nc_b)
                tn = _get_time_from_nc(nc_n)
                Vb = nc_b.variables['V']
                Vn = nc_n.variables['V']

                if Vb.shape[0] != tb.shape[0] or Vn.shape[0] != tn.shape[0]:
                    raise ValueError("Longitudes tiempo/V no coinciden en una etapa.")
                if tb.shape != tn.shape or np.any(np.abs(tb - tn) > 1e-7):
                    raise ValueError("Los tiempos de bay y nobay difieren en una etapa.")

                trim0 = False
                if t_lists_v[0]:
                    trim0 = abs(tb[0] - t_lists_v[0][-1][-1]) <= tol_time
                sl = slice(1, None) if trim0 else slice(None)
                t_use = tb[sl]

                ny_v = int(Vb.shape[-2]) - 1  # V en yp1 -> centros en y
                nx_v = int(Vb.shape[-1])
                for sidx, (i, j) in enumerate(ij):
                    if not (0 <= j < ny_v and 0 <= i < nx_v):
                        raise IndexError(
                            f"Estación fuera de rango para V centrada: "
                            f"(i,j)=({i},{j}) vs (nx,ny)=({nx_v},{ny_v})"
                        )

                    v_b = _v_center_col(Vb, int(j), int(i))
                    v_n = _v_center_col(Vn, int(j), int(i))
                    if source == "diff":
                        v_k = v_b - v_n
                    elif source == "bay":
                        v_k = v_b
                    elif source == "nobay":
                        v_k = v_n
                    else:
                        raise ValueError("source debe ser 'diff' | 'bay' | 'nobay'.")

                    t_lists_v[sidx].append(t_use)
                    V_lists[sidx].append(np.array(v_k[sl, :], dtype=float))
            finally:
                _close_mnc(nc_b)
                _close_mnc(nc_n)

        V_hov = np.full((nt, nz, nst), np.nan, dtype=np.float32)
        for sidx in range(nst):
            t_full_v = np.concatenate(t_lists_v[sidx]) if t_lists_v[sidx] else np.array([])
            if t_full_v.size and (t_full_v.shape != t_full0.shape or np.any(np.abs(t_full_v - t_full0) > 1e-7)):
                raise ValueError(f"Tiempo de V difiere de U en estación {labels[sidx]}.")

            V_full = np.concatenate(V_lists[sidx], axis=0) if V_lists[sidx] else np.empty((0, nz))
            V_cut = V_full[idx0:, :]
            if V_cut.shape != (nt, nz):
                raise ValueError(
                    f"Forma inesperada en estación {labels[sidx]}: "
                    f"V={V_cut.shape}, esperado {(nt, nz)}"
                )
            V_hov[:, :, sidx] = V_cut.astype(np.float32)

        idx_100 = int(np.argmin(np.abs(RC + 100.0)))
        idx_400 = int(np.argmin(np.abs(RC + 400.0)))

        meta = {
            "scenario": scenario_label,
            "source": source,
            "var": "UV",
            "paths": {"bay": str(outdir_bay), "nobay": str(outdir_nobay)},
            "stage_patterns": list(stage_patterns),
            "tmin_horas": float(tmin_horas),
            "idx0": int(idx0),
            "dt_seconds": float(t_cut[1] - t_cut[0]) if len(t_cut) > 1 else None,
            "interp_to_cgrid": {
                "U": "0.5*(U[..., :-1] + U[..., 1:])  # Xp1 -> X",
                "V": "0.5*(V[..., :-1, :] + V[..., 1:, :])  # Yp1 -> Y",
            },
            "stations": [{"label": labels[k], "i": int(ij[k, 0]), "j": int(ij[k, 1]), "color": colors[k]} for k in range(nst)]
        }
        meta_json = json.dumps(meta, ensure_ascii=False)

        out_path = save_dir / f"{scenario_label}_hov_ALL_UV_{source}.npz"
        np.savez_compressed(
            out_path,
            t_horas=t_h.astype(np.float64),
            RC=RC.astype(np.float64),
            U_hov=U_hov,
            V_hov=V_hov,
            labels=np.array(labels, dtype=object),
            colors=np.array(colors, dtype=object),
            ij=ij.astype(np.int32),
            idx_100=np.array(idx_100, dtype=np.int32),
            idx_400=np.array(idx_400, dtype=np.int32),
            meta_json=np.array(meta_json)
        )
        print(f"[OK] Bundle Hovmöller U/V guardado: {out_path}")
        return out_path
    finally:
        # Cierre de respaldo ante errores/salidas tempranas.
        for nc in nc_bay_list:
            _close_mnc(nc)
        for nc in nc_nobay_list:
            _close_mnc(nc)
        _close_mnc(gr_bay)

def plot_hovmoller_from_bundle(
    npz_path: str,
    station: str | int = "A",          # etiqueta ("Seiche","A","B","C") o índice 0..nst-1
    cmap=cmo.cm.balance,
    vmin=None, vmax=None,
    title_prefix=''
):
    # Carga robusta: intenta sin pickle y, si hay object arrays, reintenta con pickle
    try:
        d = np.load(npz_path, allow_pickle=False)
    except ValueError:
        d = np.load(npz_path, allow_pickle=True)

    t   = d["t_horas"]
    RC  = d["RC"]
    Z3  = d["diff"]          # (nt, nz, nst)
    meta = json.loads(str(d["meta_json"]))

    # Recoge labels/colors de forma segura
    if "labels" in d.files:
        try:
            labs = d["labels"]  # puede ser object -> permitimos pickle arriba
        except ValueError:
            labs = np.array([s["label"] for s in meta.get("stations", [])], dtype="U16")
    else:
        labs = np.array([s["label"] for s in meta.get("stations", [])], dtype="U16")

    if "colors" in d.files:
        try:
            cols = d["colors"]
        except ValueError:
            cols = np.array([s.get("color") for s in meta.get("stations", [])], dtype="U16")
    else:
        cols = np.array([s.get("color") for s in meta.get("stations", [])], dtype="U16")

    scen= meta.get("scenario","")
    var = meta.get("var","W")

    # Seleccionar estación
    if isinstance(station, int):
        k = station
        if k < 0 or k >= Z3.shape[2]:
            raise IndexError(f"Índice de estación fuera de rango: {k}")
        label = str(labs[k])
    else:
        label = str(station)
        lab_list = [str(x) for x in labs]
        if label not in lab_list:
            raise KeyError(f"Estación '{label}' no encontrada. Disponibles: {lab_list}")
        k = lab_list.index(label)

    Z = Z3[:, :, k]  # (nt, nz)

    # Escala simétrica si no se pasa
    if (vmin is None) or (vmax is None):
        m = float(np.nanmax(np.abs(Z)))
        vmin, vmax = -m, m

    # Guías
    idx_100 = int(d["idx_100"]); z100 = RC[idx_100]
    idx_400 = int(d["idx_400"]); z400 = RC[idx_400]

    fig, ax = plt.subplots(figsize=(9,5))
    pcm = ax.pcolormesh(t, RC, Z.T, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')

    ax.hlines([z100, z400], xmin=t.min(), xmax=t.max(), colors='white',
              linestyles='--', linewidth=1.0, label='≈ 100 m y ≈ 400 m')

    fig.colorbar(pcm, ax=ax, label=f"{var} (bahía – sin bahía)")
    ax.set_xlabel("Tiempo (h)")
    ax.set_ylabel("Profundidad (m)")
    ax.set_title(f"{title_prefix} Hovmöller {var}_diff — {scen} • est. {label}")
    ax.legend(loc='upper right', fontsize='small')
    plt.tight_layout()
    plt.show()
    return fig, ax

# ------------------ selección de nivel por profundidad ------------------
def _select_k_for_depth(gr, depth_m: float, prefer_RF=False):
    """depth_m>0: metros bajo la superficie. Devuelve (k_index, z_ref_negativa)."""
    target = -abs(depth_m)
    if prefer_RF and 'RF' in gr.variables:
        z = np.array(gr.variables['RF'][:], dtype=float).squeeze()
    elif 'RC' in gr.variables:
        z = np.array(gr.variables['RC'][:], dtype=float).squeeze()
    elif 'Z' in gr.variables:
        z = np.array(gr.variables['Z'][:], dtype=float).squeeze()
    else:
        raise KeyError("No encontré coordenadas verticales ('RC'/'RF'/'Z').")
    if z.ndim > 1:
        z = z[..., 0, 0]
    k = int(np.argmin(np.abs(z - target)))
    return k, float(z[k])

def guardar_zoom_2d_npz(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    varname: str = 'W',            # 'W': (t,z,y,x) ; 'Eta': (t,y,x)
    depth_m: float | None = 110.0, # para 3D
    prefer_RF: bool = False,
    subtract: bool = True,
    source: str = "diff",          # "diff" (bahía−sin), "bay", "nobay"
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern   = "grid.t*.nc",
    X_sym: int = 150,
    Y_cut: int = 420
):
    if source not in ("diff", "bay", "nobay"):
        raise ValueError("source debe ser 'diff' | 'bay' | 'nobay'.")
    save_dir.mkdir(parents=True, exist_ok=True)

    # Abre grid una sola vez (y ciérralo al final)
    gr = _open_mnc_or_raise(outdir_bay, grid_pattern)

    # --- coordenadas/máscara y slices de zoom
    # --- Coordenadas y máscara (robusto) ---
    Depth = np.array(gr.variables['Depth'][:], dtype=float)  # (ny, nx)
    ny, nx = Depth.shape

    XC = _to_1d_x(gr, nx, 'XC', 'XG')   # (nx,)
    YC = _to_1d_y(gr, ny, 'YC', 'YG')   # (ny,)

    # Slices de zoom (clamp para no salir del dominio)
    x_center = nx // 2
    x0 = max(0, x_center - X_sym)
    x1 = min(nx, x_center + X_sym)
    y0 = max(0, ny - Y_cut)
    y1 = ny

    x_slice = slice(x0, x1)
    y_slice = slice(y0, y1)

    X_zoom = XC[x_slice]                        # (nx_zoom,)
    Y_zoom = YC[y_slice]                        # (ny_zoom,)

    # Máscara (tierra donde depth==0)
    mask_2d = np.ma.getmask(np.ma.masked_values(Depth, 0.0))   # bool (ny,nx)
    mask_zoom = mask_2d[y_slice, x_slice]                      # (ny_zoom, nx_zoom)

    # --- detectar rank/level una sola vez con la PRIMERA etapa
    # buscamos el primer patrón que exista en ambas carpetas
    first_nc_b = None; first_nc_n = None
    for patt in stage_patterns:
        try:
            first_nc_b = _open_mnc_or_raise(outdir_bay, patt)
            first_nc_n = _open_mnc_or_raise(outdir_nobay, patt)
            break
        except FileNotFoundError:
            continue
    if first_nc_b is None:
        _close_mnc(gr)
        raise FileNotFoundError("No se encontró ninguna etapa en bay/nobay.")

    vb0 = first_nc_b.variables[varname]
    shp0 = tuple(getattr(vb0, "shape", ()))
    if len(shp0) not in (3,4):
        _close_mnc(first_nc_b); _close_mnc(first_nc_n); _close_mnc(gr)
        raise ValueError(f"Forma no soportada para '{varname}': {shp0}")

    var_rank = len(shp0)
    need_level = None; z_ref = None
    if var_rank == 4:
        if depth_m is None:
            _close_mnc(first_nc_b); _close_mnc(first_nc_n); _close_mnc(gr)
            raise ValueError(f"'{varname}' es 3D; debes indicar depth_m.")
        # elegir k por profundidad
        target = -abs(depth_m)
        if prefer_RF and 'RF' in gr.variables:
            z = np.array(gr.variables['RF'][:], dtype=float).squeeze()
        elif 'RC' in gr.variables:
            z = np.array(gr.variables['RC'][:], dtype=float).squeeze()
        elif 'Z' in gr.variables:
            z = np.array(gr.variables['Z'][:], dtype=float).squeeze()
        else:
            _close_mnc(first_nc_b); _close_mnc(first_nc_n); _close_mnc(gr)
            raise KeyError("No hay RC/RF/Z en grid.")
        if z.ndim > 1: z = z[..., 0, 0]
        need_level = int(np.argmin(np.abs(z - target)))
        z_ref = float(z[need_level])

    # obtener tiempo de la primera etapa y comenzar listas
    t_ref = _get_time_from_nc(first_nc_b)
    # acumuladores (listas) que iremos extendiendo
    t_list = []
    data_list = []

    # procesar etapa por etapa (streaming)
    for patt in stage_patterns:
        # abrir ambos, leer, cerrar
        try:
            nc_b = _open_mnc_or_raise(outdir_bay, patt)
            nc_n = _open_mnc_or_raise(outdir_nobay, patt)
        except FileNotFoundError:
            continue
        try:
            tb = _get_time_from_nc(nc_b)
            tn = _get_time_from_nc(nc_n)

            vb = nc_b.variables[varname]
            vn = nc_n.variables[varname]

            if var_rank == 3:
                Fb = np.array(vb[:, y_slice, x_slice], dtype=float)  # (t, y, x)
                Fn = np.array(vn[:, y_slice, x_slice], dtype=float)
            else:
                Fb = np.array(vb[:, need_level, y_slice, x_slice], dtype=float)  # (t, y, x)
                Fn = np.array(vn[:, need_level, y_slice, x_slice], dtype=float)

            if Fb.shape[0] != tb.shape[0] or Fn.shape[0] != tn.shape[0]:
                raise ValueError("Longitudes tiempo/serie no coinciden en una etapa.")

            # diferencia o caso individual
            if source == "bay":
                Dk = Fb
            elif source == "nobay":
                Dk = Fn
            else:  # "diff"
                Dk = Fb - Fn if subtract else Fb

            # append, pero evitando duplicar el primer sample si coincide con el último
            if t_list:
                if abs(tb[0] - t_list[-1][-1]) <= 1e-9:
                    tb = tb[1:]; Dk = Dk[1:, ...]
            t_list.append(tb)
            data_list.append(Dk)
        finally:
            _close_mnc(nc_b)
            _close_mnc(nc_n)

    # cerrar los primeros y el grid
    _close_mnc(first_nc_b); _close_mnc(first_nc_n); _close_mnc(gr)

    if not t_list:
        raise RuntimeError("No se pudo construir la serie temporal del zoom.")

    # concatenar
    t_full = np.concatenate(t_list)              # (nt_full,)
    D_full = np.concatenate(data_list, axis=0)   # (nt_full, y, x)

    # recorte por tmin
    tmin_seg = tmin_horas * 3600.0
    idx0 = np.where(t_full >= tmin_seg)[0][0] if np.any(t_full >= tmin_seg) else 0
    t_horas = (t_full[idx0:] / 3600.0).astype(np.float64)
    D = D_full[idx0:, :, :].astype(np.float32)

    keyz = f"z{int(round(depth_m))}m" if (var_rank == 4 and depth_m is not None) else "sfc"
    src_suffix = "" if source == "diff" else f"_{source}"
    out_path = save_dir / f"{scenario_label}_{varname}_zoom_{keyz}{src_suffix}.npz"

    meta = {
        "scenario": scenario_label,
        "var": varname,
        "subtract": bool(subtract),
        "source": source,
        "depth_m": float(depth_m) if (var_rank == 4 and depth_m is not None) else None,
        "prefer_RF": bool(prefer_RF),
        "k_index": int(need_level) if need_level is not None else None,
        "z_ref": float(z_ref) if z_ref is not None else None,
        "tmin_horas": float(tmin_horas),
        "stage_patterns": list(stage_patterns),
        "zoom": {"X_sym": int(X_sym), "Y_cut": int(Y_cut)},
        "grid_pattern": grid_pattern
    }
    meta_json = json.dumps(meta, ensure_ascii=False)

    np.savez_compressed(
        out_path,
        t_horas=t_horas,
        X_zoom=X_zoom.astype(np.float64),
        Y_zoom=Y_zoom.astype(np.float64),
        mask_zoom=mask_zoom.astype(np.uint8),
        data=D,
        meta_json=np.array(meta_json)
    )
    print(f"[OK] Guardado bundle zoom: {out_path}")
    return str(out_path)


def guardar_zoom_vel_unstagger_npz(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    varname: str = "U",          # "U"(T,Z,Y,Xp1)  o  "V"(T,Z,Yp1,X)
    depth_m: float = 110.0,
    subtract: bool = True,
    source: str = "diff",        # "diff" (bahía−sin), "bay", "nobay"
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern   = "grid.t*.nc",
    X_sym: int = 150,
    Y_cut: int = 420,
) -> str:
    """
    Guarda un zoom de velocidad U o V en malla C (unstaggered) como NPZ.

    MITgcm Arakawa C-grid:
      U(T, Z, Y, Xp1) — Xp1 tiene nx+1 caras en x.
        Unstagger: U_c[i] = 0.5 * (U[i] + U[i+1])  →  forma (T, ny_zoom, nx_zoom)
      V(T, Z, Yp1, X) — Yp1 tiene ny+1 caras en y.
        Unstagger: V_c[j] = 0.5 * (V[j] + V[j+1])  →  forma (T, ny_zoom, nx_zoom)

    Coordenadas y máscara: malla C (tracer grid) XC, YC, Depth.

    Guarda NPZ con claves idénticas a guardar_zoom_2d_npz:
      t_horas, X_zoom, Y_zoom, mask_zoom, data, meta_json
    """
    var = str(varname).upper()
    if var not in ("U", "V"):
        raise ValueError("varname debe ser 'U' o 'V'.")
    if source not in ("diff", "bay", "nobay"):
        raise ValueError("source debe ser 'diff' | 'bay' | 'nobay'.")

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Grid: tracer C-grid coords & mask ────────────────────────────
    gr = _open_mnc_or_raise(outdir_bay, grid_pattern)
    try:
        Depth = np.array(gr.variables["Depth"][:], dtype=float)   # (ny, nx)
        ny, nx = Depth.shape

        def _x1d():
            for vn in ("XC", "xC", "X", "x", "XG", "xG"):
                if vn in gr.variables:
                    arr = np.array(gr.variables[vn][:], dtype=float)
                    if arr.size > 0:
                        if arr.ndim == 2:
                            return arr[0, :]
                        return arr[:-1] if arr.size == nx + 1 else arr
            return np.arange(nx, dtype=float)

        def _y1d():
            for vn in ("YC", "yC", "Y", "y", "YG", "yG"):
                if vn in gr.variables:
                    arr = np.array(gr.variables[vn][:], dtype=float)
                    if arr.size > 0:
                        if arr.ndim == 2:
                            return arr[:, 0]
                        return arr[:-1] if arr.size == ny + 1 else arr
            return np.arange(ny, dtype=float)

        XC = _x1d()   # (nx,)  – tracer cell centers
        YC = _y1d()   # (ny,)

        # Tracer-grid zoom bounds
        x_center = nx // 2
        x0 = max(0, x_center - X_sym)
        x1 = min(nx, x_center + X_sym)
        y0 = max(0, ny - Y_cut)
        y1 = ny

        X_zoom = XC[x0:x1]                                         # (nx_zoom,)
        Y_zoom = YC[y0:y1]                                         # (ny_zoom,)
        mask_zoom = np.ma.getmask(
            np.ma.masked_values(Depth[y0:y1, x0:x1], 0.0)
        )   # True=tierra, shape (ny_zoom, nx_zoom)
        if mask_zoom is np.ma.nomask:
            mask_zoom = np.zeros((y1 - y0, x1 - x0), dtype=bool)

        # ── 2. Depth level from first available stage ────────────────────
        first_nc_b = None
        first_nc_n = None
        for patt in stage_patterns:
            try:
                first_nc_b = _open_mnc_or_raise(outdir_bay, patt)
                first_nc_n = _open_mnc_or_raise(outdir_nobay, patt)
                break
            except FileNotFoundError:
                continue
        if first_nc_b is None:
            raise FileNotFoundError("No se encontró ninguna etapa en bay/nobay.")

        try:
            if var not in first_nc_b.variables:
                raise KeyError(f"Variable '{var}' no encontrada en los archivos de estado.")
            shp0 = tuple(first_nc_b.variables[var].shape)
            if len(shp0) != 4:
                raise ValueError(f"Se esperaba '{var}' 4D (T,Z,Y,X*). Encontrado: {shp0}")

            # Vertical coordinate: prefer Z (cell center) for U/V
            z_vec = None
            for zname in ("Z", "RC", "Zl", "RF"):
                for src in (first_nc_b, gr):
                    if zname in src.variables:
                        z_vec = np.array(src.variables[zname][:], dtype=float).squeeze()
                        if z_vec.ndim > 1:
                            z_vec = z_vec[..., 0, 0]
                        break
                if z_vec is not None:
                    break
            if z_vec is None:
                raise KeyError("No se encontró coordenada vertical Z/RC/Zl/RF en grid o estado.")

            k_level = int(np.argmin(np.abs(z_vec - (-abs(depth_m)))))
            z_ref   = float(z_vec[k_level])
        finally:
            _close_mnc(first_nc_b)
            _close_mnc(first_nc_n)

        # ── 3. Staggered read slices ─────────────────────────────────────
        # U(T,Z,Y,Xp1): y on tracer grid, x has nx+1 faces.
        #   Read Xp1[x0 : x1+1] → nx_zoom+1 values → avg → nx_zoom
        # V(T,Z,Yp1,X): x on tracer grid, y has ny+1 faces.
        #   Read Yp1[y0 : y1+1] → ny_zoom+1 values → avg → ny_zoom
        if var == "U":
            rd_y = slice(y0, y1)
            rd_x = slice(x0, x1 + 1)   # x1+1 <= nx+1 always (x1 <= nx)
        else:
            rd_y = slice(y0, y1 + 1)   # y1+1 <= ny+1 always (y1 <= ny)
            rd_x = slice(x0, x1)

        # ── 4. Streaming through stage files ────────────────────────────
        t_list    = []
        data_list = []

        for patt in stage_patterns:
            try:
                nc_b = _open_mnc_or_raise(outdir_bay, patt)
                nc_n = _open_mnc_or_raise(outdir_nobay, patt)
            except FileNotFoundError:
                continue
            try:
                tb = _get_time_from_nc(nc_b)
                tn = _get_time_from_nc(nc_n)
                if tb.shape != tn.shape or np.any(np.abs(tb - tn) > 1e-4):
                    raise ValueError(f"Timestamps difieren entre bay/nobay en '{patt}'.")

                # Read staggered slice at chosen depth level
                Fb_raw = np.array(nc_b.variables[var][:, k_level, rd_y, rd_x], dtype=float)
                Fn_raw = np.array(nc_n.variables[var][:, k_level, rd_y, rd_x], dtype=float)

                # Unstagger to tracer (C) grid
                if var == "U":
                    # (T, ny_zoom, nx_zoom+1) → (T, ny_zoom, nx_zoom)
                    Fb = 0.5 * (Fb_raw[:, :, :-1] + Fb_raw[:, :, 1:])
                    Fn = 0.5 * (Fn_raw[:, :, :-1] + Fn_raw[:, :, 1:])
                else:
                    # (T, ny_zoom+1, nx_zoom) → (T, ny_zoom, nx_zoom)
                    Fb = 0.5 * (Fb_raw[:, :-1, :] + Fb_raw[:, 1:, :])
                    Fn = 0.5 * (Fn_raw[:, :-1, :] + Fn_raw[:, 1:, :])

                if source == "bay":
                    Dk = Fb
                elif source == "nobay":
                    Dk = Fn
                else:  # "diff"
                    Dk = Fb - Fn if subtract else Fb

                # Drop duplicate first sample when stages overlap
                if t_list and abs(tb[0] - t_list[-1][-1]) <= max(1e-9, 1e-4):
                    tb = tb[1:]
                    Dk = Dk[1:, ...]

                t_list.append(tb)
                data_list.append(Dk)
            finally:
                _close_mnc(nc_b)
                _close_mnc(nc_n)

    finally:
        _close_mnc(gr)

    if not t_list:
        raise RuntimeError(f"No se pudo construir la serie temporal para '{var}'.")

    # ── 5. Concatenate, trim, save ───────────────────────────────────────
    t_full = np.concatenate(t_list)
    D_full = np.concatenate(data_list, axis=0)

    tmin_seg = tmin_horas * 3600.0
    idx0 = np.where(t_full >= tmin_seg)[0][0] if np.any(t_full >= tmin_seg) else 0
    t_horas = (t_full[idx0:] / 3600.0).astype(np.float64)
    D       = D_full[idx0:].astype(np.float32)

    keyz       = f"z{int(round(depth_m))}m"
    src_suffix = "" if source == "diff" else f"_{source}"
    out_path   = save_dir / f"{scenario_label}_{var}_zoom_{keyz}{src_suffix}.npz"

    meta = {
        "scenario":       scenario_label,
        "var":            var,
        "unstagger":      True,
        "subtract":       bool(subtract),
        "source":         source,
        "depth_m":        float(depth_m),
        "k_index":        int(k_level),
        "z_ref":          float(z_ref),
        "tmin_horas":     float(tmin_horas),
        "stage_patterns": list(stage_patterns),
        "zoom":           {"X_sym": int(X_sym), "Y_cut": int(Y_cut)},
        "grid_pattern":   grid_pattern,
    }

    np.savez_compressed(
        out_path,
        t_horas   = t_horas,
        X_zoom    = X_zoom.astype(np.float64),
        Y_zoom    = Y_zoom.astype(np.float64),
        mask_zoom = mask_zoom.astype(np.uint8),
        data      = D,
        meta_json = np.array(json.dumps(meta, ensure_ascii=False)),
    )
    print(f"[OK] Guardado bundle zoom ({var}, unstaggered → C-grid): {out_path}")
    return str(out_path)



def guardar_zoom_uw_npz(*args, **kwargs) -> str:
    """Alias retrocompatible → guardar_zoom_vel_unstagger_npz()."""
    return guardar_zoom_vel_unstagger_npz(*args, **kwargs)


def guardar_zoom_uv_npz(*args, **kwargs) -> str:
    """Alias retrocompatible → guardar_zoom_vel_unstagger_npz()."""
    return guardar_zoom_vel_unstagger_npz(*args, **kwargs)


# ------------------ animación desde el NPZ guardado ------------------
def animar_zoom_desde_npz(
    npz_path: str,
    outname: str = "anim_zoom.mp4",
    vmin: float | None = None,
    vmax: float | None = None,
    fps: int = 2,
    cmap = cmo.cm.balance,
    show_mask_contour: bool = True
):
    """
    Lee el bundle .npz del zoom y genera una animación MP4.
    """
    d = np.load(npz_path, allow_pickle=True)  # etiquetas/objetos no son necesarios, pero por robustez
    t_h = d["t_horas"]        # (nt,)
    Xz  = d["X_zoom"]         # (nx_zoom,)
    Yz  = d["Y_zoom"]         # (ny_zoom,)
    Mz  = d["mask_zoom"]      # (ny_zoom, nx_zoom) 0/1
    A   = d["data"]           # (nt, ny_zoom, nx_zoom)
    meta= json.loads(str(d["meta_json"]))

    # Límites de color
    if (vmin is None) or (vmax is None):
        m = float(np.nanmax(np.abs(A)))
        vmin, vmax = -m, m

    # Figura
    fig, ax = plt.subplots(figsize=(8,6))
    mesh = ax.pcolormesh(Xz, Yz, A[0, :, :], cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    cbar = fig.colorbar(mesh, ax=ax, label=f"{meta.get('var','W')}{' (bay−nobay)' if meta.get('subtract', True) else ' (bay)'}")
    if show_mask_contour:
        ax.contour(Xz, Yz, Mz, colors='black', linewidths=1.0, linestyles='--')
    ax.set_aspect('equal')
    title_base = f"{meta.get('scenario','')} • {meta.get('var','W')} • zoom"
    if meta.get("depth_m") is not None:
        title_base += f" • {int(round(meta['depth_m']))} m"
    ax.set_title(f"{title_base} • t={t_h[0]:.2f} h")
    plt.tight_layout()

    def update(frame):
        mesh.set_array(A[frame, :, :].ravel())
        ax.set_title(f"{title_base} • t={t_h[frame]:.2f} h")
        return (mesh,)

    anim = FuncAnimation(fig, update, frames=A.shape[0], blit=True)
    anim.save(outname, writer="ffmpeg", fps=fps)
    plt.close(fig)
    print(f"[OK] Animación exportada: {outname}")
    return outname


def animar_zoom_desde_npz_normalizado_estaciones(
    npz_path: str,
    outname: str = "anim_zoom_norm_stations.mp4",
    fps: int = 2,
    cmap = cmo.cm.balance,
    show_mask_contour: bool = True,
    outdir_grid: Path | str | None = None,
    grid_pattern: str | None = None,
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    marker_size: int = 36
):
    """
    Anima un bundle zoom NPZ normalizado por max(|data|) global y superpone estaciones.

    - Normalización: A_norm = A / max(abs(A))
    - Escala color por defecto: [-1, 1]
    - Estaciones: se obtienen con make_stations_from_grid()
    """
    d = np.load(npz_path, allow_pickle=True)
    t_h = d["t_horas"]        # (nt,)
    Xz  = d["X_zoom"]         # (nx_zoom,)
    Yz  = d["Y_zoom"]         # (ny_zoom,)
    Mz  = d["mask_zoom"]      # (ny_zoom, nx_zoom) 0/1
    A   = d["data"].astype(float)   # (nt, ny_zoom, nx_zoom)
    meta = json.loads(str(d["meta_json"]))

    # Normalización global por máximo absoluto
    amax = float(np.nanmax(np.abs(A))) if A.size else 0.0
    if not np.isfinite(amax) or amax <= 0.0:
        amax = 1.0
    A_norm = A / amax
    vmin, vmax = -1.0, 1.0

    # Resolver ruta de grid para stations
    grid_dir = None
    if outdir_grid is not None:
        grid_dir = Path(outdir_grid)
    else:
        bay_path = meta.get("paths", {}).get("bay")
        if bay_path:
            grid_dir = Path(bay_path)

    stations_plot = []
    if grid_dir is not None:
        gpatt = grid_pattern or meta.get("grid_pattern", "grid.t*.nc")
        gr = _open_mnc_or_raise(grid_dir, gpatt)
        try:
            stations = make_stations_from_grid(gr, j_centro=j_centro, radio=radio)
            depth = np.array(gr.variables["Depth"][:], dtype=float)
            ny, nx = depth.shape
            XC = _as_2d_from_grid(gr, ("XC", "xC", "X", "x", "XG", "xG"), axis="x", ny=ny, nx=nx)
            YC = _as_2d_from_grid(gr, ("YC", "yC", "Y", "y", "YG", "yG"), axis="y", ny=ny, nx=nx)
            use_index_coords = _is_index_axis(Xz, nx) and _is_index_axis(Yz, ny)

            xlo, xhi = float(np.min(Xz)), float(np.max(Xz))
            ylo, yhi = float(np.min(Yz)), float(np.max(Yz))
            for st in stations:
                i, j, lab, col = st
                col = col or "k"
                i, j = int(i), int(j)
                if use_index_coords:
                    xs, ys = float(i), float(j)
                else:
                    if XC is None or YC is None:
                        continue
                    xs, ys = float(XC[j, i]), float(YC[j, i])
                if (xlo <= xs <= xhi) and (ylo <= ys <= yhi):
                    stations_plot.append((xs, ys, str(lab), col))
        finally:
            _close_mnc(gr)

    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(Xz, Yz, A_norm[0], cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
    fig.colorbar(mesh, ax=ax, label=f"{meta.get('var', 'var')} normalizado (max|.|={amax:.3e})")
    if show_mask_contour:
        ax.contour(Xz, Yz, Mz, colors="black", linewidths=1.0, linestyles="--")

    # Dibujar estaciones
    for xs, ys, lab, col in stations_plot:
        ax.scatter(xs, ys, s=marker_size, c=col, edgecolors="k", linewidths=0.6, zorder=5)
        ax.text(xs, ys, f" {lab}", color=col, fontsize=9, va="center", ha="left", zorder=6)

    ax.set_aspect("equal")
    title_base = f"{meta.get('scenario', '')} • {meta.get('var', 'var')} • zoom norm."
    if meta.get("depth_m") is not None:
        title_base += f" • {int(round(meta['depth_m']))} m"
    ax.set_title(f"{title_base} • t={t_h[0]:.2f} h")
    plt.tight_layout()

    def update(frame):
        mesh.set_array(A_norm[frame].ravel())
        ax.set_title(f"{title_base} • t={t_h[frame]:.2f} h")
        return (mesh,)

    anim = FuncAnimation(fig, update, frames=A_norm.shape[0], blit=True)
    anim.save(outname, writer="ffmpeg", fps=fps)
    plt.close(fig)
    print(f"[OK] Animación exportada (normalizada+estaciones): {outname}")
    return outname

def graficar_temp_zoom_npz_norm_estacion_c(
    pattern_npz: str = "./npz_out/*_Temp_zoom_z110m.npz",
    step_frames: int = 60,
    max_cols: int = 4,
    levels: int = 101,
    station_label_norm: str = "C",
    land_color: str = "#ece2c6",
    cmap_name: str = "RdBu_r",
    sharex: bool = True,
    sharey: bool = True,
    outdir_grid: Path | str | None = None,
    grid_pattern: str | None = None,
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    show_mask_contour: bool = False,
    marker_size: int = 24,
    figsize_panel: tuple[float, float] = (4.1, 3.2),
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    save_dir: str | Path = "./figs",
    dpi: int = 170,
    save_fig: bool = True,
    show_fig: bool = False,
):
    """
    Genera figuras (una por NPZ) con paneles cada `step_frames` usando contourf.
    Requisitos implementados:
      - normalización por max(abs(.)) en la estación C (o `station_label_norm`)
      - máscara NaN/tierra color arena
      - sharex/sharey
      - barra de color [-1, 1] en la última columna
      - estaciones delimitadas con su color de make_stations_from_grid
    """
    npz_list = sorted(glob.glob(pattern_npz))
    if not npz_list:
        raise FileNotFoundError(f"No encontré NPZ con patrón: {pattern_npz}")
    if step_frames <= 0:
        raise ValueError("step_frames debe ser > 0.")
    if max_cols <= 0:
        raise ValueError("max_cols debe ser > 0.")

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    out_paths = []
    for npz_path in npz_list:
        with np.load(npz_path, allow_pickle=True) as d:
            t_h = d["t_horas"].astype(float)
            Xz = d["X_zoom"].astype(float)
            Yz = d["Y_zoom"].astype(float)
            Mz = d["mask_zoom"].astype(bool)  # True=tierra
            A = d["data"].astype(float)
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}

        scenario = str(meta.get("scenario", Path(npz_path).name.split("_")[0]))
        depth_txt = int(round(float(meta.get("depth_m", 110.0))))

        # Coordenadas en km para los ejes
        Xz_km = Xz / 1e3
        Yz_km = Yz / 1e3

        # Resolver grid_dir y estaciones visibles en el zoom por coordenada (XC/YC).
        grid_dir = None
        if outdir_grid is not None:
            grid_dir = Path(outdir_grid)
        else:
            bay_path = meta.get("paths", {}).get("bay")
            if bay_path:
                grid_dir = Path(bay_path)
            else:
                paths = get_iwbay_paths()
                if scenario in paths["scenarios"]:
                    grid_dir = Path(paths["scenarios"][scenario]) / paths["bay_dir"]

        stations_plot = []
        use_index_coords = False
        if grid_dir is not None:
            gpatt = grid_pattern or meta.get("grid_pattern", "grid.t*.nc")
            gr = _open_mnc_or_raise(grid_dir, str(gpatt))
            try:
                stations = make_stations_from_grid(gr, j_centro=j_centro, radio=radio)
                depth = np.array(gr.variables["Depth"][:], dtype=float)
                ny, nx = depth.shape
                XC = _as_2d_from_grid(gr, ("XC", "xC", "X", "x", "XG", "xG"), axis="x", ny=ny, nx=nx)
                YC = _as_2d_from_grid(gr, ("YC", "yC", "Y", "y", "YG", "yG"), axis="y", ny=ny, nx=nx)

                use_index_coords = _is_index_axis(Xz, nx) and _is_index_axis(Yz, ny)
                xlo, xhi = float(np.min(Xz)), float(np.max(Xz))
                ylo, yhi = float(np.min(Yz)), float(np.max(Yz))
                for st in stations:
                    i, j, lab, col = st
                    col = col or "k"
                    i, j = int(i), int(j)

                    if use_index_coords:
                        xs, ys = float(i), float(j)
                    else:
                        if XC is None or YC is None:
                            continue
                        xs, ys = float(XC[j, i]), float(YC[j, i])

                    if (xlo <= xs <= xhi) and (ylo <= ys <= yhi):
                        stations_plot.append((xs, ys, str(lab), str(col)))
            finally:
                _close_mnc(gr)

        st_norm = [s for s in stations_plot if s[2] == str(station_label_norm)]
        if st_norm:
            xs_n, ys_n = float(st_norm[0][0]), float(st_norm[0][1])
            ix = int(np.argmin(np.abs(Xz - xs_n)))
            iy = int(np.argmin(np.abs(Yz - ys_n)))
        else:
            ix, iy = A.shape[2] // 2, A.shape[1] // 2

        norm_factor = float(np.nanmax(np.abs(A[:, iy, ix])))
        if not np.isfinite(norm_factor) or norm_factor <= 0.0:
            norm_factor = 1.0
        A_norm = A / norm_factor
        print(
            f"[{scenario}] Normalización est. {station_label_norm}  →  "
            f"max|ΔTemp| = {norm_factor:.4e} °C"
        )

        idxs = np.arange(0, A_norm.shape[0], step_frames, dtype=int)
        nplots = int(idxs.size)
        if nplots == 0:
            continue

        ncols = min(max_cols, nplots)
        nrows = int(math.ceil(nplots / ncols))
        fig = plt.figure(
            figsize=(figsize_panel[0] * (ncols + 0.2), figsize_panel[1] * nrows)
        )
        gs = fig.add_gridspec(
            nrows=nrows,
            ncols=ncols + 1,
            width_ratios=[1.0] * ncols + [0.05],
            wspace=0.15,
            hspace=0.22,
        )
        axes = np.empty((nrows, ncols), dtype=object)
        caxes = np.empty((nrows,), dtype=object)
        ax_ref = None
        for rr in range(nrows):
            for cc in range(ncols):
                share_ax = ax_ref if (ax_ref is not None and (sharex or sharey)) else None
                ax = fig.add_subplot(
                    gs[rr, cc],
                    sharex=share_ax if sharex else None,
                    sharey=share_ax if sharey else None,
                )
                axes[rr, cc] = ax
                if ax_ref is None:
                    ax_ref = ax
            caxes[rr] = fig.add_subplot(gs[rr, -1])

        cmap = plt.get_cmap(cmap_name).copy()
        cmap.set_bad(land_color)
        lvls = np.linspace(-1.0, 1.0, int(levels))

        # Rango visible (en km) para filtro de estaciones y límites de ejes
        if xlim is not None:
            xlo_vis, xhi_vis = float(xlim[0]) / 1e3, float(xlim[1]) / 1e3
        else:
            xlo_vis, xhi_vis = float(np.min(Xz_km)), float(np.max(Xz_km))
        if ylim is not None:
            ylo_vis, yhi_vis = float(ylim[0]) / 1e3, float(ylim[1]) / 1e3
        else:
            ylo_vis, yhi_vis = float(np.min(Yz_km)), float(np.max(Yz_km))

        row_mappables = [None] * nrows
        for k, tidx in enumerate(idxs):
            r = k // ncols
            c = k % ncols
            ax = axes[r, c]

            Z = np.where(Mz, np.nan, A_norm[tidx])
            mappable = ax.contourf(
                Xz_km, Yz_km, Z,
                levels=lvls, cmap=cmap, vmin=-1.0, vmax=1.0, extend="both"
            )
            if row_mappables[r] is None:
                row_mappables[r] = mappable
            ax.set_facecolor(land_color)

            if show_mask_contour:
                ax.contour(Xz_km, Yz_km, Mz.astype(float), colors="black",
                           linewidths=0.8, linestyles="--")

            for xs, ys, lab, col in stations_plot:
                xs_km, ys_km = xs / 1e3, ys / 1e3
                if not (xlo_vis <= xs_km <= xhi_vis and ylo_vis <= ys_km <= yhi_vis):
                    continue
                ax.scatter(xs_km, ys_km, s=marker_size, c=col,
                           edgecolors="k", linewidths=0.45, zorder=6)
                ax.text(xs_km, ys_km, f" {lab}", color=col, fontsize=8,
                        ha="left", va="center", zorder=7, clip_on=True)

            ax.set_title(f"t={t_h[tidx]:.1f} h (idx={int(tidx)})", fontsize=9)

        # Ocultar paneles vacíos
        for k in range(nplots, nrows * ncols):
            r = k // ncols
            c = k % ncols
            axes[r, c].axis("off")

        # Aplicar límites via ax_ref (se propaga a todos por sharex/sharey)
        ax_ref.set_xlim(xlo_vis, xhi_vis)
        ax_ref.set_ylim(ylo_vis, yhi_vis)

        # Etiquetas y tick visibility
        for r in range(nrows):
            for c in range(ncols):
                ax = axes[r, c]
                if r < nrows - 1:
                    ax.tick_params(labelbottom=False)
                    ax.set_xlabel("")
                else:
                    ax.set_xlabel("X (km)")
                if c == 0:
                    ax.set_ylabel("Y (km)")
                else:
                    ax.tick_params(labelleft=False)
                    ax.set_ylabel("")

        for r in range(nrows):
            cax = caxes[r]
            if row_mappables[r] is None:
                cax.axis("off")
                continue
            cb = fig.colorbar(row_mappables[r], cax=cax, ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
            cb.set_label(
                f"ΔTemp / max|ΔTemp|@{station_label_norm}\n"
                f"(max = {norm_factor:.3e} °C)",
                fontsize=8,
            )
            cb.ax.tick_params(labelsize=8)

        fig.suptitle(
            f"{scenario} • ΔTemp zoom z={depth_txt} m • "
            f"cada {step_frames} pasos • "
            f"norm: max|ΔTemp|@{station_label_norm} = {norm_factor:.3e} °C",
            fontsize=11, y=0.998,
        )
        fig.subplots_adjust(top=0.96)

        out_png = (
            save_dir
            / f"{scenario}_Temp_zoom_z{depth_txt}m_norm_{station_label_norm}_step{step_frames}_contourf.png"
        )
        if save_fig:
            fig.savefig(out_png, dpi=dpi)
            out_paths.append(str(out_png))
            print(f"[OK] Figura guardada: {out_png}")
        if show_fig:
            plt.show()
        plt.close(fig)

    return out_paths


def graficar_zoom_npz_norm_estacion_c(
    pattern_npz: str,
    var_tag: str,
    var_label: str,
    var_units: str = "",
    depth_tag: str | None = None,
    tmin_h: float | None = None,
    tmax_h: float | None = None,
    step_frames: int = 60,
    max_cols: int = 3,
    levels: int = 101,
    station_label_norm: str = "C",
    land_color: str = "#ece2c6",
    cmap_name: str = "RdBu_r",
    sharex: bool = True,
    sharey: bool = True,
    outdir_grid: Path | str | None = None,
    grid_pattern: str | None = None,
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    show_mask_contour: bool = False,
    marker_size: int = 24,
    figsize_panel: tuple[float, float] = (4.1, 3.2),
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    save_dir: str | Path = "./figs",
    dpi: int = 170,
    save_fig: bool = True,
    show_fig: bool = False,
):
    """
    Versión genérica de graficar_temp_zoom_npz_norm_estacion_c().
    - `var_tag`     se usa en el nombre de archivo: <scen>_<var_tag>_zoom_<depth_tag>_...
    - `var_label`   se usa en título y colorbar (ej. "ΔW", "Δη").
    - `depth_tag`   se usa en título/archivo. Si None, se infiere de meta.depth_m
                    (p. ej. "z110m"); si no existe, usa "sfc".
    """
    npz_list = sorted(glob.glob(pattern_npz))
    if not npz_list:
        raise FileNotFoundError(f"No encontré NPZ con patrón: {pattern_npz}")
    if step_frames <= 0:
        raise ValueError("step_frames debe ser > 0.")
    if max_cols <= 0:
        raise ValueError("max_cols debe ser > 0.")

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    out_paths = []
    for npz_path in npz_list:
        with np.load(npz_path, allow_pickle=True) as d:
            t_h = d["t_horas"].astype(float)
            Xz = d["X_zoom"].astype(float)
            Yz = d["Y_zoom"].astype(float)
            Mz = d["mask_zoom"].astype(bool)
            A = d["data"].astype(float)
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}

        scenario = str(meta.get("scenario", Path(npz_path).name.split("_")[0]))
        if depth_tag is None:
            if "depth_m" in meta:
                depth_tag_use = f"z{int(round(float(meta['depth_m'])))}m"
            else:
                depth_tag_use = "sfc"
        else:
            depth_tag_use = depth_tag

        Xz_km = Xz / 1e3
        Yz_km = Yz / 1e3

        grid_dir = None
        if outdir_grid is not None:
            grid_dir = Path(outdir_grid)
        else:
            bay_path = meta.get("paths", {}).get("bay")
            if bay_path:
                grid_dir = Path(bay_path)
            else:
                paths = get_iwbay_paths()
                if scenario in paths["scenarios"]:
                    grid_dir = Path(paths["scenarios"][scenario]) / paths["bay_dir"]

        stations_plot = []
        use_index_coords = False
        if grid_dir is not None:
            gpatt = grid_pattern or meta.get("grid_pattern", "grid.t*.nc")
            gr = _open_mnc_or_raise(grid_dir, str(gpatt))
            try:
                stations = make_stations_from_grid(gr, j_centro=j_centro, radio=radio)
                depth = np.array(gr.variables["Depth"][:], dtype=float)
                ny, nx = depth.shape
                XC = _as_2d_from_grid(gr, ("XC", "xC", "X", "x", "XG", "xG"), axis="x", ny=ny, nx=nx)
                YC = _as_2d_from_grid(gr, ("YC", "yC", "Y", "y", "YG", "yG"), axis="y", ny=ny, nx=nx)

                use_index_coords = _is_index_axis(Xz, nx) and _is_index_axis(Yz, ny)
                if use_index_coords and XC is not None and YC is not None:
                    ix_idx = np.clip(np.round(Xz).astype(int), 0, nx - 1)
                    iy_idx = np.clip(np.round(Yz).astype(int), 0, ny - 1)
                    Xz = XC[iy_idx[0], ix_idx].astype(float)
                    Yz = YC[iy_idx, ix_idx[0]].astype(float)
                    Xz_km = Xz / 1e3
                    Yz_km = Yz / 1e3
                    use_index_coords = False
                xlo, xhi = float(np.min(Xz)), float(np.max(Xz))
                ylo, yhi = float(np.min(Yz)), float(np.max(Yz))
                for st in stations:
                    i, j, lab, col = st
                    col = col or "k"
                    i, j = int(i), int(j)
                    if use_index_coords:
                        xs, ys = float(i), float(j)
                    else:
                        if XC is None or YC is None:
                            continue
                        xs, ys = float(XC[j, i]), float(YC[j, i])
                    if (xlo <= xs <= xhi) and (ylo <= ys <= yhi):
                        stations_plot.append((xs, ys, str(lab), str(col)))
            finally:
                _close_mnc(gr)

        st_norm = [s for s in stations_plot if s[2] == str(station_label_norm)]
        if st_norm:
            xs_n, ys_n = float(st_norm[0][0]), float(st_norm[0][1])
            ix = int(np.argmin(np.abs(Xz - xs_n)))
            iy = int(np.argmin(np.abs(Yz - ys_n)))
        else:
            ix, iy = A.shape[2] // 2, A.shape[1] // 2

        norm_factor = float(np.nanmax(np.abs(A[:, iy, ix])))
        if not np.isfinite(norm_factor) or norm_factor <= 0.0:
            norm_factor = 1.0
        A_norm = A / norm_factor
        print(
            f"[{scenario}] Normalización est. {station_label_norm}  →  "
            f"max|{var_label}| = {norm_factor:.4e}"
        )

        i0 = int(np.argmin(np.abs(t_h - tmin_h))) if tmin_h is not None else 0
        i1 = int(np.argmin(np.abs(t_h - tmax_h))) + 1 if tmax_h is not None else A_norm.shape[0]
        i0 = max(0, min(i0, A_norm.shape[0] - 1))
        i1 = max(i0 + 1, min(i1, A_norm.shape[0]))
        idxs = np.arange(i0, i1, step_frames, dtype=int)
        nplots = int(idxs.size)
        if nplots == 0:
            continue

        ncols = min(max_cols, nplots)
        nrows = int(math.ceil(nplots / ncols))
        fig = plt.figure(
            figsize=(figsize_panel[0] * (ncols + 0.2), figsize_panel[1] * nrows)
        )
        gs = fig.add_gridspec(
            nrows=nrows,
            ncols=ncols + 1,
            width_ratios=[1.0] * ncols + [0.05],
            wspace=0.15,
            hspace=0.22,
        )
        axes = np.empty((nrows, ncols), dtype=object)
        caxes = np.empty((nrows,), dtype=object)
        ax_ref = None
        for rr in range(nrows):
            for cc in range(ncols):
                share_ax = ax_ref if (ax_ref is not None and (sharex or sharey)) else None
                ax = fig.add_subplot(
                    gs[rr, cc],
                    sharex=share_ax if sharex else None,
                    sharey=share_ax if sharey else None,
                )
                axes[rr, cc] = ax
                if ax_ref is None:
                    ax_ref = ax
            caxes[rr] = fig.add_subplot(gs[rr, -1])

        cmap = plt.get_cmap(cmap_name).copy()
        cmap.set_bad(land_color)
        lvls = np.linspace(-1.0, 1.0, int(levels))

        if xlim is not None:
            xlo_vis, xhi_vis = float(xlim[0]) / 1e3, float(xlim[1]) / 1e3
        else:
            xlo_vis, xhi_vis = float(np.min(Xz_km)), float(np.max(Xz_km))
        if ylim is not None:
            ylo_vis, yhi_vis = float(ylim[0]) / 1e3, float(ylim[1]) / 1e3
        else:
            ylo_vis, yhi_vis = float(np.min(Yz_km)), float(np.max(Yz_km))

        row_mappables = [None] * nrows
        for k, tidx in enumerate(idxs):
            r = k // ncols
            c = k % ncols
            ax = axes[r, c]
            Z = np.where(Mz, np.nan, A_norm[tidx])
            mappable = ax.contourf(
                Xz_km, Yz_km, Z,
                levels=lvls, cmap=cmap, vmin=-1.0, vmax=1.0, extend="both"
            )
            if row_mappables[r] is None:
                row_mappables[r] = mappable
            ax.set_facecolor(land_color)
            if show_mask_contour:
                ax.contour(Xz_km, Yz_km, Mz.astype(float), colors="black",
                           linewidths=0.8, linestyles="--")
            for xs, ys, lab, col in stations_plot:
                xs_km, ys_km = xs / 1e3, ys / 1e3
                if not (xlo_vis <= xs_km <= xhi_vis and ylo_vis <= ys_km <= yhi_vis):
                    continue
                ax.scatter(xs_km, ys_km, s=marker_size, c=col,
                           edgecolors="k", linewidths=0.45, zorder=6)
                ax.text(xs_km, ys_km, f" {lab}", color=col, fontsize=8,
                        ha="left", va="center", zorder=7, clip_on=True)
            ax.set_title(f"t={t_h[tidx]:.1f} h (idx={int(tidx)})", fontsize=9)

        for k in range(nplots, nrows * ncols):
            r = k // ncols
            c = k % ncols
            axes[r, c].axis("off")

        ax_ref.set_xlim(xlo_vis, xhi_vis)
        ax_ref.set_ylim(ylo_vis, yhi_vis)

        for r in range(nrows):
            for c in range(ncols):
                ax = axes[r, c]
                if r < nrows - 1:
                    ax.tick_params(labelbottom=False)
                    ax.set_xlabel("")
                else:
                    ax.set_xlabel("X (km)")
                if c == 0:
                    ax.set_ylabel("Y (km)")
                else:
                    ax.tick_params(labelleft=False)
                    ax.set_ylabel("")

        for r in range(nrows):
            cax = caxes[r]
            if row_mappables[r] is None:
                cax.axis("off")
                continue
            cb = fig.colorbar(row_mappables[r], cax=cax, ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
            u_txt = f" {var_units}" if var_units else ""
            cb.set_label(
                f"{var_label} / max|{var_label}|@{station_label_norm}\n"
                f"(max = {norm_factor:.3e}{u_txt})",
                fontsize=8,
            )
            cb.ax.tick_params(labelsize=8)

        depth_title = depth_tag_use if depth_tag_use != "sfc" else "superficie"
        u_txt = f" {var_units}" if var_units else ""
        fig.suptitle(
            f"{scenario} • {var_label} zoom {depth_title} • "
            f"cada {step_frames} pasos • "
            f"norm: max|{var_label}|@{station_label_norm} = {norm_factor:.3e}{u_txt}",
            fontsize=11, y=0.998,
        )
        fig.subplots_adjust(top=0.96)

        win_tag = ""
        if tmin_h is not None or tmax_h is not None:
            a = f"{tmin_h:g}" if tmin_h is not None else "0"
            b = f"{tmax_h:g}" if tmax_h is not None else "end"
            win_tag = f"_t{a}-{b}h"
        out_png = (
            save_dir
            / f"{scenario}_{var_tag}_zoom_{depth_tag_use}_norm_{station_label_norm}_step{step_frames}{win_tag}_contourf.png"
        )
        if save_fig:
            fig.savefig(out_png, dpi=dpi)
            out_paths.append(str(out_png))
            print(f"[OK] Figura guardada: {out_png}")
        if show_fig:
            plt.show()
        plt.close(fig)

    return out_paths


def graficar_zoom_npz_comparativo_escenarios(
    pattern_npz: str,
    var_tag: str,
    var_label: str,
    var_units: str = "",
    depth_tag: str | None = None,
    tmin_h: float | None = None,
    tmax_h: float | None = None,
    step_frames: int = 5,
    scenario_order: tuple[str, ...] = ("2layer", "linear", "realistic"),
    station_label_norm: str = "C",
    levels: int = 101,
    land_color: str = "#ece2c6",
    cmap_name: str = "RdBu_r",
    outdir_grid: Path | str | None = None,
    grid_pattern: str | None = None,
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    marker_size: int = 22,
    figsize_panel: tuple[float, float] = (3.6, 3.0),
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    save_dir: str | Path = "./figs",
    dpi: int = 170,
    save_fig: bool = True,
    show_fig: bool = False,
    transpose: bool = False,
):
    """Una sola figura con filas=escenarios, columnas=tiempos. Norm. por escenario.

    Si ``transpose=True``: filas=tiempos, columnas=escenarios (formato vertical
    apto para insertar en documentos tipo tesis). Colorbars por escenario en la
    parte inferior.

    Cada fila tiene su propio colorbar con el factor de normalización
    `max|var|@station + var_units` para que la comparación entre escenarios sea
    explícita en patrón (norm.) y en amplitud (etiqueta del colorbar).
    """
    npz_list = sorted(glob.glob(pattern_npz))
    if not npz_list:
        raise FileNotFoundError(f"No encontré NPZ con patrón: {pattern_npz}")
    if step_frames <= 0:
        raise ValueError("step_frames debe ser > 0.")

    # Agrupar por escenario
    by_scen = {}
    for p in npz_list:
        with np.load(p, allow_pickle=True) as d:
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}
        scen = str(meta.get("scenario", Path(p).name.split("_")[0]))
        by_scen[scen] = p

    scen_used = [s for s in scenario_order if s in by_scen]
    if not scen_used:
        raise ValueError(
            f"Ningún escenario de {scenario_order} encontrado. "
            f"Disponibles: {list(by_scen)}"
        )

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # ── Carga + resolución de ejes/estaciones/norm por escenario ──────────
    per_scen = {}
    idxs_ref = None
    t_h_ref  = None
    depth_tag_final = depth_tag

    for scen in scen_used:
        with np.load(by_scen[scen], allow_pickle=True) as d:
            t_h  = d["t_horas"].astype(float)
            Xz   = d["X_zoom"].astype(float)
            Yz   = d["Y_zoom"].astype(float)
            Mz   = d["mask_zoom"].astype(bool)
            A    = d["data"].astype(float)
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}

        if depth_tag_final is None:
            depth_tag_final = (
                f"z{int(round(float(meta['depth_m'])))}m" if "depth_m" in meta else "sfc"
            )

        # Resolver grid_dir
        grid_dir = None
        if outdir_grid is not None:
            grid_dir = Path(outdir_grid)
        else:
            bay_path = meta.get("paths", {}).get("bay")
            if bay_path:
                grid_dir = Path(bay_path)
            else:
                paths = get_iwbay_paths()
                if scen in paths["scenarios"]:
                    grid_dir = Path(paths["scenarios"][scen]) / paths["bay_dir"]

        stations_plot = []
        if grid_dir is not None:
            gpatt = grid_pattern or meta.get("grid_pattern", "grid.t*.nc")
            gr = _open_mnc_or_raise(grid_dir, str(gpatt))
            try:
                stations = make_stations_from_grid(gr, j_centro=j_centro, radio=radio)
                depth = np.array(gr.variables["Depth"][:], dtype=float)
                ny, nx = depth.shape
                XC = _as_2d_from_grid(gr, ("XC","xC","X","x","XG","xG"), axis="x", ny=ny, nx=nx)
                YC = _as_2d_from_grid(gr, ("YC","yC","Y","y","YG","yG"), axis="y", ny=ny, nx=nx)

                if _is_index_axis(Xz, nx) and _is_index_axis(Yz, ny) \
                        and XC is not None and YC is not None:
                    ix_idx = np.clip(np.round(Xz).astype(int), 0, nx - 1)
                    iy_idx = np.clip(np.round(Yz).astype(int), 0, ny - 1)
                    Xz = XC[iy_idx[0], ix_idx].astype(float)
                    Yz = YC[iy_idx, ix_idx[0]].astype(float)

                xlo, xhi = float(np.min(Xz)), float(np.max(Xz))
                ylo, yhi = float(np.min(Yz)), float(np.max(Yz))
                for (i, j, lab, col) in stations:
                    i, j = int(i), int(j)
                    if XC is None or YC is None:
                        continue
                    xs, ys = float(XC[j, i]), float(YC[j, i])
                    if xlo <= xs <= xhi and ylo <= ys <= yhi:
                        stations_plot.append((xs, ys, str(lab), str(col or "k")))
            finally:
                _close_mnc(gr)

        # Estación de normalización
        st_norm = [s for s in stations_plot if s[2] == str(station_label_norm)]
        if st_norm:
            xs_n, ys_n = float(st_norm[0][0]), float(st_norm[0][1])
            ix = int(np.argmin(np.abs(Xz - xs_n)))
            iy = int(np.argmin(np.abs(Yz - ys_n)))
        else:
            ix, iy = A.shape[2] // 2, A.shape[1] // 2

        # Norm factor sobre la ventana visible
        i0 = int(np.argmin(np.abs(t_h - tmin_h))) if tmin_h is not None else 0
        i1 = int(np.argmin(np.abs(t_h - tmax_h))) + 1 if tmax_h is not None else A.shape[0]
        i0 = max(0, min(i0, A.shape[0] - 1))
        i1 = max(i0 + 1, min(i1, A.shape[0]))

        norm_factor = float(np.nanmax(np.abs(A[i0:i1, iy, ix])))
        if not np.isfinite(norm_factor) or norm_factor <= 0.0:
            norm_factor = 1.0

        idxs = np.arange(i0, i1, step_frames, dtype=int)
        if idxs.size == 0:
            raise ValueError(f"[{scen}] ventana temporal sin pasos seleccionables.")

        per_scen[scen] = dict(
            t_h=t_h, Xz=Xz, Yz=Yz, Mz=Mz, A=A,
            stations=stations_plot, norm_factor=norm_factor, idxs=idxs,
        )
        if idxs_ref is None:
            idxs_ref = idxs
            t_h_ref  = t_h
        print(
            f"[{scen}] Norm @ {station_label_norm}  →  "
            f"max|{var_label}| = {norm_factor:.4e}"
            + (f" {var_units}" if var_units else "")
        )

    # Sanity: igualar nº de pasos temporales al mínimo
    n_times = min(len(per_scen[s]["idxs"]) for s in scen_used)
    n_scen  = len(scen_used)

    if transpose:
        n_rows, n_cols = n_times, n_scen
    else:
        n_rows, n_cols = n_scen, n_times

    # ── Figura ────────────────────────────────────────────────────────────
    if transpose:
        # filas=tiempos, columnas=escenarios; colorbars en fila inferior
        fig = plt.figure(
            figsize=(figsize_panel[0] * n_cols, figsize_panel[1] * (n_rows + 0.25))
        )
        gs = fig.add_gridspec(
            nrows=n_rows + 1, ncols=n_cols,
            height_ratios=[1.0] * n_rows + [0.06],
            wspace=0.10, hspace=0.18,
        )
        axes = np.empty((n_rows, n_cols), dtype=object)
        caxes = np.empty((n_cols,), dtype=object)
        ax_ref = None
        for rr in range(n_rows):
            for cc in range(n_cols):
                ax = fig.add_subplot(
                    gs[rr, cc],
                    sharex=(ax_ref if ax_ref is not None else None),
                    sharey=(ax_ref if ax_ref is not None else None),
                )
                axes[rr, cc] = ax
                if ax_ref is None:
                    ax_ref = ax
        for cc in range(n_cols):
            caxes[cc] = fig.add_subplot(gs[-1, cc])
    else:
        fig = plt.figure(
            figsize=(figsize_panel[0] * (n_cols + 0.25), figsize_panel[1] * n_rows)
        )
        gs = fig.add_gridspec(
            nrows=n_rows, ncols=n_cols + 1,
            width_ratios=[1.0] * n_cols + [0.06],
            wspace=0.10, hspace=0.18,
        )
        axes = np.empty((n_rows, n_cols), dtype=object)
        caxes = np.empty((n_rows,), dtype=object)
        ax_ref = None
        for rr in range(n_rows):
            for cc in range(n_cols):
                ax = fig.add_subplot(
                    gs[rr, cc],
                    sharex=(ax_ref if ax_ref is not None else None),
                    sharey=(ax_ref if ax_ref is not None else None),
                )
                axes[rr, cc] = ax
                if ax_ref is None:
                    ax_ref = ax
            caxes[rr] = fig.add_subplot(gs[rr, -1])

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(land_color)
    lvls = np.linspace(-1.0, 1.0, int(levels))

    # Límites visibles
    Xz_km_any = per_scen[scen_used[0]]["Xz"] / 1e3
    Yz_km_any = per_scen[scen_used[0]]["Yz"] / 1e3
    if xlim is not None:
        xlo_vis, xhi_vis = float(xlim[0]) / 1e3, float(xlim[1]) / 1e3
    else:
        xlo_vis, xhi_vis = float(Xz_km_any.min()), float(Xz_km_any.max())
    if ylim is not None:
        ylo_vis, yhi_vis = float(ylim[0]) / 1e3, float(ylim[1]) / 1e3
    else:
        ylo_vis, yhi_vis = float(Yz_km_any.min()), float(Yz_km_any.max())

    scen_mappable = [None] * n_scen
    for s_i, scen in enumerate(scen_used):
        d = per_scen[scen]
        Xkm = d["Xz"] / 1e3
        Ykm = d["Yz"] / 1e3
        idxs = d["idxs"][:n_times]
        A_norm = d["A"] / d["norm_factor"]

        for t_i, tidx in enumerate(idxs):
            if transpose:
                rr, cc = t_i, s_i
            else:
                rr, cc = s_i, t_i
            ax = axes[rr, cc]
            Z = np.where(d["Mz"], np.nan, A_norm[tidx])
            mp = ax.contourf(
                Xkm, Ykm, Z, levels=lvls, cmap=cmap,
                vmin=-1.0, vmax=1.0, extend="both",
            )
            if scen_mappable[s_i] is None:
                scen_mappable[s_i] = mp
            ax.set_facecolor(land_color)

            for (xs, ys, lab, col) in d["stations"]:
                xs_km, ys_km = xs / 1e3, ys / 1e3
                if not (xlo_vis <= xs_km <= xhi_vis and ylo_vis <= ys_km <= yhi_vis):
                    continue
                ax.scatter(xs_km, ys_km, s=marker_size, c=col,
                           edgecolors="k", linewidths=0.4, zorder=6)
                if s_i == 0 and t_i == 0:
                    ax.text(xs_km, ys_km, f" {lab}", color=col, fontsize=7,
                            ha="left", va="center", zorder=7, clip_on=True)

            # Títulos de fila/columna superior
            if transpose:
                if t_i == 0:
                    ax.set_title(f"$\\bf{{{scen}}}$", fontsize=10)
            else:
                if s_i == 0:
                    ax.set_title(f"t={d['t_h'][tidx]:.2f} h", fontsize=9)

    ax_ref.set_xlim(xlo_vis, xhi_vis)
    ax_ref.set_ylim(ylo_vis, yhi_vis)

    # Tick visibility + labels de eje
    if transpose:
        # Primera columna lleva "t=X.XX h\nY (km)", última fila lleva "X (km)"
        t_h_first = per_scen[scen_used[0]]["t_h"]
        idxs_first = per_scen[scen_used[0]]["idxs"][:n_times]
        for rr in range(n_rows):
            for cc in range(n_cols):
                ax = axes[rr, cc]
                if rr < n_rows - 1:
                    ax.tick_params(labelbottom=False)
                    ax.set_xlabel("")
                else:
                    ax.set_xlabel("X (km)")
                if cc == 0:
                    ax.set_ylabel(
                        f"$\\bf{{t={t_h_first[idxs_first[rr]]:.2f}\\ h}}$\nY (km)",
                        fontsize=9,
                    )
                else:
                    ax.tick_params(labelleft=False)
                    ax.set_ylabel("")
    else:
        for rr in range(n_rows):
            for cc in range(n_cols):
                ax = axes[rr, cc]
                if rr < n_rows - 1:
                    ax.tick_params(labelbottom=False)
                    ax.set_xlabel("")
                else:
                    ax.set_xlabel("X (km)")
                if cc == 0:
                    ax.set_ylabel(f"$\\bf{{{scen_used[rr]}}}$\nY (km)", fontsize=9)
                else:
                    ax.tick_params(labelleft=False)
                    ax.set_ylabel("")

    # Colorbars por escenario con norm + unidades
    u_txt = f" {var_units}" if var_units else ""
    for s_i, scen in enumerate(scen_used):
        cax = caxes[s_i]
        if scen_mappable[s_i] is None:
            cax.axis("off")
            continue
        cb = fig.colorbar(
            scen_mappable[s_i], cax=cax,
            ticks=[-1.0, -0.5, 0.0, 0.5, 1.0],
            orientation=("horizontal" if transpose else "vertical"),
        )
        cb.set_label(
            f"{var_label} / max|{var_label}|@{station_label_norm}\n"
            f"(max = {per_scen[scen]['norm_factor']:.3e}{u_txt})",
            fontsize=8,
        )
        cb.ax.tick_params(labelsize=8)

    win_tag = ""
    if tmin_h is not None or tmax_h is not None:
        a = f"{tmin_h:g}" if tmin_h is not None else "0"
        b = f"{tmax_h:g}" if tmax_h is not None else "end"
        win_tag = f"_t{a}-{b}h"

    depth_title = depth_tag_final if depth_tag_final != "sfc" else "superficie"
    fig.suptitle(
        f"Comparación de escenarios • {var_label} zoom {depth_title} • "
        f"cada {step_frames} pasos • "
        f"norm: max|{var_label}|@{station_label_norm} por escenario",
        fontsize=11, y=0.998,
    )
    fig.subplots_adjust(top=0.94)

    orient_tag = "_T" if transpose else ""
    out_png = (
        save_dir
        / f"compare_{var_tag}_{depth_tag_final}_{station_label_norm}_step{step_frames}{win_tag}{orient_tag}.png"
    )
    if save_fig:
        fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
        print(f"[OK] Figura guardada: {out_png}")
    if show_fig:
        plt.show()
    plt.close(fig)
    return str(out_png)


def graficar_uv_zoom_npz_norm_estacion_c(
    pattern_npz_U: str = "./npz_out/*_U_zoom_z110m.npz",
    pattern_npz_V: str = "./npz_out/*_V_zoom_z110m.npz",
    step_frames: int = 60,
    max_cols: int = 4,
    levels: int = 101,
    station_label_norm: str = "C",
    land_color: str = "#ece2c6",
    cmap_name: str = "RdBu_r",
    sharex: bool = True,
    sharey: bool = True,
    outdir_grid: Path | str | None = None,
    grid_pattern: str | None = None,
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    show_mask_contour: bool = False,
    marker_size: int = 24,
    figsize_panel: tuple[float, float] = (4.1, 3.2),
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    save_dir: str | Path = "./figs",
    dpi: int = 170,
    save_fig: bool = True,
    show_fig: bool = False,
):
    """
    Genera figuras (una por escenario) con paneles de ΔU y ΔV normalizados.

    Layout por escenario:
      - Cada grupo de columnas representa hasta `max_cols` instantes de tiempo.
      - Para cada grupo: subfila superior = ΔU, subfila inferior = ΔV.
      - Colorbar compartida ([-1, 1]) en la columna derecha de cada subgrupofila.
      - Normalización independiente por max(|data|) en la estación `station_label_norm`.

    Parámetros
    ----------
    pattern_npz_U / pattern_npz_V : glob que encuentra los NPZ de U y V.
    step_frames   : cada cuántos índices de tiempo se toma un panel.
    max_cols      : máximo de paneles por subfila (por variable).
    station_label_norm : etiqueta de la estación usada para la normalización.
    """
    # ── 1. Cargar listas de NPZ ──────────────────────────────────────────
    u_list = sorted(glob.glob(pattern_npz_U))
    v_list = sorted(glob.glob(pattern_npz_V))
    if not u_list:
        raise FileNotFoundError(f"No encontré NPZ-U con patrón: {pattern_npz_U}")
    if not v_list:
        raise FileNotFoundError(f"No encontré NPZ-V con patrón: {pattern_npz_V}")
    if step_frames <= 0:
        raise ValueError("step_frames debe ser > 0.")
    if max_cols <= 0:
        raise ValueError("max_cols debe ser > 0.")

    # ── 2. Extraer escenario de cada NPZ ────────────────────────────────
    def _scenario_from_npz(p):
        with np.load(p, allow_pickle=True) as d:
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}
        return str(meta.get("scenario", Path(p).name.split("_")[0]))

    u_dict = {_scenario_from_npz(p): p for p in u_list}
    v_dict = {_scenario_from_npz(p): p for p in v_list}
    scenarios = sorted(set(u_dict) & set(v_dict))
    if not scenarios:
        raise ValueError("No hay escenarios en común entre los NPZ de U y V.")

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    out_paths = []

    for scenario in scenarios:
        u_path = u_dict[scenario]
        v_path = v_dict[scenario]

        # ── 4. Cargar datos ──────────────────────────────────────────────
        with np.load(u_path, allow_pickle=True) as d:
            t_h = d["t_horas"].astype(float)
            Xz = d["X_zoom"].astype(float)
            Yz = d["Y_zoom"].astype(float)
            Mz = d["mask_zoom"].astype(bool)
            AU = d["data"].astype(float)
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}

        with np.load(v_path, allow_pickle=True) as d:
            AV = d["data"].astype(float)

        # Staggered C-grid: U/V may have ny+1 or nx+1; crop to tracer-grid mask shape.
        ny_m, nx_m = Mz.shape
        AU = AU[:, :ny_m, :nx_m]
        AV = AV[:, :ny_m, :nx_m]

        # Convertir coordenadas a km para los ejes
        Xz_km = Xz / 1e3
        Yz_km = Yz / 1e3

        depth_txt = int(round(float(meta.get("depth_m", 110.0))))

        # ── 5. Resolver directorio de grid y estaciones ──────────────────
        grid_dir = None
        if outdir_grid is not None:
            grid_dir = Path(outdir_grid)
        else:
            bay_path = meta.get("paths", {}).get("bay")
            if bay_path:
                grid_dir = Path(bay_path)
            else:
                paths = get_iwbay_paths()
                if scenario in paths["scenarios"]:
                    grid_dir = Path(paths["scenarios"][scenario]) / paths["bay_dir"]

        stations_plot = []
        use_index_coords = False
        if grid_dir is not None:
            gpatt = grid_pattern or meta.get("grid_pattern", "grid.t*.nc")
            gr = _open_mnc_or_raise(grid_dir, str(gpatt))
            try:
                stations = make_stations_from_grid(gr, j_centro=j_centro, radio=radio)
                depth_grid = np.array(gr.variables["Depth"][:], dtype=float)
                ny, nx = depth_grid.shape
                XC = _as_2d_from_grid(gr, ("XC", "xC", "X", "x", "XG", "xG"), axis="x", ny=ny, nx=nx)
                YC = _as_2d_from_grid(gr, ("YC", "yC", "Y", "y", "YG", "yG"), axis="y", ny=ny, nx=nx)

                use_index_coords = _is_index_axis(Xz, nx) and _is_index_axis(Yz, ny)
                xlo, xhi = float(np.min(Xz)), float(np.max(Xz))
                ylo, yhi = float(np.min(Yz)), float(np.max(Yz))
                for st in stations:
                    i, j, lab, col = st
                    col = col or "k"
                    i, j = int(i), int(j)
                    if use_index_coords:
                        xs, ys = float(i), float(j)
                    else:
                        if XC is None or YC is None:
                            continue
                        xs, ys = float(XC[j, i]), float(YC[j, i])
                    if (xlo <= xs <= xhi) and (ylo <= ys <= yhi):
                        stations_plot.append((xs, ys, str(lab), str(col)))
            finally:
                _close_mnc(gr)

        # ── 6. Normalización por estación C ─────────────────────────────
        st_norm = [s for s in stations_plot if s[2] == str(station_label_norm)]
        if st_norm:
            xs_n, ys_n = float(st_norm[0][0]), float(st_norm[0][1])
            ix = int(np.argmin(np.abs(Xz - xs_n)))
            iy = int(np.argmin(np.abs(Yz - ys_n)))
        else:
            ix, iy = AU.shape[2] // 2, AU.shape[1] // 2

        for_U = float(np.nanmax(np.abs(AU[:, iy, ix])))
        for_V = float(np.nanmax(np.abs(AV[:, iy, ix])))
        AU_norm = AU / (for_U if np.isfinite(for_U) and for_U > 0.0 else 1.0)
        AV_norm = AV / (for_V if np.isfinite(for_V) and for_V > 0.0 else 1.0)
        print(
            f"[{scenario}] Normalización est. {station_label_norm}  →  "
            f"max|ΔU| = {for_U:.4e} m/s   |   max|ΔV| = {for_V:.4e} m/s"
        )

        # ── 7. Índices de tiempo ─────────────────────────────────────────
        idxs = np.arange(0, AU_norm.shape[0], step_frames, dtype=int)
        nplots = int(idxs.size)
        if nplots == 0:
            continue

        ncols = min(max_cols, nplots)
        ngroups = int(math.ceil(nplots / ncols))
        # Cada grupo tiene 2 subfilas: [U, V]
        nrows_total = ngroups * 2

        # ── 8. Figura y GridSpec ─────────────────────────────────────────
        fig = plt.figure(
            figsize=(figsize_panel[0] * (ncols + 0.2), figsize_panel[1] * nrows_total)
        )
        gs = fig.add_gridspec(
            nrows=nrows_total,
            ncols=ncols + 1,
            width_ratios=[1.0] * ncols + [0.05],
            wspace=0.15,
            hspace=0.22,
        )

        cmap = plt.get_cmap(cmap_name).copy()
        cmap.set_bad(land_color)
        lvls = np.linspace(-1.0, 1.0, int(levels))

        axes_U = {}   # (group, col) -> Axes
        axes_V = {}
        caxes_U = {}  # group -> cax
        caxes_V = {}
        ax_ref = None

        for g in range(ngroups):
            row_U = g * 2
            row_V = g * 2 + 1
            for c in range(ncols):
                share = ax_ref if (ax_ref is not None and (sharex or sharey)) else None
                ax = fig.add_subplot(
                    gs[row_U, c],
                    sharex=share if sharex else None,
                    sharey=share if sharey else None,
                )
                axes_U[(g, c)] = ax
                if ax_ref is None:
                    ax_ref = ax
                # V comparte directamente con ax_ref (no en cadena via U)
                # para evitar problemas con sharing transitivo en Matplotlib
                ax_v = fig.add_subplot(
                    gs[row_V, c],
                    sharex=ax_ref if sharex else None,
                    sharey=ax_ref if sharey else None,
                )
                axes_V[(g, c)] = ax_v
            caxes_U[g] = fig.add_subplot(gs[row_U, -1])
            caxes_V[g] = fig.add_subplot(gs[row_V, -1])

        # ── 9. Dibujar paneles ───────────────────────────────────────────
        if xlim is not None:
            xlo_vis, xhi_vis = float(xlim[0]) / 1e3, float(xlim[1]) / 1e3
        else:
            xlo_vis, xhi_vis = float(np.min(Xz_km)), float(np.max(Xz_km))
        if ylim is not None:
            ylo_vis, yhi_vis = float(ylim[0]) / 1e3, float(ylim[1]) / 1e3
        else:
            ylo_vis, yhi_vis = float(np.min(Yz_km)), float(np.max(Yz_km))

        row_mappable_U = [None] * ngroups
        row_mappable_V = [None] * ngroups

        for k, tidx in enumerate(idxs):
            g = k // ncols
            c = k % ncols

            for var_label, A_norm, ax, row_map in [
                ("ΔU", AU_norm, axes_U[(g, c)], row_mappable_U),
                ("ΔV", AV_norm, axes_V[(g, c)], row_mappable_V),
            ]:
                Z = np.where(Mz, np.nan, A_norm[tidx])
                mp = ax.contourf(
                    Xz_km, Yz_km, Z,
                    levels=lvls, cmap=cmap, vmin=-1.0, vmax=1.0, extend="both"
                )
                if row_map[g] is None:
                    row_map[g] = mp
                ax.set_facecolor(land_color)

                if show_mask_contour:
                    ax.contour(Xz_km, Yz_km, Mz.astype(float), colors="black",
                               linewidths=0.8, linestyles="--")

                for xs, ys, lab, col in stations_plot:
                    xs_km, ys_km = xs / 1e3, ys / 1e3
                    if not (xlo_vis <= xs_km <= xhi_vis and ylo_vis <= ys_km <= yhi_vis):
                        continue
                    ax.scatter(xs_km, ys_km, s=marker_size, c=col,
                               edgecolors="k", linewidths=0.45, zorder=6)
                    ax.text(xs_km, ys_km, f" {lab}", color=col, fontsize=8,
                            ha="left", va="center", zorder=7, clip_on=True)

                ax.set_title(f"{var_label}  t={t_h[tidx]:.1f} h (idx={int(tidx)})", fontsize=9)
                if c == 0:
                    ax.set_ylabel("Y (km)")

        # Ocultar paneles vacíos
        for k in range(nplots, ngroups * ncols):
            g = k // ncols
            c = k % ncols
            axes_U[(g, c)].axis("off")
            axes_V[(g, c)].axis("off")

        # ── Aplicar límites y ocultar tick labels compartidos ────────────
        # xlim/ylim en ax_ref: se propaga a todos los ejes vía sharex/sharey
        ax_ref.set_xlim(xlo_vis, xhi_vis)
        ax_ref.set_ylim(ylo_vis, yhi_vis)

        for g in range(ngroups):
            is_last_group = (g == ngroups - 1)
            for c in range(ncols):
                ax_u = axes_U[(g, c)]
                ax_v = axes_V[(g, c)]
                # Ocultar labels de X en todas las filas U
                ax_u.tick_params(labelbottom=False)
                ax_u.set_xlabel("")
                # Ocultar labels de X en filas V intermedias (no es el borde inferior)
                if not is_last_group:
                    ax_v.tick_params(labelbottom=False)
                    ax_v.set_xlabel("")
                else:
                    ax_v.set_xlabel("X (km)")
                # Ocultar labels de Y en columnas interiores
                if c > 0:
                    ax_u.tick_params(labelleft=False)
                    ax_u.set_ylabel("")
                    ax_v.tick_params(labelleft=False)
                    ax_v.set_ylabel("")
                else:
                    ax_u.set_ylabel("Y (km)")
                    ax_v.set_ylabel("Y (km)")

        # ── 10. Colorbars ────────────────────────────────────────────────
        for g in range(ngroups):
            for mp, cax, label in [
                (row_mappable_U[g], caxes_U[g],
                 f"ΔU / max|ΔU|@{station_label_norm}\n(max = {for_U:.3e} m/s)"),
                (row_mappable_V[g], caxes_V[g],
                 f"ΔV / max|ΔV|@{station_label_norm}\n(max = {for_V:.3e} m/s)"),
            ]:
                if mp is None:
                    cax.axis("off")
                    continue
                cb = fig.colorbar(mp, cax=cax, ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
                cb.set_label(label, fontsize=8)
                cb.ax.tick_params(labelsize=8)

        fig.suptitle(
            f"{scenario} • ΔU & ΔV zoom z={depth_txt} m • "
            f"cada {step_frames} pasos • "
            f"norm @ {station_label_norm}: max|ΔU|={for_U:.3e} m/s, "
            f"max|ΔV|={for_V:.3e} m/s",
            fontsize=11, y=0.998,
        )
        fig.subplots_adjust(top=0.96)

        out_png = (
            save_dir
            / f"{scenario}_UV_zoom_z{depth_txt}m_norm_{station_label_norm}_step{step_frames}_contourf.png"
        )
        if save_fig:
            fig.savefig(out_png, dpi=dpi)
            out_paths.append(str(out_png))
            print(f"[OK] Figura guardada: {out_png}")
        if show_fig:
            plt.show()
        plt.close(fig)

    return out_paths


def graficar_uv_zoom_npz_norm_mean_abs(
    pattern_npz_U: str = "./npz_out/*_U_zoom_z110m.npz",
    pattern_npz_V: str = "./npz_out/*_V_zoom_z110m.npz",
    step_frames: int = 60,
    max_cols: int = 4,
    levels: int = 101,
    land_color: str = "#ece2c6",
    cmap_name: str = "RdBu_r",
    sharex: bool = True,
    sharey: bool = True,
    outdir_grid: Path | str | None = None,
    grid_pattern: str | None = None,
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    show_mask_contour: bool = False,
    marker_size: int = 24,
    figsize_panel: tuple[float, float] = (4.1, 3.2),
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    save_dir: str | Path = "./figs",
    dpi: int = 170,
    save_fig: bool = True,
    show_fig: bool = False,
):
    """
    Igual que graficar_uv_zoom_npz_norm_estacion_c() pero la normalización usa
    el promedio del valor absoluto de todos los puntos del campo completo
    (mean(|ΔU|) y mean(|ΔV|)) en lugar del máximo absoluto en la estación C.

    Layout por escenario:
      - Cada grupo de columnas representa hasta `max_cols` instantes de tiempo.
      - Para cada grupo: subfila superior = ΔU, subfila inferior = ΔV.
      - Colorbar compartida ([-1, 1]) en la columna derecha de cada subgrupofila.
    """
    # ── 1. Cargar listas de NPZ ──────────────────────────────────────────
    u_list = sorted(glob.glob(pattern_npz_U))
    v_list = sorted(glob.glob(pattern_npz_V))
    if not u_list:
        raise FileNotFoundError(f"No encontré NPZ-U con patrón: {pattern_npz_U}")
    if not v_list:
        raise FileNotFoundError(f"No encontré NPZ-V con patrón: {pattern_npz_V}")
    if step_frames <= 0:
        raise ValueError("step_frames debe ser > 0.")
    if max_cols <= 0:
        raise ValueError("max_cols debe ser > 0.")

    # ── 2. Extraer escenario de cada NPZ ────────────────────────────────
    def _scenario_from_npz(p):
        with np.load(p, allow_pickle=True) as d:
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}
        return str(meta.get("scenario", Path(p).name.split("_")[0]))

    u_dict = {_scenario_from_npz(p): p for p in u_list}
    v_dict = {_scenario_from_npz(p): p for p in v_list}
    scenarios = sorted(set(u_dict) & set(v_dict))
    if not scenarios:
        raise ValueError("No hay escenarios en común entre los NPZ de U y V.")

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    out_paths = []

    for scenario in scenarios:
        u_path = u_dict[scenario]
        v_path = v_dict[scenario]

        # ── 4. Cargar datos ──────────────────────────────────────────────
        with np.load(u_path, allow_pickle=True) as d:
            t_h = d["t_horas"].astype(float)
            Xz = d["X_zoom"].astype(float)
            Yz = d["Y_zoom"].astype(float)
            Mz = d["mask_zoom"].astype(bool)
            AU = d["data"].astype(float)
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}

        with np.load(v_path, allow_pickle=True) as d:
            AV = d["data"].astype(float)

        # Staggered C-grid: U/V may have ny+1 or nx+1; crop to tracer-grid mask shape.
        ny_m, nx_m = Mz.shape
        AU = AU[:, :ny_m, :nx_m]
        AV = AV[:, :ny_m, :nx_m]

        # Convertir coordenadas a km para los ejes
        Xz_km = Xz / 1e3
        Yz_km = Yz / 1e3

        depth_txt = int(round(float(meta.get("depth_m", 110.0))))

        # ── 5. Resolver directorio de grid y estaciones ──────────────────
        grid_dir = None
        if outdir_grid is not None:
            grid_dir = Path(outdir_grid)
        else:
            bay_path = meta.get("paths", {}).get("bay")
            if bay_path:
                grid_dir = Path(bay_path)
            else:
                paths = get_iwbay_paths()
                if scenario in paths["scenarios"]:
                    grid_dir = Path(paths["scenarios"][scenario]) / paths["bay_dir"]

        stations_plot = []
        use_index_coords = False
        if grid_dir is not None:
            gpatt = grid_pattern or meta.get("grid_pattern", "grid.t*.nc")
            gr = _open_mnc_or_raise(grid_dir, str(gpatt))
            try:
                stations = make_stations_from_grid(gr, j_centro=j_centro, radio=radio)
                depth_grid = np.array(gr.variables["Depth"][:], dtype=float)
                ny, nx = depth_grid.shape
                XC = _as_2d_from_grid(gr, ("XC", "xC", "X", "x", "XG", "xG"), axis="x", ny=ny, nx=nx)
                YC = _as_2d_from_grid(gr, ("YC", "yC", "Y", "y", "YG", "yG"), axis="y", ny=ny, nx=nx)

                use_index_coords = _is_index_axis(Xz, nx) and _is_index_axis(Yz, ny)
                xlo, xhi = float(np.min(Xz)), float(np.max(Xz))
                ylo, yhi = float(np.min(Yz)), float(np.max(Yz))
                for st in stations:
                    i, j, lab, col = st
                    col = col or "k"
                    i, j = int(i), int(j)
                    if use_index_coords:
                        xs, ys = float(i), float(j)
                    else:
                        if XC is None or YC is None:
                            continue
                        xs, ys = float(XC[j, i]), float(YC[j, i])
                    if (xlo <= xs <= xhi) and (ylo <= ys <= yhi):
                        stations_plot.append((xs, ys, str(lab), str(col)))
            finally:
                _close_mnc(gr)

        # ── 6. Normalización por promedio absoluto global ────────────────
        # Se enmascaran los puntos de tierra antes de calcular el promedio.
        mask3d = np.broadcast_to(Mz[np.newaxis, :, :], AU.shape)
        AU_ocean = np.where(mask3d, np.nan, AU)
        AV_ocean = np.where(mask3d, np.nan, AV)

        for_U = float(np.nanmean(np.abs(AU_ocean)))
        for_V = float(np.nanmean(np.abs(AV_ocean)))
        AU_norm = AU / (for_U if np.isfinite(for_U) and for_U > 0.0 else 1.0)
        AV_norm = AV / (for_V if np.isfinite(for_V) and for_V > 0.0 else 1.0)
        print(
            f"[{scenario}] Normalización mean|·| global  →  "
            f"mean|ΔU| = {for_U:.4e} m/s   |   mean|ΔV| = {for_V:.4e} m/s"
        )

        # ── 7. Índices de tiempo ─────────────────────────────────────────
        idxs = np.arange(0, AU_norm.shape[0], step_frames, dtype=int)
        nplots = int(idxs.size)
        if nplots == 0:
            continue

        ncols = min(max_cols, nplots)
        ngroups = int(math.ceil(nplots / ncols))
        nrows_total = ngroups * 2

        # ── 8. Figura y GridSpec ─────────────────────────────────────────
        fig = plt.figure(
            figsize=(figsize_panel[0] * (ncols + 0.2), figsize_panel[1] * nrows_total)
        )
        gs = fig.add_gridspec(
            nrows=nrows_total,
            ncols=ncols + 1,
            width_ratios=[1.0] * ncols + [0.05],
            wspace=0.15,
            hspace=0.22,
        )

        cmap = plt.get_cmap(cmap_name).copy()
        cmap.set_bad(land_color)
        lvls = np.linspace(-1.0, 1.0, int(levels))

        axes_U = {}
        axes_V = {}
        caxes_U = {}
        caxes_V = {}
        ax_ref = None

        for g in range(ngroups):
            row_U = g * 2
            row_V = g * 2 + 1
            for c in range(ncols):
                share = ax_ref if (ax_ref is not None and (sharex or sharey)) else None
                ax = fig.add_subplot(
                    gs[row_U, c],
                    sharex=share if sharex else None,
                    sharey=share if sharey else None,
                )
                axes_U[(g, c)] = ax
                if ax_ref is None:
                    ax_ref = ax
                ax_v = fig.add_subplot(
                    gs[row_V, c],
                    sharex=ax_ref if sharex else None,
                    sharey=ax_ref if sharey else None,
                )
                axes_V[(g, c)] = ax_v
            caxes_U[g] = fig.add_subplot(gs[row_U, -1])
            caxes_V[g] = fig.add_subplot(gs[row_V, -1])

        # ── 9. Dibujar paneles ───────────────────────────────────────────
        if xlim is not None:
            xlo_vis, xhi_vis = float(xlim[0]) / 1e3, float(xlim[1]) / 1e3
        else:
            xlo_vis, xhi_vis = float(np.min(Xz_km)), float(np.max(Xz_km))
        if ylim is not None:
            ylo_vis, yhi_vis = float(ylim[0]) / 1e3, float(ylim[1]) / 1e3
        else:
            ylo_vis, yhi_vis = float(np.min(Yz_km)), float(np.max(Yz_km))

        row_mappable_U = [None] * ngroups
        row_mappable_V = [None] * ngroups

        for k, tidx in enumerate(idxs):
            g = k // ncols
            c = k % ncols

            for var_label, A_norm, ax, row_map in [
                ("ΔU", AU_norm, axes_U[(g, c)], row_mappable_U),
                ("ΔV", AV_norm, axes_V[(g, c)], row_mappable_V),
            ]:
                Z = np.where(Mz, np.nan, A_norm[tidx])
                mp = ax.contourf(
                    Xz_km, Yz_km, Z,
                    levels=lvls, cmap=cmap, vmin=-1.0, vmax=1.0, extend="both"
                )
                if row_map[g] is None:
                    row_map[g] = mp
                ax.set_facecolor(land_color)

                if show_mask_contour:
                    ax.contour(Xz_km, Yz_km, Mz.astype(float), colors="black",
                               linewidths=0.8, linestyles="--")

                for xs, ys, lab, col in stations_plot:
                    xs_km, ys_km = xs / 1e3, ys / 1e3
                    if not (xlo_vis <= xs_km <= xhi_vis and ylo_vis <= ys_km <= yhi_vis):
                        continue
                    ax.scatter(xs_km, ys_km, s=marker_size, c=col,
                               edgecolors="k", linewidths=0.45, zorder=6)
                    ax.text(xs_km, ys_km, f" {lab}", color=col, fontsize=8,
                            ha="left", va="center", zorder=7, clip_on=True)

                ax.set_title(f"{var_label}  t={t_h[tidx]:.1f} h (idx={int(tidx)})", fontsize=9)
                if c == 0:
                    ax.set_ylabel("Y (km)")

        # Ocultar paneles vacíos
        for k in range(nplots, ngroups * ncols):
            g = k // ncols
            c = k % ncols
            axes_U[(g, c)].axis("off")
            axes_V[(g, c)].axis("off")

        ax_ref.set_xlim(xlo_vis, xhi_vis)
        ax_ref.set_ylim(ylo_vis, yhi_vis)

        for g in range(ngroups):
            is_last_group = (g == ngroups - 1)
            for c in range(ncols):
                ax_u = axes_U[(g, c)]
                ax_v = axes_V[(g, c)]
                ax_u.tick_params(labelbottom=False)
                ax_u.set_xlabel("")
                if not is_last_group:
                    ax_v.tick_params(labelbottom=False)
                    ax_v.set_xlabel("")
                else:
                    ax_v.set_xlabel("X (km)")
                if c > 0:
                    ax_u.tick_params(labelleft=False)
                    ax_u.set_ylabel("")
                    ax_v.tick_params(labelleft=False)
                    ax_v.set_ylabel("")
                else:
                    ax_u.set_ylabel("Y (km)")
                    ax_v.set_ylabel("Y (km)")

        # ── 10. Colorbars ────────────────────────────────────────────────
        for g in range(ngroups):
            for mp, cax, label in [
                (row_mappable_U[g], caxes_U[g], "ΔU / mean|ΔU|"),
                (row_mappable_V[g], caxes_V[g], "ΔV / mean|ΔV|"),
            ]:
                if mp is None:
                    cax.axis("off")
                    continue
                cb = fig.colorbar(mp, cax=cax, ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
                cb.set_label(label, fontsize=8)
                cb.ax.tick_params(labelsize=8)

        fig.suptitle(
            f"{scenario} • ΔU & ΔV zoom z={depth_txt} m • "
            f"cada {step_frames} pasos • norm mean|·| global",
            fontsize=12, y=0.998,
        )
        fig.subplots_adjust(top=0.96)

        out_png = (
            save_dir
            / f"{scenario}_UV_zoom_z{depth_txt}m_norm_meanabs_step{step_frames}_contourf.png"
        )
        if save_fig:
            fig.savefig(out_png, dpi=dpi)
            out_paths.append(str(out_png))
            print(f"[OK] Figura guardada: {out_png}")
        if show_fig:
            plt.show()
        plt.close(fig)

    return out_paths


def graficar_uv_zoom_npz_quiver(
    pattern_npz_U: str = "./npz_out/*_U_zoom_z110m.npz",
    pattern_npz_V: str = "./npz_out/*_V_zoom_z110m.npz",
    pattern_npz_T: str | None = "./npz_out/*_Temp_zoom_z110m.npz",
    step_frames: int = 60,
    max_cols: int = 4,
    quiver_stride: int = 4,
    scale: float | None = None,
    scale_units: str = "inches",
    bg_speed: bool = False,
    cmap_speed: str = "viridis",
    levels_speed: int = 51,
    cmap_temp: str = "RdBu_r",
    levels_temp: int = 101,
    station_label_norm: str = "C",
    quiverkey_U: float | None = None,
    quiverkey_label: str | None = None,
    y_quiver_max_km: float | None = 0.0,
    speed_clip: float | None = None,
    land_color: str = "#ece2c6",
    arrow_color: str = "k",
    sharex: bool = True,
    sharey: bool = True,
    outdir_grid: Path | str | None = None,
    grid_pattern: str | None = None,
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    show_mask_contour: bool = False,
    marker_size: int = 24,
    figsize_panel: tuple[float, float] = (4.1, 3.2),
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    save_dir: str | Path = "./figs",
    dpi: int = 170,
    save_fig: bool = True,
    show_fig: bool = False,
):
    """
    Genera figuras (una por escenario) con paneles quiver de (ΔU, ΔV).

    Prioridad del fondo:
      1. Temperatura normalizada (pattern_npz_T no es None y hay NPZ disponible).
      2. Magnitud de velocidad (bg_speed=True).
      3. Sin fondo (solo flechas sobre color de tierra).

    Parámetros
    ----------
    pattern_npz_U / V  : glob a los NPZ de U y V.
    pattern_npz_T      : glob a los NPZ de temperatura; None desactiva la capa.
    step_frames        : cada cuántos índices de tiempo se toma un panel.
    max_cols           : máximo de paneles por fila.
    quiver_stride      : submuestreo espacial de las flechas (cada N puntos).
    scale              : parámetro ``scale`` de ax.quiver() — None = auto.
    scale_units        : ``scale_units`` de ax.quiver().
    bg_speed           : si True y no hay temperatura, dibuja |v| de fondo.
    cmap_speed         : colormap para el fondo de velocidad.
    levels_speed       : niveles del contourf de velocidad.
    cmap_temp          : colormap para el fondo de temperatura.
    levels_temp        : niveles del contourf de temperatura.
    station_label_norm : estación para normalizar ΔTemp (igual que en
                         graficar_temp_zoom_npz_norm_estacion_c).
    quiverkey_U        : valor (m/s) de la flecha de referencia; None → p90.
    quiverkey_label    : etiqueta del quiverkey; None → auto.
    arrow_color        : color de las flechas.
    """
    # ── 1. Cargar listas de NPZ ──────────────────────────────────────────
    u_list = sorted(glob.glob(pattern_npz_U))
    v_list = sorted(glob.glob(pattern_npz_V))
    if not u_list:
        raise FileNotFoundError(f"No encontré NPZ-U con patrón: {pattern_npz_U}")
    if not v_list:
        raise FileNotFoundError(f"No encontré NPZ-V con patrón: {pattern_npz_V}")
    if step_frames <= 0:
        raise ValueError("step_frames debe ser > 0.")
    if max_cols <= 0:
        raise ValueError("max_cols debe ser > 0.")

    t_list = sorted(glob.glob(pattern_npz_T)) if pattern_npz_T else []

    # ── 2. Extraer escenario de cada NPZ ────────────────────────────────
    def _scenario_from_npz(p):
        with np.load(p, allow_pickle=True) as d:
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}
        return str(meta.get("scenario", Path(p).name.split("_")[0]))

    u_dict = {_scenario_from_npz(p): p for p in u_list}
    v_dict = {_scenario_from_npz(p): p for p in v_list}
    t_dict = {_scenario_from_npz(p): p for p in t_list}
    scenarios = sorted(set(u_dict) & set(v_dict))
    if not scenarios:
        raise ValueError("No hay escenarios en común entre los NPZ de U y V.")

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    out_paths = []

    for scenario in scenarios:
        u_path = u_dict[scenario]
        v_path = v_dict[scenario]

        # ── 3. Cargar datos U / V ────────────────────────────────────────
        with np.load(u_path, allow_pickle=True) as d:
            t_h = d["t_horas"].astype(float)
            Xz = d["X_zoom"].astype(float)
            Yz = d["Y_zoom"].astype(float)
            Mz = d["mask_zoom"].astype(bool)
            AU = d["data"].astype(float)
            meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}

        with np.load(v_path, allow_pickle=True) as d:
            AV = d["data"].astype(float)

        ny_m, nx_m = Mz.shape
        AU = AU[:, :ny_m, :nx_m]
        AV = AV[:, :ny_m, :nx_m]

        mask3d = np.broadcast_to(Mz[np.newaxis, :, :], AU.shape)
        AU = np.where(mask3d, np.nan, AU)
        AV = np.where(mask3d, np.nan, AV)

        Xz_km = Xz / 1e3
        Yz_km = Yz / 1e3
        depth_txt = int(round(float(meta.get("depth_m", 110.0))))

        # ── 3b. Cargar temperatura (opcional) ───────────────────────────
        AT = None
        if scenario in t_dict:
            with np.load(t_dict[scenario], allow_pickle=True) as d:
                AT = d["data"].astype(float)
            AT = AT[:, :ny_m, :nx_m]
            AT = np.where(mask3d, np.nan, AT)

        # ── 4. Resolver estaciones ───────────────────────────────────────
        grid_dir = None
        if outdir_grid is not None:
            grid_dir = Path(outdir_grid)
        else:
            bay_path = meta.get("paths", {}).get("bay")
            if bay_path:
                grid_dir = Path(bay_path)
            else:
                paths = get_iwbay_paths()
                if scenario in paths["scenarios"]:
                    grid_dir = Path(paths["scenarios"][scenario]) / paths["bay_dir"]

        stations_plot = []
        if grid_dir is not None:
            gpatt = grid_pattern or meta.get("grid_pattern", "grid.t*.nc")
            gr = _open_mnc_or_raise(grid_dir, str(gpatt))
            try:
                stations = make_stations_from_grid(gr, j_centro=j_centro, radio=radio)
                depth_grid = np.array(gr.variables["Depth"][:], dtype=float)
                ny, nx = depth_grid.shape
                XC = _as_2d_from_grid(gr, ("XC", "xC", "X", "x", "XG", "xG"), axis="x", ny=ny, nx=nx)
                YC = _as_2d_from_grid(gr, ("YC", "yC", "Y", "y", "YG", "yG"), axis="y", ny=ny, nx=nx)
                use_index_coords = _is_index_axis(Xz, nx) and _is_index_axis(Yz, ny)
                xlo, xhi = float(np.min(Xz)), float(np.max(Xz))
                ylo, yhi = float(np.min(Yz)), float(np.max(Yz))
                for st in stations:
                    i, j, lab, col = st
                    col = col or "k"
                    i, j = int(i), int(j)
                    if use_index_coords:
                        xs, ys = float(i), float(j)
                    else:
                        if XC is None or YC is None:
                            continue
                        xs, ys = float(XC[j, i]), float(YC[j, i])
                    if (xlo <= xs <= xhi) and (ylo <= ys <= yhi):
                        stations_plot.append((xs, ys, str(lab), str(col)))
            finally:
                _close_mnc(gr)

        # ── 5. Preparar submuestreo para quiver ──────────────────────────
        s = int(quiver_stride)
        Xq_km = Xz_km[::s]
        Yq_km = Yz_km[::s]

        speed_all = np.sqrt(AU**2 + AV**2)
        ref_U = quiverkey_U
        if ref_U is None:
            ref_U = float(np.nanpercentile(speed_all, 90))
        if not np.isfinite(ref_U) or ref_U <= 0.0:
            ref_U = 0.01
        ref_label = quiverkey_label or f"{ref_U:.3g} m/s"

        # Rangos de velocidad para fondo (usados solo si bg_speed=True y AT es None)
        vmax_speed = float(np.nanpercentile(speed_all, 99))
        if not np.isfinite(vmax_speed) or vmax_speed <= 0.0:
            vmax_speed = 1.0
        lvls_speed = np.linspace(0.0, vmax_speed, int(levels_speed))

        # Normalización de temperatura por máximo en estación `station_label_norm`
        temp_norm_factor = 1.0
        lvls_temp = np.linspace(-1.0, 1.0, int(levels_temp))
        if AT is not None:
            st_norm = [s_ for s_ in stations_plot if s_[2] == str(station_label_norm)]
            if st_norm:
                xs_n, ys_n = float(st_norm[0][0]), float(st_norm[0][1])
                ix_n = int(np.argmin(np.abs(Xz - xs_n)))
                iy_n = int(np.argmin(np.abs(Yz - ys_n)))
            else:
                ix_n, iy_n = AT.shape[2] // 2, AT.shape[1] // 2
            temp_norm_factor = float(np.nanmax(np.abs(AT[:, iy_n, ix_n])))
            if not np.isfinite(temp_norm_factor) or temp_norm_factor <= 0.0:
                temp_norm_factor = 1.0
            AT_norm = AT / temp_norm_factor
            print(
                f"[{scenario}] quiver+Temp  →  ref_v = {ref_U:.4e} m/s  |  "
                f"norm_T (est. {station_label_norm}) = {temp_norm_factor:.4e} °C"
            )
        else:
            AT_norm = None
            print(
                f"[{scenario}] quiver  →  ref = {ref_U:.4e} m/s  |  "
                f"vmax_speed(p99) = {vmax_speed:.4e} m/s"
            )

        # Decidir qué se usa como fondo
        use_temp_bg = AT_norm is not None
        use_speed_bg = (not use_temp_bg) and bg_speed

        # ── 6. Índices de tiempo ─────────────────────────────────────────
        idxs = np.arange(0, AU.shape[0], step_frames, dtype=int)
        nplots = int(idxs.size)
        if nplots == 0:
            continue

        ncols = min(max_cols, nplots)
        nrows = int(math.ceil(nplots / ncols))

        # ── 7. Figura ────────────────────────────────────────────────────
        fig = plt.figure(
            figsize=(figsize_panel[0] * (ncols + 0.2), figsize_panel[1] * nrows)
        )
        gs = fig.add_gridspec(
            nrows=nrows,
            ncols=ncols + 1,
            width_ratios=[1.0] * ncols + [0.05],
            wspace=0.15,
            hspace=0.22,
        )
        axes = np.empty((nrows, ncols), dtype=object)
        caxes = np.empty((nrows,), dtype=object)
        ax_ref = None
        for rr in range(nrows):
            for cc in range(ncols):
                share_ax = ax_ref if (ax_ref is not None and (sharex or sharey)) else None
                ax = fig.add_subplot(
                    gs[rr, cc],
                    sharex=share_ax if sharex else None,
                    sharey=share_ax if sharey else None,
                )
                axes[rr, cc] = ax
                if ax_ref is None:
                    ax_ref = ax
            caxes[rr] = fig.add_subplot(gs[rr, -1])

        cmap_sp = plt.get_cmap(cmap_speed).copy()
        cmap_sp.set_bad(land_color)
        cmap_tp = plt.get_cmap(cmap_temp).copy()
        cmap_tp.set_bad(land_color)

        if xlim is not None:
            xlo_vis, xhi_vis = float(xlim[0]) / 1e3, float(xlim[1]) / 1e3
        else:
            xlo_vis, xhi_vis = float(np.min(Xz_km)), float(np.max(Xz_km))
        if ylim is not None:
            ylo_vis, yhi_vis = float(ylim[0]) / 1e3, float(ylim[1]) / 1e3
        else:
            ylo_vis, yhi_vis = float(np.min(Yz_km)), float(np.max(Yz_km))

        row_mappables = [None] * nrows

        # ── 8. Dibujar paneles ───────────────────────────────────────────
        for k, tidx in enumerate(idxs):
            r = k // ncols
            c = k % ncols
            ax = axes[r, c]
            ax.set_facecolor(land_color)

            Ut = AU[tidx]
            Vt = AV[tidx]

            # Fondo: temperatura (prioridad) o velocidad
            if use_temp_bg:
                mp = ax.contourf(
                    Xz_km, Yz_km, AT_norm[tidx],
                    levels=lvls_temp, cmap=cmap_tp,
                    vmin=-1.0, vmax=1.0, extend="both"
                )
                if row_mappables[r] is None:
                    row_mappables[r] = mp
            elif use_speed_bg:
                spd = np.sqrt(Ut**2 + Vt**2)
                mp = ax.contourf(
                    Xz_km, Yz_km, spd,
                    levels=lvls_speed, cmap=cmap_sp,
                    vmin=0.0, vmax=vmax_speed, extend="max"
                )
                if row_mappables[r] is None:
                    row_mappables[r] = mp

            # Quiver submuestreado
            Uq = Ut[::s, ::s]
            Vq = Vt[::s, ::s]
            Xq2, Yq2 = np.meshgrid(Xq_km, Yq_km)

            # Máscara: tierra + límite Y
            valid = np.isfinite(Uq) & np.isfinite(Vq)
            if y_quiver_max_km is not None:
                valid &= (Yq2 <= float(y_quiver_max_km))

            # Clamp de magnitud para no dibujar flechas gigantes
            if speed_clip is not None and float(speed_clip) > 0.0:
                spd_q = np.sqrt(Uq**2 + Vq**2)
                ratio = np.where(spd_q > float(speed_clip),
                                 float(speed_clip) / np.where(spd_q > 0, spd_q, 1.0),
                                 1.0)
                Uq = Uq * ratio
                Vq = Vq * ratio

            qv_kwargs = dict(scale_units=scale_units, angles="xy")
            if scale is not None:
                qv_kwargs["scale"] = scale
            ax.quiver(
                Xq2[valid], Yq2[valid],
                Uq[valid], Vq[valid],
                color=arrow_color, **qv_kwargs
            )

            if show_mask_contour:
                ax.contour(Xz_km, Yz_km, Mz.astype(float), colors="black",
                           linewidths=0.8, linestyles="--")

            for xs, ys, lab, col in stations_plot:
                xs_km, ys_km = xs / 1e3, ys / 1e3
                if not (xlo_vis <= xs_km <= xhi_vis and ylo_vis <= ys_km <= yhi_vis):
                    continue
                ax.scatter(xs_km, ys_km, s=marker_size, c=col,
                           edgecolors="k", linewidths=0.45, zorder=6)
                ax.text(xs_km, ys_km, f" {lab}", color=col, fontsize=8,
                        ha="left", va="center", zorder=7, clip_on=True)

            ax.set_title(f"t={t_h[tidx]:.1f} h (idx={int(tidx)})", fontsize=9)

        # Ocultar paneles vacíos
        for k in range(nplots, nrows * ncols):
            r = k // ncols
            c = k % ncols
            axes[r, c].axis("off")

        ax_ref.set_xlim(xlo_vis, xhi_vis)
        ax_ref.set_ylim(ylo_vis, yhi_vis)

        # Etiquetas y tick visibility
        for r in range(nrows):
            for c in range(ncols):
                ax = axes[r, c]
                if r < nrows - 1:
                    ax.tick_params(labelbottom=False)
                    ax.set_xlabel("")
                else:
                    ax.set_xlabel("X (km)")
                if c == 0:
                    ax.set_ylabel("Y (km)")
                else:
                    ax.tick_params(labelleft=False)
                    ax.set_ylabel("")

        # ── 9. Colorbars y quiverkey ─────────────────────────────────────
        for r in range(nrows):
            cax = caxes[r]
            if row_mappables[r] is not None:
                cb = fig.colorbar(row_mappables[r], cax=cax)
                if use_temp_bg:
                    cb.set_label(f"ΔTemp norm. (est. {station_label_norm})", fontsize=8)
                    cb.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
                else:
                    cb.set_label("|v| (m/s)", fontsize=8)
                cb.ax.tick_params(labelsize=8)
            else:
                cax.axis("off")

        # Quiverkey en el primer panel activo
        ax0 = axes[0, 0]
        qv_ref_kwargs = dict(scale_units=scale_units, angles="xy", color=arrow_color)
        if scale is not None:
            qv_ref_kwargs["scale"] = scale
        x_key = xlo_vis + 0.75 * (xhi_vis - xlo_vis)
        y_key = yhi_vis - 0.05 * (yhi_vis - ylo_vis)
        ax0.quiverkey(
            ax0.quiver(
                [x_key], [y_key], [ref_U], [0.0],
                **qv_ref_kwargs
            ),
            X=0.78, Y=0.96,
            U=ref_U,
            label=ref_label,
            labelpos="E",
            fontproperties={"size": 7},
            coordinates="axes",
        )

        bg_tag = "T" if use_temp_bg else ("spd" if use_speed_bg else "noFondo")
        fig.suptitle(
            f"{scenario} • (ΔU, ΔV) quiver z={depth_txt} m • "
            f"cada {step_frames} pasos • stride={quiver_stride} • fondo={bg_tag}",
            fontsize=12, y=0.998,
        )
        fig.subplots_adjust(top=0.96)

        out_png = (
            save_dir
            / f"{scenario}_UV_zoom_z{depth_txt}m_quiver_{bg_tag}_step{step_frames}_s{quiver_stride}.png"
        )
        if save_fig:
            fig.savefig(out_png, dpi=dpi)
            out_paths.append(str(out_png))
            print(f"[OK] Figura guardada: {out_png}")
        if show_fig:
            plt.show()
        plt.close(fig)

    return out_paths


### ─────────────────────────────────────────────────────────────

# ───────────────── helpers usados (se asume ya existen en tu módulo) ───────────────
# _open_mnc_or_raise(dirpath, pattern)   -> MNC
# _close_mnc(mnc_obj)                    -> None
# _get_time_from_nc(nc)                  -> np.ndarray (t,)

# ───────────────── coords y máscara robustos ───────────────────────────────────────
def _coords_zoom_from_grid(gr, X_sym: int, Y_cut: int):
    """Devuelve (X_zoom, Y_zoom, mask_zoom, slices, ny, nx) usando XC/XG, YC/YG, Depth."""
    Depth = np.array(gr.variables['Depth'][:], dtype=float)  # (ny,nx)
    ny, nx = Depth.shape

    XC = _to_1d_x(gr, nx, 'XC', 'XG')   # (nx,)
    YC = _to_1d_y(gr, ny, 'YC', 'YG')   # (ny,)

    x_center = nx // 2
    x0 = max(0, x_center - X_sym); x1 = min(nx, x_center + X_sym)
    y0 = max(0, ny - Y_cut);       y1 = ny

    x_slice = slice(x0, x1)
    y_slice = slice(y0, y1)

    X_zoom = XC[x_slice]
    Y_zoom = YC[y_slice]
    mask_zoom = np.ma.getmask(np.ma.masked_values(Depth, 0.0))[y_slice, x_slice]

    return X_zoom, Y_zoom, mask_zoom, (y_slice, x_slice), ny, nx

# ───────────────── guardar ETA (2D) en zoom (streaming) ────────────────────────────
def guardar_zoom_eta_npz(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    subtract: bool = True,           # bay - nobay
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern   = "grid.t*.nc",
    X_sym: int = 150,
    Y_cut: int = 420
) -> str:
    """
    Guarda un bundle .npz con Eta (o ΔEta) en una región de zoom.
    Salida: <scenario_label>_Eta_zoom_sfc.npz
    """
    save_dir.mkdir(parents=True, exist_ok=True)

    # Grid (una sola vez)
    gr = _open_mnc_or_raise(outdir_bay, grid_pattern)
    try:
        X_zoom, Y_zoom, mask_zoom, (y_slice, x_slice), _, _ = _coords_zoom_from_grid(gr, X_sym, Y_cut)

        # Busca primera etapa válida (ambas carpetas) para verificar shape 2D
        first_nc_b = None; first_nc_n = None
        for patt in stage_patterns:
            try:
                first_nc_b = _open_mnc_or_raise(outdir_bay, patt)
                first_nc_n = _open_mnc_or_raise(outdir_nobay, patt)
                break
            except FileNotFoundError:
                continue
        if first_nc_b is None:
            raise FileNotFoundError("No se encontró ninguna etapa en bay/nobay.")

        try:
            vb0 = first_nc_b.variables['Eta']
            shp0 = tuple(getattr(vb0, "shape", ()))
            if len(shp0) != 3:
                raise ValueError(f"'Eta' debería ser (t,y,x). Encontrado: {shp0}")
            # tiempo ref
            t0 = _get_time_from_nc(first_nc_b)
        finally:
            _close_mnc(first_nc_b); _close_mnc(first_nc_n)

        # Streaming por etapa
        t_list = []
        data_list = []
        for patt in stage_patterns:
            try:
                nc_b = _open_mnc_or_raise(outdir_bay, patt)
                nc_n = _open_mnc_or_raise(outdir_nobay, patt)
            except FileNotFoundError:
                continue
            try:
                tb = _get_time_from_nc(nc_b)
                tn = _get_time_from_nc(nc_n)

                vb = nc_b.variables['Eta']
                vn = nc_n.variables['Eta']

                Fb = np.array(vb[:, y_slice, x_slice], dtype=float)  # (t, y, x)
                Fn = np.array(vn[:, y_slice, x_slice], dtype=float)
                if Fb.shape[0] != tb.shape[0] or Fn.shape[0] != tn.shape[0]:
                    raise ValueError("Longitudes tiempo/serie no coinciden en una etapa.")

                Dk = Fb - Fn if subtract else Fb

                if t_list and abs(tb[0] - t_list[-1][-1]) <= 1e-9:
                    tb = tb[1:]; Dk = Dk[1:, ...]
                t_list.append(tb); data_list.append(Dk)
            finally:
                _close_mnc(nc_b); _close_mnc(nc_n)

        if not t_list:
            raise RuntimeError("No se pudo construir la serie temporal de Eta.")

        t_full = np.concatenate(t_list)            # (nt_full,)
        D_full = np.concatenate(data_list, axis=0) # (nt_full, y, x)

        # recorte por spin-up
        tmin_seg = tmin_horas * 3600.0
        idx0 = np.where(t_full >= tmin_seg)[0][0] if np.any(t_full >= tmin_seg) else 0
        t_horas = (t_full[idx0:] / 3600.0).astype(np.float64)
        D = D_full[idx0:, :, :].astype(np.float32)

        out_path = save_dir / f"{scenario_label}_Eta_zoom_sfc.npz"
        meta = {
            "scenario": scenario_label,
            "var": "Eta",
            "subtract": bool(subtract),
            "tmin_horas": float(tmin_horas),
            "stage_patterns": list(stage_patterns),
            "zoom": {"X_sym": int(X_sym), "Y_cut": int(Y_cut)},
            "grid_pattern": grid_pattern
        }
        meta_json = json.dumps(meta, ensure_ascii=False)

        np.savez_compressed(
            out_path,
            t_horas=t_horas,
            X_zoom=X_zoom.astype(np.float64),
            Y_zoom=Y_zoom.astype(np.float64),
            mask_zoom=mask_zoom.astype(np.uint8),
            data=D,
            meta_json=np.array(meta_json)
        )
        print(f"[OK] Guardado bundle zoom (Eta): {out_path}")
        return str(out_path)
    finally:
        _close_mnc(gr)

# ───────────────── animación desde el NPZ de Eta ───────────────────────────────────
def animar_zoom_eta_desde_npz(
    npz_path: str,
    outname: str = "anim_eta_zoom.mp4",
    vmin: float | None = None,
    vmax: float | None = None,
    fps: int = 3,
    cmap = cmo.cm.balance,
    show_mask_contour: bool = True
):
    """
    Lee el NPZ de Eta zoom y genera una animación MP4.
    Por defecto, si subtract=True (guardado), usa escala simétrica.
    """
    d = np.load(npz_path, allow_pickle=True)
    t_h = d["t_horas"]            # (nt,)
    Xz  = d["X_zoom"]             # (nxz,)
    Yz  = d["Y_zoom"]             # (nyz,)
    Mz  = d["mask_zoom"]          # (nyz,nxz)
    A   = d["data"]               # (nt,nyz,nxz)
    meta= json.loads(str(d["meta_json"]))
    subtract = bool(meta.get("subtract", True))

    # vmin/vmax
    if (vmin is None) or (vmax is None):
        if subtract:
            m = float(np.nanmax(np.abs(A)))
            vmin, vmax = -m, m
        else:
            vmin = float(np.nanmin(A)); vmax = float(np.nanmax(A))
            if vmin == vmax:
                m = abs(vmin) if vmin != 0 else 1e-9
                vmin, vmax = -m, m

    fig, ax = plt.subplots(figsize=(8,6))
    mesh = ax.pcolormesh(Xz, Yz, A[0], cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    cbar = fig.colorbar(mesh, ax=ax, label=f"Eta{' (bay−nobay)' if subtract else ' (bay)'}")
    if show_mask_contour:
        ax.contour(Xz, Yz, Mz, colors='black', linewidths=1.0, linestyles='--')
    ax.set_aspect('equal')
    title_base = f"{meta.get('scenario','')} • Eta • zoom"
    ax.set_title(f"{title_base} • t={t_h[0]:.2f} h")
    plt.tight_layout()

    def update(frame):
        mesh.set_array(A[frame].ravel())
        ax.set_title(f"{title_base} • t={t_h[frame]:.2f} h")
        return (mesh,)

    anim = FuncAnimation(fig, update, frames=A.shape[0], blit=True)
    anim.save(outname, writer="ffmpeg", fps=fps)
    plt.close(fig)
    print(f"[OK] Animación exportada: {outname}")
    return outname




### ----------------------------------------------------------------------------------
### VISUALIZAR VARIABLES Y DIMENSIONES

def _get_attr_safe(var, keys=("long_name","units","standard_name","description")):
    """Intenta leer atributos comunes de forma robusta para xarray, netCDF4 y MNC."""
    # xarray.DataArray
    if hasattr(var, "attrs"):
        return (
            var.attrs.get("long_name",""),
            var.attrs.get("units",""),
            var.attrs.get("standard_name",""),
            var.attrs.get("description",""),
        )
    # netCDF4.Variable
    if hasattr(var, "ncattrs"):
        def g(k):
            try: return var.getncattr(k)
            except Exception: return ""
        return (g("long_name"), g("units"), g("standard_name"), g("description"))
    # MNC variable (intenta dicts internos comunes)
    for attrdict_name in ("_attributes","attributes","__dict__"):
        ad = getattr(var, attrdict_name, None)
        if isinstance(ad, dict):
            return (
                ad.get("long_name",""),
                ad.get("units",""),
                ad.get("standard_name",""),
                ad.get("description",""),
            )
    return ("","","","")

def resumen_cualquier_dataset(ds, nombre="dataset"):
    """
    Devuelve un DataFrame con Variable, Dimensiones, Forma, Dtype y Atributos
    para:
      - xarray.Dataset
      - MITgcmutils MNC (mnc_files)
      - netCDF4.Dataset (si lo usas)
    """
    filas = []

    # Caso xarray.Dataset
    if hasattr(ds, "data_vars"):
        for varname, da in ds.data_vars.items():
            ln, un, stdn, desc = _get_attr_safe(da)
            filas.append({
                "Variable": varname,
                "Dimensiones": " x ".join(da.dims),
                "Forma": str(tuple(da.shape)),
                "Tipo de dato": str(da.dtype),
                "long_name": ln, "units": un, "standard_name": stdn, "description": desc
            })
        return pd.DataFrame(filas)

    # Caso MNC (MITgcmutils) o netCDF4-like: tienen .variables (dict-like)
    if hasattr(ds, "variables"):
        for varname, v in ds.variables.items():
            # dims
            dims = getattr(v, "dimensions", ())
            try:
                shp = getattr(v, "shape", None)
                shp = tuple(shp) if shp is not None else tuple()
            except Exception:
                # última opción: pide tamaño mínimo (costo I/O bajo)
                try:
                    shp = tuple(v[tuple(slice(0,1) for _ in range(len(dims)))].shape)
                except Exception:
                    shp = tuple()

            # dtype
            dt = getattr(v, "dtype", None)
            if dt is None:
                try:
                    dt = getattr(getattr(v, "_data", None), "dtype", "")  # por si existe
                except Exception:
                    dt = ""

            ln, un, stdn, desc = _get_attr_safe(v)
            filas.append({
                "Variable": varname,
                "Dimensiones": " x ".join(dims),
                "Forma": str(shp),
                "Tipo de dato": str(dt),
                "long_name": ln, "units": un, "standard_name": stdn, "description": desc
            })
        return pd.DataFrame(filas)

    raise TypeError(f"No reconozco el tipo de objeto para {nombre}: {type(ds)}")

# --------------------------------------------------------------------
# θ(t,z) full perfil en todas las estaciones
# --------------------------------------------------------------------
def calcular_y_guardar_theta_perfil_escenario(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    source: str = "diff",                 # "diff" (bay-nobay), "bay", "nobay"
    angle_from: str = "vertical",         # "vertical" => atan2(U,W) ; "horizontal" => atan2(W,U)
    collocate_to_tracer: bool = True,     # U_tr = 0.5*(U[i]+U[i-1]) si i>0; fallback si i==0
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern   = "grid.t*.nc",
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    incluir_gamma_omega: bool = True
) -> str:
    """
    Calcula θ(t,z) para cada estación del escenario y guarda un NPZ con:
      - t_horas (nt,), RC (nz,), theta_deg (nt,nz,nst), labels (nst,)
      - meta_json
    Definición por defecto: θ respecto a vertical = atan2(U, W).
    """
    save_dir.mkdir(parents=True, exist_ok=True)

    # GRID una sola vez
    gr = _open_mnc_or_raise(outdir_bay, grid_pattern)
    try:
        RC = _get_RC(gr)
        nz = RC.size

        # Estaciones (usa tu función ampliada)
        try:
            from MIT_DatProcesW import make_stations_from_grid as _mk
            stations_save, stations_meta, _ = _mk(
                gr, j_centro=j_centro, radio=radio, incluir_gamma_omega=incluir_gamma_omega
            )
        except Exception:
            # Fallback mínimo
            ny, nx = np.array(gr.variables['Depth'][:]).shape
            i_c = nx//2
            clamp = lambda v, lo, hi: max(lo, min(hi, v))
            j0 = clamp(j_centro, 0, ny-1)
            iS,jS = i_c, clamp(j0+20, 0, ny-1)
            iA,jA = clamp(i_c-radio,0,nx-1), j0
            iB,jB = clamp(i_c+radio,0,nx-1), j0
            iC,jC = i_c, clamp(j0-radio,0,ny-1)
            st = [(iS,jS,'Seiche'),(iA,jA,'A'),(iB,jB,'B'),(iC,jC,'C')]
            if incluir_gamma_omega:
                iG,jG = i_c, clamp(j0-2*radio,0,ny-1)
                iO,jO = i_c, j0
                st += [(iG,jG,'Gamma'),(iO,jO,'Omega')]
            stations_save = st
            stations_meta = [{"label": lab, "i": i, "j": j, "color": None} for (i,j,lab) in st]

        labels = [d["label"] if isinstance(d, dict) else lab for d,(i,j,lab) in zip(stations_meta, stations_save)]
        nst = len(stations_save)

        # Detectar rank con primera etapa
        first_nc_b = first_nc_n = None
        for patt in stage_patterns:
            try:
                first_nc_b = _open_mnc_or_raise(outdir_bay, patt)
                first_nc_n = _open_mnc_or_raise(outdir_nobay, patt)
                break
            except FileNotFoundError:
                continue
        if first_nc_b is None:
            raise FileNotFoundError("No encontré ninguna etapa en bay/nobay.")

        try:
            Ub0 = first_nc_b.variables['U']
            Wb0 = first_nc_b.variables['W']
            shpU = tuple(getattr(Ub0, "shape", ()))
            shpW = tuple(getattr(Wb0, "shape", ()))
            if len(shpU) != 4 or len(shpW) != 4:
                raise ValueError(f"Se esperan U,W 3D (t,z,y,x). U:{shpU} W:{shpW}")
            t0 = _get_time_from_nc(first_nc_b)
        finally:
            _close_mnc(first_nc_b); _close_mnc(first_nc_n)

        # Acumuladores por estación: listas de bloques (t_stage, nz)
        t_lists  = [[] for _ in range(nst)]
        U_lists  = [[] for _ in range(nst)]
        W_lists  = [[] for _ in range(nst)]

        # Recorremos etapas en streaming
        for patt in stage_patterns:
            try:
                nc_b = _open_mnc_or_raise(outdir_bay, patt)
                nc_n = _open_mnc_or_raise(outdir_nobay, patt)
            except FileNotFoundError:
                continue
            try:
                tb = _get_time_from_nc(nc_b)
                tn = _get_time_from_nc(nc_n)
                Ub = nc_b.variables['U']; Wb = nc_b.variables['W']
                Un = nc_n.variables['U']; Wn = nc_n.variables['W']

                for sidx, (i, j, lab) in enumerate(stations_save):
                    # U colocalizado (z,t) en columna (j,i)
                    if collocate_to_tracer:
                        if i > 0:
                            u_b = 0.5*(np.array(Ub[:, :, j, i], dtype=float) + np.array(Ub[:, :, j, i-1], dtype=float))
                            u_n = 0.5*(np.array(Un[:, :, j, i], dtype=float) + np.array(Un[:, :, j, i-1], dtype=float))
                        else:
                            # si i==0 y no hay i-1, usa promedio i e i+1 si existe; si no, solo i
                            try:
                                u_b = 0.5*(np.array(Ub[:, :, j, i], dtype=float) + np.array(Ub[:, :, j, i+1], dtype=float))
                                u_n = 0.5*(np.array(Un[:, :, j, i], dtype=float) + np.array(Un[:, :, j, i+1], dtype=float))
                            except Exception:
                                u_b = np.array(Ub[:, :, j, i], dtype=float)
                                u_n = np.array(Un[:, :, j, i], dtype=float)
                    else:
                        u_b = np.array(Ub[:, :, j, i], dtype=float)
                        u_n = np.array(Un[:, :, j, i], dtype=float)

                    # W ya está en puntos tracer (z,t) en (j,i)
                    w_b = np.array(Wb[:, :, j, i], dtype=float)
                    w_n = np.array(Wn[:, :, j, i], dtype=float)

                    if source == "diff":
                        u_k = u_b - u_n
                        w_k = w_b - w_n
                    elif source == "bay":
                        u_k = u_b; w_k = w_b
                    elif source == "nobay":
                        u_k = u_n; w_k = w_n
                    else:
                        raise ValueError("source debe ser 'diff', 'bay' o 'nobay'.")

                    # Evitar duplicar primer sample de la etapa
                    if t_lists[sidx]:
                        if abs(tb[0] - t_lists[sidx][-1][-1]) <= 1e-9:
                            tb_use = tb[1:]; u_k = u_k[1:,:]; w_k = w_k[1:,:]
                        else:
                            tb_use = tb
                    else:
                        tb_use = tb

                    t_lists[sidx].append(tb_use)
                    U_lists[sidx].append(u_k)   # (t_stage, nz)
                    W_lists[sidx].append(w_k)

            finally:
                _close_mnc(nc_b); _close_mnc(nc_n)

        # Construir eje de tiempo común (usamos el de la estación 0 tras recorte)
        t_full0 = np.concatenate(t_lists[0]) if t_lists[0] else np.array([])
        if t_full0.size == 0:
            raise RuntimeError("No se pudo construir el tiempo concatenado.")

        tmin_seg = tmin_horas * 3600.0
        idx0 = np.where(t_full0 >= tmin_seg)[0][0] if np.any(t_full0 >= tmin_seg) else 0
        t_horas = (t_full0[idx0:] / 3600.0).astype(np.float64)
        nt = t_horas.size

        # Salida θ(t,z,station)
        theta_deg = np.full((nt, nz, nst), np.nan, dtype=np.float32)
        eps = 1e-12

        for sidx in range(nst):
            t_full = np.concatenate(t_lists[sidx]) if t_lists[sidx] else np.array([])
            U_full = np.concatenate(U_lists[sidx], axis=0) if U_lists[sidx] else np.empty((0, nz))
            W_full = np.concatenate(W_lists[sidx], axis=0) if W_lists[sidx] else np.empty((0, nz))
            if t_full.size == 0: continue

            Uc = U_full[idx0:, :]    # (nt, nz)
            Wc = W_full[idx0:, :]    # (nt, nz)

            if angle_from == "vertical":
                # θ respecto a vertical: atan2(U, W)
                Th = np.degrees(np.arctan2(Uc, np.where(np.abs(Wc) < eps, np.sign(Wc)*eps, Wc)))
            elif angle_from == "horizontal":
                # θ respecto a horizontal: atan2(W, U)
                Th = np.degrees(np.arctan2(Wc, np.where(np.abs(Uc) < eps, np.sign(Uc)*eps, Uc)))
            else:
                raise ValueError("angle_from debe ser 'vertical' o 'horizontal'.")

            # opcional: envolver a (-180,180]
            Th = ((Th + 180.0) % 360.0) - 180.0
            theta_deg[:, :, sidx] = Th.astype(np.float32)

        # Guardar
        meta = {
            "scenario": scenario_label,
            "source": source,
            "angle_from": angle_from,
            "collocate_to_tracer": bool(collocate_to_tracer),
            "tmin_horas": float(tmin_horas),
            "stage_patterns": list(stage_patterns),
            "grid_pattern": grid_pattern,
            "stations": stations_meta
        }
        meta_json = json.dumps(meta, ensure_ascii=False)
        out_path = save_dir / f"{scenario_label}_thetaProfile_{source}_{angle_from}.npz"

        np.savez_compressed(
            out_path,
            t_horas=t_horas,
            RC=RC.astype(np.float64),
            theta_deg=theta_deg,                      # (nt, nz, nst)
            labels=np.array(labels, dtype="U16"),
            meta_json=np.array(meta_json)
        )
        print(f"[OK] θ(t,z) guardado: {out_path}")
        return str(out_path)

    finally:
        _close_mnc(gr)

# --------------------------------------------------------------------
# Plot Hovmöller de θ para una estación (desde el bundle)
# --------------------------------------------------------------------
def plot_theta_hov_from_npz(
    npz_path: str,
    station: str | int = "A",
    vmin: float | None = None,
    vmax: float | None = None,
    cmap = cmo.cm.phase,
    title_prefix: str = ""
):
    d = np.load(npz_path, allow_pickle=True)
    t   = d["t_horas"]       # (nt,)
    RC  = d["RC"]            # (nz,)
    TH3 = d["theta_deg"]     # (nt,nz,nst)
    labs= [str(x) for x in d["labels"]]
    meta= json.loads(str(d["meta_json"]))

    # Selección de estación
    if isinstance(station, int):
        k = station
        label = labs[k]
    else:
        label = str(station)
        if label not in labs:
            raise KeyError(f"Estación '{label}' no está en {labs}")
        k = labs.index(label)

    Z = TH3[:, :, k]      # (nt, nz)

    # Escala: simétrica por defecto
    if (vmin is None) or (vmax is None):
        m = float(np.nanmax(np.abs(Z)))
        vmin, vmax = -m, m
        if m == 0:
            vmin, vmax = -1.0, 1.0

    fig, ax = plt.subplots(figsize=(9,5))
    pcm = ax.pcolormesh(t, RC, Z.T, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    fig.colorbar(pcm, ax=ax, label=f"θ (°) [{meta.get('angle_from','vertical')}]")
    ax.set_xlabel("Tiempo (h)")
    ax.set_ylabel("Profundidad (m)")
    ax.set_title(f"{title_prefix} θ(t,z) — {meta.get('scenario','')} • est. {label}")
    plt.tight_layout()
    plt.show()
    return fig, ax

# --------------------------------------------------------------------
# Hovmöller de frecuencias de θ para una estación (desde el bundle)
# --------------------------------------------------------------------

def plot_omega_hov_from_npz(
    npz_path, station="A", f=6.97e-5, N=1.016e-2,
    out_units="cph", cmap="magma", vmin=0.0, vmax=None, title_prefix=""
):
    d = np.load(npz_path, allow_pickle=True)
    t   = d["t_horas"]        # (nt,)
    RC  = d["RC"]             # (nz,)
    TH3 = d["theta_deg"]      # (nt,nz,nst)
    labs= [str(x) for x in d["labels"]]
    meta= json.loads(str(d["meta_json"]))

    # estación
    if isinstance(station, int):
        k = station; label = labs[k]
    else:
        label = str(station); k = labs.index(label)

    Ztheta = TH3[:, :, k]           # (nt,nz)
    th = np.deg2rad(Ztheta)
    omega = np.sqrt((f**2)*np.sin(th)**2 + (N**2)*np.cos(th)**2)  # s^-1 (angular)

    # unidades para visualizar
    if out_units.lower() == "cph":   # ciclos por hora
        omega = omega/(2*np.pi)*3600.0
        cbar_label = r"$\omega$ (ciclos/h)"
    elif out_units.lower() == "cpd": # ciclos por día
        omega = omega/(2*np.pi)*86400.0
        cbar_label = r"$\omega$ (ciclos/día)"
    else:                            # s^-1 (angular)
        cbar_label = r"$\omega$ (s$^{-1}$)"

    if vmax is None or not np.isfinite(vmax) or vmax<=0:
        vmax = float(np.nanmax(omega))
        if vmax <= 0: vmax = 1e-6

    fig, ax = plt.subplots(figsize=(9,5))
    pcm = ax.pcolormesh(t, RC, omega.T, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(pcm, ax=ax, label=cbar_label)
    ax.set_xlabel("Tiempo (h)")
    ax.set_ylabel("Profundidad (m)")
    ax.set_title(f"{title_prefix} ω(t,z) — {meta.get('scenario','')} • est. {label}")
    plt.tight_layout(); plt.show()
    return fig, ax

# --------------------------------------------------------------------
# Histoframa de θ para una estación (desde el bundle)
# --------------------------------------------------------------------
def plot_freq_hist_at_depth(
    npz_path: str,
    depth_m: float = 110.0,
    f: float = 6.97e-5,          # s^-1
    N: float = 1.016e-2,         # s^-1
    units: str = "cph",          # 's^-1' | 'cph' | 'cpd'
    bins: int = 40,
    density: bool = True,
    separate: bool = False       # False=overlay, True=subplots
):
    d = np.load(npz_path, allow_pickle=True)
    t    = d["t_horas"]               # (nt,)
    RC   = d["RC"]                    # (nz,)
    TH3  = d["theta_deg"]             # (nt,nz,nst)
    labs = [str(x) for x in d["labels"]]

    # Elegir el nivel más cercano a 'depth_m' (tolerante a signos)
    idx = int(np.nanargmin(np.abs(np.abs(RC) - abs(depth_m))))
    depth_sel = float(RC[idx])

    # θ en radianes en ese nivel para todas las estaciones: (nt, nst)
    th = np.deg2rad(TH3[:, idx, :])

    # ω en s^-1 (frecuencia angular, ya positiva)
    omega = np.sqrt((f**2) * np.sin(th)**2 + (N**2) * np.cos(th)**2)

    # Conversión de unidades para el histograma
    if units.lower() == "cph":
        omega = omega / (2*np.pi) * 3600.0
        xlab = r"$\omega$ (ciclos/h)"
        f_ref, N_ref = f/(2*np.pi)*3600.0, N/(2*np.pi)*3600.0
    elif units.lower() == "cpd":
        omega = omega / (2*np.pi) * 86400.0
        xlab = r"$\omega$ (ciclos/día)"
        f_ref, N_ref = f/(2*np.pi)*86400.0, N/(2*np.pi)*86400.0
    else:
        xlab = r"$\omega$ (s$^{-1}$)"
        f_ref, N_ref = f, N

    # Graficar
    if not separate:
        fig, ax = plt.subplots(figsize=(9,5))
        for k, lab in enumerate(labs):
            data = omega[:, k]
            data = data[np.isfinite(data)]
            if data.size:
                ax.hist(data, bins=bins, histtype="step", density=density, label=lab)
        # referencias f y N
        ax.axvline(f_ref, ls="--", lw=1, alpha=0.7, label="f")
        ax.axvline(N_ref, ls=":",  lw=1, alpha=0.7, label="N")
        ax.set_xlabel(xlab)
        ax.set_ylabel("Densidad" if density else "Conteo")
        ax.set_title(f"Histograma de ω a {abs(depth_sel):.0f} m (nivel más cercano)")
        ax.set_xlim(left=0, right=0.8)
        ax.set_yscale("log")
        ax.legend(ncols=2)
        plt.tight_layout()
        return fig, ax, {"depth_level": depth_sel, "index": idx}
    else:
        nst = len(labs)
        cols = 3
        rows = int(np.ceil(nst / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(4.5*cols, 3.5*rows), squeeze=False)
        for k, lab in enumerate(labs):
            r, c = divmod(k, cols)
            ax = axes[r, c]
            data = omega[:, k]
            data = data[np.isfinite(data)]
            if data.size:
                ax.hist(data, bins=bins, density=density)
            ax.set_title(lab)
            ax.set_xlabel(xlab)
            ax.set_ylabel("Densidad" if density else "Conteo")
            ax.axvline(f_ref, ls="--", lw=1, alpha=0.6)
            ax.axvline(N_ref, ls=":",  lw=1, alpha=0.6)
        # quitar ejes vacíos
        for j in range(nst, rows * cols):
            r, c = divmod(j, cols)
            fig.delaxes(axes[r, c])
        fig.suptitle(f"Histograma de ω a {abs(depth_sel):.0f} m (nivel más cercano)")
        plt.tight_layout()
        return fig, axes
    



# --------------------------------------------------------------------
# Quiver θ(t,z) para una estación (desde el bundle)

def plot_quiver_from_theta_npz(
    npz_path: str,
    station: str | int = "A",
    angle_from: str | None = None,  # si None, usa el del meta del bundle
    t_stride: int = 4,
    z_stride: int = 2,
    arrow_len: float = 0.08,        # longitud fija de flecha (en unidades del eje)
    alpha: float = 0.9,
    background: str = "none",       # "none" | "angle"
    cmap = "twilight_shifted"
):
    """
    Dibuja quiver sobre (tiempo, profundidad) para una estación,
    usando el bundle de θ(t,z) guardado previamente.
    - Usa flechas UNITARIAS (solo dirección).
    - 'angle_from': "vertical" (atan2(U,W)) o "horizontal" (atan2(W,U)).
      Si es None, se lee del meta.
    - 'background="angle"' pinta un mapa de θ detrás (opcional).
    """
    d = np.load(npz_path, allow_pickle=True)
    t   = d["t_horas"]         # (nt,)
    RC  = d["RC"]              # (nz,)
    TH3 = d["theta_deg"]       # (nt, nz, nst)
    labs= [str(x) for x in d["labels"]]
    meta= json.loads(str(d["meta_json"]))

    if angle_from is None:
        angle_from = meta.get("angle_from", "vertical")

    # Selección estación
    if isinstance(station, int):
        kst = station
        label = labs[kst]
    else:
        label = str(station)
        if label not in labs:
            raise KeyError(f"Estación '{label}' no encontrada. Disponibles: {labs}")
        kst = labs.index(label)

    TH = TH3[:, :, kst]  # (nt, nz) en grados
    THr = np.radians(TH)

    # Submuestreo
    t_idx = np.arange(0, t.size, max(1, t_stride))
    z_idx = np.arange(0, RC.size, max(1, z_stride))
    Tm, Zm = np.meshgrid(t[t_idx], RC[z_idx], indexing="xy")      # (nz_sub, nt_sub)
    Ang = THr[np.ix_(t_idx, z_idx)].T                             # (nz_sub, nt_sub)

    # Componentes de dirección unitaria en (x,z)
    if angle_from == "vertical":
        Udir = np.sin(Ang)   # horizontal
        Wdir = np.cos(Ang)   # vertical
    elif angle_from == "horizontal":
        Udir = np.cos(Ang)
        Wdir = np.sin(Ang)
    else:
        raise ValueError("angle_from debe ser 'vertical' o 'horizontal'.")

    # Longitud fija de flecha (solo dirección)
    Uq = arrow_len * Udir
    Wq = arrow_len * Wdir

    fig, ax = plt.subplots(figsize=(10, 5))

    if background == "angle":
        im = ax.pcolormesh(
            t, RC, TH.T, shading="auto", cmap=cmap,
            vmin=-180, vmax=180
        )
        cbar = fig.colorbar(im, ax=ax, label="θ (°)")

    q = ax.quiver(
        Tm, Zm, Uq, Wq,
        angles="xy", scale_units="xy", scale=1.0,
        width=0.0025, headwidth=3.2, headlength=4.0,
        alpha=alpha
    )

    ax.set_xlabel("Tiempo (h)")
    ax.set_ylabel("Profundidad (m)")
    ax.set_title(f"Dirección u/w (quiver) — {meta.get('scenario','')} • est. {label} "
                 f"[ref: {angle_from}]")
    ax.invert_yaxis()  # si prefieres profundidades creciendo hacia abajo
    ax.grid(True, ls=":", alpha=0.4)
    plt.tight_layout()
    plt.show()
    return fig, ax

### ─────────────────────────────────────────────────────────────
# Componente barotrópica (promedio vertical) y componente baroclicica
def guardar_barotropico_baroclinico_escenario(
    outdir_bay: Path,
    outdir_nobay: Path,
    save_dir: Path,
    varname: str = "U",
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    depths_m=(110.0, 250.0),
    stage_patterns=("state.0000000000.t*.nc",),
    prefer_RF: bool = False
):
    """
    Calcula y guarda las componentes barotrópicas y baroclínicas
    de una variable (U o W) para un escenario dado.
    """

    # Cargar datos del MITgcm
    data = cargar_datos_mitgcm_etapas(outdir_bay, outdir_nobay, stage_patterns=stage_patterns)
    gr = data["gr_bay"]
    stations = make_stations_from_grid(gr, j_centro=j_centro, radio=radio)

    RC = _get_RC(gr)
    nz = RC.size

    # Crear carpeta de salida
    save_dir.mkdir(parents=True, exist_ok=True)

    # Procesar estaciones
    for i, j, label, *_ in stations:
        times, perfiles = [], []
        for nc in data["nc_bay_list"]:
            t = _get_time_from_nc(nc)
            arr = np.array(nc.variables[varname][:, :, j, i], dtype=float)  # (nt, nz)
            if arr.shape[1] != nz:
                raise ValueError(f"{varname}: número de niveles no coincide.")
            times.append(t)
            perfiles.append(arr)

        if not times:
            continue

        # Concatenar tiempo y datos
        T, A = _concat_time_and_data(times, perfiles)  # T:(nt,), A:(nt,nz)

        # Para cada profundidad de interés
        for depth in depths_m:
            mask_levels = np.abs(RC) <= depth
            k_sel = np.where(mask_levels)[0]
            if len(k_sel) == 0:
                continue

            A_sel = A[:, k_sel]  # (nt, nz_sel)

            # Barotrópico
            barotrop = np.nanmean(A_sel, axis=1)  # (nt,)

            # Baroclínico
            baroclin = A_sel - barotrop[:, None]  # (nt, nz_sel)

            # Guardar
            outfile = save_dir / f"{varname}_{label}_{int(depth)}m_baroclinico.npz"
            np.savez(
                outfile,
                t=T,
                RC=RC[k_sel],
                barotropica=barotrop,
                baroclinica=baroclin
            )
            print(f"Guardado: {outfile}")


def procesar_varios_escenarios_barotropico_baroclinico(
    scenarios: dict,
    subdir_bay: str,
    subdir_nobay: str,
    save_dir: Path,
    varname: str = "U",
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    depths_m=(110.0, 250.0),
    stage_patterns=("state.0000000000.t*.nc",),
    prefer_RF: bool = False
):
    """
    Ejecuta guardar_barotropico_baroclinico_escenario para varios escenarios.
    """
    for scen_name, base in scenarios.items():
        print(f"\nProcesando escenario: {scen_name}")
        outdir_bay = base / subdir_bay
        outdir_nobay = base / subdir_nobay
        scen_save_dir = save_dir / scen_name
        guardar_barotropico_baroclinico_escenario(
            outdir_bay,
            outdir_nobay,
            scen_save_dir,
            varname=varname,
            j_centro=j_centro,
            radio=radio,
            depths_m=depths_m,
            stage_patterns=stage_patterns,
            prefer_RF=prefer_RF
        )

# === viz_barotropic_baroclinic.py ===

def _to_hours(t):
    """Asume t en segundos si es muy grande; si ya está en horas, lo deja igual."""
    t = np.asarray(t, dtype=float)
    if t.size == 0:
        return t
    # Heurística: si el paso típico es ~O(10-100) s, conviene a horas
    dt = np.nanmedian(np.diff(t)) if t.size > 1 else t[0]
    if dt and dt > 120.0:  # mayor a 2 min → probablemente segundos
        return t / 3600.0
    # si ya está en horas o vacíos pequeños, no convertir
    return t

def _human_station_label(npz_path: Path):
    """Extrae estación (p.ej. Seiche, A, B, C) y profundidad del nombre del archivo."""
    # Espera algo como: U_Seiche_110m_baroclinico.npz
    m = re.search(r"_(?P<label>[^_/]+)_(?P<depth>\d+)m_", npz_path.name)
    if m:
        return m.group("label"), int(m.group("depth"))
    # Fallback
    return npz_path.stem, None

def plot_baro_npz(npz_path: Path,
                  title_prefix: str = "",
                  cmap: str = "cmo.balance" if "cmo" in plt.colormaps() else "RdBu_r",
                  vmin=None, vmax=None,
                  show=True):
    """
    Lee un .npz con claves: t, RC, barotropica, baroclinica
    y grafica:
      - (arriba) serie del barotrópico
      - (abajo) Hovmöller tiempo–profundidad del baroclínico
    """
    npz_path = Path(npz_path)
    d = np.load(npz_path, allow_pickle=False)
    t = d["t"]
    RC = d["RC"]               # niveles (negativos)
    bt = d["barotropica"]      # (nt,)
    bc = d["baroclinica"]      # (nt, nz_sel)

    th = _to_hours(t)
    # Escala para Hovmöller
    if vmin is None or vmax is None:
        m = float(np.nanmax(np.abs(bc))) if bc.size else 1.0
        vmin, vmax = -m, m

    label, depth = _human_station_label(npz_path)
    fig = plt.figure(figsize=(10, 7))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 1.4], hspace=0.25)

    # --- Serie barotrópica ---
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(th, bt, lw=1.5)
    ax1.axhline(0, color="k", lw=0.8, alpha=0.6)
    ax1.set_xlabel("Tiempo (h)")
    ax1.set_ylabel("Barotrópico")
    ttl1 = f"{title_prefix} {label} • {depth} m" if depth else f"{title_prefix} {label}"
    ax1.set_title(ttl1.strip())

    # --- Hovmöller baroclínico ---
    ax2 = fig.add_subplot(gs[1, 0])
    # Nota: RC es negativo (profundidad). Para mostrar “hacia abajo” usamos RC tal cual.
    pcm = ax2.pcolormesh(th, RC, bc.T, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
    cbar = fig.colorbar(pcm, ax=ax2, label="Baroclínico (anomalía)")
    ax2.set_xlabel("Tiempo (h)")
    ax2.set_ylabel("Profundidad (m)")
    ax2.set_title("Evolución baroclínica (tiempo–profundidad)")

    plt.tight_layout()
    if show:
        plt.show()
    return fig, (ax1, ax2)

def _find_npz(root: Path, varname: str, station: str, depth_m: int):
    """
    Busca archivo con patrón: <var>_<station>_<depth>m_baroclinico.npz
    dentro de root (recursivo).
    """
    pattern = f"{varname}_{station}_{int(depth_m)}m_baroclinico.npz"
    matches = list(root.rglob(pattern))
    if not matches:
        raise FileNotFoundError(f"No encontré {pattern} en {root}")
    # devolver el más “reciente” por nombre si hay varios
    matches.sort()
    return matches[-1]

def plot_baro_from_dir(root: Path,
                       scenario: str,
                       varname: str,
                       station: str,
                       depth_m: int,
                       cmap: str = "cmo.balance" if "cmo" in plt.colormaps() else "RdBu_r",
                       vmin=None, vmax=None,
                       show=True):
    """
    Busca el .npz de un escenario dado y grafica con plot_baro_npz.
    Espera estructura tipo:
      root/
        realistic/
          U_Seiche_110m_baroclinico.npz
        linear/
        2layer/
    """
    root = Path(root)
    npz_path = _find_npz(root / scenario, varname, station, depth_m)
    title_prefix = f"{scenario} • {varname}"
    return plot_baro_npz(npz_path, title_prefix=title_prefix, cmap=cmap, vmin=vmin, vmax=vmax, show=show)

def compare_bartropic_across_scenarios(root: Path,
                                       scenarios: list[str],
                                       varname: str,
                                       station: str,
                                       depth_m: int,
                                       show=True):
    """
    Compara la SERIE BAROTRÓPICA entre escenarios para una estación/depth/var.
    """
    root = Path(root)
    fig, ax = plt.subplots(figsize=(10, 4.2))
    for scen in scenarios:
        npz_path = _find_npz(root / scen, varname, station, depth_m)
        d = np.load(npz_path, allow_pickle=False)
        t = _to_hours(d["t"])
        bt = d["barotropica"]
        ax.plot(t, bt, lw=1.25, label=scen)
    ax.axhline(0, color="k", lw=0.8, alpha=0.6)
    ax.set_xlabel("Tiempo (h)")
    ax.set_ylabel("Barotrópico")
    ax.set_title(f"Comparación barotrópica • {varname} • {station} • {depth_m} m")
    ax.legend()
    plt.tight_layout()
    if show:
        plt.show()
    return fig, ax

# === export_pdf_escenario.py ===

def _to_hours(t):
    t = np.asarray(t, dtype=float)
    if t.size == 0: return t
    dt = np.nanmedian(np.diff(t)) if t.size > 1 else t[0]
    return t/3600.0 if (dt and dt > 120.0) else t  # si parece segundos → a horas

def _load_npz(root: Path, scenario: str, varname: str, station: str, depth_m: int):
    """
    Busca <var>_<station>_<depth>m_baroclinico.npz dentro de root/scenario y lo carga.
    Devuelve dict con claves: t, RC, barotropica, baroclinica.
    """
    root = Path(root)
    pat = f"{varname}_{station}_{int(depth_m)}m_baroclinico.npz"
    matches = list((root / scenario).rglob(pat))
    if not matches:
        raise FileNotFoundError(f"No encontré {pat} en {(root / scenario)}")
    matches.sort()
    d = np.load(matches[-1], allow_pickle=False)
    return dict(
        t=d["t"], RC=d["RC"], barotropica=d["barotropica"], baroclinica=d["baroclinica"]
    )

def _get_cmap():
    # usa cmocean si está disponible, si no cae a RdBu_r
    try:
        import cmocean as cmo
        return cmo.cm.balance
    except Exception:
        return "RdBu_r"

def _common_limits_baroclinic(bc_list):
    m = 0.0
    for bc in bc_list:
        if bc is None or bc.size == 0: continue
        m = max(m, float(np.nanmax(np.abs(bc))))
    if m == 0: m = 1.0
    return -m, m

def export_pdf_por_escenario(
    root_npz: Path,
    scenario: str,
    stations=("Seiche","A","B","C"),
    depths=(110, 250),
    varnames=("U","W"),
    outdir: Path = Path("./figs"),
    title_prefix: str | None = None
) -> Path:
    """
    Crea un PDF por escenario con comparativas por profundidad:
      - Izquierda: serie BAROTRÓPICA (todas las estaciones) para var/ profundidad.
      - Derecha: grid 2×2 de Hovmöller BAROCLÍNICO por estación (escala común).

    Estructura del PDF:
      Por cada varname en varnames:
        Por cada depth en depths:
          Página con 1 (serie) + 4 Hovmöllers (estaciones).
    """
    root_npz = Path(root_npz)
    outdir.mkdir(parents=True, exist_ok=True)
    pdf_path = outdir / f"{scenario}_comparativo_baro.pdf"

    cmap = _get_cmap()

    with PdfPages(pdf_path) as pdf:
        for var in varnames:
            for depth in depths:
                # --- Cargar todo primero para fijar escala común de baroclínico
                loaded = {}
                bc_list = []
                for st in stations:
                    try:
                        d = _load_npz(root_npz, scenario, var, st, depth)
                        loaded[st] = d
                        bc_list.append(d["baroclinica"])
                    except FileNotFoundError:
                        loaded[st] = None
                        bc_list.append(None)

                vmin, vmax = _common_limits_baroclinic(bc_list)

                # --- Crear figura: 1 fila (serie) + 2x2 Hovmöller
                fig = plt.figure(figsize=(13, 8))
                # GridSpec: 2 filas, 3 columnas: 
                #   fila 0: ax_series ocupa col 0-2 (ancho completo)
                #   fila 1: 4 ejes en 2x2 para Hovmöller
                gs = fig.add_gridspec(nrows=2, ncols=3, height_ratios=[1.0, 1.6], hspace=0.3, wspace=0.25)
                ax_series = fig.add_subplot(gs[0, :])

                # --- Serie barotrópica comparando estaciones
                any_ok = False
                for st in stations:
                    d = loaded[st]
                    if d is None: continue
                    t_h = _to_hours(d["t"])
                    ax_series.plot(t_h, d["barotropica"], lw=1.5, label=st)
                    any_ok = True
                ax_series.axhline(0, color="k", lw=0.8, alpha=0.6)
                ax_series.set_xlabel("Tiempo (h)")
                ax_series.set_ylabel("Barotrópico")
                title_left = f"{scenario} • {var} • {depth} m — Serie barotrópica (todas las estaciones)"
                ax_series.set_title(title_left)
                if any_ok:
                    ax_series.legend(ncol=len(stations), fontsize="small")

                # --- Hovmöllers 2x2
                # Orden visual fijo por estaciones
                hov_axes = [
                    fig.add_subplot(gs[1, 0]),
                    fig.add_subplot(gs[1, 1]),
                    fig.add_subplot(gs[1, 2]),
                    fig.add_axes([0.69, 0.16, 0.27, 0.35])  # último panel manual si quieres 2x2 + colorbar aparte
                ]
                # Para hacer realmente 2x2 “limpio”, usa subgridspec:
                fig.delaxes(hov_axes[-1])  # eliminamos ese axes manual
                sub = gs[1, :].subgridspec(2, 2, hspace=0.28, wspace=0.25)
                hov_axes = [
                    fig.add_subplot(sub[0,0]), fig.add_subplot(sub[0,1]),
                    fig.add_subplot(sub[1,0]), fig.add_subplot(sub[1,1]),
                ]

                for ax, st in zip(hov_axes, stations):
                    d = loaded[st]
                    ax.set_title(f"{st} — baroclínico")
                    if d is None:
                        ax.text(0.5, 0.5, "No data", ha="center", va="center")
                        ax.set_xlabel("Tiempo (h)")
                        ax.set_ylabel("Profundidad (m)")
                        continue
                    th = _to_hours(d["t"])
                    RC = d["RC"]   # niveles negativos (m)
                    bc = d["baroclinica"]  # (nt, nz_sel)
                    pm = ax.pcolormesh(th, RC, bc.T, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
                    ax.set_xlabel("Tiempo (h)")
                    ax.set_ylabel("Profundidad (m)")
                    # para no atiborrar ticks, puedes tunear aquí

                # Colorbar común (usa el último pm válido)
                # Busca el último pcolormesh creado
                last_pm = None
                for ax, st in zip(hov_axes, stations):
                    for im in ax.get_children():
                        if hasattr(im, 'get_array'):
                            last_pm = im
                if last_pm is not None:
                    cbar = fig.colorbar(last_pm, ax=hov_axes, fraction=0.03, pad=0.02)
                    cbar.set_label("Baroclínico (anomalía)")

                # Título global de página
                page_title = f"{scenario} • {var} • {depth} m"
                if title_prefix: page_title = f"{title_prefix} — {page_title}"
                fig.suptitle(page_title, y=0.98, fontsize=12, fontweight="bold")

                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

    print(f"[OK] PDF exportado: {pdf_path}")
    return pdf_path

def _get_vertical_metrics(gr, varname: str):
    """
    Devuelve (RC, drF, hFac_kji_selector) para promedios verticales.
    - RC: niveles (negativos) 1D (nz,)
    - drF: espesores (nz,) si existe; fallback = ones
    - hFac: función hfac(j,i)->(nz,), basada en varname, con fallbacks robustos
    """
    RC = _get_RC(gr)
    nz = RC.size

    # drF (espesor por nivel)
    if 'drF' in gr.variables:
        drF = np.array(gr.variables['drF'][:], dtype=float).squeeze()
        if drF.ndim > 1: drF = drF[..., 0, 0]
        if drF.size != nz: drF = np.ones(nz, float)
    else:
        drF = np.ones(nz, float)

    # Selección hFac según grid de la variable (fallbacks seguros)
    # U → HFacW (U-grid) ; W → HFacC (C-grid) (si no, usa el que exista)
    hfac_name_candidates = []
    if varname.upper().startswith('U'):
        hfac_name_candidates = ['HFacW', 'HFacC', 'HFacS']
    else:
        # W (vertical) suele cruzar interfaces; usamos HFacC como aproximación robusta
        hfac_name_candidates = ['HFacC', 'HFacW', 'HFacS']

    hFac3D = None
    for hn in hfac_name_candidates:
        if hn in gr.variables:
            hFac3D = np.array(gr.variables[hn][:], dtype=float)
            break

    def hfac_profile(j, i):
        # Devuelve (nz,) de factores [0..1] validando dimensiones
        if hFac3D is None:
            return np.ones(nz, float)
        # hFac puede venir (nz,ny,nx) o (ny,nx)→ fallback
        if hFac3D.ndim == 3:
            return hFac3D[:, j, i].astype(float, copy=False)
        elif hFac3D.ndim == 2:
            # si fuera 2D (raro), extiende a nz
            f = hFac3D[j, i]
            return np.full(nz, float(f))
        else:
            return np.ones(nz, float)

    return RC, drF, hfac_profile


def uv_to_cgrid(u, v, hFacW=None, hFacS=None, fill_value=np.nan):
    """
    Lleva U/V en mallas desfasadas (xp1/yp1) a la malla C por promedio ponderado.

    Parámetros
    ----------
    u : array-like (..., nz, ny, nx+1)   # velocidad U en U-grid (staggered en x)
    v : array-like (..., nz, ny+1, nx)   # velocidad V en V-grid (staggered en y)
    hFacW : array-like (nz, ny, nx+1) o None
        Fracción abierta en caras oeste/este. Si None, asume 1.0.
    hFacS : array-like (nz, ny+1, nx) o None
        Fracción abierta en caras sur/norte. Si None, asume 1.0.
    fill_value : valor para celdas sin apertura (denominador=0).

    Retorna
    -------
    uc, vc : arrays (..., nz, ny, nx) en malla C.
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    hw = 1.0 if hFacW is None else np.asarray(hFacW, dtype=float)
    hs = 1.0 if hFacS is None else np.asarray(hFacS, dtype=float)

    # U: promedio en x (xp1 -> x)
    u_num = u[..., 1:] * hw[..., 1:] + u[..., :-1] * hw[..., :-1]
    u_den = hw[..., 1:] + hw[..., :-1]
    uc = np.where(u_den > 0, u_num / u_den, fill_value)

    # V: promedio en y (yp1 -> y)
    v_num = v[..., 1:, :] * hs[..., 1:, :] + v[..., :-1, :] * hs[..., :-1, :]
    v_den = hs[..., 1:, :] + hs[..., :-1, :]
    vc = np.where(v_den > 0, v_num / v_den, fill_value)

    return uc, vc


def w_to_cgrid(w, hFacC=None, fill_value=np.nan):
    """
    Centra W si viene en interfaces verticales (nz+1) → malla C (nz).
    Si ya tiene shape (..., nz, ny, nx) de centros, lo devuelve tal cual.
    """
    w = np.asarray(w, dtype=float)
    if w.ndim < 4:
        raise ValueError("w debe ser (nt, nz, ny, nx) o (nt, nz+1, ny, nx).")

    nz_faces = w.shape[-3]
    nz_centers = None if hFacC is None else np.asarray(hFacC).shape[0]

    if nz_centers is not None and nz_faces == nz_centers:
        wc = w  # ya en centros
    elif nz_centers is not None and nz_faces == nz_centers + 1:
        wc = 0.5 * (w[..., 1:, :, :] + w[..., :-1, :, :])
    elif nz_centers is None:
        # sin hFac, si tiene nz+1 niveles asumimos que hay que centrar
        if nz_faces > 1:
            wc = 0.5 * (w[..., 1:, :, :] + w[..., :-1, :, :])
        else:
            wc = w
    else:
        raise ValueError(f"No puedo centrar W: niveles w={nz_faces}, hFacC={nz_centers}.")

    if hFacC is not None:
        wc = np.where(np.asarray(hFacC, dtype=float) > 0, wc, fill_value)
    return wc


def kinetic_energy_series(u, v, w, drF, area=None, dx=None, dy=None, hFacC=None, rho0=1025.0):
    """
    Energía cinética total por tiempo integrando volumen: 0.5*rho*(u^2+v^2+w^2)*vol.

    Parámetros
    ----------
    u, v, w : arrays (nt, nz, ny, nx) en malla C (usa uv_to_cgrid/w_to_cgrid antes si hace falta).
    drF     : array  (nz,) espesores verticales (p.ej. drF del grid).
    area    : opcional (ny, nx) con área de celda (p.ej. rAc). Si se pasa, ignora dx/dy.
    dx, dy  : opcionales (ny, nx) métricas horizontales; requeridas si no se pasa area.
    hFacC   : opcional (nz, ny, nx) fracción abierta en centro; si None asume 1.
    rho0    : densidad [kg/m³].

    Retorna
    -------
    ke_t : array (nt,) energía cinética integrada [J].
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    w = np.asarray(w, dtype=float)
    drF = np.asarray(drF, dtype=float)
    hfac = 1.0 if hFacC is None else np.asarray(hFacC, dtype=float)

    if area is None:
        if dx is None or dy is None:
            raise ValueError("Debes pasar 'area' (rAc) o bien 'dx' y 'dy'.")
        dx = np.asarray(dx, dtype=float)
        dy = np.asarray(dy, dtype=float)
        if dx.shape != dy.shape:
            raise ValueError("dx y dy deben tener misma forma.")
        area = dx * dy
    else:
        area = np.asarray(area, dtype=float)
    # Verificaciones mínimas
    nt, nz, ny, nx = u.shape
    drF = np.asarray(drF, dtype=float)
    if v.shape != (nt, nz, ny, nx) or w.shape != (nt, nz, ny, nx):
        raise ValueError("u, v, w deben tener misma forma (nt, nz, ny, nx).")
    if area.shape != (ny, nx):
        raise ValueError("area (o dx/dy) debe ser (ny, nx).")
    if drF.shape[0] != nz:
        raise ValueError("drF debe ser (nz,).")
    if not np.isscalar(hfac) and hfac.shape != (nz, ny, nx):
        raise ValueError("hFacC debe ser (nz, ny, nx) o None.")

    # Volumen por celda (broadcast a (nz, ny, nx))
    vol_3d = (drF[:, None, None] * area[None, :, :]).astype(np.float64, copy=False)
    if not np.isscalar(hfac):
        vol_3d = vol_3d * hfac

    # Cálculo sin crear arreglos 4D extra: einsum combina y reduce "en el aire".
    # Nota: si uv_to_cgrid dejó NaN donde no hay celda abierta, los convertimos a 0.0
    # para que no contaminen la suma (el volumen allí suele ser 0 o muy pequeño).
    u0 = np.nan_to_num(u, nan=0.0, copy=False)
    v0 = np.nan_to_num(v, nan=0.0, copy=False)
    w0 = np.nan_to_num(w, nan=0.0, copy=False)

    ke_u = np.einsum("tzyx,tzyx,zyx->t", u0, u0, vol_3d, optimize=True)
    ke_v = np.einsum("tzyx,tzyx,zyx->t", v0, v0, vol_3d, optimize=True)
    ke_w = np.einsum("tzyx,tzyx,zyx->t", w0, w0, vol_3d, optimize=True)

    return 0.5 * float(rho0) * (ke_u + ke_v + ke_w)


def _slice_plus_one(slc: slice) -> slice:
    """Devuelve un slice con stop+1 (si stop no es None)."""
    if not isinstance(slc, slice):
        raise TypeError("Se esperaba un slice.")
    if slc.stop is None:
        return slice(slc.start, None, slc.step)
    return slice(slc.start, slc.stop + 1, slc.step)


def ke_time_series_from_mnc_state(
    nc_state,
    gr,
    nz_sel: slice,
    y_slice: slice,
    x_slice: slice,
    ny_block: int = 64,
    rho0: float = 1025.0,
    var_u: str = "U",
    var_v: str = "V",
    var_w: str = "W",
    use_float32: bool = True,
):
    """
    Calcula energía cinética integrada por tiempo para un archivo MNC state (una etapa).
    Trabaja por bloques en Y para no reventar memoria.

    Retorna (t_seconds, ke_t).
    """
    # tiempo del state
    t = _get_time_from_nc(nc_state)
    nt = int(t.shape[0])

    # grid/métricas (se cortan por bloque)
    drF = np.asarray(gr.variables["drF"][nz_sel], dtype=float).squeeze()

    # rangos en y (C-grid)
    if y_slice.start is None or y_slice.stop is None:
        raise ValueError("y_slice debe tener start/stop explícitos (p.ej. slice(0,323)).")
    y0, y1 = int(y_slice.start), int(y_slice.stop)

    ke_acc = np.zeros(nt, dtype=np.float64)

    for j0 in range(y0, y1, int(ny_block)):
        j1b = min(y1, j0 + int(ny_block))

        # --- cortes staggered consistentes ---
        x_u = _slice_plus_one(x_slice) if x_slice.stop is not None else x_slice
        y_v = _slice_plus_one(slice(j0, j1b))  # yp1
        y_c = slice(j0, j1b)

        # hFac y métricas horizontales
        _hFacC_name = _pick_var(gr, ("hFacC", "HFacC"))
        _hFacW_name = _pick_var(gr, ("hFacW", "HFacW"))
        _hFacS_name = _pick_var(gr, ("hFacS", "HFacS"))
        _rA_name    = _pick_var(gr, ("rAc",   "rA"))
        hC = np.asarray(gr.variables[_hFacC_name][nz_sel, y_c, x_slice], dtype=float)
        hW = np.asarray(gr.variables[_hFacW_name][nz_sel, y_c, x_u], dtype=float)
        hS = np.asarray(gr.variables[_hFacS_name][nz_sel, y_v, x_slice], dtype=float)
        area = np.asarray(gr.variables[_rA_name][y_c, x_slice], dtype=float)

        # U/V/W
        U_blk = np.asarray(nc_state.variables[var_u][:, nz_sel, y_c, x_u])
        V_blk = np.asarray(nc_state.variables[var_v][:, nz_sel, y_v, x_slice])

        W_var = nc_state.variables[var_w]
        nz_w = int(getattr(W_var, "shape", (None, None))[1]) if len(getattr(W_var, "shape", ())) >= 2 else None
        nz_centers = int(hC.shape[0])
        if nz_w == nz_centers + 1:
            nz_w_sel = slice(nz_sel.start, nz_sel.stop + 1, nz_sel.step)
        else:
            nz_w_sel = nz_sel
        W_blk = np.asarray(W_var[:, nz_w_sel, y_c, x_slice])

        if use_float32:
            U_blk = U_blk.astype(np.float32, copy=False)
            V_blk = V_blk.astype(np.float32, copy=False)
            W_blk = W_blk.astype(np.float32, copy=False)

        uc, vc = uv_to_cgrid(U_blk, V_blk, hFacW=hW, hFacS=hS, fill_value=np.nan)
        wc = w_to_cgrid(W_blk, hFacC=hC, fill_value=np.nan)
        ke_acc += kinetic_energy_series(uc, vc, wc, drF=drF, area=area, hFacC=hC, rho0=rho0)

    return t.astype(np.float64), ke_acc


def guardar_ke_diff_npz_etapas(
    scenario_root: Path,
    scenario_label: str,
    subdir_bay: str = "run_expand",
    subdir_nobay: str = "run_expand_nobay",
    mnc_glob: str = "mnc_00*",
    stage_patterns=("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern: str = "grid.t*.nc",
    nz_sel: slice = slice(0, 34),
    y_slice: slice = slice(0, 323),
    x_slice: slice = slice(None),
    ny_block: int = 64,
    rho0: float = 1025.0,
    save_dir: Path | None = None,
    out_name: str | None = None,
    tol_time: float = 1e-7,
    recursive_search: bool = True,
):
    """
    Calcula KE(t) para bay y nobay concatenando dos (o más) etapas y guarda un NPZ
    con KE_bay, KE_nobay y KE_diff=bay−nobay.

    Guarda en `save_dir` (por defecto, carpeta del script) un archivo:
      <scenario_label>_KEdiff_nz0-34_y0-323.npz  (o `out_name`)
    """
    scenario_root = Path(scenario_root)
    if save_dir is None:
        save_dir = Path(__file__).resolve().parent
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    def _glob_any(pattern: str, recursive: bool):
        return sorted(glob.glob(str(scenario_root / pattern), recursive=recursive))

    def _open_mnc_first(patterns: list[str]):
        for patt in patterns:
            matches = _glob_any(patt, recursive=recursive_search and ("**" in patt))
            if matches:
                return MIT.mnc_files(str(scenario_root / patt))
        raise FileNotFoundError(
            "No encontré archivos MNC. Intenté patrones:\n  - "
            + "\n  - ".join(patterns)
            + f"\nEn: {scenario_root}"
        )

    # Grid desde bay (asumimos misma malla para nobay)
    gr = _open_mnc_first([
        f"{subdir_bay}/{mnc_glob}/{grid_pattern}",
        f"{subdir_bay}/{mnc_glob}/grid*.nc",
        f"{subdir_bay}/**/{grid_pattern}",
        f"{subdir_bay}/**/grid*.nc",
    ])

    tb_list, tnb_list = [], []
    keb_list, kenb_list = [], []

    for patt in stage_patterns:
        nc_b = _open_mnc_first([
            f"{subdir_bay}/{mnc_glob}/{patt}",
            f"{subdir_bay}/**/{patt}",
        ])
        nc_n = _open_mnc_first([
            f"{subdir_nobay}/{mnc_glob}/{patt}",
            f"{subdir_nobay}/**/{patt}",
        ])

        t_b, ke_b = ke_time_series_from_mnc_state(
            nc_b, gr, nz_sel=nz_sel, y_slice=y_slice, x_slice=x_slice,
            ny_block=ny_block, rho0=rho0
        )
        t_n, ke_n = ke_time_series_from_mnc_state(
            nc_n, gr, nz_sel=nz_sel, y_slice=y_slice, x_slice=x_slice,
            ny_block=ny_block, rho0=rho0
        )

        tb_list.append(t_b); keb_list.append(ke_b)
        tnb_list.append(t_n); kenb_list.append(ke_n)

    t_bay, ke_bay = _concat_time_and_data(tb_list, keb_list)
    t_nobay, ke_nobay = _concat_time_and_data(tnb_list, kenb_list)

    # Alinear tiempos si difieren levemente
    if t_bay.size == 0 or t_nobay.size == 0:
        raise ValueError("Tiempo vacío al concatenar etapas.")

    if t_bay.shape == t_nobay.shape and np.all(np.abs(t_bay - t_nobay) <= tol_time):
        t_common = t_bay
        ke_n = ke_nobay
    else:
        # Interpola nobay a tiempos de bay (típicamente deberían coincidir)
        t_common = t_bay
        ke_n = np.interp(t_common, t_nobay, ke_nobay)

    ke_diff = ke_bay - ke_n

    # nombre de salida
    if out_name is None:
        ztag = f"nz{nz_sel.start or 0}-{nz_sel.stop}"
        ytag = f"y{y_slice.start or 0}-{y_slice.stop}"
        out_name = f"{scenario_label}_KEdiff_{ztag}_{ytag}.npz"

    out_path = save_dir / out_name
    meta = {
        "scenario": scenario_label,
        "root": str(scenario_root),
        "subdir_bay": subdir_bay,
        "subdir_nobay": subdir_nobay,
        "mnc_glob": mnc_glob,
        "stage_patterns": list(stage_patterns),
        "grid_pattern": grid_pattern,
        "nz_sel": [nz_sel.start, nz_sel.stop, nz_sel.step],
        "y_slice": [y_slice.start, y_slice.stop, y_slice.step],
        "x_slice": [x_slice.start, x_slice.stop, x_slice.step],
        "ny_block": int(ny_block),
        "rho0": float(rho0),
        "tol_time": float(tol_time),
        "units": {"t": "s", "ke": "J"},
    }

    np.savez_compressed(
        out_path,
        t_seconds=t_common.astype(np.float64),
        t_horas=(t_common / 3600.0).astype(np.float64),
        ke_bay=ke_bay.astype(np.float64),
        ke_nobay=ke_n.astype(np.float64),
        ke_diff=ke_diff.astype(np.float64),
        meta_json=np.array(json.dumps(meta, ensure_ascii=False)),
    )
    print(f"[OK] Guardado KE diff: {out_path}")
    return out_path


# ─────────────────────────────────────────────────────────────
# KE usando densidad (RHOAnoma) + velocidades en malla de densidad
# ─────────────────────────────────────────────────────────────
def unstagger(ugrid, vgrid):
    """Interpola U y V a centros de celda (malla de densidad)."""
    u = np.add(ugrid[..., :-1], ugrid[..., 1:]) / 2
    v = np.add(vgrid[..., :-1, :], vgrid[..., 1:, :]) / 2
    return u[..., :-1, :], v[..., :-1]


def unstagger_w(wgrid):
    """Interpola W a centros verticales (malla de densidad)."""
    w = np.add(wgrid[..., :-1, :, :], wgrid[..., 1:, :, :]) / 2
    return w


def _normalize_slice(slc: slice, n: int) -> slice:
    if not isinstance(slc, slice):
        raise TypeError("Se esperaba un slice.")
    start = 0 if slc.start is None else int(slc.start)
    stop = n if slc.stop is None else int(slc.stop)
    step = 1 if slc.step is None else int(slc.step)
    return slice(start, stop, step)


def _pick_var(gr, names):
    for name in names:
        if name in gr.variables:
            return name
    raise KeyError(f"No encontré ninguna variable en: {names}")


def _pick_var(gr, names):
    for name in names:
        if name in gr.variables:
            return name
    raise KeyError(f"No encontré ninguna variable en: {names}")


def _match_time_indices(t1, t2, tol=1e-7):
    """
    Devuelve indices que hacen coincidir tiempos t1 y t2 dentro de tolerancia.
    Retorna (idx1, idx2, t_common).
    """
    t1 = np.asarray(t1, dtype=float)
    t2 = np.asarray(t2, dtype=float)
    i, j = 0, 0
    idx1, idx2, t_common = [], [], []
    while i < t1.size and j < t2.size:
        dt = t1[i] - t2[j]
        if abs(dt) <= tol:
            idx1.append(i); idx2.append(j); t_common.append(t1[i])
            i += 1; j += 1
        elif dt < 0:
            i += 1
        else:
            j += 1
    return np.array(idx1, dtype=int), np.array(idx2, dtype=int), np.array(t_common, dtype=float)


def _indices_to_slice(idx: np.ndarray) -> slice | None:
    """
    Convierte un array de indices consecutivos a slice. Si no es consecutivo, retorna None.
    """
    idx = np.asarray(idx, dtype=int)
    if idx.size == 0:
        return slice(0, 0, 1)
    if idx.size == 1:
        return slice(int(idx[0]), int(idx[0]) + 1, 1)
    if np.all(np.diff(idx) == 1):
        return slice(int(idx[0]), int(idx[-1]) + 1, 1)
    return None


def ke_time_series_from_state_and_rho_mnc(
    nc_state,
    nc_rho,
    gr,
    z_slice: slice,
    y_slice: slice,
    x_slice: slice,
    ny_block: int = 64,
    rho_const: float = 1025.0,
    var_rho: str = "RHOAnoma",
    var_u: str = "U",
    var_v: str = "V",
    var_w: str = "W",
    tol_time: float = 1e-7,
    use_float32: bool = True,
    rho_clip: float = 1e6,
):
    """
    Calcula KE(t) integrando 0.5*rho*(u^2+v^2+w^2)*vol usando RHOAnoma y
    velocidades, centradas en la malla de densidad con unstagger().
    """
    t_state = _get_time_from_nc(nc_state)
    t_rho = _get_time_from_nc(nc_rho)
    idx_s, idx_r, t_common = _match_time_indices(t_state, t_rho, tol=tol_time)
    if t_common.size == 0:
        raise ValueError("No hay tiempos en común entre state y diag_rho.")
    idx_s_mnc = _indices_to_slice(idx_s)
    # Mismo tratamiento para diag_rho: MITgcmutils.mnc.mnc_files no admite
    # indexado por array en el eje T (calcstrides hace `if Ellipsis in slices`,
    # que con un ndarray lanza "truth value of an array is ambiguous").
    idx_r_mnc = _indices_to_slice(idx_r)

    # métricas
    rA_name = _pick_var(gr, ("rA", "rAc"))
    hFacC_name = _pick_var(gr, ("HFacC", "hFacC"))
    rAc = np.asarray(gr.variables[rA_name][:], dtype=float)
    drF = np.asarray(gr.variables["drF"][z_slice], dtype=float).squeeze()
    ny, nx = rAc.shape
    y_slice = _normalize_slice(y_slice, ny)
    x_slice = _normalize_slice(x_slice, nx)

    # acumulador
    ke_acc = np.zeros(t_common.shape[0], dtype=np.float64)

    # rangos
    y0, y1 = int(y_slice.start), int(y_slice.stop)

    # slice para W en Zl (si corresponde)
    #
    # OJO: nz_w es la dimension COMPLETA de W en la malla (p.ej. 50 en Zl),
    # mientras que nz_centers es el numero de niveles YA RECORTADOS por z_slice
    # (p.ej. 34). Compararlos entre si hacia que ambas ramas fueran falsas y W
    # se leyera sin el nivel extra: unstagger_w devolvia 33 niveles contra los
    # 34 de rho -> "operands could not be broadcast (41,34,..) (41,33,..)".
    # Ahora se pide un nivel de mas siempre que exista, y mas abajo se decide
    # si centrar segun la forma REAL leida.
    W_var = nc_state.variables[var_w]
    nz_w = int(getattr(W_var, "shape", (None, None))[1]) if len(getattr(W_var, "shape", ())) >= 2 else None
    nz_centers = int(drF.shape[0])
    z_norm = _normalize_slice(z_slice, nz_centers if nz_w is None else nz_w)
    if nz_w is not None and (z_norm.stop + 1) <= nz_w:
        z_w = _slice_plus_one(z_norm)
    else:
        z_w = z_norm

    for j0 in range(y0, y1, int(ny_block)):
        j1b = min(y1, j0 + int(ny_block))
        y_c = slice(j0, j1b)

        # métricas y hFac
        hC = np.asarray(gr.variables[hFacC_name][z_slice, y_c, x_slice], dtype=float)
        area = np.asarray(gr.variables[rA_name][y_c, x_slice], dtype=float)

        # variables
        if idx_s_mnc is None:
            U_full = np.asarray(nc_state.variables[var_u][:, z_slice, y_c, x_slice])
            V_full = np.asarray(nc_state.variables[var_v][:, z_slice, y_c, x_slice])
            W_full = np.asarray(W_var[:, z_w, y_c, x_slice])
            U_blk = U_full[idx_s, ...]
            V_blk = V_full[idx_s, ...]
            W_blk = W_full[idx_s, ...]
        else:
            U_blk = np.asarray(nc_state.variables[var_u][idx_s_mnc, z_slice, y_c, x_slice])
            V_blk = np.asarray(nc_state.variables[var_v][idx_s_mnc, z_slice, y_c, x_slice])
            W_blk = np.asarray(W_var[idx_s_mnc, z_w, y_c, x_slice])
        if idx_r_mnc is None:
            R_full = np.asarray(nc_rho.variables[var_rho][:, z_slice, y_c, x_slice])
            R_blk = R_full[idx_r, ...]
        else:
            R_blk = np.asarray(nc_rho.variables[var_rho][idx_r_mnc, z_slice, y_c, x_slice])

        if use_float32:
            U_blk = U_blk.astype(np.float32, copy=False)
            V_blk = V_blk.astype(np.float32, copy=False)
            W_blk = W_blk.astype(np.float32, copy=False)
            R_blk = R_blk.astype(np.float32, copy=False)

        # centrar en malla de densidad. Se decide por la forma REAL leida, no
        # por las dimensiones nominales (ver nota sobre nz_w / nz_centers).
        uc, vc = unstagger(U_blk, V_blk)
        nz_read = int(W_blk.shape[1])
        if nz_read == nz_centers + 1:
            wc = unstagger_w(W_blk)          # W en caras -> centros
        elif nz_read == nz_centers:
            wc = W_blk                       # W ya en centros
        else:
            raise ValueError(
                f"W tiene {nz_read} niveles; se esperaban {nz_centers} "
                f"(centros) o {nz_centers + 1} (caras)."
            )

        # ajustar a (ny-1, nx-1)
        wc = wc[..., :-1, :-1]
        rho = (rho_const + R_blk)[..., :-1, :-1]
        hC_c = hC[..., :-1, :-1]
        area_c = area[:-1, :-1]

        if drF.shape[0] != rho.shape[1]:
            raise ValueError(
                f"drF (nz={drF.shape[0]}) no coincide con rho (nz={rho.shape[1]})."
            )

        # volumen por celda
        vol_3d = (drF[:, None, None] * area_c[None, :, :]) * hC_c

        u0 = np.nan_to_num(uc, nan=0.0, copy=False)
        v0 = np.nan_to_num(vc, nan=0.0, copy=False)
        w0 = np.nan_to_num(wc, nan=0.0, copy=False)
        if rho_clip is not None:
            rho = np.where(np.abs(rho) <= rho_clip, rho, np.nan)
        r0 = np.nan_to_num(rho, nan=0.0, copy=False)

        spd2 = u0 * u0 + v0 * v0 + w0 * w0
        ke_blk = 0.5 * np.einsum("tzyx,tzyx,zyx->t", r0, spd2, vol_3d, optimize=True)
        ke_acc += ke_blk

    return t_common.astype(np.float64), ke_acc


def ke_time_series_from_state_and_rho_chunked(
    scenario_root,
    subdir: str,
    gr,
    z_slice: slice,
    y_slice: slice,
    x_slice: slice,
    state_patterns=None,
    rho_patterns=None,
    ny_block: int = 64,
    rho_const: float = 1025.0,
    tol_time: float = 1e-7,
    verbose: bool = True,
    **kw,
):
    """KE(t) con densidad completa, robusta a chunks de longitud desigual.

    ``ke_time_series_from_state_and_rho_mnc`` recibe UN solo objeto ``mnc_files``.
    Si ese objeto abarca varios chunks de iteracion con distinto numero de
    registros temporales (la corrida de jun-2026 escribe 42,42,42,42,42,31 por
    etapa), ``MITgcmutils.mnc.mnc_files.__getitem__`` revienta con
    ``could not broadcast input array from shape (31,) into shape (42,)``: asume
    que todos los chunks tienen la misma longitud en T.

    Aqui se abre UN chunk a la vez y se llama a la funcion original por cada par
    (chunk de state, chunk de diag_rho) cuyos rangos temporales se solapan. Es
    necesario porque state y diag_rho NO comparten troceado: en jun-2026 state
    tiene 12 chunks y diag_rho 4, con fronteras distintas, asi que un chunk de
    state puede cruzar dos de diag_rho.

    Devuelve (t, ke) concatenado, ordenado y sin duplicados.
    """
    scenario_root = Path(scenario_root)
    base = scenario_root / subdir
    if state_patterns is None:
        state_patterns = discover_stage_patterns(base, var_prefix="state")
    if rho_patterns is None:
        rho_patterns = discover_stage_patterns(base, var_prefix="diag_rho")

    if verbose:
        print(f"    [chunks: {len(state_patterns)} de state, "
              f"{len(rho_patterns)} de diag_rho]")

    t_parts, ke_parts = [], []
    for sp in state_patterns:
        nc_s = _open_mnc_or_raise(base, sp)
        try:
            ts = _get_time_from_nc(nc_s)
            if ts.size == 0:
                continue
            for rp in rho_patterns:
                nc_r = _open_mnc_or_raise(base, rp)
                try:
                    tr = _get_time_from_nc(nc_r)
                    if tr.size == 0:
                        continue
                    # descartar pares sin solape temporal: evita el ValueError
                    # "No hay tiempos en comun" de la funcion original
                    if tr.max() < ts.min() - tol_time or tr.min() > ts.max() + tol_time:
                        continue
                    _, _, tc = _match_time_indices(ts, tr, tol=tol_time)
                    if tc.size == 0:
                        continue
                    t_seg, ke_seg = ke_time_series_from_state_and_rho_mnc(
                        nc_s, nc_r, gr,
                        z_slice=z_slice, y_slice=y_slice, x_slice=x_slice,
                        ny_block=ny_block, rho_const=rho_const,
                        tol_time=tol_time, **kw,
                    )
                    t_parts.append(t_seg)
                    ke_parts.append(ke_seg)
                finally:
                    _close_mnc(nc_r)
        finally:
            _close_mnc(nc_s)

    if not t_parts:
        raise ValueError(
            f"No se pudo emparejar ningun chunk de state con diag_rho en {base}."
        )

    t_all = np.concatenate(t_parts)
    ke_all = np.concatenate(ke_parts)
    order = np.argsort(t_all, kind="stable")
    t_all, ke_all = t_all[order], ke_all[order]
    keep = np.concatenate(([True], np.diff(t_all) > tol_time))
    t_all, ke_all = t_all[keep], ke_all[keep]
    if verbose:
        print(f"    [KE: {t_all.size} pasos, "
              f"t=[{t_all.min()/3600:.2f},{t_all.max()/3600:.2f}] h]")
    return t_all.astype(np.float64), ke_all


def ke_time_series_from_state_and_rho_nc(
    nc_state,
    nc_rho,
    gr,
    z_slice: slice,
    y_slice: slice,
    x_slice: slice,
    ny_block: int = 64,
    var_rho: str = "rho",
    var_u: str = "U",
    var_v: str = "V",
    var_w: str = "W",
    tol_time: float = 1e-7,
    use_float32: bool = True,
    rho_clip: float = 1e6,
):
    """
    Calcula KE(t) integrando 0.5*rho*(u^2+v^2+w^2)*vol usando rho guardado en
    NetCDF (var_rho) y velocidades, centradas en la malla de densidad.
    """
    t_state = _get_time_from_nc(nc_state)
    t_rho = _get_time_from_nc(nc_rho)
    idx_s, idx_r, t_common = _match_time_indices(t_state, t_rho, tol=tol_time)
    if t_common.size == 0:
        raise ValueError("No hay tiempos en común entre state y rho.nc.")
    idx_s_mnc = _indices_to_slice(idx_s)

    # métricas
    rA_name = _pick_var(gr, ("rA", "rAc"))
    hFacC_name = _pick_var(gr, ("HFacC", "hFacC"))
    rAc = np.asarray(gr.variables[rA_name][:], dtype=float)
    drF = np.asarray(gr.variables["drF"][z_slice], dtype=float).squeeze()
    ny, nx = rAc.shape
    y_slice = _normalize_slice(y_slice, ny)
    x_slice = _normalize_slice(x_slice, nx)

    # subdominio global
    y0, y1 = int(y_slice.start), int(y_slice.stop)
    x0, x1 = int(x_slice.start), int(x_slice.stop)

    # rho subdominio (indices locales)
    rho_var = nc_rho.variables[var_rho]
    ny_rho = int(rho_var.shape[2])
    nx_rho = int(rho_var.shape[3])
    if (y1 - y0) != ny_rho or (x1 - x0) != nx_rho:
        raise ValueError(
            f"Dimensiones rho.nc (ny={ny_rho}, nx={nx_rho}) no coinciden con "
            f"subdominio global (ny={y1 - y0}, nx={x1 - x0})."
        )

    # slice para W en Zl (si corresponde)
    W_var = nc_state.variables[var_w]
    nz_w = int(getattr(W_var, "shape", (None, None))[1]) if len(getattr(W_var, "shape", ())) >= 2 else None
    nz_centers = int(drF.shape[0])
    z_w = _slice_plus_one(z_slice) if (nz_w is not None and nz_w == nz_centers + 1) else z_slice

    # acumulador
    ke_acc = np.zeros(t_common.shape[0], dtype=np.float64)

    for j0 in range(y0, y1, int(ny_block)):
        j1b = min(y1, j0 + int(ny_block))
        y_c = slice(j0, j1b)
        y_c_local = slice(j0 - y0, j1b - y0)
        x_c_local = slice(0, x1 - x0)

        # métricas y hFac
        hC = np.asarray(gr.variables[hFacC_name][z_slice, y_c, x_slice], dtype=float)
        area = np.asarray(gr.variables[rA_name][y_c, x_slice], dtype=float)

        # variables
        if idx_s_mnc is None:
            U_full = np.asarray(nc_state.variables[var_u][:, z_slice, y_c, x_slice])
            V_full = np.asarray(nc_state.variables[var_v][:, z_slice, y_c, x_slice])
            W_full = np.asarray(W_var[:, z_w, y_c, x_slice])
            U_blk = U_full[idx_s, ...]
            V_blk = V_full[idx_s, ...]
            W_blk = W_full[idx_s, ...]
        else:
            U_blk = np.asarray(nc_state.variables[var_u][idx_s_mnc, z_slice, y_c, x_slice])
            V_blk = np.asarray(nc_state.variables[var_v][idx_s_mnc, z_slice, y_c, x_slice])
            W_blk = np.asarray(W_var[idx_s_mnc, z_w, y_c, x_slice])
        R_blk = np.asarray(rho_var[idx_r, z_slice, y_c_local, x_c_local])

        if use_float32:
            U_blk = U_blk.astype(np.float32, copy=False)
            V_blk = V_blk.astype(np.float32, copy=False)
            W_blk = W_blk.astype(np.float32, copy=False)
            R_blk = R_blk.astype(np.float32, copy=False)

        # centrar en malla de densidad
        uc, vc = unstagger(U_blk, V_blk)
        if nz_w is not None and nz_w == nz_centers:
            wc = W_blk
        else:
            wc = unstagger_w(W_blk)

        # ajustar a (ny-1, nx-1)
        wc = wc[..., :-1, :-1]
        rho = R_blk[..., :-1, :-1]
        hC_c = hC[..., :-1, :-1]
        area_c = area[:-1, :-1]

        if drF.shape[0] != rho.shape[1]:
            raise ValueError(
                f"drF (nz={drF.shape[0]}) no coincide con rho (nz={rho.shape[1]})."
            )

        # volumen por celda
        vol_3d = (drF[:, None, None] * area_c[None, :, :]) * hC_c

        u0 = np.nan_to_num(uc, nan=0.0, copy=False)
        v0 = np.nan_to_num(vc, nan=0.0, copy=False)
        w0 = np.nan_to_num(wc, nan=0.0, copy=False)
        if rho_clip is not None:
            rho = np.where(np.abs(rho) <= rho_clip, rho, np.nan)
        r0 = np.nan_to_num(rho, nan=0.0, copy=False)

        spd2 = u0 * u0 + v0 * v0 + w0 * w0
        ke_blk = 0.5 * np.einsum("tzyx,tzyx,zyx->t", r0, spd2, vol_3d, optimize=True)
        ke_acc += ke_blk

    return t_common.astype(np.float64), ke_acc


def guardar_ke_rho_diff_npz(
    scenario_root: Path,
    scenario_label: str,
    bay_dir: str = "run_expand/mnc_00*",
    nobay_dir: str = "run_expand_nobay/mnc_00*",
    state_glob: str = "state.*.t*.nc",
    rho_glob: str = "diag_rho.*.t*.nc",
    grid_glob: str = "grid.t*.nc",
    z_slice: slice = slice(0, 34),
    y_slice: slice = slice(0, 323),
    x_slice: slice = slice(None),
    ny_block: int = 64,
    rho_const: float = 999.8,
    tol_time: float = 1e-7,
    save_dir: Path | None = None,
    out_name: str | None = None,
):
    """
    Calcula KE(t) usando RHOAnoma y velocidades y guarda KE_bay, KE_nobay y KE_diff.
    """
    scenario_root = Path(scenario_root)
    if save_dir is None:
        save_dir = Path(__file__).resolve().parent
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # grid (desde bay); state/diag_rho se abren chunk por chunk dentro de
    # ke_time_series_from_state_and_rho_chunked (ver su docstring: mnc_files no
    # tolera chunks de distinta longitud en T, y state/diag_rho tienen troceados
    # distintos).
    gr = None
    try:
        gr = _open_mnc_or_raise(scenario_root, f"{bay_dir}/{grid_glob}")

        # bay_dir/nobay_dir vienen como "run_expand/mnc_00*"; el descubrimiento
        # de chunks ya mira dentro de mnc_00*/, asi que basta el subdirectorio.
        sub_b = bay_dir.split("/")[0]
        sub_n = nobay_dir.split("/")[0]

        print(f"  [bay]")
        t_b, ke_b = ke_time_series_from_state_and_rho_chunked(
            scenario_root, sub_b, gr,
            z_slice=z_slice, y_slice=y_slice, x_slice=x_slice,
            ny_block=ny_block, rho_const=rho_const, tol_time=tol_time,
        )
        print(f"  [nobay]")
        t_n, ke_n = ke_time_series_from_state_and_rho_chunked(
            scenario_root, sub_n, gr,
            z_slice=z_slice, y_slice=y_slice, x_slice=x_slice,
            ny_block=ny_block, rho_const=rho_const, tol_time=tol_time,
        )
    finally:
        if gr is not None: _close_mnc(gr)

    # alinear tiempos si difieren
    if t_b.shape == t_n.shape and np.all(np.abs(t_b - t_n) <= tol_time):
        t_common = t_b
        ke_nobay = ke_n
    else:
        t_common = t_b
        ke_nobay = np.interp(t_common, t_n, ke_n)

    ke_diff = ke_b - ke_nobay

    if out_name is None:
        ztag = f"nz{z_slice.start or 0}-{z_slice.stop}"
        ytag = f"y{y_slice.start or 0}-{y_slice.stop}"
        out_name = f"{scenario_label}_KErho_diff_{ztag}_{ytag}.npz"

    out_path = save_dir / out_name
    meta = {
        "scenario": scenario_label,
        "root": str(scenario_root),
        "bay_dir": bay_dir,
        "nobay_dir": nobay_dir,
        "state_glob": state_glob,
        "rho_glob": rho_glob,
        "grid_glob": grid_glob,
        "z_slice": [z_slice.start, z_slice.stop, z_slice.step],
        "y_slice": [y_slice.start, y_slice.stop, y_slice.step],
        "x_slice": [x_slice.start, x_slice.stop, x_slice.step],
        "rho_const": float(rho_const),
        "tol_time": float(tol_time),
        "units": {"t": "s", "ke": "J"},
        "note": "KE usa RHOAnoma + rho_const y velocidades centradas con unstagger().",
    }

    np.savez_compressed(
        out_path,
        t_seconds=t_common.astype(np.float64),
        t_horas=(t_common / 3600.0).astype(np.float64),
        ke_bay=ke_b.astype(np.float64),
        ke_nobay=ke_nobay.astype(np.float64),
        ke_diff=ke_diff.astype(np.float64),
        meta_json=np.array(json.dumps(meta, ensure_ascii=False)),
    )
    print(f"[OK] Guardado KE rho diff: {out_path}")
    return out_path


def procesar_escenarios_ke_rho(
    scenarios_dict: dict,
    bay_dir: str = "run_expand/mnc_00*",
    nobay_dir: str = "run_expand_nobay/mnc_00*",
    state_glob: str = "state.*.t*.nc",
    rho_glob: str = "diag_rho.*.t*.nc",
    grid_glob: str = "grid.t*.nc",
    z_slice: slice = slice(0, 34),
    y_slice: slice = slice(0, 323),
    x_slice: slice = slice(None),
    ny_block: int = 64,
    rho_const: float = 999.8,
    tol_time: float = 1e-7,
    save_dir: Path = Path("./npz_ke_rho"),
):
    for scen_label, base_path in scenarios_dict.items():
        print(f"\n=== Escenario (KE rho): {scen_label} ===")
        guardar_ke_rho_diff_npz(
            scenario_root=Path(base_path),
            scenario_label=scen_label,
            bay_dir=bay_dir,
            nobay_dir=nobay_dir,
            state_glob=state_glob,
            rho_glob=rho_glob,
            grid_glob=grid_glob,
            z_slice=z_slice,
            y_slice=y_slice,
            x_slice=x_slice,
            ny_block=ny_block,
            rho_const=rho_const,
            tol_time=tol_time,
            save_dir=save_dir,
        )
    print("\n[FIN] KE (con rho) para todos los escenarios.")

def _barotropic_full_column(A_tk, drF, hfac_k, nan_ok=True):
    """
    Promedio barotrópico ponderado por espesor y hFac en TODA la columna.
    - A_tk: (nt, nz)
    - drF: (nz,)
    - hfac_k: (nz,) para la columna (j,i) elegida
    Devuelve: (nt,)
    """
    w = drF * np.clip(hfac_k, 0.0, 1.0)  # pesos por nivel
    if not np.any(w > 0):
        return np.full(A_tk.shape[0], np.nan)
    Wsum = np.sum(w)
    # broadcast a (nt, nz)
    if nan_ok:
        num = np.nansum(A_tk * w[None, :], axis=1)
        den = np.nansum(np.where(np.isfinite(A_tk), w[None, :], 0.0), axis=1)
        den = np.where(den > 0, den, np.nan)
        return num / den
    else:
        return np.sum(A_tk * w[None, :], axis=1) / Wsum
    
def guardar_barotropico_baroclinico_escenario_fullcol(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    save_dir: Path = Path("./npz_baroclinico"),
    varname: str = "U",                # 'U' o 'W'
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    depths_m=(110.0, 250.0),           # solo para recortar la salida baroclínica/figura
    stage_patterns=("state.0000000000.t*.nc", "state.0000007200.t*.nc")
):
    """
    Calcula y guarda:
      - barotrópica(t) = promedio vertical EN TODA LA COLUMNA (ponderado)
      - baroclínica(t,k) = A(t,k) - barotrópica(t)
    Genera un .npz por estación y profundidad de interés (para visualizar por depth),
    PERO la barotrópica se calcula siempre con la columna completa.
    """
    data = cargar_datos_mitgcm_etapas(outdir_bay, outdir_nobay,
                                      stage_patterns=stage_patterns)  # gr, listas, time...
    gr = data["gr_bay"]
    stations = make_stations_from_grid(gr, j_centro=j_centro, radio=radio)  # [(i,j,label,color), ...]

    # métricas verticales
    RC, drF, hfac_fn = _get_vertical_metrics(gr, varname)

    # carpeta de salida
    outdir = (Path(save_dir) / scenario_label)
    outdir.mkdir(parents=True, exist_ok=True)

    # recorre estaciones
    for i, j, label, *_ in stations:
        times, perfiles = [], []
        for nc in data["nc_bay_list"]:
            t = _get_time_from_nc(nc)
            A = np.array(nc.variables[varname][:, :, j, i], dtype=float)  # (nt, nz)
            times.append(t); perfiles.append(A)

        if not times:
            continue

        T, Aall = _concat_time_and_data(times, perfiles)  # (nt,), (nt,nz)

        # barotrópica full columna (ponderada)
        hfac_k = hfac_fn(j, i)           # (nz,)
        barotrop_full = _barotropic_full_column(Aall, drF, hfac_k, nan_ok=True)  # (nt,)

        # baroclínica en toda la columna (para recortes por depth luego)
        BC_all = Aall - barotrop_full[:, None]  # (nt, nz)

        # guardar un archivo por profundidad (solo recorte vertical para gráficos/consistencia)
        for depth in depths_m:
            # subset niveles hasta esa profundidad (pero barotrópica ya es full col)
            k_sel = np.where(np.abs(RC) <= float(depth))[0]
            if k_sel.size == 0:
                continue

            outpath = outdir / f"{varname}_{label}_{int(depth)}m_baroclinico.npz"
            np.savez(
                outpath,
                t=T,            # tiempo (como venía)
                RC=RC[k_sel],   # niveles para la figura
                barotropica=barotrop_full,       # misma (full col) para todos
                baroclinica=BC_all[:, k_sel]     # recorte para esa página/figura
            )
            print(f"[OK] {scenario_label} • {varname} • {label} • {int(depth)}m  → {outpath}")

def procesar_varios_escenarios_baro_fullcol(
    scenarios: dict,
    subdir_bay: str = "run_expand/mnc_00*",
    subdir_nobay: str = "run_expand_nobay/mnc_00*",
    save_dir: Path = Path("./npz_baroclinico"),
    varnames=("U","W"),
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    depths_m=(110.0, 250.0),
    stage_patterns=("state.0000000000.t*.nc", "state.0000007200.t*.nc")
):
    for scen_label, base in scenarios.items():
        print(f"\n=== Escenario: {scen_label} ===")
        outdir_bay   = Path(base) / subdir_bay
        outdir_nobay = Path(base) / subdir_nobay
        for vn in varnames:
            guardar_barotropico_baroclinico_escenario_fullcol(
                scenario_label = scen_label,
                outdir_bay     = outdir_bay,
                outdir_nobay   = outdir_nobay,
                save_dir       = save_dir,
                varname        = vn,
                j_centro       = j_centro,
                radio          = radio,
                depths_m       = depths_m,
                stage_patterns = stage_patterns
            )

# --------------------------------------------------------------------
# Funcion para obtener las anomalias de temperatura pontencial y guardar
# en un archivo .npz
# --------------------------------------------------------------------

# ─────────────────────────── Perfil de referencia por escenario ───────────────────────────
def _tref_profile_for_scenario(scenario_label: str) -> np.ndarray:
    scen = scenario_label.strip().lower()
    if scen == "linear":
        tRef = [
            33.6, 33.5, 33.5, 33.4, 33.3, 33.2, 33.2, 33.1, 33.0,
            32.9, 32.8, 32.7, 32.5, 32.4, 32.2, 32.1, 31.9,
            31.7, 31.5, 31.3, 31.0, 30.8, 30.5, 30.2, 29.9,
            29.5, 29.1, 28.7, 28.3, 27.8, 27.3, 26.7, 26.1,
            25.4, 24.7, 23.9, 23.1, 22.2, 21.2, 20.2, 19.0,
            17.8, 16.4, 15.0, 13.4, 11.7, 9.9, 7.9, 5.8,
            3.5,
        ]
    elif scen == "realistic":
        tRef = [
            16.6, 16.6, 16.6, 16.6, 16.6, 16.6, 16.5, 16.5, 16.5,
            16.5, 16.5, 16.5, 16.5, 16.5, 16.5, 16.5, 16.4,
            16.4, 16.4, 16.3, 16.2, 16.1, 16.0, 15.8, 15.6,
            15.4, 15.1, 14.7, 14.2, 13.6, 13.0, 12.3, 11.7,
            11.1, 10.6, 10.2, 9.7, 9.3, 8.9, 8.6, 8.3,
            8.0, 7.7, 7.3, 7.0, 6.7, 6.5, 6.2, 6.0,
            5.7,
        ]
    elif scen == "2layer" or scen == "2-layer" or scen == "2 layer":
        tRef = [
            11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4,
            11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4,
            11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4,
            11.4, 11.4, 11.4, 11.4, 11.4, 7.1, 7.1, 7.1,
            7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1,
            7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1,
            7.1,
        ]
    else:
        raise ValueError(f"No tengo tRef para escenario '{scenario_label}'.")
    return np.asarray(tRef, dtype=float)

def _match_tref_to_nz(tRef: np.ndarray, nz: int) -> np.ndarray:
    """
    Re-mapea tRef (por índice) a longitud nz si difiere.
    No asume profundidades absolutas; solo “estira” o “comprime” por niveles.
    """
    if tRef.size == nz:
        return tRef.copy()
    xi_src = np.linspace(0.0, 1.0, tRef.size)
    xi_dst = np.linspace(0.0, 1.0, nz)
    return np.interp(xi_dst, xi_src, tRef)

# ───────────────────── Función principal: anomalías de Temp (Theta) ─────────────────────
def guardar_theta_anoma_npz_escenario_unico(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    varname: str = 'Temp',                 # 3D: (t, z, y, x)
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern   = "grid.t*.nc",
    source: str = "diff"                   # "diff" (bay−nobay), "bay", "nobay"
):
    """
    Construye y guarda EN UN SOLO NPZ el Hovmöller (tiempo × profundidad × estación)
    de la **anomalía** de temperatura potencial respecto a un perfil tRef por escenario.
    - Si source="diff": usa (bay − nobay) antes de restar tRef.
    - Si source="bay"/"nobay": usa ese campo antes de restar tRef.

    Salida: <save_dir>/<scenario_label>_hovThetaAnom_<source>_<varname>.npz
    Contiene:
      - t_horas: (nt,)
      - RC: (nz,)
      - anom: (nt, nz, nstations)
      - labels: (nstations,)   (unicode)
      - colors: (nstations,)   (unicode; "" si no hay color)
      - ij: (nstations, 2)     int32 con (i, j) por estación
      - idx_100, idx_400: int
      - meta_json: json con info
    """
    save_dir.mkdir(parents=True, exist_ok=True)

    # A) Carga multi-etapas (usa tus helpers existentes)
    # Debes tener estas funcs en tu módulo:
    # - cargar_datos_mitgcm_etapas(outdir_bay, outdir_nobay, stage_patterns, grid_pattern)
    # - _get_time_from_nc(nc)
    # - _concat_time_and_data(times_list, blocks_list)
    data = cargar_datos_mitgcm_etapas(
        outdir_bay, outdir_nobay,
        stage_patterns=stage_patterns,
        grid_pattern=grid_pattern
    )
    nc_bay_list   = data["nc_bay_list"]
    nc_nobay_list = data["nc_nobay_list"]
    gr_bay        = data["gr_bay"]
    time_concat   = data["time"]  # segundos

    # B) Estaciones (usa tu make_stations... que devuelva meta y/o plot list)
    # Asegúrate de que el retorno sea ((i,j,label)[, color]) o ajusta aquí.
    stations = make_stations_from_grid(gr_bay, j_centro=j_centro, radio=radio)

    # C) Eje vertical
    RC = _get_RC(gr_bay)
    nz = int(RC.shape[0])

    # D) Recorte temporal
    t_all = np.array(time_concat, dtype=float)
    tmin_seg = tmin_horas * 3600.0
    idx0 = np.where(t_all >= tmin_seg)[0][0] if np.any(t_all >= tmin_seg) else 0
    t_cut = t_all[idx0:]
    t_h   = t_cut / 3600.0
    nt = int(t_h.shape[0])

    # E) Índices guía (≈100 y ≈400 m)
    idx_100 = int(np.argmin(np.abs(RC + 100.0)))
    idx_400 = int(np.argmin(np.abs(RC + 400.0)))

    # F) Perfil de referencia para el escenario, re-mapeado a nz
    tref_raw = _tref_profile_for_scenario(scenario_label)
    tref = _match_tref_to_nz(tref_raw, nz)  # (nz,)
    tref_row = tref[None, :]                # (1, nz) para broadcast

    # G) Construcción del cubo anom (nt × nz × nstations)
    nst = len(stations)
    anom_all = np.empty((nt, nz, nst), dtype=np.float32)
    labels = []
    colors = []
    ij     = np.empty((nst, 2), dtype=np.int32)

    for k, st in enumerate(stations):
        i, j, label, color = st

        # Recolecta por etapas (tiempos y bloques nt×nz)
        times_b, blocks_b = [], []
        times_n, blocks_n = [], []
        for nc_b, nc_n in zip(nc_bay_list, nc_nobay_list):
            tb = _get_time_from_nc(nc_b)
            tn = _get_time_from_nc(nc_n)

            Xb = np.array(nc_b.variables[varname][:, :, j, i], dtype=float)  # (t_stage, nz)
            Xn = np.array(nc_n.variables[varname][:, :, j, i], dtype=float)  # (t_stage, nz)
            if Xb.shape[0] != tb.shape[0] or Xn.shape[0] != tn.shape[0]:
                raise ValueError("Longitudes tiempo/serie no coinciden en alguna etapa.")

            times_b.append(tb); blocks_b.append(Xb)
            times_n.append(tn); blocks_n.append(Xn)

        t_b, Xb_all = _concat_time_and_data(times_b, blocks_b)  # (t_full, nz)
        t_n, Xn_all = _concat_time_and_data(times_n, blocks_n)  # (t_full, nz)

        # Seleccionar fuente
        if source == "diff":
            Xsel = Xb_all - Xn_all
        elif source == "bay":
            Xsel = Xb_all
        elif source == "nobay":
            Xsel = Xn_all
        else:
            raise ValueError("source debe ser 'diff' | 'bay' | 'nobay'.")

        # Recorte temporal y anomalía respecto a tRef
        Xsel_cut = Xsel[idx0:, :]                      # (nt, nz)
        Xanom = Xsel_cut #- tref_row                    # (nt, nz)

        if Xanom.shape != (nt, nz):
            raise ValueError(f"Forma inesperada en {label}: {Xanom.shape}, esperado {(nt, nz)}.")

        anom_all[:, :, k] = Xanom.astype(np.float32)
        labels.append(str(label))
        colors.append("" if (color is None) else str(color))
        ij[k, 0] = int(i)
        ij[k, 1] = int(j)

    # Cerrar archivos NetCDF una sola vez, tras procesar TODAS las estaciones.
    # (Cerrar dentro del bucle dejaba handles HDF5 cerrados y provocaba segfault
    #  al leer la siguiente estación.)
    for nc in nc_bay_list:
        nc.close()
    for nc in nc_nobay_list:
        nc.close()

    # H) Metadatos
    meta = {
        "scenario": scenario_label,
        "var": varname,
        "paths": {"bay": str(outdir_bay), "nobay": str(outdir_nobay)},
        "stage_patterns": list(stage_patterns),
        "tmin_horas": float(tmin_horas),
        "idx0": int(idx0),
        "dt_seconds": float(t_cut[1] - t_cut[0]) if len(t_cut) > 1 else None,
        "source": source,
        "tRef_len_raw": int(tref_raw.size),
        "stations": [{"label": labels[k], "i": int(ij[k,0]), "j": int(ij[k,1]), "color": colors[k]} for k in range(nst)]
    }
    meta_json = json.dumps(meta, ensure_ascii=False)

    # I) Guardado
    out_path = save_dir / f"{scenario_label}_hovThetaAnom_{source}_{varname}.npz"
    np.savez_compressed(
        out_path,
        t_horas=t_h.astype(np.float64),
        RC=RC.astype(np.float64),
        anom=anom_all,                              # (nt, nz, nstations)
        labels=np.array(labels, dtype="U32"),
        colors=np.array(colors, dtype="U16"),
        ij=ij.astype(np.int32),
        idx_100=np.array(idx_100, dtype=np.int32),
        idx_400=np.array(idx_400, dtype=np.int32),
        meta_json=np.array(meta_json),
        tRef=tref.astype(np.float64)               # guardo el perfil usado
    )
    print(f"[OK] Bundle Hovmöller (anom. {varname}) guardado: {out_path}")
    return str(out_path)

#

#def guardar_theta_anoma__


def plot_theta_anom_hov_from_npz(
    npz_path: str,
    station: str | int = "A",
    vmin: float | None = None,
    vmax: float | None = None,
    cmap = "RdBu_r",
    title_prefix: str = "",
    symmetric: bool = False          # fuerza vmin=-m, vmax=+m
):
    d = np.load(npz_path, allow_pickle=True)
    t   = d["t_horas"]       # (nt,)
    RC  = d["RC"]            # (nz,)
    A3  = d["anom"]          # (nt,nz,nst)
    labs= [str(x) for x in d["labels"]]
    meta= json.loads(str(d["meta_json"]))

    # Selección de estación
    if isinstance(station, int):
        k = station; label = labs[k]
    else:
        label = str(station)
        if label not in labs:
            raise KeyError(f"Estación '{label}' no está en {labs}")
        k = labs.index(label)

    Z = A3[:, :, k]  # (nt, nz)

    # Límites de color
    if vmin is None or vmax is None:
        m = float(np.nanmax(np.abs(Z)))
        if m == 0:
            vmin_auto, vmax_auto = -1e-3, 1e-3
        else:
            vmin_auto, vmax_auto = (-m, m) if symmetric else (float(np.nanmin(Z)), float(np.nanmax(Z)))
        vmin = vmin if vmin is not None else vmin_auto
        vmax = vmax if vmax is not None else vmax_auto
    elif symmetric:
        # si el usuario fija vmin/vmax y symmetric=True, respetamos el mayor |.| y hacemos simétrico
        M = max(abs(vmin), abs(vmax))
        vmin, vmax = -M, M

    fig, ax = plt.subplots(figsize=(9,5))
    pcm = ax.pcolormesh(t, RC, Z.T, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    cbar = fig.colorbar(pcm, ax=ax, label=f"Anomalía {meta.get('var','Temp')} (°C)")
    ax.set_xlabel("Tiempo (h)")
    ax.set_ylabel("Profundidad (m)")
    ax.set_title(f"{title_prefix} Anomalía {meta.get('var','Temp')} — {meta.get('scenario','')} • est. {label} • fuente={meta.get('source','')}")
    plt.tight_layout()
    plt.show()
    return fig, ax

def _resolve_vertical_mode_npz_for_scenario(scenario: str, modes_dir: str | Path):
    """
    Encuentra el archivo NPZ de modos verticales para un escenario.
    Prioriza patrones explícitos y luego fallback por glob.
    """
    mdir = Path(modes_dir)
    if not mdir.exists():
        raise FileNotFoundError(f"No existe modes_dir: {mdir}")

    scenario = str(scenario).strip().lower()
    preferred = [
        mdir / f"{scenario}_modes_verticales_rho_xy30_30_t0.npz",
        mdir / f"{scenario}_mode_2layer_CR2005.npz",
    ]
    for p in preferred:
        if p.exists():
            return str(p)

    candidates = sorted(mdir.glob(f"{scenario}*.npz"))
    if not candidates:
        raise FileNotFoundError(
            f"No encontré NPZ de modos para escenario='{scenario}' en {mdir}"
        )
    return str(candidates[0])

def plot_theta_anom_hov_with_vertical_mode_from_npz(
    npz_path: str,
    station: str | int = "A",
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap = "RdBu_r",
    title_prefix: str = "",
    symmetric: bool = False,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.5"
):
    """
    Igual que plot_theta_anom_hov_from_npz, pero agrega a la izquierda
    el modo vertical del escenario (desde NPZ en modes_dir).

    n_modes_plot: cuántos modos mostrar (solo si el NPZ contiene 'modes').
    mode_index: fallback 1-based para mostrar un solo modo cuando n_modes_plot=1.
    En 2layer se usa 'phi' (modo único).
    """
    d = np.load(npz_path, allow_pickle=True)
    t   = d["t_horas"]       # (nt,)
    RC  = d["RC"]            # (nz,)
    A3  = d["anom"]          # (nt,nz,nst)
    labs= [str(x) for x in d["labels"]]
    meta= json.loads(str(d["meta_json"]))

    if isinstance(station, int):
        k = station
        label = labs[k]
    else:
        label = str(station)
        if label not in labs:
            raise KeyError(f"Estación '{label}' no está en {labs}")
        k = labs.index(label)

    Z = A3[:, :, k]  # (nt, nz)

    # Normalización por estación usando max(|anom|) de todo el Hovmöller de esa estación.
    z_scale = 1.0
    if normalize_station_max:
        z_scale = float(np.nanmax(np.abs(Z)))
        if (not np.isfinite(z_scale)) or z_scale <= 0.0:
            z_scale = 1.0
        Z = Z / z_scale

    if vmin is None or vmax is None:
        if normalize_station_max:
            # Escala fija para campo normalizado.
            vmin_auto, vmax_auto = -1.0, 1.0
        else:
            m = float(np.nanmax(np.abs(Z)))
            if m == 0:
                vmin_auto, vmax_auto = -1e-3, 1e-3
            else:
                if symmetric:
                    vmin_auto, vmax_auto = -m, m
                else:
                    vmin_auto, vmax_auto = float(np.nanmin(Z)), float(np.nanmax(Z))
        vmin = vmin if vmin is not None else vmin_auto
        vmax = vmax if vmax is not None else vmax_auto
    elif symmetric:
        M = max(abs(vmin), abs(vmax))
        vmin, vmax = -M, M

    scenario = str(meta.get("scenario", "")).strip()
    if mode_npz_path is None:
        mode_npz_path = _resolve_vertical_mode_npz_for_scenario(scenario, modes_dir)

    md = np.load(mode_npz_path, allow_pickle=True)
    mode_depth = np.array(md["depth"], dtype=float) if "depth" in md.files else np.array(RC, dtype=float)
    # Homogeneizar signo del eje vertical con RC del Hovmöller.
    if np.nanmedian(mode_depth) > 0 and np.nanmedian(RC) < 0:
        mode_depth = -mode_depth

    def _interior_zero_depths(profile, depth):
        """
        Devuelve profundidades (z) donde profile cruza 0 en el interior.
        Excluye explícitamente superficie y fondo.
        """
        p = np.asarray(profile, dtype=float)
        z = np.asarray(depth, dtype=float)
        n = p.size
        if n < 3:
            return np.array([], dtype=float)

        z0s = []
        for i in range(1, n - 2):  # pares interiores: (i, i+1), sin tocar 0 ni n-1
            p0, p1 = p[i], p[i + 1]
            z0, z1 = z[i], z[i + 1]
            if not (np.isfinite(p0) and np.isfinite(p1) and np.isfinite(z0) and np.isfinite(z1)):
                continue
            if p0 == 0.0:
                z0s.append(float(z0))
            if p0 * p1 < 0.0:
                alpha = -p0 / (p1 - p0)
                zc = z0 + alpha * (z1 - z0)
                z0s.append(float(zc))
            elif p1 == 0.0:
                z0s.append(float(z1))

        if not z0s:
            return np.array([], dtype=float)

        z0s = np.array(sorted(z0s), dtype=float)
        uniq = [z0s[0]]
        for zz in z0s[1:]:
            if abs(zz - uniq[-1]) > 1e-6:
                uniq.append(zz)
        return np.array(uniq, dtype=float)

    mode_curves = []
    zero_depths_all = []

    if "modes" in md.files:
        modes = np.array(md["modes"], dtype=float)
        c = np.array(md["c"], dtype=float) if "c" in md.files else np.array([])
        nmax = modes.shape[1]
        nshow = max(1, int(n_modes_plot))
        if nshow == 1:
            m_indices = [min(max(int(mode_index) - 1, 0), nmax - 1)]
        else:
            m_indices = list(range(min(nshow, nmax)))

        for m0 in m_indices:
            prof = modes[:, m0]
            lab = f"modo {m0 + 1}"
            if c.size > m0 and np.isfinite(c[m0]):
                lab += f" (c={c[m0]:.3f} m/s)"
            mode_curves.append((prof, lab))
            zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    elif "phi" in md.files:
        prof = np.array(md["phi"], dtype=float)
        c = float(np.array(md["c"]).squeeze()) if "c" in md.files else np.nan
        mode_label = "modo 2layer"
        if np.isfinite(c):
            mode_label += f" (c={c:.3f} m/s)"
        mode_curves.append((prof, mode_label))
        zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    else:
        raise KeyError(f"{mode_npz_path} no contiene 'modes' ni 'phi'.")

    fig, (ax_mode, ax_hov) = plt.subplots(
        1, 2, figsize=(11, 5),
        gridspec_kw={"width_ratios": [1.2, 4.8], "wspace": 0.08},
        sharey=True
    )

    # Panel izquierdo: modos verticales
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(mode_curves), 2)))
    for kk, (prof, lab) in enumerate(mode_curves):
        ax_mode.plot(prof, mode_depth, lw=1.7, color=colors[kk], label=lab)
    ax_mode.axvline(0.0, color="0.5", lw=0.8, ls="--")

    if zero_depths_all:
        zarr = np.array(sorted(zero_depths_all), dtype=float)
        zuniq = [zarr[0]]
        for zz in zarr[1:]:
            if abs(zz - zuniq[-1]) > 1e-6:
                zuniq.append(zz)
        for zz in zuniq:
            ax_mode.axhline(zz, color="0.35", lw=0.9, ls=":")

    ax_mode.set_xlabel("Modo")
    ax_mode.set_ylabel("Profundidad (m)")
    ax_mode.set_title(f"{scenario}")
    ax_mode.grid(True, alpha=0.25)
    ax_mode.legend(fontsize=8, loc="best")

    # Panel derecho: Hovmöller (contourf)
    levels = np.linspace(vmin, vmax, 61)
    pcm = ax_hov.contourf(t, RC, Z.T, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax, extend="both")
    if normalize_station_max:
        cbar_label = f"Anomalía {meta.get('var','Temp')} / max|anom| (est. {label})"
    else:
        cbar_label = f"Anomalía {meta.get('var','Temp')} (°C)"
    fig.colorbar(pcm, ax=ax_hov, label=cbar_label)
    ax_hov.set_xlabel("Tiempo (h)")
    if draw_positive_dashed_contours:
        # Contorno punteado del borde de la región positiva (nivel 0).
        Zmask = np.ma.masked_less_equal(Z.T, 0.0)
        ax_hov.contour(
            t, RC, Zmask,
            levels=[0.0],
            colors=positive_contour_color,
            linewidths=positive_contour_lw,
            linestyles="--"
        )
    if zero_depths_all:
        zarr = np.array(sorted(zero_depths_all), dtype=float)
        zuniq = [zarr[0]]
        for zz in zarr[1:]:
            if abs(zz - zuniq[-1]) > 1e-6:
                zuniq.append(zz)
        for zz in zuniq:
            ax_hov.axhline(zz, color="w", lw=0.9, ls=":", alpha=0.75)
    ax_hov.set_title(
        f"{title_prefix} Anomalía {meta.get('var','Temp')} — {scenario} • est. {label} • fuente={meta.get('source','')}"
    )

    plt.tight_layout()
    plt.show()
    return fig, (ax_mode, ax_hov)

def plot_theta_anom_hov_with_vertical_mode_multi_station_from_npz(
    npz_path: str,
    stations: list[str | int],
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap = "RdBu_r",
    title_prefix: str = "",
    symmetric: bool = False,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.5",
    fig_width: float = 12.0,
    row_height: float = 3.2
):
    """
    Figura única para varias estaciones (una fila por estación):
      - Columna izquierda: primeros modos verticales del escenario.
      - Columna derecha : Hovmöller de anomalía para la estación.
    """
    if not stations:
        raise ValueError("stations no puede estar vacío.")

    d = np.load(npz_path, allow_pickle=True)
    t   = d["t_horas"]       # (nt,)
    RC  = d["RC"]            # (nz,)
    A3  = d["anom"]          # (nt,nz,nst)
    labs= [str(x) for x in d["labels"]]
    meta= json.loads(str(d["meta_json"]))
    scenario = str(meta.get("scenario", "")).strip()

    if mode_npz_path is None:
        mode_npz_path = _resolve_vertical_mode_npz_for_scenario(scenario, modes_dir)
    md = np.load(mode_npz_path, allow_pickle=True)

    mode_depth = np.array(md["depth"], dtype=float) if "depth" in md.files else np.array(RC, dtype=float)
    if np.nanmedian(mode_depth) > 0 and np.nanmedian(RC) < 0:
        mode_depth = -mode_depth

    def _interior_zero_depths(profile, depth):
        p = np.asarray(profile, dtype=float)
        z = np.asarray(depth, dtype=float)
        n = p.size
        if n < 3:
            return np.array([], dtype=float)
        z0s = []
        for i in range(1, n - 2):
            p0, p1 = p[i], p[i + 1]
            z0, z1 = z[i], z[i + 1]
            if not (np.isfinite(p0) and np.isfinite(p1) and np.isfinite(z0) and np.isfinite(z1)):
                continue
            if p0 == 0.0:
                z0s.append(float(z0))
            if p0 * p1 < 0.0:
                alpha = -p0 / (p1 - p0)
                z0s.append(float(z0 + alpha * (z1 - z0)))
            elif p1 == 0.0:
                z0s.append(float(z1))
        if not z0s:
            return np.array([], dtype=float)
        z0s = np.array(sorted(z0s), dtype=float)
        uniq = [z0s[0]]
        for zz in z0s[1:]:
            if abs(zz - uniq[-1]) > 1e-6:
                uniq.append(zz)
        return np.array(uniq, dtype=float)

    mode_curves = []
    zero_depths_all = []
    if "modes" in md.files:
        modes = np.array(md["modes"], dtype=float)
        c = np.array(md["c"], dtype=float) if "c" in md.files else np.array([])
        nmax = modes.shape[1]
        nshow = max(1, int(n_modes_plot))
        if nshow == 1:
            m_indices = [min(max(int(mode_index) - 1, 0), nmax - 1)]
        else:
            m_indices = list(range(min(nshow, nmax)))
        for m0 in m_indices:
            prof = modes[:, m0]
            lab = f"modo {m0 + 1}"
            if c.size > m0 and np.isfinite(c[m0]):
                lab += f" (c={c[m0]:.3f} m/s)"
            mode_curves.append((prof, lab))
            zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    elif "phi" in md.files:
        prof = np.array(md["phi"], dtype=float)
        c = float(np.array(md["c"]).squeeze()) if "c" in md.files else np.nan
        mode_label = "modo 2layer"
        if np.isfinite(c):
            mode_label += f" (c={c:.3f} m/s)"
        mode_curves.append((prof, mode_label))
        zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    else:
        raise KeyError(f"{mode_npz_path} no contiene 'modes' ni 'phi'.")

    zuniq = []
    if zero_depths_all:
        zarr = np.array(sorted(zero_depths_all), dtype=float)
        zuniq = [zarr[0]]
        for zz in zarr[1:]:
            if abs(zz - zuniq[-1]) > 1e-6:
                zuniq.append(zz)

    idx_and_labels = []
    for st in stations:
        if isinstance(st, int):
            k = int(st)
            if (k < 0) or (k >= len(labs)):
                raise IndexError(f"station index {k} fuera de rango [0,{len(labs)-1}]")
            idx_and_labels.append((k, labs[k]))
        else:
            lab = str(st)
            if lab not in labs:
                raise KeyError(f"Estación '{lab}' no está en {labs}")
            idx_and_labels.append((labs.index(lab), lab))

    nrows = len(idx_and_labels)
    fig, axes = plt.subplots(
        nrows, 2,
        figsize=(fig_width, row_height * nrows),
        gridspec_kw={"width_ratios": [1.2, 4.8], "wspace": 0.08, "hspace": 0.22},
        sharex="col",
        sharey="all",
        squeeze=False
    )

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(mode_curves), 2)))
    for r, (k, label) in enumerate(idx_and_labels):
        ax_mode = axes[r, 0]
        ax_hov  = axes[r, 1]

        for kk, (prof, lab) in enumerate(mode_curves):
            ax_mode.plot(prof, mode_depth, lw=1.5, color=colors[kk], label=lab)
        ax_mode.axvline(0.0, color="0.5", lw=0.8, ls="--")
        for zz in zuniq:
            ax_mode.axhline(zz, color="0.35", lw=0.9, ls=":")
        ax_mode.set_ylabel("Profundidad (m)")
        if r == nrows - 1:
            ax_mode.set_xlabel("Modo")
        ax_mode.grid(True, alpha=0.25)
        if r == 0:
            ax_mode.set_title(f"{scenario} • modos")
            ax_mode.legend(fontsize=8, loc="best")

        Z = A3[:, :, k].astype(float)
        if normalize_station_max:
            z_scale = float(np.nanmax(np.abs(Z)))
            if (not np.isfinite(z_scale)) or z_scale <= 0.0:
                z_scale = 1.0
            Z = Z / z_scale

        if vmin is None or vmax is None:
            if normalize_station_max:
                vmin_i, vmax_i = -1.0, 1.0
            else:
                m = float(np.nanmax(np.abs(Z)))
                if m == 0:
                    vmin_i, vmax_i = -1e-3, 1e-3
                else:
                    if symmetric:
                        vmin_i, vmax_i = -m, m
                    else:
                        vmin_i, vmax_i = float(np.nanmin(Z)), float(np.nanmax(Z))
            vvmin = vmin if vmin is not None else vmin_i
            vvmax = vmax if vmax is not None else vmax_i
        else:
            vvmin, vvmax = float(vmin), float(vmax)
            if symmetric:
                M = max(abs(vvmin), abs(vvmax))
                vvmin, vvmax = -M, M

        levels_i = np.linspace(vvmin, vvmax, 61)
        pcm = ax_hov.contourf(
            t, RC, Z.T,
            levels=levels_i,
            cmap=cmap,
            vmin=vvmin,
            vmax=vvmax,
            extend="both"
        )
        if draw_positive_dashed_contours:
            Zmask = np.ma.masked_less_equal(Z.T, 0.0)
            ax_hov.contour(
                t, RC, Zmask, levels=[0.0],
                colors=positive_contour_color,
                linewidths=positive_contour_lw,
                linestyles="--"
            )
        for zz in zuniq:
            ax_hov.axhline(zz, color="w", lw=0.9, ls=":", alpha=0.75)

        if r == nrows - 1:
            ax_hov.set_xlabel("Tiempo (h)")
        ax_hov.set_ylabel("Profundidad (m)")
        ax_hov.set_title(f"{title_prefix} est. {label}")

        if normalize_station_max:
            cbar_label = f"Anomalía {meta.get('var','Temp')} / max|anom|"
        else:
            cbar_label = f"Anomalía {meta.get('var','Temp')} (°C)"
        fig.colorbar(pcm, ax=ax_hov, label=cbar_label, shrink=0.98, pad=0.02)

    fig.suptitle(f"{scenario} • {meta.get('source','')} • Hovmöller por estaciones", y=0.995)
    plt.tight_layout()
    plt.show()
    return fig, axes

def plot_theta_anom_hov_with_vertical_mode_multi_station_balance_from_npz(
    npz_path: str,
    stations: list[str | int],
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin: float | None = None,
    vmax: float | None = None,
    title_prefix: str = "",
    symmetric: bool = False,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.5",
    fig_width: float = 12.0,
    row_height: float = 3.2
):
    """
    Wrapper de conveniencia:
      modos verticales + Hovmöller de anomalía de temperatura por estación,
      usando cmocean.balance como colormap por defecto.
    """
    return plot_theta_anom_hov_with_vertical_mode_multi_station_from_npz(
        npz_path=npz_path,
        stations=stations,
        mode_npz_path=mode_npz_path,
        modes_dir=modes_dir,
        mode_index=mode_index,
        n_modes_plot=n_modes_plot,
        vmin=vmin,
        vmax=vmax,
        cmap=cmo.cm.balance,
        title_prefix=title_prefix,
        symmetric=symmetric,
        normalize_station_max=normalize_station_max,
        draw_positive_dashed_contours=draw_positive_dashed_contours,
        positive_contour_lw=positive_contour_lw,
        positive_contour_color=positive_contour_color,
        fig_width=fig_width,
        row_height=row_height
    )

def plot_theta_anom_hov_with_vertical_mode_balance_across_scenarios(
    scenario_npz_map: dict[str, str],
    stations: list[str | int],
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin: float | None = None,
    vmax: float | None = None,
    symmetric: bool = False,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.5",
    fig_width: float = 12.0,
    row_height: float = 3.2
):
    """
    Ejecuta la figura de modos + Hovmöller de anomalía de Temp para varios escenarios.
    `scenario_npz_map` ejemplo:
      {
        "2layer": "npz_out/2layer_hovThetaAnom_diff_Temp.npz",
        "linear": "npz_out/linear_hovThetaAnom_diff_Temp.npz",
        "realistic": "npz_out/realistic_hovThetaAnom_diff_Temp.npz",
      }
    Retorna dict {escenario: (fig, axes)}.
    """
    out = {}
    for scen, npz_path in scenario_npz_map.items():
        fig, axes = plot_theta_anom_hov_with_vertical_mode_multi_station_from_npz(
            npz_path=npz_path,
            stations=stations,
            mode_npz_path=None,
            modes_dir=modes_dir,
            mode_index=mode_index,
            n_modes_plot=n_modes_plot,
            vmin=vmin,
            vmax=vmax,
            cmap=cmo.cm.diff,  # colormap secuencial para anomalías de temperatura
            title_prefix=str(scen),
            symmetric=symmetric,
            normalize_station_max=normalize_station_max,
            draw_positive_dashed_contours=draw_positive_dashed_contours,
            positive_contour_lw=positive_contour_lw,
            positive_contour_color=positive_contour_color,
            fig_width=fig_width,
            row_height=row_height
        )
        out[str(scen)] = (fig, axes)
    return out

def plot_w_hov_with_vertical_mode_multi_station_from_npz(
    npz_path: str,
    stations: list[str | int],
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap = "RdBu_r",
    title_prefix: str = "",
    symmetric: bool = True,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.35",
    fig_width: float = 12.0,
    row_height: float = 3.2
):
    """
    Igual que la versión de theta multi-estación, pero para bundles de W:
      salida de guardar_hovmoller_npz_escenario_unico (clave 'diff').
    """
    if not stations:
        raise ValueError("stations no puede estar vacío.")

    d = np.load(npz_path, allow_pickle=True)
    t   = np.array(d["t_horas"], dtype=float)
    RC  = np.array(d["RC"], dtype=float)
    Z3  = np.array(d["diff"], dtype=float)    # (nt, nz, nst)
    meta= json.loads(str(d["meta_json"]))
    labs = [str(x) for x in d["labels"]] if "labels" in d.files else [str(s["label"]) for s in meta.get("stations",[])]
    scenario = str(meta.get("scenario", "")).strip()
    varname = str(meta.get("var", "W"))

    if mode_npz_path is None:
        mode_npz_path = _resolve_vertical_mode_npz_for_scenario(scenario, modes_dir)
    md = np.load(mode_npz_path, allow_pickle=True)

    mode_depth = np.array(md["depth"], dtype=float) if "depth" in md.files else np.array(RC, dtype=float)
    if np.nanmedian(mode_depth) > 0 and np.nanmedian(RC) < 0:
        mode_depth = -mode_depth

    def _interior_zero_depths(profile, depth):
        p = np.asarray(profile, dtype=float)
        z = np.asarray(depth, dtype=float)
        n = p.size
        if n < 3:
            return np.array([], dtype=float)
        z0s = []
        for i in range(1, n - 2):
            p0, p1 = p[i], p[i + 1]
            z0, z1 = z[i], z[i + 1]
            if not (np.isfinite(p0) and np.isfinite(p1) and np.isfinite(z0) and np.isfinite(z1)):
                continue
            if p0 == 0.0:
                z0s.append(float(z0))
            if p0 * p1 < 0.0:
                alpha = -p0 / (p1 - p0)
                z0s.append(float(z0 + alpha * (z1 - z0)))
            elif p1 == 0.0:
                z0s.append(float(z1))
        if not z0s:
            return np.array([], dtype=float)
        z0s = np.array(sorted(z0s), dtype=float)
        uniq = [z0s[0]]
        for zz in z0s[1:]:
            if abs(zz - uniq[-1]) > 1e-6:
                uniq.append(zz)
        return np.array(uniq, dtype=float)

    mode_curves = []
    zero_depths_all = []
    if "modes" in md.files:
        modes = np.array(md["modes"], dtype=float)
        c = np.array(md["c"], dtype=float) if "c" in md.files else np.array([])
        nmax = modes.shape[1]
        nshow = max(1, int(n_modes_plot))
        if nshow == 1:
            m_indices = [min(max(int(mode_index) - 1, 0), nmax - 1)]
        else:
            m_indices = list(range(min(nshow, nmax)))
        for m0 in m_indices:
            prof = modes[:, m0]
            lab = f"modo {m0 + 1}"
            if c.size > m0 and np.isfinite(c[m0]):
                lab += f" (c={c[m0]:.3f} m/s)"
            mode_curves.append((prof, lab))
            zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    elif "phi" in md.files:
        prof = np.array(md["phi"], dtype=float)
        c = float(np.array(md["c"]).squeeze()) if "c" in md.files else np.nan
        mode_label = "modo 2layer"
        if np.isfinite(c):
            mode_label += f" (c={c:.3f} m/s)"
        mode_curves.append((prof, mode_label))
        zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    else:
        raise KeyError(f"{mode_npz_path} no contiene 'modes' ni 'phi'.")

    zuniq = []
    if zero_depths_all:
        zarr = np.array(sorted(zero_depths_all), dtype=float)
        zuniq = [zarr[0]]
        for zz in zarr[1:]:
            if abs(zz - zuniq[-1]) > 1e-6:
                zuniq.append(zz)

    idx_and_labels = []
    for st in stations:
        if isinstance(st, int):
            k = int(st)
            if (k < 0) or (k >= len(labs)):
                raise IndexError(f"station index {k} fuera de rango [0,{len(labs)-1}]")
            idx_and_labels.append((k, labs[k]))
        else:
            lab = str(st)
            if lab not in labs:
                raise KeyError(f"Estación '{lab}' no está en {labs}")
            idx_and_labels.append((labs.index(lab), lab))

    nrows = len(idx_and_labels)
    fig, axes = plt.subplots(
        nrows, 2,
        figsize=(fig_width, row_height * nrows),
        gridspec_kw={"width_ratios": [1.2, 4.8], "wspace": 0.08, "hspace": 0.22},
        sharex="col",
        sharey="all",
        squeeze=False
    )

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(mode_curves), 2)))

    for r, (k, label) in enumerate(idx_and_labels):
        ax_mode = axes[r, 0]
        ax_hov  = axes[r, 1]

        for kk, (prof, lab) in enumerate(mode_curves):
            ax_mode.plot(prof, mode_depth, lw=1.5, color=colors[kk], label=lab)
        ax_mode.axvline(0.0, color="0.5", lw=0.8, ls="--")
        for zz in zuniq:
            ax_mode.axhline(zz, color="0.35", lw=0.9, ls=":")
        ax_mode.set_ylabel("Profundidad (m)")
        if r == nrows - 1:
            ax_mode.set_xlabel("Modo")
        ax_mode.grid(True, alpha=0.25)
        if r == 0:
            ax_mode.set_title(f"{scenario} • modos")
            ax_mode.legend(fontsize=8, loc="best")

        Z = Z3[:, :, k].astype(float)
        if normalize_station_max:
            z_scale = float(np.nanmax(np.abs(Z)))
            if (not np.isfinite(z_scale)) or z_scale <= 0.0:
                z_scale = 1.0
            Z = Z / z_scale

        if vmin is None or vmax is None:
            if normalize_station_max:
                vmin_i, vmax_i = -1.0, 1.0
            else:
                m = float(np.nanmax(np.abs(Z)))
                if m == 0:
                    vmin_i, vmax_i = -1e-9, 1e-9
                else:
                    if symmetric:
                        vmin_i, vmax_i = -m, m
                    else:
                        vmin_i, vmax_i = float(np.nanmin(Z)), float(np.nanmax(Z))
            vvmin = vmin if vmin is not None else vmin_i
            vvmax = vmax if vmax is not None else vmax_i
        else:
            vvmin, vvmax = float(vmin), float(vmax)
            if symmetric:
                M = max(abs(vvmin), abs(vvmax))
                vvmin, vvmax = -M, M

        levels_i = np.linspace(vvmin, vvmax, 61)
        pcm = ax_hov.contourf(
            t, RC, Z.T,
            levels=levels_i,
            cmap=cmap,
            vmin=vvmin,
            vmax=vvmax,
            extend="both"
        )
        if draw_positive_dashed_contours:
            Zmask = np.ma.masked_less_equal(Z.T, 0.0)
            ax_hov.contour(
                t, RC, Zmask, levels=[0.0],
                colors=positive_contour_color,
                linewidths=positive_contour_lw,
                linestyles="--"
            )
        for zz in zuniq:
            ax_hov.axhline(zz, color="w", lw=0.9, ls=":", alpha=0.75)

        if r == nrows - 1:
            ax_hov.set_xlabel("Tiempo (h)")
        ax_hov.set_ylabel("Profundidad (m)")
        ax_hov.set_title(f"{title_prefix} est. {label}")

        if normalize_station_max:
            cbar_label = f"{varname} / max|{varname}|"
        else:
            cbar_label = f"{varname} (bahía - sin bahía)"
        fig.colorbar(pcm, ax=ax_hov, label=cbar_label, shrink=0.98, pad=0.02)

    fig.suptitle(f"{scenario} • {meta.get('source','diff')} • Hovmöller por estaciones", y=0.995)
    plt.tight_layout()
    plt.show()
    return fig, axes

def plot_w_hov_with_vertical_mode_from_npz(
    npz_path: str,
    station: str | int = "A",
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap = "RdBu_r",
    title_prefix: str = "",
    symmetric: bool = True,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.35",
    fig_width: float = 11.0,
    row_height: float = 4.0
):
    """
    Wrapper de una estación para W:
    devuelve una figura con panel de modos (izq) + Hovmöller (der).
    """
    fig, axes = plot_w_hov_with_vertical_mode_multi_station_from_npz(
        npz_path=npz_path,
        stations=[station],
        mode_npz_path=mode_npz_path,
        modes_dir=modes_dir,
        mode_index=mode_index,
        n_modes_plot=n_modes_plot,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        title_prefix=title_prefix,
        symmetric=symmetric,
        normalize_station_max=normalize_station_max,
        draw_positive_dashed_contours=draw_positive_dashed_contours,
        positive_contour_lw=positive_contour_lw,
        positive_contour_color=positive_contour_color,
        fig_width=fig_width,
        row_height=row_height
    )
    return fig, (axes[0, 0], axes[0, 1])

def plot_uv_hov_with_vertical_mode_multi_station_from_npz(
    npz_path: str,
    stations: list[str | int],
    component: str = "both",  # "U", "V" o "both"
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap = "RdBu_r",
    title_prefix: str = "",
    symmetric: bool = True,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.35",
    fig_width: float = 15.0,
    row_height: float = 3.2
):
    """
    Lee un bundle UV (salida de guardar_hovmoller_uv_npz_escenario_unico) y grafica:
      - Columna izquierda: modos verticales del escenario.
      - Columna(s) derecha: Hovmöller de U y/o V por estación.

    component:
      - "U": solo U
      - "V": solo V
      - "both": U y V
    """
    if not stations:
        raise ValueError("stations no puede estar vacío.")

    comp = str(component).strip().upper()
    if comp not in {"U", "V", "BOTH"}:
        raise ValueError("component debe ser 'U', 'V' o 'both'.")

    d = np.load(npz_path, allow_pickle=True)
    t = np.array(d["t_horas"], dtype=float)
    RC = np.array(d["RC"], dtype=float)
    if "U_hov" not in d.files or "V_hov" not in d.files:
        raise KeyError(f"{npz_path} no contiene 'U_hov' y 'V_hov'.")
    U3 = np.array(d["U_hov"], dtype=float)  # (nt, nz, nst)
    V3 = np.array(d["V_hov"], dtype=float)  # (nt, nz, nst)
    meta = json.loads(str(d["meta_json"]))
    labs = [str(x) for x in d["labels"]] if "labels" in d.files else [str(s["label"]) for s in meta.get("stations", [])]

    scenario = str(meta.get("scenario", "")).strip()
    source = str(meta.get("source", "diff"))

    if mode_npz_path is None:
        mode_npz_path = _resolve_vertical_mode_npz_for_scenario(scenario, modes_dir)
    md = np.load(mode_npz_path, allow_pickle=True)

    mode_depth = np.array(md["depth"], dtype=float) if "depth" in md.files else np.array(RC, dtype=float)
    if np.nanmedian(mode_depth) > 0 and np.nanmedian(RC) < 0:
        mode_depth = -mode_depth

    def _interior_zero_depths(profile, depth):
        p = np.asarray(profile, dtype=float)
        z = np.asarray(depth, dtype=float)
        if p.size < 3:
            return np.array([], dtype=float)
        out = []
        for i in range(1, p.size - 2):
            p0, p1 = p[i], p[i + 1]
            z0, z1 = z[i], z[i + 1]
            if not (np.isfinite(p0) and np.isfinite(p1) and np.isfinite(z0) and np.isfinite(z1)):
                continue
            if p0 == 0.0:
                out.append(float(z0))
            if p0 * p1 < 0.0:
                a = -p0 / (p1 - p0)
                out.append(float(z0 + a * (z1 - z0)))
            elif p1 == 0.0:
                out.append(float(z1))
        if not out:
            return np.array([], dtype=float)
        z0s = np.array(sorted(out), dtype=float)
        uniq = [z0s[0]]
        for zz in z0s[1:]:
            if abs(zz - uniq[-1]) > 1e-6:
                uniq.append(zz)
        return np.array(uniq, dtype=float)

    mode_curves = []
    zero_depths_all = []
    if "modes" in md.files:
        modes = np.array(md["modes"], dtype=float)
        c = np.array(md["c"], dtype=float) if "c" in md.files else np.array([])
        nmax = modes.shape[1]
        nshow = max(1, int(n_modes_plot))
        if nshow == 1:
            m_indices = [min(max(int(mode_index) - 1, 0), nmax - 1)]
        else:
            m_indices = list(range(min(nshow, nmax)))
        for m0 in m_indices:
            prof = modes[:, m0]
            lab = f"modo {m0 + 1}"
            if c.size > m0 and np.isfinite(c[m0]):
                lab += f" (c={c[m0]:.3f} m/s)"
            mode_curves.append((prof, lab))
            zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    elif "phi" in md.files:
        prof = np.array(md["phi"], dtype=float)
        c = float(np.array(md["c"]).squeeze()) if "c" in md.files else np.nan
        lab = "modo 2layer"
        if np.isfinite(c):
            lab += f" (c={c:.3f} m/s)"
        mode_curves.append((prof, lab))
        zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    else:
        raise KeyError(f"{mode_npz_path} no contiene 'modes' ni 'phi'.")

    zuniq = []
    if zero_depths_all:
        zarr = np.array(sorted(zero_depths_all), dtype=float)
        zuniq = [zarr[0]]
        for zz in zarr[1:]:
            if abs(zz - zuniq[-1]) > 1e-6:
                zuniq.append(zz)

    idx_and_labels = []
    for st in stations:
        if isinstance(st, int):
            k = int(st)
            if (k < 0) or (k >= len(labs)):
                raise IndexError(f"station index {k} fuera de rango [0,{len(labs)-1}]")
            idx_and_labels.append((k, labs[k]))
        else:
            lab = str(st)
            if lab not in labs:
                raise KeyError(f"Estación '{lab}' no está en {labs}")
            idx_and_labels.append((labs.index(lab), lab))

    fields = [("U", U3)] if comp == "U" else [("V", V3)] if comp == "V" else [("U", U3), ("V", V3)]
    n_hov_cols = len(fields)
    nrows = len(idx_and_labels)
    width_ratios = [1.2] + ([4.8] if n_hov_cols == 1 else [4.2, 4.2])

    fig, axes = plt.subplots(
        nrows, 1 + n_hov_cols,
        figsize=(fig_width, row_height * nrows),
        gridspec_kw={"width_ratios": width_ratios, "wspace": 0.08, "hspace": 0.22},
        sharex="col",
        sharey="all",
        squeeze=False
    )

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(mode_curves), 2)))

    def _limits(Z):
        if vmin is None or vmax is None:
            if normalize_station_max:
                return -1.0, 1.0
            m = float(np.nanmax(np.abs(Z)))
            if m == 0:
                return -1e-9, 1e-9
            return (-m, m) if symmetric else (float(np.nanmin(Z)), float(np.nanmax(Z)))
        vvmin, vvmax = float(vmin), float(vmax)
        if symmetric:
            M = max(abs(vvmin), abs(vvmax))
            vvmin, vvmax = -M, M
        return vvmin, vvmax

    for r, (k, label) in enumerate(idx_and_labels):
        ax_mode = axes[r, 0]
        for kk, (prof, labm) in enumerate(mode_curves):
            ax_mode.plot(prof, mode_depth, lw=1.5, color=colors[kk], label=labm)
        ax_mode.axvline(0.0, color="0.5", lw=0.8, ls="--")
        for zz in zuniq:
            ax_mode.axhline(zz, color="0.35", lw=0.9, ls=":")
        ax_mode.set_ylabel("Profundidad (m)")
        if r == nrows - 1:
            ax_mode.set_xlabel("Modo")
        ax_mode.grid(True, alpha=0.25)
        if r == 0:
            ax_mode.set_title(f"{scenario} • modos")
            ax_mode.legend(fontsize=8, loc="best")

        for cidx, (fname, F3) in enumerate(fields, start=1):
            ax_hov = axes[r, cidx]
            Z = F3[:, :, k].astype(float)
            if normalize_station_max:
                z_scale = float(np.nanmax(np.abs(Z)))
                if (not np.isfinite(z_scale)) or z_scale <= 0.0:
                    z_scale = 1.0
                Z = Z / z_scale

            vvmin, vvmax = _limits(Z)
            levels_i = np.linspace(vvmin, vvmax, 61)
            pcm = ax_hov.contourf(
                t, RC, Z.T,
                levels=levels_i,
                cmap=cmap,
                vmin=vvmin,
                vmax=vvmax,
                extend="both"
            )
            if draw_positive_dashed_contours:
                Zmask = np.ma.masked_less_equal(Z.T, 0.0)
                ax_hov.contour(
                    t, RC, Zmask, levels=[0.0],
                    colors=positive_contour_color,
                    linewidths=positive_contour_lw,
                    linestyles="--"
                )
            for zz in zuniq:
                ax_hov.axhline(zz, color="w", lw=0.9, ls=":", alpha=0.75)

            if r == nrows - 1:
                ax_hov.set_xlabel("Tiempo (h)")
            ax_hov.set_ylabel("Profundidad (m)")
            ax_hov.set_title(f"{title_prefix} est. {label} • {fname}")

            if normalize_station_max:
                cbar_label = f"{fname} / max|{fname}|"
            else:
                cbar_label = f"{fname} ({source})"
            fig.colorbar(pcm, ax=ax_hov, label=cbar_label, shrink=0.98, pad=0.02)

    fig.suptitle(f"{scenario} • {source} • Hovmöller UV por estaciones", y=0.995)
    plt.tight_layout()
    plt.show()
    return fig, axes

def plot_uv_hov_with_vertical_mode_from_npz(
    npz_path: str,
    station: str | int = "A",
    component: str = "both",
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap = "RdBu_r",
    title_prefix: str = "",
    symmetric: bool = True,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.35",
    fig_width: float = 13.0,
    row_height: float = 4.0
):
    """
    Wrapper de una estación para bundles UV:
    devuelve figura con panel de modos + Hovmöller de U/V.
    """
    fig, axes = plot_uv_hov_with_vertical_mode_multi_station_from_npz(
        npz_path=npz_path,
        stations=[station],
        component=component,
        mode_npz_path=mode_npz_path,
        modes_dir=modes_dir,
        mode_index=mode_index,
        n_modes_plot=n_modes_plot,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        title_prefix=title_prefix,
        symmetric=symmetric,
        normalize_station_max=normalize_station_max,
        draw_positive_dashed_contours=draw_positive_dashed_contours,
        positive_contour_lw=positive_contour_lw,
        positive_contour_color=positive_contour_color,
        fig_width=fig_width,
        row_height=row_height
    )
    return fig, axes[0]

def plot_uvw_hov_with_vertical_mode_multi_station_from_npz(
    npz_uv_path: str,
    npz_w_path: str,
    stations: list[str | int],
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin_u: float | None = None,
    vmax_u: float | None = None,
    vmin_v: float | None = None,
    vmax_v: float | None = None,
    vmin_w: float | None = None,
    vmax_w: float | None = None,
    cmap = "RdBu_r",
    title_prefix: str = "",
    symmetric: bool = True,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.35",
    fig_width: float = 19.0,
    row_height: float = 3.2
):
    """
    Grafica modos verticales + Hovmöller de U, V y W por estación.

    Entradas:
      - npz_uv_path: salida de guardar_hovmoller_uv_npz_escenario_unico (U_hov, V_hov)
      - npz_w_path : salida de guardar_hovmoller_npz_escenario_unico para W (diff)
    """
    if not stations:
        raise ValueError("stations no puede estar vacío.")

    duv = np.load(npz_uv_path, allow_pickle=True)
    t_uv = np.array(duv["t_horas"], dtype=float)
    RC_uv = np.array(duv["RC"], dtype=float)
    if "U_hov" not in duv.files or "V_hov" not in duv.files:
        raise KeyError(f"{npz_uv_path} no contiene 'U_hov' y 'V_hov'.")
    U3 = np.array(duv["U_hov"], dtype=float)
    V3 = np.array(duv["V_hov"], dtype=float)
    meta_uv = json.loads(str(duv["meta_json"]))
    labs_uv = [str(x) for x in duv["labels"]] if "labels" in duv.files else [str(s["label"]) for s in meta_uv.get("stations", [])]

    dw = np.load(npz_w_path, allow_pickle=True)
    t_w = np.array(dw["t_horas"], dtype=float)
    RC_w = np.array(dw["RC"], dtype=float)
    if "diff" not in dw.files:
        raise KeyError(f"{npz_w_path} no contiene 'diff' (bundle W).")
    W3 = np.array(dw["diff"], dtype=float)
    meta_w = json.loads(str(dw["meta_json"]))
    labs_w = [str(x) for x in dw["labels"]] if "labels" in dw.files else [str(s["label"]) for s in meta_w.get("stations", [])]

    if t_uv.shape != t_w.shape or np.any(np.abs(t_uv - t_w) > 1e-7):
        raise ValueError("Los tiempos de UV y W no coinciden. Requieren misma malla temporal.")
    if RC_uv.shape != RC_w.shape or np.any(np.abs(RC_uv - RC_w) > 1e-7):
        raise ValueError("Las profundidades RC de UV y W no coinciden.")

    t = t_uv
    RC = RC_uv

    scenario_uv = str(meta_uv.get("scenario", "")).strip()
    scenario_w = str(meta_w.get("scenario", "")).strip()
    if scenario_uv and scenario_w and (scenario_uv != scenario_w):
        raise ValueError(f"Escenarios distintos entre UV ({scenario_uv}) y W ({scenario_w}).")
    scenario = scenario_uv if scenario_uv else scenario_w
    source = str(meta_uv.get("source", meta_w.get("source", "diff")))

    if mode_npz_path is None:
        mode_npz_path = _resolve_vertical_mode_npz_for_scenario(scenario, modes_dir)
    md = np.load(mode_npz_path, allow_pickle=True)

    mode_depth = np.array(md["depth"], dtype=float) if "depth" in md.files else np.array(RC, dtype=float)
    if np.nanmedian(mode_depth) > 0 and np.nanmedian(RC) < 0:
        mode_depth = -mode_depth

    def _interior_zero_depths(profile, depth):
        p = np.asarray(profile, dtype=float)
        z = np.asarray(depth, dtype=float)
        if p.size < 3:
            return np.array([], dtype=float)
        out = []
        for i in range(1, p.size - 2):
            p0, p1 = p[i], p[i + 1]
            z0, z1 = z[i], z[i + 1]
            if not (np.isfinite(p0) and np.isfinite(p1) and np.isfinite(z0) and np.isfinite(z1)):
                continue
            if p0 == 0.0:
                out.append(float(z0))
            if p0 * p1 < 0.0:
                a = -p0 / (p1 - p0)
                out.append(float(z0 + a * (z1 - z0)))
            elif p1 == 0.0:
                out.append(float(z1))
        if not out:
            return np.array([], dtype=float)
        z0s = np.array(sorted(out), dtype=float)
        uniq = [z0s[0]]
        for zz in z0s[1:]:
            if abs(zz - uniq[-1]) > 1e-6:
                uniq.append(zz)
        return np.array(uniq, dtype=float)

    mode_curves = []
    zero_depths_all = []
    if "modes" in md.files:
        modes = np.array(md["modes"], dtype=float)
        c = np.array(md["c"], dtype=float) if "c" in md.files else np.array([])
        nmax = modes.shape[1]
        nshow = max(1, int(n_modes_plot))
        if nshow == 1:
            m_indices = [min(max(int(mode_index) - 1, 0), nmax - 1)]
        else:
            m_indices = list(range(min(nshow, nmax)))
        for m0 in m_indices:
            prof = modes[:, m0]
            lab = f"modo {m0 + 1}"
            if c.size > m0 and np.isfinite(c[m0]):
                lab += f" (c={c[m0]:.3f} m/s)"
            mode_curves.append((prof, lab))
            zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    elif "phi" in md.files:
        prof = np.array(md["phi"], dtype=float)
        c = float(np.array(md["c"]).squeeze()) if "c" in md.files else np.nan
        lab = "modo 2layer"
        if np.isfinite(c):
            lab += f" (c={c:.3f} m/s)"
        mode_curves.append((prof, lab))
        zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    else:
        raise KeyError(f"{mode_npz_path} no contiene 'modes' ni 'phi'.")

    zuniq = []
    if zero_depths_all:
        zarr = np.array(sorted(zero_depths_all), dtype=float)
        zuniq = [zarr[0]]
        for zz in zarr[1:]:
            if abs(zz - zuniq[-1]) > 1e-6:
                zuniq.append(zz)

    idx_and_labels = []
    for st in stations:
        if isinstance(st, int):
            k = int(st)
            if (k < 0) or (k >= len(labs_uv)):
                raise IndexError(f"station index {k} fuera de rango [0,{len(labs_uv)-1}]")
            lab = labs_uv[k]
        else:
            lab = str(st)
            if lab not in labs_uv:
                raise KeyError(f"Estación '{lab}' no está en UV: {labs_uv}")
            k = labs_uv.index(lab)

        if lab not in labs_w:
            raise KeyError(f"Estación '{lab}' no está en bundle W: {labs_w}")
        kw = labs_w.index(lab)
        idx_and_labels.append((k, kw, lab))

    fig, axes = plt.subplots(
        len(idx_and_labels), 4,
        figsize=(fig_width, row_height * len(idx_and_labels)),
        gridspec_kw={"width_ratios": [1.2, 4.1, 4.1, 4.1], "wspace": 0.08, "hspace": 0.22},
        sharex="col",
        sharey="all",
        squeeze=False
    )

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(mode_curves), 2)))

    def _limits(Z, vmin_i, vmax_i):
        if vmin_i is None or vmax_i is None:
            if normalize_station_max:
                return -1.0, 1.0
            m = float(np.nanmax(np.abs(Z)))
            if m == 0:
                return -1e-9, 1e-9
            return (-m, m) if symmetric else (float(np.nanmin(Z)), float(np.nanmax(Z)))
        lo, hi = float(vmin_i), float(vmax_i)
        if symmetric:
            M = max(abs(lo), abs(hi))
            lo, hi = -M, M
        return lo, hi

    fields = [
        ("U", U3, vmin_u, vmax_u),
        ("V", V3, vmin_v, vmax_v),
        ("W", W3, vmin_w, vmax_w),
    ]

    for r, (kuv, kw, label) in enumerate(idx_and_labels):
        ax_mode = axes[r, 0]
        for kk, (prof, labm) in enumerate(mode_curves):
            ax_mode.plot(prof, mode_depth, lw=1.5, color=colors[kk], label=labm)
        ax_mode.axvline(0.0, color="0.5", lw=0.8, ls="--")
        for zz in zuniq:
            ax_mode.axhline(zz, color="0.35", lw=0.9, ls=":")
        ax_mode.set_ylabel("Profundidad (m)")
        if r == len(idx_and_labels) - 1:
            ax_mode.set_xlabel("Modo")
        ax_mode.grid(True, alpha=0.25)
        if r == 0:
            ax_mode.set_title(f"{scenario} • modos")
            ax_mode.legend(fontsize=8, loc="best")

        for cidx, (fname, F3, lo_in, hi_in) in enumerate(fields, start=1):
            ax_hov = axes[r, cidx]
            if fname in ("U", "V"):
                Z = F3[:, :, kuv].astype(float)
            else:
                Z = F3[:, :, kw].astype(float)

            if normalize_station_max:
                z_scale = float(np.nanmax(np.abs(Z)))
                if (not np.isfinite(z_scale)) or z_scale <= 0.0:
                    z_scale = 1.0
                Z = Z / z_scale

            vvmin, vvmax = _limits(Z, lo_in, hi_in)
            levels_i = np.linspace(vvmin, vvmax, 61)
            pcm = ax_hov.contourf(
                t, RC, Z.T,
                levels=levels_i,
                cmap=cmap,
                vmin=vvmin,
                vmax=vvmax,
                extend="both"
            )
            if draw_positive_dashed_contours:
                Zmask = np.ma.masked_less_equal(Z.T, 0.0)
                ax_hov.contour(
                    t, RC, Zmask, levels=[0.0],
                    colors=positive_contour_color,
                    linewidths=positive_contour_lw,
                    linestyles="--"
                )
            for zz in zuniq:
                ax_hov.axhline(zz, color="w", lw=0.9, ls=":", alpha=0.75)

            if r == len(idx_and_labels) - 1:
                ax_hov.set_xlabel("Tiempo (h)")
            ax_hov.set_ylabel("")
            ax_hov.tick_params(labelleft=False)
            ax_hov.set_title(f"{title_prefix} est. {label} • {fname}")

            if normalize_station_max:
                cbar_label = f"{fname} / max|{fname}|"
            else:
                cbar_label = f"{fname} ({source})"
            fig.colorbar(pcm, ax=ax_hov, label=cbar_label, shrink=0.98, pad=0.02)

    fig.suptitle(f"{scenario} • {source} • Hovmöller UVW por estaciones", y=0.995)
    plt.tight_layout()
    plt.show()
    return fig, axes

def plot_uvw_hov_with_vertical_mode_from_npz(
    npz_uv_path: str,
    npz_w_path: str,
    station: str | int = "A",
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin_u: float | None = None,
    vmax_u: float | None = None,
    vmin_v: float | None = None,
    vmax_v: float | None = None,
    vmin_w: float | None = None,
    vmax_w: float | None = None,
    cmap = "RdBu_r",
    title_prefix: str = "",
    symmetric: bool = True,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.35",
    fig_width: float = 17.0,
    row_height: float = 4.0
):
    """
    Wrapper de una estación para UVW: modos + Hovmöller U,V,W.
    """
    fig, axes = plot_uvw_hov_with_vertical_mode_multi_station_from_npz(
        npz_uv_path=npz_uv_path,
        npz_w_path=npz_w_path,
        stations=[station],
        mode_npz_path=mode_npz_path,
        modes_dir=modes_dir,
        mode_index=mode_index,
        n_modes_plot=n_modes_plot,
        vmin_u=vmin_u,
        vmax_u=vmax_u,
        vmin_v=vmin_v,
        vmax_v=vmax_v,
        vmin_w=vmin_w,
        vmax_w=vmax_w,
        cmap=cmap,
        title_prefix=title_prefix,
        symmetric=symmetric,
        normalize_station_max=normalize_station_max,
        draw_positive_dashed_contours=draw_positive_dashed_contours,
        positive_contour_lw=positive_contour_lw,
        positive_contour_color=positive_contour_color,
        fig_width=fig_width,
        row_height=row_height
    )
    return fig, axes[0]

def plot_uvw_theta_anom_with_vertical_mode_multi_station_from_npz(
    npz_uv_path: str,
    npz_w_path: str,
    npz_theta_path: str,
    stations: list[str | int],
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin_u: float | None = None,
    vmax_u: float | None = None,
    vmin_v: float | None = None,
    vmax_v: float | None = None,
    vmin_w: float | None = None,
    vmax_w: float | None = None,
    vmin_theta: float | None = None,
    vmax_theta: float | None = None,
    cmap_uvws = "RdBu_r",
    cmap_theta = cmo.cm.diff,
    title_prefix: str = "",
    symmetric_uvws: bool = True,
    symmetric_theta: bool = False,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.35",
    fig_width: float = 24.0,
    row_height: float = 3.2,
    save_fig: bool = False,
    save_dir: str | Path = "./figs",
    dpi: int = 170,
    save_fmt: str = "pdf",
):
    """
    Grafica modos verticales + Hovmöller de U, V, W y anomalía de temperatura por estación.

    Entradas:
      - npz_uv_path: salida de guardar_hovmoller_uv_npz_escenario_unico (U_hov, V_hov)
      - npz_w_path : salida de guardar_hovmoller_npz_escenario_unico para W (diff)
      - npz_theta_path: salida de guardar_theta_anoma_npz_escenario_unico (anom)

    Nota:
      - U/V/W usan `cmap_uvws`.
      - Theta anómala usa `cmap_theta` (default: cmo.cm.diff).
    """
    if not stations:
        raise ValueError("stations no puede estar vacío.")

    duv = np.load(npz_uv_path, allow_pickle=True)
    t_uv = np.array(duv["t_horas"], dtype=float)
    RC_uv = np.array(duv["RC"], dtype=float)
    if "U_hov" not in duv.files or "V_hov" not in duv.files:
        raise KeyError(f"{npz_uv_path} no contiene 'U_hov' y 'V_hov'.")
    U3 = np.array(duv["U_hov"], dtype=float)
    V3 = np.array(duv["V_hov"], dtype=float)
    meta_uv = json.loads(str(duv["meta_json"]))
    labs_uv = [str(x) for x in duv["labels"]] if "labels" in duv.files else [str(s["label"]) for s in meta_uv.get("stations", [])]

    dw = np.load(npz_w_path, allow_pickle=True)
    t_w = np.array(dw["t_horas"], dtype=float)
    RC_w = np.array(dw["RC"], dtype=float)
    if "diff" not in dw.files:
        raise KeyError(f"{npz_w_path} no contiene 'diff' (bundle W).")
    W3 = np.array(dw["diff"], dtype=float)
    meta_w = json.loads(str(dw["meta_json"]))
    labs_w = [str(x) for x in dw["labels"]] if "labels" in dw.files else [str(s["label"]) for s in meta_w.get("stations", [])]

    dt = np.load(npz_theta_path, allow_pickle=True)
    t_th = np.array(dt["t_horas"], dtype=float)
    RC_th = np.array(dt["RC"], dtype=float)
    if "anom" not in dt.files:
        raise KeyError(f"{npz_theta_path} no contiene 'anom' (bundle theta anómala).")
    TH3 = np.array(dt["anom"], dtype=float)
    meta_th = json.loads(str(dt["meta_json"]))
    labs_th = [str(x) for x in dt["labels"]] if "labels" in dt.files else [str(s["label"]) for s in meta_th.get("stations", [])]

    if t_uv.shape != t_w.shape or np.any(np.abs(t_uv - t_w) > 1e-7):
        raise ValueError("Los tiempos de UV y W no coinciden.")
    if t_uv.shape != t_th.shape or np.any(np.abs(t_uv - t_th) > 1e-7):
        raise ValueError("Los tiempos de UV y Theta no coinciden.")
    if RC_uv.shape != RC_w.shape or np.any(np.abs(RC_uv - RC_w) > 1e-7):
        raise ValueError("Las profundidades RC de UV y W no coinciden.")
    if RC_uv.shape != RC_th.shape or np.any(np.abs(RC_uv - RC_th) > 1e-7):
        raise ValueError("Las profundidades RC de UV y Theta no coinciden.")

    t = t_uv
    RC = RC_uv

    scen_uv = str(meta_uv.get("scenario", "")).strip()
    scen_w = str(meta_w.get("scenario", "")).strip()
    scen_th = str(meta_th.get("scenario", "")).strip()
    scen_nonempty = [s for s in (scen_uv, scen_w, scen_th) if s]
    if len(set(scen_nonempty)) > 1:
        raise ValueError(f"Escenarios distintos entre archivos: UV={scen_uv}, W={scen_w}, Theta={scen_th}.")
    scenario = scen_nonempty[0] if scen_nonempty else ""
    source = str(meta_uv.get("source", meta_w.get("source", meta_th.get("source", "diff"))))

    if mode_npz_path is None:
        mode_npz_path = _resolve_vertical_mode_npz_for_scenario(scenario, modes_dir)
    md = np.load(mode_npz_path, allow_pickle=True)

    mode_depth = np.array(md["depth"], dtype=float) if "depth" in md.files else np.array(RC, dtype=float)
    if np.nanmedian(mode_depth) > 0 and np.nanmedian(RC) < 0:
        mode_depth = -mode_depth

    def _interior_zero_depths(profile, depth):
        p = np.asarray(profile, dtype=float)
        z = np.asarray(depth, dtype=float)
        if p.size < 3:
            return np.array([], dtype=float)
        out = []
        for i in range(1, p.size - 2):
            p0, p1 = p[i], p[i + 1]
            z0, z1 = z[i], z[i + 1]
            if not (np.isfinite(p0) and np.isfinite(p1) and np.isfinite(z0) and np.isfinite(z1)):
                continue
            if p0 == 0.0:
                out.append(float(z0))
            if p0 * p1 < 0.0:
                a = -p0 / (p1 - p0)
                out.append(float(z0 + a * (z1 - z0)))
            elif p1 == 0.0:
                out.append(float(z1))
        if not out:
            return np.array([], dtype=float)
        z0s = np.array(sorted(out), dtype=float)
        uniq = [z0s[0]]
        for zz in z0s[1:]:
            if abs(zz - uniq[-1]) > 1e-6:
                uniq.append(zz)
        return np.array(uniq, dtype=float)

    mode_curves = []
    zero_depths_all = []
    if "modes" in md.files:
        modes = np.array(md["modes"], dtype=float)
        c = np.array(md["c"], dtype=float) if "c" in md.files else np.array([])
        nmax = modes.shape[1]
        nshow = max(1, int(n_modes_plot))
        if nshow == 1:
            m_indices = [min(max(int(mode_index) - 1, 0), nmax - 1)]
        else:
            m_indices = list(range(min(nshow, nmax)))
        for m0 in m_indices:
            prof = modes[:, m0]
            lab = f"modo {m0 + 1}"
            if c.size > m0 and np.isfinite(c[m0]):
                lab += f" (c={c[m0]:.3f} m/s)"
            mode_curves.append((prof, lab))
            zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    elif "phi" in md.files:
        prof = np.array(md["phi"], dtype=float)
        c = float(np.array(md["c"]).squeeze()) if "c" in md.files else np.nan
        lab = "modo 2layer"
        if np.isfinite(c):
            lab += f" (c={c:.3f} m/s)"
        mode_curves.append((prof, lab))
        zero_depths_all.extend(_interior_zero_depths(prof, mode_depth).tolist())
    else:
        raise KeyError(f"{mode_npz_path} no contiene 'modes' ni 'phi'.")

    zuniq = []
    if zero_depths_all:
        zarr = np.array(sorted(zero_depths_all), dtype=float)
        zuniq = [zarr[0]]
        for zz in zarr[1:]:
            if abs(zz - zuniq[-1]) > 1e-6:
                zuniq.append(zz)

    idx_and_labels = []
    for st in stations:
        if isinstance(st, int):
            k_uv = int(st)
            if (k_uv < 0) or (k_uv >= len(labs_uv)):
                raise IndexError(f"station index {k_uv} fuera de rango [0,{len(labs_uv)-1}]")
            lab = labs_uv[k_uv]
        else:
            lab = str(st)
            if lab not in labs_uv:
                raise KeyError(f"Estación '{lab}' no está en UV: {labs_uv}")
            k_uv = labs_uv.index(lab)

        if lab not in labs_w:
            raise KeyError(f"Estación '{lab}' no está en W: {labs_w}")
        if lab not in labs_th:
            raise KeyError(f"Estación '{lab}' no está en Theta: {labs_th}")

        k_w = labs_w.index(lab)
        k_th = labs_th.index(lab)
        idx_and_labels.append((k_uv, k_w, k_th, lab))

    fig, axes = plt.subplots(
        len(idx_and_labels), 5,
        figsize=(fig_width, row_height * len(idx_and_labels)),
        gridspec_kw={"width_ratios": [1.2, 4.0, 4.0, 4.0, 4.0], "wspace": 0.08, "hspace": 0.22},
        sharex="col",
        sharey="all",
        squeeze=False
    )

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(mode_curves), 2)))

    def _limits(Z, vmin_i, vmax_i, symmetric_local):
        if vmin_i is None or vmax_i is None:
            if normalize_station_max:
                return -1.0, 1.0
            m = float(np.nanmax(np.abs(Z)))
            if m == 0:
                return -1e-9, 1e-9
            return (-m, m) if symmetric_local else (float(np.nanmin(Z)), float(np.nanmax(Z)))
        lo, hi = float(vmin_i), float(vmax_i)
        if symmetric_local:
            M = max(abs(lo), abs(hi))
            lo, hi = -M, M
        return lo, hi

    fields = [
        ("U", U3, vmin_u, vmax_u, cmap_uvws, symmetric_uvws),
        ("V", V3, vmin_v, vmax_v, cmap_uvws, symmetric_uvws),
        ("W", W3, vmin_w, vmax_w, cmap_uvws, symmetric_uvws),
        ("Theta", TH3, vmin_theta, vmax_theta, cmap_theta, symmetric_theta),
    ]

    for r, (k_uv, k_w, k_th, label) in enumerate(idx_and_labels):
        ax_mode = axes[r, 0]
        for kk, (prof, labm) in enumerate(mode_curves):
            ax_mode.plot(prof, mode_depth, lw=1.5, color=colors[kk], label=labm)
        ax_mode.axvline(0.0, color="0.5", lw=0.8, ls="--")
        for zz in zuniq:
            ax_mode.axhline(zz, color="0.35", lw=0.9, ls=":")
        ax_mode.set_ylabel("Profundidad (m)")
        if r == len(idx_and_labels) - 1:
            ax_mode.set_xlabel("Modo")
        ax_mode.grid(True, alpha=0.25)
        if r == 0:
            ax_mode.set_title(f"{scenario} • modos")
            ax_mode.legend(fontsize=8, loc="best")

        for cidx, (fname, F3, lo_in, hi_in, cmap_local, sym_local) in enumerate(fields, start=1):
            ax_hov = axes[r, cidx]
            if fname in ("U", "V"):
                Z = F3[:, :, k_uv].astype(float)
            elif fname == "W":
                Z = F3[:, :, k_w].astype(float)
            else:
                Z = F3[:, :, k_th].astype(float)

            if normalize_station_max:
                z_scale = float(np.nanmax(np.abs(Z)))
                if (not np.isfinite(z_scale)) or z_scale <= 0.0:
                    z_scale = 1.0
                Z = Z / z_scale
                print(f"  [{scenario}] est. {label} | {fname:5s}: max|·| = {z_scale:.4e}")

            vvmin, vvmax = _limits(Z, lo_in, hi_in, sym_local)
            levels_i = np.linspace(vvmin, vvmax, 61)
            pcm = ax_hov.contourf(
                t, RC, Z.T,
                levels=levels_i,
                cmap=cmap_local,
                vmin=vvmin,
                vmax=vvmax,
                extend="both"
            )
            if draw_positive_dashed_contours:
                Zmask = np.ma.masked_less_equal(Z.T, 0.0)
                ax_hov.contour(
                    t, RC, Zmask, levels=[0.0],
                    colors=positive_contour_color,
                    linewidths=positive_contour_lw,
                    linestyles="--"
                )
            for zz in zuniq:
                ax_hov.axhline(zz, color="w", lw=0.9, ls=":", alpha=0.75)

            if r == len(idx_and_labels) - 1:
                ax_hov.set_xlabel("Tiempo (h)")
            ax_hov.set_ylabel("")
            ax_hov.tick_params(labelleft=False)
            ax_hov.set_title(f"{title_prefix} est. {label} • {fname}")

            # Colorbar solo en las dos últimas columnas de datos (W y Theta)
            if cidx >= len(fields) - 1:
                if normalize_station_max:
                    cbar_label = f"{fname} / max|{fname}|"
                else:
                    cbar_label = f"{fname} ({source})"
                cb = fig.colorbar(pcm, ax=ax_hov, label=cbar_label, shrink=0.98, pad=0.02)
                cb.ax.yaxis.set_major_formatter(
                    plt.matplotlib.ticker.FormatStrFormatter("%.1f")
                )

    fig.suptitle(f"{scenario} • {source} • Hovmöller UVW + Theta por estaciones", y=0.995)
    plt.tight_layout()
    if save_fig:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        stations_str = "_".join(str(s) for s in stations)
        out_png = save_dir / f"{scenario}_{source}_hovUVWTheta_stations_{stations_str}.{save_fmt}"
        fig.savefig(out_png, dpi=dpi, bbox_inches="tight", format=save_fmt)
        print(f"[OK] Figura guardada: {out_png}")
    plt.show()
    return fig, axes

def plot_uvw_theta_anom_with_vertical_mode_from_npz(
    npz_uv_path: str,
    npz_w_path: str,
    npz_theta_path: str,
    station: str | int = "A",
    mode_npz_path: str | None = None,
    modes_dir: str | Path = "./modes_npz",
    mode_index: int = 1,
    n_modes_plot: int = 3,
    vmin_u: float | None = None,
    vmax_u: float | None = None,
    vmin_v: float | None = None,
    vmax_v: float | None = None,
    vmin_w: float | None = None,
    vmax_w: float | None = None,
    vmin_theta: float | None = None,
    vmax_theta: float | None = None,
    cmap_uvws = "RdBu_r",
    cmap_theta = cmo.cm.diff,
    title_prefix: str = "",
    symmetric_uvws: bool = True,
    symmetric_theta: bool = False,
    normalize_station_max: bool = True,
    draw_positive_dashed_contours: bool = True,
    positive_contour_lw: float = 1.8,
    positive_contour_color: str = "0.35",
    fig_width: float = 22.0,
    row_height: float = 4.0
):
    """
    Wrapper de una estación para modos + U,V,W + Theta anómala.
    """
    fig, axes = plot_uvw_theta_anom_with_vertical_mode_multi_station_from_npz(
        npz_uv_path=npz_uv_path,
        npz_w_path=npz_w_path,
        npz_theta_path=npz_theta_path,
        stations=[station],
        mode_npz_path=mode_npz_path,
        modes_dir=modes_dir,
        mode_index=mode_index,
        n_modes_plot=n_modes_plot,
        vmin_u=vmin_u,
        vmax_u=vmax_u,
        vmin_v=vmin_v,
        vmax_v=vmax_v,
        vmin_w=vmin_w,
        vmax_w=vmax_w,
        vmin_theta=vmin_theta,
        vmax_theta=vmax_theta,
        cmap_uvws=cmap_uvws,
        cmap_theta=cmap_theta,
        title_prefix=title_prefix,
        symmetric_uvws=symmetric_uvws,
        symmetric_theta=symmetric_theta,
        normalize_station_max=normalize_station_max,
        draw_positive_dashed_contours=draw_positive_dashed_contours,
        positive_contour_lw=positive_contour_lw,
        positive_contour_color=positive_contour_color,
        fig_width=fig_width,
        row_height=row_height
    )
    return fig, axes[0]

def plot_anom_timeseries_and_spectrum_across_scenarios(
    scenario_npz_map: dict[str, str],   # {"linear": "path_linear.npz", "realistic":"...", ...}
    station: str | int = "A",
    depth_m: float = 110.0,
    detrend_series: bool = True,
    apply_hann: bool = True,
    freq_cph_max: float = 0.8,
    seiche_freq_cph: float | None = 0.303356,  # pon None si no quieres la línea
    share_ylim: bool = True
):
    """
    Para cada escenario (NPZ), extrae la serie de anomalía a 'depth_m' y:
      (a) la grafica vs tiempo,
      (b) grafica su espectro de frecuencias (PSD simple) en cph.
    Si los tiempos difieren entre escenarios, interpola a una malla común.
    """
    # 1) Cargar series
    series = {}
    dts = []
    tmins = []
    tmaxs = []
    for scen, path in scenario_npz_map.items():
        t, y, zref, lab = extract_anom_timeseries_at_depth(path, station, depth_m)
        if detrend_series:
            y = detrend(y, type='linear')
        dts.append(np.median(np.diff(t)) if t.size > 1 else np.nan)
        tmins.append(t[0]); tmaxs.append(t[-1])
        series[scen] = {"t": t, "y": y, "zref": zref, "label": lab}

    # 2) Construir malla temporal común (intersección), paso = max(median dt) para no sobre-muestrear
    t0 = max(tmins)
    t1 = min(tmaxs)
    if t1 <= t0:
        raise ValueError("Las coberturas temporales de escenarios no se traslapan.")
    dt = np.nanmax(dts)
    if not np.isfinite(dt) or dt <= 0:
        # fallback: usar el primer dt válido
        dt = next((v for v in dts if np.isfinite(v) and v > 0), 1.0)
    t_common = np.arange(t0, t1 + 1e-9, dt)

    # 3) Interpolar series a la malla común
    for scen, obj in series.items():
        obj["y_interp"] = np.interp(t_common, obj["t"], obj["y"])

    # 4) Figura: dos filas (serie y espectro)
    fig, (ax_ts, ax_sp) = plt.subplots(2, 1, figsize=(10, 8), sharex=False)

    # (a) Series temporales
    ymin, ymax = np.inf, -np.inf
    for scen, obj in series.items():
        ax_ts.plot(t_common, obj["y_interp"], label=f"{scen} (z≈{abs(obj['zref']):.1f} m)")
        ymin = min(ymin, np.nanmin(obj["y_interp"]))
        ymax = max(ymax, np.nanmax(obj["y_interp"]))
    ax_ts.set_xlabel("Tiempo (h)")
    ax_ts.set_ylabel("Anomalía Temp (°C)")
    ax_ts.set_title(f"Evolución anomalías • estación {series[list(series.keys())[0]]['label']} • profundidad ≈ {depth_m:.1f} m")
    ax_ts.grid(True, ls=":")
    ax_ts.legend()
    if share_ylim and np.isfinite(ymin) and np.isfinite(ymax):
        m = max(abs(ymin), abs(ymax))
        ax_ts.set_ylim(-m, m)

    # (b) Espectro
    for scen, obj in series.items():
        y = obj["y_interp"]
        N = y.size
        if N < 2:
            continue
        # ventana
        if apply_hann:
            w = get_window("hann", N, fftbins=True)
            yw = y * w
            U = (w**2).sum()  # corrección de energía de la ventana
        else:
            yw = y
            U = N
        dt_h = t_common[1] - t_common[0]          # horas
        dt_s = dt_h * 3600.0
        freqs = np.fft.rfftfreq(N, d=dt_s)        # Hz
        Y = np.fft.rfft(yw)
        PSD = (np.abs(Y)**2) / (U * (1.0/dt_s)) * 3600.0   # °C² cph⁻¹

        # convertir a ciclos por hora
        f_cph = freqs * 3600.0
        ax_sp.plot(f_cph, PSD, label=scen)

    ax_sp.set_xlim(0.0, freq_cph_max)
    ax_sp.set_yscale("log")
    ax_sp.set_xlabel("Frecuencia (ciclos/hora)")
    ax_sp.set_ylabel("PSD (°C² cph⁻¹)")
    ax_sp.set_title("Espectro de anomalías (serie a profundidad fija)")
    ax_sp.grid(True, which="both", ls=":")
    if seiche_freq_cph is not None:
        ax_sp.axvline(seiche_freq_cph, color="k", ls="--", lw=0.8, label="Seiche")
    ax_sp.legend()

    plt.tight_layout()
    plt.show()
    return fig, (ax_ts, ax_sp)

def plot_anom_timeseries_and_spectrum_all_stations(
    npz_path: str,
    depth_m: float = 110.0,
    detrend_series: bool = True,
    apply_hann: bool = True,
    symmetric_ylim: bool = True,
    freq_cph_max: float = 0.8,
    seiche_freq_cph: float | None = 0.303356,
    vspan: tuple[float,float] | None = None,   # (t0_h, t1_h) para sombrear una ventana
    savepath: str | None = None,               # ruta para guardar figura (png/pdf)
):
    """
    Compara TODAS las estaciones de un escenario (un NPZ) a una profundidad fija:
      (a) Serie temporal de anomalía por estación.
      (b) PSD (espectro) por estación (superpuestos).

    Parámetros clave:
      - npz_path: ruta al .npz (salida de guardar_theta_anoma_npz_escenario_unico).
      - depth_m : profundidad positiva (p.ej. 110.0).
      - detrend_series: aplica detrend lineal a cada serie.
      - apply_hann: ventana de Hann antes de la FFT.
      - symmetric_ylim: fuerza eje Y simétrico en la serie temporal.
      - freq_cph_max: límite superior del eje de frecuencia en ciclos/hora.
      - seiche_freq_cph: dibuja línea vertical; None para omitir.
      - vspan: sombrear una ventana en tiempo (t0, t1) en horas.
      - savepath: si se pasa, guarda la figura.
    """
    d = np.load(npz_path, allow_pickle=True)
    t_h   = d["t_horas"]       # (nt,)
    RC    = d["RC"]            # (nz,) negativas
    A3    = d["anom"]          # (nt, nz, nst)
    labels= [str(x) for x in d["labels"]]
    # colors es U16 o vacío; convertir blanks a None
    raw_cols = [str(x) for x in d["colors"]] if "colors" in d.files else []
    colors = [(c if c and c.lower() != "none" else None) for c in raw_cols] \
             if raw_cols else [None]*len(labels)

    meta = json.loads(str(d["meta_json"]))
    scenario = meta.get("scenario","")
    varname  = meta.get("var","Temp")

    # profundidad más cercana
    target = -abs(depth_m)
    kz = int(np.argmin(np.abs(RC - target)))
    zref = float(RC[kz])

    # series por estación
    Y = []   # lista de arrays (nt,)
    for k in range(len(labels)):
        y = A3[:, kz, k].astype(float)
        if detrend_series and y.size > 1:
            y = detrend(y, type='linear')
        Y.append(y)

    # figura
    fig, (ax_ts, ax_sp) = plt.subplots(2, 1, figsize=(11, 8), sharex=False)

    # (a) SERIE TEMPORAL
    ymin, ymax = np.inf, -np.inf
    for k, lab in enumerate(labels):
        col = colors[k] if colors[k] else None
        ax_ts.plot(t_h, Y[k], label=lab, color=col)
        ymin = min(ymin, np.nanmin(Y[k]))
        ymax = max(ymax, np.nanmax(Y[k]))

    ax_ts.set_xlabel("Tiempo (h)")
    ax_ts.set_ylabel(f"Anomalía {varname} (°C)")
    ax_ts.set_title(f"{scenario} • todas las estaciones • z≈{abs(zref):.1f} m")
    ax_ts.grid(True, ls=":", alpha=0.6)
    ax_ts.legend(ncols=min(5, len(labels)), fontsize=9)

    if vspan is not None:
        t0, t1 = vspan
        ax_ts.axvspan(t0, t1, color='k', alpha=0.08)

    if np.isfinite(ymin) and np.isfinite(ymax):
        if symmetric_ylim:
            m = max(abs(ymin), abs(ymax))
            ax_ts.set_ylim(-m, m)
        else:
            ax_ts.set_ylim(ymin, ymax)

    # (b) ESPECTRO
    # paso temporal en segundos (asumimos uniforme en este NPZ)
    if t_h.size > 1:
        dt_h = float(np.median(np.diff(t_h)))
    else:
        dt_h = 1.0
    dt_s = dt_h * 3600.0

    for k, lab in enumerate(labels):
        y = Y[k]
        N = y.size
        if N < 4:
            continue
        if apply_hann:
            w = get_window("hann", N, fftbins=True)
            yw = y * w
            U = (w**2).sum()          # corrección de energía de ventana
        else:
            yw = y
            U  = N

        freqs = np.fft.rfftfreq(N, d=dt_s)   # Hz
        Yf = np.fft.rfft(yw)
        PSD = (np.abs(Yf)**2) / (U * (1.0/dt_s)) * 3600.0  # var² cph⁻¹

        f_cph = freqs * 3600.0
        col = colors[k] if colors[k] else None
        ax_sp.plot(f_cph, PSD, label=lab, alpha=0.9, lw=1.4, color=col)

    ax_sp.set_xlim(0.0, float(freq_cph_max))
    ax_sp.set_yscale("log")
    ax_sp.set_xlabel("Frecuencia (ciclos/hora)")
    ax_sp.set_ylabel(f"PSD ({varname}² cph⁻¹)")
    ax_sp.set_title("Espectro por estación")
    ax_sp.grid(True, which="both", ls=":", alpha=0.6)
    if seiche_freq_cph is not None:
        ax_sp.axvline(seiche_freq_cph, color="k", ls="--", lw=0.9, label="Seiche")

    # mostrar leyenda compacta si aporta (evita duplicar con la de arriba)
    handles, lbls = ax_sp.get_legend_handles_labels()
    if len(lbls) and len(lbls) != len(labels)+1:  # heurística
        ax_sp.legend(ncols=min(5, len(labels)), fontsize=9)

    plt.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=200, bbox_inches="tight")
    plt.show()
    return fig, (ax_ts, ax_sp)

# -─────────────────────────────────────────────────────────────────────────────
# Funciones para animaciones de perfiles verticales
# -─────────────────────────────────────────────────────────────────────────────    

def _make_stations_from_grid_inline(gr, j_centro: int = _J_CENTRO_DEFAULT, radio: int = _RADIO_DEFAULT, incluir_gamma_omega: bool = True):
    """
    Genera estaciones con labels y colores. Devuelve stations_meta: lista de dicts {label,i,j,color}.
    """
    Depth = np.array(gr.variables['Depth'][:], dtype=float)  # (ny, nx)
    ny, nx = Depth.shape
    i_centro = nx // 2
    clamp = lambda v, lo, hi: max(lo, min(hi, v))
    j_centro = clamp(j_centro, 0, ny-1)

    iS1, jS1 = i_centro, clamp(j_centro + 30, 0, ny-1)
    iS2, jS2 = i_centro, clamp(j_centro + 10, 0, ny-1)
    iA, jA = clamp(i_centro - radio, 0, nx-1), j_centro       # A
    iB, jB = clamp(i_centro + radio, 0, nx-1), j_centro       # B
    iC, jC = i_centro, clamp(j_centro - radio, 0, ny-1)       # C
    stations = [
        {"label":"Seiche", "i":int(iS1), "j":int(jS1), "color":"#9d02d7"},
        {"label":"Seiche2", "i":int(iS2), "j":int(jS2), "color":"#ffd700"},
        {"label":"A",      "i":int(iA), "j":int(jA), "color":"#ea5f94"},
        {"label":"B",      "i":int(iB), "j":int(jB), "color":"#0000ff"},
        {"label":"C",      "i":int(iC), "j":int(jC), "color":"#fa8775"},
    ]
    if incluir_gamma_omega:
        iG, jG = i_centro, clamp(j_centro - 2*radio, 0, ny-1) # Gamma
        iO, jO = i_centro, j_centro                           # Omega
        stations += [
            {"label":"Gamma", "i":int(iG), "j":int(jG), "color":"limegreen"},
            {"label":"Omega", "i":int(iO), "j":int(jO), "color":"cyan"},
        ]
    return stations

# ---------------- Guardar sección ΔTemp (bay−nobay) en (y,z) o (x,z) ----------------
def guardar_seccion_temp_anom_npz(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    varname: str = "Temp",                # 3D: (t,z,y,x)
    orient: str = "yz",                   # "yz" (i fijo, barrido en y) | "xz" (j fijo, barrido en x)
    start_station_label: str = "Gamma",   # desde qué estación arrancar
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern   = "grid.t*.nc",
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    incluir_gamma_omega: bool = True,
    RC_cut: float | None = None,          # p.ej. -250: recorta a capas menos profundas que ese valor
    lateral_len: int | None = None,       # limita la extensión lateral desde la estación de arranque
    mask_grid_source: str = "bay",        # "bay" | "nobay" | "intersection" | "union"
    apply_mask_to_anom: bool = False,
) -> str:
    """
    Construye y guarda un NPZ con la sección vertical de ΔTemp=bay−nobay.
    Devuelve la ruta del .npz guardado.
    """
    save_dir.mkdir(parents=True, exist_ok=True)
    orient = orient.lower().strip()
    if orient not in ("yz", "xz"):
        raise ValueError("orient debe ser 'yz' o 'xz'.")

    # GRID y estaciones
    gr = _open_mnc_or_raise(outdir_bay, grid_pattern)
    gr_nobay = None
    try:
        def _get_var_any(nc_obj, candidates, label):
            for name in candidates:
                if name in nc_obj.variables:
                    return nc_obj.variables[name]
            raise KeyError(f"No encontré {label}. Probé: {candidates}")

        if mask_grid_source not in ("bay", "nobay", "intersection", "union"):
            raise ValueError("mask_grid_source debe ser 'bay', 'nobay', 'intersection' o 'union'.")
        if mask_grid_source in ("nobay", "intersection", "union"):
            gr_nobay = _open_mnc_or_raise(outdir_nobay, grid_pattern)

        # Coordenadas horizontales robustas a distintos nombres/formas
        if "XC" in gr.variables:
            XC_raw = np.asarray(gr.variables["XC"][:], dtype=float)
        elif "xC" in gr.variables:
            XC_raw = np.asarray(gr.variables["xC"][:], dtype=float)
        else:
            raise KeyError("No encontré coordenada X ('XC' o 'xC') en grid.")

        if "YC" in gr.variables:
            YC_raw = np.asarray(gr.variables["YC"][:], dtype=float)
        elif "yC" in gr.variables:
            YC_raw = np.asarray(gr.variables["yC"][:], dtype=float)
        else:
            raise KeyError("No encontré coordenada Y ('YC' o 'yC') en grid.")

        XC = XC_raw[0, :] if XC_raw.ndim == 2 else XC_raw.squeeze()
        YC = YC_raw[:, 0] if YC_raw.ndim == 2 else YC_raw.squeeze()
        if XC.ndim != 1 or YC.ndim != 1:
            raise ValueError("No pude convertir XC/YC a vectores 1D.")

        RC = np.array(_get_var_any(gr, ("RC", "rC", "Z"), "coordenada vertical (RC/rC/Z)")[:], dtype=float).squeeze()
        if RC.ndim > 1:
            RC = RC[..., 0, 0]
        HFacC_bay = np.array(_get_var_any(gr, ("HFacC", "hFacC"), "máscara HFacC/hFacC")[:], dtype=float)  # (nz,ny,nx)
        HFacC_nobay = None
        if gr_nobay is not None:
            HFacC_nobay = np.array(
                _get_var_any(gr_nobay, ("HFacC", "hFacC"), "máscara HFacC/hFacC en nobay")[:], dtype=float
            )
            if HFacC_nobay.shape != HFacC_bay.shape:
                raise ValueError(
                    f"hFacC bay y nobay tienen forma distinta: {HFacC_bay.shape} vs {HFacC_nobay.shape}"
                )

        stations_meta = _make_stations_from_grid_inline(
            gr, j_centro=j_centro, radio=radio, incluir_gamma_omega=incluir_gamma_omega
        )
        lablist = [s["label"] for s in stations_meta]
        if start_station_label not in lablist:
            raise KeyError(f"No hallé la estación '{start_station_label}' en {lablist}")
        st0 = stations_meta[lablist.index(start_station_label)]
        iG, jG = int(st0["i"]), int(st0["j"])

        # recorte lateral según orientación
        if orient == "yz":
            j0 = jG
            j1 = YC.size if lateral_len is None else min(YC.size, jG + max(1, int(lateral_len)))
            coord_full = YC[j0:j1]
            coord_name = "YC"
        else:
            i0 = iG
            i1 = XC.size if lateral_len is None else min(XC.size, iG + max(1, int(lateral_len)))
            coord_full = XC[i0:i1]
            coord_name = "XC"

        # recorte vertical
        if RC_cut is not None:
            kz_arr = np.where(RC > RC_cut)[0]
            if kz_arr.size == 0:
                raise ValueError(f"No hay niveles con RC > {RC_cut} m.")
            k0 = int(kz_arr.min())
            k1 = int(kz_arr.max()) + 1
        else:
            k0, k1 = 0, int(RC.size)
        k_slice = slice(k0, k1)
        RC_sel = RC[k_slice]

        if orient == "yz":
            wet_bay = HFacC_bay[k_slice, j0:j1, iG] > 0.0
            wet_nobay = (HFacC_nobay[k_slice, j0:j1, iG] > 0.0) if HFacC_nobay is not None else None
        else:
            wet_bay = HFacC_bay[k_slice, jG, i0:i1] > 0.0
            wet_nobay = (HFacC_nobay[k_slice, jG, i0:i1] > 0.0) if HFacC_nobay is not None else None

        if mask_grid_source == "bay":
            mask2d = wet_bay
        elif mask_grid_source == "nobay":
            mask2d = wet_nobay
        elif mask_grid_source == "intersection":
            mask2d = wet_bay & wet_nobay
        else:
            mask2d = wet_bay | wet_nobay

        # recolección por etapas
        times_b, blocks_b = [], []
        times_n, blocks_n = [], []

        for patt in stage_patterns:
            nc_b, nc_n = None, None
            try:
                nc_b = _open_mnc_or_raise(outdir_bay, patt)
                nc_n = _open_mnc_or_raise(outdir_nobay, patt)
            except FileNotFoundError:
                if nc_b is not None:
                    _close_mnc(nc_b)
                if nc_n is not None:
                    _close_mnc(nc_n)
                continue

            try:
                tb = _get_time_from_nc(nc_b)
                tn = _get_time_from_nc(nc_n)
                Xb = nc_b.variables[varname]  # (t,z,y,x)
                Xn = nc_n.variables[varname]

                if orient == "yz":
                    Xb_sec = np.array(Xb[:, k_slice, j0:j1, iG], dtype=float)
                    Xn_sec = np.array(Xn[:, k_slice, j0:j1, iG], dtype=float)
                else:
                    Xb_sec = np.array(Xb[:, k_slice, jG, i0:i1], dtype=float)
                    Xn_sec = np.array(Xn[:, k_slice, jG, i0:i1], dtype=float)

                if Xb_sec.shape[0] != tb.shape[0] or Xn_sec.shape[0] != tn.shape[0]:
                    raise ValueError(f"Tiempo/serie no coinciden en etapa '{patt}'.")

                # Validación por etapa: bay y nobay deben tener mismo eje temporal
                if tb.shape != tn.shape or np.any(np.abs(tb - tn) > 1e-7):
                    raise ValueError(f"Ejes de tiempo bay/nobay difieren en etapa '{patt}'.")

                times_b.append(tb)
                blocks_b.append(Xb_sec)
                times_n.append(tn)
                blocks_n.append(Xn_sec)
            finally:
                _close_mnc(nc_b)
                _close_mnc(nc_n)

        if not times_b or not times_n:
            raise FileNotFoundError(
                f"No se pudieron abrir etapas con patrones {stage_patterns} en:\n"
                f"  bay={outdir_bay}\n  nobay={outdir_nobay}"
            )

        t_b, Xb_all = _concat_time_and_data(times_b, blocks_b)
        t_n, Xn_all = _concat_time_and_data(times_n, blocks_n)
        if t_b.shape != t_n.shape or np.any(np.abs(t_b - t_n) > 1e-7):
            raise ValueError("Ejes de tiempo concatenados de bay y nobay no coinciden.")

        Xdiff = Xb_all - Xn_all  # (nt_full, nz_sel, nlat)

        # recorte spin-up
        tmin_s = tmin_horas * 3600.0
        idx0 = np.where(t_b >= tmin_s)[0][0] if np.any(t_b >= tmin_s) else 0
        t_cut = t_b[idx0:]
        anom = Xdiff[idx0:, :, :]
        if apply_mask_to_anom:
            anom = np.where(mask2d[None, :, :], anom, np.nan)

        idx_100 = int(np.argmin(np.abs(RC_sel + 100.0)))
        idx_400 = int(np.argmin(np.abs(RC_sel + 400.0)))

        # marcas de estaciones sobre el corte
        mark_wanted = {"Gamma", "C", "Omega", "Seiche2", "Seiche"}
        station_marks = []
        for s in stations_meta:
            lab = str(s["label"])
            if lab not in mark_wanted:
                continue
            i, j = int(s["i"]), int(s["j"])
            if orient == "yz":
                if i == iG and j0 <= j < j1:
                    station_marks.append((lab, float(YC[j])))
            else:
                if j == jG and i0 <= i < i1:
                    station_marks.append((lab, float(XC[i])))
        station_marks.sort(key=lambda tup: tup[1])
        station_labels = np.array([m[0] for m in station_marks], dtype="U32")
        station_coords = np.array([m[1] for m in station_marks], dtype=np.float64)

        meta = {
            "scenario": scenario_label,
            "var": varname,
            "orient": orient,
            "coord_name": coord_name,
            "start_station": start_station_label,
            "tmin_horas": float(tmin_horas),
            "RC_cut": None if RC_cut is None else float(RC_cut),
            "lateral_len": None if lateral_len is None else int(lateral_len),
            "idx0": int(idx0),
            "idx_100": int(idx_100),
            "idx_400": int(idx_400),
            "mask_grid_source": mask_grid_source,
            "apply_mask_to_anom": bool(apply_mask_to_anom),
        }
        out_path = save_dir / f"{scenario_label}_secThetaAnom_{orient}_{start_station_label}.npz"
        np.savez_compressed(
            out_path,
            t_horas=(t_cut / 3600.0).astype(np.float64),
            RC=RC_sel.astype(np.float64),
            coord=coord_full.astype(np.float64),
            anom=anom.astype(np.float32),
            mask2d=mask2d.astype(bool),
            station_labels=station_labels,
            station_coords=station_coords,
            station_marks=np.array(station_marks, dtype=object),
            meta_json=np.array(json.dumps(meta, ensure_ascii=False)),
        )
        print(f"[OK] Sección guardada: {out_path}")
        return str(out_path)
    finally:
        if gr_nobay is not None:
            _close_mnc(gr_nobay)
        _close_mnc(gr)

def guardar_seccion_v_anom_yz_npz(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    start_station_label: str = "Gamma",
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern: str = "grid.t*.nc",
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    incluir_gamma_omega: bool = True,
    RC_cut: float | None = None,
    lateral_len: int | None = None,
    mask_grid_source: str = "bay",     # "bay" | "nobay" | "intersection" | "union"
    apply_mask_to_anom: bool = False,
    interpolate_to_cgrid: bool = True,
) -> str:
    """
    Guarda sección vertical YZ de ΔV = V(bay) - V(nobay), con soporte para
    interpolar V(T,Z,Yp1,X) a malla C en Y:
      Vc(T,Z,Y,X) = 0.5 * [V(T,Z,Yp1_j,X) + V(T,Z,Yp1_{j+1},X)]

    Salida NPZ compatible con funciones de plot de secciones:
      - t_horas, RC, coord, anom, mask2d, station_labels, station_coords, meta_json
    """
    save_dir.mkdir(parents=True, exist_ok=True)

    def _get_var_any(nc_obj, candidates, label):
        for name in candidates:
            if name in nc_obj.variables:
                return nc_obj.variables[name]
        raise KeyError(f"No encontré {label}. Probé: {candidates}")

    def _coord_1d_from_var(var_obj, axis: str):
        arr = np.asarray(var_obj[:], dtype=float)
        if arr.ndim == 1:
            return arr
        if arr.ndim >= 2:
            return arr[0, :] if axis == "x" else arr[:, 0]
        raise ValueError("No pude convertir coordenada a 1D.")

    if mask_grid_source not in ("bay", "nobay", "intersection", "union"):
        raise ValueError("mask_grid_source debe ser 'bay', 'nobay', 'intersection' o 'union'.")

    gr = _open_mnc_or_raise(outdir_bay, grid_pattern)
    gr_nobay = None
    if mask_grid_source in ("nobay", "intersection", "union"):
        gr_nobay = _open_mnc_or_raise(outdir_nobay, grid_pattern)

    try:
        # Coordenadas: X de centros, Y de centros y Y de V (Yp1/YG)
        Xv = _coord_1d_from_var(
            _get_var_any(gr, ("X", "XC", "xC"), "coordenada X de centros"), axis="x"
        )
        Yc = _coord_1d_from_var(
            _get_var_any(gr, ("Y", "YC", "yC"), "coordenada Y de centros"), axis="y"
        )
        Yv = _coord_1d_from_var(
            _get_var_any(gr, ("Yp1", "YG", "yG", "YG"), "coordenada Y de V (Yp1/YG)"), axis="y"
        )
        RC = np.asarray(
            _get_var_any(gr, ("RC", "rC", "Z"), "coordenada vertical (RC/rC/Z)")[:], dtype=float
        ).squeeze()
        if RC.ndim > 1:
            RC = RC[..., 0, 0]

        # Máscaras
        hS_bay = np.asarray(
            _get_var_any(gr, ("hFacS", "HFacS"), "máscara hFacS/HFacS")[:], dtype=float
        )
        hC_bay = np.asarray(
            _get_var_any(gr, ("hFacC", "HFacC"), "máscara hFacC/HFacC")[:], dtype=float
        )
        hS_nobay = None
        hC_nobay = None
        if gr_nobay is not None:
            hS_nobay = np.asarray(
                _get_var_any(gr_nobay, ("hFacS", "HFacS"), "máscara hFacS/HFacS en nobay")[:], dtype=float
            )
            hC_nobay = np.asarray(
                _get_var_any(gr_nobay, ("hFacC", "HFacC"), "máscara hFacC/HFacC en nobay")[:], dtype=float
            )
            if hS_nobay.shape != hS_bay.shape:
                raise ValueError(
                    f"hFacS bay y nobay tienen forma distinta: {hS_bay.shape} vs {hS_nobay.shape}"
                )
            if hC_nobay.shape != hC_bay.shape:
                raise ValueError(
                    f"hFacC bay y nobay tienen forma distinta: {hC_bay.shape} vs {hC_nobay.shape}"
                )

        stations_meta = _make_stations_from_grid_inline(
            gr, j_centro=j_centro, radio=radio, incluir_gamma_omega=incluir_gamma_omega
        )
        labs = [s["label"] for s in stations_meta]
        if start_station_label not in labs:
            raise KeyError(f"No hallé la estación '{start_station_label}' en {labs}")
        st0 = stations_meta[labs.index(start_station_label)]
        i_fix, j_start = int(st0["i"]), int(st0["j"])

        j0 = max(0, j_start)
        if interpolate_to_cgrid:
            ny = Yc.size
            j1 = ny if lateral_len is None else min(ny, j0 + max(1, int(lateral_len)))
            coord_full = Yc[j0:j1]
            if (j1 + 1) > Yv.size:
                raise ValueError(
                    f"No puedo interpolar V a malla C: necesito Yp1 hasta {j1+1}, pero Yp1 size={Yv.size}."
                )
        else:
            nyv = Yv.size
            j1 = nyv if lateral_len is None else min(nyv, j0 + max(1, int(lateral_len)))
            coord_full = Yv[j0:j1]

        if RC_cut is not None:
            kz_arr = np.where(RC > RC_cut)[0]
            if kz_arr.size == 0:
                raise ValueError(f"No hay niveles con RC > {RC_cut} m.")
            k0, k1 = int(kz_arr.min()), int(kz_arr.max()) + 1
        else:
            k0, k1 = 0, int(RC.size)
        k_slice = slice(k0, k1)
        RC_sel = RC[k_slice]

        if interpolate_to_cgrid:
            wet_bay = hC_bay[k_slice, j0:j1, i_fix] > 0.0
            wet_nobay = (hC_nobay[k_slice, j0:j1, i_fix] > 0.0) if hC_nobay is not None else None
        else:
            wet_bay = hS_bay[k_slice, j0:j1, i_fix] > 0.0
            wet_nobay = (hS_nobay[k_slice, j0:j1, i_fix] > 0.0) if hS_nobay is not None else None
        if mask_grid_source == "bay":
            mask2d = wet_bay
        elif mask_grid_source == "nobay":
            mask2d = wet_nobay
        elif mask_grid_source == "intersection":
            mask2d = wet_bay & wet_nobay
        else:
            mask2d = wet_bay | wet_nobay

        times_b, blocks_b, times_n, blocks_n = [], [], [], []
        for patt in stage_patterns:
            nc_b, nc_n = None, None
            try:
                nc_b = _open_mnc_or_raise(outdir_bay, patt)
                nc_n = _open_mnc_or_raise(outdir_nobay, patt)
            except FileNotFoundError:
                if nc_b is not None:
                    _close_mnc(nc_b)
                if nc_n is not None:
                    _close_mnc(nc_n)
                continue

            try:
                tb = _get_time_from_nc(nc_b)
                tn = _get_time_from_nc(nc_n)
                Vb = nc_b.variables["V"]
                Vn = nc_n.variables["V"]

                if interpolate_to_cgrid:
                    # V(T,Z,Yp1,X) -> Vc(T,Z,Y,X) por promedio entre caras en Y.
                    Vb_l = np.asarray(Vb[:, k_slice, j0:j1, i_fix], dtype=float)
                    Vb_r = np.asarray(Vb[:, k_slice, j0 + 1:j1 + 1, i_fix], dtype=float)
                    Vn_l = np.asarray(Vn[:, k_slice, j0:j1, i_fix], dtype=float)
                    Vn_r = np.asarray(Vn[:, k_slice, j0 + 1:j1 + 1, i_fix], dtype=float)
                    Vb_sec = 0.5 * (Vb_l + Vb_r)
                    Vn_sec = 0.5 * (Vn_l + Vn_r)
                else:
                    Vb_sec = np.asarray(Vb[:, k_slice, j0:j1, i_fix], dtype=float)
                    Vn_sec = np.asarray(Vn[:, k_slice, j0:j1, i_fix], dtype=float)
                if Vb_sec.shape[0] != tb.shape[0] or Vn_sec.shape[0] != tn.shape[0]:
                    raise ValueError(f"Tiempo/serie no coinciden en etapa '{patt}'.")
                if tb.shape != tn.shape or np.any(np.abs(tb - tn) > 1e-7):
                    raise ValueError(f"Ejes de tiempo bay/nobay difieren en etapa '{patt}'.")

                times_b.append(tb); blocks_b.append(Vb_sec)
                times_n.append(tn); blocks_n.append(Vn_sec)
            finally:
                _close_mnc(nc_b)
                _close_mnc(nc_n)

        if not times_b or not times_n:
            raise FileNotFoundError(
                f"No se pudieron abrir etapas con patrones {stage_patterns} en:\n"
                f"  bay={outdir_bay}\n  nobay={outdir_nobay}"
            )

        t_b, Vb_all = _concat_time_and_data(times_b, blocks_b)
        t_n, Vn_all = _concat_time_and_data(times_n, blocks_n)
        if t_b.shape != t_n.shape or np.any(np.abs(t_b - t_n) > 1e-7):
            raise ValueError("Ejes de tiempo concatenados de bay y nobay no coinciden.")

        anom = Vb_all - Vn_all
        tmin_s = tmin_horas * 3600.0
        idx0 = np.where(t_b >= tmin_s)[0][0] if np.any(t_b >= tmin_s) else 0
        t_cut = t_b[idx0:]
        anom = anom[idx0:, :, :]
        if apply_mask_to_anom:
            anom = np.where(mask2d[None, :, :], anom, np.nan)

        # marcas en el transecto
        mark_wanted = {"Gamma", "C", "Omega", "Seiche2", "Seiche"}
        station_marks = []
        for s in stations_meta:
            lab = str(s["label"])
            if lab not in mark_wanted:
                continue
            i_s, j_s = int(s["i"]), int(s["j"])
            if i_s == i_fix and j0 <= j_s < j1:
                station_marks.append((lab, float(coord_full[j_s - j0])))
        station_marks.sort(key=lambda tup: tup[1])
        station_labels = np.array([m[0] for m in station_marks], dtype="U32")
        station_coords = np.array([m[1] for m in station_marks], dtype=np.float64)

        meta = {
            "scenario": scenario_label,
            "var": "V",
            "orient": "yz",
            "coord_name": "Y" if interpolate_to_cgrid else "Yp1",
            "start_station": start_station_label,
            "tmin_horas": float(tmin_horas),
            "RC_cut": None if RC_cut is None else float(RC_cut),
            "lateral_len": None if lateral_len is None else int(lateral_len),
            "idx0": int(idx0),
            "mask_grid_source": mask_grid_source,
            "apply_mask_to_anom": bool(apply_mask_to_anom),
            "interpolate_to_cgrid": bool(interpolate_to_cgrid),
        }

        out_path = save_dir / f"{scenario_label}_secVAnom_yz_{start_station_label}.npz"
        np.savez_compressed(
            out_path,
            t_horas=(t_cut / 3600.0).astype(np.float64),
            RC=RC_sel.astype(np.float64),
            coord=coord_full.astype(np.float64),
            anom=anom.astype(np.float32),
            mask2d=mask2d.astype(bool),
            station_labels=station_labels,
            station_coords=station_coords,
            station_marks=np.array(station_marks, dtype=object),
            meta_json=np.array(json.dumps(meta, ensure_ascii=False)),
        )
        print(f"[OK] Sección V anómala guardada: {out_path}")
        return str(out_path)
    finally:
        if gr_nobay is not None:
            _close_mnc(gr_nobay)
        _close_mnc(gr)

def guardar_seccion_w_anom_yz_npz(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    start_station_label: str = "Gamma",
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern: str = "grid.t*.nc",
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    incluir_gamma_omega: bool = True,
    RC_cut: float | None = None,
    lateral_len: int | None = None,
    mask_grid_source: str = "bay",     # "bay" | "nobay" | "intersection" | "union"
    apply_mask_to_anom: bool = False,
) -> str:
    """
    Guarda sección vertical YZ de ΔW = W(bay) - W(nobay).

    W vive en (T, Zl, Y, X): velocidad vertical en la cara superior de cada
    celda (interfaces Zl). Se colola a centros de celda (RC, 50 niveles) por
    promedio entre caras, conservando los 50 niveles (la celda más profunda usa
    la cara de fondo w=0):
      Wc[k] = 0.5*(W[k] + W[k+1])   k = 0..nz-2
      Wc[nz-1] = 0.5*W[nz-1]        (cara de fondo nula)

    Horizontalmente W está en malla C (mismo Y,X que Temp), así que se usa la
    coordenada de centros Yc y la máscara hFacC, igual que la sección de Temp.

    Salida NPZ compatible con funciones de plot de secciones:
      - t_horas, RC, coord, anom, mask2d, station_labels, station_coords, meta_json
    """
    save_dir.mkdir(parents=True, exist_ok=True)

    def _get_var_any(nc_obj, candidates, label):
        for name in candidates:
            if name in nc_obj.variables:
                return nc_obj.variables[name]
        raise KeyError(f"No encontré {label}. Probé: {candidates}")

    def _coord_1d_from_var(var_obj, axis: str):
        arr = np.asarray(var_obj[:], dtype=float)
        if arr.ndim == 1:
            return arr
        if arr.ndim >= 2:
            return arr[0, :] if axis == "x" else arr[:, 0]
        raise ValueError("No pude convertir coordenada a 1D.")

    def _collocate_w_keep_nz(Wslice):
        """Wslice (nt, nz, nlat) en caras Zl -> centros, conservando nz niveles."""
        Wc = np.empty_like(Wslice)
        Wc[:, :-1, :] = 0.5 * (Wslice[:, :-1, :] + Wslice[:, 1:, :])
        Wc[:, -1, :] = 0.5 * Wslice[:, -1, :]
        return Wc

    if mask_grid_source not in ("bay", "nobay", "intersection", "union"):
        raise ValueError("mask_grid_source debe ser 'bay', 'nobay', 'intersection' o 'union'.")

    gr = _open_mnc_or_raise(outdir_bay, grid_pattern)
    gr_nobay = None
    if mask_grid_source in ("nobay", "intersection", "union"):
        gr_nobay = _open_mnc_or_raise(outdir_nobay, grid_pattern)

    try:
        # Coordenadas: X de centros, Y de centros y vertical RC (centros)
        Xv = _coord_1d_from_var(
            _get_var_any(gr, ("X", "XC", "xC"), "coordenada X de centros"), axis="x"
        )
        Yc = _coord_1d_from_var(
            _get_var_any(gr, ("Y", "YC", "yC"), "coordenada Y de centros"), axis="y"
        )
        RC = np.asarray(
            _get_var_any(gr, ("RC", "rC", "Z"), "coordenada vertical (RC/rC/Z)")[:], dtype=float
        ).squeeze()
        if RC.ndim > 1:
            RC = RC[..., 0, 0]

        # Máscara de centros (hFacC)
        hC_bay = np.asarray(
            _get_var_any(gr, ("hFacC", "HFacC"), "máscara hFacC/HFacC")[:], dtype=float
        )
        hC_nobay = None
        if gr_nobay is not None:
            hC_nobay = np.asarray(
                _get_var_any(gr_nobay, ("hFacC", "HFacC"), "máscara hFacC/HFacC en nobay")[:], dtype=float
            )
            if hC_nobay.shape != hC_bay.shape:
                raise ValueError(
                    f"hFacC bay y nobay tienen forma distinta: {hC_bay.shape} vs {hC_nobay.shape}"
                )

        stations_meta = _make_stations_from_grid_inline(
            gr, j_centro=j_centro, radio=radio, incluir_gamma_omega=incluir_gamma_omega
        )
        labs = [s["label"] for s in stations_meta]
        if start_station_label not in labs:
            raise KeyError(f"No hallé la estación '{start_station_label}' en {labs}")
        st0 = stations_meta[labs.index(start_station_label)]
        i_fix, j_start = int(st0["i"]), int(st0["j"])

        j0 = max(0, j_start)
        ny = Yc.size
        j1 = ny if lateral_len is None else min(ny, j0 + max(1, int(lateral_len)))
        coord_full = Yc[j0:j1]

        if RC_cut is not None:
            kz_arr = np.where(RC > RC_cut)[0]
            if kz_arr.size == 0:
                raise ValueError(f"No hay niveles con RC > {RC_cut} m.")
            k0, k1 = int(kz_arr.min()), int(kz_arr.max()) + 1
        else:
            k0, k1 = 0, int(RC.size)
        k_slice = slice(k0, k1)
        RC_sel = RC[k_slice]

        wet_bay = hC_bay[k_slice, j0:j1, i_fix] > 0.0
        wet_nobay = (hC_nobay[k_slice, j0:j1, i_fix] > 0.0) if hC_nobay is not None else None
        if mask_grid_source == "bay":
            mask2d = wet_bay
        elif mask_grid_source == "nobay":
            mask2d = wet_nobay
        elif mask_grid_source == "intersection":
            mask2d = wet_bay & wet_nobay
        else:
            mask2d = wet_bay | wet_nobay

        times_b, blocks_b, times_n, blocks_n = [], [], [], []
        for patt in stage_patterns:
            nc_b, nc_n = None, None
            try:
                nc_b = _open_mnc_or_raise(outdir_bay, patt)
                nc_n = _open_mnc_or_raise(outdir_nobay, patt)
            except FileNotFoundError:
                if nc_b is not None:
                    _close_mnc(nc_b)
                if nc_n is not None:
                    _close_mnc(nc_n)
                continue

            try:
                tb = _get_time_from_nc(nc_b)
                tn = _get_time_from_nc(nc_n)
                Wb = nc_b.variables["W"]
                Wn = nc_n.variables["W"]

                # W(T,Zl,Y,X) en columna i_fix -> colola a centros conservando nz.
                Wb_slice = np.asarray(Wb[:, k_slice, j0:j1, i_fix], dtype=float)
                Wn_slice = np.asarray(Wn[:, k_slice, j0:j1, i_fix], dtype=float)
                Wb_sec = _collocate_w_keep_nz(Wb_slice)
                Wn_sec = _collocate_w_keep_nz(Wn_slice)

                if Wb_sec.shape[0] != tb.shape[0] or Wn_sec.shape[0] != tn.shape[0]:
                    raise ValueError(f"Tiempo/serie no coinciden en etapa '{patt}'.")
                if tb.shape != tn.shape or np.any(np.abs(tb - tn) > 1e-7):
                    raise ValueError(f"Ejes de tiempo bay/nobay difieren en etapa '{patt}'.")

                times_b.append(tb); blocks_b.append(Wb_sec)
                times_n.append(tn); blocks_n.append(Wn_sec)
            finally:
                _close_mnc(nc_b)
                _close_mnc(nc_n)

        if not times_b or not times_n:
            raise FileNotFoundError(
                f"No se pudieron abrir etapas con patrones {stage_patterns} en:\n"
                f"  bay={outdir_bay}\n  nobay={outdir_nobay}"
            )

        t_b, Wb_all = _concat_time_and_data(times_b, blocks_b)
        t_n, Wn_all = _concat_time_and_data(times_n, blocks_n)
        if t_b.shape != t_n.shape or np.any(np.abs(t_b - t_n) > 1e-7):
            raise ValueError("Ejes de tiempo concatenados de bay y nobay no coinciden.")

        anom = Wb_all - Wn_all
        tmin_s = tmin_horas * 3600.0
        idx0 = np.where(t_b >= tmin_s)[0][0] if np.any(t_b >= tmin_s) else 0
        t_cut = t_b[idx0:]
        anom = anom[idx0:, :, :]
        if apply_mask_to_anom:
            anom = np.where(mask2d[None, :, :], anom, np.nan)

        # marcas en el transecto
        mark_wanted = {"Gamma", "C", "Omega", "Seiche2", "Seiche"}
        station_marks = []
        for s in stations_meta:
            lab = str(s["label"])
            if lab not in mark_wanted:
                continue
            i_s, j_s = int(s["i"]), int(s["j"])
            if i_s == i_fix and j0 <= j_s < j1:
                station_marks.append((lab, float(coord_full[j_s - j0])))
        station_marks.sort(key=lambda tup: tup[1])
        station_labels = np.array([m[0] for m in station_marks], dtype="U32")
        station_coords = np.array([m[1] for m in station_marks], dtype=np.float64)

        meta = {
            "scenario": scenario_label,
            "var": "W",
            "orient": "yz",
            "coord_name": "Y",
            "start_station": start_station_label,
            "tmin_horas": float(tmin_horas),
            "RC_cut": None if RC_cut is None else float(RC_cut),
            "lateral_len": None if lateral_len is None else int(lateral_len),
            "idx0": int(idx0),
            "mask_grid_source": mask_grid_source,
            "apply_mask_to_anom": bool(apply_mask_to_anom),
            "collocate_z_to_center": True,
        }

        out_path = save_dir / f"{scenario_label}_secWAnom_yz_{start_station_label}.npz"
        np.savez_compressed(
            out_path,
            t_horas=(t_cut / 3600.0).astype(np.float64),
            RC=RC_sel.astype(np.float64),
            coord=coord_full.astype(np.float64),
            anom=anom.astype(np.float32),
            mask2d=mask2d.astype(bool),
            station_labels=station_labels,
            station_coords=station_coords,
            station_marks=np.array(station_marks, dtype=object),
            meta_json=np.array(json.dumps(meta, ensure_ascii=False)),
        )
        print(f"[OK] Sección W anómala guardada: {out_path}")
        return str(out_path)
    finally:
        if gr_nobay is not None:
            _close_mnc(gr_nobay)
        _close_mnc(gr)

def guardar_seccion_u_anom_yz_npz(
    scenario_label: str,
    outdir_bay: Path,
    outdir_nobay: Path,
    start_station_label: str = "Gamma",
    tmin_horas: float = 0.0,
    save_dir: Path = Path("./npz_out"),
    stage_patterns = ("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern: str = "grid.t*.nc",
    j_centro: int = _J_CENTRO_DEFAULT,
    radio: int = _RADIO_DEFAULT,
    incluir_gamma_omega: bool = True,
    RC_cut: float | None = None,
    lateral_len: int | None = None,
    mask_grid_source: str = "bay",     # "bay" | "nobay" | "intersection" | "union"
    apply_mask_to_anom: bool = False,
    interpolate_to_cgrid: bool = True,
) -> str:
    """
    Guarda sección vertical YZ de ΔU = U(bay) - U(nobay), con soporte para
    interpolar U(T,Z,Y,Xp1) a malla C en X:
      Uc(T,Z,Y,X) = 0.5 * [U(T,Z,Y,Xp1_i) + U(T,Z,Y,Xp1_{i+1})]

    Salida NPZ compatible con funciones de plot de secciones:
      - t_horas, RC, coord, anom, mask2d, station_labels, station_coords, meta_json
    """
    save_dir.mkdir(parents=True, exist_ok=True)

    def _get_var_any(nc_obj, candidates, label):
        for name in candidates:
            if name in nc_obj.variables:
                return nc_obj.variables[name]
        raise KeyError(f"No encontré {label}. Probé: {candidates}")

    def _coord_1d_from_var(var_obj, axis: str):
        arr = np.asarray(var_obj[:], dtype=float)
        if arr.ndim == 1:
            return arr
        if arr.ndim >= 2:
            return arr[0, :] if axis == "x" else arr[:, 0]
        raise ValueError("No pude convertir coordenada a 1D.")

    if mask_grid_source not in ("bay", "nobay", "intersection", "union"):
        raise ValueError("mask_grid_source debe ser 'bay', 'nobay', 'intersection' o 'union'.")

    gr = _open_mnc_or_raise(outdir_bay, grid_pattern)
    gr_nobay = None
    if mask_grid_source in ("nobay", "intersection", "union"):
        gr_nobay = _open_mnc_or_raise(outdir_nobay, grid_pattern)

    try:
        # Coordenadas: Y de centros, X de centros y X de U (Xp1/XG)
        Yc = _coord_1d_from_var(
            _get_var_any(gr, ("Y", "YC", "yC"), "coordenada Y de centros"), axis="y"
        )
        Xc = _coord_1d_from_var(
            _get_var_any(gr, ("X", "XC", "xC"), "coordenada X de centros"), axis="x"
        )
        Xu = _coord_1d_from_var(
            _get_var_any(gr, ("Xp1", "XG", "xG"), "coordenada X de U (Xp1/XG)"), axis="x"
        )
        RC = np.asarray(
            _get_var_any(gr, ("RC", "rC", "Z"), "coordenada vertical (RC/rC/Z)")[:], dtype=float
        ).squeeze()
        if RC.ndim > 1:
            RC = RC[..., 0, 0]

        # Máscaras
        hW_bay = np.asarray(
            _get_var_any(gr, ("hFacW", "HFacW"), "máscara hFacW/HFacW")[:], dtype=float
        )
        hC_bay = np.asarray(
            _get_var_any(gr, ("hFacC", "HFacC"), "máscara hFacC/HFacC")[:], dtype=float
        )
        hW_nobay = None
        hC_nobay = None
        if gr_nobay is not None:
            hW_nobay = np.asarray(
                _get_var_any(gr_nobay, ("hFacW", "HFacW"), "máscara hFacW/HFacW en nobay")[:], dtype=float
            )
            hC_nobay = np.asarray(
                _get_var_any(gr_nobay, ("hFacC", "HFacC"), "máscara hFacC/HFacC en nobay")[:], dtype=float
            )
            if hW_nobay.shape != hW_bay.shape:
                raise ValueError(
                    f"hFacW bay y nobay tienen forma distinta: {hW_bay.shape} vs {hW_nobay.shape}"
                )
            if hC_nobay.shape != hC_bay.shape:
                raise ValueError(
                    f"hFacC bay y nobay tienen forma distinta: {hC_bay.shape} vs {hC_nobay.shape}"
                )

        stations_meta = _make_stations_from_grid_inline(
            gr, j_centro=j_centro, radio=radio, incluir_gamma_omega=incluir_gamma_omega
        )
        labs = [s["label"] for s in stations_meta]
        if start_station_label not in labs:
            raise KeyError(f"No hallé la estación '{start_station_label}' en {labs}")
        st0 = stations_meta[labs.index(start_station_label)]
        i_start, j_center = int(st0["i"]), int(st0["j"])

        j0 = max(0, j_center)
        ny = Yc.size
        j1 = ny if lateral_len is None else min(ny, j0 + max(1, int(lateral_len)))
        coord_full = Yc[j0:j1]

        if RC_cut is not None:
            kz_arr = np.where(RC > RC_cut)[0]
            if kz_arr.size == 0:
                raise ValueError(f"No hay niveles con RC > {RC_cut} m.")
            k0, k1 = int(kz_arr.min()), int(kz_arr.max()) + 1
        else:
            k0, k1 = 0, int(RC.size)
        k_slice = slice(k0, k1)
        RC_sel = RC[k_slice]

        # Índice en X
        i_fix = int(i_start)
        if interpolate_to_cgrid and (i_fix + 1) > (Xu.size - 1):
            raise ValueError(
                f"No puedo interpolar U a malla C: necesito Xp1 hasta {i_fix+1}, pero Xp1 size={Xu.size}."
            )

        if interpolate_to_cgrid:
            wet_bay = hC_bay[k_slice, j0:j1, i_fix] > 0.0
            wet_nobay = (hC_nobay[k_slice, j0:j1, i_fix] > 0.0) if hC_nobay is not None else None
        else:
            wet_bay = hW_bay[k_slice, j0:j1, i_fix] > 0.0
            wet_nobay = (hW_nobay[k_slice, j0:j1, i_fix] > 0.0) if hW_nobay is not None else None

        if mask_grid_source == "bay":
            mask2d = wet_bay
        elif mask_grid_source == "nobay":
            mask2d = wet_nobay
        elif mask_grid_source == "intersection":
            mask2d = wet_bay & wet_nobay
        else:
            mask2d = wet_bay | wet_nobay

        times_b, blocks_b, times_n, blocks_n = [], [], [], []
        for patt in stage_patterns:
            nc_b, nc_n = None, None
            try:
                nc_b = _open_mnc_or_raise(outdir_bay, patt)
                nc_n = _open_mnc_or_raise(outdir_nobay, patt)
            except FileNotFoundError:
                if nc_b is not None:
                    _close_mnc(nc_b)
                if nc_n is not None:
                    _close_mnc(nc_n)
                continue

            try:
                tb = _get_time_from_nc(nc_b)
                tn = _get_time_from_nc(nc_n)
                Ub = nc_b.variables["U"]
                Un = nc_n.variables["U"]

                if interpolate_to_cgrid:
                    Ub_l = np.asarray(Ub[:, k_slice, j0:j1, i_fix], dtype=float)
                    Ub_r = np.asarray(Ub[:, k_slice, j0:j1, i_fix + 1], dtype=float)
                    Un_l = np.asarray(Un[:, k_slice, j0:j1, i_fix], dtype=float)
                    Un_r = np.asarray(Un[:, k_slice, j0:j1, i_fix + 1], dtype=float)
                    Ub_sec = 0.5 * (Ub_l + Ub_r)
                    Un_sec = 0.5 * (Un_l + Un_r)
                else:
                    Ub_sec = np.asarray(Ub[:, k_slice, j0:j1, i_fix], dtype=float)
                    Un_sec = np.asarray(Un[:, k_slice, j0:j1, i_fix], dtype=float)

                if Ub_sec.shape[0] != tb.shape[0] or Un_sec.shape[0] != tn.shape[0]:
                    raise ValueError(f"Tiempo/serie no coinciden en etapa '{patt}'.")
                if tb.shape != tn.shape or np.any(np.abs(tb - tn) > 1e-7):
                    raise ValueError(f"Ejes de tiempo bay/nobay difieren en etapa '{patt}'.")

                times_b.append(tb); blocks_b.append(Ub_sec)
                times_n.append(tn); blocks_n.append(Un_sec)
            finally:
                _close_mnc(nc_b)
                _close_mnc(nc_n)

        if not times_b or not times_n:
            raise FileNotFoundError(
                f"No se pudieron abrir etapas con patrones {stage_patterns} en:\n"
                f"  bay={outdir_bay}\n  nobay={outdir_nobay}"
            )

        t_b, Ub_all = _concat_time_and_data(times_b, blocks_b)
        t_n, Un_all = _concat_time_and_data(times_n, blocks_n)
        if t_b.shape != t_n.shape or np.any(np.abs(t_b - t_n) > 1e-7):
            raise ValueError("Ejes de tiempo concatenados de bay y nobay no coinciden.")

        anom = Ub_all - Un_all
        tmin_s = tmin_horas * 3600.0
        idx0 = np.where(t_b >= tmin_s)[0][0] if np.any(t_b >= tmin_s) else 0
        t_cut = t_b[idx0:]
        anom = anom[idx0:, :, :]
        if apply_mask_to_anom:
            anom = np.where(mask2d[None, :, :], anom, np.nan)

        mark_wanted = {"Gamma", "C", "Omega", "Seiche2", "Seiche"}
        station_marks = []
        for s in stations_meta:
            lab = str(s["label"])
            if lab not in mark_wanted:
                continue
            i_s, j_s = int(s["i"]), int(s["j"])
            if i_s == i_fix and j0 <= j_s < j1:
                station_marks.append((lab, float(coord_full[j_s - j0])))
        station_marks.sort(key=lambda tup: tup[1])
        station_labels = np.array([m[0] for m in station_marks], dtype="U32")
        station_coords = np.array([m[1] for m in station_marks], dtype=np.float64)

        meta = {
            "scenario": scenario_label,
            "var": "U",
            "orient": "yz",
            "coord_name": "Y",
            "start_station": start_station_label,
            "tmin_horas": float(tmin_horas),
            "RC_cut": None if RC_cut is None else float(RC_cut),
            "lateral_len": None if lateral_len is None else int(lateral_len),
            "idx0": int(idx0),
            "mask_grid_source": mask_grid_source,
            "apply_mask_to_anom": bool(apply_mask_to_anom),
            "interpolate_to_cgrid": bool(interpolate_to_cgrid),
        }

        out_path = save_dir / f"{scenario_label}_secUAnom_yz_{start_station_label}.npz"
        np.savez_compressed(
            out_path,
            t_horas=(t_cut / 3600.0).astype(np.float64),
            RC=RC_sel.astype(np.float64),
            coord=coord_full.astype(np.float64),
            anom=anom.astype(np.float32),
            mask2d=mask2d.astype(bool),
            station_labels=station_labels,
            station_coords=station_coords,
            station_marks=np.array(station_marks, dtype=object),
            meta_json=np.array(json.dumps(meta, ensure_ascii=False)),
        )
        print(f"[OK] Sección U anómala guardada: {out_path}")
        return str(out_path)
    finally:
        if gr_nobay is not None:
            _close_mnc(gr_nobay)
        _close_mnc(gr)

# --------------------------- Animar desde el bundle ---------------------------
def _edges_1d(x):
    x = np.asarray(x, dtype=float)
    if x.size == 1:
        dx = 1.0
        return np.array([x[0]-0.5*dx, x[0]+0.5*dx])
    dx = np.diff(x)
    # bordes internos a la mitad
    edges_mid = x[:-1] + 0.5*dx
    # extremos
    x0 = x[0] - 0.5*dx[0]
    xN = x[-1] + 0.5*dx[-1]
    return np.concatenate(([x0], edges_mid, [xN]))

def animar_seccion_temp_anom_desde_npz(
    npz_path: str,
    outname: str = "anim_thetaAnom_section.mp4",
    vmin: float = -0.01,
    vmax: float = 0.01,
    fps: int = 2,
    draw_mask_contour: bool = True,
    stations_color: str = "k",
    stations_ls: str = "--",
    stations_lw: float = 1.2
):
    d = np.load(npz_path, allow_pickle=True)
    t_h    = d["t_horas"]            # (nt,)
    RC     = d["RC"]                 # (nz,)
    coord  = d["coord"]              # (nlat,)
    A      = d["anom"]               # (nt, nz, nlat)
    mask2d = d["mask2d"]             # (nz, nlat) boolean
    marks = []
    if "station_labels" in d.files and "station_coords" in d.files:
        labs = d["station_labels"]
        xs = d["station_coords"]
        marks = list(zip(labs.tolist(), xs.tolist()))
    elif "station_marks" in d.files:
        marks = d["station_marks"]
    meta   = json.loads(str(d["meta_json"]))
    orient = meta.get("orient","yz")
    coord_name = meta.get("coord_name","YC")
    scenario = meta.get("scenario","")
    start_lab= meta.get("start_station","Gamma")

    # --- Construir bordes para pcolormesh (robusto a versiones de Matplotlib) ---
    x_edges = _edges_1d(coord)   # len = nlat+1
    z_edges = _edges_1d(RC)      # len = nz+1

    fig, ax = plt.subplots(figsize=(10, 6))

    # ¡OJO!: Z debe ser (nz, nlat). Usamos A[0] SIN transponer.
    Z0 = A[0]  # (nz, nlat)
    mesh = ax.pcolormesh(
    coord, RC, A[0],
    cmap=cmo.cm.balance, vmin=vmin, vmax=vmax, shading="nearest"
)
    cbar = fig.colorbar(mesh, ax=ax, label="Anomalía Temp (°C)")

    if draw_mask_contour and mask2d.size:
        # mask2d ya está (nz, nlat); no transponer.
        ax.contour(coord, RC, (~mask2d).astype(float), colors="black", linewidths=2, linestyles="--")

    if len(marks):
        for lab, x0 in marks:
            x0f = float(x0)
            ax.axvline(x0f, color=stations_color, lw=stations_lw, ls=stations_ls)
            ax.text(x0f, RC.min(), f" {lab}", color=stations_color, va="bottom",
                    ha="left", fontsize=9, rotation=90, rotation_mode="anchor")

    title_base = f"{scenario} • ΔTemp (bay−nobay) • corte {orient.upper()} desde {start_lab}"
    ax.set_title(f"{title_base}  |  t={t_h[0]:.2f} h")
    ax.set_xlabel(coord_name)
    ax.set_ylabel("Profundidad (m)")
    plt.tight_layout()

    # Actualización: set_array requiere (nlat*nz,)
    nlat = coord.size
    nz   = RC.size
    def update(frame):
        Zf = A[frame]  # (nz, nlat)
        mesh.set_array(A[frame].ravel())
        ax.set_title(f"{title_base}  |  t={t_h[frame]:.2f} h")
        return mesh,

    anim = FuncAnimation(fig, update, frames=np.arange(t_h.size), blit=True)
    anim.save(outname, writer="ffmpeg", fps=fps)
    plt.show()
    return anim

# intento usar cmocean; si no está, cae a 'RdBu_r'
try:
    import cmocean as cmo
    _CMAP = cmo.cm.balance
except Exception:
    _CMAP = 'RdBu_r'

def graficar_seccion_temp_anom_en_tiempo(
    npz_path: str,
    t_objetivo_horas: float | None = None,
    frame_idx: int | None = None,
    vmin: float = -0.01,
    vmax: float = 0.01,
    normalize_by_max_abs: bool = False,
    normalize_scope: str = "frame",  # "frame" | "dataset"
    draw_mask_contour: bool = False,
    stations_color: str = "k",
    station_color_map: dict | None = None,
    stations_ls: str = "--",
    stations_lw: float = 1.2,
    figsize=(10, 6),
    titulo_extra: str | None = None,
):
    """
    Grafica la sección vertical de anomalía para un tiempo dado.

    Selección de frame:
      - Si `frame_idx` no es None, usa ese índice.
      - Si no, y `t_objetivo_horas` no es None, usa el frame más cercano a ese tiempo.
      - Si ambos son None, usa el frame 0.

    Normalización:
      - Si `normalize_by_max_abs=True`, divide por max(abs(.)).
      - `normalize_scope="frame"` usa solo el frame elegido.
      - `normalize_scope="dataset"` usa todo `anom`.

    Retorna (fig, ax, idx_usado, t_horas_usado).
    """
    d = np.load(npz_path, allow_pickle=True)
    t_h = d["t_horas"]       # (nt,)
    RC = d["RC"]             # (nz,)
    coord = d["coord"]       # (nlat,)
    A = d["anom"]            # (nt, nz, nlat)
    mask2d = d["mask2d"] if "mask2d" in d.files else np.array([])

    marks = []
    if "station_labels" in d.files and "station_coords" in d.files:
        marks = list(zip(d["station_labels"].tolist(), d["station_coords"].tolist()))
    elif "station_marks" in d.files:
        marks = d["station_marks"]

    try:
        meta = json.loads(str(d["meta_json"]))
    except Exception:
        meta = {}

    orient = meta.get("orient", "yz")
    coord_name = meta.get("coord_name", "coord")
    scenario = meta.get("scenario", "")
    start_lab = meta.get("start_station", "Gamma")

    if frame_idx is not None:
        idx = int(frame_idx)
    elif t_objetivo_horas is not None:
        idx = int(np.argmin(np.abs(t_h - float(t_objetivo_horas))))
    else:
        idx = 0
    idx = max(0, min(idx, t_h.size - 1))

    if normalize_scope not in ("frame", "dataset"):
        raise ValueError("normalize_scope debe ser 'frame' o 'dataset'.")

    Z_plot = np.array(A[idx], dtype=float)
    cbar_label = "Anomalía Temp (°C)"
    vmin_plot, vmax_plot = vmin, vmax
    norm_factor = 1.0
    if normalize_by_max_abs:
        if normalize_scope == "dataset":
            norm_factor = float(np.nanmax(np.abs(A)))
        else:
            norm_factor = float(np.nanmax(np.abs(Z_plot)))
        if not np.isfinite(norm_factor) or norm_factor == 0.0:
            norm_factor = 1.0
        Z_plot = Z_plot / norm_factor
        vmin_plot, vmax_plot = vmin / norm_factor, vmax / norm_factor
        cbar_label = "Anomalía Temp normalizada (|max|=1)"

    # Aplicar máscara de tierra para rellenarla en color arena.
    if mask2d.size:
        Z_masked = np.ma.masked_where(~mask2d, Z_plot)
    else:
        Z_masked = np.ma.masked_invalid(Z_plot)

    if not np.ma.is_masked(Z_masked):
        Z_masked = np.ma.masked_invalid(Z_masked)

    cmap_plot = plt.get_cmap(_CMAP).copy()
    cmap_plot.set_bad("#ece2c6")

    if np.isclose(vmax_plot, vmin_plot):
        dv = 1e-12 if vmax_plot == 0 else abs(vmax_plot) * 0.05
        vmin_plot, vmax_plot = vmin_plot - dv, vmax_plot + dv
    levels = np.linspace(vmin_plot, vmax_plot, 121)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_facecolor("#ece2c6")
    mesh = ax.contourf(coord, RC, Z_masked, levels=levels, cmap=cmap_plot, extend="both")
    fig.colorbar(mesh, ax=ax, label=cbar_label)

    if draw_mask_contour and mask2d.size:
        ax.contour(coord, RC, (~mask2d).astype(float), colors="black", linewidths=2, linestyles="--")

    if len(marks):
        default_station_colors = {
            "Seiche": "#9d02d7",
            "Seiche2": "#ffd700",
            "A": "#ea5f94",
            "B": "#0000ff",
            "C": "#fa8775",
            "Gamma": "limegreen",
            "Omega": "cyan",
        }
        if station_color_map:
            default_station_colors.update(station_color_map)
        for lab, x0 in marks:
            x0f = float(x0)
            st_color = default_station_colors.get(str(lab), stations_color)
            ax.axvline(x0f, color=st_color, lw=stations_lw, ls=stations_ls)
            ax.text(
                x0f,
                RC.min(),
                f" {lab}",
                color=st_color,
                va="bottom",
                ha="left",
                fontsize=9,
                rotation=90,
                rotation_mode="anchor",
            )

    title = f"{scenario} • ΔTemp (bay−nobay) • corte {orient.upper()} desde {start_lab}"
    title += f" | t={t_h[idx]:.2f} h (idx={idx})"
    if normalize_by_max_abs:
        title += f" | norm={normalize_scope} (max={norm_factor:.3e})"
    if titulo_extra:
        title += f"\n{titulo_extra}"
    ax.set_title(title)
    ax.set_xlabel(coord_name)
    ax.set_ylabel("Profundidad (m)")
    plt.tight_layout()
    plt.show()
    return fig, ax, idx, float(t_h[idx])

def graficar_seccion_temp_anom_multi_tiempos(
    npz_path: str,
    tiempos_horas: list[float] | tuple[float, ...] | np.ndarray | None = None,
    frame_indices: list[int] | tuple[int, ...] | np.ndarray | None = None,
    max_cols: int = 5,
    vmin: float = -0.01,
    vmax: float = 0.01,
    normalize_by_max_abs: bool = False,
    normalize_scope: str = "dataset",  # "frame" | "dataset"
    draw_stations: bool = True,
    stations_color: str = "k",
    station_color_map: dict | None = None,
    stations_ls: str = "--",
    stations_lw: float = 1.0,
    figsize_per_panel: tuple[float, float] = (3.7, 3.0),
    sharex: bool = True,
    sharey: bool = True,
    titulo_general: str | None = None,
):
    """
    Grafica múltiples tiempos de la sección en una grilla:
    máximo `max_cols` columnas y las filas necesarias.

    Debes pasar uno de:
      - tiempos_horas: se busca el frame más cercano para cada tiempo.
      - frame_indices: índices directos de frame.

    Retorna (fig, axes, indices_usados, tiempos_horas_usados).
    """
    d = np.load(npz_path, allow_pickle=True)
    t_h = d["t_horas"]
    RC = d["RC"]
    coord = d["coord"]
    A = d["anom"]
    mask2d = d["mask2d"] if "mask2d" in d.files else np.array([])

    marks = []
    if "station_labels" in d.files and "station_coords" in d.files:
        marks = list(zip(d["station_labels"].tolist(), d["station_coords"].tolist()))
    elif "station_marks" in d.files:
        marks = d["station_marks"]

    try:
        meta = json.loads(str(d["meta_json"]))
    except Exception:
        meta = {}

    orient = meta.get("orient", "yz")
    coord_name = meta.get("coord_name", "coord")
    scenario = meta.get("scenario", "")
    start_lab = meta.get("start_station", "Gamma")

    if frame_indices is None and tiempos_horas is None:
        raise ValueError("Debes pasar 'tiempos_horas' o 'frame_indices'.")
    if frame_indices is not None and tiempos_horas is not None:
        raise ValueError("Usa solo uno: 'tiempos_horas' o 'frame_indices'.")

    if frame_indices is not None:
        idx_list = [int(i) for i in np.ravel(np.asarray(frame_indices))]
        idx_list = [max(0, min(i, t_h.size - 1)) for i in idx_list]
    else:
        t_query = [float(x) for x in np.ravel(np.asarray(tiempos_horas))]
        idx_list = [int(np.argmin(np.abs(t_h - tq))) for tq in t_query]

    if len(idx_list) == 0:
        raise ValueError("No hay tiempos/índices para graficar.")

    if normalize_scope not in ("frame", "dataset"):
        raise ValueError("normalize_scope debe ser 'frame' o 'dataset'.")

    dataset_norm = 1.0
    if normalize_by_max_abs and normalize_scope == "dataset":
        dataset_norm = float(np.nanmax(np.abs(A)))
        if not np.isfinite(dataset_norm) or dataset_norm == 0.0:
            dataset_norm = 1.0

    n_panels = len(idx_list)
    ncols = max(1, min(int(max_cols), n_panels))
    nrows = int(np.ceil(n_panels / ncols))
    fig_w = figsize_per_panel[0] * ncols
    fig_h = figsize_per_panel[1] * nrows
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(fig_w, fig_h),
        squeeze=False,
        sharex=sharex,
        sharey=sharey,
    )
    axes_flat = axes.ravel()

    cmap_plot = plt.get_cmap(_CMAP).copy()
    cmap_plot.set_bad("#ece2c6")

    default_station_colors = {
        "Seiche": "#9d02d7",
        "Seiche2": "#ffd700",
        "A": "#ea5f94",
        "B": "#0000ff",
        "C": "#fa8775",
        "Gamma": "limegreen",
        "Omega": "cyan",
    }
    if station_color_map:
        default_station_colors.update(station_color_map)

    mappable = None
    t_used = []
    for ip, idx in enumerate(idx_list):
        ax = axes_flat[ip]
        Z_plot = np.array(A[idx], dtype=float)
        norm_factor = dataset_norm

        if normalize_by_max_abs and normalize_scope == "frame":
            norm_factor = float(np.nanmax(np.abs(Z_plot)))
            if not np.isfinite(norm_factor) or norm_factor == 0.0:
                norm_factor = 1.0

        if normalize_by_max_abs:
            Z_plot = Z_plot / norm_factor
            vmin_plot, vmax_plot = vmin / norm_factor, vmax / norm_factor
        else:
            vmin_plot, vmax_plot = vmin, vmax

        if np.isclose(vmax_plot, vmin_plot):
            dv = 1e-12 if vmax_plot == 0 else abs(vmax_plot) * 0.05
            vmin_plot, vmax_plot = vmin_plot - dv, vmax_plot + dv
        levels = np.linspace(vmin_plot, vmax_plot, 121)

        if mask2d.size:
            Z_masked = np.ma.masked_where(~mask2d, Z_plot)
        else:
            Z_masked = np.ma.masked_invalid(Z_plot)

        ax.set_facecolor("#ece2c6")
        mappable = ax.contourf(coord, RC, Z_masked, levels=levels, cmap=cmap_plot, extend="both")

        if draw_stations and len(marks):
            for lab, x0 in marks:
                st_color = default_station_colors.get(str(lab), stations_color)
                x0f = float(x0)
                ax.axvline(x0f, color=st_color, lw=stations_lw, ls=stations_ls)
                ax.text(
                    x0f, RC.min(), f" {lab}", color=st_color,
                    va="bottom", ha="left", fontsize=8, rotation=90, rotation_mode="anchor"
                )

        ax.set_title(f"t={t_h[idx]:.2f} h (idx={idx})", fontsize=10)
        ax.set_xlabel(coord_name)
        ax.set_ylabel("Profundidad (m)")
        t_used.append(float(t_h[idx]))

    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    title = titulo_general or f"{scenario} • ΔTemp (bay−nobay) • corte {orient.upper()} desde {start_lab}"
    if normalize_by_max_abs:
        title += f" | norm={normalize_scope}"
    fig.suptitle(title, fontsize=12)
    fig.subplots_adjust(top=0.90, wspace=0.24, hspace=0.34)

    if mappable is not None:
        cbar_label = "Anomalía Temp normalizada (|max|=1)" if normalize_by_max_abs else "Anomalía Temp (°C)"
        fig.colorbar(mappable, ax=axes_flat.tolist(), shrink=0.95, label=cbar_label)

    plt.show()
    return fig, axes, idx_list, t_used

def graficar_seccion_norm_por_estacion(
    npz_path: str,
    station_label: str = "C",
    step_frames: int = 60,
    max_cols: int = 3,
    levels: int = 121,
    vmin: float = -1.0,
    vmax: float = 1.0,
    cmap_name: str = "RdBu_r",
    land_color: str = "#ece2c6",
    sharex: bool = True,
    sharey: bool = True,
    figsize_panel: tuple[float, float] = (4.4, 3.2),
    cbar_only_col: int | None = None,   # 1-based; ej. 3 = solo tercera columna
    cbar_label: str | None = None,
    var_display: str | None = None,
):
    """
    Grafica paneles de sección normalizando por el máximo absoluto de una estación
    (todos los tiempos y profundidades), con máscara de tierra en color arena.

    Retorna:
      fig, axes, info (dict con idx_estacion, x_estacion, norm_factor, idxs)
    """
    with np.load(npz_path, allow_pickle=True) as d:
        A = d["anom"].astype(float)  # (nt, nz, nlat)
        coord = d["coord"]
        RC = d["RC"]
        t_h = d["t_horas"]
        mask2d = d["mask2d"] if "mask2d" in d.files else np.ones(A.shape[1:], dtype=bool)
        labels = d["station_labels"] if "station_labels" in d.files else np.array([], dtype=str)
        scoords = d["station_coords"] if "station_coords" in d.files else np.array([], dtype=float)
        meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}

    if step_frames <= 0:
        raise ValueError("step_frames debe ser > 0.")
    if max_cols <= 0:
        raise ValueError("max_cols debe ser > 0.")
    if mask2d.shape != A.shape[1:]:
        raise ValueError(f"mask2d {mask2d.shape} no coincide con campo {A.shape[1:]}.")

    # Buscar coordenada de estación para normalización; fallback al centro si no existe.
    i_station = A.shape[2] // 2
    x_station = float(coord[i_station])
    if labels.size and scoords.size:
        labels_list = [str(x) for x in labels.tolist()]
        if station_label in labels_list:
            x_station = float(scoords[labels_list.index(station_label)])
            i_station = int(np.argmin(np.abs(coord - x_station)))

    norm_factor = float(np.nanmax(np.abs(A[:, :, i_station])))
    if not np.isfinite(norm_factor) or norm_factor == 0.0:
        norm_factor = 1.0
    A_norm = A / norm_factor

    idxs = np.arange(0, A_norm.shape[0], step_frames, dtype=int)
    n = len(idxs)
    ncols = min(max_cols, n)
    nrows = int(np.ceil(n / ncols))
    cbar_col_idx = None
    if cbar_only_col is not None:
        cbar_col_idx = min(max(1, int(cbar_only_col)), ncols) - 1

    coord_km = coord / 1000.0
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_panel[0] * ncols, figsize_panel[1] * nrows),
        sharex=sharex, sharey=sharey
    )
    axes = np.atleast_1d(axes).ravel()

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(land_color)
    lvls = np.linspace(vmin, vmax, int(levels))
    var_name = str(var_display if var_display is not None else meta.get("var", "Temp"))
    cbar_label_txt = cbar_label if cbar_label is not None else f"Δ{var_name} norm."

    for k, idx in enumerate(idxs):
        ax = axes[k]
        ax.set_facecolor(land_color)
        Z_nan = np.where(mask2d, A_norm[idx], np.nan)

        cf = ax.contourf(
            coord_km, RC, Z_nan,
            levels=lvls, cmap=cmap, vmin=vmin, vmax=vmax, extend="both"
        )
        show_cbar = (cbar_col_idx is None) or ((k % ncols) == cbar_col_idx)
        if show_cbar:
            cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(cbar_label_txt, fontsize=8)
            cbar.ax.tick_params(labelsize=8)

        ax.set_title(f"t={t_h[idx]:.1f} h", fontsize=10)
        ax.set_xlabel(f'{meta.get("coord_name", "YC")} (km)')
        ax.set_ylabel("Profundidad (m)")

    for ax in axes[n:]:
        ax.axis("off")

    scenario = meta.get("scenario", "")
    fig.suptitle(
        f"{scenario} • Δ{var_name} norm por estación {station_label} (max abs global)",
        y=1.01, fontsize=12
    )
    fig.tight_layout()
    plt.show()

    info = {
        "idx_estacion": int(i_station),
        "x_estacion": float(x_station),
        "norm_factor": float(norm_factor),
        "idxs": idxs,
    }
    return fig, axes, info

def graficar_seccion_velocidad_norm_por_estacion(
    npz_path: str,
    station_label: str = "C",
    step_frames: int = 60,
    max_cols: int = 3,
    levels: int = 121,
    vmin: float = -1.0,
    vmax: float = 1.0,
    cmap_name: str = "RdBu_r",
    land_color: str = "#ece2c6",
    sharex: bool = True,
    sharey: bool = True,
    figsize_panel: tuple[float, float] = (4.4, 3.2),
    cbar_only_col: int = 3,
):
    """
    Wrapper para U/V normalizados por estación.
    Muestra colorbar solo en la columna indicada (por default, columna 3).
    """
    with np.load(npz_path, allow_pickle=True) as d:
        meta = json.loads(str(d["meta_json"])) if "meta_json" in d.files else {}
    var_name = str(meta.get("var", "Vel"))

    return graficar_seccion_norm_por_estacion(
        npz_path=npz_path,
        station_label=station_label,
        step_frames=step_frames,
        max_cols=max_cols,
        levels=levels,
        vmin=vmin,
        vmax=vmax,
        cmap_name=cmap_name,
        land_color=land_color,
        sharex=sharex,
        sharey=sharey,
        figsize_panel=figsize_panel,
        cbar_only_col=cbar_only_col,
        cbar_label=f"Δ{var_name} norm.",
        var_display=var_name,
    )

def animar_anomalia_section_simple(
    npz_path: str,
    outname: str = "anim_anomalia_section.mp4",
    vmin: float = -0.01,
    vmax: float = 0.01,
    fps: int = 2,
    titulo: str | None = None
):
    """
    Anima una sección vertical de anomalías (A: shape = (nt, nz, nlat)) usando imshow.

    Espera que el NPZ contenga:
      - 'anom'  : (nt, nz, nlat)
      - 't_horas' : (nt,)
      - 'RC'    : (nz,)   (profundidades, típicamente negativas)
      - 'coord' : (nlat,) (YC o XC según corte)

    Guardará un MP4 con nombre 'outname'.
    """
    d = np.load(npz_path, allow_pickle=True)
    A  = d["anom"]        # (nt, nz, nlat)
    t  = d["t_horas"]     # (nt,)
    RC = d["RC"]          # (nz,)
    X  = d["coord"]       # (nlat,)

    # Extensión de ejes para imshow (x = coord, y = RC)
    x0, x1 = float(np.min(X)), float(np.max(X))
    y0, y1 = float(np.min(RC)), float(np.max(RC))
    extent = (x0, x1, y0, y1)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(
        A[0],                    # (nz, nlat)
        origin='lower',          # luego invertimos eje Y para ver superficie arriba
        aspect='auto',
        interpolation='nearest',
        cmap=_CMAP,
        vmin=vmin, vmax=vmax,
        extent=extent
    )
    cbar = fig.colorbar(im, ax=ax, label="Anomalía Temp (°C)")

    # Superficie arriba (profundidad menos negativa arriba)
    ax.set_ylim(y1, y0)

    # Etiquetas y título
    try:
        import json
        meta = json.loads(str(d.get("meta_json", "{}")))
        coord_name = meta.get("coord_name", "coord")
        scenario   = meta.get("scenario", "")
    except Exception:
        coord_name, scenario = "coord", ""
    ax.set_xlabel(coord_name)
    ax.set_ylabel("Profundidad (m)")
    base_title = titulo if titulo is not None else f"{scenario} • ΔTemp (bay−nobay)"
    ax.set_title(f"{base_title} | t={t[0]:.2f} h")
    plt.tight_layout()

    def update(frame):
        im.set_data(A[frame])     # misma forma (nz, nlat)
        ax.set_title(f"{base_title} | t={t[frame]:.2f} h")
        return (im,)

    anim = FuncAnimation(fig, update, frames=np.arange(A.shape[0]), blit=True)
    anim.save(outname, writer="ffmpeg", fps=fps)
    plt.show()
    return anim

# ─────────────────────────────────────────────────────────────
# FUNCIÓN PRINCIPAL: calcular desplazamiento de isoterma
# ─────────────────────────────────────────────────────────────
def calcular_desplazamiento_isoterma_escenario(
    scenario_label: str,
    npz_bay_path,
    npz_nobay_path,
    T_iso: float = 7.1,
    t_ref_index: int = 0,
    save_dir="npz_isoterma",
    nombre_base: str | None = None,
):
    """
    Calcula profundidad y desplazamiento temporal de una isoterma
    para un escenario MITgcm (bahía / sin bahía).

    Parámetros
    ----------
    scenario_label : str
        Nombre del escenario ("realista", "linear", "2layer", etc.).
    npz_bay_path : Path o str
        Ruta del archivo NPZ con Hovmöller de Temp anómala ('bay').
    npz_nobay_path : Path o str
        Ruta del archivo NPZ correspondiente ('nobay').
    T_iso : float
        Temperatura objetivo de la isoterma.
    t_ref_index : int
        Índice temporal de referencia para calcular la anomalía.
    save_dir : str o Path
        Directorio donde guardar el NPZ de salida.
    nombre_base : str
        Nombre base del archivo (si None → se genera).

    Retorna
    -------
    outpath : Path
        Ruta del archivo .npz generado.
    """
    save_dir = Path(save_dir)
    save_dir.mkdir(exist_ok=True, parents=True)

    # ──────────── Leer NPZ bay / nobay ────────────
    d_b = np.load(npz_bay_path, allow_pickle=False)
    d_n = np.load(npz_nobay_path, allow_pickle=False)

    t_horas = d_b["t_horas"]     # (nt,)
    RC      = d_b["RC"]          # (nz,)
    Ab      = d_b["anom"]        # (nt,nz,nst) bay
    An      = d_n["anom"]        # (nt,nz,nst) nobay

    labels  = d_b["labels"]
    colors  = d_b["colors"]

    nt, nz, nst = Ab.shape

    # ──────────── Reconstruir T absoluta ────────────
    tref_raw = _tref_profile_for_scenario(scenario_label)
    tref = _match_tref_to_nz(tref_raw, nz)            # (nz,)
    tref_row = tref[None, :, None]                    # (1,nz,1)

    Tb = Ab + tref_row                                # bay
    Tn = An + tref_row                                # nobay

    # ──────────── Arreglos para profundidad ────────────
    depth_b = np.full((nt, nst), np.nan, dtype=np.float32)
    depth_n = np.full((nt, nst), np.nan, dtype=np.float32)

    # ──────────── Búsqueda de isoterma ────────────
    for s in range(nst):
        for it in range(nt):
            depth_b[it, s] = _isoterma_depth_1col(Tb[it, :, s], RC, T_iso)
            depth_n[it, s] = _isoterma_depth_1col(Tn[it, :, s], RC, T_iso)

    # ──────────── Desplazamientos ────────────
    ref_b = depth_b[t_ref_index, :]       # referencia por estación
    ref_n = depth_n[t_ref_index, :]

    eta_b = depth_b - ref_b[None, :]      # bay − bay_ref
    eta_n = depth_n - ref_n[None, :]      # nobay − nobay_ref

    delta_iso = depth_b - depth_n         # bay − nobay

    # ──────────── Guardar ────────────
    if nombre_base is None:
        nombre_base = f"{scenario_label}_iso_{T_iso:.2f}C"

    outpath = save_dir / f"{nombre_base}.npz"

    meta = {
        "scenario": scenario_label,
        "T_iso": float(T_iso),
        "t_ref_index": int(t_ref_index),
        "descripcion": "Profundidad (bay,nobay), desplazamiento relativo y diferencia bay−nobay de la isoterma.",
    }

    np.savez_compressed(
        outpath,
        t_horas=t_horas,
        RC=RC,
        labels=labels,
        colors=colors,
        depth_bay=depth_b,
        depth_nobay=depth_n,
        eta_bay=eta_b,
        eta_nobay=eta_n,
        delta_iso=delta_iso,
        meta_json=json.dumps(meta),
    )

    return outpath

# ─────────────────────────────────────────────────────────────
# PLOTEO: evolución temporal de la isoterma en una estación
# ─────────────────────────────────────────────────────────────
def plot_isoterma_en_estacion(
    npz_iso_path,
    estacion_label: str = "Gamma",
    unidad_tiempo: str = "horas",
    figsize=(8, 6),
    titulo_sufijo: str | None = None,
):
    """
    Grafica la evolución temporal de la profundidad de una isoterma
    en una estación dada, para bahía / sin bahía, y la diferencia bay−nobay.

    Parámetros
    ----------
    npz_iso_path : str o Path
        Ruta al archivo .npz generado por calcular_desplazamiento_isoterma_escenario.
    estacion_label : str
        Nombre de la estación a graficar (ej. "Gamma", "Omega", "A"...).
    unidad_tiempo : {"horas","dias"}
        Eje de tiempo en horas o días.
    figsize : tuple
        Tamaño de la figura matplotlib.
    titulo_sufijo : str o None
        Texto adicional para el título (ej. "2layer, T=7.1°C").
    """
    npz_iso_path = Path(npz_iso_path)
    data = np.load(npz_iso_path, allow_pickle=False)

    t_horas   = data["t_horas"]        # (nt,)
    labels    = data["labels"]         # (nst,)
    colors    = data["colors"]         # (nst,)
    depth_bay = data["depth_bay"]      # (nt, nst)
    depth_nob = data["depth_nobay"]    # (nt, nst)
    delta_iso = data["delta_iso"]      # (nt, nst)
    meta_json = data.get("meta_json", None)

    # Buscar índice de la estación
    labels_str = [str(l) for l in labels]
    if estacion_label not in labels_str:
        raise ValueError(
            f"Estación '{estacion_label}' no encontrada. Estaciones disponibles: {labels_str}"
        )
    s = labels_str.index(estacion_label)

    color_est = colors[s] if colors is not None and len(colors) > s else "C0"

    # Eje de tiempo
    if unidad_tiempo.lower().startswith("dia"):
        t = t_horas / 24.0
        xlabel = "Tiempo [días]"
    else:
        t = t_horas
        xlabel = "Tiempo [horas]"

    # Extraer series 1D
    z_b = depth_bay[:, s]   # profundidad positiva hacia abajo [m]
    z_n = depth_nob[:, s]
    dz  = delta_iso[:, s]   # bay − nobay

    # Info extra para el título
    escenario_txt = ""
    T_iso_txt = ""
    if meta_json is not None:
        meta = json.loads(str(meta_json))
        escenario_txt = meta.get("scenario", "")
        T_iso_val = meta.get("T_iso", None)
        if T_iso_val is not None:
            T_iso_txt = f"T = {T_iso_val:.2f}°C"

    titulo_base = f"Isoterma en estación {estacion_label}"
    partes = [titulo_base]
    if escenario_txt:
        partes.append(f"Escenario: {escenario_txt}")
    if T_iso_txt:
        partes.append(T_iso_txt)
    if titulo_sufijo:
        partes.append(titulo_sufijo)
    titulo = " | ".join(partes)

    # ──────────── Figura con dos paneles ────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=figsize)

    # Panel 1: profundidad de la isoterma (bay vs sin bahía)
    ax1.plot(t, z_b, label="Bahía", linewidth=2, color=color_est)
    ax1.plot(t, z_n, label="Sin bahía", linewidth=2, linestyle="--", color="k")

    ax1.set_ylabel("Profundidad isoterma [m]")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="best")
    # Como la profundidad es positiva hacia abajo, el eje crece hacia abajo naturalmente.
    # Si quisieras que el 0 esté arriba, podrías invertirlo:
    # ax1.invert_yaxis()

    # Panel 2: diferencia bay − nobay
    ax2.plot(t, dz, label="Bahía − Sin bahía", linewidth=2, color=color_est)
    ax2.axhline(0.0, color="k", linewidth=1, linestyle=":")

    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("Δ profundidad [m]")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best")

    fig.suptitle(titulo, fontsize=11)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    return fig, (ax1, ax2)


# =====================================================================
# Cálculo y guardado de rho (EOS lineal) en NetCDF para varios escenarios
# =====================================================================
_TREF_LINEAR = 33.6
_SREF_LINEAR = None

_TREF_REALISTIC = 16.582
_SREF_REALISTIC = 33.6

_TREF_2LAYER = 11.4
_SREF_2LAYER = None 


def _pick_vertical_1d(gr):
    """
    Elige coordenada vertical 1D entre ('RC','rC','rc','Z','z'), convierte a 1D.
    """
    for name in ("RC", "rC", "rc", "Z", "z"):
        if name in gr.variables:
            arr = np.array(gr.variables[name][:], dtype=float).squeeze()
            if arr.ndim > 1:
                arr = arr[..., 0, 0]
            return arr
    raise KeyError("No encontré coordenadas verticales ('RC','rC','Z').")


def _pick_2d_from_grid(gr, candidates):
    """
    Devuelve una matriz 2D de la primera variable disponible en 'candidates'.
    Si ninguna coincide, usa la primera variable 2D que encuentre.
    """
    for name in candidates:
        if name in gr.variables:
            arr = np.array(gr.variables[name][:], dtype=float)
            if arr.ndim > 2:
                arr = arr[..., 0, 0]
            return arr
    for name, v in gr.variables.items():
        if hasattr(v, "shape") and len(v.shape) == 2:
            return np.array(v[:], dtype=float)
    raise KeyError(f"No encontré variables {candidates} ni ninguna 2D en el grid.")


def _find_sal_var(nc):
    for k in ("Salt", "S", "salt"):
        if k in nc.variables:
            return k
    return None


def _compute_rho(T, S, tRef4d, sRef4d, rhoNil, tAlpha, sBeta):
    """
    tRef4d / sRef4d pueden ser escalares o arrays 4D compatibles con T/S.
    """
    rho = rhoNil * (1.0 - tAlpha * (T - tRef4d))
    if S is not None and sBeta != 0.0 and sRef4d is not None:
        rho = rho + rhoNil * sBeta * (S - sRef4d)
    return rho


def _fit_ref_to_nz(ref_vec: np.ndarray, nz: int, name: str = "tRef"):
    """
    Ajusta (padding o truncado) un vector de referencia a longitud nz.
    Si falta 1 nivel (caso común 49 vs 50), repite el último valor.
    """
    ref = np.array(ref_vec, dtype=float).ravel()
    if ref.size == nz:
        return ref
    if ref.size < nz:
        # Rellenar con el último valor para coincidir nz
        ref_pad = np.concatenate([ref, np.full(nz - ref.size, ref[-1] if ref.size else 0.0)])
        print(f"[AVISO] {name} tenía {ref.size} niveles, se rellenó hasta nz={nz} usando el último valor.")
        return ref_pad
    # ref.size > nz
    print(f"[AVISO] {name} tenía {ref.size} niveles, se truncó a nz={nz}.")
    return ref[:nz]


def guardar_rho_lineal_netcdf(
    root: Path,
    state_pattern: str | tuple | list,
    grid_pattern: str,
    out_path: Path,
    j_slice: slice,
    i_slice: slice,
    rhoNil: float,
    tAlpha: float,
    sBeta: float,
    tRef_vec: np.ndarray | float,
    sRef_vec: np.ndarray | float | None = None,
    compress_level: int = 4,
    temp_clip: float = 1e6,
    salt_clip: float = 1e6,
    apply_hfac_mask: bool = True
):
    """
    Calcula rho(nt,nz,ny,nx) = rhoNil * [1 - tAlpha*(T - tRef) + sBeta*(S - sRef)]
    y lo guarda en NetCDF con dims (T,Z,Y,X) sobre el subdominio indicado.
    """
    # Permite tRef/sRef escalares o vectores
    tRef_is_scalar = np.isscalar(tRef_vec)
    sRef_is_scalar = np.isscalar(sRef_vec) if sRef_vec is not None else False
    if tRef_is_scalar:
        tRef_base = float(tRef_vec)
    else:
        tRef_base = np.array(tRef_vec, dtype=float)
    # Si sRef_vec es None, usar 0 (sin efecto salino) pero evita None.
    if sRef_vec is None:
        sRef_base = 0.0 if tRef_is_scalar else np.zeros_like(tRef_base, dtype=float)
    else:
        sRef_base = float(sRef_vec) if sRef_is_scalar else np.array(sRef_vec, dtype=float)
    tRef = None
    sRef = None
    tRef4d = None  # se definirán tras conocer nz real
    sRef4d = None

    # grid
    gr = _open_mnc_or_raise(root, grid_pattern)
    RC = _pick_vertical_1d(gr)
    Yfull = _pick_2d_from_grid(gr, ("YC", "yC", "Y", "YG"))
    Xfull = _pick_2d_from_grid(gr, ("XC", "xC", "X", "XG"))
    hfac_name = _pick_var(gr, ("HFacC", "hFacC"))
    hFacC = np.asarray(gr.variables[hfac_name][:], dtype=float)
    YC = Yfull[j_slice, i_slice]
    XC = Xfull[j_slice, i_slice]
    hFacC = hFacC[:, j_slice, i_slice]
    ny, nx = YC.shape
    _close_mnc(gr)

    # datos
    times, rhos = [], []
    patterns = [state_pattern] if isinstance(state_pattern, (str, Path)) else list(state_pattern)
    for patt_idx, patt in enumerate(patterns):
        nc = _open_mnc_or_raise(root, patt)
        T = np.array(nc.variables["Temp"][:, :, j_slice, i_slice], dtype=float)
        S = None
        sal_var = _find_sal_var(nc)
        if sal_var and sBeta != 0.0:
            S = np.array(nc.variables[sal_var][:, :, j_slice, i_slice], dtype=float)
        # limpiar valores corruptos/extremos
        if temp_clip is not None:
            T = np.where(np.abs(T) <= temp_clip, T, np.nan)
        if S is not None and salt_clip is not None:
            S = np.where(np.abs(S) <= salt_clip, S, np.nan)
        t = _get_time_from_nc(nc)
        nz_here = T.shape[1]

        # Definir referencias ajustadas una sola vez usando nz del primer bloque
        if tRef4d is None:
            if tRef_is_scalar:
                tRef = tRef_base
                tRef4d = tRef  # escalar
            else:
                tRef = _fit_ref_to_nz(tRef_base, nz_here, name="tRef")
                tRef4d = tRef.reshape(1, -1, 1, 1)
            if sRef_is_scalar or (sRef_vec is None and tRef_is_scalar):
                sRef = sRef_base
                sRef4d = sRef  # escalar
            else:
                sRef = _fit_ref_to_nz(sRef_base, nz_here, name="sRef")
                sRef4d = sRef.reshape(1, -1, 1, 1)

        rho = _compute_rho(T, S, tRef4d, sRef4d, rhoNil, tAlpha, sBeta).astype(np.float32)
        if apply_hfac_mask:
            rho = np.where(hFacC[None, ...] > 0.0, rho, np.nan)
        times.append(t)
        rhos.append(rho)
        _close_mnc(nc)

    time_all, rho_all = _concat_time_and_data(times, rhos)  # (nt,nz,ny,nx)

    # guardar
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(out_path, "w") as ds:
        if np.isscalar(tRef4d):
            nt = time_all.shape[0]
            nz = rho_all.shape[1]
        else:
            nt, nz = time_all.shape[0], tRef4d.shape[1]
        ds.createDimension("T", nt); ds.createDimension("Z", nz)
        ds.createDimension("Y", ny); ds.createDimension("X", nx)

        vtime = ds.createVariable("time", "f8", ("T",)); vtime.units = "s"; vtime[:] = time_all
        vZ = ds.createVariable("Z", "f4", ("Z",)); vZ.units = "m"; vZ[:] = RC
        vYC = ds.createVariable("YC", "f4", ("Y", "X")); vXC = ds.createVariable("XC", "f4", ("Y", "X"))
        vYC.units = "deg"; vXC.units = "deg"; vYC[:] = YC; vXC[:] = XC

        vrho = ds.createVariable("rho", "f4", ("T", "Z", "Y", "X"), zlib=True, complevel=compress_level)
        vrho.units = "kg m-3"; vrho.long_name = "density_linear_eos"; vrho[:] = rho_all

        ds.rhoNil = rhoNil; ds.tAlpha = tAlpha; ds.sBeta = sBeta
        if np.isscalar(tRef):
            ds.tRef = float(tRef)
        else:
            ds.tRef = tRef.tolist() if tRef is not None else tRef_base.tolist()
        if np.isscalar(sRef):
            ds.sRef = float(sRef)
        else:
            ds.sRef = sRef.tolist() if sRef is not None else sRef_base.tolist()
        ds.eos = "rho = rhoNil * (1 - tAlpha*(Temp - tRef) + sBeta*(Salt - sRef))"
        ds.note = f"Subdominio j={j_slice}, i={i_slice}"

    print(f"[OK] Guardado {out_path} con rho shape {rho_all.shape}")


def guardar_rho_para_escenarios(
    scenarios_cfg: dict,
    bay_dir: str = "run_expand/mnc_00*",
    nobay_dir: str = "run_expand_nobay/mnc_00*",
    state_patterns: str | tuple | list = "state.*.nc",
    grid_glob: str = "grid.t*.nc",
    j_slice: slice = slice(None),
    i_slice: slice = slice(None),
    out_root: Path = Path("./rho_nc")
):
    """
    Recorre escenarios y genera NetCDF para bahía y sin bahía.
    scenarios_cfg: dict con claves de escenario y campos:
      base (Path), tRef (array), sRef (array|None), tAlpha, sBeta, rhoNil
    """
    for name, cfg in scenarios_cfg.items():
        for side, subdir in (("bay", bay_dir), ("nobay", nobay_dir)):
            if isinstance(state_patterns, (str, Path)):
                pattern_state = f"{subdir}/{state_patterns}"
            else:
                pattern_state = [f"{subdir}/{p}" for p in list(state_patterns)]
            pattern_grid = f"{subdir}/{grid_glob}"
            out_path = out_root / name / f"{name}_rho_{side}.nc"
            guardar_rho_lineal_netcdf(
                root=cfg["base"],
                state_pattern=pattern_state,
                grid_pattern=pattern_grid,
                out_path=out_path,
                j_slice=j_slice,
                i_slice=i_slice,
                rhoNil=cfg.get("rhoNil", 999.8),
                tAlpha=cfg.get("tAlpha", 2e-4),
                sBeta=cfg.get("sBeta", 0.0),
                tRef_vec=cfg["tRef"],
                sRef_vec=cfg.get("sRef")
            )


# ═══════════════════════════════════════════════════════════════
# ENERGÍA: PE barotropa + APE baroclinica (EOS lineal)
# ═══════════════════════════════════════════════════════════════

def barotropic_pe_series(eta, area, rho0=1035.0, g=9.81, mask=None):
    """
    PE barotropa: 0.5 * rho0 * g * ∫∫ η² dA

    Parámetros
    ----------
    eta  : (nt, ny, nx) superficie libre [m]
    area : (ny, nx) área de celda [m²]
    rho0 : densidad de referencia [kg/m³]
    g    : gravedad [m/s²]
    mask : (ny, nx) bool — True en celdas a excluir (tierra)

    Retorna
    -------
    pe_bt : (nt,) PE barotropa integrada [J]
    """
    eta  = np.asarray(eta, dtype=float)
    area = np.asarray(area, dtype=float)
    eta2 = eta ** 2
    if mask is not None:
        eta2 = np.where(np.asarray(mask, dtype=bool)[None], 0.0, eta2)
    else:
        eta2 = np.nan_to_num(eta2, nan=0.0)
    pe_bt = np.einsum("tyx,yx->t", eta2, area, optimize=True)
    return 0.5 * float(rho0) * float(g) * pe_bt


def baroclinic_ape_series(rho_anom, z_centers, drF, area, hFacC=None, g=9.81):
    """
    APE baroclinica: g * ∫∫∫ ρ_anom · z dV

    Parámetros
    ----------
    rho_anom  : (nt, nz, ny, nx) anomalía de densidad [kg/m³]
               (= ρ − ρ_ref, con la EOS lineal)
    z_centers : (nz,) coordenadas verticales de centro [m, negativas]
    drF       : (nz,) espesores de celda [m]
    area      : (ny, nx) área de celda [m²]
    hFacC     : (nz, ny, nx) fracción abierta; si None asume 1
    g         : gravedad [m/s²]

    Retorna
    -------
    ape : (nt,) APE integrada [J]
    """
    rho_anom  = np.asarray(rho_anom, dtype=float)
    z_centers = np.asarray(z_centers, dtype=float)
    drF       = np.asarray(drF, dtype=float)
    area      = np.asarray(area, dtype=float)

    vol_3d = drF[:, None, None] * area[None, :, :]   # (nz, ny, nx)
    if hFacC is not None:
        vol_3d = vol_3d * np.asarray(hFacC, dtype=float)

    zvol = z_centers[:, None, None] * vol_3d          # (nz, ny, nx)
    rho0 = np.nan_to_num(rho_anom, nan=0.0)
    ape  = np.einsum("tzyx,zyx->t", rho0, zvol, optimize=True)
    return float(g) * ape


def N2_from_rho_ref(rho_ref, z_centers, rho0=1035.0, g=9.81, n2_floor=1e-9):
    """N²(z) = −(g/ρ₀)·dρ̄/dz a partir del perfil de referencia ρ̄(z).

    Misma construcción que energy_flux_cg.N2_from_ref(): se recorta por abajo a
    n2_floor para no dividir por cero donde la columna sale homogénea (agua bien
    mezclada, o niveles enteramente en tierra).
    """
    rho_ref   = np.asarray(rho_ref, dtype=float)
    z_centers = np.asarray(z_centers, dtype=float)
    if rho_ref.size < 2:
        return np.full_like(rho_ref, n2_floor)
    drho_dz = np.gradient(rho_ref, z_centers)
    return np.clip(-float(g) / float(rho0) * drho_dz, n2_floor, None)


def available_ape_series(rho, rho_ref, N2, drF, area, hFacC=None,
                         rho0=1035.0, g=9.81):
    """APE *disponible*: ∫∫∫ ½·g²·ρ'² / (ρ₀·N²(z)) dV     [J]

    con ρ'(x,y,z,t) = ρ(x,y,z,t) − ρ̄(z).

    Es la forma estándar (Kang & Fringer 2010; Winters et al. 1995) y la que ya
    usa energy_flux_cg.py. Se diferencia de baroclinic_ape_series() —que calcula
    g·∫ρ_anom·z dV contra una referencia CONSTANTE— en que aquélla mide la energía
    potencial TOTAL respecto a un océano homogéneo: queda dominada por la
    estratificación de fondo (estática) y sale 10⁴–10⁷ veces mayor que KE, con la
    señal de onda como residuo diminuto sobre un número enorme.

    Parámetros
    ----------
    rho     : (nt, nz, ny, nx) densidad [kg/m³]
    rho_ref : (nz,)            perfil de referencia ρ̄(z) [kg/m³]
    N2      : (nz,)            frecuencia de Brunt-Väisälä al cuadrado [s⁻²]
    drF     : (nz,)            espesores de celda [m]
    area    : (ny, nx)         área de celda [m²]
    hFacC   : (nz, ny, nx)     fracción abierta; si None asume 1
    rho0    : densidad de referencia [kg/m³]
    g       : gravedad [m/s²]

    Retorna
    -------
    ape : (nt,) APE disponible integrada [J]
    """
    rho     = np.asarray(rho, dtype=float)
    rho_ref = np.asarray(rho_ref, dtype=float)
    N2      = np.asarray(N2, dtype=float)
    drF     = np.asarray(drF, dtype=float)
    area    = np.asarray(area, dtype=float)

    vol_3d = drF[:, None, None] * area[None, :, :]      # (nz, ny, nx)
    if hFacC is not None:
        vol_3d = vol_3d * np.asarray(hFacC, dtype=float)

    # peso = ½ g² vol / (ρ₀ N²)   →  APE = Σ peso · ρ'²
    peso = 0.5 * float(g) ** 2 * vol_3d / (float(rho0) * N2[:, None, None])

    rho_p = np.nan_to_num(rho - rho_ref[None, :, None, None], nan=0.0)
    return np.einsum("tzyx,zyx->t", rho_p ** 2, peso, optimize=True)


def _rho_ref_accum_from_mnc_state(nc_state, gr,
                                  nz_sel, y_slice, x_slice,
                                  tRef_vec, sRef_vec,
                                  rhoNil, tAlpha, sBeta,
                                  ny_block=64):
    """Acumuladores para ρ̄(z): devuelve (num, den) por nivel.

    ρ̄(z) = Σ_t Σ_{x,y} ρ·A·hFac  /  Σ_t Σ_{x,y} A·hFac  → media horizontal
    ponderada por volumen y promediada en el tiempo. Se acumula por archivo para
    poder promediar sobre TODAS las etapas antes de calcular ninguna APE
    (mismo criterio que calcular_medias_run en energy_flux_IWBay.py).
    """
    Z      = np.array(gr.variables["Z"][:], dtype=float)
    nz_s   = _normalize_slice(nz_sel, Z.shape[0])

    _rA_name    = _pick_var(gr, ("rAc",   "rA"))
    _hFacC_name = _pick_var(gr, ("hFacC", "HFacC"))
    area  = np.array(gr.variables[_rA_name][y_slice, x_slice], dtype=float)
    hFacC = np.array(gr.variables[_hFacC_name][nz_s, y_slice, x_slice], dtype=float)

    nz_n = hFacC.shape[0]
    num  = np.zeros(nz_n, dtype=np.float64)
    den  = np.zeros(nz_n, dtype=np.float64)

    tRef_is_scalar = np.ndim(tRef_vec) == 0
    sRef_is_scalar = (sRef_vec is None) or (np.ndim(sRef_vec) == 0)

    ny = area.shape[0]
    for j0 in range(0, ny, ny_block):
        j1    = min(ny, j0 + ny_block)
        j_gl0 = (y_slice.start or 0) + j0
        j_blk = slice(j_gl0, (y_slice.start or 0) + j1)

        T_blk = np.array(nc_state.variables["Temp"][:, nz_s, j_blk, x_slice], dtype=float)
        S_blk = (np.array(nc_state.variables["S"][:, nz_s, j_blk, x_slice], dtype=float)
                 if sBeta != 0.0 else None)

        nz_b = T_blk.shape[1]
        tR4d = (float(tRef_vec) if tRef_is_scalar
                else _fit_ref_to_nz(tRef_vec, nz_b, "tRef").reshape(1, -1, 1, 1))
        if sRef_is_scalar or sRef_vec is None:
            sR4d = 0.0 if sRef_vec is None else float(sRef_vec)
        else:
            sR4d = _fit_ref_to_nz(sRef_vec, nz_b, "sRef").reshape(1, -1, 1, 1)

        rho_blk = _compute_rho(T_blk, S_blk, tR4d, sR4d, rhoNil, tAlpha, sBeta)
        rho_blk = np.nan_to_num(rho_blk, nan=0.0)

        w = area[j0:j1, :][None, :, :] * hFacC[:, j0:j1, :]      # (nz, ny_b, nx)
        num += np.einsum("tzyx,zyx->z", rho_blk, w, optimize=True)
        den += float(T_blk.shape[0]) * w.sum(axis=(1, 2))

    return num, den


def _pe_ape_from_mnc_state(nc_state, gr,
                            nz_sel, y_slice, x_slice,
                            tRef_vec, sRef_vec,
                            rhoNil, tAlpha, sBeta,
                            rho0=1035.0, g=9.81,
                            ny_block=64,
                            ape_method="available", rho_ref=None):
    """
    Calcula PE_BT(t) y APE(t) leyendo de archivos MNC.

    La densidad se computa desde la EOS lineal usando Temp (y Salt si sBeta≠0)
    del archivo de estado, para ser consistente con el modelo.

    ape_method : 'available' (por defecto) → ½g²ρ'²/(ρ₀N²) con ρ' contra ρ̄(z);
                 requiere pasar rho_ref (nz,) calculado sobre TODAS las etapas.
                 'total' → el cálculo antiguo g·∫ρ_anom·z dV contra ρ_nil constante.

    Retorna
    -------
    t      : (nt,) tiempo en segundos
    pe_bt  : (nt,) PE barotropa [J]
    ape    : (nt,) APE baroclinica [J]
    """
    if ape_method not in ("available", "total"):
        raise ValueError(f"ape_method debe ser 'available' o 'total', no {ape_method!r}")
    if ape_method == "available" and rho_ref is None:
        raise ValueError("ape_method='available' requiere rho_ref (nz,)")
    # ── coordenadas y grilla ──────────────────────────────────
    Z   = np.array(gr.variables["Z"][:],   dtype=float)       # (nz,)
    drF = np.array(gr.variables["drF"][:], dtype=float)       # (nz,)

    nz_all = Z.shape[0]
    nz_s   = _normalize_slice(nz_sel,  nz_all)

    Z_sel   = Z[nz_s]
    drF_sel = drF[nz_s]

    _rA_name    = _pick_var(gr, ("rAc",   "rA"))
    _hFacC_name = _pick_var(gr, ("hFacC", "HFacC"))
    area  = np.array(gr.variables[_rA_name][y_slice, x_slice], dtype=float)  # (ny, nx)
    hFacC = np.array(gr.variables[_hFacC_name][nz_s, y_slice, x_slice], dtype=float)

    # Máscara de tierra en la superficie para PE_BT
    hFac_sfc = np.array(gr.variables[_hFacC_name][0, y_slice, x_slice], dtype=float)
    mask_land = (hFac_sfc == 0.0)

    # ── tiempo ────────────────────────────────────────────────
    t = np.array(nc_state.variables["T"][:], dtype=float)
    nt = t.shape[0]

    # ── PE barotropa (no requiere bucle por bloque) ────────────
    eta = np.array(nc_state.variables["Eta"][:, y_slice, x_slice], dtype=float)
    pe_bt = barotropic_pe_series(eta, area, rho0=rho0, g=g, mask=mask_land)

    # ── APE baroclinica por bloques espaciales ─────────────────
    ny, nx = area.shape
    ape_acc = np.zeros(nt, dtype=np.float64)

    # N²(z) del perfil de referencia (solo para ape_method='available')
    N2 = (N2_from_rho_ref(rho_ref, Z_sel, rho0=rho0, g=g)
          if ape_method == "available" else None)

    # Referencias escalares o vectores por nivel
    tRef_is_scalar = np.ndim(tRef_vec) == 0
    sRef_is_scalar = (sRef_vec is None) or (np.ndim(sRef_vec) == 0)

    for j0 in range(0, ny, ny_block):
        j1      = min(ny, j0 + ny_block)
        j_gl0   = (y_slice.start or 0) + j0
        j_gl1   = (y_slice.start or 0) + j1
        j_blk   = slice(j_gl0, j_gl1)

        T_blk = np.array(nc_state.variables["Temp"][:, nz_s, j_blk, x_slice], dtype=float)
        if sBeta != 0.0:
            S_blk = np.array(nc_state.variables["S"][:, nz_s, j_blk, x_slice], dtype=float)
        else:
            S_blk = None

        # Armar referencias 4D compatibles con (nt, nz, ny_blk, nx)
        nz_b = T_blk.shape[1]
        if tRef_is_scalar:
            tR4d = float(tRef_vec)
        else:
            tR4d = _fit_ref_to_nz(tRef_vec, nz_b, "tRef").reshape(1, -1, 1, 1)

        if sRef_is_scalar or sRef_vec is None:
            sR4d = 0.0 if sRef_vec is None else float(sRef_vec)
        else:
            sR4d = _fit_ref_to_nz(sRef_vec, nz_b, "sRef").reshape(1, -1, 1, 1)

        rho_blk = _compute_rho(T_blk, S_blk, tR4d, sR4d, rhoNil, tAlpha, sBeta)

        area_blk  = area[j0:j1, :]
        hFacC_blk = hFacC[: , j0:j1, :]

        if ape_method == "available":
            ape_acc += available_ape_series(rho_blk, rho_ref, N2, drF_sel,
                                            area_blk, hFacC=hFacC_blk,
                                            rho0=rho0, g=g)
        else:
            ape_acc += baroclinic_ape_series(rho_blk - rhoNil, Z_sel, drF_sel,
                                             area_blk, hFacC=hFacC_blk, g=g)

    return t, pe_bt, ape_acc


def guardar_energy_series_npz_etapas(
    scenario_root,
    scenario_label: str,
    tRef_vec,
    sRef_vec=None,
    tAlpha: float = 2e-4,
    sBeta: float = 0.0,
    rhoNil: float = 1035.0,
    subdir_bay: str = "run_expand",
    subdir_nobay: str = "run_expand_nobay",
    mnc_glob: str = "mnc_00*",
    stage_patterns=("state.0000000000.t*.nc", "state.0000007200.t*.nc"),
    grid_pattern: str = "grid.t*.nc",
    nz_sel: slice = slice(None),
    y_slice: slice = slice(None),
    x_slice: slice = slice(None),
    ny_block: int = 64,
    rho0: float = 1035.0,
    g: float = 9.81,
    save_dir=None,
    out_name: str | None = None,
    tol_time: float = 1e-7,
    ke_npz_path=None,
    ape_method: str = "available",
    verbose: bool = True,
):
    """
    Calcula y guarda series de tiempo de KE, PE_BT y APE para bay y nobay.

    La densidad para APE, PE_BT y (opcionalmente) KE se computa mediante la EOS lineal:
      ρ = rhoNil * (1 - tAlpha*(T - tRef) + sBeta*(S - sRef))

    Parámetros
    ----------
    ke_npz_path : Path | str | None
        Si se proporciona, KE se carga desde este NPZ pre-calculado (que debe contener
        'ke_bay', 'ke_nobay', 'time_bay', 'time_nobay') en lugar de calcularse con rho0
        constante. Esto permite usar KE = 0.5*∫ρ(u²+v²+w²)dV con densidad completa
        ρ = ρ_ref + ρ_anom.  Si es None, KE se calcula con rho0 constante (Boussinesq).
    ape_method : {'available', 'total'}
        'available' (por defecto): APE = ∫½g²ρ'²/(ρ₀N²)dV con ρ' medido contra el
        perfil ρ̄(z) promediado en el tiempo sobre TODAS las etapas y en la horizontal
        sobre la región. Requiere una primera pasada de lectura para obtener ρ̄(z).
        'total': el cálculo antiguo g·∫ρ_anom·z dV contra ρ_nil constante, que mide
        energía potencial total (dominada por la estratificación de fondo).

    Archivo NPZ de salida contiene:
      t_seconds, t_horas,
      ke_bay, ke_nobay,
      pe_bt_bay, pe_bt_nobay,
      ape_bay,   ape_nobay,
      meta_json
    """
    scenario_root = Path(scenario_root)
    if save_dir is None:
        save_dir = Path(__file__).resolve().parent
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    def _open_first(patterns):
        for patt in patterns:
            full = str(scenario_root / patt)
            if glob.glob(full):
                return MIT.mnc_files(full)
        raise FileNotFoundError(
            f"No encontré archivos MNC para {scenario_label}.\n"
            f"Patrones intentados:\n  " + "\n  ".join(patterns) +
            f"\nEn: {scenario_root}"
        )

    gr = _open_first([
        f"{subdir_bay}/{mnc_glob}/{grid_pattern}",
        f"{subdir_bay}/**/{grid_pattern}",
    ])
    # Grilla PROPIA de la corrida sin bahía: no comparten batimetría. Antes se
    # usaba la de run_expand para ambas, así que el hFacC de la bahía marcaba
    # como agua unas celdas que en run_expand_nobay son tierra rellena con
    # Temp=0 → ρ(T=0)=1002.08 kg/m³ entraba en la APE como si fuera océano.
    # (KE y PE_BT no se veían afectadas porque U,V,W,Eta son 0 en tierra.)
    gr_nb = _open_first([
        f"{subdir_nobay}/{mnc_glob}/{grid_pattern}",
        f"{subdir_nobay}/**/{grid_pattern}",
    ])

    # Normalizar slices con dimensiones reales del grid
    _rA_name = _pick_var(gr, ("rAc", "rA"))
    nz_all = gr.variables["Z"][:].shape[0]
    ny_all = gr.variables[_rA_name][:].shape[0]
    nx_all = gr.variables[_rA_name][:].shape[1]
    nz_s = _normalize_slice(nz_sel, nz_all)
    y_s  = _normalize_slice(y_slice, ny_all)
    x_s  = _normalize_slice(x_slice, nx_all)

    # Centros de celda de los niveles seleccionados: se capturan aquí porque más
    # abajo se cierra `gr` y se necesitan para guardar N²(z) en el npz.
    gr_Z_sel = np.array(gr.variables["Z"][:], dtype=float)[nz_s]

    # ── PASADA 1 (solo para ape_method='available'): ρ̄(z) sobre TODAS las etapas ──
    # Tiene que ser sobre el registro completo antes de calcular ninguna APE; si se
    # hiciera por etapa, cada etapa mediría su anomalía contra un fondo distinto.
    rho_ref_bay = rho_ref_nobay = None
    if ape_method == "available":
        num_b = den_b = num_n = den_n = None
        try:
            for patt in stage_patterns:
                nc_b = _open_first([f"{subdir_bay}/{mnc_glob}/{patt}",
                                     f"{subdir_bay}/**/{patt}"])
                nc_n = _open_first([f"{subdir_nobay}/{mnc_glob}/{patt}",
                                     f"{subdir_nobay}/**/{patt}"])
                try:
                    for nc, gr_i, which in ((nc_b, gr, "bay"), (nc_n, gr_nb, "nobay")):
                        nu, de = _rho_ref_accum_from_mnc_state(
                            nc, gr_i, nz_sel=nz_s, y_slice=y_s, x_slice=x_s,
                            tRef_vec=tRef_vec, sRef_vec=sRef_vec,
                            rhoNil=rhoNil, tAlpha=tAlpha, sBeta=sBeta,
                            ny_block=ny_block,
                        )
                        if which == "bay":
                            num_b = nu if num_b is None else num_b + nu
                            den_b = de if den_b is None else den_b + de
                        else:
                            num_n = nu if num_n is None else num_n + nu
                            den_n = de if den_n is None else den_n + de
                finally:
                    nc_b.close()
                    nc_n.close()
        finally:
            pass

        def _perfil(num, den):
            # den=0 → nivel enteramente en tierra: ρ̄=rhoNil evita 0/0; el volumen
            # es cero de todos modos, así que no aporta a la APE.
            out = np.full_like(den, float(rhoNil), dtype=float)
            ok = den > 0
            out[ok] = num[ok] / den[ok]
            return out

        rho_ref_bay   = _perfil(num_b, den_b)
        rho_ref_nobay = _perfil(num_n, den_n)
        if verbose:
            for lbl, rr, dd in (("bay", rho_ref_bay, den_b),
                                ("nobay", rho_ref_nobay, den_n)):
                nwet = int((dd > 0).sum())
                if nwet == 0:
                    print(f"    ρ̄(z) {lbl:<5}: región SIN celdas mojadas → APE ≡ 0")
                else:
                    print(f"    ρ̄(z) {lbl:<5}: {nwet}/{len(dd)} niveles con agua | "
                          f"sup={rr[0]:.4f}  fondo={rr[nwet-1]:.4f}  "
                          f"Δρ={rr[nwet-1]-rr[0]:+.5f} kg/m³")

    # Acumular por etapa
    tb_list, tnb_list = [], []
    ke_b_list, ke_nb_list = [], []   # solo se usan si ke_npz_path is None
    pe_b_list, pe_nb_list = [], []
    ape_b_list, ape_nb_list = [], []

    try:
        for patt in stage_patterns:
            nc_b = _open_first([f"{subdir_bay}/{mnc_glob}/{patt}",
                                 f"{subdir_bay}/**/{patt}"])
            nc_n = _open_first([f"{subdir_nobay}/{mnc_glob}/{patt}",
                                 f"{subdir_nobay}/**/{patt}"])

            try:
                for nc, gr_i, tb_l, ke_l, pe_l, ape_l, rr in (
                    (nc_b, gr,    tb_list,  ke_b_list,  pe_b_list,  ape_b_list,  rho_ref_bay),
                    (nc_n, gr_nb, tnb_list, ke_nb_list, pe_nb_list, ape_nb_list, rho_ref_nobay),
                ):
                    t_stage, pe_bt_stage, ape_stage = _pe_ape_from_mnc_state(
                        nc, gr_i,
                        nz_sel=nz_s, y_slice=y_s, x_slice=x_s,
                        tRef_vec=tRef_vec, sRef_vec=sRef_vec,
                        rhoNil=rhoNil, tAlpha=tAlpha, sBeta=sBeta,
                        rho0=rho0, g=g, ny_block=ny_block,
                        ape_method=ape_method, rho_ref=rr,
                    )
                    tb_l.append(t_stage)
                    pe_l.append(pe_bt_stage)
                    ape_l.append(ape_stage)
                    # KE: solo calcular aquí si no se proporciona NPZ externo
                    if ke_npz_path is None:
                        _, ke_stage = ke_time_series_from_mnc_state(
                            nc, gr_i,
                            nz_sel=nz_s, y_slice=y_s, x_slice=x_s,
                            ny_block=ny_block, rho0=rho0,
                        )
                        ke_l.append(ke_stage)
            finally:
                nc_b.close()
                nc_n.close()
    finally:
        gr.close()
        gr_nb.close()

    t_bay,  pe_bay    = _concat_time_and_data(tb_list,   pe_b_list)
    t_nobay, pe_nobay = _concat_time_and_data(tnb_list,  pe_nb_list)
    _,       ape_bay  = _concat_time_and_data(tb_list,   ape_b_list)
    _,       ape_nobay= _concat_time_and_data(tnb_list,  ape_nb_list)

    # KE: desde NPZ pre-calculado (densidad completa) o desde EOS lineal con rho0
    if ke_npz_path is not None:
        _kd = np.load(ke_npz_path)
        # Dos productores escriben este npz con esquemas distintos:
        #   - test_KE.ipynb (feb-2026): 'time_bay' / 'time_nobay' (ejes separados)
        #   - guardar_ke_rho_diff_npz  : 't_seconds' (eje unico ya alineado, ver
        #     regenerar_ke_rho.py). Ahi bay y nobay comparten t_common.
        if "time_bay" in _kd.files:
            t_ke_b  = np.asarray(_kd["time_bay"],   dtype=np.float64)
            t_ke_nb = np.asarray(_kd["time_nobay"], dtype=np.float64)
        elif "t_seconds" in _kd.files:
            t_ke_b = t_ke_nb = np.asarray(_kd["t_seconds"], dtype=np.float64)
        else:
            raise KeyError(
                f"{Path(ke_npz_path).name}: no encontre eje temporal "
                f"('time_bay'/'time_nobay' ni 't_seconds'). Claves: {_kd.files}"
            )
        ke_bay   = np.interp(t_bay,   t_ke_b,  np.asarray(_kd["ke_bay"],   dtype=np.float64))
        ke_nobay = np.interp(t_nobay, t_ke_nb, np.asarray(_kd["ke_nobay"], dtype=np.float64))
        print(f"  [KE] Cargado desde {Path(ke_npz_path).name} (densidad completa ρ = ρ_ref + ρ_anom)")
    else:
        _, ke_bay   = _concat_time_and_data(tb_list,   ke_b_list)
        _, ke_nobay = _concat_time_and_data(tnb_list,  ke_nb_list)

    # Alinear tiempos bay/nobay para PE y APE
    if not (t_bay.shape == t_nobay.shape and
            np.all(np.abs(t_bay - t_nobay) <= tol_time)):
        print(f"[AVISO] Tiempos bay/nobay difieren; interpolando nobay a tiempos de bay.")
        ke_nobay  = np.interp(t_bay, t_nobay, ke_nobay)
        pe_nobay  = np.interp(t_bay, t_nobay, pe_nobay)
        ape_nobay = np.interp(t_bay, t_nobay, ape_nobay)
    t_common = t_bay

    if out_name is None:
        out_name = f"{scenario_label}_energy_series.npz"

    meta = {
        "scenario": scenario_label,
        "root": str(scenario_root),
        "subdir_bay": subdir_bay,
        "subdir_nobay": subdir_nobay,
        "stage_patterns": list(stage_patterns),
        "nz_sel": [nz_s.start, nz_s.stop, nz_s.step],
        "y_slice": [y_s.start,  y_s.stop,  y_s.step],
        "x_slice": [x_s.start,  x_s.stop,  x_s.step],
        "rho0": float(rho0), "rhoNil": float(rhoNil),
        "tAlpha": float(tAlpha), "sBeta": float(sBeta),
        "g": float(g),
        "units": {"t": "s", "ke": "J", "pe_bt": "J", "ape": "J"},
        "ke_source": (
            str(ke_npz_path) if ke_npz_path is not None
            else "rho0_constante"
        ),
        "ape_method": ape_method,
        # OJO: antes las dos ramas se concatenaban mal y la descripción de la APE
        # solo aparecía cuando ke_npz_path era None.
        "note": (
            ("KE = 0.5*∫ρ(u²+v²+w²)dV con ρ = ρ_ref + ρ_anom (de NPZ pre-calculado). "
             if ke_npz_path is not None else
             "KE usa rho0 constante (aproximación de Boussinesq). ")
            + ("APE = ∫½g²ρ'²/(ρ₀N²)dV, ρ' = ρ − ρ̄(z), ρ̄ promediada en el tiempo "
               "sobre todas las etapas y en la horizontal sobre la región (energía "
               "DISPONIBLE, comparable con KE)."
               if ape_method == "available" else
               "APE = g * sum(rho_anom * z * vol), rho_anom de EOS lineal contra "
               "ρ_nil constante (energía potencial TOTAL, no disponible).")
        ),
    }

    save_kw = {}
    if ape_method == "available":
        Z_all = np.array(gr_Z_sel, dtype=float)
        save_kw["rho_ref_bay"]   = np.asarray(rho_ref_bay,   dtype=np.float64)
        save_kw["rho_ref_nobay"] = np.asarray(rho_ref_nobay, dtype=np.float64)
        save_kw["N2_bay"]        = N2_from_rho_ref(rho_ref_bay,   Z_all, rho0=rho0, g=g)
        save_kw["N2_nobay"]      = N2_from_rho_ref(rho_ref_nobay, Z_all, rho0=rho0, g=g)
        save_kw["z_centers"]     = Z_all

    out_path = save_dir / out_name
    np.savez_compressed(
        out_path,
        t_seconds  = t_common.astype(np.float64),
        t_horas    = (t_common / 3600.0).astype(np.float64),
        ke_bay     = ke_bay.astype(np.float64),
        ke_nobay   = ke_nobay.astype(np.float64),
        pe_bt_bay  = pe_bay.astype(np.float64),
        pe_bt_nobay= pe_nobay.astype(np.float64),
        ape_bay    = ape_bay.astype(np.float64),
        ape_nobay  = ape_nobay.astype(np.float64),
        meta_json  = np.array(json.dumps(meta, ensure_ascii=False)),
        **save_kw,
    )
    print(f"[OK] Guardado: {out_path}")
    return out_path
