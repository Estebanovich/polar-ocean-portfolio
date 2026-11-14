#!/usr/bin/env python3
"""
energy_flux_IWBay.py
====================
Calcula flujos de energía de ondas internas baroclínicas para los tres
escenarios IWBay (dos capas, realista, lineal), con y sin bahía,
siguiendo el método de Nash, Alford & Kunze (2005).

  Nash, J.D., Alford, M.H., Kunze, E. (2005).
  "Estimating internal wave energy fluxes in the ocean."
  J. Atmos. Oceanic Technol., 22, 1551-1570.

─────────────────────────────────────────────────────
METODOLOGÍA (Nash et al. 2005, eqs. 2–5, 10)
─────────────────────────────────────────────────────
Las medias temporales ρ̄(z) y ū(z) se calculan sobre AMBAS etapas
(0–120 h) antes de procesar ninguna de ellas individualmente.

Paso 1  ρ̄(z)     = media temporal de ρ en la SECCIÓN DE LA BOCA, AMBAS etapas  (eq. 2)
Paso 2  ū(z,x)   = media temporal de u,v en la sección de la boca, AMBAS etapas (eq. 5)
Paso 3  ρ'(z,t)  = ρ(z,t) − ρ̄(z)                                               (eq. 2)
Paso 4  p̃(z,t)  = ρ₀·g·η + g·∫_z^0 ρ'(z') dz'                                 (eq. 3)
Paso 5  p'(z,t)  = p̃ − H⁻¹·∫p̃ dz          (remoción barotropa)               (eq. 4)
Paso 6  u'(z,t)  = [u − ū(z)] − H⁻¹·∫[u − ū] dz  (pert. barocl.)             (eq. 5)
Paso 7  F(t)     = ∫ p'(z,t)·u'(z,t) dz                                        (eq. 10)
Paso 8  ⟨F⟩     = media temporal simple de F(t) sobre la etapa

H es la profundidad LOCAL de la sección y las integrales verticales corren sólo
sobre sus niveles húmedos (en j=323 eso son los 599.9 m / 50 niveles completos;
niveles_humedos_seccion() lo verifica y aborta si la sección tocara tierra).

─────────────────────────────────────────────────────
Las etapas 1 (forzado, 0–60 h) y 2 (libre, 60–120 h) se descubren en disco:
la corrida de jun-2026 escribe 12 chunks de 1260 iteraciones, 6 por etapa, y
comparte el instante de frontera T=216000 s entre los chunks 6300 y 7200.

Salida: npz_flux/
  {escenario}_{run}_flujo.npz
  claves: t_seconds, t_horas, Fy_sec, Fx_sec,
          Fy_mean, Fx_mean, rho_ref, u_mean_sec, v_mean_sec, drF,
          dx_sec (dx real por columna — la malla es telescópica), meta_json

Variante: energy_flux_IWBay_rholocal.py hace lo mismo con una ρ̄(z,x) que
conserva variación horizontal, y escribe en {escenario}_{run}_flujo_rholocal.npz.
─────────────────────────────────────────────────────
Uso:
  /opt/anaconda3/bin/python3 energy_flux_IWBay.py
"""

import json, glob
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

# ──────────────────────────────────────────────────────────────────────
# Rutas base
# ──────────────────────────────────────────────────────────────────────
BASE_K = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification")
BASE_L = Path("/Volumes/Esteban_SSD/Maestría/MITgcm/verification")

SAVE_DIR = Path(__file__).resolve().parent / "npz_flux"

# ──────────────────────────────────────────────────────────────────────
# Configuración de escenarios
# ──────────────────────────────────────────────────────────────────────
SCENARIOS = {
    "2layer": dict(
        base   = BASE_K / "BayIW_Rectan_2layer",
        tRef   = 11.4,
        sRef   = 0.0,
        tAlpha = 2e-4,
        sBeta  = 0.0,
        rhoNil = 999.8,
    ),
    "realistic": dict(
        base   = BASE_K / "BayIW_Rectan",
        tRef   = 16.582,
        sRef   = 33.6,
        tAlpha = 2e-4,
        sBeta  = 7.4e-4,
        rhoNil = 999.8,
    ),
    "linear": dict(
        base   = BASE_L / "BayIW_Rectan_linear",
        tRef   = 33.6,
        sRef   = 0.0,
        tAlpha = 2e-4,
        sBeta  = 0.0,
        rhoNil = 999.8,
    ),
}

RUNS   = ["run_expand", "run_expand_nobay"]

MNC_GLOB = "mnc_00*"
G        = 9.81   # m/s²

# Frontera entre etapas: iteración 7200 = 216000 s = 60 h.
# La corrida de jun-2026 NO escribe 2 archivos monolíticos sino 12 chunks
# (cada 1260 iteraciones), 6 por etapa. Los patrones se descubren en disco con
# discover_stage_groups(); el antiguo STAGES fijo de 2 elementos leía sólo los
# chunks 0000000000 y 0000007200, es decir 2 de 12 → 42 de 241 pasos por etapa.
ITER_ETAPA2 = 7200

# Sección de la boca de la bahía (índices globales)
#
# Geometría real (idéntica en los 3 escenarios, verificada contra Depth de
# grid.t*.nc): el océano abierto tiene 599.9 m hasta j=323 inclusive; la bahía
# ocupa j=324..347 con 163.3 m de profundidad y sólo las columnas i=274..285.
# La BOCA es entonces la cara entre j=323 y j=324, y la sección se toma en la
# ÚLTIMA fila de océano, j=323, que es la que está pegada a esa cara.
#
# Hasta ago-2026 se usaba J_BOCA=322 (dos filas mar adentro) y X_BOCA=slice(273,285)
# (corrida una columna al oeste: metía océano abierto en i=273 y dejaba fuera i=285,
# justo el borde donde el flujo es máximo). Como frente a la costa hay un patrón de
# onda estacionaria, max|F| cambiaba hasta 3× entre filas contiguas y ⟨F⟩ cambiaba
# de signo; para el escenario 2layer j=322 caía además en un mínimo local, lo que
# hacía parecer que la bahía reducía el flujo baroclínico.
#
# Se usa j=323 y NO j=324 (ya dentro de la bahía) porque en el run sin bahía la
# ensenada está rellena de tierra: j=324 está seca ahí y el control desaparecería.
# En j=323 la geometría es idéntica en ambos runs (599.9 m, todo húmedo), así que
# la comparación bahía / sin-bahía es directa y la descomposición de Nash usa la
# columna local completa (∫p'dz = 0 sobre los 599.9 m reales).
J_BOCA   = 323              # última fila de océano, pegada a la cara 323/324
X_BOCA   = slice(274, 286)  # huella real de la bahía (12 columnas)

# Fallback si no se puede leer la grilla. OJO: desde la re-corrida de jun-2026 la
# malla es telescópica (dx va de 5000 m en el centro a 13266 m en los bordes), así
# que este escalar ya no describe el dominio. En la sección de la boca dx≈5052.6 m.
# El código usa dx_sec leído de dxF en la grilla; DX queda sólo como respaldo.
DX = 5000.0   # m


def discover_stage_groups(run_dir, mnc_glob=MNC_GLOB, iter_split=ITER_ETAPA2):
    """Descubre los chunks `state.<iter>.t*.nc` y los agrupa en las 2 etapas.

    Devuelve [globs_etapa1, globs_etapa2] ordenados por iteración, donde la
    etapa 1 (forzado) son los chunks con iter < iter_split y la etapa 2 (libre)
    los de iter >= iter_split. Funciona igual con salida monolítica (2 chunks)
    que troceada (12 chunks).
    """
    iters = set()
    for d in sorted(glob.glob(str(Path(run_dir) / mnc_glob))):
        for f in glob.glob(str(Path(d) / "state.*.t*.nc")):
            name = Path(f).name
            try:
                iters.add(int(name.split("state.")[1].split(".t")[0]))
            except (IndexError, ValueError):
                continue
    if not iters:
        raise FileNotFoundError(f"No hay archivos state.*.t*.nc en {run_dir}/{mnc_glob}/")
    g1 = [f"state.{i:010d}.t*.nc" for i in sorted(iters) if i <  iter_split]
    g2 = [f"state.{i:010d}.t*.nc" for i in sorted(iters) if i >= iter_split]
    return [g1, g2]


def _tiempos_de_glob(run_dir, stage_glob, mnc_glob=MNC_GLOB):
    """Vector T de un chunk (idéntico en todos sus tiles). None si no existe."""
    for d in sorted(glob.glob(str(Path(run_dir) / mnc_glob))):
        for f in sorted(glob.glob(str(Path(d) / stage_glob))):
            with Dataset(f) as ds:
                return ds.variables["T"][:].data.copy()
    return None


def _mascaras_sin_duplicados(run_dir, stage_globs, mnc_glob=MNC_GLOB):
    """Máscara por chunk de los registros de tiempo aún no vistos.

    Los chunks comparten el instante de frontera (p.ej. T=216000 s aparece al
    final del chunk 6300 y al inicio del 7200), así que concatenarlos a ciegas
    duplica ese paso. Devuelve una lista de máscaras booleanas alineada con
    stage_globs (None donde el chunk no existe).
    """
    vistos, masks = set(), []
    for g in stage_globs:
        T = _tiempos_de_glob(run_dir, g, mnc_glob=mnc_glob)
        if T is None:
            masks.append(None)
            continue
        keys = [round(float(t), 3) for t in T]
        masks.append(np.array([k not in vistos for k in keys], dtype=bool))
        vistos.update(keys)
    return masks


def niveles_humedos_seccion(run_dir, j_boca=J_BOCA, x_boca=X_BOCA):
    """Niveles húmedos de la sección: (nz_sec, drF_sec, profundidad).

    Se recorta a los niveles cuya cara inferior cabe en la columna MENOS profunda
    de la sección, y se aborta si alguna columna es tierra. Con J_BOCA=323 la
    sección está en océano abierto y esto devuelve los 50 niveles (599.9 m); la
    verificación existe para que mover la sección hacia la bahía (163.3 m, seca en
    el run sin bahía) no meta celdas secas (Temp=S=0 en MITgcm) en ρ' ni un H
    equivocado en la remoción barotropa de Nash (eqs. 4–6).
    """
    drF = _leer_drF(run_dir)
    tile_info, _, _ = build_tile_map(run_dir, MNC_GLOB, "grid.t*.nc")
    prof = []
    for tile in tiles_for_section(tile_info, j_boca, x_boca):
        with Dataset(tile["file"]) as dg:
            prof.append(dg.variables["Depth"][:].data[tile["j_local"],
                                                      tile["x_local_slice"]])
    if not prof:
        raise ValueError(f"Ningún tile cubre j={j_boca}, x={x_boca} en {run_dir}")
    prof = np.concatenate(prof)
    if np.any(prof <= 0):
        raise ValueError(f"La sección j={j_boca}, x={x_boca} toca tierra "
                         f"({int((prof <= 0).sum())} columnas secas) en {run_dir}")
    H_sec  = float(prof.min())
    nz_sec = int(np.count_nonzero(np.cumsum(drF) <= H_sec + 1e-6))
    return nz_sec, drF[:nz_sec], H_sec


def _leer_dx_seccion(run_dir, j_boca=J_BOCA, x_boca=X_BOCA):
    """dx [m] de cada celda de la sección de la boca, leído de dxF en la grilla.

    Necesario porque la malla es telescópica desde jun-2026: integrar con un DX
    constante sesga F_total. Devuelve (nx_sec,); si la grilla no trae dxF cae a DX.
    """
    tile_info, _, _ = build_tile_map(run_dir, MNC_GLOB, "grid.t*.nc")
    nx_sec = x_boca.stop - x_boca.start
    dx_sec = np.full(nx_sec, np.nan)
    for tile in tiles_for_section(tile_info, j_boca, x_boca):
        with Dataset(tile["file"]) as dg:
            if "dxF" not in dg.variables:
                return np.full(nx_sec, DX)
            dxF = dg.variables["dxF"][:].data
        xl       = tile["x_local_slice"]
        xg_start = tile["x_global_start"] - x_boca.start
        dx_sec[xg_start:xg_start + (xl.stop - xl.start)] = dxF[tile["j_local"], xl]
    if np.any(np.isnan(dx_sec)):
        dx_sec = np.where(np.isnan(dx_sec), DX, dx_sec)
    return dx_sec


# ══════════════════════════════════════════════════════════════════════
# 1. EOS lineal
# ══════════════════════════════════════════════════════════════════════
def eos(T, S, tRef, sRef, tAlpha, sBeta, rhoNil):
    """ρ = rhoNil·(1 − tAlpha·(T−tRef) + sBeta·(S−sRef))"""
    return rhoNil * (1.0 - tAlpha * (T - tRef) + sBeta * (S - sRef))


# ══════════════════════════════════════════════════════════════════════
# 2. Mapa de tiles
# ══════════════════════════════════════════════════════════════════════
def build_tile_map(run_dir, mnc_glob, stage_glob):
    """
    Encuentra todos los archivos de estado y construye un mapa de tiles.
    Retorna (tile_info, Y_global, X_global).
    """
    files = []
    for d in sorted(glob.glob(str(Path(run_dir) / mnc_glob))):
        files.extend(glob.glob(str(Path(d) / stage_glob)))
    if not files:
        return [], None, None

    tile_raw = []
    for fi in files:
        with Dataset(fi) as ds:
            Y = ds.variables["Y"][:].data.copy()
            X = ds.variables["X"][:].data.copy()
        tile_raw.append({"file": fi, "Y": Y, "X": X})

    Y_global = np.unique(np.concatenate([t["Y"] for t in tile_raw]))
    X_global = np.unique(np.concatenate([t["X"] for t in tile_raw]))

    tile_info = []
    for t in tile_raw:
        j0 = int(np.searchsorted(Y_global, t["Y"][0]))
        i0 = int(np.searchsorted(X_global, t["X"][0]))
        tile_info.append({
            "file": t["file"],
            "j0":   j0,
            "i0":   i0,
            "ny":   len(t["Y"]),
            "nx":   len(t["X"]),
        })

    tile_info.sort(key=lambda t: (t["j0"], t["i0"]))
    return tile_info, Y_global, X_global


def tiles_for_section(tile_info, j_global, x_global_slice):
    """
    Devuelve los tiles que contienen la fila j_global y las columnas
    x_global_slice, con los índices locales dentro de cada tile.
    """
    x_start = x_global_slice.start or 0
    x_stop  = x_global_slice.stop

    result = []
    for t in tile_info:
        j_end = t["j0"] + t["ny"] - 1
        i_end = t["i0"] + t["nx"] - 1

        if not (t["j0"] <= j_global <= j_end):
            continue

        x_ov_start = max(x_start, t["i0"])
        x_ov_end   = min(x_stop,  i_end + 1)
        if x_ov_start >= x_ov_end:
            continue

        result.append({
            "file"          : t["file"],
            "j_local"       : j_global - t["j0"],
            "x_local_slice" : slice(x_ov_start - t["i0"], x_ov_end - t["i0"]),
            "x_global_start": x_ov_start,
            "x_global_end"  : x_ov_end,
        })

    result.sort(key=lambda r: r["x_global_start"])
    return result


def _leer_drF(run_dir):
    """Lee drF (espesores verticales) del primer archivo grid disponible."""
    grid_files = []
    for d in sorted(glob.glob(str(Path(run_dir) / MNC_GLOB))):
        grid_files.extend(glob.glob(str(Path(d) / "grid.t*.nc")))
    if not grid_files:
        raise FileNotFoundError(f"No se encontraron archivos grid.t*.nc en {run_dir}")
    grid_files.sort(key=lambda p: int(Path(p).stem.split(".t")[-1]))
    with Dataset(grid_files[0]) as dg:
        return dg.variables["drF"][:].data.astype(float)  # (nz,)


def _interpolar_UV_Tpoint(U_tile, V_tile):
    """
    Interpola U, V de la grilla Arakawa-C al punto T.
    U_tile: (..., ny, nx+1)  →  U_T: (..., ny, nx)
    V_tile: (..., ny+1, nx)  →  V_T: (..., ny, nx)
    Si ya están en T-points (ny,nx) los devuelve sin modificar.
    """
    nx_T = V_tile.shape[-1]   # T-points in x  (V no tiene columna extra)
    ny_T = U_tile.shape[-2]   # T-points in y  (U no tiene fila extra)

    U_T = (0.5 * (U_tile[..., :-1] + U_tile[..., 1:])
           if U_tile.shape[-1] == nx_T + 1 else U_tile)

    V_T = (0.5 * (V_tile[..., :-1, :] + V_tile[..., 1:, :])
           if V_tile.shape[-2] == ny_T + 1 else V_tile)

    return U_T, V_T


# ══════════════════════════════════════════════════════════════════════
# 3. Paso 1–2 Nash: medias temporales sobre AMBAS etapas
# ══════════════════════════════════════════════════════════════════════
def calcular_medias_run(run_dir, stage_globs, tRef, sRef, tAlpha, sBeta, rhoNil,
                        j_boca=J_BOCA, x_boca=X_BOCA, nz_sec=None, verbose=True):
    """
    Primera pasada: calcula ρ̄(z) y ū_section(z,x) sobre AMBAS etapas.

    ρ̄(z)      = media temporal de ρ en la sección de la boca, ambas etapas
                 (Nash eq. 2: perfil de referencia LOCAL en la sección,
                  no promedio espacial del dominio completo — necesario cuando
                  hay gradientes horizontales de densidad, e.g. escenario realista)
    ū(z,x)    = media temporal de u,v en la sección de la boca, ambas etapas
                 (Nash eq. 5: ū(z) que se resta antes de remover la barotropa)

    Retorna
    -------
    Todos los perfiles se recortan a los nz_sec niveles húmedos de la sección.

    rho_ref    : (nz_sec,)           — perfil de densidad de referencia [kg/m³]
    u_mean_sec : (nz_sec, nx_sec)    — media temporal de U en la sección [m/s]
    v_mean_sec : (nz_sec, nx_sec)    — media temporal de V en la sección [m/s]
    drF        : (nz_sec,)           — espesores de capa [m]
    """
    if nz_sec is None:
        nz_sec, _, _ = niveles_humedos_seccion(run_dir, j_boca, x_boca)
    nx_sec   = x_boca.stop - x_boca.start
    rho_sum  = None   # (nz, nx_sec)
    N_rho    = 0
    u_sum    = None
    v_sum    = None
    N_uv     = 0
    drF_arr  = None

    masks = _mascaras_sin_duplicados(run_dir, stage_globs)

    for stage_glob, keep in zip(stage_globs, masks):
        tile_info, _, _ = build_tile_map(run_dir, MNC_GLOB, stage_glob)
        if not tile_info:
            if verbose:
                print(f"    [WARN] Sin archivos para {stage_glob}")
            continue
        if keep is None or not keep.any():
            continue

        if drF_arr is None:
            drF_arr = _leer_drF(run_dir)[:nz_sec]
            rho_sum = np.zeros((nz_sec, nx_sec))
            u_sum   = np.zeros((nz_sec, nx_sec))
            v_sum   = np.zeros((nz_sec, nx_sec))

        # ── ρ̄(z) y ū(z,x): acumular solo en los tiles de la sección ──
        sec_tiles = tiles_for_section(tile_info, j_boca, x_boca)
        if not sec_tiles:
            if verbose:
                print(f"    [WARN] Ningún tile cubre la sección para {stage_glob}")
            continue

        nt_stage = int(keep.sum())
        N_uv  += nt_stage
        N_rho += nt_stage

        for tile in sec_tiles:
            j_local  = tile["j_local"]
            xl       = tile["x_local_slice"]
            xg_start = tile["x_global_start"] - x_boca.start
            nx_piece = xl.stop - xl.start

            with Dataset(tile["file"]) as ds:
                # [:, :nz_sec] descarta las celdas secas bajo el fondo de la bahía
                T_tile = ds.variables["Temp"][:].data[keep][:, :nz_sec]  # (nt,nz_sec,ny,nx)
                S_tile = ds.variables["S"][:].data[keep][:, :nz_sec]
                U_tile = ds.variables["U"][:].data[keep][:, :nz_sec]     # (…, ny, nx+1)
                V_tile = ds.variables["V"][:].data[keep][:, :nz_sec]     # (…, ny+1, nx)

            # ρ̄(z): media temporal en la fila de la sección → (nz, nx_piece)
            rho_tile = eos(T_tile, S_tile, tRef, sRef, tAlpha, sBeta, rhoNil)
            rho_sum[:, xg_start:xg_start+nx_piece] += (
                rho_tile[:, :, j_local, xl].sum(axis=0)   # suma sobre tiempo
            )

            # ū(z,x): interpolar al punto T y acumular
            U_T, V_T = _interpolar_UV_Tpoint(U_tile, V_tile)
            u_sum[:, xg_start:xg_start+nx_piece] += U_T[:, :, j_local, xl].sum(axis=0)
            v_sum[:, xg_start:xg_start+nx_piece] += V_T[:, :, j_local, xl].sum(axis=0)

    if rho_sum is None:
        raise FileNotFoundError(f"No se encontraron datos de estado en {run_dir}")

    # Promedio temporal; rho_ref se colapsa a (nz,) promediando sobre x de la sección
    rho_ref    = (rho_sum / N_rho).mean(axis=1)   # (nz,)
    u_mean_sec = u_sum / N_uv                      # (nz, nx_sec)
    v_mean_sec = v_sum / N_uv

    if verbose:
        drF2 = drF_arr[:, np.newaxis]
        H    = drF_arr.sum()
        u_bar = (u_mean_sec * drF2).sum(axis=0) / H
        print(f"    ρ̄(z)      : {N_rho} pasos de tiempo | "
              f"sup={rho_ref[0]:.4f}  fondo={rho_ref[-1]:.4f}  "
              f"Δρ={rho_ref[-1]-rho_ref[0]:.5f} kg/m³")
        print(f"    ū_section : {N_uv} pasos de tiempo | "
              f"max|ū|={np.abs(u_mean_sec).max():.5f}  "
              f"max|v̄|={np.abs(v_mean_sec).max():.5f} m/s  "
              f"ū_baro={u_bar.mean():.5f} m/s")

    return rho_ref, u_mean_sec, v_mean_sec, drF_arr


# ══════════════════════════════════════════════════════════════════════
# 4. Pasos 3–5 Nash: presión baroclínica (eqs. 2–4)
# ══════════════════════════════════════════════════════════════════════
def presion_baroclínica(rho, rho_ref, Eta, drF, rhoNil):
    """
    Calcula p'(z,t) siguiendo Nash et al. eqs. 2-4.

    Nash eq. 2: ρ'(z,t) = ρ(z,t) − ρ̄(z)
    Nash eq. 3: p̃(z,t) = ρ₀·g·η + g·∫_z^0 ρ'(ẑ,t) dẑ
    Nash eq. 4: p'(z,t) = p̃ − H⁻¹·∫p̃ dz    (cond. baroclínica: ∫p' dz = 0)

    rho     : (nz_sec, ny, nx)  — densidad en un instante t, ya recortada a los
                                  niveles húmedos de la sección
    rho_ref : (nz_sec,)         — densidad de referencia ρ̄(z)
    Eta     : (ny, nx)          — elevación superficial η
    drF     : (nz_sec,)         — espesores de capa (k=0 en superficie). H=drF.sum()
                                  es la profundidad LOCAL de la sección.

    Retorna
        p_bc  : (nz, ny, nx)  [Pa]  presión baroclínica (∫p_bc dz = 0)
        p_bar : (ny, nx)      [Pa]  presión barotrópica = H⁻¹·∫p̃ dz
                                    (la parte que Nash resta; aquí la conservamos
                                     para el flujo barotrópico)
    """
    drF3 = drF[:, np.newaxis, np.newaxis]

    # Eq. 2: ρ'
    rho_prime = rho - rho_ref[:, np.newaxis, np.newaxis]

    # Eq. 3: p̃ = ρ₀gη + g·∫_z^0 ρ' dẑ
    #   int_techo[k] = Σ_{j=0}^{k-1} ρ'(j)·drF(j)  (capas ENCIMA de k)
    int_techo = np.cumsum(rho_prime * drF3, axis=0)
    int_techo = np.roll(int_techo, 1, axis=0)
    int_techo[0] = 0.0
    # Valor al centro de cada celda: sumar media celda actual
    p_tilde = (rhoNil * G * Eta[np.newaxis]
               + G * int_techo
               + G * rho_prime * drF3 / 2.0)

    # Eq. 4: p' = p̃ − H⁻¹·∫p̃ dz
    H     = drF.sum()
    p_bar = (p_tilde * drF3).sum(axis=0) / H   # (ny, nx) — presión barotrópica
    return p_tilde - p_bar[np.newaxis], p_bar   # (nz, ny, nx), (ny, nx)


# ══════════════════════════════════════════════════════════════════════
# 5. Pasos 6–7 Nash: velocidad baroclínica y flujo (eqs. 5, 10)
# ══════════════════════════════════════════════════════════════════════
def flujo_energia_seccion(U_raw, V_raw, p_bc, p_bar, drF, u_mean, v_mean, j_local, xl):
    """
    Nash eq. 5 : u'(z,t) = [u(z,t) − ū(z)] − H⁻¹·∫[u − ū] dz
    Nash eq. 10: F(t)    = ∫ p'(z,t)·u'(z,t) dz   [W/m]

    Barotrópico (para de-promediar en el tiempo aguas arriba):
        p_bt(t,x) = H⁻¹·∫p̃ dz   (= p_bar en la sección)
        u_bt(t,x) = H⁻¹·∫u dz ,  v_bt = H⁻¹·∫v dz

    U_raw  : (nz, ny, nx+1)   — un instante, grilla Arakawa-C
    V_raw  : (nz, ny+1, nx)
    p_bc   : (nz, ny, nx)     — presión baroclínica del tile completo
    p_bar  : (ny, nx)         — presión barotrópica del tile completo
    u_mean : (nz, nx_piece)   — media temporal de U en la sección (Nash ū(z))
    v_mean : (nz, nx_piece)
    j_local: fila de la sección en el tile
    xl     : slice de columnas en el tile

    Retorna
        Fx, Fy           : (nx_piece,)  [W/m]   flujo baroclínico instantáneo
        p_bt, u_bt, v_bt : (nx_piece,)          cantidades barotrópicas crudas
    """
    drF2 = drF[:, np.newaxis]   # (nz, 1)
    H    = drF.sum()

    # Interpolar U, V al punto T
    U_T, V_T = _interpolar_UV_Tpoint(U_raw, V_raw)

    # Extraer fila de la sección: (nz, nx_piece)
    u_sec = U_T[:, j_local, xl]
    v_sec = V_T[:, j_local, xl]
    p_sec = p_bc[:, j_local, xl]

    # Nash eq. 5 — dos pasos:
    #   (a) restar la media temporal ū(z) en cada nivel
    u_pert = u_sec - u_mean   # (nz, nx_piece)
    v_pert = v_sec - v_mean

    #   (b) restar la componente barotropa de la perturbación
    #       para que ∫u' dz = 0 (condición barocl., Nash eq. 6)
    u_bc = u_pert - (u_pert * drF2).sum(axis=0) / H   # (nz, nx_piece)
    v_bc = v_pert - (v_pert * drF2).sum(axis=0) / H

    # Nash eq. 10: F = ∫ p'·u' dz
    Fx = (p_sec * u_bc * drF2).sum(axis=0)   # (nx_piece,)
    Fy = (p_sec * v_bc * drF2).sum(axis=0)

    # Cantidades barotrópicas crudas (media vertical) en la sección
    u_bt = (u_sec * drF2).sum(axis=0) / H    # (nx_piece,)
    v_bt = (v_sec * drF2).sum(axis=0) / H
    p_bt = p_bar[j_local, xl]                # (nx_piece,)

    return Fx, Fy, p_bt, u_bt, v_bt


# ══════════════════════════════════════════════════════════════════════
# 6. Procesamiento de una etapa con medias pre-computadas
# ══════════════════════════════════════════════════════════════════════
def _procesar_chunk(run_dir, stage_glob, keep, rho_ref, u_mean_sec, v_mean_sec, drF,
                    tRef, sRef, tAlpha, sBeta, rhoNil, j_boca, x_boca, verbose):
    """Nash en la sección de la boca para UN chunk. Devuelve arrays crudos."""
    tile_info, _, _ = build_tile_map(run_dir, MNC_GLOB, stage_glob)
    if not tile_info:
        return None

    sec_tiles = tiles_for_section(tile_info, j_boca, x_boca)
    if not sec_tiles:
        raise ValueError(f"Ningún tile cubre j={j_boca}, x={x_boca}")

    with Dataset(sec_tiles[0]["file"]) as ds:
        t_all = ds.variables["T"][:].data.copy()
    t_all = t_all[keep]
    nt = len(t_all)
    if nt == 0:
        return None

    nx_sec = x_boca.stop - x_boca.start
    Fx_sec = np.zeros((nt, nx_sec))
    Fy_sec = np.zeros((nt, nx_sec))
    # Acumuladores barotrópicos crudos (se de-promedian en el tiempo al final,
    # sobre la etapa COMPLETA — no por chunk, que daría medias distintas)
    p_bt_sec = np.zeros((nt, nx_sec))
    u_bt_sec = np.zeros((nt, nx_sec))
    v_bt_sec = np.zeros((nt, nx_sec))

    for tile in sec_tiles:
        fi       = tile["file"]
        j_local  = tile["j_local"]
        xl       = tile["x_local_slice"]
        xg_start = tile["x_global_start"] - x_boca.start
        nx_piece = xl.stop - xl.start

        # Medias de velocidad para este trozo de sección
        u_mean_p = u_mean_sec[:, xg_start:xg_start+nx_piece]   # (nz, nx_piece)
        v_mean_p = v_mean_sec[:, xg_start:xg_start+nx_piece]

        # drF ya viene recortado a los niveles húmedos de la sección; recortar los
        # campos igual evita meter celdas secas (Temp=S=0) en ρ' y en la integral.
        nz_sec = len(drF)
        with Dataset(fi) as ds:
            T_tile   = ds.variables["Temp"][:].data[keep][:, :nz_sec]  # (nt,nz_sec,ny,nx)
            S_tile   = ds.variables["S"][:].data[keep][:, :nz_sec]
            Eta_tile = ds.variables["Eta"][:].data[keep]               # (nt, ny, nx)
            U_tile   = ds.variables["U"][:].data[keep][:, :nz_sec]     # (…, ny, nx+1)
            V_tile   = ds.variables["V"][:].data[keep][:, :nz_sec]     # (…, ny+1, nx)

        for ti in range(nt):
            rho_t = eos(T_tile[ti], S_tile[ti], tRef, sRef, tAlpha, sBeta, rhoNil)
            p_bc, p_bar = presion_baroclínica(rho_t, rho_ref, Eta_tile[ti], drF, rhoNil)
            Fx, Fy, p_bt, u_bt, v_bt = flujo_energia_seccion(
                U_tile[ti], V_tile[ti], p_bc, p_bar, drF,
                u_mean_p, v_mean_p, j_local, xl,
            )
            sl = slice(xg_start, xg_start + nx_piece)
            Fx_sec[ti, sl]   = Fx
            Fy_sec[ti, sl]   = Fy
            p_bt_sec[ti, sl] = p_bt
            u_bt_sec[ti, sl] = u_bt
            v_bt_sec[ti, sl] = v_bt

    if verbose:
        print(f"      {stage_glob.split('.')[1]}: {nt:3d} pasos, "
              f"{len(sec_tiles)} tile(s) en la sección")
    return t_all, Fx_sec, Fy_sec, p_bt_sec, u_bt_sec, v_bt_sec


def procesar_etapa(run_dir, stage_globs, rho_ref, u_mean_sec, v_mean_sec, drF,
                   tRef, sRef, tAlpha, sBeta, rhoNil,
                   j_boca=J_BOCA, x_boca=X_BOCA, verbose=True):
    """
    Aplica Nash et al. (2005) en la sección de la boca para UNA etapa,
    usando las medias temporales pre-computadas sobre AMBAS etapas.

    `stage_globs` es la LISTA de chunks que componen la etapa (6 en la corrida
    de jun-2026). Se procesan uno a uno y se concatenan; el de-promediado
    barotrópico se hace al final, sobre la etapa completa.

    Retorna dict con t_seconds, t_horas, Fx_sec, Fy_sec (arrays instantáneos).
    """
    if isinstance(stage_globs, str):
        stage_globs = [stage_globs]

    masks = _mascaras_sin_duplicados(run_dir, stage_globs)
    if verbose:
        print(f"    sección j={j_boca}, x={x_boca.start}:{x_boca.stop} | "
              f"{len(stage_globs)} chunk(s)")

    piezas = []
    for g, keep in zip(stage_globs, masks):
        if keep is None:
            continue
        res = _procesar_chunk(run_dir, g, keep, rho_ref, u_mean_sec, v_mean_sec,
                              drF, tRef, sRef, tAlpha, sBeta, rhoNil,
                              j_boca, x_boca, verbose)
        if res is not None:
            piezas.append(res)

    if not piezas:
        raise FileNotFoundError(f"No se encontraron archivos {stage_globs} en {run_dir}")

    t_all    = np.concatenate([p[0] for p in piezas])
    Fx_sec   = np.concatenate([p[1] for p in piezas], axis=0)
    Fy_sec   = np.concatenate([p[2] for p in piezas], axis=0)
    p_bt_sec = np.concatenate([p[3] for p in piezas], axis=0)
    u_bt_sec = np.concatenate([p[4] for p in piezas], axis=0)
    v_bt_sec = np.concatenate([p[5] for p in piezas], axis=0)

    idx = np.argsort(t_all, kind="stable")
    t_all    = t_all[idx]
    Fx_sec   = Fx_sec[idx]
    Fy_sec   = Fy_sec[idx]
    p_bt_sec = p_bt_sec[idx]
    u_bt_sec = u_bt_sec[idx]
    v_bt_sec = v_bt_sec[idx]

    # Flujo barotrópico: F_bt = p_bt'·u_bt'·H con perturbaciones temporales de la etapa
    H = drF.sum()
    pbt_p = p_bt_sec - p_bt_sec.mean(axis=0, keepdims=True)
    ubt_p = u_bt_sec - u_bt_sec.mean(axis=0, keepdims=True)
    vbt_p = v_bt_sec - v_bt_sec.mean(axis=0, keepdims=True)
    Fx_bt_sec = pbt_p * ubt_p * H   # (nt, nx_sec) [W/m]
    Fy_bt_sec = pbt_p * vbt_p * H

    return {
        "t_seconds":  t_all,
        "t_horas":    t_all / 3600.0,
        "Fx_sec":     Fx_sec,     # (nt, nx_sec) [W/m] baroclínico
        "Fy_sec":     Fy_sec,
        "Fx_bt_sec":  Fx_bt_sec,  # (nt, nx_sec) [W/m] barotrópico
        "Fy_bt_sec":  Fy_bt_sec,
    }


# ══════════════════════════════════════════════════════════════════════
# 7. Utilidades de análisis post-hoc
# ══════════════════════════════════════════════════════════════════════
T_SEICHE_S = 14430.0   # s — periodo observado del seiche barotropo (240.5 min)


def _media_movil(x, N):
    """Media móvil centrada de ventana N sobre array 1D. Bordes → NaN."""
    kernel = np.ones(N) / float(N)
    out    = np.convolve(x.astype(float), kernel, mode="same")
    half   = N // 2
    if half > 0:
        out[:half]  = np.nan
        out[-half:] = np.nan
    return out


def promediar_flujo_periodos(Fsec, t_seconds, T_avg_s=T_SEICHE_S):
    """
    Aplica media móvil centrada sobre el flujo instantáneo Fsec.

    Útil para suavizar series temporales en análisis y figuras.
    Nash et al. (2005) eq. (10) define ⟨F⟩ como la media temporal simple;
    esta función es un complemento para explorar la evolución temporal suavizada.

    Parámetros
    ----------
    Fsec      : (nt, nx) o (nt,)
    t_seconds : (nt,)
    T_avg_s   : float — ventana [s]

    Retorna F_avg (mismo shape), N_avg (int), dt_s (float)
    """
    dt_s  = float(np.median(np.diff(t_seconds)))
    N_avg = max(1, round(T_avg_s / dt_s))

    if Fsec.ndim == 1:
        F_avg = _media_movil(Fsec, N_avg)
    else:
        F_avg = np.apply_along_axis(_media_movil, 0, Fsec, N_avg)

    return F_avg, N_avg, dt_s


def cargar_flujo(scenario, run, etapa, npz_dir=None):
    """
    Carga un archivo .npz de flujo de energía (un solo archivo por escenario/run).

    etapa : 1 (forzado) o 2 (libre)
    Retorna dict con claves de la etapa pedida + arrays compartidos + 'meta'.
    """
    if npz_dir is None:
        npz_dir = Path(__file__).resolve().parent / "npz_flux"
    fname = Path(npz_dir) / f"{scenario}_{run}_flujo.npz"
    d    = np.load(fname, allow_pickle=True)
    e    = f"e{etapa}"
    meta_all = json.loads(str(d["meta_json"]))
    meta = meta_all.get(f"etapa{etapa}", meta_all)
    out = {
        "t_seconds" : d[f"t_seconds_{e}"],
        "t_horas"   : d[f"t_horas_{e}"],
        "Fx_sec"    : d[f"Fx_sec_{e}"],
        "Fy_sec"    : d[f"Fy_sec_{e}"],
        "Fx_mean"   : d[f"Fx_mean_{e}"],
        "Fy_mean"   : d[f"Fy_mean_{e}"],
        "rho_ref"   : d["rho_ref"],
        "u_mean_sec": d["u_mean_sec"],
        "v_mean_sec": d["v_mean_sec"],
        "drF"       : d["drF"],
        # dx real de cada columna de la sección (malla telescópica). Ausente en
        # archivos anteriores a ago-2026 → None, y quien integre cae a DX.
        "dx_sec"    : d["dx_sec"] if "dx_sec" in d.files else None,
        "meta"      : meta,
    }
    # Claves barotrópicas (presentes solo en archivos regenerados)
    for k in (f"Fx_bt_sec_{e}", f"Fy_bt_sec_{e}", f"Fx_bt_mean_{e}", f"Fy_bt_mean_{e}"):
        out[k.replace(f"_{e}", "")] = d[k] if k in d.files else None
    return out


def integrar_flujo_bahia(scenario, run, etapa, dx=None, npz_dir=None, T_avg_s=None):
    """
    Integra el flujo de energía baroclínico a través de la boca de la bahía.

        F_total(t) [W] = Σₓ Fy_sec(t, x) · dx(x)

    Positivo = flujo neto entrando a la bahía (hacia +y).

    Parámetros
    ----------
    scenario : '2layer', 'realistic', 'linear'
    run      : 'bay', 'nobay'
    etapa    : 1 (forzado) o 2 (libre)
    dx       : espaciado horizontal [m]. Por defecto usa el dx_sec guardado en el
               npz (malla telescópica); sólo cae a DX si el archivo no lo trae.
    T_avg_s  : si se indica, aplica media móvil centrada de esa ventana [s]
               antes de integrar (útil para visualización, no para el estimado medio)

    Retorna dict con t_horas, F_total, Fy_sec, rho_ref, meta.
    """
    d  = cargar_flujo(scenario, run, etapa, npz_dir=npz_dir)
    Fy = d["Fy_sec"]

    if dx is None:
        dx = d["dx_sec"] if d["dx_sec"] is not None else DX

    if T_avg_s is not None:
        Fy, _, _ = promediar_flujo_periodos(Fy, d["t_seconds"], T_avg_s=T_avg_s)

    F_total = np.nansum(Fy * np.asarray(dx), axis=1)
    filas_nan = np.any(np.isnan(Fy), axis=1)
    F_total[filas_nan] = np.nan

    return {
        "t_horas": d["t_horas"],
        "F_total": F_total,
        "Fy_sec":  Fy,
        "rho_ref": d["rho_ref"],
        "meta":    d["meta"],
    }


def resumen_flujos(dx=None, npz_dir=None):
    """Imprime resumen de flujo total medio y pico para todos los casos disponibles."""
    scenarios = ["2layer", "realistic", "linear"]
    runs      = ["bay", "nobay"]

    header = f"{'Escenario':<12} {'Run':<7} {'Etapa':<7} "
    header += f"{'<F> [kW]':>10}  {'max|F| [kW]':>12}  {'t_pico [h]':>11}"
    print(header)
    print("-" * len(header))

    for sc in scenarios:
        for ru in runs:
            for et in [1, 2]:
                try:
                    r = integrar_flujo_bahia(sc, ru, et, dx=dx, npz_dir=npz_dir)
                except FileNotFoundError:
                    continue
                F  = r["F_total"] / 1e3   # → kW
                t  = r["t_horas"]
                ip = np.nanargmax(np.abs(F))
                print(f"{sc:<12} {ru:<7} {et:<7} "
                      f"{np.nanmean(F):>10.2f}  {np.nanmax(np.abs(F)):>12.2f}  "
                      f"{t[ip]:>11.1f}")


# ══════════════════════════════════════════════════════════════════════
# 8. Bucle principal
# ══════════════════════════════════════════════════════════════════════
def main():
    SAVE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Guardando en: {SAVE_DIR}\n")

    for scen_name, cfg in SCENARIOS.items():
        base = cfg["base"]
        if not base.exists():
            print(f"[SKIP] {scen_name}: ruta no disponible ({base})\n")
            continue

        for run in RUNS:
            run_dir   = base / run
            run_label = "bay" if run == "run_expand" else "nobay"

            print(f"{'='*60}")
            print(f"Escenario: {scen_name}  |  Run: {run_label}")
            print(f"  Ruta: {run_dir}")

            try:
                stage_groups = discover_stage_groups(run_dir)
            except FileNotFoundError as e:
                print(f"  [SKIP] {e}\n")
                continue
            todos_globs = [g for grp in stage_groups for g in grp]
            print(f"  [chunks: {len(todos_globs)} → etapa1={len(stage_groups[0])}, "
                  f"etapa2={len(stage_groups[1])}]")

            dx_sec = _leer_dx_seccion(run_dir)
            print(f"  [dx en la sección: {dx_sec.min():.1f}–{dx_sec.max():.1f} m "
                  f"(DX fijo era {DX:.0f} m)]")

            try:
                nz_sec, _, H_sec = niveles_humedos_seccion(run_dir)
            except ValueError as e:
                print(f"  [SKIP] {e}\n")
                continue
            print(f"  [sección j={J_BOCA}, i={X_BOCA.start}..{X_BOCA.stop-1}: "
                  f"{H_sec:.1f} m → {nz_sec} niveles húmedos]")

            # ── Paso 1: medias sobre AMBAS etapas (0–120 h) ───────────
            print("  Paso 1: calculando ρ̄(z) y ū(z,x) sobre ambas etapas ...")
            try:
                rho_ref, u_mean, v_mean, drF = calcular_medias_run(
                    run_dir, todos_globs,
                    tRef=cfg["tRef"], sRef=cfg["sRef"],
                    tAlpha=cfg["tAlpha"], sBeta=cfg["sBeta"],
                    rhoNil=cfg["rhoNil"],
                    nz_sec=nz_sec,
                    verbose=True,
                )
            except (FileNotFoundError, ValueError) as e:
                print(f"  [SKIP] {e}\n")
                continue

            # ── Paso 2: procesar cada etapa con medias comunes ─────────
            stage_results = {}
            for etapa_n, stage_glob in enumerate(stage_groups, start=1):
                stage_label = "forzado" if etapa_n == 1 else "libre"
                print(f"\n  Etapa {etapa_n} ({stage_label}) — {len(stage_glob)} chunk(s)")
                try:
                    res = procesar_etapa(
                        run_dir, stage_glob,
                        rho_ref, u_mean, v_mean, drF,
                        tRef=cfg["tRef"], sRef=cfg["sRef"],
                        tAlpha=cfg["tAlpha"], sBeta=cfg["sBeta"],
                        rhoNil=cfg["rhoNil"],
                        verbose=True,
                    )
                except (FileNotFoundError, ValueError) as e:
                    print(f"  [SKIP] {e}")
                    continue
                stage_results[etapa_n] = (res, stage_label)

                t0 = res["t_horas"][0]; t1 = res["t_horas"][-1]
                print(f"    t=[{t0:.1f},{t1:.1f}]h  nt={len(res['t_horas'])}  "
                      f"<Fy>={float(np.nanmean(res['Fy_sec'])):.3f} W/m")

            if not stage_results:
                continue

            # ── Guardar ambas etapas en un único archivo ───────────────
            out_name = f"{scen_name}_{run_label}_flujo.npz"
            out_path = SAVE_DIR / out_name

            meta_common = dict(
                scenario      = scen_name,
                run           = run_label,
                j_boca        = J_BOCA,
                x_boca_start  = X_BOCA.start,
                x_boca_stop   = X_BOCA.stop,
                method        = "Nash, Alford & Kunze (2005) + flujo barotrópico p_bt'·u_bt'·H",
                units_F       = "W/m",
                bt_demean     = "perturbaciones barotrópicas de-promediadas por etapa",
                rho_ref_scope = f"ambas etapas (0-120 h), sección de la boca (j={J_BOCA})",
                u_mean_scope  = "ambas etapas (0-120 h), sección de la boca",
                n_chunks      = len(todos_globs),
                dx_source     = "dxF de grid.t*.nc (malla telescópica, jun-2026)",
                nz_sec        = nz_sec,
                H_sec         = H_sec,
                z_scope       = f"0-{H_sec:.1f} m ({nz_sec} niveles húmedos de la boca)",
            )
            for etapa_n, (_, stage_label) in stage_results.items():
                meta_common[f"etapa{etapa_n}"] = dict(etapa=etapa_n, stage=stage_label)

            save_kwargs = dict(
                rho_ref    = rho_ref,
                u_mean_sec = u_mean,
                v_mean_sec = v_mean,
                drF        = drF,
                dx_sec     = dx_sec,
                meta_json  = json.dumps(meta_common),
            )
            for etapa_n, (res, _) in stage_results.items():
                e = f"e{etapa_n}"
                save_kwargs[f"t_seconds_{e}"] = res["t_seconds"]
                save_kwargs[f"t_horas_{e}"]   = res["t_horas"]
                save_kwargs[f"Fx_sec_{e}"]    = res["Fx_sec"]
                save_kwargs[f"Fy_sec_{e}"]    = res["Fy_sec"]
                save_kwargs[f"Fx_mean_{e}"]   = res["Fx_sec"].mean(axis=0)
                save_kwargs[f"Fy_mean_{e}"]   = res["Fy_sec"].mean(axis=0)
                save_kwargs[f"Fx_bt_sec_{e}"] = res["Fx_bt_sec"]
                save_kwargs[f"Fy_bt_sec_{e}"] = res["Fy_bt_sec"]
                save_kwargs[f"Fx_bt_mean_{e}"]= res["Fx_bt_sec"].mean(axis=0)
                save_kwargs[f"Fy_bt_mean_{e}"]= res["Fy_bt_sec"].mean(axis=0)

            np.savez(out_path, **save_kwargs)
            print(f"\n  [OK] {out_name}  ({len(stage_results)} etapa(s))")

    print("\n[Listo] energy_flux_IWBay.py completado.")


if __name__ == "__main__":
    main()
