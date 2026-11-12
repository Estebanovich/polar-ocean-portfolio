#!/usr/bin/env python3
"""
conversion_baroclinica_2D.py
=============================
Mapa 2-D horizontal de la tasa de conversión de energía barotrópica →
baroclínica (ondas internas):

    C(x,y) = < ∫ ρ'(x,y,z,t) · g · w(x,y,z,t) dz >_t        [W/m²]

Forma depth-integrada del término de flotabilidad en la ecuación de energía
baroclínica (Kang & Fringer 2012; Omidvar et al. 2021). Es forzamiento-
agnóstica: no asume marea ni un único armónico, por lo que aplica igual al
caso IWBay, donde el "barotrópico" es el seiche de cuenca excitado por un
PULSO de viento (no forzamiento periódico sostenido — ver
make_wind_forcing_local.ipynb: pulso gaussiano cortado en t≈9000 s / 2.5 h,
NO durante toda la etapa 1 de 216000 s).

ρ'(x,y,z,t) = ρ(x,y,z,t) − ρ̄(x,y,z)   con ρ̄ = media temporal por columna
sobre TODO el run (todos los chunks descubiertos por discover_stage_patterns).

w se lee de ``W`` (definido en interfaces ``Zl``, tope de cada celda) y se
interpola al centro de celda: w_c(k) = ½[W(k) + W(k+1)], con W(fondo)=0.

La integral vertical usa el espesor real de celda drF·hFacC (respeta
celdas parciales de fondo y tierra).

Ventanas temporales (ver memoria de proyecto sobre el pulso de 2.5 h):
    pulso_viento  : [0, 9000) s    — forzamiento activo
    libre         : [9000, fin] s  — decaimiento libre completo (fusiona lo que
                    antes eran "libre_etapa1" y "etapa2": el corte en 216000 s
                    era solo el fin nominal de la etapa 1 de forzamiento, sin
                    significado físico real ya que el viento se apaga en
                    t=9000s -- ver make_wind_forcing_local.ipynb -- así que
                    separar ahí no aislaba nada distinto, solo partía en dos
                    la misma fase de decaimiento libre)

Opcional (--band-seiche): además del C(x,y) de banda ancha por ventana, se
puede calcular C_seiche(x,y) filtrando ρ' y w en la banda del seiche
(SEICHE_BAND_CPH, igual que energy_flux_cg.py) ANTES de multiplicarlos,
sobre toda la señal libre (t≥9000 s) concatenada. Esto aísla la conversión
específicamente atribuible a la oscilación coherente del seiche, separada
de transientes de banda ancha (p.ej. el propio corte del pulso de viento) o
de armónicos superiores. El filtrado FFT necesita la serie temporal COMPLETA
por columna, así que este paso procesa nivel por nivel (no columna 3-D
completa) para mantener acotado el uso de memoria (ver calcular_C_seiche).

Salida: npz_conversion/{escenario}_{run}_conversion.npz
  claves: X, Y, Depth, C_<ventana>, N_<ventana> (nº de instantes promediados),
          t_seconds, C_dom_t (serie integrada en todo el dominio, W),
          [si --band-seiche] C_seiche, N_seiche, banda_seiche_cph,
          drF, meta_json

Uso:
    /opt/anaconda3/bin/python3 conversion_baroclinica_2D.py [--scen realistic linear 2layer]
                                                             [--run run_expand run_expand_nobay]
                                                             [--band-seiche] [--plot]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from energy_flux_IWBay import SCENARIOS, MNC_GLOB, G, eos, build_tile_map
from MIT_DatProcesW import discover_stage_patterns, _pick_var
from prototipo_hovmoller_cpcg import bandpass_time

SAVE_DIR = Path(__file__).resolve().parent / "npz_conversion"
FIG_DIR = Path(__file__).resolve().parent / "figs" / "conversion_baroclinica"
RUNS = ["run_expand", "run_expand_nobay"]
# NOTA: la grilla NO es Cartesiana uniforme (está afinada cerca de la bahía),
# por lo que la integral de dominio usa el área real de celda rAc, no un DX fijo.

# Ventanas temporales [s]: (etiqueta, t_ini incl., t_fin excl. o None=hasta el final)
# t=9000 s (2.5 h) es el corte real del pulso de viento gaussiano, no el fin
# nominal de la etapa 1 (216000 s) — ver make_wind_forcing_local.ipynb. Antes
# se partía el decaimiento libre en "libre_etapa1"/"etapa2" en t=216000s, pero
# ese corte no tiene significado físico (el viento ya se apagó en t=9000s) --
# se fusionaron en una sola ventana "libre".
T_PULSO_FIN = 9000.0
WINDOWS = [
    ("pulso_viento", 0.0, T_PULSO_FIN),
    ("libre", T_PULSO_FIN, None),
]

SEICHE_BAND_CPH = (0.235, 0.290)  # igual que energy_flux_cg.py (T_seiche≈14430s≈0.25cph)


# ══════════════════════════════════════════════════════════════════════
# 1. Grilla global (drF, hFacC, Depth) ensamblada a partir de los tiles
# ══════════════════════════════════════════════════════════════════════
def leer_grid_global(run_dir):
    tile_info, Y_global, X_global = build_tile_map(run_dir, MNC_GLOB, "grid.t*.nc")
    if not tile_info:
        raise FileNotFoundError(f"No se encontraron grid.t*.nc en {run_dir}")

    ny, nx = len(Y_global), len(X_global)
    with Dataset(tile_info[0]["file"]) as ds:
        drF = ds.variables["drF"][:].data.astype(float)
    nz = len(drF)

    hFacC = np.zeros((nz, ny, nx))
    Depth = np.zeros((ny, nx))
    rAc = np.zeros((ny, nx))
    for t in tile_info:
        j0, i0, tny, tnx = t["j0"], t["i0"], t["ny"], t["nx"]
        with Dataset(t["file"]) as ds:
            hfacc_name = _pick_var(ds, ("hFacC", "HFacC"))
            ra_name = _pick_var(ds, ("rAc", "rA"))
            hFacC[:, j0:j0 + tny, i0:i0 + tnx] = ds.variables[hfacc_name][:].data
            Depth[j0:j0 + tny, i0:i0 + tnx] = ds.variables["Depth"][:].data
            rAc[j0:j0 + tny, i0:i0 + tnx] = ds.variables[ra_name][:].data

    return drF, hFacC, Depth, rAc, Y_global, X_global


# ══════════════════════════════════════════════════════════════════════
# 2. W en interfaces (Zl) → centro de celda
# ══════════════════════════════════════════════════════════════════════
def w_a_centro(W, axis):
    """w_c(k) = ½[W(k) + W(k+1)], W(fondo)=0. W definido en Zl (tope de celda)."""
    W_fondo = np.zeros_like(W)
    sl_src = [slice(None)] * W.ndim
    sl_dst = [slice(None)] * W.ndim
    sl_src[axis] = slice(1, None)
    sl_dst[axis] = slice(None, -1)
    W_fondo[tuple(sl_dst)] = W[tuple(sl_src)]
    return 0.5 * (W + W_fondo)


# ══════════════════════════════════════════════════════════════════════
# 3. Paso 1: ρ̄(x,y,z) — media temporal por columna sobre TODO el run
# ══════════════════════════════════════════════════════════════════════
def calcular_rho_media(run_dir, stage_globs, drF, ny, nx, tRef, sRef, tAlpha, sBeta,
                       rhoNil, verbose=True):
    nz = len(drF)
    rho_sum = np.zeros((nz, ny, nx))
    n_t = 0

    for stage_glob in stage_globs:
        tile_info, _, _ = build_tile_map(run_dir, MNC_GLOB, stage_glob)
        if not tile_info:
            continue
        with Dataset(tile_info[0]["file"]) as ds:
            n_t += int(ds.variables["T"].shape[0])
        for tile in tile_info:
            j0, i0, tny, tnx = tile["j0"], tile["i0"], tile["ny"], tile["nx"]
            with Dataset(tile["file"]) as ds:
                T_t = ds.variables["Temp"][:].data
                S_t = ds.variables["S"][:].data
            rho_t = eos(T_t, S_t, tRef, sRef, tAlpha, sBeta, rhoNil)  # (nt,nz,tny,tnx)
            rho_sum[:, j0:j0 + tny, i0:i0 + tnx] += rho_t.sum(axis=0)

    if n_t == 0:
        raise FileNotFoundError(f"No hay datos de estado en {run_dir}")

    rho_mean = rho_sum / n_t
    if verbose:
        print(f"    ρ̄(x,y,z): {n_t} pasos de tiempo acumulados sobre "
              f"{len(stage_globs)} chunk(s)")
    return rho_mean


# ══════════════════════════════════════════════════════════════════════
# 4. Paso 2: C(x,y,t) = g·Σ_z ρ'·w_c·(drF·hFacC), acumulado por ventana
# ══════════════════════════════════════════════════════════════════════
def procesar_run(run_dir, stage_globs, drF, hFacC, rAc, ny, nx, rho_mean,
                 tRef, sRef, tAlpha, sBeta, rhoNil, verbose=True):
    nz = len(drF)
    peso = drF[:, None, None] * hFacC  # espesor real de celda (nz,ny,nx)

    C_sum = {w[0]: np.zeros((ny, nx)) for w in WINDOWS}
    N_win = {w[0]: 0 for w in WINDOWS}
    t_chunks, Cdom_chunks = [], []

    for stage_glob in stage_globs:
        tile_info, _, _ = build_tile_map(run_dir, MNC_GLOB, stage_glob)
        if not tile_info:
            continue
        with Dataset(tile_info[0]["file"]) as ds:
            t_all = ds.variables["T"][:].data.copy()
        nt = len(t_all)

        C_stage = np.zeros((nt, ny, nx))
        for tile in tile_info:
            j0, i0, tny, tnx = tile["j0"], tile["i0"], tile["ny"], tile["nx"]
            with Dataset(tile["file"]) as ds:
                T_t = ds.variables["Temp"][:].data
                S_t = ds.variables["S"][:].data
                W_t = ds.variables["W"][:].data
            rho_t = eos(T_t, S_t, tRef, sRef, tAlpha, sBeta, rhoNil)
            rho_prime = rho_t - rho_mean[None, :, j0:j0 + tny, i0:i0 + tnx]
            w_c = w_a_centro(W_t, axis=1)
            peso_p = peso[:, j0:j0 + tny, i0:i0 + tnx]
            C_stage[:, j0:j0 + tny, i0:i0 + tnx] = (
                G * np.sum(rho_prime * w_c * peso_p[None], axis=1)
            )

        Cdom_chunks.append(np.array([np.sum(C_stage[k] * rAc) for k in range(nt)]))
        t_chunks.append(t_all)

        for label, t0, t1 in WINDOWS:
            mask = (t_all >= t0) if t1 is None else ((t_all >= t0) & (t_all < t1))
            if not np.any(mask):
                continue
            C_sum[label] += C_stage[mask].sum(axis=0)
            N_win[label] += int(mask.sum())

        if verbose:
            print(f"    {stage_glob}: nt={nt}  t=[{t_all[0]:.0f},{t_all[-1]:.0f}] s")

    C_mean = {
        label: (C_sum[label] / N_win[label] if N_win[label] > 0
                else np.full((ny, nx), np.nan))
        for label, _, _ in WINDOWS
    }

    t_seconds = np.concatenate(t_chunks)
    order = np.argsort(t_seconds)
    t_seconds = t_seconds[order]
    Cdom_t = np.concatenate(Cdom_chunks)[order]

    return C_mean, N_win, t_seconds, Cdom_t


# ══════════════════════════════════════════════════════════════════════
# 4b. C_seiche(x,y): ρ' y w filtrados en la banda del seiche antes de
#     multiplicar, sobre la señal libre completa (t≥T_PULSO_FIN)
# ══════════════════════════════════════════════════════════════════════
def construir_series_por_tile(run_dir, mnc_glob, stage_globs):
    """Agrupa por tile ESPACIAL (mismo j0,i0) los archivos de estado de
    TODOS los chunks temporales, ordenados por nombre (= orden temporal,
    ya que solo cambia el nº de iteración con ancho de dígitos fijo)."""
    grupos = {}
    for sg in stage_globs:
        tile_info, _, _ = build_tile_map(run_dir, mnc_glob, sg)
        for t in tile_info:
            key = (t["j0"], t["i0"])
            g = grupos.setdefault(key, dict(j0=t["j0"], i0=t["i0"],
                                            ny=t["ny"], nx=t["nx"], files=[]))
            g["files"].append(t["file"])
    for g in grupos.values():
        g["files"].sort()
    return list(grupos.values())


def calcular_C_seiche(run_dir, stage_globs, drF, hFacC, ny, nx, rho_mean,
                      tRef, sRef, tAlpha, sBeta, rhoNil, t_min,
                      band_cph=SEICHE_BAND_CPH, verbose=True):
    """C_seiche(x,y) = <∫ ρ'_bp·g·w_bp dz>_t, con ρ', w filtrados en la banda
    del seiche ANTES de multiplicar (aísla la generación coherente del
    ruido/transientes de banda ancha), sobre t≥t_min concatenado.

    El filtro FFT (bandpass_time) necesita la serie temporal COMPLETA por
    columna en memoria — inviable para todo el dominio 3-D a la vez (decenas
    de GB). Se procesa NIVEL por NIVEL (loop en z): por cada tile y cada
    nivel k se leen solo esas 2 rebanadas de T,S,W de cada archivo (unos
    cientos de MB como mucho), se filtra, se multiplica y se acumula.
    """
    nz = len(drF)
    grupos = construir_series_por_tile(run_dir, MNC_GLOB, stage_globs)

    C = np.zeros((ny, nx))
    n_usado = 0

    for gi, g in enumerate(grupos):
        j0, i0, tny, tnx = g["j0"], g["i0"], g["ny"], g["nx"]
        rho_mean_p = rho_mean[:, j0:j0 + tny, i0:i0 + tnx]
        hFacC_p = hFacC[:, j0:j0 + tny, i0:i0 + tnx]

        for k in range(nz):
            T_parts, S_parts, Wtop_parts, Wbot_parts, t_parts = [], [], [], [], []
            for fi in g["files"]:
                with Dataset(fi) as ds:
                    t = ds.variables["T"][:].data
                    mask = t >= t_min
                    if not np.any(mask):
                        continue
                    T_parts.append(ds.variables["Temp"][mask, k].data)
                    S_parts.append(ds.variables["S"][mask, k].data)
                    Wtop_parts.append(ds.variables["W"][mask, k].data)
                    if k + 1 < nz:
                        Wbot_parts.append(ds.variables["W"][mask, k + 1].data)
                    else:
                        Wbot_parts.append(np.zeros_like(Wtop_parts[-1]))
                    t_parts.append(t[mask])
            if not T_parts:
                continue

            T_k = np.concatenate(T_parts, axis=0)
            S_k = np.concatenate(S_parts, axis=0)
            Wtop_k = np.concatenate(Wtop_parts, axis=0)
            Wbot_k = np.concatenate(Wbot_parts, axis=0)
            t_k = np.concatenate(t_parts)
            order = np.argsort(t_k)
            T_k, S_k, Wtop_k, Wbot_k, t_k = (
                T_k[order], S_k[order], Wtop_k[order], Wbot_k[order], t_k[order])

            dt_s = float(np.median(np.diff(t_k)))
            rho_k = eos(T_k, S_k, tRef, sRef, tAlpha, sBeta, rhoNil)
            rho_prime_k = rho_k - rho_mean_p[k][None]
            w_c_k = 0.5 * (Wtop_k + Wbot_k)

            rho_bp = bandpass_time(rho_prime_k, dt_s, band_cph)
            w_bp = bandpass_time(w_c_k, dt_s, band_cph)

            peso_k = drF[k] * hFacC_p[k]
            C[j0:j0 + tny, i0:i0 + tnx] += G * (rho_bp * w_bp).mean(axis=0) * peso_k
            n_usado = len(t_k)

        if verbose:
            print(f"    tile {gi + 1}/{len(grupos)} (j0={j0},i0={i0}): "
                  f"{nz} niveles, nt={n_usado} en banda {band_cph} cph")

    return C, n_usado


# ══════════════════════════════════════════════════════════════════════
# 5. Bucle principal
# ══════════════════════════════════════════════════════════════════════
def procesar_scen_run(scen, run, band_seiche=False, verbose=True):
    p = SCENARIOS[scen]
    run_dir = p["base"] / run
    run_label = "bay" if run == "run_expand" else "nobay"

    if verbose:
        print(f"\n{'=' * 60}\nEscenario: {scen}  |  Run: {run_label}\n  Ruta: {run_dir}")

    stages = list(discover_stage_patterns(run_dir))
    if verbose:
        print(f"  {len(stages)} chunk(s) de salida descubiertos")

    drF, hFacC, Depth, rAc, Y, X = leer_grid_global(run_dir)
    ny, nx = len(Y), len(X)

    if verbose:
        print("  Paso 1: ρ̄(x,y,z) — media temporal por columna sobre todo el run ...")
    rho_mean = calcular_rho_media(
        run_dir, stages, drF, ny, nx,
        p["tRef"], p["sRef"], p["tAlpha"], p["sBeta"], p["rhoNil"], verbose=verbose,
    )

    if verbose:
        print("  Paso 2: C(x,y,t) = g·∫ρ'·w dz, acumulando por ventana temporal ...")
    C_mean, N_win, t_seconds, Cdom_t = procesar_run(
        run_dir, stages, drF, hFacC, rAc, ny, nx, rho_mean,
        p["tRef"], p["sRef"], p["tAlpha"], p["sBeta"], p["rhoNil"], verbose=verbose,
    )

    if verbose:
        for label, _, _ in WINDOWS:
            n = N_win[label]
            if n == 0:
                print(f"    [{label}] sin datos")
                continue
            C = C_mean[label]
            print(f"    [{label}] N={n}  <C>_dominio={np.nanmean(C):.4e} W/m²  "
                  f"max|C|={np.nanmax(np.abs(C)):.4e} W/m²")

    C_seiche, N_seiche = None, 0
    if band_seiche:
        if verbose:
            print(f"  Paso 3: C_seiche(x,y) — filtrado banda {SEICHE_BAND_CPH} cph "
                  f"sobre t≥{T_PULSO_FIN:.0f}s (nivel por nivel) ...")
        C_seiche, N_seiche = calcular_C_seiche(
            run_dir, stages, drF, hFacC, ny, nx, rho_mean,
            p["tRef"], p["sRef"], p["tAlpha"], p["sBeta"], p["rhoNil"],
            T_PULSO_FIN, verbose=verbose,
        )
        if verbose and N_seiche > 0:
            total_seiche_W = float(np.nansum(C_seiche * rAc))
            total_banda_ancha_W = float(np.mean(Cdom_t[t_seconds >= T_PULSO_FIN]))
            print(f"    N={N_seiche}  <C_seiche>={np.nanmean(C_seiche):.4e} W/m²  "
                  f"max|C_seiche|={np.nanmax(np.abs(C_seiche)):.4e} W/m²")
            print(f"    total banda-seiche = {total_seiche_W:.4e} W  |  "
                  f"total banda-ancha (t≥{T_PULSO_FIN:.0f}s) = {total_banda_ancha_W:.4e} W")

    return dict(
        scen=scen, run=run_label, X=X, Y=Y, Depth=Depth, drF=drF,
        C_mean=C_mean, N_win=N_win, t_seconds=t_seconds, Cdom_t=Cdom_t,
        C_seiche=C_seiche, N_seiche=N_seiche,
    )


def guardar(res):
    SAVE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = SAVE_DIR / f"{res['scen']}_{res['run']}_conversion.npz"

    meta = dict(
        scenario=res["scen"], run=res["run"],
        method="C(x,y) = <int rho'*g*w dz>_t (Kang & Fringer 2012, buoyancy-flux form)",
        units_C="W/m^2", units_Cdom="W",
        rho_ref_scope="media temporal por columna, todo el run",
        windows={label: dict(t_ini=t0, t_fin=t1, n=res["N_win"][label])
                 for label, t0, t1 in WINDOWS},
        nota_forzamiento=("pulso gaussiano de viento cortado en t=9000s (2.5h), "
                           "no forzamiento sostenido toda la etapa 1"),
    )
    if res.get("N_seiche", 0) > 0:
        meta["banda_seiche_cph"] = list(SEICHE_BAND_CPH)
        meta["metodo_seiche"] = (
            "rho' y w filtrados en banda (FFT, bandpass_time) ANTES de "
            "multiplicar, sobre t>=9000s concatenado -- aisla la conversion "
            "coherente del seiche de transientes/armonicos de banda ancha")

    save_kwargs = dict(
        X=res["X"], Y=res["Y"], Depth=res["Depth"], drF=res["drF"],
        t_seconds=res["t_seconds"], Cdom_t=res["Cdom_t"],
        meta_json=json.dumps(meta),
    )
    for label, _, _ in WINDOWS:
        save_kwargs[f"C_{label}"] = res["C_mean"][label]
        save_kwargs[f"N_{label}"] = res["N_win"][label]
    if res.get("N_seiche", 0) > 0:
        save_kwargs["C_seiche"] = res["C_seiche"]
        save_kwargs["N_seiche"] = res["N_seiche"]

    np.savez(out_path, **save_kwargs)
    print(f"  [OK] {out_path.name}")
    return out_path


# ══════════════════════════════════════════════════════════════════════
# 6. Figuras (opcional, --plot)
# ══════════════════════════════════════════════════════════════════════
def graficar(npz_path):
    import matplotlib.pyplot as plt

    with np.load(npz_path, allow_pickle=True) as d:
        X, Y, Depth = d["X"], d["Y"], d["Depth"]
        meta = json.loads(str(d["meta_json"]))
        campos = {label: d[f"C_{label}"] for label, _, _ in WINDOWS
                  if int(d[f"N_{label}"]) > 0}
        if "C_seiche" in d.files and int(d["N_seiche"]) > 0:
            campos["seiche"] = d["C_seiche"]

    if not campos:
        print(f"  [SKIP plot] {npz_path.name}: sin ventanas con datos")
        return

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    scen, run = meta["scenario"], meta["run"]

    for label, C in campos.items():
        vmax = np.nanpercentile(np.abs(C), 99) or 1e-12
        fig, ax = plt.subplots(figsize=(7, 5))
        pm = ax.pcolormesh(X / 1e3, Y / 1e3, np.ma.masked_where(Depth <= 0, C),
                           cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
        ax.contour(X / 1e3, Y / 1e3, Depth, levels=[0], colors="k", linewidths=0.8)
        ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]")
        ax.set_title(f"Conversión barotrópica→baroclínica  C(x,y)  [{scen} / {run} / {label}]")
        fig.colorbar(pm, ax=ax, label=r"$C$ [W m$^{-2}$]")
        fig.tight_layout()
        out_png = FIG_DIR / f"{scen}_{run}_{label}.png"
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print(f"    [fig] {out_png.name}")


# ══════════════════════════════════════════════════════════════════════
# 6b. Diferencia bay − nobay (aísla el efecto de la bahía)
# ══════════════════════════════════════════════════════════════════════
def graficar_diferencia_bay_nobay(scen, verbose=True):
    """C_bay(x,y) − C_nobay(x,y) por ventana (y banda-seiche si está
    disponible en ambos). Válido porque bay y nobay comparten exactamente
    la misma grilla X,Y (nobay solo rellena la bahía con tierra) -- resta
    directa punto a punto, enmascarada donde nobay es tierra (ahí no hay
    contraparte física: es literalmente donde está la bahía)."""
    import matplotlib.pyplot as plt

    f_bay = SAVE_DIR / f"{scen}_bay_conversion.npz"
    f_nobay = SAVE_DIR / f"{scen}_nobay_conversion.npz"
    if not f_bay.exists() or not f_nobay.exists():
        print(f"  [SKIP diff] {scen}: falta {f_bay.name} o {f_nobay.name} "
              "(correr ambos runs primero)")
        return None

    with np.load(f_bay, allow_pickle=True) as db, np.load(f_nobay, allow_pickle=True) as dn:
        if not (np.array_equal(db["X"], dn["X"]) and np.array_equal(db["Y"], dn["Y"])):
            raise ValueError(f"{scen}: las grillas bay/nobay no coinciden -- no se puede restar")
        X, Y = db["X"], db["Y"]
        Depth_bay, Depth_nobay = db["Depth"], dn["Depth"]
        drF = db["drF"]

        etiquetas = [label for label, _, _ in WINDOWS if int(db[f"N_{label}"]) > 0
                    and int(dn[f"N_{label}"]) > 0]
        if "C_seiche" in db.files and "C_seiche" in dn.files \
                and int(db["N_seiche"]) > 0 and int(dn["N_seiche"]) > 0:
            etiquetas.append("seiche")

        diffs = {}
        for label in etiquetas:
            key = f"C_{label}" if label != "seiche" else "C_seiche"
            diffs[label] = db[key] - dn[key]

    if not diffs:
        print(f"  [SKIP diff] {scen}: sin ventanas en común")
        return None

    # mascara comun: tierra en nobay = huella de la bahia (sin contraparte)
    comun_seco = Depth_nobay <= 0
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    run_dir_bay = SCENARIOS[scen]["base"] / "run_expand"
    _, _, _, rAc, _, _ = leer_grid_global(run_dir_bay)

    resultados = {}
    for label, Cd in diffs.items():
        Cd_masked = np.where(comun_seco, np.nan, Cd)
        total_W = float(np.nansum(np.where(comun_seco, 0.0, Cd) * rAc))
        resultados[label] = total_W
        if verbose:
            print(f"    [{scen}/{label}] efecto neto de la bahía (fuera de la bahía) "
                  f"= {total_W:.4e} W")

        vmax = np.nanpercentile(np.abs(Cd_masked), 99) or 1e-12
        fig, ax = plt.subplots(figsize=(7, 5))
        pm = ax.pcolormesh(X / 1e3, Y / 1e3, np.ma.masked_invalid(Cd_masked),
                           cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
        ax.contour(X / 1e3, Y / 1e3, Depth_bay, levels=[0], colors="k", linewidths=0.8)
        ax.contourf(X / 1e3, Y / 1e3, comun_seco, levels=[0.5, 1.5],
                   colors="0.75", alpha=0.6)
        ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]")
        ax.set_title(f"C_bay − C_nobay  [{scen} / {label}]")
        fig.colorbar(pm, ax=ax, label=r"$\Delta C$ [W m$^{-2}$]")
        fig.text(0.02, 0.02, "gris = huella de la bahía (sin contraparte en nobay)",
                 fontsize=8, color="0.4")
        fig.tight_layout()
        out_png = FIG_DIR / f"{scen}_diff_bay_menos_nobay_{label}.png"
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        if verbose:
            print(f"    [fig] {out_png.name}")

    return resultados


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scen", nargs="+", default=["2layer", "linear", "realistic"],
                    choices=["2layer", "linear", "realistic"])
    ap.add_argument("--run", nargs="+", default=RUNS, choices=RUNS)
    ap.add_argument("--plot", action="store_true", help="genera figuras PNG tras guardar")
    ap.add_argument("--band-seiche", action="store_true",
                    help="ademas calcula C_seiche(x,y) filtrando en banda antes de "
                         "multiplicar (mucho mas lento: lee cada archivo ~nz veces)")
    ap.add_argument("--diff", action="store_true",
                    help="al final, genera mapas y totales C_bay-C_nobay por escenario "
                         "a partir de los npz ya guardados (rapido, no recalcula)")
    ap.add_argument("--diff-only", action="store_true",
                    help="salta el procesamiento y solo genera las diferencias bay-nobay "
                         "a partir de npz existentes")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    if not args.diff_only:
        for scen in args.scen:
            base = SCENARIOS[scen]["base"]
            if not base.exists():
                print(f"[SKIP] {scen}: ruta no disponible ({base})\n")
                continue
            for run in args.run:
                try:
                    res = procesar_scen_run(scen, run, band_seiche=args.band_seiche,
                                            verbose=not args.quiet)
                except (FileNotFoundError, ValueError) as e:
                    print(f"  [SKIP] {scen}/{run}: {e}")
                    continue
                out_path = guardar(res)
                if args.plot:
                    graficar(out_path)

    if args.diff or args.diff_only:
        print("\n--- Diferencia bay - nobay (efecto neto de la bahia) ---")
        for scen in args.scen:
            print(f"\n  Escenario: {scen}")
            graficar_diferencia_bay_nobay(scen, verbose=not args.quiet)

    print("\n[Listo] conversion_baroclinica_2D.py completado.")


if __name__ == "__main__":
    main()
