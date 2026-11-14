#!/usr/bin/env python3
"""
ape_dos_capas.py
================
Energía potencial disponible del escenario de dos capas con el desplazamiento
de la interfaz:

    APE(t) = ½ ρ₀ g′ ∬ ζ_i(x, y, t)² dA,     g′ = g Δρ/ρ₀,

que es ½ρ₀∫N²ξ² dz (Siyanbola 2023, ec. 7) con N² concentrado en la interfaz.
La fórmula ½g²ρ′²/(ρ₀N²) de energy_calcs_IWBay.py no sirve aquí: con N² ≈ 0
dentro de las capas se dispara (1.3e12 J en la bahía).

z_i(x, y, t) = −h₁ con h₁ = Σ_k (ρ₂ − ρ_k)/(ρ₂ − ρ₁) Δz_k el espesor equivalente
de la capa superior (ρ₁, ρ₂ = capa de arriba y de abajo del estado inicial, Δz =
drF·hFacC; ρ de la EOS lineal con sBeta = 0), y ζ_i = z_i − z̄_i, con z̄_i su media
temporal por columna en todo el run. h₁ no cambia con la difusión vertical y un
desplazamiento advectivo δ lo cambia exactamente en δ. Hasta el 18-oct-2026 z_i era
el primer cruce de ½(ρ₁ + ρ₂) interpolado entre centros de celda, que con la
interfaz repartida en 1–2 celdas da ~½ de la amplitud (correlación 1.000 con h₁).
Regiones: bay (j 324–347, i 274–285) y deep_full (j 0–322), como en
energy_calcs_IWBay.py.

Salida: npz_energy/2layer_ape_interfaz.npz
Uso:    /opt/anaconda3/bin/python ape_dos_capas.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
G = 9.81


def profundidad_interfaz(rho, dz, rho1, rho2):
    """z_i = −h₁ (negativo), h₁ = Σ (ρ₂ − ρ)/(ρ₂ − ρ₁) dz el espesor equivalente
    de la capa superior; rho: (nt, nz, ...), dz: (nz, ...) = drF·hFacC (0 en
    tierra). NaN en columnas secas."""
    rho = np.asarray(rho, float)
    dz = np.asarray(dz, float)
    frac = np.where(dz[None] > 0, (rho2 - rho) / (rho2 - rho1), 0.0)
    h1 = (frac * dz[None]).sum(axis=1)
    return np.where(dz.sum(0)[None] > 0, -h1, np.nan)


def profundidad_interfaz_cruce(rho, RC, rho_obj):
    """(Método anterior, sólo para comparar.) z (negativo) del primer cruce de
    rho_obj bajando; rho: (nt, nz, ...)."""
    rho = np.asarray(rho, float)
    debajo = rho >= rho_obj
    k = np.argmax(debajo, axis=1)                       # primer nivel ≥ ρ_obj
    hay = debajo.any(axis=1) & (k > 0)
    k = np.clip(k, 1, rho.shape[1] - 1)
    r1 = np.take_along_axis(rho, k[:, None], 1)[:, 0]
    r0 = np.take_along_axis(rho, (k - 1)[:, None], 1)[:, 0]
    z1, z0 = RC[k], RC[k - 1]
    a = np.where(r1 != r0, (rho_obj - r0) / np.where(r1 != r0, r1 - r0, 1.0), 0.0)
    return np.where(hay, z0 + a * (z1 - z0), np.nan)


def ape(zeta, dA, g_prima, rho0):
    """½ρ₀g′Σζ²dA [J]; zeta: (nt, ...), dA con la forma espacial."""
    return 0.5 * rho0 * g_prima * np.nansum(np.asarray(zeta) ** 2 * dA, axis=1)


def interfaz_run(run_dir, p):
    """z_i(t, y, x) [m] para toda la corrida, y los tiempos [s]."""
    from netCDF4 import Dataset
    from energy_flux_IWBay import MNC_GLOB, build_tile_map, eos
    from MIT_DatProcesW import discover_stage_patterns
    from conversion_baroclinica_2D import leer_grid_global
    drF, hFacC, Depth, rAc, Y, X = leer_grid_global(run_dir)
    stages = list(discover_stage_patterns(run_dir))
    # ρ_obj: punto medio entre la capa de arriba y la de abajo del estado inicial
    tiles, _, _ = build_tile_map(run_dir, MNC_GLOB, stages[0])
    with Dataset(tiles[0]["file"]) as ds:
        T0 = ds.variables["Temp"][0, :, 0, 0].data
    r0 = eos(T0, 0 * T0, p["tRef"], p["sRef"], p["tAlpha"], p["sBeta"], p["rhoNil"])
    rho_obj = 0.5 * (r0[0] + r0[-1])
    dz = drF[:, None, None] * hFacC
    g_prima = G * (r0[-1] - r0[0]) / p["rhoNil"]
    Z, T = [], []
    for sg in stages:
        tiles, _, _ = build_tile_map(run_dir, MNC_GLOB, sg)
        with Dataset(tiles[0]["file"]) as ds:
            t = ds.variables["T"][:].data.copy()
        z = np.full((t.size, len(Y), len(X)), np.nan, np.float32)
        for tl in tiles:
            j0, i0, ny, nx = tl["j0"], tl["i0"], tl["ny"], tl["nx"]
            with Dataset(tl["file"]) as ds:
                Tt = ds.variables["Temp"][:].data
            rho = eos(Tt, 0 * Tt, p["tRef"], p["sRef"], p["tAlpha"], p["sBeta"], p["rhoNil"])
            z[:, j0:j0 + ny, i0:i0 + nx] = profundidad_interfaz(
                rho, dz[:, j0:j0 + ny, i0:i0 + nx], r0[0], r0[-1])
        Z.append(z); T.append(t)
        print(f"    {sg}: {t.size} pasos", flush=True)
    t = np.concatenate(T)
    o = np.argsort(t)
    _, u = np.unique(t[o], return_index=True)
    return np.concatenate(Z)[o][u], t[o][u], rAc, g_prima, rho_obj


def main():
    from energy_flux_IWBay import SCENARIOS
    p = SCENARIOS["2layer"]
    reg = {"bay": (slice(324, 348), slice(274, 286)), "deep_full": (slice(0, 323), slice(None))}
    out = {}
    for run, lab in (("run_expand", "bay"), ("run_expand_nobay", "nobay")):
        print(f"\n=== 2layer / {lab} ===", flush=True)
        z, t, rA, gp, robj = interfaz_run(p["base"] / run, p)
        zeta = z - np.nanmean(z, axis=0)[None]
        for r, (sj, si) in reg.items():
            out[f"ape_{r}_{lab}"] = ape(zeta[:, sj, si].reshape(t.size, -1),
                                        rA[sj, si].ravel(), gp, p["rhoNil"])
            print(f"  {r:9s}: APE máx {out[f'ape_{r}_{lab}'].max():.3e} J", flush=True)
        out["t_horas"] = t / 3600.0
    np.savez(BASE / "npz_energy" / "2layer_ape_interfaz.npz", **out,
             meta_json=json.dumps(dict(g_prima=gp, rho_obj=robj, metodo="½ρ₀g′∬ζ_i²dA, ζ_i de h₁ (espesor equivalente)")))
    print("[Listo] APE de dos capas", flush=True)


if __name__ == "__main__":
    main()
