#!/usr/bin/env python3
"""
fondo_kang_fringer.py
=====================
Estado de fondo único para todos los términos de energía (Kang y Fringer 2012):
ρ′ = ρ − ρ_b(z), con ρ_b el perfil inicial. Los experimentos IWBay parten del
reposo con un perfil horizontalmente uniforme (desvío < 1e-12 kg m⁻³), así que
ρ_b(z) es exactamente el estado de reposo, y N_b² = −(g/ρ₀) dρ_b/dz.

Si cada término usa otra referencia (media temporal, media de la región...) el
balance ΔE + ∇·F = C − D no cierra, porque ρ′ no es el mismo en C, F y APE.
"""
from __future__ import annotations

import numpy as np

G = 9.81


def perfil_fondo(R, mojado, tol=1e-6):
    """ρ_b(z) [kg m⁻³] = media de las celdas mojadas de cada nivel.
    R, mojado: (nz, ...). Error si algún nivel varía más de tol en la horizontal."""
    R = np.where(mojado, np.asarray(R, float), np.nan).reshape(R.shape[0], -1)
    rb = np.nanmean(R, axis=1)
    desvio = np.nanmax(np.abs(R - rb[:, None]), axis=1)
    if np.nanmax(desvio) > tol:
        raise ValueError(f"estado inicial no uniforme: desvío {np.nanmax(desvio):.3g} kg m⁻³")
    return rb


def n2_fondo(rho_b, RC, rho0):
    """N_b² [s⁻²] en los centros: −(g/ρ₀) dρ_b/dz, con z = RC (negativo hacia abajo)."""
    return -G / rho0 * np.gradient(np.asarray(rho_b, float), np.asarray(RC, float))


def leer_fondo_inicial(run_dir, p):
    """ρ_b(z) del primer registro (t = 0) de la corrida; p = SCENARIOS[esc]."""
    from netCDF4 import Dataset
    from energy_flux_IWBay import MNC_GLOB, build_tile_map, eos
    from MIT_DatProcesW import discover_stage_patterns
    sg = list(discover_stage_patterns(run_dir))[0]
    tiles, _, _ = build_tile_map(run_dir, MNC_GLOB, sg)
    R, M = [], []
    for tl in tiles:
        with Dataset(tl["file"]) as ds:
            if float(ds.variables["T"][0]) != 0.0:
                raise ValueError(f"el primer registro de {sg} no es t = 0")
            T = ds.variables["Temp"][0].data
            S = ds.variables["S"][0].data
        R.append(eos(T, S, p["tRef"], p["sRef"], p["tAlpha"], p["sBeta"],
                     p["rhoNil"]).reshape(T.shape[0], -1))
        M.append((T != 0).reshape(T.shape[0], -1))
    return perfil_fondo(np.concatenate(R, 1), np.concatenate(M, 1))
