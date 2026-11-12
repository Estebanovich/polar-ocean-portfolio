#!/usr/bin/env python3
"""
energia_viento.py
=================
Energía que el viento inyecta al dominio: la única fuente de energía del modelo
(no hay marea; el viento es un pulso de ~2.8 h).

    P(t) = ∬ τ_y(x, y, t) v_sup(x, y, t) dA      [W]     (τ_x = 0)
    E_viento = ∫ P dt                              [J]

τ_y: input/merid_local_bay_560x352.bin (>f8, 180 registros de 1200 s, igual en
los tres escenarios). MITgcm fija el registro k en t = (k + ½)·1200 s e
interpola linealmente (ciclo periódico de 216 000 s), así que τ_y sube hasta
0.099 N m⁻² a 2.5 h y vale 0 desde 2.83 h. τ_y y V están en los puntos V (cara
sur de la celda); v_sup = V del nivel 1 (la tensión entra como fuerza en la
capa de arriba) y dA = rAs·hFacS₁. V se lee cada 900 s (salida state) y P se
integra por trapecios hasta 3 h.

El parche de viento cubre la bahía: en la corrida sin bahía esa parte cae en
tierra, así que las dos corridas NO reciben la misma energía.

Salida: npz_energy/energia_viento.npz y tabla en pantalla.
Uso:    /opt/anaconda3/bin/python energia_viento.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
ARCHIVO_TAU = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/input/"
                   "merid_local_bay_560x352.bin")
NY, NX, DT_REG = 352, 560, 1200.0
T_FIN = 3 * 3600.0


def leer_tau(path=ARCHIVO_TAU):
    return np.fromfile(path, ">f8").reshape(-1, NY, NX)


def tau_en_t(R, t):
    """τ_y(t) interpolado como MITgcm: registro k en (k + ½)·DT_REG, periódico."""
    n = R.shape[0]
    x = t / DT_REG - 0.5
    k0 = int(np.floor(x))
    a = x - k0
    return (1 - a) * R[k0 % n] + a * R[(k0 + 1) % n]


def potencia(tau, v, dA):
    """P(t) = Σ τ v dA [W]; tau, v: (nt, ny, nx); dA: (ny, nx)."""
    return np.einsum("tyx,tyx,yx->t", tau, np.nan_to_num(v), dA)


def calcular(esc, run):
    from netCDF4 import Dataset
    from energy_flux_IWBay import MNC_GLOB, SCENARIOS, build_tile_map
    from MIT_DatProcesW import discover_stage_patterns
    rd = SCENARIOS[esc]["base"] / run
    tiles, _, _ = build_tile_map(rd, MNC_GLOB, "grid.t*.nc")
    dA = np.zeros((NY, NX))
    for tl in tiles:
        j0, i0, ny, nx = tl["j0"], tl["i0"], tl["ny"], tl["nx"]
        with Dataset(tl["file"]) as ds:
            hs = "HFacS" if "HFacS" in ds.variables else "hFacS"
            dA[j0:j0 + ny, i0:i0 + nx] = (ds.variables["rAs"][:ny].data
                                          * ds.variables[hs][0, :ny].data)
    sg = list(discover_stage_patterns(rd))[0]
    tiles, _, _ = build_tile_map(rd, MNC_GLOB, sg)
    with Dataset(tiles[0]["file"]) as ds:
        t = ds.variables["T"][:].data.copy()
    m = t <= T_FIN + 1
    V = np.zeros((m.sum(), NY, NX))
    for tl in tiles:
        j0, i0, ny, nx = tl["j0"], tl["i0"], tl["ny"], tl["nx"]
        with Dataset(tl["file"]) as ds:
            V[:, j0:j0 + ny, i0:i0 + nx] = ds.variables["V"][m, 0, :ny].data
    t = t[m]
    R = leer_tau()
    tau = np.stack([tau_en_t(R, tk) for tk in t])
    P = potencia(tau, V, dA)
    return dict(t=t, P=P, E=float(np.sum(0.5 * (P[1:] + P[:-1]) * np.diff(t))),
                area_forzada=float(dA[np.abs(R).max(0) > 0.01 * np.abs(R).max()].sum()))


def main():
    out = {}
    for esc in ["2layer", "linear", "realistic"]:
        for run in ["run_expand", "run_expand_nobay"]:
            r = calcular(esc, run)
            lab = "bay" if run == "run_expand" else "nobay"
            out[f"{esc}_{lab}"] = r
            print(f"{esc:9s} {lab:5s}  E_viento = {r['E'] / 1e12:6.3f} TJ   P máx = "
                  f"{r['P'].max() / 1e9:6.2f} GW   área de agua bajo el viento (>1 %) = "
                  f"{r['area_forzada'] / 1e6:8.0f} km²", flush=True)
    np.savez(BASE / "npz_energy" / "energia_viento.npz",
             **{f"{k}_{c}": v[c] for k, v in out.items() for c in ("t", "P")},
             **{f"{k}_E": v["E"] for k, v in out.items()},
             meta_json=json.dumps(dict(archivo_tau=str(ARCHIVO_TAU), t_fin_s=T_FIN)))


if __name__ == "__main__":
    main()
