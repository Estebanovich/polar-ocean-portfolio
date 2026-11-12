#!/usr/bin/env python3
"""
disipacion_viscosa_bahia.py
===========================
Disipación viscosa explícita en la bahía (filas j 324–347, columnas i 274–285,
corrida con bahía), para compararla con el residuo D del balance de Kang y
Fringer (balance_kang_fringer.py). Viscosidad laplaciana del modelo
(viscAh = 100 m² s⁻¹, viscAr = 1e-5 m² s⁻¹, paredes y fondo sin fricción:
no_slip_sides = no_slip_bottom = .FALSE., sin arrastre de fondo):

    ε = ρ₀ [ν_h (|∇_h u|² + |∇_h v|²) + ν_v (|∂_z u|² + |∂_z v|²)]   [W m⁻³]

en diferencias de la malla C, como las arma el modelo:
    ∂u/∂x, ∂v/∂y en los centros (dxF, dyF);  ∂u/∂y, ∂v/∂x en las esquinas
    (dyU, dxV), sólo entre dos caras mojadas (pared sin fricción: sin esfuerzo
    tangencial);  ∂_z en las interfaces entre niveles (drC), sólo entre dos
    caras mojadas (fondo sin fricción).

Separación como en el balance: u = ū + u′ (ū media vertical por columna).
    ε_bt = ρ₀ ν_h Σ H |∇_h ū|² dA            (parte barotrópica)
    ε_bc = ρ₀ ν_h Σ |∇_h u′|² dV + ρ₀ ν_v Σ |∂_z u′|² dV
    (el término cruzado ∇ū·∇u′ se anula al integrar en la columna con H
    uniforme; se reporta aparte como comprobación)
    ε_sup: la parte de ε_bc en la celda superficial (1 m)

Salida: npz_energy/{esc}_disipacion_viscosa_bahia.npz (series de potencia [W])
Uso:    /opt/anaconda3/bin/python disipacion_viscosa_bahia.py [--scenarios ...]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
J_BAHIA = (324, 348)
I_BAHIA = (274, 286)
NU_H, NU_V = 100.0, 1e-5


def _media(a, w):
    """media vertical ponderada; a, w: (nz, ...)."""
    W = w.sum(0)
    return (a * w).sum(0) / np.where(W > 0, W, np.inf)


def disipacion(U, V, g, rho0):
    """(ε_h_bt, ε_h_bc, ε_v_bc, cruzado, ε_sup) [W] de un instante.
    U: (nz, ny, nx+1) caras oeste; V: (nz, ny+1, nx) caras sur."""
    mW, mS = g["hW"] > 0, g["hS"] > 0
    wW, wS = g["hW"] * g["drF"][:, None, None], g["hS"] * g["drF"][:, None, None]
    Ub, Vb = _media(U, wW), _media(V, wS)
    Up, Vp = np.where(mW, U - Ub, 0.0), np.where(mS, V - Vb, 0.0)
    volC = g["hC"] * g["drF"][:, None, None] * g["rA"]

    def grad_h(u, v, vol, por_columna=False):
        """Σ (|∇u|² + |∇v|²) · peso, en centros y esquinas."""
        dux = (u[..., 1:] - u[..., :-1]) / g["dxF"]
        dvy = (v[..., 1:, :] - v[..., :-1, :]) / g["dyF"]
        s = ((dux ** 2 + dvy ** 2) * vol).sum()
        # esquinas interiores: ∂u/∂y entre filas, ∂v/∂x entre columnas
        mu = mW[..., 1:, :] & mW[..., :-1, :] if not por_columna else (g["HW"][1:] > 0) & (g["HW"][:-1] > 0)
        mv = mS[..., :, 1:] & mS[..., :, :-1] if not por_columna else (g["HS"][:, 1:] > 0) & (g["HS"][:, :-1] > 0)
        duy = (u[..., 1:, :] - u[..., :-1, :]) / g["dyU"]
        dvx = (v[..., :, 1:] - v[..., :, :-1]) / g["dxV"]
        if por_columna:
            Hz_u = np.minimum(g["HW"][1:], g["HW"][:-1])
            Hz_v = np.minimum(g["HS"][:, 1:], g["HS"][:, :-1])
            s += (np.where(mu, duy ** 2, 0) * g["rAzU"] * Hz_u).sum()
            s += (np.where(mv, dvx ** 2, 0) * g["rAzV"] * Hz_v).sum()
        else:
            dz = g["drF"][:, None, None]
            s += (np.where(mu, duy ** 2, 0) * g["rAzU"] * dz).sum()
            s += (np.where(mv, dvx ** 2, 0) * g["rAzV"] * dz).sum()
        return s

    bt = rho0 * NU_H * grad_h(Ub, Vb, g["HC"] * g["rA"], por_columna=True)
    bc_h = rho0 * NU_H * grad_h(Up, Vp, volC)
    tot_h = rho0 * NU_H * grad_h(np.where(mW, U, 0.0), np.where(mS, V, 0.0), volC)

    def grad_v(u, m, rAf):
        ok = m[1:] & m[:-1]
        du = (u[1:] - u[:-1]) / g["drC"][1:-1, None, None]
        e = np.where(ok, du ** 2, 0.0) * g["drC"][1:-1, None, None] * rAf
        return e.sum(), e[0].sum()             # total y la interfaz 1 m

    vu, vu0 = grad_v(Up, mW, g["rAw"])
    vv, vv0 = grad_v(Vp, mS, g["rAs"])
    bc_v = rho0 * NU_V * (vu + vv)
    # parte superficial: horizontal de la celda 1 + vertical de su interfaz inferior
    sup_h = rho0 * NU_H * grad_h(np.where(np.arange(U.shape[0])[:, None, None] == 0, Up, 0.0),
                                 np.where(np.arange(V.shape[0])[:, None, None] == 0, Vp, 0.0),
                                 volC)
    sup = sup_h + rho0 * NU_V * (vu0 + vv0)
    return bt, bc_h, bc_v, tot_h - bt - bc_h, sup


def _malla(run_dir):
    from energy_flux_IWBay import MNC_GLOB, build_tile_map
    from balance_kang_fringer import _recorte
    from netCDF4 import Dataset
    tiles, _, _ = build_tile_map(run_dir, MNC_GLOB, "grid.t*.nc")
    jr, ir = J_BAHIA, I_BAHIA
    jV, iU = (jr[0], jr[1] + 1), (ir[0], ir[1] + 1)
    R = lambda v, j, i, t="C": _recorte(tiles, v, None, j, i, t)
    with Dataset(tiles[0]["file"]) as ds:
        drF = ds.variables["drF"][:].data.astype(float)
        drC = ds.variables["drC"][:].data.astype(float)
        hf = "HFac" if "HFacC" in ds.variables else "hFac"   # 2layer y realista: hFac, rAc
        ra = "rA" if "rA" in ds.variables else "rAc"
    g = dict(drF=drF, drC=drC,
             hC=R(hf + "C", jr, ir), hW=R(hf + "W", jr, iU, "U"), hS=R(hf + "S", jV, ir, "V"),
             rA=R(ra, jr, ir), dxF=R("dxF", jr, ir), dyF=R("dyF", jr, ir),
             rAw=R("rAw", jr, iU, "U"), rAs=R("rAs", jV, ir, "V"))
    # esquinas: ∂u/∂y entre filas j−1, j de U (esquinas interiores en y)
    dyU = R("dyU", (jr[0] + 1, jr[1]), iU, "U")       # (Yp1, Xp1) → filas interiores
    dxV = R("dxV", jV, (ir[0] + 1, ir[1]), "U")
    rAz_u = R("rAz", (jr[0] + 1, jr[1]), iU, "U")
    rAz_v = R("rAz", jV, (ir[0] + 1, ir[1]), "U")
    g.update(dyU=dyU, dxV=dxV, rAzU=rAz_u, rAzV=rAz_v)
    dz = drF[:, None, None]
    g["HC"] = (g["hC"] * dz).sum(0)
    g["HW"] = (g["hW"] * dz).sum(0)
    g["HS"] = (g["hS"] * dz).sum(0)
    return g


def calcular(esc):
    from energy_flux_IWBay import MNC_GLOB, SCENARIOS, build_tile_map
    from MIT_DatProcesW import discover_stage_patterns
    from balance_kang_fringer import _recorte
    from netCDF4 import Dataset
    p = SCENARIOS[esc]
    run_dir = p["base"] / "run_expand"
    rho0 = p["rhoNil"]
    g = _malla(run_dir)
    jr, ir = J_BAHIA, I_BAHIA
    out = {k: [] for k in ("t", "bt", "bc_h", "bc_v", "cruzado", "sup")}
    for sg in discover_stage_patterns(run_dir):
        tiles, _, _ = build_tile_map(run_dir, MNC_GLOB, sg)
        with Dataset(tiles[0]["file"]) as ds:
            t = ds.variables["T"][:].data.copy()
        U = _recorte(tiles, "U", slice(None), jr, (ir[0], ir[1] + 1), "U")
        V = _recorte(tiles, "V", slice(None), (jr[0], jr[1] + 1), ir, "V")
        for n in range(t.size):
            r = disipacion(U[n], V[n], g, rho0)
            out["t"].append(t[n])
            for k, v in zip(("bt", "bc_h", "bc_v", "cruzado", "sup"), r):
                out[k].append(v)
    t = np.asarray(out["t"])
    o = np.argsort(t, kind="stable")
    o = o[np.concatenate(([True], np.diff(t[o]) > 1e-6))]
    return {k: np.asarray(v)[o] for k, v in out.items()}


def main():
    from figura_energia_balance_kf import _acumular, terminos
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenarios", nargs="+", default=["2layer", "linear", "realistic"])
    a = ap.parse_args()
    for esc in a.scenarios:
        r = calcular(esc)
        th = r["t"] / 3600.0
        out = BASE / "npz_energy" / f"{esc}_disipacion_viscosa_bahia.npz"
        np.savez_compressed(out, t_horas=th, **{k: v for k, v in r.items() if k != "t"},
                            meta_json=np.array(json.dumps(dict(
                                nu_h=NU_H, nu_v=NU_V, region="j 324–347, i 274–285",
                                unidades="W"))))
        b = np.load(BASE / "npz_energy" / f"{esc}_balance_kf_bahia.npz")
        T = terminos(b)
        acc = {k: _acumular(th, r[k])[-1] / 1e9 for k in ("bt", "bc_h", "bc_v", "cruzado", "sup")}
        print(f"\n=== {esc} === (0–{th[-1]:.0f} h, GJ)")
        print(f"  ε_bt          = {acc['bt']:8.3f}   D_bt = {T['bt']['D'][-1] / 1e9:8.3f}")
        print(f"  ε_bc horiz.   = {acc['bc_h']:8.3f}")
        print(f"  ε_bc vertical = {acc['bc_v']:8.3f}")
        print(f"  ε_bc total    = {acc['bc_h'] + acc['bc_v']:8.3f}   D_bc = {T['bc']['D'][-1] / 1e9:8.3f}")
        print(f"  (celda de 1 m = {acc['sup']:8.3f};  cruzado = {acc['cruzado']:.3g})")
        print(f"  [OK] {out}")


if __name__ == "__main__":
    main()
