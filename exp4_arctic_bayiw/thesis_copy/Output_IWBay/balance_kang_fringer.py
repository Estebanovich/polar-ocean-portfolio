#!/usr/bin/env python3
"""
balance_kang_fringer.py
=======================
Balance de energía barotrópica y baroclínica de la bahía (filas j 324–347,
columnas i 274–285, corrida con bahía), con las ecuaciones de Kang y Fringer
(2012, ecs. 13–17) en forma hidrostática (q = 0), de 0 a 120 h:

    ∂E_bt/∂t + ∇·(F_bt + A_bt) = P_bt − C − D_bt
    ∂E_bc/∂t + ∇·(F_bc + A_bc) = P_bc + C − D_bc

Descomposición: u = ū + u′ (ū media vertical de la columna, instantánea) y
ρ′ = ρ − ρ_b(z), con ρ_b el perfil inicial de reposo (fondo_kang_fringer.py),
el mismo en todos los términos.

    E_bt = ½ρ₀ Σ H|ū|² dA + ½ρ₀ g Σ η² dA
    E_bc = ½ρ₀ Σ |u′|² dV + Σ E′_p dV,  E′_p la APE local exacta (ec. 12):
           g∫_{z_r}^{z}[ρ − ρ_b(z′)]dz′, la misma en los tres escenarios (en dos
           capas equivale a ½ρ₀g′ζ_i²; con N_b constante, a ½g²ρ′²/(ρ₀N_b²))
    p̃    = ρ₀ g η + g ∫_z^0 ρ′ dz′,  p_bt = media vertical de p̃,  p′ = p̃ − p_bt
    F_bt = Σ_boca H v̄ p_bt dx,  F_bc = Σ_boca Σ_k v′ p′ dz dx   (trabajo de
           presión; cara sur de la fila 324, abierta en los 163.3 m de arriba;
           + = sale de la bahía)
    A_bt = Σ_boca V H E_k0 dx,  A_bc = Σ_boca Σ_k v (E′_k + E′_k0 + E′_p) dz dx
           (advección de energía, ecs. 15–16); el flujo total es F + A. A_ape es
           la parte de E′_p en A_bc; en dos capas, A_ape_interfaz es la misma con
           la APE de la interfaz (ver ape_interfaz). Los
           términos de difusión de las ecs. 15–16 (−ν_H∇E) quedan en D.
    C    = Σ ρ′ g W_bt dV,  W_bt = −∇·((H + z) ū) discreta (conversion_bt_2D.py)
    C_escalon: la misma C en la fila 323 frente a la boca (i 274–285), FUERA de
           la caja: el escalón 600 → 163 m, donde se concentra la conversión del
           dominio. Las ondas generadas ahí no cruzan la cara de la boca.
    P_bt = Σ τ_y v̄ dA,  P_bc = Σ τ_y (v_sup − v̄) dA   (τ en los puntos V; la
           cara de la boca cuenta con peso ½; τ = 0 en la etapa 2)
    D    residuo: disipación (ecs. 18–23), difusión de energía y error numérico.

El modelo es hidrostático: la KE no lleva w². Las velocidades se llevan al
centro de la celda promediando las dos caras.

Salida: npz_energy/{esc}_balance_kf_bahia.npz (series en el tiempo)
Uso:    /opt/anaconda3/bin/python balance_kang_fringer.py [--scenarios ...]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
G = 9.81
J_BAHIA = (324, 348)            # filas de la bahía (la 323 es océano)
I_BAHIA = (274, 286)
J_BOCA = 324                    # cara sur de la fila 324
T_ETAPA2 = 216000.0             # la etapa 2 corre sin archivo de viento


# ── piezas puras ──────────────────────────────────────────────────────────

def presion_tilde(rp, eta, dz, rho0):
    """p̃ en los centros: ρ₀gη + g(Σ_{m<k} ρ′dz + ½ρ′_k dz_k). rp, dz: (nz, ...)."""
    capa = rp * dz
    encima = np.cumsum(capa, axis=0) - capa
    return rho0 * G * eta[None] + G * (encima + 0.5 * capa)


def flujo_cara(v, p, dz):
    """(H v̄ p̄, Σ v′p′dz) por columna de una cara; v, p, dz: (nz, n)."""
    H = dz.sum(0)
    Hs = np.where(H > 0, H, 1.0)
    vb, pb = (v * dz).sum(0) / Hs, (p * dz).sum(0) / Hs
    return H * vb * pb, ((v - vb) * (p - pb) * dz).sum(0)


def ke_bt_bc(u, v, vol, rho0):
    """(½ρ₀ΣH|ū|²dA, ½ρ₀Σ|u′|²dV) con u, v, vol: (nz, ...) en los centros."""
    H = vol.sum(0)
    Hs = np.where(H > 0, H, 1.0)
    ub, vb = (u * vol).sum(0) / Hs, (v * vol).sum(0) / Hs
    bt = 0.5 * rho0 * ((ub ** 2 + vb ** 2) * H).sum()
    bc = 0.5 * rho0 * (((u - ub) ** 2 + (v - vb) ** 2) * vol).sum()
    return bt, bc


def pe_bt(eta, dA, rho0):
    return 0.5 * rho0 * G * float((eta ** 2 * dA).sum())


def ape_lineal(rp, n2, vol, rho0):
    """Σ ½g²ρ′²/(ρ₀N_b²) dV; n2 se difunde sobre rp (nz, ...)."""
    n2 = np.asarray(n2).reshape(np.shape(n2) + (1,) * (np.ndim(rp) - np.ndim(n2)))
    return float((0.5 * G ** 2 * rp ** 2 / (rho0 * n2) * vol).sum())


def tau_modelo(R, t):
    """τ_y que ve el modelo: interpolado como MITgcm en la etapa 1 (incluido el
    retorno hacia el registro 0 en 59.83–60 h) y 0 en la etapa 2."""
    from energia_viento import tau_en_t
    return np.zeros_like(R[0]) if t >= T_ETAPA2 else tau_en_t(R, t)


def ape_exacta(rho_t, z, rho_b, z_b):
    """APE local exacta [J m⁻³] de Kang y Fringer (2012, ec. 12; 2010):
    g∫_{z_r}^{z} [ρ_t − ρ_b(z′)] dz′, con z_r el nivel de reposo donde
    ρ_b(z_r) = ρ_t. ρ_b es lineal a trozos entre los centros z_b (negativos,
    de la superficie al fondo) y constante fuera de ellos; con ρ_b lineal da
    exactamente ½g²ρ′²/(ρ₀N²), y con un escalón, ½gΔρζ² por columna."""
    d_b = -np.asarray(z_b, float)
    rb = np.maximum.accumulate(np.asarray(rho_b, float))      # no decreciente con d
    pend = np.diff(rb) / np.diff(d_b)
    Rn = np.concatenate(([0.0], np.cumsum(0.5 * (rb[1:] + rb[:-1]) * np.diff(d_b))))

    def R(d):                                   # ∫_{d_b[0]}^{d} ρ_b, exacta
        k = np.clip(np.searchsorted(d_b, d, side="right") - 1, 0, d_b.size - 2)
        x = np.clip(d, d_b[0], d_b[-1]) - d_b[k]
        dentro = Rn[k] + rb[k] * x + 0.5 * pend[k] * x ** 2
        return (dentro + rb[0] * np.minimum(d - d_b[0], 0.0)
                + rb[-1] * np.maximum(d - d_b[-1], 0.0))

    rt = np.asarray(rho_t, float)
    d = -np.asarray(z, float) * np.ones_like(rt)
    k = np.clip(np.searchsorted(rb, rt, side="right") - 1, 0, rb.size - 2)
    sube = rb[k + 1] - rb[k]
    d_r = d_b[k] + np.where(sube > 0, (rt - rb[k]) / np.where(sube > 0, sube, 1.0), 0.0) \
        * (d_b[k + 1] - d_b[k])
    d_r = np.where(rt <= rb[0], d_b[0], np.where(rt >= rb[-1], d_b[-1], d_r))
    return G * (rt * (d_r - d) - (R(d_r) - R(d)))


def ape_interfaz(z_i, z_i0, dA, g_prima, rho0):
    """APE de un fluido de dos capas [J]: ½ρ₀g′Σζ_i²dA, ζ_i = z_i − z_i0 el
    desplazamiento de la interfaz respecto del estado inicial (Gill 1982, ec. 6.7.2,
    p. 140; su ec. 6.7.3 es el límite continuo ½ρ₀N²h²);
    columnas sin interfaz (NaN) no suman."""
    ze = np.nan_to_num(np.asarray(z_i, float) - np.asarray(z_i0, float))
    return 0.5 * rho0 * g_prima * float((ze ** 2 * dA).sum())


def adveccion_ape_interfaz(v, z, z_i, z_i0, g_prima, rho0):
    """Advección de la APE de la interfaz por una cara, por columna [W m⁻¹]:
    v(z_i)·½ρ₀g′ζ_i², con v del nivel cuyo centro está más cerca de la interfaz
    (el término Σ v E′_p dz de la ec. 16 con la APE concentrada en la interfaz);
    columnas sin interfaz dan 0. v: (nz, n), z: (nz,)."""
    z_i = np.asarray(z_i, float)
    hay = np.isfinite(z_i)
    k = np.argmin(np.abs(np.asarray(z, float)[:, None] - np.where(hay, z_i, 0.0)[None]), axis=0)
    vi = np.take_along_axis(np.asarray(v, float), k[None], 0)[0]
    ze = np.where(hay, z_i - np.asarray(z_i0, float), 0.0)
    return np.where(hay, vi * 0.5 * rho0 * g_prima * np.nan_to_num(ze) ** 2, 0.0)


def flujo_advectivo_cara(u, v, dz, ep, rho0):
    """Advección de energía por una cara normal a y (Kang y Fringer ecs. 15–16):
    (V H E_k0, Σ v (E′_k + E′_k0 + E′_p) dz) por columna; u, v, dz, ep: (nz, n)."""
    H = dz.sum(0)
    Hs = np.where(H > 0, H, 1.0)
    U, V = (u * dz).sum(0) / Hs, (v * dz).sum(0) / Hs
    up, vp = u - U, v - V
    F0 = V * H * 0.5 * rho0 * (U ** 2 + V ** 2)
    E1 = 0.5 * rho0 * (up ** 2 + vp ** 2) + rho0 * (U * up + V * vp) + ep
    return F0, (v * E1 * dz).sum(0)


def trabajo_viento(tau, v_sup, v_bar, dA):
    """(Σ τ v̄ dA, Σ τ (v_sup − v̄) dA) [W]."""
    return float((tau * v_bar * dA).sum()), float((tau * (v_sup - v_bar) * dA).sum())


# ── lectura ───────────────────────────────────────────────────────────────

def _recorte(tiles, var, m, jr, ir, tipo="C"):
    """var de los tiles en las filas jr y columnas ir globales (rangos semiabiertos);
    tipo 'U' (caras oeste, Xp1), 'V' (caras sur, Yp1) o 'C'; m = None para una
    variable sin tiempo. Devuelve ([nt,] [nz,] ny, nx)."""
    from netCDF4 import Dataset
    out = None
    for tl in tiles:
        j0, i0 = tl["j0"], tl["i0"]
        ny = tl["ny"] + (tipo == "V")
        nx = tl["nx"] + (tipo == "U")
        ja, jb = max(jr[0], j0), min(jr[1], j0 + ny)
        ia, ib = max(ir[0], i0), min(ir[1], i0 + nx)
        if ja >= jb or ia >= ib:
            continue
        with Dataset(tl["file"]) as ds:
            x = ds.variables[var]
            idx = (Ellipsis, slice(ja - j0, jb - j0), slice(ia - i0, ib - i0))
            a = x[(m,) + idx if m is not None else idx].data
        if out is None:
            out = np.zeros(a.shape[:-2] + (jr[1] - jr[0], ir[1] - ir[0]))
        out[..., ja - jr[0]:jb - jr[0], ia - ir[0]:ib - ir[0]] = a
    return out


def _malla(run_dir):
    from conversion_baroclinica_2D import leer_grid_global
    from conversion_bt_2D import leer_malla_caras
    from energy_flux_IWBay import MNC_GLOB, build_tile_map
    drF, hC, Depth, rA, Y, X = leer_grid_global(run_dir)
    hW, hS, dyG, dxG = leer_malla_caras(run_dir)
    tiles, _, _ = build_tile_map(run_dir, MNC_GLOB, "grid.t*.nc")
    rAs = _recorte(tiles, "rAs", None, (J_BOCA, J_BAHIA[1] + 1), I_BAHIA, "V")
    return drF, hC, rA, hW, hS, dyG, dxG, rAs


def calcular(esc):
    from energy_flux_IWBay import MNC_GLOB, SCENARIOS, build_tile_map, eos
    from MIT_DatProcesW import discover_stage_patterns
    from conversion_bt_2D import areas_bajo_cara, w_bt_centro
    from energia_viento import leer_tau
    from fondo_kang_fringer import leer_fondo_inicial, n2_fondo
    from ape_dos_capas import profundidad_interfaz
    p = SCENARIOS[esc]
    rho0 = p["rhoNil"]
    run_dir = p["base"] / "run_expand"
    drF, hC, rA, hW, hS, dyG, dxG, rAs = _malla(run_dir)
    (ja, jb), (ia, ib) = J_BAHIA, I_BAHIA
    RC = -(np.cumsum(drF) - 0.5 * drF)
    rho_b = leer_fondo_inicial(run_dir, p)
    dos_capas = esc == "2layer"
    n2 = n2_fondo(rho_b, RC, rho0)
    g_prima = G * (rho_b[-1] - rho_b[0]) / rho0
    z_i0 = None
    if not dos_capas and n2.min() <= 0:     # sólo para la APE lineal de comparación
        raise ValueError(f"{esc}: N_b² ≤ 0 en algún nivel; la APE lineal no está definida")
    # malla de la caja (bahía) y de la fila 323 para la presión de la boca
    dzC = drF[:, None, None] * hC[:, ja - 1:jb, ia:ib]          # filas 323…347
    vol = dzC[:, 1:] * rA[None, ja:jb, ia:ib]
    dzB = drF[:, None] * hS[:, J_BOCA, ia:ib]                   # cara de la boca
    # caras de la caja extendida a la fila 323 (el escalón de la boca, del lado del
    # océano): W_bt de la bahía no cambia, y en la fila 323 da la conversión del escalón
    dW = drF[:, None, None] * hW[:, ja - 1:jb, ia:ib + 1]
    dS = drF[:, None, None] * hS[:, ja - 1:jb + 1, ia:ib]
    A_W, A_S = areas_bajo_cara(drF, hW[:, ja - 1:jb, ia:ib + 1], hS[:, ja - 1:jb + 1, ia:ib],
                               dyG[ja - 1:jb, ia:ib + 1], dxG[ja - 1:jb + 1, ia:ib])
    rAE = rA[ja - 1:jb, ia:ib]
    vol_esc = dzC[:, 0] * rA[None, ja - 1, ia:ib]
    rAb = rA[ja:jb, ia:ib]
    peso_v = rAs * hS[0, ja:jb + 1, ia:ib]
    peso_v[0] *= 0.5                                            # cara de la boca
    tau_R = leer_tau()[:, ja:jb + 1, ia:ib]
    out = {k: [] for k in ("t", "ke_bt", "pe_bt", "ke_bc", "ape", "ape_lineal", "F_bt",
                           "F_bc", "A_bt", "A_bc", "C", "C_escalon", "P_bt", "P_bc",
                           "ape_interfaz", "A_ape", "A_ape_interfaz",
                           "F_bt_col", "F_bc_col", "A_bt_col", "A_bc_col",
                           "A_ape_col", "A_ape_interfaz_col")}
    RC3 = RC[:, None, None] * np.ones((1, jb - ja + 1, ib - ia))   # filas 323…347
    mojado = dzC > 0
    for sg in discover_stage_patterns(run_dir):
        tiles, _, _ = build_tile_map(run_dir, MNC_GLOB, sg)
        from netCDF4 import Dataset
        with Dataset(tiles[0]["file"]) as ds:
            t = ds.variables["T"][:].data.copy()
        m = slice(None)
        Tt = _recorte(tiles, "Temp", m, (ja - 1, jb), I_BAHIA)
        St = _recorte(tiles, "S", m, (ja - 1, jb), I_BAHIA)
        Et = _recorte(tiles, "Eta", m, (ja - 1, jb), I_BAHIA)
        Ut = _recorte(tiles, "U", m, (ja - 1, jb), (ia, ib + 1), "U")   # filas 323…347
        Vt = _recorte(tiles, "V", m, (ja - 1, jb + 1), I_BAHIA, "V")   # caras 323…348
        rho = eos(Tt, St, p["tRef"], p["sRef"], p["tAlpha"], p["sBeta"], rho0)
        rp = np.where(dzC[None] > 0, rho - rho_b[None, :, None, None], 0.0)
        # Ū y V̄ en las caras (para W_bt)
        HW, HS = dW.sum(0), dS.sum(0)
        Ub = np.einsum("tkyx,kyx->tyx", Ut, dW) / np.where(HW > 0, HW, np.inf)
        Vb = np.einsum("tkyx,kyx->tyx", Vt, dS) / np.where(HS > 0, HS, np.inf)
        if dos_capas:
            z_i = profundidad_interfaz(rho[:, :, 1:], dzC[:, 1:], rho_b[0], rho_b[-1])
            if z_i0 is None:
                if t[0] != 0.0:
                    raise ValueError("el primer chunk no empieza en t = 0")
                z_i0 = z_i[0]
        for n in range(t.size):
            uc = 0.5 * (Ut[n, :, 1:, :-1] + Ut[n, :, 1:, 1:])
            vc = 0.5 * (Vt[n, :, 1:-1] + Vt[n, :, 2:])
            kb, kc = ke_bt_bc(uc, vc, vol, rho0)
            pt = presion_tilde(rp[n], Et[n], dzC, rho0)
            p_boca = 0.5 * (pt[:, 0] + pt[:, 1])                 # filas 323 y 324
            fbt, fbc = flujo_cara(Vt[n, :, 1], p_boca, dzB)
            C = C_esc = 0.0
            for k in range(len(drF)):
                if not hC[k, ja - 1:jb, ia:ib].any():
                    continue
                Wc = w_bt_centro(Ub[n:n + 1], Vb[n:n + 1], A_W, A_S, rAE, k)[0]
                C += G * float((rp[n, k, 1:] * Wc[1:] * vol[k]).sum())
                C_esc += G * float((rp[n, k, 0] * Wc[0] * vol_esc[k]).sum())
            ep = np.where(mojado, ape_exacta(np.where(mojado, rho[n], rho_b[:, None, None]),
                                             RC3, rho_b, RC), 0.0)
            ape = float((ep[:, 1:] * vol).sum())
            ape_lin = np.nan if dos_capas else ape_lineal(rp[n, :, 1:], n2, vol, rho0)
            ape_i = ape_interfaz(z_i[n], z_i0, rAb, g_prima, rho0) if dos_capas else np.nan
            u_f = 0.25 * (Ut[n, :, 0, :-1] + Ut[n, :, 0, 1:] + Ut[n, :, 1, :-1] + Ut[n, :, 1, 1:])
            ep_f = 0.5 * (ep[:, 0] + ep[:, 1])
            a0, a1 = flujo_advectivo_cara(u_f, Vt[n, :, 1], dzB, ep_f, rho0)
            a_ape = (Vt[n, :, 1] * ep_f * dzB).sum(0)          # parte de E′_p en a1
            a_api = (adveccion_ape_interfaz(Vt[n, :, 1], RC, z_i[n][0], z_i0[0], g_prima,
                                            rho0) if dos_capas else np.full(ib - ia, np.nan))
            Pbt, Pbc = trabajo_viento(tau_modelo(tau_R, t[n]), Vt[n, 0, 1:], Vb[n, 1:], peso_v)
            for k, val in (("t", t[n]), ("ke_bt", kb), ("pe_bt", pe_bt(Et[n, 1:], rAb, rho0)),
                           ("ke_bc", kc), ("ape", ape), ("ape_lineal", ape_lin),
                           ("ape_interfaz", ape_i),
                           ("A_ape", -float((a_ape * dxG[J_BOCA, ia:ib]).sum())),
                           ("A_ape_interfaz", -float((a_api * dxG[J_BOCA, ia:ib]).sum())),
                           ("F_bt", -float((fbt * dxG[J_BOCA, ia:ib]).sum())),
                           ("F_bc", -float((fbc * dxG[J_BOCA, ia:ib]).sum())),
                           ("A_bt", -float((a0 * dxG[J_BOCA, ia:ib]).sum())),
                           ("A_bc", -float((a1 * dxG[J_BOCA, ia:ib]).sum())),
                           ("C", C), ("C_escalon", C_esc), ("P_bt", Pbt), ("P_bc", Pbc),
                           # por columna de la boca [W], mismo signo (+ sale de la bahía)
                           ("F_bt_col", -fbt * dxG[J_BOCA, ia:ib]),
                           ("F_bc_col", -fbc * dxG[J_BOCA, ia:ib]),
                           ("A_bt_col", -a0 * dxG[J_BOCA, ia:ib]),
                           ("A_bc_col", -a1 * dxG[J_BOCA, ia:ib]),
                           ("A_ape_col", -a_ape * dxG[J_BOCA, ia:ib]),
                           ("A_ape_interfaz_col", -a_api * dxG[J_BOCA, ia:ib])):
                out[k].append(val)
        print(f"    {sg}: {t.size} pasos", flush=True)
    t = np.array(out["t"])
    o = np.argsort(t, kind="stable")
    keep = np.concatenate(([True], np.diff(t[o]) > 1e-6))
    res = {k: np.array(v)[o][keep] for k, v in out.items()}
    res["t_horas"] = res.pop("t") / 3600.0
    return res, dict(rho_b=rho_b.tolist(), ape="exacta (Kang y Fringer 2012 ec. 12)",
                     ape_lineal="½g²ρ′²/(ρ₀N_b²), sólo como comparación (NaN en dos capas)",
                     ape_interfaz="½ρ₀g′Σζ_i²dA, sólo dos capas (NaN en los otros)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenarios", nargs="+", default=["linear", "realistic", "2layer"])
    a = ap.parse_args()
    for esc in a.scenarios:
        print(f"\n=== {esc} ===", flush=True)
        r, meta = calcular(esc)
        meta.update(marco="Kang y Fringer 2012, hidrostático y lineal",
                    caja=dict(j=list(J_BAHIA), i=list(I_BAHIA)), boca=J_BOCA)
        p = BASE / "npz_energy" / f"{esc}_balance_kf_bahia.npz"
        np.savez(p, **r, meta_json=json.dumps(meta))
        print(f"[OK] {p}", flush=True)
    print("\n[Listo] balance de Kang y Fringer")


if __name__ == "__main__":
    main()
