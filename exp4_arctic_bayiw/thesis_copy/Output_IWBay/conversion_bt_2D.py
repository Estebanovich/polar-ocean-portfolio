#!/usr/bin/env python3
"""
conversion_bt_2D.py
===================
Conversión barotrópica → baroclínica con la velocidad vertical BAROTRÓPICA
(Kang y Fringer 2012, ec. de energía baroclínica):

    C(x, y) = ⟨ Σ_k ρ′_k g W_bt,k drF_k hFacC_k ⟩_t          [W m⁻²]

conversion_baroclinica_2D.py usa la w total, lo que da el flujo de flotabilidad
total (incluye el intercambio KE ↔ APE dentro de la propia onda). Aquí W_bt es
la w que induce el flujo barotrópico sobre la topografía, W_bt = −∇·((H+z)Ū),
en su forma discreta de malla C: en la cara superior de la celda k,

    W_top,k = −[δx(Ū_W A_W,k) + δy(V̄_S A_S,k)] / rA,
    A_W,k = dyG Σ_{m≥k} drF_m hFacW_m      (ídem A_S con dxG, hFacS),

que es integrate_for_w de MITgcm con u sustituida por la Ū de la columna
(Ū_W = Σ u drF hFacW / Σ drF hFacW). En la superficie da −∇·(transporte)/rA
(= ∂η/∂t) y respeta el escalón de la boca: bajo los 163 m la cara de la boca
está cerrada y ahí el flujo que llega desde el océano sube.

ρ′ = ρ − ρ_b(z), con ρ_b el perfil inicial de reposo (fondo_kang_fringer.py),
la misma referencia que la APE y el flujo del balance de Kang y Fringer.
Sólo t ≥ 2.5 h (después del pulso de viento):

    sin filtrar, por ventana: C_total (2.5–120 h), C_temprana (2.5–27 h, régimen
              del seiche) y C_tardia (30–120 h, régimen cuasi-inercial);
    por banda (ρ′ y W_bt filtrados antes de multiplicar; W_bt es lineal en Ū,
              así que se filtra Ū una vez para todo el dominio): C_seiche
              (0.235–0.29 cph) y C_ni (0.035–0.065 cph, la de la lámina del
              flujo en el corte).

Salida: npz_conversion_bt/{escenario}_{bay|nobay}_conversion_bt.npz
Uso:    /opt/anaconda3/bin/python conversion_bt_2D.py [--scen ...] [--run ...]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

SAVE_DIR = Path(__file__).resolve().parent / "npz_conversion_bt"
RUNS = ["run_expand", "run_expand_nobay"]
T_PULSO_FIN = 9000.0
SEICHE_BAND_CPH = (0.235, 0.290)
# bandas: ρ′ y W_bt filtrados antes de multiplicar (W_bt es lineal en Ū, así que
# basta filtrar Ū); ventanas: producto sin filtrar promediado en ese intervalo
BANDAS = {"seiche": SEICHE_BAND_CPH, "ni": (0.035, 0.065)}
VENTANAS = {"total": (T_PULSO_FIN, np.inf), "temprana": (T_PULSO_FIN, 27 * 3600.0),
            "tardia": (30 * 3600.0, np.inf)}
G = 9.81


def areas_bajo_cara(drF, hFacW, hFacS, dyG, dxG):
    """A_W,k y A_S,k: área de cada cara por debajo de la cara superior del
    nivel k [m²]. hFacW: (nz, ny, nx+1); hFacS: (nz, ny+1, nx)."""
    dW = np.asarray(drF)[:, None, None] * hFacW
    dS = np.asarray(drF)[:, None, None] * hFacS
    A_W = dW[::-1].cumsum(0)[::-1] * dyG[None]
    A_S = dS[::-1].cumsum(0)[::-1] * dxG[None]
    return A_W, A_S


def w_bt_caras(Ub, Vb, A_W, A_S, rA):
    """W barotrópica en las caras superiores (nt, nz, ny, nx). Ub: (nt, ny,
    nx+1) Ū en caras oeste/este; Vb: (nt, ny+1, nx) V̄ en caras sur/norte."""
    fx = Ub[:, None] * A_W[None]
    fy = Vb[:, None] * A_S[None]
    div = (fx[..., 1:] - fx[..., :-1]) + (fy[..., 1:, :] - fy[..., :-1, :])
    return -div / rA[None, None]


def w_bt_centro(Ub, Vb, A_W, A_S, rA, k):
    """W_bt en el centro de las celdas del nivel k (nt, ny, nx):
    ½(cara k + cara k+1), con la cara del fondo = 0."""
    W = w_bt_caras(Ub, Vb, A_W[k:k + 2], A_S[k:k + 2], rA)
    return 0.5 * (W[:, 0] + (W[:, 1] if W.shape[1] > 1 else 0.0))


def _filtrar(A, dt, banda, bloque=32):
    """bandpass_time por bloques de filas (eje 1) para acotar la memoria."""
    from prototipo_hovmoller_cpcg import bandpass_time
    out = np.empty_like(A)
    for j in range(0, A.shape[1], bloque):
        out[:, j:j + bloque] = bandpass_time(A[:, j:j + bloque].astype(np.float64),
                                             dt, banda)
    return out


def leer_malla_caras(run_dir):
    """hFacW, hFacS, dxG, dyG globales (con la cara extra al este y al norte)."""
    from netCDF4 import Dataset
    from energy_flux_IWBay import MNC_GLOB, build_tile_map
    tiles, Yg, Xg = build_tile_map(run_dir, MNC_GLOB, "grid.t*.nc")
    ny, nx = len(Yg), len(Xg)
    with Dataset(tiles[0]["file"]) as ds:
        nz = ds.variables["drF"].shape[0]
    hW = np.zeros((nz, ny, nx + 1)); hS = np.zeros((nz, ny + 1, nx))
    dyG = np.zeros((ny, nx + 1)); dxG = np.zeros((ny + 1, nx))
    for t in tiles:
        j0, i0, tny, tnx = t["j0"], t["i0"], t["ny"], t["nx"]
        with Dataset(t["file"]) as ds:
            v = ds.variables
            hW[:, j0:j0 + tny, i0:i0 + tnx + 1] = v["HFacW" if "HFacW" in v else "hFacW"][:].data
            hS[:, j0:j0 + tny + 1, i0:i0 + tnx] = v["HFacS" if "HFacS" in v else "hFacS"][:].data
            dyG[j0:j0 + tny, i0:i0 + tnx + 1] = ds.variables["dyG"][:].data
            dxG[j0:j0 + tny + 1, i0:i0 + tnx] = ds.variables["dxG"][:].data
    return hW, hS, dyG, dxG


def series_ubar(run_dir, stages, drF, hW, hS, t_min):
    """Ū_W(t) y V̄_S(t) globales para t ≥ t_min (float32) y los tiempos."""
    from netCDF4 import Dataset
    from energy_flux_IWBay import MNC_GLOB, build_tile_map
    nz, ny, nxp = hW.shape
    nx = nxp - 1
    HW = np.einsum("k,kyx->yx", drF, hW)
    HS = np.einsum("k,kyx->yx", drF, hS)
    U_p, V_p, t_p = [], [], []
    for sg in stages:
        tiles, _, _ = build_tile_map(run_dir, MNC_GLOB, sg)
        if not tiles:
            continue
        with Dataset(tiles[0]["file"]) as ds:
            t = ds.variables["T"][:].data.copy()
        m = t >= t_min
        if not m.any():
            continue
        Ut = np.zeros((m.sum(), ny, nx + 1), np.float32)
        Vt = np.zeros((m.sum(), ny + 1, nx), np.float32)
        for tl in tiles:
            j0, i0, tny, tnx = tl["j0"], tl["i0"], tl["ny"], tl["nx"]
            with Dataset(tl["file"]) as ds:
                U = ds.variables["U"][m].data
                V = ds.variables["V"][m].data
            Ut[:, j0:j0 + tny, i0:i0 + tnx + 1] = np.einsum(
                "tkyx,k,kyx->tyx", U, drF, hW[:, j0:j0 + tny, i0:i0 + tnx + 1])
            Vt[:, j0:j0 + tny + 1, i0:i0 + tnx] = np.einsum(
                "tkyx,k,kyx->tyx", V, drF, hS[:, j0:j0 + tny + 1, i0:i0 + tnx])
        U_p.append(Ut / np.where(HW > 0, HW, np.inf)[None])
        V_p.append(Vt / np.where(HS > 0, HS, np.inf)[None])
        t_p.append(t[m])
        print(f"    Ū {sg}: {m.sum()} pasos", flush=True)
    t = np.concatenate(t_p)
    o = np.argsort(t)
    return (np.concatenate(U_p)[o].astype(np.float32),
            np.concatenate(V_p)[o].astype(np.float32), t[o])


def procesar(scen, run):
    from netCDF4 import Dataset
    from conversion_baroclinica_2D import leer_grid_global
    from fondo_kang_fringer import leer_fondo_inicial
    from energy_flux_IWBay import MNC_GLOB, SCENARIOS, eos
    from conversion_baroclinica_2D import construir_series_por_tile
    from MIT_DatProcesW import discover_stage_patterns
    from prototipo_hovmoller_cpcg import bandpass_time
    p = SCENARIOS[scen]
    run_dir = p["base"] / run
    lab = "bay" if run == "run_expand" else "nobay"
    print(f"\n{'=' * 60}\n{scen} / {lab}: {run_dir}", flush=True)
    stages = list(discover_stage_patterns(run_dir))
    drF, hFacC, Depth, rAc, Y, X = leer_grid_global(run_dir)
    hW, hS, dyG, dxG = leer_malla_caras(run_dir)
    A_W, A_S = areas_bajo_cara(drF, hW, hS, dyG, dxG)
    rho_b = leer_fondo_inicial(run_dir, p)
    print(f"  ρ_b(z) inicial: {rho_b[0]:.4f} … {rho_b[-1]:.4f} kg m⁻³", flush=True)
    print("  Ū(t), V̄(t) ...", flush=True)
    Ub, Vb, t = series_ubar(run_dir, stages, drF, hW, hS, T_PULSO_FIN)
    dt = float(np.median(np.diff(t)))
    Uf = {b: (_filtrar(Ub, dt, rango), _filtrar(Vb, dt, rango)) for b, rango in BANDAS.items()}
    C = {n: np.zeros((len(Y), len(X))) for n in [*VENTANAS, *BANDAS]}
    en_vent = {n: (t >= a) & (t <= b) for n, (a, b) in VENTANAS.items()}
    nz = len(drF)
    for g in construir_series_por_tile(run_dir, MNC_GLOB, stages):
        j0, i0, tny, tnx = g["j0"], g["i0"], g["ny"], g["nx"]
        sj, si = slice(j0, j0 + tny), slice(i0, i0 + tnx)
        sjp, sip = slice(j0, j0 + tny + 1), slice(i0, i0 + tnx + 1)
        args = (A_W[:, sj, sip], A_S[:, sjp, si], rAc[sj, si])
        Ut, Vt = Ub[:, sj, sip], Vb[:, sjp, si]
        for k in range(nz):
            if not hFacC[k, sj, si].any():
                continue
            Wc = w_bt_centro(Ut, Vt, *args, k)
            T_p, S_p, tt = [], [], []
            for fi in g["files"]:
                with Dataset(fi) as ds:
                    tf = ds.variables["T"][:].data
                    m = tf >= T_PULSO_FIN
                    if not m.any():
                        continue
                    T_p.append(ds.variables["Temp"][m, k].data)
                    S_p.append(ds.variables["S"][m, k].data)
                    tt.append(tf[m])
            tt = np.concatenate(tt)
            o = np.argsort(tt)
            if tt.size != t.size or np.abs(tt[o] - t).max() > 1.0:
                raise ValueError(f"tiempos del tile ({tt.size}) ≠ los de Ū ({t.size})")
            rho = eos(np.concatenate(T_p)[o], np.concatenate(S_p)[o], p["tRef"],
                      p["sRef"], p["tAlpha"], p["sBeta"], p["rhoNil"])
            rp = rho - rho_b[k]
            peso = drF[k] * hFacC[k, sj, si]
            prod = rp * Wc
            for n, m in en_vent.items():          # sin filtrar, por ventana de tiempo
                C[n][sj, si] += G * prod[m].mean(0) * peso
            for b, rango in BANDAS.items():       # ρ′ y W_bt filtrados antes del producto
                Wf = w_bt_centro(Uf[b][0][:, sj, sip], Uf[b][1][:, sjp, si], *args, k)
                C[b][sj, si] += G * (bandpass_time(rp, dt, rango) * Wf).mean(0) * peso
        print(f"    tile j0={j0} i0={i0} listo", flush=True)
    print("  " + "   ".join(f"∫C_{n} = {np.nansum(c * rAc):.4e} W" for n, c in C.items()),
          flush=True)
    SAVE_DIR.mkdir(exist_ok=True)
    meta = dict(scenario=scen, run=lab, referencia="rho_b(z) inicial",
                metodo="C = <Σ ρ' g W_bt drF hFacC>_t, "
                "W_bt = -∇·((H+z)Ū) discreta (malla C)", t_min_s=T_PULSO_FIN,
                ventanas_s={n: [a, None if np.isinf(b) else b] for n, (a, b) in VENTANAS.items()},
                bandas_cph={b: list(r) for b, r in BANDAS.items()}, nt=int(t.size))
    out = SAVE_DIR / f"{scen}_{lab}_conversion_bt.npz"
    np.savez_compressed(out, X=X, Y=Y, Depth=Depth, rAc=rAc,
                        **{f"C_{n}": c for n, c in C.items()},
                        meta_json=np.array(json.dumps(meta, ensure_ascii=False)))
    print(f"  [OK] {out}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scen", nargs="+", default=["linear", "realistic", "2layer"])
    ap.add_argument("--run", nargs="+", default=RUNS)
    a = ap.parse_args()
    for s in a.scen:
        for r in a.run:
            procesar(s, r)
    print("\n[Listo] conversión con W barotrópica")


if __name__ == "__main__":
    main()
