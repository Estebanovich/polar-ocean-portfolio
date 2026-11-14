#!/usr/bin/env python3
"""
coords_iwbay.py
===============
Origen de coordenadas del dominio IWBay: **eje de la ensenada y línea de costa**.

Las coordenadas físicas del MITgcm (`XC`/`YC` del `grid.nc`, y por herencia los
`X_zoom`/`Y_zoom` de todos los npz) NO tienen el origen en la bahía: la costa
cae en y ≈ +76.5 km y el eje de la ensenada en x ≈ +76.1 km.

La causa está en `&PARM04` de los `data`: `xgOrigin=-1400E3` e `ygOrigin=-1621E3`
son números redondos escritos a mano, pero la malla es telescópica
(`delXFile`/`delYFile`, celdas de 5 a 13.27 km) y lo que acumula desde el borde
del dominio hasta la bahía es 1476.12 km en X y 1697.49 km en Y. El sobrante es
exactamente el corrimiento. No es un error del post-proceso: reconstruir la
malla desde los binarios reproduce los `X_zoom`/`Y_zoom` de los npz con error
máximo de 4e-9 m. Es cosmético y se corrige aquí, sin re-correr el modelo.

Los offsets se DEDUCEN de la máscara, no se escriben a mano, para que sigan
siendo correctos si cambia la malla:
  - la costa abierta es la fila de agua más al norte que se repite en casi
    todas las columnas (j = 323; la tierra empieza en j = 324)
  - la ensenada son las columnas donde el agua sube muy por encima de esa fila
    (i = 274…285)
  - el origen va en la ARISTA entre la última celda de agua y la primera de
    tierra, no en el centro de una celda

`geometria_boca()` da lo mismo que `rayos_boca_difraccion.mouth_geometry()`,
que ya trabajaba en coordenadas relativas y por eso nunca sufrió el corrimiento.

Sólo depende de numpy, así que se puede importar desde cualquier entorno (los
scripts de mapas necesitan `MITgcmutils`, que sólo está en el env base).

Uso:
    from coords_iwbay import offsets_recentrado, recentrar
    x0, y0 = offsets_recentrado()          # km
    Xr, Yr = recentrar(X_km, Y_km)

    python coords_iwbay.py                 # imprime los offsets y las estaciones
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

_BASE = Path(__file__).resolve().parent
NPZ_DIR = _BASE / "npz_out"
# cualquiera de los tres sirve: los tres escenarios comparten malla y batimetría
_CANDIDATOS = ("2layer", "linear", "realistic")

_cache: dict[str, tuple] = {}


def _mascara(npz: Path | None = None):
    """(X_km, Y_km, mask) del primer npz de Eta que exista. mask True = tierra."""
    if npz is None:
        for s in _CANDIDATOS:
            p = NPZ_DIR / f"{s}_Eta_zoom_sfc_bay.npz"
            if p.exists():
                npz = p
                break
    if npz is None or not Path(npz).exists():
        raise FileNotFoundError(
            f"no encontré ningún {{scen}}_Eta_zoom_sfc_bay.npz en {NPZ_DIR}")
    with np.load(npz, allow_pickle=True) as d:
        X = np.asarray(d["X_zoom"], float) / 1e3
        Y = np.asarray(d["Y_zoom"], float) / 1e3
        M = np.asarray(d["mask_zoom"]).astype(bool)
    if X.ndim == 2:
        X = X[0, :]
    if Y.ndim == 2:
        Y = Y[:, 0]
    return X, Y, M


def geometria_boca(npz: Path | None = None) -> tuple[float, float, float]:
    """(x0, y0, a) en km: eje de la ensenada, línea de costa y ancho de la boca.

    El ancho `a` es de borde a borde (incluye la celda completa de cada extremo),
    como en `rayos_boca_difraccion.mouth_geometry`.
    """
    clave = str(npz) if npz else "_default"
    if clave in _cache:
        return _cache[clave]
    X, Y, M = _mascara(npz)
    nx = M.shape[1]
    # fila de agua más al norte de cada columna; la que más se repite es la costa
    top = np.array([np.flatnonzero(~M[:, i]).max() if (~M[:, i]).any() else -1
                    for i in range(nx)])
    if not (top >= 0).any():
        raise RuntimeError("la máscara no tiene agua")
    j_costa = int(np.bincount(top[top >= 0]).argmax())
    ens = np.flatnonzero(top > j_costa + 2)
    if ens.size == 0:
        raise RuntimeError("no encontré la ensenada en la máscara")
    dx = float(np.median(np.diff(X)))
    x0 = 0.5 * (X[ens.min()] + X[ens.max()])
    a = float(X[ens.max()] - X[ens.min()] + dx)
    y0 = 0.5 * (Y[j_costa] + Y[j_costa + 1])
    _cache[clave] = (x0, y0, a)
    return x0, y0, a


def offsets_recentrado(npz: Path | None = None) -> tuple[float, float]:
    """(x0, y0) en km a RESTAR para poner el eje de la ensenada y la costa en 0."""
    x0, y0, _ = geometria_boca(npz)
    return x0, y0


def recentrar(X_km, Y_km, npz: Path | None = None):
    """(X − x0, Y − y0) en km. Acepta arrays o escalares."""
    x0, y0 = offsets_recentrado(npz)
    return np.asarray(X_km, float) - x0, np.asarray(Y_km, float) - y0


# Rótulos de eje coherentes en todas las figuras que recentran
XLABEL = "X (km)"
YLABEL = "Y (km)"
NOTA_ORIGEN = "0 = eje de la ensenada (X) y línea de costa (Y)"


if __name__ == "__main__":
    x0, y0, a = geometria_boca()
    print(f"eje de la ensenada  x0 = {x0:9.3f} km")
    print(f"línea de costa      y0 = {y0:9.3f} km")
    print(f"ancho de la boca    a  = {a:9.3f} km")
    X, Y, _ = _mascara()
    print("\nestaciones recentradas (km):")
    for nm, (iy, ix) in {"Seiche": (345, 150), "Omega": (315, 150),
                         "A": (315, 100), "B": (315, 200),
                         "C": (265, 150)}.items():
        print(f"  {nm:7s} ({X[ix] - x0:8.1f}, {Y[iy] - y0:8.1f})")
