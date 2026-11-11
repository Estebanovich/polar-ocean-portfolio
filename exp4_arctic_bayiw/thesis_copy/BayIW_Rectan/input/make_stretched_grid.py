#!/usr/bin/env python3
"""
make_stretched_grid.py
======================
Aplica ESTIRAMIENTO GEOMETRICO a los margenes de las fronteras ABIERTAS
(sur, este, oeste) de los vectores de espaciado de malla (delXFile/delYFile),
dejando intactos el interior y la frontera NORTE (costa cerrada).

Opera SOBRE los binarios dx/dy que produce el notebook de batimetria
(bahia_rectan_impar_func.ipynb), de malla ~5 km. De este modo el interior
queda EXACTAMENTE como la malla sobre la que se interpolo la batimetria
(incluida la pequena banda de ~5053 m del centro), y solo se modifican las
celdas externas. Por eso debe ejecutarse DESPUES del notebook.

Motivacion (fronteras "realmente abiertas")
-------------------------------------------
La onda barotropica viaja a c = sqrt(gH) ~ 75 m/s (H~590 m). En la malla
uniforme de 5 km su velocidad de fase no-dimensional es

    CFL = c*dt/dx = 75 * 30 / 5000 ~ 0.45,

justo en el limite de estabilidad de la radiacion de Orlanski (CMAX < 0.5,
ver pkg/obcs/orlanski_*.F): radiacion marginal. Al engrosar las celdas hacia
el borde abierto la CFL baja muy por debajo de 0.45 (Orlanski radia con
holgura) y se crea un margen disipativo que aleja la frontera de la bahia.
Combinado con la capa esponja (data.obcs), reduce las reflexiones.

Por que es seguro (no toca el resto del setup)
----------------------------------------------
  - T y S son horizontalmente uniformes -> el estiramiento no los altera.
  - El margen sur es plano y profundo (-600 m) -> estirar dy no distorsiona.
  - Los margenes E/O tienen oceano profundo al sur y costa al norte, lejos
    de la bahia (i~274-285); estirarlos solo ensancha el margen y extiende
    la costa, sin afectar la dinamica de la bahia.
  - El forzamiento de viento es un parche centrado en la bahia, nulo en los
    margenes.
  - Se conserva Nx=560, Ny=352 (SIZE.h) y el indexado de celdas.

El factor por celda (5 % = 1.05) cumple el criterio usual de suavidad de
malla estirada (<10 % entre celdas contiguas).

Uso
---
    python3 make_stretched_grid.py      # ejecutar UNA vez tras el notebook
"""

import numpy as np

NSTRETCH = 20            # nº de celdas estiradas por margen abierto
FACTOR   = 1.05          # razon geometrica de crecimiento por celda
DT       = 30.0          # deltaT (s), solo para reportar CFL
H_DEEP   = 600.0         # profundidad del margen (m), solo para CFL
DTYPE    = np.dtype('>f8')
PREFIXES = ('bahia_01_expand', 'nobahia_01_expand')


def stretch_low(vec, n, factor):
    """Estira las primeras n celdas (borde en i=0): de FINO (interior) a
    GRUESO (borde). Ancla al primer valor interior vec[n]."""
    v = vec[n]
    for i in range(n):
        vec[i] = v * factor ** (n - i)     # i=n-1 -> v*factor ; i=0 -> v*factor**n
    return vec


def stretch_high(vec, n, factor):
    """Estira las ultimas n celdas (borde en i=N-1). Ancla al ultimo
    valor interior vec[N-1-n]."""
    N = vec.size
    v = vec[N - 1 - n]
    for k in range(n):
        vec[N - n + k] = v * factor ** (k + 1)   # k=0 -> v*factor ; k=n-1 -> v*factor**n
    return vec


def report(prefix, dx, dy):
    g = 9.81
    c = (g * H_DEEP) ** 0.5
    print(f"  [{prefix}]")
    print(f"    dx: interior={dx[NSTRETCH]:.0f} m  borde={dx[0]:.0f} m  "
          f"(CFL {c*DT/dx[NSTRETCH]:.2f} -> {c*DT/dx[0]:.2f})")
    print(f"    margen X = {dx[:NSTRETCH].sum()/1e3:.0f} km/lado ; "
          f"dominio X={dx.sum()/1e3:.0f} km Y={dy.sum()/1e3:.0f} km")
    jmax = max(abs(dx[1:]/dx[:-1] - 1).max(), abs(dy[1:]/dy[:-1] - 1).max())
    print(f"    salto max entre celdas = {jmax*100:.1f}% (criterio <10%)")


def main():
    c = (9.81 * H_DEEP) ** 0.5
    print(f"Estiramiento de margenes abiertos (S/E/O), N={NSTRETCH}, factor={FACTOR}")
    print(f"c=sqrt(gH)={c:.1f} m/s ; norte e interior intactos")
    for prefix in PREFIXES:
        dx = np.fromfile(f'{prefix}_dx.bin', DTYPE).copy()
        dy = np.fromfile(f'{prefix}_dy.bin', DTYPE).copy()
        # X: estira oeste (low) y este (high)
        stretch_low(dx, NSTRETCH, FACTOR)
        stretch_high(dx, NSTRETCH, FACTOR)
        # Y: estira solo el sur (low); norte sin tocar
        stretch_low(dy, NSTRETCH, FACTOR)
        report(prefix, dx, dy)
        dx.astype(DTYPE).tofile(f'{prefix}_dx.bin')
        dy.astype(DTYPE).tofile(f'{prefix}_dy.bin')
        print(f"    escrito {prefix}_dx.bin ({dx.size}), {prefix}_dy.bin ({dy.size})")


if __name__ == '__main__':
    main()
