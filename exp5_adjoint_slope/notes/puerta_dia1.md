# Puerta del día 1 (Task 3): resultado de GS

Fecha: 2026-11-16. Estado de partida: pickup final de G0 (30 años, configuración Dettling con relajación superficial; `v4_revision_posthoc.md` §6–7). Ventana de 10 días, J promediada en los últimos 5. Binarios Tapenade `adj`, `tlm` y `fwd` de `chan10_tap`. Corridas en ESTEBAN_MIT, copiadas a LaCie.

**Resultado: la puerta PASA.**

## 1. fc del adjunto contra el directo
`fc` = 6.32468237004301e-05 TW en el adjunto (fcref de grdchk) y en el directo (GSF). Son iguales en los 15 dígitos impresos.

## 2. grdchk (diferencias finitas centradas)

| Control | Punto (i, j, k) | eps | ∂J/∂u adjunto | DF | error rel. |
|---|---|---|---|---|---|
| T₀ (`xx_theta`) | 20, 18, 15 | 1e-1 / 1e-2 / **1e-3** | 3.850265e-4 | 4.476e-4 / 3.876e-4 / 3.850274e-4 | 16 % / 0.66 % / **2.3e-6** |
| κᵥ (`xx_diffkr`) | 20, 18, 10 | **1e-5** / 1e-6 / 1e-7 | −3.372185e-4 | −3.372216e-4 / … | **9.2e-6** / 5.2e-5 / 4.1e-4 |
| τx (`xx_fu`) | 20, 18, 1 | 1e-2 / 1e-3 / 1e-4 | 0.4749953 | 0.318 / 0.495 / 0.460 | 33 % / 4.3 % / 3.1 % |
| τx, barrido ampliado (GSU) | 20, 18, 1 | 1e-5 / **1e-6** | 0.4749953 | 0.4750223 / 0.4749819 | 5.7e-5 / **2.8e-5** |
| τx, tangente lineal (GSU_tlm) | 20, 18, 1 | 1e-4 | **0.4749953 (igual al adjunto)** | 0.460 | — |

- **τx en la fila de J.** Con el barrido de eps registrado en el plan (1e-2 a 1e-4) no se llega al 1 %.
  - J es un residuo pequeño (6e-5 TW) y el punto de τx está sobre la línea de J. Con eps = 1e-4 N m⁻², J cambia 150 %, así que la respuesta ya no es lineal.
  - El tangente lineal da el mismo gradiente que el adjunto en 7 dígitos, así que el código derivado es correcto.
  - Con eps 1e-5 y 1e-6 la DF converge: 5.7e-5 y 2.8e-5.
  - **Post hoc, declarado:** el barrido de eps se amplió después de ver el resultado. El criterio (< 1 % con el mejor eps) no cambia.
- **Otros puntos de τx** (diagnóstico, no son de la puerta):
  - j = 10, plataforma: 1.5 % con eps 1e-3 y **0.086 %** con 1e-4;
  - j = 30, mar afuera: 3.3 % y 1.5 %, todavía convergiendo.
  - Para V1 (Task 5), τx usa eps ≤ 1e-5 en todos sus puntos.

## 3. Ancla de J (GJ, 2 días)
- fc de Fortran: 4.93293e-5 TW. J de Python, desde los snapshots de cada paso: 4.93153e-5 TW.
- La diferencia absoluta es 1.4e-8 TW: 3.9e-6 del flujo bruto (|resuelto| + |bolo| + |Redi| = 0.0035 TW), y 5.8e-7 con el último registro extrapolado (el modelo no escribe snapshot del último paso).
- Respecto de J (≈ 0), la diferencia relativa es 2.8e-4. Es una cifra mal condicionada y se reporta junto con la escala del flujo bruto (`v4_revision_posthoc.md` §7).
- La J de Python y los flujos del propio modelo (`DFyE_TH`, `GM_vbT`) coinciden a 1e-15.

## 4. Lo que deja para las tareas siguientes
- **Estabilidad:** el adjunto de 10 días es estable, sin NaN, con la forma advectiva de GM y `GM_Small_Number = 1e-6`.
- **V1 (Task 5):** ≥ 10 puntos, con τx en eps ≤ 1e-5 y comparación con el tangente lineal. Se excluyen puntos sin señal (Review Focus 4).
- **J ≈ 0 en el estado medio** (V4, §7). Las sensibilidades se reportan con la partición resuelto/bolo/Redi.
