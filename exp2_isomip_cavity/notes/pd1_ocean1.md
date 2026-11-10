# PD1b — Ocean1 a 2 km: evaluación contra pd0 y Figura 1 (Exp. 2, Task 4 fase B)

Fecha: 2026-11-10. Evaluación mecánica contra `pd0_preregistro.md`, sin corridas nuevas. Script: `analysis/fig1_ocean0_validation.py` → `figures/fig1_ocean.json` y `figures/fig1_ocean.{png,pdf}`. Calibración y Ocean0: `pd1_calibracion.md` («Resultados R21»).

## Corrida

- `ocean1_2km_r21`: estado inicial COLD (AD16 Tabla 5), restauración WARM, ΓT_cal,2km = 0.030959, dt = 400 s, np = 5; 20 años en 4 etapas de 5 años con pickup. Cola de producción (`scripts/queue_production.sh`), fin 2026-11-10 06:45 CST; iteración final 1 576 800.
- Leída con `metrics.load_run` (`require_ok` en las 4 etapas, candado `melt_law`). Los 240 registros mensuales son contiguos, incluidas las 3 costuras entre etapas (registros 61, 121 y 181: el registro anterior termina donde empieza el siguiente, a 157 680 000, 315 360 000 y 473 040 000 s). **Año 20 = registros 229–240.**

## Veredictos (pd0)

| Criterio | Valor (año 20) | Rango de paso | Veredicto |
|---|---|---|---|
| V2a: melt medio de toda la plataforma | 10.34 m a⁻¹ | 4–12 | **PASA** |
| V2b: máx. \|ψ(x, z)\| en la cavidad (x < 640 km) | 0.185 Sv | 0.20–0.25 | **FALLA** |
| V2c: máx. \|ψ_bt\| en la cavidad (y > 50 km, x < 640 km) | 0.258 Sv | 0.09–0.43 | **PASA** |
| V2d: máx. \|ψ_bt\| cerca del frente (630 < x < 650 km; R13) | 0.492 Sv | 0.14–0.68 | **PASA** |
| V3 barotrópica: (m̄, C̄) = (10.34, 0.261 Sv) | banda en m̄: 0.255–0.460 | — | **PASA** |
| V3 overturning: (m̄, C̄) = (10.34, 0.186 Sv) | banda en m̄: 0.196–0.287 | — | **FALLA** |
| N1 (año 20 contra año 19) | 0.11 % | < 3 % | **PASA** |
| N2 (presupuesto del trazador, 20 años) | 1.3e-6 | < 1 % | **PASA** |
| N3 (`require_ok`, 4 etapas) | OK | — | **PASA** |

- **V3, banda.** pd0 da anclajes en m = 5 y 10 y hace evaluable [4, 12]. Como m̄ = 10.34 queda fuera de [5, 10], la banda del veredicto es la prolongación lineal de las rectas de los anclajes. Con la banda fijada en los valores del anclaje de 10 (0.25–0.45 y 0.19–0.28) el veredicto es el mismo en los dos tipos; el JSON lo registra (`variants_agree`).
- **V3, informativo (sin umbral):** ajuste con los 240 meses contra el ajuste multi-modelo impreso:
  - barotrópica: y = 0.0190 x + 0.0689, r = 0.974 (publicado: y = 0.024 x + 0.071, r = 0.89);
  - overturning: y = 0.0170 x + 0.0155, r = 0.993 (publicado: y = 0.02 x + 0.019, r = 0.97).
- **V2b y V3 overturning fallan juntos:** son la misma debilidad, una overturning ~7 % por debajo del límite inferior de V2b (0.185 contra 0.20) y ~5 % por debajo de la banda de V3. Según el plan, un V que falla va a **PD2**. No se redefine ningún criterio. **Clasificación (a)–(d): pendiente.** Hechos que la acotan, sin atribuir causa:
  - MITgcm-JPL, con el mismo modelo, ΓT = 0.0325 y otra configuración (advección 33, ajuste convectivo, hFacMin = 0.25, más celdas de plataforma), da 0.203 Sv en su año 20, dentro de V2b (comparación post hoc, `pd1_comparacion_jpl.md`).
  - La serie R20, que probaba esas opciones de una en una, está anulada por R21; no hay evidencia vigente sobre su efecto.

## Informativo (no son criterios)

- **CFL (R17):** máximo del año 20 = 0.230 (w), en 366 bloques de monitor diarios de la etapa 4 con t ≥ 19 años; máximo de los 20 años = 0.233. Ninguno supera 0.5.
- **Sobreenfriamiento (R15):**
  - estado final del año 20: 241 de 164 922 celdas mojadas bajo T_f(S, p); mínimo −0.126 °C; 10 por más de 0.1 °C;
  - media del año 20: 173 celdas; mínimo −0.122 °C; 9 por más de 0.1 °C.
  - En los dos casos todas son celdas superiores de columna, ninguna en las filas de pared. La más fría está en k = 9 (z = −190 m), hFacC = 0.46.
- **Melt profundo del año 20** (base < −300 m): 30.37 m a⁻¹. Es cercano al de Ocean0 (30.00, la meta de K1); Ocean1 restaura al mismo perfil WARM.
- **`ustar_check` del último registro:** el piso √Cd·u_tidal y la identidad γT = ΓT·u* pasan. La cota de Jensen falla por poco: razón mínima 0.9981 con tolerancia 1e-3, 3.8 % de las columnas por debajo de 1, mediana 1.004. La cota usa hFac estáticos y el modelo usa hFac con r*; con el u* corregido (menor), ese efecto pesa más que en el smoke. No es criterio y no se investiga aquí.

## Figura 1

`figures/fig1_ocean.{png,pdf}`:
- (a) melt del año 20 con la línea de apoyo y el borde de la plataforma;
- (b) serie de 20 años del melt medio, con N1;
- (c) V2 contra los rangos publicados;
- (d) diagrama melt–circulación sobre la banda de V3, etiquetada «read from Yung et al. 2026 Fig. 13 (approx.)»;
- (e) melt profundo de Ocean0 contra ΓT a 2 y 4 km, solo con las corridas de calibración R21.

El pie no presenta el acuerdo de melt entre 2 y 4 km como evidencia independiente: es la meta de la calibración (K1). Cita en su lugar ψ_bt (+125 %) y la correlación de los mapas (r = 0.91).
