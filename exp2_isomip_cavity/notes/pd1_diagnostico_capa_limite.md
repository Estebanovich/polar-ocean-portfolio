# Diagnóstico post hoc: ¿por qué nuestro ΓT_cal queda por debajo del de los otros MITgcm de ISOMIP+? (ruling R20)

Fecha: 2026-11-09. Aprobado por Esteban (R20). Escrito **antes** de lanzar cualquier corrida de este diagnóstico.

## Estatus

- Este diagnóstico es **post hoc**. **No toca K1, V1 ni ΓT_cal.**
  - K1 pasa a 2 y 4 km.
  - V1 falla (ΓT_cal,2km = 0.0099627 < 0.011), con clase (c), según `pd1_calibracion.md` §4.2–4.3.
  - ΓT_cal,2km = 0.0099627 y ΓT_cal,4km = 0.0094686 quedan fijos.
  - Ningún resultado de aquí cambia la configuración del ensamble ni la de Ocean1.
- **Objetivo:** explicar la diferencia. No busca corregirla.

## Pregunta

Nuestro ΓT_cal (≈ 0.0095–0.010, con celdas de 20 m) es menor que el de los otros MITgcm de Y26 a 20 m: MITgcm-JPL 0.0325 y MITgcm-BAS-Coupled 0.0135 (MITgcm-BAS reporta 0.011 en Ocean1). A igual ΓT, nuestro modelo funde más. Qué hacen distinto esas entradas, según Y26:
- **Advección de trazadores:** DST3 (esquema 33; final L751; preprint L873). Nosotros usamos 77.
- **Muestreo y reparto de flujos:** "sampling the temperature, salinity, and velocity within a fixed 20 m boundary layer beneath the ice. Melt fluxes are also distributed over a 20 m layer" (final L757–759; preprint L880–881).
- **Columna mínima:** MITgcm-BAS mantiene 40 m excavando el calado (final L764; preprint L887).
- **Celdas parciales:** mínimo de 25 % (JPL) y 20 % (BAS, BAS-Coupled) (final L768–769; preprint L889–890). Nosotros usamos 0.20.

## Diseño

- **Corridas:** Ocean0 a 4 km, ΓT = 0.015 (ΓS = ΓT/35), 2 años, dt = 800 s, np = 1. Cada variante tiene **exactamente un cambio** respecto a `cal4_gT0.015`.
- **Referencia:** `cal4_gT0.015`, ya corrida, con melt profundo de los meses 19–24 = 41.216 m a⁻¹.
  - Sus entradas binarias son idénticas a las que genera hoy `make_case`.
  - La única diferencia de namelist es `checkIniTemp=.FALSE.` (R18), sin efecto numérico (`pd1_calibracion.md` §6). Se comprobó con `diff` de las entradas generadas.

| Variante | Caso | Cambio único | Hipótesis |
|---|---|---|---|
| (a) | `diag4_adv33` | `tempAdvScheme = saltAdvScheme = PTRACERS_advScheme(1) = 33` (DST3), en lugar de 77 | DST3 es menos difusivo que el limitador de 2.º orden y conserva mejor la capa de agua de deshielo fría y dulce bajo el hielo. Eso aísla la base del agua cálida y baja el melt a igual ΓT. Si el efecto es grande, explica parte del ΓT mayor de las entradas MITgcm de Y26, que usaban 33. |
| (b) | `diag4_hfac025` | hFacMin = 0.25, en lugar de 0.20 (como MITgcm-JPL). `make_case` reajusta la geometría a las celdas parciales del modelo: 97 columnas cambian de calado y 102 de fondo, ≤ 4 m; 41 221 celdas mojadas contra 41 249 | Con menos celdas superiores muy delgadas (4 m), el promedio de la capa límite de Losch pondera distinto la celda delgada, y cambia el T muestreado. Se espera un efecto pequeño. |
| (c) | `diag4_bl20` | Capa fija de 20 m para muestrear y repartir flujos | **No se corre.** Ver abajo. |

**Por qué (c) no se corre:** nuestra configuración base ya usa esa capa, y las opciones de MITgcm no ofrecen otra sin código nuevo.
- `SHELFICEboundaryLayer=.TRUE.` ya está activo (`config/ocean0/input_templates/data.shelfice`). Con él, T, S y u se promedian sobre una capa de espesor drF(kTop) = 20 m bajo el hielo: la celda superior, más la fracción de la siguiente que completa los 20 m.
  - Original: `pkg/shelfice/shelfice_thermodynamics.F:341-389` (bloque `IF ( SHELFICEBoundaryLayer )`).
  - Nuestra copia: `config/ocean0/code/shelfice_thermodynamics.F:373-410`. drKp1 = drF(K)(1 − hFacC(K)), acotado por drF(K+1)·hFacC(K+1).
  - Los flujos se reparten sobre esa misma capa: `pkg/shelfice/shelfice_forcing.F:74-101, 180-206`.
  - `pkg/shelfice/SHELFICE.h:30-31` la describe como "vertical merging of cells ... for a boundary layer of drF thickness". Con delR = 20 m en todos los niveles, esa capa es de 20 m fijos.
- Las demás opciones de `shelfice_readparms.F` no definen un espesor distinto:
  - `SHI_withBL_realFWflux` (`:100, :210-211`) solo actúa con `useRealFreshWaterFlux`, que tenemos en `.FALSE.`;
  - `SHI_withBL_uStarTopDz` (`:101, :212-213`) solo cambia cómo se promedia u para u* (en los puntos u/v, `shelfice_thermodynamics.F:424-470` de nuestra copia), no el espesor de muestreo ni el de reparto.
  - No hay ningún parámetro de espesor de capa en `SHELFICE_PARM01` (`shelfice_readparms.F:44-75`).
- **Consecuencia:** la diferencia "capa fija de 20 m" de Y26 no distingue nuestra configuración de la suya, así que no puede explicar el hueco.
  - La única diferencia posible está en el detalle: en columnas donde la celda siguiente también es parcial, la capa queda por debajo de 20 m.
  - Probar otro espesor (p. ej., 40 m o el kTop completo) exigiría código nuevo. Por la regla de R20, no se hace.
  - `SHI_withBL_uStarTopDz` es un cambio disponible, pero no lo pidió R20 y **no se corre** (no se añaden variantes después de escribir esta nota).

## Lectura

Para cada variante v (con 2 años y N3 OK; si no, se reporta el hueco):
1. **Cambio del melt profundo** de los meses 19–24 a ΓT = 0.015, relativo a la referencia: Δ = (m_v − 41.216)/41.216.
2. **ΓT_cal implícito (estimación, no calibración):** ΓT_est = 0.015·(30/m_v)^(1/p), con p = 0.690.
   - p es la pendiente log–log local medida a 4 km en el tramo [0.0075, 0.015]: ln(41.216/25.541)/ln 2.
   - Supone que la variante desplaza ln⟨m⟩ sin cambiar la pendiente.
   - Con la referencia, la fórmula da 0.00947, igual al ΓT_cal,4km calibrado (0.0094686).
3. **"Explica la diferencia"** si ΓT_est ≥ 0.0135 (el ΓT de MITgcm-BAS-Coupled), es decir, m_v ≤ 32.3 m a⁻¹: una reducción del melt profundo ≥ 21.7 %.
   - Para llegar al de MITgcm-JPL (0.0325) haría falta m_v ≤ 17.6 m a⁻¹ (−57 %). Se reporta como referencia, no como criterio.
4. **Informativo:** melt medio de la plataforma, N1, N2, CFL máximo (R17) y sobreenfriamiento (R15) de cada variante.

Si ninguna variante llega a 0.0135, el resultado se reporta como "no explicado por (a) ni (b)". Las demás diferencias citadas en Y26 se dejan anotadas, sin probarlas aquí:
- convección de Cessi y Young (1996), frente a nuestro `ivdc_kappa` (final L754–756);
- la excavación de la columna mínima de MITgcm-BAS;
- la incertidumbre que el propio Y26 reconoce sobre los ΓT de BAS (final L794–812).

## Ejecución

- **Cola:** `exp2_isomip_cavity/scripts/queue_diag_bl.sh`, desacoplada. np = 1, hasta 3 simultáneas (aquí 2), con las mismas reglas de seguridad.
- **Evaluación:** con `analysis/pd1_calibration_eval.py` (N1–N3, CFL, sobreenfriamiento) y `analysis/calibrate.py --verify` (melt profundo).

## Variante (d): `SHI_withBL_uStarTopDz` (adición post hoc, aprobada por Esteban el 2026-11-09; escrita a las 06:59 CST, antes de correrla)

- **Cambio único** frente a `cal4_gT0.015`: `SHI_withBL_uStarTopDz=.TRUE.` en `data.shelfice`, nada más. El `diff` de las entradas generadas solo muestra esa línea (y los comentarios con el nombre del caso). Caso `diag4_ustartopdz`, Ocean0 a 4 km, ΓT = 0.015, 2 años, np = 1.
- **El cambio sí llega a nuestro u\*.** Lo comprobé antes de correr, en `config/ocean0/code/shelfice_thermodynamics.F`, que es la copia que enlazan `ocean0_2km` y `ocean0_4km`:
  1. El bloque de la capa límite (`IF ( SHELFICEBoundaryLayer )`, :373–421) deja velSq = u² + v² de las velocidades promediadas en el centro de la celda (:417).
  2. Después, el bloque `IF ( SHI_withBL_uStarTopDz )` (:424–481) vuelve a calcular uLoc y vLoc. Primero promedia en la vertical, en los puntos U/V, sobre la capa drF (:438–448 para U, y su análogo para V). Luego promedia en la horizontal y **reescribe velSq** (:476–478).
  3. El u\* de ISOMIP+ usa esa velSq: `ustarSq = shiCdragfld·velSq + shiUTidal²` (:509–510), y después γT = ΓT·u\* y γS = ΓS·u\*.
  4. El bloque es idéntico al del original (`pkg/shelfice/shelfice_thermodynamics.F:392`, comprobado con `diff`).
  5. La opción se activa porque `SHELFICEboundaryLayer` y `SHELFICEuseGammaFrict` están en `.TRUE.` (`pkg/shelfice/shelfice_readparms.F:212-213`), y porque `SHI_ALLOW_GAMMAFRICT` está definido.

  Por tanto, la opción no es nula en nuestra modificación. T y S muestreados no cambian; solo cambia el |u| que entra en u\*.
- **Hipótesis:** el orden del promedio, primero en la vertical en las caras U/V y luego en la horizontal, frente a primero en la horizontal en el centro, cambia u\* donde la base del hielo tiene escalones (kTop distinto entre columnas vecinas, celdas parciales en las caras). Eso cambia γT y el melt a igual ΓT. **El signo del efecto no se conoce a priori.**
- **Lectura:** la misma regla de arriba, ya fijada:
  - Δ del melt profundo de los meses 19–24 frente a 41.2156 m a⁻¹;
  - ΓT_est = 0.015·(30/m)^(1/0.690), etiquetado como estimación;
  - "explica la diferencia" si ΓT_est ≥ 0.0135, es decir, m ≤ 32.3 m a⁻¹.
- **Informativo:** melt de la plataforma, N1, N2, CFL (R17) y sobreenfriamiento (R15).
- **Cola:** `scripts/queue_diag_bl.sh`, relanzada. `--skip-ok` salta las dos variantes ya hechas.

## Resultado (post hoc; cola 2026-11-09 03:21–03:48 CST, commit e898d77)

Las dos variantes terminaron OK, con N3 aplicado por `metrics.load_run`. Lectura mecánica según la regla de arriba: p = 0.690, referencia 41.2156 m a⁻¹, umbral ΓT_est ≥ 0.0135. Valores en `figures/diag_bl_<caso>.json`.

| Variante | Caso | Melt profundo, meses 19–24 (m a⁻¹) | Δ frente a la referencia | ΓT_est (estimación) | ¿Explica la diferencia (≥ 0.0135)? |
|---|---|---|---|---|---|
| referencia | `cal4_gT0.015` | 41.2156 | — | 0.00947 | — |
| (a) DST3 (33) | `diag4_adv33` | 39.4782 | **−4.2 %** | **0.0101** | **NO** |
| (b) hFacMin 0.25 | `diag4_hfac025` | 41.1641 | **−0.12 %** | **0.0095** | **NO** |
| (c) capa fija de 20 m | — | no se corrió (ya es la base) | — | — | no aplica |

- **Informativo:**

  | Caso | Melt de la plataforma (m a⁻¹) | N1 | N2 | CFL máx. (R17) | Sobreenfriamiento (estado final, ΓT = 0.015) |
  |---|---|---|---|---|---|
  | `cal4_gT0.015` | 13.595 | 0.0014 | 4.4e-8 | 0.262 | — |
  | `diag4_adv33` | 13.233 | 0.0031 | 5.6e-8 | 0.284 | 262 celdas, mín −0.354 °C |
  | `diag4_hfac025` | 13.640 | 0.0045 | 5.1e-8 | 0.257 | 257 celdas, mín −0.300 °C |

  Las salidas de código 4 de `--verify` son esperadas: K1 no se aplica a estas corridas.
- **Conclusión mecánica:** ninguna de las dos variantes probadas explica la diferencia. DST3 baja el melt profundo un 4 %, y el ΓT_cal implícito sube ≈ 6 % (0.0095 → 0.0101). hFacMin = 0.25 no lo cambia de forma apreciable. Queda **"no explicado por (a) ni (b)"**.

### Contexto (no es veredicto)

- El ΓT corregido de MITgcm-BAS en Ocean1 es **0.011**: "MITgcm-BAS used 0.011 in the Ocean1 experiment and 0.036 in the Ocean2 experiment (Table 2)" (Y26 final, L806–807).
- Por tanto, nuestro ΓT_cal,2km = 0.00996 queda **≈ 10 % por debajo del de MITgcm-BAS** (0.00996/0.011 = 0.91) y **≈ 3× por debajo del de MITgcm-JPL** (0.0325/0.00996 = 3.3).
- **El propio Y26 no pudo verificar** qué coeficientes usaron los experimentos de MITgcm-BAS: "Due to a lack of model data traceability during the intervening years, we are unable to verify which transfer coefficients these experiments used" (Y26 final, L807–809).
- **Candidatos restantes, no probados:**
  1. `SHI_withBL_uStarTopDz`: u promediado en los puntos u/v de la capa límite, para u* (`shelfice_readparms.F:101, :212-213`).
  2. La columna mínima de 40 m con calado excavado de MITgcm-BAS (Y26 final L764). La nuestra seca las columnas de menos de 2 celdas.
  3. Otros, ya anotados arriba: la convección de Cessi y Young (Y26 final L754–756) y la incertidumbre sobre los ΓT de BAS.
- No se lanzan más corridas (R20, cerrado).
