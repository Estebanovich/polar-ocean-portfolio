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
| *V2b con la definición de Yung 2026 (máx. de \|ψ(x, z)\| en todo el dominio, año 20). **Post hoc**, no es veredicto de pd0* | *0.222 Sv (x = 766 km, z = −220 m); JPL: 0.234* | *0.20–0.25* | *dentro (post hoc)* |
| V2c: máx. \|ψ_bt\| en la cavidad (y > 50 km, x < 640 km) | 0.258 Sv | 0.09–0.43 | **PASA** |
| V2d: máx. \|ψ_bt\| cerca del frente (630 < x < 650 km; R13) | 0.492 Sv | 0.14–0.68 | **PASA** |
| V3 barotrópica: (m̄, C̄) = (10.34, 0.261 Sv) | banda en m̄: 0.255–0.460 | — | **PASA** |
| V3 overturning: (m̄, C̄) = (10.34, 0.186 Sv) | banda en m̄: 0.196–0.287 | — | **FALLA** |
| N1 (año 20 contra año 19) | 0.11 % | < 3 % | **PASA** |
| N2 (presupuesto del trazador, 20 años) | 1.3e-6 | < 1 % | **PASA** |
| N3 (`require_ok`, 4 etapas) | OK | — | **PASA** |

- **V3, banda.** pd0 da anclajes en m = 5 y 10 y hace evaluable [4, 12]. Como m̄ = 10.34 queda fuera de [5, 10], la banda del veredicto es la prolongación lineal de las rectas de los anclajes; en la Figura 1(d) esos tramos (4–5 y 10–12 m a⁻¹) van rayados como extrapolados, y el punto del año 20 cae en uno de ellos. Con la banda fijada en los valores del anclaje de 10 (0.25–0.45 y 0.19–0.28) el veredicto es el mismo en los dos tipos; el JSON lo registra (`variants_agree`).
- **V3, informativo (sin umbral):** ajuste con los 240 meses contra el ajuste multi-modelo impreso:
  - barotrópica: y = 0.0190 x + 0.0689, r = 0.974 (publicado: y = 0.024 x + 0.071, r = 0.89);
  - overturning: y = 0.0170 x + 0.0155, r = 0.993 (publicado: y = 0.02 x + 0.019, r = 0.97).
## Decisiones de PD2 (Esteban, 2026-11-10)

- **V2b: FALLA, tal como está pre-registrado** (cavidad, x < 640 km: 0.185 Sv contra 0.20–0.25; ~7 % bajo el límite inferior). No se redefine.
  - **Error (a) nuestro, en la operacionalización del criterio, no en el modelo.** El pre-registro tomó mal la fuente. Y26 L1099–1109 da los 0.20–0.25 Sv como el máximo de la overturning de la mayoría de los modelos, y describe picos tanto dentro de la cavidad como cerca de la frontera «norte» (x = 800 km). El pie de la Fig. 8 (L1182) muestra la streamfunction en todo el dominio. La restricción x < 640 km pertenece a la Fig. 13 (fuerza de la circulación, V3) y se trasladó por error a V2b.
  - **Fila post hoc**, rotulada como tal en `figures/fig1_ocean.json` (`V2b_posthoc_full_domain`), en la tabla de arriba y en la Figura 1(c) (rombo abierto y pie). Con la definición de Y26, el máximo de \|ψ(x, z)\| del año 20 en todo el dominio es **0.222 Sv** (x = 766 km, z = −220 m), dentro de 0.20–0.25. El de MITgcm-JPL, con la misma definición, es 0.234 Sv.
  - Esto no apoya que V2b y V3 overturning sean «la misma debilidad». En todo el dominio la overturning está dentro del rango publicado; solo la celda de la cavidad es débil. No se atribuye causa.
- **V3 overturning: FALLA, marginal, no resuelta dentro de la incertidumbre de lectura de la figura.** C̄ = 0.1865 Sv:
  - 4.9 % (0.0096 Sv) por debajo de la banda extendida (límite 0.1961 en m̄ = 10.34);
  - 1.8 % (0.0035 Sv) por debajo de la banda fijada en el anclaje de 10 m a⁻¹ (límite 0.19).

  Los anclajes son lecturas aproximadas de una figura y m̄ cae en el tramo extrapolado. **No se asigna clase (a)–(d).**
- Contexto post hoc (`pd1_comparacion_jpl.md`): MITgcm-JPL, con el mismo modelo, ΓT = 0.0325 y otra configuración (advección 33, ajuste convectivo, hFacMin = 0.25, más celdas de plataforma), da 0.203 Sv en la cavidad y 0.234 Sv en todo el dominio. La serie R20, que probaba esas opciones de una en una, está anulada por R21.

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
